/*
 *  Exact solver for the Capacitated Location-Routing and related problems
 *
 *  @author Ruslan Sadykov <Ruslan.Sadykov@inria.fr>, Inria Bordeaux, France (coordinator and main contributor),
 *  @author Eduardo Uchoa <eduardo.uchoa@gmail.com>, UFF, Brazil (scientific advise)
 *  @author Pedro Henrique Liguori <phliguori@gmail.com> BCG GAMMA, Brazil (separation of DCCs)
 *  @author Guillaume Marques <guillaume.marques@protonmail.com>, Atoptima, France (separation of GUBs and FCCs)
 *
 *  @date 31/03/2022
 *
 *  Dependencies:
 *  - BaPCod v.0.69
 *  - VRPSolver extension (RCSP solver) v.0.5.14
 */

#ifndef LRP_CUTSEPARATION_H
#define LRP_CUTSEPARATION_H

#include "Data.h"
#include "Parameters.h"
#include "bcModelCutConstrC.hpp"
#include <random>
#include <utility>

namespace lrp
{
    #define LRP_MAX_NUMBER_OF_CUSTOMERS 512

	class GeneralizedUpperBoundSeparationRoutine : public BcCutSeparationFunctor
	{
		const Data & data;
        const Parameters & params;
		int cutCount;
		int cutRound;

	public:
        GeneralizedUpperBoundSeparationRoutine(const Data & data_, const Parameters & params_);

		~GeneralizedUpperBoundSeparationRoutine() override;

		int operator()(BcFormulation bcForm, BcSolution &primalSol, double & maxViolation,
                       std::list<BcConstr> & cutList) override;
	};

    class DepotCoverInequalitySeparationRoutine : public BcCutSeparationFunctor
    {
        const Data & data;
        const Parameters & params;
        int cutCount;
        int cutRound;
        int printLevel;
        BapcodInit * bcInitPtr;
        Model * modelPtr;
        std::vector<ProbConfig *> probConfPts;

    public:
        DepotCoverInequalitySeparationRoutine(const Data & data_, const Parameters & params_);

        ~DepotCoverInequalitySeparationRoutine() override;

        int operator()(BcFormulation bcForm, BcSolution &primalSol, double & maxViolation,
                       std::list<BcConstr> & cutList) override;
    };

    class FenchelCutSeparationRoutine : public BcCutSeparationFunctor
    {
        const Data & data;
        const Parameters & params;
        int cutCount;
        int cutRound;
        int printLevel;
        BapcodInit * bcInitPtr;
        Model * modelPtr;
        ProbConfig * probConfPtr;
        std::vector<std::vector<int> > minimalCovers;

    public:
        FenchelCutSeparationRoutine(const Data & data_, const Parameters & params_);

        ~FenchelCutSeparationRoutine() override;

        int operator()(BcFormulation bcForm, BcSolution &primalSol, double & maxViolation,
                       std::list<BcConstr> & cutList) override;
    };

    struct DepotCapacityCutSeparationSeparatorParameters {
        int maxNumPerRound;
        double cutViolationTolerance;
        int printLevel;
        int setRcandidateMinSize;
        int setRcandidateMaxSize;
        int setRselectionMethod; /// 0 - as in Liguori's code (just shuffling and taking a certain number of first ones),
                                 /// 1 - only sets in which all corresponding y variables are positive
        int heuristicChoice; /// 0 - Liguori's GRASP-based heuristic with thresholds,
                             /// 1 - greedy construction heuristic, similar to the RCC one
        /// parameterisation for the Liguori's heuristic
        bool useDiversification;
        double alphaParamInInitialGRASP;
        double alphaParamInDiversificationGRASP;
        int maxNbSwapsInInitialGRASP;
        double nbSwapsFactorInDiversificationGRASP;
        int maxNbSwapsAfterDiversification;
        int maxInsertDeleteIterations;
        int greedyConstructionHeuristicCriterion; /// 0 - sort by the LHS change, 1 - sort by the violation change
        bool genCutsWithTwoYs;

        DepotCapacityCutSeparationSeparatorParameters() :
                maxNumPerRound(100),
                cutViolationTolerance(1e-3),
                printLevel(2),
                setRcandidateMinSize(1),
                setRcandidateMaxSize(2),
                setRselectionMethod(1),
                heuristicChoice(0),
                useDiversification(true),
                alphaParamInInitialGRASP(0.15),
                alphaParamInDiversificationGRASP(0.3),
                maxNbSwapsInInitialGRASP(10),
                nbSwapsFactorInDiversificationGRASP(0.05),
                maxNbSwapsAfterDiversification(20),
                maxInsertDeleteIterations(2),
                greedyConstructionHeuristicCriterion(0),
                genCutsWithTwoYs(true)
        {}
    };

    struct DepotCapacityCutCandidate
    {
        int depotId1;
        int depotId2;
        std::set<int> setR;
        std::vector<bool> isInS;
        double violation;

        DepotCapacityCutCandidate(int depotId1_, int depotId2_, const std::vector<int> & setR_,
                                  std::vector<bool> isInS_, double violation_) :
                depotId1(depotId1_), depotId2(depotId2_), setR(setR_.begin(), setR_.end()), isInS(std::move(isInS_)),
                violation(violation_)
        {}

        DepotCapacityCutCandidate(int depotId1_, int depotId2_, const std::vector<int> & setR_, int nbCustomers,
                                  std::vector<int> setS_) :
                depotId1(depotId1_), depotId2(depotId2_), setR(setR_.begin(), setR_.end()),
                isInS(nbCustomers + 1, false), violation(0.0)
        {
            for (auto custId : setS_)
                isInS[custId] = true;
        }
    };

    struct CompareDepotCapacityCutCandidatesBySubsets
    {
        inline bool operator()(const DepotCapacityCutCandidate & firstCandidate,
                               const DepotCapacityCutCandidate & secondCandidate) const
        {
            if (firstCandidate.depotId1 != secondCandidate.depotId1)
                return firstCandidate.depotId1 < secondCandidate.depotId1;
            if (firstCandidate.depotId2 != secondCandidate.depotId2)
                return firstCandidate.depotId2 < secondCandidate.depotId2;
            if (firstCandidate.setR != secondCandidate.setR)
                return firstCandidate.setR < secondCandidate.setR;
            return firstCandidate.isInS < secondCandidate.isInS;
        }
    };

    struct CompareDepotCapacityCutCandidatesByViolation
    {
        inline bool operator()(const DepotCapacityCutCandidate & firstCandidate,
                               const DepotCapacityCutCandidate & secondCandidate) const
        {
            if (firstCandidate.violation > secondCandidate.violation + 1e-6)
                return true;
            else if (firstCandidate.violation < secondCandidate.violation - 1e-6)
                return false;
            return CompareDepotCapacityCutCandidatesBySubsets()(firstCandidate, secondCandidate);
        }
    };

    typedef std::set<DepotCapacityCutCandidate, CompareDepotCapacityCutCandidatesBySubsets> SetOfUniqueDCCCandidates;

    struct DCCGreedyHeurCandidate
    {
        int custId;
        double lhsChange;
        double violationChange;
        bool wasGeneratedBefore;

        DCCGreedyHeurCandidate(int custId_, double lhsChange_ = 0.0, double violationChange_ = 0.0) :
                custId(custId_), lhsChange(lhsChange_), violationChange(violationChange_), wasGeneratedBefore(false)
        {}
    };

    struct CompDCCGreedyHeurCandidatesByLHSchange
    {
        bool operator()(const DCCGreedyHeurCandidate & firstCand, const DCCGreedyHeurCandidate & secondCand)
        {
            if (firstCand.wasGeneratedBefore != secondCand.wasGeneratedBefore)
                return firstCand.wasGeneratedBefore;
            if (firstCand.lhsChange > secondCand.lhsChange + 1e-6)
                return true;
            if (firstCand.lhsChange < secondCand.lhsChange - 1e-6)
                return false;
            return firstCand.custId < secondCand.custId;
        }
    };

    struct CompDCCGreedyHeurCandidatesByViolationChange
    {
        bool operator()(const DCCGreedyHeurCandidate & firstCand, const DCCGreedyHeurCandidate & secondCand)
        {
            if (firstCand.wasGeneratedBefore != secondCand.wasGeneratedBefore)
                return firstCand.wasGeneratedBefore;
            if (firstCand.violationChange < secondCand.violationChange - 1e-6)
                return true;
            if (firstCand.violationChange > secondCand.violationChange + 1e-6)
                return false;
            if (firstCand.lhsChange > secondCand.lhsChange + 1e-6)
                return true;
            if (firstCand.lhsChange < secondCand.lhsChange - 1e-6)
                return false;
            return firstCand.custId < secondCand.custId;
        }
    };


    class DepotCapacityCutSeparationRoutine : public BcCutSeparationFunctor
    {
        const Data & data;
        DepotCapacityCutSeparationSeparatorParameters params;
        int cutCount;
        int cutRound;
        std::default_random_engine randEngine;
        std::vector<std::vector<int> > setRcandidates;
        int numSetRCandidatesToConsider;
        std::vector<int> elementaryVector;
        std::vector<int> demand;
        int totalDemand;

        void buildThresholds(const std::vector<int> & setR, std::vector<int> & thresholds);
        void buildInitialFlows(const std::vector<int> & setR,
                               const std::vector<std::vector<std::vector<double> > > & xSolution,
                               std::vector<std::vector<double> > & aggXsolution, std::vector<double> & flowToS,
                               std::vector<double> & flowToScomplement, std::vector<bool> & inSetS, int & setSdemand,
                               int & sizeS, bool leaveSetSempty = false) const;
        void insertCustomerInSetS(const std::vector<std::vector<double> > & aggXsolution, int insertCustId,
                                  std::vector<bool> & inSetS, int & setSdemand, std::vector<double> & flowToS,
                                  std::vector<double> & flowToScomplement) const;
        void deleteCustomerFromSetS(const std::vector<std::vector<double> > & aggXsolution, int deleteCustId,
                                    std::vector<bool> & inSetS, int & setSdemand, std::vector<double> & flowToS,
                                    std::vector<double> & flowToScomplement) const;
        void buildSolutionWithGRASP(const std::vector<std::vector<double> > & aggXsolution, int demandThreshold,
                                    std::vector<bool> & inSetS, std::vector<double> & flowToS,
                                    std::vector<double> & flowToScomplement, double & lhs, int & setSdemand,
                                    int & sizeS);
        void localSearchWithSwaps(int maxNbSwaps, const std::vector<std::vector<double> > & aggXsolution,
                                  int demandThreshold, std::vector<bool> & inSetS, std::vector<double> & flowToS,
                                  std::vector<double> & flowToScomplement, double & lhs, int & setSdemand, int & sizeS);
        void localSearchInsertDelete(int maxNbInsertsDeletes, int baseSizeS,
                                     const std::vector<std::vector<double> > & aggXsolution, int demandThreshold,
                                     std::vector<bool> & inSetS, std::vector<double> & flowToS,
                                     std::vector<double> & flowToScomplement, double & lhs, int & setSdemand,
                                     int & sizeS);
        void diversifyWithGRASP(const std::vector<std::vector<double> > & aggXsolution, int demandThreshold,
                                std::vector<bool> & inSetS, std::vector<double> & flowToS,
                                std::vector<double> & flowToScomplement, double & lhs, int & setSdemand, int & sizeS);
        void LiguoriHeuristic(const std::vector<int> & setR, const std::vector<double> & ySolution,
                              const std::vector<std::vector<std::vector<double> > > & xSolution,
                              SetOfUniqueDCCCandidates & candidates);
        void insertCutCandidates(const std::vector<int> & setR, int setRcapacity, const std::vector<bool> & inSetS,
                                 int setSdemand, double lhs, const std::vector<double> & ySolution,
                                 SetOfUniqueDCCCandidates & candidates);
        void greedyConstructionHeuristic(const std::vector<int> & setR, const std::vector<double> & ySolution,
                                         const std::vector<std::vector<std::vector<double> > > & xSolution,
                                         SetOfUniqueDCCCandidates & cutCandidates);
        void printCut(const std::vector<double> & ySolution,
                      const std::vector<std::vector<std::vector<double> > > & xSolution,
                      const DepotCapacityCutCandidate & candidate);
        void createCut(BcConstr & newCut, const DepotCapacityCutCandidate & candidate,
                       BcVarArray & yVar, std::vector<BcVarArray> & xVar);
        double computeViolation(const DepotCapacityCutCandidate & candidate, const std::vector<double> & ySolution,
                                const std::vector<std::vector<std::vector<double> > > & xSolution,
                                bool computeOnlyLHS = false);
        void checkManualCutCandidates(const std::vector<double> & ySolution,
                                      const std::vector<std::vector<std::vector<double> > > & xSolution);

    public:
        DepotCapacityCutSeparationRoutine(const Data & data_, DepotCapacityCutSeparationSeparatorParameters params_);

        ~DepotCapacityCutSeparationRoutine() override;

        int operator()(BcFormulation bcForm, BcSolution &primalSol, double & maxViolation,
                       std::list<BcConstr> & cutList) override;
    };

}

#endif
