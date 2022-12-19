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

#include "CutSeparation.h"
#include "InputUser.h"
#include <unordered_set>
#include <bitset>

lrp::DepotCapacityCutSeparationRoutine
     ::DepotCapacityCutSeparationRoutine(const Data & data_, DepotCapacityCutSeparationSeparatorParameters params_):
        data(data_), params(params_), cutCount(0), cutRound(0), randEngine(3), setRcandidates(),
        numSetRCandidatesToConsider(0), elementaryVector(), demand(), totalDemand(0)
{
    //params.printLevel = 0;

    generateAllSubsets(1, data.nbDepots, (std::min)(params.setRcandidateMinSize, data.nbDepots),
                       (std::min)(params.setRcandidateMaxSize, data.nbDepots), setRcandidates);

    elementaryVector = std::vector<int>(setRcandidates.size());
    for (int index = 0; index < (int)elementaryVector.size(); ++index)
        elementaryVector[index] = index;

    /// this part is taken from the Liguori's code
    double candidatesFactor = 7; /// x2 in comparison with Liguori parameterisation
    if (data.nbDepots >= 20)
        candidatesFactor = 3.6; /// x2 in comparison with Liguori parameterisation
    else if (data.nbDepots >= 15)
        candidatesFactor = 5; /// x2 in comparison with Liguori parameterisation

    numSetRCandidatesToConsider = (int)ceil(candidatesFactor * data.nbDepots);

    demand.resize(data.nbCustomers + 1, 0);
    for (int custId = 1; custId <= data.nbCustomers; ++custId)
    {
        demand[custId] = data.customers[custId].demand;
        totalDemand += demand[custId];
    }
}

lrp::DepotCapacityCutSeparationRoutine::~DepotCapacityCutSeparationRoutine() = default;

void lrp::DepotCapacityCutSeparationRoutine::buildThresholds(const std::vector<int> & setR,
                                                             std::vector<int> & thresholds)
{
    /// thresholds are generated in the same way as in the Liguori's code
    /// every threshold is a value for the d(S) such that
    /// r(S,R) or r(S,R\{i_1}) or r(S,R\{i_1,i_2}) changes its value

    thresholds.clear();

    int setRcapacity = 0;
    for (auto depotId : setR)
        setRcapacity += data.depots[depotId].capacity;
    int initSetRCapacity = setRcapacity;

    int threshold = setRcapacity + 1;
    while (threshold < totalDemand)
    {
        thresholds.push_back(threshold);
        threshold += data.vehicle_cap;
    }

    setRcapacity -= data.depots[setR[0]].capacity;
    threshold = setRcapacity + 1;
    while (threshold < totalDemand)
    {
        if (threshold > initSetRCapacity)
            thresholds.push_back(threshold);
        threshold += data.vehicle_cap;
    }

    if ((int)setR.size() > 1)
    {
        setRcapacity -= data.depots[setR[1]].capacity;
        threshold = setRcapacity + 1;
        while (threshold < totalDemand)
        {
            if (threshold > initSetRCapacity)
                thresholds.push_back(threshold);
            threshold += data.vehicle_cap;
        }
    }

    std::stable_sort(thresholds.begin(), thresholds.end());
    auto iter = std::unique(thresholds.begin(), thresholds.end());
    thresholds.resize(std::distance(thresholds.begin(), iter));
}

void lrp::DepotCapacityCutSeparationRoutine
     ::buildInitialFlows(const std::vector<int> & setR,
                         const std::vector<std::vector<std::vector<double> > > & xSolution,
                         std::vector<std::vector<double> > & aggXsolution, std::vector<double> & flowToS,
                         std::vector<double> & flowToScomplement, std::vector<bool> & inSetS, int & setSdemand,
                         int & sizeS, bool leaveSetSempty) const
{
    std::vector<bool> isInSetR(data.nbDepots + 1, false);
    for (auto depotId : setR)
        isInSetR[depotId] = true;

    aggXsolution.resize(data.nbCustomers + 1, std::vector<double>(data.nbCustomers + 1, 0.0));
    flowToS.resize(data.nbCustomers + 1, 0.0);
    flowToScomplement.resize(data.nbCustomers + 1, 0.0);

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        if (isInSetR[depotId])
            continue;

        for (int firstCustId = 0; firstCustId <= data.nbCustomers; ++firstCustId)
            for (int secondCustId = firstCustId + 1; secondCustId <= data.nbCustomers; ++secondCustId)
            {
                auto value = xSolution[depotId][firstCustId][secondCustId];
                aggXsolution[firstCustId][secondCustId] += value;
                flowToScomplement[firstCustId] += value;
                flowToScomplement[secondCustId] += value;
            }
    }

    inSetS.resize(data.nbCustomers + 1, false);

    if (leaveSetSempty)
        return;

    /// initial set S includes all customers served only from depots in setR
    for (int custId = 1; custId <= data.nbCustomers; ++custId)
        if (flowToScomplement[custId] < 1e-6)
        {
            inSetS[custId] = true;
            setSdemand += demand[custId];
            sizeS += 1;
        }
}

void lrp::DepotCapacityCutSeparationRoutine
     ::insertCustomerInSetS(const std::vector<std::vector<double> > & aggXsolution, int insertCustId,
                            std::vector<bool> & inSetS, int & setSdemand, std::vector<double> & flowToS,
                            std::vector<double> & flowToScomplement) const
{
    if (inSetS[insertCustId])
        return;

    inSetS[insertCustId] = true;
    setSdemand += demand[insertCustId];
    for (int custId = 0; custId < insertCustId; ++custId)
    {
        flowToScomplement[custId] -= aggXsolution[custId][insertCustId];
        flowToS[custId] += aggXsolution[custId][insertCustId];
    }
    for (int custId = insertCustId + 1; custId <= data.nbCustomers; ++custId)
    {
        flowToScomplement[custId] -= aggXsolution[insertCustId][custId];
        flowToS[custId] += aggXsolution[insertCustId][custId];
    }
}

void lrp::DepotCapacityCutSeparationRoutine
     ::deleteCustomerFromSetS(const std::vector<std::vector<double> > & aggXsolution, int deleteCustId,
                              std::vector<bool> & inSetS, int & setSdemand, std::vector<double> & flowToS,
                              std::vector<double> & flowToScomplement) const
{
    if (!inSetS[deleteCustId])
        return;

    inSetS[deleteCustId] = false;
    setSdemand -= demand[deleteCustId];
    for (int custId = 0; custId < deleteCustId; ++custId)
    {
        flowToScomplement[custId] += aggXsolution[custId][deleteCustId];
        flowToS[custId] -= aggXsolution[custId][deleteCustId];
    }
    for (int custId = deleteCustId + 1; custId <= data.nbCustomers; ++custId)
    {
        flowToScomplement[custId] += aggXsolution[deleteCustId][custId];
        flowToS[custId] -= aggXsolution[deleteCustId][custId];
    }
}

/// build solution with a constructive GRASP heursitic
void lrp::DepotCapacityCutSeparationRoutine
     ::buildSolutionWithGRASP(const std::vector<std::vector<double> > & aggXsolution, int demandThreshold,
                              std::vector<bool> & inSetS, std::vector<double> & flowToS,
                              std::vector<double> & flowToScomplement, double & lhs, int & setSdemand, int & sizeS)
{
    while (sizeS < data.nbCustomers - 1)
    {
        std::vector<std::pair<double, int> > insertCandidates;
        insertCandidates.reserve(data.nbCustomers - sizeS);
        for (int custId = 1; custId <= data.nbCustomers; ++custId)
            if (!inSetS[custId]) {
                double lhsChange = flowToScomplement[custId] - flowToS[custId];
                insertCandidates.emplace_back(lhsChange, custId);
            }
        std::stable_sort(insertCandidates.begin(), insertCandidates.end());
        double delta =
                (insertCandidates.back().first - insertCandidates.front().first) * params.alphaParamInInitialGRASP;
        double GRASPthreshold = insertCandidates.back().first + delta;
        int candIndex = 1;
        while ((candIndex < (int) insertCandidates.size()) && (insertCandidates[candIndex].first <= GRASPthreshold))
            candIndex += 1;
        std::uniform_int_distribution<int> distribution(0, candIndex - 1);
        int randomIndex = distribution(randEngine);
        double lhsChange = insertCandidates[randomIndex].first;

        if ((lhsChange >= 0) && (setSdemand >= demandThreshold))
            break;

        int insertCustId = insertCandidates[randomIndex].second;
        insertCustomerInSetS(aggXsolution, insertCustId, inSetS, setSdemand, flowToS, flowToScomplement);
        lhs += lhsChange;
        sizeS += 1;
    }
}

/// do the local search where the neighbourhood is to exchange two customers, one in setS and another not in setS
void lrp::DepotCapacityCutSeparationRoutine
     ::localSearchWithSwaps(int maxNbSwaps, const std::vector<std::vector<double> > & aggXsolution, int demandThreshold,
                            std::vector<bool> & inSetS, std::vector<double> & flowToS,
                            std::vector<double> & flowToScomplement, double & lhs, int & setSdemand, int & sizeS)
{
    for (int numSwaps = 1; numSwaps <= maxNbSwaps; ++numSwaps)
    {
        auto bestPair = std::make_pair(0, 0);
        double bestLhsChange = 1e-6; /// we look for a non-positive LHS change
        for (int custIdInS = 1; custIdInS <= data.nbCustomers; ++custIdInS)
        {
            if (!inSetS[custIdInS])
                continue;
            for (int custIdNotInS = 1; custIdNotInS <= data.nbCustomers; ++custIdNotInS)
            {
                if (inSetS[custIdNotInS])
                    continue;
                if (setSdemand - demand[custIdInS] + demand[custIdNotInS] < demandThreshold)
                    continue;
                double lhsChange = flowToS[custIdInS] - flowToScomplement[custIdInS]
                                   + flowToScomplement[custIdNotInS] - flowToS[custIdNotInS];
                lhsChange += (custIdInS < custIdNotInS) ? 2 * aggXsolution[custIdInS][custIdNotInS]
                                                        : 2 * aggXsolution[custIdNotInS][custIdInS];
                if (lhsChange < bestLhsChange)
                {
                    bestLhsChange = lhsChange;
                    bestPair = std::make_pair(custIdInS, custIdNotInS);
                }
            }
        }
        if (bestPair.first == 0)
            break;

        lhs += bestLhsChange;
        deleteCustomerFromSetS(aggXsolution, bestPair.first, inSetS, setSdemand, flowToS, flowToScomplement);
        insertCustomerInSetS(aggXsolution, bestPair.second, inSetS, setSdemand, flowToS, flowToScomplement);
    }
}

/// do the local search where the neighbourhood is to insert a customer to setS or to delete a customer from setS
void lrp::DepotCapacityCutSeparationRoutine
     ::localSearchInsertDelete(int maxNbInsertsDeletes, int baseSizeS,
                               const std::vector<std::vector<double> > & aggXsolution, int demandThreshold,
                               std::vector<bool> & inSetS, std::vector<double> & flowToS,
                               std::vector<double> & flowToScomplement, double & lhs, int & setSdemand, int & sizeS)
{
    for (int numInsertDelete = 1; numInsertDelete <= maxNbInsertsDeletes; ++numInsertDelete)
    {
        int insertCustId = 0;
        double lhsChange = 1e-6; /// we look for a non-positive LHS change
        if (sizeS < data.nbCustomers - 1)
            for (int custId = 1; custId <= data.nbCustomers; ++custId)
                if (!inSetS[custId] && (flowToScomplement[custId] - flowToS[custId] < lhsChange))
                {
                    lhsChange = flowToScomplement[custId] - flowToS[custId];
                    insertCustId = custId;
                }
        if (insertCustId >= 1)
        {
            lhs += lhsChange;
            sizeS += 1;
            insertCustomerInSetS(aggXsolution, insertCustId, inSetS, setSdemand, flowToS, flowToScomplement);
        }

        int deleteCustId = 0;
        lhsChange = 1e-6; /// we look for a non-positive LHS change
        if (sizeS > baseSizeS)
            for (int custId = 1; custId <= data.nbCustomers; ++custId)
                if (inSetS[custId] && (flowToS[custId] - flowToScomplement[custId] < lhsChange)
                    && (setSdemand - demand[custId] >= demandThreshold))
                {
                    lhsChange = flowToS[custId] - flowToScomplement[custId];
                    deleteCustId = custId;
                }
        if (deleteCustId >= 1)
        {
            lhs += lhsChange;
            sizeS -= 1;
            deleteCustomerFromSetS(aggXsolution, deleteCustId, inSetS, setSdemand, flowToS, flowToScomplement);
        }
        if ((insertCustId == 0) && (deleteCustId == 0))
            break;
    }
}

void lrp::DepotCapacityCutSeparationRoutine
     ::diversifyWithGRASP(const std::vector<std::vector<double> > & aggXsolution,
                          int demandThreshold, std::vector<bool> & inSetS, std::vector<double> & flowToS,
                          std::vector<double> & flowToScomplement, double & lhs, int & setSdemand, int & sizeS)
{
    auto maxNbSwaps = (int)ceil(data.nbCustomers * params.nbSwapsFactorInDiversificationGRASP);
    std::vector<bool> isInTabuList(data.nbCustomers + 1, false);
    for (int numSwaps = 1; numSwaps <= maxNbSwaps; ++numSwaps)
    {
        std::vector<std::tuple<double, int, int> > swapCandidates;
        swapCandidates.reserve((data.nbCustomers - sizeS) );
        for (int custIdInS = 1; custIdInS <= data.nbCustomers; ++custIdInS)
        {
            if (!inSetS[custIdInS] || isInTabuList[custIdInS])
                continue;
            for (int custIdNotInS = 1; custIdNotInS <= data.nbCustomers; ++custIdNotInS)
            {
                if (inSetS[custIdNotInS] || isInTabuList[custIdNotInS])
                    continue;
                double lhsChange = flowToS[custIdInS] - flowToScomplement[custIdInS]
                                   + flowToScomplement[custIdNotInS] - flowToS[custIdNotInS];
                lhsChange += (custIdInS < custIdNotInS) ? 2 * aggXsolution[custIdInS][custIdNotInS]
                                                        : 2 * aggXsolution[custIdNotInS][custIdInS];
                swapCandidates.emplace_back(lhsChange, custIdInS, custIdNotInS);
            }
        }
        if (!swapCandidates.empty())
        {
            std::stable_sort(swapCandidates.begin(), swapCandidates.end());
            double delta = (std::get<0>(swapCandidates.back()) - std::get<0>(swapCandidates.front()))
                           * params.alphaParamInDiversificationGRASP;
            double GRASPthreshold = std::get<0>(swapCandidates.back()) + delta;
            int candIndex = 1;
            while ((candIndex < (int) swapCandidates.size()) &&
                   (std::get<0>(swapCandidates[candIndex]) <= GRASPthreshold))
                candIndex += 1;
            std::uniform_int_distribution<int> distribution(0, candIndex - 1);
            int randomIndex = distribution(randEngine);

            lhs += std::get<0>(swapCandidates[randomIndex]);
            int custIdInS = std::get<1>(swapCandidates[randomIndex]);
            int custIdNotInS = std::get<2>(swapCandidates[randomIndex]);
            deleteCustomerFromSetS(aggXsolution, custIdInS, inSetS,
                                   setSdemand, flowToS, flowToScomplement);
            insertCustomerInSetS(aggXsolution, custIdNotInS, inSetS,
                                 setSdemand, flowToS, flowToScomplement);
            isInTabuList[custIdNotInS] = isInTabuList[custIdInS] = true;
        }
        else
        {
            break;
        }
    }

    /// complete the setS greedily until the total demain in setS reaches demandThreshold
    while (setSdemand < demandThreshold)
    {
        int insertCustId = 0;
        double lhsChange = 2.1; /// change cannot be larger than 2
        if (sizeS < data.nbCustomers - 1)
            for (int custId = 1; custId <= data.nbCustomers; ++custId)
                if (!inSetS[custId] && (flowToScomplement[custId] - flowToS[custId] < lhsChange))
                {
                    lhsChange = flowToScomplement[custId] - flowToS[custId];
                    insertCustId = custId;
                }
        if (insertCustId >= 1)
        {
            lhs += lhsChange;
            sizeS += 1;
            insertCustomerInSetS(aggXsolution, insertCustId, inSetS, setSdemand, flowToS, flowToScomplement);
        }
        else
        {
            break;
        }
    }

    /// insert randomly another customer
    if (sizeS < data.nbCustomers - 1)
    {
        std::uniform_int_distribution<int> distribution(0, data.nbCustomers - sizeS - 1);
        int randomOrder = distribution(randEngine);
        int order = 0;
        int insertCustId = 0;
        for (int custId = 1; (custId <= data.nbCustomers) && (insertCustId == 0); ++custId)
            if (!inSetS[custId])
            {
                if (order == randomOrder)
                    insertCustId = custId;
                else
                    order += 1;
            }
        if (insertCustId == 0)
        {
            lhs += (flowToScomplement[insertCustId] - flowToS[insertCustId]);
            sizeS += 1;
            insertCustomerInSetS(aggXsolution, insertCustId, inSetS, setSdemand, flowToS, flowToScomplement);
        }
    }
}

void lrp::DepotCapacityCutSeparationRoutine
     ::insertCutCandidates(const std::vector<int> & setR, int setRcapacity, const std::vector<bool> & inSetS,
                           int setSdemand, double lhs, const std::vector<double> & ySolution,
                           SetOfUniqueDCCCandidates & candidates)
{
    /// we first find the best RHS for the Belenguer et al. cut (with $y_{i_1}$ only)
    int numVehSR = (setSdemand - setRcapacity - 1) / data.vehicle_cap + 1;
    double bestRHS = 0.0;
    int bestDepotId = 0;
    for (auto depotId : setR)
    {
        int numVehSRi1 = (setSdemand - setRcapacity + data.depots[depotId].capacity - 1) / data.vehicle_cap + 1;
        double rhs = 2 * numVehSR + ((data.depots_open) ? 0.0 : 2 * (1 - ySolution[depotId]) * (numVehSRi1 - numVehSR));
        if (rhs > bestRHS)
        {
            bestDepotId = depotId;
            bestRHS = rhs;
        }
    }
    if ((bestDepotId >= 1) && (lhs < bestRHS - params.cutViolationTolerance))
        candidates.insert(DepotCapacityCutCandidate(bestDepotId, 0, setR, inSetS,
                                                    bestRHS - lhs));

    /// we now find the best RHS for the our cut (with $y_{i_1}$ and $y_{i_2}$)
    if ((setR.size() > 1) && params.genCutsWithTwoYs && !data.depots_open)
    {
        bestRHS = 0.0;
        auto bestDepotPair = std::make_pair(0, 0);
        for (auto depotId1 : setR)
            for (auto depotId2 : setR)
                if (depotId1 != depotId2)
                {
                    int numVehSRi1 = (setSdemand - setRcapacity + data.depots[depotId1].capacity - 1)
                                     / data.vehicle_cap + 1;
                    int numVehSRi1i2 = (setSdemand - setRcapacity + data.depots[depotId1].capacity
                                        + data.depots[depotId2].capacity - 1) / data.vehicle_cap + 1;
                    double rhs = 2 * numVehSR * ySolution[depotId1] + 2 * numVehSRi1 * ySolution[depotId2]
                                 + 2 * numVehSRi1i2 * (1 - ySolution[depotId1] - ySolution[depotId2]);
                    if (rhs > bestRHS)
                    {
                        bestDepotPair = std::make_pair(depotId1, depotId2);
                        bestRHS = rhs;
                    }
                }
        if ((bestDepotPair.first >= 1) && (lhs < bestRHS - params.cutViolationTolerance))
            candidates.insert(DepotCapacityCutCandidate(bestDepotPair.first, bestDepotPair.second, setR,
                                                        inSetS, bestRHS - lhs));
    }
}

void lrp::DepotCapacityCutSeparationRoutine
     ::LiguoriHeuristic(const std::vector<int> & setR, const std::vector<double> & ySolution,
                        const std::vector<std::vector<std::vector<double> > > & xSolution,
                        SetOfUniqueDCCCandidates & candidates)
{
    int setRcapacity = 0;
    for (auto depotId : setR)
        setRcapacity += data.depots[depotId].capacity;

    std::vector<std::vector<double> > aggXsolution;
    std::vector<double> baseFlowToS;
    std::vector<double> baseFlowToScomplement;
    std::vector<bool> baseInSetS;
    int baseSetSdemand = 0;
    int baseSizeS = 0;
    buildInitialFlows(setR, xSolution, aggXsolution, baseFlowToS, baseFlowToScomplement, baseInSetS,
                      baseSetSdemand, baseSizeS);

    std::vector<int> thresholds;
    buildThresholds(setR, thresholds);

    for (auto threshold: thresholds)
    {
        std::vector<double> flowToS(baseFlowToS);
        std::vector<double> flowToScomplement(baseFlowToScomplement);
        std::vector<bool> inSetS(baseInSetS);
        double lhs = 0.0;
        int setSdemand = baseSetSdemand;
        int sizeS = baseSizeS;

        buildSolutionWithGRASP(aggXsolution, threshold, inSetS, flowToS, flowToScomplement, lhs,
                               setSdemand, sizeS);

        localSearchWithSwaps(params.maxNbSwapsInInitialGRASP, aggXsolution, threshold, inSetS, flowToS,
                             flowToScomplement, lhs, setSdemand, sizeS);

        localSearchInsertDelete(params.maxInsertDeleteIterations, baseSizeS, aggXsolution, threshold, inSetS,
                                flowToS, flowToScomplement, lhs, setSdemand, sizeS);

        if (setSdemand >= threshold)
            insertCutCandidates(setR, setRcapacity, inSetS, setSdemand, lhs, ySolution, candidates);

        if (params.useDiversification)
        {
            diversifyWithGRASP(aggXsolution, threshold, inSetS, flowToS, flowToScomplement, lhs,
                               setSdemand, sizeS);

            localSearchWithSwaps(params.maxNbSwapsAfterDiversification, aggXsolution, threshold, inSetS, flowToS,
                                 flowToScomplement, lhs, setSdemand, sizeS);

            localSearchInsertDelete(params.maxInsertDeleteIterations, baseSizeS, aggXsolution, threshold, inSetS,
                                    flowToS, flowToScomplement, lhs, setSdemand, sizeS);

            if (setSdemand >= threshold)
                insertCutCandidates(setR, setRcapacity, inSetS, setSdemand, lhs, ySolution, candidates);
        }
    } /// for (auto threshold : thresholds)
}

double computeRhs(int setSdemand, int setRcapacity, int setRi1capacity, int vehicle_cap, double yi1Value)
{
    int numVehSR = (setSdemand <= setRcapacity) ? 0 : (setSdemand - setRcapacity - 1) / vehicle_cap + 1;
    int numVehSRi1 = (setSdemand <= setRi1capacity) ? 0 : (setSdemand - setRi1capacity - 1) / vehicle_cap + 1;
    return 2 * numVehSR + 2 * (1 - yi1Value) * (numVehSRi1 - numVehSR);
}

void lrp::DepotCapacityCutSeparationRoutine
     ::greedyConstructionHeuristic(const std::vector<int> & setR, const std::vector<double> & ySolution,
                                   const std::vector<std::vector<std::vector<double> > > & xSolution,
                                   SetOfUniqueDCCCandidates & cutCandidates)
{
    if (printL(params.printLevel))
    {
        std::cout << "Launching greedy construction heuristic for set R =";
        for (auto depotId: setR)
            std::cout << " " << depotId;
        std::cout << std::endl;
    }
    int setRcapacity = 0;
    for (auto depotId : setR)
        setRcapacity += data.depots[depotId].capacity;
    int setRi1capacity = setRcapacity - data.depots[setR[0]].capacity;
    double yi1Value = ySolution[setR[0]];

    std::vector<std::vector<double> > aggXsolution;
    std::vector<double> baseFlowToS;
    std::vector<double> baseFlowToScomplement;
    std::vector<bool> baseInSetS;
    int baseSetSdemand = 0;
    int baseSizeS = 0;
    buildInitialFlows(setR, xSolution, aggXsolution, baseFlowToS, baseFlowToScomplement, baseInSetS,
                      baseSetSdemand, baseSizeS);

    /// we create adjacency lists for the fractional graph
    std::vector<std::vector<int> > adjacencyList(data.nbCustomers + 1);
    for (int firstCustId = 1; firstCustId <= data.nbCustomers; ++firstCustId)
    {
        for (int secondCustId = firstCustId + 1; secondCustId <= data.nbCustomers; ++secondCustId)
        {
            if (aggXsolution[firstCustId][secondCustId] > 1e-6)
            {
                adjacencyList[firstCustId].push_back(secondCustId);
                adjacencyList[secondCustId].push_back(firstCustId);
            }
        }
    }

    std::unordered_set<std::bitset<LRP_MAX_NUMBER_OF_CUSTOMERS> > bitsetSubsetsPool;

    for (int seedCustId = 1; seedCustId <= data.nbCustomers; ++seedCustId)
    {
        if (baseInSetS[seedCustId])
            continue;

        std::bitset<LRP_MAX_NUMBER_OF_CUSTOMERS> bitsetSubset;
        std::vector<double> flowToS(baseFlowToS);
        std::vector<double> flowToScomplement(baseFlowToScomplement);
        std::vector<bool> inSetS(baseInSetS);
        std::vector<bool> inSetSorCandidate(baseInSetS);
        int setSdemand = baseSetSdemand;
        int sizeS = baseSizeS;

        double lhs = flowToScomplement[seedCustId] - flowToS[seedCustId];
        bitsetSubset.set(seedCustId);
        sizeS += 1;
        inSetSorCandidate[seedCustId] = true;
        insertCustomerInSetS(aggXsolution, seedCustId, inSetS, setSdemand, flowToS, flowToScomplement);
        double rhs = computeRhs(setSdemand, setRcapacity, setRi1capacity, data.vehicle_cap, yi1Value);

        std::vector<DCCGreedyHeurCandidate> greedyCandidates;
        for (auto adjCustId : adjacencyList[seedCustId])
            if (!inSetSorCandidate[adjCustId])
            {
                double rhsChange = computeRhs(setSdemand + demand[adjCustId], setRcapacity, setRi1capacity,
                                              data.vehicle_cap, yi1Value) - rhs;
                double lhsChange = flowToScomplement[adjCustId] - flowToS[adjCustId];
                greedyCandidates.emplace_back(adjCustId, lhsChange, rhsChange - lhsChange);
                inSetSorCandidate[adjCustId] = true;
                bitsetSubset.set(adjCustId);
                if (bitsetSubsetsPool.find( bitsetSubset ) != bitsetSubsetsPool.end())
                    greedyCandidates.back().wasGeneratedBefore = true;
                bitsetSubset.reset(adjCustId);
            }

        while (!greedyCandidates.empty())
        {
            /// we determine the best candidate to add to the current subset : it is the last pair after sorting
            if (params.greedyConstructionHeuristicCriterion == 0)
                std::sort(greedyCandidates.begin(), greedyCandidates.end(), CompDCCGreedyHeurCandidatesByLHSchange());
            else
                std::sort(greedyCandidates.begin(), greedyCandidates.end(), CompDCCGreedyHeurCandidatesByViolationChange());
            auto bestCandidate = greedyCandidates.back();
            auto insertCustId = bestCandidate.custId;
            greedyCandidates.pop_back();

            if (bestCandidate.wasGeneratedBefore)
                break;

            for (auto adjCustId : adjacencyList[insertCustId])
                if (!inSetSorCandidate[adjCustId])
                {
                    greedyCandidates.emplace_back(adjCustId, 0.0, 0.0);
                    inSetSorCandidate[adjCustId] = true;
                }

            bitsetSubset.set(insertCustId);
            bitsetSubsetsPool.insert(bitsetSubset);
            sizeS += 1;
            lhs += flowToScomplement[insertCustId] - flowToS[insertCustId];
            insertCustomerInSetS(aggXsolution, insertCustId, inSetS, setSdemand, flowToS,
                                 flowToScomplement);
            rhs = computeRhs(setSdemand, setRcapacity, setRi1capacity, data.vehicle_cap, yi1Value);


            if (setSdemand > setRcapacity)
                insertCutCandidates(setR, setRcapacity, inSetS, setSdemand, lhs, ySolution, cutCandidates);

            bool allCandidatesWereGeneratedBefore = true;
            for (auto & candidate : greedyCandidates)
            {
                double rhsChange = computeRhs(setSdemand + demand[candidate.custId], setRcapacity, setRi1capacity,
                                              data.vehicle_cap, yi1Value) - rhs;
                candidate.lhsChange = flowToScomplement[candidate.custId] - flowToS[candidate.custId];
                candidate.violationChange = rhsChange - candidate.lhsChange;
                /// we check whether the same subset has been encountered before during the greedy construction heuristic
                bitsetSubset.set(candidate.custId);
                candidate.wasGeneratedBefore = (bitsetSubsetsPool.find( bitsetSubset ) != bitsetSubsetsPool.end());
                if (!candidate.wasGeneratedBefore)
                    allCandidatesWereGeneratedBefore = false;
                bitsetSubset.reset(candidate.custId);
            }
            if (allCandidatesWereGeneratedBefore && (sizeS + greedyCandidates.size() < data.nbCustomers))
            {
                for (int custId = 1; custId <= data.nbCustomers; ++custId)
                    if (!inSetSorCandidate[custId])
                    {
                        double rhsChange = computeRhs(setSdemand + demand[custId], setRcapacity, setRi1capacity,
                                                      data.vehicle_cap, yi1Value) - rhs;
                        double lhsChange =  flowToScomplement[custId] - flowToS[custId];
                        greedyCandidates.emplace_back(custId, lhsChange, rhsChange - lhsChange);
                        inSetSorCandidate[custId] = true;
                        bitsetSubset.set(custId);
                        if (bitsetSubsetsPool.find( bitsetSubset ) != bitsetSubsetsPool.end())
                            greedyCandidates.back().wasGeneratedBefore = true;
                        bitsetSubset.reset(custId);
                    }
            }
        }
    }
}

void lrp::DepotCapacityCutSeparationRoutine::printCut(const std::vector<double> & ySolution,
                                                      const std::vector<std::vector<std::vector<double> > > & xSolution,
                                                      const DepotCapacityCutCandidate & candidate)
{
    std::cout << "New DCC cut : setR = [";
    int setRcapacity = 0;
    for (auto depotId : candidate.setR)
    {
        setRcapacity += data.depots[depotId].capacity;
        std::cout << " " << depotId << "(cap=" << data.depots[depotId].capacity << " val=" << ySolution[depotId] << ")";
    }
    std::cout << "] tot.cap. = " << setRcapacity << " i1 = " << candidate.depotId1;
    if (candidate.depotId2 >= 1)
        std::cout << " i2 = " << candidate.depotId2;

    std::cout << " setS = [";

    int setSdemand = 0;
    for (int custId = 1; custId <= data.nbCustomers; ++custId)
        if (candidate.isInS[custId])
        {
            setSdemand += demand[custId];
            std::cout << " " << custId ;
        }
    std::cout << "] tot.dem. = " << setSdemand << " violation " << candidate.violation / 2.0
              << " computed violation = " << computeViolation(candidate, ySolution, xSolution) / 2.0 << std::endl;
}

void lrp::DepotCapacityCutSeparationRoutine::createCut(BcConstr & newCut, const DepotCapacityCutCandidate & candidate,
                                                       BcVarArray & yVar, std::vector<BcVarArray> & xVar)
{
    int setRcapacity = 0;
    std::vector<int> isInR(data.nbDepots + 1, false);
    for (auto depotId : candidate.setR)
    {
        setRcapacity += data.depots[depotId].capacity;
        isInR[depotId] = true;
    }

    std::vector<int> setS;
    std::vector<int> setScomplement;
    setS.reserve(data.nbCustomers);
    setScomplement.reserve(data.nbCustomers);

    int setSdemand = 0;
    for (int custId = 1; custId <= data.nbCustomers; ++custId)
        if (candidate.isInS[custId])
        {
            setS.push_back(custId);
            setSdemand += demand[custId];
        }
        else
        {
            setScomplement.push_back(custId);
        }

    int numVehSR = (setSdemand - setRcapacity - 1) / data.vehicle_cap + 1;
    int numVehSRi1 = (setSdemand - setRcapacity + data.depots[candidate.depotId1].capacity - 1) / data.vehicle_cap + 1;

    if (data.depots_open)
    {
        newCut >= numVehSR;
    }
    else if (candidate.depotId2 == 0)
    {
        newCut >= numVehSRi1;
        newCut += (numVehSRi1 - numVehSR) * yVar[candidate.depotId1];
    }
    else
    {
        int i1i2capacity = data.depots[candidate.depotId1].capacity + data.depots[candidate.depotId2].capacity;
        int numVehSRi1i2 = (setSdemand - setRcapacity + i1i2capacity - 1) / data.vehicle_cap + 1;
        newCut >= numVehSRi1i2;
        newCut += (numVehSRi1i2 - numVehSR) * yVar[candidate.depotId1];
        newCut += (numVehSRi1i2 - numVehSRi1) * yVar[candidate.depotId2];
    }

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        if (isInR[depotId])
            continue;

        for (int custIdInS = 1; custIdInS <= data.nbCustomers; ++custIdInS)
        {
            if (!candidate.isInS[custIdInS])
                continue;
            for (int custIdNotInS = 0; custIdNotInS <= data.nbCustomers; ++custIdNotInS)
            {
                if (candidate.isInS[custIdNotInS])
                    continue;
                if (custIdInS < custIdNotInS)
                    newCut += 0.5 * xVar[depotId][depotId][custIdInS][custIdNotInS];
                else
                    newCut += 0.5 * xVar[depotId][depotId][custIdNotInS][custIdInS];
            }
        }
    }
}

double lrp::DepotCapacityCutSeparationRoutine
       ::computeViolation(const DepotCapacityCutCandidate & candidate, const std::vector<double> & ySolution,
                          const std::vector<std::vector<std::vector<double> > > & xSolution, bool computeOnlyLHS )
{
    std::vector<std::vector<double> > aggXsolution;
    std::vector<double> flowToS;
    std::vector<double> flowToScomplement;
    std::vector<bool> inSetS;
    int setSdemand = 0;
    int sizeS = 0;
    std::vector<int> setR(candidate.setR.begin(), candidate.setR.end());
    buildInitialFlows(setR, xSolution, aggXsolution, flowToS, flowToScomplement, inSetS, setSdemand,
                      sizeS, true);
    double lhs = 0.0;

    for (int custId = 1; custId <= data.nbCustomers; ++custId)
        if (candidate.isInS[custId]) {
            insertCustomerInSetS(aggXsolution, custId, inSetS, setSdemand, flowToS, flowToScomplement);
            lhs += flowToScomplement[custId] - flowToS[custId];
            sizeS += 1;
        }

    if (computeOnlyLHS)
        return lhs;

    int setRcapacity = 0;
    for (auto depotId: setR)
        setRcapacity += data.depots[depotId].capacity;

    double rhs;
    int numVehSR = (setSdemand - setRcapacity - 1) / data.vehicle_cap + 1;
    int numVehSRi1 = (setSdemand - setRcapacity + data.depots[candidate.depotId1].capacity - 1) / data.vehicle_cap + 1;

    if (data.depots_open)
    {
        rhs = 2 * numVehSR;
    }
    else if (candidate.depotId2 == 0)
    {
        rhs = 2 * numVehSR + 2 * (1 - ySolution[candidate.depotId1]) * (numVehSRi1 - numVehSR);
    }
    else
    {
        int numVehSRi1i2 = (setSdemand - setRcapacity + data.depots[candidate.depotId1].capacity
                            + data.depots[candidate.depotId2].capacity - 1) / data.vehicle_cap + 1;
        rhs = 2 * numVehSR * ySolution[candidate.depotId1] + 2 * numVehSRi1 * ySolution[candidate.depotId2]
              + 2 * numVehSRi1i2 * (1 - ySolution[candidate.depotId1] - ySolution[candidate.depotId2]);
    }
    return rhs - lhs;
}

/// this function is for debugging (to compare with the obsolete implementation)
void lrp::DepotCapacityCutSeparationRoutine
     ::checkManualCutCandidates(const std::vector<double> & ySolution,
                                const std::vector<std::vector<std::vector<double> > > & xSolution)
{
    // 50-5-2
    //cut id =31, i1 = 51, i2 = 54, customers : 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
    //cut id =10, i1 = 55, i2 = 0, customers : 1 2 3 4 7 8 9 10 11 12 13 14 15 16 17 18 19 20 24 25 37 44
    //cut id =26, i1 = 51, i2 = 0, customers : 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
    //cut id =29, i1 = 51, i2 = 55, customers : 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50

//    std::vector<DepotCapacityCutCandidate> candidatesToCheck = {
//            DepotCapacityCutCandidate(1, 4, {1, 4}, 50,
//                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
//                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
//                                       41, 42, 43, 44, 45, 46, 47, 48, 49, 50}),
//            DepotCapacityCutCandidate(1, 5, {1, 5}, 50,
//                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
//                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
//                                       41, 42, 43, 44, 45, 46, 47, 48, 49, 50}),
//            DepotCapacityCutCandidate(1, 0, {1}, 50,
//                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
//                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
//                                       41, 42, 43, 44, 45, 46, 47, 48, 49, 50}),
//            DepotCapacityCutCandidate(5, 0, {5}, 50,
//                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
//                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
//                                       41, 42, 43, 44, 45, 46, 47, 48, 49, 50})
//    };

    // 100-10-1

    std::vector<DepotCapacityCutCandidate> candidatesToCheck = {
            DepotCapacityCutCandidate(1, 10, {1, 10}, 100, {1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 64, 65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100}),
            DepotCapacityCutCandidate(6, 10, {6, 10}, 100, {1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 64, 65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100}),
            DepotCapacityCutCandidate(3, 10, {3, 10}, 100, {1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 39, 40, 41, 43, 44, 45, 46, 48, 50, 51, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 68, 69, 70, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100}),
            DepotCapacityCutCandidate(3, 0, {3}, 100, {1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 39, 40, 41, 43, 44, 45, 46, 48, 50, 51, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 68, 69, 70, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100}),
            DepotCapacityCutCandidate(6, 0, {6}, 100, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100}),
            DepotCapacityCutCandidate(7, 5, {5, 7}, 100, {1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 69, 71, 72, 73, 74, 77, 78, 79, 80, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 95, 96, 97, 99, 100}),
            DepotCapacityCutCandidate(7, 0, {7}, 100, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100}),
            DepotCapacityCutCandidate(5, 0, {5}, 100, {2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100}),
            DepotCapacityCutCandidate(3, 0, {3}, 100, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100}),
            DepotCapacityCutCandidate(9, 0, {9}, 100, {1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 19, 22, 23, 24, 25, 26, 27, 30, 31, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 59, 60, 61, 62, 63, 64, 65, 66, 67, 69, 71, 72, 73, 77, 79, 80, 82, 83, 84, 86, 87, 88, 89, 90, 91, 92, 93, 95, 96, 99, 100}),
            DepotCapacityCutCandidate(1, 0, {1}, 100, {1, 2, 3, 4, 5, 8, 10, 12, 14, 16, 19, 23, 30, 31, 35, 36, 38, 39, 42, 44, 47, 48, 49, 50, 51, 52, 54, 60, 62, 63, 66, 67, 69, 72, 73, 77, 80, 82, 84, 86, 87, 88, 89}),
            DepotCapacityCutCandidate(2, 0, {2}, 100, {1, 2, 3, 4, 5, 8, 10, 12, 14, 16, 19, 23, 30, 31, 35, 36, 38, 39, 42, 44, 47, 48, 49, 50, 51, 52, 54, 60, 62, 63, 66, 67, 69, 72, 73, 77, 80, 82, 84, 86, 87, 88, 89})
    };


    for (auto & candidate : candidatesToCheck)
    {
        candidate.violation = computeViolation(candidate, ySolution, xSolution);
        printCut(ySolution, xSolution, candidate);
    }
}

int lrp::DepotCapacityCutSeparationRoutine::operator()(BcFormulation bcForm,
                                                       BcSolution & primalSol,
                                                       double & maxViolation,
                                                       std::list< BcConstr > & cutList)
{
    if (printL(params.printLevel))
        std::cout << "Starting DCC separation" << std::endl;

    std::vector<double> ySolution(data.nbDepots + 1, 0.0);
    if (!data.depots_open)
    {
        std::set<BcVar> yVarSet;
        primalSol.extractVar("Y", yVarSet);
        for (const auto &bcVar: yVarSet) {
            int depotId = bcVar.id().first();
            ySolution[depotId] = floor(bcVar.solVal() * 1000000) / 1000000;
        }
    }

    std::vector<std::vector<std::vector<double> > >
            xSolution(data.nbDepots + 1, std::vector<std::vector<double> >(data.nbCustomers + 1,
                                                                           std::vector<double>(data.nbCustomers + 1,
                                                                                               0.0)));
    std::set< BcVar > xVarSet;
    primalSol.extractVar("X", xVarSet);
    for (const auto & bcVar : xVarSet)
    {
        int depotId = bcVar.id().first();
        int firstCustId = bcVar.id().second();
        int secondCustId = bcVar.id().third();
        xSolution[depotId][firstCustId][secondCustId] += floor(bcVar.solVal() * 1000000) / 1000000;
        if (data.depots_open)
            ySolution[depotId] = 1.0;
    }

    cutRound++;

    SetOfUniqueDCCCandidates uniqueCandidates;

    std::vector<int> setRcandidatesIndices(elementaryVector);
    std::shuffle(setRcandidatesIndices.begin(), setRcandidatesIndices.end(), randEngine);
    if (printL(params.printLevel))
    {
        std::cout << "setRcandidatesIndices = {";
        for (auto index : setRcandidatesIndices)
        {
            std::cout << ((index != setRcandidatesIndices[0]) ? ", " : " ") << index;
        }
        std::cout << "}" << std::endl;
    }
    int numSetRcandidatesConsidered = 0;
    for (auto candIndex : setRcandidatesIndices)
    {
        std::vector<int> setR(setRcandidates[candIndex]);
        if (params.setRselectionMethod == 1)
        {
            bool allYvalsArePositive = true;
            for (auto depotId : setR)
                if (ySolution[depotId] < 1e-6 )
                {
                    allYvalsArePositive = false;
                    break;
                }
            if (!allYvalsArePositive)
                continue;
        }
        numSetRcandidatesConsidered += 1;
        if (numSetRcandidatesConsidered > numSetRCandidatesToConsider)
            break;

        /// we sort depots in setR by ratio yValue/capacity, smaller is this value, more is the chance of
        /// cut violation when the depot is used as i_1 and i_2 in the cut definition in the paper
        /// after sorting, i_1 = setR[0] and i_2 = setR[1]
        if ((int)setR.size() > 1)
            std::stable_sort(setR.begin(), setR.end(),
                             [&ySolution, this](int depot1Id, int depot2Id){
                                 return ySolution[depot1Id]/data.depots[depot1Id].capacity
                                        < ySolution[depot2Id]/data.depots[depot2Id].capacity;});

        if (params.heuristicChoice == 0)
            LiguoriHeuristic(setR, ySolution, xSolution, uniqueCandidates);
        else /// params.heuristicChoice == 1
            greedyConstructionHeuristic(setR, ySolution, xSolution, uniqueCandidates);
    } /// for (auto candIndex : setRcandidatesIndices)
    if (printL(params.printLevel))
        std::cout << "numSetRcandidatesConsidered = " << numSetRcandidatesConsidered << std::endl;

    std::vector<DepotCapacityCutCandidate> cutCandidates;
    cutCandidates.reserve(uniqueCandidates.size());
    for (const auto & candidate : uniqueCandidates)
        cutCandidates.push_back(candidate);
    std::stable_sort(cutCandidates.begin(), cutCandidates.end(), CompareDepotCapacityCutCandidatesByViolation());

    BcCutConstrArray dccCutConstr(bcForm, "DCC");
    BcVarArray yVar(bcForm, "Y");

    std::vector<BcVarArray> xVar(data.nbDepots + 1, BcVarArray());
    for (auto & spForm : bcForm.colGenSubProblemList())
        xVar[spForm.id().first()] = BcVarArray(spForm, "X");

    int numCutsToGenerate = (std::min)((int)cutCandidates.size(), params.maxNumPerRound);
    for (int candIndex = 0; candIndex < numCutsToGenerate; ++candIndex)
    {
        const auto & candidate = cutCandidates[candIndex];
        if (printL(params.printLevel))
            printCut(ySolution, xSolution, candidate);

        BcConstr newCut = dccCutConstr(cutCount++);
        createCut(newCut, candidate, yVar, xVar);
        cutList.push_back(newCut);
    }

    //checkManualCutCandidates(ySolution, xSolution);

    return numCutsToGenerate;
}
