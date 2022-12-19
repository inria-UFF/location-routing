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

#include "Parameters.h"
#include "Model.h"
#include "bcModelingLanguageC.hpp"
#include "Data.h"
#include "RCSPSolver.h"
#include "CutSeparation.h"

void lrp::generateAllSubsets(int minNumber, int maxNumber, int minSize, int maxSize,
                        std::vector<std::vector<int> > & subsets)
{
    subsets.clear();

    for (int curSize = minSize; curSize <= maxSize; ++curSize)
    {
        if (curSize > maxNumber - minNumber + 1)
            continue;

        std::vector<int> curSubset(curSize);
        for (int index = 0; index < curSize; ++index)
            curSubset[index] = minNumber + index;

        while (true)
        {
            subsets.push_back(curSubset);

            /// find maximum index strictly smaller than maxNumber - curSize + index + 1
            int baseIndex = curSize - 1;
            while ((baseIndex >= 0) && (curSubset[baseIndex] == maxNumber - curSize + baseIndex + 1))
                baseIndex -= 1;

            if (baseIndex == -1)
                break;

            int baseNumber = curSubset[baseIndex] + 1;

            for (int index = baseIndex; index < curSize; ++index)
                curSubset[index] = baseNumber + index - baseIndex;
        }
    }
}

lrp::Model::Model(const BcInitialisation& bc_init) : BcModel(bc_init, bc_init.instanceFile())
{
    BcObjective objective(*this);
    if (data.integerDistances())
        objective.setMinMaxStatus(BcObjStatus::minInt);
    else
        objective.setMinMaxStatus(BcObjStatus::minFloat);

    double upperBound = params.cutOffValue();
	if (upperBound != std::numeric_limits<double>::infinity())
	{
        if (data.integerDistances())
            objective <= upperBound + 1;
        else
            objective <= upperBound + 0.1;
	}
    if (std::abs(upperBound) < 1e6)
        objective.setArtCostValue(std::abs(upperBound));
    else
        objective.setArtCostValue(1e6);

	std::vector<int> customerIds;
	for (int custId = 1; custId <= data.nbCustomers; ++custId)
	    customerIds.push_back(custId);

	std::vector<int> depotIds;
    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        depotIds.push_back(depotId);

    int totalDemand = 0;
    for (auto custId : customerIds)
        totalDemand += data.customers[custId].demand;

    BcMaster master(*this);

    BcVarArray yVar(master, "Y");

    if (!data.depots_open)
    {
        yVar.type('B');
        yVar.priorityForMasterBranching((params.useSubsetBranching()) ? -1 : 2);
        yVar.priorityForSubproblemBranching(-1);
        yVar.defineIndexName(MultiIndexNames('d'));

        for (auto depotId: depotIds)
            objective += data.depots[depotId].cost * yVar(depotId);

        /// valid inequality to bound from below the number of depots used
        std::vector<int> depot_capacities;
        depot_capacities.reserve(depotIds.size());
        for (auto depotId : depotIds)
            depot_capacities.push_back(data.depots[depotId].capacity);
        std::sort(depot_capacities.rbegin(), depot_capacities.rend()); /// sort in reverse order
        int cumDepotCapacity = 0;
        int minNumDepots = 0;
        for (auto depot_cap : depot_capacities)
        {
            if (cumDepotCapacity >= totalDemand)
                break;
            minNumDepots += 1;
            cumDepotCapacity += depot_cap;
        }
        BcConstrArray minDepotNumConstr(master, "MDN");
        minDepotNumConstr(0) >= minNumDepots;
        for (auto depotId : depotIds)
            minDepotNumConstr[0] += yVar[depotId];
    }

    BcVarArray zVar(master, "Z");
    if (params.useMasterZvars())
    {
        zVar.type('B');
        zVar.priorityForMasterBranching(1);
        zVar.priorityForSubproblemBranching(-1);
        zVar.defineIndexName(MultiIndexNames('d', 'j'));
    }

    /// \sum_{i\in I} z_{ij} == 1 \forall j\in J
    BcConstrArray assignConstr(master, "ASS");
    for (auto custId : customerIds)
    {
        if (params.useMasterZvars())
        {
            assignConstr(custId) == 1;
            for (auto depotId: depotIds)
                assignConstr[custId] += zVar(depotId, custId);
        } else {
            assignConstr(custId) >= 2;
        }
    }

    /// sum_{e\in\delta_j}x_e^i == 2*z_{ij} \forall i\in I, j\in J
    BcConstrArray degreeConstr(master, "DEG");
    if (params.useMasterZvars())
    {
        for (auto depotId: depotIds)
            for (auto custId: customerIds)
            {
                degreeConstr(depotId, custId) >= 0;
                degreeConstr[depotId][custId] += (-2) * zVar[depotId][custId];
            }
    }

    BcConstrArray capacityConstr(master, "CPC");
    for (auto depotId : depotIds)
    {
        capacityConstr(depotId) <= ((data.depots_open) ? data.depots[depotId].capacity : 0);
        if (!data.depots_open)
            capacityConstr[depotId] += (-data.depots[depotId].capacity) * yVar[depotId];
    }

    /// valid inequality to bound from below the number of vehicles used
    int minNumVehicles = (totalDemand - 1) / data.vehicle_cap + 1;
    BcConstrArray minVehNumConstr(master, "MVN");
    minVehNumConstr(0) >= minNumVehicles;

    BcBranchingConstrArray vehNumberBranching(master, "VNB", SelectionStrategy::MostFractional, 1.0);
    vehNumberBranching(0); /// this is for the total number of vehicles
    for (auto depotId : depotIds)
        vehNumberBranching(depotId);

    BcBranchingConstrArray custNumberBranching(master, "CNB", SelectionStrategy::MostFractional, 1.0);
    for (auto depotId : depotIds)
        custNumberBranching(depotId);

    BcBranchingConstrArray edgeBranching(master, "EDGE", SelectionStrategy::MostFractional, 1.0);
    for (int firstCustId = 0; firstCustId <= data.nbCustomers; ++firstCustId )
        for (int secondCustId = firstCustId + 1; secondCustId <= data.nbCustomers; ++secondCustId)
            edgeBranching(firstCustId, secondCustId);

    BcBranchingConstrArray assignBranching(master, "ASB", SelectionStrategy::MostFractional, 1.0);
    if (!params.useMasterZvars())
    {
        for (auto depotId: depotIds)
            for (auto custId: customerIds)
                assignBranching(depotId, custId);
    }

    if (!data.depots_open && params.useSubsetBranching())
    {
        std::vector<std::vector<int> > depotSubsets;
        generateAllSubsets(1, data.nbDepots, 1, 4, depotSubsets);
        BcBranchingConstrArray subsetBranching(master, "DSS", SelectionStrategy::MostFractional, 2.0);
        for (const auto & subset : depotSubsets)
        {
            BcConstr bcConstr(nullptr);
            if (subset.size() == 1)
                bcConstr = subsetBranching(subset[0], 0, 0, 0);
            else if (subset.size() == 2)
                bcConstr = subsetBranching(subset[0], subset[1], 0, 0);
            else if (subset.size() == 3)
                bcConstr = subsetBranching(subset[0], subset[1], subset[2], 0);
            else
                bcConstr = subsetBranching(subset[0], subset[1], subset[2], subset[3]);
            for (auto depotId : subset)
                bcConstr += 1 * yVar[depotId];
        }
    }

    BcColGenSpArray depotCGSp(*this);
    for (auto depotId : depotIds)
    {
        BcFormulation spForm = depotCGSp(depotId);
        spForm <= data.depots[depotId].veh_number;

        BcVarArray xVar(spForm, "X");
        xVar.type('I');
        xVar.priorityForMasterBranching(-1);
        xVar.priorityForSubproblemBranching(-1);
        xVar.defineIndexName(MultiIndexNames('d', 'i', 'j'));

        /// customer number 0 is the depot
        for (int firstCustId = 0; firstCustId <= data.nbCustomers; ++firstCustId )
            for (int secondCustId = firstCustId + 1; secondCustId <= data.nbCustomers; ++secondCustId)
            {
                BcVar bcVar = xVar(depotId, firstCustId, secondCustId);
                if (firstCustId == 0)
                {
                    minVehNumConstr[0] += 0.5 * bcVar;
                    vehNumberBranching[0] += 0.5 * bcVar;
                    vehNumberBranching[depotId] += 0.5 * bcVar;
                    if (data.depots[depotId].veh_cost > 0.0)
                        objective += (0.5 * data.depots[depotId].veh_cost
                                        + data.getDepotToCustDistance(depotId, secondCustId)) * bcVar;
                    else
                        objective += data.getDepotToCustDistance(depotId, secondCustId) * bcVar;
                    custNumberBranching[depotId] += 0.5 * bcVar;
                }
                else
                {
                    if (params.useMasterZvars())
                    {
                        degreeConstr[depotId][firstCustId] += bcVar;
                    }
                    else
                    {
                        assignConstr[firstCustId] += bcVar;
                        assignBranching[depotId][firstCustId] += 0.5 * bcVar;
                    }
                    objective += data.getCustToCustDistance(firstCustId, secondCustId) * bcVar;
                    custNumberBranching[depotId] += bcVar;
                }
                if (params.useMasterZvars())
                {
                    degreeConstr[depotId][secondCustId] += bcVar;
                }
                else
                {
                    assignConstr[secondCustId] += bcVar;
                    assignBranching[depotId][secondCustId] += 0.5 * bcVar;
                }
                edgeBranching[firstCustId][secondCustId] += bcVar;
            }

        BcVarArray rVar(spForm, "R");
        rVar.type('C');
        rVar.priorityForMasterBranching(-1);
        rVar.priorityForSubproblemBranching(-1);
        rVar.defineIndexName(MultiIndexNames('d'));

        capacityConstr[depotId] += 1 * rVar(depotId);

        RCSPSolver solver(spForm, depotId);
        spForm.attach(solver.getOraclePtr());
    }

    std::vector<int> demands(data.nbCustomers + 1, 0);
    for (auto custId : customerIds)
        demands[custId] = data.customers[custId].demand;

    BcCapacityCutConstrArray capacityCuts(master, data.vehicle_cap, demands, true, true,
                                          -1, 3.0, 1.0);

    BcLimMemRankOneCutConstrArray limMemRank1Cuts(master);

    if (params.useRLKC())
        BcResConsumptionKnapsackCutConstrArray resConsKnapCuts(master, 1.0, 1.0);

    if (!data.depots_open && params.useGUB())
    {
        BcCutConstrArray gubCutConstr(master, "GUB", 'F', 3.0, 1.0);
        gubCutConstr.attach(new GeneralizedUpperBoundSeparationRoutine(data, params));
    }

    if (!data.depots_open && params.useCOV())
    {
        BcCutConstrArray covCutConstr(master, "COV", 'F', 3.0, 1.0);
        covCutConstr.attach(new DepotCoverInequalitySeparationRoutine(data, params));
    }

    if (!data.depots_open && params.useFCC())
    {
        BcCutConstrArray fccCutConstr(master, "FCC", 'F', 3.0, 1.0);
        fccCutConstr.attach(new FenchelCutSeparationRoutine(data, params));
    }

    if (params.useDCC())
    {
        BcCutConstrArray dccCutConstr(master, "DCC", 'F', 3.0, 1.0);
        DepotCapacityCutSeparationSeparatorParameters DCCparams;
        DCCparams.maxNumPerRound = params.DCCmaxNumPerRound();
        DCCparams.genCutsWithTwoYs = params.useDCCtwoYs();
        if (params.DCCheuristic() <= 1)
        {
            DCCparams.heuristicChoice = 0;
            DCCparams.useDiversification = params.DCCheuristic();
        }
        else
        {
            DCCparams.heuristicChoice = 1;
            DCCparams.greedyConstructionHeuristicCriterion = (params.DCCheuristic() == 2) ? 0 : 1;
        }
        dccCutConstr.attach(new DepotCapacityCutSeparationRoutine(data, DCCparams));
    }

}
