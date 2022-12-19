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


#include "RCSPSolver.h"
#include "Data.h"

lrp::RCSPSolver::RCSPSolver(BcFormulation spForm, int depotId) : spForm(std::move(spForm)), depotId(depotId)

{
	BcNetwork network(spForm, PathAttrMask, data.nbCustomers + 1,
	                  data.nbCustomers + 1);
	BcNetworkResource cap_res(network, 0);
	cap_res.setAsMainResource();

    BcVarArray rVar(spForm, "R");
    cap_res.setAssociatedVar(rVar[depotId]);

    BcNetworkResource time_res;

    if (data.type == Data::VRPTWS)
        time_res = BcNetworkResource(network, 1);

    buildVertices(depotId, network, cap_res, time_res);
	buildArcs(depotId, network, cap_res, time_res);
    buildElemSetDistanceMatrix(network);

    oraclePtr = new BcRCSPFunctor(spForm);
}

void lrp::RCSPSolver::buildVertices(int depotId, BcNetwork & network, BcNetworkResource & cap_res,
                                    BcNetworkResource & time_res)
{
    int cap_res_UB = (std::min)(data.depots[depotId].capacity, data.vehicle_cap);
	for (int custId = 0; custId <= data.nbCustomers + 1; ++custId)
	{
		BcVertex vertex = network.createVertex();

		cap_res.setVertexConsumptionLB(vertex, 0);
		cap_res.setVertexConsumptionUB(vertex, cap_res_UB);

        if (time_res.isDefined())
        {
            if (custId == 0 || custId == data.nbCustomers + 1)
            {
                time_res.setVertexConsumptionLB(vertex, data.depots[depotId].tw_start);
                time_res.setVertexConsumptionUB(vertex, data.depots[depotId].tw_end);
            }
            else
            {
                time_res.setVertexConsumptionLB(vertex, data.customers[custId].tw_start);
                time_res.setVertexConsumptionUB(vertex, data.customers[custId].tw_end);
            }
        }

        if (custId == 0)
        {
            network.setPathSource(vertex);
        }
		else if (custId <= data.nbCustomers)
        {
		    vertex.setElementaritySet(custId);
		    vertex.setPackingSet(custId);
        }
		else
        {
		    network.setPathSink(vertex);
        }
	}
}

void lrp::RCSPSolver::buildArcs(int depotId, BcNetwork & network, BcNetworkResource & cap_res,
                                BcNetworkResource & time_res)
{
	BcVarArray xVar(spForm, "X");

    for (int firstCustId = 0; firstCustId <= data.nbCustomers; ++firstCustId )
        for (int secondCustId = 1; secondCustId <= data.nbCustomers + 1; ++secondCustId)
        {
            if (firstCustId == secondCustId)
                continue;

            int minCustId = std::min(firstCustId, (secondCustId <= data.nbCustomers) ? secondCustId : 0);
            int maxCustId = std::max(firstCustId, (secondCustId <= data.nbCustomers) ? secondCustId : 0);

            if (minCustId == maxCustId)
                continue;

            BcArc arc = network.createArc(firstCustId, secondCustId, 0.0);
            arc.arcVar((BcVar)xVar[depotId][minCustId][maxCustId]);

            cap_res.setArcConsumption(arc, ((minCustId == 0) ? 0 : data.customers[minCustId].demand * 0.5)
                                           + data.customers[maxCustId].demand * 0.5);

            if (time_res.isDefined())
            {
                if (minCustId == 0)
                    time_res.setArcConsumption(arc, data.getDepotToCustDistance(depotId, maxCustId)
                                                  + ((firstCustId != 0) ? data.customers[firstCustId].service_time : 0));
                else
                    time_res.setArcConsumption(arc, data.getCustToCustDistance(firstCustId, secondCustId)
                                                   + data.customers[firstCustId].service_time);
            }
        }
}

void lrp::RCSPSolver::buildElemSetDistanceMatrix(BcNetwork & network)
{
    int nbElemSets = data.nbCustomers + 1;
    std::vector<std::vector<double> > distanceMatrix(nbElemSets, std::vector<double>(nbElemSets, 1e12));

    for (int firstCustId = 1; firstCustId <= data.nbCustomers; ++firstCustId )
        for (int secondCustId = 1; secondCustId <= data.nbCustomers; ++secondCustId)
            distanceMatrix[firstCustId][secondCustId] = data.getCustToCustDistance(firstCustId, secondCustId);

    network.setElemSetsDistanceMatrix(distanceMatrix);
}


