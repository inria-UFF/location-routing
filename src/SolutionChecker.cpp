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

#include "SolutionChecker.h"
#include "Data.h"

lrp::Route::Route(const BcSolution & solution, int id) :
        id(id), depotId(solution.formulation().id().first()), cost(solution.cost()), vertIds(), capConsumption(0.0),
        timeConsumption(0.0)
{
    const BcNetwork network(solution.formulation().network());
    vertIds.reserve(solution.orderedIds().size() + 1);

    if (solution.orderedIds().empty())
    {
        std::cout << "Error! A vertIds in a solution is empty!" << std::endl;
        return;
    }

    auto arcIdIt = solution.orderedIds().begin();
    vertIds.push_back(network.getArc(*arcIdIt).tail().ref());
    while (arcIdIt != solution.orderedIds().end())
    {
        vertIds.push_back(network.getArc(*arcIdIt).head().ref());
        ++arcIdIt;
    }
}

lrp::Solution::Solution(const BcSolution & solution) :
        routes(), cost(solution.cost()), feasible(false)
{
    BcSolution sol = solution.next(); /// skip master solution, go straight to the subproblem solutions

    int routeId = 0;
    while (sol.defined())
    {
        routes.emplace_back(sol, routeId++);
        sol = sol.next();
    }
}

bool lrp::SolutionChecker::isFeasible(const BcSolution & bcSolution, bool printSolution, bool best) const
{
    Solution solution(bcSolution);

    if (printSolution) {
        std::cout << "------------------------------------------ " << std::endl;
        if (best)
            std::cout << "Best found solution of value " << solution.cost << " : " << std::endl;
        else
            std::cout << "New solution of value " << solution.cost << " : " << std::endl;
    }

    std::vector<std::vector<int> > routeIdsByDepotId(data.nbDepots + 1);

    for (const auto & route : solution.routes)
        routeIdsByDepotId[route.depotId].push_back(route.id);

    solution.feasible = true;
    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        if (routeIdsByDepotId[depotId].empty())
            continue;

        int depotCumDemand = 0;
        int depotVehicleNumber = 0;

        for (int routeId : routeIdsByDepotId[depotId])
        {
            Route & route = solution.routes[routeId];
            if (route.vertIds.empty())
                continue;

            depotVehicleNumber += 1;
            if (printSolution)
                std::cout << "Vehicle " << depotVehicleNumber << " from depot " << route.depotId << " : ";

            int routeCumDemand = 0;
            double routeCumCost = 0.0;
            double routeCumTime = data.depots[depotId].tw_start;
            auto vertIt = route.vertIds.begin();
            int prevVertId = *vertIt;
            if (*vertIt != 0)
                solution.feasible = false; /// first vertex is not depot
            if (printSolution)
                std::cout << *vertIt;
            for (++vertIt; vertIt != route.vertIds.end(); ++vertIt)
            {
                bool vertexIsDepot = (*vertIt > data.nbCustomers);
                bool prevVertexIsDepot = (prevVertId == 0) || (prevVertId > data.nbCustomers);
                double distance = 0;
                if (vertexIsDepot)
                    distance = data.getDepotToCustDistance(depotId, prevVertId);
                else if (prevVertexIsDepot)
                    distance = data.getDepotToCustDistance(depotId, *vertIt);
                else
                    distance = data.getCustToCustDistance(prevVertId, *vertIt);

                routeCumCost += distance;
                double tw_start = (vertexIsDepot) ? data.depots[depotId].tw_start : data.customers[*vertIt].tw_start;
                double tw_end = (vertexIsDepot) ? data.depots[depotId].tw_end : data.customers[*vertIt].tw_end;
                routeCumTime = std::max(routeCumTime + distance, tw_start);

                if (printSolution)
                    std::cout << " -> " << *vertIt << "(" << routeCumTime << ")";

                if (routeCumTime > tw_end)
                {
                    solution.feasible = false;
                    if (printSolution)
                        std::cout << "!!!";
                }

                if (!vertexIsDepot)
                    routeCumDemand += data.customers[*vertIt].demand;

                prevVertId = *vertIt;
            }
            if ((prevVertId > 0) && (prevVertId <= data.nbCustomers))
                solution.feasible = false; /// last vertex is not depot
            if (printSolution)
                std::cout << ", demand = " << routeCumDemand << ", time = " << routeCumTime << ", cost = "
                          << routeCumCost << std::endl;
            route.capConsumption = routeCumDemand;
            route.timeConsumption = routeCumTime;
            if (routeCumDemand > data.vehicle_cap)
                solution.feasible = false; /// capacity is exceeded
            depotCumDemand += routeCumDemand;
        }
        if (depotCumDemand > data.depots[depotId].capacity)
            solution.feasible = false; /// depot cumulative capacity is exceeded
        if (depotVehicleNumber > data.depots[depotId].veh_number)
            solution.feasible = false; /// depot available number of vehicles is exceeded
    }
    if (printSolution)
    {
        std::cout << "Solution is " << ((solution.feasible) ? "" : "NOT ") << "feasible" << std::endl;
        std::cout << "------------------------------------------ " << std::endl;
    }

    return solution.feasible;
}

bool lrp::SolutionChecker::operator()(BcSolution new_solution) const
{
	return isFeasible(new_solution, true, false);
}
