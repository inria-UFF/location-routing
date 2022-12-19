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

#ifndef LRP_DATA_H
#define LRP_DATA_H

#include "Singleton.h"

#include <vector>
#include <string>
#include <set>
#include <limits>
#include <cmath>

namespace lrp
{

    class Customer
    {
    public:
        Customer(int id = -1, int demand = 0, double x = 0.0, double y = 0.0,
                 double tw_start = -1e12, double tw_end = 1e12, double svc_time = 0.0) :
                id(id), demand(demand), x(x), y(y), tw_start(tw_start), tw_end(tw_end), service_time(svc_time)
        {}

        int id;
        int demand;
        double x;
        double y;
        double tw_start;
        double tw_end;
        double service_time;
    };

    class Depot
    {
    public:
        Depot(int id = -1, int capacity = 0, int cost = 0, int veh_number = 10000, double veh_cost = 0.0,
              double x = 0.0, double y = 0.0, double tw_start = 0.0, double tw_end = 1e12) :
                id(id), capacity(capacity), cost(cost), veh_number(veh_number), veh_cost(veh_cost),  x(x), y(y),
                tw_start(tw_start), tw_end(tw_end)
        {}

        int id;
        int capacity;
        int cost;
        int veh_number;
        double veh_cost;
        double x;
        double y;
        double tw_start;
        double tw_end;
    };

    class Data : public Singleton<Data>
	{
		friend class Singleton<Data>;
	public:

		std::string name;

		int vehicle_cap;
		bool depots_open;

        /// customers and depots are indexed starting from 1, customers[0] and depots[0] are fictive
		int nbCustomers;
		int nbDepots;
		std::vector<Customer> customers;
        std::vector<Depot> depots;

		enum Type
		{
			LRP,
			MDCVRP,
			VRPTWS
		};

        enum RoundType
        {
            NO_ROUND,
            ROUND_CLOSEST,
            ROUND_UP,
            ROUND_DOWN_TO_ONE_DECIMAL_POINT
        };

        Type type;
        RoundType roundType;
        double distance_multiplier;

        bool integerDistances() const
        {
            return roundType == ROUND_CLOSEST || roundType == ROUND_UP;
        }

        double getCustToCustDistance(int firstCustId, int secondCustId) const
        {
            return getDistance(customers[firstCustId].x, customers[firstCustId].y,
                               customers[secondCustId].x, customers[secondCustId].y);
        }

        double getDepotToCustDistance(int depotId, int custId) const
        {
            return getDistance(depots[depotId].x, depots[depotId].y, customers[custId].x, customers[custId].y);
        }

	private:
		Data() :
		    name(), vehicle_cap(0), depots_open(false), nbCustomers(0), nbDepots(0), customers(1, Customer()),
		    depots(1, Depot()), roundType(RoundType::NO_ROUND), type(Type::LRP), distance_multiplier(1.0)
		{}

		double getDistance(double x1, double y1, double x2, double y2) const
        {
		    double distance = sqrt((x2-x1) * (x2-x1) + (y2-y1) * (y2-y1) ) * distance_multiplier;
		    if (roundType == ROUND_CLOSEST)
		        return round(distance);
		    else if (roundType == ROUND_UP)
		        return ceil(distance);
            else if (roundType == ROUND_DOWN_TO_ONE_DECIMAL_POINT)
                return floor(distance * 10) / 10;
		    return distance;
        }
	};
}

#endif
