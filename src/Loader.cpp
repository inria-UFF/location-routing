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

#include "Loader.h"

#include <iostream>
#include <algorithm>
#include <locale>
#include <cmath>

#include "Data.h"
#include "Parameters.h"

lrp::Loader::Loader() :
    data(Data::getInstance()), parameters(Parameters::getInstance())
{
}

void lrp::Loader::applyCapacitiesRatioFactor()
{
    if ((parameters.capacitiesRatioFactor() <= 0.0) || (parameters.vehCapacityFactor() != 1.0)
        && (parameters.depotCapacityFactor() != 1.0))
    {
        for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        {
            data.depots[depotId].cost = (int) ceil(data.depots[depotId].cost * parameters.depotCapacityFactor());
            data.depots[depotId].capacity = (int) ceil(data.depots[depotId].capacity * parameters.depotCapacityFactor());
        }
        data.vehicle_cap = (int)ceil(data.vehicle_cap * parameters.vehCapacityFactor());
        return;
    }

    int totalDepotCapacity = 0;
    int maxDepotCapacity = 0;
    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        totalDepotCapacity += data.depots[depotId].capacity;
        maxDepotCapacity = (std::max)(maxDepotCapacity, data.depots[depotId].capacity);
    }
    double averDepotCapacity = (double)totalDepotCapacity / data.nbDepots;
    int targetDepotCapacity = totalDepotCapacity - maxDepotCapacity;

    int totalDemand = 0;
    for (int custId = 1; custId <= data.nbCustomers; ++custId)
        totalDemand += data.customers[custId].demand;

    double tentativeDepCapacityFactor = data.vehicle_cap / averDepotCapacity / parameters.capacitiesRatioFactor();
    double depCapacityFactor = tentativeDepCapacityFactor;

    if (targetDepotCapacity * tentativeDepCapacityFactor < totalDemand)
    {
        depCapacityFactor = (double)totalDemand / targetDepotCapacity;
        data.vehicle_cap = (int)ceil(parameters.capacitiesRatioFactor() * averDepotCapacity * depCapacityFactor);
    }

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        data.depots[depotId].capacity = (int)ceil(data.depots[depotId].capacity * depCapacityFactor);
        data.depots[depotId].cost = (int)ceil(data.depots[depotId].cost * depCapacityFactor);
    }
}

bool lrp::Loader::saveProdhonTypeFile(const std::string & fileName)
{
    std::ofstream ofs(fileName.c_str(), std::ios::out);
    if (!ofs)
    {
        std::cerr << "LRP instance writer error : cannot open or create file " << fileName << std::endl;
        return false;
    }

    ofs << data.nbCustomers << std::endl << data.nbDepots << std::endl << std::endl;

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        ofs << data.depots[depotId].x << "\t" << data.depots[depotId].y << std::endl;
    ofs << std::endl;

    for (int custId = 1; custId <= data.nbCustomers; ++custId)
        ofs << data.customers[custId].x << "\t" << data.customers[custId].y << std::endl;
    ofs << std::endl;

    ofs << data.vehicle_cap << std::endl << std::endl;

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        ofs << data.depots[depotId].capacity << std::endl;
    ofs << std::endl;

    for (int custId = 1; custId <= data.nbCustomers; ++custId)
        ofs << data.customers[custId].demand << std::endl;
    ofs << std::endl;

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        ofs << data.depots[depotId].cost << std::endl;
    ofs << std::endl;

    ofs << data.depots[1].veh_cost << std::endl << std::endl << "0" << std::endl << std::endl;

    ofs.close();

    return true;
}

bool lrp::Loader::loadVRPTWTypeFile(std::ifstream & ifs)
{
    data.type = Data::VRPTWS;
    data.roundType = Data::ROUND_DOWN_TO_ONE_DECIMAL_POINT;
    data.depots_open = true;
    data.distance_multiplier = 1;

    std::string line;
    std::getline(ifs, line);

    switch (line.at(0)) {
        case 'r':
        case 'R':
        case 'c':
        case 'C':
            break;
        default:
            std::cout << "VRPTWS instance reader error "<< std::endl;
            return false;
    }
    std::getline(ifs, line);
    std::getline(ifs, line);
    std::getline(ifs, line);
    int veh_number, veh_capacity;
    if (!(ifs >> veh_number >> veh_capacity))
    {
        std::cout << "VRPTWS instance reader error : cannot read vehicles number and vehicle capacity "<< std::endl;
        return false;
    }

    std::getline(ifs, line);
    std::getline(ifs, line);
    std::getline(ifs, line);
    std::getline(ifs, line);

    double depot_xCoord, depot_yCoord, depot_tw_begin, depot_tw_end;
    int totalDemand = 0;

    while (true)
    {
        int customerId = -1, demand, readyTime, dueDate, serviceTime;
        double xCoord, yCoord;
        if (!(ifs >> customerId >> xCoord >> yCoord >> demand >> readyTime >> dueDate >> serviceTime))
        {
            std::cout << "VRPTWS instance reader error : cannot read customers data" << std::endl;
            return false;
        }
        if (customerId == 0)
        {
            depot_xCoord = xCoord;
            depot_yCoord = yCoord;
            depot_tw_begin = readyTime;
            depot_tw_end = dueDate;
        }
        else
        {
            data.customers.push_back(Customer(customerId, demand, xCoord, yCoord, readyTime, dueDate, serviceTime));
            totalDemand += demand;
        }
        if (customerId == parameters.nbCustomers())
            break;
    }

    data.nbCustomers = parameters.nbCustomers();
    data.nbDepots = 3;
    data.vehicle_cap = veh_capacity;

    int depot_cap = (int)floor((totalDemand / 3.0) * parameters.depotCapacityFactor());
    double shift_wide = floor((depot_tw_end / 3.0) * 10) / 10;

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        data.depots.push_back(Depot(depotId, depot_cap, 0, veh_number, 0.0, depot_xCoord, depot_yCoord,
                                    (depotId - 2) * shift_wide, depot_tw_end));
    }

    return true;
}

bool lrp::Loader::loadImenTypeFile(std::ifstream & ifs)
{
    data.type = Data::MDCVRP;
    data.roundType = Data::NO_ROUND;
    data.depots_open = true;
    data.distance_multiplier = 1;

    int id, temp, nb_vehicles;
    if (!(ifs >> temp >> nb_vehicles >> data.nbCustomers >> data.nbDepots))
    {
        std::cout << "MDCVRP instance reader error : cannot read number of customers and depots "<< std::endl;
        return false;
    }

    double veh_cost, veh_cap, dbl_temp, x, y;

    if (!(ifs >> veh_cost >> veh_cap))
    {
        std::cout << "MDCVRP instance reader error : cannot read vehicle cost and vehicle capacity" << std::endl;
        return false;
    }

    for (int custId = 1; custId <= data.nbCustomers; ++custId)
    {
        double dbl_demand;
        if (!(ifs >> id >> x >> y >> dbl_temp >> dbl_demand >> dbl_temp >> dbl_temp >> dbl_temp >> dbl_temp >> dbl_temp))
        {
            std::cout << "LRP instance reader error : cannot read customers data " << std::endl;
            return false;
        }
        data.customers.push_back(Customer(custId, (int)ceil(dbl_demand * 10), x, y));
    }

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        double dbl_capacity;
        if (!(ifs >> id >> x >> y >> dbl_temp >> dbl_capacity >> dbl_temp >> dbl_temp >> dbl_temp >> dbl_temp >> dbl_temp))
        {
            std::cout << "LRP instance reader error : cannot read depots data " << std::endl;
            return false;
        }
        data.depots.push_back(Depot(depotId, (int)ceil(dbl_capacity * 10), 0, nb_vehicles, veh_cost, x, y));
    }

    data.vehicle_cap = (int)ceil(veh_cap * 10);

    applyCapacitiesRatioFactor();

    return true;
}


bool lrp::Loader::loadProdhonTypeFile(std::ifstream & ifs)
{
    data.type = Data::LRP;
    data.roundType = Data::ROUND_UP;
    data.depots_open = false;
    data.distance_multiplier = 100;

    if (!(ifs >> data.nbCustomers >> data.nbDepots))
    {
        std::cout << "LRP instance reader error : cannot read number of customers and depots "<< std::endl;
        return false;
    }

    data.customers.reserve(data.nbCustomers + 1);
    data.depots.reserve(data.nbDepots + 1);

    int x, y;

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        if (!(ifs >> x >> y))
        {
            std::cout << "LRP instance reader error : cannot read coordinates of depot "<< depotId << std::endl;
            return false;
        }
        data.depots.push_back(Depot(depotId, 0, 0, 100000, 0.0, x, y));
    }

    for (int custId = 1; custId <= data.nbCustomers; ++custId)
    {
        if (!(ifs >> x >> y))
        {
            std::cout << "LRP instance reader error : cannot read coordinates of customer "<< custId << std::endl;
            return false;
        }
        data.customers.push_back(Customer(custId, 0, x, y));
    }

    if (!(ifs >> data.vehicle_cap))
    {
        std::cout << "LRP instance reader error : cannot read vehicle capacity " << std::endl;
        return false;
    }

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        if (!(ifs >> data.depots[depotId].capacity))
        {
            std::cout << "LRP instance reader error : cannot read capacity of depot "<< depotId << std::endl;
            return false;
        }
    }

    for (int custId = 1; custId <= data.nbCustomers; ++custId)
    {
        if (!(ifs >> data.customers[custId].demand))
        {
            std::cout << "LRP instance reader error : cannot read demand of customer "<< custId << std::endl;
            return false;
        }
    }

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        if (!(ifs >> data.depots[depotId].cost))
        {
            std::cout << "LRP instance reader error : cannot read cost of depot "<< depotId << std::endl;
            return false;
        }
    }

    int veh_cost;

    if (!(ifs >> veh_cost))
    {
        std::cout << "LRP instance reader error : cannot read vehicle number " << std::endl;
        return false;
    }

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        data.depots[depotId].veh_cost = veh_cost;
        data.depots[depotId].veh_number = 10000;
    }

    applyCapacitiesRatioFactor();

    return true;
}

bool lrp::Loader::loadTuzunTypeFile(std::ifstream & ifs)
{
    data.type = Data::LRP;
    data.roundType = Data::NO_ROUND;
    data.depots_open = false;

    if (!(ifs >> data.nbCustomers >> data.nbDepots))
    {
        std::cout << "LRP instance reader error : cannot read number of customers and depots "<< std::endl;
        return false;
    }

    data.customers.reserve(data.nbCustomers + 1);
    data.depots.reserve(data.nbDepots + 1);

    double x, y;

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        if (!(ifs >> x >> y))
        {
            std::cout << "LRP instance reader error : cannot read coordinates of depot "<< depotId << std::endl;
            return false;
        }
        data.depots.push_back(Depot(depotId, 0, 0, 100000, 0.0, x, y));
    }

    for (int custId = 1; custId <= data.nbCustomers; ++custId)
    {
        if (!(ifs >> x >> y))
        {
            std::cout << "LRP instance reader error : cannot read coordinates of customer "<< custId << std::endl;
            return false;
        }
        data.customers.push_back(Customer(custId, 0, x, y));
    }

    if (!(ifs >> data.vehicle_cap))
    {
        std::cout << "LRP instance reader error : cannot read vehicle capacity " << std::endl;
        return false;
    }

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        if (!(ifs >> data.depots[depotId].capacity))
        {
            std::cout << "LRP instance reader error : cannot read capacity of depot "<< depotId << std::endl;
            return false;
        }
    }

    for (int custId = 1; custId <= data.nbCustomers; ++custId)
    {
        if (!(ifs >> data.customers[custId].demand))
        {
            std::cout << "LRP instance reader error : cannot read demand of customer "<< custId << std::endl;
            return false;
        }
    }

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        if (!(ifs >> data.depots[depotId].cost))
        {
            std::cout << "LRP instance reader error : cannot read cost of depot "<< depotId << std::endl;
            return false;
        }
    }

    double veh_cost;

    if (!(ifs >> veh_cost))
    {
        std::cout << "LRP instance reader error : cannot read vehicle cost " << std::endl;
        return false;
    }

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
    {
        data.depots[depotId].veh_cost = veh_cost;
        data.depots[depotId].veh_number = 10000;
    }

    applyCapacitiesRatioFactor();

    return true;
}

bool lrp::Loader::loadSchneiderTypeFile(std::ifstream & ifs)
{
    data.type = Data::LRP;
    data.roundType = Data::ROUND_UP;
    data.depots_open = false;
    data.distance_multiplier = 100;

    json::value json = json::parse(ifs);

    auto veh_cost = json["vehicle_costs"];
    auto veh_cap = json["vehicle_capacity"];

    if (!json::is_number(veh_cost) || !json::is_number(veh_cap))
    {
        std::cout << "LRP instance reader error : cannot read the vehicle data " << std::endl;
        return false;
    }

    data.vehicle_cap = json::to_number<int>(veh_cap);

    auto customers = json["customers"];

    int custId = 0;
    for (auto customer : as_array(customers))
    {
        json::value demand = customer["demand"];
        json::value index = customer["index"];
        json::value x = customer["x"];
        json::value y = customer["y"];
        if (!json::is_number(demand) || !json::is_number(x) || !json::is_number(y))
        {
            std::cout << "LRP instance reader error : cannot read the data of customer "
                      << json::to_string(index) << std::endl;
            return false;
        }
        data.customers.push_back(Customer(++custId, json::to_number<int>(demand), json::to_number<int>(x),
                                          json::to_number<int>(y)));
    }
    data.nbCustomers = custId;

    auto depots = json["depots"];

    int depotId = 0;
    for (auto depot : as_array(depots))
    {
        json::value capacity = depot["capacity"];
        json::value cost = depot["costs"];
        json::value index = depot["index"];
        json::value x = depot["x"];
        json::value y = depot["y"];
        if (!json::is_number(capacity) || !json::is_number(cost) || !json::is_number(x) || !json::is_number(y))
        {
            std::cout << "LRP instance reader error : cannot read the data of depot "
                      << json::to_string(index) << std::endl;
            return false;
        }
        data.depots.push_back(Depot(++depotId, json::to_number<int>(capacity), json::to_number<int>(cost),
                                    100000, json::to_number<double>(veh_cost), json::to_number<int>(x),
                                    json::to_number<int>(y)));
    }
    data.nbDepots = depotId;

    applyCapacitiesRatioFactor();

    return true;
}

bool lrp::Loader::loadParameters(const std::string & file_name, int argc, char* argv[])
{
    return parameters.loadParameters(file_name, argc, argv);
}

bool lrp::Loader::loadData(const std::string & file_name)
{
    std::ifstream ifs(file_name.c_str(), std::ios::in);
	if (!ifs)
	{
		std::cerr << "LRP instance reader error : cannot open file " << file_name << std::endl;
		return false;
	}

    if (ifs.eof())
    {
        std::cout << "LRP instance reader error : empty input file " << file_name << std::endl;
        ifs.close();
        return false;
    }

    data.name = file_name;
    bool success = false;

    if (file_name.find("LRP/SetP") != std::string::npos)
        success = loadProdhonTypeFile(ifs);
    else if (file_name.find("LRP/SetT") != std::string::npos)
        success = loadTuzunTypeFile(ifs);
    else if (file_name.find("LRP/SetS") != std::string::npos)
        success = loadSchneiderTypeFile(ifs);
    else if (file_name.find("MDCVRP/SetI") != std::string::npos)
        success = loadImenTypeFile(ifs);
    else if (file_name.find("CVRPS/SetS") != std::string::npos)
        success = loadVRPTWTypeFile(ifs);
    else
        std::cerr << "LRP instance reader error : cannot determine the instance type " << file_name << std::endl;

    printInstanceStats(file_name);
    ifs.close();

    if (parameters.saveInstanceFileName() != "")
        return saveProdhonTypeFile(parameters.saveInstanceFileName());

	return success;
}

void lrp::Loader::printInstanceStats(const std::string & file_name) {
    std::cout << "Instance " << file_name << " is loaded" << std::endl;
    std::cout << "Vehicle capacity is " << data.vehicle_cap << std::endl;

    std::cout << "Depot capacities are :";
    int totalCapacity = 0;
    std::vector<int> capacities(1, 0);
    for (int depotId = 1; depotId <= data.nbDepots; ++depotId) {
        std::cout << " " << data.depots[depotId].capacity;
        totalCapacity += data.depots[depotId].capacity;
        capacities.push_back(data.depots[depotId].capacity);
    }
    std::cout << std::endl << "Total depot capacity is " << totalCapacity << std::endl;
    int averCapacity = totalCapacity / data.nbDepots;
    std::cout << "Vehicle capacity is " << ((double) data.vehicle_cap / averCapacity) * 100
              << "% from the average depot capacity " << std::endl;

    int totalDemand = 0;
    for (int custId = 1; custId <= data.nbCustomers; ++custId)
        totalDemand += data.customers[custId].demand;
    std::cout << "Total customer demand is " << totalDemand << " ("
              << ((double) totalDemand / totalCapacity) * 100 << "% from total capacity)" << std::endl;

    int minNumOfDepots = -1;
    int maxNumOfDepots = -1;
    if (totalDemand > totalCapacity) {
        std::cout << "Total capacity is not sufficient to cover the total demand!!!" << std::endl;
    } else {
        std::sort(capacities.begin(), capacities.end());

        int depotId = 1;
        int curCapacity = capacities[1];
        while (curCapacity < totalDemand) {
            depotId += 1;
            curCapacity += capacities[depotId];
        }
        maxNumOfDepots = depotId;
        depotId = data.nbDepots;
        curCapacity = capacities[depotId];
        while (curCapacity < totalDemand) {
            depotId -= 1;
            curCapacity += capacities[depotId];
        }
        minNumOfDepots = data.nbDepots - depotId + 1;
        std::cout << "Between " << minNumOfDepots << " and " << maxNumOfDepots
                  << " depots are needed to cover the total demand " << std::endl;
    }

    std::cout << file_name << "\t" << data.vehicle_cap << "\t"
              << std::fixed << std::setprecision(1) << ((double) data.vehicle_cap / averCapacity) * 100 << "%\t";

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        std::cout << ((depotId == 1) ? "" : " ") << data.depots[depotId].capacity;

    std::cout << "\t" << totalCapacity << "\t" << totalDemand << "\t"
              << std::fixed << std::setprecision(1) << ((double) totalDemand / totalCapacity) * 100 << "%\t"
              << minNumOfDepots << "\t" << maxNumOfDepots << std::defaultfloat << std::setprecision(6) << std::endl;
}