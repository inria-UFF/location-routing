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

#ifndef LRP_PARAMETERS_H
#define LRP_PARAMETERS_H

#include "Singleton.h"
#include "bcParameterParserC.hpp"

namespace lrp
{
	class Parameters : public Singleton<Parameters>, ParameterParser
	{
		friend class Singleton<Parameters>;
	public:
		Parameters();
		virtual ~Parameters() {}

		bool loadParameters(const std::string & parameterFileName, int argc, char* argv[]);

        ApplicationParameter<double> cutOffValue;
        ApplicationParameter<bool> silent;
        ApplicationParameter<std::string> saveInstanceFileName;

        ApplicationParameter<bool> useRLKC;
        ApplicationParameter<bool> useDCC;
        ApplicationParameter<bool> useDCCtwoYs;
        ApplicationParameter<bool> useFCC;
        ApplicationParameter<bool> useCOV;
        ApplicationParameter<bool> useGUB;

        ApplicationParameter<int> DCCmaxNumPerRound;
        ApplicationParameter<int> DCCheuristic; /// 0 - Liguori (GRASP) heuristic without diversification
                                                /// 1 - Liguori (GRASP) heuristic with diversification
                                                /// 2 - greedy construction heuristic with the LHS sorting
                                                /// 3 - greedy construction heuristic with the violation sorting

        ApplicationParameter<bool> useMasterZvars;
        ApplicationParameter<bool> useSubsetBranching;

        ApplicationParameter<double> depotCapacityFactor; /// multiplier for the depot capacities, needed to make
                                                          /// experiments with modified literature instances
                                                          /// also needed for VRPTWS instances
        ApplicationParameter<double> vehCapacityFactor;  /// multiplier for the vehicle capacity, needed to make
                                                         /// experiments with modified literature instances
        ApplicationParameter<double> capacitiesRatioFactor;  /// target ratio veh.cap./avg.depot cap., modifies
                                                             /// depot and vehicle capacities according to Eduardo's
                                                             /// suggestion, active only if positive, and if
                                                             /// depotCapacityFactor = 1.0 and vehCapacityFactor = 1.0
        ApplicationParameter<int> nbCustomers; /// customers number for the VRPTWS instances
    };
}

#endif
