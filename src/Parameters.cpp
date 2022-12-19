/*
 *  Exact solver for the Capacitated Location-Routing and related problems
 *
 *  @author Ruslan Sadykov <Ruslan.Sadykov@inria.fr>, Inria Bordeaux, France (coordinator and main contributor),
 *  @author Eduardo Uchoa <eduardo.uchoa@gmail.com>, UFF, Brazil (scientific advise)
 *  @author Pedro Henrique Liguori <phliguori@gmail.com> BCG GAMMA, Brazil (separation of DCCs)
 *  @author Guillaume Marques <guillaume.marques@protonmail.com>, Atoptima, France (separation of GUBs and FCCs)
 *
 *  @date 27/04/2022
 *
 *  Dependencies:
 *  - BaPCod v.0.69
 *  - VRPSolver extension (RCSP solver) v.0.5.14
 */

#include "Parameters.h"

lrp::Parameters::Parameters() :
        silent("silent", false),
        cutOffValue("cutOffValue", std::numeric_limits<double>::infinity()),
        saveInstanceFileName("saveInstanceFileName", "", ""),
        useRLKC("useRLKC", true),
        useDCC("useDCC", true),
        useDCCtwoYs("useDCCtwoYs", true),
        useFCC("useFCC", true),
        useCOV("useCOV", true),
        useGUB("useGUB", true),
        DCCmaxNumPerRound("DCCmaxNumPerRound", 100),
        DCCheuristic("DCCheuristic", 2),
        useMasterZvars("useMasterZvars", false),
        useSubsetBranching("useSubsetBranching", true),
        depotCapacityFactor("depotCapacityFactor", 1.0),
        vehCapacityFactor("vehCapacityFactor", 1.0),
        capacitiesRatioFactor("capacitiesRatioFactor", 0.0),
        nbCustomers("nbCustomers", 25)
{
}

bool lrp::Parameters::loadParameters(const std::string & parameterFileName, int argc, char* argv[])
{
    setParameterFileName(parameterFileName);
    addApplicationParameter(silent);
    addApplicationParameter(cutOffValue);
    addApplicationParameter(saveInstanceFileName);
    addApplicationParameter(useRLKC);
    addApplicationParameter(useDCC);
    addApplicationParameter(useDCCtwoYs);
    addApplicationParameter(useFCC);
    addApplicationParameter(useCOV);
    addApplicationParameter(useGUB);
    addApplicationParameter(DCCmaxNumPerRound);
    addApplicationParameter(DCCheuristic);
    addApplicationParameter(useMasterZvars);
    addApplicationParameter(useSubsetBranching);
    addApplicationParameter(depotCapacityFactor);
    addApplicationParameter(vehCapacityFactor);
    addApplicationParameter(capacitiesRatioFactor);
    addApplicationParameter(nbCustomers);
    parse(argc, argv);

    return true;
}
