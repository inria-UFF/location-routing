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
 *  - VRPSolver extension (RCSP solver) v.0.5.14a
 *  - Boost C++ library
 */


#include <iostream>
#include "Loader.h"
#include "Model.h"
#include "SolutionChecker.h"

using namespace std;

int main(int argc, char** argv)
{
	BcInitialisation bcInit(argc, argv, "", false, true);

	lrp::Loader loader;
	if (!loader.loadParameters(bcInit.configFile(), argc, argv)
	    || !loader.loadData(bcInit.instanceFile()))
		return -1;

	lrp::SolutionChecker * sol_checker = new lrp::SolutionChecker;

	lrp::Model model(bcInit);
	model.attach(sol_checker);

	BcSolution solution(model.master());
	solution = model.solve();
    bool feasibleSol = (solution.defined()) ? sol_checker->isFeasible(solution, true) : false;

    bcInit.outputBaPCodStatistics(bcInit.instanceFile());

    std::cout << ">>-!-!-<<" << std::endl;
    std::cout << "Instance,InitUB,FeasFinalSol,bcFailToSolveModel,bcCountNodeProc,bcRecRootDb,bcRecBestDb,bcRecBestInc,"
              << "rootGap,bestGap,bcCountMastSol,bcCountCg,bcCountSpSol,bcCountCol,bcCountCutInMaster,"
              << "bcCountCutCAP,bcCountCutCOV,bcCountCutDCC,bcCountCutFCC,bcCountCutGUB,bcCountCutR1C,bcCountCutRCK,"
              << "bcCountCutDCCwithOneY,bcCountCutDCCwithTwoYs,bcCountCutRLKCRounding,bcCountCutRLKCsOne2,"
              << "bcCountCutRLKCsOne3,bcCountCutRLKCsOne4,bcCountCutRLKCsOne5,bcCountCutRLKCsOne6,"
              << "bcCountCutRLKCsOne8,bcCountCutRLKCsOne10,bcCountRootActiveCutCAP,bcCountRootActiveCutCOV,"
              << "bcCountRootActiveCutDCC,bcCountRootActiveCutFCC,bcCountRootActiveCutGUB,bcCountRootActiveCutR1C,"
              << "bcCountRootActiveCutRCK,bcCountRootActiveDCCwithOneY,bcCountRootActiveDCCwithTwoYs,"
              << "bcCountRootActiveRLKCsOne2,bcCountRootActiveRLKCsOne3,bcCountRootActiveRLKCsOne4,"
              << "bcCountRootActiveRLKCsOne5,bcCountRootActiveRLKCsOne6,bcCountRootActiveRLKCsOne8,"
              << "bcCountRootActiveRLKCsOne10,bcTimeMastMPsol,bcTimeColGen,bcTimeCutSeparation,bcTimeAddCutToMaster,"
              << "bcTimeRedCostFixAndEnum,bcTimeEnumMPsol,bcTimeRootEval,bcTimeBaP" << std::endl;
    std::string fullInstanceName = bcInit.instanceFile();
    std::replace(fullInstanceName.begin(), fullInstanceName.end(), ',', '-');

    std::size_t pos = fullInstanceName.find("/Set");
    std::string instanceName = fullInstanceName.substr(pos+1);

    double rootGap = (bcInit.getStatisticValue("bcRecBestInc") - bcInit.getStatisticValue("bcRecRootDb"))
                     / bcInit.getStatisticValue("bcRecRootDb");
    double bestGap = (bcInit.getStatisticValue("bcRecBestInc") - bcInit.getStatisticValue("bcRecBestDb"))
                     / bcInit.getStatisticValue("bcRecBestDb");
    std::cout.setf(std::ios::fixed);

    if (lrp::Data::getInstance().type == lrp::Data::VRPTWS)
    {
        double capFactor = lrp::Parameters::getInstance().depotCapacityFactor();
        char abc = (capFactor == 1.05) ? 'A' : (capFactor == 1.2) ? 'B' : 'C';
        std::string shortInstanceName = fullInstanceName.substr(pos+6);
        std::cout << "SetS/" << abc << "-" << lrp::Parameters::getInstance().nbCustomers() << "-"
                  << shortInstanceName << ",";
    } else {
        std::cout << instanceName << ",";
    }

    std::cout << lrp::Parameters::getInstance().cutOffValue() << ","
              << feasibleSol << ","
              << bcInit.getStatisticValue("bcFailToSolveModel") << ","
              << bcInit.getStatisticCounter("bcCountNodeProc") << ","
              << std::setprecision(3) << bcInit.getStatisticValue("bcRecRootDb") << ","
              << std::setprecision(3) << bcInit.getStatisticValue("bcRecBestDb") << ","
              << std::setprecision(3) << bcInit.getStatisticValue("bcRecBestInc") << ","
              << std::setprecision(3) << rootGap << ","
              << std::setprecision(3) << bestGap << ","
              << bcInit.getStatisticCounter("bcCountMastSol") << ","
              << bcInit.getStatisticCounter("bcCountCg") << ","
              << bcInit.getStatisticCounter("bcCountSpSol") << ","
              << bcInit.getStatisticCounter("bcCountCol") << ","
              << bcInit.getStatisticCounter("bcCountCutInMaster") << ","
              << bcInit.getStatisticCounter("bcCountCutCAP") << ","
              << bcInit.getStatisticCounter("bcCountCutCOV") << ","
              << bcInit.getStatisticCounter("bcCountCutDCC") << ","
              << bcInit.getStatisticCounter("bcCountCutFCC") << ","
              << bcInit.getStatisticCounter("bcCountCutGUB") << ","
              << bcInit.getStatisticCounter("bcCountCutR1C") << ","
              << bcInit.getStatisticCounter("bcCountCutRCK") << ","
              << bcInit.getStatisticCounter("bcCountCutDCCwithOneY") << ","
              << bcInit.getStatisticCounter("bcCountCutDCCwithTwoYs") << ","
              << bcInit.getStatisticCounter("bcCountCutRLKCRounding") << ","
              << bcInit.getStatisticCounter("bcCountCutRLKCsOne2") << ","
              << bcInit.getStatisticCounter("bcCountCutRLKCsOne3") << ","
              << bcInit.getStatisticCounter("bcCountCutRLKCsOne4") << ","
              << bcInit.getStatisticCounter("bcCountCutRLKCsOne5") << ","
              << bcInit.getStatisticCounter("bcCountCutRLKCsOne6") << ","
              << bcInit.getStatisticCounter("bcCountCutRLKCsOne8") << ","
              << bcInit.getStatisticCounter("bcCountCutRLKCsOne10") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveCutCAP") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveCutCOV") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveCutDCC") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveCutFCC") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveCutGUB") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveCutR1C") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveCutRCK") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveDCCwithOneY") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveDCCwithTwoYs") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveRLKCsOne2") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveRLKCsOne3") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveRLKCsOne4") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveRLKCsOne5") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveRLKCsOne6") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveRLKCsOne8") << ","
              << bcInit.getStatisticCounter("bcCountRootActiveRLKCsOne10") << ","
              << bcInit.getStatisticTime("bcTimeMastMPsol") << ","
              << bcInit.getStatisticTime("bcTimeColGen") << ","
              << bcInit.getStatisticTime("bcTimeCutSeparation") << ","
              << bcInit.getStatisticTime("bcTimeAddCutToMaster") << ","
              << bcInit.getStatisticTime("bcTimeRedCostFixAndEnum") << ","
              << bcInit.getStatisticTime("bcTimeEnumMPsol") << ","
              << bcInit.getStatisticTime("bcTimeRootEval") << ","
              << bcInit.getStatisticTime("bcTimeBaP")
              << std::endl;
    std::cout << ">>-!-!-<<" << std::endl;

	return 0;
}
