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

#ifndef LRP_SOLUTIONCHECKER_H
#define LRP_SOLUTIONCHECKER_H

#include <bcModelPointerC.hpp>
#include <bcModelRCSPSolver.hpp>
#include "InputUser.h"

#include <vector>

namespace lrp
{
    struct Route
    {
        int id;
        int	depotId;
        double cost;
        std::vector<int> vertIds;
        double capConsumption;
        double timeConsumption;

        Route(const BcSolution & solution, int id);
    };

    struct Solution
    {
        std::vector<Route> routes;
        double cost;
        bool feasible;

        explicit Solution(const BcSolution & solution);
    };

	class SolutionChecker : public BcSolutionFoundCallback, public InputUser
	{
	public:
		bool isFeasible(const BcSolution& solution, bool printSolution = false, bool best = false) const;
		bool operator()(BcSolution new_solution) const override;
	};
}

#endif
