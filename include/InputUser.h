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


#ifndef LRP_INPUTUSER_H
#define LRP_INPUTUSER_H

#include "Data.h"
#include "Parameters.h"

namespace lrp
{
    void generateAllSubsets(int minNumber, int maxNumber, int minSize, int maxSize,
                            std::vector<std::vector<int> > & subsets);

	class InputUser
	{
	protected:
		InputUser() : data(Data::getInstance()), params(Parameters::getInstance()) {}

		const Data & data;
		const Parameters & params;
	};
}

#endif
