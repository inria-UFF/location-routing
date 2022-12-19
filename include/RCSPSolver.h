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

#ifndef LRP_RCSPSOLVER_H
#define LRP_RCSPSOLVER_H

#include "InputUser.h"

#include <vector>

#include "bcModelRCSPSolver.hpp"

namespace lrp
{
	class RCSPSolver : public InputUser
	{
	public:
		RCSPSolver(BcFormulation spForm, int depotId);
		virtual ~RCSPSolver() {}

        BcRCSPFunctor * getOraclePtr() { return oraclePtr; }

	private:
		BcFormulation spForm;
		int depotId;

        BcRCSPFunctor * oraclePtr;

		void buildVertices(int depotId, BcNetwork & network, BcNetworkResource & cap_res, BcNetworkResource & time_res);
		void buildArcs(int depotId, BcNetwork & network, BcNetworkResource & cap_res, BcNetworkResource & time_res);
		void buildElemSetDistanceMatrix(BcNetwork & network);
	};
}

#endif
