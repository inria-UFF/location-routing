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

#ifndef LRP_LOADER_H
#define LRP_LOADER_H

#include <string>
#include "json.h"

namespace lrp
{
	class Data;
	class Parameters;

	class Loader
	{
	public:
		Loader();

		bool loadData(const std::string & file_name);
		bool loadParameters(const std::string & file_name, int argc, char* argv[]);

	private:
		Data & data;
		Parameters & parameters;

        void applyCapacitiesRatioFactor();

        bool saveProdhonTypeFile(const std::string & fileName);
        bool loadProdhonTypeFile(std::ifstream & ifs);
        bool loadTuzunTypeFile(std::ifstream & ifs);
        bool loadSchneiderTypeFile(std::ifstream & ifs);
        bool loadImenTypeFile(std::ifstream & ifs);
        bool loadVRPTWTypeFile(std::ifstream & ifs);

        void printInstanceStats(const std::string & file_name);
	};
}

#endif
