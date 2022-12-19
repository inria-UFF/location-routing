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

#ifndef LRP_SINGLETON_H
#define LRP_SINGLETON_H

#include <stdlib.h>

namespace lrp
{
	template<typename T>
	class Singleton
	{
	public:
		virtual ~Singleton() {}

		static T& getInstance()
		{
			if (singleton == NULL)
				singleton = new T();
			return *singleton;
		}

		static inline const T& getBuilt()
		{
			return *singleton;
		}

		void destroy()
		{
			delete singleton;
		}

	protected:
		Singleton() {}
		static T* singleton;
	};

	template<typename T>
	T* Singleton<T>::singleton = NULL;
}

#endif
