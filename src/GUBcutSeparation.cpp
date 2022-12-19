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

#include "CutSeparation.h"

lrp::GeneralizedUpperBoundSeparationRoutine::GeneralizedUpperBoundSeparationRoutine(const Data & data_,
                                                                                    const Parameters & params_):
  data(data_), params(params_), cutCount(0), cutRound(0)
{
}

lrp::GeneralizedUpperBoundSeparationRoutine::~GeneralizedUpperBoundSeparationRoutine() = default;

int lrp::GeneralizedUpperBoundSeparationRoutine::operator()(BcFormulation bcForm,
                                                            BcSolution & primalSol,
                                                            double & maxViolation,
                                                            std::list< BcConstr > & cutList)
{

    std::set< BcVar > yVarSet;
    primalSol.extractVar("Y", yVarSet);

    std::set< BcVar > varSet;
    std::vector<std::vector<double> > zSolution(data.nbDepots + 1, std::vector<double>(data.nbCustomers + 1, 0.0));
    if (params.useMasterZvars())
    {
        primalSol.extractVar("Z", varSet);
        for (const auto & bcVar : varSet)
        {
            int depotId = bcVar.id().first();
            int customerId = bcVar.id().second();
            zSolution[depotId][customerId] += bcVar.solVal();
        }
    }
    else
    {
        primalSol.extractVar("X", varSet);
        for (const auto & bcVar : varSet)
        {
            int depotId = bcVar.id().first();
            int firstCustId = bcVar.id().second();
            int secondCustId = bcVar.id().third();
            if (firstCustId > 0)
                zSolution[depotId][firstCustId] += bcVar.solVal() * 0.5;
            zSolution[depotId][secondCustId] += bcVar.solVal() * 0.5;
        }
    }

    cutRound++;

    BcCutConstrArray gubCutConstr(bcForm, "GUB");
    BcVarArray yVar(bcForm, "Y");
    BcVarArray zVar(bcForm, "Z");

    std::vector<BcVarArray> xVar(data.nbDepots + 1, BcVarArray());
    if (!params.useMasterZvars())
        for (auto & spForm : bcForm.colGenSubProblemList())
            xVar[spForm.id().first()] = BcVarArray(spForm, "X");

    int initCutCount = cutCount;
    for (const auto & bcVar : yVarSet)
    {
        int depotId = bcVar.id().first();
        for (int custId = 1; custId <= data.nbCustomers; ++custId)
        {
            if (zSolution[depotId][custId] > bcVar.solVal() + 1e-6)
            {
                BcConstr newCut = gubCutConstr(cutCount++);
                newCut <= 0;
                newCut += (-1) * yVar[depotId];
                if (params.useMasterZvars())
                {
                    newCut += 1 * zVar[depotId][custId];
                }
                else
                {
                    for (int beforeCustId = 0; beforeCustId < custId; ++beforeCustId)
                        newCut += 0.5 * xVar[depotId][depotId][beforeCustId][custId];
                    for (int afterCustId = custId + 1; afterCustId <= data.nbCustomers; ++afterCustId)
                        newCut += 0.5 * xVar[depotId][depotId][custId][afterCustId];
                }
                cutList.push_back(newCut);
            }
        }
    }

    return cutCount - initCutCount;
}