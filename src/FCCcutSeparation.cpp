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
#include "bcModelingLanguageC.hpp"
#include "bcBapcodInit.hpp"
#include "bcModelC.hpp"

void minimalCoversEnumeration(const std::vector<int> & capacities, int curPosition, int remainingDemand,
                              std::vector<int> & curSolution, std::vector<std::vector<int> > & minimalCovers )
{
    if (remainingDemand <= 0)
    {
        minimalCovers.push_back(curSolution);
        return;
    }

    for (int position = curPosition + 1; position < (int)capacities.size(); ++position)
    {
        curSolution.push_back(position);
        minimalCoversEnumeration(capacities, position, remainingDemand - capacities[position], curSolution,
                                 minimalCovers);
        curSolution.pop_back();
    }
}

lrp::FenchelCutSeparationRoutine::FenchelCutSeparationRoutine(const Data & data_, const Parameters & params_) :
        data(data_), params(params_), cutCount(0), cutRound(0), printLevel(2), bcInitPtr(nullptr), modelPtr(nullptr),
        probConfPtr(nullptr), minimalCovers()
{
    int totalDemand = 0;
    for (int custId = 1; custId <= data.nbCustomers; ++custId)
        totalDemand += data.customers[custId].demand;
    std::vector<int> capacities(data.nbDepots + 1, 0);
    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        capacities[depotId] = data.depots[depotId].capacity;

    std::vector<int> curSolution;
    curSolution.reserve(data.nbDepots);
    minimalCoversEnumeration(capacities, 0, totalDemand, curSolution,minimalCovers);

    std::cout << "FCC separation inialization : " << minimalCovers.size() << " minimal covers generated." << std::endl;

    bcInitPtr = new BapcodInit();
    bcInitPtr->param().masterSolMode(SolutionMethod::lpSolver);
    bcInitPtr->param().MipSolverMultiThread(1);
    modelPtr = new Model(bcInitPtr);
    probConfPtr = modelPtr->createProbConf("FCCseparation", MultiIndex());
    BcFormulation form(probConfPtr);
    BcVarArray aVar(form, "A");
    aVar.type('C');
    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        aVar(depotId);
    BcConstrArray constrs(form, "Cstr");

    for (int coverId = 0; coverId < (int)minimalCovers.size(); ++coverId)
    {
        constrs(coverId) >= 1;
        for (auto depotId : minimalCovers[coverId])
            constrs[coverId] += 1 * aVar[depotId];
    }
}

lrp::FenchelCutSeparationRoutine::~FenchelCutSeparationRoutine()
{
    delete modelPtr;
    delete bcInitPtr;
}

int lrp::FenchelCutSeparationRoutine::operator()(BcFormulation bcForm, BcSolution &primalSol, double & maxViolation,
                                                 std::list<BcConstr> & cutList)
{
    std::set< BcVar > yVarSet;
    primalSol.extractVar("Y", yVarSet);

    BcFormulation form(probConfPtr);
    form.resetObjective();
    BcObjective objective(form);
    BcVarArray aVar(form, "A");

    std::vector<bool> objCoeffSet(data.nbDepots + 1, false);
    bool fractionalYvarValueExists = false;
    for (const auto & bcVar : yVarSet)
    {
        int depotId = bcVar.id().first();
        double yValue = bcVar.solVal();
        objective += yValue * aVar[depotId];
        objCoeffSet[depotId] = true;
        if (yValue < 1 - 1e-6)
            fractionalYvarValueExists = true;
    }

    if (!fractionalYvarValueExists)
        return 0;

    cutRound++;

    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        if (!objCoeffSet[depotId])
        {
            objective += 1e-6 * aVar[depotId];
        }

    form.update();
    bool toPrint = printL(printLevel);
    BcSolution solution = form.solve(toPrint, toPrint);
    if (toPrint)
        solution.print(std::cout);

    if (solution.cost() < 1 - 1e-6)
    {
        BcCutConstrArray fccCutConstr(bcForm, "FCC");
        BcConstr newCut = fccCutConstr(cutCount++);
        newCut >= 1;

        std::set< BcVar > aVarSet;
        solution.extractVar(aVarSet);

        BcVarArray yVar(bcForm, "Y");

        for (auto & bcVar : aVarSet)
        {
            int depotId = bcVar.id().first();
            double coeff = ceil(bcVar.solVal() * 1000000) / 1000000;
            //std::cout << "Coeff of var " << ((BcVar)yVar[depotId]).name() << " is " << coeff << std::endl;
            newCut += coeff * yVar[depotId];
        }

        cutList.push_back(newCut);
        return 1;
    }

    return 0;
}

