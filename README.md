# Presentation

This distribution includes the source code and the data sets needed to reproduce the results of the paper

"Non-Robust Strong Knapsack Cuts for Capacitated Location-Routing and Related Problems" by Pedro Henrique Liguori, A. Ridha Mahjoub, Guillaume Marques, Ruslan Sadykov, and Eduardo Uchoa

published in the Operations Research journal.

# Data sets

There are data sets for three problems:

- Capacitated Location-Routing instances
    - Data set **PPW06** is in folder *data/LRP/SetP* 
    - Data set **TB99** is in folder *data/LRP/SetT* 
    - Data set **SL19** is in folder *data/LRP/SetS* 
- VRP-CMD instances
    - Data set is in folder *data/MDCVRP/SetI*
- VRPTW with Shifts instances
    - Data set is in folder *data/CVRPS/SetS*

Some instances were not used in the paper experiments. Please refer to tables with detailed results in the paper to see which instances were used. Some instances may require additional parameters to be passed. These parameters are reviewed below.

# Requirements

In order to run the code, one needs 

- Branch-Cut-and-Price solver **BaPCod** (version 0.74 and above) which can be downloaded for free (for academic use only) on its web-site https://bapcod.math.u-bordeaux.fr/
- MILP solver **IBM ILOG Cplex** (version 12.9 and above) which can be downloaded also for free (academic edition only) on its web-site https://www.ibm.com/products/ilog-cplex-optimization-studio.

BaPCod also requires free softwares **CMake**, **Boost**, and **LEMON**. Installation instructions for BaPCod, Boost, and LEMON libraries can be found in the README file distributed with BaPCod. 

# Compiling code

1. Install *Cplex*, *CMake*, *Boost*, *LEMON*, and *BaPCod*. We suppose that *BaPCod* is installed to the folder *BapcodFramework*. 
2. Create the folder *BapcodFramework/Applications/LocationRouting* and put all files from this distribution inside this folder. For example this file should have path *BapcodFramework/Applications/LocationRouting/README.md*
3. Open file *BapcodFramework/Applications/CMakeLists.txt* and add *LocationRouting* directory inside *add_subdirectories(...)*
4. Go to folder *BapcodFramework* and type the following on Linux or Mac OS
        
        cmake -B "build"  

    On Windows, use (has been tried only with Visual Studio 2019)
        
        cmake -G "Visual Studio 16 2019" -A x64 -B "build"

    Then, compilation can be launched with 

        cmake --build build --config Release --target LocationRouting -j 8


# Running code

Go to the folder *BapcodFramework/build/Applications/LocationRouting*. On Linux or Mac OS, run 

        bin/LocationRouting -b <path_to_config_file> -a config/app.cfg -i <path_to_instance_file> --cutOffValue <initial_PB> <additional options>

On Windows use,

        .\bin\Release\LocationRouting.exe -b <path_to_config_file> -a .\config\app.cfg -i <path_to_instance_file> --cutOffValue <initial_PB> <additional options>


Configuration file <path_to_config_file> determines the time limit and use of heuristics. Available configuration files are 
- **config/bc.cfg** -- without heuristic with 12 hours time limit
- **config/bc3.cfg** -- without heuristic with 3 hours time limit
- **config/bcHeur.cfg** -- with heuristic with 12 hours time limit
- **config/bcHeur30.cfg** -- with heuristic with 30 hours time limit
- **config/bcHeur70.cfg** -- with heuristic with 70 hours time limit
- **config/bcShiftInf.cfg** -- with heuristic with 30 hours time limit (only for VRPTW with Shifts instances)

Initial primal bound <initial_PB> is the value IPB given in the tables with detailed results, or
a known feasible solution. This parameter can be skipped the the run is done without initial upper bound. Please add a small epsilon value to this parameter if the best known solution value is used. Otherwise, no solution will found if <initial_PB> exactly equals to the optimum solution value. 

By default, all families of cutting planes are used. To switch off some families one can use options 
- **--useRLKC false** -- without RLKC family
- **--useDCC false** -- without DCC family
- **--useGUB false** -- without GUB family
- **--useFCC false** -- without FC family
- **--useCOV false** -- without COV family

The modified instances *PPW06* with modified ratio $\rho$ can be used by specifying option **--capacitiesRatioFactor <$\rho$ value>**.

VRP-CMD instances should be run with additional options **--depotCapacityFactor $\rho$**, where $\rho=1.05$ (instances A), $\rho=1.2$ (instances B), and $\rho=1.5$ (instances C), and **--nbCustomers $|J|$**, where $|J|$ is 25, 50, or 100.

# Reproducing tables with detailed results

An example of the command line to reproduce the results in **Table EC.1** in the electronic companion 

        bin/LocationRouting -b config/bcHeur30.cfg -a config/app.cfg -i data/LRP/SetP/coord50-5-2b.dat --cutOffValue 67309

An example of the command line to reproduce the results in **Table EC.2** in the electronic companion 

        bin/LocationRouting -b config/bcHeur30.cfg -a config/app.cfg -i data/LRP/SetP/coord50-5-2b.dat 

An example of the command line to reproduce the results in **Table EC.3** in the electronic companion 

        bin/LocationRouting -b config/bcHeur30.cfg -a config/app.cfg -i data/LRP/SetP/coord50-5-2b.dat --cutOffValue 67341


An example of the command line to reproduce the results in **Table EC.4** in the electronic companion 

        bin/LocationRouting -b config/bcHeur30.cfg -a config/app.cfg -i data/LRP/SetT/coordP111112.dat --cutOffValue 1467.78

An example of the command line to reproduce the results in **Table EC.5** in the electronic companion 

        bin/LocationRouting -b config/bcHeur30.cfg -a config/app.cfg -i data/LRP/SetS/100-5-1c.json --cutOffValue 134517

An example of the command line to reproduce the results in **Table EC.6** in the electronic companion 

        bin/LocationRouting -b config/bc3.cfg -a config/app.cfg -i data/MDCVRP/SetI/MDCVRP-10.dat --cutOffValue 12173.41

An example of the command line to reproduce the results in **Table EC.7** in the electronic companion 

        bin/LocationRouting -b config/bcHeur30.cfg -a config/app.cfg -i data/MDCVRP/SetI/MDCVRP-87.dat --cutOffValue 10521.17

An example of the command line to reproduce the results in **Tables EC.8 - EC.10** in the electronic companion 

        bin/LocationRouting -b config/bcShiftsInf.cfg -a config/app.cfg -i data/CVRPS/SetS/c101.dat --depotCapacityFactor 1.05 --nbCustomers 25



