#include "../third_party/scTools/scTools.h"
#include "../include/WENOFV.hpp"
#include "../include/ReadConfig.hpp"

#include <iostream>
#include <fstream>
#include <mpi.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // --------------------- Two-Phase Flow Equation ---------------------//
    TwoPhaseEquation *equation = new TwoPhaseEquation();

    string inputfile = "./config/settings.cfg";

    // ---------------- Reading Configuration File - --------------//
    if (world_rank == 0)
        std::cout << "Starting to read configuration file......" << std::endl;

    CConfig config(inputfile);
    std::map<std::string, std::string> option;
    config.read(option);

    // ------------ Pass Equation and Initial Boundary Conditions to Solver -------------- //
    CWENOFV solver(equation, option);

    solver.run();

    MPI_Finalize();
    return 0;
}
