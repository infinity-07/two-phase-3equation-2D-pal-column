#include <cmath>
#include <filesystem>
#include <iomanip>
#include <map>
#include <mpi.h>
#include <vector>

#include "../include/Equations.hpp"
#include "../include/Utility.hpp"
#include "../include/WENOFV.hpp"

CWENOFV::~CWENOFV()
{
}

void CWENOFV::initializeSolver()
{
    // Note: Variables containing 'world' refer to the entire computational domain

    // Get thread rank and size
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    m_numProcessorsX = 2; // Number of processors in X direction
    m_numProcessorsY = 2; // Number of processors in Y direction
    m_procIndexX = m_rank % m_numProcessorsX;
    m_procIndexY = m_rank / m_numProcessorsX;
    const int totalProcessors = m_numProcessorsX * m_numProcessorsY;

    if (m_size != totalProcessors)
    {
        std::cerr << "This program requires exactly " << totalProcessors << " processes." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (m_rank == 0)
        std::cout << "Initializing variables..." << std::endl;

    m_testcase = (TESTCASETYPE)std::stoi(option["TESTCASE"]);
    m_rkMeth = (RKMETHOD)std::stoi(option["RKMETH"]);

    m_riemannFluxType = (RIEMANNFLUXTYPE)std::stoi(option["RIEMANNFLUX"]);
    m_reconstructionType = (RECONSTRUCTIONTYPE)std::stoi(option["RECONSTRUCTIONTYPE"]);

    // Set equation parameters based on the testcase
    equation->setEquationParameters(m_testcase);
    m_timeMoments = equation->outputTimes;
    m_varNum = equation->getVarNum();

    // Initialize scheme, PPlimiter, and CFL number from options
    m_scheme = (SCHEMETYPE)std::stoi(option["SCHEME"]);
    m_pplimiter = (SCHEMETYPE)std::stoi(option["PPlimiter"]);
    m_cfl = std::stod(option["CFL"]);

    m_ghostCellNum = 3;

    // Initialize number of elements in X and Y directions
    m_worldElemNumX = std::stoi(option["ELEMNUMX"]);
    m_worldElemNumY = std::stoi(option["ELEMNUMY"]);

    m_worldStartElemX = m_ghostCellNum;
    m_worldStartElemY = m_ghostCellNum;
    m_worldEndElemX = m_worldStartElemX + m_worldElemNumX;
    m_worldEndElemY = m_worldStartElemY + m_worldElemNumY;

    if (m_worldElemNumX % m_numProcessorsX != 0 || m_worldElemNumY % m_numProcessorsY != 0)
    {
        std::cout << m_worldElemNumX << " " << m_numProcessorsX << " " << m_worldElemNumY << " " << m_numProcessorsY << std::endl;

        std::cerr << "The total number of elements in X and Y directions must be divisible by the number of processors in X and Y directions respectively." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Define element range including ghost cells
    // Calculate the range of elements this process will handle
    m_elemNumX = m_worldElemNumX / m_numProcessorsX;
    m_elemNumY = m_worldElemNumY / m_numProcessorsY;
    m_startElemX = m_ghostCellNum;
    m_startElemY = m_ghostCellNum;
    m_endElemX = m_startElemX + m_elemNumX;
    m_endElemY = m_startElemY + m_elemNumY;
    m_totalElemNumX = m_elemNumX + 2 * m_ghostCellNum;
    m_totalElemNumY = m_elemNumY + 2 * m_ghostCellNum;

    // Calculate element size in X and Y directions
    m_deltaX = (equation->xR - equation->xL) / m_worldElemNumX; // Structured grid
    m_deltaY = (equation->yR - equation->yL) / m_worldElemNumY; // Structured grid

    const double procXLeft = equation->xL + m_procIndexX * m_elemNumX * m_deltaX;
    const double procYLeft = equation->yL + m_procIndexY * m_elemNumY * m_deltaY;

    const double procGhostXLeft = procXLeft - m_deltaX * m_ghostCellNum;
    const double procGhostYLeft = procYLeft - m_deltaY * m_ghostCellNum;

    // Set the output directory based on the current timestamp and test case types
    m_outputDir = "./output/";
    m_outputDir += getTimestamp() + "_" + TESTCASETYPE_STRINGS[m_testcase] + "_";
    m_outputDir += SCHEMETYPE_STRINGS[m_scheme] + "_";
    m_outputDir += RECONSTRUCTIONTYPE_STRINGS[m_reconstructionType] + "_";
    m_outputDir += RIEMANNFLUXTYPE_STRINGS[m_riemannFluxType] + "_";
    m_outputDir += sc_common::intToString(m_worldElemNumX) + "x" + sc_common::intToString(m_worldElemNumY) + "/";

    if (m_rank == 0)
        std::cout << "Generating grid..." << std::endl;

    // Generate computational grid
    m_grids.Resize(m_totalElemNumX, m_totalElemNumY);
    for (int ei = 0; ei < m_totalElemNumX; ++ei)
    {
        for (int ej = 0; ej < m_totalElemNumY; ++ej)
        {
            const double loc_xl = procGhostXLeft + ei * m_deltaX;
            const double loc_xr = procGhostXLeft + (ei + 1) * m_deltaX;
            const double loc_yl = procGhostYLeft + ej * m_deltaY;
            const double loc_yr = procGhostYLeft + (ej + 1) * m_deltaY;

            m_grids[ei][ej] = Cgrid(loc_xl, loc_xr, loc_yl, loc_yr);
        }
    }

    if (m_rank == 0)
    {
        m_worldGrids.Resize(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);
        for (int ei = 0; ei < m_worldElemNumX + 2 * m_ghostCellNum; ++ei)
        {
            for (int ej = 0; ej < m_worldElemNumY + 2 * m_ghostCellNum; ++ej)
            {
                const double loc_xl = equation->xL - m_ghostCellNum * m_deltaX + ei * m_deltaX;
                const double loc_xr = equation->xL - m_ghostCellNum * m_deltaX + (ei + 1) * m_deltaX;
                const double loc_yl = equation->yL - m_ghostCellNum * m_deltaY + ej * m_deltaY;
                const double loc_yr = equation->yL - m_ghostCellNum * m_deltaY + (ej + 1) * m_deltaY;

                m_worldGrids[ei][ej] = Cgrid(loc_xl, loc_xr, loc_yl, loc_yr);
            }
        }
    }

    // Resize solution matrix and initialize variables
    m_Uh.Resize(m_totalElemNumX, m_totalElemNumY);
    for (int ei = 0; ei != m_totalElemNumX; ++ei)
        for (int ej = 0; ej != m_totalElemNumY; ++ej)
            m_Uh[ei][ej].vector.Resize(m_varNum);

    if (m_rank == 0)
    {
        // Initialize arrays to store global data
        m_worldUh.Resize(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);
        for (int ei = 0; ei != m_worldElemNumX + 2 * m_ghostCellNum; ei++)
            for (int ej = 0; ej != m_worldElemNumY + 2 * m_ghostCellNum; ej++)
                m_worldUh[ei][ej].vector.Resize(m_varNum);
    }

    // Resize and initialize Un matrix for storing solutions
    m_Un.Resize(m_totalElemNumX, m_totalElemNumY);
    for (int ei = 0; ei != m_totalElemNumX; ++ei)
        for (int ej = 0; ej != m_totalElemNumY; ++ej)
            m_Un[ei][ej].vector.Resize(m_varNum);

    // Resize and initialize RHS matrix
    m_rhs.Resize(m_totalElemNumX, m_totalElemNumY);
    for (int ei = 0; ei != m_totalElemNumX; ++ei)
        for (int ej = 0; ej != m_totalElemNumY; ++ej)
            m_rhs[ei][ej].vector.Resize(m_varNum);

    // Initialize number of Gauss points and resize reference points and weights
    m_gpNum = 4;
    m_gpoints_ref.Resize(m_gpNum);
    m_gweights_ref.Resize(m_gpNum);
    sc_math::GaussLobatto_ref(m_gpNum, m_gpoints_ref, m_gweights_ref);
    // Note: The sum of m_gweights_ref is 2

    // Resize and initialize flux matrices
    m_gridFlux.Resize(m_totalElemNumX, m_totalElemNumY);
    for (int ei = 0; ei < m_totalElemNumX; ++ei)
    {
        for (int ej = 0; ej < m_totalElemNumY; ++ej)
        {
            m_gridFlux[ei][ej].bottomFlux.Resize(m_varNum, m_gpNum);
            m_gridFlux[ei][ej].topFlux.Resize(m_varNum, m_gpNum);
            m_gridFlux[ei][ej].leftFlux.Resize(m_varNum, m_gpNum);
            m_gridFlux[ei][ej].rightFlux.Resize(m_varNum, m_gpNum);
        }
    }

    if (m_rank == 0)
    {
        m_worldGridFlux.Resize(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);
        for (int ei = 0; ei < m_worldElemNumX + 2 * m_ghostCellNum; ++ei)
        {
            for (int ej = 0; ej < m_worldElemNumY + 2 * m_ghostCellNum; ++ej)
            {
                m_worldGridFlux[ei][ej].bottomFlux.Resize(m_varNum, m_gpNum);
                m_worldGridFlux[ei][ej].topFlux.Resize(m_varNum, m_gpNum);
                m_worldGridFlux[ei][ej].leftFlux.Resize(m_varNum, m_gpNum);
                m_worldGridFlux[ei][ej].rightFlux.Resize(m_varNum, m_gpNum);
            }
        }
    }

    // PP-limiter theta
    m_totalTheta.Resize(m_totalElemNumX, m_totalElemNumY);
    if (m_rank == 0)
        m_worldTotalTheta.Resize(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);

    // DeltaT for every cell
    m_totalDeltaT.Resize(m_totalElemNumX, m_totalElemNumY);
    if (m_rank == 0)
        m_worldTotalDeltaT.Resize(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);

    m_totalPr.Resize(m_totalElemNumX, m_totalElemNumY);
    if (m_rank == 0)
        m_worldTotalPr.Resize(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);

    if (m_rank == 0)
        std::cout << "Outputing config document..." << std::endl;
    if (m_rank == 0)
        backupProjectFiles();

    initializeAve();

    MPI_Barrier(MPI_COMM_WORLD);
}
void CWENOFV::initializeAve(void)
{
    const int loc_gpNum = 4; // Number of Gauss points for local integration
    double gPointX, gWeightX;
    double gPointY, gWeightY;
    Array1D<double> loc_gpoints_ref(loc_gpNum), loc_gweights_ref(loc_gpNum);
    sc_math::GaussLegendre_ref(loc_gpNum, loc_gpoints_ref, loc_gweights_ref);

    Array1D<double> Conserved_var(m_varNum);

    if (m_rank == 0)
        std::cout << "Initializing cell averages..." << std::endl;

    for (int ei = m_startElemX; ei != m_endElemX; ++ei)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ++ej)
        {
            m_Uh[ei][ej].vector.setZero();

            // Loop over Gauss points for numerical integration
            for (int gpi = 0; gpi != loc_gpNum; ++gpi)
            {
                for (int gpj = 0; gpj != loc_gpNum; ++gpj)
                {
                    // Calculate the Gauss point location within the cell
                    gPointX = m_grids[ei][ej].m_xCenter + 0.5 * m_grids[ei][ej].m_xDistance * loc_gpoints_ref[gpi];
                    gPointY = m_grids[ei][ej].m_yCenter + 0.5 * m_grids[ei][ej].m_yDistance * loc_gpoints_ref[gpj];

                    // Calculate the Gauss weights for the cell
                    gWeightX = 0.5 * m_grids[ei][ej].m_xDistance * loc_gweights_ref[gpi];
                    gWeightY = 0.5 * m_grids[ei][ej].m_yDistance * loc_gweights_ref[gpj];

                    // Get the initial conserved variables at the Gauss point
                    equation->getU0(gPointX, gPointY, Conserved_var);

                    // Accumulate the weighted conserved variables
                    for (int r = 0; r != m_varNum; ++r)
                        m_Uh[ei][ej].vector[r] += Conserved_var[r] * gWeightX * gWeightY;
                }
            }
            // Normalize by the cell area to get the average value
            for (int r = 0; r != m_varNum; ++r)
                m_Uh[ei][ej].vector[r] = m_Uh[ei][ej].vector[r] / m_grids[ei][ej].m_area;
        }
    }
    if (m_rank == 0)
        std::cout << "Cell averages initialization completed..." << std::endl;
}

void CWENOFV::run(void)
{
    ///////////////////////////////////////////////////////////////////////////////
    // Main execution function for the WENO finite volume solver
    // Handles time integration using different Runge-Kutta methods
    // Includes parallel processing support using MPI
    ///////////////////////////////////////////////////////////////////////////////

    m_solverTimer.reset("solver");
    m_mainTimer.reset("main", &m_solverTimer);
    m_wenoTimer.reset("weno", &m_mainTimer);
    m_ppLimiterTimer.reset("ppLimiter", &m_mainTimer);
    m_FVTimer.reset("FV", &m_mainTimer);
    m_outputTimer.reset("output", &m_solverTimer);
    m_initialiTimer.reset("initialization", &m_solverTimer);
    m_riemannTimer.reset("riemann", &m_mainTimer);
    m_MPITimer.reset("mpi", &m_solverTimer);

    m_solverTimer.start();

    m_initialiTimer.start();
    initializeSolver();
    m_initialiTimer.pause();

    // Initialize current time
    double now = 0;
    int count(0);

    int m_currentOutputIndex = 1; // 初始指向第一个输出时刻

    m_outputTimer.start();
    outputResults(count, now, "initial");
    m_outputTimer.pause();

    m_mainTimer.start();

    // // Main time integration loop
    // while (fabs(now - m_timeMoments.back()) > 1e-10)
    // {
    //     performTimeIntegration(count, now);

    //     // Output progress information
    //     if (m_rank == 0 && count % 10 == 0)
    //         progressTimer.updateProgress(now / m_timeMoments.back());

    //     // output results
    //     outputResults(count, now);
    // }

    // outputResults(count, now, "final");

    // 主循环改为处理多个输出时间点
    while (m_currentOutputIndex < m_timeMoments.size())
    {
        const double targetTime = m_timeMoments[m_currentOutputIndex];

        // 时间积分到当前目标时刻
        while (std::abs(now - targetTime) > 1e-10)
        {
            // 计算允许的最大时间步
            const double remainingTime = targetTime - now;
            double timestep = calculateDeltaT(now);
            timestep = std::min(timestep, remainingTime);

            performTimeIntegration(timestep);
            now += timestep;
            count++;

            // 更新进度
            if (m_rank == 0 && count % 10 == 0)
                m_solverTimer.updateProgress(now / m_timeMoments.back());

            if (m_rank == 0 && count % 1000 == 0)
            {
                std::cout << "\n"
                          << std::endl;

                m_solverTimer.printHierarchy();

                std::cout << "\n"
                          << std::endl;
            }

            m_mainTimer.pause();
            m_outputTimer.start();
            outputResults(count, now);
            m_outputTimer.pause();
            m_mainTimer.start();
        }

        // 到达目标时刻后输出
        m_mainTimer.pause();
        m_outputTimer.start();
        outputResults(count, now, "intermediate");
        m_outputTimer.pause();
        m_mainTimer.start();

        m_currentOutputIndex++;
    }

    // 最终输出
    m_mainTimer.pause();
    m_outputTimer.start();
    outputResults(count, now, "intermediate");
    m_outputTimer.pause();

    if (m_rank == 0)
        m_solverTimer.printHierarchy();
}

void CWENOFV::RunRK1(double deltaT)
{
    // Assemble the right-hand side (RHS) for the equations
    assembleRHS();

    // Update solution vector using the RK2 coefficients
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
            for (int r = 0; r < m_varNum; r++)
                m_Uh[ei][ej].vector[r] += deltaT * m_rhs[ei][ej].vector[r];
}
void CWENOFV::RunRK2(double deltaT)
{
    // Copy current solution to temporary storage
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
            for (int r = 0; r != m_varNum; ++r)
                m_Un[ei][ej].vector[r] = m_Uh[ei][ej].vector[r];

    for (int step = 0; step != 2; step++)
    {
        // Assemble the right-hand side (RHS) for the equations
        assembleRHS();

        double a, b, c;
        // Set coefficients for the second-order Runge-Kutta scheme based on the step
        switch (step)
        {
        case 0:
            a = 1.0;
            b = 0.0;
            c = 1.0;
            break;
        case 1:
            a = 0.5;
            b = 0.5;
            c = 0.5;
            break;
        default:
            std::cout << "Error: Invalid RK2 step" << std::endl;
            std::cin.get();
            exit(1);
            break;
        }

        // Update solution vector using the RK2 coefficients
        for (int ei = m_startElemX; ei != m_endElemX; ei++)
            for (int ej = m_startElemY; ej != m_endElemY; ej++)
                for (int r = 0; r < m_varNum; r++)
                    m_Uh[ei][ej].vector[r] = a * m_Un[ei][ej].vector[r] + b * m_Uh[ei][ej].vector[r] + c * deltaT * m_rhs[ei][ej].vector[r];
    }
}
void CWENOFV::RunRK3(double deltaT)
{
    // Copy current solution to temporary storage
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
            for (int r = 0; r != m_varNum; ++r)
                m_Un[ei][ej].vector[r] = m_Uh[ei][ej].vector[r];

    for (int step = 0; step != 3; step++)
    {
        // Assemble the right-hand side (RHS) for the equations
        assembleRHS();

        double a, b, c;
        // Set coefficients for the third-order Runge-Kutta scheme based on the step
        switch (step)
        {
        case 0:
            a = 1.0;
            b = 0.0;
            c = 1.0;
            break;
        case 1:
            a = 3.0 / 4.0;
            b = 1.0 / 4.0;
            c = 1.0 / 4.0;
            break;
        case 2:
            a = 1.0 / 3.0;
            b = 2.0 / 3.0;
            c = 2.0 / 3.0;
            break;
        default:
            std::cout << "Error: Invalid RK3 step" << std::endl;
            std::cin.get();
            exit(1);
            break;
        }

        // Update solution vector using the RK3 coefficients
        for (int ei = m_startElemX; ei != m_endElemX; ei++)
            for (int ej = m_startElemY; ej != m_endElemY; ej++)
                for (int r = 0; r < m_varNum; r++)
                    m_Uh[ei][ej].vector[r] = a * m_Un[ei][ej].vector[r] + b * m_Uh[ei][ej].vector[r] + c * deltaT * m_rhs[ei][ej].vector[r];
    }
}
double CWENOFV::calculateDeltaT(double now)
{
    // Initialize timestep to a large value
    double timestep(1.0);

    int localExistNan(0);

    // Temporary variables for calculations
    Array1D<double> loc_Uh(m_varNum);
    loc_Uh.setZero();

    // calculate the time step
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            // Extract average solution values for the cell
            for (int r = 0; r != m_varNum; ++r)
                loc_Uh[r] = m_Uh[ei][ej].vector[r];

            // Calculate local timestep based on CFL condition
            const double eigenvalueX = equation->getMaxEigenValue(loc_Uh, 1, 0);
            const double eigenvalueY = equation->getMaxEigenValue(loc_Uh, 0, 1);
            double tmp = m_cfl / (eigenvalueX / m_deltaX + eigenvalueY / m_deltaY);

            // tmp = 10.0 * m_deltaX * m_deltaX; // for 6-order accuracy test

            m_totalDeltaT[ei][ej] = tmp;
            if (isnan(tmp))
            {
                localExistNan = 1;
                std::cout << loc_Uh << "\n"
                          << std::endl;
            }

            // Update timestep if smaller timestep is found
            if (tmp > timestep)
                continue;
            else
                timestep = tmp;
        }
    }

    // Reduce existNan across all processes
    int globalExistNan;
    MPI_Allreduce(&localExistNan, &globalExistNan, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (globalExistNan > 0)
    {
        if (m_rank == 0)
            std::cout << "Error: exist nan DeltaT" << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);

        outputDeltaT(sc_common::doubleToString(now));
        std::cin.get();
        exit(1);
    }

    // Obtain the minimum timestep across all processes using MPI_Reduce
    const double localTimestep = timestep;
    double globalTimestep;
    MPI_Reduce(&localTimestep, &globalTimestep, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Bcast(&globalTimestep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return globalTimestep;
}
void CWENOFV::performTimeIntegration(double timestep)
{
    // // Time integration step
    // const double deltaT = calculateDeltaT(now);

    switch (m_rkMeth)
    {
    case RK1:
        RunRK1(timestep);
        break;
    case RK2:
        RunRK2(timestep);
        break;
    case RK3:
        RunRK3(timestep);
        break;
    default:
        std::cout << "Error: Invalid Runge-Kutta method" << std::endl;
        std::cin.get();
        exit(1);
    }

    // // Update current time and count iterations
    // now = now + deltaT;
    // count++;
}

void CWENOFV::assembleRHS(void)
{
    for (int ei = 0; ei != m_totalElemNumX; ei++)
        for (int ej = 0; ej != m_totalElemNumY; ej++)
            m_rhs[ei][ej].vector.setZero();

    m_wenoTimer.start();
    // Calculate numerical fluxes using WENO scheme
    switch (m_reconstructionType)
    {
    case CONVERVATIVE: // reconstruct conservative variables
        getFlux_conservative();
        break;
    case CHARACTERISTIC: // reconstruct characteristic variables
        getFlux_characteristic();
        break;
    case PRIMITIVE: // reconstruct primitive variables
        getFlux_primitive();
        break;
    case PRIMITIVE2: // reconstruct primitive variables(a imporved version)
        getFlux_primitive2();
        break;
    case PRIMITIVE3: // reconstruct primitive variables(a imporved version)
        getFlux_WENOEXPandTHINC();
        break;
    case ThincOnlyForZ1: // reconstruct primitive variables(a imporved version)
        getFlux_WENOEXPandTHINC();
        break;
    default:
        std::cerr << "Error: Invalid reconstruction type" << std::endl;
        break;
    }
    m_wenoTimer.pause();

    // useFVLimiter();

    // Apply bound-preserving limiter
    m_ppLimiterTimer.start();
    if (m_pplimiter == true)
        useBoundPreservingLimiter();
    m_ppLimiterTimer.pause();

    // Compute numerical flux using Riemann solver
    assembleBoundaryFaceTerm();

    // Assemble source term for specific test cases (e.g., gravity term for DAMBREAK)
    assembleSourceTerm();

    // Debug:
    // for (int ei = 0; ei != m_totalElemNumX; ei++)
    // {
    //     for (int ej = 0; ej != m_totalElemNumY; ej++)
    //     {
    //         for (int k = 0; k != 2; k++)
    //         {
    //             const double zkrhok = m_Uh[ei][ej].vector[k];
    //             if (zkrhok < 0)
    //             {
    //                 cout << "Cell-averaged solution of Z" << k << "RHO" << k << " is not positive !" << endl;
    //                 std::cout << scientific << zkrhok << std::endl;
    //                 std::cin.get();
    //             }
    //         }
    //     }
    // }

    // for (int ei = 0; ei != m_totalElemNumX; ei++)
    // {
    //     for (int ej = 0; ej != m_totalElemNumY; ej++)
    //     {
    //         for (int r = 0; r != m_varNum; r++)
    //         {
    //             const double uave = m_Uh[ei][ej].vector[r];
    //             if (isnan(uave))
    //             {
    //                 std::cout << scientific << uave << std::endl;
    //                 std::cin.get();
    //             }
    //         }
    //     }
    // }
}

void CWENOFV::getFlux_conservative(void)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    /////// Implement the 2D WENO scheme using a dimension by dimension approach. /////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Li, Z., Vittoz, L., Oger, G., & Le Touzé, D. (2022).
    // A simplified and efficient weakly-compressible FV-WENO scheme for immiscible two-phase flows.
    // Computers & Fluids, 244. https://doi.org/10.1016/j.compfluid.2022.105555

    // Warning: only for GaussType = Lobatto N = 4

    m_wenoTimer.pause();
    m_mainTimer.pause();

    getWorldUh();
    worldToProcUh();

    m_mainTimer.start();
    m_wenoTimer.start();

    m_totalPr.setZero();

    const int tempNum = 7; // Number of stencil points

    Array1D<double> loc_Uh(m_varNum);
    Array2D<Cvector> Uave(tempNum, tempNum);
    Array2D<Cvector> Ugp(m_gpNum, m_gpNum);
    Array2D<Cvector> Utemple(tempNum, m_gpNum);

    for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
        for (int gp = 0; gp != m_gpNum; gp++)
            Utemple[tempIndex][gp].vector.Resize(m_varNum);

    for (int i = 0; i != tempNum; ++i)
        for (int j = 0; j != tempNum; ++j)
            Uave[i][j].vector.Resize(m_varNum);

    for (int gpi = 0; gpi != m_gpNum; gpi++)
        for (int gpj = 0; gpj != m_gpNum; gpj++)
            Ugp[gpi][gpj].vector.Resize(m_varNum);

    // use WENO limiter for every cell
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            //(1)Step-1: set stencil
            for (int r = 0; r != m_varNum; ++r)
                for (int eii = 0; eii != tempNum; eii++)
                    for (int ejj = 0; ejj != tempNum; ejj++)
                        Uave[eii][ejj].vector[r] = m_Uh[ei + eii - 3][ej + ejj - 3].vector[r];

            //(3)Step-3: weno restruction in y-direction
            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        double uavemmm = Uave[tempIndex][0].vector[r];
                        double uavemm = Uave[tempIndex][1].vector[r];
                        double uavem = Uave[tempIndex][2].vector[r];
                        double uavec = Uave[tempIndex][3].vector[r];
                        double uavep = Uave[tempIndex][4].vector[r];
                        double uavepp = Uave[tempIndex][5].vector[r];
                        double uaveppp = Uave[tempIndex][6].vector[r];

                        Utemple[tempIndex][gp].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp);
                    }
                }
            }

            //(6)Step-6: weno restruction in x-direction
            for (int gpj = 0; gpj != m_gpNum; ++gpj) // y
            {
                for (int gpi = 0; gpi != m_gpNum; ++gpi) // x
                {
                    if (gpi == 0 || gpi == m_gpNum - 1 || gpj == 0 || gpj == m_gpNum - 1)
                    {
                        for (int r = 0; r != m_varNum; ++r)
                        {
                            double uavemmm = Utemple[0][gpj].vector[r];
                            double uavemm = Utemple[1][gpj].vector[r];
                            double uavem = Utemple[2][gpj].vector[r];
                            double uavec = Utemple[3][gpj].vector[r];
                            double uavep = Utemple[4][gpj].vector[r];
                            double uavepp = Utemple[5][gpj].vector[r];
                            double uaveppp = Utemple[6][gpj].vector[r];

                            Ugp[gpi][gpj].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi);
                        }
                    }
                }
            }

            // Step-6: get flux
            for (int r = 0; r != m_varNum; r++)
            {
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    m_gridFlux[ei][ej].bottomFlux[r][gp] = Ugp[gp][0].vector[r];
                    m_gridFlux[ei][ej].topFlux[r][gp] = Ugp[gp][m_gpNum - 1].vector[r];
                    m_gridFlux[ei][ej].leftFlux[r][gp] = Ugp[0][gp].vector[r];
                    m_gridFlux[ei][ej].rightFlux[r][gp] = Ugp[m_gpNum - 1][gp].vector[r];
                }
            }
        }
    }
}
void CWENOFV::getFlux_characteristic(void)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    /////// Implement the 2D WENO scheme using a dimension by dimension approach. /////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Li, Z., Vittoz, L., Oger, G., & Le Touzé, D. (2022).
    // A simplified and efficient weakly-compressible FV-WENO scheme for immiscible two-phase flows.
    // Computers & Fluids, 244. https://doi.org/10.1016/j.compfluid.2022.105555

    // Warning: only for GaussType = Lobatto N = 4
    m_wenoTimer.pause();
    m_mainTimer.pause();
    getWorldUh();
    worldToProcUh();
    m_mainTimer.start();
    m_wenoTimer.start();

    m_totalPr.setZero();

    const int tempNum = 7; // Number of stencil points

    Array1D<double> loc_Uh(m_varNum);
    Array2D<double> eigenMatrixL(m_varNum, m_varNum), eigenMatrixR(m_varNum, m_varNum);
    Array2D<Cvector> Uave(tempNum, tempNum), Uave_c(tempNum, tempNum);
    Array2D<Cvector> Ugp_c(m_gpNum, m_gpNum), Ugp(m_gpNum, m_gpNum);
    Array2D<Cvector> Utemple_c(tempNum, m_gpNum);
    Array2D<Cvector> Utemple(tempNum, m_gpNum);

    for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
    {
        for (int gp = 0; gp != m_gpNum; gp++)
        {
            Utemple[tempIndex][gp].vector.Resize(m_varNum);
            Utemple_c[tempIndex][gp].vector.Resize(m_varNum);
        }
    }

    for (int i = 0; i != tempNum; ++i)
    {
        for (int j = 0; j != tempNum; ++j)
        {
            Uave[i][j].vector.Resize(m_varNum);
            Uave_c[i][j].vector.Resize(m_varNum);
        }
    }

    for (int gpi = 0; gpi != m_gpNum; gpi++)
    {
        for (int gpj = 0; gpj != m_gpNum; gpj++)
        {
            Ugp_c[gpi][gpj].vector.Resize(m_varNum);
            Ugp[gpi][gpj].vector.Resize(m_varNum);
        }
    }

    // use WENO limiter for every cell
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            //(1)Step-1: set stencil
            for (int r = 0; r != m_varNum; ++r)
                for (int eii = 0; eii != tempNum; eii++)
                    for (int ejj = 0; ejj != tempNum; ejj++)
                        Uave[eii][ejj].vector[r] = m_Uh[ei + eii - 3][ej + ejj - 3].vector[r];

            //(2)Step-2: y-direction characteristic projection
            for (int i = 0; i != tempNum; ++i)
            {
                loc_Uh = Uave[i][3].vector;
                equation->getLEigenMatrix(loc_Uh, 0, 1, eigenMatrixL);

                for (int j = 0; j != tempNum; ++j)
                    Uave_c[i][j].vector = eigenMatrixL * Uave[i][j].vector;
            }

            //(3)Step-3: weno restruction in y-direction
            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        double uavemmm = Uave_c[tempIndex][0].vector[r];
                        double uavemm = Uave_c[tempIndex][1].vector[r];
                        double uavem = Uave_c[tempIndex][2].vector[r];
                        double uavec = Uave_c[tempIndex][3].vector[r];
                        double uavep = Uave_c[tempIndex][4].vector[r];
                        double uavepp = Uave_c[tempIndex][5].vector[r];
                        double uaveppp = Uave_c[tempIndex][6].vector[r];

                        Utemple_c[tempIndex][gp].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp);
                    }
                }
            }

            //(4)Step-4: y-direction project back to physical space
            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
            {
                loc_Uh = Uave[tempIndex][3].vector;
                equation->getREigenMatrix(loc_Uh, 0, 1, eigenMatrixR);

                for (int gp = 0; gp != m_gpNum; ++gp)
                    Utemple[tempIndex][gp].vector = eigenMatrixR * Utemple_c[tempIndex][gp].vector;
            }

            //(5)Step-5: x-direction characteristic projection compute eigenMatrix
            loc_Uh = Uave[3][3].vector;
            equation->getLEigenMatrix(loc_Uh, 1, 0, eigenMatrixL);
            equation->getREigenMatrix(loc_Uh, 1, 0, eigenMatrixR);

            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
                for (int gp = 0; gp != m_gpNum; ++gp)
                    Utemple_c[tempIndex][gp].vector = eigenMatrixL * Utemple[tempIndex][gp].vector;

            //(6)Step-6: weno restruction in x-direction
            for (int gpj = 0; gpj != m_gpNum; ++gpj) // y
            {
                for (int gpi = 0; gpi != m_gpNum; ++gpi) // x
                {
                    if (gpi == 0 || gpi == m_gpNum - 1 || gpj == 0 || gpj == m_gpNum - 1)
                    {
                        for (int r = 0; r != m_varNum; ++r)
                        {
                            double uavemmm = Utemple_c[0][gpj].vector[r];
                            double uavemm = Utemple_c[1][gpj].vector[r];
                            double uavem = Utemple_c[2][gpj].vector[r];
                            double uavec = Utemple_c[3][gpj].vector[r];
                            double uavep = Utemple_c[4][gpj].vector[r];
                            double uavepp = Utemple_c[5][gpj].vector[r];
                            double uaveppp = Utemple_c[6][gpj].vector[r];

                            Ugp_c[gpi][gpj].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi);
                        }

                        //(7)Step-7: x-direction project back to physical space
                        Ugp[gpi][gpj].vector = eigenMatrixR * Ugp_c[gpi][gpj].vector;
                    }
                }
            }

            // Step-6: get flux
            for (int r = 0; r != m_varNum; r++)
            {
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    m_gridFlux[ei][ej].bottomFlux[r][gp] = Ugp[gp][0].vector[r];
                    m_gridFlux[ei][ej].topFlux[r][gp] = Ugp[gp][m_gpNum - 1].vector[r];
                    m_gridFlux[ei][ej].leftFlux[r][gp] = Ugp[0][gp].vector[r];
                    m_gridFlux[ei][ej].rightFlux[r][gp] = Ugp[m_gpNum - 1][gp].vector[r];
                }
            }
        }
    }
}
void CWENOFV::getFlux_primitive(void)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    /////// Implement reconstructing the primitive variables /////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Li, Z., Vittoz, L., Oger, G., & Le Touzé, D. (2022).
    // A simplified and efficient weakly-compressible FV-WENO scheme for immiscible two-phase flows.
    // Computers & Fluids, 244. https://doi.org/10.1016/j.compfluid.2022.105555

    // primitive variables: z1, rho1, U, V

    m_wenoTimer.pause();
    m_mainTimer.pause();

    getWorldUh();
    worldToProcUh();

    m_mainTimer.start();
    m_wenoTimer.start();

    m_totalPr.setZero();

    const int tempNum = 7; // Number of

    Array1D<double> loc_Uh(m_varNum);

    Array2D<Cvector> conservativeVariablesAve(tempNum, tempNum);
    Array2D<Cvector> primitiveVariablesAve(tempNum, tempNum);
    Array2D<Cvector> Ugp(m_gpNum, m_gpNum);
    Array2D<Cvector> Utemple(tempNum, m_gpNum);

    for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
        for (int gp = 0; gp != m_gpNum; gp++)
            Utemple[tempIndex][gp].vector.Resize(m_varNum);

    for (int i = 0; i != tempNum; ++i)
    {
        for (int j = 0; j != tempNum; ++j)
        {
            primitiveVariablesAve[i][j].vector.Resize(m_varNum);
            conservativeVariablesAve[i][j].vector.Resize(m_varNum);
        }
    }

    for (int gpi = 0; gpi != m_gpNum; gpi++)
        for (int gpj = 0; gpj != m_gpNum; gpj++)
            Ugp[gpi][gpj].vector.Resize(m_varNum);

    // use WENO limiter for every cell
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            //(1)Step-1: set stencil
            for (int eii = 0; eii != tempNum; eii++)
                for (int ejj = 0; ejj != tempNum; ejj++)
                    for (int r = 0; r != m_varNum; ++r)
                        conservativeVariablesAve[eii][ejj].vector[r] = m_Uh[ei + eii - 3][ej + ejj - 3].vector[r];

            // Step-2: conservativeVariablesAve to primitiveVariablesAve
            for (int eii = 0; eii != tempNum; eii++)
            {
                for (int ejj = 0; ejj != tempNum; ejj++)
                {
                    loc_Uh = conservativeVariablesAve[eii][ejj].vector;

                    const double z1 = equation->getVolumeFrac(loc_Uh);
                    const double rho1 = loc_Uh[0] / z1;
                    // const double rho1 = 1000.0;
                    const double u = loc_Uh[2] / (loc_Uh[0] + loc_Uh[1]);
                    const double v = loc_Uh[3] / (loc_Uh[0] + loc_Uh[1]);

                    primitiveVariablesAve[eii][ejj].vector[0] = z1;
                    primitiveVariablesAve[eii][ejj].vector[1] = rho1;
                    primitiveVariablesAve[eii][ejj].vector[2] = u;
                    primitiveVariablesAve[eii][ejj].vector[3] = v;
                }
            }

            //(3)Step-3: weno restruction in y-direction
            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        double uavemmm = primitiveVariablesAve[tempIndex][0].vector[r];
                        double uavemm = primitiveVariablesAve[tempIndex][1].vector[r];
                        double uavem = primitiveVariablesAve[tempIndex][2].vector[r];
                        double uavec = primitiveVariablesAve[tempIndex][3].vector[r];
                        double uavep = primitiveVariablesAve[tempIndex][4].vector[r];
                        double uavepp = primitiveVariablesAve[tempIndex][5].vector[r];
                        double uaveppp = primitiveVariablesAve[tempIndex][6].vector[r];

                        Utemple[tempIndex][gp].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp);
                    }
                }
            }

            //(4)Step-4: weno restruction in x-direction
            for (int gpj = 0; gpj != m_gpNum; ++gpj) // y
            {
                for (int gpi = 0; gpi != m_gpNum; ++gpi) // x
                {
                    if (gpi == 0 || gpi == m_gpNum - 1 || gpj == 0 || gpj == m_gpNum - 1)
                    {
                        for (int r = 0; r != m_varNum; ++r)
                        {
                            double uavemmm = Utemple[0][gpj].vector[r];
                            double uavemm = Utemple[1][gpj].vector[r];
                            double uavem = Utemple[2][gpj].vector[r];
                            double uavec = Utemple[3][gpj].vector[r];
                            double uavep = Utemple[4][gpj].vector[r];
                            double uavepp = Utemple[5][gpj].vector[r];
                            double uaveppp = Utemple[6][gpj].vector[r];

                            Ugp[gpi][gpj].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi);
                        }
                    }
                }
            }

            // Step-5: primitiveVariablesAve to conservativeVariablesAve
            for (int gpi = 0; gpi != m_gpNum; gpi++)
            {
                for (int gpj = 0; gpj != m_gpNum; gpj++)
                {
                    const double z1 = Ugp[gpi][gpj].vector[0];
                    const double rho1 = Ugp[gpi][gpj].vector[1];
                    const double u = Ugp[gpi][gpj].vector[2];
                    const double v = Ugp[gpi][gpj].vector[3];

                    // // debug only
                    // const double rho1 = 1000.0;
                    // const double u = 1.0;
                    // const double v = 1.0;

                    double z2rho2 = pow(equation->g_sound1, 2) / pow(equation->g_sound2, 2) * (rho1 - equation->g_rho10) + equation->g_rho20;
                    z2rho2 = z2rho2 * (1.0 - z1);

                    const double z1rho1 = z1 * rho1;
                    const double rho = z1rho1 + z2rho2;

                    conservativeVariablesAve[gpi][gpj].vector[0] = z1rho1;
                    conservativeVariablesAve[gpi][gpj].vector[1] = z2rho2;
                    conservativeVariablesAve[gpi][gpj].vector[2] = rho * u;
                    conservativeVariablesAve[gpi][gpj].vector[3] = rho * v;
                }
            }

            // Step-6: get flux
            for (int r = 0; r != m_varNum; r++)
            {
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    m_gridFlux[ei][ej].bottomFlux[r][gp] = conservativeVariablesAve[gp][0].vector[r];
                    m_gridFlux[ei][ej].topFlux[r][gp] = conservativeVariablesAve[gp][m_gpNum - 1].vector[r];
                    m_gridFlux[ei][ej].leftFlux[r][gp] = conservativeVariablesAve[0][gp].vector[r];
                    m_gridFlux[ei][ej].rightFlux[r][gp] = conservativeVariablesAve[m_gpNum - 1][gp].vector[r];
                }
            }
        }
    }
}
void CWENOFV::getFlux_primitive2(void)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    /////// Implement reconstructing the primitive variables /////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // primitive variables: z1, prssure, U, V
    m_wenoTimer.pause();
    m_mainTimer.pause();

    getWorldUh();
    worldToProcUh();

    m_mainTimer.start();
    m_wenoTimer.start();

    m_totalPr.setZero();

    const int tempNum = 7; // Number of

    Array1D<double> loc_Uh(m_varNum);

    Array2D<Cvector> conservativeVariablesAve(tempNum, tempNum);
    Array2D<Cvector> primitiveVariablesAve(tempNum, tempNum);
    Array2D<Cvector> Ugp(m_gpNum, m_gpNum);
    Array2D<Cvector> Utemple(tempNum, m_gpNum);

    for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
        for (int gp = 0; gp != m_gpNum; gp++)
            Utemple[tempIndex][gp].vector.Resize(m_varNum);

    for (int i = 0; i != tempNum; ++i)
    {
        for (int j = 0; j != tempNum; ++j)
        {
            primitiveVariablesAve[i][j].vector.Resize(m_varNum);
            conservativeVariablesAve[i][j].vector.Resize(m_varNum);
        }
    }

    for (int gpi = 0; gpi != m_gpNum; gpi++)
        for (int gpj = 0; gpj != m_gpNum; gpj++)
            Ugp[gpi][gpj].vector.Resize(m_varNum);

    // use WENO limiter for every cell
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            //(1)Step-1: set stencil
            for (int r = 0; r != m_varNum; ++r)
                for (int eii = 0; eii != tempNum; eii++)
                    for (int ejj = 0; ejj != tempNum; ejj++)
                        conservativeVariablesAve[eii][ejj].vector[r] = m_Uh[ei + eii - 3][ej + ejj - 3].vector[r];

            // Step-2: conservativeVariablesAve to primitiveVariablesAve
            for (int eii = 0; eii != tempNum; eii++)
            {
                for (int ejj = 0; ejj != tempNum; ejj++)
                {
                    loc_Uh = conservativeVariablesAve[eii][ejj].vector;

                    const double z1 = equation->getVolumeFrac(loc_Uh);
                    double p;
                    if (z1 > 0.5)
                        p = equation->g_preb + pow(equation->g_sound1, 2) * (loc_Uh[0] / z1 - equation->g_rho10);
                    else
                        p = equation->g_preb + pow(equation->g_sound2, 2) * (loc_Uh[1] / (1.0 - z1) - equation->g_rho20);

                    const double u = loc_Uh[2] / (loc_Uh[0] + loc_Uh[1]);
                    const double v = loc_Uh[3] / (loc_Uh[0] + loc_Uh[1]);

                    primitiveVariablesAve[eii][ejj].vector[0] = z1;
                    primitiveVariablesAve[eii][ejj].vector[1] = p;
                    primitiveVariablesAve[eii][ejj].vector[2] = u;
                    primitiveVariablesAve[eii][ejj].vector[3] = v;
                }
            }

            //(3)Step-3: weno restruction in y-direction
            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        double uavemmm = primitiveVariablesAve[tempIndex][0].vector[r];
                        double uavemm = primitiveVariablesAve[tempIndex][1].vector[r];
                        double uavem = primitiveVariablesAve[tempIndex][2].vector[r];
                        double uavec = primitiveVariablesAve[tempIndex][3].vector[r];
                        double uavep = primitiveVariablesAve[tempIndex][4].vector[r];
                        double uavepp = primitiveVariablesAve[tempIndex][5].vector[r];
                        double uaveppp = primitiveVariablesAve[tempIndex][6].vector[r];

                        // 方法1
                        Utemple[tempIndex][gp].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp);
                    }
                }
            }

            //(4)Step-4: weno restruction in x-direction
            for (int gpj = 0; gpj != m_gpNum; ++gpj) // y
            {
                for (int gpi = 0; gpi != m_gpNum; ++gpi) // x
                {
                    if (gpi == 0 || gpi == m_gpNum - 1 || gpj == 0 || gpj == m_gpNum - 1)
                    {
                        for (int r = 0; r != m_varNum; ++r)
                        {
                            double uavemmm = Utemple[0][gpj].vector[r];
                            double uavemm = Utemple[1][gpj].vector[r];
                            double uavem = Utemple[2][gpj].vector[r];
                            double uavec = Utemple[3][gpj].vector[r];
                            double uavep = Utemple[4][gpj].vector[r];
                            double uavepp = Utemple[5][gpj].vector[r];
                            double uaveppp = Utemple[6][gpj].vector[r];

                            Ugp[gpi][gpj].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi);
                        }
                    }
                }
            }

            // Step-5: primitiveVariablesAve to conservativeVariablesAve
            for (int gpi = 0; gpi != m_gpNum; gpi++)
            {
                for (int gpj = 0; gpj != m_gpNum; gpj++)
                {
                    const double z1 = Ugp[gpi][gpj].vector[0];
                    const double p = Ugp[gpi][gpj].vector[1];
                    const double u = Ugp[gpi][gpj].vector[2];
                    const double v = Ugp[gpi][gpj].vector[3];

                    const double z1rho1 = ((p - equation->g_preb) / (pow(equation->g_sound1, 2)) + equation->g_rho10) * z1;
                    const double z2rho2 = ((p - equation->g_preb) / (pow(equation->g_sound2, 2)) + equation->g_rho20) * (1.0 - z1);
                    const double rho = z1rho1 + z2rho2;

                    conservativeVariablesAve[gpi][gpj].vector[0] = z1rho1;
                    conservativeVariablesAve[gpi][gpj].vector[1] = z2rho2;
                    conservativeVariablesAve[gpi][gpj].vector[2] = rho * u;
                    conservativeVariablesAve[gpi][gpj].vector[3] = rho * v;
                }
            }

            // Step-6: get flux
            for (int r = 0; r != m_varNum; r++)
            {
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    m_gridFlux[ei][ej].bottomFlux[r][gp] = conservativeVariablesAve[gp][0].vector[r];
                    m_gridFlux[ei][ej].topFlux[r][gp] = conservativeVariablesAve[gp][m_gpNum - 1].vector[r];
                    m_gridFlux[ei][ej].leftFlux[r][gp] = conservativeVariablesAve[0][gp].vector[r];
                    m_gridFlux[ei][ej].rightFlux[r][gp] = conservativeVariablesAve[m_gpNum - 1][gp].vector[r];
                }
            }
        }
    }
}
void CWENOFV::getFlux_primitive3(void)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    /////// Implement reconstructing the primitive variables /////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // primitive variables: z1, prssure, U, V

    m_mainTimer.pause();
    m_wenoTimer.pause();
    getWorldUh();
    worldToProcUh();
    m_wenoTimer.start();
    m_mainTimer.pause();

    m_totalPr.setZero();

    const int tempNum = 7; // Number of

    Array1D<double> loc_Uh(m_varNum);

    Array2D<Cvector> primitiveVariablesAve(m_totalElemNumX, m_totalElemNumY);
    for (int ei = 0; ei != m_totalElemNumX; ei++)
    {
        for (int ej = 0; ej != m_totalElemNumY; ej++)
        {
            primitiveVariablesAve[ei][ej].vector.Resize(m_varNum);

            loc_Uh = m_Uh[ei][ej].vector;

            const double z1 = equation->getVolumeFrac(loc_Uh);
            const double z2 = 1 - z1;
            double p;
            if (z1 > 0.5)
                p = equation->g_preb + pow(equation->g_sound1, 2) * (loc_Uh[0] / z1 - equation->g_rho10);
            else
                p = equation->g_preb + pow(equation->g_sound2, 2) * (loc_Uh[1] / z2 - equation->g_rho20);

            const double u = loc_Uh[2] / (loc_Uh[0] + loc_Uh[1]);
            const double v = loc_Uh[3] / (loc_Uh[0] + loc_Uh[1]);

            primitiveVariablesAve[ei][ej].vector[0] = z1;
            primitiveVariablesAve[ei][ej].vector[1] = p;
            primitiveVariablesAve[ei][ej].vector[2] = u;
            primitiveVariablesAve[ei][ej].vector[3] = v;
        }
    }

    Array2D<Cvector> Utemple(tempNum, m_gpNum);
    for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
        for (int gp = 0; gp != m_gpNum; gp++)
            Utemple[tempIndex][gp].vector.Resize(m_varNum);

    Array2D<Cvector> templeAve(tempNum, tempNum);
    for (int ei = 0; ei != tempNum; ei++)
        for (int ej = 0; ej != tempNum; ej++)
            templeAve[ei][ej].vector.Resize(m_varNum);

    Array2D<Cvector> Ugp(m_gpNum, m_gpNum);
    Array2D<Cvector> conservativeVariablesGP(m_gpNum, m_gpNum);

    for (int gpi = 0; gpi != m_gpNum; gpi++)
    {
        for (int gpj = 0; gpj != m_gpNum; gpj++)
        {
            Ugp[gpi][gpj].vector.Resize(m_varNum);
            conservativeVariablesGP[gpi][gpj].vector.Resize(m_varNum);
        }
    }

    // use WENO limiter for every cell
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            //(1)Step-1: set stencil
            for (int r = 0; r != m_varNum; ++r)
                for (int eii = 0; eii != tempNum; eii++)
                    for (int ejj = 0; ejj != tempNum; ejj++)
                        templeAve[eii][ejj].vector[r] = primitiveVariablesAve[ei + eii - 3][ej + ejj - 3].vector[r];

            //(3)Step-3: weno restruction in y-direction
            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        double uavemmm = templeAve[tempIndex][0].vector[r];
                        double uavemm = templeAve[tempIndex][1].vector[r];
                        double uavem = templeAve[tempIndex][2].vector[r];
                        double uavec = templeAve[tempIndex][3].vector[r];
                        double uavep = templeAve[tempIndex][4].vector[r];
                        double uavepp = templeAve[tempIndex][5].vector[r];
                        double uaveppp = templeAve[tempIndex][6].vector[r];

                        // Utemple[tempIndex][gp].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp);
                        if (r == 0)
                            Utemple[tempIndex][gp].vector[r] = wenoExpWithThincReconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp, m_deltaX);
                        else
                            Utemple[tempIndex][gp].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gp);
                    }
                }
            }

            //(4)Step-4: weno restruction in x-direction
            for (int gpj = 0; gpj != m_gpNum; ++gpj) // y
            {
                for (int gpi = 0; gpi != m_gpNum; ++gpi) // x
                {
                    if (gpi == 0 || gpi == m_gpNum - 1 || gpj == 0 || gpj == m_gpNum - 1)
                    {
                        for (int r = 0; r != m_varNum; ++r)
                        {
                            double uavemmm = Utemple[0][gpj].vector[r];
                            double uavemm = Utemple[1][gpj].vector[r];
                            double uavem = Utemple[2][gpj].vector[r];
                            double uavec = Utemple[3][gpj].vector[r];
                            double uavep = Utemple[4][gpj].vector[r];
                            double uavepp = Utemple[5][gpj].vector[r];
                            double uaveppp = Utemple[6][gpj].vector[r];

                            // Ugp[gpi][gpj].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi);

                            if (r == 0)
                                Ugp[gpi][gpj].vector[r] = wenoExpWithThincReconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi, m_deltaX);
                            else
                                Ugp[gpi][gpj].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gpi);
                        }
                    }
                }
            }

            // Step-5: primitiveVariablesAve to conservativeVariablesAve
            for (int gpi = 0; gpi != m_gpNum; gpi++)
            {
                for (int gpj = 0; gpj != m_gpNum; gpj++)
                {
                    const double z1 = Ugp[gpi][gpj].vector[0];
                    const double z2 = 1 - z1;
                    const double p = Ugp[gpi][gpj].vector[1];
                    const double u = Ugp[gpi][gpj].vector[2];
                    const double v = Ugp[gpi][gpj].vector[3];

                    const double z1rho1 = ((p - equation->g_preb) / (pow(equation->g_sound1, 2)) + equation->g_rho10) * z1;
                    const double z2rho2 = ((p - equation->g_preb) / (pow(equation->g_sound2, 2)) + equation->g_rho20) * z2;
                    const double rho = z1rho1 + z2rho2;

                    conservativeVariablesGP[gpi][gpj].vector[0] = z1rho1;
                    conservativeVariablesGP[gpi][gpj].vector[1] = z2rho2;
                    conservativeVariablesGP[gpi][gpj].vector[2] = rho * u;
                    conservativeVariablesGP[gpi][gpj].vector[3] = rho * v;
                }
            }

            // Step-6: get flux
            for (int r = 0; r != m_varNum; r++)
            {
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    m_gridFlux[ei][ej].bottomFlux[r][gp] = conservativeVariablesGP[gp][0].vector[r];
                    m_gridFlux[ei][ej].topFlux[r][gp] = conservativeVariablesGP[gp][m_gpNum - 1].vector[r];
                    m_gridFlux[ei][ej].leftFlux[r][gp] = conservativeVariablesGP[0][gp].vector[r];
                    m_gridFlux[ei][ej].rightFlux[r][gp] = conservativeVariablesGP[m_gpNum - 1][gp].vector[r];
                }
            }
        }
    }
}
void CWENOFV::getFlux_ThincOnlyForZ1(void)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    /////// Implement reconstructing the primitive variables /////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // primitive variables: z1, prssure, U, V

    m_mainTimer.pause();
    m_wenoTimer.pause();
    getWorldUh();
    worldToProcUh();
    m_wenoTimer.start();
    m_mainTimer.pause();

    m_totalPr.setZero();

    const int tempNum = 7; // Number of

    Array1D<double> loc_Uh(m_varNum);

    Array2D<Cvector> primitiveVariablesAve(m_totalElemNumX, m_totalElemNumY);
    for (int ei = 0; ei != m_totalElemNumX; ei++)
    {
        for (int ej = 0; ej != m_totalElemNumY; ej++)
        {
            primitiveVariablesAve[ei][ej].vector.Resize(m_varNum);

            loc_Uh = m_Uh[ei][ej].vector;

            const double z1 = equation->getVolumeFrac(loc_Uh);
            const double z2 = 1 - z1;
            double p;
            if (z1 > 0.5)
                p = equation->g_preb + pow(equation->g_sound1, 2) * (loc_Uh[0] / z1 - equation->g_rho10);
            else
                p = equation->g_preb + pow(equation->g_sound2, 2) * (loc_Uh[1] / z2 - equation->g_rho20);

            const double u = loc_Uh[2] / (loc_Uh[0] + loc_Uh[1]);
            const double v = loc_Uh[3] / (loc_Uh[0] + loc_Uh[1]);

            primitiveVariablesAve[ei][ej].vector[0] = z1;
            primitiveVariablesAve[ei][ej].vector[1] = p;
            primitiveVariablesAve[ei][ej].vector[2] = u;
            primitiveVariablesAve[ei][ej].vector[3] = v;
        }
    }

    Array2D<Cvector> Utemple(tempNum, m_gpNum);
    for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
        for (int gp = 0; gp != m_gpNum; gp++)
            Utemple[tempIndex][gp].vector.Resize(m_varNum);
    Array2D<Cvector> UtempleThinc(tempNum, m_gpNum);
    for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
        for (int gp = 0; gp != m_gpNum; gp++)
            UtempleThinc[tempIndex][gp].vector.Resize(m_varNum);

    Array2D<Cvector> templeAve(tempNum, tempNum);
    for (int ei = 0; ei != tempNum; ei++)
        for (int ej = 0; ej != tempNum; ej++)
            templeAve[ei][ej].vector.Resize(m_varNum);

    Array2D<Cvector> Ugp(m_gpNum, m_gpNum);
    Array2D<Cvector> conservativeVariablesGP(m_gpNum, m_gpNum);
    for (int gpi = 0; gpi != m_gpNum; gpi++)
    {
        for (int gpj = 0; gpj != m_gpNum; gpj++)
        {
            Ugp[gpi][gpj].vector.Resize(m_varNum);
            conservativeVariablesGP[gpi][gpj].vector.Resize(m_varNum);
        }
    }

    Array2D<Cvector> UgpThinc(m_gpNum, m_gpNum);
    Array2D<Cvector> conservativeVariablesGPThinc(m_gpNum, m_gpNum);
    for (int gpi = 0; gpi != m_gpNum; gpi++)
    {
        for (int gpj = 0; gpj != m_gpNum; gpj++)
        {
            UgpThinc[gpi][gpj].vector.Resize(m_varNum);
            conservativeVariablesGPThinc[gpi][gpj].vector.Resize(m_varNum);
        }
    }

    // use WENO limiter for every cell
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            //(1)Step-1: set stencil
            for (int r = 0; r != m_varNum; ++r)
                for (int eii = 0; eii != tempNum; eii++)
                    for (int ejj = 0; ejj != tempNum; ejj++)
                        templeAve[eii][ejj].vector[r] = primitiveVariablesAve[ei + eii - 3][ej + ejj - 3].vector[r];

            //(3)Step-3: weno restruction in y-direction
            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        double uavemmm = templeAve[tempIndex][0].vector[r];
                        double uavemm = templeAve[tempIndex][1].vector[r];
                        double uavem = templeAve[tempIndex][2].vector[r];
                        double uavec = templeAve[tempIndex][3].vector[r];
                        double uavep = templeAve[tempIndex][4].vector[r];
                        double uavepp = templeAve[tempIndex][5].vector[r];
                        double uaveppp = templeAve[tempIndex][6].vector[r];

                        // Utemple[tempIndex][gp].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp);
                        if (r == 0)
                        {
                            Utemple[tempIndex][gp].vector[r] = wenoExpWithThincReconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp, m_deltaX);
                            UtempleThinc[tempIndex][gp].vector[r] = thincReconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp, m_deltaX);
                        }
                        else
                        {
                            Utemple[tempIndex][gp].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gp);
                            UtempleThinc[tempIndex][gp].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gp);
                        }
                    }
                }
            }

            //(4)Step-4: weno restruction in x-direction
            for (int gpj = 0; gpj != m_gpNum; ++gpj) // y
            {
                for (int gpi = 0; gpi != m_gpNum; ++gpi) // x
                {
                    if (gpi == 0 || gpi == m_gpNum - 1 || gpj == 0 || gpj == m_gpNum - 1)
                    {
                        for (int r = 0; r != m_varNum; ++r)
                        {
                            double uavemmm = Utemple[0][gpj].vector[r];
                            double uavemm = Utemple[1][gpj].vector[r];
                            double uavem = Utemple[2][gpj].vector[r];
                            double uavec = Utemple[3][gpj].vector[r];
                            double uavep = Utemple[4][gpj].vector[r];
                            double uavepp = Utemple[5][gpj].vector[r];
                            double uaveppp = Utemple[6][gpj].vector[r];

                            // Ugp[gpi][gpj].vector[r] = useWENO(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi);

                            if (r == 0)
                                Ugp[gpi][gpj].vector[r] = wenoExpWithThincReconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi, m_deltaX);
                            else
                                Ugp[gpi][gpj].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gpi);
                        }

                        for (int r = 0; r != m_varNum; ++r)
                        {
                            double uavemmm = UtempleThinc[0][gpj].vector[r];
                            double uavemm = UtempleThinc[1][gpj].vector[r];
                            double uavem = UtempleThinc[2][gpj].vector[r];
                            double uavec = UtempleThinc[3][gpj].vector[r];
                            double uavep = UtempleThinc[4][gpj].vector[r];
                            double uavepp = UtempleThinc[5][gpj].vector[r];
                            double uaveppp = UtempleThinc[6][gpj].vector[r];

                            if (r == 0)
                                UgpThinc[gpi][gpj].vector[r] = wenoExpWithThincReconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi, m_deltaX);
                            else
                                UgpThinc[gpi][gpj].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gpi);
                        }
                    }
                }
            }

            // Step-5: primitiveVariablesAve to conservativeVariablesAve
            for (int gpi = 0; gpi != m_gpNum; gpi++)
            {
                for (int gpj = 0; gpj != m_gpNum; gpj++)
                {
                    const double z1 = Ugp[gpi][gpj].vector[0];
                    const double z2 = 1 - z1;
                    const double p = Ugp[gpi][gpj].vector[1];
                    const double u = Ugp[gpi][gpj].vector[2];
                    const double v = Ugp[gpi][gpj].vector[3];

                    const double z1rho1 = ((p - equation->g_preb) / (pow(equation->g_sound1, 2)) + equation->g_rho10) * z1;
                    const double z2rho2 = ((p - equation->g_preb) / (pow(equation->g_sound2, 2)) + equation->g_rho20) * z2;
                    const double rho = z1rho1 + z2rho2;

                    conservativeVariablesGP[gpi][gpj].vector[0] = z1rho1;
                    conservativeVariablesGP[gpi][gpj].vector[1] = z2rho2;
                    conservativeVariablesGP[gpi][gpj].vector[2] = rho * u;
                    conservativeVariablesGP[gpi][gpj].vector[3] = rho * v;
                }
            }

            // Step-6: get flux
            for (int r = 0; r != m_varNum; r++)
            {
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    m_gridFlux[ei][ej].bottomFlux[r][gp] = conservativeVariablesGP[gp][0].vector[r];
                    m_gridFlux[ei][ej].topFlux[r][gp] = conservativeVariablesGP[gp][m_gpNum - 1].vector[r];
                    m_gridFlux[ei][ej].leftFlux[r][gp] = conservativeVariablesGP[0][gp].vector[r];
                    m_gridFlux[ei][ej].rightFlux[r][gp] = conservativeVariablesGP[m_gpNum - 1][gp].vector[r];
                }
            }
        }
    }
}

void CWENOFV::getFlux_WENOEXPandTHINC(void)
{
    // （step1) 首先用 WENOEXP 算一版
    Array2D<CgridFlux> gridFluxForWENOEXP(m_totalElemNumX, m_totalElemNumY);
    for (int ei = 0; ei < m_totalElemNumX; ++ei)
    {
        for (int ej = 0; ej < m_totalElemNumY; ++ej)
        {
            gridFluxForWENOEXP[ei][ej].bottomFlux.Resize(m_varNum, m_gpNum);
            gridFluxForWENOEXP[ei][ej].topFlux.Resize(m_varNum, m_gpNum);
            gridFluxForWENOEXP[ei][ej].leftFlux.Resize(m_varNum, m_gpNum);
            gridFluxForWENOEXP[ei][ej].rightFlux.Resize(m_varNum, m_gpNum);
        }
    }
    getFlux_WENOEXPandTHINCzhiWENOEXP(gridFluxForWENOEXP);

    // (step2) 然后用 THINC 算一版
    Array2D<CgridFlux> gridFluxForTHINC(m_totalElemNumX, m_totalElemNumY);
    for (int ei = 0; ei < m_totalElemNumX; ++ei)
    {
        for (int ej = 0; ej < m_totalElemNumY; ++ej)
        {
            gridFluxForTHINC[ei][ej].bottomFlux.Resize(m_varNum, m_gpNum);
            gridFluxForTHINC[ei][ej].topFlux.Resize(m_varNum, m_gpNum);
            gridFluxForTHINC[ei][ej].leftFlux.Resize(m_varNum, m_gpNum);
            gridFluxForTHINC[ei][ej].rightFlux.Resize(m_varNum, m_gpNum);
        }
    }
    getFlux_WENOEXPandTHINCzhiTHINC(gridFluxForTHINC);

    // 通讯，更新边界，最临近一层的 ghost 单元的高斯点值需要同步
    getFlux_WENOEXPandTHINCzhiMPI(gridFluxForWENOEXP);
    getFlux_WENOEXPandTHINCzhiMPI(gridFluxForTHINC);

    // (step3) 使用 TVB 算法选择最恰当的一个

    // 先默认使用WENOEXP
    for (int ei = 0; ei < m_totalElemNumX; ++ei)
    {
        for (int ej = 0; ej < m_totalElemNumY; ++ej)
        {
            for (int r = 0; r != m_varNum; r++)
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    m_gridFlux[ei][ej].bottomFlux[r][gp] = gridFluxForWENOEXP[ei][ej].bottomFlux[r][gp];
                    m_gridFlux[ei][ej].topFlux[r][gp] = gridFluxForWENOEXP[ei][ej].topFlux[r][gp];
                    m_gridFlux[ei][ej].leftFlux[r][gp] = gridFluxForWENOEXP[ei][ej].leftFlux[r][gp];
                    m_gridFlux[ei][ej].rightFlux[r][gp] = gridFluxForWENOEXP[ei][ej].rightFlux[r][gp];
                }
        }
    }

    // 首先计算 TVB 值
    Array2D<double> tvbForWENOEXP(m_totalElemNumX, m_totalElemNumY);
    Array2D<double> tvbForTHINC(m_totalElemNumX, m_totalElemNumY);
    tvbForWENOEXP.setZero();
    tvbForTHINC.setZero();
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            for (int gp = 0; gp != m_gpNum; gp++)
            {
                for (int r = 0; r != m_varNum; r++)
                {
                    tvbForWENOEXP[ei][ej] += fabs(gridFluxForWENOEXP[ei][ej].rightFlux[r][gp] - gridFluxForWENOEXP[ei + 1][ej].leftFlux[r][gp]);
                    tvbForWENOEXP[ei][ej] += fabs(gridFluxForWENOEXP[ei][ej].topFlux[r][gp] - gridFluxForWENOEXP[ei][ej + 1].bottomFlux[r][gp]);
                    tvbForWENOEXP[ei][ej] += fabs(gridFluxForWENOEXP[ei][ej].leftFlux[r][gp] - gridFluxForWENOEXP[ei - 1][ej].rightFlux[r][gp]);
                    tvbForWENOEXP[ei][ej] += fabs(gridFluxForWENOEXP[ei][ej].bottomFlux[r][gp] - gridFluxForWENOEXP[ei][ej - 1].topFlux[r][gp]);

                    tvbForTHINC[ei][ej] += fabs(gridFluxForTHINC[ei][ej].rightFlux[r][gp] - gridFluxForTHINC[ei + 1][ej].leftFlux[r][gp]);
                    tvbForTHINC[ei][ej] += fabs(gridFluxForTHINC[ei][ej].topFlux[r][gp] - gridFluxForTHINC[ei][ej + 1].bottomFlux[r][gp]);
                    tvbForTHINC[ei][ej] += fabs(gridFluxForTHINC[ei][ej].leftFlux[r][gp] - gridFluxForTHINC[ei - 1][ej].rightFlux[r][gp]);
                    tvbForTHINC[ei][ej] += fabs(gridFluxForTHINC[ei][ej].bottomFlux[r][gp] - gridFluxForTHINC[ei][ej - 1].topFlux[r][gp]);
                }
            }
        }
    }

    // 如果THINC的TVB值小于WENOEXP的TVB值，则周围的一圈的cell都做好标记
    Array2D<double> kaiguan(m_totalElemNumX, m_totalElemNumY);
    kaiguan.setZero();

    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            if (tvbForTHINC[ei][ej] < tvbForWENOEXP[ei][ej])
            {
                kaiguan[ei][ej] = 1;

                // kaiguan[ei + 1][ej] = 1;
                // kaiguan[ei - 1][ej] = 1;

                // kaiguan[ei][ej + 1] = 1;
                // kaiguan[ei][ej - 1] = 1;

                // kaiguan[ei + 1][ej + 1] = 1;
                // kaiguan[ei - 1][ej - 1] = 1;

                // kaiguan[ei + 1][ej - 1] = 1;
                // kaiguan[ei - 1][ej + 1] = 1;
            }
        }
    }

    // 对于被标记的cell, 使用thinc
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            if (kaiguan[ei][ej] == 1)
            {
                for (int r = 0; r != m_varNum; r++)
                    for (int gp = 0; gp != m_gpNum; gp++)
                    {
                        m_gridFlux[ei][ej].bottomFlux[r][gp] = gridFluxForTHINC[ei][ej].bottomFlux[r][gp];
                        m_gridFlux[ei][ej].topFlux[r][gp] = gridFluxForTHINC[ei][ej].topFlux[r][gp];
                        m_gridFlux[ei][ej].leftFlux[r][gp] = gridFluxForTHINC[ei][ej].leftFlux[r][gp];
                        m_gridFlux[ei][ej].rightFlux[r][gp] = gridFluxForTHINC[ei][ej].rightFlux[r][gp];
                    }
            }
        }
    }
}
void CWENOFV::getFlux_WENOEXPandTHINCzhiWENOEXP(Array2D<CgridFlux> &gridFluxForWENOEXP)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    /////// Implement reconstructing the primitive variables /////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // primitive variables: z1, prssure, U, V
    m_wenoTimer.pause();
    m_mainTimer.pause();

    getWorldUh();
    worldToProcUh();

    m_mainTimer.start();
    m_wenoTimer.start();

    m_totalPr.setZero();

    const int tempNum = 7; // Number of

    Array1D<double> loc_Uh(m_varNum);

    Array2D<Cvector> conservativeVariablesAve(tempNum, tempNum);
    Array2D<Cvector> primitiveVariablesAve(tempNum, tempNum);
    Array2D<Cvector> Ugp(m_gpNum, m_gpNum);
    Array2D<Cvector> Utemple(tempNum, m_gpNum);

    for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
        for (int gp = 0; gp != m_gpNum; gp++)
            Utemple[tempIndex][gp].vector.Resize(m_varNum);

    for (int i = 0; i != tempNum; ++i)
    {
        for (int j = 0; j != tempNum; ++j)
        {
            primitiveVariablesAve[i][j].vector.Resize(m_varNum);
            conservativeVariablesAve[i][j].vector.Resize(m_varNum);
        }
    }

    for (int gpi = 0; gpi != m_gpNum; gpi++)
        for (int gpj = 0; gpj != m_gpNum; gpj++)
            Ugp[gpi][gpj].vector.Resize(m_varNum);

    // use WENO limiter for every cell
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            //(1)Step-1: set stencil
            for (int r = 0; r != m_varNum; ++r)
                for (int eii = 0; eii != tempNum; eii++)
                    for (int ejj = 0; ejj != tempNum; ejj++)
                        conservativeVariablesAve[eii][ejj].vector[r] = m_Uh[ei + eii - 3][ej + ejj - 3].vector[r];

            // Step-2: conservativeVariablesAve to primitiveVariablesAve
            for (int eii = 0; eii != tempNum; eii++)
            {
                for (int ejj = 0; ejj != tempNum; ejj++)
                {
                    loc_Uh = conservativeVariablesAve[eii][ejj].vector;

                    const double z1 = equation->getVolumeFrac(loc_Uh);
                    double p;
                    if (z1 > 0.5)
                        p = equation->g_preb + pow(equation->g_sound1, 2) * (loc_Uh[0] / z1 - equation->g_rho10);
                    else
                        p = equation->g_preb + pow(equation->g_sound2, 2) * (loc_Uh[1] / (1.0 - z1) - equation->g_rho20);

                    const double u = loc_Uh[2] / (loc_Uh[0] + loc_Uh[1]);
                    const double v = loc_Uh[3] / (loc_Uh[0] + loc_Uh[1]);

                    primitiveVariablesAve[eii][ejj].vector[0] = z1;
                    primitiveVariablesAve[eii][ejj].vector[1] = p;
                    primitiveVariablesAve[eii][ejj].vector[2] = u;
                    primitiveVariablesAve[eii][ejj].vector[3] = v;
                }
            }

            //(3)Step-3: weno restruction in y-direction
            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        double uavemmm = primitiveVariablesAve[tempIndex][0].vector[r];
                        double uavemm = primitiveVariablesAve[tempIndex][1].vector[r];
                        double uavem = primitiveVariablesAve[tempIndex][2].vector[r];
                        double uavec = primitiveVariablesAve[tempIndex][3].vector[r];
                        double uavep = primitiveVariablesAve[tempIndex][4].vector[r];
                        double uavepp = primitiveVariablesAve[tempIndex][5].vector[r];
                        double uaveppp = primitiveVariablesAve[tempIndex][6].vector[r];

                        // 方法1
                        if (r == 0)
                            Utemple[tempIndex][gp].vector[r] = wenoExpthreconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp, m_deltaX);
                        // debug
                        // Utemple[tempIndex][gp].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gp);

                        else
                            Utemple[tempIndex][gp].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gp);
                    }
                }
            }

            //(4)Step-4: weno restruction in x-direction
            for (int gpj = 0; gpj != m_gpNum; ++gpj) // y
            {
                for (int gpi = 0; gpi != m_gpNum; ++gpi) // x
                {
                    if (gpi == 0 || gpi == m_gpNum - 1 || gpj == 0 || gpj == m_gpNum - 1)
                    {
                        for (int r = 0; r != m_varNum; ++r)
                        {
                            double uavemmm = Utemple[0][gpj].vector[r];
                            double uavemm = Utemple[1][gpj].vector[r];
                            double uavem = Utemple[2][gpj].vector[r];
                            double uavec = Utemple[3][gpj].vector[r];
                            double uavep = Utemple[4][gpj].vector[r];
                            double uavepp = Utemple[5][gpj].vector[r];
                            double uaveppp = Utemple[6][gpj].vector[r];
                            if (r == 0)
                                Ugp[gpi][gpj].vector[r] = wenoExpthreconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi, m_deltaX);
                            // debug
                            // Ugp[gpi][gpj].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gpi);
                            else
                                Ugp[gpi][gpj].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gpi);
                        }
                    }
                }
            }

            // Step-5: primitiveVariablesAve to conservativeVariablesAve
            for (int gpi = 0; gpi != m_gpNum; gpi++)
            {
                for (int gpj = 0; gpj != m_gpNum; gpj++)
                {
                    const double z1 = Ugp[gpi][gpj].vector[0];
                    const double p = Ugp[gpi][gpj].vector[1];
                    const double u = Ugp[gpi][gpj].vector[2];
                    const double v = Ugp[gpi][gpj].vector[3];

                    const double z1rho1 = ((p - equation->g_preb) / (pow(equation->g_sound1, 2)) + equation->g_rho10) * z1;
                    const double z2rho2 = ((p - equation->g_preb) / (pow(equation->g_sound2, 2)) + equation->g_rho20) * (1.0 - z1);
                    const double rho = z1rho1 + z2rho2;

                    conservativeVariablesAve[gpi][gpj].vector[0] = z1rho1;
                    conservativeVariablesAve[gpi][gpj].vector[1] = z2rho2;
                    conservativeVariablesAve[gpi][gpj].vector[2] = rho * u;
                    conservativeVariablesAve[gpi][gpj].vector[3] = rho * v;
                }
            }

            // Step-6: get flux
            for (int r = 0; r != m_varNum; r++)
            {
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    gridFluxForWENOEXP[ei][ej].bottomFlux[r][gp] = conservativeVariablesAve[gp][0].vector[r];
                    gridFluxForWENOEXP[ei][ej].topFlux[r][gp] = conservativeVariablesAve[gp][m_gpNum - 1].vector[r];
                    gridFluxForWENOEXP[ei][ej].leftFlux[r][gp] = conservativeVariablesAve[0][gp].vector[r];
                    gridFluxForWENOEXP[ei][ej].rightFlux[r][gp] = conservativeVariablesAve[m_gpNum - 1][gp].vector[r];
                }
            }
        }
    }
}
void CWENOFV::getFlux_WENOEXPandTHINCzhiTHINC(Array2D<CgridFlux> &gridFluxForTHINC)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    /////// Implement reconstructing the primitive variables /////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // primitive variables: z1, prssure, U, V
    m_wenoTimer.pause();
    m_mainTimer.pause();

    getWorldUh();
    worldToProcUh();

    m_mainTimer.start();
    m_wenoTimer.start();

    m_totalPr.setZero();

    const int tempNum = 7; // Number of

    Array1D<double> loc_Uh(m_varNum);

    Array2D<Cvector> conservativeVariablesAve(tempNum, tempNum);
    Array2D<Cvector> primitiveVariablesAve(tempNum, tempNum);
    Array2D<Cvector> Ugp(m_gpNum, m_gpNum);
    Array2D<Cvector> Utemple(tempNum, m_gpNum);

    for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
        for (int gp = 0; gp != m_gpNum; gp++)
            Utemple[tempIndex][gp].vector.Resize(m_varNum);

    for (int i = 0; i != tempNum; ++i)
    {
        for (int j = 0; j != tempNum; ++j)
        {
            primitiveVariablesAve[i][j].vector.Resize(m_varNum);
            conservativeVariablesAve[i][j].vector.Resize(m_varNum);
        }
    }

    for (int gpi = 0; gpi != m_gpNum; gpi++)
        for (int gpj = 0; gpj != m_gpNum; gpj++)
            Ugp[gpi][gpj].vector.Resize(m_varNum);

    // use WENO limiter for every cell
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            //(1)Step-1: set stencil
            for (int r = 0; r != m_varNum; ++r)
                for (int eii = 0; eii != tempNum; eii++)
                    for (int ejj = 0; ejj != tempNum; ejj++)
                        conservativeVariablesAve[eii][ejj].vector[r] = m_Uh[ei + eii - 3][ej + ejj - 3].vector[r];

            // Step-2: conservativeVariablesAve to primitiveVariablesAve
            for (int eii = 0; eii != tempNum; eii++)
            {
                for (int ejj = 0; ejj != tempNum; ejj++)
                {
                    loc_Uh = conservativeVariablesAve[eii][ejj].vector;

                    const double z1 = equation->getVolumeFrac(loc_Uh);
                    double p;
                    if (z1 > 0.5)
                        p = equation->g_preb + pow(equation->g_sound1, 2) * (loc_Uh[0] / z1 - equation->g_rho10);
                    else
                        p = equation->g_preb + pow(equation->g_sound2, 2) * (loc_Uh[1] / (1.0 - z1) - equation->g_rho20);

                    const double u = loc_Uh[2] / (loc_Uh[0] + loc_Uh[1]);
                    const double v = loc_Uh[3] / (loc_Uh[0] + loc_Uh[1]);

                    primitiveVariablesAve[eii][ejj].vector[0] = z1;
                    primitiveVariablesAve[eii][ejj].vector[1] = p;
                    primitiveVariablesAve[eii][ejj].vector[2] = u;
                    primitiveVariablesAve[eii][ejj].vector[3] = v;
                }
            }

            //(3)Step-3: weno restruction in y-direction
            for (int tempIndex = 0; tempIndex != tempNum; tempIndex++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        double uavemmm = primitiveVariablesAve[tempIndex][0].vector[r];
                        double uavemm = primitiveVariablesAve[tempIndex][1].vector[r];
                        double uavem = primitiveVariablesAve[tempIndex][2].vector[r];
                        double uavec = primitiveVariablesAve[tempIndex][3].vector[r];
                        double uavep = primitiveVariablesAve[tempIndex][4].vector[r];
                        double uavepp = primitiveVariablesAve[tempIndex][5].vector[r];
                        double uaveppp = primitiveVariablesAve[tempIndex][6].vector[r];

                        // 方法1
                        if (r == 0)
                            Utemple[tempIndex][gp].vector[r] = thincReconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gp, m_deltaX);
                        else
                            Utemple[tempIndex][gp].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gp);
                    }
                }
            }

            //(4)Step-4: weno restruction in x-direction
            for (int gpj = 0; gpj != m_gpNum; ++gpj) // y
            {
                for (int gpi = 0; gpi != m_gpNum; ++gpi) // x
                {
                    if (gpi == 0 || gpi == m_gpNum - 1 || gpj == 0 || gpj == m_gpNum - 1)
                    {
                        for (int r = 0; r != m_varNum; ++r)
                        {
                            double uavemmm = Utemple[0][gpj].vector[r];
                            double uavemm = Utemple[1][gpj].vector[r];
                            double uavem = Utemple[2][gpj].vector[r];
                            double uavec = Utemple[3][gpj].vector[r];
                            double uavep = Utemple[4][gpj].vector[r];
                            double uavepp = Utemple[5][gpj].vector[r];
                            double uaveppp = Utemple[6][gpj].vector[r];
                            if (r == 0)
                                Ugp[gpi][gpj].vector[r] = thincReconstruction(uavemmm, uavemm, uavem, uavec, uavep, uavepp, uaveppp, gpi, m_deltaX);
                            else
                                Ugp[gpi][gpj].vector[r] = WENO5Zthreconstruction(uavemm, uavem, uavec, uavep, uavepp, gpi);
                        }
                    }
                }
            }

            // Step-5: primitiveVariablesAve to conservativeVariablesAve
            for (int gpi = 0; gpi != m_gpNum; gpi++)
            {
                for (int gpj = 0; gpj != m_gpNum; gpj++)
                {
                    const double z1 = Ugp[gpi][gpj].vector[0];
                    const double p = Ugp[gpi][gpj].vector[1];
                    const double u = Ugp[gpi][gpj].vector[2];
                    const double v = Ugp[gpi][gpj].vector[3];

                    const double z1rho1 = ((p - equation->g_preb) / (pow(equation->g_sound1, 2)) + equation->g_rho10) * z1;
                    const double z2rho2 = ((p - equation->g_preb) / (pow(equation->g_sound2, 2)) + equation->g_rho20) * (1.0 - z1);
                    const double rho = z1rho1 + z2rho2;

                    conservativeVariablesAve[gpi][gpj].vector[0] = z1rho1;
                    conservativeVariablesAve[gpi][gpj].vector[1] = z2rho2;
                    conservativeVariablesAve[gpi][gpj].vector[2] = rho * u;
                    conservativeVariablesAve[gpi][gpj].vector[3] = rho * v;
                }
            }

            // Step-6: get flux
            for (int r = 0; r != m_varNum; r++)
            {
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    gridFluxForTHINC[ei][ej].bottomFlux[r][gp] = conservativeVariablesAve[gp][0].vector[r];
                    gridFluxForTHINC[ei][ej].topFlux[r][gp] = conservativeVariablesAve[gp][m_gpNum - 1].vector[r];
                    gridFluxForTHINC[ei][ej].leftFlux[r][gp] = conservativeVariablesAve[0][gp].vector[r];
                    gridFluxForTHINC[ei][ej].rightFlux[r][gp] = conservativeVariablesAve[m_gpNum - 1][gp].vector[r];
                }
            }
        }
    }
}
void CWENOFV::getFlux_WENOEXPandTHINCzhiMPI(Array2D<CgridFlux> &gridFluxForReconstruction)
{
    Array2D<CgridFlux> worldGridFluxForReconstruction(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);

    if (m_rank == 0)
    {
        for (int ei = 0; ei < m_worldElemNumX + 2 * m_ghostCellNum; ++ei)
        {
            for (int ej = 0; ej < m_worldElemNumY + 2 * m_ghostCellNum; ++ej)
            {
                worldGridFluxForReconstruction[ei][ej].bottomFlux.Resize(m_varNum, m_gpNum);
                worldGridFluxForReconstruction[ei][ej].topFlux.Resize(m_varNum, m_gpNum);
                worldGridFluxForReconstruction[ei][ej].leftFlux.Resize(m_varNum, m_gpNum);
                worldGridFluxForReconstruction[ei][ej].rightFlux.Resize(m_varNum, m_gpNum);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // Receive data from other ranks
        for (int src_rank = 1; src_rank < m_size; ++src_rank)
        {
            MPI_Status status;
            std::vector<double> recvBuffer(m_elemNumX * m_elemNumY * m_varNum * 4 * m_gpNum); // 4 fluxes per grid element
            MPI_Recv(recvBuffer.data(), m_elemNumX * m_elemNumY * m_varNum * 4 * m_gpNum, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &status);

            // Calculate the position in worldGridFluxForReconstruction
            int offsetX = (src_rank % m_numProcessorsX) * m_elemNumX;
            int offsetY = (src_rank / m_numProcessorsX) * m_elemNumY;

            // Store received data in worldGridFluxForReconstruction
            for (int ei = 0; ei < m_elemNumX; ++ei)
            {
                for (int ej = 0; ej < m_elemNumY; ++ej)
                {
                    for (int r = 0; r < m_varNum; ++r)
                    {
                        for (int gp = 0; gp < m_gpNum; ++gp)
                        {
                            int index = ((ei * m_elemNumY + ej) * m_varNum + r) * m_gpNum + gp;
                            worldGridFluxForReconstruction[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum].bottomFlux[r][gp] = recvBuffer[index];
                            worldGridFluxForReconstruction[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum].topFlux[r][gp] = recvBuffer[index + m_elemNumX * m_elemNumY * m_varNum * m_gpNum];
                            worldGridFluxForReconstruction[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum].leftFlux[r][gp] = recvBuffer[index + 2 * m_elemNumX * m_elemNumY * m_varNum * m_gpNum];
                            worldGridFluxForReconstruction[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum].rightFlux[r][gp] = recvBuffer[index + 3 * m_elemNumX * m_elemNumY * m_varNum * m_gpNum];
                        }
                    }
                }
            }
        }

        // Process data for rank 0 itself
        for (int ei = 0; ei < m_elemNumX; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    for (int gp = 0; gp < m_gpNum; ++gp)
                    {
                        worldGridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].bottomFlux[r][gp] = gridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].bottomFlux[r][gp];
                        worldGridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].topFlux[r][gp] = gridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].topFlux[r][gp];
                        worldGridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].leftFlux[r][gp] = gridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].leftFlux[r][gp];
                        worldGridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].rightFlux[r][gp] = gridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].rightFlux[r][gp];
                    }
                }
            }
        }
    }
    else
    {
        // Serialize gridFluxForReconstruction data to a one-dimensional array
        std::vector<double> sendBuffer(m_elemNumX * m_elemNumY * m_varNum * 4 * m_gpNum); // 4 fluxes per grid element
        for (int ei = 0; ei < m_elemNumX; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    for (int gp = 0; gp < m_gpNum; ++gp)
                    {
                        int index = ((ei * m_elemNumY + ej) * m_varNum + r) * m_gpNum + gp;
                        sendBuffer[index] = gridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].bottomFlux[r][gp];
                        sendBuffer[index + m_elemNumX * m_elemNumY * m_varNum * m_gpNum] = gridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].topFlux[r][gp];
                        sendBuffer[index + 2 * m_elemNumX * m_elemNumY * m_varNum * m_gpNum] = gridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].leftFlux[r][gp];
                        sendBuffer[index + 3 * m_elemNumX * m_elemNumY * m_varNum * m_gpNum] = gridFluxForReconstruction[ei + m_ghostCellNum][ej + m_ghostCellNum].rightFlux[r][gp];
                    }
                }
            }
        }

        // Send data to rank 0
        MPI_Send(sendBuffer.data(), m_elemNumX * m_elemNumY * m_varNum * 4 * m_gpNum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    getFlux_WENOEXPandTHINCzhisetBoundaryGPs(worldGridFluxForReconstruction);

    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // Rank 0 sends data to other ranks
        for (int dest_rank = 1; dest_rank < m_size; ++dest_rank)
        {
            const int offsetX = (dest_rank % m_numProcessorsX) * m_elemNumX;
            const int offsetY = (dest_rank / m_numProcessorsX) * m_elemNumY;

            std::vector<double> sendBuffer((m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * 4 * m_gpNum); // 4 fluxes per grid element
            for (int ei = 0; ei < m_elemNumX + 2 * m_ghostCellNum; ++ei)
            {
                for (int ej = 0; ej < m_elemNumY + 2 * m_ghostCellNum; ++ej)
                {
                    for (int r = 0; r < m_varNum; ++r)
                    {
                        for (int gp = 0; gp < m_gpNum; ++gp)
                        {
                            const int index = ((ei * (m_elemNumY + 2 * m_ghostCellNum) + ej) * m_varNum + r) * m_gpNum + gp;
                            sendBuffer[index] = worldGridFluxForReconstruction[ei + offsetX][ej + offsetY].bottomFlux[r][gp];
                            sendBuffer[index + (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum] = worldGridFluxForReconstruction[ei + offsetX][ej + offsetY].topFlux[r][gp];
                            sendBuffer[index + 2 * (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum] = worldGridFluxForReconstruction[ei + offsetX][ej + offsetY].leftFlux[r][gp];
                            sendBuffer[index + 3 * (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum] = worldGridFluxForReconstruction[ei + offsetX][ej + offsetY].rightFlux[r][gp];
                        }
                    }
                }
            }

            // Send data to dest_rank
            MPI_Send(sendBuffer.data(), (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * 4 * m_gpNum, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
        }

        // Process data for rank 0 itself
        for (int ei = 0; ei < m_elemNumX + 2 * m_ghostCellNum; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY + 2 * m_ghostCellNum; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    for (int gp = 0; gp < m_gpNum; ++gp)
                    {
                        gridFluxForReconstruction[ei][ej].bottomFlux[r][gp] = worldGridFluxForReconstruction[ei][ej].bottomFlux[r][gp];
                        gridFluxForReconstruction[ei][ej].topFlux[r][gp] = worldGridFluxForReconstruction[ei][ej].topFlux[r][gp];
                        gridFluxForReconstruction[ei][ej].leftFlux[r][gp] = worldGridFluxForReconstruction[ei][ej].leftFlux[r][gp];
                        gridFluxForReconstruction[ei][ej].rightFlux[r][gp] = worldGridFluxForReconstruction[ei][ej].rightFlux[r][gp];
                    }
                }
            }
        }
    }
    else
    {
        // Other ranks receive data
        MPI_Status status;
        std::vector<double> recvBuffer((m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * 4 * m_gpNum); // 4 fluxes per grid element
        MPI_Recv(recvBuffer.data(), (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * 4 * m_gpNum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

        // Store received data in gridFluxForReconstruction
        for (int ei = 0; ei < m_elemNumX + 2 * m_ghostCellNum; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY + 2 * m_ghostCellNum; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    for (int gp = 0; gp < m_gpNum; ++gp)
                    {
                        const int index = ((ei * (m_elemNumY + 2 * m_ghostCellNum) + ej) * m_varNum + r) * m_gpNum + gp;
                        gridFluxForReconstruction[ei][ej].bottomFlux[r][gp] = recvBuffer[index];
                        gridFluxForReconstruction[ei][ej].topFlux[r][gp] = recvBuffer[index + (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum];
                        gridFluxForReconstruction[ei][ej].leftFlux[r][gp] = recvBuffer[index + 2 * (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum];
                        gridFluxForReconstruction[ei][ej].rightFlux[r][gp] = recvBuffer[index + 3 * (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum];
                    }
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
void CWENOFV::getFlux_WENOEXPandTHINCzhisetBoundaryGPs(Array2D<CgridFlux> &worldGridFluxForReconstruction)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // set the value of ghost cells by using boundary condition
        switch (equation->boundaryCondition)
        {
        case SMOOTH: // Periodic boundary
            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[r][gp] = worldGridFluxForReconstruction[m_worldEndElemX - 1][ej].rightFlux[r][gp]; // left
                        worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[r][gp] = worldGridFluxForReconstruction[m_worldStartElemX][ej].leftFlux[r][gp];           // right
                    }
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldEndElemY - 1].topFlux[r][gp]; // bottom
                        worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldStartElemY].bottomFlux[r][gp];   // top
                    }
                }
            }
            break;

        case NEUMANN: // Outflow boundary
            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[r][gp] = worldGridFluxForReconstruction[m_worldStartElemX][ej].leftFlux[r][gp]; // left
                        worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[r][gp] = worldGridFluxForReconstruction[m_worldEndElemX - 1][ej].rightFlux[r][gp];     // right
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[2][gp] *= -1;
                    worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[2][gp] *= -1;
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldStartElemY].bottomFlux[r][gp]; // bottom
                        worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldEndElemY - 1].topFlux[r][gp];     // top
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[3][gp] *= -1; // bottom
                    worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[3][gp] *= -1;    // top
                }
            }
            break;

        case FREESLIP: // Outflow boundary
            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[r][gp] = worldGridFluxForReconstruction[m_worldStartElemX][ej].leftFlux[r][gp]; // left
                        worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[r][gp] = worldGridFluxForReconstruction[m_worldEndElemX - 1][ej].rightFlux[r][gp];     // right
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[2][gp] = 0;
                    worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[2][gp] = 0;
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldStartElemY].bottomFlux[r][gp]; // bottom
                        worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldEndElemY - 1].topFlux[r][gp];     // top
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[3][gp] = 0; // bottom
                    worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[3][gp] = 0;    // top
                }
            }

            // left and right
            // for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            // {
            //     for (int r = 0; r != m_varNum; ++r)
            //     {
            //         for (int gp = 0; gp != m_gpNum; ++gp)
            //         {
            //             worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[r][gp] = worldGridFluxForReconstruction[m_worldStartElemX][ej].leftFlux[r][gp]; // left
            //             worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[r][gp] = worldGridFluxForReconstruction[m_worldEndElemX - 1][ej].rightFlux[r][gp];     // right
            //         }
            //     }

            //     for (int gp = 0; gp != m_gpNum; ++gp)
            //     {
            //         worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[2][gp] *= -1;
            //         worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[2][gp] *= -1;
            //     }
            // }

            // // bottom and top
            // for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            // {
            //     for (int r = 0; r != m_varNum; ++r)
            //     {
            //         for (int gp = 0; gp != m_gpNum; ++gp)
            //         {
            //             worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldStartElemY].bottomFlux[r][gp]; // bottom
            //             worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldEndElemY - 1].topFlux[r][gp];     // top
            //         }
            //     }

            //     for (int gp = 0; gp != m_gpNum; ++gp)
            //     {
            //         worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[3][gp] *= -1; // bottom
            //         worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[3][gp] *= -1;    // top
            //     }
            // }
            break;

        case NEUMANN2: // Outflow boundary
        {
            const int loc_gpNum = 4; // Number of Gauss points for local integration
            double gPointX, gWeightX;
            double gPointY, gWeightY;
            Array1D<double> loc_gpoints_ref(loc_gpNum), loc_gweights_ref(loc_gpNum);
            sc_math::GaussLegendre_ref(loc_gpNum, loc_gpoints_ref, loc_gweights_ref);

            Array1D<double> Conserved_var(m_varNum);

            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    // (1) worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux
                    gPointX = m_worldGrids[m_worldStartElemX - 1][ej].m_xRight;
                    gPointY = m_worldGrids[m_worldStartElemX - 1][ej].m_yCenter + 0.5 * m_worldGrids[m_worldStartElemX - 1][ej].m_yDistance * loc_gpoints_ref[gp];

                    // Get the initial conserved variables at the Gauss point
                    equation->getU0(gPointX, gPointY, Conserved_var);

                    for (int r = 0; r != m_varNum; ++r)
                        worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[r][gp] = Conserved_var[r]; // left

                    // (2) worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux
                    gPointX = m_worldGrids[m_worldEndElemX][ej].m_xLeft;
                    gPointY = m_worldGrids[m_worldEndElemX][ej].m_yCenter + 0.5 * m_worldGrids[m_worldStartElemX - 1][ej].m_yDistance * loc_gpoints_ref[gp];

                    // Get the initial conserved variables at the Gauss point
                    equation->getU0(gPointX, gPointY, Conserved_var);

                    for (int r = 0; r != m_varNum; ++r)
                        worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[r][gp] = Conserved_var[r]; // right
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    // (3) worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux
                    gPointX = m_worldGrids[ei][m_worldStartElemY - 1].m_xCenter + 0.5 * m_worldGrids[ei][m_worldStartElemY - 1].m_xDistance * loc_gpoints_ref[gp];
                    gPointY = m_worldGrids[ei][m_worldStartElemY - 1].m_yRight;

                    // Get the initial conserved variables at the Gauss point
                    equation->getU0(gPointX, gPointY, Conserved_var);

                    for (int r = 0; r != m_varNum; ++r)
                        worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[r][gp] = Conserved_var[r]; // bottom

                    // (4) worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux
                    gPointX = m_worldGrids[ei][m_worldEndElemY].m_xCenter + 0.5 * m_worldGrids[ei][m_worldEndElemY].m_xDistance * loc_gpoints_ref[gp];
                    gPointY = m_worldGrids[ei][m_worldEndElemY].m_yRight;

                    // Get the initial conserved variables at the Gauss point
                    equation->getU0(gPointX, gPointY, Conserved_var);

                    for (int r = 0; r != m_varNum; ++r)
                        worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[r][gp] = Conserved_var[r]; // top
                }
            }
            break;
        }

        case NoSlip: // Outflow boundary
            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[r][gp] = worldGridFluxForReconstruction[m_worldStartElemX][ej].leftFlux[r][gp]; // left
                        worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[r][gp] = worldGridFluxForReconstruction[m_worldEndElemX - 1][ej].rightFlux[r][gp];     // right
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[2][gp] = 0;
                    worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[2][gp] = 0;
                    worldGridFluxForReconstruction[m_worldStartElemX - 1][ej].rightFlux[3][gp] = 0;
                    worldGridFluxForReconstruction[m_worldEndElemX][ej].leftFlux[3][gp] = 0;
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldStartElemY].bottomFlux[r][gp]; // bottom
                        worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[r][gp] = worldGridFluxForReconstruction[ei][m_worldEndElemY - 1].topFlux[r][gp];     // top
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[2][gp] = 0; // bottom
                    worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[2][gp] = 0;    // top
                    worldGridFluxForReconstruction[ei][m_worldStartElemY - 1].topFlux[3][gp] = 0; // bottom
                    worldGridFluxForReconstruction[ei][m_worldEndElemY].bottomFlux[3][gp] = 0;    // top
                }
            }
            break;
        default:
            std::cout << "Error: Test case is not supported when set boundary..." << std::endl;
            std::cout << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
            std::cin.get();
            exit(1);
            break;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void CWENOFV::assembleBoundaryFaceTerm(void)
{
    double gWeightX, gWeightY;

    Array1D<double> flux_right(m_varNum);
    Array1D<double> flux_top(m_varNum);

    Array1D<double> u_right_minus(m_varNum);
    Array1D<double> u_right_plus(m_varNum);
    Array1D<double> u_top_minus(m_varNum);
    Array1D<double> u_top_plus(m_varNum);

    m_mainTimer.pause();
    // getWorldGridFlux();
    // worldToProcGridFlux();
    getFlux_WENOEXPandTHINCzhiMPI(m_gridFlux);

    m_mainTimer.start();

    // 创建通量缓存
    std::vector<std::vector<std::vector<Array1D<double>>>> flux_right_all(
        m_endElemX - m_startElemX + 1,
        std::vector<std::vector<Array1D<double>>>(m_endElemY - m_startElemY,
                                                  std::vector<Array1D<double>>(m_gpNum, Array1D<double>(m_varNum))));

    std::vector<std::vector<std::vector<Array1D<double>>>> flux_top_all(
        m_endElemX - m_startElemX,
        std::vector<std::vector<Array1D<double>>>(m_endElemY - m_startElemY + 1,
                                                  std::vector<Array1D<double>>(m_gpNum, Array1D<double>(m_varNum))));

    m_riemannTimer.start();
    // === 第一步：计算所有边上的黎曼通量 ===
    // x-direction
    for (int ei = m_startElemX - 1; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            for (int gp = 0; gp != m_gpNum; gp++)
            {
                for (int r = 0; r != m_varNum; r++)
                {
                    u_right_minus[r] = m_gridFlux[ei][ej].rightFlux[r][gp];
                    u_right_plus[r] = m_gridFlux[ei + 1][ej].leftFlux[r][gp];
                }

                equation->getRiemannFlux(u_right_minus, u_right_plus, flux_right, 1, 0, m_riemannFluxType);

                // 存储通量
                flux_right_all[ei - (m_startElemX - 1)][ej - m_startElemY][gp] = flux_right;
            }
        }
    }

    // y-direction
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY - 1; ej != m_endElemY; ej++)
        {
            for (int gp = 0; gp != m_gpNum; gp++)
            {
                for (int r = 0; r != m_varNum; r++)
                {
                    u_top_minus[r] = m_gridFlux[ei][ej].topFlux[r][gp];
                    u_top_plus[r] = m_gridFlux[ei][ej + 1].bottomFlux[r][gp];
                }

                equation->getRiemannFlux(u_top_minus, u_top_plus, flux_top, 0, 1, m_riemannFluxType);

                // 存储通量
                flux_top_all[ei - m_startElemX][ej - (m_startElemY - 1)][gp] = flux_top;
            }
        }
    }
    m_riemannTimer.pause();

    m_FVTimer.start();

    // 在 flux_right_all/flux_top_all 之后增加
    // x 方向每条边上的平均通量
    std::vector<std::vector<Array1D<double>>> flux_right_avg(
        m_endElemX - m_startElemX + 1,
        std::vector<Array1D<double>>(m_endElemY - m_startElemY, Array1D<double>(m_varNum)));

    // y 方向每条边上的平均通量
    std::vector<std::vector<Array1D<double>>> flux_top_avg(
        m_endElemX - m_startElemX,
        std::vector<Array1D<double>>(m_endElemY - m_startElemY + 1, Array1D<double>(m_varNum)));

    // === 第二步：统一更新 RHS ===
    // x-direction 平均值存储
    for (int ei = m_startElemX - 1; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            // 计算每条边（ei,ej）上各物理量的平均通量
            for (int r = 0; r < m_varNum; r++)
            {
                double ave = 0.0;
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    gWeightY = m_gweights_ref[gp] * m_grids[ei][ej].m_yDistance / 2.0;
                    ave += flux_right_all[ei - (m_startElemX - 1)][ej - m_startElemY][gp][r] * gWeightY;
                }
                // 存储平均值
                flux_right_avg[ei - (m_startElemX - 1)][ej - m_startElemY][r] = ave;
            }
        }
    }

    // y-direction 平均值存储
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY - 1; ej != m_endElemY; ej++)
        {
            for (int r = 0; r < m_varNum; r++)
            {
                double ave = 0.0;
                for (int gp = 0; gp != m_gpNum; gp++)
                {
                    gWeightX = m_gweights_ref[gp] * m_grids[ei][ej].m_xDistance / 2.0;
                    ave += flux_top_all[ei - m_startElemX][ej - (m_startElemY - 1)][gp][r] * gWeightX;
                }
                flux_top_avg[ei - m_startElemX][ej - (m_startElemY - 1)][r] = ave;
            }
        }
    }

    // x-direction 更新 RHS
    for (int ei = m_startElemX - 1; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            for (int r = 0; r < m_varNum; r++)
            {
                double ave = flux_right_avg[ei - (m_startElemX - 1)][ej - m_startElemY][r];
                m_rhs[ei][ej].vector[r] += -ave / m_grids[ei][ej].m_xDistance / m_grids[ei][ej].m_yDistance;
                m_rhs[ei + 1][ej].vector[r] += ave / m_grids[ei][ej].m_xDistance / m_grids[ei][ej].m_yDistance;
            }
        }
    }

    // y-direction 更新 RHS
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY - 1; ej != m_endElemY; ej++)
        {
            for (int r = 0; r < m_varNum; r++)
            {
                double ave = flux_top_avg[ei - m_startElemX][ej - (m_startElemY - 1)][r];
                m_rhs[ei][ej].vector[r] -= ave / m_grids[ei][ej].m_xDistance / m_grids[ei][ej].m_yDistance;
                m_rhs[ei][ej + 1].vector[r] += ave / m_grids[ei][ej].m_xDistance / m_grids[ei][ej + 1].m_yDistance;
            }
        }
    }
    m_FVTimer.pause();
}

void CWENOFV::assembleSourceTerm(void)
{
    Array1D<double> Uh(m_varNum);
    double g_g0;

    switch (m_testcase)
    {
    case DAMBREAK:
    case DAMBREAK2:
    case DAMBREAK3:
    case RTI:
        g_g0 = 9.81;
        for (int ei = m_startElemX; ei != m_endElemX; ei++)
        {
            for (int ej = m_startElemY; ej != m_endElemY; ej++)
            {
                Uh = m_Uh[ei][ej].vector;
                const double z1rho1 = Uh[0];
                const double z2rho2 = Uh[1];
                const double rho = z1rho1 + z2rho2;
                m_rhs[ei][ej].vector[3] -= rho * g_g0;
            }
        }
        break;
    case STANDING_WAVE:
        g_g0 = 2 * M_PI;
        for (int ei = m_startElemX; ei != m_endElemX; ei++)
        {
            for (int ej = m_startElemY; ej != m_endElemY; ej++)
            {
                Uh = m_Uh[ei][ej].vector;
                const double z1rho1 = Uh[0];
                const double z2rho2 = Uh[1];
                const double rho = z1rho1 + z2rho2;
                m_rhs[ei][ej].vector[3] -= rho * g_g0;
            }
        }
        break;
    default:
        break;
    }
}
double CWENOFV::useWENO(double uavemmm, double uavemm, double uavem, double uave, double uavep, double uavepp, double uaveppp, int gp)
{
    double u_hat(0);
    switch (m_scheme)
    {
    case WENO:
        u_hat = WENO5threconstruction(uavemm, uavem, uave, uavep, uavepp, gp);
        break;
    case WENOZ:
        u_hat = WENO5Zthreconstruction(uavemm, uavem, uave, uavep, uavepp, gp);
        break;
    case WENOEXP:
        u_hat = wenoExpthreconstruction(uavemmm, uavemm, uavem, uave, uavep, uavepp, uaveppp, gp, m_deltaX);
        break;
    default:
        std::cerr << "the scheme is not supported" << std::endl;
        break;
    }

    return u_hat;
}

void CWENOFV::setBoundaryAverages(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // set the value of ghost cells by using boundary condition
        switch (equation->boundaryCondition)
        {
        case PERIOD: // Periodic boundary
            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        m_worldUh[ei][m_worldStartElemY - e1].vector[r] = m_worldUh[ei][m_worldEndElemY - e1].vector[r];
                        m_worldUh[ei][m_worldEndElemY + e].vector[r] = m_worldUh[ei][m_worldStartElemY + e].vector[r];
                    }
                }
            }

            // left and right
            for (int ej = 0; ej != m_worldElemNumY + 2 * m_ghostCellNum; ej++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        m_worldUh[m_worldStartElemX - e1][ej].vector[r] = m_worldUh[m_worldEndElemX - e1][ej].vector[r];
                        m_worldUh[m_worldEndElemX + e][ej].vector[r] = m_worldUh[m_worldStartElemX + e][ej].vector[r];
                    }
                }
            }

            break;

        case NEUMANN: // Outflow boundary
            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        m_worldUh[ei][m_worldStartElemY - e1].vector[r] = m_worldUh[ei][m_worldStartElemY].vector[r]; // bottom
                        m_worldUh[ei][m_worldEndElemY + e].vector[r] = m_worldUh[ei][m_worldEndElemY - 1].vector[r];  // top
                    }

                    m_worldUh[ei][m_worldStartElemY - e1].vector[3] *= -1; // bottom
                    m_worldUh[ei][m_worldEndElemY + e].vector[3] *= -1;    // top
                }
            }

            // left and right
            for (int ej = 0; ej != m_worldElemNumY + 2 * m_ghostCellNum; ej++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        m_worldUh[m_worldStartElemX - e1][ej].vector[r] = m_worldUh[m_worldStartElemX][ej].vector[r];
                        m_worldUh[m_worldEndElemX + e][ej].vector[r] = m_worldUh[m_worldEndElemX - 1][ej].vector[r];
                    }

                    m_worldUh[m_worldStartElemX - e1][ej].vector[2] *= -1; // left
                    m_worldUh[m_worldEndElemX + e][ej].vector[2] *= -1;    // right
                }
            }
            break;

        case FREESLIP: // Outflow boundary
            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        m_worldUh[ei][m_worldStartElemY - e1].vector[r] = m_worldUh[ei][m_worldStartElemY].vector[r]; // bottom
                        m_worldUh[ei][m_worldEndElemY + e].vector[r] = m_worldUh[ei][m_worldEndElemY - 1].vector[r];  // top
                    }

                    m_worldUh[ei][m_worldStartElemY - e1].vector[3] = 0; // bottom
                    m_worldUh[ei][m_worldEndElemY + e].vector[3] = 0;    // top
                }
            }

            // left and right
            for (int ej = 0; ej != m_worldElemNumY + 2 * m_ghostCellNum; ej++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        m_worldUh[m_worldStartElemX - e1][ej].vector[r] = m_worldUh[m_worldStartElemX][ej].vector[r];
                        m_worldUh[m_worldEndElemX + e][ej].vector[r] = m_worldUh[m_worldEndElemX - 1][ej].vector[r];
                    }

                    m_worldUh[m_worldStartElemX - e1][ej].vector[2] = 0; // left
                    m_worldUh[m_worldEndElemX + e][ej].vector[2] = 0;    // right
                }
            }

            // for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            // {
            //     for (int e = 0; e != m_ghostCellNum; ++e)
            //     {
            //         int e1 = e + 1;
            //         for (int r = 0; r != m_varNum; ++r)
            //         {
            //             m_worldUh[ei][m_worldStartElemY - e1].vector[r] = m_worldUh[ei][m_worldStartElemY].vector[r]; // bottom
            //             m_worldUh[ei][m_worldEndElemY + e].vector[r] = m_worldUh[ei][m_worldEndElemY - 1].vector[r];  // top
            //         }

            //         m_worldUh[ei][m_worldStartElemY - e1].vector[3] *= -1; // bottom
            //         m_worldUh[ei][m_worldEndElemY + e].vector[3] *= -1;    // top
            //     }
            // }

            // // left and right
            // for (int ej = 0; ej != m_worldElemNumY + 2 * m_ghostCellNum; ej++)
            // {
            //     for (int e = 0; e != m_ghostCellNum; ++e)
            //     {
            //         int e1 = e + 1;
            //         for (int r = 0; r != m_varNum; ++r)
            //         {
            //             m_worldUh[m_worldStartElemX - e1][ej].vector[r] = m_worldUh[m_worldStartElemX][ej].vector[r];
            //             m_worldUh[m_worldEndElemX + e][ej].vector[r] = m_worldUh[m_worldEndElemX - 1][ej].vector[r];
            //         }

            //         m_worldUh[m_worldStartElemX - e1][ej].vector[2] *= -1; // left
            //         m_worldUh[m_worldEndElemX + e][ej].vector[2] *= -1;    // right
            //     }
            // }
            break;

        case NEUMANN2: // Outflow boundary
        {
            const int loc_gpNum = 4; // Number of Gauss points for local integration
            double gPointX, gWeightX;
            double gPointY, gWeightY;
            Array1D<double> loc_gpoints_ref(loc_gpNum), loc_gweights_ref(loc_gpNum);
            sc_math::GaussLegendre_ref(loc_gpNum, loc_gpoints_ref, loc_gweights_ref);

            Array1D<double> Conserved_var(m_varNum);

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;

                    // m_worldUh[ei][m_worldStartElemY - e1]
                    m_worldUh[ei][m_worldStartElemY - e1].vector.setZero();
                    for (int gpi = 0; gpi != loc_gpNum; ++gpi) // Loop over Gauss points for numerical integration
                    {
                        for (int gpj = 0; gpj != loc_gpNum; ++gpj)
                        {
                            // Calculate the Gauss point location within the cell
                            gPointX = m_worldGrids[ei][m_worldStartElemY - e1].m_xCenter + 0.5 * m_worldGrids[ei][m_worldStartElemY - e1].m_xDistance * loc_gpoints_ref[gpi];
                            gPointY = m_worldGrids[ei][m_worldStartElemY - e1].m_yCenter + 0.5 * m_worldGrids[ei][m_worldStartElemY - e1].m_yDistance * loc_gpoints_ref[gpj];

                            // Calculate the Gauss weights for the cell
                            gWeightX = 0.5 * m_worldGrids[ei][m_worldStartElemY - e1].m_xDistance * loc_gweights_ref[gpi];
                            gWeightY = 0.5 * m_worldGrids[ei][m_worldStartElemY - e1].m_yDistance * loc_gweights_ref[gpj];

                            // Get the initial conserved variables at the Gauss point
                            equation->getU0(gPointX, gPointY, Conserved_var);

                            // Accumulate the weighted conserved variables
                            for (int r = 0; r != m_varNum; ++r)
                                m_worldUh[ei][m_worldStartElemY - e1].vector[r] += Conserved_var[r] * gWeightX * gWeightY / m_worldGrids[ei][m_worldStartElemY - e1].m_area;
                        }
                    }

                    // m_worldUh[ei][m_worldEndElemY + e]
                    m_worldUh[ei][m_worldEndElemY + e].vector.setZero();
                    for (int gpi = 0; gpi != loc_gpNum; ++gpi) // Loop over Gauss points for numerical integration
                    {
                        for (int gpj = 0; gpj != loc_gpNum; ++gpj)
                        {
                            // Calculate the Gauss point location within the cell
                            gPointX = m_worldGrids[ei][m_worldEndElemY + e].m_xCenter + 0.5 * m_worldGrids[ei][m_worldEndElemY + e].m_xDistance * loc_gpoints_ref[gpi];
                            gPointY = m_worldGrids[ei][m_worldEndElemY + e].m_yCenter + 0.5 * m_worldGrids[ei][m_worldEndElemY + e].m_yDistance * loc_gpoints_ref[gpj];

                            // Calculate the Gauss weights for the cell
                            gWeightX = 0.5 * m_worldGrids[ei][m_worldEndElemY + e].m_xDistance * loc_gweights_ref[gpi];
                            gWeightY = 0.5 * m_worldGrids[ei][m_worldEndElemY + e].m_yDistance * loc_gweights_ref[gpj];

                            // Get the initial conserved variables at the Gauss point
                            equation->getU0(gPointX, gPointY, Conserved_var);

                            // Accumulate the weighted conserved variables
                            for (int r = 0; r != m_varNum; ++r)
                                m_worldUh[ei][m_worldEndElemY + e].vector[r] += Conserved_var[r] * gWeightX * gWeightY / m_worldGrids[ei][m_worldEndElemY + e].m_area;
                        }
                    }
                }
            }

            // left and right
            for (int ej = 0; ej != m_worldElemNumY + 2 * m_ghostCellNum; ej++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;

                    // m_worldUh[m_worldStartElemX - e1][ej]
                    m_worldUh[m_worldStartElemX - e1][ej].vector.setZero();
                    for (int gpi = 0; gpi != loc_gpNum; ++gpi) // Loop over Gauss points for numerical integration
                    {
                        for (int gpj = 0; gpj != loc_gpNum; ++gpj)
                        {
                            // Calculate the Gauss point location within the cell
                            gPointX = m_worldGrids[m_worldStartElemX - e1][ej].m_xCenter + 0.5 * m_worldGrids[m_worldStartElemX - e1][ej].m_xDistance * loc_gpoints_ref[gpi];
                            gPointY = m_worldGrids[m_worldStartElemX - e1][ej].m_yCenter + 0.5 * m_worldGrids[m_worldStartElemX - e1][ej].m_yDistance * loc_gpoints_ref[gpj];

                            // Calculate the Gauss weights for the cell
                            gWeightX = 0.5 * m_worldGrids[m_worldStartElemX - e1][ej].m_xDistance * loc_gweights_ref[gpi];
                            gWeightY = 0.5 * m_worldGrids[m_worldStartElemX - e1][ej].m_yDistance * loc_gweights_ref[gpj];

                            // Get the initial conserved variables at the Gauss point
                            equation->getU0(gPointX, gPointY, Conserved_var);

                            // Accumulate the weighted conserved variables
                            for (int r = 0; r != m_varNum; ++r)
                                m_worldUh[m_worldStartElemX - e1][ej].vector[r] += Conserved_var[r] * gWeightX * gWeightY / m_worldGrids[m_worldStartElemX - e1][ej].m_area;
                        }
                    }

                    // m_worldUh[m_worldEndElemX + e][ej]
                    m_worldUh[m_worldEndElemX + e][ej].vector.setZero();
                    for (int gpi = 0; gpi != loc_gpNum; ++gpi) // Loop over Gauss points for numerical integration
                    {
                        for (int gpj = 0; gpj != loc_gpNum; ++gpj)
                        {
                            // Calculate the Gauss point location within the cell
                            gPointX = m_worldGrids[m_worldEndElemX + e][ej].m_xCenter + 0.5 * m_worldGrids[m_worldEndElemX + e][ej].m_xDistance * loc_gpoints_ref[gpi];
                            gPointY = m_worldGrids[m_worldEndElemX + e][ej].m_yCenter + 0.5 * m_worldGrids[m_worldEndElemX + e][ej].m_yDistance * loc_gpoints_ref[gpj];

                            // Calculate the Gauss weights for the cell
                            gWeightX = 0.5 * m_worldGrids[m_worldEndElemX + e][ej].m_xDistance * loc_gweights_ref[gpi];
                            gWeightY = 0.5 * m_worldGrids[m_worldEndElemX + e][ej].m_yDistance * loc_gweights_ref[gpj];

                            // Get the initial conserved variables at the Gauss point
                            equation->getU0(gPointX, gPointY, Conserved_var);

                            // Accumulate the weighted conserved variables
                            for (int r = 0; r != m_varNum; ++r)
                                m_worldUh[m_worldEndElemX + e][ej].vector[r] += Conserved_var[r] * gWeightX * gWeightY / m_worldGrids[m_worldEndElemX + e][ej].m_area;
                        }
                    }
                }
            }

            break;
        }

        case NoSlip:
            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        m_worldUh[ei][m_worldStartElemY - e1].vector[r] = m_worldUh[ei][m_worldStartElemY].vector[r]; // bottom
                        m_worldUh[ei][m_worldEndElemY + e].vector[r] = m_worldUh[ei][m_worldEndElemY - 1].vector[r];  // top
                    }

                    m_worldUh[ei][m_worldStartElemY - e1].vector[2] = 0; // bottom
                    m_worldUh[ei][m_worldEndElemY + e].vector[2] = 0;    // top
                    m_worldUh[ei][m_worldStartElemY - e1].vector[3] = 0; // bottom
                    m_worldUh[ei][m_worldEndElemY + e].vector[3] = 0;    // top
                }
            }

            // left and right
            for (int ej = 0; ej != m_worldElemNumY + 2 * m_ghostCellNum; ej++)
            {
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        m_worldUh[m_worldStartElemX - e1][ej].vector[r] = m_worldUh[m_worldStartElemX][ej].vector[r];
                        m_worldUh[m_worldEndElemX + e][ej].vector[r] = m_worldUh[m_worldEndElemX - 1][ej].vector[r];
                    }

                    m_worldUh[m_worldStartElemX - e1][ej].vector[2] = 0; // left
                    m_worldUh[m_worldEndElemX + e][ej].vector[2] = 0;    // right
                    m_worldUh[m_worldStartElemX - e1][ej].vector[3] = 0; // left
                    m_worldUh[m_worldEndElemX + e][ej].vector[3] = 0;    // right
                }
            }
            break;

        default:
            std::cout << "Error: Test case is not supported when set boundary..." << std::endl;
            std::cout << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
            std::cin.get();
            exit(1);
            break;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
void CWENOFV::setBoundaryGPs(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // set the value of ghost cells by using boundary condition
        switch (equation->boundaryCondition)
        {
        case SMOOTH: // Periodic boundary
            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[r][gp] = m_worldGridFlux[m_worldEndElemX - 1][ej].rightFlux[r][gp]; // left
                        m_worldGridFlux[m_worldEndElemX][ej].leftFlux[r][gp] = m_worldGridFlux[m_worldStartElemX][ej].leftFlux[r][gp];           // right
                    }
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[r][gp] = m_worldGridFlux[ei][m_worldEndElemY - 1].topFlux[r][gp]; // bottom
                        m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[r][gp] = m_worldGridFlux[ei][m_worldStartElemY].bottomFlux[r][gp];   // top
                    }
                }
            }
            break;

        case NEUMANN: // Outflow boundary
            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[r][gp] = m_worldGridFlux[m_worldStartElemX][ej].leftFlux[r][gp]; // left
                        m_worldGridFlux[m_worldEndElemX][ej].leftFlux[r][gp] = m_worldGridFlux[m_worldEndElemX - 1][ej].rightFlux[r][gp];     // right
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[2][gp] *= -1;
                    m_worldGridFlux[m_worldEndElemX][ej].leftFlux[2][gp] *= -1;
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[r][gp] = m_worldGridFlux[ei][m_worldStartElemY].bottomFlux[r][gp]; // bottom
                        m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[r][gp] = m_worldGridFlux[ei][m_worldEndElemY - 1].topFlux[r][gp];     // top
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[3][gp] *= -1; // bottom
                    m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[3][gp] *= -1;    // top
                }
            }
            break;

        case FREESLIP: // Outflow boundary
            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[r][gp] = m_worldGridFlux[m_worldStartElemX][ej].leftFlux[r][gp]; // left
                        m_worldGridFlux[m_worldEndElemX][ej].leftFlux[r][gp] = m_worldGridFlux[m_worldEndElemX - 1][ej].rightFlux[r][gp];     // right
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[2][gp] = 0;
                    m_worldGridFlux[m_worldEndElemX][ej].leftFlux[2][gp] = 0;
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[r][gp] = m_worldGridFlux[ei][m_worldStartElemY].bottomFlux[r][gp]; // bottom
                        m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[r][gp] = m_worldGridFlux[ei][m_worldEndElemY - 1].topFlux[r][gp];     // top
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[3][gp] = 0; // bottom
                    m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[3][gp] = 0;    // top
                }
            }

            // left and right
            // for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            // {
            //     for (int r = 0; r != m_varNum; ++r)
            //     {
            //         for (int gp = 0; gp != m_gpNum; ++gp)
            //         {
            //             m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[r][gp] = m_worldGridFlux[m_worldStartElemX][ej].leftFlux[r][gp]; // left
            //             m_worldGridFlux[m_worldEndElemX][ej].leftFlux[r][gp] = m_worldGridFlux[m_worldEndElemX - 1][ej].rightFlux[r][gp];     // right
            //         }
            //     }

            //     for (int gp = 0; gp != m_gpNum; ++gp)
            //     {
            //         m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[2][gp] *= -1;
            //         m_worldGridFlux[m_worldEndElemX][ej].leftFlux[2][gp] *= -1;
            //     }
            // }

            // // bottom and top
            // for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            // {
            //     for (int r = 0; r != m_varNum; ++r)
            //     {
            //         for (int gp = 0; gp != m_gpNum; ++gp)
            //         {
            //             m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[r][gp] = m_worldGridFlux[ei][m_worldStartElemY].bottomFlux[r][gp]; // bottom
            //             m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[r][gp] = m_worldGridFlux[ei][m_worldEndElemY - 1].topFlux[r][gp];     // top
            //         }
            //     }

            //     for (int gp = 0; gp != m_gpNum; ++gp)
            //     {
            //         m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[3][gp] *= -1; // bottom
            //         m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[3][gp] *= -1;    // top
            //     }
            // }
            break;

        case NEUMANN2: // Outflow boundary
        {
            const int loc_gpNum = 4; // Number of Gauss points for local integration
            double gPointX, gWeightX;
            double gPointY, gWeightY;
            Array1D<double> loc_gpoints_ref(loc_gpNum), loc_gweights_ref(loc_gpNum);
            sc_math::GaussLegendre_ref(loc_gpNum, loc_gpoints_ref, loc_gweights_ref);

            Array1D<double> Conserved_var(m_varNum);

            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    // (1) m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux
                    gPointX = m_worldGrids[m_worldStartElemX - 1][ej].m_xRight;
                    gPointY = m_worldGrids[m_worldStartElemX - 1][ej].m_yCenter + 0.5 * m_worldGrids[m_worldStartElemX - 1][ej].m_yDistance * loc_gpoints_ref[gp];

                    // Get the initial conserved variables at the Gauss point
                    equation->getU0(gPointX, gPointY, Conserved_var);

                    for (int r = 0; r != m_varNum; ++r)
                        m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[r][gp] = Conserved_var[r]; // left

                    // (2) m_worldGridFlux[m_worldEndElemX][ej].leftFlux
                    gPointX = m_worldGrids[m_worldEndElemX][ej].m_xLeft;
                    gPointY = m_worldGrids[m_worldEndElemX][ej].m_yCenter + 0.5 * m_worldGrids[m_worldStartElemX - 1][ej].m_yDistance * loc_gpoints_ref[gp];

                    // Get the initial conserved variables at the Gauss point
                    equation->getU0(gPointX, gPointY, Conserved_var);

                    for (int r = 0; r != m_varNum; ++r)
                        m_worldGridFlux[m_worldEndElemX][ej].leftFlux[r][gp] = Conserved_var[r]; // right
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    // (3) m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux
                    gPointX = m_worldGrids[ei][m_worldStartElemY - 1].m_xCenter + 0.5 * m_worldGrids[ei][m_worldStartElemY - 1].m_xDistance * loc_gpoints_ref[gp];
                    gPointY = m_worldGrids[ei][m_worldStartElemY - 1].m_yRight;

                    // Get the initial conserved variables at the Gauss point
                    equation->getU0(gPointX, gPointY, Conserved_var);

                    for (int r = 0; r != m_varNum; ++r)
                        m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[r][gp] = Conserved_var[r]; // bottom

                    // (4) m_worldGridFlux[ei][m_worldEndElemY].bottomFlux
                    gPointX = m_worldGrids[ei][m_worldEndElemY].m_xCenter + 0.5 * m_worldGrids[ei][m_worldEndElemY].m_xDistance * loc_gpoints_ref[gp];
                    gPointY = m_worldGrids[ei][m_worldEndElemY].m_yRight;

                    // Get the initial conserved variables at the Gauss point
                    equation->getU0(gPointX, gPointY, Conserved_var);

                    for (int r = 0; r != m_varNum; ++r)
                        m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[r][gp] = Conserved_var[r]; // top
                }
            }
            break;
        }

        case NoSlip: // Outflow boundary
            // left and right
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[r][gp] = m_worldGridFlux[m_worldStartElemX][ej].leftFlux[r][gp]; // left
                        m_worldGridFlux[m_worldEndElemX][ej].leftFlux[r][gp] = m_worldGridFlux[m_worldEndElemX - 1][ej].rightFlux[r][gp];     // right
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[2][gp] = 0;
                    m_worldGridFlux[m_worldEndElemX][ej].leftFlux[2][gp] = 0;
                    m_worldGridFlux[m_worldStartElemX - 1][ej].rightFlux[3][gp] = 0;
                    m_worldGridFlux[m_worldEndElemX][ej].leftFlux[3][gp] = 0;
                }
            }

            // bottom and top
            for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
            {
                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[r][gp] = m_worldGridFlux[ei][m_worldStartElemY].bottomFlux[r][gp]; // bottom
                        m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[r][gp] = m_worldGridFlux[ei][m_worldEndElemY - 1].topFlux[r][gp];     // top
                    }
                }

                for (int gp = 0; gp != m_gpNum; ++gp)
                {
                    m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[2][gp] = 0; // bottom
                    m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[2][gp] = 0;    // top
                    m_worldGridFlux[ei][m_worldStartElemY - 1].topFlux[3][gp] = 0; // bottom
                    m_worldGridFlux[ei][m_worldEndElemY].bottomFlux[3][gp] = 0;    // top
                }
            }
            break;
        default:
            std::cout << "Error: Test case is not supported when set boundary..." << std::endl;
            std::cout << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
            std::cin.get();
            exit(1);
            break;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void CWENOFV::useBoundPreservingLimiter(void)
{
    double theta(1.0);

    // loop over each cell
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            // positivity-preserving limiting for partial density
            theta = getPositivityPreservingTheta(ei, ej);
            if (theta < 1 - 1e-8)
                m_totalTheta[ei][ej] = 1;
            else
                m_totalTheta[ei][ej] = 0;

            if (theta >= 0.0 && theta <= 1.0)
            {
                // if (theta < 1.0)
                //     std::cout << "111" << std::endl;

                for (int r = 0; r != m_varNum; ++r)
                {
                    for (int gp = 0; gp != m_gpNum; ++gp)
                    {
                        m_gridFlux[ei][ej].rightFlux[r][gp] = theta * (m_gridFlux[ei][ej].rightFlux[r][gp] - m_Uh[ei][ej].vector[r]) + m_Uh[ei][ej].vector[r];
                        m_gridFlux[ei][ej].leftFlux[r][gp] = theta * (m_gridFlux[ei][ej].leftFlux[r][gp] - m_Uh[ei][ej].vector[r]) + m_Uh[ei][ej].vector[r];
                        m_gridFlux[ei][ej].topFlux[r][gp] = theta * (m_gridFlux[ei][ej].topFlux[r][gp] - m_Uh[ei][ej].vector[r]) + m_Uh[ei][ej].vector[r];
                        m_gridFlux[ei][ej].bottomFlux[r][gp] = theta * (m_gridFlux[ei][ej].bottomFlux[r][gp] - m_Uh[ei][ej].vector[r]) + m_Uh[ei][ej].vector[r];
                    }
                }
            }
            else
            {
                std::cout << "Error: illegal theta in PP limiting..." << std::endl;
                std::cout << "Theta=" << theta << std::endl;
                std::cin.get();
            }
        }
    }
}
double CWENOFV::getPositivityPreservingTheta(int ei, int ej)
{
    //////////////////////////////////////////////////////////////////
    /////// warning: only for two-phase three equation model. ///////
    //////////////////////////////////////////////////////////////////

    // Positivity-preserving limiter for partial density
    double zkrhok_mloc(0);
    const double zeps = 1e-12;
    const double omega1 = 0.1;
    // const double omega1 = 0.2;
    const double omegaR = 1.0 - 2 * omega1;

    Array1D<double> loc_Uh(m_varNum);
    for (int r = 0; r != m_varNum; r++)
        loc_Uh[r] = m_Uh[ei][ej].vector[r];

    Array1D<double> thetazkrhok(2);
    const double a1 = equation->getMaxEigenValue(loc_Uh, 1, 0);
    const double a2 = equation->getMaxEigenValue(loc_Uh, 0, 1);
    const double mu1 = (a1 / m_deltaX) / (a1 / m_deltaX + a2 / m_deltaY);
    const double mu2 = (a2 / m_deltaY) / (a1 / m_deltaX + a2 / m_deltaY);

    double zkrhokRomaining(0);
    double uave;
    double eps;

    // use positivity-preserving limiter for partial density
    for (int r = 0; r != 2; r++)
    {
        thetazkrhok[r] = 1.0;

        uave = loc_Uh[r];
        // find loc minimum of rho1z1

        eps = min(zeps, uave);

        // getPR
        zkrhokRomaining = uave;
        for (int gp = 0; gp != m_gpNum; gp++)
        {
            zkrhokRomaining -= omega1 * m_gweights_ref[gp] / 2 * mu2 * m_gridFlux[ei][ej].bottomFlux[r][gp];
            zkrhokRomaining -= omega1 * m_gweights_ref[gp] / 2 * mu2 * m_gridFlux[ei][ej].topFlux[r][gp];
            zkrhokRomaining -= omega1 * m_gweights_ref[gp] / 2 * mu1 * m_gridFlux[ei][ej].leftFlux[r][gp];
            zkrhokRomaining -= omega1 * m_gweights_ref[gp] / 2 * mu1 * m_gridFlux[ei][ej].rightFlux[r][gp];
        }
        zkrhokRomaining /= omegaR;

        // get loc min
        zkrhok_mloc = zkrhokRomaining;
        for (int gp = 0; gp != m_gpNum; gp++)
        {
            zkrhok_mloc = sc_math::Min(zkrhok_mloc, m_gridFlux[ei][ej].bottomFlux[r][gp]);
            zkrhok_mloc = sc_math::Min(zkrhok_mloc, m_gridFlux[ei][ej].topFlux[r][gp]);
            zkrhok_mloc = sc_math::Min(zkrhok_mloc, m_gridFlux[ei][ej].leftFlux[r][gp]);
            zkrhok_mloc = sc_math::Min(zkrhok_mloc, m_gridFlux[ei][ej].rightFlux[r][gp]);
        }

        // find scalling factor for zk_rhok
        thetazkrhok[r] = sc_math::Min(fabs((uave - eps) / (uave - zkrhok_mloc)), thetazkrhok[r]);
    }

    const double theta = sc_math::Min(thetazkrhok[0], thetazkrhok[1]);
    return theta;
}

void CWENOFV::useFVLimiter()
{
    bool isLimiter = false;
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
    {
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
        {
            for (int gp = 0; gp != m_gpNum; gp++)
            {
                // rightFlux
                isLimiter = false;
                for (int r = 0; r != 2; r++)
                    if (fabs(m_gridFlux[ei][ej].rightFlux[r][gp] - m_Uh[ei][ej].vector[r]) > 0.9 * m_Uh[ei][ej].vector[r])
                        isLimiter = true;
                // if (isLimiter)
                //     for (int r = 0; r != m_varNum; r++)
                //         m_gridFlux[ei][ej].rightFlux[r][gp] = m_Uh[ei][ej].vector[r];
                if (isLimiter)
                {
                    // debug output
                    for (int r = 0; r != m_varNum; r++)
                        std::cout << scientific << "rightFlux: " << m_gridFlux[ei][ej].rightFlux[r][gp] << " Uh: " << m_Uh[ei][ej].vector[r] << std::endl;

                    std::cin.get();
                }

                // leftFlux
                isLimiter = false;
                for (int r = 0; r != 2; r++)
                    if (fabs(m_gridFlux[ei][ej].leftFlux[r][gp] - m_Uh[ei][ej].vector[r]) > 0.9 * m_Uh[ei][ej].vector[r])
                        isLimiter = true;
                if (isLimiter)
                    for (int r = 0; r != m_varNum; r++)
                        m_gridFlux[ei][ej].leftFlux[r][gp] = m_Uh[ei][ej].vector[r];

                // topFlux
                isLimiter = false;
                for (int r = 0; r != 2; r++)
                    if (fabs(m_gridFlux[ei][ej].topFlux[r][gp] - m_Uh[ei][ej].vector[r]) > 0.9 * m_Uh[ei][ej].vector[r])
                        isLimiter = true;
                if (isLimiter)
                    for (int r = 0; r != m_varNum; r++)
                        m_gridFlux[ei][ej].topFlux[r][gp] = m_Uh[ei][ej].vector[r];

                // bottomFlux
                isLimiter = false;
                for (int r = 0; r != 2; r++)
                    if (fabs(m_gridFlux[ei][ej].bottomFlux[r][gp] - m_Uh[ei][ej].vector[r]) > 0.9 * m_Uh[ei][ej].vector[r])
                        isLimiter = true;
                if (isLimiter)
                    for (int r = 0; r != m_varNum; r++)
                        m_gridFlux[ei][ej].bottomFlux[r][gp] = m_Uh[ei][ej].vector[r];
                isLimiter = false;
            }
        }
    }
}
void CWENOFV::getWorldUh()
{
    m_MPITimer.start();
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // Receive data from other processes
        for (int src_rank = 1; src_rank < m_size; ++src_rank)
        {
            MPI_Status status;
            std::vector<double> recvBuffer(m_elemNumX * m_elemNumY * m_varNum);
            MPI_Recv(recvBuffer.data(), m_elemNumX * m_elemNumY * m_varNum, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &status);

            // Calculate the position of the data in worldUh
            const int offsetX = src_rank % m_numProcessorsX * m_elemNumX;
            const int offsetY = src_rank / m_numProcessorsX * m_elemNumY;

            // Store the received data in worldUh
            for (int ei = 0; ei < m_elemNumX; ++ei)
            {
                for (int ej = 0; ej < m_elemNumY; ++ej)
                {
                    for (int r = 0; r < m_varNum; ++r)
                    {
                        int index = (ei * m_elemNumY * m_varNum) + (ej * m_varNum) + r;
                        m_worldUh[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum].vector[r] = recvBuffer[index];
                    }
                }
            }
        }

        // Handle data for rank 0 itself
        for (int ei = 0; ei != m_elemNumX; ei++)
            for (int ej = 0; ej != m_elemNumY; ej++)
                for (int r = 0; r != m_varNum; r++)
                    m_worldUh[ei + m_ghostCellNum][ej + m_ghostCellNum].vector[r] = m_Uh[ei + m_ghostCellNum][ej + m_ghostCellNum].vector[r];
    }
    else
    {
        // Serialize m_Uh data into a 1D array
        std::vector<double> sendBuffer(m_elemNumX * m_elemNumY * m_varNum);
        for (int ei = 0; ei < m_elemNumX; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    int index = (ei * m_elemNumY * m_varNum) + (ej * m_varNum) + r;
                    sendBuffer[index] = m_Uh[ei + m_ghostCellNum][ej + m_ghostCellNum].vector[r];
                }
            }
        }

        // Send data to rank 0
        MPI_Send(sendBuffer.data(), m_elemNumX * m_elemNumY * m_varNum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    setBoundaryAverages();
    m_MPITimer.pause();
}
void CWENOFV::worldToProcUh()
{
    m_MPITimer.start();
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // Rank 0 sends data to other ranks
        for (int dest_rank = 1; dest_rank < m_size; ++dest_rank)
        {
            int offsetX = (dest_rank % m_numProcessorsX) * m_elemNumX;
            int offsetY = (dest_rank / m_numProcessorsX) * m_elemNumY;

            std::vector<double> sendBuffer((m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum);
            for (int ei = 0; ei < m_elemNumX + 2 * m_ghostCellNum; ++ei)
            {
                for (int ej = 0; ej < m_elemNumY + 2 * m_ghostCellNum; ++ej)
                {
                    for (int r = 0; r < m_varNum; ++r)
                    {
                        int index = (ei * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum) + (ej * m_varNum) + r;
                        sendBuffer[index] = m_worldUh[ei + offsetX][ej + offsetY].vector[r];
                    }
                }
            }

            // Send data to dest_rank
            MPI_Send(sendBuffer.data(), (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
        }

        // Handle data for rank 0 itself
        for (int ei = 0; ei != m_elemNumX + 2 * m_ghostCellNum; ei++)
            for (int ej = 0; ej != m_elemNumY + 2 * m_ghostCellNum; ej++)
                for (int r = 0; r != m_varNum; r++)
                    m_Uh[ei][ej].vector[r] = m_worldUh[ei][ej].vector[r];
    }
    else
    {
        // Other ranks receive data
        MPI_Status status;
        std::vector<double> recvBuffer((m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum);
        MPI_Recv(recvBuffer.data(), (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

        // Store the received data in m_Uh
        for (int ei = 0; ei < m_elemNumX + 2 * m_ghostCellNum; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY + 2 * m_ghostCellNum; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    int index = (ei * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum) + (ej * m_varNum) + r;
                    m_Uh[ei][ej].vector[r] = recvBuffer[index];
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    m_MPITimer.pause();
}

void CWENOFV::getWorldGridFlux()
{
    m_MPITimer.start();
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // Receive data from other ranks
        for (int src_rank = 1; src_rank < m_size; ++src_rank)
        {
            MPI_Status status;
            std::vector<double> recvBuffer(m_elemNumX * m_elemNumY * m_varNum * 4 * m_gpNum); // 4 fluxes per grid element
            MPI_Recv(recvBuffer.data(), m_elemNumX * m_elemNumY * m_varNum * 4 * m_gpNum, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &status);

            // Calculate the position in m_worldGridFlux
            int offsetX = (src_rank % m_numProcessorsX) * m_elemNumX;
            int offsetY = (src_rank / m_numProcessorsX) * m_elemNumY;

            // Store received data in m_worldGridFlux
            for (int ei = 0; ei < m_elemNumX; ++ei)
            {
                for (int ej = 0; ej < m_elemNumY; ++ej)
                {
                    for (int r = 0; r < m_varNum; ++r)
                    {
                        for (int gp = 0; gp < m_gpNum; ++gp)
                        {
                            int index = ((ei * m_elemNumY + ej) * m_varNum + r) * m_gpNum + gp;
                            m_worldGridFlux[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum].bottomFlux[r][gp] = recvBuffer[index];
                            m_worldGridFlux[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum].topFlux[r][gp] = recvBuffer[index + m_elemNumX * m_elemNumY * m_varNum * m_gpNum];
                            m_worldGridFlux[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum].leftFlux[r][gp] = recvBuffer[index + 2 * m_elemNumX * m_elemNumY * m_varNum * m_gpNum];
                            m_worldGridFlux[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum].rightFlux[r][gp] = recvBuffer[index + 3 * m_elemNumX * m_elemNumY * m_varNum * m_gpNum];
                        }
                    }
                }
            }
        }

        // Process data for rank 0 itself
        for (int ei = 0; ei < m_elemNumX; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    for (int gp = 0; gp < m_gpNum; ++gp)
                    {
                        m_worldGridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].bottomFlux[r][gp] = m_gridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].bottomFlux[r][gp];
                        m_worldGridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].topFlux[r][gp] = m_gridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].topFlux[r][gp];
                        m_worldGridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].leftFlux[r][gp] = m_gridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].leftFlux[r][gp];
                        m_worldGridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].rightFlux[r][gp] = m_gridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].rightFlux[r][gp];
                    }
                }
            }
        }
    }
    else
    {
        // Serialize m_gridFlux data to a one-dimensional array
        std::vector<double> sendBuffer(m_elemNumX * m_elemNumY * m_varNum * 4 * m_gpNum); // 4 fluxes per grid element
        for (int ei = 0; ei < m_elemNumX; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    for (int gp = 0; gp < m_gpNum; ++gp)
                    {
                        int index = ((ei * m_elemNumY + ej) * m_varNum + r) * m_gpNum + gp;
                        sendBuffer[index] = m_gridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].bottomFlux[r][gp];
                        sendBuffer[index + m_elemNumX * m_elemNumY * m_varNum * m_gpNum] = m_gridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].topFlux[r][gp];
                        sendBuffer[index + 2 * m_elemNumX * m_elemNumY * m_varNum * m_gpNum] = m_gridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].leftFlux[r][gp];
                        sendBuffer[index + 3 * m_elemNumX * m_elemNumY * m_varNum * m_gpNum] = m_gridFlux[ei + m_ghostCellNum][ej + m_ghostCellNum].rightFlux[r][gp];
                    }
                }
            }
        }

        // Send data to rank 0
        MPI_Send(sendBuffer.data(), m_elemNumX * m_elemNumY * m_varNum * 4 * m_gpNum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    // setBoundaryGPs();
    getFlux_WENOEXPandTHINCzhisetBoundaryGPs(m_worldGridFlux);
    m_MPITimer.pause();
}
void CWENOFV::worldToProcGridFlux()
{
    m_MPITimer.start();

    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // Rank 0 sends data to other ranks
        for (int dest_rank = 1; dest_rank < m_size; ++dest_rank)
        {
            const int offsetX = (dest_rank % m_numProcessorsX) * m_elemNumX;
            const int offsetY = (dest_rank / m_numProcessorsX) * m_elemNumY;

            std::vector<double> sendBuffer((m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * 4 * m_gpNum); // 4 fluxes per grid element
            for (int ei = 0; ei < m_elemNumX + 2 * m_ghostCellNum; ++ei)
            {
                for (int ej = 0; ej < m_elemNumY + 2 * m_ghostCellNum; ++ej)
                {
                    for (int r = 0; r < m_varNum; ++r)
                    {
                        for (int gp = 0; gp < m_gpNum; ++gp)
                        {
                            const int index = ((ei * (m_elemNumY + 2 * m_ghostCellNum) + ej) * m_varNum + r) * m_gpNum + gp;
                            sendBuffer[index] = m_worldGridFlux[ei + offsetX][ej + offsetY].bottomFlux[r][gp];
                            sendBuffer[index + (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum] = m_worldGridFlux[ei + offsetX][ej + offsetY].topFlux[r][gp];
                            sendBuffer[index + 2 * (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum] = m_worldGridFlux[ei + offsetX][ej + offsetY].leftFlux[r][gp];
                            sendBuffer[index + 3 * (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum] = m_worldGridFlux[ei + offsetX][ej + offsetY].rightFlux[r][gp];
                        }
                    }
                }
            }

            // Send data to dest_rank
            MPI_Send(sendBuffer.data(), (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * 4 * m_gpNum, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
        }

        // Process data for rank 0 itself
        for (int ei = 0; ei < m_elemNumX + 2 * m_ghostCellNum; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY + 2 * m_ghostCellNum; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    for (int gp = 0; gp < m_gpNum; ++gp)
                    {
                        m_gridFlux[ei][ej].bottomFlux[r][gp] = m_worldGridFlux[ei][ej].bottomFlux[r][gp];
                        m_gridFlux[ei][ej].topFlux[r][gp] = m_worldGridFlux[ei][ej].topFlux[r][gp];
                        m_gridFlux[ei][ej].leftFlux[r][gp] = m_worldGridFlux[ei][ej].leftFlux[r][gp];
                        m_gridFlux[ei][ej].rightFlux[r][gp] = m_worldGridFlux[ei][ej].rightFlux[r][gp];
                    }
                }
            }
        }
    }
    else
    {
        // Other ranks receive data
        MPI_Status status;
        std::vector<double> recvBuffer((m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * 4 * m_gpNum); // 4 fluxes per grid element
        MPI_Recv(recvBuffer.data(), (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * 4 * m_gpNum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

        // Store received data in m_gridFlux
        for (int ei = 0; ei < m_elemNumX + 2 * m_ghostCellNum; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY + 2 * m_ghostCellNum; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    for (int gp = 0; gp < m_gpNum; ++gp)
                    {
                        const int index = ((ei * (m_elemNumY + 2 * m_ghostCellNum) + ej) * m_varNum + r) * m_gpNum + gp;
                        m_gridFlux[ei][ej].bottomFlux[r][gp] = recvBuffer[index];
                        m_gridFlux[ei][ej].topFlux[r][gp] = recvBuffer[index + (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum];
                        m_gridFlux[ei][ej].leftFlux[r][gp] = recvBuffer[index + 2 * (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum];
                        m_gridFlux[ei][ej].rightFlux[r][gp] = recvBuffer[index + 3 * (m_elemNumX + 2 * m_ghostCellNum) * (m_elemNumY + 2 * m_ghostCellNum) * m_varNum * m_gpNum];
                    }
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    m_MPITimer.pause();
}

void CWENOFV::getWorldTheta(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // Receive data from other processes
        for (int src_rank = 1; src_rank < m_size; ++src_rank)
        {
            MPI_Status status;
            std::vector<double> recvBuffer(m_elemNumX * m_elemNumY);
            MPI_Recv(recvBuffer.data(), m_elemNumX * m_elemNumY, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &status);

            // Calculate the position of the data in worldUh
            const int offsetX = src_rank % m_numProcessorsX * m_elemNumX;
            const int offsetY = src_rank / m_numProcessorsX * m_elemNumY;

            // Store the received data in worldUh
            for (int ei = 0; ei < m_elemNumX; ++ei)
            {
                for (int ej = 0; ej < m_elemNumY; ++ej)
                {
                    int index = (ei * m_elemNumY) + ej;
                    m_worldTotalTheta[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum] = recvBuffer[index];
                }
            }
        }

        // Handle data for rank 0 itself
        for (int ei = 0; ei != m_elemNumX; ei++)
            for (int ej = 0; ej != m_elemNumY; ej++)
                m_worldTotalTheta[ei + m_ghostCellNum][ej + m_ghostCellNum] = m_totalTheta[ei + m_ghostCellNum][ej + m_ghostCellNum];
    }
    else
    {
        // Serialize m_Uh data into a 1D array
        std::vector<double> sendBuffer(m_elemNumX * m_elemNumY);
        for (int ei = 0; ei < m_elemNumX; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY; ++ej)
            {
                int index = (ei * m_elemNumY) + ej;
                sendBuffer[index] = m_totalTheta[ei + m_ghostCellNum][ej + m_ghostCellNum];
            }
        }

        // Send data to rank 0
        MPI_Send(sendBuffer.data(), m_elemNumX * m_elemNumY, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
void CWENOFV::getWorldPr(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // Receive data from other processes
        for (int src_rank = 1; src_rank < m_size; ++src_rank)
        {
            MPI_Status status;
            std::vector<double> recvBuffer(m_elemNumX * m_elemNumY);
            MPI_Recv(recvBuffer.data(), m_elemNumX * m_elemNumY, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &status);

            // Calculate the position of the data in worldUh
            const int offsetX = src_rank % m_numProcessorsX * m_elemNumX;
            const int offsetY = src_rank / m_numProcessorsX * m_elemNumY;

            // Store the received data in worldUh
            for (int ei = 0; ei < m_elemNumX; ++ei)
            {
                for (int ej = 0; ej < m_elemNumY; ++ej)
                {
                    int index = (ei * m_elemNumY) + ej;
                    m_worldTotalPr[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum] = recvBuffer[index];
                }
            }
        }

        // Handle data for rank 0 itself
        for (int ei = 0; ei != m_elemNumX; ei++)
            for (int ej = 0; ej != m_elemNumY; ej++)
                m_worldTotalPr[ei + m_ghostCellNum][ej + m_ghostCellNum] = m_totalPr[ei + m_ghostCellNum][ej + m_ghostCellNum];
    }
    else
    {
        // Serialize m_Uh data into a 1D array
        std::vector<double> sendBuffer(m_elemNumX * m_elemNumY);
        for (int ei = 0; ei < m_elemNumX; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY; ++ej)
            {
                int index = (ei * m_elemNumY) + ej;
                sendBuffer[index] = m_totalPr[ei + m_ghostCellNum][ej + m_ghostCellNum];
            }
        }

        // Send data to rank 0
        MPI_Send(sendBuffer.data(), m_elemNumX * m_elemNumY, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
void CWENOFV::getWorldDeltaT(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_rank == 0)
    {
        // Receive data from other processes
        for (int src_rank = 1; src_rank < m_size; ++src_rank)
        {
            MPI_Status status;
            std::vector<double> recvBuffer(m_elemNumX * m_elemNumY);
            MPI_Recv(recvBuffer.data(), m_elemNumX * m_elemNumY, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &status);

            // Calculate the position of the data in worldUh
            const int offsetX = src_rank % m_numProcessorsX * m_elemNumX;
            const int offsetY = src_rank / m_numProcessorsX * m_elemNumY;

            // Store the received data in worldUh
            for (int ei = 0; ei < m_elemNumX; ++ei)
            {
                for (int ej = 0; ej < m_elemNumY; ++ej)
                {
                    int index = (ei * m_elemNumY) + ej;
                    m_worldTotalDeltaT[ei + offsetX + m_ghostCellNum][ej + offsetY + m_ghostCellNum] = recvBuffer[index];
                }
            }
        }

        // Handle data for rank 0 itself
        for (int ei = 0; ei != m_elemNumX; ei++)
            for (int ej = 0; ej != m_elemNumY; ej++)
                m_worldTotalDeltaT[ei + m_ghostCellNum][ej + m_ghostCellNum] = m_totalDeltaT[ei + m_ghostCellNum][ej + m_ghostCellNum];
    }
    else
    {
        // Serialize data into a 1D array
        std::vector<double> sendBuffer(m_elemNumX * m_elemNumY);
        for (int ei = 0; ei < m_elemNumX; ++ei)
        {
            for (int ej = 0; ej < m_elemNumY; ++ej)
            {
                int index = (ei * m_elemNumY) + ej;
                sendBuffer[index] = m_totalDeltaT[ei + m_ghostCellNum][ej + m_ghostCellNum];
            }
        }

        // Send data to rank 0
        MPI_Send(sendBuffer.data(), m_elemNumX * m_elemNumY, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

void CWENOFV::outputResults(int count, double now, std::string specialCase)
{
    MPI_Barrier(MPI_COMM_WORLD);

    if (specialCase == "initial")
    {
        if (m_rank == 0)
            std::cout << "Start iteration..." << std::endl;

        outputAve("initial");
        if (m_testcase == SMOOTH || m_testcase == VORTEX)
        {
            outputAccuracy("initial", now);
            outputAccuracyAve("initial", now);
        }

        if (m_testcase == SHARP)
            outputUPoscillation(sc_common::intToString(count) + "_" + sc_common::doubleToString(now));
    }
    else if (specialCase == "final")
    {
        if (m_rank == 0)
            std::cout << "\nTotal iterations completed. " << "Total iterations: " << count << std::endl;

        outputAve("final");
        if (m_testcase == SMOOTH || m_testcase == VORTEX)
        {
            outputAccuracy("final", now);
            outputAccuracyAve("final", now);
        }

        outputSpecialTestCasebyCell("final");
    }
    else if (specialCase == "intermediate")
    {
        if (m_rank == 0)
            std::cout << "\nTotal iterations completed. " << "Total iterations: " << count << std::endl;

        outputAve(sc_common::doubleToString(now) + "_intermediate");
    }
    else if (count == 1)
    {
        outputSpecialTestCasebyTime(now);
        outputAve(sc_common::intToString(count) + "_" + sc_common::doubleToString(now));
        // outputPerRankAveWithGhost(sc_common::intToString(count) + "_" + sc_common::doubleToString(now));

        if (m_testcase == SHARP)
            outputUPoscillation(sc_common::intToString(count) + "_" + sc_common::doubleToString(now));

        if (m_testcase == SMOOTH || m_testcase == VORTEX)
            outputAccuracy("count", now);
    }
    else
    {
        if (count % 200 == 0)
            outputSpecialTestCasebyTime(now);

        if (count % 1000 == 0)
        {
            outputAve(sc_common::intToString(count) + "_" + sc_common::doubleToString(now));
            // outputPerRankAveWithGhost(sc_common::intToString(count) + "_" + sc_common::doubleToString(now));
        }

        // if (count % 1000 == 0)
        //     outputDeltaT(sc_common::intToString(count) + "_" + sc_common::doubleToString(now));

        // if (count % 100 == 0 && m_testcase == SHARP)
        //     outputUPoscillation(sc_common::intToString(count) + "_" + sc_common::doubleToString(now));

        // if (fabs(now - m_timeMoments.front()) < 1e-10)
        // {
        //     outputAve(sc_common::intToString(count) + "_" + sc_common::doubleToString(now));
        //     m_timeMoments.erase(m_timeMoments.begin());
        // }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
void CWENOFV::outputAve(string prefix)
{
    getWorldTheta();

    m_outputTimer.pause();
    getWorldUh();
    m_outputTimer.start();

    getWorldPr();

    // output
    if (m_rank == 0)
    {
        std::cout << "\noutputing average solution..." << std::endl;

        int VitalVarNum = equation->getVitalVarNum();
        Array1D<double> VitalVar(VitalVarNum);
        Array1D<double> VitalVarAve(VitalVarNum);

        // Names of variables to be output
        Array1D<std::string> VitalVarName(VitalVarNum);
        equation->getVitalVarName(VitalVarName);

        // Output file name
        string filename = m_outputDir + "average_" + prefix + ".plt";
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            outputFile << "TITLE=FvSolution" << std::endl;
            outputFile << "VARIABLES=";
            outputFile << "X ";
            outputFile << "Y ";
            for (int k = 0; k != VitalVarNum; k++)
                outputFile << VitalVarName[k] << " ";
            outputFile << "theta" << " ";
            outputFile << "Pr" << " ";

            outputFile << std::endl;
            outputFile << "ZONE T=TA, ";
            outputFile << "I="; // Y - direction
            outputFile << m_worldElemNumY << ", ";
            outputFile << "J="; // X - direction
            outputFile << m_worldElemNumX << ", ";
            outputFile << "DATAPACKING = POINT" << std::endl;

            for (int ei = m_ghostCellNum; ei != m_worldElemNumX + m_ghostCellNum; ei++)
            {
                for (int ej = m_ghostCellNum; ej != m_worldElemNumY + m_ghostCellNum; ej++)
                {
                    Array1D<double> Uhh(m_varNum);
                    for (int r = 0; r != m_varNum; r++)
                        Uhh[r] = m_worldUh[ei][ej].vector[r];

                    // Output the values of the required variables one by one
                    equation->getVitalVarVal(Uhh, VitalVar);

                    outputFile << m_worldGrids[ei][ej].m_xCenter << " ";
                    outputFile << m_worldGrids[ei][ej].m_yCenter << " ";

                    for (int k = 0; k != VitalVarNum; k++)
                        outputFile << VitalVar[k] << " ";
                    outputFile << m_worldTotalTheta[ei][ej] << " ";
                    outputFile << m_worldTotalPr[ei][ej] << " ";
                    outputFile << std::endl;
                }
            }

            outputFile.close();
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
        std::cout << "The file " << filename << " has been output successfully..." << std::endl;
    }
}
void CWENOFV::outputPerRankAveWithGhost(string prefix)
{
    std::cout << "\noutputing average solution..." << std::endl;

    int VitalVarNum = equation->getVitalVarNum();
    Array1D<double> VitalVar(VitalVarNum);
    Array1D<double> VitalVarAve(VitalVarNum);

    // Names of variables to be output
    Array1D<std::string> VitalVarName(VitalVarNum);
    equation->getVitalVarName(VitalVarName);

    // Output file name
    string filename = m_outputDir + "average_" + prefix + "rank_" + sc_common::intToString(m_rank) + ".plt";
    std::ofstream outputFile(filename);

    if (outputFile.is_open())
    {
        outputFile << "TITLE=FvSolution" << std::endl;
        outputFile << "VARIABLES=";
        outputFile << "X ";
        outputFile << "Y ";
        for (int k = 0; k != VitalVarNum; k++)
            outputFile << VitalVarName[k] << " ";
        outputFile << "theta" << " ";
        outputFile << "Pr" << " ";

        outputFile << std::endl;
        outputFile << "ZONE T=TA, ";
        outputFile << "I="; // Y - direction
        outputFile << m_elemNumY + 2 * m_ghostCellNum << ", ";
        outputFile << "J="; // X - direction
        outputFile << m_elemNumX + 2 * m_ghostCellNum << ", ";
        outputFile << "DATAPACKING = POINT" << std::endl;

        for (int ei = 0; ei != m_totalElemNumX; ei++)
        {
            for (int ej = 0; ej != m_totalElemNumY; ej++)
            {
                Array1D<double> Uhh(m_varNum);
                for (int r = 0; r != m_varNum; r++)
                    Uhh[r] = m_Uh[ei][ej].vector[r];

                // Output the values of the required variables one by one
                equation->getVitalVarVal(Uhh, VitalVar);

                outputFile << m_grids[ei][ej].m_xCenter << " ";
                outputFile << m_grids[ei][ej].m_yCenter << " ";

                for (int k = 0; k != VitalVarNum; k++)
                    outputFile << VitalVar[k] << " ";
                outputFile << m_totalTheta[ei][ej] << " ";
                outputFile << m_totalPr[ei][ej] << " ";
                outputFile << std::endl;
            }
        }

        outputFile.close();
    }
    else
    {
        std::cout << "Unable to open file" << std::endl;
    }
    std::cout << "The file " << filename << " has been output successfully..." << std::endl;
}

void CWENOFV::outputSliceAve(string prefix)
{
    m_outputTimer.pause();
    getWorldTheta();
    getWorldUh();
    getWorldPr();
    m_outputTimer.start();

    // output
    if (m_rank == 0)
    {
        std::cout << "outputing slice solution..." << std::endl;

        int VitalVarNum = equation->getVitalVarNum();
        Array1D<double> VitalVar(VitalVarNum);
        Array1D<double> VitalVarAve(VitalVarNum);

        // Names of variables to be output
        Array1D<std::string> VitalVarName(VitalVarNum);
        equation->getVitalVarName(VitalVarName);

        // Output file name
        string filename = m_outputDir + "slice_" + prefix + ".plt";
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            outputFile << "TITLE=FvSolution_Slice" << std::endl;
            outputFile << "VARIABLES=";
            outputFile << "X ";
            for (int k = 0; k != VitalVarNum; k++)
                outputFile << VitalVarName[k] << " ";
            outputFile << "theta" << " ";
            outputFile << "Pr" << " ";

            outputFile << std::endl;
            outputFile << "ZONE T=TA, ";
            outputFile << "I="; // X - direction
            outputFile << m_worldElemNumX << ", ";
            outputFile << "DATAPACKING = POINT" << std::endl;

            // 固定 ej = 1
            const int ej = m_ghostCellNum + 1; // ej = 1 对应的实际索引

            for (int ei = m_ghostCellNum; ei != m_worldElemNumX + m_ghostCellNum; ei++)
            {
                Array1D<double> Uhh(m_varNum);
                for (int r = 0; r != m_varNum; r++)
                    Uhh[r] = m_worldUh[ei][ej].vector[r];

                // Output the values of the required variables one by one
                equation->getVitalVarVal(Uhh, VitalVar);

                outputFile << m_worldGrids[ei][ej].m_xCenter << " ";

                for (int k = 0; k != VitalVarNum; k++)
                    outputFile << VitalVar[k] << " ";
                outputFile << m_worldTotalTheta[ei][ej] << " ";
                outputFile << m_worldTotalPr[ei][ej] << " ";
                outputFile << std::endl;
            }

            outputFile.close();
            std::cout << "The slice file " << filename << " has been output successfully..." << std::endl;
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
    }
}
void CWENOFV::outputUPoscillation(string prefix)
{
    m_outputTimer.pause();
    getWorldTheta();
    getWorldUh();
    getWorldPr();
    m_outputTimer.start();

    // output
    if (m_rank == 0)
    {
        std::cout << "outputing UPoscillation..." << std::endl;

        // Output file name
        string filename = m_outputDir + "UPoscillation_" + prefix + ".plt";
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            outputFile << "TITLE=FvSolution" << std::endl;
            outputFile << "VARIABLES=";
            outputFile << "X ";
            outputFile << "Y ";
            outputFile << "Speed_oscillation" << " ";
            outputFile << "P_oscillation" << " ";
            outputFile << "U_oscillation" << " ";
            outputFile << "V_oscillation" << " ";

            outputFile << "\nZONE T=TA, ";
            outputFile << "I="; // Y - direction
            outputFile << m_worldElemNumY << ", ";
            outputFile << "J="; // X - direction
            outputFile << m_worldElemNumX << ", ";
            outputFile << "DATAPACKING = POINT" << std::endl;

            for (int ei = m_ghostCellNum; ei != m_worldElemNumX + m_ghostCellNum; ei++)
            {
                for (int ej = m_ghostCellNum; ej != m_worldElemNumY + m_ghostCellNum; ej++)
                {
                    Array1D<double> Uhh(m_varNum);
                    for (int r = 0; r != m_varNum; r++)
                        Uhh[r] = m_worldUh[ei][ej].vector[r];

                    double u = Uhh[2] / (Uhh[0] + Uhh[1]);
                    double v = Uhh[3] / (Uhh[0] + Uhh[1]);

                    outputFile << m_worldGrids[ei][ej].m_xCenter << " ";
                    outputFile << m_worldGrids[ei][ej].m_yCenter << " ";

                    outputFile << fabs(sqrt(pow(u - 1, 2) + pow(v - 1, 2))) / equation->g_sound1 << " ";
                    outputFile << fabs(equation->getPressure(Uhh)) / equation->g_rho10 / pow(equation->g_sound1, 2) << " ";
                    outputFile << fabs(u - 1) / equation->g_sound1 << " ";
                    outputFile << fabs(v - 1) / equation->g_sound1 << " ";

                    outputFile << std::endl;
                }
            }

            outputFile.close();
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
        std::cout << "The file " << filename << " has been output successfully..." << std::endl;
    }
}
void CWENOFV::outputDeltaT(string prefix)
{
    getWorldDeltaT();

    // output
    if (m_rank == 0)
    {
        std::cout << "\noutputing DeltaT..." << std::endl;

        // Output file name
        string filename = m_outputDir + "DeltaT_" + prefix + ".plt";
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            outputFile << "TITLE=FvSolution" << std::endl;
            outputFile << "VARIABLES=";
            outputFile << "X ";
            outputFile << "Y ";
            outputFile << "DeltaT" << " ";
            outputFile << "nanDeltaT" << " ";
            outputFile << std::endl;

            outputFile << "ZONE T=TA, ";
            outputFile << "I="; // Y - direction
            outputFile << m_worldElemNumY << ", ";
            outputFile << "J="; // X - direction
            outputFile << m_worldElemNumX << ", ";
            outputFile << "DATAPACKING = POINT" << std::endl;

            for (int ei = m_ghostCellNum; ei != m_worldElemNumX + m_ghostCellNum; ei++)
            {
                for (int ej = m_ghostCellNum; ej != m_worldElemNumY + m_ghostCellNum; ej++)
                {

                    outputFile << m_worldGrids[ei][ej].m_xCenter << " ";
                    outputFile << m_worldGrids[ei][ej].m_yCenter << " ";

                    if (isnan(m_worldTotalDeltaT[ei][ej]))
                        outputFile << 0 << " ";
                    else
                        outputFile << m_worldTotalDeltaT[ei][ej] << " ";

                    outputFile << isnan(m_worldTotalDeltaT[ei][ej]) << " ";

                    outputFile << std::endl;
                }
            }

            outputFile.close();
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
        std::cout << "The file " << filename << " has been output successfully..." << std::endl;
    }
}

void CWENOFV::outputUPoscillationByTime(double now)
{
    m_outputTimer.pause();
    getWorldTheta();
    getWorldUh();
    getWorldPr();
    m_outputTimer.start();

    if (m_rank == 0)
    {

        std::cout << "\noutputting UPoscillation..." << std::endl;
        const std::string filename = m_outputDir + "oscillationByTime_" + sc_common::intToString(m_worldElemNumX) + ".plt";
        const bool fileExists = std::filesystem::exists(filename);
        std::ofstream outputFile(filename, std::ios::app);

        if (!fileExists)
        {
            outputFile << "TITLE=FvSolution" << std::endl;
            outputFile << "VARIABLES=";
            outputFile << "Time ";
            outputFile << "Speed_oscillation" << " ";
            outputFile << "P_oscillation" << " ";
            outputFile << "U_oscillation" << " ";
            outputFile << "V_oscillation" << " ";
            outputFile << "\nZONE T=TA" << std::endl;
        }

        if (!outputFile.is_open())
        {
            std::cerr << "Failed to open output file." << std::endl;
            return;
        }
        if (outputFile.is_open())
        {
            // double u_oscillation(0), pre_oscillation(0);

            // for (int e = m_startElemX; e < m_endElemX; ++e)
            // {
            //     Array1D<double> cellAverage(m_varNum);
            //     for (int r = 0; r != m_varNum; r++)
            //         cellAverage[r] = m_Uh[e][r];

            //     const double u_oscillationLoc = fabs(cellAverage[2] / (cellAverage[0] + cellAverage[1]) - 1) / equation->g_sound1;
            //     const double pre_oscillationLoc = fabs(equation->getPressure(cellAverage)) / (equation->g_rho10 * pow(equation->g_sound1, 2));

            //     u_oscillation = max(u_oscillation, u_oscillationLoc);
            //     pre_oscillation = max(pre_oscillation, pre_oscillationLoc);

            double pre_oscillation_MAX(0), speed_oscillation_MAX(0), u_oscillation_MAX(0), v_oscillation_MAX(0);
            for (int ei = m_ghostCellNum; ei != m_worldElemNumX + m_ghostCellNum; ei++)
            {
                for (int ej = m_ghostCellNum; ej != m_worldElemNumY + m_ghostCellNum; ej++)
                {
                    Array1D<double> Uhh(m_varNum);
                    for (int r = 0; r != m_varNum; r++)
                        Uhh[r] = m_worldUh[ei][ej].vector[r];

                    const double pre_oscillation = fabs(equation->getPressure(Uhh)) / equation->g_rho10 / pow(equation->g_sound1, 2);
                    const double speed_oscillation = sqrt(pow(Uhh[2] / (Uhh[0] + Uhh[1]) - 1, 2) + pow(Uhh[3] / (Uhh[0] + Uhh[1]) - 1, 2)) / equation->g_sound1;
                    const double u_oscillation = (Uhh[2] / (Uhh[0] + Uhh[1]) - 1) / equation->g_sound1;
                    const double v_oscillation = (Uhh[3] / (Uhh[0] + Uhh[1]) - 1) / equation->g_sound1;

                    // get max oscillation
                    pre_oscillation_MAX = max(pre_oscillation_MAX, pre_oscillation);
                    speed_oscillation_MAX = max(speed_oscillation_MAX, speed_oscillation);
                    u_oscillation_MAX = max(u_oscillation_MAX, u_oscillation);
                    v_oscillation_MAX = max(v_oscillation_MAX, v_oscillation);
                }
            }
            outputFile << now << " " << pre_oscillation_MAX << " " << speed_oscillation_MAX << " " << u_oscillation_MAX << " " << v_oscillation_MAX << std::endl;
        }

        outputFile.close();
    }
}

void CWENOFV::outputPressureByTime(double now)
{
    m_outputTimer.pause();
    getWorldTheta();
    getWorldUh();
    getWorldPr();
    m_outputTimer.start();

    if (m_rank == 0)
    {
        std::cout << "\noutputting Pressure..." << std::endl;
        const std::string filename = m_outputDir + "pressureByTime.plt";
        const bool fileExists = std::filesystem::exists(filename);
        std::ofstream outputFile(filename, std::ios::app);

        if (!fileExists)
        {
            outputFile << "TITLE=FvSolution" << std::endl;
            outputFile << "VARIABLES=Time P_y3mm P_y15mm P_y30mm P_y80mm\n";
            outputFile << "\nZONE T=TA" << std::endl;
        }

        if (!outputFile.is_open())
        {
            std::cerr << "Failed to open output file." << std::endl;
            return;
        }

        int ei = m_ghostCellNum;

        // 需要输出的 y 位置（单位：米）
        std::vector<double> y_positions_mm = {3.0, 15.0, 30.0, 80.0};
        std::vector<double> pressures;

        for (double y_mm : y_positions_mm)
        {
            double y_m = y_mm / 1000.0;
            int ej = std::floor(y_m / m_deltaY) + m_ghostCellNum;

            // 获取5x5的模板
            Array2D<double> template_5x5(5, 5);
            for (int i = -2; i <= 2; ++i)
            {
                for (int j = -2; j <= 2; ++j)
                {
                    Array1D<double> Uhh(m_varNum);
                    for (int r = 0; r != m_varNum; r++)
                        Uhh[r] = m_worldUh[ei + i][ej + j].vector[r];

                    template_5x5[i + 2][j + 2] = equation->getPressure(Uhh) / equation->g_rho10 / 9.81 / 0.3;
                }
            }

            // x 方向上重构
            Array1D<double> pressuresAtX0(5);
            for (int i = 0; i < 5; ++i)
            {
                pressuresAtX0[i] = WENO5threconstruction(template_5x5[0][i], template_5x5[1][i], template_5x5[2][i], template_5x5[3][i], template_5x5[4][i], 0);
            }

            // y 方向上构造4个高斯点
            Array1D<double> pressuresAtY(4);
            for (int i = 0; i < 4; ++i)
            {
                pressuresAtY[i] = WENO5threconstruction(pressuresAtX0[0], pressuresAtX0[1], pressuresAtX0[2], pressuresAtX0[3], pressuresAtX0[4], i);
            }

            // 选出最接近的一个高斯点
            Array1D<double> loc_gpoints_ref(4), loc_gweights_ref(4);
            sc_math::GaussLobatto_ref(4, loc_gpoints_ref, loc_gweights_ref);

            // 计算该单元内各高斯点的物理y坐标
            // m_worldGrids[ei][ej].m_yCenter 为当前单元的中心位置，m_deltaY为单元高度
            Array1D<double> pointY(4);
            for (int gpj = 0; gpj < 4; ++gpj)
            {
                pointY[gpj] = m_worldGrids[ei][ej].m_yCenter + 0.5 * m_deltaY * loc_gpoints_ref[gpj];
            }

            // 选出最接近目标位置的高斯点
            double min_diff = std::numeric_limits<double>::max();
            int closest_gp_idx = 0;
            for (int gpj = 0; gpj < 4; ++gpj)
            {
                double diff = std::fabs(pointY[gpj] - y_m);
                if (diff < min_diff)
                {
                    min_diff = diff;
                    closest_gp_idx = gpj;
                }
            }

            // 最终取该高斯点对应的压强值
            double pre = pressuresAtY[closest_gp_idx];

            pressures.push_back(pre);
        }

        // 输出时间和各位置压强
        outputFile << now * sqrt(9.81 / 0.3);
        for (const auto &p : pressures)
            outputFile << " " << p;
        outputFile << std::endl;

        outputFile.close();
    }
}
void CWENOFV::outputThincknessByTime(double now)
{
    m_outputTimer.pause();
    getWorldUh();
    m_outputTimer.start();

    if (m_rank == 0)
    {
        std::cout << "\noutputting Thinckness..." << std::endl;
        const std::string filename = m_outputDir + "thincknessByTime.plt";
        const bool fileExists = std::filesystem::exists(filename);
        std::ofstream outputFile(filename, std::ios::app);

        if (!fileExists)
        {
            outputFile << "TITLE=FvSolution" << std::endl;
            outputFile << "VARIABLES=time thinckness\n";
            outputFile << "\nZONE T=TA" << std::endl;
        }

        if (!outputFile.is_open())
        {
            std::cerr << "Failed to open output file." << std::endl;
            return;
        }

        Array2D<double> z1World(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);

        for (int ei = 0; ei != m_worldElemNumX + 2 * m_ghostCellNum; ei++)
        {
            for (int ej = 0; ej != m_worldElemNumY + 2 * m_ghostCellNum; ej++)
            {
                Array1D<double> Uhh(m_varNum);
                for (int r = 0; r != m_varNum; r++)
                    Uhh[r] = m_worldUh[ei][ej].vector[r];

                z1World[ei][ej] = equation->getVolumeFrac(Uhh);
            }
        }

        // 计算最大z1差
        double maxDiff = 0;
        for (int ei = m_ghostCellNum; ei != m_worldElemNumX + m_ghostCellNum; ei++)
        {
            for (int ej = m_ghostCellNum; ej != m_worldElemNumY + m_ghostCellNum; ej++)
            {
                double diff1 = fabs(z1World[ei][ej] - z1World[ei + 1][ej]);
                double diff2 = fabs(z1World[ei][ej] - z1World[ei][ej + 1]);

                maxDiff = max(maxDiff, max(diff1, diff2));
            }
        }

        // 输出时间和各位置压强
        outputFile << now << " " << 1.0 / maxDiff << std::endl;

        outputFile.close();
    }
}

void CWENOFV::outputInterfaceLocationY(int count, double now)
{
    MPI_Barrier(MPI_COMM_WORLD);
    getWorldGridFlux();

    double yL = 0;
    double yR = 0;

    if (m_rank == 0)
    {
        Array1D<double> Uave(m_varNum);
        double z1(0);
        double minDistL(1e4), minDistR(1e4);

        // loop over each element
        for (int ej = m_ghostCellNum; ej != m_worldElemNumY + m_ghostCellNum; ej++)
        {
            for (int gp = 0; gp != m_gpNum; gp++)
            {
                // get Uh
                for (int r = 0; r != m_varNum; ++r)
                    Uave[r] = m_worldGridFlux[m_ghostCellNum][ej].leftFlux[r][gp];

                // get z1_L
                z1 = equation->getVolumeFrac(Uave);

                // check left interface location
                if (fabs(z1 - 0.5) <= minDistL)
                {
                    minDistL = fabs(z1 - 0.5);
                    yL = m_worldGrids[m_ghostCellNum][ej].m_yCenter + 0.5 * m_worldGrids[m_ghostCellNum][ej].m_yDistance * m_gpoints_ref[gp];
                };
            }

            //////////////////////////////////////////////////////////////////////////
            for (int gp = 0; gp != m_gpNum; gp++)
            {
                // get Uh
                for (int r = 0; r != m_varNum; ++r)
                    Uave[r] = m_worldGridFlux[m_ghostCellNum + m_worldElemNumX - 1][ej].rightFlux[r][gp];

                // get z1_R
                z1 = equation->getVolumeFrac(Uave);

                // check left interface location
                if (fabs(z1 - 0.5) <= minDistR)
                {
                    minDistR = fabs(z1 - 0.5);
                    yR = m_worldGrids[m_ghostCellNum + m_worldElemNumX - 1][ej].m_yCenter + 0.5 * m_worldGrids[m_ghostCellNum + m_worldElemNumX - 1][ej].m_yDistance * m_gpoints_ref[gp];
                };
            }
        }

        std::fstream fileout;
        string filename = m_outputDir + "interface_location_by_time.plt";
        fileout.open(filename.c_str(), ios::app);
        if (count == 1)
        {
            fileout << "TITLE=AbgrallErrorByTime" << std::endl;
            fileout << "VARIABLES=" << "time" << " , " << "yL" << " , " << "yR" << std::endl;
        }
        fileout << now << "  " << setprecision(16) << setw(20) << setiosflags(ios::scientific) << yL << " " << yR << std::endl;
        fileout.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
void CWENOFV::outputSymmetryErr(double now)
{
    m_outputTimer.pause();
    getWorldUh();
    m_outputTimer.start();

    Array1D<double> uave(m_varNum), uave_st(m_varNum);
    double z1, z1_st, err;

    double maxErr;

    if (m_rank == 0)
    {
        std::cout << "\noutputting symmetry err..." << std::endl;

        // 创建一个输出文件流
        std::ofstream outputFile(m_outputDir + "symmetry_err.plt", std::ios::app);
        if (!outputFile.is_open())
        {
            std::cerr << "Failed to open output file." << std::endl;
            return; // 或者其他错误处理
        }

        for (int ei = m_ghostCellNum; ei != m_worldElemNumX / 2 + m_ghostCellNum; ei++)
        {
            for (int ej = m_ghostCellNum; ej != m_worldElemNumY + m_ghostCellNum; ej++)
            {
                uave = m_worldUh[ei][ej].vector;
                uave_st = m_worldUh[m_worldElemNumX + 2 * m_ghostCellNum - 1 - ei][ej].vector;

                z1 = equation->getVolumeFrac(uave);
                z1_st = equation->getVolumeFrac(uave_st);

                err = fabs(z1 - z1_st);

                maxErr = max(err, maxErr);
            }
        }
        // 将误差写入文件
        outputFile << now << " " << maxErr << std::endl;
    }
}
void CWENOFV::outputBoundError(string prefix)
{
    getWorldGridFlux();

    // output
    if (m_rank == 0)
    {
        std::cout << "outputing bound error..." << std::endl;

        // Output file name
        const string filename = m_outputDir + "boundError_" + prefix + ".plt";
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            outputFile << "TITLE=FvSolution" << std::endl;
            outputFile << "VARIABLES=";
            outputFile << "X ";
            outputFile << "Y ";
            outputFile << "Error" << " ";

            outputFile << std::endl;
            outputFile << "ZONE T=TA, ";
            outputFile << "I="; // Y - direction
            outputFile << m_worldElemNumY << ", ";
            outputFile << "J="; // X - direction
            outputFile << m_worldElemNumX << ", ";
            outputFile << "DATAPACKING = POINT" << std::endl;

            for (int ei = m_ghostCellNum; ei != m_worldElemNumX + m_ghostCellNum; ei++)
            {
                for (int ej = m_ghostCellNum; ej != m_worldElemNumY + m_ghostCellNum; ej++)
                {
                    outputFile << m_worldGrids[ei][ej].m_xCenter << " ";
                    outputFile << m_worldGrids[ei][ej].m_yCenter << " ";

                    bool flag = 0;

                    for (int gp = 0; gp != m_gpNum; gp++)
                    {
                        if (m_worldGridFlux[ei][ej].bottomFlux[0][gp] < 0 || m_worldGridFlux[ei][ej].bottomFlux[1][gp] < 0)
                            flag = 1;
                        if (m_worldGridFlux[ei][ej].topFlux[0][gp] < 0 || m_worldGridFlux[ei][ej].topFlux[1][gp] < 0)
                            flag = 1;
                        if (m_worldGridFlux[ei][ej].leftFlux[0][gp] < 0 || m_worldGridFlux[ei][ej].leftFlux[1][gp] < 0)
                            flag = 1;
                        if (m_worldGridFlux[ei][ej].rightFlux[0][gp] < 0 || m_worldGridFlux[ei][ej].rightFlux[1][gp] < 0)
                            flag = 1;
                    }

                    if (flag)
                        outputFile << 1 << " ";
                    else
                        outputFile << 0 << " ";
                    outputFile << std::endl;
                }
            }

            outputFile.close();
        }
        else
            std::cout << "Unable to open file" << std::endl;

        std::cout << "The file " << filename << " has been output successfully..." << std::endl;
    }
}
void CWENOFV::outputAccuracy(string prefix, double now)
{
    m_outputTimer.pause();
    getWorldUh();
    m_outputTimer.start();

    if (m_rank == 0)
    {
        std::cout << "Verifying accuracy..." << std::endl;

        // Gauss-Legendre points and weights
        const int loc_gpNum = 4;
        double gPointX, gWeightX;
        double gPointY, gWeightY;
        Array1D<double> loc_gpoints_ref(loc_gpNum), loc_gweights_ref(loc_gpNum);
        sc_math::GaussLegendre_ref(loc_gpNum, loc_gpoints_ref, loc_gweights_ref);

        // Arrays to store exact solution and error
        Array2D<double> Uexact(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);
        Array2D<double> Uerr(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);

        double Uref(0);
        double err_1(0), err_2(0), err_inf(0);

        // Compute exact solution averages over each cell
        Uexact.setZero();
        for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
        {
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int gpi = 0; gpi != loc_gpNum; ++gpi)
                {
                    for (int gpj = 0; gpj != loc_gpNum; ++gpj)
                    {
                        gPointX = m_worldGrids[ei][ej].m_xCenter + 0.5 * m_deltaX * loc_gpoints_ref[gpi];
                        gPointY = m_worldGrids[ei][ej].m_yCenter + 0.5 * m_deltaY * loc_gpoints_ref[gpj];

                        gWeightX = 0.5 * m_deltaX * loc_gweights_ref[gpi];
                        gWeightY = 0.5 * m_deltaY * loc_gweights_ref[gpj];

                        Uref = equation->theVarExact(gPointX, gPointY, now);
                        Uexact[ei][ej] += Uref * gWeightX * gWeightY / m_deltaX / m_deltaY;
                    }
                }
            }
        }

        // Compute error norms
        Array1D<double> Uave(m_varNum);
        for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
        {
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                    Uave[r] = m_worldUh[ei][ej].vector[r];

                // Calculate error at each cell
                Uerr[ei][ej] = fabs(Uexact[ei][ej] - equation->theVarUh(Uave));

                // Update error norms
                err_1 += Uerr[ei][ej];
                err_2 += pow(Uerr[ei][ej], 2);
                if (Uerr[ei][ej] > err_inf)
                    err_inf = Uerr[ei][ej];
            }
        }

        // Calculate average error norms
        err_1 = err_1 / m_worldElemNumX / m_worldElemNumY;
        err_2 = sqrt(err_2 / m_worldElemNumX / m_worldElemNumY);

        // Output to file
        std::fstream fileout;
        string filename = m_outputDir + "accuracy_" + prefix + ".csv";

        // Open file for writing
        std::ifstream fileExists(filename.c_str());
        if (fileExists)
        {
            // File exists, append data
            fileout.open(filename.c_str(), ios::out | ios::app);
        }
        else
        {
            // File doesn't exist, create a new file and write header
            fileout.open(filename.c_str(), ios::out);
            fileout << "Cell Number, Linf-norm, L1-norm, L2-norm" << std::endl;
        }

        const string elemNum = sc_common::intToString(m_elemNumX);
        fileout << m_worldElemNumX << ", ";
        fileout << setprecision(15) << setw(20) << setiosflags(ios::scientific) << err_inf << ", ";
        fileout << setprecision(15) << setw(20) << setiosflags(ios::scientific) << err_1 << ", ";
        fileout << setprecision(15) << setw(20) << setiosflags(ios::scientific) << err_2;

        fileout << std::endl;
        fileout.close();
        std::cout << "Accuracy verification completed..." << std::endl;
    }
}
void CWENOFV::outputAccuracyAve(string prefix, double now)
{
    m_outputTimer.pause();
    getWorldUh();
    m_outputTimer.start();

    if (m_rank == 0)
    {
        std::cout << "Verifying accuracy..." << std::endl;

        // Gauss-Legendre points and weights
        const int loc_gpNum = 4;

        Array1D<double> loc_gpoints_ref(loc_gpNum), loc_gweights_ref(loc_gpNum);
        sc_math::GaussLegendre_ref(loc_gpNum, loc_gpoints_ref, loc_gweights_ref);

        // Arrays to store exact solution and error
        Array2D<double> Uexact(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);
        Array2D<double> Uerr(m_worldElemNumX + 2 * m_ghostCellNum, m_worldElemNumY + 2 * m_ghostCellNum);

        // Compute exact solution averages over each cell
        Uexact.setZero();
        for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
        {
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int gpi = 0; gpi != loc_gpNum; ++gpi)
                {
                    for (int gpj = 0; gpj != loc_gpNum; ++gpj)
                    {
                        const double gPointX = m_worldGrids[ei][ej].m_xCenter + 0.5 * m_deltaX * loc_gpoints_ref[gpi];
                        const double gPointY = m_worldGrids[ei][ej].m_yCenter + 0.5 * m_deltaY * loc_gpoints_ref[gpj];

                        const double gWeightX = 0.5 * m_deltaX * loc_gweights_ref[gpi];
                        const double gWeightY = 0.5 * m_deltaY * loc_gweights_ref[gpj];

                        const double Uref = equation->theVarExact(gPointX, gPointY, now);
                        Uexact[ei][ej] += Uref * gWeightX * gWeightY / m_deltaX / m_deltaY;
                    }
                }
            }
        }

        // Compute error norms
        Array1D<double> Uave(m_varNum);
        for (int ei = m_worldStartElemX; ei != m_worldEndElemX; ei++)
        {
            for (int ej = m_worldStartElemY; ej != m_worldEndElemY; ej++)
            {
                for (int r = 0; r != m_varNum; ++r)
                    Uave[r] = m_worldUh[ei][ej].vector[r];

                // Calculate error at each cell
                Uerr[ei][ej] = fabs(Uexact[ei][ej] - equation->theVarUh(Uave));
            }
        }

        ////////

        std::cout << "\noutputing average solution for yanzheng accuracy..." << std::endl;

        // Output file name
        string filename = m_outputDir + "average_accuracy_" + prefix + ".plt";
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            outputFile << "TITLE=FvSolution" << std::endl;
            outputFile << "VARIABLES=";
            outputFile << "X ";
            outputFile << "Y ";
            outputFile << "num" << " ";
            outputFile << "real" << " ";
            outputFile << "err" << " ";

            outputFile << std::endl;
            outputFile << "ZONE T=TA, ";
            outputFile << "I="; // Y - direction
            outputFile << m_worldElemNumY << ", ";
            outputFile << "J="; // X - direction
            outputFile << m_worldElemNumX << ", ";
            outputFile << "DATAPACKING = POINT" << std::endl;

            for (int ei = m_ghostCellNum; ei != m_worldElemNumX + m_ghostCellNum; ei++)
            {
                for (int ej = m_ghostCellNum; ej != m_worldElemNumY + m_ghostCellNum; ej++)
                {
                    outputFile << m_worldGrids[ei][ej].m_xCenter << " ";
                    outputFile << m_worldGrids[ei][ej].m_yCenter << " ";

                    for (int r = 0; r != m_varNum; ++r)
                        Uave[r] = m_worldUh[ei][ej].vector[r];
                    // Calculate error at each cell
                    Uerr[ei][ej] = fabs(Uexact[ei][ej] - equation->theVarUh(Uave));

                    outputFile << equation->theVarUh(Uave) << " ";
                    outputFile << Uexact[ei][ej] << " ";
                    outputFile << fabs(equation->theVarUh(Uave) - Uexact[ei][ej]) << " ";
                    outputFile << std::endl;
                }
            }

            outputFile.close();
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
        std::cout << "The file " << filename << " has been output successfully..." << std::endl;
    }
}
void CWENOFV::outputSpecialTestCasebyTime(double now)
{
    switch (m_testcase)
    {
    case SHARP:
        outputUPoscillationByTime(now);
        outputThincknessByTime(now);
        break;
    case RTI:
        outputSymmetryErr(now);
        break;
    case DAMBREAK:
        outputPressureByTime(now);
        break;
    default:
        break;
    }
}
void CWENOFV::outputSpecialTestCasebyCell(string prefix)
{
    switch (m_testcase)
    {
    case SHARP:
        outputUPoscillation(prefix);
        break;
    default:
        break;
    }
}
void CWENOFV::backupProjectFiles()
{
    // 创建时间戳目录
    std::filesystem::path codeOutputDir = m_outputDir + "sourceCode/";
    std::filesystem::create_directories(codeOutputDir);
    std::filesystem::create_directories(codeOutputDir / "output/");

    // 需要复制的目录列表
    const std::vector<std::pair<std::filesystem::path, std::filesystem::path>> dirsToCopy = {
        {"./bin", codeOutputDir / "bin"},
        {"./build", codeOutputDir / "build"},
        {"./config", codeOutputDir / "config"},
        {"./docs", codeOutputDir / "docs"},
        {"./include", codeOutputDir / "include"},
        {"./src", codeOutputDir / "src"},
        {"./makefile", codeOutputDir / "makefile"},
        {"./third_party", codeOutputDir / "third_party"}};

    // 执行复制操作
    for (const auto &[src, dst] : dirsToCopy)
    {
        if (std::filesystem::exists(src))
        {
            sc_common::copyDirectory(src, dst);
        }
        else
        {
            std::cerr << "警告: 源目录不存在 " << src << std::endl;
        }
    }
}