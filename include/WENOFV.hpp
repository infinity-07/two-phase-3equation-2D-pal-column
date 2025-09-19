#ifndef __WENOFV__HPP__
#define __WENOFV__HPP__

#include <iostream>
#include <fstream>
#include <map>
#include "../third_party/scTools/scTools.h"
#include "../include/Equations.hpp"
#include "../third_party/scTools/Timer.hpp"

// SSP-RK time integration methods
enum RKMETHOD
{
    RK1 = 0,
    RK2 = 1,
    RK3 = 2,
};

// WENO schemes
enum SCHEMETYPE
{
    WENO = 0,
    WENOZ = 1,
    WENOEXP = 2,
};

// WENO reconstruction types
enum RECONSTRUCTIONTYPE
{
    CONVERVATIVE = 0,
    CHARACTERISTIC = 1,
    PRIMITIVE = 2,
    PRIMITIVE2 = 3,
    PRIMITIVE3 = 4,
    ThincOnlyForZ1 = 5,
};

// 使用 inline 定义字符串数组
inline const std::string SCHEMETYPE_STRINGS[] = {
    "WENO",   // SCHEMETYPE::WENO
    "WENOZ",  // SCHEMETYPE::WENOZ
    "WENOEXP" // SCHEMETYPE::WENOEXP
};

inline const std::string RECONSTRUCTIONTYPE_STRINGS[] = {
    "CONVERVATIVE",   // RECONSTRUCTIONTYPE::CONVERVATIVE
    "CHARACTERISTIC", // RECONSTRUCTIONTYPE::CHARACTERISTIC
    "PRIMITIVE",      // RECONSTRUCTIONTYPE::PRIMITIVE
    "PRIMITIVE2",     // RECONSTRUCTIONTYPE::PRIMITIVE2
    "PRIMITIVE3"      // RECONSTRUCTIONTYPE::PRIMITIVE3
};

//
class Cgrid
{
public:
    // 成员变量
    double m_xLeft;
    double m_xRight;
    double m_xCenter;
    double m_xDistance;

    double m_yLeft;
    double m_yRight;
    double m_yCenter;
    double m_yDistance;

    double m_area;

    // 默认构造函数
    Cgrid()
        : m_xLeft(0), m_xRight(0), m_xCenter(0), m_xDistance(0),
          m_yLeft(0), m_yRight(0), m_yCenter(0), m_yDistance(0), m_area(0)
    {
    }

    // 带参数的构造函数
    Cgrid(double xLeft, double xRight, double yLeft, double yRight)
        : m_xLeft(xLeft), m_xRight(xRight), m_yLeft(yLeft), m_yRight(yRight)
    {
        m_xCenter = (m_xLeft + m_xRight) / 2;
        m_xDistance = m_xRight - m_xLeft;

        m_yCenter = (m_yLeft + m_yRight) / 2;
        m_yDistance = m_yRight - m_yLeft;

        m_area = m_xDistance * m_yDistance;
    }
};

class CgridFlux
{
public:
    // 成员变量
    Array2D<double> rightFlux;
    Array2D<double> leftFlux;
    Array2D<double> topFlux;
    Array2D<double> bottomFlux;
};

class Cvector
{
public:
    // 成员变量
    Array1D<double> vector;
};

class CWENOFV
{
public:
    CWENOFV() {}
    ~CWENOFV();

    CWENOFV(TwoPhaseEquation *equation, std::map<std::string, std::string> option)
    {
        this->equation = equation;
        this->option = option;
    }

    TwoPhaseEquation *equation;

public:
    std::map<std::string, std::string> option;

    int m_procIndexX, m_procIndexY;
    int m_rank, m_size;
    int m_worldElemNumX, m_worldElemNumY;
    int m_numProcessorsX;
    // int m_numProcessorsY;

    int m_ghostCellNum;
    int m_startElemX, m_endElemX;
    int m_startElemY, m_endElemY;

    int m_elemNumX, m_elemNumY;
    int m_totalElemNumX, m_totalElemNumY;

    int m_worldStartElemX, m_worldStartElemY;
    int m_worldEndElemX, m_worldEndElemY;
    int m_varNum;

    double m_globalXL;
    double m_globalXR;
    double m_globalYL;
    double m_globalYR;

    double m_cfl;
    double m_deltaX;
    double m_deltaY;

    int m_gpNum;
    bool m_pplimiter;
    std::vector<double> m_timeMoments;
    Array1D<double> m_gpoints_ref, m_gweights_ref;
    Array2D<Cvector> m_rhs;
    Array2D<Cvector> m_Un;

private:
    Timer m_mainTimer;
    Timer m_wenoTimer;
    Timer m_ppLimiterTimer;
    Timer m_FVTimer;
    Timer m_outputTimer;
    Timer m_solverTimer;
    Timer m_initialiTimer;
    Timer m_riemannTimer;
    Timer m_MPITimer;

private:
public:
    Array2D<Cvector> m_Uh;
    Array2D<Cvector> m_worldUh;
    Array2D<Cgrid> m_grids;
    Array2D<Cgrid> m_worldGrids;
    Array2D<CgridFlux> m_gridFlux;
    Array2D<CgridFlux> m_worldGridFlux;
    Array2D<double> m_worldTotalTheta;
    Array2D<double> m_totalTheta;
    Array2D<double> m_worldTotalDeltaT;
    Array2D<double> m_totalDeltaT;
    Array2D<double> m_worldTotalPr;
    Array2D<double> m_totalPr;

    std::string m_outputDir;

public:
    TESTCASETYPE m_testcase;
    SCHEMETYPE m_scheme;
    RKMETHOD m_rkMeth;
    RIEMANNFLUXTYPE m_riemannFluxType;
    RECONSTRUCTIONTYPE m_reconstructionType;

    void initializeSolver();
    void initializeAve();

    void performTimeIntegration(double timestep);
    void RunRK1(double deltaT);
    void RunRK2(double deltaT);
    void RunRK3(double deltaT);

    void setBoundaryGPs(void);
    // void setBoundaryAverages(void);

    // FVM discretization
    void run();
    void assembleRHS();
    double calculateDeltaT(double now);
    void getFlux_conservative(void);
    void getFlux_characteristic(void);
    void getFlux_primitive(void);
    void getFlux_primitive2(void);
    void getFlux_primitive3(void);

    void assembleBoundaryFaceTerm();
    void assembleSourceTerm();
    double useWENO(double uavemmm, double uavemm, double uavem, double uave, double uavep, double uavepp, double uaveppp, int gp);

    void useBoundPreservingLimiter();
    double getPositivityPreservingTheta(int ei, int ej);
    void useFVLimiter();

    // MPI functions
    void getWorldUh(void);
    // void getWorldGridFlux(void);
    void getWorldTheta(void);
    void getWorldPr(void);
    void getWorldDeltaT(void);
    // void worldToProcGridFlux(void);
    void worldToProcUh(void);

    void exchangeGhostCellsValue(void);
    void exchangeGhostCellsGridFlux(void);
    void setBoundaryAverages2(void);

    void setBoundaryGPs1D(void);

    // Output functions
    void outputResults(int count, double now, std::string specialCase = "");
    void backupProjectFiles(void);
    void outputInterfaceLocationY(int count, double now);
    void outputSpecialTestCasebyTime(double now);
    void outputUPoscillationByTime(double now);
    void outputPressureByTime(double now);
    void outputThincknessByTime(double now);
    void outputSymmetryErr(double now);

    void outputSpecialTestCasebyCell(string prefix);
    void outputBoundError(string prefix);
    void outputUPoscillation(string prefix);
    void outputSliceAve(string prefix);

    void outputAve(std::string filename);
    // void outputPerRankAveWithGhost(std::string filename);
    void outputDeltaT(std::string filename);
    void outputAccuracy(string prefix, double now);
    void outputAccuracyAve(string prefix, double now);
};

///////////////////////////////////////////////////////////////////////
#endif