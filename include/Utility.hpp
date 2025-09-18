#ifndef __UTILITY__HPP__
#define __UTILITY__HPP__
#include "../third_party/scTools/scTools.h"

inline double WENO5threconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp, int gp_num)
{
    double d0(0), d1(0), d2(0);
    double w0(0), w1(0), w2(0);
    double h0(0), h1(0), h2(0);

    // use linear weights
    switch (gp_num)
    {
    case 0:
        d0 = 0.3;
        d1 = 0.6;
        d2 = 0.1;
        break;
    case 1:
        d0 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 2:
        d0 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 3:
        d0 = 0.1;
        d1 = 0.6;
        d2 = 0.3;
        break;
    default:
        std::cout << gp_num << std::endl;
        std::cout << "check the gauss points for weno reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    // use nonlinear weights
    // smooth indicator for current component

    // Jiang, G.-S., & Shu, C.-W. (1996).
    // Efficient Implementation of Weighted ENO Schemes.
    // Journal of Computational Physics, 126(1), 202-228. https://doi.org/10.1006/jcph.1996.0130

    const double beta0 = 0.25 * pow(uavemm - 4 * uavem + 3 * uave, 2) + 13.0 / 12.0 * pow(uavemm - 2 * uavem + uave, 2);
    const double beta1 = 0.25 * pow(uavem - uavep, 2) + 13.0 / 12.0 * pow(uavem - 2 * uave + uavep, 2);
    const double beta2 = 0.25 * pow(3 * uave - 4 * uavep + uavepp, 2) + 13.0 / 12.0 * pow(uave - 2 * uavep + uavepp, 2);

    const double eps(1.0e-6);
    w0 = d0 / ((eps + beta0) * (eps + beta0));
    w1 = d1 / ((eps + beta1) * (eps + beta1));
    w2 = d2 / ((eps + beta2) * (eps + beta2));

    const double sum = w0 + w1 + w2;

    w0 = w0 / sum;
    w1 = w1 / sum;
    w2 = w2 / sum;

    // p(x)
    switch (gp_num)
    {
    case 0:
        h0 = -1.0 / 6.0 * uavemm + 5.0 / 6.0 * uavem + 1.0 / 3.0 * uave;
        h1 = 1.0 / 3.0 * uavem + 5.0 / 6.0 * uave - 1.0 / 6 * uavep;
        h2 = 11.0 / 6.0 * uave - 7.0 / 6.0 * uavep + 1.0 / 3.0 * uavepp;
        break;
    case 1:
        h0 = (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave;
        h1 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60 - sqrt(5.0) / 20.0) * uavep;
        h2 = (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavepp;
        break;
    case 2:
        h0 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave;
        h1 = (-1.0 / 60 - sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavep;
        h2 = (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavepp;
        break;
    case 3:
        h0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
        h1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
        h2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;
        break;
    default:
        std::cout << "Check the Gauss points for WENO reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    const double ph = w0 * h0 + w1 * h1 + w2 * h2;

    return ph;
}
inline void WENO5threconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp, Array1D<double> gpList, Array1D<double> &phList)
{
    double d0(0), d1(0), d2(0);
    double w0(0), w1(0), w2(0);
    double h0(0), h1(0), h2(0);

    int gpNum = gpList.dim();

    // use nonlinear weights
    // smooth indicator for current component

    // Jiang, G.-S., & Shu, C.-W. (1996).
    // Efficient Implementation of Weighted ENO Schemes.
    // Journal of Computational Physics, 126(1), 202-228. https://doi.org/10.1006/jcph.1996.0130

    const double beta0 = 0.25 * pow(uavemm - 4 * uavem + 3 * uave, 2) + 13.0 / 12.0 * pow(uavemm - 2 * uavem + uave, 2);
    const double beta1 = 0.25 * pow(uavem - uavep, 2) + 13.0 / 12.0 * pow(uavem - 2 * uave + uavep, 2);
    const double beta2 = 0.25 * pow(3 * uave - 4 * uavep + uavepp, 2) + 13.0 / 12.0 * pow(uave - 2 * uavep + uavepp, 2);

    for (int gpIndex = 0; gpIndex != gpNum; ++gpIndex)
    {
        int gp = gpList[gpIndex];

        // use linear weights
        switch (gp)
        {
        case 0:
            d0 = 0.3;
            d1 = 0.6;
            d2 = 0.1;
            break;
        case 1:
            d0 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
            d1 = 129.0 / 220.0;
            d2 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
            break;
        case 2:
            d0 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
            d1 = 129.0 / 220.0;
            d2 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
            break;
        case 3:
            d0 = 0.1;
            d1 = 0.6;
            d2 = 0.3;
            break;
        default:
            std::cout << gp << std::endl;
            std::cout << "check the gauss points for weno reconstruction!" << std::endl;
            std::cin.get();
            exit(1);
        }

        // p(x)
        switch (gp)
        {
        case 0:
            h0 = -1.0 / 6.0 * uavemm + 5.0 / 6.0 * uavem + 1.0 / 3.0 * uave;
            h1 = 1.0 / 3.0 * uavem + 5.0 / 6.0 * uave - 1.0 / 6 * uavep;
            h2 = 11.0 / 6.0 * uave - 7.0 / 6.0 * uavep + 1.0 / 3.0 * uavepp;
            break;
        case 1:
            h0 = (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave;
            h1 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60 - sqrt(5.0) / 20.0) * uavep;
            h2 = (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavepp;
            break;
        case 2:
            h0 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave;
            h1 = (-1.0 / 60 - sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavep;
            h2 = (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavepp;
            break;
        case 3:
            h0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
            h1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
            h2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;
            break;
        default:
            std::cout << "Check the Gauss points for WENO reconstruction!" << std::endl;
            std::cin.get();
            exit(1);
        }

        const double eps(1.0e-6);
        w0 = d0 / ((eps + beta0) * (eps + beta0));
        w1 = d1 / ((eps + beta1) * (eps + beta1));
        w2 = d2 / ((eps + beta2) * (eps + beta2));

        const double sum = w0 + w1 + w2;

        w0 = w0 / sum;
        w1 = w1 / sum;
        w2 = w2 / sum;

        phList[gpIndex] = w0 * h0 + w1 * h1 + w2 * h2;
    }
}

inline double WENO5Zthreconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp, int gp_num)
{
    // Borges, R., Carmona, M., Costa, B., & Don, W. S. (2008).
    // An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws.
    // Journal of Computational Physics, 227(6), 3191-3211. https://doi.org/10.1016/j.jcp.2007.11.038

    double d0(0), d1(0), d2(0);
    double p0(0), p1(0), p2(0), ph(0);

    // use linear weights
    switch (gp_num)
    {
    case 0:
        d0 = 0.3;
        d1 = 0.6;
        d2 = 0.1;
        break;
    case 1:
        d0 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 2:
        d0 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 3:
        d0 = 0.1;
        d1 = 0.6;
        d2 = 0.3;
        break;
    default:
        std::cout << gp_num << std::endl;
        std::cout << "check the gauss points for weno reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    // p(x)
    switch (gp_num)
    {
    case 0:
        p0 = -1.0 / 6.0 * uavemm + 5.0 / 6.0 * uavem + 1.0 / 3.0 * uave;
        p1 = 1.0 / 3.0 * uavem + 5.0 / 6.0 * uave - 1.0 / 6 * uavep;
        p2 = 11.0 / 6.0 * uave - 7.0 / 6.0 * uavep + 1.0 / 3.0 * uavepp;
        break;
    case 1:
        p0 = (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave;
        p1 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60 - sqrt(5.0) / 20.0) * uavep;
        p2 = (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavepp;
        break;
    case 2:
        p0 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave;
        p1 = (-1.0 / 60 - sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavep;
        p2 = (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavepp;
        break;
    case 3:
        p0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
        p1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
        p2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;
        break;
    default:
        std::cout << "Check the Gauss points for WENO reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    const double uave_max = std::max(std::max(std::max(uavemm, uavem), std::max(uavepp, uavep)), uave);
    const double uave_min = std::min(std::min(std::min(uavemm, uavem), std::min(uavepp, uavep)), uave);
    // dimentionless
    uavemm = (uavemm - uave_min) / (uave_max - uave_min + 1e-20);
    uavem = (uavem - uave_min) / (uave_max - uave_min + 1e-20);
    uave = (uave - uave_min) / (uave_max - uave_min + 1e-20);
    uavep = (uavep - uave_min) / (uave_max - uave_min + 1e-20);
    uavepp = (uavepp - uave_min) / (uave_max - uave_min + 1e-20);

    // wenoz
    const double beta0 = 0.25 * pow(uavemm - 4 * uavem + 3 * uave, 2) + 13.0 / 12.0 * pow(uavemm - 2 * uavem + uave, 2);
    const double beta1 = 0.25 * pow(uavem - uavep, 2) + 13.0 / 12.0 * pow(uavem - 2 * uave + uavep, 2);
    const double beta2 = 0.25 * pow(3 * uave - 4 * uavep + uavepp, 2) + 13.0 / 12.0 * pow(uave - 2 * uavep + uavepp, 2);

    const double eps(1.0e-10);
    double tau = fabs(beta0 - beta2);
    double w0 = d0 * (1 + tau / (beta0 + eps));
    double w1 = d1 * (1 + tau / (beta1 + eps));
    double w2 = d2 * (1 + tau / (beta2 + eps));

    double sum = w0 + w1 + w2;

    w0 = w0 / sum;
    w1 = w1 / sum;
    w2 = w2 / sum;

    ph = w0 * p0 + w1 * p1 + w2 * p2;

    return ph;
}
inline double wenoExpthreconstruction(double uavemmm, double uavemm, double uavem, double uave, double uavep, double uavepp, double uaveppp, int gp, double m_deltaX)
{
    // 使用lambda参数实现升阶
    // 值得一提的是，在计算矩形波平移时，
    // 由于一开始过于锐利，所以无法准确计算 lambda，因此在前期计算效果很差
    double dif5, dif3;

    if (gp >= 2)
    {
        dif5 = (-uavemm + 5.0 * uavem - 10.0 * uave + 10.0 * uavep - 5.0 * uavepp + uaveppp);
        dif3 = (uavemm - 11.0 * uavem + 28.0 * uave - 28.0 * uavep + 11.0 * uavepp - uaveppp) / 6.0;
    }
    else
    {
        dif5 = (-uavemmm + 5.0 * uavemm - 10.0 * uavem + 10.0 * uave - 5.0 * uavep + uavepp);
        dif3 = (uavemmm - 11.0 * uavemm + 28.0 * uavem - 28.0 * uave + 11.0 * uavep - uavepp) / 6.0;
    }

    double l2h2(0);

    int sw = 0;
    if (fabs(dif3) < pow(m_deltaX, 4))
    // if (fabs(dif3) < 1e-2)
    {
        sw = 1;
        l2h2 = 0;
    }
    else
    {
        l2h2 = dif5 / dif3;
        if (fabs(l2h2) > 1)
            l2h2 = 1;
    }

    // l2h2 = 0.1;

    //  计算每个模板上的重构值
    double h0, h1, h2;

    switch (gp)
    {
    case 0:
        h0 = -1.0 / 6.0 * uavemm + 5.0 / 6.0 * uavem + 1.0 / 3.0 * uave;
        h1 = 1.0 / 3.0 * uavem + 5.0 / 6.0 * uave - 1.0 / 6.0 * uavep;
        h2 = 11.0 / 6.0 * uave - 7.0 / 6.0 * uavep + 1.0 / 3.0 * uavepp;
        break;
    case 1:
        h0 = (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave;
        h1 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60 - sqrt(5.0) / 20.0) * uavep;
        h2 = (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavepp;
        break;
    case 2:
        h0 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave;
        h1 = (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavep;
        h2 = (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavepp;
        break;
    case 3:
        h0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
        h1 = -1.0 / 6.0 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
        h2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;
        break;
    default:
        std::cerr << "Check the Gauss points for WENO reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    /// 线性权重
    double d0, d1, d2;

    switch (gp)
    {
    case 0:
        d0 = 0.3 - 0.057142858 * l2h2 + 0.0071428571 * pow(l2h2, 2) - 0.00077200577 * pow(l2h2, 3);
        d1 = 0.6 + 0.078571429 * l2h2 - 0.0101190476 * pow(l2h2, 2) + 0.00112012987 * pow(l2h2, 3);
        d2 = 0.1 - 0.021428571 * l2h2 + 0.0029761905 * pow(l2h2, 2) - 0.00034812410 * pow(l2h2, 3);

        break;
    case 1:
        d0 = 0.16108042 - 0.035187607 * l2h2 + 0.0049466098 * pow(l2h2, 2) - 0.00058211525 * pow(l2h2, 3);
        d1 = 0.58636364 + 0.087692641 * l2h2 - 0.0120401515 * pow(l2h2, 2) + 0.00139436871 * pow(l2h2, 3);
        d2 = 0.25255594 - 0.052505034 * l2h2 + 0.0070935417 * pow(l2h2, 2) - 0.00081225346 * pow(l2h2, 3);

        break;
    case 2:
        d0 = 0.25255594 - 0.052505034 * l2h2 + 0.0070935417 * pow(l2h2, 2) - 0.00081225346 * pow(l2h2, 3);
        d1 = 0.58636364 + 0.087692641 * l2h2 - 0.0120401515 * pow(l2h2, 2) + 0.00139436871 * pow(l2h2, 3);
        d2 = 0.16108042 - 0.035187607 * l2h2 + 0.0049466098 * pow(l2h2, 2) - 0.00058211525 * pow(l2h2, 3);

        break;
    case 3:
        d0 = 0.1 - 0.021428571 * l2h2 + 0.0029761905 * pow(l2h2, 2) - 0.00034812410 * pow(l2h2, 3);
        d1 = 0.6 + 0.078571429 * l2h2 - 0.0101190476 * pow(l2h2, 2) + 0.00112012987 * pow(l2h2, 3);
        d2 = 0.3 - 0.057142858 * l2h2 + 0.0071428571 * pow(l2h2, 2) - 0.00077200577 * pow(l2h2, 3);

        break;
    default:
        std::cout << gp << std::endl;
        std::cout << "check the gauss points for weno reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    // round to 8 digits
    // double factor = std::pow(10, 8);
    // d0 = std::round(d0 * factor) / factor;
    // d1 = std::round(d1 * factor) / factor;
    // d2 = std::round(d2 * factor) / factor;

    // double theta = 0.25;
    // double beta0 = theta * abs((1.0 - 0) * uavemm + (2.0 * 0 - 3.0) * uavem + (2.0 - 0) * uave) + abs(uavemm - 2.0 * uavem + uave);
    // double beta1 = theta * abs((1.0 - 1.0) * uavem + (2.0 * 1.0 - 3.0) * uave + (2.0 - 1.0) * uavep) + abs(uavem - 2.0 * uave + uavep);
    // double beta2 = theta * abs((1.0 - 2.0) * uave + (2.0 * 2.0 - 3.0) * uavep + (2.0 - 2.0) * uavepp) + abs(uave - 2.0 * uavep + uavepp);

    // double beta0 = 0.25 * pow(uavemm - 4 * uavem + 3 * uave, 2) + 13.0 / 12.0 * pow(uavemm - 2 * uavem + uave, 2);
    // double beta1 = 0.25 * pow(uavem - uavep, 2) + 13.0 / 12.0 * pow(uavem - 2 * uave + uavep, 2);
    // double beta2 = 0.25 * pow(3 * uave - 4 * uavep + uavepp, 2) + 13.0 / 12.0 * pow(uave - 2 * uavep + uavepp, 2);

    const double uave_max = std::max(std::max(std::max(uavemm, uavem), std::max(uave, uavep)), uavepp);
    const double uave_min = std::min(std::min(std::min(uavemm, uavem), std::min(uave, uavep)), uavepp);
    // // dimentionless
    uavemm = (uavemm - uave_min) / (uave_max - uave_min + 1e-20);
    uavem = (uavem - uave_min) / (uave_max - uave_min + 1e-20);
    uave = (uave - uave_min) / (uave_max - uave_min + 1e-20);
    uavep = (uavep - uave_min) / (uave_max - uave_min + 1e-20);
    uavepp = (uavepp - uave_min) / (uave_max - uave_min + 1e-20);

    // WENOEXP
    const double theta = 1;
    const double beta0 = theta * 1.0 / 2.0 * abs(uavemm - 4 * uavem + 3 * uave) + abs(uavemm - 2 * uavem + uave);
    const double beta1 = theta * 1.0 / 2.0 * abs(uavem - uavep) + abs(uavem - 2 * uave + uavep);
    const double beta2 = theta * 1.0 / 2.0 * abs(3 * uave - 4 * uavep + uavepp) + abs(uave - 2 * uavep + uavepp);
    const double tau = uavemm - 4.0 * uavem + 6.0 * uave - 4.0 * uavep + uavepp;
    const double epsilon = 1e-16;
    double w0 = d0 * (1.0 + 2 * pow(tau, 2.0) / (pow(beta0, 2.0) + epsilon));
    double w1 = d1 * (1.0 + 2 * pow(tau, 2.0) / (pow(beta1, 2.0) + epsilon));
    double w2 = d2 * (1.0 + 2 * pow(tau, 2.0) / (pow(beta2, 2.0) + epsilon));

    // wenoz
    // const double beta0 = 0.25 * pow(uavemm - 4 * uavem + 3 * uave, 2) + 13.0 / 12.0 * pow(uavemm - 2 * uavem + uave, 2);
    // const double beta1 = 0.25 * pow(uavem - uavep, 2) + 13.0 / 12.0 * pow(uavem - 2 * uave + uavep, 2);
    // const double beta2 = 0.25 * pow(3 * uave - 4 * uavep + uavepp, 2) + 13.0 / 12.0 * pow(uave - 2 * uavep + uavepp, 2);
    // const double eps(1.0e-10);
    // double tau = fabs(beta0 - beta2);
    // double w0 = d0 * (1 + tau / (beta0 + eps));
    // double w1 = d1 * (1 + tau / (beta1 + eps));
    // double w2 = d2 * (1 + tau / (beta2 + eps));

    const double S = w0 + w1 + w2;

    w0 /= S;
    w1 /= S;
    w2 /= S;

    // 重构
    double ph = w0 * h0 + w1 * h1 + w2 * h2;

    // if (sw == 1)
    // {
    //     switch (gp)
    //     {
    //     case 0:
    //         ph += -1.0 / 60.0 * dif5 * pow(m_deltaX, 5);
    //         break;
    //     case 1:
    //         ph += -385.0 * sqrt(5.0) / 90000.0 * dif5 * pow(m_deltaX, 5);
    //         break;
    //     case 2:
    //         ph += +385.0 * sqrt(5.0) / 90000.0 * dif5 * pow(m_deltaX, 5);
    //         break;
    //     case 3:
    //         ph += +1.0 / 60.0 * dif5 * pow(m_deltaX, 5);
    //         break;
    //     default:
    //         std::cout << "gp error" << std::endl;
    //         std::cin.get();
    //         break;
    //     }
    // }

    return ph;
}
inline double wenoExpWithThincReconstruction(double uavemmm, double uavemm, double uavem, double uave, double uavep, double uavepp, double uaveppp, int gp, double m_deltaX)
{
    // 使用lambda参数实现升阶
    // 值得一提的是，在计算矩形波平移时，
    // 由于一开始过于锐利，所以无法准确计算 lambda，因此在前期计算效果很差
    double dif5, dif3;

    if (gp >= 2)
    {
        dif5 = (-uavemm + 5.0 * uavem - 10.0 * uave + 10.0 * uavep - 5.0 * uavepp + uaveppp);
        dif3 = (uavemm - 11.0 * uavem + 28.0 * uave - 28.0 * uavep + 11.0 * uavepp - uaveppp) / 6.0;
    }
    else
    {
        dif5 = (-uavemmm + 5.0 * uavemm - 10.0 * uavem + 10.0 * uave - 5.0 * uavep + uavepp);
        dif3 = (uavemmm - 11.0 * uavemm + 28.0 * uavem - 28.0 * uave + 11.0 * uavep - uavepp) / 6.0;
    }

    double l2h2(0);
    int flag2 = 0;
    const double L2H2MAX = 1e-3; // 震荡检测阈值

    int sw = 0;
    if (fabs(dif3) < pow(m_deltaX, 4))
    {
        sw = 1;
        flag2 = 0;
        l2h2 = 0; // 3阶导数接近0，鉴定为光滑，但是不能使用WENOEXP，使用WENOZ
    }
    else
    {
        flag2 = 0;
        l2h2 = dif5 / dif3; // 使用WENOEXP
    }

    // if (fabs(l2h2) > L2H2MAX)
    //     flag2 = 1; // 检测到震荡，使用THINC

    if ((1e-6 < uave) & (uave < 1 - 1e-6))
        flag2 = 1;

    if (flag2 == 0)

    // l2h2 = 0.1;
    {
        //  计算每个模板上的重构值
        double h0, h1, h2;

        switch (gp)
        {
        case 0:
            h0 = -1.0 / 6.0 * uavemm + 5.0 / 6.0 * uavem + 1.0 / 3.0 * uave;
            h1 = 1.0 / 3.0 * uavem + 5.0 / 6.0 * uave - 1.0 / 6.0 * uavep;
            h2 = 11.0 / 6.0 * uave - 7.0 / 6.0 * uavep + 1.0 / 3.0 * uavepp;
            break;
        case 1:
            h0 = (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave;
            h1 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60 - sqrt(5.0) / 20.0) * uavep;
            h2 = (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavepp;
            break;
        case 2:
            h0 = (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavemm + (1.0 / 30.0 - sqrt(5.0) / 5.0) * uavem + (59.0 / 60.0 + 3 * sqrt(5.0) / 20) * uave;
            h1 = (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60.0 + sqrt(5.0) / 20.0) * uavep;
            h2 = (59.0 / 60.0 - 3 * sqrt(5.0) / 20.0) * uave + (1.0 / 30.0 + sqrt(5.0) / 5.0) * uavep + (-1.0 / 60.0 - sqrt(5.0) / 20.0) * uavepp;
            break;
        case 3:
            h0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
            h1 = -1.0 / 6.0 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
            h2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;
            break;
        default:
            std::cerr << "Check the Gauss points for WENO reconstruction!" << std::endl;
            std::cin.get();
            exit(1);
        }

        /// 线性权重
        double d0, d1, d2;

        switch (gp)
        {
        case 0:
            d0 = 0.3 - 0.057142858 * l2h2 + 0.0071428571 * pow(l2h2, 2) - 0.00077200577 * pow(l2h2, 3);
            d1 = 0.6 + 0.078571429 * l2h2 - 0.0101190476 * pow(l2h2, 2) + 0.00112012987 * pow(l2h2, 3);
            d2 = 0.1 - 0.021428571 * l2h2 + 0.0029761905 * pow(l2h2, 2) - 0.00034812410 * pow(l2h2, 3);

            break;
        case 1:
            d0 = 0.16108042 - 0.035187607 * l2h2 + 0.0049466098 * pow(l2h2, 2) - 0.00058211525 * pow(l2h2, 3);
            d1 = 0.58636364 + 0.087692641 * l2h2 - 0.0120401515 * pow(l2h2, 2) + 0.00139436871 * pow(l2h2, 3);
            d2 = 0.25255594 - 0.052505034 * l2h2 + 0.0070935417 * pow(l2h2, 2) - 0.00081225346 * pow(l2h2, 3);

            break;
        case 2:
            d0 = 0.25255594 - 0.052505034 * l2h2 + 0.0070935417 * pow(l2h2, 2) - 0.00081225346 * pow(l2h2, 3);
            d1 = 0.58636364 + 0.087692641 * l2h2 - 0.0120401515 * pow(l2h2, 2) + 0.00139436871 * pow(l2h2, 3);
            d2 = 0.16108042 - 0.035187607 * l2h2 + 0.0049466098 * pow(l2h2, 2) - 0.00058211525 * pow(l2h2, 3);

            break;
        case 3:
            d0 = 0.1 - 0.021428571 * l2h2 + 0.0029761905 * pow(l2h2, 2) - 0.00034812410 * pow(l2h2, 3);
            d1 = 0.6 + 0.078571429 * l2h2 - 0.0101190476 * pow(l2h2, 2) + 0.00112012987 * pow(l2h2, 3);
            d2 = 0.3 - 0.057142858 * l2h2 + 0.0071428571 * pow(l2h2, 2) - 0.00077200577 * pow(l2h2, 3);

            break;
        default:
            std::cout << gp << std::endl;
            std::cout << "check the gauss points for weno reconstruction!" << std::endl;
            std::cin.get();
            exit(1);
        }

        // round to 8 digits
        // double factor = std::pow(10, 8);
        // d0 = std::round(d0 * factor) / factor;
        // d1 = std::round(d1 * factor) / factor;
        // d2 = std::round(d2 * factor) / factor;

        // double theta = 0.25;
        // double beta0 = theta * abs((1.0 - 0) * uavemm + (2.0 * 0 - 3.0) * uavem + (2.0 - 0) * uave) + abs(uavemm - 2.0 * uavem + uave);
        // double beta1 = theta * abs((1.0 - 1.0) * uavem + (2.0 * 1.0 - 3.0) * uave + (2.0 - 1.0) * uavep) + abs(uavem - 2.0 * uave + uavep);
        // double beta2 = theta * abs((1.0 - 2.0) * uave + (2.0 * 2.0 - 3.0) * uavep + (2.0 - 2.0) * uavepp) + abs(uave - 2.0 * uavep + uavepp);

        // double beta0 = 0.25 * pow(uavemm - 4 * uavem + 3 * uave, 2) + 13.0 / 12.0 * pow(uavemm - 2 * uavem + uave, 2);
        // double beta1 = 0.25 * pow(uavem - uavep, 2) + 13.0 / 12.0 * pow(uavem - 2 * uave + uavep, 2);
        // double beta2 = 0.25 * pow(3 * uave - 4 * uavep + uavepp, 2) + 13.0 / 12.0 * pow(uave - 2 * uavep + uavepp, 2);

        const double uave_max = std::max(std::max(std::max(uavemm, uavem), std::max(uave, uavep)), uavepp);
        const double uave_min = std::min(std::min(std::min(uavemm, uavem), std::min(uave, uavep)), uavepp);
        // // dimentionless
        uavemm = (uavemm - uave_min) / (uave_max - uave_min + 1e-20);
        uavem = (uavem - uave_min) / (uave_max - uave_min + 1e-20);
        uave = (uave - uave_min) / (uave_max - uave_min + 1e-20);
        uavep = (uavep - uave_min) / (uave_max - uave_min + 1e-20);
        uavepp = (uavepp - uave_min) / (uave_max - uave_min + 1e-20);

        // WENOEXP
        // const double theta = 1;
        // const double beta0 = theta * 1.0 / 2.0 * abs(uavemm - 4 * uavem + 3 * uave) + abs(uavemm - 2 * uavem + uave);
        // const double beta1 = theta * 1.0 / 2.0 * abs(uavem - uavep) + abs(uavem - 2 * uave + uavep);
        // const double beta2 = theta * 1.0 / 2.0 * abs(3 * uave - 4 * uavep + uavepp) + abs(uave - 2 * uavep + uavepp);
        // const double tau = uavemm - 4.0 * uavem + 6.0 * uave - 4.0 * uavep + uavepp;
        // const double epsilon = 1e-16;
        // double w0 = d0 * (1.0 + 2 * pow(tau, 2.0) / (pow(beta0, 2.0) + epsilon));
        // double w1 = d1 * (1.0 + 2 * pow(tau, 2.0) / (pow(beta1, 2.0) + epsilon));
        // double w2 = d2 * (1.0 + 2 * pow(tau, 2.0) / (pow(beta2, 2.0) + epsilon));

        // wenoz
        const double beta0 = 0.25 * pow(uavemm - 4 * uavem + 3 * uave, 2) + 13.0 / 12.0 * pow(uavemm - 2 * uavem + uave, 2);
        const double beta1 = 0.25 * pow(uavem - uavep, 2) + 13.0 / 12.0 * pow(uavem - 2 * uave + uavep, 2);
        const double beta2 = 0.25 * pow(3 * uave - 4 * uavep + uavepp, 2) + 13.0 / 12.0 * pow(uave - 2 * uavep + uavepp, 2);
        const double eps(1.0e-16);
        double tau = fabs(beta0 - beta2);
        double w0 = d0 * (1 + tau / (beta0 + eps));
        double w1 = d1 * (1 + tau / (beta1 + eps));
        double w2 = d2 * (1 + tau / (beta2 + eps));

        const double S = w0 + w1 + w2;

        w0 /= S;
        w1 /= S;
        w2 /= S;

        // 重构
        double ph = w0 * h0 + w1 * h1 + w2 * h2;

        // if (sw == 1)
        // {
        //     switch (gp)
        //     {
        //     case 0:
        //         ph += -1.0 / 60.0 * dif5 * pow(m_deltaX, 5);
        //         break;
        //     case 1:
        //         ph += -385.0 * sqrt(5.0) / 90000.0 * dif5 * pow(m_deltaX, 5);
        //         break;
        //     case 2:
        //         ph += +385.0 * sqrt(5.0) / 90000.0 * dif5 * pow(m_deltaX, 5);
        //         break;
        //     case 3:
        //         ph += +1.0 / 60.0 * dif5 * pow(m_deltaX, 5);
        //         break;
        //     default:
        //         std::cout << "gp error" << std::endl;
        //         std::cin.get();
        //         break;
        //     }
        // }

        return ph;
    }
    else
    {
        int loc_gpNum = 4;
        Array1D<double> loc_gpoints_ref(loc_gpNum), loc_gweights_ref(loc_gpNum);
        sc_math::GaussLobatto_ref(loc_gpNum, loc_gpoints_ref, loc_gweights_ref);

        double gamma = 0;
        if (uavep < uave)
            gamma = 1;
        else
            gamma = -1;

        double Uc = uave;
        double Qmin = min(uavem, uavep);
        double Qmax = max(uavem, uavep) - Qmin;
        double theta = sc_math::Sign(uavep - uave);
        double beta = 1.5;
        double C = (Uc - Qmin + 1e-20) / (Qmax + 1e-20);
        double B = std::exp(theta * beta * (2.0 * C - 1.0));
        double A = (B / std::cosh(beta) - 1.0) / std::tanh(beta);
        double xc = (std::exp(beta) - B) / (B - std::exp(-beta));
        double xcc = std::log(xc) / (2.0 * beta);
        if ((uave - uavem) * (uavep - uave) < 0)
        {
            return uave;
        }
        else
        {
            // std::cout << uavem << " " << uave << " " << uavep << " " << std::endl;
            // std::cout << xc << std::endl;
            return 0.5 * Qmax * (1 + theta * tanh(beta * ((loc_gpoints_ref[gp] + 1) / 2 - xcc))) + Qmin;
            // return uave;
        }
    }
}
inline double thincReconstruction(double uavemmm, double uavemm, double uavem, double uave, double uavep, double uavepp, double uaveppp, int gp, double m_deltaX)
{
    int loc_gpNum = 4;
    Array1D<double> loc_gpoints_ref(loc_gpNum), loc_gweights_ref(loc_gpNum);
    sc_math::GaussLobatto_ref(loc_gpNum, loc_gpoints_ref, loc_gweights_ref);

    double Uc = uave;
    double Qmin = min(uavem, uavep);
    double Qmax = max(uavem, uavep) - Qmin;
    double theta = sc_math::Sign(uavep - uavem);
    double beta = 1.3;
    double C = (Uc - Qmin + 1e-20) / (Qmax + 1e-20);
    double B = std::exp(theta * beta * (2.0 * C - 1.0));
    double A = (B / std::cosh(beta) - 1.0) / std::tanh(beta);
    double xc = (std::exp(beta) - B) / (B - std::exp(-beta));
    double xcc = std::log(xc) / (2.0 * beta);
    if ((uave - uavem) * (uavep - uave) < 0)
    {
        return uave;
    }
    else
    {
        // std::cout << uavem << " " << uave << " " << uavep << " " << std::endl;
        // std::cout << xc << std::endl;
        return 0.5 * Qmax * (1 + theta * tanh(beta * ((loc_gpoints_ref[gp] + 1) / 2 - xcc))) + Qmin;
        // return uave;
    }
}

inline double wenoZQ(double uavemm, double uavem, double uave, double uavep, double uavepp)
{
    double beta1(0);
    beta1 = 1.0 / 144.0 * pow(uavemm - 8 * uavem + 8 * uavep - uavepp, 2);
    beta1 += 1.0 / 15600.0 * pow(-11 * uavemm + 174 * uavem - 326 * uave + 174 * uavep - 11 * uavepp, 2);
    beta1 += 781.0 / 288.0 * pow(-uavemm + 2 * uavem - 2 * uavep + uavepp, 2);
    beta1 += 1421461.0 / 1310400.0 * pow(uavemm - 4 * uavem + 6 * uave - 4 * uavep + uavepp, 2);
    double beta2 = pow(uavem - uave, 2);
    double beta3 = pow(uave - uavep, 2);

    double p1 = (47.0 * uave) / 60.0 - (13.0 * uavem) / 60.0 + (9.0 * uavep) / 20.0 + uavemm / 30.0 - uavepp / 20.0;
    double p2 = (3.0 * uave) / 2.0 - uavem / 2.0;
    double p3 = uave / 2.0 + uavep / 2.0;

    double gamma1(0.98), gamma2(0.01), gamma3(0.01);

    p1 = 1.0 / gamma1 * p1 - gamma2 / gamma1 * p2 - gamma3 / gamma1 * p3;

    double tau = pow((abs(beta1 - beta2) + abs(beta1 - beta3)) / 2.0, 2);

    double eps = 1e-6;
    // 计算权重
    double w1 = gamma1 * (1 + tau / (eps + beta1));
    double w2 = gamma2 * (1 + tau / (eps + beta2));
    double w3 = gamma3 * (1 + tau / (eps + beta3));

    double S = w1 + w2 + w3;

    w1 = w1 / S;
    w2 = w2 / S;
    w3 = w3 / S;

    // 重构
    double h = w1 * p1 + w2 * p2 + w3 * p3;

    std::cout << "WENO-ZQ 有问题" << std::endl;
    std::cin.get();

    h = p2;
    // double h = gamma1 * p1 + gamma2 * p2 + gamma3 * p3;

    return h;
}
#endif