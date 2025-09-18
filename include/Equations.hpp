#ifndef __EQUATIONS__HPP__
#define __EQUATIONS__HPP__

#include "../third_party/scTools/scTools.h"

enum RIEMANNFLUXTYPE
{
	LLF = 0,
	HLLC = 1
};

enum TESTCASETYPE
{
	SMOOTH = 0,
	VORTEX = 1,
	SHARP = 2,
	DAMBREAK = 3,
	DAMBREAK2 = 4,
	DAMBREAK3 = 5,
	RTI = 6,
	STANDING_WAVE = 7,
	ZALESAK = 8,
};

enum BOUNDARYTYPE
{
	PERIOD = 0,
	NEUMANN = 1, // 反射且滑移边界条件
	NoSlip = 2,
	NEUMANN2 = 3, // 根据初值设置边界
	FREESLIP = 4  // 自由滑移边界条件
};

// 定义枚举值对应的字符串数组
inline const std::string RIEMANNFLUXTYPE_STRINGS[] = {
	"LLF", // RIEMANNFLUXTYPE::LLF
	"HLLC" // RIEMANNFLUXTYPE::HLLC
};

inline const std::string TESTCASETYPE_STRINGS[] = {
	"SMOOTH",		 // TESTCASETYPE::SMOOTH
	"VORTEX",		 // TESTCASETYPE::VORTEX
	"SHARP",		 // TESTCASETYPE::SHARP
	"DAMBREAK",		 // TESTCASETYPE::DAMBREAK
	"DAMBREAK2",	 // TESTCASETYPE::DAMBREAK2
	"DAMBREAK3",	 // TESTCASETYPE::DAMBREAK3
	"RTI",			 // TESTCASETYPE::RTI
	"STANDING_WAVE", // TESTCASETYPE::STANDING_WAVE
	"ZALESAK"		 // TESTCASETYPE::ZALESAK
};

class TwoPhaseEquation
{
	// Zhang, F., Cheng, J., & Liu, T. (2023).
	// A Physical-Constraint-Preserving Discontinuous Galerkin Method for Weakly Compressible Two-Phase Flows.
	// Journal of Scientific Computing, 96(3). https://doi.org/10.1007/s10915-023-02306-2

public:
	double xL, xR, yL, yR;
	std::vector<double> outputTimes;
	BOUNDARYTYPE boundaryCondition;
	std::function<double(Array1D<double>)> theVarUh;
	std::function<double(double, double, double)> theVarExact;
	bool u_exact_exist;

	// -------------------- 独有的参数 --------------- //
	double g_rho10, g_rho20;
	double g_sound1, g_sound2;
	double g_u0, g_v0, g_preb;
	double g_epsa;

	std::function<double(double, double)> z1_0;
	std::function<double(double, double)> p1_0;
	std::function<double(double, double)> p2_0;
	std::function<double(double, double)> u_0;
	std::function<double(double, double)> v_0;

	inline int getVarNum()
	{
		return 4;
	}

	// ------------- 方程形式 --------------//
	inline double getVolumeFrac(const Array1D<double> U)
	{
		const double z1rho1 = U[0];
		const double z2rho2 = U[1];

		const double c1_2 = g_sound1 * g_sound1;
		const double c2_2 = g_sound2 * g_sound2;

		const double q = g_rho20 * c2_2 - g_rho10 * c1_2;
		const double qs = z2rho2 * c2_2 - z1rho1 * c1_2;

		const double r = ((q - qs) + sqrt((q - qs) * (q - qs) + 4.0 * z1rho1 * c1_2 * z2rho2 * c2_2)) / (2.0 * z2rho2 * c2_2);
		const double z1 = r / (1.0 + r);

		return z1;
	}

	inline double getPressure(const Array1D<double> U)
	{
		const double z1 = getVolumeFrac(U);
		const double z2 = 1.0 - z1;

		const double z1rho1 = U[0];
		const double z2rho2 = U[1];

		const double c1_2 = g_sound1 * g_sound1;
		const double c2_2 = g_sound2 * g_sound2;

		// const double pre1 = g_preb + c1_2 * (z1rho1 / z1 - g_rho10);
		// const double pre2 = g_preb + c2_2 * (z2rho2 / z2 - g_rho20);

		const double pre1 = c1_2 * (z1rho1 / z1 - g_rho10);
		const double pre2 = c2_2 * (z2rho2 / z2 - g_rho20);

		// note: 过小的分母会导致数值不稳定
		if (z1 < z2)
			return pre2;
		else
			return pre1;
		// return z1 * pre1 + z2 * pre2;
	}

	inline double getSoundSpeed(const Array1D<double> U)
	{
		const double z1 = getVolumeFrac(U);
		const double z2 = 1.0 - z1;

		const double z1rho1 = U[0];
		const double z2rho2 = U[1];

		const double rho = z1rho1 + z2rho2;
		const double rho1 = U[0] / z1;
		const double rho2 = U[1] / z2;

		const double c1_2 = g_sound1 * g_sound1;
		const double c2_2 = g_sound2 * g_sound2;

		const double sc = sqrt(rho * (z1 / (rho1 * c1_2) + z2 / (rho2 * c2_2)));

		// const double sc = sqrt(rho * (z1 * rho2 * c2_2 + z2 * rho1 * c1_2) / (rho1 * c1_2 * rho2 * c2_2));

		return 1.0 / sc;
	}

	// eigen Values
	inline double getMaxEigenValue(const Array1D<double> &U, double nx, double ny)
	{
		const double soundspeed = getSoundSpeed(U);
		const double u = U[2] / (U[0] + U[1]);
		const double v = U[3] / (U[0] + U[1]);

		return fabs(u * nx + v * ny) + soundspeed;
	}

	inline void getPhyFlux(const Array1D<double> Uh, Array1D<double> &Flux, double nx, double ny)
	{
		const double z1rho1 = Uh[0];
		const double z2rho2 = Uh[1];

		const double rho = z1rho1 + z2rho2;

		const double u = Uh[2] / rho;
		const double v = Uh[3] / rho;
		const double pre = getPressure(Uh);

		// Wakimura, H., Takagi, S., & Xiao, F. (2022).
		// Symmetry-preserving enforcement of low-dissipation method based on boundary variation diminishing principle.
		// Computers & Fluids, 233. doi:10.1016/j.compfluid.2021.105227
		if (nx > ny)
		{
			Flux[0] = z1rho1 * u;
			Flux[1] = z2rho2 * u;
			Flux[2] = rho * pow(u, 2) + pre;
			Flux[3] = rho * u * v;
		}
		else
		{
			Flux[0] = z1rho1 * v;
			Flux[1] = z2rho2 * v;
			Flux[2] = rho * v * u;
			Flux[3] = rho * pow(v, 2) + pre;
		}
	}

	// Right eigen matrix based on {z1rho1, z2rho2, rhou, rhov}
	inline void getREigenMatrix(const Array1D<double> &U, double nx, double ny, Array2D<double> &eigMatrix)
	{
		const double z1rho1 = U[0];
		const double z2rho2 = U[1];

		const double rho = z1rho1 + z2rho2;
		const double u = U[2] / rho;
		const double v = U[3] / rho;

		const double z1 = getVolumeFrac(U);
		const double z2 = 1.0 - z1;

		const double c = getSoundSpeed(U);

		const double c1_2 = g_sound1 * g_sound1;
		const double c2_2 = g_sound2 * g_sound2;

		const double zrho1 = z1rho1 * c1_2 / (z1 * z1);
		const double zrho2 = z2rho2 * c2_2 / (z2 * z2);

		const double sc1_2 = c1_2 / z1 - c1_2 / z1 * zrho1 / (zrho1 + zrho2);
		const double sc2_2 = c2_2 / z2 - c2_2 / z2 * zrho2 / (zrho1 + zrho2);

		if (fabs(nx) > fabs(ny)) // x-direction
		{
			eigMatrix[0][0] = z1rho1 / rho;
			eigMatrix[1][0] = z2rho2 / rho;
			eigMatrix[2][0] = u - c;
			eigMatrix[3][0] = v;

			eigMatrix[0][1] = -sc2_2 / (sc1_2 - sc2_2);
			eigMatrix[1][1] = sc1_2 / (sc1_2 - sc2_2);
			eigMatrix[2][1] = u;
			eigMatrix[3][1] = v;

			eigMatrix[0][2] = 0;
			eigMatrix[1][2] = 0;
			eigMatrix[2][2] = 0;
			eigMatrix[3][2] = 1;

			eigMatrix[0][3] = z1rho1 / rho;
			eigMatrix[1][3] = z2rho2 / rho;
			eigMatrix[2][3] = u + c;
			eigMatrix[3][3] = v;
		}
		else // y-direction
		{
			eigMatrix[0][0] = z1rho1 / rho;
			eigMatrix[1][0] = z2rho2 / rho;
			eigMatrix[2][0] = u;
			eigMatrix[3][0] = v - c;

			eigMatrix[0][1] = -sc2_2 / (sc1_2 - sc2_2);
			eigMatrix[1][1] = sc1_2 / (sc1_2 - sc2_2);
			eigMatrix[2][1] = u;
			eigMatrix[3][1] = v;

			eigMatrix[0][2] = 0;
			eigMatrix[1][2] = 0;
			eigMatrix[2][2] = 1;
			eigMatrix[3][2] = 0;

			eigMatrix[0][3] = z1rho1 / rho;
			eigMatrix[1][3] = z2rho2 / rho;
			eigMatrix[2][3] = u;
			eigMatrix[3][3] = v + c;
		}
	}

	// Left eigen matrix based on {z1rho1, z2rho2, rhou, rhov}
	inline void getLEigenMatrix(const Array1D<double> &U, double nx, double ny, Array2D<double> &eigMatrix)
	{
		const double z1rho1 = U[0];
		const double z2rho2 = U[1];

		const double rho = z1rho1 + z2rho2;
		const double u = U[2] / rho;
		const double v = U[3] / rho;

		const double z1 = getVolumeFrac(U);
		const double z2 = 1.0 - z1;

		const double c = getSoundSpeed(U);
		const double c2 = pow(c, 2);

		const double c1_2 = g_sound1 * g_sound1;
		const double c2_2 = g_sound2 * g_sound2;

		const double zrho1 = z1rho1 * c1_2 / (z1 * z1);
		const double zrho2 = z2rho2 * c2_2 / (z2 * z2);

		const double sc1_2 = c1_2 / z1 - c1_2 / z1 * zrho1 / (zrho1 + zrho2);
		const double sc2_2 = c2_2 / z2 - c2_2 / z2 * zrho2 / (zrho1 + zrho2);

		if (fabs(nx) > fabs(ny)) // x-direction
		{
			eigMatrix[0][0] = (sc1_2 + c * u) / (2.0 * c2);
			eigMatrix[0][1] = (sc2_2 + c * u) / (2.0 * c2);
			eigMatrix[0][2] = -1.0 / (2.0 * c);
			eigMatrix[0][3] = 0;

			eigMatrix[1][0] = -z2rho2 * (sc1_2 - sc2_2) / (rho * c2);
			eigMatrix[1][1] = z1rho1 * (sc1_2 - sc2_2) / (rho * c2);
			eigMatrix[1][2] = 0;
			eigMatrix[1][3] = 0;

			eigMatrix[2][0] = -v;
			eigMatrix[2][1] = -v;
			eigMatrix[2][2] = 0;
			eigMatrix[2][3] = 1.0;

			eigMatrix[3][0] = (sc1_2 - c * u) / (2.0 * c2);
			eigMatrix[3][1] = (sc2_2 - c * u) / (2.0 * c2);
			eigMatrix[3][2] = 1.0 / (2.0 * c);
			eigMatrix[3][3] = 0;
		}
		else
		{
			eigMatrix[0][0] = (sc1_2 + c * v) / (2.0 * c2);
			eigMatrix[0][1] = (sc2_2 + c * v) / (2.0 * c2);
			eigMatrix[0][2] = 0;
			eigMatrix[0][3] = -1.0 / (2.0 * c);

			eigMatrix[1][0] = -z2rho2 * (sc1_2 - sc2_2) / (rho * c2);
			eigMatrix[1][1] = z1rho1 * (sc1_2 - sc2_2) / (rho * c2);
			eigMatrix[1][2] = 0;
			eigMatrix[1][3] = 0;

			eigMatrix[2][0] = -u;
			eigMatrix[2][1] = -u;
			eigMatrix[2][2] = 1.0;
			eigMatrix[2][3] = 0;

			eigMatrix[3][0] = (sc1_2 - c * v) / (2.0 * c2);
			eigMatrix[3][1] = (sc2_2 - c * v) / (2.0 * c2);
			eigMatrix[3][2] = 0;
			eigMatrix[3][3] = 1.0 / (2.0 * c);
		}
	}

	inline void getRiemannFlux(const Array1D<double> UL, const Array1D<double> UR, Array1D<double> &Flux, double nx, double ny, int riemannFluxType)
	{
		switch (riemannFluxType)
		{
		case LLF:
			getLLFRiemannFlux(UL, UR, Flux, nx, ny);
			break;
		case HLLC:
			getHLLCRiemannFlux(UL, UR, Flux, nx, ny);
			break;
		default:
			std::cerr << "Riemann flux type is not supported." << std::endl;
			break;
		}

		// debug only
		for (int r = 0; r != getVarNum(); r++)
		{
			if (isnan(Flux[r]))
			{
				std::cerr << "Riemann flux is nan." << std::endl;
				std::cerr << "UL: " << UL << std::endl;
				std::cerr << "UR: " << UR << std::endl;
				std::cerr << "nx: " << nx << std::endl;
				std::cerr << "ny: " << ny << std::endl;
				std::cerr << "riemannFluxType: " << riemannFluxType << std::endl;
				std::cerr << "r: " << r << std::endl;
				std::cerr << "Flux: " << Flux << std::endl;
				std::cin.get();
			}
		}
	}

	// LLF flux
	inline void getLLFRiemannFlux(const Array1D<double> UL, const Array1D<double> UR, Array1D<double> &Flux, double nx, double ny)
	{
		// LLF approximate Riemann flux
		const int varNum = getVarNum();
		const double ws_L = getMaxEigenValue(UL, nx, ny);
		const double ws_R = getMaxEigenValue(UR, nx, ny);
		const double ws = max(ws_L, ws_R);

		Array1D<double> FUL(varNum, 0.0), FUR(varNum, 0.0);
		getPhyFlux(UL, FUL, nx, ny);
		getPhyFlux(UR, FUR, nx, ny);

		// conservative variable
		for (int r = 0; r != varNum; ++r)
			Flux[r] = 0.5 * (FUL[r] + FUR[r] - ws * (UR[r] - UL[r]));
	}

	// HLLC flux
	inline void getHLLCRiemannFlux(const Array1D<double> &UL, const Array1D<double> &UR, Array1D<double> &Flux, double nx, double ny)
	{
		// Zhang, Y., & Zhang, F. (2025).
		// A positivity-preserving HLLC-based discontinuous Galerkin method for weakly compressible two-phase flows.
		// Journal of Computational and Applied Mathematics, 462. https://doi.org/10.1016/j.cam.2024.116467

		const int varNum = getVarNum();

		Array1D<double> FUL(varNum, 0), FUR(varNum, 0);
		Array1D<double> U(varNum, 0);

		// U_L
		const double z1rho1_L = UL[0];
		const double z2rho2_L = UL[1];
		const double rhou_L = UL[2];
		const double rhov_L = UL[3];

		const double rho_L = z1rho1_L + z2rho2_L;
		const double u_L = rhou_L / rho_L;
		const double v_L = rhov_L / rho_L;
		const double pre_L = getPressure(UL);

		const double vn_L = u_L * nx + v_L * ny;
		const double c_L = getSoundSpeed(UL);

		// U_R
		const double z1rho1_R = UR[0];
		const double z2rho2_R = UR[1];
		const double rhou_R = UR[2];
		const double rhov_R = UR[3];

		const double rho_R = z1rho1_R + z2rho2_R;
		const double u_R = rhou_R / rho_R;
		const double v_R = rhov_R / rho_R;
		const double pre_R = getPressure(UR);

		const double vn_R = u_R * nx + v_R * ny;
		const double c_R = getSoundSpeed(UR);

		// Flux_L
		FUL[0] = z1rho1_L * vn_L;
		FUL[1] = z2rho2_L * vn_L;
		FUL[2] = rhou_L * vn_L + pre_L * nx;
		FUL[3] = rhov_L * vn_L + pre_L * ny;

		// Flux_R
		FUR[0] = z1rho1_R * vn_R;
		FUR[1] = z2rho2_R * vn_R;
		FUR[2] = rhou_R * vn_R + pre_R * nx;
		FUR[3] = rhov_R * vn_R + pre_R * ny;

		// Roe-type Characteristic speed
		const double vn_roe = (vn_L * sqrt(rho_L) + vn_R * sqrt(rho_R)) / (sqrt(rho_L) + sqrt(rho_R));
		double c2_roe = (c_L * c_L * sqrt(rho_L) + c_R * c_R * sqrt(rho_R)) / (sqrt(rho_L) + sqrt(rho_R));
		c2_roe += 0.5 * sqrt(rho_L) * sqrt(rho_R) / (sqrt(rho_L) + sqrt(rho_R)) / (sqrt(rho_L) + sqrt(rho_R)) * (vn_R - vn_L) * (vn_R - vn_L);
		const double c_roe = sqrt(c2_roe);

		const double S_L = min(vn_L - c_L, vn_roe - c_roe);
		const double S_R = max(vn_R + c_R, vn_roe + c_roe);
		const double S_star = ((S_L - vn_L) * rho_L * vn_L - (S_R - vn_R) * rho_R * vn_R + pre_R - pre_L) / ((S_L - vn_L) * rho_L - (S_R - vn_R) * rho_R);

		// HLLC flux
		if (S_star >= 0)
		{
			if (S_L >= 0.0)
			{
				for (int r = 0; r != varNum; r++)
					Flux[r] = FUL[r];
			}
			else
			{
				const double k = (S_L - vn_L) / (S_L - S_star);

				U[0] = z1rho1_L;
				U[1] = z2rho2_L;
				U[2] = rho_L * (u_L + (S_star - vn_L) * nx);
				U[3] = rho_L * (v_L + (S_star - vn_L) * ny);

				for (int r = 0; r != varNum; r++)
					Flux[r] = S_L * (k * U[r] - UL[r]) + FUL[r];
			}
		}
		else
		{
			if (S_R <= 0)
			{
				for (int r = 0; r != varNum; r++)
					Flux[r] = FUR[r];
			}
			else
			{
				const double k = (S_R - vn_R) / (S_R - S_star);

				U[0] = z1rho1_R;
				U[1] = z2rho2_R;
				U[2] = rho_R * (u_R + (S_star - vn_R) * nx);
				U[3] = rho_R * (v_R + (S_star - vn_R) * ny);

				for (int r = 0; r != varNum; r++)
					Flux[r] = S_R * (k * U[r] - UR[r]) + FUR[r];
			}
		}
	}

	//------------- 初边值---------------//
	inline void setEquationParameters(TESTCASETYPE type)
	{
		switch (type)
		{
		case SMOOTH:
			xL = 0.0;
			xR = 1.0;
			yL = 0.0;
			yR = 1.0;
			outputTimes = {0, 0.05};
			boundaryCondition = PERIOD;

			g_rho10 = 1.2;
			g_rho20 = 1.0;
			g_sound1 = 1.0;
			g_sound2 = 1.0;
			g_u0 = 1.0;
			g_v0 = 1.0;
			g_preb = 0.0;
			g_epsa = 1.0e-8;

			// z1_0 = [this](double x, double y)
			// { return g_epsa + (1 - 2 * g_epsa) * pow(sin(2.0 * M_PI * (x + y)), 2); };

			//  debug
			z1_0 = [this](double x, double y)
			{ return g_epsa + (1 - 2 * g_epsa) * pow(sin(2.0 * M_PI * (x)), 2); };

			p1_0 = [this](double x, double y)
			{ return g_preb; };
			p2_0 = [this](double x, double y)
			{ return g_preb; };
			u_0 = [this](double x, double y)
			{ return g_u0; };
			v_0 = [this](double x, double y)
			{ return g_v0; };

			this->u_exact_exist = true;
			this->theVarExact = [this](double x, double y, double t)
			{ return z1_0(x - g_u0 * t, y - g_v0 * t); };
			this->theVarUh = [this](Array1D<double> Uh)
			{ return getVolumeFrac(Uh); };
			break;
		case VORTEX:
			xL = -1.0;
			xR = 1.0;
			yL = -1.0;
			yR = 1.0;
			outputTimes = {0, 0.005, 0.01};
			boundaryCondition = PERIOD;

			g_rho10 = 1000.0;
			g_rho20 = 1.0;
			g_sound1 = 30.0;
			g_sound2 = 30.0;
			g_u0 = 2;
			g_v0 = 2;
			g_preb = 0;
			g_epsa = 1.0e-4;

			z1_0 = [this](double x, double y)
			{
				while (x >= 1)
					x -= 2;
				while (x < -1)
					x += 2;
				while (y >= 1)
					y -= 2;
				while (y < -1)
					y += 2;
				const double r = sqrt(x * x + y * y);
				return (1 - 2 * g_epsa) * exp(-pow(r, 2) / pow(0.2, 2)) + g_epsa;
			};
			p1_0 = [this](double x, double y)
			{
				while (x >= 1)
					x -= 2;
				while (x < -1)
					x += 2;
				while (y >= 1)
					y -= 2;
				while (y < -1)
					y += 2;
				const double r = sqrt(x * x + y * y);
				// return 0.1 * g_rho20 * pow(g_sound2, 2) * (1 - exp(-pow(r, 2) / pow(0.2, 2)));
				return 0 * r;
			};
			p2_0 = [this](double x, double y)
			{
				while (x >= 1)
					x -= 2;
				while (x < -1)
					x += 2;
				while (y >= 1)
					y -= 2;
				while (y < -1)
					y += 2;
				const double r = sqrt(x * x + y * y);
				// return 0.1 * g_rho20 * pow(g_sound2, 2) * (1 - exp(-pow(r, 2) / pow(0.2, 2)));
				return 0 * r;
			};
			u_0 = [this](double x, double y)
			{
				while (x >= 1)
					x -= 2;
				while (x < -1)
					x += 2;
				while (y >= 1)
					y -= 2;
				while (y < -1)
					y += 2;

				const double L = 1;
				const double theta = atan2(y, x);
				const double r = sqrt(x * x + y * y);
				// only for debug
				// return g_u0 - 5 * sin(theta) * sqrt(r * exp(-(25 * r) / pow(L, 2))) / L;
				return g_u0 + r * 0.0 + theta * 0.0 + L * 0.0;
			};
			v_0 = [this](double x, double y)
			{
				while (x >= 1)
					x -= 2;
				while (x < -1)
					x += 2;
				while (y >= 1)
					y -= 2;
				while (y < -1)
					y += 2;

				const double L = 1;
				const double theta = atan2(y, x);
				const double r = sqrt(x * x + y * y);
				// only for debug
				// return g_v0 + 5 * cos(theta) * sqrt(r * exp(-(25 * r) / pow(L, 2))) / L;
				return g_v0 + r * 0.0 + theta * 0.0 + L * 0.0;
			};

			this->u_exact_exist = true;
			this->theVarExact = [this](double x, double y, double t)
			{ return z1_0(x - g_u0 * t, y - g_v0 * t); };
			this->theVarUh = [this](Array1D<double> Uh)
			{ return getVolumeFrac(Uh); };
			break;

		case SHARP:
			// Rudman, M. (1997). Volume-Tracking Methods for Interfacial Flow Calculations.
			// International Journal for Numerical Methods in Fluids, 24, 671-691.
			// https://doi.org/10.1002/(sici)1097-0363(19970415)24:7<671::Aid-fld508>3.0.Co;2-9

			xL = 0.0;
			xR = 1.0;
			yL = 0.0;
			yR = 1.0;
			outputTimes = {0, 1.0};
			boundaryCondition = PERIOD;

			g_rho10 = 1000;
			g_rho20 = 1.0;
			g_sound1 = 70.0;
			g_sound2 = 70.0;
			g_u0 = 1.0;
			g_v0 = 1.0;
			g_preb = 0.0;
			g_epsa = 1e-8;

			z1_0 = [this](double x, double y)
			{
				const double r = (xR - xL) / 10.0;
				const double xc = (xL + xR) / 2.0, yc = (yL + yR) / 2.0;

				if (fabs(x - xc) < r && fabs(y - yc) < r)
					// if (pow(x - xc, 2) + pow(y - yc, 2) < pow(r, 2))
					return 1 - g_epsa;

				else
					return g_epsa;
			};
			p1_0 = [this](double x, double y)
			{ return g_preb; };
			p2_0 = [this](double x, double y)
			{ return g_preb; };
			u_0 = [this](double x, double y)
			{ return g_u0; };
			v_0 = [this](double x, double y)
			{ return g_v0; };

			this->u_exact_exist = true;
			this->theVarExact = [this](double x, double y, double t)
			{
				const double r = (xR - xL) / 10;
				const double xc = (xL + xR) / 2, yc = (yL + yR) / 2;
				double z1;
				if (fabs(x - xc) < r && fabs(y - yc) < r)
					z1 = 1 - g_epsa;
				else
					z1 = g_epsa;

				return z1 * g_rho10;
			};
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;

		case ZALESAK:
			xL = 0.0;
			xR = 1.0;
			yL = 0.0;
			yR = 1.0;
			outputTimes = {2};
			boundaryCondition = NEUMANN2;

			g_rho10 = 1000;
			g_rho20 = 1.0;
			g_sound1 = 70.0;
			g_sound2 = 70.0;
			g_u0 = 1.0;
			g_v0 = 1.0;
			g_preb = 0.0;
			g_epsa = 1e-8;

			z1_0 = [this](double x, double y)
			{
				const double x0 = 0.5;
				const double y0 = 0.5;
				const double r = 0.15;

				// if (pow(x - x0, 2) + pow(y - y0, 2) < pow(r, 2))
				// 	if ((x < x0 - 0.025) || (x > x0 + 0.025) || (y > y0 + 0.1))
				// 		return 1 - g_epsa;
				// 	else
				// 		return g_epsa;
				// else
				// 	return g_epsa;

				if (pow(x - x0, 2) + pow(y - y0, 2) < pow(r, 2))
					return 1 - g_epsa;
				else
					return g_epsa;
			};
			p1_0 = [this](double x, double y)
			{ return g_preb; };
			p2_0 = [this](double x, double y)
			{ return g_preb; };
			u_0 = [this](double x, double y)
			{
				const double omega = M_PI;
				return omega * (0.5 - y);
			};
			v_0 = [this](double x, double y)
			{
				const double omega = M_PI;
				return omega * (x - 0.5);
			};

			this->u_exact_exist = false;
			this->theVarExact = [this](double x, double y, double t)
			{ return 0; };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;

		case STANDING_WAVE:
			xL = 0.0;
			xR = 0.5;
			yL = 0.0;
			yR = 0.5;
			outputTimes = {3.0};
			boundaryCondition = NEUMANN;

			g_rho10 = 1000.0;
			g_rho20 = 1.0;
			g_sound1 = 100;
			g_sound2 = 100;
			g_u0 = 0.0;
			g_v0 = 0.0;
			g_preb = 1e4;
			g_epsa = 1e-6;

			z1_0 = [this](double x, double y)
			{
				double hs = 0.25, A = 0.025;
				double kn = 2 * M_PI;

				if (y <= hs + A * cos(kn * x))
					return 1.0 - g_epsa;
				else
					return g_epsa;
			};
			p1_0 = [this](double x, double y)
			{ return g_preb; };
			p2_0 = [this](double x, double y)
			{ return g_preb; };
			u_0 = [this](double x, double y)
			{ return g_u0; };
			v_0 = [this](double x, double y)
			{ return g_v0; };

			this->u_exact_exist = false;
			this->theVarExact = [this](double x, double y, double t)
			{ return 0; };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;

		case RTI:
			// An, R., Gu, Z., Zhou, T., &Yu, C.(2023).
			// Numerical simulation of incompressible interfacial flows by a level set redistancing method with improved mass conservation.
			// Ocean Engineering, 290. https : // doi.org/10.1016/j.oceaneng.2023.116428

			xL = 0.0;
			xR = 1.0;
			yL = 0.0;
			yR = 4.0;
			outputTimes = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
			boundaryCondition = NEUMANN;

			g_rho10 = 1.225;
			g_rho20 = 0.1694;

			g_sound1 = 100.0;
			g_sound2 = 100.0;

			g_u0 = 0.0;
			g_v0 = 0.0;
			g_preb = 10000.0;

			g_epsa = 1.0e-8;
			z1_0 = [this](double x, double y)
			{
				const double yi = 2.0 + 0.05 * cos(2.0 * M_PI * x);
				if (y >= yi)
					return 1.0 - g_epsa;
				else
					return g_epsa;
			};
			p1_0 = [this](double x, double y)
			{ return g_preb; };
			p2_0 = [this](double x, double y)
			{ return g_preb; };
			u_0 = [this](double x, double y)
			{ return g_u0; };
			v_0 = [this](double x, double y)
			{ return g_v0; };
			break;

		case DAMBREAK:
			// Li, Z., Oger, G., & Le Touzé, D. (2020).
			// A finite volume WENO scheme for immiscible inviscid two-phase flows.
			// Journal of Computational Physics, 418. https://doi.org/10.1016/j.jcp.2020.109601
			xL = 0.0;
			xR = 1.61;
			yL = 0.0;
			yR = 0.805;
			outputTimes = {0.0, 0.1599, 0.2766, 0.3733, 0.4499, 0.5733, 0.8623, 1.0233, 1.1666};
			boundaryCondition = NEUMANN;

			// 独有变量
			g_rho10 = 1000;
			g_rho20 = 1.0;
			g_sound1 = 100.0;
			g_sound2 = 100.0;
			// g_sound1 = 300.0;
			// g_sound2 = 300.0;
			g_u0 = 0.0;
			g_v0 = 0.0;
			g_preb = 1e4;
			g_epsa = 1e-6;

			z1_0 = [this](double x, double y)
			{
				const double H0 = 0.3;
				const double L0 = 0.6;
				if (x >= xR - L0 && x <= xR && y >= yL && y <= yL + H0)
					return 1.0 - g_epsa;
				else
					return g_epsa;
			};
			p1_0 = [this](double x, double y)
			{ return g_preb; };
			p2_0 = [this](double x, double y)
			{ return g_preb; };
			u_0 = [this](double x, double y)
			{ return g_u0; };
			v_0 = [this](double x, double y)
			{ return g_v0; };

			this->u_exact_exist = false;
			this->theVarExact = [this](double x, double y, double t)
			{ return 0; };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;

		case DAMBREAK2:
			// Duy, T.-N., Nguyen, V.-T., Phan, T.-H., Kim, D.-H., & Park, W.-G. (2022).
			// A free surface flow solver based on an efficient improvement to a coupling method for interface computations.
			// Computers & Mathematics with Applications, 124, 21-41. https://doi.org/10.1016/j.camwa.2022.08.020

			xL = 0.0;
			xR = 1.61;
			yL = 0.0;
			yR = 0.6;
			outputTimes = {1.5};
			boundaryCondition = NEUMANN;

			// 独有变量
			g_rho10 = 1000;
			g_rho20 = 1.0;
			g_sound1 = 100.0;
			g_sound2 = 100.0;
			g_u0 = 0.0;
			g_v0 = 0.0;
			g_preb = 1e4;
			g_epsa = 1e-5;

			z1_0 = [this](double x, double y)
			{
				const double H0 = 0.6;
				const double L0 = 0.6;
				if (x >= xR - L0 && x <= xR && y >= yL && y <= yL + H0)
					return 1.0 - g_epsa;
				else
					return g_epsa;
			};
			p1_0 = [this](double x, double y)
			{ return g_preb; };
			p2_0 = [this](double x, double y)
			{ return g_preb; };
			u_0 = [this](double x, double y)
			{ return g_u0; };
			v_0 = [this](double x, double y)
			{ return g_v0; };

			this->u_exact_exist = false;
			this->theVarExact = [this](double x, double y, double t)
			{ return 0; };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;

		case DAMBREAK3:
			// Duy, T.-N., Nguyen, V.-T., Phan, T.-H., Kim, D.-H., & Park, W.-G. (2022).
			// A free surface flow solver based on an efficient improvement to a coupling method for interface computations.
			// Computers & Mathematics with Applications, 124, 21-41. https://doi.org/10.1016/j.camwa.2022.08.020

			xL = 0.0;
			xR = 0.8;
			yL = 0.0;
			yR = 0.6;
			outputTimes = {1.5};
			boundaryCondition = NEUMANN;

			// 独有变量
			g_rho10 = 1000;
			g_rho20 = 1.0;
			g_sound1 = 100.0;
			g_sound2 = 100.0;
			g_u0 = 0.0;
			g_v0 = 0.0;
			g_preb = 1e4;
			g_epsa = 1e-5;

			z1_0 = [this](double x, double y)
			{
				const double H0 = 0.4;
				const double L0 = 0.2;
				if (x >= xL && x <= xL + L0 && y >= yL && y <= yL + H0)
					return 1.0 - g_epsa;
				else
					return g_epsa;
			};
			p1_0 = [this](double x, double y)
			{ return g_preb; };
			p2_0 = [this](double x, double y)
			{ return g_preb; };
			u_0 = [this](double x, double y)
			{ return g_u0; };
			v_0 = [this](double x, double y)
			{ return g_v0; };

			this->u_exact_exist = false;
			this->theVarExact = [this](double x, double y, double t)
			{ return 0; };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;

		default:
			throw std::invalid_argument("Invalid example number");
		}
	}

	inline void getU0(double xP, double yP, Array1D<double> &U)
	{
		const double z_1 = z1_0(xP, yP);
		const double z_2 = 1.0 - z_1;

		const double p_1 = p1_0(xP, yP);
		const double p_2 = p2_0(xP, yP);

		const double rho10 = (p_1 - g_preb) / pow(g_sound1, 2) + g_rho10;
		const double rho20 = (p_2 - g_preb) / pow(g_sound2, 2) + g_rho20;

		const double u = u_0(xP, yP);
		const double v = v_0(xP, yP);

		const double z1rho1 = rho10 * z_1;
		const double z2rho2 = rho20 * z_2;

		const double rho = z1rho1 + z2rho2;

		U[0] = z1rho1;
		U[1] = z2rho2;
		U[2] = rho * u;
		U[3] = rho * v;
	}

	// -------------- 输出的值 ---------------//
	inline void getVitalVarVal(const Array1D<double> Uh, Array1D<double> &VitalVar)
	{
		const double z1 = getVolumeFrac(Uh);
		const double z2 = 1.0 - z1;

		// get primitive variables
		const double rho1 = Uh[0] / z1;
		const double rho2 = Uh[1] / z2;
		const double rho = (Uh[0] + Uh[1]);
		const double u = Uh[2] / (Uh[0] + Uh[1]);
		const double v = Uh[3] / (Uh[0] + Uh[1]);
		const double pre = getPressure(Uh);

		VitalVar[0] = Uh[0];
		VitalVar[1] = Uh[1];
		VitalVar[2] = z1;
		VitalVar[3] = rho1;
		VitalVar[4] = rho2;
		VitalVar[5] = rho;
		VitalVar[6] = u;
		VitalVar[7] = v;
		VitalVar[8] = pre;
	}
	inline int getVitalVarNum()
	{
		return 9;
	}
	inline void getVitalVarName(Array1D<std::string> &VitalVarName)
	{
		VitalVarName[0] = "RHO1Z1";
		VitalVarName[1] = "RHO2Z2";
		VitalVarName[2] = "Z1";
		VitalVarName[3] = "RHO1";
		VitalVarName[4] = "RHO2";
		VitalVarName[5] = "RHO";
		VitalVarName[6] = "U";
		VitalVarName[7] = "V";
		VitalVarName[8] = "PRE";
	}
};

#endif