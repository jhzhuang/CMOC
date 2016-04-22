#ifndef _AIR_H_
#define _AIR_H_

#include "common.h"
#include <math.h>

#define FIXED_CP 0
#define POLYNOMIAL_CP 1

#ifndef CP_MODEL
#define CP_MODEL FIXED_CP
#endif // !CP_MODEL

namespace CMOC {

	/*  Air  */
	typedef struct _Air {
		Real pressure;
		Real density;
		Real temperature;
	} Air;

	/*  gas constant  */
	extern Real R_Air;

	/*  initial an air object by pressure and density  */
	void Init_Air_P_Rho(Air *air, Real p, Real rho);

	/*  get heat ratio at constant pressure  */
	Real Get_Cp(const Air &air);

	/*  get heat ratio at constant volume  */
	Real Get_Cv(const Air &air);

	/*  get specific heat ratio  */
	Real Get_Gamma(const Air &air);

	/*  get sound speed  */
	Real Get_Sound_Speed(const Air &air);

#if CP_MODEL == FIXED_CP

	/*  set specific heat ratio  */
	void Set_Gamma(Real gamma);

#elif CP_MODEL == POLYNOMIAL_CP



#else

#define UNDEFINED_CP_MODEL

#endif

	inline void Init_Air_P_Rho(Air *air, Real p, Real rho) {
		air->pressure = p;
		air->density = rho;
		air->temperature = p / (rho * R_Air);
	}

	#if CP_MODEL == FIXED_CP

	/*  specific heat ratio  */
	extern Real Gamma;

	inline void Set_Gamma(Real gamma) {
		Gamma = gamma;
	}

	inline Real Get_Cp(const Air &air) {
		return Gamma * R_Air / (Gamma - 1);
	}

	inline Real Get_Cv(const Air &air) {
		return R_Air / (Gamma - 1);
	}

	inline Real Get_Gamma(const Air &air) {
		return Gamma;
	}

	inline Real Get_Sound_Speed(const Air &air) {
		return sqrt(Gamma * air.pressure / air.density);
	}

#elif CP_MODEL == POLYNOMIAL_CP



#endif


}

#endif /*    */
