#include "air.h"
#include <math.h>

namespace CMOC {

	Real Get_Sound_Speed(const Air &air) {
		return sqrt(Get_Gamma(air) * R_Air * air.temperature);
	}

}