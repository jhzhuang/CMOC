#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "common.h"

namespace CMOC {

	/*  ponit type enumerator  */
	typedef enum _Point_Type {
		Inner = 0x00,
		Boundary = 0x10,
		Axis = 0x11,
		Shock = 0x20
	} Point_Type;

	/*  point  */
	typedef struct _Point
	{
		Point_Type type;
		Real x;
		Real y;
		Real angle;
	} Point;

	/*  line  */
	typedef struct _Line {
		Point point;
		Real dx;
		Real dy;
	} Line;
}

#endif // !_GEOMETRY_H_