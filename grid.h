#ifndef _GRID_H_
#define _GRID_H_

#include "air.h"
#include "geometry.h"

namespace CMOC {

	/*  node  */
	typedef struct _Node {
		Point point;
		Air air;
		Real velocity;
		Real angle;
		int id;
	} Node;

	/*  cell type  */
	typedef enum _Cell_Type {
		Triangle = 0,
		Quadrilateral = 1
	} Cell_Type;

	/*  cell  */
	typedef struct _Cell {
		Cell_Type type;
		union {
			/*  triangle cell  */
			struct {
				Node *t_node1;
				Node *t_node2;
				Node *t_node3;
			};
			/*  quadrilateral cell  */
			struct {
				Node *q_node1;
				Node *q_node2;
				Node *q_node3;
				Node *q_node4;
			};
		};
	} Cell;

}

#endif // !_GRID_H_
