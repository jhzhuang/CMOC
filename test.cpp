#include "air.h"
#include "geometry.h"
#include "grid.h"
#include "solver.h"

int main(int argc, char **argv) {
	CMOC::Node node1, node2, node3, node4;
	CMOC::Cell cell;

	node2.air.pressure = 12474.9;
	node2.air.density = 0.1;
	node2.point.x = 0.866;
	node2.point.y = 0.5;
	node2.velocity = 733.38;
	node2.angle = 0;

	node4.air.pressure = 3409.9;
	node4.air.density = 0.0396;
	node4.point.x = 0.9357;
	node4.point.y = 0.3527;
	node4.velocity = 984.38;
	node4.angle = - 20.0 * 3.1415926535897932 / 180.0;
	
	cell.q_node1 = &node1;
	cell.q_node2 = &node2;
	cell.q_node3 = &node3;
	cell.q_node4 = &node4;

	Solve(&cell);

	return 0;
}

