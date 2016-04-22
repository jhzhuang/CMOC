#include "solver.h"

namespace CMOC {

	Real Qp(Real M, Real rho, Real v);

	Real Sx_plus(Real M, Real y, Real theta, Real alpha);

	Real Sx_minus(Real M, Real y, Real theta, Real alpha);

	Real T_plus(Real s_plus, Real q_plus, Real dx, Real p, Real theta);

	Real T_minus(Real s_minus, Real q_minus, Real dx, Real p, Real theta);

	void Intersect(Point &ipoint, Line &l1, Line &l2);

	static inline Real Qp(Real M, Real rho, Real v) {
		return sqrt(M * M - 1.0) / (rho * v * v);
	}

	static inline Real Sx_plus(Real M, Real y, Real theta, Real alpha) {
		return sin(theta) / (y * M * cos(theta + alpha));
	}

	static inline Real Sx_minus(Real M, Real y, Real theta, Real alpha) {
		return sin(theta) / (y * M * cos(theta + alpha));
	}

	static inline Real T_plus(Real s_plus, Real q_plus, Real dx, Real p, Real theta) {
		return -s_plus * dx + q_plus * p + theta;
	}

	static inline Real T_minus(Real s_minus, Real q_minus, Real dx, Real p, Real theta) {
		return -s_minus * dx + q_minus * p - theta;
	}

	static inline void Intersect(Point &ipoint, Line &l1, Line &l2) {
		ipoint.x = (l1.dx * l2.dx * (l2.point.y - l1.point.y) + l1.dx * l2.dy * l2.point.x -
			l2.dx * l1.dy * l1.point.x) / (l1.dx * l2.dy - l2.dx * l1.dy);
		ipoint.y = (l1.dy * l2.dy * (l2.point.x - l1.point.x) + l1.dy * l2.dx * l2.point.y -
			l2.dy * l1.dx * l1.point.y) / (l1.dy * l2.dx - l2.dy * l1.dx);
	}

	class Inter {
	public:
		Inter(Cell *cell):
			node2(*(cell->q_node2)), node3(*(cell->q_node3)), node4(*(cell->q_node4)) {
			Initial();
		}
	public:
		void Solve(void) {
			do
			{
				Charater_Coeff();
				Charater_Solve();
				Compatible_Coeff();
				Compatible_Solve();
			} while (err > 1.0e-4);
		}
	private:
		void Initial(void);
		void Charater_Coeff(void);
		void Charater_Solve(void);
		void Compatible_Coeff(void);
		void Compatible_Solve(void);
	private:
		Node &node2;
		Node &node3;
		Node &node4;
		Node node5;
	private:
		Real err;
		Real M2, M3, M4;
		Real alpha2, alpha3, alpha4;
		Real dx_o, dx_plus, dx_minus, dx_12;
		Real dy_o, dy_plus, dy_minus, dy_12;
		Real q_plus, q_minus;
		Real t_plus, t_minus;
		Real t_o1, t_o2;
		Real r0, A0;
	};

	void Inter::Initial(void) {
		Real w2, w4;
		Real s_plus, s_minus;
		Real a5;

		M2 = node2.velocity / Get_Sound_Speed(node2.air);
		M4 = node4.velocity / Get_Sound_Speed(node4.air);

		alpha2 = asin(1.0 / M2);
		alpha4 = asin(1.0 / M4);

		dx_o = 0.5 * (cos(node2.angle) + cos(node4.angle));
		dy_o = 0.5 * (sin(node2.angle) + sin(node4.angle));
		dx_plus = cos(node4.angle + alpha4);
		dy_plus = sin(node4.angle + alpha4);
		dx_minus = cos(node2.angle - alpha2);
		dy_minus = sin(node2.angle - alpha2);

		Charater_Solve();

		q_plus = Qp(M4, node4.air.density, node4.velocity);
		q_minus = Qp(M2, node2.air.density, node2.velocity);
		s_plus = Sx_plus(M4, node4.point.y, node4.angle, alpha4);
		s_minus = Sx_plus(M2, node2.point.y, node2.angle, alpha2);
		t_plus = T_plus(s_plus, q_plus, node3.point.x - node4.point.x, node4.air.pressure, node4.angle);
		t_minus = T_minus(s_minus, q_minus, node3.point.x - node2.point.x, node2.air.pressure, node2.angle);

		w2 = sqrt((node4.point.x - node5.point.x) * (node4.point.x - node5.point.x) +
			(node4.point.y - node5.point.y) * (node4.point.y - node5.point.y)) /
			sqrt((node4.point.x - node2.point.x) * (node4.point.x - node2.point.x) +
				(node4.point.y - node2.point.y) * (node4.point.y - node2.point.y));
		w4 = 1.0 - w2;

		node5.velocity = node2.velocity * w2 + node4.velocity * w4;
		node5.air.pressure = node2.air.pressure * w2 + node4.air.pressure * w4;
		node5.air.density = node2.air.density * w2 + node4.air.density * w4;
		a5 = Get_Sound_Speed(node5.air);

		r0 = node5.air.density * node5.velocity;
		A0 = a5 * a5;
		t_o1 = r0 * node5.velocity + node5.air.pressure;
		t_o2 = node5.air.pressure - A0 * node5.air.density;

		Compatible_Solve();
	}

	void Inter::Charater_Coeff(void) {
		dx_o = 0.5 * (cos(node5.angle) + cos(node3.angle));
		dy_o = 0.5 * (sin(node5.angle) + sin(node3.angle));
		dx_plus = 0.5 * (cos(node4.angle + alpha4) + cos(node3.angle + alpha3));
		dy_plus = 0.5 * (sin(node4.angle + alpha4) + sin(node3.angle + alpha3));
		dx_minus = 0.5 * (cos(node2.angle - alpha2) + cos(node3.angle - alpha3));
		dy_minus = 0.5 * (sin(node2.angle - alpha2) + sin(node3.angle - alpha3));
	}

	void Inter::Charater_Solve(void) {
		Line c_plus, c_minus, c_o, c_12;

		c_plus.point = node4.point;
		c_plus.dx = dx_plus;
		c_plus.dy = dy_plus;

		c_minus.point = node2.point;
		c_minus.dx = dx_minus;
		c_minus.dy = dy_minus;

		Intersect(node3.point, c_plus, c_minus);

		c_o.point = node3.point;
		c_o.dx = dx_o;
		c_o.dy = dy_o;

		c_12.point = node2.point;
		c_12.dx = node4.point.x - node2.point.x;
		c_12.dy = node4.point.y - node2.point.y;

		Intersect(node5.point, c_o, c_12);
	}

	void Inter::Compatible_Coeff(void) {
		Real w2, w4;
		Real s_plus, s_minus;
		Real a3, a5;

		a3 = Get_Sound_Speed(node3.air);
		M3 = node3.velocity / a3;
		alpha3 = asin(1.0 / M3);

		q_plus = 0.5 * (Qp(M4, node4.air.density, node4.velocity) +
			Qp(M3, node3.air.density, node3.velocity));
		q_minus = 0.5 * (Qp(M2, node2.air.density, node2.velocity) +
			Qp(M3, node3.air.density, node3.velocity));
		s_plus = 0.5 * (Sx_plus(M4, node4.point.y, node4.angle, alpha4) + 
			Sx_plus(M3, node3.point.y, node3.angle, alpha3));
		s_minus = 0.5 * (Sx_plus(M2, node2.point.y, node2.angle, alpha2) +
			Sx_plus(M3, node3.point.y, node3.angle, alpha3));
		t_plus = T_plus(s_plus, q_plus, node3.point.x - node4.point.x, node4.air.pressure, node4.angle);
		t_minus = T_minus(s_minus, q_minus, node3.point.x - node2.point.x, node2.air.pressure, node2.angle);

		w2 = sqrt((node4.point.x - node5.point.x) * (node4.point.x - node5.point.x) +
			(node4.point.y - node5.point.y) * (node4.point.y - node5.point.y)) /
			sqrt((node4.point.x - node2.point.x) * (node4.point.x - node2.point.x) +
				(node4.point.y - node2.point.y) * (node4.point.y - node2.point.y));
		w4 = 1.0 - w2;

		node5.velocity = node2.velocity * w2 + node4.velocity * w4;
		node5.air.pressure = node2.air.pressure * w2 + node4.air.pressure * w4;
		node5.air.density = node2.air.density * w2 + node4.air.density * w4;
		a5 = Get_Sound_Speed(node5.air);

		r0 = 0.5 * (node5.air.density * node5.velocity + node3.air.density * node3.velocity);
		A0 = 0.5 * (a5 * a5 + a3 * a3);
		t_o1 = r0 * node5.velocity + node5.air.pressure;
		t_o2 = node5.air.pressure - A0 * node5.air.density;
	}

	void Inter::Compatible_Solve(void) {
		Real p3, rho3, theta3, v3;

		p3 = (t_plus + t_minus) / (q_plus + q_minus);
		theta3 = t_plus - q_plus * p3;
		rho3 = (p3 - t_o2) / A0;
		v3 = (t_o1 - p3) / r0;

		err = fabs(node3.air.pressure - p3) + fabs(node3.air.density - rho3) +
			fabs(node3.angle - theta3) + fabs(node3.velocity - v3);

		node3.air.pressure = p3;
		node3.air.density = rho3;
		node3.angle = theta3;
		node3.velocity = v3;
	}

	void Solve(Cell *cell) {
		Inter solver(cell);
		solver.Solve();
	}
}