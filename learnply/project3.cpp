#define _USE_MATH_DEFINES

#include "glError.h"
#include "project3.h"
#include "polyhedron.h"
#include "iostream"
#include "GL/freeglut.h"
#include <vector>
#include <set>
#include <map>
#include <math.h>
#include "polyline2.h"

#define EPSILON 0.001
#define MIN_K 0.05
#define STEP 0.001

extern Polyhedron* poly;
extern std::vector<Polyline2> polylines;
extern std::list<Singularity> singularities;

void streamline(Polyline2& polyline, icVector3 seed, const double step) {
	streamlineFB(polyline, seed, step);
	Polyline2 lineBack;
	streamlineFB(lineBack, seed, step, false);
	polyline.merge(lineBack);
}

bool onBoundary(icVector3 nextPosition, icVector3 min, icVector3 max) {
	if (nextPosition.x == min.x || nextPosition.x == max.x) {
		return true;
	}
	else if (nextPosition.y == min.y || nextPosition.y == max.y) {
		return true;
	}
	else {
		return false;
	}
}

void streamlineFB(Polyline2& polyline, icVector3 seed, const double step, bool forward) {
	polyline.vertices.push_back(seed);
	Quad* quad = findQuad(seed);
	icVector3 min, max;
	findMinMaxField(min, max);
	icVector3 currentPosition = seed;
	double coef = 1.0;
	if (!forward) {
		coef = -1.0;
	}
	while (quad != nullptr) {
		icVector3 currentVector = getVector(quad, currentPosition);
		if (currentVector.length() < EPSILON) {
			break;
		}
		icVector3 nextPosition = currentPosition + step * currentVector * coef;
		if (onBoundary(nextPosition, min, max)) {
			polyline.vertices.push_back(nextPosition);
			break;
		}
		Quad* nextQuad = nullptr;
		//streamlineTrace(nextQuad, quad, currentPosition, currentVector * coef, step, min, max);
		nextQuad = findQuad(nextPosition);
		quad = nextQuad;
		currentPosition = nextPosition;
		polyline.vertices.push_back(currentPosition);
	}
}

// maybe done?
void findMinMaxField(icVector3& min, icVector3& max) {
	min.x = poly->vlist[0]->x;
	min.y = poly->vlist[0]->y;
	min.z = poly->vlist[0]->z;

	max = min;

	for (int i = 0; i < poly->nverts; i++) {
		if (min.x > poly->vlist[i]->x) {
			min.x = poly->vlist[i]->x;
		}
		if (min.y > poly->vlist[i]->y) {
			min.y = poly->vlist[i]->y;
		}
		if (min.z > poly->vlist[i]->z) {
			min.z = poly->vlist[i]->z;
		}

		if (max.x < poly->vlist[i]->x) {
			max.x = poly->vlist[i]->x;
		}
		if (max.y < poly->vlist[i]->y) {
			max.y = poly->vlist[i]->y;
		}
		if (max.z < poly->vlist[i]->z) {
			max.z = poly->vlist[i]->z;
		}
	}
}

bool insideQuad(const Quad* quad, const icVector3 p) {
	double v0x = quad->verts[2]->x;
	double v0y = quad->verts[2]->y;
	double v2x = quad->verts[0]->x;
	double v2y = quad->verts[0]->y;
	if (p.x >= v0x && p.x <= v2x && p.y >= v0y && p.y <= v2y) {
		return true;
	}
	else {
		return false;
	}
}

Quad* findQuad(const icVector3 p) {
	for (int i = 0; i < poly->nquads; i++) {
		Quad* tempQuad = poly->qlist[i];
		if (insideQuad(tempQuad, p)) {
			return tempQuad;
		}
	}
	return nullptr;
}

bool isZero(double x) {
	double e = std::numeric_limits<double>::epsilon();
	return std::abs(x) < e;
}

//not done
//icVector3 bilinear(icVector3 pt, Vertex v0, Vertex v1, Vertex v2, Vertex v3) {
//	double alpha0 = (pt.x - v0.x) / (v1.x - v0.x);
//	double alpha1 = (pt.x - v3.x) / (v2.x - v3.x);
//	double py = (1 - alpha0) * v0.y + (alpha0) * v1.y;
//	double qy = (1 - alpha1) * v3.y + (alpha0) * v2.y;
//	double beta = (pt.y - py) / (qy - py);
//
//	double vx = (1-alpha0) * (1-beta) * v0.vx + alpha0 * (1-beta) * v1.vx + alpha1 * beta * v2.vx + (1-alpha1) * beta * v3.vx;
//	double vy = (1-alpha0) * (1-beta) * v0.vy + alpha0 * (1-beta) * v1.vy + alpha1 * beta * v2.vy + (1-alpha1) * beta * v3.vx;
//
//	return icVector3(vx, vy, 0);
//}

// maybe done
icVector3 getVector(Quad* quad, const icVector3 p) {
	Vertex* x2y2 = quad->verts[0];
	Vertex* x1y2 = quad->verts[1];
	Vertex* x1y1 = quad->verts[2];
	Vertex* x2y1 = quad->verts[3];

	double x1 = x1y2->x;
	double x2 = x2y2->x;
	double y1 = x1y1->y;
	double y2 = x1y2->y;

	double x1y1XVector = x1y1->vx;
	double x1y1YVector = x1y1->vy;
	double x1y2XVector = x1y2->vx;
	double x1y2YVector = x1y2->vy;
	double x2y1XVector = x2y1->vx;
	double x2y1YVector = x2y1->vy;
	double x2y2XVector = x2y2->vx;
	double x2y2YVector = x2y2->vy;

	double newXVector = (((x2 - p.x) / (x2 - x1)) * ((y2 - p.y) / (y2 - y1)) * x1y1XVector) +
		(((p.x - x1) / (x2 - x1)) * ((y2 - p.y) / (y2 - y1)) * x2y1XVector) +
		(((x2 - p.x) / (x2 - x1)) * ((p.y - y1) / (y2 - y1)) * x1y2XVector) +
		(((p.x - x1) / (x2 - x1)) * ((p.y - y1) / (y2 - y1)) * x2y2XVector);

	double newYVector = (((x2 - p.x) / (x2 - x1)) * ((y2 - p.y) / (y2 - y1)) * x1y1YVector) +
		(((p.y - x1) / (x2 - x1)) * ((y2 - p.y) / (y2 - y1)) * x2y1YVector) +
		(((x2 - p.y) / (x2 - x1)) * ((p.y - y1) / (y2 - y1)) * x1y2YVector) +
		(((p.y - x1) / (x2 - x1)) * ((p.y - y1) / (y2 - y1)) * x2y2YVector);

	return icVector3(newXVector, newYVector, 0.0);
}

/*
void streamlineTrace(Quad*& nextQuad, Quad*& currentQuad, icVector3 currentPos, icVector3 currentVec, double t, const icVector3 min, const icVector3 max) {

	bool insideQuad = false;
	while (!insideQuad) {
		if ((currentPos.x < min.x || currentPos.x > max.x) || (currentPos.y < min.y || currentPos.y > max.y)) {
			nextQuad = nullptr; 
			return;
		}

		double t_ = INFINITY;
		Quad* nextQuad_ = nullptr;
		for (int i = 0; i < 4; i++) {
			Edge* edge = currentQuad->edges[i];
			Vertex* v0 = edge->verts[0];
			Vertex* v1 = edge->verts[1];
			double temp;
			if (std::abs(v0->x - v1->x) < EPSILON) {
				temp = (v0->x - currentPos.x) / currentVec.x;
			}
			else {
				temp = (v0->y - currentPos.y) / currentVec.y;
			}
			if (temp > 0 && temp < t_) {
				t_ = temp;
				if (edge->quads[0] != currentQuad && edge->quads[1] == currentQuad) {
					nextQuad_ = edge->quads[0];
				}
				else if (edge->quads[0] == currentQuad && edge->quads[1] != currentQuad) {
					nextQuad_ = edge->quads[1];
				}
			}
		}
		if (nextQuad_ == nullptr) {
			currentPos = currentPos + currentVec * t;
			nextQuad = findQuad(currentPos);
			return;
		}
		else {
			if (t_ >= t) {
				insideQuad = true;
			} 
			else {
				currentQuad = nextQuad_;
				t = t - t_;
				currentPos = currentPos + currentVec * t_;
			}
		}
	}
	nextQuad = currentQuad;
}
*/

bool quadraticRoot(double& r0, double& r1, double a, double b, double c) {
	double m = b * b - 4 * a * c;
	if (m < 0) {
		return false;
	}
	r0 = (-b - std::sqrt(m)) / (2 * a);
	r1 = (-b + std::sqrt(m)) / (2 * a);
	return true;
}

//bool singRoot(double r0, double r1, double a, double b, double c, double d) {
//	double f0 = b - a - (c + d);
//	double f1 = (c + d);
//	double f2 = a;
//	return quadraticRoot(r0, r1, f0, f1, f2);
//}

void extractSingularity() {
	singularities.clear();
	for (int i = 0; i < poly->nquads; i++) {
		icVector3 vecx1y1 = poly->qlist[i]->verts[2]->vec();
		icVector3 posx1y1 = poly->qlist[i]->verts[2]->pos();
		icVector3 vecx2y1 = poly->qlist[i]->verts[3]->vec();
		icVector3 posx2y1 = poly->qlist[i]->verts[3]->pos();
		icVector3 vecx2y2 = poly->qlist[i]->verts[0]->vec();
		icVector3 posx2y2 = poly->qlist[i]->verts[0]->pos();
		icVector3 vecx1y2 = poly->qlist[i]->verts[1]->vec();
		icVector3 posx1y2 = poly->qlist[i]->verts[1]->pos();

		icVector3 pt(0.);
		double f_11 = vecx1y1.x;
		double f_12 = vecx1y2.x;
		double f_21 = vecx2y1.x;
		double f_22 = vecx2y2.x;

		double g_11 = vecx1y1.y;
		double g_12 = vecx1y2.y;
		double g_21 = vecx2y1.y;
		double g_22 = vecx2y2.y;

		double a00 = f_11;
		double a10 = f_21 - f_11;
		double a01 = f_12 - f_11;
		double a11 = f_11 - f_21 - f_12 + f_22;
		double b00 = g_11;
		double b10 = g_21 - g_11;
		double b01 = g_12 - g_11;
		double b11 = g_11 - g_21 - g_12 + g_22;
		double c00 = a11 * b00 - a00 * b11;
		double c10 = a11 * b10 - a10 * b11;
		double c01 = a11 * b01 - a01 * b11;

		if (c01 == 0.0) {
			std::cout << "c01 is 0" << std::endl;
			continue;
		}

		double div = c10 / c01;
		double a = -a11 * div;
		double b = a10 - a01 * div;
		double c = a00 - (a01 + a11) * c00 / c01;
		double s[2];
		bool flag = quadraticRoot(s[0], s[1], a, b, c);
		if (!flag) {
			continue;
		}
		double t[2];
		t[0] = -c00 / c01 - div * s[0];
		t[1] = -c00 / c01 - div * s[1];

		for (int j = 0; j < 2; j++) {
			if (s[j] > -EPSILON && s[j] < 1 + EPSILON &&
				t[j] > -EPSILON && t[j] < 1 + EPSILON) {

				pt.x = s[j];
				pt.y = t[j];
				Singularity point;
				point.p = pt + posx1y1;
				singularities.push_back(point);
				double r0 = f_11 + pt.x * (f_21 - f_11) + pt.y * (f_12 - f_11) + pt.x * pt.y * (f_11 - f_21 - f_12 + f_22);
				std::cout << r0 << std::endl;
			}
		}
	}
}

void classifySingularityByWinding() {
	for (Singularity s : singularities) {
		icVector3 posn = s.p;
		Quad* quad = findQuad(posn);
		double winding_angle = 0;
		double angles[4];
			for (int i = 0; i < 4; i++) {
			auto vi = quad->verts[i]->vec();
			double angle = atan2(vi.entry[1], vi.entry[0]);
			if (angle < 0) {
				angle += 2 * M_PI;
			}
			angles[i] = angle;
		}

		for (int i = 0; i < 4; i++) {
			int next = i + 1;
			if (next > 3) {
				next = 0;
			}
			double diff = (angles[next] - angles[i]);
			if (diff < -M_PI) {
				diff += 2 * M_PI;
			}
			if (diff > -M_PI) {
				diff -= 2 * M_PI;
			}
			winding_angle += diff;
		}
		if (std::abs(winding_angle) < EPSILON) {
			s.type = -1;
		}
		else if (winding_angle < 0) {
			s.type = 2;
			s.rgb = icVector3(0.0, 1.0, 0.0);
		} else if (winding_angle > 0) {
			s.type = 0;
			s.rgb = icVector3(1.0, 0.0, 0.0);
		}

	}
}

void classifySingularity() {
	for (Singularity s : singularities) {
		icVector3 posn = s.p;
		Quad* quad = findQuad(posn);

		int R[4] = { 2, 3, 0, 1 };
		const icVector2 min = icVector2(quad->verts[R[0]]->x, quad->verts[R[0]]->y);
		const icVector2 max = icVector2(quad->verts[R[2]]->x, quad->verts[R[2]]->y);
		double len_x = max.x - min.x;
		double len_y = max.y - min.y;

		double dfdx = (-1 / len_x) * (max.y - posn.y) * quad->verts[R[0]]->vx
			+ (1 / len_x) * (max.y - posn.y) * quad->verts[R[1]]->vx
			+ (-1 / len_x) * (posn.y - min.y) * quad->verts[R[3]]->vx
			+ (1 / len_x) * (posn.y - min.y) * quad->verts[R[2]]->vx;

		double dfdy = (-1 / len_y) * (max.x - posn.x) * quad->verts[R[0]]->vx
			+ (-1 / len_y) * (posn.x - min.x) * quad->verts[R[1]]->vx
			+ (1 / len_y) * (max.x - posn.x) * quad->verts[R[3]]->vx
			+ (1 / len_y) * (posn.x - min.x) * quad->verts[R[2]]->vx;

		double dgdx = (-1 / len_x) * (max.y - posn.y) * quad->verts[R[0]]->vy
			+ (1 / len_x) * (max.y - posn.y) * quad->verts[R[1]]->vy
			+ (-1 / len_x) * (posn.y - min.y) * quad->verts[R[3]]->vy
			+ (1 / len_x) * (posn.y - min.y) * quad->verts[R[2]]->vy;

		double dgdy = (-1 / len_y) * (max.x - posn.x) * quad->verts[R[0]]->vy
			+ (-1 / len_y) * (posn.x - min.x) * quad->verts[R[1]]->vy
			+ (1 / len_y) * (max.x - posn.x) * quad->verts[R[3]]->vy
			+ (1 / len_y) * (posn.x - min.x) * quad->verts[R[2]]->vy;

		s.jacobi.entry[0][0] = dfdx;
		s.jacobi.entry[0][1] = dfdy;
		s.jacobi.entry[1][0] = dgdx;
		s.jacobi.entry[1][1] = dgdy;

		double tr = dfdx + dgdy;
		double det = dfdx * dgdy - dfdy * dgdx;
		double delta = tr * tr - 4 * det;
		if (delta >= 0) {
			double r1 = 0.5 * (tr + sqrt(delta));
			double r2 = 0.5 * (tr - sqrt(delta));
			if (r1 == 0 && r2 == 0) {
				s.type = -1;
			}
			else if (r1 >= 0 && r2 >= 0) {
				s.type = 0;
			}
			else if (r1 <= 0 && r2 <= 0) {
				s.type = 1;
				s.rgb = icVector3(0.0, 0.0, 1.0);
			} 
			if (r1 > 0 && r2 < 0 || r1 < 0 && r2 > 0) {
				s.type = 2;
				s.rgb = icVector3(0.0, 1.0, 0.0);
			}
			else {
				s.type = -1;
			}
		}
		else {
			if (tr == 0) {
				s.type = 3;
				s.rgb = icVector3(0.0, 1.0, 1.0);
			}
			else {
				s.type = 4;
				s.rgb = icVector3(1.0, 1.0, 0.0);
			}
		}
	}
}

void extractSeparatrix() {
	for (Singularity s : singularities) {
		// Saddle
		if (s.type == 2) {
			double a = s.jacobi.entry[0][0];
			double b = s.jacobi.entry[0][1];
			double c = s.jacobi.entry[1][0];
			double d = s.jacobi.entry[1][1];
			double yd = (a + d) / 2;
			double yr = (c - b) / 2;
			double ys = std::sqrt((a - d) * (a - d) + (b + c) * (b + c)) / 2;
			double theta = std::atan2(b + c, a - d);
			double phi = std::atan(yr / ys);
			double a_cos = std::cos(theta / 2);
			double b_sin = -std::sin(theta / 2);
			double c_sin = std::sin(theta / 2);
			// bug here? shouldn't this be negative?
			double d_cos = std::cos(theta / 2);
			double a_sin_phi = std::sqrt(std::sin(phi + M_PI / 4));
			double a_cos_phi = std::sqrt(std::cos(phi + M_PI / 4));
			icVector3 maj_v = icVector3(0.0);
			icVector3 min_v = icVector3(0.0);
			maj_v.x = a_cos * (a_sin_phi + a_cos_phi) + b_sin * (a_sin_phi - a_cos_phi);
			maj_v.y = c_sin * (a_sin_phi + a_cos_phi) + d_cos * (a_sin_phi - a_cos_phi);
			min_v.y = a_cos * (a_sin_phi - a_cos_phi) + b_sin * (a_sin_phi + a_cos_phi);
			min_v.y = c_sin * (a_sin_phi - a_cos_phi) + d_cos * (a_sin_phi + a_cos_phi);
			double k = MIN_K / maj_v.length();

			Polyline2 separatrix;
			streamlineFB(separatrix, s.p + k * maj_v, STEP);
			separatrix.rgb = icVector3(1.0, 0.0, 0.0);
			polylines.push_back(separatrix);
			separatrix.clear();

			streamlineFB(separatrix, s.p - k * maj_v, STEP);
			separatrix.rgb = icVector3(1.0, 0.0, 0.0);
			polylines.push_back(separatrix);
			separatrix.clear();

			k = MIN_K / min_v.length();

			streamlineFB(separatrix, s.p + k * min_v, STEP, false);
			separatrix.rgb = icVector3(0.0, 0.0, 1.0);
			polylines.push_back(separatrix);
			separatrix.clear();

			streamlineFB(separatrix, s.p - k * min_v, STEP, false);
			separatrix.rgb = icVector3(0.0, 0.0, 1.0);
			polylines.push_back(separatrix);

		}
	}
}

void displaySingularities() {
	CHECK_GL_ERROR();
	for (Singularity sing : singularities) {
		GLUquadric* quadric = gluNewQuadric();
		glPushMatrix();
		glTranslated(sing.p.x, sing.p.y, sing.p.z);
		glColor3f(sing.p.x, sing.p.y, sing.p.z);
		gluSphere(quadric, 0.1, 16, 16);
		glPopMatrix();
		gluDeleteQuadric(quadric);
	}
	glDisable(GL_BLEND);
}