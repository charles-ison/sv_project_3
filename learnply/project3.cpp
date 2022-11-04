#include "project3.h"
#include "polyhedron.h"
#include "iostream"
#include "GL/freeglut.h"
#include <vector>
#include <set>
#include <map>
#include "polyline2.h"

#define EPSILON 0.01
extern Polyhedron* poly;
extern std::vector<Polyline2> polylines;
extern std::list<Singularity> singularities;

int numberOfContours = 100;

double findMin() {
	double min = poly->vlist[0]->scalar;
	for (int i = 1; i < poly->nverts; i++) {
		double vertexScalar = poly->vlist[i]->scalar;
		if (vertexScalar < min) {
			min = vertexScalar;
		}
	}
	return min;
}

double findMax() {
	double max = poly->vlist[0]->scalar;
	for (int i = 1; i < poly->nverts; i++) {
		double vertexScalar = poly->vlist[i]->scalar;
		if (vertexScalar > max) {
			max = vertexScalar;
		}
	}
	return max;
}

icVector3 convertHSVToRGB(icVector3 hsv) {
	double h = hsv.x;
	double s = hsv.y;
	double v = hsv.z;
	double r, g, b;

	double c = s * v;
	double x = c * (1 - abs(fmod(h / 60.0, 2) - 1));
	double m = v - c;

	if (h >= 0 && h < 60) {
		r = c, g = x, b = 0;
	}
	else if (h >= 60 && h < 120) {
		r = x, g = c, b = 0;
	}
	else if (h >= 120 && h < 180) {
		r = 0, g = c, b = x;
	}
	else if (h >= 180 && h < 240) {
		r = 0, g = x, b = c;
	}
	else if (h >= 240 && h < 300) {
		r = x, g = 0, b = c;
	}
	else {
		r = c, g = 0, b = x;
	}

	return icVector3(r + m, g + m, b + m);
}

icVector3 convertRGBToHSV(icVector3 rgb) {
	double r = rgb.x;
	double g = rgb.y;
	double b = rgb.z;
	double h, s, v;

	double cmax = std::max(r, std::max(g, b));
	double cmin = std::min(r, std::min(g, b));
	double diff = cmax - cmin;

	if (cmax == cmin) {
		h = 0;
	} 
	else if (cmax == r) {
		h = fmod(60 * ((g - b) / diff) + 360, 360);
	}
	else if (cmax == g) {
		h = fmod(60 * ((b - r) / diff) + 120, 360);
	}
	else if (cmax == b) {
		h = fmod(60 * ((r - g) / diff) + 240, 360);
	}

	if (cmax == 0) {
		s = 0;
	}
	else {
		s = (diff / cmax);
	}

	v = cmax;

	return icVector3(h, s, v);

}

void visualizeHeight() {
	double min = findMin();
	double max = findMax();

	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double vertexScalar = vertex->scalar;
		double height = (vertexScalar - min) / (max - min);
		vertex->z = 15 * height;
	}
}

//TODO: this is a hack and should be removed in favor of constructing a polyline from vertexes in future iterations
void flattenPolyhedron() {
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		vertex->z = 0;
	}
}

void part2A() {
	polylines.clear();
	double min = findMin();
	double max = findMax();
	double interval = (max - min) / numberOfContours;
	for (int i = 0; i < numberOfContours; i++) {
		std::list<Polyline2> edges = marchingSquare(*poly, min+(i*interval));
		std::vector<Polyline2> newPolylines = makePolylineFromEdges(edges);
		for (auto polyline : newPolylines) {
			polyline.rgb = icVector3(1, 0, 0);
			polylines.push_back(polyline);
		}
	}
	glutPostOverlayRedisplay();
}

void part2B() {
	polylines.clear();
	double min = findMin();
	double max = findMax();
	double interval = (max - min) / numberOfContours;
	for (int i = 0; i < numberOfContours; i++) {
		std::list<Polyline2> edges = marchingSquare(*poly, min + (i * interval));
		std::vector<Polyline2> newPolylines = makePolylineFromEdges(edges);
		for (auto polyline : newPolylines) {
			icVector3 redRGB(1.0, 0.0, 0.0);
			icVector3 blueRGB(0.0, 0.0, 1.0);
			icVector3 redHSV = convertRGBToHSV(redRGB);
			icVector3 blueHSV = convertRGBToHSV(blueRGB);

			double doubleNumberOfContours = 1.0 * (numberOfContours - 1);
			double redScalar = (1.0 * i) / (doubleNumberOfContours);
			double blueScalar = (doubleNumberOfContours - (1.0 * i)) / doubleNumberOfContours;

			icVector3 newHSV = (redHSV * redScalar) + (blueHSV * blueScalar);
			polyline.rgb = convertHSVToRGB(newHSV);

			polylines.push_back(polyline);
		}
	}
	glutPostOverlayRedisplay();
}

void part2C() {
	polylines.clear();
	visualizeHeight();
	polylines.clear();
	double min = findMin();
	double max = findMax();
	double interval = (max - min) / numberOfContours;
	for (int i = 0; i < numberOfContours; i++) {
		std::list<Polyline2> edges = marchingSquare(*poly, min + (i * interval));
		std::vector<Polyline2> newPolylines = makePolylineFromEdges(edges);
		for (auto polyline : newPolylines) {
			polyline.rgb = icVector3(1, 0, 0);
			polylines.push_back(polyline);
		}
	}
	glutPostOverlayRedisplay();
}

void part2D() {
	polylines.clear();
	visualizeHeight();
	polylines.clear();
	double min = findMin();
	double max = findMax();
	double interval = (max - min) / numberOfContours;
	for (int i = 0; i < numberOfContours; i++) {
		std::list<Polyline2> edges = marchingSquare(*poly, min + (i * interval));
		std::vector<Polyline2> newPolylines = makePolylineFromEdges(edges);
		for (auto polyline : newPolylines) {
			icVector3 redRGB(1.0, 0.0, 0.0);
			icVector3 blueRGB(0.0, 0.0, 1.0);
			icVector3 redHSV = convertRGBToHSV(redRGB);
			icVector3 blueHSV = convertRGBToHSV(blueRGB);

			double doubleNumberOfContours = 1.0 * numberOfContours;
			double redScalar = (1.0 * i) / (doubleNumberOfContours);
			double blueScalar = (doubleNumberOfContours - (1.0 * i)) / doubleNumberOfContours;
			icVector3 newHSV = (redHSV * redScalar) + (blueHSV * blueScalar);
			polyline.rgb = convertHSVToRGB(newHSV);

			polylines.push_back(polyline);
		}
	}
	glutPostOverlayRedisplay();
}

std::list<CriticalPoint> getCriticalPoints() {
	std::list<CriticalPoint> criticalPoints;

	for (int i = 0; i < poly->nverts; i++) {
		std::set<Relationship> relationships;
		Vertex* potentialCriticalPoint = poly->vlist[i];
		if (potentialCriticalPoint->nquads != 4) {
			continue;
		}
		for (int j = 0; j < potentialCriticalPoint->nquads; j++) {
			Quad* quad = potentialCriticalPoint->quads[j];
			for (int k = 0; k < 4; k++) {
				Vertex* vertex = quad->verts[k];
				if (potentialCriticalPoint->scalar < vertex->scalar) {
					relationships.insert(min);
				}
				else if (potentialCriticalPoint->scalar > vertex->scalar) {
					relationships.insert(max);
				}				
			}
		}
		if (relationships.size() == 1) {
			Vertex criticalPoint = *potentialCriticalPoint;
			icVector3 vector = icVector3(criticalPoint.x, criticalPoint.y, criticalPoint.z);
			if (*relationships.begin() == min) {
				icVector3 color = icVector3(0.0, 0.0, 0.9);
				criticalPoints.push_back(CriticalPoint(vector, color, criticalPoint.scalar, min));
			}
			else if (*relationships.begin() == max) {
				icVector3 color = icVector3(0.9, 0.0, 0.0);
				criticalPoints.push_back(CriticalPoint(vector, color, criticalPoint.scalar, max));
			}
		}
	}

	for (int i = 0; i < poly->nquads; i++) {
		Quad* quad = poly->qlist[i];
		Vertex* x2y2 = quad->verts[0];
		Vertex* x1y2 = quad->verts[1];
		Vertex* x1y1 = quad->verts[2];
		Vertex* x2y1 = quad->verts[3];

		double x1 = x1y2->x;
		double x2 = x2y2->x;
		double y1 = x1y1->y;
		double y2 = x1y2->y;

		double x1y1Scalar = x1y1->scalar;
		double x1y2Scalar = x1y2->scalar;
		double x2y1Scalar = x2y1->scalar;
		double x2y2Scalar = x2y2->scalar;

		double x0 = (x2 * x1y1Scalar - x1 * x2y1Scalar - x2 * x1y2Scalar + x1 * x2y2Scalar) / (x1y1Scalar - x2y1Scalar - x1y2Scalar + x2y2Scalar);
		double y0 = (y2 * x1y1Scalar - y2 * x2y1Scalar - y1 * x1y2Scalar + y1 * x2y2Scalar) / (x1y1Scalar - x2y1Scalar - x1y2Scalar + x2y2Scalar);

		if (x0-EPSILON > x1 && x0+EPSILON < x2 && y0-EPSILON > y1 && y0+EPSILON < y2) {
			double zAverage = (x2y2->z + x1y2->z + x2y1->z + x2y2->z) / 4;
			double linearlyInterpolatedScalar =
				(((x2 - x0) / (x2 - x1)) * ((y2 - y0) / (y2 - y1)) * x1y1Scalar) +
				(((x0 - x1) / (x2 - x1)) * ((y2 - y0) / (y2 - y1)) * x2y1Scalar) +
				(((x2 - x0) / (x2 - x1)) * ((y0 - y1) / (y2 - y1)) * x1y2Scalar) +
				(((x0 - x1) / (x2 - x1)) * ((y0 - y1) / (y2 - y1)) * x2y2Scalar);

			icVector3 vector = icVector3(x0, y0, zAverage);
			icVector3 color = icVector3(0.70, 0.70, 0.70);
			criticalPoints.push_back(CriticalPoint(vector, color, linearlyInterpolatedScalar, saddle));
		}
	}
	return criticalPoints;
}

void addCriticalPointContours(double threshold, icVector3 color) {
	std::list<Polyline2> edges = marchingSquare(*poly, threshold);
	std::vector<Polyline2> newPolylines = makePolylineFromEdges(edges);
for (auto polyline : newPolylines) {
	polyline.rgb = color;
	polylines.push_back(polyline);
}
}

void part3B(std::list<CriticalPoint> criticalPoints) {
	for (CriticalPoint criticalPoint : criticalPoints) {
		if (criticalPoint.relationship == saddle) {
			addCriticalPointContours(criticalPoint.scalar, criticalPoint.color);
		}
	}
}

// Actual Project 3 Stuff
void streamline(Polyline2 polyline, icVector3 seed, const double step) {
	streamlineFB(polyline, seed, step);
	Polyline2 lineBack;
	streamlineFB(lineBack, seed, step, false);
	polyline.merge(lineBack);
}

void streamlineFB(Polyline2 polyline, icVector3 seed, const double step, bool forward) {
	polyline.vertices.push_back(seed);
	Quad quad = findQuad(seed);
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
		icVector3 nextPosition = currPos + step * currVec * coef;
		if (sinp2Boundary(nextPos, min, max)) {
			polyline.vertices.push_back(nextPosition);
			break;
		}
		Quad nextQuad = nullptr;
		streamlineTrace(nextQuad, quad, currentPosition, currentVector * coef, step, min, max);
		quad = nextQuad;
		currentPosition = nextPosition;
		polyline.vertices.push_back(currentPosition);
	}
}

// not done
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
	}
}

bool insideQuad(const Quad quad, const icVector3 p) {
	double v0x = quad.verts[2]->x;
	double v0y = quad.verts[2]->y;
	double v2x = quad.verts[0]->x;
	double v2y = quad.verts[0]->y;
	if (p.x >= v0x && p.x <= v2x && p.y >= v0y && p.y <= v2y) {
		return true;
	}
	else {
		return false;
	}
}

Quad findQuad(const icVector3 p) {
	for (int i = 0; i < poly->nquads; i++) {
		Quad tempQuad = poly->qlist[i];
		if (insideQuad(tempQuad, p)) {
			return tempQuad;
		}
	}
	return nullptr;
}

void streamlineTrace(Quad nextQuad, icVector3 nextPos, icVector3 nextVec, Quad currentQuad, icVector3 currentPos,
	icVector3 currentVec, double t, const icVector3 min, const icVector3 max) {

	bool insideQuad = false;
	while (!insideQuad) {
		if (currentPos.x < min.x || currentPos.x > max.x) ||
			currentPos.y < min.y || currentPos.y > max.y) {
			nextQuad = nullptr; 
			return;
		}

		double t_ = INFINITY;
		Quad nextQuad = nullptr;
		for (int i = 0; i < 4; i++) {
			Edge edge = currentQuad.edges[i];
			Vertex v0 = edge.verts[0];
			Vertex v1 = edge.verts[1];
			double temp;
			if (std::abs(v0.x - v1.x) < EPSILON) {
				temp = (v0.x - currentPos.x) / currentVec.x;
			}
			else {
				temp = (v0.y - currentPos.y) / currentVec.y;
			}
			if (temp > 0 && temp < t) {
				t_ = temp;
				if (edge.quads[0] != currentQuad && edge.quads[1] == currentQuad) {
					nextQuad = edge.quads[0];
				}
				else if (edge.quads[0] == currQuad && edge.quads[1] != currentQuad) {
					nextQuad = edge.quads[1];
				}
			}
		}
		if (nextQuad == nullptr) {
			currentPos = currentPos + currentVec * t;
			nextQuad = findQuad(currentPos);
			return;
		}
		else {
			if (t_ >= t) {
				insideQuad = true;
			} 
			else {
				currentQuad = nextQuad;
				t = t - t_;
				currentPos = currentPos + currentVec * t_;
			}
		}
	}
	nextQuad = currentQuad;
}

bool singRoot(double r0, double r1, double a, double b, double c, double d) {
	double f0 = b - a - (c + d);
	double f1 = (c + d);
	double f2 = a;
	return quadraticRoot(r0, r1, f0, f1, f2);
}

// not done
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
		double b01 - g_12 - g_11;
		double b11 = g_11 - g_21 - g_12 + g_22;
		double c00 = a11 * b00 - a00 * b11;
		double c10 = a11 * b10 - a10 * b11;
		double c01 = a11 * b01 - a01 * b11;

		if (c01 == 0.0) {
			std::cout << "c01 is 0" << std::endl;
			continue;
		}
	}
}