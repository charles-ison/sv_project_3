#pragma once
#include "polyline2.h"
#include "GL/glew.h"
#include "iostream"
#include <vector>
#define EPSILON 1.0e-100

bool Polyline2::isNeighbor(const Polyline2 polyline) {
	return (vertices.front() - polyline.vertices.front()).length() < EPSILON ||
		(vertices.front() - polyline.vertices.back()).length() < EPSILON ||
		(vertices.back() - polyline.vertices.front()).length() < EPSILON ||
		(vertices.back() - polyline.vertices.back()).length() < EPSILON;
}

void Polyline2::merge(const Polyline2 polyline) {
	Polyline2 polylineCopy = polyline;
	if ((vertices.front() - polyline.vertices.front()).length() < EPSILON) {
		polylineCopy.vertices.pop_front();
		for (auto i = polylineCopy.vertices.begin(); i != polylineCopy.vertices.end(); i++) {
			vertices.push_front(*i);
		}
	}
	else if ((vertices.front() - polyline.vertices.back()).length() < EPSILON) {
		polylineCopy.vertices.pop_back();
		polylineCopy.vertices.reverse();
		for (auto i = polylineCopy.vertices.begin(); i != polylineCopy.vertices.end(); i++) {
			vertices.push_front(*i);
		}
	}
	else if ((vertices.back() - polyline.vertices.front()).length() < EPSILON) {
		polylineCopy.vertices.pop_front();
		for (auto i = polylineCopy.vertices.begin(); i != polylineCopy.vertices.end(); i++) {
			vertices.push_back(*i);
		}
	}
	else if ((vertices.back() - polyline.vertices.back()).length() < EPSILON) {
		polylineCopy.vertices.pop_back();
		polylineCopy.vertices.reverse();
		for (auto i = polylineCopy.vertices.begin(); i != polylineCopy.vertices.end(); i++) {
			vertices.push_back(*i);
		}
	}
}

void displayPolylines(std::vector<Polyline2> polylines) {
	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	for (Polyline2 polyline : polylines) {
		glLineWidth(polyline.weight);
		glColor3f(polyline.rgb.entry[0], polyline.rgb.entry[1], polyline.rgb.entry[2]);
		glBegin(GL_LINE_STRIP);

		for (auto it = polyline.vertices.begin(); it != polyline.vertices.end(); it++) {
			glVertex3d(it->entry[0], it->entry[1], it->entry[2]);
		}
		glEnd();
	}

	glDisable(GL_BLEND);
	glLineWidth(1);
}


Vertex lineInterpolateByScalar(const Vertex v0, const Vertex v1, const double threshold) {
	double f0_f1 = v0.scalar - v1.scalar;
	Vertex r(0.0, 0.0, 0.0);
	if (std::abs(f0_f1) < EPSILON) {
		r.x = (v0.x + v1.x) / 2;
		r.y = (v0.y + v1.y) / 2;
		r.z = (v0.z + v1.z) / 2;
	}
	else {
		double t = std::abs((v0.scalar - threshold) / ((v0.scalar - threshold) - (v1.scalar - threshold)));
		r.x = v0.x + t * (v1.x - v0.x);
		r.y = v0.y + t * (v1.y - v0.y);
		r.z = v0.z + t * (v1.z - v0.z);
	}
	return r;
}

void lookUpTable(std::vector<Vertex>& r, const Vertex& v0, const Vertex& v1, const Vertex& v2, const Vertex& v3, const double threshold) {
	r.reserve(2);
	int id = 0;
	if (v0.scalar <= threshold + EPSILON) {
		id += 1;
	}
	if (v1.scalar <= threshold + EPSILON) {
		id += 2;
	}
	if (v2.scalar <= threshold + EPSILON) {
		id += 4;
	}
	if (v3.scalar <= threshold + EPSILON) {
		id += 8;
	}
	double center = 0;
	switch (id) {
		case 0:
			break;
		case 1:
			r.push_back(lineInterpolateByScalar(v0, v1, threshold));
			r.push_back(lineInterpolateByScalar(v0, v3, threshold));
			break;
		case 2:
			r.push_back(lineInterpolateByScalar(v0, v1, threshold));
			r.push_back(lineInterpolateByScalar(v1, v2, threshold));
			break;
		case 3:
			r.push_back(lineInterpolateByScalar(v1, v2, threshold));
			r.push_back(lineInterpolateByScalar(v0, v3, threshold));
			break;
		case 4:
			r.push_back(lineInterpolateByScalar(v1, v2, threshold));
			r.push_back(lineInterpolateByScalar(v2, v3, threshold));
			break;
		case 5:
			center = v0.scalar + v1.scalar + v2.scalar + v3.scalar;
			center /= 4;
			if (center <= threshold) {
				r.push_back(lineInterpolateByScalar(v0, v1, threshold));
				r.push_back(lineInterpolateByScalar(v1, v2, threshold));
				r.push_back(lineInterpolateByScalar(v2, v3, threshold));
				r.push_back(lineInterpolateByScalar(v0, v3, threshold));
			}
			else {
				r.push_back(lineInterpolateByScalar(v0, v1, threshold));
				r.push_back(lineInterpolateByScalar(v0, v3, threshold));
				r.push_back(lineInterpolateByScalar(v1, v2, threshold));
				r.push_back(lineInterpolateByScalar(v2, v3, threshold));
			}
			break;
		case 6:
			r.push_back(lineInterpolateByScalar(v0, v1, threshold));
			r.push_back(lineInterpolateByScalar(v2, v3, threshold));
			break;
		case 7:
			r.push_back(lineInterpolateByScalar(v2, v3, threshold));
			r.push_back(lineInterpolateByScalar(v0, v3, threshold));
			break;
		case 8:
			r.push_back(lineInterpolateByScalar(v2, v3, threshold));
			r.push_back(lineInterpolateByScalar(v0, v3, threshold));
			break;
		case 9:
			r.push_back(lineInterpolateByScalar(v0, v1, threshold));
			r.push_back(lineInterpolateByScalar(v2, v3, threshold));
			break;
		case 10:
			center = v0.scalar + v1.scalar + v2.scalar + v3.scalar;
			center /= 4;
			if (center <= threshold) {
				r.push_back(lineInterpolateByScalar(v0, v1, threshold));
				r.push_back(lineInterpolateByScalar(v0, v3, threshold));
				r.push_back(lineInterpolateByScalar(v1, v2, threshold));
				r.push_back(lineInterpolateByScalar(v2, v3, threshold));
			}
			else {
				r.push_back(lineInterpolateByScalar(v0, v1, threshold));
				r.push_back(lineInterpolateByScalar(v1, v2, threshold));
				r.push_back(lineInterpolateByScalar(v2, v3, threshold));
				r.push_back(lineInterpolateByScalar(v0, v3, threshold));
			}
			break;
		case 11:
			r.push_back(lineInterpolateByScalar(v1, v2, threshold));
			r.push_back(lineInterpolateByScalar(v2, v3, threshold));
			break;
		case 12:
			r.push_back(lineInterpolateByScalar(v1, v2, threshold));
			r.push_back(lineInterpolateByScalar(v0, v3, threshold));
			break;
		case 13:
			r.push_back(lineInterpolateByScalar(v0, v1, threshold));
			r.push_back(lineInterpolateByScalar(v1, v2, threshold));
			break;
		case 14:
			r.push_back(lineInterpolateByScalar(v0, v1, threshold));
			r.push_back(lineInterpolateByScalar(v0, v3, threshold));
			break;
		case 15:
			break;
	}
}

std::list<Polyline2> marchingSquare(const Polyhedron& polyhedron, const double threshold) {
	std::list<Polyline2> edges;
	for (int i = 0; i < polyhedron.nquads; i++) {
		std::vector<Vertex> r;
		lookUpTable(
			r,
			*polyhedron.qlist[i]->verts[0],
			*polyhedron.qlist[i]->verts[1],
			*polyhedron.qlist[i]->verts[2],
			*polyhedron.qlist[i]->verts[3],
			threshold);

		if (r.size() > 0) {
			for (int j = 0; j < r.size()/2; j++) {
				Polyline2 polyline;
				auto v0 = icVector3(r[j * 2].x, r[j * 2].y, r[j * 2].z);
				auto v1 = icVector3(r[j * 2 + 1].x, r[j * 2 + 1].y, r[j * 2 + 1].z);
				polyline.vertices.push_back(v0);
				polyline.vertices.push_back(v1);
				edges.push_back(polyline);
			}
		}
	}
	return edges;
}

std::vector<Polyline2> makePolylineFromEdges(std::list<Polyline2> edges) {

	std::vector<Polyline2> newPolylines;
	newPolylines.reserve(edges.size());
	std::list<Polyline2> edgesTemp(edges);
	while (edgesTemp.size() > 0) {
		newPolylines.push_back(edgesTemp.front());
		edgesTemp.erase(edgesTemp.begin());
		int initSize = 0;
		while (initSize != edgesTemp.size()) {
			initSize = edgesTemp.size();
			for (auto i = edgesTemp.begin(); i != edgesTemp.end();) {
				if (newPolylines.back().isNeighbor(*i)) {
					newPolylines.back().merge(*i);
					i = edgesTemp.erase(i);
				}
				else {
					i++;
				}
			}
		}
	}
	return newPolylines;
}
