#pragma once

#include <list>
#include <vector>
#include "polyline2.h"
#include "polyhedron.h"
#include "icMatrix.H"

bool isZero(double x);

void findMinMaxField(icVector3& min, icVector3& max);
icVector3 getVector(Quad* quad, const icVector3 p);

bool insideQuad(const Quad* quad, const icVector3 p);
Quad* findQuad(const icVector3 p);

void streamlineFB(Polyline2& polyline, icVector3 seed, const double step, bool forward = true);
void streamline(Polyline2& polyline, icVector3 seed, const double step);

struct Singularity {
	int type = -1;
	icVector3 p;
	icVector3 rgb = icVector3(0.0);
	icMatrix2x2 jacobi;
};

void displaySingularities(std::list<Singularity> singularities);
void extractSingularity();
void classifySingularity();
void classifySingularityByWinding();
void extractSeparatrix();

