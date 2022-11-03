#include "icVector.H"
#include <list>
#include <vector>
#include "polyhedron.h"

class Polyline2 {
	public:
		std::list<icVector3> vertices;
		icVector3 rgb = icVector3(1.0, 0.0, 0.0);
		double weight = 1.0;
		double scalar = 0.0;
		bool isNeighbor(const Polyline2 polyline);
		void merge(const Polyline2 polyline);
};

std::list<Polyline2> marchingSquare(const Polyhedron& polyhedron, const double threshold);
void displayPolylines(std::vector<Polyline2> polylines);
std::vector<Polyline2> makePolylineFromEdges(std::list<Polyline2> edges);
