#include <list>
#include <vector>
#include "polyhedron.h"

void part2A();
void part2B();
void part2C();
void part2D();

enum Relationship { min, max, saddle };

class CriticalPoint {
public:
	icVector3 vector;
	icVector3 color;
	double scalar;
	Relationship relationship;
public:
	CriticalPoint(icVector3 newVector, icVector3 newColor, double newScalar, Relationship newRelationship) {
		color = newColor;
		vector = newVector;
		scalar = newScalar;
		relationship = newRelationship;
	}
};

std::list<CriticalPoint> getCriticalPoints();
void part3B(std::list<CriticalPoint> criticalPoints);
void flattenPolyhedron();