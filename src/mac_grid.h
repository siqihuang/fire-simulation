#ifndef MACGrid_H_
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)

#include "open_gl_headers.h" // PETER KUTZ.

#include <glm/glm.hpp>

#include "grid_data.h"
#include "grid_data_matrix.h" // PETER KUTZ.
#include "spectrum.h"

class Camera;

class MACGrid
{

public:
	MACGrid();
	~MACGrid();
	MACGrid(const MACGrid& orig);
	MACGrid& operator=(const MACGrid& orig);

	void reset();

	void draw(const Camera& c);
	void updateSources();
	void advectVelocity(double dt,bool tag);
	void addExternalForces(double dt,bool tag);
	void project(double dt);
	void advectTemperature(double dt,bool tag);
	void advectDensity(double dt);
	double getUgradientH(int i, int j, int k);

	//fire simulation

	void initBoundary();
	void advectBoundary(double dt,bool tag);
	void advectY(double dt);
	void findBoundary(bool tag);
	void findBoundaryGradient();
	void reconditioning();
	void insideBoundary();
	bool isOnBorder(int i,int j,int k);

	//fire simulation

protected:

	// Setup
	void initialize();

	// Simulation
	void computeBouyancy(double dt);
	void computeVorticityConfinement(double dt);

	// Rendering
	struct Cube { glm::dvec3 pos; glm::dvec4 color; double dist; };
	void drawWireGrid();
	void drawSmokeCubes(const Camera& c);
	void drawSmoke(const Camera& c);
	void drawCube(const MACGrid::Cube& c);
	void drawFace(const MACGrid::Cube& c);
	void drawVelocities();

	glm::dvec4 getRenderColor(int i, int j, int k);
	glm::dvec4 getRenderColor(const glm::dvec3& pt);
	glm::dvec3 getRenderRGB(double t);

	void drawZSheets(bool backToFront);
	void drawXSheets(bool backToFront);

	// GridData accessors
	enum Direction { X, Y, Z };
	glm::dvec3 getVelocityf(const glm::dvec3& pt);
	glm::dvec3 getVelocityh(const glm::dvec3& pt);
	glm::dvec3 getCenter(int i, int j, int k);
	bool isValidCell(int i, int j, int k);

	 // Sets up A matrix for calculation
	void calculateAMatrix();
	void calculateAMatrixOut();

	// Conjugate Gradient stuff
	bool conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
	double dotProduct(const GridData & vector1, const GridData & vector2);
	void add(const GridData & vector1, const GridData & vector2, GridData & result);
	void subtract(const GridData & vector1, const GridData & vector2, GridData & result);
	void multiply(const double scalar, const GridData & vector, GridData & result);
	double maxMagnitude(const GridData & vector);
	void apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result);


	// TODO : Fill in the necessary data structures to maintain velocity, pressure and density
	GridDataX mUh; 
	GridDataY mVh; 
	GridDataZ mWh;
	GridDataX mUf;
	GridDataY mVf;
	GridDataZ mWf;
	GridData mPressure;  
	GridData mDensityh;
	GridData mDensityf;
	GridData mTemp;
	GridData mInside;
	GridDataMatrix AMatrix;
	
	//use for the color rendering
	Spectrum spect;


	//fire simulation
	GridData mY;//record the temperature information
	GridData mPhi;//the bundary of the blue core surface
	GridData mN_x;//bundary normal x
	GridData mN_y;//bundary normal y
	GridData mN_z;//bundary normal z
	GridData mPhi_State;//indicate whether the point is inside or outside
	GridData mNeighborX;
	GridData mNeighborY;
	GridData mNeighborZ;

	int unTagedNum;

public:

	enum RenderMode { CUBES, SHEETS };
	static RenderMode theRenderMode;
	static bool theDisplayVel;
	
	void saveSmoke(const char* fileName);
};

#endif
