#include "grid_data.h"

GridData::GridData() 
{
	// TODO : GridData constructor
	mDflt = 0.0;
	mBound = glm::dvec3(0.0, 0.0, 0.0);
}

GridData::GridData(const GridData& orig) 
{
	// TODO : GridData copy constructor
	mDflt = orig.mDflt; 
	mBound = orig.mBound; 
	mData = orig.mData; 
}

GridData::~GridData() 
{
  // TODO : GridData destructor
}

GridData& GridData::operator=(const GridData& orig) 
{
	// TODO : Override GridData '=' operator with copy 
	if(this == &orig)
	{
		return *this;
	}
	mDflt = orig.mDflt;
	mData = orig.mData;
	mBound = orig.mBound;
	return *this;
}

void GridData::initialize(double dfltValue) 
{
	// TODO : Initialize the grid to a default value
	mDflt = dfltValue;	
	for(int i = 0; i < 3; ++i)
	{
		mBound[i] = theCellSize * (theDim[i]);
	}
	mData.resize(theDim[0] * theDim[1] * theDim[2], false);
	std::fill(mData.begin(), mData.end(), mDflt);
}

bool GridData::isValidIndex(const int index, const int axis) const
{
	return axis <=2 && index >= 0 && index < theDim[axis];
}

double& GridData::operator()(int i, int j, int k) 
{
	static double dflt = 0;
	dflt = mDflt;  
	// TODO : Grid accessor that allows for client to access and set cell data
	if(!isValidIndex(i, 0) || !isValidIndex(j, 1) || !isValidIndex(k, 2))
	{
		return dflt;
	}
	//col + row + stack
	return mData[i + k * theDim[0] + j * theDim[0] * theDim[2]];
}

const double GridData::operator()(int i, int j, int k) const 
{
	static double dflt = 0;
	dflt = mDflt;  
	// TODO : Grid accessor
	if(!isValidIndex(i, 0) || !isValidIndex(j, 1) || !isValidIndex(k, 2))
	{
		return dflt;
	}
	//col + row + stack
	return mData[i + k * theDim[0] + j * theDim[0] * theDim[2]];
}

double GridData::interpolate(const glm::dvec3& pt) 
{
	// TODO : Given a point, interpolate the value in the grid at that point
	//linear 
	glm::dvec3 pos = worldToSelf(pt);
	int i = (int) (pos[0]/theCellSize);
	int j = (int) (pos[1]/theCellSize);
	int k = (int) (pos[2]/theCellSize);

	double s = 1.0 / theCellSize;  
	double fx,fy,fz;
	fx = s*(pos[0] - i * theCellSize);
	fy = s*(pos[1] - j * theCellSize);
	fz = s*(pos[2] - k * theCellSize);

	double tmp[9];
	int idx = 0;
	//first
	for(int kk = k; kk <= k + 1; ++kk)
	{
		for(int ii = i; ii <= i + 1; ++ii)
		{
			for(int jj = j; jj <= j + 1; ++jj)
			{
				tmp[++idx] = (*this)(ii,jj,kk); 
			}
		}
	}
	double tmp2[4];
	idx = 0;
	//second
	for(int ii = 1; ii < 9; ii +=2 )
	{
		tmp2[idx++] = LERP(tmp[ii], tmp[ii+1], fy); 
	}
	//third
	double a = LERP (tmp2[0], tmp2[1], fx);
	double b = LERP (tmp2[2], tmp2[3], fx);
	return LERP(a, b, fz);
}

std::vector<double>& GridData::data() 
{
	// TODO : Return underlying data structure (you may change the method header
	// to fit whatever design you choose).
	return mData;
}

void GridData::getCell(const glm::dvec3& pt, int& i, int& j, int& k) 
{
	// TODO : Given a point in world coordinates, return the cell index
	// corresponding to it.
	glm::dvec3 tmp = worldToSelf(pt);
	i = (int)(tmp.x / theCellSize); 
	j = (int)(tmp.y / theCellSize); 
	k= (int)(tmp.z / theCellSize); 
}

glm::dvec3 GridData::worldToSelf(const glm::dvec3& pt) const 
{
	// TODO : Given a point, returns the cell index that the grid uses in its own
	// space.
	glm::dvec3 res;
	for(int i = 0; i < 3; ++i)
	{
		res[i] = min(max(0.0, pt[i] - theCellSize * 0.5), mBound[i]);
	}
	return res; 
}

GridDataX::GridDataX() : GridData() {
}

GridDataX::~GridDataX() {
}

void GridDataX::initialize(double dfltValue) 
{
	// TODO : Intialize GridDataX
	GridData::initialize(dfltValue);
	mBound[0] = theCellSize * (theDim[0]+1);
	mBound[1] = theCellSize * theDim[1];
	mBound[2] = theCellSize * theDim[2];
	mData.resize((theDim[0] + 1) * theDim[1] * theDim[2], false);
	std::fill(mData.begin(), mData.end(), mDflt);
}

double& GridDataX::operator()(int i, int j, int k) 
{
  // TODO : GridX accessor
	static double dflt = 0;
	dflt = mDflt;  
	if (i < 0 || i > theDim[0]) 
		return dflt;
	j = max(j, 0);
	j = min(j, theDim[1]-1);
	k = max(k, 0);
	k = min(k, theDim[2]-1);
	int col = i;
	int row = k * (theDim[0]+1);
	int stack = j * (theDim[0]+1) * theDim[2];
	return mData[stack + row + col];
}

const double GridDataX::operator()(int i, int j, int k) const 
{
	// TODO : GridX accessor
	return 0.0;
	static double dflt = 0;
	dflt = mDflt;  
	if (i < 0 || i > theDim[0]) 
	{
		return dflt;
	}
	j = max(j, 0);
	j = min(j, theDim[1]-1);
	k = max(k, 0);
	k = min(k, theDim[2]-1);
	int col = i;
	int row = k * (theDim[0]+1);
	int stack = j * (theDim[0]+1) * theDim[2];
	return mData[stack + row + col];
}

glm::dvec3 GridDataX::worldToSelf(const glm::dvec3& pt) const 
{
  // TODO : Given a point, returns the cell index that the grid uses in its own space
   glm::dvec3 res;
   res[0] = min(max(0.0, pt[0]), mBound[0]);
   res[1] = min(max(0.0, pt[1] - theCellSize * 0.5), mBound[1]);
   res[2] = min(max(0.0, pt[2] - theCellSize * 0.5), mBound[2]);
   return res;
}

GridDataY::GridDataY() : GridData() {
}

GridDataY::~GridDataY() {
}

void GridDataY::initialize(double dfltValue) 
{
  // TODO : Intialize GridDataY
	GridData::initialize(dfltValue);
	mBound[0] = theCellSize * theDim[0];
	mBound[1] = theCellSize * (theDim[1] + 1);
	mBound[2] = theCellSize * theDim[2];
	mData.resize(theDim[0] * (theDim[1]+1) * theDim[2], false);
	std::fill(mData.begin(), mData.end(), mDflt);
}

double& GridDataY::operator()(int i, int j, int k) 
{
	// TODO : GridY accessor
    static double dflt = 0;
	dflt = mDflt;  

	if (j < 0 || j > theDim[1])
	{
		return dflt;
	}

	if (i < 0) i = 0;
	i = max(i, 0);
	i = min(i, theDim[0] - 1);
	k = max(k, 0);
	k = min(k, theDim[2] - 1);

	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * theDim[2];
	return mData[stack + row + col];
}

const double GridDataY::operator()(int i, int j, int k) const 
{
	// TODO : GridY accessor
	static double dflt = 0;
	dflt = mDflt;  

	if (j < 0 || j > theDim[1])
	{
		return dflt;
	}

	if (i < 0) i = 0;
	i = max(i, 0);
	i = min(i, theDim[0] - 1);
	k = max(k, 0);
	k = min(k, theDim[2] - 1);

	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * theDim[2];
	return mData[stack + row + col];
}

glm::dvec3 GridDataY::worldToSelf(const glm::dvec3& pt) const 
{
  // TODO : Given a point, returns the cell index that the grid uses in its own space
   glm::dvec3 res;
   res[0] = min(max(0.0, pt[0] - theCellSize * 0.5), mBound[0]);
   res[1] = min(max(0.0, pt[1]), mBound[1]);
   res[2] = min(max(0.0, pt[2] - theCellSize * 0.5), mBound[2]);
   return res;
}

GridDataZ::GridDataZ() : GridData() {
}

GridDataZ::~GridDataZ() {
}

void GridDataZ::initialize(double dfltValue) 
{
	// TODO : Intialize GridDataZ
	GridData::initialize(dfltValue);
	mBound[0] = theCellSize * theDim[0];
	mBound[1] = theCellSize * theDim[1];
	mBound[2] = theCellSize * (theDim[2] + 1);
	mData.resize(theDim[0] * theDim[1] * (theDim[2] + 1), false);
	std::fill(mData.begin(), mData.end(), mDflt);
}

double& GridDataZ::operator()(int i, int j, int k) 
{
	// TODO : GridZ accessor
	static double dflt = 0;
	dflt = mDflt;  
	if (k < 0 || k > theDim[2])
	{
		return dflt;
	}

	if (i < 0) i = 0;
	i = max(i, 0);
	i = min(i, theDim[0] - 1);
	j = max(j, 0);
	j = min(j, theDim[1] - 1);

	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * (theDim[2] + 1);
	return mData[stack + row + col];
}

const double GridDataZ::operator()(int i, int j, int k) const 
{
	// TODO : GridY accessor
	static double dflt = 0;
	dflt = mDflt;  
	if (k < 0 || k > theDim[2])
	{
		return dflt;
	}
	if (i < 0) i = 0;
	i = max(i, 0);
	i = min(i, theDim[0] - 1);
	j = max(j, 0);
	j = min(j, theDim[1] - 1);
	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * (theDim[2] + 1);
	return mData[stack + row + col];
}

glm::dvec3 GridDataZ::worldToSelf(const glm::dvec3& pt) const 
{
	// TODO : Given a point, returns the cell index that the grid uses in its own space
   glm::dvec3 res;
   res[0] = min(max(0.0, pt[0] - theCellSize * 0.5), mBound[0]);
   res[1] = min(max(0.0, pt[1] - theCellSize * 0.5), mBound[1]);
   res[2] = min(max(0.0, pt[2]), mBound[2]);
   return res;
}
