#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"

#include <math.h>
#include <map>
#include <stdio.h>

#undef max
#undef min
#include <fstream>

#define LEN 0.5
// Globals
MACGrid target;

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;//true

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; ++k)  \
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; ++k) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]+1; ++i) 

MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
	// TODO : Copy constructor for MAC Grid 
	mUf = orig.mUf;
	mVf = orig.mVf;
	mWf = orig.mWf;
	mUh = orig.mUh;
	mVh = orig.mVh;
	mWh = orig.mWh;
	mPressure = orig.mPressure;
	mDensityf = orig.mDensityf;
	mDensityh = orig.mDensityh;
	mTemp = orig.mTemp;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
	// TODO : Copy constructor for MAC Grid 
	if (&orig == this)
	{
      return *this;
	}
	mUf = orig.mUf;
	mVf = orig.mVf;
	mWf = orig.mWf;
	mUh = orig.mUh;
	mVh = orig.mVh;
	mWh = orig.mWh;
	mPressure = orig.mPressure;
	mDensityf = orig.mDensityf;
	mDensityh = orig.mDensityh;
	mTemp = orig.mTemp;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
	// TODO : Initialize the MAC Grid.
	mUf.initialize();
	mVf.initialize();
	mWf.initialize();
	mUh.initialize();
	mVh.initialize();
	mWh.initialize();
	mPressure.initialize(0);
	mDensityf.initialize();
	mDensityh.initialize();
	mTemp.initialize(0.0);
	mY.initialize();
	mN_x.initialize();
	mN_y.initialize();
	mN_z.initialize();
	mPhi.initialize();
	mPhi_State.initialize(0);
	initBoundary();
	findBoundaryGradient();
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
	// TODO: Set initial values for density, mTemp, velocity	
	FOR_EACH_CELL
	{
		//core
		if(mPhi(i,j,k) < 0)
		{
			/*
			if(mPhi(i,j,k) == 0)
			{
				mTemp(i,j,k) = Tmax; 
			}
			else
			{
				mTemp(i,j,k) = Tignition; 
			}
			*/
			mUf(i,j,k)=mN_x(i,j,k)*1;
			mVf(i,j,k)=mN_y(i,j,k)*1;
			mWf(i,j,k)=mN_z(i,j,k)*1;
			
			mTemp(i,j,k) = Tignition; 
			//mDensityf(i,j,k)=1;
			//mDensityf(i,j,k) = -1*mPhi(i,j,k)/radius; 
			mDensityf(i,j,k)=sqrt(-mPhi(i,j,k)/radius);
			mY(i,j,k)=1.0;
		}
		/*
		else if(mPhi(i,j,k) == 0)
		{
			mTemp(i,j,k) = Tmax; 
		}
		*/
	}
}

void MACGrid::advectVelocity(double dt,bool tag)
{
    // TODO: Calculate new velocities and store in target
	target.mUf = mUf;
    target.mVf = mVf;
    target.mWf = mWf;
	target.mUh = mUh;
	target.mVh = mVh;
	target.mWh = mWh;
	for(int k = 0; k < theDim[MACGrid::Z]; ++k)
	{
      for(int j = 0; j < theDim[MACGrid::Y]; ++j)
	  {
         for(int i = 0; i < theDim[MACGrid::X]+1; ++i)
		 {
			 glm::dvec3 posG = glm::dvec3(i*theCellSize, (j+LEN)*theCellSize,(k+LEN)*theCellSize);
			 glm::dvec3 velocityf=getVelocityf(posG);
			 glm::dvec3 velocityh=getVelocityh(posG);
			 glm::dvec3 normal=glm::dvec3(mN_x.interpolate(posG),mN_y.interpolate(posG),mN_z.interpolate(posG));
			 glm::dvec3 oldPosG;
			 if(glm::length(normal)<1e-5) normal=glm::dvec3(0,0,0);
			 else normal=glm::normalize(normal);
			 if(mPhi(i,j,k)>0) oldPosG=posG-velocityh*dt;
			 else oldPosG=posG-velocityf*dt;
			 double oldPhi=mPhi.interpolate(oldPosG);
			 float vNormal;
			 if(oldPhi<0)
			 {//inside
				 target.mUf(i,j,k)=mUf.interpolate(oldPosG);
				 //fuel
				 vNormal=glm::dot(velocityf,normal)+(fuelDensity/gasDensity-1.0)*S;
				 target.mUh(i,j,k)=vNormal*mN_x.interpolate(posG)+(velocityf-glm::dot(velocityf,normal)*normal).x;		
			 }
			 else
			 {//outside
				 target.mUh(i,j,k)=mUh.interpolate(oldPosG);
				 //hgp
				 vNormal=glm::dot(velocityh,normal)-(fuelDensity/gasDensity-1.0)*S;
				 target.mUf(i,j,k)=vNormal*mN_x.interpolate(posG)+(velocityh-glm::dot(velocityh,normal)*normal).x;
			 }
			// target.mUf(i,j,k)=0;
			// target.mUh(i,j,k)=0;
		 }
	  }
	}
	for(int k = 0; k < theDim[MACGrid::Z]; ++k) 
	{
      for(int j = 0; j < theDim[MACGrid::Y]+1; ++j) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 
		 {
			 glm::dvec3 posG = glm::dvec3((i+LEN)*theCellSize,j*theCellSize,(k+LEN)*theCellSize);
			 double newPhi=mPhi.interpolate(posG);
			 glm::dvec3 velocityf=getVelocityf(posG);
			 glm::dvec3 velocityh=getVelocityh(posG);
			 glm::dvec3 normal=glm::dvec3(mN_x.interpolate(posG),mN_y.interpolate(posG),mN_z.interpolate(posG));
			 glm::dvec3 oldPosG;
			 if(glm::length(normal)<1e-5) normal=glm::dvec3(0,0,0);
			 else normal=glm::normalize(normal);
			 if(mPhi(i,j,k)>0) oldPosG=posG-velocityh*dt;
			 else oldPosG=posG-velocityf*dt;
			 double oldPhi=mPhi.interpolate(oldPosG);
			 float vNormal;
			 if(oldPhi<0){//inside
				 target.mVf(i,j,k)=mVf.interpolate(oldPosG);
				 //fuel
				 vNormal=glm::dot(velocityf,normal)+(fuelDensity/gasDensity-1.0)*S;
				 target.mVh(i,j,k)=vNormal*mN_y.interpolate(posG)+(velocityf-glm::dot(velocityf,normal)*normal).y;
			 }
			 else{//outside
				 target.mVh(i,j,k)=mVh.interpolate(oldPosG);
				 //hgp
				 vNormal=glm::dot(velocityh,normal)-(fuelDensity/gasDensity-1.0)*S;
				 target.mVf(i,j,k)=vNormal*mN_y.interpolate(posG)+(velocityh-glm::dot(velocityh,normal)*normal).y;
			 }
			 //target.mVf(i,j,k)=0;
			// target.mVh(i,j,k)=0;
		 }
	  }
	}

	for(int k = 0; k < theDim[MACGrid::Z]+1; ++k) 
	{
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 
		 {
			 glm::dvec3 posG = glm::dvec3((i+LEN)*theCellSize, (j+LEN)*theCellSize,k*theCellSize);
			 double newPhi=mPhi.interpolate(posG);
			 glm::dvec3 velocityf=getVelocityf(posG);
			 glm::dvec3 velocityh=getVelocityh(posG);
			 glm::dvec3 normal=glm::dvec3(mN_x.interpolate(posG),mN_y.interpolate(posG),mN_z.interpolate(posG));
			 glm::dvec3 oldPosG;
			 if(glm::length(normal)<1e-5) normal=glm::dvec3(0,0,0);
			 else normal=glm::normalize(normal);
			 if(mPhi(i,j,k)>0) oldPosG=posG-velocityh*dt;
			 else oldPosG=posG-velocityf*dt;
			 double oldPhi=mPhi.interpolate(oldPosG);
			 float vNormal;
			 if(oldPhi<0){//inside
				 target.mWf(i,j,k)=mWf.interpolate(oldPosG);
				 //fuel
				 vNormal=glm::dot(velocityf,normal)+(fuelDensity/gasDensity-1.0)*S;
				 target.mWh(i,j,k)=vNormal*mN_z.interpolate(posG)+(velocityf-glm::dot(velocityf,normal)*normal).z;
			 }
			 else{//outside
				 target.mWh(i,j,k)=mWh.interpolate(oldPosG);
				 //hgp
				 vNormal=glm::dot(velocityh,normal)-(fuelDensity/gasDensity-1)*S;
				 target.mWf(i,j,k)=vNormal*mN_z.interpolate(posG)+(velocityh-glm::dot(velocityh,normal)*normal).z;
			 }
			 //target.mWf(i,j,k)=0;
			 //target.mWh(i,j,k)=0;
		 }
	  }
	}
    mUf = target.mUf;
    mVf = target.mVf;
    mWf = target.mWf;
	mUh = target.mUh;
	mVh = target.mVh;
	mWh = target.mWh;
}

double MACGrid::getUgradientH(int i, int j, int k)
{
	double udiv=0.0;
	if(i!=0){
		if(mPhi(i,j,k)<0) udiv-=mUf(i,j,k);
		else udiv-=mUh(i,j,k);
	}
	if(j!=0){
		if(mPhi(i,j,k)<0) udiv-=mVf(i,j,k);
		else udiv-=mVh(i,j,k);
	}
	if(k!=0){
		if(mPhi(i,j,k)<0) udiv-=mWf(i,j,k);
		else udiv-=mWh(i,j,k);
	}
	if(i!=theDim[0]-1){
		if(mPhi(i+1,j,k)<0) udiv+=mUf(i+1,j,k);
		else udiv+=mUh(i+1,j,k);
	}
	if(j!=theDim[1]-1){
		if(mPhi(i,j+1,k)<0) udiv+=mVf(i,j+1,k);
		else udiv+=mVh(i,j+1,k);
	}
	if(k!=theDim[2]-1){
		if(mPhi(i,j,k+1)<0) udiv+=mWf(i,j,k+1);
		else udiv+=mWh(i,j,k+1);
	}
	return udiv;
}

void MACGrid::advectTemperature(double dt,bool tag)
{
	target.mTemp = mTemp;
	 FOR_EACH_CELL
	  {
			glm::dvec3 pos = getCenter(i,j,k);
			double T;
			double Y = mY(i,j,k);
			//make temp rise 
			if(Y >= Y_threshold)
			{
				T = Tmax - (Y - Y_threshold)/ (1.0 - Y_threshold) * Trise;
			}
			//temp is cooling
			else  
			{
				glm::dvec3 vel = getVelocityh(pos);
				glm::dvec3 oldPos = pos - vel * dt;	
				double oldTemp = mTemp.interpolate(oldPos); 
				double tmp = (oldTemp - Tair) / (Tmax - Tair); 
				double udiv = getUgradientH(i,j,k); 
				T = oldTemp + (-udiv * oldTemp - CoolT * pow(tmp, 4)) * dt; 
				T = std::min(T, Tmax);
				//T = std::max(T, 0.0);
				//printf("%lf ", T);
				//T = std::max(0.0, Tmax - (1.0 - Y) * Tmax);
				//printf("dis:%lf ", oldTemp -T);
			} 
			/*
			if(mPhi(i,j,k) < 0)
			{
				T = Tignition;
			}
			else
			{
				glm::dvec3 vel = getVelocityh(pos);
				glm::dvec3 oldPos = pos - vel * dt;
				double oldPhi = mPhi.interpolate(oldPos);
				//printf("o:%lf ", oldPhi);
				T = std::max(Tignition, Tmax - oldPhi); 
			}
			*/
			target.mTemp(i,j,k)=T;
			//debug
			//target.mTemp(i,j,k) = mTemp.interpolate(oldPos); 
	  }
	 mTemp = target.mTemp;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target
	target.mDensityf = mDensityf;
	//target.mDensityh = mDensityh;
	FOR_EACH_CELL
	{
		glm::dvec3 posG = getCenter(i,j,k);
		if(mPhi(i,j,k) < 0)
		{
			target.mDensityf(i,j,k) = mDensityf.interpolate(posG - getVelocityf(posG)*dt);
		}
		else
		{
			target.mDensityf(i,j,k) = mDensityf.interpolate(posG - getVelocityh(posG)*dt);
			//temp falling
			if(Y < Y_threshold)
			{
				double old = target.mDensityf(i,j,k); 
				target.mDensityf(i,j,k) = old + (-getUgradientH(i,j,k) * old) * dt; 
			}
		}	
		//if(target.mDensityf(i,j,k)==0) n++;
	}
    mDensityf = target.mDensityf;
	//mDensityh = target.mDensityh;
}

void MACGrid::computeBouyancy(double dt)
{
  	// TODO: Calculate bouyancy and store in target
	//target.mVf = mVf;
	target.mVh = mVh;
	for(int k = 0; k < theDim[MACGrid::Z]; ++k) 
	{
      for(int j = 1; j < theDim[MACGrid::Y]; ++j) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 
		 {
			 double s = (mDensityf(i,j,k) + mDensityf(i,j-1,k))/2;
			 double temp = (mTemp(i,j,k)+mTemp(i,j-1,k))/2;
			 double f = -theBuoyancyAlpha * s + theBuoyancyBeta * temp;
			 target.mVh(i,j,k) += dt * f;
			 //target.mVf(i,j,k) += dt * f;
			 //if(s!=0) std::cout<<s<<std::endl;
		 }
	  }
	}
    //mVf = target.mVf;
	mVh = target.mVh;
}

void MACGrid::computeVorticityConfinement(double dt)
{
   // TODO: Calculate vorticity confinement forces
   // Apply the forces to the current velocity and store the result in target
	target.mUf = mUf;
	target.mVf = mVf;
	target.mWf = mWf;
	target.mUh = mUh;
	target.mVh = mVh;
	target.mWh = mWh;
	GridData velu,velv,velw,wu,wv,ww,wlen,fu,fv,fw;
	velu.initialize();
	velv.initialize();
	velw.initialize();
	wu.initialize();
	wv.initialize();
	ww.initialize();
	wlen.initialize();
	fu.initialize();
	fv.initialize();
	fw.initialize();

	FOR_EACH_CELL 
	{
		glm::dvec3 pos = getCenter(i,j,k);
		glm::dvec3 vel = getVelocityh(pos);
		velu(i,j,k) = vel[0];
		velv(i,j,k) = vel[1];
		velw(i,j,k) = vel[2];
    }

	FOR_EACH_CELL 
	{
		double u,v,w; 
		u = v = w = 0;
		if(j!=theDim[1]) 
		{ 
			u += velw(i,j+1,k); 
			w -= velu(i,j+1,k);
		}
		if(j!=0) 
		{ 
			u -= velw(i,j-1,k); 
			w += velu(i,j-1,k);
		}
		if(k!=theDim[2]) 
		{ 
			u -= velv(i,j,k+1); 
			v += velu(i,j,k+1);
		}
		if(k!=0) 
		{ 
			u += velu(i,j,k-1); 
			v -= velu(i,j,k-1);
		}
		if(i!=theDim[0]) 
		{ 
			v -= velw(i+1,j,k); 
			w += velv(i+1,j,k);
		}
		if(i!=0) 
		{ 
			v += velw(i-1,j,k); 
			w -= velv(i-1,j,k);
		}
		wu(i,j,k) = u/(2*theCellSize);
		wv(i,j,k) = v/(2*theCellSize);
		ww(i,j,k) = w/(2*theCellSize);
		glm::dvec3 wvec = glm::dvec3(wu(i,j,k),wv(i,j,k),ww(i,j,k));
		wlen(i,j,k) = glm::length(wvec); 
	}

	FOR_EACH_CELL 
	{
		double u,v,w; 
		u = v = w = 0;
		if(j!=theDim[1]) 
		{ 
			v += wlen(i,j+1,k);
		}
		if(j!=0) 
		{
			v -= wlen(i,j-1,k);
		}
		if(k!=theDim[2]) 
		{
			w += wlen(i,j,k+1);
		}
		if(k!=0) 
		{
			w -= wlen(i,j,k-1);
		}
		if(i!=theDim[0]) 
		{
			u += wlen(i+1,j,k);
		}
		if(i!=0)
		{
			u -= wlen(i-1,j,k);
		}
		double wgradlen=0.0;
		glm::dvec3 wgradvec = glm::dvec3(u,v,w)/(2*theCellSize);
		wgradlen = glm::length(wgradvec);
		glm::dvec3 N = wgradvec/(wgradlen + 1e-20);
		glm::dvec3 wvec = glm::dvec3(wu(i,j,k),wv(i,j,k),ww(i,j,k));
		double ff; 
		if(mPhi(i,j,k) < 0)
		{
			ff = 5.0;
		}
		else
		{
			ff = 2.0;
		}
		glm::dvec3 fvec = ff * theCellSize *(glm::cross(N, wvec));
		fu(i,j,k) = fvec[0];
		fv(i,j,k) = fvec[1];
		fw(i,j,k) = fvec[2];
	}
	for(int k = 0; k < theDim[MACGrid::Z]; ++k) 
	{
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) 
	  {
         for(int i = 1; i < theDim[MACGrid::X]; ++i) 
		 {
			 target.mUf(i,j,k) += dt*(fu(i,j,k)+fu(i-1,j,k))/2;
			 target.mUh(i,j,k) += dt*(fu(i,j,k)+fu(i-1,j,k))/2;
		 }
	  }
	}
	for(int k = 0; k < theDim[MACGrid::Z]; ++k) 
	{
      for(int j = 1; j < theDim[MACGrid::Y]; ++j) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 
		 {
			 target.mVf(i,j,k) += dt*(fv(i,j,k)+fu(i,j-1,k))/2;
			 target.mVh(i,j,k) += dt*(fv(i,j,k)+fu(i,j-1,k))/2;
		 }
	  }
	}
	for(int k = 1; k < theDim[MACGrid::Z]; ++k) 
	{
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 
		 {
			 target.mWf(i,j,k) += dt*(fw(i,j,k)+fu(i,j,k-1))/2;
			 target.mWh(i,j,k) += dt*(fw(i,j,k)+fu(i,j,k-1))/2;
		 }
	  }
	}
	mUf = target.mUf;
	mVf = target.mVf;
	mWf = target.mWf;
	mUh = target.mUh;
	mVh = target.mVh;
	mWh = target.mWh;
}

bool MACGrid::isOnBorder(int i,int j,int k){
	int sign1,sign2;
	if(mPhi(i,j,k)>0) sign1=1;
	else sign1=-1;

	if(sign1==1){
		if(i!=0&&mPhi(i-1,j,k)<=0) return true;
		if(i!=theDim[0]-1&&mPhi(i+1,j,k)<=0) return true;
		if(j!=0&&mPhi(i,j-1,k)<=0) return true;
		if(j!=theDim[1]-1&&mPhi(i,j+1,k)<=0) return true;
		if(k!=0&&mPhi(i,j,k-1)<=0) return true;
		if(k!=theDim[2]-1&&mPhi(i,j,k+1)<=0) return true;
	}
	else{
		if(i!=0&&mPhi(i-1,j,k)>0) return true;
		if(i!=theDim[0]-1&&mPhi(i+1,j,k)>0) return true;
		if(j!=0&&mPhi(i,j-1,k)>0) return true;
		if(j!=theDim[1]-1&&mPhi(i,j+1,k)>0) return true;
		if(k!=0&&mPhi(i,j,k-1)>0) return true;
		if(k!=theDim[2]-1&&mPhi(i,j,k+1)>0) return true;
	}
	return false;
}

void MACGrid::addExternalForces(double dt,bool tag)
{
	computeBouyancy(dt);
	computeVorticityConfinement(dt);
}

void MACGrid::project(double dt)
{
   // TODO: Solve Ax = b for pressure
   // 1. Contruct b
   // 2. Construct A 
   // 3. Solve for p
   // Subtract pressure from our velocity and save in target

	target.mPressure=mPressure;
	target.mUf = mUf;
	target.mVf = mVf;
	target.mWf = mWf;

	target.mUh = mUh;
	target.mVh = mVh;
	target.mWh = mWh;
	
	calculateAMatrix();
	GridData d;
	d.initialize();
	
	FOR_EACH_CELL 
	{
		double udiv=0.0;
		if(i!=0){
			if(mPhi(i,j,k)<0) udiv-=mUf(i,j,k);
			else udiv-=mUh(i,j,k);
		}
		if(j!=0){
			if(mPhi(i,j,k)<0) udiv-=mVf(i,j,k);
			else udiv-=mVh(i,j,k);
		}
		if(k!=0){
			if(mPhi(i,j,k)<0) udiv-=mWf(i,j,k);
			else udiv-=mWh(i,j,k);
		}
		/*
		if(mPhi(i,j,k)<0) udiv+=mUf(i,j,k);
		else udiv+=mUh(i,j,k);
		if(mPhi(i,j,k)<0) udiv+=mVf(i,j,k);
		else udiv+=mVh(i,j,k);
		if(mPhi(i,j,k)<0) udiv+=mWf(i,j,k);
		else udiv+=mWh(i,j,k);
		*/
		if(i!=theDim[0]-1){
			if(mPhi(i+1,j,k)<0) udiv+=mUf(i+1,j,k);
			else udiv+=mUh(i+1,j,k);
		}
		if(j!=theDim[1]-1){
			if(mPhi(i,j+1,k)<0) udiv+=mVf(i,j+1,k);
			else udiv+=mVh(i,j+1,k);
		}
		if(k!=theDim[2]-1){
			if(mPhi(i,j,k+1)<0) udiv+=mWf(i,j,k+1);
			else udiv+=mWh(i,j,k+1);
		}
		
		if(mPhi(i,j,k)>0) d(i,j,k)=-theCellSize/dt * udiv*gasDensity;
		else d(i,j,k)=-theCellSize/dt * udiv*gasDensity;

		/*if(d(i,j,k)!=0){
			std::cout<<udiv<<","<<i<<","<<j<<","<<k<<std::endl;
			getchar();
		}*/
	}

	conjugateGradient(AMatrix,mPressure,d,3000,0.01);
	
	float C=S*S*fuelDensity*(fuelDensity/gasDensity-1);

	for(int k = 0; k < theDim[MACGrid::Z]; ++k) 
	{
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) 
	  {
         for(int i = 1; i < theDim[MACGrid::X]; ++i) 
		 {
			 //target.mUf(i,j,k) -= dt*(target.mPressure(i,j,k)-target.mPressure(i-1,j,k))/theCellSize;
			 float p1,p2,p3,p4;
			 if(isOnBorder(i,j,k)&&isOnBorder(i-1,j,k)){
				 p1=target.mPressure(i,j,k);
				 if(mPhi(i,j,k)>0){
					 p3=target.mPressure(i,j,k)+C;
					 if(mPhi(i-1,j,k)>0){
						 p2=target.mPressure(i-1,j,k);
						 p4=target.mPressure(i-1,j,k)+C;
					 }
					 else{
						 p2=target.mPressure(i-1,j,k)-C;
						 p4=target.mPressure(i-1,j,k);
					 }
					 target.mUh(i,j,k)-=dt*(p1-p2)/theCellSize;
					 target.mUf(i,j,k)-=dt*(p3-p4)/theCellSize;
				 }
				 else{
					 p3=target.mPressure(i,j,k)-C;
					 if(mPhi(i-1,j,k)>0){
						 p2=target.mPressure(i-1,j,k)+C;
						 p4=target.mPressure(i-1,j,k);
					 }
					 else{
						 p2=target.mPressure(i-1,j,k);
						 p4=target.mPressure(i-1,j,k)-C;
					 }
					 target.mUf(i,j,k)-=dt*(p1-p2)/theCellSize;
					 target.mUh(i,j,k)-=dt*(p3-p4)/theCellSize;
				 }
			 }
			 else{
				 p1=target.mPressure(i,j,k);
				 p2=target.mPressure(i-1,j,k);
				 if(mPhi(i,j,k)>0) target.mUh(i,j,k)-=dt*(p1-p2)/theCellSize;
				 else target.mUf(i,j,k)-=dt*(p1-p2)/theCellSize;
			 }
		 }
	  }
	}
	for(int k = 0; k < theDim[MACGrid::Z]; ++k) 
	{
      for(int j = 1; j < theDim[MACGrid::Y]; ++j) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 
		 {
			 //target.mVf(i,j,k) -= dt*(target.mPressure(i,j,k)-target.mPressure(i,j-1,k))/theCellSize;
			float p1,p2,p3,p4;
			 if(isOnBorder(i,j,k)&&isOnBorder(i,j-1,k)){
				 p1=target.mPressure(i,j,k);
				 if(mPhi(i,j,k)>0){
					 p3=target.mPressure(i,j,k)+C;
					 if(mPhi(i,j-1,k)>0){
						 p2=target.mPressure(i,j-1,k);
						 p4=target.mPressure(i,j-1,k)+C;
					 }
					 else{
						 p2=target.mPressure(i,j-1,k)-C;
						 p4=target.mPressure(i,j-1,k);
					 }
					 target.mVh(i,j,k)-=dt*(p1-p2)/theCellSize;
					 target.mVf(i,j,k)-=dt*(p3-p4)/theCellSize;
				 }
				 else{
					 p3=target.mPressure(i,j,k)-C;
					 if(mPhi(i-1,j,k)>0){
						 p2=target.mPressure(i,j-1,k)+C;
						 p4=target.mPressure(i,j-1,k);
					 }
					 else{
						 p2=target.mPressure(i,j-1,k);
						 p4=target.mPressure(i,j-1,k)-C;
					 }
					 target.mVf(i,j,k)-=dt*(p1-p2)/theCellSize;
					 target.mVh(i,j,k)-=dt*(p3-p4)/theCellSize;
				 }
			 }
			 else{
				 p1=target.mPressure(i,j,k);
				 p2=target.mPressure(i,j-1,k);
				 if(mPhi(i,j,k)>0) target.mVh(i,j,k)-=dt*(p1-p2)/theCellSize;
				 else target.mVf(i,j,k)-=dt*(p1-p2)/theCellSize;
			 }
		 }
	  }
	}
	for(int k = 1; k < theDim[MACGrid::Z]; ++k) 
	{
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 
		 {
			 //target.mWf(i,j,k) -= dt*(target.mPressure(i,j,k)-target.mPressure(i,j,k-1))/theCellSize;
			float p1,p2,p3,p4;
			 if(isOnBorder(i,j,k)&&isOnBorder(i,j,k-1)){
				 p1=target.mPressure(i,j,k);
				 if(mPhi(i,j,k)>0){
					 p3=target.mPressure(i,j,k)+C;
					 if(mPhi(i,j,k-1)>0){
						 p2=target.mPressure(i,j,k-1);
						 p4=target.mPressure(i,j,k-1)+C;
					 }
					 else{
						 p2=target.mPressure(i,j,k-1)-C;
						 p4=target.mPressure(i,j,k-1);
					 }
					 target.mWh(i,j,k)-=dt*(p1-p2)/theCellSize;
					 target.mWf(i,j,k)-=dt*(p3-p4)/theCellSize;
				 }
				 else{
					 p3=target.mPressure(i,j,k)-C;
					 if(mPhi(i,j,k-1)>0){
						 p2=target.mPressure(i,j,k-1)+C;
						 p4=target.mPressure(i,j,k-1);
					 }
					 else{
						 p2=target.mPressure(i,j,k-1);
						 p4=target.mPressure(i,j,k-1)-C;
					 }
					 target.mWf(i,j,k)-=dt*(p1-p2)/theCellSize;
					 target.mWh(i,j,k)-=dt*(p3-p4)/theCellSize;
				 }
			 }
			 else{
				 p1=target.mPressure(i,j,k);
				 p2=target.mPressure(i,j,k-1);
				 if(mPhi(i,j,k)>0) target.mWh(i,j,k)-=dt*(p1-p2)/theCellSize;
				 else target.mWf(i,j,k)-=dt*(p1-p2)/theCellSize;
			 }
		 }
	  }
	}
	mPressure = target.mPressure;
	mUf = target.mUf;
	mVf = target.mVf;
	mWf = target.mWf;

	mUh = target.mUh;
	mVh = target.mVh;
	mWh = target.mWh;

	// IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
}

glm::dvec3 MACGrid::getVelocityf(const glm::dvec3& pt) 
{
	// TODO : Given a point in space, give the 3D velocity field at the point
	glm::dvec3 v;
	v[0] = mUf.interpolate(pt); 
	v[1] = mVf.interpolate(pt); 
	v[2] = mWf.interpolate(pt); 
	return v;
}

glm::dvec3 MACGrid::getVelocityh(const glm::dvec3& pt) 
{
	// TODO : Given a point in space, give the 3D velocity field at the point
	glm::dvec3 v;
	v[0] = mUh.interpolate(pt); 
	v[1] = mVh.interpolate(pt); 
	v[2] = mWh.interpolate(pt); 
	return v;
}

glm::dvec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return glm::dvec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::calculateAMatrix() 
{
	FOR_EACH_CELL 
	{
		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

void MACGrid::calculateAMatrixOut() {

	FOR_EACH_CELL 
	{
		if(mPhi(i,j,k)<=0) AMatrix.diag(i,j,k)=1;
		else{
			int numFluidNeighbors = 0;
			if (i-1 >= 0) {
				AMatrix.plusI(i-1,j,k) = -1;
				numFluidNeighbors++;
			}
			if (i+1 < theDim[MACGrid::X]) {
				AMatrix.plusI(i,j,k) = -1;
				numFluidNeighbors++;
			}
			if (j-1 >= 0) {
				AMatrix.plusJ(i,j-1,k) = -1;
				numFluidNeighbors++;
			}
			if (j+1 < theDim[MACGrid::Y]) {
				AMatrix.plusJ(i,j,k) = -1;
				numFluidNeighbors++;
			}
			if (k-1 >= 0) {
				AMatrix.plusK(i,j,k-1) = -1;
				numFluidNeighbors++;
			}
			if (k+1 < theDim[MACGrid::Z]) {
				AMatrix.plusK(i,j,k) = -1;
				numFluidNeighbors++;
			}
			// Set the diagonal:
			AMatrix.diag(i,j,k) = numFluidNeighbors;
		}
	}
}


bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL 
	{
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.
	GridData z; 
	z.initialize();
	// TODO : Apply preconditioner; for now, bypass the preconditioner
	z = r;
	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) 
	{
		double rho = sigma; // According to Aline. Here???

		apply(A, s, z); // z = applyA(s);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; 
		alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; 
		alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) 
		{
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true; //return p;
		}

    // TODO : Apply preconditioner; for now, bypass the preconditioner
		z = r;		
		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	std::cout<<maxMagnitude(r)<<std::endl;
	return false;

}

//fire simUflation

void MACGrid::initBoundary()
{
	int n,m;
	n = m = 0;
	FOR_EACH_CELL
	{
		glm::dvec3 dis=glm::dvec3(i,j,k)-fireCenter;
		//float distance=theCellSize*glm::length(dis)-radius;
		double distance=std::sqrt(dis.x*dis.x/(a*a)+dis.y*dis.y/(b*b)+dis.z*dis.z/(c*c))*radius-radius;
		mPhi(i,j,k) = distance;
		mTemp(i,j,k) = Tignition;	

		if(mPhi(i,j,k)<0){
			//mDensityf(i,j,k)=sqrt(-mPhi(i,j,k)/radius);
			//mDensityf(i,j,k) = -1*mPhi(i,j,k)/radius; 
			mTemp(i,j,k)=Tignition;
		}
	}
}

void MACGrid::findBoundaryGradient(){
	double dx,dy,dz,length,v1,v2;
	FOR_EACH_CELL{
		if(i==theDim[0]-1){
			dx=mPhi(i,j,k)-mPhi(i-1,j,k);
		}
		else if(i==0){
			dx=mPhi(i+1,j,k)-mPhi(i,j,k);
		}
		else{
			dx=(mPhi(i+1,j,k)-mPhi(i-1,j,k))/2;
		}
		if(j==theDim[1]-1){
			dy=mPhi(i,j,k)-mPhi(i,j-1,k);
		}
		else if(j==0){
			dy=mPhi(i,j+1,k)-mPhi(i,j,k);
		}
		else{
			dy=(mPhi(i,j+1,k)-mPhi(i,j-1,k))/2;
		}

		if(k==theDim[2]-1){
			dz=mPhi(i,j,k)-mPhi(i,j,k-1);
		}
		else if(k==0){
			dz=mPhi(i,j,k+1)-mPhi(i,j,k);
		}
		else{
			dz=(mPhi(i,j,k+1)-mPhi(i,j,k-1))/2;
		}

		length=std::sqrt(dx*dx+dy*dy+dz*dz);
		if(length>1e-5){
			mN_x(i,j,k)=dx/length;
			mN_y(i,j,k)=dy/length;
			mN_z(i,j,k)=dz/length;
		}
		else{
			mN_x(i,j,k)=0;
			mN_y(i,j,k)=0;
			mN_z(i,j,k)=0;
		}
	}//FOR_EACH_CELL
}

void MACGrid::advectBoundary(double dt,bool tag){
	FOR_EACH_CELL{
		glm::dvec3 w,uf,n;
		uf=getVelocityf(getCenter(i,j,k));
		n=glm::dvec3(mN_x.interpolate(getCenter(i,j,k)),mN_y.interpolate(getCenter(i,j,k)),mN_z.interpolate(getCenter(i,j,k)));
		w=uf+S*n;

		mPhi(i,j,k)=mPhi(i,j,k)-dt*(w.x*mN_x(i,j,k)+w.y*mN_y(i,j,k)+w.z*mN_z(i,j,k));
	}
}

void MACGrid::findBoundary(bool tag){//find the new boundary
	mPhi_State.initialize(0);
	mNeighborX.initialize(-1);
	mNeighborY.initialize(-1);
	mNeighborZ.initialize(-1);
	unTagedNum=theDim[0]*theDim[1]*theDim[2];
	FOR_EACH_CELL{
		if(i!=0){
			/*if(tag){
				std::cout<<mPhi(i,j,k)<<std::endl;
				getchar();
			}*/
			/*if(mPhi(i,j,k)<0){
				std::cout<<mPhi(i,j,k)<<std::endl;
				getchar();
			}*/
			if(mPhi(i,j,k)*mPhi(i-1,j,k)<=0){
				if(mPhi_State(i,j,k)==0){
					mPhi_State(i,j,k)=1;
					mNeighborX(i,j,k)=i;
					mNeighborY(i,j,k)=j;
					mNeighborZ(i,j,k)=k;
					unTagedNum--;
				}
				if(mPhi_State(i-1,j,k)==0){
					mPhi_State(i-1,j,k)=1;
					mNeighborX(i-1,j,k)=i-1;
					mNeighborY(i,j,k)=j;
					mNeighborZ(i,j,k)=k;
					unTagedNum--;
				}
			}
		}
		if(j!=0){
			if(mPhi(i,j,k)*mPhi(i,j-1,k)<=0){
				if(mPhi_State(i,j,k)==0){
					mPhi_State(i,j,k)=1;
					mNeighborX(i,j,k)=i;
					mNeighborY(i,j,k)=j;
					mNeighborZ(i,j,k)=k;
					unTagedNum--;
				}
				if(mPhi_State(i,j-1,k)==0){
					mPhi_State(i,j-1,k)=1;
					mNeighborX(i,j,k)=i;
					mNeighborY(i,j-1,k)=j-1;
					mNeighborZ(i,j,k)=k;
					unTagedNum--;
				}
			}
		}
		if(k!=0){
			if(mPhi(i,j,k)*mPhi(i,j,k-1)<=0){
				if(mPhi_State(i,j,k)==0){
					mPhi_State(i,j,k)=1;
					mNeighborX(i,j,k)=i;
					mNeighborY(i,j,k)=j;
					mNeighborZ(i,j,k)=k;
					unTagedNum--;
				}
				if(mPhi_State(i,j,k-1)==0){
					mPhi_State(i,j,k-1)=1;
					mNeighborX(i,j,k)=i;
					mNeighborY(i,j,k)=j;
					mNeighborZ(i,j,k-1)=k-1;
					unTagedNum--;
				}
			}
		}
	}
	std::cout<<unTagedNum<<std::endl;
}

void MACGrid::insideBoundary(){
	mInside.initialize(0);
	FOR_EACH_CELL{
		int posCount=0,negCount=0;
		if(i!=0){
			if(mPhi(i-1,j,k)>=0) posCount++;
			else negCount++;
		}
		if(i!=theDim[0]-1){
			if(mPhi(i+1,j,k)>=0) posCount++;
			else negCount++;
		}
		if(j!=0){
			if(mPhi(i,j-1,k)>=0) posCount++;
			else negCount++;
		}
		if(j!=theDim[1]-1){
			if(mPhi(i,j+1,k)>=0) posCount++;
			else negCount++;
		}
		if(k!=0){
			if(mPhi(i,j,k-1)>=0) posCount++;
			else negCount++;
		}
		if(k!=theDim[2]-1){
			if(mPhi(i,j,k+1)>=0) posCount++;
			else negCount++;
		}
		if(posCount>=negCount) mInside(i,j,k)=1;
		else mInside(i,j,k)=0;
	}//FOR EACH CELL
}

void MACGrid::reconditioning(){
	if(unTagedNum==theDim[0]*theDim[1]*theDim[2]){
		std::cout<<"!"<<std::endl;
		getchar();
	}
	while(unTagedNum>0&&unTagedNum!=theDim[0]*theDim[1]*theDim[2]){
		FOR_EACH_CELL{
			if(mPhi_State(i,j,k)==0){
				float minDis=1e7;
				if(i!=0&&mPhi_State(i-1,j,k)==1){
					glm::dvec3 target=glm::dvec3(mNeighborX(i-1,j,k),mNeighborY(i-1,j,k),mNeighborZ(i-1,j,k));
					glm::dvec3 pos=glm::dvec3(i,j,k);
					float dis=theCellSize*glm::length(target-pos);
					if(fabs(dis)<fabs(mPhi(i-1,j,k))){
						mPhi_State(i-1,j,k)=0;
						unTagedNum++;
					}
					mPhi_State(i,j,k)=1;
					if(fabs(dis)<fabs(minDis)){
						minDis=dis;
						mNeighborX(i,j,k)=mNeighborX(i-1,j,k);
						mNeighborY(i,j,k)=mNeighborY(i-1,j,k);
						mNeighborZ(i,j,k)=mNeighborZ(i-1,j,k);
					}
				}
				if(i!=theDim[0]-1&&mPhi_State(i+1,j,k)==1){
					glm::dvec3 target=glm::dvec3(mNeighborX(i+1,j,k),mNeighborY(i+1,j,k),mNeighborZ(i+1,j,k));
					glm::dvec3 pos=glm::dvec3(i,j,k);
					float dis=theCellSize*glm::length(target-pos);
					if(fabs(dis)<fabs(mPhi(i+1,j,k))){
						mPhi_State(i+1,j,k)=0;
						unTagedNum++;
					}
					mPhi_State(i,j,k)=1;
					if(fabs(dis)<fabs(minDis)){
						minDis=dis;
						mNeighborX(i,j,k)=mNeighborX(i+1,j,k);
						mNeighborY(i,j,k)=mNeighborY(i+1,j,k);
						mNeighborZ(i,j,k)=mNeighborZ(i+1,j,k);
					}
				}
				if(j!=0&&mPhi_State(i,j-1,k)==1){
					glm::dvec3 target=glm::dvec3(mNeighborX(i,j-1,k),mNeighborY(i,j-1,k),mNeighborZ(i,j-1,k));
					glm::dvec3 pos=glm::dvec3(i,j,k);
					float dis=theCellSize*glm::length(target-pos);
					if(fabs(dis)<fabs(mPhi(i,j-1,k))){
						mPhi_State(i,j-1,k)=0;
						unTagedNum++;
					}
					mPhi_State(i,j,k)=1;
					if(fabs(dis)<fabs(minDis)){
						minDis=dis;
						mNeighborX(i,j,k)=mNeighborX(i,j-1,k);
						mNeighborY(i,j,k)=mNeighborY(i,j-1,k);
						mNeighborZ(i,j,k)=mNeighborZ(i,j-1,k);
					}
				}
				if(j!=theDim[1]-1&&mPhi_State(i,j+1,k)==1){
					glm::dvec3 target=glm::dvec3(mNeighborX(i,j+1,k),mNeighborY(i,j+1,k),mNeighborZ(i,j+1,k));
					glm::dvec3 pos=glm::dvec3(i,j,k);
					float dis=theCellSize*glm::length(target-pos);
					if(fabs(dis)<fabs(mPhi(i,j+1,k))){
						mPhi_State(i,j+1,k)=0;
						unTagedNum++;
					}
					mPhi_State(i,j,k)=1;
					if(fabs(dis)<fabs(minDis)){
						minDis=dis;
						mNeighborX(i,j,k)=mNeighborX(i,j+1,k);
						mNeighborY(i,j,k)=mNeighborY(i,j+1,k);
						mNeighborZ(i,j,k)=mNeighborZ(i,j+1,k);
					}
				}
				if(k!=0&&mPhi_State(i,j,k-1)==1){
					glm::dvec3 target=glm::dvec3(mNeighborX(i,j,k-1),mNeighborY(i,j,k-1),mNeighborZ(i,j,k-1));
					glm::dvec3 pos=glm::dvec3(i,j,k);
					float dis=theCellSize*glm::length(target-pos);
					if(fabs(dis)<fabs(mPhi(i,j,k-1))){
						mPhi_State(i,j,k-1)=0;
						unTagedNum++;
					}
					mPhi_State(i,j,k)=1;
					if(fabs(dis)<fabs(minDis)){
						minDis=dis;
						mNeighborX(i,j,k)=mNeighborX(i,j,k-1);
						mNeighborY(i,j,k)=mNeighborY(i,j,k-1);
						mNeighborZ(i,j,k)=mNeighborZ(i,j,k-1);
					}
				}
				if(k!=theDim[2]-1&&mPhi_State(i,j,k+1)==1){
					glm::dvec3 target=glm::dvec3(mNeighborX(i,j,k+1),mNeighborY(i,j,k+1),mNeighborZ(i,j,k+1));
					glm::dvec3 pos=glm::dvec3(i,j,k);
					float dis=theCellSize*glm::length(target-pos);
					if(fabs(dis)<fabs(mPhi(i,j,k+1))){
						mPhi_State(i,j,k+1)=0;
						unTagedNum++;
					}
					mPhi_State(i,j,k)=1;
					if(fabs(dis)<fabs(minDis)){
						minDis=dis;
						mNeighborX(i,j,k)=mNeighborX(i,j,k+1);
						mNeighborY(i,j,k)=mNeighborY(i,j,k+1);
						mNeighborZ(i,j,k)=mNeighborZ(i,j,k+1);
					}
				}
				if(mPhi_State(i,j,k)==1){
					unTagedNum--;
					if(mInside(i,j,k)==1) mPhi(i,j,k)=minDis;
					else mPhi(i,j,k)=-minDis;
				}
			}//if
		}//FOR_EACH_CELL
	}
	if(unTagedNum==theDim[0]*theDim[1]*theDim[2]){
		std::cout<<"!"<<std::endl;
		getchar();
	}
}

void MACGrid::advectY(double dt) 
{       
	target.mY = mY; 
    // TODO: Calculate new temp and store in target.
	FOR_EACH_CELL
	{
		glm::dvec3 pos = getCenter(i,j,k);
		double Ynew;
		//core
		if(mPhi(i,j,k) < 0)
		{
			/*
			glm::dvec3 vf = getVelocityf(pos);
			glm::dvec3 oldPos = pos - vf *dt;
			//double Ymid = mY.interpolate(oldPos) ;
			//Ynew = -AdvectYk * dt + Ymid; 
			Ynew = mY.interpolate(oldPos) ;
			*/
			Ynew = 1.0;
		}
		else //make Y falling
		{
			glm::dvec3 vh = getVelocityh(pos);
			glm::dvec3 oldPos = pos - vh *dt;
			double Ymid = mY.interpolate(oldPos) ;
			Ynew = -AdvectYk * dt + Ymid; 
		}
		target.mY(i,j,k) = Ynew;
	}	
    mY = target.mY;
}

//fire simUflation

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) 
{
	double result = 0.0;
	FOR_EACH_CELL 
	{
		result = std::max(result, abs(vector(i,j,k)));
	}
	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);
		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}
}

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mDensityf(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0);

   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         glm::dvec3 pos = getCenter(i,j,k);
         glm::dvec3 vel = getVelocityf(pos);
         if (glm::length(vel) > 0.0001)
         {
           //vel.Normalize(); // PETER KUTZ.
           vel *= theCellSize/2.0;
           vel += pos;
		       glColor4f(1.0, 1.0, 0.0, 1.0);

           GLdouble doublePos[3];
           doublePos[0] = pos.x, doublePos[1] = pos.y, doublePos[2] = pos.z;
           glVertex3dv(doublePos);
		       
           GLdouble doubleVel[3];
           glColor4f(0.0, 1.0, 0.0, 1.0);
           doubleVel[0] = vel.x, doubleVel[1] = vel.y, doubleVel[2] = vel.z;
           glVertex3dv(doubleVel);
         }
      }
   glEnd();
}

glm::dvec4 MACGrid::getRenderColor(int i, int j, int k)
{
	double tvalue = mTemp(i, j, k); 
	if(tvalue < CoolT)
	{
		return glm::dvec4(0,0,0,0);
	}
	glm::dvec3 rgb = getRenderRGB(tvalue);
	double dvalue;
	if(mPhi(i,j,k) >= 0)
	{
		dvalue = mDensityh(i, j, k); 
	}
	else
	{
		dvalue = mDensityf(i, j, k); 
	}
    return glm::dvec4(rgb[0], rgb[1], rgb[2], dvalue);
}

glm::dvec3 MACGrid::getRenderRGB(double t)
{
	return spect.getRGB(t);	
}

glm::dvec4 MACGrid::getRenderColor(const glm::dvec3& pt)
{
	double tvalue = mTemp.interpolate(pt); 	
	tvalue = std::min(tvalue, Tmax);
	
	double phi = mPhi.interpolate(pt);
	double dvalue;
	dvalue = mDensityf.interpolate(pt);	
	glm::dvec3 rgb; 
	if(tvalue < Tignition)
	{
		rgb = getRenderRGB(Tignition) * tvalue / Tignition;
	}
	else
	{
		rgb = getRenderRGB(tvalue);
	}
	/*
	if(tvalue > Tmax)
	{
		printf("%lf %lf\n", phi, tvalue);
	}
	*/

	//glm::dvec3 test=getRenderRGB(40);
	//return glm::dvec4(test[0],test[1],test[2],1);

	//if(fabs(phi)<1) return glm::dvec4(1-fabs(phi),0,0,1);
	//else return glm::dvec4(0,0,0,1);

	if(phi < 0.0) 
	{
		dvalue /= fuelDensity;
		//printf("%lf %lf %lf\n", tvalue, phi, dvalue);
		return glm::dvec4(rgb[0], rgb[1], rgb[2], dvalue);
		//return glm::dvec4(0,1,0,1);
		
	}
	else
	{
		return glm::dvec4(rgb[0], rgb[1], rgb[2], dvalue);
	}
	/*
	if(phi > 0)
	{
		dvalue = mDensityf.interpolate(pt);
		return glm::dvec4(rgb[0], rgb[1], rgb[2], dvalue);
	}
    //return glm::dvec4(rgb[0], rgb[1], rgb[2], dvalue);
	*/
	
	/*
	double v=mPhi.interpolate(pt);
	if(v>1) 
	//if(v >= 0)
		return glm::dvec4(0,1,0,1);
	else if(v<-1) 
	//else if(v < 0)
			return glm::dvec4(1,1,1,1);
	else 
		return glm::dvec4(1-fabs(v),0,0,1);
	*/
}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            glm::dvec3 pos1 = glm::dvec3(i,j,k); 
            glm::dvec3 pos2 = glm::dvec3(i, j+stepsize, k); 

            glm::dvec4 color1 = getRenderColor(pos1);
            glm::dvec4 color2 = getRenderColor(pos2);

            glColor4dv(glm::value_ptr(color1));
            glVertex3dv(glm::value_ptr(pos1));

            glColor4dv(glm::value_ptr(color2));
            glVertex3dv(glm::value_ptr(pos2));
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            glm::dvec3 pos1 = glm::dvec3(i,j,k); 
            glm::dvec3 pos2 = glm::dvec3(i, j+stepsize, k); 

            glm::dvec4 color1 = getRenderColor(pos1);
            glm::dvec4 color2 = getRenderColor(pos2);

            glColor4dv(glm::value_ptr(color1));
            glVertex3dv(glm::value_ptr(pos1));

            glColor4dv(glm::value_ptr(color2));
            glVertex3dv(glm::value_ptr(pos2));
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            glm::dvec3 pos1 = glm::dvec3(i,j,k); 
            glm::dvec3 pos2 = glm::dvec3(i, j+stepsize, k); 

            glm::dvec4 color1 = getRenderColor(pos1);
            glm::dvec4 color2 = getRenderColor(pos2);

            glColor4dv(glm::value_ptr(color1));
            glVertex3dv(glm::value_ptr(pos1));

            glColor4dv(glm::value_ptr(color2));
            glVertex3dv(glm::value_ptr(pos2));
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            glm::dvec3 pos1 = glm::dvec3(i,j,k); 
            glm::dvec3 pos2 = glm::dvec3(i, j+stepsize, k); 

            glm::dvec4 color1 = getRenderColor(pos1);
            glm::dvec4 color2 = getRenderColor(pos2);

            glColor4dv(glm::value_ptr(color1));
            glVertex3dv(glm::value_ptr(pos1));

            glColor4dv(glm::value_ptr(color2));
            glVertex3dv(glm::value_ptr(pos2));
         } 
         glEnd();
         j+=stepsize;
      }
   }
}


void MACGrid::drawSmoke(const Camera& c)
{
   glm::dvec3 eyeDir = c.getBackward();
   double zresult = fabs(glm::dot(eyeDir, glm::dvec3(1,0,0)));
   double xresult = fabs(glm::dot(eyeDir, glm::dvec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, glm::dvec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = glm::length((cube.pos - c.getPosition()));
      cubes.insert(std::make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; ++i)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; ++i)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(glm::value_ptr(cube.color));
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(glm::value_ptr(cube.color));
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
