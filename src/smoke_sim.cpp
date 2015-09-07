#include "smoke_sim.h"
#include "constants.h"
#include "open_gl_headers.h"
#include "stb_image_write.h"
#include "custom_output.h"
#include "basic_math.h"

SmokeSim::SmokeSim() : mFrameNum(0), mTotalFrameNum(0), mRecordEnabled(false)
{
   reset();
}

SmokeSim::~SmokeSim()
{
}

void SmokeSim::reset()
{
   mGrid.reset();
	mTotalFrameNum = 0;
}

void SmokeSim::step()
{
	double dt = 0.03;
	// TODO : Write the simulation step
	mGrid.updateSources();
	mGrid.findBoundaryGradient();
	mGrid.advectBoundary(dt,mTotalFrameNum==1);
	mGrid.findBoundary(mTotalFrameNum==-1);
	mGrid.insideBoundary();
	mGrid.reconditioning();
	std::cout<<"-----"<<mTotalFrameNum<<"-----"<<std::endl;

	mGrid.advectVelocity(dt,false);
	mGrid.addExternalForces(dt,false);
	mGrid.project(dt);
	mGrid.advectY(dt);
	mGrid.advectTemperature(dt,false);

	mGrid.advectDensity(dt);
	mTotalFrameNum++;
}

void SmokeSim::setRecording(bool on, int width, int height)
{
   if (on && ! mRecordEnabled)  // reset counter
   {
      mFrameNum = 0;
   }
   mRecordEnabled = on;
	
	recordWidth = width;
	recordHeight = height;
}

bool SmokeSim::isRecording()
{
   return mRecordEnabled;
}

void SmokeSim::draw(const Camera& c)
{
   drawAxes(); 
   mGrid.draw(c);
   if (mRecordEnabled) grabScreen();
}

void SmokeSim::drawAxes()
{
	glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
		glDisable(GL_LIGHTING);

		glLineWidth(2.0); 
		glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(1.0, 0.0, 0.0);

			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 1.0, 0.0);

			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 1.0);
		glEnd();
	glPopAttrib();
}

void SmokeSim::grabScreen()  // Code adapted from asst#1 // OVERHAULED BY PETER KUTZ. USING STB_IMAGE_WRITE INSTEAD OF DEVIL.
{	
	if (mFrameNum > 9999) exit(0);
	

	// TODO: Un-comment this to save each frame of smoke in CIS 460 volumetric text file format.
	// Save the smoke densities:
	char smoke_filename[2048];
	sprintf_s(smoke_filename, 2048, "%d.txt", mFrameNum); // Use snprintf if your compiler supports C99.
	mGrid.saveSmoke(smoke_filename);
	
	// Save an image:

	unsigned char* bitmapData = new unsigned char[3 * recordWidth * recordHeight];

	for (int i=0; i<recordHeight; i++) 
	{
		glReadPixels(0,i,recordWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
			bitmapData + (recordWidth * 3 * ((recordHeight-1)-i)));
	}

	char anim_filename[2048];
	sprintf_s(anim_filename, 2048, "smoke_%04d.png", mFrameNum); 
	
	stbi_write_png(anim_filename, recordWidth, recordHeight, 3, bitmapData, recordWidth * 3);
	
	delete [] bitmapData;
	
	
	mFrameNum++;
}

int SmokeSim::getTotalFrames() {
	return mTotalFrameNum;
}
