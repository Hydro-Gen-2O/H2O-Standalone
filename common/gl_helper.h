#ifndef GL_HELPER
	#define GL_HELPER

	#include "common_defs.h"

	#include <gl/glee.h>
	#include <gl/glext.h>	
	#include <gl/glut.h>
	#include <iostream>

	extern void checkOpenGL ();
	extern void drawText ( int x, int y, char* msg);
	extern void drawGrid ();

	extern int			tm_cnt;
	extern float		tm_fps;
	
	extern void checkFrameBuffers ();

	extern GLuint glSphere;
	extern float  glRadius;
	extern void setSphereRadius ( float f );
	extern void drawSphere ();

#endif