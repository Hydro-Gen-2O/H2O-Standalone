#include "common_defs.h"
#include "gl_helper.h"

#include <math.h>

// Shadow Light
float light_proj[16];
float light_x, light_y, light_z;
float light_tox, light_toy, light_toz;
float light_mfov;

// Fonts
void *font = GLUT_BITMAP_8_BY_13;
void *fonts[] = {GLUT_BITMAP_9_BY_15,
				 GLUT_BITMAP_TIMES_ROMAN_10,
				GLUT_BITMAP_TIMES_ROMAN_24};
// Timing
mint::Time	tm_last;
int			tm_cnt;
float		tm_fps;

GLuint glSphere = 65535;
float  glRadius = 0.0;

void setSphereRadius ( float r )
{
	if ( glRadius == r ) return;
	glRadius = r;

	// GL sphere
	if ( glSphere != 65535 ) glDeleteLists ( glSphere, 1 );
	glSphere = glGenLists ( 1 );
	float x, y, z, x1, y1, z1;	
	glNewList ( glSphere, GL_COMPILE );
		glBegin ( GL_TRIANGLE_STRIP );
		for ( float tilt=-90; tilt <= 90; tilt += 10.0) {
			for ( float ang=0; ang <= 360; ang += 30.0) {
				x = sin ( ang*DEGtoRAD) * cos ( tilt*DEGtoRAD );
				y = cos ( ang*DEGtoRAD) * cos ( tilt*DEGtoRAD );
				z = sin ( tilt*DEGtoRAD ) ;
				x1 = sin ( ang*DEGtoRAD) * cos ( (tilt+10.0)*DEGtoRAD ) ;
				y1 = cos ( ang*DEGtoRAD) * cos ( (tilt+10.0)*DEGtoRAD ) ;
				z1 = sin ( (tilt+10.0)*DEGtoRAD );
				glNormal3f ( x, y, z );		glVertex3f ( x*r, y*r, z*r );		
				glNormal3f ( x1, y1, z1 );	glVertex3f ( x1*r, y1*r, z1*r );
			}
		}
		glEnd ();
	glEndList ();
}

void drawSphere ()
{
	if ( glRadius == 0.0 ) setSphereRadius ( 1.0 );
	glCallList ( glSphere );
}

// Check if there have been any openGL problems
void checkOpenGL ()
{
	GLenum errCode = glGetError();
	if (errCode != GL_NO_ERROR) {
		const GLubyte* errString = gluErrorString(errCode);
		fprintf( stderr, "OpenGL error: %s\n", errString );
	}
}

void drawText ( int x, int y, char* msg)
{
  int len, i;
  glRasterPos2f(x, y);
  len = (int) strlen(msg);
  for (i = 0; i < len; i++) 
    glutBitmapCharacter(font, msg[i]);  
}

void drawGrid ()
{
	glColor3f ( 0.3, 0.3, 0.3 );
	glBegin ( GL_LINES );
	for (float x=-40; x<=40.0; x+=10.0 ) {
		glVertex3f ( x, -40.0, 0 );
		glVertex3f ( x,  40.0, 0 );
	}
	for (float y=-40; y<=40.0; y+=10.0 ) {
		glVertex3f ( -40.0, y, 0 );
		glVertex3f (  40.0, y, 0 );
	}
	glEnd ();
}

void measureFPS ()
{
	// Measure FPS
	mint::Time tm_elaps;	
	if ( ++tm_cnt > 5 ) {		
		tm_elaps.SetSystemTime ( ACC_NSEC );			// get current sytem time - accurate to 1 ns
		tm_elaps = tm_elaps - tm_last;					// get elapsed time from 5 frames ago
		tm_fps = 5.0 * 1000.0 / tm_elaps.GetMSec ();	// compute fps
		tm_cnt = 0;										// reset frame counter
		tm_last.SetSystemTime ( ACC_NSEC );
	}
}

void checkFrameBuffers ()
{                                                            
	GLenum status;                                             
	status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);  
	switch(status) {                                          
	case GL_FRAMEBUFFER_COMPLETE_EXT: printf ( "FBO complete\n" ); break;                                                
	case GL_FRAMEBUFFER_UNSUPPORTED_EXT: printf ( "FBO format unsupported\n"); break;                                                
	default:  printf ( "Unknown FBO error\n");
	}
}