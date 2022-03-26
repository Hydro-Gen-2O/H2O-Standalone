/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/
#include <iostream>
extern "C" { FILE __iob_func[3] = { *stdin,*stdout,*stderr }; }

#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "common_defs.h"

#include "fluid_system.h"
#include "gl_helper.h"

#include <gl/glut.h>
// Globals
FluidSystem			psys;

float window_width  = 1024;
float window_height = 768;

Vector3DF	cam_from, cam_angs, cam_to;			// Camera stuff
float		cam_fov;	

int		psys_demo = -1; // CHANGE TO CHANGE  SIM TYPE
int		psys_nmax = 4096;

bool	bHelp = false;						// Toggles
int		iShade = 1;			
int		iClrMode = 0;
bool	bPntDraw = false;
bool    bPause = false;

// View matricies
float view_matrix[16];					// View matrix (V)
float model_matrix[16];					// Model matrix (M)
float proj_matrix[16];					// Projective matrix

// Different things we can move around

// Mouse control
#define DRAG_OFF		0				// mouse states
#define DRAG_LEFT		1
#define DRAG_RIGHT		2
int		last_x = -1, last_y = -1;		// mouse vars
int		dragging = 0;

void drawScene ( float* viewmat, bool bShade )
{
	if ( iShade <= 1 && bShade ) {		
		glEnable ( GL_LIGHT0 );
		GLfloat diff[4];
		GLfloat spec[4];
		GLfloat shininess = 60.0;
		
		diff[0] = 0.8f; diff[1] = 0.8f; diff[2] = 0.8f; diff[3] = 1.0f;
		spec[0] = 1.0f; spec[1] = 1.0f; spec[2] = 1.0f; spec[3] = 1.0f;
		glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, &diff[0]);
		glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, &spec[0]);
		glMaterialfv (GL_FRONT_AND_BACK, GL_SHININESS, &shininess);
		glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );

		glColor3f ( 1, 1, 1 );
		glLoadMatrixf ( viewmat );
		glBegin ( GL_LINES );
		for (float n=-100; n <= 100; n += 20.0 ) {
			glVertex3f ( -100, n, 0.1 );
			glVertex3f ( 100, n, 0.1 );
			glVertex3f ( n, -100, 0.1 );
			glVertex3f ( n,  100, 0.1 );
		}
		glEnd ();

		psys.Draw ( &viewmat[0], 0.8 );				// Draw particles		
		psys.SPH_DrawDomain();
	} else {
		glDisable ( GL_LIGHTING );
		psys.Draw ( &viewmat[0], 0.55 );			// Draw particles
	}
}

void computeFromPositions ()
{
	cam_from.x = cam_to.x + sin( cam_angs.x * DEGtoRAD) * sin( cam_angs.y * DEGtoRAD) * cam_angs.z;
	cam_from.y = cam_to.y + -cos( cam_angs.x * DEGtoRAD) * sin( cam_angs.y * DEGtoRAD) * cam_angs.z;
	cam_from.z = cam_to.z + cos( cam_angs.y * DEGtoRAD) * cam_angs.z;	
}

void computeProjection ()
{
	// ---- Create projection matrix for eye-coordinate transformations 
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective ( cam_fov, window_width / ( float ) window_height, 10.0, 800.0 );
	glPushMatrix (); 
	glGetFloatv ( GL_MODELVIEW_MATRIX, proj_matrix ); 
	glPopMatrix ();
}

void computeView ()
{
	glMatrixMode ( GL_MODELVIEW );
	glLoadIdentity ();
	gluLookAt ( cam_from.x, cam_from.y, cam_from.z, cam_to.x, cam_to.y, cam_to.z, 0, 0, 1 );
	glPushMatrix (); 
	glGetFloatv ( GL_MODELVIEW_MATRIX, view_matrix ); 
	glPopMatrix ();
}	

void display () 
{
	// Do simulation!
	if ( !bPause ) psys.Run ();

	glEnable ( GL_DEPTH_TEST );

	// Clear frame buffer
	glClearColor ( 0, 0, 0, 0 );
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable ( GL_CULL_FACE );
	glShadeModel ( GL_SMOOTH );

	// Compute camera view
	computeFromPositions ();
	computeProjection ();
	computeView ();

	// Draw 3D
	glEnable ( GL_LIGHTING );  
	glLoadMatrixf ( view_matrix );	
	drawScene ( view_matrix, true );
 
	// Swap buffers
	glutSwapBuffers();  
	glutPostRedisplay();
}

void reshape ( int width, int height ) 
{
  window_width  = (float) width;
  window_height = (float) height;
  glViewport( 0, 0, width, height );  
}

void keyboard_func ( unsigned char key, int x, int y )
{
	// m: run sim with double particles; n run sim with half
	switch( key ) {
	case 'x': case 'X':
		if ( ++iClrMode > 2) iClrMode = 0;
		psys.SetParam ( CLR_MODE, iClrMode );
		break;
	case 'd': case 'D': {
		int d = psys.GetParam ( PNT_DRAWMODE ) + 1;
		if ( d > 2 ) d = 0;
		psys.SetParam ( PNT_DRAWMODE, d );
		} break;
	case 's': case 'S':	if ( ++iShade > 2 ) iShade = 0;		break;
	case ' ':
		bPause = !bPause;	break;
	default:
	break;
  }
}


void mouse_click_func ( int button, int state, int x, int y )
{
  if( state == GLUT_DOWN ) {
    if ( button == GLUT_LEFT_BUTTON )		dragging = DRAG_LEFT;
    else if ( button == GLUT_RIGHT_BUTTON ) dragging = DRAG_RIGHT;	
    last_x = x;
    last_y = y;	
  } else {
    dragging = DRAG_OFF;
  }
}

void mouse_move_func ( int x, int y )
{
	int dx = x - last_x;
	int dy = y - last_y;

	if ( dragging == DRAG_LEFT ) {
		cam_angs.x += dx;
		cam_angs.y += dy;
		if ( cam_angs.x >= 360.0 )	cam_angs.x -= 360.0;
		if ( cam_angs.x < 0 )		cam_angs.x += 360.0;
		if ( cam_angs.y >= 180.0 )	cam_angs.y = 180.0;
		if ( cam_angs.y <= -180.0 )	cam_angs.y = -180.0;
		//printf ( "Cam Ang: %f %f %f\n", cam_angs.x, cam_angs.y, cam_angs.z );
		//printf ( "Cam To:  %f %f %f\n", cam_to.x, cam_to.y, cam_to.z );
		//printf ( "Cam FOV: %f\n", cam_fov);
	} else if ( dragging == DRAG_RIGHT ) {
		cam_angs.z += dy*.15;
		if ( cam_angs.z < 0)		cam_angs.z = 0;
		//printf ( "Cam Ang: %f %f %f\n", cam_angs.x, cam_angs.y, cam_angs.z );
		//printf ( "Cam To:  %f %f %f\n", cam_to.x, cam_to.y, cam_to.z );
		//printf ( "Cam FOV: %f\n", cam_fov );
	}

	if ( x < 10 || y < 10 || x > 1000 || y > 700 ) {
		glutWarpPointer ( 1024/2, 768/2 );
		last_x = 1024/2;
		last_y = 768/2;
	} else {
		last_x = x;
		last_y = y;
	}
}


void idle_func ()
{
}

void init ()
{
	glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);	
	glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

	srand ( time ( 0x0 ) );

	glClearColor( 0.49, 0.49, 0.49, 1.0 );
	glShadeModel( GL_SMOOTH );

	glEnable ( GL_COLOR_MATERIAL );
	glEnable (GL_DEPTH_TEST);  
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
	glDepthMask ( 1 );
	glEnable ( GL_TEXTURE_2D );

	// callbacks
	glutDisplayFunc( display );
	glutReshapeFunc( reshape );
	glutKeyboardFunc( keyboard_func );
	glutMouseFunc( mouse_click_func );  
	glutMotionFunc( mouse_move_func );
	glutIdleFunc( idle_func );
	glutSetCursor ( GLUT_CURSOR_NONE );

	cam_angs.x = 29;		cam_angs.y = 75;		cam_angs.z = 80.0;
	cam_to.x = 0;		cam_to.y = 0;		cam_to.z = 5;
	cam_fov = 35.0;

	psys.SPH_CreateExample ( psys_demo, psys_nmax );

	psys.SetParam ( PNT_DRAWMODE, int(bPntDraw ? 1:0) );
	psys.SetParam ( CLR_MODE, iClrMode );	
}


int main ( int argc, char **argv )
{
	glutInit( &argc, &argv[0] ); 
	glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowPosition( 100, 100 );
	glutInitWindowSize( (int) window_width, (int) window_height );
	glutCreateWindow ( "Fluids v.1 (c) 2008, R. Hoetzlein (ZLib)" );

//	glutFullScreen ();
 
	// initialize parameters
	init();

	// wait for something to happen
	glutMainLoop();

  return 0;
}
