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

#include <vector.h>
#include "geomx.h"
//#include "mdebug.h"

#include <assert.h>

GeomX::GeomX ()
{
}

void GeomX::FreeBuffers ()
{
	for (int n=0; n < (int) mBuf.size(); n++) {
		free ( mBuf[n].data );
		mBuf[n].data = 0x0;
	}
	mBuf.clear();
}

int GeomX::AddBuffer ( uchar typ, ushort stride, int max )
{
	GeomBuf buf;
	buf.dtype = typ;
	buf.stride = stride;
	buf.max = max;
	buf.num = 0;
	buf.size = 0;
	buf.data = (char*) malloc ( buf.max * buf.stride );
	mBuf.push_back ( buf );
	return (int) mBuf.size()-1;
}
void GeomX::ResetBuffer ( uchar b, int n )
{
	mBuf[b].max = n;		

	if ( mBuf[b].data != 0x0 ) free ( mBuf[b].data );

	char* new_data = (char*) malloc ( mBuf[b].max * mBuf[b].stride );
	
	mBuf[b].data = new_data;
	
	mBuf[b].num = 0;
	mBuf[b].size = mBuf[b].num*mBuf[b].stride;
}

char* GeomX::AddElem ( uchar b, href& ndx )
{
	if ( mBuf[b].num >= mBuf[b].max ) {
		if ( long(mBuf[b].max) * 2 > ELEM_MAX ) {
			std::cout << " Max elemnts " << std::endl;
			exit(-1);
		}
		mBuf[b].max *= 2;		
		char* new_data = (char*) malloc ( mBuf[b].max * mBuf[b].stride );
		memcpy ( new_data, mBuf[b].data, mBuf[b].num*mBuf[b].stride );
		free ( mBuf[b].data );
		mBuf[b].data = new_data;
	}
	mBuf[b].num++;
	mBuf[b].size += mBuf[b].stride; 
	ndx = mBuf[b].num-1;
	return mBuf[b].data + ndx*mBuf[b].stride;
}

char* GeomX::RandomElem ( uchar b, href& ndx )
{
	ndx = mBuf[b].num * rand() / RAND_MAX;
	return mBuf[b].data + ndx*mBuf[b].stride;
}