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

#ifndef DEF_GEOM
	#define DEF_GEOM

	#include <vector>

	#define	HEAP_MAX			2147483640	// largest heap size (range of hpos)

	#define	ELEM_MAX			2147483640	// largest number of elements in a buffer (range of hval)
	//#define ELEM_MAX			32768		// largest number of elements in a buffer (range of hval)

	#define BUF_UNDEF			255

	#define FPOS				2			// free position offsets
	typedef unsigned char		uchar;
	typedef unsigned short		ushort;
	typedef signed int			hpos;		// pointers into heap
	typedef signed int			hval;		// values in heap	
	typedef hval				href;		// values are typically references 
	struct hList {
		ushort		cnt;
		ushort		max;
		hpos		pos;
	};

	class GeomBuf {
	public:
		GeomBuf()	{ dtype = 0; num = 0; max = 0; stride = 0; data = 0x0; }		
		uchar		dtype;
		hval		num;
		hval		max;
		long		size;
		ushort		stride;		
		char*		data;
	};

	class GeomX {
	public:
		GeomX ();
		// Basic geometry setup	
		void FreeBuffers ();
		void ResetBuffer ( uchar b, int n );
		int AddBuffer ( uchar typ, ushort stride, int max );
		
		int NumElem ( uchar b )				{ if ( b==BUF_UNDEF) return 0; else return mBuf[b].num; }		
		char* GetElem ( uchar b, int n )	{ return mBuf[b].data + n*mBuf[b].stride; }
		char* RandomElem(uchar b, href& ndx);
		char* AddElem(uchar b, href& pos);
		
	protected:
		std::vector< GeomBuf > mBuf;
	};

#endif
