/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
*/

#include <AliMUONTUBE.h>
#include <TTUBE.h>
#include "TVirtualPad.h"
ClassImp(AliMUONTUBE)

AliMUONTUBE::AliMUONTUBE():TTUBE()
{;}

//_____________________________________________________________________________

AliMUONTUBE::AliMUONTUBE(Text_t *name, Text_t *title, Text_t *material, Float_t rmin, 
			 Float_t rmax, Float_t dz,Float_t aspect):
    TTUBE(name, title, material, rmin, rmax, dz, aspect)
{}

//______________________________________________________________________________
 AliMUONTUBE::AliMUONTUBE(Text_t *name, Text_t *title, Text_t *material, Float_t rmax, 
			  Float_t dz)
      : TTUBE(name, title, material, rmax, dz)
{
}

//______________________________________________________________________________
 void AliMUONTUBE::Paint(Option_t *option)
{
//*-*-*-*-*-*-*-*Paint this 3-D shape with its current attributes*-*-*-*-*-*-*-*
//*-*            ================================================

   Int_t i, j;
   Int_t n = GetNumberOfDivisions();
   const Int_t numpoints = 4*n;

//*-* Allocate memory for points *-*

   Float_t *points = new Float_t[3*numpoints];
   if (!points) return;

   SetPoints(points);

   if (gPad->GetView3D()) PaintGLPoints(points);

//==   for (i = 0; i < numpoints; i++)
//==            gNode->Local2Master(&points[3*i],&points[3*i]);

    X3DBuffer *buff = new X3DBuffer;
    if (buff) {
        buff->numPoints = numpoints;
        buff->numSegs   = n*6;
	if (strstr(option, "x3d"))  buff->numSegs   = n*8;
        buff->numPolys  = n*4;
    }

//*-* Allocate memory for points *-*

    buff->points = points;

    Int_t c = ((GetLineColor() % 8) - 1) * 4;     // Basic colors: 0, 1, ... 7
    if (c < 0) c = 0;

//*-* Allocate memory for segments *-*

    buff->segs = new Int_t[buff->numSegs*3];
    if (buff->segs) {
        for (i = 0; i < 4; i++) {
            for (j = 0; j < n; j++) {
                buff->segs[(i*n+j)*3  ] = c;
                buff->segs[(i*n+j)*3+1] = i*n+j;
                buff->segs[(i*n+j)*3+2] = i*n+j+1;
            }
            buff->segs[(i*n+j-1)*3+2] = i*n;
        }
	for (i = 4; i < 6; i++) {
	    for (j = 0; j < n; j++) {
		buff->segs[(i*n+j)*3  ] = c+1;
		buff->segs[(i*n+j)*3+1] = (i-4)*n+j;
		buff->segs[(i*n+j)*3+2] = (i-2)*n+j;
	    }
	}
	if (strstr(option, "x3d")) 
	{
	    for (i = 6; i < 8; i++) {
		for (j = 0; j < n; j++) {
		    buff->segs[(i*n+j)*3  ] = c;
		    buff->segs[(i*n+j)*3+1] = 2*(i-6)*n+j;
		    buff->segs[(i*n+j)*3+2] = (2*(i-6)+1)*n+j;
		}
	    }
	}
    }
//*-* Allocate memory for polygons *-*

    Int_t indx = 0;

    buff->polys = new Int_t[buff->numPolys*6];
    if (buff->polys) {
        for (i = 0; i < 2; i++) {
            for (j = 0; j < n; j++) {
                indx = 6*(i*n+j);
                buff->polys[indx  ] = c;
                buff->polys[indx+1] = 4;
                buff->polys[indx+2] = i*n+j;
                buff->polys[indx+3] = (4+i)*n+j;
                buff->polys[indx+4] = (2+i)*n+j;
                buff->polys[indx+5] = (4+i)*n+j+1;
            }
            buff->polys[indx+5] = (4+i)*n;
        }
        for (i = 2; i < 4; i++) {
            for (j = 0; j < n; j++) {
                indx = 6*(i*n+j);
                buff->polys[indx  ] = c+(i-2)*2+1;
                buff->polys[indx+1] = 4;
                buff->polys[indx+2] = (i-2)*2*n+j;
                buff->polys[indx+3] = (4+i)*n+j;
                buff->polys[indx+4] = ((i-2)*2+1)*n+j;
                buff->polys[indx+5] = (4+i)*n+j+1;
            }
            buff->polys[indx+5] = (4+i)*n;
        }
    }

    //*-* Paint in the pad
    Bool_t rangeView = strcmp(option,"range")==0 ? kTRUE : kFALSE;
    PaintShape(buff,rangeView);

    if (strstr(option, "x3d")) {
        if(buff && buff->points && buff->segs)
            FillX3DBuffer(buff);
        else {
            gSize3D.numPoints -= buff->numPoints;
            gSize3D.numSegs   -= buff->numSegs;
            gSize3D.numPolys  -= buff->numPolys;
        }
    }

    if (buff->points)   delete [] buff->points;
    if (buff->segs)     delete [] buff->segs;
    if (buff->polys)    delete [] buff->polys;
    if (buff)           delete    buff;
}

