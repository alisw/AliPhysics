// -*- C++ -*-
// 
// 1998/10/19
// ---------------------------------------------------------------------------
//
// AliGTRD1 Class
//
// This file is part of the ALICE Geometry Database .
//
// Author:  Joana E. Santo & company...
//
// ---------------------------------------------------------------------------

#include <TView.h>
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <iostream.h>
#include <TGLKernelABC.h>
#include "AliGTRD1.h"
#include <TROOT.h>

ClassImp(AliGTRD1)

//-------------------------------------------------------------------------

AliGTRD1::AliGTRD1(Text_t *name, Text_t *title, Float_t dx1, Float_t dx2, Float_t dy, Float_t dz) : AliGShape(name, title)
{
//*-*-*-*-*-*-*-*-*-*-*-*-*AliGTRD1 shape normal constructor*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                      =============================

    fDx2 = dx2;
    fDx1 = dx1;
    fDy  = dy;
    fDz  = dz;
}

//-------------------------------------------------------------------------   

AliGTRD1::AliGTRD1()
{
    /* Default Constructor */

    fDx2   = 0;        
    fName  = "";
    fTitle = "";
}

//-------------------------------------------------------------------------

AliGTRD1::AliGTRD1( AliGTRD1* trd1 )
{
    /* Copy Constructor */
    fColor     = trd1->fColor;
    fDx1       = trd1->fDx1;
    fDx2       = trd1->fDx2;
    fDy        = trd1->fDy; 
    fDz        = trd1->fDz;
    fName      = trd1->fName;
    fTitle     = trd1->fTitle;
}

//----------------------------------------------------------------------------

void AliGTRD1::SetPoints(Float_t *buff)
{
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*Create AliGTRD1 points*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                            ==================

    Float_t dx1, dx2, dy, dz;

    dx1 = fDx1;
    dx2 = fDx2;
    dy  = fDy;
    dz  = fDz;
    
    if (buff) {
        buff[ 0] = -dx1;  buff[ 1] = -dy;  buff[ 2] = -dz;
        buff[ 3] =  dx1;  buff[ 4] = -dy;  buff[ 5] = -dz;
        buff[ 6] =  dx1;  buff[ 7] =  dy;  buff[ 8] = -dz;
        buff[ 9] = -dx1;  buff[10] =  dy;  buff[11] = -dz;
        buff[12] = -dx2;  buff[13] = -dy;  buff[14] =  dz;
        buff[15] =  dx2;  buff[16] = -dy;  buff[17] =  dz;
        buff[18] =  dx2;  buff[19] =  dy;  buff[20] =  dz;
        buff[21] = -dx2;  buff[22] =  dy;  buff[23] =  dz;
    }
}

//-------------------------------------------------------------------------

void AliGTRD1::Paint(Option_t *option)
{
    SetLineColor( GetCol() );

    const Int_t numpoints = 8;

    //*-* Allocate memory for points *-*

    Float_t *points = new Float_t[3*numpoints];
    if (!points) return;

    SetPoints(points);

    if (gPad->GetView3D()) PaintGLPoints(points);

    //==  for (Int_t i = 0; i < numpoints; i++)
    //            gNode->Local2Master(&points[3*i],&points[3*i]);

    Int_t c = ((GetLineColor() % 8) - 1) * 4;     // Basic colors: 0, 1, ... 7
    if (c < 0) c = 0;

    //*-* Allocate memory for segments *-*

    X3DBuffer *buff = new X3DBuffer;
    if (buff) {
        buff->numPoints = 8;
        buff->numSegs   = 12;
        buff->numPolys  = 6;
    }

//*-* Allocate memory for points *-*

    buff->points = points;
    buff->segs = new Int_t[buff->numSegs*3];

    if (buff->segs) {
        buff->segs[ 0] = c;    buff->segs[ 1] = 0;    buff->segs[ 2] = 1;
        buff->segs[ 3] = c+1;  buff->segs[ 4] = 1;    buff->segs[ 5] = 2;
        buff->segs[ 6] = c+1;  buff->segs[ 7] = 2;    buff->segs[ 8] = 3;
        buff->segs[ 9] = c;    buff->segs[10] = 3;    buff->segs[11] = 0;
        buff->segs[12] = c+2;  buff->segs[13] = 4;    buff->segs[14] = 5;
        buff->segs[15] = c+2;  buff->segs[16] = 5;    buff->segs[17] = 6;
        buff->segs[18] = c+3;  buff->segs[19] = 6;    buff->segs[20] = 7;
        buff->segs[21] = c+3;  buff->segs[22] = 7;    buff->segs[23] = 4;
        buff->segs[24] = c;    buff->segs[25] = 0;    buff->segs[26] = 4;
        buff->segs[27] = c+2;  buff->segs[28] = 1;    buff->segs[29] = 5;
        buff->segs[30] = c+1;  buff->segs[31] = 2;    buff->segs[32] = 6;
        buff->segs[33] = c+3;  buff->segs[34] = 3;    buff->segs[35] = 7;
    }

//*-* Allocate memory for polygons *-*

    buff->polys = new Int_t[buff->numPolys*6];

    if (buff->polys) {
        buff->polys[ 0] = c;   buff->polys[ 1] = 4;  buff->polys[ 2] = 0;
        buff->polys[ 3] = 9;   buff->polys[ 4] = 4;  buff->polys[ 5] = 8;
        buff->polys[ 6] = c+1; buff->polys[ 7] = 4;  buff->polys[ 8] = 1;
        buff->polys[ 9] = 10;  buff->polys[10] = 5;  buff->polys[11] = 9;
        buff->polys[12] = c;   buff->polys[13] = 4;  buff->polys[14] = 2;
        buff->polys[15] = 11;  buff->polys[16] = 6;  buff->polys[17] = 10;
        buff->polys[18] = c+1; buff->polys[19] = 4;  buff->polys[20] = 3;
        buff->polys[21] = 8;   buff->polys[22] = 7;  buff->polys[23] = 11;
        buff->polys[24] = c+2; buff->polys[25] = 4;  buff->polys[26] = 0;
        buff->polys[27] = 3;   buff->polys[28] = 2;  buff->polys[29] = 1;
        buff->polys[30] = c+3; buff->polys[31] = 4;  buff->polys[32] = 4;
        buff->polys[33] = 5;   buff->polys[34] = 6;  buff->polys[35] = 7;
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

    delete [] points;
    if (buff->segs)     delete [] buff->segs;
    if (buff->polys)    delete [] buff->polys;
    if (buff)           delete    buff;
}

