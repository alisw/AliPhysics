// -*- C++ -*-
// 
// 1998/10/19
// ---------------------------------------------------------------------------
//
// AliGBox Class
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
#include "AliGBox.h"
#include <TROOT.h>

ClassImp(AliGBox)


//-------------------------------------------------------------------------

AliGBox::AliGBox(Text_t* name,Text_t* title, Float_t Dx, Float_t Dy, Float_t Dz) : AliGShape(name, title)
{
    /* Constructor */
    fDx = Dx;
    fDy = Dy;
    fDz = Dz;
}

//-------------------------------------------------------------------------   

AliGBox::AliGBox( AliGBox* box )
{
    if( box ) {
        /* Copy Constructor */
        fDx    = box->fDx;
        fDy    = box->fDy; 
        fDz    = box->fDz;
        fName  = box->fName;
        fTitle = box->fTitle;
        fColor = box->fColor;
    }
    else{
        /* Default Constructor */
        fDx    = 0;
        fDy    = 0;
        fDz    = 0;
        fName  = "";
        fTitle = "";
    }
}

//-------------------------------------------------------------------------

AliGBox* AliGBox::operator=(const AliGBox *box)
{
    /* Operator = */
    if( this == box) return this; // special case.

    fDx    = box->fDx;
    fDy    = box->fDy; 
    fDz    = box->fDz; 
    fName  = box->fName;
    fTitle = box->fTitle;
    fColor = box->fColor;

    return this;
}

//-------------------------------------------------------------------------

void AliGBox::DrawShape(Option_t *option)
{
    Draw(option);
    gPad->Update();
}

//-------------------------------------------------------------------------

void AliGBox::Paint(Option_t *option)
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




// ---------------------------------------------------------------------------

void AliGBox::PaintGLPoints(Float_t *vertex)
{
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*Paint BRIK via OpenGL *-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                            =====================
    gGLKernel->PaintBrik(vertex);
}


// ---------------------------------------------------------------------------

void AliGBox::SetPoints(Float_t *buff)
{
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*Create AliGBox points*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                            =====================
 
    if (buff) {
        buff[ 0] = -fDx;    buff[ 1] = -fDy;    buff[ 2] = -fDz;
        buff[ 3] =  fDx;    buff[ 4] = -fDy;    buff[ 5] = -fDz;
        buff[ 6] =  fDx;    buff[ 7] =  fDy;    buff[ 8] = -fDz;
        buff[ 9] = -fDx;    buff[10] =  fDy;    buff[11] = -fDz;
        buff[12] = -fDx;    buff[13] = -fDy;    buff[14] =  fDz;
        buff[15] =  fDx;    buff[16] = -fDy;    buff[17] =  fDz;
        buff[18] =  fDx;    buff[19] =  fDy;    buff[20] =  fDz;
        buff[21] = -fDx;    buff[22] =  fDy;    buff[23] =  fDz;
    }
}

//-------------------------------------------------------------------------

void AliGBox::Sizeof3D() const
{
//*-*-*-*-*-*-*Return total X3D size of this shape with its attributes*-*-*-*-*-*
//*-*          =======================================================

  //cout << " Entra en AliGBox::Sizeof3D() " << endl;

    gSize3D.numPoints += 8;
    gSize3D.numSegs   += 12;
    gSize3D.numPolys  += 6;
}








