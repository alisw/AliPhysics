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

// -*- C++ -*-
// 
// 1998/10/19
// ---------------------------------------------------------------------------
//
// AliTube Class
//
// This file is part of the ALICE Geometry Database .
//
// By:  Joana E. Santo & David Collados
//
// ---------------------------------------------------------------------------

#include <TVirtualPad.h>
#include <TView.h>
#include <TCanvas.h>
#include <iostream.h>
#include <TGLKernelABC.h>
#include "TROOT.h"
#include "AliGTube.h"
#include "AliGShape.h"

ClassImp(AliGTube)

AliGTube::AliGTube()
{
    /* Default Constructor */
    fCoTab = NULL;    // Table of cos(fPhi1) .... cos(fPhil+fDphi1)
    fAspectRatio = 1; // defines  (the ellipse semi-axis in Y)/(the ellipse semi-axis in X)
    fDz    = 0.;      // half length in z
    fName  = "";
    fNdiv  = 0;       // number of segments (precision)
    fRmax  = 0.;      // ellipse  semi-axis   in  X outside
    fRmin  = 0.;      // ellipse  semi-axis   in  X inside
    fSiTab = NULL;    // Table of sin(fPhi1) .... sin(fPhil+fDphi1)
    fTitle = "";
}

//-------------------------------------------------------------------------

AliGTube::AliGTube( Text_t *name, Text_t *title, Float_t rmin, Float_t rmax, Float_t dz, Float_t aspect ) : AliGShape(name, title)
{
    /* Constructor */
    fCoTab       = NULL; // Table of cos(fPhi1) .... cos(fPhil+fDphi1)
    fAspectRatio = aspect; // defines  (the ellipse semi-axis in Y)/(the ellipse semi-axis in X)
    fDz          = dz;   // half length in z
    fNdiv        = 0;    // number of segments (precision)
    fRmax        = rmax; // ellipse  semi-axis   in  X outside
    fRmin        = rmin; // ellipse  semi-axis   in  X inside
    fSiTab       = NULL; // Table of sin(fPhi1) .... sin(fPhil+fDphi1)

    MakeTableOfCoSin();
}


//-------------------------------------------------------------------------

AliGTube::AliGTube( AliGTube *tube ) 
{
    /* Copy Constructor */
    fCoTab       = NULL; // Table of cos(fPhi1) .... cos(fPhil+fDphi1)
    fAspectRatio = tube->fAspectRatio; // defines  (the ellipse semi-axis in Y)/(the ellipse semi-axis in X)
    fDz          = tube->fDz;   // half length in z
    fColor       = tube->fColor;
    fName        = tube->fName;
    fNdiv        = 0;    // number of segments (precision)
    fRmax        = tube->fRmax; // ellipse  semi-axis   in  X outside
    fRmin        = tube->fRmin; // ellipse  semi-axis   in  X inside
    fSiTab       = NULL; // Table of sin(fPhi1) .... sin(fPhil+fDphi1)
    fTitle       = tube->fTitle;

    MakeTableOfCoSin();
}

//-------------------------------------------------------------------------

AliGTube::AliGTube(Text_t *name, Text_t *title, Float_t rmax, Float_t dz) : AliGShape(name, title)
{
    /* Tube simplified constructor */
    fCoTab       = NULL; // Table of cos(fPhi1) .... cos(fPhil+fDphi1)
    fAspectRatio = 1;    // defines  (the ellipse semi-axis in Y)/(the ellipse semi-axis in X)
    fDz          = dz;   // half length in z
    fNdiv        = 0;    // number of segments (precision)
    fRmax        = rmax; // ellipse  semi-axis   in  X outside
    fRmin        = 0.;   // ellipse  semi-axis   in  X inside
    fSiTab       = NULL; // Table of sin(fPhi1) .... sin(fPhil+fDphi1)

    MakeTableOfCoSin();
}

//-------------------------------------------------------------------------

AliGTube::~AliGTube()
{
    /* Destructor */
    delete [] fCoTab;
    delete [] fSiTab;
}

//-------------------------------------------------------------------------

void AliGTube::MakeTableOfCoSin()
{
    const Double_t PI    = TMath::ATan(1) * 4.0;
    const Double_t TWOPI = 2*PI;

    Int_t j;
    Int_t n = GetNumberOfDivisions ();

    if (fCoTab)
        delete [] fCoTab; // Delete the old tab if any
    fCoTab = new Double_t [n];

    if (!fCoTab ) {
        Error("MakeTableOfCoSin()","No cos table done");
        return;
    }

    if (fSiTab)
        delete [] fSiTab; // Delete the old tab if any
    fSiTab = new Double_t [n];

    if (!fSiTab )
    {
        Error("MakeTableOfCoSin()","No sin table done");
        return;
    }

    Double_t range   = TWOPI;
    Double_t angstep = range/n;

    Double_t ph = 0;
    for (j = 0; j < n; j++)
    {
        ph = j*angstep;
        fCoTab[j] = TMath::Cos(ph);
        fSiTab[j] = TMath::Sin(ph);
    }

}

//-------------------------------------------------------------------------

void AliGTube::DrawShape(Option_t *option)
{
    Draw(option);
    gPad->Update();
}

//-------------------------------------------------------------------------

void AliGTube::Draw(Option_t *option)
{
    //cout << " Entra en " << this->GetName() << "::Draw " << endl;
    TString opt = option;
    opt.ToLower();

    if( !gPad ) {
        gPad = new TCanvas("AliGTube","AliGTube",0,0,400,300);
        gPad->Range(0,0,1,1);
        gPad->SetFillColor(32); // Light Green
        gPad->SetBorderSize(3);
        gPad->SetBorderMode(0); // -1 (down) 0 (no) 1 (up)
    }
    else {
        if( !opt.Contains("same") ) {
            gPad->Clear();
            gPad->SetName("AliGTube");
            gPad->SetTitle("AliGTube");
        }
        else {
            gPad->SetName("AliShapes");
            gPad->SetTitle("AliShapes"); 
        }
    }

    AppendPad(option);
    TView *view = gPad->GetView();

    if (!view)
        view = new TView(1);

    view->SetAutoRange(kTRUE);
    Paint(option);
    view->SetAutoRange(kFALSE);

    cout << " Sale de " << this->GetName() << "::Draw " << endl;
}

//-------------------------------------------------------------------------

void AliGTube::Paint(Option_t *option)
{
//*-*-*-*-*-*-*-*Paint this 3-D shape with its current attributes*-*-*-*-*-*-*-*
//*-*            ================================================

    SetLineColor( GetCol() );
        
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
        buff->numSegs   = n*8;
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
        for (i = 6; i < 8; i++) {
            for (j = 0; j < n; j++) {
                buff->segs[(i*n+j)*3  ] = c;
                buff->segs[(i*n+j)*3+1] = 2*(i-6)*n+j;
                buff->segs[(i*n+j)*3+2] = (2*(i-6)+1)*n+j;
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
    //*-* Paint in the pad
    Bool_t rangeView = strcmp(option,"range")==0 ? kTRUE : kFALSE;
   
    PaintShape(buff,rangeView);
    //PaintShape(buff);

    if (strstr(option, "x3d")) {
        if(buff && buff->points && buff->segs)
            FillX3DBuffer(buff);
        else {
            gSize3D.numPoints -= buff->numPoints;
            gSize3D.numSegs   -= buff->numSegs;
            gSize3D.numPolys  -= buff->numPolys;
        }
    }
    
   if( points ) delete [] points;
  
    //if (buff->points)   delete [] buff->points;
    
    if (buff->segs)     delete [] buff->segs;
    
    if (buff->polys)    delete [] buff->polys;
    
    if (buff)           delete    buff;
      
}

//-------------------------------------------------------------------------

void AliGTube::PaintGLPoints(Float_t *vertex)
{
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*Paint BRIK via OpenGL *-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                            =====================
    gGLKernel->PaintCone(vertex,GetNumberOfDivisions(),2);
}

//-------------------------------------------------------------------------

void AliGTube::SetPoints(Float_t *buff)
{
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*Create TUBE points*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                          ==================

    Float_t dz;
    Int_t j, n;

    n = GetNumberOfDivisions();

    dz = fDz;

    Int_t indx = 0;

    if (buff) {
//*-* We've to checxk whether the table does exist and create it
//*-* since fCoTab/fSiTab are not saved with any TShape::Streamer function
        if (!fCoTab)   MakeTableOfCoSin();

        for (j = 0; j < n; j++) {
            buff[indx+6*n] = buff[indx] = fRmin * fCoTab[j];
            indx++;
            buff[indx+6*n] = buff[indx] = fAspectRatio*fRmin * fSiTab[j];
            indx++;
            buff[indx+6*n] = dz;
            buff[indx]     =-dz;
            indx++;
        }
        for (j = 0; j < n; j++) {
            buff[indx+6*n] = buff[indx] = fRmax * fCoTab[j];
            indx++;
            buff[indx+6*n] = buff[indx] = fAspectRatio*fRmax * fSiTab[j];
            indx++;
            buff[indx+6*n]= dz;
            buff[indx]    =-dz;
            indx++;
        }
    }
}

//-------------------------------------------------------------------------

void AliGTube::Sizeof3D() const
{
//*-*-*-*-*-*Return total X3D size of this shape with its attributes*-*-*-*-*-*-*
//*-*        =======================================================

  cout << " Entra en AliGTube::Sizeof3D() " << endl;

    Int_t n = GetNumberOfDivisions();

    gSize3D.numPoints += n*4;
    gSize3D.numSegs   += n*8;
    gSize3D.numPolys  += n*4;
}

//-------------------------------------------------------------------------

void AliGTube::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliGTube.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliGShape::Streamer(R__b);
      R__b.ReadArray(fCoTab); //
      R__b >> fAspectRatio;
      R__b >> fDz;
      R__b >> fNdiv;
      R__b >> fRmax;
      R__b >> fRmin;
      R__b.ReadArray(fSiTab); //
   } else {
      R__b.WriteVersion(AliGTube::IsA());
      AliGShape::Streamer(R__b);
      R__b.WriteArray(fCoTab, GetNumberOfDivisions()); //
      R__b << fAspectRatio;
      R__b << fDz;
      R__b << fNdiv;
      R__b << fRmax;
      R__b << fRmin;
      R__b.WriteArray(fSiTab, GetNumberOfDivisions()); //
   }
}

//-------------------------------------------------------------------------

