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
// AliGPCone Class
//
// This file is part of the ALICE Geometry Database .
//
// By:  Joana E. Santo & David Collados
//
// ---------------------------------------------------------------------------

#include <TView.h>
#include <TVirtualPad.h>
#include <iostream.h>
#include <TGLKernelABC.h>
#include <TCanvas.h>
#include "AliGPCone.h"
#include "TROOT.h"

ClassImp(AliGPCone)

AliGPCone::AliGPCone() : AliGShape()
{
    /* Default Constructor */
    
    fPhi1  = 0.;
    fDphi1 = 0.;
    fNz    = 0;
    fNdiv  = 0;
    fRmin  = NULL;
    fRmax  = NULL;
    fDz    = NULL;
    fCoTab = NULL;
    fSiTab = NULL;
}

//-------------------------------------------------------------------------

AliGPCone::AliGPCone( AliGPCone* pcone )
{
    /* Copy Constructor */
    if( pcone->fNz < 2 ) {
        Error( pcone->GetName(), "number of z planes for %s must be at least two !", pcone->GetName() );
        return;
    }

    fName  = pcone->fName;
    fTitle = pcone->fTitle;

    fColor = pcone->fColor;
    fPhi1  = pcone->fPhi1;
    fDphi1 = pcone->fDphi1;
    fNz    = pcone->fNz;
    fNdiv  = pcone->fNdiv;

    fRmin  = new Float_t [fNz+1];
    fRmax  = new Float_t [fNz+1];
    fDz    = new Float_t [fNz+1];

    fCoTab = NULL;
    fSiTab = NULL;

    while (fDphi1 > 360) fDphi1 -= 360;

    MakeTableOfCoSin();

    for( int k=0; k<fNz; k++ ) {
        fDz[k]   = pcone->fDz[k];
        fRmin[k] = pcone->fRmin[k];
        fRmax[k] = pcone->fRmax[k];
    }

}

//-------------------------------------------------------------------------

AliGPCone::AliGPCone( Text_t *name, Text_t *title,  Float_t phi1, Float_t dphi1, Int_t nz ) : AliGShape(name, title)
{
    /* Constructor */
    
    if (nz < 2 ) {
        Error(name, "number of z planes for %s must be at least two !", name);
        return;
    }

    fPhi1  = phi1;
    fDphi1 = dphi1;
    fNz    = nz;
    fNdiv  = 0;
    fRmin  = new Float_t [fNz+1];
    fRmax  = new Float_t [fNz+1];
    fDz    = new Float_t [fNz+1];

    fCoTab = NULL;
    fSiTab = NULL;

    while (fDphi1 > 360) fDphi1 -= 360;

    MakeTableOfCoSin();
}

//---------------------------------------------------------------------------

AliGPCone::AliGPCone( Text_t *name, Text_t* title,  Float_t *upar, Int_t np) : AliGShape(name, title)
{
    if (upar[2] < 2 ) {
        Error(name, "number of z planes for %s must be at least two !", name);
        return;
    }

    fPhi1  = upar[0];
    fDphi1 = upar[1];
    fNz    = (Int_t) upar[2];
    fNdiv  = 0;
    fRmin  = new Float_t [fNz+1];
    fRmax  = new Float_t [fNz+1];
    fDz    = new Float_t [fNz+1];
    fCoTab = NULL;
    fSiTab = NULL;

    while (fDphi1 > 360) fDphi1 -= 360;

    MakeTableOfCoSin();

    for( int j=3, k=0; k<fNz; k++ ) {
        fDz[k]   = upar[j];
        fRmin[k] = upar[j+1];
        fRmax[k] = upar[j+2];
        j+=3;
    }
}

//-------------------------------------------------------------------------

AliGPCone::AliGPCone( Text_t *name, Text_t* title,  Float_t phi1, Float_t dphi1, Int_t nz, Float_t Z[10], Float_t RMIN[10], Float_t RMAX[10] ) : AliGShape(name, title)
{       
     cout << " ENTRA EN EL CONSTRUCTOR CHUNGO !!!!!! " << endl;
    
     //AliGPCone(name, title,  phi1, dphi1, nz);
     
     //for (int i= 0; i<nz; i++)
         //DefineSection(i, Z[i], RMIN[i], RMAX[i]);
}

//---------------------------------------------------------------------------

AliGPCone::~AliGPCone() 
{
    /* Destructor */
    
    if( fRmin  ) delete [] fRmin;
    if( fRmax  ) delete [] fRmax;
    if( fDz    ) delete [] fDz;
    //if( fSiTab ) delete fSiTab;
    //if( fCoTab ) delete fCoTab;

    fRmin  = NULL;
    fRmax  = NULL;
    fDz    = NULL;
    //fCoTab = NULL;
    //fSiTab = NULL;
}

//-------------------------------------------------------------------------

void AliGPCone::DrawShape(Option_t *option)
{
    Draw(option);
    gPad->Update();
}

//-------------------------------------------------------------------------

void AliGPCone::Draw(Option_t *option)
{
    //cout << " Entra en " << this->GetName() << "::Draw " << endl;

    TString opt = option;
    opt.ToLower();

    if( !gPad ) {
      //TCanvas* Cone = new TCanvas("AliGPCone","AliGPCone",0,0,400,300);
        gPad = new TCanvas("AliGPCone","AliGPCone",0,0,400,300);
        gPad->Range(0,0,1,1);
        gPad->SetFillColor(32); // Light Green
        gPad->SetBorderSize(3);
        gPad->SetBorderMode(0); // -1 (down) 0 (no) 1 (up)
    }
    else {
        if( !opt.Contains("same") ) {
            gPad->Clear();
            gPad->SetName("AliGPCone");
            gPad->SetTitle("AliGPCone");
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
    //cout << " Sale de " << this->GetName() << "::Draw " << endl;
}

//-------------------------------------------------------------------------

void AliGPCone::SetPoints(Float_t *buff)
{
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*Create PCON points*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                            ==================

    Int_t i, j;
    Int_t indx = 0;

    if (buff) {

        Int_t n            = GetNumberOfDivisions()+1;

//*-* We've to check whether the table does exist and create it
//*-* since fCoTab/fSiTab are not saved with any TShape::Streamer function
        if (!fCoTab)   MakeTableOfCoSin();

        for (i = 0; i < fNz; i++)
        {
            for (j = 0; j < n; j++)
            {
                buff[indx++] = fRmin[i] * fCoTab[j];
                buff[indx++] = fRmin[i] * fSiTab[j];
                buff[indx++] = fDz[i];
            }
            for (j = 0; j < n; j++)
            {
                buff[indx++] = fRmax[i] * fCoTab[j];
                buff[indx++] = fRmax[i] * fSiTab[j];
                buff[indx++] = fDz[i];
            }
        }
    }
}

//-------------------------------------------------------------------------

void AliGPCone::Paint(Option_t *option)
{

//*-*-*-*-*-*-*-*Paint this 3-D shape with its current attributes*-*-*-*-*-*-*-*
//*-*            ================================================

    Int_t i, j;
    if (fNz < 2) return;
    const Int_t n = GetNumberOfDivisions()+1;

    SetLineColor( GetCol() );
    
    Int_t numpoints = fNz*2*n;
    if (numpoints <= 0) return;
    //*-* Allocate memory for points *-*

    Float_t *points = new Float_t[3*numpoints];
    if (!points) return;
    SetPoints(points);

    Bool_t rangeView = strcmp(option,"range")==0 ? kTRUE : kFALSE;
    if (!rangeView && gPad->GetView3D()) PaintGLPoints(points);

    //==  for (i = 0; i < numpoints; i++)
    //==          gNode->Local2Master(&points[3*i],&points[3*i]);

    Bool_t specialCase = kFALSE;

    if (fDphi1 == 360)           //mark this as a very special case, when
        specialCase = kTRUE;     //we have to draw this PCON like a TUBE

    X3DBuffer *buff = new X3DBuffer;

    if (buff) {
        buff->numPoints = numpoints;
        buff->numSegs   = 4*(fNz*n-1+(specialCase == kTRUE));
        buff->numPolys  = 2*(fNz*n-1+(specialCase == kTRUE));
    }

    //*-* Allocate memory for points *-*

    buff->points = points;

    Int_t c = ((GetLineColor() % 8) - 1) * 4;     // Basic colors: 0, 1, ... 7
    if (c < 0) c = 0;

    //*-* Allocate memory for segments *-*

    Int_t indx, indx2, k;
    indx = indx2 = 0;

    buff->segs = new Int_t[buff->numSegs*3];
    if (buff->segs) {

        //inside & outside circles, number of segments: 2*fNz*(n-1)
        //             special case number of segments: 2*fNz*n
        for (i = 0; i < fNz*2; i++) {
            indx2 = i*n;
            for (j = 1; j < n; j++) {
                buff->segs[indx++] = c;
                buff->segs[indx++] = indx2+j-1;
                buff->segs[indx++] = indx2+j;
            }
            if (specialCase) {
                buff->segs[indx++] = c;
                buff->segs[indx++] = indx2+j-1;
                buff->segs[indx++] = indx2;
            }
        }

        //bottom & top lines, number of segments: 2*n
        for (i = 0; i < 2; i++) {
            indx2 = i*(fNz-1)*2*n;
            for (j = 0; j < n; j++) {
                buff->segs[indx++] = c;
                buff->segs[indx++] = indx2+j;
                buff->segs[indx++] = indx2+n+j;
            }
        }

        //inside & outside cilindres, number of segments: 2*(fNz-1)*n
        for (i = 0; i < (fNz-1); i++) {

            //inside cilinder
            indx2 = i*n*2;
            for (j = 0; j < n; j++) {
                buff->segs[indx++] = c+2;
                buff->segs[indx++] = indx2+j;
                buff->segs[indx++] = indx2+n*2+j;
            }
            //outside cilinder
            indx2 = i*n*2+n;
            for (j = 0; j < n; j++) {
                buff->segs[indx++] = c+3;
                buff->segs[indx++] = indx2+j;
                buff->segs[indx++] = indx2+n*2+j;
            }
        }

        //left & right sections, number of segments: 2*(fNz-2)
        //          special case number of segments: 0
        if (!specialCase) {
            for (i = 1; i < (fNz-1); i++) {
                for (j = 0; j < 2; j++) {
                    buff->segs[indx++] = c;
                    buff->segs[indx++] =  2*i    * n + j*(n-1);
                    buff->segs[indx++] = (2*i+1) * n + j*(n-1);
                }
            }
        }
    }


    Int_t m = n - 1 + (specialCase == kTRUE);

//*-* Allocate memory for polygons *-*

    indx = 0;

    buff->polys = new Int_t[buff->numPolys*6];

    if (buff->polys) {

        //bottom & top, number of polygons: 2*(n-1)
        // special case number of polygons: 2*n
        for (i = 0; i < 2; i++) {
            for (j = 0; j < n-1; j++) {
                buff->polys[indx++] = c+3;
                buff->polys[indx++] = 4;
                buff->polys[indx++] = 2*fNz*m+i*n+j;
                buff->polys[indx++] = i*(fNz*2-2)*m+m+j;
                buff->polys[indx++] = 2*fNz*m+i*n+j+1;
                buff->polys[indx++] = i*(fNz*2-2)*m+j;
            }
            if (specialCase) {
                buff->polys[indx++] = c+3;
                buff->polys[indx++] = 4;
                buff->polys[indx++] = 2*fNz*m+i*n+j;
                buff->polys[indx++] = i*(fNz*2-2)*m+m+j;
                buff->polys[indx++] = 2*fNz*m+i*n;
                buff->polys[indx++] = i*(fNz*2-2)*m+j;
            }
        }


        //inside & outside, number of polygons: (fNz-1)*2*(n-1)
        for (k = 0; k < (fNz-1); k++) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < n-1; j++) {
                    buff->polys[indx++] = c+i;
                    buff->polys[indx++] = 4;
                    buff->polys[indx++] = (2*k+i*1)*m+j;
                    buff->polys[indx++] = fNz*2*m+(2*k+i*1+2)*n+j;
                    buff->polys[indx++] = (2*k+i*1+2)*m+j;
                    buff->polys[indx++] = fNz*2*m+(2*k+i*1+2)*n+j+1;
                }
                if (specialCase) {
                    buff->polys[indx++] = c+i;
                    buff->polys[indx++] = 4;
                    buff->polys[indx++] = (2*k+i*1)*m+j;
                    buff->polys[indx++] = fNz*2*m+(2*k+i*1+2)*n+j;
                    buff->polys[indx++] = (2*k+i*1+2)*m+j;
                    buff->polys[indx++] = fNz*2*m+(2*k+i*1+2)*n;
                }
            }
        }


        //left & right sections, number of polygons: 2*(fNz-1)
        //          special case number of polygons: 0
        if (!specialCase) {
            indx2 = fNz*2*(n-1);
            for (k = 0; k < (fNz-1); k++) {
                for (i = 0; i < 2; i++) {
                    buff->polys[indx++] = c+2;
                    buff->polys[indx++] = 4;
                    buff->polys[indx++] = k==0 ? indx2+i*(n-1) : indx2+2*fNz*n+2*(k-1)+i;
                    buff->polys[indx++] = indx2+2*(k+1)*n+i*(n-1);
                    buff->polys[indx++] = indx2+2*fNz*n+2*k+i;
                    buff->polys[indx++] = indx2+(2*k+3)*n+i*(n-1);
                }
            }
            buff->polys[indx-8] = indx2+n;
            buff->polys[indx-2] = indx2+2*n-1;
        }
    }

    //*-* Paint in the pad
    
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
    
//------------------------------------------------------------------------- 

void AliGPCone::MakeTableOfCoSin()
{
    const Double_t PI  = TMath::ATan(1) * 4.0;
    const Double_t ragrad  = PI/180.0;

    Int_t j;
    Int_t n = GetNumberOfDivisions () + 1;
    if (fCoTab)
        delete [] fCoTab; // Delete the old tab if any
        fCoTab = new Double_t [n];
    if (!fCoTab ) return;

    if (fSiTab)
        delete [] fSiTab; // Delete the old tab if any
    fSiTab = new Double_t [n];
    if (!fSiTab ) return;

    Double_t range   = Double_t(fDphi1 * ragrad);
    Double_t phi1    = Double_t(fPhi1  * ragrad);
    Double_t angstep = range/(n-1);

    Double_t ph = phi1;
    for (j = 0; j < n; j++)
    {
        ph = phi1 + j*angstep;
        fCoTab[j] = TMath::Cos(ph);
        fSiTab[j] = TMath::Sin(ph);
    }

}
 
//----------------------------------------------------    
    
void AliGPCone::SetNumberOfDivisions (Int_t p)
{
    if (GetNumberOfDivisions () == p) return;
    fNdiv=p;
    MakeTableOfCoSin();
}   
    
//-------------------------------------------------------   

void AliGPCone::DefineSection(Int_t secNum, Float_t z, Float_t rmin, Float_t rmax )
{

//*-*-*-*-*-*-*-*-*-*Defines section secNum of the polycone*-*-*-*-*-*-*-*-*-*-*
//*-*                ======================================
//
//     - rmin  radius of the inner circle in the cross-section
//
//     - rmax  radius of the outer circle in the cross-section
//
//     - z     z coordinate of the section

    if ((secNum < 0) || (secNum >= fNz)) return;

    fRmin[secNum] = rmin;
    fRmax[secNum] = rmax;
    fDz[secNum]   = z;
}

//-------------------------------------------------------   

void AliGPCone::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliGPCone.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliGShape::Streamer(R__b);
      R__b.ReadArray(fSiTab); //
      R__b.ReadArray(fCoTab); //
      R__b >> fPhi1;
      R__b >> fDphi1;
      R__b >> fNz;
      R__b.ReadArray(fRmin); //
      R__b.ReadArray(fRmax); //
      R__b.ReadArray(fDz); //
      R__b >> fNdiv;
   } else {
      R__b.WriteVersion(AliGPCone::IsA());
      AliGShape::Streamer(R__b);
      R__b.WriteArray(fSiTab, GetNumberOfDivisions()+1); //
      R__b.WriteArray(fCoTab, GetNumberOfDivisions()+1); //
      R__b << fPhi1;
      R__b << fDphi1;
      R__b << fNz;
      R__b.WriteArray(fRmin, fNz+1); //
      R__b.WriteArray(fRmax, fNz+1); //
      R__b.WriteArray(fDz, fNz+1); //
      R__b << fNdiv;
   }
}

//-------------------------------------------------------   


