   // -*- C++ -*-
// 
// 1998/07/23
// ---------------------------------------------------------------------------
//
// AliGSphere Class
//
// This file is part of the ALICE Geometry Database .
//
// Author:  Joana E. Santo
//
// ---------------------------------------------------------------------------

#include <TView.h>
#include <TVirtualPad.h>
#include <iostream.h>
#include <TCanvas.h>
#include <TGLKernelABC.h>
#include "TROOT.h"
#include "AliGSphere.h"


ClassImp(AliGSphere)

AliGSphere::AliGSphere()
{
    /* Default Constructor */
    faX          = 1.;  // Coeff along Ox
    faY          = 1.;  // Coeff along Oy
    faZ          = 1.;  // Coeff along Oz
    fAspectRatio = 1.;  // Relation between asumth and grid size (by default 1.0)
    fCoTab       = NULL;// Table of cos(fPhimin) .... cos(Phi)
    fCoThetaTab  = NULL;// Table of sin(gThemin) .... cos(Theta)
    fName        = "";
    fNdiv        = 0;   // number of divisions
    fNz          = 0;   // number of sections
    fPhimax      = 0.;  // maximum phi
    fPhimin      = 0.;  // minimum phi
    fRmax        = 0.;  // maximum radius
    fRmin        = 0.;  // minimum radius
    fSiTab       = NULL;// Table of sin(fPhimin) .... sin(Phi)
    fThemax      = 0.;  // maximum theta
    fThemin      = 0.;  // minimum theta
    fTitle       = "";
}

// ---------------------------------------------------------------------------

AliGSphere::AliGSphere(Text_t *name, Text_t *title, Float_t rmin, Float_t rmax, Float_t themin, Float_t themax, Float_t phimin, Float_t phimax) : AliGShape(name, title)
{
    /* Constructor */
    faX          = 1.;    // Coeff along Ox
    faY          = 1.;    // Coeff along Oy
    faZ          = 1.;    // Coeff along Oz
    fAspectRatio = 1.;    // Relation between asumth and grid size (by default 1.0)
    fCoTab       = NULL;  // Table of cos(fPhimin) .... cos(Phi)
    fCoThetaTab  = NULL;  // Table of sin(gThemin) .... cos(Theta)
    fNdiv        = 0;     // number of divisions
    fNz          = 0;     // number of sections
    fPhimax      = phimax;// maximum phi
    fPhimin      = phimin;// minimum phi
    fRmax        = rmax;  // maximum radius
    fRmin        = rmin;  // minimum radius
    fSiTab       = NULL;  // Table of sin(fPhimin) .... sin(Phi)
    fThemax      = themax;// maximum theta
    fThemin      = themin;// minimum theta

    SetNumberOfDivisions (20);
}

// ---------------------------------------------------------------------------

AliGSphere::AliGSphere(AliGSphere *sphere) 
{
    /* Copy Constructor */
    faX          = 1.;    // Coeff along Ox
    faY          = 1.;    // Coeff along Oy
    faZ          = 1.;    // Coeff along Oz
    fAspectRatio = 1.;    // Relation between asumth and grid size (by default 1.0)
    fColor       = sphere->fColor;
    fCoTab       = NULL;  // Table of cos(fPhimin) .... cos(Phi)
    fCoThetaTab  = NULL;  // Table of sin(gThemin) .... cos(Theta)
    fNdiv        = 0;     // number of divisions
    fNz          = 0;     // number of sections
    fPhimax      = sphere->fPhimax;// maximum phi
    fPhimin      = sphere->fPhimin;// minimum phi
    fRmax        = sphere->fRmax;  // maximum radius
    fRmin        = sphere->fRmin;  // minimum radius
    fSiTab       = NULL;  // Table of sin(fPhimin) .... sin(Phi)
    fThemax      = sphere->fThemax;// maximum theta
    fThemin      = sphere->fThemin;// minimum theta

    SetNumberOfDivisions (20);
}
// ---------------------------------------------------------------------------

AliGSphere::AliGSphere(Text_t *name, Text_t *title, Float_t rmax) : AliGShape(name, title)
{
    /* Simplified Constructor */
    faX          = 1.;   // Coeff along Ox
    faY          = 1.;   // Coeff along Oy
    faZ          = 1.;   // Coeff along Oz
    fAspectRatio = 1.;   // Relation between asumth and grid size (by default 1.0)
    fCoTab       = NULL; // Table of cos(fPhimin) .... cos(Phi)
    fCoThetaTab  = NULL; // Table of sin(gThemin) .... cos(Theta)
    fNdiv        = 0;    // number of divisions
    fNz          = 0;    // number of sections
    fPhimax      = 360.; // maximum phi
    fPhimin      = 0.;   // minimum phi
    fRmax        = rmax; // maximum radius
    fRmin        = 0.;   // minimum radius
    fSiTab       = NULL; // Table of sin(fPhimin) .... sin(Phi)
    fThemax      = 180.; // maximum theta
    fThemin      = 0.;   // minimum theta

    SetNumberOfDivisions (20);
}

// ---------------------------------------------------------------------------

AliGSphere::~AliGSphere() {
    /* Destructor */
    delete [] fCoThetaTab;  // Table of sin(gThemin) .... cos(Theta)
    delete [] fCoTab;
    delete [] fSiTab;
}

// ---------------------------------------------------------------------------

void AliGSphere::DrawShape(Option_t *option)
{
    Draw(option);
    gPad->Update();
}

//-------------------------------------------------------------------------

//void AliGSphere::Draw()
void AliGSphere::Draw(Option_t *option)
{
    TString opt = option;
    opt.ToLower();

    if( !gPad ) {
        gPad = new TCanvas("AliGSphere","AliGSphere",0,0,400,300);
        gPad->Range(0,0,1,1);
        gPad->SetFillColor(32); // Light Green
        gPad->SetBorderSize(3);
        gPad->SetBorderMode(0); // -1 (down) 0 (no) 1 (up)
    }
    else {
        if( !opt.Contains("same") ) {
            gPad->Clear();
            gPad->SetName("AliGSphere");
            gPad->SetTitle("AliGSphere");
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
}

// ---------------------------------------------------------------------------

void AliGSphere::Paint(Option_t *option)
{
//*-*-*-*-*-*-*-*Paint this 3-D shape with its current attributes*-*-*-*-*-*-*-*
//*-*            ================================================

    SetLineColor( GetCol() );

    Int_t i, j;
    const Int_t n = GetNumberOfDivisions()+1;
    Int_t nz = fNz+1;
    Int_t numpoints = 2*n*nz;
    if (nz < 2) return;

    if (numpoints <= 0) return;
    //*-* Allocate memory for points *-*

    Float_t *points = new Float_t[3*numpoints];
    if (!points) return;
    SetPoints(points);

    if (gPad->GetView3D()) PaintGLPoints(points);

 //==  for (i = 0; i < numpoints; i++)
 //==          gNode->Local2Master(&points[3*i],&points[3*i]);

    Bool_t specialCase = kFALSE;

    if (TMath::Abs(TMath::Sin(2*(fPhimax - fPhimin))) <= 0.01)  //mark this as a very special case, when
        specialCase = kTRUE;                                  //we have to draw this PCON like a TUBE

    X3DBuffer *buff = new X3DBuffer;
    if (buff) {
        buff->numPoints = numpoints;
        buff->numSegs   = 4*(nz*n-1+(specialCase == kTRUE));
        buff->numPolys  = 2*(nz*n-1+(specialCase == kTRUE));
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

        //inside & outside spheres, number of segments: 2*nz*(n-1)
        //             special case number of segments: 2*nz*n
        for (i = 0; i < nz*2; i++) {
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
            indx2 = i*(nz-1)*2*n;
            for (j = 0; j < n; j++) {
                buff->segs[indx++] = c;
                buff->segs[indx++] = indx2+j;
                buff->segs[indx++] = indx2+n+j;
            }
        }

        //inside & outside spheres, number of segments: 2*(nz-1)*n
        for (i = 0; i < (nz-1); i++) {

            //inside sphere
            indx2 = i*n*2;
            for (j = 0; j < n; j++) {
                buff->segs[indx++] = c+2;
                buff->segs[indx++] = indx2+j;
                buff->segs[indx++] = indx2+n*2+j;
            }
            //outside sphere
            indx2 = i*n*2+n;
            for (j = 0; j < n; j++) {
                buff->segs[indx++] = c+3;
                buff->segs[indx++] = indx2+j;
                buff->segs[indx++] = indx2+n*2+j;
            }
        }

        //left & right sections, number of segments: 2*(nz-2)
        //          special case number of segments: 0
        /*if (!specialCase) {
            for (i = 1; i < (nz-1); i++) {
                for (j = 0; j < 2; j++) {
                    buff->segs[indx++] = c;
                    buff->segs[indx++] =  2*i    * n + j*(n-1);
                    buff->segs[indx++] = (2*i+1) * n + j*(n-1);
                }
            }
        }*/
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
                buff->polys[indx++] = 2*nz*m+i*n+j;
                buff->polys[indx++] = i*(nz*2-2)*m+m+j;
                buff->polys[indx++] = 2*nz*m+i*n+j+1;
                buff->polys[indx++] = i*(nz*2-2)*m+j;
            }
            if (specialCase) {
                buff->polys[indx++] = c+3;
                buff->polys[indx++] = 4;
                buff->polys[indx++] = 2*nz*m+i*n+j;
                buff->polys[indx++] = i*(nz*2-2)*m+m+j;
                buff->polys[indx++] = 2*nz*m+i*n;
                buff->polys[indx++] = i*(nz*2-2)*m+j;
            }
        }


        //inside & outside, number of polygons: (nz-1)*2*(n-1)
        for (k = 0; k < (nz-1); k++) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < n-1; j++) {
                    buff->polys[indx++] = c+i;
                    buff->polys[indx++] = 4;
                    buff->polys[indx++] = (2*k+i*1)*m+j;
                    buff->polys[indx++] = nz*2*m+(2*k+i*1+2)*n+j;
                    buff->polys[indx++] = (2*k+i*1+2)*m+j;
                    buff->polys[indx++] = nz*2*m+(2*k+i*1+2)*n+j+1;
                }
                if (specialCase) {
                    buff->polys[indx++] = c+i;
                    buff->polys[indx++] = 4;
                    buff->polys[indx++] = (2*k+i*1)*m+j;
                    buff->polys[indx++] = nz*2*m+(2*k+i*1+2)*n+j;
                    buff->polys[indx++] = (2*k+i*1+2)*m+j;
                    buff->polys[indx++] = nz*2*m+(2*k+i*1+2)*n;
                }
            }
        }

        //left & right sections, number of polygons: 2*(nz-1)
        //          special case number of polygons: 0
        /*if (!specialCase) {
            indx2 = nz*2*(n-1);
            for (k = 0; k < (nz-1); k++) {
                for (i = 0; i < 2; i++) {
                    buff->polys[indx++] = c+2;
                    buff->polys[indx++] = 4;
                    buff->polys[indx++] = k==0 ? indx2+i*(n-1) : indx2+2*nz*n+2*(k-1)+i;
                    buff->polys[indx++] = indx2+2*(k+1)*n+i*(n-1);
                    buff->polys[indx++] = indx2+2*nz*n+2*k+i;
                    buff->polys[indx++] = indx2+(2*k+3)*n+i*(n-1);
                }
            }
            buff->polys[indx-8] = indx2+n;
            buff->polys[indx-2] = indx2+2*n-1;
        }*/
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

    delete [] points;
    if (buff->segs)     delete [] buff->segs;
    if (buff->polys)    delete [] buff->polys;
    if (buff)           delete    buff;

    

}

// ---------------------------------------------------------------------------

void AliGSphere::SetEllipse(Float_t *factors)
{
  if (factors[0] > 0) faX = factors[0];
  if (factors[1] > 0) faY = factors[1];
  if (factors[2] > 0) faZ = factors[2];
  MakeTableOfCoSin();
}

// ---------------------------------------------------------------------------

void AliGSphere::SetNumberOfDivisions (Int_t p)
{
    if (GetNumberOfDivisions () == p)
        return;
    fNdiv = p;
    fNz   = Int_t(fAspectRatio*fNdiv*(fThemax - fThemin )/(fPhimax - fPhimin )) + 1;
    MakeTableOfCoSin();
}

// ---------------------------------------------------------------------------

void AliGSphere::SetPoints(Float_t *buff)
{
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*Create SPHE points*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                            ==================
    Int_t i, j;
    Int_t indx = 0;

    if (buff) {
        Int_t n = GetNumberOfDivisions()+1;

//*-* We've to check whether the table does exist and create it
//*-* since fCoTab/fSiTab are not saved with any TShape::Streamer function
        if (!fCoTab)   MakeTableOfCoSin();

        Float_t z;
        for (i = 0; i < fNz+1; i++)
        {
            z = fRmin * fCoThetaTab[i]; // fSinPhiTab[i];
            Float_t sithet = TMath::Sqrt(TMath::Abs(1-fCoThetaTab[i]*fCoThetaTab[i]));
            Float_t zi = fRmin*sithet;
            for (j = 0; j < n; j++)
            {
                buff[indx++] = zi * fCoTab[j];
                buff[indx++] = zi * fSiTab[j];
                buff[indx++] = z;
            }
            z = fRmax * fCoThetaTab[i];
            zi = fRmax*sithet;
            for (j = 0; j < n; j++)
            {
                buff[indx++] = zi * fCoTab[j];
                buff[indx++] = zi * fSiTab[j];
                buff[indx++] = z;
            }
        }
    }
}

// ---------------------------------------------------------------------------

void AliGSphere::MakeTableOfCoSin()
{
    const Double_t PI  = TMath::ATan(1) * 4.0;
    const Double_t ragrad  = PI/180.0;

    Float_t dphi = fPhimax - fPhimin;
    while (dphi > 360) dphi -= 360;

    Float_t dtet = fThemax - fThemin;
    while (dtet > 180) dtet -= 180;

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

    Double_t range   = Double_t(dphi * ragrad);
    Double_t phi1    = Double_t(fPhimin  * ragrad);
    Double_t angstep = range/(n-1);

    Double_t ph = phi1;
    for (j = 0; j < n; j++)
    {
        ph = phi1 + j*angstep;
        fCoTab[j] = TMath::Cos(ph);
        fSiTab[j] = TMath::Sin(ph);
    }

    n  = fNz + 1;

    if (fCoThetaTab)
        delete [] fCoThetaTab; // Delete the old tab if any
    fCoThetaTab = new Double_t [n];
    if (!fCoThetaTab ) return;

    range   = Double_t(dtet * ragrad);
    phi1    = Double_t(fThemin  * ragrad);
    angstep = range/(n-1);

    ph = phi1;
    for (j = 0; j < n; j++)
    {
        fCoThetaTab[n-j-1] = TMath::Cos(ph);
        ph += angstep;
    }

}

// ---------------------------------------------------------------------------

void AliGSphere::PaintGLPoints(Float_t *vertex)
{
    gGLKernel->PaintCone(vertex,-(GetNumberOfDivisions()+1),fNz+1);
}

// ---------------------------------------------------------------------------

void AliGSphere::Sizeof3D() const
{
//*-*-*-*-*-*-*Return total X3D size of this shape with its attributes*-*-*-*-*-*
//*-*          =======================================================

    cout << " Entra en AliGSphere::Sizeof3D() " << endl;

    Int_t n;

    n = GetNumberOfDivisions()+1;
    Int_t nz = fNz+1;

    //cout << " n = " << n << "   y   nz = " << nz << endl;

    Bool_t specialCase = kFALSE;

    if (TMath::Abs(TMath::Sin(2*(fPhimax - fPhimin))) <= 0.01)  //mark this as a very special case, when
          specialCase = kTRUE;                                  //we have to draw this PCON like a TUBE

    gSize3D.numPoints += 2*n*nz;
    gSize3D.numSegs   += 4*(nz*n-1+(specialCase == kTRUE));
    gSize3D.numPolys  += 2*(nz*n-1+(specialCase == kTRUE));
}

// ---------------------------------------------------------------------------

void AliGSphere::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliGSphere.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliGShape::Streamer(R__b);
      R__b >> fAspectRatio;
      R__b.ReadArray(fCoTab); //
      R__b.ReadArray(fCoThetaTab); //
      R__b >> fNdiv;
      R__b >> fNz;
      R__b.ReadArray(fSiTab); //
      R__b >> faX;
      R__b >> faY;
      R__b >> faZ;
      R__b >> fPhimax;
      R__b >> fPhimin;
      R__b >> fRmax;
      R__b >> fRmin;
      R__b >> fThemax;
      R__b >> fThemin;
   } else {
      R__b.WriteVersion(AliGSphere::IsA());
      AliGShape::Streamer(R__b);
      R__b << fAspectRatio;
      R__b.WriteArray(fCoTab, GetNumberOfDivisions()+1); //
      R__b.WriteArray(fCoThetaTab, fNz+1); //
      R__b << fNdiv;
      R__b << fNz;
      R__b.WriteArray(fSiTab, GetNumberOfDivisions()+1); //
      R__b << faX;
      R__b << faY;
      R__b << faZ;
      R__b << fPhimax;
      R__b << fPhimin;
      R__b << fRmax;
      R__b << fRmin;
      R__b << fThemax;
      R__b << fThemin;
   }
}

// ---------------------------------------------------------------------------

