// -*- C++ -*-
// 
// 1998/10/19
// ---------------------------------------------------------------------------
//
// AliGCone Class
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
#include "AliGCone.h"
#include "TROOT.h"

ClassImp(AliGCone)

AliGCone::AliGCone() : AliGTube()
{
    /* Default Constructor */
    fRmax2 = 0.;      // outside radius at the high z limit
    fRmin2 = 0.;      // inside radius at the high z limit

    //SetLineColor(5);  // Yellow
}

//-------------------------------------------------------------------------

AliGCone::AliGCone(Text_t *name, Text_t *title, Float_t dz, Float_t rmin1, Float_t rmax1, Float_t rmin2, Float_t rmax2) : AliGTube(name, title, rmin1, rmax1, dz)
{
    /* Constructor */
    fRmax2 = rmax2;
    fRmin2 = rmin2;

    //SetLineColor(5); // Yellow
}

//-------------------------------------------------------------------------

AliGCone::AliGCone(AliGCone *cone) 
{
    /* Copy Constructor */
    fRmax2 = cone->fRmax2;
    fRmin2 = cone->fRmin2;
    fColor = cone->fColor;
    
    //SetLineColor(5); // Yellow
}

//-------------------------------------------------------------------------

AliGCone::AliGCone(Text_t *name, Text_t *title, Float_t dz, Float_t rmax1, Float_t rmax2) : AliGTube(name, title, 0, rmax1, dz)
{
    /* Simplified Constructor */
    fRmin2 = 0;
    fRmax2 = rmax2;

    //SetLineColor(5); // Yellow
}

//-------------------------------------------------------------------------

AliGCone::~AliGCone() {
    /* Destructor */
}

//-------------------------------------------------------------------------

void AliGCone::DrawShape(Option_t *option)
{
    Draw(option);
    gPad->Update();
}

//-------------------------------------------------------------------------

void AliGCone::Draw(Option_t *option)
{
    cout << " Entra en " << this->GetName() << "::Draw " << endl;
    TString opt = option;
    opt.ToLower();

    if( !gPad ) {
      //TCanvas* Cone = new TCanvas("AliGCone","AliGCone",0,0,400,300);
        gPad = new TCanvas("AliGCone","AliGCone",0,0,400,300);
        gPad->Range(0,0,1,1);
        gPad->SetFillColor(32); // Light Green
        gPad->SetBorderSize(3);
        gPad->SetBorderMode(0); // -1 (down) 0 (no) 1 (up)
    }
    else {
        if( !opt.Contains("same") ) {
            gPad->Clear();
            gPad->SetName("AliGCone");
            gPad->SetTitle("AliGCone");
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

void AliGCone::SetPoints(Float_t *buff)
{
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*Create CONE points*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                            ==================

    //cout << " Entra en " << this->GetName() << "::Paint " << endl;
    SetLineColor( GetCol() );

    Float_t rmin1, rmax1, dz;
    Int_t j, n;

    n = GetNumberOfDivisions();

    rmin1 = AliGTube::fRmin;
    rmax1 = AliGTube::fRmax;
    dz    = AliGTube::fDz;

    Int_t indx = 0;

//*-* We've to checxk whether the table does exist and create it
//*-* since fCoTab/fSiTab are not saved with any TShape::Streamer function
    if (!fCoTab)   MakeTableOfCoSin();

    if (buff) {
        for (j = 0; j < n; j++) {
            buff[indx++] = rmin1 * fCoTab[j];
            buff[indx++] = rmin1 * fSiTab[j];
            buff[indx++] = -dz;
        }
        for (j = 0; j < n; j++) {

            buff[indx++] = rmax1 * fCoTab[j];
            buff[indx++] = rmax1 * fSiTab[j];
            buff[indx++] = -dz;
        }

        for (j = 0; j < n; j++) {
            buff[indx++] = fRmin2 * fCoTab[j];
            buff[indx++] = fRmin2 * fSiTab[j];
            buff[indx++] = dz;
        }

        for (j = 0; j < n; j++) {
            buff[indx++] = fRmax2 * fCoTab[j];
            buff[indx++] = fRmax2 * fSiTab[j];
            buff[indx++] = dz;
        }
    }

    //cout << " Sale de " << this->GetName() << "::Paint " << endl;   
}

//-------------------------------------------------------------------------


