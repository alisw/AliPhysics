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
Revision 1.9  2001/01/23 18:58:19  hristov
Initialisation of some pointers

Revision 1.8  2000/10/06 09:09:01  morsch
Pad colour according to z-position (slats).

Revision 1.7  2000/10/02 21:28:09  fca
Removal of useless dependecies via forward declarations

Revision 1.6  2000/07/03 11:54:57  morsch
AliMUONSegmentation and AliMUONHitMap have been replaced by AliSegmentation and AliHitMap in STEER
The methods GetPadIxy and GetPadXxy of AliMUONSegmentation have changed name to GetPadI and GetPadC.

Revision 1.5  2000/06/28 15:16:35  morsch
(1) Client code adapted to new method signatures in AliMUONSegmentation (see comments there)
to allow development of slat-muon chamber simulation and reconstruction code in the MUON
framework. The changes should have no side effects (mostly dummy arguments).
(2) Hit disintegration uses 3-dim hit coordinates to allow simulation
of chambers with overlapping modules (MakePadHits, Disintegration).

Revision 1.4  2000/06/27 09:46:57  morsch
kMAXZOOM global constant now in AliMUONConstants

Revision 1.3  2000/06/26 14:02:38  morsch
Add class AliMUONConstants with MUON specific constants using static memeber data and access methods.

Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.15  2000/06/14 14:37:53  morsch
method Trigger() modified

Revision 1.1.2.14  2000/06/09 21:57:09  morsch
Bug in color scale diplay corrected.
Most coding rule violations corrected.

Revision 1.1.2.13  2000/05/02 11:57:27  morsch
Coding rules RN3, RN13, RN17 violations corrected.

Revision 1.1.2.12  2000/04/26 12:27:33  morsch
Mods for trigger display (P. Crochet):
- color code versus time for pad hits in trigger chambers
- call to TriggerDecision corrected

Revision 1.1.2.11  2000/04/26 09:04:46  morsch
Obsolete cathode correlation related code removed.

Revision 1.1.2.10  2000/04/19 19:43:47  morsch
change NCH to kNCH as in AliMUON.h
no more TreeC related methods

Revision 1.1.2.9  2000/03/20 18:10:33  morsch
Trigger method for "online" trigger decission added

Revision 1.1.2.8  2000/02/23 10:12:01  morsch
Dont't try to draw reconstructed hit coordinates for Trigger Chambers.
General clean-up of commented code.

Revision 1.1.2.7  2000/02/17 14:36:55  morsch
Display of Trigger hits and clusters added.
Displacement between clusters and hits has to be investigated and corrected ! (A.M.)

Revision 1.1.2.6  2000/02/15 10:19:42  morsch
Previous log messages included

Revision 1.1.2.5  2000/02/15 10:09:09  morsch
Log Message added

Revision 1.1.2.4  2000/02/08 09:17:16  gosset    
One more improvement of MUON display:
same zoom for both cathode planes in a chamber

Revision 1.1.2.3  2000/02/07 15:37:21  gosset
A few improvements of the MUON display:
new buttons to change either chamber or cathode,
added to the more complicated way
(right mouse click and explicit filling of chamber and cathode)

Revision 1.1.2.2  2000/02/04 10:57:34  gosset
Z position of the chambers:
it was the Z position of the stations;
it is now really the Z position of the chambers.
   !!!! WARNING: THE CALLS TO "AliMUONChamber::SetZPOS"
   !!!!                   AND "AliMUONChamber::ZPosition"
   !!!! HAVE TO BE CHANGED TO "AliMUONChamber::"SetZ"
   !!!!                   AND "AliMUONChamber::Z"                                                           
*/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliDisplay                                                           //
//                                                                      //
// Utility class to display ALICE outline, tracks, hits,..              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TTree.h>
#include <TButton.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TView.h>
#include <TText.h>
#include <TPolyMarker3D.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TList.h>
#include <TDiamond.h>
#include <TNode.h>
#include <TArc.h>
#include <TTUBE.h>
#include <TSlider.h>
#include <TSliderBox.h>
#include <TGaxis.h>
#include <TVirtualX.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TGeometry.h>
#include <X3DBuffer.h>
#include <TMarker3DBox.h>

#include "AliRun.h"
#include "AliDetector.h"
#include "AliMUON.h"
#include "AliMUONDisplay.h"
#include "AliMUONPoints.h"
#include "TParticle.h"
#include "AliMUONTriggerDecision.h"

#include "AliMUONHit.h"
#include "AliMUONPadHit.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"

#include "AliSegmentation.h"
#include "AliMUONResponse.h"
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
// to manage the same zoom on both cathodes



ClassImp(AliMUONDisplay)


//_____________________________________________________________________________
AliMUONDisplay::AliMUONDisplay()
{
// Constructor
    fPoints = 0;
    fPhits = 0;
    fRpoints = 0;
    fR2points = 0;
    fCpoints = 0;
    fCanvas = 0;
    fNextCathode = kFALSE; 
    fColPad = 0;
}

//_____________________________________________________________________________
AliMUONDisplay::AliMUONDisplay(Int_t size)
{
// Create an event display object.
// A canvas named "edisplay" is created with a vertical size in pixels
//
//    A QUICK Overview of the Event Display functions
//    ===============================================
//
//  The event display can ve invoked by executing the macro "display.C"
// A canvas like in the picture below will appear.
//
//  On the left side of the canvas, the following buttons appear:
//   *Next*       to move to the next event
//   *Previous*   to move to the previous event

//   *Pick*       Select this option to be able to point on a track with the
//                mouse. Once on the track, use the right button to select
//                an action. For example, select SetMarkerAttributes to
//                change the marker type/color/size for the track.
//   *Zoom*       Select this option (default) if you want to zoom.
//                To zoom, simply select the selected area with the left button.
//   *UnZoom*     To revert to the previous picture size.
//
//   slider R     On the left side, the vertical slider can be used to
//                set the default picture size.
//
//    When you are in Zoom mode, you can click on the black part of the canvas
//  to select special options with the right mouse button.

//
//  When you are in pick mode, you can "Inspect" the object pointed by the mouse.
//  When you are on a track, select the menu item "InspectParticle"
//  to display the current particle attributes.
//
//  You can activate the Root browser by selecting the Inspect menu
//  in the canvas tool bar menu. Then select "Start Browser"
//  This will open a new canvas with the browser. At this point, you may want
//  to display some histograms (from the Trees). Go to the "File" menu
//  of the browser and click on "New canvas".
//  In the browser, click on item "ROOT files" in the left pane.
//  Click on galice.root.
//  Click on TH
//  Click on TPC for example
//  Click on any variable (eg TPC.fX) to histogram the variable.
//
//   If you are lost, you can click on HELP in any Root canvas or browser.
//Begin_Html
/*
<img src="gif/AliMUONDisplay.gif">
*/
//End_Html


    fPad = 0;
    
    gAlice->SetDisplay(this);
   
   // Initialize display default parameters
    SetRange(200,2000);
   // Set front view by default
    fTheta =   0;
    fPhi   = -90;
    fPsi   =   0;
    fChamber = 1;
    fCathode = 1;
    //   fRzone   = 1.e10;
    fDrawClusters  = kTRUE;
    fDrawCoG       = kTRUE;
    fDrawCoG  = kTRUE;
    fZoomMode      = 1;
    fZooms         = 0;
    fClustersCuts  = 0;
    fPoints        = 0;
    fPhits         = 0;
    fRpoints       = 0;
    fR2points = 0;
    fCpoints = 0;
    // Create colors
    CreateColors();
    // Create display canvas
    Int_t ysize = size;
    if (ysize < 100) ysize = 750;
    Int_t xsize = Int_t(size*830./ysize);
    fCanvas = new TCanvas("Canvas", "MUON Clusters Display",14,47,xsize,ysize);
    fCanvas->ToggleEventStatus();
    
   // Create main display pad
    fPad = new TPad("viewpad", "MUON display",0.15,0,0.9,1);
    fPad->Draw();
    fPad->Modified();
    fPad->SetFillColor(30);
    fPad->SetBorderSize(2);

    fCanvas->cd();

   // Create colors pad
    fColPad = new TPad("colpad", "Colors pad",0.9,0,1,1);
    fColPad->Draw();
    fColPad->SetFillColor(17);
    fColPad->SetBorderSize(2);
    fColPad->cd();
    DisplayColorScale();

    fCanvas->cd();
   // Create user interface control pad
    DisplayButtons();
    fCanvas->cd();

   // Create Range and mode pad
    Float_t dxtr     = 0.15;
    Float_t dytr     = 0.45;
    fTrigPad = new TPad("trigger", "range and mode pad",0,0,dxtr,dytr);
    fTrigPad->SetEditable(kFALSE);
    fTrigPad->Draw();
    fTrigPad->cd();
    fTrigPad->SetFillColor(22);
    fTrigPad->SetBorderSize(2);
    fRangeSlider = new TSlider("range","range",0.7,0.42,0.9,0.98);
    fRangeSlider->SetObject(this);
    char pickmode[] = "gAlice->Display()->SetPickMode()";
    Float_t db = 0.09;
    fPickButton = new TButton("Pick",pickmode,0.05,0.32,0.65,0.32+db);
    fPickButton->SetFillColor(38);
    fPickButton->Draw();
    char zoommode[] = "gAlice->Display()->SetZoomMode()";
    fZoomButton = new TButton("Zoom",zoommode,0.05,0.21,0.65,0.21+db);
    fZoomButton->SetFillColor(38);
    fZoomButton->Draw();
    fArcButton = new TArc(.8,fZoomButton->GetYlowNDC()+0.5*db,0.33*db);
    fArcButton->SetFillColor(kGreen);
    fArcButton->Draw();
    char butUnzoom[] = "gAlice->Display()->UnZoom()";
    TButton *button = new TButton("UnZoom",butUnzoom,0.05,0.05,0.95,0.15);
    button->SetFillColor(38);
    button->Draw();
    AppendPad(); // append display object as last object to force selection
    
    fCanvas->cd();
    fCanvas->Update();
    fNextCathode = kFALSE; 
}

AliMUONDisplay::AliMUONDisplay(const AliMUONDisplay & display)
{
// Dummy copy constructor    
    ;
}



//_____________________________________________________________________________
AliMUONDisplay::~AliMUONDisplay()
{
  // Delete space point structure
    if (fPoints) fPoints->Delete();
    delete fPoints;
    fPoints     = 0;
    //
    if (fPhits) fPhits->Delete();
    delete fPhits;
    fPhits     = 0;
    //
    if (fRpoints) fRpoints->Delete();
    delete fRpoints;
    fRpoints     = 0;
//
    if (fR2points) fR2points->Delete();
    delete fR2points;
    fR2points     = 0;
//
    if (fCpoints) fCpoints->Delete();
    delete fCpoints;
    fCpoints     = 0;
}

//_____________________________________________________________________________
void AliMUONDisplay::Clear(Option_t *)
{
//    Delete graphics temporary objects
}

//_____________________________________________________________________________
void AliMUONDisplay::DisplayButtons()
{
//    Create the user interface buttons


    fButtons = new TPad("buttons", "newpad",0,0.45,0.15,1);
    fButtons->SetEditable(kFALSE);
    fButtons->Draw();
    fButtons->SetFillColor(38);
    fButtons->SetBorderSize(2);
    fButtons->cd();
    
//   Int_t butcolor = 33;
    Float_t dbutton = 0.08;
    Float_t y  = 0.96;
    Float_t dy = 0.014;
    Float_t x0 = 0.05;
    Float_t x1 = 0.95;
    
    TButton *button;
    char but1[] = "gAlice->Display()->ShowNextEvent(1)";
    button = new TButton("Event +", but1, x0, y - dbutton, x1, y);
    button->SetFillColor(38);
    button->Draw();
    
    y -= dbutton + dy;
    char but2[] = "gAlice->Display()->ShowNextEvent(-1)";
    button = new TButton("Event -", but2, x0, y - dbutton, x1, y);
    button->SetFillColor(38);
    button->Draw();
   
    y -= dbutton + dy;
    char but3[] = "((AliMUONDisplay*)(gAlice->Display()))->NextChamber(1)";
    button = new TButton("Chamber +", but3, x0, y - dbutton, x1, y);
    button->SetFillColor(38);
    button->Draw();
    
    y -= dbutton + dy;
    char but4[] = "((AliMUONDisplay*)(gAlice->Display()))->NextChamber(-1)";
    button = new TButton("Chamber -", but4, x0, y - dbutton, x1, y);
    button->SetFillColor(38);
    button->Draw();
    
    y -= dbutton + dy;
    char but5[] = "((AliMUONDisplay*)(gAlice->Display()))->SetChamberAndCathode(1,1)";
    button = new TButton("Chamber 1", but5, x0, y - dbutton, x1, y);
    button->SetFillColor(38);
    button->Draw();
   
    y -= dbutton + dy;
    char but6[] = "((AliMUONDisplay*)(gAlice->Display()))->NextCathode()";
    button = new TButton("Cathode <>", but6, x0, y - dbutton, x1, y);
    button->SetFillColor(38);
    button->Draw();

    y -= dbutton + dy;
    char but7[] = "((AliMUONDisplay*)(gAlice->Display()))->Trigger()";
    button = new TButton("Trigger", but7, x0, y - dbutton, x1, y);
    button->SetFillColor(38);
    button->Draw();
    
   // display logo
    TDiamond *diamond = new TDiamond(0.05,0.015,0.95,0.22);
    diamond->SetFillColor(50);
    diamond->SetTextAlign(22);
    diamond->SetTextColor(5);
    diamond->SetTextSize(0.11);
    diamond->Draw();
    diamond->AddText(".. ");
    diamond->AddText("ROOT");
    diamond->AddText("MUON");
    diamond->AddText("... ");
    diamond->AddText(" ");
}

//_____________________________________________________________________________
void AliMUONDisplay::CreateColors()
{
//    Create the colors palette used to display clusters

    Int_t k,i;
    Int_t color;
    Float_t r,g,b;
    
    for (k=1;k<=5;k++) {
	switch(k) {
	case 1:
	    for (i=1;i<=5;i++) {
		r=1.;
		g=i*0.2;  
		b=0.;
		color=i;
		color=260+23-color;
		new TColor(color,r,g,b);
	    } 
	    break;
	case 2:
	    for (i=1;i<=4;i++) {
		r=1.1-i*0.2;
		g=1.;  
		b=0.;
		color=i+5;
		color=260+23-color;
		new TColor(color,r,g,b);
	    } 
	    break;
	case 3:
	    for (i=1;i<=4;i++) {
		r=0.;
		g=1.;  
		b=i*0.2+0.2;
		color=i+9;
		color=260+23-color;
		new TColor(color,r,g,b);
	    } 
	    break;
	case 4:
	    for (i=1;i<=4;i++) {
		r=0.;
		g=1.1-i*0.2;  
		b=1.;
		color=i+13;
		color=260+23-color;
		new TColor(color,r,g,b);
	    } 
	    break;
	case 5:
	    for (i=1;i<=5;i++) {
		r=i*0.2;
		g=0.;  
		b=1.;
		color=i+17;
		color=260+23-color;
		new TColor(color,r,g,b);
	    } 
	    break;
	}
    }
}

//_____________________________________________________________________________
void AliMUONDisplay::DisplayColorScale()
{
// Display pulse height color scale
    Int_t i;
    Int_t color;
    Float_t xlow, ylow, xup, yup, hs;
    Float_t x1, y1, x2, y2;
    x1 = y1 = 0;
    x2 = y2 = 1.0;
    
    TText *text = new TText(0,0,"");
    text->SetTextFont(61);
    text->SetTextSize(0.2);
    text->SetTextAlign(22);
    
    AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
    AliMUONChamber *iChamber = &(pMUON->Chamber(fChamber-1));
    AliMUONResponse * response=iChamber->ResponseModel();
    Int_t adcmax=1024;
    if (response) adcmax= (Int_t) response->MaxAdc();
    

    TBox *box;
    char label[8];
//*-* draw colortable boxes
    hs = (y2-y1)/Float_t(22);
    xlow=x1+.05;
    xup=x2-0.5;
    for (i=0;i<22;i++) {
	ylow = y1 + hs*(Float_t(i));
	yup  = y1 + hs*(Float_t(i+1));
	color = 261+i;
	Double_t logscale=Double_t(i+1)*(TMath::Log(adcmax)/22);
	Int_t scale=(Int_t)TMath::Exp(logscale);
	sprintf(label,"%d",scale);
	box = new TBox(xlow, ylow, xup, yup);
	box->Draw();
	box->SetFillColor(color);
	text->DrawText(xlow+0.7, 0.5*(ylow+yup),label);
    }
}

//______________________________________________________________________________
Int_t AliMUONDisplay::DistancetoPrimitive(Int_t px, Int_t)
{
// Compute distance from point px,py to objects in event

    gPad->SetCursor(kCross);
    
    if (gPad == fTrigPad) return 9999;
    
    const Int_t kBig = 9999;
    Int_t dist   = kBig;
    Float_t xmin = gPad->GetX1();
    Float_t xmax = gPad->GetX2();
    Float_t dx   = 0.02*(xmax - xmin);
    Float_t x    = gPad->AbsPixeltoX(px);
    if (x < xmin+dx || x > xmax-dx) return dist;
    
    if (fZoomMode) return 0;
    else           return 7;
}

//_____________________________________________________________________________
void AliMUONDisplay::Draw(Option_t *)
{
//    Display current event

    fPad->cd();

    DrawView(fTheta, fPhi, fPsi);   
    // Display the event number and title
    fPad->cd();
    DrawTitle();
}

void AliMUONDisplay::DrawSegmentation()
{
// Draw graphical representation of segmenatation
// Attention: still experimental code
    Int_t icat=1;
    
    AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
    AliMUONChamber*   iChamber;
    AliSegmentation*  seg;
    iChamber = &(pMUON->Chamber(fChamber));
    seg=iChamber->SegmentationModel(icat);
    Float_t zpos=iChamber->Z();
    Float_t r=iChamber->ROuter();
    
    TMarker3DBox *marker;
    if (icat == 1) {
	for (Int_t j=0; j<seg->Npy(); j++) {
	    Float_t y0;
	    y0=j*seg->Dpy()-seg->Dpy()/2.;
	    for (seg->FirstPad(0.,y0,0,300,0.); 
		 seg->MorePads();
		 seg->NextPad())
	    {
		if (seg->ISector()==0) continue;
		Float_t x,y,z;
		seg->GetPadC(seg->Ix(), seg->Iy(), x, y, z);
		Float_t dpx=seg->Dpx(seg->ISector())/2;
		Float_t dpy=seg->Dpy(seg->ISector())/2;
		marker=new TMarker3DBox(x,y,zpos,dpx,dpy,0,0,0);
		marker->SetLineColor(seg->ISector()+1);
		marker->SetFillStyle(1001);
		marker->SetFillColor(0);
		marker->Draw();
	    }
	}
    } else {
	for (Int_t j=0; j<250; j++) {
	    Float_t x0=j*seg->Dpx();
	    Float_t y0=TMath::Sqrt(r*r-x0*x0);
	    
	    for (seg->FirstPad(x0,0,0,0,y0); 
		 seg->MorePads();
		 seg->NextPad())
	    {
		if (seg->ISector()==0) continue;
		
		Float_t x,y,z;
		seg->GetPadC(seg->Ix(), seg->Iy(), x, y, z);
		Float_t dpx=seg->Dpx(seg->ISector())/2;
		Float_t dpy=seg->Dpy(seg->ISector())/2;
		marker=new TMarker3DBox(x,y,zpos,dpx,dpy,0,0,0);
		marker->SetLineColor(seg->ISector()+1);
		marker->SetFillStyle(1001);
		marker->SetFillColor(0);
		marker->Draw();
	    }
	}
    }
}

//_____________________________________________________________________________
void AliMUONDisplay::DrawClusters()
{
//    Draw clusters for MUON chambers

    Int_t ndigits, digit;
    TObjArray *points;
    AliMUONPoints *pm;

      
    fClustersCuts = 0;
    points = Points();
    if (!points) return;
    ndigits = points->GetEntriesFast();
    for (digit=0;digit<ndigits;digit++){
	pm = (AliMUONPoints*)points->UncheckedAt(digit);
	if (!pm) continue;
	Float_t *pxyz;
	pxyz=pm->GetP();
	for (Int_t im=0;im<3;im++) {
	    TMarker3DBox *marker=pm->GetMarker(im);
	    if (marker)
		marker->Draw();
	}
	pm->Draw();
	fClustersCuts +=pm->GetN();
    }
}

//_____________________________________________________________________________
void AliMUONDisplay::DrawHits()
{
//    Draw hits for MUON chambers

    LoadHits(fChamber);

    Int_t ntracks, track;
    TObjArray *points;
    AliMUONPoints *pm;
    
    fHitsCuts = 0;
    points = Phits();
    if (!points) return;
    ntracks = points->GetEntriesFast();
    for (track=0;track<ntracks;track++) {
	pm = (AliMUONPoints*)points->UncheckedAt(track);
	if (!pm) continue;
	pm->Draw();
	fHitsCuts += pm->GetN();
    }
}


//_____________________________________________________________________________
void AliMUONDisplay::DrawCoG()
{
//    Draw hits for MUON chambers
    if (!fDrawCoG) return;
    if (fChamber > 10) return;
    LoadCoG(fChamber,fCathode);
    
    Int_t ncog, icog;
    TObjArray *points;
    AliMUONPoints *pm;

    points = Rpoints();
    if (!points) return;
    ncog = points->GetEntriesFast();
    for (icog=0;icog<ncog;icog++) {
	pm = (AliMUONPoints*)points->UncheckedAt(icog);
	if (!pm) continue;
	pm->Draw();
    }
}
void AliMUONDisplay::DrawCoG2()
{
//    Draw hits for MUON chambers

    if (!fDrawCoG) return;
    if (fChamber > 10) return;  

    if (fCathode==1) {
	LoadCoG2(fChamber,2);
    } else if (fCathode==2) {
	LoadCoG2(fChamber,1);
    }

    Int_t ncog, icog;
    TObjArray *points;
    AliMUONPoints *pm;
    
    points = R2points();
    if (!points) return;
    ncog = points->GetEntriesFast();
    for (icog=0;icog<ncog;icog++) {
	pm = (AliMUONPoints*)points->UncheckedAt(icog);
	if (!pm) continue;
	pm->Draw();
    }
}
//_____________________________________________________________________________

void AliMUONDisplay::DrawTitle(Option_t *option)
{
//    Draw the event title

    Float_t xmin = gPad->GetX1();
    Float_t xmax = gPad->GetX2();
    Float_t ymin = gPad->GetY1();
    Float_t ymax = gPad->GetY2();
    Float_t dx   = xmax-xmin;
    Float_t dy   = ymax-ymin;
    
    if (strlen(option) == 0) {
	TPaveText *title = new TPaveText(xmin +0.01*dx, ymax-0.09*dy, xmin +0.5*dx, ymax-0.01*dy);
//      title->SetTextSize(0.023932);
	title->SetTextSize(0.02);
	title->SetBit(kCanDelete);
	title->SetFillColor(42);
	title->Draw();
	char ptitle[100];
	sprintf(ptitle, "Alice event:%d Run:%d Chamber:%d Cathode:%d",
		gAlice->GetHeader()->GetEvent(),
		gAlice->GetHeader()->GetRun(),
		fChamber,
		fCathode);
	title->AddText(ptitle);
	Int_t nparticles = gAlice->Particles()->GetEntriesFast();
	sprintf(ptitle,"Nparticles = %d Nhits = %d Npads fired = %d",
		nparticles, fHitsCuts,fClustersCuts);
	title->AddText(ptitle);
    } else {
	TPaveLabel *label = new TPaveLabel(xmin +0.01*dx, ymax-0.07*dy, xmin +0.2*dx, ymax-0.01*dy,option);
	label->SetBit(kCanDelete);
	label->SetFillColor(42);
	label->Draw();
    }
}

//_____________________________________________________________________________
void AliMUONDisplay::DrawView(Float_t theta, Float_t phi, Float_t psi)
{
//    Draw a view of MUON clusters
    printf("\n Draw View");
    
    gPad->SetCursor(kWatch);
    // gPad->SetFillColor(39);
    gPad->SetFillColor(1);
    gPad->Clear();
    // gPad->SetFillColor(39);
    gPad->SetFillColor(1);
    

    Int_t iret=0;
    TView *view = new TView(1);
    
    Float_t range = fRrange*fRangeSlider->GetMaximum();
    view->SetRange(-range,-range,-range,range, range, range);
    // zoom back to full scale only if DrawView not called from NextCathode
    if (!fNextCathode) {
	fZoomX0[0] = -1;
	fZoomY0[0] = -1;
	fZoomX1[0] =  1;
	fZoomY1[0] =  1;
	fZooms = 0;
    }

// Display MUON Chamber Geometry
    char nodeName[7];
    sprintf(nodeName,"MUON%d",100+fChamber);
    printf("\n Node name %s %p", nodeName, gAlice->GetGeometry());
    
    TNode *node1=gAlice->GetGeometry()->GetNode(nodeName);
    if (node1) node1->Draw("same");  
//add clusters to the pad
    DrawClusters();
    printf("Node name %s", nodeName);   
    DrawHits();
    DrawCoG();
    DrawCoG2();
//     DrawSegmentation();
    // add itself to the list (must be last)
    AppendPad();
    view->SetView(phi, theta, psi, iret);
}


//______________________________________________________________________________
void AliMUONDisplay::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
//  Execute action corresponding to the mouse event

    static Float_t x0, y0, x1, y1;
    
    static Int_t pxold, pyold;
    static Int_t px0, py0;
    static Int_t linedrawn;
    Float_t temp;
    
    if (px == 0 && py == 0) { //when called by sliders
	if (event == kButton1Up) {
	    Draw();
	}
	return;
    }
    if (!fZoomMode && gPad->GetView()) {
	gPad->GetView()->ExecuteRotateView(event, px, py);
	return;
    }

    // something to zoom ?
    gPad->SetCursor(kCross);
    
    switch (event) {

    case kButton1Down:
	gVirtualX->SetLineColor(-1);
	gPad->TAttLine::Modify();  //Change line attributes only if necessary
	x0 = gPad->AbsPixeltoX(px);
	y0 = gPad->AbsPixeltoY(py);
	px0   = px; py0   = py;
	pxold = px; pyold = py;
	linedrawn = 0;
	return;
	
    case kButton1Motion:
	if (linedrawn) gVirtualX->DrawBox(px0, py0, pxold, pyold, TVirtualX::kHollow);
	pxold = px;
	pyold = py;
	linedrawn = 1;
	gVirtualX->DrawBox(px0, py0, pxold, pyold,  TVirtualX::kHollow);
	return;
	
    case kButton1Up:
	gPad->GetCanvas()->FeedbackMode(kFALSE);
	if (px == px0) return;
	if (py == py0) return;
	x1 = gPad->AbsPixeltoX(px);
	y1 = gPad->AbsPixeltoY(py);
	
	if (x1 < x0) {temp = x0; x0 = x1; x1 = temp;}
	if (y1 < y0) {temp = y0; y0 = y1; y1 = temp;}
	gPad->Range(x0,y0,x1,y1);
	if (fZooms < AliMUONConstants::MaxZoom()-1) {
	    fZooms++;
	    fZoomX0[fZooms] = x0;
	    fZoomY0[fZooms] = y0;
	    fZoomX1[fZooms] = x1;
	    fZoomY1[fZooms] = y1;
	}
	gPad->Modified(kTRUE);
	return;
    }
}
 
//___________________________________________
void AliMUONDisplay::LoadDigits(Int_t chamber, Int_t cathode)
{
// Read digits info and store x,y,z info in arrays fPoints
// Loop on all detectors

    if (chamber > 14) return;
    printf(" chamber %d \n",chamber);
    fChamber=chamber;
    fCathode=cathode;
    
    ResetPoints();
    
    AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
    AliMUONChamber*       iChamber;
    AliSegmentation*      segmentation;
    AliMUONResponse*      response;

    TClonesArray *muonDigits  = pMUON->DigitsAddress(chamber-1);
    if (muonDigits == 0) return;

    gAlice->ResetDigits();
    Int_t nent = 0;
 
   if (gAlice->TreeD()) {
	nent=(Int_t)gAlice->TreeD()->GetEntries();
	gAlice->TreeD()->GetEvent(nent-2+cathode-1);
    }
    
    Int_t ndigits = muonDigits->GetEntriesFast();
    if (ndigits == 0) return;
    if (fPoints == 0) fPoints = new TObjArray(ndigits);
    
    iChamber = &(pMUON->Chamber(chamber-1));

    segmentation = iChamber->SegmentationModel(cathode);
    response     = iChamber->ResponseModel();
    Float_t zpos=iChamber->Z();  
    AliMUONDigit  *mdig;
    AliMUONPoints *points = 0;
    TMarker3DBox *marker=0;
    //
    //loop over all digits and store their position
    
    Int_t npoints=1;
    Float_t adcmax=1024;
    if (response) adcmax= response->MaxAdc();

    for (Int_t digit=0;digit<ndigits;digit++) {
        mdig    = (AliMUONDigit*)muonDigits->UncheckedAt(digit);
        //
        // First get all needed parameters
        //
        Int_t charge=mdig->fSignal;
        Int_t index=Int_t(TMath::Log(charge)/(TMath::Log(adcmax)/22));
        Int_t color=261+index;
	Int_t colorTrigger=2;   
        if (color>282) color=282;

	if (chamber > 10) { // trigger chamber 
	    Int_t sumCharge=0;
	    for (Int_t icharge=0; icharge<10; icharge++) {
		sumCharge=sumCharge+mdig->fTcharges[icharge];
	    }
	    Int_t testCharge=sumCharge-(Int_t(sumCharge/10))*10;
	    if(sumCharge<=10||testCharge>0) {
		colorTrigger=color; 
	    } else {
		colorTrigger=5; 
	    }
	}

	// get the center of the pad - add on x and y half of pad size
	Float_t xpad, ypad, zpad;
	segmentation->GetPadC(mdig->fPadX, mdig->fPadY,xpad, ypad, zpad);
	
        Int_t isec=segmentation->Sector(mdig->fPadX, mdig->fPadY);
        Float_t dpx=segmentation->Dpx(isec)/2;
        Float_t dpy=segmentation->Dpy(isec)/2;
	Int_t nPara, offset;
        segmentation->GetNParallelAndOffset(mdig->fPadX,mdig->fPadY,
		&nPara,&offset);
	//
	// Then set the objects
	//
        points = new AliMUONPoints(npoints);
	fPoints->AddAt(points,digit);
	if (chamber > 10) {
	    points->SetMarkerColor(colorTrigger);
	} else {  
	    points->SetMarkerColor(color);
	}
        points->SetMarkerStyle(21);
        points->SetMarkerSize(0.5);
        points->SetParticle(-1);
        points->SetHitIndex(-1);
        points->SetTrackIndex(-1);
        points->SetDigitIndex(digit);
        points->SetPoint(0,xpad,ypad,zpos);	
	for (Int_t imark=0;imark<nPara; imark++)
	{
	    Int_t lineColor = (zpad-zpos > 0) ? 2:3;
	    segmentation->GetPadC(mdig->fPadX + imark*offset, mdig->fPadY,xpad, ypad, zpad);
	    marker=new TMarker3DBox(xpad,ypad,zpos,dpx,dpy,0,0,0);
	    marker->SetLineColor(lineColor);
	    marker->SetFillStyle(1001);
	    marker->SetFillColor(color);
	    marker->SetRefObject((TObject*)points);
	    points->Set3DMarker(imark, marker);
	}
    }
}
//___________________________________________
void AliMUONDisplay::LoadCoG(Int_t chamber, Int_t cathode)
{
// Read raw clusters info and store x,y,z info in arrays fRpoints
// Loop on all detectors

    if (chamber > 10) return;
    
    ResetRpoints();
    
    AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
    AliMUONChamber*  iChamber;
    
    TClonesArray *muonRawClusters  = pMUON->RawClustAddress(chamber-1);
    if (muonRawClusters == 0) return;

    pMUON->ResetRawClusters();

    Int_t nent = 0;
    if (gAlice->TreeR()) {
	nent=(Int_t)gAlice->TreeR()->GetEntries();
	gAlice->TreeR()->GetEvent(nent-2+cathode-1);
    }
    
    Int_t nrawcl = muonRawClusters->GetEntriesFast();
    if (nrawcl == 0) return;
    if (fRpoints == 0) fRpoints = new TObjArray(nrawcl);
    
    iChamber = &(pMUON->Chamber(chamber-1));
    Float_t zpos=iChamber->Z();  
    AliMUONRawCluster  *mRaw;
    AliMUONPoints *points = 0;
    //
    //loop over all raw clusters and store their position
    points = new AliMUONPoints(nrawcl);
    for (Int_t iraw=0;iraw<nrawcl;iraw++) {
  	mRaw   = (AliMUONRawCluster*)muonRawClusters->UncheckedAt(iraw);
	fRpoints->AddAt(points,iraw);
        points->SetMarkerColor(51);
        points->SetMarkerStyle(2);
        points->SetMarkerSize(1.);
        points->SetParticle(-1);
        points->SetHitIndex(-1);
        points->SetTrackIndex(-1);
        points->SetDigitIndex(-1);
        points->SetPoint(iraw,mRaw->fX[0],mRaw->fY[0],zpos);
    }
}
//___________________________________________
void AliMUONDisplay::LoadCoG2(Int_t chamber, Int_t cathode)
{
// Read raw clusters info and store x,y,z info in arrays fRpoints
// Loop on all detectors

    if (chamber > 10) return;

    ResetR2points();

    AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
    AliMUONChamber*  iChamber;
    
    TClonesArray *muonRawClusters  = pMUON->RawClustAddress(chamber-1);
    if (muonRawClusters == 0) return;
    
    pMUON->ResetRawClusters();
    
    Int_t nent = 0;
    if (gAlice->TreeR()) {
	nent=(Int_t)gAlice->TreeR()->GetEntries();
	gAlice->TreeR()->GetEvent(nent-2+cathode-1);
    }
    
    Int_t nrawcl = muonRawClusters->GetEntriesFast();
    if (nrawcl == 0) return;
    if (fR2points == 0) fR2points = new TObjArray(nrawcl);
    
    iChamber = &(pMUON->Chamber(chamber-1));
    Float_t zpos=iChamber->Z();  
    AliMUONRawCluster  *mRaw;
    AliMUONPoints *points = 0;
    //
    //loop over all raw clusters and store their position
    points = new AliMUONPoints(nrawcl);
    for (Int_t iraw=0;iraw<nrawcl;iraw++) {
  	mRaw   = (AliMUONRawCluster*)muonRawClusters->UncheckedAt(iraw);
	fR2points->AddAt(points,iraw);
        points->SetMarkerColor(51);
        points->SetMarkerStyle(4);
        points->SetMarkerSize(1.3);
        points->SetParticle(-1);
        points->SetHitIndex(-1);
        points->SetTrackIndex(-1);
        points->SetDigitIndex(-1);
        points->SetPoint(iraw,mRaw->fX[0],mRaw->fY[0],zpos);
   }
}
//___________________________________________
void AliMUONDisplay::LoadHits(Int_t chamber)
{
// Read hits info and store x,y,z info in arrays fPhits
// Loop on all detectors

    if (chamber > 14) return;
    Int_t track;

    fChamber=chamber;
 
    ResetPhits();
    
    AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
    AliMUONChamber*  iChamber;

    iChamber = &(pMUON->Chamber(chamber-1));
    Float_t zpos=iChamber->Z();

    Int_t ntracks = (Int_t)gAlice->TreeH()->GetEntries();
    Int_t nthits=0;
    for (track=0; track<ntracks;track++) {
	gAlice->ResetHits();
	gAlice->TreeH()->GetEvent(track);
	TClonesArray *muonHits  = pMUON->Hits();
	if (muonHits == 0) return;
	nthits += muonHits->GetEntriesFast();
    } 
    if (fPhits == 0) fPhits = new TObjArray(nthits);
    Int_t nhold=0;
    for (track=0; track<ntracks;track++) {
	gAlice->ResetHits();
	gAlice->TreeH()->GetEvent(track);
	TClonesArray *muonHits  = pMUON->Hits();
	if (muonHits == 0) return;
	Int_t nhits = muonHits->GetEntriesFast();
	if (nhits == 0) continue;
	AliMUONHit *mHit;
	AliMUONPoints *points = 0;
	Int_t npoints=1;
	for (Int_t hit=0;hit<nhits;hit++) {
            mHit = (AliMUONHit*)muonHits->UncheckedAt(hit);
            Int_t nch  = mHit->fChamber;              // chamber number
            if (nch != chamber) continue;
	    //
	    // Retrieve info and set the objects
	    //
	    points = new AliMUONPoints(npoints);
	    fPhits->AddAt(points,nhold+hit);
            points->SetMarkerColor(kRed);
            points->SetMarkerStyle(5);
            points->SetMarkerSize(1.);
            points->SetParticle(mHit->Track());
            points->SetHitIndex(hit);
            points->SetTrackIndex(track);
            points->SetDigitIndex(-1);
	    points->SetPoint(0,mHit->X(),mHit->Y(),zpos);
	}
	nhold+=nhits;
    }
}

//_____________________________________________________________________________
void AliMUONDisplay::Paint(Option_t *)
{
//    Paint miscellaneous items
}

//_____________________________________________________________________________
void AliMUONDisplay::SetPickMode()
{
// Set parameters for pick mode.
// 
    fZoomMode = 0;

    fArcButton->SetY1(fPickButton->GetYlowNDC()+0.5*fPickButton->GetHNDC());
    fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliMUONDisplay::SetZoomMode()
{
//  Set parameters for zoom mode
    fZoomMode = 1;
    
    fArcButton->SetY1(fZoomButton->GetYlowNDC()+0.5*fZoomButton->GetHNDC());
    fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliMUONDisplay::NextChamber(Int_t delta)
{
  // to go from chamber to next chamber if delta = 1
  // or previous chamber otherwise
    if (delta == 1) {
	if (fChamber < AliMUONConstants::NCh()) fChamber++;
    } else {
	if (fChamber > 1) fChamber--;
    }
    if (!fPad) return;
    fPad->Clear();
    LoadDigits(fChamber, fCathode);
    Draw();
}

//_____________________________________________________________________________
void AliMUONDisplay::NextCathode()
{
    // to switch to other cathode plane
    if (!fPad) return;
    fPad->Clear();
    if (fCathode == 1) {
	LoadDigits(fChamber, 2);	
    } else {
	LoadDigits(fChamber, 1);
    }
    fNextCathode = kTRUE; // to keep the same zoom
    Draw();
    fNextCathode = kFALSE;
    TPad *pad = (TPad*)gPad->GetPadSave();
    pad->Range(fZoomX0[fZooms], fZoomY0[fZooms],
	       fZoomX1[fZooms], fZoomY1[fZooms]);
    pad->Modified();
    fPad->cd();
    DrawTitle();
}

//_____________________________________________________________________________
void AliMUONDisplay::Trigger()
{
  // returns Trigger Decision for current event
  AliMUONTriggerDecision* decision= new AliMUONTriggerDecision(1);
  decision->Trigger(); 
}


//_____________________________________________________________________________
void AliMUONDisplay::SetChamberAndCathode(Int_t chamber, Int_t cathode)
{
// Set chamber and cathode number
   fChamber = chamber;
   fCathode = cathode;

   if (!fPad) return;
   fPad->Clear();
   LoadDigits(chamber,cathode);
   Draw();
}

void AliMUONDisplay::SetEvent(Int_t newevent)
{
// Chose event 
    gAlice->GetEvent(newevent);
    fEvent=newevent;
    if (!gAlice->TreeD()) return; 
    if (!fPad) return;
    fPad->Clear();
    LoadDigits(fChamber,fCathode);
    Draw();
}

//_____________________________________________________________________________
void AliMUONDisplay::SetRange(Float_t rrange, Float_t zrange)
{
// Set view range along R and Z
    fRrange = rrange;
    fZrange = zrange;

    if (!fPad) return;
    fPad->Clear();
    Draw();
}
   
//_____________________________________________________________________________
void AliMUONDisplay::SetView(Float_t theta, Float_t phi, Float_t psi)
{
//  change viewing angles for current event

    fPad->cd();
    fPhi   = phi;
    fTheta = theta;
    fPsi   = psi;
    Int_t iret = 0;
    
    TView *view = gPad->GetView();
    if (view) view->SetView(fPhi, fTheta, fPsi, iret);
    else      Draw();
    gPad->Modified();
}

//_____________________________________________________________________________
void AliMUONDisplay::ShowNextEvent(Int_t delta)
{
//  Display (current event_number + delta)
//    delta =  1  shown next event
//    delta = -1 show previous event
    if (delta) {
	gAlice->Clear();
	Int_t currentEvent = gAlice->GetHeader()->GetEvent();
	Int_t newEvent     = currentEvent + delta;
	gAlice->GetEvent(newEvent);
	fEvent=newEvent;
	if (!gAlice->TreeD()) return; 
    }
    
    LoadDigits(fChamber, fCathode);
    fPad->cd(); 
    Draw();
}

//______________________________________________________________________________
void AliMUONDisplay::UnZoom()
{
// Unzoom 
    if (fZooms <= 0) return;
    fZooms--;
    TPad *pad = (TPad*)gPad->GetPadSave();
    pad->Range(fZoomX0[fZooms],fZoomY0[fZooms], fZoomX1[fZooms],fZoomY1[fZooms]);
    pad->Modified();
}

//_____________________________________________________________________________
void AliMUONDisplay::ResetPoints()
{
    //
    // Reset array of points
    //
    if (fPoints) {
	fPoints->Delete();
	delete fPoints;
	fPoints = 0;
    }
}
//_____________________________________________________________________________
void AliMUONDisplay::ResetPhits()
{
    //
    // Reset array of points
    //
    if (fPhits) {
	fPhits->Delete();
	delete fPhits;
	fPhits = 0;
    }
}
//_____________________________________________________________________________
void AliMUONDisplay::ResetRpoints()
{
  //
  // Reset array of points
  //
    if (fRpoints) {
	fRpoints->Delete();
	delete fRpoints;
	fRpoints = 0;
    }
}
//_____________________________________________________________________________
void AliMUONDisplay::ResetR2points()
{
  //
  // Reset array of points
  //
    if (fR2points) {
	fR2points->Delete();
	delete fR2points;
	fR2points = 0;
    }
}
//_____________________________________________________________________________
void AliMUONDisplay::ResetCpoints()
{
    //
    // Reset array of points
    //
  if (fCpoints) {
      fCpoints->Delete();
      delete fCpoints;
      fCpoints = 0;
  }
}


AliMUONDisplay & AliMUONDisplay::operator = (const AliMUONDisplay &)
{
// Dummy assignment operator
    return *this;
}









