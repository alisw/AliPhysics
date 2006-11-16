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

/* $Id$ */

/// \class AliMUONDisplay
/// Create an event display object.
/// A canvas named "edisplay" is created with a vertical size in pixels  \n
///
///    A QUICK Overview of the Event Display functions                   \n
///    ===============================================                   \n
///
/// The event display can ve invoked by executing the macro "display.C"  \n
/// A canvas like in the picture below will appear.
///
///  On the left side of the canvas, the following buttons appear:
///  - *Next*       to move to the next event
///  - *Previous*   to move to the previous event
///  - *Pick*       Select this option to be able to point on a track with the
///                 mouse. Once on the track, use the right button to select
///                 an action. For example, select SetMarkerAttributes to
///                 change the marker type/color/size for the track.
///  - *Zoom*       Select this option (default) if you want to zoom.
///                 To zoom, simply select the selected area with the left button.
///  - *UnZoom*     To revert to the previous picture size.
///
///  - slider R     On the left side, the vertical slider can be used to
///                 set the default picture size.
///
///  When you are in Zoom mode, you can click on the black part of the canvas
///  to select special options with the right mouse button.
///
/// When you are in pick mode, you can "Inspect" the object pointed by the mouse.
/// When you are on a track, select the menu item "InspectParticle"
/// to display the current particle attributes.
///
/// You can activate the Root browser by selecting the Inspect menu
/// in the canvas tool bar menu. Then select "Start Browser"
/// This will open a new canvas with the browser. At this point, you may want
/// to display some histograms (from the Trees). Go to the "File" menu
/// of the browser and click on "New canvas".
/// In the browser, click on item "ROOT files" in the left pane.
/// Click on galice.root.
/// Click on TH
/// Click on TPC for example
/// Click on any variable (eg TPC.fX) to histogram the variable.
///
/// If you are lost, you can click on HELP in any Root canvas or browser.

//Begin_Html
/*
<img src="gif/AliMUONDisplay.gif">
*/
//End_Html

#include "AliMUONDisplay.h"
#include "AliMUON.h"
#include "AliMUONPoints.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONHit.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONSegmentation.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONConstants.h"
#include "AliMUONTriggerSegmentation.h"

#include "AliMpDEIterator.h"
#include "AliMpSegmentation.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSector.h"
#include "AliMpTriggerSegmentation.h"
#include "AliMpTrigger.h"
#include "AliMpStationType.h"
#include "AliMpDEManager.h"

#include "AliMC.h"
#include "AliLog.h"
#include "AliRun.h"
#include "AliHeader.h"

#include <TButton.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TView.h>
#include <TText.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TDiamond.h>
#include <TNode.h>
#include <TArc.h>
#include <TSlider.h>
#include <TVirtualX.h>
#include <TMath.h>
#include <TGeometry.h>
#include <TMarker3DBox.h>
#include <TParticle.h>
#include <TPolyLine3D.h>
#include <TBox.h>

/// \cond CLASSIMP
ClassImp(AliMUONDisplay)
/// \endcond

//_____________________________________________________________________________
AliMUONDisplay::AliMUONDisplay()
  : AliDisplay(),
    fEvent(0),
    fChamber(0),
    fCathode(0),
    fDrawClusters(kTRUE),
    fDrawCoG(kTRUE),
    fDrawTracks(kFALSE),
    fClustersCuts(0),
    fColPad(0),
    fPoints(0),
    fPhits(0),
    fRpoints(0),
    fNextCathode(0),
    fLoader(0), 
    fMUONData(0)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMUONDisplay::AliMUONDisplay(Int_t size, AliLoader * loader)
  : AliDisplay(),
    fEvent(0),
    fChamber(1),
    fCathode(1),
    fDrawClusters(kTRUE),
    fDrawCoG(kTRUE),
    fDrawTracks(kFALSE),
    fClustersCuts(0),
    fColPad(0),
    fPoints(0),
    fPhits(0),
    fRpoints(0),
    fNextCathode(kFALSE),
    fLoader(loader), 
    fMUONData(0)
  
{
/// Standard constructor to create an event display object.

    fPad = 0;
    
    gAlice->SetDisplay(this);
   
   // Initialize display default parameters
    SetRange(200,2000);

   // Set front view by default
    fTheta =   0;
    fPhi   = -90;
    fPsi   =   0;
    fZoomMode      = 1;
    fZooms         = 0;

    // Create colors
    CreateColors();
    // Create display canvas
    Int_t ysize = size;
    if (ysize < 100) ysize = 750;
    Int_t xsize = Int_t(size*830./ysize);
    fCanvas = new TCanvas("Canvas", "MUON Display",14,47,xsize,ysize);
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
    fTrigPad->SetEditable(kFALSE);
    fButtons->SetEditable(kFALSE);
    fCanvas->Update();
 
    // initialize container
    if(fLoader) 
      fMUONData = new AliMUONData(fLoader,"MUON","MUON");
    else
      fMUONData =0x0;
}

//_____________________________________________________________________________
AliMUONDisplay::~AliMUONDisplay()
{
/// Destructor

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
}

//_____________________________________________________________________________
void AliMUONDisplay::Clear(Option_t *)
{
/// Delete graphics temporary objects
}

//_____________________________________________________________________________
void AliMUONDisplay::DisplayButtons()
{
/// Create the user interface buttons


    fButtons = new TPad("buttons", "newpad",0,0.45,0.15,1);
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
    
    y -= dbutton + dy;
    char but8[] = "((AliMUONDisplay*)(gAlice->Display()))->DrawReco()";
    button = new TButton("Tracking", but8, x0, y - dbutton, x1, y);
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
void AliMUONDisplay::CreateColors() const
{
/// Create the colors palette used to display clusters

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
/// Display pulse height color scale

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
    

    Int_t adcmax=4096; // default 12 bits ADC

    

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
/// Compute distance from point px,py to objects in event

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
/// Display current event

    if (!fDrawTracks) 
      DrawChamber();   
    else 
      DrawReco();
    
}
//_____________________________________________________________________________
void AliMUONDisplay::DrawChamber()
{
/// Display current event

  fDrawTracks = kFALSE;
  fPad->cd();
  DrawView(fTheta, fPhi, fPsi);   
  // Display the event number and title
  fPad->cd();
  DrawTitle();
 
}
//_____________________________________________________________________________
void AliMUONDisplay::DrawReco(Option_t *)
{
/// Display current event

  fDrawTracks = kTRUE;
  // print kinematics of generated particles
    PrintKinematics();
    // Draw global view of muon system
    fPad->cd();
    DrawGlobalView(135, -50, -140); 
  
    // Display the event number and title
    fPad->cd();
    DrawTitle();
}

//_____________________________________________________________________________
void AliMUONDisplay::PrintKinematics()
{
/// Print kinematic tree

  AliRunLoader * runLoader;
  TParticle *particle = new TParticle();
  Int_t nPart;
  Float_t vertex[3], momentum[3];

  if (fLoader)
    runLoader = fLoader->GetRunLoader();
  else
    runLoader = 0x0;
  
  printf("******  Event # %d ******\n",runLoader->GetEventNumber());
  runLoader->TreeK()->GetBranch("Particles")->SetAddress(&particle);
  nPart = (Int_t)runLoader->TreeK()->GetEntries();
  for(Int_t iPart = 0; iPart < nPart; iPart++) {
    runLoader->TreeK()->GetEvent(iPart);
    vertex[0] = particle->Vx();
    vertex[1] = particle->Vy();
    vertex[2] = particle->Vz();
    momentum[0] = particle->Px();
    momentum[1] = particle->Py();
    momentum[2] = particle->Pz();
    
    printf("===================================================\n");
    printf(" Generated particle # %d \n",iPart);
    printf(" name: %s \n",particle->GetName());
    printf(" vertex x,y,z (cm): %f %f %f \n",vertex[0],vertex[1],vertex[2]); 
    printf(" momentum Px,Py,Pz (GeV/c): %f %f %f \n",momentum[0],momentum[1],momentum[2]);
  }
  delete particle;    
}

//_____________________________________________________________________________
void AliMUONDisplay::DrawSegmentation()
{
/// \todo to be re-written for new seg
/// Draw graphical representation of segmenatation
/// Attention: still experimental code
//     Int_t icat=1;
    
//     AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
//     AliMUONChamber*   iChamber;

//     AliSegmentation*  seg;
//     iChamber = &(pMUON->Chamber(fChamber));
//     seg=iChamber->SegmentationModel(icat);

//     Float_t zpos=iChamber->Z();
//     Float_t r=iChamber->ROuter();
    
//     TMarker3DBox *marker;
//     if (icat == 1) {
// 	for (Int_t j=0; j<seg->Npy(); j++) {
// 	    Float_t y0;
// 	    y0=j*seg->Dpy()-seg->Dpy()/2.;
// 	    for (seg->FirstPad(0.,y0,0,300,0.); 
// 		 seg->MorePads();
// 		 seg->NextPad())
// 	    {
// 		if (seg->ISector()==0) continue;
// 		Float_t x,y,z;
// 		seg->GetPadC(seg->Ix(), seg->Iy(), x, y, z);
// 		Float_t dpx=seg->Dpx(seg->ISector())/2;
// 		Float_t dpy=seg->Dpy(seg->ISector())/2;
// 		marker=new TMarker3DBox(x,y,zpos,dpx,dpy,0,0,0);
// 		marker->SetLineColor(seg->ISector()+1);
// 		marker->SetFillStyle(1001);
// 		marker->SetFillColor(0);
// 		marker->Draw();
// 	    }
// 	}
//     } else {
// 	for (Int_t j=0; j<250; j++) {
// 	    Float_t x0=j*seg->Dpx();
// 	    Float_t y0=TMath::Sqrt(r*r-x0*x0);
	    
// 	    for (seg->FirstPad(x0,0,0,0,y0); 
// 		 seg->MorePads();
// 		 seg->NextPad())
// 	    {
// 		if (seg->ISector()==0) continue;
		
// 		Float_t x,y,z;
// 		seg->GetPadC(seg->Ix(), seg->Iy(), x, y, z);
// 		Float_t dpx=seg->Dpx(seg->ISector())/2;
// 		Float_t dpy=seg->Dpy(seg->ISector())/2;
// 		marker=new TMarker3DBox(x,y,zpos,dpx,dpy,0,0,0);
// 		marker->SetLineColor(seg->ISector()+1);
// 		marker->SetFillStyle(1001);
// 		marker->SetFillColor(0);
// 		marker->Draw();
// 	    }
// 	}
//     }
 }

//_____________________________________________________________________________
void AliMUONDisplay::DrawClusters()
{
/// Draw clusters for MUON chambers

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
	fClustersCuts += pm->GetN();
    }
}

//_____________________________________________________________________________
void AliMUONDisplay::DrawHits()
{
/// Draw hits for MUON chambers

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
/// Draw hits for MUON chambers

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
//_____________________________________________________________________________
void AliMUONDisplay::DrawTracks()
{
/// Draw tracks

    if (!fDrawTracks) return;
    LoadTracks();
    
    Int_t nTrack, iTrack;
    TObjArray *points;
    TPolyLine3D *pm;

    points = Rpoints();
    if (!points) return;
    nTrack = points->GetEntriesFast();
    for ( iTrack = 0; iTrack < nTrack; iTrack++) {
	pm = (TPolyLine3D*)points->UncheckedAt(iTrack);
	if (!pm) continue;
	pm->Draw();
    }
}
//_____________________________________________________________________________

void AliMUONDisplay::DrawTitle(Option_t *option)
{
/// Draw the event title

    Float_t xmin = gPad->GetX1();
    Float_t xmax = gPad->GetX2();
    Float_t ymin = gPad->GetY1();
    Float_t ymax = gPad->GetY2();
    Float_t dx   = xmax-xmin;
    Float_t dy   = ymax-ymin;
    
    AliRunLoader * runLoader;
    if (fLoader)
      runLoader = fLoader->GetRunLoader();
    else
      runLoader = 0x0;


    if (strlen(option) == 0) {
	TPaveText *title = new TPaveText(xmin +0.01*dx, ymax-0.09*dy, xmin +0.5*dx, ymax-0.01*dy);
//      title->SetTextSize(0.023932);
	title->SetTextSize(0.02);
	title->SetBit(kCanDelete);
	title->SetFillColor(42);
	title->Draw();
	char ptitle[100];
	sprintf(ptitle, "Alice event:%d Run:%d Chamber:%d Cathode:%d",
		runLoader->GetEventNumber(),
		gAlice->GetHeader()->GetRun(),
		fChamber,
		fCathode);
	title->AddText(ptitle);
	Int_t nparticles = gAlice->GetMCApp()->Particles()->GetEntriesFast();
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
/// Draw a view of MUON clusters

  AliInfo(" Draw View");
  
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
  
  Float_t xg1, xg2, yg1, yg2, zg1, zg2;
  
  // Recovering the chamber 
  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
  
  const AliMUONGeometryTransformer* kGeomTransformer 
    = pMUON->GetGeometryTransformer();
  
  AliMUONSegmentation* segmentation = pMUON->GetSegmentation();
  
  // Display MUON Chamber Geometry
  char nodeName[7];
  sprintf(nodeName,"MUON%d",100+fChamber);
  printf(">>>> chamber is %d\n",fChamber);
  
  if(fChamber < 5) {
    AliMpDEIterator it;
    for ( it.First(fChamber-1); ! it.IsDone(); it.Next() ) {
      
      Int_t detElemId = it.CurrentDE();
      AliMpSectorSegmentation * seg =   
        (AliMpSectorSegmentation *) AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, 0);
      const AliMpSector * sector = seg->GetSector();
      
      // get sector measurements
      TVector2 position  = sector->Position(); 
      TVector2 dimension = sector->Dimensions(); // half length
      
      Float_t xlocal1 =  position.Px(); // FIXME: not really needed as it's 0 ?
      Float_t ylocal1 =  position.Py(); // FIXME: not really needed as it's 0 ?
      Float_t xlocal2 =  dimension.Px() * 2.;
      Float_t ylocal2 =  dimension.Px() * 2.;
      
      kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
      kGeomTransformer->Local2Global(detElemId, xlocal2, ylocal2, 0, xg2, yg2, zg2);
      
      // drawing 
      TPolyLine3D* poly = new  TPolyLine3D();
      Int_t nPoint = 0;
      
      poly->SetPoint(nPoint++, xg1, yg1, 0.);
      for (Float_t d = 0; d < TMath::Pi()/2.; d+= 0.01) {
        Float_t x = xg1 + xg2 * TMath::Cos(d);
        Float_t y = yg1 + yg2 * TMath::Sin(d);
        poly->SetPoint(nPoint++, x, y, 0.);
      }
      poly->SetPoint(nPoint++, xg1, yg1, 0.);
      
      poly->SetLineColor(2);
      poly->Draw("s");
    }
    
  }
  
  if (fChamber>4) 
  {
    AliMpDEIterator it;
    for ( it.First(fChamber-1); ! it.IsDone(); it.Next() ) 
    {
      Int_t detElemId = it.CurrentDE();
      AliMpStationType stationType = AliMpDEManager::GetStationType(detElemId);

      if (  segmentation->HasDE(detElemId) ) 
      {
        const AliMpVSegmentation* seg 
	  = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, 0);
        if (seg) 
        {  
          Float_t deltax = seg->Dimensions().X();
          Float_t deltay = seg->Dimensions().Y();
          Float_t xlocal1 =  -deltax;
          Float_t ylocal1 =  -deltay;
          Float_t xlocal2 =  +deltax;
          Float_t ylocal2 =  +deltay;
          kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
          kGeomTransformer->Local2Global(detElemId, xlocal2, ylocal2, 0, xg2, yg2, zg2);
          
          // drawing slat active volumes
          Float_t xCenter = (xg1 + xg2)/2.;
          Float_t yCenter = (yg1 + yg2)/2.;
          
          TMarker3DBox* box = new TMarker3DBox(xCenter,yCenter,0,xlocal1,ylocal2,0,0,0);
          
          box->SetFillStyle(0);
          box->SetLineColor( stationType == kStationTrigger ? 4 : 2);
          box->Draw("s");
          
          if ( stationType == kStation345 )
          {
            // drawing inner circle + disc
            TPolyLine3D* poly  = new  TPolyLine3D();
            TPolyLine3D* poly1 = new  TPolyLine3D();
            
            Int_t nPoint = 0;
            Int_t nPoint1 = 0;
            for (Float_t d = 0; d < 6.24; d+= 0.005) 
            {
              Double_t x = AliMUONConstants::Dmin((fChamber-1)/2) * TMath::Cos(d)/2.;
              Double_t y = AliMUONConstants::Dmin((fChamber-1)/2) * TMath::Sin(d)/2.;
              if (nPoint % 2 == 0) poly->SetPoint(nPoint++, 0., 0., 0.);
              poly->SetPoint(nPoint++, x, y, 0.);
              poly1->SetPoint(nPoint1++, x, y, 0.);            
            }
            poly->SetLineColor(1);
            poly->Draw("s");
            poly1->SetLineColor(2);
            poly1->Draw("s");
          }
        }  
      }
    }
  }
  
  //add clusters to the pad
  DrawClusters();
  DrawHits();
  DrawCoG();
  //     DrawSegmentation();
  // add itself to the list (must be last)
  AppendPad();
  view->SetView(phi, theta, psi, iret);
}

//_____________________________________________________________________________
void AliMUONDisplay::DrawGlobalView(Float_t theta, Float_t phi, Float_t psi)
{
/// Draw a view of muons chambers with tracks
    
    gPad->SetCursor(kWatch);
    // gPad->SetFillColor(39);
    gPad->SetFillColor(1);
    gPad->Clear();
    // gPad->SetFillColor(39);
    gPad->SetFillColor(1);
    

    Int_t iret=0;
    TView *view = new TView(1);
    
    Float_t range = fRrange*fRangeSlider->GetMaximum()*3.;
    view->SetRange(-range,-range,-range,range,range,range);

// Display all MUON Chambers segmentation
    char nodeName[7];
    TNode *node1;
    sprintf(nodeName,"alice");
    
    node1=gAlice->GetGeometry()->GetNode(nodeName);
    if (node1) node1->Draw("same"); 
     
      
// Draw clusters for all chambers
    Int_t chamberSave = fChamber;
    for (fChamber = 1; fChamber <= 10; fChamber++){
      DrawCoG();
    }
    fChamber = chamberSave;
// Draw reconstructed tracks
    DrawTracks();

    AppendPad();

    Float_t zoom = 2.;
    Float_t shift = 0.9;
    Float_t x0 = (-1+shift)/zoom;
    Float_t y0 = (-1+shift)/zoom;
    Float_t x1 = (1+shift)/zoom;
    Float_t y1 = (1+shift)/zoom;
    gPad->Range(x0,y0,x1,y1);
    view->SetView(phi, theta, psi, iret);

}

//______________________________________________________________________________
void AliMUONDisplay::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
/// Execute action corresponding to the mouse event

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
/// Read digits info and store x,y,z info in arrays fPoints.
/// Loop on all detectors

    if (chamber > 14) return;
    fChamber = chamber;
    fCathode = cathode;
    
    ResetPoints();
    
    AliMUON *pMUON  =     (AliMUON*)gAlice->GetModule("MUON");

    GetMUONData()->SetTreeAddress("D");

    TClonesArray *muonDigits  = GetMUONData()->Digits(chamber-1);
    if (muonDigits == 0) return;

    gAlice->ResetDigits();
    Int_t nent = 0;
 
    if (GetLoader()->TreeD()) {
	nent = (Int_t) GetLoader()->TreeD()->GetEntries();
	//     gAlice->TreeD()->GetEvent(nent-2+cathode-1);
	GetMUONData()->GetDigits();
    }
    
    Int_t ndigits = muonDigits->GetEntriesFast();    
    if (ndigits == 0) return;
    if (fPoints == 0) fPoints = new TObjArray(ndigits);
    
    Float_t zpos = AliMUONConstants::DefaultChamberZ(chamber-1);

    AliMUONDigit  *mdig;
    AliMUONPoints *points  = 0;
    TMarker3DBox  *marker  = 0;

    Int_t npoints  = 1;
    Float_t adcmax = 1024; // default
    if (chamber<11) adcmax = 4096;

// check if trigger is using new or old segmentation
    Bool_t old = true;
    AliMUONSegmentation* segmentation = pMUON->GetSegmentation();
    const AliMUONVGeometryDESegmentation* kdeSegmentation 
      = segmentation->GetDESegmentation(1100, cathode-1);
    if ( dynamic_cast<const AliMUONTriggerSegmentation*>(kdeSegmentation) ) old = false;

    if ( old  && chamber > 10) {
	if (chamber > 10) printf(">>> old segmentation for trigger \n");
	else printf(">>> old segmentation for tracking \n");

	for (Int_t digit = 0; digit < ndigits; digit++) {
	    mdig    = (AliMUONDigit*)muonDigits->UncheckedAt(digit);
	    if (mdig->Cathode() != cathode-1) continue;
	    
	    //
	    // First get all needed parameters
	    //
	    Float_t charge = mdig->Signal();
	    Int_t index  = Int_t(TMath::Log(charge)/(TMath::Log(adcmax)/22));
	    Int_t color  = 261+index;
	    Int_t colorTrigger = 2;   
	    if (color > 282) color = 282;
	    
	    if (chamber > 10) { // trigger chamber 
		
		Float_t sumCharge = 0;
		for (Int_t icharge = 0; icharge < 10; icharge++) {
		    sumCharge = sumCharge+mdig->TrackCharge(icharge);
		}
		Float_t testCharge = sumCharge-(Int_t(sumCharge/10))*10;
		if(sumCharge <= 10 || testCharge > 0) {
		    colorTrigger = color; 
		} else {
		    colorTrigger = 5; 
		}
	    }
	    
	    // get the center of the pad - add on x and y half of pad size
	    Float_t xpad, ypad, zpad;
	    Int_t isec;
	    Float_t dpx, dpy;
	    
	    Int_t detElemId = mdig->DetElemId();
	    AliMUONGeometrySegmentation* segmentation2
	      = pMUON->GetSegmentation()->GetModuleSegmentationByDEId(detElemId, cathode-1);
	    segmentation2->GetPadC(detElemId, mdig->PadX(), mdig->PadY(), xpad, ypad, zpad);
	    isec = segmentation2->Sector(detElemId, mdig->PadX(), mdig->PadY());
	    dpx = segmentation2->Dpx(detElemId, isec)/2;
	    dpy = segmentation2->Dpy(detElemId, isec)/2;
	    
	    // Then set the objects
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
	    
	    Int_t lineColor = (zpad-zpos > 0) ? 2:3;
	    marker=new TMarker3DBox(xpad,ypad,zpos,dpx,dpy,0,0,0);
	    
	    
	    marker->SetLineColor(lineColor);
	    marker->SetFillStyle(1001);
	    marker->SetFillColor(color);
	    marker->SetRefObject((TObject*)points);
	    points->Set3DMarker(0, marker);
	} // end loop on digits	    
	
    } else {	
	if (chamber > 10) printf(">>> new segmentation for trigger \n");
	else printf(">>> new segmentation for tracking \n");
	
	const AliMUONGeometryTransformer* kGeomTransformer 
	    = pMUON->GetGeometryTransformer();
	
	//loop over all digits and store their position
	for (Int_t digit = 0; digit < ndigits; digit++) {
	    mdig    = (AliMUONDigit*)muonDigits->UncheckedAt(digit);
	    if (mdig->Cathode() != cathode-1) continue;
	    
	    // get all needed parameters
	    Int_t ix=mdig->PadX();
	    Int_t iy=mdig->PadY();
	    Int_t detElemId=mdig->DetElemId();      
	    Float_t charge = mdig->Signal();
	    Int_t index  = Int_t(TMath::Log(charge)/(TMath::Log(adcmax)/22));
	    Int_t color  = 261+index;
	    Int_t colorTrigger = 2;   
	    if (color > 282) color = 282;
	    
	    const AliMpVSegmentation* seg = 
		AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,cathode-1);
	    
	    AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy),kTRUE);
	    
	    if (chamber > 10) { // trigger chamber
		Float_t sumCharge = 0;
    Int_t n = mdig->Ntracks();
		for (Int_t icharge = 0; icharge < n; icharge++) {
		    sumCharge = sumCharge+mdig->TrackCharge(icharge);
		}
		Float_t testCharge = sumCharge-(Int_t(sumCharge/n))*n;
		if(sumCharge <= n || testCharge > 0) {
		    colorTrigger = color; 
		} else {
		    colorTrigger = 5; 
		}
	    }
	    
	    // get the pad position and dimensions
	    Float_t xlocal1 = pad.Position().X();
	    Float_t ylocal1 = pad.Position().Y();
	    Float_t xlocal2 = pad.Dimensions().X();
	    Float_t ylocal2 = pad.Dimensions().Y();
	    
	    Float_t xg1, xg2, yg1, yg2, zg1;
	    
	    kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
	    // (no transformation for pad dimensions)
	    xg2 = xlocal2;
	    yg2 = ylocal2;
	    
	    // Then set the objects
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
	    points->SetPoint(0,xg1,yg1,zpos);
	    
	    Int_t lineColor = (zg1-zpos > 0) ? 2:3;
	    marker=new TMarker3DBox(xg1,yg1,zpos,xg2,yg2,0,0,0);
	    
	    marker->SetLineColor(lineColor);
	    marker->SetFillStyle(1001);
	    marker->SetFillColor(color);
	    marker->SetRefObject((TObject*)points);
	    points->Set3DMarker(0, marker);
	    
	} // end loop on digits
    } // end of new segmentation    
}
//___________________________________________
void AliMUONDisplay::LoadCoG(Int_t chamber, Int_t /*cathode*/)
{
/// Read raw clusters info and store x,y,z info in arrays fRpoints.
/// Loop on all detectors

    if (chamber > 10) return;
    
    ResetRpoints();
    
    GetMUONData()->SetTreeAddress("RC");
    TClonesArray *muonRawClusters  = GetMUONData()->RawClusters(chamber-1);

    if (muonRawClusters == 0) return;

    Int_t nent = 0;
    if (GetMUONData()->TreeR()) {
	nent=(Int_t) GetMUONData()->TreeR()->GetEntries();
	GetMUONData()->TreeR()->GetEvent(0);
    }
    
    Int_t nrawcl = muonRawClusters->GetEntriesFast();
    if (nrawcl == 0) return;
    if (fRpoints == 0) fRpoints = new TObjArray(nrawcl);
    
    Float_t zpos = AliMUONConstants::DefaultChamberZ(chamber-1);  
    AliMUONRawCluster  *mRaw;
    AliMUONPoints *points = 0;
    //
    //loop over all raw clusters and store their position
    points = new AliMUONPoints(nrawcl);
    for (Int_t iraw=0;iraw<nrawcl;iraw++) {
  	mRaw   = (AliMUONRawCluster*)muonRawClusters->UncheckedAt(iraw);
        points->SetMarkerColor(51);
        points->SetMarkerStyle(2);
        points->SetMarkerSize(1.);
        points->SetParticle(-1);
        points->SetHitIndex(-1);
        points->SetTrackIndex(-1);
        points->SetDigitIndex(-1);
        points->SetPoint(iraw,mRaw->GetX(0),mRaw->GetY(0),zpos);
	fRpoints->AddAt(points,iraw);
	//	printf("%f and %f and %f\n",mRaw->GetX(0),mRaw->GetY(0),mRaw->GetZ(0));
    }
}

//___________________________________________
void AliMUONDisplay::LoadTracks()
{
/// Load tracks

  AliMUONTrack* recTrack = 0;
  AliMUONTrackParam* trackParam = 0;
  TClonesArray *  trackParamAtHit = 0;   

  ResetRpoints();

  GetMUONData()->SetTreeAddress("RT");
  TClonesArray* recTracksArray = GetMUONData()->RecTracks();
  if (recTracksArray == NULL) return;
  GetMUONData()->GetRecTracks();

  Int_t nRecTracks = 0;
  if (recTracksArray)
    nRecTracks = (Int_t) recTracksArray->GetEntriesFast();


  if (fRpoints == 0) fRpoints = new TObjArray(nRecTracks);

  for (Int_t iRecTracks = 0; iRecTracks <  nRecTracks;  iRecTracks++) {
    // reading info from tracks
    recTrack = (AliMUONTrack*) recTracksArray->At(iRecTracks);

    Int_t nTrackHits = recTrack->GetNTrackHits();

    if (nTrackHits == 0) continue;

    Int_t iPoint = 0;
    TPolyLine3D *points = new TPolyLine3D(nTrackHits+1); 
    points->SetLineColor(6);
    points->SetLineWidth(1);
    fRpoints->AddAt(points,iRecTracks);

    Float_t xRec=0;
    Float_t yRec=0;
    Float_t zRec=0;

    trackParam = recTrack->GetTrackParamAtVertex(); 
    xRec  = trackParam->GetNonBendingCoor();
    yRec  = trackParam->GetBendingCoor();
    zRec  = trackParam->GetZ();
    points->SetPoint(iPoint,xRec,yRec,zRec);
    iPoint++;	

    for (Int_t iHit = 0; iHit < nTrackHits; iHit++){
      trackParamAtHit = recTrack->GetTrackParamAtHit();
      trackParam = (AliMUONTrackParam*) trackParamAtHit->At(iHit); 
      xRec  = trackParam->GetNonBendingCoor();
      yRec  = trackParam->GetBendingCoor();
      zRec  = trackParam->GetZ();
      points->SetPoint(iPoint,xRec,yRec,zRec);
      iPoint++;	
    } // end loop rec. hits
    PrintTrack(iRecTracks,recTrack);
  } // end loop tracks
  

}

//___________________________________________
void AliMUONDisplay::PrintTrack(Int_t iRecTracks, AliMUONTrack *recTrack)
{
/// Print reconstructed track

  AliMUONTrackParam *trackParam;
  Float_t vertex[3], momentum[3];
  Float_t pYZ, bendingSlope, nonBendingSlope, chi2dof;
  Int_t charge;

  trackParam = recTrack->GetTrackParamAtVertex();
  vertex[0] = trackParam->GetNonBendingCoor();
  vertex[1] = trackParam->GetBendingCoor();
  vertex[2] = trackParam->GetZ();
  pYZ =  1./TMath::Abs(trackParam->GetInverseBendingMomentum());
  bendingSlope = trackParam->GetBendingSlope();
  nonBendingSlope = trackParam->GetNonBendingSlope();
  momentum[2] = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);
  momentum[0] = momentum[2] * nonBendingSlope;
  momentum[1] = momentum[2] * bendingSlope;
  charge = Int_t(TMath::Sign(1.,trackParam->GetInverseBendingMomentum()));
  chi2dof = recTrack->GetFitFMin()/(2.0 * recTrack->GetNTrackHits() - 5.);
  
  printf("===================================================\n");
  printf(" Reconstructed track # %d \n",iRecTracks);
  printf(" charge: %d \n",charge);
  printf(" vertex x,y,z (cm): %f %f %f \n",vertex[0],vertex[1],vertex[2]); 
  printf(" momentum Px,Py,Pz (GeV/c): %f %f %f \n",momentum[0],momentum[1],momentum[2]);
  printf(" track chi2/dof: %f \n",chi2dof); 
  
}

//___________________________________________
void AliMUONDisplay::LoadHits(Int_t chamber)
{
/// Read hits info and store x,y,z info in arrays fPhits.
/// Loop on all detectors

    if (chamber > 14) return;
    Int_t track;

    fChamber=chamber;
 
    ResetPhits();
    
    Float_t zpos=AliMUONConstants::DefaultChamberZ(chamber-1);

    if (GetMUONData()->TreeH()) {
      GetMUONData()->SetTreeAddress("H");
      Int_t ntracks = (Int_t)GetMUONData()->TreeH()->GetEntries(); //skowron
      Int_t nthits  = 0;
      for (track = 0; track < ntracks; track++) {
	GetMUONData()->ResetHits();
	GetMUONData()->GetTrack(track);//skowron
	TClonesArray *muonHits  = GetMUONData()->Hits();
	if (muonHits == 0) return;
	nthits += muonHits->GetEntriesFast();
      } 
      if (fPhits == 0) fPhits = new TObjArray(nthits);
      Int_t nhold=0;
      for (track=0; track<ntracks;track++) {
	GetMUONData()->ResetHits();
	GetMUONData()->GetTrack(track);//skowron
	TClonesArray *muonHits  = GetMUONData()->Hits();
	if (muonHits == 0) return;
	Int_t nhits = muonHits->GetEntriesFast();
	if (nhits == 0) continue;
	AliMUONHit *mHit;
	AliMUONPoints *points = 0;
	Int_t npoints=1;
	for (Int_t hit=0;hit<nhits;hit++) {
	  mHit = (AliMUONHit*)muonHits->UncheckedAt(hit);
	  Int_t nch  = mHit->Chamber();              // chamber number
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
	  //	    printf("%f and %f and %f\n",mHit->X(),mHit->Y(),mHit->Z());
	}
	nhold+=nhits;
      }
    }
}

//_____________________________________________________________________________
void AliMUONDisplay::Paint(Option_t *)
{
/// Paint miscellaneous items
}

//_____________________________________________________________________________
void AliMUONDisplay::SetPickMode()
{
/// Set parameters for pick mode.
 
    fZoomMode = 0;

    fArcButton->SetY1(fPickButton->GetYlowNDC()+0.5*fPickButton->GetHNDC());
    fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliMUONDisplay::SetZoomMode()
{
/// Set parameters for zoom mode

    fZoomMode = 1;
    
    fArcButton->SetY1(fZoomButton->GetYlowNDC()+0.5*fZoomButton->GetHNDC());
    fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliMUONDisplay::NextChamber(Int_t delta)
{
/// To go from chamber to next chamber if delta = 1
/// or previous chamber otherwise
    if (delta == 1) {
	if (fChamber < AliMUONConstants::NCh()) fChamber++;
    } else {
	if (fChamber > 1) fChamber--;
    }
    if (!fPad) return;
    fPad->Clear();
    LoadDigits(fChamber, fCathode);
    DrawChamber();
}

//_____________________________________________________________________________
void AliMUONDisplay::NextCathode()
{
/// To switch to other cathode plane

    if (!fPad) return;
    fPad->Clear();
    if (fCathode == 1) {
	LoadDigits(fChamber, 2);	
    } else {
	LoadDigits(fChamber, 1);
    }
    fNextCathode = kTRUE; // to keep the same zoom
    DrawChamber();
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
/// Print global trigger output

  AliMUONGlobalTrigger* globalTrig;

  GetMUONData()->SetTreeAddress("GLT");
  GetMUONData()->GetTriggerD();

  globalTrig =  (AliMUONGlobalTrigger*)GetMUONData()->GlobalTrigger()->UncheckedAt(0);
  if (globalTrig == 0) return;

  globalTrig->Print("full");

 //  // returns Trigger Decision for current event
//   AliMUONTriggerDecision* decision= new AliMUONTriggerDecision(GetLoader(),1);

//   //  AliMUONTriggerDecision* decision= new AliMUONTriggerDecision(1);
//   AliMUONData* muonData = decision->GetMUONData();
//   muonData->SetTreeAddress("D");
//   decision->Trigger(); 
}
//_____________________________________________________________________________
void AliMUONDisplay::SetChamberAndCathode(Int_t chamber, Int_t cathode)
{
/// Set chamber and cathode number

   fChamber = chamber;
   fCathode = cathode;

   if (!fPad) return;
   fPad->Clear();
   LoadDigits(chamber,cathode);
   DrawChamber();
}

//_____________________________________________________________________________
void AliMUONDisplay::SetEvent(Int_t newevent)
{
/// Chose event 

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
/// Set view range along R and Z

    fRrange = rrange;
    fZrange = zrange;

    if (!fPad) return;
    fPad->Clear();
    Draw();
}
   
//_____________________________________________________________________________
void AliMUONDisplay::SetView(Float_t theta, Float_t phi, Float_t psi)
{
/// Change viewing angles for current event

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
/// Display (current event_number + delta)
///  -  delta =  1  shown next event
///  -  delta = -1 show previous event

  AliRunLoader * runLoader;
  if (fLoader)
    runLoader = fLoader->GetRunLoader();
  else
    runLoader = 0x0;
    
    if (delta) {
      //runLoader->CleanDetectors();
      //runLoader->CleanKinematics();
      Int_t currentEvent = runLoader->GetEventNumber();
      Int_t newEvent     = currentEvent + delta;
      runLoader->GetEvent(newEvent);
      fEvent=newEvent;
    }
    LoadDigits(fChamber, fCathode);
    fPad->cd(); 
    Draw();
}

//______________________________________________________________________________
void AliMUONDisplay::UnZoom()
{
/// Unzoom 

    if (fZooms <= 0) return;
    fZooms--;
    TPad *pad = (TPad*)gPad->GetPadSave();
    pad->Range(fZoomX0[fZooms],fZoomY0[fZooms], fZoomX1[fZooms],fZoomY1[fZooms]);
    pad->Modified();
}

//_____________________________________________________________________________
void AliMUONDisplay::ResetPoints()
{
/// Reset array of points

    if (fPoints) {
	fPoints->Delete();
	delete fPoints;
	fPoints = 0;
    }
}
//_____________________________________________________________________________
void AliMUONDisplay::ResetPhits()
{
/// Reset array of points

    if (fPhits) {
	fPhits->Delete();
	delete fPhits;
	fPhits = 0;
    }
}
//_____________________________________________________________________________
void AliMUONDisplay::ResetRpoints()
{
/// Reset array of points

    if (fRpoints) {
	fRpoints->Clear();
	//	delete fRpoints;
	fRpoints = 0;
    }
}
