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
  Revision 1.7  2000/10/02 21:28:12  fca
  Removal of useless dependecies via forward declarations

  Revision 1.6  2000/10/02 15:46:38  jbarbosa
  Fixed forward declarations.

  Revision 1.5  2000/06/30 16:49:34  dibari
  Different call for ring drawing.

  Revision 1.4  2000/06/12 15:21:08  jbarbosa
  Cleaned up version.

  Revision 1.3  2000/06/09 14:52:08  jbarbosa
  New tentative ellipse drawing routine

  Revision 1.1  2000/04/19 13:07:45  morsch
  Digits, clusters and reconstruction results added.

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
#include <TAtt3D.h>
#include <TAttLine.h>
#include <TPolyMarker.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TList.h>
#include <TDiamond.h>
#include <TNode.h>
#include <TArc.h>
#include <TTUBE.h>
#include <TSlider.h>
#include <TGeometry.h>
#include <TSliderBox.h>
#include <TGaxis.h>
#include <TVirtualX.h>
#include <TMath.h>
#include <TRandom.h>
#include <X3DBuffer.h>
#include <TParticle.h>

#include "AliRun.h"
#include "AliPDG.h"
#include "AliDetector.h"
#include "AliRICH.h"
#include "AliRICHConst.h"
#include "AliRICHDisplay.h"
#include "AliRICHPoints.h"

#include "AliRICHHit.h"
#include "AliRICHCerenkov.h"
#include "AliRICHPadHit.h"
#include "AliRICHDigit.h"
#include "AliRICHRawCluster.h"
#include "AliRICHRecHit.h"
#include "AliRICHEllipse.h"

ClassImp(AliRICHDisplay)
    

//____________________________________________________________________________
AliRICHDisplay::AliRICHDisplay()
{ 

// default constructor

    fPoints = 0;
    fPhits = 0;
    fPCerenkovs = 0;
    fCanvas = 0;
    fRpoints = 0;
    fRecpoints = 0; 
}

//_____________________________________________________________________________
AliRICHDisplay::AliRICHDisplay(Int_t size)
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
  <img src="gif/AliRICHDisplay.gif">
*/
//End_Html
    
    
    fPad = 0;
    
    gAlice->SetDisplay(this);
    
    // Initialize display default parameters
    SetRange();
    
    // Set front view by default
    fTheta = 90;
    fPhi   = 90;
    fPsi   = 0;
    fChamber = 1;
    fCathode = 1;
    //   fRzone   = 1.e10;
    fDrawClusters  = kTRUE;
    fDrawCoG       = kTRUE;
    fDrawRecHits   = kTRUE;
    fZoomMode      = 1;
    fZooms         = 0;
    fClustersCuts  = 0;
    fPoints        = 0;
    fPCerenkovs    = 0;
    fPhits         = 0;
    fRpoints       = 0;   
    fRecpoints       = 0;
    // Create colors
    CreateColors();
    // Create display canvas
    Int_t ysize = size;
    if (ysize < 100) ysize = 750;
    Int_t xsize = Int_t(size*830./ysize);
    fCanvas = new TCanvas("Canvas", "RICH Clusters Display",14,47,xsize,ysize);
    fCanvas->ToggleEventStatus();
    
    // Create main display pad
    fPad = new TPad("viewpad", "RICH display",0.15,0,0.9,1);
    fPad->Draw();
    fPad->Modified();
    fPad->SetFillColor(1);
    fPad->SetBorderSize(2);
    
    fCanvas->cd();
    
    // Create colors pad
    fColPad = new TPad("colpad", "Colors pad",0.9,0,1,1);
    fColPad->Draw();
    fColPad->Modified();
    fColPad->SetFillColor(19);
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
}


//_____________________________________________________________________________
AliRICHDisplay::~AliRICHDisplay()
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
    if (fRecpoints) fRecpoints->Delete();
    delete fRecpoints;
    fRecpoints     = 0;
    //
    if (fPCerenkovs) fPCerenkovs->Delete();
    delete fPCerenkovs;
    fPCerenkovs     = 0;
}

//_____________________________________________________________________________
void AliRICHDisplay::Clear(Option_t *)
{
//    Delete graphics temporary objects
}

//_____________________________________________________________________________
void AliRICHDisplay::DisplayButtons()
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
    button = new TButton("Next",but1,x0,y-dbutton,x1,y);
    button->SetFillColor(38);
    button->Draw();
    
    y -= dbutton +dy;
    char but2[] = "gAlice->Display()->ShowNextEvent(-1)";
    button = new TButton("Previous",but2,x0,y-dbutton,x1,y);
    button->SetFillColor(38);
    button->Draw();

    y -= dbutton +dy;
    char but7[] = "gAlice->Display()->DrawViewGL()";
    button = new TButton("OpenGL",but7,x0,y-dbutton,x1,y);
    button->SetFillColor(38);
    button->Draw();
    
    y -= dbutton +dy;
    char but8[] = "gAlice->Display()->DrawViewX3D()";
    button = new TButton("X3D",but8,x0,y-dbutton,x1,y);
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
    diamond->AddText("RICH");
    diamond->AddText("... ");
    diamond->AddText(" ");
}

//_____________________________________________________________________________
void AliRICHDisplay::CreateColors()
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
		color=700+23-color;
		new TColor(color,r,g,b);
	    } 
	    break;
	case 2:
	    for (i=1;i<=4;i++) {
		r=1.1-i*0.2;
		g=1.;  
		b=0.;
		color=i+5;
		color=700+23-color;
		new TColor(color,r,g,b);
	    } 
	    break;
    case 3:
	for (i=1;i<=4;i++) {
	    r=0.;
	    g=1.;  
	    b=i*0.2+0.2;
	    color=i+9;
	    color=700+23-color;
	    new TColor(color,r,g,b);
	} 
	break;
	case 4:
	    for (i=1;i<=4;i++) {
		r=0.;
		g=1.1-i*0.2;  
		b=1.;
		color=i+13;
		color=700+23-color;
		new TColor(color,r,g,b);
	    } 
	    break;
	case 5:
	    for (i=1;i<=5;i++) {
		r=i*0.2;
		g=0.;  
		b=1.;
		color=i+17;
		color=700+23-color;
		new TColor(color,r,g,b);
	    } 
	    break;
	}
    }
}

//_____________________________________________________________________________
void AliRICHDisplay::DisplayColorScale()
{

// Draw the color scale in the RICH display canvas
    
    Int_t i;
    Int_t color;
    Float_t xlow, ylow, xup, yup, hs;
    Float_t x1, y1, x2, y2;
    x1 = y1 = 0;
    x2 = y2 = 20;
    
    gPad->SetFillColor(0);
    gPad->Clear();
    gPad->Range(x1,y1,x2,y2);
    TText *text = new TText(0,0,"");
    text->SetTextFont(61);
    text->SetTextSize(0.03);
    text->SetTextAlign(22);
    
    TBox *box;
    char label[8];
//*-* draw colortable boxes
    hs = (y2-y1)/Float_t(22);
    xlow=x1+1;
    xup=x2-9;
    for (i=0;i<22;i++) {
	ylow = y1 + hs*(Float_t(i));
	yup  = y1 + hs*(Float_t(i+1));
	color = 701+i;
	Double_t logscale=Double_t(i+1)*(TMath::Log(adc_satm)/22);
	Int_t scale=(Int_t)TMath::Exp(logscale);
	sprintf(label,"%d",scale);
	box = new TBox(xlow, ylow, xup, yup);
	box->SetFillColor(color);
	box->Draw();
	text->DrawText(xup+4, 0.5*(ylow+yup),label);
    }
}

//______________________________________________________________________________
Int_t AliRICHDisplay::DistancetoPrimitive(Int_t px, Int_t)
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
void AliRICHDisplay::Draw(Option_t *)
{
//    Display current event

    fPad->cd();
    
    DrawView(fTheta, fPhi, fPsi);      // see how to draw PGON+inner frames
    
    // Display the event number and title
    fPad->cd();
    DrawTitle();
}

//_____________________________________________________________________________
void AliRICHDisplay::DrawCoG()
{
//    Draw hits for RICH chambers

    if (!fDrawCoG) return;
    ResetRpoints();
    for (Int_t chamber=0;chamber<kNCH;chamber++) {
	LoadCoG(chamber,1);
    }
    
    Int_t ncog, icog;
    TObjArray *points;
    AliRICHPoints *pm;
    points = fRpoints;
    if (!points) return;
    ncog = points->GetEntriesFast();
    for (icog=0; icog < ncog; icog++) {
	pm = (AliRICHPoints*)points->UncheckedAt(icog);
	if (!pm) continue;
	pm->Draw();
    }
}

void AliRICHDisplay::DrawRecHits()
{
//    Draw rec hits for RICH chambers

    if (!fDrawRecHits) return;
    //ResetRecpoints();
    for (Int_t chamber=0;chamber<kNCH;chamber++) {
	LoadRecHits(chamber,1);
    }
    
    Int_t nrec, irec;
    TObjArray *points;
    AliRICHPoints *pm;
    points = fRecpoints;
    if (!points) return;
    nrec = points->GetEntriesFast();
    printf("Nrec %d\n",nrec);
    for (irec=0; irec < nrec; irec++) {
	pm = (AliRICHPoints*)points->UncheckedAt(irec);
	if (!pm) continue;
	pm->Draw();
    }
}

//_____________________________________________________________________________

void AliRICHDisplay::DrawCerenkovs()
{
//    Draw cerenkovs hits for RICH chambers
    
    LoadCerenkovs(fChamber);
    //printf("\nDrawCerenkovs\n");
    
    Int_t ntracks, track;
    TObjArray *cpoints;
    AliRICHPoints *pm;
    
    fHitsCuts = 0;
    cpoints = fPCerenkovs;
    //printf ("Cpoints: %p",cpoints);
    if (!cpoints) return;
    ntracks = cpoints->GetEntriesFast();
    //printf("DrawCerenkovs - ntracks %d \n",ntracks);
    for (track=0;track<ntracks;track++) {
	pm = (AliRICHPoints*)cpoints->UncheckedAt(track);
	if (!pm) continue;
	pm->Draw();
	fHitsCuts += pm->GetN();
    }
}

//_____________________________________________________________________________

void AliRICHDisplay::DrawClusters()
{
//    Draw clusterss for RICH chambers
    
    Int_t ndigits, digit;
    TObjArray *points;
    AliRICHPoints *pm;
    
    fClustersCuts = 0;
    points = fPoints;
    if (!points) return;
    ndigits = points->GetEntriesFast();
    for (digit=0;digit<ndigits;digit++){
	pm = (AliRICHPoints*)points->UncheckedAt(digit);
	if (!pm) continue;
	pm->Draw();
	Float_t *pxyz;
	pxyz=pm->GetP();
	for (Int_t im=0;im<3;im++) {
	    TMarker3DBox *marker=pm->GetMarker(im);
	    if (marker)
		marker->Draw();
	}
	fClustersCuts +=pm->GetN();
    }
}

//_____________________________________________________________________________
void AliRICHDisplay::DrawHits()
{
//    Draw hits for RICH chambers
    
    LoadHits(fChamber);
    
    Int_t ntracks, track;
    TObjArray *points;
    AliRICHPoints *pm;
    
    fHitsCuts = 0;
    points = Phits();
    if (!fPhits) return;
//    ntracks = points->GetEntriesFast();
    ntracks = fPhits->GetEntriesFast();
    
    //printf("DrawHits - ntracks %d \n",ntracks);
    for (track=0;track<ntracks;track++) {
	pm = (AliRICHPoints*)fPhits->UncheckedAt(track);
	if (!pm) continue;
	pm->Draw();
	fHitsCuts += pm->GetN();
    }
}


//_____________________________________________________________________________
void AliRICHDisplay::DrawTitle(Option_t *option)
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
	title->SetBit(kCanDelete);
	title->SetFillColor(42);
	title->Draw();
	char ptitle[100];
	sprintf(ptitle,"Alice event: %d, Run:%d",
		gAlice->GetHeader()->GetEvent(), gAlice->GetHeader()->GetRun());
	title->AddText(ptitle);
	Int_t nparticles = gAlice->Particles()->GetEntriesFast();
	sprintf(ptitle,"Nparticles = %d Nhits = %d Npads fired = %d",
		nparticles, fHitsCuts,fClustersCuts);
	title->AddText(ptitle);
    } else {
	TPaveLabel *label = 
	    new TPaveLabel(xmin +0.01*dx, ymax-0.07*dy, xmin +0.2*dx, ymax-0.01*dy,option);
	label->SetBit(kCanDelete);
	label->SetFillColor(42);
	label->Draw();
    }
}

//_____________________________________________________________________________
void AliRICHDisplay::DrawView(Float_t theta, Float_t phi, Float_t psi)
{
//    Draw a view of RICH clusters
    
   gPad->SetCursor(kWatch);
   gPad->SetFillColor(1);
   gPad->Clear();
   
   Int_t iret;
   TView *view = new TView(1);
   Float_t range = fRrange*fRangeSlider->GetMaximum();
   view->SetRange(-range,-range,-range,range, range, range);
   fZoomX0[0] = -1;
   fZoomY0[0] = -1;
   fZoomX1[0] =  1;
   fZoomY1[0] =  1;
   fZooms = 0;
   
   //Display RICH Chamber Geometry
   gAlice->GetGeometry()->Draw("same");
   
   //add clusters to the pad
   DrawClusters();
   DrawHits();
//   DrawCerenkovs();
   printf("Calling DrawCoG\n");
   DrawCoG();
   printf("Calling DrawRecHits\n");
   DrawRecHits();
   /*for (Int_t i=0;i<7;i++)
     LoadRecHits(i,1);*/
   
   // add itself to the list (must be last)
   AppendPad();
   
   view->SetView(phi, theta, psi, iret);
}

//______________________________________________________________________________
void AliRICHDisplay::ExecuteEvent(Int_t event, Int_t px, Int_t py)
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
       gVirtualX->DrawBox(px0, py0, pxold, pyold, TVirtualX::kHollow);
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
       if (fZooms < kMAXZOOM-1) {
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
void AliRICHDisplay::LoadCoG(Int_t chamber, Int_t cathode)
{
// Read raw clusters info and store x,y,z info in arrays fRpoints
// Loop on all detectors

   if (chamber > 6) return;

   printf("Entering LoadCoG\n");


   AliRICH *pRICH  = (AliRICH*)gAlice->GetModule("RICH");
   AliRICHChamber*  iChamber;

   TClonesArray *pRICHrawclust  = pRICH->RawClustAddress(chamber);
   printf ("Chamber:%d has adress:%p\n", chamber, pRICHrawclust );
   if (pRICHrawclust == 0) return;

   pRICH->ResetRawClusters();


   Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
   gAlice->TreeR()->GetEvent(nent-1+cathode-1);
   Int_t nrawcl = pRICHrawclust->GetEntriesFast();
   printf ("nrawcl:%d\n",nrawcl);
   if (nrawcl == 0) return;
   if (fRpoints == 0) fRpoints = new TObjArray(nrawcl);
   
   iChamber = &(pRICH->Chamber(chamber));
   AliRICHRawCluster  *mRaw;
   AliRICHPoints *points = 0;
   //
   //loop over all raw clusters and store their position
   points = new AliRICHPoints(nrawcl);
   for (Int_t iraw=0;iraw<nrawcl;iraw++) {
       mRaw   = (AliRICHRawCluster*)pRICHrawclust->UncheckedAt(iraw);
       fRpoints->AddAt(points,iraw);
       points->SetMarkerColor(3);
       points->SetMarkerStyle(3);
       points->SetMarkerSize(1.);
       points->SetParticle(-1);
       points->SetHitIndex(-1);
       points->SetTrackIndex(-1);
       points->SetDigitIndex(-1);
       Float_t  vectorLoc[3]={mRaw->fX,6.276,mRaw->fY};
       Float_t  vectorGlob[3];
       iChamber->LocaltoGlobal(vectorLoc,vectorGlob);
       points->SetPoint(iraw,vectorGlob[0],vectorGlob[1],vectorGlob[2]);
   }
}
//___________________________________________
void AliRICHDisplay::LoadRecHits(Int_t chamber, Int_t cathode)
{
// Read rec. hits info 
// Loop on all detectors

   if (chamber > 6) return;

   printf("Entering LoadRecHits\n");


   AliRICH *pRICH  = (AliRICH*)gAlice->GetModule("RICH");
   AliRICHChamber*  iChamber;

   TClonesArray *pRICHrechits  = pRICH->RecHitsAddress(chamber);
   printf ("Chamber:%d\n", chamber);
   if (pRICHrechits == 0) return;

   //RICH->ResetRecHits();


   Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
   gAlice->TreeR()->GetEvent(nent-1+cathode-1);
   Int_t nrechits = pRICHrechits->GetEntriesFast();
   printf ("nrechits:%d\n",nrechits);
   if (nrechits == 0) return;
   if (fRecpoints == 0) fRecpoints = new TObjArray(50);
   
   iChamber = &(pRICH->Chamber(chamber));
   AliRICHRecHit  *mRec;
   AliRICHPoints *points = 0;
   //AliRICHEllipse *ellipse = 0;
   //
   //loop over all rechits and store their position  

   points = new AliRICHPoints(nrechits);
   for (Int_t irec=0;irec<nrechits;irec++) {
       mRec   = (AliRICHRecHit*)pRICHrechits->UncheckedAt(irec);
       fRecpoints->AddAt(points,irec);
       points->SetMarkerColor(38);
       points->SetMarkerStyle(8);
       points->SetMarkerSize(1.);
       points->SetParticle(-1);
       points->SetHitIndex(-1);
       points->SetTrackIndex(-1);
       points->SetDigitIndex(-1);
       Float_t  vectorLoc[3]={mRec->fX,6.276,mRec->fY};
       Float_t  vectorGlob[3];
       iChamber->LocaltoGlobal(vectorLoc,vectorGlob);
       points->SetPoint(irec,vectorGlob[0],vectorGlob[1],vectorGlob[2]);
       //Float_t theta = iChamber->GetRotMatrix()->GetTheta();
       //Float_t phi   = iChamber->GetRotMatrix()->GetPhi();	   
       //ellipse=new TEllipse(vectorGlob[0],vectorGlob[2],10,10,0,360,phi);
       printf("Generating ellipse %d\n",irec);
       AliRICHEllipse *ellipse=new AliRICHEllipse(mRec->fX,mRec->fY,mRec->fOmega,mRec->fTheta,mRec->fPhi,mRec->fEmissPoint);
       ellipse->CerenkovRingDrawing(chamber,irec);
       //ellipse->SetFillStyle(1001);
       ellipse->SetMarkerColor(38);
       ellipse->Draw();
       //marker->SetRefObject((TObject*)points);
       //points->Set3DMarker(0, marker); 
   }
}
//___________________________________________
void AliRICHDisplay::LoadDigits()
{
// Read digits info and store x,y,z info in arrays fPoints
// Loop on all detectors
    
   ResetPoints();
   AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
   AliRICHChamber*       iChamber;
   AliSegmentation*      segmentation;
   Int_t nAllDigits=0;
   Int_t ich;
   
   for (ich=0; ich<kNCH; ich++) {
       TClonesArray *pRICHdigits  = pRICH->DigitsAddress(ich);
       if (pRICHdigits == 0) continue;
       gAlice->ResetDigits();
       gAlice->TreeD()->GetEvent(1);
       Int_t ndigits = pRICHdigits->GetEntriesFast();
        nAllDigits+=ndigits;
   }
   if (fPoints == 0) fPoints = new TObjArray(nAllDigits);   
   Int_t counter=0;
   for (ich=0; ich<kNCH; ich++) {
       TClonesArray *pRICHdigits  = pRICH->DigitsAddress(ich);
       if (pRICHdigits == 0) continue;
       gAlice->ResetDigits();
       gAlice->TreeD()->GetEvent(1);
       Int_t ndigits = pRICHdigits->GetEntriesFast();
       if (ndigits == 0) continue;
       iChamber = &(pRICH->Chamber(ich));
       segmentation=iChamber->GetSegmentationModel();
       Float_t dpx  = segmentation->Dpx();
       Float_t dpy  = segmentation->Dpy();

       //printf("Dpx:%d, Dpy:%d\n",dpx,dpy);
       
       AliRICHDigit  *mdig;
       AliRICHPoints *points = 0;
       TMarker3DBox  *marker = 0;
       //
       //loop over all digits and store their position
       Int_t npoints=1;
       
       for (Int_t digit=0;digit<ndigits;digit++) {
	   mdig    = (AliRICHDigit*)pRICHdigits->UncheckedAt(digit);
	   points = new AliRICHPoints(npoints);
	   fPoints->AddAt(points,counter);
	   counter++;
	   Int_t charge=mdig->fSignal;
	   Int_t index=Int_t(TMath::Log(charge)/(TMath::Log(adc_satm)/22));
	   Int_t color=701+index;
	   if (color>722) color=722;
	   points->SetMarkerColor(color);
	   points->SetMarkerStyle(21);
	   points->SetMarkerSize(0.5);
	   Float_t xpad, ypad, zpad;
	   segmentation->GetPadC(mdig->fPadX, mdig->fPadY,xpad, ypad, zpad);
	   Float_t vectorLoc[3]={xpad,6.276,ypad};
	   Float_t  vectorGlob[3];
	   iChamber->LocaltoGlobal(vectorLoc,vectorGlob);
	   points->SetParticle(-1);
	   points->SetHitIndex(-1);
	   points->SetTrackIndex(-1);
	   points->SetDigitIndex(digit);
	   points->SetPoint(0,vectorGlob[0],vectorGlob[1],vectorGlob[2]);
	   
	   segmentation->GetPadC(mdig->fPadX, mdig->fPadY, xpad, ypad, zpad);
	   Float_t theta = iChamber->GetRotMatrix()->GetTheta();
	   Float_t phi   = iChamber->GetRotMatrix()->GetPhi();	   
	   marker=new TMarker3DBox(vectorGlob[0],vectorGlob[1],vectorGlob[2],
				   dpy/2,0,dpx/2,theta,phi);
	   marker->SetLineColor(2);
	   marker->SetFillStyle(1001);
	   marker->SetFillColor(color);
	   marker->SetRefObject((TObject*)points);
	   points->Set3DMarker(0, marker);
       } // loop over digits
   } // loop over chambers 
}


//___________________________________________
void AliRICHDisplay::LoadHits(Int_t chamber)
{
// Read hits info and store x,y,z info in arrays fPhits
// Loop on all detectors
    
    
    fChamber=chamber; 
    ResetPhits();
    
    AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
    AliRICHChamber*  iChamber;
    
    iChamber = &(pRICH->Chamber(chamber-1));
    Int_t ntracks = (Int_t)gAlice->TreeH()->GetEntries();
    Int_t track;
    
    if (fPhits == 0) fPhits = new TObjArray(ntracks);
    //TVector *xp = new TVector(1000);
    //TVector *yp = new TVector(1000);
    //TVector *zp = new TVector(1000);
    //TVector *ptrk = new TVector(1000);
    //TVector *phit = new TVector(1000);
    Int_t nAllHits=0;
    for (track=0; track<ntracks;track++) {
	gAlice->ResetHits();
	gAlice->TreeH()->GetEvent(track);
	TClonesArray *pRICHhits  = pRICH->Hits();
	if (pRICHhits == 0) return;
	Int_t nhits = pRICHhits->GetEntriesFast();
	nAllHits+=nhits;
    }

    fPhits = new TObjArray(nAllHits);

    Int_t npoints=0;
    for (track=0; track<ntracks;track++) {
	gAlice->ResetHits();
	gAlice->TreeH()->GetEvent(track);
	TClonesArray *pRICHhits  = pRICH->Hits();
	if (pRICHhits == 0) return;
	Int_t nhits = pRICHhits->GetEntriesFast();
	if (nhits == 0) continue;
	AliRICHHit *mHit;
	AliRICHPoints *points = 0;
	for (Int_t hit=0;hit<nhits;hit++) {
	    points = new AliRICHPoints(1);
	    fPhits->AddAt(points,npoints);
            mHit = (AliRICHHit*)pRICHhits->UncheckedAt(hit);
	    TParticle *current = 
		(TParticle*)(*gAlice->Particles())[mHit->Track()];
	    if (current->GetPdgCode() == 50000050) {
		points->SetMarkerColor(kBlue);
	    } else if (current->GetPdgCode() == 50000051) {
		points->SetMarkerColor(kYellow);
	    } else {
		points->SetMarkerColor(kRed);
	    }
            points->SetMarkerStyle(5);
            points->SetMarkerSize(1.);
            points->SetParticle(mHit->Track());
	    points->SetHitIndex(hit);
            points->SetTrackIndex(track);
            points->SetDigitIndex(-1);
            points->SetPoint(hit,mHit->X(), mHit->Y(), mHit->Z());
	    npoints++;
	}
    }
}

//_____________________________________________________________________________

void AliRICHDisplay::LoadCerenkovs(Int_t chamber)
{
// Read cerenkov hits info and store x,y,z info in array fPCerenkovs
// Loop on all detectors
    
    fChamber=chamber; 
    ResetPCerenkovs();
    
    AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
    AliRICHChamber*  iChamber;
    
    iChamber = &(pRICH->Chamber(chamber-1));
    
    pRICH->SetTreeAddress();
    Int_t ntracks = (Int_t)gAlice->TreeH()->GetEntries();
    
    if (fPCerenkovs == 0) fPCerenkovs = new TObjArray(ntracks);
    TVector *xp = new TVector(1000);
    TVector *yp = new TVector(1000);
    TVector *zp = new TVector(1000);
    TVector *ptrk = new TVector(1000);
    TVector *phit = new TVector(1000);
    for (Int_t track=0; track<ntracks;track++) {
	gAlice->ResetHits();
	gAlice->TreeH()->GetEvent(track);
	TClonesArray *pRICHCerenkovs  = pRICH->Cerenkovs();
	if (pRICHCerenkovs == 0) return;
	Int_t nhits = pRICHCerenkovs->GetEntriesFast();
	if (nhits == 0) continue;
	AliRICHCerenkov *mCerenkov;
	AliRICHPoints *cpoints = 0;
	Int_t npoints=0;
	
	
//Display Cerenkov hits in blue
	
	for (Int_t hit=0;hit<nhits;hit++) {
            mCerenkov = (AliRICHCerenkov*)pRICHCerenkovs->UncheckedAt(hit);
	    (*xp)(npoints)=mCerenkov->X();
            (*yp)(npoints)=mCerenkov->Y();
            (*zp)(npoints)=mCerenkov->Z();
            (*ptrk)(npoints)=Float_t(mCerenkov->GetTrack());
            (*phit)(npoints)=Float_t(hit);
            npoints++;
	}
	if (npoints == 0) continue;
	cpoints = new AliRICHPoints(npoints);
	for (Int_t p=0;p<npoints;p++) {
	    cpoints->SetMarkerColor(kBlue);
	    cpoints->SetMarkerStyle(5);
	    cpoints->SetMarkerSize(1.);
	    cpoints->SetParticle(Int_t((*ptrk)(p)));
	    cpoints->SetHitIndex(Int_t((*phit)(p)));
	    cpoints->SetTrackIndex(track);
	    cpoints->SetDigitIndex(-1);
	    cpoints->SetPoint(p,(*xp)(p),(*yp)(p),(*zp)(p));
	}
	xp->Zero();
	yp->Zero();
	ptrk->Zero();
	phit->Zero();
	fPCerenkovs->AddAt(cpoints,track);
    }
}

//_____________________________________________________________________________
void AliRICHDisplay::Paint(Option_t *)
{
//    Paint miscellaneous items

}

//_____________________________________________________________________________
void AliRICHDisplay::SetPickMode()
{

// Toggle pick mode

    fZoomMode = 0;
    
    fArcButton->SetY1(fPickButton->GetYlowNDC()+0.5*fPickButton->GetHNDC());
    fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliRICHDisplay::SetZoomMode()
{

// Toggle Zoom mode

    fZoomMode = 1;
    
    fArcButton->SetY1(fZoomButton->GetYlowNDC()+0.5*fZoomButton->GetHNDC());
    fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliRICHDisplay::SetChamberAndCathode(Int_t chamber, Int_t cathode)
{
// Set chamber and cathode number
    fChamber = chamber;
    fCathode = cathode;
    
    printf("SetChamberAndCathode - fChamber fCathode %d %d\n",fChamber,fCathode);
    if (!fPad) return;
    fPad->Clear();
    LoadDigits();
    Draw();
}

//_____________________________________________________________________________
void AliRICHDisplay::SetRange(Float_t rrange, Float_t zrange)
{
// Set view range along R and Z
    fRrange = rrange;
    fZrange = zrange;
    
   if (!fPad) return;
   fPad->Clear();
   Draw();
}

//_____________________________________________________________________________
void AliRICHDisplay::SetView(Float_t theta, Float_t phi, Float_t psi)
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
void AliRICHDisplay::ShowNextEvent(Int_t delta)
{
//  Display (current event_number+delta)
//    delta =  1  shown next event
//    delta = -1 show previous event
    
    if (delta) {
	gAlice->Clear();
	Int_t currentEvent = gAlice->GetHeader()->GetEvent();
	Int_t newEvent     = currentEvent + delta;
	gAlice->GetEvent(newEvent);
	if (!gAlice->TreeD()) return; 
    }
    LoadDigits();
    DrawClusters();
    fPad->cd(); 
    Draw();

  
}

//______________________________________________________________________________
void AliRICHDisplay::UnZoom()
{

// Return to previous zoom factor

    if (fZooms <= 0) return;
    fZooms--;
    TPad *pad = (TPad*)gPad->GetPadSave();
    pad->Range(fZoomX0[fZooms],fZoomY0[fZooms], fZoomX1[fZooms],fZoomY1[fZooms]);
    pad->Modified();
}

//_____________________________________________________________________________
void AliRICHDisplay::ResetPoints()
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
void AliRICHDisplay::ResetRpoints()
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
void AliRICHDisplay::ResetRecpoints()
{
  //
  // Reset array of points
  //
  if (fRecpoints) {
    fRecpoints->Delete();
    delete fRecpoints;
    fRecpoints = 0;
  }
}
//_____________________________________________________________________________
void AliRICHDisplay::ResetPhits()
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
void AliRICHDisplay::ResetPCerenkovs()
{
    //
    // Reset array of points
    //
    if (fPCerenkovs) {
	fPCerenkovs->Delete();
	delete fPCerenkovs;
	fPCerenkovs = 0;
    }
}

//_____________________________________________________________________________

void AliRICHDisplay::DrawViewGL()
{
//    Draw current view using OPENGL

   TPad *pad = (TPad*)gPad->GetPadSave();
   pad->cd();
   TView *view = pad->GetView();
   if (!view) return;
   pad->x3d("OPENGL");
}

//_____________________________________________________________________________
void AliRICHDisplay::DrawViewX3D()
{
//    Draw current view using X3D

   TPad *pad = (TPad*)gPad->GetPadSave();
   pad->cd();
   TView *view = pad->GetView();
   if (!view) return;
   pad->x3d();
}
