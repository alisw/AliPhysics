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
Revision 1.10  1999/11/10 07:37:04  fca
Pads do not inherit editability from canvas any more

Revision 1.9  1999/11/09 07:38:51  fca
Changes for compatibility with version 2.23 of ROOT

Revision 1.8  1999/09/29 09:24:23  fca
Introduction of the Copyright and cvs Log

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
#include <TStyle.h>
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
#include <TH2.h>
#include <TObjArray.h>
#include <X3DBuffer.h>

#include "AliRun.h"
#include "AliDetector.h"
#include "AliITS.h"
#include "AliITSdigitNew.h"
#include "AliITSMap.h"
#include "AliITSresponseV0.h"
#include "AliITSdisplay.h"
#include "AliITSpoints.h"
#include "TParticle.h"


const  Int_t kMAXHIST = 20; 

static Int_t  sModule=0; 


ClassImp(AliITSdisplay)


//_____________________________________________________________________________
AliITSdisplay::AliITSdisplay()
{
   fPoints = 0;
   fPhits = 0;
   fRpoints = 0;
   fR2points = 0;
   fCpoints = 0;
   fCanvas = 0;

   fMap = 0;

}

//_____________________________________________________________________________
AliITSdisplay::AliITSdisplay(Int_t size)
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
<img src="gif/aliITSdisplay.gif">
*/
//End_Html


   fPad = 0;

   gAlice->SetDisplay(this);
   
   // Set module view by default
   fTheta = -90;
   fPhi   =  90;
   fPsi   =  180;

   fYmodule = 0;

   fModule = 393;
   fNmodules = 0;

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   AliITSgeom *gm=ITS->GetITSgeom();
   gm->GetModuleId(fModule,fLayer,fLadder,fDetector);   

   //   fRzone   = 1.e10;
   fDrawClusters  = kTRUE;
   fDrawCoG       = kTRUE;
   fDrawCathCor  = kFALSE;
   fDrawHist  = kFALSE;

   fMap = 0;


   fZoomMode      = 1;
   fZooms         = 0;
   fClustersCuts  = 0;
   fPoints        = 0;
   fPhits         = 0;
   fRpoints       = 0;
   fR2points = 0;
   fCpoints = 0;

   /*
   // Initialize display default parameters
   SetRange(3.5,3.5);
   // Create colors
   CreateColors();
   // Create display canvas
   Int_t ysize = size;
   if (ysize < 100) ysize = 750;
   Int_t xsize = Int_t(size*830./ysize);
   fCanvas = new TCanvas("ModuleCanvas", "ITS Clusters Display",14,47,xsize,ysize);
   fCanvas->ToggleEventStatus();
   
   // Create main display pad
   fPad = new TPad("mviewpad", "ITS display",0.15,0,0.9,1);
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
   fTrigPad = new TPad("mtrigger", "range and mode pad",0,0,dxtr,dytr);
   fTrigPad->SetEditable(kFALSE);
   fTrigPad->Draw();
   fTrigPad->cd();
   fTrigPad->SetFillColor(22);
   fTrigPad->SetBorderSize(2);
   fRangeSlider = new TSlider("mrange","range",0.7,0.42,0.9,0.98);
   fRangeSlider->SetObject(this);
   //char pickmode[] = "gAlice->Display()->SetPickMode()";
   char pickmode[] = "((AliITSdisplay*)(gAlice->Display()))->SetPickMode()";
   Float_t db = 0.09;
   fPickButton = new TButton("Pick",pickmode,0.05,0.32,0.65,0.32+db);
   fPickButton->SetFillColor(38);
   fPickButton->Draw();
   //char zoommode[] = "gAlice->Display()->SetZoomMode()";
   char zoommode[] = "((AliITSdisplay*)(gAlice->Display()))->SetZoomMode()";
   fZoomButton = new TButton("Zoom",zoommode,0.05,0.21,0.65,0.21+db);
   fZoomButton->SetFillColor(38);
   fZoomButton->Draw();
   fArcButton = new TArc(.8,fZoomButton->GetYlowNDC()+0.5*db,0.33*db);
   fArcButton->SetFillColor(kGreen);
   fArcButton->Draw();
   //char butUnzoom[] = "gAlice->Display()->UnZoom()";
   char butUnzoom[] = "((AliITSdisplay*)(gAlice->Display()))->UnZoom()";
   TButton *button = new TButton("UnZoom",butUnzoom,0.05,0.05,0.95,0.15);
   button->SetFillColor(38);
   button->Draw();
   AppendPad(); // append display object as last object to force selection

   fCanvas->cd();
   fCanvas->Update();
   */

   CreateModuleCanvas(size);

}

//_____________________________________________________________________________
void AliITSdisplay::CreateModuleCanvas(Int_t size)
{

   // Initialize display default parameters
   SetRange(3.5,3.5);
   // Create colors
   CreateColors();
   // Create display canvas
   Int_t ysize = size;
   if (ysize < 100) ysize = 750;
   Int_t xsize = Int_t(size*830./ysize);
   fCanvas = new TCanvas("ModuleCanvas", "ITS Clusters Display",14,47,xsize,ysize);
   fCanvas->ToggleEventStatus();

    Int_t wtopx, wtopy;
    UInt_t ww,wh;
    fCanvas->GetCanvasPar(wtopx,wtopy,ww,wh);

      Float_t cx = gStyle->GetScreenFactor();

    printf("CreateCanvas: cx wtopx wtopy ww wh %f %d %d %d %d\n", cx,wtopx,wtopy,ww,wh);
 
      Int_t dum1,dum2;
      UInt_t fCw,fCh;

      Int_t fCanvasID=fCanvas->GetCanvasID();
      gVirtualX->GetGeometry(fCanvasID, dum1, dum2, fCw, fCh);

    printf("CreateCanvas: fCw, fCh %d %d \n", fCw, fCh );


   
   // Create main display pad
   fPad = new TPad("mviewpad", "ITS display",0.15,0,0.9,1);
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
   fTrigPad = new TPad("mtrigger", "range and mode pad",0,0,dxtr,dytr);
   fTrigPad->SetEditable(kFALSE);
   fTrigPad->Draw();
   fTrigPad->cd();
   fTrigPad->SetFillColor(22);
   fTrigPad->SetBorderSize(2);
   fRangeSlider = new TSlider("mrange","range",0.7,0.42,0.9,0.98);
   fRangeSlider->SetObject(this);

   // to use  "gAlice->Display()->SetPickMode()" KEEP "trigger" instead
   // of "mtrigger" and "range" instead of "mrange" !!!! - do NOT change
   // the pad names if mix methods from AliDisplay and AliITSdisplay !

   //char pickmode[] = "gAlice->Display()->SetPickMode()";
   char pickmode[] = "((AliITSdisplay*)(gAlice->Display()))->SetPickMode()";
   Float_t db = 0.09;
   fPickButton = new TButton("Pick",pickmode,0.05,0.32,0.65,0.32+db);
   fPickButton->SetFillColor(38);
   fPickButton->Draw();
   //char zoommode[] = "gAlice->Display()->SetZoomMode()";
   char zoommode[] = "((AliITSdisplay*)(gAlice->Display()))->SetZoomMode()";
   fZoomButton = new TButton("Zoom",zoommode,0.05,0.21,0.65,0.21+db);
   fZoomButton->SetFillColor(38);
   fZoomButton->Draw();
   fArcButton = new TArc(.8,fZoomButton->GetYlowNDC()+0.5*db,0.33*db);
   fArcButton->SetFillColor(kGreen);
   fArcButton->Draw();
   //char butUnzoom[] = "gAlice->Display()->UnZoom()";
   char butUnzoom[] = "((AliITSdisplay*)(gAlice->Display()))->UnZoom()";
   TButton *button = new TButton("UnZoom",butUnzoom,0.05,0.05,0.95,0.15);
   button->SetFillColor(38);
   button->Draw();
   AppendPad(); // append display object as last object to force selection

   fCanvas->cd();
   fCanvas->Update();
}


//_____________________________________________________________________________
AliITSdisplay::~AliITSdisplay()
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
void AliITSdisplay::Clear(Option_t *)
{
//    Delete graphics temporary objects
}

//_____________________________________________________________________________
void AliITSdisplay::DisplayButtons()
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
   char but1[] = "((AliITSdisplay*)(gAlice->Display()))->ShowNextEvent(1)";
   button = new TButton("Next",but1,x0,y-dbutton,x1,y);
   button->SetFillColor(38);
   button->Draw();

   y -= dbutton +dy;
   //char but2[] = "gAlice->Display()->ShowNextEvent(-1)";
   char but2[] = "((AliITSdisplay*)(gAlice->Display()))->ShowNextEvent(-1)";
   button = new TButton("Previous",but2,x0,y-dbutton,x1,y);
   button->SetFillColor(38);
   button->Draw();

   y -= dbutton + dy;
   char but3[] = "((AliITSdisplay*)(gAlice->Display()))->NextModule(1)";
   button = new TButton("Module +", but3, x0, y - dbutton, x1, y);
   button->SetFillColor(38);
   button->Draw();
    
   y -= dbutton + dy;
   char but4[] = "((AliITSdisplay*)(gAlice->Display()))->NextModule(-1)";
   button = new TButton("Module -", but4, x0, y - dbutton, x1, y);
   button->SetFillColor(38);
   button->Draw();

   y -= dbutton + dy;
   char but5[] = "((AliITSdisplay*)(gAlice->Display()))->DrawHistograms()";
   button = new TButton("DrawHist", but5, x0, y - dbutton, x1, y);
   button->SetFillColor(38);
   button->Draw();
   /*
   y -= dbutton +dy;
   char but6[] = "((AliITSdisplay*)(gAlice->Display()))->Trigger()";
   button = new TButton("Trigger",but6,x0,y-dbutton,x1,y);
   button->SetFillColor(38);
   button->Draw();
   */
   // display logo
   TDiamond *diamond = new TDiamond(0.05,0.015,0.95,0.22);
   diamond->SetFillColor(50);
   diamond->SetTextAlign(22);
   diamond->SetTextColor(5);
   diamond->SetTextSize(0.11);
   diamond->Draw();
   diamond->AddText(".. ");
   diamond->AddText("ROOT");
   diamond->AddText("ITS");
   diamond->AddText("... ");
   diamond->AddText(" ");
}

//_____________________________________________________________________________
void AliITSdisplay::CreateColors()
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
void AliITSdisplay::DisplayColorScale()
{

   Int_t i;
   Int_t color;
   Float_t xlow, ylow, xup, yup, hs;
   Float_t x1, y1, x2, y2;
   x1 = y1 = 0;
   x2 = y2 = 20;

   printf("DisplayColorScale - gPad %p\n",gPad);

   /*
   gPad->SetFillColor(0);
   gPad->Clear();
   gPad->Range(x1,y1,x2,y2);
   */

   fColPad->SetFillColor(0);
   fColPad->Clear();
   fColPad->Range(x1,y1,x2,y2);

   TText *text = new TText(0,0,"");
   text->SetTextFont(61);
   text->SetTextSize(0.2);
   text->SetTextAlign(22);

   /*
   AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
   AliMUONchamber *iChamber = &(MUON->Chamber(fChamber-1));
   AliMUONresponse * response=iChamber->GetResponseModel();
   Int_t adcmax= (Int_t) response->MaxAdc();
   */

   Int_t adcmax=1024;
   

   TBox *box;
   char label[8];
//*-* draw colortable boxes
   hs = (y2-y1)/Float_t(22);
   xlow=x1+1;
   xup=x2-9;
   for (i=0;i<22;i++) {
       ylow = y1 + hs*(Float_t(i));
       yup  = y1 + hs*(Float_t(i+1));
         color = 261+i;
	 Double_t logscale=Double_t(i+1)*(TMath::Log(adcmax)/22);
         Int_t scale=(Int_t)TMath::Exp(logscale);
	 sprintf(label,"%d",scale);
         box = new TBox(xlow, ylow, xup, yup);
         box->SetFillColor(color);
         box->Draw();
         text->DrawText(xup+4, 0.5*(ylow+yup),label);
   }
}

//______________________________________________________________________________
Int_t AliITSdisplay::DistancetoPrimitive(Int_t px, Int_t)
{
// Compute distance from point px,py to objects in event

   gPad->SetCursor(kCross);
   
   if (gPad == fTrigPad) return 9999;

   const Int_t big = 9999;
   Int_t dist   = big;
   Float_t xmin = gPad->GetX1();
   Float_t xmax = gPad->GetX2();
   Float_t dx   = 0.02*(xmax - xmin);
   Float_t x    = gPad->AbsPixeltoX(px);
   if (x < xmin+dx || x > xmax-dx) return dist;

   if (fZoomMode) return 0;
   else           return 7;
}

//_____________________________________________________________________________
void AliITSdisplay::Draw(Option_t *)
{
//    Display current event

   printf("Draw\n");

   fPad->cd();

   DrawView(fTheta, fPhi, fPsi);   
   // Display the event number and title
   fPad->cd();
   DrawTitle();
}


//_____________________________________________________________________________
void AliITSdisplay::DrawClusters()
{
//    Draw clusters for ITS modules

   if (!fDrawClusters) return;

   Int_t ndigits, digit;
   TObjArray *points;
   AliITSpoints *pm;

   LoadDigits(fModule);
      
   fClustersCuts = 0;
      points = Points();
      if (!points) return;
      ndigits = points->GetEntriesFast();
      for (digit=0;digit<ndigits;digit++){
         pm = (AliITSpoints*)points->UncheckedAt(digit);
         if (!pm) continue;
         Float_t *pxyz;
         pxyz=pm->GetP();
	 TMarker3DBox *marker=pm->GetMarker(0);
	 if (marker) marker->Draw();
      	 pm->Draw();
         fClustersCuts +=pm->GetN();

      }

    sModule=fModule;
}

//_____________________________________________________________________________
void AliITSdisplay::DrawHits()
{
//    Draw hits for ITS modules

   LoadHits(fModule);

   Int_t ntracks, track;
   TObjArray *points;
   AliITSpoints *pm;

   fHitsCuts = 0;
      points = Phits();
      if (!points) return;
      ntracks = points->GetEntriesFast();
      for (track=0;track<ntracks;track++) {
         pm = (AliITSpoints*)points->UncheckedAt(track);
	 if (!pm) continue;
         pm->Draw();
         fHitsCuts += pm->GetN();
      }
}


//_____________________________________________________________________________
void AliITSdisplay::DrawCoG()
{
//    Draw rec hits for ITS module

  /*
    if (!fDrawCoG) return;
   LoadCoG(fChamber,fCathode);

   Int_t ncog, icog;
   TObjArray *points;
   AliITSpoints *pm;

      points = Rpoints();
      if (!points) return;
      ncog = points->GetEntriesFast();
      for (icog=0;icog<ncog;icog++) {
         pm = (AliITSpoints*)points->UncheckedAt(icog);
	 if (!pm) continue;
         pm->Draw();
      }
  */
}
//_____________________________________________________________________________
void AliITSdisplay::DrawCoG2()
{
//    Draw rec hits for ITS module

  /*

   if (!fDrawCoG) return;
  

  if (fCathode==1) {
     LoadCoG2(fChamber,2);
  } else if (fCathode==2) {
     LoadCoG2(fChamber,1);
  }

   Int_t ncog, icog;
   TObjArray *points;
   AliITSpoints *pm;

      points = R2points();
      if (!points) return;
      ncog = points->GetEntriesFast();
      for (icog=0;icog<ncog;icog++) {
         pm = (AliITSpoints*)points->UncheckedAt(icog);
	 if (!pm) continue;
         pm->Draw();
      }
  */
}
//_____________________________________________________________________________
void AliITSdisplay::DrawCathCor()
{
//    Draw hits for ITS chambers

  /*
    
   if (!fDrawCathCor) return;

   LoadCathCor(fChamber);

   Int_t ncog, icog;
   TObjArray *points;
   AliITSpoints *pm;

      points = Cpoints();
      if (!points) return;
      ncog = points->GetEntriesFast();
      for (icog=0;icog<ncog;icog++) {
         pm = (AliITSpoints*)points->UncheckedAt(icog);
	 if (!pm) continue;
         pm->Draw();
      }
  */
}

//_____________________________________________________________________________
void AliITSdisplay::DrawTitle(Option_t *option)
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
      sprintf(ptitle,"Alice event: %d, Run:%d Module:%d",gAlice->GetHeader()->GetEvent(), gAlice->GetHeader()->GetRun(),fModule);
      title->AddText(ptitle);
      Int_t nparticles = gAlice->Particles()->GetEntriesFast();
      sprintf(ptitle,"Nparticles = %d Nhits = %d Npads fired = %d",nparticles, fHitsCuts,fClustersCuts);
      title->AddText(ptitle);
   } else {
      TPaveLabel *label = new TPaveLabel(xmin +0.01*dx, ymax-0.07*dy, xmin +0.2*dx, ymax-0.01*dy,option);
      label->SetBit(kCanDelete);
      label->SetFillColor(42);
      label->Draw();
   }
}

//_____________________________________________________________________________
void AliITSdisplay::DrawView(Float_t theta, Float_t phi, Float_t psi)
{
//    Draw a view of ITS clusters

  printf("DrawView - fPad gPad %p %p\n",fPad,gPad);

   gPad->SetCursor(kWatch);
   gPad->SetFillColor(1);
   gPad->Clear();
   gPad->SetFillColor(1);

   Int_t iret=0;
   TView *view = new TView(1);

   printf("DrawView - view %p \n",view);

   Float_t range = fRrange*fRangeSlider->GetMaximum();

   printf("DrawView - range fRrange fRangeSlider %f %f %f \n",range,fRrange,fRangeSlider->GetMaximum());
   view->SetRange(-range,-range,-range,range, range, range);
   fZoomX0[0] = -1;
   fZoomY0[0] = -1;
   fZoomX1[0] =  1;
   fZoomY1[0] =  1;
   fZooms = 0;

// Display ITS Chamber Geometry
// gAlice->GetGeometry()->Draw("same");
   char NodeName[7];
   sprintf(NodeName,"ITS%d",100+fModule);
   printf("Node name %s\n", NodeName);
   
   TNode *node1=gAlice->GetGeometry()->GetNode(NodeName);
   if (node1) node1->Draw("same");  

//add clusters to the pad
   DrawClusters();
   DrawHits();
   //DrawCoG();
   // DrawCoG2();
   //DrawCathCor();

  printf("DrawView - before append fPad gPad %p %p\n",fPad,gPad);
   // add itself to the list (must be last)
   AppendPad();
   view->SetView(phi, theta, psi, iret);
}


//______________________________________________________________________________
void AliITSdisplay::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
//  Execute action corresponding to the mouse event

   static Float_t x0, y0, x1, y1;

   static Int_t pxold, pyold;
   static Int_t px0, py0;
   static Int_t linedrawn;
   Float_t temp;


   /*
   //printf("ExecuteEvent - px py %d %d\n",px,py);

   int px1 = gPad->GetEventX();
   int py1 = gPad->GetEventY();
   float uxmin = gPad->GetUxmin();
   float uxmax = gPad->GetUxmax();
   int pxmin = gPad->XtoAbsPixel(uxmin);
   int pxmax = gPad->XtoAbsPixel(uxmax);

   //printf("ExecuteEvent - px1 py1 uxmin uxmax pxmin pxmax %d %d %f %f %d %d\n",px1,py1,uxmin,uxmax,pxmin,pxmax);

   Float_t x = gPad->AbsPixeltoX(px);
   Float_t y = gPad->AbsPixeltoY(py);
   //printf("x=%.3g, y=%.3g \n",gPad->PadtoX(x),gPad->PadtoY(y));

   TObject *select = gPad->GetSelected();
   if(!select) {printf("no select \n"); return;}
   if (select->InheritsFrom("AliITSdisplay")) printf("Inherits from AliITSdisplay \n");

   if (select->InheritsFrom("TCanvas")) printf("Inherits from TCanvas \n");
   if (select->InheritsFrom("TView")) printf("Inherits from TView \n");

    Float_t x,y;

   // x,y in the -1,1 scale
      x = gPad->AbsPixeltoX(px);
      y = gPad->AbsPixeltoY(py);
     
   printf("abspixeltoX: px py %d %d x=%.3g, y=%.3g \n",px,py,x,y);
   printf("x=%.3g, y=%.3g \n",gPad->PadtoX(x),gPad->PadtoY(y));

   // it should be smth like ? - no ! 
      x = gPad->AbsPixeltoX(px)*fRrange*fRangeSlider->GetMaximum();
      y = gPad->AbsPixeltoY(py)*fRrange*fRangeSlider->GetMaximum();
      // now find out how to convert these x,y to x,y in detector system

      Int_t dum1,dum2;
      UInt_t fCw, fCh;
      Float_t fYsizeReal, fXsizeReal;

      Int_t fCanvasID=fCanvas->GetCanvasID();
      gVirtualX->GetGeometry(fCanvasID, dum1, dum2, fCw, fCh);

    printf("Exec: fCw, fCh %d %d \n", fCw, fCh );

   if (fCw < fCh) {
      fYsizeReal = 20;
      fXsizeReal = fYsizeReal*Float_t(fCw)/Float_t(fCh);
   }
   else {
      fXsizeReal = 20;
      fYsizeReal = fXsizeReal*Float_t(fCh)/Float_t(fCw);
   }

    printf("Exec: fYsizeReal, fXsizeReal %f %f \n",fYsizeReal, fXsizeReal);

      fXsizeReal = fXsizeReal*Float_t(px)/Float_t(fCw);
      fYsizeReal = fYsizeReal*Float_t(py)/Float_t(fCh);

    printf("Exec: fYsizeReal, fXsizeReal %f %f \n",fYsizeReal, fXsizeReal);

   */

    if (!gPad) return;

    Float_t scale[3],center[3];
    Int_t irep;
    gPad->GetView()->FindScope(scale,center,irep);
    for(int i=0;i<3;i++) {
      //printf("Exec: scale center irep %f %f %d \n",scale[i],center[i],irep);
    }
    Float_t x,y;
   // x,y in the -1,1 scale
    x = gPad->AbsPixeltoX(px);
    y = gPad->AbsPixeltoY(py);
    Float_t xhit=x*scale[0]/(fRrange*fRangeSlider->GetMaximum());
    Float_t yhit=y*scale[0]/(fRrange*fRangeSlider->GetMaximum());
    //printf("Exec: x*scale y*scale xhit yhit%f %f %f %f  \n",x*scale[0],y*scale[1],xhit,yhit);
    
    Int_t anode, timebin;
    GetPadIxy(anode,timebin,x*scale[0],y*scale[0]);
    //printf("Exec: anode timebin %d %d  \n",anode,timebin);


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
void AliITSdisplay::LoadDigits(Int_t module)
{
// Read digits info and store x,y,z info in arrays fPoints

   printf("LoadDigits - fZooms %d\n",fZooms);

   if (module > fNmodules ) return;

   fModule=module;

   ResetPoints();

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   AliITSgeom *geom = ITS->GetITSgeom();

   // get detector type
   AliITSDetType *iDetType = ITS->DetType(fLayer); 
   Int_t id=(Int_t)((fLayer-1)/2);

   // for the moment do only SDD
   if (id != 1) return;

   TClonesArray *ITSdigits  = ITS->DigitsAddress(id);
   if (ITSdigits == 0) return;

   Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
   printf("Entries in TreeD - nent %d\n",nent);
   Int_t lastSPD=geom->GetLastSPD();
   Int_t lastSDD=geom->GetLastSDD();
   if (nent < lastSDD) {
      printf("You have not done the digitisation for all det types !!!");
      module=fModule-lastSPD;
   }


   ITS->ResetDigits();
   gAlice->TreeD()->GetEvent(module);
   //gAlice->TreeD()->GetEvent(nent-nofmodules+module-1);
   Int_t ndigits = ITSdigits->GetEntriesFast();
   printf("Found %d digits for module %d in det type %d \n",ndigits,fModule,id+1);
   if (ndigits == 0) return;
   if (fPoints == 0) fPoints = new TObjArray(ndigits);
   
   //AliITSresponseV0* resp = (AliITSresponseV0*)iDetType->GetResponseModel();
   AliITSresponse* resp = (AliITSresponse*)iDetType->GetResponseModel();
   resp->SetDetSize(38400,75200);
   resp->SetDetParam(376,256,40.,6.);

   if (!fMap) fMap = new AliITSMapA2(resp); 
   if (sModule != fModule) fMap->ClearMap();


   // if parameters in response where changed during digitisation set these
   // parameters here as well !! resp->SetDetParam(....), SetDetSize(...)
   // or pass in the constructor of display the pointer to a response 
   // instantiation done in ITSdisplay.C !

   AliITSdigitSDD  *mdig;
   AliITSpoints *points = 0;
   TMarker3DBox *marker=0;
   Float_t xl[3], dx[3];
   //
   //loop over all digits and store their position
   Int_t adcsatm=1024;
   Int_t npoints=1;
   for (Int_t digit=0;digit<ndigits;digit++) {
        mdig    = (AliITSdigitSDD*)ITSdigits->UncheckedAt(digit);
        //
        // First get all needed parameters
        //
        Int_t charge=0;
        Int_t chargecolor=0;
        Int_t signal=mdig->fSignal;
        for (int i=0;i<3;i++) {charge += mdig->fTcharges[i];}
	//printf("charge %d\n",charge);
        Int_t chargeindex=Int_t(TMath::Log(charge)/(TMath::Log(adcsatm)/22));
        Int_t index=Int_t(TMath::Log(signal)/(TMath::Log(adcsatm)/22));
        Int_t color=261+index;
        if (color>282) color=282;
        if (charge) chargecolor=261+chargeindex;
        if (chargecolor>282) chargecolor=282;
        Int_t anode=mdig->fCellX ;  
        Int_t timebin=mdig->fCellY; 
        fMap->SetHit(anode,timebin,(Double_t)signal);
	// get the center of the pad - add on x and y half of pad size
        GetPadCxy(anode,timebin,xl,dx);
      
        Float_t xpad=xl[0];
        Float_t ypad=xl[1];
        Float_t zpad=xl[2];

        Float_t dpx=dx[0];
        Float_t dpy=dx[1];
        Float_t dpz=dx[2];

        //printf("anode timebin charge signal %d %d %d %d \n",anode, timebin,charge,signal);
        //printf("xpad ypad zpad %f %f %f\n",xpad,ypad,zpad);
	//printf(" dpx dpz %f %f \n",dpx,dpz);
	//
	// Then set the objects
	//
        points = new AliITSpoints(npoints);
	fPoints->AddAt(points,digit);

        points->SetMarkerColor(color);
        points->SetMarkerStyle(21);
        //points->SetMarkerSize(0.5);
        points->SetMarkerSize(0.7);
        points->SetParticle(-1);
        points->SetDigitIndex(digit);
        points->SetPoint(0,xpad,ypad,zpad);	
	marker=new TMarker3DBox(xpad,ypad,zpad,dpx,dpy,dpz,0,0);
	//marker->SetLineColor(2);
	marker->SetLineColor(chargecolor);
	marker->SetFillStyle(1001);
	marker->SetFillColor(color);
	marker->SetRefObject((TObject*)points);
	points->Set3DMarker(0, marker);

   }
}
//___________________________________________
void AliITSdisplay::GetPadCxy(Int_t anode, Int_t timebin, Float_t *xl, Float_t *dx  )

{
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   // get detector type
   AliITSDetType *iDetType = ITS->DetType(fLayer); 
   AliITSresponse* resp = (AliITSresponse*)iDetType->GetResponseModel();


    Int_t na,ns;
    Float_t clock, speed;
    resp->GetDetParam(na,ns,clock,speed);
    //printf("resp na ns clock speed %p %d %d %f %f\n",resp,na,ns,clock,speed);

    Float_t timestep=1000./clock;
    Float_t pitch=resp->Pitch();
    Float_t width=resp->Width();
    Float_t length=resp->Length();
    //printf("timestep pitch width length %f %f %f %f \n",timestep,pitch,width,length);

    const Float_t kconv=10000;  // um->cm


    //Float_t driftpath=timebin*timestep*speed;
    Float_t driftpath=(timebin+1)*timestep*speed;
    if (anode >= na) xl[0]=(length-driftpath)/kconv;
    else xl[0]= -(length-driftpath)/kconv;
    if (anode >= na) anode-=na;
    //xl[2]=(anode*pitch-width/2)/kconv;
    xl[2]=((anode+1)*pitch-width/2)/kconv;
    xl[1]=fYmodule;
    //printf("anode timebin %d %d \n",anode,timebin);
    //printf("anode*pitch timebin*timestep*speed %f %f\n",anode*pitch,timebin*timestep*speed);

    // keep these for drawing
    //dx[0]=timestep*speed/kconv/2.1;
    //dx[1]=0;
    //dx[2]=pitch/kconv/2.05;

    dx[0]=timestep*speed/kconv/2;
    dx[1]=0;
    dx[2]=pitch/kconv/2;
}

//___________________________________________
void AliITSdisplay::GetPadIxy(Int_t &anode, Int_t &timebin,Float_t x,Float_t z)
{
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   // get detector type
   AliITSDetType *iDetType = ITS->DetType(fLayer); 
   AliITSresponse* resp = (AliITSresponse*)iDetType->GetResponseModel();


    Int_t na,ns;
    Float_t clock, speed;
    resp->GetDetParam(na,ns,clock,speed);

    Float_t timestep=1000./clock;
    Float_t pitch=resp->Pitch();
    Float_t width=resp->Width();
    Float_t length=resp->Length();

    const Float_t kconv=10000;  // cm->um

    Float_t driftpath=length-TMath::Abs(kconv*x);
    timebin=(Int_t)(driftpath/speed/timestep);
    anode=kconv*z/pitch + na/2;
    if (x > 0) anode += na;

    timebin+=1;
    anode+=1;
    //printf("anode timebin %d %d\n",anode,timebin);

}
//___________________________________________
void AliITSdisplay::LoadCoG(Int_t chamber, Int_t cathode)
{
// Read raw clusters info and store x,y,z info in arrays fRpoints
// Loop on all detectors

  /*

   if (chamber > 10) return;

   ResetRpoints();

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   AliITSchamber*  iChamber;

   TClonesArray *ITSrawclust  = ITS->RawClustAddress(chamber-1);
   if (ITSrawclust == 0) return;

   ITS->ResetRawClusters();


   Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
   gAlice->TreeR()->GetEvent(nent-2+cathode-1);
   //gAlice->TreeR()->GetEvent(cathode);
   Int_t nrawcl = ITSrawclust->GetEntriesFast();
   if (nrawcl == 0) return;
   if (fRpoints == 0) fRpoints = new TObjArray(nrawcl);
	  
   iChamber = &(ITS->Chamber(chamber-1));
   Float_t zpos=iChamber->ZPosition();  
   AliITSRawCluster  *mRaw;
   AliITSpoints *points = 0;
   //
   //loop over all raw clusters and store their position
   points = new AliITSpoints(nrawcl);
   for (Int_t iraw=0;iraw<nrawcl;iraw++) {
  	mRaw   = (AliITSRawCluster*)ITSrawclust->UncheckedAt(iraw);
	fRpoints->AddAt(points,iraw);
        points->SetMarkerColor(51);
        points->SetMarkerStyle(2);
        points->SetMarkerSize(1.);
        points->SetParticle(-1);
        points->SetHitIndex(-1);
        points->SetTrackIndex(-1);
        points->SetDigitIndex(-1);
        points->SetPoint(iraw,mRaw->fX,mRaw->fY,zpos);
   }
  */
}
//___________________________________________
void AliITSdisplay::LoadCoG2(Int_t chamber, Int_t cathode)
{
// Read raw clusters info and store x,y,z info in arrays fRpoints
// Loop on all detectors
  /*

   if (chamber > 10) return;

   ResetR2points();

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   AliMUONchamber*  iChamber;

   TClonesArray *MUONrawclust  = MUON->RawClustAddress(chamber-1);
   if (MUONrawclust == 0) return;

   MUON->ResetRawClusters();

   Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
   gAlice->TreeR()->GetEvent(nent-2+cathode-1);
   //gAlice->TreeR()->GetEvent(cathode);
   Int_t nrawcl = MUONrawclust->GetEntriesFast();
   if (nrawcl == 0) return;
   if (fR2points == 0) fR2points = new TObjArray(nrawcl);
	  
   iChamber = &(MUON->Chamber(chamber-1));
   Float_t zpos=iChamber->ZPosition();  
   AliMUONRawCluster  *mRaw;
   AliMUONpoints *points = 0;
   //
   //loop over all raw clusters and store their position
   points = new AliMUONpoints(nrawcl);
   for (Int_t iraw=0;iraw<nrawcl;iraw++) {
  	mRaw   = (AliMUONRawCluster*)MUONrawclust->UncheckedAt(iraw);
	fR2points->AddAt(points,iraw);
        points->SetMarkerColor(51);
        points->SetMarkerStyle(4);
        points->SetMarkerSize(1.3);
        points->SetParticle(-1);
        points->SetHitIndex(-1);
        points->SetTrackIndex(-1);
        points->SetDigitIndex(-1);
        points->SetPoint(iraw,mRaw->fX,mRaw->fY,zpos);
   }
  */
}
//___________________________________________
void AliITSdisplay::LoadCathCor(Int_t chamber)
{
// Read correlation info and store x,y,z info in arrays fCpoints
// Loop on all detectors

  /*
     if (chamber > 10) return;
     fChamber=chamber;
     ResetCpoints();

     AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
     AliMUONchamber*  iChamber;
     iChamber = &(MUON->Chamber(chamber-1));
     Float_t zpos=iChamber->ZPosition();  // check with Andreas


   //new
     MUON->GetTreeC(fEvent);
     TTree *TC=MUON->TreeC();
     if (!TC) return;
     //     Int_t nent=(Int_t)TC->GetEntries();
 
     TClonesArray *MUONcorrel  = MUON->CathCorrelAddress(chamber-1);
     if (MUONcorrel == 0) return;

     MUON->ResetCorrelation();
     TC->GetEvent();

     Int_t ncor = MUONcorrel->GetEntries();
     if (!ncor) return;
     if (!fCpoints) fCpoints = new TObjArray(ncor);
	  
     AliMUONcorrelation  *mCor;
     AliMUONpoints *points = 0;
   //
   //loop over all raw clusters and store their position
     points = new AliMUONpoints(ncor);
     for (Int_t icor=0;icor<ncor;icor++) {
          mCor   = (AliMUONcorrelation*)MUONcorrel->UncheckedAt(icor);
	  fCpoints->AddAt(points,icor);
          points->SetMarkerColor(4);
          points->SetMarkerStyle(4);
          points->SetMarkerSize(0.8);
          points->SetParticle(-1);
          points->SetHitIndex(-1);
          points->SetTrackIndex(-1);
          points->SetDigitIndex(-1);
          points->SetPoint(icor,mCor->fX[0],mCor->fY[0],zpos);
    }
  */
}

//___________________________________________
void AliITSdisplay::LoadHits(Int_t module)
{

// Read hits info and store x,y,z info in arrays fPhits
// Loop on all detectors


   if (module > fNmodules || fNmodules==0) return;

   fModule=module;
   printf("LoadModuleHits - fModule module %d %d\n",fModule,module);
 
   ResetPhits();

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   AliITSgeom *geom = ITS->GetITSgeom();

   // try fITSmodules first (memory) - check also TreeH with the condition
   // if (module != fModule) continue; - see the speed

    TObjArray *mods = ITS->GetModules();
    Int_t nentries=(Int_t)mods->GetEntriesFast();
    printf("nentries %d\n",nentries);

    AliITSmodule *mod = (AliITSmodule *)mods->UncheckedAt(fModule);
    TObjArray *fHits = mod->GetHits();
    Int_t nhits = fHits->GetEntriesFast();
    if (fPhits == 0) fPhits = new TObjArray(nhits);
    printf("LoadModuleHits - fPhits fModule nhits %p %d %d\n",fPhits,fModule,nhits);
    AliITSpoints *points = 0;
    AliITShit *mHit = 0;
    Int_t idtrack, idhit;
    for (Int_t hit=0;hit<nhits;hit++) {
        mHit = (AliITShit*)fHits->UncheckedAt(hit);
        mod->GetHitTrackAndHitIndex(hit,idtrack,idhit);
        printf("hit eloss %d %f\n",hit,mHit->GetIonization());
        printf("layer ladder det %d %d %d\n",mHit->GetLayer(),mHit->GetLadder(),mHit->GetDetector());
        if (!mHit->GetIonization()) continue;
        //fYmodule=mHit->fY;
        Float_t xl[3], xg[3];
        xg[0]=mHit->fX;
        xg[1]=mHit->fY;
        xg[2]=mHit->fZ;

	geom->GtoL(fModule,xg,xl);
        printf("hit xhit yhit zhit %d %f %f %f\n",hit,xl[0],xl[1],xl[2]);
        Float_t x,y,z,tof;
        mHit->GetPositionL(x,y,z,tof);
        printf("hit x y z %d %f %f %f\n",hit,x,y,z);

        fYmodule=xl[1];
        //
	// Retrieve info and set the objects
	//
	points = new AliITSpoints();
	points->SetMarkerColor(kRed);
	points->SetMarkerStyle(5);
	points->SetMarkerSize(1.);
	points->SetParticle(mHit->fTrack);
	points->SetHitIndex(idhit);
	points->SetTrackIndex(idtrack);
	//points->SetPoint(0,mHit->fX,mHit->fY,mHit->fZ);
	points->SetPoint(0,xl[0],xl[1],xl[2]);
	fPhits->AddAt(points,hit);
    }

}

//_____________________________________________________________________________
void AliITSdisplay::NextModule(Int_t delta)
{
  // to go from module to next module if delta = 1
  // or previous module otherwise
    if (delta == 1) {
	if (fModule < fNmodules) fModule++;
    } else {
	if (fModule > 1) fModule--;
    }
    if (!fPad) return;
    fPad->Clear();
    Draw();
}

//_____________________________________________________________________________
void AliITSdisplay::Paint(Option_t *)
{
//    Paint miscellaneous items

}

//_____________________________________________________________________________
void AliITSdisplay::DrawHistograms()
{

  // better write a class for this - see TInspectCanvas - with a methode
  // called xxxx and here just do smth like it's done for Inspect() in TObject
  //gROOT->ProcessLine(Form("TInspectCanvas::Inspector((TObject *)0x%lx);",(Long_t)this));
  // with xxxx instead of Inspector - so one can create buttons on the 
  // new canvas and have more choices


   printf("DrawHistograms: fZooms fZoomMode  %d %d ! \n",fZooms,fZoomMode);

   if (!fZooms) return;

    Float_t scale[3],center[3];
    Int_t irep;
    gPad->GetView()->FindScope(scale,center,irep);

    Int_t amin, tmin, amax, tmax;
    GetPadIxy(amin,tmin, fZoomX0[fZooms]*scale[0], fZoomY0[fZooms]*scale[0]);
    printf("DrawHist: amin tmin %d %d  \n",amin,tmin);
    GetPadIxy(amax,tmax, fZoomX1[fZooms]*scale[0], fZoomY1[fZooms]*scale[0]);
    printf("DrawHist: amax tmax %d %d  \n",amax,tmax);

    if ((fZoomX0[fZooms] < 0 && fZoomX1[fZooms] > 0) || (fZoomX0[fZooms] > 0 && fZoomX1[fZooms] < 0) ) tmax=256;

    if( amax-amin > kMAXHIST) {
        Error("DrawHistograms","Too many histograms %d! Zoom a smaller area and try again!",amax-amin);
        return;
    }
    Int_t nbinx = amax-amin; 
    Int_t nbiny = tmax-tmin;

   fDrawHist = kTRUE;
   TObjArray *fHis=0;

    //create or set the new canvas c2
   TVirtualPad *padsav = gPad;
   TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
   /*
   if(c2) {
     if (fHis) {
       fHis->Delete(); 
       delete fHis;     
     }                
   }
   */
   if(c2) delete c2->GetPrimitive("h2map");
   else     c2=new TCanvas("c2","c2",14,47,830,750);
   c2->cd();

   TPad *hpad = new TPad("hpad","",0.01,0.01,0.99,0.99);
   hpad->SetFillColor(0);
   //hpad->Draw();
   hpad->Divide(2,nbinx/2,0.,0.);
   hpad->Draw();

   char title[30], name[10];
   sprintf(title,"Signal map for module %d",fModule);
    TH2F *fTH2F = new TH2F("h2map",title,nbinx,amin,amax,nbiny,tmin,tmax);
    fHis=new TObjArray(nbinx);

   //TH1F *h1=new TH1F(" ","",nbiny,(float)tmin,(float)tmax);

   for (int i=amin;i<amax;i++) {
     //h1->Reset();
      printf("i, i-amin  %d %d \n",i,i-amin);
      (*fHis)[i-amin] = new TH1F(" "," ",nbiny,(float)tmin,(float)tmax);
      //sprintf(title,"Projection on anode=%d",i);
      sprintf(title,"anode=%d",i);
      sprintf(name,"h%d",i);
      ((TH1F*)((*fHis)[i-amin]))->SetTitle(title);
      ((TH1F*)((*fHis)[i-amin]))->SetTitleSize(2.);
      ((TH1F*)((*fHis)[i-amin]))->SetTitleOffset(2.);
      ((TH1F*)((*fHis)[i-amin]))->SetLabelSize(0.05);
      //((TH1F*)((*fHis)[i-amin]))->SetYTitle(title);
      ((TH1F*)((*fHis)[i-amin]))->SetName(name);
      ((TH1F*)((*fHis)[i-amin]))->SetStats(0);
      //h1->SetTitle(title);
      //h1->SetTitleSize(0.04);
      //h1->SetYTitle(title);
      //h1->SetName(title);
      hpad->cd(i-amin+1);
      for (int j=tmin;j<tmax;j++) {
	((TH1F*)((*fHis)[i-amin]))->Fill((float)j,(float)fMap->GetSignal(i,j));
	//h1->Fill((float)j,(float)fMap->GetSignal(i,j));
          fTH2F->Fill((float)i,(float)j,(float)fMap->GetSignal(i,j));
      }
      ((TH1F*)((*fHis)[i-amin]))->Smooth(1);
      //((TH1F*)((*fHis)[i-amin]))->Fit("gaus","ql");
      ((TH1F*)((*fHis)[i-amin]))->Draw();
      //h1->Fit("gaus","ql");
      //h1->Draw();
      //hpad->Update();
   }
   c2->Update();
   //c2->SetFillColor(1);
   //fTH2F->Draw("LEGO2");

   padsav->cd();
   
   fDrawHist=kFALSE;

   /*
  if (fHis) {
     fHis->Delete(); 
     delete fHis;     
  }                
   
   */


   /*
   //create or set the new canvas c2
   TVirtualPad *padsav = gPad;
   TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
   if(c2) delete c2->GetPrimitive("Projection");
   if(c2) delete c2->GetPrimitive("h2map");
   else     c2=new TCanvas("c2","c2",14,47,830,750);
   c2->cd();

   c2->SetFillColor(1);
   fTH2F->Draw("LEGO2");

   TPad *hpad = new TPad("hpad","",0.01,0.01,0.99,0.99);
   hpad->SetFillColor(11);
   hpad->Divide(1,nbinx/2,0.,0.);
   hpad->Draw();
   // -empty projections !!

   //draw slice corresponding to mouse position
   for(int i=amin;i<amax;i++) {
     hpad->cd(i-amin+1);
     //Int_t biny = fTH2F->GetYaxis()->FindBin(y);
     //TH1D *hp = fTH2F->ProjectionX("",biny,biny);
     TH1D *hp = fTH2F->ProjectionY("",i,i);
     char title[80];
     sprintf(title,"Projection on anode=%d",i);
     hp->SetName("Projection");
     hp->SetTitle(title);
     //hp->Fit("gaus","ql");
     hp->Draw();
     hpad->Update();
   }
     c2->Update();
     padsav->cd();
   
   fDrawHist=kFALSE;
   */
}

//_____________________________________________________________________________
void AliITSdisplay::SetPickMode()
{
   fZoomMode = 0;

   fArcButton->SetY1(fPickButton->GetYlowNDC()+0.5*fPickButton->GetHNDC());
   fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliITSdisplay::SetZoomMode()
{
   fZoomMode = 1;

   fArcButton->SetY1(fZoomButton->GetYlowNDC()+0.5*fZoomButton->GetHNDC());
   fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliITSdisplay::SetModule(Int_t layer, Int_t ladder, Int_t detector)
{
// Set chamber and cathode number
   fLayer = layer;
   fLadder = ladder;
   fDetector = ladder;

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   AliITSgeom *gm=ITS->GetITSgeom();
   fModule=gm->GetModuleIndex(fLayer,fLadder,fDetector);

   if (!fPad) return;
   fPad->Clear();
   Draw();
}

//_____________________________________________________________________________
void AliITSdisplay::SetModuleNumber(Int_t module)
{
// Set chamber and cathode number
   fModule=module;

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   AliITSgeom *gm=ITS->GetITSgeom();
   gm->GetModuleId(fModule,fLayer,fLadder,fDetector);   

   if (!fPad) return;
   fPad->Clear();
   Draw();
}
//_____________________________________________________________________________
void AliITSdisplay::SetRange(Float_t rrange, Float_t zrange)
{
// Set view range along R and Z
   fRrange = rrange;
   fZrange = zrange;

   if (!fPad) return;
   fPad->Clear();
   Draw();
}
   
//_____________________________________________________________________________
void AliITSdisplay::SetView(Float_t theta, Float_t phi, Float_t psi)
{
//  change viewing angles for current event

   fPad->cd();
   fPhi   = phi;
   fTheta = theta;
   fPsi   = psi;
   Int_t iret = 0;

   TView *view = gPad->GetView();
   printf("SetView view %p \n",view);
   if (view) view->SetView(fPhi, fTheta, fPsi, iret);
   Draw();
   //else      Draw();

   gPad->Modified();
}

//_____________________________________________________________________________
void AliITSdisplay::ShowNextEvent(Int_t delta)
{
//  Display (current event_number+delta)
//    delta =  1  shown next event
//    delta = -1 show previous event
  if (delta) {
     gAlice->Clear();
     Int_t current_event = gAlice->GetHeader()->GetEvent();
     Int_t new_event     = current_event + delta;
     gAlice->GetEvent(new_event);
     fEvent=new_event;
     if (!gAlice->TreeD() || !gAlice->TreeD() ) return; 
   }
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  ITS->ClearModules();
  ITS->InitModules(-1,fNmodules); 
  ITS->FillModules(fEvent,0,fNmodules," "," "); //memory problems ?

  fPad->cd(); 
  Draw();
}

//_____________________________________________________________________________

void AliITSdisplay::SetEvent(Int_t newevent)
{
    gAlice->GetEvent(newevent);
    fEvent=newevent;
    if (!gAlice->TreeD() || !gAlice->TreeD() ) return; 
    if (!fPad) return;
    AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
    ITS->ClearModules();
    ITS->InitModules(-1,fNmodules); 
    ITS->FillModules(fEvent,0,fNmodules," "," "); //memory problems ?
    fPad->cd(); 
    fPad->Clear();
    Draw();

}
//______________________________________________________________________________
void AliITSdisplay::UnZoom()
{
   printf("UnZoom - fZooms %d\n",fZooms);
   if (fZooms <= 0) return;
   fZooms--;
   TPad *pad = (TPad*)gPad->GetPadSave();
   printf("UnZoom - pad fPad gPad %p %p %p\n",pad,fPad,gPad);

   //pad->cd();

   pad->Range(fZoomX0[fZooms],fZoomY0[fZooms], fZoomX1[fZooms],fZoomY1[fZooms]);
   pad->Modified();
}

//_____________________________________________________________________________
void AliITSdisplay::ResetPoints()
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
void AliITSdisplay::ResetPhits()
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
void AliITSdisplay::ResetRpoints()
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
void AliITSdisplay::ResetR2points()
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
void AliITSdisplay::ResetCpoints()
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











