
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
#include <TPolyMarker.h>
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
#include <TGXW.h>
#include <TMath.h>
#include <X3DBuffer.h>

#include "AliRun.h"
#include "AliDetector.h"
#include "AliMUON.h"
#include "AliMUONConst.h"
#include "AliMUONdisplay.h"
#include "AliMUONpoints.h"
#include "GParticle.h"


ClassImp(AliMUONdisplay)


//_____________________________________________________________________________
AliMUONdisplay::AliMUONdisplay()
{
   fPoints = 0;
   fPhits = 0;
   fCanvas = 0;
}

//_____________________________________________________________________________
AliMUONdisplay::AliMUONdisplay(Int_t size)
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
<img src="picts/aliMUONdisplay.gif">
*/
//End_Html


   fPad = 0;

   gAlice->SetDisplay(this);
   
   // Initialize display default parameters
   SetRange();

   // Set front view by default
   fTheta = 0;
   fPhi   = -90;
   fPsi   = 0;
   fChamber = 1;
   fCathode = 1;
   //   fRzone   = 1.e10;
   fDrawClusters  = kTRUE;
   fZoomMode      = 1;
   fZooms         = 0;
   fClustersCuts  = 0;
   fPoints        = 0;
   fPhits        = 0;
   // Create colors
   CreateColors();
   // Create display canvas
   Int_t ysize = size;
   if (ysize < 100) ysize = 750;
   Int_t xsize = Int_t(size*830./ysize);
   fCanvas = new TCanvas("Canvas", "MUON Clusters Display",14,47,xsize,ysize);
   fCanvas->SetEditable(kIsNotEditable);
   fCanvas->ToggleEventStatus();
   
   // Create main display pad
   fPad = new TPad("viewpad", "MUON display",0.15,0,0.9,1);
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
AliMUONdisplay::~AliMUONdisplay()
{
  // Delete space point structure
  if (fPoints) fPoints->Delete();
  delete fPoints;
  fPoints     = 0;
  //
  if (fPhits) fPhits->Delete();
  delete fPhits;
  fPhits     = 0;
}

//_____________________________________________________________________________
void AliMUONdisplay::Clear(Option_t *)
{
//    Delete graphics temporary objects
}

//_____________________________________________________________________________
void AliMUONdisplay::DisplayButtons()
{
//    Create the user interface buttons


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
   button = new TButton("Next",but1,x0,y-dbutton,x1,y);
   button->SetFillColor(38);
   button->Draw();

   y -= dbutton +dy;
   char but2[] = "gAlice->Display()->ShowNextEvent(-1)";
   button = new TButton("Previous",but2,x0,y-dbutton,x1,y);
   button->SetFillColor(38);
   button->Draw();
   /*
   y -= dbutton +dy;
   char but3[] = "gAlice->Display()->SetChamberAndCathode(1,1)";
   button = new TButton("Cham&Cath",but3,x0,y-dbutton,x1,y);
   button->SetFillColor(butcolor);
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
   diamond->AddText("MUON");
   diamond->AddText("... ");
   diamond->AddText(" ");
}

//_____________________________________________________________________________
void AliMUONdisplay::CreateColors()
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
        color=50+23-color;
	//        printf("CreateColors - i, k, color %d %d %d \n",i,k,color);
        new TColor(color,r,g,b);
      } 
      //      printf("I'm out of case %d \n",k);
      break;
    case 2:
      for (i=1;i<=4;i++) {
        r=1.1-i*0.2;
        g=1.;  
        b=0.;
        color=i+5;
        color=50+23-color;
	//        printf("CreateColors - i, k, color %d %d %d \n",i,k,color);
        new TColor(color,r,g,b);
      } 
      //      printf("I'm out of case %d \n",k);
      break;
    case 3:
      for (i=1;i<=4;i++) {
        r=0.;
        g=1.;  
        b=i*0.2+0.2;
        color=i+9;
        color=50+23-color;
	//        printf("CreateColors - i, k, color %d %d %d \n",i,k,color);
        new TColor(color,r,g,b);
      } 
      //      printf("I'm out of case %d \n",k);
      break;
    case 4:
      for (i=1;i<=4;i++) {
        r=0.;
        g=1.1-i*0.2;  
        b=1.;
        color=i+13;
        color=50+23-color;
	//        printf("CreateColors - i, k, color %d %d %d \n",i,k,color);
        new TColor(color,r,g,b);
      } 
      //      printf("I'm out of case %d \n",k);
      break;
    case 5:
      for (i=1;i<=5;i++) {
        r=i*0.2;
        g=0.;  
        b=1.;
        color=i+17;
        color=50+23-color;
	//        printf("CreateColors - i, k, color %d %d %d \n",i,k,color);
        new TColor(color,r,g,b);
      } 
      //      printf("I'm out of case %d \n",k);
      break;
    }
    
  }

}

//_____________________________________________________________________________
void AliMUONdisplay::DisplayColorScale()
{

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
         color = 51+i;
	 //	 Int_t scale=(i+1)*(Int_t)adc_satm/22;
	 //	 sprintf(label,"%d",scale);
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
Int_t AliMUONdisplay::DistancetoPrimitive(Int_t px, Int_t)
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
void AliMUONdisplay::Draw(Option_t *)
{
//    Display current event

   fPad->cd();

   DrawView(fTheta, fPhi, fPsi);      // see how to draw PGON+inner frames

   // Display the event number and title
   fPad->cd();
   DrawTitle();
}


//_____________________________________________________________________________
void AliMUONdisplay::DrawClusters()
{
//    Draw clusterss for MUON chambers

   Int_t ndigits, digit;
   TObjArray *points;
   AliMUONpoints *pm;

   //   rmax=fRzone;
      
   fClustersCuts = 0;
      points = Points();
      if (!points) return;
      ndigits = points->GetEntriesFast();
      printf("DrawClusters - ndigits %d \n",ndigits);
      for (digit=0;digit<ndigits;digit++){
         pm = (AliMUONpoints*)points->UncheckedAt(digit);
         if (!pm) continue;
	 //        Float_t *pxyz;
	 //        pxyz=pm->GetP();
	 //        printf("DrawClusters - pxyz[0],pxyz[1] %f %f \n",pxyz[0],pxyz[1]);
         pm->Draw();
	 //	          Int_t n=pm->GetN();
	 //	          printf("DrawClusters - n %d \n",n);
         fClustersCuts +=pm->GetN();
      }
}

//_____________________________________________________________________________
void AliMUONdisplay::DrawHits()
{
//    Draw hits for MUON chambers

   LoadHits(fChamber);

   Int_t ntracks, track;
   TObjArray *points;
   AliMUONpoints *pm;

   fHitsCuts = 0;
      points = Phits();
      if (!points) return;
      ntracks = points->GetEntriesFast();
      printf("DrawHits - ntracks %d \n",ntracks);
      for (track=0;track<ntracks;track++) {
         pm = (AliMUONpoints*)points->UncheckedAt(track);
	 if (!pm) continue;
         pm->Draw();
         fHitsCuts += pm->GetN();
      }
}


//_____________________________________________________________________________
void AliMUONdisplay::DrawTitle(Option_t *option)
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
      sprintf(ptitle,"Alice event: %d, Run:%d",gAlice->GetHeader()->GetEvent(), gAlice->GetHeader()->GetRun());
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
void AliMUONdisplay::DrawView(Float_t theta, Float_t phi, Float_t psi)
{
//    Draw a view of MUON clusters

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
   
   // Display MUON Chamber Geometry
   //   gAlice->GetGeometry()->Draw("same");
   
   //add clusters to the pad
   DrawClusters();
   DrawHits();

    // add itself to the list (must be last)
   AppendPad();
   
   view->SetView(phi, theta, psi, iret);
}


//______________________________________________________________________________
void AliMUONdisplay::ExecuteEvent(Int_t event, Int_t px, Int_t py)
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
      gGXW->SetLineColor(-1);
      gPad->TAttLine::Modify();  //Change line attributes only if necessary
      x0 = gPad->AbsPixeltoX(px);
      y0 = gPad->AbsPixeltoY(py);
      px0   = px; py0   = py;
      pxold = px; pyold = py;
      linedrawn = 0;
      return;

   case kButton1Motion:
      if (linedrawn) gGXW->DrawBox(px0, py0, pxold, pyold, TGXW::kHollow);
      pxold = px;
      pyold = py;
      linedrawn = 1;
      gGXW->DrawBox(px0, py0, pxold, pyold, TGXW::kHollow);
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
void AliMUONdisplay::LoadDigits(Int_t chamber, Int_t cathode)
{
// Read digits info and store x,y,z info in arrays fPoints
// Loop on all detectors

   fChamber=chamber;
   fCathode=cathode;
 
   ResetPoints();

   AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
   AliMUONchamber*  iChamber;
   AliMUONsegmentation*  segmentation;

   TClonesArray *MUONdigits  = MUON->DigitsAddress(chamber-1);
   if (MUONdigits == 0) return;

   gAlice->ResetDigits();
   gAlice->TreeD()->GetEvent(cathode);
   Int_t ndigits = MUONdigits->GetEntriesFast();
   if (ndigits == 0) return;
   if (fPoints == 0) fPoints = new TObjArray(ndigits);
   printf("Found %d digits for cathode %d in chamber %d \n",ndigits,cathode,chamber);
	  
   iChamber = &(MUON->Chamber(chamber-1));
   printf("LoadPoints - iChamber %p \n",iChamber);
   segmentation=iChamber->GetSegmentationModel(cathode);
   printf("LoadPoints - segmentation %p \n",segmentation);
   Float_t dpx  = segmentation->Dpx();
   Float_t dpy  = segmentation->Dpy();
   printf("LoadPoints - dpx, dpy %f %f \n",dpx,dpy);
   Float_t zpos=iChamber->ZPosition();
   printf("LoadPoint - zpos %f \n",zpos);
   AliMUONdigit  *mdig;
   AliMUONpoints *points = 0;
   //
   //loop over all digits and store their position
   //    points = new AliMUONpoints(ndigits);
   Int_t npoints=1;
   for (Int_t digit=0;digit<ndigits;digit++) {
       mdig    = (AliMUONdigit*)MUONdigits->UncheckedAt(digit);
        points = new AliMUONpoints(npoints);
	fPoints->AddAt(points,digit);
        Int_t charge=mdig->fSignal;
        // set the color according to the color scale
	//        Int_t scale=(Int_t)adc_satm/22;
	//        Int_t index=(Int_t)charge/scale;
        Int_t index=Int_t(TMath::Log(charge)/(TMath::Log(adc_satm)/22));
        Int_t color=51+index;
        if (color>72) color=72;
        points->SetMarkerColor(color);
        points->SetMarkerStyle(21);
        points->SetMarkerSize(0.5);
	// get the center of the pad - add on x and y half of pad size
	Float_t xpad, ypad;
	segmentation->GetPadCxy(mdig->fPadX, mdig->fPadY,xpad, ypad);
        points->SetParticle(-1);
        points->SetHitIndex(-1);
        points->SetTrackIndex(-1);
        points->SetDigitIndex(digit);
        points->SetPoint(0,xpad,ypad,zpos);
	//        Float_t *pxyz;
	//       pxyz=points->GetP();
	//        printf("pxyz[0],pxyz[1],color %f %f %d \n",pxyz[0],pxyz[1],color);
	//        Int_t np=points->GetN();
	//        printf("np %d \n",np);
   }
}


//___________________________________________
void AliMUONdisplay::LoadHits(Int_t chamber)
{
// Read hits info and store x,y,z info in arrays fPhits
// Loop on all detectors

   fChamber=chamber;
 
   ResetPhits();

   AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
   AliMUONchamber*  iChamber;

   iChamber = &(MUON->Chamber(chamber-1));
   Float_t zpos=iChamber->ZPosition();
   printf("LoadHits - zpos %f \n",zpos);

   Int_t ntracks = (Int_t)gAlice->TreeH()->GetEntries();
   printf("ntracks %d\n",ntracks);
   Int_t ntrks = gAlice->GetNtrack();
   printf("ntrks %d\n",ntrks);

   if (fPhits == 0) fPhits = new TObjArray(ntracks);
   //   if (fPhits == 0) fPhits = new TObjArray(ntrks);

    TVector *xp = new TVector(10);
    TVector *yp = new TVector(10);
    TVector *ptrk = new TVector(10);
    TVector *phit = new TVector(10);
    for (Int_t track=0; track<ntracks;track++) {
      gAlice->ResetHits();
      gAlice->TreeH()->GetEvent(track);
         TClonesArray *MUONhits  = MUON->Hits();
	 //	 printf("MUONhits %p\n",MUONhits);
         if (MUONhits == 0) return;
         Int_t nhits = MUONhits->GetEntriesFast();
         if (nhits == 0) continue;
	 //	 printf("nhits %d \n",nhits);
         AliMUONhit *mHit;
         AliMUONpoints *points = 0;
	 //         Int_t trko=-99, trk;
	 //         points = new AliPoints(nhits);
         Int_t npoints=0;
         for (Int_t hit=0;hit<nhits;hit++) {
            mHit = (AliMUONhit*)MUONhits->UncheckedAt(hit);
            Int_t nch  = mHit->fChamber;              // chamber number
            if (nch != chamber) continue;
            (*xp)(npoints)=mHit->fX;
            (*yp)(npoints)=mHit->fY;
            (*ptrk)(npoints)=Float_t(mHit->GetTrack());
            (*phit)(npoints)=Float_t(hit);
	    //	    printf("track, trk %d %d\n",track,mHit->GetTrack());
            npoints++;
	 }
         if (npoints == 0) continue;
	 //	 printf("npoints %d \n",npoints);
	 points = new AliMUONpoints(npoints);
         for (Int_t p=0;p<npoints;p++) {
            points->SetMarkerColor(kRed);
            points->SetMarkerStyle(5);
            points->SetMarkerSize(1.);
            points->SetParticle(Int_t((*ptrk)(p)));
	    //            Int_t index=points->GetIndex();
	    //	    printf("index %d \n",index);
            points->SetHitIndex(Int_t((*phit)(p)));
            points->SetTrackIndex(track);
            points->SetDigitIndex(-1);
            points->SetPoint(p,(*xp)(p),(*yp)(p),zpos);
	 }
      xp->Zero();
      yp->Zero();
      ptrk->Zero();
      phit->Zero();
      fPhits->AddAt(points,track);
      //            Int_t np=points->GetN();
      //            printf("np %d \n",np);
   }

}

//_____________________________________________________________________________
void AliMUONdisplay::Paint(Option_t *)
{
//    Paint miscellaneous items

}

//_____________________________________________________________________________
void AliMUONdisplay::SetPickMode()
{
   fZoomMode = 0;

   fArcButton->SetY1(fPickButton->GetYlowNDC()+0.5*fPickButton->GetHNDC());
   fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliMUONdisplay::SetZoomMode()
{
   fZoomMode = 1;

   fArcButton->SetY1(fZoomButton->GetYlowNDC()+0.5*fZoomButton->GetHNDC());
   fTrigPad->Modified();
}

//_____________________________________________________________________________
void AliMUONdisplay::SetChamberAndCathode(Int_t chamber, Int_t cathode)
{
// Set chamber and cathode number
   fChamber = chamber;
   fCathode = cathode;

   //   Pad();
   //   printf("SetChamberAndCathode - fPad %p \n",fPad);
   printf("SetChamberAndCathode - fChamber fCathode %d %d\n",fChamber,fCathode);
   if (!fPad) return;
   fPad->Clear();
   LoadDigits(chamber,cathode);
   Draw();
}

//_____________________________________________________________________________
void AliMUONdisplay::SetRange(Float_t rrange, Float_t zrange)
{
// Set view range along R and Z
   fRrange = rrange;
   fZrange = zrange;

   if (!fPad) return;
   fPad->Clear();
   Draw();
}
   
//_____________________________________________________________________________
void AliMUONdisplay::SetView(Float_t theta, Float_t phi, Float_t psi)
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
void AliMUONdisplay::ShowNextEvent(Int_t delta)
{
//  Display (current event_number+delta)
//    delta =  1  shown next event
//    delta = -1 show previous event

  if (delta) {
     gAlice->Clear();
     Int_t current_event = gAlice->GetHeader()->GetEvent();
     Int_t new_event     = current_event + delta;
     gAlice->GetEvent(new_event);
     if (!gAlice->TreeD()) return; 
   }
  LoadDigits(fChamber,fCathode);
  fPad->cd(); 
  Draw();
}

//______________________________________________________________________________
void AliMUONdisplay::UnZoom()
{
   if (fZooms <= 0) return;
   fZooms--;
   TPad *pad = (TPad*)gPad->GetPadSave();
   pad->Range(fZoomX0[fZooms],fZoomY0[fZooms], fZoomX1[fZooms],fZoomY1[fZooms]);
   pad->Modified();
}

//_____________________________________________________________________________
void AliMUONdisplay::ResetPoints()
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
void AliMUONdisplay::ResetPhits()
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
