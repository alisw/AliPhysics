
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
#include <TGXW.h>
#include <TMath.h>
#include <TMatrix.h>
#include <X3DBuffer.h>

#include "AliRun.h"
#include "AliDetector.h"
#include "AliMUON.h"
#include "AliMUONConst.h"
#include "AliMUONdisplay.h"
#include "AliMUONpoints.h"
#include "TParticle.h"



ClassImp(AliMUONdisplay)


//_____________________________________________________________________________
AliMUONdisplay::AliMUONdisplay()
{
   fPoints = 0;
   fPhits = 0;
   fRpoints = 0;
   fR2points = 0;
   fCpoints = 0;
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
<img src="gif/aliMUONdisplay.gif">
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
   fDrawCathCor  = kTRUE;
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
   text->SetTextSize(0.2);
   text->SetTextAlign(22);
   
   AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
   AliMUONchamber *iChamber = &(MUON->Chamber(fChamber-1));
   AliMUONresponse * response=iChamber->GetResponseModel();
   Int_t adcmax= (Int_t) response->MaxAdc();
   

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

   DrawView(fTheta, fPhi, fPsi);   
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

      
   fClustersCuts = 0;
      points = Points();
      if (!points) return;
      ndigits = points->GetEntriesFast();
      for (digit=0;digit<ndigits;digit++){
         pm = (AliMUONpoints*)points->UncheckedAt(digit);
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
      for (track=0;track<ntracks;track++) {
         pm = (AliMUONpoints*)points->UncheckedAt(track);
	 if (!pm) continue;
         pm->Draw();
         fHitsCuts += pm->GetN();
      }
}


//_____________________________________________________________________________
void AliMUONdisplay::DrawCoG()
{
//    Draw hits for MUON chambers
    if (!fDrawCoG) return;
   LoadCoG(fChamber,fCathode);

   Int_t ncog, icog;
   TObjArray *points;
   AliMUONpoints *pm;

      points = Rpoints();
      if (!points) return;
      ncog = points->GetEntriesFast();
      for (icog=0;icog<ncog;icog++) {
         pm = (AliMUONpoints*)points->UncheckedAt(icog);
	 if (!pm) continue;
         pm->Draw();
      }
}
void AliMUONdisplay::DrawCoG2()
{
//    Draw hits for MUON chambers

   if (!fDrawCoG) return;
  

  if (fCathode==1) {
     LoadCoG2(fChamber,2);
  } else if (fCathode==2) {
     LoadCoG2(fChamber,1);
  }

   Int_t ncog, icog;
   TObjArray *points;
   AliMUONpoints *pm;

      points = R2points();
      if (!points) return;
      ncog = points->GetEntriesFast();
      for (icog=0;icog<ncog;icog++) {
         pm = (AliMUONpoints*)points->UncheckedAt(icog);
	 if (!pm) continue;
         pm->Draw();
      }

}
//_____________________________________________________________________________
void AliMUONdisplay::DrawCathCor()
{
//    Draw hits for MUON chambers
    
   if (!fDrawCathCor) return;

   LoadCathCor(fChamber);

   Int_t ncog, icog;
   TObjArray *points;
   AliMUONpoints *pm;

      points = Cpoints();
      if (!points) return;
      ncog = points->GetEntriesFast();
      for (icog=0;icog<ncog;icog++) {
         pm = (AliMUONpoints*)points->UncheckedAt(icog);
	 if (!pm) continue;
         pm->Draw();
      }
}
//_____________________________________________________________________________
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

   Int_t iret=0;
   TView *view = new TView(1);

   Float_t range = fRrange*fRangeSlider->GetMaximum();
   view->SetRange(-range,-range,-range,range, range, range);
   fZoomX0[0] = -1;
   fZoomY0[0] = -1;
   fZoomX1[0] =  1;
   fZoomY1[0] =  1;
   fZooms = 0;

// Display MUON Chamber Geometry
// gAlice->GetGeometry()->Draw("same");
   char NodeName[7];
   sprintf(NodeName,"MUON%d",100+fChamber);
   printf("Node name %s", NodeName);
   
   TNode *node1=gAlice->GetGeometry()->GetNode(NodeName);
   if (node1) node1->Draw("same");  
// ok if I rotate the chamber in a proper way

//add clusters to the pad
   DrawClusters();
   DrawHits();
   DrawCoG();
   DrawCoG2();
   DrawCathCor();
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

   if (chamber > 10) return;

   fChamber=chamber;
   fCathode=cathode;

   ResetPoints();

   AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
   AliMUONchamber*  iChamber;
   AliMUONsegmentation*  segmentation;

   TClonesArray *MUONdigits  = MUON->DigitsAddress(chamber-1);
   if (MUONdigits == 0) return;

   gAlice->ResetDigits();

   Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
   gAlice->TreeD()->GetEvent(nent-2+cathode-1);
  //gAlice->TreeD()->GetEvent(cathode);
   Int_t ndigits = MUONdigits->GetEntriesFast();
   if (ndigits == 0) return;
   if (fPoints == 0) fPoints = new TObjArray(ndigits);
	  
   iChamber = &(MUON->Chamber(chamber-1));
   segmentation=iChamber->GetSegmentationModel(cathode);
   //Float_t dpxbig  = segmentation->Dpx()/2.;
   //   Float_t dpy  = segmentation->Dpy()/2.;
   Float_t zpos=iChamber->ZPosition();  
   AliMUONdigit  *mdig;
   AliMUONpoints *points = 0;
   TMarker3DBox *marker=0;
   //   TMatrix *matrix;
   //
   //loop over all digits and store their position
   //    points = new AliMUONpoints(ndigits);
   Int_t npoints=1;
   for (Int_t digit=0;digit<ndigits;digit++) {
        mdig    = (AliMUONdigit*)MUONdigits->UncheckedAt(digit);
        //
        // First get all needed parameters
        //
        Int_t charge=mdig->fSignal;
        Int_t index=Int_t(TMath::Log(charge)/(TMath::Log(adc_satm)/22));
        Int_t color=261+index;
        if (color>282) color=282;
	// get the center of the pad - add on x and y half of pad size
	Float_t xpad, ypad;
	segmentation->GetPadCxy(mdig->fPadX, mdig->fPadY,xpad, ypad);
	//	printf("xpad,ypad,zpos,dpx,dpy %f %f %f %f %f\n",xpad,ypad,zpos,dpx,dpy);
        Int_t isec=segmentation->Sector(mdig->fPadX, mdig->fPadY);
	//        printf(" isec %d \n",isec);
        Float_t dpx=segmentation->Dpx(isec)/2;
        Float_t dpy=segmentation->Dpy(isec)/2;
	//        printf(" dpx %f \n",dpx);
	Int_t nPara, offset;
        segmentation->GetNParallelAndOffset(mdig->fPadX,mdig->fPadY,
		&nPara,&offset);
	//
	// Then set the objects
	//
        points = new AliMUONpoints(npoints);
	fPoints->AddAt(points,digit);

        points->SetMarkerColor(color);
        points->SetMarkerStyle(21);
        points->SetMarkerSize(0.5);
        points->SetParticle(-1);
        points->SetHitIndex(-1);
        points->SetTrackIndex(-1);
        points->SetDigitIndex(digit);
        points->SetPoint(0,xpad,ypad,zpos);	
	for (Int_t imark=0;imark<nPara; imark++)
	  {
 	  segmentation->GetPadCxy(mdig->fPadX + imark*offset, mdig->fPadY,xpad, ypad);
	  marker=new TMarker3DBox(xpad,ypad,zpos,dpx,dpy,0,0,0);
	  marker->SetLineColor(2);
	  marker->SetFillStyle(1001);
          marker->SetFillColor(color);
          marker->SetRefObject((TObject*)points);
          points->Set3DMarker(imark, marker);

	  }
   }
}
//___________________________________________
void AliMUONdisplay::LoadCoG(Int_t chamber, Int_t cathode)
{
// Read raw clusters info and store x,y,z info in arrays fRpoints
// Loop on all detectors

   if (chamber > 10) return;

   ResetRpoints();

   AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
   AliMUONchamber*  iChamber;

   TClonesArray *MUONrawclust  = MUON->RawClustAddress(chamber-1);
   if (MUONrawclust == 0) return;

   MUON->ResetRawClusters();


   Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
   gAlice->TreeR()->GetEvent(nent-2+cathode-1);
   //gAlice->TreeR()->GetEvent(cathode);
   Int_t nrawcl = MUONrawclust->GetEntriesFast();
   if (nrawcl == 0) return;
   if (fRpoints == 0) fRpoints = new TObjArray(nrawcl);
	  
   iChamber = &(MUON->Chamber(chamber-1));
   Float_t zpos=iChamber->ZPosition();  
   AliMUONRawCluster  *mRaw;
   AliMUONpoints *points = 0;
   //
   //loop over all raw clusters and store their position
   points = new AliMUONpoints(nrawcl);
   for (Int_t iraw=0;iraw<nrawcl;iraw++) {
  	mRaw   = (AliMUONRawCluster*)MUONrawclust->UncheckedAt(iraw);
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
}
//___________________________________________
void AliMUONdisplay::LoadCoG2(Int_t chamber, Int_t cathode)
{
// Read raw clusters info and store x,y,z info in arrays fRpoints
// Loop on all detectors

   if (chamber > 10) return;

   ResetR2points();

   AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
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
}
//___________________________________________
void AliMUONdisplay::LoadCathCor(Int_t chamber)
{
// Read correlation info and store x,y,z info in arrays fCpoints
// Loop on all detectors

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
}

//___________________________________________
void AliMUONdisplay::LoadHits(Int_t chamber)
{
// Read hits info and store x,y,z info in arrays fPhits
// Loop on all detectors

   if (chamber > 10) return;

   fChamber=chamber;
 
   ResetPhits();

   AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
   AliMUONchamber*  iChamber;

   iChamber = &(MUON->Chamber(chamber-1));
   Float_t zpos=iChamber->ZPosition();

   Int_t ntracks = (Int_t)gAlice->TreeH()->GetEntries();
   //Int_t ntrks = gAlice->GetNtrack();


   Int_t nthits=0;
    for (Int_t track=0; track<ntracks;track++) {
      gAlice->ResetHits();
      gAlice->TreeH()->GetEvent(track);
         TClonesArray *MUONhits  = MUON->Hits();
         if (MUONhits == 0) return;
         nthits += MUONhits->GetEntriesFast();
   } 
   if (fPhits == 0) fPhits = new TObjArray(nthits);
   //printf("nthits %d \n",nthits);

   // old stuff 
   //
   //if (fPhits == 0) fPhits = new TObjArray(ntracks);
   /*
    TVector *xp = new TVector(20);
    TVector *yp = new TVector(20);
    //    TVector *zp = new TVector(20);
    TVector *ptrk = new TVector(20);
    TVector *phit = new TVector(20);
   */
   // end old stuff

    Int_t nhold=0;
    for (Int_t track=0; track<ntracks;track++) {
      gAlice->ResetHits();
      gAlice->TreeH()->GetEvent(track);
         TClonesArray *MUONhits  = MUON->Hits();
         if (MUONhits == 0) return;
         Int_t nhits = MUONhits->GetEntriesFast();
         if (nhits == 0) continue;
         AliMUONhit *mHit;
         AliMUONpoints *points = 0;
         Int_t npoints=1;
         for (Int_t hit=0;hit<nhits;hit++) {
            mHit = (AliMUONhit*)MUONhits->UncheckedAt(hit);
            Int_t nch  = mHit->fChamber;              // chamber number
            if (nch != chamber) continue;
	    //
	    // Retrieve info and set the objects
	    //
	    points = new AliMUONpoints(npoints);
	    fPhits->AddAt(points,nhold+hit);
            points->SetMarkerColor(kRed);
            points->SetMarkerStyle(5);
            points->SetMarkerSize(1.);
            points->SetParticle(mHit->fTrack);
            points->SetHitIndex(hit);
            points->SetTrackIndex(track);
            points->SetDigitIndex(-1);
	    points->SetPoint(0,mHit->fX,mHit->fY,zpos);
	 }
	 nhold+=nhits;


	 // old stuff
	 /*
         Int_t npoints=0;
         for (Int_t hit=0;hit<nhits;hit++) {
            mHit = (AliMUONhit*)MUONhits->UncheckedAt(hit);
            Int_t nch  = mHit->fChamber;              // chamber number
            if (nch != chamber) continue;

            (*xp)(npoints)=mHit->fX;
            (*yp)(npoints)=mHit->fY;
	    //            (*zp)(npoints)=mHit->fZ;
            (*ptrk)(npoints)=Float_t(mHit->GetTrack());
            //(*ptrk)(npoints)=Float_t(mHit->fTrack);
	    (*phit)(npoints)=Float_t(hit);
	     printf("hit,(*phit)(npoints), track, trk, fTrack ipart %d %f %d %d %d %f\n",hit,(*phit)(npoints),track,mHit->GetTrack(),mHit->fTrack,mHit->fParticle);
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
            Int_t index=points->GetIndex();
            points->SetHitIndex(Int_t((*phit)(p)));
            points->SetTrackIndex(track);
	    	    printf("p, index, Int_t((*ptrk)(p)), hit, track  %d %d %d %d %d \n",p, index,Int_t((*ptrk)(p)),Int_t((*phit)(p)),track);
            points->SetDigitIndex(-1);
	    points->SetPoint(p,(*xp)(p),(*yp)(p),zpos);
	    //            points->SetPoint(p,(*xp)(p),(*yp)(p),(*zp)(p));
	 }
	 xp->Zero();
	 yp->Zero();
	 //      zp->Zero();
	 ptrk->Zero();
	 phit->Zero();
	 fPhits->AddAt(points,track);
	 //            Int_t np=points->GetN();
	 //            printf("np %d \n",np);


	 */
	 // end old stuff

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
     fEvent=new_event;
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
//_____________________________________________________________________________
void AliMUONdisplay::ResetRpoints()
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
void AliMUONdisplay::ResetR2points()
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
void AliMUONdisplay::ResetCpoints()
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











