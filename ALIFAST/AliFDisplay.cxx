
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFDisplay                                                          //
//                                                                      //
// Utility class to display ALICE outline, tracks, clusters, jets,..    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TView.h>
#include <TArc.h>
#include <TText.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TList.h>
#include <TDiamond.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TMath.h>
#include <X3DBuffer.h>

#include "AliFDisplay.h"
#include "AliFFruit.h"
#include "AliFParticle.h"
#include "AliFast.h"
#include "AliFMCMaker.h"


ClassImp(AliFDisplay)


//_____________________________________________________________________________
AliFDisplay::AliFDisplay() : AliFVirtualDisplay()
{
   fParticle  = 0;
   fFruits    = 0;
}

//_____________________________________________________________________________
AliFDisplay::AliFDisplay(const char *title) : AliFVirtualDisplay()
{

   gAliFast->SetDisplay(this);

   // Initialize display default parameters
   SetPTcut();
   SetPTcutEGMUNU();
   SetGeometry();

   // Set front view by default
   fTheta = 0;
   fPhi   = -90;
   fDrawAllViews  = kFALSE;
   fDrawParticles = kTRUE;

   // Create display canvas
   fCanvas = new TCanvas("Canvas", (char*)title,14,47,740,650);
   fCanvas->SetEditable(kIsNotEditable);

   // Create main display pad
   fPad = new TPad("viewpad", "AliFast display",0.15,0,1,1);
   fPad->Draw();
   fPad->Modified();
   fPad->SetFillColor(1);
   fPad->SetBorderSize(2);

   // Create user interface control pad
   DisplayButtons();
   fCanvas->cd();

   // Create trigger view pad
   Float_t dxtr = 0.15;
   Float_t dytr = 0.45;
   Float_t xt   = 0.3*dxtr;
   Float_t yt   = 0.8*dytr;
   Float_t dyt  = 0.07*dytr;
   Float_t xarc = 0.7*dxtr;
   Float_t rarc = 0.3*dyt;
   fTrigPad = new TPad("TrigPad", "trigger pad",0,0,dxtr,dytr);
   fTrigPad->Range(0,0,dxtr,dytr);
   fTrigPad->Draw();
   fTrigPad->cd();
   fTrigPad->SetFillColor(22);
   fTrigPad->SetBorderSize(2);

   TText *t = new TText();
   t->SetTextFont(61);
   t->SetTextSize(0.2);
   t->SetTextAlign(22);
   t->DrawText(0.5*dxtr, 0.93*dytr,"Trigger");
   t->SetTextSize(0.14);
   t->SetTextAlign(22);
   t->DrawText(xt,yt,      "EM1");
   t->DrawText(xt,yt-dyt,  "PH1");
   t->DrawText(xt,yt-2*dyt,"EM2");
   t->DrawText(xt,yt-3*dyt,"MU1");
   t->DrawText(xt,yt-4*dyt,"MU2");
   t->DrawText(xt,yt-5*dyt,"EMU");
   t->DrawText(xt,yt-6*dyt,"JT1");
   t->DrawText(xt,yt-7*dyt,"JT2");
   t->DrawText(xt,yt-8*dyt,"JT3");
   t->DrawText(xt,yt-9*dyt,"ALL");
   AppendPad(); // append display object as last object to force selection

   fTubin = new TTUBE("tubin","inner tube"," ", fRin, fRin+5, fZin);
   fNodin = new TNode("nodin","ALIAS outline","tubin",0,0,0," ");
   fNodin->SetLineColor(7);
         

   // Create list to support list of fruits
   fFruits = new TList();

   // Create particle manager
   fParticle = new AliFParticle("particle_manager");

   fCanvas->cd();
   fCanvas->Update();

}


//_____________________________________________________________________________
AliFDisplay::~AliFDisplay()
{
   delete fParticle;
   if (fFruits) fFruits->Delete();
   delete fFruits;
}

//_____________________________________________________________________________
void AliFDisplay::Clear(Option_t *)
{
//    Delete graphics temporary objects

   fFruits->Delete();

}

//_____________________________________________________________________________
void AliFDisplay::DisplayButtons()
{
//    Create the user interface buttons

   fButtons = new TPad("buttons", "newpad",0,0.45,0.15,1);
   fButtons->Draw();
   fButtons->SetFillColor(38);
   fButtons->SetBorderSize(2);
   fButtons->cd();

   Int_t butcolor = 33;
   Float_t dbutton = 0.08;
   Float_t y  = 0.96;
   Float_t dy = 0.014;
   Float_t x0 = 0.05;
   Float_t x1 = 0.95;

   TButton *button;
   char *but1 = "gAliFast->Display()->ShowNextEvent(1)";
   button = new TButton("Next",but1,x0,y-dbutton,x1,y);
   button->SetFillColor(38);
   button->Draw();

   y -= dbutton +dy;
   char *but2 = "gAliFast->Display()->ShowNextEvent(-1)";
   button = new TButton("Previous",but2,x0,y-dbutton,x1,y);
   button->SetFillColor(38);
   button->Draw();

   y -= dbutton +dy;
   char *but3 = "gAliFast->Display()->SetView(90,-90)";
   button = new TButton("Top View",but3,x0,y-dbutton,x1,y);
   button->SetFillColor(butcolor);
   button->Draw();

   y -= dbutton +dy;
   char *but4 = "gAliFast->Display()->SetView(90,0)";
   button = new TButton("Side View",but4,x0,y-dbutton,x1,y);
   button->SetFillColor(butcolor);
   button->Draw();

   y -= dbutton +dy;
   char *but5 = "gAliFast->Display()->SetView(0,-90)";
   button = new TButton("Front View",but5,x0,y-dbutton,x1,y);
   button->SetFillColor(butcolor);
   button->Draw();

   y -= dbutton +dy;
   char *but6 = "gAliFast->Display()->DrawAllViews()";
   button = new TButton("All Views",but6,x0,y-dbutton,x1,y);
   button->SetFillColor(butcolor);
   button->Draw();

   y -= dbutton +dy;
   char *but7 = "gAliFast->Display()->DrawViewGL()";
   button = new TButton("OpenGL",but7,x0,y-dbutton,x1,y);
   button->SetFillColor(38);
   button->Draw();

   y -= dbutton +dy;
   char *but8 = "gAliFast->Display()->DrawViewX3D()";
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
   diamond->AddText("AliFAST");
   diamond->AddText("... ");
   diamond->AddText(" ");
}

//_____________________________________________________________________________

Int_t AliFDisplay::DistancetoPrimitive(Int_t px, Int_t py)
{
// Compute distance from point px,py to objects in event

   if (gPad == fTrigPad) {gPad->SetCursor(kCross); return 0;}

   const Int_t big = 9999;
   Int_t dist = big;
   Float_t xmin = gPad->GetX1();
   Float_t xmax = gPad->GetX2();
   Float_t dx   = 0.05*(xmax - xmin);
   Float_t x    = gPad->AbsPixeltoX(px);
   if (x < xmin+dx || x > xmax-dx) return dist;

    // scan list of particles
   dist = fParticle->DistancetoPrimitive(px, py);
   if (dist <= 0) return 0;

    // scan list of fruits
   TIter nextf(fFruits);
   AliFFruit *fruit;
   while((fruit=(AliFFruit*)nextf())) {
      dist = fruit->DistancetoPrimitive(px, py);
      if (dist < 5) {
         gPad->SetSelected(fruit->Fruit());
         gPad->SetCursor(kCross);
         return 0;
      }
   }

    // scan list of detectors (currently only one tube)
   dist = fNodin->DistancetoPrimitive(px, py);
   if (gPad->GetCanvas()->GetSelected() == gPad->GetView()) {
      gPad->SetSelected(this);
   }
   return 0;
}

//_____________________________________________________________________________
void AliFDisplay::Draw(Option_t *)
{
//    Insert current event in graphics pad list

   if (fDrawAllViews) {
      DrawAllViews();
      return;
   }

   fPad->cd();

   DrawView(fTheta, fPhi);

   // Display the event number and title
   fPad->cd();
   DrawTitle();
}

//_____________________________________________________________________________
void AliFDisplay::DrawAllViews()
{
//    Draw front,top,side and 30 deg views

   fDrawAllViews = kTRUE;
   fPad->cd();
   fPad->SetFillColor(15);
   fPad->Clear();
   fPad->Divide(2,2);

   // draw 30 deg view
   fPad->cd(1);
   DrawView(30, 30);
   DrawTitle();

   // draw front view
   fPad->cd(2);
   DrawView(0, -90);
   DrawTitle("Front");

   // draw top view
   fPad->cd(3);
   DrawView(90, -90);
   DrawTitle("Top");

   // draw side view
   fPad->cd(4);
   DrawView(90, 0);
   DrawTitle("Side");

   fPad->cd(2);
}

//_____________________________________________________________________________
void AliFDisplay::DrawTitle(Option_t *option)
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
      sprintf(ptitle,"Pythia event: %d, Run:%d",gAliFast->Event(), gAliFast->Run());
      title->AddText(ptitle);
      sprintf(ptitle,"Pythia Mode: %s",gAliFast->MCMaker()->GetTitle());
      title->AddText(ptitle);
   } else {
      TPaveLabel *label = new TPaveLabel(xmin +0.01*dx, ymax-0.07*dy, xmin +0.2*dx, ymax-0.01*dy,option);
      label->SetBit(kCanDelete);
      label->SetFillColor(42);
      label->Draw();
   }
}

//_____________________________________________________________________________
void AliFDisplay::DrawView(Float_t theta, Float_t phi)
{
//    Draw a view of ALIAS

   gPad->SetFillColor(1);
   // Display ALIAS outline
   gPad->Clear();

   Int_t iret;
   TView *view = new TView(1);
   view->SetRange(-fRin, -fRin, -fZin, fRin, fRin, fZin);

   fNodin->Draw("same");

    // add itself to the list
   AppendPad();
   
   //Loop on all makers to add their products to the pad
   TIter next(gAliFast->Makers());
   AliFMaker *maker;
   while ((maker = (AliFMaker*)next())) {
      maker->Draw();
   }
   view->SetView(phi, theta, 0, iret);
}

//_____________________________________________________________________________
void AliFDisplay::DrawViewGL()
{
//    Draw current view using OPENGL

   TPad *pad = (TPad*)gPad->GetPadSave();
   pad->cd();
   TView *view = pad->GetView();
   if (!view) return;
   pad->x3d("OPENGL");
}

//_____________________________________________________________________________
void AliFDisplay::DrawViewX3D()
{
//    Draw current view using X3D

   TPad *pad = (TPad*)gPad->GetPadSave();
   pad->cd();
   TView *view = pad->GetView();
   if (!view) return;
   pad->x3d();
}

//______________________________________________________________________________
void AliFDisplay::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
//*-*-*-*-*-*-*-*-*-*-*Execute action corresponding to one event*-*-*-*
//*-*                  =========================================

   if (gPad->GetView()) {
      gPad->GetView()->ExecuteRotateView(event, px, py);
   }
}

//_____________________________________________________________________________
void AliFDisplay::GetEvent(Int_t event)
{
//    Read event in memory

   gAliFast->GetTreeEvent(event);

   Draw();
}

//_____________________________________________________________________________
void AliFDisplay::Paint(Option_t *)
{
//    Paint miscellaneous items

}

//_____________________________________________________________________________
void AliFDisplay::PaintFruit(TObject *obj, Float_t eta, Float_t phi, Float_t
pt, Int_t type, Option_t *option)
{
//    Display fruit from obj

   AliFFruit *fruit = new AliFFruit(obj, eta, phi, pt, type);
   fFruits->Add(fruit);
   fruit->Paint(option);
}

//_____________________________________________________________________________
void AliFDisplay::PaintParticles(Option_t *option)
{
   if (fDrawParticles) fParticle->Paint(option);
}

//_____________________________________________________________________________
void AliFDisplay::SetGeometry(Float_t rin)
{
//  Set ALIAS in/out outline parameters

   fRin  = rin;
   fRout = 1.2*rin;
   fZin  = 600;
   fZout = 680;
}

//_____________________________________________________________________________
void AliFDisplay::SetPTcut(Float_t ptcut)
{
   fPTcut = ptcut;

   if (fDrawAllViews) {
      fPad->cd(1); gPad->Modified();
      fPad->cd(2); gPad->Modified();
      fPad->cd(3); gPad->Modified();
      fPad->cd(4); gPad->Modified();
      fPad->cd();
   }
}

//_____________________________________________________________________________
void AliFDisplay::SetPTcutEGMUNU(Float_t ptcut)
{
   fPTcutEGMUNU = ptcut;

   if (fDrawAllViews) {
      fPad->cd(1); gPad->Modified();
      fPad->cd(2); gPad->Modified();
      fPad->cd(3); gPad->Modified();
      fPad->cd(4); gPad->Modified();
      fPad->cd();
   }
}

//_____________________________________________________________________________
void AliFDisplay::SetView(Float_t theta, Float_t phi)
{
//  change viewing angles for current event

   fPad->cd();
   fDrawAllViews = kFALSE;
   fPhi   = phi;
   fTheta = theta;
   Int_t iret;

   TView *view = gPad->GetView();
   if (view) view->SetView(fPhi, fTheta, 0, iret);
   else      Draw();

   gPad->Modified();
}

//_____________________________________________________________________________
void AliFDisplay::ShowNextEvent(Int_t delta)
{
//  Display (current event_number+delta)
//    delta =  1  shown next event
//    delta = -1 show previous event

  if (delta) {
     gAliFast->Clear();
     Int_t current_event = gAliFast->Event();
     Int_t new_event     = current_event + delta;
     gAliFast->GetTreeEvent(new_event); 
   }
  fPad->cd(); 
  Draw();
}

//______________________________________________________________________________
void AliFDisplay::SizeFruit() const
{
   const Int_t npoints = 2;
   gSize3D.numPoints += npoints;
   gSize3D.numSegs   += (npoints-1);
   gSize3D.numPolys  += 0;
}

//______________________________________________________________________________
void AliFDisplay::SizeParticles() const
{
   if (fDrawParticles)  fParticle->SizeParticles();
}










