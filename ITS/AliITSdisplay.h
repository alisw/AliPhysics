#ifndef AliITSdisplay_H
#define AliITSdisplay_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliDisplay                                                           //
//                                                                      //
// Utility class to display ALice outline, tracks, hits,..              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//#ifndef ROOT_TObject
#include <TObject.h>
//#endif
#include "AliDisplay.h"

class TCanvas;
class TPad;
class TList;
class TSlider;
class TButton;
class TArc;
class AliITSMapA2;

const Int_t kMAXZOOM = 20;

class AliITSdisplay : public AliDisplay {

private:
   Int_t             fEvent;

   Int_t             fLayer;
   Int_t             fLadder;
   Int_t             fDetector;
   Int_t             fModule;
   Int_t             fNmodules;

   Float_t             fYmodule;
   Option_t            *fOption;   // option for type od det ?

   Int_t             fPdgParticle;
   Int_t             fPdgParent;

   Int_t             fZoomMode;             //=1 if in zoom mode

   Bool_t            fDrawHist;             //Flag True if draw histograms
   AliITSMapA2      *fMap;                

   Bool_t            fDrawClusters;         //Flag True if Clusters to be drawn
   Bool_t            fDrawCoG;              //Flag True if CoG to be drawn
   Bool_t            fDrawCathCor;          //Flag True if correlated point 
                                            //to be drawn
   Float_t           fTheta;                //Viewing angle theta
   Float_t           fPhi;                  //Viewing angle phi
   Float_t           fPsi;                  //Viewving angle psi (rotation on display)
  //   Float_t           fRzone;                //
   Float_t           fRrange;               //Size of view in R
   Float_t           fZrange;               //Size of view along Z
   Float_t           fZoomX0[20];           //Low x range of zoom number i
   Float_t           fZoomY0[20];           //Low y range of zoom number i
   Float_t           fZoomX1[20];           //High x range of zoom number i
   Float_t           fZoomY1[20];           //High y range of zoom number i
   Int_t             fZooms;                //Number of zooms

   Int_t             fHitsCuts;             //Number of hits surviving cuts
   Int_t             fClustersCuts;         //Number of clusters surviving cuts

   TCanvas          *fCanvas;               //Pointer to the display canvas
   TPad             *fTrigPad;              //Pointer to the trigger pad 
   TPad             *fColPad;               //Pointer to the colors pad 
   TPad             *fButtons;              //Pointer to the buttons pad
   TPad             *fPad;                  //Pointer to the event display main pad
   TSlider          *fRangeSlider;          //Range slider
   TButton          *fPickButton;           //Button to activate Pick mode
   TButton          *fZoomButton;           //Button to activate Zoom mode
   TArc             *fArcButton;            //Gren/Red button to show Pick/Zoom mode
   TObjArray        *fPoints;               //Array of points for each cathode
   TObjArray        *fPhits;                //Array of hit points for each chamber
   TObjArray        *fRpoints;              //Array of cog points for each cathode
   TObjArray        *fR2points;              //Array of cog points for each cathode
   TObjArray        *fCpoints;              //Array of correlated points for each first cathode


public:
                     AliITSdisplay();
                     AliITSdisplay(Int_t size);
   virtual          ~AliITSdisplay();
   virtual void      Clear(Option_t *option="");
   virtual void      CreateModuleCanvas(Int_t size);
   virtual void      DisplayButtons();
   virtual void      CreateColors();
   virtual void      DisplayColorScale();
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
   virtual void      Draw(Option_t *option="");
   virtual void      DrawClusters();
   virtual void      DrawHits();
   virtual void      DrawCoG();
   virtual void      DrawCoG2();
   virtual void      DrawCathCor();
   
   virtual void      DrawTitle(Option_t *option="");
   virtual void      DrawView(Float_t theta, Float_t phi, Float_t psi=0);
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Int_t             GetZoomMode() {return fZoomMode;}
   Int_t             GetLayer() {return fLayer;}
   Int_t             GetLadder() {return fLadder;}
   Int_t             GetDetector() {return fDetector;}
   Int_t             GetModule() {return fModule;}
   Int_t             GetEvent() {return fEvent;}
   Float_t           GetYmodule() {return fYmodule;}
   Bool_t            GetDrawHistOpt() {return fDrawHist;}
   AliITSMapA2      *GetMap() {return fMap;}
   void              GetPadCxy(Int_t,Int_t,Float_t *,Float_t *);
   void              GetPadIxy(Int_t &,Int_t &,Float_t,Float_t);
   virtual void      LoadDigits(Int_t module);
   virtual void      LoadHits(Int_t module);
   virtual void      LoadCoG(Int_t chamber, Int_t cathode);
   virtual void      LoadCoG2(Int_t chamber, Int_t cathode);
   virtual void      LoadCathCor(Int_t chamber);
   TPad             *Pad() {return fPad;}
   TObjArray        *Points() {return fPoints;}
   TObjArray        *Phits() {return fPhits;}
   TObjArray        *Rpoints() {return fRpoints;}
   TObjArray        *R2points() {return fR2points;}
   TObjArray        *Cpoints() {return fCpoints;}
   virtual void      Paint(Option_t *option="");
   virtual void      SetModule(Int_t layer=1, Int_t ladder=1, Int_t detector=1); // *MENU*
   virtual void      SetModuleNumber(Int_t module=393); // *MENU*
   virtual void      SetParticle(Int_t code=-211) {fPdgParticle=code;}// *MENU*
   virtual void      SetParent(Int_t code=0) {fPdgParent=code;}// *MENU* 
   virtual void      SetEvent(Int_t newevent=0); // *MENU*   
   virtual void      SetDrawHist(Bool_t draw=kFALSE) {fDrawHist=draw;}   // *MENU*
   virtual void      SetDrawClusters(Bool_t draw=kTRUE) {fDrawClusters=draw;}   // *MENU*
   virtual void      SetDrawCoG(Bool_t draw=kTRUE) {fDrawCoG=draw;}   // *MENU*
   virtual void      SetDrawCathCor(Bool_t draw=kFALSE) {fDrawCathCor=draw;} // *MENU*
   virtual void      SetRange(Float_t rrange=250., Float_t zrange=1050.); // *MENU*
   virtual void      SetView(Float_t theta=-90, Float_t phi=90, Float_t psi=180); // *MENU* 
   virtual void      SetPickMode();
   virtual void      SetZoomMode();
   virtual void      ShowNextEvent(Int_t delta=1);
   virtual void      UnZoom(); // *MENU*
   virtual void      DrawHistograms();
   virtual void      NextModule(Int_t delta=1);
   virtual void      ResetPoints();
   virtual void      ResetPhits();
   virtual void      ResetRpoints();
   virtual void      ResetR2points();
   virtual void      ResetCpoints();

   ClassDef(AliITSdisplay, 0)   //Utility class to display ITS clusters...
};

#endif
