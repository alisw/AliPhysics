#ifndef ALIDISPLAY_H
#define ALIDISPLAY_H
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

#include <TObject.h>

class TCanvas;
class TPad;
class TList;
class TSlider;
class TButton;
class TArc;

static const Int_t kMAXZOOMS = 20;

class AliDisplay : public TObject {
   
public:
                     AliDisplay();
                     AliDisplay(Int_t size);
                     AliDisplay(const AliDisplay &disp);
   virtual          ~AliDisplay();
   virtual Bool_t    AllViews() {return fDrawAllViews;}
   virtual void      Clear(Option_t *option="");
   TSlider          *CutSlider() {return fCutSlider;}
   virtual void      ShowTrack(Int_t trackNumber); // *MENU*
   virtual void      HideTrack(Int_t trackNumber); // *MENU*
           void      Copy(AliDisplay &disp) const;
   virtual void      DisableDetector(const char *name); // *MENU*
   virtual void      DisplayButtons();
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
   virtual void      Draw(Option_t *option="");
   virtual void      DrawAllViews();
   virtual void      DrawHits();
   virtual void      DrawTitle(Option_t *option="");
   virtual void      DrawView(Float_t theta, Float_t phi, Float_t psi=0);
   virtual void      DrawViewGL();
   virtual void      DrawViewX3D();
   virtual void      EnableDetector(const char *name); // *MENU*
   TSlider          *EtaSlider() {return fEtaSlider;}
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Int_t             GetZoomMode() {return fZoomMode;}
   virtual void      LoadPoints();
   TPad             *Pad() {return fPad;}
   virtual void      Paint(Option_t *option="");
   Float_t           PTcut() {return fPTcut;}
   virtual void      SetDrawHits(Bool_t draw=kTRUE) {fDrawHits=draw;}   // *MENU*
   virtual void      SetDrawParticles(Bool_t draw=kTRUE) {fDrawParticles=draw;} // *MENU*
   virtual void      SetPTcut(Float_t ptcut=1.5); // *MENU*
   virtual void      SetRange(Float_t rrange=350, Float_t zrange=350); // *MENU*
   virtual void      SetView(Float_t theta, Float_t phi, Float_t psi=0);
   virtual void      SetPickMode();
   virtual void      SetZoomMode();
   virtual void      ShowNextEvent(Int_t delta=1);
   virtual void      UnZoom(); // *MENU*
   AliDisplay&       operator= (const AliDisplay &disp);
   
protected:
   Int_t             fZoomMode;             //=1 if in zoom mode
   Bool_t            fDrawAllViews;         //Flag True if AllViews selected
   Bool_t            fDrawParticles;        //Flag True if particles to be drawn
   Bool_t            fDrawHits;             //Flag True if Hits to be drawn
   Float_t           fPTcut;                //PT cut to display objects
   Float_t           fTheta;                //Viewing angle theta
   Float_t           fPhi;                  //Viewing angle phi
   Float_t           fPsi;                  //Viewving angle psi (rotation on display)
   Float_t           fRrange;               //Size of view in R
   Float_t           fZrange;               //Size of view along Z
   Float_t           fZoomX0[20];           //Low x range of zoom number i
   Float_t           fZoomY0[20];           //Low y range of zoom number i
   Float_t           fZoomX1[20];           //High x range of zoom number i
   Float_t           fZoomY1[20];           //High y range of zoom number i
   Int_t             fZooms;                //Number of zooms
   Int_t             fHitsCuts;             //Number of hits surviving cuts
   TCanvas          *fCanvas;               //Pointer to the display canvas
   TPad             *fTrigPad;              //Pointer to the trigger pad 
   TPad             *fCutPad;               //Pointer to the momentum cut slider pad 
   TPad             *fEtaPad;               //Pointer to the rapidity cut slider pad 
   TPad             *fButtons;              //Pointer to the buttons pad
   TPad             *fPad;                  //Pointer to the event display main pad
   TSlider          *fCutSlider;            //Momentum cut slider
   TSlider          *fEtaSlider;            //Rapidity slider
   TSlider          *fRangeSlider;          //Range slider
   TButton          *fPickButton;           //Button to activate Pick mode
   TButton          *fZoomButton;           //Button to activate Zoom mode
   TArc             *fArcButton;            //Gren/Red button to show Pick/Zoom mode
   TList            *fFruits;               //List for fruits

   ClassDef(AliDisplay, 0)   //Utility class to display ALICE outline, tracks, hits,..
};

#endif
