#ifndef ALIRICHDISPLAY_H
#define ALIRICHDISPLAY_H

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
#include <TPolyMarker3D.h>
#include <TAttLine.h>
#include <TAtt3D.h>

//#endif
#include "AliDisplay.h"


class TCanvas;
class TPad;
class TList;
class TSlider;
class TButton;
class TArc;

class AliRICHEllipse;

const Int_t kMAXZOOM = 20;

class AliRICHDisplay : public AliDisplay {
 public:
    AliRICHDisplay();
    AliRICHDisplay(Int_t size);
    virtual          ~AliRICHDisplay();
    virtual void      Clear(Option_t *option="");
    virtual void      DisplayButtons();
    virtual void      CreateColors();
    virtual void      DisplayColorScale();
    virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
    virtual void      Draw(Option_t *option="");
    virtual void      DrawCoG();
    virtual void      DrawRecHits();
    virtual void      DrawCerenkovs();
    virtual void      DrawClusters();
    virtual void      DrawHits();
    virtual void      DrawTitle(Option_t *option="");
    virtual void      DrawView(Float_t theta, Float_t phi, Float_t psi=0);
    virtual void      DrawViewGL();
    virtual void      DrawViewX3D();
    virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
    Int_t             GetZoomMode() {return fZoomMode;}
    Int_t             GetChamber() {return fChamber;}
    Int_t             GetCathode() {return fCathode;}
    virtual void      LoadDigits();
    virtual void      LoadRecHits(Int_t chamber, Int_t cathode);
    virtual void      LoadCoG(Int_t chamber, Int_t cathode);
    virtual void      LoadCerenkovs(Int_t chamber);
    virtual void      LoadHits(Int_t chamber);
    TPad             *Pad() {return fPad;}
    TObjArray        *Points() {return fPoints;}
    TObjArray        *Phits() {return fPhits;}
    TObjArray        *Rpoints() {return fRpoints;}
    TObjArray        *PCerenkovs() {return fPCerenkovs;}
    virtual void      Paint(Option_t *option="");
    virtual void      SetDrawClusters(Bool_t draw=kTRUE) {fDrawClusters=draw;}   // *MENU*
    virtual void      SetDrawCoG(Bool_t draw=kTRUE) {fDrawCoG=draw;}   // *MENU*
    virtual void      SetChamberAndCathode(Int_t chamber=1, Int_t cathode=1); // *MENU*
    virtual void      SetRange(Float_t rrange=250, Float_t zrange=1050); // *MENU*
    virtual void      SetView(Float_t theta, Float_t phi, Float_t psi=0);
    virtual void      SetPickMode();
    virtual void      SetZoomMode();
    virtual void      ShowNextEvent(Int_t delta=1);
    virtual void      UnZoom(); // *MENU*
    virtual void      ResetPoints();
    virtual void      ResetRpoints();
    virtual void      ResetRecpoints();
    virtual void      ResetPhits();
    virtual void      ResetPCerenkovs();
 private:
    Int_t             fChamber;              //Chamber number
    Int_t             fCathode;              //Cathode number
    Int_t             fZoomMode;             //=1 if in zoom mode
    
    Bool_t            fDrawClusters;         //Flag True if Clusters to be drawn
    Bool_t            fDrawCoG;              //Flag True if CoG to be drawn
    Bool_t            fDrawRecHits;          //Flag True if rec hits to be drawn
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
    TObjArray        *fPCerenkovs;           //Array of cerenkov hits for each chamber
    TObjArray        *fRpoints;              //Array of cog points for each cathode 
    TObjArray        *fRecpoints;            //Array of rec points for each cathode 
    ClassDef(AliRICHDisplay, 0)   //Utility class to display RICH clusters...
	
};
#endif
	
