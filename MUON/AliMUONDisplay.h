#ifndef ALIMUONDISPLAY_H
#define ALIMUONDISPLAY_H

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

class AliMUONDisplay : /*splaypublic TObject,*/ public AliDisplay {
public:
                     AliMUONDisplay();
                     AliMUONDisplay(Int_t size);
		     AliMUONDisplay(const AliMUONDisplay& display);
		     
   virtual          ~AliMUONDisplay();
   virtual void      Clear(Option_t *option="");
   virtual void      DisplayButtons();
   virtual void      CreateColors();
   virtual void      DisplayColorScale();
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
   virtual void      Draw(Option_t *option="");
   virtual void      DrawClusters();
   virtual void      DrawHits();
   virtual void      DrawCoG();
   virtual void      DrawCoG2();
   virtual void      DrawSegmentation();
   virtual void      DrawTitle(Option_t *option="");
   virtual void      DrawView(Float_t theta, Float_t phi, Float_t psi=0);
   virtual void      DrawP(Float_t,Float_t,Float_t,Float_t,Float_t,Int_t){}
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Int_t             GetZoomMode() {return fZoomMode;}
   Int_t             GetChamber() {return fChamber;}
   Int_t             GetCathode() {return fCathode;}
   virtual void      LoadDigits(Int_t chamber, Int_t cathode);
   virtual void      LoadHits(Int_t chamber);
   virtual void      LoadCoG(Int_t chamber, Int_t cathode);
   virtual void      LoadCoG2(Int_t chamber, Int_t cathode);
   TPad             *Pad() {return fPad;}
   TObjArray        *Points() {return fPoints;}
   TObjArray        *Phits() {return fPhits;}
   TObjArray        *Rpoints() {return fRpoints;}
   TObjArray        *R2points() {return fR2points;}
   TObjArray        *Cpoints() {return fCpoints;}
   virtual void      Paint(Option_t *option="");
   virtual void      SetDrawClusters(Bool_t draw=kTRUE) {fDrawClusters=draw;}   // *MENU*
   virtual void      SetChamberAndCathode(Int_t chamber=1, Int_t cathode=1);    // *MENU*
   virtual void      SetDrawCoG(Bool_t draw=kTRUE) {fDrawCoG=draw;}             // *MENU*
   virtual void      SetRange(Float_t rrange=250., Float_t zrange=1050.);       // *MENU*
   virtual void      SetEvent(Int_t newevent=0);                                // *MENU*   
   virtual void      SetView(Float_t theta=0, Float_t phi=-90, Float_t psi=0);
   virtual void      SetPickMode();
   virtual void      SetZoomMode();
   virtual void      ShowNextEvent(Int_t delta=1);
   virtual void      UnZoom();                                                   // *MENU* 
   virtual void      ResetPoints();
   virtual void      ResetPhits();
   virtual void      ResetRpoints();
   virtual void      ResetR2points();
   virtual void      ResetCpoints();
   virtual void      NextChamber(Int_t delta=1);
   virtual void      NextCathode();
           void      Trigger();
	   AliMUONDisplay& operator = (const AliMUONDisplay& rhs);
	   

private:
   Int_t             fEvent;                // Current event
   Int_t             fChamber;              // Current Chamber
   Int_t             fCathode;              // Current cathode plane
   Int_t             fZoomMode;             // =1 if in zoom mode

   Bool_t            fDrawClusters;         //Flag True if Clusters to be drawn
   Bool_t            fDrawCoG;              //Flag True if CoG to be drawn
                                            //to be drawn
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
   TObjArray        *fRpoints;              //Array of cog points for each cathode
   TObjArray        *fR2points;             //Array of cog points for each cathode
   TObjArray        *fCpoints;              //Array of correlated points for each first cathode
   Int_t            fNextCathode;           //Flagging next cathode
   ClassDef(AliMUONDisplay, 0)   //Utility class to display MUON clusters...
};

#endif








