#ifndef ALIMUONDISPLAY_H
#define ALIMUONDISPLAY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliDisplay                                                           //
//                                                                      //
// Utility class to display ALice outline, tracks, hits,..              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliDisplay.h"

class AliLoader;
class AliMUONData;
class AliMUONTrack;
class TCanvas;
class TPad;
class TList;
class TSlider;
class TButton;
class TArc;

class AliMUONDisplay : public AliDisplay 
{
public:
                     AliMUONDisplay();
                     AliMUONDisplay(Int_t size, AliLoader * loader=0x0);
		     
   virtual          ~AliMUONDisplay();
   virtual void      Clear(Option_t *option="");
   virtual void      DisplayButtons();
   virtual void      CreateColors() const;
   virtual void      DisplayColorScale();
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
   virtual void      DrawReco(Option_t *option="");
   virtual void      PrintTrack(Int_t iRecTracks, AliMUONTrack *recTrack); 
   virtual void      PrintKinematics(); 
   virtual void      Draw(Option_t *option="");
   virtual void      DrawChamber();
   virtual void      DrawClusters();
   virtual void      DrawHits();
   virtual void      DrawCoG();
   virtual void      DrawTracks();
   virtual void      DrawSegmentation();
   virtual void      DrawTitle(Option_t *option="");
   virtual void      DrawView(Float_t theta, Float_t phi, Float_t psi=0);
   virtual void      DrawGlobalView(Float_t theta, Float_t phi, Float_t psi=0);
   virtual void      DrawP(Float_t,Float_t,Float_t,Float_t,Float_t,Int_t){}
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Int_t             GetZoomMode() const {return fZoomMode;}
   Int_t             GetChamber() const {return fChamber;}
   Int_t             GetCathode() const {return fCathode;}

   AliMUONData*      GetMUONData() {return fMUONData;}
   AliLoader*        GetLoader()  {return fLoader;}

   virtual void      LoadDigits(Int_t chamber, Int_t cathode);
   virtual void      LoadHits(Int_t chamber);
   virtual void      LoadCoG(Int_t chamber, Int_t cathode);
   virtual void      LoadTracks();
   TPad             *Pad() {return fPad;}
   TObjArray        *Points() {return fPoints;}
   TObjArray        *Phits() {return fPhits;}
   TObjArray        *Rpoints() {return fRpoints;}
   virtual void      Paint(Option_t *option="");
   virtual void      SetDrawClusters(Bool_t draw=kTRUE) {fDrawClusters=draw;}   // *MENU*
   virtual void      SetChamberAndCathode(Int_t chamber=1, Int_t cathode=1);    // *MENU*
   virtual void      SetDrawCoG(Bool_t draw=kTRUE) {fDrawCoG=draw;}             // *MENU*
   virtual void      SetDrawTracks(Bool_t draw=kTRUE) {fDrawTracks=draw;}       // *MENU*
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
   virtual void      NextChamber(Int_t delta=1);
   virtual void      NextCathode();
           void      Trigger();

protected:
		     AliMUONDisplay(const AliMUONDisplay& display);
   AliMUONDisplay&   operator = (const AliMUONDisplay& rhs);
	   

private:
   Int_t             fEvent;                // Current event
   Int_t             fChamber;              // Current Chamber
   Int_t             fCathode;              // Current cathode plane
   Bool_t            fDrawClusters;         //Flag True if Clusters to be drawn
   Bool_t            fDrawCoG;              //Flag True if CoG to be drawn
   Bool_t            fDrawTracks;           //Flag True if tracks to be drawn
                                           
   Int_t             fClustersCuts;         //Number of clusters surviving cuts
   TPad             *fColPad;               //Pointer to the colors pad 
   TObjArray        *fPoints;               //Array of points for each cathode
   TObjArray        *fPhits;                //Array of hit points for each chamber
   TObjArray        *fRpoints;              //Array of cog points for each cathode
   Int_t            fNextCathode;           //Flagging next cathode

   AliLoader*       fLoader;                //! MUON loader to get data
   AliMUONData*     fMUONData;              //! Data container for MUON subsystem 

   ClassDef(AliMUONDisplay, 0)   //Utility class to display MUON events
};

#endif








