#ifndef MUON_RECODISPLAY
#define MUON_RECODISPLAY
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

// Authors : M.Gheata, A.Gheata 09/10/00

#include <TApplication.h>
#include <TROOT.h>
#include <TFile.h>
#include <TPolyLine3D.h>
#include <TParticle.h>
#include <AliDisplay.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <AliDetector.h>
#include "AliMUONHit.h"


class AliMUONRecoDisplay:public AliDisplay {

private:
//methods
   Int_t              GetBestMatch(Int_t indr, Float_t tolerance=3.0);
   TClonesArray*      MakePolyLines3D(TClonesArray *tracklist);
   void               MapEvent(Int_t nevent);
   Bool_t             IsReconstructible(Int_t track);
//data members
   Int_t              fEvent;                   // current event number
   AliMUONRecoEvent  *fEvGen;                   // Geant event
   AliMUONRecoEvent  *fEvReco;                  // reconstructed event
   TFile             *fFile;                    // file with reco. event tree
   TTree             *fTree;                    // tree with reco. events
   TClonesArray      *fPolyRecoList;            // list of TPolyLine3D's for reco. tracks
   TClonesArray      *fPolyGenList;             // list of TPolyLine3D's for generated tracks
   TClonesArray      *fRecoTracks;              // list of reco tracks
   TClonesArray      *fGenTracks;               // list of GEANT tracks
   Int_t              fHighlited;               // index of current highlited track
   Double_t           fMinMomentum;             // min. cut of momentum
   Double_t           fMaxMomentum;             // max. cut of momentum
   Bool_t             fPrinted;			// tracks info switch
   Bool_t             fEmpty;                   // true if current reco. event empty

public:
   AliMUONRecoDisplay(Int_t nevent=0);
   virtual             ~AliMUONRecoDisplay();
   virtual void       DrawHits();
   virtual void       DrawView(Float_t theta, Float_t phi, Float_t psi = 0);
   Bool_t             Event(Int_t nevent);
   virtual void       SetDrawHits(Bool_t hits = kTRUE); 	// *MENU*
   virtual void       ShowNextEvent(Int_t delta = 1);
   void               ListTracks();      			// *MENU*
   void               Highlight(Int_t track=0); 		// *MENU*
   void               UnHighlight();               		// *MENU*
   void               CutMomentum(Double_t min=0, Double_t max=999); 	// *MENU*
   void               PolyLineInfo(TClonesArray *line3Dlist);
   void               RecoEfficiency(Int_t first=0, Int_t last=10000);  // *MENU*
   void               XYPlot();                                 // *MENU*
   
   ClassDef(AliMUONRecoDisplay,0)	// MUON reco. event display
};

#endif
