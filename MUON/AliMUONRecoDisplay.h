// Authors : M.Gheata, A.Gheata 09/10/00
#ifndef MUON_RECDISPLAY
#define MUON_RECDISPLAY

//////////////////////////////////////////////////////////////////////
//                                                                  //
// AliMUONRecoDisplay						    //
//								    //
// This class subclasses AliDisplay and provides display of         //
// reconstructed tracks with following functionality : 		    //
//	- front/top/side/3D display of MUON reconstructed tracks    //
//        as polylines ;                                            //
//	- context menu activated when the main pad is right-clicked //
//	The context menu contains following functions :		    //
//	* SetDrawHits()	- switches on or off Geant hits ;	    //
//	* CutMomentum()	- displays only tracks within Pmin - Pmax   //
//	* ListTracks()	- prints ID and momentum info. for all	    //
//	tracks within momentum range Pmin,Pmax ;		    //
//	* Highlight()	- shows only one selected reco. track	    //
//	and its best matching Geant track;			    //
//	* UnHighlight()	- self explaining;			    //
//	* RecoEfficiency() - compute reco. efficiency for all events//
//        from galice.root file; also fake track percentage; make   //
//        plots for momentum precision                              //
//      * XYPlot()      - make X-Y plots of reconstructed and       //
//        generated tracks in all chambers                          //
//								    //
//      Starting : generate and reconstruct events, then use the    //
//                 MUONrecodisplay.C macro                          //
//                                                                  //
//////////////////////////////////////////////////////////////////////

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

public:
   AliMUONRecoDisplay(Int_t nevent=0);
   virtual             ~AliMUONRecoDisplay();
   virtual void       DrawHits();
   virtual void       DrawView(Float_t theta, Float_t phi, Float_t psi = 0);
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
