/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

/*****************************************************************
  AliFlowEventSimple: A simple event 
  for flow analysis                  
                                     
  origin: Naomi van der Kolk (kolk@nikhef.nl)           
          Ante Bilandzic     (anteb@nikhef.nl)         
          Raimond Snellings  (Raimond.Snellings@nikhef.nl)    
  mods:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#ifndef ALIFLOWEVENTSIMPLE_H
#define ALIFLOWEVENTSIMPLE_H

#include "TObject.h"
#include "TParameter.h"
#include "TMath.h"
class TTree;
class TF1;
class AliFlowVector;
class AliFlowTrackSimple;
class AliFlowTrackSimpleCuts;

class AliFlowEventSimple: public TObject {

 public:

  enum ConstructionMethod {kEmpty,kGenerate};

  AliFlowEventSimple();
  AliFlowEventSimple( Int_t nParticles,
                      ConstructionMethod m=kEmpty,
                      TF1* ptDist=NULL,
                      Double_t phiMin=0.0,
                      Double_t phiMax=TMath::TwoPi(),
                      Double_t etaMin=-1.0,
                      Double_t etaMax= 1.0 );
  AliFlowEventSimple(TTree* anInput, const AliFlowTrackSimpleCuts* rpCuts, const AliFlowTrackSimpleCuts* poiCuts);
  AliFlowEventSimple(const AliFlowEventSimple& anEvent);
  AliFlowEventSimple& operator=(const AliFlowEventSimple& anEvent);
  virtual  ~AliFlowEventSimple();

  Bool_t  IsFolder() const {return kTRUE;};
  void    Browse(TBrowser *b); 
  void    Print(Option_t* option = "") const;      //method to print stats
  
  Int_t    NumberOfTracks() const                   { return fNumberOfTracks; }
  Int_t    GetReferenceMultiplicity() const         { return fReferenceMultiplicity; }
  void     SetReferenceMultiplicity( Int_t m )      { fReferenceMultiplicity = m; }
  Int_t    GetEventNSelTracksRP() const             { return fNumberOfRPs; } 
  void     SetEventNSelTracksRP(Int_t nr)           { fNumberOfRPs = nr; }  
  Int_t    GetEventNSelTracksPOI() const            { return fNumberOfPOIs; } 
  void     SetEventNSelTracksPOI(Int_t np)          { fNumberOfPOIs = np; }  
  Int_t    GetNumberOfRPs() const                   { return fNumberOfRPs; }
  void     SetNumberOfRPs( Int_t nr )               { fNumberOfRPs=nr; }
  Int_t    GetNumberOfPOIs() const                  { return fNumberOfPOIs; }
  void     SetNumberOfPOIs( Int_t np )              { fNumberOfPOIs=np; }

  void     SetUseGlauberMCSymmetryPlanes()          { fUseGlauberMCSymmetryPlanes = kTRUE; }
  void     SetUseExternalSymmetryPlanes(TF1 *gPsi1Psi3 = 0x0,
					TF1 *gPsi2Psi4 = 0x0,
					TF1 *gPsi3Psi5 = 0x0);
  void     SetPsi1(Double_t gPsi1)                  { fPsi1 = gPsi1; }
  void     SetPsi2(Double_t gPsi2)                  { fPsi2 = gPsi2; }
  void     SetPsi3(Double_t gPsi3)                  { fPsi3 = gPsi3; }
  void     SetPsi4(Double_t gPsi4)                  { fPsi4 = gPsi4; }
  void     SetPsi5(Double_t gPsi5)                  { fPsi5 = gPsi5; }
  Double_t GetPsi1() const                          { return fPsi1; }
  Double_t GetPsi2() const                          { return fPsi2; }
  Double_t GetPsi3() const                          { return fPsi3; }
  Double_t GetPsi4() const                          { return fPsi4; }
  Double_t GetPsi5() const                          { return fPsi5; }

  Double_t GetMCReactionPlaneAngle() const          { return fMCReactionPlaneAngle; }
  void     SetMCReactionPlaneAngle(Double_t fPhiRP) { fMCReactionPlaneAngle=fPhiRP; fMCReactionPlaneAngleIsSet=kTRUE; }
  Bool_t   IsSetMCReactionPlaneAngle() const        { return fMCReactionPlaneAngleIsSet; }
  void     SetAfterBurnerPrecision(Double_t p)      { fAfterBurnerPrecision=p; }
  Double_t GetAfterBurnerPrecision() const          { return fAfterBurnerPrecision; }
  void     SetUserModified(Bool_t s=kTRUE)          { fUserModified=s; }
  Bool_t   IsUserModified() const                   { return fUserModified; }
  void     SetShuffleTracks(Bool_t b)               {fShuffleTracks=b;}
  void     ShuffleTracks();

  void ResolutionPt(Double_t res);
  void TagSubeventsInEta(Double_t etaMinA, Double_t etaMaxA, Double_t etaMinB, Double_t etaMaxB );
  void TagSubeventsByCharge();
  void TagRP(const AliFlowTrackSimpleCuts* cuts );
  void TagPOI(const AliFlowTrackSimpleCuts* cuts );
  void CloneTracks(Int_t n);
  void AddV1( Double_t v1 );
  void AddV2( Double_t v2 );
  void AddV3( Double_t v3 );
  void AddV4( Double_t v4 );
  void AddV5( Double_t v5 );
  void AddFlow( Double_t v1, Double_t v2, Double_t v3, Double_t v4, Double_t v5 );
  void AddFlow(Double_t v1, Double_t v2, Double_t v3, Double_t v4, Double_t v5,
               Double_t rp1, Double_t rp2, Double_t rp3, Double_t rp4, Double_t rp5 );
  void AddV2( TF1* ptDepV2 );
  void DefineDeadZone( Double_t etaMin, Double_t etaMax, Double_t phiMin, Double_t phiMax );
  Int_t CleanUpDeadTracks();
  void ClearFast();
 
  static TF1* SimplePtSpectrum();
  static TF1* SimplePtDepV2();

  AliFlowTrackSimple* GetTrack(Int_t i);
  void AddTrack( AliFlowTrackSimple* track ); 
  AliFlowTrackSimple* MakeNewTrack();
 
  AliFlowVector GetQ(Int_t n=2, TList *weightsList=NULL, Bool_t usePhiWeights=kFALSE, Bool_t usePtWeights=kFALSE, Bool_t useEtaWeights=kFALSE);
  void Get2Qsub(AliFlowVector* Qarray, Int_t n=2, TList *weightsList=NULL, Bool_t usePhiWeights=kFALSE, Bool_t usePtWeights=kFALSE, Bool_t useEtaWeights=kFALSE);  

 protected:
  virtual void Generate( Int_t nParticles,
                         TF1* ptDist=NULL,
                         Double_t phiMin=0.0,
                         Double_t phiMax=TMath::TwoPi(),
                         Double_t etaMin=-1.0,
                         Double_t etaMax= 1.0 );

  //data members
  TObjArray*              fTrackCollection;           //-> collection of tracks
  Int_t                   fReferenceMultiplicity;     // reference multiplicity
  Int_t                   fNumberOfTracks;            // number of tracks
  Int_t                   fNumberOfRPs;               // number of tracks that have passed the RP selection
  Int_t                   fNumberOfPOIs;              // number of tracks that have passed the POI selection
  Bool_t                  fUseGlauberMCSymmetryPlanes;// Use symmetry planes (Glauber MC)
  Bool_t                  fUseExternalSymmetryPlanes; // Use symmetry planes (external)
  Double_t                fPsi1;                      // Psi_1
  Double_t                fPsi2;                      // Psi_2
  Double_t                fPsi3;                      // Psi_3
  Double_t                fPsi4;                      // Psi_4
  Double_t                fPsi5;                      // Psi_5
  TF1*                    fPsi1Psi3;                  // Correlation between Psi_1 and Psi_3
  TF1*                    fPsi2Psi4;                  // Correlation between Psi_2 and Psi_4
  TF1*                    fPsi3Psi5;                  // Correlation between Psi_3 and Psi_5
  Double_t                fMCReactionPlaneAngle;      // the angle of the reaction plane from the MC truth
  Bool_t                  fMCReactionPlaneAngleIsSet; // did we set it from MC?
  Double_t                fAfterBurnerPrecision;      // iteration precision in afterburner
  Bool_t                  fUserModified;              // did we modify the event in any way (afterburner etc) ?
  TParameter<Int_t>*      fNumberOfTracksWrap;        //! number of tracks in TBrowser
  TParameter<Int_t>*      fNumberOfRPsWrap;           //! number of tracks that have passed the RP selection in TBrowser
  TParameter<Int_t>*      fNumberOfPOIsWrap;          //! number of tracks that have passed the POI selection in TBrowser
  TParameter<Double_t>*   fMCReactionPlaneAngleWrap;  //! the angle of the reaction plane from the MC truth in TBrowser
  Int_t*                  fShuffledIndexes;           //! placeholder for randomized indexes
  Bool_t                  fShuffleTracks;             // do we shuffle tracks on get?

  ClassDef(AliFlowEventSimple,1)
};

#endif


