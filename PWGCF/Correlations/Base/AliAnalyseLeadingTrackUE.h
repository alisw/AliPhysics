//-*- Mode: C++ -*-
#ifndef ALIANALYSELEADINGTRACKUE_H
#define ALIANALYSELEADINGTRACKUE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class  for transverse regions analysis w.r.t leading track
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---
#include <TObject.h> 

class AliAODEvent;
class AliAODInputHandler;
class AliESDEvent;
class AliESDtrackCuts;
class AliAODTrack;
class AliESDTrack;
class AliGenPythiaEventHeader;
class AliInputEventHandler;
class AliLog;
class AliMCEvent;
class AliStack;
class AliVParticle;
class TClonesArray;
class TObject;
class TROOT;
class TVector3;
class AliVTrack;
class AliHelperPID;
class TFormula;

class AliAnalyseLeadingTrackUE : public TObject {

 public: 

  AliAnalyseLeadingTrackUE();                                                   //constructor
  AliAnalyseLeadingTrackUE & operator = (const AliAnalyseLeadingTrackUE & g);   //assignment operator
  virtual ~AliAnalyseLeadingTrackUE();                                          //virtual destructor

  // Setters
  void  SetParticleSelectionCriteria(Int_t filterbit, Bool_t onlyhadrons, Double_t etacut, Double_t etacutMin=-1., Double_t ptMin = 0) { fFilterBit   = filterbit;
  											     fOnlyHadrons = onlyhadrons;
										             fTrackEtaCut = etacut;     
										             fTrackEtaCutMin = etacutMin;     
										             fTrackPtMin = ptMin;}

  void SetEventSelection(UInt_t bits) { fEventSelection = bits; }
  void  SetDebug(Int_t debug) { fDebug = debug; } 
  Bool_t         ApplyCuts(TObject* track);                       // Reproduces the cuts of the corresponding bit in the ESD->AOD filtering
  void           DefineESDCuts(Int_t filterbit);                                 // Emulate filterbit
  TObjArray*     FindLeadingObjects(TObject* obj);                                 // Looks for leading track or MC particle
  TObjArray*     GetMinMaxRegion(TList* transv1, TList* transv2);                  // Sorts the transverse regions in MIN and MAX
  Int_t          NParticles(TObject *obj);                                         // Counts tracks or MC particles
  AliVParticle*  ParticleWithCuts(TObject* obj, Int_t ipart, Bool_t onlyprimaries = kTRUE, Int_t particleSpecies = -1, Bool_t onlyCharged = kTRUE);                     // Returns track or MC particle at position "ipart" if passes selection criteria
  void  	 QSortTracks(TObjArray &a, Int_t first, Int_t last);               // Sort by pT an array of AliVParticles 
  TObjArray*     SortRegions(const AliVParticle* leading, TObject* obj, TObject* arrayMC, Bool_t onlyprimaries = kTRUE); // Assign particles to towards, away or transverse regions
  TObjArray*     GetAcceptedParticles(TObject* obj, TObject* arrayMC, Bool_t onlyprimaries = kTRUE, Int_t particleSpecies = -1, Bool_t useEtaPtCuts = kFALSE, Bool_t speciesOnTracks = kTRUE, Double_t evtPlane = -999., Bool_t onlyCharged = kTRUE, ULong64_t generatorIndexMask = 0ull); 
  TObjArray*     GetFakeParticles(TObject* obj, TObject* arrayMC, Bool_t onlyprimaries, Int_t particleSpecies, Bool_t useEtaPtCuts);
  Bool_t         TriggerSelection(const TObject* obj);                                   // Select good triggers with AliPhysicsSelection class
  Bool_t         VertexSelection(const TObject* obj, Int_t ntracks, Double_t zed);       // Vertex selection: see implementation
  void 		 RemoveInjectedSignals(TObjArray* tracks, TObject* arrayMC, Int_t maxLabel);
  void 		 RemoveWeakDecays(TObjArray* tracks, TObject* mcObj);
  AliHelperPID*  GetHelperPID() { return fHelperPID; }
  void 		 SetHelperPID(AliHelperPID* pid) { fHelperPID = pid; }
  void		 SetTrackStatus(UInt_t status) { fTrackStatus = status; }
  UInt_t	 GetTrackStatus() { return fTrackStatus; }
  void		 SetCheckMotherPDG(Bool_t checkpdg) { fCheckMotherPDG = checkpdg; }
  Bool_t	 GetCheckMotherPDG() { return fCheckMotherPDG; }
  void		 NextEvent() { fEventCounter++; }
  UInt_t	 GetEventCounter() { return fEventCounter; }
  void		 SetDCAXYCut(TFormula* value) { fDCAXYCut = value; }
  void		 SetDCAZCut(TFormula* value) { fDCAZCut = value; }
  void		 SetSharedClusterCut(Double_t value) { fSharedClusterCut = value; }
  void		 SetCrossedRowsCut(Int_t value) { fCrossedRowsCut = value; }
  void		 SetFoundFractionCut(Double_t value) { fFoundFractionCut = value; }
  void           SetParticlePhiCutEventPlane(Double_t min, Double_t max) { fTrackPhiCutEvPlMin = min; fTrackPhiCutEvPlMax = max; }

protected:
  Bool_t CheckTrack(AliVParticle * part);

private:
  AliAnalyseLeadingTrackUE(const AliAnalyseLeadingTrackUE & g);                 //copy constructor, not implemented
  Int_t          fDebug;             // debug flag
  Int_t          fFilterBit;         // track selection cuts
  UInt_t         fTrackStatus;       // if non-0, the bits set in this variable are required for each track
  Bool_t         fOnlyHadrons;       // consider only charged Pions, Protons and Kaons 
  Bool_t         fCheckMotherPDG;     // Check the PDG code of mother for secondaries 
  Double_t       fTrackEtaCut;       // pseudo-rapidity limit of transverse regions     
  Double_t       fTrackEtaCutMin;       // minimum of the pseudo-rapidity limit of transverse regions     
  Double_t       fTrackPhiCutEvPlMin;   // Minimum Phi cut on particles with respect to the Event Plane (values between 0 and Pi/2)
  Double_t       fTrackPhiCutEvPlMax;   // Maximum Phi cut on particles with respect to the Event Plane (values between 0 and Pi/2)
  Double_t       fTrackPtMin;        // pt limit for selecting particles
  UInt_t         fEventSelection;    // bit for physics selection
  TFormula*      fDCAXYCut;          // additional pt dependent cut on DCA XY (only for AOD)
  TFormula*      fDCAZCut;           // additional pt dependent cut on DCA Z (only for AOD)
  Double_t       fSharedClusterCut;  // cut on shared clusters (only for AOD)
  Int_t		 fCrossedRowsCut;   // cut on crossed rows (only for AOD)
  Double_t	 fFoundFractionCut;     // cut on crossed rows/findable clusters (only for AOD)
  
  AliESDtrackCuts *fEsdTrackCuts;    // set of cuts when reading ESD
  AliESDtrackCuts *fEsdTrackCutsExtra1;    // set of cuts when reading ESD
  AliESDtrackCuts *fEsdTrackCutsExtra2;    // set of cuts when reading ESD
  AliHelperPID   *fHelperPID;    // PID Helper object
  UInt_t fEventCounter;		  //! processed events (counter)
  
  ClassDef(AliAnalyseLeadingTrackUE,0)
};
#endif
