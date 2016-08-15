#ifndef AliAnaTaskV0EffDecomposition_H
#define AliAnaTaskV0EffDecomposition_H

// ROOT includes
#include <TList.h>
#include <TH1.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 

class AliAnaTaskV0EffDecomposition : public AliAnalysisTaskSE {
 public:
  enum AnalysisMode { kInvalid = -1, kGlobalTrk = 0x1, kTPCTrk = 0x2 }; 
  AliAnaTaskV0EffDecomposition();
  AliAnaTaskV0EffDecomposition(const char *name);

  virtual ~AliAnaTaskV0EffDecomposition();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     

  virtual void  SetTrigger(UInt_t ktriggerInt) {fTrigBit = ktriggerInt;}
  virtual void  SetTrackFilterBit(UInt_t trackF) {fTrackFilterBit = trackF;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
  virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
  virtual void  SetPdgV0(Int_t pdg) {fPdgV0 = pdg;}
  virtual void  SetPdgPos(Int_t pdg) {fPdgPos = pdg;}
  virtual void  SetPdgNeg(Int_t pdg) {fPdgNeg = pdg;}
  
 private:
  virtual void ProcessMCTruthAOD();
  AliAODMCParticle* FindPrimaryMotherAOD(AliAODMCParticle* startParticle, Int_t& nSteps);
  virtual void AnalyzeV0AOD();
  AliAODMCParticle* ValidateTrack(AliAODTrack* track, 
				  Int_t pdgDaughter);
  virtual void AnalyzeDaughtersAOD();
  virtual void AnalyzeRecMothersAOD();

  AliAODEvent* fAOD;                  //! AOD object
  //  AliMCEvent*  fMC;               //! MC object
  TClonesArray* fMCArray;             //! MC array for AOD
  UInt_t fTrackFilterBit;             //  Track Filter, old cuts 2010

  //
  // Cuts and options
  //
  UInt_t       fTrigBit;
  Double_t     fVtxCut;    // Vtx cut on z position in cm
  Double_t     fEtaCut;    // Eta cut used to select particles
  Float_t      fMinCent;   // minimum centrality
  Float_t      fMaxCent;   // maximum centrality
  Int_t        fPdgV0;     // pdg of mother 
  Int_t        fPdgPos;    // pdg of pos daughter 
  Int_t        fPdgNeg;    // pdg of neg daughter 

  //
  // Output objects
  //
  TList*        fListOfObjects;     //! Output list of objects
  TH1D*         hV0Gen;             //! No of V0s generated
  TH1D*         hV0Rec;             //! No of V0s reconstructed
  TH1D*         hDaughterRec;       //! No of daughter pairs reconstructed
  TH1D*         hV0Ghost;           //! No of multi-rec V0s
  TH1D*         hTrackGhost;        //! No of multi-rec tracks
  TH1D*         hV0ButNoTracks;     //! No of V0s rec where daughters not found

  ClassDef(AliAnaTaskV0EffDecomposition, 1);    //Analysis task for v0 eff decomposition 
};

#endif