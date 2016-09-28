#ifndef AliAnaTaskV0EffDecomposition_H
#define AliAnaTaskV0EffDecomposition_H

#include <TList.h>
#include <TH1.h>
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 


/* /\* #include <TRandom.h> *\/ */
/* #include <TObject.h> */
/* #include <AliESDEvent.h> */
/* #include <AliMCEvent.h> */
/* #include <AliESDtrackCuts.h> */
/* #include "DebugClassesMultESA2013.h" */


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
  //  Double_t GetLowPtFraction()  { return fLowPtFraction; }
  Double_t GetMassCut()  { return fMassCut; }
  Double_t GetMinPtV0()  { return fMinPtV0; }
  Double_t GetDecayRmax()  { return fDecayRmax; }
  Double_t GetDecayRmin()  { return fDecayRmin; }
  Double_t GetDcaDaugh()  { return fDcaDaugh; }
  Double_t GetV0pndca()  { return fV0pndca; }
  Double_t GetCospt()  { return fCospt; }
  Double_t GetPdca()  { return fPdca; }
  Double_t GetNdca()  { return fNdca; }
  Double_t GetDcaV0()  { return fDcaV0; }
  Double_t GetCt()  { return fCt; }
  Double_t GetInvMass()  { return fInvMass; }
  Double_t GetNcl()  { return fNcl; }
  Double_t GetChi2perNDF()  { return fChi2perNDF; }


  virtual void  SetTrigger(UInt_t ktriggerInt) {fTrigBit = ktriggerInt;}
  virtual void  SetTrackFilterBit(UInt_t trackF) {fTrackFilterBit = trackF;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
  virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
  virtual void  SetPdgV0(Int_t pdg) {fPdgV0 = pdg;}
  virtual void  SetPdgPos(Int_t pdg) {fPdgPos = pdg;}
  virtual void  SetPdgNeg(Int_t pdg) {fPdgNeg = pdg;}

  virtual void  SetLowPtFraction(Double_t value) {fLowPtFraction = value;} 
  virtual void  SetMassCut(Double_t massCut){fMassCut = massCut;} 
  virtual void SetPdca(Double_t Pdca){fPdca = Pdca;}
  virtual void SetNdca(Double_t Ndca){fNdca = Ndca;}
  virtual void SetDcaV0(Double_t DcaV0){fDcaV0 = DcaV0;}
  virtual void  SetMinPtV0(Double_t value) {fMinPtV0 = value;}
  virtual void SetDecayRmax(Double_t DecayRmax){fDecayRmax = DecayRmax; }
  virtual void SetDecayRmin(Double_t DecayRmin){fDecayRmin = DecayRmin; }
  virtual void SetDcaDaugh(Double_t DcaDaugh){fDcaDaugh = DcaDaugh; }
  virtual void SetV0pndca(Double_t V0pndca){fV0pndca = V0pndca; }
  virtual void SetCospt(Double_t Cospt){fCospt = Cospt; }

   virtual void SetCt(Double_t Ct){fCt = Ct;}
   virtual void SetInvMass(Double_t InvMass){fInvMass = InvMass;}
   virtual void SetNcl(Double_t Ncl){fNcl = Ncl;}
   virtual void SetChi2perNDF(Double_t Chi2perNDF){fChi2perNDF = Chi2perNDF;}


  
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
  Double_t fLowPtFraction; 
  Double_t fMassCut; 
  Double_t fMinPtV0;
  Double_t fDecayRmax;
  Double_t fDecayRmin;
  Double_t fDcaDaugh;
  Double_t fV0pndca;
  Double_t fCospt;
  Double_t fPdca;
  Double_t fNdca;
  Double_t fDcaV0;
  Double_t fCt;
  Double_t fInvMass;
  Double_t fNcl;
  Double_t fChi2perNDF;


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

  TH1D*         hCentr;                //! 
  TH1D*         hCt;                //! 
  TH1D*         hDCAdaugh;          //! DCA between daughters
  TH1D*         hDecayR;            //! V0 decay radius
  TH1D*         hCosPA;             //! cosine of pointing angle
  TH1D*         hInvMass;           //!
  TH1D*         hNcl;    //!
  TH1D*         hChi2perNDF;  //!
  TH1D*         hPdca;              //!
  TH1D*         hNdca;              //!
  TH1D*         hDcaV0;              //!

  TH1D*         hCtHigh;                //! 
  TH1D*         hDCAdaughHigh;          //! DCA between daughters
  TH1D*         hDecayRHigh;            //! V0 decay radius
  TH1D*         hCosPAHigh;             //! cosine of pointing angle
  TH1D*         hInvMassHigh;           //!
  TH1D*         hNclHigh;    //!
  TH1D*         hChi2perNDFHigh;  //!
  TH1D*         hPdcaHigh;              //!
  TH1D*         hNdcaHigh;              //!
  TH1D*         hDcaV0High;              //!

  ClassDef(AliAnaTaskV0EffDecomposition, 1);    //Analysis task for v0 eff decomposition 
};

#endif
