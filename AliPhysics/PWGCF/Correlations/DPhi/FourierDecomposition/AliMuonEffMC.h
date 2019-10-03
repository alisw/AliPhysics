// Author : Saehanseul Oh

#ifndef AliMuonEffMC_h
#define AliMuonEffMC_h

#include "AliAnalysisTaskSE.h"
#include "AliPool.h"
#include "AliMuonTrackCuts.h"

class TH1F;
class TH1D;
class TH2F;
class TH3F;
class THn;
class TList;
class TObjArray;
class TObject;
class AliStack;
class AliAODEvent;
class AliESDEvent;
class AliESDtrackCuts;
class AliEvtPoolManager;
class AliMCEvent;
class AliESDMuonTrack;
class AliAODTrack;
class AliMCParticle;

class AliMuonEffMC : public AliAnalysisTaskSE {
 public:
  AliMuonEffMC();
  AliMuonEffMC(const char *name);
  virtual ~AliMuonEffMC();
 
  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);

  void         SetMcAna(Bool_t IsMc)           { fIsMc = IsMc;               }
  void         SetPlotMode(Int_t kPlotMode)    { fPlotMode = kPlotMode;      }
  void         SetCentEstimator(TString Cent)  { fCentralityEstimator = Cent;}
  void         SetNEtaBins(Int_t NEtaBins)     { fNEtaBins = NEtaBins;       }
  void         SetNpTBins(Int_t NpTBins)       { fNpTBins = NpTBins;         }
  void         SetNCentBins(Int_t NCentBins)   { fNCentBins = NCentBins;     }
  void         SetNZvtxBins(Int_t NZvtxBins)   { fNZvtxBins = NZvtxBins;     }
  void         SetNPhiBins(Int_t NPhiBins)     { fNPhiBins = NPhiBins;       }
 
  Double_t     deltaphi(Double_t phi);
  Int_t        GetFirstMother(Int_t muonlabel);
  Int_t        GetFirstPPMother(Int_t muonlabel);
  Double_t     GetSpecies(Int_t PdgCode);
  void         SetMuonCutMask(UInt_t mask = AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca | AliMuonTrackCuts::kMuMatchApt ) { fMuonCutMask = mask; }

 protected:
  Bool_t       VertexOk(TObject* obj) const;
  Bool_t       IsGoodMUONtrack(AliESDMuonTrack &track);
  Bool_t       IsGoodMUONtrack(AliAODTrack &track);
  Double_t     GetMuonTrackType(AliMCParticle &track);

 private:
  AliESDEvent *fESD;               //! ESD object
  AliAODEvent *fAOD;               //! AOD object
  AliMCEvent  *fMC;                //! MC object
  AliStack    *fStack;             //! MC stack
  Double_t     fCentrality;        //! Of current event
  Double_t     fZVertex;           //! Of current event
  TList       *fOutputList;        //! Output list
  TH1D        *fHEventStat;        //! statistics histo
  TH2F        *fHEvt;              //! Cent, vtx

  Bool_t       fIsMc;              //
  Int_t        fPlotMode;          //
  UInt_t       fMuonCutMask;       //  muon cut mask for above
  AliMuonTrackCuts *fMuonTrackCuts;//  muon track cut object
  TString      fCentralityEstimator;//
  Int_t        fNEtaBins;          // number of eta bins
  Int_t        fNpTBins;           // number of p_T bins
  Int_t        fNCentBins;         // number of centrality bins
  Int_t        fNZvtxBins;         // number of Z-vertex bins
  Int_t        fNPhiBins;          // number of phi bins

  THn         *fHFM;               //! first mother 
  THn         *fHPP;               //! physical primary mother
  THn         *fHMuFM;             //! first mother of MUON track
  THn         *fHMuFMrec;          //! reconstructed MUON track with first mother species 
  THn         *fHMuPP;             //! physical primary mother of MUON track 
  THn         *fHMuPPrec;          //! reconstructed MUON track with physical primary mother species
  THn         *fHMuRec;            //! reconstructed MUON track
  
  AliMuonEffMC(const AliMuonEffMC&);            // not implemented
  AliMuonEffMC &operator=(const AliMuonEffMC&); // not implemented

  ClassDef(AliMuonEffMC, 15);
};

#endif
