// MUON tracking efficiency + (mother hadron) : (daughter muon) kinematic relation 
// Author : Saehanseul Oh


#ifndef AliMuonEffMC_h
#define AliMuonEffMC_h

class TH1F;
class TH1D;
class TH2F;
class TH3F;
class THn;
class TList;
class TObjArray;
class TObject;

class AliAODEvent;
class AliESDEvent;
class AliESDtrackCuts;
class AliEvtPoolManager;
class AliMCEvent;
class AliESDMuonTrack;
class AliAODTrack;
#include "AliAnalysisTaskSE.h"


class AliMuonEffMC : public AliAnalysisTaskSE {
 public:
  AliMuonEffMC();
  AliMuonEffMC(const char *name);
  virtual ~AliMuonEffMC();
 
  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);

  void         SetMcAna(Bool_t IsMc)               { fIsMc = IsMc;               }
  void         SetMDProcess(Bool_t MDProcess)      { fMDProcess = MDProcess;     }
  void         SetCentEstimator(TString Cent)      { fCentralityEstimator = Cent;}
  void         SetNEtaBins(Int_t NEtaBins)         { fNEtaBins = NEtaBins;       }
  void         SetNpTBins(Int_t NpTBins)           { fNpTBins = NpTBins;         }
  void         SetNCentBins(Int_t NCentBins)       { fNCentBins = NCentBins;     }
  void         SetNZvtxBins(Int_t NZvtxBins)       { fNZvtxBins = NZvtxBins;     }
  void         SetNPhiBins(Int_t NPhiBins)         { fNPhiBins = NPhiBins;       }

 protected:
  Bool_t       VertexOk(TObject* obj) const;
  Bool_t       IsGoodMUONtrack(AliESDMuonTrack &track);
  Bool_t       IsGoodMUONtrack(AliAODTrack &track);

 private:
  AliESDEvent *fESD;               //! ESD object
  AliAODEvent *fAOD;               //! AOD object 
  AliMCEvent  *fMC;                //! MC object
  Double_t     fCentrality;        //! Of current event
  Double_t     fZVertex;           //! Of current event
  TList       *fOutputList;        //! Output list
  TH1D        *fHEventStat;        //! statistics histo
  TH2F        *fHEvt;              //! Cent, vtx

  Bool_t       fIsMc;              //
  Bool_t       fMDProcess;         // (mother hadron) : (daughter muon) QA 

  TString      fCentralityEstimator;// 
  Int_t        fNEtaBins;          // number of eta bins
  Int_t        fNpTBins;           // number of p_T bins
  Int_t        fNCentBins;         // number of centrality bins
  Int_t        fNZvtxBins;         // number of Z-vertex bins
  Int_t        fNPhiBins;          // number of phi bins

  THn         *fHMuonParGen;       //! truth muon track eta, p_T, Centrality, Z-vertex, phi 
  THn         *fHMuonDetGen;       //! detector level muon track generated eta, p_T, Centrality, Z-vertex, phi 
  THn         *fHMuonDetRec;       //! reconstructed muon track eta, p_T, Centrality, Z-vertex, phi 
  THn         *fHEtcDetRec;        //! particle reconstructed at MUON detector, but not muon

  TH2F        *fHMuMotherGenPt[4]; //! particle-level muon p_T vs. mother p_T
  TH2F        *fHMuMotherRecPt[4]; //! detector-level muon p_T vs. mother p_T
  TH2F        *fHMuMotherGenPhi[4];//! particle-level muon phi vs. mother phi
  TH2F        *fHMuMotherRecPhi[4];//! detector-level muon phi vs. mother phi
  TH2F        *fHMuMotherGenEta[4];//! particle-level muon eta vs. mother eta
  TH2F        *fHMuMotherRecEta[4];//! detector-level muon eta vs. mother eta
  TH1F        *fHMuDCA[4];         //! muon DCA

  AliMuonEffMC(const AliMuonEffMC&);            // not implemented
  AliMuonEffMC &operator=(const AliMuonEffMC&); // not implemented

  ClassDef(AliMuonEffMC, 3);
};

#endif
