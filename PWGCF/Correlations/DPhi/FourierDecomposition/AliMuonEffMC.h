// MUON tracking efficiency + (mother hadron):(daughter muon) kinematic relation + x_F analysis
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
class AliStack;
class AliAODEvent;
class AliESDEvent;
class AliESDtrackCuts;
class AliEvtPoolManager;
class AliMCEvent;
class AliESDMuonTrack;
class AliAODTrack;
class AliMCParticle;

#include "AliAnalysisTaskSE.h"

class AliMuonEffMC : public AliAnalysisTaskSE {
 public:
  AliMuonEffMC();
  AliMuonEffMC(const char *name);
  virtual ~AliMuonEffMC();
 
  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);
  virtual Bool_t    Notify();
 
  void         SetMcAna(Bool_t IsMc)           { fIsMc = IsMc;               }
  void         SetIsPYTHIA(Bool_t IsPythia)    { fIsPythia = IsPythia;       }
  void         SetMDProcess(Bool_t kMDProcess) { fMDProcess = kMDProcess;    }
  void         SetPlotMode(Int_t kPlotMode)    { fPlotMode = kPlotMode;      }
  void         SetCentEstimator(TString Cent)  { fCentralityEstimator = Cent;}
  void         SetNEtaBins(Int_t NEtaBins)     { fNEtaBins = NEtaBins;       }
  void         SetNpTBins(Int_t NpTBins)       { fNpTBins = NpTBins;         }
  void         SetNCentBins(Int_t NCentBins)   { fNCentBins = NCentBins;     }
  void         SetNZvtxBins(Int_t NZvtxBins)   { fNZvtxBins = NZvtxBins;     }
  void         SetNPhiBins(Int_t NPhiBins)     { fNPhiBins = NPhiBins;       }

  void         MDProcess(Int_t isprimary, Int_t cutNum, Double_t trackpt, Double_t trackphi, Double_t tracketa, Double_t motherpt, Double_t motherphi, Double_t mothereta);
  Double_t     deltaphi(Double_t phi);
  Int_t        GetFirstPrimaryMother(Int_t muonlabel);
  Int_t        GetFirstPPMother(Int_t muonlabel);
  Double_t     GetSpecies(Int_t PdgCode);

 protected:
  Bool_t       VertexOk(TObject* obj) const;
  Int_t        GetMUONCutType(AliESDMuonTrack &track);
  Int_t        GetMUONCutType(AliAODTrack &track);
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
  TH1F        *fHXsec;             //! sum of cross section value for hard bins
  TH1F        *fHTrials;           //! N_trial
  TH2F        *fHEvt;              //! Cent, vtx

  Bool_t       fIsMc;              //
  Bool_t       fIsPythia;          //
  Bool_t       fMDProcess;         // (mother hadron) : (daughter muon) QA
  Int_t        fPlotMode;          //

  TString      fCentralityEstimator;//
  Int_t        fNEtaBins;          // number of eta bins
  Int_t        fNpTBins;           // number of p_T bins
  Int_t        fNCentBins;         // number of centrality bins
  Int_t        fNZvtxBins;         // number of Z-vertex bins
  Int_t        fNPhiBins;          // number of phi bins

  THn         *fHFPM;              //! first primary mother of all charged particles with centrality and z-vertex bin
  THn         *fHPP;               //! physical primary mothers with centrality and z-vertex bin
  THn         *fHDetRecMu[2];      //! reconstructed MUON track, detector level with centrality and z-vertex bin
  THn         *fHDetRecMuFPM[2];   //! reconstructed MUON track, detector level with FPM species
  THn         *fHDetRecMuPP[2];    //! reconstructed MUON track, detector level with PP species
  THn         *fHMuFPM[2];         //! reconstructed MUON track's first primary mother
  THn         *fHMuPP[2];          //! reconstructed MUON track's first physical primary mother

  TH2F        *fHMuMotherRecPt[5][3]; //! detector-level muon p_T vs. mother p_T
  TH2F        *fHMuMotherRecPhi[5][3];//! detector-level muon phi vs. mother phi
  TH2F        *fHMuMotherRecEta[5][3];//! detector-level muon eta vs. mother eta
  TH1F        *fHMuMohterPtDifRec[5][3][3]; //!
  TH1F        *fHMuMohterPhiDifRec[5][3][3]; //!
  TH1F        *fHMuMohterEtaDifRec[5][3][3]; //!

  TH2F        *fHZvRv[3];         //! muon decay Vertex Z-R for primary muon 
  TH2F        *fHXvYv[3];         //! muon decay Vertex x-y for primary muon 

  AliMuonEffMC(const AliMuonEffMC&);            // not implemented
  AliMuonEffMC &operator=(const AliMuonEffMC&); // not implemented

  ClassDef(AliMuonEffMC, 14);
};

#endif
