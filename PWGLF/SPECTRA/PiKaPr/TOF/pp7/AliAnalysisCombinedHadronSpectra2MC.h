#ifndef ALIANALYSISTASKCHARGEDHADRONSPECTRAMC_H
#define ALIANALYSISTASKCHARGEDHADRONSPECTRAMC_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// This analysis extracts pT-spectra of charged kaons, protons, and pions.  //
// It is based on particles identifation via the dE/dx signal of the TPC.   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TH1;
class TH1F;
class TH2F;
class TH3F;
class TList;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliHeader;
class AliESDpid;
#include "AliTOFT0v1.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include <TTree.h>

class AliAnalysisFilter;
class AliCFContainer;
class TDatabasePDG;

#include "AliAnalysisTask.h"
#include "AliESDVertex.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTOFT0v1.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"



#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class AliAnalysisCombinedHadronSpectra2MC : public AliAnalysisTaskSE {
 public:
  AliAnalysisCombinedHadronSpectra2MC(const char *name);
  AliAnalysisCombinedHadronSpectra2MC();
  virtual ~AliAnalysisCombinedHadronSpectra2MC() {}
  //
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //
  //
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  //void           SetAlephParameters(const Double_t * parameters){for(Int_t j=0;j<5;j++) fAlephParameters[j] = parameters[j]; Initialize();};
  Int_t           Mult();
  Int_t multiplicity;
  Int_t vert;
  //
  
 private:
  //
  //void  BinLogAxis(const THnSparse *h, Int_t axisNumber);
  
  //
  AliESDEvent *fESD;                  //! ESD object
  TList       *fListHist;             //! list for histograms
  //
  AliESDtrackCuts * fESDtrackCuts;    // basic cut variables
  AliESDpid       * fESDpid;          // basic TPC object for n-sigma cuts
  Bool_t        fMCtrue;              // flag if real data or MC is processed
  Double_t      fAlephParameters[5];  // Aleph Parameters for Bethe-Bloch
  //
 
  TList* TOFCheck;

  Bool_t calibrateESD;
  Bool_t correctTExp;
  Bool_t useT0TOF;
  Double_t timeResolution; 
  Bool_t tuneTOFMC;
  TTree *fTreeTrack;
  TTree *fTreeEv;
  Bool_t fLoadOCDB;
  Int_t frunOld;
  Int_t frun; 
  AliTOFcalib *tofCalib;
  AliTOFT0maker *t0maker;
  Int_t fMCtracks;       // n MC trk 
  Int_t fMCPrimaries;    // MC primaries 
  Double_t fT0TOF0;      // best t0
  Double_t fT0TOF1;      // sigma best t0 in ps
  Double_t fT0TOF2;      // t0 fill
  Double_t fT0TOF3;      // n TOF tracks
  Double_t fT0TOF4;      // TOF t0
  Double_t fT0TOF5;      // TOF t0 sigma
  Double_t fT0TOF6;      // sigma t0 fill
  Double_t fT0TOF7;      // n TOF tracks used for T0
  Double_t XPrimVertex;
  Double_t YPrimVertex;
  Double_t ZPrimVertex;
  Int_t NContrPrimVertex;
  Double_t rapidityMC;
  Float_t fDCAXY;
  Float_t fDCAZ;
  Int_t fcut;
  Int_t fTOFout;
  Int_t ftrdout;
  Int_t ftime;
  Int_t ftpcclust;
  Double_t flength;
  Int_t fsign;
  Double_t ftimetof;
  Int_t ftofchan;
  Double_t feta;
  Double_t fphi;
  Double_t fmomtrasv;
  Double_t sigmapi;
  Double_t sigmaka;
  Double_t sigmapr;
  Float_t fTot;
  Double_t r1[5];
  Double_t fmom;
  Double_t fexptimepi;
  Double_t fexptimeka;
  Double_t fexptimepr;
  Double_t ftofz; // local z  of track's impact on the TOF pad  
  Double_t ftofx;// local x  of track's impact on the TOF pad 
  Float_t t0track;
  Double_t TPCSignal;
  Float_t TPCSigmaPI;
  Float_t TPCSigmaKA;
  Float_t TPCSigmaPR;
  Int_t fmatch;
  Double_t fPhiout;
  Double_t fXout;
  Double_t fYout;
  Double_t fZout;
  Int_t  fTimeZeroType;      // flag to select timeZero type 
  
  Float_t spdCorr;
  Double_t treeMCP;
  Double_t treeMCPt;
  Double_t treeMCEta;
  Double_t treeMCPhi;
  Int_t treeMCPdg;
  Double_t treeMCPBis; 
  Double_t treeMCPtBis; 
  Double_t treeMCEtaBis; 
  Double_t treeMCPhiBis; 
  Int_t treeMCPdgBis; 
  Float_t t0trackSigma;
  Double_t fptMC;
  Double_t fphiMC;
  Double_t fetaMC;
  Int_t fPdgcode;


  TH2D* pad;
  TH1D* resx; TH1D * resz; TH1D * tofres; TH1D * tofresTOF; TH1D * tofresgood; 
  TH1F *hNumMatch; 
  TH1F* hNumMatchPos; TH1F*  hNumMatchNeg; TH1F*  hDenMatch; 

  TH1F*hNumMatchPip;  TH1F*hNumMatchPim; TH1F*hNumMatchKap;  TH1F*hNumMatchKam; TH1F*hNumMatchPrp;  TH1F*hNumMatchPrm; TH1F*hDenMatchPip;  TH1F*hDenMatchPim; TH1F*hDenMatchKap;  TH1F*hDenMatchKam; TH1F*hDenMatchPrp;  TH1F*hDenMatchPrm;

 
TH1F*  hDenMatchPos; TH1F*  hDenMatchNeg; TH1F*  hNumMatchEta; TH1F*  hNumMatchPosEta; TH1F*  hNumMatchNegEta; TH1F*  hDenMatchEta; TH1F*  hDenMatchPosEta; TH1F*   hDenMatchNegEta; TH1F*  hNumMatchphiOut; TH1F*  hNumMatchPosphiOut; TH1F*  hNumMatchNegphiOut; TH1F*  hDenMatchphiOut; TH1F*  hDenMatchPosphiOut; TH1F*  hDenMatchNegphiOut; TH1F*  hNumMatchEtaPtMa; TH1F*  hNumMatchPosEtaPtMa; TH1F*  hNumMatchNegEtaPtMa; TH1F*  hDenMatchEtaPtMa; TH1F*  hDenMatchPosEtaPtMa; TH1F*  hDenMatchNegEtaPtMa; TH1F*  hNumMatchphiOutPtMa; TH1F*  hNumMatchPosphiOutPtMa; TH1F*  hNumMatchNegphiOutPtMa; TH1F*  hDenMatchphiOutPtMa; TH1F*  hDenMatchPosphiOutPtMa; TH1F*  hDenMatchNegphiOutPtMa; TH1F* hNumMatchTRDOut; TH1F* hNumMatchPosTRDOut; TH1F* hNumMatchNegTRDOut; TH1F* hDenMatchTRDOut; TH1F* hDenMatchPosTRDOut; TH1F* hDenMatchNegTRDOut; TH1F* hNumMatchNoTRDOut; TH1F* hNumMatchPosNoTRDOut; TH1F* hNumMatchNegNoTRDOut; TH1F* hDenMatchNoTRDOut; TH1F* hDenMatchPosNoTRDOut; TH1F* hDenMatchNegNoTRDOut; TH1F* hNumMatchTPCpip; TH1F* hNumMatchTPCkap; TH1F* hNumMatchTPCprp; TH1F* hDenMatchTPCpip; TH1F* hDenMatchTPCkap; TH1F* hDenMatchTPCprp; TH1F* hNumMatchTPCpim; TH1F* hNumMatchTPCkam; TH1F* hNumMatchTPCprm; TH1F* hDenMatchTPCpim; TH1F* hDenMatchTPCkam; TH1F* hDenMatchTPCprm;

  TH1F *hNumEv;
  TH1F* hNumMatchMultTrkTRDOut[7][6];
  TH1F* hDenMatchMultTrkTRDOut[7][6];
  TH1F* hDenTrkMultTrkTRDOut[7][6];
  TH1F* hNumMatchMultTrkNoTRDOut[7][6];
  TH1F* hDenMatchMultTrkNoTRDOut[7][6];
  TH1F* hDenTrkMultTrkNoTRDOut[7][6];

  TH1F *hNumMatchMultTrk[7][6];
  TH1F *hDenMatchMultTrk[7][6]; 
  TH1F *hDenTrkMultTrk[7][6];
  TH1F *hNumMatchMultSPD[7][6];
  TH1F *hDenMatchMultSPD[7][6]; 
  TH1F *hDenTrkMultSPD[7][6];
  
  TH1F *hNumMatchMultTrkInc[7][2];
  TH1F *hDenMatchMultTrkInc[7][2];
  TH1F *hNumMatchMultSPDInc[7][2];
  TH1F *hDenMatchMultSPDInc[7][2];
  TH1F* hDenTrkVertMultTrk[6];
   
  TH1F* hDenTrkTriggerMultTrk[6];
 
  //
  AliAnalysisCombinedHadronSpectra2MC(const AliAnalysisCombinedHadronSpectra2MC&); 
  AliAnalysisCombinedHadronSpectra2MC & operator=(const AliAnalysisCombinedHadronSpectra2MC&); 

  ClassDef(AliAnalysisCombinedHadronSpectra2MC, 1); 
};

#endif