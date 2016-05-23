#ifndef AliAnalysisTaskElectronEfficiency_h
#define AliAnalysisTaskElectronEfficiency_h
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//###########################################################
//#                                                         #
//#             Single Electron Efficiency Task             #
//#        also, mainly for cross-checks:                   #
//#             Pair Efficiency Task                        #
//#             Pair-Prefilter Efficiency Task              #
//#                                                         #
//#  Authors:                                               #
//#   Patrick Reichelt, Uni Ffm / Patrick.Reichelt@cern.ch  #
//#   Theo Broeker, Uni Ffm / Theo.Broeker@cern.ch          #
//#   Markus Koehler, GSI / M.Koehler@gsi.de                #
//#                                                         #
//###########################################################
/**
 Fills 3D histograms (mcPt, mcEta, mcPhi) for generated and reconstructed electrons.
 Ratios of these for given acceptance regions give 3D track efficiencies, to be then used in a pair efficiency generator.
 Cut instances are defined by adding AliAnalysisFilters via a Config file. Therefore also an LMEECutLib is supported.
 Additional functionality:
 Determination of random electron rejection efficiency due to pair-prefiltering (used for photon conversion + Dalitz rejection).
 It is estimated by pairing primary, non-injected, charged pions with the selected electrons (so the pair has no real correlation)
 and applying the prefilter pair cuts to these random pairs. All and rejected pions are stored in 3D histograms.
 HOWEVER: this can be done in data as well. For this purpose the task "AliAnalysisTaskRandomRejection" is available now.
 ---
 As examples to set up the task, see PWGDQ/dielectron/macrosLMEE/Config_reichelt_ElectronEfficiency.C or Config_tbroeker_ElectronEfficiency.C
**/


#include "AliAnalysisTaskSE.h"

#include "THn.h"
#include <vector>
#include "AliAnalysisCuts.h"
#include "AliAnalysisFilter.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TTreeStream.h"//why?

class TString;
class TObject;
class TObjArray;
class AliPIDResponse;
class TH1F;
class TH2F;
class TH3F;
class TList;
class AliESDEvent;
class AliMCEvent;
class AliDielectronSignalMC;


class AliAnalysisTaskElectronEfficiency : public AliAnalysisTaskSE {
 public:
  
  AliAnalysisTaskElectronEfficiency();
  /// default constructor is mandatory for local LEGO train:
  /// W-TBufferFile::WriteObjectAny: since AliAnalysisTaskElectronEfficiency has no public constructor
  /// which can be called without argument, objects of this class can not be read with the current library.
  /// You will need to add a default constructor before attempting to read it.
  /// default constructor is mandatory to call in the AddTask: 'task->SetTriggerMask(triggerNames);'
  /// dont understand how this is inherited from: void AliEventTagCuts::SetTriggerMask(ULong64_t trmask)
  /// and it can crash, so instead, we implement the function and checks within this task, as done in 'AliAnalysisTaskMultiDielectron'.
  AliAnalysisTaskElectronEfficiency(const char *name); //const char *name = "AliAnalysisTaskElectronEfficiency"
  virtual ~AliAnalysisTaskElectronEfficiency();
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(const Option_t*);
  
  void          SetDoPairing(Bool_t b=kTRUE)                        {fDoPairing=b;}
  void          SetCalcResolution(Bool_t b=kTRUE)                   {fCalcResolution=b;}
  void          SetResolutionCuts(AliAnalysisFilter *cuts)          {fResolutionCuts=cuts;}
  void          SetKineTrackCuts(AliAnalysisFilter *cuts)           {fKineTrackCuts=cuts;}
  //void        SetPairCuts(AliAnalysisFilter *cuts)                {fPairCuts=cuts;}
  void          UsePhysicsSelection(Bool_t phy=kTRUE)               {fSelectPhysics=phy;}  // from AliAnalysisTaskMultiDielectron
  void          SetTriggerMask(ULong64_t mask)                      {fTriggerMask=mask;}   // from AliAnalysisTaskMultiDielectron
  void          SetEventFilter(AliAnalysisCuts * const filter)      {fEventFilter=filter;} // from Mahmuts AliAnalysisTaskSingleElectron
  void          SetCentralityRange(Double_t min, Double_t max)      {fCentMin=min; fCentMax=max;}
  void          SetCutInjectedSignal(Bool_t bCutInj)                { return; } //obsolete.
  void          SetNminEleInEventForRej(UInt_t nEle)                {fNminEleInEventForRej=nEle;}
  void          SetEtaRangeGEN(Double_t min, Double_t max)          {fEtaMinGEN=min; fEtaMaxGEN=max;}
  void          SetPtRangeGEN(Double_t min, Double_t max)           {fPtMinGEN=min; fPtMaxGEN=max;}
  void          SetSupportedCutInstance(Int_t supp)                 {fSupportedCutInstance=supp;}
  void          SetWriteTree(Bool_t write)                          {fWriteTree=write;}
  void          SetPIDResponse(AliPIDResponse *fPIDRespIn)          {fPIDResponse=fPIDRespIn;}
  void          SetRandomizeDaughters(Bool_t random=kTRUE)          {fRandomizeDaughters=random;}
  void          SetResolutionP(TObjArray *resArr)                   {fPResArr=resArr;}
  void          SetResolutionPt(TObjArray *resArr)                  {fPtResArr=resArr;}
  void          SetResolutionTheta(TObjArray *resArr)               {fThetaResArr=resArr;}
  void          SetResolutionEta(TObjArray *resArr)                 {fEtaResArr=resArr;}
  void          SetResolutionPhi(TObjArray *resArr)                 {fPhiResArr=resArr;}
  void          SetCalcEfficiencyRec(Bool_t b)                      {fCalcEfficiencyRec=b;}
  void          SetCalcEfficiencyPoslabel(Bool_t b)                 {fCalcEfficiencyPoslabel=b;}
  
  void          AddSignalMC(AliDielectronSignalMC* signal);   // use the functionality from AliDielectronSignalMC & AliDielectronMC to choose electron sources.
  void          SetCentroidCorrFunction(TObject *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void          SetWidthCorrFunction(TObject *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void          SetCentroidCorrFunctionITS(TObject *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void          SetWidthCorrFunctionITS(TObject *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  
  void          SetBins(Int_t Nptbins, Double_t *PtBins, Int_t Netabins, Double_t *EtaBins, Int_t Nphibins, Double_t *PhiBins, Int_t Nmeebins=0, Double_t *Meebins=0x0, Int_t Npteebins=0, Double_t *Pteebins=0x0) {
    /**/          fPtBins=PtBins;   fEtaBins=EtaBins;   fPhiBins=PhiBins;   fMeeBins=Meebins; fPteeBins=Pteebins;
    /**/          fNptBins=Nptbins; fNetaBins=Netabins; fNphiBins=Nphibins; fNmeeBins=Nmeebins; fNpteeBins=Npteebins;
    /**/        }
  void          SetRunBins(TString runs)                      { fsRunBins=runs; }
  void          AttachTrackCuts(AliAnalysisFilter *cuts)      { fvTrackCuts.push_back(cuts); }
  void          AttachExtraTrackCuts(AliAnalysisFilter *cuts) { fvExtraTrackCuts.push_back(cuts); }
  void          AttachDoPrefilterEff(Bool_t doPref)           { fvDoPrefilterEff.push_back(doPref); }
  void          AttachRejCutMee(Double_t rejcut)              { fvRejCutMee.push_back(rejcut); }
  void          AttachRejCutTheta(Double_t rejcut)            { fvRejCutTheta.push_back(rejcut); }
  void          AttachRejCutPhiV(Double_t rejcut)             { fvRejCutPhiV.push_back(rejcut); }
  void          SetPairCutMee(Double_t cut)                   { fPairCutMee=cut; }
  void          SetPairCutTheta(Double_t cut)                 { fPairCutTheta=cut; }
  void          SetPairCutPhiV(Double_t cut)                  { fPairCutPhiV=cut; }
  
  
  virtual void  CreateHistograms(TString names, Int_t cutInstance);
  void          CreateHistoGen();
  void          CreateSupportHistos();
  
  UInt_t        GetNCutsets() const { return fvReco_Ele.size(); } // 'cutset' and 'cutInstance' are used as synonymes in this task!
  //AliPIDResponse* GetPIDResponse() { return fPIDResponse; }
  
 private:
  void          CalcPrefilterEff(AliMCEvent* mcEventLocal, const std::vector< std::vector<Int_t> > & vvEleCand, const std::vector<Bool_t> & vbEleExtra);
  Double_t      PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2);
  const char*   GetParticleName(Int_t pdg) {
    /**/          TParticlePDG* p1 = TDatabasePDG::Instance()->GetParticle(pdg);
    /**/          if(p1) return p1->GetName();
    /**/          return Form("%d", pdg);
    /**/        }
  TVectorD*     GetPDGcodes();
  Double_t      GetSmearing(TObjArray *arr, Double_t x);
  
  AliESDEvent*      fESD;
  AliMCEvent*       mcEvent;
  AliPIDResponse*   fPIDResponse;
  TH1*              fPostPIDCntrdCorrTPC;     // post pid correction object for centroids in TPC
  TH1*              fPostPIDWdthCorrTPC;      // post pid correction object for widths in TPC
  TH1*              fPostPIDCntrdCorrITS;     // post pid correction object for centroids in ITS
  TH1*              fPostPIDWdthCorrITS;      // post pid correction object for widths in ITS
  TBits*            fUsedVars;                // used variables by AliDielectronVarManager
  TObjArray*        fSignalsMC;               // array of AliDielectronSignalMC

  Bool_t            fDoPairing;
  Bool_t            fSelectPhysics;           // Whether to use physics selection
  UInt_t            fTriggerMask;             // Event trigger mask
  AliAnalysisCuts*  fEventFilter;             // event filter
  Bool_t            fRequireVtx;
  UInt_t            fNminEleInEventForRej;    // should be fNminEleInEventForRej=2, because if there are less than 2 electrons, then no random ee-pair would be possible.
  Int_t             fSupportedCutInstance;    // for debugging
  //Int_t             fEventcount;
  Bool_t            fRandomizeDaughters;      // shuffle daughters at pair creation (sorted according to pt by default, which affects PhivPair at least for Like Sign)
  TRandom3          fRandom3;
  
  Double_t          fMaxVtxZ;
  Double_t          fCentMin;                 // should be fCentMin=-1 for pp and p-Pb
  Double_t          fCentMax;
  Double_t          fEtaMinGEN;
  Double_t          fEtaMaxGEN;
  Double_t          fPtMinGEN;
  Double_t          fPtMaxGEN;
  Int_t             fNptBins;
  Int_t             fNetaBins;
  Int_t             fNphiBins;
  Double_t*         fPtBins;                  //! ("!" to avoid streamer error)
  Double_t*         fEtaBins;                 //! ("!" to avoid streamer error)
  Double_t*         fPhiBins;                 //! ("!" to avoid streamer error)
  TString           fsRunBins;                // for run dependency histograms
  Bool_t            fCalcEfficiencyRec;
  Bool_t            fCalcEfficiencyPoslabel;
  
  
  //Cut Settings
  std::vector<AliAnalysisFilter*> fvTrackCuts;
  std::vector<AliAnalysisFilter*> fvExtraTrackCuts; // used in prefilter cutsets for global electron cuts to find relevant events for prefilter efficiency determination. // has to be a subset of 'fvTrackCuts', otherwise the treatment is incorrect!
  std::vector<Bool_t>             fvDoPrefilterEff;
  std::vector<Double_t>           fvRejCutMee;
  std::vector<Double_t>           fvRejCutTheta;
  std::vector<Double_t>           fvRejCutPhiV;
  //Efficiency Histograms
  TH3D*                           fNgen_Ele;
  std::vector<TH3D*>              fvReco_Ele;           // store reconstructed electrons (N vs pT, eta, phi) per cutset.
  std::vector<TH3D*>              fvReco_Ele_poslabel;  // store also result when using only tracks with positive label, for systematic checks.
  TH3D*                           fNgen_Pos;
  std::vector<TH3D*>              fvReco_Pos;           // store reconstructed positrons (N vs pT, eta, phi) per cutset.
  std::vector<TH3D*>              fvReco_Pos_poslabel;  // store also result when using only tracks with positive label, for systematic checks.
  std::vector<TH3D*>              fvAllPionsForRej;     // testparticles for prefilter efficiency determination.
  std::vector<TH3D*>              fvPionsRejByAllSigns;
  std::vector<TH3D*>              fvPionsRejByUnlike;
  //std::vector<TH3F*>             fvReco_Pio; // be really careful if you need to implement this (see comments in UserExec).
  //std::vector<TH3F*>             fvReco_Kao; // be really careful if you need to implement this (see comments in UserExec).
  //std::vector<TH3F*>             fvReco_Pro; // be really careful if you need to implement this (see comments in UserExec).
  
  // histograms with reconstructed observables
  TH3D*                           fNgen_Rec_Ele;
  TH3D*                           fNgen2_Rec_Ele;
  std::vector<TH3D*>              fvReco_Rec_Ele;           // store reconstructed electrons (N vs pT, eta, phi) per cutset.
  std::vector<TH3D*>              fvReco_Rec_Ele_poslabel;  // store also result when using only tracks with positive label, for systematic checks.
  TH3D*                           fNgen_Rec_Pos;
  TH3D*                           fNgen2_Rec_Pos;
  std::vector<TH3D*>              fvReco_Rec_Pos;           // store reconstructed positrons (N vs pT, eta, phi) per cutset.
  std::vector<TH3D*>              fvReco_Rec_Pos_poslabel;
  
  Int_t                           fNmeeBins;
  Int_t                           fNpteeBins;
  Double_t*                       fMeeBins;   //! ("!" to avoid streamer error)
  Double_t*                       fPteeBins;  //! ("!" to avoid streamer error)
  
  TH2D*                           fNgenPairs_sameMother;
  std::vector<TH2D*>              fvRecoPairs_sameMother;
  std::vector<TH2D*>              fvRecoPairs_poslabel_sameMother;
  TH2D*                           fNgenPairs_diffMothers;
  std::vector<TH2D*>              fvRecoPairs_diffMothers;
  std::vector<TH2D*>              fvRecoPairs_poslabel_diffMothers;
  
  
  // resolutions
  Bool_t                          fCalcResolution;
  
  TH1D*                           fDeltaPhiAll;
  TH1D*                           fDeltaPhi;
  TH2D*                           fDeltaPhi_alpha;
  TH2D*                           fDeltaPhi_pt;
  TH2D*                           fDeltaPhi_eta;
  TH2D*                           fDeltaPhi_MCcharge;
  TH2D*                           fDeltaPhi_charge;
  
  TH2D*                           fPGen_DeltaP;
  TH2D*                           fPtGen_DeltaPt;
  TH2D*                           fEtaGen_DeltaEta;
  TH2D*                           fThetaGen_DeltaTheta;
  TH2D*                           fPhiGen_DeltaPhi;
  TH2D*                           fOpeningAngleGen_DeltaOpeningAngleUS;
  TH2D*                           fOpeningAngleGen_DeltaOpeningAngleLS;
  
  TH2D*                           fPGen_DeltaP_pions;
  TH2D*                           fPtGen_DeltaPt_pions;
  TH2D*                           fEtaGen_DeltaEta_pions;
  TH2D*                           fThetaGen_DeltaTheta_pions;
  TH2D*                           fPhiGen_DeltaPhi_pions;
  TH2D*                           fOpeningAngleGen_DeltaOpeningAngleUS_pions;
  TH2D*                           fOpeningAngleGen_DeltaOpeningAngleLS_pions;  
  
  THnSparseF*                     fMgen_PtGen_mRes_ptRes;
  
  // external resolutions
  TObjArray*                      fPResArr; 
  TObjArray*                      fPtResArr; 
  TObjArray*                      fThetaResArr;
  TObjArray*                      fEtaResArr;
  TObjArray*                      fPhiResArr;
  TH1D*                           fSmearing_Ele[6];
  TH1D*                           fSmearing_Pos[6];
  
  AliAnalysisFilter*              fResolutionCuts;
  AliAnalysisFilter*              fKineTrackCuts;   // used for MC track acceptance cuts in pair efficiency calculation.
  Double_t                        fPairCutMee;      // used for pair cuts in pair efficiency calculation.
  Double_t                        fPairCutTheta;    // ''
  Double_t                        fPairCutPhiV;     // ''
  //AliAnalysisFilter*            fPairCuts;        // would be nicer, but cannot get AliDielectronPair working on MC :-(
  
  TList*                          fOutputList; // ! output data container
  TList*                          fOutputListSupportHistos; // ! output data container   
  TH1D*                           fEventStat;               // ! Histogram with event statistics
  
  //Output Tree with Tracks
  TTree*                          tracksT;
  Bool_t                          fWriteTree;
  
  //Track Variables
  // these are for data and MC:
  Float_t   pxESD;
  Float_t   pyESD;
  Float_t   pzESD;
  Float_t   pTPC; // determined via GetTPCInnerParam()->P() or GetInnerParam()->P(), see .cxx file.
  Int_t     chargeT;
  Float_t   signalITS;
  Float_t   signalTPC;
  Float_t   beta;
  Float_t   kchi2ITS;
  Int_t     kNclsITS;
  Float_t   kITSchi2Cl; // redundant
  Int_t     kNclsTPC;
  Float_t   kTPCchi2Cl;
  Int_t     kNclsTPCdEdx;
  Float_t   kNFclsTPCr;
  Float_t   kNFclsTPCfCross;
  Int_t     kNtrkltsTRD;
  Int_t     kNtrkltsTRDPID;
  Float_t   sigmaEleITS;
  Float_t   sigmaEleTPC;
  Float_t   sigmaEleTOF;
  //Float_t   sigmaEleTRD;
  Float_t   probEleTRD;
  Float_t   sigmaPioITS;
  Float_t   sigmaPioTPC;
  Float_t   sigmaPioTOF;
  Float_t   sigmaKaoITS;
  Float_t   sigmaKaoTPC;
  Float_t   sigmaProITS;
  Float_t   sigmaProTPC;
  Bool_t    isGlobalT;
  Bool_t    isGlobalSDD;
  // these are only for MC:
  Int_t     labelT;
  Int_t     pdgT;
  Int_t     labelmotherT;
  Int_t     pdgmotherT;
  Int_t     labelgrandmotherT;
  Int_t     pdggrandmotherT;
  //Float_t   pMC = 0;
  Float_t   pxMC;
  Float_t   pyMC;
  Float_t   pzMC;
  UInt_t    fSelectedByCut; // bit mask
  UInt_t    fSelectedByExtraCut; // bit mask
  
  //protected:
  enum {kAllEvents=0, kPhysicsSelectionEvents, kFilteredEvents , kEventStatBins};

  AliAnalysisTaskElectronEfficiency(const AliAnalysisTaskElectronEfficiency&); // not implemented
  AliAnalysisTaskElectronEfficiency& operator=(const AliAnalysisTaskElectronEfficiency&); // not implemented
  
  ClassDef(AliAnalysisTaskElectronEfficiency, 5);
};

#endif
