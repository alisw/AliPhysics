#ifndef ALIANALYSISTASKFRAGMENTATIONFUNCTION_H
#define ALIANALYSISTASKFRAGMENTATIONFUNCTION_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"

class AliJetHeader;
class AliAODJet;
class TProfile;
class TH1F;

class AliAnalysisTaskFragmentationFunction : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskFragmentationFunction();
    AliAnalysisTaskFragmentationFunction(const char* name);
    virtual ~AliAnalysisTaskFragmentationFunction() {;}
    // Implementation of interface methods
    virtual Bool_t Notify();
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    virtual void SetAODInput(Bool_t b){fUseAODInput = b;}
    virtual void SetLimitGenJetEta(Bool_t b){fLimitGenJetEta = b;}
    virtual void SetRecEtaWindow(Float_t f){fRecEtaWindow = f;}
    virtual void SetAnalysisType(Int_t i){fAnalysisType = i;}
    virtual void SetBranchGen(const char* c){fBranchGen = c;}
    virtual void SetBranchRec(const char* c){fBranchRec = c;}
    virtual void SetFilterMask(UInt_t i){fFilterMask = i;}

    virtual void FillMonoJetH(Int_t goodBin, AliAODJet* jet, TClonesArray* Tracks);
    virtual void DefineJetH();

//    virtual void DeleteHists();
    virtual void SetProperties(TH1* h,const char* x, const char* y);

    // Helper
    //

    // we have different cases
    // AOD reading -> MC from AOD
    // ESD reading -> MC from Kinematics
    // this has to match with our selection of input events
    enum {kTrackUndef = 0, kTrackAOD, kTrackKineAll,kTrackKineCharged, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance};
    enum {kAnaMC =  0x1, kAnaMCESD = 0x2};
    enum {kMaxJets = 4};
    enum {kMaxCorrelation =  3};
    
    // 
    // 0 all jets
    // 1 all jet in eta window
    // 2 all jets with partner
    // 3 all jets in eta window with partner
    // 4 all jets with partner in eta window
    enum {kStep0 = 0, kStep1, kStep2, kStep3, kStep4,kMaxStep};


 private:
  AliAnalysisTaskFragmentationFunction(const AliAnalysisTaskFragmentationFunction &det);
  AliAnalysisTaskFragmentationFunction &operator=(const AliAnalysisTaskFragmentationFunction &det);
  
  void MakeJetContainer();
 
 private:
  AliJetHeader *fJetHeaderRec;
  AliJetHeader *fJetHeaderGen;
  AliAODEvent  *fAOD; // wherewe take the jets from can be input or output AOD
  //  THnSparseF   *fhnJetContainer[kMaxStep*2];   // like particle container in corrfw with different steps need AliCFContainer with Scale(), and clone() to do the same
  //  THnSparseF   *fhnCorrelation;           // response matrix for unfolding 
 
  TString       fBranchRec;  // AOD branch name for reconstructed
  TString       fBranchGen;  // AOD brnach for genereated

  Bool_t        fUseAODInput;           // use AOD input
  Bool_t        fUseAODJetInput;        // use AOD input
  Bool_t        fUseAODTrackInput;      // take track from input AOD not from ouptu AOD
  Bool_t        fUseAODMCInput;         // take MC from input AOD not from ouptu AOD
  Bool_t        fUseGlobalSelection;    // Limit the eta of the generated jets
  Bool_t        fUseExternalWeightOnly; // use only external weight
  Bool_t        fLimitGenJetEta;        // Limit the eta of the generated jets
  UInt_t        fFilterMask;            // filter bit for slecected tracks
  Int_t         fAnalysisType;          // Analysis type 
  Int_t         fTrackTypeRec;          // type of tracks used for FF 
  Int_t         fTrackTypeGen;          // type of tracks used for FF 
  Float_t       fAvgTrials;             // Average nimber of trials
  Float_t       fExternalWeight;        // external weight
  Float_t       fRecEtaWindow;          // eta window used for corraltion plots between rec and gen 

  Double_t      fR;
  Double_t      fdRdNdxi;
  Double_t      fPartPtCut;
  Double_t      fEfactor;
  Int_t         fNff;
  Int_t         fNim;
  TList*        fList;

  Int_t         fGlobVar;

  Bool_t        fCDFCut;

  // INTERVALS
  Int_t     fnEBin;       // Number of energy bins
  Double_t  fEmin;
  Double_t  fEmax;
  Int_t     fnEInterval;

  Int_t     fnRBin;       // Number of radius bins
  Double_t  fRmin;
  Double_t  fRmax;
  Int_t     fnRInterval;

  // HISTOGRAMS LIMITS
  Int_t     fnEtaHBin;
  Double_t  fEtaHBinMin;
  Double_t  fEtaHBinMax;

  Int_t     fnPhiHBin;
  Double_t  fPhiHBinMin;
  Double_t  fPhiHBinMax;

  Int_t     fnPtHBin;
  Double_t  fPtHBinMin;
  Double_t  fPtHBinMax;

  Int_t     fnEHBin;
  Double_t  fEHBinMin;
  Double_t  fEHBinMax;

  Int_t     fnXiHBin;
  Double_t  fXiHBinMax;
  Double_t  fXiHBinMin;

  Int_t     fnPthadHBin;
  Double_t  fPthadHBinMin;
  Double_t  fPthadHBinMax;

  Int_t     fnZHBin;
  Double_t  fZHBinMin;
  Double_t  fZHBinMax;

  Int_t     fnThetaHBin;
  Double_t  fThetaHBinMin;
  Double_t  fThetaHBinMax;

  Int_t     fnCosThetaHBin;
  Double_t  fcosThetaHBinMin;
  Double_t  fcosThetaHBinMax;

  Int_t     fnkTHBin;
  Double_t  fkTHBinMin;
  Double_t  fkTHBinMax;

  Int_t     fnRHBin;
  Double_t  fRHBinMin;
  Double_t  fRHBinMax;

  Int_t fnPtTrigBin;

  //HISTOGRAMS
  TH1F**        fEtaMonoJet1H;
  TH1F**        fPhiMonoJet1H;
  TH1F**        fPtMonoJet1H;
  TH1F**        fEMonoJet1H;

  TH1F***        fdNdXiMonoJet1H;
  TH1F***        fdNdPtMonoJet1H;
  TH1F***        fdNdZMonoJet1H;
  TH1F***        fdNdThetaMonoJet1H;
  TH1F***        fdNdcosThetaMonoJet1H;
  TH1F***        fdNdkTMonoJet1H;
  TH1F***        fdNdpTvsZMonoJet1H;
  TH1F***        fShapeMonoJet1H;
  TH1F***        fNMonoJet1sH;

  TH2F***        fThetaPtPartMonoJet1H;
  TH2F***        fcosThetaPtPartMonoJet1H;
  TH2F***        fkTPtPartMonoJet1H;
  TH2F***        fThetaPtJetMonoJet1H;
  TH2F***        fcosThetaPtJetMonoJet1H;
  TH2F***        fkTPtJetMonoJet1H;
  TH2F***        fpTPtJetMonoJet1H;

  //ARRAYS
  Double_t*        farrayEmin; //!
  Double_t*        farrayEmax; //!
  Double_t*        farrayRadii; //!
  Double_t*        farrayPtTrigmin; //!
  Double_t*        farrayPtTrigmax; //!

  // TRACK CONTROL PLOTS
  TH1F* fptAllTracks; //!
  TH1F* fetaAllTracks; //!
  TH1F* fphiAllTracks; //!
  TH2F* fetaphiptAllTracks; //!
  TH2F* fetaphiAllTracks; //!
  TH1F* fptAllTracksCut; //!
  TH1F* fetaAllTracksCut; //!
  TH1F* fphiAllTracksCut; //!
  TH2F* fetaphiptAllTracksCut; //!
  TH2F* fetaphiAllTracksCut; //!

  TH1F** fptTracks; //!
  TH1F** fetaTracks; //!
  TH1F** fphiTracks; //!
  TH1F** fdetaTracks; //!
  TH1F** fdphiTracks; //!
  TH2F** fetaphiptTracks; //!
  TH2F** fetaphiTracks; //!
  TH2F** fdetadphiTracks; //!
  TH1F** fptTracksCut; //!
  TH1F** fetaTracksCut; //!
  TH1F** fphiTracksCut; //!
  TH1F** fdetaTracksCut; //!
  TH1F** fdphiTracksCut; //!
  TH2F** fetaphiptTracksCut; //!
  TH2F** fetaphiTracksCut; //!
  TH2F** fdetadphiTracksCut; //!
  TH1F** fNPtTrig;
  TH1F** fNPtTrigCut;

  TH2F* fvertexXY; //!
  TH1F* fvertexZ; //!
  TH1F* fEvtMult; //!
  TH2F* fEvtMultvsJetPt; //!
  TH2F* fPtvsEtaJet; //!
  TH2F* fNpvsEtaJet; //!
  TH2F* fNpevtvsEtaJet; //!
  TH2F* fPtvsPtJet; //!
  TH2F* fNpvsPtJet; //!
  TH2F* fNpevtvsPtJet; //!
  TH1F* fPtvsPtJet1D; //!
  TH1F* fNpvsPtJet1D; //!
  TH1F* fNpevtvsPtJet1D; //!
  TH1F* fptLeadingJet; //!
  TH1F* fetaLeadingJet; //!
  TH1F* fphiLeadingJet; //!
  TH1F* fptJet; //!
  TH1F* fetaJet; //!
  TH1F* fphiJet; //!


  TList*        fHistList; //! Output list

  Int_t fNBadRuns; //!
  TH1F* fNBadRunsH; //!

  ClassDef(AliAnalysisTaskFragmentationFunction, 1) // Analysis task for standard jet analysis
};
 
#endif
