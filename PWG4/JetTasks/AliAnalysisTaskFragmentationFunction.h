#ifndef ALIANALYSISTASKFRAGMENTATIONFUNCTION_H
#define ALIANALYSISTASKFRAGMENTATIONFUNCTION_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//
// Task for fragmentation
// 

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

  Double_t      fR;                     // tmp
  Double_t      fdRdNdxi;               // tmp
  Double_t      fPartPtCut;             // tmp
  Double_t      fEfactor;               // tmp
  Int_t         fNff;                   // tmp
  Int_t         fNim;                   // tmp
  TList*        fList;                  // tmp

  Int_t         fGlobVar;               // tmp
 
  Bool_t        fCDFCut;                // tmp

  // INTERVALS
  Int_t     fnEBin;       // Number of energy bins
  Double_t  fEmin; // tmp
  Double_t  fEmax; // tmp
  Int_t     fnEInterval; // tmp

  Int_t     fnRBin;       // Number of radius bins
  Double_t  fRmin; // tmp
  Double_t  fRmax; // tmp
  Int_t     fnRInterval; // tmp
 
  // HISTOGRAMS LIMITS
  Int_t     fnEtaHBin; // tmp
  Double_t  fEtaHBinMin; // tmp
  Double_t  fEtaHBinMax; // tmp
 
  Int_t     fnPhiHBin; // tmp
  Double_t  fPhiHBinMin; // tmp
  Double_t  fPhiHBinMax; // tmp 

  Int_t     fnPtHBin; // tmp
  Double_t  fPtHBinMin; // tmp
  Double_t  fPtHBinMax; // tmp

  Int_t     fnEHBin; // tmp
  Double_t  fEHBinMin; // tmp
  Double_t  fEHBinMax; // tmp 

  Int_t     fnXiHBin; // tmp
  Double_t  fXiHBinMax; // tmp
  Double_t  fXiHBinMin; // tmp

  Int_t     fnPthadHBin; // tmp
  Double_t  fPthadHBinMin; // tmp
  Double_t  fPthadHBinMax; // tmp

  Int_t     fnZHBin; // tmp
  Double_t  fZHBinMin; // tmp
  Double_t  fZHBinMax; // tmp

  Int_t     fnThetaHBin; // tmp 
  Double_t  fThetaHBinMax; // tmp

  Int_t     fnCosThetaHBin; // tmp
  Double_t  fcosThetaHBinMin; // tmp
  Double_t  fcosThetaHBinMax; // tmp

  Int_t     fnkTHBin; // tmp
  Double_t  fkTHBinMin; // tmp
  Double_t  fkTHBinMax; // tmp

  Int_t     fnRHBin; // tmp
  Double_t  fRHBinMin; // tmp
  Double_t  fRHBinMax; // tmp

  Int_t fnPtTrigBin; // tmp

  //HISTOGRAMS
  TH1F**        fEtaMonoJet1H; // tmp
  TH1F**        fPhiMonoJet1H; // tmp 
  TH1F**        fPtMonoJet1H; // tmp
  TH1F**        fEMonoJet1H; // tmp

  TH1F***        fdNdXiMonoJet1H; // tmp
  TH1F***        fdNdPtMonoJet1H; // tmp
  TH1F***        fdNdZMonoJet1H; // tmp
  TH1F***        fdNdThetaMonoJet1H; // tmp
  TH1F***        fdNdcosThetaMonoJet1H; // tmp
  TH1F***        fdNdkTMonoJet1H; // tmp
  TH1F***        fdNdpTvsZMonoJet1H; // tmp
  TH1F***        fShapeMonoJet1H; // tmp
  TH1F***        fNMonoJet1sH; // tmp

  TH2F***        fThetaPtPartMonoJet1H; // tmp
  TH2F***        fcosThetaPtPartMonoJet1H; // tmp
  TH2F***        fkTPtPartMonoJet1H; // tmp
  TH2F***        fThetaPtJetMonoJet1H; // tmp
  TH2F***        fcosThetaPtJetMonoJet1H; // tmp
  TH2F***        fkTPtJetMonoJet1H; // tmp
  TH2F***        fpTPtJetMonoJet1H; // tmp

  //ARRAYS
  Double_t*        farrayEmin; //! tmp
  Double_t*        farrayEmax; //! tmp
  Double_t*        farrayRadii; //!  tmp
  Double_t*        farrayPtTrigmin; //! tmp
  Double_t*        farrayPtTrigmax; //! tmp

  // TRACK CONTROL PLOTS
  TH1F* fptAllTracks; //! tmp
  TH1F* fetaAllTracks; //! tmp
  TH1F* fphiAllTracks; //! tmp
  TH2F* fetaphiptAllTracks; //! tmp
  TH2F* fetaphiAllTracks; //! tmp
  TH1F* fptAllTracksCut; //! tmp
  TH1F* fetaAllTracksCut; //! tmp 
  TH1F* fphiAllTracksCut; //! tmp
  TH2F* fetaphiptAllTracksCut; //! tmp
  TH2F* fetaphiAllTracksCut; //! tmp

  TH1F** fptTracks; //! tmp
  TH1F** fetaTracks; //! tmp
  TH1F** fphiTracks; //!  tmp
  TH1F** fdetaTracks; //! tmp
  TH1F** fdphiTracks; //! tmp 
  TH2F** fetaphiptTracks; //! tmp 
  TH2F** fetaphiTracks; //! tmp
  TH2F** fdetadphiTracks; //! tmp 
  TH1F** fptTracksCut; //! tmp
  TH1F** fetaTracksCut; //! tmp
  TH1F** fphiTracksCut; //! tmp
  TH1F** fdetaTracksCut; //! tmp
  TH1F** fdphiTracksCut; //! tmp
  TH2F** fetaphiptTracksCut; //! tmp
  TH2F** fetaphiTracksCut; //! tmp
  TH2F** fdetadphiTracksCut; //! tmp
  TH1F** fNPtTrig; // tmp
  TH1F** fNPtTrigCut; // tmp

  TH2F* fvertexXY; //! tmp
  TH1F* fvertexZ; //! tmp
  TH1F* fEvtMult; //! tmp
  TH2F* fEvtMultvsJetPt; //! tmp
  TH2F* fPtvsEtaJet; //! tmp
  TH2F* fNpvsEtaJet; //! tmp
  TH2F* fNpevtvsEtaJet; //! tmp
  TH2F* fPtvsPtJet; //! tmp
  TH2F* fNpvsPtJet; //! tmp
  TH2F* fNpevtvsPtJet; //! tmp
  TH1F* fPtvsPtJet1D; //! tmp
  TH1F* fNpvsPtJet1D; //! tmp
  TH1F* fNpevtvsPtJet1D; //! tmp
  TH1F* fptLeadingJet; //! tmp
  TH1F* fetaLeadingJet; //! tmp
  TH1F* fphiLeadingJet; //! // 
  TH1F* fptJet; //! // tmp
  TH1F* fetaJet; //! // tmp
  TH1F* fphiJet; //! // tmp


  TList*        fHistList; //! Output list

  Int_t fNBadRuns; //!
  TH1F* fNBadRunsH; //!

  ClassDef(AliAnalysisTaskFragmentationFunction, 1) // Analysis task for standard jet analysis
};
 
#endif
