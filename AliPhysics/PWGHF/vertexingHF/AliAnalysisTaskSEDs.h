#ifndef ALIANALYSISTASKDS_H
#define ALIANALYSISTASKDS_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
///                                                                      //
/// \class Analysis task to produce Ds candidates mass spectra           //
/// \author Origin: F.Prino, Torino, prino@to.infn.it                    //
///                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliLog.h"
#include "AliExternalBDT.h"

class AliNormalizationCounter;


class AliAnalysisTaskSEDs : public AliAnalysisTaskSE
{
 public:

  enum {kpp, kPbPb, kUpgr};

  AliAnalysisTaskSEDs();
  AliAnalysisTaskSEDs(const char *name, AliRDHFCutsDstoKKpi* analysiscuts, Int_t fillNtuple=0);
  virtual ~AliAnalysisTaskSEDs();
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetWriteOnlySignalInNtuple(Bool_t opt=kTRUE){
    if(fReadMC) fWriteOnlySignal=opt;
    else AliError("fReadMC has to be kTRUE");
  }
  void SetFillNtuple(Int_t fill=0){fFillNtuple=fill;}
  void SetFillNSparse(Bool_t fill=kTRUE){fFillSparse=fill;}
  void SetFillNSparseDplus(Bool_t fill=kTRUE){fFillSparseDplus=fill;if(fill)fFillSparse=fill;}
  void SetFillNSparseImpPar(Bool_t fill=kTRUE){fFillImpParSparse=fill;}
  void SetFillNSparseAcceptanceLevel(Bool_t fill=kTRUE){fFillAcceptanceLevel=fill;}
  void SetMassRange(Double_t rang=0.4){fMassRange=rang;}
  void SetDoCutVarHistos(Bool_t opt=kTRUE) {fDoCutVarHistos=opt;}
  void SetUseSelectionBit(Bool_t opt=kFALSE){ fUseSelectionBit=opt;}
  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
  void SetUseRotBkg(Bool_t flag=kFALSE) {fDoRotBkg=flag;}
  void SetUseBkgFromPhiSB(Bool_t flag=kFALSE) {fDoBkgPhiSB=flag;}
  void SetPhiMassRange4RotBkg(Double_t range) {fMaxDeltaPhiMass4Rot=range;}
  void SetUseCutV0multVsTPCout(Bool_t flag) {fDoCutV0multTPCout=flag;}
  Bool_t CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau);
  Bool_t GetUseWeight() const {return fUseWeight;}
  void FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader, Double_t nTracklets);
  void GenerateRotBkg(AliAODRecoDecayHF3Prong *d, Int_t dec, Int_t iPtBin);
  void CreateCutVarsAndEffSparses();
  void CreateImpactParameterSparses();
  Float_t GetTrueImpactParameterDstoPhiPi(const AliAODMCHeader *mcHeader, TClonesArray* arrayMC, const AliAODMCParticle *partDs) const;

  void SetPtWeightsFromFONLL5anddataoverLHC16i2a();
  void SetPtWeightsFromFONLL5overLHC16i2abc();
  void SetPtWeightsFromFONLL5andBAMPSoverLHC16i2abc();
  void SetPtWeightsFromFONLL5andTAMUoverLHC16i2abc();
  void SetPtWeightsFromFONLL13overLHC17c3a12();
  void SetPtWeightsFromFONLL5overLHC18a4a2();

  void SetInvMassBinSize(Double_t binsiz=0.002){fMassBinSize=binsiz;}
  void SetPtBins(Int_t n, Float_t* lim);
  void SetAnalysisCuts(AliRDHFCutsDstoKKpi* cuts){fAnalysisCuts=cuts;}
  void SetSystem(Int_t system){fSystem = system;}

  void SetUseFinePtBinsForSparse(Bool_t usefinebins=kTRUE) {fUseFinPtBinsForSparse=usefinebins;} //use only in case of few candidates (e.g. MC signal only)

  Double_t GetPtWeightFromHistogram(Double_t pt);

  void SetKeepOnlyBkgFromHIJING(Bool_t keeponlyhijing=true) {fKeepOnlyBkgFromHIJING = keeponlyhijing;}
  void SetFillBkgSparse(Bool_t dofill=true) {fFillBkgSparse = dofill;}

  /// methods for ML application
  void SetDoMLApplication(Bool_t flag = kTRUE) {fApplyML = flag;}
  void SetMLConfigFile(TString path = ""){fConfigPath = path;}
  void EnablePIDMLHistos(Bool_t enable=kTRUE) {fEnablePIDMLHistos=enable;}
  void SetMLBinsForSparse(Int_t nbins = 300, Double_t min = 0.85, Double_t max = 1.) { fNMLBins = nbins; fMLOutputMin = min; fMLOutputMax = max;}

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*4;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*4+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*4+2;}
  Int_t GetReflSignalHistoIndex(Int_t iPtBin) const { return iPtBin*4+3;}
  /// methods for ML application
  Bool_t SetMLVariables(TString path);
  std::string GetFile(const std::string path);
  double CombineNsigmaDiffDet(double nsigmaTPC, double nsigmaTOF);

  enum {kMaxPtBins=36,knVarForSparse=14,knVarForSparseAcc=2,kVarForImpPar=3};

  AliAnalysisTaskSEDs(const AliAnalysisTaskSEDs &source);
  AliAnalysisTaskSEDs& operator=(const AliAnalysisTaskSEDs& source);

  TList*  fOutput;                        //!<! list send on output slot 0
  TH1F*   fHistNEvents;                   //!<! hist. for No. of events
  TH1F*   fHistoPtWeight;                 //-> user-defined histogram to calculate the Pt weights
  TH1F*   fChanHist[4];                   //!<! hist. with KKpi and piKK candidates (sig,bkg,tot)
  TH1F*   fMassHist[4*kMaxPtBins];        //!<! hist. of mass spectra (sig,bkg,tot)
  TH1F*   fMassHistPhi[4*kMaxPtBins];     //!<! hist. of mass spectra via phi (sig,bkg,tot)
  TH1F*   fMassHistK0st[4*kMaxPtBins];    //!<! hist. of mass spectra via K0* (sig,bkg,tot)
  TH1F*   fMassHistKK[kMaxPtBins];        //!<! hist. of mass spectra of KK
  TH1F*   fMassHistKpi[kMaxPtBins];       //!<! hist. of mass spectra of Kpi
  TH2F*   fMassHistKKVsKKpi[kMaxPtBins];  //!<! hist. of mass spectra of KK vs. mass spectra of KKpi
  TH2F*   fMassHistKpiVsKKpi[kMaxPtBins]; //!<! hist. of mass spectra of KK vs. mass spectra of KKpi
  TH1F*   fMassRotBkgHistPhi[kMaxPtBins]; //!<! hist. of bkg generated from rot. of the pion
  TH1F*   fMassLSBkgHistPhi[kMaxPtBins];  //!<! hist. of bkg generated from left phi sideband + pion
  TH1F*   fMassRSBkgHistPhi[kMaxPtBins];  //!<! hist. of bkg generated from right phi sideband + pion
  TH1F*   fCosPHist[4*kMaxPtBins];        //!<! hist. of cos pointing angle (sig,bkg,tot)
  TH1F*   fCosPxyHist[4*kMaxPtBins];      //!<! hist. of cosXY pointing angle (sig,bkg,tot)
  TH1F*   fDLenHist[4*kMaxPtBins];        //!<! hist. of decay length (sig,bkg,tot)
  TH1F*   fDLenxyHist[4*kMaxPtBins];      //!<! hist. of norm decay length XY (sig,bkg,tot)
  TH1F*   fNDLenxyHist[4*kMaxPtBins];     //!<! hist. of decay length XY (sig,bkg,tot)
  TH1F*   fSumd02Hist[4*kMaxPtBins];      //!<! hist. for sum d02 (Prod Cuts)
  TH1F*   fSigVertHist[4*kMaxPtBins];     //!<! hist. for sigVert (Prod Cuts)
  TH1F*   fPtMaxHist[4*kMaxPtBins];       //!<! hist. for Pt Max (Prod Cuts)
  TH1F*   fPtCandHist[4*kMaxPtBins];      //!<! hist. for Pt Max (Prod Cuts)
  TH1F*   fDCAHist[4*kMaxPtBins];         //!<! hist. for DCA (Prod Cuts)
  TH1F*   fNormIPHist[4*kMaxPtBins];      //!<! hist. for topomatic variable
  TH1F*   fCosPiDsHist[4*kMaxPtBins];     //!<! hist. for CosPiDs
  TH1F*   fCosPiKPhiHist[4*kMaxPtBins];   //!<! hist. for CosPiKPhi
  TH1F*   fPtProng0Hist[4*kMaxPtBins];    //!<! hist. for Pt Max (Prod Cuts)
  TH1F*   fPtProng1Hist[4*kMaxPtBins];    //!<! hist. for DCA (Prod Cuts)
  TH1F*   fPtProng2Hist[4*kMaxPtBins];    //!<! hist. for DCA (Prod Cuts)
  TH2F*   fDalitz[4*kMaxPtBins];          //!<! dalitz plot (sig,bkg,tot)
  TH2F*   fDalitzPhi[4*kMaxPtBins];       //!<! dalitz plot via phi (sig,bkg,tot)
  TH2F*   fDalitzK0st[4*kMaxPtBins];      //!<! dalitz plot via K0* (sig,bkg,tot)
  TH2F*   fPtVsMass;                      //!<! hist. of pt vs. mass (prod. cuts)
  TH2F*   fPtVsMassPhi;                   //!<! hist. of pt vs. mass (phi selection)
  TH2F*   fPtVsMassK0st;                  //!<! hist. of pt vs. mass (K0* selection)
  TH2F*   fYVsPt;                         //!<! hist. of Y vs. Pt (prod. cuts)
  TH2F*   fYVsPtSig;                      //!<! hist. of Y vs. Pt (MC, only sig, prod. cuts)
  TH2F*   fHistAllV0multNTPCout;          //!<! histo for V0mult vs #tracks TPCout (all)
  TH2F*   fHistSelV0multNTPCout;          //!<! histo for V0mult vs #tracks TPCout (sel)
  TH1F*   fHistCentrality[3];             //!<!hist. for cent distr (all,sel ev, )
  TH2F*   fHistCentralityMult[3];         //!<!hist. for cent distr vs mult (all,sel ev, )
  TH3F*   fCosPHist3D;                    //!<! cosP vs Ds mass vs pt
  TH3F*   fCosPxyHist3D;                  //!<! cosPxy vs Ds mass vs pt
  TH3F*   fDLenHist3D;                    //!<! Dlen vs Ds mass vs pt
  TH3F*   fDLenxyHist3D;                  //!<! Dlenxy vs Ds mass vs pt
  TH3F*   fNDLenxyHist3D;                 //!<! NDlenxy vs Ds mass vs pt
  TH3F*   fSigVertHist3D;                 //!<! SigVert vs Ds mass vs pt
  TH3F*   fDCAHist3D;                     //!<! DCA vs Ds mass vs pt
  TH3F*   fNormIPHist3D;                  //!<! nIP vs Ds mass vs pt
  TH3F*   fCosPiDsHist3D;                 //!<! cosPiDs vs Ds mass vs pt
  TH3F*   fCosPiKPhiHist3D;               //!<! cosPiKPhi vs Ds mass vs pt
  TH3F*   fPtProng0Hist3D;                //!<! Pt prong0 vs Ds mass vs pt
  TH3F*   fPtProng1Hist3D;                //!<! Pt prong1 vs Ds mass vs pt
  TH3F*   fPtProng2Hist3D;                //!<! Pt prong2 vs Ds mass vs pt
  TNtuple *fNtupleDs;                     //!<! output ntuple
  Int_t fFillNtuple;                      /// 0 not to fill ntuple
                                          /// 1 for filling ntuple for events through Phi
                                          /// 2 for filling ntuple for events through K0Star
                                          /// 3 for filling all

  Int_t   fSystem;                        /// 0 = pp, 1 = pPb,PbPb
  Bool_t  fReadMC;                        ///  flag for access to MC
  Bool_t  fWriteOnlySignal;               ///  flag to control ntuple writing in MC
  Bool_t  fDoCutVarHistos;                ///  flag to create and fill histos with distributions of cut variables
  Bool_t  fUseSelectionBit;               /// flag for usage of HasSelectionBit
  Bool_t  fFillSparse;                    /// flag for usage of THnSparse
  Bool_t  fFillSparseDplus;               /// flag for usage of THnSparse
  Bool_t  fFillImpParSparse;              /// flag for usage of sparse for imp. parameter
  Bool_t  fFillAcceptanceLevel;           /// flag for filling true reconstructed Ds at acceptance level (see FillMCGenAccHistos)
  Bool_t fDoRotBkg;                       ///flag to create rotational bkg (rotating pi track)
  Bool_t fDoBkgPhiSB;                     ///flag to create bkg from phi sidebands
  Bool_t fDoCutV0multTPCout;              ///flag to activate cut on V0mult vs #tracks TPCout
  Bool_t fUseWeight;                      /// flag to decide whether to use pt-weights != 1 when filling the container or not
  Int_t fAODProtection;                   /// flag to activate protection against AOD-dAOD mismatch.
                                          /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  UChar_t fNPtBins;                       /// number of Pt bins
  TList *fListCuts;                       /// list of cuts
  Float_t fPtLimits[kMaxPtBins+1];        ///  limits for pt bins
  Double_t fMassRange;                    /// range for mass histogram
  Double_t fMassBinSize;                  /// bin size for inv. mass histo
  Double_t fminMass;
  Double_t fmaxMass;
  Double_t fMaxDeltaPhiMass4Rot;          ///flag to set mass window of phi meson (when using pion rotation to create bkg)

  AliNormalizationCounter *fCounter;      //!<!Counter for normalization
  AliRDHFCutsDstoKKpi *fAnalysisCuts;     /// Cuts for Analysis

  THnSparseF *fnSparse;                   ///!<!THnSparse for candidates on data
  THnSparseF *fnSparseMC[5];              ///!<!THnSparse for MC
  ///[0]: Acc step prompt Ds
  ///[1]: Acc step FD Ds
  ///[2]: Selected prompt Ds
  ///[3]: Selected FD Ds
  ///[4]: Selected bkg candidates
  THnSparseF *fnSparseMCDplus[4];         ///!<!THnSparse for MC for D+->kkpi
  THnSparseF *fImpParSparse;              ///!<!THnSparse for imp. par. on data
  THnSparseF *fImpParSparseMC[4];         ///!<!THnSparse for imp. par. on MC

  TString fMultSelectionObjectName;       /// name of the AliMultSelection object to be considered
  Bool_t fUseFinPtBinsForSparse;          ///flag to fill pt axis of sparse with 0.1 GeV/c wide bins

  /// variables for ML application
  Bool_t fApplyML;                        /// flag to enable ML application
  TString fConfigPath;                    /// path to ML config file
  Int_t fNumVars;                         /// number of variables used in the model
  std::vector<std::string> fModelPaths;   /// vector of model paths
  std::vector<double> fModelOutputCuts;   /// vector of thresholds on model output
  std::vector<double> fPtBinsModel;       /// vector of pt bin lims
  std::vector<AliExternalBDT> fModels;    //!<! vector of ML models (BDTs for now)
  TH3F* fHistNsigmaPIDVsML[3][6];         //!<! histograms with PID Nsigma variables vs ML output
  Bool_t fEnablePIDMLHistos;              /// flag to enable control histograms for PID with ML
  Int_t fNMLBins;                         /// number of bins for ML output axis in THnSparse
  Double_t fMLOutputMin;                  /// min for ML output axis in THnSparse
  Double_t fMLOutputMax;                  /// max for ML output axis in THnSparse

  Bool_t fFillBkgSparse;                  /// flag to fill bkg sparse
  Bool_t fKeepOnlyBkgFromHIJING;          /// flag to keep the background from HIJING only

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEDs,34);    ///  AliAnalysisTaskSE for Ds mass spectra
  /// \endcond
};

#endif


