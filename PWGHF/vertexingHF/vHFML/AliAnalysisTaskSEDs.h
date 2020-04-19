#ifndef ALIANALYSISTASKDS_H
#define ALIANALYSISTASKDS_H

/* Copyright(c) 2007-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////
//                                                               //
// Analysis task to produce Ds candidates mass spectra           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
// Mantainers: F. Catalano, fabio.catalano@cern.ch               //
//             F. Grosa, fabrizio.grosa@cern.ch                  //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <THnSparse.h>
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"
#include "AliHFMLVarHandlerDstoKKpi.h"

class TH1F;
class TH2F;
class TH3F;
class TTree;
class TClonesArray;

class AliAODMCHeader;
class AliAODMCParticle;
class AliNormalizationCounter;
class AliAODRecoDecayHF3Prong;
class AliRDHFCutsDstoKKpi;
class AliHFMLResponseDstoKKpi;


class AliAnalysisTaskSEDs : public AliAnalysisTaskSE
{
 public:

  enum {kpp, kPbPb, kUpgr};

  AliAnalysisTaskSEDs();
  AliAnalysisTaskSEDs(const char *name, AliRDHFCutsDstoKKpi* analysiscuts, Bool_t createMLtree);
  virtual ~AliAnalysisTaskSEDs();
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetFillNSparse(Bool_t fill=kTRUE){fFillSparse=fill;}
  void SetFillNSparseDplus(Bool_t fill=kTRUE){fFillSparseDplus=fill;if(fill)fFillSparse=fill;}
  void SetFillNSparseMCAccWoQuark(Bool_t fill=kTRUE){fFillSparseAccWoQuark=fill;}
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
  void SetInvMassBinSize(Double_t binsiz=0.002){fMassBinSize=binsiz;}
  void SetPtBins(Int_t n, Float_t* lim);
  void SetAnalysisCuts(AliRDHFCutsDstoKKpi* cuts){fAnalysisCuts=cuts;}
  void SetSystem(Int_t system){fSystem = system;}
  void SetUseFinePtBinsForSparse(Bool_t usefinebins=kTRUE) {fUseFinPtBinsForSparse=usefinebins;} //use only in case of few candidates (e.g. MC signal only)
  void SetKeepOnlyBkgFromHIJING(Bool_t keeponlyhijing=kTRUE) {fKeepOnlyBkgFromHIJING = keeponlyhijing;}
  void SetFillBkgSparse(Bool_t dofill=kTRUE) {fFillBkgSparse = dofill;}
  /// methods for ML application
  void SetDoMLApplication(Bool_t flag = kTRUE) {fApplyML = flag;}
  void SetMLConfigFile(TString path = ""){fConfigPath = path;}
  void SetMLBinsForSparse(Int_t nbins = 300, Double_t min = 0.85, Double_t max = 1.) { fNMLBins = nbins; fMLOutputMin = min; fMLOutputMax = max;}
  void EnablePIDMLSparses(Bool_t enable=kTRUE) {fEnablePIDMLSparses=enable;}
  /// methods for ML tree creation
  void SetCreateMLTree(Bool_t flag = kTRUE) {fCreateMLtree = flag;}
  void SetMLTreePIDopt(int opt) {fPIDopt = opt;} // default AliHFMLVarHandlerDstoKKpi::kNsigmaDetAndCombPID
  void SetMLTreeAddTrackVar(Bool_t flag = kTRUE) {fAddSingleTrackVar = flag;}
  void SetFillOnlySignalInMLtree(Bool_t opt = kTRUE) {
    if(fReadMC) fFillOnlySignal = opt;
    else {
      if(opt)
        AliError("fReadMC has to be kTRUE");
    }
  }
  void EnableMLTreeEvtSampling(Float_t fractokeep, ULong_t seed) {
    fEnableEvtSampling = kTRUE;
    fFracEvtToKeep = fractokeep;
    fSeedSampling = seed;
  }
  void EnableMLTreeCandSampling(Float_t fractokeep, Float_t maxptsampling) {
    fEnableCandSampling = kTRUE;
    fFracCandToKeep = fractokeep;
    fMaxCandPtSampling = maxptsampling;
  }
  
  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:
  enum {kMaxPtBins=36,knVarForSparse=14,knVarForSparseAcc=2,knVarForSparseAccFD=3,kVarForImpPar=3,knVarPID=14,knVarPIDcomb=8};

  AliAnalysisTaskSEDs(const AliAnalysisTaskSEDs &source);
  AliAnalysisTaskSEDs& operator=(const AliAnalysisTaskSEDs& source);

  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*4;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*4+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*4+2;}
  Int_t GetReflSignalHistoIndex(Int_t iPtBin) const { return iPtBin*4+3;}
  Bool_t CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau);
  Float_t GetTrueImpactParameterDstoPhiPi(const AliAODMCHeader *mcHeader, TClonesArray* arrayMC, const AliAODMCParticle *partDs) const;
  void FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader);
  void GenerateRotBkg(AliAODRecoDecayHF3Prong *d, Int_t dec, Int_t iPtBin);
  void CreateCutVarsAndEffSparses();
  void CreateImpactParameterSparses();
  void CreatePIDMLSparses();

  TList*  fOutput = nullptr;                    //!<! list send on output slot 0
  TH1F*   fHistNEvents = nullptr;               //!<! hist. for No. of events
  TH1F*   fChanHist[4] = {};                    //!<! hist. with KKpi and piKK candidates (sig,bkg,tot)
  TH1F*   fMassHist[4*kMaxPtBins] = {};         //!<! hist. of mass spectra (sig,bkg,tot)
  TH1F*   fMassHistPhi[4*kMaxPtBins] = {};      //!<! hist. of mass spectra via phi (sig,bkg,tot)
  TH1F*   fMassHistK0st[4*kMaxPtBins] = {};     //!<! hist. of mass spectra via K0* (sig,bkg,tot)
  TH1F*   fMassHistKK[kMaxPtBins] = {};         //!<! hist. of mass spectra of KK
  TH1F*   fMassHistKpi[kMaxPtBins] = {};        //!<! hist. of mass spectra of Kpi
  TH2F*   fMassHistKKVsKKpi[kMaxPtBins] = {};   //!<! hist. of mass spectra of KK vs. mass spectra of KKpi
  TH2F*   fMassHistKpiVsKKpi[kMaxPtBins] = {};  //!<! hist. of mass spectra of KK vs. mass spectra of KKpi
  TH1F*   fMassRotBkgHistPhi[kMaxPtBins] = {};  //!<! hist. of bkg generated from rot. of the pion
  TH1F*   fMassLSBkgHistPhi[kMaxPtBins] = {};   //!<! hist. of bkg generated from left phi sideband + pion
  TH1F*   fMassRSBkgHistPhi[kMaxPtBins] = {};   //!<! hist. of bkg generated from right phi sideband + pion
  TH1F*   fCosPHist[4*kMaxPtBins] = {};         //!<! hist. of cos pointing angle (sig,bkg,tot)
  TH1F*   fCosPxyHist[4*kMaxPtBins] = {};       //!<! hist. of cosXY pointing angle (sig,bkg,tot)
  TH1F*   fDLenHist[4*kMaxPtBins] = {};         //!<! hist. of decay length (sig,bkg,tot)
  TH1F*   fDLenxyHist[4*kMaxPtBins] = {};       //!<! hist. of norm decay length XY (sig,bkg,tot)
  TH1F*   fNDLenxyHist[4*kMaxPtBins] = {};      //!<! hist. of decay length XY (sig,bkg,tot)
  TH1F*   fSumd02Hist[4*kMaxPtBins] = {};       //!<! hist. for sum d02 (Prod Cuts)
  TH1F*   fSigVertHist[4*kMaxPtBins] = {};      //!<! hist. for sigVert (Prod Cuts)
  TH1F*   fPtMaxHist[4*kMaxPtBins] = {};        //!<! hist. for Pt Max (Prod Cuts)
  TH1F*   fPtCandHist[4*kMaxPtBins] = {};       //!<! hist. for Pt Max (Prod Cuts)
  TH1F*   fDCAHist[4*kMaxPtBins] = {};          //!<! hist. for DCA (Prod Cuts)
  TH1F*   fNormIPHist[4*kMaxPtBins] = {};       //!<! hist. for topomatic variable
  TH1F*   fCosPiDsHist[4*kMaxPtBins] = {};      //!<! hist. for CosPiDs
  TH1F*   fCosPiKPhiHist[4*kMaxPtBins] = {};    //!<! hist. for CosPiKPhi
  TH1F*   fPtProng0Hist[4*kMaxPtBins] = {};     //!<! hist. for Pt Max (Prod Cuts)
  TH1F*   fPtProng1Hist[4*kMaxPtBins] = {};     //!<! hist. for DCA (Prod Cuts)
  TH1F*   fPtProng2Hist[4*kMaxPtBins] = {};     //!<! hist. for DCA (Prod Cuts)
  TH2F*   fDalitz[4*kMaxPtBins] = {};           //!<! dalitz plot (sig,bkg,tot)
  TH2F*   fDalitzPhi[4*kMaxPtBins] = {};        //!<! dalitz plot via phi (sig,bkg,tot)
  TH2F*   fDalitzK0st[4*kMaxPtBins] = {};       //!<! dalitz plot via K0* (sig,bkg,tot)
  TH2F*   fPtVsMass = nullptr;                  //!<! hist. of pt vs. mass (prod. cuts)
  TH2F*   fPtVsMassPhi = nullptr;               //!<! hist. of pt vs. mass (phi selection)
  TH2F*   fPtVsMassK0st = nullptr;              //!<! hist. of pt vs. mass (K0* selection)
  TH2F*   fYVsPt = nullptr;                     //!<! hist. of Y vs. Pt (prod. cuts)
  TH2F*   fYVsPtSig = nullptr;                  //!<! hist. of Y vs. Pt (MC, only sig, prod. cuts)
  TH2F*   fHistAllV0multNTPCout = nullptr;      //!<! histo for V0mult vs #tracks TPCout (all)
  TH2F*   fHistSelV0multNTPCout = nullptr;      //!<! histo for V0mult vs #tracks TPCout (sel)
  TH1F*   fHistCentrality[3] = {};              //!<! hist. for cent distr (all,sel ev, )
  TH2F*   fHistCentralityMult[3] = {};          //!<! hist. for cent distr vs mult (all,sel ev, )
  TH3F*   fCosPHist3D = nullptr;                //!<! cosP vs Ds mass vs pt
  TH3F*   fCosPxyHist3D = nullptr;              //!<! cosPxy vs Ds mass vs pt
  TH3F*   fDLenHist3D = nullptr;                //!<! Dlen vs Ds mass vs pt
  TH3F*   fDLenxyHist3D = nullptr;              //!<! Dlenxy vs Ds mass vs pt
  TH3F*   fNDLenxyHist3D = nullptr;             //!<! NDlenxy vs Ds mass vs pt
  TH3F*   fSigVertHist3D = nullptr;             //!<! SigVert vs Ds mass vs pt
  TH3F*   fDCAHist3D = nullptr;                 //!<! DCA vs Ds mass vs pt
  TH3F*   fNormIPHist3D = nullptr;              //!<! nIP vs Ds mass vs pt
  TH3F*   fCosPiDsHist3D = nullptr;             //!<! cosPiDs vs Ds mass vs pt
  TH3F*   fCosPiKPhiHist3D = nullptr;           //!<! cosPiKPhi vs Ds mass vs pt
  TH3F*   fPtProng0Hist3D = nullptr;            //!<! Pt prong0 vs Ds mass vs pt
  TH3F*   fPtProng1Hist3D = nullptr;            //!<! Pt prong1 vs Ds mass vs pt
  TH3F*   fPtProng2Hist3D = nullptr;            //!<! Pt prong2 vs Ds mass vs pt
  Int_t   fSystem = kpp;                        /// 0 = pp, 1 = pPb,PbPb
  Bool_t  fReadMC = kFALSE;                     /// flag for access to MC
  Bool_t  fDoCutVarHistos = kTRUE;              /// flag to create and fill histos with distributions of cut variables
  Bool_t  fUseSelectionBit = kFALSE;            /// flag for usage of HasSelectionBit
  Bool_t  fFillSparse = kTRUE;                  /// flag for usage of THnSparse
  Bool_t  fFillSparseDplus = kFALSE;            /// flag for usage of THnSparse
  Bool_t  fFillImpParSparse = kFALSE;           /// flag for usage of THnSparse for imp. parameter
  Bool_t  fFillSparseAccWoQuark = kFALSE;       /// flag for usage of THnSparse for generated Ds w/o quark
  Bool_t  fFillAcceptanceLevel = kTRUE;         /// flag for filling true reconstructed Ds at acceptance level (see FillMCGenAccHistos)
  Bool_t  fDoRotBkg = kFALSE;                   /// flag to create rotational bkg (rotating pi track)
  Bool_t  fDoBkgPhiSB = kFALSE;                 /// flag to create bkg from phi sidebands
  Bool_t  fDoCutV0multTPCout = kFALSE;          /// flag to activate cut on V0mult vs #tracks TPCout
  Int_t   fAODProtection = 1;                   /// flag to activate protection against AOD-dAOD mismatch.
                                                /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  UChar_t fNPtBins = 0;                         /// number of Pt bins
  TList *fListCuts = nullptr;                   /// list of cuts
  Float_t fPtLimits[kMaxPtBins+1] = {};         /// limits for pt bins
  Double_t fMassRange = 0.8;                    /// range for mass histogram
  Double_t fMassBinSize = 0.002;                /// bin size for inv. mass histo
  Double_t fminMass = 1.6;
  Double_t fmaxMass = 2.5;
  Double_t fMaxDeltaPhiMass4Rot = 0.010;        /// flag to set mass window of phi meson (when using pion rotation to create bkg)
  AliNormalizationCounter *fCounter = nullptr;  //!<! Counter for normalization
  AliRDHFCutsDstoKKpi *fAnalysisCuts = nullptr; /// Cuts for Analysis
  /// variables for Sparse
  THnSparseF *fnSparse = nullptr;               //!<! THnSparse for candidates on data
  THnSparseF *fnSparseMC[7] = {};               //!<! THnSparse for MC
                                                ///[0]: Acc step prompt Ds
                                                ///[1]: Acc step FD Ds
                                                ///[2]: Selected prompt Ds
                                                ///[3]: Selected FD Ds
                                                ///[4]: Selected bkg candidates
                                                ///[5]: Acc step prompt Ds w/o quark
                                                ///[6]: Acc step FD Ds w/o quark
  THnSparseF *fnSparseMCDplus[4] = {};          //!<! THnSparse for MC for D+->kkpi
  THnSparseF *fImpParSparse = nullptr;          //!<! THnSparse for imp. par. on data
  THnSparseF *fImpParSparseMC[4] = {};          //!<! THnSparse for imp. par. on MC
  Bool_t fUseFinPtBinsForSparse = kFALSE;       ///flag to fill pt axis of sparse with 0.1 GeV/c wide bins
  Bool_t fFillBkgSparse = kFALSE;               /// flag to fill bkg sparse
  Bool_t fKeepOnlyBkgFromHIJING = kFALSE;       /// flag to keep the background from HIJING only
  /// variables for ML application
  Bool_t fApplyML = kFALSE;                     /// flag to enable ML application
  TString fConfigPath = "";                     /// path to ML config file
  AliHFMLResponseDstoKKpi* fMLResponse = nullptr;  //!<! object to handle ML response
  Int_t fNMLBins = 300;                         /// number of bins for ML output axis in THnSparse
  Double_t fMLOutputMin = 0.85;                 /// min for ML output axis in THnSparse
  Double_t fMLOutputMax = 1.0;                  /// max for ML output axis in THnSparse
  Bool_t fEnablePIDMLSparses = kFALSE;          /// flag to enable control histograms for PID with ML
  THnSparseF* fnSparseNsigmaPIDVsML[2] = {};    //!<! histograms with PID Nsigma variables vs ML output
  /// variables for ML tree creation
  Bool_t fCreateMLtree = kFALSE;
  AliHFMLVarHandlerDstoKKpi* fMLhandler = nullptr;  //!<! object to handle ML tree creation and filling
  TTree* fMLtree = nullptr;                     //!<! tree with candidates for ML
  int fPIDopt = AliHFMLVarHandlerDstoKKpi::kNsigmaDetAndCombPID; /// option for PID variables
  Bool_t fAddSingleTrackVar = kFALSE;           /// option to store single track variables
  Bool_t fFillOnlySignal = kFALSE;              /// option to store only signal when using MC
  Bool_t fEnableEvtSampling = kFALSE;           /// flag to apply event sampling
  Float_t fFracEvtToKeep = 1.1;                 /// fraction of events to be kept by event sampling
  ULong_t fSeedSampling = 0;                    /// seed for sampling
  Bool_t fEnableCandSampling = kFALSE;          /// flag to apply candidate sampling
  Float_t fFracCandToKeep = 1.1;                /// fraction of candidates to be kept by sampling
  Float_t fMaxCandPtSampling = 0.;              /// maximun candidate pt to apply sampling

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEDs, 40);       /// AliAnalysisTaskSE for Ds mass spectra
  /// \endcond
};

#endif


