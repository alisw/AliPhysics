#ifndef ALIANALYSISEBYERATIOS_H
#define ALIANALYSISEBYERATIOS_H

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//                        Analysis for Jet Hadrochemistry                       //
//                                                                              //
//    This analysis extracts pT-spectra of charged kaons, protons, and pions    //
//                      for the inclusive event and in jets.                    //
//   It is based on particles identification via the dE/dx signal of the TPC.   //
//                                                                              //
// Author: Sierra Weyhmiller <sierra.lisa.weyhmiller@cern.ch>, Yale University  //
//      Author: Mesut Arslandok <mesut.arslandok@cern.ch>, Yale University      //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

class TH1;
class THn;
class TH1F;
class TH2D;
class TH3D;
class TList;
class TTree;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliPIDResponse;
class AliHeader;
class AliESDpid;
//class AliESDtools;
class AliPIDCombined;
class AliJetContainer;
class AliAnalysisTaskRho;


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliPIDCombined.h"
#include "AliTPCdEdxInfo.h"
#include "AliESDv0KineCuts.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "THn.h"
#include "TVectorF.h"
#include "TCutG.h"
#include "TTreeStream.h"
#include "AliESDv0Cuts.h"
#include "AliEventCuts.h"

class AliAnalysisJetHadro : public AliAnalysisTaskEmcalJet {
public:

  AliEventCuts fEventCuts;     /// Event cuts

  // ---------------------------------------------------------------------------------
  //                           Constructor and Destructor
  // ---------------------------------------------------------------------------------

  AliAnalysisJetHadro(const char *name);
  AliAnalysisJetHadro();
  virtual ~AliAnalysisJetHadro();

  enum trackCutBit {
    kNCrossedRowsTPC60=0,
    kNCrossedRowsTPC80=1,
    kNCrossedRowsTPC100=2,
    kMaxChi2PerClusterTPCSmall=3,
    kMaxChi2PerClusterTPC=4,
    kMaxChi2PerClusterTPCLarge=5,
    kMaxDCAToVertexXYPtDepSmall=6,
    kMaxDCAToVertexXYPtDep=7,
    kMaxDCAToVertexXYPtDepLarge=8,
    kVertexZSmall=9,
    kVertexZ=10,
    kVertexZLarge=11,
    kEventVertexZSmall=12,
    kEventVertexZ=13,
    kEventVertexZLarge=14,
    kRequireITSRefit=15,
    kPixelRequirementITS=16,
    kNewITSCut=17,
    kActiveZoneSmall=18,
    kActiveZone=19,
    kActiveZoneLarge=20,
    kTPCSignalNSmall=21,
    kTPCSignalN=22,
    kTPCSignalNLarge=23,
    kCleanPrTOF=24,
    kCleanKaTOF=25,
    kCleanKaTOFTRD=26,
    kTrackProbKaTOF=27,
    kTrackProbPrTOF=28,
    kCleanDeTOF=29,
    kEventVertexZALICE=30,
    kEventVertexZALICETight=31,
  };

  enum centEst {
    kV0M=0,
    kCL0=1,
    kCL1=2,
  };

  enum kPDGpart{
    kPDGel=11,
    kPDGpi=211,
    kPDGka=321,
    kPDGpr=2212,
    kPDGde=1000010020,
    kPDGmu=13,
    kPDGla=3122,
  };

  /*
  kV0M=0,           // Centrality from V0A+V0C
  kCL0=1,           // Centrality from Clusters in layer 0
  kCL1=2,           // Centrality from Clusters in layer 1
  kTRK=3,           // Centrality from tracks
  kTKL=4,           // Centrality from tracklets
  kV0MvsFMD=5,      // Centrality from V0 vs FMD
  kTKLvsV0M=6,      // Centrality from tracklets vs V0
  kZEMvsZDC=7,      // Centrality from ZEM vs ZDC
  kV0A=8,           // Centrality from V0A
  kV0C=9,           // Centrality from V0C
  kZNA=10,          // Centrality from ZNA
  kZNC=11,          // Centrality from ZNC
  kZPA=12,          // Centrality from ZPA
  kZPC=13,          // Centrality from ZPC
  kCND=14,          // Centrality from tracks (candle condition)
  kFMD=15,          // Centrality from FMD
  kNPA=16,          // Centrality from Npart (MC)
  kV0A0=17,         // Centrality from V0A0
  kV0A123=18,       // Centrality from V0A123
  kV0A23=19,        // Centrality from V0A23
  kV0C01=20,        // Centrality from V0C01
  kV0S=21,          // Centrality from V0S
  kV0MEq=22,        // Centrality from V0A+V0C equalized channel
  kV0AEq=23,        // Centrality from V0A equalized channel
  kV0CEq=24,        // Centrality from V0C equalized channel
  kSPDClusters=25,  // Centrality from SPD Clusters
  kSPDTracklets=26, // Centrality from SPD Tracklets
  */
  //
  // ---------------------------------------------------------------------------------
  //                                    Methods
  // ---------------------------------------------------------------------------------
  //
  virtual void   UserCreateOutputObjects();            // create output objects
  virtual Bool_t Run();                                // Equivalent to UserExec --> tuned in
  virtual void   Terminate(Option_t *);                // run only once and terminate

  // ---------------------------------------------------------------------------------
  //                                    Settings
  // ---------------------------------------------------------------------------------

  void   Initialize();
  void   PrintNumInBinary(UInt_t num);
  void   SetESDtrackCuts(AliESDtrackCuts * trackCuts)                 {fESDtrackCuts        = trackCuts;};
  void   SetIsMCtrue(Bool_t isMCdata = kTRUE)                         {fMCtrue              = isMCdata;};
  void   SetYear(const Int_t ifYear = 0)                              {fYear                = ifYear;}
  void   SetPeriodName(const TString ifPeriodName = "")               {fPeriodName          = ifPeriodName;}
  void   SetPassIndex(const Int_t ifPassIndex = 0)                    {fPassIndex           = ifPassIndex;}
  // Some boolian settings
  void   SetRunOnGrid(const Bool_t ifRunOnGrid = kTRUE)               {fRunOnGrid           = ifRunOnGrid;}
  void   SetIncludeITScuts(const Bool_t ifITSCuts = kTRUE)            {fIncludeITS          = ifITSCuts;}
  void   SetFillArmPodTree(const Bool_t ifArmpodTree = kTRUE)         {fFillArmPodTree      = ifArmpodTree;}
  void   SetFillBGJetsFJTree(const Bool_t ifBGJetsFJTree = kTRUE)     {fFillBGJetsFJTree    = ifBGJetsFJTree;}
  void   SetFilldscaledTree(const Bool_t ifdscaledTree = kTRUE)       {fFilldscaledTree     = ifdscaledTree;}
  void   SetFillFastJet(const Bool_t ifFastJet = kTRUE)               {fFillFastJet         = ifFastJet;}
  void   SetDeDxCheck(const Bool_t ifDeDxCheck = kFALSE)              {fDEdxCheck           = ifDeDxCheck;}
  void   SetEffMatrix(const Bool_t ifEffMatrix = kFALSE)              {fEffMatrix           = ifEffMatrix;}
  void   SetFillAllCutVariables(const Bool_t ifAllCuts = kFALSE)      {fFillTracks          = ifAllCuts;}
  void   SetFillOnlyHists(const Bool_t ifFillOnlyHists = kFALSE)      {fFillOnlyHists       = ifFillOnlyHists;}
  void   SetFillEffLookUpTable(const Bool_t ifEffLookUpTable = kFALSE){fFillEffLookUpTable  = ifEffLookUpTable;}
  void   SetFillHigherMomentsMCclosure(const Bool_t ifHigherMomentsMCclosure = kFALSE){fFillHigherMomentsMCclosure  = ifHigherMomentsMCclosure;}
  void   SetRunFastSimulation(const Bool_t ifFastSimul = kFALSE)      {fRunFastSimulation   = ifFastSimul;}
  void   SetRunFastHighMomentCal(const Bool_t ifFastHighMom = kFALSE) {fRunFastHighMomentCal= ifFastHighMom;}
  void   SetFillDistributions(const Bool_t ifGenDistributions = kFALSE) {fFillDistributions= ifGenDistributions;}
  void   SetFillTreeMC(const Bool_t ifTreeMC = kFALSE)                {fFillTreeMC= ifTreeMC;}

  void   SetDefaultTrackCuts(const Bool_t ifDefaultTrackCuts = kFALSE){fDefaultTrackCuts= ifDefaultTrackCuts;}
  void   SetDefaultEventCuts(const Bool_t ifDefaultEventCuts = kFALSE){fDefaultEventCuts= ifDefaultEventCuts;}
  void   SetFillNudynFastGen(const Bool_t ifNudynFastGen = kFALSE)    {fFillNudynFastGen= ifNudynFastGen;}
  void   SetUsePtCut(const Int_t ifUsePtCut = 1)                      {fUsePtCut            = ifUsePtCut;}
  void   SetMCTrackOriginType(const Int_t ifTrackOriginOnlyPrimary = 0) {fTrackOriginOnlyPrimary     = ifTrackOriginOnlyPrimary;}
  void   SetRapidityType(const Int_t ifRapidityType = 0)              {fRapidityType        = ifRapidityType;}
  void   SetSisterCheck(const Int_t ifSisterCheck = 0)                {fSisterCheck         = ifSisterCheck;}
  void   SetFillDnchDeta(const Bool_t ifDnchDetaCal = kFALSE)         {fFillDnchDeta        = ifDnchDetaCal;}
  void   SetIncludeTOF(const Bool_t ifIncludeTOF = kFALSE)            {fIncludeTOF          = ifIncludeTOF;}
  void   SetUseThnSparse(const Bool_t ifUseThnSparse = kFALSE)        {fUseThnSparse        = ifUseThnSparse;}
  void   SetUseCouts(const Bool_t ifUseCouts = kFALSE)                {fUseCouts            = ifUseCouts;}
  void   SetWeakAndMaterial(const Bool_t ifWeakAndMaterial = kFALSE)  {fWeakAndMaterial     = ifWeakAndMaterial;}
  void   SetPercentageOfEvents(const Int_t nPercentageOfEvents = 0)   {fPercentageOfEvents = nPercentageOfEvents;}
  void   SetNSettings(const Int_t nSettings = 22)                     {fNSettings = nSettings;}

  //
  Bool_t GetRunOnGrid() const { return fRunOnGrid; }

  // Setters for the systematic uncertainty checks
  void   SetSystDCAxy(const Int_t systDCAxy = 0)                  {fSystDCAxy           = systDCAxy;}
  void   SetSystNCrossedRows(const Int_t systNCrossedRows = 0)    {fSystCrossedRows     = systNCrossedRows;}
  void   SetSystTPCChi2(const Int_t systTPCChi2 = 0)              {fSystChi2            = systTPCChi2;}
  void   SetSystVz(const Int_t systVz = 0)                        {fSystVz              = systVz;}

  // Setters for the eta momentum dEdx and centrality bins
  void   SetSampleDeDxUpperEdge(const Float_t dEdxCleanUp = 200.) {fDEdxCleanUp         = dEdxCleanUp;}
  void   SetDeDxBinWidth(const Float_t dEdxBinWidth = 2.5)        {fDEdxBinWidth        = dEdxBinWidth;}
  void   SetDeDxLowerEdge(const Float_t dEdxLowerEdge = 20.)      {fDEdxDown            = dEdxLowerEdge;}
  void   SetDeDxUpperEdge(const Float_t dEdxUpperEdge = 1020.)    {fDEdxUp              = dEdxUpperEdge;}

  void   SetEtaLowerEdge(const Float_t etaLowerEdge = -0.8)      {fEtaDown             = etaLowerEdge;}
  void   SetEtaUpperEdge(const Float_t etaUpperEdge = 0.8)       {fEtaUp               = etaUpperEdge;}
  void   SetNEtabins(const Int_t nEtaBins = 20)                   {fNEtaBins            = nEtaBins;}
  void   SetMomLowerEdge(const Float_t momLowerEdge = 0.)         {fMomDown             = momLowerEdge;}
  void   SetMomUpperEdge(const Float_t momUpperEdge = 12.)        {fMomUp               = momUpperEdge;}
  void   SetNMomBins(const Int_t nMombins = 600)                  {fNMomBins            = nMombins;}
  void   SetNGenprotonBins(const Int_t nGenprotonBins = 100)      {fGenprotonBins       = nGenprotonBins;}
  //function to add jet task
  void   AddJet(AliJetContainer* jet = 0)                         {fJetContainer        = jet; }

  // Set the binning of centrality
  void SetCentralityBinning(const Int_t tmpCentbins, Float_t tmpfxCentBins[])
  {
    // Create the histograms to be used in the binning of eta, cent and momentum
    std::cout << " Info::siweyhmi: !!!!!! Centrality binning is being set !!!!!!! " << std::endl;
    fHistCent =  new TH1F("fHistCent","Centrality Bins",tmpCentbins-1    ,tmpfxCentBins );
    fHistPhi  =  new TH1F("fHistPhi" ,"Phi Bins"       ,36               ,-TMath::Pi(), TMath::Pi());
    // ==========================================
    // prepare real data centrality bins
    fNCentbinsData = tmpCentbins;
    fNCentBinsMC   = tmpCentbins-1;
    fxCentBins.resize(fNCentbinsData);
    for (Int_t i=0; i<fNCentbinsData; i++) fxCentBins[i] =  tmpfxCentBins[i];
    fcentDownArr.resize(fNCentBinsMC);
    fcentUpArr.resize(fNCentBinsMC);
    for (Int_t i=0; i<fNCentbinsData-1; i++) fcentDownArr[i] =  tmpfxCentBins[i];
    for (Int_t i=1; i<fNCentbinsData; i++)   fcentUpArr[i-1] =  tmpfxCentBins[i];
  }

  void SetMCEtaScanArray(const Int_t tmpEtaBinsMC, Float_t tmpetaDownArr[], Float_t tmpetaUpArr[])
  {
    // set MC eta values to scan
    std::cout << " Info::siweyhmi: !!!!!! SetMCEtaScanArray is being set !!!!!!! " << std::endl;
    fNEtaWinBinsMC = tmpEtaBinsMC;
    fetaDownArr.resize(fNEtaWinBinsMC);
    fetaUpArr.resize(fNEtaWinBinsMC);
    for (Int_t i=0; i<fNEtaWinBinsMC; i++) {
      fetaDownArr[i] =  tmpetaDownArr[i];
      fetaUpArr[i]   =  tmpetaUpArr[i];
    }
  }

  void SetMCResonanceArray(const Int_t tmpNRes, TString tmpResArr[])
  {
    // set MC eta values to scan
    std::cout << " Info::siweyhmi: !!!!!! SetMCResonanceArray is being set !!!!!!! " << std::endl;
    fNResBins = tmpNRes;
    fResonances.resize(fNResBins);
    for (Int_t i=0; i<fNResBins; i++) fResonances[i] = tmpResArr[i];

  }

  void SetMCBaryonArray(const Int_t tmpNBar, Int_t tmpBarArr[])
  {
    // set MC eta values to scan
    std::cout << " Info::siweyhmi: !!!!!! SetMCBaryonArray is being set !!!!!!! " << std::endl;
    fNBarBins = tmpNBar;
    fBaryons.resize(fNBarBins);
    for (Int_t i=0; i<fNBarBins; i++) fBaryons[i] = tmpBarArr[i];
  }

  void SetMCMomScanArray(const Int_t tmpMomBinsMC, Float_t tmppDownArr[], Float_t tmppUpArr[])
  {
    // set MC momentum values to scan
    std::cout << " Info::siweyhmi: !!!!!! SetMCMomScanArray is being set !!!!!!! " << std::endl;
    fNMomBinsMC = tmpMomBinsMC;
    fpDownArr.resize(fNMomBinsMC);
    fpUpArr.resize(fNMomBinsMC);
    for (Int_t i=0; i<fNMomBinsMC; i++) {
      fpDownArr[i] =  tmppDownArr[i];
      fpUpArr[i]   =  tmppUpArr[i];
    }
  }

private:

  AliAnalysisJetHadro(const AliAnalysisJetHadro&);
  AliAnalysisJetHadro& operator=(const AliAnalysisJetHadro&);

  // ---------------------------------------------------------------------------------
  //                                   Functions
  // ---------------------------------------------------------------------------------

  void FindJetsEMC();                          // Find and Fill Jets with EMCAL framework
  void FindJetsFJ();                          // Find and Fill Jets with FJ framework
  void FillTPCdEdxReal();                   // Main function to fill all info + TIden
  void FillTrackVariables(AliESDtrack *track);
  void FillMCFull();                        // Fill all info + TIdenMC from MC to do MC closure test
  void WeakAndMaterial();                   // Look full acceptance, weak decay and material
  void FillEffMatrix();                     // Prepare efficiency matrix
  void FillCleanSamples();                    // Fill Clean Pions
  void SelectCleanSamplesFromV0s(AliESDv0 *v0, AliESDtrack *track0, AliESDtrack *track1);
  void SetSpecialV0Cuts(AliESDv0KineCuts* cuts);
  void BinLogAxis(TH1 *h);
  void GetExpecteds(AliESDtrack *track, Float_t closestPar[3]);
  void FindJetsFJGen();
  void FillEventTree();


  //
  Int_t CountEmptyEvents(Int_t counterBin);  // Just count if there is empty events
  UInt_t SetCutBitsAndSomeTrackVariables(AliESDtrack *track, Int_t particleType);
  Bool_t CheckIfFromResonance(Int_t mcType, AliMCParticle *trackMCgen, Int_t trackIndex, Bool_t parInterest, Double_t ptot, Double_t eta, Double_t cent, Bool_t fillTree);
  Bool_t CheckIfFromAnyResonance(AliMCParticle *trackMCgen, Float_t etaLow, Float_t etaUp, Float_t pDown, Float_t pUp);
  Bool_t ApplyDCAcutIfNoITSPixel(AliESDtrack *track);
  // ---------------------------------------------------------------------------------
  //                                   Members
  // ---------------------------------------------------------------------------------

  AliPIDResponse   * fPIDResponse;            //! PID response object
  AliESDEvent      * fESD;                    //! ESD object
  TList            * fListHist;               //! list for histograms
  AliESDtrackCuts  * fESDtrackCuts;           //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit96;     //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit96_spd;     //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit96_sdd;     //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit128;    //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit768;    //! basic cut variables
  AliESDtrackCuts  * fESDtrackCutsLoose;      //! basic cut variables for debugging
  AliESDv0Cuts     * fESDtrackCutsV0;         //! basic cut variables for V0
  AliESDtrackCuts  * fESDtrackCutsCleanSamp;  //! basic cut variables for clean pion and electron form V0s
  AliPIDCombined   * fPIDCombined;            //! combined PID object
  AliTPCdEdxInfo   * fTPCdEdxInfo;            //! detailed dEdx info
  AliStack         * fMCStack;                //! stack object to get Mc info
  AliESDv0KineCuts * fV0OpenCuts;             // v0 strong filter for tagged V0s
  AliESDv0KineCuts * fV0StrongCuts;           // v0 strong filter for tagged V0s
  AliAnalysisCuts  * fK0sPionCuts;            // filter for pions from K0s
  AliAnalysisCuts  * fLambdaProtonCuts;       // filter for protons from Lambda
  AliAnalysisCuts  * fLambdaPionCuts;         // filter for pions from Lambda
  AliAnalysisCuts  * fGammaElectronCuts;      // filter for electrons from gamma conversions
  const AliESDVertex * fVertex;               // primary vertex
//  AliESDtools      * fESDtool;                 // tools to calculate derived variables from the ESD

  TTree            * fArmPodTree;             // Tree for clean pion and proton selection
  TTreeSRedirector * fTreeSRedirector;        /// temp tree to dump output
  TTree            * fTreejetsEMCconst;       // tree for EMCal signal jet constituents
  TTree            * fTreeMC;                 // tree for mc samples
  TTree            * fTreejetsFJ;             // tree for fastjet signal jets
  TTree            * fTreeBGjetsFJ;           // tree for fastjet background jets
  TTree            * fTreeCuts;               // tree to save all variables for control plots
  TTree            * fTreejetsFJconst;          // tree for fastjet signal jet constituents
  TTree            * fTreeResonance;          // tree with full acceptance filled with MC
  TTree            * fTreeEvents;
  TTree            * fTreeDScaled;
  TTree            * fTreejetsFJGen;
  TTree            * fTreejetEMC;            // tree for EMCal signal jets
  TRandom          fRandom;


  TString            fPeriodName;
  Int_t              fYear;
  Int_t              fPassIndex;
  UInt_t             fPileUpBit;
  TH1F             * fHistCent;               // helper histogram for TIdentity tree
  TH1F             * fHistPhi;
  TH1F             * fHistGenMult;
  TH2F             * fHistRapDistFullAccPr;
  TH2F             * fHistRapDistFullAccAPr;
  TH1F             * fHistInvK0s;             // helper histogram for TIdentity tree
  TH1F             * fHistInvLambda;          // helper histogram for TIdentity tree
  TH1F             * fHistInvAntiLambda;      // helper histogram for TIdentity tree
  TH1F             * fHistInvPhoton;          // helper histogram for TIdentity tree
  //
  THnSparseF       * fHndEdx;                 // histogram which hold all dEdx info
  THnSparseF       * fHnExpected;         // histogram which hold all PIDresponse info

  TString           fChunkName;

  UInt_t            fTrackCutBits;           // integer which hold all cut variations as bits
  Int_t             fSystClass;
  Double_t          fEtaDown;
  Double_t          fEtaUp;
  Int_t             fNEtaBins;
  Int_t             fPercentageOfEvents;     // when only a fPercentageOfEvents is enough

  Bool_t            fRunOnGrid;              // flag if real data or MC is processed
  Bool_t            fMCtrue;                 // flag if real data or MC is processed
  Bool_t            fWeakAndMaterial;        // flag for the Weak and Material analysis
  Bool_t            fEffMatrix;              // flag for efficiency matrix filling
  Bool_t            fDEdxCheck;              // flag to check only the dEdx performance
  Bool_t            fIncludeITS;             // decide whether to use ITS or not
  Bool_t            fFillTracks;             // switch whether to fill all cut variables
  Bool_t            fFillOnlyHists;          //
  Bool_t            fFillEffLookUpTable;     //
  Bool_t            fFillHigherMomentsMCclosure;
  Bool_t            fFillArmPodTree;         // switch whether to fill clean sample tree
  Bool_t            fFillBGJetsFJTree;         // switch whether to fill BG Jets FJ tree
  Bool_t            fFilldscaledTree;         // switch whether to fill dscaled tree
  Bool_t            fFillFastJet;         // switch whether to fill dscaled tree
  Bool_t            fRunFastSimulation;      // when running over galice.root do not fill other objects
  Bool_t            fRunFastHighMomentCal;   // when running over galice.root do not fill other objects
  Bool_t            fFillDistributions;   // when running over galice.root do not fill other objects
  Bool_t            fFillTreeMC;
  Bool_t            fDefaultTrackCuts;
  Bool_t            fDefaultEventCuts;
  Bool_t            fFillNudynFastGen;
  Int_t             fUsePtCut;
  Int_t             fTrackOriginOnlyPrimary;
  Int_t             fRapidityType;
  Int_t             fSisterCheck;           // 0: reject the mother anyways, 1: if both girls are in acceptance rejet mother

  Bool_t            fFillDnchDeta;         // switch on calculation of the dncdeta for fastgens
  Bool_t            fIncludeTOF;             // Include TOF information to investigate the efficiency loss effects on observable
  Bool_t            fUseThnSparse;           // in case thnsparse is filled
  Bool_t            fUseCouts;               // for debugging

  Int_t             fNSettings;
  Int_t             fNMomBins;               // number of mombins --> for 20MeV slice 150 and 10MeV 300
  Float_t           fMomDown;                // bottom limit for the momentum range (default 0.2)
  Float_t           fMomUp;                  // uppper limit for the momentum range (default 3.2)
  Float_t           fDEdxBinWidth;           // bin width for the dEdx histograms (default 2.5)
  Float_t           fDEdxUp;                 // bottom limit for dEdx histogram (default 20)
  Float_t           fDEdxDown;               // upper limit for dEdx histogram (default 1020)
  Float_t           fDEdxCleanUp;            // upper limit for dEdx histogram of clean kaons and electrons (default 140)

  Float_t           fArmPodTPCSignal;
  Float_t           fArmPodptot;
  Float_t           fArmPodEta;
  Float_t           fArmPodCentrality;
  Float_t           fQt;
  Float_t           fAlfa;
  Float_t           fNSigmasElTOF;           // TOF N sigma for Electron
  Float_t           fNSigmasPiTOF;           // TOF N sigma for Pion
  Float_t           fNSigmasKaTOF;           // TOF N sigma for Kaon
  Float_t           fNSigmasPrTOF;           // TOF N sigma for Proton
  Float_t           fNSigmasDeTOF;           // TOF N sigma for Deuteron

  Float_t           fDEdxEl;                 // Expected Electron dEdx
  Float_t           fDEdxKa;                 // Expected Kaon dEdx
  Float_t           fDEdxPi;                 // Expected Pion dEdx
  Float_t           fDEdxPr;                 // Expected Proton dEdx
  Float_t           fDEdxDe;                 // Expected Deuteron dEdx

  Float_t           fSigmaEl;                // Expected Electron sigma
  Float_t           fSigmaKa;                // Expected Kaon sigma
  Float_t           fSigmaPi;                // Expected Pion sigma
  Float_t           fSigmaPr;                // Expected Proton sigma
  Float_t           fSigmaDe;                // Expected Deuteron sigma

  Float_t           fNSigmasElTPC;           // TOF N sigma for Electron
  Float_t           fNSigmasPiTPC;           // TOF N sigma for Pion
  Float_t           fNSigmasKaTPC;           // TOF N sigma for Kaon
  Float_t           fNSigmasPrTPC;           // TOF N sigma for Proton
  Float_t           fNSigmasDeTPC;           // TOF N sigma for Deuteron

  Float_t           fTPCSignalMC;
  Float_t           fPtotMC;
  Float_t           fPtotMCtruth;
  Float_t           fPtMC;
  Float_t           fEtaMC;
  Int_t             fSignMC;

  Float_t           fPxMC;                     // x component of momentum
  Float_t           fPyMC;                     // y component of momentum
  Float_t           fPzMC;                     // z component of momentum

  Float_t           fElMC;
  Float_t           fPiMC;
  Float_t           fKaMC;
  Float_t           fPrMC;
  Float_t           fDeMC;
  Float_t           fMuMC;
  Float_t           fLaMC;

  Double_t          fMCImpactParameter;
  Int_t             fNHardScatters;            // Number of hard scatterings
  Int_t             fNProjectileParticipants;  // Number of projectiles participants
  Int_t             fNTargetParticipants;      // Number of target participants
  Int_t             fNNColl;                   // Number of N-N collisions
  Int_t             fNNwColl;                  // Number of N-Nwounded collisions
  Int_t             fNwNColl;                  // Number of Nwounded-N collisons
  Int_t             fNwNwColl;                 // Number of Nwounded-Nwounded collisions


  Float_t           fElMCgen;
  Float_t           fPiMCgen;
  Float_t           fKaMCgen;
  Float_t           fPrMCgen;
  Float_t           fDeMCgen;
  Float_t           fMuMCgen;
  Float_t           fLaMCgen;
  Float_t           fBaMCgen;


  Float_t           fPx;                     // x component of momentum
  Float_t           fPy;                     // y component of momentum
  Float_t           fPz;                     // z component of momentum
  Float_t           fPtot;                   // TPC momentum
  Float_t           fPVertex;                // TPC momentum
  Float_t           fPt;                     // Transverse momentum
  Float_t           fY;                      // rapidity

  Int_t              fMultiplicity;           // Multiplicity in case of PbPb
  Int_t              fMultiplicityMC;
  Float_t            fCentrality;             // centrality information
  Float_t            fCentImpBin;
  Double_t           fVz;                     // Vertex position
  ULong64_t          fEventGID;               // global Event Id
  Int_t              fEventGIDMC;             // global MC event id
  Int_t              fEventCountInFile;       // event count per job
  Int_t              fEvent;                  // Event counter for Christian
  Int_t              fEventMC;                // Event id for MC data
  Int_t              fEventMCgen;             // Event id for MC generated

  Float_t            fTPCSignal;              // Measured dE/dx
  Float_t            fEta;                    // pseudo rapidity
  Float_t            fNContributors;          // Ntracks
  Float_t            fTheta;                  // theta
  Float_t            fPhi;                    // azimuthal angle
  Int_t              fSign;                   // sign of the particle
  Int_t              fTPCShared;              // number of shared clusters
  Int_t              fNcl;                    // number of points used for dEdx
  Int_t              fNclCorr;                // number of points used for dEdx

  Int_t              fNResBins;
  Int_t              fNBarBins;
  Int_t              fNEtaWinBinsMC;
  Int_t              fNMomBinsMC;
  Int_t              fNCentBinsMC;
  Int_t              fGenprotonBins;
  Int_t              fNResModeMC;
  Int_t              fNCentbinsData;
  Float_t            fMissingCl;
  Int_t              fTPCMult;
  Int_t              fEventMult;
  Double_t           fTimeStamp;
  Float_t            fIntRate;
  Int_t              fRunNo;
  Float_t            fBField;
  TString            fBeamType;
  AliJetContainer*   fJetContainer;
  Double_t           fJetPt;
  Double_t           fJetEta;
  Double_t           fJetPhi;
  Float_t            fjetRhoVal;
  Float_t            frhoFJ;
  Bool_t             fhasAcceptedFJjet;
  Bool_t             fhasAcceptedEMCjet;
  Bool_t             fhasRealFJjet;
  Bool_t             fhasRealEMCjet;
  // Cut variables
  Double_t fTrackProbElTPC;
  Double_t fTrackProbPiTPC;
  Double_t fTrackProbKaTPC;
  Double_t fTrackProbPrTPC;
  Bool_t   fTrackProbDeTPC;
  Double_t fTrackProbElTOF;
  Double_t fTrackProbPiTOF;
  Double_t fTrackProbKaTOF;
  Double_t fTrackProbPrTOF;
  Bool_t   fTrackProbDeTOF;

  Float_t fTrackTPCCrossedRows;
  Float_t fTrackChi2TPC;
  Float_t fTrackChi2TPCcorr;
  Float_t fTrackDCAxy;
  Float_t fTrackDCAz;
  Float_t fTrackLengthInActiveZone;
  Float_t fTrackTPCSignalN;
  Bool_t  fTrackIsFirstITSlayer;
  Bool_t  fTrackIsSecondITSlayer;
  Bool_t  fTrackNewITScut;
  Bool_t  fTrackRequireITSRefit;

  // Additional cuts from marian
  Bool_t             fIsITSpixel01;           // if track has hits in innermost 2 pixels of ITS
  Int_t              fNITSclusters;           // number of ITS clusters
  Float_t            fPrimRestriction;        // prim vertex cut recommended by marian
  Float_t            fTPCvZ;                  // TPC vertex
  Float_t            fSPDvZ;                  // SPD vertex

  //   CleanSample cuts
  Bool_t             fCleanPionsFromK0;
  Bool_t             fCleanPion0FromK0;
  Bool_t             fCleanPion1FromK0;
  Bool_t             fCleanPion0FromLambda;
  Bool_t             fCleanPion1FromLambda;
  Bool_t             fCleanProton0FromLambda;
  Bool_t             fCleanProton1FromLambda;
  Bool_t             fHasTrack0FirstITSlayer;
  Bool_t             fHasTrack1FirstITSlayer;
  Bool_t             fHasV0FirstITSlayer;

  //  Variables for systematic uncertainty checks
  //  B field configurations -->  use default settings and analyse the following set of runs
  //  ------------------------------------------------
  //  Field (++)  --> run interval is [137161, 138275]
  //  Field (--)  --> run interval is [138364, 139510]
  //  ------------------------------------------------
  Int_t              fSystCrossedRows;       // 0 -->  80     ||| -1 -->   60   ||| +1 -->  100
  Int_t              fSystDCAxy;             // 0 --> default ||| -1 --> -sigma ||| +1 --> +sigma
  Int_t              fSystChi2;              // 0 -->  4      ||| -1 -->    3   ||| +1 -->   5
  Int_t              fSystVz;                // 0 -->  10     ||| -1 -->    8   ||| +1 -->   12

  std::vector<float>  fetaDownArr;
  std::vector<float>  fetaUpArr;
  std::vector<float>  fcentDownArr;
  std::vector<float>  fcentUpArr;
  std::vector<float>  fpDownArr;
  std::vector<float>  fpUpArr;
  std::vector<float>  fxCentBins;
  std::vector<int>    fBaryons;
  std::vector<std::string> fResonances;
  //
  // control and QA histograms
  //
  TGraph           * fEventInfo_LumiGraph;           // grap for the interaction rate info for a run
  THnF             * fHistPosEffMatrixRec;       // histogram efficiency matrix --> reconstructed traks
  THnF             * fHistNegEffMatrixRec;       // histogram efficiency matrix --> generated traks
  THnF             * fHistPosEffMatrixGen;       // histogram efficiency matrix --> reconstructed pions
  THnF             * fHistNegEffMatrixGen;       // histogram efficiency matrix --> generated pions
  THnF             * fHistPosEffMatrixScanRec;   // histogram efficiency matrix --> reconstructed traks
  THnF             * fHistNegEffMatrixScanRec;   // histogram efficiency matrix --> generated traks
  THnF             * fHistPosEffMatrixScanGen;   // histogram efficiency matrix --> reconstructed pions
  THnF             * fHistNegEffMatrixScanGen;   // histogram efficiency matrix --> generated pions
  TH1F             * fHistEmptyEvent;            // control histogram for empty event
  TH1F             * fHistCentrality;            // control histogram for centrality
  TH1F             * fHistCentralityImpPar;      // control histogram for centrality
  TH1F             * fHistImpParam;              // control histogram for impact parameter
  TH1F             * fHistVertex;                // control histogram for vertexZ
  THnF             * fHistdEdxTPC;               // 5D hist of dEdx from all TPC
  TH2F             * fHistArmPod;                // control histogram for Armanteros Podolanski plot
  //
  // Counters for Marian
  //
  AliEventCuts* fPileUpTightnessCut4;
  AliEventCuts* fPileUpTightnessCut3;
  AliEventCuts* fPileUpTightnessCut2;
  AliEventCuts* fPileUpTightnessCut1;

  ClassDef(AliAnalysisJetHadro, 6);

};

#endif
