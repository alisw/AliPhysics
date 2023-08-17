#ifndef ALIANALYSISJETHADRO_H
#define ALIANALYSISJETHADRO_H

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
class TH3F;
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
#include "THnSparse.h"
#include "TClonesArray.h"
#include "THn.h"
#include "TVectorF.h"
#include "TCutG.h"
#include "TTreeStream.h"
#include "AliEventCuts.h"
#include "TRandom3.h"

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
    kNCrossedRowsTPC70=0,
    kNCrossedRowsTPC80=1,
    kNCrossedRowsTPC90=2,
    kMaxChi2PerClusterTPCSmall=3,
    kMaxChi2PerClusterTPC=4,
    kMaxChi2PerClusterTPCLarge=5,
    kMaxDCAToVertexXYPtDep=6,
    kMaxDCAToVertexXYPtDepLarge=7,
    kVertexZSmall=8,
    kVertexZ=9,
    kEventVertexZ=10,
    kEventVertexZLarge=11,
    kActiveZone=12,
    kTPCSignalNSmall=13,
    kTPCSignalN=14,
    kTPCSignalNLarge=16,
    kCleanPrTOF=17,
    kCleanKaTOF=18,
    kCleanKaTOFTRD=19,
    kTrackProbKaTOF=20,
    kTrackProbPrTOF=21,
    kCleanDeTOF=22,
    kPileup=23,
    kPileupLoose=24,
    kSharedCls=25,
    kSharedClsLoose=26,
    kFindableCls=27,
    kFindableClsTight=28,
    kFindableClsLoose=29,
    kBFieldPos=30,
    kBFieldNeg=31
  };

  enum cutSettings {
    kCutReference=0,
    kCutCrossedRowsTPC70=1,
    kCutCrossedRowsTPC90=2,
    kCutActiveZone=3,
    kCutMaxChi2PerClusterTPCSmall=4,
    kCutMaxChi2PerClusterTPCLarge=5,
    kCutMaxDCAToVertexXYPtDepLarge=6,
    kCutVertexZSmall=7,
    kCutEventVertexZLarge=8,
    kCutSharedCls=9,
    kCutFindableClsTight=10,
    kCutFindableClsLoose=11,
    kCutPileupLoose=12,
    kCutBFieldPos=13,
    kCutBFieldNeg=14,
    kCutTPCSignalNSmall=15,
    kCutTPCSignalNLarge=16
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
  void   SetSmallOut(const Bool_t ifSmallOut = kTRUE)               {fSmallOut           = ifSmallOut;}
  void   SetIncludeITScuts(const Bool_t ifITSCuts = kTRUE)            {fIncludeITS          = ifITSCuts;}
  void   SetFilljetsFJBGTree(const Bool_t ifjetsFJBGTree = kTRUE)     {fFilljetsFJBGTree    = ifjetsFJBGTree;}
  void   SetFillJetsFJBGConst(const Bool_t ifJetsFJBGConst = kTRUE)     {fFillJetsFJBGConst     = ifJetsFJBGConst;}
  void   SetFilldscaledTree(const Bool_t ifdscaledTree = kTRUE)       {fFilldscaledTree     = ifdscaledTree;}
  void   SetDoFastJet(const Bool_t ifFastJet = kTRUE)               {fDoFastJet         = ifFastJet;}
  void   SetDoEMCJet(const Bool_t ifEMCJet = kTRUE)               {fDoEMCJet         = ifEMCJet;}
  void   SetFillFastJet(const Bool_t ifFastJet = kTRUE)               {fFillFastJet         = ifFastJet;}
  void   SetFillEMCJet(const Bool_t ifEMCJet = kTRUE)               {fFillEMCJet         = ifEMCJet;}
  void   SetFillJetsEMCConst(const Bool_t ifJetsEMCConst = kTRUE)     {fFillJetsEMCConst     = ifJetsEMCConst;}
  void   SetFillJetsFJConst(const Bool_t ifJetsFJConst = kTRUE)     {fFillJetsFJConst     = ifJetsFJConst;}
  void   SetFillJetsEMCBG(const Bool_t ifJetsEMCBG = kTRUE)     {fFillJetsEMCBG     = ifJetsEMCBG;}
  void   SetFillJetsEMCBGConst(const Bool_t ifJetsEMCBGConst = kTRUE)     {fFillJetsEMCBGConst     = ifJetsEMCBGConst;}
  void   SetJetMinPtSub(const Double_t jetminptsub = -1000.0)         {fjetMinPtSub         = jetminptsub;}
  void   SetJetMinArea(const Double_t jetminarea = -1000.0)         {fjetMinArea         = jetminarea;}
  void   SetMinCent(const Double_t mincent = 0.0)         {fcent_min         = mincent;}
  void   SetMaxCent(const Double_t maxcent = 100.0)         {fcent_max         = maxcent;}

  void   SetDeDxCheck(const Bool_t ifDeDxCheck = kFALSE)              {fDEdxCheck           = ifDeDxCheck;}
  void   SetFillOnlyHists(const Bool_t ifFillOnlyHists = kFALSE)      {fFillOnlyHists       = ifFillOnlyHists;}
  void   SetFillEffLookUpTable(const Bool_t ifEffLookUpTable = kFALSE){fFillEffLookUpTable  = ifEffLookUpTable;}
  void   SetRunFastSimulation(const Bool_t ifFastSimul = kFALSE)      {fRunFastSimulation   = ifFastSimul;}
  void   SetFillTreeMC(const Bool_t ifTreeMC = kFALSE)                {fFillTreeMC= ifTreeMC;}
  void   SetFillIncTracks(const Bool_t ifIncTracks = kTRUE)           {fFillIncTracks       = ifIncTracks;}
  void   SetFillpTPC(const Bool_t ifpTPC = kTRUE)           {fFillpTPC       = ifpTPC;}
  void   SetFillp(const Bool_t ifp = kTRUE)           {fFillp       = ifp;}
  void   SetFillpT(const Bool_t ifpT = kTRUE)           {fFillpT       = ifpT;}

  void   SetDefaultTrackCuts(const Bool_t ifDefaultTrackCuts = kFALSE){fDefaultTrackCuts= ifDefaultTrackCuts;}
  void   SetDefaultEventCuts(const Bool_t ifDefaultEventCuts = kFALSE){fDefaultEventCuts= ifDefaultEventCuts;}
  void   SetCorrectForMissCl(const Int_t ifCorrectForMissCl = kFALSE)    {fCorrectForMissCl= ifCorrectForMissCl;}
  void   SetUsePtCut(const Int_t ifUsePtCut = 1)                      {fUsePtCut            = ifUsePtCut;}
  void   SetMCTrackOriginType(const Int_t ifTrackOriginOnlyPrimary = 0) {fTrackOriginOnlyPrimary     = ifTrackOriginOnlyPrimary;}
  void   SetRapidityType(const Int_t ifRapidityType = 0)              {fRapidityType        = ifRapidityType;}
  void   SetSisterCheck(const Int_t ifSisterCheck = 0)                {fSisterCheck         = ifSisterCheck;}
  void   SetIncludeTOF(const Bool_t ifIncludeTOF = kFALSE)            {fIncludeTOF          = ifIncludeTOF;}
  void   SetUseCouts(const Bool_t ifUseCouts = kFALSE)                {fUseCouts            = ifUseCouts;}
  void   SetFillEventInfo(const Bool_t ifEventInfo = kFALSE)          {fEventInfo           = ifEventInfo;}
  void   SetPercentageOfEvents(const Int_t nPercentageOfEvents = 0)   {fPercentageOfEvents = nPercentageOfEvents;}
  void   SetNSettings(const Int_t nSettings = 22)                     {fNSettings = nSettings;}
  void   SetRunNumberForExpecteds(const Int_t ifRunNumberForExpecteds = 0)    {fRunNumberForExpecteds = ifRunNumberForExpecteds;}

  //
  Bool_t GetRunOnGrid() const { return fRunOnGrid; }

  // Setters for the systematic uncertainty checks
  void   SetSystCentEstimator(const Int_t systCentEstimator = 0)  {fSystCentEstimatetor = systCentEstimator;}

  // Setters for the eta momentum dEdx and centrality bins
  void   SetSampleDeDxUpperEdge(const Float_t dEdxCleanUp = 200.) {fDEdxCleanUp         = dEdxCleanUp;}
  void   SetNDeDxBins(const Float_t ndEdxBins = 2000)              {fNdEdxBins        = ndEdxBins;}
  void   SetDeDxLowerEdge(const Float_t dEdxLowerEdge = 0.)      {fDEdxDown            = dEdxLowerEdge;}
  void   SetDeDxUpperEdge(const Float_t dEdxUpperEdge = 2000.)    {fDEdxUp              = dEdxUpperEdge;}

  void   SetEtaLowerEdge(const Float_t etaLowerEdge = -0.8)      {fEtaDown             = etaLowerEdge;}
  void   SetEtaUpperEdge(const Float_t etaUpperEdge = 0.8)       {fEtaUp               = etaUpperEdge;}
  void   SetNEtabins(const Int_t nEtaBins = 20)                   {fNEtaBins            = nEtaBins;}
  void   SetMomLowerEdge(const Float_t momLowerEdge = 0.)         {fMomDown             = momLowerEdge;}
  void   SetMomUpperEdge(const Float_t momUpperEdge = 12.)        {fMomUp               = momUpperEdge;}
  void   SetNMomBins(const Int_t nMombins = 600)                  {fNMomBins            = nMombins;}
  void   SetNGenprotonBins(const Int_t nGenprotonBins = 100)      {fGenprotonBins       = nGenprotonBins;}

  void   SetMomExpec_LowerEdge(const Float_t momLowerEdge = 0.)         {fMomExpec_Low             = momLowerEdge;}
  void   SetMomExpec_UpperEdge(const Float_t momUpperEdge = 20.0)        {fMomExpec_High               = momUpperEdge;}
  void   SetMomExpec_NBins(const Int_t nMombins = 2000)                  {fMomExpec_NBins            = nMombins;}

  void   SetEtaExpec_LowerEdge(const Float_t EtaLowerEdge = 0.)         {fEtaExpec_Low             = EtaLowerEdge;}
  void   SetEtaExpec_UpperEdge(const Float_t EtaUpperEdge = 20.0)        {fEtaExpec_High               = EtaUpperEdge;}
  void   SetEtaExpec_NBins(const Int_t nEtabins = 2000)                  {fEtaExpec_NBins            = nEtabins;}

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
  void FillTreeMC();
  void FastGen();                           // Run over galice.root for Fastgen 2nd moments
  void GetExpecteds(AliESDtrack *track, Double_t closestPar[3]);
  void FillEventTree();

  //
  Int_t CountEmptyEvents(Int_t counterBin);  // Just count if there is empty events
  Int_t CacheTPCEventInformation();
  UInt_t SetCutBitsAndSomeTrackVariables(AliESDtrack *track);
  Bool_t CheckIfFromResonance(Int_t mcType, AliMCParticle *trackMCgen, Int_t trackIndex, Bool_t parInterest, Double_t ptot, Double_t eta, Double_t cent, Bool_t fillTree);
  Bool_t CheckIfFromAnyResonance(AliMCParticle *trackMCgen, Float_t etaLow, Float_t etaUp, Float_t pDown, Float_t pUp);
  Bool_t ApplyDCAcutIfNoITSPixel(AliESDtrack *track);
  Bool_t GetSystematicClassIndex(UInt_t cut,Int_t syst);
  // ---------------------------------------------------------------------------------
  //                                   Members
  // ---------------------------------------------------------------------------------

  AliPIDResponse   * fPIDResponse;            //! PID response object
  AliESDEvent      * fESD;                    //! ESD object
  TList            * fListHist;               //! list for histograms
  AliESDtrackCuts  * fESDtrackCuts;           //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit96;     //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit128;    //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit768;    //! basic cut variables
  AliESDtrackCuts  * fESDtrackCutsLoose;      //! basic cut variables for debugging
  AliESDtrackCuts  * fESDtrackCutsCleanSamp;  //! basic cut variables for clean pion and electron form V0s
  AliPIDCombined   * fPIDCombined;            //! combined PID object
  AliTPCdEdxInfo   * fTPCdEdxInfo;            //! detailed dEdx info
  AliStack         * fMCStack;                //! stack object to get Mc info
  AliAnalysisCuts  * fK0sPionCuts;            // filter for pions from K0s
  AliAnalysisCuts  * fLambdaProtonCuts;       // filter for protons from Lambda
  AliAnalysisCuts  * fLambdaPionCuts;         // filter for pions from Lambda
  AliAnalysisCuts  * fGammaElectronCuts;      // filter for electrons from gamma conversions
  const AliESDVertex * fVertex;               // primary vertex
//  AliESDtools      * fESDtool;                 // tools to calculate derived variables from the ESD

  TTreeSRedirector * fTreeSRedirector;        /// temp tree to dump output
  TTree            * fTreeMC;                 // tree for mc samples
  TTree            * fTreeCuts;               // tree to save all variables for control plots
  TTree            * fTreejetsEMCconst;       // tree for EMCal signal jet constituents
  TTree            * fTreejetsEMCBGconst;       // tree for EMCal background jet constituents
  TTree            * fTreejetsFJ;             // tree for fastjet signal jets
  TTree            * fTreejetsFJBG;           // tree for fastjet background jets
  TTree            * fTreejetsFJconst;          // tree for fastjet signal jet constituents
  TTree            * fTreejetsFJBGconst;          // tree for fastjet signal jet constituents
  TTree            * fTreejetResonance;          // tree with full acceptance filled with MC
  TTree            * fTreejetEvents;
  TTree            * fTreejetsEMC;            // tree for EMCal signal jets
  TTree            * fTreejetsEMCBG;            // tree for EMCal background jets
  TRandom3         fRandom;


  TString            fPeriodName;
  Int_t              fYear;
  Int_t              fPassIndex;
  UInt_t             fPileUpBit;
  TH1F             * fHistCent;               // helper histogram for TIdentity tree
  TH1F             * fHistPhi;
  TString           fChunkName;

  UInt_t            fTrackCutBits;           // integer which hold all cut variations as bits
  Int_t             fSystClass;
  Double_t          fEtaDown;
  Double_t          fEtaUp;
  Int_t             fNEtaBins;
  Int_t             fPercentageOfEvents;     // when only a fPercentageOfEvents is enough

  Bool_t            fRunOnGrid;              // flag if real data or MC is processed
  Bool_t            fSmallOut;              // flag for small output
  Bool_t            fMCtrue;                 // flag if real data or MC is processed
  Bool_t            fEventInfo;              // flag if event info and downscaled track tree is filled
  Bool_t            fDEdxCheck;              // flag to check only the dEdx performance
  Bool_t            fIncludeITS;             // decide whether to use ITS or not
  Bool_t            fFillIncTracks;          // switch whether to fill tracks tree
  Bool_t            fFillpTPC;               // switch whether to fill histos with pTPC
  Bool_t            fFillp;               // switch whether to fill histos with p
  Bool_t            fFillpT;               // switch whether to fill histos with pT
  Bool_t            fFillOnlyHists;          //
  Bool_t            fFillEffLookUpTable;     //
  Bool_t            fFillArmPodTree;         // switch whether to fill clean sample tree
  Bool_t            fFilljetsFJBGTree;         // switch whether to fill BG Jets FJ tree
  Bool_t            fFillJetsFJBGConst;        // switch whether to fill jetsEMCBG constituent tree
  Bool_t            fFilldscaledTree;         // switch whether to fill dscaled tree
  Bool_t            fFillFastJet;         // switch whether to fill FJ tree
  Bool_t            fFillEMCJet;         // switch whether to fill EMC tree
  Bool_t            fDoFastJet;         // switch whether to use FJ jets
  Bool_t            fDoEMCJet;         // switch whether to use EMC jets
  Bool_t            fFillJetsEMCConst;        // switch whether to fill jetsEMC constituent tree
  Bool_t            fFillJetsFJConst;        // switch whether to fill jetsFJ constituent tree
  Bool_t            fFillJetsEMCBG;        // switch whether to fill jetsEMCBG tree
  Bool_t            fFillJetsEMCBGConst;        // switch whether to fill jetsEMCBG constituent tree
  Double_t          fjetMinPtSub;            // minimium jet pt after subtraction to keep jet
  Double_t          fjetMinArea;            // minimium jet pt after subtraction to keep jet
  Float_t           fcent_min;            // minimium centrality cut
  Float_t           fcent_max;            // maximium centrality cut

  Bool_t            fRunFastSimulation;      // when running over galice.root do not fill other objects
  Bool_t            fFillDistributions;   // when running over galice.root do not fill other objects
  Bool_t            fFillTreeMC;
  Bool_t            fDefaultTrackCuts;
  Bool_t            fDefaultEventCuts;
  Int_t             fCorrectForMissCl;       // 0; defaults crows, 1; ncls used wo correction, 2; ncls used with correction
  Int_t             fUsePtCut;
  Int_t             fTrackOriginOnlyPrimary;
  Int_t             fRapidityType;
  Int_t             fSisterCheck;           // 0: reject the mother anyways, 1: if both girls are in acceptance rejet mother

  Bool_t            fIncludeTOF;             // Include TOF information to investigate the efficiency loss effects on observable
  Bool_t            fUseCouts;               // for debugging
  Int_t             fRunNumberForExpecteds;  // Run number in which to fill the expecteds tree
  Bool_t            fFillExpecteds;

  Int_t             fNSettings;
  Int_t             fNMomBins;               // number of mombins --> for 20MeV slice 150 and 10MeV 300
  Float_t           fMomDown;                // bottom limit for the momentum range (default 0.2)
  Float_t           fMomUp;                  // uppper limit for the momentum range (default 3.2)
  Float_t           fNdEdxBins;           // bin width for the dEdx histograms (default 2.5)
  Float_t           fDEdxUp;                 // bottom limit for dEdx histogram (default 20)
  Float_t           fDEdxDown;               // upper limit for dEdx histogram (default 1020)
  Float_t           fDEdxCleanUp;            // upper limit for dEdx histogram of clean kaons and electrons (default 140)
  Int_t             fMomExpec_NBins;         // number of mom bins for expecteds (default 2000)
  Float_t           fMomExpec_Low;           // bottom limit for the momentum range for expecteds (default 0.0)
  Float_t           fMomExpec_High;          // uppper limit for the momentum range for expecteds (default 20.0)
  Int_t             fEtaExpec_NBins;         // number of absolute eta bins for expecteds (default 9)
  Float_t           fEtaExpec_Low;           // bottom limit for the absolute eta range for expecteds (default 0.0)
  Float_t           fEtaExpec_High;          // uppper limit for the absolute eta range for expecteds (default 0.9)

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
  Int_t              fTPCFindable;            // number of findable clusters
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
  Bool_t             fIsMCPileup;

  AliJetContainer*   fJetContainer;
  AliJetContainer*   fbgJetContainer;
  Double_t           fJetPt;
  Double_t           fJetEta;
  Double_t           fJetPhi;
  Float_t            fjetRhoVal;
  Float_t            frhoFJ;
  Bool_t             fisGoodIncEvent;
  Bool_t             fhasAcceptedFJjet;
  Bool_t             fhasAcceptedEMCjet;
  Bool_t             fhasRealFJjet;
  Bool_t             fhasRealEMCjet;
  Int_t              fNumRealJets;
  Double_t           ftotalJetArea;
  Int_t              ftotalNumRealJets;
  Int_t              ftotalNumRealJetEvents;
  Int_t              ftotalNumIncEvents;

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

  //  Variables for systematic uncertainty checks
  //  B field configurations -->  use default settings and analyse the following set of runs
  //  ------------------------------------------------
  //  Field (++)  --> run interval is [137161, 138275]
  //  Field (--)  --> run interval is [138364, 139510]
  //  ------------------------------------------------
  Int_t              fSystCentEstimatetor;   // 0 --> "V0M"   ||| -1 -->  "TRK" ||| +1 --> "CL1"

  std::vector<float>  fetaDownArr;
  std::vector<float>  fetaUpArr;
  std::vector<float>  fcentDownArr;
  std::vector<float>  fcentUpArr;
  std::vector<float>  fpDownArr;
  std::vector<float>  fpUpArr;
  std::vector<float>  fxCentBins;
  std::vector<std::string> fResonances;
  std::vector<int>    fBaryons;
  //
  // control and QA histograms
  //
  //
  TH1F             * fHistEmptyEvent;            // control histogram for empty event
  TH1F             * fHistCentrality;            // control histogram for centrality
  TH1F             * fHistCentralityImpPar;      // control histogram for centrality
  TH1F             * fHistImpParam;              // control histogram for impact parameter
  TH1F             * fHistVertex;                // control histogram for vertexZ
  TH3F             * fHistIncTracks_dEdx;        // histogram for inclusive tracks dEdx all eta v ptpc
  TH3F             * fHistIncTracks_dEdx_p;      // histogram for inclusive tracks dEdx all eta v p
  TH3F             * fHistIncTracks_dEdx_pT;     // histogram for inclusive tracks dEdx all eta v pT
  TH2F             * fHistIncTracks_moms;        // histogram for inclusive tracks ptpc to pT
  TH2F             * fHistIncTracks_moms_p;        // histogram for inclusive tracks p to pT
  TH3F             * fHistIncTracks_kin;        // histogram for inclusive tracks dEdx all eta

  TH3F             * fHistJetTracks_dEdx;        // histogram for jet tracks dEdx all eta v ptpc
  TH3F             * fHistJetTracks_dEdx_p;        // histogram for jet tracks dEdx all eta v p
  TH3F             * fHistJetTracks_dEdx_pT;        // histogram for jet tracks dEdx all eta v pt
  TH2F             * fHistJetTracks_moms;        // histogram for jet tracks ptpc to pT
  TH2F             * fHistJetTracks_moms_p;        // histogram for jet tracks p to Pt
  TH3F             * fHistJetTracks_kin;         // histogram for jet tracks dEdx all eta
/*
  TH2D             * fHistIncTracks_mpi_small;        // histogram for inclusive tracks dEdx expected pion mean
  TH2D             * fHistIncTracks_spi_small;        // histogram for inclusive tracks dEdx expected pion sigma
  TH2D             * fHistIncTracks_mel_small;        // histogram for inclusive tracks dEdx expected electron mean
  TH2D             * fHistIncTracks_sel_small;        // histogram for inclusive tracks dEdx expected electron sigma
  TH2D             * fHistIncTracks_mka_small;        // histogram for inclusive tracks dEdx expected kaon mean
  TH2D             * fHistIncTracks_ska_small;        // histogram for inclusive tracks dEdx expected kaon sigma
  TH2D             * fHistIncTracks_mpr_small;        // histogram for inclusive tracks dEdx expected proton mean
  TH2D             * fHistIncTracks_spr_small;        // histogram for inclusive tracks dEdx expected proton sigma
*/
  TH3F             * fHistIncTracks_mpi;        // intermediate histogram for inclusive tracks dEdx expected pion mean
  TH3F             * fHistIncTracks_spi;        // intermediate histogram for inclusive tracks dEdx expected pion sigma
  TH3F             * fHistIncTracks_mel;        // intermediate histogram for inclusive tracks dEdx expected electron mean
  TH3F             * fHistIncTracks_sel;        // intermediate histogram for inclusive tracks dEdx expected electron sigma
  TH3F             * fHistIncTracks_mka;        // intermediate histogram for inclusive tracks dEdx expected kaon mean
  TH3F             * fHistIncTracks_ska;        // intermediate histogram for inclusive tracks dEdx expected kaon sigma
  TH3F             * fHistIncTracks_mpr;        // intermediate histogram for inclusive tracks dEdx expected proton mean
  TH3F             * fHistIncTracks_spr;        // intermediate histogram for inclusive tracks dEdx expected proton sigma

  TH2F             * fHistJet_ptsub_v_area;     // histogram for before any cuts, jet pt after bg subtraction vs jet area
  TH3F             * fHistJet_kin;     // histogram for jet ptsub, eta, phi after area cut
  TH2F             * fHistJet_moms;     // histogram for jet pt v jet ptsub after area cut
  //
  // Counters for Marian
  //
  TVectorF         * fEventInfo_CentralityEstimates;
  TGraph           * fEventInfo_LumiGraph;           // grap for the interaction rate info for a run
  static const char*  fEventInfo_centEstStr[];              //!centrality types

  AliEventCuts* fPileUpTightnessCut4;
  AliEventCuts* fPileUpTightnessCut3;
  AliEventCuts* fPileUpTightnessCut2;
  AliEventCuts* fPileUpTightnessCut1;
  Double_t fEffMatrixNSigmasTOF;

  ClassDef(AliAnalysisJetHadro, 11);

};

#endif
