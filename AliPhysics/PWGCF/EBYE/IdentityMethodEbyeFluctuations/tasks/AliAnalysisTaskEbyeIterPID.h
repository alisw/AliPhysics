#ifndef ALIANALYSISEBYERATIOS_H
#define ALIANALYSISEBYERATIOS_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// This analysis extracts pT-spectra of charged kaons, protons, and pions.  //
// It is based on particles identifation via the dE/dx signal of the TPC.   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

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
class fPIDCombined;


#include "AliAnalysisTaskSE.h"
#include "AliPIDCombined.h"
#include "AliTPCdEdxInfo.h"
#include "AliESDv0KineCuts.h"
#include "THnSparse.h"
#include "THn.h"
#include "TVectorF.h"
#include "TCutG.h"
#include "TTreeStream.h"
#include "AliESDv0Cuts.h"
#include "AliEventCuts.h"

// class AliAnalysisTaskPIDetaTreeElectrons : public AliAnalysisTaskPIDV0base {
class AliAnalysisTaskEbyeIterPID : public AliAnalysisTaskSE {
 public:

  AliEventCuts fEventCuts;     /// Event cuts

  // ---------------------------------------------------------------------------------
  //                           Constructor and Destructor
  // ---------------------------------------------------------------------------------

  AliAnalysisTaskEbyeIterPID(const char *name);
  AliAnalysisTaskEbyeIterPID();
  virtual ~AliAnalysisTaskEbyeIterPID();

  enum momentType {kPi=0,kKa=1,kPr=2,kPiPi=3,kKaKa=4,kPrPr=5,kPiKa=6,kPiPr=7,kKaPr=8,kLa=9,kLaLa=10,kCh=11,kChCh=12,kBa=13,kBaBa=14};

  enum momentTypeUnlike {
    kPiPosPiNeg=0,
    kKaPosKaNeg=1,
    kPrPosPrNeg=2,
    kPiPosKaNeg=3,
    kPiPosPrNeg=4,
    kKaPosPiNeg=5,
    kKaPosPrNeg=6,
    kPrPosPiNeg=7,
    kPrPosKaNeg=8,
    kLaPosLaNeg=9,
    kChPosChNeg=10,
    kBaPosBaNeg=11,
  };

  enum trackCutBit {
    kNCrossedRowsTPC60=0,
    kNCrossedRowsTPC80=1,
    kNCrossedRowsTPC100=2,
    kMaxChi2PerClusterTPC3=3,
    kMaxChi2PerClusterTPC4=4,
    kMaxChi2PerClusterTPC5=5,
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
    kTrackProbPiTPC=24,
    kTrackProbKaTPC=25,
    kTrackProbPrTPC=26,
    kTrackProbDeTPC=27,
    kTrackProbPiTOF=28,
    kTrackProbKaTOF=29,
    kTrackProbPrTOF=30,
    kTrackProbDeTOF=31,
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

  enum kNetMoments{
    kA=0,
    kB=1,
    kAA=2,
    kBB=3,
    kAB=4,
    kAAA=5,
    kBBB=6,
    kAAB=7,
    kBBA=8,
    kABBB=9,
    kAABB=10,
    kAAAB=11,
    kAAAA=12,
    kBBBB=13,
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
  virtual void   UserExec(Option_t *option);           // run over event-by-event and fill output objects
  virtual void   Terminate(Option_t *);                // run only once and terminate

  // ---------------------------------------------------------------------------------
  //                                    Settings
  // ---------------------------------------------------------------------------------

  void   SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  void   SetIsMCtrue(Bool_t isMCdata = kTRUE){fMCtrue = isMCdata;};
  void   Initialize();
  void   PrintNumInBinary(UInt_t num);

  // Some boolian settings
  void   SetRunOnGrid(const Bool_t ifRunOnGrid = kTRUE)               {fRunOnGrid           = ifRunOnGrid;}
  void   SetIncludeITScuts(const Bool_t ifITSCuts = kTRUE)            {fIncludeITS          = ifITSCuts;}
  void   SetFillArmPodTree(const Bool_t ifArmpodTree = kTRUE)         {fFillArmPodTree      = ifArmpodTree;}
  void   SetDeDxCheck(const Bool_t ifDeDxCheck = kFALSE)              {fDEdxCheck           = ifDeDxCheck;}
  void   SetEffMatrix(const Bool_t ifEffMatrix = kFALSE)              {fEffMatrix           = ifEffMatrix;}
  void   SetFillAllCutVariables(const Bool_t ifAllCuts = kFALSE)      {fFillTracks          = ifAllCuts;}
  void   SetFillOnlyHists(const Bool_t ifFillOnlyHists = kFALSE)      {fFillOnlyHists       = ifFillOnlyHists;}
  void   SetFillEffLookUpTable(const Bool_t ifEffLookUpTable = kFALSE){fFillEffLookUpTable  = ifEffLookUpTable;}
  void   SetFillHigherMomentsMCclosure(const Bool_t ifHigherMomentsMCclosure = kFALSE){fFillHigherMomentsMCclosure  = ifHigherMomentsMCclosure;}
  void   SetRunFastSimulation(const Bool_t ifFastSimul = kFALSE)      {fRunFastSimulation   = ifFastSimul;}
  void   SetRunFastHighMomentCal(const Bool_t ifFastHighMom = kFALSE) {fRunFastHighMomentCal= ifFastHighMom;}
  void   SetFillGenDistributions(const Bool_t ifGenDistributions = kFALSE) {fFillGenDistributions= ifGenDistributions;}
  void   SetFillTreeMC(const Bool_t ifTreeMC = kFALSE)                {fFillTreeMC= ifTreeMC;}

  void   SetDefaultTrackCuts(const Bool_t ifDefaultTrackCuts = kFALSE){fDefaultTrackCuts= ifDefaultTrackCuts;}
  void   SetDefaultEventCuts(const Bool_t ifDefaultEventCuts = kFALSE){fDefaultEventCuts= ifDefaultEventCuts;}
  void   SetFillNudynFastGen(const Bool_t ifNudynFastGen = kFALSE)    {fFillNudynFastGen= ifNudynFastGen;}
  void   SetUsePtCut(const Int_t ifUsePtCut = 1)                      {fUsePtCut            = ifUsePtCut;}
  void   SetTrackOriginType(const Int_t ifTrackOriginType = 0)        {fTrackOriginType     = ifTrackOriginType;}
  void   SetRapidityType(const Int_t ifRapidityType = 0)              {fRapidityType        = ifRapidityType;}
  void   SetFillDnchDeta(const Bool_t ifDnchDetaCal = kFALSE)         {fFillDnchDeta        = ifDnchDetaCal;}
  void   SetIncludeTOF(const Bool_t ifIncludeTOF = kFALSE)            {fIncludeTOF          = ifIncludeTOF;}
  void   SetUseThnSparse(const Bool_t ifUseThnSparse = kFALSE)        {fUseThnSparse        = ifUseThnSparse;}
  void   SetUseCouts(const Bool_t ifUseCouts = kFALSE)                {fUseCouts            = ifUseCouts;}
  void   SetWeakAndMaterial(const Bool_t ifWeakAndMaterial = kFALSE)  {fWeakAndMaterial     = ifWeakAndMaterial;}
  void   SetFillEventInfo(const Bool_t ifEventInfo = kFALSE)          {fEventInfo           = ifEventInfo;}
  void   SetPercentageOfEvents(const Int_t nPercentageOfEvents = 0)   {fPercentageOfEvents = nPercentageOfEvents;}
  //
  Bool_t GetRunOnGrid() const { return fRunOnGrid; }


  // Setters for the systematic uncertainty checks
  void   SetSystCentEstimator(const Int_t systCentEstimator = 0)  {fSystCentEstimatetor = systCentEstimator;}
  void   SetSystDCAxy(const Int_t systDCAxy = 0)                  {fSystDCAxy           = systDCAxy;}
  void   SetSystNCrossedRows(const Int_t systNCrossedRows = 0)    {fSystCrossedRows     = systNCrossedRows;}
  void   SetSystTPCChi2(const Int_t systTPCChi2 = 0)              {fSystChi2            = systTPCChi2;}
  void   SetSystVz(const Int_t systVz = 0)                        {fSystVz              = systVz;}

  // Setters for the eta momentum dEdx and centrality bins
  void   SetSampleDeDxUpperEdge(const Float_t dEdxCleanUp = 200.) {fDEdxCleanUp         = dEdxCleanUp;}
  void   SetDeDxBinWidth(const Float_t dEdxBinWidth = 2.5)        {fDEdxBinWidth        = dEdxBinWidth;}
  void   SetDeDxLowerEdge(const Float_t dEdxLowerEdge = 20.)      {fDEdxDown            = dEdxLowerEdge;}
  void   SetDeDxUpperEdge(const Float_t dEdxUpperEdge = 1020.)    {fDEdxUp              = dEdxUpperEdge;}

  void   SetEtaLowerEdge(const Float_t etaLowerEdge = -1.)        {fEtaDown             = etaLowerEdge;}
  void   SetEtaUpperEdge(const Float_t etaUpperEdge = 1.)         {fEtaUp               = etaUpperEdge;}
  void   SetNEtabins(const Int_t nEtaBins = 20)                   {fNEtaBins            = nEtaBins;}
  void   SetMomLowerEdge(const Float_t momLowerEdge = 0.)         {fMomDown             = momLowerEdge;}
  void   SetMomUpperEdge(const Float_t momUpperEdge = 12.)        {fMomUp               = momUpperEdge;}
  void   SetNMomBins(const Int_t nMombins = 600)                  {fNMomBins            = nMombins;}
  void   SetNGenprotonBins(const Int_t nGenprotonBins = 100)      {fGenprotonBins       = nGenprotonBins;}



  // Set the binning of centrality
  void SetCentralityBinning(const Int_t tmpCentbins, Float_t tmpfxCentBins[])
  {
    // Create the histograms to be used in the binning of eta, cent and momentum
    std::cout << " Info::marsland: !!!!!! Centrality binning is being set !!!!!!! " << std::endl;
    fHistCent =  new TH1F("fHistCent","Centrality Bins",tmpCentbins-1    ,tmpfxCentBins );
    fHistPhi  =  new TH1F("fHistPhi" ,"Phi Bins"       ,36               ,-TMath::Pi(), TMath::Pi());
    // ==========================================
    // prepare real data centrality bins
    fNCentbinsData = tmpCentbins;
    fNCentBinsMC   = tmpCentbins-1;
    fxCentBins = new Float_t[fNCentbinsData];
    for (Int_t i=0; i<fNCentbinsData; i++) fxCentBins[i] =  tmpfxCentBins[i];
    fcentDownArr = new Float_t[fNCentBinsMC];
    fcentUpArr   = new Float_t[fNCentBinsMC];
    for (Int_t i=0; i<fNCentbinsData-1; i++) fcentDownArr[i] =  tmpfxCentBins[i];
    for (Int_t i=1; i<fNCentbinsData; i++)   fcentUpArr[i-1] =  tmpfxCentBins[i];
  }

  void SetMCEtaScanArray(const Int_t tmpEtaBinsMC, Float_t tmpetaDownArr[], Float_t tmpetaUpArr[])
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetMCEtaScanArray is being set !!!!!!! " << std::endl;
    fNEtaWinBinsMC = tmpEtaBinsMC;
    fetaDownArr = new Float_t[fNEtaWinBinsMC];
    fetaUpArr   = new Float_t[fNEtaWinBinsMC];
    for (Int_t i=0; i<fNEtaWinBinsMC; i++) {
      fetaDownArr[i] =  tmpetaDownArr[i];
      fetaUpArr[i]   =  tmpetaUpArr[i];
    }
  }

  void SetMCResonanceArray(const Int_t tmpNRes, TString tmpResArr[])
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetMCResonanceArray is being set !!!!!!! " << std::endl;
    fNResBins = tmpNRes;
    fResonances = new TString[fNResBins];
    for (Int_t i=0; i<fNResBins; i++) fResonances[i] = tmpResArr[i];

  }

  void SetMCBaryonArray(const Int_t tmpNBar, Int_t tmpBarArr[])
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetMCBaryonArray is being set !!!!!!! " << std::endl;
    fNBarBins = tmpNBar;
    fBaryons = new Int_t[fNBarBins];
    for (Int_t i=0; i<fNBarBins; i++) fBaryons[i] = tmpBarArr[i];
  }

  void SetMCMomScanArray(const Int_t tmpMomBinsMC, Float_t tmppDownArr[], Float_t tmppUpArr[])
  {
    // set MC momentum values to scan
    std::cout << " Info::marsland: !!!!!! SetMCMomScanArray is being set !!!!!!! " << std::endl;
    fNMomBinsMC = tmpMomBinsMC;
    fpDownArr = new Float_t[fNMomBinsMC];
    fpUpArr   = new Float_t[fNMomBinsMC];
    for (Int_t i=0; i<fNMomBinsMC; i++) {
      fpDownArr[i] =  tmppDownArr[i];
      fpUpArr[i]   =  tmppUpArr[i];
    }
  }

  void SetLookUpTableFirstMoments(TTree *lookUpTree, Int_t partType, Float_t pArr[],Float_t centArr[],Float_t etaArr[],const Int_t tmpMomBinsMC, const Int_t tmpCentbins, const Int_t tmpEtaBinsMC)
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetLookUpTableFirstMoments is being set !!!!!!!   " << std::endl;
    //
    // fill arrays from lookup table
    TH1D *h=NULL, *h1=NULL;
    for (Int_t imom=0; imom<tmpMomBinsMC; imom++){
      for (Int_t icent=0; icent<tmpCentbins; icent++){
        for (Int_t ieta=0; ieta<tmpEtaBinsMC; ieta++){
          //
          // with resonances
          lookUpTree->Draw(Form("momentPos.fElements[%d]-momentNeg.fElements[%d]",partType,partType),Form("abs(etaUp-%f)<0.01&&abs(pDown-%f)<0.01&&abs(centDown-%f)<0.01",etaArr[ieta],pArr[imom],centArr[icent]),"goff");
          h= (TH1D*)lookUpTree->GetHistogram()->Clone(); h-> SetName("Res");
          if (partType==0)  fNetPiFirstMoments[0][imom][icent][ieta] = h->GetMean();
          if (partType==1)  fNetKaFirstMoments[0][imom][icent][ieta] = h->GetMean();
          if (partType==2)  fNetPrFirstMoments[0][imom][icent][ieta] = h->GetMean();
          if (partType==9)  fNetLaFirstMoments[0][imom][icent][ieta] = h->GetMean();
          if (partType==11) fNetChFirstMoments[0][imom][icent][ieta] = h->GetMean();
          delete h;
          //
          // without resonances
          lookUpTree->Draw(Form("noResmomentPos.fElements[%d]-noResmomentNeg.fElements[%d]",partType,partType),Form("abs(etaUp-%f)<0.01&&abs(pDown-%f)<0.01&&abs(centDown-%f)<0.01",etaArr[ieta],pArr[imom],centArr[icent]),"goff");
          h1= (TH1D*)lookUpTree->GetHistogram()->Clone(); h1-> SetName("noRes");
          if (partType==0)  fNetPiFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          if (partType==1)  fNetKaFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          if (partType==2)  fNetPrFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          if (partType==9)  fNetLaFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          if (partType==11) fNetChFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          delete h1;

        }
      }
    }

  }

  void SetLookUpTableEfficiencyCorrection(TTree *lookUpTree, Int_t partType, Float_t pDownArr[], Float_t pUpArr[], Float_t centArr[],Float_t etaArr[],const Int_t tmpMomBinsMC, const Int_t tmpCentbins, const Int_t tmpEtaBinsMC)
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetLookUpTableEfficiencyCorrection is being set !!!!!!!   " << std::endl;
    //
    // fill arrays from lookup table
    TH1D *hGenNet=NULL,   *hRecNet=NULL;
    TH1D *hGenCross=NULL, *hRecCross=NULL;
    TH1D *hGenPos=NULL,   *hRecPos=NULL;
    TH1D *hGenNeg=NULL,   *hRecNeg=NULL;
    for (Int_t imom=0; imom<tmpMomBinsMC; imom++){
      for (Int_t ieta=0; ieta<tmpEtaBinsMC; ieta++){
        for (Int_t icent=0; icent<tmpCentbins; icent++){
          //
          // ----------------------------
          // NET Particle cumulants
          // ----------------------------
          //
          TString etaMomCentCut = Form("abs(etaUp-%f)<0.01 && abs(pDown-%f)<0.01 && abs(pUp-%f)<0.01 && abs(centDown-%f)<0.01",etaArr[ieta],pDownArr[imom],pUpArr[imom],centArr[icent]);
          // generated net particles
          TString genStr = Form("momentPosGen.fElements[%d]-momentNegGen.fElements[%d]",partType,partType);
          lookUpTree->Draw(genStr,etaMomCentCut,"goff");
          hGenNet = (TH1D*)lookUpTree->GetHistogram()->Clone(); hGenNet-> SetName("netProtonGen");
          if (partType==0)  fNetPiFirstMomentsGen[imom][icent][ieta] = hGenNet->GetMean();
          if (partType==1)  fNetKaFirstMomentsGen[imom][icent][ieta] = hGenNet->GetMean();
          if (partType==2)  fNetPrFirstMomentsGen[imom][icent][ieta] = hGenNet->GetMean();
          //
          // reconstructed net particles
          TString recStr = Form("momentPosRec.fElements[%d]-momentNegRec.fElements[%d]",partType,partType);
          lookUpTree->Draw(recStr,etaMomCentCut,"goff");
          hRecNet= (TH1D*)lookUpTree->GetHistogram()->Clone(); hRecNet-> SetName("netProtonRec");
          if (partType==0)  fNetPiFirstMomentsRec[imom][icent][ieta] = hRecNet->GetMean();
          if (partType==1)  fNetKaFirstMomentsRec[imom][icent][ieta] = hRecNet->GetMean();
          if (partType==2)  fNetPrFirstMomentsRec[imom][icent][ieta] = hRecNet->GetMean();
          //
          // ----------------------------
          // Cross cumulants
          // ----------------------------
          //
          // generated net particles
          TString momentCrossGenStr = Form("momentCrossGen.fElements[%d]",partType);
          lookUpTree->Draw(momentCrossGenStr,etaMomCentCut,"goff");
          hGenCross = (TH1D*)lookUpTree->GetHistogram()->Clone(); hGenCross-> SetName("crossProtonGen");
          if (partType==0)  fCrossPiFirstMomentsGen[imom][icent][ieta] = hGenCross->GetMean();
          if (partType==1)  fCrossKaFirstMomentsGen[imom][icent][ieta] = hGenCross->GetMean();
          if (partType==2)  fCrossPrFirstMomentsGen[imom][icent][ieta] = hGenCross->GetMean();
          //
          // reconstructed net particles
          TString momentCrossRecStr = Form("momentCrossRec.fElements[%d]",partType);
          lookUpTree->Draw(momentCrossRecStr,etaMomCentCut,"goff");
          hRecCross= (TH1D*)lookUpTree->GetHistogram()->Clone(); hRecCross-> SetName("crossProtonRec");
          if (partType==0)  fCrossPiFirstMomentsRec[imom][icent][ieta] = hRecCross->GetMean();
          if (partType==1)  fCrossKaFirstMomentsRec[imom][icent][ieta] = hRecCross->GetMean();
          if (partType==2)  fCrossPrFirstMomentsRec[imom][icent][ieta] = hRecCross->GetMean();
          //
          // ----------------------------
          // Single  Particle cumulants
          // ----------------------------
          //
          // generated particles
          TString momentPosGenStr = Form("momentPosGen.fElements[%d]",partType);
          lookUpTree->Draw(momentPosGenStr,etaMomCentCut,"goff");
          hGenPos = (TH1D*)lookUpTree->GetHistogram()->Clone(); hGenPos-> SetName("protonPosGen");
          if (partType==0)  fPiFirstMomentsGen[0][imom][icent][ieta] = hGenPos->GetMean();
          if (partType==1)  fKaFirstMomentsGen[0][imom][icent][ieta] = hGenPos->GetMean();
          if (partType==2)  fPrFirstMomentsGen[0][imom][icent][ieta] = hGenPos->GetMean();
          //
          TString momentNegGenStr = Form("momentNegGen.fElements[%d]",partType);
          lookUpTree->Draw(momentNegGenStr,etaMomCentCut,"goff");
          hGenNeg = (TH1D*)lookUpTree->GetHistogram()->Clone(); hGenNeg-> SetName("protonNegGen");
          if (partType==0)  fPiFirstMomentsGen[1][imom][icent][ieta] = hGenNeg->GetMean();
          if (partType==1)  fKaFirstMomentsGen[1][imom][icent][ieta] = hGenNeg->GetMean();
          if (partType==2)  fPrFirstMomentsGen[1][imom][icent][ieta] = hGenNeg->GetMean();
          //
          // reconstruced particles
          TString momentPosRecStr = Form("momentPosRec.fElements[%d]",partType);
          lookUpTree->Draw(momentPosRecStr,etaMomCentCut,"goff");
          hRecPos = (TH1D*)lookUpTree->GetHistogram()->Clone(); hRecPos-> SetName("protonPosRec");
          if (partType==0)  fPiFirstMomentsRec[0][imom][icent][ieta] = hRecPos->GetMean();
          if (partType==1)  fKaFirstMomentsRec[0][imom][icent][ieta] = hRecPos->GetMean();
          if (partType==2)  fPrFirstMomentsRec[0][imom][icent][ieta] = hRecPos->GetMean();
          //
          // reconstruced negative particles
          TString momentNegRecStr = Form("momentNegRec.fElements[%d]",partType);
          lookUpTree->Draw(momentNegRecStr,etaMomCentCut,"goff");
          hRecNeg = (TH1D*)lookUpTree->GetHistogram()->Clone(); hRecNeg-> SetName("protonNegRec");
          if (partType==0)  fPiFirstMomentsRec[1][imom][icent][ieta] = hRecNeg->GetMean();
          if (partType==1)  fKaFirstMomentsRec[1][imom][icent][ieta] = hRecNeg->GetMean();
          if (partType==2)  fPrFirstMomentsRec[1][imom][icent][ieta] = hRecNeg->GetMean();
          //
          // delete pointers to tmp histograms
          delete hGenNet;
          delete hRecNet;
          delete hGenPos;
          delete hGenNeg;
          delete hRecPos;
          delete hRecNeg;
          delete hGenCross;
          delete hRecCross;


          //

        }
      }
    }

  }

private:

  AliAnalysisTaskEbyeIterPID(const AliAnalysisTaskEbyeIterPID&);
  AliAnalysisTaskEbyeIterPID& operator=(const AliAnalysisTaskEbyeIterPID&);

  // ---------------------------------------------------------------------------------
  //                                   Functions
  // ---------------------------------------------------------------------------------

  void FillTPCdEdxReal();                   // Main function to fill all info + TIden
  void FillTPCdEdxCheck();                  // Quick check for the TPC dEdx
  void FillMCFull();                     // Fill all info + TIdenMC from MC to do MC closure test
  void FillTreeMC();
  void FillMCFull_NetParticles();
  void FastGen();                           // Run over galice.root for Fastgen
  void FastGenHigherMoments();     // Run over galice.root for Fastgen and calculate higher moments
  void MCclosureHigherMoments();   // Calculate higher moments for REC and GEN
  void WeakAndMaterial();                   // Look full acceptance, weak decay and material
  void FillDnchDeta();                      // Fill dnch/deta values for each cent and eta bin
  void FillEffMatrix();            // Prepare efficiency matrix
  void FillCleanSamples();                    // Fill Clean Pions
  void SelectCleanSamplesFromV0s(AliESDv0 *v0, AliESDtrack *track0, AliESDtrack *track1);
  void SetSpecialV0Cuts(AliESDv0KineCuts* cuts);
  void BinLogAxis(TH1 *h);
  void CalculateEventVariables();
  void SetCutBitsAndSomeTrackVariables(AliESDtrack *track);
  void DumpDownScaledTree();
  void GetExpecteds(AliESDtrack *track, Float_t closestPar[3]);
  void DumpEventVariables();
  Bool_t ApplyDCAcutIfNoITSPixel(AliESDtrack *track);
  Bool_t GetSystematicClassIndex(UInt_t cut,Int_t syst);
  Int_t CountEmptyEvents(Int_t counterBin);  // Just count if there is empty events
  Int_t CacheTPCEventInformation();
  Bool_t CheckIfFromResonance(Int_t mcType, AliMCParticle *trackMCgen, Int_t trackIndex, Bool_t parInterest, Double_t ptot, Double_t eta, Double_t cent, Bool_t fillTree);
  Bool_t CheckIfFromAnyResonance(AliMCParticle *trackMCgen);
  void FillGenDistributions();

  // ---------------------------------------------------------------------------------
  //                                   Members
  // ---------------------------------------------------------------------------------

  AliPIDResponse   * fPIDResponse;            //! PID response object
  AliESDEvent      * fESD;                    //! ESD object
  TList            * fListHist;               //! list for histograms
  AliESDtrackCuts  * fESDtrackCuts;           //! basic cut variables
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

  TTree            * fArmPodTree;             // Tree for clean pion and proton selection
  TTreeSRedirector * fTreeSRedirector;        //! temp tree to dump output
  TTree            * fTreeMCFull;             // tree for reconstructed moments
  TTree            * fTreeMCgen;              // tree for reconstructed moments
  TTree            * fTreeDnchDeta;           // tree for dnch/deta calculation
  TTree            * fTreeMC;                 // tree for mc samples
  TTree            * fTreedEdxCheck;          // tree to check dEdx performance for a small data sample
  TTree            * fTreeCuts;               // tree to save all variables for control plots
  TTree            * fTreeMCFullAcc;          // tree with full acceptance filled with MC
  TTree            * fTreeResonance;          // tree with full acceptance filled with MC
  TTree            * fTreeMCgenMoms;          // tree with higher moment calculations
  TTree            * fTreeEvents;
  TTree            * fTreeDScaled;
  TTree            * fTreeMCEffCorr;

  TH1F             * fHistCent;               // helper histogram for TIdentity tree
  TH1F             * fHistPhi;
  TH1F             * fHistGenMult;
  TH1F             * fHistInvK0s;             // helper histogram for TIdentity tree
  TH1F             * fHistInvLambda;          // helper histogram for TIdentity tree
  TH1F             * fHistInvAntiLambda;      // helper histogram for TIdentity tree
  TH1F             * fHistInvPhoton;          // helper histogram for TIdentity tree
  //
  TH1F             * fHistPhiTPCcounterA;     // helper histogram for TIdentity tree
  TH1F             * fHistPhiTPCcounterC;     // helper histogram for TIdentity tree
  TH1F             * fHistPhiTPCcounterAITS;  // helper histogram for TIdentity tree
  TH1F             * fHistPhiTPCcounterCITS;  // helper histogram for TIdentity tree
  TH1F             * fHistPhiITScounterA;     // helper histogram for TIdentity tree
  TH1F             * fHistPhiITScounterC;     // helper histogram for TIdentity tree

  THnSparseF       * fHndEdx;                 // histogram which hold all dEdx info
  THnSparseF       * fHnExpected[20];         // histogram which hold all PIDresponse info
  THnSparseF       * fHnCleanKa;              // histogram which hold Clean Kaons
  THnSparseF       * fHnCleanDe;              // histogram which hold Clean Deuterons

  TString           fChunkName;

  UInt_t            fTrackCutBits;           // integer which hold all cut variations as bits
  Int_t             fSystClass;
  Double_t          fEtaDown;
  Double_t          fEtaUp;
  Int_t             fNEtaBins;
  Int_t             fPercentageOfEvents;     // when only a fPercentageOfEvents is enough

  Bool_t            fRunOnGrid;              // flag if real data or MC is processed
  Bool_t            fMCtrue;                 // flag if real data or MC is processed
  Bool_t            fEventInfo;              // flag if event info and downscaled track tree is filled
  Bool_t            fWeakAndMaterial;        // flag for the Weak and Material analysis
  Bool_t            fEffMatrix;              // flag for efficiency matrix filling
  Bool_t            fDEdxCheck;              // flag to check only the dEdx performance
  Bool_t            fIncludeITS;             // decide whether to use ITS or not
  Bool_t            fFillTracks;             // switch whether to fill all cut variables
  Bool_t            fFillOnlyHists;          //
  Bool_t            fFillEffLookUpTable;     //
  Bool_t            fFillHigherMomentsMCclosure;
  Bool_t            fFillArmPodTree;         // switch whether to fill clean sample tree
  Bool_t            fRunFastSimulation;      // when running over galice.root do not fill other objects
  Bool_t            fRunFastHighMomentCal;   // when running over galice.root do not fill other objects
  Bool_t            fFillGenDistributions;   // when running over galice.root do not fill other objects
  Bool_t            fFillTreeMC;
  Bool_t            fDefaultTrackCuts;
  Bool_t            fDefaultEventCuts;
  Bool_t            fFillNudynFastGen;
  Int_t             fUsePtCut;
  Int_t             fTrackOriginType;
  Int_t             fRapidityType;


  Bool_t            fFillDnchDeta;           // switch on calculation of the dncdeta for fastgens
  Bool_t            fIncludeTOF;             // Include TOF information to investigate the efficiency loss effects on observable
  Bool_t            fUseThnSparse;           // in case thnsparse is filled
  Bool_t            fUseCouts;               // for debugging

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
  Float_t           fNSigmasElTOF;           // TOF N sigma for Pion
  Float_t           fNSigmasPiTOF;           // TOF N sigma for Pion
  Float_t           fNSigmasKaTOF;           // TOF N sigma for Pion
  Float_t           fNSigmasPrTOF;           // TOF N sigma for Pion
  Float_t           fNSigmasDeTOF;           // TOF N sigma for Proton

  Float_t           fDEdxEl;                 // Expected Electron dEdx
  Float_t           fDEdxKa;                 // Expected Kaon dEdx
  Float_t           fDEdxPi;                 // Expected Pion dEdx
  Float_t           fDEdxPr;                 // Expected Proton dEdx
  Float_t           fDEdxDe;                 // Expected Proton dEdx

  Float_t           fSigmaEl;                // Expected Electron sigma
  Float_t           fSigmaKa;                // Expected Kaon sigma
  Float_t           fSigmaPi;                // Expected Pion sigma
  Float_t           fSigmaPr;                // Expected Proton sigma
  Float_t           fSigmaDe;                // Expected Proton sigma

  Float_t           fNSigmasElTPC;           // TOF N sigma for Electron
  Float_t           fNSigmasPiTPC;           // TOF N sigma for Pion
  Float_t           fNSigmasKaTPC;           // TOF N sigma for Kaon
  Float_t           fNSigmasPrTPC;           // TOF N sigma for Proton
  Float_t           fNSigmasDeTPC;           // TOF N sigma for Proton

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
  Float_t            fPhi;                    // azimut angle
  Int_t              fSign;                   // sign of the particle
  Int_t              fTPCShared;              // number of shared clusters
  Int_t              fNcl;                    // number of points used for dEdx

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

  // Cut variables
  Double_t fTrackProbElTPC;
  Double_t fTrackProbPiTPC;
  Double_t fTrackProbKaTPC;
  Double_t fTrackProbPrTPC;
  Bool_t fTrackProbDeTPC;
  Double_t fTrackProbElTOF;
  Double_t fTrackProbPiTOF;
  Double_t fTrackProbKaTOF;
  Double_t fTrackProbPrTOF;
  Bool_t fTrackProbDeTOF;

  Float_t fTrackTPCCrossedRows;
  Float_t fTrackChi2TPC;
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
  //  ***********************************************
  //  Field (++)  --> run interval is [137161, 138275]
  //  Field (--)  --> run interval is [138364, 139510]
  //  ***********************************************
  Int_t              fSystCentEstimatetor;   // 0 --> "V0M"   ||| -1 -->  "TRK" ||| +1 --> "CL1"
  Int_t              fSystCrossedRows;       // 0 -->  80     ||| -1 -->   60   ||| +1 -->  100
  Int_t              fSystDCAxy;             // 0 --> default ||| -1 --> -sigma ||| +1 --> +sigma
  Int_t              fSystChi2;              // 0 -->  4      ||| -1 -->    3   ||| +1 -->   5
  Int_t              fSystVz;                // 0 -->  10     ||| -1 -->    8   ||| +1 -->   12
  Float_t            fNetPiFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetKaFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetPrFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetLaFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetChFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]

  Float_t            fNetPiFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetKaFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetPrFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetPiFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetKaFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetPrFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]

  Float_t            fCrossPiFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossKaFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossPrFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossPiFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossKaFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossPrFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]

  Float_t            fPiFirstMomentsGen[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fKaFirstMomentsGen[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fPrFirstMomentsGen[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fPiFirstMomentsRec[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fKaFirstMomentsRec[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fPrFirstMomentsRec[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]


  Float_t            *fetaDownArr;           //[fNEtaWinBinsMC]
  Float_t            *fetaUpArr;             //[fNEtaWinBinsMC]
  Float_t            *fcentDownArr;          //[fNCentBinsMC]
  Float_t            *fcentUpArr;            //[fNCentBinsMC]
  Float_t            *fpDownArr;             //[fNMomBinsMC]
  Float_t            *fpUpArr;               //[fNMomBinsMC]
  Float_t            *fxCentBins;            //[fNCentbinsData]
  TString            *fResonances;           //[fNResBins]
  Int_t              *fBaryons;              //[fNBarBins]
  //
  // control and QA histograms
  //
  THnF             * fHistPosEffMatrixRec;       // histogram efficiency matrix --> reconstructed traks
  THnF             * fHistNegEffMatrixRec;       // histogram efficiency matrix --> generated traks
  THnF             * fHistPosEffMatrixGen;       // histogram efficiency matrix --> reconstructed pions
  THnF             * fHistNegEffMatrixGen;       // histogram efficiency matrix --> generated pions
  TH1F             * fHistEmptyEvent;         // control histogram for empty event
  TH1F             * fHistCentrality;         // control histogram for centrality
  TH1F             * fHistCentralityImpPar;         // control histogram for centrality
  TH1F             * fHistImpParam;           // control histogram for impact parameter
  TH1F             * fHistVertex;             // control histogram for vertexZ
  THnF             * fHistdEdxTPC;            // 5D hist of dEdx from all TPC
  TH2F             * fHistArmPod;             // control histogram for Armanteros Podolanski plot
  //
  // Counters for Marian
  //
  TVectorF         * fPhiTPCdcarA;  // track counter
  TVectorF         * fPhiTPCdcarC; // dedx info counter
  TVectorF         * fCacheTrackCounters;  // track counter
  TVectorF         * fCacheTrackdEdxRatio; // dedx info counter
  TVectorF         * fCacheTrackNcl;       // ncl counter
  TVectorF         * fCacheTrackChi2;      // chi2 counter
  TVectorF         * fCacheTrackMatchEff;  // matchEff counter
  TVectorF         * fCentralityEstimates;
  TGraph           * fLumiGraph;           // grap for the interaction rate info for a run
  TH1F             * fHisTPCVertexA;
  TH1F             * fHisTPCVertexC;
  TH1F             * fHisTPCVertexACut;
  TH1F             * fHisTPCVertexCCut;
  TH1F             * fHisTPCVertex;
  TVectorF         * fCacheTrackTPCCountersZ; // track counter with DCA z cut
  static const char*  centEstStr[];              //!centrality types



  ClassDef(AliAnalysisTaskEbyeIterPID, 4);

};

#endif
