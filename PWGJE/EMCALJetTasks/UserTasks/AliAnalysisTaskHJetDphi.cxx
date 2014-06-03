#include <TCanvas.h>
#include <TChain.h>
#include <TFormula.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile2D.h>
#include <THnSparse.h>
#include <TROOT.h>
#include <TTree.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include "TKey.h"
#include "TList.h"
#include "TSystem.h"
#include "AliFJWrapper.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskHJetDphi.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliAODMCParticle.h"
#include "AliNamedArrayI.h"
#include "AliNamedString.h"
#include "AliPicoTrack.h"
#include "AliAnalysisTaskFastEmbedding.h"
#include "AliEmcalJet.h"
#include "AliAODJet.h"
#include "AliAODJetEventBackground.h"
#include "AliAnalysisHelperJetTasks.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

ClassImp(AliAnalysisTaskHJetDphi)

const Double_t pi = TMath::Pi();
const Double_t kSector = pi/9;

//________________________________________________________________________
AliAnalysisTaskHJetDphi::AliAnalysisTaskHJetDphi() : 
  AliAnalysisTaskSE(), 
  fVerbosity(0), fIsEmbedding(kFALSE), fAnaType(0), fPeriod("lhc11h"), fCollisionSystem("PbPb"),
  fIsMC(kFALSE), fAnalyzeMCTruth(kFALSE), fMC(0), 
  fEvent(0x0), fESD(0x0), fAODIn(0x0), fAODOut(0x0), fAODExtension(0x0),
  fOfflineTrgMask(AliVEvent::kAny), fTriggerType(-1), fCentrality(-1), fMaxVtxZ(10),
  fEsdTrkCut(0x0), fEsdHybCut(0x0), fFilterMask(0), fRequireITSRefit(kTRUE), fRequireSharedClsCut(kTRUE),
  fIsInit(kFALSE), fNonStdFile(""), fMcParticleArrName(""), fMcParticleArray(0x0),  fMcParticlelMap(0x0), 
  fEmbTrkArrName(""), fEmbTrkArray(0x0), fTrackArrName(""), fTrackArray(0x0), 
  fTriggerTrkIndex(-1), fTriggerTrkPt(-1), fSwitchOnAvoidTpcHole(kFALSE), fAvoidTpcHole(0), fCutTPCBoundary(kFALSE), fDistToTPCBoundary(0.), 
  fMinTrkPt(0.15), fMaxTrkPt(1e4), fMinTrkEta(-0.9), fMaxTrkEta(0.9), fMinTrkPhi(0), fMaxTrkPhi(2*pi), 
  fRadius(0.4), fJetArrName(""), fPLJetArrName(""), fDLJetArrName(""), fJetArray(0x0), fPLJetArray(0x0), fDLJetArray(0x0),
  fRhoName(""), fRho(0x0), fRhoValue(0), fEvtBkg(0x0), fPtHardBinName(0x0), fPtHardBin(-1),
  fRunTrkQA(kFALSE), fRunJetQA(kFALSE), fRunSingleInclHJet(kFALSE),  fTTtype(0), fTTMinPt(9), fTTMaxPt(10), fJetPtMin(10), 
  fRunPLHJet(kFALSE), fRunDLHJet(kFALSE), fRunLeadTrkQA(kFALSE), fStudyKtEffects(kFALSE), fKtValue(0), fRandom(0), 
  fRunBkgFlow(kFALSE),
  fOutputList(0x0), fhEventStat(0x0), fhNTrials(0x0), fhPtHardBins(0x0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  //DefineOutput(1, TList::Class());
  // Output slot #0 id reserved by the base class for AOD

  for(Int_t i=0; i<4; i++)
    {
      fhVtxZ[i]                  = 0x0;
      fhCentrality[i]            = 0x0;
      fhRhoVsCent[i]             = 0x0;

      fhTrkPt[i]                 = 0x0;
      fhTrkQA[i]                 = 0x0;
      fhTrkPtRes[i]              = 0x0;
      fhTrkPhiRes[i]             = 0x0;

      fhNumberOfTT[i]            = 0x0;
      for(Int_t j=0; j<3; j++)
	{
	  fhJetPt[i][j]          = 0x0;
	  fhJetArea[i][j]        = 0x0;
	  fhJetQA[i][j]          = 0x0;
	  
	  fhTTPt[i][j]           = 0x0;
	  fHJetPhiCorr[i][j]     = 0x0;
	}
      fHJetPhiCorrUp[i]          = 0x0;
      fHJetPhiCorrDown[i]        = 0x0;

      fhLeadTrkQA[i]             = 0x0;
      fhKtEffects[i]             = 0x0;
    }

  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
}

//________________________________________________________________________
AliAnalysisTaskHJetDphi::AliAnalysisTaskHJetDphi(const char *name) : 
  AliAnalysisTaskSE(name), 
  fVerbosity(0), fIsEmbedding(kFALSE), fAnaType(0), fPeriod("lhc11h"), fCollisionSystem("PbPb"),
  fIsMC(kFALSE), fAnalyzeMCTruth(kFALSE), fMC(0), 
  fEvent(0x0), fESD(0x0), fAODIn(0x0), fAODOut(0x0), fAODExtension(0x0),
  fOfflineTrgMask(AliVEvent::kAny), fTriggerType(-1), fCentrality(-1), fMaxVtxZ(10),
  fEsdTrkCut(0x0), fEsdHybCut(0x0), fFilterMask(0), fRequireITSRefit(kTRUE), fRequireSharedClsCut(kTRUE),
  fIsInit(kFALSE), fNonStdFile(""), fMcParticleArrName(""), fMcParticleArray(0x0),  fMcParticlelMap(0x0), 
  fEmbTrkArrName(""), fEmbTrkArray(0x0), fTrackArrName(""), fTrackArray(0x0), 
  fTriggerTrkIndex(-1), fTriggerTrkPt(-1), fSwitchOnAvoidTpcHole(kFALSE), fAvoidTpcHole(0), fCutTPCBoundary(kFALSE), fDistToTPCBoundary(0.), 
  fMinTrkPt(0.15), fMaxTrkPt(1e4), fMinTrkEta(-0.9), fMaxTrkEta(0.9), fMinTrkPhi(0), fMaxTrkPhi(2*pi), 
  fRadius(0.4), fJetArrName(""), fPLJetArrName(""), fDLJetArrName(""), fJetArray(0x0), fPLJetArray(0x0), fDLJetArray(0x0),
  fRhoName(""), fRho(0x0), fRhoValue(0), fEvtBkg(0x0), fPtHardBinName(0x0), fPtHardBin(-1),
  fRunTrkQA(kFALSE), fRunJetQA(kFALSE), fRunSingleInclHJet(kFALSE),  fTTtype(0), fTTMinPt(9), fTTMaxPt(10), fJetPtMin(10), 
  fRunPLHJet(kFALSE), fRunDLHJet(kFALSE), fRunLeadTrkQA(kFALSE), fStudyKtEffects(kFALSE), fKtValue(0), fRandom(0), 
  fRunBkgFlow(kFALSE),
  fOutputList(0x0), fhEventStat(0x0), fhNTrials(0x0), fhPtHardBins(0x0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  // Output slot #0 id reserved by the base class for AOD

  for(Int_t i=0; i<4; i++)
    {
      fhVtxZ[i]                  = 0x0;
      fhCentrality[i]            = 0x0;
      fhRhoVsCent[i]             = 0x0;

      fhTrkPt[i]                 = 0x0;
      fhTrkQA[i]                 = 0x0;
      fhTrkPtRes[i]              = 0x0;
      fhTrkPhiRes[i]             = 0x0;

      fhNumberOfTT[i]            = 0x0;
      for(Int_t j=0; j<3; j++)
	{
	  fhJetPt[i][j]          = 0x0;
	  fhJetArea[i][j]        = 0x0;
	  fhJetQA[i][j]          = 0x0;
	  
	  fhTTPt[i][j]           = 0x0;
	  fHJetPhiCorr[i][j]     = 0x0;
	}

      fHJetPhiCorrUp[i]          = 0x0;
      fHJetPhiCorrDown[i]        = 0x0;

      fhLeadTrkQA[i]             = 0x0;
      fhKtEffects[i]             = 0x0;
    }

  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
}


//________________________________________________________________________
AliAnalysisTaskHJetDphi::~AliAnalysisTaskHJetDphi()
{
  //Destructor
  if(fRandom)      delete fRandom;
  if(fEsdTrkCut)   delete fEsdTrkCut;
  if(fEsdHybCut)   delete fEsdHybCut;
  if(fOutputList) { fOutputList->Delete(); delete fOutputList;}
}

//________________________________________________________________________
void AliAnalysisTaskHJetDphi::UserCreateOutputObjects()
{
  // Create histograms

  const Int_t nTrkPtBins = 100;
  const Float_t lowTrkPtBin=0, upTrkPtBin=100;
  const Int_t nJetPtBins = 300;
  const Float_t lowJetPtBin=-100, upJetPtBin=200;

  // track QA
  const Int_t dimTrkqa = 4;
  const Int_t nBinsTrkqa[dimTrkqa]     = {nTrkPtBins/5,  36,  40, 10};
  const Double_t lowBinTrkqa[dimTrkqa] = {lowTrkPtBin,   0,  -1,  0};
  const Double_t hiBinTrkqa[dimTrkqa]  = {upTrkPtBin,    360, 1,  10};

  const Int_t dimTrkRes = 5;
  const Int_t nBinsTrkRes[dimTrkRes]     = {nTrkPtBins,    50,  50,  3, 10};
  const Double_t lowBinTrkRes[dimTrkRes] = {lowTrkPtBin,   0,   0,   0,  0};
  const Double_t hiBinTrkRes[dimTrkRes]  = {upTrkPtBin,    0.5, 0.5, 3,  10};

  const Int_t dimPhiRes = 4;
  const Int_t nBinsPhiRes[dimPhiRes]     = {nTrkPtBins,    200,     3, 10};
  const Double_t lowBinPhiRes[dimPhiRes] = {lowTrkPtBin,   -0.00995, 0, 0};
  const Double_t hiBinPhiRes[dimPhiRes]  = {upTrkPtBin,    0.01005,  3, 10};

  // jet QA
  const Int_t dimJetpt = 4;
  const Int_t nBinsJetpt[dimJetpt]     = {nJetPtBins,    300, 10, 10};
  const Double_t lowBinJetpt[dimJetpt] = {lowJetPtBin,   0,   0,  0};
  const Double_t hiBinJetpt[dimJetpt]  = {upJetPtBin,    300, 10, 10};

  const Int_t dimJetA = 4;
  const Int_t nBinsJetA[dimJetA]     = {nJetPtBins,    100, 10, 10};
  const Double_t lowBinJetA[dimJetA] = {lowJetPtBin,   0,   0,  0};
  const Double_t hiBinJetA[dimJetA]  = {upJetPtBin,    1,   10, 10};

  const Int_t dimJetqa = 7;
  const Int_t nBinsJetqa[dimJetqa]     = {nJetPtBins/5, 36,  24,  6,   100, 10, 11};
  const Double_t lowBinJetqa[dimJetqa] = {lowJetPtBin,   0,  -0.6, 0,   0,   0,  0};
  const Double_t hiBinJetqa[dimJetqa]  = {upJetPtBin,    360, 0.6, 1.2, 500, 10, 11};

  // h-jet analysis
  const Int_t dimTT = 4;
  const Int_t nBinsTT[dimTT]     = {nTrkPtBins,  10,  11, 10};
  const Double_t lowBinTT[dimTT] = {lowTrkPtBin, 0,   0,   0};
  const Double_t hiBinTT[dimTT]  = {upTrkPtBin,  100, 11,  10}; 
  
  const Int_t dimCor = 8;
  const Int_t nBinsCor[dimCor]     = {nTrkPtBins, nJetPtBins,  140,     6,   10, 40,    11, 10};
  const Double_t lowBinCor[dimCor] = {lowTrkPtBin,lowJetPtBin, pi-4.95, 0,   0,  -1.95, 0,  0};
  const Double_t hiBinCor[dimCor]  = {upTrkPtBin, upJetPtBin,  pi+2.05, 1.2, 100, 2.05, 11, 10};

  // Leading track QA
  const Int_t dimLeadTrkqa = 5;
  const Int_t nBinsLeadTrkqa[dimLeadTrkqa]     = {nTrkPtBins,  200,  80,   55, 10};
  const Double_t lowBinLeadTrkqa[dimLeadTrkqa] = {lowTrkPtBin,   0,  -0.4, 0,   0};
  const Double_t hiBinLeadTrkqa[dimLeadTrkqa]  = {upTrkPtBin,  200,  0.4,  1.1, 100};

  // kt effects
  const Int_t dimKt = 5;
  const Int_t nBinsKt[dimKt]     = {nTrkPtBins,   20,  81,    10, 11};
  const Double_t lowBinKt[dimKt] = {lowTrkPtBin,   0,  -40.5, 0,  0};
  const Double_t hiBinKt[dimKt]  = {upTrkPtBin,  200,  40.5,  10, 11};

   // Called once
   OpenFile(1);
   if(!fOutputList) fOutputList = new TList;
   fOutputList->SetOwner(kTRUE);
 
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   if(fAnaType==1) TH1::AddDirectory(kFALSE);

  fhEventStat = new TH1F("fhEventStat","Event statistics for jet analysis",8,0,8);
  fhEventStat->GetXaxis()->SetBinLabel(1,"All");
  fhEventStat->GetXaxis()->SetBinLabel(2,"PS");
  fhEventStat->GetXaxis()->SetBinLabel(3,"Vtx");
  fhEventStat->GetXaxis()->SetBinLabel(4,"Vtx+10cm");
  fhEventStat->GetXaxis()->SetBinLabel(5,"kMB");
  fhEventStat->GetXaxis()->SetBinLabel(6,"kCentral");
  fhEventStat->GetXaxis()->SetBinLabel(7,"kSemiCentral");
  fhEventStat->GetXaxis()->SetBinLabel(8,"kJetService");
  fOutputList->Add(fhEventStat);

  fhNTrials = new TH1F("fhNTrials","Number of trials",1,0,1);
  fOutputList->Add(fhNTrials);

  fhPtHardBins = new TH1F("fhPtHardBins","Number of events in each pT hard bin",11,0,11);
  fOutputList->Add(fhPtHardBins);

  const char *triggerName[4] = {"kMB","kEGA","kEJE","MC"};
  const char *dataType[3]    = {"", "_DL","_PL"};
  const char *dataName[3]    = {"Data","DL","PL"};
  
  Double_t newbins[7] = {0,0.07,0.2,0.4,0.6,0.8,1};

  for(Int_t i=0; i<4; i++)
    {
      if( fAnalyzeMCTruth )
	{
	  if(i!=3) continue;
	}
      else
	{
	  if( fPeriod=="lhc11a" && i>1 ) continue;
	  if( fPeriod=="lhc10h" && i!=0 ) continue;
	  if( fPeriod=="lhc11h" && i!=0 ) continue;
	  if( fPeriod.Contains("lhc12a15a") && i!=0 ) continue;
	  if( fPeriod.Contains("lhc12a15e") && i!=0 ) continue;
	}
 
      fhVtxZ[i] = new TH1F(Form("%s_fhVtxZ",triggerName[i]),Form("%s: z distribution of event vertexz;z(cm)",triggerName[i]),400,-20,20);
      fOutputList->Add(fhVtxZ[i]);

      fhCentrality[i] = new TH1F(Form("%s_fhCentrality",triggerName[i]),Form("%s: Event centrality;centrality",triggerName[i]),100,0,100);
      fOutputList->Add(fhCentrality[i]);

      fhRhoVsCent[i] = new TH2F(Form("%s_fhRhoVsCent",triggerName[i]),Form("%s: Rho vs centrality (R=%1.1f);centrality;Rho",triggerName[i],fRadius),100,0,100,300,0,300);
      fOutputList->Add(fhRhoVsCent[i]);

      if(fRunTrkQA)
	{
	  fhTrkPt[i] = new TH2F(Form("%s_fhTrkPt",triggerName[i]),Form("%s: Track p_{T} vs centrality;p_{T}^{track} (GeV/c);Centrality",triggerName[i]),nTrkPtBins,lowTrkPtBin,upTrkPtBin,10,0,100);
	  fOutputList->Add(fhTrkPt[i]);

	  fhTrkQA[i] = new THnSparseF(Form("%s_fhTrkQA",triggerName[i]),Form("%s: track p_{T} vs #phi vs #eta vs centrality;p_{T,track} (GeV/c);#phi;#eta;centrality",triggerName[i]),dimTrkqa,nBinsTrkqa,lowBinTrkqa,hiBinTrkqa);
	  fOutputList->Add(fhTrkQA[i]);

	  fhTrkPtRes[i] = new THnSparseF(Form("%s_fhTrkPtRes",triggerName[i]),Form("%s: track p_{T} vs resolution vs (p_{T}^{gen}-p_{T}^{rec})/p_{T}^{gen} vs type vs centrality;p_{T,track} (GeV/c);#sigma(p_{T})/p_{T};type;centrality",triggerName[i]),dimTrkRes,nBinsTrkRes,lowBinTrkRes,hiBinTrkRes);
	  fOutputList->Add(fhTrkPtRes[i]);

	  fhTrkPhiRes[i] = new THnSparseF(Form("%s_fhTrkPhiRes",triggerName[i]),Form("%s: track p_{T} vs #varphi^{gen}-#varphi^{rec} vs type vs centrality;p_{T,track} (GeV/c);#Delta#varphi;type;centrality",triggerName[i]),dimPhiRes,nBinsPhiRes,lowBinPhiRes,hiBinPhiRes);
	  fOutputList->Add(fhTrkPhiRes[i]);
	}

      for(Int_t j=0; j<3; j++)
	{
	  if(!fRunDLHJet && j==1) continue;
	  if(!fRunPLHJet && j==2) continue;
	  if(fRunJetQA)
	    {
	      fhJetPt[i][j] = new THnSparseF(Form("%s_fhJetPt%s",triggerName[i],dataType[j]),Form("%s-%s: jet p_{T} vs raw jet p_{T} vs centrality vs pt hard bin (R=%1.1f);p_{T,jet}^{ch} (GeV/c);p_{T,jet}^{raw} (GeV/c);centrality;ptHardBin",dataName[j],triggerName[i],fRadius),dimJetpt,nBinsJetpt,lowBinJetpt,hiBinJetpt);
	      fOutputList->Add(fhJetPt[i][j]);

	      fhJetArea[i][j] = new THnSparseF(Form("%s_fhJetArea%s",triggerName[i],dataType[j]),Form("%s-%s: jet p_{T} vs area vs centrality vs pt hard bin (R=%1.1f);p_{T,jet}^{ch} (GeV/c);area;centrality;ptHardBin",dataName[j],triggerName[i],fRadius),dimJetA,nBinsJetA,lowBinJetA,hiBinJetA);
	      fOutputList->Add(fhJetArea[i][j]);

	      fhJetQA[i][j] = new THnSparseF(Form("%s_fhJetQA%s",triggerName[i],dataType[j]),Form("%s-%s: jet p_{T} vs #phi vs #eta vs area vs # of constituents vs centrality vs pt hard bin (R=%1.1f);p_{T,jet}^{ch} (GeV/c);#phi;#eta;area;# of constituents;centrality;ptHardBin",dataName[j],triggerName[i],fRadius),dimJetqa,nBinsJetqa,lowBinJetqa,hiBinJetqa);
	      fhJetQA[i][j]->SetBinEdges(3,newbins);
	      fOutputList->Add(fhJetQA[i][j]);
	    }

	  if(fRunSingleInclHJet)
	    {
	      fhTTPt[i][j] = new THnSparseF(Form("%s_fhTTPt%s",triggerName[i],dataType[j]),Form("%s-%s: TT p_{T} vs centrality vs pT hard bin;p_{T,TT}^{ch} (GeV/c);centrality;ptHardBin",dataName[j],triggerName[i]),dimTT,nBinsTT,lowBinTT,hiBinTT);
	      fOutputList->Add(fhTTPt[i][j]);

	      fHJetPhiCorr[i][j] = new THnSparseF(Form("%s_fHJetPhiCorr%s",triggerName[i],dataType[j]),Form("%s-%s: single inclusive TT p_{T} vs recoil jet p_{T} vs #Delta#varphi vs area vs centrality vs #Delta#eta vs pt hard bin vs event bin (R=%1.1f);TT p_{T} (GeV/c);p_{T,jet}^{ch} (GeV/c);#Delta#varphi;area;centrality;#Delta#eta;ptHardBin;EventBin",dataName[j],triggerName[i],fRadius),dimCor,nBinsCor,lowBinCor,hiBinCor);
	      fHJetPhiCorr[i][j]->SetBinEdges(3,newbins);
	      fOutputList->Add(fHJetPhiCorr[i][j]);
	    }
	}


      if(fRunSingleInclHJet)
	{
	  fhNumberOfTT[i] = new TH1F(Form("%s_fhNumberOfTT",triggerName[i]), Form("%s: number of TT",triggerName[i]),6,0,6);
	  fOutputList->Add(fhNumberOfTT[i]);


	  if(fRunBkgFlow)
	    {
	      fHJetPhiCorrUp[i] = new THnSparseF(Form("%s_fHJetPhiCorrUp",triggerName[i]),Form("%s: single inclusive TT p_{T} vs recoil jet p_{T} vs #Delta#varphi vs area vs centrality vs #Delta#eta vs pt hard bin vs event bin (R=%1.1f);TT p_{T} (GeV/c);p_{T,jet}^{ch} (GeV/c);#Delta#varphi;area;centrality;#Delta#eta;ptHardBin;EventBin",triggerName[i],fRadius),dimCor,nBinsCor,lowBinCor,hiBinCor);
	      fOutputList->Add(fHJetPhiCorrUp[i]);
	      
	      fHJetPhiCorrDown[i] = new THnSparseF(Form("%s_fHJetPhiCorrDown",triggerName[i]),Form("%s: single inclusive TT p_{T} vs recoil jet p_{T} vs #Delta#varphi vs area vs centrality vs #Delta#eta vs pt hard bin vs event bin (R=%1.1f);TT p_{T} (GeV/c);p_{T,jet}^{ch} (GeV/c);#Delta#varphi;area;centrality;#Delta#eta;ptHardBin;EventBin",triggerName[i],fRadius),dimCor,nBinsCor,lowBinCor,hiBinCor);
	      fOutputList->Add(fHJetPhiCorrDown[i]);
	    }
	}

      if(fRunLeadTrkQA)
	{
	  fhLeadTrkQA[i] = new THnSparseF(Form("%s_fhLeadTrkQA",triggerName[i]),Form("%s: p_{T,trk}^{leading} vs p_{T,jet}^{full} vs #Delta#varphi vs z vs centrality;p_{T,trk}^{leading} (GeV/c); p_{T,jet}^{full} (GeV/c);#Delta#varphi;z;centrality",triggerName[i]),dimLeadTrkqa,nBinsLeadTrkqa,lowBinLeadTrkqa,hiBinLeadTrkqa);
	  fOutputList->Add(fhLeadTrkQA[i]);
	}

      if(fStudyKtEffects)
	{
	  fhKtEffects[i] = new THnSparseF(Form("%s_fhKtEffects",triggerName[i]),Form("%s: TT p_{T} vs recoil jet p_{T} vs k_{t} vs centrality vs pt hard bin (R=%1.1f);TT p_{T} (GeV/c);p_{T,jet}^{ch} (GeV/c);k_{t} (GeV/c);centrality;ptHardBin",triggerName[i],fRadius),dimKt,nBinsKt,lowBinKt,hiBinKt);
	  fOutputList->Add(fhKtEffects[i]);
	}
    }

  //error calculation in THnSparse
  Int_t nObj = fOutputList->GetEntries();
  for(Int_t i=0; i<nObj; i++)
    {
      TObject *obj = (TObject*) fOutputList->At(i);
      if (obj->IsA()->InheritsFrom( "THnSparse" ))
      	{
      	  THnSparseF *hn = (THnSparseF*)obj;
      	  hn->Sumw2();
      	}
    }

  if(fRunTrkQA)
    {
      fEsdTrkCut = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
      fEsdTrkCut->SetMaxDCAToVertexXY(2.4);
      fEsdTrkCut->SetMaxDCAToVertexZ(3.2);
      fEsdTrkCut->SetDCAToVertex2D(kTRUE);
      fEsdTrkCut->SetMaxChi2TPCConstrainedGlobal(36);
      fEsdTrkCut->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

      fEsdHybCut = new AliESDtrackCuts(*fEsdTrkCut);
      fEsdHybCut->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
    }

  fRandom = new TRandom3(0);

  PrintConfig();

  if(fAnaType==1) TH1::AddDirectory(oldStatus);
  PostData(1, fOutputList);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetDphi::UserNotify()
{
  AliInfo("User Nofity");

  Int_t runNumber = InputEvent()->GetRunNumber();

  fAvoidTpcHole = 0;
  if(fSwitchOnAvoidTpcHole)
    {
      Int_t runs_iroc[28] = {169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309};
      Int_t runs_oroc[23] = {169591, 169590, 169588, 169587, 169586, 169584, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169418, 169099, 169040, 169045, 169044};

      for(Int_t i=0; i<28; i++)
	{
	  if(runNumber==runs_iroc[i])
	    {
	      fAvoidTpcHole = 1;
	      break;
	    }
	}
      for(Int_t i=0; i<23; i++)
	{
	  if(runNumber==runs_oroc[i])
	    {
	      fAvoidTpcHole = 2;
	      break;
	    }
	}
    }

  if(fIsMC)
    {
      TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
      TFile *currFile = tree->GetCurrentFile();
      TString fileName(currFile->GetName());
      if(fileName.Contains("root_archive.zip#"))
	{
	  Ssiz_t pos = fileName.Index("#",0,TString::kExact);
	  fileName.Replace(pos+1,20,"");
	}
      else
	{
	  fileName.ReplaceAll(gSystem->BaseName(fileName.Data()),"");
	}

      TFile *fxsec = TFile::Open(Form("%s%s",fileName.Data(),"pyxsec_hists.root"),"read");
      if(fxsec)
	{
	  TKey *key = (TKey*)fxsec->GetListOfKeys()->At(0);
	  if(key)
	    {
	      TList *list = dynamic_cast<TList*>(key->ReadObj());
	      if(list)
		{
		  fhNTrials->Fill(0.5, ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1));
		}
	      else
		return kFALSE;
	    }
	  else
	    return kFALSE;
	}
      else
	{
	  fxsec = TFile::Open(Form("%s%s",fileName.Data(),"pyxsec.root"),"read");
	  TTree *xtree = (TTree*)fxsec->Get("Xsection");
	  if(xtree)
	    {
	      UInt_t ntrials = 0;
	      xtree->SetBranchAddress("ntrials",&ntrials);
	      xtree->GetEntry(0);
	      fhNTrials->Fill(0.5, ntrials);
	    }
	  else
	    return kFALSE;
	}
      fxsec->Close();
    }
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskHJetDphi::UserExec(Option_t *) 
{  
  // Main loop, called for each event.

  fTriggerType = -1;

  // Get pointers to input events
  fEvent = InputEvent();
  if (!fEvent) 
    {
      AliError("Input event not available");
      return;
    }

  if(fIsMC)  
    {
      fMC  = MCEvent();
      if (!fMC)
	{
	  AliError("MC event available");
	  return;
	}
    }

  fESD=dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) 
    {
      fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
    }
  else
    {
      if(fAnaType==1) fEvent = AODEvent();
    }
  if(fAnaType==1)
    {
      fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());
      if(fNonStdFile.Length()!=0)
	{
	  // case that we have an AOD extension we need can fetch the jets from the extended output
	  AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
	  fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
	  if(!fAODExtension)
	    {
	      if(fVerbosity>1) printf("AODExtension not found for %s",fNonStdFile.Data());
	    }
	}
    }

  fhEventStat->Fill(0.5);

  if(fVerbosity>1)
    {
      TList *list = 0x0;
      if(fAnaType==0) list = fEvent->GetList();
      else            list = fAODOut->GetList();
      for(Int_t i=0; i<list->GetEntries(); i++)
	{
	  TObject *obj = (TObject*)list->At(i);
	  cout<<i<<": "<<obj->GetName()<<" : "<<obj->ClassName()<<endl;
	}
    }
  // Retrieve arraies from memory
  if(!fIsInit) fIsInit = RetrieveArraies();
  if(!fIsInit) return;

  if(fIsEmbedding)
    {
      if(fAnaType==0)
	{
	  TString fileName = fPtHardBinName->GetString();
	  fileName.Remove(0,51);
	  fileName.Remove(fileName.Index("/"));
	  fPtHardBin = fileName.Atoi();
	}
      if(fAnaType==1)
	{
	  Double_t pthard = AliAnalysisTaskFastEmbedding::GetPtHard();
	  fPtHardBin = GetPtHardBin(pthard);
	}
    }

  // physics selection
  UInt_t trigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(fOfflineTrgMask==AliVEvent::kAny) fTriggerType = 0;
  else
    {
      if(trigger & fOfflineTrgMask) fTriggerType=0;

      if(fPeriod.Contains("lhc11a",TString::kIgnoreCase))
	{
	  if (trigger & AliVEvent::kMB)     { fTriggerType=0; }
	  if (trigger & AliVEvent::kEMC1)   { fTriggerType=1; }
	}  
    }
  if(fTriggerType==-1) return;
  if(fAnalyzeMCTruth) fTriggerType = 3;
  fhEventStat->Fill(1.5);
 
  // Vertex cut 
  const AliVVertex* vtx = fEvent->GetPrimaryVertex();
  if (!vtx || vtx->GetNContributors()<1) return;
  fhEventStat->Fill(2.5);
  fhVtxZ[fTriggerType]->Fill(vtx->GetZ());
  if (TMath::Abs(vtx->GetZ())>fMaxVtxZ)  return;
  fhEventStat->Fill(3.5);
  
  if(fTriggerType==0)
    {
      if(trigger & AliVEvent::kCentral) fhEventStat->Fill(5.5);
      else if (trigger & AliVEvent::kCentral) fhEventStat->Fill(6.5);
      else fhEventStat->Fill(4.5);
    }

  if(!AliAnalysisHelperJetTasks::Selected()) return;
  fhEventStat->Fill(7.5);
  
    // GetCentrality
  if(fCollisionSystem=="PbPb")
    {
      AliCentrality *centrality = fEvent->GetCentrality();
      if (centrality) fCentrality = centrality->GetCentralityPercentile("V0M");
      else            fCentrality = 99;
    }
  else if(fCollisionSystem=="pp")
    fCentrality = 0;
  fhCentrality[fTriggerType]->Fill(fCentrality);

  // Get Rho value
  if(fCollisionSystem=="PbPb")
    {
      if(fAnaType==0) fRhoValue = fRho->GetVal();
      if(fAnaType==1) fRhoValue = fEvtBkg->GetBackground(0);
    }
  else if(fCollisionSystem=="pp")
    {
      fRhoValue = 0;
    }
  fhRhoVsCent[fTriggerType]->Fill(fCentrality,fRhoValue);
  fhPtHardBins->Fill(fPtHardBin);

  if(fRunSingleInclHJet) 
    {
      Double_t trigPt = -1, trigPhi = -999, trigEta = -999;
      Int_t nTrig = FindSingleIncTrigger(fTrackArray, trigPt, trigPhi, trigEta, fIsEmbedding);
      if(nTrig>0) fhNumberOfTT[fTriggerType]->Fill(nTrig);
      RunSingleInclHJetCorr(trigPt, trigPhi, trigEta, fJetArray, fRhoValue, fhTTPt[fTriggerType][0], fHJetPhiCorr[fTriggerType][0]);
      if(fRunBkgFlow)
	{
	  RunSingleInclHJetCorr(trigPt, trigPhi, trigEta, fJetArray, fRhoValue+1.8, 0x0, fHJetPhiCorrUp[fTriggerType]);
	  RunSingleInclHJetCorr(trigPt, trigPhi, trigEta, fJetArray, fRhoValue-1.8, 0x0, fHJetPhiCorrDown[fTriggerType]);
	}

      if(fIsEmbedding)
	{
	  if(fRunDLHJet) 
	    {
	      FindSingleIncTrigger(fTrackArray, trigPt, trigPhi, trigEta, 1);
	      RunSingleInclHJetCorr(trigPt, trigPhi, trigEta, fDLJetArray, 0, fhTTPt[fTriggerType][1], fHJetPhiCorr[fTriggerType][1]);
	    }

	  if(fRunPLHJet) 
	    {
	      FindSingleIncTrigger(fMcParticleArray, trigPt, trigPhi, trigEta, 2);
	      RunSingleInclHJetCorr(trigPt, trigPhi, trigEta, fPLJetArray, 0, fhTTPt[fTriggerType][2], fHJetPhiCorr[fTriggerType][2]);
	    }
	}
    }
  
  if(fRunJetQA) 
    {
      RunJetQA(fJetArray, fRhoValue, fhJetPt[fTriggerType][0], fhJetArea[fTriggerType][0], fhJetQA[fTriggerType][0]);
      if(fIsEmbedding)
	{
	  if(fRunDLHJet) RunJetQA(fDLJetArray, 0, fhJetPt[fTriggerType][1], fhJetArea[fTriggerType][1], fhJetQA[fTriggerType][1]);
	  if(fRunPLHJet)  RunJetQA(fPLJetArray, 0, fhJetPt[fTriggerType][2], fhJetArea[fTriggerType][2], fhJetQA[fTriggerType][2]);
	}
    }
  
  if(!fIsEmbedding)
    {
      if(fRunTrkQA) RunTrackQA();
      if(fRunLeadTrkQA) RunLeadTrkQA();
    }

  if(fStudyKtEffects) StudyKtEffects();

  PostData(1, fOutputList);
  return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskHJetDphi::FindSingleIncTrigger(const TClonesArray *trackArray, Double_t &trigPt, Double_t &trigPhi, Double_t &trigEta, const Int_t arrayType)
{
  Int_t trigIndex = -1;
  trigPt = -1;
  trigPhi = -999; 
  trigEta = -999;

  // arrayType: 0 -> data, 1 -> embedding, 2 -> MC
  if(!trackArray) return 0;

  Int_t nTT = 0, counter = 0; 
  Int_t ntracks = trackArray->GetEntries();
  TArrayI arr;
  arr.Set(ntracks);
  for(Int_t it=0; it<ntracks; it++)
    {
      AliVParticle *t = dynamic_cast<AliVParticle*>(trackArray->At(it));
      if(!t || t->Charge()==0 || !AcceptTrack(t)) continue;

      if(fAnaType==0 && arrayType==1 && t->GetLabel()==0) continue; 

      if(fAnaType==1 && arrayType<2 && !IsGoodAODtrack(t) ) continue;

      Double_t pt = t->Pt();
      if(fTTtype==0)
	{
	  if (pt<fTTMaxPt && pt>=fTTMinPt)
	    {
	      nTT++;
	      arr.AddAt(it,counter);
	      counter++;
	    }
	}
    }
  arr.Set(counter);

  if(fTTtype==0)
    {
      if(counter==0) trigIndex = -1;
      else if(counter==1) trigIndex = arr.At(0);
      else
	{
	  fRandom->SetSeed(arr.At(0)); //make this random selection reproducible
	  Double_t pro = fRandom->Uniform() * counter;
	  trigIndex = arr.At(TMath::FloorNint(pro));
	}
      arr.Reset();
    }

  if(trigIndex>-1)
    {
      AliVParticle *tt = (AliVParticle*) fTrackArray->At(trigIndex);
      trigPt = tt->Pt();
      trigPhi = tt->Phi();
      trigEta = tt->Eta();

      if(fSwitchOnAvoidTpcHole)
	{
	  if(fAvoidTpcHole==1 && !(trigPhi>3.89 && trigPhi<5.53)) trigIndex = -1;
	  if(fAvoidTpcHole==2 && !(trigPhi>2.45 && trigPhi<3.44)) trigIndex = -1;
	}

      if(fCutTPCBoundary)
	{
	  Double_t phiDist = trigPhi - TMath::FloorNint(trigPhi/kSector)*kSector;
	  if(phiDist<fDistToTPCBoundary || phiDist>kSector-fDistToTPCBoundary)
	    {
	      trigIndex = -1;
	    }
	}
    }

  if(trigIndex==-1) { trigPt = -1; trigPhi = -999; trigEta = -999;}
  if(arrayType<2) fTriggerTrkIndex = trigIndex; 
  return nTT;
}

//________________________________________________________________________
void AliAnalysisTaskHJetDphi::RunSingleInclHJetCorr(Double_t trigPt, Double_t trigPhi, Double_t trigEta, const TClonesArray *jetArray, Double_t rho, THnSparse *hTT, THnSparse *hn)
{
  if(trigPt<0 || !fJetArray) return;

  if(hTT)
    {
      Double_t fillTT[] = {trigPt, fCentrality, (Double_t)fPtHardBin,static_cast<Double_t>(Entry()%10)};
      hTT->Fill(fillTT);
    }

  Int_t nJets = jetArray->GetEntries();
  Double_t jetPt = 0, jetEta = 0, jetPhi = 0, jetArea = 0;
  for(Int_t ij=0; ij<nJets; ij++)
    {
      if(fAnaType==0)
	{
	  AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(jetArray->At(ij));
	  jetPt   = jet->Pt();
	  jetEta  = jet->Eta();
	  jetPhi  = jet->Phi();
	  jetArea = jet->Area();
	}
      else if(fAnaType==1)
	{
	  AliAODJet* jet = dynamic_cast<AliAODJet*>(jetArray->At(ij));
	  jetPt   = jet->Pt();
	  jetEta  = jet->Eta();
	  jetPhi  = jet->Phi();
	  jetArea = jet->EffectiveAreaCharged();
	}
      else
	return;
      if(!IsGoodJet(jetEta)) continue; // eta cut
      Double_t dPhi = CalculateDPhi(trigPhi,jetPhi);

      Double_t jetPtCorr = jetPt-rho*jetArea;
      if(jetPtCorr>fJetPtMin)
	{
	  Double_t fill[] = {trigPt,jetPtCorr,dPhi,jetArea,fCentrality,trigEta-jetEta, (Double_t)fPtHardBin,static_cast<Double_t>(Entry()%10)};
	  hn->Fill(fill);
	}
    }
}


//________________________________________________________________________
void AliAnalysisTaskHJetDphi::RunTrackQA()
{
  if(fIsEmbedding) return;
  if(!fTrackArray) return;

  Int_t ntracks = fTrackArray->GetEntries();
  for(Int_t it=0; it<ntracks; it++)
    {
      AliVParticle *t = dynamic_cast<AliVParticle*>(fTrackArray->At(it));
      if(!t || t->Charge()==0 || !AcceptTrack(t)) continue;
      if(fAnaType==1 && !IsGoodAODtrack(t)) continue;
      fhTrkPt[fTriggerType]->Fill(t->Pt(),fCentrality);
      Double_t phi = t->Phi();
      if(phi<0) phi += 2*pi;
      Double_t fill[] = {t->Pt(), phi*TMath::RadToDeg(), t->Eta(), fCentrality};
      fhTrkQA[fTriggerType]->Fill(fill);
    }

  if(fESD)
    {
      Int_t ntrack = fESD->GetNumberOfTracks();
      for(Int_t itr=0; itr<ntrack; itr++)
	{
	  AliESDtrack *esdtrack = fESD->GetTrack(itr);
	  if(!esdtrack ||  TMath::Abs(esdtrack->Eta())>0.9 )continue;
	  Int_t type = -1;
	  if(fEsdTrkCut->AcceptTrack(esdtrack))
	    type = 0;
	  else if(fEsdHybCut->AcceptTrack(esdtrack))
	    type = 1;
	  else
	    continue;

	  Double_t resolution = -1;
	  Int_t label = esdtrack->GetLabel();
	  if(label>0 && fMC)
	    {
	      AliMCParticle *part = (AliMCParticle*)fMC->GetTrack(label);
	      TParticle *tpart = (TParticle*)part->Particle();
	      if(TMath::Abs(tpart->GetPdgCode())==211)
		{
		  Double_t fillPhiRes[] = {esdtrack->Pt(),part->Phi()-esdtrack->Phi(),(Double_t)type,fCentrality};
		  fhTrkPhiRes[fTriggerType]->Fill(fillPhiRes);
		}
	      resolution = (part->Pt()-esdtrack->Pt())/part->Pt();
	    }
	  Double_t fillPtRes[] = {esdtrack->Pt(),esdtrack->Pt()*TMath::Sqrt(esdtrack->GetSigma1Pt2()),resolution,(Double_t)type,fCentrality};
	  fhTrkPtRes[fTriggerType]->Fill(fillPtRes);
	}
    }
  else if(fAODIn)
    {
      ntracks = fAODIn->GetNumberOfTracks();
      Int_t type = -1;
      for(Int_t itrack=0; itrack<ntracks; itrack++)
	{
	  AliAODTrack *aodtrack = (AliAODTrack*)fAODIn->GetTrack(itrack);
	  if(!aodtrack || !AcceptTrack(dynamic_cast<AliVParticle*>(aodtrack)) ) continue;
	  if (fAODfilterBits[0] < 0) 
	    {
	      if (aodtrack->IsHybridGlobalConstrainedGlobal())
		type = 3;
	      else
		continue;
	    }
	  else 
	    {
	      if (aodtrack->TestFilterBit(fAODfilterBits[0]))  type = 0;
	      else if (aodtrack->TestFilterBit(fAODfilterBits[1]) && (aodtrack->GetStatus()&AliESDtrack::kITSrefit)!=0)	type = 1;
	      else continue;
	    }
	  Double_t sigma1Pt2 = GetAODTrackPtRes(aodtrack);
	  Double_t resolution = -1;
	  Int_t label = aodtrack->GetLabel();
	  if(label>0 && fMC)
	    {
	      AliMCParticle *part = (AliMCParticle*)fMC->GetTrack(label);
	      resolution = (part->Pt()-aodtrack->Pt())/part->Pt();
	      Double_t fillPhiRes[] = {aodtrack->Pt(),part->Phi()-aodtrack->Phi(),(Double_t)type,fCentrality};
	      fhTrkPhiRes[fTriggerType]->Fill(fillPhiRes);
	    }
	  if(sigma1Pt2>0)
	    {
	      Double_t fillPtRes[5] = {aodtrack->Pt(),aodtrack->Pt()*TMath::Sqrt(sigma1Pt2),resolution,(Double_t)type,fCentrality};
	      fhTrkPtRes[fTriggerType]->Fill(fillPtRes);
	    }
	}
    }
}

//________________________________________________________________________
void AliAnalysisTaskHJetDphi::RunJetQA(const TClonesArray *jetArray, const Double_t rho, THnSparse *hJetPt, THnSparse *hJetArea, THnSparse *hJetQA)
{
  Int_t nJets = jetArray->GetEntries();
  Double_t jetPt = 0, jetEta = 0, jetPhi = 0, jetArea = 0, jetPtCorr = 0;
  Int_t nCons = 0;
  for(Int_t ij=0; ij<nJets; ij++)
    {
      if(fAnaType==0)
	{
	  AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(jetArray->At(ij));
	  jetPt   = jet->Pt();
	  jetEta  = jet->Eta();
	  jetPhi  = jet->Phi();
	  jetArea = jet->Area();
	  nCons   = jet->GetNumberOfConstituents();
	}
      else if(fAnaType==1)
	{
	  AliAODJet* jet = dynamic_cast<AliAODJet*>(jetArray->At(ij));
	  jetPt   = jet->Pt();
	  jetEta  = jet->Eta();
	  jetPhi  = jet->Phi();
	  jetArea = jet->EffectiveAreaCharged();
	  nCons   = jet->GetRefTracks()->GetEntriesFast();
	}
      else
	return;
      if(!IsGoodJet(jetEta)) continue; // eta cut

      jetPtCorr = jetPt-rho*jetArea;
      Double_t fillPt[] = {jetPtCorr, jetPt, fCentrality, (Double_t)fPtHardBin};
      hJetPt->Fill(fillPt);

      Double_t fillA[] = {jetPtCorr, jetArea, fCentrality, (Double_t)fPtHardBin};
      hJetArea->Fill(fillA);

      Double_t fillQA[] = {jetPtCorr, jetPhi*TMath::RadToDeg(), jetEta, jetArea, (Double_t)nCons, fCentrality, (Double_t)fPtHardBin};
      hJetQA->Fill(fillQA);
    }
}


//________________________________________________________________________
void AliAnalysisTaskHJetDphi::RunLeadTrkQA()
{
  if(fIsEmbedding || fTriggerTrkIndex<0) return;
  Double_t jetPt = -1, jetPhi = -999;

  if(fAnaType==0)
    {
      Int_t nJets = fJetArray->GetEntries();
      for(Int_t ij=0; ij<nJets; ij++)
	{
	  AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fJetArray->At(ij));
	  if(!IsGoodJet(jet->Eta())) continue; // eta cut
	  for (Int_t i = 0; i < jet->GetNumberOfTracks(); i++) 
	    {
	      if(jet->TrackAt(i)==fTriggerTrkIndex)
		{
		  jetPt = jet->Pt();
		  jetPhi = jet->Phi();
		  break;
		}
	    }
	}
    }
  if(fAnaType==1)
    {
      jetPt = -1;
    }
  

  if(jetPt<=0) return;
  AliVParticle *tt = (AliVParticle*) fTrackArray->At(fTriggerTrkIndex);
  Double_t fill[] = {tt->Pt(), jetPt, tt->Phi()-jetPhi, tt->Pt()/jetPt, fCentrality};
  fhLeadTrkQA[fTriggerType]->Fill(fill);
}


//________________________________________________________________________
void AliAnalysisTaskHJetDphi::StudyKtEffects()
{
  if(fAnaType==1) return;
  if(fTriggerTrkIndex<0) return;

  fKtValue = 0; // dummy


  AliPicoTrack *tt = (AliPicoTrack*) fTrackArray->At(fTriggerTrkIndex);
  Double_t triggerPhi = tt->GetTrackPhiOnEMCal();

  Int_t nJets = fDLJetArray->GetEntries();
  for(Int_t ij=0; ij<nJets; ij++)
    {
      AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fDLJetArray->At(ij));
      if(!IsGoodJet(jet->Eta())) continue; // eta cut
      Double_t jetPhi = jet->Phi();
      Double_t jetPt  = jet->Pt();
      Double_t dPhi = CalculateDPhi(triggerPhi,jetPhi);
      if(dPhi<pi+0.6 && dPhi>pi-0.6)
	{
	  Double_t fill[] = {tt->Pt(), jetPt, jetPt*TMath::Tan(tt->GetTrackPtOnEMCal()), fCentrality, (Double_t)fPtHardBin};
	  fhKtEffects[fTriggerType]->Fill(fill);
	}
    }
}


//________________________________________________________________________
Bool_t AliAnalysisTaskHJetDphi::RetrieveArraies()
{
  if(fAnaType==0)
    {
      // Get mc particles
      if (!fMcParticleArrName.IsNull()) 
	{
	  AliInfo(Form("Retrieve mc particles %s!", fMcParticleArrName.Data()));
	  fMcParticleArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fMcParticleArrName));
	  if (!fMcParticleArray) 
	    {
	      AliError(Form("Could not retrieve mc particles %s!", fMcParticleArrName.Data())); 
	      return kFALSE;
	    }

	  TString fMcParticleMapArrName = fMcParticleArrName + "_Map";
	  fMcParticlelMap = dynamic_cast<AliNamedArrayI*>(fEvent->FindListObject(fMcParticleMapArrName));
	  if (!fMcParticlelMap) 
	    {
	      AliWarning(Form("%s: Could not retrieve map for MC particles %s! Will assume MC labels consistent with indexes...", GetName(), fMcParticleArrName.Data())); 
	      return kFALSE;
	    }
	}

      // Get track collection
      if (!fTrackArrName.IsNull()) 
	{
	  AliInfo(Form("Retrieve tracks %s!", fTrackArrName.Data()));
	  fTrackArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fTrackArrName));
	  if (!fTrackArray) 
	    {
	      AliError(Form("Could not retrieve tracks %s!", fTrackArrName.Data())); 
	      return kFALSE;
	    }
	}

      // Get mixed track collection
      if (!fEmbTrkArrName.IsNull()) 
	{
	  AliInfo(Form("Retrieve PYTHIA+PbPb tracks %s!", fEmbTrkArrName.Data()));
	  fEmbTrkArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fEmbTrkArrName));
	  if (!fEmbTrkArray) 
	    {
	      AliError(Form("Could not retrieve PYTHIA+PbPb tracks %s!", fEmbTrkArrName.Data())); 
	      return kFALSE;
	    }
	}
  
      // Get Rho value
      if(fCollisionSystem=="PbPb")
	{
	  if(!fRhoName.IsNull())
	    {
	      AliInfo(Form("Retrieve rho %s!", fRhoName.Data()));
	      fRho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoName));
	      if(!fRho)
		{
		  AliError(Form("Could not retrieve rho %s!",fRhoName.Data()));
		  return kFALSE;
		}
	    }
	}
      
      // Get jet collection
      if (!fJetArrName.IsNull())
	{
	  AliInfo(Form("Retrieve jets %s!", fJetArrName.Data()));
	  fJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fJetArrName));
	  if (!fJetArray)
	    {
	      AliError(Form("%s: Could not retrieve jets %s!", GetName(), fJetArrName.Data()));
	      return kFALSE;
	    }
	}

      // Get DL jet collection
      if (!fDLJetArrName.IsNull())
	{
	  AliInfo(Form("Retrieve DL jets %s!", fDLJetArrName.Data()));
	  fDLJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fDLJetArrName));
	  if (!fDLJetArray)
	    {
	      AliError(Form("%s: Could not retrieve DL jets %s!", GetName(), fDLJetArrName.Data()));
	      return kFALSE;
	    }
	}

      // Get PL jet collection
      if (!fPLJetArrName.IsNull())
	{
	  AliInfo(Form("Retrieve PL jets %s!", fPLJetArrName.Data()));
	  fPLJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fPLJetArrName));
	  if (!fPLJetArray)
	    {
	      AliError(Form("%s: Could not retrieve PL jets %s!", GetName(), fPLJetArrName.Data()));
	      return kFALSE;
	    }
	}

      if(fIsEmbedding)
	{
	  if(fAnaType==0 && !fPtHardBinName)
	    {
	      // Get embedded pt hard bin number
	      fPtHardBinName = static_cast<AliNamedString*>(fEvent->FindListObject("AODEmbeddingFile"));
	      if(!fPtHardBinName)
		{
		  AliError("The object for pt hard bin information is not available!");
		  return kFALSE;
		}
	    }
	}
    }

  if(fAnaType==1)
    {
      // Get mc particles
      if (!fMcParticleArrName.IsNull()) 
	{
	  AliInfo(Form("Retrieve mc particles %s!", fMcParticleArrName.Data()));
	  if(fAODOut        && !fMcParticleArray) fMcParticleArray = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fMcParticleArrName));
	  if(fAODIn         && !fMcParticleArray) fMcParticleArray = dynamic_cast<TClonesArray*>(fAODIn ->FindListObject(fMcParticleArrName));
	  if(fAODExtension  && !fMcParticleArray) fMcParticleArray = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fMcParticleArrName));
	  if (!fMcParticleArray) 
	    {
	      AliError(Form("Could not retrieve mc particles %s!", fMcParticleArrName.Data())); 
	      return kFALSE;
	    }
	}

      // Get tracks
      if (!fTrackArrName.IsNull()) 
	{
	  AliInfo(Form("Retrieve tracks %s!", fTrackArrName.Data()));
	  if(fAODIn         && !fTrackArray) fTrackArray = dynamic_cast<TClonesArray*>(fAODIn ->FindListObject(fTrackArrName));
	  if(fAODOut        && !fTrackArray) fTrackArray = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fTrackArrName));
	  if (!fTrackArray) 
	    {
	      AliError(Form("Could not retrieve tracks %s!", fTrackArrName.Data())); 
	      return kFALSE;
	    }
	}

      // Get PYTHIA+PbPb tracks
      if (!fEmbTrkArrName.IsNull()) 
	{
	  AliInfo(Form("Retrieve PYTHIA+PbPb tracks %s!", fEmbTrkArrName.Data()));
	  if(fAODOut        && !fEmbTrkArray) fEmbTrkArray = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fEmbTrkArrName));
	  if(fAODIn         && !fEmbTrkArray) fEmbTrkArray = dynamic_cast<TClonesArray*>(fAODIn ->FindListObject(fEmbTrkArrName));
	  if(fAODExtension  && !fEmbTrkArray) fEmbTrkArray = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fEmbTrkArrName));
	  if (!fTrackArray) 
	    {
	      AliError(Form("Could not retrieve PYTHIA+PbPb tracks %s!", fTrackArrName.Data())); 
	      return kFALSE;
	    }
	}

      // Get jet collection
      if (!fJetArrName.IsNull())
	{
	  AliInfo(Form("Retrieve jets %s!", fJetArrName.Data()));
	  if(fAODOut        && !fJetArray) fJetArray = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fJetArrName));
	  if(fAODExtension  && !fJetArray) fJetArray = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetArrName));
	  if(fAODIn         && !fJetArray) fJetArray = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fJetArrName));
	  if (!fJetArray)
	    {
	      AliError(Form("%s: Could not retrieve jets %s!", GetName(), fJetArrName.Data()));
	      return kFALSE;
	    }
	}

      // Get DL jet collection
      if (!fDLJetArrName.IsNull())
	{
	  AliInfo(Form("Retrieve DL jets %s!", fDLJetArrName.Data()));
	  if(fAODOut        && !fDLJetArray) fDLJetArray = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fDLJetArrName));
	  if(fAODExtension  && !fDLJetArray) fDLJetArray = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fDLJetArrName));
	  if(fAODIn         && !fDLJetArray) fDLJetArray = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fDLJetArrName));
	  if (!fDLJetArray)
	    {
	      AliError(Form("%s: Could not retrieve DL jets %s!", GetName(), fDLJetArrName.Data()));
	      return kFALSE;
	    }
	}

      // Get PL jet collection
      if (!fPLJetArrName.IsNull())
	{
	  AliInfo(Form("Retrieve PL jets %s!", fPLJetArrName.Data()));
	  if(fAODOut        && !fPLJetArray) fPLJetArray = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fPLJetArrName));
	  if(fAODExtension  && !fPLJetArray) fPLJetArray = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fPLJetArrName));
	  if(fAODIn         && !fPLJetArray) fPLJetArray = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fPLJetArrName));
	  if (!fPLJetArray)
	    {
	      AliError(Form("%s: Could not retrieve PL jets %s!", GetName(), fPLJetArrName.Data()));
	      return kFALSE;
	    }
	}

      // Get Rho 
      if(fCollisionSystem=="PbPb" && !fRhoName.IsNull())
	{
	  AliInfo(Form("Retrieve rho %s!", fRhoName.Data()));
	  if(fAODOut       && !fEvtBkg ) fEvtBkg = (AliAODJetEventBackground*)(fAODOut->FindListObject(fRhoName));
	  if(fAODExtension && !fEvtBkg ) fEvtBkg = (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fRhoName));
	  if(fAODIn        && !fEvtBkg ) fEvtBkg = (AliAODJetEventBackground*)(fAODIn->FindListObject(fRhoName));
	  if(!fEvtBkg)
	    {
	      AliError(Form("Could not retrieve rho %s!",fRhoName.Data()));
	      return kFALSE;
	    }
	}
    }

  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetDphi::CalculatePhi(const Double_t py, const Double_t px)
{
  Double_t phi = TMath::ATan2(py,px);
  if(phi<0) phi += 2*pi;
  return phi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetDphi::CalculateDPhi(const Double_t phi1, const Double_t phi2)
{
  Double_t dPhi = phi1-phi2;
  if(dPhi>2*pi)  dPhi -= 2*pi;
  if(dPhi<-2*pi) dPhi += 2*pi;
  if(dPhi<-0.5*pi) dPhi += 2*pi;
  if(dPhi>1.5*pi)  dPhi -= 2*pi;
  return dPhi;
}

//________________________________________________________________________
Int_t AliAnalysisTaskHJetDphi::LocateToTPCHole(const Double_t phi)
{
  if(phi<4-pi/2 && phi>5.5-1.5*pi) return 1; // away-side
  else return 0;   //near-side
}


//________________________________________________________________________
Bool_t AliAnalysisTaskHJetDphi::AcceptTrack(AliVParticle *track)
{
  if(track->Pt()<fMinTrkPt || track->Pt()>fMaxTrkPt) return kFALSE;
  if(track->Eta()<fMinTrkEta || track->Eta()>fMaxTrkEta) return kFALSE;
  if(track->Phi()<fMinTrkPhi || track->Phi()>fMaxTrkPhi) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetDphi::IsGoodAODtrack(AliVParticle *track)
{
  AliAODTrack *aodtrack = static_cast<AliAODTrack*>(track);
  if( fFilterMask>0)
    {
      if(!aodtrack->TestFilterBit(fFilterMask) ) return kFALSE;
    }
  else
    {
      if(!aodtrack->IsHybridGlobalConstrainedGlobal()) return kFALSE;
    }
  if( fRequireITSRefit && (aodtrack->GetStatus()&AliESDtrack::kITSrefit)==0 ) return kFALSE;
  if (fRequireSharedClsCut)
    {
      Double_t frac = Double_t(aodtrack->GetTPCnclsS())/Double_t(aodtrack->GetTPCncls());
      if (frac > 0.4) return kFALSE;
    }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetDphi::IsGoodJet(Double_t jetEta)
{
  Double_t etaCut = (0.9-fRadius>0.5)?0.5:0.9-fRadius;
  if(TMath::Abs(jetEta)>etaCut) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Int_t AliAnalysisTaskHJetDphi::GetPtHardBin(Double_t ptHard)
{
  const Int_t nBins = 10;
  Double_t binLimits[nBins] = { 5., 11., 21., 36., 57., 84., 117., 156., 200., 249. }; // lower limits
  Int_t bin = -1;
  while(bin<nBins-1 && binLimits[bin+1]<ptHard)
    {
      bin++;
    }
  return bin;
}

//________________________________________________________________________
Int_t AliAnalysisTaskHJetDphi::GetParticleType(Int_t pdg_input)
{
  Int_t type = -1;
  Int_t pdg = TMath::Abs(pdg_input);
  if(pdg==211) type = 0;
  else if(pdg==321) type = 1;
  else if(pdg==2212) type = 2;
  else if(pdg==11) type = 3;
  else if(pdg>=3122 && pdg<=3334) type = 4;
  else type = 9;
  return type;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetDphi::GetAODTrackPtRes(AliAODTrack *aodtrack)
{
  Double_t sigma1Pt2 = -1;
  Double_t cov[21] = {0,}, pxpypz[3] = {0,}, xyz[3] = {0,};
  AliExternalTrackParam *exParam = new  AliExternalTrackParam();
  aodtrack->GetCovMatrix(cov);
  aodtrack->PxPyPz(pxpypz);
  aodtrack->GetXYZ(xyz);
  exParam->Set(xyz,pxpypz,cov,aodtrack->Charge());
  sigma1Pt2 = exParam->GetSigma1Pt2();
  delete exParam;
  return sigma1Pt2;
}

//
//________________________________________________________________________
//
void AliAnalysisTaskHJetDphi::PrintConfig()
{
  const char *decision[2] = {"no","yes"};
  printf("\n\n===== h-jet analysis configuration =====\n");
  printf("Input event type: %s - %s\n",fCollisionSystem.Data(),fPeriod.Data());
  printf("Is this MC data: %s\n",decision[fIsMC]);
  printf("Run on particle level: %s\n",decision[fAnalyzeMCTruth]);
  printf("Event type selection: %d\n",fOfflineTrgMask);
  printf("Is embedding? %s\n",decision[fIsEmbedding]);
  printf("Track filter mask: %d\n",fFilterMask);
  printf("Require track to have ITS refit? %s\n",decision[fRequireITSRefit]);
  printf("Require to cut on fraction of shared TPC clusters? %s\n",decision[fRequireSharedClsCut]);
  printf("Track pt range: %2.2f < pt < %2.2f\n",fMinTrkPt, fMaxTrkPt);
  printf("Track eta range: %2.1f < eta < %2.1f\n",fMinTrkEta, fMaxTrkEta);
  printf("Track phi range: %2.0f < phi < %2.0f\n",fMinTrkPhi*TMath::RadToDeg(),fMaxTrkPhi*TMath::RadToDeg());
  printf("Cut TT away from boundary: %s with distance = %2.2f\n",decision[fCutTPCBoundary],fDistToTPCBoundary);
  printf("Avoid TPC holes: %s\n", decision[fSwitchOnAvoidTpcHole]);
  printf("Jet cone size R = %1.1f, and jet pt > %1.0f GeV/c \n",fRadius,fJetPtMin);
  printf("Run track QA: %s\n",decision[fRunTrkQA]);
  printf("Run jet QA: %s\n",decision[fRunJetQA]);
  printf("Run single inclusive h+jet analysis: %s\n",decision[fRunSingleInclHJet]);
  printf("TT interval:    %2.0f < pt < %2.0f\n",fTTMinPt, fTTMaxPt);
  printf("Run leading track QA: %s\n",decision[fRunLeadTrkQA]);
  printf("Run kT effects study: %s\n",decision[fStudyKtEffects]);
  printf("Run background flow: %s\n",decision[fRunBkgFlow]);
  printf("=======================================\n\n");
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetDphi::GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz)
{
  return (trkPx*jetPx+trkPy*jetPy+trkPz*jetPz)/(jetPx*jetPx+jetPy*jetPy+jetPz*jetPz);
}

//________________________________________________________________________
void AliAnalysisTaskHJetDphi::Terminate(Option_t *) 
{
  // Called once at the end of the query
}
