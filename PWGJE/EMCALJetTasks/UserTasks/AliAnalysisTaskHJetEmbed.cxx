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
#include <TParameter.h>

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskHJetEmbed.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliNamedString.h"

#include "AliEmcalJet.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

ClassImp(AliAnalysisTaskHJetEmbed)

const Double_t pi = TMath::Pi();
//const Double_t areaCut[4] = {0.1, 0.23, 0.4, 0.63};

//________________________________________________________________________
AliAnalysisTaskHJetEmbed::AliAnalysisTaskHJetEmbed() : 
  AliAnalysisTaskSE(), 
  fVerbosity(0), fAnaType(1), fPeriod("lhc11h"), fCollisionSystem("PbPb"),
  fEvent(0), fTriggerType(-1), fCentrality(-1), fMaxVtxZ(10),
  fMCParticleArrName(""), fMCParticleArray(0x0), 
  fTrackArrName(""), fTrackArray(0x0), fTriggerTrkIndex(-1), 
  fMinTrkPt(0.15), fMaxTrkPt(1e4), fMinTrkEta(-0.9), fMaxTrkEta(0.9), fMinTrkPhi(0), fMaxTrkPhi(2*pi), 
  fTTtype(0), 
  fRadius(0.4), fJetArrName(""), fPLJetArrName(""), fDLJetArrName(""),
  fJetArray(0x0), fPLJetArray(0x0), fDLJetArray(0x0),
  fRhoName(""), fRho(0x0), fRhoValue(0), 
  fPtHardBinName(0x0), fPtHardBin(-1), fRandom(0x0),
  fRunQA(kTRUE), fRunHJet(kTRUE), fRunMatch(kTRUE),
  fRunPL(kFALSE), fRunDL(kTRUE),
  fOutputList(0x0), fhEventStat(0x0), fhPtHardBins(0x0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  //DefineOutput(1, TList::Class());
  // Output slot #0 id reserved by the base class for AOD

  for(Int_t i=0; i<kNTrig; i++)
    {
      fhVtxZ[i]                  = 0x0;
      fhCentrality[i]            = 0x0;
      fhRhoVsCent[i]             = 0x0;

      fhDLJetPtVsCent[i]         = 0x0;
      fhPLJetPtVsCent[i]         = 0x0;

      fhPLTT[i]                  = 0x0;
      fhDLTT[i]                   = 0x0;
      fhPLHJet[i]                = 0x0;
      fhDLHJet[i]                = 0x0;
      fhTTPtQA[i]                = 0x0;
      fhTTPt[i]                  = 0x0;
      fhHJet[i]                  = 0x0;

      fhJetPtGeoMatch[i]         = 0x0;
      fhJetPtEnMatch[i]          = 0x0;
      fhJetPhiGeoMatch[i]        = 0x0;
      fhJetPhiEnMatch[i]         = 0x0;
    }
  for(Int_t i=0; i<kNTT; i++)
    {
      if(i==0)
	{ fMinTTPt[i] = 19; fMaxTTPt[i] = 25; }
      else
	{ fMinTTPt[i] = -1; fMaxTTPt[i] = -1; }
    }
}

//________________________________________________________________________
AliAnalysisTaskHJetEmbed::AliAnalysisTaskHJetEmbed(const char *name) : 
  AliAnalysisTaskSE(name), 
  fVerbosity(0), fAnaType(1), fPeriod("lhc11h"), fCollisionSystem("PbPb"),
  fEvent(0), fTriggerType(-1), fCentrality(-1), fMaxVtxZ(10),
  fMCParticleArrName(""), fMCParticleArray(0x0), 
  fTrackArrName(""), fTrackArray(0x0), fTriggerTrkIndex(-1), 
  fMinTrkPt(0.15), fMaxTrkPt(1e4), fMinTrkEta(-0.9), fMaxTrkEta(0.9), fMinTrkPhi(0), fMaxTrkPhi(2*pi), 
  fTTtype(0), 
  fRadius(0.4), fJetArrName(""), fPLJetArrName(""), fDLJetArrName(""),
  fJetArray(0x0), fPLJetArray(0x0), fDLJetArray(0x0),
  fRhoName(""), fRho(0x0), fRhoValue(0), 
  fPtHardBinName(0x0), fPtHardBin(-1), fRandom(0x0),
  fRunQA(kTRUE), fRunHJet(kTRUE), fRunMatch(kTRUE),
  fRunPL(kFALSE), fRunDL(kTRUE),
  fOutputList(0x0), fhEventStat(0x0), fhPtHardBins(0x0)
{
  // Constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

  for(Int_t i=0; i<kNTrig; i++)
    {
      fhVtxZ[i]                  = 0x0;
      fhCentrality[i]            = 0x0;
      fhRhoVsCent[i]             = 0x0;

      fhDLJetPtVsCent[i]         = 0x0;
      fhPLJetPtVsCent[i]         = 0x0;

      fhPLTT[i]                  = 0x0;
      fhDLTT[i]                   = 0x0;
      fhPLHJet[i]                = 0x0;
      fhDLHJet[i]                = 0x0;
      fhTTPtQA[i]                = 0x0;
      fhTTPt[i]                  = 0x0;
      fhHJet[i]                  = 0x0;

      fhJetPtGeoMatch[i]         = 0x0;
      fhJetPtEnMatch[i]          = 0x0;
      fhJetPhiGeoMatch[i]        = 0x0;
      fhJetPhiEnMatch[i]         = 0x0;
    }
  for(Int_t i=0; i<kNTT; i++)
    {
      if(i==0)
	{ fMinTTPt[i] = 19; fMaxTTPt[i] = 25; }
      else
	{ fMinTTPt[i] = -1; fMaxTTPt[i] = -1; }
    }
}
//________________________________________________________________________
AliAnalysisTaskHJetEmbed::~AliAnalysisTaskHJetEmbed()
{
  //Destructor
  if(fRho)         delete fRho;
  if(fOutputList) { fOutputList->Delete(); delete fOutputList;}
}

//________________________________________________________________________
void AliAnalysisTaskHJetEmbed::UserCreateOutputObjects()
{
  // Create histograms

  const Int_t nJetPtBins = 40;
  const Float_t lowJetPtBin=-50, upJetPtBin=150;

  const Int_t nTrkPtBins = 100;
  const Float_t lowTrkPtBin=0, upTrkPtBin=100;


  // QA
  const Int_t dimJetqa = 3;
  const Int_t nBinsJetqa[dimJetqa]     = {nJetPtBins,  30, 11};
  const Double_t lowBinJetqa[dimJetqa] = {lowJetPtBin, 0,  0};
  const Double_t hiBinJetqa[dimJetqa]  = {upJetPtBin,  30, 11};

  // h+jet
  const Int_t dimTT = 3;
  const Int_t nBinsTT[dimTT]     = {nTrkPtBins,  30, 11};
  const Double_t lowBinTT[dimTT] = {lowTrkPtBin, 0,  0};
  const Double_t hiBinTT[dimTT]  = {upTrkPtBin,  30, 11};  

  const Int_t dimHJet = 6;
  const Int_t nBinsHJet[dimHJet]     = {nTrkPtBins,  nJetPtBins,  140,     8,   30, 11};
  const Double_t lowBinHJet[dimHJet] = {lowTrkPtBin, lowJetPtBin, pi-4.95, 0,   0,  0};
  const Double_t hiBinHJet[dimHJet]  = {upTrkPtBin,  upJetPtBin,  pi+2.05, 0.8, 30, 11};  

  // Match
  const Int_t dimMthPt = 6;
  const Int_t nBinsMthPt[dimMthPt]     = {nJetPtBins,  nJetPtBins,  20,    20, 30, 11};
  const Double_t lowBinMthPt[dimMthPt] = {lowJetPtBin, lowJetPtBin, -0.95, 0,  0,  0};
  const Double_t hiBinMthPt[dimMthPt]  = {upJetPtBin,  upJetPtBin,  1.05,  1,  30, 11};

  const Int_t dimMthPhi = 7;
  const Int_t nBinsMthPhi[dimMthPhi]     = {nTrkPtBins,  nJetPtBins/2,  70,      70,      10, 1,  11};
  const Double_t lowBinMthPhi[dimMthPhi] = {lowTrkPtBin, lowJetPtBin,   pi-4.95, pi-4.95, 0,  0,  0};
  const Double_t hiBinMthPhi[dimMthPhi]  = {upTrkPtBin,  upJetPtBin,    pi+2.05, pi+2.05, 1,  10, 11};

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fhEventStat = new TH1F("fhEventStat","Event statistics for jet analysis",9,0,9);
  fhEventStat->GetXaxis()->SetBinLabel(1,"All");
  fhEventStat->GetXaxis()->SetBinLabel(2,"PS");
  fhEventStat->GetXaxis()->SetBinLabel(3,"Vtx");
  fhEventStat->GetXaxis()->SetBinLabel(4,"Vtx+10cm");
  fhEventStat->GetXaxis()->SetBinLabel(5,"kMB");
  fhEventStat->GetXaxis()->SetBinLabel(6,"kCentral");
  fhEventStat->GetXaxis()->SetBinLabel(7,"kSemiCentral");
  fhEventStat->GetXaxis()->SetBinLabel(8,"kEMCEGA");
  fhEventStat->GetXaxis()->SetBinLabel(9,"kEMCEJE");
  fOutputList->Add(fhEventStat);

  fhPtHardBins = new TH1F("fhPtHardBins","Number of events in each pT hard bin",11,0,11);
  fOutputList->Add(fhPtHardBins);

  const char *triggerName[kNTrig] = {"kMB","kEGA","kEJE"};
  
  for(Int_t i=0; i<kNTrig; i++)
    {
      if( i!=0 ) continue;
 
      fhVtxZ[i] = new TH1F(Form("%s_fhVtxZ",triggerName[i]),Form("%s: z distribution of event vertexz;z(cm)",triggerName[i]),400,-20,20);
      fOutputList->Add(fhVtxZ[i]);

      fhCentrality[i] = new TH1F(Form("%s_fhCentrality",triggerName[i]),Form("%s: Event centrality;centrality",triggerName[i]),100,0,100);
      fOutputList->Add(fhCentrality[i]);

      fhRhoVsCent[i] = new TH2F(Form("%s_fhRhoVsCent",triggerName[i]),Form("%s: Rho vs centrality (R=%1.1f);centrality;Rho",triggerName[i],fRadius),100,0,100,300,0,300);
      fOutputList->Add(fhRhoVsCent[i]);

      if(fRunQA)
	{
	  fhPLJetPtVsCent[i] = new THnSparseF(Form("%s_fhPLJetPtVsCent",triggerName[i]),Form("PYTHIA: jet p_{T} vs centrality vs pT hard bin (particle-level,R=%1.1f);p_{T,jet}^{ch} (GeV/c);centrality;pT hard bin",fRadius),dimJetqa,nBinsJetqa,lowBinJetqa,hiBinJetqa);
	  fOutputList->Add(fhPLJetPtVsCent[i]);

	  fhDLJetPtVsCent[i] = new THnSparseF(Form("%s_fhDLJetPtVsCent",triggerName[i]),Form("PYTHIA: jet p_{T} vs centrality vs pT hard bin (detector-level,R=%1.1f);p_{T,jet}^{ch} (GeV/c);centrality;pT hard bin",fRadius),dimJetqa,nBinsJetqa,lowBinJetqa,hiBinJetqa);
	  fOutputList->Add(fhDLJetPtVsCent[i]);
	}

      if(fRunHJet)
	{
	  fhPLTT[i] = new THnSparseF(Form("%s_fhPLTT",triggerName[i]),Form("PYTHIA: TT p_{T} vs centrality vs pT hard bin (particle-level);p_{T,TT}^{ch} (GeV/c);centrality;pT hard bin"),dimTT,nBinsTT,lowBinTT,hiBinTT);
	  fOutputList->Add(fhPLTT[i]);

	  fhDLTT[i] = new THnSparseF(Form("%s_fhDLTT",triggerName[i]),Form("PYTHIA: TT p_{T} vs centrality vs pT hard bin (detector-level);p_{T,TT}^{ch} (GeV/c);centrality;pT hard bin"),dimTT,nBinsTT,lowBinTT,hiBinTT);
	  fOutputList->Add(fhDLTT[i]);

	  fhTTPt[i] = new THnSparseF(Form("%s_fhTTPt",triggerName[i]),Form("Embedded: TT p_{T} vs centrality vs pT hard bin;p_{T,TT}^{ch} (GeV/c);centrality;pT hard bin"),dimTT,nBinsTT,lowBinTT,hiBinTT);
	  fOutputList->Add(fhTTPt[i]);

	  fhPLHJet[i] = new THnSparseF(Form("%s_fhPLHJet",triggerName[i]),Form("PYTHIA: TT p_{T} vs jet p_{T} vs #Delta#varphi vs jet area vs centrality vs pT hard bin (particle-level, R=%1.1f);p_{T,TT}^{ch} (GeV/c);p_{T,jet}^{ch} (GeV/c);#Delta#varphi;Area;centrality;pT hard bin",fRadius),dimHJet,nBinsHJet,lowBinHJet,hiBinHJet);
	  fOutputList->Add(fhPLHJet[i]);

	  fhDLHJet[i] = new THnSparseF(Form("%s_fhDLHJet",triggerName[i]),Form("PYTHIA: TT p_{T} vs jet p_{T} vs #Delta#varphi vs jet area vs centrality vs pT hard bin (detector-level, R=%1.1f);p_{T,TT}^{ch} (GeV/c);p_{T,jet}^{ch} (GeV/c);#Delta#varphi;Area;centrality;pT hard bin",fRadius),dimHJet,nBinsHJet,lowBinHJet,hiBinHJet);
	  fOutputList->Add(fhDLHJet[i]);

	  fhHJet[i] = new THnSparseF(Form("%s_fhHJet",triggerName[i]),Form("Embedded: TT p_{T} vs jet p_{T} vs #Delta#varphi vs jet area vs centrality vs pT hard bin (R=%1.1f);p_{T,TT}^{ch} (GeV/c);p_{T,jet}^{ch} (GeV/c);#Delta#varphi;Area;centrality;pT hard bin",fRadius),dimHJet,nBinsHJet,lowBinHJet,hiBinHJet);
	  fOutputList->Add(fhHJet[i]);
	}

      if(fRunMatch)
	{
	  fhJetPtGeoMatch[i] = new THnSparseF(Form("%s_fhJetPtGeoMatch",triggerName[i]),Form("Embed: generated p_{T,jet} vs reconstructed p_{T,jet} vs jet p_{T} difference vs dR vs centrality vs pT hard bin (R=%1.1f);p_{T,jet}^{gen} (GeV/c);p_{T,jet}^{rec} (GeV/c);(p_{T,jet}^{rec}-p_{T,jet}^{gen})/p_{T,jet}^{gen};dR;centrality;pT hard bin",fRadius),dimMthPt,nBinsMthPt,lowBinMthPt,hiBinMthPt);
	  fOutputList->Add(fhJetPtGeoMatch[i]);

	  fhJetPtEnMatch[i] = new THnSparseF(Form("%s_fhJetPtEnMatch",triggerName[i]),Form("Embed: generated p_{T,jet} vs reconstructed p_{T,jet} vs jet p_{T} difference vs dR vs centrality vs pT hard bin (R=%1.1f);p_{T,jet}^{gen} (GeV/c);p_{T,jet}^{rec} (GeV/c);(p_{T,jet}^{rec}-p_{T,jet}^{gen})/p_{T,jet}^{gen};dR;centrality;pT hard bin",fRadius),dimMthPt,nBinsMthPt,lowBinMthPt,hiBinMthPt);
	  fOutputList->Add(fhJetPtEnMatch[i]);

	  fhJetPhiGeoMatch[i] = new THnSparseF(Form("%s_fhJetPhiGeoMatch",triggerName[i]),Form("Embed: p_{T,TT} vs p_{T,jet}^{det} vs #Delta#varphi_{TT} vs #Delta#varphi_{jet} vs dR vs centrality vs pT hard bin (R=%1.1f);p_{T,TT} (GeV/c);p_{T,jet}^{det} (GeV/c);#Delta#varphi_{TT};#Delta#varphi_{jet};dR;centrality;pThard bin",fRadius),dimMthPhi,nBinsMthPhi,lowBinMthPhi,hiBinMthPhi);
	  fOutputList->Add(fhJetPhiGeoMatch[i]);

	  fhJetPhiEnMatch[i] = new THnSparseF(Form("%s_fhJetPhiEnMatch",triggerName[i]),Form("Embed: p_{T,TT} vs p_{T,jet}^{det} vs #Delta#varphi_{TT} vs #Delta#varphi_{jet} vs dR vs centrality vs pT hard bin (R=%1.1f);p_{T,TT} (GeV/c);p_{T,jet}^{det} (GeV/c);#Delta#varphi_{TT};#Delta#varphi_{jet};dR;centrality;pThard bin",fRadius),dimMthPhi,nBinsMthPhi,lowBinMthPhi,hiBinMthPhi);
	  fOutputList->Add(fhJetPhiEnMatch[i]);
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

  fRandom = new TRandom3();

  PrintConfig();
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHJetEmbed::UserExec(Option_t *) 
{  
  // Main loop, called for each event.

  AliDebug(5,"Entering UserExec");
  fTriggerType = -1;
  fEvent = InputEvent();
  if (!fEvent) 
    {
      AliError("Input event not available");
      return;
    }
  AliDebug(5,"Got the input event");
  if(!fPtHardBinName)
    {
      // Get embedded pt hard bin number
      fPtHardBinName = static_cast<AliNamedString*>(fEvent->FindListObject("AODEmbeddingFile"));
      if(!fPtHardBinName)
  	{
  	  AliError("The object for pt hard bin information is not available!");
  	  return;
  	}
    }
  TString fileName = fPtHardBinName->GetString();
  fileName.Remove(0,50);
  fileName.Remove(fileName.Index("/"));
  fPtHardBin = fileName.Atoi();
  fhEventStat->Fill(0.5);
  UInt_t trigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  if(fPeriod.Contains("lhc11h",TString::kIgnoreCase))
    {
      if (trigger & AliVEvent::kAnyINT)      { fTriggerType=0; }
      else if (trigger & AliVEvent::kCentral)     { fTriggerType=0; }
      else if (trigger & AliVEvent::kSemiCentral) { fTriggerType=0; }
      else if (trigger & AliVEvent::kEMCEGA)      { fTriggerType=1; }
      else if (trigger & AliVEvent::kEMCEJE)      { fTriggerType=2; }
    }
  else if(fPeriod.Contains("lhc10h",TString::kIgnoreCase))
    {
      if (trigger & AliVEvent::kAnyINT)   { fTriggerType=0; }
    }
  else if(fPeriod.Contains("lhc12a15a",TString::kIgnoreCase))
    {
      fTriggerType=0;
    }
  
  if(fTriggerType==-1) return;
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
  if(fTriggerType==1) fhEventStat->Fill(7.5);
  if(fTriggerType==2) fhEventStat->Fill(8.5);
  
  // GetCentrality
  if(fCollisionSystem=="PbPb")
    {
      AliCentrality *centrality = fEvent->GetCentrality();
      if (centrality)
      	fCentrality = centrality->GetCentralityPercentile("V0M");
      else 
	fCentrality = 99;
    }
  else if(fCollisionSystem=="pp")
    {
      fCentrality = 0;
    }

  // Get track collection and run QA
  if (!fTrackArray) 
    {
      fTrackArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fTrackArrName));
      if (!fTrackArray) 
	{
	  AliError(Form("Could not retrieve tracks %s!", fTrackArrName.Data())); 
	  return;
	}
      if (!fTrackArray->GetClass()->GetBaseClass("AliVParticle")) 
	{
	  AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTrackArrName.Data())); 
	  fTrackArray = 0;
	  return;
	}
    }

  // Get MC particle array
  if (fRunPL && !fMCParticleArray) 
    {
      fMCParticleArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fMCParticleArrName));
      if (!fMCParticleArray) 
	{
	  AliError(Form("Could not retrieve tracks %s!", fMCParticleArrName.Data())); 
	  return;
	}
    }

  // Get Rho value
  if(fCollisionSystem=="PbPb")
    {
      if(!fRho && !fRhoName.IsNull())
	{
	  fRho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoName));
	  if(!fRho)
	    {
	      AliError(Form("Could not retrieve rho %s!",fRhoName.Data()));
	      return;
	    }
	}
      if(fRho) fRhoValue = fRho->GetVal();
    }
  else if(fCollisionSystem=="pp")
    {
      fRhoValue = 0;
    }

  // Get jet collection
  if (!fJetArray && !fJetArrName.IsNull())
    {
      fJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fJetArrName));
      if (!fJetArray)
	{
	  AliError(Form("%s: Could not retrieve jets %s!", GetName(), fJetArrName.Data()));
	  return;
	}
      if (!fJetArray->GetClass()->GetBaseClass("AliEmcalJet")) 
	{
	  AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJetArrName.Data())); 
	  fJetArray = 0;
	  return;
	}
    }
  // Get particle-level jet array
  if (fRunPL && !fPLJetArray && !fPLJetArrName.IsNull())
    {
      fPLJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fPLJetArrName));
      if (!fPLJetArray)
	{
	  AliError(Form("%s: Could not retrieve jets %s!", GetName(), fPLJetArrName.Data()));
	  return;
	}
      if (!fPLJetArray->GetClass()->GetBaseClass("AliEmcalJet")) 
	{
	  AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fPLJetArrName.Data())); 
	  fPLJetArray = 0;
	  return;
	}
    }
 // Get detector-level jet array
  if (fRunDL && !fDLJetArray && !fDLJetArrName.IsNull())
    {
      fDLJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fDLJetArrName));
      if (!fDLJetArray)
	{
	  AliError(Form("%s: Could not retrieve jets %s!", GetName(), fDLJetArrName.Data()));
	  return;
	}
      if (!fDLJetArray->GetClass()->GetBaseClass("AliEmcalJet")) 
	{
	  AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fDLJetArrName.Data())); 
	  fDLJetArray = 0;
	  return;
	}
    }

  fhCentrality[fTriggerType]->Fill(fCentrality);
  fhRhoVsCent[fTriggerType]->Fill(fCentrality,fRhoValue);
  fhPtHardBins->Fill(fPtHardBin);

  if(fRunQA) RunQA();
  if(fRunHJet) 
    {
      for(Int_t i=0; i<kNTT; i++)
	RunHJet(fMinTTPt[i],fMaxTTPt[i]);
    }

  PostData(1, fOutputList);
  return;
}


//________________________________________________________________________
void AliAnalysisTaskHJetEmbed::RunHJet(const Double_t minPt, const Double_t maxPt)
{
  TArrayI arr;
  Int_t counter = 0;
  Int_t indexPL = -1;

  if(fRunPL)
    {
      // Find trigger track on particle level
      const Int_t nParticles = fMCParticleArray->GetEntries();
      Double_t maxPLPt = -1;
      arr.Set(nParticles);
      counter = 0;
      for(Int_t iPart=0; iPart<nParticles; iPart++)
	{
	  AliVParticle *t = static_cast<AliVParticle*>(fMCParticleArray->At(iPart));
	  //if(!t || t->Charge()==0) continue;
	  //if(!AcceptTrack(t)) continue;
	  Double_t pt = t->Pt();
	  if(fTTtype==0) // single inclusive triggers
	    {
	      if (pt<maxPt && pt>=minPt)
		{
		  arr.AddAt(iPart,counter);
		  counter++;
		}
	    }
	  else if(fTTtype==1) // leading triggers
	    {
	      if(maxPLPt<pt)
		{
		  maxPLPt = pt;
		  indexPL = iPart;
		}
	    }
	}
      arr.Set(counter);
      if(fTTtype==0)
	{
	  if(counter==0) indexPL = -1;
	  else if(counter==1) indexPL = arr.At(0);
	  else
	    {
	      Double_t pro = fRandom->Uniform() * counter;
	      indexPL = arr.At(TMath::FloorNint(pro));
	    }
	}
      arr.Reset();
    }


  // Find trigger track on detector level and after embedding
  const Int_t Ntracks = fTrackArray->GetEntries();
  Double_t maxDLPt = 0;
  Int_t indexDL = -1;
  arr.Set(Ntracks);
  counter = 0;
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) 
    {
      AliVParticle *t = static_cast<AliVParticle*>(fTrackArray->At(iTracks));
      if(!t || t->Charge()==0) continue;
      if(!AcceptTrack(t)) continue;
      if(t->GetLabel()!=0) 
	{
	  //cout<<iTracks<<"  "<<t->Pt()<<" "<<t->GetLabel()<<endl;
	  Double_t pt = t->Pt();
	  if(fTTtype==0) 
	    {
	      if (pt<maxPt && pt>=minPt)
		{
		  arr.AddAt(iTracks,counter);
		  counter++;
		}
	    }
	  else if(fTTtype==1)
	    {
	      if(maxDLPt<pt)
		{
		  maxDLPt = pt;
		  indexDL = iTracks;
		}
	    }
	}
    }
  arr.Set(counter);
  if(fTTtype==0)
    {
      if(counter==0) indexDL = -1;
      else if(counter==1) indexDL = arr.At(0);
      else
	{
	  Double_t pro = fRandom->Uniform() * counter;
	  indexDL = arr.At(TMath::FloorNint(pro));
	}
    }
  arr.Reset();

  AliDebug(2,Form("TT indices: PL=%d, DL=%d\n",indexPL,indexDL));

  // Run h+jet
  if(fRunPL)  FillHJetCor(fMCParticleArray, indexPL, fPLJetArray, fhPLTT[fTriggerType], fhPLHJet[fTriggerType], kFALSE);
  if(fRunDL)  FillHJetCor(fTrackArray,      indexDL, fDLJetArray, fhDLTT[fTriggerType], fhDLHJet[fTriggerType], kFALSE);
  FillHJetCor(fTrackArray,      indexDL, fJetArray,   fhTTPt[fTriggerType], fhHJet[fTriggerType],   kTRUE);  

  if(fRunMatch) RunMatch(fTrackArray, indexDL);
}

//________________________________________________________________________
void  AliAnalysisTaskHJetEmbed::FillHJetCor(const TClonesArray *tracks, const Int_t leadingIndex, const TClonesArray *jetArray, THnSparse *hTT, THnSparse *hn, Bool_t isBkg)
{
  if(leadingIndex<0) return;

  AliVParticle *tt = (AliVParticle*) tracks->At(leadingIndex);
  Double_t triggerPt = tt->Pt();
  Double_t fill1[] = {triggerPt, fCentrality, static_cast<Double_t>(fPtHardBin) };
  hTT->Fill(fill1);
  AliDebug(2,Form("Found a trigger with pt = %2.2f",triggerPt));

  Double_t triggerPhi = tt->Phi();
  if(triggerPhi<0) triggerPhi += 2*pi;
  Int_t nJets = jetArray->GetEntries();
  for(Int_t ij=0; ij<nJets; ij++)
    {
      AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(jetArray->At(ij));
      if(!jet) continue;
      if(!IsGoodJet(jet)) continue; // eta cut
      Double_t jetPhi = jet->Phi();
      Double_t jetPt  = jet->Pt();
      Double_t jetArea = jet->Area();
      Double_t dPhi = CalculateDPhi(triggerPhi,jetPhi);
      Double_t fill[] = {triggerPt,jetPt-jetArea*fRhoValue,dPhi,jetArea,fCentrality,static_cast<Double_t>(fPtHardBin)};
      if(!isBkg) fill[1] = jetPt; 
      AliDebug(10,"Fill the histograms");
      hn->Fill(fill);
    }
}

//________________________________________________________________________
void AliAnalysisTaskHJetEmbed::RunMatch(const TClonesArray *tracks, const Int_t leadingIndex)
{
  if(leadingIndex<0) return;

  if(!fDLJetArray || !fJetArray)
    {
      AliWarning("Jet array is not available.");
      return;
    }

  AliVParticle *tt = (AliVParticle*) tracks->At(leadingIndex);
  Double_t dR = 999, fraction = -1;
  Int_t nJets = fDLJetArray->GetEntries();
  for(Int_t ij=0; ij<nJets; ij++)
    {
      AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fDLJetArray->At(ij));
      if(!jet) continue;
      if(!IsGoodJet(jet)) continue; // eta cut
      Double_t jetPt = jet->Pt();
      if(jetPt<10) continue;
      
      // energy matching
      Int_t mthJetIndexEn = FindEnergyMatchedJet(jet,fJetArray,dR,fraction);
      if(mthJetIndexEn>-1 && fraction>0.5)
	{
	  AliEmcalJet* jetMthEn = dynamic_cast<AliEmcalJet*>(fJetArray->At(mthJetIndexEn));
	  if(jetMthEn)
	    {
	      Double_t fill[] = {tt->Pt(),jetPt,CalculateDPhi(tt->Phi(),jet->Phi()),CalculateDPhi(jetMthEn->Phi(),jet->Phi()),dR,fCentrality,static_cast<Double_t>(fPtHardBin)};
	      fhJetPhiEnMatch[fTriggerType]->Fill(fill);
	    }
	}
    }
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskHJetEmbed::FindGeoMatchedJet(const AliEmcalJet* jet, const TClonesArray *jetArray, Double_t &dR)
{
  dR = 999;
  if(!jetArray) return -1;
  
  Int_t index = -1;
  Int_t nJets = jetArray->GetEntries();
  Double_t dRMax = 1;
  for(Int_t ij=0; ij<nJets; ij++)
    {
      AliEmcalJet* jetTmp = dynamic_cast<AliEmcalJet*>(jetArray->At(ij)); 
      if(!jetTmp) continue;
      if(TMath::Abs(jetTmp->Eta())>1) continue; // Generous eta cut
      Double_t dPhi = GetDPhi(jet->Phi(),jetTmp->Phi());
      Double_t dEta = jet->Eta()-jetTmp->Eta();
      Double_t dRTmp = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
      if(dRTmp<dRMax)
	{
	  dRMax = dRTmp;
	  index = ij;
	}
    }
  dR = dRMax;
  return index;
}

//________________________________________________________________________
Int_t AliAnalysisTaskHJetEmbed::FindEnergyMatchedJet(const AliEmcalJet* jet, const TClonesArray *jetArray, Double_t &dR, Double_t &fraction)
{
  dR = 999;
  fraction=-1;
  if(!jetArray || !jet) return -1;
  
  Int_t index = -1;
  Int_t nJets = jetArray->GetEntries();
  Double_t maxFrac = 0;
  Int_t nJetC = (Int_t)jet->GetNumberOfConstituents();
  Double_t jetPt = jet->Pt();
  for(Int_t ij=0; ij<nJets; ij++)
    {
      AliEmcalJet* jetTmp = dynamic_cast<AliEmcalJet*>(jetArray->At(ij)); 
      if(!jetTmp) continue;
      if(TMath::Abs(jetTmp->Eta())>1) continue; // Generous eta cut
      if(GetJetDistance(jet,jetTmp)>1) continue;

      Int_t nc = (Int_t)jetTmp->GetNumberOfConstituents();
      Double_t sumPt = 0;
      for(Int_t ic=0; ic<nc; ic++)
	{
	  for(Int_t ijc=0; ijc<nJetC; ijc++)
	    {
	      if(jetTmp->TrackAt(ic)==jet->TrackAt(ijc))
		{
		  AliVParticle *part = (AliVParticle*)jet->TrackAt(ijc,fTrackArray);
		  sumPt += part->Pt();
		}
	    }
	}
      Double_t frac = sumPt/jetPt;
      if(frac>maxFrac)
	{
	  maxFrac = frac;
	  index = ij;
	}
    }
  fraction = maxFrac;

  if(index>0)
    {
      AliEmcalJet* jetTmp = dynamic_cast<AliEmcalJet*>(jetArray->At(index)); 
      if(jetTmp)
	dR = GetJetDistance(jet,jetTmp);
    }
  return index;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetEmbed::CalculateDPhi(const Double_t phi1, const Double_t phi2)
{
  Double_t dPhi = phi1-phi2;
  if(dPhi>2*pi)  dPhi -= 2*pi;
  if(dPhi<-2*pi) dPhi += 2*pi;
  if(dPhi<-0.5*pi) dPhi += 2*pi;
  if(dPhi>1.5*pi)  dPhi -= 2*pi;
  return dPhi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetEmbed::GetDPhi(const Double_t phi1, const Double_t phi2)
{
  Double_t dPhi = TMath::Abs(phi1-phi2);
  if(dPhi>2*pi) dPhi -= 2*pi;
  if(dPhi>pi)   dPhi = 2*pi - dPhi;
  return dPhi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetEmbed::GetJetDistance(const AliEmcalJet *jet1, const AliEmcalJet* jet2)
{
  Double_t dPhi = GetDPhi(jet1->Phi(),jet2->Phi());
  Double_t dEta = jet1->Eta()-jet2->Eta();
  return TMath::Sqrt(dPhi*dPhi+dEta*dEta);
}

//________________________________________________________________________
void AliAnalysisTaskHJetEmbed::RunQA()
{
  if(!fPLJetArray)
    {
      AliWarning(Form("Particle-level jet array is not available: %s\n",fPLJetArrName.Data()));
    }
  else
    {
      Int_t nPLJets = fPLJetArray->GetEntries();
      for(Int_t ij=0; ij<nPLJets; ij++)
	{
	  AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fPLJetArray->At(ij));
	  if(!jet) continue;
	  if(!IsGoodJet(jet)) continue; // eta cut
	  Double_t jetPt = jet->Pt();
	  Double_t fill[] = {jetPt, fCentrality, static_cast<Double_t>(fPtHardBin)};
	  fhPLJetPtVsCent[fTriggerType]->Fill(fill);
	  AliDebug(5, Form("PL jet %d has (pt,eta,phi) = (%2.2f,%2.2f,%2.2f)",ij,jetPt,jet->Eta(),jet->Phi()));
	}
    }

  if(!fDLJetArray)
    {
      AliWarning(Form("Detector-level jet array is not available: %s\n",fDLJetArrName.Data()));
    }
  else
    {
      Int_t nDLJets = fDLJetArray->GetEntries();
      for(Int_t ij=0; ij<nDLJets; ij++)
	{
	  AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fDLJetArray->At(ij));
	  if(!jet) continue;
	  if(!IsGoodJet(jet)) continue; // eta cut
	  Double_t jetPt = jet->Pt();
	  Double_t fill[] = {jetPt, fCentrality, static_cast<Double_t>(fPtHardBin)};
	  fhDLJetPtVsCent[fTriggerType]->Fill(fill);
	  AliDebug(5, Form("DL jet %d has pt = %2.2f",ij,jetPt));
	}
    }
}


//________________________________________________________________________
Bool_t AliAnalysisTaskHJetEmbed::AcceptTrack(const AliVParticle *track)
{
  if(track->Pt()<fMinTrkPt || track->Pt()>fMaxTrkPt) return kFALSE;
  if(track->Eta()<fMinTrkEta || track->Eta()>fMaxTrkEta) return kFALSE;
  if(track->Phi()<fMinTrkPhi || track->Phi()>fMaxTrkPhi) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetEmbed::IsGoodJet(const AliEmcalJet* jet)
{
  Double_t etaCut = (0.9-fRadius>0.5)?0.5:0.9-fRadius;
  if(TMath::Abs(jet->Eta())>etaCut) return kFALSE;
  return kTRUE;
}

//
//________________________________________________________________________
//
void AliAnalysisTaskHJetEmbed::PrintConfig()
{
  const char *decision[2] = {"no","yes"};
  const char *TTtype[2] = {"Single inclusive","Leading"};
  printf("\n\n===== h-jet analysis configuration =====\n");
  printf("Input event type: %s - %s\n",fCollisionSystem.Data(),fPeriod.Data());
  printf("Track pt range: %2.2f < pt < %2.2f\n",fMinTrkPt, fMaxTrkPt);
  printf("Track eta range: %2.1f < eta < %2.1f\n",fMinTrkEta, fMaxTrkEta);
  printf("Track phi range: %2.0f < phi < %2.0f\n",fMinTrkPhi*TMath::RadToDeg(),fMaxTrkPhi*TMath::RadToDeg());
  printf("TT type: %s\n", TTtype[fTTtype]);
  for(Int_t i=0; i<kNTT; i++)
    printf("TT range %d:  %2.0f < pt < %2.0f\n", i+1, fMinTTPt[i], fMaxTTPt[i]);
  printf("Run QA: %s\n",decision[fRunQA]);
  printf("Run particle level: %s\n",decision[fRunPL]);
  printf("Run detector level: %s\n",decision[fRunDL]);
  printf("Run h+jet: %s\n",decision[fRunHJet]);
  printf("Run matching: %s\n",decision[fRunMatch]);
  printf("=======================================\n\n");
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetEmbed::GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz)
{
  return (trkPx*jetPx+trkPy*jetPy+trkPz*jetPz)/(jetPx*jetPx+jetPy*jetPy+jetPz*jetPz);
}

//________________________________________________________________________
void AliAnalysisTaskHJetEmbed::Terminate(Option_t *) 
{
  // Called once at the end of the query
}
