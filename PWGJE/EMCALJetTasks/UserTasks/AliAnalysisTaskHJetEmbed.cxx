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

#include "AliEmcalJet.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

ClassImp(AliAnalysisTaskHJetEmbed)

const Double_t pi = TMath::Pi();
const Double_t areaCut[4] = {0.1, 0.23, 0.4, 0.63};

//________________________________________________________________________
AliAnalysisTaskHJetEmbed::AliAnalysisTaskHJetEmbed() : 
  AliAnalysisTaskSE(), 
  fVerbosity(0), fAnaType(1), fPeriod("lhc11h"), fCollisionSystem("PbPb"),
  fEvent(0), fTriggerType(-1), fCentrality(-1), fMaxVtxZ(10),
  fMCParticleArrName(""), fMCParticleArray(0x0), 
  fTrackArrName(""), fTrackArray(0x0), fTriggerTrkIndex(-1), 
  fMinTrkPt(0.15), fMaxTrkPt(1e4), fMinTrkEta(-0.9), fMaxTrkEta(0.9), fMinTrkPhi(0), fMaxTrkPhi(2*pi), 
  fRadius(0.4), fJetArrName(""), fPLJetArrName(""), fDLJetArrName(""),
  fJetArray(0x0), fPLJetArray(0x0), fDLJetArray(0x0),
  fRhoName(""), fRho(0x0), fRhoValue(0), 
  fPtHardBinParam(0), fPtHardBin(-1), 
  fRunQA(kTRUE), fRunHJet(kTRUE), fRunMatch(kTRUE),
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
}

//________________________________________________________________________
AliAnalysisTaskHJetEmbed::AliAnalysisTaskHJetEmbed(const char *name) : 
  AliAnalysisTaskSE(name), 
  fVerbosity(0), fAnaType(1), fPeriod("lhc11h"), fCollisionSystem("PbPb"),
  fEvent(0), fTriggerType(-1), fCentrality(-1), fMaxVtxZ(10),
  fMCParticleArrName(""), fMCParticleArray(0x0), 
  fTrackArrName(""), fTrackArray(0x0), fTriggerTrkIndex(-1), 
  fMinTrkPt(0.15), fMaxTrkPt(1e4), fMinTrkEta(-0.9), fMaxTrkEta(0.9), fMinTrkPhi(0), fMaxTrkPhi(2*pi), 
  fRadius(0.4), fJetArrName(""), fPLJetArrName(""), fDLJetArrName(""),
  fJetArray(0x0), fPLJetArray(0x0), fDLJetArray(0x0),
  fRhoName(""), fRho(0x0), fRhoValue(0), 
  fPtHardBinParam(0), fPtHardBin(-1), 
  fRunQA(kTRUE), fRunHJet(kTRUE), fRunMatch(kTRUE),
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
  const Int_t dimTTqa = 5;
  const Int_t nBinsTTqa[dimTTqa]     = {nTrkPtBins, nTrkPtBins, nTrkPtBins, 30, 11};
  const Double_t lowBinTTqa[dimTTqa] = {lowTrkPtBin,lowTrkPtBin,lowTrkPtBin, 0,  0};
  const Double_t hiBinTTqa[dimTTqa]  = {upTrkPtBin, upTrkPtBin, upTrkPtBin, 30, 11};  

  const Int_t dimTT = 3;
  const Int_t nBinsTT[dimTT]     = {nTrkPtBins,  30, 11};
  const Double_t lowBinTT[dimTT] = {lowTrkPtBin, 0,  0};
  const Double_t hiBinTT[dimTT]  = {upTrkPtBin,  30, 11};  

  const Int_t dimHJet = 6;
  const Int_t nBinsHJet[dimHJet]     = {nTrkPtBins,  nJetPtBins,  144,            8,   30, 11};
  const Double_t lowBinHJet[dimHJet] = {lowTrkPtBin, lowJetPtBin, -0.5*pi+pi/144, 0,   0,  0};
  const Double_t hiBinHJet[dimHJet]  = {upTrkPtBin,  upJetPtBin,  1.5*pi+pi/144,  0.8, 30, 11};  

  // Match
  const Int_t dimMthPt = 6;
  const Int_t nBinsMthPt[dimMthPt]     = {nJetPtBins,  nJetPtBins,  20,    20, 30, 11};
  const Double_t lowBinMthPt[dimMthPt] = {lowJetPtBin, lowJetPtBin, -0.95, 0,  0,  0};
  const Double_t hiBinMthPt[dimMthPt]  = {upJetPtBin,  upJetPtBin,  1.05,  1,  30, 11};

  const Int_t dimMthPhi = 5;
  const Int_t nBinsMthPhi[dimMthPhi]     = {nJetPtBins,  181,                       20, 30, 11};
  const Double_t lowBinMthPhi[dimMthPhi] = {lowJetPtBin, -pi/2-TMath::DegToRad()/2, 0,  0,  0};
  const Double_t hiBinMthPhi[dimMthPhi]  = {upJetPtBin,  pi/2+TMath::DegToRad()/2,  1,  30, 11};

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

	  fhTTPtQA[i] = new THnSparseF(Form("%s_fhTTPtQA",triggerName[i]),Form("PL p_{T} vs DL p_{T} vs embed p_{T} vs centrality vs pT hard bin;p_{T,TT}^{PL} (GeV/c);p_{T,TT}^{DL} (GeV/c);p_{T,TT}^{embed} (GeV/c);centrality;pT hard bin"),dimTTqa,nBinsTTqa,lowBinTTqa,hiBinTTqa);
	  fOutputList->Add(fhTTPtQA[i]);

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

	  fhJetPhiGeoMatch[i] = new THnSparseF(Form("%s_fhJetPhiGeoMatch",triggerName[i]),Form("Embed: generated p_{T,jet} vs #Delta#varphi vs dR vs centrality vs pT hard bin (R=%1.1f);p_{T,jet}^{gen} (GeV/c);#Delta#varphi;dR;centrality;pT hard bin",fRadius),dimMthPhi,nBinsMthPhi,lowBinMthPhi,hiBinMthPhi);
	  fOutputList->Add(fhJetPhiGeoMatch[i]);

	  fhJetPhiEnMatch[i] = new THnSparseF(Form("%s_fhJetPhiEnMatch",triggerName[i]),Form("Embed: generated p_{T,jet} vs #Delta#varphi vs dR vs centrality vs pT hard bin (R=%1.1f);p_{T,jet}^{gen} (GeV/c);#Delta#varphi;dR;centrality;pT hard bin",fRadius),dimMthPhi,nBinsMthPhi,lowBinMthPhi,hiBinMthPhi);
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

  PrintConfig();
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHJetEmbed::UserExec(Option_t *) 
{  
  // Main loop, called for each event.

  fTriggerType = -1;

  fEvent = InputEvent();
  if (!fEvent) 
    {
      AliError("Input event not available");
      return;
    }

  if(!fPtHardBinParam)
    {
      // Get embedded pt hard bin number
      fPtHardBinParam = static_cast<TParameter<int>*>(fEvent->FindListObject("PYTHIAPtHardBin"));
      if(!fPtHardBinParam)
	{
	  AliError("The object for pt hard bin information is not available!");
	  return;
	}
    }
  fPtHardBin = fPtHardBinParam->GetVal();
  AliDebug(2,Form("Embed pt hard bin: %d\n",fPtHardBin));
  if(fPtHardBin<0) return;

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

  const Int_t Ntracks = fTrackArray->GetEntries();
  Int_t nLabel = 0;
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) 
    {
      AliVParticle *t = static_cast<AliVParticle*>(fTrackArray->At(iTracks));
      if (!t) continue;
      if(t->GetLabel()!=0) nLabel ++;
    }

  // Get MC particle array
  if (!fMCParticleArray) 
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
  if (!fPLJetArray && !fPLJetArrName.IsNull())
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
  if (!fDLJetArray && !fDLJetArrName.IsNull())
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
  if(fRunHJet) RunHJet();
  if(fRunMatch) RunMatch();

  PostData(1, fOutputList);
  return;
}


//________________________________________________________________________
void AliAnalysisTaskHJetEmbed::RunHJet()
{
  // Find trigger track on particle level
  const Int_t nParticles = fMCParticleArray->GetEntries();
  Double_t maxPLPt = -1;
  Int_t indexPL = -1;
  for(Int_t iPart=0; iPart<nParticles; iPart++)
    {
      AliVParticle *t = static_cast<AliVParticle*>(fMCParticleArray->At(iPart));
      if(!t || t->Charge()==0) continue;
      if(!AcceptTrack(t)) continue;
      if(maxPLPt<t->Pt())
	{
	  maxPLPt = t->Pt();
	  indexPL = iPart;
	}
    }
  Double_t fill[] = {maxPLPt, fCentrality, fPtHardBin };
  fhPLTT[fTriggerType]->Fill(fill);
  

  // Find trigger track on detector level and after embedding
  const Int_t Ntracks = fTrackArray->GetEntries();
  Double_t maxDLPt = 0, maxEmbPt = 0;
  Int_t indexDL = -1, indexEmb = -1;
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) 
    {
      AliVParticle *t = static_cast<AliVParticle*>(fTrackArray->At(iTracks));
      if(!t || t->Charge()==0) continue;
      if(!AcceptTrack(t)) continue;

      if(t->GetLabel()!=0) 
	{
	  if(maxDLPt<t->Pt())
	    {
	      maxDLPt = t->Pt();
	      indexDL = iTracks;
	    }
	}
      else
	{
	  if(maxEmbPt<t->Pt())
	    {
	      maxEmbPt = t->Pt();
	      indexEmb = iTracks;
	    }
	}
    }
  Double_t fill1[] = {maxDLPt, fCentrality, fPtHardBin };
  fhDLTT[fTriggerType]->Fill(fill1);
  Double_t fill2[] = {maxEmbPt, fCentrality, fPtHardBin };
  fhTTPt[fTriggerType]->Fill(fill2);
  Double_t fill3[] = {maxPLPt, maxDLPt, maxEmbPt, fCentrality, fPtHardBin };
  fhTTPtQA[fTriggerType]->Fill(fill3);
  AliDebug(5,Form("Leading indices: PL=%d, DL=%d, Emb=%d\n",indexPL,indexDL,indexEmb));

  // Run h+jet
  if(indexPL>-1) FillHJetCor(fMCParticleArray, indexPL, fPLJetArray, fhPLHJet[fTriggerType], kFALSE);
  if(indexDL>-1) 
    {
      FillHJetCor(fTrackArray, indexDL, fDLJetArray, fhDLHJet[fTriggerType], kFALSE);
      FillHJetCor(fTrackArray, indexDL, fJetArray, fhHJet[fTriggerType], kTRUE);  
    }
}

//________________________________________________________________________
void  AliAnalysisTaskHJetEmbed::FillHJetCor(const TClonesArray *tracks, const Int_t leadingIndex, const TClonesArray *jetArray, THnSparse *hn, Bool_t isBkg)
{
  if(leadingIndex<0) return;

  AliVParticle *tt = (AliVParticle*) tracks->At(leadingIndex);
  Double_t triggerPt = tt->Pt();
  Double_t triggerPhi = tt->Phi();
  if(triggerPhi<0) triggerPhi += 2*pi;

  Int_t nJets = jetArray->GetEntries();
  for(Int_t ij=0; ij<nJets; ij++)
    {
      AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(jetArray->At(ij));
      if(!IsGoodJet(jet)) continue; // eta cut
      Double_t jetPhi = jet->Phi();
      Double_t jetPt  = jet->Pt();
      Double_t jetArea = jet->Area();
      Double_t dPhi = CalculateDPhi(triggerPhi,jetPhi);
      Double_t fill[] = {triggerPt,jetPt-jetArea*fRhoValue,dPhi,jetArea,fCentrality,fPtHardBin};
      if(!isBkg) fill[1] = jetPt; 
      hn->Fill(fill);
    }
}

//________________________________________________________________________
void AliAnalysisTaskHJetEmbed::RunMatch()
{
  if(!fPLJetArray || !fJetArray)
    {
      AliWarning("Jet array is not available.");
      return;
    }
  Double_t dR = 999;
  Int_t nJets = fJetArray->GetEntries();
  for(Int_t ij=0; ij<nJets; ij++)
    {
      AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fJetArray->At(ij));
      if(!IsGoodJet(jet)) continue; // eta cut
      Double_t jetPt = jet->Pt();
      
      // Geometrial matching
      Int_t mthJetIndexGeo = FindGeoMatchedJet(jet,fPLJetArray,dR);
      if(mthJetIndexGeo>-1)
	{
	  AliEmcalJet* jetMthGeo = dynamic_cast<AliEmcalJet*>(fPLJetArray->At(mthJetIndexGeo));
	  Int_t dataJetIndex = FindGeoMatchedJet(jetMthGeo,fJetArray,dR);
	  if(dataJetIndex==ij) // one-to-one match
	    {
	      Double_t jetPtMthGeo = jetMthGeo->Pt();
	      Double_t fill1[] = {jetPtMthGeo,jetPt,(jetPtMthGeo-jetPt)/jetPtMthGeo,dR,fCentrality, fPtHardBin};
	      fhJetPtGeoMatch[fTriggerType]->Fill(fill1);
	      Double_t fill2[] = {jetPtMthGeo,jetMthGeo->Phi()-jet->Phi(),dR,fCentrality, fPtHardBin};
	      fhJetPhiGeoMatch[fTriggerType]->Fill(fill2);
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
Int_t AliAnalysisTaskHJetEmbed::FindEnergyMatchedJet(const AliEmcalJet* jet, const TClonesArray *jetArray, Double_t &dR)
{
  dR = 999;
  if(!jetArray) return -1;
  
  Int_t index = -1;
  Int_t nJets = jetArray->GetEntries();
  Double_t fMin = 0;
  Int_t nJetC = (Int_t)jet->GetNumberOfConstituents();
  for(Int_t ij=0; ij<nJets; ij++)
    {
      AliEmcalJet* jetTmp = dynamic_cast<AliEmcalJet*>(jetArray->At(ij)); 
      if(TMath::Abs(jetTmp->Eta())>1) continue; // Generous eta cut
      Double_t jetPt = jet->Pt();
      Int_t nc = (Int_t)jetTmp->GetNumberOfConstituents();
      Double_t sumPt = 0;
      for(Int_t ic=0; ic<nc; ic++)
	{
	  AliVParticle *part = jetTmp->TrackAt(ic, fMCParticleArray);
	  for(Int_t ijc=0; ijc<nJetC; ijc++)
	    {
	      AliVParticle *track = jet->TrackAt(ijc, fTrackArray);
	      if(track->GetLabel()==part->GetLabel())
		sumPt += part->Pt();
	    }
	}
      Double_t frac = sumPt/jetPt;
      if(frac>fMin)
	{
	  fMin = frac;
	  index = ij;
	}
    }

  if(index>0)
    {
      AliEmcalJet* jetTmp = dynamic_cast<AliEmcalJet*>(jetArray->At(index)); 
      Double_t dPhi = GetDPhi(jet->Phi(),jetTmp->Phi());
      Double_t dEta = jet->Eta()-jetTmp->Eta();
      dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
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
	  if(!IsGoodJet(jet)) continue; // eta cut
	  Double_t jetPt = jet->Pt();
	  Double_t fill[] = {jetPt, fCentrality, fPtHardBin};
	  fhPLJetPtVsCent[fTriggerType]->Fill(fill);
	  AliDebug(5, Form("PL jet %d has (pt,eta,phi) = (%2.2f,%2.2f,%2.2f)",ij,jetPt,jet->Eta(),jet->Phi()));
	}
    }

  if(!fDLJetArray)
    {
      AliWarning(Form("Particle-level jet array is not available: %s\n",fDLJetArrName.Data()));
    }
  else
    {
      Int_t nDLJets = fDLJetArray->GetEntries();
      for(Int_t ij=0; ij<nDLJets; ij++)
	{
	  AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fDLJetArray->At(ij));
	  if(!IsGoodJet(jet)) continue; // eta cut
	  Double_t jetPt = jet->Pt();
	  Double_t fill[] = {jetPt, fCentrality, fPtHardBin};
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
  printf("\n\n===== h-jet analysis configuration =====\n");
  printf("Input event type: %s - %s\n",fCollisionSystem.Data(),fPeriod.Data());
  printf("Track pt range: %2.2f < pt < %2.2f\n",fMinTrkPt, fMaxTrkPt);
  printf("Track eta range: %2.1f < eta < %2.1f\n",fMinTrkEta, fMaxTrkEta);
  printf("Track phi range: %2.0f < phi < %2.0f\n",fMinTrkPhi*TMath::RadToDeg(),fMaxTrkPhi*TMath::RadToDeg());
  printf("Run QA: %s\n",decision[fRunQA]);
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
