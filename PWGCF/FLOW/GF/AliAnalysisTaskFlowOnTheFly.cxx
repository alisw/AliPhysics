/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskFlowOnTheFly.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"

class AliAnalysisTaskFlowOnTheFly;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskFlowOnTheFly) // classimp: necessary for root

AliAnalysisTaskFlowOnTheFly::AliAnalysisTaskFlowOnTheFly() : AliAnalysisTaskSE(), 
    fFilterbit(96),
    fFilterbitDefault(96),
    fEtaCut(0.8),
    fVtxCut(10.0),
    fVtxCutDefault(10.0),
    fMinPt(0.2),
    fMaxPt(3.0),
	fUseImpactXaxis(false),
	fEventWeightSetToOne(false),
	fAddTPCPileupCuts(false),
	fESDvsTPConlyLinearCut(15000.),
	fUseCL1Centrality(0),
    fTrigger(0),
    fListOfObjects(0),

    fTrackEfficiency(0),
    hTrackEfficiencyRun(0),

    fFlowRunByRunWeights(false),
    fFlowPeriodWeights(false),
    fFlowUse3Dweights(false),
    fFlowWeightsList(nullptr),

    hPhiWeight(0),

    hEventCount(0),
    fCentralityDis(0),

	hPhi(0),
	hPhiBefore(0),
	fEtaDis(0),
	fPtDis(0),

    rand(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskFlowOnTheFly::AliAnalysisTaskFlowOnTheFly(const char* name) : AliAnalysisTaskSE(name),
	fFilterbit(96),
	fEtaCut(0.8),
	fVtxCut(10.0),
	fMinPt(0.2),
	fMaxPt(3.0),
	fUseImpactXaxis(false),
	fEventWeightSetToOne(false),
	fAddTPCPileupCuts(false),
	fESDvsTPConlyLinearCut(15000.),
	fUseCL1Centrality(0),
	fTrigger(0),
	fAliTrigger(0),
	fListOfObjects(0),

	fTrackEfficiency(0),
	hTrackEfficiencyRun(0),

        fFlowRunByRunWeights(false),
        fFlowPeriodWeights(false),
        fFlowUse3Dweights(false),
        fFlowWeightsList(nullptr),

	hPhiWeight(0),

	hEventCount(0),
	fCentralityDis(0),

	hPhi(0),
	hPhiBefore(0),
	fEtaDis(0),
	fPtDis(0),
	rand(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    
	
	//the number of DefineInput should be in line with NUE,NUA
	//use NUA only--> One DefineInput 
	//use NUE and NUA --> Two DefineInput
	DefineInput(1, TFile::Class());
	DefineInput(2, TFile::Class());

    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskFlowOnTheFly::~AliAnalysisTaskFlowOnTheFly()
{
    // destructor
    if (fListOfObjects)
		delete fListOfObjects;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowOnTheFly::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fListOfObjects = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fListOfObjects->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


	if(fUseImpactXaxis){
		nn = 1000;
		for (int i = 0; i <= 1000; i++) {
		xbins[i] = 30.0/nn*i;
		}
	}
	else{
		vector<double> cent = {0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0};
		if(centralitymap.empty()) {
			vector<double> b = {0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,13.50,14.51,100.0};
			for(size_t i(0); i<b.size(); ++i) centralitymap[b[i]]=cent[i];
		}
		nn = 10;
		for (int i = 0; i <= nn; i++) {
			xbins[i] = cent[i];
		}
	}

    fIP = new TH1D("fIP","Impact parameter",1000,0.0,30.0);
	fListOfObjects->Add(fIP);

    hEventCount = new TH1D("hEventCount", "; centrality;;", 7, 0, 7);
	fListOfObjects->Add(hEventCount);

	fCentralityDis = new TH1F("fCentralityDis", "centrality distribution; centrality; Counts", 100, 0, 100);
	fListOfObjects->Add(fCentralityDis);

	hPhiBefore = new TH1D("hPhiBefore", "phi distribution before the weight correction", 60, 0, 2*3.1415926);
  	fListOfObjects->Add(hPhiBefore);
  	hPhi  = new TH1D("hPhi", "phi distribution after the weight correction", 60, 0, 2*3.1415926);
  	fListOfObjects->Add(hPhi);
	fEtaDis = new TH1D("hEtaDis", "eta distribution", 100, -2, 2);
	fListOfObjects->Add(fEtaDis);
	fPtDis = new TH1D("hPtDis", "pt distribution", 100, 0, 5);
	fListOfObjects->Add(fPtDis);

    // Int_t inSlotCounter=1;
	// if(fNUA) {
    //             fFlowWeightsList = (TList*) GetInputData(inSlotCounter);
	// 	inSlotCounter++;
	// };
	// if(fNUE) {
	// 	fTrackEfficiency = (TList*)GetInputData(inSlotCounter);
	// 	inSlotCounter++;
	// };

    // Physics profiles
	//	NL response
	InitProfile(multProfile, "");
    for (int i = 0; i < 30; i++) InitProfile(multProfile_bin[i], Form("_%d", i));
                                        // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!
    
    PostData(1, fListOfObjects);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
    // Post output data.
	
}

//_________________________________________________________________
void AliAnalysisTaskFlowOnTheFly::NotifyRun() {
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowOnTheFly::UserExec(Option_t *)
{
	
	//hEventCount->Fill("after fEventCuts", 1.);
	hEventCount->GetXaxis()->SetBinLabel(1,"Loop Number");
    hEventCount->Fill(0.5);

	// user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain   
  
    bootstrap_value = rand.Integer(30);	

	//..all charged particles
	//AnalyzeAOD calculate Q vector and Profile。
	//All the Profiles are stored in fListOfObjects 
	//Output in PostData
	ProcessOnTheFly();

	hEventCount->GetXaxis()->SetBinLabel(2,"after ProcessOnTheFly");
	hEventCount->Fill(1.5);

    PostData(1, fListOfObjects);                          // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
    // Post output data.
	
}
AliMCEvent *AliAnalysisTaskFlowOnTheFly::getMCEvent() {
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!ev) { AliFatal("MC event not found!"); return 0; }
  AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!"); return 0; }
  AliCollisionGeometry* headerH;
  TString genName;
  TList *ltgen = (TList*)ev->GetCocktailList();
  if (ltgen) {
  for(auto&& listObject: *ltgen){
    genName = Form("%s",listObject->GetName());
    if (genName.Contains("Hijing")) {
      headerH = dynamic_cast<AliCollisionGeometry*>(listObject);
      break;
      }
    }
  }
  else headerH = dynamic_cast<AliCollisionGeometry*>(ev->GenEventHeader());
  if(headerH){
      fImpactParameterMC = headerH->ImpactParameter();
  }
  return ev;
}
double AliAnalysisTaskFlowOnTheFly::getAMPTCentrality()
{
  vector<double> b;
  if(centralitymap.empty()) AliFatal("Centralitymap is empty!");
  for (auto const& element : centralitymap) b.push_back(element.first);
  vector<double>::iterator it = upper_bound(b.begin(),b.end(),fImpactParameterMC);
  double l_cent = (fImpactParameterMC<0)?-1.0:(centralitymap[b[it-b.begin()]]+centralitymap[b[it-b.begin()-1]])/2.0;
  return l_cent;
}
void AliAnalysisTaskFlowOnTheFly::ProcessOnTheFly(){
	fMCEvent = getMCEvent();
	fIP->Fill(fImpactParameterMC);
	//Double_t l_Cent 
	fCurrCentrality = getAMPTCentrality();
	fCentralityDis->Fill(fCurrCentrality);
	Int_t nTracks = fMCEvent->GetNumberOfPrimaries();
	if(nTracks < 1) { return; }

	double NtrksBefore = 0;
    NtrksAfter = 0;
	NtrksAfterGap0M = 0;
	NtrksAfterGap0P = 0;
	NtrksAfterGap2M = 0;
	NtrksAfterGap2P = 0;
	NtrksAfterGap4M = 0;
	NtrksAfterGap4P = 0;
	NtrksAfterGap6M = 0;
	NtrksAfterGap6P = 0;
	NtrksAfterGap8M = 0;
	NtrksAfterGap8P = 0;	
	NtrksAfterGap10M = 0;
	NtrksAfterGap10P = 0;
	NtrksAfterGap12M = 0;
	NtrksAfterGap12P = 0;
	NtrksAfterGap14M = 0;
	NtrksAfterGap14P = 0;
	NtrksAfter3subL = 0;
	NtrksAfter3subM = 0;
	NtrksAfter3subR = 0;

	double Qcos[20][20] = {0};
	double Qsin[20][20] = {0};
	double QcosGap0M[20][20] = {0};
	double QsinGap0M[20][20] = {0};
	double QcosGap0P[20][20] = {0};
	double QsinGap0P[20][20] = {0};
	double QcosGap2M[20][20] = {0};
	double QsinGap2M[20][20] = {0};
	double QcosGap2P[20][20] = {0};
	double QsinGap2P[20][20] = {0};
	double QcosGap4M[20][20] = {0};
	double QsinGap4M[20][20] = {0};
	double QcosGap4P[20][20] = {0};
	double QsinGap4P[20][20] = {0};
	double QcosGap6M[20][20] = {0};
	double QsinGap6M[20][20] = {0};
	double QcosGap6P[20][20] = {0};
	double QsinGap6P[20][20] = {0};
	double QcosGap8M[20][20] = {0};
	double QsinGap8M[20][20] = {0};
	double QcosGap8P[20][20] = {0};
	double QsinGap8P[20][20] = {0};
	double QcosGap10M[20][20] = {0};
	double QsinGap10M[20][20] = {0};
	double QcosGap10P[20][20] = {0};
	double QsinGap10P[20][20] = {0};
	double QcosGap12M[20][20] = {0};
	double QsinGap12M[20][20] = {0};
	double QcosGap12P[20][20] = {0};
	double QsinGap12P[20][20] = {0};
	double QcosGap14M[20][20] = {0};
	double QsinGap14M[20][20] = {0};
	double QcosGap14P[20][20] = {0};
	double QsinGap14P[20][20] = {0};
	double QcosSubLeft[20][20] = {0};
	double QsinSubLeft[20][20] = {0};
	double QcosSubMiddle[20][20] = {0};
	double QsinSubMiddle[20][20] = {0};
	double QcosSubRight[20][20] = {0};
	double QsinSubRight[20][20] = {0};

	for(Int_t i=0;i<nTracks;i++){
		AliMCParticle* lPart = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(i));
		if(!lPart) { continue; };
		if(!lPart->IsPhysicalPrimary()) continue;
		Double_t l_pt=lPart->Pt();
		Double_t l_phi=lPart->Phi();
		Double_t l_eta=lPart->Eta();

		//Pt, Eta cut
		if(l_pt < fMinPt) continue;
		if(l_pt > fMaxPt) continue;
		if(TMath::Abs(l_eta) > fEtaCut) continue;
		if (lPart->Charge() == 0) continue;

		//..get phi-weight for NUA correction
		double weight = 1;
		double weightPt = 1;

		hPhiBefore->Fill(lPart->Phi());
    	fPtDis->Fill(lPart->Pt());
    	fEtaDis->Fill(lPart->Eta());
    	hPhi->Fill(lPart->Phi(), weight*weightPt);

		//..calculate Q-vectors
		//..no eta gap
		for(int iharm=0; iharm<20; iharm++)
		{
			for(int ipow=0; ipow<20; ipow++)
			{
				Qcos[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
				Qsin[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
			}
		}

		//..Gap > 0.
		if(lPart->Eta() < 0.)
		{
			NtrksAfterGap0M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() > 0.)
		{
			NtrksAfterGap0P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}

		//..Gap > 0.2
		if(lPart->Eta() < -0.1)
		{
			NtrksAfterGap2M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() > 0.1)
		{
			NtrksAfterGap2P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}

		//..Gap > 0.4
		if(lPart->Eta() < -0.2)
		{
			NtrksAfterGap4M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() > 0.2)
		{
			NtrksAfterGap4P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}

		//..Gap > 0.6
		if(lPart->Eta() < -0.3)
		{
			NtrksAfterGap6M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap6M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap6M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() > 0.3)
		{
			NtrksAfterGap6P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap6P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap6P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}

		//..Gap > 0.8
		if(lPart->Eta() < -0.4)
		{
			NtrksAfterGap8M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() > 0.4)
		{
			NtrksAfterGap8P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}

		//..Gap > 1.0
		if(lPart->Eta() < -0.5)
		{
			NtrksAfterGap10M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap10M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap10M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() > 0.5)
		{
			NtrksAfterGap10P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap10P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap10P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}

		//..Gap > 1.2
		if(lPart->Eta() < -0.6)
		{
			NtrksAfterGap12M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap12M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap12M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() > 0.6)
		{
			NtrksAfterGap12P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap12P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap12P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}

		//..Gap > 1.4
		if(lPart->Eta() < -0.7)
		{
			NtrksAfterGap14M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() > 0.7)
		{
			NtrksAfterGap14P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		//..3-subevent method
		if(lPart->Eta() < -0.4)
		{//..left part
			NtrksAfter3subL += 1;
			for(int iharm=0; iharm<10; iharm++)
			{
				for(int ipow=0; ipow<10; ipow++)
				{
					QcosSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() >= -0.4 && lPart->Eta() <= 0.4)
		{//..middle part
			NtrksAfter3subM += 1;
			for(int iharm=0; iharm<10; iharm++)
			{
				for(int ipow=0; ipow<10; ipow++)
				{
					QcosSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
		if(lPart->Eta() > 0.4)
		{//..right part
			NtrksAfter3subR += 1;
			for(int iharm=0; iharm<10; iharm++)
			{
				for(int ipow=0; ipow<10; ipow++)
				{
					QcosSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*lPart->Phi());
					QsinSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*lPart->Phi());
				}
			}
		}
	}// end loop of all track

	//............................
	//..GENERIC FRAMEWORK RP
	//............................

	//..calculate Q-vector for each harmonics n and power p
        correlator.FillQVector(correlator.Qvector, Qcos, Qsin); 
		correlator.FillQVector(correlator.Qvector0M, QcosGap0M, QsinGap0M); 
        correlator.FillQVector(correlator.Qvector0P, QcosGap0P, QsinGap0P); 
		correlator.FillQVector(correlator.Qvector2M, QcosGap2M, QsinGap2M); 
        correlator.FillQVector(correlator.Qvector2P, QcosGap2P, QsinGap2P);
		correlator.FillQVector(correlator.Qvector4M, QcosGap4M, QsinGap4M); 
        correlator.FillQVector(correlator.Qvector4P, QcosGap4P, QsinGap4P);
		correlator.FillQVector(correlator.Qvector6M, QcosGap6M, QsinGap6M); 
        correlator.FillQVector(correlator.Qvector6P, QcosGap6P, QsinGap6P);
		correlator.FillQVector(correlator.Qvector8M, QcosGap8M, QsinGap8M); 
        correlator.FillQVector(correlator.Qvector8P, QcosGap8P, QsinGap8P);
        correlator.FillQVector(correlator.Qvector10M, QcosGap10M, QsinGap10M); 
        correlator.FillQVector(correlator.Qvector10P, QcosGap10P, QsinGap10P); 
		correlator.FillQVector(correlator.Qvector12M, QcosGap12M, QsinGap12M); 
        correlator.FillQVector(correlator.Qvector12P, QcosGap12P, QsinGap12P); 
        correlator.FillQVector(correlator.Qvector14M, QcosGap14M, QsinGap14M); 
        correlator.FillQVector(correlator.Qvector14P, QcosGap14P, QsinGap14P); 
        correlator.FillQVector(correlator.QvectorSubLeft, QcosSubLeft, QsinSubLeft); 
        correlator.FillQVector(correlator.QvectorSubRight, QcosSubRight, QsinSubRight); 
        correlator.FillQVector(correlator.QvectorSubMiddle, QcosSubMiddle, QsinSubMiddle); 

	// CalculateProfile(centProfile, cent);
	if (fUseImpactXaxis) {
	    CalculateProfile(multProfile, fImpactParameterMC);
	    CalculateProfile(multProfile_bin[bootstrap_value], fImpactParameterMC);
        } else {
	    CalculateProfile(multProfile, fCurrCentrality);
	    CalculateProfile(multProfile_bin[bootstrap_value], fCurrCentrality);
        }
}

//_____________________________________________________________________________
void AliAnalysisTaskFlowOnTheFly::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//=============================================================================
void AliAnalysisTaskFlowOnTheFly::InitProfile(PhysicsProfileFlowOnTheFly& multProfile, TString label) {

	
	for(int h=0; h<5; h++)
	{
		//Fill 2 correlation. 
		//fChcn2[i]: v_(i+2){2}
		multProfile.fChcn2[h] = new TProfile(Form("fChc%d{2}%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2[h]);

		multProfile.fChcn2_Gap0[h] = new TProfile(Form("fChc%d{2}_Gap0%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap0[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap0[h]);

		multProfile.fChcn2_Gap2[h] = new TProfile(Form("fChc%d{2}_Gap2%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap2[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap2[h]);

		multProfile.fChcn2_Gap4[h] = new TProfile(Form("fChc%d{2}_Gap4%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap4[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap4[h]);

		multProfile.fChcn2_Gap6[h] = new TProfile(Form("fChc%d{2}_Gap6%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap6[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap6[h]);

		multProfile.fChcn2_Gap8[h] = new TProfile(Form("fChc%d{2}_Gap8%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap8[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap8[h]);

		multProfile.fChcn2_Gap10[h] = new TProfile(Form("fChc%d{2}_Gap10%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap10[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap10[h]);

		multProfile.fChcn2_Gap12[h] = new TProfile(Form("fChc%d{2}_Gap12%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap12[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap12[h]);

		multProfile.fChcn2_Gap14[h] = new TProfile(Form("fChc%d{2}_Gap14%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap14[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap14[h]);

		multProfile.fChcn2_3subLM[h] = new TProfile(Form("fChc%d{2}_3subLM%s", h+2, label.Data()), "<<2>> Re for 3-subevent method, left+middle; # of tracks", nn, xbins);
		multProfile.fChcn2_3subLM[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_3subLM[h]);

		multProfile.fChcn2_3subRM[h] = new TProfile(Form("fChc%d{2}_3subRM%s", h+2, label.Data()), "<<2>> Re for 3-subevent method, right+middle; # of tracks", nn, xbins);
		multProfile.fChcn2_3subRM[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_3subRM[h]);

		multProfile.fChcn2_3subLR[h] = new TProfile(Form("fChc%d{2}_3subLR%s", h+2, label.Data()), "<<2>> Re for 3-subevent method, left+right; # of tracks", nn, xbins);
		multProfile.fChcn2_3subLR[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_3subLR[h]);

		multProfile.fChcn4[h] = new TProfile(Form("fChc%d{4}%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4[h]);

		multProfile.fChcn4_Gap0[h] = new TProfile(Form("fChc%d{4}_Gap0%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4_Gap0[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_Gap0[h]);

		multProfile.fChcn4_Gap2[h] = new TProfile(Form("fChc%d{4}_Gap2%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4_Gap2[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_Gap2[h]);

		multProfile.fChcn4_Gap4[h] = new TProfile(Form("fChc%d{4}_Gap4%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4_Gap4[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_Gap4[h]);

		multProfile.fChcn4_Gap6[h] = new TProfile(Form("fChc%d{4}_Gap6%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4_Gap6[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_Gap6[h]);

		multProfile.fChcn4_Gap8[h] = new TProfile(Form("fChc%d{4}_Gap8%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4_Gap8[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_Gap8[h]);

		multProfile.fChcn4_Gap10[h] = new TProfile(Form("fChc%d{4}_Gap10%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4_Gap10[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_Gap10[h]);

		multProfile.fChcn4_Gap12[h] = new TProfile(Form("fChc%d{4}_Gap12%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4_Gap12[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_Gap12[h]);

		multProfile.fChcn4_3subLLMR[h] = new TProfile(Form("fChc%d{4}_3subLLMR%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
		multProfile.fChcn4_3subLLMR[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_3subLLMR[h]);

		multProfile.fChcn4_3subRRML[h] = new TProfile(Form("fChc%d{4}_3subRRML%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
		multProfile.fChcn4_3subRRML[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_3subRRML[h]);

		multProfile.fChcn4_3subMMLR[h] = new TProfile(Form("fChc%d{4}_3subMMLR%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
		multProfile.fChcn4_3subMMLR[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_3subMMLR[h]);

		multProfile.fChcn4_3subGap2[h] = new TProfile(Form("fChc%d{4}_3subGap2%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
		multProfile.fChcn4_3subGap2[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_3subGap2[h]);

		multProfile.fChcn6[h] = new TProfile(Form("fChc%d{6}%s", h+2, label.Data()), "<<6>> Re; # of tracks", nn, xbins);
		multProfile.fChcn6[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn6[h]);

		multProfile.fChcn6_Gap0[h] = new TProfile(Form("fChc%d{6}_Gap0%s", h+2, label.Data()), "<<6>> Re; # of tracks", nn, xbins);
		multProfile.fChcn6_Gap0[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn6_Gap0[h]);

		multProfile.fChcn6_Gap10[h] = new TProfile(Form("fChc%d{6}_Gap10%s", h+2, label.Data()), "<<6>> Re; # of tracks", nn, xbins);
		multProfile.fChcn6_Gap10[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn6_Gap10[h]);

		multProfile.fChcn8[h] = new TProfile(Form("fChc%d{8}%s", h+2, label.Data()), "<<8>> Re; # of tracks", nn, xbins);
		multProfile.fChcn8[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn8[h]);

		multProfile.fChcn8_Gap0[h] = new TProfile(Form("fChc%d{8}_Gap0%s", h+2, label.Data()), "<<8>> Re; # of tracks", nn, xbins);
		multProfile.fChcn8_Gap0[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn8_Gap0[h]);

	} // harmonics

	//Three Correlator

	multProfile.fChc422 = new TProfile(Form("fChc422%s", label.Data()), "", nn, xbins);
	multProfile.fChc422->Sumw2();
	fListOfObjects->Add(multProfile.fChc422);

	multProfile.fChc532 = new TProfile(Form("fChc532%s", label.Data()), "", nn, xbins);
	multProfile.fChc532->Sumw2();
	fListOfObjects->Add(multProfile.fChc532);
	// Gap0
	multProfile.fChc422_Gap0A = new TProfile(Form("fChc422_Gap0A%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap0A->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap0A);

	multProfile.fChc422_Gap0B = new TProfile(Form("fChc422_Gap0B%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap0B->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap0B);

	multProfile.fChc532_Gap0A = new TProfile(Form("fChc532_Gap0A%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap0A->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap0A);

	multProfile.fChc532_Gap0B = new TProfile(Form("fChc532_Gap0B%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap0B->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap0B);

	// Gap2
	multProfile.fChc422_Gap2A = new TProfile(Form("fChc422_Gap2A%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap2A->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap2A);

	multProfile.fChc422_Gap2B = new TProfile(Form("fChc422_Gap2B%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap2B->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap2B);

	multProfile.fChc532_Gap2A = new TProfile(Form("fChc532_Gap2A%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap2A->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap2A);

	multProfile.fChc532_Gap2B = new TProfile(Form("fChc532_Gap2B%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap2B->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap2B);

	// Gap4
	multProfile.fChc422_Gap4A = new TProfile(Form("fChc422_Gap4A%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap4A->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap4A);

	multProfile.fChc422_Gap4B = new TProfile(Form("fChc422_Gap4B%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap4B->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap4B);

	multProfile.fChc532_Gap4A = new TProfile(Form("fChc532_Gap4A%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap4A->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap4A);

	multProfile.fChc532_Gap4B = new TProfile(Form("fChc532_Gap4B%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap4B->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap4B);

	// Gap6
	multProfile.fChc422_Gap6A = new TProfile(Form("fChc422_Gap6A%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap6A->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap6A);

	multProfile.fChc422_Gap6B = new TProfile(Form("fChc422_Gap6B%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap6B->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap6B);

	multProfile.fChc532_Gap6A = new TProfile(Form("fChc532_Gap6A%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap6A->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap6A);

	multProfile.fChc532_Gap6B = new TProfile(Form("fChc532_Gap6B%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap6B->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap6B);

	// Gap8
	multProfile.fChc422_Gap8A = new TProfile(Form("fChc422_Gap8A%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap8A->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap8A);

	multProfile.fChc422_Gap8B = new TProfile(Form("fChc422_Gap8B%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap8B->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap8B);

	multProfile.fChc532_Gap8A = new TProfile(Form("fChc532_Gap8A%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap8A->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap8A);

	multProfile.fChc532_Gap8B = new TProfile(Form("fChc532_Gap8B%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap8B->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap8B);

	// Gap10
	multProfile.fChc422_Gap10A = new TProfile(Form("fChc422_Gap10A%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap10A->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap10A);

	multProfile.fChc422_Gap10B = new TProfile(Form("fChc422_Gap10B%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap10B->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap10B);

	multProfile.fChc532_Gap10A = new TProfile(Form("fChc532_Gap10A%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap10A->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap10A);

	multProfile.fChc532_Gap10B = new TProfile(Form("fChc532_Gap10B%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap10B->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap10B);

	// Gap12
	multProfile.fChc422_Gap12A = new TProfile(Form("fChc422_Gap12A%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap12A->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap12A);

	multProfile.fChc422_Gap12B = new TProfile(Form("fChc422_Gap12B%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap12B->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap12B);

	multProfile.fChc532_Gap12A = new TProfile(Form("fChc532_Gap12A%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap12A->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap12A);

	multProfile.fChc532_Gap12B = new TProfile(Form("fChc532_Gap12B%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap12B->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap12B);

	// SC(n,m): SC(3,2)
	multProfile.fChsc3232 = new TProfile(Form("fChsc3232%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232);

	multProfile.fChsc3232_Gap0 = new TProfile(Form("fChsc3232_Gap0%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_Gap0->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_Gap0);

	multProfile.fChsc3232_Gap2 = new TProfile(Form("fChsc3232_Gap2%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_Gap2->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_Gap2);

	multProfile.fChsc3232_Gap4 = new TProfile(Form("fChsc3232_Gap4%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_Gap4->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_Gap4);

	multProfile.fChsc3232_Gap6 = new TProfile(Form("fChsc3232_Gap6%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_Gap6->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_Gap6);

	multProfile.fChsc3232_Gap8 = new TProfile(Form("fChsc3232_Gap8%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_Gap8->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_Gap8);

	multProfile.fChsc3232_Gap10 = new TProfile(Form("fChsc3232_Gap10%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_Gap10->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_Gap10);

	multProfile.fChsc3232_Gap12 = new TProfile(Form("fChsc3232_Gap12%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_Gap12->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_Gap12);

	multProfile.fChsc3232_3subMMLRA = new TProfile(Form("fChsc3232_3subMMLRA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subMMLRA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subMMLRA);

	multProfile.fChsc3232_3subMMLRB = new TProfile(Form("fChsc3232_3subMMLRB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subMMLRB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subMMLRB);

	multProfile.fChsc3232_3subLLMRA = new TProfile(Form("fChsc3232_3subLLMRA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subLLMRA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subLLMRA);

	multProfile.fChsc3232_3subLLMRB = new TProfile(Form("fChsc3232_3subLLMRB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subLLMRB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subLLMRB);

	multProfile.fChsc3232_3subRRMLA = new TProfile(Form("fChsc3232_3subRRMLA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subRRMLA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subRRMLA);

	multProfile.fChsc3232_3subRRMLB = new TProfile(Form("fChsc3232_3subRRMLB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subRRMLB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subRRMLB);

	// SC(n,m): SC(4,2)
	multProfile.fChsc4242 = new TProfile(Form("fChsc4242%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242);

	multProfile.fChsc4242_Gap0 = new TProfile(Form("fChsc4242_Gap0%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_Gap0->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_Gap0);

	multProfile.fChsc4242_Gap2 = new TProfile(Form("fChsc4242_Gap2%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_Gap2->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_Gap2);

	multProfile.fChsc4242_Gap4 = new TProfile(Form("fChsc4242_Gap4%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_Gap4->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_Gap4);

	multProfile.fChsc4242_Gap6 = new TProfile(Form("fChsc4242_Gap6%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_Gap6->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_Gap6);

	multProfile.fChsc4242_Gap8 = new TProfile(Form("fChsc4242_Gap8%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_Gap8->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_Gap8);

	multProfile.fChsc4242_Gap10 = new TProfile(Form("fChsc4242_Gap10%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_Gap10->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_Gap10);

	multProfile.fChsc4242_Gap12 = new TProfile(Form("fChsc4242_Gap12%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_Gap12->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_Gap12);

	multProfile.fChsc4242_3subMMLRA = new TProfile(Form("fChsc4242_3subMMLRA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subMMLRA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subMMLRA);

	multProfile.fChsc4242_3subMMLRB = new TProfile(Form("fChsc4242_3subMMLRB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subMMLRB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subMMLRB);

	multProfile.fChsc4242_3subLLMRA = new TProfile(Form("fChsc4242_3subLLMRA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subLLMRA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subLLMRA);

	multProfile.fChsc4242_3subLLMRB = new TProfile(Form("fChsc4242_3subLLMRB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subLLMRB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subLLMRB);

	multProfile.fChsc4242_3subRRMLA = new TProfile(Form("fChsc4242_3subRRMLA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subRRMLA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subRRMLA);

	multProfile.fChsc4242_3subRRMLB = new TProfile(Form("fChsc4242_3subRRMLB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subRRMLB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subRRMLB);

	// //5,6 correlation
	// multProfile.fChc5_A42222 = new TProfile("fChc5_A42222", "<<5>> Re; # of tracks", nn, xbins);
	// multProfile.fChc5_A42222->Sumw2();
	// fListOfObjects->Add(multProfile.fChc5_A42222);

	// multProfile.fChc5_A52322 = new TProfile("fChc5_A52322", "<<5>> Re; # of tracks", nn, xbins);
	// multProfile.fChc5_A52322->Sumw2();
	// fListOfObjects->Add(multProfile.fChc5_A52322);

	// multProfile.fChc6_222222 = new TProfile("fChc6_222222", "<<6>> Re; # of tracks", nn, xbins);
	// multProfile.fChc6_222222->Sumw2();
	// fListOfObjects->Add(multProfile.fChc6_222222);

	// multProfile.fChc6_322322 = new TProfile("fChc6_322322", "<<6>> Re; # of tracks", nn, xbins);
	// multProfile.fChc6_322322->Sumw2();
	// fListOfObjects->Add(multProfile.fChc6_322322);

	// //Additional 3,4 correlation
	// multProfile.fChsc6222_Gap0 = new TProfile(Form("fChsc6222_Gap0%s", label.Data()), "# of tracks", nn, xbins);
	// multProfile.fChsc6222_Gap0->Sumw2();
	// fListOfObjects->Add(multProfile.fChsc6222_Gap0);

	// multProfile.fChsc6222_Gap10 = new TProfile(Form("fChsc6222_Gap10%s", label.Data()), "# of tracks", nn, xbins);
	// multProfile.fChsc6222_Gap10->Sumw2();
	// fListOfObjects->Add(multProfile.fChsc6222_Gap10);

	// multProfile.fChsc633_Gap0A = new TProfile(Form("fChsc633_Gap0A%s", label.Data()), "# of tracks", nn, xbins);
	// multProfile.fChsc633_Gap0A->Sumw2();
	// fListOfObjects->Add(multProfile.fChsc633_Gap0A);

	// multProfile.fChsc633_Gap10A = new TProfile(Form("fChsc633_Gap10A%s", label.Data()), "# of tracks", nn, xbins);
	// multProfile.fChsc633_Gap10A->Sumw2();
	// fListOfObjects->Add(multProfile.fChsc633_Gap10A);

}

void AliAnalysisTaskFlowOnTheFly::CalculateProfile(PhysicsProfileFlowOnTheFly& profile, double Ntrks) {
	//..calculate 2-particle correlations
	//..................................
	double Dn2 = correlator.Two(0, 0).Re();
	double Dn2Gap0 = correlator.TwoGap0(0, 0).Re();
	double Dn2Gap2 = correlator.TwoGap2(0, 0).Re();
	double Dn2Gap4 = correlator.TwoGap4(0, 0).Re();
	double Dn2Gap6 = correlator.TwoGap6(0, 0).Re();
	double Dn2Gap8 = correlator.TwoGap8(0, 0).Re();
	double Dn2Gap10 = correlator.TwoGap10(0, 0).Re();
	double Dn2Gap12 = correlator.TwoGap12(0, 0).Re();
	double Dn2Gap14 = correlator.TwoGap14(0, 0).Re();
	double Dn2_3subLM = correlator.Two_3SubLM(0, 0).Re();
	double Dn2_3subRM = correlator.Two_3SubRM(0, 0).Re();
	double Dn2_3subLR = correlator.Two_3SubLR(0, 0).Re();
	double Weight_Dn2 = 1.;
	double Weight_Dn2Gap0 = 1.;
	double Weight_Dn2Gap2 = 1.;
	double Weight_Dn2Gap4 = 1.;
	double Weight_Dn2Gap6 = 1.;
	double Weight_Dn2Gap8 = 1.;
	double Weight_Dn2Gap10 = 1.;
	double Weight_Dn2Gap12 = 1.;
	double Weight_Dn2Gap14 = 1.;
	double Weight_Dn2_3subLM = 1.;
	double Weight_Dn2_3subRM = 1.;
	double Weight_Dn2_3subLR = 1.;
	if(!fEventWeightSetToOne){
		Weight_Dn2=Dn2;
		Weight_Dn2Gap0 = Dn2Gap0;
		Weight_Dn2Gap2 = Dn2Gap2;
		Weight_Dn2Gap4 = Dn2Gap4;
		Weight_Dn2Gap6 = Dn2Gap6;
		Weight_Dn2Gap8 = Dn2Gap8;
		Weight_Dn2Gap10 = Dn2Gap10;
		Weight_Dn2Gap12 = Dn2Gap12;
		Weight_Dn2Gap14 = Dn2Gap14;
		Weight_Dn2_3subLM = Dn2_3subLM;
		Weight_Dn2_3subRM = Dn2_3subRM;
		Weight_Dn2_3subLR = Dn2_3subLR;
	}

	//calculate no eta-gap, gap1.0, gap1.4, 3subevent v2
	if(NtrksAfter > 1 && Dn2 != 0) {
		//..v2{2} = <cos2(phi1 - phi2)>
		TComplex v22 = correlator.Two(2, -2);
		double v22Re = v22.Re()/Dn2;
		//v22Re is the Real of v2{2}
		//Dn2 is the denominator of Generic framework
		//Profile histogram use Dn2 as weight
		profile.fChcn2[0]->Fill(Ntrks, v22Re, Weight_Dn2);

		//..v3{2} = <cos3(phi1 - phi2)>
		TComplex v32 = correlator.Two(3, -3);
		double v32Re = v32.Re()/Dn2;
		profile.fChcn2[1]->Fill(Ntrks, v32Re, Weight_Dn2);

		//..v4{2} = <cos4(phi1 - phi2)>
		TComplex v42 = correlator.Two(4, -4);
		double v42Re = v42.Re()/Dn2;
		profile.fChcn2[2]->Fill(Ntrks, v42Re, Weight_Dn2);

		
		//v5{2}=<cos4(phi1 - phi2)>
		TComplex v52 = correlator.Two(5, -5);
		double v52Re = v52.Re()/Dn2;
		profile.fChcn2[3]->Fill(Ntrks, v52Re, Weight_Dn2);

		//v6{2}=<cos4(phi1 - phi2)>
		TComplex v62 = correlator.Two(6, -6);
		double v62Re = v62.Re()/Dn2;
		profile.fChcn2[4]->Fill(Ntrks, v62Re, Weight_Dn2);
		

	}

	if(NtrksAfterGap0M > 0 && NtrksAfterGap0P > 0 && Dn2Gap0 != 0)
	{
		//..v2{2} with eta Gap > 0
		TComplex v22Gap0 = correlator.TwoGap0(2, -2);
		double v22ReGap0 = v22Gap0.Re()/Dn2Gap0;
		profile.fChcn2_Gap0[0]->Fill(Ntrks, v22ReGap0, Weight_Dn2Gap0);

		//..v3{2} with eta Gap > 0
		TComplex v32Gap0 = correlator.TwoGap0(3, -3);
		double v32ReGap0 = v32Gap0.Re()/Dn2Gap0;
		profile.fChcn2_Gap0[1]->Fill(Ntrks, v32ReGap0, Weight_Dn2Gap0);

		//..v4{2} with eta Gap > 0
		TComplex v42Gap0 = correlator.TwoGap0(4, -4);
		double v42ReGap0 = v42Gap0.Re()/Dn2Gap0;
		profile.fChcn2_Gap0[2]->Fill(Ntrks, v42ReGap0, Weight_Dn2Gap0);

		//..v5{2} with eta Gap > 0
		TComplex v52Gap0 = correlator.TwoGap0(5, -5);
		double v52ReGap0 = v52Gap0.Re()/Dn2Gap0;
		profile.fChcn2_Gap0[3]->Fill(Ntrks, v52ReGap0, Weight_Dn2Gap0);

		//..v6{2} with eta Gap > 0
		TComplex v62Gap0 = correlator.TwoGap0(6, -6);
		double v62ReGap0 = v62Gap0.Re()/Dn2Gap0;
		profile.fChcn2_Gap0[4]->Fill(Ntrks, v62ReGap0, Weight_Dn2Gap0);
	}

	if(NtrksAfterGap2M > 0 && NtrksAfterGap2P > 0 && Dn2Gap2 != 0)
	{
		//..v2{2} with eta Gap > 0.2
		TComplex v22Gap2 = correlator.TwoGap2(2, -2);
		double v22ReGap2 = v22Gap2.Re()/Dn2Gap2;
		profile.fChcn2_Gap2[0]->Fill(Ntrks, v22ReGap2, Weight_Dn2Gap2);

		//..v3{2} with eta Gap > 0.2
		TComplex v32Gap2 = correlator.TwoGap2(3, -3);
		double v32ReGap2 = v32Gap2.Re()/Dn2Gap2;
		profile.fChcn2_Gap2[1]->Fill(Ntrks, v32ReGap2, Weight_Dn2Gap2);

		//..v4{2} with eta Gap > 0.2
		TComplex v42Gap2 = correlator.TwoGap2(4, -4);
		double v42ReGap2 = v42Gap2.Re()/Dn2Gap2;
		profile.fChcn2_Gap2[2]->Fill(Ntrks, v42ReGap2, Weight_Dn2Gap2);

		//..v5{2} with eta Gap > 0.2
		TComplex v52Gap2 = correlator.TwoGap2(5, -5);
		double v52ReGap2 = v52Gap2.Re()/Dn2Gap2;
		profile.fChcn2_Gap2[3]->Fill(Ntrks, v52ReGap2, Weight_Dn2Gap2);

		//..v6{2} with eta Gap > 0.2
		TComplex v62Gap2 = correlator.TwoGap2(6, -6);
		double v62ReGap2 = v62Gap2.Re()/Dn2Gap2;
		profile.fChcn2_Gap2[4]->Fill(Ntrks, v62ReGap2, Weight_Dn2Gap2);
	}

	if(NtrksAfterGap4M > 0 && NtrksAfterGap4P > 0 && Dn2Gap4 != 0)
	{
		//..v2{2} with eta Gap > 0.4
		TComplex v22Gap4 = correlator.TwoGap4(2, -2);
		double v22ReGap4 = v22Gap4.Re()/Dn2Gap4;
		profile.fChcn2_Gap4[0]->Fill(Ntrks, v22ReGap4, Weight_Dn2Gap4);

		//..v3{2} with eta Gap > 0.4
		TComplex v32Gap4 = correlator.TwoGap4(3, -3);
		double v32ReGap4 = v32Gap4.Re()/Dn2Gap4;
		profile.fChcn2_Gap4[1]->Fill(Ntrks, v32ReGap4, Weight_Dn2Gap4);

		//..v4{2} with eta Gap > 0.4
		TComplex v42Gap4 = correlator.TwoGap4(4, -4);
		double v42ReGap4 = v42Gap4.Re()/Dn2Gap4;
		profile.fChcn2_Gap4[2]->Fill(Ntrks, v42ReGap4, Weight_Dn2Gap4);

		//..v5{2} with eta Gap > 0.4
		TComplex v52Gap4 = correlator.TwoGap4(5, -5);
		double v52ReGap4 = v52Gap4.Re()/Dn2Gap4;
		profile.fChcn2_Gap4[3]->Fill(Ntrks, v52ReGap4, Weight_Dn2Gap4);

		//..v6{2} with eta Gap > 0.4
		TComplex v62Gap4 = correlator.TwoGap4(6, -6);
		double v62ReGap4 = v62Gap4.Re()/Dn2Gap4;
		profile.fChcn2_Gap4[4]->Fill(Ntrks, v62ReGap4, Weight_Dn2Gap4);
	}

	if(NtrksAfterGap6M > 0 && NtrksAfterGap6P > 0 && Dn2Gap6 != 0)
	{
		//..v2{2} with eta Gap > 0.6
		TComplex v22Gap6 = correlator.TwoGap6(2, -2);
		double v22ReGap6 = v22Gap6.Re()/Dn2Gap6;
		profile.fChcn2_Gap6[0]->Fill(Ntrks, v22ReGap6, Weight_Dn2Gap6);

		//..v3{2} with eta Gap > 0.6
		TComplex v32Gap6 = correlator.TwoGap6(3, -3);
		double v32ReGap6 = v32Gap6.Re()/Dn2Gap6;
		profile.fChcn2_Gap6[1]->Fill(Ntrks, v32ReGap6, Weight_Dn2Gap6);

		//..v4{2} with eta Gap > 0.6
		TComplex v42Gap6 = correlator.TwoGap6(4, -4);
		double v42ReGap6 = v42Gap6.Re()/Dn2Gap6;
		profile.fChcn2_Gap6[2]->Fill(Ntrks, v42ReGap6, Weight_Dn2Gap6);

		//..v5{2} with eta Gap > 0.6
		TComplex v52Gap6 = correlator.TwoGap6(5, -5);
		double v52ReGap6 = v52Gap6.Re()/Dn2Gap6;
		profile.fChcn2_Gap6[3]->Fill(Ntrks, v52ReGap6, Weight_Dn2Gap6);

		//..v6{2} with eta Gap > 0.6
		TComplex v62Gap6 = correlator.TwoGap6(6, -6);
		double v62ReGap6 = v62Gap6.Re()/Dn2Gap6;
		profile.fChcn2_Gap6[4]->Fill(Ntrks, v62ReGap6, Weight_Dn2Gap6);
	}

	if(NtrksAfterGap8M > 0 && NtrksAfterGap8P > 0 && Dn2Gap8 != 0)
	{
		//..v2{2} with eta Gap > 0.8
		TComplex v22Gap8 = correlator.TwoGap8(2, -2);
		double v22ReGap8 = v22Gap8.Re()/Dn2Gap8;
		profile.fChcn2_Gap8[0]->Fill(Ntrks, v22ReGap8, Weight_Dn2Gap8);

		//..v3{2} with eta Gap > 0.8
		TComplex v32Gap8 = correlator.TwoGap8(3, -3);
		double v32ReGap8 = v32Gap8.Re()/Dn2Gap8;
		profile.fChcn2_Gap8[1]->Fill(Ntrks, v32ReGap8, Weight_Dn2Gap8);

		//..v4{2} with eta Gap > 0.8
		TComplex v42Gap8 = correlator.TwoGap8(4, -4);
		double v42ReGap8 = v42Gap8.Re()/Dn2Gap8;
		profile.fChcn2_Gap8[2]->Fill(Ntrks, v42ReGap8, Weight_Dn2Gap8);

		//..v5{2} with eta Gap > 0.8
		TComplex v52Gap8 = correlator.TwoGap8(5, -5);
		double v52ReGap8 = v52Gap8.Re()/Dn2Gap8;
		profile.fChcn2_Gap8[3]->Fill(Ntrks, v52ReGap8, Weight_Dn2Gap8);

		//..v6{2} with eta Gap > 0.8
		TComplex v62Gap8 = correlator.TwoGap8(6, -6);
		double v62ReGap8 = v62Gap8.Re()/Dn2Gap8;
		profile.fChcn2_Gap8[4]->Fill(Ntrks, v62ReGap8, Weight_Dn2Gap8);
	}

	if(NtrksAfterGap10M > 0 && NtrksAfterGap10P > 0 && Dn2Gap10 != 0)
	{
		//..v2{2} with eta Gap > 1.0
		TComplex v22Gap10 = correlator.TwoGap10(2, -2);
		double v22ReGap10 = v22Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[0]->Fill(Ntrks, v22ReGap10, Weight_Dn2Gap10);

		//..v3{2} with eta Gap > 1.0
		TComplex v32Gap10 = correlator.TwoGap10(3, -3);
		double v32ReGap10 = v32Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[1]->Fill(Ntrks, v32ReGap10, Weight_Dn2Gap10);

		//..v4{2} with eta Gap > 1.0
		TComplex v42Gap10 = correlator.TwoGap10(4, -4);
		double v42ReGap10 = v42Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[2]->Fill(Ntrks, v42ReGap10, Weight_Dn2Gap10);

		//..v5{2} with eta Gap > 1.0
		TComplex v52Gap10 = correlator.TwoGap10(5, -5);
		double v52ReGap10 = v52Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[3]->Fill(Ntrks, v52ReGap10, Weight_Dn2Gap10);

		//..v6{2} with eta Gap > 1.0
		TComplex v62Gap10 = correlator.TwoGap10(6, -6);
		double v62ReGap10 = v62Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[4]->Fill(Ntrks, v62ReGap10, Weight_Dn2Gap10);
	}

	if(NtrksAfterGap12M > 0 && NtrksAfterGap12P > 0 && Dn2Gap12 != 0)
	{
		//..v2{2} with eta Gap > 1.2
		TComplex v22Gap12 = correlator.TwoGap12(2, -2);
		double v22ReGap12 = v22Gap12.Re()/Dn2Gap12;
		profile.fChcn2_Gap12[0]->Fill(Ntrks, v22ReGap12, Weight_Dn2Gap12);

		//..v3{2} with eta Gap > 1.2
		TComplex v32Gap12 = correlator.TwoGap12(3, -3);
		double v32ReGap12 = v32Gap12.Re()/Dn2Gap12;
		profile.fChcn2_Gap12[1]->Fill(Ntrks, v32ReGap12, Weight_Dn2Gap12);

		//..v4{2} with eta Gap > 1.2
		TComplex v42Gap12 = correlator.TwoGap12(4, -4);
		double v42ReGap12 = v42Gap12.Re()/Dn2Gap12;
		profile.fChcn2_Gap12[2]->Fill(Ntrks, v42ReGap12, Weight_Dn2Gap12);

		//..v5{2} with eta Gap > 1.2
		TComplex v52Gap12 = correlator.TwoGap12(5, -5);
		double v52ReGap12 = v52Gap12.Re()/Dn2Gap12;
		profile.fChcn2_Gap12[3]->Fill(Ntrks, v52ReGap12, Weight_Dn2Gap12);

		//..v6{2} with eta Gap > 1.2
		TComplex v62Gap12 = correlator.TwoGap12(6, -6);
		double v62ReGap12 = v62Gap12.Re()/Dn2Gap12;
		profile.fChcn2_Gap12[4]->Fill(Ntrks, v62ReGap12, Weight_Dn2Gap12);
	}

	if(NtrksAfterGap14M > 0 && NtrksAfterGap14P > 0 && Dn2Gap14 != 0)
	{
		//..v2{2} with eta Gap > 1.4
		TComplex v22Gap14 = correlator.TwoGap14(2, -2);
		double v22ReGap14 = v22Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[0]->Fill(Ntrks, v22ReGap14, Weight_Dn2Gap14);

		//..v3{2} with eta Gap > 1.4
		TComplex v32Gap14 = correlator.TwoGap14(3, -3);
		double v32ReGap14 = v32Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[1]->Fill(Ntrks, v32ReGap14, Weight_Dn2Gap14);

		//..v4{2} with eta Gap > 1.4
		TComplex v42Gap14 = correlator.TwoGap14(4, -4);
		double v42ReGap14 = v42Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[2]->Fill(Ntrks, v42ReGap14, Weight_Dn2Gap14);

		//..v5{2} with eta Gap > 1.4
		TComplex v52Gap14 = correlator.TwoGap14(5, -5);
		double v52ReGap14 = v52Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[3]->Fill(Ntrks, v52ReGap14, Weight_Dn2Gap14);

		//..v6{2} with eta Gap > 1.4
		TComplex v62Gap14 = correlator.TwoGap14(6, -6);
		double v62ReGap14 = v62Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[4]->Fill(Ntrks, v62ReGap14, Weight_Dn2Gap14);
	}

	//..for 3-subevent method, Gap0
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && Dn2_3subLM != 0)
	{//..left+middle
		TComplex v22_3subLM = correlator.Two_3SubLM(2, -2);
		double v22Re_3subLM = v22_3subLM.Re()/Dn2_3subLM;
		profile.fChcn2_3subLM[0]->Fill(Ntrks, v22Re_3subLM, Weight_Dn2_3subLM);

		TComplex v32_3subLM = correlator.Two_3SubLM(3, -3);
		double v32Re_3subLM = v32_3subLM.Re()/Dn2_3subLM;
		profile.fChcn2_3subLM[1]->Fill(Ntrks, v32Re_3subLM, Weight_Dn2_3subLM);

		TComplex v42_3subLM = correlator.Two_3SubLM(4, -4);
		double v42Re_3subLM = v42_3subLM.Re()/Dn2_3subLM;
		profile.fChcn2_3subLM[2]->Fill(Ntrks, v42Re_3subLM, Weight_Dn2_3subLM);
	}

	if(NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn2_3subRM != 0)
	{//..right+middle
		TComplex v22_3subRM = correlator.Two_3SubRM(2, -2);
		double v22Re_3subRM = v22_3subRM.Re()/Dn2_3subRM;
		profile.fChcn2_3subRM[0]->Fill(Ntrks, v22Re_3subRM, Weight_Dn2_3subRM);

		TComplex v32_3subRM = correlator.Two_3SubRM(3, -3);
		double v32Re_3subRM = v32_3subRM.Re()/Dn2_3subRM;
		profile.fChcn2_3subRM[1]->Fill(Ntrks, v32Re_3subRM, Weight_Dn2_3subRM);

		TComplex v42_3subRM = correlator.Two_3SubRM(4, -4);
		double v42Re_3subRM = v42_3subRM.Re()/Dn2_3subRM;
		profile.fChcn2_3subRM[2]->Fill(Ntrks, v42Re_3subRM, Weight_Dn2_3subRM);
	}

	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && Dn2_3subLR != 0)
	{//..right+middle
		TComplex v22_3subLR = correlator.Two_3SubLR(2, -2);
		double v22Re_3subLR = v22_3subLR.Re()/Dn2_3subLR;
		profile.fChcn2_3subLR[0]->Fill(Ntrks, v22Re_3subLR, Weight_Dn2_3subLR);

		TComplex v32_3subLR = correlator.Two_3SubLR(3, -3);
		double v32Re_3subLR = v32_3subLR.Re()/Dn2_3subLR;
		profile.fChcn2_3subLR[1]->Fill(Ntrks, v32Re_3subLR, Weight_Dn2_3subLR);

		TComplex v42_3subLR = correlator.Two_3SubLR(4, -4);
		double v42Re_3subLR = v42_3subLR.Re()/Dn2_3subLR;
		profile.fChcn2_3subLR[2]->Fill(Ntrks, v42Re_3subLR, Weight_Dn2_3subLR);
	}

	//..calculate 3-particle correlations
	//................................
	double Dn3 = correlator.Three(0, 0, 0).Re();
	double Dn3Gap0A = correlator.ThreeGap0A(0, 0, 0).Re();
	double Dn3Gap2A = correlator.ThreeGap2A(0, 0, 0).Re();
	double Dn3Gap4A = correlator.ThreeGap4A(0, 0, 0).Re();
	double Dn3Gap6A = correlator.ThreeGap6A(0, 0, 0).Re();
	double Dn3Gap8A = correlator.ThreeGap8A(0, 0, 0).Re();
	double Dn3Gap10A = correlator.ThreeGap10A(0, 0, 0).Re();
	double Dn3Gap12A = correlator.ThreeGap12A(0, 0, 0).Re();
	double Dn3Gap0B = correlator.ThreeGap0B(0, 0, 0).Re();
	double Dn3Gap2B = correlator.ThreeGap2B(0, 0, 0).Re();
	double Dn3Gap4B = correlator.ThreeGap4B(0, 0, 0).Re();
	double Dn3Gap6B = correlator.ThreeGap6B(0, 0, 0).Re();
	double Dn3Gap8B = correlator.ThreeGap8B(0, 0, 0).Re();
	double Dn3Gap10B = correlator.ThreeGap10B(0, 0, 0).Re();
	double Dn3Gap12B = correlator.ThreeGap12B(0, 0, 0).Re();
	double Weight_Dn3 = 1.;
	double Weight_Dn3Gap0A = 1.;
	double Weight_Dn3Gap2A = 1.;
	double Weight_Dn3Gap4A = 1.;
	double Weight_Dn3Gap6A = 1.;
	double Weight_Dn3Gap8A = 1.;
	double Weight_Dn3Gap10A = 1.;
	double Weight_Dn3Gap12A = 1.;
	double Weight_Dn3Gap0B = 1.;
	double Weight_Dn3Gap2B = 1.;
	double Weight_Dn3Gap4B = 1.;
	double Weight_Dn3Gap6B = 1.;
	double Weight_Dn3Gap8B = 1.;
	double Weight_Dn3Gap10B = 1.;
	double Weight_Dn3Gap12B = 1.;
	if(!fEventWeightSetToOne){
		Weight_Dn3 = Dn3;
		Weight_Dn3Gap0A = Dn3Gap0A;
		Weight_Dn3Gap2A = Dn3Gap2A;
		Weight_Dn3Gap4A = Dn3Gap4A;
		Weight_Dn3Gap6A = Dn3Gap6A;
		Weight_Dn3Gap8A = Dn3Gap8A;
		Weight_Dn3Gap10A = Dn3Gap10A;
		Weight_Dn3Gap12A = Dn3Gap12A;
		Weight_Dn3Gap0B = Dn3Gap0B;
		Weight_Dn3Gap2B = Dn3Gap2B;
		Weight_Dn3Gap4B = Dn3Gap4B;
		Weight_Dn3Gap6B = Dn3Gap6B;
		Weight_Dn3Gap8B = Dn3Gap8B;
		Weight_Dn3Gap10B = Dn3Gap10B;
		Weight_Dn3Gap12B = Dn3Gap12B;
	}


	if(NtrksAfter > 2 && Dn3 != 0 )
	{
		//..v4{psi2}
		TComplex v422 = correlator.Three(4, -2, -2);
		double v422Re = v422.Re()/Dn3;
		profile.fChc422->Fill(Ntrks, v422Re, Weight_Dn3);

		//..v5{psi32}
		TComplex v532 = correlator.Three(5, -3, -2);
		double v532Re = v532.Re()/Dn3;
		profile.fChc532->Fill(Ntrks, v532Re, Weight_Dn3 );
	
	}
	// Gap 0
    // A-type
	if(NtrksAfterGap0M > 0 && NtrksAfterGap0P > 1 && Dn3Gap0A != 0)
	{

		TComplex v422Gap0A = correlator.ThreeGap0A(4, -2, -2);
		double v422Gap0ARe = v422Gap0A.Re()/Dn3Gap0A;
		profile.fChc422_Gap0A->Fill(Ntrks, v422Gap0ARe, Weight_Dn3Gap0A);

		TComplex v532Gap0A = correlator.ThreeGap0A(5, -3, -2);
		double v532Gap0ARe = v532Gap0A.Re()/Dn3Gap0A;
		profile.fChc532_Gap0A->Fill(Ntrks, v532Gap0ARe, Weight_Dn3Gap0A);

		// TComplex v633Gap0A = correlator.ThreeGap0A(6, -3, -3);
		// double v633Gap0ARe = v633Gap0A.Re()/Dn3Gap0A;
		// profile.fChsc633_Gap0A->Fill(Ntrks, v633Gap0ARe, Dn3Gap0A);
	}

	// B-type
	if(NtrksAfterGap0P > 0 && NtrksAfterGap0M > 1 && Dn3Gap0B != 0)
	{

		TComplex v422Gap0B = correlator.ThreeGap0B(4, -2, -2);
		double v422Gap0BRe = v422Gap0B.Re()/Dn3Gap0B;
		profile.fChc422_Gap0B->Fill(Ntrks, v422Gap0BRe, Weight_Dn3Gap0B);

		TComplex v532Gap0B = correlator.ThreeGap0B(5, -3, -2);
		double v532Gap0BRe = v532Gap0B.Re()/Dn3Gap0B;
		profile.fChc532_Gap0B->Fill(Ntrks, v532Gap0BRe, Weight_Dn3Gap0B);
	}

	// Gap 2
    // A-type
	if(NtrksAfterGap2M > 0 && NtrksAfterGap2P > 1 && Dn3Gap2A != 0)
	{

		TComplex v422Gap2A = correlator.ThreeGap2A(4, -2, -2);
		double v422Gap2ARe = v422Gap2A.Re()/Dn3Gap2A;
		profile.fChc422_Gap2A->Fill(Ntrks, v422Gap2ARe, Weight_Dn3Gap2A);

		TComplex v532Gap2A = correlator.ThreeGap2A(5, -3, -2);
		double v532Gap2ARe = v532Gap2A.Re()/Dn3Gap2A;
		profile.fChc532_Gap2A->Fill(Ntrks, v532Gap2ARe, Weight_Dn3Gap2A);
	}

	// B-type
	if(NtrksAfterGap2P > 0 && NtrksAfterGap2M > 1 && Dn3Gap2B != 0)
	{

		TComplex v422Gap2B = correlator.ThreeGap2B(4, -2, -2);
		double v422Gap2BRe = v422Gap2B.Re()/Dn3Gap2B;
		profile.fChc422_Gap2B->Fill(Ntrks, v422Gap2BRe, Weight_Dn3Gap2B);

		TComplex v532Gap2B = correlator.ThreeGap2B(5, -3, -2);
		double v532Gap2BRe = v532Gap2B.Re()/Dn3Gap2B;
		profile.fChc532_Gap2B->Fill(Ntrks, v532Gap2BRe, Weight_Dn3Gap2B);
	}

	// Gap 4
    // A-type
	if(NtrksAfterGap4M > 0 && NtrksAfterGap4P > 1 && Dn3Gap4A != 0)
	{

		TComplex v422Gap4A = correlator.ThreeGap4A(4, -2, -2);
		double v422Gap4ARe = v422Gap4A.Re()/Dn3Gap4A;
		profile.fChc422_Gap4A->Fill(Ntrks, v422Gap4ARe, Weight_Dn3Gap4A);

		TComplex v532Gap4A = correlator.ThreeGap4A(5, -3, -2);
		double v532Gap4ARe = v532Gap4A.Re()/Dn3Gap4A;
		profile.fChc532_Gap4A->Fill(Ntrks, v532Gap4ARe, Weight_Dn3Gap4A);
	}

	// B-type
	if(NtrksAfterGap4P > 0 && NtrksAfterGap4M > 1 && Dn3Gap4B != 0)
	{

		TComplex v422Gap4B = correlator.ThreeGap4B(4, -2, -2);
		double v422Gap4BRe = v422Gap4B.Re()/Dn3Gap4B;
		profile.fChc422_Gap4B->Fill(Ntrks, v422Gap4BRe, Weight_Dn3Gap4B);

		TComplex v532Gap4B = correlator.ThreeGap4B(5, -3, -2);
		double v532Gap4BRe = v532Gap4B.Re()/Dn3Gap4B;
		profile.fChc532_Gap4B->Fill(Ntrks, v532Gap4BRe, Weight_Dn3Gap4B);
	}

	// Gap 6
    // A-type
	if(NtrksAfterGap6M > 0 && NtrksAfterGap6P > 1 && Dn3Gap6A != 0)
	{

		TComplex v422Gap6A = correlator.ThreeGap6A(4, -2, -2);
		double v422Gap6ARe = v422Gap6A.Re()/Dn3Gap6A;
		profile.fChc422_Gap6A->Fill(Ntrks, v422Gap6ARe, Weight_Dn3Gap6A);

		TComplex v532Gap6A = correlator.ThreeGap6A(5, -3, -2);
		double v532Gap6ARe = v532Gap6A.Re()/Dn3Gap6A;
		profile.fChc532_Gap6A->Fill(Ntrks, v532Gap6ARe, Weight_Dn3Gap6A);
	}

	// B-type
	if(NtrksAfterGap6P > 0 && NtrksAfterGap6M > 1 && Dn3Gap6B != 0)
	{

		TComplex v422Gap6B = correlator.ThreeGap6B(4, -2, -2);
		double v422Gap6BRe = v422Gap6B.Re()/Dn3Gap6B;
		profile.fChc422_Gap6B->Fill(Ntrks, v422Gap6BRe, Weight_Dn3Gap6B);

		TComplex v532Gap6B = correlator.ThreeGap6B(5, -3, -2);
		double v532Gap6BRe = v532Gap6B.Re()/Dn3Gap6B;
		profile.fChc532_Gap6B->Fill(Ntrks, v532Gap6BRe, Weight_Dn3Gap6B);
	}

	// Gap 8
    // A-type
	if(NtrksAfterGap8M > 0 && NtrksAfterGap8P > 1 && Dn3Gap8A != 0)
	{

		TComplex v422Gap8A = correlator.ThreeGap8A(4, -2, -2);
		double v422Gap8ARe = v422Gap8A.Re()/Dn3Gap8A;
		profile.fChc422_Gap8A->Fill(Ntrks, v422Gap8ARe, Weight_Dn3Gap8A);

		TComplex v532Gap8A = correlator.ThreeGap8A(5, -3, -2);
		double v532Gap8ARe = v532Gap8A.Re()/Dn3Gap8A;
		profile.fChc532_Gap8A->Fill(Ntrks, v532Gap8ARe, Weight_Dn3Gap8A);
	}

	// B-type
	if(NtrksAfterGap8P > 0 && NtrksAfterGap8M > 1 && Dn3Gap8B != 0)
	{

		TComplex v422Gap8B = correlator.ThreeGap8B(4, -2, -2);
		double v422Gap8BRe = v422Gap8B.Re()/Dn3Gap8B;
		profile.fChc422_Gap8B->Fill(Ntrks, v422Gap8BRe, Weight_Dn3Gap8B);

		TComplex v532Gap8B = correlator.ThreeGap8B(5, -3, -2);
		double v532Gap8BRe = v532Gap8B.Re()/Dn3Gap8B;
		profile.fChc532_Gap8B->Fill(Ntrks, v532Gap8BRe, Weight_Dn3Gap8B);
	}

	//Gap 10
        // A-type
	if(NtrksAfterGap10M > 0 && NtrksAfterGap10P > 1 && Dn3Gap10A != 0)
	{

		TComplex v422Gap10A = correlator.ThreeGap10A(4, -2, -2);
		double v422Gap10ARe = v422Gap10A.Re()/Dn3Gap10A;
		profile.fChc422_Gap10A->Fill(Ntrks, v422Gap10ARe, Weight_Dn3Gap10A);

		TComplex v532Gap10A = correlator.ThreeGap10A(5, -3, -2);
		double v532Gap10ARe = v532Gap10A.Re()/Dn3Gap10A;
		profile.fChc532_Gap10A->Fill(Ntrks, v532Gap10ARe, Weight_Dn3Gap10A);

		// TComplex v633Gap10A = correlator.ThreeGap10A(6, -3, -3);
		// double v633Gap10ARe = v633Gap10A.Re()/Dn3Gap10A;
		// profile.fChsc633_Gap10A->Fill(Ntrks, v633Gap10ARe, Dn3Gap10A);
	}

        // B-type
	if(NtrksAfterGap10P > 0 && NtrksAfterGap10M > 1 && Dn3Gap10B != 0)
	{

		TComplex v422Gap10B = correlator.ThreeGap10B(4, -2, -2);
		double v422Gap10BRe = v422Gap10B.Re()/Dn3Gap10B;
		profile.fChc422_Gap10B->Fill(Ntrks, v422Gap10BRe, Weight_Dn3Gap10B);

		TComplex v532Gap10B = correlator.ThreeGap10B(5, -3, -2);
		double v532Gap10BRe = v532Gap10B.Re()/Dn3Gap10B;
		profile.fChc532_Gap10B->Fill(Ntrks, v532Gap10BRe, Weight_Dn3Gap10B);
	}

	// Gap 12
    // A-type
	if(NtrksAfterGap12M > 0 && NtrksAfterGap12P > 1 && Dn3Gap12A != 0)
	{

		TComplex v422Gap12A = correlator.ThreeGap12A(4, -2, -2);
		double v422Gap12ARe = v422Gap12A.Re()/Dn3Gap12A;
		profile.fChc422_Gap12A->Fill(Ntrks, v422Gap12ARe, Weight_Dn3Gap12A);

		TComplex v532Gap12A = correlator.ThreeGap12A(5, -3, -2);
		double v532Gap12ARe = v532Gap12A.Re()/Dn3Gap12A;
		profile.fChc532_Gap12A->Fill(Ntrks, v532Gap12ARe, Weight_Dn3Gap12A);
	}

	// B-type
	if(NtrksAfterGap12P > 0 && NtrksAfterGap12M > 1 && Dn3Gap12B != 0)
	{

		TComplex v422Gap12B = correlator.ThreeGap12B(4, -2, -2);
		double v422Gap12BRe = v422Gap12B.Re()/Dn3Gap12B;
		profile.fChc422_Gap12B->Fill(Ntrks, v422Gap12BRe, Weight_Dn3Gap12B);

		TComplex v532Gap12B = correlator.ThreeGap12B(5, -3, -2);
		double v532Gap12BRe = v532Gap12B.Re()/Dn3Gap12B;
		profile.fChc532_Gap12B->Fill(Ntrks, v532Gap12BRe, Weight_Dn3Gap12B);
	}

	//..calculate 4-particle correlations
	//................................
	double Dn4 = correlator.Four(0, 0, 0, 0).Re();
	double Dn4Gap0 = correlator.FourGap0(0, 0, 0, 0).Re();
	double Dn4Gap2 = correlator.FourGap2(0, 0, 0, 0).Re();
	double Dn4Gap4 = correlator.FourGap4(0, 0, 0, 0).Re();
	double Dn4Gap6 = correlator.FourGap6(0, 0, 0, 0).Re();
	double Dn4Gap8 = correlator.FourGap8(0, 0, 0, 0).Re();
	double Dn4Gap10 = correlator.FourGap10(0, 0, 0, 0).Re();
	double Dn4Gap12 = correlator.FourGap12(0, 0, 0, 0).Re();
	double Dn4_3subMMLR = correlator.Four_3SubMMLR(0, 0, 0, 0).Re();
	double Dn4_3subLLMR = correlator.Four_3SubLLMR(0, 0, 0, 0).Re();
	double Dn4_3subRRML = correlator.Four_3SubRRML(0, 0, 0, 0).Re();
	double Weight_Dn4 = 1.;
	double Weight_Dn4Gap0 = 1.;
	double Weight_Dn4Gap2 = 1.;
	double Weight_Dn4Gap4 = 1.;
	double Weight_Dn4Gap6 = 1.;
	double Weight_Dn4Gap8 = 1.;
	double Weight_Dn4Gap10 = 1.;
	double Weight_Dn4Gap12 = 1.;
	double Weight_Dn4_3subMMLR = 1.;
	double Weight_Dn4_3subLLMR = 1.;
	double Weight_Dn4_3subRRML = 1.;
	if(!fEventWeightSetToOne){
		Weight_Dn4 = Dn4;
		Weight_Dn4Gap0 = Dn4Gap0;
		Weight_Dn4Gap2 = Dn4Gap2;
		Weight_Dn4Gap4 = Dn4Gap4;
		Weight_Dn4Gap6 = Dn4Gap6;
		Weight_Dn4Gap8 = Dn4Gap8;
		Weight_Dn4Gap10 = Dn4Gap10;
		Weight_Dn4Gap12 = Dn4Gap12;
		Weight_Dn4_3subMMLR = Dn4_3subMMLR;
		Weight_Dn4_3subLLMR = Dn4_3subLLMR;
		Weight_Dn4_3subRRML = Dn4_3subRRML;
	}

	if(NtrksAfter > 3 && Dn4 != 0)
	{

		TComplex v24 = correlator.Four(2, 2, -2, -2);
		double v24Re = v24.Re()/Dn4;
		profile.fChcn4[0]->Fill(Ntrks, v24Re, Weight_Dn4);
		// fcn4Ntrks1bin[0][fBin]->Fill(NtrksAfter, v24Re, Dn4);

		TComplex v34 = correlator.Four(3, 3, -3, -3);
		double v34Re = v34.Re()/Dn4;
		profile.fChcn4[1]->Fill(Ntrks, v34Re, Weight_Dn4);
		// fcn4Ntrks1bin[1][fBin]->Fill(NtrksAfter, v34Re, Dn4);

		TComplex v44 = correlator.Four(4, 4, -4, -4);
		double v44Re = v44.Re()/Dn4;
		profile.fChcn4[2]->Fill(Ntrks, v44Re, Weight_Dn4);
		// fcn4Ntrks1bin[2][fBin]->Fill(NtrksAfter, v44Re, Dn4);

		//..SC(3,2,-3,-2)
		TComplex sc3232 = correlator.Four(3, 2, -3, -2);
		double sc3232Re = sc3232.Re()/Dn4;
		profile.fChsc3232->Fill(Ntrks, sc3232Re, Weight_Dn4);

		//..SC(4,2,-4,-2)
		TComplex sc4242 = correlator.Four(4, 2, -4, -2);
		double sc4242Re = sc4242.Re()/Dn4;
		profile.fChsc4242->Fill(Ntrks, sc4242Re, Weight_Dn4);

	}
	// Gap 0
	if(NtrksAfterGap0M > 1 && NtrksAfterGap0P > 1 && Dn4Gap0 !=0)
	{
		TComplex v24Gap0 = correlator.FourGap0(2, 2, -2, -2);
		double v24Gap0Re = v24Gap0.Re()/Dn4Gap0;
		profile.fChcn4_Gap0[0]->Fill(Ntrks, v24Gap0Re, Weight_Dn4Gap0);

		TComplex v34Gap0 = correlator.FourGap0(3, 3, -3, -3);
		double v34Gap0Re = v34Gap0.Re()/Dn4Gap0;
		profile.fChcn4_Gap0[1]->Fill(Ntrks, v34Gap0Re, Weight_Dn4Gap0);

		TComplex v44Gap0 = correlator.FourGap0(4, 4, -4, -4);
		double v44Gap0Re = v44Gap0.Re()/Dn4Gap0;
		profile.fChcn4_Gap0[2]->Fill(Ntrks, v44Gap0Re, Weight_Dn4Gap0);

		TComplex sc3232Gap0 = correlator.FourGap0(3, 2, -3, -2);
		double sc3232Gap0Re = sc3232Gap0.Re()/Dn4Gap0;
		profile.fChsc3232_Gap0->Fill(Ntrks, sc3232Gap0Re, Weight_Dn4Gap0);

		TComplex sc4242Gap0 = correlator.FourGap0(4, 2, -4, -2);
		double sc4242Gap0Re = sc4242Gap0.Re()/Dn4Gap0;
		profile.fChsc4242_Gap0->Fill(Ntrks, sc4242Gap0Re, Weight_Dn4Gap0);

		// TComplex sc6222Gap0 = correlator.FourGap0(6, -2, -2, -2);
		// double sc6222Gap0Re = sc6222Gap0.Re()/Dn4Gap0;
		// profile.fChsc6222_Gap0->Fill(Ntrks, sc6222Gap0Re, Dn4Gap0);
	}

	// Gap 2
	if(NtrksAfterGap2M > 1 && NtrksAfterGap2P > 1 && Dn4Gap2 !=0)
	{
		TComplex v24Gap2 = correlator.FourGap2(2, 2, -2, -2);
		double v24Gap2Re = v24Gap2.Re()/Dn4Gap2;
		profile.fChcn4_Gap2[0]->Fill(Ntrks, v24Gap2Re, Weight_Dn4Gap2);

		TComplex v34Gap2 = correlator.FourGap2(3, 3, -3, -3);
		double v34Gap2Re = v34Gap2.Re()/Dn4Gap2;
		profile.fChcn4_Gap2[1]->Fill(Ntrks, v34Gap2Re, Weight_Dn4Gap2);

		TComplex v44Gap2 = correlator.FourGap2(4, 4, -4, -4);
		double v44Gap2Re = v44Gap2.Re()/Dn4Gap2;
		profile.fChcn4_Gap2[2]->Fill(Ntrks, v44Gap2Re, Weight_Dn4Gap2);

		TComplex sc3232Gap2 = correlator.FourGap2(3, 2, -3, -2);
		double sc3232Gap2Re = sc3232Gap2.Re()/Dn4Gap2;
		profile.fChsc3232_Gap2->Fill(Ntrks, sc3232Gap2Re, Weight_Dn4Gap2);

		TComplex sc4242Gap2 = correlator.FourGap2(4, 2, -4, -2);
		double sc4242Gap2Re = sc4242Gap2.Re()/Dn4Gap2;
		profile.fChsc4242_Gap2->Fill(Ntrks, sc4242Gap2Re, Weight_Dn4Gap2);
	}

	// Gap 4
	if(NtrksAfterGap4M > 1 && NtrksAfterGap4P > 1 && Dn4Gap4 !=0)
	{
		TComplex v24Gap4 = correlator.FourGap4(2, 2, -2, -2);
		double v24Gap4Re = v24Gap4.Re()/Dn4Gap4;
		profile.fChcn4_Gap4[0]->Fill(Ntrks, v24Gap4Re, Weight_Dn4Gap4);

		TComplex v34Gap4 = correlator.FourGap4(3, 3, -3, -3);
		double v34Gap4Re = v34Gap4.Re()/Dn4Gap4;
		profile.fChcn4_Gap4[1]->Fill(Ntrks, v34Gap4Re, Weight_Dn4Gap4);

		TComplex v44Gap4 = correlator.FourGap4(4, 4, -4, -4);
		double v44Gap4Re = v44Gap4.Re()/Dn4Gap4;
		profile.fChcn4_Gap4[2]->Fill(Ntrks, v44Gap4Re, Weight_Dn4Gap4);

		TComplex sc3232Gap4 = correlator.FourGap4(3, 2, -3, -2);
		double sc3232Gap4Re = sc3232Gap4.Re()/Dn4Gap4;
		profile.fChsc3232_Gap4->Fill(Ntrks, sc3232Gap4Re, Weight_Dn4Gap4);

		TComplex sc4242Gap4 = correlator.FourGap4(4, 2, -4, -2);
		double sc4242Gap4Re = sc4242Gap4.Re()/Dn4Gap4;
		profile.fChsc4242_Gap4->Fill(Ntrks, sc4242Gap4Re, Weight_Dn4Gap4);
	}

	// Gap 6
	if(NtrksAfterGap6M > 1 && NtrksAfterGap6P > 1 && Dn4Gap6 !=0)
	{
		TComplex v24Gap6 = correlator.FourGap6(2, 2, -2, -2);
		double v24Gap6Re = v24Gap6.Re()/Dn4Gap6;
		profile.fChcn4_Gap6[0]->Fill(Ntrks, v24Gap6Re, Weight_Dn4Gap6);

		TComplex v34Gap6 = correlator.FourGap6(3, 3, -3, -3);
		double v34Gap6Re = v34Gap6.Re()/Dn4Gap6;
		profile.fChcn4_Gap6[1]->Fill(Ntrks, v34Gap6Re, Weight_Dn4Gap6);

		TComplex v44Gap6 = correlator.FourGap6(4, 4, -4, -4);
		double v44Gap6Re = v44Gap6.Re()/Dn4Gap6;
		profile.fChcn4_Gap6[2]->Fill(Ntrks, v44Gap6Re, Weight_Dn4Gap6);

		TComplex sc3232Gap6 = correlator.FourGap6(3, 2, -3, -2);
		double sc3232Gap6Re = sc3232Gap6.Re()/Dn4Gap6;
		profile.fChsc3232_Gap6->Fill(Ntrks, sc3232Gap6Re, Weight_Dn4Gap6);

		TComplex sc4242Gap6 = correlator.FourGap6(4, 2, -4, -2);
		double sc4242Gap6Re = sc4242Gap6.Re()/Dn4Gap6;
		profile.fChsc4242_Gap6->Fill(Ntrks, sc4242Gap6Re, Weight_Dn4Gap6);
	}

	// Gap 8
	if(NtrksAfterGap8M > 1 && NtrksAfterGap8P > 1 && Dn4Gap8 !=0)
	{
		TComplex v24Gap8 = correlator.FourGap8(2, 2, -2, -2);
		double v24Gap8Re = v24Gap8.Re()/Dn4Gap8;
		profile.fChcn4_Gap8[0]->Fill(Ntrks, v24Gap8Re, Weight_Dn4Gap8);

		TComplex v34Gap8 = correlator.FourGap8(3, 3, -3, -3);
		double v34Gap8Re = v34Gap8.Re()/Dn4Gap8;
		profile.fChcn4_Gap8[1]->Fill(Ntrks, v34Gap8Re, Weight_Dn4Gap8);

		TComplex v44Gap8 = correlator.FourGap8(4, 4, -4, -4);
		double v44Gap8Re = v44Gap8.Re()/Dn4Gap8;
		profile.fChcn4_Gap8[2]->Fill(Ntrks, v44Gap8Re, Weight_Dn4Gap8);

		TComplex sc3232Gap8 = correlator.FourGap8(3, 2, -3, -2);
		double sc3232Gap8Re = sc3232Gap8.Re()/Dn4Gap8;
		profile.fChsc3232_Gap8->Fill(Ntrks, sc3232Gap8Re, Weight_Dn4Gap8);

		TComplex sc4242Gap8 = correlator.FourGap8(4, 2, -4, -2);
		double sc4242Gap8Re = sc4242Gap8.Re()/Dn4Gap8;
		profile.fChsc4242_Gap8->Fill(Ntrks, sc4242Gap8Re, Weight_Dn4Gap8);
	}

	// Gap 10
	if(NtrksAfterGap10M > 1 && NtrksAfterGap10P > 1 && Dn4Gap10 !=0)
	{
		TComplex v24Gap10 = correlator.FourGap10(2, 2, -2, -2);
		double v24Gap10Re = v24Gap10.Re()/Dn4Gap10;
		profile.fChcn4_Gap10[0]->Fill(Ntrks, v24Gap10Re, Weight_Dn4Gap10);

		TComplex v34Gap10 = correlator.FourGap10(3, 3, -3, -3);
		double v34Gap10Re = v34Gap10.Re()/Dn4Gap10;
		profile.fChcn4_Gap10[1]->Fill(Ntrks, v34Gap10Re, Weight_Dn4Gap10);

		TComplex v44Gap10 = correlator.FourGap10(4, 4, -4, -4);
		double v44Gap10Re = v44Gap10.Re()/Dn4Gap10;
		profile.fChcn4_Gap10[2]->Fill(Ntrks, v44Gap10Re, Weight_Dn4Gap10);

		TComplex sc3232Gap10 = correlator.FourGap10(3, 2, -3, -2);
		double sc3232Gap10Re = sc3232Gap10.Re()/Dn4Gap10;
		profile.fChsc3232_Gap10->Fill(Ntrks, sc3232Gap10Re, Weight_Dn4Gap10);

		TComplex sc4242Gap10 = correlator.FourGap10(4, 2, -4, -2);
		double sc4242Gap10Re = sc4242Gap10.Re()/Dn4Gap10;
		profile.fChsc4242_Gap10->Fill(Ntrks, sc4242Gap10Re, Weight_Dn4Gap10);

		// TComplex sc6222Gap10 = correlator.FourGap10(6, -2, -2, -2);
		// double sc6222Gap10Re = sc6222Gap10.Re()/Dn4Gap10;
		// profile.fChsc6222_Gap10->Fill(Ntrks, sc6222Gap10Re, Dn4Gap10);
	}

	// Gap 12
	if(NtrksAfterGap12M > 1 && NtrksAfterGap12P > 1 && Dn4Gap12 !=0)
	{
		TComplex v24Gap12 = correlator.FourGap12(2, 2, -2, -2);
		double v24Gap12Re = v24Gap12.Re()/Dn4Gap12;
		profile.fChcn4_Gap12[0]->Fill(Ntrks, v24Gap12Re, Weight_Dn4Gap12);

		TComplex v34Gap12 = correlator.FourGap12(3, 3, -3, -3);
		double v34Gap12Re = v34Gap12.Re()/Dn4Gap12;
		profile.fChcn4_Gap12[1]->Fill(Ntrks, v34Gap12Re, Weight_Dn4Gap12);

		TComplex v44Gap12 = correlator.FourGap12(4, 4, -4, -4);
		double v44Gap12Re = v44Gap12.Re()/Dn4Gap12;
		profile.fChcn4_Gap12[2]->Fill(Ntrks, v44Gap12Re, Weight_Dn4Gap12);

		TComplex sc3232Gap12 = correlator.FourGap12(3, 2, -3, -2);
		double sc3232Gap12Re = sc3232Gap12.Re()/Dn4Gap12;
		profile.fChsc3232_Gap12->Fill(Ntrks, sc3232Gap12Re, Weight_Dn4Gap12);

		TComplex sc4242Gap12 = correlator.FourGap12(4, 2, -4, -2);
		double sc4242Gap12Re = sc4242Gap12.Re()/Dn4Gap12;
		profile.fChsc4242_Gap12->Fill(Ntrks, sc4242Gap12Re, Weight_Dn4Gap12);
	}


	//..3-subevent method
	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && NtrksAfter3subM > 1 && Dn4_3subMMLR != 0)
	{
		TComplex v24_3sub = correlator.Four_3SubMMLR(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3subMMLR;
		profile.fChcn4_3subMMLR[0]->Fill(Ntrks, v24_3subRe, Weight_Dn4_3subMMLR);

		TComplex v34_3sub = correlator.Four_3SubMMLR(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3subMMLR;
		profile.fChcn4_3subMMLR[1]->Fill(Ntrks, v34_3subRe, Weight_Dn4_3subMMLR);

		TComplex v44_3sub = correlator.Four_3SubMMLR(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3subMMLR;
		profile.fChcn4_3subMMLR[2]->Fill(Ntrks, v44_3subRe, Weight_Dn4_3subMMLR);

		TComplex sc3232_3subA = correlator.Four_3SubMMLR(3, 2, -3, -2);
		double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subMMLR;
		profile.fChsc3232_3subMMLRA->Fill(Ntrks, sc3232_3subARe, Weight_Dn4_3subMMLR);

		TComplex sc3232_3subB = correlator.Four_3SubMMLR(3, 2, -2, -3);
		double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subMMLR;
		profile.fChsc3232_3subMMLRB->Fill(Ntrks, sc3232_3subBRe, Weight_Dn4_3subMMLR);

		TComplex sc4242_3subA = correlator.Four_3SubMMLR(4, 2, -4, -2);
		double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subMMLR;
		profile.fChsc4242_3subMMLRA->Fill(Ntrks, sc4242_3subARe, Weight_Dn4_3subMMLR);

		TComplex sc4242_3subB = correlator.Four_3SubMMLR(4, 2, -2, -4);
		double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subMMLR;
		profile.fChsc4242_3subMMLRB->Fill(Ntrks, sc4242_3subBRe, Weight_Dn4_3subMMLR);
	}

	//..3-subevent method
	if(NtrksAfter3subL > 1 && NtrksAfter3subR > 0 && NtrksAfter3subM > 0 && Dn4_3subLLMR != 0)
	{
		TComplex v24_3sub = correlator.Four_3SubLLMR(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3subLLMR;
		profile.fChcn4_3subLLMR[0]->Fill(Ntrks, v24_3subRe, Weight_Dn4_3subLLMR);

		TComplex v34_3sub = correlator.Four_3SubLLMR(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3subLLMR;
		profile.fChcn4_3subLLMR[1]->Fill(Ntrks, v34_3subRe, Weight_Dn4_3subLLMR);

		TComplex v44_3sub = correlator.Four_3SubLLMR(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3subLLMR;
		profile.fChcn4_3subLLMR[2]->Fill(Ntrks, v44_3subRe, Weight_Dn4_3subLLMR);

		TComplex sc3232_3subA = correlator.Four_3SubLLMR(3, 2, -3, -2);
		double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subLLMR;
		profile.fChsc3232_3subLLMRA->Fill(Ntrks, sc3232_3subARe, Weight_Dn4_3subLLMR);

		TComplex sc3232_3subB = correlator.Four_3SubLLMR(3, 2, -2, -3);
		double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subLLMR;
		profile.fChsc3232_3subLLMRB->Fill(Ntrks, sc3232_3subBRe, Weight_Dn4_3subLLMR);

		TComplex sc4242_3subA = correlator.Four_3SubLLMR(4, 2, -4, -2);
		double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subLLMR;
		profile.fChsc4242_3subLLMRA->Fill(Ntrks, sc4242_3subARe, Weight_Dn4_3subLLMR);

		TComplex sc4242_3subB = correlator.Four_3SubLLMR(4, 2, -2, -4);
		double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subLLMR;
		profile.fChsc4242_3subLLMRB->Fill(Ntrks, sc4242_3subBRe, Weight_Dn4_3subLLMR);
	}

	//..3-subevent method
	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 1 && NtrksAfter3subM > 0 && Dn4_3subRRML != 0)
	{
		TComplex v24_3sub = correlator.Four_3SubRRML(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3subRRML;
		profile.fChcn4_3subRRML[0]->Fill(Ntrks, v24_3subRe, Weight_Dn4_3subRRML);

		TComplex v34_3sub = correlator.Four_3SubRRML(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3subRRML;
		profile.fChcn4_3subRRML[1]->Fill(Ntrks, v34_3subRe, Weight_Dn4_3subRRML);

		TComplex v44_3sub = correlator.Four_3SubRRML(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3subRRML;
		profile.fChcn4_3subRRML[2]->Fill(Ntrks, v44_3subRe, Weight_Dn4_3subRRML);

		TComplex sc3232_3subA = correlator.Four_3SubRRML(3, 2, -3, -2);
		double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subRRML;
		profile.fChsc3232_3subRRMLA->Fill(Ntrks, sc3232_3subARe, Weight_Dn4_3subRRML);

		TComplex sc3232_3subB = correlator.Four_3SubRRML(3, 2, -2, -3);
		double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subRRML;
		profile.fChsc3232_3subRRMLB->Fill(Ntrks, sc3232_3subBRe, Weight_Dn4_3subRRML);

		TComplex sc4242_3subA = correlator.Four_3SubRRML(4, 2, -4, -2);
		double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subRRML;
		profile.fChsc4242_3subRRMLA->Fill(Ntrks, sc4242_3subARe, Weight_Dn4_3subRRML);

		TComplex sc4242_3subB = correlator.Four_3SubRRML(4, 2, -2, -4);
		double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subRRML;
		profile.fChsc4242_3subRRMLB->Fill(Ntrks, sc4242_3subBRe, Weight_Dn4_3subRRML);
	}

	// //..calculate 5-particle correlations
	// //................................
	// double Dn5 = correlator.Five(0, 0, 0, 0, 0).Re();
	// double Dn5Gap10A = correlator.FiveGap10A(0, 0, 0, 0, 0).Re();
	// double Dn5Gap10B = correlator.FiveGap10B(0, 0, 0, 0, 0).Re();

	// // A type
	// if(NtrksAfterGap10M > 1 && NtrksAfterGap10P > 2 && Dn5Gap10A != 0)
	// {

	// 	TComplex v5_A42222 = correlator.FiveGap10A(4, 2, -2, -2, -2);
	// 	double v5_A42222Re = v5_A42222.Re()/Dn5Gap10A;
	// 	profile.fChc5_A42222->Fill(Ntrks, v5_A42222Re, Dn5Gap10A);

	// 	TComplex v5_A52322 = correlator.FiveGap10A(5, 2, -3, -2, -2);
	// 	double v5_A52322Re = v5_A52322.Re()/Dn5Gap10A;
	// 	profile.fChc5_A52322->Fill(Ntrks, v5_A52322Re, Dn5Gap10A);
	// }

	//..calculate 6-particle correlations
	//................................
	double Dn6 = correlator.Six(0, 0, 0, 0, 0, 0).Re();
	double Dn6Gap0 = correlator.SixGap0(0, 0, 0, 0, 0, 0).Re();
	double Dn6Gap10 = correlator.SixGap10(0, 0, 0, 0, 0, 0).Re();
	double Weight_Dn6 = 1.;
	double Weight_Dn6Gap0 = 1.;
	double Weight_Dn6Gap10 = 1.;
	if(!fEventWeightSetToOne){
		Weight_Dn6 = Dn6;
		Weight_Dn6Gap0 = Dn6Gap0;
		Weight_Dn6Gap10 = Dn6Gap10;
	}

	if(NtrksAfter > 5 && Dn6 != 0)
	{
		// v2{6}
		TComplex v26 = correlator.Six(2, 2, 2, -2, -2, -2);
		double v26Re = v26.Re()/Dn6;
		profile.fChcn6[0]->Fill(Ntrks, v26Re, Weight_Dn6);
	}

	if(NtrksAfterGap0M > 2 && NtrksAfterGap0P > 2 && Dn6Gap0 !=0){
		TComplex v26Gap0 = correlator.SixGap0(2, 2, 2, -2, -2, -2);
		double v26ReGap0 = v26Gap0.Re()/Dn6Gap0;
		profile.fChcn6_Gap0[0]->Fill(Ntrks, v26ReGap0, Weight_Dn6Gap0);
	}

	if(NtrksAfterGap10M > 2 && NtrksAfterGap10P > 2 && Dn6Gap10 !=0)
	{
		//..v2{6} with eta Gap > 1.0
		TComplex v26Gap10 = correlator.SixGap10(2, 2, 2, -2, -2, -2);
		double v26ReGap10 = v26Gap10.Re()/Dn6Gap10;
		profile.fChcn6_Gap10[0]->Fill(Ntrks, v26ReGap10, Weight_Dn6Gap10);

		// TComplex v6_222222 = correlator.SixGap10(2, 2, 2, -2, -2, -2);
		// double v6_222222Re = v6_222222.Re()/Dn6Gap10;
		// profile.fChc6_222222->Fill(Ntrks, v6_222222Re, Dn6Gap10);

		// TComplex v6_322322 = correlator.SixGap10(3, 2, 2, -3, -2, -2);
		// double v6_322322Re = v6_322322.Re()/Dn6Gap10;
		// profile.fChc6_322322->Fill(Ntrks, v6_322322Re, Dn6Gap10);
	}

	// 8-Particles Correlation
	double Dn8 = correlator.Eight(0, 0, 0, 0, 0, 0, 0, 0).Re();
	// Gap0 isn't calculated in this version
	double Dn8Gap0 = correlator.EightGap0(0, 0, 0, 0, 0, 0, 0, 0).Re();
	double Weight_Dn8 = 1.;
	// Gap0 isn't calculated in this version
	double Weight_Dn8Gap0 = 1.;
	if(!fEventWeightSetToOne){
		Weight_Dn8=Dn8;
		Weight_Dn8Gap0=Dn8Gap0;
	}

	if(NtrksAfter > 7 && Dn8 != 0)
	{
		// v2{8}
		TComplex v28 = correlator.Eight(2, 2, 2, 2, -2, -2, -2, -2);
		double v28Re = v28.Re()/Dn8;
		profile.fChcn8[0]->Fill(Ntrks, v28Re, Weight_Dn8);
	}

	if(NtrksAfterGap0M > 3 && NtrksAfterGap0P > 3 && Dn8Gap0 !=0){
		//..v2{8} with eta Gap > 0.
		TComplex v28Gap0 = correlator.EightGap0(2, 2, 2, 2, -2, -2, -2, -2);
		double v28ReGap0 = v28Gap0.Re()/Dn8Gap0;
		profile.fChcn8_Gap0[0]->Fill(Ntrks, v28ReGap0, Weight_Dn8Gap0);
	}


}


//_____________________________________________________________________________
ClassImp(PhysicsProfileFlowOnTheFly);
PhysicsProfileFlowOnTheFly::PhysicsProfileFlowOnTheFly() :
		fChsc4242(nullptr),
		fChsc4242_Gap0(nullptr),
		fChsc4242_Gap2(nullptr),
		fChsc4242_Gap4(nullptr),
		fChsc4242_Gap6(nullptr),
		fChsc4242_Gap8(nullptr),      
		fChsc4242_Gap10(nullptr),    
		fChsc4242_Gap12(nullptr),    
		fChsc4242_3sub(nullptr),	
		fChsc4242_3subMMLRA(nullptr),
		fChsc4242_3subMMLRB(nullptr),							
		fChsc4242_3subLLMRA(nullptr),							
		fChsc4242_3subLLMRB(nullptr),							
		fChsc4242_3subRRMLA(nullptr),							
		fChsc4242_3subRRMLB(nullptr),							
		fChsc4224_3sub(nullptr),							
		fChsc4242_3subGap2(nullptr),					
		fChsc4224_3subGap2(nullptr),					
		fChsc3232(nullptr),									
		fChsc3232_Gap0(nullptr),							
		fChsc3232_Gap2(nullptr),							
		fChsc3232_Gap4(nullptr),                            
		fChsc3232_Gap6(nullptr),                            
		fChsc3232_Gap8(nullptr),                            
		fChsc3232_Gap10(nullptr),                            
		fChsc3232_Gap12(nullptr),                            
		fChsc3232_3sub(nullptr),							
		fChsc3232_3subMMLRA(nullptr),							
		fChsc3232_3subMMLRB(nullptr),							
		fChsc3232_3subLLMRA(nullptr),							
		fChsc3232_3subLLMRB(nullptr),							
		fChsc3232_3subRRMLA(nullptr),							
		fChsc3232_3subRRMLB(nullptr),							
		fChsc3223_3sub(nullptr),							
		fChsc3232_3subGap2(nullptr),					
		fChsc3223_3subGap2(nullptr),					
		fChc422(nullptr), 
		fChc532(nullptr),
		fChc422_Gap0A(nullptr),   
		fChc422_Gap0B(nullptr),   
		fChc532_Gap0A(nullptr),   
		fChc532_Gap0B(nullptr),   
		fChc422_Gap2A(nullptr),   
		fChc422_Gap2B(nullptr),   
		fChc532_Gap2A(nullptr),   
		fChc532_Gap2B(nullptr),   
		fChc422_Gap4A(nullptr),   
		fChc422_Gap4B(nullptr),   
		fChc532_Gap4A(nullptr),   
		fChc532_Gap4B(nullptr),   
		fChc422_Gap6A(nullptr),   
		fChc422_Gap6B(nullptr),   
		fChc532_Gap6A(nullptr),   
		fChc532_Gap6B(nullptr),   
		fChc422_Gap8A(nullptr),   
		fChc422_Gap8B(nullptr),   
		fChc532_Gap8A(nullptr),   
		fChc532_Gap8B(nullptr),   
		fChc422_Gap10A(nullptr), 
		fChc422_Gap10B(nullptr), 
		fChc532_Gap10A(nullptr), 
		fChc532_Gap10B(nullptr),
		fChc532_Gap12A(nullptr), 
		fChc532_Gap12B(nullptr)
		// fChc5_A42222(nullptr),
		// fChc5_A52322(nullptr),
		// fChc6_222222(nullptr),
		// fChc6_322322(nullptr),
		// fChsc6222_Gap0(nullptr),
		// fChsc6222_Gap10(nullptr),
		// fChsc633_Gap10A(nullptr)
{
		//Initial TProfile memory
		memset(fChcn2, 0, sizeof(fChcn2));
		memset(fChcn2_Gap0, 0, sizeof(fChcn2_Gap0));
		memset(fChcn2_Gap2, 0, sizeof(fChcn2_Gap2));
		memset(fChcn2_Gap4, 0, sizeof(fChcn2_Gap4));
		memset(fChcn2_Gap6, 0, sizeof(fChcn2_Gap6));
		memset(fChcn2_Gap8, 0, sizeof(fChcn2_Gap8));
		memset(fChcn2_Gap10, 0, sizeof(fChcn2_Gap10));
		memset(fChcn2_Gap12, 0, sizeof(fChcn2_Gap12));
		memset(fChcn2_Gap14, 0, sizeof(fChcn2_Gap14));
		memset(fChcn2_Gap16, 0, sizeof(fChcn2_Gap16));
		memset(fChcn2_Gap18, 0, sizeof(fChcn2_Gap18));

		memset(fChcn2_3subLM, 0, sizeof(fChcn2_3subLM));
		memset(fChcn2_3subRM, 0, sizeof(fChcn2_3subRM));
		memset(fChcn2_3subLR, 0, sizeof(fChcn2_3subLR));
		memset(fChcn2_3subGap2LM, 0, sizeof(fChcn2_3subGap2LM));
		memset(fChcn2_3subGap2RM, 0, sizeof(fChcn2_3subGap2RM));

		memset(fChcn4, 0, sizeof(fChcn4));
		memset(fChcn4_Gap0, 0, sizeof(fChcn4_Gap0));
		memset(fChcn4_Gap2, 0, sizeof(fChcn4_Gap2));
		memset(fChcn4_Gap4, 0, sizeof(fChcn4_Gap4));
		memset(fChcn4_Gap6, 0, sizeof(fChcn4_Gap6));
		memset(fChcn4_Gap8, 0, sizeof(fChcn4_Gap8));
		memset(fChcn4_Gap10, 0, sizeof(fChcn4_Gap10));
		memset(fChcn4_Gap12, 0, sizeof(fChcn4_Gap12));
		memset(fChcn4_3subMMLR, 0, sizeof(fChcn4_3subMMLR));
		memset(fChcn4_3subLLMR, 0, sizeof(fChcn4_3subLLMR));
		memset(fChcn4_3subRRML, 0, sizeof(fChcn4_3subRRML));
		memset(fChcn4_3subGap2, 0, sizeof(fChcn4_3subGap2));

		// memset(fChcn6, 0, sizeof(fChcn6));
		// memset(fChcn6_Gap0, 0, sizeof(fChcn6_Gap0));
		// memset(fChcn6_Gap10, 0, sizeof(fChcn6_Gap10));

		// memset(fChcn8, 0, sizeof(fChcn8));
		// memset(fChcn8_Gap0, 0, sizeof(fChcn8_Gap0));
}