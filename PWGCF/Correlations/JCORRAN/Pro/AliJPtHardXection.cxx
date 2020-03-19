#include "TChain.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TArrayI.h" 
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TRandom3.h"

#include "AliLog.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliInputEventHandler.h"
#include "AliGenPythiaEventHeader.h"
#include "AliMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
#include "AliJPtHardXection.h"
#include "AliHeader.h" //KF//
#include <iostream>

using std::cout;
using std::endl;

ClassImp(AliJPtHardXection)

//Filip Krizek 1st March 2013

//---------------------------------------------------------------------
AliJPtHardXection::AliJPtHardXection() :
AliAnalysisTaskSE(),
fOutputList(NULL),
fh1Xsec(0x0),
fh1Trials(0x0),
fh1PtHard(0x0),
fh1PtHardNoW(0x0),  
fh1PtHardTrials(0x0),
fAvgTrials(1)
{
   // default Constructor
}

//---------------------------------------------------------------------

AliJPtHardXection::AliJPtHardXection(const char *name) :
AliAnalysisTaskSE(name),
fOutputList(NULL),
fh1Xsec(0x0),
fh1Trials(0x0), 
fh1PtHard(0x0),
fh1PtHardNoW(0x0),  
fh1PtHardTrials(0x0),
fAvgTrials(1)
{
// Constructor
   DefineOutput(1, TList::Class());
}

//--------------------------------------------------------------
AliJPtHardXection::AliJPtHardXection(const AliJPtHardXection& a):
AliAnalysisTaskSE(a.GetName()),
fOutputList(a.fOutputList),
fh1Xsec(a.fh1Xsec),
fh1Trials(a.fh1Trials),
fh1PtHard(a.fh1PtHard),
fh1PtHardNoW(a.fh1PtHardNoW),  
fh1PtHardTrials(a.fh1PtHardTrials),
fAvgTrials(a.fAvgTrials)
{
   //Copy Constructor
}
//--------------------------------------------------------------

AliJPtHardXection& AliJPtHardXection::operator = (const AliJPtHardXection& a){
  // assignment operator
  this->~AliJPtHardXection();
  new(this) AliJPtHardXection(a);
  return *this;
}
//--------------------------------------------------------------

AliJPtHardXection::~AliJPtHardXection()
{
   //Destructor 
   delete fOutputList; // ????
}

//--------------------------------------------------------------


Bool_t AliJPtHardXection::Notify()
{
	return kTRUE;
}
//--------------------------------------------------------------

void AliJPtHardXection::Init()
{

}

//--------------------------------------------------------------

void AliJPtHardXection::UserCreateOutputObjects()
{
	// Create histograms and initilize variables
	OpenFile(1);
	if(!fOutputList) fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);

	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);

	const Int_t nBinPt = 150;
	Double_t binLimitsPt[nBinPt+1];
	for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
		if(iPt == 0){
			binLimitsPt[iPt] = -50.0;
		}else{// 1.0
			binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 2.;
		}
	}

	fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
	fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
	fOutputList->Add(fh1Xsec);
	fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
	fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
	fOutputList->Add(fh1Trials);
	fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
	fOutputList->Add(fh1PtHard);
	fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
	fOutputList->Add(fh1PtHardNoW);
	fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);
	fOutputList->Add(fh1PtHardTrials);      


	// =========== Switch on Sumw2 for all histos ===========
	for(Int_t i=0; i<fOutputList->GetEntries(); i++){
		TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
		if(h1){
			h1->Sumw2();
			continue;
		}
		THnSparse *hn = dynamic_cast<THnSparse*>(fOutputList->At(i));
		if(hn){
			hn->Sumw2();
		}	  
	}
	TH1::AddDirectory(oldStatus);

	PostData(1, fOutputList);
}

//--------------------------------------------------------------------

void AliJPtHardXection::UserExec(Option_t *)
{
	//User Exec
	//Event loop
	Double_t eventW  = 1.0;
	Double_t ptHard  = 0.0;
	Double_t xsection  = 0.0;
	Double_t nTrials = 1.0; // Trials for MC trigger

	AliMCEvent *mcEvent = MCEvent();

	AliGenPythiaEventHeader *pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
	if(pythiaGenHeader){
		nTrials = pythiaGenHeader->Trials();
		xsection = pythiaGenHeader->GetXsection();
		ptHard  = pythiaGenHeader->GetPtHard();
		fh1Xsec->Fill("<#sigma>",xsection);
      		fh1Trials->Fill("#sum{ntrials}",nTrials);
		fh1PtHard->Fill(ptHard,eventW);
		fh1PtHardNoW->Fill(ptHard,1);
		fh1PtHardTrials->Fill(ptHard,nTrials);
	}

	PostData(1, fOutputList);
}

//----------------------------------------------------------------------------
void AliJPtHardXection::Terminate(const Option_t *)
{
	// Draw result to the screen
	// Called once at the end of the query

	if(fDebug) printf("AliJPtHardXection DONE\n");
	if(!GetOutputData(1)) return;
}
