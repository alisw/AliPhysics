#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TList.h"
#include "TH2F.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDVertex.h"

#include "AliProtonFeedDownAnalysis.h"
#include "AliProtonAnalysisBase.h"
#include "AliProtonFeedDownAnalysisTask.h"
#include <Riostream.h>

ClassImp(AliProtonFeedDownAnalysisTask)
  
//________________________________________________________________________ 
AliProtonFeedDownAnalysisTask::AliProtonFeedDownAnalysisTask()
  : AliAnalysisTask(), fESD(0), fAOD(0), fMC(0),
    fList(0), fProtonAnalysis(0),fStatHist(0)
 {
  //Dummy constructor
  
}

//________________________________________________________________________
AliProtonFeedDownAnalysisTask::AliProtonFeedDownAnalysisTask(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fAOD(0), fMC(0),
    fList(0), fProtonAnalysis(0),fStatHist(0)
    {
	// Constructor
	
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TList container
	DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliProtonFeedDownAnalysisTask::ConnectInputData(Option_t *) 
{
	// Connect ESD or AOD here
	// Called once
	TString gAnalysisLevel = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisLevel(); 
	
	TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
	if (!tree) 
	{
		Printf("ERROR: Could not read chain from input slot 0");
	} 
	else 
	{
		if(gAnalysisLevel == "ESD") 
		{
			AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());   
			if (!esdH) 
			{
				Printf("ERROR: Could not get ESDInputHandler");
			} 
			else
				fESD = esdH->GetEvent();
			AliMCEventHandler* mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
			if (!mcH) 
			{
				Printf("ERROR: Could not retrieve MC event handler");
			}
			else
				fMC = mcH->MCEvent();
		}
		else if(gAnalysisLevel == "AOD") 
		{
			AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());	
			if (!aodH) 
			{
				Printf("ERROR: Could not get AODInputHandler");
			}
			 else
				fAOD = aodH->GetEvent();
		}
		else
			Printf("Wrong analysis type: Only ESD, AOD  types are allowed!");
	}
}
//________________________________________________________________________
void AliProtonFeedDownAnalysisTask::CreateOutputObjects() 
{
	// Create output objects
	// Called once
	fList = new TList();
	fList->Add(fProtonAnalysis->GetProtonContainer());
	fList->Add(fProtonAnalysis->GetAntiProtonContainer());
	fList->Add(fProtonAnalysis->GetLambdaHist());
	fList->Add(fProtonAnalysis->GetLambdaweightedHist());
	fList->Add(fProtonAnalysis->GetAntiLambdaHist());
	fList->Add(fProtonAnalysis->GetAntiLambdaweightedHist());
	fStatHist=new TH1F("StatsHist","StatsHist",10,-0.5,9.5);
	fList->Add(fStatHist);
	fStatHist->GetXaxis()->SetBinLabel(1,"level1cut");
	fStatHist->GetXaxis()->SetBinLabel(2,"level2cut");
	fStatHist->GetXaxis()->SetBinLabel(3,"level3cut");
	fStatHist->GetXaxis()->SetBinLabel(4,"level4cut");
}
//________________________________________________________________________
void AliProtonFeedDownAnalysisTask::Exec(Option_t *) 
{
	// Main loop
	// Called for each event
	TString gAnalysisLevel = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisLevel(); 
	//TString gAnalysisLevel = (fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisLevel(); 
	if(gAnalysisLevel == "ESD") 
	{
		if (!fESD) 
		{
			Printf("ERROR: fESD not available");
			return;
		}
		fStatHist->Fill(0);
		if(dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->IsEventTriggered(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetTriggerMode())) 
		{
			fStatHist->Fill(1);
			const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
			if(vertex) 
			{
				fStatHist->Fill(2);
				Printf("Proton ESD analysis task: There are %d tracks in this event", fESD->GetNumberOfTracks());
				AliStack* stack1=0x0;
				if(!fMC)
					return;			
				stack1 = fMC->Stack();
				if(!stack1)
					return;
				fStatHist->Fill(3);	
				fProtonAnalysis->Analyze(fESD,vertex,stack1);
				fProtonAnalysis->Analyze(stack1);

			}//reconstructed vertex
		}//triggered event
	}//ESD analysis              
	
	else if(gAnalysisLevel == "AOD") 
	{
		if (!fAOD) 
		{
			Printf("ERROR: fAOD not available");
			return;
		}
		Printf("Proton AOD analysis task: There are %d tracks in this event", fAOD->GetNumberOfTracks());
		fProtonAnalysis->Analyze(fAOD);
	}//AOD analysis
	                 
	
	// Post output data.
	PostData(0, fList);
}    
//__________________________________________________________________________________________
void AliProtonFeedDownAnalysisTask::Terminate(Option_t *) 
{

}


