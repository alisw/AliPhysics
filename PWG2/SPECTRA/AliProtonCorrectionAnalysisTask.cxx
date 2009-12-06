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
#include "AliProtonCorrectionAnalysisTask.h"

#include "AliProtonAbsorptionCorrection.h"
#include "AliProtonSpectraCorrection.h"

#include <Riostream.h>

ClassImp(AliProtonCorrectionAnalysisTask)
  
//________________________________________________________________________ 
AliProtonCorrectionAnalysisTask::AliProtonCorrectionAnalysisTask()
  : AliAnalysisTask(), fESD(0), fAOD(0), fMC(0),
    fList(0),fProtonAnalysisBase(0),fProtonAbsorptionCorrection(0), fProtonFeedDownAnalysis(0),fProtonSpectraCorrection(0),fStatHist(0),
    fIsOn_AliProtonAbsorptionCorrection(0),fIsOn_AliProtonFeedDownAnalysis(0),fIsOn_AliProtonSpectraCorrection(0)
 {
  //Dummy constructor
  
}

//________________________________________________________________________
AliProtonCorrectionAnalysisTask::AliProtonCorrectionAnalysisTask(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fAOD(0), fMC(0),
    fList(0),fProtonAnalysisBase(0),fProtonAbsorptionCorrection(0), fProtonFeedDownAnalysis(0),fProtonSpectraCorrection(0),fStatHist(0),
    fIsOn_AliProtonAbsorptionCorrection(0),fIsOn_AliProtonFeedDownAnalysis(0),fIsOn_AliProtonSpectraCorrection(0)
    {
	// Constructor
	
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TList container
	DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliProtonCorrectionAnalysisTask::ConnectInputData(Option_t *) 
{
	// Connect ESD or AOD here
	// Called once
	TString gAnalysisLevel = fProtonAnalysisBase->GetAnalysisLevel(); 
	
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
void AliProtonCorrectionAnalysisTask::CreateOutputObjects() 
{
	// Create output objects
	// Called once
	fList = new TList();
	if(fIsOn_AliProtonAbsorptionCorrection)
	{
		fProtonAbsorptionCorrection->GetProtonContainer()->SetName("containerProtonsAbsorptionCorrection");
		fList->Add(fProtonAbsorptionCorrection->GetProtonContainer());
		fProtonAbsorptionCorrection->GetAntiProtonContainer()->SetName("containerAntiProtonsAbsorptionCorrection");
		fList->Add(fProtonAbsorptionCorrection->GetAntiProtonContainer());
	}
	if(fIsOn_AliProtonFeedDownAnalysis)
	{
		fProtonFeedDownAnalysis->GetProtonContainer()->SetName("containerProtonsFeedDown");
		fList->Add(fProtonFeedDownAnalysis->GetProtonContainer());
		fProtonFeedDownAnalysis->GetAntiProtonContainer()->SetName("containerAntiProtonsFeedDown");
		fList->Add(fProtonFeedDownAnalysis->GetAntiProtonContainer());
		fList->Add(fProtonFeedDownAnalysis->GetLambdaHist());
		fList->Add(fProtonFeedDownAnalysis->GetLambdaweightedHist());
		fList->Add(fProtonFeedDownAnalysis->GetAntiLambdaHist());
		fList->Add(fProtonFeedDownAnalysis->GetAntiLambdaweightedHist());
	}
	if(fIsOn_AliProtonSpectraCorrection)
	{
		fProtonSpectraCorrection->GetProtonContainer()->SetName("containerProtonsSpectraCorrection");
		fList->Add(fProtonSpectraCorrection->GetProtonContainer());
		fProtonSpectraCorrection->GetProtonContainer()->SetName("containerAntiProtonsSpectraCorrection");
		fList->Add(fProtonSpectraCorrection->GetAntiProtonContainer());
	}
	fStatHist=new TH1F("StatsHist","StatsHist",10,-0.5,9.5);
	fList->Add(fStatHist);
	fStatHist->GetXaxis()->SetBinLabel(1,"level1cutESD");
	fStatHist->GetXaxis()->SetBinLabel(2,"level2cutTrigger");
	fStatHist->GetXaxis()->SetBinLabel(3,"level3cutVerstex");
	fStatHist->GetXaxis()->SetBinLabel(4,"level4cutMC");
}
//________________________________________________________________________
void AliProtonCorrectionAnalysisTask::Exec(Option_t *) 
{
	// Main loop
	// Called for each event
	TString gAnalysisLevel =fProtonAnalysisBase->GetAnalysisLevel(); 
	//TString gAnalysisLevel = (fProtonFeedDownAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisLevel(); 
	if(gAnalysisLevel == "ESD") 
	{
		if (!fESD) 
		{
			Printf("ERROR: fESD not available");
			return;
		}
		fStatHist->Fill(0);
		if(fProtonAnalysisBase->IsEventTriggered(fESD,fProtonAnalysisBase->GetTriggerMode())) 
		{
			fStatHist->Fill(1);
			const AliESDVertex *vertex = fProtonAnalysisBase->GetVertex(fESD,fProtonAnalysisBase->GetAnalysisMode(),fProtonAnalysisBase->GetVxMax(),fProtonAnalysisBase->GetVyMax(),fProtonAnalysisBase->GetVzMax());
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
				if(fIsOn_AliProtonAbsorptionCorrection)
					fProtonAbsorptionCorrection->FillAbsorptionMaps(fESD,vertex,fMC);
				if(fIsOn_AliProtonFeedDownAnalysis)
				{	
					fProtonFeedDownAnalysis->Analyze(fESD,vertex,stack1);
					fProtonFeedDownAnalysis->Analyze(stack1);
				}	
				if(fIsOn_AliProtonSpectraCorrection)
					fProtonSpectraCorrection->FillCorrectionMaps(fESD,vertex,fMC);

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
		if(fIsOn_AliProtonAbsorptionCorrection)
			fProtonAbsorptionCorrection->FillAbsorptionMaps(fAOD);
		if(fIsOn_AliProtonFeedDownAnalysis)
			fProtonFeedDownAnalysis->Analyze(fAOD);
		if(fIsOn_AliProtonSpectraCorrection)
			fProtonSpectraCorrection->FillCorrectionMaps(fAOD);
	}//AOD analysis
	                 
	
	// Post output data.
	PostData(0, fList);
}    
//__________________________________________________________________________________________
void AliProtonCorrectionAnalysisTask::Terminate(Option_t *) 
{

}
//___________________________________________________________________________________________
void AliProtonCorrectionAnalysisTask::SetAnalysisObjectAbsorptionCorrection(AliProtonAbsorptionCorrection *const analysis) 
{
	if (analysis&&fProtonAnalysisBase)
	{
		Printf("Absorption Correection ON\n");
		fIsOn_AliProtonAbsorptionCorrection=kTRUE;
		fProtonAbsorptionCorrection = analysis;
		fProtonAbsorptionCorrection->SetBaseAnalysis(fProtonAnalysisBase);
		fProtonAbsorptionCorrection->InitAnalysisHistograms(fProtonAnalysisBase->GetNBinsX(),fProtonAnalysisBase->GetMinX(),fProtonAnalysisBase->GetMaxX(),fProtonAnalysisBase->GetNBinsY(),fProtonAnalysisBase->GetMinY(),fProtonAnalysisBase->GetMaxY());
	}
	else
		fIsOn_AliProtonAbsorptionCorrection=kFALSE;
}
//___________________________________________________________________________________________
void AliProtonCorrectionAnalysisTask::SetAnalysisObjectFeedDown(AliProtonFeedDownAnalysis *const analysis) 
{
	if (analysis&&fProtonAnalysisBase)
	{
		Printf("Feed Down ON\n");
		 fIsOn_AliProtonFeedDownAnalysis=kTRUE;
		fProtonFeedDownAnalysis = analysis;
		fProtonFeedDownAnalysis->SetBaseAnalysis(fProtonAnalysisBase);
		fProtonFeedDownAnalysis->InitAnalysisHistograms(fProtonAnalysisBase->GetNBinsX(),fProtonAnalysisBase->GetMinX(),fProtonAnalysisBase->GetMaxX(),fProtonAnalysisBase->GetNBinsY(),fProtonAnalysisBase->GetMinY(),fProtonAnalysisBase->GetMaxY());
	}
	else
		 fIsOn_AliProtonFeedDownAnalysis=kFALSE;
}
//___________________________________________________________________________________________
void AliProtonCorrectionAnalysisTask::SetAnalysisObjectSpectraCorrection(AliProtonSpectraCorrection *const analysis) 
{ 
	if (analysis&&fProtonAnalysisBase)
	{
		Printf("Spectra Correection ON\n");
		fIsOn_AliProtonSpectraCorrection=kTRUE;
		fProtonSpectraCorrection= analysis;
		fProtonSpectraCorrection->SetBaseAnalysis(fProtonAnalysisBase);
		fProtonSpectraCorrection->InitAnalysisHistograms(fProtonAnalysisBase->GetNBinsX(),fProtonAnalysisBase->GetMinX(),fProtonAnalysisBase->GetMaxX(),fProtonAnalysisBase->GetNBinsY(),fProtonAnalysisBase->GetMinY(),fProtonAnalysisBase->GetMaxY());
	}
	else
		fIsOn_AliProtonSpectraCorrection=kFALSE;
}
  



