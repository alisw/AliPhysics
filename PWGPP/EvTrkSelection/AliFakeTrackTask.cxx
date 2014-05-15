#include "AliFakeTrackTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"


//#include "AliESDtrack.h"

#include "AliInputEventHandler.h"
#include "AliStack.h"
//#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TH3F.h"
//#include "TMCProcess.h"
#include "AliVEvent.h"

#include "AliESDtrackCuts.h"
//#include "AliESDpidCuts.h"
//#include "AliESDpid.h"
#include "AliCentrality.h"
#include <iostream>
#include "TChain.h"

using namespace std;

ClassImp(AliFakeTrackTask)

AliFakeTrackTask::AliFakeTrackTask(const char *name ):
AliAnalysisTaskSE(name),fESD(0),
fptvsTPCsignalvsITSsignalAll(0),fptvsTPCsignalvsITSsignalGlobalgood(0), 
fptvsTPCsignalvsITSsignalGlobalfake(0),fptvsTPCsignalvsITSsignalTPCfake(0),
fptvsTPCsignalvsITSsignalITSfake(0),ffakestat(0),
ftrackcuts(0),flistout(0)	
{
	 DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
}
//___________________________________________________________________________
AliFakeTrackTask:: ~AliFakeTrackTask()
{
	if(fptvsTPCsignalvsITSsignalAll)
		delete fptvsTPCsignalvsITSsignalAll;
 	if(fptvsTPCsignalvsITSsignalGlobalgood)
 		delete fptvsTPCsignalvsITSsignalGlobalgood;
	if(fptvsTPCsignalvsITSsignalGlobalfake)
		delete fptvsTPCsignalvsITSsignalGlobalfake;
 	if(fptvsTPCsignalvsITSsignalTPCfake)
		delete  fptvsTPCsignalvsITSsignalTPCfake;
	if(fptvsTPCsignalvsITSsignalITSfake)
		delete fptvsTPCsignalvsITSsignalITSfake;
 	if(ffakestat)
		delete ffakestat;	
	if(ftrackcuts)
	 	delete  ftrackcuts;
 	if(flistout)
		delete flistout;

}
//____________________________________________________________________
void AliFakeTrackTask::UserCreateOutputObjects()
{
	flistout=new TList();
	flistout->SetOwner(kTRUE);
	fptvsTPCsignalvsITSsignalAll=new TH3F("ptvsTPCsignalvsITSsignalAll","ptvsTPCsignalvsITSsignalAll",50,0,5.0,200,0,1000,200,0,1000);
 	fptvsTPCsignalvsITSsignalGlobalgood=new TH3F("ptvsTPCsignalvsITSsignalGlobalgood","ptvsTPCsignalvsITSsignalGlobalgood",50,0,5.0,200,0,1000,200,0,1000);
 	fptvsTPCsignalvsITSsignalGlobalfake=new TH3F("ptvsTPCsignalvsITSsignalGlobalfake","ptvsTPCsignalvsITSsignalGlobalfake",50,0,5.0,200,0,1000,200,0,1000);
 	fptvsTPCsignalvsITSsignalTPCfake=new TH3F("ptvsTPCsignalvsITSsignalTPCfake","ptvsTPCsignalvsITSsignalTPCfake",50,0,5.0,200,0,1000,200,0,1000);
 	fptvsTPCsignalvsITSsignalITSfake=new TH3F("ptvsTPCsignalvsITSsignalITSfake","ptvsTPCsignalvsITSsignalITSfake",50,0,5.0,200,0,1000,200,0,1000);
	ffakestat= new TH1F("fake stats","fake stats",8,-0.5,7.5);
	flistout->Add(fptvsTPCsignalvsITSsignalAll);
	flistout->Add(fptvsTPCsignalvsITSsignalGlobalgood);
	flistout->Add(fptvsTPCsignalvsITSsignalGlobalfake);
	flistout->Add(fptvsTPCsignalvsITSsignalTPCfake);
	flistout->Add(fptvsTPCsignalvsITSsignalITSfake);
	flistout->Add(ffakestat);
	PostData(1,  flistout);	

}
//_______________________________________________________________________
void AliFakeTrackTask::UserExec(Option_t *)
{
	fESD = dynamic_cast<AliESDEvent*> (InputEvent());
	if (!fESD) 
	{
		Printf("ERROR: fESD not available");
		return;
	}
	//UInt_t isSelected=0;
/*	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())))
		isSelected=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
     	if(!(isSelected&AliVEvent::kMB))
     	{
		 PostData(1,  flistout);
		 return;
	}*/

	const AliESDVertex *vertex = 0x0;
	vertex = fESD->GetPrimaryVertexTracks();
	if(vertex->GetNContributors()<1) 
	{
		// SPD vertex
		vertex = fESD->GetPrimaryVertexSPD();
		if(vertex->GetNContributors()<1) 
		{
			//Printf("No good  Vertex.........\n");
			PostData(1,  flistout);
			return;			
		}	
	}
	if(fESD->IsPileupFromSPDInMultBins())
	{
		PostData(1,  flistout);
		return;
	}		
	if(TMath::Abs(vertex ->GetZ())>10.0)
	{
			PostData(1,  flistout);
			return;
	}	

	if(!fPIDResponse) 
  	{
    		AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   		AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    		fPIDResponse = inputHandler->GetPIDResponse();
  	}

	AliStack* stack=0x0;
	AliMCEvent* mcEvent  = (AliMCEvent*) MCEvent();
	if(mcEvent)	
		stack = mcEvent->Stack();
	
	Int_t nTracks=fESD->GetNumberOfTracks();
	for (int i=0;i<nTracks;i++)
	{
		AliESDtrack* esdtrack=fESD->GetTrack(i);
		//Printf("%d %d",i,nTracks);
		if(!ftrackcuts->AcceptTrack(esdtrack))
			continue;
		Int_t nSDDSSD=0;
		for (int j=0;j<4;j++)
		{
			if(esdtrack->HasPointOnITSLayer(i))
				nSDDSSD++;
		}	 
		//Printf("%d",nSDDSSD);

		if(nSDDSSD<3)
			continue;
		Float_t pt=esdtrack->Pt();
		Float_t signalTPC=esdtrack->GetTPCsignal();
		Float_t singalITS=esdtrack->GetITSsignal();
		Int_t faketype=0;
		fptvsTPCsignalvsITSsignalAll->Fill(pt,signalTPC,singalITS);
		if(!stack)
			continue;
		if(esdtrack->GetLabel()>0)
			fptvsTPCsignalvsITSsignalGlobalgood->Fill(pt,signalTPC,singalITS);
		else
		{
			fptvsTPCsignalvsITSsignalGlobalfake->Fill(pt,signalTPC,singalITS);
			faketype+=1;

		}	
		if(esdtrack->GetITSLabel()<0)
		{
			fptvsTPCsignalvsITSsignalITSfake->Fill(pt,signalTPC,singalITS);
			faketype+=2;
		}
		if(esdtrack->GetTPCLabel()<0)
		{
			fptvsTPCsignalvsITSsignalTPCfake->Fill(pt,signalTPC,singalITS);

			faketype+=4;
		}
		ffakestat->Fill(faketype);
		

	}	
	PostData(1,  flistout);

} 
