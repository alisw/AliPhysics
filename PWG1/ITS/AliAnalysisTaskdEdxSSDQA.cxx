#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include "AliAnalysisTaskdEdxSSDQA.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"


#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "Riostream.h"

#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliTrackPointArray.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "TMath.h"


using namespace std;
ClassImp(AliAnalysisTaskdEdxSSDQA)

//________________________________________________________________________
AliAnalysisTaskdEdxSSDQA::AliAnalysisTaskdEdxSSDQA(const char *name) 
: AliAnalysisTaskSE(name), fHist1(0),fHist2(0),fListOfHistos(0),fPcut(0.0)
{
  // Constructor
	// Define input and output slots here
  // Input slot #1 works with a TChain
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskdEdxSSDQA::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  	
	fListOfHistos = new TList();
	fHist1 =new TH2F("QAChargeRatio","QAChargeRatio;Module;CR",1698,-0.5,1697.5,80,-1.0,1.0);
	fListOfHistos->Add(fHist1);
	fHist2=new TH2F("QACharge","QACharge;Module;Q",1698,-0.5,1697.5,150,0,300);
	fListOfHistos->Add(fHist2);
}
//________________________________________________________________________
void AliAnalysisTaskdEdxSSDQA::LocalInit() 
{
	
	Printf("end of LocalInit");
}

//________________________________________________________________________

void AliAnalysisTaskdEdxSSDQA::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
    AliESDEvent* esd = dynamic_cast<AliESDEvent*> (InputEvent());
    if(!esd) 
    {
	Printf("ERROR: Input ESD Event not available");
	return;
    }
    
    Int_t nTracks=esd->GetNumberOfTracks();
    
    

	if(!ESDfriend())
	{
		Printf("problem with friend");
		return;
	}
	
	AliTrackPointArray*  trackar=0x0;
	Bool_t l5;
	Bool_t l6;
	Int_t npoints;
	AliTrackPoint point;
	Int_t nlayer=0;
	Int_t id=0;
	Float_t chargeratio=0;
	for(int itr=0;itr<nTracks;itr++)
	{
		AliESDtrack* track= esd->GetTrack(itr);

		if(TMath::Abs(track->Eta())>0.9)
			continue;
		if(track->GetP()>10.0)
			continue;	
		if (track->IsOn(AliESDtrack::kITSrefit)&&track->IsOn(AliESDtrack::kTPCrefit))
		{
			 l5=track->HasPointOnITSLayer(4);
			 l6=track->HasPointOnITSLayer(5);
			if (!(l5||l6))//only tracks with SSD point
				continue;
			Double_t tmpQESD[4]={-1.0,-1.0,-1.0,-1.0};
			track->GetITSdEdxSamples(tmpQESD);
			trackar=(AliTrackPointArray*)track->GetTrackPointArray(); 
			if(!trackar)
				continue;
			npoints=trackar->GetNPoints();			
			for(int itnp=0;itnp<npoints;itnp++)		
			{									
				if(trackar->GetPoint(point,itnp))
					nlayer=AliGeomManager::VolUIDToLayer(point.GetVolumeID(), id);//layer number
				else
					continue;		
				if(nlayer==5||nlayer==6)
				{
					if(point.GetCharge()>0.0&&point.IsExtra()==kFALSE)//do not use additional clusters
					{		
						 chargeratio=point.GetChargeRatio();						 
						if(nlayer==5&&tmpQESD[2]>0.0)
						{
							fHist1->Fill(id,chargeratio);
							if(track->GetP()>fPcut)
								fHist2->Fill(id,tmpQESD[2]);
						}	
						if(nlayer==6&&tmpQESD[3]>0.0)
						{
							fHist1->Fill(id+748,chargeratio);
							if(track->GetP()>fPcut)
								fHist2->Fill(id+748,tmpQESD[3]);
						}		
					}	
				}
			}	
		}
	}
	
	// Post output data.
	PostData(1,  fListOfHistos);
}      

//________________________________________________________________________
void AliAnalysisTaskdEdxSSDQA::Terminate(Option_t *) 
{

	Printf("end of Terminate");
}
//_____________________________________________________________________________



