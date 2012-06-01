#include <Riostream.h>
#include <TFile.h>
#include "AliHFAssociatedTrackCuts.h"
#include <TClonesArray.h>
#include <TParameter.h>
//#include "AliAODPidHF.h"

/* $Id$ */

void makeInputHFCorrelation(){
	
	
	AliHFAssociatedTrackCuts* HFCorrelationCuts=new AliHFAssociatedTrackCuts();
	HFCorrelationCuts->SetName("AssociatedCuts");
	HFCorrelationCuts->SetTitle("Cuts for associated track");
	Float_t eta = 0.9;

	//______________________________ set ESD track cuts
	AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();

	esdTrackCuts->SetRequireTPCRefit(kTRUE);
	esdTrackCuts->SetRequireITSRefit(kTRUE);
	esdTrackCuts->SetMinNClustersTPC(80);
	esdTrackCuts->SetMinNClustersITS(2); 
	//esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); // set requirement on Pixels
	// default is kBoth, otherwise kAny
	esdTrackCuts->SetPtRange(0.3,1.e10);
	esdTrackCuts->SetEtaRange(-eta,eta);
	HFCorrelationCuts->AddTrackCuts(esdTrackCuts);
	
	//______________________________ set kinematics cuts for AOD track 
	const int nofcuts = 4;
	Float_t* trackcutsarray;
	trackcutsarray=new Float_t[nofcuts];
	trackcutsarray[0] = 0.3;//track min pt
	trackcutsarray[1] = 100.;//track max pt
	trackcutsarray[2] = 0.;//track min impact parameter
	trackcutsarray[3] = 100.;//track max impact parameter
	
	HFCorrelationCuts->SetNVarsTrack(nofcuts);
	HFCorrelationCuts->SetAODTrackCuts(trackcutsarray);
	
	
	
	//______________________________ set kinematics cuts for AOD v0 
	
	const int nofcuts2 = 7;
	
	Float_t* vzerocutsarray;
	vzerocutsarray=new Float_t[nofcuts2];
	vzerocutsarray[0] = 0.2; // max dca between two daugters (cm)
	vzerocutsarray[1] = 2; //  max chi square
	vzerocutsarray[2] = 2.; // min decay length (cm) 
	vzerocutsarray[3] = 15; // max decay length (cm)
	vzerocutsarray[4] = 0.2; // min opening angle between two daugters
	vzerocutsarray[5] = 0; // min pt of k0 (GeV/c)
	vzerocutsarray[6] = 0.9; // set eta acceptance

	HFCorrelationCuts->SetNVarsVzero(nofcuts2);
	HFCorrelationCuts->SetAODvZeroCuts(vzerocutsarray);
	
	//______________________________ set PID
	 
	Int_t mode =1;
	AliAODPidHF* pidObj=new AliAODPidHF();
	pidObj->SetMatch(mode);
	pidObj->SetSigma(0,2); // TPC
	pidObj->SetSigma(3,3); // TOF
	pidObj->SetTPC(kTRUE);
	pidObj->SetTOF(kTRUE);
	pidObj->SetCompat(kTRUE);
	
	HFCorrelationCuts->SetPidHF(pidObj);
	
	
    //______________________________ save to *.root file
	HFCorrelationCuts->PrintAll();
	
	TFile* fout=new TFile("AssocPartCuts.root","recreate");   //set this!! 
	fout->cd();
	HFCorrelationCuts->Write();
	fout->Close();


}
