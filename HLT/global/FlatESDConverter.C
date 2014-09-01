/**
 * >> Testing Macro to fill FlatESDEvent from ALiESDEvent <<
 **
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli
 *
 * Usage:
 *  aliroot -b -l -q LoadLibs.C FlatESDConverter.C++
 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "../TPC/Rec/AliTPCseed.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "./AliFlatESDEvent.h"
#include "./AliFlatESDTrack.h"
#include "./AliFlatESDTrigger.h"
#include "./AliFlatTPCCluster.h"
#include "./AliFlatExternalTrackParam.h"
#include "Riostream.h"
#include "AliSysInfo.h"
#endif   

void FlatESDConverter(const char* filename="AliESDs.root", const char* filenameFriends="AliESDfriends.root",const char* filenameOut="out.dat", Bool_t useESDFriends = kTRUE, Bool_t useHLTtree = kFALSE,Int_t verbose = 0) {

	if(useESDFriends) Printf("using friends");
	if(useHLTtree) Printf("using HLT tree");
	
	
  // -- Convert AliESDEvent to AliFlatESDEvent
/*
  if ( access( filename, F_OK ) == -1 ){
	Printf("input file not readable!");
	return;
  }
*/

  TFile *file    = new TFile(Form("%s", filename));
	
	
  ofstream outFile(Form("%s",filenameOut), std::ifstream::binary | std::ifstream::out);
  //ofstream outFile("outFlatESD.dat");

  
  TTree *esdTree = useHLTtree? dynamic_cast<TTree*>(file->Get("HLTesdTree")) : dynamic_cast<TTree*>(file->Get("esdTree"));
  
  // -- Connect ESD
  AliESDEvent *esd = new AliESDEvent;
  esd->ReadFromTree(esdTree);
  
  // -- Connect ESD friend
  AliESDfriend *esdFriend = NULL; 
  if (useESDFriends && !esdTree->FindBranch("ESDfriend.")) {
    esdTree->AddFriend("esdFriendTree", Form("%s", filenameFriends));
    esdTree->SetBranchStatus("ESDfriend.", 1);
    esdFriend = (AliESDfriend*)esd->FindListObject("AliESDfriend");
    if (esdFriend)
      esdTree->SetBranchAddress("ESDfriend.", &esdFriend);
  } // if (!esdTree->FindBranch("ESDfriend.")) {
  ;
  AliFlatESDEvent *flatEsd = NULL;

  // -- Event Loop
  for (Int_t idxEvent = 0; idxEvent < esdTree->GetEntries(); idxEvent++) {
    Printf("Processing event nr %d", idxEvent);
  //  esd->SaveAs("esdTemp.root");
   // TFile  fTmp = TFile("esdTemp.root");
    Int_t sizeIn = 1;//fTmp.GetSize();

    AliSysInfo::AddStamp("getEntry",0,0,idxEvent);

    esdTree->GetEntry(idxEvent);
    // -- Book memory for AliFlatESDEvent
    // -- TEST >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	Int_t size = AliFlatESDEvent::EstimateSize(esd, kTRUE);
    Byte_t *mem = new Byte_t[size];

    flatEsd = reinterpret_cast<AliFlatESDEvent*>(mem);
    new (flatEsd) AliFlatESDEvent;
	
	
	
    AliSysInfo::AddStamp("DoEvent.Start",0,0,idxEvent);
    // -- Fill AliFlatESDEvent

  Int_t err=  flatEsd->SetFromESD( size,  esd, kTRUE ); 
  
  if(err) Printf("!!! Error while filling flatESD event  %d!!!", err);
  Printf("trigger classes: %d size: %d , = %d",flatEsd->GetNumberOfTriggerClasses(), sizeof(AliFlatESDTrigger) , flatEsd->GetNumberOfTriggerClasses() *sizeof(AliFlatESDTrigger) );      
    AliSysInfo::AddStamp("DoEvent.Stop",sizeIn,flatEsd->GetSize(),idxEvent);

		if(useESDFriends){
      Printf("ESD : Event %d || V0s %d || Tracks %d | FRIEND Tracks %d || estimated size %llu", 
	     idxEvent, esd->GetNumberOfV0s(),esd->GetNumberOfTracks(), esdFriend->GetNumberOfTracks(), 
	     AliFlatESDEvent::EstimateSize(esd, useESDFriends));
      Printf("FLAT: Event %d || V0s %d || Tracks %d | FRIEND Tracks %d || estimated size %llu", 
	     idxEvent, flatEsd->GetNumberOfV0s(),flatEsd->GetNumberOfTracks(), esdFriend->GetNumberOfTracks(), flatEsd->GetSize());
		}
		else{
      Printf("ESD : Event %d || V0s %d || Tracks %d || estimated size %llu ", 
	     idxEvent,esd->GetNumberOfV0s(), esd->GetNumberOfTracks(), 
	     AliFlatESDEvent::EstimateSize(esd, useESDFriends) );
      Printf("FLAT: Event %d || V0s %d || Tracks %d || estimated size %llu ", 
	     idxEvent, flatEsd->GetNumberOfV0s(),flatEsd->GetNumberOfTracks(),  flatEsd->GetSize());
		}
      AliFlatESDTrack *track = const_cast<AliFlatESDTrack*> (flatEsd->GetTracks());
      for (Int_t idxTrack = 0; idxTrack < flatEsd->GetNumberOfTracks(); ++idxTrack) {      
	AliESDtrack *esdTrack = esd->GetTrack(idxTrack);      
	AliESDfriendTrack *friendTrack = useESDFriends ?  esdFriend->GetTrack(idxTrack) :NULL;
       	if (track && !esdTrack) {
	  Printf("ERROR THIS SHOULD NOT HAPPEN AT ALL !!! TRACK %d HAS NO ESD TRACK!!!", idxTrack);
	  return;
	}
    if (verbose) {
	if (track) {
	  const AliFlatExternalTrackParam* exp1 = track->GetFlatTrackParamCp();
	   const AliFlatExternalTrackParam* exp2 = track->GetFlatTrackParamIp();
	   const AliFlatExternalTrackParam* exp3 = track->GetFlatTrackParamTPCInner();
	   const AliFlatExternalTrackParam* exp4 = track->GetFlatTrackParamOp();

	  Float_t alphaFLAT[4] = {0., 0., 0., 0.};
	  if (exp1) alphaFLAT[0] = exp1->GetAlpha();
	  if (exp2) alphaFLAT[1] = exp2->GetAlpha();
	  if (exp3) alphaFLAT[2] = exp3->GetAlpha();
	  if (exp4) alphaFLAT[3] = exp4->GetAlpha();

	  Float_t alphaOLD[4]  = {0., 0., 0., 0.};
	  if (esdTrack->GetConstrainedParam())  alphaOLD[0] = esdTrack->GetConstrainedParam()->GetAlpha();
	  if (esdTrack->GetInnerParam())        alphaOLD[1] = esdTrack->GetInnerParam()->GetAlpha();
	  if (esdTrack->GetTPCInnerParam())     alphaOLD[2] = esdTrack->GetTPCInnerParam()->GetAlpha();
	  if (esdTrack->GetOuterParam())        alphaOLD[3] = esdTrack->GetOuterParam()->GetAlpha();

	  Printf("TEST: FlatTrack %d > FlatExternalTrackParam > %p %p %p %p", idxTrack, exp1, exp2, exp3, exp4);
	  Printf("TEST: FlatTrack %d > Alpha %f %f %f %f", idxTrack, alphaFLAT[0], alphaFLAT[1], alphaFLAT[2], alphaFLAT[3]);
	  Printf("TEST: Old Track %d > Alpha %f %f %f %f", idxTrack, alphaOLD[0], alphaOLD[1], alphaOLD[2], alphaOLD[3]);
	  Printf("TEST: Diff      %d > Alpha %f %f %f %f", idxTrack, 
		 alphaFLAT[0]-alphaOLD[0], alphaFLAT[1]-alphaOLD[1], alphaFLAT[2]-alphaOLD[2], alphaFLAT[3]-alphaOLD[3]);

	  Int_t nCl = track->GetNumberOfTPCClusters();
	  Printf("TEST: FlatTrack %d has %d FlatClusters", idxTrack, nCl);
	  
#if 0
	  if(nCl && useESDFriends && verbose > 1){
	    TObject* calibObject = NULL;
	    AliTPCseed* seed = NULL;
	    for (Int_t idx = 0; (calibObject = friendTrack->GetCalibObject(idx)); ++idx) {
	      if ((seed = dynamic_cast<AliTPCseed*>(calibObject))) break;
	    }
	    // -- Fill cluster
	    if (seed) {
    	      Int_t idxRow2=0;
	      for (Int_t idxRow = 0; idxRow <  nCl; idxRow++){
		AliFlatTPCCluster * cl = track->GetTPCCluster(idxRow);
		cout << " idx fX fY fZ  fSigmaY2 fSigmaZ2 fCharge fQMax fPadRow" << endl;
		if(cl){
		  cout << idxRow << " " << cl->GetX() << " " << cl->GetY() << " " << cl->GetZ() << " " << cl->GetSigmaY2() << " " << cl->GetSigmaZ2() << " " << cl->GetCharge() << " " << cl->GetQMax() << " " << cl->GetPadRow() << endl;
		}
		else{
		  cout << idxRow << "---------------------------------" << endl << endl;		
		}
		AliTPCclusterMI* cl2 = NULL; 
		while(!cl2 && idxRow2<160){
		  cl2 = seed->GetClusterPointer(idxRow2++);
		}
		if (cl2) {
		  //cout<<" normalCl fX fY fZ fPadRow fSigmaY2 fSigmaZ2 fCharge fQMax" <<endl;
		  cout << idxRow << " " << cl2->GetX() << " " << cl2->GetY() << " " << cl2->GetZ() << " " << cl2->GetSigmaY2() << " " << cl2->GetSigmaZ2() << " " << cl2->GetQ() << " " << cl2->GetMax() << " " << cl2->GetRow() << endl << endl;
		}
		else
		  cout << idxRow << "---------------------------------" << endl << endl;	
	      }
	    }
	  }
#endif

	}
	track = const_cast<AliFlatESDTrack*> (track->GetNextTrack() );
      }
    }

    outFile.write(reinterpret_cast<char*>(mem), flatEsd->GetSize());
    delete[] mem; 

  } // for (Int_t idxEvent = 1; idxEvent < 2; idxEvent++) {

  outFile.close();

  return;
}    
