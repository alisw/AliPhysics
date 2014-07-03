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
#include "./AliFlatTPCCluster.h"
#include "./AliFlatExternalTrackParam.h"
#include "Riostream.h"
#endif   

void FlatESDConverter(const char* filename="AliESDs.root", const char* filenameFriends="AliESDfriends.root",const char* filenameOut="out.dat", Bool_t useESDFriends = kTRUE, Bool_t useHLTtree = kFALSE,Int_t verbose = 0) {
  // -- Convert AliESDEvent to AliFlatESDEvent

  ofstream outFile(Form("%s",filenameOut), std::ifstream::binary | std::ifstream::out);
  //ofstream outFile("outFlatESD.dat");

  TFile *file    = new TFile(Form("%s", filename));
	
	TTree *esdTree = useHLTtree? dynamic_cast<TTree*>(file->Get("HLTesdTree")) : dynamic_cast<TTree*>(file->Get("esdTree"));

  
  // -- Connect ESD
  AliESDEvent *esd = new AliESDEvent;
  esd->ReadFromTree(esdTree);
  
  // -- Connect ESD friend
  AliESDfriend *esdFriend = NULL; 
  if (useESDFriends && !esdTree->FindBranch("ESDfriend.")) {
    esdTree->AddFriend("esdFriendTree", Form("%s", filenameFriends));
    esdTree->SetBranchStatus("ESDfriend.", 1);

    esdFriend = dynamic_cast<AliESDfriend*>((const_cast<AliESDEvent*>(esd))->FindListObject("AliESDfriend"));
    if (esdFriend)
      esdTree->SetBranchAddress("ESDfriend.", &esdFriend);
  } // if (!esdTree->FindBranch("ESDfriend.")) {
  
  AliFlatESDEvent *flatEsd = NULL;

  // -- Event Loop
  for (Int_t idxEvent = 0; idxEvent < esdTree->GetEntries(); idxEvent++) {
  
  AliSysInfo::AddStamp("getEntry",0,0,edxEvent);
    esdTree->GetEntry(idxEvent);
    // -- Book memory for AliFlatESDEvent
   // -- TEST >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    Byte_t *mem = new Byte_t[AliFlatESDEvent::EstimateSize(esd, useESDFriends)];
    
AliSysInfo::AddStamp("DoEvent.Start",0,0,edxEvent);
    flatEsd = reinterpret_cast<AliFlatESDEvent*>(mem);    
	new (flatEsd) AliFlatESDEvent(1);

    // -- Fill AliFlatESDEvent
    flatEsd->Fill(esd, useESDFriends);  
    
    
AliSysInfo::AddStamp("DoEvent.Start",0,0,edxEvent);
if(verbose){
     Printf("TEST: Event %d || Tracks %d | FRIEND Tracks %d || estimated size %llu || sizeof(AliFlatESDEvent) %llu", 
	   idxEvent, esd->GetNumberOfTracks(), esdFriend->GetNumberOfTracks(), 
	   AliFlatESDEvent::EstimateSize(esd, useESDFriends), flatEsd->GetSize());


    AliFlatESDTrack *track = flatEsd->GetTracks();
    for (Int_t idxTrack = 0; idxTrack < flatEsd->GetNumberOfTracks(); ++idxTrack) {      
       AliESDtrack *esdTrack = esd->GetTrack(idxTrack);      
       AliESDfriendTrack *friendTrack = esdFriend->GetTrack(idxTrack) ;
       
       if (track && !esdTrack) {
	 Printf("ERROR THIS SHOULD NOT HAPPEN AT ALL !!! TRACK %d HAS NO ESD TRACK!!!", idxTrack);
	 return;
       }

      if (track ) {
	
	/*
	AliFlatExternalTrackParam* exp1 = track->GetTrackParamCp();
	AliFlatExternalTrackParam* exp2 = track->GetTrackParamIp();
	AliFlatExternalTrackParam* exp3 = track->GetTrackParamTPCInner();
	AliFlatExternalTrackParam* exp4 = track->GetTrackParamOp();

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

     	Printf("  TEST: FlatTrack %d > FlatExternalTrackParam > %p %p %p %p", idxTrack, exp1, exp2, exp3, exp4);
	Printf("  TEST: FlatTrack %d > Alpha %f %f %f %f", idxTrack, alphaFLAT[0], alphaFLAT[1], alphaFLAT[2], alphaFLAT[3]);
	Printf("  TEST: Old Track %d > Alpha %f %f %f %f", idxTrack, alphaOLD[0], alphaOLD[1], alphaOLD[2], alphaOLD[3]);
	Printf("  TEST: Diff      %d > Alpha %f %f %f %f", idxTrack, 
	       alphaFLAT[0]-alphaOLD[0], alphaFLAT[1]-alphaOLD[1], alphaFLAT[2]-alphaOLD[2], alphaFLAT[3]-alphaOLD[3]);
*/

Int_t nCl = track->GetNumberOfTPCClusters();
	Printf("  TEST: FlatTrack %d has %d FlatClusters", idxTrack,  nCl );
	if(nCl && verbose > 1){
	
	  TObject* calibObject = NULL;
    AliTPCseed* seed = NULL;

    for (Int_t idx = 0; (calibObject = friendTrack->GetCalibObject(idx)); ++idx) {
      if ((seed = dynamic_cast<AliTPCseed*>(calibObject))) break;
    }

    // -- Fill cluster
    if (seed) {
      for (Int_t idxRow = 0; idxRow <  nCl; idxRow++){
      AliFlatTPCCluster * cl = track->GetTPCCluster(idxRow);

			 	cout<<" idx fX fY fZ  fSigmaY2 fSigmaZ2 fCharge fQMax fPadRow" <<endl;
				cout<< idxRow <<" "<< cl->GetX()<<" "<< cl->GetY()<<" "<< cl->GetZ()<<" "<< cl->GetSigmaY2()<<" "<< cl->GetSigmaZ2()<<" "<< cl->GetCharge()<<" "<< cl->GetQMax() <<" "<< cl->GetPadRow()<<endl;
				
				AliTPCclusterMI* cl2 =seed->GetClusterPointer( cl->GetPadRow());
      
	if (cl2) {
	  //cout<<" normalCl fX fY fZ fPadRow fSigmaY2 fSigmaZ2 fCharge fQMax" <<endl;
				cout<< idxRow <<" "<< cl2->GetX()<<" "<< cl2->GetY()<<" "<< cl2->GetZ()<<" "<< cl2->GetSigmaY2()<<" "<< cl2->GetSigmaZ2()<<" "<< cl2->GetQ()<<" "<< cl2->GetMax()<<" "<< cl2->GetRow()<<" "<< cl2->GetPad() <<endl<<endl;
	}
		else
		  	cout<<idxRow<<"---------------------------------"<<endl<<endl;	
      }
    }
	}
	
	
	
/*
	for (Int_t idxCluster = 0; idxCluster < track->GetNumberOfTPCClusters(); ++idxCluster){
	  
	  }
	*/  
      }
      track = track->GetNextTrack();
    }
}

    outFile.write(reinterpret_cast<char*>(mem), flatEsd->GetSize());

    delete[] mem; 

  } // for (Int_t idxEvent = 1; idxEvent < 2; idxEvent++) {

  outFile.close();

  return;
}    
