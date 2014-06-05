/**
 * >> Testing Macro to compare FlatESDEvents from output files <<
 **
 * Primary Authors : Steffen Weber
 *
 * Usage:
 *  aliroot -b -l -q LoadLibs.C CompareFlatESDs.C++
 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "./AliFlatESDEvent.h"
#include "./AliFlatESDTrack.h"
#include "./AliFlatTPCCluster.h"
#include "./AliFlatExternalTrackParam.h"
#include "Riostream.h"
#endif   

void CompareFlatESDs(const char* filename1="outFlatESD1.root",const char* filename2="outFlatESD2.root") {
  
  
  // Create output histograms
  
    TH2F* hNTracks;
  
  TString outputFilename = "$PWD/compare.root";
  
	cout<< "creating histograms"<<endl;
	hNTracks = new TH2F("nTracks","number of tracks", 100,0,100, 100,0,100);
  
  

  ifstream is1(filename1, std::ifstream::binary | std::ifstream::in);
  ifstream is2(filename2, std::ifstream::binary | std::ifstream::in);
  if (is1 && is2 ){
    is1.seekg (0, is1.end);
    int length1 = is1.tellg();
    is1.seekg (0, is1.beg);
    char * buffer1 = new char [length1];
    
    std::cout << "Reading " << length1 << " characters... ";
    
    is1.read (buffer1,length1);
    if (is1)
      std::cout << "all characters read successfully." << endl;
    else
      std::cout << "error: only " << is1.gcount() << " could be read";
    is1.close();
	
	
    is2.seekg (0, is2.end);
    int length2 = is2.tellg();
    is2.seekg (0, is2.beg);
    char * buffer2 = new char [length2];
    
    std::cout << "Reading " << length2 << " characters... ";
    
    is2.read (buffer2,length2);
    if (is2)
      std::cout << "all characters read successfully." << endl;
    else
      std::cout << "error: only " << is2.gcount() << " could be read";
    is2.close();
    
	
	
    // ...buffer contains the entire file...
    
    char *curr1 = buffer1;
    char *endBuff1 = buffer1+length1;
	
    char *curr2 = buffer2;
    char *endBuff2 = buffer2+length2;
	
    int iEvent = 0;
    
	while( curr1 < endBuff1  && curr2 < endBuff2 ){
      AliFlatESDEvent *flatEsd1 = reinterpret_cast<AliFlatESDEvent *>(curr1);
      AliFlatESDEvent *flatEsd2 = reinterpret_cast<AliFlatESDEvent *>(curr2);
	  
	  if(flatEsd1->GetNumberOfTracks()==0  ||  flatEsd2->GetNumberOfTracks() ==0){
		if(flatEsd1->GetNumberOfTracks()==0){
		  curr1=curr1+ flatEsd1->GetSize();
		}
		if(flatEsd2->GetNumberOfTracks()==0){
		  curr2=curr2+ flatEsd2->GetSize();
		}
		continue;
	  }
	  
	  
      cout<<endl<<"Reading event "<<iEvent<<":"<<endl;
	  
	  
	  if( (Bool_t) flatEsd1->GetPrimaryVertexSPD() != (Bool_t) flatEsd2->GetPrimaryVertexSPD()  ) cout<<"DIFFERENCE!: ";
      cout<<"vtx SPD: "<<(Bool_t) flatEsd1->GetPrimaryVertexSPD() << " | " << (Bool_t) flatEsd2->GetPrimaryVertexSPD() << endl;
	 
	  if( (Bool_t) flatEsd1->GetPrimaryVertexTracks() != (Bool_t) flatEsd2->GetPrimaryVertexTracks()  ) cout<<"DIFFERENCE!: ";
	  cout<<" vtx tracks: "<<(Bool_t) flatEsd1->GetPrimaryVertexTracks()<< " | " <<	(Bool_t) flatEsd2->GetPrimaryVertexTracks()<< endl;
	  
	  if(  flatEsd1->GetNumberOfTracks() != flatEsd2->GetNumberOfTracks() ) cout<<"DIFFERENCE!: ";
	  cout<<" ntracks: "<<flatEsd1->GetNumberOfTracks()<< " | " <<flatEsd2->GetNumberOfTracks()<< endl;
	  
	  if(  flatEsd1->GetNumberOfV0s() != flatEsd2->GetNumberOfV0s() ) cout<<"DIFFERENCE!: ";
	  cout<<" nV0's: "<<flatEsd1->GetNumberOfV0s()<< " | " <<flatEsd2->GetNumberOfV0s()<< endl;
	  
	  
	  
	  // compare tracks
	  
	  
	  AliFlatESDTrack *track1 = flatEsd1->GetTracks();
	  AliFlatESDTrack *track2 = flatEsd2->GetTracks();
    for (Int_t idxTrack = 0; idxTrack < flatEsd1->GetNumberOfTracks(); ++idxTrack) { 

      if (track1 && track2) {
	AliFlatExternalTrackParam* exp11 = track1->GetTrackParamCp();
	AliFlatExternalTrackParam* exp21 = track1->GetTrackParamIp();
	AliFlatExternalTrackParam* exp31 = track1->GetTrackParamTPCInner();
	AliFlatExternalTrackParam* exp41 = track1->GetTrackParamOp();
	
	AliFlatExternalTrackParam* exp12 = track2->GetTrackParamCp();
	AliFlatExternalTrackParam* exp22 = track2->GetTrackParamIp();
	AliFlatExternalTrackParam* exp32 = track2->GetTrackParamTPCInner();
	AliFlatExternalTrackParam* exp42 = track2->GetTrackParamOp();

	Float_t alphaFLAT1[4] = {-99., -99., -99., -99.};
	if (exp11) alphaFLAT1[0] = exp11->GetAlpha();
	if (exp21) alphaFLAT1[1] = exp21->GetAlpha();
	if (exp31) alphaFLAT1[2] = exp31->GetAlpha();
	if (exp41) alphaFLAT1[3] = exp41->GetAlpha();
	
	Float_t alphaFLAT2[4] = {-99.,-99., -99., -99.};
	if (exp12) alphaFLAT2[0] = exp12->GetAlpha();
	if (exp22) alphaFLAT2[1] = exp22->GetAlpha();
	if (exp32) alphaFLAT2[2] = exp32->GetAlpha();
	if (exp42) alphaFLAT2[3] = exp42->GetAlpha();
	
	
	Float_t pFLAT1[4] = {-99., -99., -99., -99.};
	if (exp11) pFLAT1[0] = exp11->GetSigned1Pt();
	if (exp21) pFLAT1[1] = exp21->GetSigned1Pt();
	if (exp31) pFLAT1[2] = exp31->GetSigned1Pt();
	if (exp41) pFLAT1[3] = exp41->GetSigned1Pt();
	
	Float_t pFLAT2[4] = {-99., -99., -99., -99.};
	if (exp12) pFLAT2[0] = exp12->GetSigned1Pt();
	if (exp22) pFLAT2[1] = exp22->GetSigned1Pt();
	if (exp32) pFLAT2[2] = exp32->GetSigned1Pt();
	if (exp42) pFLAT2[3] = exp42->GetSigned1Pt();
	
	
	
	  if( alphaFLAT1[0] != alphaFLAT2[0] ) cout<<"DIFFERENCE!: ";
	  cout<<" alpha0: "<<alphaFLAT1[0]<< " | " <<alphaFLAT2[0]<< endl;
	  if( alphaFLAT1[1] != alphaFLAT2[1] ) cout<<"DIFFERENCE!: ";
	  cout<<" alpha1: "<<alphaFLAT1[1]<< " | " <<alphaFLAT2[1]<< endl;
	  if( alphaFLAT1[2] != alphaFLAT2[2] ) cout<<"DIFFERENCE!: ";
	  cout<<" alpha2: "<<alphaFLAT1[2]<< " | " <<alphaFLAT2[2]<< endl;
	  if( alphaFLAT1[3] != alphaFLAT2[3] ) cout<<"DIFFERENCE!: ";
	  cout<<" alpha3: "<<alphaFLAT1[3]<< " | " <<alphaFLAT2[3]<< endl;
	
	
	
	  if( pFLAT1[0] != pFLAT2[0] ) cout<<"DIFFERENCE!: ";
	  cout<<" p0: "<<pFLAT1[0]<< " | " <<pFLAT2[0]<< endl;
	  if( pFLAT1[1] != pFLAT2[1] ) cout<<"DIFFERENCE!: ";
	  cout<<" p1: "<<pFLAT1[1]<< " | " <<pFLAT2[1]<< endl;
	  if( pFLAT1[2] != pFLAT2[2] ) cout<<"DIFFERENCE!: ";
	  cout<<" p2: "<<pFLAT1[2]<< " | " <<pFLAT2[2]<< endl;
	  if( pFLAT1[3] != pFLAT2[3] ) cout<<"DIFFERENCE!: ";
	  cout<<" p3: "<<pFLAT1[3]<< " | " <<pFLAT2[3]<< endl;
	  
	  
	  if( track1->GetNumberOfTPCClusters() != track2->GetNumberOfTPCClusters() ) cout<<"DIFFERENCE!: ";
	  cout<<" nTPCclusters: "<<track1->GetNumberOfTPCClusters()<< " | " <<track2->GetNumberOfTPCClusters()<< endl;
	  
	  if( track1->GetNumberOfITSClusters() != track2->GetNumberOfITSClusters() ) cout<<"DIFFERENCE!: ";
	  cout<<" nITSclusters: "<<track1->GetNumberOfITSClusters()<< " | " <<track2->GetNumberOfITSClusters()<< endl;
	
/*

	for (Int_t idxCluster = 0; idxCluster < track1->GetNumberOfTPCClusters(); ++idxCluster){
	  Printf(" TEST: FlatTrack %d > FlatCluster %d has row %d", idxTrack, idxCluster, track1->GetTPCCluster(idxCluster).GetPadRow());
      }
  */    
	  }
     
      track1 = track1->GetNextTrack();
      track2 = track2->GetNextTrack();
	  
	  
	}
	  
	  
	  
	  hNTracks->Fill(flatEsd1->GetNumberOfTracks(),flatEsd2->GetNumberOfTracks());
	  
	  
      curr1=curr1+ flatEsd1->GetSize();
      curr2=curr2+ flatEsd2->GetSize();
      iEvent++;
    }

    delete[] buffer1;
    delete[] buffer2;
  }
  else {
    cout << "File could not be read" << endl;
  }
  
	TList histosList;
	histosList.Add(hNTracks);
  histosList.SaveAs(outputFilename);
  
  return;
}
