/**
 * >> Testing Macro to read FlatESDEvent from output file <<
 **
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli
 *
 * Usage:
 *  aliroot -b -l -q LoadLibs.C ReadFlatESD.C++
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

void ReadFlatESD(const char* filename="outFlatESD.dat", Bool_t verbose=kFALSE) {

  ifstream is(filename, std::ifstream::binary | std::ifstream::in);
  if (is){
    std::cout << "ifstream available"<<endl;
    is.seekg (0, is.end);
    int length = is.tellg();
    is.seekg (0, is.beg);
    std::cout << "length"<<length<<endl;
    char * buffer = new char [length];
    
    std::cout << "Reading " << length << " characters... ";
    
    is.read (buffer,length);
    if (is)
      std::cout << "all characters read successfully." << endl;
    else
      std::cout << "error: only " << is.gcount() << " could be read";
    is.close();
    
    // ...buffer contains the entire file...
    
    char *curr = buffer;
    char *endBuff = buffer+length;
    int iEvent = 0;
    while( curr < endBuff ){
      cout<<"Reading event "<<iEvent<<":"<<endl;
      AliFlatESDEvent *flatEsd = reinterpret_cast<AliFlatESDEvent *>(curr);
			new (flatEsd) AliFlatESDEvent(1);

      cout<<"vtx SPD: "<<(Bool_t) flatEsd->GetPrimaryVertexSPD()
	  <<" vtx tracks: "<<(Bool_t) flatEsd->GetPrimaryVertexTracks()	
	  <<" ntracks: "<<flatEsd->GetNumberOfTracks()
	  <<" nV0's: "<<flatEsd->GetNumberOfV0s()
	  <<endl;




// compare tracks
if(verbose){
	static const int nExt = 4;
	  AliFlatESDTrack *track = flatEsd->GetTracks();
    for (Int_t idxTrack = 0; idxTrack < flatEsd->GetNumberOfTracks() && track; ++idxTrack) { 

		AliFlatExternalTrackParam* ext[nExt] ={
			
				track->GetTrackParamRefitted(),
				track->GetTrackParamIp(),
				track->GetTrackParamTPCInner(),
				track->GetTrackParamOp(),
		
		};
	
     	//Printf("  TEST: FlatTrack1 %d > FlatExternalTrackParam1 > %p %p %p %p", idxTrack, exp11, exp21, exp31, exp41);
     	//Printf("  TEST: FlatTrack2 %d > FlatExternalTrackParam2 > %p %p %p %p", idxTrack, exp12, exp22, exp32, exp42);


	for(int iExt=0; iExt<nExt; ++iExt){
cout<<endl<<iExt<<endl;		
		if(!ext[iExt]){
		//	cout<<"DIFFERENCE!: ";
	 		cout<<" ext"<<iExt<<" not set"<<endl;
		}	


		cout<<" alpha"<<iExt<<" :"  << (ext[iExt] ? ext[iExt]->GetAlpha() : -9999) <<endl;
			

		
cout<<" GetX"<<iExt<<" :"  << (ext[iExt] ? ext[iExt]->GetX(): -9999) <<endl;


		
	cout<<" 1/pt"<<iExt<<" :"  <<  (ext[iExt] ? ext[iExt]->GetSigned1Pt(): -9999)  <<endl;
			

}
	
		cout<<" nTPCclusters: "<<track->GetNumberOfTPCClusters()<< endl;

	 	cout<<" nITSclusters: "<<track->GetNumberOfITSClusters()<< endl;


// compare clusters
	if(  track->GetNumberOfTPCClusters()){
		for (Int_t idxCluster = 0; idxCluster < track->GetNumberOfTPCClusters(); ++idxCluster){
			AliFlatTPCCluster * cl = track->GetTPCCluster(idxCluster);

			 	cout<<" clusterNr fX fY fZ fPadRow fSigmaY2 fSigmaZ2 fCharge fQMax" <<endl;
				cout<< idxCluster<<" "<< cl->GetX()<<" "<< cl->GetY()<<" "<< cl->GetZ()<<" "<< cl->GetPadRow()<<" "<< cl->GetSigmaY2()<<" "<< cl->GetSigmaZ2()<<" "<< cl->GetCharge()<<" "<< cl->GetQMax() <<endl;
				
		}
	  }
     
      track = track->GetNextTrack();
	  
	  
	  }
}











      curr=curr+ flatEsd->GetSize();
      iEvent++;
    }











    delete[] buffer;
  }
  else {
    cout << "File "<<filename<<" could not be read" << endl;
  }
  return;
}
