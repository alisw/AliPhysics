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
#include "./AliFlatESDTrigger.h"
#include "./AliFlatESDV0.h"
#include "Riostream.h"
#endif   

void ReadFlatESD(const char* filename="outFlatESD.dat", Int_t verbose=0) {

  ifstream is(filename, std::ifstream::binary | std::ifstream::in);
  if (is){
    std::cout << "ifstream available"<<endl;
    is.seekg (0, is.end);
    int length = is.tellg();
    is.seekg (0, is.beg);
    std::cout << "length "<<length<<endl;
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
      cout<<endl<<"Reading event "<<iEvent<<":"<<endl;
 //   Printf("curr: %p \t endBuff: %p \t diff %p ", curr, endBuff, endBuff-curr);
      AliFlatESDEvent *flatEsd = reinterpret_cast<AliFlatESDEvent *>(curr);
			flatEsd->Reinitialize();

      
cout<<"vtx SPD: "<<(Bool_t) flatEsd->GetFlatPrimaryVertexSPD()
	  <<" vtx tracks: "<<(Bool_t) flatEsd->GetFlatPrimaryVertexTracks()	
	  <<" ntracks: "<<flatEsd->GetNumberOfTracks()
	  <<" nV0's: "<<flatEsd->GetNumberOfV0s()
	  <<endl;

// compare tracks
if(verbose){/*
	static const int nExt = 4;
	  AliFlatESDTrack *track = const_cast<AliFlatESDTrack*> ( flatEsd->GetTracks() );
	  //new (track)AliFlatESDTrack(1);
    for (Int_t idxTrack = 0; idxTrack < flatEsd->GetNumberOfTracks() && track; ++idxTrack) { 

	cout<<"track nr "<<idxTrack<<endl;

		const AliFlatExternalTrackParam* ext[nExt] ={
			
				track->GetFlatTrackParamRefitted(),
				track->GetFlatTrackParamIp(),
				track->GetFlatTrackParamTPCInner(),
				track->GetFlatTrackParamOp(),
		
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
      track = const_cast<AliFlatESDTrack*> (track->GetNextTrack());
	  
	  
	  }
	  
*/
	  
	  
	  // compare triggers
	  
	  /*
		cout<<"------------------\ntriggers\n------------------\n";
    AliFlatESDTrigger * trigger =const_cast<AliFlatESDTrigger*>(flatEsd->GetTriggerClasses() ) ;
    for( Int_t i = 0; i < flatEsd->GetNumberOfTriggerClasses() ; i++ ){
      cout<<"\nnew trigger\n";
			cout<<"AliFlatESDTrigger::GetSize"<<trigger->GetSize()<<endl;
			cout<<"AliFlatESDTrigger::GetTriggerIndex"<<trigger->GetTriggerIndex()<<endl;
			cout<< "AliFlatESDTrigger::GetTriggerClassName"<<trigger->GetTriggerClassName()<<endl;
			
      trigger= trigger->GetNextTriggerNonConst();
    }
	  **/
	  
	  // compare v0s

	  
	if(flatEsd->GetNumberOfV0s()  ){
		cout<<"------------------\nv0s\n------------------\n";
		
    AliFlatESDV0 * v0 = const_cast<AliFlatESDV0*>(flatEsd->GetV0s() ) ;
    for( Int_t i = 0; i < flatEsd->GetNumberOfV0s(); i++ ){
      cout<<"\nnew v0\n";
			cout<<"AliFlatESDV0::GetSize "<<v0->GetSize()<<endl; 
			cout<<"AliFlatESDV0::GetNegTrackID "<<v0->GetNegTrackID()<<endl ; 
			cout<<"AliFlatESDV0::GetPosTrackID "<<v0->GetPosTrackID()<<endl; 
			
      v0 = v0->GetNextV0NonConst();
    }
	}
	  
	  
	  
	  
	  
}

  //  Printf("curr: %p \t + %d = %p , diff:%p", curr, flatEsd->GetSize() ,curr+ flatEsd->GetSize(), endBuff-(curr+ flatEsd->GetSize())   );
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
