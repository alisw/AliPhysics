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

void ReadFlatESD(const char* filename="outFlatESD.dat") {

  ifstream is(filename, std::ifstream::binary | std::ifstream::in);
  if (is){
    is.seekg (0, is.end);
    int length = is.tellg();
    is.seekg (0, is.beg);
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
      AliFlatESDEvent *flatEsd = reinterpret_cast<AliFlatESDEvent *>(curr);
			new (flatEsd) AliFlatESDEvent(1);

      cout<<"Reading event "<<iEvent<<":"<<endl;
      cout<<"vtx SPD: "<<(Bool_t) flatEsd->GetPrimaryVertexSPD()
	  <<" vtx tracks: "<<(Bool_t) flatEsd->GetPrimaryVertexTracks()	
	  <<" ntracks: "<<flatEsd->GetNumberOfTracks()
	  <<" nV0's: "<<flatEsd->GetNumberOfV0s()
	  <<endl;
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
