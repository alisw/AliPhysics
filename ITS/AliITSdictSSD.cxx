
#include <stdio.h>
#include <TObjArray.h>

#include "AliITSsegmentationSSD.h"
#include "AliITSresponseSSD.h"
#include "AliITSsimulationSSD.h"
#include "AliITSdictSSD.h"
#include "AliITSdcsSSD.h"
#include "AliITS.h"
#include "AliRun.h"


void AliITSdictSSD::AddTrack(Int_t track) {
  // add track
    if (fTracks > 9) return;
    Int_t exist = 0,i;
    
    for(i=0; i<10; i++) 
     {
       if(track == fTrack[i]) 
        {
	 exist = 1;
	}
      
      } 
    if (!exist) 
     {
     fTrack[fTracks++] = track;     
     }

}


Int_t AliITSdictSSD::GetTrack(Int_t index) {
  // get track
    if(index > fTracks) return 0;
    else return fTrack[index];
}


//****************************************************************************
