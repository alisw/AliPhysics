////////////////////////////////////////////////
//  Reconstructed point class for set:ITS     //
////////////////////////////////////////////////


#include "AliITSRecPoint.h"
ClassImp(AliITSRecPoint)

AliITSRecPoint::AliITSRecPoint() {
    // default creator
    fTracks[0]=fTracks[1]=fTracks[2]=-3; 
    fX=fZ=fQ=fdEdX=0.;
    fSigmaX2=fSigmaZ2=0.;
   }
   

