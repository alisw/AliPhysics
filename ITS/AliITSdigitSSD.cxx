#include "AliITSdigitSSD.h"

ClassImp(AliITSdigitSSD);

AliITSdigitSSD::AliITSdigitSSD 
    (Int_t *tracks, Int_t *digits, Int_t strNo, Int_t s, Bool_t p)
    :AliITSdigit(tracks,digits) {
 
    fSignal = s;
    fStripNumber = strNo;
    fSide = p;
}

