/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
*/

#include "AliITSmodule.h"
#include "AliRun.h"

ClassImp(AliITSmodule)

//_______________________________________________________________________
//
// Impementation of class AliITSmodule
//
// created by: A. Bouchm, W. Peryt, S. Radomski, P. Skowronski
//             R. Barbers, B. Batyunia, B. S. Nilsen
// ver 1.0     CERN 16.09.1999
//_______________________________________________________________________
//

//________________________________________________________________________
// 
// Constructors and deconstructor
//________________________________________________________________________
//


AliITSmodule::AliITSmodule() {

    fIndex  = 0;
    fHitsM  = new TObjArray();
    fNhitsM = 0;
    fITS    = (AliITS*)(gAlice->GetDetector("ITS"));
}


//_________________________________________________________________________

AliITSmodule::AliITSmodule(Int_t index) {

    fIndex  = index;
    fHitsM  = new TObjArray();
    fNhitsM = 0;
    fITS    = (AliITS*)(gAlice->GetDetector("ITS"));
}


//__________________________________________________________________________

AliITSmodule::~AliITSmodule() {
    if(fHitsM) delete fHitsM;
    fNhitsM = 0;
}

//_________________________________________________________________________
// 
// Hits menagement
//__________________________________________________________________________


Int_t AliITSmodule::AddHit(AliITShit* hit ) {
    fHitsM->AddLast(hit);
    fNhitsM = fHitsM->GetEntriesFast();
    return fNhitsM;
}

//___________________________________________________________________________
void AliITSmodule::Streamer(TBuffer &R__b){
   // Stream an object of class AliITSmodule.

    return;
// This class is not to be written out to any file.
//   if (R__b.IsReading()) {
//      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
//      TObject::Streamer(R__b);
//      R__b >> fITS;
//      R__b >> fIndex;
//      R__b >> fHitsM;
//      R__b >> fNhitsM;
//      R__b >> fIDigits;
//      R__b >> fNdigits;
//      R__b >> fIPoints;
//      R__b >> fNpoints;
//   } else {
//      R__b.WriteVersion(AliITSmodule::IsA());
//      TObject::Streamer(R__b);
//      R__b << fITS;
//      R__b << fIndex;
//      R__b << fHitsM;
//      R__b << fNhitsM;
//      R__b << fIDigits;
//      R__b << fNdigits;
//      R__b << fIPoints;
//      R__b << fNpoints;
//   }
}
//______________________________________________________________________
