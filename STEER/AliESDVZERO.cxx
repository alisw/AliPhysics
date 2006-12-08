#include "AliESDVZERO.h"

ClassImp(AliESDVZERO)

AliESDVZERO::AliESDVZERO():TObject(),
   fNbPMV0A(0),
   fNbPMV0C(0),
   fMTotV0A(0),
   fMTotV0C(0)
{   
   // Default constructor 
   for(Int_t j=0; j<4; j++){ 
       fMRingV0A[j] = 0;
       fMRingV0C[j] = 0;   }
}

AliESDVZERO::AliESDVZERO(const AliESDVZERO &o):TObject(o),
   fNbPMV0A(o.fNbPMV0A),
   fNbPMV0C(o.fNbPMV0C),
   fMTotV0A(o.fMTotV0A),
   fMTotV0C(o.fMTotV0C)
{   
   // Default constructor 
   for(Int_t j=0; j<4; j++){ 
       fMRingV0A[j] = o.fMRingV0A[j];
       fMRingV0C[j] = o.fMRingV0C[j];   }
}

AliESDVZERO::AliESDVZERO(Int_t NbPMV0A, Int_t NbPMV0C, Int_t MTotV0A ,
             Int_t MTotV0C, Int_t* MRingV0A, Int_t* MRingV0C) :TObject(),
   fNbPMV0A(NbPMV0A),
   fNbPMV0C(NbPMV0C),
   fMTotV0A(MTotV0A),
   fMTotV0C(MTotV0C)

{
   // Constructor
   for(Int_t j=0; j<4; j++){ 
       fMRingV0A[j] = MRingV0A[j];
       fMRingV0C[j] = MRingV0C[j]; }   
}

AliESDVZERO& AliESDVZERO::operator=(const AliESDVZERO& o)
{
  // Assignment operator
  fNbPMV0A=o.fNbPMV0A;
  fNbPMV0C=o.fNbPMV0C;
  fMTotV0A=o.fMTotV0A;
  fMTotV0C=o.fMTotV0C;

  for(Int_t j=0; j<4; j++){ 
      fMRingV0A[j] = o.fMRingV0A[j];
      fMRingV0C[j] = o.fMRingV0C[j];   }

  return *this;
}

