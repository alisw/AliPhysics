#include "AliESDVZERO.h"

ClassImp(AliESDVZERO)

AliESDVZERO::AliESDVZERO():TObject(),
   fMTotV0A(0),
   fMTotV0C(0),
   fNbPMV0A(0),
   fNbPMV0C(0)

{   
   // Default constructor 
   for(Int_t j=0; j<4; j++){ 
       fMRingV0A[j] = fMRingV0C[j] = 0;   
   }
}

AliESDVZERO::AliESDVZERO(const AliESDVZERO &o)
  :TObject(o),
   fMTotV0A(o.fMTotV0A),
   fMTotV0C(o.fMTotV0C),
   fNbPMV0A(o.fNbPMV0A),
   fNbPMV0C(o.fNbPMV0C)
{   
   // Default constructor 
   for(Int_t j=0; j<4; j++){ 
       fMRingV0A[j] = o.fMRingV0A[j];
       fMRingV0C[j] = o.fMRingV0C[j];   
   }
}

AliESDVZERO::AliESDVZERO(Int_t NbPMV0A, Int_t NbPMV0C, Int_t MTotV0A ,
             Int_t MTotV0C, Int_t* MRingV0A, Int_t* MRingV0C) 
  :TObject(),
   fMTotV0A(MTotV0A),
   fMTotV0C(MTotV0C),
   fNbPMV0A(NbPMV0A),
   fNbPMV0C(NbPMV0C)
{
   // Constructor
   for(Int_t j=0; j<4; j++){ 
       fMRingV0A[j] = MRingV0A[j];
       fMRingV0C[j] = MRingV0C[j]; }   
}

AliESDVZERO& AliESDVZERO::operator=(const AliESDVZERO& o)
{

  if(this==&o)return *this;
  TObject::operator=(o);
  // Assignment operator
  fMTotV0A=o.fMTotV0A;
  fMTotV0C=o.fMTotV0C;
  fNbPMV0A=o.fNbPMV0A;
  fNbPMV0C=o.fNbPMV0C;


  for(Int_t j=0; j<4; j++){ 
      fMRingV0A[j] = o.fMRingV0A[j];
      fMRingV0C[j] = o.fMRingV0C[j];   }

  return *this;
}

