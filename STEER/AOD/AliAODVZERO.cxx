/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     Container class for AOD VZERO data
//     Author: Cvetan Cheshkov
//     cvetan.cheshkov@cern.ch 2/02/2011
//-------------------------------------------------------------------------

#include "AliAODVZERO.h"
#include "AliLog.h"

ClassImp(AliAODVZERO)

//__________________________________________________________________________
AliAODVZERO::AliAODVZERO()
  :AliVVZERO(),
   fBBtriggerV0A(0),
   fBGtriggerV0A(0),
   fBBtriggerV0C(0),
   fBGtriggerV0C(0),
   fV0ATime(-1024),
   fV0CTime(-1024),
   fV0ADecision(kV0Invalid),
   fV0CDecision(kV0Invalid),
   fTriggerChargeA(0),
   fTriggerChargeC(0),
   fTriggerBits(0)
{   
   // Default constructor 
   for(Int_t j=0; j<64; j++){ 
      fMultiplicity[j] = 0.0;   
      fBBFlag[j]= kFALSE;
      fBGFlag[j]= kFALSE;  
   }
}

//__________________________________________________________________________
AliAODVZERO::AliAODVZERO(const AliAODVZERO &source)
  :AliVVZERO(source),
   fBBtriggerV0A(source.fBBtriggerV0A),
   fBGtriggerV0A(source.fBGtriggerV0A),
   fBBtriggerV0C(source.fBBtriggerV0C),
   fBGtriggerV0C(source.fBGtriggerV0C),
   fV0ATime(source.fV0ATime),
   fV0CTime(source.fV0CTime),
   fV0ADecision(source.fV0ADecision),
   fV0CDecision(source.fV0CDecision),
   fTriggerChargeA(source.fTriggerChargeA),
   fTriggerChargeC(source.fTriggerChargeC),
   fTriggerBits(source.fTriggerBits)
{   
   // Default constructor 
   for(Int_t j=0; j<64; j++) {
       fMultiplicity[j] = source.fMultiplicity[j];
       fBBFlag[j] = source.fBBFlag[j];
       fBGFlag[j] = source.fBGFlag[j];
   }
}

//__________________________________________________________________________
AliAODVZERO::AliAODVZERO(const AliVVZERO &source)
  :AliVVZERO(source),
   fBBtriggerV0A(0),
   fBGtriggerV0A(0),
   fBBtriggerV0C(0),
   fBGtriggerV0C(0),
   fV0ATime(source.GetV0ATime()),
   fV0CTime(source.GetV0CTime()),
   fV0ADecision(source.GetV0ADecision()),
   fV0CDecision(source.GetV0CDecision()),
   fTriggerChargeA(source.GetTriggerChargeA()),
   fTriggerChargeC(source.GetTriggerChargeC()),
   fTriggerBits(source.GetTriggerBits())
{   
   // Default constructor 
   for(Int_t j=0; j<64; j++) {
     fMultiplicity[j] = source.GetMultiplicity(j);
     fBBFlag[j] = source.GetBBFlag(j);
     fBGFlag[j] = source.GetBGFlag(j);
   }

   for(Int_t j=0; j<32; j++) {
     if (source.BBTriggerV0A(j)) fBBtriggerV0A |= (1 << j);
     if (source.BGTriggerV0A(j)) fBGtriggerV0A |= (1 << j);
     if (source.BBTriggerV0C(j)) fBBtriggerV0C |= (1 << j);
     if (source.BGTriggerV0C(j)) fBGtriggerV0C |= (1 << j);
   }
}

//__________________________________________________________________________
AliAODVZERO& AliAODVZERO::operator=(const AliAODVZERO& source)
{
  // Assignment operator
  //
  if(this==&source) return *this;
  AliVVZERO::operator=(source);
  // Assignment operator
  fBBtriggerV0A=source.fBBtriggerV0A;
  fBGtriggerV0A=source.fBGtriggerV0A;
  fBBtriggerV0C=source.fBBtriggerV0C;
  fBGtriggerV0C=source.fBGtriggerV0C;

  fV0ATime = source.fV0ATime;
  fV0CTime = source.fV0CTime;
  fV0ADecision = source.fV0ADecision;
  fV0CDecision = source.fV0CDecision;
  fTriggerChargeA = source.fTriggerChargeA;
  fTriggerChargeC = source.fTriggerChargeC;
  fTriggerBits = source.fTriggerBits;

   for(Int_t j=0; j<64; j++) {
       fMultiplicity[j] = source.fMultiplicity[j];
       fBBFlag[j] = source.fBBFlag[j];
       fBGFlag[j] = source.fBGFlag[j];
   }
  return *this;
}

//__________________________________________________________________________
AliAODVZERO& AliAODVZERO::operator=(const AliVVZERO& source)
{
  // Assignment operator
  // used in esd->aod filter
  if(this==&source) return *this;
  AliVVZERO::operator=(source);

  fV0ATime = source.GetV0ATime();
  fV0CTime = source.GetV0CTime();
  fV0ADecision = source.GetV0ADecision();
  fV0CDecision = source.GetV0CDecision();
  fTriggerChargeA = source.GetTriggerChargeA();
  fTriggerChargeC = source.GetTriggerChargeC();
  fTriggerBits = source.GetTriggerBits();

  for(Int_t j=0; j<64; j++) {
    fMultiplicity[j] = source.GetMultiplicity(j);
    fBBFlag[j] = source.GetBBFlag(j);
    fBGFlag[j] = source.GetBGFlag(j);
  }

  fBBtriggerV0A = fBGtriggerV0A = fBBtriggerV0C = fBGtriggerV0C = 0;
  for(Int_t j=0; j<32; j++) {
    if (source.BBTriggerV0A(j)) fBBtriggerV0A |= (1 << j);
    if (source.BGTriggerV0A(j)) fBGtriggerV0A |= (1 << j);
    if (source.BBTriggerV0C(j)) fBBtriggerV0C |= (1 << j);
    if (source.BGTriggerV0C(j)) fBGtriggerV0C |= (1 << j);
  }

  return *this;

}

//__________________________________________________________________________
Short_t AliAODVZERO::GetNbPMV0A() const
{
  // Returns the number of
  // fired PM in V0A
  Short_t n=0;
  for(Int_t i=32;i<64;i++) 
    if (fMultiplicity[i]>0) n++;
  return n;
}

//__________________________________________________________________________
Short_t AliAODVZERO::GetNbPMV0C() const
{
  // Returns the number of
  // fired PM in V0C
  Short_t n=0;
  for(Int_t i=0;i<32;i++) 
    if (fMultiplicity[i]>0) n++;
  return n;
}

//__________________________________________________________________________
Float_t AliAODVZERO::GetMTotV0A() const
{
  // returns total multiplicity
  // in V0A
  Float_t mul=0.0;
  for(Int_t i=32;i<64;i++) 
    mul+=  fMultiplicity[i];
  return mul;
}

//__________________________________________________________________________
Float_t AliAODVZERO::GetMTotV0C() const
{
  // returns total multiplicity
  // in V0C
  Float_t mul=0.0;
  for(Int_t i=0;i<32;i++) 
    mul+=  fMultiplicity[i];
  return mul;
}

//__________________________________________________________________________
Float_t AliAODVZERO::GetMRingV0A(Int_t ring) const
{ 
  // returns multiplicity in a
  // given ring of V0A
  if (OutOfRange(ring, "AliAODVZERO:::GetMRingV0A",4)) return -1;
  Float_t mul =0.0;

  if (ring == 0) for(Int_t i=32;i<40;i++) mul +=  fMultiplicity[i];
  if (ring == 1) for(Int_t i=40;i<48;i++) mul +=  fMultiplicity[i];
  if (ring == 2) for(Int_t i=48;i<56;i++) mul +=  fMultiplicity[i];
  if (ring == 3) for(Int_t i=56;i<64;i++) mul +=  fMultiplicity[i];
  return mul ;
}

//__________________________________________________________________________
Float_t AliAODVZERO::GetMRingV0C(Int_t ring) const
{ 
  // returns multiplicity in a
  // given ring of V0C
  if (OutOfRange(ring, "AliAODVZERO:::GetMRingV0C",4)) return -1;
  Float_t mul =0.0;

  if (ring == 0) for(Int_t i=0;i<8;i++)   mul +=  fMultiplicity[i];
  if (ring == 1) for(Int_t i=8;i<16;i++)  mul +=  fMultiplicity[i];
  if (ring == 2) for(Int_t i=16;i<24;i++) mul +=  fMultiplicity[i];
  if (ring == 3) for(Int_t i=24;i<32;i++) mul +=  fMultiplicity[i];
  return mul ;
}

//__________________________________________________________________________
Float_t AliAODVZERO::GetMultiplicity(Int_t i) const

{
  // returns multiplicity in a
  // given cell of V0
  if (OutOfRange(i, "AliAODVZERO::GetMultiplicity:",64)) return -1;
  return fMultiplicity[i];
}

//__________________________________________________________________________
Float_t AliAODVZERO::GetMultiplicityV0A(Int_t i) const

{
  // returns multiplicity in a
  // given cell of V0A
  if (OutOfRange(i, "AliAODVZERO::GetMultiplicityV0A:",32)) return -1;
  return fMultiplicity[32+i];
}

//__________________________________________________________________________
Float_t AliAODVZERO::GetMultiplicityV0C(Int_t i) const

{
  // returns multiplicity in a
  // given cell of V0C
  if (OutOfRange(i, "AliAODVZERO::GetMultiplicityV0C:",32)) return -1;
  return fMultiplicity[i];
}

//__________________________________________________________________________
Bool_t AliAODVZERO::BBTriggerV0A(Int_t i) const
{
  // returns offline beam-beam flags in V0A
  // one bit per cell
  if (OutOfRange(i, "AliAODVZERO:::BBTriggerV0A",32)) return kFALSE;
  UInt_t test = 1;
  return ( fBBtriggerV0A & (test << i) ? kTRUE : kFALSE );
}

//__________________________________________________________________________
Bool_t AliAODVZERO::BGTriggerV0A(Int_t i) const
{
  // returns offline beam-gas flags in V0A
  // one bit per cell
  if (OutOfRange(i, "AliAODVZERO:::BGTriggerV0A",32)) return kFALSE;
  UInt_t test = 1;
  return ( fBGtriggerV0A & (test << i) ? kTRUE : kFALSE );
}

//__________________________________________________________________________
Bool_t AliAODVZERO::BBTriggerV0C(Int_t i) const
{
  // returns offline beam-beam flags in V0C
  // one bit per cell
  if (OutOfRange(i, "AliAODVZERO:::BBTriggerV0C",32)) return kFALSE;
  UInt_t test = 1;
  return ( fBBtriggerV0C & (test << i) ? kTRUE : kFALSE );
}

//__________________________________________________________________________
Bool_t AliAODVZERO::BGTriggerV0C(Int_t i) const
{
  // returns offline beam-gasflags in V0C
  // one bit per cell
  if (OutOfRange(i, "AliAODVZERO:::BGTriggerV0C",32)) return kFALSE;
  UInt_t test = 1;
  return ( fBGtriggerV0C & (test << i) ? kTRUE : kFALSE );
}

//__________________________________________________________________________
Bool_t AliAODVZERO::GetBBFlag(Int_t i) const

{
  // returns online beam-beam flag in V0
  // one boolean per cell
  if (OutOfRange(i, "AliAODVZERO::GetBBFlag:",64)) return kFALSE;
  return fBBFlag[i];
}

//__________________________________________________________________________
Bool_t AliAODVZERO::GetBGFlag(Int_t i) const

{
  // returns online beam-gas flag in V0
  // one boolean per cell
  if (OutOfRange(i, "AliAODVZERO::GetBGFlag:",64)) return kFALSE;
  return fBGFlag[i];
}
