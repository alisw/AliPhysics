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
//     Container class for AOD AD data
//     Author: Michal Broz
//     michal.broz@cern.ch
//-------------------------------------------------------------------------

#include "AliAODAD.h"
#include "AliLog.h"

ClassImp(AliAODAD)

//__________________________________________________________________________
AliAODAD::AliAODAD()
  :AliVAD(),
   fBBtriggerADA(0),
   fBGtriggerADA(0),
   fBBtriggerADC(0),
   fBGtriggerADC(0),
   fADATime(-1024),
   fADCTime(-1024),
   fADADecision(kADInvalid),
   fADCDecision(kADInvalid),
   fTriggerChargeA(0),
   fTriggerChargeC(0),
   fTriggerBits(0)
{   
   // Default constructor 
   for(Int_t j=0; j<16; j++){ 
      fMultiplicity[j] = 0.0;   
      fBBFlag[j]= kFALSE;
      fBGFlag[j]= kFALSE;  
   }
}

//__________________________________________________________________________
AliAODAD::AliAODAD(const AliAODAD &source)
  :AliVAD(source),
   fBBtriggerADA(source.fBBtriggerADA),
   fBGtriggerADA(source.fBGtriggerADA),
   fBBtriggerADC(source.fBBtriggerADC),
   fBGtriggerADC(source.fBGtriggerADC),
   fADATime(source.fADATime),
   fADCTime(source.fADCTime),
   fADADecision(source.fADADecision),
   fADCDecision(source.fADCDecision),
   fTriggerChargeA(source.fTriggerChargeA),
   fTriggerChargeC(source.fTriggerChargeC),
   fTriggerBits(source.fTriggerBits)
{   
   // Default constructor 
   for(Int_t j=0; j<16; j++) {
       fMultiplicity[j] = source.fMultiplicity[j];
       fBBFlag[j] = source.fBBFlag[j];
       fBGFlag[j] = source.fBGFlag[j];
   }
}

//__________________________________________________________________________
AliAODAD::AliAODAD(const AliVAD &source)
  :AliVAD(source),
   fBBtriggerADA(0),
   fBGtriggerADA(0),
   fBBtriggerADC(0),
   fBGtriggerADC(0),
   fADATime(source.GetADATime()),
   fADCTime(source.GetADCTime()),
   fADADecision(source.GetADADecision()),
   fADCDecision(source.GetADCDecision()),
   fTriggerChargeA(source.GetTriggerChargeA()),
   fTriggerChargeC(source.GetTriggerChargeC()),
   fTriggerBits(source.GetTriggerBits())
{   
   // Default constructor 
   for(Int_t j=0; j<16; j++) {
     fMultiplicity[j] = source.GetMultiplicity(j);
     fBBFlag[j] = source.GetBBFlag(j);
     fBGFlag[j] = source.GetBGFlag(j);
   }

   for(Int_t j=0; j<8; j++) {
     if (source.BBTriggerADA(j)) fBBtriggerADA |= (1 << j);
     if (source.BGTriggerADA(j)) fBGtriggerADA |= (1 << j);
     if (source.BBTriggerADC(j)) fBBtriggerADC |= (1 << j);
     if (source.BGTriggerADC(j)) fBGtriggerADC |= (1 << j);
   }
}

//__________________________________________________________________________
AliAODAD& AliAODAD::operator=(const AliAODAD& source)
{
  // Assignment operator
  //
  if(this==&source) return *this;
  AliVAD::operator=(source);
  // Assignment operator
  fBBtriggerADA=source.fBBtriggerADA;
  fBGtriggerADA=source.fBGtriggerADA;
  fBBtriggerADC=source.fBBtriggerADC;
  fBGtriggerADC=source.fBGtriggerADC;

  fADATime = source.fADATime;
  fADCTime = source.fADCTime;
  fADADecision = source.fADADecision;
  fADCDecision = source.fADCDecision;
  fTriggerChargeA = source.fTriggerChargeA;
  fTriggerChargeC = source.fTriggerChargeC;
  fTriggerBits = source.fTriggerBits;

   for(Int_t j=0; j<16; j++) {
       fMultiplicity[j] = source.fMultiplicity[j];
       fBBFlag[j] = source.fBBFlag[j];
       fBGFlag[j] = source.fBGFlag[j];
   }
  return *this;
}

//__________________________________________________________________________
AliAODAD& AliAODAD::operator=(const AliVAD& source)
{
  // Assignment operator
  // used in esd->aod filter
  if(this==&source) return *this;
  AliVAD::operator=(source);

  fADATime = source.GetADATime();
  fADCTime = source.GetADCTime();
  fADADecision = source.GetADADecision();
  fADCDecision = source.GetADCDecision();
  fTriggerChargeA = source.GetTriggerChargeA();
  fTriggerChargeC = source.GetTriggerChargeC();
  fTriggerBits = source.GetTriggerBits();

  for(Int_t j=0; j<16; j++) {
    fMultiplicity[j] = source.GetMultiplicity(j);
    fBBFlag[j] = source.GetBBFlag(j);
    fBGFlag[j] = source.GetBGFlag(j);
  }

  fBBtriggerADA = fBGtriggerADA = fBBtriggerADC = fBGtriggerADC = 0;
  for(Int_t j=0; j<8; j++) {
    if (source.BBTriggerADA(j)) fBBtriggerADA |= (1 << j);
    if (source.BGTriggerADA(j)) fBGtriggerADA |= (1 << j);
    if (source.BBTriggerADC(j)) fBBtriggerADC |= (1 << j);
    if (source.BGTriggerADC(j)) fBGtriggerADC |= (1 << j);
  }

  return *this;

}

//__________________________________________________________________________
Short_t AliAODAD::GetNbPMADA() const
{
  // Returns the number of
  // fired PM in ADA
  Short_t n=0;
  for(Int_t i=8;i<16;i++) 
    if (fMultiplicity[i]>0) n++;
  return n;
}

//__________________________________________________________________________
Short_t AliAODAD::GetNbPMADC() const
{
  // Returns the number of
  // fired PM in ADC
  Short_t n=0;
  for(Int_t i=0;i<8;i++) 
    if (fMultiplicity[i]>0) n++;
  return n;
}

//__________________________________________________________________________
Float_t AliAODAD::GetMTotADA() const
{
  // returns total multiplicity
  // in ADA
  Float_t mul=0.0;
  for(Int_t i=8;i<16;i++) 
    mul+=  fMultiplicity[i];
  return mul;
}

//__________________________________________________________________________
Float_t AliAODAD::GetMTotADC() const
{
  // returns total multiplicity
  // in ADC
  Float_t mul=0.0;
  for(Int_t i=0;i<8;i++) 
    mul+=  fMultiplicity[i];
  return mul;
}

//__________________________________________________________________________
Float_t AliAODAD::GetMultiplicity(Int_t i) const

{
  // returns multiplicity in a
  // given cell of AD
  if (OutOfRange(i, "AliAODAD::GetMultiplicity:",16)) return -1;
  return fMultiplicity[i];
}

//__________________________________________________________________________
Float_t AliAODAD::GetMultiplicityADA(Int_t i) const

{
  // returns multiplicity in a
  // given cell of ADA
  if (OutOfRange(i, "AliAODAD::GetMultiplicityADA:",8)) return -1;
  return fMultiplicity[8+i];
}

//__________________________________________________________________________
Float_t AliAODAD::GetMultiplicityADC(Int_t i) const

{
  // returns multiplicity in a
  // given cell of ADC
  if (OutOfRange(i, "AliAODAD::GetMultiplicityADC:",8)) return -1;
  return fMultiplicity[i];
}

//__________________________________________________________________________
Bool_t AliAODAD::BBTriggerADA(Int_t i) const
{
  // returns offline beam-beam flags in ADA
  // one bit per cell
  if (OutOfRange(i, "AliAODAD:::BBTriggerADA",8)) return kFALSE;
  UInt_t test = 1;
  return ( fBBtriggerADA & (test << i) ? kTRUE : kFALSE );
}

//__________________________________________________________________________
Bool_t AliAODAD::BGTriggerADA(Int_t i) const
{
  // returns offline beam-gas flags in ADA
  // one bit per cell
  if (OutOfRange(i, "AliAODAD:::BGTriggerADA",8)) return kFALSE;
  UInt_t test = 1;
  return ( fBGtriggerADA & (test << i) ? kTRUE : kFALSE );
}

//__________________________________________________________________________
Bool_t AliAODAD::BBTriggerADC(Int_t i) const
{
  // returns offline beam-beam flags in ADC
  // one bit per cell
  if (OutOfRange(i, "AliAODAD:::BBTriggerADC",8)) return kFALSE;
  UInt_t test = 1;
  return ( fBBtriggerADC & (test << i) ? kTRUE : kFALSE );
}

//__________________________________________________________________________
Bool_t AliAODAD::BGTriggerADC(Int_t i) const
{
  // returns offline beam-gasflags in ADC
  // one bit per cell
  if (OutOfRange(i, "AliAODAD:::BGTriggerADC",8)) return kFALSE;
  UInt_t test = 1;
  return ( fBGtriggerADC & (test << i) ? kTRUE : kFALSE );
}

//__________________________________________________________________________
Bool_t AliAODAD::GetBBFlag(Int_t i) const

{
  // returns online beam-beam flag in AD
  // one boolean per cell
  if (OutOfRange(i, "AliAODAD::GetBBFlag:",16)) return kFALSE;
  return fBBFlag[i];
}

//__________________________________________________________________________
Bool_t AliAODAD::GetBGFlag(Int_t i) const

{
  // returns online beam-gas flag in AD
  // one boolean per cell
  if (OutOfRange(i, "AliAODAD::GetBGFlag:",16)) return kFALSE;
  return fBGFlag[i];
}
