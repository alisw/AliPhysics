#ifndef ALIFLATESDVZERO_H
#define ALIFLATESDVZERO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/**
 * Flat structure representing a ESD VZERO 
 */

#include "Rtypes.h"
#include "AliVMisc.h"
#include "AliVVZERO.h"
#include "AliESDVZERO.h"

class AliFlatESDVZERO
{
 public:
  // -- Constructor / Destructors
  
  AliFlatESDVZERO();
 ~AliFlatESDVZERO() {}

  // constructor and method for reinitialisation of virtual table
  AliFlatESDVZERO( AliVConstructorReinitialisationFlag );
  void Reinitialize() const {} // no virtual table - do nothing

  void SetFromESDVZERO(const AliESDVZERO &v );
  void GetESDVZERO( AliESDVZERO &v ) const;

  static size_t GetSize() { return sizeof(AliFlatESDVZERO); }
  
private: 

  UInt_t  fBBtriggerV0A;     // bit mask for Beam-Beam trigger in V0A
  UInt_t  fBGtriggerV0A;     // bit mask for Beam-Gas trigger in V0A
  UInt_t  fBBtriggerV0C;     // bit mask for Beam-Beam trigger in V0C
  UInt_t  fBGtriggerV0C;     // bit mask for Beam-Gas trigger in V0C

  Float_t fMultiplicity[64]; //  multiplicity for each channel
  Float_t fAdc[64];          //  adc for each channel
  Float_t fTime[64];         //  time for each channel
  Float_t fWidth[64];        //  time width for each channel
  Bool_t  fBBFlag[64];       //  BB Flags from Online V0 Electronics
  Bool_t  fBGFlag[64];       //  BG Flags from Online V0 Electronics

  Float_t fV0ATime;          // Average time in V0A
  Float_t fV0CTime;          // Average time in V0C
  Float_t fV0ATimeError;     // Error in the average time in V0A
  Float_t fV0CTimeError;     // Error in the average time in V0C
  
  AliVVZERO::Decision fV0ADecision;     // V0A final decision based on average time of channels
  AliVVZERO::Decision fV0CDecision;     // V0C final decision based on average time of channels

  UShort_t fTriggerChargeA;  // Sum of the trigger (clock=10) charge on A side
  UShort_t fTriggerChargeC;  // Sum of the trigger (clock=10) charge on C side
  UShort_t fTriggerBits;     // V0 trigger bits as defined in the firmware  
};

inline AliFlatESDVZERO::AliFlatESDVZERO() :
   fBBtriggerV0A(0),
   fBGtriggerV0A(0),
   fBBtriggerV0C(0),
   fBGtriggerV0C(0),
   fV0ATime(-1024),
   fV0CTime(-1024),
   fV0ATimeError(0),
   fV0CTimeError(0),
   fV0ADecision(AliVVZERO::kV0Invalid),
   fV0CDecision(AliVVZERO::kV0Invalid),
   fTriggerChargeA(0),
   fTriggerChargeC(0),
   fTriggerBits(0)
{   
  // Default constructor 
  for(Int_t j=0; j<64; j++){ 
    fMultiplicity[j] = 0.0;   
    fAdc[j]   = 0.0;   
    fTime[j]  = 0.0; 
    fWidth[j] = 0.0; 
    fBBFlag[j]= kFALSE;
    fBGFlag[j]= kFALSE;  
  }
}

#pragma GCC diagnostic ignored "-Weffc++" 
inline AliFlatESDVZERO::AliFlatESDVZERO( AliVConstructorReinitialisationFlag ){}  // do nothing
#pragma GCC diagnostic warning "-Weffc++" 

inline void AliFlatESDVZERO::SetFromESDVZERO(const AliESDVZERO &v )
{
  // Set from ESD VZERO
  fBBtriggerV0A = v.GetBBTriggerV0A();
  fBGtriggerV0A = v.GetBGTriggerV0A();
  fBBtriggerV0C = v.GetBBTriggerV0C();
  fBGtriggerV0C = v.GetBGTriggerV0C();
  fV0ATime = v.GetV0ATime();
  fV0CTime = v.GetV0CTime();
  fV0ATimeError = v.GetV0ATimeError();
  fV0CTimeError = v.GetV0CTimeError();
  fV0ADecision = v.GetV0ADecision();
  fV0CDecision = v.GetV0CDecision();
  fTriggerChargeA = v.GetTriggerChargeA();
  fTriggerChargeC = v.GetTriggerChargeC();
  fTriggerBits = v.GetTriggerBits();

  for(Int_t j=0; j<64; j++){ 
    fMultiplicity[j] = v.GetMultiplicity(j);
    fAdc[j]   = v.GetAdc(j);
    fTime[j]  = v.GetTime(j);    
    fWidth[j] = v.GetWidth(j);
    fBBFlag[j]= v.GetBBFlag(j);
    fBGFlag[j]= v.GetBGFlag(j);  
  }
}

inline void AliFlatESDVZERO::GetESDVZERO( AliESDVZERO &v ) const
{
  // copy content to ESD VZERO object

  v.SetBBtriggerV0A( fBBtriggerV0A );
  v.SetBGtriggerV0A( fBGtriggerV0A );
  v.SetBBtriggerV0C( fBBtriggerV0C );
  v.SetBGtriggerV0C( fBGtriggerV0C );

  v.SetMultiplicity( fMultiplicity );
  v.SetADC( fAdc );
  v.SetTime( fTime );
  v.SetWidth( fWidth );
  v.SetBBFlag( fBBFlag );
  v.SetBGFlag( fBGFlag );

  v.SetV0ATime( fV0ATime );
  v.SetV0CTime( fV0CTime );
  v.SetV0ATimeError( fV0ATimeError );
  v.SetV0CTimeError( fV0CTimeError );

  v.SetV0ADecision( fV0ADecision );
  v.SetV0CDecision( fV0CDecision );

  v.SetTriggerChargeA( fTriggerChargeA );
  v.SetTriggerChargeC( fTriggerChargeC );
  v.SetTriggerBits( fTriggerBits );
}

#endif
