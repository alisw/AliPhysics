#ifndef AliESDVZERO_H
#define AliESDVZERO_H

#include <TObject.h>

class AliESDVZERO : public TObject 
{
public:
  AliESDVZERO();
  AliESDVZERO(const AliESDVZERO&);
  AliESDVZERO(UInt_t BBtriggerV0A, UInt_t BGtriggerV0A,
	      UInt_t BBtriggerV0C, UInt_t BGtriggerV0C,
	      Short_t *Multiplicity);

  virtual ~AliESDVZERO() {};
  
  virtual void SetBBtriggerV0A(UInt_t BBtrigger) {fBBtriggerV0A=BBtrigger;}
  virtual void SetBGtriggerV0A(UInt_t BGtrigger) {fBGtriggerV0A=BGtrigger;}
  virtual void SetBBtriggerV0C(UInt_t BBtrigger) {fBBtriggerV0C=BBtrigger;}
  virtual void SetBGtriggerV0C(UInt_t BGtrigger) {fBGtriggerV0C=BGtrigger;}
  virtual void SetMultiplicity(Short_t Multiplicity[64])
    {for(Int_t i=0;i<64;i++) fMultiplicity[i]=Multiplicity[i];}

  // Getters  
  Short_t GetNbPMV0A();
  Short_t GetNbPMV0C();
  Int_t   GetMTotV0A();
  Int_t   GetMTotV0C(); 
  Int_t*  GetMRingV0A();
  Int_t*  GetMRingV0C();
  Int_t   GetMRingV0A(Int_t ring);
  Int_t   GetMRingV0C(Int_t ring);

  Int_t   GetMultiplicity(Int_t i);
  Bool_t  BBTriggerV0A(Int_t i);
  Bool_t  BGTriggerV0A(Int_t i);
  Bool_t  BBTriggerV0C(Int_t i);
  Bool_t  BGTriggerV0C(Int_t i);

  Bool_t  OutOfRange(Int_t i, const char *s, Int_t upper) const;
  AliESDVZERO &operator=(const AliESDVZERO& source);
    
protected:

  UInt_t  fBBtriggerV0A; // bit mask for Beam-Beam trigger in V0A
  UInt_t  fBGtriggerV0A; // bit mask for Beam-Gas trigger in V0A
  UInt_t  fBBtriggerV0C; // bit mask for Beam-Beam trigger in V0C
  UInt_t  fBGtriggerV0C; // bit mask for Beam-Gas trigger in V0C

  UShort_t fMultiplicity[64]; // multiplicity per PMT

  ClassDef(AliESDVZERO,3)
};

#endif
