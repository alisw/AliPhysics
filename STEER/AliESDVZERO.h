#ifndef AliESDVZERO_H
#define AliESDVZERO_H

#include <TObject.h>

class AliESDVZERO : public TObject 
{
public:
  AliESDVZERO();
  AliESDVZERO(const AliESDVZERO&);
  AliESDVZERO(Int_t NbPMV0A, Int_t NbPMV0C, Int_t MTotV0A, Int_t MTotV0C, 
              Int_t *MRingV0A, Int_t *MRingV0C);
  virtual ~AliESDVZERO() {};
  
// Setters
  virtual void  SetNbPMV0A(Int_t NbPMV0A)  {fNbPMV0A = NbPMV0A;}
  virtual void  SetNbPMV0C(Int_t NbPMV0C)  {fNbPMV0C = NbPMV0C;}
  virtual void  SetMTotV0A(Int_t MTotV0A)  {fMTotV0A = MTotV0A;}
  virtual void  SetMTotV0C(Int_t MTotV0C)  {fMTotV0C = MTotV0C;}					      
  virtual void  SetMRingV0A(Int_t MRingV0A[4]){for(Int_t j=0; j<4; j++){  
                                                fMRingV0A[j] = MRingV0A[j];} }
  virtual void  SetMRingV0C(Int_t MRingV0C[4]){for(Int_t j=0; j<4; j++){  
                                                fMRingV0C[j] = MRingV0C[j];} }
  
// Getters  
  Int_t GetNbPMV0A()  const {return fNbPMV0A;}
  Int_t GetNbPMV0C()  const {return fNbPMV0C;}
  Int_t GetMTotV0A()  const {return fMTotV0A;}
  Int_t GetMTotV0C()  const {return fMTotV0C;}
  Int_t* GetMRingV0A() const {return (int*) fMRingV0A;}
  Int_t* GetMRingV0C() const {return (int*) fMRingV0C;}

  AliESDVZERO &operator=(const AliESDVZERO& source);
    
protected:
  Int_t fNbPMV0A;     // Number of PMs fired in V0A - out of 32
  Int_t fNbPMV0C;     // Number of PMs fired in V0C - out of 32
  Int_t fMTotV0A;     // Total multiplicity in V0A
  Int_t fMTotV0C;     // Total multiplicity in V0C
  Int_t fMRingV0A[4]; // Multiplicity per ring in V0A - 4 rings
  Int_t fMRingV0C[4]; // Multiplicity per ring in V0C - 4 rings

  ClassDef(AliESDVZERO,1)
};

#endif
