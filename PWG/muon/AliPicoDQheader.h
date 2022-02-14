#ifndef ALIPICODQHEADER_H
#define ALIPICODQHEADER_H

#include <TNamed.h>

class AliPicoDQheader : public TNamed {

 public:

  AliPicoDQheader(const TString s="AliPicoDQheader");
  AliPicoDQheader(const AliPicoDQheader &a);
  AliPicoDQheader &operator=(const AliPicoDQheader &a);
  virtual ~AliPicoDQheader();
//=============================================================================

  UInt_t  PSmask()            const { return fPSmask;            }
  UInt_t  TriggerInputs()     const { return fTriggerInputs;     }
  TString FiredTriggerClass() const { return fFiredTriggerClass; }
//=============================================================================

  Double_t Vz() const { return fVz; }
  Float_t MultSPDtrkls()  const { return fMultSPDtrkls; }
  Float_t MultV0M()  const { return fMultV0M; }
  Float_t MultV0C()  const { return fMultV0C; }
  Int_t   NSPDtracklets() const { return fNSPDtrkls; }
//=============================================================================

  void Reset();
  void SetEventInfo(UInt_t   wMask,
                    UInt_t   wTrgIn,
                    TString  sTrgCls,
                    Double_t dVz,
                    Float_t  dMultSPDtrkls,
                    Float_t  dMultV0M,
                    Float_t  dMultV0C,
                    Int_t    nTrkls);
//=============================================================================

 private :

  UInt_t  fPSmask;            // Physics selection mask
  UInt_t  fTriggerInputs;     // Trigger inputs (ID)
  TString fFiredTriggerClass; // Trigger class

  Double_t fVz;     //
  Float_t fMultSPDtrkls;    //
  Float_t fMultV0M;    //
  Float_t fMultV0C;    //
  Int_t fNSPDtrkls; //
//=============================================================================

 ClassDef(AliPicoDQheader, 3)
};

#endif

