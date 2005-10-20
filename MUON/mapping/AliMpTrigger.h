/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$
// $MpId$

/// \ingroup trigger
/// \class AliMpTrigger
/// \brief A trigger slat
/// 
/// A trigger 'slat' object. It is to be viewed as a superposition of  
/// virtual layers of AliMpSlat objects. The need for more than one layer  
/// arise from the fact that a given local board deals with strips  
/// located in different detelem. So a given strip (pad) can have several  
/// "locations".
/// Author: Laurent Aphecetche

#ifndef ALI_MP_TRIGGER_H
#define ALI_MP_TRIGGER_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

#ifndef ROOT_TObjArray
#  include "TObjArray.h"
#endif

#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif

#ifndef ALI_MP_PLANE_TYPE
#  include "AliMpPlaneType.h"
#endif

class AliMpPCB;
class AliMpSlat;
class TArrayI;

class AliMpTrigger : public TObject
{
public:
  AliMpTrigger();
  AliMpTrigger(const char* slatType, AliMpPlaneType bendingOrNonBending);
  virtual ~AliMpTrigger();
  
  Bool_t AdoptLayer(AliMpSlat* slat);
    
  void GetAllLocalBoardNumbers(TArrayI& lbn) const;
  
  const char* GetID() const;
  
  const char* GetName() const;

  Double_t DX() const;
  Double_t DY() const;
  
  TVector2 Position() const;
  
  AliMpSlat* GetLayer(int layer) const;
  
  Int_t GetNofPadsX() const;
  
  Int_t GetMaxNofPadsY() const;
  
  /// Returns the number of layers.
  Int_t GetSize() const;
  
  void Print(Option_t* option="") const;

private:
    
  Bool_t IsLayerValid(int layer) const;
  
  TString fId;
  AliMpPlaneType fPlaneType;
  TObjArray fSlats;
  Int_t fMaxNofPadsY;
  Double_t fDX;
  Double_t fDY;
  
  ClassDef(AliMpTrigger,1) // Slat for trigger
};

#endif
