/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpTrigger.h,v 1.5 2006/05/24 13:58:27 ivana Exp $

/// \ingroup trigger
/// \class AliMpTrigger
/// \brief A trigger slat
///
//  Author: Laurent Aphecetche

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
class AliMpVSegmentation;
class TArrayI;

class AliMpTrigger : public TObject
{
public:
  AliMpTrigger();
  AliMpTrigger(const char* slatType, AliMp::PlaneType bendingOrNonBending);
  virtual ~AliMpTrigger();
  
  Bool_t AdoptLayer(AliMpSlat* slat);
    
  void GetAllLocalBoardNumbers(TArrayI& lbn) const;
  
  const char* GetID() const;
  
  const char* GetName() const;

  Double_t DX() const;
  Double_t DY() const;
  
  TVector2 Position() const;
  
  AliMpSlat* GetLayer(int layer) const;
  AliMpVSegmentation* GetLayerSegmentation(int layer) const;
  
  Int_t GetNofPadsX() const;
  
  Int_t GetMaxNofPadsY() const;
  
  /// Returns the number of layers.
  Int_t GetSize() const;
  
  void Print(Option_t* option="") const;

  AliMp::PlaneType PlaneType() const;
  
  TVector2 Dimensions() const;
  
private:
    
  Bool_t IsLayerValid(int layer) const;
  
  TString fId; ///< name of that slat
  AliMp::PlaneType fPlaneType; ///< bending or non-bending
  TObjArray fSlats; ///< virtual slat composing this trigger slat
  TObjArray fSlatSegmentations; ///< segmentations for virtual slats
  Int_t fMaxNofPadsY; ///< max number of pads in y direction
  Double_t fDX; ///< half-size in x (cm)
  Double_t fDY; ///< half-size in y (cm)
  
  ClassDef(AliMpTrigger,4) // Slat for trigger
};

#endif
