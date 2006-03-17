#ifndef ALI_MP_SEG_FACTORY_H
#define ALI_MP_SEG_FACTORY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 
// $MpId: AliMpSegFactory.h,v 1.5 2006/03/17 11:35:58 ivana Exp $ 

/// \ingroup management
/// \class AliMpSegFactory
/// \brief The factory for building mapping segmentations
//
/// The factory does not delete the created segmentation objects.
/// They have to be deleted in the client code.
/// As the same segmentation objects can be shared with more detection elements,
/// the class provides Clear() method for a safe deleting.
///
/// Authors: Ivana Hrivnacova, IPN Orsay

#ifndef ROOT_TObject
#  include <TObject.h>
#endif

#ifndef ALI_MP_STRING_OBJ_MAP_H
#  include "AliMpStringObjMap.h"
#endif

class AliMpExMap;
class AliMpVSegmentation;

class AliMpSegFactory : public  TObject {

  public:
    AliMpSegFactory();
    virtual ~AliMpSegFactory();
    
    // methods
    AliMpVSegmentation* CreateMpSegmentation(
                              Int_t detElemId, Int_t cath);

    AliMpVSegmentation* CreateMpSegmentationByElectronics(
                              Int_t detElemId, Int_t elCardID);

    void DeleteSegmentations();

  protected:
    AliMpSegFactory(const AliMpSegFactory& rhs);
    AliMpSegFactory& operator=(const AliMpSegFactory& rhs);

  private:
    AliMpExMap* FillMpMap(Int_t detElemId);
  
  private:

    AliMpStringObjMap  fMpSegmentations;// Map of mapping segmentations to DE names
    AliMpExMap*        fMpMap;          // Map of el. cards IDs to segmentations
      
  ClassDef(AliMpSegFactory,0)  // MUON Factory for Chambers and Segmentation
};

#endif //ALI_MP_SEG_FACTORY_H















