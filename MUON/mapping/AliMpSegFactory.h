#ifndef ALI_MP_SEG_FACTORY_H
#define ALI_MP_SEG_FACTORY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 
// $MpId: AliMpSegFactory.h,v 1.7 2006/05/24 13:58:16 ivana Exp $ 

/// \ingroup management
/// \class AliMpSegFactory
/// \brief The factory for building mapping segmentations
//
/// The factory does not delete the created segmentation objects.
/// They have to be deleted in the client code.
/// As the same segmentation objects can be shared with more detection elements,
/// the class provides Clear() method for a safe deleting.
///
/// \author Ivana Hrivnacova, IPN Orsay

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

    static Int_t NumberOfInstances() { return fgNumberOfInstances; }
    
  private:
    AliMpSegFactory(const AliMpSegFactory& rhs);
    AliMpSegFactory& operator=(const AliMpSegFactory& rhs);

    AliMpExMap* FillMpMap(Int_t detElemId);

    AliMpStringObjMap  fMpSegmentations;///< Map of mapping segmentations to DE names
    AliMpExMap*        fMpMap;          ///< Map of el. cards IDs to segmentations
      
    static Int_t fgNumberOfInstances; ///< number of AliMpSegFactory objects...
    
  ClassDef(AliMpSegFactory,0)  // The factory for building mapping segmentations
};

#endif //ALI_MP_SEG_FACTORY_H















