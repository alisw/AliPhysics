#ifndef ALI_MP_SEGMENTATION_H
#define ALI_MP_SEGMENTATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 
// $MpId: AliMpSegmentation.h,v 1.7 2006/05/24 13:58:16 ivana Exp $ 

/// \ingroup management
/// \class AliMpSegmentation
/// \brief Singleton container class for mapping segmentations
///
/// It provides access to mapping segmentations based on the
/// AliMpVSegmentation interface.                                        \n
/// Mapping segmentations for all detection elements
/// are created at the first call to AliMpSegmentation::Instance().
/// The class is a singleton, it has all constructors
/// private, except for the special constructor for Root I/O.
///
/// \author Ivana Hrivnacova, IPN Orsay; Laurent Aphecetche, SUBATECH

#ifndef ROOT_TObject
#  include <TObject.h>
#endif

#ifndef ALI_MP_STRING_OBJ_MAP_H
#  include "AliMpStringObjMap.h"
#endif

#ifndef ALI_MP_EX_MAP_H
#  include "AliMpExMap.h"
#endif

class AliMpVSegmentation;
class AliMpSegmentation;

class TRootIOCtor;

class AliMpSegmentation : public  TObject {

  public:
    AliMpSegmentation(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpSegmentation();
    
    // static methods
    static AliMpSegmentation* Instance();

    // methods
    const AliMpVSegmentation* GetMpSegmentation(
                                 Int_t detElemId, Int_t cath, 
			         Bool_t warn = true) const;

    const AliMpVSegmentation* GetMpSegmentationByElectronics(
                                 Int_t detElemId, Int_t elCardID, 
			         Bool_t warn = true) const;
    
  private:
    AliMpSegmentation();
    AliMpSegmentation(const AliMpSegmentation& rhs);
    AliMpSegmentation& operator=(const AliMpSegmentation& rhs);

    AliMpVSegmentation* CreateMpSegmentation(
                              Int_t detElemId, Int_t cath);

    AliMpExMap* FillElCardsMap(Int_t detElemId);

    // static data members
    static AliMpSegmentation* fgInstance; ///< Singleton instance

    // data members
    AliMpStringObjMap  fMpSegmentations;///< Map of mapping segmentations to DE seg names
    AliMpExMap         fElCardsMap;     ///< Map of el. cards IDs to segmentations
      
    
  ClassDef(AliMpSegmentation,1)  // The factory for building mapping segmentations
};

#endif //ALI_MP_SEGMENTATION_H















