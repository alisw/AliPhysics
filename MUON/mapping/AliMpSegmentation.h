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
/// AliMpVSegmentation interface and to segmentations underlying
/// detectors (sector, slat, trigger).                                  \n
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

#ifndef ALI_MP_CATHOD_TYPE_H
#  include "AliMpCathodType.h"
#endif

#ifndef ALI_MP_V_SEGMENTATION_H
#  include "AliMpVSegmentation.h"
#endif

class AliMpDEStore;
class AliMpSegmentation;
class AliMpSlatMotifMap;
class AliMpDataStreams;
class AliMpSector;
class AliMpSlat;
class AliMpTrigger;
class TRootIOCtor;

class AliMpSegmentation : public  TObject {

  public:
    AliMpSegmentation(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpSegmentation();
    
    // static methods
    static AliMpSegmentation* Instance(Bool_t warn = true);
    static AliMpSegmentation* ReadData(const AliMpDataStreams& dataStreams, 
                                       Bool_t warn = true);

    // methods
    
    //
    // Access to segmentations
    //
    
    const AliMpVSegmentation* GetMpSegmentation(
                                 Int_t detElemId, AliMp::CathodType cath, 
			         Bool_t warn = true) const;

    const AliMpVSegmentation* GetMpSegmentationByElectronics(
                                 Int_t detElemId, Int_t elCardID, 
			         Bool_t warn = true) const;
    //
    // Access to sectors, slats, triggers
    //

    const AliMpSector*  GetSector(const AliMpVSegmentation* kSegmentation, 
                                  Bool_t warn = true) const;
    const AliMpSector*  GetSector(Int_t detElemId, AliMp::CathodType cath, 
                                  Bool_t warn = true) const;
    const AliMpSector*  GetSectorByElectronics(Int_t detElemId, Int_t elCardID, 
                                  Bool_t warn = true) const;
  
    const AliMpSlat*    GetSlat(const AliMpVSegmentation* kSegmentation, 
                                  Bool_t warn = true) const;
    const AliMpSlat*    GetSlat(Int_t detElemId, AliMp::CathodType cath, 
                                  Bool_t warn = true) const;
    const AliMpSlat*    GetSlatByElectronics(Int_t detElemId, Int_t elCardID, 
                                  Bool_t warn = true) const;
  
    const AliMpTrigger* GetTrigger(const AliMpVSegmentation* kSegmentation, 
                                  Bool_t warn = true) const;
    const AliMpTrigger* GetTrigger(Int_t detElemId, AliMp::CathodType cath, 
                                  Bool_t warn = true) const;
    const AliMpTrigger* GetTriggerByElectronics(Int_t detElemId, Int_t elCardID, 
                                  Bool_t warn = true) const;
  
  private:
    AliMpSegmentation(const AliMpDataStreams& dataStreams);
    /// Not implemented
    AliMpSegmentation(const AliMpSegmentation& rhs);
    /// Not implemented
    AliMpSegmentation& operator=(const AliMpSegmentation& rhs);

    AliMpVSegmentation* CreateMpSegmentation(
                              Int_t detElemId, AliMp::CathodType cath);
                              

    AliMpExMap* FillElCardsMap(Int_t detElemId);

    // static data members
    static AliMpSegmentation* fgInstance; ///< Singleton instance

    // data members
    const AliMpDataStreams&  fDataStreams;  //!< Data streams
    AliMpDEStore*      fDetElements;    ///< Detection element store
    AliMpStringObjMap  fMpSegmentations;///< Map of mapping segmentations to DE seg names
    AliMpExMap         fElCardsMap;     ///< Map of el. cards IDs to segmentations
    AliMpSlatMotifMap* fSlatMotifMap;   ///< Map of motif, motifTypes to avoid duplications and allow proper deletion
    
  ClassDef(AliMpSegmentation,3)  // The factory for building mapping segmentations
};

#endif //ALI_MP_SEGMENTATION_H















