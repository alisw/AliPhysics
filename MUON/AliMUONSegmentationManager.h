#ifndef ALIMUONSEGMENTATIONMANAGER_H
#define ALIMUONSEGMENTATIONMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup base
/// \class AliMUONSegmentationManager
/// \brief Segmentation manager

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ALI_MP_EX_MAP_H
#  include "AliMpExMap.h"
#endif

#ifndef ALI_MP_PLANE_TYPE
#  include "AliMpPlaneType.h"
#endif

#ifndef ALI_MP_STATION_TYPE
#  include "AliMpStationType.h"
#endif

class AliMpSlat;
class AliMpTriggerSegmentation;
class AliMpVSegmentation;
class TList;

class AliMUONSegmentationManager : public TObject
{
public:
  AliMUONSegmentationManager();
  virtual ~AliMUONSegmentationManager();

  static Bool_t IsValidDetElemId(Int_t detElemId);

  static AliMpVSegmentation* Segmentation(Int_t detElemId,
                                          AliMpPlaneType planeType);
  
  static TList* SegmentationList(Int_t localBoardNumber);

  static const char* DetElemName(Int_t detElemId);
  
  static AliMpStationType StationType(Int_t detElemId);
  
private:

    static void FillLocalBoardMap(AliMpTriggerSegmentation* seg);
  
  static const char* SlatType(Int_t detElemId);
  
  static bool ReadDetElemIdToName(AliMpStationType stationType);
    
  static AliMpVSegmentation* ReadSegmentation(Int_t detElemId,
                                              AliMpPlaneType planeType);
  
  static AliMpExMap fgDetElemIdToNameMap; // map of int to TObjString
  
  static AliMpExMap fgMap; // map of int to TPair<AliMpVSegmentation*, AliMpVSegmentation*>
  
  static AliMpExMap fgLocalBoardMap; // map of int to TList* of AliMpVSegmentation*
  
  ClassDef(AliMUONSegmentationManager,1) // Holder for various segmentations
};

#endif
