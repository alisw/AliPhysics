#ifndef ALIMUONSEGMENTATIONMANAGER_H
#define ALIMUONSEGMENTATIONMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup base
/// \class AliMUONSegmentationManager
/// \brief Segmentation manager

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TExMap
#include "TExMap.h"
#endif

#include "AliMpPlaneType.h"

class AliMpSlat;
class AliMpVSegmentation;

class AliMUONSegmentationManager : public TObject
{
public:
  AliMUONSegmentationManager();
  virtual ~AliMUONSegmentationManager();

  static Bool_t IsValidDetElemId(Int_t detElemId);

  static AliMpVSegmentation* Segmentation(Int_t detElemId,
                                          AliMpPlaneType planeType);

private:

  static const char* SlatType(Int_t detElemId);
  
  static AliMpSlat* ReadSlat(Int_t detElemId, AliMpPlaneType planeType);
  
  static bool ReadDetElemIdToSlatType();
    
  static AliMpVSegmentation* ReadSegmentation(Int_t detElemId,
                                              AliMpPlaneType planeType);
  
  static TExMap fgDetElemIdToSlatTypeMap; // map of int to TObjString
  
  static TExMap fgMap; // map of int to TPair<AliMpVSegmentation*, AliMpVSegmentation*>
  
  ClassDef(AliMUONSegmentationManager,1) // Holder for various segmentations
};

#endif
