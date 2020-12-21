///
/// \file AliFemtoXiSharedDaughterCut.h
///


/// \class AliFemtoXiSharedDaughterCut
///
/// \This macro is parallels, and is based on AliFemtoV0SharedDaughterCut,
/// \ which was written by Dominik Arominski
/// \As indicated by the class name, this is for Xi particles


#pragma once

#ifndef ALIFEMTOXISHAREDDAUGHTERCUT_H
#define ALIFEMTOXISHAREDDAUGHTERCUT_H

#include "AliFemtoXiTrackCut.h"
#include "AliFemtoV0SharedDaughterCut.h"

class AliFemtoXiSharedDaughterCut {
public:

  AliFemtoXiSharedDaughterCut();
  virtual ~AliFemtoXiSharedDaughterCut();
  virtual AliFemtoXiCollection AliFemtoXiSharedDaughterCutCollection(AliFemtoXiCollection *XiCollection, AliFemtoXiTrackCut *pCut);



#ifdef __ROOT__
  ClassDef(AliFemtoXiSharedDaughterCut, 1);
#endif

};

#endif
