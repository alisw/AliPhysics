///
/// AliFemtoEventReaderAlt.h
///

#pragma once

#ifndef ALIFEMTOEVENTREADERALT_H
#define ALIFEMTOEVENTREADERALT_H

#include "AliFemtoEventReaderAODMultSelection.h"


class TRandom3;

/// \class AliFemtoEventReaderAlt
/// \brief Alternative Event Reader
///
/// AliFemtoEventReaderAODMultSelection
///
class AliFemtoEventReaderAlt : public AliFemtoEventReaderAODMultSelection {
public:

  AliFemtoEventReaderAlt();
  ~AliFemtoEventReaderAlt();

  void SetEnhanceSmearing(double n);
  double GetEnhanceSmearing() const
    { return fEnhanceSmearing; }

protected:
  AliFemtoEventReaderAlt(const AliFemtoEventReaderAlt&);
  AliFemtoEventReaderAlt operator=(const AliFemtoEventReaderAlt&);

  virtual AliFemtoTrack* CopyAODtoFemtoTrack(AliAODTrack *src);
  void CopyPIDtoFemtoTrack(AliAODTrack *, AliFemtoTrack *);

  TRandom3 *fRng;
  double fEnhanceSmearing;
};


#endif
