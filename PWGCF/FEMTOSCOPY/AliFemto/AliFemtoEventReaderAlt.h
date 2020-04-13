///
/// \file AliFemtoEventReaderAlt.h
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
  virtual ~AliFemtoEventReaderAlt();

  /// Randomly smear particles' momentum components with weights
  /// from gaussian distribution of mean 1.0 & the sigma provided
  void SetEnhanceSmearing(double sigma);
  double GetEnhanceSmearing() const
    { return fEnhanceSmearing; }

  /// randomly distribute particles around gaussian approximation of
  /// event shape
  void SetShouldDistribute(bool val=true);

protected:
  AliFemtoEventReaderAlt(const AliFemtoEventReaderAlt&);
  AliFemtoEventReaderAlt operator=(const AliFemtoEventReaderAlt&);

  virtual AliFemtoEvent* CopyAODtoFemtoEvent();

  virtual AliFemtoTrack* CopyAODtoFemtoTrack(AliAODTrack *src);

  virtual void CopyPIDtoFemtoTrack(AliAODTrack *, AliFemtoTrack *);

  //
  void RandomlyDistributeParticles(AliFemtoEvent &);

  TRandom3 *fRng;
  double fEnhanceSmearing;

  bool fDistributeMCParticles;
};


#endif
