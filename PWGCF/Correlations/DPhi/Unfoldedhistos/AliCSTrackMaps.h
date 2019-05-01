#ifndef ALICSTRACKMAPS_H
#define ALICSTRACKMAPS_H

/// \file AliCSTrackMaps.h
/// \brief Track selection cuts support for the correlations studies tasks
///

#include <TNamed.h>
#include <TObjArray.h>
#include "AliAODTrack.h"

class AliVTrack;

/// \class AliCSTrackMaps
/// \brief Class which implements the track mapping needed for AOD events
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM, GSI
/// \date Oct 15, 2017

class AliCSTrackMaps : public TNamed {
public:
                                     AliCSTrackMaps();
                                     AliCSTrackMaps(const char *name, const char * title);
  virtual                           ~AliCSTrackMaps();

  virtual void                       NotifyEvent();

  static AliVTrack                  *GetOriginalTrack(AliAODTrack *trk);

private:
  static TObjArray                   fgAODTracksIdMap;                   ///< the map for AOD tracks ID

  /// Copy constructor
  /// Not allowed. Forced private.
  AliCSTrackMaps(const AliCSTrackMaps&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliCSTrackMaps& operator=(const AliCSTrackMaps&);

  /// \cond CLASSIMP
  ClassDef(AliCSTrackMaps,1);
  /// \endcond
};

/// \brief Returns the associated original track for a potentially constrained AOD track
/// \param trk the potentially constrained AOD track
/// \return the original track associated to the potentially constrained AOD track
inline AliVTrack *AliCSTrackMaps::GetOriginalTrack(AliAODTrack *trk)
{
  Int_t id = trk->GetID();
  if (id < 0) id = -1 -id;

  return (AliVTrack *) AliCSTrackMaps::fgAODTracksIdMap.At(id);
}



#endif /* ALICSTRACKMAPS_H */
