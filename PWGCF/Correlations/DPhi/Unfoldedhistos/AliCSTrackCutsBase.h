#ifndef ALICSTRACKCUTSBASE_H
#define ALICSTRACKCUTSBASE_H

/// \file AliCSTrackCutsBase.h
/// \brief Base track cuts support for the correlations studies tasks
///

#include "AliCSAnalysisCutsBase.h"

class AliVTrack;

/// \class AliCSTrackCutsBase
/// \brief Track cuts base class for correlation studies analysis
///
/// Abstract class support for different kind of tracks cuts.
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date March 05, 2017

class AliCSTrackCutsBase : public AliCSAnalysisCutsBase {
public:

                                AliCSTrackCutsBase();
                                AliCSTrackCutsBase(Int_t nCuts, Int_t nParams, const char *name="CS TrackCuts", const char * title="CS TrackCuts");
  virtual                      ~AliCSTrackCutsBase();

                                /// Initializes the cuts
                                /// Pure virtual function
                                /// \param name the name to assign to the histograms list
  virtual void                  InitCuts(const char *name) = 0;
                                /// Processes a potential change in the run number
                                /// Pure virtual function
  virtual void                  NotifyRun() = 0;

                                /// Is the track accepted by the set of cuts?
                                /// Pure virtual function
                                /// \param trk the track to be accepted or rejected
                                /// \return kTRUE if the track is accepted kFALSE otherwise
  virtual Bool_t                IsTrackAccepted(AliVTrack *trk) = 0;
                                /// Is the corresponding true track accepted by the set of cuts?
                                /// Pure virtual function
                                /// \param trk the track whose corresponding true track is to be accepted or rejected
                                /// \return kTRUE if the associated true track is accepted kFALSE otherwise
  virtual Bool_t                IsTrueTrackAccepted(AliVTrack *trk) = 0;
                                /// Is the true track accepted by the set of cuts?
                                /// Pure virtual function
                                /// \param itrk the stack track index of the true track to be accepted or rejected
                                /// \return kTRUE if the true track is accepted kFALSE otherwise
  virtual Bool_t                IsTrueTrackAccepted(Int_t itrk) = 0;


private:
  /// Copy constructor
  /// Not allowed. Forced private.
  AliCSTrackCutsBase(const AliCSTrackCutsBase&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliCSTrackCutsBase& operator=(const AliCSTrackCutsBase&);

  /// \cond CLASSIMP
  ClassDef(AliCSTrackCutsBase,1);
  /// \endcond
};

#endif /* ALICSTRACKCUTSBASE_H */
