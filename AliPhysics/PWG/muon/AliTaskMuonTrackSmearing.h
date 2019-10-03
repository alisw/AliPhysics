#ifndef ALITASKMUONTRACKSMEARING_H
#define ALITASKMUONTRACKSMEARING_H

/* $Id$ */ 

///
/// \class AliTaskMuonTrackSmearing
/// \brief Task to smear the muon track parameter according to resolution
///
/// \author Diego Stocco <stocco@subatech.in2p3.fr>, Subatech
/// \date Jan 26, 2017

#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackSmearing.h"
class TRootIOCtor;

class AliTaskMuonTrackSmearing : public AliAnalysisTaskSE {
 public:
  AliTaskMuonTrackSmearing ( TRootIOCtor* ioCtor );
  AliTaskMuonTrackSmearing ( const char *name, Int_t chosenFunc );
  virtual ~AliTaskMuonTrackSmearing();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  /// Get the muon parameter resolution smearer
  AliMuonTrackSmearing& GetMuonTrackSmearing () { return fMuonTrackSmearing; }

 private:
  AliTaskMuonTrackSmearing(const AliTaskMuonTrackSmearing&);
  AliTaskMuonTrackSmearing& operator=(const AliTaskMuonTrackSmearing&);

  AliMuonTrackSmearing fMuonTrackSmearing; ///< Class to smear track resolution

  ClassDef(AliTaskMuonTrackSmearing, 1); // Muon pair analysis
};
#endif
