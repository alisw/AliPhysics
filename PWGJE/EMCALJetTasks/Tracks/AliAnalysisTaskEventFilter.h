/*
 * AliAnalysisTaskEventFilter.h
 *
 *  Created on: Jan 29, 2016
 *      Author: markus
 */

#ifndef ALIANALYSISTASKEVENTFILTER_H
#define ALIANALYSISTASKEVENTFILTER_H

#include <AliAnalysisTaskSE.h>

#include <vector>

class TArray;
class THistManager;
class TList;

class AliAnalysisUtils;
class AliESDtrackCuts;
class AliVTrack;
class AliVEvent;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEventFilter: public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEventFilter();
  AliAnalysisTaskEventFilter(const char *name);
  virtual ~AliAnalysisTaskEventFilter();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);

protected:
  Bool_t FakeVertexSelection2013pA(const AliVEvent * const inputevent) const;
  void CreatePtBinning(TArrayD& binning) const;

  UChar_t FilterEvent() const;
  std::vector<const AliVTrack *> FilterTracks() const;
  void FillEvent(const char *filterstep, double vz);
  void FillTracks(const char *filterstep, const std::vector<const AliVTrack *> &tracks);

private:
  AliAnalysisUtils          *fAnalysisUtils;
  AliESDtrackCuts           *fTrackCuts;
  THistManager              *fHistos;

  AliAnalysisTaskEventFilter(const AliAnalysisTaskEventFilter&);
  AliAnalysisTaskEventFilter& operator=(const AliAnalysisTaskEventFilter&);

  ClassDef(AliAnalysisTaskEventFilter, 1);
};

}

}

#endif /* ALIANALYSISTASKEVENTFILTER_H */
