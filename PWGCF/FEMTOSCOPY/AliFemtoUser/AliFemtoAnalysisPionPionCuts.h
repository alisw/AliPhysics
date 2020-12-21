///
/// \file AliFemtoAnalysisPionPionCuts.h
///

#pragma once

#ifndef ALIFEMTOANALYSISPIONPIONCUTS_H
#define ALIFEMTOANALYSISPIONPIONCUTS_H

#include "AliFemtoEventCut.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoPairCut.h"

#include "AliFemtoCutAttrEvent.h"
#include "AliFemtoCutAttrTrack.h"
#include "AliFemtoCutAttrPairTrack.h"

#include "AliFemtoModelHiddenInfo.h"

#include <TList.h>


namespace pwgfemto {

// event-cut
typedef AddEventCutAttrs< EventCutAttrCentrality,
        AddEventCutAttrs< EventCutAttrVertexZ,
        AddEventCutAttrs< EventCutAttrTrigger,
                          EventCutAttrZdcParticipants
                          > > >
        EventCutAttrsAK;

// track-cut
typedef TrackSelectionCut<
          // track-selection properties
          AddTrackCutAttrs< TrackCutAttrStatus,
          AddTrackCutAttrs< TrackCutAttrRemoveNegLabel,
          AddTrackCutAttrs< TrackCutAttrMinNclsTPC,
          AddTrackCutAttrs< TrackCutAttrMinNclsITS,
                            TrackCutAttrCharge
                            > > > >,

          // physical-cuts
          AddTrackCutAttrs< TrackCutAttrRejectTpcElectron,
          AddTrackCutAttrs< TrackCutAttrImpact,
          AddTrackCutAttrs< TrackCutAttrPt,
          AddTrackCutAttrs< TrackCutAttrEta,

          // PID cuts
          AddTrackCutAttrs< TrackCutAttrChi2TPC,
          AddTrackCutAttrs< TrackCutAttrChi2ITS,
          AddTrackCutAttrs< TrackCutAttrRejectKaonTofSigma,
          AddTrackCutAttrs< TrackCutAttrRejectProtonTofSigma,
          AddTrackCutAttrs< TrackCutAttrTofSigmaPion,
                            TrackCutAttrTpcSigmaPion
                            > > > > > > > > >
        > TrackCutAttrsAK;

// paircut
typedef AddPairCutAttrs< PairCutTrackAttrSameLabel,
        AddPairCutAttrs< PairCutTrackAttrPt,
        AddPairCutAttrs< PairCutTrackAttrRemoveEE,
                         PairCutTrackAttrShareQuality
                         > > >
        PairCutAttrsBaseAK;

typedef AddPairCutAttrs< PairCutAttrsBaseAK,
                         PairCutTrackAttrAvgSep
                         >
        PairCutAttrsAvgSepAK;

typedef AddPairCutAttrs< PairCutAttrsBaseAK,
                         PairCutTrackAttrDetaDphiStar
                         >
        PairCutAttrsDphiDetaAK;

}  // namespace pwgfemto


/// \class AliFemtoEventCutPionPionAK
/// \brief Andrew Kubera's event cut for identical-pion analysis
///
///
class AliFemtoEventCutPionPionAK : public AliFemtoEventCutAttr<AliFemtoEventCutPionPionAK, pwgfemto::EventCutAttrsAK> {
public:

  typedef AliFemtoEventCutAttr<AliFemtoEventCutPionPionAK, pwgfemto::EventCutAttrsAK> Super;
  typedef pwgfemto::EventCutAttrsAK CutAttrs;

  AliFemtoEventCutPionPionAK()
    : Super()
    {}

  AliFemtoEventCutPionPionAK(AliFemtoConfigObject &cfg)
    : Super(cfg)
    {}

  virtual ~AliFemtoEventCutPionPionAK()
    {}

  virtual void AppendSettings(TCollection &) const;

  virtual const char* ClassName() const
    { return "AliFemtoEventCutPionPionAK"; }
};


/// \class AliFemtoTrackCutPionPionAK
/// \brief Andrew Kubera's pion-track cut
///
class AliFemtoTrackCutPionPionAK : public AliFemtoTrackCutAttr<AliFemtoTrackCutPionPionAK, pwgfemto::TrackCutAttrsAK> {
public:
  typedef AliFemtoTrackCutAttr<AliFemtoTrackCutPionPionAK, pwgfemto::TrackCutAttrsAK> Super;
  typedef pwgfemto::TrackCutAttrsAK CutAttrs;

  AliFemtoTrackCutPionPionAK()
    : fNumPass(0)
    , fNumFail(0)
    {
      SetMass(0.139570);
    }

  AliFemtoTrackCutPionPionAK(AliFemtoConfigObject &cfg)
    : Super(cfg)
    , fNumPass(0)
    , fNumFail(0)
    {
      SetMass(0.139570);
    }

  virtual ~AliFemtoTrackCutPionPionAK()
    { }

  virtual void AppendSettings(TCollection &, TString prefix="") const;

  virtual bool Pass(const AliFemtoTrack *track)
    {
      bool passes = CutAttrs::PassSelection(*track);
      if (passes) {
        passes = CutAttrs::PassCut(*track);
        ++(passes ? fNumPass : fNumFail);
      }
      return passes;
    }

  virtual const char* ClassName() const
    { return "AliFemtoTrackCutPionPionAK"; }

  ULong_t fNumPass,
          fNumFail;
};

/// Cut on MonteCarlo PDG code to select *only* pions
///
class AliFemtoTrackCutPionPionIdealAK : public AliFemtoTrackCutPionPionAK {
public:

  int ideal_pid;

  AliFemtoTrackCutPionPionIdealAK()
    : AliFemtoTrackCutPionPionAK()
    , ideal_pid(211)
    {}

  AliFemtoTrackCutPionPionIdealAK(AliFemtoConfigObject &cfg)
    : AliFemtoTrackCutPionPionAK(cfg)
    , ideal_pid(cfg.pop_num("ideal_pid", 211))
    {}

  virtual ~AliFemtoTrackCutPionPionIdealAK() {}

  virtual bool Pass(const AliFemtoTrack *track)
    {
      const AliFemtoModelHiddenInfo *info = static_cast<AliFemtoModelHiddenInfo*>(track->GetHiddenInfo());
      const Int_t pid = info->GetPDGPid();

      if (pid != std::copysign(ideal_pid, charge)) {
        return false;
      }
      return AliFemtoTrackCutPionPionAK::Pass(track);
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      AliFemtoTrackCutPionPionAK::FillConfiguration(cfg);
      cfg.insert("ideal_pid", ideal_pid);
    }

  virtual const char* ClassName() const
    { return "AliFemtoTrackCutPionPionIdealAK"; }
};


/// Cut on MonteCarlo PDG code to select *only* misidentified pions
///
class AliFemtoTrackCutPionPionMisidentAK : public AliFemtoTrackCutPionPionAK {
public:

  AliFemtoTrackCutPionPionMisidentAK()
    : AliFemtoTrackCutPionPionAK()
    {}

  AliFemtoTrackCutPionPionMisidentAK(AliFemtoConfigObject &cfg)
    : AliFemtoTrackCutPionPionAK(cfg)
    {}

  virtual ~AliFemtoTrackCutPionPionMisidentAK() {}

  virtual bool Pass(const AliFemtoTrack *track)
    {
      const AliFemtoModelHiddenInfo *info = static_cast<AliFemtoModelHiddenInfo*>(track->GetHiddenInfo());
      const Int_t pid = info->GetPDGPid();
      if (pid == std::copysign(211, charge)) {
        return false;
      }
      return AliFemtoTrackCutPionPionAK::Pass(track);
    }

  virtual const char* ClassName() const
    { return "AliFemtoTrackCutPionPionMisidentAK"; }
};


/// \class AliFemtoPairCutPionPionAKAvgSep
/// \brief Andrew Kubera's average separation pair cut
///
class AliFemtoPairCutPionPionAKAvgSep : public AliFemtoPairCutAttrTracks<AliFemtoPairCutPionPionAKAvgSep,
                                                                         pwgfemto::PairCutAttrsAvgSepAK> {
public:

  typedef AliFemtoPairCutAttrTracks<AliFemtoPairCutPionPionAKAvgSep, pwgfemto::PairCutAttrsAvgSepAK> Super;

  virtual ~AliFemtoPairCutPionPionAKAvgSep();

  AliFemtoPairCutPionPionAKAvgSep()
    : Super()
    {}

  AliFemtoPairCutPionPionAKAvgSep(AliFemtoConfigObject &cfg)
    : Super(cfg)
    {}

  /// user-written method to return string describing cuts
  virtual AliFemtoString Report();

  virtual void AppendSettings(TCollection &) const;

  virtual const char* ClassName() const
    { return "AliFemtoPairCutPionPionAKAvgSep"; }
};


/// \class AliFemtoPairCutPionPionAKDetaDphi
/// \brief Andrew Kubera's Deta-Dphi pair cut
///
class AliFemtoPairCutPionPionAKDetaDphi : public AliFemtoPairCutAttrTracks<AliFemtoPairCutPionPionAKDetaDphi,
                                                                           pwgfemto::PairCutAttrsDphiDetaAK> {
public:

  typedef AliFemtoPairCutAttrTracks<AliFemtoPairCutPionPionAKDetaDphi, pwgfemto::PairCutAttrsDphiDetaAK> Super;

  virtual ~AliFemtoPairCutPionPionAKDetaDphi()
    { }

  AliFemtoPairCutPionPionAKDetaDphi()
    : Super()
    { }

  AliFemtoPairCutPionPionAKDetaDphi(AliFemtoConfigObject &cfg)
    : Super(cfg)
    { }

  virtual void EventBegin(const AliFemtoEvent *ev);

  /// user-written method to return string describing cuts
  virtual AliFemtoString Report();

  virtual void AppendSettings(TCollection &) const;

  virtual const char* ClassName() const
    { return "AliFemtoPairCutPionPionAKDetaDphi"; }
};



template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoEventCutPionPionAK &cut);

template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoTrackCutPionPionAK &cut);

template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoTrackCutPionPionIdealAK &cut);

template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoPairCutPionPionAKAvgSep &cut);

template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoTrackCutPionPionMisidentAK &cut);

template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoPairCutPionPionAKDetaDphi &cut);

template <>
AliFemtoEventCutPionPionAK* AliFemtoConfigObject::Into<AliFemtoEventCutPionPionAK>(bool);

template <>
AliFemtoTrackCutPionPionAK* AliFemtoConfigObject::Into<AliFemtoTrackCutPionPionAK>(bool);

template <>
AliFemtoTrackCutPionPionIdealAK* AliFemtoConfigObject::Into<AliFemtoTrackCutPionPionIdealAK>(bool);

template <>
AliFemtoTrackCutPionPionMisidentAK* AliFemtoConfigObject::Into<AliFemtoTrackCutPionPionMisidentAK>(bool);

template <>
AliFemtoPairCutPionPionAKAvgSep* AliFemtoConfigObject::Into<AliFemtoPairCutPionPionAKAvgSep>(bool);

template <>
AliFemtoPairCutPionPionAKDetaDphi* AliFemtoConfigObject::Into<AliFemtoPairCutPionPionAKDetaDphi>(bool);




#endif
