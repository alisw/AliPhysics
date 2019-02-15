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

#include <TList.h>


namespace pwgfemto {

// event-cut
typedef AddEventCutAttrs<
          EventCutAttrMultiplicty,
          // AddEventCutAttrs<
          //   EventCutAttrVertexZ,
            EventCutAttrTrigger
          // >
        > EventCutAttrsAK;

// track-cut
typedef TrackSelectionCut<
          // non-physics properties
          AddTrackCutAttrs<
            TrackCutAttrStatus,
            AddTrackCutAttrs<
              TrackCutAttrRemoveKinks,
              TrackCutAttrRemoveFakeITS> >,

          // physics-cuts
          AddTrackCutAttrs<
            TrackCutAttrImpact,
            TrackCutAttrCharge>

        > TrackCutAttrsAK;

// paircut
typedef AddPairCutAttrs<
          PairCutTrackAttrSameLabel,
            AddPairCutAttrs<
              PairCutTrackAttrPt,
              PairCutTrackAttrShareQuality> > PairCutAttrsBaseAK;

typedef AddPairCutAttrs<
          PairCutAttrsBaseAK,
          PairCutTrackAttrAvgSep> PairCutAttrsAvgSepAK;

typedef AddPairCutAttrs<
          PairCutAttrsBaseAK,
          PairCutTrackAttrDetaDphi> PairCutAttrsDphiDetaAK;

}  // namespace pwgfemto


/// \class AliFemtoEventCutAttr
/// \brief Bridge from AliFemtoPairCut to a metaclass of PairCut-Attrs
///
template <typename CRTP, typename CutAttrType>
class AliFemtoEventCutAttr : public AliFemtoEventCut, public CutAttrType {
public:

  typedef CutAttrType CutAttrs;

  virtual bool Pass(AliFemtoEvent *ev)
    {
      return CutAttrs::Pass(*ev);
    }

  virtual AliFemtoString Report()
    {
      return "AliFemtoEventCutAttr Report\n";
    }

  virtual TList* ListSettings() const
    {
      TList* list = new TList();
      AppendSettings(*list);
      return list;
    }

  virtual void AppendSettings(TCollection &) const = 0;
  virtual ~AliFemtoEventCutAttr() = 0;
};


/// \class AliFemtoEventCutPionPionAK
/// \brief Andrew Kubera's event cut for identical-pion analysis
///
///
class AliFemtoEventCutPionPionAK : public AliFemtoEventCutAttr<AliFemtoEventCutPionPionAK, pwgfemto::EventCutAttrsAK> {
public:

  typedef pwgfemto::EventCutAttrsAK CutAttrs;

  virtual ~AliFemtoEventCutPionPionAK()
    { }

  virtual void AppendSettings(TCollection &) const;

};


/// \class AliFemtoEventCutAttr
/// \brief Bridge from AliFemtoPairCut to a metaclass of PairCut-Attrs
///
template <typename CRTP, typename CutAttrType>
class AliFemtoTrackCutAttr : public AliFemtoTrackCut, public CutAttrType {
public:

  typedef CutAttrType CutAttrs;

  virtual bool Pass(AliFemtoTrack *ev)
    {
      return CutAttrs::Pass(*ev);
    }

  virtual AliFemtoString Report()
    {
      return "AliFemtoTrackCutAttr Report\n";
    }

  virtual TList* ListSettings() const
    {
      TList* list = new TList();
      AppendSettings(*list);
      return list;
    }

  virtual void AppendSettings(TCollection &) const = 0;
  virtual ~AliFemtoTrackCutAttr() = 0;
};


/// \class AliFemtoTrackCutPionPionAK
/// \brief Andrew Kubera's pion-track cut
///
class AliFemtoTrackCutPionPionAK : public AliFemtoTrackCutAttr<AliFemtoTrackCutPionPionAK, pwgfemto::TrackCutAttrsAK> {
public:

  typedef pwgfemto::TrackCutAttrsAK CutAttrs;

  virtual ~AliFemtoTrackCutPionPionAK()
    { }

  virtual void AppendSettings(TCollection &) const;

  virtual bool Pass(const AliFemtoTrack *track)
    {
      bool passes = CutAttrs::PassSelection(*track);
      if (passes) {
        passes = CutAttrs::PassCut(*track);
        ++(passes ? fNumPass : fNumFail);
      }
      return passes;
    }

  ULong_t fNumPass,
          fNumFail;
};


/// \class AliFemtoPairCutAttr
/// \brief Bridge from AliFemtoPairCut to a metaclass of PairCut-Attrs
///
/// Subclass and implement your method:
///  `const char* GetName() const`
///
///
template <typename CRTP, typename CutAttrType>
class AliFemtoPairCutAttr : public AliFemtoPairCut, public CutAttrType {
public:

  typedef CutAttrType CutAttrs;

  virtual ~AliFemtoPairCutAttr()
    { }

  AliFemtoPairCutAttr();

  /// user-written method to return string describing cuts
  virtual AliFemtoString Report()
    { return ""; }

  /// Return a TList of settings
  virtual TList* ListSettings()
    {
      TList* list = new TList();
      AppendSettings(*list);
      return list;
    }

  virtual void AppendSettings(TCollection &) const = 0;

  virtual bool Pass(const AliFemtoPair *pair)
    {
      return CutAttrs::Pass(*pair->Track1()->Track(), *pair->Track2()->Track());
    }

  void StoreConfiguration(AliFemtoConfigObject &cfg) const
    {
      CutAttrs::FillConfiguration(cfg);
      cfg.insert("class", static_cast<CRTP*>(this)->GetName());
    }

};


/// \class AliFemtoPairCutPionPionAKAvgSep
/// \brief Andrew Kubera's average separation pair cut
///
class AliFemtoPairCutPionPionAKAvgSep : public AliFemtoPairCutAttr<AliFemtoPairCutPionPionAKAvgSep, pwgfemto::PairCutAttrsAvgSepAK> {
public:

  typedef AliFemtoPairCutAttr<AliFemtoPairCutPionPionAKAvgSep, pwgfemto::PairCutAttrsAvgSepAK> Super;

  virtual ~AliFemtoPairCutPionPionAKAvgSep();

  AliFemtoPairCutPionPionAKAvgSep();

  /// user-written method to return string describing cuts
  virtual AliFemtoString Report();

  virtual void AppendSettings(TCollection &) const;

  const char* GetName() const
    { return "AliFemtoPairCutPionPIonAkAvgSep"; }
};


/// \class AliFemtoPairCutPionPionAKDetaDphi
/// \brief Andrew Kubera's Deta-Dphi pair cut
///
class AliFemtoPairCutPionPionAKDetaDphi : public AliFemtoPairCutAttr<AliFemtoPairCutPionPionAKDetaDphi, pwgfemto::PairCutAttrsDphiDetaAK> {
public:

  typedef AliFemtoPairCutAttr<AliFemtoPairCutPionPionAKDetaDphi, pwgfemto::PairCutAttrsDphiDetaAK> Super;

  virtual ~AliFemtoPairCutPionPionAKDetaDphi()
    { }

  AliFemtoPairCutPionPionAKDetaDphi();

  virtual void EventBegin(const AliFemtoEvent *ev);

  /// user-written method to return string describing cuts
  virtual AliFemtoString Report();

  virtual void AppendSettings(TCollection &) const;

  const char* GetName() const
    { return "AliFemtoPairCutPionPionAKDetaDphi"; }
};



template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoEventCutPionPionAK &cut);

template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoTrackCutPionPionAK &cut);

template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoPairCutPionPionAKAvgSep &cut);

template <>
AliFemtoConfigObject AliFemtoConfigObject::From(const AliFemtoPairCutPionPionAKDetaDphi &cut);

template <>
AliFemtoEventCutPionPionAK* AliFemtoConfigObject::Into<AliFemtoEventCutPionPionAK>(bool);

template <>
AliFemtoTrackCutPionPionAK* AliFemtoConfigObject::Into<AliFemtoTrackCutPionPionAK>(bool);

template <>
AliFemtoPairCutPionPionAKAvgSep* AliFemtoConfigObject::Into<AliFemtoPairCutPionPionAKAvgSep>(bool);

template <>
AliFemtoPairCutPionPionAKDetaDphi* AliFemtoConfigObject::Into<AliFemtoPairCutPionPionAKDetaDphi>(bool);



#endif
