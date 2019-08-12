///
/// \class AliFemtoUser/AliFemtoCutAttrEvent.h
///

#pragma once

#ifndef ALIFEMTOCUTATTREVENT_H
#define ALIFEMTOCUTATTREVENT_H

#include "AliFemtoConfigObject.h"
#include "AliFemtoEvent.h"

#include <utility>


namespace pwgfemto {

template <typename T1, typename T2>
struct AddEventCutAttrs : public T1 , public T2 {

  bool Pass(const AliFemtoEvent &ev)
    {
      return T1::Pass(ev) && T2::Pass(ev);
    }

  AddEventCutAttrs()
    : T1()
    , T2()
    {}

  AddEventCutAttrs(AliFemtoConfigObject &cfg)
    : T1(cfg)
    , T2(cfg)
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      T1::FillConfiguration(cfg);
      T2::FillConfiguration(cfg);
    }

  virtual ~AddEventCutAttrs() {}
};


/// Cut on Reaction-Plane
///
struct EventCutAttrEpPsi {
  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> ep_psi_range;

  bool Pass(const AliFemtoEvent &ev) const
    {
      const double phi = ev.ReactionPlaneAngle();
      return ep_psi_range.first <= phi && phi < ep_psi_range.second;
    }

  EventCutAttrEpPsi()
    : ep_psi_range(DEFAULT)
    {}

  EventCutAttrEpPsi(AliFemtoConfigObject &cfg)
    : ep_psi_range(cfg.pop_range("ep_psi_range", DEFAULT))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("ep_psi_range", ep_psi_range);
    }

  virtual ~EventCutAttrEpPsi() {}
};

/// cut on event multiplicty
struct EventCutAttrMultiplicty {
  static const std::pair<int, int> DEFAULT;

  std::pair<int, int> mult_range;

  bool Pass(const AliFemtoEvent &ev) const
    {
      const int mult = ev.UncorrectedNumberOfPrimaries();
      return mult_range.first <= mult && mult < mult_range.second;
    }

  EventCutAttrMultiplicty()
    : mult_range(DEFAULT)
    {}

  EventCutAttrMultiplicty(AliFemtoConfigObject &cfg)
    : mult_range(cfg.pop_range("mult_range", DEFAULT))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("mult_range", mult_range);
    }

  virtual ~EventCutAttrMultiplicty() {}
};

/// Cut on event centrality
struct EventCutAttrCentrality {
  static const std::pair<double, double> DEFAULT;
  std::pair<double, double> cent_range;

  bool Pass(const AliFemtoEvent &ev) const
    {
      const double cent = ev.CentralityV0();
      return cent_range.first <= cent && cent < cent_range.second;
    }

  EventCutAttrCentrality()
    : cent_range(DEFAULT)
    {}

  EventCutAttrCentrality(AliFemtoConfigObject &cfg)
    : cent_range(cfg.pop_range("cent_range", DEFAULT))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("cent_range", cent_range);
    }

  virtual ~EventCutAttrCentrality() {}
};


/// Cut on z-position of vertex
///
struct EventCutAttrVertexZ {
  std::pair<double, double> zvert_range;

  bool Pass(const AliFemtoEvent &ev) const
    {
      const double vertex_z = ev.PrimVertPos().z();
      return zvert_range.first <= vertex_z && vertex_z < zvert_range.second;
    }

  EventCutAttrVertexZ()
    : zvert_range(-100.0, 100.0)
    {}

  EventCutAttrVertexZ(AliFemtoConfigObject &cfg)
    : zvert_range(cfg.pop_range("zvert_range", std::make_pair(-100.0, 100.0)))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("zvert_range", zvert_range);
    }

  virtual ~EventCutAttrVertexZ() {}
};

/// Cut bad vertex based on ZDC participants
/// (default requires at least two 2 participants)
///
struct EventCutAttrZdcParticipants {
  unsigned int zdc_participants_min;

  bool Pass(const AliFemtoEvent &ev) const
    {
      return zdc_participants_min <= ev.ZDCParticipants();
    }

  EventCutAttrZdcParticipants()
    : zdc_participants_min(0)
    {}

  EventCutAttrZdcParticipants(AliFemtoConfigObject &cfg)
    : zdc_participants_min(cfg.pop_uint("zdc_participants_min", 0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      if (zdc_participants_min == 0) {
        return;
      }
      cfg.insert("zdc_participants_min", (Long64_t)zdc_participants_min);
    }

  virtual ~EventCutAttrZdcParticipants() {}
};


/// Trigger cut
struct EventCutAttrTrigger {
  unsigned char trigger;

  bool Pass(const AliFemtoEvent &ev) const
    {
      return trigger == 0 || ev.TriggerCluster() == trigger;
    }

  EventCutAttrTrigger()
    : trigger(0)
    {}

  EventCutAttrTrigger(AliFemtoConfigObject &cfg)
    : trigger(cfg.pop_uint("trigger", 0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("trigger", trigger);
    }

  virtual ~EventCutAttrTrigger() {}
};


} // namespace pwgfemto


#include "AliFemtoEventCut.h"


/// \class AliFemtoEventCutAttr
/// \brief Bridge from AliFemtoEventCut to a metaclass of EventCut-Attrs
///
template <typename CRTP, typename CutAttrType>
class AliFemtoEventCutAttr : public AliFemtoEventCut, public CutAttrType {
public:

  typedef CutAttrType CutAttrs;

  AliFemtoEventCutAttr()
    {}

  AliFemtoEventCutAttr(AliFemtoConfigObject &cfg)
    : AliFemtoEventCut()
    , CutAttrType(cfg)
    {}

  virtual bool Pass(const AliFemtoEvent *ev)
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

  void FillConfig(AliFemtoConfigObject &cfg) const
    {
      CutAttrs::FillConfiguration(cfg);
    }

  AliFemtoConfigObject GetConfiguration() const
    {
      AliFemtoConfigObject cfg = AliFemtoConfigObject::BuildMap()
                                  ("_class", static_cast<const CRTP*>(this)->ClassName());
      FillConfig(cfg);
      return cfg;
    }

  virtual void AppendSettings(TCollection &) const = 0;
  virtual ~AliFemtoEventCutAttr()
    {}
};




#endif
