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

  AddEventCutAttrs()
    : T1()
    , T2()
    {}

  AddEventCutAttrs(AliFemtoConfigObject &cfg)
    : T1(cfg)
    , T2(cfg)
    {}

  bool Pass(const AliFemtoEvent &ev)
    {
      return T1::Pass(ev) && T2::Pass(ev);
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      T1::FillConfiguration(cfg);
      T2::FillConfiguration(cfg);
    }

  virtual ~AddEventCutAttrs() = 0;
};

struct EventCutAttrEpPsi {
  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> ep_psi_range;

  EventCutAttrEpPsi()
    : ep_psi_range(-1000.0, 1000.0)
    {}

  EventCutAttrEpPsi(AliFemtoConfigObject &cfg)
    : ep_psi_range(cfg.pop_range("ep_psi_range", DEFAULT))
    {}

  bool Pass(const AliFemtoEvent &ev)
    {
      const double epvzero = ev.ReactionPlaneAngle();
      return ep_psi_range.first <= epvzero && epvzero < ep_psi_range.second;
    }

  virtual ~EventCutAttrEpPsi() = 0;
};

/// cut on event multiplicty
struct EventCutAttrMultiplicty {
  static const std::pair<int, int> DEFAULT;
  std::pair<int, int> mult_range;

  EventCutAttrMultiplicty()
    : mult_range(0, 100000)
    {}

  EventCutAttrMultiplicty(AliFemtoConfigObject &cfg)
    : mult_range(cfg.pop_range("mult_range", std::make_pair(0, 100000)))
    {}

  bool Pass(const AliFemtoEvent &ev)
    {
      const int mult = ev.UncorrectedNumberOfPrimaries();
      return mult_range.first <= mult && mult < mult_range.second;
    }

  virtual ~EventCutAttrMultiplicty() = 0;
};

/// Cut on event centrality
struct EventCutAttrCentrality {
  std::pair<double, double> cent_range;

  EventCutAttrCentrality()
    : cent_range(0.0, 100.0)
    {}

  EventCutAttrCentrality(AliFemtoConfigObject &cfg)
    : cent_range(cfg.pop_range("cent_range", std::make_pair(0.0, 100.0)))
    {}

  bool Pass(const AliFemtoEvent &ev)
    {
      const double cent = ev.CentralityV0();
      return cent_range.first <= cent && cent < cent_range.second;
    }

  virtual ~EventCutAttrCentrality() = 0;
};

struct EventCutAttrVertexZ {
  std::pair<double, double> zvert_range;

  EventCutAttrVertexZ()
    : zvert_range(-100.0, 100.0)
    {}

  EventCutAttrVertexZ(AliFemtoConfigObject &cfg)
    : zvert_range(cfg.pop_range("zvert_range", std::make_pair(-100.0, 100.0)))
    {}

  bool Pass(const AliFemtoEvent &ev)
    {
      const double vertex_z = ev.PrimVertPos().z();
      return zvert_range.first <= vertex_z && vertex_z < zvert_range.second;
    }

  virtual ~EventCutAttrVertexZ() = 0;
};

/// Cut bad vertex based on ZDC participants
/// (default requires at least two 2 participants)
///
struct EventCutAttrZdcParticipants {
  unsigned int min_zdc_participants;

  EventCutAttrZdcParticipants()
    : min_zdc_participants(2)
    {}

  EventCutAttrZdcParticipants(AliFemtoConfigObject &cfg)
    : min_zdc_participants(cfg.pop_uint("min_zdc_participants", 2))
    {}

  bool Pass(const AliFemtoEvent &ev)
    {
      return ev.ZDCParticipants() >= min_zdc_participants;
    }

  virtual ~EventCutAttrZdcParticipants() = 0;
};


/// Trigger cut
struct EventCutAttrTrigger {
  unsigned char trigger;

  EventCutAttrTrigger()
    : trigger(0)
    {}

  EventCutAttrTrigger(AliFemtoConfigObject &cfg)
    : trigger(cfg.pop_uint("trigger", 0))
    {}

  bool Pass(const AliFemtoEvent &ev)
    {
      return trigger == 0 || ev.TriggerCluster() == trigger;
    }

  virtual ~EventCutAttrTrigger() = 0;
};




} // namespace pwgfemto

#endif
