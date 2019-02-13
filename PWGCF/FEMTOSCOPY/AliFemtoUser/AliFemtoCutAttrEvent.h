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

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      T1::FillConfiguration(cfg);
      T2::FillConfiguration(cfg);
    }

  virtual ~AddEventCutAttrs() = 0;
};

struct EventCutAttrEpPsi {
  std::pair<double, double> ep_psi_range;

  EventCutAttrEpPsi()
    : ep_psi_range(-1000.0, 1000.0)
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
  std::pair<int, int> mult_range;

  EventCutAttrMultiplicty()
    : mult_range(0, 100000)
    {}

  bool Pass(const AliFemtoEvent &ev)
    {
      const int mult = ev.UncorrectedNumberOfPrimaries();
      return mult_range.first <= mult && mult < mult_range.second;
    }

  virtual ~EventCutAttrMultiplicty() = 0;
};

/// Cut on event centrality
struct EventCutAttrCent {
  std::pair<double, double> cent_range;

  bool Pass(const AliFemtoEvent &ev)
    {
      const double cent = ev.CentralityV0();
      return cent_range.first <= cent && cent < cent_range.second;
    }

  virtual ~EventCutAttrCent() = 0;
};

struct EventCutAttrVertexZ {
  std::pair<double, double> zvert_range;

  EventCutAttrVertexZ()
    : zvert_range(-100.0, 100.0)
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

  bool Pass(const AliFemtoEvent &ev)
    {
      return trigger == 0 || ev.TriggerCluster() == trigger;
    }

  virtual ~EventCutAttrTrigger() = 0;
};




} // namespace pwgfemto

#endif
