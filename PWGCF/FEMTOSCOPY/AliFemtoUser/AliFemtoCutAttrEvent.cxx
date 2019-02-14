///
/// \class AliFemtoUser/AliFemtoCutAttrEvent.cxx
///

#include "AliFemtoCutAttrEvent.h"

#ifndef ALIFEMTOCUTATTREVENT_H

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
};

/// Cut on event centrality
struct EventCutAttrCent {
  std::pair<double, double> cent_range;

  bool Pass(const AliFemtoEvent &ev)
    {
      const double cent = ev.CentralityV0();
      return cent_range.first <= cent && cent < cent_range.second;
    }
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
};




} // namespace pwgfemto

#endif
