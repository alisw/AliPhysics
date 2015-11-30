///
/// \file AliFemtoEventCutCentrality.h
///

#ifndef ALIFEMTOEVENTCUTCENTRALITY_H
#define ALIFEMTOEVENTCUTCENTRALITY_H

#pragma once

#include "AliFemtoEventCut.h"

/// \class AliFemtoEventCutCentrality
/// \brief Event cut based on the determined event centrality
///
/// Cuts cuts on event multiplicity and z-vertex position
///
class AliFemtoEventCutCentrality : public AliFemtoEventCut {
public:

  /**
   * Enumerated type used to select the centrality algorithm
   */
  enum CentralityType {
    kV0
  , kV0A
  , kV0C
  , kZNA
  , kZNC
  , kCL1
  , kCL0
  , kTKL
  , kFMD
  , kTrk
  , kCND
  , kNPA
  , kSPD1
  };

public:
  /**
   * Default Constructor
   *
   * Uses the 'V0' event centrality calcuation.
   */
  AliFemtoEventCutCentrality();

  /**
   * Copy Constructor
   */
  AliFemtoEventCutCentrality(const AliFemtoEventCutCentrality& c);

  /** Destructor */
  virtual ~AliFemtoEventCutCentrality();

  /** Assignment Operator */
  AliFemtoEventCutCentrality& operator=(const AliFemtoEventCutCentrality& c);

  void SetCentralityRange(const int lo,const int hi); ///< Set min and max acceptable event centrality
  void SetZPosRange(const float lo, const float hi);  ///< Set min and max acceptable vertex z-coordinate
  void SetEPVZERO(const float lo, const float hi);    ///< Set the min and max allowed event reaction plane angle

  int NEventsPassed() const;  ///< Number of events passed
  int NEventsFailed() const;  ///< Number of events failed

  void SetCentralityType(const CentralityType);
  void SetCentralityType(TString typestr);

  void SetTriggerSelection(int trig);  ///< Set the trigger cluster


  virtual TList* AppendSettings(TList*, const TString& prefix="") const;
  virtual AliFemtoString Report();

  virtual bool Pass(const AliFemtoEvent* event);

  bool PassCentrality(const AliFemtoEvent* event) const;
  bool PassVertex(const AliFemtoEvent* event) const;
  bool PassEventPlane(const AliFemtoEvent* event) const;
  bool PassTrigger(const AliFemtoEvent* event) const;

  AliFemtoEventCutCentrality* Clone() const;

  /**
   * Return the centrality of the event based on algorithm selected in the
   * event cut
   */
  float GetCentrality(const AliFemtoEvent*) const;

protected:

  float fEventCentrality[2];      ///< range of centrality
  CentralityType fCentralityType; ///< Selects which centrality calculation algorithm to use
  float fVertZPos[2];             ///< range of z-position of vertex
  float fPsiEP[2];                ///< range of event plane angle
  int fSelectTrigger;             ///< If set, only given triggers will be selected

  long fNEventsPassed;  ///< Number of events checked by this cut that passed
  long fNEventsFailed;  ///< Number of events checked by this cut that failed

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventCutCentrality, 0);
  /// \endcond
#endif

};

inline void AliFemtoEventCutCentrality::SetCentralityRange(const int lo, const int hi)
{
  fEventCentrality[0] = lo;
  fEventCentrality[1] = hi;
}

inline void AliFemtoEventCutCentrality::SetCentralityType(const CentralityType type)
{
  fCentralityType = type;
}

inline void AliFemtoEventCutCentrality::SetCentralityType(TString typestr)
{
  typestr.ToLower();
  const TString &t = typestr;

  if (t == "v0") {
    SetCentralityType(kV0);
  } else if (t == "v0a") {
    SetCentralityType(kV0A);
  } else if (t == "V0Collection") {
    SetCentralityType(kV0C);
  } else if (t == "zna") {
    SetCentralityType(kZNA);
  } else if (t == "znc") {
    SetCentralityType(kZNC);
  } else if (t == "cl0") {
    SetCentralityType(kCL0);
  } else if (t == "cl1") {
    SetCentralityType(kCL1);
  } else if (t == "tkl") {
    SetCentralityType(kTKL);
  } else if (t == "fmd") {
    SetCentralityType(kFMD);
  } else if (t == "trk") {
    SetCentralityType(kTrk);
  } else if (t == "cnd") {
    SetCentralityType(kCND);
  } else if (t == "npa") {
    SetCentralityType(kNPA);
  } else if (t == "spd1") {
    SetCentralityType(kSPD1);
  } else {

  }
}

inline void AliFemtoEventCutCentrality::SetZPosRange(const float lo, const float hi)
{
  fVertZPos[0] = lo;
  fVertZPos[1] = hi;
}

inline void AliFemtoEventCutCentrality::SetEPVZERO(const float lo, const float hi)
{
  fPsiEP[0] = lo;
  fPsiEP[1] = hi;
}

inline int AliFemtoEventCutCentrality::NEventsPassed() const
{
  return fNEventsPassed;
}

inline int AliFemtoEventCutCentrality::NEventsFailed() const
{
  return fNEventsFailed;
}

inline void AliFemtoEventCutCentrality::SetTriggerSelection(int trig)
{
  fSelectTrigger = trig;
}

inline AliFemtoEventCutCentrality* AliFemtoEventCutCentrality::Clone() const
{
  return new AliFemtoEventCutCentrality(*this);
}

inline AliFemtoEventCutCentrality::AliFemtoEventCutCentrality(const AliFemtoEventCutCentrality& c):
  AliFemtoEventCut(c)
  , fCentralityType(c.fCentralityType)
  , fSelectTrigger(c.fSelectTrigger)
  , fNEventsPassed(0)
  , fNEventsFailed(0)
{
  fEventCentrality[0] = c.fEventCentrality[0];
  fEventCentrality[1] = c.fEventCentrality[1];
  fVertZPos[0] = c.fVertZPos[0];
  fVertZPos[1] = c.fVertZPos[1];
  fPsiEP[0] = c.fPsiEP[0];
  fPsiEP[1] = c.fPsiEP[1];
}

inline AliFemtoEventCutCentrality& AliFemtoEventCutCentrality::operator=(const AliFemtoEventCutCentrality& c)
{
  if (this != &c) {
    return *this;
  }

  AliFemtoEventCut::operator=(c);
  fSelectTrigger = c.fSelectTrigger;
  fEventCentrality[0] = c.fEventCentrality[0];
  fEventCentrality[1] = c.fEventCentrality[1];
  fCentralityType = c.fCentralityType;
  fVertZPos[0] = c.fVertZPos[0];
  fVertZPos[1] = c.fVertZPos[1];
  fPsiEP[0] = c.fPsiEP[0];
  fPsiEP[1] = c.fPsiEP[1];

  return *this;
}

inline
float AliFemtoEventCutCentrality::GetCentrality(const AliFemtoEvent *ev) const
{
  switch (fCentralityType) {
  case kV0: return ev->CentralityV0();
  case kV0A: return ev->CentralityV0A();
  case kV0C: return ev->CentralityV0C();
  case kZNA: return ev->CentralityZNA();
  case kZNC: return ev->CentralityZNC();
  case kCL1: return ev->CentralityCL1();
  case kCL0: return ev->CentralityCL0();
  case kTKL: return ev->CentralityTKL();
  case kFMD: return ev->CentralityFMD();
  case kTrk: return ev->CentralityTrk();
  case kCND: return ev->CentralityCND();
  case kNPA: return ev->CentralityNPA();
  case kSPD1: return ev->CentralitySPD1();
  default:
    return -1.0;
  }
}

inline
bool AliFemtoEventCutCentrality::PassCentrality(const AliFemtoEvent* event) const
{
  const float cent = GetCentrality(event);
  return (fEventCentrality[0] <= cent) && (cent < fEventCentrality[1]);
}

inline
bool AliFemtoEventCutCentrality::PassVertex(const AliFemtoEvent* event) const
{
  const float vertex_z = event->PrimVertPos().z();
  return (fVertZPos[0] < vertex_z) && (vertex_z < fVertZPos[1]);
}

inline
bool AliFemtoEventCutCentrality::PassEventPlane(const AliFemtoEvent* event) const
{
  const double epvzero = event->ReactionPlaneAngle();
  return (fPsiEP[0] < epvzero && epvzero < fPsiEP[1]);
}

inline
bool AliFemtoEventCutCentrality::PassTrigger(const AliFemtoEvent* event) const
{
  return fSelectTrigger && (event->TriggerCluster() == fSelectTrigger);
}

#endif
