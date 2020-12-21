///
/// \file AliFemtoEventCutCentrality.h
///

#pragma once

#ifndef ALIFEMTOEVENTCUTCENTRALITY_H
#define ALIFEMTOEVENTCUTCENTRALITY_H

#include "AliFemtoEventCut.h"
#include "AliLog.h"

#include <utility> // std::pair
#include "AliFemtoConfigObject.h"


/// \class AliFemtoEventCutCentrality
/// \brief Event cut based on the determined event centrality
///
/// Cuts cuts on event centrality, z-vertex position, and event plane angle (Ψ-EP)
///
/// \author Andrew Kubera, The Ohio State University <andrew.kubera@cern.ch>
///
class AliFemtoEventCutCentrality : public AliFemtoEventCut {
public:

  /// Enumerated type used to select the centrality algorithm.
  /// Look at :class:`AliFemtoEvent` for more information.
  ///
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
  , kUNKNOWN = -9999
  };

  struct Parameters {

    CentralityType centrality_type;
    std::pair<double, double> centrality_range;
    std::pair<double, double> vertex_z_range;
    std::pair<double, double> psi_ep_range;
    int select_trigger;

    Parameters();
    Parameters(AliFemtoConfigObject&);
    // Parameters(const AliFemtoConfigObject &);

    operator AliFemtoEventCutCentrality() const
    {
      return AliFemtoEventCutCentrality(*this);
    }

  };

public:

  /// Default Constructor
  ///
  /// Uses the 'V0' event centrality calculation and wide ranges to
  /// allow all events.
  ///
  AliFemtoEventCutCentrality();

  /// Copy Constructor - Copies the parameters but NOT number of events
  /// which have passed
  ///
  AliFemtoEventCutCentrality(const AliFemtoEventCutCentrality& c);

  /// Parameter Constructor
  AliFemtoEventCutCentrality(const Parameters &);

  /// Build AliFemtoConfigObject
  AliFemtoConfigObject GetConfigObject() const;
  AliFemtoConfigObject* GetConfigObjectPtr() const;

  /// Assignment Operator - Copies the parameters and RESETS the number
  /// of events to 0
  ///
  AliFemtoEventCutCentrality& operator=(const AliFemtoEventCutCentrality& c);

  /// Reset object by applying parameters number of events passed/fail
  ///
  void ResetWithParameters(const Parameters &);

  /// Set min and max acceptable event centrality
  void SetCentralityRange(const float lo, const float hi);

  /// Set min and max acceptable vertex z-coordinate
  void SetZPosRange(const float lo, const float hi);

  /// Set the min and max allowed event reaction plane angle
  void SetEPVZERO(const float lo, const float hi);

  /// Number of events passed
  ULong_t NEventsPassed() const;

  /// Number of events failed
  ULong_t NEventsFailed() const;

  /// Set centrality type using the enum
  void SetCentralityType(const CentralityType);

  /// Set centrality type by string identification (case-insensitive).
  /// If string is unknown, a warning will be printed and no modification
  /// to the cut will be made.
  ///
  void SetCentralityType(const TString& typestr);

  void SetTriggerSelection(int trig);  ///< Set the trigger cluster

  virtual TList* AppendSettings(TList*, const TString& prefix="") const;
  virtual AliFemtoString Report();

  virtual bool Pass(const AliFemtoEvent* event);

  bool PassCentrality(const AliFemtoEvent* event) const;
  bool PassVertex(const AliFemtoEvent* event) const;
  bool PassEventPlane(const AliFemtoEvent* event) const;
  bool PassTrigger(const AliFemtoEvent* event) const;

  virtual AliFemtoEventCut* Clone() const;

  /// Return the centrality of the event based on whatever algorithm is
  /// selected by the CentralityType member.
  ///
  float GetCentrality(const AliFemtoEvent*) const;

  /// Function returnting whether first parameter is within the bounds
  /// set by the second parameter
  template <typename RangeType>
  static bool within_range(float, const RangeType&);

  ///
  static CentralityType CentralityTypeFromName(const TString &);

  std::string CentralityTypeName(const CentralityType) const;

protected:
  typedef std::pair<float, float> Range_t;

  CentralityType fCentralityType;   ///< Selects which centrality calculation algorithm to use
  Range_t fEventCentrality;  ///< range of centrality
  Range_t fVertZPos;         ///< range of z-position of vertex
  Range_t fPsiEP;            ///< range of event plane angle
  int fSelectTrigger;        ///< If set, only given triggers will be selected

  ULong_t fNEventsPassed;  ///< Number of events checked by this cut that passed
  ULong_t fNEventsFailed;  ///< Number of events checked by this cut that failed

private:
  /// Return the name of the class - required for use with AliWarning
  TString ClassName() { return "AliFemtoEventCutCentrality"; }
};

inline void AliFemtoEventCutCentrality::SetCentralityType(const CentralityType type)
{
  fCentralityType = type;
}

inline void AliFemtoEventCutCentrality::SetCentralityType(const TString& typestr)
{
  auto type = CentralityTypeFromName(typestr);
  SetCentralityType(type);
}

inline AliFemtoEventCutCentrality::CentralityType AliFemtoEventCutCentrality::CentralityTypeFromName(const TString &typestr)
{
  TString t = typestr;
  t.ToLower();

  CentralityType type = t == "v0" ? kV0
                      : t == "v0m" ? kV0A
                      : t == "v0a" ? kV0A
                      : t == "v0c" ? kV0C
                      : t == "zna" ? kZNA
                      : t == "znc" ? kZNC
                      : t == "cl0" ? kCL0
                      : t == "cl1" ? kCL1
                      : t == "tkl" ? kTKL
                      : t == "fmd" ? kFMD
                      : t == "trk" ? kTrk
                      : t == "cnd" ? kCND
                      : t == "npa" ? kNPA
                      : t == "spd1" ? kSPD1
                      : kUNKNOWN;

  // if (type == BAD_STRING) {
  //   AliWarning("Bad centrality type string: '" + typestr + "'");
  // }

  return type;
}

inline std::string AliFemtoEventCutCentrality::CentralityTypeName(const CentralityType type) const
{
  #define CASE(val, name) case val: return name;

  switch (type) {
    CASE(kV0, "V0")
    CASE(kV0A, "V0A")
    CASE(kV0C, "V0C")
    CASE(kZNA, "ZNA")
    CASE(kZNC, "ZNC")
    CASE(kCL0, "CL0")
    CASE(kCL1, "CL1")
    CASE(kTKL, "TKL")
    CASE(kFMD, "FMD")
    CASE(kTrk, "TRK")
    CASE(kCND, "CND")
    CASE(kNPA, "NPA")
    CASE(kSPD1, "SPD1")
    default: break;
  }

  #undef CASE

  return "UNKNOWN";
}

inline void AliFemtoEventCutCentrality::SetCentralityRange(const float lo, const float hi)
{
  fEventCentrality = std::make_pair(lo, hi);
}

inline void AliFemtoEventCutCentrality::SetZPosRange(const float lo, const float hi)
{
  fVertZPos = std::make_pair(lo, hi);
}

inline void AliFemtoEventCutCentrality::SetEPVZERO(const float lo, const float hi)
{
  fPsiEP = std::make_pair(lo, hi);
}

inline ULong_t AliFemtoEventCutCentrality::NEventsPassed() const
{
  return fNEventsPassed;
}

inline ULong_t AliFemtoEventCutCentrality::NEventsFailed() const
{
  return fNEventsFailed;
}

inline void AliFemtoEventCutCentrality::SetTriggerSelection(int trig)
{
  fSelectTrigger = trig;
}

inline AliFemtoEventCut* AliFemtoEventCutCentrality::Clone() const
{
  return new AliFemtoEventCutCentrality(*this);
}

inline AliFemtoEventCutCentrality::AliFemtoEventCutCentrality(const AliFemtoEventCutCentrality& c):
  AliFemtoEventCut(c)
  , fCentralityType(c.fCentralityType)
  , fEventCentrality(c.fEventCentrality)
  , fVertZPos(c.fVertZPos)
  , fPsiEP(c.fPsiEP)
  , fSelectTrigger(c.fSelectTrigger)
  , fNEventsPassed(0)
  , fNEventsFailed(0)
{
}

inline AliFemtoEventCutCentrality& AliFemtoEventCutCentrality::operator=(const AliFemtoEventCutCentrality& c)
{
  if (this != &c) {
     AliFemtoEventCut::operator=(c);
     fCentralityType = c.fCentralityType;
     fSelectTrigger = c.fSelectTrigger;
     fEventCentrality = c.fEventCentrality;
     fVertZPos = c.fVertZPos;
     fPsiEP = c.fPsiEP;
     fNEventsPassed = 0;
     fNEventsFailed = 0;
  }

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

template <typename RangeType> inline
bool AliFemtoEventCutCentrality::within_range(float val, const RangeType &range) {
    return (std::get<0>(range) <= val) && (val < std::get<1>(range));
}

inline
bool AliFemtoEventCutCentrality::PassCentrality(const AliFemtoEvent* event) const
{
  return within_range(GetCentrality(event), fEventCentrality);
}

inline
bool AliFemtoEventCutCentrality::PassVertex(const AliFemtoEvent* event) const
{
  return within_range(event->PrimVertPos().z(), fVertZPos);
}

inline
bool AliFemtoEventCutCentrality::PassEventPlane(const AliFemtoEvent* event) const
{
  return within_range(event->ReactionPlaneAngle(), fPsiEP);
}

inline
bool AliFemtoEventCutCentrality::PassTrigger(const AliFemtoEvent* event) const
{
  return (fSelectTrigger == 0) || (event->TriggerCluster() == fSelectTrigger);
}

#endif
