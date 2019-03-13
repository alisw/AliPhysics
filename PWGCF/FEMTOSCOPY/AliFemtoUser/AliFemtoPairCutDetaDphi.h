///
/// \file AliFemtoPairCutDetaDphi.h
///

#ifndef ALIFEMTOPAIRCUTDETADPHI_H_
#define ALIFEMTOPAIRCUTDETADPHI_H_

#pragma once

#include "AliFemtoShareQualityPairCut.h"
#include <TList.h>

#include <cmath>


/// \class AliFemtoPairCutDetaDphi
/// \brief A pair cut which cuts on the Δη Δφ of the pair.
///
/// The difference in phi is calculated by examining the tracks' azimuthal angle at a particular radius
/// of all tracks. The default value for this radius is 1.6 (inside the TPC) but this may be changed via
/// the SetR method.
///
///
///    \Delta \phi_{min}* = \phi_1 - \phi_2
///                       + arcsin \left( \frac{ z_1 \cdot B_z \cdot R}{2 p_{T1}} \right)
///                       - arcsin \left( \frac{ z_2 \cdot B_z \cdot R}{2 p_{T2}} \right)
///
///
/// This class encapsulates two methods for cutting pairs: the "simple"
/// method of independently cutting each one, and "quadrature" method
/// which cuts Δη greater than the first parameter and \sqrt{Δη^2 + Δφ^2}
/// greater than the second parameter. The technique may be chosen by
/// using the SetCutTechnique method.
///
/// \author Andrew Kubera <andrew.kubera@cern.ch>
///
class AliFemtoPairCutDetaDphi : public AliFemtoShareQualityPairCut {
public:
  /// Cut technique for
  enum Technique {
    /// Simple cut testing each parameter independently
    Simple,
    /// Tests if the second parameter is within \sqrt{Δη^2 + Δφ^2}
    Quad,
  };

public:
  /// Default Constructor
  AliFemtoPairCutDetaDphi();
  AliFemtoPairCutDetaDphi(Float_t min_eta, Float_t min_phi, Float_t radius=1.6f);

  virtual bool Pass(const AliFemtoPair *);
  virtual bool Pass(AliFemtoPair *);

  /// Returns the simple
  bool PassesSimple(float deta, float dphi_star) const;

  /// Returns we
  bool PassesQuad(float deta, float dphi_star) const;

  /// Returns an empty string
  virtual AliFemtoString Report();

  /// Return a list of TObjStrings describing the configuration settings for this cut
  virtual TList* ListSettings();
  virtual TList* AppendSettings(TList *, const TString &prefix="");

  /// Called every time a new event is processed. This method updates the fCurrentMagneticField variable
  virtual void EventBegin(const AliFemtoEvent *event);

  /// Sets which method to use
  void SetCutTechnique(Technique t);
  Bool_t GetUsesQuadratureTechnique() const;

  /// Set the radius for $\Delta\phi^{*}$ calculations
  void SetR(Float_t r);
  Float_t GetRadius() const;

  /// Set the minimum allowed Delta eta value
  void SetMinEta(Float_t eta);
  Float_t GetMinDeltaEta() const;

  /// Set the minimum allowed Delta phi value
  void SetMinPhi(Float_t phi);
  Float_t GetMinDeltaPhi() const;

  /// Calculate \Delta\eta between two particles.
  /// \param a Momentum of first particle
  /// \param b Momentum of second particle
  ///
  /// The calculation returns $\Delta\eta = \eta_2 - \eta_1$, with each eta caluclated via the
  /// PseudoRapidity method of AliFemtoThreeVector.
  ///
  static Float_t CalculateDEta(const AliFemtoThreeVector& a, const AliFemtoThreeVector& b);

  /// Calculate the \Delta\phi^{*} between two particles.
  /// \param p_a momentum of first particle
  /// \param charge_a charge of the first particle
  /// \param p_b momentum of second particle
  /// \param charge_b charge of the second particle
  /// \param radius_in_meters Radial distance at which the angle should be taken [Meters]
  /// \param magnetic_field_in_tesla Strength and direction of the magnetic field in the detector [Tesla]
  ///
  static Float_t CalculateDPhiStar(const AliFemtoThreeVector& p_a,
                                   const Short_t charge_a,
                                   const AliFemtoThreeVector& p_b,
                                   const Short_t charge_b,
                                   const Double_t radius_in_meters,
                                   const Double_t magnetic_field_in_tesla);

  /// Calculate \Delta\eta between two particles.
  /// \param a Momentum of first particle
  /// \param b Momentum of second particle
  /// \param minRad Radial distance at which the eta value should be taken?
  static Float_t CalculateDEtaStar(const AliFemtoThreeVector& a, const AliFemtoThreeVector& b, const Double_t minRad);

protected:

  /// The magnetic field of the most recent event, updated on EventUpdate
  Float_t fCurrentMagneticField;

  /// Radius at which we calculate φ*
  Float_t fR;

  /// Minimum acceptable Δη
  Float_t fDeltaEtaMin;

  /// Minimum acceptable Δφ*
  Float_t fDeltaPhiMin;

  /// which method to use to cut pairs
  Technique fCutTechnique;
};

inline
void AliFemtoPairCutDetaDphi::SetCutTechnique(Technique t)
{
  // no change
  if (fCutTechnique == t) {
    return;
  }

  // change meaning of fDeltaPhiMin - either limit or square of limit
  switch (t) {
  case Simple:
    fDeltaPhiMin = TMath::Sqrt(fDeltaPhiMin);
    break;
  case Quad:
    fDeltaPhiMin = fDeltaPhiMin * fDeltaPhiMin;
    break;
  }

  fCutTechnique = t;
}

inline
Bool_t AliFemtoPairCutDetaDphi::GetUsesQuadratureTechnique() const
{
  return fCutTechnique == Quad;
}

inline
AliFemtoPairCutDetaDphi::AliFemtoPairCutDetaDphi():
  AliFemtoShareQualityPairCut()
  , fCurrentMagneticField(0.0)
  , fR(1.6f)
  , fDeltaEtaMin(0.0)
  , fDeltaPhiMin(0.0)
  , fCutTechnique(Simple)
{
  // no-op
}


inline
AliFemtoPairCutDetaDphi::AliFemtoPairCutDetaDphi(Float_t min_eta, Float_t min_phi, Float_t radius):
  AliFemtoShareQualityPairCut()
  , fCurrentMagneticField(0.0)
  , fR(radius)
  , fDeltaEtaMin(min_eta)
  , fDeltaPhiMin(min_phi)
  , fCutTechnique(Simple)
{
  // no-op
}


inline
void AliFemtoPairCutDetaDphi::SetR(Float_t R)
{
  fR = R;
}

inline
Float_t AliFemtoPairCutDetaDphi::GetRadius() const
{
  return fR;
}


inline
void AliFemtoPairCutDetaDphi::SetMinEta(Float_t eta)
{
  fDeltaEtaMin = eta;
}


inline
Float_t AliFemtoPairCutDetaDphi::GetMinDeltaEta() const
{
  return fDeltaEtaMin;
}


inline
void AliFemtoPairCutDetaDphi::SetMinPhi(Float_t phi)
{
  fDeltaPhiMin = phi;
  if (fCutTechnique == Quad) {
    fDeltaPhiMin *= fDeltaPhiMin;
  }
}

inline
Float_t AliFemtoPairCutDetaDphi::GetMinDeltaPhi() const
{
  return (fCutTechnique == Quad)
       ? std::sqrt(fDeltaPhiMin)
       : fDeltaPhiMin;
}

inline
bool AliFemtoPairCutDetaDphi::Pass(AliFemtoPair *pair)
{
  return Pass(const_cast<const AliFemtoPair*>(pair));
}


inline
bool AliFemtoPairCutDetaDphi::Pass(const AliFemtoPair *pair)
{
  const AliFemtoTrack *track1 = pair->Track1()->Track(),
                      *track2 = pair->Track2()->Track();

  const AliFemtoThreeVector &p1 = track1->P(),
                            &p2 = track2->P();

  const short charge1 = track1->Charge(),
              charge2 = track2->Charge();

  const Float_t dEta = fabs(CalculateDEta(p1, p2)),
                dPhi = fabs(CalculateDPhiStar(p1, charge1, p2, charge2, fR, fCurrentMagneticField));

  // Initialize to false - so if fCutTechnique is invalid, all pairs
  // will fail and user will eventually find this comment.
  bool passes = false;
  switch (fCutTechnique) {
  case Simple:
    passes = PassesSimple(dEta, dPhi);
    break;
  case Quad:
    passes = PassesQuad(dEta, dPhi);
    break;
  }

  passes = passes && AliFemtoShareQualityPairCut::Pass(pair);

  return passes;
}


inline
bool AliFemtoPairCutDetaDphi::PassesSimple(float delta_eta, float delta_phi_star) const
{
  return (fDeltaEtaMin <= delta_eta) || (fDeltaPhiMin <= delta_phi_star);
}

inline
bool AliFemtoPairCutDetaDphi::PassesQuad(float delta_eta, float delta_phi_star) const
{
  if (fDeltaEtaMin <= delta_eta) {
    return true;
  }
  const float quad = delta_eta * delta_eta + delta_phi_star * delta_phi_star;
  // fDeltaPhiMin stores limit^2 - so we don't have to take the squareroot
  return fDeltaPhiMin <= quad;
}

inline
Float_t AliFemtoPairCutDetaDphi::CalculateDEta(const AliFemtoThreeVector& a,
                                               const AliFemtoThreeVector& b)
{
  const double eta1 = a.PseudoRapidity(),
               eta2 = b.PseudoRapidity();
  const double deta = eta2 - eta1;
  return deta;
}


inline
Float_t AliFemtoPairCutDetaDphi::CalculateDEtaStar(
  const AliFemtoThreeVector& a,
  const AliFemtoThreeVector& b,
  const Double_t radius_in_meters)
{
  const double
    PI_OVER_2 = TMath::Pi() / 2.0,
    RADIUS_CM = radius_in_meters * 100.0,

    thetas1 = PI_OVER_2 - TMath::ATan(a.z() / RADIUS_CM),
    thetas2 = PI_OVER_2 - TMath::ATan(b.z() / RADIUS_CM),

    etas1 = -TMath::Log( TMath::Tan(thetas1 / 2.0) ),
    etas2 = -TMath::Log( TMath::Tan(thetas2 / 2.0) ),

    delta_eta_star = TMath::Abs(etas1 - etas2);

  return delta_eta_star;
}


inline
Float_t AliFemtoPairCutDetaDphi::CalculateDPhiStar(
  const AliFemtoThreeVector& p_a,
  const Short_t charge_a,
  const AliFemtoThreeVector& p_b,
  const Short_t charge_b,
  const Double_t radius_in_meters,
  const Double_t magnetic_field)
{
  // phi shift at radius R in magnetic field B, with charge and momentum q & p_T is:
  //    φ = arcsin( q * B * R / 2 p_T )
  //
  // Unit analysis:
  //   q * B * R : [Coulomb][Tesla][Meter] = 1.60218e-19 [e][Tesla][Meter]
  //   p_T : [Joule][Second]/[Meter] = 1.60218e-10 [GeV][Second]/[Meter] = 1.60218e-10 / c [GeV/c]
  //
  //  q * B * R / p_T = (1.60218e-19 * c / 1.60218e-10) [e] [Tesla] [Meters] / [GeV/c]
  //                 ~ 0.3 [e] [Tesla] [Meters] / [GeV/c]
  //

  const double
    UNIT_FACTOR = 0.299792458,

    phi_a = p_a.Phi(),
    phi_b = p_b.Phi(),

    prefix = -0.5 * UNIT_FACTOR * magnetic_field * radius_in_meters,

    shift_a = TMath::ASin(prefix * charge_a / p_a.Perp()),
    shift_b = TMath::ASin(prefix * charge_b / p_b.Perp());

  const Double_t delta_phi_star = (phi_b + shift_b) - (phi_a + shift_a);

  return delta_phi_star;
}


inline
TList* AliFemtoPairCutDetaDphi::ListSettings()
{
  TList* settings = AliFemtoShareQualityPairCut::ListSettings();
  AppendSettings(settings);

  return settings;
};

inline
void
AliFemtoPairCutDetaDphi::EventBegin(const AliFemtoEvent *event)
{
  fCurrentMagneticField = event->MagneticField();

  // Correct AliFemto units error for magnetic field (back to Tesla)
  // TODO: Fix this bug in AliFemtoEventReaderAOD::CopyAODtoFemtoEvent
  if (fabs(fCurrentMagneticField) < 1e-10) {
    fCurrentMagneticField *= 1e13;
  }

  AliFemtoCutMonitorHandler::EventBegin(event);
}

inline
TList*
AliFemtoPairCutDetaDphi::AppendSettings(TList *setting_list, const TString &prefix)
{
  setting_list->AddVector(
    new TObjString(prefix + Form("AliFemtoPairCutDetaDphi.radius=%g", fR)),
    new TObjString(prefix + Form("AliFemtoPairCutDetaDphi.min_delta_eta=%g", GetMinDeltaEta())),
    new TObjString(prefix + Form("AliFemtoPairCutDetaDphi.min_delta_phi=%g", GetMinDeltaPhi())),
  nullptr);

  return setting_list;
}

inline
AliFemtoString AliFemtoPairCutDetaDphi::Report()
{
  return "";
}


#endif
