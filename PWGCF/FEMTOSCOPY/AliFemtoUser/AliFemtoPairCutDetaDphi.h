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
/// \author Andrew Kubera <andrew.kubera@cern.ch>
///
class AliFemtoPairCutDetaDphi : public AliFemtoShareQualityPairCut {
public:
  /// Default Constructor
  AliFemtoPairCutDetaDphi();
  AliFemtoPairCutDetaDphi(Float_t min_eta, Float_t min_phi, Float_t radius=1.6f);

  virtual bool Pass(const AliFemtoPair *);
  virtual bool Pass(AliFemtoPair *);

  /// Returns an empty string
  virtual AliFemtoString Report();

  /// Return a list of TObjStrings describing the configuration settings for this cut
  virtual TList* ListSettings();
  virtual TList* AppendSettings(TList *, const TString &prefix="");

  /// Called every time a new event is processed. This method updates the fCurrentMagneticField variable
  virtual void EventBegin(const AliFemtoEvent *event);

  /// Set the radius for $\Delta\phi^{*}$ calculations
  void SetR(Float_t r);

  /// Set the minimum allowed Delta eta value
  void SetMinEta(Float_t eta);

  /// Set the minimum allowed Delta phi value
  void SetMinPhi(Float_t phi);

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
};


inline
AliFemtoPairCutDetaDphi::AliFemtoPairCutDetaDphi():
  AliFemtoShareQualityPairCut()
  , fR(1.6f)
  , fDeltaEtaMin(0.0)
  , fDeltaPhiMin(0.0)
{
  // no-op
}


inline
AliFemtoPairCutDetaDphi::AliFemtoPairCutDetaDphi(Float_t min_eta, Float_t min_phi, Float_t radius):
  AliFemtoShareQualityPairCut()
  , fR(radius)
  , fDeltaEtaMin(min_eta)
  , fDeltaPhiMin(min_phi)
{
  // no-op
}


inline
void AliFemtoPairCutDetaDphi::SetR(Float_t R)
{
  fR = R;
}


inline
void AliFemtoPairCutDetaDphi::SetMinEta(Float_t eta)
{
  fDeltaEtaMin = eta;
}


inline
void AliFemtoPairCutDetaDphi::SetMinPhi(Float_t phi)
{
  fDeltaPhiMin = phi;
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

  const bool within_cut_range = (dEta < fDeltaEtaMin) && (dPhi < fDeltaPhiMin);

  bool passes = !within_cut_range && AliFemtoShareQualityPairCut::Pass(pair);

  return passes;
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
  const double PI_OVER_2 = TMath::Pi() / 2.0,
              PI_TIMES_2 = TMath::Pi() * 2.0,
               RADIUS_CM = radius_in_meters * 100.0;

  double thetas1 = PI_OVER_2 - TMath::ATan(a.z() / RADIUS_CM);
  double thetas2 = PI_OVER_2 - TMath::ATan(b.z() / RADIUS_CM);
  double etas1 = -TMath::Log( TMath::Tan(thetas1/2.) );
  double etas2 = -TMath::Log( TMath::Tan(thetas2/2.) );
  double delta_eta_star = TMath::Abs(etas1 - etas2);
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
  //  q * B * R / pT = (1.60218e-19 * c / 1.60218e-10) [e] [Tesla] [Meters] / [GeV/c]
  //                 ~ 0.3 [e] [Tesla] [Meters] / [GeV/c]
  //

  const Double_t unit_factor = 0.299792458 / 2.0;

  const Double_t phi_a = p_a.Phi(),
                 phi_b = p_b.Phi(),

                 prefix = -1.0 * unit_factor * magnetic_field * radius_in_meters,

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

}

inline
TList*
AliFemtoPairCutDetaDphi::AppendSettings(TList *setting_list, const TString &prefix)
{
  setting_list->AddVector(
    new TObjString(
      prefix + TString::Format("AliFemtoPairCutDetaDphi.radius=%f", fR)
    ),

    new TObjString(
      prefix + TString::Format("AliFemtoPairCutDetaDphi.min_delta_eta=%f", fDeltaEtaMin)
    ),

    new TObjString(
      prefix + TString::Format("AliFemtoPairCutDetaDphi.min_delta_phi=%f", fDeltaPhiMin)
    ),

  NULL);

  return setting_list;
}

inline
AliFemtoString AliFemtoPairCutDetaDphi::Report()
{
  return "";
}


#endif
