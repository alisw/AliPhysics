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
  AliFemtoPairCutDetaDphi();
  AliFemtoPairCutDetaDphi(Float_t min_eta, Float_t min_phi, Float_t radius=1.6f);

  virtual bool Pass(const AliFemtoPair *);
  virtual bool Pass(AliFemtoPair *);

  virtual AliFemtoString Report() {return "";};
  virtual TList* ListSettings();
  virtual TList* AppendSettings(TList *, const TString &prefix="");

  void SetR(Float_t r);
  void SetMinEta(Float_t eta);
  void SetMinPhi(Float_t phi);

  static Float_t CalculateDEta(const AliFemtoThreeVector& a, const AliFemtoThreeVector& b);

  static Float_t CalculateDPhiStar(const AliFemtoThreeVector& a, const AliFemtoThreeVector& b, const Double_t minRad);
  static Float_t CalculateDEtaStar(const AliFemtoThreeVector& a, const AliFemtoThreeVector& b, const Double_t minRad);

protected:

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

  bool passes = AliFemtoShareQualityPairCut::Pass(pair)
             && fDeltaEtaMin <= fabs(CalculateDEta(p1, p2))
             && fDeltaPhiMin <= fabs(CalculateDPhiStar(p1, p2, fR));

  // cout << "> " << passes << " " << CalculateDEta(p1, p2) << " " << CalculateDPhiStar(p1, p2, fR) << "\n";

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
  const AliFemtoThreeVector& a,
  const AliFemtoThreeVector& b,
  const Double_t radius_in_meters)
{
  // phi shift at radius R in magnetic field B, with charge and momentum q & pT is:
  //    φ = arcsin( q * B * R / 2 pT )
  //
  // Unit analysis:
  //   q * B * R : [Coulomb][Tesla][Meter] = 1.60218e-19 [e][Tesla][Meter]
  //   pT : [Joule][Second]/[Meter] = 1.60218e-10 [GeV][Second]/[Meter] = 1.60218e-10 / c [GeV/c]
  //
  //  q * B * R / pT = (1.60218e-19 * c / 1.60218e-10) [e] [Tesla] [Meters] / [GeV/c]
  //                 ~ 0.3 [e] [Tesla] [Meters] / [GeV/c]
  //

  const Double_t unit_factor = 0.299792458 / 2.0;
  const Double_t b_field = 0.5006670488586,
                  charge = 1.0;
  // ((AliAODInputHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->GetEvent()->GetMagneticField();

  const Double_t phi_a = a.Phi(),
                 phi_b = b.Phi(),

                 shift_a = TMath::ASin(- unit_factor * charge * b_field * radius_in_meters / a.Perp()),
                 shift_b = TMath::ASin(- unit_factor * charge * b_field * radius_in_meters / b.Perp()),

                 delta_phi_star = (phi_a + shift_a) - (phi_b + shift_b);

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


#endif
