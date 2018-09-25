#ifndef Lifetimes_MiniV0_h
#define Lifetimes_MiniV0_h

#include <cmath>
#include <cstdlib>
#include <limits>

#include "Utils.h"

namespace Lifetimes {

template<int n>
class MiniV0 {
 public:
  float GetV0radius() const { return std::abs(fV0radius); }
  float GetV0pt() const { return std::abs(fV0pt); }
  float GetV0eta() const { return fV0eta; }
  float GetCandidateInvMass(int i) const { return fInvMass[i]; }
  float GetDistOverP() const { return fDistOverTotMom; }
  float GetV0CosPA() const;
  float GetV0chi2() const;
  float GetNegProngPvDCA() const;
  float GetPosProngPvDCA() const;
  float GetProngsDCA() const;
  float GetArmenterosAlpha() const;
  float GetArmenterosPt() const;
  unsigned char GetLeastNumberOfXedRows() const { return fLeastNxedRows; }
  float GetLeastXedRowsOverFindable() const;
  float GetMaxChi2perCluster() const;
  float GetNegProngTPCnsigmaPion() const;
  float GetNegProngTPCnsigmaProton() const;
  float GetPosProngTPCnsigmaPion() const;
  float GetPosProngTPCnsigmaProton() const;
  float GetNegProngEta() const;
  float GetPosProngEta() const;
  bool IsCowboy() const { return fChi2V0 & (1 << 7); }
  bool IsLikeSign() const { return fV0radius < 0.; } //TODO: switch to signbit with ROOT6
  bool IsFake() const { return fV0pt < 0.; }         //TODO: switch to signbit with ROOT6

  void SetV0radiusAndLikeSign(float r, bool ls = false) { fV0radius = ls ? -r : r; }
  void SetV0ptAndFake(float pt, bool fake) { fV0pt = fake ? -pt : pt; }
  void SetV0eta(float eta) { fV0eta = eta; }
  void SetInvMass(int i, float m) { fInvMass[i] = m; }
  void SetDistOverP(float distOnP) { fDistOverTotMom = distOnP; }
  void SetV0CosPA(float cospa);
  void SetV0Chi2andCowBoy(float chi2, bool cow);
  void SetProngsPvDCA(float posD, float negD);
  void SetProngsDCA(float dca);
  void SetArmenterosVariables(float alpha, float pt);
  void SetLeastNumberOfXedRows(unsigned char xedrows);
  void SetLeastXedRowsOverFindable(float ratio);
  void SetMaxChi2perCluster(float chi2);
  void SetProngsTPCnsigmas(float pPi, float pP, float nPi, float nP);
  void SetProngsEta(float posEta, float negEta);

  static const int fgkV0cosPA_n;
  static const float fgkV0cosPA_f;
  static const float fgkV0cosPA_l;
  static const float fgkV0cosPA_w;

  static const int fgkV0chi2_n;
  static const float fgkV0chi2_f;
  static const float fgkV0chi2_l;
  static const float fgkV0chi2_w;

  static const int fgkDCAProng2PV_n;
  static const float fgkDCAProng2PV_f;
  static const float fgkDCAProng2PV_l;
  static const float fgkDCAProng2PV_w;

  static const int fgkDCAProngs_n;
  static const float fgkDCAProngs_f;
  static const float fgkDCAProngs_l;
  static const float fgkDCAProngs_w;

  static const int fgkArmAlpha_n;
  static const float fgkArmAlpha_f;
  static const float fgkArmAlpha_l;
  static const float fgkArmAlpha_w;

  static const int fgkArmPt_n;
  static const float fgkArmPt_f;
  static const float fgkArmPt_l;
  static const float fgkArmPt_w;

  static const int fgkXedOverFindable_n;
  static const float fgkXedOverFindable_f;
  static const float fgkXedOverFindable_l;
  static const float fgkXedOverFindable_w;

  static const int fgkChi2xCluster_n;
  static const float fgkChi2xCluster_f;
  static const float fgkChi2xCluster_l;
  static const float fgkChi2xCluster_w;

  static const int fgkTPCsigma_n;
  static const float fgkTPCsigma_f;
  static const float fgkTPCsigma_l;
  static const float fgkTPCsigma_w;

  static const int fgkEta_n;
  static const float fgkEta_f;
  static const float fgkEta_l;
  static const float fgkEta_w;

 private:
  float fV0radius;                      // V0 decay vertex radius (negarive -> LikeSign V0)
  float fV0pt;                          // V0 transverse momentum (in MC if negative -> fake V0)
  float fV0eta;                         // V0 pseudorapidity
  float fInvMass[n];                    // Invariant mass for the candidate
  float fDistOverTotMom;                // L/p
  unsigned short fV0CosPA;              // V0 cosine of pointing angle
  unsigned char fChi2V0;                // V0 fit chi2
  unsigned char fDcaNeg2PrimaryVertex;  // DCA of the negative prong to the PV
  unsigned char fDcaPos2PrimaryVertex;  // DCA of the positive prong to the PV
  unsigned char fDcaV0daughters;        // DCA between the two prongs
  unsigned char fV0armAlpha;            // Armenteros alpha
  unsigned char fV0armPt;               // Armenteros pt
  unsigned char fLeastNxedRows;         // Min number of xed roads
  unsigned char fLeastXedOverFindable;  // Min xed roads/findable clusters
  unsigned char fMaxChi2PerCluster;     // Max chi2 per cluster in TPC
  unsigned char fNsigmaPos;  // # sigma TPC pion/proton for the positive prong
  unsigned char fNsigmaNeg;  // # sigma TPC pion/proton for the negative prong
  unsigned char fEtaPos;     // Pseudorapidity of the positive prong
  unsigned char fEtaNeg;     // Pseudorapidity of the negative prong
};

template<int n>
inline float MiniV0<n>::GetV0CosPA() const {
  return getBinCenter(fV0CosPA, fgkV0cosPA_w, fgkV0cosPA_f);
}
template<int n>
inline float MiniV0<n>::GetV0chi2() const {
  return getBinCenter(fChi2V0 & 0x7F, fgkV0chi2_w, fgkV0chi2_f);
}
template<int n>
inline float MiniV0<n>::GetNegProngPvDCA() const {
  return getBinCenter(fDcaNeg2PrimaryVertex, fgkDCAProng2PV_w, fgkDCAProng2PV_f,
                      false, true);
}
template<int n>
inline float MiniV0<n>::GetPosProngPvDCA() const {
  return getBinCenter(fDcaPos2PrimaryVertex, fgkDCAProng2PV_w, fgkDCAProng2PV_f,
                      false, true);
}
template<int n>
inline float MiniV0<n>::GetProngsDCA() const {
  return getBinCenter(fDcaV0daughters, fgkDCAProngs_w, fgkDCAProngs_f,
                      false, true);
}
template<int n>
inline float MiniV0<n>::GetArmenterosAlpha() const {
  return getBinCenter(fV0armAlpha, fgkArmAlpha_w, fgkArmAlpha_f,
                      true, true);
}
template<int n>
inline float MiniV0<n>::GetArmenterosPt() const {
  return getBinCenter(fV0armPt, fgkArmPt_w, fgkArmPt_f, false, true);
}
template<int n>
inline float MiniV0<n>::GetLeastXedRowsOverFindable() const {
  return getBinCenter(fLeastXedOverFindable, fgkXedOverFindable_w,
                      fgkXedOverFindable_f);
}
template<int n>
inline float MiniV0<n>::GetMaxChi2perCluster() const {
  return getBinCenter(fMaxChi2PerCluster, fgkChi2xCluster_w, fgkChi2xCluster_f,
                      false, true);
}

template<int n>
inline float MiniV0<n>::GetNegProngEta() const {
  return getBinCenter(fEtaNeg, fgkEta_w, fgkEta_f, true, true);
}
template<int n>
inline float MiniV0<n>::GetPosProngEta() const {
  return getBinCenter(fEtaPos, fgkEta_w, fgkEta_f, true, true);
}

template<int n>
inline float MiniV0<n>::GetNegProngTPCnsigmaPion() const {
  return getBinCenter((fNsigmaNeg & 0xF0) >> 4, fgkTPCsigma_w, fgkTPCsigma_f,
                      false, true);
}
template<int n>
inline float MiniV0<n>::GetNegProngTPCnsigmaProton() const {
  return getBinCenter((fNsigmaNeg & 0x0F), fgkTPCsigma_w, fgkTPCsigma_f,
                      false, true);
}
template<int n>
inline float MiniV0<n>::GetPosProngTPCnsigmaPion() const {
  return getBinCenter((fNsigmaPos & 0xF0) >> 4, fgkTPCsigma_w, fgkTPCsigma_f,
                      false, true);
}
template<int n>
inline float MiniV0<n>::GetPosProngTPCnsigmaProton() const {
  return getBinCenter((fNsigmaPos & 0x0F), fgkTPCsigma_w, fgkTPCsigma_f,
                      false, true);
}

template<int n>
inline void MiniV0<n>::SetV0CosPA(float cospa) {
  fV0CosPA = getBinnedValue<unsigned short>(cospa, fgkV0cosPA_w, fgkV0cosPA_f,
                                            fgkV0cosPA_l);
}

template<int n>
inline void MiniV0<n>::SetV0Chi2andCowBoy(float chi2, bool cow) {
  fChi2V0 = getBinnedValue<char>(chi2, fgkV0chi2_w, fgkV0chi2_l);
  if (cow) fChi2V0 += 1 << 7;
}

template<int n>
inline void MiniV0<n>::SetProngsPvDCA(float posD, float negD) {
  fDcaPos2PrimaryVertex =
      getBinnedValue<unsigned char>(posD, fgkDCAProng2PV_w, fgkDCAProng2PV_l);
  fDcaNeg2PrimaryVertex =
      getBinnedValue<unsigned char>(negD, fgkDCAProng2PV_w, fgkDCAProng2PV_l);
}

template<int n>
inline void MiniV0<n>::SetProngsDCA(float dca) {
  fDcaV0daughters =
      getBinnedValue<unsigned char>(dca, fgkDCAProngs_w, fgkDCAProngs_l);
}

template<int n>
inline void MiniV0<n>::SetArmenterosVariables(float alpha, float pt) {
  fV0armAlpha = getBinnedValue<unsigned char>(alpha, fgkArmAlpha_w,
                                              fgkArmAlpha_f, fgkArmAlpha_l);
  fV0armPt = getBinnedValue<unsigned char>(pt, fgkArmPt_w, fgkArmPt_l);
}

template<int n>
inline void MiniV0<n>::SetLeastNumberOfXedRows(unsigned char xedrows) {
  fLeastNxedRows = xedrows;
}

template<int n>
inline void MiniV0<n>::SetLeastXedRowsOverFindable(float ratio) {
  fLeastXedOverFindable =
      getBinnedValue<unsigned char>(ratio, fgkXedOverFindable_w, fgkXedOverFindable_l);
}

template<int n>
inline void MiniV0<n>::SetMaxChi2perCluster(float chi2) {
  fMaxChi2PerCluster = getBinnedValue<unsigned char>(chi2, fgkChi2xCluster_w);
}

template<int n>
inline void MiniV0<n>::SetProngsTPCnsigmas(float pPi, float pP, float nPi,
                                        float nP) {
  fNsigmaNeg = (std::abs(nP) > fgkTPCsigma_l)
                   ? 0x0F
                   : getBinnedValue<unsigned char>(std::abs(nP), fgkTPCsigma_w);
  fNsigmaNeg +=
      ((std::abs(nPi) > fgkTPCsigma_l)
           ? 0x0F
           : getBinnedValue<unsigned char>(std::abs(nPi), fgkTPCsigma_w))
      << 4;
  fNsigmaPos = (std::abs(pP) > fgkTPCsigma_l)
                   ? 0x0F
                   : getBinnedValue<unsigned char>(std::abs(pP), fgkTPCsigma_w);
  fNsigmaPos +=
      ((std::abs(pPi) > fgkTPCsigma_l)
           ? 0x0F
           : getBinnedValue<unsigned char>(std::abs(pPi), fgkTPCsigma_w))
      << 4;
}

template<int n>
inline void MiniV0<n>::SetProngsEta(float posEta, float negEta) {
  fEtaNeg = getBinnedValue<unsigned char>(negEta, fgkEta_w, fgkEta_f, fgkEta_l);
  fEtaPos = getBinnedValue<unsigned char>(posEta, fgkEta_w, fgkEta_f, fgkEta_l);
}

template<int n> const int MiniV0<n>::fgkV0cosPA_n = 50000;
template<int n> const float MiniV0<n>::fgkV0cosPA_f = 0.9f;
template<int n> const float MiniV0<n>::fgkV0cosPA_l = 1.f;
template<int n> const float MiniV0<n>::fgkV0cosPA_w = \
    (MiniV0<n>::fgkV0cosPA_l - MiniV0<n>::fgkV0cosPA_f) / \
    MiniV0<n>::fgkV0cosPA_n;

template<int n> const int MiniV0<n>::fgkV0chi2_n = 100;
template<int n> const float MiniV0<n>::fgkV0chi2_f = 0.f;
template<int n> const float MiniV0<n>::fgkV0chi2_l = 10.f;
template<int n> const float MiniV0<n>::fgkV0chi2_w = \
    (MiniV0<n>::fgkV0chi2_l - MiniV0<n>::fgkV0chi2_f) / MiniV0<n>::fgkV0chi2_n;

template<int n> const int MiniV0<n>::fgkDCAProng2PV_n = 250;
template<int n> const float MiniV0<n>::fgkDCAProng2PV_f = 0.f;
template<int n> const float MiniV0<n>::fgkDCAProng2PV_l = 0.25f;
template<int n> const float MiniV0<n>::fgkDCAProng2PV_w = \
    (MiniV0<n>::fgkDCAProng2PV_l - MiniV0<n>::fgkDCAProng2PV_f) / \
    MiniV0<n>::fgkDCAProng2PV_n;

template<int n> const int MiniV0<n>::fgkDCAProngs_n = 250;
template<int n> const float MiniV0<n>::fgkDCAProngs_f = 0.f;
template<int n> const float MiniV0<n>::fgkDCAProngs_l = 2.f;
template<int n> const float MiniV0<n>::fgkDCAProngs_w = \
    (MiniV0<n>::fgkDCAProngs_l - MiniV0<n>::fgkDCAProngs_f) / \
    MiniV0<n>::fgkDCAProngs_n;

template<int n> const int MiniV0<n>::fgkArmAlpha_n = 250;
template<int n> const float MiniV0<n>::fgkArmAlpha_f = -1.f;
template<int n> const float MiniV0<n>::fgkArmAlpha_l = 1.f;
template<int n> const float MiniV0<n>::fgkArmAlpha_w = \
    (MiniV0<n>::fgkArmAlpha_l - MiniV0<n>::fgkArmAlpha_f) / \
    MiniV0<n>::fgkArmAlpha_n;

template<int n> const int MiniV0<n>::fgkArmPt_n = 254;
template<int n> const float MiniV0<n>::fgkArmPt_f = 0.f;
template<int n> const float MiniV0<n>::fgkArmPt_l = 0.254f;
template<int n> const float MiniV0<n>::fgkArmPt_w = \
    (MiniV0<n>::fgkArmPt_l - MiniV0<n>::fgkArmPt_f) / MiniV0<n>::fgkArmPt_n;

template<int n> const int MiniV0<n>::fgkXedOverFindable_n = 250;
template<int n> const float MiniV0<n>::fgkXedOverFindable_f = 0.f;
template<int n> const float MiniV0<n>::fgkXedOverFindable_l = 1.f;
template<int n> const float MiniV0<n>::fgkXedOverFindable_w = \
    (MiniV0<n>::fgkXedOverFindable_l - MiniV0<n>::fgkXedOverFindable_f) / \
    MiniV0<n>::fgkXedOverFindable_n;

template<int n> const int MiniV0<n>::fgkChi2xCluster_n = 250;
template<int n> const float MiniV0<n>::fgkChi2xCluster_f = 0.f;
template<int n> const float MiniV0<n>::fgkChi2xCluster_l = 10.f;
template<int n> const float MiniV0<n>::fgkChi2xCluster_w = \
    (MiniV0<n>::fgkChi2xCluster_l - MiniV0<n>::fgkChi2xCluster_f) / \
    MiniV0<n>::fgkChi2xCluster_n;

template<int n> const int MiniV0<n>::fgkTPCsigma_n = 12;
template<int n> const float MiniV0<n>::fgkTPCsigma_f = 0.f;
template<int n> const float MiniV0<n>::fgkTPCsigma_l = 6.f;
template<int n> const float MiniV0<n>::fgkTPCsigma_w = \
    (MiniV0<n>::fgkTPCsigma_l - MiniV0<n>::fgkTPCsigma_f) / \
    MiniV0<n>::fgkTPCsigma_n;

template<int n> const int MiniV0<n>::fgkEta_n = 200;
template<int n> const float MiniV0<n>::fgkEta_f = -1.f;
template<int n> const float MiniV0<n>::fgkEta_l = 1.f;
template<int n> const float MiniV0<n>::fgkEta_w = \
    (MiniV0<n>::fgkEta_l - MiniV0<n>::fgkEta_f) / MiniV0<n>::fgkEta_n;

}  // namespace Lifetimes

#endif
