#ifndef Lifetimes_HyperTriton2Body_h
#define Lifetimes_HyperTriton2Body_h

#include <cmath>
#include <cstdlib>
#include <limits>
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "Rtypes.h"
#include "Utils.h"
#include <TMath.h>
#include <vector>

namespace Lifetimes
{

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;
constexpr double Sq(double x) { return x * x; }

class HyperTriton2Body
{
public:
  float GetV0radius() const { return std::abs(fV0radius); }
  float GetV0pt() const { return std::abs(fV0pt); }
  float GetV0eta() const { return fV0eta; }
  float GetCandidateInvMass() const { return fInvMass; }
  float GetDistOverP() const { return fDistOverTotMom; }
  Double32_t GetV0CosPA() const { return fV0CosPA; }
  Double32_t GetV0chi2() const { return fChi2V0; }
  Double32_t GetNegProngPvDCA() const { return fDcaNeg2PrimaryVertex; };
  Double32_t GetPosProngPvDCA() const { return fDcaPos2PrimaryVertex; }
  Double32_t GetProngsDCA() const { return fDcaV0daughters; };
  Double32_t GetArmenterosAlpha() const { return fV0armAlpha; }
  Double32_t GetArmenterosPt() const { return fV0armPt; }
  unsigned char GetLeastNumberOfXedRows() const { return fLeastNxedRows; }
  Double32_t GetLeastXedRowsOverFindable() const { return fLeastXedOverFindable; }
  Double32_t GetMaxChi2perCluster() const { return fMaxChi2PerCluster; }
  Double32_t GetNegProngTPCnsigmaPion() const { return fNsigmaPionNeg; }
  Double32_t GetNegProngTPCnsigmaHe3() const { return fNsigmaHe3Neg; }
  Double32_t GetPosProngTPCnsigmaPion() const { return fNsigmaPionPos; }
  Double32_t GetPosProngTPCnsigmaHe3() const { return fNsigmaHe3Pos; }
  Double32_t GetNegProngEta() const { return fEtaNeg; }
  Double32_t GetPosProngEta() const { return fEtaPos; }
  Double32_t GetNegProngPt() const { return fPtNeg; }
  Double32_t GetPosProngPt() const { return fPtPos; }
  Double32_t GetNegProngPhi() const { return fPhiNeg; }
  Double32_t GetPosProngPhi() const { return fPhiPos; }
  bool IsCowboy() const { return fFlags & kCowboySailor; }
  bool IsLikeSign() const { return fV0radius < 0.; } //TODO: switch to signbit with ROOT6
  bool IsFake() const { return fV0pt < 0.; }         //TODO: switch to signbit with ROOT6
  bool NegativeProngHasTOF() const { return fFlags & kNegativeTOF; }
  bool PositiveProngHasTOF() const { return fFlags & kPositiveTOF; }
  bool OneProngHasTOF() const { return NegativeProngHasTOF() || PositiveProngHasTOF(); }
  bool NegativeProngHasITSrefit() const { return fITSInfo & (1 << 7); }
  bool PositiveProngHasITSrefit() const { return fITSInfo & (1 << 6); }
  bool NegativeProngHasSPDcluster() const { return fITSInfo & (1 << 5); }
  bool PositiveProngHasSPDcluster() const { return fITSInfo & (1 << 4); }
  int GetLeastNumberOfITSclusters() const { return fITSInfo & 0x07; }

  void SetV0radiusAndLikeSign(float r, bool ls = false) { fV0radius = ls ? -r : r; }
  void SetV0ptAndFake(float pt, bool fake) { fV0pt = fake ? -pt : pt; }
  void SetV0eta(float eta) { fV0eta = eta; }
  void SetInvMass(float m) { fInvMass = m; }
  void SetDistOverP(float distOnP) { fDistOverTotMom = distOnP; }
  void SetV0CosPA(float cospa) { fV0CosPA = cospa; }
  void SetV0Chi2(float chi2) { fChi2V0 = chi2; }
  void SetProngsPvDCA(float posD, float negD);
  void SetProngsDCA(float dca) { fDcaV0daughters = dca; };
  void SetArmenterosVariables(float alpha, float pt);
  void SetLeastNumberOfXedRows(unsigned char xedrows) { fLeastNxedRows = xedrows; }
  void SetLeastXedRowsOverFindable(float ratio) { fLeastXedOverFindable = ratio; }
  void SetMaxChi2perCluster(float chi2) { fMaxChi2PerCluster = chi2; }
  void SetProngsTPCnsigmas(float pPi, float pP, float nPi, float nP);
  void SetProngsEtaTOF(float posEta, float negEta);
  void SetProngsPt(float posPt, float negPt);
  void SetProngsPhi(float posPhi, float negPhi);
  void SetITSinformation(bool, bool, bool, bool, int);
  void SetCowboyAndSailor(bool cs) { fFlags = flipBits(fFlags, static_cast<unsigned char>(kCowboySailor), cs); }
  void SetTOFbits(bool pTOF, bool nTOF);
  void SetOptimalParameters(bool opt) { fFlags = flipBits(fFlags, static_cast<unsigned char>(kOptimalParams), opt); }
  static HyperTriton2Body FillHyperTriton2Body(AliESDv0 *v0, AliESDtrack *pTrack, AliESDtrack *nTrack, float nsigmaposhe3, float nsigmaneghe3, float nsigmapospion, float nsigmanegpion, float magneticField, double primaryVertex[3], bool fake);
  LVector_t GetV0LorentzVector(AliESDtrack *negTrack, AliESDtrack *posTrack, double alpha);

private:
  enum Flags
  {
    kCowboySailor = 1,
    kNegativeTOF = 1 << 1,
    kPositiveTOF = 1 << 2,
    kOptimalParams = 1 << 3
  };

  float fV0radius;                  // V0 decay vertex radius (negarive -> LikeSign V0)
  float fV0pt;                      // V0 transverse momentum (in MC if negative -> fake V0)
  float fV0eta;                     // V0 pseudorapidity
  float fInvMass;                   // Invariant mass for the candidate
  float fDistOverTotMom;            // L/p
  Double32_t fV0CosPA;              //[0.9,1.0,16] V0 cosine of pointing angle
  Double32_t fChi2V0;               //[0.0,10.24,8] V0 fit chi2
  Double32_t fDcaNeg2PrimaryVertex; //[0.0,0.256,8] DCA of the negative prong to the PV
  Double32_t fDcaPos2PrimaryVertex; //[0.0,0.256,8]  DCA of the positive prong to the PV
  Double32_t fDcaV0daughters;       //[0.0,2.56,8] DCA between the two prongs
  Double32_t fV0armAlpha;           //[-1.28,1.28,8] Armenteros alpha
  Double32_t fV0armPt;              //[0.0,0.256,8] Armenteros pt
  Double32_t fLeastXedOverFindable; //[0.36,1.0,7] Min xed roads/findable clusters
  Double32_t fMaxChi2PerCluster;    //[0,6.4,7] Max chi2 per cluster in TPC
  Double32_t fNsigmaPionPos;        //[-8.0,8.0,5] number of sigmas TPC pion for the positive prong
  Double32_t fNsigmaHe3Pos;         //[-8.0,8.0,5] number of sigmas TPC 3He for the positive prong
  Double32_t fNsigmaPionNeg;        //[-8.0,8.0,5] number of sigmas TPC pion for the negative prong
  Double32_t fNsigmaHe3Neg;         //[-8.0,8.0,5] number of sigmas TPC 3He for the negative prong
  Double32_t fPtPos;                //[0,12.8,8]
  Double32_t fPtNeg;                //[0,12.8,8]
  Double32_t fPhiPos;               //[0,2*pi,8]
  Double32_t fPhiNeg;               //[0,2*pi,8]
  Double32_t fEtaPos;               //[-1.0,1.0,7] Pseudorapidity of the positive prong. MSB is the TOF bit.
  Double32_t fEtaNeg;               //[-1.0,1.0,7] Pseudorapidity of the negative prong. MSB is the TOF bit.
  unsigned char fLeastNxedRows;     // Min number of xed roads
  unsigned char fITSInfo;           // Starting from the MSB: kITSrefit for neg and pos, kSPDany for neg and pos, least number of ITS clusters (last 4 bits)
  unsigned char fFlags;             // Cowboy&Saylor, TOF bits for neg and pos, optimal tracking parameters
};
/// Cowboy&saylor, 2xTOF bit

inline void HyperTriton2Body::SetProngsPvDCA(float posD, float negD)
{
  fDcaPos2PrimaryVertex = posD;
  fDcaNeg2PrimaryVertex = negD;
}

inline void HyperTriton2Body::SetArmenterosVariables(float alpha, float pt)
{
  fV0armAlpha = alpha;
  fV0armPt = pt;
}

inline void HyperTriton2Body::SetProngsTPCnsigmas(float pPi, float pHe, float nPi,
                                                  float nHe)
{
  fNsigmaPionPos = pPi;
  fNsigmaHe3Pos = pHe;
  fNsigmaPionNeg = nPi;
  fNsigmaHe3Neg = nHe;
}

inline void HyperTriton2Body::SetProngsEtaTOF(float posEta, float negEta)
{
  fEtaNeg = negEta;
  fEtaPos = posEta;
}

inline void HyperTriton2Body::SetProngsPt(float posPt, float negPt)
{
  fPtNeg = negPt;
  fPtPos = posPt;
}

inline void HyperTriton2Body::SetProngsPhi(float posPhi, float negPhi)
{
  fPhiNeg = negPhi;
  fPhiPos = posPhi;
}

inline void HyperTriton2Body::SetITSinformation(bool negRefit, bool posRefit, bool negSPD, bool posSPD, int nITS)
{
  fITSInfo = 0u;
  if (negRefit)
    fITSInfo |= 1 << 7;
  if (posRefit)
    fITSInfo |= 1 << 6;
  if (negSPD)
    fITSInfo |= 1 << 5;
  if (posSPD)
    fITSInfo |= 1 << 4;
  fITSInfo |= 0x07 & nITS;
}

inline LVector_t HyperTriton2Body::GetV0LorentzVector(AliESDtrack *negTrack, AliESDtrack *posTrack, double alpha)
{
  constexpr AliPID::EParticleType children[1][2]{
      {AliPID::kHe3, AliPID::kPion},
  };
  int posIndex = int(alpha < 0);
  int negIndex = int(alpha >= 0);
  double posMass = AliPID::ParticleMass(children[0][posIndex]);
  double negMass = AliPID::ParticleMass(children[0][negIndex]);
  int posCharge = AliPID::ParticleCharge(children[0][posIndex]);
  int negCharge = AliPID::ParticleCharge(children[0][negIndex]);
  double posMom[3], negMom[3];
  posTrack->GetPxPyPz(posMom);
  negTrack->GetPxPyPz(negMom);
  LVector_t posLvec{posMom[0] * posCharge, posMom[1] * posCharge, posCharge * posMom[2], posMass};
  LVector_t negLvec{negMom[0] * negCharge, negMom[1] * negCharge, negCharge * negMom[2], negMass};

  posLvec += negLvec;
  return posLvec;
}

inline void HyperTriton2Body::SetTOFbits(bool pTOF, bool nTOF)
{
  flipBits(fFlags, static_cast<unsigned char>(kPositiveTOF), pTOF);
  flipBits(fFlags, static_cast<unsigned char>(kNegativeTOF), nTOF);
}

inline HyperTriton2Body HyperTriton2Body::FillHyperTriton2Body(AliESDv0 *v0, AliESDtrack *pTrack, AliESDtrack *nTrack, float nsigmaposhe3, float nsigmaneghe3, float nsigmapospion, float nsigmanegpion, float magneticField, double primaryVertex[3], bool fake)
{

  HyperTriton2Body miniHyper;
  double decayVtx[3];
  v0->GetXYZ(decayVtx[0], decayVtx[1], decayVtx[2]);
  double v0Radius = std::hypot(decayVtx[0], decayVtx[1]);
  double alpha = v0->AlphaV0();
  auto lvector = miniHyper.GetV0LorentzVector(nTrack, pTrack, alpha);
  double mass = lvector.M();
  double v0Pt = lvector.Pt();
  double lV0TotalMomentum = lvector.P();
  float distOverP = std::sqrt(Sq(decayVtx[0] - primaryVertex[0]) +
                              Sq(decayVtx[1] - primaryVertex[1]) +
                              Sq(decayVtx[2] - primaryVertex[2])) /
                    (lV0TotalMomentum + 1e-16);
  double PtPos;
  double PtNeg;
  if (alpha < 0)
  {
    PtPos = pTrack->Pt();
    PtNeg = 2 * (nTrack->Pt());
  }
  else
  {
    PtPos = 2 * (pTrack->Pt());
    PtNeg = nTrack->Pt();
  }
  unsigned char posXedRows = pTrack->GetTPCClusterInfo(2, 1);
  unsigned char negXedRows = nTrack->GetTPCClusterInfo(2, 1);
  float posChi2PerCluster =
      pTrack->GetTPCchi2() / (pTrack->GetTPCNcls() + 1.e-16);
  float negChi2PerCluster =
      nTrack->GetTPCchi2() / (nTrack->GetTPCNcls() + 1.e-16);
  float posXedRowsOverFindable = float(posXedRows) / pTrack->GetTPCNclsF();
  float negXedRowsOverFindable = float(negXedRows) / nTrack->GetTPCNclsF();
  unsigned char minXedRows =
      posXedRows < negXedRows ? posXedRows : negXedRows;
  float minXedRowsOverFindable =
      posXedRowsOverFindable < negXedRowsOverFindable
          ? posXedRowsOverFindable
          : negXedRowsOverFindable;
  float maxChi2PerCluster = posChi2PerCluster > negChi2PerCluster
                                ? posChi2PerCluster
                                : negChi2PerCluster;

  double cosPA = v0->GetV0CosineOfPointingAngle(primaryVertex[0], primaryVertex[1], primaryVertex[2]);

  double dcaPosToPrimVertex = std::abs(
      pTrack->GetD(primaryVertex[0], primaryVertex[1], magneticField));

  double dcaNegToPrimVertex = std::abs(
      nTrack->GetD(primaryVertex[0], primaryVertex[1], magneticField));

  bool negTOF = nTrack->GetTOFsignal() * 1.e-3 < 100; // in ns, loose cut on TOF beta (<~0.2)
  bool posTOF = pTrack->GetTOFsignal() * 1.e-3 < 100; // in ns, loose cut on TOF beta (<~0.2)

  bool posITSrefit = pTrack->GetStatus() & AliESDtrack::kITSrefit;
  bool negITSrefit = nTrack->GetStatus() & AliESDtrack::kITSrefit;
  bool posSPDany = pTrack->HasPointOnITSLayer(0) || pTrack->HasPointOnITSLayer(1);
  bool negSPDany = nTrack->HasPointOnITSLayer(0) || nTrack->HasPointOnITSLayer(1);
  int ITSnCl = (nTrack->GetITSclusters(0) > pTrack->GetITSclusters(0)) ? pTrack->GetITSclusters(0) : nTrack->GetITSclusters(0);

  double momPos[3];
  v0->GetPPxPyPz(momPos[0], momPos[1], momPos[2]);
  double momNeg[3];
  v0->GetNPxPyPz(momNeg[0], momNeg[1], momNeg[2]);

  double lVecProd = momPos[0] * momNeg[1] - momPos[1] * momNeg[0];
  bool isCowboy = lVecProd * magneticField < 0;

  miniHyper.SetV0ptAndFake(v0Pt, fake);
  miniHyper.SetV0eta(v0->Eta());
  miniHyper.SetLeastNumberOfXedRows(minXedRows);
  miniHyper.SetDistOverP(distOverP);
  miniHyper.SetInvMass(mass);
  miniHyper.SetArmenterosVariables(v0->AlphaV0(), v0->PtArmV0());
  miniHyper.SetV0CosPA(cosPA);
  miniHyper.SetV0Chi2(v0->GetChi2V0());
  miniHyper.SetProngsPt(PtPos, PtNeg);
  miniHyper.SetProngsPhi(pTrack->Phi(), nTrack->Phi());
  miniHyper.SetProngsDCA(v0->GetDcaV0Daughters());
  miniHyper.SetProngsPvDCA(dcaPosToPrimVertex, dcaNegToPrimVertex);
  miniHyper.SetV0radiusAndLikeSign(v0Radius);
  miniHyper.SetLeastXedRowsOverFindable(minXedRowsOverFindable);
  miniHyper.SetMaxChi2perCluster(maxChi2PerCluster);
  miniHyper.SetProngsEtaTOF(pTrack->Eta(), nTrack->Eta());
  miniHyper.SetProngsTPCnsigmas(nsigmapospion, nsigmaposhe3,
                                nsigmanegpion, nsigmaneghe3);
  miniHyper.SetITSinformation(negITSrefit, posITSrefit, negSPDany, posSPDany, ITSnCl);
  miniHyper.SetTOFbits(posTOF, negTOF);
  miniHyper.SetCowboyAndSailor(isCowboy);
  return miniHyper;
}

} // namespace Lifetimes

#endif
