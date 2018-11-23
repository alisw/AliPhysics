#ifndef Lifetimes_HyperTriton2Body_h
#define Lifetimes_HyperTriton2Body_h

#include <cmath>
#include <cstdlib>
#include <limits>

#include "Rtypes.h"
#include "Utils.h"


namespace Lifetimes {

class HyperTriton2Body {
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
  Double32_t GetProngsDCA() const {return fDcaV0daughters;};
  Double32_t GetArmenterosAlpha() const{ return fV0armAlpha;}
  Double32_t GetArmenterosPt() const { return fV0armPt;}
  unsigned char GetLeastNumberOfXedRows() const { return fLeastNxedRows; }
  Double32_t GetLeastXedRowsOverFindable() const { return fLeastXedOverFindable;}
  Double32_t GetMaxChi2perCluster() const{ return fMaxChi2PerCluster;}
  Double32_t GetNegProngTPCnsigmaPion() const { return fNsigmaPionNeg;}
  Double32_t GetNegProngTPCnsigmaHe3() const { return fNsigmaHe3Neg;}
  Double32_t GetPosProngTPCnsigmaPion() const {return fNsigmaPionPos;}
  Double32_t GetPosProngTPCnsigmaHe3() const {return fNsigmaHe3Pos;}
  Double32_t GetNegProngEta() const { return fEtaNeg;}
  Double32_t GetPosProngEta() const { return fEtaPos;}
  Double32_t GetNegProngPt() const { return fPtNeg;}
  Double32_t GetPosProngPt() const { return fPtPos;}
  Double32_t GetNegProngPhi() const { return fPhiNeg;}
  Double32_t GetPosProngPhi() const { return fPhiPos;}  
  bool IsCowboy() const { return fChi2V0 & (1 << 7); }
  bool IsLikeSign() const { return fV0radius < 0.; } //TODO: switch to signbit with ROOT6
  bool IsFake() const { return fV0pt < 0.; }         //TODO: switch to signbit with ROOT6
  bool NegativeProngHasTOF() const { return fEtaNeg & (1 << 7); }
  bool PositiveProngHasTOF() const { return fEtaPos & (1 << 7); }
  bool OneProngHasTOF() const { return NegativeProngHasTOF() || PositiveProngHasTOF(); }
  bool NegativeProngHasITSrefit() const { return fITSInfo & (1 << 7); }
  bool PositiveProngHasITSrefit() const { return fITSInfo & (1 << 6); }
  bool NegativeProngHasSPDcluster() const { return fITSInfo & (1 << 5); }
  bool PositiveProngHasSPDcluster() const { return fITSInfo & (1 << 4); }
  int  GetLeastNumberOfITSclusters() const { return fITSInfo & 0x07; }
 
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
  void HyperTriton2Body::SetProngsPt(float posPt, float negPt);
  void HyperTriton2Body::SetProngsPhi(float posPhi, float negPhi);
  void SetITSinformation(bool, bool, bool, bool, int);
  void SetCowboyAndSailor(bool cs) { fFlags = flipBits(fFlags, kCowboySailor, cs); }
  void SetTOFbits(bool pTOF, bool nTOF);
  void SetOptimalParameters(bool opt) { fFlags = flipBits(fFlags, kOptimalParams, opt); }

 private:
  enum {
    kCowboySailor = 1,
    kNegativeTOF = 1 << 1,
    kPositiveTOF = 1 << 2,
    kOptimalParams = 1 << 3
  };

  float fV0radius;                      // V0 decay vertex radius (negarive -> LikeSign V0)
  float fV0pt;                          // V0 transverse momentum (in MC if negative -> fake V0)
  float fV0eta;                         // V0 pseudorapidity
  float fInvMass;                    // Invariant mass for the candidate
  float fDistOverTotMom;                // L/p
  Double32_t fV0CosPA;                  //[0.9,1.0,16] V0 cosine of pointing angle
  Double32_t fChi2V0;                   //[0.0,10.24,8] V0 fit chi2
  Double32_t fDcaNeg2PrimaryVertex;     //[0.0,0.256,8] DCA of the negative prong to the PV
  Double32_t fDcaPos2PrimaryVertex;     //[0.0,0.256,8]  DCA of the positive prong to the PV
  Double32_t fDcaV0daughters;           //[0.0,2.56,8] DCA between the two prongs
  Double32_t fV0armAlpha;               //[-1.28,1.28,8] Armenteros alpha
  Double32_t fV0armPt;                  //[0.0,0.256,8] Armenteros pt
  unsigned char fLeastNxedRows;         // Min number of xed roads
  Double32_t fLeastXedOverFindable;     //[0.0,1.0,8] Min xed roads/findable clusters
  Double32_t fMaxChi2PerCluster;        //[0,12.8,8] Max chi2 per cluster in TPC
  Double32_t fNsigmaPionPos;            //[0.0,8.0,4] # sigma TPC pion for the positive prong
  Double32_t fNsigmaHe3Pos;  
  Double32_t fNsigmaPionNeg;            //[0.0,8.0,4] # sigma TPC pion for the positive prong
  Double32_t fNsigmaHe3Neg;          //[0.0,8.0,4] # sigma TPC proton for the positive prong
  Double32_t fPtPos;
  Double32_t fPtNeg;
  Double32_t fPhiPos;
  Double32_t fPhiNeg;
  Double32_t fEtaPos;                   //[-1.0,1.0,7] Pseudorapidity of the positive prong. MSB is the TOF bit.
  Double32_t fEtaNeg;                   //[-1.0,1.0,7] Pseudorapidity of the negative prong. MSB is the TOF bit.
  unsigned char fITSInfo;               // Starting from the MSB: kITSrefit for neg and pos, kSPDany for neg and pos, least number of ITS clusters (last 4 bits)
  unsigned char fFlags;                 // Cowboy&Saylor, TOF bits for neg and pos, optimal tracking parameters
};
/// Cowboy&saylor, 2xTOF bit


inline void HyperTriton2Body::SetProngsPvDCA(float posD, float negD) {
  fDcaPos2PrimaryVertex = posD;
  fDcaNeg2PrimaryVertex = negD;
}

inline void HyperTriton2Body::SetArmenterosVariables(float alpha, float pt) {
  fV0armAlpha = alpha;
  fV0armPt = pt;
}

inline void HyperTriton2Body::SetProngsTPCnsigmas(float pPi, float pHe, float nPi,
                                        float nHe) {
  fNsigmaPionPos = pPi;
  fNsigmaHe3Pos = pHe;
  fNsigmaPionNeg = nPi;
  fNsigmaHe3Neg = nHe;
}

inline void HyperTriton2Body::SetProngsEtaTOF(float posEta, float negEta) {
  fEtaNeg = negEta;
  fEtaPos = posEta;
}

inline void HyperTriton2Body::SetProngsPt(float posPt, float negPt) {
  fPtNeg = negPt;
  fPtPos = posPt;
}

inline void HyperTriton2Body::SetProngsPhi(float posPhi, float negPhi) {
  fPhiNeg = negPhi;
  fPhiPos = posPhi;
}

inline void HyperTriton2Body::SetITSinformation(bool negRefit, bool posRefit, bool negSPD, bool posSPD, int nITS) {
  fITSInfo = 0u;
  if (negRefit) fITSInfo |= 1 << 7;
  if (posRefit) fITSInfo |= 1 << 6;
  if (negSPD) fITSInfo |= 1 << 5;
  if (posSPD) fITSInfo |= 1 << 4;
  fITSInfo |= 0x07 & nITS;
}


inline void HyperTriton2Body::SetTOFbits(bool pTOF, bool nTOF) {
  flipBits(fFlags, kPositiveTOF, pTOF);
  flipBits(fFlags, kNegativeTOF, nTOF);
}

}  // namespace Lifetimes

#endif
