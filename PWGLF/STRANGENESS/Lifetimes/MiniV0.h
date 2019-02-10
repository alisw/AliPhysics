#ifndef Lifetimes_MiniV0_h
#define Lifetimes_MiniV0_h

#include <cmath>
#include <cstdlib>
#include <limits>
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "Rtypes.h"
#include "Utils.h"



namespace Lifetimes {

constexpr double Squ(double x) { return x * x; }
class MiniV0 {
 public:
  float GetV0radius() const { return std::abs(fV0radius); }
  float GetV0pt() const { return std::abs(fV0pt); }
  float GetV0eta() const { return fV0eta; }
  float GetKInvMass() const { return fKInvMass; }
  float GetLambdaInvMass() const { return fLambdaInvMass; }  
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
  Double32_t GetNegProngTPCnsigmaProton() const { return fNsigmaProtonNeg;}
  Double32_t GetPosProngTPCnsigmaPion() const {return fNsigmaPionPos;}
  Double32_t GetPosProngTPCnsigmaProton() const {return fNsigmaProtonPos;}
  Double32_t GetNegProngEta() const { return fEtaNeg;}
  Double32_t GetPosProngEta() const { return fEtaPos;}
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
  int  GetLeastNumberOfITSclusters() const { return fITSInfo & 0x07; }

  void SetV0radiusAndLikeSign(float r, bool ls = false) { fV0radius = ls ? -r : r; }
  void SetV0ptAndFake(float pt, bool fake) { fV0pt = fake ? -pt : pt; }
  void SetV0eta(float eta) { fV0eta = eta; }
  void SetKInvMass(float m) { fKInvMass = m; }
  void SetLambdaInvMass(float m2){fLambdaInvMass= m2;}
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
  void SetProngsEta(float posEta, float negEta);
  void SetITSinformation(bool, bool, bool, bool, int);
  void SetCowboyAndSailor(bool cs) { fFlags = flipBits(fFlags, static_cast<unsigned char>(kCowboySailor), cs); }
  void SetTOFbits(bool pTOF, bool nTOF);
  void SetOptimalParameters(bool opt) { fFlags = flipBits(fFlags, static_cast<unsigned char>(kOptimalParams), opt); }
  static MiniV0 FillMiniV0(AliESDv0 *v0, AliESDtrack *pTrack , AliESDtrack *nTrack, float nsigmaposproton
,float nsigmanegproton, float nsigmapospion,float nsigmanegpion, float magneticField , double primaryVertex[3],bool fake);

 private:
  enum Flags {
    kCowboySailor = 1,
    kNegativeTOF = 1 << 1,
    kPositiveTOF = 1 << 2,
    kOptimalParams = 1 << 3
  };
  float fV0radius;                      // V0 decay vertex radius (negarive -> LikeSign V0)
  float fV0pt;                          // V0 transverse momentum (in MC if negative -> fake V0)
  float fV0eta;                         // V0 pseudorapidity
  float fKInvMass;                      // k invariant mass for the candidate
  float fLambdaInvMass;                 // lambda invariant mass for the candidate
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
  Double32_t fNsigmaProtonPos;          //[0.0,8.0,4] # sigma TPC proton for the positive prong
  Double32_t fNsigmaPionNeg;            //[0.0,8.0,4] # sigma TPC pion for the negative prong
  Double32_t fNsigmaProtonNeg;          //[0.0,8.0,4] # sigma TPC proton for the negative prong
  Double32_t fEtaPos;                   //[-1.0,1.0,7] Pseudorapidity of the positive prong. MSB is the TOF bit.
  Double32_t fEtaNeg;                   //[-1.0,1.0,7] Pseudorapidity of the negative prong. MSB is the TOF bit.
  unsigned char fITSInfo;               // Starting from the MSB: kITSrefit for neg and pos, kSPDany for neg and pos, least number of ITS clusters (last 4 bits)
  unsigned char fFlags;                 // Cowboy&Saylor, TOF bits for neg and pos, optimal tracking parameters
};


inline void MiniV0::SetProngsPvDCA(float posD, float negD) {
  fDcaPos2PrimaryVertex = posD;
  fDcaNeg2PrimaryVertex = negD;
}

inline void MiniV0::SetArmenterosVariables(float alpha, float pt) {
  fV0armAlpha = alpha;
  fV0armPt = pt;
}

inline void MiniV0::SetProngsTPCnsigmas(float pPi, float pP, float nPi,
                                        float nP) {
  fNsigmaPionNeg = pPi;
  fNsigmaProtonNeg = pP;
  fNsigmaPionPos = pPi;
  fNsigmaProtonPos = pP;
}

inline void MiniV0::SetProngsEta(float posEta, float negEta) {
  fEtaNeg = negEta;
  fEtaPos = posEta;
}

inline void MiniV0::SetITSinformation(bool negRefit, bool posRefit, bool negSPD, bool posSPD, int nITS) {
  fITSInfo = 0u;
  if (negRefit) fITSInfo |= 1 << 7;
  if (posRefit) fITSInfo |= 1 << 6;
  if (negSPD) fITSInfo |= 1 << 5;
  if (posSPD) fITSInfo |= 1 << 4;
  fITSInfo |= 0x07 & nITS;
}

inline void MiniV0::SetTOFbits(bool pTOF, bool nTOF) {
  flipBits(fFlags, static_cast<unsigned char>(kPositiveTOF), pTOF);
  flipBits(fFlags, static_cast<unsigned char>(kNegativeTOF), nTOF);
}

inline MiniV0 MiniV0::FillMiniV0(AliESDv0 *v0, AliESDtrack *pTrack , AliESDtrack *nTrack, float nsigmaposproton
,float nsigmanegproton, float nsigmapospion,float nsigmanegpion, float magneticField , double primaryVertex[3], bool fake){

MiniV0 miniV0;
double decayVtx[3];
v0->GetXYZ(decayVtx[0], decayVtx[1], decayVtx[2]);
double v0Radius = std::hypot(decayVtx[0], decayVtx[1]);
double masses[2];
v0->ChangeMassHypothesis(310);
masses[0] = v0->GetEffMass();
v0->ChangeMassHypothesis(3122);
masses[1] = v0->GetEffMass();

double v0Pt = v0->Pt();
double tV0mom[3];
v0->GetPxPyPz(tV0mom[0], tV0mom[1], tV0mom[2]);
double lV0TotalMomentum = std::sqrt(tV0mom[0] * tV0mom[0] + tV0mom[1] * tV0mom[1] + tV0mom[2] * tV0mom[2]);
float distOverP = std::sqrt(Squ(decayVtx[0] - primaryVertex[0]) +
                            Squ(decayVtx[1] - primaryVertex[1]) +
                            Squ(decayVtx[2] - primaryVertex[2])) /
                  (lV0TotalMomentum + 1e-16); 
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


miniV0.SetV0ptAndFake(v0Pt, fake);
miniV0.SetV0eta(v0->Eta());
miniV0.SetLeastNumberOfXedRows(minXedRows);
miniV0.SetDistOverP(distOverP);
miniV0.SetKInvMass(masses[0]);
miniV0.SetLambdaInvMass(masses[1]);
miniV0.SetArmenterosVariables(v0->AlphaV0(), v0->PtArmV0());
miniV0.SetV0CosPA(cosPA);
miniV0.SetV0Chi2(v0->GetChi2V0());
miniV0.SetProngsDCA(v0->GetDcaV0Daughters());
miniV0.SetProngsPvDCA(dcaPosToPrimVertex, dcaNegToPrimVertex);
miniV0.SetV0radiusAndLikeSign(v0Radius);
miniV0.SetLeastXedRowsOverFindable(minXedRowsOverFindable);
miniV0.SetMaxChi2perCluster(maxChi2PerCluster);
miniV0.SetProngsEta(pTrack->Eta(), nTrack->Eta());
miniV0.SetProngsTPCnsigmas(nsigmapospion, nsigmaposproton,
                                 nsigmanegpion, nsigmanegproton);
miniV0.SetITSinformation(negITSrefit, posITSrefit, negSPDany, posSPDany, ITSnCl);
miniV0.SetTOFbits(posTOF, negTOF);
miniV0.SetCowboyAndSailor(isCowboy);
return miniV0;

}




}  // namespace Lifetimes

#endif
