#include "AliAnalysisCODEX.h"
#include <vector>

namespace AliAnalysisCODEX {

  int BitNo(BitMask v) {
    // From https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
    static const unsigned int b[] = {
      0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0,0xFF00FF00, 0xFFFF0000
    };
    unsigned int r = (v & b[0]) != 0;
    r |= ((v & b[4]) != 0) << 4;
    r |= ((v & b[3]) != 0) << 3;
    r |= ((v & b[2]) != 0) << 2;
    r |= ((v & b[1]) != 0) << 1;
    return r;
  }

  void Track::SetTPCsigma(BitMask bit, float s) {
    const int iB = BitNo(bit);
    if (iB < 8) {
      if (fabs(s) <= 10.f)
        TPCsigmas[iB] = char(round(s / kTPCsigmaBinWidth));
      else
        TPCsigmas[iB] = s > 0 ? 127 : -127;
    }
  }

  int Track::ITSnClusters() const {
    return
      bool(ITSmap & 1) + bool(ITSmap & 2) + \
      bool(ITSmap & 4) + bool(ITSmap & 8) + \
      bool(ITSmap & 16) + bool(ITSmap & 32);
  }

  TOFpidLite::TOFpidLite(Header* head) :
    mMomBins{0.3f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.f, 1.2f, 1.5f, 2.f, 3.f},
    //mParams{0.008, 0.008, 0.002, 40},
    mParams{0.004f, 0.f, 0.009f, 15.f},
    mTOFres{80.f},
    mTOFtail{0.75f},
    mHeader{head} {}

  float TOFpidLite::GetExpectedSignal(const Track& track, float mass) const {
    const float& pTOF = track.GetTOFmomentum();
    const float energy = sqrt((mass * mass) + (pTOF * pTOF));
    return track.GetIntegratedLength() * energy / (kCtof * pTOF);
  }

  float TOFpidLite::GetExpectedSigma(const Track& track, float mass) const {
    const float& p = track.P();
    const float time = GetExpectedSignal(track, mass);
    const float dpp = mParams[0] + (mParams[1] * p) + (mParams[2] * mass / p);      //mean relative pt resolution;
    const float sigma = dpp * time / (1. + (p * p / (mass * mass)));
    const float t0res = mHeader->mT0resolution[GetMomBin(p)];

    return sqrt((sigma * sigma) + (mParams[3] * mParams[3] / (p * p)) + (mTOFres * mTOFres) + (t0res * t0res));
  }

  void TOFpidLite::SetParams(const float params[4]) {
    for (int i = 0; i < 4; i++) {
      mParams[i] = params[i];
    }
  }

  float TOFpidLite::GetNumberOfSigmas(const Track &t, float mass)  const {
    const float sigma = GetExpectedSigma(t, mass);
    const float mu = GetExpectedSignal(t, mass);
    return (t.GetTOFsignal() - mu) / sigma;
  }

  int TOFpidLite::GetMomBin(float p) const {
    for (int i = 0; i < 10; ++i)
      if (p < mMomBins[i + 1] && p > mMomBins[i]) return i;
    return 10;
  }

  int TOFpidLite::GetTOFchannel(const Track &t) const {
    const float etaAbs = fabs(t.Eta());
    int channel = int(4334.09 - 4758.36 * etaAbs -1989.71 * etaAbs*etaAbs + 1957.62*etaAbs*etaAbs*etaAbs);
    if (channel < 1 || etaAbs > 1.f) channel = 1;
    return channel;
  }
}
