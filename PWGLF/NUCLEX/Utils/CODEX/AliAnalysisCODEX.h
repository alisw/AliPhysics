#ifndef ALIANALYSISCODEX_H
#define ALIANALYSISCODEX_H

#include <Rtypes.h>
#include <cmath>
#include <climits>
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include <assert.h>

using std::string;
using std::vector;
using ROOT::Math::XYZVectorF;

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> FourVector_t;

namespace AliAnalysisCODEX {
  /// Constants
  const float kElectronMass = 5.10998909999999971e-04;
  const float kPionMass = 1.39569997787475586e-01f;
  const float kKaonMass = 4.93676990270614624e-01;
  const float kProtonMass = 9.38271999359130859e-01f;
  const float kDeuteronMass = 1.87561297416687012f;
  const float kH3Mass = 2.80925011634826660f;
  const float kHe3Mass = 2.80923008918762207f;
  const float kHe4Mass = 3.72737908363342285f;
  const float kMasses[8] = {
    kElectronMass, kPionMass, kKaonMass, kProtonMass,
    kDeuteronMass,kH3Mass,kHe3Mass,kHe4Mass
  };
  const string kNames[8] = {
    "electron","pion","kaon","proton",
    "deuteron","triton","helium-3","helium-4"
  };
  const string kShortNames[8] = {
    "el","pi","k","p",
    "d","t","he3","he4"
  };
  const char kPosNeg[2] = {'P','N'};

  const float kC = 2.99792457999999984e8; // m / s
  const float kCtof = 2.99792457999999984e-02; // cm / ps

  const float kTPCsigmaBinWidth = 0.08f;
  const float kTPCChi2binWidth = 0.05f;
  const float kITSChi2binWidth = 0.5f;
  const float kGoldenChi2binWidth = 0.5f;
  const float kDCAbinWidth = 0.0004f;
  const float kZvertConv = 1000.f;

  /// Enums
  enum BitMask {
    kEl  = BIT(0),
    kPi  = BIT(1),
    kKa  = BIT(2),
    kPr  = BIT(3),
    kDe  = BIT(4),
    kH3  = BIT(5),
    kHe3 = BIT(6),
    kHe4 = BIT(7),
    kTPCrefit = BIT(8),
    kIsPrimary = BIT(9),
    kIsSecondaryFromMaterial = BIT(10),
    kIsReconstructed = BIT(11),
    kIsFake = BIT(12),
    kTOFmismatch = BIT(13),
    kIsKink = BIT(14),
    kTRDrefit = BIT(15)
  };

  enum EventMask {
    kMCevent = BIT(0),
    kNegativeB = BIT(1),
    kInelGt0 = BIT(2),
    kTriggerClasses = BIT(3)
  };

  enum ITSbits {
    kSPD0 = BIT(0),
    kSPD1 = BIT(1),
    kSDD0 = BIT(2),
    kSDD1 = BIT(3),
    kSSD0 = BIT(4),
    kSSD1 = BIT(5),
    kAny  = 63,
    kITSrefit = BIT(7)
  };

  /// Helper functions
  int BitNo(BitMask v);

  class Header {
    public:
      float GetCentrality()      const  { return mCentrality + 0.5f; }
      float GetVertexZposition() const  { return mZvert / kZvertConv; }
      void  SetCentrality(float cent)   { mCentrality = char(cent); }
      void  SetVertexZposition(float z) { mZvert = round(z * kZvertConv); }

      short  mZvert;                               /// Binned, 100um resolution in the [-30,30] cm range
      char   mCentrality;                          /// Binned, 1% resolution in [0,100%]
      char   mEventMask;                           /// Bitmask with misc info
      char   mT0mask[10];                          /// TOF T0 mask
      float  mT0event[10];                         /// Event T0
      float  mT0resolution[10];                    /// Event T0 resolution

  };

  /// Structs and classes
  class Track {
    public:
      /// Kinematics
      float P()                        const { return fabs(pT) * cosh(eta); }
      float Pt()                       const { return fabs(pT); }
      float Px()                       const { return fabs(pT) * cos(phi); }
      float Py()                       const { return fabs(pT) * sin(phi); }
      float Pz()                       const { return fabs(pT) * sinh(eta); }
      float Charge()                   const { return pT > 0.f ? 1.f : -1.f; }
      float Eta()                      const { return eta; }
      float Phi()                      const { return phi; }

      float E(float m)                 const { return sqrt(m * m + P() * P()); }
      float Mt(float m)                const { return sqrt(m * m + Pt() * Pt()); }
      float Y(float m)                 const { return asinh(Pt() / Mt(m) * sinh(Eta())); }

      /// Other getters
      float GetTPCmomentum()           const { return pTPC; }
      float GetTPCsignal()             const { return TPCsignal; }
      int   GetTPCNcls()               const { return TPCnClusters; }
      int   GetTPCsignalN()            const { return TPCnSignal; }

      /// TOF business
      float GetIntegratedLength()      const { return length; }
      int   GetTOFchannel()            const { return wildcard; }
      float GetTOFmomentum()           const { return pTOF; }
      int   GetTOFNcls()                     { return TOFnClusters; }
      float GetTOFsignal()             const { return TOFsignal; }
      bool  HasTOF()                   const { return TOFsignal > -1 && length > 350.f; }
      void  SetIntegratedLength(float len)   { length = len; }
      void  SetTOFchannel(int ch)            { wildcard = ch; }
      void  SetTOFmomentum(float ptof)       { pTOF = ptof; }
      void  SetTOFsignal(float signal)       { TOFsignal = signal; }
      void  SetTOFNcls(int tofcls)           { TOFnClusters = tofcls; }
      float TOFbeta()                  const { return (TOFsignal < 1.e-16) ? -1.f : length / (TOFsignal * kCtof); }
      float TOFgamma()                 const { return (TOFsignal < 1.e-16) ? -1.f : (TOFbeta() > 1. - 1.e-16) ? 1.e16 : 1. / sqrt(1. - TOFbeta() * TOFbeta()); }
      float TOFmass()                  const { return (TOFsignal < 1.e-16) ? -1.f : pTPC / (TOFbeta() * TOFgamma()); }

      bool  Is(BitMask bit)            const { return mask & bit; }
      bool  TPCrefit()                 const { return Is(kTPCrefit); }
      bool  ITSrefit()                 const { return ITSmap & kITSrefit; }
      bool  Primary()                  const { return Is(kIsPrimary); }
      bool  SecondaryFromMaterial()    const { return Is(kIsSecondaryFromMaterial); }
      bool  SecondaryFromDecay()       const { return !Primary() && !SecondaryFromMaterial(); }

      int   ITSnClusters()             const;
      bool  HasPointOnITSLayer(int i)  const { return TESTBIT(ITSmap,i); }

      float GetNumberOfSigmasTPC(int species)   const { return (species < 8 && species >= 0) ? TPCsigmas[species] * kTPCsigmaBinWidth : -999.f; }
      float GetNumberOfSigmasTPC(BitMask bit)   const { return GetNumberOfSigmasTPC(BitNo(bit)); }
      void  SetTPCsigma(BitMask bit, float s);

      float GetTPCChi2NDF() const { return TPCchi2NDF * kTPCChi2binWidth;}
      void  SetTPCChi2NDF(float c) { TPCchi2NDF = (c >= kTPCChi2binWidth * 250.f) ? 255u : round(c / kTPCChi2binWidth); }

      float GetITSChi2NDF() const { return ITSchi2NDF * kITSChi2binWidth;}
      void  SetITSChi2NDF(float c) { ITSchi2NDF = (c >= kITSChi2binWidth * 250.f) ? 255u : round(c / kITSChi2binWidth); }

      float GetGoldenChi2NDF() const { return GoldenChi2 * kGoldenChi2binWidth;}
      void  SetGoldenChi2NDF(float c) { GoldenChi2 = (c >= kGoldenChi2binWidth * 250.f) ? 255u : round(c / kGoldenChi2binWidth); }

      float GetDCAxy() const { return DCAxy * kDCAbinWidth;}
      void  SetDCAxy(float c) { DCAxy = round(c / kDCAbinWidth); }
      float GetDCAz() const { return DCAz * kDCAbinWidth;}
      void  SetDCAz(float c) { DCAz = round(c / kDCAbinWidth); }

      unsigned char GetTRDnTracklets() const { return TRDnTracklets; }
      void SetTRDnTracklets(unsigned char nTRDtrklt) { TRDnTracklets = nTRDtrklt; }

      /// Templates
      template<typename F>void  P(F p[3]) const { p[0] = Px(); p[1] = Py(); p[2] = Pz(); }

      float            eta;          /// Eta of the track at the primary vertex
      float            phi;          /// Phi of the track at the primary vertex
      float            pT;           /// Transverse momentum
      float            pTPC;         /// Momentum at the innerwall of the TPC
      float            pTOF;         /// p used to compute the TOF expected time
      float            TPCsignal;    /// TPC dE / dx (a.u.)
      float            TOFsignal;    /// TOF time (T0 already subtracted)
      int              wildcard;     /// In the MC: index of the mother, in the data: TOF channel
      float            length;       /// Track length
      char             TPCsigmas[8]; /// TPC sigmas. If TPC PID is not available ITS sigms is stored.
      unsigned short   mask;         /// Mask (see ne enumerator above for the meaning)
      short            DCAxy;        /// DCAxy (binned)
      short            DCAz;         /// DCAz (binned)
      unsigned char    TPCnClusters; /// Number of clusters in TPC
      unsigned char    TPCnClustersF;/// Number of findable clusters in TPC
      unsigned char    TPCnXedRows;  /// Number of crossed rows in TPC
      unsigned char    TPCnSignal;   /// Number of clusters with PID info in TPC
      unsigned char    TOFnClusters; /// Number of clusters in TOF
      unsigned char    ITSmap;       /// Map of ITS clusters
      unsigned short   ITSSignal[4]; /// Energy loss (binned) in the sensitive layers of the ITS
      unsigned char    TPCchi2NDF;   /// Chi2/ndf (binned) in TPC
      unsigned char    ITSchi2NDF;   /// Chi2/ndf (binned) in ITS
      unsigned char    GoldenChi2;   /// Golden Chi2 defined as Constrained Global TPC chi2 (binned)
      unsigned char    ActiveLength; /// Active length in the TPCchi2NDF
      unsigned char    TRDnTracklets;/// Number of TRD tracklets associated with the track

  };

  class TOFpidLite {
    public:
      TOFpidLite(Header* head = 0x0);
      void  ConnectHeader(Header* head) { mHeader = head; }
      float GetExpectedSignal(const Track &t, float mass) const;
      float GetExpectedSigma(const Track &t, float mass)  const;
      float GetNumberOfSigmas(const Track &t, int species) const { return (species >=0 && species < 8) ? GetNumberOfSigmas(t, kMasses[species]) : -999.f; }
      float GetNumberOfSigmas(const Track &t, float mass)  const;
      void  SetParams(const float params[4]);
      int   GetTOFchannel(const Track& t) const;
      void  SetTOFresolution(const float TOFresolution) { mTOFres = TOFresolution; }
      float GetTOFresolution() const { return mTOFres; }
      void  SetTOFtail(const float TOFtail) { mTOFtail = TOFtail; };
      float GetTOFtail() const { return mTOFtail; }
      float GetT0event(const Track &t) const { return mHeader->mT0event[GetMomBin(t.GetTPCmomentum())]; }
      char  GetT0source(const Track &t) const { return mHeader->mT0mask[GetMomBin(t.GetTPCmomentum())]; }

  // fTOFpid
    private:
  // fTOFpid
      int   GetMomBin(float mom) const;
  // fTOFpid

      float   mMomBins[11];
      float   mParams[4];
      float   mTOFres;
      float   mTOFtail;
      Header* mHeader;
  };

  // Template that stores tracks from different events used for background estimation
  // with Mixed-Event technique, essentially tracks are stored in a 4D matrix
  // depending on the value of centrality and vertex of the event
  template <int centr, int vert, int part> class EventMixingPool : public TObject {

    public:
      EventMixingPool(int maxcent = 100, float maxvtz = 10., int depth=10) :
        mPool(),
        mCentralityBins(centr),
        mVertexBins(vert),
        mNparticles(part),
        mDepth(depth),
        mLevel(),
        mCWBin(maxcent / centr),
        mVWBin(2.0f * maxvtz / (float)vert),
        mCMax(maxcent),
        mVMax(maxvtz),
        mPartMass() {
          ResetLevels();
        }

      // fills the vector of tracks in the correct matrix element
      // depending on the value of centrality and vertex the event
      void FillEvent(vector<FourVector_t> &track, char c, float v, int p) {
        if (track.size() > 0 && fabs(v) < mVMax && c < mCMax) {
          int cbin = c / mCWBin;
          int vbin = ( v + mVMax ) / mVWBin;
          int index = mLevel[cbin][vbin][p] < mDepth ? mLevel[cbin][vbin][p] : 0;
          vector<FourVector_t> vv3;
          vv3.resize(track.size());
          for (size_t iT = 0; iT < track.size(); iT++) {
            FourVector_t ntr = track[iT];
            FourVector_t vec = {ntr.Pt(), ntr.Eta(), ntr.Phi(), ntr.M()};
            vv3[iT] = vec;
          }
          mPool[cbin][vbin][p][index] = vv3;
          mLevel[cbin][vbin][p] = index + 1;
        }
      }

      // sets to 0 every occupation level and clear all mPool elements
      void ResetLevels() {
        for (int i = 0; i < mCentralityBins; i++) {
          for (int j = 0; j < mVertexBins; j++) {
            for (int k = 0; k < mNparticles; k++) {
              mLevel[i][j][k] = 0;
              mPool[i][j][k].resize(mDepth);
            }
          }
        }
      }

      // returns the vector of FourVector_t stored in the mPool
      vector<FourVector_t> GetVectorFV(int c, float v, int p, int d) {
        if (fabs(v) < mVMax && c < mCMax) {
          int cbin = c / mCWBin;
          int vbin = ( v + mVMax ) / mVWBin;
          return mPool[cbin][vbin][p][d];
        }
      }

      // returns centrality bin of mPool for a given value of centrality
      int GetCentralityBin(const int c) const              { return int(c / mCWBin); }
      // returns vertex bin of mPool for a given value of vertex
      int GetVertexBin(const float v) const                { return int(( v + mVMax ) / mVWBin); }

      // returns number of level occupied in mPool[c][v][p]
      int   GetNLevelOccupied(int c, float v, int p) {
        if (fabs(v) < mVMax && c < mCMax) {
          int cbin = c / mCWBin;
          int vbin = ( v + mVMax ) / mVWBin;
          return mLevel[cbin][vbin][p];
        } else {
          return 0;
        }
      }

      void  SetmCMax(const float cmax)                      { mCMax = cmax; }
      float GetmCMax() const                                { return mCMax; }

      void  SetmVMax(const float vmax)                      { mVMax = vmax; }
      float GetmVMax() const                                { return mVMax; }

      void  SetPartMassI(const int parti, const float mass) { mPartMass[parti] = mass; }
      float GetPartMassI(const int parti) const             { return mPartMass[parti]; }

    private:

      vector<vector<FourVector_t> > mPool[centr][vert][part];
      const int mCentralityBins;
      const int mVertexBins;
      const int mNparticles;
      const int mDepth;

      int       mLevel[centr][vert][part];

      int       mCWBin;
      int       mCMax;
      float     mVWBin;
      float     mVMax;
      float     mPartMass[part];
  };
}
#endif
