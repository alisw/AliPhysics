///
/// \file AliFemtoTrack.h
///

#pragma once

#ifndef ALIFEMTOTRACK_H
#define ALIFEMTOTRACK_H

#include "AliFemtoTypes.h"
#include "AliFmPhysicalHelixD.h"
#include "TBits.h"
#include "AliFemtoHiddenInfo.h"

#include <algorithm>


/// \class AliFemtoTrack
/// \brief Main class holding track information (before identification)
///
/// AliFemtoTrack holds all the necessary information about a track that is
/// required during femtoscopic analysis. This class is filled with information
/// from the input stream by the reader. A particle has a link back to the Track
/// it was created from, so we do not copy the information.
///
class AliFemtoTrack {
public:
  AliFemtoTrack();
  AliFemtoTrack(const AliFemtoTrack& aTrack);
  ~AliFemtoTrack();
  AliFemtoTrack& operator=(const AliFemtoTrack& aTrack);

  short Charge() const;
  float PidProbElectron() const;
  float PidProbPion() const;
  float PidProbKaon() const;
  float PidProbProton() const;
  float PidProbMuon() const;

  const AliFemtoThreeVector& P() const;
  float Pt() const;
  float InnerMomentum() const;

  const AliFmPhysicalHelixD& Helix() const;
  int TrackId() const;
  long int Flags() const;
  int Label() const;
  float ImpactD() const;

  float ImpactDprim() const;
  float ImpactDweak() const;
  float ImpactDmat() const;

  float ImpactZ() const;
  float Cdd() const;
  float Cdz() const;
  float Czz() const;

  float ITSchi2perNDF() const;
  float ITSchi2() const;
  int   ITSncls() const;
  float TPCchi2() const;

  /// Calculate reduced chi-squared (χ²/NDoF) for TPC
  ///
  /// Calculation comes from AliAODTrack::GetTPCchi2perNDF
  ///
  float TPCchi2perNDF() const;
  int   TPCncls() const;
  short TPCnclsF() const;
  float TPCsignal() const;
  short TPCsignalN() const;
  float TPCsignalS() const;

  //new PID
  float NSigmaTPCPi() const;
  float NSigmaTPCK() const;
  float NSigmaTPCP() const;
  float NSigmaTPCE() const;

  /***********************/
  //
  float NSigmaTPCD() const;
  float NSigmaTPCT() const;
  float NSigmaTPCH() const;
  float NSigmaTPCA() const;
  //
  /***********************/

  float VTOF() const;
  float TOFsignal() const;
  float NSigmaTOFPi() const;
  float NSigmaTOFK() const;
  float NSigmaTOFP() const;
  float NSigmaTOFE() const;

  /**************************/
  //
  float NSigmaTOFD() const;
  float NSigmaTOFT() const;
  float NSigmaTOFH() const;
  float NSigmaTOFA() const;

  float MassTOF() const;
  //
  /**************************/

  float TOFpionTime() const;
  float TOFkaonTime() const;
  float TOFprotonTime() const;

  /***************************/
  //
  float TOFdeuteronTime() const;
  float TOFtritonTime() const;
  float TOFhe3Time() const;
  float TOFalphaTime() const;
  //
  /******************************/

  double XatDCA() const;
  double YatDCA() const;
  double ZatDCA() const;

  float CorrectionPion() const;
  float CorrectionKaon() const;
  float CorrectionProton() const;
  float CorrectionPionMinus() const;
  float CorrectionKaonMinus() const;
  float CorrectionProtonMinus() const;
  float CorrectionAll() const;

  /**************************************/
  //
  float CorrectionDeuteron() const;
  float CorrectionTriton() const;
  float CorrectionHe3() const;
  float CorrectionAlpha() const;
  float CorrectionDeuteronMinus() const;
  float CorrectionTritonMinus() const;
  float CorrectionHe3Minus() const;
  float CorrectionAlphaMinus() const;
  //
  /**************************************/

  const TBits& TPCclusters() const;
  const TBits& TPCsharing()  const;

  void SetCharge(const short& s);
  void SetPidProbElectron(const float& x);
  void SetPidProbPion(const float& x);
  void SetPidProbKaon(const float& x);
  void SetPidProbProton(const float& x);
  void SetPidProbMuon(const float& x);
  void SetTofExpectedTimes(const float& tpi, const float& tkn, const float& tpr, const float& ttof);


  void SetP(const AliFemtoThreeVector& p);
  void SetPt(const float& x);
  void SetInnerMomentum(const float& x);
  void SetHelix(const AliFmPhysicalHelixD& h);
  void SetTrackId(const int& s);
  void SetFlags(const long int& i);
  void SetLabel(const int& i);
  void SetImpactD(const float& x);

  void SetImpactDprim(const float& x);
  void SetImpactDweak(const float& x);
  void SetImpactDmat(const float& x);

  void SetImpactZ(const float& x);
  void SetCdd(const float& x);
  void SetCdz(const float& x);
  void SetCzz(const float& x);

  void SetITSchi2(const float& x);
  void SetITSncls(const int& i);
  void SetTPCchi2(const float& x);
  void SetTPCncls(const int& i);
  void SetTPCnclsF(const short& s);
  void SetTPCsignal(const float& s);
  void SetTPCsignalN(const short& s);
  void SetTPCsignalS(const float& x);

  //new PID
  void SetNSigmaTPCPi(const float& x);
  void SetNSigmaTPCK(const float& x);
  void SetNSigmaTPCP(const float& x);
  void SetNSigmaTPCE(const float& x);
  void SetVTOF(const float& x);
  void SetTOFsignal(const float& x);
  void SetNSigmaTOFPi(const float& x);
  void SetNSigmaTOFK(const float& x);
  void SetNSigmaTOFP(const float& x);
  void SetNSigmaTOFE(const float& x);

  /***************************************************/
  //
  void SetNSigmaTPCD(const float& x);
  void SetNSigmaTPCT(const float& x);
  void SetNSigmaTPCH(const float& x);
  void SetNSigmaTPCA(const float& x);
  void SetNSigmaTOFD(const float& x);
  void SetNSigmaTOFT(const float& x);
  void SetNSigmaTOFH(const float& x);
  void SetNSigmaTOFA(const float& x);

  void SetMassTOF(const float& x);
  //
  /**************************************************/

  void SetTPCcluster(const short& aNBit, const Bool_t& aValue);
  void SetTPCshared(const short& aNBit, const Bool_t& aValue);

  void SetTPCClusterMap(const TBits& aBits);
  void SetTPCSharedMap(const TBits& aBits);

  void SetKinkIndexes(const int points[3]);
  int  KinkIndex(int aIndex) const;
  void SetITSHitOnLayer(int i, bool val);
  bool HasPointOnITSLayer(int aIndex) const; ///< i: 0-5, for 6 layers

  void SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo);
  bool ValidHiddenInfo() const;
  AliFemtoHiddenInfo* GetHiddenInfo() const;

  const AliFemtoThreeVector& NominalTpcExitPoint() const;
  const AliFemtoThreeVector& NominalTpcPoint(int i) const;
  const AliFemtoThreeVector& NominalTpcEntrancePoint() const;
  const AliFemtoThreeVector& NominalTpcPointShifted() const;

  void SetNominalTPCEntrancePoint(const AliFemtoThreeVector& aXTPC);
  void SetNominalTPCEntrancePoint(double *aXTPC);

  void SetNominalTPCPoints(const AliFemtoThreeVector * const);
  void SetNominalTPCPoints(double **aXTPC);

  void SetNominalTPCExitPoint(const AliFemtoThreeVector& aXTPC);
  void SetNominalTPCExitPoint(double *aXTPC);

  void SetNominalTPCPointShifted(const AliFemtoThreeVector& aXTPC);
  void SetNominalTPCPointShifted(double *aXTPC);

  void SetSigmaToVertex(const float& Sigma);
  float SigmaToVertex() const;

  void SetXatDCA(const double& x);
  void SetYatDCA(const double& x);
  void SetZatDCA(const double& x);

  void SetCorrectionPion(const double& x);
  void SetCorrectionKaon(const double& x);
  void SetCorrectionProton(const double& x);
  void SetCorrectionPionMinus(const double& x);
  void SetCorrectionKaonMinus(const double& x);
  void SetCorrectionProtonMinus(const double& x);
  void SetCorrectionAll(const double& x);

  /**********************************************/
  //
  void SetCorrectionDeuteron(const double& x);
  void SetCorrectionTriton(const double& x);
  void SetCorrectionHe3(const double& x);
  void SetCorrectionAlpha(const double& x);
  void SetCorrectionDeuteronMinus(const double& x);
  void SetCorrectionTritonMinus(const double& x);
  void SetCorrectionHe3Minus(const double& x);
  void SetCorrectionAlphaMinus(const double& x);
  //
  /********************************************/

  void SetTrueMomentum(AliFemtoThreeVector *aMom);
  void SetTrueMomentum(const AliFemtoThreeVector& aMom);
  void SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz);
  void SetEmissionPoint(AliFemtoLorentzVector *aPos);
  void SetEmissionPoint(const AliFemtoLorentzVector& aPos);
  void SetEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT);
  void SetPDGPid(Int_t aPid);
  void SetMass(Double_t aMass);

  AliFemtoThreeVector   *GetTrueMomentum() const;
  const AliFemtoLorentzVector *GetEmissionPoint() const;
  Int_t                  GetPDGPid() const;
  Double_t               GetMass() const;

  const AliFemtoThreeVector   *GetGlobalEmissionPoint() const;
  void                   SetGlobalEmissionPoint(const AliFemtoThreeVector& aPos);
  void                   SetGlobalEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz);


  void SetPrimaryVertex(const double *vertex);
  void GetPrimaryVertex(double *fisvertex) const;

  int Multiplicity() const;
  double Zvtx() const;
  void SetMultiplicity(int mult);
  void SetZvtx(double vtx);

  //Alice stuff
  enum {
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kRICHpid=0x20000,
    kTRDbackup=0x80000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000
  };

 private:
  char  fCharge;          ///< track charge
  float fPidProbElectron; ///< electron pid
  float fPidProbPion;     ///< pion pid
  float fPidProbKaon;     ///< kaon pid
  float fPidProbProton;   ///< proton pid
  float fPidProbMuon;     ///< muon pid
  int   fTrackId;         ///< track unique id
  float fTofPionTime;     ///< TOF time - pion expected time
  float fTofKaonTime;     ///< TOF time - kaon expected time
  float fTofProtonTime;   ///< TOF time - proton expected time

  /**************************************************************/
  //
  float fTofDeuteronTime; ///< TOF time - deuteron expected time
  float fTofTritonTime;   ///< TOF time - triton expected time
  float fTofHe3Time;      ///< TOF time - he3 expected time
  float fTofAlphaTime;    ///< TOF time - alpha expected time
  //
  /**************************************************************/

  AliFemtoThreeVector fP; ///< track momentum
  float fPt;              ///< transverse momenta
  float fInnerMomentum;   ///< *total* momentum at the *inner* wall of the TPC
  int fMultiplicity;
  double fZvtx;

  AliFmPhysicalHelixD fHelix; ///< track helix
  //alice stuff
  //Long_t  fFlags;
  long int  fFlags;       ///< Reconsruction status flags
  int       fLabel;       ///< Track label
  float     fImpactD;     ///< Impact parameter in xy plane

  float fImpactDprim;     ///< Impact parameter in xy plane
  float fImpactDweak;     ///< Impact parameter in xy plane
  float fImpactDmat;      ///< Impact parameter in xy plane

  float fImpactZ;         ///< impacct parameter in z
  float fCdd,fCdz,fCzz;   ///< covariance matrix of the impact parameters

  // ITS related track information
  float fITSchi2;         ///< chi2 in the ITS
  int   fITSncls;         ///< number of clusters assigned in the ITS

  // TPC related track information
  float fTPCchi2;         ///< chi2 in the TPC
  int   fTPCncls;         ///< number of clusters assigned in the TPC
  short fTPCnclsF;        ///< number of findable clusters in the TPC
  float fTPCsignal;       ///< dEdx TPC value
  short fTPCsignalN;      ///< number of points used for dEdx
  float fTPCsignalS;      ///< RMS of dEdx measurement (not used in AOD files)

  float fVTOF;            ///< v=length/TOF
  float fTOFsignal;            ///< TOF signal
  float fNSigmaTPCPi;     ///< nsigma TPC for pion
  float fNSigmaTPCK;      ///< nsigma TPC for K
  float fNSigmaTPCP;      ///< nsigma TPC for P
  float fNSigmaTPCE;      ///< nsigma TPC for electron
  float fNSigmaTOFPi;     ///< nsigma TPC for pion
  float fNSigmaTOFK;      ///< nsigma TPC for K
  float fNSigmaTOFP;      ///< nsigma TPC for P
  float fNSigmaTOFE;      ///< nsigma TPC for electron

  /*******************************************************/
  //
  float fNSigmaTPCD;      ///< nsigma TPC for deuteron
  float fNSigmaTPCT;      ///< nsigma TPC for triton
  float fNSigmaTPCH;      ///< nsigma TPC for he3
  float fNSigmaTPCA;      ///< nsigma TPC for alpha
  float fNSigmaTOFD;      ///< nsigma TPC for deuteron
  float fNSigmaTOFT;      ///< nsigma TPC for triton
  float fNSigmaTOFH;      ///< nsigma TPC for he3
  float fNSigmaTOFA;      ///< nsigma TPC for alpha
  float fMassTOF;         ///<
  //
  /******************************************************/

  float fSigmaToVertex;   ///< Distance from track to vertex in sigmas
  TBits fClusters;        ///< Cluster per padrow map
  TBits fShared;          ///< Sharing per padrow map

  AliFemtoThreeVector fNominalTpcEntrancePoint;  ///< Nominal track entrance point into TPC
  AliFemtoThreeVector fNominalTpcPoints[9];      ///< Nominal track points in TCP
  AliFemtoThreeVector fNominalTpcExitPoint;      ///< Nominal track exit point from TPC
  AliFemtoThreeVector fNominalTpcPointShifted;   ///< Nominal track at given point from TPC

  int   fKinkIndexes[3];    ///< Kink Index list
  bool  fHasPointOnITS[6];  ///< if track has hit on the ITS layer (6 layers: 2 x 3 (SPD, SSD, SDD))

  double fXatDCA;
  double fYatDCA;
  double fZatDCA;

  AliFemtoHiddenInfo* fHiddenInfo;              //!< hidden info containing MC data

  AliFemtoThreeVector   *fTrueMomentum;         ///< True (simulated) momentumfse
  AliFemtoLorentzVector *fEmissionPoint;        ///< Emission point coordinates
  Int_t                  fPDGPid;               ///< True PID of the particle
  Double_t               fMass;                 ///< True particle mass
  AliFemtoThreeVector   *fGlobalEmissionPoint;  ///< Global emission point

  double fVertex[3];

  //Corrections related information
  float fCorrPi;  ///< corrections for pion hypothesis
  float fCorrK;   ///< corrections for kaon hypothesis
  float fCorrP;   ///< corrections for proton hypothesis

  float fCorrPiMinus;  ///< corrections for pion hypothesis
  float fCorrKMinus;   ///< corrections for kaon hypothesis
  float fCorrPMinus;   ///< corrections for proton hypothesis

  float fCorrAll;  ///< corrections for particles without PID

  /***********************************************************/
  //
  float fCorrD;  ///< corrections for deuteron hypothesis
  float fCorrT;  ///< corrections for triton hypothesis
  float fCorrH;  ///< corrections for he3 hypothesis
  float fCorrA;  ///< corrections for alpha hypothesis

  float fCorrDMinus;  ///< corrections for deuteron hypothesis
  float fCorrTMinus;  ///< corrections for triton hypothesis
  float fCorrHMinus;  ///< corrections frr he3 hypothesis
  float fCorrAMinus;  ///< corrections for alpha hypothesis
  //
  /**************************************************************/

};

//inline const float* AliFemtoTrack::NSigma() const
//{ return &mNSigmaElectron; } // Fab private
inline float AliFemtoTrack::PidProbElectron() const { return fPidProbElectron; }
inline float AliFemtoTrack::PidProbPion() const { return fPidProbPion; }
inline float AliFemtoTrack::PidProbKaon() const { return fPidProbKaon; }
inline float AliFemtoTrack::PidProbProton() const { return fPidProbProton; }
inline float AliFemtoTrack::PidProbMuon() const { return fPidProbMuon; }

inline void AliFemtoTrack::SetPidProbElectron(const float& x) { fPidProbElectron = x; }
inline void AliFemtoTrack::SetPidProbPion(const float& x) { fPidProbPion = x; }
inline void AliFemtoTrack::SetPidProbKaon(const float& x) { fPidProbKaon = x; }
inline void AliFemtoTrack::SetPidProbProton(const float& x) { fPidProbProton = x; }
inline void AliFemtoTrack::SetPidProbMuon(const float& x) { fPidProbMuon = x; }

inline void AliFemtoTrack::SetMultiplicity(int mult) { fMultiplicity = mult; }
inline int AliFemtoTrack::Multiplicity() const { return fMultiplicity; }

inline void AliFemtoTrack::SetZvtx(double vtx) { fZvtx=vtx; }
inline double AliFemtoTrack::Zvtx() const { return fZvtx; }

inline float AliFemtoTrack::TPCchi2perNDF() const
{
  Int_t ndof = 2 * fTPCncls - 5;
  return __builtin_expect(ndof > 0, 1) ? fTPCchi2 / ndof : 9999;
}

inline float AliFemtoTrack::ITSchi2perNDF() const
{
  Int_t ndof = 2 * fITSncls - 5;
  return __builtin_expect(ndof > 0, 1) ? fITSchi2 / ndof : 9999;
}

inline void AliFemtoTrack::SetCharge(const short& ch) { fCharge = ch; }
inline short AliFemtoTrack::Charge() const { return fCharge; }

inline void AliFemtoTrack::SetP(const AliFemtoThreeVector& p) { fP = p; }
inline const AliFemtoThreeVector& AliFemtoTrack::P() const { return fP; }

inline void AliFemtoTrack::SetPt(const float& pt) { fPt = pt; }
inline float AliFemtoTrack::Pt() const { return fPt; }

inline void AliFemtoTrack::SetInnerMomentum(const float& x) { fInnerMomentum = x; }
inline float AliFemtoTrack::InnerMomentum() const { return fInnerMomentum; }

inline void AliFemtoTrack::SetHelix(const AliFmPhysicalHelixD& h) { fHelix = h; }
inline const AliFmPhysicalHelixD& AliFemtoTrack::Helix() const { return fHelix; }

inline void AliFemtoTrack::SetTrackId(const int & id) { fTrackId = id;}
inline int AliFemtoTrack::TrackId() const { return fTrackId; }

inline void AliFemtoTrack::SetFlags(const long int &flags) { fFlags = flags; }
inline long int AliFemtoTrack::Flags() const { return fFlags; }

inline void AliFemtoTrack::SetLabel(const int &label) { fLabel = label; }
inline int AliFemtoTrack::Label() const { return fLabel; }

inline void AliFemtoTrack::SetImpactD(const float& aImpactD) { fImpactD = aImpactD; }
inline float AliFemtoTrack::ImpactD() const { return fImpactD; }

inline void AliFemtoTrack::SetImpactDprim(const float& aImpactDprim) { fImpactDprim = aImpactDprim; }
inline float AliFemtoTrack::ImpactDprim() const { return fImpactDprim; }

inline void AliFemtoTrack::SetImpactDweak(const float& aImpactDweak) { fImpactDweak = aImpactDweak; }
inline float AliFemtoTrack::ImpactDweak() const { return fImpactDweak; }

inline void AliFemtoTrack::SetImpactDmat(const float& aImpactDmat) { fImpactDmat = aImpactDmat; }
inline float AliFemtoTrack::ImpactDmat() const { return fImpactDmat; }

inline void AliFemtoTrack::SetImpactZ(const float& aImpactZ) { fImpactZ = aImpactZ; }
inline void AliFemtoTrack::SetCdd(const float& aCdd) { fCdd = aCdd; }
inline void AliFemtoTrack::SetCdz(const float& aCdz) { fCdz = aCdz; }
inline void AliFemtoTrack::SetCzz(const float& aCzz) { fCzz = aCzz; }
inline void AliFemtoTrack::SetITSchi2(const float& aITSchi2) { fITSchi2 = aITSchi2; }
inline void AliFemtoTrack::SetITSncls(const int& aITSncls) { fITSncls = aITSncls; }
inline void AliFemtoTrack::SetTPCchi2(const float& aTPCchi2) { fTPCchi2 = aTPCchi2; }
inline void AliFemtoTrack::SetTPCncls(const int& aTPCncls) { fTPCncls = aTPCncls; }
inline void AliFemtoTrack::SetTPCnclsF(const short& aTPCnclsF) { fTPCnclsF = aTPCnclsF; }
inline void AliFemtoTrack::SetTPCsignal(const float& aTPCsig) { fTPCsignal = aTPCsig; }
inline void AliFemtoTrack::SetTPCsignalN(const short& aTPCsignalN) { fTPCsignalN = aTPCsignalN; }
inline void AliFemtoTrack::SetTPCsignalS(const float& aTPCsignalS) { fTPCsignalS = aTPCsignalS; }
inline void AliFemtoTrack::SetVTOF(const float& aVTOF) { fVTOF = aVTOF; }
inline void AliFemtoTrack::SetTOFsignal(const float& aTOFsignal) { fTOFsignal = aTOFsignal; }
inline void AliFemtoTrack::SetNSigmaTPCPi(const float& aNSigmaTPCPi) { fNSigmaTPCPi = aNSigmaTPCPi; }
inline void AliFemtoTrack::SetNSigmaTPCK(const float& aNSigmaTPCK) { fNSigmaTPCK = aNSigmaTPCK; }
inline void AliFemtoTrack::SetNSigmaTPCP(const float& aNSigmaTPCP) { fNSigmaTPCP = aNSigmaTPCP; }
inline void AliFemtoTrack::SetNSigmaTPCE(const float& aNSigmaTPCE) { fNSigmaTPCE = aNSigmaTPCE; }
inline void AliFemtoTrack::SetNSigmaTPCD(const float& aNSigmaTPCD) { fNSigmaTPCD = aNSigmaTPCD; }
inline void AliFemtoTrack::SetNSigmaTPCT(const float& aNSigmaTPCT) { fNSigmaTPCT = aNSigmaTPCT; }
inline void AliFemtoTrack::SetNSigmaTPCH(const float& aNSigmaTPCH) { fNSigmaTPCH = aNSigmaTPCH; }
inline void AliFemtoTrack::SetNSigmaTPCA(const float& aNSigmaTPCA) { fNSigmaTPCA = aNSigmaTPCA; }
inline void AliFemtoTrack::SetMassTOF(const float& aMassTOF) { fMassTOF = aMassTOF; }

inline void AliFemtoTrack::SetNSigmaTOFPi(const float& aNSigmaTOFPi) { fNSigmaTOFPi = aNSigmaTOFPi; }
inline void AliFemtoTrack::SetNSigmaTOFK(const float& aNSigmaTOFK) { fNSigmaTOFK = aNSigmaTOFK; }
inline void AliFemtoTrack::SetNSigmaTOFP(const float& aNSigmaTOFP) { fNSigmaTOFP = aNSigmaTOFP; }
inline void AliFemtoTrack::SetNSigmaTOFE(const float& aNSigmaTOFE) { fNSigmaTOFE = aNSigmaTOFE; }
inline void AliFemtoTrack::SetNSigmaTOFD(const float& aNSigmaTOFD) { fNSigmaTOFD = aNSigmaTOFD; }
inline void AliFemtoTrack::SetNSigmaTOFT(const float& aNSigmaTOFT) { fNSigmaTOFT = aNSigmaTOFT; }
inline void AliFemtoTrack::SetNSigmaTOFH(const float& aNSigmaTOFH) { fNSigmaTOFH = aNSigmaTOFH; }
inline void AliFemtoTrack::SetNSigmaTOFA(const float& aNSigmaTOFA) { fNSigmaTOFA = aNSigmaTOFA; }
inline void AliFemtoTrack::SetSigmaToVertex(const float& aSigma) { fSigmaToVertex = aSigma; }

inline void AliFemtoTrack::SetXatDCA(const double& x) { fXatDCA = x; }
inline void AliFemtoTrack::SetYatDCA(const double& x) { fYatDCA = x; }
inline void AliFemtoTrack::SetZatDCA(const double& x) { fZatDCA = x; }


inline void AliFemtoTrack::SetCorrectionPion(const double& x) { fCorrPi = x; }
inline void AliFemtoTrack::SetCorrectionKaon(const double& x) { fCorrK = x; }
inline void AliFemtoTrack::SetCorrectionProton(const double& x) { fCorrP = x; }
inline void AliFemtoTrack::SetCorrectionDeuteron(const double& x) { fCorrD = x; }
inline void AliFemtoTrack::SetCorrectionTriton(const double& x) { fCorrT = x; }
inline void AliFemtoTrack::SetCorrectionHe3(const double& x) { fCorrH = x; }
inline void AliFemtoTrack::SetCorrectionAlpha(const double& x) { fCorrA = x; }
inline void AliFemtoTrack::SetCorrectionPionMinus(const double& x) { fCorrPiMinus = x; }
inline void AliFemtoTrack::SetCorrectionKaonMinus(const double& x) { fCorrKMinus = x; }
inline void AliFemtoTrack::SetCorrectionProtonMinus(const double& x) { fCorrPMinus = x; }
inline void AliFemtoTrack::SetCorrectionDeuteronMinus(const double& x) { fCorrDMinus = x; }
inline void AliFemtoTrack::SetCorrectionTritonMinus(const double& x) { fCorrTMinus = x; }
inline void AliFemtoTrack::SetCorrectionHe3Minus(const double& x) { fCorrHMinus = x; }
inline void AliFemtoTrack::SetCorrectionAlphaMinus(const double& x) { fCorrAMinus = x; }
inline void AliFemtoTrack::SetCorrectionAll(const double& x) { fCorrAll = x; }

inline float AliFemtoTrack::ImpactZ() const { return fImpactZ; }
inline float AliFemtoTrack::Cdd() const { return fCdd; }
inline float AliFemtoTrack::Cdz() const { return fCdz; }
inline float AliFemtoTrack::Czz() const { return fCzz; }
inline float AliFemtoTrack::ITSchi2() const { return fITSchi2; }
inline int   AliFemtoTrack::ITSncls() const { return fITSncls; }
inline float AliFemtoTrack::TPCchi2() const { return fTPCchi2; }
inline int   AliFemtoTrack::TPCncls() const { return fTPCncls; }
inline short AliFemtoTrack::TPCnclsF() const { return fTPCnclsF; }
inline float AliFemtoTrack::TPCsignal() const { return fTPCsignal; }
inline short AliFemtoTrack::TPCsignalN() const { return fTPCsignalN; }
inline float AliFemtoTrack::TPCsignalS() const { return fTPCsignalS; }
inline float AliFemtoTrack::VTOF() const { return fVTOF; }
inline float AliFemtoTrack::TOFsignal() const { return fTOFsignal; }
inline float AliFemtoTrack::NSigmaTPCPi() const { return fNSigmaTPCPi; }
inline float AliFemtoTrack::NSigmaTPCK() const { return fNSigmaTPCK; }
inline float AliFemtoTrack::NSigmaTPCP() const { return fNSigmaTPCP; }
inline float AliFemtoTrack::NSigmaTPCE() const { return fNSigmaTPCE; }
inline float AliFemtoTrack::NSigmaTPCD() const { return fNSigmaTPCD; }
inline float AliFemtoTrack::NSigmaTPCT() const { return fNSigmaTPCT; }
inline float AliFemtoTrack::NSigmaTPCH() const { return fNSigmaTPCH; }
inline float AliFemtoTrack::NSigmaTPCA() const { return fNSigmaTPCA; }
inline float AliFemtoTrack::NSigmaTOFPi() const { return fNSigmaTOFPi; }
inline float AliFemtoTrack::NSigmaTOFK() const { return fNSigmaTOFK; }
inline float AliFemtoTrack::NSigmaTOFP() const { return fNSigmaTOFP; }
inline float AliFemtoTrack::NSigmaTOFE() const { return fNSigmaTOFE; }
inline float AliFemtoTrack::NSigmaTOFD() const { return fNSigmaTOFD; }
inline float AliFemtoTrack::NSigmaTOFT() const { return fNSigmaTOFT; }
inline float AliFemtoTrack::NSigmaTOFH() const { return fNSigmaTOFH; }
inline float AliFemtoTrack::NSigmaTOFA() const { return fNSigmaTOFA; }
inline float AliFemtoTrack::MassTOF() const { return fMassTOF; }
inline float AliFemtoTrack::SigmaToVertex() const { return fSigmaToVertex; }
inline float AliFemtoTrack::TOFpionTime() const { return fTofPionTime; }
inline float AliFemtoTrack::TOFkaonTime() const { return fTofKaonTime; }
inline float AliFemtoTrack::TOFprotonTime() const { return fTofProtonTime; }
inline float AliFemtoTrack::TOFdeuteronTime() const { return fTofDeuteronTime; }
inline float AliFemtoTrack::TOFtritonTime() const { return fTofTritonTime; }
inline float AliFemtoTrack::TOFhe3Time() const { return fTofHe3Time; }
inline float AliFemtoTrack::TOFalphaTime() const { return fTofAlphaTime; }

//corrections
inline float AliFemtoTrack::CorrectionPion() const { return fCorrPi; }
inline float AliFemtoTrack::CorrectionKaon() const { return fCorrK; }
inline float AliFemtoTrack::CorrectionProton() const { return fCorrP; }
inline float AliFemtoTrack::CorrectionDeuteron() const { return fCorrD; }
inline float AliFemtoTrack::CorrectionTriton() const { return fCorrT; }
inline float AliFemtoTrack::CorrectionHe3() const { return fCorrH; }
inline float AliFemtoTrack::CorrectionAlpha() const { return fCorrA; }
inline float AliFemtoTrack::CorrectionPionMinus() const { return fCorrPiMinus; }
inline float AliFemtoTrack::CorrectionKaonMinus() const { return fCorrKMinus; }
inline float AliFemtoTrack::CorrectionProtonMinus() const { return fCorrPMinus; }
inline float AliFemtoTrack::CorrectionDeuteronMinus() const { return fCorrDMinus; }
inline float AliFemtoTrack::CorrectionTritonMinus() const { return fCorrTMinus; }
inline float AliFemtoTrack::CorrectionHe3Minus() const { return fCorrHMinus; }
inline float AliFemtoTrack::CorrectionAlphaMinus() const { return fCorrAMinus; }
inline float AliFemtoTrack::CorrectionAll() const { return fCorrAll; }

inline double AliFemtoTrack::XatDCA() const { return fXatDCA; }
inline double AliFemtoTrack::YatDCA() const { return fYatDCA; }
inline double AliFemtoTrack::ZatDCA() const { return fZatDCA; }

inline void AliFemtoTrack::SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo) {fHiddenInfo=aHiddenInfo;}
inline bool AliFemtoTrack::ValidHiddenInfo() const { return (fHiddenInfo != NULL); }
inline AliFemtoHiddenInfo* AliFemtoTrack::GetHiddenInfo() const { return fHiddenInfo; }

inline AliFemtoThreeVector* AliFemtoTrack::GetTrueMomentum() const { return fTrueMomentum; }
inline const AliFemtoLorentzVector* AliFemtoTrack::GetEmissionPoint() const { return fEmissionPoint; }
inline const AliFemtoThreeVector* AliFemtoTrack::GetGlobalEmissionPoint() const { return fGlobalEmissionPoint; }

inline void AliFemtoTrack::SetPDGPid(Int_t aPid) { fPDGPid = aPid; }
inline Int_t AliFemtoTrack::GetPDGPid() const { return fPDGPid; }

inline void AliFemtoTrack::SetMass(Double_t aMass) { fMass = aMass; }
inline Double_t AliFemtoTrack::GetMass() const { return fMass; }


inline void AliFemtoTrack::SetPrimaryVertex(const double* vertex)
{
  fVertex[0] = vertex[0];
  fVertex[1] = vertex[1];
  fVertex[2] = vertex[2];
}


inline void AliFemtoTrack::GetPrimaryVertex(double* vertex) const
{
  vertex[0] = fVertex[0];
  vertex[1] = fVertex[1];
  vertex[2] = fVertex[2];
}

inline void AliFemtoTrack::SetKinkIndexes(const int points[3])
{
  // Transfer the Kink indices
  fKinkIndexes[0] = points[0];
  fKinkIndexes[1] = points[1];
  fKinkIndexes[2] = points[2];
}

inline int AliFemtoTrack::KinkIndex(int aIndex) const
{
  // Return Kink index
  if ((aIndex < 3) && (aIndex >= 0))
    return fKinkIndexes[aIndex];
  else
    return 0;
}

inline const AliFemtoThreeVector& AliFemtoTrack::NominalTpcExitPoint() const { return fNominalTpcExitPoint; }
inline const AliFemtoThreeVector& AliFemtoTrack::NominalTpcPointShifted() const { return fNominalTpcPointShifted; }
inline const AliFemtoThreeVector& AliFemtoTrack::NominalTpcEntrancePoint() const { return fNominalTpcEntrancePoint; }

inline const AliFemtoThreeVector& AliFemtoTrack::NominalTpcPoint(int i) const
{
  if(i<0)
    return fNominalTpcPoints[0];
  if(i>8)
    return fNominalTpcPoints[8];
  return fNominalTpcPoints[i];
}

inline const TBits& AliFemtoTrack::TPCclusters() const {return fClusters;}
inline const TBits& AliFemtoTrack::TPCsharing()  const {return fShared;}

inline void AliFemtoTrack::SetTPCcluster(const short& aNBit, const Bool_t& aValue) { fClusters.SetBitNumber(aNBit, aValue); }
inline void AliFemtoTrack::SetTPCshared(const short& aNBit, const Bool_t& aValue) { fShared.SetBitNumber(aNBit, aValue); }

inline void AliFemtoTrack::SetTPCClusterMap(const TBits& aBits) { fClusters = aBits; }
inline void AliFemtoTrack::SetTPCSharedMap(const TBits& aBits) { fShared = aBits; }


inline void AliFemtoTrack::SetITSHitOnLayer(int i, bool val)
{
  // Transfer ITS hit
  fHasPointOnITS[i] = val;
}


inline bool AliFemtoTrack::HasPointOnITSLayer(int aIndex) const
{
  // Return if i-th ITS layer had a hit for this track
  if ((aIndex < 6) && (aIndex >= 0))
    return fHasPointOnITS[aIndex];
  else
    return false;
}

inline void AliFemtoTrack::SetNominalTPCEntrancePoint(const AliFemtoThreeVector& aXTPC) { fNominalTpcEntrancePoint = aXTPC; }
inline void AliFemtoTrack::SetNominalTPCEntrancePoint(double *aXTPC)
{
  // Store the nominal TPC entrance point
  fNominalTpcEntrancePoint.SetX(aXTPC[0]);
  fNominalTpcEntrancePoint.SetY(aXTPC[1]);
  fNominalTpcEntrancePoint.SetZ(aXTPC[2]);
}

inline void AliFemtoTrack::SetNominalTPCPoints(const AliFemtoThreeVector* const points)
{
  std::copy_n(points, 9, fNominalTpcPoints);
}
inline void AliFemtoTrack::SetNominalTPCPoints(double **aXTPC)
{
  // Store the nominal TPC points
  for(int i=0;i<9;i++)
    {
      fNominalTpcPoints[i].SetX(aXTPC[i][0]);
      fNominalTpcPoints[i].SetY(aXTPC[i][1]);
      fNominalTpcPoints[i].SetZ(aXTPC[i][2]);
    }
}

inline void AliFemtoTrack::SetNominalTPCExitPoint(const AliFemtoThreeVector& aXTPC)
{
  fNominalTpcExitPoint = aXTPC;
}
inline void AliFemtoTrack::SetNominalTPCExitPoint(double *aXTPC)
{
  // Store the nominal TPC exit point
  fNominalTpcExitPoint.SetX(aXTPC[0]);
  fNominalTpcExitPoint.SetY(aXTPC[1]);
  fNominalTpcExitPoint.SetZ(aXTPC[2]);
}

inline void AliFemtoTrack::SetNominalTPCPointShifted(const AliFemtoThreeVector& aXTPC)
{
  fNominalTpcPointShifted = aXTPC;
}

inline void AliFemtoTrack::SetNominalTPCPointShifted(double *aXTPC)
{
  fNominalTpcPointShifted.SetX(aXTPC[0]);
  fNominalTpcPointShifted.SetY(aXTPC[1]);
  fNominalTpcPointShifted.SetZ(aXTPC[2]);
}

inline void AliFemtoTrack::SetTrueMomentum(AliFemtoThreeVector *aMom)
{
  // Set momentum from vector
  if (fTrueMomentum) {
    *fTrueMomentum = *aMom;
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector(*aMom);
  }
}

inline void AliFemtoTrack::SetTrueMomentum(const AliFemtoThreeVector& aMom)
{
  // Set momentum from vector
  if (fTrueMomentum) {
    *fTrueMomentum = aMom;
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector(aMom);
  }
}

inline void AliFemtoTrack::SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz)
{
  // Set momentum from components
  if (fTrueMomentum) {
    fTrueMomentum->SetX(aPx);
    fTrueMomentum->SetY(aPy);
    fTrueMomentum->SetZ(aPz);
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector(aPx, aPy, aPz);
  }
}

inline void AliFemtoTrack::SetEmissionPoint(AliFemtoLorentzVector *aPos)
{
  // Set position from vector
  if (fEmissionPoint) {
    *fEmissionPoint = *aPos;
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(*aPos);
  }
}

inline void AliFemtoTrack::SetEmissionPoint(const AliFemtoLorentzVector& aPos)
{
  // Set position from vector
  if (fEmissionPoint) {
    *fEmissionPoint = aPos;
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(aPos);
  }
}

inline void AliFemtoTrack::SetEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT)
{
  // Set position from components
  if (fEmissionPoint) {
    fEmissionPoint->SetX(aRx);
    fEmissionPoint->SetY(aRy);
    fEmissionPoint->SetZ(aRz);
    fEmissionPoint->SetT(aT);
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(aRx, aRy, aRz, aT);
  }
}

inline void AliFemtoTrack::SetGlobalEmissionPoint(const AliFemtoThreeVector& aPos)
{
  // set position from vector
  if (fGlobalEmissionPoint) {
    *fGlobalEmissionPoint = aPos;
  } else {
    fGlobalEmissionPoint = new AliFemtoThreeVector(aPos);
  }
}

inline void AliFemtoTrack::SetGlobalEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz)
{
  // Set position from components
  if (fGlobalEmissionPoint) {
    fGlobalEmissionPoint->SetX(aRx);
    fGlobalEmissionPoint->SetY(aRy);
    fGlobalEmissionPoint->SetZ(aRz);
  }
  else {
    fGlobalEmissionPoint = new AliFemtoThreeVector(aRx, aRy, aRz);
  }
}

inline void AliFemtoTrack::SetTofExpectedTimes(const float& tpi, const float& tkn, const float& tpr, const float& ttof)
{
  fTofPionTime = tpi;
  fTofKaonTime = tkn;
  fTofProtonTime = tpr;
  fTofDeuteronTime = ttof;
  fTofTritonTime = ttof;
  fTofHe3Time = ttof;
  fTofAlphaTime = ttof;
}

#endif
