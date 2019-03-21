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

  AliFemtoThreeVector P() const;
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
  AliFemtoLorentzVector *GetEmissionPoint() const;
  Int_t                  GetPDGPid() const;
  Double_t               GetMass() const;

  AliFemtoThreeVector   *GetGlobalEmissionPoint() const;
  void                   SetGlobalEmissionPoint(const AliFemtoThreeVector& aPos);
  void                   SetGlobalEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz);


  void SetPrimaryVertex(const double *vertex);
  void GetPrimaryVertex(double *fisvertex);

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
//{return &mNSigmaElectron;} // Fab private
inline float AliFemtoTrack::PidProbElectron() const {return fPidProbElectron;}
inline float AliFemtoTrack::PidProbPion() const {return fPidProbPion;}
inline float AliFemtoTrack::PidProbKaon() const {return fPidProbKaon;}
inline float AliFemtoTrack::PidProbProton() const {return fPidProbProton;}
inline float AliFemtoTrack::PidProbMuon() const {return fPidProbMuon;}
inline int AliFemtoTrack::Multiplicity() const { return fMultiplicity;}
inline double AliFemtoTrack::Zvtx() const { return fZvtx;}

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

#endif
