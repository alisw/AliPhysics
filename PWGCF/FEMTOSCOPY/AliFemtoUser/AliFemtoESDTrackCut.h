///
/// \file AliFemtoESDTrackCut.h
///


#ifndef ALIFEMTOESDTRACKCUT_H
#define ALIFEMTOESDTRACKCUT_H

#include "AliESDtrackCuts.h" //for enum with ITS layers
#include "AliFemtoTrackCut.h"

#include <utility>


/// \class AliFemtoESDTrackCut
/// \brief A basic track cut that used information from ALICE ESD to accept or reject the track.
///
class AliFemtoESDTrackCut : public AliFemtoTrackCut {
public:

  enum PIDMethodType {knSigma=0, kContour=1};
  typedef enum PIDMethodType ReadPIDMethodType;
  AliFemtoESDTrackCut();
  virtual ~AliFemtoESDTrackCut();

  virtual bool Pass(const AliFemtoTrack* aTrack);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtTrack;}

  void SetPt(const float& lo, const float& hi);
  void SetRapidity(const float& lo, const float& hi);
  void SetEta(const float& lo, const float& hi);
  void SetCharge(const int& ch);
  void SetPidProbElectron(const float& lo, const float& hi);
  void SetPidProbPion(const float& lo, const float& hi);
  void SetPidProbKaon(const float& lo, const float& hi);
  void SetPidProbProton(const float& lo, const float& hi);
  void SetPidProbMuon(const float& lo, const float& hi);

  void SetLabel(const bool& flag);
  void SetStatus(const ULong64_t w);
  void SetminTPCclsF(const short& s);
  void SetminTPCncls(const short& s);
  void SetminITScls(const int& s);
  void SetRemoveKinks(const bool& flag);
  void SetRemoveITSFake(const bool& flag);
  void SetMaxITSChiNdof(const float& maxchi);
  void SetMaxTPCChiNdof(const float& maxchi);
  void SetMaxSigmaToVertex(const float& maxsig);
  void SetMaxImpactXY(const float& maximpxy);
  void SetMinImpactXY(const float& minimpxy);
  void SetMaxImpactZ(const float& maximpz);
  void SetMaxImpactXYPtDep(const float& maxoff, const float& maxnrm, const float& maxpow);
  void SetMostProbablePion();
  void SetMostProbableKaon();
  void SetMostProbableProton();
  void SetLeastProbableProton();
  void SetNoMostProbable();
   //
  void SetMostProbableDeuteron();
  void SetMostProbableTriton();
  void SetMostProbableHe3();
  void SetMostProbableAlpha();
  //

  void SetMostProbable(const int& num);
  void SetPIDMethod(ReadPIDMethodType newMethod);
  void SetNsigmaTPCTOF(Bool_t);
  void SetNsigmaTPConly(Bool_t);
  void SetNsigma(Double_t);
  void SetNsigmaMass(Double_t);
  void SetClusterRequirementITS(AliESDtrackCuts::Detector det, AliESDtrackCuts::ITSClusterRequirement req = AliESDtrackCuts::kOff);

  void SetMomRangeTOFpidIs(const float& minp, const float& maxp);
  void SetMomRangeTPCpidIs(const float& minp, const float& maxp);
  void SetMomRangeITSpidIs(const float& minp, const float& maxp);
  void SetElectronRejection(Bool_t);

  int GetCharge() const { return fCharge; }
  std::pair<float, float> GetPt() const { return std::make_pair(fPt[0], fPt[1]); }
  std::pair<float, float> GetRapidity() const { return std::make_pair(fRapidity[0], fRapidity[1]); }
  std::pair<float, float> GetEta() const { return std::make_pair(fEta[0], fEta[1]); }
  std::pair<float, float> GetProbElectron() const { return std::make_pair(fPidProbElectron[0], fPidProbElectron[1]); }
  std::pair<float, float> GetProbPion() const { return std::make_pair(fPidProbPion[0], fPidProbPion[1]); }
  std::pair<float, float> GetProbKaon() const { return std::make_pair(fPidProbKaon[0], fPidProbKaon[1]); }
  std::pair<float, float> GetProbProton() const { return std::make_pair(fPidProbProton[0], fPidProbProton[1]); }
  std::pair<float, float> GetProbMuon() const { return std::make_pair(fPidProbMuon[0], fPidProbMuon[1]); }

  bool GetLabel() const { return fLabel; }
  ULong64_t GetStatus() const { return fStatus; }
  int GetPIDmethod() const { return fPIDMethod; }

  int GetMinFindableClustersTPC() const { return fminTPCclsF; }
  int GetMinNClustersTPC() const { return fminTPCncls; }
  int GetMinNClustersITS() const { return fminITScls; }

  float GetMaxITSchiNdof() const { return fMaxITSchiNdof; }
  float GetMaxTPCchiNdof() const { return fMaxTPCchiNdof; }
  float GetMaxSigmaToVertex() const { return fMaxSigmaToVertex; }



  /// Use TPC & TOF information
  bool GetDualNsigma() const { return fNsigmaTPCTOF; }
  bool GetNsigmaTPConly() const { return fNsigmaTPConly; }
  double GetNsigma() const { return fNsigma; }
  double GetNsigmaMass() const { return fNsigmaMass; }
  bool GetRemoveKinks() const { return fRemoveKinks; }
  bool GetRemoveITSFake() const { return fRemoveITSFake; }
  int GetMostProbable() const { return fMostProbable; }

  float GetMinImpactXY() const { return fMinImpactXY; }
  float GetMaxImpactXY() const { return fMaxImpactXY; }
  float GetMaxImpactZ() const { return fMaxImpactZ; }

  float GetMaxImpactXyPtOff() const { return fMaxImpactXYPtOff; }   ///< Max XY DCA Pt dependent offset
  float GetMaxImpactXYPtNrm() const { return fMaxImpactXYPtNrm; }  ///< Max XY DCA Pt dependent normalization
  float GetMaxImpactXYPtPow() const { return fMaxImpactXYPtPow; }  ///< Max XY DCA Pt dependent power

  std::pair<float, float> GetTOFpidMomentumRange() const { return std::make_pair(fMinPforTOFpid, fMaxPforTOFpid); }
  std::pair<float, float> GetTPCpidMomentumRange() const { return std::make_pair(fMinPforTPCpid, fMaxPforTPCpid); }
  std::pair<float, float> GetITSpidMomentumRange() const { return std::make_pair(fMinPforITSpid, fMaxPforITSpid); }

  bool GetElectronRejection() const { return fElectronRejection; }

protected:   // here are the quantities I want to cut on...

  int               fCharge;             ///< particle charge
  float             fPt[2];              ///< bounds for transverse momentum
  float             fRapidity[2];        ///< bounds for rapidity
  float             fEta[2];             ///< bounds for pseudorapidity
  float             fPidProbElectron[2]; ///< bounds for electron probability
  float             fPidProbPion[2];     ///< bounds for pion probability
  float             fPidProbKaon[2];     ///< bounds for kaon probability
  float             fPidProbProton[2];   ///< bounds for proton probability
  float             fPidProbMuon[2];     ///< bounds for muon probability



  
  AliESDtrackCuts::ITSClusterRequirement fCutClusterRequirementITS[3];  ///< detailed ITS cluster requirements for (SPD, SDD, SSD) - from AliESDtrackcuts!
  bool              fLabel;              ///< if true label<0 will not pass throught
  ULong64_t         fStatus;             ///< staus flag
  ReadPIDMethodType fPIDMethod;          ///< which PID mehod to use. 0 - nsgima, 1 - contour
  Bool_t            fNsigmaTPCTOF;       ///< true if squared nsigma from TPC and TOF, false if separately from TPC and TOF
  Bool_t            fNsigmaTPConly;      ///< true if nsigma from TPC only
  Double_t          fNsigma;             ///< number of sigmas - 3 by default
  Double_t          fNsigmaMass;             ///< number of sigmas in mass condition (deuterons)

  short             fminTPCclsF;         ///< min number of findable clusters in the TPC
  short             fminTPCncls;         ///< min number of clusters in the TPC
  int               fminITScls;          ///< min number of clusters assigned in the ITS
  float             fMaxITSchiNdof;      ///< maximum allowed chi2/ndof for ITS clusters
  float             fMaxTPCchiNdof;      ///< maximum allowed chi2/ndof for TPC clusters
  float             fMaxSigmaToVertex;   ///< maximum allowed sigma to primary vertex
  long              fNTracksPassed;      ///< passed tracks count
  long              fNTracksFailed;      ///< failed tracks count
  bool              fRemoveKinks;        ///< if true particles with any kink label will not pass
  bool              fRemoveITSFake;      ///< if true particles with ITS fake flag will not pass
  int               fMostProbable;       ///< this particle type is required to be most probable

  float             fMaxImpactXY;        ///< Max XY impact parameter
  float             fMinImpactXY;        ///< Max XY impact parameter
  float             fMaxImpactZ;         ///< Max Z impact parameter

  float             fMaxImpactXYPtOff;   ///< Max XY DCA Pt dependent offset
  float             fMaxImpactXYPtNrm;   ///< Max XY DCA Pt dependent normalization
  float             fMaxImpactXYPtPow;   ///< Max XY DCA Pt dependent power

  float             fMinPforTOFpid;  ///< momentum from which TOF PID is requested
  float             fMaxPforTOFpid;  ///< momentum till which TOF PID is requested
  float             fMinPforTPCpid;  ///< momentum from which TPC PID is requested
  float             fMaxPforTPCpid;  ///< momentum till which TPC PID is requested
  float             fMinPforITSpid;  ///< momentum from which ITS PID is requested
  float             fMaxPforITSpid;  ///< momentum till which ITS PID is requested
  bool fElectronRejection;

  float PidFractionElectron(float mom) const;
  float PidFractionPion(float mom) const;
  float PidFractionKaon(float mom) const;
  float PidFractionProton(float mom) const;

  
  bool IsPionTPCdEdx(float mom, float dEdx);
  bool IsKaonTPCdEdx(float mom, float dEdx);
  bool IsProtonTPCdEdx(float mom, float dEdx);
  bool IsDeuteronTPCdEdx(float mom, float dEdx);

  
  bool IsPionTOFTime(float mom, float ttof);
  bool IsKaonTOFTime(float mom, float ttof);
  bool IsProtonTOFTime(float mom, float ttof);

  //
  bool IsDeuteronTOFTime(float mom, float ttof);
  bool IsTritonTOFTime(float mom, float ttof);
  bool IsHe3TOFTime(float mom, float ttof);
  bool IsAlphaTOFTime(float mom, float ttof);
  //



  
  bool IsKaonTPCdEdxNSigma(float mom, float nsigma);
  bool IsKaonTOFNSigma(float mom, float nsigma);
  bool IsKaonNSigma(float mom, float nsigmaTPC, float nsigmaTOF);
  bool IsPionNSigma(float mom, float nsigmaTPC, float nsigmaTOF);
  bool IsProtonNSigma(float mom, float nsigmaTPC, float nsigmaTOF);
  bool IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP);
  //
  bool IsDeuteronNSigma(float mom,float massTOFDPG, float nsigmaTPC, float nsigmaTOF);
  bool IsTritonNSigma(float mom, float nsigmaTPC, float nsigmaTOF);
  bool IsHe3NSigma(float mom, float nsigmaTPC, float nsigmaTOF);
  bool IsAlphaNSigma(float mom, float nsigmaTPC, float nsigmaTOF);
  //



  
  Bool_t CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2); //the same as in AliESDtrackCuts


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoESDTrackCut, 3);
  /// \endcond
#endif
};


inline void AliFemtoESDTrackCut::SetPt(const float& lo, const float& hi){fPt[0]=lo; fPt[1]=hi;}
inline void AliFemtoESDTrackCut::SetRapidity(const float& lo,const float& hi){fRapidity[0]=lo; fRapidity[1]=hi;}
inline void AliFemtoESDTrackCut::SetEta(const float& lo,const float& hi){fEta[0]=lo; fEta[1]=hi;}
inline void AliFemtoESDTrackCut::SetCharge(const int& ch){fCharge = ch;}
inline void AliFemtoESDTrackCut::SetPidProbElectron(const float& lo,const float& hi){fPidProbElectron[0]=lo; fPidProbElectron[1]=hi;}
inline void AliFemtoESDTrackCut::SetPidProbPion(const float& lo,const float& hi){fPidProbPion[0]=lo; fPidProbPion[1]=hi;}
inline void AliFemtoESDTrackCut::SetPidProbKaon(const float& lo,const float& hi){fPidProbKaon[0]=lo; fPidProbKaon[1]=hi;}
inline void AliFemtoESDTrackCut::SetPidProbProton(const float& lo,const float& hi){fPidProbProton[0]=lo; fPidProbProton[1]=hi;}
inline void AliFemtoESDTrackCut::SetPidProbMuon(const float& lo,const float& hi){fPidProbMuon[0]=lo; fPidProbMuon[1]=hi;}


inline void AliFemtoESDTrackCut::SetLabel(const bool& flag){fLabel=flag;}
inline void AliFemtoESDTrackCut::SetStatus(const ULong64_t status){fStatus=status;}
inline void AliFemtoESDTrackCut::SetminTPCclsF(const short& minTPCclsF){fminTPCclsF=minTPCclsF;}
inline void AliFemtoESDTrackCut::SetminTPCncls(const short& s){fminTPCncls=s;}
inline void AliFemtoESDTrackCut::SetminITScls(const int& minITScls){fminITScls=minITScls;}
inline void AliFemtoESDTrackCut::SetMostProbablePion() { fMostProbable = 2; }
inline void AliFemtoESDTrackCut::SetMostProbableKaon() { fMostProbable = 3; }
inline void AliFemtoESDTrackCut::SetMostProbableProton() { fMostProbable = 4; }
inline void AliFemtoESDTrackCut::SetLeastProbableProton() { fMostProbable = 5; }
//
inline void AliFemtoESDTrackCut::SetMostProbableDeuteron() { fMostProbable = 13; }
inline void AliFemtoESDTrackCut::SetMostProbableTriton() { fMostProbable = 14; }
inline void AliFemtoESDTrackCut::SetMostProbableHe3() { fMostProbable = 15; }
inline void AliFemtoESDTrackCut::SetMostProbableAlpha() { fMostProbable = 16; }
//

inline void AliFemtoESDTrackCut::SetNoMostProbable() { fMostProbable = 0; }
inline void AliFemtoESDTrackCut::SetMostProbable(const int& num) {  fMostProbable =  num; }
inline void AliFemtoESDTrackCut::SetMaxITSChiNdof(const float& maxchi) { fMaxITSchiNdof = maxchi; }
inline void AliFemtoESDTrackCut::SetMaxTPCChiNdof(const float& maxchi) { fMaxTPCchiNdof = maxchi; }
inline void AliFemtoESDTrackCut::SetMaxSigmaToVertex(const float& maxsig) { fMaxSigmaToVertex = maxsig; }
inline void AliFemtoESDTrackCut::SetMaxImpactXY(const float& maximpxy) { fMaxImpactXY = maximpxy; }
inline void AliFemtoESDTrackCut::SetMinImpactXY(const float& minimpxy) { fMinImpactXY = minimpxy; }
inline void AliFemtoESDTrackCut::SetMaxImpactXYPtDep(const float& maxoff, const float& maxnrm, const float& maxpow) { fMaxImpactXYPtOff = maxoff; fMaxImpactXYPtNrm = maxnrm; fMaxImpactXYPtPow = maxpow; }
inline void AliFemtoESDTrackCut::SetMaxImpactZ(const float& maximpz) { fMaxImpactZ = maximpz; }
inline void AliFemtoESDTrackCut::SetElectronRejection(Bool_t setE) { fElectronRejection = setE; }

#endif


