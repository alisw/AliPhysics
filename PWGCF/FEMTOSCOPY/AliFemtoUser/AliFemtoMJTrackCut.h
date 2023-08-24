///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoMJTrackCut: A basic track cut that used information from     //
// ALICE ESD to accept or reject the track.                             ////////////////////////////////////////////////
#include "AliESDtrackCuts.h"

#ifndef ALIFEMTOMJTRACKCUT_H
#define ALIFEMTOMJTRACKCUT_H

//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliESDtrackCuts.h" //for enum with ITS layers
#include "AliFemtoTrackCut.h"


class AliFemtoMJTrackCut : public AliFemtoTrackCut
{
  public:

  enum PIDMethodType {knSigma=0, kContour=1};
  typedef enum PIDMethodType ReadPIDMethodType;

  AliFemtoMJTrackCut();
  virtual ~AliFemtoMJTrackCut();

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
  void SetStatus(const long& w);
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
  void SetMaxImpactZPtDep(const float& maxoff, const float& maxnrm, const float& maxpow);
  void SetMostProbablePion();
  void SetMostProbableKaon();
  void SetMostProbableProton();
  void SetMostProbableDeuteron();
  void SetLeastProbableProton();
  void SetNoMostProbable();
  void SetMostProbable(const int& num);
  void SetPIDMethod(ReadPIDMethodType newMethod);
  void SetNsigmaTPCTOF(Bool_t);
  void SetNsigmaTPConly(Bool_t);
  void SetNsigma(Double_t);
  void SetNsigma2(Double_t);
  void SetNsigmaRejection(Double_t);
  void SetNsigmaAccept(Double_t);
  void SetClusterRequirementITS(AliESDtrackCuts::Detector det, AliESDtrackCuts::ITSClusterRequirement req = AliESDtrackCuts::kOff);

  void SetMomRangeTOFpidIs(const float& minp, const float& maxp);
  void SetMomRangeTPCpidIs(const float& minp, const float& maxp);
  void SetMomRangeITSpidIs(const float& minp, const float& maxp);
  void SetElectronRejection(Bool_t);

 private:   // here are the quantities I want to cut on...

  int               fCharge;             // particle charge
  float             fPt[2];              // bounds for transverse momentum
  float             fRapidity[2];        // bounds for rapidity
  float             fEta[2];             // bounds for pseudorapidity
  float             fPidProbElectron[2]; // bounds for electron probability
  float             fPidProbPion[2];     // bounds for pion probability
  float             fPidProbKaon[2];     // bounds for kaon probability
  float             fPidProbProton[2];   // bounds for proton probability
  float             fPidProbMuon[2];     // bounds for muon probability

  AliESDtrackCuts::ITSClusterRequirement fCutClusterRequirementITS[3];  // detailed ITS cluster requirements for (SPD, SDD, SSD) - from AliESDtrackcuts!
  bool              fLabel;              // if true label<0 will not pass throught
  long              fStatus;             // staus flag
  ReadPIDMethodType fPIDMethod;          // which PID mehod to use. 0 - nsgima, 1 - contour
  Bool_t            fNsigmaTPCTOF;       // true if squared nsigma from TPC and TOF, false if separately from TPC and TOF
  Bool_t            fNsigmaTPConly;      // true if nsigma from TPC only
  Double_t          fNsigma;             // number of sigmas - 3 by default
  Double_t          fNsigma2;             // number of sigmas - 3 by default
  Double_t          fNsigmaRejection;     // number of sigmas for rejection - 3 by default
  Double_t          fNsigmaAccept;        // number of sigmas for rejection - 3 by default

  short             fminTPCclsF;         // min number of findable clusters in the TPC
  short             fminTPCncls;         // min number of clusters in the TPC
  int               fminITScls;          // min number of clusters assigned in the ITS
  float             fMaxITSchiNdof;      // maximum allowed chi2/ndof for ITS clusters
  float             fMaxTPCchiNdof;      // maximum allowed chi2/ndof for TPC clusters
  float             fMaxSigmaToVertex;   // maximum allowed sigma to primary vertex
  long              fNTracksPassed;      // passed tracks count
  long              fNTracksFailed;      // failed tracks count
  bool              fRemoveKinks;        // if true particles with any kink label will not pass
  bool              fRemoveITSFake;      // if true particles with ITS fake flag will not pass
  int               fMostProbable;       // this particle type is required to be most probable

  float             fMaxImpactXY;        // Max XY impact parameter
  float             fMinImpactXY;        // Max XY impact parameter
  float             fMaxImpactZ;         // Max Z impact parameter

  float             fMaxImpactXYPtOff;   // Max XY DCA Pt dependent offset
  float             fMaxImpactXYPtNrm;   // Max XY DCA Pt dependent normalization
  float             fMaxImpactXYPtPow;   // Max XY DCA Pt dependent power

  float             fMaxImpactZPtOff;   // Max Z DCA Pt dependent offset
  float             fMaxImpactZPtNrm;   // Max Z DCA Pt dependent normalization
  float             fMaxImpactZPtPow;   // Max Z DCA Pt dependent power
  
  float             fMinPforTOFpid;  // momentum from which TOF PID is requested
  float             fMaxPforTOFpid;  // momentum till which TOF PID is requested
  float             fMinPforTPCpid;  // momentum from which TPC PID is requested
  float             fMaxPforTPCpid;  // momentum till which TPC PID is requested
  float             fMinPforITSpid;  // momentum from which ITS PID is requested
  float             fMaxPforITSpid;  // momentum till which ITS PID is requested
  bool fElectronRejection;

  float PidFractionElectron(float mom) const;
  float PidFractionPion(float mom) const;
  float PidFractionKaon(float mom) const;
  float PidFractionProton(float mom) const;

  bool IsPionTPCdEdx(float mom, float dEdx);
  bool IsKaonTPCdEdx(float mom, float dEdx);
  bool IsProtonTPCdEdx(float mom, float dEdx);

  bool IsPionTOFTime(float mom, float ttof);
  bool IsKaonTOFTime(float mom, float ttof);
  bool IsProtonTOFTime(float mom, float ttof);

  bool IsKaonTPCdEdxNSigma(float mom, float nsigma);
  bool IsKaonTOFNSigma(float mom, float nsigma);
  bool IsKaonNSigma(float mom, float nsigmaTPC, float nsigmaTOF);
  bool IsPionNSigma(float mom, float nsigmaTPC, float nsigmaTOF);
  bool IsProtonNSigma(float mom, float nsigmaTPC, float nsigmaTOF);
  bool IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP);

  bool IsKaonNSigmaRejection(float mom, float nsigmaTPC, float nsigmaTOF, float TOFtime);
  bool IsPionNSigmaRejection(float mom, float nsigmaTPC, float nsigmaTOF, float TOFtime);
  bool IsProtonNSigmaRejection(float mom, float nsigmaTPC, float nsigmaTOF, float TOFtime);

  bool IsKaonNSigmaAccept(float mom, float nsigmaTPC, float nsigmaTOF, float TOFtime);
  bool IsPionNSigmaAccept(float mom, float nsigmaTPC, float nsigmaTOF, float TOFtime);
  bool IsProtonNSigmaAccept(float mom, float nsigmaTPC, float nsigmaTOF, float TOFtime);

  bool IsDeuteronTPCNSigma(float mom, float nsigmaTPC);



  Bool_t CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2); //the same as in AliESDtrackCuts


#ifdef __ROOT__
  ClassDef(AliFemtoMJTrackCut, 1)
#endif
    };


inline void AliFemtoMJTrackCut::SetPt(const float& lo, const float& hi){fPt[0]=lo; fPt[1]=hi;}
inline void AliFemtoMJTrackCut::SetRapidity(const float& lo,const float& hi){fRapidity[0]=lo; fRapidity[1]=hi;}
inline void AliFemtoMJTrackCut::SetEta(const float& lo,const float& hi){fEta[0]=lo; fEta[1]=hi;}
inline void AliFemtoMJTrackCut::SetCharge(const int& ch){fCharge = ch;}
inline void AliFemtoMJTrackCut::SetPidProbElectron(const float& lo,const float& hi){fPidProbElectron[0]=lo; fPidProbElectron[1]=hi;}
inline void AliFemtoMJTrackCut::SetPidProbPion(const float& lo,const float& hi){fPidProbPion[0]=lo; fPidProbPion[1]=hi;}
inline void AliFemtoMJTrackCut::SetPidProbKaon(const float& lo,const float& hi){fPidProbKaon[0]=lo; fPidProbKaon[1]=hi;}
inline void AliFemtoMJTrackCut::SetPidProbProton
(const float& lo,const float& hi){fPidProbProton[0]=lo; fPidProbProton[1]=hi;}
inline void AliFemtoMJTrackCut::SetPidProbMuon(const float& lo,const float& hi){fPidProbMuon[0]=lo; fPidProbMuon[1]=hi;}
inline void AliFemtoMJTrackCut::SetLabel(const bool& flag){fLabel=flag;}
inline void AliFemtoMJTrackCut::SetStatus(const long& status){fStatus=status;}
inline void AliFemtoMJTrackCut::SetminTPCclsF(const short& minTPCclsF){fminTPCclsF=minTPCclsF;}
inline void AliFemtoMJTrackCut::SetminTPCncls(const short& s){fminTPCncls=s;}
inline void AliFemtoMJTrackCut::SetminITScls(const int& minITScls){fminITScls=minITScls;}
inline void AliFemtoMJTrackCut::SetMostProbablePion() { fMostProbable = 2; }
inline void AliFemtoMJTrackCut::SetMostProbableKaon() { fMostProbable = 3; }
inline void AliFemtoMJTrackCut::SetMostProbableProton() { fMostProbable = 4; }
inline void AliFemtoMJTrackCut::SetMostProbableDeuteron() { fMostProbable = 30; }
inline void AliFemtoMJTrackCut::SetLeastProbableProton() { fMostProbable = 5; }
inline void AliFemtoMJTrackCut::SetNoMostProbable() { fMostProbable = 0; }
inline void AliFemtoMJTrackCut::SetMostProbable(const int& num) {  fMostProbable =  num; }
inline void AliFemtoMJTrackCut::SetMaxITSChiNdof(const float& maxchi) { fMaxITSchiNdof = maxchi; }
inline void AliFemtoMJTrackCut::SetMaxTPCChiNdof(const float& maxchi) { fMaxTPCchiNdof = maxchi; }
inline void AliFemtoMJTrackCut::SetMaxSigmaToVertex(const float& maxsig) { fMaxSigmaToVertex = maxsig; }
inline void AliFemtoMJTrackCut::SetMaxImpactXY(const float& maximpxy) { fMaxImpactXY = maximpxy; }
inline void AliFemtoMJTrackCut::SetMinImpactXY(const float& minimpxy) { fMinImpactXY = minimpxy; }
inline void AliFemtoMJTrackCut::SetMaxImpactXYPtDep(const float& maxoff, const float& maxnrm, const float& maxpow) { fMaxImpactXYPtOff = maxoff; fMaxImpactXYPtNrm = maxnrm; fMaxImpactXYPtPow = maxpow; }
inline void AliFemtoMJTrackCut::SetMaxImpactZ(const float& maximpz) { fMaxImpactZ = maximpz; }
inline void AliFemtoMJTrackCut::SetMaxImpactZPtDep(const float& maxoff, const float& maxnrm, const float& maxpow) { fMaxImpactZPtOff = maxoff; fMaxImpactZPtNrm = maxnrm; fMaxImpactZPtPow = maxpow; }
inline void AliFemtoMJTrackCut::SetElectronRejection(Bool_t setE) { fElectronRejection = setE; }

#endif
