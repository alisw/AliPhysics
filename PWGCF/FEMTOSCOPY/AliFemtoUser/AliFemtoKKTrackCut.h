///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoESDTrackCut: A basic track cut that used information from     //
// ALICE ESD to accept or reject the track.                              //  
// Enables the selection on charge, transverse momentum, rapidity,       //
// pid probabilities, number of ITS and TPC clusters                     //
// Author: Marek Chojnacki (WUT), mchojnacki@knf.pw.edu.pl               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#include "AliESDtrackCuts.h"

#ifndef ALIFEMTOKKTRACKCUT_H
#define ALIFEMTOKKTRACKCUT_H

//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliESDtrackCuts.h" //for enum with ITS layers
#include "AliFemtoTrackCut.h"


class AliFemtoKKTrackCut : public AliFemtoTrackCut 
{
  public:

  enum PIDMethodType {knSigma=0, kContour=1};
  typedef enum PIDMethodType ReadPIDMethodType; 

  AliFemtoKKTrackCut();
  virtual ~AliFemtoKKTrackCut();

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
  void SetMaxImpactZ(const float& maximpz);
  void SetMaxImpactXYPtDep(const float& maxoff, const float& maxnrm, const float& maxpow);
  void SetMostProbablePion();
  void SetMostProbableKaon();
  void SetMostProbableProton();
  void SetNoMostProbable(); 
  void SetPIDMethod(ReadPIDMethodType newMethod);
  void SetClusterRequirementITS(AliESDtrackCuts::Detector det, AliESDtrackCuts::ITSClusterRequirement req = AliESDtrackCuts::kOff);

  void SetMomRangeTOFpidIs(const float& minp, const float& maxp);
  void SetMomRangeTPCpidIs(const float& minp, const float& maxp);
  void SetMomRangeITSpidIs(const float& minp, const float& maxp);

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
  float             fMaxImpactZ;         // Max Z impact parameter

  float             fMaxImpactXYPtOff;   // Max XY DCA Pt dependent offset
  float             fMaxImpactXYPtNrm;   // Max XY DCA Pt dependent normalization
  float             fMaxImpactXYPtPow;   // Max XY DCA Pt dependent power

  float             fMinPforTOFpid;  // momentum from which TOF PID is requested
  float             fMaxPforTOFpid;  // momentum till which TOF PID is requested
  float             fMinPforTPCpid;  // momentum from which TPC PID is requested
  float             fMaxPforTPCpid;  // momentum till which TPC PID is requested
  float             fMinPforITSpid;  // momentum from which ITS PID is requested
  float             fMaxPforITSpid;  // momentum till which ITS PID is requested

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

  Bool_t CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2); //the same as in AliESDtrackCuts


#ifdef __ROOT__ 
  ClassDef(AliFemtoKKTrackCut, 1)
#endif
    };


inline void AliFemtoKKTrackCut::SetPt(const float& lo, const float& hi){fPt[0]=lo; fPt[1]=hi;}
inline void AliFemtoKKTrackCut::SetRapidity(const float& lo,const float& hi){fRapidity[0]=lo; fRapidity[1]=hi;}
inline void AliFemtoKKTrackCut::SetEta(const float& lo,const float& hi){fEta[0]=lo; fEta[1]=hi;}
inline void AliFemtoKKTrackCut::SetCharge(const int& ch){fCharge = ch;}
inline void AliFemtoKKTrackCut::SetPidProbElectron(const float& lo,const float& hi){fPidProbElectron[0]=lo; fPidProbElectron[1]=hi;}
inline void AliFemtoKKTrackCut::SetPidProbPion(const float& lo,const float& hi){fPidProbPion[0]=lo; fPidProbPion[1]=hi;}
inline void AliFemtoKKTrackCut::SetPidProbKaon(const float& lo,const float& hi){fPidProbKaon[0]=lo; fPidProbKaon[1]=hi;}
inline void AliFemtoKKTrackCut::SetPidProbProton(const float& lo,const float& hi){fPidProbProton[0]=lo; fPidProbProton[1]=hi;}
inline void AliFemtoKKTrackCut::SetPidProbMuon(const float& lo,const float& hi){fPidProbMuon[0]=lo; fPidProbMuon[1]=hi;}
inline void AliFemtoKKTrackCut::SetLabel(const bool& flag){fLabel=flag;}
inline void AliFemtoKKTrackCut::SetStatus(const long& status){fStatus=status;}
inline void AliFemtoKKTrackCut::SetminTPCclsF(const short& minTPCclsF){fminTPCclsF=minTPCclsF;}
inline void AliFemtoKKTrackCut::SetminTPCncls(const short& s){fminTPCncls=s;}
inline void AliFemtoKKTrackCut::SetminITScls(const int& minITScls){fminITScls=minITScls;}
inline void AliFemtoKKTrackCut::SetMostProbablePion() { fMostProbable = 2; }
inline void AliFemtoKKTrackCut::SetMostProbableKaon() { fMostProbable = 3; }
inline void AliFemtoKKTrackCut::SetMostProbableProton() { fMostProbable = 4; }
inline void AliFemtoKKTrackCut::SetNoMostProbable() { fMostProbable = 0; }
inline void AliFemtoKKTrackCut::SetMaxITSChiNdof(const float& maxchi) { fMaxITSchiNdof = maxchi; }
inline void AliFemtoKKTrackCut::SetMaxTPCChiNdof(const float& maxchi) { fMaxTPCchiNdof = maxchi; }
inline void AliFemtoKKTrackCut::SetMaxSigmaToVertex(const float& maxsig) { fMaxSigmaToVertex = maxsig; }
inline void AliFemtoKKTrackCut::SetMaxImpactXY(const float& maximpxy) { fMaxImpactXY = maximpxy; }
inline void AliFemtoKKTrackCut::SetMaxImpactXYPtDep(const float& maxoff, const float& maxnrm, const float& maxpow) { fMaxImpactXYPtOff = maxoff; fMaxImpactXYPtNrm = maxnrm; fMaxImpactXYPtPow = maxpow; }
inline void AliFemtoKKTrackCut::SetMaxImpactZ(const float& maximpz) { fMaxImpactZ = maximpz; }


#endif

