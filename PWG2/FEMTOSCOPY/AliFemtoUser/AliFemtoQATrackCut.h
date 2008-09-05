///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoQATrackCut: A basic track cut that used information from     //
// ALICE ESD to accept or reject the track.                              //  
// Enables the selection on charge, transverse momentum, rapidity,       //
// pid probabilities, number of ITS and TPC clusters                     //
// Author: Marek Chojnacki (WUT), mchojnacki@knf.pw.edu.pl               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoQATrackCut_H
#define AliFemtoQATrackCut_H

//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoTrackCut.h"

class AliFemtoQATrackCut : public AliFemtoTrackCut 
{

 public:
  AliFemtoQATrackCut();
  virtual ~AliFemtoQATrackCut();

  virtual bool Pass(const AliFemtoTrack* aTrack);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtTrack;}

  void SetPt(const float& lo, const float& hi);
  void SetRapidity(const float& lo, const float& hi);
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
  void SetminTPCchiNdof(const float& s);
  void SetRemoveKinks(const bool& flag);
  void SetMaxTPCncls(const short& s);
  void SetMaxITSChiNdof(const float& maxchi);
  void SetMaxTPCChiNdof(const float& maxchi);
  void SetMaxSigmaToVertex(const float& maxsig);
  void SetMostProbablePion();
  void SetMostProbableKaon();
  void SetMostProbableProton();
  void SetNoMostProbable(); 
  
  void SetTPCnclsExclusionSwitch(const bool& excluSwitch);
  void SetTPCnclsExclusion(const short& lo, const short& hi);
  void SetTPCchiNdofExclusionSwitch(const bool& excluSwitch);
  void SetTPCchiNdofExclusion(const float& lo, const float& hi);

 private:   // here are the quantities I want to cut on...

  int               fCharge;             // particle charge
  float             fPt[2];              // bounds for transverse momentum
  float             fRapidity[2];        // bounds for rapidity
  float             fPidProbElectron[2]; // bounds for electron probability
  float             fPidProbPion[2];     // bounds for pion probability
  float             fPidProbKaon[2];     // bounds for kaon probability
  float             fPidProbProton[2];   // bounds for proton probability
  float             fPidProbMuon[2];     // bounds for muon probability 
  bool              fLabel;              // if true label<0 will not pass throught 
  long              fStatus;             // staus flag

  short             fminTPCclsF;         // min number of findable clusters in the TPC
  short             fminTPCncls;         // min number of clusters in the TPC
  int               fminITScls;          // min number of clusters assigned in the ITS 
  float             fminTPCchiNdof;      // min allowed chi2/ndof for TPC clusters
  short             fMaxTPCncls;         // maximum allowed clusters in the TPC
  float             fMaxITSchiNdof;      // maximum allowed chi2/ndof for ITS clusters
  float             fMaxTPCchiNdof;      // maximum allowed chi2/ndof for TPC clusters
  float             fMaxSigmaToVertex;   // maximum allowed sigma to primary vertex
  long              fNTracksPassed;      // passed tracks count
  long              fNTracksFailed;      // failed tracks count
  bool              fRemoveKinks;        // if true particles with any kink label will not pass
  int               fMostProbable;       // this particle type is required to be most probable
  
  bool    fTPCnclsExclusionSwitch;       // turn on/off TPCncls exclusion zone (true=on)
  short   fTPCnclsExclusion[2];          // lower and upper limit of TPCncls QA exclusion zone
  bool    fTPCchiNdofExclusionSwitch;        // turn on/off TPCchi exclusion zone (true=on)
  float   fTPCchiNdofExclusion[2];           // lower and upper limit of TPCchi QA exclusion zone

  float PidFractionElectron(float mom) const;
  float PidFractionPion(float mom) const;
  float PidFractionKaon(float mom) const;
  float PidFractionProton(float mom) const;

#ifdef __ROOT__ 
  ClassDef(AliFemtoQATrackCut, 1)
#endif
    };


inline void AliFemtoQATrackCut::SetPt(const float& lo, const float& hi){fPt[0]=lo; fPt[1]=hi;}
inline void AliFemtoQATrackCut::SetRapidity(const float& lo,const float& hi){fRapidity[0]=lo; fRapidity[1]=hi;}
inline void AliFemtoQATrackCut::SetCharge(const int& ch){fCharge = ch;}
inline void AliFemtoQATrackCut::SetPidProbElectron(const float& lo,const float& hi){fPidProbElectron[0]=lo; fPidProbElectron[1]=hi;}
inline void AliFemtoQATrackCut::SetPidProbPion(const float& lo,const float& hi){fPidProbPion[0]=lo; fPidProbPion[1]=hi;}
inline void AliFemtoQATrackCut::SetPidProbKaon(const float& lo,const float& hi){fPidProbKaon[0]=lo; fPidProbKaon[1]=hi;}
inline void AliFemtoQATrackCut::SetPidProbProton(const float& lo,const float& hi){fPidProbProton[0]=lo; fPidProbProton[1]=hi;}
inline void AliFemtoQATrackCut::SetPidProbMuon(const float& lo,const float& hi){fPidProbMuon[0]=lo; fPidProbMuon[1]=hi;}
inline void AliFemtoQATrackCut::SetLabel(const bool& flag){fLabel=flag;}
inline void AliFemtoQATrackCut::SetStatus(const long& status){fStatus=status;}
inline void AliFemtoQATrackCut::SetminTPCclsF(const short& minTPCclsF){fminTPCclsF=minTPCclsF;}
inline void AliFemtoQATrackCut::SetminTPCncls(const short& s){fminTPCncls=s;}
inline void AliFemtoQATrackCut::SetminITScls(const int& minITScls){fminITScls=minITScls;}
inline void AliFemtoQATrackCut::SetminTPCchiNdof(const float& s){fminTPCchiNdof = s;}
inline void AliFemtoQATrackCut::SetMostProbablePion() { fMostProbable = 2; }
inline void AliFemtoQATrackCut::SetMostProbableKaon() { fMostProbable = 3; }
inline void AliFemtoQATrackCut::SetMostProbableProton() { fMostProbable = 4; }
inline void AliFemtoQATrackCut::SetNoMostProbable() { fMostProbable = 0; }
inline void AliFemtoQATrackCut::SetMaxTPCncls(const short& s){fMaxTPCncls=s;}
inline void AliFemtoQATrackCut::SetMaxITSChiNdof(const float& maxchi) { fMaxITSchiNdof = maxchi; }
inline void AliFemtoQATrackCut::SetMaxTPCChiNdof(const float& maxchi) { fMaxTPCchiNdof = maxchi; }
inline void AliFemtoQATrackCut::SetMaxSigmaToVertex(const float& maxsig) { fMaxSigmaToVertex = maxsig; }

inline void AliFemtoQATrackCut::SetTPCnclsExclusionSwitch(const bool& excluSwitch) { fTPCnclsExclusionSwitch = excluSwitch; }
inline void AliFemtoQATrackCut::SetTPCnclsExclusion(const short& lo, const short& hi) {fTPCnclsExclusion[0] = lo; fTPCnclsExclusion[1] = hi;}
inline void AliFemtoQATrackCut::SetTPCchiNdofExclusionSwitch(const bool& excluSwitch) { fTPCchiNdofExclusionSwitch = excluSwitch; }
inline void AliFemtoQATrackCut::SetTPCchiNdofExclusion(const float& lo, const float& hi) {fTPCchiNdofExclusion[0] = lo; fTPCchiNdofExclusion[1] = hi;}

#endif

