///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoAODTrackCut: A basic track cut that used information from     //
// ALICE AOD to accept or reject the track.                              //  
// Enables the selection on charge, transverse momentum, rapidity,       //
// pid probabilities, number of ITS and TPC clusters                     //
// Author: Adam Kisiel (WUT, OSU), Adam.Kisiel@cern.ch                   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOAODTRACKCUT_H
#define ALIFEMTOAODTRACKCUT_H

//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoTrackCut.h"

class AliFemtoAODTrackCut : public AliFemtoTrackCut 
{

 public:
  AliFemtoAODTrackCut();
  virtual ~AliFemtoAODTrackCut();

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
  void SetMaxChiNdof(const float& maxchi);
  void SetMaxSigmaToVertex(const float& maxsig);
  void SetMostProbablePion();
  void SetMostProbableKaon();
  void SetMostProbableProton();
  void SetNoMostProbable(); 

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
  float             fMaxchiNdof;         // maximum allowed chi2/ndof for TPC clusters
  float             fMaxSigmaToVertex;   // maximum allowed sigma to primary vertex
  long              fNTracksPassed;      // passed tracks count
  long              fNTracksFailed;      // failed tracks count
  int               fMostProbable;       // this particle type is required to be most probable

  float PidFractionElectron(float mom) const;
  float PidFractionPion(float mom) const;
  float PidFractionKaon(float mom) const;
  float PidFractionProton(float mom) const;

#ifdef __ROOT__ 
  ClassDef(AliFemtoAODTrackCut, 1)
#endif
    };


inline void AliFemtoAODTrackCut::SetPt(const float& lo, const float& hi){fPt[0]=lo; fPt[1]=hi;}
inline void AliFemtoAODTrackCut::SetRapidity(const float& lo,const float& hi){fRapidity[0]=lo; fRapidity[1]=hi;}
inline void AliFemtoAODTrackCut::SetCharge(const int& ch){fCharge = ch;}
inline void AliFemtoAODTrackCut::SetPidProbElectron(const float& lo,const float& hi){fPidProbElectron[0]=lo; fPidProbElectron[1]=hi;}
inline void AliFemtoAODTrackCut::SetPidProbPion(const float& lo,const float& hi){fPidProbPion[0]=lo; fPidProbPion[1]=hi;}
inline void AliFemtoAODTrackCut::SetPidProbKaon(const float& lo,const float& hi){fPidProbKaon[0]=lo; fPidProbKaon[1]=hi;}
inline void AliFemtoAODTrackCut::SetPidProbProton(const float& lo,const float& hi){fPidProbProton[0]=lo; fPidProbProton[1]=hi;}
inline void AliFemtoAODTrackCut::SetPidProbMuon(const float& lo,const float& hi){fPidProbMuon[0]=lo; fPidProbMuon[1]=hi;}
inline void AliFemtoAODTrackCut::SetLabel(const bool& flag){fLabel=flag;}
inline void AliFemtoAODTrackCut::SetMostProbablePion() { fMostProbable = 2; }
inline void AliFemtoAODTrackCut::SetMostProbableKaon() { fMostProbable = 3; }
inline void AliFemtoAODTrackCut::SetMostProbableProton() { fMostProbable = 4; }
inline void AliFemtoAODTrackCut::SetNoMostProbable() { fMostProbable = 0; }
inline void AliFemtoAODTrackCut::SetMaxChiNdof(const float& maxchi) { fMaxchiNdof = maxchi; }
inline void AliFemtoAODTrackCut::SetMaxSigmaToVertex(const float& maxsig) { fMaxSigmaToVertex = maxsig; }

#endif

