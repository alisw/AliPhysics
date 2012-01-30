///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoMCTrackCut: A basic track cut that used information from     //
// ALICE MC to accept or reject the track.                              //  
// Enables the selection on charge, transverse momentum, rapidity,       //
// and PDG of the particle										         //
// Authors: Malgorzata Janik (WUT)    majanik@cern.ch 					//
//			Lukasz Graczykowski (WUT) lgraczyk@cern.ch					 //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMCTRACKCUT_H
#define ALIFEMTOMCTRACKCUT_H


#include "AliFemtoTrackCut.h"

class AliFemtoMCTrackCut : public AliFemtoTrackCut 
{

 public:
  AliFemtoMCTrackCut();
  virtual ~AliFemtoMCTrackCut();

  virtual bool Pass(const AliFemtoTrack* aTrack);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtTrack;}

  void SetPt(const float& lo, const float& hi);
  void SetRapidity(const float& lo, const float& hi);
  void SetEta(const float& lo, const float& hi);
  void SetCharge(const int& ch);
  void SetPDG(const int& pdg);
  void SetLabel(const bool& flag);


 private:   // here are the quantities I want to cut on...

  int               fCharge;             // particle charge
  float             fPt[2];              // bounds for transverse momentum
  float             fRapidity[2];        // bounds for rapidity
  float             fEta[2];             // bounds for pseudorapidity
  bool              fLabel;              // if true label<0 will not pass throught 
  int               fPDGcode;            // PDG code of the particle
 
  long              fNTracksPassed;      // passed tracks count
  long              fNTracksFailed;      // failed tracks count



#ifdef __ROOT__ 
  ClassDef(AliFemtoMCTrackCut, 1)
#endif
    };


inline void AliFemtoMCTrackCut::SetPt(const float& lo, const float& hi){fPt[0]=lo; fPt[1]=hi;}
inline void AliFemtoMCTrackCut::SetRapidity(const float& lo,const float& hi){fRapidity[0]=lo; fRapidity[1]=hi;}
inline void AliFemtoMCTrackCut::SetEta(const float& lo,const float& hi){fEta[0]=lo; fEta[1]=hi;}
inline void AliFemtoMCTrackCut::SetCharge(const int& ch){fCharge = ch;}
inline void AliFemtoMCTrackCut::SetPDG(const int& pdg){fPDGcode = pdg;}
inline void AliFemtoMCTrackCut::SetLabel(const bool& flag){fLabel=flag;}


#endif

