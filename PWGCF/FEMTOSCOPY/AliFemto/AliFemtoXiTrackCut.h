///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoESDTrackCut: A basic track cut that used information from     //
// ALICE ESD to accept or reject the track.                              //  
// Enables the selection on charge, transverse momentum, rapidity,       //
// pid probabilities, number of ITS and TPC clusters                     //
// Author: Marek Chojnacki (WUT), mchojnacki@knf.pw.edu.pl               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOXITRACKCUT_H
#define ALIFEMTOXITRACKCUT_H

#include "AliFemtoV0TrackCut.h"

class AliFemtoXiTrackCut : public AliFemtoV0TrackCut 
{
  public:
  enum XiType {kXiMinus = 0, kXiPlus=1};
  typedef enum XiType AliFemtoXiType;


  AliFemtoXiTrackCut();
  virtual ~AliFemtoXiTrackCut();

  virtual bool Pass(const AliFemtoXi* aXi);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtXi;}
 
  
 private:   // here are the quantities I want to cut on...



#ifdef __ROOT__ 
  ClassDef(AliFemtoXiTrackCut, 1)
#endif

};


#endif

