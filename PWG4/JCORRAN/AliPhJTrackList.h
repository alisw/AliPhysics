#ifndef ALIPHJTRACKLIST_H
#define ALIPHJTRACKLIST_H

////////////////////////////////////////////////////
/*!
  \file AliPhJTrackList.hh
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.4 $
  \date $Date: 2008/05/08 13:44:45 $
*/
////////////////////////////////////////////////////

// $Id: AliPhJTrackList.h,v 1.4 2008/05/08 13:44:45 djkim Exp $

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "TClonesArray.h"
#include <iostream>

#include "JConst.h"
#include "AliPhJBaseTrack.h"
#include "AliJTrack.h"

class AliJTrack;
class AliPhJBaseTrack;
class TClonesArray;

class AliPhJTrackList : public TObject {

public:
  AliPhJTrackList();
  AliPhJTrackList(expName exp);
  AliPhJTrackList(const AliPhJTrackList& a);
  virtual ~AliPhJTrackList();

  void Reset();
  //getters
  unsigned short GetNTracks() const { return fTracks; }
  AliPhJBaseTrack* GetTrack(const unsigned int itrk); 
  AliJTrack*       GetAliJTrack(const unsigned int itrk); // ALICE getter
  //setters
  void SetNTracks(const unsigned short ntrk) { fTracks = ntrk; }
  int  SetTClonesArraySize(const unsigned int ntrk);
  // add tracks
  void AddAliJTrack(const unsigned int itrk);     // ALICE add

  AliPhJTrackList& operator=(const AliPhJTrackList& list);

protected:
  TClonesArray *GetList() const { return fTrackList; }
  TClonesArray *fTrackList; // track list
  unsigned short fTracks; //number of tracks in the list

private:
  ClassDef(AliPhJTrackList,1)

};

#endif
