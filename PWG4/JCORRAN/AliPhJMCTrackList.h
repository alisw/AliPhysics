// $Id: AliPhJMCTrackList.h,v 1.4 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliPhJMCTrackList.hh
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.4 $
  \date $Date: 2008/05/08 13:44:45 $
*/
////////////////////////////////////////////////////

#ifndef ALIPHJMCTRACKLIST_H
#define ALIPHJMCTRACKLIST_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "TClonesArray.h"
#include <iostream>

#include "JConst.h"
#include "AliJMCTrack.h"

class AliJMCTrack;
class TClonesArray;

class AliPhJMCTrackList : public TObject {

public:
  AliPhJMCTrackList();
  AliPhJMCTrackList(expName exp);
  AliPhJMCTrackList(const AliPhJMCTrackList& a);  
  virtual ~AliPhJMCTrackList();

  void Reset();
  //getters
  unsigned short GetNTracks() const { return fTracks; }
  AliJMCTrack*        GetTrack(const unsigned int itrk); 
  //setters
  void SetNTracks(const unsigned short ntrk) { fTracks = ntrk; }
  int  SetTClonesArraySize(const unsigned int ntrk);
  // add tracks
  void AddJMCTrack(const unsigned int itrk);     // MC add

  AliPhJMCTrackList& operator=(const AliPhJMCTrackList&  list);

protected:
  TClonesArray *GetList() const { return fMcTrackList; }
  TClonesArray *fMcTrackList;    //list of MC tracks
  unsigned short fTracks;        //number of objects in the list

private:
  ClassDef(AliPhJMCTrackList,1)

};

#endif
