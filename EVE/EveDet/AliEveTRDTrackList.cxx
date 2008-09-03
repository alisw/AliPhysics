#include "AliEveTRDTrackList.h"

ClassImp(AliEveTRDTrackList)

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTrackList ////////////////////////
///////////////////////////////////////////////////////////
AliEveTRDTrackList::AliEveTRDTrackList(const Text_t* n, const Text_t* t, Bool_t doColor):
  TEveElementList(n, t, doColor),
  macroList(0),
  macroSelList(0)
{
  SetChildClass(AliEveTRDTrack::Class());

  macroList = new TList();
  macroSelList = new TList();
}
