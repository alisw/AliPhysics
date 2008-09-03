#ifndef AliEveTRDTrackList_H
#define AliEveTRDTrackList_H

#include <TEveElement.h>
#include <EveDet/AliEveTRDData.h>
#include <TList.h>


class AliEveTRDTrackList: public TEveElementList
{
  friend class AliEveTRDTrackListEditor;

public:
  AliEveTRDTrackList(const Text_t* n = "AliEveTRDTrackList", const Text_t* t = "", Bool_t doColor = kFALSE);
 
protected:
  TList* macroList;                 // List of (process) macros
  TList* macroSelList;              // List of (selection) macros

private:
  AliEveTRDTrackList(const AliEveTRDTrackList&);            // Not implemented
  AliEveTRDTrackList& operator=(const AliEveTRDTrackList&); // Not implemented

  ClassDef(AliEveTRDTrackList, 0);  // Class containing a list of tracks
};

#endif
