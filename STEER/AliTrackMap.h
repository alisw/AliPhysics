#ifndef ALITRACKMAP_H
#define ALITRACKMAP_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliTrackMap.h
// description:
//   contains a relation between track label and it's index
//   in a TreeH.
//   See http://AliSoft.cern.ch/people/chudoba/classes/AliTrackMap.html
//                  
//  Author: Jiri Chudoba (CERN), 2002
//
////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

typedef enum { kOutOfBounds = -2, kNoEntry} MapConsts_t;

class AliTrackMap: public TNamed {

public:
  AliTrackMap();
  AliTrackMap(Int_t size, Int_t *array);
  AliTrackMap(const AliTrackMap& trm);
  AliTrackMap& operator=(const AliTrackMap& trm)
    {trm.Copy(*this); return(*this);}
  ~AliTrackMap();
  Int_t At(Int_t label) const;
  Int_t Size() const {return fSize;}
  void SetEventNr(Int_t eventNr);
  void PrintValues() const;

private:
  void Copy(TObject &trm) const;

  Int_t fSize;             // size of the array
  Int_t *fArray;           //[fSize] actual map

  ClassDef(AliTrackMap,1)  // connection between track label and TreeH indeces
};

#endif // ALITRACKMAP_H





