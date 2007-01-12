#ifndef ALIMPSLATMOTIFMAP_H
#define ALIMPSLATMOTIFMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup motif
/// \class AliMpSlatMotifMap
/// \brief A container to keep track of allocated motifs and motifTypes for slats
/// (both St345 and trigger ones).
///
//  Author Laurent Aphecetche, Subatech

#ifndef ROOT_TMap
#  include "TMap.h"
#endif

class AliMpVMotif;
class AliMpMotifType;
class TString;

class AliMpSlatMotifMap : public TObject
{
public:
  AliMpSlatMotifMap();
  virtual ~AliMpSlatMotifMap();
  
  AliMpVMotif* FindMotif(const TString& id) const;
  AliMpMotifType* FindMotifType(const TString& id) const;
  
  Bool_t AddMotif(AliMpVMotif* motif, Bool_t warn=kTRUE);
  Bool_t AddMotifType(AliMpMotifType* motifType, Bool_t warn=kTRUE);
        
  void Print(Option_t* opt="") const;
  
  void Reset();
  
private:
  TMap fMotifs; //< collection of motifs
  TMap fMotifTypes; //< collection of motifTypes
  
  ClassDef(AliMpSlatMotifMap,1) // 
};

#endif
