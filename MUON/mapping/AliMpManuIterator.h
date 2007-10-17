#ifndef ALIMPMANUITERATOR_H
#define ALIMPMANUITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup management
/// \class AliMpManuIterator
/// \brief Class to loop over all manus of MUON Tracker
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMpBusPatch;
class TExMapIter;

class AliMpManuIterator : public TObject
{
public:
  
  AliMpManuIterator();
  virtual ~AliMpManuIterator();
  
  Bool_t Next(Int_t& detElemId, Int_t& manuId);
  
  void Reset();
  
private:

    AliMpBusPatch* NextBusPatch() const;
    
private:

    TExMapIter* fIterator; ///< internal iterator
    AliMpBusPatch* fCurrentBusPatch; ///< current bus patch
    Int_t fCurrentManuIndex; ///< current manu index in current bus patch
    
  ClassDef(AliMpManuIterator,1) // Iterator on MUON tracker manus
};

#endif
