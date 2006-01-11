/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup management
/// \class AliMpDEIterator
/// \brief The iterator over valid detection element IDs
///
/// The valid detection element Ids are defined in the files denames*.dat.    \n
/// It can iterate 
/// - over all valid detection elements, if started with First() function; 
/// - or over detection elements in a selected module, if started with
///   First(Int_t moduleId) function                                          \n 
/// 
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MP_DE_ITERATOR_H
#define ALI_MP_DE_ITERATOR_H

#include "AliMpStationType.h"

#include <TObject.h>
#include <TArrayI.h>

class TString;

class AliMpDEIterator : public  TObject {

  public:
    AliMpDEIterator();
    AliMpDEIterator(const AliMpDEIterator& rhs);
    virtual ~AliMpDEIterator();

    // Operators
    AliMpDEIterator& operator=(const AliMpDEIterator& rhs);
    
    // Methods for iterating over DE elements
    // 
    void First();
    void First(Int_t moduleId);
    void Next();
    Bool_t IsDone() const;
    Int_t CurrentDE() const;

  private:
    // static methods
    static Bool_t ReadDEIds(AliMpStationType station);
    static void   ReadData();

    // static data members	
    static const Int_t  fgkMaxNofDetElements; // Maximum number of DEs
    static TArrayI      fgDetElemIds;         // DE Ids	
    static Int_t        fgNofDetElemIds;      // Number of DE Ids	

    // data members	
    Int_t  fIndex;    // Current DE index
    Int_t  fModuleId; // The iterated module 

  ClassDef(AliMpDEIterator,0)  // MUON Factory for Chambers and Segmentation
};

#endif //ALI_MP_DE_ITERATOR_H















