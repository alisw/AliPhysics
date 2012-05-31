/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 
// $MpId: AliMpDEIterator.h,v 1.5 2006/05/24 13:58:16 ivana Exp $ 

/// \ingroup management
/// \class AliMpDEIterator
/// \brief The iterator over detection elements
///
/// It can iterate 
/// - over all detection elements, if started with First() function; 
/// - or over detection elements in a selected chamber, if started with
///   First(Int_t chamberId) function                                          \n 
/// 
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MP_DE_ITERATOR_H
#define ALI_MP_DE_ITERATOR_H

#include <TObject.h>

#include <TArrayI.h>

class AliMpDetElement;
class TIterator;
class TString;

class AliMpDEIterator : public  TObject {

  public:
    AliMpDEIterator();
    //AliMpDEIterator(const AliMpDEIterator& rhs);
    virtual ~AliMpDEIterator();

    // Operators
    //AliMpDEIterator& operator=(const AliMpDEIterator& rhs);
    
    // Methods for iterating over DE elements
    // 
    void First();
    void First(Int_t chamberId);
    void Next();
    Bool_t IsDone() const;
    
    AliMpDetElement* CurrentDE() const;
    Int_t CurrentDEId() const;

  private:
    /// Not implemented
    AliMpDEIterator(const AliMpDEIterator& rhs);
    /// Not implemented
    AliMpDEIterator& operator=(const AliMpDEIterator& rhs);

    // data members	
    AliMpDetElement* fCurrentDE; ///< current element in iteration
    TIterator*       fIterator;  ///< iterator
    Int_t            fChamberId; ///< The iterated chamber 

  ClassDef(AliMpDEIterator,0)  // The iterator over valid detection element IDs
};

#endif //ALI_MP_DE_ITERATOR_H















