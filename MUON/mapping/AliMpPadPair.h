/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPadPair.h,v 1.9 2006/05/23 13:07:29 ivana Exp $

/// \ingroup basic
/// \class AliMpPadPair
/// \brief Wrap up for std::pair<AliMpPad, AliMpPad>
/// to avoid problems with CINT.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PAD_PAIR_H
#define ALI_MP_PAD_PAIR_H

#include <TObject.h>

#include "AliMpPad.h"

class AliMpPadPair : public TObject
{
  public:
    AliMpPadPair(const AliMpPad& pad1, const AliMpPad& pad2);      
    AliMpPadPair(const AliMpPadPair& pair);      
    AliMpPadPair();
    virtual ~AliMpPadPair();

    // operators    
    Bool_t operator == (const AliMpPadPair& right) const;
    Bool_t operator != (const AliMpPadPair& right) const;
    AliMpPadPair& operator = (const AliMpPadPair& right);

    // methods
    AliMpPad GetFirst() const;  
    AliMpPad GetSecond() const;  

  private:
    // data members
    AliMpPad  fPadFirst;  ///<  first pad
    AliMpPad  fPadSecond; ///<  second pad
    
    
  ClassDef(AliMpPadPair,1) //utility class for the motif type
};

// inline functions

inline AliMpPad AliMpPadPair::GetFirst() const  { return fPadFirst; } 
inline AliMpPad AliMpPadPair::GetSecond() const { return fPadSecond; } 


#endif //ALI_MP_PAD_PAIR_H
