#ifndef ALIMUONERRORCOUNTER_H
#define ALIMUONERRORCOUNTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONErrorCounter
/// \brief 
/// 
//  Author Alberto Baldisseri, JL Charvet 

#include <TNamed.h>

class AliMUONErrorCounter : public TNamed
{
  public :
    AliMUONErrorCounter(Int_t bp = 0, Int_t manu = 0, Int_t ev = 1);
    
    ///
    void Increment() {fEvents++;}
    ///
    Int_t BusPatch() const {return fBusPatch;}
    ///
    Int_t ManuId() const {return fManuId;}
    ///
    Int_t Events() const {return fEvents;}
    Int_t Compare(const TObject* obj) const;
    void Print(const Option_t* option="") const;
    void Print_uncal(const Option_t* option="") const;

  private :
    Int_t fBusPatch; ///< Buspath ID
    Int_t fManuId;   ///< Manu ID
    Int_t fEvents;   ///< Events with error in this buspatch
};

#endif //ALIMUONERRORCOUNTER_H 


