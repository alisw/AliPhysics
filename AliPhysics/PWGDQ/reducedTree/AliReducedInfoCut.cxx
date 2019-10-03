/*
***********************************************************
  Implementation of AliReducedInfoCut class.
  Contact: iarsene@cern.ch
  2015/09/10
  *********************************************************
*/

#ifndef ALIREDUCEDINFOCUT_H
#include "AliReducedInfoCut.h"
#endif

ClassImp(AliReducedInfoCut)

//____________________________________________________________________________
AliReducedInfoCut::AliReducedInfoCut() :
  TNamed()
{
  //
  // default constructor
  //
}

//____________________________________________________________________________
AliReducedInfoCut::AliReducedInfoCut(const Char_t* name, const Char_t* title) :
  TNamed(name, title)
{
  //
  // named constructor
  //
}

//____________________________________________________________________________
AliReducedInfoCut::~AliReducedInfoCut() {
  //
  // destructor
  //
}
