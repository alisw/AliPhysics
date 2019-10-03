/*
***********************************************************
  FMD data class
  Contact: jacobus.onderwaater@cern.ch
  2015/04/08
  *********************************************************
*/

#include "AliReducedFMDInfo.h"
#include <TClass.h>

ClassImp(AliReducedFMDInfo)


//_______________________________________________________________________________
AliReducedFMDInfo::AliReducedFMDInfo() :
  fMultiplicity(0.0),
  fId(0)
{
  AliReducedFMDInfo::Class()->IgnoreTObjectStreamer();
  //
  // Constructor
  //
}


//_______________________________________________________________________________
AliReducedFMDInfo::AliReducedFMDInfo(const AliReducedFMDInfo &c) :
  TObject(c),
  fMultiplicity(c.Multiplicity()),
  fId(c.Id())
{
  //
  // copy constructor
  //
}


//_______________________________________________________________________________
AliReducedFMDInfo::~AliReducedFMDInfo()
{
  //
  // De-Constructor
  //
}
