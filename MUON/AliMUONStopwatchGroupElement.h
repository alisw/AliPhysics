#ifndef ALIMUONSTOPWATCHGROUPELEMENT_H
#define ALIMUONSTOPWATCHGROUPELEMENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONStopwatchGroupElement
/// \brief A class to group timers by name
/// 
// Author Laurent Aphecetche

#ifndef ALIMUONSTOPWATCHGROUP_H
#  include "AliMUONStopwatchGroup.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONStopwatchGroupElement
{
public:
  AliMUONStopwatchGroupElement(AliMUONStopwatchGroup* group, const char* a, const char *b)
  : fGroup(group), fA(a), fB(b)
  { group->Start(a,b); }
  AliMUONStopwatchGroupElement(const AliMUONStopwatchGroupElement& rhs) : fGroup(0),fA(),fB()
  { fGroup = rhs.fGroup; fA = rhs.fA; fB=rhs.fB ; }
  AliMUONStopwatchGroupElement& operator=(const AliMUONStopwatchGroupElement& rhs)
  { if ( this != &rhs ) { fGroup = rhs.fGroup; fA = rhs.fA; fB=rhs.fB ; } return *this; }
  
  ~AliMUONStopwatchGroupElement()
  { fGroup->Stop(fA.Data(),fB.Data()); }
  
private:
  AliMUONStopwatchGroup* fGroup; // the group for which we're just a proxy
  TString fA; // first parameter
  TString fB; // second parameter
};

#endif
