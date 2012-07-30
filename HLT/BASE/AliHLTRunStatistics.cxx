//-*- Mode: C++ -*-
// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTRunStatistics.cxx
    @author Jochen Thaeder
    @date   
    @brief  Base class for run statistics, for all detectors
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTRunStatistics.h"
#include <iostream>
#include <cerrno>

using std::cout;

ClassImp(AliHLTRunStatistics)
    
AliHLTRunStatistics::AliHLTRunStatistics()
  : TNamed("HLT", "HLT Run Statistics")
  , fNEvents(0)
  , fMyObjects()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTRunStatistics::AliHLTRunStatistics (const AliHLTRunStatistics& src)
  : TNamed(src)
  , fNEvents(src.fNEvents)
  , fMyObjects()
{
  // see header file for class documentation
  for (int i=0; i<src.fMyObjects.GetEntriesFast(); i++) {
    fMyObjects.Add(src.fMyObjects.Clone());
  }
  fMyObjects.SetOwner(kTRUE);
}

AliHLTRunStatistics& AliHLTRunStatistics::operator= (const AliHLTRunStatistics& src)
{
  // see header file for class documentation
  this->~AliHLTRunStatistics();
  new (this) AliHLTRunStatistics(src);
  return *this;
}

AliHLTRunStatistics::~AliHLTRunStatistics()
{
  // see header file for class documentation
  fMyObjects.Delete();
}

void AliHLTRunStatistics::Print(Option_t* option) const
{
  // see header file for class documentation
  cout << "============ " << GetTitle() << " ============" << endl;
  cout << "\t" << GetNEvents() << " event(s)" << endl;
  for (int i=0; i<fMyObjects.GetEntriesFast(); i++) {
    fMyObjects.Print(option);
  }
  cout << "==============================================" << endl;
}

void AliHLTRunStatistics::Copy(TObject &object) const
{
  // copy this to the specified object

  AliHLTRunStatistics* pStatistics=dynamic_cast<AliHLTRunStatistics*>(&object);
  if (pStatistics) {
    // copy members if target is a AliHLTTriggerDecision
    *pStatistics=*this;
  }

  // copy the base class
  TObject::Copy(object);
}

TObject *AliHLTRunStatistics::Clone(const char */*newname*/) const
{
  // create a new clone, classname is ignored

  return new AliHLTRunStatistics(*this);
}

void  AliHLTRunStatistics::Clear(Option_t * /*option*/)
{
  // clear the content
  fNEvents=0;
  fMyObjects.Clear();
}

int AliHLTRunStatistics::Add(const TObject* pObject)
{
  // Add clone of object

  if (!pObject) return -EINVAL;
  fMyObjects.Add(pObject->Clone());
  return 0;
}
