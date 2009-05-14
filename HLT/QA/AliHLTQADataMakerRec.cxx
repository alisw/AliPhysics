// $Id:$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTQADataMakerRec.cxx
    @author Matthias Richter
    @date   2009-05-14
    @brief  Container for the HLT offline QA
*/
#include "AliHLTQADataMakerRec.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTQADataMakerRec)

AliHLTQADataMakerRec::AliHLTQADataMakerRec()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTQADataMakerRec::~AliHLTQADataMakerRec()
{
  // see header file for class documentation
}

void AliHLTQADataMakerRec::StartOfDetectorCycle()
{
  // see header file for class documentation
}

void AliHLTQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray** /*list*/)
{
  // see header file for class documentation
}
