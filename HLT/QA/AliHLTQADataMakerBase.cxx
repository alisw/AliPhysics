// $Id$

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

/** @file   AliHLTQADataMakerBase.cxx
    @author Matthias Richter
    @date   2010-03-10
    @brief  
*/
#include "AliHLTQADataMakerBase.h"
#include "AliESDEvent.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTQADataMakerBase)

AliHLTQADataMakerBase::AliHLTQADataMakerBase()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTQADataMakerBase::~AliHLTQADataMakerBase()
{
  // see header file for class documentation
}

void AliHLTQADataMakerBase::Exec(AliQAv1::TASKINDEX_t task, TObject * data) 
{ 
  // special handling for esds
  if ( task == AliQAv1::kESDS ) {
    AliESDEvent * esd = NULL;
    AliESDEvent * hltesd = NULL;
    if (data->IsA() == AliESDEvent::Class()) {
      // silently skip this. Currently HLT QA is still called as
      // part of AliQAManager::RunOneEvent with the esd
      return;
    }
    if (data->InheritsFrom("TObjArray")) {
      TObjArray* array=dynamic_cast<TObjArray*>(data);
      if (array && array->GetEntriesFast()>0) {
	esd = dynamic_cast<AliESDEvent *>(array->At(0)) ;
      }
      if (array && array->GetEntriesFast()>1) {
	hltesd = dynamic_cast<AliESDEvent *>(array->At(1)) ;
      }
    } else {
      esd = static_cast<AliESDEvent *>(data) ; 
    }

    if (esd && strcmp(esd->ClassName(), "AliESDEvent") == 0) {
      if (hltesd) {
	MakeESDs(esd, hltesd);
      } else {
	AliError(Form("HLT ESD missing or wrong class type (%p), skipping HLT QA task kESDs", hltesd));
      }
    } else {
      AliError(Form("ESD missing or wrong class type (%p), skipping HLT QA task kESDSs", esd));
    }
  } else {
    // forward for all other types
    AliQADataMakerRec::Exec(task, data);
  }
}

void AliHLTQADataMakerBase::MakeESDs(AliESDEvent * esd)
{
  // see header file for class documentation
  
  // as an extension in the QA interface also the hlt esd can be sent
  // in order to preserve backward compatibility, a new function has been
  // introduced.
  //
  // NOTE: This function is not the place for HLT QA
  if (!esd) return;
}

void AliHLTQADataMakerBase::MakeESDs(AliESDEvent * esd, AliESDEvent* hltesd)
{
  // HLT QA on ESDs
  if (!esd || !hltesd) {
    AliError("invalid parameter: missing ESDs");
    return;
  }

  // nothing to do in the base class, QA implemented in the overloaded
  // child function
}
