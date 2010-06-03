// @(#) $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTLoaderPublisherComponent.cxx
    @author Matthias Richter
    @date   
    @brief  A general tree publisher component for the AliLoader.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTLoaderPublisherComponent.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "TTree.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTLoaderPublisherComponent)

AliHLTLoaderPublisherComponent::AliHLTLoaderPublisherComponent()
  :
  fMaxSize(0),
  fLoaderType(),
  fTreeType("digits"),
  fVerbose(kFALSE),
  fDataType(kAliHLTAnyDataType),
  fSpecification(kAliHLTVoidDataSpec),
  fpLoader(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTLoaderPublisherComponent::~AliHLTLoaderPublisherComponent()
{
  // see header file for class documentation
}

const char* AliHLTLoaderPublisherComponent::GetComponentID()
{
  // see header file for class documentation
  return "AliLoaderPublisher";
}

AliHLTComponentDataType AliHLTLoaderPublisherComponent::GetOutputDataType()
{
  // see header file for class documentation
  return fDataType;
}

void AliHLTLoaderPublisherComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=fMaxSize;
  inputMultiplier=1;
}

AliHLTComponent* AliHLTLoaderPublisherComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTLoaderPublisherComponent;
}

int AliHLTLoaderPublisherComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  // scan arguments
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -loader
    if (argument.CompareTo("-loader")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fLoaderType=argv[i];

      // -tree
    } else if (argument.CompareTo("-tree")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fTreeType=argv[i];

      // -verbose
    } else if (argument.CompareTo("-verbose")==0) {
      fVerbose=kTRUE;

      // -datatype
    } else if (argument.CompareTo("-datatype")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&fDataType.fID, argv[i], TMath::Min(kAliHLTComponentDataTypefIDsize, (Int_t)strlen(argv[i])));
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&fDataType.fOrigin, argv[i], TMath::Min(kAliHLTComponentDataTypefOriginSize, (Int_t)strlen(argv[i])));

      // -dataspec
    } else if (argument.CompareTo("-dataspec")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fSpecification=(AliHLTUInt32_t)parameter.Atoi();
      } else if (parameter.BeginsWith("0x") &&
		 parameter.Replace(0,2,"",0).IsHex()) {
	sscanf(parameter.Data(),"%x", &fSpecification);
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
    } else {
      HLTError("unknown argument %s", argument.Data());
      iResult=-EINVAL;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (iResult<0) return iResult;

  if (fLoaderType.IsNull()) {
    AliErrorStream() << "loader type required, use \'-loader\' option" << endl;
    return -EINVAL;
  }

  // fetch runLoader instance from interface
  AliRunLoader* pRunLoader=GetRunLoader();
  if (pRunLoader) {

    // get the specific loader for the module
    fpLoader=pRunLoader->GetLoader(fLoaderType.Data());
    if (fpLoader) {
      // prepare the loader
      fpLoader->LoadDigits("read");

      // scan trough all events and estimate the size of the digits
      for (int i=0; i<pRunLoader->GetNumberOfEvents(); i++) {
	pRunLoader->GetEvent(i);
	TTree* pTree=GetTree();
	if (pTree) {
	  int size=EstimateObjectSize(pTree);
	  if (size>fMaxSize) fMaxSize=size;
	  if (fVerbose) {
	    AliInfoStream() << "event " << i << " " 
			    << fTreeType <<" size " << size 
			    << " count " << pTree->GetEntries() << endl;
	  }
	} else {
	  AliWarningStream() << "no " << fTreeType << " tree for event " << i << endl;
 	}
      }
    } else {
      AliErrorStream() << "can not get loader of type " << fLoaderType << endl;
      iResult=-EFAULT;
    }

    // prepare the RunLoader to provide the kinematics tree, some components
    // (e.g. the offline ITS Clusterfinder) need the kinematics tree to
    // propagate the mc information
    pRunLoader->LoadKinematics("READ");
  } else {
    AliErrorStream() << "can not get runLoader" << endl;
    iResult=-EFAULT;
  }
  return iResult;
}

int AliHLTLoaderPublisherComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  if (fpLoader) {
    fpLoader->UnloadDigits();
  }
  fpLoader=NULL;
  return iResult;
}

int AliHLTLoaderPublisherComponent::GetEvent(const AliHLTComponentEventData& /*evtData*/,
					     AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation

  // process data events only
  if (!IsDataEvent()) return 0;

  int iResult=0;
  // fetch runLoader instance from interface
  AliRunLoader* pRunLoader=GetRunLoader();
  if (pRunLoader && fpLoader) {
    pRunLoader->GetEvent(GetEventCount());
    TTree* pTree=GetTree();
    if (pTree) {
      PushBack(pTree, fDataType);
    } else {
      AliWarningStream() << "no " << fTreeType << " tree for event " << GetEventCount() << endl;
    }
  } else {
    AliErrorStream() << "component not initialized" << endl;
    iResult=-EFAULT;
  }
  return iResult;
}

TTree* AliHLTLoaderPublisherComponent::GetTree()
{
  // see header file for class documentation
  TTree* pTree=NULL;
  if (fpLoader) {
    if (fTreeType.CompareTo("digits")==0)
      pTree=fpLoader->TreeD();
    else if (fTreeType.CompareTo("clusters")==0) {
      pTree=fpLoader->TreeR();
    }
  } else {
  }
  return pTree;
}
