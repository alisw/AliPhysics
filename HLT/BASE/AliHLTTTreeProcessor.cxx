// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors:                                                       *
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

/// @file   AliHLTTTreeProcessor.cxx
/// @author 
/// @date   05.07.2010
/// @brief  Generic component for data collection in a TTree

#include "AliHLTTTreeProcessor.h"
#include "TString.h"
#include "TTree.h"
#include "TH1.h"
#include "TDatime.h"
#include "TDirectory.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTTreeProcessor)

AliHLTTTreeProcessor::AliHLTTTreeProcessor()
  : AliHLTProcessor()
  , fTree(NULL)
  , fMaxEntries(1000)
  , fPublishInterval(5)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTTreeProcessor::~AliHLTTTreeProcessor()
{
  // see header file for class documentation
}

AliHLTComponentDataType AliHLTTTreeProcessor::GetOutputDataType()
{
  // get the component output data type
  return kAliHLTDataTypeHistogram;
}

void AliHLTTTreeProcessor::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // get the output size estimator
  //
  // TODO: the output size most likely fixed and independent of the input data size.
  // calculate constBase from the list of histogram definitions
}

int AliHLTTTreeProcessor::DoInit( int argc, const char** argv )
{
  // init component

  // ask child to create the tree, passing of arguments needs to be discussed
  fTree=CreateTree(0, NULL);
  if (fTree) {
    // memory resident tree, TODO check if this is ok
    fTree->SetDirectory(NULL);
    fTree->SetCircular(fMaxEntries);
  }
}

int AliHLTTTreeProcessor::DoDeinit()
{
  // cleanup component
  delete fTree;
  fTree=NULL;

  // TODO: cleanup definition table
}

int AliHLTTTreeProcessor::DoEvent( const AliHLTComponentEventData& evtData,
				   AliHLTComponentTriggerData& trigData )
{
  // process input data blocks and fill the tree
  int iResult=0;

  TDatime time;
  static unsigned lasttime=0;
  if (lasttime-time.Get()>fPublishInterval) {
    // TODO: loop over histogram definitions, create and publish
    // check for correct cleanup of the created histograms

    // TODO: think about data origin for the histograms and how to define it
    // - automatically from origin of the input blocks -> working for
    //   detector histogramming components
    // - optional data type field in the histogram definition?

    // {
    //   TH1* histo=CreateHistogram(...);
    //   iResult=PushBack(...);
    //   if (iResult==-ENOSPC) {
    // 	// update size estimator
    // 	... += GetLastObjectSize();
    //   }
    // }

    lasttime=time.Get();
  }
}

int AliHLTTTreeProcessor::ScanConfigurationArgument(int argc, const char** argv)
{
  // scan one argument and its parameters from the list
  // return number of processed entries
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -maxentries
  if (argument.CompareTo("-maxentries")==0) {
    if (argc<2) return -EPROTO;
    TString parameter=argv[++i];
    fMaxEntries=parameter.Atoi();
    return ++i;
  }    

  // -interval
  if (argument.CompareTo("-interval")==0) {
    if (argc<2) return -EPROTO;
    TString parameter=argv[++i];
    fPublishInterval=parameter.Atoi();
    return ++i;
  }

  // TODO add scanning of histogram definitions

  return 0;
}

TH1* AliHLTTTreeProcessor::CreateHistogram(TTree* pTree,
					   const char* name,
					   const char* expression, 
					   const char* selection,
					   const char* option)
{
  // create a histogram from the tree
  // TODO: check if Project() works correctly if the tree has no
  // directory assigned
  if (!pTree || !name || !expression) return NULL;
  pTree->Project(name, expression, selection, option);
  TObject* histogram=gDirectory->Get(name);
  return dynamic_cast<TH1*>(histogram);
}
