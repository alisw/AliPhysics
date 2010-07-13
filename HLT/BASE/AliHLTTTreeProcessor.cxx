// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timur Pocheptsov <Timur.Pocheptsov@cern.ch>           *
//*                  Matthias Richter <Matthias.Richter@cern.ch>
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
/// @author Timur Pocheptsov, Matthias Richter
/// @date   05.07.2010
/// @brief  Generic component for data collection in a TTree

#include <cerrno>
#include <memory>

#include "AliHLTTTreeProcessor.h"
#include "TDirectory.h"
#include "TDatime.h"
#include "TString.h"
#include "TTree.h"
#include "TH1.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTTreeProcessor)

AliHLTTTreeProcessor::AliHLTTTreeProcessor()
                        : AliHLTProcessor(), 
                          fDefinitions(),
                          fTree(0),
                          fMaxEntries(kMaxEntries),
                          fPublishInterval(kInterval),
                          fLastTime(0)
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

void AliHLTTTreeProcessor::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // get the output size estimator
  //
  if (!fDefinitions.size()) {
    HLTError("Can not calculate output data size, no histogram definitions were provided");
    return;
  }

  constBase = 0;
  for (list_const_iterator i = fDefinitions.begin(); i != fDefinitions.end(); ++i)
    constBase += i->GetSize();

  inputMultiplier = 1.;
}

int AliHLTTTreeProcessor::DoInit(int argc, const char** argv)
{
  // init component
  // ask child to create the tree.
  int iResult = 0;

  if (!fTree) {
    std::auto_ptr<TTree> ptr(CreateTree(argc, argv));
    if (ptr.get()) {
      //Stage 1: default initialization.
      ptr->SetDirectory(0);
      //"Default" (for derived component) histograms.
      FillHistogramDefinitions();
      //Default values.
      fMaxEntries = kMaxEntries;
      fPublishInterval = kInterval;
      fLastTime = 0;
      //Stage 2: OCDB.
      TString cdbPath("HLT/ConfigHLT/");
      cdbPath += GetComponentID();
      //
      iResult = ConfigureFromCDBTObjString(cdbPath);
      //
      if (iResult < 0)
        return iResult;
      //Stage 3: command line arguments.
      if (argc && (iResult = ConfigureFromArgumentString(argc, argv)) <= 0)
        return iResult;

      ptr->SetCircular(fMaxEntries);
      fTree = ptr.release();
    } else //No way to process error correctly - error is unknown here.
      return -EINVAL;
  } else {
    HLTError("fTree pointer must be null before DoInit call");
    return -EINVAL;
  }

  return iResult;
}

int AliHLTTTreeProcessor::DoDeinit()
{
  // cleanup component
  delete fTree;
  fTree = 0;
  fDefinitions.clear();
  return 0;
}

int AliHLTTTreeProcessor::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData)
{
  //Process event and publish histograms.
  //I'm pretty sure, that if fTree == 0 (DoInit failed) DoEvent is not called.
  //But interface itself does not force you to call DoInit before DoEvent, so,
  //I make this check explicit.
  if (!fTree) {
    HLTError("fTree is a null pointer, try to call AliHLTTTreeProcessor::DoInit first.");
    return -EINVAL;//-ENULLTREE? :)
  }

  // process input data blocks and fill the tree
  const int iResult = FillTree(fTree, evtData, trigData);

  if (iResult < 0)
    return iResult;

  const TDatime time;

  if (fLastTime - time.Get() > fPublishInterval) {
    for (list_const_iterator i = fDefinitions.begin(); i != fDefinitions.end(); ++i) {
      if (TH1* h = CreateHistogram(*i)) {
        //I do not care about errors here - since I'm not able
        //to rollback changes.
        PushBack(h, GetOriginDataType(), GetDataSpec());
      }
    }

    fLastTime = time.Get();
  }

  return iResult;
}

int AliHLTTTreeProcessor::ScanConfigurationArgument(int argc, const char** argv)
{
  // scan one argument and its parameters from the list
  // return number of processed entries.
  // possible arguments: 
  // -maxentries number
  // -interval number
  // -histogram name -size number -expression expression [-cut expression ][-opt option]
  // As soon as "-histogram" found, -size and -expression and -outtype are required, 
  // cut and option can be omitted.
  if (argc <= 0)
    return 0;

  std::list<AliHLTHistogramDefinition> newDefs;
  AliHLTHistogramDefinition def;

  int i = 0;
  int maxEntries = 0;

  while (i < argc) {
    const TString argument(argv[i]);

    if (argument.CompareTo("-maxentries") == 0) { //1. Max entries argument for TTree.
      if (i + 1 == argc) {
        HLTError("Numeric value for '-maxentries' is expected");
        return -EPROTO;
      }
      //Next must be a number.
      //TString returns 0 (number) even if string contains non-numeric symbols.
      maxEntries = TString(argv[i + 1]).Atoi();
      if (maxEntries <= 0) {
        HLTError("Bad value for '-maxentries': %d", maxEntries);
        return -EPROTO;
      }
  
      i += 2;
    } else if (argument.CompareTo("-interval") == 0) { //2. Interval argument for publishing.
      if (i + 1 == argc) {
        HLTError("Numeric value for '-interval' is expected");
        return -EPROTO;
      }

      const Int_t interval = TString(argv[i + 1]).Atoi();
      if (interval <= 0) {
        HLTError("Bad value for '-interval' argument: %d", interval);
        return -EPROTO;
      }

      fPublishInterval = interval;

      i += 2;
    } else if (argument.CompareTo("-histogram") == 0) { //3. Histogramm definition.
      const int nParsed = ParseHistogramDefinition(argc, argv, i, def);
      if (!nParsed)
        return -EPROTO;

      newDefs.push_back(def);

      i += nParsed;   
    } else {
      HLTError("Unknown argument %s", argument.Data());
      return -EPROTO;
    }
  }

  if (maxEntries != fMaxEntries) {
    fMaxEntries = maxEntries;
    if (fTree) {
      fTree->Reset();
      fTree->SetCircular(fMaxEntries);
    }
  }

  if (newDefs.size())
    fDefinitions.swap(newDefs);

  return i;
}

TH1* AliHLTTTreeProcessor::CreateHistogram(const AliHLTHistogramDefinition& d)
{
  // create a histogram from the tree
  if (!fTree) {
    HLTError("fTree is a null pointer, try to call AliHLTTTreeProcessor::DoInit first.");
    return 0;
  }

  const Long64_t rez = fTree->Project(d.GetName().Data(), d.GetExpression().Data(), d.GetCut().Data(), d.GetDrawOption().Data()); 

  if (rez == -1) {
    HLTError("TTree::Project failed");
    return 0;
  }

  return dynamic_cast<TH1*>(gDirectory->Get(d.GetName().Data()));
}

int AliHLTTTreeProcessor::ParseHistogramDefinition(int argc, const char** argv, int pos, AliHLTHistogramDefinition& dst)const
{
  //Histogram-definition:
  //    -histogram name -size number -expression expression [-cut expression][-opt option]

  //at pos we have '-histogram', at pos + 1 must be the name.
  if (pos + 1 == argc) {
    HLTError("Bad histogram definition, histogram name is expected");
    return 0;
  }

  dst.SetName(argv[pos + 1]);
  pos += 2;
  
  //At pos must be '-size', and number at pos + 1.
  if (pos == argc || TString(argv[pos]).CompareTo("-size")) {
    HLTError("Bad histogram definition, '-size' is expected");
    return 0;
  }

  if (pos + 1 == argc) {
    HLTError("Bad histogram definition, size is expected");
    return 0;
  }

  dst.SetSize(TString(argv[pos + 1]).Atoi());
  if (dst.GetSize() <= 0) {
    HLTError("Bad histogram definition, positive size is required");
    return 0;
  }

  pos += 2;
  //At pos must be '-expression', and expression at pos + 1. 
  if (pos == argc || TString(argv[pos]).CompareTo("-expression")) {
    HLTError("Bad histogram definition, '-expression' is expected");
    return 0;
  }

  if (pos + 1 == argc) {
    HLTError("Bad histogram definition, expression is expected");
    return 0;
  }

  dst.SetExpression(argv[pos + 1]);
  pos += 2;

  int processed = 6;
  dst.SetCut("");
  dst.SetDrawOption("");

  //remaining options can be the cut and Draw option.
  //cut must be first.
  if (pos + 1 >= argc)
    return processed;

  if (TString(argv[pos]).CompareTo("-cut") == 0) {
    dst.SetCut(argv[pos + 1]);
    pos += 2;
    processed += 2;
  }

  if (pos + 1 >= argc)
    return processed;

  if (TString(argv[pos]).CompareTo("-opt") == 0) {
    dst.SetDrawOption(argv[pos + 1]);
    processed += 2;
  }

  return processed;
}
