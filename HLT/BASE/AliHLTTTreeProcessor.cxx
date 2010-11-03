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
#include "TStopwatch.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTTreeProcessor)

AliHLTTTreeProcessor::AliHLTTTreeProcessor()
                        : AliHLTProcessor(), 
                          fDefinitions(),
                          fTree(0),
                          fMaxEntries(kMaxEntries),
                          fPublishInterval(kInterval),
                          fLastTime(0),
                          fpEventTimer(NULL),
                          fpCycleTimer(NULL),
                          fMaxEventTime(0),
                          fNofEventsForce(0),
                          fForcedEventsCount(0),
                          fSkippedEventsCount(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

const AliHLTUInt32_t AliHLTTTreeProcessor::fgkTimeScale=1000000; // ticks per second

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
      if (argc && (iResult = ConfigureFromArgumentString(argc, argv)) < 0)
        return iResult;

      ptr->SetCircular(fMaxEntries);
      fTree = ptr.release();
    } else //No way to process error correctly - error is unknown here.
      return -EINVAL;
  } else {
    HLTError("fTree pointer must be null before DoInit call");
    return -EINVAL;
  }

  if (iResult>=0 && fMaxEventTime>0) {
    fpEventTimer=new TStopwatch;
    if (fpEventTimer) {
      fpEventTimer->Reset();
    }
    fpCycleTimer=new TStopwatch;
    if (fpCycleTimer) {
      fpCycleTimer->Reset();
    }
  }
  fSkippedEventsCount=0;

  return iResult;
}

int AliHLTTTreeProcessor::DoDeinit()
{
  // cleanup component
  delete fTree;
  fTree = 0;
  fDefinitions.clear();

  if (fpEventTimer) delete fpEventTimer;
  fpEventTimer=NULL;
  if (fpCycleTimer) delete fpCycleTimer;
  fpCycleTimer=NULL;

  return 0;
}

int AliHLTTTreeProcessor::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData)
{
  //Process event and publish histograms.
  AliHLTUInt32_t eventType=0;
  if (!IsDataEvent(&eventType) && eventType!=gkAliEventTypeEndOfRun) return 0;

  //I'm pretty sure, that if fTree == 0 (DoInit failed) DoEvent is not called.
  //But interface itself does not force you to call DoInit before DoEvent, so,
  //I make this check explicit.
  if (!fTree) {
    HLTError("fTree is a null pointer, try to call AliHLTTTreeProcessor::DoInit first.");
    return -EINVAL;//-ENULLTREE? :)
  }

  AliHLTUInt32_t averageEventTime=0;
  AliHLTUInt32_t averageCycleTime=0;

  AliHLTUInt32_t proctime=0;
  bool bDoFilling=true;
  if (fpEventTimer && fpCycleTimer) {
    averageEventTime=(fpEventTimer->RealTime()*fgkTimeScale)/(GetEventCount()+1);
    proctime=fpEventTimer->RealTime()*fgkTimeScale;
    fpEventTimer->Start(kFALSE);
    fpCycleTimer->Stop();
    averageCycleTime=(fpCycleTimer->RealTime()*fgkTimeScale)/(GetEventCount()+1);
    // adapt processing to 3/4 of the max time
    bDoFilling=4*averageEventTime<3*fMaxEventTime || averageEventTime<averageCycleTime;
    if (fNofEventsForce>0 && fForcedEventsCount<fNofEventsForce) {
      fForcedEventsCount++;
      bDoFilling=true;
    }
  }

  // process input data blocks and fill the tree
  int iResult = 0;
  if (eventType!=gkAliEventTypeEndOfRun) {
    if (bDoFilling) iResult=FillTree(fTree, evtData, trigData);
    else fSkippedEventsCount++;
  }

  if (iResult < 0)
    return iResult;

  const TDatime time;

  if ( time.Get() - fLastTime > fPublishInterval ||
      eventType==gkAliEventTypeEndOfRun) {
    for (list_const_iterator i = fDefinitions.begin(); i != fDefinitions.end(); ++i) {
      if (TH1* h = CreateHistogram(*i)) {
        //I do not care about errors here - since I'm not able
        //to rollback changes.
	// TODO: in case of -ENOSPC et the size of the last object by calling
	// GetLastObjectSize() and accumulate the necessary output buffer size
        PushBack(h, GetOriginDataType(), GetDataSpec());
      }
    }

    fLastTime = time.Get();
  }

  if (fpEventTimer) {
    fpEventTimer->Stop();
    proctime=fpEventTimer->RealTime()*fgkTimeScale-proctime;
    averageEventTime=(fpEventTimer->RealTime()*fgkTimeScale)/(GetEventCount()+1);

    // info output once every 5 seconds
    static UInt_t lastTime=0;
    if (time.Get()-lastTime>5 ||
      eventType==gkAliEventTypeEndOfRun) {
      lastTime=time.Get();
      unsigned eventcount=GetEventCount();
      HLTBenchmark("event time %d us, average time %d us, cycle time %d us, accumulated %d of %d events (%.1f%%)", proctime, averageEventTime, averageCycleTime, eventcount-fSkippedEventsCount, eventcount, eventcount>0?(100*float(eventcount-fSkippedEventsCount)/eventcount):0);
    }
  }
  if (fpCycleTimer) {
    fpCycleTimer->Start(kFALSE);
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
    } else if (argument.CompareTo("-maxeventtime") == 0) { // max average processing time in us
      if (i + 1 == argc) {
        HLTError("Numeric value for '-maxeventtime' is expected");
        return -EPROTO;
      }

      const Int_t time = TString(argv[i + 1]).Atoi();
      if (time <= 0) {
        HLTError("Bad value for '-maxeventtime' argument: %d", time);
        return -EPROTO;
      }

      fMaxEventTime = time;

      i += 2;
    } else if (argument.CompareTo("-forced-events") == 0) { // number of forced events
      if (i + 1 == argc) {
        HLTError("Numeric value for '-forced-events' is expected");
        return -EPROTO;
      }

      const Int_t count = TString(argv[i + 1]).Atoi();
      if (count <= 0) {
        HLTError("Bad value for '-forced-events' argument: %d", count);
        return -EPROTO;
      }

      fNofEventsForce = count;
      fForcedEventsCount=0;

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

  TString histName(d.GetName());
  if (!histName.Contains("(")) {
    //Without number of bins, the histogram will be "fixed"
    //and most of values can go to underflow/overflow bins,
    //since kCanRebin will be false.
    histName += TString::Format("(%d)", Int_t(kDefaultNBins));
  }

  const Long64_t rez = fTree->Project(histName.Data(), d.GetExpression().Data(), d.GetCut().Data(), d.GetDrawOption().Data());

  if (rez == -1) {
    HLTError("TTree::Project failed");
    return 0;
  }

  //Now, cut off the binning part of a name
  histName = histName(0, histName.Index("("));
  TH1 * hist = dynamic_cast<TH1*>(gDirectory->Get(histName.Data()));
  if (!hist) {
    const TString msg(Form("Hist %s is a null pointer, selection was %s, strange name or hist's type\n", histName.Data(), d.GetExpression().Data()));
    HLTError(msg.Data());
  }else if (d.GetDrawOption().Length()) {
    hist->SetOption(d.GetDrawOption().Data());
  }

  return hist;
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
