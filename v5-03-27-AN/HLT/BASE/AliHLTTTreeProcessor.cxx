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
#include "AliHLTErrorGuard.h"
#include "TDirectory.h"
#include "TDatime.h"
#include "TString.h"
#include "TTree.h"
#include "TH1.h"
#include "TStopwatch.h"
#include "TUUID.h"
#include "TSystem.h"
#include "TRandom3.h"

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
                          fMaxMemory(700000),
                          fMaxEventTime(0),
                          fNofEventsForce(0),
                          fForcedEventsCount(0),
                          fSkippedEventsCount(0),
                          fNewEventsCount(0),
                          fUniqueId(0),
                          fIgnoreCycleTime(10),
                          fCycleTimeFactor(1.0)
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

  // component configuration
  //Stage 1: default initialization.
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

  // calculating a unique id from the hostname and process id
  // used for identifying output of multiple components
  TUUID guid = GenerateGUID();
  union
  {
    UChar_t buf[16];
    UInt_t bufAsInt[4];
  };
  guid.GetUUID(buf);
  fUniqueId = bufAsInt[0];
  
  if (!fTree) {
    // originally foreseen to pass the arguments to the function, however
    // this is not appropriate. Argument scan via overloaded function
    // ScanConfigurationArgument
    std::auto_ptr<TTree> ptr(CreateTree(0, NULL));
    if (ptr.get()) {
      ptr->SetDirectory(0);
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

  int fillingtime=0;
  int publishtime=0;
  bool bDoFilling=true;
  bool bDoPublishing=false;
  const int cycleResetInterval=1000;
  if (fpEventTimer && fpCycleTimer) {
    averageEventTime=AliHLTUInt32_t(fpEventTimer->RealTime()*fgkTimeScale)/(GetEventCount()+1);
    fillingtime=int(fpEventTimer->RealTime()*fgkTimeScale);
    publishtime=fillingtime;
    fpEventTimer->Start(kFALSE);
    fpCycleTimer->Stop();
    averageCycleTime=AliHLTUInt32_t(fpCycleTimer->RealTime()*fgkTimeScale)/((GetEventCount()%cycleResetInterval)+1);
    // adapt processing to 3/4 of the max time
    bDoFilling=4*averageEventTime<3*fMaxEventTime ||
      (averageEventTime<fCycleTimeFactor*averageCycleTime && fpCycleTimer->RealTime()>fIgnoreCycleTime);
    if (fNofEventsForce>0 && fForcedEventsCount<fNofEventsForce) {
      fForcedEventsCount++;
      bDoFilling=true;
    }
  }

  // FIXME: there is still an unclear increase in memory consumption, even if the number of entries
  // in the tree is restricted. Valgrind studies did not show an obvious memory leak. This is likely
  // to be caused by something deep in the Root TTree functionality and needs to be studied in detail.
  ProcInfo_t ProcInfo;
  gSystem->GetProcInfo(&ProcInfo);
  if (ProcInfo.fMemResident>fMaxMemory) bDoFilling=false;

  // process input data blocks and fill the tree
  int iResult = 0;
  if (eventType!=gkAliEventTypeEndOfRun) {
    if (bDoFilling) {iResult=FillTree(fTree, evtData, trigData); fNewEventsCount++;}
    else fSkippedEventsCount++;
  }
  if (fpEventTimer) {
    fpEventTimer->Stop();
    fillingtime=int(fpEventTimer->RealTime()*fgkTimeScale)-fillingtime;
    if (fillingtime<0) fillingtime=0;
    fpEventTimer->Start(kFALSE);
  }

  if (iResult < 0) {
    ALIHLTERRORGUARD(5, "FillTree failed with %d, first event %d", iResult, GetEventCount());
    return iResult;
  }

  const TDatime time;

  if (( time.Get() - fLastTime > fPublishInterval && fNewEventsCount>0) ||
      eventType==gkAliEventTypeEndOfRun) {
    if ((bDoPublishing=fLastTime>0)) { // publish earliest after the first interval but set the timer

    for (list_const_iterator i = fDefinitions.begin(); i != fDefinitions.end(); ++i) {
      if (TH1* h = CreateHistogram(*i)) {
        //I do not care about errors here - since I'm not able
        //to rollback changes.
	// TODO: in case of -ENOSPC et the size of the last object by calling
	// GetLastObjectSize() and accumulate the necessary output buffer size
        PushBack(h, GetOriginDataType(), GetDataSpec());
	delete h;
      }
    }
    unsigned eventcount=GetEventCount()+1;
    HLTBenchmark("publishing %d histograms, %d entries in tree, %d new events since last publishing, accumulated %d of %d events (%.1f%%)", fDefinitions.size(), fTree->GetEntriesFast(), fNewEventsCount, eventcount-fSkippedEventsCount, eventcount, eventcount>0?(100*float(eventcount-fSkippedEventsCount)/eventcount):0);
    fNewEventsCount=0;
    HLTBenchmark("current memory usage %d %d", ProcInfo.fMemResident, ProcInfo.fMemVirtual);
    }

    fLastTime=time.Get();
    if (fLastTime==0) {
      // choose a random offset at beginning to equalize traffic for multiple instances
      // of the component
      gRandom->SetSeed(fUniqueId);
      fLastTime-=gRandom->Integer(fPublishInterval);
    }
  }

  if (fpEventTimer) {
    fpEventTimer->Stop();
    publishtime=int(fpEventTimer->RealTime()*fgkTimeScale)-publishtime;
    if (publishtime>fillingtime) publishtime-=fillingtime;
    else publishtime=0;

    averageEventTime=AliHLTUInt32_t(fpEventTimer->RealTime()*fgkTimeScale)/(GetEventCount()+1);

    // info output once every 5 seconds
    static UInt_t lastTime=0;
    if (time.Get()-lastTime>5 ||
	eventType==gkAliEventTypeEndOfRun ||
	bDoPublishing) {
      lastTime=time.Get();
      unsigned eventcount=GetEventCount()+1;
      HLTBenchmark("filling time %d us, publishing time %d, average total processing time %d us, cycle time %d us, accumulated %d of %d events (%.1f%%)", fillingtime, publishtime, averageEventTime, averageCycleTime, eventcount-fSkippedEventsCount, eventcount, eventcount>0?(100*float(eventcount-fSkippedEventsCount)/eventcount):0);
    }
  }
  if (fpCycleTimer) {
    bool bReset=(GetEventCount()%cycleResetInterval)==0;
    fpCycleTimer->Start(bReset);
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
  // -histogram name -size number -expression expression [-title expression ] -cut expression ][-opt option]
  // As soon as "-histogram" found, -size and -expression and -outtype are required, 
  // cut and option can be omitted.
  if (argc <= 0)
    return 0;

  std::list<AliHLTHistogramDefinition> newDefs;
  AliHLTHistogramDefinition def;

  int i = 0;
  int maxEntries = fMaxEntries;

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
      if (interval < 0) {
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
      if (time < 0) {
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
      if (count < 0) {
        HLTError("Bad value for '-forced-events' argument: %d", count);
        return -EPROTO;
      }

      fNofEventsForce = count;
      fForcedEventsCount=0;

      i += 2;
    } else if (argument.CompareTo("-ignore-cycletime") == 0) { // ignore cycle time for n sec
      if (i + 1 == argc) {
        HLTError("Numeric value for '-ignore-cycletime' is expected");
        return -EPROTO;
      }

      const Int_t time = TString(argv[i + 1]).Atoi();
      if (time < 0) {
        HLTError("Bad value for '-ignore-cycletime' argument: %d", time);
        return -EPROTO;
      }

      fIgnoreCycleTime = time;
      i += 2;
    } else if (argument.CompareTo("-maxmemory") == 0) { // maximum of memory in kByte to be used by the component
      if (i + 1 == argc) {
        HLTError("Numeric value for '-maxmemory' is expected");
        return -EPROTO;
      }

      const Int_t mem = TString(argv[i + 1]).Atoi();
      if (mem < 0) {
        HLTError("Bad value for '-maxmemory' argument: %d", time);
        return -EPROTO;
      }

      fMaxMemory = mem;
      i += 2;
    } else if (argument.CompareTo("-cycletime-factor") == 0) { // weight factor for cycle time
      if (i + 1 == argc) {
        HLTError("Numeric value for '-cycletime-factor' is expected");
        return -EPROTO;
      }

      const Float_t factor = TString(argv[i + 1]).Atof();
      if (factor < 0) {
        HLTError("Bad value for '-cycletime-factor' argument: %f", factor);
        return -EPROTO;
      }

      fCycleTimeFactor = factor;
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

  //Reformatting the histogram name
  TString str2=d.GetCut().Data();
  str2.ReplaceAll("Track_", "");
  str2.ReplaceAll("&&", " ");
  str2 = histName+" "+str2;
  hist->SetTitle(str2);

  if(d.GetTitle().Length()){
  
    //removing underscore
    TString axis=d.GetTitle().Data();
    axis.ReplaceAll("_{T}", "underscore{T}");
    axis.ReplaceAll("_", " ");
    axis.ReplaceAll("underscore{T}", "_{T}");
  
    hist->SetXTitle(axis);
    hist->GetXaxis()->CenterTitle();
  }
  return hist;
}

int AliHLTTTreeProcessor::ParseHistogramDefinition(int argc, const char** argv, int pos, AliHLTHistogramDefinition& dst)const
{
  //Histogram-definition:
  //    -histogram name -size number -expression expression [-title expression][-cut expression][-opt option]

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
  dst.SetTitle("");
  dst.SetCut("");
  dst.SetDrawOption("");

  //remaining options can be the title, cut and Draw option.
  //title must be first
  if (pos + 1 >= argc){
    return processed;
  }
  if (TString(argv[pos]).CompareTo("-title") == 0) {
    dst.SetTitle(argv[pos + 1]);
    pos += 2;
    processed += 2;
  }

  //cut must be second.
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
