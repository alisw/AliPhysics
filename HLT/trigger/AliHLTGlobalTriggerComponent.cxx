// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
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

/// @file   AliHLTGlobalTriggerComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   26 Nov 2008
/// @brief  Implementation of the AliHLTGlobalTriggerComponent component class.
///
/// The AliHLTGlobalTriggerComponent class applies the global HLT trigger to all
/// trigger information produced by components deriving from AliHLTTrigger.

#include "AliHLTGlobalTriggerComponent.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTGlobalTrigger.h"
#include "AliHLTGlobalTriggerConfig.h"
#include "AliHLTTriggerMenu.h"
#include "AliHLTCTPData.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliRawDataHeader.h"
#include "TUUID.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TInterpreter.h"
#include "TDatime.h"
#include "TClass.h"
#include "TNamed.h"
#include "TRandom3.h"
#include <fstream>
#include <cerrno>
#include <cassert>
#include <vector>
#include <algorithm>

ClassImp(AliHLTGlobalTriggerComponent)

const char* AliHLTGlobalTriggerComponent::fgkTriggerMenuCDBPath = "HLT/ConfigHLT/HLTGlobalTrigger";


namespace
{
  /**
   * This method is used as a comparison functor with the STL sort routines.
   */
  bool AliHLTDescendingNumbers(UInt_t a, UInt_t b)
  {
    return a > b;
  }
} // end of namespace


AliHLTGlobalTriggerComponent::AliHLTGlobalTriggerComponent() :
	AliHLTTrigger(),
	fTrigger(NULL),
	fDebugMode(false),
	fRuntimeCompile(true),
	fDeleteCodeFile(false),
	fMakeSoftwareTriggers(true),
	fCodeFileName(),
	fClassName(),
	fCTPDecisions(NULL),
	fBufferSizeConst(2*(sizeof(AliHLTGlobalTriggerDecision) + sizeof(AliHLTReadoutList))),
	fBufferSizeMultiplier(1.),
	fIncludePaths(TObjString::Class()),
	fIncludeFiles(TObjString::Class()),
	fLibStateAtLoad(),
	fBits(0),
	fDataEventsOnly(true),
	fMonitorPeriod(-1),
	fUniqueID(0),
	fSoftwareTrigger(true, "SOFTWARE"),
	fTotalEventCounter(0),
	fCDH(NULL)
{
  // Default constructor.
  
  ClearInfoForNewEvent(false);
}


AliHLTGlobalTriggerComponent::~AliHLTGlobalTriggerComponent()
{
  // Default destructor.
  
  if (fTrigger != NULL) delete fTrigger;

  if (fCTPDecisions) {
    fCTPDecisions->Delete();
    delete fCTPDecisions;
  }
}


void AliHLTGlobalTriggerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list) const
{
  // Returns the kAliHLTDataTypeGlobalTrigger type as output.
  list.push_back(kAliHLTDataTypeGlobalTrigger);
}


void AliHLTGlobalTriggerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // Returns the output data size estimate.

  constBase = fBufferSizeConst;
  inputMultiplier = fBufferSizeMultiplier;
}


Int_t AliHLTGlobalTriggerComponent::DoInit(int argc, const char** argv)
{
  // Initialises the global trigger component.
  
  fDebugMode = false;
  fClassName = "";
  fCodeFileName = "";
  fDeleteCodeFile = false;
  fMakeSoftwareTriggers = true;
  const char* configFileName = NULL;
  const char* codeFileName = NULL;
  fIncludePaths.Clear();
  fIncludeFiles.Clear();
  SetBit(kIncludeInput);
  fDataEventsOnly = true;
  UInt_t randomSeed = 0;
  bool randomSeedSet = false;
  fCDH = NULL;
  
  for (int i = 0; i < argc; i++)
  {
    if (strcmp(argv[i], "-config") == 0)
    {
      if (configFileName != NULL)
      {
        HLTWarning("Trigger configuration macro was already specified."
                   " Will replace previous value given by -config."
        );
      }
      if (argc <= i+1)
      {
        HLTError("The trigger configuration macro filename was not specified." );
        return -EINVAL;
      }
      configFileName = argv[i+1];
      i++;
      continue;
    }
    
    if (strcmp(argv[i], "-includepath") == 0)
    {
      if (argc <= i+1)
      {
        HLTError("The include path was not specified." );
        return -EINVAL;
      }
      try
      {
        new (fIncludePaths[fIncludePaths.GetEntriesFast()]) TObjString(argv[i+1]);
      }
      catch (const std::bad_alloc&)
      {
        HLTError("Could not allocate more memory for the fIncludePaths array.");
        return -ENOMEM;
      }
      i++;
      continue;
    }
    
    if (strcmp(argv[i], "-include") == 0)
    {
      if (argc <= i+1)
      {
        HLTError("The include file name was not specified." );
        return -EINVAL;
      }
      try
      {
        new (fIncludeFiles[fIncludeFiles.GetEntriesFast()]) TObjString(argv[i+1]);
      }
      catch (const std::bad_alloc&)
      {
        HLTError("Could not allocate more memory for the fIncludeFiles array.");
        return -ENOMEM;
      }
      i++;
      continue;
    }
    
    if (strcmp(argv[i], "-debug") == 0)
    {
      if (fDebugMode == true)
      {
        HLTWarning("The debug flag was already specified. Ignoring this instance.");
      }
      fDebugMode = true;
      continue;
    }
    
    if (strcmp(argv[i], "-cint") == 0)
    {
      fRuntimeCompile = false;
      continue;
    }
    
    if (strcmp(argv[i], "-usecode") == 0)
    {
      if (codeFileName != NULL)
      {
        HLTWarning("Custom trigger code file was already specified."
                   " Will replace previous value given by -usecode."
        );
      }
      if (argc <= i+1)
      {
        HLTError("The custom trigger code filename was not specified." );
        return -EINVAL;
      }
      codeFileName = argv[i+1];
      if (argc <= i+2)
      {
        HLTError("The custom trigger class name was not specified." );
        return -EINVAL;
      }
      fClassName = argv[i+2];
      i += 2;
      continue;
    }

    if (strcmp(argv[i], "-skipctp") == 0)
    {
      HLTInfo("Skipping CTP counters in trigger decision");
      SetBit(kSkipCTP);
      continue;
    }

    if (strcmp(argv[i], "-forward-input") == 0)
    {
      HLTInfo("Forwarding input objects and trigger decisions");
      SetBit(kForwardInput);
      SetBit(kIncludeShort);
      SetBit(kIncludeInput, false);
      continue;
    }

    if (strstr(argv[i], "-include-input") == argv[i])
    {
      SetBit(kForwardInput,false);
      TString param=argv[i];
      param.ReplaceAll("-include-input", "");
      if (param.CompareTo("=none")==0) 
      {
        HLTInfo("skipping objects and trigger decisions");
        SetBit(kIncludeShort, false);
        SetBit(kIncludeInput, false);
      }
      else if (param.CompareTo("=short")==0) 
      {
        HLTInfo("including short info on objects and trigger decisions");
        SetBit(kIncludeShort);
        SetBit(kIncludeInput, false);
      }
      else if (param.CompareTo("=both")==0) 
      {
        HLTInfo("including input objects, trigger decisions and short info");
        SetBit(kIncludeShort);
        SetBit(kIncludeInput);
      }
      else if (param.CompareTo("=objects")==0 || param.IsNull())
      {
        HLTInfo("including input objects and trigger decisions");
        SetBit(kIncludeShort, false);
        SetBit(kIncludeInput);
      }
      else
      {
        HLTError("unknown parameter '%s' for argument '-include-input'", param.Data());
      }
      continue;
    }

    if (strcmp(argv[i], "-process-all-events") == 0)
    {
      fDataEventsOnly = false;
      continue;
    }

    if (strstr(argv[i], "-monitoring") == argv[i])
    {
      TString param=argv[i];
      param.ReplaceAll("-monitoring", "");
      if (param.IsNull()) 
      {
	// -monitoring
	// enable monitoring trigger for all events
	fMonitorPeriod=0;
      } else {
	// -monitoring=n
	// enable monitoring trigger once every n seconds
	param.ReplaceAll("=", "");
	if (param.IsDigit()) {
	  fMonitorPeriod=param.Atoi();
	} else {
	  HLTError("expecting number as parameter for argument '-monitoring=', got %s, skipping monitoring trigger", param.Data());
	}
      }
      continue;
    }

    if (strcmp(argv[i], "-dont-make-software-triggers") == 0)
    {
      fMakeSoftwareTriggers = false;
      continue;
    }
    if (strcmp(argv[i], "-randomseed") == 0)
    {
      if (randomSeedSet)
      {
        HLTWarning("The random seed was already specified previously with the option -randomseed."
                   "Will replace the previous value of %d with the new one.",
                   randomSeed
        );
      }
      if (argc <= i+1)
      {
        HLTError("The number to use as the seed was not specified for the -randomseed option." );
        return -EINVAL;
      }
      TString numstr = argv[i+1];
      if (not numstr.IsDigit())
      {
        HLTError("The number specified in the -randomseed option is not a valid decimal integer." );
        return -EINVAL;
      }
      randomSeed = numstr.Atoi();
      randomSeedSet = true;
      i++;
      continue;
    }
    
    HLTError("Unknown option '%s'.", argv[i]);
    return -EINVAL;
  } // for loop
  
  const AliHLTTriggerMenu* menu = NULL;
  if (configFileName != NULL)
  {
    TString cmd = ".x ";
    cmd += configFileName;
    gROOT->ProcessLine(cmd);
    menu = AliHLTGlobalTriggerConfig::Menu();
  }
  
  // Try load the trigger menu from the CDB if it is not loaded yet with the
  // -config option
  int result = -ENOENT;
  if (menu == NULL)
  {
    result = LoadTriggerMenu(fgkTriggerMenuCDBPath, menu);
  }
  if (menu == NULL)
  {
    HLTError("No trigger menu configuration found or specified.");
    return result;
  }
  
  if (codeFileName == NULL)
  {
    result = GenerateTrigger(menu, fClassName, fCodeFileName, fIncludePaths, fIncludeFiles);
    if (result == 0) fDeleteCodeFile = true;
  }
  else
  {
    result = LoadTriggerClass(codeFileName, fIncludePaths);
    if (result == 0) fCodeFileName = codeFileName;
  }
  if (result != 0) return result;
  
  try
  {
    fTrigger = AliHLTGlobalTrigger::CreateNew(fClassName.Data());
  }
  catch (const std::bad_alloc&)
  {
    HLTError("Could not allocate memory for the AliHLTGlobalTrigger instance.");
    return -ENOMEM;
  }
  if (fTrigger == NULL)
  {
    HLTError("Could not create a new instance of '%s'.", fClassName.Data());
    return -EIO;
  }
  
  fTrigger->FillFromMenu(*menu);
  if (fTrigger->CallFailed()) return -EPROTO;

  // setup the CTP accounting in AliHLTComponent
  SetupCTPData();

  // Set the default values from the trigger menu.
  SetDescription(menu->DefaultDescription());
  SetTriggerDomain(menu->DefaultTriggerDomain());
  
  // Initialise the random number generator seed value.
  // NOTE: The GenerateTrigger method called above will set the fUniqueID value
  // with a random value based on a GUID that should be unique across the system.
  // This is then used as the seed to the random number generator if the -randomseed
  // option is not used.
  if (not randomSeedSet) randomSeed = fUniqueID;
  gRandom->SetSeed(randomSeed);
  
  fTotalEventCounter = 0;
  return 0;
}


Int_t AliHLTGlobalTriggerComponent::DoDeinit()
{
  // Cleans up the global trigger component.
  
  if (fTrigger != NULL)
  {
    delete fTrigger;
    fTrigger = NULL;
  }
  
  if (fCTPDecisions) {
    fCTPDecisions->Delete();
    delete fCTPDecisions;
  }
  fCTPDecisions=NULL;
  fCDH = NULL;
  
  Int_t result = UnloadTriggerClass(fCodeFileName);
  if (result != 0) return result;
  
  if (fDeleteCodeFile and !fCodeFileName.IsNull() && gSystem->AccessPathName(fCodeFileName)==0 && !fDebugMode) {
    fCodeFileName.ReplaceAll(".cxx", "*");
    TString command="rm "; command+=fCodeFileName;
    gSystem->Exec(command);
  }
  fCodeFileName="";
  fDeleteCodeFile=false;

  return 0;
}


AliHLTComponent* AliHLTGlobalTriggerComponent::Spawn()
{
  // Creates a new object instance.
  AliHLTComponent* comp = NULL;
  try
  {
    comp = new AliHLTGlobalTriggerComponent;
  }
  catch (const std::bad_alloc&)
  {
    HLTError("Could not allocate memory for a new instance of AliHLTGlobalTriggerComponent.");
    return NULL;
  }
  return comp;
}


int AliHLTGlobalTriggerComponent::DoTrigger()
{
  // This method will apply the global trigger decision.

  if (fTrigger == NULL)
  {
    HLTFatal("Global trigger implementation object is NULL!");
    return -EIO;
  }

  AliHLTUInt32_t eventType=0;
  if (!IsDataEvent(&eventType)) {
    if (fDataEventsOnly)
    {
      if (eventType == gkAliEventTypeEndOfRun) PrintStatistics(fTrigger, kHLTLogImportant);
      IgnoreEvent();  // dont generate any trigger decision.
      return 0;
    }
  }
  
  fCDH = NULL;  // have to reset this in case ExtractTriggerData fails.
  ExtractTriggerData(*GetTriggerData(), NULL, NULL, &fCDH, NULL);

  // Copy the trigger counters in case we need to set them back to their original
  // value because the PushBack method fails with ENOSPC.
  TArrayL64 originalCounters = fTrigger->GetCounters();
  if (fTrigger->CallFailed()) return -EPROTO;
  
  fTrigger->NewEvent();
  if (fTrigger->CallFailed()) return -EPROTO;
  
  // Fill in the input data.
  const TObject* obj = GetFirstInputObject();
  while (obj != NULL)
  {
    fTrigger->Add(obj, GetDataType(), GetSpecification());
    if (fTrigger->CallFailed()) return -EPROTO;
    obj = GetNextInputObject();
  }

  // add trigger decisions for every CTP class
  const AliHLTCTPData* pCTPData=CTPData();
  if (pCTPData) {
    AddCTPDecisions(fTrigger, pCTPData, GetTriggerData());
  }
  
  bool softwareTriggerIsValid = FillSoftwareTrigger();
  if (softwareTriggerIsValid)
  {
    fTrigger->Add(&fSoftwareTrigger, kAliHLTDataTypeTriggerDecision, kAliHLTVoidDataSpec);
  }
  
  // Calculate the global trigger result and trigger domain, then create and push
  // back the new global trigger decision object.
  TString description;
  AliHLTTriggerDomain triggerDomain;
  bool triggerResult = false;
  bool matchedItems = fTrigger->CalculateTriggerDecision(triggerResult, triggerDomain, description);
  if (fTrigger->CallFailed()) return -EPROTO;
  AliHLTGlobalTriggerDecision decision(
      triggerResult,
      // The following will cause the decision to be generated with default values
      // (set in fTriggerDomain and fDescription) if the trigger result is false.
      (matchedItems == true) ? triggerDomain : GetTriggerDomain(),
      (matchedItems == true) ? description.Data() : GetDescription()
    );

  decision.SetUniqueID(fUniqueID);
  decision.SetCounters(fTrigger->GetCounters(), fTotalEventCounter+1);
  if (fTrigger->CallFailed()) return -EPROTO;
  
  TClonesArray shortInfo(TNamed::Class(), GetNumberOfInputBlocks());
  
  // Add the input objects used to make the global decision.
  obj = GetFirstInputObject();
  while (obj != NULL)
  {
    const AliHLTTriggerDecision* intrig = dynamic_cast<const AliHLTTriggerDecision*>(obj);
    
    if (TestBit(kForwardInput)) Forward(obj);
    
    if (TestBit(kIncludeInput))
    {
      if (intrig != NULL)
      {
         decision.AddTriggerInput(*intrig);
      }
      else
      {
        // The const_cast should be safe in this case because the list of inputObjects
        // will be added to the global decision with AddInputObjectRef, which only
        // modifies the kCanDelete bit and nothing else.
        // This is necessary since GetFirstInputObject only returns const objects.
        decision.AddInputObjectRef( const_cast<TObject*>(obj) );
      }
    }
    
    if (TestBit(kIncludeShort))
    {
      int entries = shortInfo.GetEntriesFast();
      try
      {
        new (shortInfo[entries]) TNamed(obj->GetName(), obj->GetTitle());
      }
      catch (const std::bad_alloc&)
      {
        HLTError("Could not allocate more memory for the short list of input objects.");
        return -ENOMEM;
      }
      if (intrig != NULL)
      {
        shortInfo[entries]->SetBit(BIT(16)); // indicate that this is a trigger decision
        shortInfo[entries]->SetBit(BIT(15), intrig->Result());
      }
    }

    obj = GetNextInputObject();
  }
  if (TestBit(kIncludeShort)) decision.AddInputObjectRef(&shortInfo);
  
  // The const_cast should be safe in this case because AddInputObjectRef just
  // modifies the kCanDelete bit and nothing else.
  if (!TestBit(kSkipCTP) && CTPData()) decision.AddInputObjectRef(const_cast<AliHLTCTPData*>(CTPData()));
  
  if (softwareTriggerIsValid and TestBit(kIncludeInput)) decision.AddTriggerInput(fSoftwareTrigger);
  
  static UInt_t lastTime=0;
  TDatime time;
  if (time.Get()-lastTime>60)
  {
    lastTime=time.Get();
    PrintStatistics(fTrigger, kHLTLogImportant);
  }
  else if (eventType==gkAliEventTypeEndOfRun)
  {
    PrintStatistics(fTrigger, kHLTLogImportant);
  }
  
  // add readout filter to event done data
  CreateEventDoneReadoutFilter(decision.TriggerDomain(), 3);

  // add monitoring filter to event done data if enabled by setting
  // a monitoring period >=0: -1 means off, 0 means for every event
  // configured by argument '-monitoring[=n]'
  if (fMonitorPeriod>=0) {
    static UInt_t lastMonitorEvent=0;

    AliHLTTriggerDomain monitoringFilter(decision.TriggerDomain());
    if (decision.Result() &&
	int(time.Get()-lastMonitorEvent)>fMonitorPeriod) {
      lastMonitorEvent=time.Get();
      // add monitoring event command for triggered events
      CreateEventDoneReadoutFilter(decision.TriggerDomain(), 5);
    } else {
      // empty filter list if events are not triggered
      // or within the monitoring interval
      monitoringFilter.Clear();
    }
    // add monitoring filter list
    CreateEventDoneReadoutFilter(monitoringFilter, 4);
  }
  
  // Mask the readout list according to the CTP trigger
  // if the classes have been initialized (mask non-zero).
  // If we are dealing with a software trigger on the other hand then
  // mask with the participating detector list.
  // In both cases we must make sure that HLT is part of the readout mask.
  if (CTPData() != NULL and CTPData()->Mask() != 0x0)
  {
    AliHLTReadoutList readoutlist = decision.ReadoutList();
    AliHLTReadoutList ctpreadout = CTPData()->ReadoutList(*GetTriggerData());
    ctpreadout.Enable(AliHLTReadoutList::kHLT);
    readoutlist.AndEq(ctpreadout);
    decision.ReadoutList(readoutlist); // override the readout list with the masked one.
  }
  else if (softwareTriggerIsValid)
  {
    assert(fCDH != NULL);
    AliHLTReadoutList readoutlist = decision.ReadoutList();
    UInt_t detectors = fCDH->GetSubDetectors();
    AliHLTReadoutList softwareReadout(Int_t(detectors | AliHLTReadoutList::kHLT));
    readoutlist.AndEq(softwareReadout);
    decision.ReadoutList(readoutlist); // override the readout list with the masked one.
  }

  if (TriggerEvent(&decision, kAliHLTDataTypeGlobalTrigger) == -ENOSPC)
  {
    // Increase the estimated buffer space required if the PushBack methods in TriggerEvent
    // returned the "no buffer space" error code. Also remember to set the trigger counters
    // back to what they were, otherwise triggers will be double counted when we try to reprocess
    // this event with larger buffers.
    fBufferSizeConst += 1024*1024;
    fBufferSizeMultiplier *= 2.;
    fTrigger->SetCounters(originalCounters);
    if (fTrigger->CallFailed()) return -EPROTO;
    return -ENOSPC;
  }
  
  ++fTotalEventCounter;
  return 0;
}


int AliHLTGlobalTriggerComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // Reconfigure the component by loading the trigger menu and recreating the
  // trigger logic class.
  
  const char* path = fgkTriggerMenuCDBPath;
  const char* id = "(unknown)";
  if (cdbEntry != NULL) path = cdbEntry;
  if (chainId != NULL and chainId[0] != '\0') id = chainId;
  HLTInfo("Reconfiguring from '%s' for chain component '%s'.", path, id);
  
  const AliHLTTriggerMenu* menu = NULL;
  int result = LoadTriggerMenu(path, menu);
  if (result != 0) return result;
  
  TString className;
  TString codeFileName;
  result = GenerateTrigger(menu, className, codeFileName, fIncludePaths, fIncludeFiles);
  if (result != 0) return result;
  
  AliHLTGlobalTrigger* trigger = NULL;
  try
  {
    trigger = AliHLTGlobalTrigger::CreateNew(className.Data());
  }
  catch (const std::bad_alloc&)
  {
    HLTError("Could not allocate memory for the AliHLTGlobalTrigger instance.");
    return -ENOMEM;
  }
  if (trigger == NULL)
  {
    HLTError("Could not create a new instance of '%s'.", className.Data());
    // Make sure to cleanup after the new code file.
    UnloadTriggerClass(codeFileName);
    if (not codeFileName.IsNull() and gSystem->AccessPathName(codeFileName)==0 and not fDebugMode)
    {
      codeFileName.ReplaceAll(".cxx", "*");
      TString command="rm "; command+=codeFileName;
      gSystem->Exec(command);
    }
    return -EIO;
  }
  
  if (fTrigger != NULL)
  {
    delete fTrigger;
    fTrigger = NULL;
  }
  
  fTrigger = trigger;
  fTrigger->FillFromMenu(*menu);
  if (fTrigger->CallFailed()) return -EPROTO;

  // Set the default values from the trigger menu.
  SetDescription(menu->DefaultDescription());
  SetTriggerDomain(menu->DefaultTriggerDomain());
  
  // Cleanup the old class code.
  UnloadTriggerClass(fCodeFileName);
  if (fDeleteCodeFile and not fCodeFileName.IsNull() and gSystem->AccessPathName(fCodeFileName)==0 and not fDebugMode)
  {
    fCodeFileName.ReplaceAll(".cxx", "*");
    TString command="rm "; command+=fCodeFileName;
    gSystem->Exec(command);
  }
  fCodeFileName = codeFileName;
  fDeleteCodeFile = true;  // Since we generated a new class.
  
  return 0;
}


int AliHLTGlobalTriggerComponent::LoadTriggerMenu(const char* cdbPath, const AliHLTTriggerMenu*& menu)
{
  // Loads the trigger menu object from the CDB path.
  
  HLTDebug("Trying to load trigger menu from '%s'.", cdbPath);
  TObject* obj = LoadAndExtractOCDBObject(cdbPath);
  if (obj == NULL)
  {
    HLTError("Configuration object for \"%s\" is missing.", cdbPath);
    return -ENOENT;
  }
  if (obj->IsA() != AliHLTTriggerMenu::Class())
  {
    HLTError("Wrong type for configuration object in \"%s\". Found a %s but we expect a AliHLTTriggerMenu.",
             cdbPath, obj->ClassName()
    );
    return -EPROTO;
  }
  menu = static_cast<AliHLTTriggerMenu*>(obj);
  return 0;
}


void AliHLTGlobalTriggerComponent::GenerateFileName(TString& name, TString& filename)
{
  // Creates a unique file name for the generated code.
  
  TUUID guid = GenerateGUID();
  union
  {
    UChar_t buf[16];
    UInt_t bufAsInt[4];
  };
  guid.GetUUID(buf);
  fUniqueID = bufAsInt[0];
  TString guidstr = guid.AsString();
  // Replace '-' with '_' in string.
  for (int i = 0; i < guidstr.Length(); ++i)
  {
    if (guidstr[i] == '-') guidstr[i] = '_';
  }
  name = "AliHLTGlobalTriggerImpl_";
  name += guidstr;
  filename = name + ".cxx";

  // For good measure, check that the file names are not used. If they are then regenerate.
  fstream file(filename.Data(), ios_base::in);
  if (file.good())
  {
    file.close();
    GenerateFileName(name, filename);
  }
}


int AliHLTGlobalTriggerComponent::GenerateTrigger(
    const AliHLTTriggerMenu* menu, TString& name, TString& filename,
    const TClonesArray& includePaths, const TClonesArray& includeFiles
  )
{
  // Generates the global trigger class that will implement the specified trigger menu.
  // See header for more details.
  
  GenerateFileName(name, filename);
  HLTDebug("Generating custom HLT trigger class named %s, in file %s, using trigger menu %p.",
    name.Data(), filename.Data(), ((void*)menu)
  );
  
  // Open a text file to write the code and generate the new class.
  fstream code(filename.Data(), ios_base::out | ios_base::trunc);
  if (not code.good())
  {
    HLTError("Could not open file '%s' for writing.", filename.Data());
    return -EIO;
  }
  
  TClonesArray symbols(AliHLTTriggerMenuSymbol::Class());
  int result = BuildSymbolList(menu, symbols);
  if (result != 0) return result;
  
  code << "#if !defined(__CINT__) || defined(__MAKECINT__)" << endl;
  code << "#include <cstring>" << endl;
  code << "#include \"TClass.h\"" << endl;
  code << "#include \"TString.h\"" << endl;
  code << "#include \"TClonesArray.h\"" << endl;
  code << "#include \"TRandom3.h\"" << endl;
  code << "#include \"AliHLTLogging.h\"" << endl;
  code << "#include \"AliHLTGlobalTrigger.h\"" << endl;
  code << "#include \"AliHLTGlobalTriggerDecision.h\"" << endl;
  code << "#include \"AliHLTDomainEntry.h\"" << endl;
  code << "#include \"AliHLTTriggerDomain.h\"" << endl;
  code << "#include \"AliHLTReadoutList.h\"" << endl;
  code << "#include \"AliHLTTriggerMenu.h\"" << endl;
  code << "#include \"AliHLTTriggerMenuItem.h\"" << endl;
  code << "#include \"AliHLTTriggerMenuSymbol.h\"" << endl;
  
  // Add any include files that were specified on the command line.
  for (Int_t i = 0; i < includeFiles.GetEntriesFast(); i++)
  {
    TString file = static_cast<TObjString*>(includeFiles.UncheckedAt(i))->String();
    code << "#include \"" << file.Data() << "\"" << endl;
  }
  
  if (fDebugMode)
  {
    code << "#else" << endl;
    code << "const char* gFunctionName = \"???\";" << endl;
    code << "#define HLTDebug(msg) if (CheckFilter(kHLTLogDebug) && CheckGroup(Class_Name())) SendMessage(kHLTLogDebug, Class_Name(), ::gFunctionName, __FILE__, __LINE__, msg)" << endl;
  }
  code << "#endif" << endl;
  
  code << "class " << name << " :" << endl;
  // Add appropriate #ifdef sections since we need to prevent inheritance from
  // AliHLTGlobalTrigger. CINT does not seem to support multiple inheritance nor
  // multiple levels of inheritance. Neither of the following schemes worked:
  //
  // 1)  class AliHLTGlobalTrigger : public AliHLTLogging {};
  //     class AliHLTGlobalTriggerImpl_xyz : public AliHLTGlobalTrigger {};
  //
  // 2)  class AliHLTGlobalTrigger {};
  //     class AliHLTGlobalTriggerImpl_xyz : public AliHLTGlobalTrigger, public AliHLTLogging {};
  //
  // Thus, we are forced to just inherit from AliHLTLogging when running in the CINT
  // interpreter. But we anyway have to call the global trigger implementation class
  // through the AliHLTGlobalTriggerWrapper so this is not such a problem.
  code << "#if !defined(__CINT__) || defined(__MAKECINT__)" << endl;
  code << "  public AliHLTGlobalTrigger," << endl;
  code << "#endif" << endl;
  code << "  public AliHLTLogging" << endl;
  code << "{" << endl;
  code << "public:" << endl;
  
  // Generate constructor method.
  code << "  " << name << "() :" << endl;
  code << "#if !defined(__CINT__) || defined(__MAKECINT__)" << endl;
  code << "    AliHLTGlobalTrigger()," << endl;
  code << "#endif" << endl;
  code << "    AliHLTLogging()";
  // Write the symbols in the trigger menu in the initialisation list.
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    code << "," << endl;
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    code << "    " << symbol->Name() << "()," << endl;
    if (strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0)
    {
      code << "    " << symbol->Name() << "TriggerDomain()," << endl;
    }
    code << "    " << symbol->Name() << "DomainEntry(kAliHLTAnyDataType)";
  }
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    code << "," << endl << "    fMenuItemDescription" << i << "()";
  }
  code << endl << "  {" << endl;
  if (fDebugMode)
  {
    code << "#ifdef __CINT__" << endl;
    code << "    gFunctionName = \"" << name.Data() <<"\";" << endl;
    code << "#endif" << endl;
    code << "    SetLocalLoggingLevel(kHLTLogAll);" << endl;
    code << "    HLTDebug(Form(\"Creating new instance at %p.\", this));" << endl;
  }
  code << "  }" << endl;
  
  code << "  virtual ~" << name << "() {" << endl;
  if (fDebugMode)
  {
    code << "#ifdef __CINT__" << endl;
    code << "    gFunctionName = \"~" << name.Data() << "\";" << endl;
    code << "#endif" << endl;
    code << "    HLTDebug(Form(\"Deleting instance at %p.\", this));" << endl;
  }
  code << "  }" << endl;
  
  // Generate the FillFromMenu method.
  code << "  virtual void FillFromMenu(const AliHLTTriggerMenu& menu) {" << endl;
  if (fDebugMode)
  {
    code << "#ifdef __CINT__" << endl;
    code << "    gFunctionName = \"FillFromMenu\";" << endl;
    code << "#endif" << endl;
    code << "    HLTDebug(Form(\"Filling description entries from trigger menu for global trigger %p.\", this));" << endl;
  }
  code << "    fCounter.Set(menu.NumberOfItems());" << endl;
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    code << "    fMenuItemDescription" << i << " = (menu.Item(" << i
         << ") != NULL) ? menu.Item(" << i << ")->Description() : \"\";" << endl;
  }
  if (fDebugMode)
  {
    code << "    HLTDebug(Form(\"Finished filling description entries from trigger menu.\"));" << endl;
    code << "    HLTDebug(Form(\"Filling domain entries from trigger menu symbols for global trigger %p.\", this));" << endl;
  }
  code << "    for (Int_t i = 0; i < menu.SymbolArray().GetEntriesFast(); i++) {" << endl;
  // 30 Oct 2009 - CINT sometimes evaluates the dynamic_cast incorrectly.
  // Have to use the TClass system for extra protection.
  code << "      if (menu.SymbolArray().UncheckedAt(i) == NULL) continue;" << endl;
  code << "      if (menu.SymbolArray().UncheckedAt(i)->IsA() != AliHLTTriggerMenuSymbol::Class()) continue;" << endl;
  code << "      const AliHLTTriggerMenuSymbol* symbol = dynamic_cast<const"
           " AliHLTTriggerMenuSymbol*>(menu.SymbolArray().UncheckedAt(i));" << endl;
  code << "      if (symbol == NULL) continue;" << endl;
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    code << "      if (strcmp(symbol->RealName(), \"" << symbol->RealName() << "\") == 0) {" << endl;
    if (fDebugMode)
    {
      code << "        HLTDebug(Form(\"Assinging domain entry value corresponding with symbol '%s' to '%s'.\","
              " symbol->RealName(), symbol->BlockType().AsString().Data()));" << endl;
    }
    code << "        " << symbol->Name() << "DomainEntry = symbol->BlockType();" << endl;
    code << "        continue;" << endl;
    code << "      }" << endl;
  }
  code << "    }" << endl;
  // The following is an optimisation where symbols without any assignment operators
  // are treated as constant and only initialised in FillFromMenu rather than reseting
  // them in the NewEvent method.
  // Note: we putting this initialisation into the constructor can lead to seg faults
  // under CINT interpretation. Thus we must put it into the FillFromMenu method instead.
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    if (TString(symbol->AssignExpression()) != "") continue;
    if (strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0) continue;
    // CINT has problems with the implicit equals operator for complex types, so if
    // the type has an equals operater we need to write the operator call explicitly.
    TClass* clas = TClass::GetClass(symbol->Type());
    if (clas != NULL and clas->GetMethodAny("operator=") != NULL)
    {
      code << "    " << symbol->Name() << ".operator = (" << symbol->DefaultValue() << ");" << endl;
    }
    else
    {
      code << "    " << symbol->Name() << " = " << symbol->DefaultValue() << ";" << endl;
    }
  }
  if (fDebugMode)
  {
    code << "    HLTDebug(Form(\"Finished filling domain entries from trigger menu symbols.\"));" << endl;
  }
  code << "  }" << endl;
  
  // Generate the NewEvent method.
  code << "  virtual void NewEvent() {" << endl;
  if (fDebugMode)
  {
    code << "#ifdef __CINT__" << endl;
    code << "    gFunctionName = \"NewEvent\";" << endl;
    code << "#endif" << endl;
    code << "    HLTDebug(Form(\"New event for global trigger object %p, initialising variables to default values.\", this));" << endl;
  }
  // Write code to initialise the symbols in the trigger menu to their default values.
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    // The following is an optimisation. If the symbol does not have an assignment expression
    // then it is effectively a constant symbol and can be initialised earlier and only once.
    // In this case we initialise it in the FillFromMenu method instead.
    if (TString(symbol->AssignExpression()) == "") continue;
    // CINT has problems with the implicit equals operator for complex types, so if
    // the type has an equals operater we need to write the operator call explicitly.
    TClass* clas = TClass::GetClass(symbol->Type());
    if (clas != NULL and clas->GetMethodAny("operator=") != NULL)
    {
      code << "    " << symbol->Name() << ".operator = (" << symbol->DefaultValue() << ");" << endl;
    }
    else
    {
      code << "    " << symbol->Name() << " = " << symbol->DefaultValue() << ";" << endl;
    }
    if (strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0)
    {
      code << "    " << symbol->Name() << "TriggerDomain.Clear();" << endl;
    }
  }
  if (fDebugMode)
  {
    code << "    HLTDebug(Form(\"Finished initialising variables.\"));" << endl;
  }
  code << "  }" << endl;
  
  // Generate the Add method.
  bool haveAssignments = false;
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    // First check if we have any symbols with assignment expressions.
    // This is needed to get rid of the on the fly compilation warning about '_object_' not being used.
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    TString expr = symbol->AssignExpression();
    if (expr == "") continue; // Skip entries that have no assignment expression.
    haveAssignments = true;
    break;
  }
  if (haveAssignments or fDebugMode)
  {
    code << "  virtual void Add(const TObject* _object_, const AliHLTComponentDataType& _type_, AliHLTUInt32_t _spec_) {" << endl;
  }
  else
  {
    code << "  virtual void Add(const TObject* /*_object_*/, const AliHLTComponentDataType& _type_, AliHLTUInt32_t _spec_) {" << endl;
  }
  if (fDebugMode)
  {
    code << "#ifdef __CINT__" << endl;
    code << "    gFunctionName = \"Add\";" << endl;
    code << "#endif" << endl;
  }
  code << "    AliHLTDomainEntry _type_spec_(_type_, _spec_);" << endl;
  if (fDebugMode)
  {
    code << "    HLTDebug(Form(\"Adding TObject %p, with class name '%s' from data block"
            " '%s', to global trigger object %p\", _object_, _object_->ClassName(),"
            " _type_spec_.AsString().Data(), this));" << endl;
    code << "    _object_->Print();" << endl;
    code << "    bool _object_assigned_ = false;" << endl;
  }
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    // Write code to check if the block type, specification and class name is correct.
    // Then write code to assign the variable from the input object.
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    TString expr = symbol->AssignExpression();
    if (expr == "") continue; // Skip entries that have no assignment expression.
    bool isTrigDecision = strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0;
    if (fDebugMode)
    {
      if (isTrigDecision)
      {
        code << "    HLTDebug(Form(\"Trying to match input object to class '"
             << symbol->ObjectClass() << "', trigger name '" << symbol->RealName()
             << "' and block type '%s'\", " << symbol->Name()
             << "DomainEntry.AsString().Data()));" << endl;
      }
      else
      {
        code << "    HLTDebug(Form(\"Trying to match input object to class '"
             << symbol->ObjectClass() << "' and block type '%s'\", "
             << symbol->Name() << "DomainEntry.AsString().Data()));" << endl;
      }
    }
    // 30 Oct 2009 - CINT sometimes evaluates the dynamic_cast incorrectly.
    // Have to use the TClass system for extra protection.
    code << "    const " << symbol->ObjectClass() << "* " << symbol->Name() << "_object_ = NULL;" << endl;
    code << "    if (_object_->InheritsFrom(" << symbol->ObjectClass() << "::Class())) " << symbol->Name()
         << "_object_ = dynamic_cast<const " << symbol->ObjectClass()
         << "*>(_object_);" << endl;
    code << "    if (" << symbol->Name() << "_object_ != NULL && ";
    if (isTrigDecision)
    {
      code << "strcmp(" << symbol->Name() << "_object_->Name(), \""
           << symbol->RealName() << "\") == 0 && ";
    }
    code << symbol->Name() << "DomainEntry == _type_spec_) {" << endl;
    TString fullname = symbol->Name();
    fullname += "_object_";
    expr.ReplaceAll("this", fullname);
    code << "      this->" << symbol->Name() << " = " << expr.Data() << ";" << endl;
    if (isTrigDecision)
    {
      code << "      this->" << symbol->Name() << "TriggerDomain = "
           << fullname.Data() << "->TriggerDomain();" << endl;
    }
    if (fDebugMode)
    {
      code << "      HLTDebug(Form(\"Added TObject %p with class name '%s' to variable "
           << symbol->Name() << "\", _object_, _object_->ClassName()));" << endl;
      code << "      _object_assigned_ = true;" << endl;
    }
    code << "    }" << endl;
  }
  if (fDebugMode)
  {
    code << "    if (! _object_assigned_) {" << endl;
    code << "      HLTDebug(Form(\"Did not assign TObject %p"
            " with class name '%s' to any variable.\", _object_, _object_->ClassName()));"
         << endl;
    code << "    }" << endl;
  }
  code << "  }" << endl;
  
  // Generate the CalculateTriggerDecision method.
  // This requires code to be generated that checks which items in the trigger menu
  // have their conditions asserted and then the trigger domain is generated from
  // those fired items.
  // The processing will start from the highest priority trigger group and stop
  // after at least one trigger from the current priority group being processed
  // is positive. For each priority group all the trigger menu items are checked.
  // Their combined trigger condition expression must be true for the trigger priority
  // group to be triggered positive. The full condition expression is formed by
  // concatenating the individual condition expressions. If no trailing operators are
  // used in the individual expressions then the default condition operator is placed
  // between two concatenated condition expressions.
  // If a trigger priority group has at least one trigger fired then the trigger domain
  // is calculated such that it will give the same result as the concatenated trigger
  // domain merging expressions for all the individual trigger menu items with
  // positive results. Again, if no trailing operators are used in the individual
  // merging expressions then the default domain operator is placed between two
  // expression fragments.
  code << "  virtual bool CalculateTriggerDecision(bool& _trigger_result_, AliHLTTriggerDomain& _domain_, TString& _description_) {" << endl;
  if (fDebugMode)
  {
    code << "#ifdef __CINT__" << endl;
    code << "    gFunctionName = \"CalculateTriggerDecision\";" << endl;
    code << "#endif" << endl;
    code << "    HLTDebug(Form(\"Calculating global HLT trigger result with trigger object at %p.\", this));" << endl;
  }
  
  // Build a list of priorities used in the trigger menu.
  std::vector<UInt_t> priorities;
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    const AliHLTTriggerMenuItem* item = menu->Item(i);
    bool priorityNotInList = std::find(priorities.begin(), priorities.end(), item->Priority()) == priorities.end();
    if (priorityNotInList) priorities.push_back(item->Priority());
  }
  std::sort(priorities.begin(), priorities.end(), AliHLTDescendingNumbers);
  // From the priority list, build the priority groups in the correct order,
  // i.e. highest priority first.
  // The priority group is a list of vectors of integers. The integers are the
  // index numbers into the trigger menu item list for the items which form part
  // of the priority group.
  std::vector<std::vector<Int_t> > priorityGroup;
  priorityGroup.insert(priorityGroup.begin(), priorities.size(), std::vector<Int_t>());
  for (size_t n = 0; n < priorities.size(); n++)
  {
    UInt_t priority = priorities[n];
    for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
    {
      const AliHLTTriggerMenuItem* item = menu->Item(i);
      if (item->Priority() == priority) priorityGroup[n].push_back(i);
    }
  }
  
  for (size_t n = 0; n < priorityGroup.size(); n++)
  {
    if (fDebugMode)
    {
      code << "    HLTDebug(Form(\"Processing trigger priority group " << priorities[n] << "\"));" << endl;
    }
    code << "    ";
    if (n == 0) code << "UInt_t ";
    code << "_previous_match_ = 0xFFFFFFFF;" << endl;
    code << "    ";
    if (n == 0) code << "bool ";
    code << "_trigger_matched_ = false;" << endl;
    code << "    ";
    if (n == 0) code << "bool ";
    code << "_group_result_ = false;" << endl;
    std::vector<TString> conditionOperator;
    conditionOperator.insert(conditionOperator.begin(), priorityGroup[n].size(), TString(""));
    std::vector<TString> domainOperator;
    domainOperator.insert(domainOperator.begin(), priorityGroup[n].size(), TString(""));
    for (size_t m = 0; m < priorityGroup[n].size(); m++)
    {
      UInt_t i = priorityGroup[n][m];
      const AliHLTTriggerMenuItem* item = menu->Item(i);
      TString triggerCondition = item->TriggerCondition();
      TString mergeExpr = item->MergeExpression();
      // Replace the symbols found in the trigger condition and merging expressions
      // with appropriate compilable versions.
      for (Int_t j = 0; j < symbols.GetEntriesFast(); j++)
      {
        AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(j) );
        bool symbolNamesDifferent = strcmp(symbol->RealName(), symbol->Name()) != 0;
        if (strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0)
        {
          TString newname = symbol->Name();
          newname += "TriggerDomain";
          mergeExpr.ReplaceAll(symbol->RealName(), newname);
        }
        else
        {
          if (symbolNamesDifferent) mergeExpr.ReplaceAll(symbol->RealName(), symbol->Name());
        }
        if (symbolNamesDifferent) triggerCondition.ReplaceAll(symbol->RealName(), symbol->Name());
      }
      // We allow the trigger conditions and merging expressions to have trailing operators.
      // Thus, we need to extract the operators and cleanup the expressions so that they will
      // compile. This means that we silently ignore the trailing operator if not needed.
      if (ExtractedOperator(triggerCondition, conditionOperator[m]))
      {
        // If the trailing operator is the same as the default operator then reset
        // the value in the operator list so that the default is used in the generated
        // code. This creates more compact code.
        if (conditionOperator[m] == menu->DefaultConditionOperator()) conditionOperator[m] = "";
      }
      if (ExtractedOperator(mergeExpr, domainOperator[m]))
      {
        if (domainOperator[m] == menu->DefaultDomainOperator()) domainOperator[m] = "";
      }
      if (fDebugMode)
      {
        code << "    HLTDebug(Form(\"Trying trigger condition " << i
             << " (Description = '%s').\", fMenuItemDescription" << i << ".Data()));"
             << endl;
      }
      code << "    ";
      if (n == 0 and m == 0) code << "bool ";
      code << "_item_result_ = false;" << endl;
      code << "    if (" << triggerCondition << ") {" << endl;
      code << "      ++fCounter[" << i << "];" << endl;
      const char* indentation = "";
      // Generate the code to handle the prescalar and scale-down
      bool havePrescalar = item->PreScalar() != 0;
      bool haveScaledown = item->ScaleDown() < 1;
      if (havePrescalar or haveScaledown) 
      {
        indentation = "  ";
        code << "      if (";
        if (havePrescalar) code << "(fCounter[" << i << "] % " << item->PreScalar() << ") == 1";
        if (havePrescalar and haveScaledown) code << " && ";
        if (haveScaledown)
        {
          std::streamsize oldprecision = code.precision(17);
          code << "gRandom->Rndm() < " << item->ScaleDown();
          code.precision(oldprecision);
        }
        code << ") {" << endl;
      }
      code << indentation << "      _item_result_ = true;" << endl;
      if (fDebugMode)
      {
        code << indentation << "      HLTDebug(Form(\"Matched trigger condition " << i
             << " (Description = '%s').\", fMenuItemDescription" << i << ".Data()));" << endl;
      }
      if (havePrescalar or haveScaledown)
      {
        code << "      }" << endl;
      }
      code << "    }" << endl;
      if (m == 0)
      {
        // Since this is the first item of the trigger group,
        // the generated trigger logic can be simplified a little.
        code << "    _group_result_ = _item_result_;" << endl;
        code << "    if (_item_result_) {" << endl;
        code << "      _domain_ = " << mergeExpr.Data() << ";" << endl;
        code << "      _description_ = fMenuItemDescription" << i << ";" << endl;
        code << "      _previous_match_ = " << i << ";" << endl;
        code << "      _trigger_matched_ = true;" << endl;
        code << "    }" << endl;
      }
      else
      {
        if (conditionOperator[m-1] == "")
        {
          code << "    _group_result_ = _group_result_ "
               << menu->DefaultConditionOperator() << " _item_result_;" << endl;
        }
        else
        {
          code << "    _group_result_ = _group_result_ "
               << conditionOperator[m-1] << " _item_result_;" << endl;
        }
        code << "    if (_item_result_) {" << endl;
        code << "      if (_trigger_matched_) {" << endl;
        bool switchWillBeEmpty = true;
        for (size_t k = 0; k < m; k++)
        {
          if (domainOperator[k] == "") continue;
          switchWillBeEmpty = false;
        }
        if (switchWillBeEmpty)
        {
          code << "        _domain_ = _domain_ " << menu->DefaultDomainOperator() << " "
               << mergeExpr.Data() << ";" << endl;
        }
        else
        {
          code << "        switch(_previous_match_) {" << endl;
          for (size_t k = 0; k < m; k++)
          {
            if (domainOperator[k] == "") continue;
            code << "        case " << k << ": _domain_ = _domain_ "
                 << domainOperator[k] << " " << mergeExpr.Data() << "; break;" << endl;
          }
          code << "        default: _domain_ = _domain_ "
               << menu->DefaultDomainOperator() << " " << mergeExpr.Data() << ";" << endl;
          code << "        }" << endl;
        }
        code << "        _description_ += \",\";" << endl;
        code << "        _description_ += fMenuItemDescription" << i << ";" << endl;
        code << "      } else {" << endl;
        code << "        _domain_ = " << mergeExpr.Data() << ";" << endl;
        code << "        _description_ = fMenuItemDescription" << i << ";" << endl;
        code << "      }" << endl;
        code << "      _previous_match_ = " << i << ";" << endl;
        code << "      _trigger_matched_ = true;" << endl;
        code << "    }" << endl;
      }
    }
    code << "    if (_group_result_) {" << endl;
    if (fDebugMode)
    {
      if (n < priorities.size() - 1)
      {
        code << "      HLTDebug(Form(\"Matched triggers in trigger priority group " << priorities[n]
             << ". Stopping processing here because all other trigger groups have lower priority.\"));" << endl;
      }
      else
      {
        code << "      HLTDebug(Form(\"Matched triggers in trigger priority group " << priorities[n] << ".\"));" << endl;
      }
    }
    bool methodReturnResult = true;
    if (priorityGroup[n].size() > 0)
    {
      const AliHLTTriggerMenuItem* item = menu->Item(priorityGroup[n][0]);
      methodReturnResult = item->DefaultResult();
    }
    // Check to see if the items of the group all have the same default result.
    // If not then warn the user since only the first item's value will be used.
    for (size_t m = 1; m < priorityGroup[n].size(); m++)
    {
      const AliHLTTriggerMenuItem* item = menu->Item(priorityGroup[n][m]);
      if (item->DefaultResult() != methodReturnResult)
      {
        HLTWarning("Found items with different default results set for priority group %d."
                   "Will only use the value from the first item.",
                   item->Priority()
                  );
        break;
      }
    }
    code << "      _trigger_result_ = " << (methodReturnResult ? "true" : "false") << ";" << endl;
    code << "      return true;" << endl;
    code << "    }" << endl;
  }
  code << "    _domain_.Clear();" << endl;
  code << "    _description_ = \"\";" << endl;
  code << "    _trigger_result_ = " << (menu->DefaultResult() ? "true" : "false") << ";" << endl;
  code << "    return false;" << endl;
  code << "  }" << endl;
  
  // Generate getter and setter methods for the counters.
  code << "  const TArrayL64& GetCounters() const { return fCounter; }" << endl;
  code << "  void SetCounters(const TArrayL64& counters) { fCounter = counters; }" << endl;
  
  code << "private:" << endl;
  // Add the symbols in the trigger menu to the list of private variables.
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    code << "  " << symbol->Type() << " " << symbol->Name() << ";" << endl;
    if (strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0)
    {
      code << "  AliHLTTriggerDomain " << symbol->Name() << "TriggerDomain;" << endl;
    }
    code << "  AliHLTDomainEntry " << symbol->Name() << "DomainEntry;" << endl;
  }
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    code << "  TString fMenuItemDescription" << i << ";" << endl;
  }
  code << "  TArrayL64 fCounter;" << endl;
  code << "#if !defined(__CINT__) || defined(__MAKECINT__)" << endl;
  code << "  ClassDef(" << name.Data() << ", 0)" << endl;
  code << "#else" << endl;
  code << "  virtual const char* Class_Name() const { return \"" << name.Data() << "\"; }" << endl;
  code << "#endif" << endl;
  code << "};" << endl;
  code << "#if !defined(__CINT__) || defined(__MAKECINT__)" << endl;
  code << "ClassImp(" << name.Data() << ")" << endl;
  code << "#endif" << endl;
  
  code.close();
  
  // Now we need to compile and load the new class.
  result = LoadTriggerClass(filename, includePaths);
  return result;
}


int AliHLTGlobalTriggerComponent::LoadTriggerClass(
    const char* filename, const TClonesArray& includePaths
  )
{
  // Loads the code for a custom global trigger class implementation on the fly.
  
  HLTDebug("Loading HLT trigger class from file '%s'.", filename);
  
  TString compiler = gSystem->GetBuildCompilerVersion();
  if (fRuntimeCompile && (compiler.Contains("gcc") or compiler.Contains("icc")))
  {
    TString includePath;
#if defined(PKGINCLUDEDIR)
    // this is especially for the HLT build system where the package is installed
    // in a specific directory including proper treatment of include files
    includePath.Form("-I%s", PKGINCLUDEDIR);
#else
    // the default AliRoot behavior, all include files can be found in the
    // $ALICE_ROOT subfolders
    includePath = "-I${ALICE_ROOT}/include -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/trigger";
#endif
    // Add any include paths that were specified on the command line.
    for (Int_t i = 0; i < includePaths.GetEntriesFast(); i++)
    {
      TString path = static_cast<TObjString*>(includePaths.UncheckedAt(i))->String();
      includePath += " -I";
      includePath += path;
    }
    HLTDebug("using include settings: %s", includePath.Data());
    gSystem->SetIncludePath(includePath);
    gSystem->SetFlagsOpt("-O3 -DNDEBUG");
    gSystem->SetFlagsDebug("-g3 -DDEBUG -D__DEBUG");
    
    int result = kTRUE;
    if (fDebugMode)
    {
      result = gSystem->CompileMacro(filename, "g");
    }
    else
    {
      result = gSystem->CompileMacro(filename, "O");
    }
    if (result != kTRUE)
    {
      HLTFatal("Could not compile and load global trigger menu implementation.");
      return -ENOENT;
    }
  }
  else
  {
    // Store the library state to be checked later in UnloadTriggerClass.
    fLibStateAtLoad = gSystem->GetLibraries();
    
    // If we do not support the compiler then try interpret the class instead.
    TString cmd = ".L ";
    cmd += filename;
    Int_t errorcode = TInterpreter::kNoError;
    gROOT->ProcessLine(cmd, &errorcode);
    if (errorcode != TInterpreter::kNoError)
    {
      HLTFatal("Could not load interpreted global trigger menu implementation"
               " (Interpreter error code = %d).",
               errorcode
      );
      return -ENOENT;
    }
  }
  
  return 0;
}


int AliHLTGlobalTriggerComponent::UnloadTriggerClass(const char* filename)
{
  // Unloads the code previously loaded by LoadTriggerClass.
  
  HLTDebug("Unloading HLT trigger class in file '%s'.", filename);
  
  TString compiler = gSystem->GetBuildCompilerVersion();
  if (fRuntimeCompile && (compiler.Contains("gcc") or compiler.Contains("icc")))
  {
    // Generate the library name.
    TString libname = filename;
    Ssiz_t dotpos = libname.Last('.');
    if (0 <= dotpos and dotpos < libname.Length()) libname[dotpos] = '_';
    libname += ".";
    libname += gSystem->GetSoExt();
    
    // This is a workaround for a problem with unloading shared libraries in ROOT.
    // If the trigger logic library is loaded before the libAliHLTHOMER.so library
    // or any other library is loaded afterwards, then during the gInterpreter->UnloadFile
    // call all the subsequent libraries get unloded. This means that any objects created
    // from classes implemented in the libAliHLTHOMER.so library will generate segfaults
    // since the executable code has been unloaded.
    // We need to check if there are any more libraries loaded after the class we
    // are unloading and in that case don't unload the class.
    TString libstring = gSystem->GetLibraries();
    TString token, lastlib;
    Ssiz_t from = 0;
    Int_t numOfLibs = 0, posOfLib = -1;
    while (libstring.Tokenize(token, from, " "))
    {
      ++numOfLibs;
      lastlib = token;
      if (token.Contains(libname)) posOfLib = numOfLibs;
    }
    if (numOfLibs != posOfLib)
    {
      HLTWarning(Form("ROOT limitation! Cannot properly cleanup and unload the shared"
          " library '%s' since another library '%s' was loaded afterwards. Trying to"
          " unload this library will remove the others and lead to serious memory faults.",
          libname.Data(), lastlib.Data()
      ));
      return 0;
    }
    
    char* path = NULL;
    int result = 0;
    if ((path = gSystem->DynamicPathName(libname)) != NULL)
    {
      result = gInterpreter->UnloadFile(path);
      delete [] path;
    }
    if (result != TInterpreter::kNoError) return -ENOENT;
  }
  else
  {
    // This is again a workaround for the problem with unloading files in ROOT.
    // If the trigger logic class is loaded before the libAliHLTHOMER.so library
    // or any other library is loaded afterwards, then during the gInterpreter->UnloadFile
    // call all the subsequent libraries get unloded.
    // We need to check if the list of loaded libraries has changed since the last
    // call to LoadTriggerClass. If it has then don't unload the class.
    if (fLibStateAtLoad != gSystem->GetLibraries())
    {
      TString libstring = gSystem->GetLibraries();
      TString token;
      Ssiz_t from = 0;
      while (libstring.Tokenize(token, from, " "))
      {
        if (not fLibStateAtLoad.Contains(token)) break;
      }
      HLTWarning(Form("ROOT limitation! Cannot properly cleanup and unload the file"
          " '%s' since another library '%s' was loaded afterwards. Trying to unload"
          " this file will remove the other library and lead to serious memory faults.",
          filename, token.Data()
      ));
      return 0;
    }
  
    // If we did not compile the trigger logic then remove the interpreted class.
    TString cmd = ".U ";
    cmd += filename;
    Int_t errorcode = TInterpreter::kNoError;
    gROOT->ProcessLine(cmd, &errorcode);
    if (errorcode != TInterpreter::kNoError)
    {
      HLTFatal("Could not unload interpreted global trigger menu implementation"
               " (Interpreter error code = %d).",
               errorcode
      );
      return -ENOENT;
    }
  }
  
  return 0;
}


int AliHLTGlobalTriggerComponent::FindSymbol(const char* name, const TClonesArray& list)
{
  // Searches for the named symbol in the given list.
  // See header for more details.
  
  for (int i = 0; i < list.GetEntriesFast(); i++)
  {
    const AliHLTTriggerMenuSymbol* symbol = dynamic_cast<const AliHLTTriggerMenuSymbol*>( list.UncheckedAt(i) );
    if (symbol == NULL) continue;
    if (strcmp(symbol->Name(), name) == 0) return i;
  }
  return -1;
}


namespace
{
  /**
   * Helper routine to compare two trigger menu symbols to see if b is a subset of a.
   * \returns true if b is a subset of a or they are unrelated.
   */
  bool AliHLTCheckForContainment(const AliHLTTriggerMenuSymbol* a, const AliHLTTriggerMenuSymbol* b)
  {
    TString bstr = b->Name();
    return bstr.Contains(a->Name());
  }

} // end of namespace


int AliHLTGlobalTriggerComponent::BuildSymbolList(const AliHLTTriggerMenu* menu, TClonesArray& list)
{
  // Builds the list of symbols to use in the custom global trigger menu
  // implementation class.
  // See header for more details.
  
  // Note: when we build the symbol list we must use the symbol name as returned
  // by the Name() method and not the RealName() method when using FindSymbol.
  // This is so that we avoid problems with the generated code not compiling
  // because names like "abc-xyz" and "abc_xyz" are synonymous.
  // Name() returns the converted C++ symbol name as used in the generated code.
  
  for (UInt_t i = 0; i < menu->NumberOfSymbols(); i++)
  {
    const AliHLTTriggerMenuSymbol* symbol = menu->Symbol(i);
    if (FindSymbol(symbol->Name(), list) != -1)
    {
      HLTError("Multiple symbols with the name '%s' defined in the trigger menu.", symbol->Name());
      return -EIO;
    }
    try
    {
      new (list[list.GetEntriesFast()]) AliHLTTriggerMenuSymbol(*symbol);
    }
    catch (const std::bad_alloc&)
    {
      HLTError("Could not allocate more memory for the symbols list when adding a trigger menu symbol.");
      return -ENOMEM;
    }
  }
  Int_t initialEntryCount = list.GetEntriesFast();
  
  // Note: the \\. must not be the first element in the character class, otherwise
  // it is interpreted as an "any character" dot symbol.
  TRegexp exp("[_a-zA-Z][-\\._a-zA-Z0-9]*");
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    const AliHLTTriggerMenuItem* item = menu->Item(i);
    TString str = item->TriggerCondition();
    Ssiz_t start = 0;
    do
    {
      Ssiz_t length = 0;
      Ssiz_t pos = exp.Index(str, &length, start);
      if (pos == kNPOS) break;
      start = pos+length;
      
      // Check if there is a numerical character before the found
      // regular expression. If so, then the symbol is not a valid one
      // and should be skipped.
      if (pos > 0)
      {
        bool notValid = false;
        switch (str[pos-1])
        {
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9':
          notValid = true;
          break;
        default:
          notValid = false;
          break;
        }
        if (notValid) continue;
      }
      TString s = str(pos, length);
      
      if (s == "and" or s == "and_eq" or s == "bitand" or s == "bitor" or
          s == "compl" or s == "not" or s == "not_eq" or s == "or" or
          s == "or_eq" or s == "xor" or s == "xor_eq" or s == "true" or
          s == "false"
         )
      {
        // Ignore iso646.h and other keywords.
        continue;
      }
      
      // We need to handle the special case where the symbol contains a dot.
      // In C++ this is a dereferencing operator. So we need to check if the
      // current symbol we are handling starts with the same string as any of
      // the existing symbols defined manually in the symbols table.
      // If we do find such a case then revert to treating the dot as an operator
      // rather than part of the symbol name. i.e. skip adding the automatic symbol.
      bool dereferencedSymbol = false;
      for (int j = 0; j < initialEntryCount; j++)
      {
        const AliHLTTriggerMenuSymbol* symbol = dynamic_cast<const AliHLTTriggerMenuSymbol*>( list.UncheckedAt(j) );
        if (symbol == NULL) continue;
        TString symstr = symbol->Name();
        symstr += ".";
        if (s.BeginsWith(symstr))
        {
          dereferencedSymbol = true;
          break;
        }
      }
      if (dereferencedSymbol) continue;

      // Need to create the symbols first and check if its name is in the list
      // before actually adding it to the symbols list.
      AliHLTTriggerMenuSymbol newSymbol;
      newSymbol.Name(s.Data());
      newSymbol.Type("bool");
      newSymbol.ObjectClass("AliHLTTriggerDecision");
      newSymbol.AssignExpression("this->Result()");
      newSymbol.DefaultValue("false");
      if (FindSymbol(newSymbol.Name(), list) == -1)
      {
        try
        {
          new (list[list.GetEntriesFast()]) AliHLTTriggerMenuSymbol(newSymbol);
        }
        catch (const std::bad_alloc&)
        {
          HLTError("Could not allocate more memory for the symbols list when adding a trigger name symbol.");
          return -ENOMEM;
        }
      }
    }
    while (start < str.Length());
  }
  
  // This last part is necessary to make sure that symbols are replaced in the
  // trigger condition and domain merging expressions in a greedy manner.
  // I.e. we need to make sure that if one symbol's string representation is
  // contained inside another (string subset) that the longer symbol name is
  // always first in the symbols list.
  // This will work because the symbol table is traversed from first to last
  // element and TString::ReplaceAll is used to replace the substrings inside
  // the AliHLTGlobalTriggerComponent::GenerateTrigger method.
  std::vector<AliHLTTriggerMenuSymbol*> orderedList;
  for (Int_t i = 0; i < list.GetEntriesFast(); i++)
  {
    orderedList.push_back( static_cast<AliHLTTriggerMenuSymbol*>(list.UncheckedAt(i)) );
  }
  std::sort(orderedList.begin(), orderedList.end(), AliHLTCheckForContainment);
  //std::sort(orderedList.begin(), orderedList.end());
  // Now swap values around according to the orderedList.
  for (Int_t i = 0; i < list.GetEntriesFast(); i++)
  {
    AliHLTTriggerMenuSymbol* target = static_cast<AliHLTTriggerMenuSymbol*>(list.UncheckedAt(i));
    AliHLTTriggerMenuSymbol tmp = *target;
    *target = *orderedList[i];
    *orderedList[i] = tmp;
  }
  
  return 0;
}


bool AliHLTGlobalTriggerComponent::ExtractedOperator(TString& expr, TString& op)
{
  // Extracts the trailing operator from the expression.
  
  Ssiz_t i = 0;
  // First skip the trailing whitespace.
  bool whitespace = true;
  for (i = expr.Length()-1; i >= 0 and whitespace; i--)
  {
    switch (expr[i])
    {
    case ' ': case '\t': case '\r': case '\n':
      whitespace = true;
      break;
    default:
      whitespace = false;
    }
  }
  if (i < 0 or whitespace) return false;
  
  // Now find the first whitespace character before the trailing symbol.
  bool nonwhitespace = true;
  for (; i >= 0 and nonwhitespace; i--)
  {
    switch (expr[i])
    {
    case ' ': case '\t': case '\r': case '\n':
      nonwhitespace = false;
      break;
    default:
      nonwhitespace = true;
    }
  }
  if (i < 0 or nonwhitespace) return false;
  
  // Extract the last symbols and check if it is a valid operator.
  TString s = expr;
  s.Remove(0, i+2);
  if (s == "and" or s == "and_eq" or s == "bitand" or s == "bitor" or
      s == "compl" or s == "not" or s == "not_eq" or s == "or" or
      s == "or_eq" or s == "xor" or s == "xor_eq" or s == "&&" or
      s == "&=" or s == "&" or s == "|" or s == "~" or s == "!" or
      s == "!=" or s == "||" or s == "|=" or s == "^" or s == "^=" or
      s == "==" or s == "+" or s == "-" or s == "*" or s == "/" or
      s == "%" or s == ">" or s == "<" or s == ">=" or s == "<="
     )
  {
    expr.Remove(i+1);
    op = s;
    return true;
  }
  
  return false;
}


bool AliHLTGlobalTriggerComponent::FillSoftwareTrigger()
{
  // Fills the fSoftwareTrigger structure.
  
  if (fCDH == NULL) return false;
  UChar_t l1msg = fCDH->GetL1TriggerMessage();
  if ((l1msg & 0x1) == 0x0) return false;  // skip physics events.
  // From here on everything must be a software trigger.
  if (((l1msg >> 2) & 0xF) == 0xE)
  {
    fSoftwareTrigger.Name("START_OF_DATA");
    fSoftwareTrigger.Description("Generated internal start of data trigger.");
  }
  else if (((l1msg >> 2) & 0xF) == 0xF)
  {
    fSoftwareTrigger.Name("END_OF_DATA");
    fSoftwareTrigger.Description("Generated internal end of data trigger.");
  }
  else if (((l1msg >> 6) & 0x1) == 0x1)
  {
    fSoftwareTrigger.Name("CALIBRATION");
    fSoftwareTrigger.Description("Generated internal calibration trigger.");
  }
  else
  {
    fSoftwareTrigger.Name("SOFTWARE");
    fSoftwareTrigger.Description("Generated internal software trigger.");
  }
  UInt_t detectors = fCDH->GetSubDetectors();
  fSoftwareTrigger.ReadoutList( AliHLTReadoutList(Int_t(detectors)) );
  return true;
}


int AliHLTGlobalTriggerComponent::PrintStatistics(const AliHLTGlobalTrigger* pTrigger, AliHLTComponentLogSeverity level, int offset) const
{
  // print some statistics
  int totalEvents=fTotalEventCounter+offset;
  const TArrayL64& counters = pTrigger->GetCounters();
  if (pTrigger->CallFailed()) return -EPROTO;
  TString msg;
  for (int i = 0; i < counters.GetSize(); i++) {
    ULong64_t count = counters[i];
    float ratio=0;
    if (totalEvents>0) ratio=100*(float)count/totalEvents;
    if (i != 0) msg += "\n";
    msg += Form("Item %d: total events: %d - counted events: %llu (%.1f%%)", i, totalEvents, count, ratio);
  }
  HLTLog(level, msg.Data());
  return 0;
}

int AliHLTGlobalTriggerComponent::AddCTPDecisions(AliHLTGlobalTrigger* pTrigger, const AliHLTCTPData* pCTPData, const AliHLTComponentTriggerData* trigData)
{
  // add trigger decisions for the valid CTP classes
  if (!pCTPData || !pTrigger) return 0;

  AliHLTUInt64_t triggerMask=pCTPData->Mask();
  AliHLTUInt64_t bit0=0x1;
  if (!fCTPDecisions) {
    try
    {
      fCTPDecisions=new TClonesArray(AliHLTTriggerDecision::Class(), gkNCTPTriggerClasses);
    }
    catch (const std::bad_alloc&)
    {
      HLTError("Could not allocate memory for the CTP decisions array.");
      return -ENOMEM;
    }
    if (!fCTPDecisions) return -ENOMEM;

    try
    {
      fCTPDecisions->ExpandCreate(gkNCTPTriggerClasses);
    }
    catch (const std::bad_alloc&)
    {
      HLTError("Could not allocate more memory for the CTP decisions array.");
      return -ENOMEM;
    }
    for (int i=0; i<gkNCTPTriggerClasses; i++) {
      const char* name=pCTPData->Name(i);
      if (triggerMask&(bit0<<i) && name) {
	AliHLTTriggerDecision* pDecision=dynamic_cast<AliHLTTriggerDecision*>(fCTPDecisions->At(i));
	assert(pDecision);
	if (!pDecision) {
	  delete fCTPDecisions;
	  fCTPDecisions=NULL;
	  return -ENOENT;
	}
	pDecision->Name(name);
      }
    }
  }

  for (int i=0; i<gkNCTPTriggerClasses; i++) {
    const char* name=pCTPData->Name(i);
    if ((triggerMask&(bit0<<i))==0 || name==NULL) continue;
    AliHLTTriggerDecision* pDecision=dynamic_cast<AliHLTTriggerDecision*>(fCTPDecisions->At(i));
    HLTDebug("updating CTP trigger decision %d %s (%p casted %p)", i, name, fCTPDecisions->At(i), pDecision);
    if (!pDecision) return -ENOENT;

    bool result=false;
    // 13 March 2010 - Optimisation:
    // Dont use the EvaluateCTPTriggerClass method, which uses slow TFormula objects.
    AliHLTUInt64_t triggers = 0;
    if (trigData) triggers = pCTPData->ActiveTriggers(*trigData);
    else triggers = pCTPData->Triggers();
    result = (triggers&((AliHLTUInt64_t)0x1<<i)) ? true : false;
    //if (trigData) result=pCTPData->EvaluateCTPTriggerClass(name, *trigData);
    //else result=pCTPData->EvaluateCTPTriggerClass(name);
    pDecision->Result(result);
    pDecision->TriggerDomain().Clear();
    if (trigData) pDecision->TriggerDomain().Add(pCTPData->ReadoutList(*trigData));
    else pDecision->TriggerDomain().Add(pCTPData->ReadoutList());

    pTrigger->Add(fCTPDecisions->At(i), kAliHLTDataTypeTriggerDecision, kAliHLTVoidDataSpec);
    if (pTrigger->CallFailed()) return -EPROTO;
  }

  return 0;
}
