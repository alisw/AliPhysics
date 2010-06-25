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

//  @file   AliHLTCTPData.cxx
//  @author Matthias Richter
//  @date   2009-08-20
//  @brief  Container for CTP trigger classes and counters
//  @note

#include "AliHLTCTPData.h"
#include "TClass.h"
#include "TObjString.h"
#include "TFormula.h"
#include "AliHLTComponent.h"
#include "AliRawDataHeader.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTCTPData)

AliHLTCTPData::AliHLTCTPData()
  : TNamed("AliHLTCTPData", "HLT counters for the CTP")
  , AliHLTLogging()
  , fMask(0)
  , fTriggers(0)
  , fClassIds(AliHLTReadoutList::Class(), gkNCTPTriggerClasses)
  , fCounters(gkNCTPTriggerClasses)
  , fMap()
{
  // constructor
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTCTPData::AliHLTCTPData(const char* parameter)
  : TNamed("AliHLTCTPData", "HLT counters for the CTP")
  , AliHLTLogging()
  , fMask(0)
  , fTriggers(0)
  , fClassIds(AliHLTReadoutList::Class(), gkNCTPTriggerClasses)
  , fCounters(gkNCTPTriggerClasses)
  , fMap()
{
  // constructor, init the CTP trigger classes
  InitCTPTriggerClasses(parameter);
}

AliHLTCTPData::~AliHLTCTPData()
{
  // destructor
  fClassIds.Delete();
}

AliHLTCTPData::AliHLTCTPData(const AliHLTCTPData& src)
  : TNamed(src.GetName(), src.GetTitle())
  , AliHLTLogging()
  , fMask(src.Mask())
  , fTriggers(src.fTriggers)
  , fClassIds(src.fClassIds)
  , fCounters(src.Counters())
  , fMap()
{
  // copy constructor
  ReadMap();
}

AliHLTCTPData& AliHLTCTPData::operator=(const AliHLTCTPData& src)
{
  // assignment operator, clone content
  if (this!=&src) {
    SetName(src.GetName());
    SetTitle(src.GetTitle());
    fMask=src.Mask();
    fClassIds.Delete();
    fClassIds.ExpandCreate(gkNCTPTriggerClasses);
    for (int i=0; i<gkNCTPTriggerClasses; i++) {
      if (i>src.fClassIds.GetLast()) break;
      ((TNamed*)fClassIds.At(i))->SetName(src.fClassIds.At(i)->GetName());
      ((TNamed*)fClassIds.At(i))->SetTitle(src.fClassIds.At(i)->GetTitle());
    }
    fCounters=src.Counters();
  }

  ReadMap();
  return *this;
}

int AliHLTCTPData::Add(const AliHLTCTPData& src, int factor, int &skipped)
{
  // see header file for class documentation
  
  skipped=0;
  for (int i=0; i<gkNCTPTriggerClasses; i++) {
    TString c;
    c=fClassIds.At(i)->GetName();
    if (c.IsNull()) continue;
    if (c.CompareTo(src.fClassIds.At(i)->GetName())==0) {
      fCounters[i]+=factor*src.Counter(i);
    } else {
      skipped++;
    }
  }
  return 0;
}

AliHLTCTPData& AliHLTCTPData::operator += (const AliHLTCTPData& src)
{
  // see header file for class documentation
  
  int nofInconsistencies=0;
  Add(src, 1, nofInconsistencies);
  if (nofInconsistencies>0) {
    HLTError("Inconsistent operants: skipping %d of %d CTP classes for operation", nofInconsistencies, gkNCTPTriggerClasses);
  }
  return *this;
}

AliHLTCTPData& AliHLTCTPData::operator -= (const AliHLTCTPData& src)
{
  // see header file for class documentation
  
  int nofInconsistencies=0;
  Add(src, -1, nofInconsistencies);
  if (nofInconsistencies>0) {
    HLTError("Inconsistent operants: skipping %d of %d CTP classes for operation", nofInconsistencies, gkNCTPTriggerClasses);
  }
  return *this;
}

AliHLTCTPData AliHLTCTPData::operator + (const AliHLTCTPData& src) const
{
  // see header file for class documentation

  AliHLTCTPData result(*this);
  result+=src;
  return result;
}

AliHLTCTPData AliHLTCTPData::operator - (const AliHLTCTPData& src) const
{
  // see header file for class documentation

  AliHLTCTPData result(*this);
  result-=src;
  return result;
}

int AliHLTCTPData::InitCTPTriggerClasses(const char* ctpString)
{
  // see header file for function documentation
  if (!ctpString) return -EINVAL;

  HLTImportant("Parameter: %s", ctpString);

  fMask=0;
  fClassIds.Delete();
  fClassIds.ExpandCreate(gkNCTPTriggerClasses);

  // general format of the CTP_TRIGGER_CLASS parameter
  // <bit position>:<Trigger class identifier string>:<detector-id-nr>-<detector-id-nr>-...,<bit position>:<Trigger class identifier string>:<detector-id-nr>-<detector-id-nr>-...,...
  HLTDebug(": %s", ctpString);
  TString string=ctpString;
  if (string.BeginsWith("CTP_TRIGGER_CLASS=")) string.ReplaceAll("CTP_TRIGGER_CLASS=", "");
  TObjArray* classEntries=string.Tokenize(",");
  if (classEntries) {
    enum {kBit=0, kName, kDetectors};
    for (int i=0; i<classEntries->GetEntriesFast(); i++) {
      TString entry=((TObjString*)classEntries->At(i))->GetString();
      TObjArray* entryParams=entry.Tokenize(":");
      if (entryParams) {
	if (entryParams->GetEntriesFast()==3 &&
	    (((TObjString*)entryParams->At(kBit))->GetString()).IsDigit()) {
	  int index=(((TObjString*)entryParams->At(kBit))->GetString()).Atoi();
	  if (index<gkNCTPTriggerClasses) {
	    AliHLTReadoutList* pCTPClass=dynamic_cast<AliHLTReadoutList*>(fClassIds.At(index));
	    if (pCTPClass) {
	      fMask|=(AliHLTUInt64_t)0x1 << index;
	      pCTPClass->SetTitle("CTP Class");
	      pCTPClass->SetName((((TObjString*)entryParams->At(kName))->GetString()).Data());
	      TObjArray* detectors=(((TObjString*)entryParams->At(kDetectors))->GetString()).Tokenize("-");
	      if (detectors) {
		for (int dix=0; dix<detectors->GetEntriesFast(); dix++) {
		  if (!(((TObjString*)detectors->At(dix))->GetString()).IsDigit()) {
		    HLTError("invalid detector list format: trigger class entry %s", entry.Data());
		    break;
		  }
		  // see AliHLTReadoutList::EDetectorId for defines of detectors
		  pCTPClass->Enable(0x1<<(((TObjString*)detectors->At(dix))->GetString()).Atoi());
		}
		delete detectors;
	      }
	    } else {
	    }
	  } else {
	    // the trigger bitfield is fixed to 50 bits (gkNCTPTriggerClasses)
	    HLTError("invalid trigger class entry %s, index width of trigger bitfield exceeded (%d)", entry.Data(), gkNCTPTriggerClasses);
	  }
	} else {
	  HLTError("invalid trigger class entry %s", entry.Data());
	}
	delete entryParams;
      }
    }
    delete classEntries;
  }

  ResetCounters();
  ReadMap();

  return 0;
}

AliHLTUInt64_t AliHLTCTPData::ActiveTriggers(const AliHLTComponentTriggerData& trigData)
{
  // extract active triggers from the trigger data

  const AliRawDataHeader* cdh = NULL;
  if (AliHLTComponent::ExtractTriggerData(trigData, NULL, NULL, &cdh, NULL) != 0) return (AliHLTUInt64_t)0;
  // trigger mask is 50 bit wide and is stored in word 5 and 6 of the CDH
  AliHLTUInt64_t triggerMask = cdh->GetTriggerClasses();
  return triggerMask;
}

bool AliHLTCTPData::EvaluateCTPTriggerClass(const char* expression, const AliHLTComponentTriggerData& trigData) const
{
  // see header file for function documentation
  
  const AliRawDataHeader* cdh = NULL;
  if (AliHLTComponent::ExtractTriggerData(trigData, NULL, NULL, &cdh, NULL, true) != 0) return false;
  // trigger mask is 50 bit wide and is stored in word 5 and 6 of the CDH
  AliHLTUInt64_t triggerMask = cdh->GetTriggerClasses();

  if (fMask!=0 && (triggerMask & fMask)==0) {
    AliHLTEventTriggerData* evtData=reinterpret_cast<AliHLTEventTriggerData*>(trigData.fData);
    HLTWarning("invalid trigger mask 0x%llx, unknown CTP trigger, initialized 0x%llx", triggerMask, fMask);
    for (int i=0; i<8; i++) HLTWarning("\t CDH[%d]=0x%lx", i, evtData->fCommonHeader[i]);
    return false;
  }

  return EvaluateCTPTriggerClass(expression, triggerMask);
}

bool AliHLTCTPData::EvaluateCTPTriggerClass(const char* expression, AliHLTUInt64_t triggerMask) const
{
  // see header file for function documentation

  // use a TFormula to interprete the expression
  // all classname are replaced by '[n]' which means the n'th parameter in the formula
  // the parameters are set to 0 or 1 depending on the bit in the trigger mask
  const vector<unsigned> *pMap=&fMap;
  vector<unsigned> tmp;
  if (fMap.size()==0 && fClassIds.GetLast()>=0) {
    // read map into temporary array and use it
    ReadMap(tmp);
    pMap=&tmp;
    static bool suppressWarning=false;
    if (!suppressWarning) HLTWarning("map not yet initialized, creating local map (slow), suppressing further warnings");
    suppressWarning=true;
  }
  vector<Double_t> par;
  TString condition=expression;
  for (unsigned index=0; index<pMap->size(); index++) {
    const char* className=Name((*pMap)[index]);
    if (className && strlen(className)>0) {
      //HLTDebug("checking trigger class %s", className.Data());
      if (condition.Contains(className)) {
	TString replace; replace.Form("[%d]", par.size());
	//HLTDebug("replacing %s with %s in \"%s\"", className.Data(), replace.Data(), condition.Data());
	condition.ReplaceAll(className, replace);
	if (triggerMask&((AliHLTUInt64_t)0x1<<(*pMap)[index])) par.push_back(1.0);
	else par.push_back(0.0);
      }
    }
  }

  TFormula form("trigger expression", condition);
  if (form.Compile()!=0) {
    HLTError("invalid expression %s", expression);
    return false;
  }
  if (form.EvalPar(&par[0], &par[0])>0.5) return true;
  return false;
}

void AliHLTCTPData::ResetCounters()
{
  // see header file for function documentation
  fCounters.Set(gkNCTPTriggerClasses);
  fCounters.Reset();
}

int AliHLTCTPData::Index(const char* name) const
{
  // see header file for function documentation
  TObject* obj=fClassIds.FindObject(name);
  return obj!=NULL?fClassIds.IndexOf(obj):-1;
}

int AliHLTCTPData::CheckTrigger(const char* name) const
{
  // check status of a trigger class
  int index=Index(name);
  if (index<0) return index;
  return (fTriggers&(0x1<<index))>0?1:0;
}

void AliHLTCTPData::Increment(const char* classIds)
{
  // see header file for function documentation
  TString string=classIds;
  TObjArray* classEntries=string.Tokenize(",");
  if (classEntries) {
    for (int i=0; i<classEntries->GetEntriesFast(); i++) {
      int index=Index(((TObjString*)classEntries->At(i))->GetString().Data());
      if (index>=0 && index<fCounters.GetSize()) fCounters[index]++;
    }
    delete classEntries;
  }
}

void AliHLTCTPData::Increment(AliHLTUInt64_t triggerPattern)
{
  // see header file for function documentation
  AliHLTUInt64_t pattern=triggerPattern&fMask;
  for (int i=0; i<fCounters.GetSize(); i++) {
    if ((pattern&((AliHLTUInt64_t)0x1<<i))==0) continue;
    fCounters[i]++;    
  }
}

void AliHLTCTPData::Increment(int classIdx)
{
  // see header file for function documentation
  if (classIdx<fCounters.GetSize() &&
      (fMask&((AliHLTUInt64_t)0x1<<classIdx))) {
    fCounters[classIdx]++;    
  }
  
}

int AliHLTCTPData::Increment(AliHLTComponentTriggerData& trigData)
{
  // see header file for function documentation
  const AliRawDataHeader* cdh = NULL;
  int result = AliHLTComponent::ExtractTriggerData(trigData, NULL, NULL, &cdh, NULL, true);
  if (result != 0) return result;
  // trigger mask is 50 bit wide and is stored in word 5 and 6 of the CDH
  AliHLTUInt64_t triggerMask = cdh->GetTriggerClasses();

  if (fMask!=0 && (triggerMask & fMask)==0) {
    AliHLTEventTriggerData* evtData=reinterpret_cast<AliHLTEventTriggerData*>(trigData.fData);
    HLTWarning("invalid trigger mask 0x%llx, unknown CTP trigger, initialized 0x%llx", triggerMask, fMask);
    for (int i=0; i<8; i++) HLTWarning("\t CDH[%d]=0x%lx", i, evtData->fCommonHeader[i]);
  }
  Increment(triggerMask);
  return 0;
}

AliHLTUInt64_t AliHLTCTPData::Counter(int index) const
{
  // see header file for function documentation
  if (index>=0 && index<Counters().GetSize()) return Counters()[index];
  return 0;
}

AliHLTUInt64_t AliHLTCTPData::Counter(const char* classId) const
{
  // see header file for function documentation
  return Counter(Index(classId));
}

const char* AliHLTCTPData::Name(int index) const
{
  // see header file for function documentation
  if (index>fClassIds.GetLast()) return NULL;
  return fClassIds.At(index)->GetName();
}

int AliHLTCTPData::ReadMap(vector<unsigned> &map) const
{
  // read the index map for class names
  // for nested class names (e.g. 'myclass' is contained in
  // 'myclassA') the longer names is added first to the map.
  for (int index=0; index<=fClassIds.GetLast(); index++) {
    vector<unsigned>::iterator element=map.begin();
    for (; element!=map.end(); element++) {
      TString name=Name(index);
      if (name.Contains(Name(*element))) {
	// current name contains another one already in the map
	// -> add before and go to next entry
	element=map.insert(element, index);
	break;
      }
    }

    if (element==map.end()) {
      // unique class name, append to map
      map.push_back(index);
    }
  }
  return 0;
}


AliHLTReadoutList AliHLTCTPData::ReadoutList(const AliHLTComponentTriggerData& trigData) const
{
  // see header file for function documentation

  const AliRawDataHeader* cdh = NULL;
  if (AliHLTComponent::ExtractTriggerData(trigData, NULL, NULL, &cdh, NULL, true) != 0) return AliHLTReadoutList();
  // trigger mask is 50 bit wide and is stored in word 5 and 6 of the CDH
  AliHLTUInt64_t triggerMask = cdh->GetTriggerClasses();

  if (fMask!=0 && (triggerMask & fMask)==0) {
    AliHLTEventTriggerData* evtData=reinterpret_cast<AliHLTEventTriggerData*>(trigData.fData);
    HLTWarning("invalid trigger mask 0x%llx, unknown CTP trigger, initialized 0x%llx", triggerMask, fMask);
    for (int i=0; i<8; i++) HLTWarning("\t CDH[%d]=0x%lx", i, evtData->fCommonHeader[i]);
  }

  return ReadoutList(triggerMask);
}

AliHLTReadoutList AliHLTCTPData::ReadoutList(AliHLTUInt64_t  triggerMask) const
{
  // take an 'OR' of all active trigger classes 
  AliHLTReadoutList list;
  for (int i=0; i<gkNCTPTriggerClasses; i++) {
    if (i>fClassIds.GetLast()) break;
    if ((triggerMask&((AliHLTUInt64_t)0x1<<i))==0) continue;
    AliHLTReadoutList* tcrl=(AliHLTReadoutList*)fClassIds.At(i);
    list.OrEq(*tcrl);
  }

  return list;
}


void AliHLTCTPData::Print(Option_t* /*option*/) const
{
  // see header file for function documentation
  cout << GetTitle() << endl;
  cout << "\tactive trigger mask: 0x" << hex << fTriggers << dec << endl;
  int count=0;
  for (int i=0; i<gkNCTPTriggerClasses; i++) {
    if (i>=Counters().GetSize()) break;
    if (i>fClassIds.GetLast()) break;
    if ((fMask&((AliHLTUInt64_t)0x1<<i))==0) continue;
    count++;
    cout << "\t" << i << "\t" << Name(i) << "\t" << Counter(i) << endl;
  }
  if (count==0) cout << "\t(none)" << endl;
}
