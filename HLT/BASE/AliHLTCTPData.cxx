// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
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

/** @file   AliHLTCTPData.cxx
    @author Matthias Richter
    @date   2009-08-20
    @brief  Container for CTP trigger classes and counters
*/

#include "AliHLTCTPData.h"
#include "AliHLTReadoutList.h"
#include "TClass.h"
#include "TObjString.h"
#include "TFormula.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTCTPData)

AliHLTCTPData::AliHLTCTPData()
  : TNamed("AliHLTCTPData", "HLT counters for the CTP")
  , AliHLTLogging()
  , fMask(0)
  , fClassIds(AliHLTReadoutList::Class(), gkNCTPTriggerClasses)
  , fCounters()
{
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
  , fClassIds(AliHLTReadoutList::Class(), gkNCTPTriggerClasses)
  , fCounters()
{
  // see header file for class documentation
  InitCTPTriggerClasses(parameter);
}

AliHLTCTPData::~AliHLTCTPData()
{
  // see header file for class documentation
  fClassIds.Delete();
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
    for (int i=0; i<classEntries->GetEntries(); i++) {
      TString entry=((TObjString*)classEntries->At(i))->GetString();
      TObjArray* entryParams=entry.Tokenize(":");
      if (entryParams) {
	if (entryParams->GetEntries()==3 &&
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

  return 0;
}

bool AliHLTCTPData::EvaluateCTPTriggerClass(const char* expression, AliHLTComponentTriggerData& trigData) const
{
  // see header file for function documentation
  if (trigData.fDataSize != sizeof(AliHLTEventTriggerData)) {
    HLTError("invalid trigger data size: %d expected %d", trigData.fDataSize, sizeof(AliHLTEventTriggerData));
    return false;
  }

  // trigger mask is 50 bit wide and is stored in word 5 and 6 of the CDH
  AliHLTEventTriggerData* evtData=reinterpret_cast<AliHLTEventTriggerData*>(trigData.fData);
  AliHLTUInt64_t triggerMask=evtData->fCommonHeader[6]&0x3ffff;
  triggerMask<<=32;
  triggerMask|=evtData->fCommonHeader[5];

  if (fMask!=0 && (triggerMask & fMask)==0) {
    HLTWarning("invalid trigger mask 0x%llx, unknown CTP trigger, initialized 0x%llx", triggerMask, fMask);
    for (int i=0; i<8; i++) HLTWarning("\t CDH[%d]=0x%lx", i, evtData->fCommonHeader[i]);
    return false;
  }

  // use a TFormula to interprete the expression
  // all classname are replaced by '[n]' which means the n'th parameter in the formula
  // the parameters are set to 0 or 1 depending on the bit in the trigger mask
  //
  // TODO: this will most likely fail for class names like 'base', 'baseA', 'baseB'
  // the class names must be fully unique, none must be contained as substring in
  // another class name. Probably not needed for the moment but needs to be extended.
  vector<Double_t> par;
  TString condition=expression;
  for (int i=0; i<gkNCTPTriggerClasses; i++) {
    const char* className=Name(i);
    if (className && strlen(className)>0) {
      //HLTDebug("checking trigger class %s", className.Data());
      if (condition.Contains(className)) {
	TString replace; replace.Form("[%d]", par.size());
	//HLTDebug("replacing %s with %s in \"%s\"", className.Data(), replace.Data(), condition.Data());
	condition.ReplaceAll(className, replace);
	if (triggerMask&((AliHLTUInt64_t)0x1<<i)) par.push_back(1.0);
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

void AliHLTCTPData::Increment(const char* classIds)
{
  // see header file for function documentation
  TString string=classIds;
  TObjArray* classEntries=string.Tokenize(",");
  if (classEntries) {
    for (int i=0; i<classEntries->GetEntries(); i++) {
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
  if (trigData.fDataSize != sizeof(AliHLTEventTriggerData)) {
    HLTError("invalid trigger data size: %d expected %d", trigData.fDataSize, sizeof(AliHLTEventTriggerData));
    return -EBADF;
  }

  // trigger mask is 50 bit wide and is stored in word 5 and 6 of the CDH
  AliHLTEventTriggerData* evtData=reinterpret_cast<AliHLTEventTriggerData*>(trigData.fData);
  AliHLTUInt64_t triggerMask=evtData->fCommonHeader[6]&0x3ffff;
  triggerMask<<=32;
  triggerMask|=evtData->fCommonHeader[5];

  if (fMask!=0 && (triggerMask & fMask)==0) {
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

AliHLTEventDDL AliHLTCTPData::ReadoutList(const AliHLTComponentTriggerData& trigData) const
{
  // see header file for function documentation
  if (trigData.fDataSize != sizeof(AliHLTEventTriggerData)) {
    HLTError("invalid trigger data size: %d expected %d", trigData.fDataSize, sizeof(AliHLTEventTriggerData));
    AliHLTEventDDL dummy;
    memset(&dummy, 0, sizeof(AliHLTEventDDL));
    return dummy;
  }

  // trigger mask is 50 bit wide and is stored in word 5 and 6 of the CDH
  AliHLTEventTriggerData* evtData=reinterpret_cast<AliHLTEventTriggerData*>(trigData.fData);
  AliHLTUInt64_t triggerMask=evtData->fCommonHeader[6]&0x3ffff;
  triggerMask<<=32;
  triggerMask|=evtData->fCommonHeader[5];

  if (fMask!=0 && (triggerMask & fMask)==0) {
    HLTWarning("invalid trigger mask 0x%llx, unknown CTP trigger, initialized 0x%llx", triggerMask, fMask);
    for (int i=0; i<8; i++) HLTWarning("\t CDH[%d]=0x%lx", i, evtData->fCommonHeader[i]);
  }

  // take an 'OR' of all active trigger classes 
  AliHLTReadoutList list;
  for (int i=0; i<gkNCTPTriggerClasses; i++) {
    if (i>fClassIds.GetLast()) break;
    if ((triggerMask&((AliHLTUInt64_t)0x1<<i))==0) continue;
    AliHLTReadoutList* tcrl=(AliHLTReadoutList*)fClassIds.At(i);
    // 2009-08-27: this is a temorary bugfix:
    // the operator functions of the AliHLTReadoutList class did not work
    // when running online on the HLT cluster. The fix for the moment is
    // to send out the readout list only for the first matching trigger
    // class instead of merging the list. This is sufficient for the
    // current trigger setups but needs to be corrected
    return *tcrl;
    list|=*tcrl;
  }

  return list;
}

void AliHLTCTPData::Print(Option_t* /*option*/) const
{
  // see header file for function documentation
  cout << "CTP counters:" << endl;
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
