// $Id$

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

/** @file   testAliHLTComponent_CTPTrigger.C
    @author Matthias Richter
    @date   
    @brief  Test program for default data types
 */

#ifndef __CINT__
#include "TDatime.h"
#include "TRandom.h"
#include "AliHLTDataTypes.h"
#include "AliHLTProcessor.h"
#endif

class AliHLTTestComponent : public AliHLTProcessor
{
public:
  AliHLTTestComponent() {}
  ~AliHLTTestComponent() {}

  const char* GetComponentID() {return "TestComponent";};
  void GetInputDataTypes(AliHLTComponentDataTypeList& list) {list.clear();}
  AliHLTComponentDataType GetOutputDataType() {return kAliHLTAnyDataType;}
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier) {constBase=0; inputMultiplier=0;}
  AliHLTComponent* Spawn() {return new AliHLTTestComponent;}

  int InitTest(const char* param) {
    return ScanECSParam(param);
    //return InitCTPTriggerClasses(param);
  }

  bool Check(const char* expression, AliHLTComponentTriggerData* data) {
    return EvaluateCTPTriggerClass(expression, *data);
  }

  class AliHLTCTPTriggerClass {
  public:
    AliHLTCTPTriggerClass() : fBit(~(unsigned)0), fClassName(""), fTrigger(false) {}
    ~AliHLTCTPTriggerClass() {}

    bool        Valid() {return fBit!=~(unsigned)0;}
    unsigned    Bit() {return fBit;}
    void        Bit(unsigned bit) {fBit=bit;}
    bool        Trigger() {return fTrigger;}
    void        Trigger(bool trigger) {fTrigger=trigger;}
    const char* ClassName() {return fClassName.c_str();}    
    void        ClassName(const char* classname) {fClassName=classname;}    
  private:
    unsigned fBit;
    string   fClassName;
    bool     fTrigger;
  };
protected:
  int DoInit(int /*argc*/, const char** /*argv*/) {
    return 0;
  }

  int DoDeinit() {
    return 0;
  }

  int DoEvent( const AliHLTComponentEventData& /*evtData*/,
	       AliHLTComponentTriggerData& /*trigData*/) {
    return 0;
  }
};

class AliHLTTriggerDataAccess
{
public:
  AliHLTTriggerDataAccess()
    : fData(NULL)
    , fEventData(NULL)
    , fCDH(NULL)
    , fMine(NULL)
  {
    unsigned size=sizeof(AliHLTComponentTriggerData) + sizeof(AliHLTEventTriggerData);
    fMine=new Byte_t[size];
    memset(fMine, 0, size);
    AliHLTComponentTriggerData* data=reinterpret_cast<AliHLTComponentTriggerData*>(fMine);
    data->fData=fMine+sizeof(AliHLTComponentTriggerData);
    Set(data);
  }

  AliHLTTriggerDataAccess(AliHLTComponentTriggerData* pData)
    : fData(NULL)
    , fEventData(NULL)
    , fCDH(NULL)
    , fMine(NULL)
  {
    if (fMine) delete [] fMine;
    fMine=NULL;
    Set(pData);
  }

  ~AliHLTTriggerDataAccess(){
    if (fMine) delete [] fMine;
    fMine=NULL;
    fData=NULL;
    fEventData=NULL;
    fCDH=NULL;
  }

  AliHLTComponentTriggerData* Data() {return fData;}

  Long64_t              TriggerMask() {
    Long64_t mask=0;
    if (fCDH) {
      mask=fCDH[6];
      mask<<=32;
      mask|=fCDH[5];
    }
    return mask;
  }

  int Set(AliHLTComponentTriggerData* data) {
    fData=data;
    fData->fDataSize=sizeof(AliHLTEventTriggerData);
    fEventData=reinterpret_cast<AliHLTEventTriggerData*>(fData->fData);
    fCDH=fEventData->fCommonHeader;
    return 0;
  }

  int ResetCDH() {
    if (fCDH) {
      memset(fCDH, 0, 32);
    }
    return 0;
  }

  int TriggerBit(unsigned bit, bool set) {
    if ((int)bit>=gkNCTPTriggerClasses) return -EINVAL;
    if (!fCDH) return -ENODATA;

    int word=5;
    if (bit>=32) {
      word++;
      bit-=32;
    }
    if (set)
      fCDH[word]|=(UInt_t)0x1<<bit;
    else
      fCDH[word]&=~((UInt_t)0x1<<bit);
      
    return bit;
  }

private:
  AliHLTTriggerDataAccess(const AliHLTTriggerDataAccess&);
  AliHLTTriggerDataAccess& operator=(const AliHLTTriggerDataAccess&);

  AliHLTComponentTriggerData* fData;
  AliHLTEventTriggerData*     fEventData;
  AliHLTUInt32_t*             fCDH;
  Byte_t*                     fMine;
};

/**
 * Get a random number in the given range.
 */
int GetRandom(int min, int max)
{
  if (max-min<2) return min;
  static TRandom rand;
  static bool seedSet=false;
  if (!seedSet) {
    TDatime dt;
    rand.SetSeed(dt.Get());
    seedSet=true;
  }
  return min+rand.Integer(max-min);
}

/**
 * Generate a random name of given length
 */
string GenerateTriggerClassName(int length)
{
  string tcn;
  for (int i=0; i<length; i++) {
    unsigned char c=GetRandom(48, 83);
    if (c>57) c+=7;
    tcn+=c;
  }
  return tcn;
}

/**
 * Generate an array of trigger classes.
 * The array has the specified size but the number antries which are actually
 * filled is randomized.
 */
int GenerateTriggerClasses(int size, vector<AliHLTTestComponent::AliHLTCTPTriggerClass>& classes)
{
  classes.clear();
  classes.resize(size);
  unsigned count=GetRandom(4, size>16?size/2:size);
  for (unsigned i=0; i<count; i++) {
    int bit=0;
    do {
      bit=GetRandom(0, size);
    } while (classes[bit].Valid());
    classes[bit].Bit(bit);
    classes[bit].ClassName((GenerateTriggerClassName(GetRandom(5,15))).c_str());
  }
  return classes.size();
}

int testAliHLTComponent_CTPTrigger()
{
  int iResult=0;
  vector<AliHLTTestComponent::AliHLTCTPTriggerClass> triggerClasses;
  if (GenerateTriggerClasses(GetRandom(5,gkNCTPTriggerClasses), triggerClasses)<=0) {
    return -1;
  }

  TString parameter="CONFIGURE;CTP_TRIGGER_CLASS=";
  vector<AliHLTTestComponent::AliHLTCTPTriggerClass>::iterator element=triggerClasses.begin();
  while (element!=triggerClasses.end()) {
    if (!element->Valid()) {
      element=triggerClasses.erase(element);
      continue;
    }
    if (!parameter.EndsWith("=")) parameter+=",";
    parameter+=element->Bit(); parameter+=":";
    parameter+=element->ClassName(); parameter+=":";
    parameter+=0;
    element++;
  }
  parameter+=";HLT_MODE=A;RUN_NO=0";

  AliHLTTestComponent component;
  component.SetGlobalLoggingLevel(kHLTLogDefault);
  if ((iResult=component.InitTest(parameter.Data()))<0) return iResult;

  AliHLTTriggerDataAccess trigData;
  for (int cycle=0; cycle<500 && iResult>=0; cycle++) {
    for (element=triggerClasses.begin();
	 element!=triggerClasses.end(); element++) {
      element->Trigger(GetRandom(0,100)>50);
    }

    vector<AliHLTTestComponent::AliHLTCTPTriggerClass> shuffle;
    shuffle.assign(triggerClasses.begin(), triggerClasses.end());
    for (unsigned int trial=0; trial<2*triggerClasses.size() && iResult>=0; trial++) {
      random_shuffle(shuffle.begin(), shuffle.end());

      bool result=0;
      bool trigger=0;
      TString expression;
      trigData.ResetCDH();
      for (element=shuffle.begin();
	   element!=shuffle.end(); element++) {
	trigData.TriggerBit(element->Bit(), element->Trigger());
      }

      // single class
      for (element=shuffle.begin();
	   element!=shuffle.end() && iResult>=0 && trial<3;
	   element++) {
	// is
	result=element->Trigger();
	expression=element->ClassName();
	trigger=component.Check(expression.Data(), trigData.Data());
	if (trigger!=result) {
	  cout << expression << ": " << element->Trigger()
	       << "->" << trigger 
	       << std::hex << "   (" << trigData.TriggerMask() << ")"
	       << endl;
	  cerr << "trigger does not match, expected " << result << endl;
	  iResult=-1;
	  break;
	}

	// is not
	expression="!";
	expression+=element->ClassName();
	result=!result;
	trigger=component.Check(expression.Data(), trigData.Data());
	if (trigger!=result) {
	  cout << expression << ": " << element->Trigger()
	       << "->" << trigger 
	       << std::hex << "   (" << trigData.TriggerMask() << ")"
	       << endl;
	  cerr << "trigger does not match, expected " << result << endl;
	  iResult=-1;
	  break;
	}
      }

      // OR
      result=shuffle[0].Trigger() || shuffle[1].Trigger() || shuffle[2].Trigger();
      expression.Form("%s || %s || %s",
		      shuffle[0].ClassName(), shuffle[1].ClassName(), shuffle[2].ClassName());
      trigger=component.Check(expression.Data(), trigData.Data());
      if (trigger!=result) {
	cout << expression << ": " << shuffle[0].Trigger() << shuffle[1].Trigger() << shuffle[2].Trigger() 
	     << "->" << trigger 
	     << std::hex << "   (" << trigData.TriggerMask() << ")"
	     << endl;
	cerr << "trigger does not match, expected " << result << endl;
	iResult=-1;
	break;
      }

      // AND
      result=shuffle[0].Trigger() && shuffle[1].Trigger() && shuffle[2].Trigger();
      expression.Form("%s && %s && %s",
		      shuffle[0].ClassName(), shuffle[1].ClassName(), shuffle[2].ClassName());

      trigger=component.Check(expression.Data(), trigData.Data());
      if (trigger!=result) {
	cout << expression << ": " << shuffle[0].Trigger() << shuffle[1].Trigger() << shuffle[2].Trigger() 
	     << "->" << trigger 
	     << std::hex << "   (" << trigData.TriggerMask() << ")"
	     << endl;
	cerr << "trigger does not match, expected " << result << endl;
	iResult=-1;
	break;
      }

      // mixed OR/AND
      result=shuffle[0].Trigger() && (shuffle[1].Trigger() || shuffle[2].Trigger());
      expression.Form("%s && (%s || %s)",
		      shuffle[0].ClassName(), shuffle[1].ClassName(), shuffle[2].ClassName());

      trigger=component.Check(expression.Data(), trigData.Data());
      if (trigger!=result) {
	cout << expression << ": " << shuffle[0].Trigger() << shuffle[1].Trigger() << shuffle[2].Trigger() 
	     << "->" << trigger 
	     << std::hex << "   (" << trigData.TriggerMask() << ")"
	     << endl;
	cerr << "trigger does not match, expected " << result << endl;
	iResult=-1;
	break;
      }
    }
  }
  if (iResult<0) {
    cerr << "check failed, dumping info" << endl;
    cerr << "ECS param: " << parameter << endl;
    for (element=triggerClasses.begin();
	 element!=triggerClasses.end(); element++) {
      cerr << element->Trigger() << " " << element->Bit() << ": " << element->ClassName() << endl;
    }
  }
  return iResult;
}

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  iResult=testAliHLTComponent_CTPTrigger();
  return iResult;
}
