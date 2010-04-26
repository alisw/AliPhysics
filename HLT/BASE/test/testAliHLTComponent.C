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

// @file   testAliHLTComponent.C
// @author Matthias Richter
// @date   
// @brief  Test program for the AliHLTComponent base class
// 

#ifndef __CINT__
#include "TDatime.h"
#include "TRandom.h"
#include "AliHLTDataTypes.h"
#include "AliHLTProcessor.h"
#include "AliHLTPluginBase.h"
#include "AliHLTSystem.h"
#include "AliHLTMisc.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTConfiguration.h"
#include "algorithm"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#endif

using namespace std;

class AliHLTTestComponent : public AliHLTProcessor
{
public:
  AliHLTTestComponent() :
    fArguments(),
    fCurrentArgument(fArguments.begin())
  {}

  ~AliHLTTestComponent() {}

  const char* GetComponentID() {return "TestComponent";};
  void GetInputDataTypes(AliHLTComponentDataTypeList& list) {list.clear();}
  AliHLTComponentDataType GetOutputDataType() {return kAliHLTAnyDataType;}
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier) {constBase=0; inputMultiplier=0;}
  AliHLTComponent* Spawn() {return new AliHLTTestComponent;}

  class AliHLTConfigurationArgument {
  public:
    AliHLTConfigurationArgument(const char* argument) :
      fElements(NULL)
    {
      TString work=argument;
      fElements=work.Tokenize(" ");
    }

    AliHLTConfigurationArgument(const AliHLTConfigurationArgument& src) :
      fElements(NULL)
    {
      if (src.fElements) {
	fElements=dynamic_cast<TObjArray*>(src.fElements->Clone());
      }
    }

    ~AliHLTConfigurationArgument() {
      if (fElements) delete fElements;
    }

    AliHLTConfigurationArgument& operator=(const AliHLTConfigurationArgument& src) {
      delete fElements;
      if (src.fElements) {
	fElements=dynamic_cast<TObjArray*>(src.fElements->Clone());
      }
      return *this;
    }

    const char* Argument() {
      if (!fElements) return NULL;
      return ((TObjString*)fElements->At(0))->GetString().Data();
    }

    int NofParameters() {
      if (!fElements) return 0;
      return fElements->GetEntriesFast()-1;
    }

    const char* Parameter(int i) {
      if (!fElements ||
	  fElements->GetEntriesFast()<=i+1) return NULL;
      return ((TObjString*)fElements->At(i+1))->GetString().Data();      
    }

    bool operator==(const char* string) {
      if (!fElements) return 0;
      return (((TObjString*)fElements->At(0))->GetString().CompareTo(string))==0;
    }

    bool operator!=(const char* string) {
      return !operator==(string);
    }

    void Print() {
      if (!fElements) {
	cout << "#############  empty ############" << endl;
	return;
      }
      cout << "   Print: " << Argument() << " with " << NofParameters() << " parameter(s)";
      for (int n=0; n<NofParameters(); n++) cout << " " << Parameter(n);
      cout << endl;
    }

  private:
    TObjArray* fElements;
  };

  int ScanConfigurationArgument(int argc, const char** argv) {
    if (fCurrentArgument==fArguments.end()) return 0;
    int count=0;

    // check whether it is an argument at all
    if (*(argv[count])!='-') {
      cerr << "not a recognized argument: " << argv[count] << endl;
      return -EINVAL;
    }

    // check whether the argument matches
    //fCurrentArgument->Print();
    if (*fCurrentArgument!=argv[count]) {
      cerr << "argument sequence does not match: got " << argv[count] << "  expected " << fCurrentArgument->Argument() << endl;
      return -EINVAL;
    }

    count++;
    // read parameters
    if (fCurrentArgument->NofParameters()>0) {
      if (argc<=count) {
	cerr << "missing parameter" << endl;
	return -EPROTO;
      }

      // implement more checks here
      count+=fCurrentArgument->NofParameters();
    }
    fCurrentArgument++;
    return count;
  }

  int FillArgumentVector(const char** arguments, vector<AliHLTConfigurationArgument>& list) {
    list.clear();
    for (const char** iter=arguments; *iter!=NULL; iter++) {
      list.push_back(AliHLTConfigurationArgument(*iter));
    }
    return list.size();
  }

  int FillArgv(vector<AliHLTConfigurationArgument>& list, vector<const char*>& argv) {
    argv.clear();
    for (vector<AliHLTConfigurationArgument>::iterator argument=list.begin();
	 argument!=list.end(); argument++) {
      argv.push_back(argument->Argument());
      for (int n=0; n<argument->NofParameters(); n++) {
	argv.push_back(argument->Parameter(n));
      }
    }
    return argv.size();
  }

  void PrintArgv(int argc, const char** argv) {
    for (int n=0; n<argc; n++) {
      cout << "   " << n << " : " << argv[n] << endl;
    }
  }

  int CheckSequence(const char* sequence[], int mode=0) {
    int iResult=0;
    if ((iResult=FillArgumentVector(sequence, fArguments))<0) {
      cerr << "failed to fill argument vector" << endl;
      return iResult;
    }
    vector<const char*> argv;
    if (mode==0) {
      if ((iResult=FillArgv(fArguments, argv))<0) {
	cerr << "failed to fill argument array" << endl;
	return iResult;
      }
    } else {
      for (const char** element=sequence; *element!=NULL; element++)
	argv.push_back(*element);
    }
    fCurrentArgument=fArguments.begin();
    //PrintArgv(argv.size(), &argv[0]);
    if ((iResult=ConfigureFromArgumentString(argv.size(), &argv[0]))<0) {
      cerr << "ConfigureFromArgumentString failed " << endl;
      return iResult;
    }

    return iResult;
  }

  int CheckConfigure() {
    int iResult=0;
    const char* sequence1[]={"-sequence1","-argument2 5", NULL};
    if ((iResult=CheckSequence(sequence1))<0) {
      cerr << "failed checking sequence " << sequence1[0] << endl;
      return iResult;
    }

    const char* sequence2[]={"-sequence2","-argument2 5 8", "-argument3 test", NULL};
    if ((iResult=CheckSequence(sequence2))<0) {
      cerr << "failed checking sequence in mode 0: " << sequence2[0] << endl;
      return iResult;
    }

    if ((iResult=CheckSequence(sequence2, 1))<0) {
      cerr << "failed checking sequence in mode 1: " << sequence2[0] << endl;
      return iResult;
    }

    const char* sequence3[]={"-solenoidBz 5", NULL};
    if ((iResult=CheckSequence(sequence3))<0) {
      cerr << "failed checking sequence " << sequence3[0] << endl;
      return iResult;
    }
    return iResult;
  }

  int InitCTPTest(const char* param) {
    // this quick test needs to be the functions of the base class to be
    // defined 'protected'
    SetupCTPData();
    //return ScanECSParam(param);
    //return InitCTPTriggerClasses(param);
    return -ENOSYS;
  }

  bool CheckCTP(const char* expression, AliHLTComponentTriggerData* data) {
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
    SetupCTPData();
    return 0;
  }

  int DoDeinit() {
    return 0;
  }

  int DoEvent( const AliHLTComponentEventData& /*evtData*/,
	       AliHLTComponentTriggerData& /*trigData*/) {
    if (!IsDataEvent()) return 0;

    cout << "DoEvent: run no " << GetRunNo() << endl;
    if ((int)GetRunNo()!=AliHLTMisc::Instance().GetCDBRunNo()) {
      return -1;
    }
    return 0;
  }
private:
  vector<AliHLTConfigurationArgument> fArguments;
  vector<AliHLTConfigurationArgument>::iterator fCurrentArgument;
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

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// setup of the CTP test

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

/**
 * Test the CTP trigger tools
 * The base class is initialized with an ECS string of randomly defined trigger
 * classes consisting of random bits and names. Than a couple of expressions is
 * test for various random bit patterns.
 */
int testCTPTrigger()
{
  cout << "checking CTP functionality of the base class" << endl;
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
    if (element->Bit()<10) parameter+="0";
    parameter+=element->Bit(); 
    parameter+=":";
    parameter+=element->ClassName(); parameter+=":";
    parameter+="05-01-06"; // just a test pattern for the detector ids, ignored for the moment
    element++;
  }
  parameter+=";HLT_MODE=A;RUN_NO=0";

  AliHLTTestComponent component;
  component.SetGlobalLoggingLevel(kHLTLogDefault);
  if ((iResult=component.InitCTPTest(parameter.Data()))<0) {
    cerr << "InitCTPTest failed :" << iResult << endl;
    return iResult;
  }
  cout << "init ECS parameter: " << parameter << endl;

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
	trigger=component.CheckCTP(expression.Data(), trigData.Data());
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
	trigger=component.CheckCTP(expression.Data(), trigData.Data());
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
      trigger=component.CheckCTP(expression.Data(), trigData.Data());
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

      trigger=component.CheckCTP(expression.Data(), trigData.Data());
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

      trigger=component.CheckCTP(expression.Data(), trigData.Data());
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

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// setup of the Configure test
int testConfigure() 
{
  cout << "checking common configuration tools of the base class" << endl;
  AliHLTTestComponent component;
  return component.CheckConfigure();
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// test event precessing sequence
int testEventProcessing()
{
  int iResult=0;
  cout << "checking AliHLTComponent processing sequence" << endl;
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  if (!pHLT) {
    cerr << "error: failed to create HLTSystem" << endl;
    return -1;
  }

  // add the test component to the handler
  pHLT->fpComponentHandler->AddComponent(new AliHLTTestComponent);

  // configurattion with just the TestComponent
  AliHLTConfiguration testcomponent("testcomponent", "TestComponent", "", "");
  pHLT->BuildTaskList("testcomponent");

  // setup the OCDB and init a random runno
  AliHLTMisc::Instance().InitCDB("local:///tmp");
  int runNo=GetRandom(0,1000);
  AliHLTMisc::Instance().SetCDBRunNo(runNo);

  // run one event
  if ((iResult=pHLT->Run(1))<0) {
    cerr << "error: event processing failed" << endl;
    return iResult;
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// main functions

int testAliHLTComponent()
{
  int iResult=0;
  //if ((iResult=testCTPTrigger())<0) return iResult;
  if ((iResult=testConfigure())<0) return iResult;
  if ((iResult=testEventProcessing())<0) return iResult;
  return iResult;
}

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  iResult=testAliHLTComponent();
  return iResult;
}
