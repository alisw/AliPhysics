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

/** @file   testAliHLTCTPData.C
    @author Matthias Richter
    @date   
    @brief  Test program for the AliHLTCTPData class
 */

#ifndef __CINT__
#include "TDatime.h"
#include "TRandom.h"
#include "AliHLTDataTypes.h"
#include "algorithm"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "AliHLTDAQ.h"
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cerrno>
#include "AliHLTCTPData.h"
#include "AliHLTReadoutList.h"
#endif

using namespace std;

class AliHLTCTPDataTest {
public:
  AliHLTCTPDataTest() {}
  ~AliHLTCTPDataTest() {}

  class AliHLTCTPTriggerClass {
  public:
    AliHLTCTPTriggerClass() : fBit(~(unsigned)0), fClassName(""), fTrigger(false), fDetectorParam(), fDetectors() {}
    ~AliHLTCTPTriggerClass() {}

    bool        Valid() {return fBit!=~(unsigned)0;}
    unsigned    Bit() {return fBit;}
    void        Bit(unsigned bit) {fBit=bit;}
    bool        Trigger() {return fTrigger;}
    void        Trigger(bool trigger) {fTrigger=trigger;}
    const char* ClassName() {return fClassName.c_str();}
    void        ClassName(const char* classname) {fClassName=classname;}
    const char* DetectorParam() {
      if (fDetectorParam.IsNull() && fDetectors.size()>0) {
	for (unsigned  i=0; i<fDetectors.size(); i++) {
	  if (!fDetectorParam.IsNull()) fDetectorParam+="-";
	  if (fDetectors[i]<10) fDetectorParam+="0";
	  fDetectorParam+=fDetectors[i];
	}
      }
      return fDetectorParam.Data();
    }
    void        DetectorParam(const char* detectorparam) {fDetectorParam=detectorparam;}
    int         AddDetector(unsigned short detector) {
      fDetectorParam.Clear();
      for (vector<unsigned short>::iterator element=fDetectors.begin();
	   element!=fDetectors.end(); element++) {
	if (detector<*element) {
	  fDetectors.insert(element, detector);
	  return 0;
	}
      }
      fDetectors.push_back(detector);
      return 0;
    }
    bool        HasDetector(int detector) const {
      for (unsigned i=0; i<fDetectors.size(); i++) {
	if (fDetectors[i]==detector) {
	  return true;
	}
      }
      return false;
    }
    bool        HasDetector(const char* id) const {
      TString detectorid=id;
      for (unsigned i=0; i<fDetectors.size(); i++) {
	if (detectorid.CompareTo(AliHLTDAQ::OnlineName(fDetectors[i]))==0) {
	  return true;
	}
      }
      return false;
    }

  private:
    unsigned fBit;
    string   fClassName;
    bool     fTrigger;
    TString  fDetectorParam;
    vector<unsigned short> fDetectors;
  };

protected:
private:
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
    data->fStructSize=sizeof(AliHLTComponentTriggerData);
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
    fEventData->fCommonHeaderWordCnt=gkAliHLTCommonHeaderCount;
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
  bool toggle=false;
  for (int i=0; i<length; i++) {
    // add random number of '-' characters, but not at beginning or end
    if (i>0 && i<length-1 && GetRandom(5,10)==7 && (toggle=(!toggle))) {
      tcn+='-';
      continue;
    }
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
int GenerateTriggerClasses(int size, vector<AliHLTCTPDataTest::AliHLTCTPTriggerClass>& classes)
{
  classes.clear();
  classes.resize(size);
  unsigned count=GetRandom(4, size>16?size/2:size);
  int last=-1;
  for (unsigned i=0; i<count; i++) {
    int bit=0;
    do {
      bit=GetRandom(0, size);
    } while (classes[bit].Valid());
    classes[bit].Bit(bit);
    if (last>=0 && GetRandom(5,10)==7) {
      // generate similar class names by just adding or truncating characters
      string name=classes[last].ClassName();
      if (GetRandom(0,100)%2) {
	name+=GenerateTriggerClassName(GetRandom(3,7)).c_str();
	cout << "using appended name " << name << " (" << classes[last].ClassName() << ")" << endl;
      } else {
	int num=GetRandom(1,name.length()/2);
	name.erase(name.length()-num, num);
	cout << "using truncated name " << name << " (" << classes[last].ClassName() << ")" << endl;
      }
      classes[bit].ClassName(name.c_str());
    } else {
    classes[bit].ClassName((GenerateTriggerClassName(GetRandom(5,15))).c_str());
    }
    last=bit;
    unsigned nofdetectors=GetRandom(1, 10);
    unsigned short detector=17;
    for (unsigned k=0; k<nofdetectors; k++) {
      detector=GetRandom(nofdetectors-k-1, detector-1);
      classes[bit].AddDetector(detector);
    }
  }
  return classes.size();
}

/**
 * Test the CTP trigger tools
 * The base class is initialized with an ECS string of randomly defined trigger
 * classes consisting of random bits and names. Than a couple of expressions is
 * test for various random bit patterns.
 */
int testAliHLTCTPData()
{
  cout << "checking AliHLTCTPData class" << endl;
  int iResult=0;
  vector<AliHLTCTPDataTest::AliHLTCTPTriggerClass> triggerClasses;
  if (GenerateTriggerClasses(GetRandom(5,gkNCTPTriggerClasses), triggerClasses)<=0) {
    return -1;
  }

  TString ecs_parameter="CONFIGURE";
  TString ctp_parameter="CTP_TRIGGER_CLASS=";
  vector<AliHLTCTPDataTest::AliHLTCTPTriggerClass>::iterator element=triggerClasses.begin();
  while (element!=triggerClasses.end()) {
    if (!element->Valid()) {
      element=triggerClasses.erase(element);
      continue;
    }
    if (!ctp_parameter.EndsWith("=")) ctp_parameter+=",";
    if (element->Bit()<10) ctp_parameter+="0";
    ctp_parameter+=element->Bit(); 
    ctp_parameter+=":";
    ctp_parameter+=element->ClassName(); ctp_parameter+=":";
    ctp_parameter+=element->DetectorParam();
    element++;
  }

  vector<const char*> parameters;
  parameters.push_back("HLT_MODE=B");
  parameters.push_back("RUN_NO=0");
  parameters.push_back(ctp_parameter.Data());
  random_shuffle(parameters.begin(), parameters.end());
  for (vector<const char*>::iterator par=parameters.begin();
       par!=parameters.end(); par++) {
    ecs_parameter+=";";
    ecs_parameter+=*par;
  }

  AliHLTCTPData ctpdata;
  ctpdata.SetGlobalLoggingLevel(kHLTLogDefault);
  if ((iResult=ctpdata.InitCTPTriggerClasses(ctp_parameter.Data()))<0) {
    cerr << "InitCTPTriggerClasses failed :" << iResult << endl;
    return iResult;
  }

  for (element=triggerClasses.begin();
       element!=triggerClasses.end();
       element++) {
    int index=ctpdata.Index(element->ClassName());
    if (index<0) {
      cerr << "error: can not find CTP trigger class name " << element->ClassName() << endl;
      return -1;
    }
    if (element->Bit()!=(unsigned)index) {
      cerr << "error: wrong index for CTP trigger class " << element->ClassName() << ": expected " << element->Bit() << " - got " << index << endl;
      return -1;      
    }
  }

  AliHLTTriggerDataAccess trigData;

  // check the readout lists for the different trigger classes
  for (element=triggerClasses.begin();
       element!=triggerClasses.end(); element++) {
    trigData.ResetCDH();
    trigData.TriggerBit(element->Bit(), 1);
    AliHLTReadoutList detectorReadout(ctpdata.ReadoutList(*trigData.Data()));
    for (int detectorid=0; detectorid<17; detectorid++) {
	if ((detectorReadout.DetectorEnabled(0x1<<detectorid) && !element->HasDetector(detectorid)) ||
	    (!detectorReadout.DetectorEnabled(0x1<<detectorid) && element->HasDetector(detectorid))) {
	  cerr << "readout list does not match trigger class " << element->Bit() << ":" << element->ClassName() << ":" << element->DetectorParam() << endl;
	  detectorReadout.Print();
	  return -1;
	}
     }
  }

  const int nofCycles=500;
  for (int cycle=0; cycle<nofCycles && iResult>=0; cycle++) {
    bool bHaveTrigger=false;
    for (element=triggerClasses.begin();
	 element!=triggerClasses.end(); element++) {
      element->Trigger(GetRandom(0,100)>50);
      bHaveTrigger|=element->Trigger();
    }
    if (!bHaveTrigger) {
      // trigger at least one class
      (triggerClasses.begin())->Trigger(1);
    }

    vector<AliHLTCTPDataTest::AliHLTCTPTriggerClass> shuffle;
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
	//cout << " " << element->Bit() << ":" << element->Trigger();
      }
      //cout << endl;
      ctpdata.SetTriggers(*(trigData.Data()));

      // single class
      for (element=shuffle.begin();
	   element!=shuffle.end() && iResult>=0 && trial<3;
	   element++) {
	int check=ctpdata.CheckTrigger(element->ClassName());
	if (check<0) {
	  cerr << "error avaluating CheckTrigger: class name " << element->ClassName() << " not found" << endl;
	  iResult=-1;
	  break;	  
	}
	if (check!=element->Trigger()) {
	  cerr << "error avaluating CheckTrigger, class name " << element->ClassName() << ": expecting " << element->Trigger() << " - got " << check << endl;
	  ctpdata.Print();
	  iResult=-1;
	  break;	  
	}

	// is
	result=element->Trigger();
	expression=element->ClassName();
	trigger=ctpdata.EvaluateCTPTriggerClass(expression.Data(), *trigData.Data());
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
	trigger=ctpdata.EvaluateCTPTriggerClass(expression.Data(), *trigData.Data());
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
      trigger=ctpdata.EvaluateCTPTriggerClass(expression.Data(), *trigData.Data());
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

      trigger=ctpdata.EvaluateCTPTriggerClass(expression.Data(), *trigData.Data());
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

      trigger=ctpdata.EvaluateCTPTriggerClass(expression.Data(), *trigData.Data());
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
    // readout list
    AliHLTReadoutList detectorReadout(ctpdata.ReadoutList(*trigData.Data()));
    for (int detectorid=0; detectorid<17 && iResult>=0; detectorid++) {
	if (detectorReadout.DetectorEnabled(0x1<<detectorid)) {
	// detector is included in the readout, find at least one active trigger
	// class with this detector
	for (element=triggerClasses.begin();
	     element!=triggerClasses.end(); element++) {
	  if (element->Trigger() && element->HasDetector(detectorid)) break;
	}
	if (element==triggerClasses.end()) {
	  cerr << "can not find any active trigger class for detector " << detectorid << " enabled in the readout" << endl;
	  iResult=-1;
	}
      } else {
	// check that this detector is not part of any of the active trigger classes
	for (element=triggerClasses.begin();
	     element!=triggerClasses.end(); element++) {
	  if (element->Trigger() && element->HasDetector(detectorid)) {
	    cerr << "detector " << detectorid << " not enabled in the readout but enabled in active trigger class " << element->ClassName() << endl;
	    iResult=-1;
	  }
	}
      }
      if (iResult<0) {
	detectorReadout.Print();
      }
    }
    if ((cycle+1)%(nofCycles/5)==0) cout << " " << (100*(cycle+1))/nofCycles << " % done" << endl;
  }

  if (iResult<0) {
    cerr << "check failed, dumping info" << endl;
    cerr << "CTP param: " << ctp_parameter << endl;
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
// main functions

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  iResult=testAliHLTCTPData();
  return iResult;
}
