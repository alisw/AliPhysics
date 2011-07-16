#ifndef ALISYSINFO_H
#define ALISYSINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
// This is the class which is to be used during the writing of
// simulated raw data (DDL files format).
// It is using the root functionality in order to deal correctly
// with little/big endian issue. By convention the detector raw
// data payload is stored always with little endian (this corresponds
// to the real life situation when the detector data is coming from
// the hardware).
//-------------------------------------------------------------------------

#include <TObject.h>
class TStopwatch;
class TTree;
class TMemStatManager;

class AliSysInfo : public TObject {
public:
  AliSysInfo();
  static AliSysInfo * Instance();
  static void AddStamp(const char *sname, Int_t id0=-1, Int_t id1=-1, Int_t id2=-1);
  static TTree * MakeTree(const char *lname);
  static void OpenMemStat();
  static void CloseMemStat();
  static Bool_t Contain(const char * str1, const char * str2);
  typedef void (*StampCallback_t)(const Char_t * desription);
  static  void AddCallBack(StampCallback_t callback);
  //
  // Object size function
  static Double_t EstimateObjectSize(TObject* object);
  static  TTree* Test();
private:
  AliSysInfo(const AliSysInfo& source);
  AliSysInfo& operator= (const AliSysInfo& rec);

  fstream         *fSysWatch;       // system watch - Memory and CPU usage 
  TStopwatch      *fTimer;          // timer
  TMemStatManager *fMemStat;      
  static AliSysInfo *   fInstance; //instance pointer
  StampCallback_t *fCallBackFunc; // call back functions
  Int_t           fNCallBack;        // number of call back functions
  ClassDef(AliSysInfo,0)
};

#endif
