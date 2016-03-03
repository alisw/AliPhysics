#ifndef ALISYSINFO_H
#define ALISYSINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <TObject.h>
class TStopwatch;
class TTree;
class TMemStatManager;
using std::fstream;

class AliSysInfo : public TObject {
public:
  AliSysInfo();
  static AliSysInfo * Instance();
  static void AddStamp(const char *sname, Int_t id0=-1, Int_t id1=-1, Int_t id2=-1, Int_t id3=-1);
  static TTree * MakeTree(const char *lname, const char * fout=0);
  static void OpenMemStat();
  static void CloseMemStat();
  static Bool_t Contain(const char * str1, const char * str2);
  typedef void (*StampCallback_t)(const Char_t * desription);
  static  void AddCallBack(StampCallback_t callback);
  //
  // Object size function
  static Double_t EstimateObjectSize(TObject* object);
  static  TTree* Test();
  
  static void SetVerbose(Bool_t v=kFALSE)   {fgVerbose = v;}
  static Bool_t GetVerbose()                {return fgVerbose;} 
  static void SetDisabled(Bool_t v=kTRUE)   {fgDisabled = v; fgVerbose=kFALSE;}
  static Bool_t IsDisabled()               {return fgDisabled;}

private:
  AliSysInfo(const AliSysInfo& source);
  AliSysInfo& operator= (const AliSysInfo& rec);

  fstream         *fSysWatch;       // system watch - Memory and CPU usage 
  TStopwatch      *fTimer;          // timer
  TMemStatManager *fMemStat;      
  static AliSysInfo *   fInstance; //instance pointer
  StampCallback_t *fCallBackFunc; // call back functions
  Int_t           fNCallBack;        // number of call back functions
  static Bool_t   fgVerbose;      // do we want actually to write the stamps ?
  static Bool_t   fgDisabled;      // do not even open the file
  ClassDef(AliSysInfo,0)
};

#endif
