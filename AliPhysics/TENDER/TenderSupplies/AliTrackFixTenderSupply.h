#ifndef ALITRACKFIXTENDERSUPPLY_H
#define ALITRACKFIXTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Apply on-the-fly fix to tracks                                    //
//                                                                    //
//  19/06/2012: RS: Add 1/pt shift from AODB to TPC and TPC-ITS       //
//                  Optionally correct also track coordinate          //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include "AliTenderSupply.h"


class AliESDVertex;
class AliExternalTrackParam;
class AliOADBContainer;
class AliESDtrack;
class AliOADBTrackFix;

class AliTrackFixTenderSupply: public AliTenderSupply {
  
public:
  //
  AliTrackFixTenderSupply();
  AliTrackFixTenderSupply(const char *name, const AliTender *tender=NULL);
  virtual ~AliTrackFixTenderSupply();
  virtual  void ProcessEvent();
  virtual  void Init() {}
  //
  Double_t GetSideAFraction(const AliESDtrack* track) const;
  void     CorrectTrackPtInv(AliExternalTrackParam* trc, int mode, double sideAfraction, double phi) const;
  Bool_t   GetRunCorrections(int run);
  Bool_t   LoadOADBObjects();
  //
  void     SetOADBObjPath(const char* path)        { fOADBObjPath = path; }
  void     SetOADBObjName(const char* name)        { fOADBObjName = name; }
  TString& GetOADBObjPath()                 const  { return (TString&)fOADBObjPath; }
  TString& GetOADBObjName()                 const  { return (TString&)fOADBObjName; }
  //
  void     SetDebugLevel(Int_t l=1)                {fDebug = l;}
  Int_t    GetDebugLevel()                  const  {return fDebug;}
  //
private:
  
  AliTrackFixTenderSupply(const AliTrackFixTenderSupply&c);
  AliTrackFixTenderSupply& operator= (const AliTrackFixTenderSupply&c);
  //
  Int_t             fDebug;                  // Debug level
  Double_t          fBz;                     // mag field from ESD
  AliOADBTrackFix*  fParams;                 // parameters for current run
  TString           fOADBObjPath;            // path of file with parameters to use, starting from OADB dir
  TString           fOADBObjName;            // name of the corrections object in the OADB container
  AliOADBContainer* fOADBCont;               // OADB container with parameters collection
  //
  ClassDef(AliTrackFixTenderSupply, 1);  // track fixing tender task 
};


#endif

