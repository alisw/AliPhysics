/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


//
// Origin:  marian.ivanov@cern.ch
//
//  Make   a log file for the CPU and Memory usage
//  
//  The typical usage:
//  Make a set of stamps in the code in the place of interest
//  e.g. 
//
//  AliSysInfo::AddStamp("Start");
//
//  loader->UnloadRecPoints();
//    AliSysInfo::AddStamp(Form("LRec%s_%d",fgkDetectorName[iDet],eventNr), iDet,1,eventNr);
//

// The log file can be transformed to the tree - to make a visualization
// See $ALICE_ROOT/macros/PlotSysInfo.C as an example


#include <Riostream.h>
#include "AliLog.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "TTimeStamp.h"
#include "AliSysInfo.h"

//#include "TMemStatManager.h"  //USE IFDEF


ClassImp(AliSysInfo)

AliSysInfo* AliSysInfo::fInstance=0;

AliSysInfo::AliSysInfo():
    TObject(),
    fSysWatch(0),
    fTimer(0),
    fMemStat(0),
    fCallBackFunc(0),
    fNCallBack(0)
{
  fTimer = new TStopwatch;
  fSysWatch  = new fstream("syswatch.log", ios_base::out|ios_base::trunc);

  //hname/C:sname/C:sec/I:mI.fMemUsed/F:mI.fSwapUsed/F:pI.fMemResident/F:pI.fMemVirtual/F:cI.fUser/F:cI.fSys/F:cI.fCpuUser/F:pI.fCpuSys/F


  (*fSysWatch) <<"hname"<<"/C:"               // hname - hostname  
               <<"sname"<<"/C:"              // stamp name
               <<"id0"<<"/I:"                // 0 id
               <<"id1"<<"/I:"                // 1 id
               <<"id2"<<"/I:"                // 1 id
               <<"first"<<"/I:"              // first stamp
    //
	       <<"stampSec"<<"/I:"         // time  - time stamp in seconds
	       <<"mi.fMemUsed"<<"/F:"       // system info 
	       <<"mi.fSwapUsed"<<"/F:"      //
	       <<"cI.fUser"<<"/F:"         //
	       <<"cI.fSys"<<"/F:"         //
    // 
	       <<"pI.fMemResident"<<"/F:"  // process info
	       <<"pI.fMemVirtual"<<"/F:"   //    
	       <<"pI.fCpuUser"<<"/F:"      //
	       <<"pI.fCpuSys"<<"/F:"       //
    //
    	       <<"stampOldSec"<<"/I:"         // time  - time stamp in seconds
	       <<"miOld.fMemUsed"<<"/F:"       // system info - previous
	       <<"miOld.fSwapUsed"<<"/F:"      //
	       <<"cIOld.fUser"<<"/F:"         //
	       <<"cIOld.fSys"<<"/F:"         //
    // 
	       <<"pIOld.fMemResident"<<"/F:"  // process info -previous
	       <<"pIOld.fMemVirtual"<<"/F:"   //    
	       <<"pIOld.fCpuUser"<<"/F:"      //
	       <<"pIOld.fCpuSys"<<"/F"       //
	       << endl;
  
}




AliSysInfo * AliSysInfo::Instance(){
  //
  //
  //
  if (!fInstance){
    fInstance = new AliSysInfo;
  }
  return fInstance;
}


void AliSysInfo::AddStamp(const char *sname, Int_t id0, Int_t id1, Int_t id2){
  //
  // 
  //
  //
  TTimeStamp stamp;
  CpuInfo_t  cpuInfo;
  MemInfo_t  memInfo;
  ProcInfo_t procInfo;  
  gSystem->GetCpuInfo(&cpuInfo, 10);
  gSystem->GetMemInfo(&memInfo);
  gSystem->GetProcInfo(&procInfo);
  procInfo.fMemVirtual*=0.001;  //size in MBy
  procInfo.fMemResident*=0.001;  //size in MBy

  const char * hname = gSystem->HostName();

  static Int_t entry=0;
  static Int_t first=stamp.GetSec();
  //
  static TTimeStamp stampOld;
  static CpuInfo_t  cpuInfoOld;
  static MemInfo_t  memInfoOld;
  static ProcInfo_t procInfoOld;  


  (*(Instance()->fSysWatch)) << hname   <<"\t"               // hname - hostname  
               << sname    <<"\t"              // stamp name
               << id0      <<"\t"
               << id1      <<"\t"
               << id2      <<"\t"
               << first    <<"\t"              // first stamp               
               //
	       << stamp.GetSec()<<"\t"         // time  - time stamp in seconds
	       << memInfo.fMemUsed<<"\t"       // system info 
	       << memInfo.fSwapUsed<<"\t"      //
	       << cpuInfo.fUser <<"\t"         //
	       << cpuInfo.fSys  <<"\t"         //
               // 
	       << procInfo.fMemResident<<"\t"  // process info
	       << procInfo.fMemVirtual<<"\t"   //    
	       << procInfo.fCpuUser<<"\t"      //
	       << procInfo.fCpuSys<<"\t"       //
    //
    	       << stampOld.GetSec()<<"\t"         // time  - time stamp in seconds
	       << memInfoOld.fMemUsed<<"\t"       // system info - previous
	       << memInfoOld.fSwapUsed<<"\t"      //
	       << cpuInfoOld.fUser <<"\t"         //
	       << cpuInfoOld.fSys  <<"\t"         //
               // 
	       << procInfoOld.fMemResident<<"\t"  // process info -previous
	       << procInfoOld.fMemVirtual<<"\t"   //    
	       << procInfoOld.fCpuUser<<"\t"      //
	       << procInfoOld.fCpuSys<<"\t"       //
	       << endl;
  stampOld   = stamp;
  cpuInfoOld = cpuInfo;
  memInfoOld = memInfo;
  procInfoOld= procInfo;

  //  if (fInstance->fMemStat) fInstance->fMemStat->AddStamps(sname);
  for (Int_t icallback=0; icallback<Instance()->fNCallBack; icallback++){
    Instance()->fCallBackFunc[icallback](sname);
  }
  entry++;
}


TTree * AliSysInfo::MakeTree(const char *lname){
  // char * lname = "syswatch.log"
  TTree * tree = new TTree;
  tree->ReadFile(lname);
  tree->SetAlias("deltaT","stampSec-stampOldSec");
  tree->SetAlias("T","stampSec-first");
  tree->SetAlias("deltaVM","(pI.fMemVirtual-pIOld.fMemVirtual)");
  tree->SetAlias("VM","pI.fMemVirtual");
  return tree;
}


Bool_t AliSysInfo::Contain(const char * str1, const char * str2){
  //
  //
  //
  TString str(str1);
  return str.Contains(str2);
}



void AliSysInfo::OpenMemStat(){
  //
  //
  //
  //USE IFDEF if MEMSTAT ENABLED  
  //  Instance()->fMemStat = TMemStatManager::GetInstance();
  //   Instance()->fMemStat->SetAutoStamp(10000000, 10000000,1000000);
  //   Instance()->fMemStat->Enable();  
}

void AliSysInfo::CloseMemStat(){
  //
  //
  //
  //USE IFDEF if MEMSTAT ENABLED
  //if (Instance()->fMemStat  == TMemStatManager::GetInstance()) Instance()->fMemStat->Close();
  //Instance()->fMemStat=0;
}



void AliSysInfo::AddCallBack(StampCallback_t callback){
  //
  // add cal back function
  //
  AliSysInfo *info =  Instance();
  if (!info->fCallBackFunc)
    info->fCallBackFunc = new StampCallback_t[100];
  info->fCallBackFunc[info->fNCallBack]=callback;
  info->fNCallBack++;
}
