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
//  Principles:
//  Snapshots of the system information are writen to the text log files.
//  Text files were chosen in order to get the ouptu also in case code 
//  is crashing. 
//  Following information is stored in the log file:
//  TTimeStamp stamp;
//  CpuInfo_t  cpuInfo;
//  MemInfo_t  memInfo;
//  ProcInfo_t procInfo;
// 
//  Root TSystem is used to retrieve this information:
//  gSystem->GetCpuInfo(&cpuInfo, 10);
//  gSystem->GetMemInfo(&memInfo);
//  gSystem->GetProcInfo(&procInfo);
//  for details see:
//  http://root.cern.ch/root/html/TUnixSystem.html
//  http://root.cern.ch/root/html/ProcInfo_t.html
//  http://root.cern.ch/root/html/MemInfo_t.html
//  http://root.cern.ch/root/html/CpuInfo_t.html
//  -------------------------------------------------------------------
//  class CpuInfo_t
//   Float_t	fIdle	cpu idle percentage
//   Float_t	fLoad15m	cpu load average over 15 m
//   Float_t	fLoad1m	cpu load average over 1 m
//   Float_t	fLoad5m	cpu load average over 5 m
//   Float_t	fSys	cpu sys load in percentage
//   Float_t	fTotal	cpu user+sys load in percentage
//   Float_t	fUser	cpu user load in percentage

//  -------------------------------------------------------------------
//  class ProcInfo_t:
//  Float_t	fCpuSys	system time used by this process in seconds
//  Float_t	fCpuUser	user time used by this process in seconds
//  Long_t	fMemResident	resident memory used by this process in KB
//  Long_t	fMemVirtual	virtual memory used by this process in KB
//  -------------------------------------------------------------------


//  The information from the AliSysInfo can be used as measurement
//  of the code quality. Be aware of the limitation induced by
//  using only system info described in the  AliSysInfo::Test() function
// 
//  The example usage of the AliSysInfo is shown in the
//  AliSysInfo::Test() example.
//  
//  
//
//  
//  The typical usage in the AliRoot code:
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
//#include "AliLog.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"

#include "TTimeStamp.h"
#include "AliSysInfo.h"
#include "TBufferFile.h"
#include "TTreePlayer.h"

//#include "TMemStatManager.h"  //USE IFDEF


using std::endl;
using std::cout;
using std::ios_base;
using std::setprecision;
ClassImp(AliSysInfo)

AliSysInfo* AliSysInfo::fInstance=0;
Bool_t AliSysInfo::fgVerbose = kTRUE;
Bool_t AliSysInfo::fgDisabled = kFALSE;

AliSysInfo::AliSysInfo():
    TObject(),
    fSysWatch(0),
    fTimer(0),
    fMemStat(0),
    fCallBackFunc(0),
    fNCallBack(0)
{
  if (fgDisabled) return;
  fTimer = new TStopwatch;
  fSysWatch  = new fstream("syswatch.log", ios_base::out|ios_base::trunc);

  //hname/C:sname/C:sec/D:mI.fMemUsed/F:mI.fSwapUsed/F:pI.fMemResident/F:pI.fMemVirtual/F:cI.fUser/F:cI.fSys/F:cI.fCpuUser/F:pI.fCpuSys/F


  (*fSysWatch) <<"hname"<<"/C:"               // hname - hostname  
               <<"sname"<<"/C:"              // stamp name
               <<"id0"<<"/I:"                // 0 id
               <<"id1"<<"/I:"                // 1 id
               <<"id2"<<"/I:"                // 1 id
               <<"id3"<<"/I:"                // 1 id
               <<"first"<<"/D:"              // first stamp
    //
	       <<"stampSec"<<"/D:"         // time  - time stamp in seconds
	       <<"mi.fMemUsed"<<"/D:"       // system info 
	       <<"mi.fSwapUsed"<<"/D:"      //
	       <<"cI.fUser"<<"/D:"         //
	       <<"cI.fSys"<<"/D:"         //
	       <<"cI.fLoad1m"<<"/D:"         //
	       <<"cI.fLoad5m"<<"/D:"         //
	       <<"cI.fLoad15m"<<"/D:"         //
    // 
	       <<"pI.fMemResident"<<"/D:"  // process info
	       <<"pI.fMemVirtual"<<"/D:"   //    
	       <<"pI.fCpuUser"<<"/D:"      //
	       <<"pI.fCpuSys"<<"/D:"       //
    //
    	       <<"stampOldSec"<<"/D:"         // time  - time stamp in seconds
	       <<"miOld.fMemUsed"<<"/D:"       // system info - previous
	       <<"miOld.fSwapUsed"<<"/D:"      //
	       <<"cIOld.fUser"<<"/D:"         //
	       <<"cIOld.fSys"<<"/D:"         //
    // 
	       <<"pIOld.fMemResident"<<"/D:"  // process info -previous
	       <<"pIOld.fMemVirtual"<<"/D:"   //    
	       <<"pIOld.fCpuUser"<<"/D:"      //
	       <<"pIOld.fCpuSys"<<"/D:"       //
    // 
	       <<"fileBytesRead"<<"/D:"       // file IO information
	       <<"fileBytesWritten"<<"/D:"    //    
	       <<"fileCounter"<<"/D:"         //
	       <<"fileReadCalls"<<"/D"        //
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


void AliSysInfo::AddStamp(const char *sname, Int_t id0, Int_t id1, Int_t id2, Int_t id3){
  //
  // 
  if (!fgVerbose) return;
  //
  //
  TTimeStamp stamp;
  CpuInfo_t  cpuInfo;
  MemInfo_t  memInfo;
  ProcInfo_t procInfo;  
  gSystem->GetCpuInfo(&cpuInfo, 10);
  gSystem->GetMemInfo(&memInfo);
  gSystem->GetProcInfo(&procInfo);
  //  procInfo.fMemVirtual/=1024;  //size in MBy
  //procInfo.fMemResident/=1024;  //size in MBy

  const char * hname = gSystem->HostName();

  static Int_t entry=0;
  static Double_t  first=stamp.GetSec()+stamp.GetNanoSec()/1000000000.;
  //
  static TTimeStamp stampOld;
  static CpuInfo_t  cpuInfoOld;
  static MemInfo_t  memInfoOld;
  static ProcInfo_t procInfoOld;  
  Double_t fileBytesRead    = TFile::GetFileBytesRead();
  Double_t fileBytesWritten = TFile::GetFileBytesWritten();
  Double_t fileCounter      = TFile::GetFileCounter();
  Double_t fileReadCalls    = TFile::GetFileReadCalls();


  (*(Instance()->fSysWatch)) 
    << hname   <<"\t"               // hname - hostname  
    << sname    <<"\t"              // stamp name
    << id0      <<"\t"
    << id1      <<"\t"
    << id2      <<"\t"
    << id3      <<"\t"
    <<setprecision(15)<< first    <<"\t"              // first stamp               
    //
    <<setprecision(15)<< stamp.GetSec()+stamp.GetNanoSec()/1000000000.<<"\t"         // time  - time stamp in seconds
    << memInfo.fMemUsed<<"\t"       // system info 
    << memInfo.fSwapUsed<<"\t"      //
    << cpuInfo.fUser <<"\t"         //
    << cpuInfo.fSys  <<"\t"         //
    << cpuInfo.fLoad1m  <<"\t"         //
    << cpuInfo.fLoad5m  <<"\t"         //
    << cpuInfo.fLoad15m  <<"\t"         //
    // 
    <<setprecision(15)<< procInfo.fMemResident/1024.<<"\t"  // process info
    <<setprecision(15)<< procInfo.fMemVirtual/1024.<<"\t"   //    
    << procInfo.fCpuUser<<"\t"      //
    << procInfo.fCpuSys<<"\t"       //
    //
    <<setprecision(15)<< stampOld.GetSec()+stampOld.GetNanoSec()/1000000000.<<"\t"         // time  - time stamp in seconds
    << memInfoOld.fMemUsed<<"\t"       // system info - previous
    << memInfoOld.fSwapUsed<<"\t"      //
    << cpuInfoOld.fUser <<"\t"         //
    << cpuInfoOld.fSys  <<"\t"         //
    // 
    <<setprecision(15)<< procInfoOld.fMemResident/1024.<<"\t"  // process info -previous
    <<setprecision(15)<< procInfoOld.fMemVirtual/1024.<<"\t"   //    
    << procInfoOld.fCpuUser<<"\t"      //
    << procInfoOld.fCpuSys<<"\t"       //
    //
    <<fileBytesRead<<"\t"           // file IO information
    <<fileBytesWritten<<"\t"        //    
    <<fileCounter<<"\t"             //
    <<fileReadCalls<<"\t"            //
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


TTree * AliSysInfo::MakeTree(const char *lname, const char * fout){
  // char * lname = "syswatch.log"
  TTree * tree = new TTree;
  tree->ReadFile(lname);
  tree->SetAlias("deltaT","stampSec-stampOldSec");
  tree->SetAlias("T","stampSec-first");
  tree->SetAlias("deltaVM","(pI.fMemVirtual-pIOld.fMemVirtual)");
  tree->SetAlias("VM","pI.fMemVirtual");
  tree->SetAlias("deltaRM","(pI.fMemResident-pIOld.fMemResident)");
  tree->SetAlias("RM","pI.fMemResident");
  if (fout!=0){
    TFile * f = TFile::Open(fout,"recreate");
    f->cd();
    tree->Write("AliSysInfo");
    delete f;
  }
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



TTree*  AliSysInfo::Test(){
  //
  // Test example for AliSysInfo:
  // 1. Make huge memory leak
  // 2. Slow down execution
  /*
    To use the test:
    TTree * tree = AliSysInfo::Test();
    // Make alias what we set as input
    tree->SetAlias("deltaVMIn","(id0*100000+id1*10000+id2*1000)/1000000.")
    tree->SetAlias("deltaVM","(pI.fMemVirtual-pIOld.fMemVirtual)");
    tree->SetAlias("deltaTIn","(id1+id0*10)");
    tree->SetAlias("deltaT","stampSec-stampOldSec");
    //
    tree->SetMarkerStyle(23); tree->SetMarkerSize(0.5);
    // Memory usage
    tree->Draw("deltaVM:deltaVMIn","Entry$>0");
    // or alternative
    tree->Draw("deltaVM:deltaVMIn","Entry$>0","prof");
    //
    // draw time usage
    tree->Draw("deltaT:deltaTIn","Entry$>0"); 
  */
  //
  // The delta of VM as obtained from the AliSysInfo starts to be proportional
  // to  the input allocation after 0.12 MBy (and it is system dependent) 
  // Bellow these limit the deltaVM can be used only in mean.
  // (compare first and  profile histogram)
  for (Int_t id0=0; id0<5; id0++)
    for (Int_t id1=1; id1<10; id1++)
      for (Int_t id2=0; id2<20; id2++){
	new Char_t[id2*1000+id1*10000+id0*100000];  // emulate memory leak
	gSystem->Sleep(id1+id0*10);         // emulate CPU usage 
	AliSysInfo::AddStamp("Leak",id0,id1,id2);
      }
  TTree * tree = AliSysInfo::MakeTree("syswatch.log");
  return tree;  
}

Double_t AliSysInfo::EstimateObjectSize(TObject* object){
  //
  // Estimate size of object as represented in the memory size in bytes
  // Warnings:
  //  1. Only data memebrs which are persistent are counted
  //  2. Local copy of the object is temporary created in memory
  //  3. Do not use it in standard programs, time and memory consument procedure
  //
  if (!object) return 0;
  TBufferFile * file = new TBufferFile(TBuffer::kWrite);
  file->WriteObject(object);
  Double_t size=file->Length();
  delete file;
  return size;
}


void  AliSysInfo::PrintJiraTable(TTree * tree, const char *var, const char *cut, const char *format, const char *outputTable){
  //
  // Print Tree query table as a JIRA table
  // Generic implementation using "custom format" in the future, e.g to print html tables
  // 
  /*
    Example:
    AliSysInfo::PrintJiraTable(tree, "deltaT:sname", "deltaT>0", "col=10:100", "syswatch.table");
    AliSysInfo::PrintJiraTable(tree, "deltaT:sname", "deltaT>0", "col=10:100", 0);
    AliSysInfo::PrintJiraTable(tree, "deltaT:sname:fileReadCalls:fileBytesRead/1000000", "deltaT>0", "col=10:100:25:15", 0);
  */
  tree->SetScanField(tree->GetEntries());
  ((TTreePlayer*)(tree->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(tree->GetPlayer()))->SetScanFileName("AliSysInfo.PrintJiraTable.tmp0");
  tree->Scan(var,cut,format);
  ((TTreePlayer*)(tree->GetPlayer()))->SetScanRedirect(0);
  gSystem->Exec("cat AliSysInfo.PrintJiraTable.tmp0 | sed s_*_\"|\"_g | grep -v \"|||\" > AliSysInfo.PrintJiraTable.tmp1");
  if (outputTable==0){
    gSystem->Exec("cat AliSysInfo.PrintJiraTable.tmp1");
  }else{
    gSystem->Exec(TString::Format("mv AliSysInfo.PrintJiraTable.tmp1 %s", outputTable).Data());
  }

}


TTree * AliSysInfo::MakeDUTree(const char *lname, const char * fout){
  //
  // Provide formatted du - "disk usage"  tree 
  // Input du.txt files assumed to be in 2 collumn format  
  // 1.) process du.txt - Input du.txt files assumed to be in 2 collumn format   - data volume (kBy) dirname
  //     a.) do sorting
  //     b.) append deep
  // 2.) Return tree with defined aliases
  //
  // Example usage:
  //  TTree * tree  =  AliSysInfo::MakeDUTree("du.txt","du.tree");
  /*
    // 1.) Find hotspot
    tree->Scan("sizeTB:dir","depth==3","col=10:70")
    **************************************************************************************************
    *    Row   *     sizeTB *                                                                    dir *
    **************************************************************************************************
    *        3 *    21.1578 *                                   ./SpaceChargeDistortion/data/ATO-108 *
    *        8 *    11.4985 *                     ./reconstruction/dataProductionPreparation/ATO-240 *
    *        9 *    11.4672 *                                                      ./QA/ATO-102/data *
    *       11 *    5.73582 *                                    ./reconstruction/distortionFit/data *
    *       12 *    1.30158 *                ./reconstruction/dataProductionPreparation/ALIROOT-6252 *
    *       15 *    0.74605 *                                                    ./JIRA/ATO-238/data *
    *       16 *   0.474617 *                                                       ./QA/ATO-102/sim *
    // 2. Investigate hotspots  individually 
    tree->Scan("sizeTB:depth:dir","strstr(dir,\"ATO-108\")","col=10:5:70")
    **********************************************************************************************************
    *    Row   *     sizeTB * depth *                                                                    dir *
    **********************************************************************************************************
    *        3 *    21.1578 *     3 *                                   ./SpaceChargeDistortion/data/ATO-108 *
    *        4 *    21.0026 *     4 *                             ./SpaceChargeDistortion/data/ATO-108/alice *
    *        5 *    20.9637 *     5 *                        ./SpaceChargeDistortion/data/ATO-108/alice/data *
    *        6 *    20.9637 *     6 *                   ./SpaceChargeDistortion/data/ATO-108/alice/data/2015 *
    *       15 *    6.62019 *     7 *   ./SpaceChargeDistortion/data/ATO-108/alice/data/2015/LHC15o020212015 *
    *       16 *    6.37346 *     7 *        ./SpaceChargeDistortion/data/ATO-108/alice/data/2015/LHC15o1002 *
    *       25 *    3.61322 *     7 *    ./SpaceChargeDistortion/data/ATO-108/alice/data/2015/LHC15o30012015 *
    *       35 *    1.42999 *     7 *        ./SpaceChargeDistortion/data/ATO-108/alice/data/2015/LHC15o2701 *
    
    or dump JIRA table
    AliSysInfo::PrintJiraTable(tree,"sizeGB:dir","depth==3","col=10:50");
  */

  if (lname==NULL) lname="du.txt";
  if (fout==NULL) fout="du.tree";
  gSystem->GetFromPipe(TString::Format("cat %s | gawk '{print $1\" \"$2\" \"(split($2,a,\"/\")-1)}' | sort -g -r   > %s",lname,fout).Data());
  //
  TTree * tree = new TTree("du","du");
  tree->ReadFile(fout,"size/D:dir/C:depth/d");
  tree->SetAlias("sizeMB","size/(1024)");
  tree->SetAlias("sizeGB","size/(1024*1024)");
  tree->SetAlias("sizeTB","size/(1024*1024*1024)");
  return tree;  
}
