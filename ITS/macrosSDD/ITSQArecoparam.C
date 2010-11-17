#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TAlienFile.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <Riostream.h>
#include <TObjArray.h>
#include <TKey.h>
#include <TStyle.h>
#include <TTree.h>
#include <TROOT.h>
#include <TPluginManager.h>
#include <TClass.h>
#include <iostream>
#include <stdio.h>
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"
#include "AliITSQADataMakerRec.h"
#include "AliITSQASDDDataMakerRec.h"
#include "AliITSQAChecker.h"
#include "AliQAChecker.h"
#include "AliITSQASDDChecker.h"
#include "AliReconstructor.h"
#include "AliCDBManager.h"
#include "AliQAv1.h"
#include "AliGeomManager.h"
#include "AliITSInitGeometry.h"
#include "AliITSgeom.h"
#include "AliRecoParam.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliRecoParam.h"
#include "AliDetectorRecoParam.h"
#include "AliITSReconstructor.h"
#include "AliITSRecPointContainer.h"
#include "AliLog.h"
#endif
void ITSQArecoparam(char *iFile,Int_t idet=2,Int_t FirstEvt=0, Int_t MaxEvts=1000000)
{
  TString namefile(iFile);
  if(namefile.Contains("alien"))
    {
      TGrid::Connect("alien://"); 
      if(!gGrid) {
	printf("gGrid not found! exit macro\n");
	return;
      }
    }
  gStyle->SetPalette(1);
  Int_t ic=0;
  AliRawReader *rd=NULL; 
  if(strstr(iFile,".root")!=0){rd = new AliRawReaderRoot(iFile,FirstEvt);}
  else{rd=new AliRawReaderDate(iFile,FirstEvt);}
  Int_t runNumber = rd->GetRunNumber();
  cout << "ITS Quality Assurance Prototype" << endl; 
  //TStopwatch mytimer;
  //TString namefile(iFile);
  // Set OCDB if needed
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    printf("Setting a local default storage and run number \n");
    if(namefile.Contains("alien")){
      man->SetDefaultStorage("raw://");
    }
    else{man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");}
    man->SetRun(runNumber);
    AliQAv1::SetQARefStorage("local://$ALICE_ROOT/QARef") ;
  }

  AliITSQADataMakerRec *itsQAdm = new AliITSQADataMakerRec(kTRUE,idet,0); //online kTRUE
  itsQAdm->SetWriteExpert() ;
  itsQAdm->SetRunNumber(runNumber);  
  //________________________For the RecPoints____________________________________
  /************************************************/
  TPluginManager* pluginManager=NULL;
  TPluginHandler* pluginHandler=NULL; 
  AliReconstructor* reconstructor = NULL;
  AliITSRecPointContainer* rpcont=NULL;
  AliGeomManager::LoadGeometry("geometry.root");
  AliGeomManager::ApplyAlignObjsFromCDB("ITS");
  //   ITS initializations
  
  AliITSInitGeometry initgeom;
  AliITSgeom *geom = initgeom.CreateAliITSgeom();
  printf("Geometry name: %s\n",(initgeom.GetGeometryName()).Data());
  
  printf("Loading reconstruction parameter objects for detector ITS\n");
  AliRecoParam fRecoParam; 
  AliCDBPath path("ITS","Calib","RecoParam");
  AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
  Bool_t cacheStatus = AliCDBManager::Instance()->GetCacheFlag();
  if(!entry){printf("Couldn't find RecoParam entry in OCDB for detector ITS");entry=NULL;}
  else {
    TObject *recoParamObj = entry->GetObject();
    if (dynamic_cast<TObjArray*>(recoParamObj)) {
      printf("RecoParam TObjArray\n");
      fRecoParam.AddDetRecoParamArray(0,dynamic_cast<TObjArray*>(recoParamObj));
    }
    else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
      printf("RecoParam AliDetectorRecoParam\n");
      printf("Single set of reconstruction parameters found for detector ITS");
      (dynamic_cast<AliDetectorRecoParam*>(recoParamObj))->SetAsDefault();
      fRecoParam.AddDetRecoParam(0,(dynamic_cast<AliDetectorRecoParam*>(recoParamObj)));
    }
    else {printf("Error: No valid RecoParam object found in the OCDB for detector ITS");}
    entry->SetOwner(0);
  }
  if(!cacheStatus)entry->SetObject(NULL);
  if(!cacheStatus){ delete entry;}
  
  // load the reconstructor object
 pluginManager = gROOT->GetPluginManager();
  TString detName = "ITS";
  TString recName = "Ali" + detName + "Reconstructor";  
  
pluginHandler = pluginManager->FindHandler("AliReconstructor", "ITS");
  // if not, add a plugin for it
  if (!pluginHandler) {
    printf("defining plugin for ITS\n");
    TString libs = gSystem->GetLibraries();
    if (libs.Contains("lib" + detName + "base.so") ||
	(gSystem->Load("lib" + detName + "base.so") >= 0)) {pluginManager->AddHandler("AliReconstructor", detName,recName, detName + "rec", recName + "()");}
    else {pluginManager->AddHandler("AliReconstructor", detName,recName, detName, recName + "()");}
    pluginHandler = pluginManager->FindHandler("AliReconstructor", detName);
  }
  if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {reconstructor = (AliReconstructor*) pluginHandler->ExecPlugin(0);}
  if (fRecoParam.GetDetRecoParamArray(0) && !AliReconstructor::GetRecoParam(0)) {
    const AliDetectorRecoParam *par = fRecoParam.GetDetRecoParam(0);
    reconstructor->Init();
    reconstructor->SetRecoParam(par);
  }
  
  /*AliITSRecPointContainer**/ rpcont=AliITSRecPointContainer::Instance();
  rpcont->PrepareToRead();

  Int_t cycleLength = 5;  
  //cout << "Processing Run " << runNumber << endl;
  cout << "Init: " << AliQAv1::kRAWS << endl;
  
  TObjArray **objArray= itsQAdm->Init(AliQAv1::kRAWS, cycleLength);
  cout<<"raw tobjarray :"<<objArray<<"\n"<<endl;
  for(Int_t spec = 0 ; spec < 5 ; spec++){
    if(spec==1){
      itsQAdm->SetEventSpecie(AliRecoParam::kLowMult);
      AliQAv1::Instance()->SetEventSpecie(AliRecoParam::kLowMult);
      itsQAdm->InitRaws();
    }
    else{continue;}
  }
  itsQAdm->StartOfCycle(AliQAv1::kRAWS,runNumber,kFALSE);
  /*********************************************************************/

    cout << "Init: " << AliQAv1::kRECPOINTS << endl;
    TObjArray **objArray1=  itsQAdm->Init(AliQAv1::kRECPOINTS, cycleLength);
    cout<<"recpoint tobjarray :"<<objArray1<<"\n"<<endl;
    for(Int_t spec = 0 ; spec < 5 ; spec++){
      if(spec==1){
	itsQAdm->SetEventSpecie(AliRecoParam::kLowMult);
	AliQAv1::Instance()->SetEventSpecie(AliRecoParam::kLowMult);
	itsQAdm->InitRecPoints();
      }
      else{continue;}
    }
    itsQAdm->StartOfCycle(AliQAv1::kRECPOINTS,runNumber,kTRUE);

  /*********************************************************************/
  Int_t iev = 0;
  while(rd->NextEvent() && iev < MaxEvts ) {
    cout<<">>>>>>>   Processing event number: "<<++iev<<endl;
    /*******************************************************/
    if(itsQAdm->IsCycleDone()) {
      cout << "end of cycle" << endl;
      AliQAChecker::Instance()->SetRunNumber(AliCDBManager::Instance()->GetRun());	
      itsQAdm->EndOfCycle(AliQAv1::kRAWS);
      itsQAdm->StartOfCycle(AliQAv1::kRAWS,ic++,kFALSE);
    }
    /******************************************************/
    /*************************************************/

      if(itsQAdm->IsCycleDone()) {
	cout << "end of cycle" << endl;
	AliQAChecker::Instance()->SetRunNumber(AliCDBManager::Instance()->GetRun());
	itsQAdm->EndOfCycle(AliQAv1::kRECPOINTS);
	itsQAdm->StartOfCycle(AliQAv1::kRECPOINTS,ic,kTRUE);
      } 
 
    /*************************************************/
    cout<<"Beginning Exec"<<endl;
    cout<<"AliQAv1::kRAWS   "<<AliQAv1::kRAWS<<endl;
    itsQAdm->SetEventSpecie(AliRecoParam::kLowMult);
    AliQAv1::Instance()->SetEventSpecie(AliRecoParam::kLowMult);
    itsQAdm->Exec(AliQAv1::kRAWS,rd);
    /***************************************************/
    //  if(kRecPoints){
      cout<<"AliQAv1::kRECPOINTS   "<<AliQAv1::kRECPOINTS<<endl;
      cout << "DigitsToRecPoints" << endl;
      TTree* fTreeR = new TTree("TreeR", "Reconstructed Points Container"); //make a tree
      Char_t option[5];

      if(idet==0)sprintf(option,"ALL");
      else if(idet==1)sprintf(option,"SPD");
      else if(idet==2)sprintf(option,"SDD");
      else if(idet==3)sprintf(option,"SSD");
      printf("\t\t===========>option is %s\n",option);

      rpcont->PrepareToRead();
      reconstructor->Reconstruct(rd,fTreeR);

      itsQAdm->SetEventSpecie(AliRecoParam::kLowMult);
      AliQAv1::Instance()->SetEventSpecie(AliRecoParam::kLowMult)  ; 
      itsQAdm->Exec(AliQAv1::kRECPOINTS,fTreeR);
      cout<<"Finishing Exec"<<endl;
    
      ((AliITSReconstructor*)reconstructor)->ResetRecPoints();
      delete fTreeR; 
      fTreeR=NULL; 
      /*****************************************************/
    }

  cout << "end RAWS cycle: " << AliQAv1::kRAWS << endl;
  cout << "refStorage: " << AliQAv1::GetQARefStorage() << endl;
  cout << "end of cycle 2" << endl;
  AliQAChecker::Instance()->SetRunNumber(AliCDBManager::Instance()->GetRun());
  itsQAdm->EndOfCycle(AliQAv1::kRAWS);   
  cout << "Raws QA completed for " << iev << " events" << endl;
  /*******************************************************************/

    AliQAChecker::Instance()->SetRunNumber(AliCDBManager::Instance()->GetRun());
    itsQAdm->EndOfCycle(AliQAv1::kRECPOINTS);
    cout << "RecPoints QA completed for " << iev << " events" << endl;

  /*******************************************************************/
  itsQAdm->Finish(); // write to the output File

  cout << "Call AliITSQASDDDataMakerRec destructor" << endl;
  delete itsQAdm;
  itsQAdm=NULL;

}

