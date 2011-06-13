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
#include <TGeoGlobalMagField.h>
#include <TMap.h>
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
#include "AliMagF.h"
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
#include "AliGRPObject.h"
#include "AliRunInfo.h"
#endif
void ITSQArecoparam(char *iFile, Int_t runNb=150000, Int_t idet=2, Int_t FirstEvt=0, Int_t MaxEvts=1000000)
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
  //Int_t runNumber = rd->GetRunNumber();
  cout << "ITS Quality Assurance Prototype for run "<< runNb << endl; 
  //TStopwatch mytimer;
  //TString namefile(iFile);
  // Set OCDB if needed
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    printf("Setting a local default storage and run number \n");
    //  if(namefile.Contains("alien")){
      man->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
      //}
      //else{man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");}
    man->SetRun(runNb);
    AliQAv1::SetQARefStorage("local://$ALICE_ROOT/QARef") ;
  }
  AliCDBManager::Instance()->GetAll(Form("ITS/Calib/*"));
  

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");

  AliGRPObject *fGRPData=NULL;
  
  TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (m) {
       printf("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
       m->Print();
       fGRPData = new AliGRPObject();
       fGRPData->ReadValuesFromMap(m);
    }

    else {
       printf("Found an AliGRPObject in GRP/GRP/Data, reading it");
       fGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0);
    }


 if (!fGRPData) {
     printf("Error   No GRP entry found in OCDB!");
     return;
  }


  TString lhcState = fGRPData->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    printf("Error  GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = fGRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    printf("Error   GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = fGRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    printf("Error   GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }

  TString runType = fGRPData->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    printf("Error  GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  Int_t activeDetectors = fGRPData->GetDetectorMask();
  if (activeDetectors==AliGRPObject::GetInvalidUInt()) {
    printf("Error  GRP/GRP/Data entry:  missing value for the detector mask ! Using 1074790399");
    activeDetectors = 1074790399;
  }

  AliRunInfo *fRunInfo = new AliRunInfo(lhcState, beamType, beamEnergy, runType, activeDetectors);
  fRunInfo->Dump();

 if ( TGeoGlobalMagField::Instance()->IsLocked() ) {
    if (TGeoGlobalMagField::Instance()->GetField()->TestBit(AliMagF::kOverrideGRP)) {
      printf("ExpertMode!!! GRP information will be ignored !");
      printf("ExpertMode!!! Running with the externally locked B field !");
    }
    else {
      printf("Destroying existing B field instance!");
      delete TGeoGlobalMagField::Instance();
    }    
  }
  if ( !TGeoGlobalMagField::Instance()->IsLocked() ) {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = fGRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      printf("Error :  GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }
    
    Char_t l3Polarity = fGRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      printf("Error:   GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

   // Dipole
    Float_t diCurrent = fGRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      printf("Error  GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }

    Char_t diPolarity = fGRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      printf("Error GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }


   Int_t  polConvention = fGRPData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
   Bool_t uniformB = fGRPData->IsUniformBMap();
   
   if (ok) { 
     AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1), 
					    TMath::Abs(diCurrent) * (diPolarity ? -1:1), 
					    polConvention,uniformB,beamEnergy, beamType.Data());
     if (fld) {
       TGeoGlobalMagField::Instance()->SetField( fld );
       TGeoGlobalMagField::Instance()->Lock();
       printf("Running with the B field constructed out of GRP !");
     }
     else printf("Fatal Failed to create a B field map !");
   }
   else printf("Fatal  B field is neither set nor constructed from GRP ! Exitig...");
  }
  



  AliITSQADataMakerRec *itsQAdm = new AliITSQADataMakerRec(kTRUE,idet,0); //online kTRUE
  itsQAdm->SetWriteExpert() ;
  itsQAdm->SetRunNumber(runNb);  
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
  AliCDBEntry *entry2=AliCDBManager::Instance()->Get(path.GetPath());
  Bool_t cacheStatus = AliCDBManager::Instance()->GetCacheFlag();
  if(!entry2){printf("Couldn't find RecoParam entry in OCDB for detector ITS");entry2=NULL;}
  else {
    TObject *recoParamObj = entry2->GetObject();
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
    entry2->SetOwner(0);
  }
  if(!cacheStatus)entry2->SetObject(NULL);
  if(!cacheStatus){ delete entry2;}
  
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
      itsQAdm->SetEventSpecie(AliRecoParam::kCosmic);
      AliQAv1::Instance()->SetEventSpecie(AliRecoParam::kCosmic);
      itsQAdm->InitRaws();
    }
    else{continue;}
  }
  itsQAdm->StartOfCycle(AliQAv1::kRAWS,runNb,kFALSE);
  /*********************************************************************/

    cout << "Init: " << AliQAv1::kRECPOINTS << endl;
    TObjArray **objArray1=  itsQAdm->Init(AliQAv1::kRECPOINTS, cycleLength);
    cout<<"recpoint tobjarray :"<<objArray1<<"\n"<<endl;
    for(Int_t spec = 0 ; spec < 5 ; spec++){
      if(spec==1){
	itsQAdm->SetEventSpecie(AliRecoParam::kCosmic);
	AliQAv1::Instance()->SetEventSpecie(AliRecoParam::kCosmic);
	itsQAdm->InitRecPoints();
      }
      else{continue;}
    }
    itsQAdm->StartOfCycle(AliQAv1::kRECPOINTS,runNb,kTRUE);

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
    itsQAdm->SetEventSpecie(AliRecoParam::kCosmic);
    AliQAv1::Instance()->SetEventSpecie(AliRecoParam::kCosmic);
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

      itsQAdm->SetEventSpecie(AliRecoParam::kCosmic);
      AliQAv1::Instance()->SetEventSpecie(AliRecoParam::kCosmic)  ; 
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

void ITSQArecoparam(Int_t runNb=150000,Int_t year=2011,Char_t period[10]="LHC11a",Int_t idet=2,Int_t FirstEvt=0, Int_t MaxEvts=1000000)
{
  //TString namefile(iFile);
  //if(namefile.Contains("alien"))
  //{
      TGrid::Connect("alien://");
      if(!gGrid) {
        printf("gGrid not found! exit macro\n");
        return;
      }
      //}

  gSystem->AddIncludePath("-I. -I$ALICE_ROOT/include");

  char path2[200];

  sprintf(path2,"/alice/data/");
  printf("path %s\n",path2);
  char rawfilename[200];
  sprintf(rawfilename,"%04i/%s/%09i/raw/%02i%09i*.*.root",year,period,runNb,year-2000,runNb);


  Int_t number=0;

  TGridResult *gr = gGrid->Query(path2,rawfilename);
  if (gr->GetEntries() < 1) {
    printf("In this run there are not raws files: Exit macro\n");
    return;
  }
  //if(gr->GetEntries() > 10) number=3;
  //else number=3;
  Bool_t chunkok=kFALSE;
  TString filenamedef; //= gr->GetKey(number,"turl");
  //printf("-> FILE %s \n",filename);
  while(chunkok==kFALSE || number<gr->GetEntries())
    {
      const char* filename = gr->GetKey(number,"turl");
      printf("-> FILE %s \n",filename);
      TString namefile(filename);
      if(namefile.Contains(".10.root")==kFALSE){chunkok=kTRUE; filenamedef.Form("%s",filename); break;} 
      else{
	chunkok=kFALSE;
	number++;
      }
      if(number==gr->GetEntries()-1) {chunkok=kTRUE; filenamedef.Form("%s",filename); break;}       
    }

  char *filetouse=(char*)filenamedef.Data();//=filenamedef.Data(); //=NULL;
  //sprintf(filetouse,"%s",(char*)filenamedef.Data());

  printf("File to use = %s\n\n",filetouse);

  Int_t runn=runNb;

  printf("Run: %d \n",runn);

  ITSQArecoparam(filetouse,runn,idet,FirstEvt,MaxEvts);

}
