/*
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");

gSystem->Load("libXrdClient.so");
gSystem->Load("libNetx.so");
if (!gGrid) TGrid::Connect("alien://",0,0,"t");


.L CalibEnv.C++
Init();
CalibEnv("listAll.txt");
TFile f("dcsTime.root")


*/

#include <iostream>
#include <fstream>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliLog.h>
#include <AliMagF.h>
#include "AliTPCcalibDB.h"
#include "AliTPCAltroMapping.h"
#include "AliTPCExB.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCSensorTempArray.h"
#include "AliGRPObject.h"
#include "AliTPCTransform.h"
#include "TFile.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalibPulser.h"
#include "AliTPCCalibPedestal.h"
#include "AliTPCCalibCE.h"
#include "AliTPCExBFirst.h"
#include "TTreeStream.h"
#include "AliTPCTempMap.h"

void ProcessGoofie( AliDCSSensorArray* goofieArray, TVectorD & vecEntries, TVectorD & vecMedian, TVectorD &vecMean, TVectorD &vecRMS);

void Init(){
  //
  //
  //
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Parameters","local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("GRP/GRP/Data","local:///lustre_alpha/alice/alien/alice/data/2008/LHC08d/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Temperature","local:///lustre_alpha/alice/alien/alice/data/2008/LHC08d/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/HighVoltage","local:///lustre_alpha/alice/alien/alice/data/2008/LHC08d/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Goofie","local:///lustre_alpha/alice/alien/alice/data/2008/LHC08d/OCDB/");
  AliCDBManager::Instance()->SetRun(1);
}


void InitAlien(const char *path="LHC08b"){
  //
  //
  //
  TString alpath="alien://folder=/alice/data/2008/";
  alpath+=path;
  alpath+="/OCDB";
    
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Parameters","local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("GRP/GRP/Data",alpath.Data());
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Temperature",alpath.Data());
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Goofie",alpath.Data());
  AliCDBManager::Instance()->SetRun(1);
}


void CalibEnv(const char * runList){
  //
  //
  //  
  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  ifstream in;
  in.open(runList);
  Int_t irun=0;
  TTreeSRedirector *pcstream = new TTreeSRedirector("dcsTime.root");
  //  for (Int_t irun=startRun; irun<stopRun; irun++){
  while(in.good()) {
    in >> irun;
    if (irun==0) continue;
    printf("Processing run %d ...\n",irun);
    AliDCSSensor * sensorPressure = AliTPCcalibDB::Instance()->GetPressureSensor(irun);
    if (!sensorPressure) continue;
    AliTPCSensorTempArray * tempArray = AliTPCcalibDB::Instance()->GetTemperatureSensor(irun);
    AliTPCTempMap * tempMap = new AliTPCTempMap(tempArray);
    AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(irun);
    //
    Int_t startTime = sensorPressure->GetStartTime();
    Int_t endTime = sensorPressure->GetEndTime();
    Int_t dtime = TMath::Max((endTime-startTime)/20,10*60);
    for (Int_t itime=startTime; itime<endTime; itime+=dtime){
      //
      TTimeStamp tstamp(itime); 
      Float_t valuePressure  = calibDB->GetPressure(tstamp,irun,0);
      Float_t valuePressure2 = calibDB->GetPressure(tstamp,irun,1);

      TLinearFitter * fitter = 0;
      TVectorD vecTemp[10];
      for (Int_t itype=0; itype<5; itype++)
	for (Int_t iside=0; iside<2; iside++){	  
	  fitter= tempMap->GetLinearFitter(itype,iside,tstamp);	  
	  if (!fitter) continue;
	  fitter->Eval(); 
	  fitter->GetParameters(vecTemp[itype+iside*5]);
	  delete fitter;
	} 
      
      
      TVectorD vecGoofie, vecEntries, vecMean, vecMedian,vecRMS;
      if (goofieArray){	
	vecGoofie.ResizeTo(goofieArray->NumSensors());
	ProcessGoofie(goofieArray, vecEntries ,vecMedian, vecMean, vecRMS);
	//
	for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
	  AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
	  if (gsensor){
	    vecGoofie[isensor] = gsensor->GetValue(tstamp);
	  }
	}
      }
      Double_t ptrelative0 = AliTPCcalibDB::GetPTRelative(tstamp,irun,0);
      Double_t ptrelative1 = AliTPCcalibDB::GetPTRelative(tstamp,irun,1);
      //
      Double_t voltagesIROC[36]; 
      Double_t voltagesOROC[36]; 
      for(Int_t j=1; j<36; j++) voltagesIROC[j-1] = AliTPCcalibDB::Instance()->GetChamberHighVoltage(startTime, irun, j);
      for(Int_t j=36; j<72; j++) voltagesOROC[j-36] = AliTPCcalibDB::Instance()->GetChamberHighVoltage(startTime, irun, j);
      Double_t voltIROC = TMath::Median(36, voltagesIROC);
      Double_t voltOROC = TMath::Median(36, voltagesOROC);

      //tempMap->GetLinearFitter(0,0,itime);
      (*pcstream)<<"dcs"<<
	"run="<<irun<<
	"time="<<itime<<
	"voltageIROC="<<voltIROC<<
	"voltageOROC="<<voltOROC<<
	"ptrel0="<<ptrelative0<<
	"ptrel1="<<ptrelative1<<
	"goofie.="<<&vecGoofie<<
	"goofieE.="<<&vecEntries<<
	"goofieMean.="<<&vecMean<<
	"goofieMedian.="<<&vecMedian<<
	"goofieRMS.="<<&vecRMS<<
	"press="<<valuePressure<<
	"press2="<<valuePressure2<<
	"temp00.="<<&vecTemp[0]<<
	"temp10.="<<&vecTemp[1]<<
	"temp20.="<<&vecTemp[2]<<
	"temp30.="<<&vecTemp[3]<<
	"temp40.="<<&vecTemp[4]<<
	"temp01.="<<&vecTemp[5]<<
	"temp11.="<<&vecTemp[6]<<
	"temp21.="<<&vecTemp[7]<<
	"temp31.="<<&vecTemp[8]<<
	"temp41.="<<&vecTemp[9]<<
	"\n";
    }
  }
  delete pcstream;
}


void ProcessGoofie( AliDCSSensorArray* goofieArray, TVectorD & vecEntries, TVectorD & vecMedian, TVectorD &vecMean, TVectorD &vecRMS){
  /*
    
  1       TPC_ANODE_I_A00_STAT
  2       TPC_DVM_CO2
  3       TPC_DVM_DriftVelocity
  4       TPC_DVM_FCageHV
  5       TPC_DVM_GainFar
  6       TPC_DVM_GainNear
  7       TPC_DVM_N2
  8       TPC_DVM_NumberOfSparks
  9       TPC_DVM_PeakAreaFar
  10      TPC_DVM_PeakAreaNear
  11      TPC_DVM_PeakPosFar
  12      TPC_DVM_PeakPosNear
  13      TPC_DVM_PickupHV
  14      TPC_DVM_Pressure
  15      TPC_DVM_T1_Over_P
  16      TPC_DVM_T2_Over_P
  17      TPC_DVM_T_Over_P
  18      TPC_DVM_TemperatureS1
   */
  //
  //
  //  TVectorD  vecMedian; TVectorD  vecEntries; TVectorD  vecMean; TVectorD  vecRMS;
  Double_t kEpsilon=0.0000000001;
  Double_t kBig=100000000000.;
  Int_t nsensors = goofieArray->NumSensors();
  vecEntries.ResizeTo(nsensors);
  vecMedian.ResizeTo(nsensors);
  vecMean.ResizeTo(nsensors);
  vecRMS.ResizeTo(nsensors);
  TVectorF values;
  for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
    AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
    if (gsensor &&  gsensor->GetGraph()){
      Int_t npoints = gsensor->GetGraph()->GetN();
      // filter zeroes
      values.ResizeTo(npoints);
      Int_t nused =0;
      for (Int_t ipoint=0; ipoint<npoints; ipoint++){
	if (TMath::Abs(gsensor->GetGraph()->GetY()[ipoint])>kEpsilon && 
	   TMath::Abs(gsensor->GetGraph()->GetY()[ipoint])<kBig ){
	  values[nused]=gsensor->GetGraph()->GetY()[ipoint];
	  nused++;
	}
      }
      //
      vecEntries[isensor]= nused;      
      if (nused>1){
	vecMedian[isensor] = TMath::Median(nused,values.GetMatrixArray());
	vecMean[isensor]   = TMath::Mean(nused,values.GetMatrixArray());
	vecRMS[isensor]    = TMath::RMS(nused,values.GetMatrixArray());
      }
    }
  }
}


void FilterMag(const char * runList){
  //
  //
  //
  //  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  ifstream in;
  in.open(runList);
  Int_t irun=0;
  while(in.good()) {
    in >> irun;
    if (irun==0) continue;
    AliGRPObject *grp = AliTPCcalibDB::GetGRP(irun);
    Float_t current = -1;
    Float_t bz      = -1;
    Float_t press   =  0;
    if (grp){
      current = grp->GetL3Current((AliGRPObject::Stats)0);
      bz = 5*current/30000.;
      printf("Run%d\tL3 current%f\tBz\t%f\n",irun,current,bz);
    }
    else{
      printf("Run%d\tL3 current%f\tBz\t%f\n",irun,current,bz);
    }
  }

}

/*

AliDCSSensor * sensorPressure = AliTPCcalibDB::Instance()->GetPressureSensor(62084);
entry = AliCDBManager::Instance()->Get("TPC/Calib/Temperature",run);
AliTPCSensorTempArray * tempArray = (AliTPCSensorTempArray *)entry->GetObject();
AliTPCSensorTempArray * tempArray = (AliTPCSensorTempArray *)AliTPCcalibDB::Instance()->GetTemperatureSensor(62084)
AliTPCTempMap * tempMap = new AliTPCTempMap(tempArray);
TLinearFitter * fitter = tempMap->GetLinearFitter(0,0,tempArray->GetStartTime());

AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(62084);

*/
