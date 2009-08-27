/*
.x ~/NimStyle.C
.x ~/rootlogon.C

gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");



.L $ALICE_ROOT/TPC/CalibMacros/CalibEnv.C+
Init();
CalibEnv("listAll.txt");
GetTree();

TFile f("dcsTime.root")

//
// if you want to use alien OCDB
// 
gSystem->Load("libXrdClient.so");
gSystem->Load("libNetx.so");
if (!gGrid) TGrid::Connect("alien://",0,0,"t");



*/

#include <iostream>
#include <fstream>
#include <stdio.h>
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
#include "AliTPCROC.h"
#include "AliTPCParam.h"
#include "AliTPCCalibPulser.h"
#include "AliTPCCalibPedestal.h"
#include "AliTPCCalibCE.h"
#include "AliTPCExBFirst.h"
#include "TTreeStream.h"
#include "AliTPCTempMap.h"
#include "TVectorD.h"
#include "TMatrixD.h"


TTree * dcsTree=0;

void ProcessGoofie( AliDCSSensorArray* goofieArray, TVectorD & vecEntries, TVectorD & vecMedian, TVectorD &vecMean, TVectorD &vecRMS);
void ProcessCEdata(const char* fitFormula, TVectorD &fitResultsA, TVectorD &fitResultsC);
void ProcessNoiseData(TVectorD &vNoiseMean, TVectorD &vNoiseMeanSenRegions,
                      TVectorD &vNoiseRMS, TVectorD &vNoiseRMSSenRegions,
                     Int_t &nonMaskedZero);
void ProcessPulser(Int_t &nMasked, Int_t &nonMaskedZero);
void GetProductionInfo(Int_t run, Int_t &nalien, Int_t &nRawAlien, Int_t &nlocal, Int_t &nRawLocal);
  
void Init(){
  //
  //
  //
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Parameters","local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("GRP/GRP/Data","local:///lustre/alice/alien/alice/data/2009/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Temperature","local:///lustre/alice/alien/alice/data/2009/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/HighVoltage","local:///lustre/alice/alien/alice/data/2009/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Goofie","local:///lustre/alice/alien/alice/data/2009/OCDB/");
  AliCDBManager::Instance()->SetRun(1);
}


void InitAlien(const char *path="LHC08b"){
  //
  //
  //
  TString alpath="alien://folder=/alice/data/2008/";
  alpath+=path;
  alpath+="/OCDB";
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Parameters","local://$ALICE_ROOT/OCDB");
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
    AliTPCcalibDB::Instance()->SetRun(irun);
    AliDCSSensor * sensorPressure = AliTPCcalibDB::Instance()->GetPressureSensor(irun);
    if (!sensorPressure) continue;
    AliTPCSensorTempArray * tempArray = AliTPCcalibDB::Instance()->GetTemperatureSensor(irun);
    AliTPCTempMap * tempMap = new AliTPCTempMap(tempArray);
    AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(irun);
    //
    Int_t startTime = sensorPressure->GetStartTime();
    Int_t endTime = sensorPressure->GetEndTime();
    Int_t startTimeGRP = AliTPCcalibDB::GetGRP(irun)->GetTimeStart();
    Int_t stopTimeGRP  = AliTPCcalibDB::GetGRP(irun)->GetTimeEnd();
    Int_t dtime = TMath::Max((endTime-startTime)/20,10*60);
    //CE data processing - see ProcessCEdata function for description of the results
    TVectorD fitResultsA, fitResultsC;
    ProcessCEdata("gx++gy++lx++lx**2",fitResultsA,fitResultsC);
    //noise data Processing - see ProcessNoiseData function for description of the results
    TVectorD vNoiseMean, vNoiseMeanSenRegions, vNoiseRMS, vNoiseRMSSenRegions;
    Int_t nonMaskedZero=0;
    ProcessNoiseData(vNoiseMean, vNoiseMeanSenRegions, vNoiseRMS, vNoiseRMSSenRegions, nonMaskedZero);
    //L3 data
    Float_t bz=AliTPCcalibDB::GetBz(irun);
    Char_t  l3pol=AliTPCcalibDB::GetL3Polarity(irun);
    //calibration Pulser data processing
    Int_t nMasked=0;
    Int_t nOffChannels=0;
    ProcessPulser(nMasked,nOffChannels);
    //production information
    Int_t nalien=0,nRawAlien=0,nlocal=0,nRawLocal=0;
    GetProductionInfo(irun, nalien, nRawAlien, nlocal,nRawLocal);
    //run type
    TObjString runType(AliTPCcalibDB::GetRunType(irun).Data());
    
    for (Int_t itime=startTime; itime<endTime; itime+=dtime){
      //
      TTimeStamp tstamp(itime);
      Float_t valuePressure  = calibDB->GetPressure(tstamp,irun,0);
      Float_t valuePressure2 = calibDB->GetPressure(tstamp,irun,1);
      //temperature fits
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
      //measured skirt temperatures
      TVectorD vecSkirtTempA(18);
      TVectorD vecSkirtTempC(18);
      Int_t nsenTemp=tempArray->NumSensors();
      for (Int_t isenTemp=0;isenTemp<nsenTemp;++isenTemp){
        AliTPCSensorTemp *senTemp=(AliTPCSensorTemp*)tempArray->GetSensorNum(isenTemp);
        if (senTemp->GetType()!=3) continue;
        if (TMath::Sqrt(senTemp->GetX()*senTemp->GetX()+senTemp->GetY()*senTemp->GetY())<100) continue; //only skirt, outer FC vessel
        Double_t val=senTemp->GetValue(tstamp);
        if (senTemp->GetSide()==0)
          vecSkirtTempA[senTemp->GetSector()]=val;
        else
          vecSkirtTempC[senTemp->GetSector()]=val;
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
      TVectorD voltagesIROC(36);
      TVectorD voltagesOROC(36);
      for(Int_t j=1; j<36; j++) voltagesIROC[j-1] = AliTPCcalibDB::Instance()->GetChamberHighVoltage(irun, j,tstamp);
      for(Int_t j=36; j<72; j++) voltagesOROC[j-36] = AliTPCcalibDB::Instance()->GetChamberHighVoltage(irun, j,tstamp);
      Double_t voltIROC = TMath::Median(36, voltagesIROC.GetMatrixArray());
      Double_t voltOROC = TMath::Median(36, voltagesOROC.GetMatrixArray());
      //
      Float_t  coverIA=AliTPCcalibDB::GetCoverVoltage(irun,0,itime);
      Float_t  coverIC=AliTPCcalibDB::GetCoverVoltage(irun,18,itime);
      Float_t  coverOA=AliTPCcalibDB::GetCoverVoltage(irun,36,itime);
      Float_t  coverOC=AliTPCcalibDB::GetCoverVoltage(irun,54,itime);
      Float_t  skirtA=AliTPCcalibDB::GetSkirtVoltage(irun,0,itime);
      Float_t  skirtC=AliTPCcalibDB::GetSkirtVoltage(irun,18,itime);
      Float_t  ggOffA=AliTPCcalibDB::GetGGoffsetVoltage(irun,0,itime);
      Float_t  ggOffC=AliTPCcalibDB::GetGGoffsetVoltage(irun,18,itime);
      
      
      
      //tempMap->GetLinearFitter(0,0,itime);
      TTreeStream &mistream = (*pcstream)<<"dcs"<<
        "run="<<irun<<
        "time="<<itime<<
        "startTimeGRP="<<startTimeGRP<<
        "stopTimeGRP="<<stopTimeGRP<<
        //run type
        "runType.="<<&runType<<
        // voltage setting
        "VIROC.="<<&voltagesIROC<<
        "VOROC.="<<&voltagesOROC<<
        "medianVIROC="<<voltIROC<<
        "medianVOROC="<<voltOROC<<
        "coverIA=" << coverIA <<
        "coverIC=" << coverIC <<
        "coverOA=" << coverOA <<
        "coverOC=" << coverOC <<
        "skirtA=" << skirtA <<
        "skirtC=" << skirtC <<
        "ggOffA=" << ggOffA <<
        "ggOffC=" << ggOffC <<
        //
        "ptrel0="<<ptrelative0<<  // deltaTP/TP  - A side
        "ptrel1="<<ptrelative1<<  // deltaTP/TPC - C side
        "goofie.="<<&vecGoofie<<
        "goofieE.="<<&vecEntries<<
        "goofieMean.="<<&vecMean<<
        "goofieMedian.="<<&vecMedian<<
        "goofieRMS.="<<&vecRMS<<
        //
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
        "tempSkirtA.="<<&vecSkirtTempA<<
        "tempSkirtC.="<<&vecSkirtTempC;
	
        //noise data
      mistream<<
        "meanNoise.="<<&vNoiseMean<<
        "meanNoiseSen.="<<&vNoiseMeanSenRegions<<
        "rmsNoise.="<<&vNoiseRMS<<
        "rmsNoiseSen.="<<&vNoiseRMSSenRegions<<
        "zeroNoise="<<nonMaskedZero<<
        //pulser data
        "nMasked="<<nMasked<< //should perhaps go to altro data
        "nOffPulser="<<nOffChannels<<
        //ce data
        "CEfitA.="<<&fitResultsA<<
        "CEfitC.="<<&fitResultsC<<
        // b field
        "Bz="<< bz <<
        "L3polarity="<<l3pol<<
        // production information
        "nalien="<<nalien<<
        "nRawAlien="<<nRawAlien<<
        "nlocal="<<nlocal<<
        "nRawLocal="<<nRawLocal<<
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

void ProcessCEdata(const char* fitFormula, TVectorD &fitResultsA, TVectorD &fitResultsC)
{
  //
  // Process the CE data for this run
  // the return TVectorD arrays contian the results of the fit
  //
  const Float_t irocToffset=0.2;
  const Float_t tMaxLimit=1.2;
  //retrieve CE and ALTRO data
  AliTPCCalPad *cet0=AliTPCcalibDB::Instance()->GetCETmean();
  if (!cet0){
    TString fitString(fitFormula);
    fitString.ReplaceAll("++","#");
    Int_t ndim=fitString.CountChar('#')+1;
    fitResultsA.ResizeTo(ndim);
    fitResultsC.ResizeTo(ndim);
    return;
  }
  AliTPCCalPad padT0(*cet0);
  AliTPCCalPad *padMasked=AliTPCcalibDB::Instance()->GetALTROMasked();
  //create outlier map
  AliTPCCalPad out("out","out");
  //loop over all channels
  for (UInt_t iroc=0;iroc<padT0.kNsec;++iroc){
    AliTPCCalROC *rocData=padT0.GetCalROC(iroc);
    AliTPCCalROC *rocMasked=padMasked->GetCalROC(iroc);
    AliTPCCalROC *rocOut=out.GetCalROC(iroc);
    if (!rocData) continue;
    //add time offset to IROCs
    if (iroc<AliTPCROC::Instance()->GetNInnerSector())
      rocData->Add(irocToffset);
    //select outliers
    for (UInt_t ichannel=0;ichannel<rocData->GetNchannels();++ichannel){
      if (rocMasked && rocMasked->GetValue(ichannel)) rocOut->SetValue(ichannel,1);
      Float_t valTmean=rocData->GetValue(ichannel);
      if (valTmean==0) rocOut->SetValue(ichannel,1); //exclude values that are exactly 0
      if (TMath::Abs(valTmean)>tMaxLimit) rocOut->SetValue(ichannel,1); // exclude channels with too large variations
    }
  }
  //perform fit
  TMatrixD dummy;
  Float_t chi2A,chi2C;
  padT0.GlobalSidesFit(&out,fitFormula,fitResultsA,fitResultsC,dummy,dummy,chi2A,chi2C);
}

void ProcessNoiseData(TVectorD &vNoiseMean, TVectorD &vNoiseMeanSenRegions,
                      TVectorD &vNoiseRMS, TVectorD &vNoiseRMSSenRegions,
                      Int_t &nonMaskedZero)
{
  //
  // process noise data
  // vNoiseMean/RMS contains the Mean/RMS noise of the complete TPC [0], IROCs only [1],
  //    OROCs small pads [2] and OROCs large pads [3]
  // vNoiseMean/RMSsenRegions constains the same information, but only for the sensitive regions (edge pads, corners, IROC spot)
  // 

  //set proper size and reset
  const UInt_t infoSize=4;
  vNoiseMean.ResizeTo(infoSize);
  vNoiseMeanSenRegions.ResizeTo(infoSize);
  vNoiseRMS.ResizeTo(infoSize);
  vNoiseRMSSenRegions.ResizeTo(infoSize);
  vNoiseMean.Zero();
  vNoiseMeanSenRegions.Zero();
  vNoiseRMS.Zero();
  vNoiseRMSSenRegions.Zero();
  //counters
  TVectorD c(infoSize);
  TVectorD cs(infoSize);
  //tpc parameters
  AliTPCParam par;
  par.Update();
  //retrieve noise and ALTRO data
  AliTPCCalPad *noise=AliTPCcalibDB::Instance()->GetPadNoise();
  AliTPCCalPad *padMasked=AliTPCcalibDB::Instance()->GetALTROMasked();
  //create IROC, OROC1, OROC2 and sensitive region masks
  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
    AliTPCCalROC *noiseROC=noise->GetCalROC(isec);
    UInt_t nrows=noiseROC->GetNrows();
    for (UInt_t irow=0;irow<nrows;++irow){
      UInt_t npads=noiseROC->GetNPads(irow);
      for (UInt_t ipad=0;ipad<npads;++ipad){
        Int_t masked=(Int_t)padMasked->GetCalROC(isec)->GetValue(irow,ipad);
        Float_t noiseVal=noiseROC->GetValue(irow,ipad);
        if (masked) continue; // don't use inactive pads
        //check if noise==0
        if (noiseVal==0) ++nonMaskedZero;
        Int_t cpad=(Int_t)ipad-(Int_t)npads/2;
        Int_t masksen=1; // sensitive pards are not masked (0)
        if (ipad<2||npads-ipad-1<2) masksen=0; //don't mask edge pads (sensitive)
        if (isec<AliTPCROC::Instance()->GetNInnerSector()){
          //IROCs
          if (irow>19&&irow<46){
            if (TMath::Abs(cpad)<7) masksen=0; //IROC spot
          }
          Int_t type=1;
          vNoiseMean[type]+=noiseVal;
          vNoiseRMS[type]+=noiseVal*noiseVal;
          ++c[type];
          if (!masksen){
            vNoiseMeanSenRegions[type]+=noiseVal;
            vNoiseRMSSenRegions[type]+=noiseVal*noiseVal;
            ++cs[type];
          }          
        } else {
          //OROCs
          //define sensive regions
          if ((nrows-irow-1)<3) masksen=0; //last three rows in OROCs are sensitive
          if ( irow>75 ){
            Int_t padEdge=(Int_t)TMath::Min(ipad,npads-ipad);
            if (padEdge<((((Int_t)irow-76)/4+1))*2) masksen=0; //OROC outer corners are sensitive
          }
          if ((Int_t)irow<par.GetNRowUp1()){
            //OROC1
            Int_t type=2;
            vNoiseMean[type]+=noiseVal;
            vNoiseRMS[type]+=noiseVal*noiseVal;
            ++c[type];
            if (!masksen){
              vNoiseMeanSenRegions[type]+=noiseVal;
              vNoiseRMSSenRegions[type]+=noiseVal*noiseVal;
              ++cs[type];
            }
          }else{
            //OROC2
            Int_t type=3;
            vNoiseMean[type]+=noiseVal;
            vNoiseRMS[type]+=noiseVal*noiseVal;
            ++c[type];
            if (!masksen){
              vNoiseMeanSenRegions[type]+=noiseVal;
              vNoiseRMSSenRegions[type]+=noiseVal*noiseVal;
              ++cs[type];
            }
          }
        }
        //whole tpc
        Int_t type=0;
        vNoiseMean[type]+=noiseVal;
        vNoiseRMS[type]+=noiseVal*noiseVal;
        ++c[type];
        if (!masksen){
          vNoiseMeanSenRegions[type]+=noiseVal;
          vNoiseRMSSenRegions[type]+=noiseVal*noiseVal;
          ++cs[type];
        }
      }//end loop pads
    }//end loop rows
  }//end loop sectors (rocs)
  
  //calculate mean and RMS
  const Double_t verySmall=0.0000000001;
  for (UInt_t i=0;i<infoSize;++i){
    Double_t mean=0;
    Double_t rms=0;
    Double_t meanSen=0;
    Double_t rmsSen=0;
    
    if (c[i]>verySmall){
      mean=vNoiseMean[i]/c[i];
      rms=vNoiseRMS[i];
      rms=TMath::Sqrt(TMath::Abs(rms/c[i]-mean*mean));
    }
    vNoiseMean[i]=mean;
    vNoiseRMS[i]=rms;
    
    if (cs[i]>verySmall){
      meanSen=vNoiseMeanSenRegions[i]/cs[i];
      rmsSen=vNoiseRMSSenRegions[i];
      rmsSen=TMath::Sqrt(TMath::Abs(rmsSen/cs[i]-meanSen*meanSen));
    }
    vNoiseMeanSenRegions[i]=meanSen;
    vNoiseRMSSenRegions[i]=rmsSen;
  }
}

void ProcessPulser(Int_t &nMasked, Int_t &nonMaskedZero)
{
  //
  // Process the Pulser information
  //

  //reset counters
  nonMaskedZero=0;
  nMasked=0;
  //retrieve pulser and ALTRO data
  AliTPCCalPad *pulserTmean=AliTPCcalibDB::Instance()->GetPulserTmean();
  if (!pulserTmean) return;
//   AliTPCCalPad *pulserTrms=AliTPCcalibDB::Instance()->GetPulserTrms();
//   AliTPCCalPad *pulserQmean=AliTPCcalibDB::Instance()->GetPulserQmean();
  AliTPCCalPad *padMasked=AliTPCcalibDB::Instance()->GetALTROMasked();
  //create IROC, OROC1, OROC2 and sensitive region masks
  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
    AliTPCCalROC *tmeanROC=pulserTmean->GetCalROC(isec);
    if (!tmeanROC) continue;
//     AliTPCCalROC *trmsROC=pulserTrms->GetCalROC(isec);
//     AliTPCCalROC *qmeanROC=pulserQmean->GetCalROC(isec);
    Float_t tmeanMedian=tmeanROC->GetMedian();
    UInt_t nrows=tmeanROC->GetNrows();
    for (UInt_t irow=0;irow<nrows;++irow){
      UInt_t npads=tmeanROC->GetNPads(irow);
      for (UInt_t ipad=0;ipad<npads;++ipad){
        Int_t masked=(Int_t)padMasked->GetCalROC(isec)->GetValue(irow,ipad);
        Float_t tmeanVal=tmeanROC->GetValue(irow,ipad);
//         Float_t trmsVal =trmsROC->GetValue(irow,ipad);
//         Float_t qmeanVal=qmeanROC->GetValue(irow,ipad);
        if (masked){
          ++nMasked;
          continue; // don't use inactive pads
        }
        if ( TMath::Abs(tmeanVal-tmeanMedian)>1.5 ) ++nonMaskedZero;
      }
    }
  }
}

void GetProductionInfo(Int_t run, Int_t &nalien, Int_t &nRawAlien, Int_t &nlocal, Int_t &nRawLocal){
  //
  // find number of ESDs in central and local production for run
  //

  nalien=0;
  nRawAlien=0;
  nlocal=0;
  nRawLocal=0;
  TString sNlines;
  FILE *pipe = 0x0;
  //find number of ESDs in alien
  TString command="alien_find /alice/data/2009 ";
  command += Form("%09d",run);
  command += " | grep AliESDs.root | wc -l";
  sNlines = gSystem->GetFromPipe(command.Data());
  nalien=sNlines.Atoi();
  //find number of raw files on alien
  command="alien_find /alice/data/2009 ";
  command += Form("%09d",run);
  command += " | grep raw | grep -v tag | wc -l";
  sNlines = gSystem->GetFromPipe(command.Data());
  nRawAlien=sNlines.Atoi();
  //find number of ESDs local
  command="find /lustre/alice/alien/alice/data/2009 -name AliESDs.root | grep ";
  command += Form("%09d",run);
  command += " | wc -l";
  sNlines = gSystem->GetFromPipe(command.Data());
  nlocal=sNlines.Atoi();
  //find number of local raw data files
  command="find /lustre/alice/alien/alice/data/2009 -name \"*.root\" | grep ";
  command += Form("%09d",run);
  command += " | grep raw | grep -v tag | wc -l";
  sNlines = gSystem->GetFromPipe(command.Data());
  nRawLocal=sNlines.Atoi();
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
//     Float_t press   =  0;
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


void GetTree(){
  TFile *fdcs = new TFile("dcsTime.root");
  dcsTree  = (TTree*)fdcs->Get("dcs");
  //
  // mean temp A
  
  dcsTree->Draw("temp30.fElements[0]");
  
}

void GetNominalValues(){
  //
  if (!dcsTree) return;
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

/*

void PlotPressureResol(){
  //
  // Example
  //
  dcs->Draw("100000*(press-press2-4.782)/press/sqrt(2.)>>his(100,-50,50)","run>61400","")
  his->SetXTitle("#sigma_{P/P_{0}}(10^{-5})");
  gPad->SaveAs("picDCS/deltaPoverP.eps");
  gPad->SaveAs("picDCS/deltaPoverP.gif");

}
void PlotTresol(){
  //
  // T resolution example
  // plot difference of the temperature from A and C side
  // Supposing the error is independent - (division by sqrt(2))
  dcs->Draw("100000*(temp30.fElements[0]-temp31.fElements[0]+0.00509)/(temp31.fElements[0]+273.15)/sqrt(2.)>>his(100,-5,5)","run>61400","");
  his->SetXTitle("#sigma_{T/T_{0}}(10^{-5})");
  gPad->SaveAs("picDCS/deltaToverT.eps");
  gPad->SaveAs("picDCS/deltaToverT.gif");
}
*/
