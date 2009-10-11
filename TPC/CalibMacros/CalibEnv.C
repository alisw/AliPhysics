/*
//Make a tree dump of TPC calibration:


gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
.L $ALICE_ROOT/TPC/CalibMacros/CalibEnv.C+


CalibEnv("run.list");
TFile f("dcsTime.root")
*/
 
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliLog.h>
#include <AliMagF.h>
#include "AliTPCcalibDB.h"
#include "AliTPCcalibDButil.h"
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
#include "AliTPCCalibRaw.h"
#include "AliSplineFit.h"


TTree * dcsTree=0;
TString refFile="dummy.root";
void GetProductionInfo(Int_t run, Int_t &nalien, Int_t &nRawAlien, Int_t &nlocal, Int_t &nRawLocal);
  

void CalibEnv(const char * runList, Int_t first=1, Int_t last=-1){
  //
  // runList - listOfRuns to process
  // first   - first run to process
  // last    - last  to process
  // 

  refFile=gSystem->Getenv("REF_DATA_FILE");
  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  //
  // make list of runs
  //
  ifstream in;
  in.open(runList);
  Int_t irun=0;
  TArrayI runArray(100000);
  Int_t indexes[100000];
  Int_t nruns=0;
  while(in.good()) {
    in >> irun;
    if (in.eof()) break;
    if (irun<first) continue;  // process only subset of list
    if (last>0 && irun>=last) continue;  // process only subset of list
    runArray[nruns]=irun;
    nruns++;
  }
  TMath::Sort(nruns, runArray.fArray, indexes,kFALSE);

  TTreeSRedirector *pcstream = new TTreeSRedirector("dcsTime.root");
  AliTPCcalibDButil dbutil;
  dbutil.SetRefFile(refFile.Data());
  Int_t startTime = 0;
  Int_t endTime = 0;
  Int_t startTimeGRP = 0;
  Int_t stopTimeGRP  = 0;

  AliSplineFit *fitVdrift=0x0;
  for (Int_t run=0; run<nruns; run++){
    Int_t irun=runArray[indexes[run]];
    printf("Processing run %d ...\n",irun);
    AliTPCcalibDB::Instance()->SetRun(irun);
    dbutil.UpdateFromCalibDB();
    //
    AliDCSSensorArray *arrHV=calibDB->GetVoltageSensors(irun);
    if (!arrHV) continue;
    for  (Int_t isenHV=0; isenHV<arrHV->NumSensors(); ++isenHV){
      AliDCSSensor *senHV=arrHV->GetSensorNum(isenHV);
      if (!senHV) continue;
      startTime=senHV->GetStartTime();
      endTime  =senHV->GetEndTime();
      if (startTime>0&&endTime>0) break;
    }
    if (calibDB->GetGRP(irun)){
      startTimeGRP = AliTPCcalibDB::GetGRP(irun)->GetTimeStart();
      stopTimeGRP  = AliTPCcalibDB::GetGRP(irun)->GetTimeEnd();
    }
    //    AliDCSSensor * sensorPressure = AliTPCcalibDB::Instance()->GetPressureSensor(irun);
//     if (!sensorPressure) continue;
//     Int_t startTime = sensorPressure->GetStartTime();
//     Int_t endTime = sensorPressure->GetEndTime();
//     Int_t startTimeGRP = AliTPCcalibDB::GetGRP(irun)->GetTimeStart();
//     Int_t stopTimeGRP  = AliTPCcalibDB::GetGRP(irun)->GetTimeEnd();
    AliTPCSensorTempArray * tempArray = AliTPCcalibDB::Instance()->GetTemperatureSensor(irun);
    AliTPCTempMap * tempMap = new AliTPCTempMap(tempArray);
    AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(irun);
    //
    Int_t dtime = TMath::Max((endTime-startTime)/20,10*60);
    //Goofie statistical data
    TVectorD vecEntries, vecMean, vecMedian,vecRMS;
    dbutil.ProcessGoofie(vecEntries ,vecMedian, vecMean, vecRMS);
    //CE data processing - see ProcessCEdata function for description of the results
    TVectorD fitResultsA, fitResultsC;
    Int_t nmaskedCE;
    Double_t chi2ACE=0,chi2CCE=0;
//     dbutil.ProcessCEdata("(sector<36)++gx++gy++lx++lx**2",fitResultsA,fitResultsC,nmaskedCE);
    dbutil.ProcessCEdata("(sector<36)++gy++gx++(lx-134)++(sector<36)*(lx-134)++(ly/lx)^2",fitResultsA,fitResultsC,nmaskedCE,chi2ACE,chi2CCE);
    TVectorD vecTEntries, vecTMean, vecTRMS, vecTMedian, vecQEntries, vecQMean, vecQRMS, vecQMedian;
    Float_t driftTimeA, driftTimeC;
    dbutil.ProcessCEgraphs(vecTEntries, vecTMean, vecTRMS, vecTMedian,
                           vecQEntries, vecQMean, vecQRMS, vecQMedian,
                           driftTimeA, driftTimeC );
    //drift velocity using tracks
//     fitVdrift=calibDB->GetVdriftSplineFit("ALISPLINEFIT_MEAN_VDRIFT_COSMICS_ALL",irun);
    fitVdrift=calibDB->CreateVdriftSplineFit("TGRAPHERRORS_MEAN_VDRIFT_COSMICS_ALL",irun);
    //noise data Processing - see ProcessNoiseData function for description of the results
    TVectorD vNoiseMean, vNoiseMeanSenRegions, vNoiseRMS, vNoiseRMSSenRegions;
    Int_t nonMaskedZero=0;
    dbutil.ProcessNoiseData(vNoiseMean, vNoiseMeanSenRegions, vNoiseRMS, vNoiseRMSSenRegions, nonMaskedZero);
    // comparisons
    TVectorF pedestalDeviations;
    TVectorF noiseDeviations;
    TVectorF pulserQdeviations;
    Float_t varQMean;
    Int_t npadsOutOneTB;
    Int_t npadsOffAdd;
    dbutil.ProcessPedestalVariations(pedestalDeviations);
    dbutil.ProcessNoiseVariations(noiseDeviations);
    dbutil.ProcessPulserVariations(pulserQdeviations,varQMean,npadsOutOneTB,npadsOffAdd);
    
    //L3 data 
    Float_t bz=AliTPCcalibDB::GetBz(irun);
    Char_t  l3pol=AliTPCcalibDB::GetL3Polarity(irun);
    //calibration Pulser data processing
    Int_t nOffChannels=0;
    TVectorD vTimePulser;
    nOffChannels=dbutil.GetNPulserOutliers();
    dbutil.ProcessPulser(vTimePulser);
    //ALTRO data
    Int_t nMasked=0;
    dbutil.ProcessALTROConfig(nMasked);
    //Calib RAW data
    Int_t nFailL1=-1;
    if (calibDB->GetCalibRaw())
      nFailL1=calibDB->GetCalibRaw()->GetNFailL1Phase();
    //production information
    Int_t nalien=0,nRawAlien=0,nlocal=0,nRawLocal=0;
//     GetProductionInfo(irun, nalien, nRawAlien, nlocal,nRawLocal);
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
      //goofie data
      TVectorD vecGoofie;
      if (goofieArray){
        vecGoofie.ResizeTo(goofieArray->NumSensors());
        for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
          AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
          if (gsensor){
            vecGoofie[isensor] = gsensor->GetValue(tstamp);
          }
        }
      } else {
        vecGoofie.ResizeTo(19);
      }
      Double_t ptrelative0 = AliTPCcalibDB::GetPTRelative((UInt_t)itime,irun,0);
      Double_t ptrelative1 = AliTPCcalibDB::GetPTRelative((UInt_t)itime,irun,1);
      //
      TVectorD voltagesIROC(36);
      TVectorD voltagesOROC(36);
      for(Int_t j=1; j<36; j++) voltagesIROC[j-1] = AliTPCcalibDB::Instance()->GetChamberHighVoltage(irun, j,itime);
      for(Int_t j=36; j<72; j++) voltagesOROC[j-36] = AliTPCcalibDB::Instance()->GetChamberHighVoltage(irun, j,itime);
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
      //drift velocity
      Float_t dvCorr=-5;
      if (fitVdrift) dvCorr=fitVdrift->Eval(itime);
      
      
      //tempMap->GetLinearFitter(0,0,itime);
      (*pcstream)<<"dcs"<<
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
      
      (*pcstream)<<"dcs"<<	
        //noise data
        "meanNoise.="<<&vNoiseMean<<
        "meanNoiseSen.="<<&vNoiseMeanSenRegions<<
        "rmsNoise.="<<&vNoiseRMS<<
        "rmsNoiseSen.="<<&vNoiseRMSSenRegions<<
        "zeroNoise="<<nonMaskedZero<<
        //pulser data
        "timePulser.=" << &vTimePulser <<
        "nOffPulser="<<nOffChannels<<
        //altro data
        "nMasked="<<nMasked<<
        //ce data
        "CEfitA.="<<&fitResultsA<<
        "CEfitC.="<<&fitResultsC<<
        "nmaskedCE="<<nmaskedCE<<
	"chi2ACE="<<chi2ACE<<
	"chi2CCE="<<chi2CCE<<
        //ce graph data
        "CEgrTEntries.="<<&vecTEntries<<
        "CEgrTMean.="<<&vecTMean<<
        "CEgrTRMS.="<<&vecTRMS<<
        "CEgrTMedian.="<<&vecTMedian<<
        "CEgrQEntries.="<<&vecQEntries<<
        "CEgrQMean.="<<&vecQMean<<
        "CEgrQRMS.="<<&vecQRMS<<
        "CEgrQMedian.="<<&vecQMedian<<
        "CEgrDriftA="<<driftTimeA<<
        "CEgrDriftC="<<driftTimeC<<
        //calib raw data
        "nFailL1="<<nFailL1<<
        // b field
        "Bz="<< bz <<
        "L3polarity="<<l3pol<<
        // production information
        "nalien="<<nalien<<
        "nRawAlien="<<nRawAlien<<
        "nlocal="<<nlocal<<
        "nRawLocal="<<nRawLocal<<
        //comparisons with ref data
        "pedestalDeviations.="<<&pedestalDeviations<<
        "noiseDeviations.="<<&noiseDeviations<<
        "pulserQdeviations.="<<&pulserQdeviations<<
//         "pulserVarQMean="<<varQMean<<
        "pulserNpadsOutOneTB="<<npadsOutOneTB<<
        "pulserNpadsOffAdd="<<npadsOffAdd<<
        "driftCorrCosmAll="<<dvCorr<<
      
        "\n";
    }//end run loop
//     delete fitVdrift;
//     fitVdrift=0;
  }
  delete pcstream;
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

