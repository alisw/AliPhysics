/*
// Make a summary information of calibration.
// Store results in the summary trees
// OCDB configuration

Example usage:

gSystem->Load("libANALYSIS");
gSystem->Load("libTPCcalib");

Int_t irun=119037;
gROOT->LoadMacro("$ALICE_ROOT/TPC/scripts/OCDBscan/ConfigOCDB.C");
ConfigOCDB(irun)

AliTPCcalibSummary *calibSummary = new AliTPCcalibSummary;
calibSummary->ProcessRun(irun);
delete calibSummary;

*/

#include <TROOT.h>
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
#include "TGraphErrors.h"
#include <AliCTPTimeParams.h>
#include <AliTPCcalibSummary.h>
#include <TStatToolkit.h>
#include <TCut.h> 
#include "AliTPCCalibGlobalMisalignment.h"
#include "AliTPCExBTwist.h"
#include "AliTPCComposedCorrection.h"
//
//
//
AliTPCcalibSummary::AliTPCcalibSummary():
  TNamed(),
  fCalibDB(0),
  fDButil(0),
  fPcstream(0)
{
  //
  // default constructor
  // OCDB have to be setupe before - not part of class
  // usualy ConfigOCDB.C macro used
  // 
  fPcstream = new TTreeSRedirector("dcsTime.root");
  fCalibDB = AliTPCcalibDB::Instance();
  fDButil= new AliTPCcalibDButil;
}

AliTPCcalibSummary::~AliTPCcalibSummary(){
  //
  // destructor  - close streamer
  //
  delete fPcstream;  
}

void AliTPCcalibSummary::Process(const char * runList, Int_t first, Int_t last){
  //
  // runList - listOfRuns to process
  // first   - first run to process
  // last    - last  to process
  // 
  //
  // make list of runs
  //

  ifstream inputFile;
  inputFile.open("run.list");
  Int_t irun=0;
  TArrayI runArray(100000);
  Int_t indexes[100000];
  Int_t nruns=0;
  printf("Runs to process:\n");
  if (!inputFile.is_open()) { 
    printf("Problem to open file %s\n",runList);
  }
  while( inputFile.good() ) {
    inputFile >> irun;
    printf("Run \t%d\n",irun);
    if (irun<first) continue;  // process only subset of list
    if (last>0 && irun>last) continue;  // process only subset of list
    runArray[nruns]=irun;
    nruns++;
  }


  TMath::Sort(nruns, runArray.fArray, indexes,kFALSE);
  Int_t startTime = 0;
  Int_t endTime   = 0;
  for (Int_t run=0; run<nruns; run++){
    irun=runArray[indexes[run]];
    printf("Processing run %d ...\n",irun);
    fCalibDB->SetRun(irun);
    fDButil->UpdateFromCalibDB();
    fDButil->SetReferenceRun(irun);
    fDButil->UpdateRefDataFromOCDB();
    fCalibDB->CreateGUITree("calPads.root");
    fDButil->CreateGUIRefTree("calPadsRef.root");
    //
    AliDCSSensorArray *arrHV=fCalibDB->GetVoltageSensors(irun);
    if (!arrHV) continue;
    for  (Int_t isenHV=0; isenHV<arrHV->NumSensors(); ++isenHV){
      AliDCSSensor *senHV=arrHV->GetSensorNum(isenHV);
      if (!senHV) {
	printf("Not interesting OCDB info\n");
	continue;
      }
      startTime=senHV->GetStartTime();
      endTime  =senHV->GetEndTime();
      if (startTime>0&&endTime>0) break;
    }
    AliDCSSensorArray* goofieArray = fCalibDB->GetGoofieSensors(irun);	
    if (goofieArray) fDButil->FilterGoofie(goofieArray,0.5,4.,4,10,fPcstream);
    // don't filter goofie for the moment
    ProcessRun(irun, startTime,endTime);
  }
}


void AliTPCcalibSummary::ProcessRun(Int_t irun, Int_t startTime, Int_t endTime){
  //
  // Process run irun
  // 
  fCalibDB->SetRun(irun);
  fDButil->UpdateFromCalibDB();
  fDButil->SetReferenceRun(irun);
  fDButil->UpdateRefDataFromOCDB();
  fCalibDB->CreateGUITree("calPads.root");
  fDButil->CreateGUIRefTree("calPadsRef.root");

  //
  AliSplineFit *fitVdrift=0x0;
  Int_t startTimeGRP=0, stopTimeGRP=0;
  if (fCalibDB->GetGRP(irun)){
    startTimeGRP = AliTPCcalibDB::GetGRP(irun)->GetTimeStart();
    stopTimeGRP  = AliTPCcalibDB::GetGRP(irun)->GetTimeEnd();
  }
  if (startTime==0){
    startTime=startTimeGRP;
    endTime=stopTimeGRP;
  }
  AliTPCSensorTempArray * tempArray = fCalibDB->GetTemperatureSensor(irun);
  AliTPCTempMap * tempMap = new AliTPCTempMap(tempArray);
  AliDCSSensorArray* goofieArray = fCalibDB->GetGoofieSensors(irun);
  //
  Int_t dtime = TMath::Max((endTime-startTime)/20,10);
  //
  //Goofie statistical data
  //
  TVectorD vecEntries, vecMean, vecMedian,vecRMS;
  fDButil->ProcessGoofie(vecEntries ,vecMedian, vecMean, vecRMS);
  //
  //CE data processing - see ProcessCEdata function for description of the results
  //
  TVectorD fitResultsA, fitResultsC;
  Int_t nmaskedCE;
  Double_t chi2ACE=0,chi2CCE=0;
  //     fDButil->ProcessCEdata("(sector<36)++gx++gy++lx++lx**2",fitResultsA,fitResultsC,nmaskedCE);
  fDButil->ProcessCEdata("(sector<36)++gx++gy++(lx-134)++(sector<36)*(lx-134)++(ly/lx)^2",fitResultsA,fitResultsC,nmaskedCE,chi2ACE,chi2CCE);

  TVectorD fitCEResultsA(7), fitCEResultsC(7);
  Int_t    noutCE;
  Double_t chi2CEA=0,chi2CEC=0;
  AliTPCCalPad *time0 = fDButil->CreatePadTime0CE(fitCEResultsA, fitCEResultsC, noutCE, chi2CEA, chi2CEC);
  delete time0;
  //  
  //
  TVectorD vecTEntries, vecTMean, vecTRMS, vecTMedian, vecQEntries, vecQMean, vecQRMS, vecQMedian;
  Float_t driftTimeA, driftTimeC;
  fDButil->ProcessCEgraphs(vecTEntries, vecTMean, vecTRMS, vecTMedian,
			  vecQEntries, vecQMean, vecQRMS, vecQMedian,
			  driftTimeA, driftTimeC );
  //
  //
  //
  //drift velocity using tracks
  //
  //     fitVdrift=fCalibDB->GetVdriftSplineFit("ALISPLINEFIT_MEAN_VDRIFT_COSMICS_ALL",irun);
  fitVdrift=fCalibDB->CreateVdriftSplineFit("TGRAPHERRORS_MEAN_VDRIFT_COSMICS_ALL",irun);
  //noise data Processing - see ProcessNoiseData function for description of the results
  TVectorD vNoiseMean, vNoiseMeanSenRegions, vNoiseRMS, vNoiseRMSSenRegions;
  Int_t nonMaskedZero=0, nNaN=0;
  fDButil->ProcessNoiseData(vNoiseMean, vNoiseMeanSenRegions, vNoiseRMS, vNoiseRMSSenRegions, nonMaskedZero, nNaN);
  //
  // comparisons
  //
  TVectorF pedestalDeviations;
  TVectorF noiseDeviations;
  TVectorF pulserQdeviations;
  Float_t varQMean;
  Int_t npadsOutOneTB;
  Int_t npadsOffAdd;
  fDButil->ProcessPedestalVariations(pedestalDeviations);
  fDButil->ProcessNoiseVariations(noiseDeviations);
  fDButil->ProcessPulserVariations(pulserQdeviations,varQMean,npadsOutOneTB,npadsOffAdd);
  //
  //L3 data 
  //
  Float_t bz=AliTPCcalibDB::GetBz(irun);               
  Char_t  l3pol=AliTPCcalibDB::GetL3Polarity(irun);    
  //
  //QA data processing
  //
  TVectorD vQaOcc;
  TVectorD vQaQtot;
  TVectorD vQaQmax;
  fDButil->ProcessQAData(vQaOcc, vQaQtot, vQaQmax);
  //
  //calibration Pulser data processing
  //
  Int_t nOffChannels=0;
  TVectorD vTimePulser;
  nOffChannels=fDButil->GetNPulserOutliers();
  fDButil->ProcessPulser(vTimePulser);
  //
  //ALTRO data
  //
  Int_t nMasked=0;
  fDButil->ProcessALTROConfig(nMasked);
  //
  //Calib RAW data
  //
  Int_t nFailL1=-1;
  if (fCalibDB->GetCalibRaw()) nFailL1=fCalibDB->GetCalibRaw()->GetNFailL1Phase();
  //
  //production information
  //
  Int_t nalien=0,nRawAlien=0,nlocal=0,nRawLocal=0;
  //run type
  TObjString runType(AliTPCcalibDB::GetRunType(irun).Data());
  //
  //
  //
  
  for (Int_t itime=startTime; itime<endTime; itime+=dtime){
    //
    TTimeStamp tstamp(itime);
    Float_t valuePressure  = fCalibDB->GetPressure(tstamp,irun,0);
    Float_t valuePressure2 = fCalibDB->GetPressure(tstamp,irun,1);
    Double_t ptrelative0 = AliTPCcalibDB::GetPTRelative((UInt_t)itime,irun,0);
    Double_t ptrelative1 = AliTPCcalibDB::GetPTRelative((UInt_t)itime,irun,1);
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
    //
    //measured skirt temperatures
    //
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
    //
    //goofie data
    //
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
    //
    static TVectorF voltagesIROC(36);
    static TVectorF voltagesIROCMedian(36);
    static TVectorF voltagesIROCNominal(36);
    static TVectorF voltagesIROCCurrentNominal(36);
    static TVectorF voltagesIROCStatus(36);
    static TVectorF voltagesIROCGoodFraction(36);
    //
    static TVectorF voltagesOROC(36);
    static TVectorF voltagesOROCMedian(36);
    static TVectorF voltagesOROCNominal(36);
    static TVectorF voltagesOROCCurrentNominal(36);
    static TVectorF voltagesOROCStatus(36);
    static TVectorF voltagesOROCGoodFraction(36);
    
    for(Int_t j=0; j<36; j++){
      voltagesIROC[j]               = fCalibDB->GetChamberHighVoltage(irun, j,itime);
      voltagesIROCMedian[j]         = fCalibDB->GetChamberHighVoltageMedian(j);
      voltagesIROCNominal[j]        = fCalibDB->GetParameters()->GetNominalVoltage(j);
      voltagesIROCCurrentNominal[j] = fCalibDB->GetChamberCurrentNominalHighVoltage(j);
      voltagesIROCStatus[j]         = fCalibDB->GetChamberHVStatus(j);
      voltagesIROCGoodFraction[j]   = fCalibDB->GetChamberGoodHighVoltageFraction(j);
    }
    
    for(Int_t j=36; j<72; j++) {
      voltagesOROC[j-36]               = fCalibDB->GetChamberHighVoltage(irun, j,itime);
      voltagesOROCMedian[j-36]         = fCalibDB->GetChamberHighVoltageMedian(j);
      voltagesOROCNominal[j-36]        = fCalibDB->GetParameters()->GetNominalVoltage(j);
      voltagesOROCCurrentNominal[j-36] = fCalibDB->GetChamberCurrentNominalHighVoltage(j);
      voltagesOROCStatus[j-36]         = fCalibDB->GetChamberHVStatus(j);
      voltagesOROCGoodFraction[j-36]   = fCalibDB->GetChamberGoodHighVoltageFraction(j);
    }
    
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
    //data taking active
    Bool_t dataTakingActive=fCalibDB->IsDataTakingActive((time_t)itime);
    
    //tempMap->GetLinearFitter(0,0,itime);
    (*fPcstream)<<"dcs"<<
      "run="<<irun<<
      "time="<<itime<<
      "startTimeGRP="<<startTimeGRP<<
      "stopTimeGRP="<<stopTimeGRP<<
      "dataTakingActive="<<dataTakingActive<<
      //run type
      "runType.="<<&runType<<
      // voltage setting
      "VIROC.="               << &voltagesIROC<<
      "VIROCMedian.="         << &voltagesIROCMedian<<
      "VIROCNominal.="        << &voltagesIROCNominal <<
      "VIROCCurrentNominal.=" << &voltagesIROCCurrentNominal <<
      "VIROCGoodHVFraction.=" << &voltagesIROCGoodFraction <<
      "VIROCStatus.="         << &voltagesIROCStatus <<
      //
      "VOROC.="               << &voltagesOROC<<
      "VOROCMedian.="         << &voltagesOROCMedian<<
      "VOROCNominal.="        << &voltagesOROCNominal <<
      "VOROCCurrentNominal.=" << &voltagesOROCCurrentNominal <<
      "VOROCGoodHVFraction.=" << &voltagesOROCGoodFraction <<
      "VOROCStatus.="         << &voltagesOROCStatus <<
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

    ProcessDrift(irun, itime);
    ProcessDriftCE(irun,itime);
    ProcessDriftAll(irun,itime);
    //    ProcessKryptonTime(irun,itime);
    ProcessCTP(irun,itime);
    ProcessAlign(irun,itime);
    ProcessGain(irun,itime);
    //ProcessDriftCERef();
    //ProcessPulserRef();
    //ProcessCurrent(irun,itime);


    (*fPcstream)<<"dcs"<<	
      //noise data
      "meanNoise.="<<&vNoiseMean<<
      "meanNoiseSen.="<<&vNoiseMeanSenRegions<<
      "rmsNoise.="<<&vNoiseRMS<<
      "rmsNoiseSen.="<<&vNoiseRMSSenRegions<<
      "zeroNoise="<<nonMaskedZero<<
      "nNaN="<<nNaN<<  
      //QA data
      "occQA.="  << &vQaOcc  <<
      "qQA.="    << &vQaQtot <<
      "qmaxQA.=" << &vQaQmax <<
      //pulser data
      "timePulser.=" << &vTimePulser <<
      "nOffPulser="<<nOffChannels<<
      //altro data
      "nMasked="<<nMasked<<
      //ce data -Jens version
      "CEfitA.="<<&fitResultsA<<
      "CEfitC.="<<&fitResultsC<<
      "nmaskedCE="<<nmaskedCE<<
      "chi2ACE="<<chi2ACE<<
      "chi2CCE="<<chi2CCE<<
      //ce data new - MI version
      "CEfitAMI.="<<&fitCEResultsA<<
      "CEfitCMI.="<<&fitCEResultsC<<
      "chi2CEA="<<chi2CEA<<
      "chi2CEC="<<chi2CEC<<
      //
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
}






void AliTPCcalibSummary::ProcessDrift(Int_t run, Int_t timeStamp){
  //
  // dump drift calibration data to the tree
  //
  TObjArray *array =fCalibDB->GetTimeVdriftSplineRun(run);
  TGraphErrors *laserA[3]={0,0,0};
  TGraphErrors *laserC[3]={0,0,0};
  TGraphErrors *cosmicAll=0;
  static Double_t     vlaserA[3]={0,0,0};
  static Double_t     vlaserC[3]={0,0,0};
  static Double_t     vcosmicAll=0;
  static Double_t     vdrift1=0;
  vdrift1=fCalibDB->GetVDriftCorrectionTime(timeStamp,run,0,1);

  if (array){
    laserA[0]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DELAY_LASER_ALL_A");
    laserA[1]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_A");
    laserA[2]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_A");
    laserC[0]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DELAY_LASER_ALL_C");
    laserC[1]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_C");
    laserC[2]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_C");
    cosmicAll =(TGraphErrors*)array->FindObject("TGRAPHERRORS_MEAN_VDRIFT_COSMICS_ALL");
  }
  if (laserA[0]) vlaserA[0]= AliTPCcalibDButil::EvalGraphConst(laserA[0],timeStamp);
  if (laserA[1]) vlaserA[1]= AliTPCcalibDButil::EvalGraphConst(laserA[1],timeStamp);
  if (laserA[2]) vlaserA[2]= AliTPCcalibDButil::EvalGraphConst(laserA[2],timeStamp);
  if (laserC[0]) vlaserC[0]= AliTPCcalibDButil::EvalGraphConst(laserC[0],timeStamp);
  if (laserC[1]) vlaserC[1]= AliTPCcalibDButil::EvalGraphConst(laserC[1],timeStamp);
  if (laserC[2]) vlaserC[2]= AliTPCcalibDButil::EvalGraphConst(laserC[2],timeStamp);
  if (cosmicAll) vcosmicAll= AliTPCcalibDButil::EvalGraphConst(cosmicAll,timeStamp); 
  (*fPcstream)<<"dcs"<<
    "vlaserA0="<<vlaserA[0]<<
    "vlaserA1="<<vlaserA[1]<<
    "vlaserA2="<<vlaserA[2]<<
    "vlaserC0="<<vlaserC[0]<<
    "vlaserC1="<<vlaserC[1]<<
    "vlaserC2="<<vlaserC[2]<<
    "vcosmicAll="<<vcosmicAll<<
    "vdrift1="<<vdrift1;

  //
  // define distance to measurement
  //
  static Double_t dlaserA=0; 
  static Double_t dlaserC=0; 
  static Double_t dcosmic=0; 
  static Double_t slaserA=0; 
  static Double_t slaserC=0; 
  static Double_t scosmic=0; 
  static Double_t  vclaserA[3]={0,0,0};
  static Double_t  vclaserC[3]={0,0,0};
  static Double_t  vccosmicAll=0;
  for (Int_t i=0;i<3;i++){
    if (laserA[i]) AliTPCcalibDButil::GetNearest(laserA[i],timeStamp,dlaserA,vclaserA[i]);
    if (laserC[i]) AliTPCcalibDButil::GetNearest(laserC[i],timeStamp,dlaserC,vclaserC[i]);
  }  
  if (cosmicAll) AliTPCcalibDButil::GetNearest(cosmicAll,timeStamp,dcosmic,vccosmicAll);
  (*fPcstream)<<"dcs"<<
    "vclaserA0="<<vclaserA[0]<<
    "vclaserA1="<<vclaserA[1]<<
    "vclaserA2="<<vclaserA[2]<<
    "vclaserC0="<<vclaserC[0]<<
    "vclaserC1="<<vclaserC[1]<<
    "vclaserC2="<<vclaserC[2]<<
    "vccosmicAll="<<vccosmicAll<<
    "dlaserA="<<dlaserA<<
    "dlaserC="<<dlaserC<<
    "dcosmic="<<dcosmic<<
    "slaserA="<<slaserA<<
    "slaserC="<<slaserC<<
    "scosmic="<<scosmic;

  static TGeoMatrix * matrixAlign=0;
  static Double_t twistX=0;
  static Double_t twistY=0;
  if (matrixAlign==0){
    AliTPCComposedCorrection * corr =  (AliTPCComposedCorrection *)array->FindObject("FitCorrectionTime");
    if (!corr) {
      matrixAlign=new TGeoHMatrix;
      
    }
    if (corr){
       AliTPCCalibGlobalMisalignment *align = (AliTPCCalibGlobalMisalignment*)corr->GetCorrections()->FindObject("FitAlignTime");
       AliTPCExBTwist *twist  = (AliTPCExBTwist*)corr->GetCorrections()->FindObject("FitExBTwistTime");
       if (twist){
	 twistX=twist->GetXTwist();
	 twistY=twist->GetYTwist();
	 //delete twist;
       }
       if (align && align->GetAlignGlobal()){
	 matrixAlign =  (TGeoMatrix*) (align->GetAlignGlobal()->Clone());	 
	 //delete align;
       }
    }    
  }
  (*fPcstream)<<"dcs"<<
    "alignTime.="<<matrixAlign<<
    "twistX="<<twistX<<
    "twistY="<<twistY;
}

void AliTPCcalibSummary::ProcessDriftCE(Int_t run,Int_t timeStamp){
  //
  // dump drift calibration data CE
  //
  TObjArray *arrT=fCalibDB->GetCErocTtime();
  AliTPCParam *param=fCalibDB->GetParameters();
  static TVectorD tdriftCE(74);
  static TVectorD tndriftCE(74);
  static TVectorD vdriftCE(74);
  static TVectorD tcdriftCE(74);
  static TVectorD tddriftCE(74);
  static Double_t ltime0A=0.;
  static Double_t ltime0C=0.;
  //
  //
  //
  ltime0A  = fDButil->GetLaserTime0(run,timeStamp,36000,0);
  ltime0C  = fDButil->GetLaserTime0(run,timeStamp,36000,1);
  //
  for (Int_t i=0; i<arrT->GetEntries();i++){
    tdriftCE[i]=0;
    vdriftCE[i]=0;
    TGraph *graph = (TGraph*)arrT->At(i);
    if (!graph) continue;
    tdriftCE[i]=AliTPCcalibDButil::EvalGraphConst(graph,timeStamp);
    Double_t deltaT,gry;
    AliTPCcalibDButil::GetNearest(graph,timeStamp,deltaT,gry);
    tndriftCE[i]=graph->GetN();
    tcdriftCE[i]=gry;	       
    tddriftCE[i]=deltaT;	       
    if (i%36<18){
      vdriftCE[i] =(param->GetZLength(i)/(tdriftCE[i]*param->GetTSample()*(1.-ltime0A)-param->GetL1Delay()))/param->GetDriftV();
    }else{
      vdriftCE[i] =(param->GetZLength(i)/(tdriftCE[i]*param->GetTSample()*(1.-ltime0A)-param->GetL1Delay()))/param->GetDriftV();
    }
  }
  // export values
  (*fPcstream)<<"dcs"<<  
    "tdriftCE.="<<&tdriftCE<<    // arrival time 
    "vdriftCE.="<<&vdriftCE<<    // derived drift velocity per chamber 
    "tndriftCE.="<<&tndriftCE<<  // number of points (chambers)
    "tcdriftCE.="<<&tcdriftCE<<  // constant evaluation - nearest point used
    "tddriftCE.="<<&tddriftCE<<  // distance to closest measuement
    "ltime0A="<<ltime0A<<        // laser offset expected in reconstruction
    "ltime0C="<<ltime0C;         // laser offset expected in reconstruction 
   }


void AliTPCcalibSummary::ProcessDriftAll(Int_t run,Int_t timeStamp){
  //
  // dump drift calibration data  all calibrations form DB util
  // test of utils
  static Double_t vdriftCEA=0, vdriftCEC=0, vdriftCEM=0;
  static Double_t vdriftLTA=0, vdriftLTC=0, vdriftLTM=0;
  static Double_t vdriftLTAon=0, vdriftLTCon=0, vdriftLTMon=0;
  static Double_t vdriftITS=0;
  static Double_t vdriftP=0;
  static Double_t dcea=0, dcec=0, dcem=0,  dla=0,dlc=0,dlm=0, dlaon=0,dlcon=0,dlmon=0, dp=0;
  static Double_t dits=0;
  static Double_t ltime0A=0.;
  static Double_t ltime0C=0.;
  static Double_t ctime0=0.;
  static Double_t vdrift1=0;
  vdrift1= fCalibDB->GetVDriftCorrectionTime(timeStamp,run,0,1);
  vdriftP = fDButil->GetVDriftTPC(dp, run, timeStamp, 86400, 3600,0);
  ctime0 = AliTPCcalibDButil::GetTriggerOffsetTPC(run,timeStamp, 36000, 3600,0);
  //
  vdriftCEA= fDButil->GetVDriftTPCCE(dcea,run,timeStamp,36000,0);
  vdriftCEC= fDButil->GetVDriftTPCCE(dcec,run,timeStamp,36000,1);
  vdriftCEM= fDButil->GetVDriftTPCCE(dcem,run,timeStamp,36000,2);
  //
  vdriftLTA= fDButil->GetVDriftTPCLaserTracks(dla,run,timeStamp,36000,0);
  vdriftLTC= fDButil->GetVDriftTPCLaserTracks(dlc,run,timeStamp,36000,1);
  vdriftLTM= fDButil->GetVDriftTPCLaserTracks(dlm,run,timeStamp,36000,2);
  //
  vdriftLTAon= fDButil->GetVDriftTPCLaserTracksOnline(dlaon,run,timeStamp,36000,0);
  vdriftLTCon= fDButil->GetVDriftTPCLaserTracksOnline(dlcon,run,timeStamp,36000,1);
  vdriftLTMon= fDButil->GetVDriftTPCLaserTracksOnline(dlmon,run,timeStamp,36000,2);
  //
  vdriftITS= fDButil->GetVDriftTPCITS(dits, run,timeStamp);
  //
  ltime0A  = fDButil->GetLaserTime0(run,timeStamp,36000,0);
  ltime0C  = fDButil->GetLaserTime0(run,timeStamp,36000,1);

  (*fPcstream)<<"dcs"<<  
    //
    "vdriftCEA="<<vdriftCEA<<   // drift velocity CE
    "vdriftCEC="<<vdriftCEC<<
    "vdriftCEM="<<vdriftCEM<<
    "dcea="<<dcea<<
    "dcec="<<dcec<<
    "dcem="<<dcem<<
    "vdriftLTA="<<vdriftLTA<<   // drift velocity laser tracks
    "vdriftLTC="<<vdriftLTC<<
    "vdriftLTM="<<vdriftLTM<<
    "dla="<<dla<<
    "dlc="<<dlc<<
    "dlm="<<dlm<<
    "vdriftLTAon="<<vdriftLTAon<<   // drift velocity laser tracks and CE from online algorithm
    "vdriftLTCon="<<vdriftLTCon<<
    "vdriftLTMon="<<vdriftLTMon<<
    "dlaOn="<<dlaon<<
    "dlcOn="<<dlcon<<
    "dlmOn="<<dlmon<<
    //
    //
    "vdriftITS="<<vdriftITS<<
    "dits="<<dits<<
    "ctime0="<<ctime0<<
    "vdriftP="<<vdriftP<<       // drift velocity comsic 
    "dp="<<dp<<
    "vdrift1="<<vdrift1;        // combined drift velocity

}



void AliTPCcalibSummary::ProcessKryptonTime(Int_t run, Int_t timeStamp){
  //
  // Dumping  krypton calibration results
  //
  static TObjArray * krArray=0;
  if (!krArray) {
    AliCDBEntry* entry = AliCDBManager::Instance()->Get("TPC/Calib/TimeGainKrypton", run);
    if (entry){
      krArray = (TObjArray*)entry->GetObject();
    }
  }
  static TVectorD krMean(74);
  static TVectorD krErr(74);
  static TVectorD krDist(74);
  TGraphErrors *gr=0;
  Double_t deltaT=0,gry=0;
  if (krArray){
    for (Int_t isec=0; isec<72; isec++){
      krMean[isec]=0;
      krDist[isec]=0;
      krErr[isec]=0;
      gr=(TGraphErrors*)krArray->At(isec);
      if (gr) {
	krMean[isec]=AliTPCcalibDButil::EvalGraphConst(gr,timeStamp);
	AliTPCcalibDButil::GetNearest(gr, timeStamp,deltaT,gry);
	krDist[isec]=deltaT;
      }     
      if (72+isec<krArray->GetEntries()) {
	gr=(TGraphErrors*)krArray->At(72+isec);
	if (gr) krErr[isec]=AliTPCcalibDButil::EvalGraphConst(gr,timeStamp);
      }
    }
    krMean[72]= TMath::Median(36,krMean.GetMatrixArray());
    krMean[73]= TMath::Median(36,&(krMean.GetMatrixArray()[36]));
    krErr[72]= TMath::Median(36,krErr.GetMatrixArray());
    krErr[73]= TMath::Median(36,&(krErr.GetMatrixArray()[36]));
  }
  (*fPcstream)<<"dcs"<<
    "krMean.="<<&krMean<<
    "krErr.="<<&krErr<<
    "krDist.="<<&krDist;
}

void AliTPCcalibSummary::ProcessCTP(Int_t irun, Int_t timeStamp){
  //
  // 
  //
  static TClonesArray *pcarray = new TClonesArray("AliCTPInputTimeParams",1);
  static AliCTPTimeParams *ctpParams =0;
  ctpParams = fCalibDB->GetCTPTimeParams(); //
  const TObjArray        *arr = ctpParams->GetInputTimeParams();
  pcarray->ExpandCreateFast(TMath::Max(arr->GetEntries(),1));
  for (Int_t i=0; i<arr->GetEntries(); i++){
    new ((*pcarray)[i]) AliCTPInputTimeParams(*((AliCTPInputTimeParams*)arr->At(i)));
  }
  (*fPcstream)<<"ctp"<<
    "run="<<irun<<
    "time="<<timeStamp<<
    "ctpP.="<<ctpParams<<
    "ctpArr="<<pcarray<<
    "\n";
  (*fPcstream)<<"dcs"<<
    "ctpP.="<<ctpParams<<
    "ctpArr="<<pcarray;
}

void AliTPCcalibSummary::ProcessGain(Int_t irun, Int_t timeStamp){
  //
  // Dump gain related information to the tree
  //
  static Float_t  gainCosmic = 0;
  static Float_t  gainMIP = 0;
  static Float_t  attachMIP = 0;
  static Double_t dMIP=0; 
  Double_t dummy=0;
  static TVectorD vGainGraphIROC(36);
  static TVectorD vGainGraphOROCmed(36);
  static TVectorD vGainGraphOROClong(36);
  static TVectorD vGainGraphIROCErr(36);
  static TVectorD vGainGraphOROCmedErr(36);
  static TVectorD vGainGraphOROClongErr(36);
  
  vGainGraphIROC.Zero();
  vGainGraphOROCmed.Zero();
  vGainGraphOROClong.Zero();
  vGainGraphIROCErr.Zero();
  vGainGraphOROCmedErr.Zero();
  vGainGraphOROClongErr.Zero();
  
  TGraphErrors grDummy;
  TObjArray * gainSplines = fCalibDB->GetTimeGainSplinesRun(irun);
  if (gainSplines) {
    TGraphErrors * graphMIP = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL");
    TGraphErrors * graphCosmic = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL");
    TGraphErrors * graphAttach = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_ATTACHMENT_BEAM_ALL");
    //
    TGraphErrors * graphGainIROC       = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_CHAMBERGAIN_SHORT_BEAM_ALL");
    TGraphErrors * graphGainOROCMedium = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_CHAMBERGAIN_MEDIUM_BEAM_ALL");
    TGraphErrors * graphGainOROCLong   = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_CHAMBERGAIN_LONG_BEAM_ALL");

    if (graphGainIROC && graphGainOROCMedium && graphGainOROCLong) {
      Double_t x=0,y=0;
      for (Int_t i=0; i<36; ++i){
        graphGainIROC->GetPoint(i,x,y);
        vGainGraphIROC(i)=y;
        graphGainOROCMedium->GetPoint(i,x,y);
        vGainGraphOROCmed(i)=y;
        graphGainOROCLong->GetPoint(i,x,y);
        vGainGraphOROClong(i)=y;
        //errors
        vGainGraphIROCErr(i)     = graphGainIROC->GetEY()[i];
        vGainGraphOROCmedErr(i)  = graphGainOROCMedium->GetEY()[i];
        vGainGraphOROClongErr(i) = graphGainOROCLong->GetEY()[i];
      }
    }
    
    if (graphMIP) gainMIP = AliTPCcalibDButil::EvalGraphConst(graphMIP,timeStamp);
    if (graphCosmic) gainCosmic = AliTPCcalibDButil::EvalGraphConst(graphCosmic,timeStamp);
    if (graphAttach) attachMIP = AliTPCcalibDButil::EvalGraphConst(graphAttach,timeStamp);
    if (graphMIP)  AliTPCcalibDButil::GetNearest(graphMIP, timeStamp, dMIP,dummy);    
  }
    
  // time dependence of gain
  (*fPcstream)<<"dcs"<<
    "rocGainIROC.="            << &vGainGraphIROC        <<
    "rocGainOROCMedium.="      << &vGainGraphOROCmed     <<
    "rocGainOROCLong.="        << &vGainGraphOROClong    <<
    "rocGainErrIROC.="         << &vGainGraphIROCErr     <<
    "rocGainErrOROCMedium.="   << &vGainGraphOROCmedErr  <<
    "rocGainErrOROCLong.="     << &vGainGraphOROClongErr <<
    "gainMIP="                 << gainMIP                <<
    "attachMIP="               << attachMIP              <<
    "dMIP="                    << dMIP                   <<
    "gainCosmic="              << gainCosmic;
}


void AliTPCcalibSummary::ProcessAlign(Int_t run, Int_t timeStamp){
  //
  // Proccess alignment 
  //
  TString   grnames[12]={"ALIGN_ITS", "ALIGN_ITSP", "ALIGN_ITSM", "ALIGN_ITSB",
			 "ALIGN_TRD", "ALIGN_TRDP", "ALIGN_TRDM","ALIGN_TRDB",
			 "ALIGN_TOF", "ALIGN_TOFP", "ALIGN_TOFM","ALIGN_TOFB"};
  TString   grpar[9]={"DELTAPSI", "DELTATHETA", "DELTAPHI",
		      "DELTAX", "DELTAY", "DELTAZ",
		      "DRIFTVD", "T0", "VDGY"};
  static Double_t values[12*9];
  static Double_t errs[12*9];

  TObjArray *array =fCalibDB->GetTimeVdriftSplineRun(run);
  TGraphErrors *gr=0;
  for (Int_t idet=0; idet<12; idet++){
    for (Int_t ipar=0; ipar<9; ipar++){
      TString grName=grnames[idet];
      grName+="_TPC_";
      grName+=grpar[ipar];
      if (array){
	gr = (TGraphErrors*)array->FindObject(grName.Data());
      }
      values[9*idet+ipar]=0;
      errs[9*idet+ipar]=0;
      if (gr) values[9*idet+ipar]=AliTPCcalibDButil::EvalGraphConst(gr,timeStamp);
      (*fPcstream)<<"dcs"<<
	Form("%s=",grName.Data())<<values[9*idet+ipar];
      (*fPcstream)<<"align"<<
	Form("%s=",grName.Data())<<values[9*idet+ipar];
    }
  }
  (*fPcstream)<<"align"<<
    "time="<<timeStamp<<
    "run="<<run<<
    "\n";
}


void AliTPCcalibSummary::ProcessDriftCERef(){
  //
  // Get fit of residuals if CE in respect with reference
  // data
  //
  static TVectorD  sec(72);
  static TVectorD vec0(72);
  static TVectorD vecLy(72);
  static TVectorD vecLx(72);
  static TVectorD vecChi2(72);
  static TVectorD vecN(72);
  //
  static TVectorD vecA0(72);
  static TVectorD vecALy(72);
  static TVectorD vecALx(72);
  static TVectorD vecAChi2(72);
  //
  static TVectorD vecASide(4);
  static TVectorD vecCSide(4);
  static Bool_t isCalc=kFALSE;
  
  TFile f("calPads.root");
  TFile fref("calPadsRef.root");
  TTree * tree = (TTree*)f.Get("calPads");
  TTree * treeRef = (TTree*)fref.Get("calPads");
  tree->AddFriend(treeRef,"R");
  tree->SetAlias("inCE","((CEQmean.fElements>35)&&abs(CETmean.fElements)<1.5&&abs(CETrms.fElements/1.2-1)<0.2)");  // outlyerTrms
  tree->SetAlias("inCER","((R.CEQmean.fElements>35)&&abs(R.CETmean.fElements)<1.5&&abs(R.CETrms.fElements/1.2-1)<0.2)");  // outlyerTrms
  //
  if (!isCalc){
    // make fits only once
    TStatToolkit toolkit;
    Double_t chi2=0;
    Int_t    npoints=0;
    TVectorD param;
    TMatrixD covar;
    tree->SetAlias("dt","CETmean.fElements-R.CETmean.fElements");
    TCut cutAll ="inCE&&inCER&&abs(CETmean.fElements-R.CETmean.fElements)<0.5"; 
    TString  fstringG="";              // global part
    fstringG+="ly.fElements++";
    fstringG+="(lx.fElements-134.)++";  
    for (Int_t isec=0; isec<72; isec++){
      TStatToolkit::FitPlane(tree,"0.264*dt", fstringG.Data(),Form("sector==%d",isec)+cutAll, chi2,npoints,param,covar,-1,0, 10000000, kFALSE);
      if (npoints<3) continue;
      printf("Sector=%d\n",isec);
      vec0[isec]=param[0];
      vecLy[isec]=param[1];
      vecLx[isec]=param[2];
      sec[isec]=isec;
      vecN[isec]=npoints;
      vecChi2[isec]=TMath::Sqrt(chi2/npoints);

      TStatToolkit::FitPlane(tree,"0.264*CETmean.fElements", fstringG.Data(),Form("sector==%d",isec)+cutAll, chi2,npoints,param,covar,-1,0, 10000000, kFALSE);
      if (npoints<3) continue;
      printf("Sector=%d\n",isec);
      vecA0[isec]=param[0];
      vecALy[isec]=param[1];
      vecALx[isec]=param[2];
      vecAChi2[isec]=TMath::Sqrt(chi2/npoints);
    }
    isCalc=kTRUE;
    //
    TString  fstringRef="";              // global fit
    fstringRef+="gx.fElements++";
    fstringRef+="gy.fElements++";  
    fstringRef+="lx.fElements++";
    TStatToolkit::FitPlane(tree,"0.264*dt", fstringG.Data(),"(sector%36)<18"+cutAll, chi2,npoints,vecASide,covar,-1,0, 10000000, kFALSE);
    TStatToolkit::FitPlane(tree,"0.264*dt", fstringG.Data(),"(sector%36)>=18"+cutAll, chi2,npoints,vecCSide,covar,-1,0, 10000000, kFALSE);

  }
  (*fPcstream)<<"dcs"<<     // CE information
    "CETSector.="<<&sec<<    // sector numbers
    "CETRefA.="<<&vecASide<<   // diff to reference A side
    "CETRefC.="<<&vecCSide<<   // diff to reference C side    
    //                      // fit in respect to reference data
    "CETRef0.="<<&vec0<<    // offset change
    "CETRefY.="<<&vecLy<<   // slope y change - rad
    "CETRefX.="<<&vecLx<<   // slope x change - rad
    "CETRefChi2.="<<&vecChi2<<    // chi2 (rms in cm)
    "CETRefN.="<<&vecN<<     //number of accepted points
    //                       // fit in respect per mean per side
    "CET0.="<<&vecA0<<       // offset change
    "CETY.="<<&vecALy<<      // slope y change - rad
    "CETX.="<<&vecALx<<      // slope x change - rad
    "CETChi2.="<<&vecAChi2;  // chi2 (rms in cm)
}

void AliTPCcalibSummary::ProcessPulserRef(){
  //
  // Get fit of residuals if Pulser in respect with reference
  // data
  //
  static TVectorD  sec(72);
  static TVectorD vec0(72);
  static TVectorD vecLy(72);
  static TVectorD vecLx(72);
  static TVectorD vecChi2(72);
  static TVectorD vecN(72);
  //
  static TVectorD vecA0(72);
  static TVectorD vecALy(72);
  static TVectorD vecALx(72);
  static TVectorD vecAChi2(72);
  static Bool_t isCalc=kFALSE;
  
  TFile f("calPads.root");
  TFile fref("calPadsRef.root");
  TTree * tree = (TTree*)f.Get("calPads");
  TTree * treeRef = (TTree*)fref.Get("calPads");
  tree->AddFriend(treeRef,"R");
  
  tree->SetAlias("inPulser","(abs(PulserTmean.fElements-PulserTmean_Median)<1.5&&abs(PulserTrms.fElements-PulserTrms_Median)<0.2)");  // outlyerTrms
  tree->SetAlias("inPulserR","(abs(R.PulserTmean.fElements-R.PulserTmean_Median)<1.5&&abs(R.PulserTrms.fElements-R.PulserTrms_Median)<0.2)");  // outlyerTrms
  //
  if (!isCalc){
    // make fits only once
    TStatToolkit toolkit;
    Double_t chi2=0;
    Int_t    npoints=0;
    TVectorD param;
    TMatrixD covar;
    tree->SetAlias("dt","PulserTmean.fElements-R.PulserTmean.fElements");
    TCut cutAll ="inPulser&&inPulserR"; 
    TString  fstringG="";              // global part
    fstringG+="ly.fElements++";
    fstringG+="(lx.fElements-134.)++";  
    for (Int_t isec=0; isec<72; isec++){
      TStatToolkit::FitPlane(tree,"dt", fstringG.Data(),Form("sector==%d",isec)+cutAll, chi2,npoints,param,covar,-1,0, 10000000, kFALSE);
      if (npoints<3) continue;
      printf("Setor=%d\n",isec);
      vec0[isec]=param[0];
      vecLy[isec]=param[1];
      vecLx[isec]=param[2];
      sec[isec]=isec;
      vecN[isec]=npoints;

      TStatToolkit::FitPlane(tree,"PulserTmean.fElements", fstringG.Data(),Form("sector==%d",isec)+cutAll, chi2,npoints,param,covar,-1,0, 10000000, kFALSE);
      if (npoints<3) continue;
      printf("Setor=%d\n",isec);
      vecA0[isec]=param[0];
      vecALy[isec]=param[1];
      vecALx[isec]=param[2];
      vecAChi2[isec]=TMath::Sqrt(chi2/npoints);
    }

    isCalc=kTRUE;
  }
  (*fPcstream)<<"dcs"<<     // Pulser information
    "PulserTSector.="<<&sec<<    // sector numbers
    //                      // fit in respect to reference
    "PulserTRef0.="<<&vec0<<    // offset change
    "PulserTRefY.="<<&vecLy<<   // slope y change - rad
    "PulserTRefX.="<<&vecLx<<   // slope x change - rad
    "PulserTRefChi2.="<<&vecChi2<<    // chi2 (rms in cm)
    "PulserTRefN.="<<&vecN<<     //number of accepted points
    //                       // fit in respect per mean per side
    "PulserT0.="<<&vecA0<<       // offset change
    "PulserTY.="<<&vecALy<<      // slope y change - rad
    "PulserTX.="<<&vecALx<<      // slope x change - rad
    "PulserTChi2.="<<&vecAChi2;  // chi2 (rms in cm)
}

void AliTPCcalibSummary::ProcessCurrent(Int_t irun, Int_t itime){
  //
  // Dump current 
  //
  //variables to export
  //
  static TObjArray *currentArray=new TObjArray(72);   // current graphs
  static TObjArray *currentArray2=new TObjArray(72);  // current graphs to export
  //
  static TVectorD currentIROC(36);                    // current snapshots
  static TVectorD currentOROC(36); 
  static TVectorF sector(72);                         //
  static Double_t medcurIROC = 0;
  static Double_t medcurOROC = 0;
  //
  static TVectorF minROC(72);                         // current mean +-5 minutes
  static TVectorF maxROC(72);
  static TVectorF meanROC(72);
  static TVectorF medianROC(72);
  static Double_t meanIIROC=0;
  static Double_t meanIOROC=0;
  static Double_t medianIIROC=0;
  static Double_t medianIOROC=0;
  //
  AliDCSSensorArray* voltageArray = AliTPCcalibDB::Instance()->GetVoltageSensors(irun);
  //
  for(Int_t j=1; j<36; j++) currentIROC[j-1] = fCalibDB->GetChamberHighVoltage(irun, j,itime,-1,kTRUE);
  for(Int_t j=36; j<72; j++) currentOROC[j-36] = fCalibDB->GetChamberHighVoltage(irun, j,itime,-1,kTRUE);
  medcurIROC = TMath::Median(36, currentIROC.GetMatrixArray());
  medcurOROC = TMath::Median(36, currentOROC.GetMatrixArray());


  if (currentArray->At(0)==0){
    for (Int_t isec=0; isec<72; isec++){
      TString sensorName="";
      const char* sideName=(isec%36<18) ? "A":"C";
      if (isec<36){
	//IROC
	sensorName=Form("TPC_ANODE_I_%s%02d_IMEAS",sideName,isec%18);
      }else{
	//OROC
	sensorName=Form("TPC_ANODE_O_%s%02d_0_IMEAS",sideName,isec%18);
      }      
    
      AliDCSSensor *sensor = 0;
      if (voltageArray) sensor= voltageArray->GetSensor(sensorName);   
      TGraph *gr=0;
      if (!sensor) gr=new TGraph(1);
      else{
	if (!sensor->GetGraph()) gr=new TGraph(1);
	else{
	  gr=sensor->GetGraph();
	  Double_t startTime=sensor->GetStartTime();
	  Double_t * time = new Double_t[gr->GetN()];
	  for (Int_t ip=0; ip<gr->GetN(); ip++){ time[ip]= (gr->GetX()[ip]*3600.)+startTime;}	  
	  gr=new TGraph(gr->GetN(), time, gr->GetY());	
	  delete [] time;
	}      
      }
      gr->Sort();
      currentArray->AddAt(gr, isec);
      currentArray->AddAt(gr->Clone(), isec);
    }
  }


  for (Int_t isec=0; isec<72; isec++){
    sector[isec]=isec;
    TGraph * gr = (TGraph*)currentArray->At(isec);
    TGraph * graph2 = (TGraph*)currentArray2->At(isec);    
    Int_t firstBin= TMath::BinarySearch(gr->GetN(), gr->GetX(), itime-300.)-2;
    Int_t lastBin= TMath::BinarySearch(gr->GetN(), gr->GetX(), itime+300.)+2;
    if (firstBin<0) firstBin=0;
    if (lastBin>=gr->GetN()) lastBin=gr->GetN()-1;
    //
    if (firstBin<lastBin){
      //
      minROC[isec]=TMath::MinElement(lastBin-firstBin, &(gr->GetY()[firstBin]));
      maxROC[isec]=TMath::MaxElement(lastBin-firstBin, &(gr->GetY()[firstBin]));
      meanROC[isec]=TMath::Mean(lastBin-firstBin, &(gr->GetY()[firstBin]));
      medianROC[isec]=TMath::Median(lastBin-firstBin, &(gr->GetY()[firstBin]));       
      graph2 = new TGraph(lastBin-firstBin, &(gr->GetX()[firstBin]), &(gr->GetY()[firstBin]));
      delete currentArray2->At(isec);
      currentArray2->AddAt(graph2,isec);
    }
    (*fPcstream)<<"dcs"<<     // current information
      Form("current%d.=",isec)<<graph2;
  }     
  meanIIROC=TMath::Mean(36, &(meanROC.GetMatrixArray()[0]));
  meanIOROC=TMath::Mean(36, &(meanROC.GetMatrixArray()[36]));
  medianIIROC=TMath::Median(36, &(meanROC.GetMatrixArray()[0]));
  medianIOROC=TMath::Median(36, &(meanROC.GetMatrixArray()[36]));
  //
  (*fPcstream)<<"dcs"<<     // current information
    "isec.="<<&sector<<                       //sector number
    "IIROC.="<<&currentIROC<<               // current sample at given moment
    "IOROC.="<<&currentOROC<<               // current sample at given moment
    "medianIIROC="<<medcurIROC<<            // median at given moment 
    "medianIOROC="<<medcurOROC<<            // median at given moment
    //
    "minIROC.="<<&minROC<<                  // minimum in +-5 min 
    "maxIROC.="<<&maxROC<<                  // maximum in +-5 min
    "meanIROC.="<<&meanROC<<                // mean in +-5 min
    "medianIROC.="<<&medianROC<<              // median in +-5 min
    "meanIIROC5="<<meanIIROC<<               // mean current in IROC +-5 minutes 
    "meanIOROC5="<<meanIOROC<<               // mean current in OROC 
    "medianIIROC5="<<medianIIROC<<           // median current in IROC 
    "medianIOROC5="<<medianIOROC;           // medianan current in OROC 
   

  (*fPcstream)<<"current"<<     // current information
    "time="<<itime<<
    "isec.="<<&sector<<                       //sector number
    "IIROC.="<<&currentIROC<<               // current sample at given moment
    "IOROC.="<<&currentOROC<<               // current sample at given moment
    "medianIIROC="<<medcurIROC<<            // median at given moment 
    "medianIOROC="<<medcurOROC<<            // median at given moment
    //
    "minIROC.="<<&minROC<<                  // minimum in +-5 min 
    "maxIROC.="<<&maxROC<<                  // maximum in +-5 min
    "meanIROC.="<<&meanROC<<                // mean in +-5 min
    "medianIROC.="<<&medianROC<<              // median in +-5 min
    "meanIIROC5="<<meanIIROC<<               // mean current in IROC +-5 minutes 
    "meanIOROC5="<<meanIOROC<<               // mean current in OROC 
    "medianIIROC5="<<medianIIROC<<           // median current in IROC 
    "medianIOROC5="<<medianIOROC<< // medianan current in OROC 
    "\n";

}




// TCanvas * DrawCEDiff(TTree * tree){
  
//   TCanvas *canvasIO = new TCanvas("canvasCEIO","canvasCEIO");
//   canvasIO->Divide(6,6);
//   for (Int_t isec=0; isec<36; isec++){
//     canvasIO->cd(isec+1);
//     dcs->Draw(Form("CET0.fElements[%d]-CET0.fElements[%d]",isec+36,isec),Form("abs(CETRef0.fElements[%d])<0.3",isec),"");
//     printf("%d\t%f\t%f\n",isec,dcs->GetHistogram()->GetMean(),dcs->GetHistogram()->GetRMS());
//   }

// }
