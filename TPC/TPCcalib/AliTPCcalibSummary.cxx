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
#include "AliLHCData.h"

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
  // ===| LHC data |===========================================================
  //
  //
  TVectorF vecMeanLHCBckgAlice(AliLHCData::kNBGs);
  GetAverageLHCData(vecMeanLHCBckgAlice);
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
    static TVectorF gainCorrIROCHVandPT(36);
    //
    static TVectorF voltagesOROC(36);
    static TVectorF voltagesOROCMedian(36);
    static TVectorF voltagesOROCNominal(36);
    static TVectorF voltagesOROCCurrentNominal(36);
    static TVectorF voltagesOROCStatus(36);
    static TVectorF voltagesOROCGoodFraction(36);
    static TVectorF gainCorrOROCHVandPT(36);

    for(Int_t j=0; j<36; j++){
      voltagesIROC[j]               = fCalibDB->GetChamberHighVoltage(irun, j,itime);
      voltagesIROCMedian[j]         = fCalibDB->GetChamberHighVoltageMedian(j);
      voltagesIROCNominal[j]        = fCalibDB->GetParameters()->GetNominalVoltage(j);
      voltagesIROCCurrentNominal[j] = fCalibDB->GetChamberCurrentNominalHighVoltage(j);
      voltagesIROCStatus[j]         = fCalibDB->GetChamberHVStatus(j);
      voltagesIROCGoodFraction[j]   = fCalibDB->GetChamberGoodHighVoltageFraction(j);
      gainCorrIROCHVandPT[j]        = fCalibDB->GetGainCorrectionHVandPT(itime, irun, j, 5, 1);
    }
    
    for(Int_t j=36; j<72; j++) {
      voltagesOROC[j-36]               = fCalibDB->GetChamberHighVoltage(irun, j,itime);
      voltagesOROCMedian[j-36]         = fCalibDB->GetChamberHighVoltageMedian(j);
      voltagesOROCNominal[j-36]        = fCalibDB->GetParameters()->GetNominalVoltage(j);
      voltagesOROCCurrentNominal[j-36] = fCalibDB->GetChamberCurrentNominalHighVoltage(j);
      voltagesOROCStatus[j-36]         = fCalibDB->GetChamberHVStatus(j);
      voltagesOROCGoodFraction[j-36]   = fCalibDB->GetChamberGoodHighVoltageFraction(j);
      gainCorrOROCHVandPT[j-36]        = fCalibDB->GetGainCorrectionHVandPT(itime, irun, j, 5, 1);
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
    // parameter description is as follows: // <name>;<unit>;<group>;<array type/array info>
    //      <name> should also include the unit; e.g. IROC anode voltage (V))
    // where group is e.g.
    //       o HV:              for High Voltage information
    //       o Environment:     for environmental information (temperature; pressure)
    //       o Pulser:          Calibartion pulser information
    //       o CE:              Central electrode inforamtion
    //       o Noise Pedestals: Noise and Pedestal information
    //       o ALTRO:           ALTRO configuration information
    //
    // <array type> describes the information stored in the elements of a TVectorT<> variable
    // this can e.g. be
    //       o Sector:   The array has 36 entries; one per sector (A00-A17 and C00 to C17)
    //       o Sector-A: The array has 18 entries; one per sector on the A-Side (A00-A17)
    //       o Sector-C: The array has 18 entries; one per sector on the C-Side (C00-C17)
    //       o ROC:      The array has 72 entries; one per ROC (IA00-IA17, IC00-IC17, OA00-OA17, OC00-OC17)
    //       o Parameters from a fit: in this case the description per parameter should be given
    (*fPcstream)<<"dcs"<<
      "run="<<irun<< // Run number
      "time="<<itime<< // Time stamp of calibration entry
      "startTimeGRP="<<startTimeGRP<< // Start time of run from GRP
      "stopTimeGRP="<<stopTimeGRP<< // Stop time of run from GRP
      "dataTakingActive="<<dataTakingActive<< // If data taking is active
      //run type
      "runType.="<<&runType<< // Run Type; e.g. PHYSICS; LASER; COSMIC; PEDESTAL; PULSER
      // voltage setting
      "VIROC.="               << &voltagesIROC               << // IROC anode voltage [calib interval] (V);HV;Sector
      "VIROCMedian.="         << &voltagesIROCMedian         << // IROC anode voltage [Median of run] (V);HV;Sector
      "VIROCNominal.="        << &voltagesIROCNominal        << // IROC anode voltage [global nominal] (V);HV;Sector
      "VIROCCurrentNominal.=" << &voltagesIROCCurrentNominal << // IROC anode voltage [current nominal] (V);HV;Sector
      "VIROCGoodHVFraction.=" << &voltagesIROCGoodFraction   << // IROC anode voltage [fraction of good settings];-;HV;Sector
      "VIROCStatus.="         << &voltagesIROCStatus         << // IROC HV status;-;HV;Sector
      "gainCorrIROCHVandPT.=" << &gainCorrIROCHVandPT        << // IROC gain correction factor using HV, P and T
      //
      "VOROC.="               << &voltagesOROC               << // OROC anode voltage [calib interval] (V);HV;Sector
      "VOROCMedian.="         << &voltagesOROCMedian         << // OROC anode voltage [Median of run] (V);HV;Sector
      "VOROCNominal.="        << &voltagesOROCNominal        << // OROC anode voltage [global nominal] (V);HV;Sector
      "VOROCCurrentNominal.=" << &voltagesOROCCurrentNominal << // OROC anode voltage [current nominal] (V);HV;Sector
      "VOROCGoodHVFraction.=" << &voltagesOROCGoodFraction   << // OROC anode voltage [fraction of good settings];-;HV;Sector
      "VOROCStatus.="         << &voltagesOROCStatus         << // OROC HV status;-;HV;Sector
      "gainCorrOROCHVandPT.=" << &gainCorrOROCHVandPT        << // OROC gain correction factor using HV, P and T
      //
      "medianVIROC="          << voltIROC                    << // IROC anode voltage [median of all IROCs] (V);HV
      "medianVOROC="          << voltOROC                    << // OROC anode voltage [median of all OROCs] (V);HV
      "coverIA="              << coverIA                     << // Cover voltage IROC A-Side (V);HV
      "coverIC="              << coverIC                     << // Cover voltage IROC C-Side (V);HV
      "coverOA="              << coverOA                     << // Cover voltage OROC A-Side (V);HV
      "coverOC="              << coverOC                     << // Cover voltage OROC C-Side (V);HV
      "skirtA="               << skirtA                      << // Skirt voltage IROC A-Side (V);HV
      "skirtC="               << skirtC                      << // Skirt voltage IROC C-Side (V);HV
      "ggOffA="               << ggOffA                      << // Gating grid offset voltage A-Side (V);HV
      "ggOffC="               << ggOffC                      << // Gating grid offset voltage C-Side (V);HV
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
      "temp00.="<<&vecTemp[0]<< // T-Fit ROC A-Side;Environment;Mean Temp (#circC);dT/dgx (K/cm);dT/dgy (K/cm)
      "temp10.="<<&vecTemp[1]<< // T-Fit OFC A-Side;Environment;Mean Temp (#circC);dT/dz (K/cm);dT/d#phi (K/rad)
      "temp20.="<<&vecTemp[2]<< // T-Fit IFC+TS A-Side;Environment;Mean Temp (#circC);dT/dz (K/cm);dT/d#phi (K/rad)
      "temp30.="<<&vecTemp[3]<< // T-Fit Skirt A-Side;Environment;Mean Temp (#circC);dT/dgx (K/cm);dT/dgy (K/cm)
      "temp40.="<<&vecTemp[4]<< // T-Fit IFC A-Side;Environment;Mean Temp (#circC);dT/dz (K/cm);dT/d#phi (K/rad)
      "temp01.="<<&vecTemp[5]<< // T-Fit ROC C-Side;Environment;Mean Temp (#circC);dT/dgx (K/cm);dT/dgy (K/cm)
      "temp11.="<<&vecTemp[6]<< // T-Fit OFC C-Side;Environment;Mean Temp (#circC);dT/dz (K/cm);dT/d#phi (K/rad)
      "temp21.="<<&vecTemp[7]<< // T-Fit IFC+TS C-Side;Environment;Mean Temp (#circC);dT/dz (K/cm);dT/d#phi (K/rad)
      "temp31.="<<&vecTemp[8]<< // T-Fit Skirt C-Side;Environment;Mean Temp (#circC);dT/dgx (K/cm);dT/dgy (K/cm)
      "temp41.="<<&vecTemp[9]<< // T-Fit IFC C-Side;Environment;Mean Temp (#circC);dT/dz (K/cm);dT/d#phi (K/rad)
      //
      "tempSkirtA.="<<&vecSkirtTempA<< // T Skirt A-Side;Environment;Sector-A
      "tempSkirtC.="<<&vecSkirtTempC;  // T Skirt C-Side;Environment;Sector-C

    ProcessGas(irun, itime);
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
    ProcessLHCData(irun, itime);


    (*fPcstream)<<"dcs"<<	
      //noise data
      "meanNoise.="    << &vNoiseMean           << // Mean Noise;Noise Pedestals;All Pads;IROCs;OROCs small pads;OROCs large pads
      "meanNoiseSen.=" << &vNoiseMeanSenRegions << // Mean Noise in sensitive regions;Noise Pedestals;All Pads;IROCs;OROCs small pads;OROCs large pads
      "rmsNoise.="     << &vNoiseRMS            << // RMS Noise;Noise Pedestals;All Pads;IROCs;OROCs small pads;OROCs large pads
      "rmsNoiseSen.="  << &vNoiseRMSSenRegions  << // RMS Noise in sensitive regions;Noise Pedestals;All Pads;IROCs;OROCs small pads;OROCs large pads
      "zeroNoise="     << nonMaskedZero         << // Pads with zero noise;Noise Pedestals
      "nNaN="          << nNaN                  << // Pads with NaN noise;Noise Pedestals
      //QA data
      "occQA.="  << &vQaOcc  <<
      "qQA.="    << &vQaQtot <<
      "qmaxQA.=" << &vQaQmax <<
      //pulser data
      "timePulser.=" << &vTimePulser <<
      "nOffPulser="<<nOffChannels<<
      //altro data
      "nMasked="<< nMasked << // Number of masked pads;ALTRO
      //
      //ce data -Jens version
      //
      "CEfitA.="        << &fitResultsA  << // CE-Fit A-Side;CE;Offset (timebins);IROC/OROC Offset (timebins);dt/dgx (timebins/cm);dt/dgy (timebins/cm);dt/dlx common (timebins/cm);dt/dlx IROCs (timebins/cm)
      "CEfitC.="        << &fitResultsC  << // CE Fit C-Side;CE;Offset (timebins);IROC/OROC Offset (timebins);dt/dgx (timebins/cm);dt/dgy (timebins/cm);dt/dlx common (timebins/cm);dt/dlx IROCs (timebins/cm)
      "nmaskedCE="      << nmaskedCE     << // CE Number of outliers;CE
      "chi2ACE="        << chi2ACE       << // CE-Fit Chi^{2} A-Side;CE
      "chi2CCE="        << chi2CCE       << // CE-Fit Chi^{2} C-Side;CE
      //
      //ce data new - MI version
      //
      "CEfitAMI.="      << &fitCEResultsA<< // CE-Fit A-Side [MI];CE;Offset (timebins);IROC/OROC Offset (timebins);dt/dgx (timebins/cm);dt/dgy (timebins/cm);dt/dlx common (timebins/cm);dt/dlx IROCs (timebins/cm)
      "CEfitCMI.="      << &fitCEResultsC<< // CE-Fit C-Side [MI];CE;Offset (timebins);IROC/OROC Offset (timebins);dt/dgx (timebins/cm);dt/dgy (timebins/cm);dt/dlx common (timebins/cm);dt/dlx IROCs (timebins/cm)
      "chi2CEA="        << chi2CEA       << // CE-Fit Chi^{2} A-Side [MI];CE
      "chi2CEC="        << chi2CEC       << // CE-Fit Chi^{2} C-Side [MI];CE
      //
      //ce graph data
      //
      "CEgrTEntries.="   << &vecTEntries << // CE-graph drift time - entries;CE;ROC
      "CEgrTMean.="      << &vecTMean    << // CE-graph mean drift time;CE;ROC
      "CEgrTRMS.="       << &vecTRMS     << // CE-graph RMS of drift time;CE;ROC
      "CEgrTMedian.="    << &vecTMedian  << // CE-graph median drift time;CE;ROC
      "CEgrQEntries.="   << &vecQEntries << // CE-graph charge - entries;CE;ROC
      "CEgrQMean.="      << &vecQMean    << // CE-graph mean charge;CE;ROC
      "CEgrQRMS.="       << &vecQRMS     << // CE-graph RMS charge;CE;ROC
      "CEgrQMedian.="    << &vecQMedian  << // CE-graph median charge;CE;ROC
      "CEgrDriftA="      << driftTimeA   << // CE median drift time A-Side;CE
      "CEgrDriftC="      << driftTimeC   << // CE median drift time C-Side;CE
      //
      //calib raw data
      //
      "nFailL1="         << nFailL1      << // RCU synchonisation failures;ALTRO
      // b field
      "Bz="              << bz           << // Magnetic Field (T);Environment
      "L3polarity="      << l3pol        << // L3 polarity;Environment
      // production information
      "nalien="          << nalien       << // obsolete
      "nRawAlien="       << nRawAlien    << // obsolete
      "nlocal="          << nlocal       << // obsolete
      "nRawLocal="       <<nRawLocal     << // obsolete
      //
      // comparisons with ref data
      //
      "pedestalDeviations.=" << &pedestalDeviations << // Pedestal variation to ref (fraction);Noise Pedestals;>#pm 0.5 ADC;>#pm 1 ADC;>#pm 1.5 ADC;>#pm 2.0 ADC
      "noiseDeviations.="    << &noiseDeviations    << // Noise var to ref (fraction);Noise Pedestals;>5%;>10%;>15%;>20%
      "pulserQdeviations.="  << &pulserQdeviations  << // Pulser-Q var to ref (fraction);Pulser;>0.5%;>1%;>5%;>10%
      //         "pulserVarQMean="<<varQMean<<
      "pulserNpadsOutOneTB=" << npadsOutOneTB       << // Number of pads with Pulser time var >#pm 1 tb to ROC mean;Pulser
      "pulserNpadsOffAdd="   << npadsOffAdd         << // Number of pads without signal but signal in ref;Pulser
      "driftCorrCosmAll="    << dvCorr              <<
      //
      // LHCData
      //
      "meanBckgAlice.="       << &vecMeanLHCBckgAlice <<
      "\n";
  }//end run loop
}



void AliTPCcalibSummary::ProcessGas(Int_t run, Int_t timeStamp)
{
  const Int_t nsensors=Int_t(AliTPCcalibDB::kNGasSensor);
  static TVectorF gasValues(nsensors);

  for (Int_t isen=0; isen<nsensors; ++isen) {
    gasValues(isen) = fCalibDB->GetGasSensorValue((AliTPCcalibDB::EDcsGasSensor)isen);
  }

  (*fPcstream)<<"dcs"<<
    "gasValues.=" << &gasValues;
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
  static Double_t     vdrift1=0;                                // TODO: repeated below, obsolete?
  vdrift1=fCalibDB->GetVDriftCorrectionTime(timeStamp,run,0,1); // TODO: repeated below, obsolete?

  if (array){
    laserA[0]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DELAY_LASER_ALL_A");
    laserA[1]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_A");
    laserA[2]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_A");
    laserC[0]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DELAY_LASER_ALL_C");
    laserC[1]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_C");
    laserC[2]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_C");
    cosmicAll =(TGraphErrors*)array->FindObject("TGRAPHERRORS_MEAN_VDRIFT_COSMICS_ALL");
  }

  //
  // TODO: the information stored in vlaserXX, vcosmicAll and vclaserXX,vccosmicAll 
  //       seem to be redundant information do we need to keep vlaserXX vcosmicAll
  //
  if (laserA[0]) vlaserA[0]= AliTPCcalibDButil::EvalGraphConst(laserA[0],timeStamp);
  if (laserA[1]) vlaserA[1]= AliTPCcalibDButil::EvalGraphConst(laserA[1],timeStamp);
  if (laserA[2]) vlaserA[2]= AliTPCcalibDButil::EvalGraphConst(laserA[2],timeStamp);
  if (laserC[0]) vlaserC[0]= AliTPCcalibDButil::EvalGraphConst(laserC[0],timeStamp);
  if (laserC[1]) vlaserC[1]= AliTPCcalibDButil::EvalGraphConst(laserC[1],timeStamp);
  if (laserC[2]) vlaserC[2]= AliTPCcalibDButil::EvalGraphConst(laserC[2],timeStamp);
  if (cosmicAll) vcosmicAll= AliTPCcalibDButil::EvalGraphConst(cosmicAll,timeStamp);
  (*fPcstream)<<"dcs"<<
    "vlaserA0="   << vlaserA[0] << // Laser offset A-Side;Drift            //TODO: Obsolete
    "vlaserA1="   << vlaserA[1] << // Laser drift correction A-Side;Drift  //TODO: Obsolete
    "vlaserA2="   << vlaserA[2] << // Laser gy correction A-Side;Drift     //TODO: Obsolete
    "vlaserC0="   << vlaserC[0] << // Laser offset C-Side;Drift            //TODO: Obsolete
    "vlaserC1="   << vlaserC[1] << // Laser drift correction C-Side;Drift  //TODO: Obsolete
    "vlaserC2="   << vlaserC[2] << // Laser gy correction C-Side;Drift     //TODO: Obsolete
    "vcosmicAll=" << vcosmicAll << // Cosmic drift corrrection;Drift       //TODO: Obsolete
    //
    "vdrift1="    << vdrift1;      // Combined drift correction ;Drift      // TODO: repeated below, obsolete?

  //
  // define distance to measurement
  //
  static Double_t dlaserA=0; 
  static Double_t dlaserC=0; 
  static Double_t dcosmic=0; 
  static Double_t slaserA=0; //TODO: Obsolete?
  static Double_t slaserC=0; //TODO: Obsolete?
  static Double_t scosmic=0; //TODO: Obsolete?
  static Double_t  vclaserA[3]={0,0,0};
  static Double_t  vclaserC[3]={0,0,0};
  static Double_t  vccosmicAll=0;
  for (Int_t i=0;i<3;i++){
    if (laserA[i]) AliTPCcalibDButil::GetNearest(laserA[i],timeStamp,dlaserA,vclaserA[i]);
    if (laserC[i]) AliTPCcalibDButil::GetNearest(laserC[i],timeStamp,dlaserC,vclaserC[i]);
  }  
  if (cosmicAll) AliTPCcalibDButil::GetNearest(cosmicAll,timeStamp,dcosmic,vccosmicAll);
  (*fPcstream)<<"dcs"<<
    "vclaserA0="   << vclaserA[0]<< // Laser offset A-Side;Drift
    "vclaserA1="   << vclaserA[1]<< // Laser drift correction A-Side;Drift
    "vclaserA2="   << vclaserA[2]<< // Laser gy correction A-Side;Drift
    "vclaserC0="   << vclaserC[0]<< // Laser offset A-Side;Drift
    "vclaserC1="   << vclaserC[1]<< // Laser drift correction A-Side;Drift
    "vclaserC2="   << vclaserC[2]<< // Laser gy correction A-Side;Drift
    "vccosmicAll=" << vccosmicAll<< // Cosmic drift corrrection;Drift
    "dlaserA="     << dlaserA    << // Distance to laser measurement A-Side;Drift
    "dlaserC="     << dlaserC    << // Distance to laser measurement C-Side;Drift
    "dcosmic="     << dcosmic    << // Distance to cosmics measurement A-Side;Drift
    "slaserA="     << slaserA    << //TODO: Obsolete?
    "slaserC="     << slaserC    << //TODO: Obsolete?
    "scosmic="     << scosmic;      //TODO: Obsolete?

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
  // TODO: are the values tdriftCE and tcdriftCE redundant and need both be kept?
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
    // TODO: use tdriftCE or tndriftCE
    if (i%36<18){
      vdriftCE[i] =(param->GetZLength(i)/(tdriftCE[i]*param->GetTSample()*(1.-ltime0A)-param->GetL1Delay()))/param->GetDriftV();
    }else{
      vdriftCE[i] =(param->GetZLength(i)/(tdriftCE[i]*param->GetTSample()*(1.-ltime0A)-param->GetL1Delay()))/param->GetDriftV();
    }
  }
  // export values
  (*fPcstream)<<"dcs"<<  
    "tdriftCE.="  << &tdriftCE  << // CE arrival time;CE;Sector                  // TODO: obsolete, redundant?
    "vdriftCE.="  << &vdriftCE  << // CE derived drift velocity;CE;Sector
    "tndriftCE.=" << &tndriftCE << // CE number of points;CE;Sector
    "tcdriftCE.=" << &tcdriftCE << // CE arrival time - nearest point;CE;Sector
    "tddriftCE.=" << &tddriftCE << // CE distance to closest measuement;CE;Sector
    "ltime0A="    << ltime0A    << // CE laser offset A-Side;CE
    "ltime0C="    << ltime0C     ; // CE laser offset C-Side;CE
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
  ltime0A  = fDButil->GetLaserTime0(run,timeStamp,36000,0); // TODO: not used, needed?
  ltime0C  = fDButil->GetLaserTime0(run,timeStamp,36000,1); // TODO: not used, needed?

  (*fPcstream)<<"dcs"<<  
    //
    "vdriftCEA="   << vdriftCEA   << // CE drift correction A-Side;Drift
    "vdriftCEC="   << vdriftCEC   << // CE drift correction C-Side;Drift
    "vdriftCEM="   << vdriftCEM   << // CE drift correction Mean;Drift
    "dcea="        << dcea        << // CE distance to closest measurement A-Side;Drift
    "dcec="        << dcec        << // CE distance to closest measurement C-Side;Drift
    "dcem="        << dcem        << // CE distance to closest measurement Mean;Drift
    "vdriftLTA="   << vdriftLTA   << // Offline Laser track vdrift correction A-Side;Drift
    "vdriftLTC="   << vdriftLTC   << // Offline Laser track vdrift correction C-Side;Drift
    "vdriftLTM="   << vdriftLTM   << // Offline Laser track vdrift correction Mean;Drift
    "dla="         << dla         << // Offline Laser track distance to closest measurement A-Side;Drift
    "dlc="         << dlc         << // Offline Laser track distance to closest measurement C-Side;Drift
    "dlm="         << dlm         << // Offline Laser track distance to closest measurement Mean;Drift
    "vdriftLTAon=" << vdriftLTAon << // Online Laser track vdrift correction A-Side;Drift
    "vdriftLTCon=" << vdriftLTCon << // Online Laser track vdrift correction C-Side;Drift
    "vdriftLTMon=" << vdriftLTMon << // Online Laser track vdrift correction Mean;Drift
    "dlaOn="       << dlaon       << // Online Laser track distance to closest measurement A-Side;Drift
    "dlcOn="       << dlcon       << // Online Laser track distance to closest measurement C-Side;Drift
    "dlmOn="       << dlmon       << // Online Laser track distance to closest measurement Mean;Drift
    //
    //
    "vdriftITS="   << vdriftITS   << // TPC-ITS vdrift correction;Drift
    "dits="        << dits        << // TPC-ITS vdrift correction distance to closest measurement;Drift
    "ctime0="      << ctime0      << // Trigger offset correction;Drift
    "vdriftP="     << vdriftP     << // Cosmics vdrift correction;Drift
    "dp="          << dp          << // Cosmics vdrift correction distance to closest measurement;Drift
    "vdrift1="     << vdrift1      ; // combined drift velocity;Drift

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
  //
  static TVectorD vGainQMaxGraphRegion(3);
  static TVectorD vGainQTotGraphRegion(3);
  //
  static TGraphErrors ggrPadEqualMax(36);
  static TGraphErrors ggrPadEqualTot(36);
  //
  static TGraphErrors ggrDipAngleMaxShort;
  static TGraphErrors ggrDipAngleMaxMedium;
  static TGraphErrors ggrDipAngleMaxLong;
  static TGraphErrors ggrDipAngleMaxAbsolute;
  //
  static TGraphErrors ggrDipAngleTotShort;
  static TGraphErrors ggrDipAngleTotMedium;
  static TGraphErrors ggrDipAngleTotLong;
  static TGraphErrors ggrDipAngleTotAbsolute;
  //
  static TGraphErrors ggrMultiplicityTot;
  static TGraphErrors ggrMultiplicityMax;
  //
  static TVectorD vFitDipAngleParMaxShort(3);
  static TVectorD vFitDipAngleParMaxMedium(3);
  static TVectorD vFitDipAngleParMaxLong(3);
  static TVectorD vFitDipAngleParMaxAbsolute(3);
  //
  static TVectorD vFitDipAngleParTotShort(3);
  static TVectorD vFitDipAngleParTotMedium(3);
  static TVectorD vFitDipAngleParTotLong(3);
  static TVectorD vFitDipAngleParTotAbsolute(3);
  static Int_t gainMIPTime0=0;
  static Int_t gainMIPTime1=0;
   
  vGainGraphIROC.Zero();
  vGainGraphOROCmed.Zero();
  vGainGraphOROClong.Zero();
  vGainGraphIROCErr.Zero();
  vGainGraphOROCmedErr.Zero();
  vGainGraphOROClongErr.Zero();
  vGainQMaxGraphRegion.Zero();
  vGainQTotGraphRegion.Zero();
  TGraphErrors grDummy;
  TObjArray * gainSplines = fCalibDB->GetTimeGainSplinesRun(irun);
  if (gainSplines) {
    TGraphErrors * graphMIP = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL");
    TGraphErrors * graphCosmic = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL");
    TGraphErrors * graphAttach = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_ATTACHMENT_BEAM_ALL");
    //
    TGraphErrors * grPadEqualQMax = (TGraphErrors * ) gainSplines->FindObject("TGRAPHERRORS_MEANQMAX_PADREGIONGAIN_BEAM_ALL");
    TGraphErrors * grPadEqualQTot = (TGraphErrors * ) gainSplines->FindObject("TGRAPHERRORS_MEANQTOT_PADREGIONGAIN_BEAM_ALL");
    //
    TGraphErrors * graphGainIROC       = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_CHAMBERGAIN_SHORT_BEAM_ALL");
    TGraphErrors * graphGainOROCMedium = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_CHAMBERGAIN_MEDIUM_BEAM_ALL");
    TGraphErrors * graphGainOROCLong   = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEAN_CHAMBERGAIN_LONG_BEAM_ALL");
    //
    //
    TF1*  funDipAngleMax[4]={0x0,0x0,0x0,0x0};
    TF1*  funDipAngleTot[4]={0x0,0x0,0x0,0x0};
    TGraphErrors*  grDipAngleMax[4]={0x0,0x0,0x0,0x0};
    TGraphErrors*  grDipAngleTot[4]={0x0,0x0,0x0,0x0};
    const char* names[4]={"SHORT","MEDIUM","LONG","ABSOLUTE"};
    for (Int_t iPadRegion=0; iPadRegion<4; ++iPadRegion) {
      funDipAngleMax[iPadRegion]=(TF1*) gainSplines->FindObject(Form("TF1_QMAX_DIPANGLE_%s_BEAM_ALL",names[iPadRegion]));
      funDipAngleTot[iPadRegion]=(TF1*) gainSplines->FindObject(Form("TF1_QTOT_DIPANGLE_%s_BEAM_ALL",names[iPadRegion]));
      grDipAngleMax[iPadRegion]= (TGraphErrors*) gainSplines->FindObject(Form("TGRAPHERRORS_QMAX_DIPANGLE_%s_BEAM_ALL",names[iPadRegion]));
      grDipAngleTot[iPadRegion]= (TGraphErrors*) gainSplines->FindObject(Form("TGRAPHERRORS_QTOT_DIPANGLE_%s_BEAM_ALL",names[iPadRegion]));
    }
    //
    for(Int_t iPar=0; iPar < 3; iPar++) {
      if (funDipAngleMax[0]) vFitDipAngleParMaxShort(iPar)    = funDipAngleMax[0]->GetParameter(iPar);
      if (funDipAngleMax[1]) vFitDipAngleParMaxMedium(iPar)   = funDipAngleMax[1]->GetParameter(iPar);
      if (funDipAngleMax[2]) vFitDipAngleParMaxLong(iPar)     = funDipAngleMax[2]->GetParameter(iPar);
      if (funDipAngleMax[3]) vFitDipAngleParMaxAbsolute(iPar) = funDipAngleMax[3]->GetParameter(iPar);
      //
      if (funDipAngleTot[0]) vFitDipAngleParTotShort(iPar)    = funDipAngleTot[0]->GetParameter(iPar);
      if (funDipAngleTot[1]) vFitDipAngleParTotMedium(iPar)   = funDipAngleTot[1]->GetParameter(iPar);
      if (funDipAngleTot[2]) vFitDipAngleParTotLong(iPar)     = funDipAngleTot[2]->GetParameter(iPar);
      if (funDipAngleTot[3]) vFitDipAngleParTotAbsolute(iPar) = funDipAngleTot[3]->GetParameter(iPar);
    }
    //
    if (grDipAngleMax[0]) ggrDipAngleMaxShort    = * grDipAngleMax[0];
    if (grDipAngleMax[1]) ggrDipAngleMaxMedium   = * grDipAngleMax[1];
    if (grDipAngleMax[2]) ggrDipAngleMaxLong     = * grDipAngleMax[2];
    if (grDipAngleMax[3]) ggrDipAngleMaxAbsolute = * grDipAngleMax[3];
    //
    if (grDipAngleTot[0]) ggrDipAngleTotShort    = * grDipAngleTot[0];
    if (grDipAngleTot[1]) ggrDipAngleTotMedium   = * grDipAngleTot[1];
    if (grDipAngleTot[2]) ggrDipAngleTotLong     = * grDipAngleTot[2];
    if (grDipAngleTot[3]) ggrDipAngleTotAbsolute = * grDipAngleTot[3];
    //
    //
    TGraphErrors *grPadEqualMax = (TGraphErrors * ) gainSplines->FindObject("TGRAPHERRORS_MEANQMAX_PADREGIONGAIN_BEAM_ALL");
    TGraphErrors *grPadEqualTot = (TGraphErrors * ) gainSplines->FindObject("TGRAPHERRORS_MEANQTOT_PADREGIONGAIN_BEAM_ALL");
    if (grPadEqualMax) ggrPadEqualMax = *grPadEqualMax;
    if (grPadEqualTot) ggrPadEqualTot = *grPadEqualTot;
    //
    TGraphErrors * grMultiplicityTot  = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEANQTOT_MULTIPLICITYDEPENDENCE_BEAM_ALL");
    TGraphErrors * grMultiplicityMax  = (TGraphErrors *) gainSplines->FindObject("TGRAPHERRORS_MEANQMAX_MULTIPLICITYDEPENDENCE_BEAM_ALL");
    if (grMultiplicityTot) ggrMultiplicityTot = *grMultiplicityTot;
    if (grMultiplicityMax) ggrMultiplicityTot = *grMultiplicityMax;

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
       for (Int_t i=0; i<3; ++i){
	 vGainQMaxGraphRegion[i]=grPadEqualQMax->GetY()[i];
	 vGainQTotGraphRegion[i]=grPadEqualQTot->GetY()[i];
       }
    }
    
    if (graphMIP) {
      gainMIP = AliTPCcalibDButil::EvalGraphConst(graphMIP,timeStamp);
      gainMIPTime0=graphMIP->GetX()[0];
      gainMIPTime1=graphMIP->GetX()[graphMIP->GetN()-1];
    }
    if (graphCosmic) gainCosmic = AliTPCcalibDButil::EvalGraphConst(graphCosmic,timeStamp);
    if (graphAttach) attachMIP = AliTPCcalibDButil::EvalGraphConst(graphAttach,timeStamp);
    if (graphMIP)  AliTPCcalibDButil::GetNearest(graphMIP, timeStamp, dMIP,dummy);
  }
    
  // time dependence of gain 
  (*fPcstream)<<"dcs"<<
    "grPadEqualMax.="             << &ggrPadEqualMax             <<
    "grPadEqualTot.="             << &ggrPadEqualTot             <<
    "rocGainIROC.="               << &vGainGraphIROC             <<
    "rocGainOROCMedium.="         << &vGainGraphOROCmed          <<
    "rocGainOROCLong.="           << &vGainGraphOROClong         <<
    "rocGainErrIROC.="            << &vGainGraphIROCErr          <<
    "rocGainErrOROCMedium.="      << &vGainGraphOROCmedErr       <<
    "rocGainErrOROCLong.="        << &vGainGraphOROClongErr      <<
    "vGainQMaxGraphRegion.="      << &vGainQMaxGraphRegion       <<
    "vGainQTotGraphRegion.="      << &vGainQTotGraphRegion       <<
    //
    "vFitDipAngleParMaxShort.="   << &vFitDipAngleParMaxShort    <<
    "vFitDipAngleParMaxMedium.="  << &vFitDipAngleParMaxMedium   <<
    "vFitDipAngleParMaxLong.="    << &vFitDipAngleParMaxLong     <<
    "vFitDipAngleParMaxAbsolute.="<< &vFitDipAngleParMaxAbsolute <<
    //
    "vFitDipAngleParTotShort.="   << &vFitDipAngleParTotShort    <<
    "vFitDipAngleParTotMedium.="  << &vFitDipAngleParTotMedium   <<
    "vFitDipAngleParTotLong.="    << &vFitDipAngleParTotLong     <<
    "vFitDipAngleParTotAbsolute.="<< &vFitDipAngleParTotAbsolute <<
    //
    "grDipAngleMaxShort.="        << &ggrDipAngleMaxShort        <<
    "grDipAngleMaxMedium.="       << &ggrDipAngleMaxMedium       <<
    "grDipAngleMaxLong.="         << &ggrDipAngleMaxLong         <<
    "grDipAngleMaxAbsolute.="     << &ggrDipAngleMaxAbsolute     <<
    //
    "grDipAngleTotShort.="        << &ggrDipAngleTotShort        <<
    "grDipAngleTotMedium.="       << &ggrDipAngleTotMedium       <<
    "grDipAngleTotLong.="         << &ggrDipAngleTotLong         <<
    "grDipAngleTotAbsolute.="     << &ggrDipAngleTotAbsolute     <<
    //
    "grMultiplicityTot.="         << &ggrMultiplicityTot         <<
    "grMultiplicityMax.="         << &ggrMultiplicityMax         <<
    //
    "gainMIP="                    << gainMIP                     <<   // gain normalization parameter
    "gainMIPTime0="               << gainMIPTime0                <<   // first bin for time calibration
    "gainMIPTime1="               << gainMIPTime1                <<   // last bin for gain calibration
    "attachMIP="                  << attachMIP                   <<
    "dMIP="                       << dMIP                        <<
    "gainCosmic="                 << gainCosmic                   ;
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
    "CETSector.="  << &sec      << // sector numbers
    "CETRefA.="    << &vecASide << // diff to reference A side
    "CETRefC.="    << &vecCSide << // diff to reference C side
    //                      // fit in respect to reference data
    "CETRef0.="    << &vec0     << // offset change
    "CETRefY.="    << &vecLy    << // slope y change - rad
    "CETRefX.="    << &vecLx    << // slope x change - rad
    "CETRefChi2.=" << &vecChi2  << // chi2 (rms in cm)
    "CETRefN.="    << &vecN     << //number of accepted points
    //                       // fit in respect per mean per side
    "CET0.="       << &vecA0    << // offset change
    "CETY.="       << &vecALy   << // slope y change - rad
    "CETX.="       << &vecALx   << // slope x change - rad
    "CETChi2.="    << &vecAChi2  ; // chi2 (rms in cm)
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
    "PulserTSector.="  << &sec     << // sector numbers
    //                      // fit in respect to reference
    "PulserTRef0.="    << &vec0    << // offset change
    "PulserTRefY.="    << &vecLy   << // slope y change - rad
    "PulserTRefX.="    << &vecLx   << // slope x change - rad
    "PulserTRefChi2.=" << &vecChi2 << // chi2 (rms in cm)
    "PulserTRefN.="    << &vecN    << //number of accepted points
    //                       // fit in respect per mean per side
    "PulserT0.="       << &vecA0   << // offset change
    "PulserTY.="       << &vecALy  << // slope y change - rad
    "PulserTX.="       << &vecALx  << // slope x change - rad
    "PulserTChi2.="    << &vecAChi2 ; // chi2 (rms in cm)
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
    "isec.="        << &sector      << // sector number
    "IIROC.="       << &currentIROC << // current sample at given moment
    "IOROC.="       << &currentOROC << // current sample at given moment
    "medianIIROC="  << medcurIROC   << // median at given moment
    "medianIOROC="  << medcurOROC   << // median at given moment
    //
    "minIROC.="     << &minROC      << // minimum in +-5 min
    "maxIROC.="     << &maxROC      << // maximum in +-5 min
    "meanIROC.="    << &meanROC     << // mean in +-5 min
    "medianIROC.="  << &medianROC   << // median in +-5 min
    "meanIIROC5="   << meanIIROC    << // mean current in IROC +-5 minutes
    "meanIOROC5="   << meanIOROC    << // mean current in OROC
    "medianIIROC5=" << medianIIROC  << // median current in IROC
    "medianIOROC5=" << medianIOROC   ; // medianan current in OROC
   

  (*fPcstream)<<"current"<<     // current information
    "time="         << itime        <<
    "isec.="        << &sector      << // sector number
    "IIROC.="       << &currentIROC << // current sample at given moment
    "IOROC.="       << &currentOROC << // current sample at given moment
    "medianIIROC="  << medcurIROC   << // median at given moment
    "medianIOROC="  << medcurOROC   << // median at given moment
    //
    "minIROC.="     << &minROC      << // minimum in +-5 min
    "maxIROC.="     << &maxROC      << // maximum in +-5 min
    "meanIROC.="    << &meanROC     << // mean in +-5 min
    "medianIROC.="  << &medianROC   << // median in +-5 min
    "meanIIROC5="   << meanIIROC    << // mean current in IROC +-5 minutes
    "meanIOROC5="   << meanIOROC    << // mean current in OROC
    "medianIIROC5=" << medianIIROC  << // median current in IROC
    "medianIOROC5=" << medianIOROC  << // medianan current in OROC
    "\n";

}

void AliTPCcalibSummary::ProcessLHCData(Int_t irun, Int_t itime)
{
  // itime
  static TVectorF vecLHCBckgAlice(AliLHCData::kNBGs);

  // ---| fill data if lhcData object exists |---
  const AliLHCData *lhcData=GetLHCdata();

  for (Int_t i=0; i<Int_t(AliLHCData::kNBGs); ++i) {
    vecLHCBckgAlice[i]=lhcData?lhcData->GetBckgAlice(i,Double_t(itime)):-1.f;
  }

  (*fPcstream)<<"dcs"<<     // current information
  "bckgAlice.=" << &vecLHCBckgAlice;

}

void AliTPCcalibSummary::AddMetadataRawQA(TTree * treeRawQATPC){
  //
  // Make Aliases and description for Raw QA trending
  //
  treeRawQATPC->SetAlias("isIROC","(Iteration$<36)");
  treeRawQATPC->SetAlias("isOROC","(Iteration$>=36)");
  treeRawQATPC->SetAlias("occOK0","(occQA.fElements>0)");
  treeRawQATPC->SetAlias("occIROC","Sum$(occQA.fElements*(isIROC*occOK0))/Sum$(isIROC*occOK0)");
  treeRawQATPC->SetAlias("occOROC","Sum$(occQA.fElements*(isOROC*occOK0))/Sum$(isOROC*occOK0)");

  TStatToolkit::AddMetadata(treeRawQATPC,"occIROC.AxisTitle","Occupancy IROC ");
  TStatToolkit::AddMetadata(treeRawQATPC,"occIROC.Title","Sum$(occQA.fElements*(isIROC*occOK0))/Sum$(isIROC*occOK0)");
  TStatToolkit::AddMetadata(treeRawQATPC,"occIROC.Legend","IROC occ.");
  TStatToolkit::AddMetadata(treeRawQATPC,"occIROC.Comment","Digits occupancy in IROC  #(A>thr)/#All as obtained in AMORE QA");
  TStatToolkit::AddMetadata(treeRawQATPC,"occOROC.AxisTitle","Occupancy IROC ");
  TStatToolkit::AddMetadata(treeRawQATPC,"occOROC.Title","Sum$(occQA.fElements*(isOROC*occOK0))/Sum$(isOROC*occOK0)");
  TStatToolkit::AddMetadata(treeRawQATPC,"occOROC.Legend","OROC occ.");
  TStatToolkit::AddMetadata(treeRawQATPC,"occOROC.Comment","Digits occupanncy in OROC  #(A>thr)/#All as obtained in AMORE QA");
}

void AliTPCcalibSummary::AddMetadataGain(TTree * treeRawQATPC){
  //
  // Define aliases and valriable description for some important gain calibration parameters
  //
  treeRawQATPC->SetAlias("gainDefined","gainMIP!=0&&abs(dits)<1800");
  TStatToolkit::AddMetadata(treeRawQATPC,"gainMIP.AxisTitle","gain conversion (50/MIP)");
  TStatToolkit::AddMetadata(treeRawQATPC,"gainMIP.Title","gainMIP");
  TStatToolkit::AddMetadata(treeRawQATPC,"gainMIP.Legend","gain(t)");
  TStatToolkit::AddMetadata(treeRawQATPC,"gainMIP.Comment","Gain normalization coeficient used for time dependenta gain calibation. Last correction to move combined dEdx for MIP to channel 50");

}


void AliTPCcalibSummary::AddMetadata(TTree * tree){
  //
  // add metadata for automatic documentiation 
  //
  AddMetadataRawQA(tree);
  AddMetadataGain(tree);
  
}

AliLHCData* AliTPCcalibSummary::GetLHCdata()
{
  static AliCDBEntry *cdbLHCData=0x0;
  static AliLHCData *lhcData=0x0;

  AliCDBManager *man=AliCDBManager::Instance();

  // ===| get LHCData OCDB object |============================================

  // ===| check if we already have the latest object |=========================
  const Int_t currentRun=man->GetRun();
  if (cdbLHCData && cdbLHCData->GetId().GetAliCDBRunRange().Comprises(AliCDBRunRange(currentRun, currentRun))) {
    return lhcData;
  }

  // ---| if not, load the object |---------------------------------------------
  cdbLHCData=man->Get("GRP/GRP/LHCData");

  if (!cdbLHCData) {
    AliError("Could not get LHCData");
    return lhcData;
  }

  lhcData=static_cast<AliLHCData*>(cdbLHCData->GetObject());
  return lhcData;
}

void AliTPCcalibSummary::GetAverageLHCData(TVectorF &valsBckgAlice)
{
  // average Alice background
  const AliLHCData *lhcData=GetLHCdata();
  valsBckgAlice.ResizeTo(AliLHCData::kNBGs);

  // ===| reset to -1 |---
  for (Int_t i=0; i<Int_t(AliLHCData::kNBGs); ++i) {
    valsBckgAlice[i]=-1.f;
  }

  // ===| exit if lhcData does not exist |===
  if (!lhcData) {
    return;
  }

  // ===| build simple average per background estimator |===
  for (Int_t ibckg=0; ibckg<Int_t(AliLHCData::kNBGs); ++ibckg) {
    const Int_t nValues=lhcData->GetNBckgAlice(ibckg);
    if (nValues<=0) { continue; }
    Double_t sum=0.;
    for (Int_t ivalue=0; ivalue<nValues; ++ivalue ) {
      sum+=lhcData->GetBckgAlice(ibckg,ivalue)->GetValue();
    }
    sum/=Double_t(nValues);
    valsBckgAlice[ibckg]=sum;
  }

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
