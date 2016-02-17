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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// 	************** Class for ZDC reconstruction      **************      //
//                  Author: Chiara.Oppedisano@to.infn.it		     //
//                                                                           //
// NOTATIONS ADOPTED TO IDENTIFY DETECTORS (used in different ages!):	     //
//   (ZNC,ZPC) or (ZNC, ZPC) or RIGHT refers to side C (RB26)		     //
//   (ZNA,ZPA) or (ZNA, ZPA) or LEFT refers to side A (RB24)		     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TH2F.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TMap.h>

#include "AliRawReader.h"
#include "AliESDEvent.h"
#include "AliESDZDC.h"
#include "AliZDCDigit.h"
#include "AliZDCRawStream.h"
#include "AliZDCReco.h"
#include "AliZDCReconstructor.h"
#include "AliZDCPedestals.h"
#include "AliZDCEnCalib.h"
#include "AliZDCSaturationCalib.h"
#include "AliZDCTowerCalib.h"
#include "AliZDCMBCalib.h"
#include "AliZDCTDCCalib.h"
#include "AliZDCChMap.h"
#include "AliZDCRecoParam.h"
#include "AliZDCRecoParampp.h"
#include "AliZDCRecoParamPbPb.h"
#include "AliRunInfo.h"
#include "AliLHCClockPhase.h"


ClassImp(AliZDCReconstructor)
AliZDCRecoParam *AliZDCReconstructor::fgRecoParam=0;  //reconstruction parameters
AliZDCMBCalib *AliZDCReconstructor::fgMBCalibData=0;  //calibration parameters for A-A reconstruction

//_____________________________________________________________________________
AliZDCReconstructor:: AliZDCReconstructor() :
  fPedData(GetPedestalData()),
  fEnCalibData(GetEnergyCalibData()),
  fSatCalibData(GetSaturationCalibData()),
  fTowCalibData(GetTowerCalibData()),
  fTDCCalibData(GetTDCCalibData()),
  fMapping(GetMapping()),
  fRecoMode(0),
  fBeamEnergy(0.),
  fNRun(0),
  fIsCalibrationMB(kFALSE),
  fSignalThreshold(7),
  fMeanPhase(0),
  fESDZDC(NULL)
{
  // **** Default constructor
}


//_____________________________________________________________________________
AliZDCReconstructor::~AliZDCReconstructor()
{
// destructor
//   if(fgRecoParam)    delete fgRecoParam;
  AliCDBManager * man = AliCDBManager::Instance();
  if (man && !man->GetCacheFlag()) { // CDB objects must NOT be deleted if cache is active!
    if(fPedData)      delete fPedData;    
    if(fEnCalibData)  delete fEnCalibData;
    if(fSatCalibData) delete fSatCalibData;
    if(fTowCalibData) delete fTowCalibData;
    if(fgMBCalibData) delete fgMBCalibData;
    if(fTDCCalibData) delete fTDCCalibData;
    if(fMapping)      delete fMapping;
  }
  if(fESDZDC)       delete fESDZDC;
}

//____________________________________________________________________________
void AliZDCReconstructor::Init()
{
  // Setting reconstruction parameters
    
  TString runType = GetRunInfo()->GetRunType();
  if((runType.CompareTo("CALIBRATION_MB")) == 0){
    fIsCalibrationMB = kTRUE;
  }
    
  TString beamType = GetRunInfo()->GetBeamType();
  // To allow reconstruction in tests without beam
  if(((beamType.CompareTo("UNKNOWN"))==0) && 
     ((runType.CompareTo("PHYSICS"))==0 || (runType.CompareTo("CALIBRATION_BC"))==0)){
    fRecoMode=1;
  }
    
  fBeamEnergy = GetRunInfo()->GetBeamEnergy();
  if(fBeamEnergy<0.01){
     AliWarning(" Beam energy value missing -> setting it to 1380 GeV ");
     fBeamEnergy = 1380.;
  }
  
  if(((beamType.CompareTo("pp"))==0) || ((beamType.CompareTo("p-p"))==0)
     ||((beamType.CompareTo("PP"))==0) || ((beamType.CompareTo("P-P"))==0)){
    fRecoMode=1;
  }
  else if(((beamType.CompareTo("p-A"))==0) || ((beamType.CompareTo("A-p"))==0)
     ||((beamType.CompareTo("P-A"))==0) || ((beamType.CompareTo("A-P"))==0)){
    fRecoMode=1;
  }
  else if((beamType.CompareTo("A-A")) == 0 || (beamType.CompareTo("AA")) == 0){
    fRecoMode=2;
    /*if(!fgRecoParam) fgRecoParam = const_cast<AliZDCRecoParam*>(GetRecoParam());
    if(fgRecoParam){
      fgRecoParam->SetGlauberMCDist(fBeamEnergy);	
    }*/ 
  }

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase"); 
  if (!entry) AliFatal("LHC clock-phase shift is not found in OCDB !");
  else{
    AliLHCClockPhase *phaseLHC = (AliLHCClockPhase*)entry->GetObject();
    // 4/2/2011 According to A. Di Mauro BEAM1 measurement is more reliable 
    // than BEAM2 and therefore also than the average of the 2
    fMeanPhase = phaseLHC->GetMeanPhaseB1();
  }  
  if(fIsCalibrationMB==kFALSE)  
    AliInfo(Form("\n\n ***** ZDC reconstruction initialized for %s @ %1.0f + %1.0f GeV *****\n\n",
    	beamType.Data(), fBeamEnergy, fBeamEnergy));
  
  // if EMD calibration run NO ENERGY CALIBRATION should be performed
  // pp-like reconstruction must be performed (E calib. coeff. = 1)
  if((runType.CompareTo("CALIBRATION_EMD")) == 0) fRecoMode=1; 
  
  AliInfo(Form("  ZDC reconstruction mode %d (1 -> p-p/p-A, 2-> A-A)\n\n",fRecoMode));
  
  fESDZDC = new AliESDZDC();
}


//____________________________________________________________________________
void AliZDCReconstructor::Init(TString beamType, Float_t beamEnergy)
{
  // Setting reconstruction mode
  // Needed to work in the HLT framework
  
  fIsCalibrationMB = kFALSE;
     
  fBeamEnergy = beamEnergy;
  
  if(((beamType.CompareTo("pp"))==0) || ((beamType.CompareTo("p-p"))==0)
     ||((beamType.CompareTo("PP"))==0) || ((beamType.CompareTo("P-P"))==0)){
    fRecoMode=1;
  }
  else if(((beamType.CompareTo("p-A"))==0) || ((beamType.CompareTo("A-p"))==0)
     ||((beamType.CompareTo("P-A"))==0) || ((beamType.CompareTo("A-P"))==0)){
    fRecoMode=1;
  }
  else if((beamType.CompareTo("A-A")) == 0 || (beamType.CompareTo("AA")) == 0){
    fRecoMode=2;
    //if(!fgRecoParam) fgRecoParam = const_cast<AliZDCRecoParam*>(GetRecoParam());
    //if( fgRecoParam ) fgRecoParam->SetGlauberMCDist(fBeamEnergy);	
  }    

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase"); 
  if (!entry) AliFatal("LHC clock-phase shift is not found in OCDB !");
  else{
    AliLHCClockPhase *phaseLHC = (AliLHCClockPhase*)entry->GetObject();
    fMeanPhase = phaseLHC->GetMeanPhase();
  }
  fESDZDC = new AliESDZDC();
  
  AliInfo(Form("\n\n ***** ZDC reconstruction initialized for %s @ %1.0f + %1.0f GeV *****\n\n",
    	beamType.Data(), fBeamEnergy, fBeamEnergy));
  
}

//_____________________________________________________________________________
void AliZDCReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  // *** Local ZDC reconstruction for digits
  // Works on the current event
    
  // Retrieving calibration data  
  // Parameters for mean value pedestal subtraction
  int const kNch = 24;
  Float_t meanPed[2*kNch];    
  for(Int_t jj=0; jj<2*kNch; jj++) meanPed[jj] = fPedData->GetMeanPed(jj);
  // Parameters pedestal subtraction through correlation with out-of-time signals
  Float_t corrCoeff0[2*kNch], corrCoeff1[2*kNch];
  for(Int_t jj=0; jj<2*kNch; jj++){
     corrCoeff0[jj] = fPedData->GetPedCorrCoeff0(jj);
     corrCoeff1[jj] = fPedData->GetPedCorrCoeff1(jj);
  }
  Bool_t testPedSubBit = fPedData->TestPedModeBit();
  
  Bool_t chPedSubMode[kNch] = {0,};
  if(testPedSubBit){
       for(int i=0; i<kNch; i++) chPedSubMode[i] = fPedData->GetUseCorrFit(i);
  }

  int adcInTime[kNch][2], adcOutOfTime[kNch][2];
  float adcCorr[kNch][2];
  int signalCode[2*kNch]={0,};
  for(int ich=0; ich<24; ich++){
     for(int ig=0; ig<=1; ig++){
        adcInTime[ich][ig] = adcOutOfTime[ich][ig] = adcCorr[ich][ig] = -1;
     }
  }
  int tdcCabling[7]={-1,-1,-1,-1,-1,-1,-1};

  // get digits
  AliZDCDigit digit;
  AliZDCDigit* pdigit = &digit;
  digitsTree->SetBranchAddress("ZDC", &pdigit);
  //printf("\n\t # of digits in tree: %d\n",(Int_t) digitsTree->GetEntries());

  // loop over digits  
  Int_t digNentries = digitsTree->GetEntries();
  for(Int_t iDigit=0; iDigit<digNentries; iDigit++){
     digitsTree->GetEntry(iDigit);
     if(pdigit){       
       // first kNch channels are in-time ADC values
       if(iDigit<kNch){
  	  for(int igain=0; igain<=1; igain++) adcInTime[iDigit][igain] = digit.GetADCValue(igain);
       }
       else if(iDigit>=kNch && iDigit<2*kNch){
  	  for(int igain=0; igain<=1; igain++) adcOutOfTime[iDigit-kNch][igain] = digit.GetADCValue(igain);
       }
       else AliWarning(" Looking for wrong index in digit tree: index >kNch !!!\n");
       //
       // Writing cabled signal code for digits!!!!!!!!!!!!!!!!!!!!!!!!
       int det = digit.GetSector(0);
       int quad = digit.GetSector(1);
       if(iDigit<kNch){  // first Nch digits are for in time signals
         signalCode[iDigit] = GetChannelSignal(det, quad, kTRUE);
       } 
       if(testPedSubBit){
	 if(iDigit>=kNch && iDigit<2*kNch){  // 2nd Nch digits are for out of time signals
            signalCode[iDigit] = GetChannelSignal(det, quad, kFALSE);
	  }
      } // if(testPedSubBit)
    }
  } // Loop over digits
  
  // PEDESTAL subtraction
  // Nov 2015: if PedSubMode==kTRUE but coefficients are null, mean value must be subtracted!!!!!
  for(int ich=0; ich<24; ich++){
     if(chPedSubMode[ich]==kTRUE && (TMath::Abs(corrCoeff1[ich])>0. && TMath::Abs(corrCoeff0[ich])>0.)){
       // Pedestal subtraction from correlation ------------------------------------------------
       for(int igain=0; igain<=1; igain++) 
       	   adcCorr[ich][igain] = (float) (adcInTime[ich][igain] - (corrCoeff1[ich+igain*kNch]*adcOutOfTime[ich][igain]+corrCoeff0[ich+igain*kNch]));
     }
     else if((chPedSubMode[ich]==kFALSE) || (chPedSubMode[ich]==kTRUE && ((TMath::Abs(corrCoeff1[ich]))<1e5 && TMath::Abs(corrCoeff0[ich])<1e5))){
        // Pedestal subtraction from mean value ------------------------------------------------
        for(int igain=0; igain<=1; igain++){
	   adcCorr[ich][igain] = (float) (adcInTime[ich][igain] - meanPed[ich+igain*kNch]);
	}
      }
  // Ch. debug
  //printf("ZDC rec.ADC: ch.%d (code %d) rawADC %d subMode %d meanPed %f pedfromcorr = %f  -> corrADC HR %f \n", ich, adcSignal[ich], adcInTime[ich][0], chPedSubMode[ich], meanPed[ich], corrCoeff1[ich]*adcOutOfTime[ich][0]+corrCoeff0[ich], adcCorr[ich][0]);
  }
  
   // Ch. debug
   //printf("\n ---- AliZDCReconstructor: rec from digits\n");
   //for(int ich=0; ich<kNch; ich++) if(adcInTime[ich][0]>0.) printf(" ch.%d signalcode %d rawADC HR %d subMode %d  corrADC HR %f \n", ich, signalCode[ich], adcInTime[ich][0], chPedSubMode[ich], adcCorr[ich][0]);
   
  UInt_t counts[32];
  Int_t  tdc[32][4];
  for(Int_t jj=0; jj<32; jj++){
    counts[jj]=0;
    for(Int_t ii=0; ii<4; ii++) tdc[jj][ii]=0;
  }
  
  Int_t  evQualityBlock[4] = {1,0,0,0};
  Int_t  triggerBlock[4] = {0,0,0,0};
  Int_t  chBlock[3] = {0,0,0};
  UInt_t puBits=0;
  
  // reconstruct the event
  if(fRecoMode==1)
    ReconstructEventpp(clustersTree, adcCorr, signalCode, tdcCabling, kFALSE, counts, tdc,
      evQualityBlock,  triggerBlock,  chBlock, puBits);
  else if(fRecoMode==2)
    ReconstructEventPbPb(clustersTree, adcCorr, signalCode, tdcCabling, kFALSE, counts, tdc,
      evQualityBlock,  triggerBlock,  chBlock, puBits);    
}

//_____________________________________________________________________________
void AliZDCReconstructor::Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const
{
  // *** ZDC raw data reconstruction
  // Works on the current event

  // Retrieving calibration data  
  // Parameters for pedestal subtraction
  int const kNch = 24;
  Float_t meanPed[2*kNch];    
  for(Int_t jj=0; jj<2*kNch; jj++) meanPed[jj] = fPedData->GetMeanPed(jj);
  // Parameters pedestal subtraction through correlation with out-of-time signals
  Float_t corrCoeff0[2*kNch], corrCoeff1[2*kNch];
  for(Int_t jj=0; jj<2*kNch; jj++){
     corrCoeff0[jj] =  fPedData->GetPedCorrCoeff0(jj);
     corrCoeff1[jj] =  fPedData->GetPedCorrCoeff1(jj);
     //printf("  %d   %1.4f  %1.4f\n", jj,corrCoeff0[jj],corrCoeff1[jj]);
  }
  Bool_t testPedSubBit = fPedData->TestPedModeBit();
  //ch. debug
  //printf("   pedSubMode from OCDB object (test bit) %d (FALSE = mean value, TRUE = from corr. w. o.o.t)\n", testPedSubBit);
  
  // Reading mapping from OCDB 
  //fMapping->Print("");
  //
  Int_t tdcSignalCode[32] = {0};
  // **** Adding a fix since kL0 in RUN1 was not set in the mapping (23/7/2015)
  Bool_t iskL0set = kFALSE;
  for(int i=0; i<32; i++){
     tdcSignalCode[i] = fMapping->GetTDCSignalCode(i);
     // Ch. debug
     //printf(" TDC ch.%d  signal %d \n",i, tdcSignalCode[i]);
  
     if(fMapping->GetTDCSignalCode(i) == kL0) iskL0set = kTRUE;
  }
  // if kL0 is not set (RUN1) it is manually set to ch.15
  if(!iskL0set) tdcSignalCode[15] = kL0;
  // Ch. debug
  //printf(" iskL0set %d  tdcSignalCode[15] %d \n",iskL0set, tdcSignalCode[15]);
  
  Bool_t isScalerOn=kFALSE;
  Int_t jsc=0, itdc=0, ihittdc=0;
  Int_t iprevtdc[32];
  UInt_t scalerData[32]={0,};
  Int_t tdcData[32][4];	
  for(Int_t k=0; k<32; k++){
    iprevtdc[k]=-1;
    for(Int_t i=0; i<4; i++) tdcData[k][i]=0;
  }
  
  
  Int_t  evQualityBlock[4] = {1,0,0,0};
  Int_t  triggerBlock[4] = {0,0,0,0};
  Int_t  chBlock[3] = {0,0,0};
  UInt_t puBits=0;

  int adcInTime[kNch][2], adcOutOfTime[kNch][2];
  float adcCorr[kNch][2];
  for(int ich=0; ich<kNch; ich++){
    for(int ig=0; ig<=1; ig++) adcInTime[ich][ig] = adcOutOfTime[ich][ig] = adcCorr[ich][ig] = 0;
  }
  
  Bool_t chPedSubMode[kNch] = {0,};
  if(testPedSubBit){
    for(int i=0; i<kNch; i++) chPedSubMode[i] = fPedData->GetUseCorrFit(i);
  }
  
  int adcSignal[kNch]={0};
  int tdcCabling[7]={-1,-1,-1,-1,-1,-1,-1};

  
  // loop over raw data
  //rawReader->Reset();
  AliZDCRawStream rawData(rawReader);
  while(rawData.Next()){
   
   // ***************************** Reading ADCs
   if((rawData.GetADCModule()>=kFirstADCGeo) && (rawData.GetADCModule()<=kLastADCGeo)){    
    //printf(" **** Reading ADC raw data from module %d **** \n",rawData.GetADCModule());
    //
    if((rawData.IsADCDataWord()) && (rawData.GetNChannelsOn()<48))    chBlock[0] = kTRUE;
    if((rawData.IsADCDataWord()) && (rawData.IsOverflow() == kTRUE))  chBlock[1] = kTRUE;
    if((rawData.IsADCDataWord()) && (rawData.IsUnderflow() == kTRUE)) chBlock[2] = kTRUE;
    if((rawData.IsADCDataWord()) && (rawData.IsADCEventGood() == kTRUE)) evQualityBlock[0] = kTRUE;
    
    if((rawData.IsADCDataWord()) && (rawData.IsUnderflow()==kFALSE) 
        && (rawData.IsOverflow()==kFALSE)){
        //&& (rawData.IsOverflow()==kFALSE) && (rawData.IsADCEventGood()==kTRUE)){
     
      Int_t adcMod = rawData.GetADCModule();
      Int_t adcCh = rawData.GetADCChannel();
      Int_t det = rawData.GetSector(0);
      Int_t quad = rawData.GetSector(1);
      Int_t gain = rawData.GetADCGain();
      Int_t value = rawData.GetADCValue();
      //
      // NB -> index is needed to cpompute which of the 24 signal I want to reconstruct I am dealing with!!!!!!!!!
      // The order is: ZNC(PMC,1,2,3,4) ZPC(PMC,1,2,3,4) ZEM1 ZEM2 ZNA(PMC,1,2,3,4) ZPA(PMC,1,2,3,4) PMREF1 PMREF2
      // that MUST correspond to the oredr with which are written the pedestals!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      int index = -1;
      if(quad!=5){
        if(det==1)	index = quad;
        else if(det==2) index = 5+quad;
        else if(det==3) index = 9+quad; // giusto 9!!!!
        else if(det==4) index = 12+quad;
        else if(det==5) index = 17+quad;
      }
      else index = (det-1)/3+22;     
      //Ch.debug
      //printf("\t AliZDCReconstructor: ADC mod.%d ch.%d det.%d quad.%d, gain.%d   value %d\n",adcMod, adcCh,det,quad,gain,value);
      
      if((adcMod==0 || adcMod==1) && index!=-1){
         adcInTime[index][gain] = value; // in time signals     
         adcSignal[index] = GetChannelSignal(det, quad, kTRUE);
         //Ch.debug
         //printf("   mod.%d ch.%d sigcode %d det.%d  sec.%d  adcInTime[%d][%d] = %d\n\n",adcMod, adcCh, adcSignal[index], det,quad,index, gain, value);
      }
      else if((adcMod==2 || adcMod==3) && index!=-1){
         adcOutOfTime[index][gain] = value; // out of time signals
      }

    }//IsADCDataWord
   }// ADC DATA
   // ***************************** Reading Scaler
   else if(rawData.GetADCModule()==kScalerGeo){
     if(rawData.IsScalerWord()==kTRUE){
       isScalerOn = kTRUE;
       scalerData[jsc] = rawData.GetTriggerCount();
       // Ch. debug
       //printf("   Reconstructed VME Scaler: %d %d  ",jsc,scalerData[jsc]);
       //
       jsc++;
     }
   }// VME SCALER DATA
   // ***************************** Reading ZDC TDC
   else if(rawData.GetADCModule()==kZDCTDCGeo && rawData.IsZDCTDCDatum()==kTRUE){
       itdc = rawData.GetChannel(); 
       // Setting signal code read from OCDB
       rawData.SetCabledSignal(tdcSignalCode[itdc]);
       // setting TDC channels from raw data (=0 if reconstructing from digits!!!!)
       // NB -> I want to store the ch. in tdcCabling variable ONLY if it is not yet stored!!!
       // NB -> THE ORDER MUST BE THE SAME AS THE ONE IN ZDCMAPPINGDA.cxx (7/7/2015)
       //	ZEM1 ZEM2 ZNC ZPC ZNA ZPA
       // otherwise calibrated TDC won't be centered around zero (wrong offset subtracted)
       if(tdcSignalCode[itdc]==kZEM1D && tdcCabling[0]<0) tdcCabling[0] = itdc;
       else if(tdcSignalCode[itdc]==kZEM2D && tdcCabling[1]<0) tdcCabling[1] = itdc;
       else if(tdcSignalCode[itdc]==kZNCD && tdcCabling[2]<0)  tdcCabling[2] = itdc;
       else if(tdcSignalCode[itdc]==kZPCD && tdcCabling[3]<0)  tdcCabling[3] = itdc;
       else if(tdcSignalCode[itdc]==kZNAD && tdcCabling[4]<0)  tdcCabling[4] = itdc;
       else if(tdcSignalCode[itdc]==kZPAD && tdcCabling[5]<0)  tdcCabling[5] = itdc;
       else if(tdcSignalCode[itdc]==kL0 && tdcCabling[6]<0)    tdcCabling[6] = itdc;
       //
       iprevtdc[itdc]++;
       if(iprevtdc[itdc]<4) tdcData[itdc][iprevtdc[itdc]] = rawData.GetZDCTDCDatum();
       // Ch. debug
       //printf("   TDCch.%d hit %d signal %d value %d\n",itdc, iprevtdc[itdc], rawData.GetCabledSignal(),rawData.GetZDCTDCDatum());
   }// ZDC TDC DATA
   // ***************************** Reading PU
   else if(rawData.GetADCModule()==kPUGeo){
     puBits = rawData.GetDetectorPattern();
   }
    // ***************************** Reading trigger history
   else if(rawData.IstriggerHistoryWord()==kTRUE){
     triggerBlock[0] = rawData.IsCPTInputEMDTrigger();
     triggerBlock[1] = rawData.IsCPTInputSemiCentralTrigger();
     triggerBlock[2] = rawData.IsCPTInputCentralTrigger();
     triggerBlock[3] = rawData.IsCPTInputMBTrigger();
   }
  
  }//loop on raw data
  // Ch. debug
  //printf("\n  TDC channels  ZEM1 %d ZEM2 %d ZNC %d ZPC %d ZNA %d  ZPA %d L0 %d\n\n", tdcCabling[0],tdcCabling[1],tdcCabling[2],tdcCabling[3],tdcCabling[4],tdcCabling[5],tdcCabling[6]);
  
  // PEDESTAL subtraction
  // Jul 2015: if PedSubMode==kTRUE but coefficients are null, mean value must be subtracted!!!!!
  for(int ich=0; ich<24; ich++){
     if(chPedSubMode[ich]==kTRUE && (TMath::Abs(corrCoeff1[ich])>0. && TMath::Abs(corrCoeff0[ich])>0.)){
       // Pedestal subtraction from correlation ------------------------------------------------
       for(int igain=0; igain<=1; igain++) 
       	   adcCorr[ich][igain] = (float) (adcInTime[ich][igain] - (corrCoeff1[ich+igain*kNch]*adcOutOfTime[ich][igain]+corrCoeff0[ich+igain*kNch]));
     }
     else if((chPedSubMode[ich]==kFALSE) || (chPedSubMode[ich]==kTRUE && ((TMath::Abs(corrCoeff1[ich]))<1e5 && TMath::Abs(corrCoeff0[ich])<1e5))){
        // Pedestal subtraction from mean value ------------------------------------------------
        for(int igain=0; igain<=1; igain++){
	   adcCorr[ich][igain] = (float) (adcInTime[ich][igain] - meanPed[ich+igain*kNch]);
	}
      }
  // Ch. debug
  //printf("ZDC rec.ADC: ch.%d (code %d) rawADC %d subMode %d meanPed %f pedfromcorr = %f  -> corrADC HR %f \n", ich, adcSignal[ich], adcInTime[ich][0], chPedSubMode[ich], meanPed[ich], corrCoeff1[ich]*adcOutOfTime[ich][0]+corrCoeff0[ich], adcCorr[ich][0]);
  }
  
  
    
  if(fRecoMode==1) // p-p data
    ReconstructEventpp(clustersTree, adcCorr, adcSignal, tdcCabling, isScalerOn, scalerData, tdcData,
      evQualityBlock, triggerBlock, chBlock, puBits);
  else if(fRecoMode==2) // Pb-Pb data
      ReconstructEventPbPb(clustersTree, adcCorr, adcSignal, tdcCabling, isScalerOn, scalerData,  tdcData,
      evQualityBlock, triggerBlock, chBlock, puBits);
}

//_____________________________________________________________________________
void AliZDCReconstructor::ReconstructEventpp(TTree *clustersTree, 
	Float_t adc[24][2], Int_t signalCodeADC[24], Int_t tdcCabling[7], Bool_t isScalerOn, UInt_t* scaler, 
	Int_t tdc[32][4], const Int_t* const evQualityBlock, 
	const Int_t* const triggerBlock, const Int_t* const chBlock, UInt_t puBits) const
{
  // ****************** Reconstruct one event ******************

  int const kNch = 24;
  Float_t corrADCZNC[10] = {0,}, corrADCZPC[10] = {0,};
  Float_t corrADCZNA[10] = {0,}, corrADCZPA[10] = {0,};
  Float_t corrADCZEM1[2] = {0,0}, corrADCZEM2[2] = {0,0};
  Float_t sPMRef1[2] = {0,0}, sPMRef2[2] = {0,0};
  
  for(int i=0; i<kNch; i++){
    if(signalCodeADC[i]==kZNAC){
      for(int igain=0; igain<=1; igain++) corrADCZNA[0+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNA1){
      for(int igain=0; igain<=1; igain++) corrADCZNA[1+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNA2){
      for(int igain=0; igain<=1; igain++) corrADCZNA[2+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNA3){
      for(int igain=0; igain<=1; igain++) corrADCZNA[3+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNA4){
      for(int igain=0; igain<=1; igain++) corrADCZNA[4+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPAC){
      for(int igain=0; igain<=1; igain++) corrADCZPA[0+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPA1){
      for(int igain=0; igain<=1; igain++) corrADCZPA[1+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPA2){
      for(int igain=0; igain<=1; igain++) corrADCZPA[2+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPA3){
      for(int igain=0; igain<=1; igain++) corrADCZPA[3+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPA4){
      for(int igain=0; igain<=1; igain++) corrADCZPA[4+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZEM1){
      for(int igain=0; igain<=1; igain++) corrADCZEM1[igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZEM2){  
      for(int igain=0; igain<=1; igain++) corrADCZEM2[igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNCC){ 
      for(int igain=0; igain<=1; igain++) corrADCZNC[0+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNC1){
      for(int igain=0; igain<=1; igain++) corrADCZNC[1+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNC2){
      for(int igain=0; igain<=1; igain++) corrADCZNC[2+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNC3){
      for(int igain=0; igain<=1; igain++) corrADCZNC[3+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNC4){
      for(int igain=0; igain<=1; igain++) corrADCZNC[4+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPCC){ 
      for(int igain=0; igain<=1; igain++) corrADCZPC[0+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPC1){
      for(int igain=0; igain<=1; igain++) corrADCZPC[1+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPC2){
      for(int igain=0; igain<=1; igain++) corrADCZPC[2+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPC3){
      for(int igain=0; igain<=1; igain++) corrADCZPC[3+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPC4){
      for(int igain=0; igain<=1; igain++) corrADCZPC[4+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZDCCMon){
      for(int igain=0; igain<=1; igain++) sPMRef1[igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZDCAMon){
      for(int igain=0; igain<=1; igain++) sPMRef2[igain] = adc[i][igain];
    }
    
    // Ch. Debug
    //printf("  ADC ch. %d  cabled signal %d\n",i, signalCodeADC[i]);
  }
  
  // CH. debug
  /*printf("\n*************************************************\n");
  printf(" ReconstructEventpp -> values after pedestal subtraction:\n");
  for(int ich=0; ich<24; ich++) printf(" ch.%d ADC hg %1.2f  lg %1.2f\n", ich, adc[ich][0],adc[ich][1]);
  printf("*************************************************\n");
  // CH. debug
  printf("\n*************************************************\n");
  printf(" ADCZNC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
        corrADCZNC[0],corrADCZNC[1],corrADCZNC[2],corrADCZNC[3],corrADCZNC[4]);
  printf(" ADCZPC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
        corrADCZPC[0],corrADCZPC[1],corrADCZPC[2],corrADCZPC[3],corrADCZPC[4]);
  printf(" ADCZNA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
        corrADCZNA[0],corrADCZNA[1],corrADCZNA[2],corrADCZNA[3],corrADCZNA[4]);
  printf(" ADCZPA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
        corrADCZPA[0],corrADCZPA[1],corrADCZPA[2],corrADCZPA[3],corrADCZPA[4]);
  printf(" ADCZEM1 [%1.2f] ADCZEM2 [%1.2f] \n",corrADCZEM1[0],corrADCZEM2[0]);
  printf("*************************************************\n"); */
  
  // ---------------------- Setting reco flags for ESD
  UInt_t rFlags[32];
  for(Int_t ifl=0; ifl<32; ifl++) rFlags[ifl]=0;
  
  if(evQualityBlock[0] == 1) rFlags[31] = 0x0;
  else rFlags[31] = 0x1;
  //
  if(evQualityBlock[1] == 1) rFlags[30] = 0x1;
  if(evQualityBlock[2] == 1) rFlags[29] = 0x1;
  if(evQualityBlock[3] == 1) rFlags[28] = 0x1;

  if(triggerBlock[0] == 1) rFlags[27] = 0x1;
  if(triggerBlock[1] == 1) rFlags[26] = 0x1;
  if(triggerBlock[2] == 1) rFlags[25] = 0x1;
  if(triggerBlock[3] == 1) rFlags[24] = 0x1;
  
  if(chBlock[0] == 1) rFlags[18] = 0x1;
  if(chBlock[1] == 1) rFlags[17] = 0x1;
  if(chBlock[2] == 1) rFlags[16] = 0x1;
    
  rFlags[13] = puBits & 0x00000020;
  rFlags[12] = puBits & 0x00000010;
  rFlags[11] = puBits & 0x00000080;
  rFlags[10] = puBits & 0x00000040;
  rFlags[9]  = puBits & 0x00000020;
  rFlags[8]  = puBits & 0x00000010;
  
  if(corrADCZPC[0]>fSignalThreshold)  rFlags[5] = 0x1;
  if(corrADCZNC[0]>fSignalThreshold)  rFlags[4] = 0x1;
  if(corrADCZEM2[0]>fSignalThreshold) rFlags[3] = 0x1;
  if(corrADCZEM1[0]>fSignalThreshold) rFlags[2] = 0x1;
  if(corrADCZPA[0]>fSignalThreshold)  rFlags[1] = 0x1;
  if(corrADCZNA[0]>fSignalThreshold)  rFlags[0] = 0x1;
  
  UInt_t recoFlag = rFlags[31] << 31 | rFlags[30] << 30 | rFlags[29] << 29 | rFlags[28] << 28 |
             rFlags[27] << 27 | rFlags[26] << 26 | rFlags[25] << 25 | rFlags[24] << 24 |
	     0x0 << 23 | 0x0 << 22 | 0x0 << 21 | 0x0 << 20 |
	     0x0 << 19 | rFlags[18] << 18 |  rFlags[17] << 17 |  rFlags[16] << 16 |
	     0x0 << 15 | 0x0 << 14 | rFlags[13] << 13 | rFlags[12] << 12 | 
             rFlags[11] << 11 |rFlags[10] << 10 | rFlags[9] << 9 | rFlags[8] << 8 |
	     0x0 << 7 | 0x0 << 6 | rFlags[5] << 5 | rFlags[4] << 4 | 
	     rFlags[3] << 3 | rFlags[2] << 2 | rFlags[1] << 1 | rFlags[0];
  // --------------------------------------------------

  // ******	Retrieving calibration data 
  // --- Equalization coefficients ---------------------------------------------
  Float_t equalCoeffZNC[5], equalCoeffZPC[5], equalCoeffZNA[5], equalCoeffZPA[5];
  for(Int_t ji=0; ji<5; ji++){
     equalCoeffZNC[ji] = fTowCalibData->GetZNCEqualCoeff(ji);
     equalCoeffZPC[ji] = fTowCalibData->GetZPCEqualCoeff(ji); 
     equalCoeffZNA[ji] = fTowCalibData->GetZNAEqualCoeff(ji); 
     equalCoeffZPA[ji] = fTowCalibData->GetZPAEqualCoeff(ji); 
  }
  // --- Energy calibration factors ------------------------------------
  Float_t calibEne[6], calibSatZNA[4], calibSatZNC[4];
  // **** Energy calibration coefficient set to 1 
  // **** (no trivial way to calibrate in p-p runs)
  for(Int_t ij=0; ij<6; ij++) calibEne[ij] = fEnCalibData->GetEnCalib(ij);
  //fEnCalibData->Print("");
  for(Int_t ij=0; ij<4; ij++){
    calibSatZNA[ij] = fSatCalibData->GetZNASatCalib(ij);
    calibSatZNC[ij] = fSatCalibData->GetZNCSatCalib(ij);
  }
  
  // ******	Equalization of detector responses ************************
  Float_t equalTowZNC[10], equalTowZNA[10], equalTowZPC[10], equalTowZPA[10];
  for(Int_t gi=0; gi<10; gi++){
     if(gi<5){
       equalTowZNC[gi] = corrADCZNC[gi]*equalCoeffZNC[gi];
       equalTowZPC[gi] = corrADCZPC[gi]*equalCoeffZPC[gi];
       equalTowZNA[gi] = corrADCZNA[gi]*equalCoeffZNA[gi];
       equalTowZPA[gi] = corrADCZPA[gi]*equalCoeffZPA[gi];
     }
     else{
       equalTowZNC[gi] = corrADCZNC[gi]*equalCoeffZNC[gi-5];
       equalTowZPC[gi] = corrADCZPC[gi]*equalCoeffZPC[gi-5];
       equalTowZNA[gi] = corrADCZNA[gi]*equalCoeffZNA[gi-5];
       equalTowZPA[gi] = corrADCZPA[gi]*equalCoeffZPA[gi-5];
     }
  }
  // Ch. debug
  /*printf("\n ------------- EQUALIZATION -------------\n");
  printf(" ADCZNC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZNC[0],equalTowZNC[1],equalTowZNC[2],equalTowZNC[3],equalTowZNC[4]);
  printf(" ADCZPC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZPC[0],equalTowZPC[1],equalTowZPC[2],equalTowZPC[3],equalTowZPC[4]);
  printf(" ADCZNA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZNA[0],equalTowZNA[1],equalTowZNA[2],equalTowZNA[3],equalTowZNA[4]);
  printf(" ADCZPA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZPA[0],equalTowZPA[1],equalTowZPA[2],equalTowZPA[3],equalTowZPA[4]);
  printf(" ----------------------------------------\n");*/
  
  //  *** p-A RUN 2013 -> new calibration object
  //      to take into account saturation in ZN PMC
  //   -> 5th order pol. fun. to be applied BEFORE en. calibration 
  equalTowZNC[0] = equalTowZNC[0] + calibSatZNC[0]*equalTowZNC[0]*equalTowZNC[0] +
  	calibSatZNC[1]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0] +
	calibSatZNC[2]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0] +
	calibSatZNC[3]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0];
  equalTowZNA[0] = equalTowZNA[0] + calibSatZNA[0]*equalTowZNA[0]*equalTowZNA[0] +
  	calibSatZNA[1]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0] +
	calibSatZNA[2]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0] +
	calibSatZNA[3]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0];

  // Ch. debug
  /*printf("\n ------------- SATURATION CORRECTION -------------\n");
  printf(" ZNC PMC %1.2f\n", equalTowZNC[0]);
  printf(" ZNA PMC %1.2f\n", equalTowZNA[0]);
  printf(" ----------------------------------------\n");*/
  
  // ******	Summed response for hadronic calorimeter (SUMMED and then CALIBRATED!)
  Float_t calibSumZNC[2]={0,0}, calibSumZNA[2]={0,0}, calibSumZPC[2]={0,0}, calibSumZPA[2]={0,0};
  for(Int_t gi=0; gi<5; gi++){
       calibSumZNC[0] += equalTowZNC[gi];
       calibSumZPC[0] += equalTowZPC[gi];
       calibSumZNA[0] += equalTowZNA[gi];
       calibSumZPA[0] += equalTowZPA[gi];
       //
       calibSumZNC[1] += equalTowZNC[gi+5];
       calibSumZPC[1] += equalTowZPC[gi+5];
       calibSumZNA[1] += equalTowZNA[gi+5];
       calibSumZPA[1] += equalTowZPA[gi+5];
  }
  // High gain chain
  calibSumZNC[0] = calibSumZNC[0]*calibEne[0];
  calibSumZPC[0] = calibSumZPC[0]*calibEne[1];
  calibSumZNA[0] = calibSumZNA[0]*calibEne[2];
  calibSumZPA[0] = calibSumZPA[0]*calibEne[3];
  // Low gain chain
  calibSumZNC[1] = calibSumZNC[1]*calibEne[0];
  calibSumZPC[1] = calibSumZPC[1]*calibEne[1];
  calibSumZNA[1] = calibSumZNA[1]*calibEne[2];
  calibSumZPA[1] = calibSumZPA[1]*calibEne[3];
  
  // ******	Energy calibration of detector responses
  Float_t calibTowZNC[10], calibTowZNA[10], calibTowZPC[10], calibTowZPA[10];
  for(Int_t gi=0; gi<5; gi++){
     // High gain chain
     calibTowZNC[gi] = equalTowZNC[gi]*calibEne[0];
     calibTowZPC[gi] = equalTowZPC[gi]*calibEne[1];
     calibTowZNA[gi] = equalTowZNA[gi]*calibEne[2];
     calibTowZPA[gi] = equalTowZPA[gi]*calibEne[3];
     // Low gain chain
     calibTowZNC[gi+5] = equalTowZNC[gi+5]*calibEne[0];
     calibTowZPC[gi+5] = equalTowZPC[gi+5]*calibEne[1];
     calibTowZNA[gi+5] = equalTowZNA[gi+5]*calibEne[2];
     calibTowZPA[gi+5] = equalTowZPA[gi+5]*calibEne[3];
  }
  //
  Float_t calibZEM1[2]={0,0}, calibZEM2[2]={0,0};
  calibZEM1[0] = corrADCZEM1[0]*calibEne[4];
  calibZEM1[1] = corrADCZEM1[1]*calibEne[4];
  calibZEM2[0] = corrADCZEM2[0]*calibEne[5];
  calibZEM2[1] = corrADCZEM2[1]*calibEne[5];
  // Ch. debug
  /*printf("\n ------------- CALIBRATION -------------\n");
  printf(" ADCZNC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZNC[0],calibTowZNC[1],calibTowZNC[2],calibTowZNC[3],calibTowZNC[4]);
  printf(" ADCZPC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZPC[0],calibTowZPC[1],calibTowZPC[2],calibTowZPC[3],calibTowZPC[4]);
  printf(" ADCZNA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZNA[0],calibTowZNA[1],calibTowZNA[2],calibTowZNA[3],calibTowZNA[4]);
  printf(" ADCZPA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZPA[0],calibTowZPA[1],calibTowZPA[2],calibTowZPA[3],calibTowZPA[4]);
  printf(" ADCZEM1 [%1.2f] ADCZEM2 [%1.2f] \n",calibZEM1[0],calibZEM2[0]);
  printf(" ----------------------------------------\n");*/
  
  //  ******	No. of spectator and participants nucleons
  //  Variables calculated to comply with ESD structure
  //  *** N.B. -> They have a meaning only in Pb-Pb!!!!!!!!!!!!
  Int_t nDetSpecNLeft=0, nDetSpecPLeft=0, nDetSpecNRight=0, nDetSpecPRight=0;
  Int_t nGenSpec=0, nGenSpecLeft=0, nGenSpecRight=0;
  Int_t nPart=0, nPartTotLeft=0, nPartTotRight=0;
  Double_t impPar=0., impPar1=0., impPar2=0.;
  
  Bool_t energyFlag = kFALSE;
  // create the output tree
  AliZDCReco* reco = new AliZDCReco(calibSumZNC, calibSumZPC, calibSumZNA, calibSumZPA, 
  		   	calibTowZNC, calibTowZPC, calibTowZNA, calibTowZPA, 
		   	calibZEM1, calibZEM2, sPMRef1, sPMRef2,
		   	nDetSpecNLeft, nDetSpecPLeft, nDetSpecNRight, nDetSpecPRight, 
		   	nGenSpec, nGenSpecLeft, nGenSpecRight, 
		   	nPart, nPartTotLeft, nPartTotRight, 
		   	impPar, impPar1, impPar2,
		   	recoFlag, energyFlag, isScalerOn, scaler, tdc, tdcCabling);
		  
  const Int_t kBufferSize = 4000;
  clustersTree->Branch("ZDC", "AliZDCReco", &reco, kBufferSize);
  // write the output tree
  clustersTree->Fill();
  delete reco;
}

//_____________________________________________________________________________
void AliZDCReconstructor::ReconstructEventPbPb(TTree *clustersTree, 
	Float_t adc[24][2], Int_t signalCodeADC[24], Int_t tdcCabling[7], Bool_t isScalerOn, UInt_t* scaler, 
	Int_t tdc[32][4], const Int_t* const evQualityBlock, 
	const Int_t* const triggerBlock, const Int_t* const chBlock, UInt_t puBits) const
{
  // ****************** Reconstruct one event ******************

  int const kNch = 24;
  Float_t corrADCZNC[10] = {0,}, corrADCZPC[10] = {0,};
  Float_t corrADCZNA[10] = {0,}, corrADCZPA[10] = {0,};
  Float_t corrADCZEM1[2] = {0,0}, corrADCZEM2[2] = {0,0};
  Float_t sPMRef1[2] = {0,0}, sPMRef2[2] = {0,0};

  for(int i=0; i<kNch; i++){
     if(signalCodeADC[i]==kZNAC){
      for(int igain=0; igain<=1; igain++) corrADCZNA[0+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNA1){
      for(int igain=0; igain<=1; igain++) corrADCZNA[1+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNA2){
      for(int igain=0; igain<=1; igain++) corrADCZNA[2+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNA3){
      for(int igain=0; igain<=1; igain++) corrADCZNA[3+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNA4){
      for(int igain=0; igain<=1; igain++) corrADCZNA[4+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPAC){
      for(int igain=0; igain<=1; igain++) corrADCZPA[0+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPA1){
      for(int igain=0; igain<=1; igain++) corrADCZPA[1+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPA2){
      for(int igain=0; igain<=1; igain++) corrADCZPA[2+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPA3){
      for(int igain=0; igain<=1; igain++) corrADCZPA[3+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPA4){
      for(int igain=0; igain<=1; igain++) corrADCZPA[4+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZEM1){
      for(int igain=0; igain<=1; igain++) corrADCZEM1[igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZEM2){  
      for(int igain=0; igain<=1; igain++) corrADCZEM2[igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNCC){ 
      for(int igain=0; igain<=1; igain++) corrADCZNC[0+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNC1){
      for(int igain=0; igain<=1; igain++) corrADCZNC[1+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNC2){
      for(int igain=0; igain<=1; igain++) corrADCZNC[2+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNC3){
      for(int igain=0; igain<=1; igain++) corrADCZNC[3+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZNC4){
      for(int igain=0; igain<=1; igain++) corrADCZNC[4+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPCC){ 
      for(int igain=0; igain<=1; igain++) corrADCZPC[0+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPC1){
      for(int igain=0; igain<=1; igain++) corrADCZPC[1+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPC2){
      for(int igain=0; igain<=1; igain++) corrADCZPC[2+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPC3){
      for(int igain=0; igain<=1; igain++) corrADCZPC[3+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZPC4){
      for(int igain=0; igain<=1; igain++) corrADCZPC[4+5*igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZDCCMon){
      for(int igain=0; igain<=1; igain++) sPMRef1[igain] = adc[i][igain];
    }
    else if(signalCodeADC[i]==kZDCAMon){
      for(int igain=0; igain<=1; igain++) sPMRef2[igain] = adc[i][igain];
    }
        
    // Ch. Debug
    //printf("  ADC ch. %d  cabled signal %d\n",i, signalCodeADC[i]);
  }
  
  // CH. debug
  /*printf("\n*************************************************\n");
  printf(" ReconstructEventPbPb -> values after pedestal subtraction:\n");
  for(int ich=0; ich<24; ich++) printf(" ch.%d ADC hg %1.2f  lg %1.2f\n", ich, adc[ich][0],adc[ich][1]);
  printf("*************************************************\n");*/
  // CH. debug
  /*printf("\n*************************************************\n");
  printf(" ADCZNC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
        corrADCZNC[0],corrADCZNC[1],corrADCZNC[2],corrADCZNC[3],corrADCZNC[4]);
  printf(" ADCZPC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
        corrADCZPC[0],corrADCZPC[1],corrADCZPC[2],corrADCZPC[3],corrADCZPC[4]);
  printf(" ADCZPC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
        corrADCZPC[0],corrADCZPC[1],corrADCZPC[2],corrADCZPC[3],corrADCZPC[4]);
  printf(" ADCZPA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
        corrADCZPA[0],corrADCZPA[1],corrADCZPA[2],corrADCZPA[3],corrADCZPA[4]);
  printf(" ADCZEM1 [%1.2f] ADCZEM2 [%1.2f] \n",corrADCZEM1[0],corrADCZEM2[0]);
  printf("*************************************************\n"); */

  // ---------------------- Setting reco flags for ESD
  UInt_t rFlags[32];
  for(Int_t ifl=0; ifl<32; ifl++) rFlags[ifl]=0;
  
  if(evQualityBlock[0] == 1) rFlags[31] = 0x0;
  else rFlags[31] = 0x1;
  //
  if(evQualityBlock[1] == 1) rFlags[30] = 0x1;
  if(evQualityBlock[2] == 1) rFlags[29] = 0x1;
  if(evQualityBlock[3] == 1) rFlags[28] = 0x1;

  if(triggerBlock[0] == 1) rFlags[27] = 0x1;
  if(triggerBlock[1] == 1) rFlags[26] = 0x1;
  if(triggerBlock[2] == 1) rFlags[25] = 0x1;
  if(triggerBlock[3] == 1) rFlags[24] = 0x1;
  
  if(chBlock[0] == 1) rFlags[18] = 0x1;
  if(chBlock[1] == 1) rFlags[17] = 0x1;
  if(chBlock[2] == 1) rFlags[16] = 0x1;
  
  rFlags[13] = puBits & 0x00000020;
  rFlags[12] = puBits & 0x00000010;
  rFlags[11] = puBits & 0x00000080;
  rFlags[10] = puBits & 0x00000040;
  rFlags[9]  = puBits & 0x00000020;
  rFlags[8]  = puBits & 0x00000010;  
  
  if(corrADCZPC[0]>fSignalThreshold)  rFlags[5] = 0x1;
  if(corrADCZNC[0]>fSignalThreshold)  rFlags[4] = 0x1;
  if(corrADCZEM2[0]>fSignalThreshold) rFlags[3] = 0x1;
  if(corrADCZEM1[0]>fSignalThreshold) rFlags[2] = 0x1;
  if(corrADCZPA[0]>fSignalThreshold)  rFlags[1] = 0x1;
  if(corrADCZNA[0]>fSignalThreshold)  rFlags[0] = 0x1;

  UInt_t recoFlag = rFlags[31] << 31 | rFlags[30] << 30 | rFlags[29] << 29 | rFlags[28] << 28 |
             rFlags[27] << 27 | rFlags[26] << 26 | rFlags[25] << 25 | rFlags[24] << 24 |
	     0x0 << 23 | 0x0 << 22 | 0x0 << 21 | 0x0 << 20 |
	     0x0 << 19 | rFlags[18] << 18 |  rFlags[17] << 17 |  rFlags[16] << 16 |
	     0x0 << 15 | 0x0 << 14 | rFlags[13] << 13 | rFlags[12] << 12 | 
             rFlags[11] << 11 |rFlags[10] << 10 | rFlags[9] << 9 | rFlags[8] << 8 |
	     0x0 << 7 | 0x0 << 6 | rFlags[5] << 5 | rFlags[4] << 4 | 
	     rFlags[3] << 3 | rFlags[2] << 2 | rFlags[1] << 1 | rFlags[0];
  // --------------------------------------------------
  
  
  // ******	Retrieving calibration data 
  // --- Equalization coefficients ---------------------------------------------
  Float_t equalCoeffZNC[5], equalCoeffZPC[5], equalCoeffZNA[5], equalCoeffZPA[5];
  for(Int_t ji=0; ji<5; ji++){
     equalCoeffZNC[ji] = fTowCalibData->GetZNCEqualCoeff(ji);
     equalCoeffZPC[ji] = fTowCalibData->GetZPCEqualCoeff(ji); 
     equalCoeffZNA[ji] = fTowCalibData->GetZNAEqualCoeff(ji); 
     equalCoeffZPA[ji] = fTowCalibData->GetZPAEqualCoeff(ji); 
  }
  // --- Energy calibration factors ------------------------------------
  Float_t calibEne[6], calibSatZNA[4], calibSatZNC[4];
  // **** Energy calibration coefficient set to 1 
  // **** (no trivial way to calibrate in p-p runs)
  for(Int_t ij=0; ij<6; ij++) calibEne[ij] = fEnCalibData->GetEnCalib(ij);
  for(Int_t ij=0; ij<4; ij++){
    calibSatZNA[ij] = fSatCalibData->GetZNASatCalib(ij);
    calibSatZNC[ij] = fSatCalibData->GetZNCSatCalib(ij);
  }
  
  // ******	Equalization of detector responses
  Float_t equalTowZNC[10], equalTowZNA[10], equalTowZPC[10], equalTowZPA[10];
  for(Int_t gi=0; gi<10; gi++){
     if(gi<5){
       equalTowZNC[gi] = corrADCZNC[gi]*equalCoeffZNC[gi];
       equalTowZPC[gi] = corrADCZPC[gi]*equalCoeffZPC[gi];
       equalTowZNA[gi] = corrADCZNA[gi]*equalCoeffZNA[gi];
       equalTowZPA[gi] = corrADCZPA[gi]*equalCoeffZPA[gi];
     }
     else{
       equalTowZNC[gi] = corrADCZNC[gi]*equalCoeffZNC[gi-5];
       equalTowZPC[gi] = corrADCZPC[gi]*equalCoeffZPC[gi-5];
       equalTowZNA[gi] = corrADCZNA[gi]*equalCoeffZNA[gi-5];
       equalTowZPA[gi] = corrADCZPA[gi]*equalCoeffZPA[gi-5];
     }
  }
  
  // Ch. debug
  /*printf("\n ------------- EQUALIZATION -------------\n");
  printf(" ADCZNC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZNC[0],equalTowZNC[1],equalTowZNC[2],equalTowZNC[3],equalTowZNC[4]);
  printf(" ADCZPC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZPC[0],equalTowZPC[1],equalTowZPC[2],equalTowZPC[3],equalTowZPC[4]);
  printf(" ADCZNA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZNA[0],equalTowZNA[1],equalTowZNA[2],equalTowZNA[3],equalTowZNA[4]);
  printf(" ADCZPA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZPA[0],equalTowZPA[1],equalTowZPA[2],equalTowZPA[3],equalTowZPA[4]);
  printf(" ----------------------------------------\n");*/
  
  //  *** p-A RUN 2013 -> new calibration object
  //      to take into account saturation in ZN PMC
  //   -> 5th order pol. fun. to be applied BEFORE en. calibration 
  equalTowZNC[0] = equalTowZNC[0] + calibSatZNC[0]*equalTowZNC[0]*equalTowZNC[0] +
  	calibSatZNC[1]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0] +
	calibSatZNC[2]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0] +
	calibSatZNC[3]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0]*equalTowZNC[0];
  equalTowZNA[0] = equalTowZNA[0] + calibSatZNA[0]*equalTowZNA[0]*equalTowZNA[0] +
  	calibSatZNA[1]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0] +
	calibSatZNA[2]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0] +
	calibSatZNA[3]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0]*equalTowZNA[0];
  
  // ******	Summed response for hadronic calorimeter (SUMMED and then CALIBRATED!)
  Float_t calibSumZNC[]={0,0}, calibSumZNA[]={0,0}, calibSumZPC[]={0,0}, calibSumZPA[]={0,0};
  for(Int_t gi=0; gi<5; gi++){
       calibSumZNC[0] += equalTowZNC[gi];
       calibSumZPC[0] += equalTowZPC[gi];
       calibSumZNA[0] += equalTowZNA[gi];
       calibSumZPA[0] += equalTowZPA[gi];
       //
       calibSumZNC[1] += equalTowZNC[gi+5];
       calibSumZPC[1] += equalTowZPC[gi+5];
       calibSumZNA[1] += equalTowZNA[gi+5];
       calibSumZPA[1] += equalTowZPA[gi+5];
  }
   
  // High gain chain
  // NB -> The calibration factor is extracted from *8 chain
  //       i.e. it is equal to 1.38 TeV/100 ch. ADC -> must be *8. for HG chain!!!!
  calibSumZNC[0] = calibSumZNC[0]*calibEne[0]*8.;
  calibSumZPC[0] = calibSumZPC[0]*calibEne[1]*8.;
  calibSumZNA[0] = calibSumZNA[0]*calibEne[2]*8.;
  calibSumZPA[0] = calibSumZPA[0]*calibEne[3]*8.;
  // Low gain chain
  calibSumZNC[1] = calibSumZNC[1]*calibEne[0];
  calibSumZPC[1] = calibSumZPC[1]*calibEne[1];
  calibSumZNA[1] = calibSumZNA[1]*calibEne[2];
  calibSumZPA[1] = calibSumZPA[1]*calibEne[3];
  //
  Float_t calibZEM1[2]={0,0}, calibZEM2[2]={0,0};
  calibZEM1[0] = corrADCZEM1[0]*calibEne[4];
  calibZEM1[1] = corrADCZEM1[1]*calibEne[4];
  calibZEM2[0] = corrADCZEM2[0]*calibEne[5];
  calibZEM2[1] = corrADCZEM2[1]*calibEne[5];
    
  // ******	Energy calibration of detector responses
  Float_t calibTowZNC[10], calibTowZNA[10], calibTowZPC[10], calibTowZPA[10];
  for(Int_t gi=0; gi<5; gi++){
     // High gain chain
     calibTowZNC[gi] = equalTowZNC[gi]*2*calibEne[0]*8.;
     calibTowZPC[gi] = equalTowZPC[gi]*2*calibEne[1]*8.;
     calibTowZNA[gi] = equalTowZNA[gi]*2*calibEne[2]*8.;
     calibTowZPA[gi] = equalTowZPA[gi]*2*calibEne[3]*8.;
     // Low gain chain
     calibTowZNC[gi+5] = equalTowZNC[gi+5]*2*calibEne[0];
     calibTowZPC[gi+5] = equalTowZPC[gi+5]*2*calibEne[1];
     calibTowZNA[gi+5] = equalTowZNA[gi+5]*2*calibEne[2];
     calibTowZPA[gi+5] = equalTowZPA[gi+5]*2*calibEne[3];
  }

  // Ch. debug
  /*printf("\n ------------- CALIBRATION -------------\n");
  printf(" ADCZNC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZNC[0],calibTowZNC[1],calibTowZNC[2],calibTowZNC[3],calibTowZNC[4]);
  printf(" ADCZPC [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZPC[0],calibTowZPC[1],calibTowZPC[2],calibTowZPC[3],calibTowZPC[4]);
  printf(" ADCZNA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZNA[0],calibTowZNA[1],calibTowZNA[2],calibTowZNA[3],calibTowZNA[4]);
  printf(" ADCZPA [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZPA[0],calibTowZPA[1],calibTowZPA[2],calibTowZPA[3],calibTowZPA[4]);
  printf(" ADCZEM1 [%1.2f] ADCZEM2 [%1.2f] \n",calibZEM1[0],calibZEM2[0]);
  printf(" ----------------------------------------\n");*/
  
  //  ******	Number of detected spectator nucleons
  Int_t nDetSpecNLeft=0, nDetSpecPLeft=0, nDetSpecNRight=0, nDetSpecPRight=0;
  if(fBeamEnergy>0.01){
    nDetSpecNLeft = (Int_t) (calibSumZNC[0]/fBeamEnergy);
    nDetSpecPLeft = (Int_t) (calibSumZPC[0]/fBeamEnergy);
    nDetSpecNRight = (Int_t) (calibSumZNA[0]/fBeamEnergy);
    nDetSpecPRight = (Int_t) (calibSumZPA[0]/fBeamEnergy);
  }
  else AliWarning(" ATTENTION!!! fBeamEnergy=0 -> N_spec will be ZERO!!! \n");
  
  Int_t nGenSpec=0, nGenSpecA=0, nGenSpecC=0;
  Int_t nPart=0, nPartA=0, nPartC=0;
  Double_t b=0., bA=0., bC=0.;
  
  /*if(fIsCalibrationMB == kFALSE){
   // ******   Reconstruction parameters ------------------ 
   if(!fgRecoParam) fgRecoParam = const_cast<AliZDCRecoParam*>(GetRecoParam());
   if(!fgRecoParam){  
     AliError("  RecoParam object not retrieved correctly: not reconstructing ZDC event!!!");
     return;
   }
   TH1D* hNpartDist = fgRecoParam->GethNpartDist();
   TH1D*   hbDist = fgRecoParam->GethbDist();	 
   Float_t  fClkCenter = fgRecoParam->GetClkCenter();
   if(!hNpartDist || !hbDist){  
      AliError("Something wrong in Glauber MC histos got from AliZDCREcoParamPbPb: NO EVENT RECO FOR ZDC DATA!!!\n\n");
      //return;
   }
   else{  
    if(!fgMBCalibData) fgMBCalibData = const_cast<AliZDCMBCalib*>(GetMBCalibData()); 
    TH2F *hZDCvsZEM  = fgMBCalibData->GethZDCvsZEM();
    TH2F *hZDCCvsZEM = fgMBCalibData->GethZDCCvsZEM();
    TH2F *hZDCAvsZEM = fgMBCalibData->GethZDCAvsZEM();
    //
    Double_t xHighEdge = hZDCvsZEM->GetXaxis()->GetXmax();
    Double_t origin = xHighEdge*fClkCenter;
    // Ch. debug
    //printf("\n\n  xHighEdge %1.2f, origin %1.4f \n", xHighEdge, origin);
    //
    // ====> Summed ZDC info (sideA+side C)
    TF1 *line = new TF1("line","[0]*x+[1]",0.,xHighEdge);
    Float_t y = (calibSumZNC[0]+calibSumZPC[0]+calibSumZNA[0]+calibSumZPA[0])/1000.;
    Float_t x = (calibZEM1[0]+calibZEM2[0])/1000.;
    line->SetParameter(0, y/(x-origin));
    line->SetParameter(1, -origin*y/(x-origin));
    // Ch. debug
    //printf("  ***************** Summed ZDC info (sideA+side C) \n");
    //printf("  E_{ZEM} %1.4f, E_{ZDC} %1.2f, TF1: %1.2f*x + %1.2f   ", x, y,y/(x-origin),-origin*y/(x-origin));
    //
    Double_t countPerc=0;
    Double_t xBinCenter=0, yBinCenter=0;
    for(Int_t nbinx=1; nbinx<=hZDCvsZEM->GetNbinsX(); nbinx++){
      for(Int_t nbiny=1; nbiny<=hZDCvsZEM->GetNbinsY(); nbiny++){
         xBinCenter = hZDCvsZEM->GetXaxis()->GetBinCenter(nbinx);
         yBinCenter = hZDCvsZEM->GetYaxis()->GetBinCenter(nbiny);
         //
	 if(line->GetParameter(0)>0){
           if(yBinCenter < (line->GetParameter(0)*xBinCenter + line->GetParameter(1))){
             countPerc += hZDCvsZEM->GetBinContent(nbinx,nbiny);
             // Ch. debug
             //printf(" xBinCenter  %1.3f, yBinCenter %1.0f,  countPerc %1.0f\n", 
	 	//xBinCenter, yBinCenter, countPerc);
           }
   	 }
   	 else{
   	   if(yBinCenter > (line->GetParameter(0)*xBinCenter + line->GetParameter(1))){
   	     countPerc += hZDCvsZEM->GetBinContent(nbinx,nbiny);
             // Ch. debug
             //printf(" xBinCenter  %1.3f, yBinCenter %1.0f,  countPerc %1.0f\n", 
	 	//xBinCenter, yBinCenter, countPerc);
   	   }
   	 }
      }
    }
    //
    Double_t xSecPerc = 0.;
    if(hZDCvsZEM->GetEntries()!=0){ 
      xSecPerc = countPerc/hZDCvsZEM->GetEntries();
    }
    else{
      AliWarning("  Histogram hZDCvsZEM from OCDB has no entries!!!");
    }
    // Ch. debug
    //printf("  xSecPerc %1.4f  \n", xSecPerc);

    // ====> side C
    TF1 *lineC = new TF1("lineC","[0]*x+[1]",0.,xHighEdge);
    Float_t yC = (calibSumZNC[0]+calibSumZPC[0])/1000.;
    lineC->SetParameter(0, yC/(x-origin));
    lineC->SetParameter(1, -origin*yC/(x-origin));
    // Ch. debug
    //printf("  ***************** Side C \n");
    //printf("  E_{ZEM} %1.4f, E_{ZDCC} %1.2f, TF1: %1.2f*x + %1.2f   ", x, yC,yC/(x-origin),-origin*yC/(x-origin));
    //
    Double_t countPercC=0;
    Double_t xBinCenterC=0, yBinCenterC=0;
    for(Int_t nbinx=1; nbinx<=hZDCCvsZEM->GetNbinsX(); nbinx++){
      for(Int_t nbiny=1; nbiny<=hZDCCvsZEM->GetNbinsY(); nbiny++){
    	 xBinCenterC = hZDCCvsZEM->GetXaxis()->GetBinCenter(nbinx);
    	 yBinCenterC = hZDCCvsZEM->GetYaxis()->GetBinCenter(nbiny);
    	 if(lineC->GetParameter(0)>0){
           if(yBinCenterC < (lineC->GetParameter(0)*xBinCenterC + lineC->GetParameter(1))){
             countPercC += hZDCCvsZEM->GetBinContent(nbinx,nbiny);
           }
    	 }
    	 else{
    	   if(yBinCenterC > (lineC->GetParameter(0)*xBinCenterC + lineC->GetParameter(1))){
    	     countPercC += hZDCCvsZEM->GetBinContent(nbinx,nbiny);
    	   }
    	 }
      }
    }
    //
    Double_t xSecPercC = 0.;
    if(hZDCCvsZEM->GetEntries()!=0){ 
      xSecPercC = countPercC/hZDCCvsZEM->GetEntries();
    }
    else{
      AliWarning("  Histogram hZDCCvsZEM from OCDB has no entries!!!");
    }
    // Ch. debug
    //printf("  xSecPercC %1.4f  \n", xSecPercC);
    
    // ====> side A
    TF1 *lineA = new TF1("lineA","[0]*x+[1]",0.,xHighEdge);
    Float_t yA = (calibSumZNA[0]+calibSumZPA[0])/1000.;
    lineA->SetParameter(0, yA/(x-origin));
    lineA->SetParameter(1, -origin*yA/(x-origin));
    //
    // Ch. debug
    //printf("  ***************** Side A \n");
    //printf("  E_{ZEM} %1.4f, E_{ZDCA} %1.2f, TF1: %1.2f*x + %1.2f   ", x, yA,yA/(x-origin),-origin*yA/(x-origin));
    //
    Double_t countPercA=0;
    Double_t xBinCenterA=0, yBinCenterA=0;
    for(Int_t nbinx=1; nbinx<=hZDCAvsZEM->GetNbinsX(); nbinx++){
      for(Int_t nbiny=1; nbiny<=hZDCAvsZEM->GetNbinsY(); nbiny++){
    	 xBinCenterA = hZDCAvsZEM->GetXaxis()->GetBinCenter(nbinx);
    	 yBinCenterA = hZDCAvsZEM->GetYaxis()->GetBinCenter(nbiny);
    	 if(lineA->GetParameter(0)>0){
           if(yBinCenterA < (lineA->GetParameter(0)*xBinCenterA + lineA->GetParameter(1))){
             countPercA += hZDCAvsZEM->GetBinContent(nbinx,nbiny);
           }
   	 }
   	 else{
   	   if(yBinCenterA > (lineA->GetParameter(0)*xBinCenterA + lineA->GetParameter(1))){
   	     countPercA += hZDCAvsZEM->GetBinContent(nbinx,nbiny);
   	   }
   	 }
      }
    }
    //
    Double_t xSecPercA = 0.;
    if(hZDCAvsZEM->GetEntries()!=0){ 
      xSecPercA = countPercA/hZDCAvsZEM->GetEntries();
    }
    else{
      AliWarning("  Histogram hZDCAvsZEM from OCDB has no entries!!!");
    }
    // Ch. debug
    //printf("  xSecPercA %1.4f  \n", xSecPercA);
    
    //  ******    Number of participants (from E_ZDC vs. E_ZEM correlation)
    Double_t nPartFrac=0., nPartFracC=0., nPartFracA=0.;
    for(Int_t npbin=1; npbin<hNpartDist->GetNbinsX(); npbin++){
      nPartFrac += (hNpartDist->GetBinContent(npbin))/(hNpartDist->GetEntries());
      if((1.-nPartFrac) < xSecPerc){
    	nPart = (Int_t) hNpartDist->GetBinLowEdge(npbin);
        // Ch. debug
        //printf("  ***************** Summed ZDC info (sideA+side C) \n");
        //printf("  nPartFrac %1.4f, nPart %d\n", nPartFrac, nPart);
    	break;
      }
    }
    if(nPart<0) nPart=0;
    //
    for(Int_t npbin=1; npbin<hNpartDist->GetNbinsX(); npbin++){
      nPartFracC += (hNpartDist->GetBinContent(npbin))/(hNpartDist->GetEntries());
      if((1.-nPartFracC) < xSecPercC){
    	nPartC = (Int_t) hNpartDist->GetBinLowEdge(npbin);
        // Ch. debug
        //printf("  ***************** Side C \n");
        //printf("  nPartFracC %1.4f, nPartC %d\n", nPartFracC, nPartC);
    	break;
    }
    }
    if(nPartC<0) nPartC=0;
    //
    for(Int_t npbin=1; npbin<hNpartDist->GetNbinsX(); npbin++){
      nPartFracA += (hNpartDist->GetBinContent(npbin))/(hNpartDist->GetEntries());
      if((1.-nPartFracA) < xSecPercA){
    	nPartA = (Int_t) hNpartDist->GetBinLowEdge(npbin);
        // Ch. debug
        //printf("  ***************** Side A \n");
        //printf("  nPartFracA %1.4f, nPartA %d\n\n", nPartFracA, nPartA);
        break;
      }
    }
    if(nPartA<0) nPartA=0;
    
    //  ******    Impact parameter (from E_ZDC vs. E_ZEM correlation)
    Double_t bFrac=0., bFracC=0., bFracA=0.;
    for(Int_t ibbin=1; ibbin<hbDist->GetNbinsX(); ibbin++){
      bFrac += (hbDist->GetBinContent(ibbin))/(hbDist->GetEntries());
      if(bFrac > xSecPerc){
   	b = hbDist->GetBinLowEdge(ibbin);
   	break;
      }
    }
    //
    for(Int_t ibbin=1; ibbin<hbDist->GetNbinsX(); ibbin++){
      bFracC += (hbDist->GetBinContent(ibbin))/(hbDist->GetEntries());
      if(bFracC > xSecPercC){
   	bC = hbDist->GetBinLowEdge(ibbin);
   	break;
      }
    }
    //
    for(Int_t ibbin=1; ibbin<hbDist->GetNbinsX(); ibbin++){
      bFracA += (hbDist->GetBinContent(ibbin))/(hbDist->GetEntries());
      if(bFracA > xSecPercA){
   	bA = hbDist->GetBinLowEdge(ibbin);
   	break;
      }
    }

    //  ******	Number of spectator nucleons 
    nGenSpec = 416 - nPart;
    nGenSpecC = 416 - nPartC;
    nGenSpecA = 416 - nPartA;
    if(nGenSpec>416) nGenSpec=416; if(nGenSpec<0) nGenSpec=0;
    if(nGenSpecC>416) nGenSpecC=416; if(nGenSpecC<0) nGenSpecC=0;
    if(nGenSpecA>416) nGenSpecA=416; if(nGenSpecA<0) nGenSpecA=0;    
    
    delete line; 
    delete lineC;  delete lineA;
   }
  }*/ // ONLY IF fIsCalibrationMB==kFALSE
  
  Bool_t energyFlag = kTRUE;  
  AliZDCReco* reco = new AliZDCReco(calibSumZNC, calibSumZPC, calibSumZNA, calibSumZPA, 
  	  	  calibTowZNC, calibTowZPC, calibTowZNA, calibTowZPA, 
		  calibZEM1, calibZEM2, sPMRef1, sPMRef2,
		  nDetSpecNLeft, nDetSpecPLeft, nDetSpecNRight, nDetSpecPRight, 
		  nGenSpec, nGenSpecA, nGenSpecC, 
		  nPart, nPartA, nPartC, b, bA, bC,
		  recoFlag, energyFlag, isScalerOn, scaler, tdc, tdcCabling);
		    
  const Int_t kBufferSize = 4000;
  clustersTree->Branch("ZDC", "AliZDCReco", &reco, kBufferSize);
  //reco->Print("");
  // write the output tree
  clustersTree->Fill();
  delete reco;
}


//_____________________________________________________________________________
void AliZDCReconstructor::FillZDCintoESD(TTree *clustersTree, AliESDEvent* esd) const
{
  // fill energies and number of participants to the ESD

  // Retrieving TDC calibration data  
  // Parameters for TDC centering around zero
  int const knTDC = 6;
  Float_t tdcOffset[knTDC];
  for(Int_t jj=0; jj<knTDC; jj++) tdcOffset[jj] = fTDCCalibData->GetMeanTDC(jj);
  //fTDCCalibData->Print("");

  AliZDCReco reco;
  AliZDCReco* preco = &reco;
  clustersTree->SetBranchAddress("ZDC", &preco);
  clustersTree->GetEntry(0);

  fESDZDC->SetZDC(reco.GetZN1HREnergy(), reco.GetZP1HREnergy(), 
  		reco.GetZEM1HRsignal(), reco.GetZEM2HRsignal(), 
		reco.GetZN2HREnergy(), reco.GetZP2HREnergy(), 
		reco.GetNParticipants(), reco.GetNPartSideA(), reco.GetNPartSideC(),
		reco.GetImpParameter(), reco.GetImpParSideA(), reco.GetImpParSideC(),
		reco.GetRecoFlag());

  Float_t tZNCEne[5], tZNAEne[5], tZPCEne[5], tZPAEne[5];
  Float_t tZNCEneLR[5], tZNAEneLR[5], tZPCEneLR[5], tZPAEneLR[5];
  for(Int_t i=0; i<5; i++){
     tZNCEne[i] = reco.GetZN1HREnTow(i);
     tZNAEne[i] = reco.GetZN2HREnTow(i);
     tZPCEne[i] = reco.GetZP1HREnTow(i);
     tZPAEne[i] = reco.GetZP2HREnTow(i);
     //
     tZNCEneLR[i] = reco.GetZN1LREnTow(i);
     tZNAEneLR[i] = reco.GetZN2LREnTow(i);
     tZPCEneLR[i] = reco.GetZP1LREnTow(i);
     tZPAEneLR[i] = reco.GetZP2LREnTow(i);
  }
  //
  fESDZDC->SetZN1TowerEnergy(tZNCEne);
  fESDZDC->SetZN2TowerEnergy(tZNAEne);
  fESDZDC->SetZP1TowerEnergy(tZPCEne);
  fESDZDC->SetZP2TowerEnergy(tZPAEne);
  //
  fESDZDC->SetZN1TowerEnergyLR(tZNCEneLR);
  fESDZDC->SetZN2TowerEnergyLR(tZNAEneLR);
  fESDZDC->SetZP1TowerEnergyLR(tZPCEneLR);
  fESDZDC->SetZP2TowerEnergyLR(tZPAEneLR);
  // 
  // Writing ZDC scaler for cross section calculation
  // ONLY IF the scaler has been read during the event
  if(reco.IsScalerOn()==kTRUE){
    UInt_t counts[32];
    for(Int_t jk=0; jk<32; jk++) counts[jk] = reco.GetZDCScaler(jk);
    fESDZDC->SetZDCScaler(counts);
  }    
  
  Int_t *tdcCabling = reco.GetTDCchCabling();
  //
  // Ch. debug
  //printf("\n  FillZDCintoESD: TDC channels  ZEM1 %d  ZEM2 %d ZNC %d ZPC %d ZNA %d ZPA %d L0 %d\n\n",  tdcCabling[0],tdcCabling[1],tdcCabling[2],tdcCabling[3],tdcCabling[4],tdcCabling[5],tdcCabling[6]);
  
  fESDZDC->SetZEM1TDChit(kFALSE);
  fESDZDC->SetZEM2TDChit(kFALSE);
  fESDZDC->SetZNCTDChit(kFALSE);
  fESDZDC->SetZPCTDChit(kFALSE);
  fESDZDC->SetZNATDChit(kFALSE);
  fESDZDC->SetZPATDChit(kFALSE);
  
  Int_t tdcValues[32][4] = {{0,}}; 
  Float_t tdcCorrected[32][4] = {{999.,}};
  for(Int_t jk=0; jk<32; jk++){
    for(Int_t lk=0; lk<4; lk++){
      tdcValues[jk][lk] = reco.GetZDCTDCData(jk, lk);
      //
      // NB -> THE ORDER MUST BE THE SAME AS THE ONE IN ZDCMAPPINGDA.cxx (7/7/2015)
      // otherwise calibrated TDC won't be centered around zero (wrong offset subtracted)
      if(tdcCabling[0]>0. && jk==tdcCabling[0] && TMath::Abs(tdcValues[jk][lk])>1e-09)      fESDZDC->SetZEM1TDChit(kTRUE);
      else if(tdcCabling[1]>0. && jk==tdcCabling[1] && TMath::Abs(tdcValues[jk][lk])>1e-09) fESDZDC->SetZEM2TDChit(kTRUE);
      else if(tdcCabling[2]>0. && jk==tdcCabling[2] && TMath::Abs(tdcValues[jk][lk])>1e-09) fESDZDC->SetZNCTDChit(kTRUE);
      else if(tdcCabling[3]>0. && jk==tdcCabling[3] && TMath::Abs(tdcValues[jk][lk])>1e-09) fESDZDC->SetZPCTDChit(kTRUE);
      else if(tdcCabling[4]>0. && jk==tdcCabling[4] && TMath::Abs(tdcValues[jk][lk])>1e-09) fESDZDC->SetZNATDChit(kTRUE);
      else if(tdcCabling[5]>0. && jk==tdcCabling[5] && TMath::Abs(tdcValues[jk][lk])>1e-09) fESDZDC->SetZPATDChit(kTRUE);
    }
  }
  //printf("\n  FillZDCintoESD: hit TDC         ZEM1 %d  ZEM2 %d ZNC %d ZPC %d ZNA %d ZPA %d \n",  fESDZDC->IsZEM1hit(),fESDZDC->IsZEM2hit(),fESDZDC->IsZNChit(),fESDZDC->IsZPChit(),fESDZDC->IsZNAhit(),fESDZDC->IsZPAhit());
  
  // Writing TDC data into ZDC ESDs
  // 4/2/2011 -> Subtracting L0 (tdcValues[15]) instead of ADC gate 
  // we try to keep the TDC oscillations as low as possible!
  for(Int_t jk=0; jk<32; jk++){
    for(Int_t lk=0; lk<4; lk++){
      if(TMath::Abs(tdcValues[jk][lk])>1e-09){
        // Feb2013 -> TDC corrected entry filled ONLY IF tdc has a hit && L0 TDC ch. is defined
        if(TMath::Abs(tdcValues[jk][lk])>1e-09 && tdcCabling[6]>0.){
	   tdcCorrected[jk][lk] = 0.025*(tdcValues[jk][lk]-tdcValues[tdcCabling[6]][0])+fMeanPhase;
           // Detector channels are centered around zero using the OCDB object
	   for(int idch=0; idch<6; idch++){
	      if(jk==tdcCabling[idch]){
	        tdcCorrected[jk][lk] -= tdcOffset[idch];
                // Ch. debug
		//printf("ch.%d -> %d (TDC raw - L0)  %f  offset %f  TDC corrected %f\n", jk,idch, 0.025*(tdcValues[jk][lk]-tdcValues[tdcCabling[6]][0]), tdcOffset[idch], tdcCorrected[jk][lk]);
              }
	   }
	}
      }
    }
  }
  
  for(int idch=0; idch<7; idch++){
     if(tdcCabling[idch]>0.) fESDZDC->SetZDCTDCChannel(idch, tdcCabling[idch]);
  }
  

  fESDZDC->SetZDCTDCData(tdcValues);
  fESDZDC->SetZDCTDCCorrected(tdcCorrected);
  fESDZDC->AliESDZDC::SetBit(AliESDZDC::kTDCcablingSet, kTRUE);
  fESDZDC->AliESDZDC::SetBit(AliESDZDC::kCorrectedTDCFilled, kTRUE);
  fESDZDC->AliESDZDC::SetBit(AliESDZDC::kEnergyCalibratedSignal, reco.GetEnergyFlag());
  //fESDZDC->Print("");
  
  if(esd) esd->SetZDCData(fESDZDC);
}

//_____________________________________________________________________________
AliCDBStorage* AliZDCReconstructor::SetStorage(const char *uri) 
{
  // Setting the storage

  Bool_t deleteManager = kFALSE;
  
  AliCDBManager *manager = AliCDBManager::Instance();
  AliCDBStorage *defstorage = manager->GetDefaultStorage();
  
  if(!defstorage || !(defstorage->Contains("ZDC"))){ 
     AliWarning("No default storage set or default storage doesn't contain ZDC!");
     manager->SetDefaultStorage(uri);
     deleteManager = kTRUE;
  }
 
  AliCDBStorage *storage = manager->GetDefaultStorage();

  if(deleteManager){
    AliCDBManager::Instance()->UnsetDefaultStorage();
    defstorage = 0;   // the storage is killed by AliCDBManager::Instance()->Destroy()
  }

  return storage; 
}

//_____________________________________________________________________________
AliZDCPedestals* AliZDCReconstructor::GetPedestalData() const
{

  // Getting pedestal calibration object for ZDC set
  AliZDCPedestals *calibdata = 0x0;
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/Pedestals");
  if(!entry) AliFatal("No calibration data loaded!");
  else{
    //entry->SetOwner(kFALSE);

    calibdata = dynamic_cast<AliZDCPedestals*>  (entry->GetObject());
    if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  }
  return calibdata;
}

//_____________________________________________________________________________
AliZDCEnCalib* AliZDCReconstructor::GetEnergyCalibData() const
{

  // Getting energy and equalization calibration object for ZDC set
  AliZDCEnCalib *calibdata = 0x0;
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/EnergyCalib");
  if(!entry) AliFatal("No calibration data loaded!");  
  else{
    entry->SetOwner(kFALSE);

    calibdata = dynamic_cast<AliZDCEnCalib*> (entry->GetObject());
    if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");
  }
  return calibdata;
}

//_____________________________________________________________________________
AliZDCSaturationCalib* AliZDCReconstructor::GetSaturationCalibData() const
{

  // Getting energy and equalization calibration object for ZDC set
  AliZDCSaturationCalib *calibdata = 0x0;
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/SaturationCalib");
  if(!entry) AliFatal("No calibration data loaded!");  
  else{
    entry->SetOwner(kFALSE);

    calibdata = dynamic_cast<AliZDCSaturationCalib*> (entry->GetObject());
    if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");
  }
  return calibdata;
}

//_____________________________________________________________________________
AliZDCTowerCalib* AliZDCReconstructor::GetTowerCalibData() const
{

  // Getting energy and equalization calibration object for ZDC set
  AliZDCTowerCalib *calibdata = 0x0;
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/TowerCalib");
  if(!entry) AliFatal("No calibration data loaded!");  
  else{
    entry->SetOwner(kFALSE);

    calibdata = dynamic_cast<AliZDCTowerCalib*> (entry->GetObject());
    if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");
  }
  return calibdata;
}

//_____________________________________________________________________________
AliZDCMBCalib* AliZDCReconstructor::GetMBCalibData() const
{

  // Getting energy and equalization calibration object for ZDC set
  AliZDCMBCalib *calibdata = 0x0;
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/MBCalib");
  if(!entry) AliFatal("No calibration data loaded!");  
  else{
    entry->SetOwner(kFALSE);

    calibdata = dynamic_cast<AliZDCMBCalib*> (entry->GetObject());
    if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");
  }
  return calibdata;
}

//_____________________________________________________________________________
AliZDCTDCCalib* AliZDCReconstructor::GetTDCCalibData() const
{

  // Getting TDC object for ZDC 
  AliZDCTDCCalib *calibdata = 0x0;
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/TDCCalib");
  if(!entry) AliFatal("No calibration data loaded!");  
  else{
    entry->SetOwner(kFALSE);

    calibdata = dynamic_cast<AliZDCTDCCalib*> (entry->GetObject());
    if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  }
  return calibdata;
}

//_____________________________________________________________________________
AliZDCChMap* AliZDCReconstructor::GetMapping() const
{

  // Getting TDC object for ZDC 
  AliZDCChMap *calibdata = 0x0;
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/ChMap");
  if(!entry) AliFatal("No calibration data loaded!");  
  else{
    entry->SetOwner(kFALSE);

    calibdata = dynamic_cast<AliZDCChMap*> (entry->GetObject());
    if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  }
  return calibdata;
}


//_____________________________________________________________________________
int AliZDCReconstructor::GetChannelSignal(int det, int quad, Bool_t intime) const
{
  if(intime){
   if(quad!=5){ // No ref. PMT
      if(det==1){
   	if(quad==0) return kZNCC;
   	else if(quad==1) return kZNC1;
   	else if(quad==2) return kZNC2;
   	else if(quad==3) return kZNC3;
   	else if(quad==4) return kZNC4;
      }
      else if(det==2){
   	if(quad==0) return kZPCC;
   	else if(quad==1) return kZPC1;
   	else if(quad==2) return kZPC2;
   	else if(quad==3) return kZPC3;
   	else if(quad==4) return kZPC4;
      }
      else if(det==3){
   	if(quad==1) return kZEM1;
   	else if(quad==2) return kZEM2;
      }
      else if(det==4){
   	if(quad==0) return kZNAC;
   	else if(quad==1) return kZNA1;
   	else if(quad==2) return kZNA2;
   	else if(quad==3) return kZNA3;
   	else if(quad==4) return kZNA4;
      }
      else if(det==5){
   	if(quad==0) return kZPAC;
   	else if(quad==1) return kZPA1;
   	else if(quad==2) return kZPA2;
   	else if(quad==3) return kZPA3;
   	else if(quad==4) return kZPA4;
      }
   }
   else{ // Ref. PMT
     if(quad==1) return kZDCCMon;
     else if(quad==4)  return kZDCAMon;
   }
  }
  else{
      if(quad!=5){ // No ref. PMT
         if(det==1){
           if(quad==0) return kZNCCoot;
           else if(quad==1) return kZNC1oot;
           else if(quad==2) return kZNC2oot;
           else if(quad==3) return kZNC3oot;
           else if(quad==4) return kZNC4oot;
         }
         else if(det==2){
           if(quad==0) return kZPCCoot;
           else if(quad==1) return kZPC1oot;
           else if(quad==2) return kZPC2oot;
           else if(quad==3) return kZPC3oot;
           else if(quad==4) return kZPC4oot;
         }
         else if(det==3){
           if(quad==1) return kZEM1oot;
           else if(quad==2) return kZEM2oot;
         }
         else if(det==4){
           if(quad==0) return kZNACoot;
           else if(quad==1) return kZNA1oot;
           else if(quad==2) return kZNA2oot;
           else if(quad==3) return kZNA3oot;
           else if(quad==4) return kZNA4oot;
         }
         else if(det==5){
           if(quad==0) return kZPACoot;
           else if(quad==1) return kZPA1oot;
           else if(quad==2) return kZPA2oot;
           else if(quad==3) return kZPA3oot;
           else if(quad==4) return kZPA4oot;
         }
      }
      else{
        if(quad==1) return kZDCCMon;
        else if(quad==4)  return kZDCAMon;
      }
    }
}
