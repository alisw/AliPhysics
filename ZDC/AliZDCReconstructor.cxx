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
//   (ZN1,ZP1) or (ZNC, ZPC) or RIGHT refers to side C (RB26)		     //
//   (ZN2,ZP2) or (ZNA, ZPA) or LEFT refers to side A (RB24)		     //
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
#include "AliZDCTowerCalib.h"
#include "AliZDCMBCalib.h"
#include "AliZDCTDCCalib.h"
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
  fTowCalibData(GetTowerCalibData()),
  fTDCCalibData(GetTDCCalibData()),
  fRecoMode(0),
  fBeamEnergy(0.),
  fNRun(0),
  fIsCalibrationMB(kFALSE),
  fPedSubMode(0),
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
   if(fPedData)      delete fPedData;    
   if(fEnCalibData)  delete fEnCalibData;
   if(fTowCalibData) delete fTowCalibData;
   if(fgMBCalibData) delete fgMBCalibData;
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
  // This is a temporary solution to allow reconstruction in tests without beam
  if(((beamType.CompareTo("UNKNOWN"))==0) && 
     ((runType.CompareTo("PHYSICS"))==0 || (runType.CompareTo("CALIBRATION_BC"))==0)){
    fRecoMode=1;
  }
  /*else if((beamType.CompareTo("UNKNOWN"))==0){
    AliError("\t UNKNOWN beam type\n");
    return;
  }*/
    
  fBeamEnergy = GetRunInfo()->GetBeamEnergy();
  if(fBeamEnergy<0.01){
     AliWarning(" Beam energy value missing -> setting it to 1380 GeV ");
     fBeamEnergy = 1380.;
  }
  
  if(((beamType.CompareTo("pp"))==0) || ((beamType.CompareTo("p-p"))==0)
     ||((beamType.CompareTo("PP"))==0) || ((beamType.CompareTo("P-P"))==0)){
    fRecoMode=1;
  }
  else if((beamType.CompareTo("A-A")) == 0 || (beamType.CompareTo("AA")) == 0){
    fRecoMode=2;
    if(!fgRecoParam) fgRecoParam = const_cast<AliZDCRecoParam*>(GetRecoParam());
    if(fgRecoParam){
      fgRecoParam->SetGlauberMCDist(fBeamEnergy);	
    } 
  }

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase"); 
  if (!entry) AliFatal("LHC clock-phase shift is not found in OCDB !");
  AliLHCClockPhase *phaseLHC = (AliLHCClockPhase*)entry->GetObject();
  // 4/2/2011 According to A. Di Mauro BEAM1 measurement is more reliable 
  // than BEAM2 and therefore also than the average of the 2
  fMeanPhase = phaseLHC->GetMeanPhaseB1();
    
  if(fIsCalibrationMB==kFALSE)  
    AliInfo(Form("\n\n ***** ZDC reconstruction initialized for %s @ %1.0f + %1.0f GeV *****\n\n",
    	beamType.Data(), fBeamEnergy, fBeamEnergy));
  
  // if EMD calibration run NO ENERGY CALIBRATION should be performed
  // pp-like reconstruction must be performed (E cailb. coeff. = 1)
  if((runType.CompareTo("CALIBRATION_EMD")) == 0){
    fRecoMode=1; 
    fBeamEnergy = 1380.;
  }
  
  AliInfo(Form("\n   ZDC reconstruction mode %d (1 -> p-p, 2-> A-A)\n\n",fRecoMode));
  
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
  else if((beamType.CompareTo("A-A")) == 0 || (beamType.CompareTo("AA")) == 0){
    fRecoMode=2;
    if(!fgRecoParam) fgRecoParam = const_cast<AliZDCRecoParam*>(GetRecoParam());
    if( fgRecoParam ) fgRecoParam->SetGlauberMCDist(fBeamEnergy);	
  }    

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase"); 
  if (!entry) AliFatal("LHC clock-phase shift is not found in OCDB !");
  AliLHCClockPhase *phaseLHC = (AliLHCClockPhase*)entry->GetObject();
  fMeanPhase = phaseLHC->GetMeanPhase();
  
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

  // get digits
  AliZDCDigit digit;
  AliZDCDigit* pdigit = &digit;
  digitsTree->SetBranchAddress("ZDC", &pdigit);
  //printf("\n\t # of digits in tree: %d\n",(Int_t) digitsTree->GetEntries());

  // loop over digits
  Float_t tZN1Corr[10], tZP1Corr[10], tZN2Corr[10], tZP2Corr[10]; 
  Float_t dZEM1Corr[2], dZEM2Corr[2], sPMRef1[2], sPMRef2[2]; 
  for(Int_t i=0; i<10; i++){
     tZN1Corr[i] = tZP1Corr[i] = tZN2Corr[i] = tZP2Corr[i] = 0.;
     if(i<2) dZEM1Corr[i] = dZEM2Corr[i] = sPMRef1[i] = sPMRef2[i] = 0.;
  }  
  
  Int_t digNentries = digitsTree->GetEntries();
  Float_t ootDigi[kNch]; Int_t i=0;
  // -- Reading out-of-time signals (last kNch entries) for current event
  if(fPedSubMode==1){
    for(Int_t iDigit=kNch; iDigit<digNentries; iDigit++){
       if(i<=kNch) ootDigi[i] = digitsTree->GetEntry(iDigit);
       else AliWarning(" Can't read more out of time values: index>kNch !!!\n");
       i++;
    }
  }
  
  for(Int_t iDigit=0; iDigit<(digNentries/2); iDigit++) {
   digitsTree->GetEntry(iDigit);
   if (!pdigit) continue;
   //  
   Int_t det = digit.GetSector(0);
   Int_t quad = digit.GetSector(1);
   Int_t pedindex = -1;
   Float_t ped2SubHg=0., ped2SubLg=0.;
   if(quad!=5){
     if(det==1)      pedindex = quad;
     else if(det==2) pedindex = quad+5;
     else if(det==3) pedindex = quad+9;
     else if(det==4) pedindex = quad+12;
     else if(det==5) pedindex = quad+17;
   }
   else pedindex = (det-1)/3+22;
   //
   if(fPedSubMode==0){
     ped2SubHg = meanPed[pedindex];
     ped2SubLg = meanPed[pedindex+kNch];
   }
   else if(fPedSubMode==1){
     ped2SubHg = corrCoeff1[pedindex]*ootDigi[pedindex]+corrCoeff0[pedindex];
     ped2SubLg = corrCoeff1[pedindex+kNch]*ootDigi[pedindex+kNch]+corrCoeff0[pedindex+kNch];
   }
      
   if(quad != 5){ // ZDC (not reference PTMs!)
    if(det == 1){ // *** ZNC
       tZN1Corr[quad] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       tZN1Corr[quad+5] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
    }
    else if(det == 2){ // *** ZP1
       tZP1Corr[quad] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       tZP1Corr[quad+5] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
    }
    else if(det == 3){
       if(quad == 1){	    // *** ZEM1  
         dZEM1Corr[0] += (Float_t) (digit.GetADCValue(0)-ped2SubHg); 
         dZEM1Corr[1] += (Float_t) (digit.GetADCValue(1)-ped2SubLg); 
       }
       else if(quad == 2){  // *** ZEM2
         dZEM2Corr[0] += (Float_t) (digit.GetADCValue(0)-ped2SubHg); 
         dZEM2Corr[1] += (Float_t) (digit.GetADCValue(1)-ped2SubLg); 
       }
    }
    else if(det == 4){  // *** ZN2
       tZN2Corr[quad] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       tZN2Corr[quad+5] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
   }
    else if(det == 5){  // *** ZP2 
       tZP2Corr[quad] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       tZP2Corr[quad+5] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
    }
   }
   else{ // Reference PMs
     if(det == 1){
       sPMRef1[0] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       sPMRef1[1] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
     }
     else if(det == 4){
       sPMRef2[0] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       sPMRef2[1] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
     }
   }

   // Ch. debug
   /*printf("AliZDCReconstructor: digit #%d det %d quad %d pedHG %1.0f pedLG %1.0f\n",
   	 iDigit, det, quad, ped2SubHg, ped2SubLg);
   printf(" -> pedindex %d\n", pedindex);
   printf("   HGChain -> RawDig %d DigCorr %1.2f", 
   	digit.GetADCValue(0), digit.GetADCValue(0)-ped2SubHg); 
   printf("   LGChain -> RawDig %d DigCorr %1.2f\n", 
   	digit.GetADCValue(1), digit.GetADCValue(1)-ped2SubLg);*/ 
   
  }//digits loop
 
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
    ReconstructEventpp(clustersTree, tZN1Corr, tZP1Corr, tZN2Corr, tZP2Corr, 
      dZEM1Corr, dZEM2Corr, sPMRef1, sPMRef2, 
      kFALSE, counts, tdc,
      evQualityBlock,  triggerBlock,  chBlock, puBits);
  else if(fRecoMode==2)
    ReconstructEventPbPb(clustersTree, tZN1Corr, tZP1Corr, tZN2Corr, tZP2Corr, 
      dZEM1Corr, dZEM2Corr, sPMRef1, sPMRef2, 
      kFALSE, counts, tdc,
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

  Int_t adcZN1[5], adcZN1oot[5], adcZN1lg[5], adcZN1ootlg[5];
  Int_t adcZP1[5], adcZP1oot[5], adcZP1lg[5], adcZP1ootlg[5];
  Int_t adcZN2[5], adcZN2oot[5], adcZN2lg[5], adcZN2ootlg[5];
  Int_t adcZP2[5], adcZP2oot[5], adcZP2lg[5], adcZP2ootlg[5];
  Int_t adcZEM[2], adcZEMoot[2], adcZEMlg[2], adcZEMootlg[2];
  Int_t pmRef[2], pmRefoot[2], pmReflg[2], pmRefootlg[2];
  for(Int_t ich=0; ich<5; ich++){
    adcZN1[ich] = adcZN1oot[ich] = adcZN1lg[ich] = adcZN1ootlg[ich] = 0;
    adcZP1[ich] = adcZP1oot[ich] = adcZP1lg[ich] = adcZP1ootlg[ich] = 0;
    adcZN2[ich] = adcZN2oot[ich] = adcZN2lg[ich] = adcZN2ootlg[ich] = 0;
    adcZP2[ich] = adcZP2oot[ich] = adcZP2lg[ich] = adcZP2ootlg[ich] = 0;
    if(ich<2){
      adcZEM[ich] = adcZEMoot[ich] = adcZEMlg[ich] = adcZEMootlg[ich] = 0;
      pmRef[ich] = pmRefoot[ich] = pmReflg[ich] = pmRefootlg[ich] = 0;
    }
  }
  
  Float_t tZN1Corr[10], tZP1Corr[10], tZN2Corr[10], tZP2Corr[10]; 
  Float_t dZEM1Corr[2], dZEM2Corr[2], sPMRef1[2], sPMRef2[2]; 
  for(Int_t i=0; i<10; i++){
     tZN1Corr[i] = tZP1Corr[i] = tZN2Corr[i] = tZP2Corr[i] = 0.;
     if(i<2) dZEM1Corr[i] = dZEM2Corr[i] = sPMRef1[i] = sPMRef2[i] = 0.;
  }  

  Bool_t isScalerOn=kFALSE;
  Int_t jsc=0, itdc=0, iprevtdc=-1, ihittdc=0;
  UInt_t scalerData[32];
  Int_t tdcData[32][4];	
  for(Int_t k=0; k<32; k++){
    scalerData[k]=0;
    for(Int_t i=0; i<4; i++) tdcData[k][i]=0;
  }
  
  
  Int_t  evQualityBlock[4] = {1,0,0,0};
  Int_t  triggerBlock[4] = {0,0,0,0};
  Int_t  chBlock[3] = {0,0,0};
  UInt_t puBits=0;

  Int_t kFirstADCGeo=0, kLastADCGeo=3, kScalerGeo=8, kZDCTDCGeo=4, kPUGeo=29;
  //Int_t kTrigScales=30, kTrigHistory=31;

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
        && (rawData.IsOverflow()==kFALSE) && (rawData.IsADCEventGood()==kTRUE)){
     
      Int_t adcMod = rawData.GetADCModule();
      Int_t det = rawData.GetSector(0);
      Int_t quad = rawData.GetSector(1);
      Int_t gain = rawData.GetADCGain();
      Int_t pedindex=0;
      //
      // Mean pedestal value subtraction -------------------------------------------------------
      if(fPedSubMode == 0){
       //  **** Pb-Pb data taking 2010 -> subtracting some ch. from correlation ****
       // Not interested in o.o.t. signals (ADC modules 2, 3)
       //if(adcMod == 2 || adcMod == 3) continue;
       //  **** Pb-Pb data taking 2011 -> subtracting only ZEM from correlation ****
       if(det==3){
	 if(adcMod==0 || adcMod==1){
	   if(gain==0) adcZEM[quad-1] = rawData.GetADCValue();
           else adcZEMlg[quad-1] = rawData.GetADCValue();
	 }
	 else if(adcMod==2 || adcMod==3){ 
	   if(gain==0) adcZEMoot[quad-1] = rawData.GetADCValue();
           else adcZEMootlg[quad-1] = rawData.GetADCValue();
	 }
       }
       // When oot values are read the ADC modules 2, 3 can be skipped!!!
       if(adcMod == 2 || adcMod == 3) continue;
       
       // *************************************************************************
       if(quad != 5){ // ZDCs (not reference PTMs)
        if(det==1){    
          pedindex = quad;
          if(gain == 0) tZN1Corr[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else tZN1Corr[quad+5]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
        }
        else if(det==2){ 
          pedindex = quad+5;
          if(gain == 0) tZP1Corr[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else tZP1Corr[quad+5]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
        }
        /*else if(det == 3){ 
          pedindex = quad+9;
          if(quad==1){     
            if(gain == 0) dZEM1Corr[0] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
            else dZEM1Corr[1] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
          }
          else if(quad==2){ 
            if(gain == 0) dZEM2Corr[0] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
            else dZEM2Corr[1] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
          }
        }*/
        else if(det == 4){	 
          pedindex = quad+12;
          if(gain == 0) tZN2Corr[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else tZN2Corr[quad+5]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
        }
        else if(det == 5){
          pedindex = quad+17;
          if(gain == 0) tZP2Corr[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else tZP2Corr[quad+5]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
        }
       }
       else{ // reference PM
         pedindex = (det-1)/3 + 22;
         if(det == 1){
           if(gain==0) sPMRef1[0] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]);
	 else sPMRef1[1] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]);
         }
         else if(det == 4){
           if(gain==0) sPMRef2[0] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]);
      	   else sPMRef2[1] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]);
         }
       }
       // Ch. debug
       /*if(gain==0){
         printf(" AliZDCReconstructor: det %d quad %d res %d -> Pedestal[%d] %1.0f", 
           det,quad,gain, pedindex, meanPed[pedindex]);
         printf("   RawADC %d ADCCorr %1.0f\n", 
           rawData.GetADCValue(), rawData.GetADCValue()-meanPed[pedindex]);
       }*/ 
      }// mean pedestal subtraction
      // Pedestal subtraction from correlation ------------------------------------------------
      else if(fPedSubMode == 1){
       // In time signals
       if(adcMod==0 || adcMod==1){
         if(quad != 5){ // signals from ZDCs
           if(det == 1){
	     if(gain==0) adcZN1[quad] = rawData.GetADCValue();
             else adcZN1lg[quad] = rawData.GetADCValue();
	   }
	   else if(det == 2){
	     if(gain==0) adcZP1[quad] = rawData.GetADCValue();
             else adcZP1lg[quad] = rawData.GetADCValue();
	   }
	   else if(det == 3){
	     if(gain==0) adcZEM[quad-1] = rawData.GetADCValue();
             else adcZEMlg[quad-1] = rawData.GetADCValue();
	   }
	   else if(det == 4){
	     if(gain==0) adcZN2[quad] = rawData.GetADCValue();
             else adcZN2lg[quad] = rawData.GetADCValue();
	   }
	   else if(det == 5){
	     if(gain==0) adcZP2[quad] = rawData.GetADCValue();
             else adcZP2lg[quad] = rawData.GetADCValue();
	   }
	 }
	 else{ // signals from reference PM
	    if(gain==0) pmRef[quad-1] = rawData.GetADCValue();
            else pmReflg[quad-1] = rawData.GetADCValue();
	 }
       }
       // Out-of-time pedestals
       else if(adcMod==2 || adcMod==3){
         if(quad != 5){ // signals from ZDCs
           if(det == 1){
	     if(gain==0) adcZN1oot[quad] = rawData.GetADCValue();
             else adcZN1ootlg[quad] = rawData.GetADCValue();
	   }
	   else if(det == 2){
	     if(gain==0) adcZP1oot[quad] = rawData.GetADCValue();
             else adcZP1ootlg[quad] = rawData.GetADCValue();
	   }
	   else if(det == 3){
	     if(gain==0) adcZEMoot[quad-1] = rawData.GetADCValue();
             else adcZEMootlg[quad-1] = rawData.GetADCValue();
	   }
	   else if(det == 4){
	     if(gain==0) adcZN2oot[quad] = rawData.GetADCValue();
             else adcZN2ootlg[quad] = rawData.GetADCValue();
	   }
	   else if(det == 5){
	     if(gain==0) adcZP2oot[quad] = rawData.GetADCValue();
             else adcZP2ootlg[quad] = rawData.GetADCValue();
	   }
	 }
	 else{ // signals from reference PM
	    if(gain==0) pmRefoot[quad-1] = rawData.GetADCValue();
            else pmRefootlg[quad-1] = rawData.GetADCValue();
	 }
       }
      } // pedestal subtraction from correlation
      // Ch. debug
      /*printf("\t AliZDCReconstructor: det %d quad %d res %d -> Ped[%d] = %1.0f\n", 
        det,quad,gain, pedindex, meanPed[pedindex]);*/
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
       if(itdc==iprevtdc) ihittdc++;
       else ihittdc=0;
       iprevtdc=itdc;
       if(ihittdc<4) tdcData[itdc][ihittdc] = rawData.GetZDCTDCDatum();
       // Ch. debug
       //if(ihittdc==0) printf("   TDC%d %d  ",itdc, tdcData[itdc][ihittdc]);
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
  
  if(fPedSubMode==1){
    for(Int_t t=0; t<5; t++){
       tZN1Corr[t] = adcZN1[t] - (corrCoeff1[t]*adcZN1oot[t]+corrCoeff0[t]);
       tZN1Corr[t+5] = adcZN1lg[t] - (corrCoeff1[t+kNch]*adcZN1ootlg[t]+corrCoeff0[t+kNch]);
       //
       tZP1Corr[t] = adcZP1[t] - (corrCoeff1[t+5]*adcZP1oot[t]+corrCoeff0[t+5]);
       tZP1Corr[t+5] = adcZP1lg[t] - (corrCoeff1[t+5+kNch]*adcZP1ootlg[t]+corrCoeff0[t+5+kNch]);
       //
       tZN2Corr[t] = adcZN2[t] - (corrCoeff1[t+12]*adcZN2oot[t]+corrCoeff0[t+12]);
       tZN2Corr[t+5] = adcZN2lg[t] - (corrCoeff1[t+12+kNch]*adcZN2ootlg[t]+corrCoeff0[t+12+kNch]);
       //
       tZP2Corr[t] = adcZP2[t] - (corrCoeff1[t+17]*adcZP2oot[t]+corrCoeff0[t+17]);
       tZP2Corr[t+5] = adcZP2lg[t] - (corrCoeff1[t+17+kNch]*adcZP2ootlg[t]+corrCoeff0[t+17+kNch]);
    }
    dZEM1Corr[0] = adcZEM[0]   - (corrCoeff1[10]*adcZEMoot[0]+corrCoeff0[10]);
    dZEM1Corr[1] = adcZEMlg[0] - (corrCoeff1[10+kNch]*adcZEMootlg[0]+corrCoeff0[10+kNch]);
    dZEM2Corr[0] = adcZEM[1]   - (corrCoeff1[11]*adcZEMoot[1]+corrCoeff0[11]);
    dZEM2Corr[1] = adcZEMlg[1] - (corrCoeff1[11+kNch]*adcZEMootlg[1]+corrCoeff0[11+kNch]);
    //
    sPMRef1[0] = pmRef[0]   - (corrCoeff1[22]*pmRefoot[0]+corrCoeff0[22]);
    sPMRef1[1] = pmReflg[0] - (corrCoeff1[22+kNch]*pmRefootlg[0]+corrCoeff0[22+kNch]);
    sPMRef2[0] = pmRef[0]   - (corrCoeff1[23]*pmRefoot[1]+corrCoeff0[23]);
    sPMRef2[1] = pmReflg[0] - (corrCoeff1[23+kNch]*pmRefootlg[1]+corrCoeff0[23+kNch]);
  }
  if(fPedSubMode==0 && fRecoMode==2){
    //  **** Pb-Pb data taking 2011 -> subtracting some ch. from correlation ****
    //tZN1Corr[0] = adcZN1[0] - (corrCoeff1[0]*adcZN1oot[0]+corrCoeff0[0]);
    //tZN1Corr[5] = adcZN1lg[0] - (corrCoeff1[kNch]*adcZN1ootlg[0]+corrCoeff0[kNch]);
    // Ch. debug
    //printf(" adcZN1 %d  adcZN1oot %d tZN1Corr %1.2f \n", adcZN1[0],adcZN1oot[0],tZN1Corr[0]);
    //printf(" adcZN1lg %d  adcZN1ootlg %d tZN1Corrlg %1.2f \n", adcZN1lg[0],adcZN1ootlg[0],tZN1Corr[5]);
    //
    //tZP1Corr[2] = adcZP1[2] - (corrCoeff1[2+5]*adcZP1oot[2]+corrCoeff0[2+5]);
    //tZP1Corr[2+5] = adcZP1lg[2] - (corrCoeff1[2+5+kNch]*adcZP1ootlg[2]+corrCoeff0[2+5+kNch]);
    //
    dZEM1Corr[0] = adcZEM[0]   - (corrCoeff1[10]*adcZEMoot[0]+corrCoeff0[10]);
    dZEM1Corr[1] = adcZEMlg[0] - (corrCoeff1[10+kNch]*adcZEMootlg[0]+corrCoeff0[10+kNch]);
    dZEM2Corr[0] = adcZEM[1]   - (corrCoeff1[11]*adcZEMoot[1]+corrCoeff0[11]);
    dZEM2Corr[1] = adcZEMlg[1] - (corrCoeff1[11+kNch]*adcZEMootlg[1]+corrCoeff0[11+kNch]);
    // *************************************************************************
  }
  else if(fPedSubMode==0 && fRecoMode==1){
    //  **** p-p data taking 2011 -> temporary patch to overcome DA problem ****
    tZN1Corr[0] = adcZN1[0] - meanPed[0];
    tZN1Corr[5] = adcZN1lg[0] - meanPed[kNch];
    //
    dZEM1Corr[0] = adcZEM[0]   - meanPed[10];
    dZEM1Corr[1] = adcZEMlg[0] - meanPed[10+kNch];
    dZEM2Corr[0] = adcZEM[1]   - meanPed[11];
    dZEM2Corr[1] = adcZEMlg[1] - meanPed[11+kNch];
        // *************************************************************************
  }
    
  if(fRecoMode==1) // p-p data
    ReconstructEventpp(clustersTree, tZN1Corr, tZP1Corr, tZN2Corr, tZP2Corr, 
      dZEM1Corr, dZEM2Corr, sPMRef1, sPMRef2, 
      isScalerOn, scalerData, tdcData,
      evQualityBlock, triggerBlock, chBlock, puBits);
  else if(fRecoMode==2) // Pb-Pb data
      ReconstructEventPbPb(clustersTree, tZN1Corr, tZP1Corr, tZN2Corr, tZP2Corr, 
      dZEM1Corr, dZEM2Corr, sPMRef1, sPMRef2, 
      isScalerOn, scalerData,  tdcData,
      evQualityBlock, triggerBlock, chBlock, puBits);
}

//_____________________________________________________________________________
void AliZDCReconstructor::ReconstructEventpp(TTree *clustersTree, 
	const Float_t* const corrADCZN1, const Float_t* const corrADCZP1, 
	const Float_t* const corrADCZN2, const Float_t* const corrADCZP2,
	const Float_t* const corrADCZEM1, const Float_t* const corrADCZEM2,
	Float_t* sPMRef1, Float_t* sPMRef2, Bool_t isScalerOn, UInt_t* scaler, 
	Int_t tdcData[32][4], const Int_t* const evQualityBlock, 
	const Int_t* const triggerBlock, const Int_t* const chBlock, UInt_t puBits) const
{
  // ****************** Reconstruct one event ******************
  
  // CH. debug
  /*printf("\n*************************************************\n");
  printf(" ReconstructEventpp -> values after pedestal subtraction:\n");
  printf(" ADCZN1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	corrADCZN1[0],corrADCZN1[1],corrADCZN1[2],corrADCZN1[3],corrADCZN1[4]);
  printf(" ADCZP1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	corrADCZP1[0],corrADCZP1[1],corrADCZP1[2],corrADCZP1[3],corrADCZP1[4]);
  printf(" ADCZN2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	corrADCZN2[0],corrADCZN2[1],corrADCZN2[2],corrADCZN2[3],corrADCZN2[4]);
  printf(" ADCZP2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	corrADCZP2[0],corrADCZP2[1],corrADCZP2[2],corrADCZP2[3],corrADCZP2[4]);
  printf(" ADCZEM1 [%1.2f] ADCZEM2 [%1.2f] \n",corrADCZEM1[0],corrADCZEM2[0]);
  printf("*************************************************\n");*/
    
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
  
  if(corrADCZP1[0]>fSignalThreshold)  rFlags[5] = 0x1;
  if(corrADCZN1[0]>fSignalThreshold)  rFlags[4] = 0x1;
  if(corrADCZEM2[0]>fSignalThreshold) rFlags[3] = 0x1;
  if(corrADCZEM1[0]>fSignalThreshold) rFlags[2] = 0x1;
  if(corrADCZP2[0]>fSignalThreshold)  rFlags[1] = 0x1;
  if(corrADCZN2[0]>fSignalThreshold)  rFlags[0] = 0x1;

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
  Float_t equalCoeffZN1[5], equalCoeffZP1[5], equalCoeffZN2[5], equalCoeffZP2[5];
  for(Int_t ji=0; ji<5; ji++){
     equalCoeffZN1[ji] = fTowCalibData->GetZN1EqualCoeff(ji);
     equalCoeffZP1[ji] = fTowCalibData->GetZP1EqualCoeff(ji); 
     equalCoeffZN2[ji] = fTowCalibData->GetZN2EqualCoeff(ji); 
     equalCoeffZP2[ji] = fTowCalibData->GetZP2EqualCoeff(ji); 
  }
  // --- Energy calibration factors ------------------------------------
  Float_t calibEne[6];
  // **** Energy calibration coefficient set to 1 
  // **** (no trivial way to calibrate in p-p runs)
  for(Int_t ij=0; ij<6; ij++) calibEne[ij] = fEnCalibData->GetEnCalib(ij);
  
  // ******	Equalization of detector responses
  Float_t equalTowZN1[10], equalTowZN2[10], equalTowZP1[10], equalTowZP2[10];
  for(Int_t gi=0; gi<10; gi++){
     if(gi<5){
       equalTowZN1[gi] = corrADCZN1[gi]*equalCoeffZN1[gi];
       equalTowZP1[gi] = corrADCZP1[gi]*equalCoeffZP1[gi];
       equalTowZN2[gi] = corrADCZN2[gi]*equalCoeffZN2[gi];
       equalTowZP2[gi] = corrADCZP2[gi]*equalCoeffZP2[gi];
     }
     else{
       equalTowZN1[gi] = corrADCZN1[gi]*equalCoeffZN1[gi-5];
       equalTowZP1[gi] = corrADCZP1[gi]*equalCoeffZP1[gi-5];
       equalTowZN2[gi] = corrADCZN2[gi]*equalCoeffZN2[gi-5];
       equalTowZP2[gi] = corrADCZP2[gi]*equalCoeffZP2[gi-5];
     }
  }
  // Ch. debug
  /*printf("\n ------------- EQUALIZATION -------------\n");
  printf(" ADCZN1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZN1[0],equalTowZN1[1],equalTowZN1[2],equalTowZN1[3],equalTowZN1[4]);
  printf(" ADCZP1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZP1[0],equalTowZP1[1],equalTowZP1[2],equalTowZP1[3],equalTowZP1[4]);
  printf(" ADCZN2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZN2[0],equalTowZN2[1],equalTowZN2[2],equalTowZN2[3],equalTowZN2[4]);
  printf(" ADCZP2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZP2[0],equalTowZP2[1],equalTowZP2[2],equalTowZP2[3],equalTowZP2[4]);
  printf(" ----------------------------------------\n");*/
  
  // ******	Summed response for hadronic calorimeter (SUMMED and then CALIBRATED!)
  Float_t calibSumZN1[]={0,0}, calibSumZN2[]={0,0}, calibSumZP1[]={0,0}, calibSumZP2[]={0,0};
  for(Int_t gi=0; gi<5; gi++){
       calibSumZN1[0] += equalTowZN1[gi];
       calibSumZP1[0] += equalTowZP1[gi];
       calibSumZN2[0] += equalTowZN2[gi];
       calibSumZP2[0] += equalTowZP2[gi];
       //
       calibSumZN1[1] += equalTowZN1[gi+5];
       calibSumZP1[1] += equalTowZP1[gi+5];
       calibSumZN2[1] += equalTowZN2[gi+5];
       calibSumZP2[1] += equalTowZP2[gi+5];
  }
  // High gain chain
  calibSumZN1[0] = calibSumZN1[0]*calibEne[0];
  calibSumZP1[0] = calibSumZP1[0]*calibEne[1];
  calibSumZN2[0] = calibSumZN2[0]*calibEne[2];
  calibSumZP2[0] = calibSumZP2[0]*calibEne[3];
  // Low gain chain
  calibSumZN1[1] = calibSumZN1[1]*calibEne[0];
  calibSumZP1[1] = calibSumZP1[1]*calibEne[1];
  calibSumZN2[1] = calibSumZN2[1]*calibEne[2];
  calibSumZP2[1] = calibSumZP2[1]*calibEne[3];
  
  // ******	Energy calibration of detector responses
  Float_t calibTowZN1[10], calibTowZN2[10], calibTowZP1[10], calibTowZP2[10];
  for(Int_t gi=0; gi<5; gi++){
     // High gain chain
     calibTowZN1[gi] = equalTowZN1[gi]*calibEne[0];
     calibTowZP1[gi] = equalTowZP1[gi]*calibEne[1];
     calibTowZN2[gi] = equalTowZN2[gi]*calibEne[2];
     calibTowZP2[gi] = equalTowZP2[gi]*calibEne[3];
     // Low gain chain
     calibTowZN1[gi+5] = equalTowZN1[gi+5]*calibEne[0];
     calibTowZP1[gi+5] = equalTowZP1[gi+5]*calibEne[1];
     calibTowZN2[gi+5] = equalTowZN2[gi+5]*calibEne[2];
     calibTowZP2[gi+5] = equalTowZP2[gi+5]*calibEne[3];
  }
  //
  Float_t sumZEM[]={0,0}, calibZEM1[]={0,0}, calibZEM2[]={0,0};
  calibZEM1[0] = corrADCZEM1[0]*calibEne[4];
  calibZEM1[1] = corrADCZEM1[1]*calibEne[4];
  calibZEM2[0] = corrADCZEM2[0]*calibEne[5];
  calibZEM2[1] = corrADCZEM2[1]*calibEne[5];
  for(Int_t k=0; k<2; k++) sumZEM[k] = calibZEM1[k] + calibZEM2[k];
  // Ch. debug
  /*printf("\n ------------- CALIBRATION -------------\n");
  printf(" ADCZN1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZN1[0],calibTowZN1[1],calibTowZN1[2],calibTowZN1[3],calibTowZN1[4]);
  printf(" ADCZP1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZP1[0],calibTowZP1[1],calibTowZP1[2],calibTowZP1[3],calibTowZP1[4]);
  printf(" ADCZN2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZN2[0],calibTowZN2[1],calibTowZN2[2],calibTowZN2[3],calibTowZN2[4]);
  printf(" ADCZP2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZP2[0],calibTowZP2[1],calibTowZP2[2],calibTowZP2[3],calibTowZP2[4]);
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
  AliZDCReco* reco = new AliZDCReco(calibSumZN1, calibSumZP1, calibSumZN2, calibSumZP2, 
  		   calibTowZN1, calibTowZP1, calibTowZN2, calibTowZP2, 
		   calibZEM1, calibZEM2, sPMRef1, sPMRef2,
		   nDetSpecNLeft, nDetSpecPLeft, nDetSpecNRight, nDetSpecPRight, 
		   nGenSpec, nGenSpecLeft, nGenSpecRight, 
		   nPart, nPartTotLeft, nPartTotRight, 
		   impPar, impPar1, impPar2,
		   recoFlag, energyFlag, isScalerOn, scaler, tdcData);
		  
  const Int_t kBufferSize = 4000;
  clustersTree->Branch("ZDC", "AliZDCReco", &reco, kBufferSize);
  // write the output tree
  clustersTree->Fill();
  delete reco;
}

//_____________________________________________________________________________
void AliZDCReconstructor::ReconstructEventPbPb(TTree *clustersTree, 
	const Float_t* const corrADCZN1, const Float_t* const corrADCZP1, 
	const Float_t* const corrADCZN2, const Float_t* const corrADCZP2,
	const Float_t* const corrADCZEM1, const Float_t* const corrADCZEM2,
	Float_t* sPMRef1, Float_t* sPMRef2, Bool_t isScalerOn, UInt_t* scaler, 
	Int_t tdcData[32][4], const Int_t* const evQualityBlock, 
	const Int_t* const triggerBlock, const Int_t* const chBlock, UInt_t puBits) const
{
  // ****************** Reconstruct one event ******************
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
  
  if(corrADCZP1[0]>fSignalThreshold)  rFlags[5] = 0x1;
  if(corrADCZN1[0]>fSignalThreshold)  rFlags[4] = 0x1;
  if(corrADCZEM2[0]>fSignalThreshold) rFlags[3] = 0x1;
  if(corrADCZEM1[0]>fSignalThreshold) rFlags[2] = 0x1;
  if(corrADCZP2[0]>fSignalThreshold)  rFlags[1] = 0x1;
  if(corrADCZN2[0]>fSignalThreshold)  rFlags[0] = 0x1;

  UInt_t recoFlag = rFlags[31] << 31 | rFlags[30] << 30 | rFlags[29] << 29 | rFlags[28] << 28 |
             rFlags[27] << 27 | rFlags[26] << 26 | rFlags[25] << 25 | rFlags[24] << 24 |
	     0x0 << 23 | 0x0 << 22 | 0x0 << 21 | 0x0 << 20 |
	     0x0 << 19 | rFlags[18] << 18 |  rFlags[17] << 17 |  rFlags[16] << 16 |
	     0x0 << 15 | 0x0 << 14 | rFlags[13] << 13 | rFlags[12] << 12 | 
             rFlags[11] << 11 |rFlags[10] << 10 | rFlags[9] << 9 | rFlags[8] << 8 |
	     0x0 << 7 | 0x0 << 6 | rFlags[5] << 5 | rFlags[4] << 4 | 
	     rFlags[3] << 3 | rFlags[2] << 2 | rFlags[1] << 1 | rFlags[0];
  // --------------------------------------------------
  
  
  // CH. debug
/*  printf("\n*************************************************\n");
  printf(" ReconstructEventPbPb -> values after pedestal subtraction:\n");
  printf(" ADCZN1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	corrADCZN1[0],corrADCZN1[1],corrADCZN1[2],corrADCZN1[3],corrADCZN1[4]);
  printf(" ADCZP1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	corrADCZP1[0],corrADCZP1[1],corrADCZP1[2],corrADCZP1[3],corrADCZP1[4]);
  printf(" ADCZN2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	corrADCZN2[0],corrADCZN2[1],corrADCZN2[2],corrADCZN2[3],corrADCZN2[4]);
  printf(" ADCZP2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	corrADCZP2[0],corrADCZP2[1],corrADCZP2[2],corrADCZP2[3],corrADCZP2[4]);
  printf(" ADCZEM1 [%1.2f] ADCZEM2 [%1.2f] \n",corrADCZEM1[0],corrADCZEM2[0]);
  printf("*************************************************\n");
*/
  // ******	Retrieving calibration data 
  // --- Equalization coefficients ---------------------------------------------
  Float_t equalCoeffZN1[5], equalCoeffZP1[5], equalCoeffZN2[5], equalCoeffZP2[5];
  for(Int_t ji=0; ji<5; ji++){
     equalCoeffZN1[ji] = fTowCalibData->GetZN1EqualCoeff(ji);
     equalCoeffZP1[ji] = fTowCalibData->GetZP1EqualCoeff(ji); 
     equalCoeffZN2[ji] = fTowCalibData->GetZN2EqualCoeff(ji); 
     equalCoeffZP2[ji] = fTowCalibData->GetZP2EqualCoeff(ji); 
  }
  // --- Energy calibration factors ------------------------------------
  Float_t calibEne[6];
  // The energy calibration object already takes into account of E_beam 
  // -> the value from the OCDB can be directly used (Jul 2010)
  for(Int_t ij=0; ij<6; ij++) calibEne[ij] = fEnCalibData->GetEnCalib(ij);
  
  // ******	Equalization of detector responses
  Float_t equalTowZN1[10], equalTowZN2[10], equalTowZP1[10], equalTowZP2[10];
  for(Int_t gi=0; gi<10; gi++){
     if(gi<5){
       equalTowZN1[gi] = corrADCZN1[gi]*equalCoeffZN1[gi];
       equalTowZP1[gi] = corrADCZP1[gi]*equalCoeffZP1[gi];
       equalTowZN2[gi] = corrADCZN2[gi]*equalCoeffZN2[gi];
       equalTowZP2[gi] = corrADCZP2[gi]*equalCoeffZP2[gi];
     }
     else{
       equalTowZN1[gi] = corrADCZN1[gi]*equalCoeffZN1[gi-5];
       equalTowZP1[gi] = corrADCZP1[gi]*equalCoeffZP1[gi-5];
       equalTowZN2[gi] = corrADCZN2[gi]*equalCoeffZN2[gi-5];
       equalTowZP2[gi] = corrADCZP2[gi]*equalCoeffZP2[gi-5];
     }
  }
  
  // Ch. debug
/*  printf("\n ------------- EQUALIZATION -------------\n");
  printf(" ADCZN1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZN1[0],equalTowZN1[1],equalTowZN1[2],equalTowZN1[3],equalTowZN1[4]);
  printf(" ADCZP1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZP1[0],equalTowZP1[1],equalTowZP1[2],equalTowZP1[3],equalTowZP1[4]);
  printf(" ADCZN2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZN2[0],equalTowZN2[1],equalTowZN2[2],equalTowZN2[3],equalTowZN2[4]);
  printf(" ADCZP2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	equalTowZP2[0],equalTowZP2[1],equalTowZP2[2],equalTowZP2[3],equalTowZP2[4]);
  printf(" ----------------------------------------\n");
*/
  
  // ******	Summed response for hadronic calorimeter (SUMMED and then CALIBRATED!)
  Float_t calibSumZN1[]={0,0}, calibSumZN2[]={0,0}, calibSumZP1[]={0,0}, calibSumZP2[]={0,0};
  for(Int_t gi=0; gi<5; gi++){
       calibSumZN1[0] += equalTowZN1[gi];
       calibSumZP1[0] += equalTowZP1[gi];
       calibSumZN2[0] += equalTowZN2[gi];
       calibSumZP2[0] += equalTowZP2[gi];
       //
       calibSumZN1[1] += equalTowZN1[gi+5];
       calibSumZP1[1] += equalTowZP1[gi+5];
       calibSumZN2[1] += equalTowZN2[gi+5];
       calibSumZP2[1] += equalTowZP2[gi+5];
  }
  //
  //fEnCalibData->Print("");
  
  // High gain chain
  calibSumZN1[0] = calibSumZN1[0]*calibEne[0]*8.;
  calibSumZP1[0] = calibSumZP1[0]*calibEne[1]*8.;
  calibSumZN2[0] = calibSumZN2[0]*calibEne[2]*8.;
  calibSumZP2[0] = calibSumZP2[0]*calibEne[3]*8.;
  // Low gain chain
  calibSumZN1[1] = calibSumZN1[1]*calibEne[0];
  calibSumZP1[1] = calibSumZP1[1]*calibEne[1];
  calibSumZN2[1] = calibSumZN2[1]*calibEne[2];
  calibSumZP2[1] = calibSumZP2[1]*calibEne[3];
  //
  Float_t sumZEM[]={0,0}, calibZEM1[]={0,0}, calibZEM2[]={0,0};
  calibZEM1[0] = corrADCZEM1[0]*calibEne[4]*8.;
  calibZEM1[1] = corrADCZEM1[1]*calibEne[4];
  calibZEM2[0] = corrADCZEM2[0]*calibEne[5]*8.;
  calibZEM2[1] = corrADCZEM2[1]*calibEne[5];
  for(Int_t k=0; k<2; k++) sumZEM[k] = calibZEM1[k] + calibZEM2[k];
    
  // ******	Energy calibration of detector responses
  Float_t calibTowZN1[10], calibTowZN2[10], calibTowZP1[10], calibTowZP2[10];
  for(Int_t gi=0; gi<5; gi++){
     // High gain chain
     calibTowZN1[gi] = equalTowZN1[gi]*2*calibEne[0]*8.;
     calibTowZP1[gi] = equalTowZP1[gi]*2*calibEne[1]*8.;
     calibTowZN2[gi] = equalTowZN2[gi]*2*calibEne[2]*8.;
     calibTowZP2[gi] = equalTowZP2[gi]*2*calibEne[3]*8.;
     // Low gain chain
     calibTowZN1[gi+5] = equalTowZN1[gi+5]*2*calibEne[0];
     calibTowZP1[gi+5] = equalTowZP1[gi+5]*2*calibEne[1];
     calibTowZN2[gi+5] = equalTowZN2[gi+5]*2*calibEne[2];
     calibTowZP2[gi+5] = equalTowZP2[gi+5]*2*calibEne[3];
  }

  // Ch. debug
/*  printf("\n ------------- CALIBRATION -------------\n");
  printf(" ADCZN1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZN1[0],calibTowZN1[1],calibTowZN1[2],calibTowZN1[3],calibTowZN1[4]);
  printf(" ADCZP1 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZP1[0],calibTowZP1[1],calibTowZP1[2],calibTowZP1[3],calibTowZP1[4]);
  printf(" ADCZN2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZN2[0],calibTowZN2[1],calibTowZN2[2],calibTowZN2[3],calibTowZN2[4]);
  printf(" ADCZP2 [%1.2f %1.2f %1.2f %1.2f %1.2f]\n",
  	calibTowZP2[0],calibTowZP2[1],calibTowZP2[2],calibTowZP2[3],calibTowZP2[4]);
  printf(" ADCZEM1 [%1.2f] ADCZEM2 [%1.2f] \n",calibZEM1[0],calibZEM2[0]);
  printf(" ----------------------------------------\n");
*/  
  //  ******	Number of detected spectator nucleons
  Int_t nDetSpecNLeft=0, nDetSpecPLeft=0, nDetSpecNRight=0, nDetSpecPRight=0;
  if(fBeamEnergy>0.01){
    nDetSpecNLeft = (Int_t) (calibSumZN1[0]/fBeamEnergy);
    nDetSpecPLeft = (Int_t) (calibSumZP1[0]/fBeamEnergy);
    nDetSpecNRight = (Int_t) (calibSumZN2[0]/fBeamEnergy);
    nDetSpecPRight = (Int_t) (calibSumZP2[0]/fBeamEnergy);
  }
  else AliWarning(" ATTENTION!!! fBeamEnergy=0 -> N_spec will be ZERO!!! \n");
  /*printf("\n\t AliZDCReconstructor -> fBeamEnergy %1.0f: nDetSpecNsideA %d, nDetSpecPsideA %d,"
    " nDetSpecNsideC %d, nDetSpecPsideC %d\n",fBeamEnergy,nDetSpecNLeft, nDetSpecPLeft, 
    nDetSpecNRight, nDetSpecPRight);*/
  
  Int_t nGenSpec=0, nGenSpecA=0, nGenSpecC=0;
  Int_t nPart=0, nPartA=0, nPartC=0;
  Double_t b=0., bA=0., bC=0.;
  
  if(fIsCalibrationMB == kFALSE){
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
    Float_t y = (calibSumZN1[0]+calibSumZP1[0]+calibSumZN2[0]+calibSumZP2[0])/1000.;
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
    Float_t yC = (calibSumZN1[0]+calibSumZP1[0])/1000.;
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
    Float_t yA = (calibSumZN2[0]+calibSumZP2[0])/1000.;
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
  } // ONLY IF fIsCalibrationMB==kFALSE
  
  Bool_t energyFlag = kTRUE;  
  AliZDCReco* reco = new AliZDCReco(calibSumZN1, calibSumZP1, calibSumZN2, calibSumZP2, 
  	  	  calibTowZN1, calibTowZP1, calibTowZN2, calibTowZP2, 
		  calibZEM1, calibZEM2, sPMRef1, sPMRef2,
		  nDetSpecNLeft, nDetSpecPLeft, nDetSpecNRight, nDetSpecPRight, 
		  nGenSpec, nGenSpecA, nGenSpecC, 
		  nPart, nPartA, nPartC, b, bA, bC,
		  recoFlag, energyFlag, isScalerOn, scaler, tdcData);
		    
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
  //
  Float_t tZN1Ene[5], tZN2Ene[5], tZP1Ene[5], tZP2Ene[5];
  Float_t tZN1EneLR[5], tZN2EneLR[5], tZP1EneLR[5], tZP2EneLR[5];
  for(Int_t i=0; i<5; i++){
     tZN1Ene[i] = reco.GetZN1HREnTow(i);
     tZN2Ene[i] = reco.GetZN2HREnTow(i);
     tZP1Ene[i] = reco.GetZP1HREnTow(i);
     tZP2Ene[i] = reco.GetZP2HREnTow(i);
     //
     tZN1EneLR[i] = reco.GetZN1LREnTow(i);
     tZN2EneLR[i] = reco.GetZN2LREnTow(i);
     tZP1EneLR[i] = reco.GetZP1LREnTow(i);
     tZP2EneLR[i] = reco.GetZP2LREnTow(i);
  }
  //
  fESDZDC->SetZN1TowerEnergy(tZN1Ene);
  fESDZDC->SetZN2TowerEnergy(tZN2Ene);
  fESDZDC->SetZP1TowerEnergy(tZP1Ene);
  fESDZDC->SetZP2TowerEnergy(tZP2Ene);
  //
  fESDZDC->SetZN1TowerEnergyLR(tZN1EneLR);
  fESDZDC->SetZN2TowerEnergyLR(tZN2EneLR);
  fESDZDC->SetZP1TowerEnergyLR(tZP1EneLR);
  fESDZDC->SetZP2TowerEnergyLR(tZP2EneLR);
  // 
  Int_t nPart  = reco.GetNParticipants();
  Int_t nPartA = reco.GetNPartSideA();
  Int_t nPartC = reco.GetNPartSideC();
  Double_t b  = reco.GetImpParameter();
  Double_t bA = reco.GetImpParSideA();
  Double_t bC = reco.GetImpParSideC();
  UInt_t recoFlag = reco.GetRecoFlag();
  
  fESDZDC->SetZDC(reco.GetZN1HREnergy(), reco.GetZP1HREnergy(), 
  	reco.GetZEM1HRsignal(), reco.GetZEM2HRsignal(), 
	reco.GetZN2HREnergy(), reco.GetZP2HREnergy(), 
	nPart, nPartA, nPartC, b, bA, bC, recoFlag);
  
  // Writing ZDC scaler for cross section calculation
  // ONLY IF the scaler has been read during the event
  if(reco.IsScalerOn()==kTRUE){
    UInt_t counts[32];
    for(Int_t jk=0; jk<32; jk++) counts[jk] = reco.GetZDCScaler(jk);
    fESDZDC->SetZDCScaler(counts);
  }    
  
  Int_t tdcValues[32][4]; 
  Float_t tdcCorrected[32][4];
  for(Int_t jk=0; jk<32; jk++){
    for(Int_t lk=0; lk<4; lk++){
      tdcValues[jk][lk] = reco.GetZDCTDCData(jk, lk);
      //Ch debug
      //if((jk>=8 && jk<=13 && lk==0) || jk==15) printf(" *** ZDC: tdc%d =  %d = %f ns \n",jk,tdcValues[jk][lk],0.025*tdcValues[jk][lk]);
    }
  }
  
  // Writing TDC data into ZDC ESDs
  // 4/2/2011 -> Subtracting L0 (tdcValues[15]) instead of ADC gate 
  // we try to keep the TDC oscillations as low as possible!
  for(Int_t jk=0; jk<32; jk++){
    for(Int_t lk=0; lk<4; lk++){
      if(TMath::Abs(tdcValues[jk][lk])>1e-10){
        tdcCorrected[jk][lk] = 0.025*(tdcValues[jk][lk]-tdcValues[15][0])+fMeanPhase;
        // Sep 2011: TDC ch. from 8 to 13 centered around 0 using OCDB 
	if(jk>=8 && jk<=13) tdcCorrected[jk][lk] =  tdcCorrected[jk][lk] - tdcOffset[jk-8];
	//Ch. debug
	//if(jk>=8 && jk<=13) printf(" *** tdcOffset%d %f  tdcCorr%d %f \n",jk,tdcOffset[jk-8],tdcCorrected[jk][lk]);
   
      }
    }
  }

  fESDZDC->SetZDCTDCData(tdcValues);
  fESDZDC->SetZDCTDCCorrected(tdcCorrected);
  fESDZDC->AliESDZDC::SetBit(AliESDZDC::kCorrectedTDCFilled, reco.GetEnergyFlag());
  fESDZDC->AliESDZDC::SetBit(AliESDZDC::kEnergyCalibratedSignal, kTRUE);
  
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

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/Pedestals");
  if(!entry) AliFatal("No calibration data loaded!");
  entry->SetOwner(kFALSE);

  AliZDCPedestals *calibdata = dynamic_cast<AliZDCPedestals*>  (entry->GetObject());
  if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}

//_____________________________________________________________________________
AliZDCEnCalib* AliZDCReconstructor::GetEnergyCalibData() const
{

  // Getting energy and equalization calibration object for ZDC set

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/EnergyCalib");
  if(!entry) AliFatal("No calibration data loaded!");  
  entry->SetOwner(kFALSE);

  AliZDCEnCalib *calibdata = dynamic_cast<AliZDCEnCalib*> (entry->GetObject());
  if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}

//_____________________________________________________________________________
AliZDCTowerCalib* AliZDCReconstructor::GetTowerCalibData() const
{

  // Getting energy and equalization calibration object for ZDC set

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/TowerCalib");
  if(!entry) AliFatal("No calibration data loaded!");  
  entry->SetOwner(kFALSE);

  AliZDCTowerCalib *calibdata = dynamic_cast<AliZDCTowerCalib*> (entry->GetObject());
  if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}

//_____________________________________________________________________________
AliZDCMBCalib* AliZDCReconstructor::GetMBCalibData() const
{

  // Getting energy and equalization calibration object for ZDC set

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/MBCalib");
  if(!entry) AliFatal("No calibration data loaded!");  
  entry->SetOwner(kFALSE);

  AliZDCMBCalib *calibdata = dynamic_cast<AliZDCMBCalib*> (entry->GetObject());
  if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}

//_____________________________________________________________________________
AliZDCTDCCalib* AliZDCReconstructor::GetTDCCalibData() const
{

  // Getting TDC object for ZDC 

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/TDCCalib");
  if(!entry) AliFatal("No calibration data loaded!");  
  entry->SetOwner(kFALSE);

  AliZDCTDCCalib *calibdata = dynamic_cast<AliZDCTDCCalib*> (entry->GetObject());
  if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}
