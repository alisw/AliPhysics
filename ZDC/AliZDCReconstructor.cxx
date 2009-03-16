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
// class for ZDC reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TF1.h>
#include <TMap.h>

#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliGRPObject.h"
#include "AliESDEvent.h"
#include "AliESDZDC.h"
#include "AliZDCDigit.h"
#include "AliZDCRawStream.h"
#include "AliZDCReco.h"
#include "AliZDCReconstructor.h"
#include "AliZDCPedestals.h"
#include "AliZDCCalib.h"
#include "AliZDCRecoParam.h"
#include "AliZDCRecoParampp.h"
#include "AliZDCRecoParamPbPb.h"


ClassImp(AliZDCReconstructor)
AliZDCRecoParam *AliZDCReconstructor::fRecoParam=0;  //reconstruction parameters

//_____________________________________________________________________________
AliZDCReconstructor:: AliZDCReconstructor() :
  fPedData(GetPedData()),
  fECalibData(GetECalibData()),
  fRecoMode(0),
  fBeamEnergy(0.),
  fPedSubMode(0)
{
  // **** Default constructor
  SetRecoMode();

}


//_____________________________________________________________________________
AliZDCReconstructor::~AliZDCReconstructor()
{
// destructor
   if(fRecoParam)  delete fRecoParam;
   if(fPedData)    delete fPedData;    
   if(fECalibData) delete fECalibData;
}

//____________________________________________________________________________
void AliZDCReconstructor::SetRecoMode()
{
  // Setting reconstruction mode

  // Initialization of the GRP entry 
  AliCDBEntry*  entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* grpData = 0x0;
  if(entry){
    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry
    if(m){
       m->Print();
       grpData = new AliGRPObject();
       grpData->ReadValuesFromMap(m);
    }
    else{
       grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0);
    }
    AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
  }
  if(!grpData) AliError("No GRP entry found in OCDB!");

  TString beamType = grpData->GetBeamType();
  if(beamType==AliGRPObject::GetInvalidString()){
    AliError("GRP/GRP/Data entry:  missing value for the beam energy !");
    AliError("\t ZDC does not reconstruct event 4 UNKNOWN beam type\n");
    return;
  }
  //
  if((beamType.CompareTo("p-p")) == 0){
    fRecoMode=0;
    fRecoParam = (AliZDCRecoParampp*) AliZDCRecoParampp::GetppRecoParam();
  }
  else if((beamType.CompareTo("A-A")) == 0){
    fRecoMode=1;
    fRecoParam = (AliZDCRecoParamPbPb*) AliZDCRecoParamPbPb::GetPbPbRecoParam();
  }
  
  fBeamEnergy = grpData->GetBeamEnergy();
  if(fBeamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    fBeamEnergy = 0.;
  }
  
  printf("\n ***** ZDC reconstruction initialized for %s @ %1.3f GeV\n\n",beamType.Data(), fBeamEnergy);
}

//_____________________________________________________________________________
void AliZDCReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  // *** Local ZDC reconstruction for digits
  // Works on the current event
    
  // Retrieving calibration data  
  // Parameters for mean value pedestal subtraction
  Float_t meanPed[48];    
  for(Int_t jj=0; jj<48; jj++) meanPed[jj] = fPedData->GetMeanPed(jj);
  // Parameters pedestal subtraction through correlation with out-of-time signals
  Float_t corrCoeff0[48], corrCoeff1[48];
  for(Int_t jj=0; jj<48; jj++){
     corrCoeff0[jj] =  fPedData->GetPedCorrCoeff0(jj);
     corrCoeff1[jj] =  fPedData->GetPedCorrCoeff1(jj);
  }

  // get digits
  AliZDCDigit digit;
  AliZDCDigit* pdigit = &digit;
  digitsTree->SetBranchAddress("ZDC", &pdigit);
  //printf("\n\t # of digits in tree: %d\n",(Int_t) digitsTree->GetEntries());

  // loop over digits
  Float_t tZN1Corr[10], tZP1Corr[10], tZN2Corr[10], tZP2Corr[10]; 
  Float_t dZEM1Corr[2], dZEM2Corr[2], PMRef1[2], PMRef2[2]; 
  for(Int_t i=0; i<10; i++){
     tZN1Corr[i] = tZP1Corr[i] = tZN2Corr[i] = tZP2Corr[i] = 0.;
     if(i<2) dZEM1Corr[i] = dZEM2Corr[i] = PMRef1[i] = PMRef2[i] = 0.;
  }  
  
  Int_t digNentries = digitsTree->GetEntries();
  int const kNch = 24;
  Float_t ootDigi[kNch];
  // -- Reading out-of-time signals (last kNch entries) for current event
  if(fPedSubMode==1){
    for(Int_t iDigit=kNch; iDigit<digNentries; iDigit++){
       ootDigi[iDigit] = digitsTree->GetEntry(iDigit);
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

   //printf("\n\t Digit #%d det %d quad %d", iDigit, det, quad);
      
   if(quad != 5){ // ZDC (not reference PTMs!)
    if(det == 1){ // *** ZNC
       tZN1Corr[quad] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       tZN1Corr[quad+5] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
       if(tZN1Corr[quad]<0.) tZN1Corr[quad] = 0.;
       if(tZN1Corr[quad+5]<0.) tZN1Corr[quad+5] = 0.;
       //printf("\t pedindex %d tZN1Corr[%d] = %1.0f tZN1Corr[%d] = %1.0f", 
       //	pedindex, quad, tZN1Corr[quad], quad+5, tZN1Corr[quad+5]);
    }
    else if(det == 2){ // *** ZP1
       tZP1Corr[quad] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       if(tZP1Corr[quad]<0.) tZP1Corr[quad] = 0.;
       tZP1Corr[quad+5] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
       if(tZP1Corr[quad+5]<0.) tZP1Corr[quad+5] = 0.;
       //printf("\t pedindex %d tZP1Corr[%d] = %1.0f tZP1Corr[%d] = %1.0f", 
       //	pedindex, quad, tZP1Corr[quad], quad+5, tZP1Corr[quad+5]);
    }
    else if(det == 3){
       if(quad == 1){	    // *** ZEM1  
         dZEM1Corr[0] += (Float_t) (digit.GetADCValue(0)-ped2SubHg); 
         if(dZEM1Corr[0]<0.) dZEM1Corr[0] = 0.;
         dZEM1Corr[1] += (Float_t) (digit.GetADCValue(1)-ped2SubLg); 
         if(dZEM1Corr[1]<0.) dZEM1Corr[1] = 0.;
         //printf("\t pedindex %d tZEM1Corr[%d] = %1.0f tZEM1Corr[%d] = %1.0f", 
         //	pedindex, quad, tZEM1Corr[quad], quad+1, tZEM1Corr[quad+1]);
       }
       else if(quad == 2){  // *** ZEM2
         dZEM2Corr[0] += (Float_t) (digit.GetADCValue(0)-ped2SubHg); 
         if(dZEM2Corr[0]<0.) dZEM2Corr[0] = 0.;
         dZEM2Corr[1] += (Float_t) (digit.GetADCValue(1)-ped2SubLg); 
         if(dZEM2Corr[1]<0.) dZEM2Corr[1] = 0.;
         //printf("\t pedindex %d tZEM2Corr[%d] = %1.0f tZEM2Corr[%d] = %1.0f", 
         //	pedindex, quad, tZEM2Corr[quad], quad+1, tZEM2Corr[quad+1]);
       }
    }
    else if(det == 4){  // *** ZN2
       tZN2Corr[quad] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       if(tZN2Corr[quad]<0.) tZN2Corr[quad] = 0.;
       tZN2Corr[quad+5] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
       if(tZN2Corr[quad+5]<0.) tZN2Corr[quad+5] = 0.;
       //printf("\t pedindex %d tZN2Corr[%d] = %1.0f tZN2Corr[%d] = %1.0f\n", 
       //	pedindex, quad, tZN2Corr[quad], quad+5, tZN2Corr[quad+5]);
    }
    else if(det == 5){  // *** ZP2 
       tZP2Corr[quad] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       if(tZP2Corr[quad]<0.) tZP2Corr[quad] = 0.;
       tZP2Corr[quad+5] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
       if(tZP2Corr[quad+5]<0.) tZP2Corr[quad+5] = 0.;
       //printf("\t pedindex %d tZP2Corr[%d] = %1.0f tZP2Corr[%d] = %1.0f\n", 
       //	pedindex, quad, tZP2Corr[quad], quad+5, tZP2Corr[quad+5]);
    }
   }
   else{ // Reference PMs
     if(det == 1){
       PMRef1[0] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       if(PMRef1[0]<0.) PMRef1[0] = 0.;
       PMRef1[1] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
       if(PMRef2[1]<0.) PMRef1[1] = 0.;
     }
     else if(det == 4){
       PMRef2[0] = (Float_t) (digit.GetADCValue(0)-ped2SubHg);
       if(PMRef2[0]<0.) PMRef2[0] = 0.;
       PMRef2[1] = (Float_t) (digit.GetADCValue(1)-ped2SubLg);
       if(PMRef2[1]<0.) PMRef2[1] = 0.;
     }
   }
  }

  // reconstruct the event
  if(fRecoMode==0)
    ReconstructEventpp(clustersTree, tZN1Corr, tZP1Corr, tZN2Corr, tZP2Corr, 
    	dZEM1Corr, dZEM2Corr, PMRef1, PMRef2);
  else if(fRecoMode==1)
    ReconstructEventPbPb(clustersTree, tZN1Corr, tZP1Corr, tZN2Corr, tZP2Corr, 
    	dZEM1Corr, dZEM2Corr, PMRef1, PMRef2);

}

//_____________________________________________________________________________
void AliZDCReconstructor::Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const
{
  // *** ZDC raw data reconstruction
  // Works on the current event
  
  // Retrieving calibration data  
  // Parameters for mean value pedestal subtraction
  Float_t meanPed[48];    
  for(Int_t jj=0; jj<48; jj++) meanPed[jj] = fPedData->GetMeanPed(jj);
  // Parameters pedestal subtraction through correlation with out-of-time signals
  Float_t corrCoeff0[48], corrCoeff1[48];
  for(Int_t jj=0; jj<48; jj++){
     corrCoeff0[jj] =  fPedData->GetPedCorrCoeff0(jj);
     corrCoeff1[jj] =  fPedData->GetPedCorrCoeff1(jj);
  }

  rawReader->Reset();
  
  // loop over raw data
  Float_t tZN1Corr[10], tZP1Corr[10], tZN2Corr[10], tZP2Corr[10]; 
  Float_t dZEM1Corr[2], dZEM2Corr[2], PMRef1[2], PMRef2[2]; 
  for(Int_t i=0; i<10; i++){
     tZN1Corr[i] = tZP1Corr[i] = tZN2Corr[i] = tZP2Corr[i] = 0.;
     if(i<2) dZEM1Corr[i] = dZEM2Corr[i] = PMRef1[i] = PMRef2[i] = 0.;
  }  
  //
  AliZDCRawStream rawData(rawReader);
  Int_t const kNch = 24;
  while(rawData.Next()) {
    if(rawData.IsADCDataWord()){
     Int_t det = rawData.GetSector(0);
     Int_t quad = rawData.GetSector(1);
     Int_t gain = rawData.GetADCGain();
     Int_t pedindex=0;
     //
     if(quad !=5){ // ZDCs (not reference PTMs)
      if(det == 1){    
        pedindex = quad;
        if(gain == 0) tZN1Corr[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else tZN1Corr[quad+5]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
      }
      else if(det == 2){ 
        pedindex = quad+5;
        if(gain == 0) tZP1Corr[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else tZP1Corr[quad+5]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
      }
      else if(det == 3){ 
        pedindex = quad+9;
        if(quad==1){	 
          if(gain == 0) dZEM1Corr[0] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else dZEM1Corr[1] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
        }
        else if(quad==2){ 
          if(gain == 0) dZEM2Corr[0] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else dZEM2Corr[1] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+kNch]); 
        }
      }
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
      //printf("\t AliZDCReconstructor - det %d quad %d res %d -> Ped[%d] = %1.0f\n", 
      //  det,quad,gain, pedindex, meanPed[pedindex]);
     }
     else{ // reference PM
       pedindex = (det-1)/3 + 22;
       if(det == 1){
         if(gain==0) PMRef1[0] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]);
	 else PMRef1[1] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]);
       }
       else if(det ==4){
         if(gain==0) PMRef2[0] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]);
	 else PMRef2[1] += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]);
       }
     }
    }//IsADCDataWord
  }
    
  // reconstruct the event
  if(fRecoMode==0)
    ReconstructEventpp(clustersTree, tZN1Corr, tZP1Corr, tZN2Corr, tZP2Corr, 
    	dZEM1Corr, dZEM2Corr, PMRef1, PMRef2);
  else if(fRecoMode==1)
    ReconstructEventPbPb(clustersTree, tZN1Corr, tZP1Corr, tZN2Corr, tZP2Corr, 
    	dZEM1Corr, dZEM2Corr, PMRef1, PMRef2);

}

//_____________________________________________________________________________
void AliZDCReconstructor::ReconstructEventpp(TTree *clustersTree, Float_t* ZN1ADCCorr, 
	Float_t* ZP1ADCCorr, Float_t* ZN2ADCCorr, Float_t* ZP2ADCCorr,
	Float_t* ZEM1ADCCorr, Float_t* ZEM2ADCCorr, Float_t* PMRef1, Float_t* PMRef2) const
{
  // ***** Reconstruct one event
  
  // *** RECONSTRUCTION FROM "REAL" DATA
  //
  // Retrieving calibration data
  // --- Equalization coefficients ---------------------------------------------
  Float_t equalCoeffZN1[5], equalCoeffZP1[5], equalCoeffZN2[5], equalCoeffZP2[5];
  for(Int_t ji=0; ji<5; ji++){
     equalCoeffZN1[ji] = fECalibData->GetZN1EqualCoeff(ji);
     equalCoeffZP1[ji] = fECalibData->GetZP1EqualCoeff(ji); 
     equalCoeffZN2[ji] = fECalibData->GetZN2EqualCoeff(ji); 
     equalCoeffZP2[ji] = fECalibData->GetZP2EqualCoeff(ji); 
  }
  // --- Energy calibration factors ------------------------------------
  Float_t calibEne[4];
  // **** Energy calibration coefficient set to 1 
  // **** (no trivial way to calibrate in p-p runs)
  //for(Int_t ij=0; ij<4; ij++) calibEne[ij] = fECalibData->GetEnCalib(ij);
  for(Int_t ij=0; ij<4; ij++) calibEne[ij] = 1.;
  
  // Equalization of detector responses
  Float_t equalTowZN1[10], equalTowZN2[10], equalTowZP1[10], equalTowZP2[10];
  for(Int_t gi=0; gi<5; gi++){
     equalTowZN1[gi] = ZN1ADCCorr[gi]*equalCoeffZN1[gi];
     equalTowZN1[gi+5] = ZN1ADCCorr[gi+5]*equalCoeffZN1[gi];
     equalTowZP1[gi] = ZP1ADCCorr[gi]*equalCoeffZP1[gi];
     equalTowZP1[gi+5] = ZP1ADCCorr[gi+5]*equalCoeffZP1[gi];
     equalTowZN2[gi] = ZN2ADCCorr[gi]*equalCoeffZN2[gi];
     equalTowZN2[gi+5] = ZN2ADCCorr[gi+5]*equalCoeffZN2[gi];
     equalTowZP2[gi] = ZP2ADCCorr[gi]*equalCoeffZP2[gi];
     equalTowZP2[gi+5] = ZP2ADCCorr[gi+5]*equalCoeffZP2[gi];
  }
  
  // Energy calibration of detector responses
  Float_t calibTowZN1[10], calibTowZN2[10], calibTowZP1[10], calibTowZP2[10];
  Float_t calibSumZN1[]={0,0}, calibSumZN2[]={0,0}, calibSumZP1[]={0,0}, calibSumZP2[]={0,0};
  for(Int_t gi=0; gi<10; gi++){
     calibTowZN1[gi] = equalTowZN1[gi]*calibEne[0];
     calibTowZP1[gi] = equalTowZP1[gi]*calibEne[1];
     calibTowZN2[gi] = equalTowZN2[gi]*calibEne[2];
     calibTowZP2[gi] = equalTowZP2[gi]*calibEne[3];
     //
     if(gi<5){
       calibSumZN1[0] += calibTowZN1[gi];
       calibSumZP1[0] += calibTowZP1[gi];
       calibSumZN2[0] += calibTowZN2[gi];
       calibSumZP2[0] += calibTowZP2[gi];
     }
     //
     else{
       calibSumZN1[1] += calibTowZN1[gi];
       calibSumZP1[1] += calibTowZP1[gi];
       calibSumZN2[1] += calibTowZN2[gi];
       calibSumZP2[1] += calibTowZP2[gi];
     }
  }
  
  //  ---      Number of detected spectator nucleons
  //  *** N.B. -> It works only in Pb-Pb!!!!!!!!!!!!
  //  Variables calculated to comply with ESD structure
  Int_t nDetSpecNLeft=0, nDetSpecPLeft=0, nDetSpecNRight=0, nDetSpecPRight=0;
  if(fBeamEnergy!=0){
   nDetSpecNLeft = (Int_t) (calibSumZN1[0]/fBeamEnergy);
   nDetSpecPLeft = (Int_t) (calibSumZP1[0]/fBeamEnergy);
   nDetSpecNRight = (Int_t) (calibSumZN2[0]/fBeamEnergy);
   nDetSpecPRight = (Int_t) (calibSumZP2[0]/fBeamEnergy);
  }
  else AliWarning(" ATTENTION -> fBeamEnergy = 0\n");
  /*printf("\n\t AliZDCReconstructor -> nDetSpecNLeft %d, nDetSpecPLeft %d,"
    " nDetSpecNRight %d, nDetSpecPRight %d\n",nDetSpecNLeft, nDetSpecPLeft, 
    nDetSpecNRight, nDetSpecPRight);*/

  //  ---      Number of generated spectator nucleons (from HIJING parameterization)
  Int_t nGenSpecNLeft=0, nGenSpecPLeft=0, nGenSpecLeft=0;
  Int_t nGenSpecNRight=0, nGenSpecPRight=0, nGenSpecRight=0;
  Int_t nPartTotLeft=0, nPartTotRight=0;
  Double_t impPar=0.;
  
  // create the output tree
  AliZDCReco reco(calibSumZN1, calibSumZP1, calibSumZN2, calibSumZP2, 
  		  calibTowZN1, calibTowZP1, calibTowZN2, calibTowZP2, 
		  ZEM1ADCCorr, ZEM2ADCCorr, PMRef1, PMRef2,
		  nDetSpecNLeft, nDetSpecPLeft, nDetSpecNRight, nDetSpecPRight, 
		  nGenSpecNLeft, nGenSpecPLeft, nGenSpecLeft, nGenSpecNRight, 
		  nGenSpecPRight, nGenSpecRight, nPartTotLeft, nPartTotRight, impPar);
		  
  AliZDCReco* preco = &reco;
  const Int_t kBufferSize = 4000;
  clustersTree->Branch("ZDC", "AliZDCReco", &preco, kBufferSize);

  // write the output tree
  clustersTree->Fill();
}

//_____________________________________________________________________________
void AliZDCReconstructor::ReconstructEventPbPb(TTree *clustersTree, Float_t* ZN1ADCCorr, 
	Float_t* ZP1ADCCorr, Float_t* ZN2ADCCorr, Float_t* ZP2ADCCorr,
	Float_t* ZEM1ADCCorr, Float_t* ZEM2ADCCorr, Float_t* PMRef1, Float_t* PMRef2) const
{
  // ***** Reconstruct one event
  
  // *** RECONSTRUCTION FROM "REAL" DATA
  //
  // Retrieving calibration data
  // --- Equalization coefficients ---------------------------------------------
  Float_t equalCoeffZN1[5], equalCoeffZP1[5], equalCoeffZN2[5], equalCoeffZP2[5];
  for(Int_t ji=0; ji<5; ji++){
     equalCoeffZN1[ji] = fECalibData->GetZN1EqualCoeff(ji);
     equalCoeffZP1[ji] = fECalibData->GetZP1EqualCoeff(ji); 
     equalCoeffZN2[ji] = fECalibData->GetZN2EqualCoeff(ji); 
     equalCoeffZP2[ji] = fECalibData->GetZP2EqualCoeff(ji); 
  }
  // --- Energy calibration factors ------------------------------------
  Float_t calibEne[4];
  for(Int_t ij=0; ij<4; ij++) calibEne[ij] = fECalibData->GetEnCalib(ij);
  
  // Equalization of detector responses
  Float_t equalTowZN1[10], equalTowZN2[10], equalTowZP1[10], equalTowZP2[10];
  for(Int_t gi=0; gi<5; gi++){
     equalTowZN1[gi] = ZN1ADCCorr[gi]*equalCoeffZN1[gi];
     equalTowZN1[gi+5] = ZN1ADCCorr[gi+5]*equalCoeffZN1[gi];
     equalTowZP1[gi] = ZP1ADCCorr[gi]*equalCoeffZP1[gi];
     equalTowZP1[gi+5] = ZP1ADCCorr[gi+5]*equalCoeffZP1[gi];
     equalTowZN2[gi] = ZN2ADCCorr[gi]*equalCoeffZN2[gi];
     equalTowZN2[gi+5] = ZN2ADCCorr[gi+5]*equalCoeffZN2[gi];
     equalTowZP2[gi] = ZP2ADCCorr[gi]*equalCoeffZP2[gi];
     equalTowZP2[gi+5] = ZP2ADCCorr[gi+5]*equalCoeffZP2[gi];
  }
  
  // Energy calibration of detector responses
  Float_t calibTowZN1[10], calibTowZN2[10], calibTowZP1[10], calibTowZP2[10];
  Float_t calibSumZN1[]={0,0}, calibSumZN2[]={0,0}, calibSumZP1[]={0,0}, calibSumZP2[]={0,0};
  for(Int_t gi=0; gi<10; gi++){
     calibTowZN1[gi] = equalTowZN1[gi]*calibEne[0];
     calibTowZP1[gi] = equalTowZP1[gi]*calibEne[1];
     calibTowZN2[gi] = equalTowZN2[gi]*calibEne[2];
     calibTowZP2[gi] = equalTowZP2[gi]*calibEne[3];
     //
     if(gi<5){
       calibSumZN1[0] += calibTowZN1[gi];
       calibSumZP1[0] += calibTowZP1[gi];
       calibSumZN2[0] += calibTowZN2[gi];
       calibSumZP2[0] += calibTowZP2[gi];
     }
     //
     else{
       calibSumZN1[1] += calibTowZN1[gi];
       calibSumZP1[1] += calibTowZP1[gi];
       calibSumZN2[1] += calibTowZN2[gi];
       calibSumZP2[1] += calibTowZP2[gi];
     }
  }

  //
  // --- Reconstruction parameters ------------------ 
  if(!fRecoParam)  fRecoParam = (AliZDCRecoParamPbPb*) AliZDCRecoParamPbPb::GetPbPbRecoParam();
  //
  Float_t endPointZEM = fRecoParam->GetZEMEndValue();
  Float_t cutFractionZEM = fRecoParam->GetZEMCutFraction();
  Float_t dZEMSup = fRecoParam->GetDZEMSup();
  Float_t dZEMInf = fRecoParam->GetDZEMInf();
  //
  Float_t cutValueZEM = endPointZEM*cutFractionZEM;
  Float_t supValueZEM = cutValueZEM+(endPointZEM*dZEMSup);
  Float_t infValueZEM = cutValueZEM-(endPointZEM*dZEMInf);
  //
  Float_t maxValEZN1  = fRecoParam->GetEZN1MaxValue();
  Float_t maxValEZP1  = fRecoParam->GetEZP1MaxValue();
  Float_t maxValEZDC1 = fRecoParam->GetEZDC1MaxValue();
  Float_t maxValEZN2  = fRecoParam->GetEZN2MaxValue();
  Float_t maxValEZP2  = fRecoParam->GetEZP2MaxValue();
  Float_t maxValEZDC2 = fRecoParam->GetEZDC2MaxValue();
  //
  //printf("\n\t AliZDCReconstructor -> ZEMEndPoint %1.0f, ZEMCutValue %1.0f,"
  //   " ZEMSupValue %1.0f, ZEMInfValue %1.0f\n",endPointZEM,cutValueZEM,supValueZEM,infValueZEM);
  
  //  ---      Number of detected spectator nucleons
  //  *** N.B. -> It works only in Pb-Pb
  Int_t nDetSpecNLeft=0, nDetSpecPLeft=0, nDetSpecNRight=0, nDetSpecPRight=0;
  if(fBeamEnergy!=0){
    nDetSpecNLeft = (Int_t) (calibSumZN1[0]/fBeamEnergy);
    nDetSpecPLeft = (Int_t) (calibSumZP1[0]/fBeamEnergy);
    nDetSpecNRight = (Int_t) (calibSumZN2[0]/fBeamEnergy);
    nDetSpecPRight = (Int_t) (calibSumZP2[0]/fBeamEnergy);
  }
  else AliWarning(" ATTENTION -> fBeamEnergy = 0\n");
  /*printf("\n\t AliZDCReconstructor -> nDetSpecNLeft %d, nDetSpecPLeft %d,"
    " nDetSpecNRight %d, nDetSpecPRight %d\n",nDetSpecNLeft, nDetSpecPLeft, 
    nDetSpecNRight, nDetSpecPRight);*/

  //  ---      Number of generated spectator nucleons (from HIJING parameterization)
  Int_t nGenSpecNLeft=0, nGenSpecPLeft=0, nGenSpecLeft=0;
  Int_t nGenSpecNRight=0, nGenSpecPRight=0, nGenSpecRight=0;
  Double_t impPar=0.;
  //
  Float_t corrADCZEMHG = ZEM1ADCCorr[0] + ZEM2ADCCorr[0];
  //
  if(corrADCZEMHG > supValueZEM){
    nGenSpecNLeft  = (Int_t) ((fRecoParam->GetfZNCen())->Eval(calibSumZN1[0]));
    nGenSpecPLeft  = (Int_t) ((fRecoParam->GetfZPCen())->Eval(calibSumZP1[0]));
    nGenSpecLeft   = (Int_t) ((fRecoParam->GetfZDCCen())->Eval(calibSumZN1[0]+calibSumZP1[0]));
    nGenSpecNRight = (Int_t) ((fRecoParam->GetfZNCen())->Eval(calibSumZN2[0]));
    nGenSpecPRight = (Int_t) ((fRecoParam->GetfZNCen())->Eval(calibSumZP2[0]));
    nGenSpecRight  = (Int_t) ((fRecoParam->GetfZNCen())->Eval(calibSumZN2[0]+calibSumZP2[0]));
    impPar  = (fRecoParam->GetfbCen())->Eval(calibSumZN1[0]+calibSumZP1[0]);
  }
  else if(corrADCZEMHG < infValueZEM){
    nGenSpecNLeft = (Int_t) ((fRecoParam->GetfZNPer())->Eval(calibSumZN1[0])); 
    nGenSpecPLeft = (Int_t) ((fRecoParam->GetfZPPer())->Eval(calibSumZP1[0]));
    nGenSpecLeft  = (Int_t) ((fRecoParam->GetfZDCPer())->Eval(calibSumZN1[0]+calibSumZP1[0]));
    impPar   = (fRecoParam->GetfbPer())->Eval(calibSumZN1[0]+calibSumZP1[0]);
  }
  else if(corrADCZEMHG >= infValueZEM && corrADCZEMHG <= supValueZEM){
    nGenSpecNLeft = (Int_t) ((fRecoParam->GetfZEMn())->Eval(corrADCZEMHG));
    nGenSpecPLeft = (Int_t) ((fRecoParam->GetfZEMp())->Eval(corrADCZEMHG));
    nGenSpecLeft  = (Int_t)((fRecoParam->GetfZEMsp())->Eval(corrADCZEMHG));
    impPar   =  (fRecoParam->GetfZEMb())->Eval(corrADCZEMHG);
  }
  // 
  if(calibSumZN1[0]/maxValEZN1>1.)  nGenSpecNLeft = (Int_t) ((fRecoParam->GetfZEMn())->Eval(corrADCZEMHG));
  if(calibSumZP1[0]/maxValEZP1>1.)  nGenSpecPLeft = (Int_t) ((fRecoParam->GetfZEMp())->Eval(corrADCZEMHG));
  if((calibSumZN1[0]+calibSumZP1[0]/maxValEZDC1)>1.){
     nGenSpecLeft = (Int_t)((fRecoParam->GetfZEMsp())->Eval(corrADCZEMHG));
     impPar = (fRecoParam->GetfZEMb())->Eval(corrADCZEMHG);
  }
  if(calibSumZN2[0]/maxValEZN2>1.)  nGenSpecNRight = (Int_t) ((fRecoParam->GetfZEMn())->Eval(corrADCZEMHG));
  if(calibSumZP2[0]/maxValEZP2>1.)  nGenSpecPRight = (Int_t) ((fRecoParam->GetfZEMp())->Eval(corrADCZEMHG));
  if((calibSumZN2[0]+calibSumZP2[0]/maxValEZDC2)>1.) nGenSpecRight = (Int_t)((fRecoParam->GetfZEMsp())->Eval(corrADCZEMHG));
  //
  if(nGenSpecNLeft>125)    nGenSpecNLeft=125;
  else if(nGenSpecNLeft<0) nGenSpecNLeft=0;
  if(nGenSpecPLeft>82)     nGenSpecPLeft=82;
  else if(nGenSpecPLeft<0) nGenSpecPLeft=0;
  if(nGenSpecLeft>207)     nGenSpecLeft=207;
  else if(nGenSpecLeft<0)  nGenSpecLeft=0;
  
  //  ---      Number of generated participants (from HIJING parameterization)
  Int_t nPart, nPartTotLeft, nPartTotRight;
  nPart = 207-nGenSpecNLeft-nGenSpecPLeft;
  nPartTotLeft = 207-nGenSpecLeft;
  nPartTotRight = 207-nGenSpecRight;
  if(nPart<0) nPart=0;
  if(nPartTotLeft<0) nPartTotLeft=0;
  if(nPartTotRight<0) nPartTotRight=0;
  //
  // *** DEBUG ***
  /*printf("\n\t AliZDCReconstructor -> calibSumZN1[0] %1.0f, calibSumZP1[0] %1.0f,"
      "  calibSumZN2[0] %1.0f, calibSumZP2[0] %1.0f, corrADCZEMHG %1.0f\n", 
      calibSumZN1[0],calibSumZP1[0],calibSumZN2[0],calibSumZP2[0],corrADCZEMHG);
  printf("\t AliZDCReconstructor -> nGenSpecNLeft %d, nGenSpecPLeft %d, nGenSpecLeft %d\n"
      "\t\t nGenSpecNRight %d, nGenSpecPRight %d, nGenSpecRight %d\n", 
      nGenSpecNLeft, nGenSpecPLeft, nGenSpecLeft, 
      nGenSpecNRight, nGenSpecPRight, nGenSpecRight);
  printf("\t AliZDCReconstructor ->  NpartL %d,  NpartR %d,  b %1.2f fm\n\n",nPartTotLeft, nPartTotRight, impPar);
  */
  
  // create the output tree
  AliZDCReco reco(calibSumZN1, calibSumZP1, calibSumZN2, calibSumZP2, 
  		  calibTowZN1, calibTowZP1, calibTowZN2, calibTowZP2, 
		  ZEM1ADCCorr, ZEM2ADCCorr, PMRef1, PMRef2,
		  nDetSpecNLeft, nDetSpecPLeft, nDetSpecNRight, nDetSpecPRight, 
		  nGenSpecNLeft, nGenSpecPLeft, nGenSpecLeft, nGenSpecNRight, 
		  nGenSpecPRight, nGenSpecRight, nPartTotLeft, nPartTotRight, impPar);
		  
  AliZDCReco* preco = &reco;
  const Int_t kBufferSize = 4000;
  clustersTree->Branch("ZDC", "AliZDCReco", &preco, kBufferSize);

  // write the output tree
  clustersTree->Fill();
}

//_____________________________________________________________________________
void AliZDCReconstructor::FillZDCintoESD(TTree *clustersTree, AliESDEvent* esd) const
{
  // fill energies and number of participants to the ESD

  AliZDCReco reco;
  AliZDCReco* preco = &reco;
  clustersTree->SetBranchAddress("ZDC", &preco);

  clustersTree->GetEntry(0);
  //
  AliESDZDC * esdzdc = esd->GetESDZDC();
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
  esdzdc->SetZN1TowerEnergy(tZN1Ene);
  esdzdc->SetZN2TowerEnergy(tZN2Ene);
  esdzdc->SetZP1TowerEnergy(tZP1Ene);
  esdzdc->SetZP2TowerEnergy(tZP2Ene);
  esdzdc->SetZN1TowerEnergyLR(tZN1EneLR);
  esdzdc->SetZN2TowerEnergyLR(tZN2EneLR);
  esdzdc->SetZP1TowerEnergyLR(tZP1EneLR);
  esdzdc->SetZP2TowerEnergyLR(tZP2EneLR);
  // 
  esd->SetZDC(reco.GetZN1HREnergy(), reco.GetZP1HREnergy(), reco.GetZEM1HRsignal(), 
  	      reco.GetZEM2HRsignal(), reco.GetZN2HREnergy(), reco.GetZP2HREnergy(), 
	      reco.GetNPartLeft(), reco.GetNPartRight());
  //
  
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
AliZDCPedestals* AliZDCReconstructor::GetPedData() const
{

  // Getting pedestal calibration object for ZDC set

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/Pedestals");
  if(!entry) AliFatal("No calibration data loaded!");  

  AliZDCPedestals *calibdata = dynamic_cast<AliZDCPedestals*>  (entry->GetObject());
  if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}

//_____________________________________________________________________________
AliZDCCalib* AliZDCReconstructor::GetECalibData() const
{

  // Getting energy and equalization calibration object for ZDC set

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/EMDCalib");
  if(!entry) AliFatal("No calibration data loaded!");  

  AliZDCCalib *calibdata = dynamic_cast<AliZDCCalib*>  (entry->GetObject());
  if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}

