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

#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliESDEvent.h"
#include "AliZDCDigit.h"
#include "AliZDCRawStream.h"
#include "AliZDCReco.h"
#include "AliZDCReconstructor.h"
#include "AliZDCCalibData.h"


ClassImp(AliZDCReconstructor)


//_____________________________________________________________________________
AliZDCReconstructor:: AliZDCReconstructor() :

  fZNCen(new TF1("fZNCen", 
	"(-2.287920+sqrt(2.287920*2.287920-4*(-0.007629)*(11.921710-x)))/(2*(-0.007629))",0.,164.)),
  fZNPer(new TF1("fZNPer",
      "(-37.812280-sqrt(37.812280*37.812280-4*(-0.190932)*(-1709.249672-x)))/(2*(-0.190932))",0.,164.)),
  fZPCen(new TF1("fZPCen",
       "(-1.321353+sqrt(1.321353*1.321353-4*(-0.007283)*(3.550697-x)))/(2*(-0.007283))",0.,60.)),
  fZPPer(new TF1("fZPPer",
      "(-42.643308-sqrt(42.643308*42.643308-4*(-0.310786)*(-1402.945615-x)))/(2*(-0.310786))",0.,60.)),
  fZDCCen(new TF1("fZDCCen",
      "(-1.934991+sqrt(1.934991*1.934991-4*(-0.004080)*(15.111124-x)))/(2*(-0.004080))",0.,225.)),
  fZDCPer(new TF1("fZDCPer",
      "(-34.380639-sqrt(34.380639*34.380639-4*(-0.104251)*(-2612.189017-x)))/(2*(-0.104251))",0.,225.)),
  fbCen(new TF1("fbCen","-0.056923+0.079703*x-0.0004301*x*x+0.000001366*x*x*x",0.,220.)),
  fbPer(new TF1("fbPer","17.943998-0.046846*x+0.000074*x*x",0.,220.)),
  fZEMn(new TF1("fZEMn","126.2-0.05399*x+0.000005679*x*x",0.,4000.)),
  fZEMp(new TF1("fZEMp","82.49-0.03611*x+0.00000385*x*x",0.,4000.)),
  fZEMsp(new TF1("fZEMsp","208.7-0.09006*x+0.000009526*x*x",0.,4000.)),
  fZEMb(new TF1("fZEMb",
	"16.06-0.01633*x+1.44e-5*x*x-6.778e-9*x*x*x+1.438e-12*x*x*x*x-1.112e-16*x*x*x*x*x",0.,4000.)),
  //
  fCalibData(GetCalibData())

{
  // **** Default constructor

}


//_____________________________________________________________________________
AliZDCReconstructor::~AliZDCReconstructor()
{
// destructor

  delete fZNCen;
  delete fZNPer;
  delete fZPCen;
  delete fZPPer;
  delete fZDCCen;
  delete fZDCPer;
  delete fbCen;
  delete fbPer;
  delete fZEMn;
  delete fZEMp;
  delete fZEMsp;
  delete fZEMb;

}


//_____________________________________________________________________________
void AliZDCReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  // *** Local ZDC reconstruction for digits
  // Works on the current event
    
  // Retrieving calibration data  
  Float_t meanPed[47];
  for(Int_t jj=0; jj<47; jj++) meanPed[jj] = fCalibData->GetMeanPed(jj);

  // get digits
  AliZDCDigit digit;
  AliZDCDigit* pdigit = &digit;
  digitsTree->SetBranchAddress("ZDC", &pdigit);

  // loop over digits
  Float_t tZN1CorrHG[]={0.,0.,0.,0.,0.}, tZP1CorrHG[]={0.,0.,0.,0.,0.}; 
  Float_t dZEMCorrHG=0.; 
  Float_t tZN2CorrHG[]={0.,0.,0.,0.,0.}, tZP2CorrHG[]={0.,0.,0.,0.,0.};
  Float_t tZN1CorrLG[]={0.,0.,0.,0.,0.}, tZP1CorrLG[]={0.,0.,0.,0.,0.};
  Float_t dZEMCorrLG=0.; 
  Float_t tZN2CorrLG[]={0.,0.,0.,0.,0.}, tZP2CorrLG[]={0.,0.,0.,0.,0.};
  
  //printf("\n\t # of digits in tree: %d\n",(Int_t) digitsTree->GetEntries());
  for (Int_t iDigit = 0; iDigit < (digitsTree->GetEntries()/2); iDigit++) {
    digitsTree->GetEntry(iDigit);
    if (!pdigit) continue;
    //pdigit->Print("");
    //  
    Int_t det = digit.GetSector(0);
    Int_t quad = digit.GetSector(1);
    Int_t pedindex = -1;
    //printf("\n\t #%d det %d quad %d", iDigit, det, quad);
    //
    if(det == 1){ // *** ZN1
       pedindex = quad;
       tZN1CorrHG[quad] = (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]);
       if(tZN1CorrHG[quad]<0.) tZN1CorrHG[quad] = 0.;
       tZN1CorrLG[quad] = (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+5]);
       if(tZN1CorrLG[quad]<0.) tZN1CorrLG[quad] = 0.;
       //printf("\t pedindex %d tZN1CorrHG[%d] = %1.0f tZN1CorrLG[%d] = %1.0f", 
       //	pedindex, quad, tZN1CorrHG[quad], quad, tZN1CorrLG[quad]);
    }
    else if(det == 2){ // *** ZP1
       pedindex = quad+10;
       tZP1CorrHG[quad] = (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]);
       if(tZP1CorrLG[quad]<0.) tZP1CorrLG[quad] = 0.;
       tZP1CorrLG[quad] = (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+5]);
       if(tZP1CorrHG[quad]<0.) tZP1CorrHG[quad] = 0.;
       //printf("\t pedindex %d tZP1CorrHG[%d] = %1.0f tZP1CorrLG[%d] = %1.0f", 
       //	pedindex, quad, tZP1CorrHG[quad], quad, tZP1CorrLG[quad]);
    }
    else if(det == 3){
       if(quad == 1){	    // *** ZEM1  
    	 pedindex = quad+19;
         dZEMCorrHG += (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]); 
         if(dZEMCorrHG<0.) dZEMCorrHG = 0.;
         dZEMCorrLG += (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+2]); 
         if(dZEMCorrLG<0.) dZEMCorrLG = 0.;
         //printf("\t pedindex %d ADC(0) = %d ped = %1.0f ADCCorr = %1.0f\n", 
	 //	pedindex, digit.GetADCValue(0), meanPed[pedindex], dZEMCorrHG);
         //printf("\t pedindex %d ADC(1) = %d ped = %1.0f ADCCorr = %1.0f\n", 
	 //	pedindex+2, digit.GetADCValue(1), meanPed[pedindex+2], dZEMCorrLG);
         ////printf("\t pedindex %d dZEMCorrHG = %1.0f dZEMCorrLG = %1.0f\n", pedindex, dZEMCorrHG, dZEMCorrLG);
       }
       else if(quad == 2){  // *** ZEM1
    	 pedindex = quad+19;
         dZEMCorrHG += (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]); 
         if(dZEMCorrHG<0.) dZEMCorrHG = 0.;
         dZEMCorrLG += (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+2]); 
         if(dZEMCorrLG<0.) dZEMCorrLG = 0.;
         //printf("\t pedindex %d ADC(0) = %d ped = %1.0f ADCCorr = %1.0f\n", 
	 //	pedindex, digit.GetADCValue(0), meanPed[pedindex], dZEMCorrHG);
         //printf("\t pedindex %d ADC(1) = %d ped = %1.0f ADCCorr = %1.0f\n", 
	 //	pedindex+2, digit.GetADCValue(1),meanPed[pedindex+2], dZEMCorrLG);
         ////printf("\t pedindex %d dZEMCorrHG = %1.0f dZEMCorrLG = %1.0f\n", pedindex, dZEMCorrHG, dZEMCorrLG);
       }
    }
    else if(det == 4){  // *** ZN2
       pedindex = quad+24;
       tZN2CorrHG[quad] = (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]);
       if(tZN2CorrHG[quad]<0.) tZN2CorrHG[quad] = 0.;
       tZN2CorrLG[quad] = (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+5]);
       if(tZN2CorrLG[quad]<0.) tZN2CorrLG[quad] = 0.;
       //printf("\t pedindex %d tZN2CorrHG[%d] = %1.0f tZN2CorrLG[%d] = %1.0f\n", 
       //	pedindex, quad, tZN2CorrHG[quad], quad, tZN2CorrLG[quad]);
    }
    else if(det == 5){  // *** ZP2 
       pedindex = quad+34;
       tZP2CorrHG[quad] = (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]);
       if(tZP2CorrHG[quad]<0.) tZP2CorrHG[quad] = 0.;
       tZP2CorrLG[quad] = (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+5]);
       if(tZP2CorrLG[quad]<0.) tZP2CorrLG[quad] = 0.;
       //printf("\t pedindex %d tZP2CorrHG[%d] = %1.0f tZP2CorrLG[%d] = %1.0f\n", 
       //	pedindex, quad, tZP2CorrHG[quad], quad, tZP2CorrLG[quad]);
    }
  }

  // reconstruct the event
    ReconstructEvent(clustersTree, tZN1CorrHG, tZP1CorrHG, tZN2CorrHG, 
    	tZP2CorrHG, tZN1CorrLG, tZP1CorrLG, tZN2CorrLG, 
    	tZP2CorrLG, dZEMCorrHG);

}

//_____________________________________________________________________________
void AliZDCReconstructor::Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const
{
  // *** ZDC raw data reconstruction
  // Works on the current event
  
  // Retrieving calibration data  
  Float_t meanPed[47];
  for(Int_t jj=0; jj<47; jj++) meanPed[jj] = fCalibData->GetMeanPed(jj);

  rawReader->Reset();

  // loop over raw data rawDatas
  Float_t tZN1CorrHG[]={0.,0.,0.,0.,0.}, tZP1CorrHG[]={0.,0.,0.,0.,0.};
  Float_t dZEMCorrHG=0.;
  Float_t tZN2CorrHG[]={0.,0.,0.,0.,0.}, tZP2CorrHG[]={0.,0.,0.,0.,0.};
  Float_t tZN1CorrLG[]={0.,0.,0.,0.,0.}, tZP1CorrLG[]={0.,0.,0.,0.,0.};
  Float_t dZEMCorrLG=0.; 
  Float_t tZN2CorrLG[]={0.,0.,0.,0.,0.}, tZP2CorrLG[]={0.,0.,0.,0.,0.};
  //
  AliZDCRawStream rawData(rawReader);
  while (rawData.Next()) {
    if(rawData.IsADCDataWord()){
      Int_t det = rawData.GetSector(0);
      Int_t quad = rawData.GetSector(1);
      Int_t gain = rawData.GetADCGain();
      Int_t pedindex;
      //
      if(det == 1){    
        pedindex = quad;
        if(gain == 0) tZN1CorrHG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else tZN1CorrLG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+5]); 
      }
      else if(det == 2){ 
        pedindex = quad+10;
        if(gain == 0) tZP1CorrHG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else tZP1CorrLG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+5]); 
      }
      else if(det == 3){ 
        if(quad==1){	 
          pedindex = quad+20;
          if(gain == 0) dZEMCorrHG += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else dZEMCorrLG += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+2]); 
        }
        else if(quad==2){ 
          pedindex = rawData.GetSector(1)+21;
          if(gain == 0) dZEMCorrHG += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else dZEMCorrLG += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+2]); 
        }
      }
      else if(det == 4){       
        pedindex = rawData.GetSector(1)+24;
        if(gain == 0) tZN2CorrHG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else tZN2CorrLG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+2]); 
      }
      else if(det == 5){
        pedindex = rawData.GetSector(1)+34;
        if(gain == 0) tZP2CorrHG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else tZP2CorrLG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+5]); 
      }
    }
  }
    
  // reconstruct the event
    ReconstructEvent(clustersTree, tZN1CorrHG, tZP1CorrHG, tZN2CorrHG, 
    	tZP2CorrHG, tZN1CorrLG, tZP1CorrLG, tZN2CorrLG, 
    	tZP2CorrLG, dZEMCorrHG);

}

//_____________________________________________________________________________
void AliZDCReconstructor::ReconstructEvent(TTree *clustersTree, 
		Float_t* ZN1ADCCorrHG, Float_t* ZP1ADCCorrHG, 
		Float_t* ZN2ADCCorrHG, Float_t* ZP2ADCCorrHG, 
		Float_t* ZN1ADCCorrLG, Float_t* ZP1ADCCorrLG, 
		Float_t* ZN2ADCCorrLG, Float_t* ZP2ADCCorrLG, 
		Float_t corrADCZEMHG) const
{
  // ***** Reconstruct one event
  
  // *** RECONSTRUCTION FROM SIMULATED DATA
  // It passes trhough the no. of phe which is known from simulations
  //  ---      ADCchannel -> photoelectrons
  // NB-> PM gain = 10^(5), ADC resolution = 6.4*10^(-7)
  // Move to V965 (E.S.,15/09/04) NB-> PM gain = 10^(5), ADC resolution = 8*10^(-7)
  //Float_t zn1phe, zp1phe, zemphe, zn2phe, zp2phe, convFactor = 0.08;
  //zn1phe  = ZN1Corr/convFactor;
  //zp1phe  = ZP1Corr/convFactor;
  //zemphe = ZEMCorr/convFactor;
  //zn2phe  = ZN2Corr/convFactor;
  //zp2phe  = ZP2Corr/convFactor;
  ////if AliDebug(1,Form("\n    znphe = %f, zpphe = %f, zemphe = %f\n",znphe, zpphe, zemphe);
  //
  ////  ---      Energy calibration
  //// Conversion factors for hadronic ZDCs goes from phe yield to TRUE 
  //// incident energy (conversion from GeV to TeV is included); while for EM 
  //// calos conversion is from light yield to detected energy calculated by
  //// GEANT NB -> ZN and ZP conversion factors are constant since incident
  //// spectators have all the same energy, ZEM energy is obtained through a
  //// fit over the whole range of incident particle energies 
  //// (obtained with full HIJING simulations) 
  //Float_t zn1energy, zp1energy, zemenergy, zdc1energy, zn2energy, zp2energy, zdc2energy;
  //Float_t zn1phexTeV=329., zp1phexTeV=369., zn2phexTeV=329., zp2phexTeV=369.;
  //zn1energy  = zn1phe/zn1phexTeV;
  //zp1energy  = zp1phe/zp1phexTeV;
  //zdc1energy = zn1energy+zp1energy;
  //zn2energy  = zn2phe/zn2phexTeV;
  //zp2energy  = zp2phe/zp2phexTeV;
  //zdc2energy = zn2energy+zp2energy;
  //zemenergy = -4.81+0.3238*zemphe;
  //if(zemenergy<0) zemenergy=0;
  ////  if AliDebug(1,Form("    znenergy = %f TeV, zpenergy = %f TeV, zdcenergy = %f GeV, "
  ////			   "\n		zemenergy = %f TeV\n", znenergy, zpenergy, 
  ////			   zdcenergy, zemenergy);
  ////  if(zdcenergy==0)
  ////    if AliDebug(1,Form("\n\n	###	ATTENZIONE!!! -> ev# %d: znenergy = %f TeV, zpenergy = %f TeV, zdcenergy = %f GeV, "
  ////			     " zemenergy = %f TeV\n\n", fMerger->EvNum(), znenergy, zpenergy, zdcenergy, zemenergy); 
  
  //
  // *** RECONSTRUCTION FROM "REAL" DATA
  //
  // Retrieving calibration data
  Float_t equalCoeffZN1[5], equalCoeffZP1[5], equalCoeffZN2[5], equalCoeffZP2[5];
  for(Int_t ji=0; ji<5; ji++){
     equalCoeffZN1[ji] = fCalibData->GetZN1EqualCoeff(ji);
     equalCoeffZP1[ji] = fCalibData->GetZP1EqualCoeff(ji); 
     equalCoeffZN2[ji] = fCalibData->GetZN2EqualCoeff(ji); 
     equalCoeffZP2[ji] = fCalibData->GetZP2EqualCoeff(ji); 
  }
  //
  Float_t calibEne[4];
  for(Int_t ij=0; ij<4; ij++) calibEne[ij] = fCalibData->GetEnCalib(ij);
  //
  Float_t endPointZEM = fCalibData->GetZEMEndValue();
  Float_t cutFractionZEM = fCalibData->GetZEMCutFraction();
  Float_t dZEMSup = fCalibData->GetDZEMSup();
  Float_t dZEMInf = fCalibData->GetDZEMInf();
  //
  Float_t cutValueZEM = endPointZEM*cutFractionZEM;
  Float_t supValueZEM = cutValueZEM+(endPointZEM*dZEMSup);
  Float_t infValueZEM = cutValueZEM-(endPointZEM*dZEMInf);
  //
  Float_t maxValEZN1 = fCalibData->GetEZN1MaxValue();
  Float_t maxValEZP1 = fCalibData->GetEZP1MaxValue();
  Float_t maxValEZDC1 = fCalibData->GetEZDC1MaxValue();
  Float_t maxValEZN2 = fCalibData->GetEZN1MaxValue();
  Float_t maxValEZP2 = fCalibData->GetEZP1MaxValue();
  Float_t maxValEZDC2 = fCalibData->GetEZDC1MaxValue();
  //
  //printf("\n\t AliZDCReconstructor -> ZEMEndPoint %1.0f, ZEMCutValue %1.0f,"
  //   " ZEMSupValue %1.0f, ZEMInfValue %1.0f\n",endPointZEM,cutValueZEM,supValueZEM,infValueZEM);
  
  // Equalization of detector responses
  Float_t equalTowZN1HG[5], equalTowZN2HG[5], equalTowZP1HG[5], equalTowZP2HG[5];
  Float_t equalTowZN1LG[5], equalTowZN2LG[5], equalTowZP1LG[5], equalTowZP2LG[5];
  for(Int_t gi=0; gi<5; gi++){
     equalTowZN1HG[gi] = ZN1ADCCorrHG[gi]*equalCoeffZN1[gi];
     equalTowZP1HG[gi] = ZP1ADCCorrHG[gi]*equalCoeffZP1[gi];
     equalTowZN2HG[gi] = ZN2ADCCorrHG[gi]*equalCoeffZN2[gi];
     equalTowZP2HG[gi] = ZP2ADCCorrHG[gi]*equalCoeffZP2[gi];
     //
     equalTowZN1LG[gi] = ZN1ADCCorrLG[gi]*equalCoeffZN1[gi];
     equalTowZP1LG[gi] = ZP1ADCCorrLG[gi]*equalCoeffZP1[gi];
     equalTowZN2LG[gi] = ZN2ADCCorrLG[gi]*equalCoeffZN2[gi];
     equalTowZP2LG[gi] = ZP2ADCCorrLG[gi]*equalCoeffZP2[gi];
  }
  
  // Energy calibration of detector responses
  Float_t calibTowZN1HG[5], calibTowZN2HG[5], calibTowZP1HG[5], calibTowZP2HG[5];
  Float_t calibSumZN1HG=0., calibSumZN2HG=0., calibSumZP1HG=0., calibSumZP2HG=0.;
  Float_t calibTowZN1LG[5], calibTowZN2LG[5], calibTowZP1LG[5], calibTowZP2LG[5];
  Float_t calibSumZN1LG=0., calibSumZN2LG=0., calibSumZ12LG=0., calibSumZP2LG=0.;
  for(Int_t gi=0; gi<5; gi++){
     calibTowZN1HG[gi] = equalTowZN1HG[gi]*calibEne[0];
     calibTowZP1HG[gi] = equalTowZP1HG[gi]*calibEne[1];
     calibTowZN2HG[gi] = equalTowZN2HG[gi]*calibEne[2];
     calibTowZP2HG[gi] = equalTowZP2HG[gi]*calibEne[3];
     calibSumZN1HG += calibTowZN1HG[gi];
     calibSumZP1HG += calibTowZP1HG[gi];
     calibSumZN2HG += calibTowZN2HG[gi];
     calibSumZP2HG += calibTowZP2HG[gi];
     //
     calibTowZN1LG[gi] = equalTowZN1LG[gi]*calibEne[0];
     calibTowZP1LG[gi] = equalTowZP1LG[gi]*calibEne[1];
     calibTowZN2LG[gi] = equalTowZN2LG[gi]*calibEne[2];
     calibTowZP2LG[gi] = equalTowZP2LG[gi]*calibEne[3];
     calibSumZN1LG += calibTowZN1LG[gi];
     calibSumZ12LG += calibTowZP1LG[gi];
     calibSumZN2LG += calibTowZN2LG[gi];
     calibSumZP2LG += calibTowZP2LG[gi];
  }
  
  //  ---      Number of detected spectator nucleons
  //  *** N.B. -> It works only in Pb-Pb
  Int_t nDetSpecNLeft, nDetSpecPLeft, nDetSpecNRight, nDetSpecPRight;
  nDetSpecNLeft = (Int_t) (calibSumZN1HG/2.760);
  nDetSpecPLeft = (Int_t) (calibSumZP1HG/2.760);
  nDetSpecNRight = (Int_t) (calibSumZN2HG/2.760);
  nDetSpecPRight = (Int_t) (calibSumZP2HG/2.760);
  /*printf("\n\t AliZDCReconstructor -> nDetSpecNLeft %d, nDetSpecPLeft %d,"
    " nDetSpecNRight %d, nDetSpecPRight %d\n",nDetSpecNLeft, nDetSpecPLeft, 
    nDetSpecNRight, nDetSpecPRight);*/

  //  ---      Number of generated spectator nucleons (from HIJING parameterization)
  Int_t nGenSpecNLeft=0, nGenSpecPLeft=0, nGenSpecLeft=0;
  Int_t nGenSpecNRight=0, nGenSpecPRight=0, nGenSpecRight=0;
  Double_t impPar=0.;
  //
  // *** RECONSTRUCTION FROM SIMULATED DATA
  // Cut value for Ezem (GeV)
  // ### Results from production  -> 0<b<18 fm (Apr 2002)
  /*Float_t eZEMCut = 420.;
  Float_t deltaEZEMSup = 690.; 
  Float_t deltaEZEMInf = 270.; 
  if(zemenergy > (eZEMCut+deltaEZEMSup)){
    nGenSpecNLeft  = (Int_t) (fZNCen->Eval(ZN1CalibSum));
    nGenSpecPLeft  = (Int_t) (fZPCen->Eval(ZP1CalibSum));
    nGenSpecLeft   = (Int_t) (fZDCCen->Eval(ZN1CalibSum+ZP1CalibSum));
    nGenSpecNRight = (Int_t) (fZNCen->Eval(ZN2CalibSum));
    nGenSpecPRight = (Int_t) (fZNCen->Eval(ZP2CalibSum));
    nGenSpecRight  = (Int_t) (fZNCen->Eval(ZN2CalibSum+ZP2CalibSum));
    impPar  = fbCen->Eval(ZN1CalibSum+ZP1CalibSum);
  }
  else if(zemenergy < (eZEMCut-deltaEZEMInf)){
    nGenSpecNLeft = (Int_t) (fZNPer->Eval(ZN1CalibSum)); 
    nGenSpecPLeft = (Int_t) (fZPPer->Eval(ZP1CalibSum));
    nGenSpecLeft  = (Int_t) (fZDCPer->Eval(ZN1CalibSum+ZP1CalibSum));
    impPar   = fbPer->Eval(ZN1CalibSum+ZP1CalibSum);
  }
  else if(zemenergy >= (eZEMCut-deltaEZEMInf) && zemenergy <= (eZEMCut+deltaEZEMSup)){
    nGenSpecNLeft = (Int_t) (fZEMn->Eval(zemenergy));
    nGenSpecPLeft = (Int_t) (fZEMp->Eval(zemenergy));
    nGenSpecLeft  = (Int_t)(fZEMsp->Eval(zemenergy));
    impPar   =  fZEMb->Eval(zemenergy);
  }
  // ### Results from production  -> 0<b<18 fm (Apr 2002)
  if(ZN1CalibSum>162.)  nGenSpecNLeft = (Int_t) (fZEMn->Eval(zemenergy));
  if(ZP1CalibSum>59.75)  nGenSpecPLeft = (Int_t) (fZEMp->Eval(zemenergy));
  if(ZN1CalibSum+ZP1CalibSum>221.5) nGenSpecLeft  = (Int_t)(fZEMsp->Eval(zemenergy));
  if(ZN1CalibSum+ZP1CalibSum>220.)  impPar    =  fZEMb->Eval(zemenergy);
  */
  //
  //
  // *** RECONSTRUCTION FROM REAL DATA
  //
  if(corrADCZEMHG > supValueZEM){
    nGenSpecNLeft  = (Int_t) (fZNCen->Eval(calibSumZN1HG));
    nGenSpecPLeft  = (Int_t) (fZPCen->Eval(calibSumZP1HG));
    nGenSpecLeft   = (Int_t) (fZDCCen->Eval(calibSumZN1HG+calibSumZP1HG));
    nGenSpecNRight = (Int_t) (fZNCen->Eval(calibSumZN2HG));
    nGenSpecPRight = (Int_t) (fZNCen->Eval(calibSumZP2HG));
    nGenSpecRight  = (Int_t) (fZNCen->Eval(calibSumZN2HG+calibSumZP2HG));
    impPar  = fbCen->Eval(calibSumZN1HG+calibSumZP1HG);
  }
  else if(corrADCZEMHG < infValueZEM){
    nGenSpecNLeft = (Int_t) (fZNPer->Eval(calibSumZN1HG)); 
    nGenSpecPLeft = (Int_t) (fZPPer->Eval(calibSumZP1HG));
    nGenSpecLeft  = (Int_t) (fZDCPer->Eval(calibSumZN1HG+calibSumZP1HG));
    impPar   = fbPer->Eval(calibSumZN1HG+calibSumZP1HG);
  }
  else if(corrADCZEMHG >= infValueZEM && corrADCZEMHG <= supValueZEM){
    nGenSpecNLeft = (Int_t) (fZEMn->Eval(corrADCZEMHG));
    nGenSpecPLeft = (Int_t) (fZEMp->Eval(corrADCZEMHG));
    nGenSpecLeft  = (Int_t)(fZEMsp->Eval(corrADCZEMHG));
    impPar   =  fZEMb->Eval(corrADCZEMHG);
  }
  // 
  if(calibSumZN1HG/maxValEZN1>1.)  nGenSpecNLeft = (Int_t) (fZEMn->Eval(corrADCZEMHG));
  if(calibSumZP1HG/maxValEZP1>1.)  nGenSpecPLeft = (Int_t) (fZEMp->Eval(corrADCZEMHG));
  if((calibSumZN1HG+calibSumZP1HG/maxValEZDC1)>1.){
     nGenSpecLeft = (Int_t)(fZEMsp->Eval(corrADCZEMHG));
     impPar = fZEMb->Eval(corrADCZEMHG);
  }
  if(calibSumZN2HG/maxValEZN2>1.)  nGenSpecNRight = (Int_t) (fZEMn->Eval(corrADCZEMHG));
  if(calibSumZP2HG/maxValEZP2>1.)  nGenSpecPRight = (Int_t) (fZEMp->Eval(corrADCZEMHG));
  if((calibSumZN2HG+calibSumZP2HG/maxValEZDC2)>1.) nGenSpecRight = (Int_t)(fZEMsp->Eval(corrADCZEMHG));
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
  //
  /*printf("\n\t AliZDCReconstructor -> nGenSpecNLeft %d, nGenSpecPLeft %d, nGenSpecLeft %d,"
      " nGenSpecNRight %d, nGenSpecPRight %d, nGenSpecRight %d\n", 
      nGenSpecNLeft, nGenSpecPLeft, nGenSpecLeft, nGenSpecNRight, 
      nGenSpecPRight, nGenSpecRight);*/

  // create the output tree
  AliZDCReco reco(calibSumZN1HG, calibSumZP1HG, calibSumZN2HG, calibSumZP2HG, 
  		  calibTowZN1LG, calibTowZN2LG, calibTowZP1LG, calibTowZP2LG, 
		  corrADCZEMHG, 
		  nDetSpecNLeft, nDetSpecPLeft, nDetSpecNRight, nDetSpecPRight, 
		  nGenSpecNLeft, nGenSpecPLeft, nGenSpecLeft, nGenSpecNRight, 
		  nGenSpecPRight, nGenSpecRight,
		  nPartTotLeft, nPartTotRight, impPar);
		  
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
  esd->SetZDC(reco.GetZN1Energy(), reco.GetZP1Energy(), reco.GetZEMsignal(),
	      reco.GetZN2Energy(), reco.GetZP2Energy(), 
	      reco.GetNPartLeft());
  //
  /*Double_t tZN1Ene[4], tZN2Ene[4];
  for(Int_t i=0; i<4; i++){
     tZN1Ene[i] = reco.GetZN1EnTow(i);
     tZN2Ene[i] = reco.GetZN2EnTow(i);
  }*/
  //esd->SetZN1TowerEnergy(tZN1Ene);
  //esd->SetZN2TowerEnergy(tZN2Ene);
  
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
AliZDCCalibData* AliZDCReconstructor::GetCalibData() const
{

  // Getting calibration object for ZDC set

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/Data");
  if(!entry) AliFatal("No calibration data loaded!");  

  AliZDCCalibData *calibdata = dynamic_cast<AliZDCCalibData*>  (entry->GetObject());
  if(!calibdata)  AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}
