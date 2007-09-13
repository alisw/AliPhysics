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
    Float_t ZN1TowCorrHG[5], ZP1TowCorrHG[5], ZEMCorrHG=0., 
    	    ZN2TowCorrHG[5], ZP2TowCorrHG[5];
    Float_t ZN1TowCorrLG[5], ZP1TowCorrLG[5], ZEMCorrLG=0., 
            ZN2TowCorrLG[5], ZP2TowCorrLG[5];
    
  for (Int_t iDigit = 0; iDigit < digitsTree->GetEntries(); iDigit++) {
    digitsTree->GetEntry(iDigit);
    if (!pdigit) continue;
      
    Int_t det = digit.GetSector(0);
    Int_t quad = digit.GetSector(1);
    Int_t pedindex;
    //
    if(det == 1){ // *** ZN1
       pedindex = quad;
       ZN1TowCorrHG[quad] = (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]);
       ZN1TowCorrLG[quad] = (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+5]);
    }
    else if(det == 2){ // *** ZP1
       pedindex = quad+10;
       ZP1TowCorrHG[quad] = (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]);
       ZP1TowCorrLG[quad] = (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+5]);
    }
    else if(det == 3){
       if(quad == 1){	    // *** ZEM1  
    	 pedindex = quad+20;
         ZEMCorrHG += (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]); 
         ZEMCorrLG += (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+2]); 
       }
       else if(quad == 2){  // *** ZEM1
    	 pedindex = quad+21;
         ZEMCorrHG += (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]); 
         ZEMCorrLG += (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+2]); 
       }
    }
    else if(det == 4){  // *** ZN2
       pedindex = quad+24;
       ZN2TowCorrHG[quad] = (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]);
       ZN2TowCorrLG[quad] = (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+5]);
    }
    else if(det == 5){  // *** ZP2 
       pedindex = quad+34;
       ZP2TowCorrHG[quad] = (Float_t) (digit.GetADCValue(0)-meanPed[pedindex]);
       ZP2TowCorrLG[quad] = (Float_t) (digit.GetADCValue(1)-meanPed[pedindex+5]);
    }
  }

  // reconstruct the event
    ReconstructEvent(clustersTree, ZN1TowCorrHG, ZP1TowCorrHG, ZN2TowCorrHG, 
    	ZP2TowCorrHG, ZN1TowCorrLG, ZP1TowCorrLG, ZN2TowCorrLG, 
    	ZP2TowCorrLG, ZEMCorrHG);

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
  Float_t ZN1TowCorrHG[5], ZP1TowCorrHG[5], ZEMCorrHG=0., 
  	  ZN2TowCorrHG[5], ZP2TowCorrHG[5];
  Float_t ZN1TowCorrLG[5], ZP1TowCorrLG[5], ZEMCorrLG=0., 
  	  ZN2TowCorrLG[5], ZP2TowCorrLG[5];
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
        if(gain == 0) ZN1TowCorrHG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else ZN1TowCorrLG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+5]); 
      }
      else if(det == 2){ 
        pedindex = quad+10;
        if(gain == 0) ZP1TowCorrHG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else ZP1TowCorrLG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+5]); 
      }
      else if(det == 3){ 
        if(quad==1){	 
          pedindex = quad+20;
          if(gain == 0) ZEMCorrHG += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else ZEMCorrLG += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+2]); 
        }
        else if(quad==2){ 
          pedindex = rawData.GetSector(1)+21;
          if(gain == 0) ZEMCorrHG += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
          else ZEMCorrLG += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+2]); 
        }
      }
      else if(det == 4){       
        pedindex = rawData.GetSector(1)+24;
        if(gain == 0) ZN2TowCorrHG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else ZN2TowCorrLG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+2]); 
      }
      else if(det == 5){
        pedindex = rawData.GetSector(1)+34;
        if(gain == 0) ZP2TowCorrHG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex]); 
        else ZP2TowCorrLG[quad]  += (Float_t) (rawData.GetADCValue()-meanPed[pedindex+5]); 
      }
    }
  }
    
  // reconstruct the event
    ReconstructEvent(clustersTree, ZN1TowCorrHG, ZP1TowCorrHG, ZN2TowCorrHG, 
    	ZP2TowCorrHG, ZN1TowCorrLG, ZP1TowCorrLG, ZN2TowCorrLG, 
    	ZP2TowCorrLG, ZEMCorrHG);

}

//_____________________________________________________________________________
void AliZDCReconstructor::ReconstructEvent(TTree *clustersTree, 
		Float_t* ZN1ADCCorrHG, Float_t* ZP1ADCCorrHG, 
		Float_t* ZN2ADCCorrHG, Float_t* ZP2ADCCorrHG, 
		Float_t* ZN1ADCCorrLG, Float_t* ZP1ADCCorrLG, 
		Float_t* ZN2ADCCorrLG, Float_t* ZP2ADCCorrLG, 
		Float_t ZEMADCCorrHG) const
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
  Float_t ZN1EqualCoeff[5], ZP1EqualCoeff[5], ZN2EqualCoeff[5], ZP2EqualCoeff[5];
  for(Int_t ji=0; ji<5; ji++){
     ZN1EqualCoeff[ji] = fCalibData->GetZN1EqualCoeff(ji);
     ZP1EqualCoeff[ji] = fCalibData->GetZP1EqualCoeff(ji); 
     ZN2EqualCoeff[ji] = fCalibData->GetZN2EqualCoeff(ji); 
     ZP2EqualCoeff[ji] = fCalibData->GetZP2EqualCoeff(ji); 
  }
  //
  Float_t CalibEne[4];
  for(Int_t ij=0; ij<4; ij++) CalibEne[ij] = fCalibData->GetEnCalib(ij);
  //
  Float_t ZEMEndPoint = fCalibData->GetZEMEndValue();
  Float_t ZEMCutFraction = fCalibData->GetZEMCutFraction();
  Float_t DZEMSup = fCalibData->GetDZEMSup();
  Float_t DZEMInf = fCalibData->GetDZEMInf();
  //
  Float_t ZEMCutValue = ZEMEndPoint*ZEMCutFraction;
  Float_t ZEMSupValue = ZEMCutValue+(ZEMEndPoint*DZEMSup);
  Float_t ZEMInfValue = ZEMCutValue-(ZEMEndPoint*DZEMInf);
  //
  Float_t EZN1MaxVal = fCalibData->GetEZN1MaxValue();
  Float_t EZP1MaxVal = fCalibData->GetEZP1MaxValue();
  Float_t EZDC1MaxVal = fCalibData->GetEZDC1MaxValue();
  Float_t EZN2MaxVal = fCalibData->GetEZN1MaxValue();
  Float_t EZP2MaxVal = fCalibData->GetEZP1MaxValue();
  Float_t EZDC2MaxVal = fCalibData->GetEZDC1MaxValue();
  
  // Equalization of detector responses
  Float_t ZN1EqualTowHG[5], ZN2EqualTowHG[5], ZP1EqualTowHG[5], ZP2EqualTowHG[5];
  Float_t ZN1EqualTowLG[5], ZN2EqualTowLG[5], ZP1EqualTowLG[5], ZP2EqualTowLG[5];
  for(Int_t gi=0; gi<5; gi++){
     ZN1EqualTowHG[gi] = ZN1ADCCorrHG[gi]*ZN1EqualCoeff[gi];
     ZP1EqualTowHG[gi] = ZP1ADCCorrHG[gi]*ZP1EqualCoeff[gi];
     ZN2EqualTowHG[gi] = ZN2ADCCorrHG[gi]*ZN2EqualCoeff[gi];
     ZP2EqualTowHG[gi] = ZP2ADCCorrHG[gi]*ZP2EqualCoeff[gi];
     //
     ZN1EqualTowLG[gi] = ZN1ADCCorrLG[gi]*ZN1EqualCoeff[gi];
     ZP1EqualTowLG[gi] = ZP1ADCCorrLG[gi]*ZP1EqualCoeff[gi];
     ZN2EqualTowLG[gi] = ZN2ADCCorrLG[gi]*ZN2EqualCoeff[gi];
     ZP2EqualTowLG[gi] = ZP2ADCCorrLG[gi]*ZP2EqualCoeff[gi];
  }
  
  // Energy calibration of detector responses
  Float_t ZN1CalibTowHG[5], ZN2CalibTowHG[5], ZP1CalibTowHG[5], ZP2CalibTowHG[5];
  Float_t ZN1CalibSumHG=0., ZN2CalibSumHG=0., ZP1CalibSumHG=0., ZP2CalibSumHG=0.;
  Float_t ZN1CalibTowLG[5], ZN2CalibTowLG[5], ZP1CalibTowLG[5], ZP2CalibTowLG[5];
  Float_t ZN1CalibSumLG=0., ZN2CalibSumLG=0., ZP1CalibSumLG=0., ZP2CalibSumLG=0.;
  for(Int_t gi=0; gi<5; gi++){
     ZN1CalibTowHG[gi] = ZN1EqualTowHG[gi]*CalibEne[0];
     ZP1CalibTowHG[gi] = ZP1EqualTowHG[gi]*CalibEne[1];
     ZN2CalibTowHG[gi] = ZN2EqualTowHG[gi]*CalibEne[2];
     ZP2CalibTowHG[gi] = ZP2EqualTowHG[gi]*CalibEne[3];
     ZN1CalibSumHG += ZN1CalibTowHG[gi];
     ZP1CalibSumHG += ZP1CalibTowHG[gi];
     ZN2CalibSumHG += ZN2CalibTowHG[gi];
     ZP2CalibSumHG += ZP2CalibTowHG[gi];
     //
     ZN1CalibTowLG[gi] = ZN1EqualTowLG[gi]*CalibEne[0];
     ZP1CalibTowLG[gi] = ZP1EqualTowLG[gi]*CalibEne[1];
     ZN2CalibTowLG[gi] = ZN2EqualTowLG[gi]*CalibEne[2];
     ZP2CalibTowLG[gi] = ZP2EqualTowLG[gi]*CalibEne[3];
     ZN1CalibSumLG += ZN1CalibTowLG[gi];
     ZP1CalibSumLG += ZP1CalibTowLG[gi];
     ZN2CalibSumLG += ZN2CalibTowLG[gi];
     ZP2CalibSumLG += ZP2CalibTowLG[gi];
  }
  
  //  ---      Number of detected spectator nucleons
  //  *** N.B. -> It works only in Pb-Pb
  Int_t nDetSpecNLeft, nDetSpecPLeft, nDetSpecNRight, nDetSpecPRight;
  nDetSpecNLeft = (Int_t) (ZN1CalibSumHG/2.760);
  nDetSpecPLeft = (Int_t) (ZP1CalibSumHG/2.760);
  nDetSpecNRight = (Int_t) (ZN2CalibSumHG/2.760);
  nDetSpecPRight = (Int_t) (ZP2CalibSumHG/2.760);

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
  if(ZEMADCCorrHG > ZEMSupValue){
    nGenSpecNLeft  = (Int_t) (fZNCen->Eval(ZN1CalibSumHG));
    nGenSpecPLeft  = (Int_t) (fZPCen->Eval(ZP1CalibSumHG));
    nGenSpecLeft   = (Int_t) (fZDCCen->Eval(ZN1CalibSumHG+ZP1CalibSumHG));
    nGenSpecNRight = (Int_t) (fZNCen->Eval(ZN2CalibSumHG));
    nGenSpecPRight = (Int_t) (fZNCen->Eval(ZP2CalibSumHG));
    nGenSpecRight  = (Int_t) (fZNCen->Eval(ZN2CalibSumHG+ZP2CalibSumHG));
    impPar  = fbCen->Eval(ZN1CalibSumHG+ZP1CalibSumHG);
  }
  else if(ZEMADCCorrHG < ZEMInfValue){
    nGenSpecNLeft = (Int_t) (fZNPer->Eval(ZN1CalibSumHG)); 
    nGenSpecPLeft = (Int_t) (fZPPer->Eval(ZP1CalibSumHG));
    nGenSpecLeft  = (Int_t) (fZDCPer->Eval(ZN1CalibSumHG+ZP1CalibSumHG));
    impPar   = fbPer->Eval(ZN1CalibSumHG+ZP1CalibSumHG);
  }
  else if(ZEMADCCorrHG >= ZEMInfValue && ZEMADCCorrHG <= ZEMSupValue){
    nGenSpecNLeft = (Int_t) (fZEMn->Eval(ZEMADCCorrHG));
    nGenSpecPLeft = (Int_t) (fZEMp->Eval(ZEMADCCorrHG));
    nGenSpecLeft  = (Int_t)(fZEMsp->Eval(ZEMADCCorrHG));
    impPar   =  fZEMb->Eval(ZEMADCCorrHG);
  }
  // 
  if(ZN1CalibSumHG/EZN1MaxVal>1.)  nGenSpecNLeft = (Int_t) (fZEMn->Eval(ZEMADCCorrHG));
  if(ZP1CalibSumHG/EZP1MaxVal>1.)  nGenSpecPLeft = (Int_t) (fZEMp->Eval(ZEMADCCorrHG));
  if((ZN1CalibSumHG+ZP1CalibSumHG/EZDC1MaxVal)>1.){
     nGenSpecLeft = (Int_t)(fZEMsp->Eval(ZEMADCCorrHG));
     impPar = fZEMb->Eval(ZEMADCCorrHG);
  }
  if(ZN2CalibSumHG/EZN2MaxVal>1.)  nGenSpecNRight = (Int_t) (fZEMn->Eval(ZEMADCCorrHG));
  if(ZP2CalibSumHG/EZP2MaxVal>1.)  nGenSpecPRight = (Int_t) (fZEMp->Eval(ZEMADCCorrHG));
  if((ZN2CalibSumHG+ZP2CalibSumHG/EZDC2MaxVal)>1.) nGenSpecRight = (Int_t)(fZEMsp->Eval(ZEMADCCorrHG));
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

  // create the output tree
  AliZDCReco reco(ZN1CalibSumHG, ZP1CalibSumHG, ZN2CalibSumHG, ZP2CalibSumHG, 
  		  ZN1CalibTowLG, ZN2CalibTowLG, ZP1CalibTowLG, ZP2CalibTowLG, 
		  ZEMADCCorrHG, 
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
