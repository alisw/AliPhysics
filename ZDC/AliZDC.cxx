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
//  			Zero Degree Calorimeter			             //
//  	     This class contains the basic functions for the ZDCs;           //
//            functions specific to one particular geometry are              //
//                      contained in the derived classes                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// --- ROOT system
#include <TBRIK.h>
#include <TClonesArray.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TRandom.h>

// --- AliRoot header files
#include "AliDetector.h"
#include "AliRawDataHeaderSim.h"
#include "AliRawReader.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliLog.h"
#include "AliDAQ.h"
#include "AliZDC.h"
#include "AliZDCHit.h"
#include "AliZDCSDigit.h"
#include "AliZDCDigit.h"
#include "AliZDCDigitizer.h"
#include "AliZDCRawStream.h"
#include "AliZDCPedestals.h"
#include "AliZDCCalib.h"
#include "AliFstream.h"

 
ClassImp(AliZDC)

//_____________________________________________________________________________
AliZDC::AliZDC() :
  AliDetector(),
  fNoShower(0),
  fPedCalib(0),
  fCalibData(0),
  fZDCCalibFName("")
{
  //
  // Default constructor for the Zero Degree Calorimeter base class
  //
  
  fIshunt = 1;
  fNhits  = 0;
  fHits = 0;
  fDigits = 0;
  fNdigits = 0;
  
}
 
//_____________________________________________________________________________
AliZDC::AliZDC(const char *name, const char *title) : 
  AliDetector(name,title),
  fNoShower  (0),
  fPedCalib(0),
  fCalibData(0),
  fZDCCalibFName("")
{
  //
  // Standard constructor for the Zero Degree Calorimeter base class
  //
    
  fIshunt = 1;
  fNhits  = 0;
  fDigits = 0;
  fNdigits = 0;
 
  fHits = new TClonesArray("AliZDCHit",1000);
  gAlice->GetMCApp()->AddHitList(fHits);
  
  char sensname[5],senstitle[25];
  sprintf(sensname,"ZDC");
  sprintf(senstitle,"ZDC dummy");
  SetName(sensname); SetTitle(senstitle);

}

//____________________________________________________________________________ 
AliZDC::~AliZDC()
{
  //
  // ZDC destructor
  //

  fIshunt = 0;
  delete fPedCalib;
  delete fCalibData;

}

//_____________________________________________________________________________
AliZDC::AliZDC(const AliZDC& ZDC) :
AliDetector("ZDC","ZDC"),
fNoShower(ZDC.fNoShower),
fPedCalib(ZDC.fPedCalib),
fCalibData(ZDC.fCalibData),
fZDCCalibFName(ZDC.fZDCCalibFName)
{
  // copy constructor
}

//_____________________________________________________________________________
AliZDC& AliZDC::operator=(const AliZDC& ZDC)
{
  // assignement operator
  if(this!=&ZDC){
    fNoShower = ZDC.fNoShower;
    fPedCalib = ZDC.fPedCalib;
    fCalibData = ZDC.fCalibData;
    fZDCCalibFName = ZDC.fZDCCalibFName;
  } return *this;
}

//_____________________________________________________________________________
void AliZDC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // 		Add a ZDC hit to the hit list.
  
  static Float_t primKinEn=0., xImpact=0., yImpact=0., sFlag=0.;
  static Int_t   pcPDGcode;

  AliZDCHit *newquad, *curprimquad;
  newquad = new AliZDCHit(fIshunt, track, vol, hits);
  TClonesArray &lhits = *fHits;
  
  if(fNhits==0){
      // First hit -> setting flag for primary or secondary particle
      Int_t primary = gAlice->GetMCApp()->GetPrimary(track);     
      //
      if(track != primary){
        newquad->SetSFlag(1);  // SECONDARY particle entering the ZDC
      }
      else if(track == primary){
        newquad->SetSFlag(0);  // PRIMARY particle entering the ZDC
      }  
      sFlag 	= newquad->GetSFlag();
      primKinEn = newquad->GetPrimKinEn();
      xImpact 	= newquad->GetXImpact();
      yImpact 	= newquad->GetYImpact();
      pcPDGcode	= newquad->GetPDGCode();
   }
   else{       
      newquad->SetPrimKinEn(primKinEn);
      newquad->SetXImpact(xImpact);
      newquad->SetYImpact(yImpact);
      newquad->SetSFlag(sFlag);
      newquad->SetPDGCode(pcPDGcode);
   }
 
  Int_t j;
  for(j=0; j<fNhits; j++){
    // If hits are equal (same track, same volume), sum them.
     curprimquad = (AliZDCHit*) lhits[j];
     if(*curprimquad == *newquad){
        *curprimquad = *curprimquad+*newquad;
        // CH. debug
        /*if(newquad->GetEnergy() != 0. || newquad->GetLightPMC() != 0. || 
	   newquad->GetLightPMQ() != 0.){
	  printf("\n\t --- Equal hits found\n");
	  curprimquad->Print("");
	  newquad->Print("");
          printf("\t --- Det. %d, Quad. %d: X = %f, E = %f, LightPMC = %f, LightPMQ = %f\n",
          curprimquad->GetVolume(0),curprimquad->GetVolume(1),curprimquad->GetXImpact(),
          curprimquad->GetEnergy(), curprimquad->GetLightPMC(), curprimquad->GetLightPMQ());
	}*/
	//
	delete newquad;
	return;
     } 
  }

    //Otherwise create a new hit
    new(lhits[fNhits]) AliZDCHit(*newquad);
    fNhits++;
    // CH. debug
    /*printf("\n\t New ZDC hit added! fNhits = %d\n", fNhits);
    printf("\t Det. %d, Quad.t %d: X = %f, E = %f, LightPMC = %f, LightPMQ = %f\n",
    newquad->GetVolume(0),newquad->GetVolume(1),newquad->GetXImpact(),
    newquad->GetEnergy(), newquad->GetLightPMC(), newquad->GetLightPMQ());
    */
    delete newquad;
}

//____________________________________________________________________________
Float_t AliZDC::ZMin(void) const
{
  // Minimum dimension of the ZDC module in z
  return -11600.;
}

//____________________________________________________________________________
Float_t AliZDC::ZMax(void) const
{
  // Maximum dimension of the ZDC module in z
  return 11750.;
}
  

//_____________________________________________________________________________
void AliZDC::MakeBranch(Option_t *opt)
{
  //
  // Create Tree branches for the ZDC
  //

  char branchname[10];
  sprintf(branchname,"%s",GetName());

  const char *cH = strstr(opt,"H");
  
  if(cH && fLoader->TreeH()) {
    if (fHits) {
      fHits->Clear();
      fNhits = 0;
    }
    else {
      fHits   = new TClonesArray("AliZDCHit",1000); 
      if (gAlice && gAlice->GetMCApp())
	gAlice->GetMCApp()->AddHitList(fHits);
    }
  }
  
  AliDetector::MakeBranch(opt);
}

//_____________________________________________________________________________
void AliZDC::Hits2SDigits()
{
  // Create summable digits from hits
  
  AliDebug(1,"\n	Entering AliZDC::Hits2SDigits() ");
  
  fLoader->LoadHits("read");
  fLoader->LoadSDigits("recreate");
  AliRunLoader* runLoader = fLoader->GetRunLoader();
  AliZDCSDigit sdigit;
  AliZDCSDigit* psdigit = &sdigit;

  // Event loop
  for(Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    Float_t pmCZNC=0, pmCZPC=0, pmCZNA=0, pmCZPA=0, pmZEM1 = 0, pmZEM2 = 0;
    Float_t pmQZNC[4], pmQZPC[4], pmQZNA[4], pmQZPA[4];
    for(Int_t i = 0; i < 4; i++) pmQZNC[i] = pmQZPC[i] =  pmQZNA[i] = pmQZPA[i] = 0;

    runLoader->GetEvent(iEvent);
    TTree* treeH = fLoader->TreeH();
    Int_t ntracks = (Int_t) treeH->GetEntries();
    ResetHits();

    // Tracks loop
    Int_t sector[2];
    for(Int_t itrack = 0; itrack < ntracks; itrack++) {
      treeH->GetEntry(itrack);
      for(AliZDCHit* zdcHit = (AliZDCHit*)FirstHit(-1); zdcHit;
                      zdcHit = (AliZDCHit*)NextHit()) { 
		      
	sector[0] = zdcHit->GetVolume(0);
	sector[1] = zdcHit->GetVolume(1);
	if((sector[1] < 1) || (sector[1] > 5)) {
	  Error("Hits2SDigits", "sector[0] = %d, sector[1] = %d", 
		sector[0], sector[1]);
	  continue;
	}
	Float_t lightQ = zdcHit->GetLightPMQ();
	Float_t lightC = zdcHit->GetLightPMC();
     
	if(sector[0] == 1) { //ZNC 
	  pmCZNC += lightC;
	  pmQZNC[sector[1]-1] += lightQ;
	} 
	else if(sector[0] == 2) { //ZPC 
	  pmCZPC += lightC;
	  pmQZPC[sector[1]-1] += lightQ;
	} 
	else if(sector[0] == 3) { //ZEM 
	  if(sector[1] == 1) pmZEM1 += lightC;
	  else pmZEM2 += lightQ;
	}
	if(sector[0] == 4) { //ZNA 
	  pmCZNA += lightC;
	  pmQZNA[sector[1]-1] += lightQ;
	} 
	else if(sector[0] == 5) { //ZPA 
	  pmCZPA += lightC;
	  pmQZPA[sector[1]-1] += lightQ;
	} 
      }//Hits loop
    }

    // create the output tree
    fLoader->MakeTree("S");
    TTree* treeS = fLoader->TreeS();
    const Int_t kBufferSize = 4000;
    treeS->Branch(GetName(), "AliZDCSDigit", &psdigit, kBufferSize);

    // Create sdigits for ZNC
    sector[0] = 1; // Detector = ZNC
    sector[1] = 0; // Common PM ADC
    new(psdigit) AliZDCSDigit(sector, pmCZNC);
    if(pmCZNC > 0) treeS->Fill();
    for(Int_t j = 0; j < 4; j++) {
      sector[1] = j+1; // Towers PM ADCs
      new(psdigit) AliZDCSDigit(sector, pmQZNC[j]);
      if(pmQZNC[j] > 0) treeS->Fill();
    }
  
    // Create sdigits for ZPC
    sector[0] = 2; // Detector = ZPC
    sector[1] = 0; // Common PM ADC
    new(psdigit) AliZDCSDigit(sector, pmCZPC);
    if(pmCZPC > 0) treeS->Fill();
    for(Int_t j = 0; j < 4; j++) {
      sector[1] = j+1; // Towers PM ADCs
      new(psdigit) AliZDCSDigit(sector, pmQZPC[j]);
      if(pmQZPC[j] > 0) treeS->Fill();
    }

    // Create sdigits for ZEM
    sector[0] = 3; 
    sector[1] = 1; // Detector = ZEM1
    new(psdigit) AliZDCSDigit(sector, pmZEM1);
    if(pmZEM1 > 0) treeS->Fill();
    sector[1] = 2; // Detector = ZEM2
    new(psdigit) AliZDCSDigit(sector, pmZEM2);
    if(pmZEM2 > 0) treeS->Fill();

    // Create sdigits for ZNA
    sector[0] = 4; // Detector = ZNA
    sector[1] = 0; // Common PM ADC
    new(psdigit) AliZDCSDigit(sector, pmCZNA);
    if(pmCZNA > 0) treeS->Fill();
    for(Int_t j = 0; j < 4; j++) {
      sector[1] = j+1; // Towers PM ADCs
      new(psdigit) AliZDCSDigit(sector, pmQZNA[j]);
      if(pmQZNA[j] > 0) treeS->Fill();
    }
  
    // Create sdigits for ZPA
    sector[0] = 5; // Detector = ZPA
    sector[1] = 0; // Common PM ADC
    new(psdigit) AliZDCSDigit(sector, pmCZPA);
    if(pmCZPA > 0) treeS->Fill();
    for(Int_t j = 0; j < 4; j++) {
      sector[1] = j+1; // Towers PM ADCs
      new(psdigit) AliZDCSDigit(sector, pmQZPA[j]);
      if(pmQZPA[j] > 0) treeS->Fill();
    }

    // write the output tree
    fLoader->WriteSDigits("OVERWRITE");
  }

  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
}

//_____________________________________________________________________________
AliDigitizer* AliZDC::CreateDigitizer(AliRunDigitizer* manager) const
{
  // Create the digitizer for ZDC

  return new AliZDCDigitizer(manager);
}

//_____________________________________________________________________________
void AliZDC::Digits2Raw()
{
  // Convert ZDC digits to raw data

  // Format: 24 int values -> ZN1(C+Q1-4), ZP1(C+Q1-4), ZEM1, ZEM2, ZN(C+Q1-4), ZP2(C+Q1-4), 2 Ref PMs
  // 	     + 24 int values for the corresponding out of time channels
  // For the CAEN module V965 we have an Header, the Data Words and an End Of Block
  // 	12 channels x 2 gain chains read from 1st ADC module
  // 	12 channels x 2 gain chains read from 2nd ADC module
  // 	12 channels x 2 gain chains read from 3rd ADC module (o.o.t.)
  // 	12 channels x 2 gain chains read from 4rth ADC module (o.o.t.)
  //
  const int knADCData1=24, knADCData2=24; // In principle the 2 numbers can be different!
  UInt_t lADCHeader1; 
  UInt_t lADCHeader2; 
  UInt_t lADCData1[knADCData1];
  UInt_t lADCData2[knADCData2];
  UInt_t lADCData3[knADCData1];
  UInt_t lADCData4[knADCData2];
  //
  UInt_t lADCEndBlock;

  // load the digits
  fLoader->LoadDigits("read");
  AliZDCDigit digit;
  AliZDCDigit* pdigit = &digit;
  TTree* treeD = fLoader->TreeD();
  if(!treeD) return;
  treeD->SetBranchAddress("ZDC", &pdigit);
  //printf("\t AliZDC::Digits2Raw -> TreeD has %d entries\n",(Int_t) treeD->GetEntries());

  // Fill data array
  // ADC header
  UInt_t lADCHeaderGEO = 0;
  UInt_t lADCHeaderCRATE = 0;
  UInt_t lADCHeaderCNT1 = knADCData1;
  UInt_t lADCHeaderCNT2 = knADCData2;
    
  lADCHeader1 = lADCHeaderGEO << 27 | 0x1 << 25 | lADCHeaderCRATE << 16 |
               lADCHeaderCNT1 << 8 ;
  lADCHeader2 = lADCHeaderGEO << 27 | 0x1 << 25 | lADCHeaderCRATE << 16 |
               lADCHeaderCNT2 << 8 ;
      
  // ADC data word
  UInt_t lADCDataGEO = lADCHeaderGEO;
  //
  UInt_t lADCDataValue1[knADCData1];
  UInt_t lADCDataValue2[knADCData2];
  UInt_t lADCDataValue3[knADCData1];
  UInt_t lADCDataValue4[knADCData2];
  //
  UInt_t lADCDataOvFlw1[knADCData1];
  UInt_t lADCDataOvFlw2[knADCData2];
  UInt_t lADCDataOvFlw3[knADCData1];
  UInt_t lADCDataOvFlw4[knADCData2];
  //
  for(Int_t i=0; i<knADCData1 ; i++){
    lADCDataValue1[i] = 0;
    lADCDataOvFlw1[i] = 0;
    lADCDataValue3[i] = 0;
    lADCDataOvFlw3[i] = 0;
  }
  for(Int_t i=0; i<knADCData2 ; i++){
    lADCDataValue2[i] = 0;
    lADCDataOvFlw2[i] = 0;
    lADCDataValue4[i] = 0;
    lADCDataOvFlw4[i] = 0;
  }
  //
  UInt_t lADCDataChannel = 0;
  
  // loop over digits
  for(Int_t iDigit=0; iDigit<treeD->GetEntries(); iDigit++){
    treeD->GetEntry(iDigit);
    if(!pdigit) continue;
    //digit.Print("");
    
    // *** ADC data
    Int_t index=0;
    if(digit.GetSector(1)!=5){ // ZDC signal channels
      // *** ADC1 (ZN1, ZP1, ZEM1,2) or ADC3 (ZN1, ZP1, ZEM1,2 o.o.t.)
      if(digit.GetSector(0)==1 || digit.GetSector(0)==2 || digit.GetSector(0)==3){
        if(digit.GetSector(0)==1 || digit.GetSector(0)==2){
          index = (digit.GetSector(0)-1) + 4*digit.GetSector(1); // ZN1 or ZP1
          lADCDataChannel = 8*(digit.GetSector(0)-1) + digit.GetSector(1);
        }
        else if(digit.GetSector(0)==3){ // ZEM 1,2
          index = 20 + (digit.GetSector(1)-1);
          lADCDataChannel = 5 + 8*(digit.GetSector(1)-1);
        }
        //
        /*printf("\t AliZDC::Digits2Raw -> idig%d det %d quad %d index %d, ADCch %d ADCVal[%d, %d]\n",
		iDigit,digit.GetSector(0),digit.GetSector(1),index,lADCDataChannel,
		digit.GetADCValue(0),digit.GetADCValue(1));// Ch. debug
        */
        //
        if(iDigit<knADCData1){ // *** In-time signals
          lADCDataValue1[index] = digit.GetADCValue(0);   // High gain ADC ch.    
          if(lADCDataValue1[index] > 2047) lADCDataOvFlw1[index] = 1; 
          lADCDataValue1[index+2] = digit.GetADCValue(1); // Low gain ADC ch.
          if(lADCDataValue1[index+2] > 2047) lADCDataOvFlw1[index+2] = 1; 
        
          lADCData1[index] = lADCDataGEO << 27 | 0x1 << 24 | lADCDataChannel << 17 | 
        		  lADCDataOvFlw1[index] << 12 | (lADCDataValue1[index] & 0xfff); 
          lADCData1[index+2] = lADCDataGEO << 27 | 0x1 << 24  | lADCDataChannel << 17 | 0x1 << 16 |
        		  lADCDataOvFlw1[index+2] << 12 | (lADCDataValue1[index+2] & 0xfff);  
        }
      	else{ // *** Out-of-time signals
          lADCDataValue3[index] = digit.GetADCValue(0);   // High gain ADC ch.    
          if(lADCDataValue3[index] > 2047) lADCDataOvFlw3[index] = 1; 
     	  lADCDataValue3[index+2] = digit.GetADCValue(1); // Low gain ADC ch.
     	  if(lADCDataValue3[index+2] > 2047) lADCDataOvFlw3[index+2] = 1; 
      
     	  lADCData3[index] = lADCDataGEO << 27 | lADCDataChannel << 17 | 
     			  lADCDataOvFlw3[index] << 12 | (lADCDataValue3[index] & 0xfff); 
     	  lADCData3[index+2] = lADCDataGEO << 27 | lADCDataChannel << 17 | 0x1 << 16 |
     			  lADCDataOvFlw3[index+2] << 12 | (lADCDataValue3[index+2] & 0xfff);  
     	}		   
      }
      // *** ADC2 (ZN2, ZP2) or ADC4 (ZN2, ZP2 o.o.t.)
      else if(digit.GetSector(0)==4 || digit.GetSector(0)==5){
     	index = (digit.GetSector(0)-4) + 4*digit.GetSector(1); // ZN2 or ZP2
     	lADCDataChannel = 8*(digit.GetSector(0)-4) + digit.GetSector(1);
     	//
        /*printf("\t AliZDC::Digits2Raw -> idig%d det %d quad %d index %d, ADCch %d ADCVal[%d, %d]\n",
		iDigit,digit.GetSector(0),digit.GetSector(1),index,lADCDataChannel,
		digit.GetADCValue(0),digit.GetADCValue(1));// Ch. debug
        */
        //
        if(iDigit<knADCData2){ // *** In-time signals
          lADCDataValue2[index] = digit.GetADCValue(0);
          if(lADCDataValue2[index] > 2047) lADCDataOvFlw2[index] = 1; 
          lADCDataValue2[index+2] = digit.GetADCValue(1);
          if(lADCDataValue2[index+2] > 2047) lADCDataOvFlw2[index+2] = 1; 
          //
          lADCData2[index] =   lADCDataGEO << 27 | lADCDataChannel << 17 | 
        		  lADCDataOvFlw2[index] << 12 | (lADCDataValue2[index] & 0xfff); 
          lADCData2[index+2] = lADCDataGEO << 27 | lADCDataChannel << 17 | 0x1 << 16 |
        		  lADCDataOvFlw2[index+2] << 12 | (lADCDataValue2[index+2] & 0xfff);   
        }		  
      	else{ // *** Out-of-time signals
          lADCDataValue4[index] = digit.GetADCValue(0);
          if(lADCDataValue4[index] > 2047) lADCDataOvFlw4[index] = 1; 
          lADCDataValue4[index+2] = digit.GetADCValue(1);
          if(lADCDataValue4[index+2] > 2047) lADCDataOvFlw4[index+2] = 1; 
          //
          lADCData4[index] =   lADCDataGEO << 27 | lADCDataChannel << 17 | 
                        lADCDataOvFlw4[index] << 12 | (lADCDataValue4[index] & 0xfff); 
          lADCData4[index+2] = lADCDataGEO << 27 | lADCDataChannel << 17 | 0x1 << 16 |
                        lADCDataOvFlw4[index+2] << 12 | (lADCDataValue4[index+2] & 0xfff);   
        }                 
      }
    }
    // *** ADC2 (Reference PTMs) or ADC4 (Reference PTMs o.o.t.)
    else if(digit.GetSector(1)==5){
      index = 20 + (digit.GetSector(0)-1)/3; 
      lADCDataChannel = 5 + 8*(digit.GetSector(0)-1)/3;
      //
      /*printf("\t AliZDC::Digits2Raw -> idig%d det %d quad %d index %d, ADCch %d ADCVal[%d, %d]\n",
		iDigit,digit.GetSector(0),digit.GetSector(1),index,lADCDataChannel,
		digit.GetADCValue(0),digit.GetADCValue(1));// Ch. debug
      */
      //
      if(iDigit<knADCData2){ // *** In-time signals
        lADCDataValue2[index] = digit.GetADCValue(0);
        if(lADCDataValue2[index] > 2047) lADCDataOvFlw2[index] = 1; 
        lADCDataValue2[index+2] = digit.GetADCValue(1);
        if(lADCDataValue2[index+2] > 2047) lADCDataOvFlw2[index+2] = 1; 
        //
        lADCData2[index] =   lADCDataGEO << 27 | lADCDataChannel << 17 | 
        		lADCDataOvFlw2[index] << 12 | (lADCDataValue2[index] & 0xfff); 
        lADCData2[index+2] = lADCDataGEO << 27 | lADCDataChannel << 17 | 0x1 << 16 |
        		lADCDataOvFlw2[index+2] << 12 | (lADCDataValue2[index+2] & 0xfff);   
      } 		
      else{ // *** Out-of-time signals
        lADCDataValue4[index] = digit.GetADCValue(0);
        if(lADCDataValue4[index] > 2047) lADCDataOvFlw4[index] = 1; 
        lADCDataValue4[index+2] = digit.GetADCValue(1);
        if(lADCDataValue4[index+2] > 2047) lADCDataOvFlw4[index+2] = 1; 
        //
        lADCData4[index] =   lADCDataGEO << 27 | lADCDataChannel << 17 | 
        	      lADCDataOvFlw4[index] << 12 | (lADCDataValue4[index] & 0xfff); 
        lADCData4[index+2] = lADCDataGEO << 27 | lADCDataChannel << 17 | 0x1 << 16 |
        	      lADCDataOvFlw4[index+2] << 12 | (lADCDataValue4[index+2] & 0xfff);   
      } 		
           
    }
    if((index<0) || (index>23)) {
      Error("Digits2Raw", "sector[0] = %d, sector[1] = %d", 
	    digit.GetSector(0), digit.GetSector(1));
      continue;
    }
    
    
  }
  //
  /*
  for(Int_t i=0;i<knADCData1;i++) printf("\t ADCData1[%d] = %x\n",i,lADCData1[i]);
  for(Int_t i=0;i<knADCData2;i++) printf("\t ADCData2[%d] = %x\n",i,lADCData2[i]);
  for(Int_t i=0;i<knADCData1;i++) printf("\t ADCData3[%d] = %x\n",i,lADCData3[i]);
  for(Int_t i=0;i<knADCData2;i++) printf("\t ADCData4[%d] = %x\n",i,lADCData4[i]);
  */
 
  // End of Block
  UInt_t lADCEndBlockGEO = lADCHeaderGEO;
  UInt_t lADCEndBlockEvCount = gAlice->GetEventNrInRun();
  //  
  lADCEndBlock = lADCEndBlockGEO << 27 | 0x1 << 26 | lADCEndBlockEvCount;
  //printf("\t AliZDC::Digits2Raw -> ADCEndBlock = %d\n",lADCEndBlock);


  // open the output file
  char fileName[30];
  strcpy(fileName,AliDAQ::DdlFileName("ZDC",0));

  AliFstream* file = new AliFstream(fileName);

  // write the DDL data header
  AliRawDataHeaderSim header;
  header.fSize = sizeof(header) + 
                 sizeof(lADCHeader1) + sizeof(lADCData1) + sizeof(lADCEndBlock) +
  		 sizeof(lADCHeader2) + sizeof(lADCData2) + sizeof(lADCEndBlock) +
                 sizeof(lADCHeader1) + sizeof(lADCData3) + sizeof(lADCEndBlock) +
  		 sizeof(lADCHeader2) + sizeof(lADCData4) + sizeof(lADCEndBlock);
  //
  /*printf("sizeof header = %d, ADCHeader1 = %d, ADCData1 = %d, ADCEndBlock = %d\n",
          sizeof(header),sizeof(lADCHeader1),sizeof(lADCData1),sizeof(lADCEndBlock));
  printf("sizeof header = %d, ADCHeader2 = %d, ADCData2 = %d, ADCEndBlock = %d\n",
          sizeof(header),sizeof(lADCHeader2),sizeof(lADCData2),sizeof(lADCEndBlock));
  */
  //	 
  header.SetAttribute(0);  // valid data
  file->WriteBuffer((char*)(&header), sizeof(header));

  // write the raw data and close the file
  file->WriteBuffer((char*) &lADCHeader1, sizeof (lADCHeader1));
  file->WriteBuffer((char*)(lADCData1), sizeof(lADCData1));
  file->WriteBuffer((char*) &lADCEndBlock, sizeof(lADCEndBlock));
  file->WriteBuffer((char*) &lADCHeader2, sizeof (lADCHeader2));
  file->WriteBuffer((char*)(lADCData2), sizeof(lADCData2));
  file->WriteBuffer((char*) &lADCEndBlock, sizeof(lADCEndBlock));
  file->WriteBuffer((char*) &lADCHeader1, sizeof (lADCHeader1));
  file->WriteBuffer((char*)(lADCData3), sizeof(lADCData3));
  file->WriteBuffer((char*) &lADCEndBlock, sizeof(lADCEndBlock));
  file->WriteBuffer((char*) &lADCHeader2, sizeof (lADCHeader2));
  file->WriteBuffer((char*)(lADCData4), sizeof(lADCData4));
  file->WriteBuffer((char*) &lADCEndBlock, sizeof(lADCEndBlock));
  delete file;

  // unload the digits
  fLoader->UnloadDigits();
}

//_____________________________________________________________________________
Bool_t AliZDC::Raw2SDigits(AliRawReader* rawReader)
{
  // Convert ZDC raw data to Sdigits
  
  AliLoader* loader = (AliRunLoader::GetRunLoader())->GetLoader("ZDCLoader");
  if(!loader) {
    AliError("no ZDC loader found");
    return kFALSE;
  }

//  // Event loop
  Int_t iEvent = 0;
  while(rawReader->NextEvent()){
    (AliRunLoader::GetRunLoader())->GetEvent(iEvent++);
    // Create the output digit tree
    TTree* treeS = loader->TreeS();
    if(!treeS){
      loader->MakeTree("S");
      treeS = loader->TreeS();
    }
    //
    AliZDCSDigit sdigit;
    AliZDCSDigit* psdigit = &sdigit;
    const Int_t kBufferSize = 4000;
    treeS->Branch("ZDC", "AliZDCSDigit",  &psdigit, kBufferSize);
    //
    AliZDCRawStream rawStream(rawReader);
    Int_t sector[2], resADC, rawADC, corrADC, nPheVal;
    Int_t jcount = 0;
    while(rawStream.Next()){
      if(rawStream.IsADCDataWord()){
        //For the moment only in-time SDigits are foreseen (1st 48 raw values)
        if(jcount < 48){ 
          for(Int_t j=0; j<2; j++) sector[j] = rawStream.GetSector(j);
	  rawADC = rawStream.GetADCValue();
	  resADC = rawStream.GetADCGain();
	  //printf("\t RAw2SDigits raw%d ->  RawADC[%d, %d, %d] read\n",
	  //	jcount, sector[0], sector[1], rawADC);
	  //
	  corrADC = rawADC - Pedestal(sector[0], sector[1], resADC);
	  if(corrADC<0) corrADC=0;
	  nPheVal = ADCch2Phe(sector[0], sector[1], corrADC, resADC);
          //
	  //printf("\t \t ->  SDigit[%d, %d, %d] created\n",
	  //	sector[0], sector[1], nPheVal);
	  //
          new(psdigit) AliZDCSDigit(sector, (Float_t) nPheVal);
          treeS->Fill();
          jcount++;
        }
      }//IsADCDataWord
    }//rawStream.Next
    // write the output tree
    fLoader->WriteSDigits("OVERWRITE");
    fLoader->UnloadSDigits();
  }//Event loop 
   
  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliZDC::Pedestal(Int_t Det, Int_t Quad, Int_t Res) const
{
  // Returns a pedestal for detector det, PM quad, channel with res.
  //
  // Getting calibration object for ZDC set
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry  *entry = man->Get("ZDC/Calib/Pedestals");
  AliZDCPedestals *calibPed = (AliZDCPedestals*) entry->GetObject();
  //
  if(!calibPed){
    printf("\t No calibration object found for ZDC!");
    return -1;
  }
  //
  Int_t index=0, kNch=24;
  if(Quad!=5){
    if(Det==1)        index = Quad+kNch*Res;	 // ZN1
    else if(Det==2)   index = Quad+5+kNch*Res;  	 // ZP1
    else if(Det==3)   index = Quad+9+kNch*Res; // ZEM
    else if(Det==4)   index = Quad+12+kNch*Res; // ZN2
    else if(Det==5)   index = Quad+17+kNch*Res; // ZP2
  }
  else index = (Det-1)/3+22+kNch*Res; // Reference PMs
  //
  //
  Float_t meanPed = calibPed->GetMeanPed(index);
  Float_t pedWidth = calibPed->GetMeanPedWidth(index);
  Float_t pedValue = gRandom->Gaus(meanPed,pedWidth);
  //
  //printf("\t AliZDC::Pedestal - det(%d, %d) - Ped[%d] = %d\n",Det, Quad, index,(Int_t) pedValue); // Chiara debugging!
  
  

  return (Int_t) pedValue;
}


//_____________________________________________________________________________
Int_t AliZDC::ADCch2Phe(Int_t Det, Int_t Quad, Int_t ADCVal, Int_t Res) const
{
  // Evaluation of the no. of phe produced
  Float_t pmGain[6][5];
  Float_t resADC[2];
  for(Int_t j = 0; j < 5; j++){
    pmGain[0][j] = 50000.;
    pmGain[1][j] = 100000.;
    pmGain[2][j] = 100000.;
    pmGain[3][j] = 50000.;
    pmGain[4][j] = 100000.;
    pmGain[5][j] = 100000.;
  }
  // ADC Caen V965
  resADC[0] = 0.0000008; // ADC Resolution high gain: 200 fC/adcCh
  resADC[1] = 0.0000064; // ADC Resolution low gain:  25  fC/adcCh
  //
  Int_t nPhe = (Int_t) (ADCVal * pmGain[Det-1][Quad] * resADC[Res]);
  //
  //printf("\t AliZDC::ADCch2Phe -> det(%d, %d) - ADC %d  phe %d\n",Det,Quad,ADCVal,nPhe);

  return nPhe;
}

//______________________________________________________________________
void AliZDC::SetTreeAddress(){

  // Set branch address for the Trees.
  if(fLoader->TreeH() && (fHits == 0x0))
    fHits   = new TClonesArray("AliZDCHit",1000);
      
  AliDetector::SetTreeAddress();
}
