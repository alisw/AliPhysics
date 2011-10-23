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
#include <TClonesArray.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TParticle.h>

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
#include "AliZDCEnCalib.h"
#include "AliZDCTowerCalib.h"
#include "AliFstream.h"

 
ClassImp(AliZDC)

//_____________________________________________________________________________
AliZDC::AliZDC() :
  AliDetector(),
  fNoShower(0),
  fPedCalib(0),
  fEnCalibData(0),
  fTowCalibData(0),
  fZDCCalibFName(""),
  fSpectatorTracked(1)
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
  fEnCalibData(0),
  fTowCalibData(0),
  fZDCCalibFName(""),
  fSpectatorTracked(1)
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
  
  SetName("ZDC"); SetTitle("ZDC");

}

//____________________________________________________________________________ 
AliZDC::~AliZDC()
{
  //
  // ZDC destructor
  //

  fIshunt = 0;
  if(fPedCalib) delete fPedCalib;
  if(fEnCalibData) delete fEnCalibData;
  if(fEnCalibData) delete fEnCalibData;

}

//_____________________________________________________________________________
AliZDC::AliZDC(const AliZDC& ZDC) :
AliDetector("ZDC","ZDC"),
fNoShower(ZDC.fNoShower),
fPedCalib(ZDC.fPedCalib),
fEnCalibData(ZDC.fEnCalibData),
fTowCalibData(ZDC.fTowCalibData),
fZDCCalibFName(ZDC.fZDCCalibFName),
fSpectatorTracked(ZDC.fSpectatorTracked)
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
    fEnCalibData = ZDC.fEnCalibData;
    fTowCalibData = ZDC.fTowCalibData;
    fZDCCalibFName = ZDC.fZDCCalibFName;
  } return *this;
}

//_____________________________________________________________________________
void AliZDC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // 		Add a ZDC hit to the hit list.
  
  static Float_t trackTime=0., primKinEn=0., xImpact=0., yImpact=0., sFlag=0.;
  static Int_t   pcPDGcode, motPDGcode;

  AliZDCHit *newquad, *curprimquad;
  newquad = new AliZDCHit(fIshunt, track, vol, hits);
  TClonesArray &lhits = *fHits;
  
  if(fNhits==0){
      // First hit -> setting flag for primary or secondary particle
      TParticle * p = gAlice->GetMCApp()->Particle(track);
      Int_t imo = p->GetFirstMother();
      //
      if(track != imo){
        newquad->SetSFlag(1);  // SECONDARY particle entering the ZDC
      }
      else if(track == imo){
        newquad->SetSFlag(0);  // PRIMARY particle entering the ZDC
      }
      //  
      sFlag 	 = newquad->GetSFlag();
      primKinEn  = newquad->GetPrimKinEn();
      xImpact 	 = newquad->GetXImpact();
      yImpact 	 = newquad->GetYImpact();
      pcPDGcode	 = newquad->GetPDGCode();
      motPDGcode = newquad->GetMotherPDGCode();
      trackTime  = newquad->GetTrackTOF();
   }
   else{       
      newquad->SetPrimKinEn(primKinEn);
      newquad->SetXImpact(xImpact);
      newquad->SetYImpact(yImpact);
      newquad->SetSFlag(sFlag);
      newquad->SetPDGCode(pcPDGcode);
      newquad->SetMotherPDGCode(motPDGcode);
      newquad->SetTrackTOF(trackTime);
   }
 
  Int_t j;
  for(j=0; j<fNhits; j++){
    // If hits are equal (same track, same volume), sum them.
     curprimquad = (AliZDCHit*) lhits[j];
     if(*curprimquad == *newquad){
        *curprimquad = *curprimquad+*newquad;
        // Ch. debug
        //printf("\n\t Summing hits **************** \n", fNhits);
        //curprimquad->Print("");
	//
	delete newquad;
	return;
     } 
  }

    //Otherwise create a new hit
    new(lhits[fNhits]) AliZDCHit(*newquad);
    fNhits++;
    // Ch. debug
    //printf("\n\t New ZDC hit added! fNhits = %d\n", fNhits);
    //newquad->Print("");
    
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
  snprintf(branchname, 10, "%s", GetName());

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
  
  AliDebug(1,"\n	AliZDC::Hits2SDigits() ");
  
  fLoader->LoadHits("read");
  fLoader->LoadSDigits("recreate");
  AliRunLoader* runLoader = fLoader->GetRunLoader();
  AliZDCSDigit sdigit;
  AliZDCSDigit* psdigit = &sdigit;

  // Event loop
  for(Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    Float_t pmZNC[5], pmZPC[5], pmZNA[5], pmZPA[5], pmZEM1=0., pmZEM2=0.;
    for(Int_t i=0; i<5; i++) pmZNC[i] = pmZPC[i] =  pmZNA[i] = pmZPA[i] = 0;

    runLoader->GetEvent(iEvent);
    TTree* treeH = fLoader->TreeH();
    Int_t ntracks = (Int_t) treeH->GetEntries();
    ResetHits();

    // Tracks loop
    Int_t sector[2]; Float_t trackTime = 0.;
    for(Int_t itrack = 0; itrack < ntracks; itrack++) {
      treeH->GetEntry(itrack);
      for(AliZDCHit* zdcHit = (AliZDCHit*)FirstHit(-1); zdcHit;
          zdcHit = (AliZDCHit*)NextHit()) { 
		      
	sector[0] = zdcHit->GetVolume(0);
	sector[1] = zdcHit->GetVolume(1);
	if((sector[1] < 1) || (sector[1]>5)) {
	  Error("Hits2SDigits", "sector[0] = %d, sector[1] = %d", sector[0], sector[1]);
	  continue;
	}
	Float_t lightQ = zdcHit->GetLightPMQ();
	Float_t lightC = zdcHit->GetLightPMC();
	trackTime = zdcHit->GetTrackTOF();
	// Signals from ZEM are delayed to arrive in time with ZDC signals
	if(sector[0] == 3)  trackTime += 320;
	// Ch. debug
	//printf("\t det %d vol %d trackTOF %f lightQ %1.0f lightC %1.0f\n",
	//	sector[0], sector[1], trackTime, lightQ, lightC);
     
	if(sector[0] == 1) { //ZNC 
	  pmZNC[0] += lightC;
	  pmZNC[sector[1]] += lightQ;
	} 
	else if(sector[0] == 2) { //ZPC 
	  pmZPC[0] += lightC;
	  pmZPC[sector[1]] += lightQ;
	} 
	else if(sector[0] == 3) { //ZEM 
	  if(sector[1] == 1) pmZEM1 += lightC;
	  else pmZEM2 += lightQ;
	}
	if(sector[0] == 4) { //ZNA 
	  pmZNA[0] += lightC;
	  pmZNA[sector[1]] += lightQ;
	} 
	else if(sector[0] == 5) { //ZPA 
	  pmZPA[0] += lightC;
	  pmZPA[sector[1]] += lightQ;
	} 
      }//Hits loop
    }//Tracks loop

    // create the output tree
    fLoader->MakeTree("S");
    TTree* treeS = fLoader->TreeS();
    const Int_t kBufferSize = 4000;
    treeS->Branch(GetName(), "AliZDCSDigit", &psdigit, kBufferSize);

    // Create sdigits for ZNC
    sector[0] = 1; // Detector = ZNC
    for(Int_t j = 0; j < 5; j++) {
      sector[1] = j; 
      if(pmZNC[j]>0){
        new(psdigit) AliZDCSDigit(sector, pmZNC[j], trackTime);
        treeS->Fill();
	// Ch. debug
	//printf("\t SDigit created: det %d quad %d pmZNC[%d] %1.0f trackTOF %f\n",
	//	sector[0], sector[1], j, pmZNC[j], trackTime);
      }
    }
  
    // Create sdigits for ZPC
    sector[0] = 2; // Detector = ZPC
    for(Int_t j = 0; j < 5; j++) {
      sector[1] = j; // Towers PM ADCs
      if(pmZPC[j]>0){
        new(psdigit) AliZDCSDigit(sector, pmZPC[j], trackTime);
        treeS->Fill();
	// Ch. debug
	//printf("\t SDigit created: det %d quad %d pmZPC[%d] %1.0f trackTOF %f\n",
	//	sector[0], sector[1], j, pmZPC[j], trackTime);
      }
    }

    // Create sdigits for ZEM
    sector[0] = 3; 
    sector[1] = 1; // Detector = ZEM1
    if(pmZEM1>0){ 
      new(psdigit) AliZDCSDigit(sector, pmZEM1, trackTime);
      treeS->Fill();
      // Ch. debug
      //printf("\t SDigit created: det %d quad %d pmZEM1 %1.0f trackTOF %f\n",
      //	sector[0], sector[1], pmZEM1, trackTime);
    }
    sector[1] = 2; // Detector = ZEM2
    if(pmZEM2>0){
      new(psdigit) AliZDCSDigit(sector, pmZEM2, trackTime);
      treeS->Fill();
      // Ch. debug
      //printf("\t SDigit created: det %d quad %d pmZEM2 %1.0f trackTOF %f\n",
      //	sector[0], sector[1], pmZEM2, trackTime);
    }

    // Create sdigits for ZNA
    sector[0] = 4; // Detector = ZNA
    for(Int_t j = 0; j < 5; j++) {
      sector[1] = j; // Towers PM ADCs
      if(pmZNA[j]>0){
        new(psdigit) AliZDCSDigit(sector, pmZNA[j], trackTime);
        treeS->Fill();
	// Ch. debug
	//printf("\t SDigit created: det %d quad %d pmZNA[%d] %1.0f trackTOF %f\n",
	//	sector[0], sector[1], j, pmZNA[j], trackTime);
      }
    }
  
    // Create sdigits for ZPA
    sector[0] = 5; // Detector = ZPA
    sector[1] = 0; // Common PM ADC
    for(Int_t j = 0; j < 5; j++) {
      sector[1] = j; // Towers PM ADCs
      if(pmZPA[j]>0){
        new(psdigit) AliZDCSDigit(sector, pmZPA[j], trackTime);
        treeS->Fill();
	// Ch. debug
	//printf("\t SDigit created: det %d quad %d pmZPA[%d] %1.0f trackTOF %f\n",
	//	sector[0], sector[1], j, pmZPA[j], trackTime);
      }
    }

    // write the output tree
    fLoader->WriteSDigits("OVERWRITE");
  }

  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
}

//_____________________________________________________________________________
AliDigitizer* AliZDC::CreateDigitizer(AliDigitizationInput* digInput) const
{
  // Create the digitizer for ZDC
  AliZDCDigitizer *zdcDigitizer = new AliZDCDigitizer(digInput);
  if(fSpectatorTracked==0) zdcDigitizer->SetSpectators2Track();
  //printf("\n**************************ZDC digitizer created with Spectators2Track = %d\n\n", fSpectatorTracked);
  return zdcDigitizer;
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
  const int knADCData1=12, knADCData2=12; 
  const int knADCData3=12, knADCData4=12; 
  //
  UInt_t lADCHeader1; 
  UInt_t lADCHeader2; 
  UInt_t lADCHeader3; 
  UInt_t lADCHeader4; 
  //
  UInt_t lADCData1[2*knADCData1];
  UInt_t lADCData2[2*knADCData2];
  UInt_t lADCData3[2*knADCData3];
  UInt_t lADCData4[2*knADCData4];
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

  // Reading channel map
  //printf("\n\t Reading ADC mapping from OCDB\n");
  AliZDCChMap * chMap = GetChMap();
  const int nCh = knADCData1+knADCData2+knADCData3+knADCData4;
  Int_t  mapADC[nCh][4]; 
  for(Int_t i=0; i<nCh; i++){
    mapADC[i][0] = chMap->GetADCModule(i);
    mapADC[i][1] = chMap->GetADCChannel(i);
    mapADC[i][2] = chMap->GetDetector(i);
    mapADC[i][3] = chMap->GetSector(i);
    // Ch. debug
    //printf("  mapADC[%d] = (%d %d %d %d)\n", i,
    //	mapADC[i][0],mapADC[i][1],mapADC[i][2],mapADC[i][3]);
  }

  // *** Fill data array
  // ** ADC header
  UInt_t lADCHeaderGEO1 = 0;
  UInt_t lADCHeaderGEO2 = 1;
  UInt_t lADCHeaderGEO3 = 2;
  UInt_t lADCHeaderGEO4 = 3;
  UInt_t lADCHeaderCRATE = 0;
  UInt_t lADCHeaderCNT1 = knADCData1;
  UInt_t lADCHeaderCNT2 = knADCData2;
  UInt_t lADCHeaderCNT3 = knADCData3;
  UInt_t lADCHeaderCNT4 = knADCData4;
    
  lADCHeader1 = lADCHeaderGEO1 << 27 | 0x1 << 25 | lADCHeaderCRATE << 16 |
               lADCHeaderCNT1 << 8 ;
  lADCHeader2 = lADCHeaderGEO2 << 27 | 0x1 << 25 | lADCHeaderCRATE << 16 |
               lADCHeaderCNT2 << 8 ;
  lADCHeader3 = lADCHeaderGEO3 << 27 | 0x1 << 25 | lADCHeaderCRATE << 16 |
               lADCHeaderCNT3 << 8 ;
  lADCHeader4 = lADCHeaderGEO4 << 27 | 0x1 << 25 | lADCHeaderCRATE << 16 |
               lADCHeaderCNT4 << 8 ;
      
  // ** ADC data word
  UInt_t lADCDataGEO = 0;
  //
  UInt_t lADCDataValue1[2*knADCData1];
  UInt_t lADCDataValue2[2*knADCData2];
  UInt_t lADCDataValue3[2*knADCData3];
  UInt_t lADCDataValue4[2*knADCData4];
  //
  UInt_t lADCDataOvFlwHG = 0;
  UInt_t lADCDataOvFlwLG = 0;
  //
  for(Int_t i=0; i<2*knADCData1 ; i++) lADCDataValue1[i] = 0;
  for(Int_t i=0; i<2*knADCData2 ; i++) lADCDataValue2[i] = 0;
  for(Int_t i=0; i<2*knADCData3 ; i++) lADCDataValue3[i] = 0;
  for(Int_t i=0; i<2*knADCData4 ; i++) lADCDataValue4[i] = 0;
  //
  UInt_t lADCDataChannel = 0;
  
  Int_t indADC0=0, indADC1=0, indADC2=0, indADC3=0;
  
  // loop over digits
  for(Int_t iDigit=0; iDigit<(Int_t) (treeD->GetEntries()); iDigit++){
    treeD->GetEntry(iDigit);
    if(!pdigit) continue;
    //digit.Print("");
   
    // *** ADC data
    // Scan of the map to assign the correct ADC module-channel
    for(Int_t k=0; k<nCh; k++){
      if(iDigit<knADCData1+knADCData2){ 
       if(digit.GetSector(0)==mapADC[k][2] && digit.GetSector(1)==mapADC[k][3]){
    	 lADCDataGEO = (UInt_t) mapADC[k][0];
    	 lADCDataChannel = (UInt_t) mapADC[k][1];
    	 break;
       } 
      }
      else{
       if(digit.GetSector(0)==mapADC[k][2] && digit.GetSector(1)==mapADC[k][3]){
    	 lADCDataGEO = (UInt_t) mapADC[k][0];
    	 lADCDataChannel = (UInt_t) mapADC[k][1];
    	 if(k>knADCData1+knADCData2) break;
       } 
      }
    }
    // Ch. debug
    //printf("iDigit %d det %d sec %d -> lADCDataGEO %d  lADCDataChannel %d\n",
    //	iDigit,digit.GetSector(0),digit.GetSector(1),lADCDataGEO,lADCDataChannel);
     
    if(lADCDataGEO==0){ 
      if(indADC0>=knADCData1){
        AliWarning(" Problem with digit index 4 ADC0\n");
	return;
      }
      Int_t indLG = indADC0+knADCData1;
      // High gain ADC ch.	 
      if(digit.GetADCValue(0) > 2047) lADCDataOvFlwHG = 1; 
      lADCDataValue1[indADC0] = digit.GetADCValue(0);    
      lADCData1[indADC0] = lADCDataGEO << 27 |  lADCDataChannel << 17 | 
        	    lADCDataOvFlwHG << 12 | (lADCDataValue1[indADC0] & 0xfff); 
      // Low gain ADC ch.
      if(digit.GetADCValue(1) > 2047) lADCDataOvFlwLG = 1; 
      lADCDataValue1[indLG] = digit.GetADCValue(1); 
      lADCData1[indLG] = lADCDataGEO << 27 |  lADCDataChannel << 17 | 0x1 << 16 |
        	    lADCDataOvFlwLG << 12 | (lADCDataValue1[indLG] & 0xfff);  
      // Ch. debug
      //printf(" lADCDataGEO %d  ADCdataHG[%d] %d  ADCdataLG[%d] %d\n", 
      //  lADCDataGEO,indADC0,lADCDataValue1[indADC0],indLG,lADCDataValue1[indLG]);
		    
      indADC0++;
    }
    else if(lADCDataGEO==1){ 
      if(indADC1>=knADCData2){
         AliWarning(" Problem with digit index 4 ADC1\n");
	 return;
      }
      Int_t indLG = indADC1+knADCData2;
      // High gain ADC ch.	 
      if(digit.GetADCValue(0) > 2047) lADCDataOvFlwHG = 1; 
      lADCDataValue2[indADC1] = digit.GetADCValue(0);	 
      lADCData2[indADC1] = lADCDataGEO << 27 | lADCDataChannel << 17 | 
        	    lADCDataOvFlwHG << 12 | (lADCDataValue2[indADC1] & 0xfff); 
      // Low gain ADC ch.
      if(digit.GetADCValue(1) > 2047) lADCDataOvFlwLG = 1; 
      lADCDataValue2[indLG] = digit.GetADCValue(1); 
      lADCData2[indLG] = lADCDataGEO << 27 |  lADCDataChannel << 17 | 0x1 << 16 |
        	    lADCDataOvFlwLG << 12 | (lADCDataValue2[indLG] & 0xfff);  
      // Ch. debug
      //printf(" lADCDataGEO %d  ADCdataHG[%d] %d  ADCdataLG[%d] %d\n", 
      //  lADCDataGEO,indADC1,lADCDataValue2[indADC1],indLG,lADCDataValue2[indLG]);
        	  
      indADC1++;
    }
    else if(lADCDataGEO==2){ 
      if(indADC2>=knADCData3){
        AliWarning(" Problem with digit index 4 ADC2\n");
	return;
      }
      Int_t indLG = indADC2+knADCData3;
      // High gain ADC ch.	 
      if(digit.GetADCValue(0) > 2047) lADCDataOvFlwHG = 1; 
      lADCDataValue3[indADC1] = digit.GetADCValue(0);    
      lADCData3[indADC1] = lADCDataGEO << 27 | lADCDataChannel << 17 | 
        	    lADCDataOvFlwHG << 12 | (lADCDataValue3[indADC2] & 0xfff); 
      // Low gain ADC ch.
      if(digit.GetADCValue(1) > 2047) lADCDataOvFlwLG = 1; 
      lADCDataValue3[indLG] = digit.GetADCValue(1); 
      lADCData3[indLG] = lADCDataGEO << 27 |  lADCDataChannel << 17 | 0x1 << 16 |
        	    lADCDataOvFlwLG << 12 | (lADCDataValue3[indLG] & 0xfff);  
      // Ch. debug
      //printf(" lADCDataGEO %d  ADCdataHG[%d] %d  ADCdataLG[%d] %d\n", 
      //  lADCDataGEO,indADC2,lADCDataValue3[indADC2],indLG,lADCDataValue3[indLG]);
        	  
      indADC2++;
    }
    else if(lADCDataGEO==3){ 
      if(indADC3>=knADCData4){
         AliWarning(" Problem with digit index 4 ADC2\n");
	 return;
      }
      Int_t indLG = indADC3+knADCData4;
      // High gain ADC ch.	 
      if(digit.GetADCValue(0) > 2047) lADCDataOvFlwHG = 1; 
      lADCDataValue4[indADC3] = digit.GetADCValue(0);    
      lADCData4[indADC3] = lADCDataGEO << 27 | lADCDataChannel << 17 | 
        	    lADCDataOvFlwHG << 12 | (lADCDataValue4[indADC3] & 0xfff); 
      // Low gain ADC ch.
      if(digit.GetADCValue(1) > 2047) lADCDataOvFlwLG = 1; 
      lADCDataValue4[indLG] = digit.GetADCValue(1); 
      lADCData4[indLG] = lADCDataGEO << 27 |  lADCDataChannel << 17 | 0x1 << 16 |
        	    lADCDataOvFlwLG << 12 | (lADCDataValue4[indLG] & 0xfff);  
      // Ch. debug
      //printf(" lADCDataGEO %d  ADCdataHG[%d] %d  ADCdataLG[%d] %d\n", 
      //  lADCDataGEO,indADC3,lADCDataValue4[indADC3],indLG,lADCDataValue4[indLG]);
        	  
      indADC3++;
    }		  

  }
  //
  /*for(Int_t i=0;i<2*knADCData1;i++) printf("\t ADCData1[%d] = %x\n",i,lADCData1[i]);
  for(Int_t i=0;i<2*knADCData2;i++) printf("\t ADCData2[%d] = %x\n",i,lADCData2[i]);
  for(Int_t i=0;i<2*knADCData3;i++) printf("\t ADCData3[%d] = %x\n",i,lADCData3[i]);
  for(Int_t i=0;i<2*knADCData4;i++) printf("\t ADCData4[%d] = %x\n",i,lADCData4[i]);*/
   
  // End of Block
  UInt_t lADCEndBlockGEO = 0;
  // Event counter in ADC EOB -> getting no. of events in run from AliRunLoader
  // get run loader
  AliRunLoader* runLoader = fLoader->GetRunLoader(); 
  UInt_t lADCEndBlockEvCount = runLoader->GetEventNumber();
  //  
  lADCEndBlock = lADCEndBlockGEO << 27 | 0x1 << 26 | lADCEndBlockEvCount;
  //printf("\t AliZDC::Digits2Raw -> ADCEndBlock = %d\n",lADCEndBlock);

  // open the output file
  TString fileName;
  fileName.Form("%s",AliDAQ::DdlFileName("ZDC",0)); 

  AliFstream* file = new AliFstream(fileName.Data());

  // write the DDL data header
  AliRawDataHeaderSim header;
  header.fSize = sizeof(header) + 
                 sizeof(lADCHeader1) + sizeof(lADCData1) + sizeof(lADCEndBlock) +
  		 sizeof(lADCHeader2) + sizeof(lADCData2) + sizeof(lADCEndBlock) +
                 sizeof(lADCHeader3) + sizeof(lADCData3) + sizeof(lADCEndBlock) +
  		 sizeof(lADCHeader4) + sizeof(lADCData4) + sizeof(lADCEndBlock);
  //
  /*printf("sizeof header = %d, ADCHeader1 = %d, ADCData1 = %d, ADCEndBlock = %d\n",
          sizeof(header),sizeof(lADCHeader1),sizeof(lADCData1),sizeof(lADCEndBlock));
  printf("sizeof header = %d, ADCHeader2 = %d, ADCData2 = %d, ADCEndBlock = %d\n",
          sizeof(header),sizeof(lADCHeader2),sizeof(lADCData2),sizeof(lADCEndBlock));
  printf("sizeof header = %d, ADCHeader3 = %d, ADCData3 = %d, ADCEndBlock = %d\n",
          sizeof(header),sizeof(lADCHeader1),sizeof(lADCData1),sizeof(lADCEndBlock));
  printf("sizeof header = %d, ADCHeader4 = %d, ADCData4 = %d, ADCEndBlock = %d\n",
          sizeof(header),sizeof(lADCHeader2),sizeof(lADCData2),sizeof(lADCEndBlock));*/
  
  header.SetAttribute(0);  // valid data
  file->WriteBuffer((char*)(&header), sizeof(header));
  // write the raw data and close the file
  file->WriteBuffer((char*) &lADCHeader1,  sizeof (lADCHeader1));
  file->WriteBuffer((char*) &lADCData1,   sizeof(lADCData1));
  file->WriteBuffer((char*) &lADCEndBlock, sizeof(lADCEndBlock));
  file->WriteBuffer((char*) &lADCHeader2,  sizeof (lADCHeader2));
  file->WriteBuffer((char*) (lADCData2),   sizeof(lADCData2));
  file->WriteBuffer((char*) &lADCEndBlock, sizeof(lADCEndBlock));
  file->WriteBuffer((char*) &lADCHeader3,  sizeof (lADCHeader3));
  file->WriteBuffer((char*) (lADCData3),   sizeof(lADCData3));
  file->WriteBuffer((char*) &lADCEndBlock, sizeof(lADCEndBlock));
  file->WriteBuffer((char*) &lADCHeader4,  sizeof (lADCHeader4));
  file->WriteBuffer((char*) (lADCData4),   sizeof(lADCData4));
  file->WriteBuffer((char*) &lADCEndBlock, sizeof(lADCEndBlock));
  delete file;

  // unload the digits
  fLoader->UnloadDigits();
}

//_____________________________________________________________________________
Bool_t AliZDC::Raw2SDigits(AliRawReader* rawReader)
{
  // Convert ZDC raw data to Sdigits
  const int kNch = 48;
  AliLoader* loader = (AliRunLoader::Instance())->GetLoader("ZDCLoader");
  if(!loader) {
    AliError("no ZDC loader found");
    return kFALSE;
  }

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
      if(jcount < kNch){ 
  	for(Int_t j=0; j<2; j++) sector[j] = rawStream.GetSector(j);
        rawADC = rawStream.GetADCValue();
        resADC = rawStream.GetADCGain();
        //printf("\t RAw2SDigits raw%d ->  RawADC[%d, %d, %d] read\n",
        //    jcount, sector[0], sector[1], rawADC);
        //
        corrADC = rawADC - Pedestal(sector[0], sector[1], resADC);
        if(corrADC<0) corrADC=0;
        nPheVal = ADCch2Phe(sector[0], sector[1], corrADC, resADC);
  	//
        //printf("\t \t ->  SDigit[%d, %d, %d] created\n",
        //    sector[0], sector[1], nPheVal);
        //
  	new(psdigit) AliZDCSDigit(sector, (Float_t) nPheVal, 0.);
  	treeS->Fill();
  	jcount++;
      }
    }//IsADCDataWord
  }//rawStream.Next
  // write the output tree
  fLoader->WriteSDigits("OVERWRITE");
  fLoader->UnloadSDigits();
   
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
  if(!entry) AliFatal("No calibration data loaded!");  
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
  Int_t nPhe = (Int_t) (ADCVal / (pmGain[Det-1][Quad] * resADC[Res]));
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

//_____________________________________________________________________________
AliZDCChMap* AliZDC::GetChMap() const
{

  // Getting calibration object for ZDC

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/ChMap");
  if(!entry) AliFatal("No calibration data loaded!");  

  AliZDCChMap *calibdata = dynamic_cast<AliZDCChMap*> (entry->GetObject());
  if(!calibdata) AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}
