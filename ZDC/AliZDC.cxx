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
#include <TGeometry.h>
#include <TNode.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>

// --- AliRoot header files
#include "AliDetector.h"
#include "AliZDC.h"
#include "AliZDCHit.h"
#include "AliZDCSDigit.h"
#include "AliZDCDigit.h"
#include "AliZDCDigitizer.h"
#include "AliZDCRawStream.h"
#include "AliZDCCalibData.h"

#include "AliRawDataHeader.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliMC.h"

 
ClassImp(AliZDC)

AliZDC *gZDC;
 
//_____________________________________________________________________________
AliZDC::AliZDC()
{
  //
  // Default constructor for the Zero Degree Calorimeter base class
  //
  
  fIshunt     = 1;
  fNoShower   = 0;

  fHits       = 0;
  fNhits      = 0;

  fDigits     = 0;
  fNdigits    = 0;

}
 
//_____________________________________________________________________________
AliZDC::AliZDC(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Standard constructor for the Zero Degree Calorimeter base class
  //

  fIshunt   = 1;
  fNoShower = 0;

  // Allocate the hits array  
  fHits   = new TClonesArray("AliZDCHit",1000);
  gAlice->GetMCApp()->AddHitList(fHits);
  
  char sensname[5],senstitle[25];
  sprintf(sensname,"ZDC");
  sprintf(senstitle,"ZDC dummy");
  SetName(sensname); SetTitle(senstitle);

  fDigits     = 0;
  fNdigits    = 0;
  
  gZDC=this;

}
//____________________________________________________________________________ 
AliZDC::~AliZDC()
{
  //
  // ZDC destructor
  //

  fIshunt   = 0;
  gZDC=0;

}
//_____________________________________________________________________________
void AliZDC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // 		Add a ZDC hit to the hit list.
  // -> We make use of 2 array of hits:
  // [1]  fHits (the usual one) that contains hits for each PRIMARY
  // [2]  fStHits that contains hits for each EVENT and is used to
  //	  obtain digits at the end of each event
  //
  
  static Float_t primKinEn, xImpact, yImpact, sFlag;

  AliZDCHit *newquad, *curprimquad;
  newquad = new AliZDCHit(fIshunt, track, vol, hits);
  TClonesArray &lhits = *fHits;
  
  if(fNhits==0){
      // First hit -> setting flag for primary or secondary particle
      Int_t primary = gAlice->GetMCApp()->GetPrimary(track);     
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
   }
   else{       
      newquad->SetPrimKinEn(primKinEn);
      newquad->SetXImpact(xImpact);
      newquad->SetYImpact(yImpact);
      newquad->SetSFlag(sFlag);
   }
 
  Int_t j;
  for(j=0; j<fNhits; j++){
    // If hits are equal (same track, same volume), sum them.
     curprimquad = (AliZDCHit*) lhits[j];
     if(*curprimquad == *newquad){
        *curprimquad = *curprimquad+*newquad;
	delete newquad;
	return;
     } 
  }

    //Otherwise create a new hit
    new(lhits[fNhits]) AliZDCHit(newquad);
    fNhits++;
    
    delete newquad;
}

//_____________________________________________________________________________
void AliZDC::BuildGeometry()
{
  //
  // Build the ROOT TNode geometry for event display 
  // in the Zero Degree Calorimeter
  // This routine is dummy for the moment
  //

  TNode *node, *top;
  TBRIK *brik;
  const int kColorZDC  = kBlue;
  
  //
  top=gAlice->GetGeometry()->GetNode("alice");
  
  // ZDC
    brik = new TBRIK("S_ZDC","ZDC box","void",300,300,5);
    top->cd();
    node = new TNode("ZDC","ZDC","S_ZDC",0,0,600,"");
    node->SetLineColor(kColorZDC);
    fNodes->Add(node);
}

//_____________________________________________________________________________
Int_t AliZDC::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance from the mouse to the Zero Degree Calorimeter
  // Dummy routine
  //
  return 9999;
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
  return -11750.;
}
  

//_____________________________________________________________________________
void AliZDC::MakeBranch(Option_t *opt, const char * /*file*/)
{
  //
  // Create Tree branches for the ZDC
  //

  char branchname[10];
  sprintf(branchname,"%s",GetName());

  const char *cH = strstr(opt,"H");
  
  if (cH && fLoader->TreeH())
   fHits   = new TClonesArray("AliZDCHit",1000); 
  
  AliDetector::MakeBranch(opt);
}

//_____________________________________________________________________________
void AliZDC::Hits2SDigits()
{
  // Create summable digits from hits
  
  if (GetDebug()) printf("\n	Entering AliZDC::Hits2Digits() ");
  
  fLoader->LoadHits("read");
  fLoader->LoadSDigits("recreate");
  AliRunLoader* runLoader = fLoader->GetRunLoader();
  AliZDCSDigit sdigit;
  AliZDCSDigit* psdigit = &sdigit;

  // Event loop
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    Float_t pmCZN = 0, pmCZP = 0, pmQZN[4], pmQZP[4], pmZEM1 = 0, pmZEM2 = 0;
    for (Int_t i = 0; i < 4; i++) pmQZN[i] = pmQZP[i] = 0;

    runLoader->GetEvent(iEvent);
    TTree* treeH = fLoader->TreeH();
    Int_t ntracks = (Int_t) treeH->GetEntries();
    ResetHits();

    // Tracks loop
    Int_t sector[2];
    for (Int_t itrack = 0; itrack < ntracks; itrack++) {
      treeH->GetEntry(itrack);
      for (AliZDCHit* zdcHit = (AliZDCHit*)FirstHit(-1); zdcHit;
                      zdcHit = (AliZDCHit*)NextHit()) { 
		      
	sector[0] = zdcHit->GetVolume(0);
	sector[1] = zdcHit->GetVolume(1);
	if ((sector[1] < 1) || (sector[1] > 4)) {
	  Error("Hits2SDigits", "sector[0] = %d, sector[1] = %d", 
		sector[0], sector[1]);
	  continue;
	}
	Float_t lightQ = zdcHit->GetLightPMQ();
	Float_t lightC = zdcHit->GetLightPMC();
     
	if (sector[0] == 1) {          //ZN 
	  pmCZN += lightC;
	  pmQZN[sector[1]-1] += lightQ;
	} else if (sector[0] == 2) {   //ZP 
	  pmCZP += lightC;
	  pmQZP[sector[1]-1] += lightQ;
	} else if (sector[0] == 3) {   //ZEM 
	  if (sector[1] == 1) pmZEM1 += lightC;
	  else                pmZEM2 += lightQ;
	}
      }//Hits loop
    }

    // create the output tree
    fLoader->MakeTree("S");
    TTree* treeS = fLoader->TreeS();
    const Int_t kBufferSize = 4000;
    treeS->Branch(GetName(), "AliZDCSDigit", &psdigit, kBufferSize);

    // Create sdigits for ZN
    sector[0] = 1; // Detector = ZN
    sector[1] = 0; // Common PM ADC
    new(psdigit) AliZDCSDigit(sector, pmCZN);
    if (pmCZN > 0) treeS->Fill();
    for (Int_t j = 0; j < 4; j++) {
      sector[1] = j+1; // Towers PM ADCs
      new(psdigit) AliZDCSDigit(sector, pmQZN[j]);
      if (pmQZN[j] > 0) treeS->Fill();
    }
  
    // Create sdigits for ZP
    sector[0] = 2; // Detector = ZP
    sector[1] = 0; // Common PM ADC
    new(psdigit) AliZDCSDigit(sector, pmCZP);
    if (pmCZP > 0) treeS->Fill();
    for (Int_t j = 0; j < 4; j++) {
      sector[1] = j+1; // Towers PM ADCs
      new(psdigit) AliZDCSDigit(sector, pmQZP[j]);
      if (pmQZP[j] > 0) treeS->Fill();
    }

    // Create sdigits for ZEM
    sector[0] = 3; 
    sector[1] = 1; // Detector = ZEM1
    new(psdigit) AliZDCSDigit(sector, pmZEM1);
    if (pmZEM1 > 0) treeS->Fill();
    sector[1] = 2; // Detector = ZEM2
    new(psdigit) AliZDCSDigit(sector, pmZEM2);
    if (pmZEM2 > 0) treeS->Fill();

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

  // preliminary format: 12 interger values (ZNC, ZNQ1-4, ZPC, ZPQ1-4, ZEM1,2)
  // For the CAEN module V965 we have an header, the Data Words and an End Of Block
  UInt_t ADCHeader; 
  UInt_t ADCData[24];
  UInt_t ADCEndBlock;

  // load the digits
  fLoader->LoadDigits("read");
  AliZDCDigit digit;
  AliZDCDigit* pdigit = &digit;
  TTree* treeD = fLoader->TreeD();
  if (!treeD) return;
  treeD->SetBranchAddress("ZDC", &pdigit);

  // Fill data array
  // ADC header
  UInt_t ADCHeaderGEO = 0;
  UInt_t ADCHeaderCRATE = 0;
  UInt_t ADCHeaderCNT = (UInt_t) treeD->GetEntries();
    
  ADCHeader = ADCHeaderGEO << 27 | 0x1 << 25 | ADCHeaderCRATE << 16 |
              ADCHeaderCNT << 8 ;

  //printf("ADCHeader = %d\n",ADCHeader);
      
  // ADC data word
  UInt_t ADCDataGEO = ADCHeaderGEO;
  UInt_t ADCDataValue[24];
  UInt_t ADCDataOvFlw[24];
  for(Int_t i = 0; i < 24; i++){
    ADCDataValue[i] = 0;
    ADCDataOvFlw[i] = 0;
  }
  UInt_t ADCDataChannel = 0;
  
  // loop over digits
  for (Int_t iDigit = 0; iDigit < treeD->GetEntries(); iDigit++) {
    treeD->GetEntry(iDigit);
    if (!pdigit) continue;

    //ADC data
    Int_t index = 0;
    if(digit.GetSector(0)!=3){
      index = (digit.GetSector(0)-1) + digit.GetSector(1)*4;
      ADCDataChannel = (digit.GetSector(0)-1)*8 + digit.GetSector(1);
    }
    else {
      index = 19 + digit.GetSector(1);
      ADCDataChannel = 5 + digit.GetSector(1)*8;
    }
     
    if ((index < 0) || (index >= 22)) {
      Error("Digits2Raw", "sector[0] = %d, sector[1] = %d", 
	    digit.GetSector(0), digit.GetSector(1));
      continue;
    }
    
    ADCDataValue[index] = digit.GetADCValue(0);
    if (ADCDataValue[index] > 2047) ADCDataOvFlw[index] = 1; 
    ADCDataValue[index+2] = digit.GetADCValue(1);
    if (ADCDataValue[index+2] > 2047) ADCDataOvFlw[index+2] = 1; 
    
    ADCData[index] =   ADCDataGEO << 27 | ADCDataChannel << 17 | 
                       ADCDataOvFlw[index] << 12 | (ADCDataValue[index] & 0xfff); 
    ADCData[index+2] = ADCDataGEO << 27 | ADCDataChannel << 17 | 0x1 << 16 |
                       ADCDataOvFlw[index+2] << 12 | (ADCDataValue[index+2] & 0xfff);                    
  }
  //for (Int_t i=0;i<24;i++)printf("ADCData[%d] = %d\n",i,ADCData[i]);
  
  // End of Block
  UInt_t ADCEndBlockGEO = ADCHeaderGEO;
  UInt_t ADCEndBlockEvCount = gAlice->GetEventNrInRun();
  
  ADCEndBlock = ADCEndBlockGEO << 27 | 0x1 << 26 | ADCEndBlockEvCount;
  
  //printf("ADCEndBlock = %d\n",ADCEndBlock);


  // open the output file
  char fileName[30];
  sprintf(fileName, "ZDC_%d.ddl", AliZDCRawStream::kDDLOffset);
#ifndef __DECCXX
  ofstream file(fileName, ios::binary);
#else
  ofstream file(fileName);
#endif

  // write the DDL data header
  AliRawDataHeader header;
  header.fSize = sizeof(header) + sizeof(ADCHeader) + 
                 sizeof(ADCData) + sizeof(ADCEndBlock);
  //printf("sizeof header = %d, ADCHeader = %d, ADCData = %d, ADCEndBlock = %d\n",
  //        sizeof(header),sizeof(ADCHeader),sizeof(ADCData),sizeof(ADCEndBlock));
  header.SetAttribute(0);  // valid data
  file.write((char*)(&header), sizeof(header));

  // write the raw data and close the file
  file.write((char*) &ADCHeader, sizeof (ADCHeader));
  file.write((char*)(ADCData), sizeof(ADCData));
  file.write((char*) &ADCEndBlock, sizeof(ADCEndBlock));
  file.close();

  // unload the digits
  fLoader->UnloadDigits();
}

//______________________________________________________________________
void AliZDC::SetTreeAddress(){
  // Set branch address for the Trees.
  // Inputs:
  //      none.
  // Outputs:
  //      none.
  // Return:
  //      none.
  if (fLoader->TreeH() && (fHits == 0x0))
    fHits   = new TClonesArray("AliZDCHit",1000);
      
  AliDetector::SetTreeAddress();
}
 
 
//Calibration methods (by Alberto Colla)
 
 
//________________________________________________________________
void AliZDC::CreateCalibData()
{
  // 
  //if (fCalibData) delete fCalibData; // delete previous version
  fCalibData = new AliZDCCalibData(GetName());
}
//________________________________________________________________
void AliZDC::WriteCalibData(Int_t option)
{
  //
  const int kCompressLevel = 9;
  const char defname[] = "$(ALICE)/AliRoot/data/AliZDCCalib.root";
  char* fnam = gAlice->GetZDCCalibFName();
  if (!fnam || fnam[0]=='\0') {
    fnam = gSystem->ExpandPathName(defname);
    Warning("WriteCalibData","No File Name is provided, using default %s",fnam);
  }
  TFile* cdfile = TFile::Open(fnam,"UPDATE","",kCompressLevel);

  // Writes Calibration Data to current directory. 
  // User MUST take care of corresponding file opening and ->cd()... !!!
  // By default, the object is overwritten. Use 0 option for opposite.
  if (option) option = TObject::kOverwrite;
  if (fCalibData) fCalibData->Write(0,option);
  else if (fCalibData) fCalibData->Write(0,option);

  cdfile->Close();
  delete cdfile;
}

//________________________________________________________________
void AliZDC::LoadCalibData()
{
  //
  char* fnam = gAlice->GetZDCCalibFName();
  if (!fnam || fnam[0]=='\0') return; 
  if (!gAlice->IsFileAccessible(fnam)) {
    Error("LoadCalibData","ZDC Calibration Data file is not accessible, %s",fnam);
    exit(1);
  }
  TFile* cdfile = TFile::Open(fnam);

  // Loads Calibration Data from current directory. 
  // User MUST take care of corresponding file opening and ->cd()...!!!
  //
  if (fCalibData) delete fCalibData; // delete previous version
  TString dtname = "Calib_";
  dtname += GetName();
  fCalibData = (AliZDCCalibData*) gDirectory->Get(dtname.Data());
  if (!fCalibData) { 
    Error("LoadCalibData","No Calibration data found for %s",GetName());
    exit(1);
  }

  cdfile->Close();
  delete cdfile;
}


//Calibration methods (by Alberto Colla)
