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

/* $Id: $ */

// Objects of this class contain basis for reference calibrations
//

#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliEMCALCalibReference.h"

using namespace std;

ClassImp(AliEMCALCalibReference)

//____________________________________________________________________________
AliEMCALCalibReference::AliEMCALCalibReference(const int nSM) : 
  fNSuperModule(nSM),
  fSuperModuleData()
{
  //Default constructor.
  for (int i=0; i<fNSuperModule; i++) {
    fSuperModuleData.Add(new AliEMCALSuperModuleCalibReference(i));
  }
  fSuperModuleData.Compress(); // compress the TObjArray
  fSuperModuleData.SetOwner(kTRUE); 
}

//____________________________________________________________________________
void AliEMCALCalibReference::ReadTextCalibReferenceInfo(Int_t nSM, const TString &txtFileName,
					    Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALCalibReference::ReadCalibReferenceInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t id = 0;

  // list of values to be read
  // first: overall values for the whole SuperModule
  Int_t iReferenceTime = 0; 
  // second: additional info for LED Reference and SM temperature
  Float_t rLEDRefAmp = 0;
  Float_t rLEDRefAmpRMS = 0;
  Int_t iLEDRefHighLow = 0;
  Float_t temperature = 0;
  Float_t temperatureRMS = 0;
  // third: info for each tower
  Int_t iHighLow = 0; // 
  Float_t rLEDAmp = 0; // low gain eq. amplitude
  Float_t rLEDAmpRMS = 0; //
  // end - all values

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibReference * t = (AliEMCALSuperModuleCalibReference*) fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALCalibReference::ReadCalibReferenceInfo - Error while reading input file; likely EOF..\n");
      return;
    }
    inputFile >> iSM;
    t->SetSuperModuleNum(iSM);

    // first: overall values for the whole SuperModule
    inputFile >> iReferenceTime;
    t->SetReferenceTime(iReferenceTime);

    // second: additional info for LED Reference and SM temperature
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      inputFile >> id >> iLEDRefHighLow >> rLEDRefAmp >> rLEDRefAmpRMS;
      if (id<0 || id>(AliEMCALGeoParams::fgkEMCALLEDRefs-1) ) {
	printf("AliEMCALCalibReference::ReadCalibReferenceInfo - Error while reading input file; LEDRef j %d id %d\n", j, id);
	return;
      }
      t->SetLEDRefHighLow(id, iLEDRefHighLow);
      t->SetLEDRefAmp(id, rLEDRefAmp);
      t->SetLEDRefAmpRMS(id, rLEDRefAmpRMS);
    }

    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      inputFile >> id >> temperature >> temperatureRMS;
      if (id<0 || id>(AliEMCALGeoParams::fgkEMCALTempSensors-1) ) {
	printf("AliEMCALCalibReference::ReadCalibReferenceInfo - Error while reading input file; TempSensor j %d id %d\n", j, id);
	return;
      }
      t->SetTemperature(id, temperature);
      t->SetTemperatureRMS(id, temperatureRMS);
    }

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow 
		>> iHighLow >> rLEDAmp >> rLEDAmpRMS;

      // check that input values are not out bounds
      if (iCol<0 || iCol>(AliEMCALGeoParams::fgkEMCALCols-1) ||
	  iRow<0 || iRow>(AliEMCALGeoParams::fgkEMCALRows-1) ) {
	printf("AliEMCALCalibReference::ReadCalibReferenceInfo - Error while reading input file; j %d iCol %d iRow %d\n", j, iCol, iRow);
      return;
      }

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      AliEMCALCalibReferenceVal * v = t->GetAPDVal(iCol, iRow);

      v->SetHighLow(iHighLow);
      v->SetLEDAmp(rLEDAmp);
      v->SetLEDAmpRMS(rLEDAmpRMS);
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibReference::WriteTextCalibReferenceInfo(const TString &txtFileName,
					     Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALCalibReference::WriteCalibReferenceInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibReference * t = (AliEMCALSuperModuleCalibReference*) fSuperModuleData[i];

    // first: overall values for the whole SuperModule
    outputFile << t->GetSuperModuleNum() << endl;
    outputFile << t->GetReferenceTime() << endl;

    // second: additional info for LED Reference and SM temperature
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      outputFile << j << " " << t->GetLEDRefHighLow(j) 
		 << " " << t->GetLEDRefAmp(j) << " " << t->GetLEDRefAmpRMS(j) 
		 << endl;
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      outputFile << j << " " << t->GetTemperature(j) << " " << t->GetTemperatureRMS(j) << endl;
    }

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      AliEMCALCalibReferenceVal * v = t->GetAPDVal(iCol, iRow);

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow 
		 << " " << v->GetHighLow() 
		 << " " << v->GetLEDAmp() 
		 << " " << v->GetLEDAmpRMS() << endl;
    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibReference::ReadRootCalibReferenceInfo(const TString &rootFileName,
					    Bool_t swapSides)
{
  //Read data from root file. ; coordinates given on SuperModule basis
  TFile inputFile(rootFileName, "read");  

  TTree *tree = (TTree*) inputFile.Get("tree");

  ReadTreeCalibReferenceInfo(tree, swapSides);

  inputFile.Close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibReference::ReadTreeCalibReferenceInfo(TTree *tree,
					    Bool_t swapSides)
{
  // how many SuperModule's worth of info do we have?
  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  fNSuperModule = tree->GetEntries();

  Int_t iSM = 0; // SuperModule index
  // list of values to be read
  // first: overall values for the whole SuperModule
  Int_t iReferenceTime= 0; 
  // second: additional info for LED Reference and SM temperature
  Float_t rLEDRefAmp[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t rLEDRefAmpRMS[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Int_t iLEDRefHighLow[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t temperature[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  Float_t temperatureRMS[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  // third: info for each tower
  Int_t iHighLow[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Float_t rLEDAmp[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Float_t rLEDAmpRMS[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(rLEDRefAmp, 0, sizeof(rLEDRefAmp)); 
  memset(rLEDRefAmpRMS, 0, sizeof(rLEDRefAmpRMS)); 
  memset(iLEDRefHighLow, 0, sizeof(iLEDRefHighLow)); 
  memset(temperature, 0, sizeof(temperature)); 
  memset(temperatureRMS, 0, sizeof(temperatureRMS)); 
  memset(iHighLow, 0, sizeof(iHighLow)); 
  memset(rLEDAmp, 0, sizeof(rLEDAmp)); 
  memset(rLEDAmpRMS, 0, sizeof(rLEDAmpRMS)); 

  // declare the branches
  tree->SetBranchAddress("iSM", &iSM);
  tree->SetBranchAddress("ReferenceTime", &iReferenceTime);
  //
  tree->SetBranchAddress("LEDRefAmp", rLEDRefAmp);
  tree->SetBranchAddress("LEDRefAmpRMS", rLEDRefAmpRMS);
  tree->SetBranchAddress("LEDRefHighLow", iLEDRefHighLow);
  tree->SetBranchAddress("Temperature", temperature);
  tree->SetBranchAddress("TemperatureRMS", temperatureRMS);
  //
  tree->SetBranchAddress("HighLow", iHighLow);
  tree->SetBranchAddress("LEDAmp", rLEDAmp);
  tree->SetBranchAddress("LEDAmpRMS", rLEDAmpRMS);

  // indices for looping over the towers
  Int_t iCol = 0;
  Int_t iRow = 0;

  for (int ient=0; ient<tree->GetEntries(); ient++) {
    tree->GetEntry(ient);

    // assume the index SuperModules come in order: i=iSM
    AliEMCALSuperModuleCalibReference * t = (AliEMCALSuperModuleCalibReference*) fSuperModuleData[iSM];

    t->SetSuperModuleNum(iSM);
    // first, overall values
    t->SetReferenceTime(iReferenceTime);

    // second: additional info for LED references and SM temperatures
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      t->SetLEDRefAmp(j, rLEDRefAmp[j]);
      t->SetLEDRefAmpRMS(j, rLEDRefAmpRMS[j]);
      t->SetLEDRefHighLow(j, iLEDRefHighLow[j]);
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      t->SetTemperature(j, temperature[j]);
      t->SetTemperatureRMS(j, temperatureRMS[j]);
    }

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      // help variables: possibly modified or swapped indices
      int iColMod = iCol;
      int iRowMod = iRow;
      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iColMod = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRowMod = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      AliEMCALCalibReferenceVal * v = t->GetAPDVal(iColMod, iRowMod);

      v->SetHighLow(iHighLow[iCol][iRow]);
      v->SetLEDAmp(rLEDAmp[iCol][iRow]);
      v->SetLEDAmpRMS(rLEDAmpRMS[iCol][iRow]);
    }

  } // loop over entries

  return;
}

//____________________________________________________________________________
void AliEMCALCalibReference::WriteRootCalibReferenceInfo(const TString &rootFileName,
					     Bool_t swapSides)
{
  // write data to root file. ; coordinates given on SuperModule basis
  TFile destFile(rootFileName, "recreate");  
  if (destFile.IsZombie()) {
    return;
  }  
  destFile.cd();

  TTree *tree = new TTree("tree","");

  // variables for filling the TTree
  Int_t iSM = 0; // SuperModule index
  // list of values to be written
  // first: overall values for the whole SuperModule
  Int_t iReferenceTime = 0; 
  // second: additional info for LED Reference and SM temperature
  Float_t rLEDRefAmp[AliEMCALGeoParams::fgkEMCALLEDRefs] = {0};
  Float_t rLEDRefAmpRMS[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Int_t iLEDRefHighLow[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t temperature[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  Float_t temperatureRMS[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  // third: info for each tower
  Int_t iHighLow[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Float_t rLEDAmp[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Float_t rLEDAmpRMS[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(rLEDRefAmp, 0, sizeof(rLEDRefAmp)); 
  memset(rLEDRefAmpRMS, 0, sizeof(rLEDRefAmpRMS)); 
  memset(iLEDRefHighLow, 0, sizeof(iLEDRefHighLow)); 
  memset(temperature, 0, sizeof(temperature)); 
  memset(temperatureRMS, 0, sizeof(temperatureRMS)); 
  memset(iHighLow, 0, sizeof(iHighLow)); 
  memset(rLEDAmp, 0, sizeof(rLEDAmp)); 
  memset(rLEDAmpRMS, 0, sizeof(rLEDAmpRMS)); 

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  // for looping over towers
  Int_t iCol = 0;
  Int_t iRow = 0;

  // declare the branches
  // first
  tree->Branch("iSM", &iSM, "iSM/I");
  tree->Branch("ReferenceTime", &iReferenceTime, "ReferenceTime/I");
  // second  
  tree->Branch( "LEDRefAmp", &rLEDRefAmp, Form("LEDRefAmp[%d]/F", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "LEDRefAmpRMS", &rLEDRefAmpRMS, Form("LEDRefAmpRMS[%d]/F", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "LEDRefHighLow", &iLEDRefHighLow, Form("LEDRefHighLow[%d]/I", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "Temperature", &temperature, Form("Temperature[%d]/F", AliEMCALGeoParams::fgkEMCALTempSensors) );
  tree->Branch( "TemperatureRMS", &temperatureRMS, Form("TemperatureRMS[%d]/F", AliEMCALGeoParams::fgkEMCALTempSensors) );
  // third: info for each tower; see if a 2D array works OK or if we'll have to use 1D arrays instead 
  tree->Branch( "HighLow", &iHighLow, Form("HighLow[%d][%d]/I", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );
  tree->Branch( "LEDAmp", &rLEDAmp, Form("LEDAmp[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );
  tree->Branch( "LEDAmpRMS", &rLEDAmpRMS, Form("LEDAmpRMS[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );

  for (iSM = 0; iSM < fNSuperModule; iSM++) {
    AliEMCALSuperModuleCalibReference * t = (AliEMCALSuperModuleCalibReference*) fSuperModuleData[iSM];

    iSM = t->GetSuperModuleNum();
    // first, overall values
    iReferenceTime = t->GetReferenceTime();

    // second: additional info for LED references and SM temperatures
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      rLEDRefAmp[j] = t->GetLEDRefAmp(j);
      rLEDRefAmpRMS[j] = t->GetLEDRefAmpRMS(j);
      iLEDRefHighLow[j] = t->GetLEDRefHighLow(j);
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      temperature[j] = t->GetTemperature(j);
      temperatureRMS[j] = t->GetTemperatureRMS(j);
    }

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      // help variables: possibly modified or swapped indices
      int iColMod = iCol;
      int iRowMod = iRow;
      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iColMod = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRowMod = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      AliEMCALCalibReferenceVal * v = t->GetAPDVal(iCol, iRow);

      iHighLow[iColMod][iRowMod] = v->GetHighLow();
      rLEDAmp[iColMod][iRowMod] = v->GetLEDAmp();
      rLEDAmpRMS[iColMod][iRowMod] = v->GetLEDAmpRMS();
    }

    tree->Fill();
  } // i, SuperModule

  tree->Write();
  destFile.Close();

  return;
}

//____________________________________________________________________________
AliEMCALCalibReference::~AliEMCALCalibReference()
{
  fSuperModuleData.Delete();
}

//____________________________________________________________________________
AliEMCALSuperModuleCalibReference * AliEMCALCalibReference::GetSuperModuleCalibReferenceNum(Int_t supModIndex)const
{ // getter via index
  for (int i=0; i<fNSuperModule; i++) {
    AliEMCALSuperModuleCalibReference * t = (AliEMCALSuperModuleCalibReference*) fSuperModuleData[i];
    if (t->GetSuperModuleNum() == supModIndex) {
      return t;
    }
  }

  // if we arrived here, then nothing was found.. just return a NULL pointer 
  return NULL;
}

