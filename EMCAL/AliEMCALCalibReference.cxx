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
  Int_t ReferenceTime; 
  // second: additional info for LED Reference and SM temperature
  Float_t LEDRefAmp;
  Float_t LEDRefAmpRMS;
  Int_t LEDRefHighLow;
  Float_t Temperature;
  Float_t TemperatureRMS;
  // third: info for each tower
  Int_t HighLow; // 
  Float_t LEDAmp; // low gain eq. amplitude
  Float_t LEDAmpRMS; //
  // end - all values

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibReference * t = (AliEMCALSuperModuleCalibReference*) fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALCalibReference::ReadCalibReferenceInfo - Error while reading input file; likely EOF..");
      return;
    }
    inputFile >> iSM;
    t->SetSuperModuleNum(iSM);

    // first: overall values for the whole SuperModule
    inputFile >> ReferenceTime;
    t->SetReferenceTime(ReferenceTime);

    // second: additional info for LED Reference and SM temperature
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      inputFile >> id >> LEDRefHighLow >> LEDRefAmp >> LEDRefAmpRMS;
      t->SetLEDRefHighLow(id, LEDRefHighLow);
      t->SetLEDRefAmp(id, LEDRefAmp);
      t->SetLEDRefAmpRMS(id, LEDRefAmpRMS);
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      inputFile >> id >> Temperature >> TemperatureRMS;
      t->SetTemperature(id, Temperature);
      t->SetTemperatureRMS(id, TemperatureRMS);
    }

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow 
		>> HighLow >> LEDAmp >> LEDAmpRMS;

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      AliEMCALCalibReferenceVal * v = t->GetAPDVal(iCol, iRow);

      v->SetHighLow(HighLow);
      v->SetLEDAmp(LEDAmp);
      v->SetLEDAmpRMS(LEDAmpRMS);
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
  Int_t ReferenceTime= {0}; 
  // second: additional info for LED Reference and SM temperature
  Float_t LEDRefAmp[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t LEDRefAmpRMS[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Int_t LEDRefHighLow[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t Temperature[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  Float_t TemperatureRMS[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  // third: info for each tower
  Int_t HighLow[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Float_t LEDAmp[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Float_t LEDAmpRMS[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(LEDRefAmp, 0, sizeof(LEDRefAmp)); 
  memset(LEDRefAmpRMS, 0, sizeof(LEDRefAmpRMS)); 
  memset(LEDRefHighLow, 0, sizeof(LEDRefHighLow)); 
  memset(Temperature, 0, sizeof(Temperature)); 
  memset(TemperatureRMS, 0, sizeof(TemperatureRMS)); 
  memset(HighLow, 0, sizeof(HighLow)); 
  memset(LEDAmp, 0, sizeof(LEDAmp)); 
  memset(LEDAmpRMS, 0, sizeof(LEDAmpRMS)); 

  // declare the branches
  tree->SetBranchAddress("iSM", &iSM);
  tree->SetBranchAddress("ReferenceTime", &ReferenceTime);
  //
  tree->SetBranchAddress("LEDRefAmp", LEDRefAmp);
  tree->SetBranchAddress("LEDRefAmpRMS", LEDRefAmpRMS);
  tree->SetBranchAddress("LEDRefHighLow", LEDRefHighLow);
  tree->SetBranchAddress("Temperature", Temperature);
  tree->SetBranchAddress("TemperatureRMS", TemperatureRMS);
  //
  tree->SetBranchAddress("HighLow", HighLow);
  tree->SetBranchAddress("LEDAmp", LEDAmp);
  tree->SetBranchAddress("LEDAmpRMS", LEDAmpRMS);

  // indices for looping over the towers
  Int_t iCol = 0;
  Int_t iRow = 0;

  for (int ient=0; ient<tree->GetEntries(); ient++) {
    tree->GetEntry(ient);

    // assume the index SuperModules come in order: i=iSM
    AliEMCALSuperModuleCalibReference * t = (AliEMCALSuperModuleCalibReference*) fSuperModuleData[iSM];

    t->SetSuperModuleNum(iSM);
    // first, overall values
    t->SetReferenceTime(ReferenceTime);

    // second: additional info for LED references and SM temperatures
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      t->SetLEDRefAmp(j, LEDRefAmp[j]);
      t->SetLEDRefAmpRMS(j, LEDRefAmpRMS[j]);
      t->SetLEDRefHighLow(j, LEDRefHighLow[j]);
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      t->SetTemperature(j, Temperature[j]);
      t->SetTemperatureRMS(j, TemperatureRMS[j]);
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

      v->SetHighLow(HighLow[iCol][iRow]);
      v->SetLEDAmp(LEDAmp[iCol][iRow]);
      v->SetLEDAmpRMS(LEDAmpRMS[iCol][iRow]);
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
  Int_t ReferenceTime = 0; 
  // second: additional info for LED Reference and SM temperature
  Float_t LEDRefAmp[AliEMCALGeoParams::fgkEMCALLEDRefs] = {0};
  Float_t LEDRefAmpRMS[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Int_t LEDRefHighLow[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t Temperature[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  Float_t TemperatureRMS[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  // third: info for each tower
  Int_t HighLow[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Float_t LEDAmp[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Float_t LEDAmpRMS[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(LEDRefAmp, 0, sizeof(LEDRefAmp)); 
  memset(LEDRefAmpRMS, 0, sizeof(LEDRefAmpRMS)); 
  memset(LEDRefHighLow, 0, sizeof(LEDRefHighLow)); 
  memset(Temperature, 0, sizeof(Temperature)); 
  memset(TemperatureRMS, 0, sizeof(TemperatureRMS)); 
  memset(HighLow, 0, sizeof(HighLow)); 
  memset(LEDAmp, 0, sizeof(LEDAmp)); 
  memset(LEDAmpRMS, 0, sizeof(LEDAmpRMS)); 

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  // for looping over towers
  Int_t iCol = 0;
  Int_t iRow = 0;

  // declare the branches
  // first
  tree->Branch("iSM", &iSM, "iSM/I");
  tree->Branch("ReferenceTime", &ReferenceTime, "ReferenceTime/I");
  // second  
  tree->Branch( "LEDRefAmp", &LEDRefAmp, Form("LEDRefAmp[%d]/F", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "LEDRefAmpRMS", &LEDRefAmpRMS, Form("LEDRefAmpRMS[%d]/F", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "LEDRefHighLow", &LEDRefHighLow, Form("LEDRefHighLow[%d]/I", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "Temperature", &Temperature, Form("Temperature[%d]/F", AliEMCALGeoParams::fgkEMCALTempSensors) );
  tree->Branch( "TemperatureRMS", &TemperatureRMS, Form("TemperatureRMS[%d]/F", AliEMCALGeoParams::fgkEMCALTempSensors) );
  // third: info for each tower; see if a 2D array works OK or if we'll have to use 1D arrays instead 
  tree->Branch( "HighLow", &HighLow, Form("HighLow[%d][%d]/I", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );
  tree->Branch( "LEDAmp", &LEDAmp, Form("LEDAmp[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );
  tree->Branch( "LEDAmpRMS", &LEDAmpRMS, Form("LEDAmpRMS[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );

  for (iSM = 0; iSM < fNSuperModule; iSM++) {
    AliEMCALSuperModuleCalibReference * t = (AliEMCALSuperModuleCalibReference*) fSuperModuleData[iSM];

    iSM = t->GetSuperModuleNum();
    // first, overall values
    ReferenceTime = t->GetReferenceTime();

    // second: additional info for LED references and SM temperatures
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      LEDRefAmp[j] = t->GetLEDRefAmp(j);
      LEDRefAmpRMS[j] = t->GetLEDRefAmpRMS(j);
      LEDRefHighLow[j] = t->GetLEDRefHighLow(j);
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      Temperature[j] = t->GetTemperature(j);
      TemperatureRMS[j] = t->GetTemperatureRMS(j);
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

      HighLow[iColMod][iRowMod] = v->GetHighLow();
      LEDAmp[iColMod][iRowMod] = v->GetLEDAmp();
      LEDAmpRMS[iColMod][iRowMod] = v->GetLEDAmpRMS();
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
{
  for (int i=0; i<fNSuperModule; i++) {
    AliEMCALSuperModuleCalibReference * t = (AliEMCALSuperModuleCalibReference*) fSuperModuleData[i];
    if (t->GetSuperModuleNum() == supModIndex) {
      return t;
    }
  }

  // if we arrived here, then nothing was found.. just return a NULL pointer 
  return NULL;
}

