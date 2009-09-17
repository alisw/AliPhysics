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

// Objects of this class contain basis for absolute calibrations
//

#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliEMCALCalibAbs.h"

using namespace std;

ClassImp(AliEMCALCalibAbs)

//____________________________________________________________________________
AliEMCALCalibAbs::AliEMCALCalibAbs(const int nSM) : 
  fNSuperModule(nSM),
  fSuperModuleData()
{
  //Default constructor.
  for (int i=0; i<fNSuperModule; i++) {
    fSuperModuleData.Add(new AliEMCALSuperModuleCalibAbs(i));
  }
  fSuperModuleData.Compress(); // compress the TObjArray
}

//____________________________________________________________________________
void AliEMCALCalibAbs::ReadTextCalibAbsInfo(Int_t nSM, const TString &txtFileName,
					    Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALCalibAbs::ReadCalibAbsInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t id = 0;

  // list of values to be read
  // first: overall values for the whole SuperModule
  Int_t CalibMethod; 
  Int_t CalibPass; 
  Int_t CalibTime; 
  Float_t AbsoluteGain; 
  // second: additional info for LED Reference and SM temperature
  Float_t LEDRefAmp;
  Float_t LEDRefAmpRMS;
  Float_t LEDRefHighLowRatio;
  Int_t LEDRefHighLow;
  Float_t Temperature;
  Float_t TemperatureRMS;
  // third: info for each tower
  Float_t RelativeGain; // (ADC>GeV relative gain/conversion), value around 1
  Float_t HighLowRatio; // value around 16 or so
  Int_t HighLow; // 
  Float_t LEDAmp; // low gain eq. amplitude
  Float_t LEDAmpRMS; //
  // end - all values

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALCalibAbs::ReadCalibAbsInfo - Error while reading input file; likely EOF..");
      return;
    }
    inputFile >> iSM;
    t->SetSuperModuleNum(iSM);

    // first: overall values for the whole SuperModule
    inputFile >> CalibMethod >> CalibPass >> CalibTime >> AbsoluteGain;
    t->SetCalibMethod(CalibMethod);
    t->SetCalibPass(CalibPass);
    t->SetCalibTime(CalibTime);
    t->SetAbsoluteGain(AbsoluteGain);

    // second: additional info for LED Reference and SM temperature
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      inputFile >> id >> LEDRefAmp >> LEDRefAmpRMS >> LEDRefHighLowRatio >> LEDRefHighLow;
      t->SetLEDRefAmp(id, LEDRefAmp);
      t->SetLEDRefAmpRMS(id, LEDRefAmpRMS);
      t->SetLEDRefHighLowRatio(id, LEDRefHighLowRatio);
      t->SetLEDRefHighLow(id, LEDRefHighLow);
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      inputFile >> id >> Temperature >> TemperatureRMS;
      t->SetTemperature(id, Temperature);
      t->SetTemperatureRMS(id, TemperatureRMS);
    }

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow 
		>> RelativeGain >> HighLowRatio >> HighLow >> LEDAmp >> LEDAmpRMS;

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      AliEMCALCalibAbsVal * v = t->GetAPDVal(iCol, iRow);

      v->SetRelativeGain(RelativeGain);
      v->SetHighLowRatio(HighLowRatio);
      v->SetHighLow(HighLow);
      v->SetLEDAmp(LEDAmp);
      v->SetLEDAmpRMS(LEDAmpRMS);
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibAbs::WriteTextCalibAbsInfo(const TString &txtFileName,
					     Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALCalibAbs::WriteCalibAbsInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[i];

    // first: overall values for the whole SuperModule
    outputFile << t->GetSuperModuleNum() << endl;
    outputFile << t->GetCalibMethod() << " " 
	       << t->GetCalibPass() << " " 
	       << t->GetCalibTime() << " " 
	       << t->GetAbsoluteGain() << endl;

    // second: additional info for LED Reference and SM temperature
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      outputFile << j << " " << t->GetLEDRefAmp(j) << " " << t->GetLEDRefAmpRMS(j) 
		 << " " << t->GetLEDRefHighLowRatio(j) << " " << t->GetLEDRefHighLow(j) 
		 << endl;
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      outputFile << j << " " << t->GetTemperature(j) << " " << t->GetTemperatureRMS(j) << endl;
    }

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      AliEMCALCalibAbsVal * v = t->GetAPDVal(iCol, iRow);

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow 
		 << " " << v->GetRelativeGain()
		 << " " << v->GetHighLowRatio() 
		 << " " << v->GetHighLow() 
		 << " " << v->GetLEDAmp() 
		 << " " << v->GetLEDAmpRMS() << endl;
    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibAbs::ReadRootCalibAbsInfo(const TString &rootFileName,
					    Bool_t swapSides)
{
  //Read data from root file. ; coordinates given on SuperModule basis
  TFile inputFile(rootFileName, "read");  

  TTree *tree = (TTree*) inputFile.Get("tree");

  ReadTreeCalibAbsInfo(tree, swapSides);

  inputFile.Close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibAbs::ReadTreeCalibAbsInfo(TTree *tree,
					    Bool_t swapSides)
{
  // how many SuperModule's worth of entries / APDs do we have?
  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  fNSuperModule = tree->GetEntries() / nAPDPerSM;

  Int_t iSM = 0; // SuperModule index
  // list of values to be read
  // first: overall values for the whole SuperModule
  Int_t CalibMethod; 
  Int_t CalibPass= {0}; 
  Int_t CalibTime= {0}; 
  Float_t AbsoluteGain= {0}; 
  // second: additional info for LED Reference and SM temperature
  Float_t LEDRefAmp[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t LEDRefAmpRMS[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t LEDRefHighLowRatio[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Int_t LEDRefHighLow[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t Temperature[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  Float_t TemperatureRMS[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  // third: info for each tower
  Float_t RelativeGain[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  Float_t HighLowRatio[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  Int_t HighLow[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  Float_t LEDAmp[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  Float_t LEDAmpRMS[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(LEDRefAmp, 0, sizeof(LEDRefAmp)); 
  memset(LEDRefAmpRMS, 0, sizeof(LEDRefAmpRMS)); 
  memset(LEDRefHighLowRatio, 0, sizeof(LEDRefHighLowRatio)); 
  memset(LEDRefHighLow, 0, sizeof(LEDRefHighLow)); 
  memset(Temperature, 0, sizeof(Temperature)); 
  memset(TemperatureRMS, 0, sizeof(TemperatureRMS)); 
  memset(RelativeGain, 0, sizeof(RelativeGain)); 
  memset(HighLowRatio, 0, sizeof(HighLowRatio)); 
  memset(HighLow, 0, sizeof(HighLow)); 
  memset(LEDAmp, 0, sizeof(LEDAmp)); 
  memset(LEDAmpRMS, 0, sizeof(LEDAmpRMS)); 

  // declare the branches
  tree->SetBranchAddress("iSM", &iSM);
  tree->SetBranchAddress("CalibMethod", &CalibMethod);
  tree->SetBranchAddress("CalibPass", &CalibPass);
  tree->SetBranchAddress("CalibTime", &CalibTime);
  tree->SetBranchAddress("AbsoluteGain", &AbsoluteGain);
  //
  tree->SetBranchAddress("LEDRefAmp", LEDRefAmp);
  tree->SetBranchAddress("LEDRefAmpRMS", LEDRefAmpRMS);
  tree->SetBranchAddress("LEDRefHighLowRatio", LEDRefHighLowRatio);
  tree->SetBranchAddress("LEDRefHighLow", LEDRefHighLow);
  tree->SetBranchAddress("Temperature", Temperature);
  tree->SetBranchAddress("TemperatureRMS", TemperatureRMS);
  //
  tree->SetBranchAddress("RelativeGain", RelativeGain);
  tree->SetBranchAddress("HighLowRatio", HighLowRatio);
  tree->SetBranchAddress("HighLow", HighLow);
  tree->SetBranchAddress("LEDAmp", LEDAmp);
  tree->SetBranchAddress("LEDAmpRMS", LEDAmpRMS);

  // indices for looping over the towers
  Int_t iCol = 0;
  Int_t iRow = 0;

  for (int ient=0; ient<tree->GetEntries(); ient++) {
    tree->GetEntry(ient);

    // assume the index SuperModules come in order: i=iSM
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[iSM];

    t->SetSuperModuleNum(iSM);
    // first, overall values
    t->SetCalibMethod(CalibMethod);
    t->SetCalibPass(CalibPass);
    t->SetCalibTime(CalibTime);
    t->SetAbsoluteGain(AbsoluteGain);

    // second: additional info for LED references and SM temperatures
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      t->SetLEDRefAmp(j, LEDRefAmp[j]);
      t->SetLEDRefAmpRMS(j, LEDRefAmpRMS[j]);
      t->SetLEDRefHighLowRatio(j, LEDRefHighLowRatio[j]);
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

      AliEMCALCalibAbsVal * v = t->GetAPDVal(iCol, iRow);

      v->SetRelativeGain(RelativeGain[iCol][iRow]);
      v->SetHighLowRatio(HighLowRatio[iCol][iRow]);
      v->SetHighLow(HighLow[iCol][iRow]);
      v->SetLEDAmp(LEDAmp[iCol][iRow]);
      v->SetLEDAmpRMS(LEDAmpRMS[iCol][iRow]);
    }

  } // loop over entries

  return;
}

//____________________________________________________________________________
void AliEMCALCalibAbs::WriteRootCalibAbsInfo(const TString &rootFileName,
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
  Int_t CalibMethod = 0; 
  Int_t CalibPass = 0; 
  Int_t CalibTime = 0; 
  Float_t AbsoluteGain = 0; 
  // second: additional info for LED Reference and SM temperature
  Float_t LEDRefAmp[AliEMCALGeoParams::fgkEMCALLEDRefs] = {0};
  Float_t LEDRefAmpRMS[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t LEDRefHighLowRatio[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Int_t LEDRefHighLow[AliEMCALGeoParams::fgkEMCALLEDRefs]= {0};
  Float_t Temperature[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  Float_t TemperatureRMS[AliEMCALGeoParams::fgkEMCALTempSensors]= {0};
  // third: info for each tower
  Float_t RelativeGain[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  Float_t HighLowRatio[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  Int_t HighLow[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  Float_t LEDAmp[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  Float_t LEDAmpRMS[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]= {0}; 
  // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(LEDRefAmp, 0, sizeof(LEDRefAmp)); 
  memset(LEDRefAmpRMS, 0, sizeof(LEDRefAmpRMS)); 
  memset(LEDRefHighLowRatio, 0, sizeof(LEDRefHighLowRatio)); 
  memset(LEDRefHighLow, 0, sizeof(LEDRefHighLow)); 
  memset(Temperature, 0, sizeof(Temperature)); 
  memset(TemperatureRMS, 0, sizeof(TemperatureRMS)); 
  memset(RelativeGain, 0, sizeof(RelativeGain)); 
  memset(HighLowRatio, 0, sizeof(HighLowRatio)); 
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
  tree->Branch("CalibMethod", &CalibMethod, "CalibMethod/I");
  tree->Branch("CalibPass", &CalibPass, "CalibPass/I");
  tree->Branch("CalibTime", &CalibTime, "CalibTime/I");
  tree->Branch("AbsoluteGain", &AbsoluteGain, "AbsoluteGain/F");
  // second  
  tree->Branch( "LEDRefAmp", &LEDRefAmp, Form("LEDRefAmp[%d]/F", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "LEDRefAmpRMS", &LEDRefAmpRMS, Form("LEDRefAmpRMS[%d]/F", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "LEDRefHighLowRatio", &LEDRefHighLowRatio, Form("LEDRefHighLowRatio[%d]/F", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "LEDRefHighLow", &LEDRefHighLow, Form("LEDRefHighLow[%d]/I", AliEMCALGeoParams::fgkEMCALLEDRefs) );
  tree->Branch( "Temperature", &Temperature, Form("Temperature[%d]/F", AliEMCALGeoParams::fgkEMCALTempSensors) );
  tree->Branch( "TemperatureRMS", &TemperatureRMS, Form("TemperatureRMS[%d]/F", AliEMCALGeoParams::fgkEMCALTempSensors) );
  // third: info for each tower; see if a 2D array works OK or if we'll have to use 1D arrays instead 
  tree->Branch( "RelativeGain", &RelativeGain, Form("RelativeGain[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );
  tree->Branch( "HighLowRatio", &HighLowRatio, Form("HighLowRatio[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );
  tree->Branch( "HighLow", &HighLow, Form("HighLow[%d][%d]/I", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );
  tree->Branch( "LEDAmp", &LEDAmp, Form("LEDAmp[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );
  tree->Branch( "LEDAmpRMS", &LEDAmpRMS, Form("LEDAmpRMS[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );

  for (iSM = 0; iSM < fNSuperModule; iSM++) {
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[iSM];

    iSM = t->GetSuperModuleNum();
    // first, overall values
    CalibMethod = t->GetCalibMethod();
    CalibPass = t->GetCalibPass();
    CalibTime = t->GetCalibTime();
    AbsoluteGain = t->GetAbsoluteGain();

    // second: additional info for LED references and SM temperatures
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      LEDRefAmp[j] = t->GetLEDRefAmp(j);
      LEDRefAmpRMS[j] = t->GetLEDRefAmpRMS(j);
      LEDRefHighLowRatio[j] = t->GetLEDRefHighLowRatio(j);
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

      AliEMCALCalibAbsVal * v = t->GetAPDVal(iCol, iRow);

      RelativeGain[iCol][iRow] = v->GetRelativeGain();
      HighLowRatio[iCol][iRow] = v->GetHighLowRatio();
      HighLow[iCol][iRow] = v->GetHighLow();
      LEDAmp[iCol][iRow] = v->GetLEDAmp();
      LEDAmpRMS[iCol][iRow] = v->GetLEDAmpRMS();
    }

    tree->Fill();
  } // i, SuperModule

  tree->Write();
  destFile.Close();

  return;
}

//____________________________________________________________________________
AliEMCALCalibAbs::~AliEMCALCalibAbs()
{
  fSuperModuleData.Delete();
}

//____________________________________________________________________________
AliEMCALSuperModuleCalibAbs * AliEMCALCalibAbs::GetSuperModuleCalibAbsNum(Int_t supModIndex)const
{
  for (int i=0; i<fNSuperModule; i++) {
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[i];
    if (t->GetSuperModuleNum() == supModIndex) {
      return t;
    }
  }

  // if we arrived here, then nothing was found.. just return a NULL pointer 
  return NULL;
}

