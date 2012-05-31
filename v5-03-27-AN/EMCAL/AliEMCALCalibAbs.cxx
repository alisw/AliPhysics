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
  fSuperModuleData.SetOwner(kTRUE); 
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

  // list of values to be read
  // first: overall values for the whole SuperModule
  Int_t iCalibMethod = 0; 
  Int_t iCalibPass = 0; 
  Float_t absoluteCalib = 0; 
  // third: info for each tower
  Float_t relativeCalib = 0; // (ADC>GeV relative gain/conversion), value around 1
  // end - all values

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALCalibAbs::ReadCalibAbsInfo - Error while reading input file; likely EOF..\n");
      return;
    }
    inputFile >> iSM;
    t->SetSuperModuleNum(iSM);

    // first: overall values for the whole SuperModule
    inputFile >> iCalibMethod >> iCalibPass >> absoluteCalib;
    t->SetCalibMethod(iCalibMethod);
    t->SetCalibPass(iCalibPass);
    t->SetAbsoluteCalib(absoluteCalib);

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow >> relativeCalib;

      // check that input values are not out bounds
      if (iCol<0 || iCol>(AliEMCALGeoParams::fgkEMCALCols-1) ||
	  iRow<0 || iRow>(AliEMCALGeoParams::fgkEMCALRows-1) ) {
	printf("AliEMCALCalibAbs::ReadCalibAbsInfo - Error while reading input file; j %d iCol %d iRow %d\n", j, iCol, iRow);
      return;
      }

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      t->SetRelativeCalib(iCol, iRow, relativeCalib);
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
  Float_t relativeCalib = 0;
  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[i];

    // first: overall values for the whole SuperModule
    outputFile << t->GetSuperModuleNum() << endl;
    outputFile << t->GetCalibMethod() << " " 
	       << t->GetCalibPass() << " " 
	       << t->GetAbsoluteCalib() << endl;

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      relativeCalib = t->GetRelativeCalib(iCol, iRow);

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow 
		 << " " << relativeCalib << endl;
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
  // how many SuperModule's worth of info do we have?
  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  fNSuperModule = tree->GetEntries();

  Int_t iSM = 0; // SuperModule index
  // list of values to be read
  // first: overall values for the whole SuperModule
  Int_t iCalibMethod = 0; 
  Int_t iCalibPass = 0; 
  Float_t absoluteCalib = 0; 
  // third: info for each tower
  Float_t relativeCalib[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(relativeCalib, 0, sizeof(relativeCalib)); 

  // declare the branches
  tree->SetBranchAddress("iSM", &iSM);
  tree->SetBranchAddress("CalibMethod", &iCalibMethod);
  tree->SetBranchAddress("CalibPass", &iCalibPass);
  tree->SetBranchAddress("AbsoluteCalib", &absoluteCalib);
  //
  tree->SetBranchAddress("RelativeCalib", relativeCalib);

  // indices for looping over the towers
  Int_t iCol = 0;
  Int_t iRow = 0;

  for (int ient=0; ient<tree->GetEntries(); ient++) {
    tree->GetEntry(ient);

    // assume the index SuperModules come in order: i=iSM
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[iSM];

    t->SetSuperModuleNum(iSM);
    // first, overall values
    t->SetCalibMethod(iCalibMethod);
    t->SetCalibPass(iCalibPass);
    t->SetAbsoluteCalib(absoluteCalib);

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

      t->SetRelativeCalib(iColMod, iRowMod, relativeCalib[iCol][iRow]);
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
  Int_t iCalibMethod = 0; 
  Int_t iCalibPass = 0; 
  Float_t absoluteCalib = 0; 
  // third: info for each tower
  Float_t relativeCalib[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(relativeCalib, 0, sizeof(relativeCalib)); 

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  // for looping over towers
  Int_t iCol = 0;
  Int_t iRow = 0;

  // declare the branches
  // first
  tree->Branch("iSM", &iSM, "iSM/I");
  tree->Branch("CalibMethod", &iCalibMethod, "CalibMethod/I");
  tree->Branch("CalibPass", &iCalibPass, "CalibPass/I");
  tree->Branch("AbsoluteCalib", &absoluteCalib, "AbsoluteCalib/F");
  // third: info for each tower; see if a 2D array works OK or if we'll have to use 1D arrays instead 
  tree->Branch( "RelativeCalib", &relativeCalib, Form("RelativeCalib[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );

  for (iSM = 0; iSM < fNSuperModule; iSM++) {
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[iSM];

    iSM = t->GetSuperModuleNum();
    // first, overall values
    iCalibMethod = t->GetCalibMethod();
    iCalibPass = t->GetCalibPass();
    absoluteCalib = t->GetAbsoluteCalib();

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

      relativeCalib[iColMod][iRowMod] = t->GetRelativeCalib(iCol, iRow);
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
{ // getter via index
  for (int i=0; i<fNSuperModule; i++) {
    AliEMCALSuperModuleCalibAbs * t = (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[i];
    if (t->GetSuperModuleNum() == supModIndex) {
      return t;
    }
  }

  // if we arrived here, then nothing was found.. just return a NULL pointer 
  return NULL;
}

