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

// Objects of this class contain info on APD bias settings/voltages
//

#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliEMCALBiasAPD.h"

using namespace std;

ClassImp(AliEMCALBiasAPD)

//____________________________________________________________________________
AliEMCALBiasAPD::AliEMCALBiasAPD(const int nSM) : 
  fNSuperModule(nSM), // make space for everyone 
  fSuperModuleData()
{
  //Default constructor.
  for (int i=0; i<fNSuperModule; i++) {
    fSuperModuleData.Add(new AliEMCALSuperModuleBiasAPD(i));
  }
  fSuperModuleData.Compress(); // compress the TObjArray
  fSuperModuleData.SetOwner(kTRUE); 
}

//____________________________________________________________________________
void AliEMCALBiasAPD::ReadTextBiasAPDInfo(Int_t nSM, const TString &txtFileName,
					  Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALBiasAPD::ReadBiasAPDInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t iElecId = 0;
  Int_t iDAC = 0;
  Float_t voltage = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleBiasAPD * t = (AliEMCALSuperModuleBiasAPD*) fSuperModuleData[i];

    if (!inputFile) {
      printf("AliEMCALBiasAPD::ReadBiasAPDInfo - Error while reading input file; likely EOF..\n");
      return;
    }
    inputFile >> iSM;
    t->SetSuperModuleNum(iSM);

    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow >> iElecId >> iDAC >> voltage;

      // check that input values are not out bounds
      if (iCol<0 || iCol>(AliEMCALGeoParams::fgkEMCALCols-1) ||
	  iRow<0 || iRow>(AliEMCALGeoParams::fgkEMCALRows-1) ) {
	printf("AliEMCALBiasAPD::ReadBiasAPDInfo - Error while reading input file; j %d iCol %d iRow %d\n", j, iCol, iRow);
      return;
      }

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      t->SetElecId(iCol, iRow, iElecId);
      t->SetDAC(iCol, iRow, iDAC);
      t->SetVoltage(iCol, iRow, voltage);
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALBiasAPD::WriteTextBiasAPDInfo(const TString &txtFileName,
					   Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALBiasAPD::WriteBiasAPDInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t iElecId = 0;
  Int_t iDAC = 0;
  Float_t voltage = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleBiasAPD * t = (AliEMCALSuperModuleBiasAPD*) fSuperModuleData[i];
    outputFile << t->GetSuperModuleNum() << endl;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      iElecId = t->GetElecId(iCol, iRow);
      iDAC = t->GetDAC(iCol, iRow);
      voltage = t->GetVoltage(iCol, iRow);

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow << " " 
		 << iElecId << " " << iDAC << " "
		 << voltage << endl;
    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALBiasAPD::ReadRootBiasAPDInfo(const TString &rootFileName,
					  Bool_t swapSides)
{
  //Read data from root file. ; coordinates given on SuperModule basis
  TFile inputFile(rootFileName, "read");  

  TTree *tree = (TTree*) inputFile.Get("tree");

  ReadTreeBiasAPDInfo(tree, swapSides);

  inputFile.Close();

  return;
}

//____________________________________________________________________________
void AliEMCALBiasAPD::ReadTreeBiasAPDInfo(TTree *tree,
					  Bool_t swapSides)
{
 // how many SuperModule's worth of entries / APDs do we have?
  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  fNSuperModule = tree->GetEntries() / nAPDPerSM;

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  // list of values to be read
  Int_t iElecId = 0;
  Int_t iDAC = 0;
  Float_t voltage = 0;     
  // end - all values

  // declare the branches
  tree->SetBranchAddress("iSM", &iSM);
  tree->SetBranchAddress("iCol", &iCol);
  tree->SetBranchAddress("iRow", &iRow);
  tree->SetBranchAddress("iElecId", &iElecId);
  tree->SetBranchAddress("iDAC", &iDAC);
  tree->SetBranchAddress("voltage", &voltage);

  for (int ient=0; ient<tree->GetEntries(); ient++) {
    tree->GetEntry(ient);

    // assume the index SuperModules come in order: i=iSM
    AliEMCALSuperModuleBiasAPD * t = (AliEMCALSuperModuleBiasAPD*) fSuperModuleData[iSM];
    t->SetSuperModuleNum(iSM);

    // assume that this info is already swapped and done for this basis?
    if (swapSides) {
      // C side, oriented differently than A side: swap is requested
      iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
      iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
    }

    t->SetElecId(iCol, iRow, iElecId);
    t->SetDAC(iCol, iRow, iDAC);
    t->SetVoltage(iCol, iRow, voltage);

  } // 

  return;
}

//____________________________________________________________________________
void AliEMCALBiasAPD::WriteRootBiasAPDInfo(const TString &rootFileName,
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
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t iElecId = 0;
  Int_t iDAC = 0;
  Float_t voltage = 0;
  // declare the branches
  tree->Branch("iSM", &iSM, "iSM/I");
  tree->Branch("iCol", &iCol, "iCol/I");
  tree->Branch("iRow", &iRow, "iRow/I");
  tree->Branch("iElecId", &iElecId, "iElecId/I");
  tree->Branch("iDAC", &iDAC, "iDAC/I");
  tree->Branch("voltage", &voltage, "voltage/F");

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (iSM = 0; iSM < fNSuperModule; iSM++) {
    AliEMCALSuperModuleBiasAPD * t = (AliEMCALSuperModuleBiasAPD*) fSuperModuleData[iSM];

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      iElecId = t->GetElecId(iCol, iRow);
      iDAC = t->GetDAC(iCol, iRow);
      voltage = t->GetVoltage(iCol, iRow);

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      tree->Fill();
    }

  } // i, SuperModule

  tree->Write();
  destFile.Close();

  return;
}

//____________________________________________________________________________
AliEMCALBiasAPD::~AliEMCALBiasAPD()
{
  fSuperModuleData.Delete();
}

//____________________________________________________________________________
AliEMCALSuperModuleBiasAPD * AliEMCALBiasAPD::GetSuperModuleBiasAPDNum(Int_t supModIndex)const
{ // getter via index
  for (int i=0; i<fNSuperModule; i++) {
    AliEMCALSuperModuleBiasAPD * t = (AliEMCALSuperModuleBiasAPD*) fSuperModuleData[i];
    if (t->GetSuperModuleNum() == supModIndex) {
      return t;
    }
  }

  // if we arrived here, then nothing was found.. just return a NULL pointer 
  return NULL;
}

