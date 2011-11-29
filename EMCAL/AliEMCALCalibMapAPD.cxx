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

// Objects of this class contain info on APD calibration and map info
//

#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliEMCALCalibMapAPD.h"

using namespace std;

ClassImp(AliEMCALCalibMapAPD)

//____________________________________________________________________________
AliEMCALCalibMapAPD::AliEMCALCalibMapAPD(const int nSM) : 
  fNSuperModule(nSM),
  fSuperModuleData()
{
  //Default constructor.
  for (int i=0; i<fNSuperModule; i++) {
    fSuperModuleData.Add(new AliEMCALSuperModuleCalibMapAPD(i));
  }
  fSuperModuleData.Compress(); // compress the TObjArray
  fSuperModuleData.SetOwner(kTRUE); 
}

//____________________________________________________________________________
void AliEMCALCalibMapAPD::ReadTextCalibMapAPDInfo(Int_t nSM, const TString &txtFileName,
						  Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALCalibMapAPD::ReadCalibMapAPDInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  // list of values to be read
  Int_t iHW = 0;
  Int_t iAPDNum = 0;
  Float_t v30 = 0;     
  Float_t par[3] = {0};   
  Float_t parErr[3] = {0}; 
  Int_t iBreakDown = 0;
  Float_t darkCurrent = 0; 
  // end - all values

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibMapAPD * t = (AliEMCALSuperModuleCalibMapAPD*) fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALCalibMapAPD::ReadCalibMapAPDInfo - Error while reading input file; likely EOF..\n");
      return;
    }
    inputFile >> iSM;
    t->SetSuperModuleNum(iSM);

    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow >> iHW 
		>> iAPDNum >> v30 
		>> par[0] >> par[1] >> par[2]
		>> parErr[0] >> parErr[1] >> parErr[2]
		>> iBreakDown >> darkCurrent;

      // check that input values are not out bounds
      if (iCol<0 || iCol>(AliEMCALGeoParams::fgkEMCALCols-1) ||
	  iRow<0 || iRow>(AliEMCALGeoParams::fgkEMCALRows-1) ) {
	printf("AliEMCALCalibMapAPD::ReadCalibMapAPDInfo - Error while reading input file; j %d iCol %d iRow %d\n", j, iCol, iRow);
      return;
      }

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      AliEMCALCalibMapAPDVal * v = t->GetAPDVal(iCol, iRow);

      v->SetHardWareId(iHW);
      v->SetAPDNum(iAPDNum);
      v->SetV30(v30);
      v->SetPar(0, par[0]);
      v->SetPar(1, par[1]);
      v->SetPar(2, par[2]);
      v->SetParErr(0, parErr[0]);
      v->SetParErr(1, parErr[1]);
      v->SetParErr(2, parErr[2]);
      v->SetBreakDown(iBreakDown);
      v->SetDarkCurrent(darkCurrent);
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibMapAPD::WriteTextCalibMapAPDInfo(const TString &txtFileName,
						   Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALCalibMapAPD::WriteCalibMapAPDInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibMapAPD * t = (AliEMCALSuperModuleCalibMapAPD*) fSuperModuleData[i];
    outputFile << t->GetSuperModuleNum() << endl;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      AliEMCALCalibMapAPDVal * v = t->GetAPDVal(iCol, iRow);

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow << " " << v->GetHardWareId() 
		 << " " << v->GetAPDNum() << " " << v->GetV30() 
		 << " " << v->GetPar(0) << " " << v->GetPar(1) << " " << v->GetPar(2)
		 << " " << v->GetParErr(0) << " " << v->GetParErr(1) << " " << v->GetParErr(2)
		 << " " << v->GetBreakDown() << " " << v->GetDarkCurrent() << endl;
    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibMapAPD::ReadRootCalibMapAPDInfo(const TString &rootFileName,
						  Bool_t swapSides)
{
  //Read data from root file. ; coordinates given on SuperModule basis
  TFile inputFile(rootFileName, "read");  

  TTree *tree = (TTree*) inputFile.Get("tree");

  ReadTreeCalibMapAPDInfo(tree, swapSides);

  inputFile.Close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibMapAPD::ReadTreeCalibMapAPDInfo(TTree *tree,
						  Bool_t swapSides)
{
  // how many SuperModule's worth of entries / APDs do we have?
  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  fNSuperModule = tree->GetEntries() / nAPDPerSM;

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  // list of values to be read
  Int_t iHW = 0;
  Int_t iAPDNum = 0;
  Float_t v30 = 0;     
  Float_t par[3] = {0};   
  Float_t parErr[3] = {0}; 
  Int_t iBreakDown = 0;
  Float_t darkCurrent = 0; 
  // end - all values

  // declare the branches
  tree->SetBranchAddress("iSM", &iSM);
  tree->SetBranchAddress("iCol", &iCol);
  tree->SetBranchAddress("iRow", &iRow);
  tree->SetBranchAddress("iHW", &iHW);
  tree->SetBranchAddress("APDNum", &iAPDNum);
  tree->SetBranchAddress("V30", &v30);
  tree->SetBranchAddress("Par", par);
  tree->SetBranchAddress("ParErr", parErr);
  tree->SetBranchAddress("BreakDown", &iBreakDown);
  tree->SetBranchAddress("DarkCurrent", &darkCurrent);

  for (int ient=0; ient<tree->GetEntries(); ient++) {
    tree->GetEntry(ient);

    // assume the index SuperModules come in order: i=iSM
    AliEMCALSuperModuleCalibMapAPD * t = (AliEMCALSuperModuleCalibMapAPD*) fSuperModuleData[iSM];
    t->SetSuperModuleNum(iSM);

    // assume that this info is already swapped and done for this basis?
    if (swapSides) {
      // C side, oriented differently than A side: swap is requested
      iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
      iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
    }

    AliEMCALCalibMapAPDVal * v = t->GetAPDVal(iCol, iRow);

    v->SetHardWareId(iHW);
    v->SetAPDNum(iAPDNum);
    v->SetV30(v30);
    v->SetPar(0, par[0]);
    v->SetPar(1, par[1]);
    v->SetPar(2, par[2]);
    v->SetParErr(0, parErr[0]);
    v->SetParErr(1, parErr[1]);
    v->SetParErr(2, parErr[2]);
    v->SetBreakDown(iBreakDown);
    v->SetDarkCurrent(darkCurrent);
  } // 

  return;
}

//____________________________________________________________________________
void AliEMCALCalibMapAPD::WriteRootCalibMapAPDInfo(const TString &rootFileName,
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
  Int_t iHW = 0;
  Int_t iAPDNum = 0;
  Float_t v30 = 0;     
  Float_t par[3] = {0};   
  Float_t parErr[3] = {0}; 
  Int_t iBreakDown = 0;
  Float_t darkCurrent = 0; 
  //
  Int_t iCol = 0;
  Int_t iRow = 0;
  // declare the branches
  tree->Branch("iSM", &iSM, "iSM/I");
  tree->Branch("iCol", &iCol, "iCol/I");
  tree->Branch("iRow", &iRow, "iRow/I");
  tree->Branch("iHW", &iHW, "iHW/I");
  tree->Branch("APDNum", &iAPDNum, "APDNum/I");
  tree->Branch("V30", &v30, "V30/F");
  tree->Branch("Par", &par, "Par[3]/F");
  tree->Branch("ParErr", &parErr, "ParErr[3]/F");
  tree->Branch("BreakDown", &iBreakDown, "BreakDown/I");
  tree->Branch("DarkCurrent", &darkCurrent, "DarkCurrent/F");

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (iSM = 0; iSM < fNSuperModule; iSM++) {
    AliEMCALSuperModuleCalibMapAPD * t = (AliEMCALSuperModuleCalibMapAPD *) fSuperModuleData[iSM];

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      AliEMCALCalibMapAPDVal * v = t->GetAPDVal(iCol, iRow);

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      iHW = v->GetHardWareId(); 
      iAPDNum = v->GetAPDNum();
      v30 = v->GetV30();
      for (int k=0; k<3; k++) {
	par[k] = v->GetPar(k);
	parErr[k] = v->GetParErr(k);
      } 
      iBreakDown = v->GetBreakDown();
      darkCurrent = v->GetDarkCurrent();

      tree->Fill();
    }

  } // i, SuperModule

  tree->Write();
  destFile.Close();

  return;
}

//____________________________________________________________________________
AliEMCALCalibMapAPD::~AliEMCALCalibMapAPD()
{
  fSuperModuleData.Delete();
}

//____________________________________________________________________________
AliEMCALSuperModuleCalibMapAPD * AliEMCALCalibMapAPD::GetSuperModuleCalibMapAPDNum(Int_t supModIndex)const
{ // getter via index
  for (int i=0; i<fNSuperModule; i++) {
    AliEMCALSuperModuleCalibMapAPD * t = (AliEMCALSuperModuleCalibMapAPD*) fSuperModuleData[i];
    if (t->GetSuperModuleNum() == supModIndex) {
      return t;
    }
  }

  // if we arrived here, then nothing was found.. just return a NULL pointer 
  return NULL;
}

