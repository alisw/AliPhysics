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
/*
$Log$
Revision 1.3  2000/06/25 16:47:43  pcrochet
pow replaced by TMath::Power

*/

#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerLut.h"
#include "TTree.h"
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONPoints.h"
#include "TMath.h"
#include "TFile.h"
#include "TH3.h"
#include <iostream.h>

ClassImp(AliMUONTriggerLut)

//----------------------------------------------------------------------
AliMUONTriggerLut::AliMUONTriggerLut() {
// constructor
  fLptPlus = fLptMinu = fLptUnde = 0;
  fHptPlus = fHptMinu = fHptUnde = 0;
  fAptPlus = fAptMinu = fAptUnde = 0;
}
//----------------------------------------------------------------------
AliMUONTriggerLut::~AliMUONTriggerLut() {
// destructor
  delete fLptPlus;  
  delete fLptMinu;
  delete fLptUnde;
  delete fHptPlus;  
  delete fHptMinu;
  delete fHptUnde;
  delete fAptPlus;  
  delete fAptMinu;
  delete fAptUnde;
  fLptPlus = fLptMinu = fLptUnde = 0;
  fHptPlus = fHptMinu = fHptUnde = 0;  
  fAptPlus = fAptMinu = fAptUnde = 0;    
}

//----------------------------------------------------------------------
AliMUONTriggerLut::AliMUONTriggerLut (const AliMUONTriggerLut& MUONTriggerLut)
{
// Dummy copy constructor
}

//----------------------------------------------------------------------
AliMUONTriggerLut & AliMUONTriggerLut::operator=(const AliMUONTriggerLut& MUONTriggerLut)
{
// Dummy assignment operator
    return *this;
}

//----------------------------------------------------------------------
void AliMUONTriggerLut::GetLutOutput(Int_t circuit, Int_t xstrip, Int_t idev,
				     Int_t ystrip, Int_t lutLpt[2], 
				     Int_t lutHpt[2], Int_t lutApt[2]){
// return output of LuT for corresponding TH3S  

  static TFile *fileLut;
  static Bool_t first=kTRUE;  
  if(first) {
    cout << " opening MUONTriggerLut.root " << "\n";
    fileLut = new TFile("$(ALICE_ROOT)/MUON/MUONTriggerLut.root","READ");
    first=kFALSE;
  }
  fileLut->cd();

// get the pointers to the TH3S objects of the file
  TH3S *lptPlus = (TH3S*)gROOT->FindObject("LptPlus");  
  TH3S *lptMinu = (TH3S*)gROOT->FindObject("LptMinu");
  TH3S *lptUnde = (TH3S*)gROOT->FindObject("LptUnde");
  TH3S *hptPlus = (TH3S*)gROOT->FindObject("HptPlus");  
  TH3S *hptMinu = (TH3S*)gROOT->FindObject("HptMinu");
  TH3S *hptUnde = (TH3S*)gROOT->FindObject("HptUnde");
  TH3S *aptPlus = (TH3S*)gROOT->FindObject("AptPlus");  
  TH3S *aptMinu = (TH3S*)gROOT->FindObject("AptMinu");
  TH3S *aptUnde = (TH3S*)gROOT->FindObject("AptUnde");

  Int_t bin;
  Short_t binc; 
  Int_t mask = GetMask(ystrip);        // get ystrip mask
  
  // Low pt.............................................. 
  bin    =          lptPlus->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)lptPlus->GetBinContent(bin);
  if ((binc & mask)!=0) lutLpt[1]=1;

  bin    =          lptMinu->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)lptMinu->GetBinContent(bin);
  if ((binc & mask)!=0) lutLpt[0]=1;
  
  bin    =          lptUnde->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)lptUnde->GetBinContent(bin);
  if ((binc & mask)!=0) lutLpt[0]=lutLpt[1]=1;

  // High pt.............................................
  bin    =          hptPlus->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)hptPlus->GetBinContent(bin);
  if ((binc & mask)!=0) lutHpt[1]=1;

  bin    =          hptMinu->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)hptMinu->GetBinContent(bin);
  if ((binc & mask)!=0) lutHpt[0]=1;

  bin    =          hptUnde->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)hptUnde->GetBinContent(bin);
  if ((binc & mask)!=0) lutHpt[0]=lutHpt[1]=1;

  // All pts.............................................
  bin    =          aptPlus->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)aptPlus->GetBinContent(bin);
  if ((binc & mask)!=0) lutApt[1]=1;

  bin    =          aptMinu->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)aptMinu->GetBinContent(bin);
  if ((binc & mask)!=0) lutApt[0]=1;

  bin    =          aptUnde->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)aptUnde->GetBinContent(bin);
  if ((binc & mask)!=0) lutApt[0]=lutApt[1]=1;

// get back to the first file
  TTree *tK = gAlice->TreeK();
  TFile *file1 = 0;
  if (tK) file1 = tK->GetCurrentFile();
  file1->cd();
}

//----------------------------------------------------------------------
Int_t AliMUONTriggerLut::GetMask(Int_t ystrip){
// returns the mask corresponding to ystrip
  Int_t tabMask[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t mask=0;
  tabMask[ystrip]=1;
  for (Int_t i=0; i<16; i++) {          
    mask=mask+Int_t(tabMask[i]*TMath::Power(2,i));   
  }
  return mask;
}

//----------------------------------------------------------------------
void AliMUONTriggerLut::LoadLut(){
// !!!!!!! This is dummy version of the LoadLut method !!!!!!!
// !!!!!!!         calibration to be done              !!!!!!!
// 1) Loop on circuit/Xstrip1/deviation/Ystrip
// 2) get corresponding ptCal from AliMUONTriggerCircuit
// 3) fill histos with cuts on deviation, ptLow and ptHigh 
// 4) store histos in a file

  char fileName[60];
  sprintf(fileName,"$(ALICE_ROOT)/MUON/MUONTriggerLut.root");
  cout << " file name is " << fileName << "\n";

// open output file containing histos  
  TFile *hfile = new TFile(fileName,"RECREATE","Trigger Look Up Table");

  //..........................................circuit/stripX/deviation
  TH3S *fLptPlus=new TH3S("LptPlus","LptPlus",234,0,234,31,0,31,31,0,31);
  TH3S *fLptMinu=new TH3S("LptMinu","LptMinu",234,0,234,31,0,31,31,0,31);
  TH3S *fLptUnde=new TH3S("LptUnde","LptUnde",234,0,234,31,0,31,31,0,31);

  TH3S *fHptPlus=new TH3S("HptPlus","HptPlus",234,0,234,31,0,31,31,0,31);
  TH3S *fHptMinu=new TH3S("HptMinu","HptMinu",234,0,234,31,0,31,31,0,31);
  TH3S *fHptUnde=new TH3S("HptUnde","HptUnde",234,0,234,31,0,31,31,0,31);

  TH3S *fAptPlus=new TH3S("AptPlus","AptPlus",234,0,234,31,0,31,31,0,31);
  TH3S *fAptMinu=new TH3S("AptMinu","AptMinu",234,0,234,31,0,31,31,0,31);
  TH3S *fAptUnde=new TH3S("AptUnde","AptUnde",234,0,234,31,0,31,31,0,31);
  
  Float_t lptTreshold=0.75;
  Float_t hptTreshold=1.75;
  
  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");  
  AliMUONTriggerCircuit* triggerCircuit;

  for (Int_t icirc=0; icirc<234; icirc++) {
    cout << " Loading LuT for circuit " << icirc << " of 234 " << "\n";
    triggerCircuit = &(pMUON->TriggerCircuit(icirc));  	      

    for (Int_t istripX=0; istripX<31; istripX++) {
      for (Int_t idev=0; idev<31; idev++) {
	
	Short_t iLptPlus, iLptMinu, iLptUnde;
	Short_t iHptPlus, iHptMinu, iHptUnde;
	Short_t iAptPlus, iAptMinu, iAptUnde;
	iLptPlus = iLptMinu = iLptUnde = 0;
	iHptPlus = iHptMinu = iHptUnde = 0;
	iAptPlus = iAptMinu = iAptUnde = 0;
	
	for (Int_t istripY=0; istripY<16; istripY++) {
	  Float_t pt=triggerCircuit->PtCal(istripX,idev,istripY);
	  
	  if (pt>lptTreshold) {
	    if (idev<15)       iLptMinu=iLptMinu+Int_t(TMath::Power(2,istripY));
	    else if (idev==15) iLptUnde=iLptUnde+Int_t(TMath::Power(2,istripY));
	    else if (idev>15)  iLptPlus=iLptPlus+Int_t(TMath::Power(2,istripY));
	  }
	  if (pt>hptTreshold) {
	    if (idev<15)       iHptMinu=iHptMinu+Int_t(TMath::Power(2,istripY));
	    else if (idev==15) iHptUnde=iHptUnde+Int_t(TMath::Power(2,istripY));
	    else if (idev>15)  iHptPlus=iHptPlus+Int_t(TMath::Power(2,istripY));
	  }
	  if (idev<15) 	     iAptMinu=iAptMinu+Int_t(TMath::Power(2,istripY));
	  else if (idev==15) iAptUnde=iAptUnde+Int_t(TMath::Power(2,istripY));
	  else if (idev>15)  iAptPlus=iAptPlus+Int_t(TMath::Power(2,istripY));

	} // loop on istripY

	Int_t bin; 
	
	bin = fLptMinu->GetBin(icirc,istripX,idev);
	fLptMinu->SetBinContent(bin,iLptMinu);
	bin = fLptUnde->GetBin(icirc,istripX,idev);
	fLptUnde->SetBinContent(bin,iLptUnde);
	bin = fLptPlus->GetBin(icirc,istripX,idev);
	fLptPlus->SetBinContent(bin,iLptPlus);

	bin = fHptMinu->GetBin(icirc,istripX,idev);
	fHptMinu->SetBinContent(bin,iHptMinu);
	bin = fHptUnde->GetBin(icirc,istripX,idev);
	fHptUnde->SetBinContent(bin,iHptUnde);
	bin = fHptPlus->GetBin(icirc,istripX,idev);
	fHptPlus->SetBinContent(bin,iHptPlus);
	
	bin = fAptMinu->GetBin(icirc,istripX,idev);
	fAptMinu->SetBinContent(bin,iAptMinu);
	bin = fAptUnde->GetBin(icirc,istripX,idev);
	fAptUnde->SetBinContent(bin,iAptUnde);
	bin = fAptPlus->GetBin(icirc,istripX,idev);
	fAptPlus->SetBinContent(bin,iAptPlus);
	  
      } // loop on idev
    } // loop on istripX
  } // loop on circuit

  hfile->Write();
  hfile->Close();
}







