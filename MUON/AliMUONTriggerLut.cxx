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

#include "AliMUONTriggerLut.h"

#include "AliLog.h"

#include "TFile.h"
#include "TH3.h"

ClassImp(AliMUONTriggerLut)

//----------------------------------------------------------------------
AliMUONTriggerLut::AliMUONTriggerLut() 
  : TNamed()
{
// constructor
  fLptPlus = fLptMinu = fLptUnde = 0;
  fHptPlus = fHptMinu = fHptUnde = 0;
  fAptPlus = fAptMinu = fAptUnde = 0;
}
//----------------------------------------------------------------------
AliMUONTriggerLut::~AliMUONTriggerLut() 
{
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
}

//----------------------------------------------------------------------
AliMUONTriggerLut::AliMUONTriggerLut (const AliMUONTriggerLut& theMUONTriggerLut)
  : TNamed(theMUONTriggerLut)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//----------------------------------------------------------------------
AliMUONTriggerLut & 
AliMUONTriggerLut::operator=(const AliMUONTriggerLut& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal( "Not implemented.");
    
  return *this;  
}

void
AliMUONTriggerLut::ReadFromFile(const char* filename)
{
  TFile f(filename);
  
  if ( f.IsZombie() )
  {
    AliFatal(Form("Could not open file %s",filename));
  }
  
  fLptPlus = (TH3*)(f.Get("LptPlus")->Clone());  
  fLptMinu = (TH3*)(f.Get("LptMinu")->Clone());
  fLptUnde = (TH3*)(f.Get("LptUnde")->Clone());
  fHptPlus = (TH3*)(f.Get("HptPlus")->Clone());  
  fHptMinu = (TH3*)(f.Get("HptMinu")->Clone());
  fHptUnde = (TH3*)(f.Get("HptUnde")->Clone());
  fAptPlus = (TH3*)(f.Get("AptPlus")->Clone());  
  fAptMinu = (TH3*)(f.Get("AptMinu")->Clone());
  fAptUnde = (TH3*)(f.Get("AptUnde")->Clone());

  // insure we "detach" those histograms from file f
  fLptPlus->SetDirectory(0);
  fLptMinu->SetDirectory(0);
  fLptUnde->SetDirectory(0);
  fHptPlus->SetDirectory(0);
  fHptMinu->SetDirectory(0);
  fHptUnde->SetDirectory(0);
  fAptPlus->SetDirectory(0);
  fAptMinu->SetDirectory(0);
  fAptUnde->SetDirectory(0);
}

//----------------------------------------------------------------------
void AliMUONTriggerLut::GetLutOutput(Int_t circuit, Int_t xstrip, Int_t idev,
				     Int_t ystrip, Int_t lutLpt[2], 
				     Int_t lutHpt[2], Int_t lutApt[2])
{
  // return output of LuT for corresponding TH3S  

  if ( !fLptPlus )
  {
    ReadFromFile("$(ALICE_ROOT)/MUON/data/MUONTriggerLut.root");
  }
  
  Int_t bin;
  Short_t binc; 
  Int_t mask = GetMask(ystrip);        // get ystrip mask
  
  // Low pt.............................................. 
  bin    =          fLptPlus->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)fLptPlus->GetBinContent(bin);
  if ((binc & mask)!=0) lutLpt[1]=1;

  bin    =          fLptMinu->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)fLptMinu->GetBinContent(bin);
  if ((binc & mask)!=0) lutLpt[0]=1;
  
  bin    =          fLptUnde->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)fLptUnde->GetBinContent(bin);
  if ((binc & mask)!=0) lutLpt[0]=lutLpt[1]=1;

  // High pt.............................................
  bin    =          fHptPlus->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)fHptPlus->GetBinContent(bin);
  if ((binc & mask)!=0) lutHpt[1]=1;

  bin    =          fHptMinu->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)fHptMinu->GetBinContent(bin);
  if ((binc & mask)!=0) lutHpt[0]=1;

  bin    =          fHptUnde->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)fHptUnde->GetBinContent(bin);
  if ((binc & mask)!=0) lutHpt[0]=lutHpt[1]=1;

  // All pts.............................................
  bin    =          fAptPlus->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)fAptPlus->GetBinContent(bin);
  if ((binc & mask)!=0) lutApt[1]=1;

  bin    =          fAptMinu->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)fAptMinu->GetBinContent(bin);
  if ((binc & mask)!=0) lutApt[0]=1;

  bin    =          fAptUnde->GetBin(circuit,xstrip,idev);
  binc   = (Short_t)fAptUnde->GetBinContent(bin);
  if ((binc & mask)!=0) lutApt[0]=lutApt[1]=1;

}

//----------------------------------------------------------------------
Int_t AliMUONTriggerLut::GetMask(Int_t ystrip)
{
  // returns the mask corresponding to ystrip
  Int_t tabMask[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t mask=0;
  tabMask[ystrip]=1;
  for (Int_t i=0; i<16; i++) 
  {          
    mask += tabMask[i]<<i; 
  }
  return mask;
}





