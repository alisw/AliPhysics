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

/// -----------------------
/// Class AliMUONTriggerLut
/// -----------------------
/// Local Trigger Look Up Table
/// reading interface LUT data is stored into TH3S histograms and readout 
/// from the Local Trigger algorithm
/// Author: Philippe Crochet

#include "AliMUONTriggerLut.h"

#include "AliLog.h"

#include "TFile.h"
#include "TH3.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerLut)
/// \endcond

//----------------------------------------------------------------------
AliMUONTriggerLut::AliMUONTriggerLut() 
    : TNamed(),
      fLptPlus(0),
      fLptMinu(0),
      fLptUnde(0),
      fHptPlus(0),
      fHptMinu(0),
      fHptUnde(0),
      fAptPlus(0),
      fAptMinu(0),
      fAptUnde(0)
{
    //ctor
}

//----------------------------------------------------------------------
AliMUONTriggerLut::~AliMUONTriggerLut() 
{
/// Destructor
  
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

void
AliMUONTriggerLut::ReadFromFile(const char* filename)
{
/// Return output of LuT for corresponding TH3S  

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
				     Int_t lutHpt[2])
{
/// Return output of LuT for corresponding TH3S  

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
/*
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
*/
}

//----------------------------------------------------------------------
Int_t AliMUONTriggerLut::GetMask(Int_t ystrip)
{
/// Return the mask corresponding to ystrip

  Int_t tabMask[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t mask=0;
  tabMask[ystrip]=1;
  for (Int_t i=0; i<16; i++) 
  {          
    mask += tabMask[i]<<i; 
  }
  return mask;
}





