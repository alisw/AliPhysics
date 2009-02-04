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

//-----------------------------------------------------------------------------
/// \class AliMUONTriggerLut
/// 
/// Local Trigger Look Up Table
/// reading interface LUT data is stored into TH3S histograms and readout 
/// from the Local Trigger algorithm
///
/// Histograms structure is :
/// X 234 bins, 1 to 235   = local board number
/// Y  31 bins, 0 to  31   = x strip
/// Z  31 bins, 0 to  31   = x deviation
/// content = Short_t      = y strip mask
///
///  overflow bin is used !
///
/// \author Philippe Crochet
//-----------------------------------------------------------------------------

#include "AliMUONTriggerLut.h"

#include "AliLog.h"

#include <TFile.h>
#include <TH3.h>
#include <TMap.h>
#include <TObjString.h>

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
      fAptUnde(0),
  fMap(0x0)
{
/// ctor
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
  delete fMap;
}

//----------------------------------------------------------------------
Int_t
AliMUONTriggerLut::Compare(TH3* h1, TH3* h2) const
{
/// Return 0 if both histograms are strictly equal (at the bin-by-bin level)

  AliDebug(1,Form("h1 %s h2 %s",h1 ? h1->GetName() : "null", h2 ? h2->GetName() : "null"));

  if (!h1 || !h2) 
  {
    return 0;
  }
  Int_t bin;
  for ( Int_t ix = 0; ix <= h1->GetNbinsX()+1; ix++ ) 
    for ( Int_t iy = 0; iy <= h1->GetNbinsY()+1; iy++ ) 
      for ( Int_t iz = 0; iz <= h1->GetNbinsZ()+1; iz++ ) 
	{
	  {
	    {
	      bin = h1->GetBin(ix,iy,iz);
	      Double_t x1 = h1->GetBinContent(bin);
	      Double_t x2 = h2->GetBinContent(bin);
	      if ( x1 != x2 ) return 0;
	    }
	  }
	}

  AliDebug(1,"same");
  
  return 1;
}

//----------------------------------------------------------------------
Int_t 
AliMUONTriggerLut::Compare(const TObject* object) const
{
/// Return 0 if the two luts are strictly equal

  const AliMUONTriggerLut* lut = static_cast<const AliMUONTriggerLut*>(object);
  
  Int_t rvLpt(0);
  
  rvLpt += Compare(fLptPlus,lut->fLptPlus);
  rvLpt += Compare(fLptMinu,lut->fLptMinu);
  rvLpt += Compare(fLptUnde,lut->fLptUnde);

  Int_t rvHpt(0);
  
  rvHpt += Compare(fHptPlus,lut->fHptPlus);
  rvHpt += Compare(fHptMinu,lut->fHptMinu);
  rvHpt += Compare(fHptUnde,lut->fHptUnde);
  
  Int_t rv(0);
  
  rv += Compare(fAptPlus,lut->fAptPlus);
  rv += Compare(fAptMinu,lut->fAptMinu);
  rv += Compare(fAptUnde,lut->fAptUnde);
  
  AliDebug(1,Form("Same Lpt %d Hpt %d Apt %d",rvLpt,rvHpt,rv));
  
  if ( rvLpt == 3 && rvHpt == 3 ) 
  {
    return 0;
  }

  return 1;
}

//----------------------------------------------------------------------
void AliMUONTriggerLut::ReadFromFile(const char* filename)
{
/// Return output of LuT for corresponding TH3S  

  TFile f(filename);
  
  if ( f.IsZombie() )
  {
    AliFatal(Form("Could not open file %s",filename));
  }
  
  AliDebug(1,Form("filename=%s",filename));
  
  fLptPlus = (TH3*)(f.Get("LptPlus"));  
  fLptMinu = (TH3*)(f.Get("LptMinu"));
  fLptUnde = (TH3*)(f.Get("LptUnde"));
  fHptPlus = (TH3*)(f.Get("HptPlus"));  
  fHptMinu = (TH3*)(f.Get("HptMinu"));
  fHptUnde = (TH3*)(f.Get("HptUnde"));
  fAptPlus = (TH3*)(f.Get("AptPlus"));  
  fAptMinu = (TH3*)(f.Get("AptMinu"));
  fAptUnde = (TH3*)(f.Get("AptUnde"));
  
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
  
  RegisterHistos();
}

//----------------------------------------------------------------------
void
AliMUONTriggerLut::RegisterHistos()
{
/// Add histos to our internal map

  Add(fLptPlus);
  Add(fLptMinu);
  Add(fLptUnde);

  Add(fHptPlus);
  Add(fHptMinu);
  Add(fHptUnde);

  Add(fAptPlus);
  Add(fAptMinu);
  Add(fAptUnde);
}

//----------------------------------------------------------------------
void
AliMUONTriggerLut::Add(TH3* h)
{
  /// Update internal map
  if (!fMap)
  {
    fMap = new TMap;
    fMap->SetOwner(kTRUE);
  }
  
  if (h) fMap->Add(new TObjString(h->GetName()),h);
}

//----------------------------------------------------------------------
void 
AliMUONTriggerLut::SetContent(const char* hname, Int_t icirc, UChar_t istripX, 
                              UChar_t idev, Short_t value)
{
  /// Set the content of one bin of one histogram
  
  if (!fMap)
  {
    //..........................................circuit/stripX/deviation
    fLptPlus = new TH3S("LptPlus","LptPlus",234,0,234,31,0,31,31,0,31);
    fLptMinu = new TH3S("LptMinu","LptMinu",234,0,234,31,0,31,31,0,31);
    fLptUnde = new TH3S("LptUnde","LptUnde",234,0,234,31,0,31,31,0,31);
    
    fHptPlus = new TH3S("HptPlus","HptPlus",234,0,234,31,0,31,31,0,31);
    fHptMinu = new TH3S("HptMinu","HptMinu",234,0,234,31,0,31,31,0,31);
    fHptUnde = new TH3S("HptUnde","HptUnde",234,0,234,31,0,31,31,0,31);
    
    fAptPlus = new TH3S("AptPlus","AptPlus",234,0,234,31,0,31,31,0,31);
    fAptMinu = new TH3S("AptMinu","AptMinu",234,0,234,31,0,31,31,0,31);
    fAptUnde = new TH3S("AptUnde","AptUnde",234,0,234,31,0,31,31,0,31);
    
    RegisterHistos();
  }
  
  TH3* h = static_cast<TH3*>(fMap->GetValue(hname));
  
  Int_t bin = h->GetBin(icirc,istripX,idev);
  h->SetBinContent(bin,value);
}

//----------------------------------------------------------------------
void AliMUONTriggerLut::GetLutOutput(Int_t circuit, Int_t xstrip, Int_t idev,
                                     Int_t ystrip, Int_t lutLpt[2], 
                                     Int_t lutHpt[2]) const
{
/// Return output of LuT for corresponding TH3S  

  if ( !fLptPlus )
  {
    AliError("LUT not initialized");
//    ReadFromFile("$(ALICE_ROOT)/MUON/data/MUONTriggerLut.root");
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
Int_t AliMUONTriggerLut::GetMask(Int_t ystrip) const
{
/// Return the mask corresponding to ystrip

    return (Int_t)(1<<ystrip);
}





