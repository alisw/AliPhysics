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

// $Id$

#include "AliLog.h"
#include "AliMpConstants.h"

#include "TH1F.h"
#include "TFile.h"

#include <fstream>
#include <cassert>

#include "AliMUONTriggerEfficiencyCells.h"


//-----------------------------------------------------------------------------
/// \class AliMUONTriggerEfficiencyCells
/// A class to store and give access to the numerator and denominator 
/// histograms for the trigger chamber efficiency calculation.
///
/// \author Diego Stocco; Subatech, Nantes
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONTriggerEfficiencyCells)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells()
:
TObject(),
fCountHistoList(0x0),
fNoCountHistoList(0x0),
fFiredStrips(0x0)
{
///  Default constructor.
  CheckConstants();
  ResetHistos();
}

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells(const Char_t* filename,
							     const Char_t* listname)
:
TObject(),
fCountHistoList(0x0),
fNoCountHistoList(0x0),
fFiredStrips(0x0)
{
///  Constructor using an ASCII file.
  CheckConstants();
  ResetHistos();
  ReadFile(filename, listname);
}

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells(TList *countHistoList)
:
TObject(),
fCountHistoList(countHistoList),
fNoCountHistoList(0x0),
fFiredStrips(0x0)
{
///  Constructor using a list of histograms with counts.
  CheckConstants();
  ResetHistos();
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells(const AliMUONTriggerEfficiencyCells& other)
:
TObject(other),
fCountHistoList(other.fCountHistoList),
fNoCountHistoList(other.fNoCountHistoList),
fFiredStrips(other.fFiredStrips)
{
/// Copy constructor

  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    fBoardEfficiency[chCath] = other.fBoardEfficiency[chCath];
    fSlatEfficiency[chCath] = other.fSlatEfficiency[chCath];
  }
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells& AliMUONTriggerEfficiencyCells::operator=(const AliMUONTriggerEfficiencyCells& other)
{
  /// Asignment operator
  // check assignement to self
  if (this == &other)
    return *this;

  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    fBoardEfficiency[chCath] = other.fBoardEfficiency[chCath];
    fSlatEfficiency[chCath] = other.fSlatEfficiency[chCath];
  }

  fCountHistoList = other.fCountHistoList;
  fNoCountHistoList = other.fNoCountHistoList;
  fFiredStrips = other.fFiredStrips;
    
  return *this;
}

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::~AliMUONTriggerEfficiencyCells()
{
///  Destructor.
  //delete [] fBoardEfficiency;
  //delete [] fSlatEfficiency;
  ResetHistos(kTRUE);
  delete fCountHistoList;
  delete fNoCountHistoList;
  delete fFiredStrips;
}


//_____________________________________________________________________________
void AliMUONTriggerEfficiencyCells::CheckConstants() const
{
/// Check consistence of redefined constants 

  assert(fgkNcathodes == AliMpConstants::NofCathodes());    
  assert(fgkNchambers == AliMpConstants::NofTriggerChambers());    
  assert(fgkNplanes == AliMpConstants::NofTriggerChambers() * fgkNcathodes);    
}


//__________________________________________________________________________
void
AliMUONTriggerEfficiencyCells::ResetHistos(Bool_t deleteObjects)
{
///  Sets our internal array contents to zero.

  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    if ( deleteObjects ) {
      delete fBoardEfficiency[chCath];
      delete fSlatEfficiency[chCath];
    }
    fBoardEfficiency[chCath] = 0x0;
    fSlatEfficiency[chCath] = 0x0;
  }
}


//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::ReadFile(const Char_t* filename, const Char_t* listname)
{
///  Structure of file (.root) containing local board efficency
    TFile *file = new TFile(filename, "read");
    if(!file || !file->IsOpen()) {
      AliError(Form("Can't read file %s",filename));
      return;
    }

    if ( ! fCountHistoList ) {
      fCountHistoList = new TList();
      fCountHistoList->SetOwner();
    }

    TH1F *histo = 0x0;
    const Char_t* histoName;

    TList* listInFile = 0x0;
    TString listNameString(listname);
    if ( ! listNameString.IsNull() )
      listInFile = (TList*)file->FindObjectAny(listname);

    for ( Int_t ide=0; ide<=kHchamberCount; ide++){
      for(Int_t ich=0; ich<fgkNchambers; ich++){

	// Efficiency per chamber is provided by 1 histogram only
	if ( ide == kHchamberCount ) ich = fgkNchambers;

	for(Int_t hType=0; hType<kNcounts; hType++){
	  histoName = GetHistoName(ide, hType, ich);
	  histo = ( listInFile ) ? (TH1F*)listInFile->FindObject(histoName) : (TH1F*)file->Get(histoName);
	  if ( ! histo ) {
	    AliWarning(Form("Cannot find %s in file. Skip histogram", histoName));
	    continue;
	  }
	  histo->SetDirectory(0);
	  fCountHistoList->Add(histo);

	  AliDebug(5,Form("Adding histogram %s\n",histoName));

	  // Do not fill efficiency per chamber histos
	  if ( ide == kHchamberCount )
	    continue;

	  // Fill old historgrams for consistency
	  if ( hType < kBothPlanesEff ){
	    TString newName = histoName;
	    newName.ReplaceAll("Counts","Eff");
	    TH1F* auxHisto = (TH1F*)histo->Clone(newName.Data());
	    auxHisto->SetDirectory(0);
	    if ( ide == kHboardCount )
	      fBoardEfficiency[fgkNchambers*hType+ich] = auxHisto;
	    else if ( ide == kHslatCount )
	      fSlatEfficiency[fgkNchambers*hType+ich] = auxHisto;

	    AliDebug(5,Form("Creating histogram %s\n",auxHisto->GetName()));
	  }
	  else if ( hType == kAllTracks ){
	    for ( Int_t icath=0; icath<2; icath++){
	      if ( ide == kHboardCount )
		fBoardEfficiency[fgkNchambers*icath+ich]->Divide(histo);
	      else if ( ide == kHslatCount )
		fSlatEfficiency[fgkNchambers*icath+ich]->Divide(histo);
	    }
	  }
	} // loop on count types
      } // loop on chambers
    } // loop on detection element type

    file->Close();
}


//__________________________________________________________________________
const Char_t*
AliMUONTriggerEfficiencyCells::GetHistoName(Int_t histoType, Int_t countType, 
					    Int_t chamber)
{
  //
  /// Return the name of the histogram for efficiency calculation
  //
  TString histoTypeName[kNcounts] = {"bendPlane", "nonBendPlane", "bothPlanes", "allTracks"};

  switch ( histoType ) {
  case kHchamberCount:
    return Form("%sCountChamber", histoTypeName[countType].Data());
  case kHslatCount:
    return Form("%sCountSlatCh%i", histoTypeName[countType].Data(), 11+chamber);
  case kHboardCount:
    return Form("%sCountBoardCh%i", histoTypeName[countType].Data(), 11+chamber);
  }

  return 0x0;
}

//__________________________________________________________________________
TH1F* AliMUONTriggerEfficiencyCells::GetOldEffHisto(Int_t histoType,
						    Int_t ich, Int_t icath) const
{
  //
  /// Compatibility with old class
  /// Gets the efficiency from the array
  /// (which are empty in the new implementation)

  switch ( histoType ) {
  case kHboardCount:
    return fBoardEfficiency[fgkNchambers*icath+ich];
  case kHslatCount:
    return fSlatEfficiency[fgkNchambers*icath+ich];
  }
  return 0x0;
}
