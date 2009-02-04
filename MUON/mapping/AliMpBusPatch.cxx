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
// $MpId: AliMpBusPatch.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $

//-----------------------------------------------------------------------------
// Class AliMpBusPatch
// --------------------
// The class defines the properties of BusPatch
// Author: Ivana Hrivnacova, IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpBusPatch.h"

#include "AliDAQ.h"
#include "AliMpConstants.h"
#include "AliMpDEManager.h"
#include "AliMpSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpPCB.h"
#include "AliMpMotifPosition.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpBusPatch)
/// \endcond

const Int_t  AliMpBusPatch::fgkOffset = 100;
//
// static methods
//

//____________________________________________________________________
Int_t AliMpBusPatch::GetGlobalBusID(Int_t localID, Int_t ddlID)
{
  /// return global bus id from local bus and ddl id

  return ddlID*fgkOffset + localID;

}
//____________________________________________________________________
Int_t AliMpBusPatch::GetLocalBusID(Int_t globalID, Int_t ddlID)
{
  /// return local bus id from local bus id

  return globalID - ddlID*fgkOffset;

}

//______________________________________________________________________________
AliMpBusPatch::AliMpBusPatch(Int_t id, Int_t detElemId, Int_t ddlId)
  : TObject(),
    fId(id),
    fDEId(detElemId),
    fDdlId(ddlId),
    fManus(false),
    fNofManusPerModule(false),
    fCableLength(-1),
    fCableLabel(),
    fTranslatorLabel(),
    fFrtId(0)
{
/// Standard constructor
}

//______________________________________________________________________________
AliMpBusPatch::AliMpBusPatch(TRootIOCtor* /*ioCtor*/)
  : TObject(),
    fId(),
    fDEId(),
    fDdlId(),
    fManus(false),
    fNofManusPerModule(false),
    fCableLength(-1),
    fCableLabel(),
    fTranslatorLabel(),
    fFrtId(0)
{
/// Root IO constructor
}

//______________________________________________________________________________
AliMpBusPatch::~AliMpBusPatch()
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
Bool_t AliMpBusPatch::AddManu(Int_t manuId)
{
/// Add detection element with given detElemId.
/// Return true if the detection element was added

  if ( HasManu(manuId) ) {
    AliWarningStream() 
      << "Manu with manuId=" << manuId << " already present."
      << endl;
    return false;
  }    

  fManus.Add(manuId);
  return true;
}   

//______________________________________________________________________________
Bool_t AliMpBusPatch::SetNofManusPerModule(Int_t manuNumber)
{
/// Set the number of manus per patch module (PCB):
/// - for stations 1 all manus are connected to one PCB,
/// - for stations 2 there maximum two PCBs per buspatch,
/// - for slat stations there are maximum three PCBs per buspatch

  if ( AliMpDEManager::GetStation12Type(fDEId) == AliMq::kStation1) {

    // simply fill the number of manus, no bridge for station 1
       
    fNofManusPerModule.Add(GetNofManus());
    return true;
  }

 if ( AliMpDEManager::GetStation12Type(fDEId) == AliMq::kStation2) {

    // there is max two patch modules per buspatch
       
    fNofManusPerModule.Add(manuNumber);
    if (manuNumber != GetNofManus())
	fNofManusPerModule.Add(GetNofManus() - manuNumber);

    return true;
  }

  if ( AliMpDEManager::GetStationType(fDEId) == AliMp::kStation345 ) {
  
    const AliMpSlat* kSlat0 
	= AliMpSegmentation::Instance()->GetSlat(fDEId, AliMp::kCath0);

    const AliMpSlat* kSlat1 
	= AliMpSegmentation::Instance()->GetSlat(fDEId, AliMp::kCath1);
       
    Int_t iPcb = 0;
    Int_t iPcbPrev = -1;
    Int_t manuPerPcb = 0;

    Double_t x = 0.;
    Double_t length = 0.;

    // Loop over manu
    for (Int_t iManu = 0; iManu < GetNofManus(); ++iManu) {
      Int_t manuId = GetManuId(iManu);
      AliMpMotifPosition* motifPos0 = kSlat0->FindMotifPosition(manuId);
      AliMpMotifPosition* motifPos1 = kSlat1->FindMotifPosition(manuId);	  
      
      if ( !motifPos0 && !motifPos1 ) {
        // should never happen
        AliErrorStream() 
          << "Motif position for manuId = " << manuId << "not found" << endl;
        return false;
      }

      // find PCB id
      if ( motifPos0 ) {
        x = motifPos0->Position().X();
        length = kSlat0->GetPCB(0)->DX()*2.;
      }
      if ( motifPos1 ) {
        x = motifPos1->Position().X();
        length = kSlat1->GetPCB(0)->DX()*2.;
      }
      
      iPcb = Int_t(x/length + AliMpConstants::LengthTolerance());

      // check when going to next PCB
      if ( iPcb == iPcbPrev )
        manuPerPcb++;
      else if ( iPcbPrev != -1 ) {
        //vec.Set(vec.GetSize()+1);
        //vec[vec.GetSize()-1] = manuPerPcb+1;
        fNofManusPerModule.Add(manuPerPcb+1);
        manuPerPcb = 0;
      }
      iPcbPrev = iPcb;
    }
   
    // store last PCB
    //vec.Set(vec.GetSize()+1);
    //vec[vec.GetSize()-1] = manuPerPcb+1;
    fNofManusPerModule.Add(manuPerPcb+1);
    return true;  
  } 
  
  return false; 
}     

//______________________________________________________________________________
void  AliMpBusPatch::RevertReadout()
{
/// Revert order of manus

  fManus.Revert();
}

//______________________________________________________________________________
void  AliMpBusPatch::ResetReadout()
{
/// Revert order of manus

  fManus.Reset();
}

//______________________________________________________________________________
Int_t AliMpBusPatch::GetNofManus() const
{  
/// Return the number of detection elements connected to this DDL

  return fManus.GetSize(); 
}

//______________________________________________________________________________
Int_t  AliMpBusPatch::GetManuId(Int_t index) const
{  
/// Return the detection element by index (in loop)

  return fManus.GetValue(index); 
}

//______________________________________________________________________________
Bool_t  AliMpBusPatch::HasManu(Int_t manuId) const
{  
/// Return true if bus patch has manu with given manuId

  return fManus.HasValue(manuId); 
}

//______________________________________________________________________________
Int_t  AliMpBusPatch::GetNofPatchModules() const
{
/// Return the number of patch modules (PCB) connected to this bus patch.

  return fNofManusPerModule.GetSize();
}  
  
//______________________________________________________________________________
TString
AliMpBusPatch::GetFRTPosition() const
{
  /// Return CRXX-Y-Z where XX is the Crocus number, Y the FRT number
  /// and Z the local bus patch number.
  return Form("CR%2d-%d-%d",fDdlId,fFrtId+1,GetLocalBusID(fId,fDdlId));
}

//______________________________________________________________________________
Int_t  AliMpBusPatch::GetNofManusPerModule(Int_t patchModule) const
{
/// Return the number of manus per patch module (PCB)

  if ( patchModule < 0 || patchModule >= GetNofPatchModules() ) {
    AliErrorStream() << "Invalid patch module number = " << patchModule << endl;
    return 0;
  }
  
  return fNofManusPerModule.GetValue(patchModule);
}     

//______________________________________________________________________________
void 
AliMpBusPatch::Print(Option_t* opt) const
{
  /// Printout
  
  cout << Form("BusPatch %04d DDL %d : %s <> %s / %s",
               fId,
               AliDAQ::DdlID("MUONTRK",fDdlId),
               GetFRTPosition().Data(),
               fCableLabel.Data(),
               fTranslatorLabel.Data()) << endl;

  TString sopt(opt);
  sopt.ToUpper();
  
  if ( sopt.Contains("FULL") ) 
  {
    cout << Form("Nof of PCBs (i.e. patch modules) = %d",fNofManusPerModule.GetSize()) << endl;
    
    for ( Int_t i = 0; i < fNofManusPerModule.GetSize(); ++i ) 
    {
      cout << Form("\t\t %d manus in patch module %d",fNofManusPerModule.GetValue(i),i) << endl;
    }
    
    if ( sopt.Contains("MANU") )
    {
      cout << "Manus of that buspatch=" << endl;
      
      for ( Int_t i = 0; i < fManus.GetSize(); ++i ) 
      {
        cout << Form("%4d,",fManus.GetValue(i));
      }
      cout << endl;
    }
  }
  
//  Int_t        fId;     ///< Identifier (unique)
//  Int_t        fDEId;   ///< Detection element to which this bus patch is connected
//  Int_t        fDdlId;  ///< DDL to which this bus patch is connected
//  AliMpArrayI  fManus;  ///< Manu Ids connected to this bus patch
//  AliMpArrayI  fNofManusPerModule; ///< Nof Manus per patch modules (PCBs)
//  Float_t      fCableLength;       ///< length of the buspatch cable
//  TString      fCableLabel;        ///< label of the buspatch cable
//  TString      fTranslatorLabel;   ///< label of the translator board
//  Int_t        fFrtId;               ///< FRT Ids connected to this bus patch
  
}
