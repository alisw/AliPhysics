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
// $MpId: AliMpDDL.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $
// Category: management

//-----------------------------------------------------------------------------
// Class AliMpDDL
// --------------------
// The class defines electronics properties of DDL
// Authors: Ivana Hrivnacova, IPN Orsay
//          Christian Finck, SUBATECH Nantes
//-----------------------------------------------------------------------------

#include "AliMpDDL.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"

#include "AliLog.h"

#include <Riostream.h>


using std::endl;
/// \cond CLASSIMP
ClassImp(AliMpDDL)
/// \endcond

//______________________________________________________________________________
AliMpDDL::AliMpDDL(Int_t id)
  : TObject(),
    fId(id),
    fDEIds(),
    fFrtIds(false),
    fBusPatchIds(),
    fTriggerCrateIds(false)

{
/// Standard constructor
}

//______________________________________________________________________________
AliMpDDL::AliMpDDL(TRootIOCtor* /*ioCtor*/)
  : TObject(),
    fId(0),
    fDEIds(),
    fFrtIds(false),    
    fBusPatchIds(),
    fTriggerCrateIds()
{
/// Root IO constructor
}

//______________________________________________________________________________
AliMpDDL::~AliMpDDL()
{
/// Destructor
}

//
// private methods
//

//______________________________________________________________________________
void AliMpDDL::FillBusPatchIds()
{
/// Fill array with bus patch Ids

  for ( Int_t i=0; i<GetNofDEs(); i++ ) {
    AliMpDetElement* detElement 
      = AliMpDEManager::GetDetElement(GetDEId(i));
    
    for ( Int_t j=0; j<detElement->GetNofBusPatches(); j++ )
      fBusPatchIds.Add(detElement->GetBusPatchId(j));
  }
}      

//
// public methods
//

//______________________________________________________________________________
Bool_t AliMpDDL::AddDE(Int_t detElemId)
{
/// Add detection element with given detElemId.
/// Return true if the detection element was added

  if ( ! AliMpDEManager::IsValidDetElemId(detElemId) ) return false;
 
  if ( HasDEId(detElemId) ) {
    AliWarningStream() 
      << "Detection element Id = " << detElemId << " already present."
      << endl;
    return false;
  }    

  AliDebugStream(3) << "Adding detElemId " << detElemId << endl;

  fDEIds.Add(detElemId);
  return true;
}   

//______________________________________________________________________________
Bool_t AliMpDDL::AddTriggerCrate(Int_t crateId)
{
/// Add trigger crate with given crateId.
/// Return true if the trigger crate was added

  if ( HasTriggerCrateId(crateId) ) {
    AliWarningStream() 
	<< "Trigger crate Id = " << crateId << " already present."
	<< endl;
    return false;
  }    
  
  fTriggerCrateIds.Add(crateId);

  return true;
}      

//______________________________________________________________________________
Bool_t AliMpDDL::AddFrt(Int_t frtId)
{
/// Add FRT with given frtId.
/// Return true if the FRT was added

  if ( HasFrtId(frtId) ) {
    AliWarningStream() 
	<< "FRT Id = " << frtId << " already present."
	<< endl;
    return false;
  }    
  
  fFrtIds.Add(frtId);

  return true;
}      


//______________________________________________________________________________
Int_t AliMpDDL::GetNofDEs() const
{  
/// Return the number of detection elements connected to this DDL

  return fDEIds.GetSize(); 
}

//______________________________________________________________________________
Int_t  AliMpDDL::GetDEId(Int_t index) const
{  
/// Return the detection element by index (in loop)

  return fDEIds.GetValue(index); 
}

//______________________________________________________________________________
Bool_t  AliMpDDL::HasDEId(Int_t detElemId) const
{  
/// Return true if the detection element Id is present

  return fDEIds.HasValue(detElemId);; 
}

//______________________________________________________________________________
Int_t AliMpDDL::GetNofFrts() const
{  
/// Return the number of FRT connected to this DDL

  return fFrtIds.GetSize(); 
}

//______________________________________________________________________________
Int_t  AliMpDDL::GetFrtId(Int_t index) const
{  
/// Return the FRT by index (in loop)

  return fFrtIds.GetValue(index); 
}

//______________________________________________________________________________
Bool_t  AliMpDDL::HasFrtId(Int_t frtId) const
{  
/// Return true if the FRT Id is present

  return fFrtIds.HasValue(frtId);; 
}


//______________________________________________________________________________
Int_t AliMpDDL::GetNofBusPatches() const
{  
/// Return the number of detection elements connected to this DDL

  return fBusPatchIds.GetSize(); 
}

//______________________________________________________________________________
Int_t  AliMpDDL::GetBusPatchId(Int_t index) const
{  
/// Return the detection element by index (in loop)

  return fBusPatchIds.GetValue(index); 
}

//______________________________________________________________________________
Bool_t  AliMpDDL::HasBusPatchId(Int_t busPatchId) const
{  
/// Return true if the detection element Id is present

  return fBusPatchIds.HasValue(busPatchId);; 
}

//______________________________________________________________________________
Int_t AliMpDDL::GetNofTriggerCrates() const
{  
/// Return the number of trigger crate connected to this DDL

  return fTriggerCrateIds.GetSize(); 
}

//______________________________________________________________________________
Int_t  AliMpDDL::GetTriggerCrateId(Int_t index) const
{  
/// Return the trigger crate by index (in loop)

  return fTriggerCrateIds.GetValue(index); 
}

//______________________________________________________________________________
Bool_t  AliMpDDL::HasTriggerCrateId(Int_t triggerCrateId) const
{  
/// Return true if the trigger crate Id is present

  return fTriggerCrateIds.HasValue(triggerCrateId);
}

//____________________________________________________________________
Int_t AliMpDDL::GetMaxDsp() const
{
/// calculates the number of DSP 

  Int_t iBusPerBlk = fBusPatchIds.GetSize()/2; //per block

  Int_t iDspMax =  iBusPerBlk/5; //number max of DSP per block
  if (iBusPerBlk % 5 != 0)
    iDspMax += 1;

  return iDspMax;
}

//____________________________________________________________________
void AliMpDDL::GetBusPerDsp(Int_t* iBusPerDSP) 
const
{
/// calculates buspatch per block

  Int_t iBusPerBlk = fBusPatchIds.GetSize()/2; //per block

  for (Int_t i = 0; i < GetMaxDsp(); i++) {
    if ((iBusPerBlk -= 5) > 0) 
      iBusPerDSP[i] = 5;
    else 
      iBusPerDSP[i] = iBusPerBlk + 5;
  }
}

