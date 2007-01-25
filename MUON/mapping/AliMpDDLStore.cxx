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
// $MpId: AliMpDDLStore.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $
// Category: management
//
// Class AliMpDDLStore
// --------------------
// The top container class for DDLs, det elements and bus patched
// It provides acces to DDL, det element and bus patches objects
// via various characteristics.
// Authors: Ivana Hrivnacova, IPN Orsay
//          Christian Finck, SUBATECH Nantes

#include "AliMpDDLStore.h"
#include "AliMpDEStore.h"
#include "AliMpDDL.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpBusPatch.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpDDLStore)
/// \endcond

AliMpDDLStore* AliMpDDLStore::fgInstance = 0;
const Int_t    AliMpDDLStore::fgkNofDDLs = 20;

//
// static methods
//

//______________________________________________________________________________
AliMpDDLStore* AliMpDDLStore::Instance()
{
/// Create the DDL store if it does not yet exist
/// and return its instance

  if ( ! fgInstance )
    fgInstance = new AliMpDDLStore();
    
  return fgInstance;
}    

//
// ctors, dtor
//

//______________________________________________________________________________
AliMpDDLStore::AliMpDDLStore()
: TObject(),
  fDDLs(fgkNofDDLs),
  fDetElements(AliMpDEStore::Instance()),
  fBusPatches(true)
{  
/// Standard constructor

  AliDebug(1,"");
  fDDLs.SetOwner(true);
  fBusPatches.SetOwner(true);
  fBusPatches.SetSize(900);

  // Create all detection elements
  ReadDDLs();
}

//______________________________________________________________________________
AliMpDDLStore::AliMpDDLStore(TRootIOCtor* /*ioCtor*/)
: TObject(),
  fDDLs(),
  fDetElements(0),
  fBusPatches()
{  
/// Constructor for IO

  AliDebug(1,"");

  fgInstance = this;
}


//______________________________________________________________________________
AliMpDDLStore::~AliMpDDLStore()
{
/// Destructor

  AliDebug(1,"");

  delete fDetElements;

  // DDL objects are deleted with fDDLs 
  // Bus patches objects are deleted with fBusPatches 
  
  fgInstance = 0;
}

//
// private methods
//

//______________________________________________________________________________
Bool_t AliMpDDLStore::ReadDDLs()
{
/// Read ddl <-> bus patch file
 
  TString infile = AliMpFiles::BusPatchFilePath();

  ifstream in(infile, ios::in);
  if (!in) {
    AliErrorStream() << "Data file " << infile << " not found.";
    return false;
  }  
      
  char line[80];

  while ( in.getline(line,80) ) {

     if ( line[0] == '#' ) continue;

     TString tmp(AliMpHelper::Normalize(line));

     Int_t blankPos  = tmp.First(' ');
     Int_t blankPos1 = tmp.Last(' ');

     TString sDE(tmp(0, blankPos));

     Int_t idDE = atoi(sDE.Data());

     TString sDDL(tmp(blankPos1 + 1, tmp.Length()-blankPos1));

     Int_t iDDL = atoi(sDDL.Data());

     // always working local DDL number... for the moment.

     // not really needed remove for stand alone purpose (Ch.F)
     //      if (iDDL >= AliDAQ::DdlIDOffset("MUONTRK"))
     //        iDDL -= AliDAQ::DdlIDOffset("MUONTRK");


     TString busPatch(tmp(blankPos + 1,blankPos1-blankPos-1));
     AliDebugStream(3)
       << "idDE " << idDE << " buspatch " << busPatch.Data() << " iDDL " << iDDL 
       << endl;
     
     if ( iDDL < 0 || iDDL >= fgkNofDDLs ) {
       AliErrorStream() << "DDL id "<< iDDL << " outside limits." << endl;
       return false;
     }  

     if ( ! AliMpDEManager::IsValidDetElemId(idDE, false) ) {
       AliErrorStream() << "DetElemId "<< idDE << " not valid." << endl;
       return false;
     }  
         

     AliMpDDL* ddl = GetDDL(iDDL, false);
     if ( !ddl) {
       ddl = new AliMpDDL(iDDL);
       fDDLs.AddAt(ddl, iDDL);
     }  
     ddl->AddDE(idDE);
     

     TArrayI busPatchList;
     // decoding range of buspatch
     AliMpHelper::DecodeName(busPatch,';',busPatchList);
     
     // Update DE
     AliMpDetElement* de = AliMpDEManager::GetDetElement(idDE);
     de->SetDdlId(iDDL);
     // filling buspatch -> idDE
     for (Int_t i = 0; i < busPatchList.GetSize(); i++) {
       fBusPatches.Add(busPatchList[i], 
                       new AliMpBusPatch(busPatchList[i], idDE, iDDL));
       de->AddBusPatch(busPatchList[i]);
     }  
   }
   
   // Fill bus patch Ids array in DDLs now
   for ( Int_t i=0; i<fDDLs.GetEntriesFast(); i++ ) {
     AliMpDDL* ddl = (AliMpDDL*) fDDLs.At(i);
     ddl->FillBusPatchIds();
   }  
   
   in.close();
   return true;
}

//
// public methods
//


//______________________________________________________________________________
AliMpDDL* AliMpDDLStore::GetDDL(Int_t ddlId, Bool_t warn) const
{
/// Return DDL for given ddlId

  AliMpDDL* ddl 
    = (AliMpDDL*)fDDLs.At(ddlId);
    
  if ( ! ddl && warn ) {  
    AliErrorStream() 
        << "DDL with Id = " << ddlId << " not defined." << endl;
  }	

  return ddl;
}    

//______________________________________________________________________________
AliMpDetElement*  AliMpDDLStore::GetDetElement(Int_t detElemId, Bool_t warn) const
{
/// Return detection element with given detElemId

  return fDetElements->GetDetElement(detElemId, warn);
}  

//______________________________________________________________________________
AliMpBusPatch* AliMpDDLStore::GetBusPatch(Int_t busPatchId, Bool_t warn) const
{
/// Return bus patch with given Id

  AliMpBusPatch* busPatch
    = (AliMpBusPatch*) fBusPatches.GetValue(busPatchId);
    
  if ( ! busPatch && warn ) {  
    AliErrorStream() 
        << "Bus patch with Id = " << busPatchId << " not defined." << endl;
  }	

  return busPatch;
}    

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetDEfromBus(Int_t busPatchId) const
{
/// Return detection element Id for given busPatchId

  AliMpBusPatch* busPatch = GetBusPatch(busPatchId);
  
  if ( ! busPatch ) {
    AliErrorStream() 
        << "Bus patch with Id = " << busPatchId << " not defined." << endl;
    return 0;    
  }	

  return busPatch->GetDEId();
}  

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetDDLfromBus(Int_t busPatchId) const
{
/// Return DDL Id for given busPatchId

  AliMpBusPatch* busPatch = GetBusPatch(busPatchId);
  
  if ( ! busPatch ) {
    AliErrorStream() 
        << "Bus patch with Id = " << busPatchId << " not defined." << endl;
    return 0;    
  }	

  return busPatch->GetDdlId();
}  

//______________________________________________________________________________
AliMpIntPair  AliMpDDLStore::GetDetElemIdManu(Int_t manuSerial) const
{
/// Return the detElemId and manuId for given serial manu number

  return fDetElements->GetDetElemIdManu(manuSerial);
}  
