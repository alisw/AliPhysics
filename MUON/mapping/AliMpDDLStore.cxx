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
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEStore.h"
#include "AliMpDDL.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpBusPatch.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TList.h>

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
  fBusPatches(true),
  fManuList12()
{  
/// Standard constructor

  AliDebug(1,"");
  fDDLs.SetOwner(true);
  fBusPatches.SetOwner(true);
  fBusPatches.SetSize(900);

  // Create all detection elements
  ReadDDLs();
  SetManus();
  SetPatchModules();
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
Int_t  AliMpDDLStore::GetManuListIndex(Int_t detElemId) const
{
/// Return the index of the manu list for given detElemId

  return AliMpDEManager::GetChamberId(detElemId)*4 + (detElemId % 100);
}  
  

//______________________________________________________________________________
Int_t AliMpDDLStore::GetBusPatchIndex(Int_t detElemId, Int_t manuId) const
{
/// Calculate the index of the buspatch from manuId

    Int_t pos = 0;
    AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
    static Int_t manuMask = AliMpConstants::ManuMask(AliMp::kNonBendingPlane) - 1;

    if( stationType == AliMp::kStation345) {
      pos = (manuId & manuMask)/100;
    } else {
      Int_t idx = GetManuListIndex(detElemId);

      // using array defined from DetElemIdToBusPatch.dat file
      for (pos = fManuList12[idx].GetSize()-1; pos >= 0; --pos)
	  if ( manuId >= fManuList12[idx].At(pos)) break;
    }

    return pos;
}


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
      
  char line[255];

  while ( in.getline(line,255) ) {

     if ( line[0] == '#' ) continue;

     TString tmp(AliMpHelper::Normalize(line));

     Int_t blankPos  = tmp.First(' ');
     Int_t blankPos1 = tmp.Last(' ');
     Int_t length = 0;

     TString sDE(tmp(0, blankPos));

     Int_t idDE = atoi(sDE.Data());

     // reading 1st manu Id for each bus patch (station 1 & 2)
     if(AliMpDEManager::GetStationType(idDE) != AliMp::kStation345) {

       TString sManu(tmp(blankPos1 + 1, tmp.Length()-blankPos1));
       AliMpHelper::DecodeName(sManu,',',fManuList12[GetManuListIndex(idDE)]);

       TString tmp1(tmp(blankPos + 1, blankPos1 -  blankPos));
       blankPos1 = blankPos + tmp1.First(' ') + 1;
       length    = tmp.Last(' ') - blankPos1; 

     } else {
       length = tmp.Length()-blankPos1;
     }

     TString sDDL(tmp(blankPos1 + 1, length));
     Int_t  iDDL = atoi(sDDL.Data());
           

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
//______________________________________________________________________________
Bool_t AliMpDDLStore::SetManus()
{
/// Set manus for each bus patch

    Int_t manuMask = AliMpConstants::ManuMask(AliMp::kNonBendingPlane) - 1;

    // loop over DDL 
    for (Int_t iDDL = 0; iDDL < fgkNofDDLs; ++iDDL) {

      AliDebug(3, Form("DDL # %d\n", iDDL));

      AliMpDDL* ddl = GetDDL(iDDL);

      // loop over DE in the given DDL
      for (Int_t detElemIdx = 0; detElemIdx < ddl->GetNofDEs(); ++detElemIdx) {

	Int_t detElemId = ddl->GetDEId(detElemIdx);

	AliMpDetElement* detElement = GetDetElement(detElemId);

	AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);


	// list of manu per DE on both cathode
	TList manuList;
	for ( Int_t cath = 0; cath < 2 ; ++cath ) {
	  const AliMpVSegmentation* seg 
	      = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));
        
	  AliMp::PlaneType planeType = detElement->GetPlaneType(AliMp::GetCathodType(cath));

	  TArrayI manus;
	  seg->GetAllElectronicCardIDs(manus);
          
	  // filling TList manu
	  for ( Int_t im = 0; im < manus.GetSize(); ++im ) {

	    AliMpIntPair* manu = 0x0;
	    if( stationType == AliMp::kStation345) 
		manu = new AliMpIntPair((manus[im] & manuMask), planeType, kTRUE); //remove offset for NB
	    else 
		manu = new AliMpIntPair(manus[im], planeType, kTRUE); //keep offset for NB
	    
	    manuList.Add(manu);
	  }        
	}// cathode

	manuList.Sort(); // sort regardless B or NB plane

	// filling manu to the corresponding buspatch
	for (Int_t iEntry = 0; iEntry < manuList.GetEntries(); ++iEntry) {

	  AliMpIntPair* manuPtr = (AliMpIntPair*)manuList.At(iEntry);

	  Int_t manuId = manuPtr->GetFirst();
	  Int_t pos    = GetBusPatchIndex(detElemId, manuId);

	  if (pos > detElement->GetNofBusPatches()) {
	    AliError(Form("pos greater %d than size %d manuId %d detElemId %d \n", 
			  pos, detElement->GetNofBusPatches(), manuId, detElemId));
	    return false;
	  }
	  
	  // get buspatch and fill manus
	  Int_t busPatchId = detElement->GetBusPatchId(pos);
	  AliMpBusPatch* busPatch = GetBusPatch(busPatchId);

	  if( stationType == AliMp::kStation345) {

	    if (manuPtr->GetSecond()) 
		busPatch->AddManu(manuId+manuMask+1); // add offset again after sorted
	    else 
		busPatch->AddManu(manuId); 

	  } else {
	
	    busPatch->AddManu(manuId); 

	  }
	}

	manuList.Delete();

	if (AliDebugLevel() == 3) {

	  // print out for checking
	  for(Int_t pos = 0; pos < detElement->GetNofBusPatches(); ++pos) {
	    Int_t busPatchId = detElement->GetBusPatchId(pos);
	    AliMpBusPatch* busPatch = GetBusPatch(busPatchId);
	    printf("BusPatch: %d\n", busPatch->GetId());
	    for (Int_t iEntry = 0; iEntry < busPatch->GetNofManus(); ++iEntry) 
		printf("manu Id: %d\n", busPatch->GetManuId(iEntry));
	  }
	}

      } // detection element loop
    }// DDL loop

    return true;
}

//______________________________________________________________________________
Bool_t AliMpDDLStore::SetPatchModules()
{
/// Compute the number of manu per PCB for each buspatch 

  Bool_t result = true;

  for (Int_t i = 0; i < fBusPatches.GetSize(); ++i) {
    AliMpBusPatch* busPatch = (AliMpBusPatch*)fBusPatches.GetObject(i);
    Bool_t newResult = busPatch->SetNofManusPerModule();
    result = result && newResult;

    if (AliDebugLevel() == 3) {
      // print out for checking
      printf("\nbus patch %d\n", busPatch->GetId());
      for (Int_t i = 0; i < busPatch->GetNofPatchModules(); ++i) 
        printf("manu per %dth pcb %d\n", i, busPatch->GetNofManusPerModule(i));
    }    
  }  
  
  return result;
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
Int_t AliMpDDLStore::GetBusPatchId(Int_t detElemId, Int_t manuId) const
{
/// Return bus patch for a given manuId

  AliMpDetElement* detElement = GetDetElement(detElemId);
  Int_t pos = GetBusPatchIndex(detElemId, manuId);

  if ( pos > detElement->GetNofBusPatches() ) {
    AliErrorStream() 
      << "Pos = " << pos 
      << " greater than the size = " <<  detElement->GetNofBusPatches()
      << " for detElemId = " << detElemId 
      << " manuId = " << manuId << endl;
    return -1;
  }

  return detElement->GetBusPatchId(pos);
}    

//______________________________________________________________________________
AliMpIntPair  AliMpDDLStore::GetDetElemIdManu(Int_t manuSerial) const
{
/// Return the detElemId and manuId for given serial manu number

  return fDetElements->GetDetElemIdManu(manuSerial);
}  

//______________________________________________________________________________
void AliMpDDLStore::PrintAllManu() const
{
/// Print all manu Ids and their serial numbers sorted by detection element
/// and bus patch.                                                            \n
/// As serial manu numbers are filled in a different way than manu Ids this
/// printing allows to check that both ways are consistent

  // Loop over DE
  AliMpDEIterator it;
  for ( it.First(); ! it.IsDone(); it.Next() ) {
     AliMpDetElement* de = it.CurrentDE();
     cout << "DE: " << de->GetId() << endl; 
     
     // Loop over bus patches in this DE
     for ( Int_t i=0; i< de->GetNofBusPatches(); ++i ) {

       AliMpBusPatch* busPatch = GetBusPatch(de->GetBusPatchId(i));    
       cout << "  busPatch: " << busPatch->GetId() << endl; 

       cout << "    Manu       : ";
       for ( Int_t j=0; j<busPatch->GetNofManus(); ++j ) {
         cout << std::setw(6) << busPatch->GetManuId(j) << " ";
       }
       cout << endl;
       
       cout << "    Manu serial: ";
       for ( Int_t k=0; k<busPatch->GetNofManus(); ++k ) {
         cout << std::setw(6) << de->GetManuSerialFromId(busPatch->GetManuId(k)) << " ";
       }        
       cout << endl;
     }
  }
}  


