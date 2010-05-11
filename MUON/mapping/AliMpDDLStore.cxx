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

//-----------------------------------------------------------------------------
// Class AliMpDDLStore
// --------------------
// The top container class for DDLs, det elements and bus patched
// It provides acces to DDL, det element and bus patches objects
// via various characteristics.
// Authors: Ivana Hrivnacova, IPN Orsay
//          Christian Finck, SUBATECH Nantes
//-----------------------------------------------------------------------------

#include <cstdlib>
#include "AliMpDDLStore.h"
#include "AliMpExMapIterator.h"
#include "AliMpConstants.h"
#include "AliMpDEStore.h"
#include "AliMpFrtCrocusConstants.h"
#include "AliMpDDL.h"
#include "AliMpFiles.h"
#include "AliMpDataStreams.h"
#include "AliMpHelper.h"
#include "AliMpDEManager.h"
#include "AliMpManuStore.h"
#include "AliMpDetElement.h"
#include "AliMpBusPatch.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpStringObjMap.h"
#include "AliMpEncodePair.h"
#include "AliMpIntPair.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TList.h>
#include <TObjArray.h>
#include <TString.h>
#include <TObjString.h>
#include <TClass.h>

/// \cond CLASSIMP
ClassImp(AliMpDDLStore)
/// \endcond

AliMpDDLStore* AliMpDDLStore::fgInstance = 0;
const Int_t    AliMpDDLStore::fgkNofDDLs = 20;
const Int_t    AliMpDDLStore::fgkNofTriggerDDLs = 2;
const TString  AliMpDDLStore::fgkRevertKeyword = "REVERT"; 
const TString  AliMpDDLStore::fgkExplicitKeyword = "EXPLICIT";

//
// static methods
//

//______________________________________________________________________________
AliMpDDLStore* AliMpDDLStore::Instance(Bool_t warn) 
{
    /// Create the DDL store if it does not yet exist
    /// and return its instance

    if ( ! fgInstance && warn  ) {
        AliWarningClass("DDL Store has not been loaded");
    }

    return fgInstance;
}

//______________________________________________________________________________
AliMpDDLStore* AliMpDDLStore::ReadData(const AliMpDataStreams& dataStreams,
                                       Bool_t warn) 
{
    /// Load the DDL store from ASCII data files
    /// and return its instance

    if ( fgInstance ) {
        if ( warn )
            AliWarningClass("DDL Store has been already loaded");
        return fgInstance;
    }

    if ( dataStreams.GetReadFromFiles() )
      AliInfoClass("Reading DDL Store from ASCII files.");

    fgInstance = new AliMpDDLStore(dataStreams);
    return fgInstance;
}

//
// ctors, dtor
//

//______________________________________________________________________________
AliMpDDLStore::AliMpDDLStore(const AliMpDataStreams& dataStreams)
        : TObject(),
        fkDataStreams(dataStreams),
        fDDLs(fgkNofDDLs+fgkNofTriggerDDLs), // FIXEME
        fBusPatches(),
        fManuList12(),
        fManuBridge2(),
        fRegionalTrigger()
{
  /// Standard constructor
  
  AliDebug(1,"");
  fDDLs.SetOwner(true);
  fBusPatches.SetOwner(true);
  fBusPatches.SetSize(900);

  // Load segmentation & DE store data
  if ( ! AliMpSegmentation::Instance(false) )
    AliMpSegmentation::ReadData(dataStreams, true);
    
  // Create all detection elements
  ReadDDLs();
  ReadTrigger();
  SetTriggerDDLs();
  SetManus();
  ReadBusPatchSpecial();
  SetPatchModules();
  ReadBusPatchInfo();
}

//______________________________________________________________________________
AliMpDDLStore::AliMpDDLStore(TRootIOCtor* ioCtor)
        : TObject(),
        fkDataStreams(ioCtor),
        fDDLs(),
        fBusPatches(ioCtor),
        fRegionalTrigger(ioCtor)
{
    /// Constructor for I0

    AliDebug(1,"");

    fgInstance = this;
}


//______________________________________________________________________________
AliMpDDLStore::~AliMpDDLStore() 
{
    /// Destructor

    AliDebug(1,"");

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
            if ( manuId >= fManuList12[idx].At(pos))
                break;
    }

    return pos;
}

//______________________________________________________________________________
Bool_t AliMpDDLStore::ReadDDLs() 
{
    /// Read ddl <-> bus patch file

    istream& in 
      = fkDataStreams.
          CreateDataStream(AliMpFiles::BusPatchFilePath());

    char line[255];

    while ( in.getline(line,255) ) {

        if ( line[0] == '#' )
            continue;

        TString tmp(AliMpHelper::Normalize(line));

        TObjArray* stringList = tmp.Tokenize(TString(" "));

        TString sDE = ((TObjString*)stringList->At(0))->GetString();
        Int_t idDE  = atoi(sDE.Data());

        if ( ! AliMpDEManager::IsValidDetElemId(idDE, false) ) {
            AliErrorStream() << "DetElemId "<< idDE << " not valid." << endl;
            return false;
        }

        TString busPatch = ((TObjString*)stringList->At(1))->GetString();


        TString sDDL = ((TObjString*)stringList->At(2))->GetString();
        Int_t  iDDL  = atoi(sDDL.Data());

        if ( iDDL < 0 || iDDL >= fgkNofDDLs ) {
            AliErrorStream() << "DDL id "<< iDDL << " outside limits." << endl;
            return false;
        }

        AliDebugStream(3)
        << "idDE " << idDE << " buspatch " << busPatch.Data() << " iDDL " << iDDL
        << endl;

        // reading 1st manu Id for each bus patch (station 1 & 2)
        if(AliMpDEManager::GetStationType(idDE) != AliMp::kStation345) {

            TString sManu = ((TObjString*)stringList->At(3))->GetString();
            AliMpHelper::DecodeName(sManu,',',fManuList12[GetManuListIndex(idDE)]);

            if(AliMpDEManager::GetStation12Type(idDE) == AliMq::kStation2) {
                TString sManuBridge = ((TObjString*)stringList->At(4))->GetString();
                AliMpHelper::DecodeName(sManuBridge,',',fManuBridge2[GetManuListIndex(idDE)]);
            }

        }

        delete stringList;

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

    delete &in;
    return true;
}

//______________________________________________________________________________
Bool_t  AliMpDDLStore::ReadTrigger() 
{
    /// create trigger DDL object and Global crate object
  
  if ( ! fRegionalTrigger.ReadData(fkDataStreams) ) return false;

  return true;
}

//______________________________________________________________________________
Bool_t  
AliMpDDLStore::SetTriggerDDLs() 
{
  /// Create trigger DDLs and set DDL Ids in the regional trigger

  Int_t iDDL = -1;
  TIter next(fRegionalTrigger.CreateCrateIterator());
  AliMpTriggerCrate* crate;
  
  while ( ( crate = static_cast<AliMpTriggerCrate*>(next()) ) )
  {
    TString crateName = crate->GetName();

    // determine ddl number vs crate side
    if (crateName.Contains("R"))
      iDDL = fgkNofDDLs; // starts where tracker ends
    else
      iDDL = fgkNofDDLs + 1;

    // Create DDL if it does not yet exist and set it to the crate
    AliMpDDL* ddl = (AliMpDDL*)fDDLs.At(iDDL);
    if ( !ddl) {
      ddl = new AliMpDDL(iDDL);
      fDDLs.AddAt(ddl, iDDL);
    }
    crate->SetDdlId(iDDL);
    
    
    // Add trigger crate number for given ddl if not present
    if ( !ddl->HasTriggerCrateId(crate->GetId()) )
      ddl->AddTriggerCrate(crate->GetId());
    
    
    // Loop over local boards in this crate

    for ( Int_t j=0; j<crate->GetNofLocalBoards(); ++j ) 
    {
      Int_t localBoardId = crate->GetLocalBoardId(j);
      AliMpLocalBoard* localBoard 
        = fRegionalTrigger.FindLocalBoard(localBoardId);
      if (!localBoard ) {
        AliFatalClass("Cannot find local board.");
        return kFALSE;
      }   
      
      // Loop over DEs in this localBoard 
        
      for ( Int_t k=0; k<localBoard->GetNofDEs(); ++k ) 
      {
    
        Int_t deId = localBoard->GetDEId(k);
        AliMpDetElement* de = AliMpDEManager::GetDetElement(deId);

        if ( de->GetDdlId() == -1 ) de->SetDdlId(iDDL);

        if ( ! ddl->HasDEId(deId) ) ddl->AddDE(deId);
      }   
    }
  }
  return kTRUE;
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

                    detElement->AddManu(manus[im]);
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
Bool_t AliMpDDLStore::ReadBusPatchSpecial()
{
/// Read file with bus patches with a special order of manus
/// and reset the manus arrays filled via SetManu function

  istream& in 
    = fkDataStreams.
        CreateDataStream(AliMpFiles::BusPatchSpecialFilePath());

  char line[255];

  while ( in.getline(line,255) ) {

    if ( line[0] == '#' ) continue;

    TString tmp(AliMpHelper::Normalize(line));
    TObjArray* stringList = tmp.Tokenize(TString(" "));

    TString sKey = ((TObjString*)stringList->At(0))->GetString();
    
    TString sDDL = ((TObjString*)stringList->At(1))->GetString();
    TArrayI ddlList;
    AliMpHelper::DecodeName(sDDL,';',ddlList);

    TString sBusPatch = ((TObjString*)stringList->At(2))->GetString();
    TArrayI busPatchList;
    AliMpHelper::DecodeName(sBusPatch,',',busPatchList);
    
    // Loop over DDL and Bus Patch
    for (Int_t iDDL = 0; iDDL < ddlList.GetSize(); ++iDDL ) {
      for (Int_t iBusPatch = 0; iBusPatch < busPatchList.GetSize(); ++iBusPatch) {
        // Global bus patch ID
        Int_t busPatchID 
          = AliMpBusPatch::GetGlobalBusID(
              busPatchList.At(iBusPatch), ddlList.At(iDDL));
        
        // Get this bus patch
        AliMpBusPatch* busPatch = GetBusPatch(busPatchID);
        if ( ! busPatch ) {
          AliErrorStream() << "Bus patch " << busPatchID << " does not exist." << endl;
          return kFALSE;
        }
     
        if ( sKey == fgkRevertKeyword ) {
          AliDebugStream(3)
            << "Reverting readout of bus patch " << busPatchID  << endl;
        
          // Now revert the manus in this bus patch
          busPatch->RevertReadout();     
        }
        else if ( sKey == fgkExplicitKeyword ) {
        
          busPatch->ResetReadout();

          TString sManus = ((TObjString*)stringList->At(3))->GetString();
          TArrayI manuList;
          AliMpHelper::DecodeName(sManus,',',manuList);

          AliDebugStream(3)
            << "Reseting readout of bus patch " << busPatchID  
            << "  manus: " << sManus.Data() << endl;

          for (Int_t i = 0; i < manuList.GetSize(); i++) {
            busPatch->AddManu(manuList.At(i));
          }
        }
        else {             
          AliErrorStream() << "Unrecognized key." << endl;
          return kFALSE;
        }
      }
    }
    delete stringList;
  }
  
  delete &in;
  
  return kTRUE;
}    
 

//______________________________________________________________________________
Bool_t AliMpDDLStore::SetPatchModules() 
{
    /// Compute the number of manu per PCB for each buspatch

    AliMpDEIterator it;
    Bool_t result = true;

    for ( it.First(); !it.IsDone(); it.Next() ) {

        AliMpDetElement* detElement = it.CurrentDE();

        for (Int_t i = 0; i < detElement->GetNofBusPatches(); ++i) {
            AliMpBusPatch* busPatch = GetBusPatch(detElement->GetBusPatchId(i));
            Bool_t newResult = false;
            Int_t idDE = busPatch->GetDEId();

            if (AliMpDEManager::GetStation12Type(idDE) == AliMq::kStation2 )
                newResult = busPatch->SetNofManusPerModule(fManuBridge2[GetManuListIndex(idDE)].At(i));
            else
                newResult = busPatch->SetNofManusPerModule();
        }
    }

    return result;
}

//______________________________________________________________________________
Bool_t AliMpDDLStore::ReadBusPatchInfo() 
{
    /// read the buspatch info file and set buspatch info

    istream& in 
      = fkDataStreams.
          CreateDataStream(AliMpFiles::BusPatchInfoFilePath());

    char line[255];

    for (Int_t iDDL = 0; iDDL < fgkNofDDLs; ++iDDL ) {
        AliMpDDL* ddl = GetDDL(iDDL);

        for (Int_t iBusPatch = 0; iBusPatch < ddl->GetNofBusPatches(); ++iBusPatch) {

            do {
                if (!in.getline(line,255)) {
                    AliWarning(Form("Wrong size in bus patch length file; index %d DDL %d",
                                    iBusPatch, iDDL));
                    return false;
                }
            } while(line[0] == '#');

            TString tmp(AliMpHelper::Normalize(line));

            TObjArray* stringList = tmp.Tokenize(TString(" "));

            // Crocus label
            TString crLabel    = ((TObjString*)stringList->At(0))->GetString();
            Int_t pos          = crLabel.First('-');
            tmp                = crLabel(pos-2, crLabel.Length()-pos+2);
            TArrayI list;
            AliMpHelper::DecodeName(tmp.Data(), '-', list);
            
            Int_t localDDLId  = list[0];
            Int_t frtId       = list[1] - 1; // begin at zero ! 
            Int_t localBusId  = list[2];

            // Add FRT number for given ddl if not present
            if ( !ddl->HasFrtId(frtId) )
              ddl->AddFrt(frtId);

            // BP & translator label
            TString label      = ((TObjString*)stringList->At(1))->GetString();
            TString transLabel = ((TObjString*)stringList->At(2))->GetString();

            // BP length
            TString sLength    = ((TObjString*)stringList->At(3))->GetString();
            Float_t length     = sLength.Atof();

            delete stringList;
                       
            if (localBusId != iBusPatch + 1)
               AliWarning(Form("Wrong local buspatch id %d instead of %d", iBusPatch+1, localBusId));
               
            if(localDDLId != ddl->GetId()+1)
                AliWarning(Form("Wrong local DDL id %d instead of %d", ddl->GetId()+1, localDDLId));

            Int_t busPatchId = ddl->GetBusPatchId(iBusPatch);
            AliMpBusPatch* busPatch = GetBusPatch(busPatchId);
            busPatch->SetCableLength(length);
            busPatch->SetCableLabel(label);
            busPatch->SetTranslatorLabel(transLabel);
            busPatch->SetFrtId(frtId);

        }
    }
    
    delete &in;

    return true;
}


//________________________________________________________________
Int_t AliMpDDLStore::GetLocalBoardId(TString name) const {
    /// return the first board with a given side and line

  TIter next(fRegionalTrigger.CreateLocalBoardIterator());
  AliMpLocalBoard* local;
  
  while ( ( local = static_cast<AliMpLocalBoard*>(next()) ) )
  {
        TString tmp(&local->GetName()[4], 2);
        if (name.Contains(tmp))
            if (name[0] == local->GetName()[0])
                return local->GetId();
    }

    return 0;
}

//
// public methods
//


//______________________________________________________________________________
AliMpDDL* AliMpDDLStore::GetDDL(Int_t ddlId, Bool_t warn) const {
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
AliMpDetElement*  AliMpDDLStore::GetDetElement(Int_t detElemId, Bool_t warn) const {
    /// Return detection element with given detElemId

    if ( ! AliMpDEStore::Instance() ) {
        AliFatal("DE Store has not been loaded.");
        return 0;
    }

    return AliMpDEStore::Instance()->GetDetElement(detElemId, warn);
}

//______________________________________________________________________________
AliMpBusPatch* AliMpDDLStore::GetBusPatch(Int_t busPatchId, Bool_t warn) const {
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
AliMpLocalBoard* AliMpDDLStore::GetLocalBoard(Int_t localBoardId, Bool_t warn) const {
    /// Return bus patch with given Id

    return fRegionalTrigger.FindLocalBoard(localBoardId, warn);
}

//______________________________________________________________________________
AliMpTriggerCrate* AliMpDDLStore::GetTriggerCrate(TString name, Bool_t warn) const  {
    /// Return trigger crate with given name

    return fRegionalTrigger.FindTriggerCrate(name, warn);
}

//______________________________________________________________________________
AliMpTriggerCrate* AliMpDDLStore::GetTriggerCrate(Int_t ddlId, Int_t index, Bool_t warn) const  {
    /// Return trigger crate with given ddl and index crate

    if (ddlId == 0 || ddlId == 1)
        ddlId += fgkNofDDLs;

    AliMpDDL* ddl = GetDDL(ddlId, warn);
    if ( ! ddl )
        return 0;

    if ( index >= ddl->GetNofTriggerCrates() ) {
        AliError(Form("crate id %d greater than array[%d]", index, ddl->GetNofTriggerCrates()));
        return 0;
    }

    TString name = AliMpTriggerCrate::GenerateName(index, ddlId, fgkNofDDLs);

    return GetTriggerCrate(name, warn);
}

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetDEfromBus(Int_t busPatchId) const {
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
Int_t  AliMpDDLStore::GetDEfromLocalBoard(Int_t localBoardId, Int_t chamberId) const {
    /// Return detElemId for local board Id and chamber id.

    AliMpLocalBoard* localBoard = GetLocalBoard(localBoardId);

    if ( ! localBoard ) {
        AliErrorStream()
        << "Loacl board with Id = " << localBoardId << " not defined." << endl;
        return 0;
    }

    return localBoard->GetDEIdByChamber(chamberId);
}

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetDDLfromBus(Int_t busPatchId) const {
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
Int_t AliMpDDLStore::GetBusPatchId(Int_t detElemId, Int_t manuId) const {
    /// Return bus patch for a given manuId

    AliMpDetElement* detElement = GetDetElement(detElemId);
    Int_t pos = GetBusPatchIndex(detElemId, manuId);

    if ( pos >= detElement->GetNofBusPatches() ) 
    {
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
Long_t AliMpDDLStore::GetLinkPortId(Int_t busPatchId) const {

    /// Get link port and DSP from busPatch id.
    /// Return -1 if the value is not valid 
    /// (the validity has to be tested in the client code)

    AliMpBusPatch* busPatch = GetBusPatch(busPatchId);
    Int_t ddlId = busPatch->GetDdlId();
        
    Int_t localBusPatchId = AliMpBusPatch::GetLocalBusID(busPatchId, ddlId) - 1; // begin at zero

    Int_t pos = (localBusPatchId % AliMpFrtCrocusConstants::GetNofBusPatches()); 
    
    return AliMpFrtCrocusConstants::GetLinkPortId(pos);

}

//______________________________________________________________________________
void AliMpDDLStore::PrintAllManu() const {
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

            if ( AliMpManuStore::Instance(kFALSE) ) {
              cout << "    Manu serial: ";
              for ( Int_t k=0; k<busPatch->GetNofManus(); ++k ) {
                cout << std::setw(6) 
                     << AliMpManuStore::Instance()
                        ->GetManuSerial(de->GetId(), busPatch->GetManuId(k)) << " ";
              }
              cout << endl;
            }  
        }
    }
}

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetNextDEfromLocalBoard(Int_t localBoardId, Int_t chamberId ) const {
    /// return the next detection element in line

    AliMpLocalBoard* localBoard  =  GetLocalBoard(localBoardId);

    TString name(localBoard->GetName());

    Int_t line = AliMp::PairFirst(localBoard->GetPosition());
    ++line;

    name.Replace(4,1,Form("%d", line));

    Int_t nextLocalId;
    if ((nextLocalId = GetLocalBoardId(name)))
        return GetDEfromLocalBoard(nextLocalId, chamberId);
    else
        return 0;

    return 0;
}

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetPreviousDEfromLocalBoard(Int_t localBoardId, Int_t chamberId) const {
    /// return the previous detection element in line

    AliMpLocalBoard* localBoard  =  GetLocalBoard(localBoardId);

    TString name(localBoard->GetName());

    Int_t line = AliMp::PairFirst(localBoard->GetPosition());
    --line;

    name.Replace(4,1,Form("%d", line));

    Int_t prevLocalId;
    if ((prevLocalId = GetLocalBoardId(name)))
        return GetDEfromLocalBoard(prevLocalId, chamberId);
    else
        return 0;

}

//______________________________________________________________________________
void AliMpDDLStore::SetRegionalTrigger(const AliMpRegionalTrigger& regionalTrigger)
{
/// Replace the existing regional trigger with the given one

  fRegionalTrigger = regionalTrigger;
  
  // Remove the existing trigger DDLsf
  fDDLs.RemoveAt(fgkNofDDLs+1);
  fDDLs.RemoveAt(fgkNofDDLs);
  
  // Set new trigger DDLs from new regional trigger
  SetTriggerDDLs();
}  


//______________________________________________________________________________
TIterator* 
AliMpDDLStore::CreateBusPatchIterator() const
{
/// Create the iterator over bus patches

  return fBusPatches.CreateIterator();
}
