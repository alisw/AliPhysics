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

/// \class AliMUONData
///
/// Class containing MUON data: hits, digits, rawclusters, globaltrigger, localtrigger, etc ..
/// The classe makes the lik between the MUON data lists and the event trees from loaders
///
/// \author Gines Martinez, Subatech,  September 2003
///

#include "AliMUONData.h"
#include "AliMUONDataIterator.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONRawCluster.h"

#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliLog.h"

#include <TString.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <Riostream.h>
#include <TFile.h>

/// \cond CLASSIMP
ClassImp(AliMUONData)
/// \endcond
 
//_____________________________________________________________________________
  AliMUONData::AliMUONData():
    TNamed(),
    fRunLoader(0x0),
    fLoader(0x0),
    fSDigits(0x0),
    fDigits(0x0),
    fGlobalTrigger(0x0),
    fLocalTrigger(0x0),
    fRegionalTrigger(0x0),
    fNSdigits(0x0),
    fNdigits(0x0),
    fNglobaltrigger(0),
    fNlocaltrigger(0),
    fNregionaltrigger(0),
    fSplitLevel(0),
    fCurrentEvent(-1)
{
/// Default constructor
}
//_____________________________________________________________________________
AliMUONData::AliMUONData(AliLoader * loader, const char* name, const char* title):
  TNamed(name,title),
    fRunLoader(0x0),
    fLoader(loader),
    fSDigits(0x0),
    fDigits(0x0),
    fGlobalTrigger(0x0),
    fLocalTrigger(0x0),
    fRegionalTrigger(0x0),
    fNSdigits(0x0),
    fNdigits(0x0),
    fNglobaltrigger(0),
    fNlocaltrigger(0),
    fNregionaltrigger(0),
    fSplitLevel(0),
    fCurrentEvent(-1)
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMUONData::AliMUONData(const char* galiceFile, const char* folderName):
  TNamed("MUON", "MUON"),
    fRunLoader(0x0),
    fLoader(0x0),
    fSDigits(0x0),
    fDigits(0x0),
    fGlobalTrigger(0x0),
    fLocalTrigger(0x0),
    fRegionalTrigger(0x0),
    fNSdigits(0x0),
    fNdigits(0x0),
    fNglobaltrigger(0),
    fNlocaltrigger(0),
    fNregionaltrigger(0),
    fSplitLevel(0),
    fCurrentEvent(-1)
{
/// Constructor for loading data from gAlice file

  fRunLoader = AliRunLoader::Open(galiceFile, folderName, "READ");
  if (!fRunLoader) {
    AliError(Form("Error opening %s file \n", galiceFile));
    return;
  }  

  fLoader = fRunLoader->GetLoader("MUONLoader");
  if ( ! fLoader ) {
    AliError(Form("Could get MUONLoader"));
    return;
  }  
}

//_____________________________________________________________________________
AliMUONData::~AliMUONData()
{
/// Destructor for AliMUONData
  
  if (fSDigits) {
    fSDigits->Delete();
    delete fSDigits;
  }

  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
  }
  if (fGlobalTrigger){
    fGlobalTrigger->Delete();
    delete fGlobalTrigger;
  }  
  if (fRegionalTrigger){
    fRegionalTrigger->Delete();
    delete fRegionalTrigger;
  }
  if (fLocalTrigger){
    fLocalTrigger->Delete();
    delete fLocalTrigger;
  }

  if (fRunLoader) {
    fRunLoader->UnloadAll();
    delete fRunLoader;
  }  
}
//_____________________________________________________________________________
void AliMUONData::AddSDigit(Int_t id, const AliMUONDigit& Sdigit)
{
/// Add a MUON Sdigit to the list of SDigits of the detection plane id

  TClonesArray &lSdigits = * SDigits(id) ; 
  new(lSdigits[fNSdigits[id]++]) AliMUONDigit(Sdigit);
}
//_____________________________________________________________________________
void AliMUONData::AddDigit(Int_t id, const AliMUONDigit& digit)
{
/// Add a MUON digit to the list of Digits of the detection plane id

  TClonesArray &ldigits = * Digits(id) ; 
  new(ldigits[fNdigits[id]++]) AliMUONDigit(digit);
}

//_____________________________________________________________________________
void AliMUONData::AddGlobalTrigger(const AliMUONGlobalTrigger& trigger )
{
/// Add a MUON Global Trigger to the list (only one GlobalTrigger per event !);

  TClonesArray &globalTrigger = *fGlobalTrigger;
  new(globalTrigger[fNglobaltrigger++]) AliMUONGlobalTrigger(trigger);
}

//____________________________________________________________________________
void AliMUONData::AddRegionalTrigger(const  AliMUONRegionalTrigger& trigger)
{
/// add a MUON regional Trigger to the list
  TClonesArray &regionalTrigger = *fRegionalTrigger;
  new(regionalTrigger[fNregionaltrigger++]) AliMUONRegionalTrigger(trigger);
}
//____________________________________________________________________________
void AliMUONData::AddLocalTrigger(const  AliMUONLocalTrigger& trigger)
{
/// add a MUON Local Trigger to the list

  TClonesArray &localTrigger = *fLocalTrigger;
  new(localTrigger[fNlocaltrigger++]) AliMUONLocalTrigger(trigger);
}

//____________________________________________________________________________
TClonesArray*  AliMUONData::SDigits(Int_t DetectionPlane) const
{
/// Getting List of SDigits

  if (fSDigits)
    return ( (TClonesArray*) fSDigits->At(DetectionPlane) );
  else
    return NULL;
}
//____________________________________________________________________________
TClonesArray*  AliMUONData::Digits(Int_t DetectionPlane) const
{
/// Getting List of Digits

  if (fDigits)
    return ( (TClonesArray*) fDigits->At(DetectionPlane) );
  else
    return NULL;
}
//____________________________________________________________________________
Bool_t   AliMUONData::IsDigitsBranchesInTree()
{
/// Checking if there are Digits Branches In TreeD

  if (TreeD()==0x0) {
    AliError("No treeD in memory");
    return kFALSE;
  }
  else {
     char branchname[30];
     sprintf(branchname,"%sDigits1",GetName());
     TBranch * branch = 0x0;
     branch = TreeD()->GetBranch(branchname);
     if (branch)  return kTRUE;
     else return kFALSE;    
  }
}
//____________________________________________________________________________
Bool_t   AliMUONData::IsTriggerBranchesInTreeD()
{
/// Checking if there are Trigger Branches In TreeD
 if (TreeD()==0x0) {
    AliError("No treeD in memory");
    return kFALSE;
  }
  else {
     char branchname[30];
     sprintf(branchname,"%sLocalTrigger",GetName());
     TBranch * branch = 0x0;
     branch = TreeD()->GetBranch(branchname);
     if (branch)  return kTRUE;
     else return kFALSE;    
  }
}

//____________________________________________________________________________
void AliMUONData::Fill(Option_t* option)
{
/// Method to fill the trees

  const char *cS   = strstr(option,"S");   // SDigits branches in TreeS
  const char *cD   = strstr(option,"D");   // Digits branches in TreeD
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger branches in TreeD
  
  char branchname[30];
  TBranch * branch = 0x0;

  // Filling TreeS
  if ( TreeS() && cS) 
  {
    TreeS()->Fill();
  }

  // Filling TreeD

  if ( TreeD() && cD && cGLT )
  {
    // Writing digits and (global+local) trigger at once.
    TreeD()->Fill();
  }
  else
  {
    if ( TreeD() && cD ) 
    {
      if ( IsTriggerBranchesInTreeD() ) 
      {
        for (int i=0; i<AliMUONConstants::NCh(); i++) 
        {
          sprintf(branchname,"%sDigits%d",GetName(),i+1);
          branch = TreeD()->GetBranch(branchname);
          branch->Fill();
        }
      } 
      else
      {
        TreeD()->Fill();
      }
    }
    
    if ( TreeD() && cGLT ) 
    {
      if ( IsDigitsBranchesInTree() ) 
      {
        sprintf(branchname,"%sLocalTrigger",GetName());
        branch = TreeD()->GetBranch(branchname); 
        branch->Fill();
	sprintf(branchname,"%sRegionalTrigger",GetName());
        branch = TreeD()->GetBranch(branchname);
        branch->Fill();
        sprintf(branchname,"%sGlobalTrigger",GetName());
        branch = TreeD()->GetBranch(branchname);
        branch->Fill();

      } 
      else
      {
        TreeD()->Fill();
      }
    }
  } // end of TreeD() handling.
}

//_____________________________________________________________________________
void AliMUONData::MakeBranch(Option_t* option)
{
/// Create Tree branches for the MUON.

  const Int_t kBufferSize = 4000;
  char branchname[30];
  
  //Setting Data Container
  SetDataContainer(option);  

  const char *cS   = strstr(option,"S");   // Digits branches in TreeS
  const char *cD   = strstr(option,"D");   // Digits branches in TreeD
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger branches in TreeD
  
  TBranch * branch = 0x0;
  
  //Creating Branches for SDigits
  if (TreeS() && cS ) {
    // one branch for Sdigits per chamber
    for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
      sprintf(branchname,"%sSDigits%d",GetName(),iDetectionPlane+1);
      branch = 0x0;
      branch = TreeS()->GetBranch(branchname);
      if (branch) {  
        AliInfo(Form("Branch %s is already in tree.",branchname));
        return;
      }
      TClonesArray * sdigits = SDigits(iDetectionPlane); 
      branch = TreeS()->Branch(branchname, &sdigits, kBufferSize,1);
      //Info("MakeBranch","Making Branch %s for sdigits in detection plane %d\n",branchname,iDetectionPlane+1);
    }
  }

  //Creating Branches for Digits
  TTree* treeD = 0x0;
  if ( cD || cGLT )
  {
    treeD = TreeD();
  }

  if ( treeD && cD ) 
  {
    // one branch for digits per chamber
    for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) 
    {
      sprintf(branchname,"%sDigits%d",GetName(),iDetectionPlane+1);
      branch = treeD->GetBranch(branchname);
      if (branch) 
      {  
        AliInfo(Form("Branch %s is already in tree.",branchname));
        return;
      }
      TClonesArray * digits = Digits(iDetectionPlane); 
      branch = treeD->Branch(branchname, &digits, kBufferSize,1);
    }
  }
  
  if ( treeD && cGLT ) 
  {
    //
    // one branch for global trigger
    //
    sprintf(branchname,"%sGlobalTrigger",GetName());
    branch = treeD->GetBranch(branchname);
    if (branch) 
    {  
      AliInfo(Form("Branch GlobalTrigger is already in treeD."));
      return ;
    }
    branch = treeD->Branch(branchname, &fGlobalTrigger, kBufferSize);

  //
    // one branch for regional trigger
    //  
    sprintf(branchname,"%sRegionalTrigger",GetName());
    branch = 0x0;
    branch = treeD->GetBranch(branchname);
    if (branch) 
    {  
      AliInfo(Form("Branch RegionalTrigger is already in treeD."));
      return;
    }
    branch = treeD->Branch(branchname, &fRegionalTrigger, kBufferSize);
  

    //
    // one branch for local trigger
    //  
    sprintf(branchname,"%sLocalTrigger",GetName());
    branch = 0x0;
    branch = treeD->GetBranch(branchname);
    if (branch) 
    {  
      AliInfo(Form("Branch LocalTrigger is already in treeD."));
      return;
    }
    branch = treeD->Branch(branchname, &fLocalTrigger, kBufferSize);
  }
}

//____________________________________________________________________________
TClonesArray*  
AliMUONData::LocalTrigger() const
{
/// Getting local trigger

  return fLocalTrigger;
}

//____________________________________________________________________________
TClonesArray*  
AliMUONData::RegionalTrigger() const
{
/// Getting regional trigger

  return fRegionalTrigger;
}

//____________________________________________________________________________
void
AliMUONData::GetDigits() const 
{
/// Load the digits from TreeD for the current event.

  Int_t event = fLoader->GetRunLoader()->GetEventNumber();
  if ( fCurrentEvent != event )
  {
    if (fLoader->TreeD()) {
      fLoader->TreeD()->GetEvent(0);
      fCurrentEvent = event;
    }
  }
}

//____________________________________________________________________________
TClonesArray*  
AliMUONData::GlobalTrigger() const
{
/// Return the global trigger 

  return fGlobalTrigger;
}

//____________________________________________________________________________
void AliMUONData::ResetSDigits()
{
/// Reset number of Sdigits and the Sdigits array for this detector

    if (fSDigits == 0x0) return;
    for ( int i=0;i<AliMUONConstants::NCh();i++ ) {
      if ((*fSDigits)[i])    ((TClonesArray*)fSDigits->At(i))->Clear();
      if (fNSdigits)  fNSdigits[i]=0;
    }
}
//____________________________________________________________________________
void AliMUONData::ResetDigits()
{
/// Reset number of digits and the digits array for this detector

    if (fDigits == 0x0) return;
    for ( int i=0;i<AliMUONConstants::NCh();i++ ) {
      if ((*fDigits)[i])    ((TClonesArray*)fDigits->At(i))->Clear("C");
      if (fNdigits)  fNdigits[i]=0;
    }
}
//_______________________________________________________________________________
void AliMUONData::ResetTrigger()
{
/// Reset Local and Global Trigger 

  fNglobaltrigger = 0;
  if (fGlobalTrigger) fGlobalTrigger->Clear();
  fNregionaltrigger = 0;
  if (fRegionalTrigger) fRegionalTrigger->Clear();
  fNlocaltrigger = 0;
  if (fLocalTrigger) fLocalTrigger->Clear();

}
//____________________________________________________________________________
void AliMUONData::SetDataContainer(Option_t* option)
{
/// Setting data containers of muon data

  const char *cS   = strstr(option,"S");   // SDigits
  const char *cD   = strstr(option,"D");   // Digits
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger

  AliDebug(1,Form("option=%s",option));
  
  //
  // Container for Sdigits
  if (cS) {
    if (fSDigits == 0x0) { 
      AliDebug(1,"Creating fSDigits TObjArray");
      fSDigits = new TObjArray(AliMUONConstants::NCh());
      fNSdigits= new Int_t[AliMUONConstants::NCh()];
      for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
	TClonesArray* a = new TClonesArray("AliMUONDigit",10000);
	a->SetOwner();
	fSDigits->AddAt(a,i);
	AliDebug(1,Form("fSDigits[%d]=%p",i,a));
        fNSdigits[i]=0;
      }
    }
    else {
      AliDebug(1,Form("fSDigits already there = %p",fSDigits));
    }
    ResetSDigits();
  }  

  //
  // ObjArray of ClonesArrays for Digits
  if ( cD ) {      
    if (fDigits == 0x0 ) {
      fDigits = new TObjArray(AliMUONConstants::NCh());
      fNdigits= new Int_t[AliMUONConstants::NCh()];
      for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
	TClonesArray * tca = new TClonesArray("AliMUONDigit",10000);
	tca->SetOwner();
        fDigits->AddAt(tca,i); 
        fNdigits[i]=0;
      }
    } 
    else {
      AliDebug(1,Form("fDigits already there = %p",fDigits));
    }
    ResetDigits();
  }

  //
  // ClonesArrays for Trigger
  if ( cGLT ) { 
    if (fLocalTrigger == 0x0) {
      fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);
    }
    if (fRegionalTrigger == 0x0) {
      fRegionalTrigger  = new TClonesArray("AliMUONRegionalTrigger",16);
    }
    if (fGlobalTrigger== 0x0) {
      fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1); 
    }
    ResetTrigger();
  }
}

//____________________________________________________________________________
void AliMUONData::SetTreeAddress(Option_t* option)
{
  // Setting Data containers
  SetDataContainer(option);

/// Setting Addresses to the events trees

  const char *cS   = strstr(option,"S");   // SDigits branches in TreeS
  const char *cD   = strstr(option,"D");   // Digits branches in TreeD
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger branches in TreeD
  
  // Set branch address for the Hits, Digits, RawClusters, GlobalTrigger and LocalTrigger Tree.
  char branchname[30];
  TBranch * branch = 0x0;
  
  AliDebug(1,Form("option=%s",option));
  
  //
  // Branch address for Sdigit tree
  if (TreeS() && fSDigits && cS) {
    AliDebug(1,"Setting branch addresses");
    for (int i=0; i<AliMUONConstants::NCh(); i++) {
      sprintf(branchname,"%sSDigits%d",GetName(),i+1);
      if (fSDigits) {
        AliDebug(1,Form("TreeS=%p for ich=%d branchname=%s",
                        TreeS(),i,branchname));
        branch = TreeS()->GetBranch(branchname);
        TClonesArray * sdigits = SDigits(i);
        if (branch) branch->SetAddress( &sdigits );
        else AliWarning(Form("(%s) Failed for SDigits Detection plane %d. Can not find branch in tree.",GetName(),i));
      }
    }
  }

  //
  // Branch address for digit tree
  if (TreeD() && fDigits && cD) {
    for (int i=0; i<AliMUONConstants::NCh(); i++) {
      sprintf(branchname,"%sDigits%d",GetName(),i+1);
      if (fDigits) {
        branch = TreeD()->GetBranch(branchname);
        TClonesArray * digits = Digits(i);
        if (branch) {
          branch->SetAddress( &digits );
        }
        else AliWarning(Form("(%s) Failed for Digits Detection plane %d. Can not find branch in tree.",GetName(),i));
      }
    }
  }
  if ( TreeD()  && fLocalTrigger && cGLT) {
    sprintf(branchname,"%sLocalTrigger",GetName());
    branch = TreeD()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fLocalTrigger);
    else AliWarning(Form("(%s) Failed for LocalTrigger. Can not find branch in treeD.",GetName()));
  }
 if ( TreeD()  && fRegionalTrigger && cGLT) {
    sprintf(branchname,"%sRegionalTrigger",GetName());
    branch = TreeD()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fRegionalTrigger);
    else AliWarning(Form("(%s) Failed for RegionalTrigger. Can not find branch in treeD.",GetName()));
  }
  if ( TreeD() && fGlobalTrigger && cGLT) {
    sprintf(branchname,"%sGlobalTrigger",GetName());
    branch = TreeD()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fGlobalTrigger);
    else AliWarning(Form("(%s) Failed for GlobalTrigger. Can not find branch in treeD.",GetName()));
  }
}

//_____________________________________________________________________________
void
AliMUONData::Print(Option_t* opt) const
{
/// Dump object on screen

  TString options(opt);
  options.ToUpper();
  
  if ( options.Contains("D") )
  {
    for ( Int_t ich = 0; ich < AliMUONConstants::NCh(); ++ich)
    {
      TClonesArray* digits = Digits(ich);
      Int_t ndigits = digits->GetEntriesFast();
      for ( Int_t id = 0; id < ndigits; ++id )
      {
        AliMUONDigit* digit = 
          static_cast<AliMUONDigit*>(digits->UncheckedAt(id));
        digit->Print();
      }
    }
  }

  if ( options.Contains("S") )
  {
    for ( Int_t ich = 0; ich < AliMUONConstants::NCh(); ++ich)
    {
      TClonesArray* digits = SDigits(ich);
      Int_t ndigits = digits->GetEntriesFast();
      for ( Int_t id = 0; id < ndigits; ++id )
      {
        AliMUONDigit* digit = 
        static_cast<AliMUONDigit*>(digits->UncheckedAt(id));
        digit->Print();
      }
    }
  }  
}

//_____________________________________________________________________________
void 
AliMUONData::DumpSDigits(Int_t event2Check, Option_t* opt)
{
/// Dump SDigits

  fLoader->LoadSDigits("READ");
  
  // Event loop
  Int_t nevents = fRunLoader->GetNumberOfEvents();
  for (Int_t ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);

    // Getting event ievent
    fRunLoader->GetEvent(ievent);
    SetTreeAddress("S");
    GetSDigits();

    // Loop on chambers
    Int_t nchambers = AliMUONConstants::NCh(); ;
    for (Int_t ichamber=0; ichamber<nchambers; ichamber++) {
      TClonesArray* digits = SDigits(ichamber);

      // Loop on Sdigits
      Int_t ndigits = (Int_t)digits->GetEntriesFast();
      for (Int_t idigit=0; idigit<ndigits; idigit++) {
        AliMUONDigit* mDigit = static_cast<AliMUONDigit*>(digits->At(idigit));
        mDigit->Print(opt);
      }
    }
    ResetSDigits();
    if (event2Check!=0) ievent=nevents;
  }
  fLoader->UnloadSDigits();
}
//_____________________________________________________________________________
void 
AliMUONData::DumpDigits(Int_t event2Check, Option_t* opt)
{
/// Dump digits

  fLoader->LoadDigits("READ");
  
  // Event loop
  Int_t firstEvent = 0;
  Int_t lastEvent = fRunLoader->GetNumberOfEvents()-1;
  if ( event2Check != 0 ) {
    firstEvent = event2Check;
    lastEvent = event2Check;
  }  
  
  for ( Int_t ievent = firstEvent; ievent <= lastEvent; ++ievent ) {
    printf(">>> Event %d \n",ievent);
    fRunLoader->GetEvent(ievent);

    AliMUONDataIterator it(this, "digit", AliMUONDataIterator::kTrackingChambers);
    AliMUONDigit* digit;
 
     while ( ( digit = (AliMUONDigit*)it.Next() ) )
     {
       digit->Print(opt);
     }
  } 
  fLoader->UnloadDigits();
}
