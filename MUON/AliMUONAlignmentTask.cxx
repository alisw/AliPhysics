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

//-----------------------------------------------------------------------------
/// \class AliMUONAlignmentTask
/// AliAnalysisTask to align the MUON spectrometer.
/// The Task reads as input ESDs and feeds the MUONTracks to AliMUONAlignment.
/// The alignment itself is performed by AliMillepede.
/// A OCDB entry is written with the alignment parameters.
/// A list of graph are written to monitor the alignment parameters.
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN
//-----------------------------------------------------------------------------

#include <fstream>

#include <TString.h>
#include <TError.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <Riostream.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliMagF.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliGeomManager.h"

#include "AliMpCDB.h"
#include "AliMUONAlignment.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONESDInterface.h"

#include "AliMUONAlignmentTask.h"

///\cond CLASSIMP  
ClassImp(AliMUONAlignmentTask)
///\endcond

// //________________________________________________________________________
// AliMUONAlignmentTask::AliMUONAlignmentTask(const char *name) 
//   : AliAnalysisTask(name, ""),
//     fESD(0x0),
//     fAlign(0x0),
//     fGeoFilename(""),
//     fTransform(0x0),
//     fTrackTot(0),
//     fTrackOk(0),
//     fMSDEx(0x0), 
//     fMSDEy(0x0), 
//     fMSDEz(0x0),
//     fMSDEp(0x0),
//     fList(0x0)
// {
//   /// Default Constructor
//   // Define input and output slots here
//   // Input slot #0 works with a TChain
//   DefineInput(0, TChain::Class());
//   // Output slot #0 writes NTuple/histos into a TList
//   DefineOutput(0, TList::Class());  

//   // initialize parameters ...
//   for(Int_t k=0;k<4*156;k++) {
//     fParameters[k]=0.;
//     fErrors[k]=0.;
//     fPulls[k]=0.;
//   }

//   fAlign = new AliMUONAlignment();
//   fTransform = new AliMUONGeometryTransformer();  
// }

//________________________________________________________________________
AliMUONAlignmentTask::AliMUONAlignmentTask(const char *name, const char *geofilename, const char *defaultocdb, const char *misalignocdb) 
  : AliAnalysisTask(name, ""),
    fESD(0x0),
    fAlign(0x0),
    fGeoFilename(geofilename),
    fMisAlignOCDB(misalignocdb),
    fDefaultOCDB(defaultocdb),
    fTransform(0x0),
    fTrackTot(0),
    fTrackOk(0),
    fLastRunNumber(-1),
    fMSDEx(0x0), 
    fMSDEy(0x0), 
    fMSDEz(0x0),
    fMSDEp(0x0),
    fList(0x0)
{
  /// Default Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes NTuple/histos into a TList
  DefineOutput(0, TList::Class());  

  // initialize parameters ...
  for(Int_t k=0;k<4*156;k++) {
    fParameters[k]=0.;
    fErrors[k]=0.;
    fPulls[k]=0.;
  }

  fAlign = new AliMUONAlignment();
  fTransform = new AliMUONGeometryTransformer();  
}


//________________________________________________________________________
AliMUONAlignmentTask::AliMUONAlignmentTask(const AliMUONAlignmentTask& obj) 
  : AliAnalysisTask(obj),
    fESD(0x0),
    fAlign(0x0),
    fGeoFilename(""),
    fMisAlignOCDB(""),
    fDefaultOCDB(""),
    fTransform(0x0),
    fTrackTot(0),
    fTrackOk(0),
    fLastRunNumber(-1),
    fMSDEx(0x0), 
    fMSDEy(0x0), 
    fMSDEz(0x0),
    fMSDEp(0x0),
    fList(0x0)
{
  /// Copy constructor
  fESD = obj.fESD;
  fAlign = obj.fAlign;
  fGeoFilename = obj.fGeoFilename;
  fTransform = obj.fTransform;
  fTrackTot = obj.fTrackTot;  
  fTrackOk = obj.fTrackOk;  
  fLastRunNumber = obj.fLastRunNumber;
  fMSDEx = obj.fMSDEx; 
  fMSDEy = obj.fMSDEy; 
  fMSDEz = obj.fMSDEz;
  fMSDEp = obj.fMSDEp;
  fList = obj.fList;  
}

//________________________________________________________________________
AliMUONAlignmentTask& AliMUONAlignmentTask::operator=(const AliMUONAlignmentTask& other) 
{
  /// Assignment
  AliAnalysisTask::operator=(other);
  fESD = other.fESD;
  fAlign = other.fAlign;
  fGeoFilename = other.fGeoFilename;
  fMisAlignOCDB = other.fMisAlignOCDB;
  fDefaultOCDB = other.fDefaultOCDB;
  fTransform = other.fTransform;
  fTrackTot = other.fTrackTot;  
  fTrackOk = other.fTrackOk;  
  fLastRunNumber = other.fLastRunNumber;
  fMSDEx = other.fMSDEx; 
  fMSDEy = other.fMSDEy; 
  fMSDEz = other.fMSDEz;
  fMSDEp = other.fMSDEp;
  fList = other.fList;  
  
  return *this;
}

//________________________________________________________________________
AliMUONAlignmentTask::~AliMUONAlignmentTask() 
{ 
  /// Destructor
  if (fAlign) delete fAlign;
  if (fTransform) delete fTransform;
}

//________________________________________________________________________
void AliMUONAlignmentTask::LocalInit() 
{
  /// Local initialization, called once per task on the client machine 
  /// where the analysis train is assembled
  fLastRunNumber = 0; 
  //  Prepare(fGeoFilename.Data(),fDefaultOCDB.Data(),fMisAlignOCDB.Data());
  Prepare(fGeoFilename.Data(),"local://$ALICE_ROOT/OCDB",fMisAlignOCDB.Data());
  fLastRunNumber = -1;

  // Set initial values here, good guess may help convergence
  // St 1 
  //  Int_t iPar = 0;
  //  fParameters[iPar++] =  0.010300 ;  fParameters[iPar++] =  0.010600 ;  fParameters[iPar++] =  0.000396 ;  


  fAlign->InitGlobalParameters(fParameters);


  fTransform->LoadGeometryData();
  fAlign->SetGeometryTransformer(fTransform);

  // Do alignment with magnetic field off
  fAlign->SetBFieldOn(kFALSE);
  
  // Set tracking station to use
  //  Bool_t bStOnOff[5] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  Bool_t bChOnOff[10] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kFALSE,kTRUE,kTRUE,kTRUE,kTRUE};

  // Set degrees of freedom to align (see AliMUONAlignment)
  fAlign->AllowVariations(bChOnOff);

  // Fix parameters or add constraints here
  //   for (Int_t iSt=0; iSt<5; iSt++)
  //     if (!bStOnOff[iSt]) fAlign->FixStation(iSt+1);
  for (Int_t iCh=0; iCh<10; iCh++)
    if (!bChOnOff[iCh]) fAlign->FixChamber(iCh+1);

  // Left and right sides of the detector are independent, one can choose to align 
  // only one side
  Bool_t bSpecLROnOff[2] = {kTRUE,kTRUE};
  fAlign->FixHalfSpectrometer(bChOnOff,bSpecLROnOff);

  fAlign->SetChOnOff(bChOnOff);
  fAlign->SetSpecLROnOff(bChOnOff);

    // Here we can fix some detection elements
  fAlign->FixDetElem(908);
  fAlign->FixDetElem(1020);

  // Set predifined global constrains: X, Y, P, XvsZ, YvsZ, PvsZ, XvsY, YvsY, PvsY
//   Bool_t bVarXYT[9] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
//   Bool_t bDetTLBR[4] = {kFALSE,kTRUE,kFALSE,kTRUE};
  //  fAlign->AddConstraints(bChOnOff,bVarXYT,bDetTLBR,bSpecLROnOff);

}

//________________________________________________________________________
void AliMUONAlignmentTask::ConnectInputData(Option_t *) 
{
  /// Connect ESD here. Called on each input data change.

  // Connect ESD here
  TTree* esdTree = dynamic_cast<TTree*> (GetInputData(0));
  if (!esdTree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());   
    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } 
    else {
      fESD = esdH->GetEvent();
    }
  }
}

//________________________________________________________________________
void AliMUONAlignmentTask::CreateOutputObjects()
{
  /// Executed once on each worker (machine actually running the analysis code)
  //
  // This method has to be called INSIDE the user redefined CreateOutputObjects
  // method, before creating each object corresponding to the output containers
  // that are to be written to a file. This need to be done in general for the big output
  // objects that may not fit memory during processing. 
  //  OpenFile(0); 

  // Creating graphs 
  fMSDEx = new TGraphErrors(156);
  fMSDEy = new TGraphErrors(156);
  fMSDEz = new TGraphErrors(156);
  fMSDEp = new TGraphErrors(156);

  // Add Ntuples to the list
  fList = new TList();
  fList->Add(fMSDEx);
  fList->Add(fMSDEy);
  fList->Add(fMSDEz);
  fList->Add(fMSDEp);
}

//________________________________________________________________________
void AliMUONAlignmentTask::Exec(Option_t *) 
{
  /// Main loop, called for each event
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  Double_t trackParams[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
  if (fESD->GetRunNumber()!=fLastRunNumber){
    fLastRunNumber = fESD->GetRunNumber();
    Prepare(fGeoFilename.Data(),fDefaultOCDB.Data(),fMisAlignOCDB.Data());
  }
 
  Int_t nTracks = Int_t(fESD->GetNumberOfMuonTracks());
  //  if (!event%100) cout << " there are " << nTracks << " tracks in event " << event << endl;
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliESDMuonTrack* esdTrack = fESD->GetMuonTrack(iTrack);
    if (!esdTrack->ClustersStored()) continue;
    Double_t invBenMom = esdTrack->GetInverseBendingMomentum();
//     fInvBenMom->Fill(invBenMom);
//     fBenMom->Fill(1./invBenMom);
    if (TMath::Abs(invBenMom)<=1.04) {
      AliMUONTrack track;
      AliMUONESDInterface::ESDToMUON(*esdTrack, track);
      fAlign->ProcessTrack(&track);
      fAlign->LocalFit(fTrackOk++,trackParams,0);
    }
    fTrackTot++;
    cout << "Processed " << fTrackTot << " Tracks." << endl;
  }
  
  // Post final data. Write histo list to a file with option "RECREATE"
  PostData(0,fList);
}      

//________________________________________________________________________
void AliMUONAlignmentTask::FinishTaskOutput()
{
  /// Called once per task on the client machine at the end of the analysis.

  cout << "Processed " << fTrackTot << " Tracks." << endl;
  // Perform global fit
  fAlign->GlobalFit(fParameters,fErrors,fPulls);

  cout << "Done with GlobalFit " << endl;

//   // Update pointers reading them from the output slot
//   fList = (TList*)GetOutputData(0);
//   fMSDEx = (TGraphErrors*)fList->At(0);
//   fMSDEy = (TGraphErrors*)fList->At(1);
//   fMSDEz = (TGraphErrors*)fList->At(2);
//   fMSDEp = (TGraphErrors*)fList->At(3);

  // Store results
  Double_t DEid[156] = {0};
  Double_t MSDEx[156] = {0};
  Double_t MSDEy[156] = {0};
  Double_t MSDEz[156] = {0};
  Double_t MSDEp[156] = {0};
  Double_t DEidErr[156] = {0};
  Double_t MSDExErr[156] = {0};
  Double_t MSDEyErr[156] = {0};
  Double_t MSDEzErr[156] = {0};
  Double_t MSDEpErr[156] = {0};
  Int_t lNDetElem = 4*2+4*2+18*2+26*2+26*2;
  Int_t lNDetElemCh[10] = {4,4,4,4,18,18,26,26,26,26};
  // Int_t lSNDetElemCh[10] = {4,8,12,16,34,52,78,104,130,156};
  Int_t idOffset = 0; // 400
  Int_t lSDetElemCh = 0;
  for(Int_t iDE=0; iDE<lNDetElem; iDE++){
    DEidErr[iDE] = 0.;
    DEid[iDE] = idOffset+100;
    DEid[iDE] += iDE; 
    lSDetElemCh = 0;
    for(Int_t iCh=0; iCh<9; iCh++){
      lSDetElemCh += lNDetElemCh[iCh];
      if (iDE>=lSDetElemCh) {
	DEid[iDE] += 100;
	DEid[iDE] -= lNDetElemCh[iCh];
      }
    }
    MSDEx[iDE]=fParameters[3*iDE+0];
    MSDEy[iDE]=fParameters[3*iDE+1];
    MSDEz[iDE]=fParameters[3*iDE+3];
    MSDEp[iDE]=fParameters[3*iDE+2];
    MSDExErr[iDE]=(Double_t)fAlign->GetParError(3*iDE+0);
    MSDEyErr[iDE]=(Double_t)fAlign->GetParError(3*iDE+1);
    MSDEzErr[iDE]=(Double_t)fAlign->GetParError(3*iDE+3);
    MSDEpErr[iDE]=(Double_t)fAlign->GetParError(3*iDE+2);
    fMSDEx->SetPoint(iDE,DEid[iDE],fParameters[3*iDE+0]);
    fMSDEx->SetPoint(iDE,DEidErr[iDE],(Double_t)fAlign->GetParError(3*iDE+0));
    fMSDEy->SetPoint(iDE,DEid[iDE],fParameters[3*iDE+1]);
    fMSDEy->SetPoint(iDE,DEidErr[iDE],(Double_t)fAlign->GetParError(3*iDE+1));
    fMSDEp->SetPoint(iDE,DEid[iDE],fParameters[3*iDE+2]);
    fMSDEp->SetPoint(iDE,DEidErr[iDE],(Double_t)fAlign->GetParError(3*iDE+2));
    fMSDEz->SetPoint(iDE,DEid[iDE],fParameters[3*iDE+3]);
    fMSDEz->SetPoint(iDE,DEidErr[iDE],(Double_t)fAlign->GetParError(3*iDE+3));
  }

  // Post final data. Write histo list to a file with option "RECREATE"
  PostData(0,fList);

  // Re Align
  AliMUONGeometryTransformer *newTransform = fAlign->ReAlign(fTransform,fParameters,true); 
  newTransform->WriteTransformations("transform2ReAlign.dat");
  
  // Generate realigned data in local cdb
  const TClonesArray* array = newTransform->GetMisAlignmentData();

  // 100 mum residual resolution for chamber misalignments?
  fAlign->SetAlignmentResolution(array,-1,0.01,0.01,0.004,0.003);
   
  // CDB manager
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage(fDefaultOCDB.Data());
  cdbManager->SetSpecificStorage("MUON/Align/Data",fMisAlignOCDB.Data());
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON alignment objects with residual misalignment");
  AliCDBId id("MUON/Align/Data", 0, AliCDBRunRange::Infinity()); 
  cdbManager->Put(const_cast<TClonesArray*>(array), id, cdbData);

}

//________________________________________________________________________
void AliMUONAlignmentTask::Terminate(const Option_t*)
{
  /// Called once per task on the client machine at the end of the analysis.

}

//-----------------------------------------------------------------------
void AliMUONAlignmentTask::Prepare(const char* geoFilename, const char* defaultOCDB, const char* misAlignOCDB)
{
  /// Set the geometry, the magnetic field, the mapping and the reconstruction parameters
      
  // Load mapping
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(defaultOCDB);
  man->SetSpecificStorage("MUON/Align/Data",misAlignOCDB);
  man->Print();
  man->SetRun(fLastRunNumber);
  if ( ! AliMpCDB::LoadDDLStore() ) {
    Error("MUONRefit","Could not access mapping from OCDB !");
    exit(-1);
  }

  // Import TGeo geometry (needed by AliMUONTrackExtrap::ExtrapToVertex)
  if (!gGeoManager) {
    AliGeomManager::LoadGeometry(geoFilename);
    if (!gGeoManager) {
      Error("AliMUONReAlignTask", "getting geometry from file %s failed", "generated/galice.root");
      return;
    }
  }

  // set mag field
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    printf("Loading field map...\n");
    AliGRPManager *grpMan = new AliGRPManager();
    grpMan->ReadGRPEntry();
    grpMan->SetMagField();
    delete grpMan;
  }
  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField();
  
}
