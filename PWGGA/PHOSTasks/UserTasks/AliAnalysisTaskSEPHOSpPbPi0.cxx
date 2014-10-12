/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the gamma and pi0 from pPb collision analysis
//
// Author: H-S. Zhu, hongsheng.zhu@cern.ch
//                   hszhu@iopp.ccnu.edu.cn
///////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TH2I.h>
#include <TList.h>
#include <TMath.h>
#include <TArray.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TGeoManager.h>
#include <TClonesArray.h>

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODHeader.h"
#include "AliESDHeader.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliAODCaloCells.h"
#include "AliESDCaloCells.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSGeoUtils.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSAodCluster.h"
#include "AliPHOSEsdCluster.h"
#include "AliOADBContainer.h"
#include "AliPHOSpPbPi0Header.h"
#include "AliCaloClusterInfo.h"
#include "AliAnalysisTaskSEPHOSpPbPi0.h"

ClassImp(AliAnalysisTaskSEPHOSpPbPi0)

Bool_t   AliAnalysisTaskSEPHOSpPbPi0::fgRemovePileup    = kFALSE;
Bool_t   AliAnalysisTaskSEPHOSpPbPi0::fgUseFiducialCut  = kFALSE;
Bool_t   AliAnalysisTaskSEPHOSpPbPi0::fgUseTOFCut       = kFALSE;
Double_t AliAnalysisTaskSEPHOSpPbPi0::fgCuts[5]         = { 0.3, 2., 0.2, 2.5, 1e-7 };

//________________________________________________________________________
AliAnalysisTaskSEPHOSpPbPi0::AliAnalysisTaskSEPHOSpPbPi0():
  AliAnalysisTaskSE(), fIsMC(kFALSE), fCentralityBin(10), fBufferSize(10), fRunNumber(-1),
  fPHOSGeo(0), fOutputListQA(0), fOutputListRD(0), fOutputListMC(0), fHeader(0), fCaloClArr(0)
{
  //
  // Default constructor
  //
  for (Int_t i=0; i<10; i++) { for (Int_t j=0; j<10; j++) fEventList[i][j] = 0; }
}
//________________________________________________________________________
AliAnalysisTaskSEPHOSpPbPi0::AliAnalysisTaskSEPHOSpPbPi0(const char *name):
  AliAnalysisTaskSE(name), fIsMC(kFALSE), fCentralityBin(10), fBufferSize(10), fRunNumber(-1),
  fPHOSGeo(0), fOutputListQA(0), fOutputListRD(0), fOutputListMC(0), fHeader(0), fCaloClArr(0)
{
  // Constructor
  for (Int_t i=0; i<10; i++) { for (Int_t j=0; j<10; j++) fEventList[i][j] = 0; }

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEPHOSpPbPi0::~AliAnalysisTaskSEPHOSpPbPi0()
{
  //
  // Default destructor
  //
  if (fOutputListQA)  { delete fOutputListQA;  fOutputListQA = NULL; }
  if (fOutputListRD)  { delete fOutputListRD;  fOutputListRD = NULL; }
  if (fOutputListMC)  { delete fOutputListMC;  fOutputListMC = NULL; }
  if (fHeader)        { delete fHeader;        fHeader       = NULL; }
  if (fCaloClArr)     { delete fCaloClArr;     fCaloClArr    = NULL; }
}
//________________________________________________________________________
void AliAnalysisTaskSEPHOSpPbPi0::UserCreateOutputObjects()
{
  // Create the output container

  AliPHOSpPbPi0Header::SetIsMC(fIsMC);
  AliPHOSpPbPi0Header::SetUseFiducialCut(fgUseFiducialCut);

  if (!fHeader)       fHeader       = new AliPHOSpPbPi0Header();
  if (!fCaloClArr)    fCaloClArr    = new TClonesArray("AliCaloClusterInfo", 0);
  if (!fOutputListQA) fOutputListQA = new TList();   fOutputListQA->SetOwner(kTRUE);
  if (!fOutputListRD) fOutputListRD = new TList();   fOutputListRD->SetOwner(kTRUE);
  if (!fOutputListMC) fOutputListMC = new TList();   fOutputListMC->SetOwner(kTRUE);

  fHeader->SetNCent(fCentralityBin.GetSize()-1);
  fHeader->CreateHistograms(fOutputListQA, fOutputListRD, fOutputListMC);

  // Post output data.
  PostData(1, fOutputListQA);
  PostData(2, fOutputListRD);
  if (fIsMC) PostData(3, fOutputListMC);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEPHOSpPbPi0::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (fIsMC) {
    if (MCEvent()) {
      if (MCEvent()->GetNumberOfTracks()<=0)
           { AliError("MC event not found. Nothing done!"); return; }
    } else { AliError("MC event not found. Nothing done!"); return; }
  }

  AliAODEvent *aod   = 0x0;
  AliESDEvent *esd   = 0x0;

  if (((TString)fInputEvent->IsA()->GetName())=="AliAODEvent") {
    aod = dynamic_cast<AliAODEvent*>(fInputEvent);
    if (!aod) { AliError("AOD event not found. Nothing done!");   return; }
    if (!fIsMC && (aod->GetHeader()->GetEventType()!=7))          return; // check event type; should be PHYSICS = 7 for data and 0 for MC
  } else {
    esd = dynamic_cast<AliESDEvent*>(fInputEvent);
    if (!esd) { AliError("ESD event not found. Nothing done!");   return; }
    if (!fIsMC && (esd->GetHeader()->GetEventType()!=7))          return; // check event type; should be PHYSICS = 7 for data and 0 for MC
  }

  // Fill Event info
  fHeader->SetEventInfo(fInputHandler);
  fHeader->FillHistosEvent(fOutputListQA);

  // PHOS Geometry and Misalignment initialization at the first time it runs
  if(fRunNumber != fInputEvent->GetRunNumber()) {
    fRunNumber = fInputEvent->GetRunNumber();
    PHOSInitialize();
  }

  // Event Selection
  if (!fHeader->IsSelected())                     return;
  if (fgRemovePileup && fHeader->IsPileup())      return;

  // Fill PHOS cells QA histograms
  fHeader->FillHistosCaloCellsQA(fOutputListQA, fInputEvent->GetPHOSCells(), fPHOSGeo);

  // Fill PHOS cluster Clones Array
  FillCaloClusterInfo();
  aod = 0x0;
  esd = 0x0;

  // vertex bining and centrality bining
  Int_t zvtx = (Int_t)((fHeader->Vz() + 10.)/2.);   if (zvtx<0) zvtx = 0;   if (zvtx>9) zvtx = 9;
  Int_t cent = TMath::BinarySearch<Float_t>(fCentralityBin.GetSize()-1, fCentralityBin.GetArray(), fHeader->Centrality());
  if (!fEventList[zvtx][cent]) fEventList[zvtx][cent] = new TList(); 
  TList *eventList = fEventList[zvtx][cent];

  // Fill MC info JUST FOR ESD
  AliStack *stack = 0x0;
  if (fIsMC) {
    if (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()) {
      if (static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent())
        stack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
    }
    fHeader->FillHistosMC(fOutputListMC, stack, fCaloClArr, fPHOSGeo, cent);
  }

  if (!fCaloClArr->GetEntriesFast()) return;

  // Fill cluster histograms
  fHeader->FillHistosCaloCluster(fOutputListQA, fCaloClArr, cent);

  // Fill pi0 histograms
  fHeader->FillHistosPi0(fOutputListRD, fCaloClArr, cent);

  // Fill mixed pi0 histograms
  fHeader->FillHistosMixPi0(fOutputListRD, fCaloClArr, eventList, cent);

  // Fill event list for mixing
  if (fCaloClArr->GetEntriesFast()>0) {
    eventList->AddFirst(fCaloClArr);   fCaloClArr = 0x0;
    if (eventList->GetSize()>fBufferSize[cent]) { // Remove redundant events
      TClonesArray *tmp = static_cast<TClonesArray*>(eventList->Last());
      eventList->RemoveLast();
      delete tmp;
    }
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEPHOSpPbPi0::Terminate(Option_t *) 
{
  // Terminate analysis

  // add the correction matrix

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEPHOSpPbPi0::PHOSInitialize()
{
  // Initialize PHOS Geometry ,misalignment and calibration 

  fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP");
  AliOADBContainer geomContainer("phosGeo");
  geomContainer.InitFromFile("$ALICE_ROOT/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
  TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
  for (Int_t mod=0; mod<5; mod++) {
    if (!matrixes->At(mod)) continue;
    else fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod);
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEPHOSpPbPi0::FillCaloClusterInfo()
{
  // Fill calo cluster info
  if (fCaloClArr) fCaloClArr->Clear();
  else fCaloClArr = new TClonesArray("AliCaloClusterInfo", 0);

  Int_t    nclsts       = fInputEvent->GetNumberOfCaloClusters();
  Int_t    countN       = fCaloClArr->GetEntriesFast();
  Int_t    relID[4]     = {0,0,0,0 }; // module = relID[0]; cellX = relID[2]; cellZ = relID[3];
  Float_t  position[3]  = {0.,0.,0.};
  Double_t vtx[3]       = {0.,0.,0.};   fHeader->GetXYZ(vtx);   TVector3 vtxVector(vtx); 

  TClonesArray &caloRef = *fCaloClArr;
  TLorentzVector momentum;
  AliVCluster        *clust       = 0x0;
  AliCaloClusterInfo *caloCluster = 0x0;
  for (Int_t iclst=0; iclst<nclsts; iclst++) {  // loop over all clusters
    clust = fInputEvent->GetCaloCluster(iclst);
    if (!(clust && clust->IsPHOS() && clust->E()>fgCuts[0]))                 { clust=0; continue; }
    if (!(clust->GetNCells()>(Int_t)fgCuts[1] && clust->GetM02()>fgCuts[2])) { clust=0; continue; } // To remove exotic clusters
    if (!(clust->GetDistanceToBadChannel()>fgCuts[3]))                       { clust=0; continue; }
    if (fgUseTOFCut && !(TMath::Abs(clust->GetTOF())<fgCuts[4]))             { clust=0; continue; } // remove clusters from pileup or same benches
    clust->GetPosition(position); TVector3 global(position); fPHOSGeo->GlobalPos2RelId(global,relID);
    if (relID[0] == 2)                                                       { clust=0; continue; } // !remove module 2

    caloCluster = new AliCaloClusterInfo(clust, relID);
//  if (fgUseFiducialCut && !caloCluster->IsInFiducialRegion(relID[2], relID[3])) { delete caloCluster; caloCluster=0; clust=0; continue; } // Fiducial cut 

    clust->GetMomentum(momentum, vtx);
    caloCluster->SetLorentzVector(momentum);
    if (fIsMC) caloCluster->SetLabels(clust->GetNLabels(), clust->GetLabels());

    clust = 0;

    new(caloRef[countN++]) AliCaloClusterInfo(*caloCluster);
    delete caloCluster; caloCluster=0;
  }  // end loop of all clusters
}
