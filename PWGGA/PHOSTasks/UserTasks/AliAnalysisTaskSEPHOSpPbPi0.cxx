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

#include <TH2I.h>
#include <TList.h>
#include <TMath.h>
#include <TArray.h>
#include <TClonesArray.h>

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDHeader.h"
#include "AliAODHeader.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliPHOSGeoUtils.h"
#include "AliPHOSGeometry.h"
#include "AliOADBContainer.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliPHOSpPbPi0Header.h"
#include "AliCaloClusterInfo.h"
#include "AliAnalysisTaskSEPHOSpPbPi0.h"

ClassImp(AliAnalysisTaskSEPHOSpPbPi0)

Int_t    AliAnalysisTaskSEPHOSpPbPi0::fgMinNCells        = 2;
Double_t AliAnalysisTaskSEPHOSpPbPi0::fgMinClusterEnergy = 0.3;
Double_t AliAnalysisTaskSEPHOSpPbPi0::fgMinM02           = 0.2;
Double_t AliAnalysisTaskSEPHOSpPbPi0::fgMinDistToBad     = 2.5;

//________________________________________________________________________
AliAnalysisTaskSEPHOSpPbPi0::AliAnalysisTaskSEPHOSpPbPi0():
  AliAnalysisTaskSE(), fIsMC(kFALSE), fCentralityBin(10), fBufferSize(10),
  fRunNumber(-1), fPHOSGeo(0), fList(0), fHeader(0), fCaloClArr(0)
{
  //
  // Default constructor
  //
  for (Int_t i=0; i<10; i++) { for (Int_t j=0; j<10; j++) fEventList[i][j] = 0; }

  // Set bad channel map
  for (Int_t i=0; i<5; i++)
    fPHOSBadMap[i] = new TH2I(Form("PHOS_BadMap_mod%d", i), Form("PHOS_BadMap_mod%d", i), 64, 0., 64., 56, 0., 56.);
}
//________________________________________________________________________
AliAnalysisTaskSEPHOSpPbPi0::AliAnalysisTaskSEPHOSpPbPi0(const char *name):
  AliAnalysisTaskSE(name), fIsMC(kFALSE), fCentralityBin(10), fBufferSize(10),
  fRunNumber(-1), fPHOSGeo(0), fList(0), fHeader(0), fCaloClArr(0)
{
  // Constructor
  for (Int_t i=0; i<10; i++) { for (Int_t j=0; j<10; j++) fEventList[i][j] = 0; }

  // Set bad channel map
  for (Int_t i=0; i<5; i++)
    fPHOSBadMap[i] = new TH2I(Form("PHOS_BadMap_mod%d", i), Form("PHOS_BadMap_mod%d", i), 64, 0., 64., 56, 0., 56.);

  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEPHOSpPbPi0::~AliAnalysisTaskSEPHOSpPbPi0()
{
  //
  // Default destructor
  //
  if (fList)        { delete fList;        fList      = NULL; }
  if (fHeader)      { delete fHeader;      fHeader    = NULL; }
  if (fCaloClArr)   { delete fCaloClArr;   fCaloClArr = NULL; }
}

//________________________________________________________________________
void AliAnalysisTaskSEPHOSpPbPi0::UserCreateOutputObjects()
{
  // Create the output container
  // Initialize the PHOS geometry
//fPHOSGeo = new AliPHOSGeoUtils("PHOSGeo");
//fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP");

  if (!fHeader)    fHeader    = new AliPHOSpPbPi0Header();
  if (!fCaloClArr) fCaloClArr = new TClonesArray("AliCaloClusterInfo", 0);
  if (!fList)      fList      = new TList();

  fHeader->SetIsMC(fIsMC);
  fHeader->SetNCent(fCentralityBin.GetSize());
  fHeader->CreateHistograms(fList);

  // Post output data.
  PostData(1, fList);
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

  Int_t       nclsts = 0;
  AliAODEvent *aod   = 0;
  AliESDEvent *esd   = 0;

  if (((TString)fInputEvent->IsA()->GetName())=="AliAODEvent") {
    aod = dynamic_cast<AliAODEvent*>(fInputEvent);
    if (!aod) { AliError("AOD event not found. Nothing done!");   return; }
    if (!fIsMC && (aod->GetHeader()->GetEventType()!=7))          return; // check event type; should be PHYSICS = 7 for data and 0 for MC
    if (!fHeader->IspAVertexOK(aod))                              return; // check p-A collision vertex
    nclsts = aod->GetNumberOfCaloClusters();   if (nclsts<1)      return;
  } else {
    esd = dynamic_cast<AliESDEvent*>(fInputEvent);
    if (!esd) { AliError("ESD event not found. Nothing done!");   return; }
    if (!fIsMC && (esd->GetHeader()->GetEventType()!=7))          return; // check event type; should be PHYSICS = 7 for data and 0 for MC
    if (!fHeader->IspAVertexOK(esd))                              return; // check p-A collision vertex
    nclsts = esd->GetNumberOfCaloClusters();   if (nclsts<1)      return;
  }

  // Fill Event info
  fHeader->SetEventInfo(fInputHandler);
  fHeader->FillHistosEvnH(fList);

  // PHOS Geometry and Misalignment initialization at the first time it runs
  if(fRunNumber != fInputEvent->GetRunNumber()) {
    fRunNumber = fInputEvent->GetRunNumber();
    PHOSInitialize(esd);
  }

  // Fill PHOS cells QA histograms
  fHeader->FillHistosCaloCellsQA(fList, fInputEvent->GetPHOSCells(), fPHOSGeo);
//aod = 0;

  // Fill PHOS cluster Clones Array
  FillCaloClusterInfo(nclsts, esd);

  if (!fCaloClArr->GetEntriesFast()) return;
  Int_t zvtx = (Int_t)((fHeader->Vz() + 10.)/2.);    if (zvtx<0) zvtx = 0; if (zvtx>9) zvtx = 9;
  Int_t cent = TMath::BinarySearch<Double_t>(fCentralityBin.GetSize()-1, fCentralityBin.GetArray(), fHeader->Centrality());
  if (!fEventList[zvtx][cent]) fEventList[zvtx][cent] = new TList(); 
  TList *eventList = fEventList[zvtx][cent];

  // Fill cluster histograms
  fHeader->FillHistosCaloCluster(fList, fCaloClArr, cent);

  // Fill pi0 histograms
  fHeader->FillHistosPi0(fList, fCaloClArr, cent);

  // Fill mixed pi0 histograms
  fHeader->FillHistosMixPi0(fList, fCaloClArr, eventList, cent);

  // Fill MC info
  AliStack *stack = 0;
  if (fIsMC) {
    if (esd) {
      if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
      if(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent())
        stack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
      }
      fHeader->FillHistosMC(fList, stack, fPHOSGeo, cent);
    } else
      fHeader->FillHistosMC(fList, MCEvent(), fPHOSGeo, cent);
  }
//esd = 0;

  // Fill event list for mixing
  if(fCaloClArr->GetEntriesFast()>0) {
    eventList->AddFirst(fCaloClArr);   fCaloClArr = 0;
    //fCaloClArr->Clear();
    if(eventList->GetSize()>fBufferSize[cent]) { // Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(eventList->Last()) ;
      eventList->RemoveLast();
      delete tmp ;
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
void AliAnalysisTaskSEPHOSpPbPi0::PHOSInitialize(AliESDEvent* const esd)
{
  // Initialize PHOS Geometry and misalignment
  // PHOS Geometry
  AliOADBContainer geomContainer("phosGeo");
  geomContainer.InitFromFile("$ALICE_ROOT/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
  TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
  fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP");
  for (Int_t mod=0; mod<5; mod++) {
    if (!matrixes->At(mod)) continue;
    else fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod);
  }

  // sets the PHOS Misalignment vertex if ESD
  if (esd) {
    for (Int_t mod=0; mod<5; mod++) {
      const TGeoHMatrix* modMatrix = fInputEvent->GetPHOSMatrix(mod);
      if (!modMatrix) continue;
      else fPHOSGeo->SetMisalMatrix(modMatrix, mod);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEPHOSpPbPi0::FillCaloClusterInfo(Int_t nclsts, AliESDEvent* const esd)
{
  // Fill calo cluster info
  if (fCaloClArr) fCaloClArr->Clear();
  else fCaloClArr = new TClonesArray("AliCaloClusterInfo", 0);

  TClonesArray       &caloRef     = *fCaloClArr;
  Int_t              countN       = 0;
  Int_t              relId[4]     = {0,0,0,0}; // module = relId[0]; cellX = relId[2]; cellZ = relId[3];
  Float_t            position[3]  = {0,0,0};
  Double_t           vtx[3]       = {0,0,0}; fHeader->GetXYZ(vtx);
  TLorentzVector     momentum;
  AliVCluster        *clust       = 0;
  AliCaloClusterInfo *caloCluster = 0;
  for (Int_t iclst=0; iclst<nclsts; iclst++) {  // loop over all clusters
    clust = fInputEvent->GetCaloCluster(iclst);
    if (!(clust && clust->IsPHOS() && clust->E()>fgMinClusterEnergy))     { clust=0; continue; }
    if (!(clust->GetNCells()>fgMinNCells && clust->GetM02()>fgMinM02))    { clust=0; continue; } // To remove exotic clusters
    if (!(clust->GetDistanceToBadChannel()>fgMinDistToBad))               { clust=0; continue; }
    clust->GetPosition(position); TVector3 global(position);
    fPHOSGeo->GlobalPos2RelId(global,relId);
    if (!IsGoodCaloCluster(relId[0], relId[2], relId[3]))                 { clust=0; continue; } // calo cluster selection 
    if (relId[0] == 2)                                                    { clust=0; continue; } // !remove module 2

    caloCluster = new AliCaloClusterInfo(clust, esd, fPHOSGeo, vtx);
//  if (esd && !(caloCluster->LorentzVector().E()>fgMinClusterEnergy)) { delete caloCluster; caloCluster=0; clust=0; continue; } // check again for ESD
    if (caloCluster->TestCPV(fHeader->MagneticField())) caloCluster->SetPIDBit(BIT(0));          // set CPV bit

    clust = 0;

    new(caloRef[countN++]) AliCaloClusterInfo(*caloCluster);

    delete caloCluster; caloCluster=0;
  }  // end loop of all clusters
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEPHOSpPbPi0::IsGoodCaloCluster(Int_t iMod, Int_t cellX, Int_t cellZ)
{
  // !Check whether this cluster is not in bad channel

  if (!(iMod>0 && iMod<5 && fPHOSBadMap[iMod]))          return kTRUE;   //No bad maps for this Module
  if (fPHOSBadMap[iMod]->GetBinContent(cellX,cellZ)>0)   return kFALSE;

  return kTRUE;
}
