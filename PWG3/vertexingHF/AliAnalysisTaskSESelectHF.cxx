/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the selection of heavy flavor
// decay candidates and creation a stand-alone AOD.
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSESelectHF.h"

ClassImp(AliAnalysisTaskSESelectHF)


//________________________________________________________________________
AliAnalysisTaskSESelectHF::AliAnalysisTaskSESelectHF():
AliAnalysisTaskSE(),
fVerticesHFTClArr(0),
fD0toKpiTClArr(0),
fVHF(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSESelectHF::AliAnalysisTaskSESelectHF(const char *name):
AliAnalysisTaskSE(name),
fVerticesHFTClArr(0),
fD0toKpiTClArr(0),
fVHF(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSESelectHF::~AliAnalysisTaskSESelectHF()
{
  // Destructor

  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  }

}  

//________________________________________________________________________
void AliAnalysisTaskSESelectHF::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSESelectHF::Init() \n");

  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  fVHF->PrintStatus();

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSESelectHF::UserCreateOutputObjects() \n");

  fVerticesHFTClArr = new TClonesArray("AliAODVertex", 0);
  fVerticesHFTClArr->SetName("VerticesHF");
  AddAODBranch("TClonesArray", &fVerticesHFTClArr);

  fD0toKpiTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
  fD0toKpiTClArr->SetName("D0toKpi");
  AddAODBranch("TClonesArray", &fD0toKpiTClArr);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates selection and histograms
  
  AliAODEvent *aodIn = dynamic_cast<AliAODEvent*> (InputEvent());

  TClonesArray *inputArrayD0toKpi = 0;

  if(!aodIn && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodIn = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      // load D0 candidates                                                   
      inputArrayD0toKpi=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
    }
  } else if(aodIn) {
    // load D0 candidates                                                   
    inputArrayD0toKpi=(TClonesArray*)aodIn->GetList()->FindObject("D0toKpi");
  }

  if(!inputArrayD0toKpi || !aodIn) {
    printf("AliAnalysisTaskSESelectHF::UserExec: D0toKpi branch not found!\n");
    return;
  }

  //print event info
  //aodIn->GetHeader()->Print();
  
  // primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aodIn->GetPrimaryVertex();
  //vtx1->Print();
    
  // make trkIDtoEntry register (temporary)
  Int_t trkIDtoEntry[100000];
  for(Int_t it=0;it<aodIn->GetNumberOfTracks();it++) {
    AliAODTrack *track = aodIn->GetTrack(it);
    trkIDtoEntry[track->GetID()]=it;
  }

  Int_t iOutVerticesHF=0,iOutD0toKpi=0;
  fVerticesHFTClArr->Delete();
  iOutVerticesHF = fVerticesHFTClArr->GetEntriesFast();
  TClonesArray &verticesHFRef = *fVerticesHFTClArr;
  fD0toKpiTClArr->Delete();
  iOutD0toKpi = fD0toKpiTClArr->GetEntriesFast();
  TClonesArray &aodD0toKpiRef = *fD0toKpiTClArr;


  // loop over D0->Kpi candidates
  Int_t nInD0toKpi = inputArrayD0toKpi->GetEntriesFast();
  printf("Number of D0->Kpi: %d\n",nInD0toKpi);
  
  for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
    AliAODRecoDecayHF2Prong *dIn = (AliAODRecoDecayHF2Prong*)inputArrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!dIn->GetOwnPrimaryVtx()) {
      dIn->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    
    //Int_t okD0=0,okD0bar=0; 
    //if(dIn->SelectD0(fVHF->GetD0toKpiCuts(),okD0,okD0bar)) {
      // get daughter AOD tracks
      AliAODTrack *trk0 = (AliAODTrack*)dIn->GetDaughter(0);
      AliAODTrack *trk1 = (AliAODTrack*)dIn->GetDaughter(1);
      if(!trk0 || !trk1) {
	trk0=aodIn->GetTrack(trkIDtoEntry[dIn->GetProngID(0)]);
	trk1=aodIn->GetTrack(trkIDtoEntry[dIn->GetProngID(1)]);
      }
      printf("pt of positive track: %f\n",trk0->Pt());
      printf("pt of negative track: %f\n",trk1->Pt());
      // HERE ONE COULD RECALCULATE THE VERTEX USING THE KF PACKAGE

      // clone candidate for output AOD
      AliAODVertex *v = new(verticesHFRef[iOutVerticesHF++]) 
	AliAODVertex(*(dIn->GetSecondaryVtx()));
      AliAODRecoDecayHF2Prong *dOut=new(aodD0toKpiRef[iOutD0toKpi++]) 
	AliAODRecoDecayHF2Prong(*dIn);
      dOut->SetSecondaryVtx(v);
      dOut->SetOwnPrimaryVtx((AliAODVertex*)((dIn->GetOwnPrimaryVtx())->Clone()));
      v->SetParent(dOut);
    
    
    if(unsetvtx) dIn->UnsetOwnPrimaryVtx();
  } // end loop on D0->Kpi

  printf("Number of selected D0->Kpi: %d\n",iOutD0toKpi);


  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSESelectHF: Terminate() \n");
}

