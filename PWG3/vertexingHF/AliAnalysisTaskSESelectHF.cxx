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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the selection of heavy flavor
// decay candidates and creation a stand-alone AOD.
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>

#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSESelectHF.h"

ClassImp(AliAnalysisTaskSESelectHF)


//________________________________________________________________________
AliAnalysisTaskSESelectHF::AliAnalysisTaskSESelectHF():
AliAnalysisTaskSE(),
fVerticesHFTClArr(0),
fD0toKpiTClArr(0)
{
  // Default constructor
  SetD0toKpiCuts();
}

//________________________________________________________________________
AliAnalysisTaskSESelectHF::AliAnalysisTaskSESelectHF(const char *name):
AliAnalysisTaskSE(name),
fVerticesHFTClArr(0),
fD0toKpiTClArr(0)
{
  // Default constructor
  SetD0toKpiCuts();
}

//________________________________________________________________________
AliAnalysisTaskSESelectHF::~AliAnalysisTaskSESelectHF()
{
  // Destructor
}  

//________________________________________________________________________
void AliAnalysisTaskSESelectHF::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSESelectHF::Init() \n");

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

  // load D0->Kpi candidates                                                   
  TClonesArray *inputArrayD0toKpi =
    (TClonesArray*)aodIn->GetList()->FindObject("D0toKpi");
  if(!inputArrayD0toKpi) {
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
    Int_t okD0=0,okD0bar=0; 
    if(dIn->SelectD0(fD0toKpiCuts,okD0,okD0bar)) {
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
    }
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

//________________________________________________________________________
void AliAnalysisTaskSESelectHF::SetD0toKpiCuts(Double_t cut0,Double_t cut1,
				   Double_t cut2,Double_t cut3,Double_t cut4,
				   Double_t cut5,Double_t cut6,
				   Double_t cut7,Double_t cut8) 
{
  // Set the cuts for D0 selection
  // cuts[0] = inv. mass half width [GeV]   
  // cuts[1] = dca [cm]
  // cuts[2] = cosThetaStar 
  // cuts[3] = pTK [GeV/c]
  // cuts[4] = pTPi [GeV/c]
  // cuts[5] = d0K [cm]   upper limit!
  // cuts[6] = d0Pi [cm]  upper limit!
  // cuts[7] = d0d0 [cm^2]
  // cuts[8] = cosThetaPoint

  fD0toKpiCuts[0] = cut0;
  fD0toKpiCuts[1] = cut1;
  fD0toKpiCuts[2] = cut2;
  fD0toKpiCuts[3] = cut3;
  fD0toKpiCuts[4] = cut4;
  fD0toKpiCuts[5] = cut5;
  fD0toKpiCuts[6] = cut6;
  fD0toKpiCuts[7] = cut7;
  fD0toKpiCuts[8] = cut8;

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF::SetD0toKpiCuts(const Double_t cuts[9]) 
{
  // Set the cuts for D0 selection
  // cuts[0] = inv. mass half width [GeV]   
  // cuts[1] = dca [cm]
  // cuts[2] = cosThetaStar 
  // cuts[3] = pTK [GeV/c]
  // cuts[4] = pTPi [GeV/c]
  // cuts[5] = d0K [cm]   upper limit!
  // cuts[6] = d0Pi [cm]  upper limit!
  // cuts[7] = d0d0 [cm^2]
  // cuts[8] = cosThetaPoint

  for(Int_t i=0; i<9; i++) fD0toKpiCuts[i] = cuts[i];

  return;
}
