/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TList.h>
#include <TChain.h>
#include <TH1F.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCuts.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSECleanupVertexingHF.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSECleanupVertexingHF)


//________________________________________________________________________
AliAnalysisTaskSECleanupVertexingHF::AliAnalysisTaskSECleanupVertexingHF():
AliAnalysisTaskSE()//,
//fNentries(0)
{
  // Default constructor

}
//_______________________________________________________
AliAnalysisTaskSECleanupVertexingHF::AliAnalysisTaskSECleanupVertexingHF(const char *name):
  AliAnalysisTaskSE(name)//,
  //fNentries(0)
{
  // Default constructor
  //DefineInput(0,TChain::Class());
  // DefineOutput(0,TH1F::Class());  //My private output
}
//_______________________________________________________
AliAnalysisTaskSECleanupVertexingHF::~AliAnalysisTaskSECleanupVertexingHF()
{
  //if (fNentries){
  //  delete fNentries;
  //  fNentries = 0;
  // }
}
//________________________________________________________
void AliAnalysisTaskSECleanupVertexingHF::UserCreateOutputObjects()
{

  //const char* nameoutput=GetOutputSlot(1)->GetContainer()->GetName();
  //fNentries=new TH1F(nameoutput, "Integral(1,2) = number of AODs ", 7,-0.5,6.5);

  //fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  //fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  //PostData(1,fNentries);

}
//________________________________________________________________________
void AliAnalysisTaskSECleanupVertexingHF::UserExec(Option_t */*option*/)
{

  //3prong
  TClonesArray *array3Prong = 0;
  //2prong D0 TClones Array
  TClonesArray *inputArrayD0=0;
  //DStar TClonesArray 
  TClonesArray *arrayDStartoD0pi=0;
  //Cascade TClonesArray
   TClonesArray *arrayCascade=0;

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      inputArrayD0=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      arrayDStartoD0pi=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar"); 
      arrayCascade=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
    }
  } else if(aod) {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    inputArrayD0=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    arrayDStartoD0pi=(TClonesArray*)aod->GetList()->FindObject("Dstar");  
    arrayCascade=(TClonesArray*)aod->GetList()->FindObject("CascadesHF");
  }
  if(!aod){
    printf("AliAnalysisTaskSECleanupVertexingHF::UserExec: aod not found!\n");
    return;
  }
  if(!array3Prong || !inputArrayD0 || !arrayDStartoD0pi) {
   printf("AliAnalysisTaskSECleanupVertexingHF::UserExec: an input branch not found!\n");
  }
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  // DStar
  for(Int_t iDStartoD0pi = 0; iDStartoD0pi<arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) { // D* candidates and D0 from D*
    AliAODRecoCascadeHF* recdstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
    if(!recdstarD0pi)continue;
    if(recdstarD0pi->GetIsFilled()==2){ 
    AliAODVertex *vtxDS = (AliAODVertex*)recdstarD0pi->GetSecondaryVtx();
    if(vtxDS) {delete vtxDS; vtxDS=0;}
    }else continue;
   }

  // D0toKpi
  for (Int_t iD0toKpi = 0; iD0toKpi < inputArrayD0->GetEntriesFast(); iD0toKpi++) {
    AliAODRecoDecayHF2Prong *recd = (AliAODRecoDecayHF2Prong*)inputArrayD0->UncheckedAt(iD0toKpi);
    if(!recd)continue;
    if(recd->GetIsFilled()==2){
    AliAODVertex *vtx2 = (AliAODVertex*)recd->GetSecondaryVtx();
    if(vtx2){delete vtx2;vtx2=0;}
  }else continue;
  }
  //3Prong
  for (Int_t i3prong = 0; i3prong < array3Prong->GetEntriesFast(); i3prong++) {
    AliAODRecoDecayHF3Prong *recd3 = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3prong);
    if(!recd3)continue;
    if(recd3->GetIsFilled()==2){
    AliAODVertex *vtx3 = (AliAODVertex*)recd3->GetSecondaryVtx();
    if(vtx3){delete vtx3;vtx3=0;}
  }else continue;
   }
  //Cascade
  if(arrayCascade){
  for(Int_t iCasc=0; iCasc<arrayCascade->GetEntriesFast(); iCasc++){
  AliAODRecoCascadeHF* recCasc = dynamic_cast<AliAODRecoCascadeHF*>(arrayCascade->At(iCasc));
  if(!recCasc)continue;
   if(recCasc->GetIsFilled()==2){
   AliAODVertex *vtxCasc = (AliAODVertex*)recCasc->GetSecondaryVtx();
   if(vtxCasc){delete vtxCasc; vtxCasc=0;}
   }else continue;
  }
  }

  return;
}
//_________________________
void AliAnalysisTaskSECleanupVertexingHF::Terminate(Option_t *option)
{
  // Terminate analysis
  //
  //fNentries = dynamic_cast<TH1F*>(GetOutputData(0));
  //if(!fNentries){
  //  printf("ERROR: fNEntries not available\n");
  return;
}



