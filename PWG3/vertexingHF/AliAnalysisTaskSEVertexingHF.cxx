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
// AliAnalysisTaskSE for the reconstruction of heavy flavor
// decays, using the class AliAnalysisVertexingHF.
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TString.h>

#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSEVertexingHF.h"

ClassImp(AliAnalysisTaskSEVertexingHF)


//________________________________________________________________________
AliAnalysisTaskSEVertexingHF::AliAnalysisTaskSEVertexingHF():
AliAnalysisTaskSE(),
fVHF(0),
fListOfCuts(0),
fDeltaAODFileName("AliAOD.VertexingHF.root"),
fVerticesHFTClArr(0),
fD0toKpiTClArr(0),
fJPSItoEleTClArr(0),
fCharm3ProngTClArr(0),
fCharm4ProngTClArr(0),
fDstarTClArr(0),
fLikeSign2ProngTClArr(0),
fLikeSign3ProngTClArr(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEVertexingHF::AliAnalysisTaskSEVertexingHF(const char *name):
AliAnalysisTaskSE(name),
fVHF(0),
fListOfCuts(0),
fDeltaAODFileName("AliAOD.VertexingHF.root"),
fVerticesHFTClArr(0),
fD0toKpiTClArr(0),
fJPSItoEleTClArr(0),
fCharm3ProngTClArr(0),
fCharm4ProngTClArr(0),
fDstarTClArr(0),
fLikeSign2ProngTClArr(0),
fLikeSign3ProngTClArr(0)
{
  // Standard constructor

  DefineOutput(1,TList::Class()); // analysis cuts
}

//________________________________________________________________________
AliAnalysisTaskSEVertexingHF::~AliAnalysisTaskSEVertexingHF()
{
  // Destructor

  if(fListOfCuts) {
    delete fListOfCuts;
    fListOfCuts=NULL;
  }

}  

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHF::Init()
{
  // Initialization
  // Instanciates vHF and loads its parameters

  if(fDebug > 1) printf("AnalysisTaskSEVertexingHF::Init() \n");

  if(gROOT->LoadMacro("ConfigVertexingHF.C")) {
    printf("AnalysisTaskSEVertexingHF::Init() \n Using $ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C\n");
    gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C");
  }

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  fVHF->PrintStatus();


  // write the objects AliRDHFCuts to a list to store in the output

  fListOfCuts = fVHF->FillListOfCuts();

  PostData(1,fListOfCuts);

  AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile(fDeltaAODFileName.Data());

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHF::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEVertexingHF::UserCreateOutPutData() \n");
  // Support both the case when the AOD + deltaAOD are produced in an ESD
  // analysis or if the deltaAOD is produced on an analysis on AOD's. (A.G. 27/04/09)
  if(!AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()) {
    Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
    return;
  }   
  TString filename = fDeltaAODFileName;
  // When running on standard AOD to produce deltas, IsStandardAOD is never set,
  // If AODEvent is NULL, new branches have to be added to the new file(s) (A.G. 15/01/10)
  if(!IsStandardAOD() && AODEvent()) filename = "";
  if(!fVHF) {
    printf("AnalysisTaskSEVertexingHF::UserCreateOutPutData() \n ERROR! no fvHF!\n");
    return;
  }

  fVerticesHFTClArr = new TClonesArray("AliAODVertex", 0);
  fVerticesHFTClArr->SetName("VerticesHF");
  AddAODBranch("TClonesArray", &fVerticesHFTClArr, filename);

  if(fVHF->GetD0toKpi()) {
    fD0toKpiTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fD0toKpiTClArr->SetName("D0toKpi");
    AddAODBranch("TClonesArray", &fD0toKpiTClArr, filename);
  }

  if(fVHF->GetJPSItoEle()) {
    fJPSItoEleTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fJPSItoEleTClArr->SetName("JPSItoEle");
    AddAODBranch("TClonesArray", &fJPSItoEleTClArr, filename);
  }

  if(fVHF->Get3Prong()) {
    fCharm3ProngTClArr = new TClonesArray("AliAODRecoDecayHF3Prong", 0);
    fCharm3ProngTClArr->SetName("Charm3Prong");
    AddAODBranch("TClonesArray", &fCharm3ProngTClArr, filename);
  }

  if(fVHF->Get4Prong()) {
    fCharm4ProngTClArr = new TClonesArray("AliAODRecoDecayHF4Prong", 0);
    fCharm4ProngTClArr->SetName("Charm4Prong");
    AddAODBranch("TClonesArray", &fCharm4ProngTClArr, filename);
  }

  if(fVHF->GetDstar()) {
    fDstarTClArr = new TClonesArray("AliAODRecoCascadeHF", 0);
    fDstarTClArr->SetName("Dstar");
    AddAODBranch("TClonesArray", &fDstarTClArr, filename);
  }

  if(fVHF->GetLikeSign()) {                      
    fLikeSign2ProngTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fLikeSign2ProngTClArr->SetName("LikeSign2Prong");
    AddAODBranch("TClonesArray", &fLikeSign2ProngTClArr, filename);
  }

  if(fVHF->GetLikeSign() && fVHF->Get3Prong()) {                      
    fLikeSign3ProngTClArr = new TClonesArray("AliAODRecoDecayHF3Prong", 0);
    fLikeSign3ProngTClArr->SetName("LikeSign3Prong");
    AddAODBranch("TClonesArray", &fLikeSign3ProngTClArr, filename);
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHF::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor vertexing
  
  AliVEvent *event = dynamic_cast<AliVEvent*> (InputEvent());
  // In case there is an AOD handler writing a standard AOD, use the AOD 
  // event in memory rather than the input (ESD) event. (A.G. 27/04/09)
  if (AODEvent() && IsStandardAOD()) event = dynamic_cast<AliVEvent*> (AODEvent());

  // heavy flavor vertexing
  fVHF->FindCandidates(event,
		       fVerticesHFTClArr,
		       fD0toKpiTClArr,
		       fJPSItoEleTClArr,
		       fCharm3ProngTClArr,
		       fCharm4ProngTClArr,
		       fDstarTClArr,
                       fLikeSign2ProngTClArr,
                       fLikeSign3ProngTClArr);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHF::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEVertexingHF: Terminate() \n");
}
