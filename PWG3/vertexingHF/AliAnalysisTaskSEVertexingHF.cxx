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

#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEVertexingHF.h"

ClassImp(AliAnalysisTaskSEVertexingHF)


//________________________________________________________________________
AliAnalysisTaskSEVertexingHF::AliAnalysisTaskSEVertexingHF():
AliAnalysisTaskSE(),
fVHF(0),
fVerticesHFTClArr(0),
fD0toKpiTClArr(0),
fJPSItoEleTClArr(0),
fCharm3ProngTClArr(0),
fCharm4ProngTClArr(0),
fDstarTClArr(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEVertexingHF::AliAnalysisTaskSEVertexingHF(const char *name):
AliAnalysisTaskSE(name),
fVHF(0),
fVerticesHFTClArr(0),
fD0toKpiTClArr(0),
fJPSItoEleTClArr(0),
fCharm3ProngTClArr(0),
fCharm4ProngTClArr(0),
fDstarTClArr(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEVertexingHF::~AliAnalysisTaskSEVertexingHF()
{
  // Destructor
}  

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHF::Init()
{
  // Initialization
  // Instanciates vHF and loads its parameters

  if(fDebug > 1) printf("AnalysisTaskSEVertexingHF::Init() \n");

  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  fVHF->PrintStatus();

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHF::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEVertexingHF::UserCreateOutPutData() \n");

  if(!fVHF) {
    printf("AnalysisTaskSEVertexingHF::UserCreateOutPutData() \n ERROR! no fvHF!\n");
    return;
  }

  fVerticesHFTClArr = new TClonesArray("AliAODVertex", 0);
  fVerticesHFTClArr->SetName("VerticesHF");
  AddAODBranch("TClonesArray", &fVerticesHFTClArr);

  if(fVHF->GetD0toKpi()) {
    fD0toKpiTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fD0toKpiTClArr->SetName("D0toKpi");
    AddAODBranch("TClonesArray", &fD0toKpiTClArr);
  }

  if(fVHF->GetJPSItoEle()) {
    fJPSItoEleTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fJPSItoEleTClArr->SetName("JPSItoEle");
    AddAODBranch("TClonesArray", &fJPSItoEleTClArr);
  }

  if(fVHF->Get3Prong()) {
    fCharm3ProngTClArr = new TClonesArray("AliAODRecoDecayHF3Prong", 0);
    fCharm3ProngTClArr->SetName("Charm3Prong");
    AddAODBranch("TClonesArray", &fCharm3ProngTClArr);
  }

  if(fVHF->Get4Prong()) {
    fCharm4ProngTClArr = new TClonesArray("AliAODRecoDecayHF4Prong", 0);
    fCharm4ProngTClArr->SetName("Charm4Prong");
    AddAODBranch("TClonesArray", &fCharm4ProngTClArr);
  }

  if(fVHF->GetDstar()) {
    fDstarTClArr = new TClonesArray("AliAODRecoCascadeHF", 0);
    fDstarTClArr->SetName("Dstar");
    AddAODBranch("TClonesArray", &fDstarTClArr);
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHF::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor vertexing
  
  AliVEvent *event = dynamic_cast<AliVEvent*> (InputEvent());

  // heavy flavor vertexing
  fVHF->FindCandidates(event,
		       fVerticesHFTClArr,
		       fD0toKpiTClArr,
		       fJPSItoEleTClArr,
		       fCharm3ProngTClArr,
		       fCharm4ProngTClArr,
		       fDstarTClArr);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHF::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEVertexingHF: Terminate() \n");
}
