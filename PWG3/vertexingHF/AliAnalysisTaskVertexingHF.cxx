/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
//     AliAnalysisTask for the reconstruction of heavy flavor
// decays, using the class AliAnalysisVertexingHF.
//
// Author: J.Faivre, julien.faivre@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TSystem.h>
#include "Riostream.h"

#include "AliAnalysisTask.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskVertexingHF.h"

ClassImp(AliAnalysisTaskVertexingHF)


//________________________________________________________________________
AliAnalysisTaskVertexingHF::AliAnalysisTaskVertexingHF(const char *name) :
AliAnalysisTask(name,""), 
fESD(0), 
fChain(0), 
fVHF(0), 
fTrees(0)
{
  //Constructor

  //Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  //Output slots 0 to 3 write into a TTree
  DefineOutput(0, TTree::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}
//________________________________________________________________________
AliAnalysisTaskVertexingHF::~AliAnalysisTaskVertexingHF()
{
  // Destructor
  if (fTrees) delete [] fTrees;
}  
//________________________________________________________________________
void AliAnalysisTaskVertexingHF::ConnectInputData(Option_t *)
{
  //Implementation of AliAnalysisTask::ConnectInputData

  Info("ConnectInputData","ConnectInputData of task %s\n",GetName());
  fChain = (TChain*)GetInputData(0);

  if (!fChain) {
    printf("ERROR: Could not read chain from input slot 0\n");
    return;
  } 

  //fESD = new AliESDEvent();
  //fESD->ReadFromTree(fChain); 
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  
  if (!esdH) {
    printf("ERROR: Could not get ESDInputHandler\n");
  } else {
    fESD = esdH->GetEvent();
    if (!fESD) printf("ERROR: Could not get ESD object\n"); 
 }
 
  return;
}
//________________________________________________________________________
void AliAnalysisTaskVertexingHF::CreateOutputObjects()
{
  //Implementation of AliAnalysisTask::CreateOutputObjects
  //4 output trees (D0 in 2-prongs, J/Psi to e+e-, 3-prongs (D+, Ds, /\c), D0 in 4-prongs)

  fTrees = new TTree*[4];
  AliAODRecoDecayHF2Prong *rd2=0;
  AliAODRecoDecayHF3Prong *rd3=0;
  AliAODRecoDecayHF4Prong *rd4=0;
  fTrees[0] = new TTree("NameD0toKpi","TitleD0toKpi");
  fTrees[0]->Branch("D0toKpi","AliAODRecoDecayHF2Prong",&rd2);
  fTrees[1] = new TTree("NameJPSItoEle", "TitleJPSItoEle");
  fTrees[1]->Branch("JPSItoEle","AliAODRecoDecayHF2Prong",&rd2);
  fTrees[2] = new TTree("NameCharmto3Prong", "TitleCharmto3Prong");
  fTrees[2]->Branch("Charmto3Prong","AliAODRecoDecayHF3Prong",&rd3);
  fTrees[3] = new TTree("NameD0to4Prong","TitleD0to4Prong");
  fTrees[3]->Branch("D0to4Prong","AliAODRecoDecayHF4Prong",&rd4);

  return;
}
//________________________________________________________________________
void AliAnalysisTaskVertexingHF::LocalInit()
{
  //Instanciates vHF and loads its parameters

  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  fVHF->PrintStatus();

  return;
}
//________________________________________________________________________
void AliAnalysisTaskVertexingHF::Exec(Option_t *)
{
  //Performs heavy flavor vertexing
  
  //Heavy flavor vertexing :
  //fVHF->FindCandidates(fESD,fTrees);
  
  printf(" NOTHING DONE, THIS CLASS IS OBSOLETE, PLEASE USE AliAnalysisTaskSEVertexingHF\n");

  //Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fTrees[0]);
  PostData(1, fTrees[1]);
  PostData(2, fTrees[2]);
  PostData(3, fTrees[3]);

  return;
}      
//________________________________________________________________________
void AliAnalysisTaskVertexingHF::Terminate(Option_t *)
{
  //Implementation of AliAnalysisTask::Terminate

  return;
}









