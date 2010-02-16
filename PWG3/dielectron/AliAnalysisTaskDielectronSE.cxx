/*************************************************************************
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                        Basic Analysis Task                            //
//                      for Dielectron Analysis                          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>

#include <AliCFContainer.h>
#include <AliVEvent.h>

#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronCF.h"
#include "AliAnalysisTaskDielectronSE.h"

ClassImp(AliAnalysisTaskDielectronSE)

//_________________________________________________________________________________
AliAnalysisTaskDielectronSE::AliAnalysisTaskDielectronSE() :
  AliAnalysisTaskSE(),
  fDielectron(0)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskDielectronSE::AliAnalysisTaskDielectronSE(const char *name) :
  AliAnalysisTaskSE(name),
  fDielectron(0)
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1, THashList::Class());
  DefineOutput(2, AliCFContainer::Class());
}

//_________________________________________________________________________________
void AliAnalysisTaskDielectronSE::UserCreateOutputObjects()
{
  //
  // Initialise the framework objects
  //
  if (!fDielectron){
    AliError("No Dielectron framework object set !!!");
    return;
  }
  fDielectron->Init();
}

//_________________________________________________________________________________
void AliAnalysisTaskDielectronSE::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //

  if (!fDielectron) return;
  
  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );
  
  fDielectron->Process(InputEvent());
  fDielectron->FillHistograms();

  if (fDielectron->GetHistogramList()){
    PostData(1, const_cast<THashList*>(fDielectron->GetHistogramList()));
  }
  if (fDielectron->GetCFManagerPair()){
    PostData(2, const_cast<AliCFContainer*>(fDielectron->GetCFManagerPair()->GetContainer()));
  }
}

