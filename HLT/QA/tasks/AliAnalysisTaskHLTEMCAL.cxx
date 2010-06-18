// $Id: AliAnalysisTaskHLTEMCAL.cxx 40285 2010-04-09 14:04:51Z kkanaki $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Zhongbao Yin <zbyin@mail.ccnu.edu.cn>,                *
//*                  Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
//*                  Svein Lindal <svein.lindal@gmail.com>                 *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliAnalysisTaskHLTEMCAL.cxx
    @author Zhongbao Yin, Kalliopi Kanaki, Svein Lindal
    @date
    @brief
*/

#include <iostream>

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TString.h"
#include "TObjArray.h"
#include "TFile.h"

#include "AliESDEvent.h"
#include "AliESDRun.h"
#include "AliESDInputHandler.h"
#include "AliESDCaloCluster.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskHLTEMCAL.h"


ClassImp(AliAnalysisTaskHLTEMCAL)


//===========================================================================================

AliAnalysisTaskHLTEMCAL::AliAnalysisTaskHLTEMCAL(const char *name) : AliAnalysisTaskHLTCalo(name)
{
  // Constructor
}

void AliAnalysisTaskHLTEMCAL::CreateSpecificStuff(TList * fOutputList){
  return;
}


void AliAnalysisTaskHLTEMCAL::DoSpecificStuff(AliESDEvent * evESD, AliESDEvent * evHLTESD) {
  return;
}

Int_t AliAnalysisTaskHLTEMCAL::GetClusters(AliESDEvent * event, TRefArray * clusters) {
  return event->GetEMCALClusters(clusters);
}

Bool_t AliAnalysisTaskHLTEMCAL::IsThisDetector(AliESDCaloCluster * cluster) {
  return cluster->IsEMCAL();
}
