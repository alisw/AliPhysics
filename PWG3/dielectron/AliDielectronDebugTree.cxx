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
//                Dielectron DebugTree                                  //
//                                                                       //
//                                                                       //
/*
register variables from the variable manager. The output will be written
to a tree

NOTE: Please use with extream care! Only for debugging and test purposes!!!

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTreeStream.h>
#include <TString.h>

#include <AliAnalysisManager.h>

#include "AliDielectronPair.h"

#include "AliDielectronDebugTree.h"

ClassImp(AliDielectronDebugTree)

AliDielectronDebugTree::AliDielectronDebugTree() :
  TNamed(),
  fFileName("jpsi_debug.root"),
  fNVars(0),
  fNVarsLeg(0),
  fStreamer(0x0)
{
  //
  // Default Constructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues;++i){
    fVariables[i]=0;
    fVariablesLeg[i]=0;
  }
}

//______________________________________________
AliDielectronDebugTree::AliDielectronDebugTree(const char* name, const char* title) :
  TNamed(name, title),
  fFileName("jpsi_debug.root"),
  fNVars(0),
  fNVarsLeg(0),
  fStreamer(0x0)
{
  //
  // Named Constructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues;++i){
    fVariables[i]=0;
    fVariablesLeg[i]=0;
  }
}

//______________________________________________
AliDielectronDebugTree::~AliDielectronDebugTree()
{
  //
  // Default Destructor
  //
  if (fStreamer){
    fStreamer->GetFile()->Write();
    delete fStreamer;
  }
}

//______________________________________________
void AliDielectronDebugTree::Fill(AliDielectronPair *pair)
{
  //
  // Fill configured variables to the tree
  //

  //is there anything to fill
  if (fNVars==0&&fNVarsLeg==0) return;
  
  //only in local mode!!!
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (man && man->GetAnalysisType()!=AliAnalysisManager::kLocalAnalysis) return;
  
  if (!fStreamer) fStreamer=new TTreeSRedirector(fFileName.Data());

  Int_t var=0;
  Double_t values[AliDielectronVarManager::kNMaxValues];
  Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues];
  Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues];
// fill pair values
  if (fNVars>0){
    AliDielectronVarManager::Fill(pair,values);

    for (Int_t i=0; i<fNVars; ++i){
      var=fVariables[i];
      (*fStreamer) << "Pair"
                   << Form("%s=",AliDielectronVarManager::GetValueName(var))
                   << values[var];
    }
  }

  if (fNVarsLeg>0){
    //leg1
    AliDielectronVarManager::Fill(pair->GetFirstDaughter(),valuesLeg1);
    //leg2
    AliDielectronVarManager::Fill(pair->GetSecondDaughter(),valuesLeg2);
    
    for (Int_t i=0; i<fNVarsLeg; ++i){
      var=fVariablesLeg[i];
      (*fStreamer) << "Pair"
                   << Form("Leg1_%s=",AliDielectronVarManager::GetValueName(var))
                   << valuesLeg1[var]
                   << Form("Leg2_%s=",AliDielectronVarManager::GetValueName(var))
                   << valuesLeg2[var];
    }
    
  }
  
  (*fStreamer) << "Pair" << "\n";
    
  
}

//______________________________________________
void AliDielectronDebugTree::DeleteStreamer()
{
  //
  // delete the streamer
  //
  if (!fStreamer) return;
  delete fStreamer;
  fStreamer=0x0;

}

//______________________________________________
void AliDielectronDebugTree::WriteTree()
{
  //
  // Write file
  //
  if (!fStreamer) return;
  fStreamer->GetFile()->Write();
}
