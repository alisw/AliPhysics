/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


// $Id$

// AliAnalysisMuMuConfig : class to hold various configuration
// options for the AliAnalysisMuMu and AliAnalysisMuMuEvolution classes
// like the list of triggers to consider, the fit to be performed, etc...
// both for real data and for simulations (which might differ in e.g.
// the naming of the triggers)
//
// author: Laurent Aphecetche, Subatech
//
//
// TODO : make it readeable/writeable from/to a simple ASCII file ?
//

#include "AliAnalysisMuMuConfig.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "Riostream.h"
#include "TColor.h"
#include "TStyle.h"

ClassImp(AliAnalysisMuMuConfig)

//_____________________________________________________________________________
AliAnalysisMuMuConfig::AliAnalysisMuMuConfig(const char* beamYear) : TObject(),
fLists(new TObjArray),
fOCDBPath("raw://"),
fIsCompactGraphs(kFALSE)
{
  // ctor
  
  fLists->SetOwner(kTRUE);
  
  fLists->Add(new TObjArray); // list for real data
  fLists->Add(new TObjArray); // list for simulations

  DefineDefaults(beamYear);  
}

//_____________________________________________________________________________
AliAnalysisMuMuConfig::~AliAnalysisMuMuConfig()
{
  // dtor
  delete fLists;
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::DefineDefaults(const char* beamYear)
{
  // define some sensible defaults
  
  TString sbeam(beamYear);
  
  SetList(kDimuonTriggerList,kTRUE,"CMULLO-B-NOPF-MUON");
  SetList(kMuonTriggerList,kTRUE,"CMSNGL-B-NOPF-MUON");
  SetList(kEventSelectionList,kTRUE,"ALL");
  SetList(kPairSelectionList,kTRUE,"pRABSETAMATCHLOWPAIRY"); //pRABSETAPDCAMATCHLOWPAIRYPAIRPTIN0.0-15.0 //pRABSETAMATCHLOWPAIRYPAIRPTIN0.0-15.0
  SetList(kCentralitySelectionList,kTRUE,"V0A");
SetList(kFitTypeList,kTRUE,"func=PSICB2:histoType=minvJPsi:rebin=2:range=2.2;3.9,func=PSINA60NEW:histoType=minvJPsi:rebin=2:range=1.5;4.2,func=PSICB2:histoType=minvPsiP:rebin=2:range=2.0;4.2,func=PSICOUNT:histoType=minv");
//  SetList(kFitTypeList,kTRUE,"func=PSICB2:histoType=minvJPsi:rebin=2:range=2.2;3.9,func=PSICOUNT:histoType=minv");
  
  //,func=PSINA60NEW:histoType=minvJPsi:rebin=2:range=1.5;4.2
  //func=PSINA60NEW:histoType=minvJPsi:range=1.6;4.0:rebin=2,func=PSINA60NEW:histoType=minvPsiP:range=1.4;5.0:rebin=2
  //:alJPsi=1.0469:nlJPsi=4.1687:auJPsi=2.2517:nuJPsi=3.0778 (JPsi tails from raw spectra)
  //:alPsiP=1.0289:nlPsiP=3.86131:auPsiP=2.2737:nuPsiP=2.8995 (JPsi tails from AccxEff corr spectra to apply to PsiP spectra)
  
  if (sbeam=="pPb2013" || sbeam=="Pbp2013")
  {
    SetList(kDimuonTriggerList,kFALSE,"CMUL7-B-NOPF-MUON");
    SetList(kMuonTriggerList,kFALSE,"CMSL7-B-NOPF-MUON");
    SetList(kMinbiasTriggerList,kFALSE,"CINT7-B-NOPF-ALLNOTRD");
    SetList(kEventSelectionList,kFALSE,"PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00,PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50,PSALLHASSPD,PSALL");
    SetList(kPairSelectionList,kFALSE,"pRABSETAPDCAMATCHLOWPAIRYPAIRPTIN0.0-15.0");
    SetList(kCentralitySelectionList,kFALSE,"V0A");
 
//    SetList(kFitTypeList,kFALSE,"func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=0.9,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.1,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=0.9,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.1,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt");
    
SetList(kFitTypeList,kFALSE,"func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt");
    
//    SetList(kFitTypeList,kFALSE,"func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=0.9,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.1,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=0.9,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWVWG:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWVWG:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWVWG:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWVWG:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWVWG:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWVWG:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=0.9,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.1,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=0.9,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt");
    
    //,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt //Removed because bad results for pPb mean pt

//    SetList(kFitTypeList,kFALSE,"func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=0.9,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.1,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=0.9,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWVWG:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWVWG:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWVWG:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWVWG:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWVWG:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWVWG:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=0.9,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.1,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=0.9,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt");
    
    //,func=PSIPSIPRIMENA60NEWVWG:rebin=2:histoType=minv:tails=mctails
    //func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctailsJPsi&PsiP, func=PSIPSIPRIMECB2VWG;MPT2CB2VWGPOL2:rebin=2:histoType=minv&mpt:tails=mctailsJsi&PsiP(We can think about somth like this to make the combined fits),func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:alJPsi=0.984:nlJPsi=5.839:auJPsi=1.972:nuJPsi=3.444
    //func=PSIPSIPRIMECB2VWGINDEPTAILS:rebin=2:tails=mctails:histoType=minv,
    //Tails key migth be unneccesary since we have already different fitting function names(think about it)
  }
  else if (sbeam=="pp2012_7")
  {
    //    SetList(kDimuonTriggerList,kFALSE,"CMUL7-S-NOPF-MUON");
    SetList(kDimuonTriggerList,kFALSE,"CMUL7-S-NOPF-MUON");
    SetList(kMuonTriggerList,kFALSE,"CMSL7-S-NOPF-MUON");
    SetList(kMinbiasTriggerList,kFALSE,"CINT7-S-NOPF-ALLNOTRD");
    SetList(kEventSelectionList,kFALSE,"PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00");
    SetList(kPairSelectionList,kFALSE,"pRABSETAPDCAMATCHLOWPAIRYPAIRPTIN0.0-15.0");
    SetList(kCentralitySelectionList,kFALSE,"V0A");
    
    SetList(kFitTypeList,kFALSE,"func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt");//:sigmapsip=1.0
    SetList(kFitTypeList,kFALSE,"func=PSIPSIPRIMECB2POL4EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt");

    SetList(kFitTypeList,kFALSE,"func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=0.9,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.1,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=0.9,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.,func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWVWG:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWVWG:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWVWG:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWVWG:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWVWG:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWVWG:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=0.9,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.0;5.0:fsigmapsip=1.1,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=0.9,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.,func=PSIPSIPRIMECB2POL2EXP:rebin=2:histoType=minv:tails=mctails:range=2.2;4.7:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.0;5.0:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=0.9,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.,func=PSIPSIPRIMENA60NEWPOL2EXP:range=2.2;4.7:rebin=2:histoType=minv:tails=mctails:fsigmapsip=1.1,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2:range=2.2;4.7:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2EXP:range=2.0;5.:rebin=2:histoType=mpt,func=MPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2EXP:range=2.2;4.7:rebin=2:histoType=mpt");
    
//SetList(kFitTypeList,kFALSE,"func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:rebin=2:histoType=mpt");
    //MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2
    //func=PSIPSIPRIMECB2VWGINDEPTAILS:rebin=2:tails=mctails:histoType=minv,
    //func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctailsJPsi&PsiP, func=PSIPSIPRIMECB2VWG;MPT2CB2VWGPOL2:rebin=2:histoType=minv&mpt:tails=mctailsJsi&PsiP(We can think about somth like this to make the combined fits),func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:alJPsi=0.984:nlJPsi=5.839:auJPsi=1.972:nuJPsi=3.444
    //Tails key migth be unneccesary since we have already different fitting function names(think about it)
  }
  
  else if (sbeam=="pp2012_8")
  {
    //    SetList(kDimuonTriggerList,kFALSE,"CMUL7-S-NOPF-MUON");
    SetList(kDimuonTriggerList,kFALSE,"CMUL8-S-NOPF-MUON");
    SetList(kMuonTriggerList,kFALSE,"CMSL7-8-NOPF-MUON");
    SetList(kMinbiasTriggerList,kFALSE,"CINT8-S-NOPF-ALLNOTRD");
    SetList(kEventSelectionList,kFALSE,"PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00");
    SetList(kPairSelectionList,kFALSE,"pRABSETAPDCAMATCHLOWPAIRYPAIRPTIN0.0-15.0");
    SetList(kCentralitySelectionList,kFALSE,"V0A");
    SetList(kFitTypeList,kFALSE,"func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctails,func=MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2:rebin=2:histoType=mpt");
    //MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2
    //func=PSIPSIPRIMECB2VWGINDEPTAILS:rebin=2:tails=mctails:histoType=minv,
    //func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:tails=mctailsJPsi&PsiP, func=PSIPSIPRIMECB2VWG;MPT2CB2VWGPOL2:rebin=2:histoType=minv&mpt:tails=mctailsJsi&PsiP(We can think about somth like this to make the combined fits),func=PSIPSIPRIMECB2VWG:rebin=2:histoType=minv:alJPsi=0.984:nlJPsi=5.839:auJPsi=1.972:nuJPsi=3.444
    //Tails key migth be unneccesary since we have already different fitting function names(think about it)
  }

}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::GetList(ETypeList type, Bool_t simulation) const
{
  /// Get the value for a given list. If simulation=true and the list is not
  /// there for that type, the list from real data is returned
  TObjArray* array = static_cast<TObjArray*>(fLists->At(simulation));
  TObjString* str = static_cast<TObjString*>(array->At(type));
  if ( !str && simulation )
  {
    return GetList(type,kFALSE);
  }
  return str->String();
}

//_____________________________________________________________________________
TObjArray* AliAnalysisMuMuConfig::GetListElements(ETypeList type, Bool_t simulation) const
{
  /// Get list as an array (to be deleted by the user)
  TString list = GetList(type,simulation);
  return list.Tokenize(",");
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::SetList(ETypeList type, Bool_t simulation, const char* list)
{
  /// Set the list of a given type
  TObjArray* array = static_cast<TObjArray*>(fLists->At(simulation));
  TObjString* str = static_cast<TObjString*>(array->At(type));
  if (!str)
  {
    str = new TObjString;
    array->AddAt(str,type);
  }
  str->String() = list;
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::ShowLists(const char* title, ETypeList type, const char separator, const TString& opt) const
{
  /// Show the real and sim list of a given type
  
  std::cout << title << std::endl;
  
  TString list;
  
  if ( opt.Contains("REAL",TString::kIgnoreCase) )
  {
    list = GetList(type,kFALSE);
    ShowList("real",list,separator);
  }
  if ( opt.Contains("SIM",TString::kIgnoreCase) )
  {
    list = GetList(type,kTRUE);
    ShowList("sim",list,separator);
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::ShowList(const char* title, const TString& list, const char separator) const
{
  /// Show the list content
  
  TObjArray* parts = list.Tokenize(separator);

  TIter next(parts);
  TObjString* str;
  
  std::cout << "   " << title << " (" << parts->GetEntries() << ")" << std::endl;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    std::cout << "       " << str->String().Data() << std::endl;
  }
  
  delete parts;
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::Print(Option_t* opt) const
{
  /// printout
  /// Use opt = "REAL" to show only things relevant to real data
  /// Use opt = "SIM" to show only things relevant to simulation
  /// Use opt = "REAL SIM" or "" to show everything
  
  TString sopt(opt);
  sopt.ToUpper();
  
  if (sopt.Length()==0)
  {
    sopt = "REAL SIM";
  }
  
  ShowLists("Dimuon triggers",kDimuonTriggerList,',',sopt.Data());
  ShowLists("Muon triggers",kMuonTriggerList,',',sopt.Data());
  ShowLists("MB triggers",kMinbiasTriggerList,',',sopt.Data());
  ShowLists("Event selection",kEventSelectionList,',',sopt.Data());
  ShowLists("Pair selection",kPairSelectionList,',',sopt.Data());
  ShowLists("Centrality selection",kCentralitySelectionList,',',sopt.Data());
  ShowLists("Fit types",kFitTypeList,',',sopt.Data());
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::SetColorScheme()
{
  /// Set a few custom colors
  
  new TColor(AliAnalysisMuMuConfig::kBlue,4/255.0,44/255.0,87/255.0,"my blue");
  new TColor(AliAnalysisMuMuConfig::kOrange,255/255.0,83/255.0,8/255.0,"my orange");
  new TColor(AliAnalysisMuMuConfig::kGreen,152/255.0,202/255.0,52/255.0,"my green");
  
  gStyle->SetGridColor(AliAnalysisMuMuConfig::kBlue);
  
  gStyle->SetFrameLineColor(AliAnalysisMuMuConfig::kBlue);
  gStyle->SetAxisColor(AliAnalysisMuMuConfig::kBlue,"xyz");
  gStyle->SetLabelColor(AliAnalysisMuMuConfig::kBlue,"xyz");
  
  gStyle->SetTitleColor(AliAnalysisMuMuConfig::kBlue);
  gStyle->SetTitleTextColor(AliAnalysisMuMuConfig::kBlue);
  gStyle->SetLabelColor(AliAnalysisMuMuConfig::kBlue);
  gStyle->SetStatTextColor(AliAnalysisMuMuConfig::kBlue);
  
  gStyle->SetOptStat(0);
}
