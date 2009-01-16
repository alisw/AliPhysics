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

/*
$Log$
*/ 

#include "AliFlowCommonConstants.h" 
#include "TMath.h" 

// AliFlowCommonConstants:
//
// Constants for the common histograms in the flow analysis
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)

//ClassImp(AliFlowCommonConstants)

Double_t AliFlowCommonConstants::fgMultMin =  0.;            
Double_t AliFlowCommonConstants::fgMultMax = 10000.;
Double_t AliFlowCommonConstants::fgPtMin   =  0.;	     
Double_t AliFlowCommonConstants::fgPtMax   = 10.;
Double_t AliFlowCommonConstants::fgPhiMin  =  0.;	     
Double_t AliFlowCommonConstants::fgPhiMax  =  TMath::TwoPi();
Double_t AliFlowCommonConstants::fgEtaMin  = -2.;	     
Double_t AliFlowCommonConstants::fgEtaMax  =  2.;	     
Double_t AliFlowCommonConstants::fgQMin    =  0.;	     
Double_t AliFlowCommonConstants::fgQMax    =  3.;

//getters
Int_t AliFlowCommonConstants::GetNbinsMult() { return kNbinsMult; }
Int_t AliFlowCommonConstants::GetNbinsPt()   { return kNbinsPt; }
Int_t AliFlowCommonConstants::GetNbinsPhi()  { return kNbinsPhi; }
Int_t AliFlowCommonConstants::GetNbinsEta()  { return kNbinsEta; }
Int_t AliFlowCommonConstants::GetNbinsQ()    { return kNbinsQ; }
 
//getters
Double_t AliFlowCommonConstants::GetMultMin() { return fgMultMin; }
Double_t AliFlowCommonConstants::GetMultMax() { return fgMultMax; }
Double_t AliFlowCommonConstants::GetPtMin()   { return fgPtMin; }
Double_t AliFlowCommonConstants::GetPtMax()   { return fgPtMax; }
Double_t AliFlowCommonConstants::GetPhiMin()  { return fgPhiMin; }
Double_t AliFlowCommonConstants::GetPhiMax()  { return fgPhiMax; }
Double_t AliFlowCommonConstants::GetEtaMin()  { return fgEtaMin; }
Double_t AliFlowCommonConstants::GetEtaMax()  { return fgEtaMax; }
Double_t AliFlowCommonConstants::GetQMin()    { return fgQMin; }
Double_t AliFlowCommonConstants::GetQMax()    { return fgQMax; }
  
