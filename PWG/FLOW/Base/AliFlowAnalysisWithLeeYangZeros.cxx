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

/////////////////////////////////////////////////////////////////////////////////////////
//Description: Maker to analyze Flow using the LeeYangZeros method  
//             Equation numbers are from Big Paper (BP): Nucl. Phys. A 727, 373 (2003)
//             Practical Guide (PG):    J. Phys. G: Nucl. Part. Phys. 30, S1213 (2004)  
//             Adapted from StFlowLeeYangZerosMaker.cxx           
//             by Markus Oldenberg and Art Poskanzer, LBNL        
//             with advice from Jean-Yves Ollitrault and Nicolas Borghini   
//
//Author: Naomi van der Kolk (kolk@nikhef.nl)
/////////////////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TObject.h" 
#include "TMath.h" 
#include "TFile.h"
#include "TList.h"
#include "TVector2.h"
#include "TH1F.h"
#include "TComplex.h"
#include "TProfile.h"
#include "TDirectoryFile.h"

#include "AliFlowCommonConstants.h"
#include "AliFlowLYZConstants.h"
#include "AliFlowAnalysisWithLeeYangZeros.h"
#include "AliFlowLYZHist1.h"
#include "AliFlowLYZHist2.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowVector.h"

ClassImp(AliFlowAnalysisWithLeeYangZeros)

  //-----------------------------------------------------------------------
 
  AliFlowAnalysisWithLeeYangZeros::AliFlowAnalysisWithLeeYangZeros():
    fQsum(NULL),
    fQ2sum(0),
    fEventNumber(0),
    fFirstRun(kTRUE),
    fUseSum(kTRUE),
    fDoubleLoop(kFALSE),
    fDebug(kFALSE),
    fHistList(NULL),
    fFirstRunList(NULL),
    fHistVtheta(NULL),
    fHistProVetaRP(NULL),  
    fHistProVetaPOI(NULL), 
    fHistProVPtRP(NULL),   
    fHistProVPtPOI(NULL),  
    fHistR0theta(NULL),
    fHistProReDenom(NULL),
    fHistProImDenom(NULL),
    fHistReDtheta(NULL),
    fHistImDtheta(NULL),
    fHistQsumforChi(NULL),
    fCommonHists(NULL),
    fCommonHistsRes(NULL)
  
{
  //default constructor
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::AliFlowAnalysisWithLeeYangZeros default constructor****"<<endl;

  fHistList = new TList();
  fFirstRunList = new TList();

  for(Int_t i = 0;i<5;i++)
    {
      fHist1[i]=0;
      fHist2RP[i]=0;
      fHist2POI[i]=0;
    }

  fQsum = new TVector2();

}

//-----------------------------------------------------------------------


 AliFlowAnalysisWithLeeYangZeros::~AliFlowAnalysisWithLeeYangZeros() 
 {
   //default destructor
   if (fDebug) cout<<"****~AliFlowAnalysisWithLeeYangZeros****"<<endl;
   delete fQsum;
   delete fHistList;
   delete fFirstRunList;
   
 }
 
//-----------------------------------------------------------------------

void AliFlowAnalysisWithLeeYangZeros::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file

  TFile *output = new TFile(outputFileName->Data(),"RECREATE");
  if (GetFirstRun()) {
    //output->WriteObject(fHistList, "cobjLYZ1","SingleKey");
    if (fUseSum) { fHistList->SetName("cobjLYZ1SUM");}
    else {fHistList->SetName("cobjLYZ1PROD");}
    fHistList->SetOwner(kTRUE);
    fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
  }
  else {
    //output->WriteObject(fHistList, "cobjLYZ2","SingleKey");
    if (fUseSum) { fHistList->SetName("cobjLYZ2SUM"); }
    else { fHistList->SetName("cobjLYZ2PROD"); }
    fHistList->SetOwner(kTRUE);
    fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
  }
  delete output;
}

//-----------------------------------------------------------------------

void AliFlowAnalysisWithLeeYangZeros::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file

  TFile *output = new TFile(outputFileName.Data(),"RECREATE");
  if (GetFirstRun()) {
    //output->WriteObject(fHistList, "cobjLYZ1","SingleKey");
    if (fUseSum) { fHistList->SetName("cobjLYZ1SUM");}
    else {fHistList->SetName("cobjLYZ1PROD");}
    fHistList->SetOwner(kTRUE);
    fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
  }
  else {
    //output->WriteObject(fHistList, "cobjLYZ2","SingleKey");
    if (fUseSum) { fHistList->SetName("cobjLYZ2SUM"); }
    else { fHistList->SetName("cobjLYZ2PROD"); }
    fHistList->SetOwner(kTRUE);
    fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
  }
  delete output;
}

//-----------------------------------------------------------------------

void AliFlowAnalysisWithLeeYangZeros::WriteHistograms(TDirectoryFile *outputFileName)
{
 //store the final results in output .root file
 if (GetFirstRun()) {
   if (fUseSum) { fHistList->SetName("cobjLYZ1SUM");}
   else {fHistList->SetName("cobjLYZ1PROD");}
   fHistList->SetOwner(kTRUE);
   outputFileName->Add(fHistList);
   outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey); 
 }
 else {
  if (fUseSum) { fHistList->SetName("cobjLYZ2SUM"); }
  else { fHistList->SetName("cobjLYZ2PROD"); }
  fHistList->SetOwner(kTRUE);
  outputFileName->Add(fHistList);
  outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey); 
 }
}

//-----------------------------------------------------------------------

Bool_t AliFlowAnalysisWithLeeYangZeros::Init() 
{
  //init method 
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::Init()****"<<endl;

  //save old value and prevent histograms from being added to directory
  //to avoid name clashes in case multiple analaysis objects are used
  //in an analysis
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
 
  // Book histograms
  Int_t iNtheta = AliFlowLYZConstants::GetMaster()->GetNtheta();
  Int_t iNbinsPt = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
  Int_t iNbinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();

  Double_t  dPtMin = AliFlowCommonConstants::GetMaster()->GetPtMin();	     
  Double_t  dPtMax = AliFlowCommonConstants::GetMaster()->GetPtMax();
  Double_t  dEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();	     
  Double_t  dEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();
    
  //for control histograms
  if (fFirstRun){ 
    if (fUseSum) {fCommonHists = new AliFlowCommonHist("AliFlowCommonHistLYZ1SUM");}
    else {fCommonHists = new AliFlowCommonHist("AliFlowCommonHistLYZ1PROD");}
    fHistList->Add(fCommonHists);
    if (fUseSum) {fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsLYZ1SUM");}
    else {fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsLYZ1PROD");}
    fHistList->Add(fCommonHistsRes);
  }
  else { 
    if (fUseSum) {fCommonHists = new AliFlowCommonHist("AliFlowCommonHistLYZ2SUM");}
    else {fCommonHists = new AliFlowCommonHist("AliFlowCommonHistLYZ2PROD");}
    fHistList->Add(fCommonHists);
    if (fUseSum) {fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsLYZ2SUM");}
    else {fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsLYZ2PROD");}
    fHistList->Add(fCommonHistsRes); 
  }
  
  TString nameChiHist;
  if (fUseSum) {nameChiHist = "Flow_QsumforChi_LYZSUM";}
  else {nameChiHist = "Flow_QsumforChi_LYZPROD";}
  fHistQsumforChi = new TH1F(nameChiHist.Data(),nameChiHist.Data(),3,-1.,2.);
  fHistQsumforChi->SetXTitle("Qsum.X , Qsum.Y, Q2sum");
  fHistQsumforChi->SetYTitle("value");
  fHistList->Add(fHistQsumforChi);

  //for first loop over events 
  if (fFirstRun){
    TString nameR0Hist;
    if (fUseSum) {nameR0Hist = "First_Flow_r0theta_LYZSUM";}
    else {nameR0Hist = "First_Flow_r0theta_LYZPROD";}
    fHistR0theta  = new TH1D(nameR0Hist.Data(),nameR0Hist.Data(),iNtheta,-0.5,iNtheta-0.5);
    fHistR0theta->SetXTitle("#theta");
    fHistR0theta->SetYTitle("r_{0}^{#theta}");
    fHistList->Add(fHistR0theta);

    TString nameVHist;
    if (fUseSum) {nameVHist = "First_Flow_Vtheta_LYZSUM";}
    else {nameVHist = "First_Flow_Vtheta_LYZPROD";}
    fHistVtheta  = new TH1D(nameVHist.Data(),nameVHist.Data(),iNtheta,-0.5,iNtheta-0.5);
    fHistVtheta->SetXTitle("#theta");
    fHistVtheta->SetYTitle("V_{n}^{#theta}");        
    fHistList->Add(fHistVtheta);

    //class AliFlowLYZHist1 defines the histograms: fHistProGtheta, fHistProReGtheta, fHistProImGtheta
    for (Int_t theta=0;theta<iNtheta;theta++) {  
      TString nameHist1;
      if (fUseSum) {nameHist1 = "AliFlowLYZHist1_";}
      else {nameHist1 = "AliFlowLYZHist1_";}
      nameHist1 += theta;
      fHist1[theta]=new AliFlowLYZHist1(theta, nameHist1, fUseSum);
      fHistList->Add(fHist1[theta]);
    }
         
  }
  //for second loop over events 
  else {
    TString nameReDenomHist;
    if (fUseSum) {nameReDenomHist = "Second_FlowPro_ReDenom_LYZSUM";}
    else {nameReDenomHist = "Second_FlowPro_ReDenom_LYZPROD";}
    fHistProReDenom = new TProfile(nameReDenomHist.Data(),nameReDenomHist.Data(), iNtheta, -0.5, iNtheta-0.5);
    fHistProReDenom->SetXTitle("#theta");
    fHistProReDenom->SetYTitle("Re(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");
    fHistList->Add(fHistProReDenom);

    TString nameImDenomHist;
    if (fUseSum) {nameImDenomHist = "Second_FlowPro_ImDenom_LYZSUM";}
    else {nameImDenomHist = "Second_FlowPro_ImDenom_LYZPROD";}
    fHistProImDenom = new TProfile(nameImDenomHist.Data(),nameImDenomHist.Data(), iNtheta, -0.5, iNtheta-0.5);
    fHistProImDenom->SetXTitle("#theta");
    fHistProImDenom->SetYTitle("Im(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");
    fHistList->Add(fHistProImDenom);

    TString nameVetaRPHist;
    if (fUseSum) {nameVetaRPHist = "Second_FlowPro_VetaRP_LYZSUM";}
    else {nameVetaRPHist = "Second_FlowPro_VetaRP_LYZPROD";}
    fHistProVetaRP = new TProfile(nameVetaRPHist.Data(),nameVetaRPHist.Data(),iNbinsEta,dEtaMin,dEtaMax);
    fHistProVetaRP->SetXTitle("rapidity");
    fHistProVetaRP->SetYTitle("v_{2}(#eta) for RP selection");
    fHistList->Add(fHistProVetaRP);

    TString nameVetaPOIHist;
    if (fUseSum) {nameVetaPOIHist = "Second_FlowPro_VetaPOI_LYZSUM";}
    else {nameVetaPOIHist = "Second_FlowPro_VetaPOI_LYZPROD";}
    fHistProVetaPOI = new TProfile(nameVetaPOIHist.Data(),nameVetaPOIHist.Data(),iNbinsEta,dEtaMin,dEtaMax);
    fHistProVetaPOI->SetXTitle("rapidity");
    fHistProVetaPOI->SetYTitle("v_{2}(#eta) for POI selection");
    fHistList->Add(fHistProVetaPOI);

    TString nameVPtRPHist;
    if (fUseSum) {nameVPtRPHist = "Second_FlowPro_VPtRP_LYZSUM";}
    else {nameVPtRPHist = "Second_FlowPro_VPtRP_LYZPROD";}
    fHistProVPtRP = new TProfile(nameVPtRPHist.Data(),nameVPtRPHist.Data(),iNbinsPt,dPtMin,dPtMax);
    fHistProVPtRP->SetXTitle("Pt");
    fHistProVPtRP->SetYTitle("v_{2}(p_{T}) for RP selection");
    fHistList->Add(fHistProVPtRP);

    TString nameVPtPOIHist;
    if (fUseSum) {nameVPtPOIHist = "Second_FlowPro_VPtPOI_LYZSUM";}
    else {nameVPtPOIHist = "Second_FlowPro_VPtPOI_LYZPROD";}
    fHistProVPtPOI = new TProfile(nameVPtPOIHist.Data(),nameVPtPOIHist.Data(),iNbinsPt,dPtMin,dPtMax);
    fHistProVPtPOI->SetXTitle("p_{T}");
    fHistProVPtPOI->SetYTitle("v_{2}(p_{T}) for POI selection");
    fHistList->Add(fHistProVPtPOI);

    if (fUseSum){
      fHistReDtheta = new TH1D("Second_Flow_ReDtheta_LYZSUM","Second_Flow_ReDtheta_LYZSUM",iNtheta, -0.5, (double)iNtheta-0.5);
      fHistReDtheta->SetXTitle("#theta");
      fHistReDtheta->SetYTitle("Re(D^{#theta})");
      fHistList->Add(fHistReDtheta);

      fHistImDtheta = new TH1D("Second_Flow_ImDtheta_LYZSUM","Second_Flow_ImDtheta_LYZSUM",iNtheta, -0.5, (double)iNtheta-0.5);
      fHistImDtheta->SetXTitle("#theta");
      fHistImDtheta->SetYTitle("Im(D^{#theta})");
      fHistList->Add(fHistImDtheta);
    }

    //class AliFlowLYZHist2 defines the histograms: 
    for (Int_t theta=0;theta<iNtheta;theta++)  {  
      TString nameRP = "AliFlowLYZHist2RP_";
      nameRP += theta;
      fHist2RP[theta]=new AliFlowLYZHist2(theta, "RP", nameRP, fUseSum);
      fHistList->Add(fHist2RP[theta]);

      TString namePOI = "AliFlowLYZHist2POI_";
      namePOI += theta;
      fHist2POI[theta]=new AliFlowLYZHist2(theta, "POI", namePOI, fUseSum);
      fHistList->Add(fHist2POI[theta]);
    }
     
    //read histogram fHistR0theta from the first run list
    if (fFirstRunList) {
      if (fUseSum) { fHistR0theta  = (TH1D*)fFirstRunList->FindObject("First_Flow_r0theta_LYZSUM");}
      else{ fHistR0theta  = (TH1D*)fFirstRunList->FindObject("First_Flow_r0theta_LYZPROD");}
      if (!fHistR0theta) {cout<<"fHistR0theta has a NULL pointer!"<<endl;}
      fHistList->Add(fHistR0theta);
    } else { cout<<"list is NULL pointer!"<<endl; }

  }
   

  if (fDebug) cout<<"****Histograms initialised****"<<endl;
    
  fEventNumber = 0; //set event counter to zero
  
  //resore old stuff
  TH1::AddDirectory(oldHistAddStatus);

  return kTRUE; 
}
 
 //-----------------------------------------------------------------------
 
Bool_t AliFlowAnalysisWithLeeYangZeros::Make(AliFlowEventSimple* anEvent) 
{
  //make method
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::Make()****"<<endl;
        
  //get tracks from event
  if (anEvent) {
    if (fFirstRun){
      fCommonHists->FillControlHistograms(anEvent);
      FillFromFlowEvent(anEvent);
    }
    else {
      fCommonHists->FillControlHistograms(anEvent);
      SecondFillFromFlowEvent(anEvent);
    }
  }
 
  else {
    cout<<"##### FlowLeeYangZero: Stack pointer null"<<endl;
    return kFALSE;
  }
  //  cout<<"^^^^read event "<<fEventNumber<<endl;
  fEventNumber++;
  
     
  return kTRUE; 
}

   //-----------------------------------------------------------------------     
void AliFlowAnalysisWithLeeYangZeros::GetOutputHistograms(TList *outputListHistos) {
 // get the pointers to all output histograms before calling Finish()
 
  const Int_t iNtheta = AliFlowLYZConstants::GetMaster()->GetNtheta();
  
  if (outputListHistos) {

    //define histograms for first and second run
    AliFlowCommonHist *pCommonHist = NULL;
    AliFlowCommonHistResults *pCommonHistResults = NULL;
    TH1D*     pHistR0theta = NULL;
    TH1D*     pHistVtheta  = NULL;
    TProfile* pHistProReDenom = NULL;
    TProfile* pHistProImDenom = NULL;
    TProfile* pHistProVetaRP  = NULL;
    TProfile* pHistProVetaPOI = NULL;
    TProfile* pHistProVPtRP   = NULL;
    TProfile* pHistProVPtPOI  = NULL;
    TH1F* pHistQsumforChi = NULL;
    //AliFlowLYZHist1 *pLYZHist1[iNtheta] = {NULL};      //array of pointers to AliFlowLYZHist1
    //AliFlowLYZHist2 *pLYZHist2RP[iNtheta] = {NULL};    //array of pointers to AliFlowLYZHist2
    //AliFlowLYZHist2 *pLYZHist2POI[iNtheta] = {NULL};   //array of pointers to AliFlowLYZHist2
    AliFlowLYZHist1 **pLYZHist1 = new AliFlowLYZHist1*[iNtheta];      //array of pointers to AliFlowLYZHist1
    AliFlowLYZHist2 **pLYZHist2RP = new AliFlowLYZHist2*[iNtheta];    //array of pointers to AliFlowLYZHist2
    AliFlowLYZHist2 **pLYZHist2POI = new AliFlowLYZHist2*[iNtheta];   //array of pointers to AliFlowLYZHist2
    for (Int_t i=0; i<iNtheta; i++)
    {
      pLYZHist1[i] = NULL;
      pLYZHist2RP[i] = NULL;
      pLYZHist2POI[i] = NULL;
    }

    if (GetFirstRun()) { //first run
      //Get the common histograms from the output list
      if (GetUseSum()){
	pCommonHist = dynamic_cast<AliFlowCommonHist*> 
	  (outputListHistos->FindObject("AliFlowCommonHistLYZ1SUM")); 
	pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
	  (outputListHistos->FindObject("AliFlowCommonHistResultsLYZ1SUM"));
	//Get the histograms from the output list
	for(Int_t theta = 0;theta<iNtheta;theta++){
	  TString name = "AliFlowLYZHist1_"; 
	  name += theta;
	  pLYZHist1[theta] = dynamic_cast<AliFlowLYZHist1*> 
	    (outputListHistos->FindObject(name));
	}
	pHistVtheta = dynamic_cast<TH1D*> 
	  (outputListHistos->FindObject("First_Flow_Vtheta_LYZSUM"));

	pHistR0theta = dynamic_cast<TH1D*> 
	  (outputListHistos->FindObject("First_Flow_r0theta_LYZSUM"));

	pHistQsumforChi = dynamic_cast<TH1F*> 
	  (outputListHistos->FindObject("Flow_QsumforChi_LYZSUM"));

	//Set the histogram pointers and call Finish()
	if (pCommonHist && pCommonHistResults && pLYZHist1[0] && 
	    pHistVtheta && pHistR0theta && pHistQsumforChi ) {
	  this->SetCommonHists(pCommonHist);
	  this->SetCommonHistsRes(pCommonHistResults);
	  this->SetHist1(pLYZHist1);
	  this->SetHistVtheta(pHistVtheta);
	  this->SetHistR0theta(pHistR0theta);
	  this->SetHistQsumforChi(pHistQsumforChi);
	} 
	else { 
	  cout<<"WARNING: Histograms needed to run Finish() firstrun (SUM) are not accessable!"<<endl; 
	}
      }
      else {
	pCommonHist = dynamic_cast<AliFlowCommonHist*> 
	  (outputListHistos->FindObject("AliFlowCommonHistLYZ1PROD"));
	pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
	  (outputListHistos->FindObject("AliFlowCommonHistResultsLYZ1PROD"));
	//Get the histograms from the output list
	for(Int_t theta = 0;theta<iNtheta;theta++){
	  TString name = "AliFlowLYZHist1_"; 
	  name += theta;
	  pLYZHist1[theta] = dynamic_cast<AliFlowLYZHist1*> 
	    (outputListHistos->FindObject(name));
	}
	pHistVtheta = dynamic_cast<TH1D*> 
	  (outputListHistos->FindObject("First_Flow_Vtheta_LYZPROD"));

	pHistR0theta = dynamic_cast<TH1D*> 
	  (outputListHistos->FindObject("First_Flow_r0theta_LYZPROD"));

	pHistQsumforChi = dynamic_cast<TH1F*> 
	  (outputListHistos->FindObject("Flow_QsumforChi_LYZPROD"));

	//Set the histogram pointers and call Finish()
	if (pCommonHist && pCommonHistResults && pLYZHist1[0] && 
	    pHistVtheta && pHistR0theta && pHistQsumforChi ) {
	  this->SetCommonHists(pCommonHist);
	  this->SetCommonHistsRes(pCommonHistResults);
	  this->SetHist1(pLYZHist1);
	  this->SetHistVtheta(pHistVtheta);
	  this->SetHistR0theta(pHistR0theta);
	  this->SetHistQsumforChi(pHistQsumforChi);
	} else { 
	  cout<<"WARNING: Histograms needed to run Finish() firstrun (PROD) are not accessable!"<<endl; 
	}
      }
    }
    else { //second run
      //Get the common histograms from the output list
      if (GetUseSum()){
	pCommonHist = dynamic_cast<AliFlowCommonHist*> 
	  (outputListHistos->FindObject("AliFlowCommonHistLYZ2SUM"));
	pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
	  (outputListHistos->FindObject("AliFlowCommonHistResultsLYZ2SUM"));

	pHistR0theta = dynamic_cast<TH1D*> 
	  (outputListHistos->FindObject("First_Flow_r0theta_LYZSUM"));

	pHistQsumforChi = dynamic_cast<TH1F*> 
	  (outputListHistos->FindObject("Flow_QsumforChi_LYZSUM"));

	//Get the histograms from the output list
	for(Int_t theta = 0;theta<iNtheta;theta++){
	  TString nameRP = "AliFlowLYZHist2RP_"; 
	  nameRP += theta;
	  pLYZHist2RP[theta] = dynamic_cast<AliFlowLYZHist2*> 
	    (outputListHistos->FindObject(nameRP));
	  TString namePOI = "AliFlowLYZHist2POI_"; 
	  namePOI += theta;
	  pLYZHist2POI[theta] = dynamic_cast<AliFlowLYZHist2*> 
	    (outputListHistos->FindObject(namePOI));
	}
	pHistProReDenom = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_ReDenom_LYZSUM"));
	pHistProImDenom = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_ImDenom_LYZSUM"));
	
	TH1D* pHistReDtheta = dynamic_cast<TH1D*> 
	  (outputListHistos->FindObject("Second_Flow_ReDtheta_LYZSUM"));
	TH1D* pHistImDtheta = dynamic_cast<TH1D*> 
	  (outputListHistos->FindObject("Second_Flow_ImDtheta_LYZSUM"));
	
	pHistProVetaRP = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_VetaRP_LYZSUM"));
	pHistProVetaPOI = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_VetaPOI_LYZSUM"));
	pHistProVPtRP = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_VPtRP_LYZSUM"));
	pHistProVPtPOI = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_VPtPOI_LYZSUM"));


	//Set the histogram pointers and call Finish()
	if (pCommonHist && pCommonHistResults && 
	    pLYZHist2RP[0] && pLYZHist2POI[0] && 
	    pHistR0theta && 
	    pHistProReDenom && pHistProImDenom && 
	    pHistReDtheta && pHistImDtheta && 
	    pHistProVetaRP && pHistProVetaPOI && 
	    pHistProVPtRP && pHistProVPtPOI && 
	    pHistQsumforChi) {
	  this->SetCommonHists(pCommonHist);
	  this->SetCommonHistsRes(pCommonHistResults);
	  this->SetHist2RP(pLYZHist2RP);
	  this->SetHist2POI(pLYZHist2POI);
	  this->SetHistR0theta(pHistR0theta);
	  this->SetHistProReDenom(pHistProReDenom);
	  this->SetHistProImDenom(pHistProImDenom);
	  this->SetHistReDtheta(pHistReDtheta);
	  this->SetHistImDtheta(pHistImDtheta);
	  this->SetHistProVetaRP(pHistProVetaRP);
	  this->SetHistProVetaPOI(pHistProVetaPOI);
	  this->SetHistProVPtRP(pHistProVPtRP);
	  this->SetHistProVPtPOI(pHistProVPtPOI);
	  this->SetHistQsumforChi(pHistQsumforChi);
	} 
	else { 
	  cout<<"WARNING: Histograms needed to run Finish() secondrun (SUM) are not accessable!"<<endl; 
	}
      }
      else {
	pCommonHist = dynamic_cast<AliFlowCommonHist*> 
	  (outputListHistos->FindObject("AliFlowCommonHistLYZ2PROD"));
	pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
	  (outputListHistos->FindObject("AliFlowCommonHistResultsLYZ2PROD"));
      

	pHistR0theta = dynamic_cast<TH1D*> 
	  (outputListHistos->FindObject("First_Flow_r0theta_LYZPROD"));

	pHistQsumforChi = dynamic_cast<TH1F*> 
	  (outputListHistos->FindObject("Flow_QsumforChi_LYZPROD"));

	//Get the histograms from the output list
	for(Int_t theta = 0;theta<iNtheta;theta++){
	  TString nameRP = "AliFlowLYZHist2RP_"; 
	  nameRP += theta;
	  pLYZHist2RP[theta] = dynamic_cast<AliFlowLYZHist2*> 
	    (outputListHistos->FindObject(nameRP));
	  TString namePOI = "AliFlowLYZHist2POI_"; 
	  namePOI += theta;
	  pLYZHist2POI[theta] = dynamic_cast<AliFlowLYZHist2*> 
	    (outputListHistos->FindObject(namePOI));
	}
	pHistProReDenom = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_ReDenom_LYZPROD"));
	pHistProImDenom = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_ImDenom_LYZPROD"));
	
	pHistProVetaRP = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_VetaRP_LYZPROD"));
	pHistProVetaPOI = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_VetaPOI_LYZPROD"));
	pHistProVPtRP = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_VPtRP_LYZPROD"));
	pHistProVPtPOI = dynamic_cast<TProfile*> 
	  (outputListHistos->FindObject("Second_FlowPro_VPtPOI_LYZPROD"));

	//Set the histogram pointers and call Finish()
	if (pCommonHist && pCommonHistResults && pLYZHist2RP[0] && pLYZHist2POI[0] && 
	    pHistR0theta && pHistProReDenom && pHistProImDenom && pHistProVetaRP && 
	    pHistProVetaPOI && pHistProVPtRP && pHistProVPtPOI && pHistQsumforChi) {
	  this->SetCommonHists(pCommonHist);
	  this->SetCommonHistsRes(pCommonHistResults);
	  this->SetHist2RP(pLYZHist2RP);
	  this->SetHist2POI(pLYZHist2POI);
	  this->SetHistR0theta(pHistR0theta);
	  this->SetHistProReDenom(pHistProReDenom);
	  this->SetHistProImDenom(pHistProImDenom);
	  this->SetHistProVetaRP(pHistProVetaRP);
	  this->SetHistProVetaPOI(pHistProVetaPOI);
	  this->SetHistProVPtRP(pHistProVPtRP);
	  this->SetHistProVPtPOI(pHistProVPtPOI);
	  this->SetHistQsumforChi(pHistQsumforChi);
	} 
	else { 
	  cout<<"WARNING: Histograms needed to run Finish() secondrun (PROD) are not accessable!"<<endl; 
	}
      }
    } //secondrun
    //outputListHistos->Print(); 
    delete [] pLYZHist1;
    delete [] pLYZHist2RP;
    delete [] pLYZHist2POI;
  } //listhistos
  else { 
    cout << "histogram list pointer is empty in method AliFlowAnalysisWithLeeYangZeros::GetOutputHistograms() " << endl;}
}

  //-----------------------------------------------------------------------     
 Bool_t AliFlowAnalysisWithLeeYangZeros::Finish() 
{
  //finish method
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::Finish()****"<<endl; 
  
  //define variables for both runs
  Double_t  dJ01 = 2.405; 
  Int_t iNtheta = AliFlowLYZConstants::GetMaster()->GetNtheta();

  //set the event number
  SetEventNumber((int)fCommonHists->GetHistQ()->GetEntries());
  
  //Get multiplicity for RP selection
  Double_t dMultRP = fCommonHists->GetHistMultRP()->GetMean();
  if (fDebug) cout<<"The average multiplicity is "<<dMultRP<<endl;  

  //set the sum of Q vectors
  fQsum->Set(fHistQsumforChi->GetBinContent(1),fHistQsumforChi->GetBinContent(2));
  SetQ2sum(fHistQsumforChi->GetBinContent(3));  
    
  if (fFirstRun){

    //define variables for the first run
    Double_t  dR0 = 0;
    Double_t  dVtheta = 0; 
    Double_t  dv = 0;
    Double_t  dV = 0; 

    //reset histograms in case of merged output files
    fHistR0theta->Reset();
    fHistVtheta->Reset();
    
    for (Int_t theta=0;theta<iNtheta;theta++)
      {
	//get the first minimum r0
	fHist1[theta]->FillGtheta();
	dR0 = fHist1[theta]->GetR0();
		    	   	   
	//calculate integrated flow
	if (!TMath::AreEqualAbs(dR0, 0., 1e-100)) { dVtheta = dJ01/dR0; }
	else { cout<<"r0 is not found! Leaving LYZ analysis."<<endl; return kFALSE; }

	//for estimating systematic error resulting from d0
	Double_t dBinsize =0.;
	if (fUseSum){ dBinsize = (AliFlowLYZConstants::GetMaster()->GetMaxSUM())/(AliFlowLYZConstants::GetMaster()->GetNbins());}
	else { dBinsize = (AliFlowLYZConstants::GetMaster()->GetMaxPROD())/(AliFlowLYZConstants::GetMaster()->GetNbins());}
	Double_t dVplus = -1.;
	Double_t dVmin  = -1.;
	if (!TMath::AreEqualAbs(dR0+dBinsize, 0., 1e-100)) {dVplus = dJ01/(dR0+dBinsize);}
	if (!TMath::AreEqualAbs(dR0-dBinsize, 0., 1e-100)) {dVmin = dJ01/(dR0-dBinsize);}
	//convert V to v 
	Double_t dvplus = -1.;
	Double_t dvmin= -1.;
	if (fUseSum){
	  //for SUM: V=v because the Q-vector is scaled by 1/M
	  dv = dVtheta;
	  dvplus = dVplus;
	  dvmin = dVmin; }
	else {
	  //for PRODUCT: v=V/M
	  if (!TMath::AreEqualAbs(dMultRP, 0., 1e-100)){
	    dv = dVtheta/dMultRP;
	    dvplus = dVplus/dMultRP;
	    dvmin = dVmin/dMultRP; }}

	if (fDebug) cout<<"dv = "<<dv<<" and dvplus = "<<dvplus<< " and dvmin = "<<dvmin<<endl;
	     
	//fill the histograms
	fHistR0theta->SetBinContent(theta+1,dR0);
	fHistR0theta->SetBinError(theta+1,0.0);
	fHistVtheta->SetBinContent(theta+1,dVtheta);
	fHistVtheta->SetBinError(theta+1,0.0);
	//get average value of fVtheta = fV
	dV += dVtheta;
      } //end of loop over theta

    //get average value of fVtheta = fV
    dV /=iNtheta;
    if (!fUseSum) { if (dMultRP!=0.){dV /=dMultRP;}} //scale with multiplicity for PRODUCT
    
    //sigma2 and chi 
    Double_t  dSigma2 = 0;
    Double_t  dChi= 0;
    if (fEventNumber!=0) {
      *fQsum /= fEventNumber;
      //cout<<"fQsum is "<<fQsum->X()<<" "<<fQsum->Y()<<endl; 
      fQ2sum /= fEventNumber;
      //cout<<"fQ2sum is "<<fQ2sum<<endl; 
      dSigma2 = fQ2sum - TMath::Power(fQsum->X(),2.) - TMath::Power(fQsum->Y(),2.) - TMath::Power(dV,2.);  //BP eq. 62
      //cout<<"dSigma2 is "<<dSigma2<<endl; 
      if (dSigma2>0) dChi = dV/TMath::Sqrt(dSigma2);
      else dChi = -1.;
      fCommonHistsRes->FillChi(dChi);
  
      cout<<"*************************************"<<endl;
      cout<<"*************************************"<<endl;
      cout<<"      Integrated flow from           "<<endl;
      if (fUseSum) {
	cout<<"       Lee-Yang Zeroes SUM          "<<endl;}
      else {
	cout<<"     Lee-Yang Zeroes PRODUCT      "<<endl;}
      cout<<endl;
      cout<<"Chi = "<<dChi<<endl;
      //cout<<endl;
    }
	   
    // recalculate statistical errors on integrated flow
    //combining 5 theta angles to 1 relative error BP eq. 89
    Double_t dRelErr2comb = 0.;
    Int_t iEvts = fEventNumber; 
    if (iEvts!=0) {
      for (Int_t theta=0;theta<iNtheta;theta++){
	Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
				       TMath::Cos(dTheta));
	Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
				      TMath::Cos(dTheta));
	dRelErr2comb += (1/(2*iEvts*(dJ01*dJ01)*TMath::BesselJ1(dJ01)*
			  TMath::BesselJ1(dJ01)))*
	  (dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2)) + 
	   dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2)));
      }
      dRelErr2comb /= iNtheta;
    }
    Double_t dRelErrcomb = TMath::Sqrt(dRelErr2comb);

    //copy content of profile into TH1D and add error
    Double_t dv2pro = dV;   //in the case that fv is equal to fV
    Double_t dv2Err = dv2pro*dRelErrcomb ; 
    cout<<"dV = "<<dv2pro<<" +- "<<dv2Err<<endl;
    cout<<endl;
    cout<<"*************************************"<<endl;
    cout<<"*************************************"<<endl;
    fCommonHistsRes->FillIntegratedFlow(dv2pro, dv2Err);  


    if (fDebug) cout<<"****histograms filled****"<<endl;  
    
    return kTRUE;
  }  //firstrun
   
 
  else {  //second run to calculate differential flow
   
    //declare variables for the second run
    TComplex cDenom, cNumerRP, cNumerPOI, cDtheta;
    Int_t m = 1;
    TComplex i = TComplex::I();
    Double_t dBesselRatio[3] = {1., 1.202, 2.69};
    Int_t iNbinsPt = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
    Int_t iNbinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
         
    Double_t dR0 = 0.; 
    Double_t dVtheta = 0.;
    Double_t dV = 0.;
    Double_t dReDenom = 0.;
    Double_t dImDenom = 0.; 

    //reset histograms in case of merged output files
    if (fUseSum) { 
      fHistReDtheta->Reset();
      fHistImDtheta->Reset();
    }
    fHistProVetaRP ->Reset(); 
    fHistProVetaPOI->Reset(); 
    fHistProVPtRP  ->Reset(); 
    fHistProVPtPOI ->Reset(); 

    //scale fHistR0theta by the number of merged files to undo the merging
    if (!fHistR0theta) { cout<<"Hist pointer R0theta in file does not exist"<<endl; }
    else {
      Int_t iEntries = (int)fHistR0theta->GetEntries();
      if (iEntries > iNtheta){
	//for each individual file fHistR0theta has iNtheta entries
	Int_t iFiles = iEntries/iNtheta;
	cout<<iFiles<<" files were merged!"<<endl;
	fHistR0theta->Scale(1./iFiles);
      }

      //loop over theta
      for (Int_t theta=0;theta<iNtheta;theta++)  { 
      	dR0 = fHistR0theta->GetBinContent(theta+1); //histogram starts at bin 1
	if (fDebug) cerr<<"dR0 = "<<dR0<<endl;
	if (!TMath::AreEqualAbs(dR0, 0., 1e-100)) dVtheta = dJ01/dR0;
	dV += dVtheta; 
	// BP Eq. 9  -> Vn^theta = j01/r0^theta
	if (!fHistProReDenom || !fHistProImDenom) {
	  cout << "Hist pointer fDenom in file does not exist" <<endl;
	  cout<< "Leaving LYZ second pass analysis!" <<endl;
	  return kFALSE;
	} else {
	  dReDenom = fHistProReDenom->GetBinContent(theta+1);
	  dImDenom = fHistProImDenom->GetBinContent(theta+1);
	  cDenom(dReDenom,dImDenom);
	  
	  //for new method and use by others (only with the sum generating function):
	  if (fUseSum) {
	    cDtheta = dR0*cDenom/dJ01;
	    fHistReDtheta->SetBinContent(theta+1,cDtheta.Re()); 
	    fHistReDtheta->SetBinError(theta+1,0.0);
	    fHistImDtheta->SetBinContent(theta+1,cDtheta.Im());
	    fHistImDtheta->SetBinError(theta+1,0.0);
	  }
	 
	  cDenom *= TComplex::Power(i, m-1);
	  
	  //v as a function of eta for RP selection
	  for (Int_t be=1;be<=iNbinsEta;be++)  {
	    Double_t dReRatioRP = 0.0;
	    Double_t dEtaRP = fHist2RP[theta]->GetBinCenter(be);
	    cNumerRP = fHist2RP[theta]->GetNumerEta(be);
	    if (cNumerRP.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cNumerRP is zero in Finish()"<<endl;
	    }
	    else if (TMath::AreEqualAbs(cDenom.Rho(), 0, 1e-100)) {
	      if (fDebug) cerr<<"WARNING: modulus of cDenom is zero"<<endl;
	    }
	    else {
	      dReRatioRP = (cNumerRP/cDenom).Re();
	    }
	    Double_t dVetaRP = dBesselRatio[m-1]*dReRatioRP*dVtheta; //BP eq. 12
	    fHistProVetaRP->Fill(dEtaRP,dVetaRP);
	  } //loop over bins be
	 

	  //v as a function of eta for POI selection
	  for (Int_t be=1;be<=iNbinsEta;be++)  {
	    Double_t dReRatioPOI = 0.0;
	    Double_t dEtaPOI = fHist2POI[theta]->GetBinCenter(be);
	    cNumerPOI = fHist2POI[theta]->GetNumerEta(be);
	    if (cNumerPOI.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cNumerPOI is zero in Finish()"<<endl;
	    }
	    else if (cDenom.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cDenom is zero"<<endl;
	    }
	    else {
	      dReRatioPOI = (cNumerPOI/cDenom).Re();
	    }
	    Double_t dVetaPOI = dBesselRatio[m-1]*dReRatioPOI*dVtheta; //BP eq. 12
	    fHistProVetaPOI->Fill(dEtaPOI,dVetaPOI);
	  } //loop over bins be


	  //v as a function of Pt for RP selection
	  for (Int_t bp=1;bp<=iNbinsPt;bp++)  {
	    Double_t dReRatioRP = 0.0;
	    Double_t dPtRP = fHist2RP[theta]->GetBinCenterPt(bp);
	    cNumerRP = fHist2RP[theta]->GetNumerPt(bp);
	    if (cNumerRP.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cNumerRP is zero"<<endl;
	    }
	    else if (cDenom.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cDenom is zero"<<endl;
	    }
	    else {
	      dReRatioRP = (cNumerRP/cDenom).Re();
	    }
	    Double_t dVPtRP = dBesselRatio[m-1]*dReRatioRP*dVtheta; //BP eq. 12
	    fHistProVPtRP->Fill(dPtRP,dVPtRP);
	  } //loop over bins bp


	  //v as a function of Pt for POI selection
	  for (Int_t bp=1;bp<=iNbinsPt;bp++)  {
	    Double_t dReRatioPOI = 0.0;
	    Double_t dPtPOI = fHist2POI[theta]->GetBinCenterPt(bp);
	    cNumerPOI = fHist2POI[theta]->GetNumerPt(bp);
	    if (cNumerPOI.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cNumerPOI is zero"<<endl;
	    }
	    else if (cDenom.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cDenom is zero"<<endl;
	    }
	    else {
	      dReRatioPOI = (cNumerPOI/cDenom).Re();
	    }
	    Double_t dVPtPOI = dBesselRatio[m-1]*dReRatioPOI*dVtheta; //BP eq. 12
	    fHistProVPtPOI->Fill(dPtPOI,dVPtPOI);
	  } //loop over bins bp

	}
      }

    }//end of loop over theta

    //calculate the average of fVtheta = fV
    dV /= iNtheta;
    if (!fUseSum) { if (dMultRP!=0.) { dV /=dMultRP; }} //scale by the multiplicity for PRODUCT
    if (TMath::AreEqualAbs(dV, 0., 1e-100)) { cout<<"dV = 0! Leaving LYZ analysis."<<endl; return kFALSE; }

    //sigma2 and chi (for statistical error calculations)
    Double_t  dSigma2 = 0;
    Double_t  dChi= 0;
    if (fEventNumber!=0) {
      *fQsum /= fEventNumber;
      //cerr<<"fQsum.X() = "<<fQsum.X()<<endl;
      //cerr<<"fQsum.Y() = "<<fQsum.Y()<<endl;
      fQ2sum /= fEventNumber;
      //cout<<"fQ2sum = "<<fQ2sum<<endl;
      dSigma2 = fQ2sum - TMath::Power(fQsum->X(),2.) - TMath::Power(fQsum->Y(),2.) - TMath::Power(dV,2.);  //BP eq. 62
      if (dSigma2>0) dChi = dV/TMath::Sqrt(dSigma2);
      else dChi = -1.;
      fCommonHistsRes->FillChi(dChi);

      // recalculate statistical errors on integrated flow
      //combining 5 theta angles to 1 relative error BP eq. 89
      Double_t dRelErr2comb = 0.;
      Int_t iEvts = fEventNumber; 
      if (iEvts!=0) {
	for (Int_t theta=0;theta<iNtheta;theta++){
	Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
				       TMath::Cos(dTheta));
	Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
				      TMath::Cos(dTheta));
	dRelErr2comb += (1/(2*iEvts*(dJ01*dJ01)*TMath::BesselJ1(dJ01)*
			  TMath::BesselJ1(dJ01)))*
	  (dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2)) + 
	   dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2)));
	}
	dRelErr2comb /= iNtheta;
      }
      Double_t dRelErrcomb = TMath::Sqrt(dRelErr2comb);
      Double_t dVErr = dV*dRelErrcomb ; 
      fCommonHistsRes->FillIntegratedFlow(dV,dVErr); 

      cout<<"*************************************"<<endl;
      cout<<"*************************************"<<endl;
      cout<<"      Integrated flow from           "<<endl;
      if (fUseSum) {
	cout<<"       Lee-Yang Zeroes SUM          "<<endl;}
      else {
	cout<<"     Lee-Yang Zeroes PRODUCT      "<<endl;}
      cout<<endl;
      cout<<"Chi = "<<dChi<<endl;
      cout<<"dV = "<<dV<<" +- "<<dVErr<<endl;
      //cout<<endl;
    }
	     
    //copy content of profile into TH1D and add error and fill the AliFlowCommonHistResults

    //v as a function of eta for RP selection
    for(Int_t b=0;b<iNbinsEta;b++) {
      Double_t dv2pro  = fHistProVetaRP->GetBinContent(b);
      Double_t dNprime = fCommonHists->GetEntriesInEtaBinRP(b);  
      Double_t dErrdifcomb = 0.;  //set error to zero
      Double_t dErr2difcomb = 0.; //set error to zero
      //calculate error
      if (!TMath::AreEqualAbs(dNprime, 0., 1e-100)) { 
	for (Int_t theta=0;theta<iNtheta;theta++) {
	  Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	  Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
					   TMath::Cos(dTheta));
	  Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
					  TMath::Cos(dTheta));
	  dErr2difcomb += (TMath::Cos(dTheta)/(4*dNprime*TMath::BesselJ1(dJ01)*
						 TMath::BesselJ1(dJ01)))*
	    ((dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2))) - 
	     (dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2))));
	} //loop over theta
      } 
      
      if (!TMath::AreEqualAbs(dErr2difcomb, 0., 1e-100)) {
	dErr2difcomb /= iNtheta;
	dErrdifcomb = TMath::Sqrt(dErr2difcomb);
      }
      else {dErrdifcomb = 0.;}
      //fill TH1D
      fCommonHistsRes->FillDifferentialFlowEtaRP(b, dv2pro, dErrdifcomb); 
    } //loop over bins b


    //v as a function of eta for POI selection
    for(Int_t b=0;b<iNbinsEta;b++) {
      Double_t dv2pro  = fHistProVetaPOI->GetBinContent(b);
      Double_t dNprime = fCommonHists->GetEntriesInEtaBinPOI(b);  
      Double_t dErrdifcomb = 0.;  //set error to zero
      Double_t dErr2difcomb = 0.; //set error to zero
      //calculate error
      if (dNprime!=0.) { 
	for (Int_t theta=0;theta<iNtheta;theta++) {
	  Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	  Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
					   TMath::Cos(dTheta));
	  Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
					  TMath::Cos(dTheta));
	  dErr2difcomb += (TMath::Cos(dTheta)/(4*dNprime*TMath::BesselJ1(dJ01)*
					       TMath::BesselJ1(dJ01)))*
	    ((dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2))) - 
	     (dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2))));
	} //loop over theta
      } 
      
      if (dErr2difcomb!=0.) {
	dErr2difcomb /= iNtheta;
	dErrdifcomb = TMath::Sqrt(dErr2difcomb);
      }
      else {dErrdifcomb = 0.;}
      //fill TH1D
      fCommonHistsRes->FillDifferentialFlowEtaPOI(b, dv2pro, dErrdifcomb); 
    } //loop over bins b
    
    
    
    //v as a function of Pt for RP selection
    TH1F* fHistPtRP = fCommonHists->GetHistPtRP(); //for calculating integrated flow
    Double_t dVRP = 0.;
    Double_t dSum = 0.;
    Double_t dErrV =0.;
    for(Int_t b=0;b<iNbinsPt;b++) {
      Double_t dv2pro  = fHistProVPtRP->GetBinContent(b);
      Double_t dNprime = fCommonHists->GetEntriesInPtBinRP(b);  
      Double_t dErrdifcomb = 0.;  //set error to zero
      Double_t dErr2difcomb = 0.; //set error to zero
      //calculate error
      if (dNprime!=0.) { 
	for (Int_t theta=0;theta<iNtheta;theta++) {
	  Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	  Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
					   TMath::Cos(dTheta));
	  Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
					  TMath::Cos(dTheta));
	  dErr2difcomb += (TMath::Cos(dTheta)/(4*dNprime*TMath::BesselJ1(dJ01)*
					       TMath::BesselJ1(dJ01)))*
	    ((dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2))) - 
	     (dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2))));
	} //loop over theta
      } 
      
      if (dErr2difcomb!=0.) {
	dErr2difcomb /= iNtheta;
	dErrdifcomb = TMath::Sqrt(dErr2difcomb);
	//cerr<<"dErrdifcomb = "<<dErrdifcomb<<endl;
      }
      else {dErrdifcomb = 0.;}
      
      //fill TH1D
      fCommonHistsRes->FillDifferentialFlowPtRP(b, dv2pro, dErrdifcomb); 
      //calculate integrated flow for RP selection
      if (fHistPtRP){
	Double_t dYieldPt = fHistPtRP->GetBinContent(b);
	dVRP += dv2pro*dYieldPt;
	dSum +=dYieldPt;
	dErrV += dYieldPt*dYieldPt*dErrdifcomb*dErrdifcomb;
      } else { cout<<"fHistPtRP is NULL"<<endl; }
    } //loop over bins b

    if (!TMath::AreEqualAbs(dSum, 0., 1e-100)) {
      dVRP /= dSum; //the pt distribution should be normalised
      dErrV /= (dSum*dSum);
      dErrV = TMath::Sqrt(dErrV);
    }
    cout<<"dV(RP) = "<<dVRP<<" +- "<<dErrV<<endl;
    //cout<<endl;
    fCommonHistsRes->FillIntegratedFlowRP(dVRP,dErrV);

             
    //v as a function of Pt for POI selection 
    TH1F* fHistPtPOI = fCommonHists->GetHistPtPOI(); //for calculating integrated flow
    Double_t dVPOI = 0.;
    dSum = 0.;
    dErrV =0.;

    for(Int_t b=0;b<iNbinsPt;b++) {
      Double_t dv2pro = fHistProVPtPOI->GetBinContent(b);
      Double_t dNprime = fCommonHists->GetEntriesInPtBinPOI(b);   
      Double_t dErrdifcomb = 0.;  //set error to zero
      Double_t dErr2difcomb = 0.; //set error to zero
      //calculate error
      if (dNprime!=0.) { 
	for (Int_t theta=0;theta<iNtheta;theta++) {
	  Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	  Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
					   TMath::Cos(dTheta));
	  Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
					  TMath::Cos(dTheta));
	  dErr2difcomb += (TMath::Cos(dTheta)/(4*dNprime*TMath::BesselJ1(dJ01)*
						 TMath::BesselJ1(dJ01)))*
	    ((dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2))) - 
	     (dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2))));
	} //loop over theta
      } 
      
      if (dErr2difcomb!=0.) {
	dErr2difcomb /= iNtheta;
	dErrdifcomb = TMath::Sqrt(dErr2difcomb);
	//cerr<<"dErrdifcomb = "<<dErrdifcomb<<endl;
      }
      else {dErrdifcomb = 0.;}
	  
      //fill TH1D
      fCommonHistsRes->FillDifferentialFlowPtPOI(b, dv2pro, dErrdifcomb); 
      //calculate integrated flow for POI selection
      if (fHistPtPOI){
	Double_t dYieldPt = fHistPtPOI->GetBinContent(b);
	dVPOI += dv2pro*dYieldPt;
	dSum +=dYieldPt;
	dErrV += dYieldPt*dYieldPt*dErrdifcomb*dErrdifcomb;
      } else { cout<<"fHistPtPOI is NULL"<<endl; }
    } //loop over bins b

    if (dSum != 0.) {
      dVPOI /= dSum; //the pt distribution should be normalised
      dErrV /= (dSum*dSum);
      dErrV = TMath::Sqrt(dErrV);
    }
    cout<<"dV(POI) = "<<dVPOI<<" +- "<<dErrV<<endl;
    cout<<endl;
    cout<<"*************************************"<<endl;
    cout<<"*************************************"<<endl;
    fCommonHistsRes->FillIntegratedFlowPOI(dVPOI,dErrV);

  } //secondrun
   
  //cout<<"----LYZ analysis finished....----"<<endl<<endl;

  return kTRUE;
}


//-----------------------------------------------------------------------
 
 Bool_t AliFlowAnalysisWithLeeYangZeros::FillFromFlowEvent(AliFlowEventSimple* anEvent) 
{ 
  // Get event quantities from AliFlowEvent for all particles

  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::FillFromFlowEvent()****"<<endl;
   
  if (!anEvent){
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }
   
  //define variables
  TComplex cExpo, cGtheta, cGthetaNew, cZ;
  Int_t iNtheta = AliFlowLYZConstants::GetMaster()->GetNtheta();
  Int_t iNbins = AliFlowLYZConstants::GetMaster()->GetNbins();
  
  //calculate flow
  Double_t dOrder = 2.;
      
  //get the Q vector 
  AliFlowVector vQ = anEvent->GetQ();
  //weight by the multiplicity
  Double_t dQX = 0;
  Double_t dQY = 0;
  if (!TMath::AreEqualAbs(vQ.GetMult(), 0., 1e-100)) {
    dQX = vQ.X()/vQ.GetMult(); 
    dQY = vQ.Y()/vQ.GetMult();
  }
  vQ.Set(dQX,dQY);

  //for chi calculation:
  *fQsum += vQ;
  fHistQsumforChi->SetBinContent(1,fQsum->X());
  fHistQsumforChi->SetBinContent(2,fQsum->Y());
  fQ2sum += vQ.Mod2();
  fHistQsumforChi->SetBinContent(3,fQ2sum);
  
  for (Int_t theta=0;theta<iNtheta;theta++) {
    Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi()/dOrder; 
	  
    //calculate dQtheta = cos(dOrder*(fPhi-dTheta);the projection of the Q vector on the reference direction dTheta
    Double_t dQtheta = GetQtheta(vQ, dTheta);
	     	   
    for (Int_t bin=1;bin<=iNbins;bin++) {
      Double_t dR = fHist1[theta]->GetBinCenter(bin); //bincentre of bins in histogram  //FIXED???
      if (fUseSum) {
	//calculate the sum generating function
	cExpo(0.,dR*dQtheta); //Re=0 ; Im=dR*dQtheta
	cGtheta = TComplex::Exp(cExpo);
      }
      else {
	//calculate the product generating function
	cGtheta = GetGrtheta(anEvent, dR, dTheta);  
	if (cGtheta.Rho2() > 100.) break;
      }
      //fill real and imaginary part of cGtheta
      fHist1[theta]->Fill(dR,cGtheta);    
    } //loop over bins
  } //loop over theta 
    
  return kTRUE;
          
}

 //-----------------------------------------------------------------------   
 Bool_t AliFlowAnalysisWithLeeYangZeros::SecondFillFromFlowEvent(AliFlowEventSimple* anEvent) 
{ 
  //for differential flow

  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::SecondFillFromFlowEvent()****"<<endl;
    
  if (!anEvent){
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }
   
  //define variables
  TComplex cExpo, cDenom, cNumerRP, cNumerPOI, cCosTermComplex;
  Double_t dR0 = 0.;
  Double_t dCosTermRP = 0.;
  Double_t dCosTermPOI = 0.;
  Double_t m = 1.;
  Double_t dOrder = 2.;
  Int_t iNtheta = AliFlowLYZConstants::GetMaster()->GetNtheta();
     
  //get the Q vector 
  AliFlowVector vQ = anEvent->GetQ();
  //weight by the multiplicity
  Double_t dQX = 0.;
  Double_t dQY = 0.;
  if (vQ.GetMult() != 0) {
    dQX = vQ.X()/vQ.GetMult(); 
    dQY = vQ.Y()/vQ.GetMult();
  }
  vQ.Set(dQX,dQY); 
              
  //for chi calculation:
  *fQsum += vQ;
  fHistQsumforChi->SetBinContent(1,fQsum->X());
  fHistQsumforChi->SetBinContent(2,fQsum->Y());
  fQ2sum += vQ.Mod2();
  fHistQsumforChi->SetBinContent(3,fQ2sum);

  for (Int_t theta=0;theta<iNtheta;theta++)
    {
      Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi()/dOrder;   

      //calculate dQtheta = cos(dOrder*(dPhi-dTheta);the projection of the Q vector on the reference direction dTheta	  
      Double_t dQtheta = GetQtheta(vQ, dTheta);
        	 
      //denominator for differential v
      if (fHistR0theta) {
	dR0 = fHistR0theta->GetBinContent(theta+1);
      }
      else { cout <<"pointer fHistR0theta does not exist" << endl;
      }
      
      if (fUseSum) //sum generating function
	{
	  cExpo(0.,dR0*dQtheta);
	  cDenom = dQtheta*(TComplex::Exp(cExpo)); //BP eq 12
	  //loop over tracks in event
	  Int_t iNumberOfTracks = anEvent->NumberOfTracks();
	  for (Int_t i=0;i<iNumberOfTracks;i++)  {
	    AliFlowTrackSimple*  pTrack = anEvent->GetTrack(i);
	    if (pTrack) {
	      Double_t dEta = pTrack->Eta();
	      Double_t dPt = pTrack->Pt();
	      Double_t dPhi = pTrack->Phi();
	      if (pTrack->InRPSelection()) { // RP selection
		dCosTermRP = cos(m*dOrder*(dPhi-dTheta));
		cNumerRP = dCosTermRP*(TComplex::Exp(cExpo));
		if (cNumerRP.Rho()==0) { cerr<<"WARNING: modulus of cNumerRP is zero in SecondFillFromFlowEvent"<<endl;}
		if (fDebug) { cerr<<"modulus of cNumerRP is "<<cNumerRP.Rho()<<endl;}
		if (fHist2RP[theta]) {
		  fHist2RP[theta]->Fill(dEta,dPt,cNumerRP); 
		}
	      }
	      if (pTrack->InPOISelection()) { //POI selection
		dCosTermPOI = cos(m*dOrder*(dPhi-dTheta));
		cNumerPOI = dCosTermPOI*(TComplex::Exp(cExpo));
		if (cNumerPOI.Rho()==0) { cerr<<"WARNING: modulus of cNumerPOI is zero in SecondFillFromFlowEvent"<<endl;}
		if (fDebug) { cerr<<"modulus of cNumerPOI is "<<cNumerPOI.Rho()<<endl;}
		if (fHist2POI[theta]) {
		  fHist2POI[theta]->Fill(dEta,dPt,cNumerPOI); 
		}
	      }
	    } //if track
	    else {cerr << "no particle!!!"<<endl;}
	  } //loop over tracks
	} //sum
      else {    //product generating function
	cDenom = GetDiffFlow(anEvent, dR0, theta); 
      }//product
      
      if (fHistProReDenom && fHistProImDenom) {
	fHistProReDenom->Fill(theta,cDenom.Re());               //fill the real part of fDenom
	fHistProImDenom->Fill(theta,cDenom.Im());               //fill the imaginary part of fDenom
      }
      else { cout << "Pointers to cDenom  mising" << endl;}
            
    }//end of loop over theta
  
  return kTRUE;
    
}
 //-----------------------------------------------------------------------   
 Double_t AliFlowAnalysisWithLeeYangZeros::GetQtheta(AliFlowVector aQ, Double_t aTheta) 
{
  //calculate Qtheta. Qtheta is the sum over all particles of cos(dOrder*(dPhi-dTheta)) BP eq. 3
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::GetQtheta()****"<<endl;

  Double_t dQtheta = 0.;
  Double_t dOrder = 2.;
  
  dQtheta = aQ.X()*cos(dOrder*aTheta)+aQ.Y()*sin(dOrder*aTheta);

  return dQtheta;
 
}
 

//-----------------------------------------------------------------------   
TComplex AliFlowAnalysisWithLeeYangZeros::GetGrtheta(AliFlowEventSimple* const anEvent, Double_t aR, Double_t aTheta) 
{
  // Product Generating Function for LeeYangZeros method
  // PG Eq. 3 (J. Phys. G Nucl. Part. Phys 30 S1213 (2004))
  
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::GetGrtheta()****"<<endl;
  
  
  TComplex cG = TComplex::One();
  Double_t dOrder =  2.;
  Double_t dWgt = 1.;
  //Double_t dWgt = 1./anEvent->GetEventNSelTracksRP(); //weight with the multiplicity
    
  Int_t iNumberOfTracks = anEvent->NumberOfTracks();
  
  for (Int_t i=0;i<iNumberOfTracks;i++) //loop over tracks in event
    {
      AliFlowTrackSimple* pTrack = anEvent->GetTrack(i) ; 
      if (pTrack){
	if (pTrack->InRPSelection()) {
	  Double_t dPhi = pTrack->Phi();
	  Double_t dGIm = aR * dWgt*cos(dOrder*(dPhi - aTheta));
	  TComplex cGi(1., dGIm);
	  cG *= cGi;     //product over all tracks
	}
      }
      else {cerr << "no particle pointer !!!"<<endl;}
    }//loop over tracks
  
  return cG;
  
} 


//-----------------------------------------------------------------------   
TComplex AliFlowAnalysisWithLeeYangZeros::GetDiffFlow(AliFlowEventSimple* const anEvent, Double_t aR0, Int_t theta) 
{
  // Sum for the denominator for diff. flow for the Product Generating Function for LeeYangZeros method
  // PG Eq. 9 (J. Phys. G Nucl. Part. Phys 30 S1213 (2004))
  // Also for v1 mixed harmonics: DF Eq. 5
  // It is the deriverative of Grtheta at r0 divided by Grtheta at r0
  
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::GetGrtheta()****"<<endl;
  
  TComplex cG = TComplex::One();
  TComplex cdGr0(0.,0.);
  Double_t dOrder =  2.;
  Double_t dWgt = 1.;
  //Double_t dWgt = 1./anEvent->GetEventNSelTracksRP(); //weight with the multiplicity

  Int_t iNumberOfTracks = anEvent->NumberOfTracks();
  
  Int_t iNtheta = AliFlowLYZConstants::GetMaster()->GetNtheta();
  Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi()/dOrder;
  
  //for the denominator (use all RP selected particles)
  for (Int_t i=0;i<iNumberOfTracks;i++) //loop over tracks in event
    {
      AliFlowTrackSimple* pTrack = anEvent->GetTrack(i) ;  
      if (pTrack){
	if (pTrack->InRPSelection()) {
	  Double_t dPhi = pTrack->Phi();
	  Double_t dCosTerm = dWgt*cos(dOrder*(dPhi - dTheta));
	  //GetGr0theta
	  Double_t dGIm = aR0 * dCosTerm;
	  TComplex cGi(1., dGIm);
	  TComplex cCosTermComplex(1., aR0*dCosTerm);
	  cG *= cGi;     //product over all tracks
	  //GetdGr0theta
	  cdGr0 +=(dCosTerm / cCosTermComplex);  //sum over all tracks
	}
      } //if particle
      else {cerr << "no particle!!!"<<endl;}
    }//loop over tracks
  
  //for the numerator
  for (Int_t i=0;i<iNumberOfTracks;i++) 
    {
      AliFlowTrackSimple* pTrack = anEvent->GetTrack(i) ;  
      if (pTrack){
	Double_t dEta = pTrack->Eta();
	Double_t dPt = pTrack->Pt();
	Double_t dPhi = pTrack->Phi();
	Double_t dCosTerm = cos(dOrder*(dPhi-dTheta));
	TComplex cCosTermComplex(1.,aR0*dCosTerm);
	//RP selection
	if (pTrack->InRPSelection()) {
	  TComplex cNumerRP = cG*dCosTerm/cCosTermComplex;  //PG Eq. 9
	  fHist2RP[theta]->Fill(dEta,dPt,cNumerRP);  
	}
	//POI selection
	if (pTrack->InPOISelection()) {
	  TComplex cNumerPOI = cG*dCosTerm/cCosTermComplex;  //PG Eq. 9
	  fHist2POI[theta]->Fill(dEta,dPt,cNumerPOI);  
	}
      } //if particle
      else {cerr << "no particle pointer!!!"<<endl;}
    }//loop over tracks
  
  TComplex cDenom = cG*cdGr0;  
  return cDenom;
  
} 

//----------------------------------------------------------------------- 

