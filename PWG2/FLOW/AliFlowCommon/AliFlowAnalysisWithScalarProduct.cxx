/*************************************************************************
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

#define AliFlowAnalysisWithScalarProduct_cxx
 
#include "Riostream.h"
#include "TFile.h"      
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TH1D.h"
#include "TH2D.h"

#include "AliFlowCommonConstants.h"
#include "AliFlowEventSimple.h"
#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowAnalysisWithScalarProduct.h"

//////////////////////////////////////////////////////////////////////////////
// AliFlowAnalysisWithScalarProduct:
// Description: 
// Maker to analyze Flow with the Scalar product method.
//
// authors: N. van der Kolk (kolk@nikhef.nl), A. Bilandzic (anteb@nikhef.nl)
//////////////////////////////////////////////////////////////////////////////

ClassImp(AliFlowAnalysisWithScalarProduct)

  //-----------------------------------------------------------------------
  AliFlowAnalysisWithScalarProduct::AliFlowAnalysisWithScalarProduct():
   fEventNumber(0),
   fDebug(kFALSE),
   fApplyCorrectionForNUA(kFALSE),
   fHarmonic(2),
   fTotalQvector(NULL),
   fRelDiffMsub(1.),
   fWeightsList(NULL),
   fUsePhiWeights(kFALSE),
   fPhiWeightsSub0(NULL),
   fPhiWeightsSub1(NULL),
   fHistProFlags(NULL),
   fHistProUQetaRP(NULL),
   fHistProUQetaPOI(NULL),
   fHistProUQetaAllEventsPOI(NULL),
   fHistProUQPtRP(NULL),
   fHistProUQPtPOI(NULL),
   fHistProUQPtAllEventsPOI(NULL),
   fHistProQNorm(NULL),
   fHistProQaQb(NULL),
   fHistProQaQbVsM(NULL),
   fHistProQaQbNorm(NULL),
   fHistProQaQbReImNorm(NULL),
   fHistProNonIsotropicTermsQ(NULL),
   fHistSumOfLinearWeights(NULL),
   fHistSumOfQuadraticWeights(NULL),
   fHistProUQQaQbPtRP(NULL),
   fHistProUQQaQbEtaRP(NULL),
   fHistProUQQaQbPtPOI(NULL),
   fHistProUQQaQbEtaPOI(NULL),
   fCommonHistsSP(NULL),
   fCommonHistsResSP(NULL),
   fCommonHistsmuQ(NULL),
   fHistQNorm(NULL),
   fHistQaQb(NULL),
   fHistQaQbNorm(NULL),
   fHistQNormvsQaQbNorm(NULL),
   fHistQaQbCos(NULL),
   fHistResolution(NULL),
   fHistQaNorm(NULL),
   fHistQaNormvsMa(NULL),
   fHistQbNorm(NULL),
   fHistQbNormvsMb(NULL),
   fHistMavsMb(NULL),
   fHistList(NULL)
{
  // Constructor.
  fWeightsList = new TList();
  fHistList = new TList();
  fHistList->SetOwner(kTRUE);
  
  // Total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb"
  fTotalQvector = new TString("QaQb");
  
  // Initialize arrays:
  for(Int_t i=0;i<3;i++)
  {
   fHistSumOfWeightsPtRP[i] = NULL;
   fHistSumOfWeightsEtaRP[i] = NULL;
   fHistSumOfWeightsPtPOI[i] = NULL;
   fHistSumOfWeightsEtaPOI[i] = NULL;
  }
  for(Int_t rp=0;rp<2;rp++)
  {
   for(Int_t pe=0;pe<2;pe++)
   {
    for(Int_t sc=0;sc<2;sc++)
    {
     fHistProNonIsotropicTermsU[rp][pe][sc] = NULL;
    }
   } 
  }
}
 //-----------------------------------------------------------------------
 AliFlowAnalysisWithScalarProduct::~AliFlowAnalysisWithScalarProduct() 
 {
   //destructor
   delete fWeightsList;
   delete fHistList;
 }
 
//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file

  TFile *output = new TFile(outputFileName->Data(),"RECREATE");
  //output->WriteObject(fHistList, "cobjSP","SingleKey");
  fHistList->SetName("cobjSP");
  fHistList->SetOwner(kTRUE);
  fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
  delete output;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file

  TFile *output = new TFile(outputFileName.Data(),"RECREATE");
  //output->WriteObject(fHistList, "cobjSP","SingleKey");
  fHistList->SetName("cobjSP");
  fHistList->SetOwner(kTRUE);
  fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
  delete output;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::WriteHistograms(TDirectoryFile *outputFileName)
{
 //store the final results in output .root file
 fHistList->SetName("cobjSP");
 fHistList->SetOwner(kTRUE);
 outputFileName->Add(fHistList);
 outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::Init() {

  //Define all histograms
  cout<<"---Analysis with the Scalar Product Method--- Init"<<endl;

  //save old value and prevent histograms from being added to directory
  //to avoid name clashes in case multiple analaysis objects are used
  //in an analysis

  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
 
  Int_t iNbinsPt   = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
  Double_t dPtMin  = AliFlowCommonConstants::GetMaster()->GetPtMin();	     
  Double_t dPtMax  = AliFlowCommonConstants::GetMaster()->GetPtMax();
  Int_t iNbinsEta  = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
  Double_t dEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();	     
  Double_t dEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();

  fHistProFlags = new TProfile("FlowPro_Flags_SP","Flow_Flags_SP",1,0,1,"s");
  fHistProFlags->GetXaxis()->SetBinLabel(1,"fApplyCorrectionForNUA");
  fHistList->Add(fHistProFlags);
  
  fHistProUQetaRP = new TProfile("FlowPro_UQetaRP_SP","Flow_UQetaRP_SP",iNbinsEta,dEtaMin,dEtaMax,"s");
  fHistProUQetaRP->SetXTitle("{eta}");
  fHistProUQetaRP->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQetaRP);

  fHistProUQetaPOI = new TProfile("FlowPro_UQetaPOI_SP","Flow_UQetaPOI_SP",iNbinsEta,dEtaMin,dEtaMax,"s");
  fHistProUQetaPOI->SetXTitle("{eta}");
  fHistProUQetaPOI->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQetaPOI);

  fHistProUQetaAllEventsPOI = new TProfile("FlowPro_UQetaAllEventsPOI_SP","FlowPro_UQetaAllEventsPOI_SP",iNbinsEta,dEtaMin,dEtaMax);
  fHistProUQetaAllEventsPOI->SetXTitle("{eta}");
  fHistProUQetaAllEventsPOI->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQetaAllEventsPOI);

  fHistProUQPtRP = new TProfile("FlowPro_UQPtRP_SP","Flow_UQPtRP_SP",iNbinsPt,dPtMin,dPtMax,"s");
  fHistProUQPtRP->SetXTitle("p_t (GeV)");
  fHistProUQPtRP->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQPtRP);

  fHistProUQPtPOI = new TProfile("FlowPro_UQPtPOI_SP","Flow_UQPtPOI_SP",iNbinsPt,dPtMin,dPtMax,"s");
  fHistProUQPtPOI->SetXTitle("p_t (GeV)");
  fHistProUQPtPOI->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQPtPOI);

  fHistProUQPtAllEventsPOI = new TProfile("FlowPro_UQPtAllEventsPOI_SP","FlowPro_UQPtAllEventsPOI_SP",iNbinsPt,dPtMin,dPtMax);
  fHistProUQPtAllEventsPOI->SetXTitle("p_t (GeV)");
  fHistProUQPtAllEventsPOI->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQPtAllEventsPOI);

  fHistProQNorm = new TProfile("FlowPro_QNorm_SP","FlowPro_QNorm_SP", 1, 0.5, 1.5,"s");
  fHistProQNorm ->SetYTitle("<|Qa+Qb|>");
  fHistList->Add(fHistProQNorm); 

  fHistProQaQb = new TProfile("FlowPro_QaQb_SP","FlowPro_QaQb_SP", 1, 0.5, 1.5,"s");
  fHistProQaQb->SetYTitle("<QaQb>");
  fHistList->Add(fHistProQaQb); 

  fHistProQaQbVsM = new TProfile("FlowPro_QaQbVsM_SP","FlowPro_QaQbVsM_SP",1000,0.,10000.); // to be improved (hardwired numbers)
  fHistProQaQbVsM->SetYTitle("#LTQ_{a}Q_{b}^{*}#GT");
  fHistProQaQbVsM->SetXTitle("M");
  fHistProQaQbVsM->Sumw2();
  fHistList->Add(fHistProQaQbVsM); 

  fHistProQaQbNorm = new TProfile("FlowPro_QaQbNorm_SP","FlowPro_QaQbNorm_SP", 1, 0.5, 1.5,"s");
  fHistProQaQbNorm->SetYTitle("<QaQb/MaMb>");
  fHistList->Add(fHistProQaQbNorm);
  
  fHistProQaQbReImNorm = new TProfile("FlowPro_QaQbReImNorm_SP","FlowPro_QaQbReImNorm_SP", 4, 0.5, 4.5,"s");
  fHistProQaQbReImNorm->GetXaxis()->SetBinLabel(1,"#LT#LTsin(#phi_{a})#GT#GT");
  fHistProQaQbReImNorm->GetXaxis()->SetBinLabel(2,"#LT#LTcos(#phi_{a})#GT#GT");
  fHistProQaQbReImNorm->GetXaxis()->SetBinLabel(3,"#LT#LTsin(#phi_{b})#GT#GT");
  fHistProQaQbReImNorm->GetXaxis()->SetBinLabel(4,"#LT#LTcos(#phi_{b})#GT#GT");
  fHistList->Add(fHistProQaQbReImNorm); 
  
  fHistProNonIsotropicTermsQ = new TProfile("FlowPro_NonIsotropicTermsQ_SP","FlowPro_NonIsotropicTermsQ_SP", 2, 0.5, 2.5,"s");
  if(!strcmp(fTotalQvector->Data(),"QaQb"))
  {
   fHistProNonIsotropicTermsQ->GetXaxis()->SetBinLabel(1,"#LT#LTsin(#phi_{a+b})#GT#GT");
   fHistProNonIsotropicTermsQ->GetXaxis()->SetBinLabel(2,"#LT#LTcos(#phi_{a+b})#GT#GT");
  } else if(!strcmp(fTotalQvector->Data(),"Qa"))
    {
     fHistProNonIsotropicTermsQ->GetXaxis()->SetBinLabel(1,"#LT#LTsin(#phi_{a})#GT#GT");
     fHistProNonIsotropicTermsQ->GetXaxis()->SetBinLabel(2,"#LT#LTcos(#phi_{a})#GT#GT");
    } else if(!strcmp(fTotalQvector->Data(),"Qb"))
      {
       fHistProNonIsotropicTermsQ->GetXaxis()->SetBinLabel(1,"#LT#LTsin(#phi_{b})#GT#GT");
       fHistProNonIsotropicTermsQ->GetXaxis()->SetBinLabel(2,"#LT#LTcos(#phi_{b})#GT#GT");    
      }
  fHistList->Add(fHistProNonIsotropicTermsQ); 
  
  TString rpPoi[2] = {"RP","POI"};
  TString ptEta[2] = {"Pt","Eta"};
  TString sinCos[2] = {"sin","cos"};
  Int_t nBinsPtEta[2] = {iNbinsPt,iNbinsEta};
  Double_t minPtEta[2] = {dPtMin,dEtaMin};
  Double_t maxPtEta[2] = {dPtMax,dEtaMax};
  for(Int_t rp=0;rp<2;rp++)
  {
   for(Int_t pe=0;pe<2;pe++)
   {
    for(Int_t sc=0;sc<2;sc++)
    {  
     fHistProNonIsotropicTermsU[rp][pe][sc] = new TProfile(Form("FlowPro_NonIsotropicTerms_%s_%s_%s_SP",rpPoi[rp].Data(),ptEta[pe].Data(),sinCos[sc].Data()),Form("FlowPro_NonIsotropicTerms_%s_%s_%s_SP",rpPoi[rp].Data(),ptEta[pe].Data(),sinCos[sc].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
     fHistList->Add(fHistProNonIsotropicTermsU[rp][pe][sc]);
    } 
   }
  } 
   
  fHistSumOfLinearWeights = new TH1D("Flow_SumOfLinearWeights_SP","Flow_SumOfLinearWeights_SP",1,-0.5, 0.5);
  fHistSumOfLinearWeights -> SetYTitle("sum (*)");
  fHistSumOfLinearWeights -> SetXTitle("sum (Ma*Mb)");
  fHistList->Add(fHistSumOfLinearWeights);
  
  fHistSumOfQuadraticWeights = new TH1D("Flow_SumOfQuadraticWeights_SP","Flow_SumOfQuadraticWeights_SP",1,-0.5, 0.5);
  fHistSumOfQuadraticWeights -> SetYTitle("sum (*)");
  fHistSumOfQuadraticWeights -> SetXTitle("sum (Ma*Mb)^2");
  fHistList->Add(fHistSumOfQuadraticWeights);
  
  fHistProUQQaQbPtRP = new TProfile("FlowPro_UQQaQbPtRP_SP","FlowPro_UQQaQbPtRP_SP",iNbinsPt,dPtMin,dPtMax);
  fHistProUQQaQbPtRP -> SetYTitle("<*>");
  fHistProUQQaQbPtRP -> SetXTitle("<Qu QaQb>");
  fHistList->Add(fHistProUQQaQbPtRP);
  
  fHistProUQQaQbEtaRP = new TProfile("FlowPro_UQQaQbEtaRP_SP","FlowPro_UQQaQbEtaRP_SP",iNbinsEta,dEtaMin,dEtaMax);
  fHistProUQQaQbEtaRP -> SetYTitle("<*>");
  fHistProUQQaQbEtaRP -> SetXTitle("<Qu QaQb>");
  fHistList->Add(fHistProUQQaQbEtaRP);
  
  fHistProUQQaQbPtPOI = new TProfile("FlowPro_UQQaQbPtPOI_SP","FlowPro_UQQaQbPtPOI_SP",iNbinsPt,dPtMin,dPtMax);
  fHistProUQQaQbPtPOI -> SetYTitle("<*>");
  fHistProUQQaQbPtPOI -> SetXTitle("<Qu QaQb>");
  fHistList->Add(fHistProUQQaQbPtPOI);
  
  fHistProUQQaQbEtaPOI = new TProfile("FlowPro_UQQaQbEtaPOI_SP","FlowPro_UQQaQbEtaPOI_SP",iNbinsEta,dEtaMin,dEtaMax);
  fHistProUQQaQbEtaPOI -> SetYTitle("<*>");
  fHistProUQQaQbEtaPOI -> SetXTitle("<Qu QaQb>");
  fHistList->Add(fHistProUQQaQbEtaPOI);
   
  TString weightFlag[3] = {"w_Qu_","w_Qu^2_","w_QuQaQb_"}; 
  for(Int_t i=0;i<3;i++)
  {
   fHistSumOfWeightsPtRP[i] = new TH1D(Form("Flow_SumOfWeights%sPtRP_SP",weightFlag[i].Data()),
                              Form("Flow_SumOfWeights%sPtRP_SP",weightFlag[i].Data()),iNbinsPt,dPtMin,dPtMax);
   fHistSumOfWeightsPtRP[i] -> SetYTitle("sum (*)");
   fHistSumOfWeightsPtRP[i] -> SetXTitle("p_{T}");
   fHistList->Add(fHistSumOfWeightsPtRP[i]);
 
   fHistSumOfWeightsEtaRP[i] = new TH1D(Form("Flow_SumOfWeights%sEtaRP_SP",weightFlag[i].Data()),
                               Form("Flow_SumOfWeights%sEtaRP_SP",weightFlag[i].Data()),iNbinsEta,dEtaMin,dEtaMax);
   fHistSumOfWeightsEtaRP[i] -> SetYTitle("sum (*)");
   fHistSumOfWeightsEtaRP[i] -> SetXTitle("#eta");
   fHistList->Add(fHistSumOfWeightsEtaRP[i]);
  
   fHistSumOfWeightsPtPOI[i] = new TH1D(Form("Flow_SumOfWeights%sPtPOI_SP",weightFlag[i].Data()),
                               Form("Flow_SumOfWeights%sPtPOI_SP",weightFlag[i].Data()),iNbinsPt,dPtMin,dPtMax);
   fHistSumOfWeightsPtPOI[i] -> SetYTitle("sum (*)");
   fHistSumOfWeightsPtPOI[i] -> SetXTitle("p_{T}");
   fHistList->Add(fHistSumOfWeightsPtPOI[i]);
 
   fHistSumOfWeightsEtaPOI[i] = new TH1D(Form("Flow_SumOfWeights%sEtaPOI_SP",weightFlag[i].Data()),
                                Form("Flow_SumOfWeights%sEtaPOI_SP",weightFlag[i].Data()),iNbinsEta,dEtaMin,dEtaMax);
   fHistSumOfWeightsEtaPOI[i] -> SetYTitle("sum (*)");
   fHistSumOfWeightsEtaPOI[i] -> SetXTitle("#eta");
   fHistList->Add(fHistSumOfWeightsEtaPOI[i]);
  }
      
  fCommonHistsSP = new AliFlowCommonHist("AliFlowCommonHistSP");
  fHistList->Add(fCommonHistsSP);
  fCommonHistsResSP = new AliFlowCommonHistResults("AliFlowCommonHistResultsSP","",fHarmonic);
  fHistList->Add(fCommonHistsResSP);  
  fCommonHistsmuQ = new AliFlowCommonHist("AliFlowCommonHistmuQ");
  fHistList->Add(fCommonHistsmuQ);

  (fCommonHistsSP->GetHarmonic())->Fill(0.5,fHarmonic); // store harmonic 
  (fCommonHistsmuQ->GetHarmonic())->Fill(0.5,fHarmonic); // store harmonic 

  fHistQNorm = new TH1D("Flow_QNorm_SP","Flow_QNorm_SP",110,0.,1.1);
  fHistQNorm -> SetYTitle("dN/d(|(Qa+Qb)/(Ma+Mb)|)");
  fHistQNorm -> SetXTitle("|(Qa+Qb)/(Ma+Mb)|");
  fHistList->Add(fHistQNorm);

  fHistQaQb = new TH1D("Flow_QaQb_SP","Flow_QaQb_SP",200,-100.,100.);
  fHistQaQb -> SetYTitle("dN/dQaQb");
  fHistQaQb -> SetXTitle("QaQb");
  fHistList->Add(fHistQaQb);

  fHistQaQbNorm = new TH1D("Flow_QaQbNorm_SP","Flow_QaQbNorm_SP",44,-1.1,1.1);
  fHistQaQbNorm -> SetYTitle("dN/d(QaQb/MaMb)");
  fHistQaQbNorm -> SetXTitle("QaQb/MaMb");
  fHistList->Add(fHistQaQbNorm);

  fHistQNormvsQaQbNorm = new TH2D("Flow_QNormvsQaQbNorm_SP","Flow_QNormvsQaQbNorm_SP",88,-1.1,1.1,22,0.,1.1);
  fHistQNormvsQaQbNorm -> SetYTitle("|Q/Mq|");
  fHistQNormvsQaQbNorm -> SetXTitle("QaQb/MaMb");
  fHistList->Add(fHistQNormvsQaQbNorm);

  fHistQaQbCos = new TH1D("Flow_QaQbCos_SP","Flow_QaQbCos_SP",63,0.,TMath::Pi());
  fHistQaQbCos -> SetYTitle("dN/d(#phi)");
  fHistQaQbCos -> SetXTitle("#phi");
  fHistList->Add(fHistQaQbCos);

  fHistResolution = new TH1D("Flow_resolution_SP","Flow_resolution_SP",100,-1.0,1.0);
  fHistResolution -> SetYTitle("dN/d(cos(2(#phi_a - #phi_b))");
  fHistResolution -> SetXTitle("cos(2*(#phi_a - #phi_b))");
  fHistList->Add(fHistResolution);

  fHistQaNorm = new TH1D("Flow_QaNorm_SP","Flow_QaNorm_SP",22,0.,1.1);
  fHistQaNorm -> SetYTitle("dN/d(|Qa/Ma|)");
  fHistQaNorm -> SetXTitle("|Qa/Ma|");
  fHistList->Add(fHistQaNorm);

  fHistQaNormvsMa = new TH2D("Flow_QaNormvsMa_SP","Flow_QaNormvsMa_SP",100,0.,100.,22,0.,1.1);
  fHistQaNormvsMa -> SetYTitle("|Qa/Ma|");
  fHistQaNormvsMa -> SetXTitle("Ma");
  fHistList->Add(fHistQaNormvsMa);

  fHistQbNorm = new TH1D("Flow_QbNorm_SP","Flow_QbNorm_SP",22,0.,1.1);
  fHistQbNorm -> SetYTitle("dN/d(|Qb/Mb|)");
  fHistQbNorm -> SetXTitle("|Qb/Mb|");
  fHistList->Add(fHistQbNorm);

  fHistQbNormvsMb = new TH2D("Flow_QbNormvsMb_SP","Flow_QbNormvsMb_SP",100,0.,100.,22,0.,1.1);
  fHistQbNormvsMb -> SetYTitle("|Qb/Mb|");
  fHistQbNormvsMb -> SetXTitle("|Mb|");
  fHistList->Add(fHistQbNormvsMb);

  fHistMavsMb = new TH2D("Flow_MavsMb_SP","Flow_MavsMb_SP",100,0.,100.,100,0.,100.);
  fHistMavsMb -> SetYTitle("Ma");
  fHistMavsMb -> SetXTitle("Mb");
  fHistList->Add(fHistMavsMb);


  //weights
  if(fUsePhiWeights) {
    if(!fWeightsList) {
      cout<<"WARNING: fWeightsList is NULL in the Scalar Product method."<<endl;
      exit(0);  
    }
    if(fWeightsList->FindObject("phi_weights_sub0"))  {
      fPhiWeightsSub0 = dynamic_cast<TH1F*>
	(fWeightsList->FindObject("phi_weights_sub0"));
      fHistList->Add(fPhiWeightsSub0);
    } else {
      cout<<"WARNING: histogram with phi weights is not accessible in Scalar Product"<<endl;
      exit(0);
    }
    if(fWeightsList->FindObject("phi_weights_sub1"))  {
      fPhiWeightsSub1 = dynamic_cast<TH1F*>
	(fWeightsList->FindObject("phi_weights_sub1"));
      fHistList->Add(fPhiWeightsSub1);
    } else {
      cout<<"WARNING: histogram with phi weights is not accessible in Scalar Product"<<endl;
      exit(0);
    }

  } // end of if(fUsePhiWeights)

  fEventNumber = 0;  //set number of events to zero 
  
  //store all boolean flags needed in Finish():
  this->StoreFlags();   

  TH1::AddDirectory(oldHistAddStatus);
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::Make(AliFlowEventSimple* anEvent) {


  if (anEvent) {

  //Calculate muQ (for comparing pp and PbPb)
  FillmuQ(anEvent);

  //Calculate flow based on  <QaQb/MaMb> = <v^2>
  FillSP(anEvent);

  }
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::FillSP(AliFlowEventSimple* anEvent) {

  //Calculate flow based on  <QaQb/MaMb> = <v^2>

  //Fill histograms
  if (anEvent) {

    //get Q vectors for the eta-subevents
    AliFlowVector* vQarray = new AliFlowVector[2];
    if (fUsePhiWeights) {
      anEvent->Get2Qsub(vQarray,fHarmonic,fWeightsList,kTRUE);
    } else {
      anEvent->Get2Qsub(vQarray,fHarmonic);
    }
    //subevent a
    AliFlowVector vQa = vQarray[0];
    //subevent b
    AliFlowVector vQb = vQarray[1];

    //For calculating v2 only events should be taken where both subevents are not empty
    //check that the subevents are not empty:
    Double_t dMa = vQa.GetMult();
    Double_t dMb = vQb.GetMult();
    if (dMa > 0. && dMb > 0.) {
      
      //request that the subevent multiplicities are not too different
      //fRelDiffMsub can be set from the configuration macro
      Double_t dRelDiff = TMath::Abs((dMa - dMb)/(dMa + dMb));
      if (dRelDiff < fRelDiffMsub) {

	//fill control histograms 	   
	if (fUsePhiWeights) {
	  fCommonHistsSP->FillControlHistograms(anEvent,fWeightsList,kTRUE);
	} else {
	  fCommonHistsSP->FillControlHistograms(anEvent);
	}

	//fill some SP control histograms
	fHistProQaQb -> Fill(1.,vQa*vQb,1.); //Fill with weight 1 -> Weight with MaMb????
	fHistProQaQbVsM->Fill(anEvent->GetNumberOfRPs()+0.5,(vQa*vQb)/(dMa*dMb),dMa*dMb);  
	fHistQaQbCos ->Fill(TMath::ACos((vQa/vQa.Mod())*(vQb/vQb.Mod())));  //Acos(Qa*Qb) = angle
	fHistResolution -> Fill(TMath::Cos( vQa.Phi()- vQb.Phi() ));  //vQa.Phi() returns 2*phi
	fHistQaQb -> Fill(vQa*vQb);
	fHistMavsMb -> Fill(dMb,dMa);

	//get total Q vector = the sum of subevent a and subevent b
	AliFlowVector vQ;
	Double_t dMQ = 0.; // multiplicity in Q-vector 
	if(!strcmp(fTotalQvector->Data(),"QaQb"))
	{
	 vQ = vQa + vQb;
	 dMQ = dMa + dMb;
	} else if(!strcmp(fTotalQvector->Data(),"Qa"))
	  {
	   vQ = vQa; 
   	   dMQ = dMa;
	  } else if(!strcmp(fTotalQvector->Data(),"Qb"))
	    {
	     vQ = vQb; 
     	     dMQ = dMb;
	    }

	//needed to correct for non-uniform acceptance:
	fHistProNonIsotropicTermsQ->Fill(1.,vQ.Y()/dMQ,dMQ);
	fHistProNonIsotropicTermsQ->Fill(2.,vQ.X()/dMQ,dMQ);

	//weight the Q vectors for the subevents by the multiplicity
	//Note: Weight Q only in the particle loop when it is clear 
	//if it should be (m-1) or M
	Double_t dQXa = vQa.X()/dMa; 
	Double_t dQYa = vQa.Y()/dMa;
	vQa.Set(dQXa,dQYa);
	
	Double_t dQXb = vQb.X()/dMb; 
	Double_t dQYb = vQb.Y()/dMb;
	vQb.Set(dQXb,dQYb);
        
	//scalar product of the two subevents
	Double_t dQaQb = (vQa*vQb);
	fHistProQaQbNorm -> Fill(1.,dQaQb,dMa*dMb);  //Fill (QaQb/MaMb) with weight (MaMb). 
	//needed for the error calculation:
	fHistSumOfLinearWeights -> Fill(0.,dMa*dMb);
	fHistSumOfQuadraticWeights -> Fill(0.,pow(dMa*dMb,2.));
	//needed for correcting non-uniform acceptance: 
	fHistProQaQbReImNorm->Fill(1.,dQYa,dMa); // to get <<sin(phi_a)>>
	fHistProQaQbReImNorm->Fill(2.,dQXa,dMa); // to get <<cos(phi_a)>>
	fHistProQaQbReImNorm->Fill(3.,dQYb,dMb); // to get <<sin(phi_b)>>
	fHistProQaQbReImNorm->Fill(4.,dQXb,dMb); // to get <<cos(phi_b)>>
	
	//fill some SP control histograms
	fHistQaQbNorm ->Fill(vQa*vQb);
	fHistQaNorm ->Fill(vQa.Mod());
	fHistQaNormvsMa->Fill(dMa,vQa.Mod());
	fHistQbNorm ->Fill(vQb.Mod());
	fHistQbNormvsMb->Fill(dMb,vQb.Mod());
	
	//loop over the tracks of the event
	AliFlowTrackSimple*   pTrack = NULL; 
	Int_t iNumberOfTracks = anEvent->NumberOfTracks(); 
	Double_t dMq =  vQ.GetMult();
	
	for (Int_t i=0;i<iNumberOfTracks;i++) 
	  {
	    pTrack = anEvent->GetTrack(i) ; 
	    if (pTrack){
	      Double_t dPhi = pTrack->Phi();
	      Double_t dWeightUQ = 1.; // weight for u*Q 	    
	      //calculate vU
	      TVector2 vU;
	      //do not need to use weight for v as the length will be made 1
	      Double_t dUX = TMath::Cos(fHarmonic*dPhi);
	      Double_t dUY = TMath::Sin(fHarmonic*dPhi);
	      vU.Set(dUX,dUY);
	      Double_t dModulus = vU.Mod();
	      if (dModulus > 0.) vU.Set(dUX/dModulus,dUY/dModulus);  // make length 1
	      else cerr<<"dModulus is zero!"<<endl;
	    
	      //redefine the Q vector and devide by its multiplicity
	      TVector2 vQm;
	      Double_t dQmX = 0.;
	      Double_t dQmY = 0.;
	      //subtract particle from the flowvector if used to define it
	      if (pTrack->InSubevent(0) || pTrack->InSubevent(1)) { 
		//set default phi weight to 1
		Double_t dW = 1.; 
		//if phi weights are used
		if(fUsePhiWeights && fPhiWeightsSub0 && fPhiWeightsSub1) 
		{
		  if(strcmp(fTotalQvector->Data(),"QaQb"))
		  {
		   printf("\n WARNING (SP): If you use phi-weights total Q-vector has to be Qa+Qb in the current implementation!!!! \n");
		   exit(0);
		  }
		  //value of the center of the phi bin
		  Double_t dPhiCenter = 0.;  
		  if (pTrack->InSubevent(0) ) {
		    Int_t iNBinsPhiSub0 = fPhiWeightsSub0->GetNbinsX();
		    Int_t phiBin = 1+(Int_t)(TMath::Floor(dPhi*iNBinsPhiSub0/TMath::TwoPi()));
		    dW = fPhiWeightsSub0->GetBinContent(phiBin); 
		    dPhiCenter = fPhiWeightsSub0->GetBinCenter(phiBin);
		    dQmX = (vQ.X() - dW*(pTrack->Weight())* TMath::Cos(fHarmonic*dPhiCenter) )/(dMq-dW*pTrack->Weight());
		    dQmY = (vQ.Y() - dW*(pTrack->Weight())* TMath::Sin(fHarmonic*dPhiCenter) )/(dMq-dW*pTrack->Weight());
		    
		    vQm.Set(dQmX,dQmY);
		  }

		  else if ( pTrack->InSubevent(1)) { 
		    Int_t iNBinsPhiSub1 = fPhiWeightsSub1->GetNbinsX();
		    Int_t phiBin = 1+(Int_t)(TMath::Floor(dPhi*iNBinsPhiSub1/TMath::TwoPi()));
		    dW = fPhiWeightsSub1->GetBinContent(phiBin);
		    dPhiCenter = fPhiWeightsSub1->GetBinCenter(phiBin);
		    dQmX = (vQ.X() - dW*(pTrack->Weight())* TMath::Cos(fHarmonic*dPhiCenter) )/(dMq-dW*pTrack->Weight());
		    dQmY = (vQ.Y() - dW*(pTrack->Weight())* TMath::Sin(fHarmonic*dPhiCenter) )/(dMq-dW*pTrack->Weight());
		    
		    vQm.Set(dQmX,dQmY);
		  }
		  //bin = 1 + value*nbins/range
		  //TMath::Floor rounds to the lower integer
		}     
		// if no phi weights are used
		else 
		{
		 if(!strcmp(fTotalQvector->Data(),"QaQb"))
		 {
		  dQmX = (vQ.X() - (pTrack->Weight())*dUX)/(dMq-pTrack->Weight());
		  dQmY = (vQ.Y() - (pTrack->Weight())*dUY)/(dMq-pTrack->Weight());
		  dWeightUQ = dMq-pTrack->Weight();
		  vQm.Set(dQmX,dQmY);
		 } else if((!strcmp(fTotalQvector->Data(),"Qa") && pTrack->InSubevent(0)) ||
		           (!strcmp(fTotalQvector->Data(),"Qb") && pTrack->InSubevent(1)))
		   {
		    dQmX = (vQ.X() - (pTrack->Weight())*dUX)/(dMq-pTrack->Weight());
		    dQmY = (vQ.Y() - (pTrack->Weight())*dUY)/(dMq-pTrack->Weight());
		    dWeightUQ = dMq-pTrack->Weight();
		    vQm.Set(dQmX,dQmY);
		   } else if((!strcmp(fTotalQvector->Data(),"Qa") && pTrack->InSubevent(1)) ||
		             (!strcmp(fTotalQvector->Data(),"Qb") && pTrack->InSubevent(0)))
		     {
   		      dQmX = vQ.X()/dMq;
		      dQmY = vQ.Y()/dMq;
		      dWeightUQ = dMq;
		      vQm.Set(dQmX,dQmY);
		     }
		}
			      
		//dUQ = scalar product of vU and vQm
		Double_t dUQ = (vU * vQm);
		Double_t dPt = pTrack->Pt();
		Double_t dEta = pTrack->Eta();
		
		//fill the profile histograms
		if (pTrack->InRPSelection()) {
		  fHistProUQetaRP -> Fill(dEta,dUQ,dWeightUQ); //Fill (Qu/(Mq-1)) with weight (Mq-1) 
		  //needed for the error calculation:
		  fHistProUQQaQbEtaRP -> Fill(dEta,dUQ*dQaQb,dWeightUQ*dMa*dMb); //Fill [Qu/(Mq-1)]*[QaQb/MaMb] with weight (Mq-1)MaMb	    
		  fHistProUQPtRP -> Fill(dPt,dUQ,dWeightUQ);                     //Fill (Qu/(Mq-1)) with weight (Mq-1)
		  fHistProUQQaQbPtRP -> Fill(dPt,dUQ*dQaQb,dWeightUQ*dMa*dMb);   //Fill [Qu/(Mq-1)]*[QaQb/MaMb] with weight (Mq-1)MaMb	
		  
		  fHistSumOfWeightsEtaRP[0]->Fill(dEta,dWeightUQ);        // sum of Mq-1     
		  fHistSumOfWeightsEtaRP[1]->Fill(dEta,pow(dWeightUQ,2.));// sum of (Mq-1)^2     
		  fHistSumOfWeightsEtaRP[2]->Fill(dEta,dWeightUQ*dMa*dMb);// sum of (Mq-1)*MaMb     
		  fHistSumOfWeightsPtRP[0]->Fill(dPt,dWeightUQ);          // sum of Mq-1     
		  fHistSumOfWeightsPtRP[1]->Fill(dPt,pow(dWeightUQ,2.));  // sum of (Mq-1)^2     
		  fHistSumOfWeightsPtRP[2]->Fill(dPt,dWeightUQ*dMa*dMb);  // sum of (Mq-1)*MaMb   
		  //nonisotropic terms:
		  fHistProNonIsotropicTermsU[0][0][0]->Fill(dPt,dUY,1.);
		  fHistProNonIsotropicTermsU[0][0][1]->Fill(dPt,dUX,1.);
		  fHistProNonIsotropicTermsU[0][1][0]->Fill(dEta,dUY,1.);
		  fHistProNonIsotropicTermsU[0][1][1]->Fill(dEta,dUX,1.);
		}
		if (pTrack->InPOISelection()) {
		  fHistProUQetaPOI -> Fill(dEta,dUQ,dWeightUQ);//Fill (Qu/(Mq-1)) with weight (Mq-1)
		  //needed for the error calculation:
		  fHistProUQQaQbEtaPOI -> Fill(dEta,dUQ*dQaQb,dWeightUQ*dMa*dMb); //Fill [Qu/(Mq-1)]*[QaQb/MaMb] with weight (Mq-1)MaMb	    
		  fHistProUQPtPOI -> Fill(dPt,dUQ,dWeightUQ);                     //Fill (Qu/(Mq-1)) with weight (Mq-1)
		  fHistProUQQaQbPtPOI -> Fill(dPt,dUQ*dQaQb,dWeightUQ*dMa*dMb);   //Fill [Qu/(Mq-1)]*[QaQb/MaMb] with weight (Mq-1)MaMb	    
		  
		  fHistSumOfWeightsEtaPOI[0]->Fill(dEta,dWeightUQ);        // sum of Mq-1     
		  fHistSumOfWeightsEtaPOI[1]->Fill(dEta,pow(dWeightUQ,2.));// sum of (Mq-1)^2     
		  fHistSumOfWeightsEtaPOI[2]->Fill(dEta,dWeightUQ*dMa*dMb);// sum of (Mq-1)*MaMb     
		  fHistSumOfWeightsPtPOI[0]->Fill(dPt,dWeightUQ);          // sum of Mq-1     
		  fHistSumOfWeightsPtPOI[1]->Fill(dPt,pow(dWeightUQ,2.)); // sum of (Mq-1)^2     
		  fHistSumOfWeightsPtPOI[2]->Fill(dPt,dWeightUQ*dMa*dMb); // sum of (Mq-1)*MaMb   
		  //nonisotropic terms:
		  fHistProNonIsotropicTermsU[1][0][0]->Fill(dPt,dUY,1.);
		  fHistProNonIsotropicTermsU[1][0][1]->Fill(dPt,dUX,1.);
		  fHistProNonIsotropicTermsU[1][1][0]->Fill(dEta,dUY,1.);
		  fHistProNonIsotropicTermsU[1][1][1]->Fill(dEta,dUX,1.);		  	     
		}  
		
	      } else { //do not subtract the particle from the flowvector
		dQmX = vQ.X()/dMq;
		dQmY = vQ.Y()/dMq;
		vQm.Set(dQmX,dQmY);

		//fill histograms with vQm
		fHistProQNorm->Fill(1.,vQm.Mod(),dMq);
		fHistQNorm->Fill(vQm.Mod());
		fHistQNormvsQaQbNorm->Fill(vQa*vQb ,vQm.Mod()); 
	      
		//dUQ = scalar product of vU and vQm
		Double_t dUQ = (vU * vQm);
		Double_t dPt = pTrack->Pt();
		Double_t dEta = pTrack->Eta();
		
		//fill the profile histograms
		if (pTrack->InRPSelection()) {
		  fHistProUQetaRP -> Fill(dEta,dUQ,dMq);                   //Fill (Qu/Mq) with weight Mq 
		  //needed for the error calculation:
		  fHistProUQQaQbEtaRP -> Fill(dEta,dUQ*dQaQb,dMq*dMa*dMb); //Fill [Qu/Mq]*[QaQb/MaMb] with weight Mq*MaMb	    
		  fHistProUQPtRP -> Fill(dPt,dUQ,dMq);                     //Fill (Qu/Mq) with weight Mq 
		  fHistProUQQaQbPtRP -> Fill(dPt,dUQ*dQaQb,dMq*dMa*dMb);   //Fill [Qu/Mq]*[QaQb/MaMb] with weight Mq*MaMb	    
		  
		  fHistSumOfWeightsEtaRP[0]->Fill(dEta,dMq);        // sum of Mq     
		  fHistSumOfWeightsEtaRP[1]->Fill(dEta,pow(dMq,2.));// sum of Mq^2     
		  fHistSumOfWeightsEtaRP[2]->Fill(dEta,dMq*dMa*dMb);// sum of Mq*MaMb     
		  fHistSumOfWeightsPtRP[0]->Fill(dPt,dMq);          // sum of Mq     
		  fHistSumOfWeightsPtRP[1]->Fill(dPt,pow(dMq,2.));  // sum of Mq^2     
		  fHistSumOfWeightsPtRP[2]->Fill(dPt,dMq*dMa*dMb);  // sum of Mq*MaMb   
		  //nonisotropic terms:
		  fHistProNonIsotropicTermsU[0][0][0]->Fill(dPt,dUY,1.);
		  fHistProNonIsotropicTermsU[0][0][1]->Fill(dPt,dUX,1.);
		  fHistProNonIsotropicTermsU[0][1][0]->Fill(dEta,dUY,1.);
		  fHistProNonIsotropicTermsU[0][1][1]->Fill(dEta,dUX,1.);  
		}
		if (pTrack->InPOISelection()) {
		  fHistProUQetaPOI -> Fill(dEta,dUQ,dMq); //Fill (Qu/Mq) with weight Mq 
		  //needed for the error calculation:
		  fHistProUQQaQbEtaPOI -> Fill(dEta,dUQ*dQaQb,dMq*dMa*dMb); //Fill [Qu/Mq]*[QaQb/MaMb] with weight Mq*MaMb	    
		  fHistProUQPtPOI -> Fill(dPt,dUQ,dMq);                     //Fill (Qu/Mq) with weight Mq 
		  fHistProUQQaQbPtPOI -> Fill(dPt,dUQ*dQaQb,dMq*dMa*dMb);   //Fill [Qu/Mq]*[QaQb/MaMb] with weight Mq*MaMb	    
		  
		  fHistSumOfWeightsEtaPOI[0]->Fill(dEta,dMq);        // sum of Mq     
		  fHistSumOfWeightsEtaPOI[1]->Fill(dEta,pow(dMq,2.));// sum of Mq^2     
		  fHistSumOfWeightsEtaPOI[2]->Fill(dEta,dMq*dMa*dMb);// sum of Mq*MaMb     
		  fHistSumOfWeightsPtPOI[0]->Fill(dPt,dMq);          // sum of Mq     
		  fHistSumOfWeightsPtPOI[1]->Fill(dPt,pow(dMq,2.));  // sum of Mq^2     
		  fHistSumOfWeightsPtPOI[2]->Fill(dPt,dMq*dMa*dMb);  // sum of Mq*MaMb     
		  //nonisotropic terms:
		  fHistProNonIsotropicTermsU[1][0][0]->Fill(dPt,dUY,1.);
		  fHistProNonIsotropicTermsU[1][0][1]->Fill(dPt,dUX,1.);
		  fHistProNonIsotropicTermsU[1][1][0]->Fill(dEta,dUY,1.);
		  fHistProNonIsotropicTermsU[1][1][1]->Fill(dEta,dUX,1.);	
		}  
	      }//track not in subevents
	      
	    }//track
	    
	  }//loop over tracks
	
	fEventNumber++;

      } //difference Ma and Mb

    }// subevents not empty 
    delete [] vQarray;

  } //event

}//end of FillSP()

//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::FillmuQ(AliFlowEventSimple* anEvent) {

  if (anEvent) {

    //get Q vectors for the eta-subevents
    AliFlowVector* vQarray = new AliFlowVector[2];
    if (fUsePhiWeights) {
      anEvent->Get2Qsub(vQarray,fHarmonic,fWeightsList,kTRUE);
    } else {
      anEvent->Get2Qsub(vQarray,fHarmonic);
    }
    //subevent a
    AliFlowVector vQa = vQarray[0];
    //subevent b
    AliFlowVector vQb = vQarray[1];

    //get total Q vector = the sum of subevent a and subevent b
    AliFlowVector vQ;
    if(!strcmp(fTotalQvector->Data(),"QaQb"))
    {
     if(vQa.GetMult() > 0 || vQb.GetMult() > 0) 
     {
      vQ = vQa + vQb;
     } else {return;}         
    } else if(!strcmp(fTotalQvector->Data(),"Qa"))
      {
       if(vQa.GetMult() > 0)
       {
        vQ = vQa;
       } else {return;}
      } else if(!strcmp(fTotalQvector->Data(),"Qb"))
        {
         if(vQb.GetMult() > 0)
         {
          vQ = vQb;
         } else {return;}
        }
      
    //For calculating uQ for comparison all events should be taken also if one of the subevents is empty
    //check if the total Q vector is not empty
    Double_t dMq =  vQ.GetMult();
    if (dMq > 0.) {
                  
      //Fill control histograms
      if (fUsePhiWeights) {
	fCommonHistsmuQ->FillControlHistograms(anEvent,fWeightsList,kTRUE);
      } else {
	fCommonHistsmuQ->FillControlHistograms(anEvent);
      }

      //loop over all POI tracks and fill uQ
      AliFlowTrackSimple*   pTrack = NULL; 
      for (Int_t i=0;i<anEvent->NumberOfTracks();i++) {
	pTrack = anEvent->GetTrack(i) ; 
	if (pTrack){

	  if (pTrack->InPOISelection()) {

	    Double_t dPhi = pTrack->Phi();
	    //weights do not need to be used as the length of vU will be set to 1
	    	    
	    //calculate vU
	    TVector2 vU;
	    Double_t dUX = TMath::Cos(fHarmonic*dPhi);
	    Double_t dUY = TMath::Sin(fHarmonic*dPhi);
	    vU.Set(dUX,dUY);
	    Double_t dModulus = vU.Mod();
	    // make length 1
	    if (dModulus!=0.) vU.Set(dUX/dModulus,dUY/dModulus);  
	    else cerr<<"dModulus is zero!"<<endl;
	    
	    //redefine the Q vector 
	    TVector2 vQm;
	    Double_t dQmX = 0.;
	    Double_t dQmY = 0.;
	    //subtract particle from the flowvector if used to define it
	    if (pTrack->InSubevent(0) || pTrack->InSubevent(1)) { 
	      //the number of tracks contributing to vQ must be more than 1
	      if (dMq > 1) { 
		//set default phi weight to 1
		Double_t dW = 1.; 
		//if phi weights are used
		if(fUsePhiWeights && fPhiWeightsSub0 && fPhiWeightsSub1) 
		{
		 if(strcmp(fTotalQvector->Data(),"QaQb"))
		 {
	              printf("\n WARNING (SP): If you use phi-weights total Q-vector has to be Qa+Qb in the current implementation!!!! \n");
	              exit(0);
	             }

		  //value of the center of the phi bin
		  Double_t dPhiCenter = 0.;  
		  if (pTrack->InSubevent(0) ) {
		    Int_t iNBinsPhiSub0 = fPhiWeightsSub0->GetNbinsX();
		    Int_t phiBin = 1+(Int_t)(TMath::Floor(dPhi*iNBinsPhiSub0/TMath::TwoPi()));
		    dW = fPhiWeightsSub0->GetBinContent(phiBin); 
		    dPhiCenter = fPhiWeightsSub0->GetBinCenter(phiBin);
		    dQmX = (vQ.X() - dW*(pTrack->Weight())* TMath::Cos(fHarmonic*dPhiCenter) );
		    dQmY = (vQ.Y() - dW*(pTrack->Weight())* TMath::Sin(fHarmonic*dPhiCenter) );
		    
		    vQm.Set(dQmX,dQmY);
		  }
		
		  else if ( pTrack->InSubevent(1)) { 
		    Int_t iNBinsPhiSub1 = fPhiWeightsSub1->GetNbinsX();
		    Int_t phiBin = 1+(Int_t)(TMath::Floor(dPhi*iNBinsPhiSub1/TMath::TwoPi()));
		    dW = fPhiWeightsSub1->GetBinContent(phiBin);
		    dPhiCenter = fPhiWeightsSub1->GetBinCenter(phiBin);
		    dQmX = (vQ.X() - dW*(pTrack->Weight())* TMath::Cos(fHarmonic*dPhiCenter) );
		    dQmY = (vQ.Y() - dW*(pTrack->Weight())* TMath::Sin(fHarmonic*dPhiCenter) );
		    
		    vQm.Set(dQmX,dQmY);
		  }
		  //bin = 1 + value*nbins/range
		  //TMath::Floor rounds to the lower integer
		}     
		// if no phi weights are used
		else 
		{
		 if(!strcmp(fTotalQvector->Data(),"QaQb"))
		 {
		  dQmX = (vQ.X() - (pTrack->Weight())*dUX);
		  dQmY = (vQ.Y() - (pTrack->Weight())*dUY);
		  vQm.Set(dQmX,dQmY);
		 } else if((!strcmp(fTotalQvector->Data(),"Qa") && pTrack->InSubevent(0)) ||
		           (!strcmp(fTotalQvector->Data(),"Qb") && pTrack->InSubevent(1)))
		   {
		    //printf("\n A \n");exit(0);
		    dQmX = (vQ.X() - (pTrack->Weight())*dUX);
		    dQmY = (vQ.Y() - (pTrack->Weight())*dUY);
		    vQm.Set(dQmX,dQmY);
		   } else if((!strcmp(fTotalQvector->Data(),"Qa") && pTrack->InSubevent(1)) ||
		             (!strcmp(fTotalQvector->Data(),"Qb") && pTrack->InSubevent(0)))
		     {
   		      //printf("\n B \n");exit(0);
		      dQmX = vQ.X();
		      dQmY = vQ.Y();
		      vQm.Set(dQmX,dQmY);
		     }
		}

		//dUQ = scalar product of vU and vQm
		Double_t dUQ = (vU * vQm);
		Double_t dPt = pTrack->Pt();
		Double_t dEta = pTrack->Eta();
		//fill the profile histograms
		fHistProUQetaAllEventsPOI -> Fill(dEta,dUQ);   //Fill (Qu)
		fHistProUQPtAllEventsPOI -> Fill(dPt,dUQ);     //Fill (Qu)
	      
	      } //dMq > 1
	    } 
	    else { //do not subtract the particle from the flowvector

	      dQmX = vQ.X();
	      dQmY = vQ.Y();
	      vQm.Set(dQmX,dQmY);
	   
	      //dUQ = scalar product of vU and vQm
	      Double_t dUQ = (vU * vQm);
	      Double_t dPt = pTrack->Pt();
	      Double_t dEta = pTrack->Eta();
	      //fill the profile histograms
	      fHistProUQetaAllEventsPOI -> Fill(dEta,dUQ);   //Fill (Qu)
	      fHistProUQPtAllEventsPOI -> Fill(dPt,dUQ);     //Fill (Qu)
	       
	    }

	  } //in POI selection
	} //track valid
      } //end of loop over tracks
    } //Q vector is not empty
           
  } //anEvent valid
  
} //end of FillmuQ

//--------------------------------------------------------------------  
void AliFlowAnalysisWithScalarProduct::GetOutputHistograms(TList *outputListHistos){
  
  //get pointers to all output histograms (called before Finish())

  if (outputListHistos) {
  //Get the common histograms from the output list
    AliFlowCommonHist *pCommonHistSP = dynamic_cast<AliFlowCommonHist*> 
      (outputListHistos->FindObject("AliFlowCommonHistSP"));
    AliFlowCommonHistResults *pCommonHistResultsSP = dynamic_cast<AliFlowCommonHistResults*> 
      (outputListHistos->FindObject("AliFlowCommonHistResultsSP"));
    AliFlowCommonHist *pCommonHistmuQ = dynamic_cast<AliFlowCommonHist*> 
      (outputListHistos->FindObject("AliFlowCommonHistmuQ"));

    TProfile* pHistProQNorm    = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_QNorm_SP"));
    TProfile* pHistProQaQb     = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_QaQb_SP"));
    TProfile* pHistProQaQbNorm = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_QaQbNorm_SP"));
    TProfile* pHistProQaQbReImNorm = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_QaQbReImNorm_SP"));
    TProfile* pHistProNonIsotropicTermsQ = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_NonIsotropicTermsQ_SP"));
    TH1D*     pHistSumOfLinearWeights    = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_SumOfLinearWeights_SP"));
    TH1D*     pHistSumOfQuadraticWeights = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_SumOfQuadraticWeights_SP"));

    TProfile* pHistProFlags    = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_Flags_SP"));
    TProfile* pHistProUQetaRP  = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_UQetaRP_SP"));
    TProfile* pHistProUQetaPOI = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_UQetaPOI_SP"));
    TProfile* pHistProUQPtRP   = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_UQPtRP_SP"));
    TProfile* pHistProUQPtPOI  = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_UQPtPOI_SP"));
    TProfile* pHistProUQQaQbPtRP    = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_UQQaQbPtRP_SP"));
    TProfile* pHistProUQQaQbEtaRP   = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_UQQaQbEtaRP_SP"));
    TProfile* pHistProUQQaQbPtPOI   = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_UQQaQbPtPOI_SP"));
    TProfile* pHistProUQQaQbEtaPOI  = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_UQQaQbEtaPOI_SP"));
    TString weightFlag[3] = {"w_Qu_","w_Qu^2_","w_QuQaQb_"}; 

   
    TH1D* pHistSumOfWeightsPtRP[3] = {NULL};                    
    TH1D* pHistSumOfWeightsEtaRP[3] = {NULL};                    
    TH1D* pHistSumOfWeightsPtPOI[3] = {NULL};                    
    TH1D* pHistSumOfWeightsEtaPOI[3] = {NULL}; 
    
    for(Int_t i=0;i<3;i++) {
      pHistSumOfWeightsPtRP[i]   = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("Flow_SumOfWeights%sPtRP_SP",weightFlag[i].Data())));
      pHistSumOfWeightsEtaRP[i]  = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("Flow_SumOfWeights%sEtaRP_SP",weightFlag[i].Data())));
      pHistSumOfWeightsPtPOI[i]  = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("Flow_SumOfWeights%sPtPOI_SP",weightFlag[i].Data())));
      pHistSumOfWeightsEtaPOI[i] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("Flow_SumOfWeights%sEtaPOI_SP",weightFlag[i].Data())));
    }   
    
    TString rpPoi[2] = {"RP","POI"};
    TString ptEta[2] = {"Pt","Eta"};
    TString sinCos[2] = {"sin","cos"};
    TProfile *pHistProNonIsotropicTermsU[2][2][2] = {{{NULL}}};
    for(Int_t rp=0;rp<2;rp++) {
      for(Int_t pe=0;pe<2;pe++)	{
	for(Int_t sc=0;sc<2;sc++) {      
	  pHistProNonIsotropicTermsU[rp][pe][sc] = dynamic_cast<TProfile*>(outputListHistos->FindObject(Form("FlowPro_NonIsotropicTerms_%s_%s_%s_SP",rpPoi[rp].Data(),ptEta[pe].Data(),sinCos[sc].Data())));   
	} 
      }
    }   
 
    TH1D*     pHistQNorm    = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QNorm_SP"));
    TH1D*     pHistQaQb     = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QaQb_SP"));
    TH1D*     pHistQaQbNorm = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QaQbNorm_SP"));
    TH2D*     pHistQNormvsQaQbNorm = dynamic_cast<TH2D*>(outputListHistos->FindObject("Flow_QNormvsQaQbNorm_SP"));
    TH1D*     pHistQaQbCos  = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QaQbCos_SP"));
    TH1D*     pHistResolution = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_resolution_SP"));
    TH1D*     pHistQaNorm   = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QaNorm_SP"));
    TH2D*     pHistQaNormvsMa   = dynamic_cast<TH2D*>(outputListHistos->FindObject("Flow_QaNormvsMa_SP"));
    TH1D*     pHistQbNorm   = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QbNorm_SP"));
    TH2D*     pHistQbNormvsMb   = dynamic_cast<TH2D*>(outputListHistos->FindObject("Flow_QbNormvsMb_SP"));
    TH2D*     pHistMavsMb = dynamic_cast<TH2D*>(outputListHistos->FindObject("Flow_MavsMb_SP"));

    //pass the pointers to the task
    if (pCommonHistSP && 
	pCommonHistResultsSP && 
	pCommonHistmuQ &&
	pHistProQNorm && 
	pHistProQaQb && 
	pHistProQaQbNorm && 
	pHistProQaQbReImNorm && 
	pHistProNonIsotropicTermsQ &&
	pHistSumOfLinearWeights && 
	pHistSumOfQuadraticWeights && 
	pHistProFlags &&
	pHistProUQetaRP	&& 
	pHistProUQetaPOI && 
	pHistProUQPtRP && 
	pHistProUQPtPOI &&  
	pHistProUQQaQbPtRP && 
	pHistProUQQaQbEtaRP && 
	pHistProUQQaQbPtPOI && 
	pHistProUQQaQbEtaPOI &&
	pHistSumOfWeightsPtRP[0] && pHistSumOfWeightsPtRP[1] && pHistSumOfWeightsPtRP[2] &&
	pHistSumOfWeightsEtaRP[0] && pHistSumOfWeightsEtaRP[1] && pHistSumOfWeightsEtaRP[2] &&
	pHistSumOfWeightsPtPOI[0] && pHistSumOfWeightsPtPOI[1] && pHistSumOfWeightsPtPOI[2] &&
	pHistSumOfWeightsEtaPOI[0] && pHistSumOfWeightsEtaPOI[1] && pHistSumOfWeightsEtaPOI[2] && 
	pHistProNonIsotropicTermsU[0][0][0] && pHistProNonIsotropicTermsU[1][0][0] && pHistProNonIsotropicTermsU[0][1][0] && pHistProNonIsotropicTermsU[0][0][1] && 
	pHistProNonIsotropicTermsU[1][1][0] && pHistProNonIsotropicTermsU[1][0][1] && pHistProNonIsotropicTermsU[0][1][1] && pHistProNonIsotropicTermsU[1][1][1] &&
	pHistQNorm && 
	pHistQaQb && 
	pHistQaQbNorm && 
	pHistQNormvsQaQbNorm &&
	pHistQaQbCos && 
	pHistResolution &&
	pHistQaNorm && 
	pHistQaNormvsMa && 
	pHistQbNorm && 
	pHistQbNormvsMb && 
	pHistMavsMb 
	) {

      this -> SetCommonHistsSP(pCommonHistSP);
      this -> SetCommonHistsResSP(pCommonHistResultsSP);
      this -> SetCommonHistsmuQ(pCommonHistmuQ);
      this -> SetHistProQNorm(pHistProQNorm);
      this -> SetHistProQaQb(pHistProQaQb);
      this -> SetHistProQaQbNorm(pHistProQaQbNorm);
      this -> SetHistProQaQbReImNorm(pHistProQaQbReImNorm);      
      this -> SetHistProNonIsotropicTermsQ(pHistProNonIsotropicTermsQ);
      this -> SetHistSumOfLinearWeights(pHistSumOfLinearWeights);
      this -> SetHistSumOfQuadraticWeights(pHistSumOfQuadraticWeights); 
      this -> SetHistProFlags(pHistProFlags);
      this -> SetHistProUQetaRP(pHistProUQetaRP);
      this -> SetHistProUQetaPOI(pHistProUQetaPOI);
      this -> SetHistProUQPtRP(pHistProUQPtRP);
      this -> SetHistProUQPtPOI(pHistProUQPtPOI);
      this -> SetHistProUQQaQbPtRP(pHistProUQQaQbPtRP);
      this -> SetHistProUQQaQbEtaRP(pHistProUQQaQbEtaRP);
      this -> SetHistProUQQaQbPtPOI(pHistProUQQaQbPtPOI);
      this -> SetHistProUQQaQbEtaPOI(pHistProUQQaQbEtaPOI); 
      for(Int_t i=0;i<3;i++) {
	if(pHistSumOfWeightsPtRP[i]) this -> SetHistSumOfWeightsPtRP(pHistSumOfWeightsPtRP[i],i);      
	if(pHistSumOfWeightsEtaRP[i]) this -> SetHistSumOfWeightsEtaRP(pHistSumOfWeightsEtaRP[i],i);      
	if(pHistSumOfWeightsPtPOI[i]) this -> SetHistSumOfWeightsPtPOI(pHistSumOfWeightsPtPOI[i],i);      
	if(pHistSumOfWeightsEtaPOI[i]) this -> SetHistSumOfWeightsEtaPOI(pHistSumOfWeightsEtaPOI[i],i);      
      }
      for(Int_t rp=0;rp<2;rp++)  {
	for(Int_t pe=0;pe<2;pe++) {
	  for(Int_t sc=0;sc<2;sc++) {
	    if(pHistProNonIsotropicTermsU[rp][pe][sc]) {
	      this->SetHistProNonIsotropicTermsU(pHistProNonIsotropicTermsU[rp][pe][sc],rp,pe,sc);
	    }
	  }
	}
      }        
      this -> SetHistQNorm(pHistQNorm);
      this -> SetHistQaQb(pHistQaQb);
      this -> SetHistQaQbNorm(pHistQaQbNorm);
      this -> SetHistQNormvsQaQbNorm(pHistQNormvsQaQbNorm);
      this -> SetHistQaQbCos(pHistQaQbCos);
      this -> SetHistResolution(pHistResolution);
      this -> SetHistQaNorm(pHistQaNorm);
      this -> SetHistQaNormvsMa(pHistQaNormvsMa);
      this -> SetHistQbNorm(pHistQbNorm);
      this -> SetHistQbNormvsMb(pHistQbNormvsMb);
      this -> SetHistMavsMb(pHistMavsMb);

    } else {
      cout<<"WARNING: Histograms needed to run Finish() in SP are not accessable!"<<endl; }
         
  } // end of if(outputListHistos)
}            

//--------------------------------------------------------------------            
void AliFlowAnalysisWithScalarProduct::Finish() {
   
  //calculate flow and fill the AliFlowCommonHistResults
  if (fDebug) cout<<"AliFlowAnalysisWithScalarProduct::Finish()"<<endl;
  
  // access harmonic:
  if(fCommonHistsSP->GetHarmonic())
  {
   fHarmonic = (Int_t)(fCommonHistsSP->GetHarmonic())->GetBinContent(1); 
  }
  
  //access all boolean flags needed in Finish():
  this->AccessFlags();

  cout<<"*************************************"<<endl;
  cout<<"*************************************"<<endl;
  cout<<"      Integrated flow from           "<<endl;
  cout<<"         Scalar product              "<<endl;
  cout<<endl;
  
  Int_t iNbinsPt  = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
  Int_t iNbinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
   
  //Calculate the event plane resolution
  //----------------------------------
  Double_t dCos2phi = fHistResolution->GetMean();
  if (dCos2phi > 0.0){
    Double_t dResolution = TMath::Sqrt(2*dCos2phi); 
    cout<<"An estimate of the event plane resolution is: "<<dResolution<<endl;
    cout<<endl;
  }

  //Calculate reference flow (noname)
  //----------------------------------
  //weighted average over (QaQb/MaMb) with weight (MaMb)
  Double_t dQaQb  = fHistProQaQbNorm->GetBinContent(1);
  Double_t dSpreadQaQb = fHistProQaQbNorm->GetBinError(1);
  Double_t dEntriesQaQb = fHistProQaQbNorm->GetEntries();
  
  //non-isotropic terms:  
  Double_t dImQa = fHistProQaQbReImNorm->GetBinContent(1);
  Double_t dReQa = fHistProQaQbReImNorm->GetBinContent(2);
  Double_t dImQb = fHistProQaQbReImNorm->GetBinContent(3);
  Double_t dReQb = fHistProQaQbReImNorm->GetBinContent(4);

  if(fApplyCorrectionForNUA) 
  {
   dQaQb = dQaQb - dImQa*dImQb - dReQa*dReQb; 
  }
  
  if (dEntriesQaQb > 0.) {
    cout<<"QaQb = "<<dQaQb<<" +- "<<(dSpreadQaQb/TMath::Sqrt(dEntriesQaQb))<<endl;
    cout<<endl;
  }

  Double_t dV = -999.; 
  if(dQaQb>=0.)
  {
   dV = TMath::Sqrt(dQaQb); 
  }
  //statistical error of dQaQb: 
  //  statistical error = term1 * spread * term2:
  //  term1 = sqrt{sum_{i=1}^{N} w^2}/(sum_{i=1}^{N} w)
  //  term2 = 1/sqrt(1-term1^2) 
  Double_t dSumOfLinearWeights = fHistSumOfLinearWeights->GetBinContent(1);
  Double_t dSumOfQuadraticWeights = fHistSumOfQuadraticWeights->GetBinContent(1);
  Double_t dTerm1 = 0.;
  Double_t dTerm2 = 0.;
  if(dSumOfLinearWeights) {
    dTerm1 = pow(dSumOfQuadraticWeights,0.5)/dSumOfLinearWeights;
  }
  if(1.-pow(dTerm1,2.) > 0.) {
    dTerm2 = 1./pow(1-pow(dTerm1,2.),0.5);
  }
  Double_t dStatErrorQaQb = dTerm1 * dSpreadQaQb * dTerm2;
  //calculate the statistical error
  Double_t dVerr = 0.;
  if(dQaQb > 0.) { 
    dVerr = (1./(2.*pow(dQaQb,0.5)))*dStatErrorQaQb;
  } 
  fCommonHistsResSP->FillIntegratedFlow(dV,dVerr);
  cout<<Form("v%i(subevents) = ",fHarmonic)<<dV<<" +- "<<dVerr<<endl;
	
  //Calculate differential flow and integrated flow (RP, POI)
  //---------------------------------------------------------
  //v as a function of eta for RP selection
  for(Int_t b=1;b<iNbinsEta+1;b++) {
    Double_t duQpro = fHistProUQetaRP->GetBinContent(b);
    if(fApplyCorrectionForNUA) {
      duQpro = duQpro 
	- fHistProNonIsotropicTermsU[0][1][1]->GetBinContent(b)*fHistProNonIsotropicTermsQ->GetBinContent(2)
	- fHistProNonIsotropicTermsU[0][1][0]->GetBinContent(b)*fHistProNonIsotropicTermsQ->GetBinContent(1);  
    }
    Double_t dv2pro = -999.;
    if (dV!=0.) { dv2pro = duQpro/dV; }
    //calculate the statistical error
    Double_t dv2ProErr = CalculateStatisticalError(b, dStatErrorQaQb, fHistProUQetaRP, fHistProUQQaQbEtaRP, fHistSumOfWeightsEtaRP);
    //fill TH1D
    fCommonHistsResSP->FillDifferentialFlowEtaRP(b, dv2pro, dv2ProErr);   
  } //loop over bins b


  //v as a function of eta for POI selection
  for(Int_t b=1;b<iNbinsEta+1;b++) {
    Double_t duQpro = fHistProUQetaPOI->GetBinContent(b);
    if(fApplyCorrectionForNUA)  {
      duQpro = duQpro 
	- fHistProNonIsotropicTermsU[1][1][1]->GetBinContent(b)*fHistProNonIsotropicTermsQ->GetBinContent(2)
	- fHistProNonIsotropicTermsU[1][1][0]->GetBinContent(b)*fHistProNonIsotropicTermsQ->GetBinContent(1); 
    }    
    Double_t dv2pro = -999.;
    if (dV!=0.) { dv2pro = duQpro/dV; }
    //calculate the statistical error
    Double_t dv2ProErr = CalculateStatisticalError(b, dStatErrorQaQb, fHistProUQetaPOI, fHistProUQQaQbEtaPOI, fHistSumOfWeightsEtaPOI);
   
    //fill TH1D
    fCommonHistsResSP->FillDifferentialFlowEtaPOI(b, dv2pro, dv2ProErr); 
  } //loop over bins b
  

  //v as a function of Pt for RP selection
  TH1F* fHistPtRP = fCommonHistsSP->GetHistPtRP(); //for calculating integrated flow
  Double_t dVRP = 0.;
  Double_t dSumRP = 0.;
  Double_t dErrVRP =0.;
  
  for(Int_t b=1;b<iNbinsPt+1;b++) {
    Double_t duQpro = fHistProUQPtRP->GetBinContent(b);
    if(fApplyCorrectionForNUA) {
      duQpro = duQpro 
	- fHistProNonIsotropicTermsU[0][0][1]->GetBinContent(b)*fHistProNonIsotropicTermsQ->GetBinContent(2)
	- fHistProNonIsotropicTermsU[0][0][0]->GetBinContent(b)*fHistProNonIsotropicTermsQ->GetBinContent(1);  
    }
    Double_t dv2pro = -999.;
    if (dV!=0.) { dv2pro = duQpro/dV; }
    //calculate the statistical error
    Double_t dv2ProErr = CalculateStatisticalError(b, dStatErrorQaQb, fHistProUQPtRP, fHistProUQQaQbPtRP, fHistSumOfWeightsPtRP);
              
    //fill TH1D
    fCommonHistsResSP->FillDifferentialFlowPtRP(b, dv2pro, dv2ProErr);

    //calculate integrated flow for RP selection
    if (fHistPtRP){
      Double_t dYieldPt = fHistPtRP->GetBinContent(b);
      dVRP += dv2pro*dYieldPt;
      dSumRP +=dYieldPt;
      dErrVRP += dYieldPt*dYieldPt*dv2ProErr*dv2ProErr;
    } else { cout<<"fHistPtRP is NULL"<<endl; }
  } //loop over bins b
  
  if (dSumRP != 0.) {
    dVRP /= dSumRP; //the pt distribution should be normalised
    dErrVRP /= (dSumRP*dSumRP);
    dErrVRP = TMath::Sqrt(dErrVRP);
  }
  fCommonHistsResSP->FillIntegratedFlowRP(dVRP,dErrVRP);
  cout<<Form("v%i(RP) = ",fHarmonic)<<dVRP<<" +- "<<dErrVRP<<endl;
  

  //v as a function of Pt for POI selection 
  TH1F* fHistPtPOI = fCommonHistsSP->GetHistPtPOI(); //for calculating integrated flow
  Double_t dVPOI = 0.;
  Double_t dSumPOI = 0.;
  Double_t dErrVPOI =0.;
  
  for(Int_t b=1;b<iNbinsPt+1;b++) {
    Double_t duQpro = fHistProUQPtPOI->GetBinContent(b);
    if(fApplyCorrectionForNUA)  {
     duQpro = duQpro  
       - fHistProNonIsotropicTermsU[1][0][1]->GetBinContent(b)*fHistProNonIsotropicTermsQ->GetBinContent(2)
       - fHistProNonIsotropicTermsU[1][0][0]->GetBinContent(b)*fHistProNonIsotropicTermsQ->GetBinContent(1);  
    }    
    Double_t dv2pro = -999.;
    if (dV!=0.) { dv2pro = duQpro/dV; }
    //calculate the statistical error
    Double_t dv2ProErr = CalculateStatisticalError(b, dStatErrorQaQb, fHistProUQPtPOI, fHistProUQQaQbPtPOI, fHistSumOfWeightsPtPOI);
        
    //fill TH1D
    fCommonHistsResSP->FillDifferentialFlowPtPOI(b, dv2pro, dv2ProErr); 
    
    //calculate integrated flow for POI selection
    if (fHistPtPOI){
      Double_t dYieldPt = fHistPtPOI->GetBinContent(b);
      dVPOI += dv2pro*dYieldPt;
      dSumPOI +=dYieldPt;
      dErrVPOI += dYieldPt*dYieldPt*dv2ProErr*dv2ProErr;
    } else { cout<<"fHistPtPOI is NULL"<<endl; }
  } //loop over bins b
  
  if (dSumPOI > 0.) {
    dVPOI /= dSumPOI; //the pt distribution should be normalised
    dErrVPOI /= (dSumPOI*dSumPOI);
    dErrVPOI = TMath::Sqrt(dErrVPOI);
  }
  fCommonHistsResSP->FillIntegratedFlowPOI(dVPOI,dErrVPOI);
  cout<<Form("v%i(POI) = ",fHarmonic)<<dVPOI<<" +- "<<dErrVPOI<<endl;

  cout<<endl;
  cout<<"*************************************"<<endl;
  cout<<"*************************************"<<endl;   	  
     
  //cout<<".....finished"<<endl;
}


//--------------------------------------------------------------------            
Double_t AliFlowAnalysisWithScalarProduct::CalculateStatisticalError(Int_t b, Double_t aStatErrorQaQb, TProfile* pHistProUQ, TProfile* pHistProUQQaQb, TH1D** pHistSumOfWeights) {
  //calculate the statistical error for differential flow for bin b
  Double_t duQproSpread = pHistProUQ->GetBinError(b);
  Double_t sumOfMq = pHistSumOfWeights[0]->GetBinContent(b);
  Double_t sumOfMqSquared = pHistSumOfWeights[1]->GetBinContent(b);
  Double_t dQaQb = fHistProQaQbNorm->GetBinContent(1);  
  //non-isotropic terms:  
  Double_t dImQa = fHistProQaQbReImNorm->GetBinContent(1);
  Double_t dReQa = fHistProQaQbReImNorm->GetBinContent(2);
  Double_t dImQb = fHistProQaQbReImNorm->GetBinContent(3);
  Double_t dReQb = fHistProQaQbReImNorm->GetBinContent(4);
  if(fApplyCorrectionForNUA) 
  {
   dQaQb = dQaQb - dImQa*dImQb - dReQa*dReQb; 
  }  
  Double_t dTerm1 = 0.;
  Double_t dTerm2 = 0.;
  if(sumOfMq) {
    dTerm1 = (pow(sumOfMqSquared,0.5)/sumOfMq);
  } 
  if(1.-pow(dTerm1,2.)>0.) {
    dTerm2 = 1./pow(1.-pow(dTerm1,2.),0.5); 
  }
  Double_t duQproErr = dTerm1*duQproSpread*dTerm2;
  // covariances:
  Double_t dTerm1Cov = pHistSumOfWeights[2]->GetBinContent(b);
  Double_t dTerm2Cov = fHistSumOfLinearWeights->GetBinContent(1);
  Double_t dTerm3Cov = sumOfMq;
  Double_t dWeightedCovariance = 0.;
  if(dTerm2Cov*dTerm3Cov>0.) {
    Double_t dDenominator = 1.-dTerm1Cov/(dTerm2Cov*dTerm3Cov);
    Double_t dPrefactor = dTerm1Cov/(dTerm2Cov*dTerm3Cov);
    if(dDenominator!=0) {
      Double_t dCovariance = (pHistProUQQaQb->GetBinContent(b)-dQaQb*pHistProUQ->GetBinContent(b))/dDenominator;            
      dWeightedCovariance = dCovariance*dPrefactor; 
    }
  }
  
  Double_t dv2ProErr = 0.; // final statitical error: 
  if(dQaQb>0.) {
    Double_t dv2ProErrorSquared = (1./4.)*pow(dQaQb,-3.)*
      (pow(pHistProUQ->GetBinContent(b),2.)*pow(aStatErrorQaQb,2.)
       + 4.*pow(dQaQb,2.)*pow(duQproErr,2.)
       - 4.*dQaQb*pHistProUQ->GetBinContent(b)*dWeightedCovariance);
    if(dv2ProErrorSquared>0.) dv2ProErr = pow(dv2ProErrorSquared,0.5);
  } 
   
  return dv2ProErr;
}


//--------------------------------------------------------------------     

void AliFlowAnalysisWithScalarProduct::StoreFlags()
{
 // Store all boolean flags needed in Finish() in profile fHistProFlags.

 // Apply correction for non-uniform acceptance or not:
 fHistProFlags->Fill(0.5,fApplyCorrectionForNUA);

} 

//-------------------------------------------------------------------- 

void AliFlowAnalysisWithScalarProduct::AccessFlags()
{
 // Access all boolean flags needed in Finish() from profile fHistProFlags.

 // Apply correction for non-uniform acceptance or not:
 fApplyCorrectionForNUA = (Bool_t) fHistProFlags->GetBinContent(1);

} 

//--------------------------------------------------------------------     
