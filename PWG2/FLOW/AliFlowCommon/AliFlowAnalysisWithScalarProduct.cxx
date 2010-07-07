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
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowAnalysisWithScalarProduct.h"
#include "AliFlowVector.h"

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
   fRelDiffMsub(1.),
   fWeightsList(NULL),
   fUsePhiWeights(kFALSE),
   fPhiWeightsSub0(NULL),
   fPhiWeightsSub1(NULL),
   fHistList(NULL),
   fHistProUQetaRP(NULL),
   fHistProUQetaPOI(NULL),
   fHistProUQPtRP(NULL),
   fHistProUQPtPOI(NULL),
   fHistProQaQb(NULL),
   fHistProQaQbNorm(NULL),
   fHistSumOfLinearWeights(NULL),
   fHistSumOfQuadraticWeights(NULL),
   fHistProUQQaQbPtRP(NULL),
   fHistProUQQaQbEtaRP(NULL),
   fHistProUQQaQbPtPOI(NULL),
   fHistProUQQaQbEtaPOI(NULL),
   fCommonHists(NULL),
   fCommonHistsRes(NULL),
   fHistQaQb(NULL),
   fHistQaQbNorm(NULL),
   fHistQaQbCos(NULL),
   fHistResolution(NULL),
   fHistQaNorm(NULL),
   fHistQaNormvsMa(NULL),
   fHistQbNorm(NULL),
   fHistQbNormvsMb(NULL),
   fHistMavsMb(NULL)

{
  // Constructor.
  fWeightsList = new TList();
  fHistList = new TList();
  
  // Initialize arrays:
  for(Int_t i=0;i<3;i++)
  {
   fHistSumOfWeightsPtRP[i] = NULL;
   fHistSumOfWeightsEtaRP[i] = NULL;
   fHistSumOfWeightsPtPOI[i] = NULL;
   fHistSumOfWeightsEtaPOI[i] = NULL;
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

  fHistProUQetaRP = new TProfile("Flow_UQetaRP_SP","Flow_UQetaRP_SP",iNbinsEta,dEtaMin,dEtaMax,"s");
  fHistProUQetaRP->SetXTitle("{eta}");
  fHistProUQetaRP->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQetaRP);

  fHistProUQetaPOI = new TProfile("Flow_UQetaPOI_SP","Flow_UQetaPOI_SP",iNbinsEta,dEtaMin,dEtaMax,"s");
  fHistProUQetaPOI->SetXTitle("{eta}");
  fHistProUQetaPOI->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQetaPOI);

  fHistProUQPtRP = new TProfile("Flow_UQPtRP_SP","Flow_UQPtRP_SP",iNbinsPt,dPtMin,dPtMax,"s");
  fHistProUQPtRP->SetXTitle("p_t (GeV)");
  fHistProUQPtRP->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQPtRP);

  fHistProUQPtPOI = new TProfile("Flow_UQPtPOI_SP","Flow_UQPtPOI_SP",iNbinsPt,dPtMin,dPtMax,"s");
  fHistProUQPtPOI->SetXTitle("p_t (GeV)");
  fHistProUQPtPOI->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQPtPOI);

  fHistProQaQb = new TProfile("Flow_QaQb_SP","Flow_QaQb_SP", 1, 0.5, 1.5,"s");
  fHistProQaQb->SetYTitle("<QaQb>");
  fHistList->Add(fHistProQaQb);

  fHistProQaQbNorm = new TProfile("FlowPro_QaQbNorm_SP","FlowPro_QaQbNorm_SP", 1, 0.5, 1.5,"s");
  fHistProQaQbNorm->SetYTitle("<QaQb/MaMb>");
  fHistList->Add(fHistProQaQbNorm);

  fHistSumOfLinearWeights = new TH1D("Flow_SumOfLinearWeights_SP","Flow_SumOfLinearWeights_SP",1,-0.5, 0.5);
  fHistSumOfLinearWeights -> SetYTitle("sum (*)");
  fHistSumOfLinearWeights -> SetXTitle("sum (Ma*Mb)");
  fHistList->Add(fHistSumOfLinearWeights);
  
  fHistSumOfQuadraticWeights = new TH1D("Flow_SumOfQuadraticWeights_SP","Flow_SumOfQuadraticWeights_SP",1,-0.5, 0.5);
  fHistSumOfQuadraticWeights -> SetYTitle("sum (*)");
  fHistSumOfQuadraticWeights -> SetXTitle("sum (Ma*Mb)^2");
  fHistList->Add(fHistSumOfQuadraticWeights);
  
  fHistProUQQaQbPtRP = new TProfile("Flow_UQQaQbPtRP_SP","Flow_UQQaQbPtRP_SP",iNbinsPt,dPtMin,dPtMax);
  fHistProUQQaQbPtRP -> SetYTitle("<*>");
  fHistProUQQaQbPtRP -> SetXTitle("<Qu QaQb>");
  fHistList->Add(fHistProUQQaQbPtRP);
  
  fHistProUQQaQbEtaRP = new TProfile("Flow_UQQaQbEtaRP_SP","Flow_UQQaQbEtaRP_SP",iNbinsEta,dEtaMin,dEtaMax);
  fHistProUQQaQbEtaRP -> SetYTitle("<*>");
  fHistProUQQaQbEtaRP -> SetXTitle("<Qu QaQb>");
  fHistList->Add(fHistProUQQaQbEtaRP);
  
  fHistProUQQaQbPtPOI = new TProfile("Flow_UQQaQbPtPOI_SP","Flow_UQQaQbPtPOI_SP",iNbinsPt,dPtMin,dPtMax);
  fHistProUQQaQbPtPOI -> SetYTitle("<*>");
  fHistProUQQaQbPtPOI -> SetXTitle("<Qu QaQb>");
  fHistList->Add(fHistProUQQaQbPtPOI);
  
  fHistProUQQaQbEtaPOI = new TProfile("Flow_UQQaQbEtaPOI_SP","Flow_UQQaQbEtaPOI_SP",iNbinsEta,dEtaMin,dEtaMax);
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
    
  fCommonHists = new AliFlowCommonHist("AliFlowCommonHistSP");
  fHistList->Add(fCommonHists);
  fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsSP");
  fHistList->Add(fCommonHistsRes);  

  fHistQaQb = new TH1D("Flow_QaQb_SP","Flow_QaQb_SP",200,-100.,100.);
  fHistQaQb -> SetYTitle("dN/dQaQb");
  fHistQaQb -> SetXTitle("QaQb");
  fHistList->Add(fHistQaQb);

  fHistQaQbNorm = new TH1D("Flow_QaQbNorm_SP","Flow_QaQbNorm_SP",44,-1.1,1.1);
  fHistQaQbNorm -> SetYTitle("dN/d(QaQb/MaMb)");
  fHistQaQbNorm -> SetXTitle("QaQb/MaMb");
  fHistList->Add(fHistQaQbNorm);

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

  TH1::AddDirectory(oldHistAddStatus);
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::Make(AliFlowEventSimple* anEvent) {

  //Calculate flow based on  <QaQb/MaMb> = <v^2>

  //Fill histograms
  if (anEvent) {

    //get Q vectors for the eta-subevents
    AliFlowVector* vQarray = new AliFlowVector[2];
    if (fUsePhiWeights) {
      anEvent->Get2Qsub(vQarray,2,fWeightsList,kTRUE);
    } else {
      anEvent->Get2Qsub(vQarray);
    }
    //subevent a
    AliFlowVector vQa = vQarray[0];
    //subevent b
    AliFlowVector vQb = vQarray[1];

    //check that the subevents are not empty:
    Double_t dMa = vQa.GetMult();
    Double_t dMb = vQb.GetMult();
    if (dMa != 0 && dMb != 0) {
      

      //request that the subevent multiplicities are not too different
      Double_t dRelDiff = TMath::Abs((dMa - dMb)/(dMa + dMb));
      if (dRelDiff < fRelDiffMsub) {

	//fill control histograms     
	fCommonHists->FillControlHistograms(anEvent);

	//fill some SP control histograms
	fHistProQaQb -> Fill(1.,vQa*vQb,1.); //Fill with weight 1 -> Weight with MaMb????
	fHistQaQbCos ->Fill(TMath::ACos((vQa/vQa.Mod())*(vQb/vQb.Mod())));  //Acos(Qa*Qb) = angle
	fHistResolution -> Fill(TMath::Cos( vQa.Phi()- vQb.Phi() ));  //vQa.Phi() returns 2*phi
	fHistQaQb -> Fill(vQa*vQb);
	fHistMavsMb -> Fill(dMb,dMa);

	//get total Q vector = the sum of subevent a and subevent b
	AliFlowVector vQ = vQa + vQb;
	//weight the Q vectors for the subevents by the multiplicity
	//Note: Weight Q only in the particle loop when it is clear if it should be (m-1) or M
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
	      //set default phi weight to 1
	      Double_t dW = 1.; 
	      //phi weight of pTrack
	      if(fUsePhiWeights && fPhiWeightsSub0 && fPhiWeightsSub1) {
		if (pTrack->InSubevent(0) ) {
		  Int_t iNBinsPhiSub0 = fPhiWeightsSub0->GetNbinsX();
		  dW = fPhiWeightsSub0->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*iNBinsPhiSub0/TMath::TwoPi())));  
		}
		else if ( pTrack->InSubevent(1)) { 
		  Int_t iNBinsPhiSub1 = fPhiWeightsSub1->GetNbinsX();
		  dW = fPhiWeightsSub1->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*iNBinsPhiSub1/TMath::TwoPi())));
		}
   		//bin = 1 + value*nbins/range
		//TMath::Floor rounds to the lower integer
	      }     
	    
	      //calculate vU
	      TVector2 vU;
	      Double_t dUX = TMath::Cos(2*dPhi);
	      Double_t dUY = TMath::Sin(2*dPhi);
	      vU.Set(dUX,dUY);
	      Double_t dModulus = vU.Mod();
	      if (dModulus!=0.) vU.Set(dUX/dModulus,dUY/dModulus);  // make length 1
	      else cerr<<"dModulus is zero!"<<endl;
	    
	      //redefine the Q vector and devide by its multiplicity
	      TVector2 vQm;
	      Double_t dQmX = 0.;
	      Double_t dQmY = 0.;
	      //subtract particle from the flowvector if used to define it
	      if (pTrack->InSubevent(0) || pTrack->InSubevent(1)) { 
		dQmX = (vQ.X() - dW*(pTrack->Weight())*dUX)/(dMq-1);
		dQmY = (vQ.Y() - dW*(pTrack->Weight())*dUY)/(dMq-1);
		vQm.Set(dQmX,dQmY);
	      
		//dUQ = scalar product of vU and vQm
		Double_t dUQ = (vU * vQm);
		Double_t dPt = pTrack->Pt();
		Double_t dEta = pTrack->Eta();
		//fill the profile histograms
		if (pTrack->InRPSelection()) {
		  fHistProUQetaRP -> Fill(dEta,dUQ,(dMq-1)); //Fill (Qu/(Mq-1)) with weight (Mq-1) 
		  //needed for the error calculation:
		  fHistProUQQaQbEtaRP -> Fill(dEta,dUQ*dQaQb,(dMq-1)*dMa*dMb); //Fill [Qu/(Mq-1)]*[QaQb/MaMb] with weight (Mq-1)MaMb	    
		  fHistProUQPtRP -> Fill(dPt,dUQ,(dMq-1));                     //Fill (Qu/(Mq-1)) with weight (Mq-1)
		  fHistProUQQaQbPtRP -> Fill(dPt,dUQ*dQaQb,(dMq-1)*dMa*dMb);   //Fill [Qu/(Mq-1)]*[QaQb/MaMb] with weight (Mq-1)MaMb	
		  
		  fHistSumOfWeightsEtaRP[0]->Fill(dEta,(dMq-1));        // sum of Mq-1     
		  fHistSumOfWeightsEtaRP[1]->Fill(dEta,pow((dMq-1),2.));// sum of (Mq-1)^2     
		  fHistSumOfWeightsEtaRP[2]->Fill(dEta,(dMq-1)*dMa*dMb);// sum of (Mq-1)*MaMb     
		  fHistSumOfWeightsPtRP[0]->Fill(dPt,(dMq-1));          // sum of Mq-1     
		  fHistSumOfWeightsPtRP[1]->Fill(dPt,pow((dMq-1),2.));  // sum of (Mq-1)^2     
		  fHistSumOfWeightsPtRP[2]->Fill(dPt,(dMq-1)*dMa*dMb);  // sum of (Mq-1)*MaMb     
		}
		if (pTrack->InPOISelection()) {
		  fHistProUQetaPOI -> Fill(dEta,dUQ,(dMq-1));//Fill (Qu/(Mq-1)) with weight (Mq-1)
		  //needed for the error calculation:
		  fHistProUQQaQbEtaPOI -> Fill(dEta,dUQ*dQaQb,(dMq-1)*dMa*dMb); //Fill [Qu/(Mq-1)]*[QaQb/MaMb] with weight (Mq-1)MaMb	    
		  fHistProUQPtPOI -> Fill(dPt,dUQ,(dMq-1));                     //Fill (Qu/(Mq-1)) with weight (Mq-1)
		  fHistProUQQaQbPtPOI -> Fill(dPt,dUQ*dQaQb,(dMq-1)*dMa*dMb);   //Fill [Qu/(Mq-1)]*[QaQb/MaMb] with weight (Mq-1)MaMb	    
		  
		  fHistSumOfWeightsEtaPOI[0]->Fill(dEta,(dMq-1));        // sum of Mq-1     
		  fHistSumOfWeightsEtaPOI[1]->Fill(dEta,pow((dMq-1),2.));// sum of (Mq-1)^2     
		  fHistSumOfWeightsEtaPOI[2]->Fill(dEta,(dMq-1)*dMa*dMb);// sum of (Mq-1)*MaMb     
		  fHistSumOfWeightsPtPOI[0]->Fill(dPt,(dMq-1));          // sum of Mq-1     
		  fHistSumOfWeightsPtPOI[1]->Fill(dPt,pow((dMq-1),2.)); // sum of (Mq-1)^2     
		  fHistSumOfWeightsPtPOI[2]->Fill(dPt,(dMq-1)*dMa*dMb); // sum of (Mq-1)*MaMb     
		}  
		
	      } else { //do not subtract the particle from the flowvector
		dQmX = vQ.X()/dMq;
		dQmY = vQ.Y()/dMq;
		vQm.Set(dQmX,dQmY);
	      
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
		}  
	      }//track not in subevents
	      
	    }//track
	    
	  }//loop over tracks
	
	fEventNumber++;

      } //difference Ma and Mb

    }// subevents not empty 
    delete [] vQarray;

  } //event

}//end of Make()

//--------------------------------------------------------------------  
void AliFlowAnalysisWithScalarProduct::GetOutputHistograms(TList *outputListHistos){
  
  //get pointers to all output histograms (called before Finish())

  if (outputListHistos) {
  //Get the common histograms from the output list
    AliFlowCommonHist *pCommonHist = dynamic_cast<AliFlowCommonHist*> 
      (outputListHistos->FindObject("AliFlowCommonHistSP"));
    AliFlowCommonHistResults *pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
      (outputListHistos->FindObject("AliFlowCommonHistResultsSP"));
    TProfile* pHistProQaQb     = dynamic_cast<TProfile*>(outputListHistos->FindObject("Flow_QaQb_SP"));
    TProfile* pHistProQaQbNorm = dynamic_cast<TProfile*>(outputListHistos->FindObject("FlowPro_QaQbNorm_SP"));
    TH1D*     pHistSumOfLinearWeights    = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_SumOfLinearWeights_SP"));
    TH1D*     pHistSumOfQuadraticWeights = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_SumOfQuadraticWeights_SP"));

    TProfile* pHistProUQetaRP  = dynamic_cast<TProfile*>(outputListHistos->FindObject("Flow_UQetaRP_SP"));
    TProfile* pHistProUQetaPOI = dynamic_cast<TProfile*>(outputListHistos->FindObject("Flow_UQetaPOI_SP"));
    TProfile* pHistProUQPtRP   = dynamic_cast<TProfile*>(outputListHistos->FindObject("Flow_UQPtRP_SP"));
    TProfile* pHistProUQPtPOI  = dynamic_cast<TProfile*>(outputListHistos->FindObject("Flow_UQPtPOI_SP"));
    TProfile* pHistProUQQaQbPtRP    = dynamic_cast<TProfile*>(outputListHistos->FindObject("Flow_UQQaQbPtRP_SP"));
    TProfile* pHistProUQQaQbEtaRP   = dynamic_cast<TProfile*>(outputListHistos->FindObject("Flow_UQQaQbEtaRP_SP"));
    TProfile* pHistProUQQaQbPtPOI   = dynamic_cast<TProfile*>(outputListHistos->FindObject("Flow_UQQaQbPtPOI_SP"));
    TProfile* pHistProUQQaQbEtaPOI  = dynamic_cast<TProfile*>(outputListHistos->FindObject("Flow_UQQaQbEtaPOI_SP"));
    TString weightFlag[3] = {"w_Qu_","w_Qu^2_","w_QuQaQb_"}; 
    TH1D* pHistSumOfWeightsPtRP[3] = {NULL};                    
    TH1D* pHistSumOfWeightsEtaRP[3] = {NULL};                    
    TH1D* pHistSumOfWeightsPtPOI[3] = {NULL};                    
    TH1D* pHistSumOfWeightsEtaPOI[3] = {NULL};                    
    for(Int_t i=0;i<3;i++)
    {
     pHistSumOfWeightsPtRP[i]   = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("Flow_SumOfWeights%sPtRP_SP",weightFlag[i].Data())));
     pHistSumOfWeightsEtaRP[i]  = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("Flow_SumOfWeights%sEtaRP_SP",weightFlag[i].Data())));
     pHistSumOfWeightsPtPOI[i]  = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("Flow_SumOfWeights%sPtPOI_SP",weightFlag[i].Data())));
     pHistSumOfWeightsEtaPOI[i] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("Flow_SumOfWeights%sEtaPOI_SP",weightFlag[i].Data())));
    }

    TH1D*     pHistQaQb     = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QaQb_SP"));
    TH1D*     pHistQaQbNorm = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QaQbNorm_SP"));
    TH1D*     pHistQaQbCos  = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QaQbCos_SP"));
    TH1D*     pHistResolution = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_resolution_SP"));
    TH1D*     pHistQaNorm   = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QaNorm_SP"));
    TH2D*     pHistQaNormvsMa   = dynamic_cast<TH2D*>(outputListHistos->FindObject("Flow_QaNormvsMa_SP"));
    TH1D*     pHistQbNorm   = dynamic_cast<TH1D*>(outputListHistos->FindObject("Flow_QbNorm_SP"));
    TH2D*     pHistQbNormvsMb   = dynamic_cast<TH2D*>(outputListHistos->FindObject("Flow_QbNormvsMb_SP"));
    TH2D*     pHistMavsMb = dynamic_cast<TH2D*>(outputListHistos->FindObject("Flow_MavsMb_SP"));

    //pass the pointers to the task
    if (pCommonHist && pCommonHistResults && pHistProQaQb && pHistProQaQbNorm &&
	pHistSumOfLinearWeights && pHistSumOfQuadraticWeights && 
	pHistQaQb && pHistQaQbNorm && pHistQaQbCos && pHistResolution &&
	pHistQaNorm && pHistQaNormvsMa && pHistQbNorm && pHistQbNormvsMb && 
	pHistMavsMb &&
	pHistProUQetaRP	&& pHistProUQetaPOI && 
	pHistProUQPtRP && pHistProUQPtPOI &&  
	pHistProUQQaQbPtRP && pHistProUQQaQbEtaRP && 
	pHistProUQQaQbPtPOI && pHistProUQQaQbEtaPOI) {
      this -> SetCommonHists(pCommonHist);
      this -> SetCommonHistsRes(pCommonHistResults);
      this -> SetHistProQaQb(pHistProQaQb);
      this -> SetHistProQaQbNorm(pHistProQaQbNorm);
      this -> SetHistSumOfLinearWeights(pHistSumOfLinearWeights);
      this -> SetHistSumOfQuadraticWeights(pHistSumOfQuadraticWeights); 
      this -> SetHistQaQb(pHistQaQb);
      this -> SetHistQaQbNorm(pHistQaQbNorm);
      this -> SetHistQaQbCos(pHistQaQbCos);
      this -> SetHistResolution(pHistResolution);
      this -> SetHistQaNorm(pHistQaNorm);
      this -> SetHistQaNormvsMa(pHistQaNormvsMa);
      this -> SetHistQbNorm(pHistQbNorm);
      this -> SetHistQbNormvsMb(pHistQbNormvsMb);
      this -> SetHistMavsMb(pHistMavsMb);
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

    } else {
      cout<<"WARNING: Histograms needed to run Finish() in SP are not accessable!"<<endl; }
         
  } // end of if(outputListHistos)
}            

//--------------------------------------------------------------------            
void AliFlowAnalysisWithScalarProduct::Finish() {
   
  //calculate flow and fill the AliFlowCommonHistResults
  if (fDebug) cout<<"AliFlowAnalysisWithScalarProduct::Finish()"<<endl;

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
  Double_t dV = -999.; 
  if(dQaQb>=0.)
  {
   dV = TMath::Sqrt(dQaQb); 
  }
  //statistical error of dQaQb: 
  //  statistical error = term1 * spread * term2:
  //  term1 = sqrt{sum_{i=1}^{N} w^2}/(sum_{i=1}^{N} w)
  //  term2 = 1/sqrt(1-term1^2) 
  Double_t dSpreadQaQb = fHistProQaQbNorm->GetBinError(1);
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
  fCommonHistsRes->FillIntegratedFlow(dV,dVerr);
  cout<<"dV = "<<dV<<" +- "<<dVerr<<endl;
	
  //Calculate differential flow and integrated flow (RP, POI)
  //---------------------------------------------------------
  //v as a function of eta for RP selection
  for(Int_t b=1;b<iNbinsEta+1;b++) {
    Double_t duQpro = fHistProUQetaRP->GetBinContent(b);
    Double_t dv2pro = -999.;
    if (dV!=0.) { dv2pro = duQpro/dV; }
    //calculate the statistical error
    Double_t dv2ProErr = CalculateStatisticalError(b, dStatErrorQaQb, fHistProUQetaRP, fHistProUQQaQbEtaRP, fHistSumOfWeightsEtaRP);
    //fill TH1D
    fCommonHistsRes->FillDifferentialFlowEtaRP(b, dv2pro, dv2ProErr);   
  } //loop over bins b
  

  //v as a function of eta for POI selection
  for(Int_t b=1;b<iNbinsEta+1;b++) {
    Double_t duQpro = fHistProUQetaPOI->GetBinContent(b);
    Double_t dv2pro = -999.;
    if (dV!=0.) { dv2pro = duQpro/dV; }
    //calculate the statistical error
    Double_t dv2ProErr = CalculateStatisticalError(b, dStatErrorQaQb, fHistProUQetaPOI, fHistProUQQaQbEtaPOI, fHistSumOfWeightsEtaPOI);
   
    //fill TH1D
    fCommonHistsRes->FillDifferentialFlowEtaPOI(b, dv2pro, dv2ProErr); 
  } //loop over bins b
  

  //v as a function of Pt for RP selection
  TH1F* fHistPtRP = fCommonHists->GetHistPtRP(); //for calculating integrated flow
  Double_t dVRP = 0.;
  Double_t dSumRP = 0.;
  Double_t dErrVRP =0.;
  
  for(Int_t b=1;b<iNbinsPt+1;b++) {
    Double_t duQpro = fHistProUQPtRP->GetBinContent(b);
    Double_t dv2pro = -999.;
    if (dV!=0.) { dv2pro = duQpro/dV; }
    //calculate the statistical error
    Double_t dv2ProErr = CalculateStatisticalError(b, dStatErrorQaQb, fHistProUQPtRP, fHistProUQQaQbPtRP, fHistSumOfWeightsPtRP);
              
    //fill TH1D
    fCommonHistsRes->FillDifferentialFlowPtRP(b, dv2pro, dv2ProErr);

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
  fCommonHistsRes->FillIntegratedFlowRP(dVRP,dErrVRP);
  cout<<"dV(RP) = "<<dVRP<<" +- "<<dErrVRP<<endl;
  

  //v as a function of Pt for POI selection 
  TH1F* fHistPtPOI = fCommonHists->GetHistPtPOI(); //for calculating integrated flow
  Double_t dVPOI = 0.;
  Double_t dSumPOI = 0.;
  Double_t dErrVPOI =0.;
  
  for(Int_t b=1;b<iNbinsPt+1;b++) {
    Double_t duQpro = fHistProUQPtPOI->GetBinContent(b);
    Double_t dv2pro = -999.;
    if (dV!=0.) { dv2pro = duQpro/dV; }
    //calculate the statistical error
    Double_t dv2ProErr = CalculateStatisticalError(b, dStatErrorQaQb, fHistProUQPtPOI, fHistProUQQaQbPtPOI, fHistSumOfWeightsPtPOI);
        
    //fill TH1D
    fCommonHistsRes->FillDifferentialFlowPtPOI(b, dv2pro, dv2ProErr); 
    
    //calculate integrated flow for POI selection
    if (fHistPtPOI){
      Double_t dYieldPt = fHistPtPOI->GetBinContent(b);
      dVPOI += dv2pro*dYieldPt;
      dSumPOI +=dYieldPt;
      dErrVPOI += dYieldPt*dYieldPt*dv2ProErr*dv2ProErr;
    } else { cout<<"fHistPtPOI is NULL"<<endl; }
  } //loop over bins b
  
  if (dSumPOI != 0.) {
    dVPOI /= dSumPOI; //the pt distribution should be normalised
    dErrVPOI /= (dSumPOI*dSumPOI);
    dErrVPOI = TMath::Sqrt(dErrVPOI);
  }
  fCommonHistsRes->FillIntegratedFlowPOI(dVPOI,dErrVPOI);
  cout<<"dV(POI) = "<<dVPOI<<" +- "<<dErrVPOI<<endl;

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
