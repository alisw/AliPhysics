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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */
#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectoryFile.h"

#include "TAxis.h"

#include "TChain.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskMomentumFlowCorrelation.h"

#include "TComplex.h"
#include "TMath.h"

#include "AliMultSelection.h"



class AliAnalysisTaskMomentumFlowCorrelation;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskMomentumFlowCorrelation) // classimp: necessary for root

AliAnalysisTaskMomentumFlowCorrelation::AliAnalysisTaskMomentumFlowCorrelation() : AliAnalysisTaskSE(), 
    fAOD(0), 
    fOutputList(0), 
    fProfileRho(0),
    fProfileCov(0),
    fProfileVar(0),
    fProfileC2(0),
    fProfileTwoParCorr(0),
    fProfileThreeParCorr(0),
    fProfileFourParCorr(0),
    fProfileMeanPt(0),
    fProfileTwoParPtCorr(0),
    fProfileTwoParPtCorrCent(0),
    fProfileOneParPtCorrCent(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskMomentumFlowCorrelation::AliAnalysisTaskMomentumFlowCorrelation(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), 
    fOutputList(0), 
    fProfileRho(0),
    fProfileCov(0),
    fProfileVar(0),
    fProfileC2(0),
    fProfileTwoParCorr(0),
    fProfileThreeParCorr(0),
    fProfileFourParCorr(0),
    fProfileMeanPt(0),
    fProfileTwoParPtCorr(0),
    fProfileTwoParPtCorrCent(0),
    fProfileOneParPtCorrCent(0)
{
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events             
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 

}
//_____________________________________________________________________________
AliAnalysisTaskMomentumFlowCorrelation::~AliAnalysisTaskMomentumFlowCorrelation()
{
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskMomentumFlowCorrelation::UserCreateOutputObjects()
{
    // create output objects
    fOutputList = new TList();          
    fOutputList->SetOwner(kTRUE);      
    
    const int BINWIDTH  = 5;
    const int BINMIN    = 0;
    const int BINMAX    = 100;
    const int NBINS = (int)((BINMAX-BINMIN)/BINWIDTH);
    double edges[NBINS+1];

    for (int i = 0; i <= NBINS; i++)
    {
        edges[i] = i*BINWIDTH;
    }

    // Create Profiles for centrality dependent datas
    fProfileRho = new TProfile("Rho", "Rho", NBINS, edges);
    fProfileVar = new TProfile("Var", "Var(v_2{2}^2)", NBINS, edges);
    fProfileCov = new TProfile("Cov", "Cov(v_2^2, [pt])", NBINS, edges);
    fProfileC2  = new TProfile("C2", "Dynamic pt variance", NBINS, edges);

    fProfileTwoParCorr       = new TProfile("Two_par_v2^2", "Subevent method", NBINS, edges);
    fProfileThreeParCorr     = new TProfile("Three_par_meanpt_vn", "Correlation of v2^2 and [pT]", NBINS, edges);
    fProfileFourParCorr      = new TProfile("Four_par_v2^4", "Sub event method", NBINS, edges);

    fProfileMeanPt  = new TProfile("MeanPt", "Mean Pt", NBINS, edges);
    fProfileTwoParPtCorr  = new TProfile("2par_pt_corr", "Two particle pt correlation", NBINS, edges);
    fProfileTwoParPtCorrCent  = new TProfile("2par_pt_corr_cent", "Central two particle pt correlation", NBINS, edges);
    fProfileOneParPtCorrCent  = new TProfile("1par_pt_corr_cent", "Central one particle pt correlation", NBINS, edges);

    // Save output
    fOutputList->Add(fProfileRho);
    fOutputList->Add(fProfileVar);
    fOutputList->Add(fProfileCov);
    fOutputList->Add(fProfileC2);
    fOutputList->Add(fProfileMeanPt);
    fOutputList->Add(fProfileTwoParPtCorr);
    fOutputList->Add(fProfileTwoParPtCorrCent);
    fOutputList->Add(fProfileOneParPtCorrCent);

    fOutputList->Add(fProfileTwoParCorr);
    fOutputList->Add(fProfileThreeParCorr);
    fOutputList->Add(fProfileFourParCorr);

    PostData(1, fOutputList);         



}
//_____________________________________________________________________________
void AliAnalysisTaskMomentumFlowCorrelation::UserExec(Option_t *)
{
    // user exec this function is called once for each event
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    

    if(!fAOD) return;                                 

    Double_t vertexZ(fAOD->GetPrimaryVertex()->GetZ());
    if(!(std::abs(vertexZ) < 10.)) return;

    Double_t centrality(0);
    AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(multSelection) centrality = multSelection->GetMultiplicityPercentile("V0M");

    TComplex QSubA = TComplex(0., 0.);
    TComplex QSubAFour = TComplex(0., 0.);

    TComplex QSubC = TComplex(0., 0.);
    TComplex QSubCFour = TComplex(0., 0.);

    TComplex PSubB = TComplex(0., 0.);

    Double_t sumWeightSubA = 0;
    Double_t sumWeightSubA_sq = 0;

    Double_t w1 = 0;
    Double_t w2 = 0;

    Double_t sumWeightSubC = 0;
    Double_t sumWeightSubC_sq = 0;

    Double_t p11 = 0;
    Double_t p22 = 0;

    Int_t iTracks(fAOD->GetNumberOfTracks());      

    //Double_t meanPt = GetMeanPt(centrality);

    for(Int_t i(0); i < iTracks; i++) {                 
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track || !track->TestFilterBit(96)) continue;                            

        // Kinematic cuts
        Double_t pt = track->Pt();
        if(!(0.2 < pt && pt < 3.0)) continue;

        Double_t eta = track->Eta();
        if(!(std::abs(eta) < 0.8)) continue;

        Double_t phi = track->Phi();
        
        Double_t ptWeight = 1.0;
        Double_t phiWeight = 1.0;

        if( std::abs(eta) < 0.4)            
        {
            // Fill sub-event B for mean pt
            p11 += pt*ptWeight;
            p22 += TMath::Power(pt*ptWeight, 2.);
            w1 += ptWeight;    
            w2 += ptWeight*ptWeight;    
        }
        else if( 0.4 < eta && eta < 0.8)
        {
            // Fill sub-event C for vn
            QSubC += phiWeight*TComplex(TMath::Cos(2.*phi), TMath::Sin(2.*phi));
            QSubCFour += phiWeight*phiWeight*TComplex(TMath::Cos(4.*phi), TMath::Sin(4.*phi));
            sumWeightSubC += phiWeight;
            sumWeightSubC_sq += phiWeight*phiWeight;
        }
        else if( -0.8 < eta && eta < -0.4)
        {
            // Fill sub-event A for vn
            QSubA += phiWeight*TComplex(TMath::Cos(2.*phi), TMath::Sin(2.*phi));
            QSubAFour += phiWeight*phiWeight*TComplex(TMath::Cos(4.*phi), TMath::Sin(4.*phi));

            sumWeightSubA += phiWeight;
            sumWeightSubA_sq += phiWeight*phiWeight;
        }
 
    }
    if(w1 < 4) return;
    if(sumWeightSubA < 4) return;
    if(sumWeightSubC < 4) return;

    PSubB = TComplex(p11, 0.);
    
    Double_t tau = w2/TMath::Power(w1, 2.);

    // Dynamic fluctuations non-central to central
    Double_t meanPt = (p11/w1);

    Double_t oneParPtCentralCorrelation = (p11/w1);
    Double_t twoParPtCentralCorrelation = (p22/w2) - 2.*(p11/w1)*meanPt + meanPt*meanPt;
    Double_t c2_central = TMath::Power(oneParPtCentralCorrelation, 2.) - twoParPtCentralCorrelation;
    c2_central *= 1./(1-tau);

    Double_t oneParPtCorrelation = p11/w1;
    Double_t twoParPtCorrelation = (p11*p11 - p22)/(w1*w1-w2);
    Double_t c2 = - twoParPtCorrelation + oneParPtCorrelation*oneParPtCorrelation;

    fProfileC2->Fill(centrality, TMath::Sqrt(c2));
    fProfileMeanPt->Fill(centrality, p11/w1);
    fProfileTwoParPtCorr->Fill(centrality, twoParPtCorrelation);
    fProfileOneParPtCorrCent->Fill(centrality, oneParPtCentralCorrelation);
    fProfileTwoParPtCorrCent->Fill(centrality, twoParPtCentralCorrelation);


    // Internal subevent correlation -> Useed for four particle correlation later;
    TComplex twoParCorrSubA = QSubA*QSubA - QSubAFour;
    TComplex twoParCorrSubC = TComplex::Conjugate(QSubC)*TComplex::Conjugate(QSubC) - TComplex::Conjugate(QSubCFour);

    twoParCorrSubA *= 1./(sumWeightSubA*sumWeightSubA - sumWeightSubA_sq);
    twoParCorrSubC *= 1./(sumWeightSubC*sumWeightSubC - sumWeightSubC_sq);

    // Multi-particle correlation
    Double_t twoParCorr     = (QSubA*TComplex::Conjugate(QSubC)).Re();                       // <2>
    Double_t threeParCorr   = (QSubA*TComplex::Conjugate(QSubC) * PSubB).Re();               // <v2^2*pT>
    Double_t fourParCorr    = (twoParCorrSubA * twoParCorrSubC).Re();                        // <4>

    // Update with event weight
    twoParCorr      *= 1./(sumWeightSubA * sumWeightSubC);
    threeParCorr    *= 1./(sumWeightSubA * w1 * sumWeightSubC);

    Double_t twoParCum      = twoParCorr;
    Double_t threeParCum    = threeParCorr - twoParCorr*oneParPtCorrelation;
    Double_t fourParCum     = fourParCorr - 2.*twoParCorr*twoParCorr;
    Double_t var            = twoParCum*twoParCum + fourParCum;                               // v2{2}^4 -  v2{4}^4 

    fProfileVar->Fill(centrality, var);
    fProfileCov->Fill(centrality, threeParCum);

    fProfileTwoParCorr->Fill(centrality, twoParCorr );
    fProfileThreeParCorr->Fill(centrality, threeParCorr );
    fProfileFourParCorr->Fill(centrality, fourParCorr );

   
    if(var < 0 || c2 < 0) return;

    Double_t rho = threeParCum/(TMath::Sqrt(var)*TMath::Sqrt(c2));

    fProfileRho->Fill(centrality, rho);


    
    PostData(1, fOutputList);                           
                                                       
}
//_____________________________________________________________________________
void AliAnalysisTaskMomentumFlowCorrelation::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________




    // printf("\n twoParCum: %f \n", twoParCum);
    // printf("fourParCum: %f \n", fourParCum);
    // printf("threeParCorr: %f \n", threeParCorr);
    // printf("var: %f \n ", var);
    // printf("c2: %f \n ", c2);
    // printf("RHO: %f \n ", rho);
    // Double_t twoParCorrSubA = (TMath::Power(TComplex::Abs(QSubA), 2) - sumWeightSubA_sq)/(Double_t)(sumWeightSubA*sumWeightSubA - sumWeightSubC_sq); 
    // Double_t twoParCorrSubC = (TMath::Power(TComplex::Abs(QSubC), 2) - sumWeightSubC_sq)/(Double_t)(sumWeightSubC*sumWeightSubC - sumWeightSubC_sq); 
    