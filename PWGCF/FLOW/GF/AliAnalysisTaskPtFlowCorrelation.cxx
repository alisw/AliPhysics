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
#include <iostream>

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
#include "AliAnalysisTaskPtFlowCorrelation.h"

#include "TComplex.h"
#include "TMath.h"

#include "AliMultSelection.h"



class AliAnalysisTaskPtFlowCorrelation;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskPtFlowCorrelation) // classimp: necessary for root

AliAnalysisTaskPtFlowCorrelation::AliAnalysisTaskPtFlowCorrelation() : AliAnalysisTaskSE(), 
    fAOD(0), 
    fOutputList(0), 
    fProfileRho(0),
    fProfileCov(0),
    fProfileVar(0),
    fProfileC2(0),
    inputFile(0),
    inputProfileMeanPt(0),
    inputAxisMeanPt(0),
    inputHistPhi(0),
    inputAxisPhi(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskPtFlowCorrelation::AliAnalysisTaskPtFlowCorrelation(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), 
    fOutputList(0), 
    fProfileRho(0),
    fProfileCov(0),
    fProfileVar(0),
    fProfileC2(0),
    inputFile(0),
    inputProfileMeanPt(0),
    inputAxisMeanPt(0),
    inputHistPhi(0),
    inputAxisPhi(0)
{
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events             
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 

}
//_____________________________________________________________________________
AliAnalysisTaskPtFlowCorrelation::~AliAnalysisTaskPtFlowCorrelation()
{
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskPtFlowCorrelation::UserCreateOutputObjects()
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

    // Save output
    fOutputList->Add(fProfileRho);
    fOutputList->Add(fProfileVar);
    fOutputList->Add(fProfileCov);
    fOutputList->Add(fProfileC2);

    PostData(1, fOutputList);         


        
    //inputFile = TFile::Open("alien:///alice/cern.ch/user/f/frjensen/Stat-Pass2-2015o.root", "READ");  
    inputFile = new TFile("Stat-Pass2-2015o.root", "READ"); 

    TDirectoryFile* statDirectory = (TDirectoryFile*)inputFile->Get("GetPassStatistic");
    if( statDirectory == nullptr ) printf("Could not load directory \n");

    TList*  statList = (TList*)statDirectory->Get("MyOutputContainer");
    if( statList == nullptr ) printf("\n Could not load TList \n");

    inputProfileMeanPt = (TProfile*)statList->FindObject("fProfileMeanPt");
    if( inputProfileMeanPt == nullptr ) printf("\n Could not load pt profile \n");

    inputAxisMeanPt =  (TAxis*)inputProfileMeanPt->GetXaxis();
    if( inputAxisMeanPt == nullptr ) printf("\n Could not load axis \n");
  

    // For calculations of phi-weight
    inputHistPhi = (TH1D*)statList->FindObject("fHistPhi");
    if( inputHistPhi == nullptr ) printf("\n Could not load pt profile \n");

    inputAxisPhi =  (TAxis*)inputHistPhi->GetXaxis();
    if( inputAxisPhi == nullptr ) printf("\n Could not load axis \n");
    N_max = inputHistPhi->GetBinContent( inputHistPhi->GetMaximumBin() );

}
//_____________________________________________________________________________
void AliAnalysisTaskPtFlowCorrelation::UserExec(Option_t *)
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

    Double_t sumWeightSubB = 0;
    Double_t sumWeightSubB_sq = 0;
    Double_t sumWeightSubC = 0;
    Double_t sumWeightSubC_sq = 0;


    Double_t p11 = 0;
    Double_t p22 = 0;

    Int_t iTracks(fAOD->GetNumberOfTracks());      

    Double_t meanPt = GetMeanPt(centrality);

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
        Double_t phiWeight = GetPhiWeight(phi);

        if( std::abs(eta) < 0.4)            
        {
            // Fill sub-event B for mean pt
            p11 += ptWeight*pt;
            p22 += TMath::Power(pt*ptWeight, 2.);
            sumWeightSubB += ptWeight;    
            sumWeightSubB_sq += ptWeight*ptWeight;    

        }
        else if( 0.4 < eta && eta < 0.8)
        {
            // Fill sub-event C for vn
            QSubC += phiWeight*TComplex(TMath::Cos(2.*phi), TMath::Sin(2.*phi));
            QSubCFour += phiWeight*TComplex(TMath::Cos(4.*phi), TMath::Sin(4.*phi));
            sumWeightSubC += phiWeight;
            sumWeightSubC_sq += phiWeight*phiWeight;

        }
        else if( -0.8 < eta && eta < -0.4)
        {
            // Fill sub-event A for vn
            QSubA += phiWeight*TComplex(TMath::Cos(2.*phi), TMath::Sin(2.*phi));
            QSubAFour += phiWeight*TComplex(TMath::Cos(4.*phi), TMath::Sin(4.*phi));
            sumWeightSubA += phiWeight;
            sumWeightSubA_sq += phiWeight*phiWeight;
        }
 
    }

    if(sumWeightSubA < 2) return;
    if(sumWeightSubB < 2) return;
    if(sumWeightSubC < 2) return;

    Double_t tau = sumWeightSubB_sq/TMath::Power(sumWeightSubB, 2.);

    p11 *= 1./sumWeightSubB; 
    p22 *= 1./TMath::Power(sumWeightSubB, 2.);

    PSubB = TComplex(p11, 0.);

    //Dynamic fluctuations non-central to central
    Double_t oneParPtCentralCorrelation = p11 - meanPt;
    Double_t twoParPtCentralCorrelation = p22 - 2.*p11*meanPt + meanPt*meanPt;
    Double_t c2 = 1+(TMath::Power(oneParPtCentralCorrelation, 2.) - twoParPtCentralCorrelation)/(1.-tau);

    // Internal subevent correlation -> Useed for four particle correlation later;
    Double_t twoParCorrSubA = ((QSubA*TComplex::Conjugate(QSubA)).Re() - sumWeightSubA_sq)/(Double_t)(sumWeightSubA*sumWeightSubA - sumWeightSubC_sq); 
    Double_t twoParCorrSubC = ((QSubC*TComplex::Conjugate(QSubC)).Re() - sumWeightSubC_sq)/(Double_t)(sumWeightSubC*sumWeightSubC - sumWeightSubC_sq); 
    
    // Multi-particle correlation
    Double_t twoParCorr     = (QSubA*TComplex::Conjugate(QSubC)).Re();             // v2^2
    Double_t threeParCorr   = (QSubA * TComplex::Conjugate(QSubC) * PSubB).Re();     // v2^2*pT
    Double_t fourParCorr    = twoParCorrSubA * TComplex::Conjugate(twoParCorrSubC);  // v2^4

    // // Update with event weight
    twoParCorr      *= 1./(Double_t)(sumWeightSubA * sumWeightSubC);
    threeParCorr    *= 1./(Double_t)(sumWeightSubA * sumWeightSubB * sumWeightSubC);

    Double_t twoParCum      = twoParCorr;
    Double_t threeParCum    = threeParCorr - twoParCorr*meanPt;
    Double_t fourParCum     = fourParCorr - 2.*twoParCorr;
    Double_t var            = twoParCum*twoParCum + fourParCum;

    if(var < 0 || c2 < 0) return;

    Double_t rho = threeParCum/(TMath::Sqrt(var)*TMath::Sqrt(c2));

    fProfileRho->Fill(centrality, rho);
    fProfileVar->Fill(centrality, var);
    fProfileCov->Fill(centrality, threeParCum);
    fProfileC2->Fill(centrality, c2);

    
    PostData(1, fOutputList);                           
                                                       
}
//_____________________________________________________________________________
void AliAnalysisTaskPtFlowCorrelation::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________




Double_t AliAnalysisTaskPtFlowCorrelation::GetMeanPt(Double_t centrality)
{
    Int_t binx = inputAxisMeanPt->FindBin(centrality);
    return inputProfileMeanPt->GetBinContent(binx);
}


Double_t AliAnalysisTaskPtFlowCorrelation::GetPhiWeight(Double_t phi)
{
    Int_t binx = inputAxisPhi->FindBin(phi);
    Double_t w = N_max/inputHistPhi->GetBinContent(binx);
    
    return w;
}
