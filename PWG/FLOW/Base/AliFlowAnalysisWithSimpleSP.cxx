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

#define ALIFLOWANALYSISWITHSIMPLESP_CXX

#include "TFile.h"      
#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"

#include "AliFlowCommonConstants.h"
#include "AliFlowEventSimple.h"
#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowAnalysisWithSimpleSP.h"

using std::cout;
using std::endl;


ClassImp(AliFlowAnalysisWithSimpleSP)

    //-----------------------------------------------------------------------
    AliFlowAnalysisWithSimpleSP::AliFlowAnalysisWithSimpleSP():
        fDebug(0),
        fMinimalBook(kFALSE),
        fUsePhiWeights(0),
        fApplyCorrectionForNUA(0),
        fHarmonic(2),
        fWeights(kFALSE),
        fScaling(kTRUE),
        fNormalizationType(1),
        fV0SanityCheck(0),
        fExternalResolution(-999.),
        fExternalResErr(0.),
        fTotalQvector(3),
        fWeightsList(NULL),
        fHistList(NULL),
        fHistProConfig(NULL),
        fHistProQaQbNorm(NULL),
        fHistSumOfWeights(NULL),
        fHistProNUAq(NULL),
        fCentralityWeight(1.),
        fHistProQNorm(NULL),
        fHistProQaQb(NULL),
        fHistProQaQbM(NULL),
        fHistMaMb(NULL),
        fHistQNormQaQbNorm(NULL),
        fHistQaNormMa(NULL),
        fHistQbNormMb(NULL),
        fResolution(NULL),
        fHistQaQb(NULL),
        fHistQaQc(NULL),
        fHistQbQc(NULL),
        fHistQaQbCos(NULL),
        fHistNumberOfSubtractedDaughters(NULL),
        fCommonHists(NULL),
        fCommonHistsuQ(NULL),
        fCommonHistsRes(NULL)
{
    //ctor
    for(int i=0; i!=2; ++i) {
        fPhiWeightsSub[i] = NULL;
        for(int j=0; j!=2; ++j) {
            fHistProUQ[i][j] = NULL;
            fHistProUQQaQb[i][j] = NULL;
            fHistProNUAu[i][j][0] = NULL;
            fHistProNUAu[i][j][1] = NULL;
            for(int k=0; k!=3; ++k)
                fHistSumOfWeightsu[i][j][k] = NULL;
        }
    }
}
//-----------------------------------------------------------------------
AliFlowAnalysisWithSimpleSP::~AliFlowAnalysisWithSimpleSP() 
{
    //destructor
    delete fWeightsList;
    delete fHistList;
}
//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleSP::Init() {
    //Define all histograms
    //printf("---Analysis with the Scalar Product Method--- Init\n");
    //printf("--- fNormalizationType %d ---\n", fNormalizationType);
    //save old value and prevent histograms from being added to directory
    //to avoid name clashes in case multiple analaysis objects are used
    //in an analysis
    Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    fHistList = new TList();
    fHistList->SetName("cobjSP");
    fHistList->SetOwner();

    TList *uQRelated = new TList();
    uQRelated->SetName("uQ");
    uQRelated->SetOwner();

    TList *nuaRelated = new TList();
    nuaRelated->SetName("NUA");
    nuaRelated->SetOwner();

    TList *errorRelated = new TList();
    errorRelated->SetName("error");
    errorRelated->SetOwner();

    TList *tQARelated = new TList();
    tQARelated->SetName("QA");
    tQARelated->SetOwner();

    fCommonHists = new AliFlowCommonHist("AliFlowCommonHist_SP","AliFlowCommonHist",fMinimalBook, fHarmonic);
    (fCommonHists->GetHarmonic())->Fill(0.5,fHarmonic); // store harmonic 
    fHistList->Add(fCommonHists);

    //  fCommonHistsuQ = new AliFlowCommonHist("AliFlowCommonHist_uQ");
    //  (fCommonHistsuQ->GetHarmonic())->Fill(0.5,fHarmonic); // store harmonic 
    //  fHistList->Add(fCommonHistsuQ);

    fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResults_SP","",fHarmonic);
    fHistList->Add(fCommonHistsRes);

    fHistProConfig = new TProfile("FlowPro_Flags_SP","Flow_Flags_SP",4,0.5,4.5,"s");
    fHistProConfig->GetXaxis()->SetBinLabel(1,"fApplyCorrectionForNUA");
    fHistProConfig->GetXaxis()->SetBinLabel(2,"fNormalizationType");
    fHistProConfig->GetXaxis()->SetBinLabel(3,"fUsePhiWeights");
    fHistProConfig->GetXaxis()->SetBinLabel(4,"fHarmonic");
    fHistProConfig->Fill(1,fApplyCorrectionForNUA);
    fHistProConfig->Fill(2,fNormalizationType);
    fHistProConfig->Fill(3,fUsePhiWeights);
    fHistProConfig->Fill(4,fHarmonic);
    fHistList->Add(fHistProConfig);

    fHistProQaQbNorm = new TProfile("FlowPro_QaQbNorm_SP","FlowPro_QaQbNorm_SP", 10000, -1000, 1000);
    fHistProQaQbNorm->SetYTitle("<QaQb/Na/Nb>");
    errorRelated->Add(fHistProQaQbNorm);

    fHistProNUAq = new TProfile("FlowPro_NUAq_SP","FlowPro_NUAq_SP", 6, 0.5, 6.5,"s");
    fHistProNUAq->GetXaxis()->SetBinLabel( 1,"<<sin(#Phi_{a})>>");
    fHistProNUAq->GetXaxis()->SetBinLabel( 2,"<<cos(#Phi_{a})>>");
    fHistProNUAq->GetXaxis()->SetBinLabel( 3,"<<sin(#Phi_{b})>>");
    fHistProNUAq->GetXaxis()->SetBinLabel( 4,"<<cos(#Phi_{b})>>");
    fHistProNUAq->GetXaxis()->SetBinLabel( 5,"<<sin(#Phi_{t})>>");
    fHistProNUAq->GetXaxis()->SetBinLabel( 6,"<<cos(#Phi_{t})>>");
    nuaRelated->Add(fHistProNUAq);

    fHistSumOfWeights = new TH1D("Flow_SumOfWeights_SP","Flow_SumOfWeights_SP",2,0.5, 2.5);
    fHistSumOfWeights->GetXaxis()->SetBinLabel( 1,"#Sigma Na*Nb");
    fHistSumOfWeights->GetXaxis()->SetBinLabel( 2,"#Sigma (Na*Nb)^2");
    errorRelated->Add(fHistSumOfWeights);

    TString sPOI[2] = {"RP","POI"}; // backward compatibility
    TString sEta[2] = {"Pt","eta"}; // backward compatibility
    TString sTitle[2] = {"p_{T} [GeV]","#eta"};
    TString sWeights[3] = {"uQ","uQuQ","uQQaQb"};
    Int_t iNbins[2];
    Double_t dMin[2], dMax[2];
    iNbins[0] = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
    iNbins[1] = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
    dMin[0]   = AliFlowCommonConstants::GetMaster()->GetPtMin();
    dMin[1]   = AliFlowCommonConstants::GetMaster()->GetEtaMin();
    dMax[0]   = AliFlowCommonConstants::GetMaster()->GetPtMax();
    dMax[1]   = AliFlowCommonConstants::GetMaster()->GetEtaMax();
    for(Int_t iPOI=0; iPOI!=2; ++iPOI) 
        for(Int_t iSpace=0; iSpace!=2; ++iSpace) {
            // uQ
            fHistProUQ[iPOI][iSpace] = new TProfile( Form( "FlowPro_UQ_%s%s_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ),
                    Form( "FlowPro_UQ%s%s_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ),
                    iNbins[iSpace], dMin[iSpace], dMax[iSpace], "s");
            fHistProUQ[iPOI][iSpace]->SetXTitle(sTitle[iSpace].Data());
            fHistProUQ[iPOI][iSpace]->SetYTitle("<uQ>");
            uQRelated->Add(fHistProUQ[iPOI][iSpace]);

            // NUAu
            fHistProNUAu[iPOI][iSpace][0] = new TProfile( Form("FlowPro_NUAu_%s%s_IM_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ),
                    Form("FlowPro_NUAu_%s%s_IM_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ),
                    iNbins[iSpace], dMin[iSpace], dMax[iSpace]);
            fHistProNUAu[iPOI][iSpace][0]->SetXTitle(sTitle[iSpace].Data());
            nuaRelated->Add(fHistProNUAu[iPOI][iSpace][0]);
            fHistProNUAu[iPOI][iSpace][1] = new TProfile( Form("FlowPro_NUAu_%s%s_RE_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ),
                    Form("FlowPro_NUAu_%s%s_RE_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ),
                    iNbins[iSpace], dMin[iSpace], dMax[iSpace]);
            fHistProNUAu[iPOI][iSpace][1]->SetXTitle(sTitle[iSpace].Data());
            nuaRelated->Add(fHistProNUAu[iPOI][iSpace][1]);

            // uQ QaQb
            fHistProUQQaQb[iPOI][iSpace] = new TProfile( Form("FlowPro_UQQaQb_%s%s_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ),
                    Form("FlowPro_UQQaQb_%s%s_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ),
                    iNbins[iSpace], dMin[iSpace], dMax[iSpace]);
            fHistProUQQaQb[iPOI][iSpace]->SetXTitle(sTitle[iSpace].Data());
            fHistProUQQaQb[iPOI][iSpace]-> SetYTitle("<Qu QaQb>");
            errorRelated->Add(fHistProUQQaQb[iPOI][iSpace]);

            // uWeights
            for(Int_t i=0; i!=3; ++i) {
                fHistSumOfWeightsu[iPOI][iSpace][i] = new TH1D( Form("Flow_SumOfWeights_%s%s_%s_SP",sWeights[i].Data(),sPOI[iPOI].Data(),sEta[iSpace].Data()),
                        Form("Flow_SumOfWeights_%s%s_%s_SP",sWeights[i].Data(),sPOI[iPOI].Data(),sEta[iSpace].Data()),
                        iNbins[iSpace], dMin[iSpace], dMax[iSpace]);
                fHistSumOfWeightsu[iPOI][iSpace][i]->SetYTitle(sWeights[i].Data());
                fHistSumOfWeightsu[iPOI][iSpace][i]->SetXTitle(sTitle[iSpace].Data());
                errorRelated->Add(fHistSumOfWeightsu[iPOI][iSpace][i]);
            }
        }
    //weights
    if(fUsePhiWeights) {
        if(!fWeightsList) {
            printf( "WARNING: fWeightsList is NULL in the Scalar Product method.\n" );
            exit(0);
        }
        fPhiWeightsSub[0] = dynamic_cast<TH1F*>
            (fWeightsList->FindObject("phi_weights_sub0"));
        if(!fPhiWeightsSub[0]) {
            printf( "WARNING: phi_weights_sub0 not found in the Scalar Product method.\n" );
            exit(0);  
        }
        nuaRelated->Add( fPhiWeightsSub[0] );
        fPhiWeightsSub[1] = dynamic_cast<TH1F*>
            (fWeightsList->FindObject("phi_weights_sub1"));
        if(!fPhiWeightsSub[1]) {
            printf( "WARNING: phi_weights_sub1 not found in the Scalar Product method.\n" );
            exit(0);  
        }
        nuaRelated->Add( fPhiWeightsSub[1] );
    }

    if(!fMinimalBook) {
        fHistProQNorm = new TProfile("FlowPro_QNorm_SP","FlowPro_QNorm_SP",       1,0.5,1.5,"s");
        fHistProQNorm->SetYTitle("<|Qa+Qb|>");
        tQARelated->Add(fHistProQNorm);

        fHistProQaQb  = new TH1D("FlowPro_QaQb_SP","FlowPro_QaQb_SP", 10000,-100,100);
        fHistProQaQb->SetYTitle("<QaQb>");
        fHistProQaQb->StatOverflows(kTRUE); 
        tQARelated->Add(fHistProQaQb);

        fHistProQaQbM = new TH1D("FlowPro_QaQbvsM_SP","FlowPro_QaQbvsM_SP",1000,0.0,10000);
        fHistProQaQbM->SetYTitle("<QaQb>");
        fHistProQaQbM->SetXTitle("M");
        fHistProQaQbM->Sumw2();
        tQARelated->Add(fHistProQaQbM);

        fHistMaMb = new TH2D("Flow_MavsMb_SP","Flow_MavsMb_SP",100,0.,100.,100,0.,100.);
        fHistMaMb->SetYTitle("Ma");
        fHistMaMb->SetXTitle("Mb");
        tQARelated->Add(fHistMaMb);

        fHistQNormQaQbNorm = new TH2D("Flow_QNormvsQaQbNorm_SP","Flow_QNormvsQaQbNorm_SP",88,-1.1,1.1,22,0.,1.1);
        fHistQNormQaQbNorm->SetYTitle("|Q/Mq|");
        fHistQNormQaQbNorm->SetXTitle("QaQb/MaMb");
        tQARelated->Add(fHistQNormQaQbNorm);

        fHistQaNormMa = new TH2D("Flow_QaNormvsMa_SP","Flow_QaNormvsMa_SP",100,0.,100.,22,0.,1.1);
        fHistQaNormMa->SetYTitle("|Qa/Ma|");
        fHistQaNormMa->SetXTitle("Ma");
        tQARelated->Add(fHistQaNormMa);

        fHistQbNormMb = new TH2D("Flow_QbNormvsMb_SP","Flow_QbNormvsMb_SP",100,0.,100.,22,0.,1.1);
        fHistQbNormMb->SetYTitle("|Qb/Mb|");
        fHistQbNormMb->SetXTitle("Mb");
        tQARelated->Add(fHistQbNormMb);

        fResolution = new TH1D("Flow_resolution_SP","Flow_resolution_SP",100,-1.0,1.0);
        fResolution->SetYTitle("dN/d(Cos2(#phi_a - #phi_b))");
        fResolution->SetXTitle("Cos2(#phi_a - #phi_b)");
        tQARelated->Add(fResolution);

        fHistQaQb = new TH1D("Flow_QaQb_SP","Flow_QaQb_SP",20000,-1000.,1000.);
        fHistQaQb->SetYTitle("dN/dQaQb");
        fHistQaQb->SetXTitle("dQaQb");
        fHistQaQb->StatOverflows(kTRUE);
        tQARelated->Add(fHistQaQb);

        fHistQaQc = new TH1D("Flow_QaQc_SP","Flow_QaQc_SP",20000,-1000.,1000.);
        fHistQaQc->SetYTitle("dN/dQaQc");
        fHistQaQc->SetXTitle("dQaQc");
        fHistQaQc->StatOverflows(kTRUE);
        tQARelated->Add(fHistQaQc);

        fHistQbQc = new TH1D("Flow_QbQc_SP","Flow_QbQc_SP",20000,-1000.,1000.);
        fHistQbQc->SetYTitle("dN/dQbQc");
        fHistQbQc->SetXTitle("dQbQc");
        fHistQbQc->StatOverflows(kTRUE);
        tQARelated->Add(fHistQbQc);

        fHistQaQbCos = new TH1D("Flow_QaQbCos_SP","Flow_QaQbCos_SP",63,0.,TMath::Pi());
        fHistQaQbCos->SetYTitle("dN/d#phi");
        fHistQaQbCos->SetXTitle("#phi");
        tQARelated->Add(fHistQaQbCos);

        fHistNumberOfSubtractedDaughters = new TH1I("NumberOfSubtractedDaughtersInQ",";daughters;counts;Number of daughters subtracted from Q",5,0.,5.);
        tQARelated->Add(fHistNumberOfSubtractedDaughters);
    }

    fHistList->Add(uQRelated);
    fHistList->Add(nuaRelated);
    fHistList->Add(errorRelated);
    fHistList->Add(tQARelated);

    TH1::AddDirectory(oldHistAddStatus);
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleSP::Make(AliFlowEventSimple* anEvent) {
    // Scalar Product method
    if (!anEvent) return; // for coverity

    fCentralityWeight = anEvent->GetPsi5();
    // Get Q vectors for the subevents
    AliFlowVector* vQarray = new AliFlowVector[2];
    if (fUsePhiWeights)
        anEvent->Get2Qsub(vQarray,fHarmonic,fWeightsList,kTRUE);
    else
        anEvent->Get2Qsub(vQarray,fHarmonic);
    // Subevent a
    AliFlowVector vQa = vQarray[0];
    // Subevent b
    AliFlowVector vQb = vQarray[1];
    delete [] vQarray;

    // multiplicity here corresponds to the V0 equalized multiplicity
    Double_t dMa = vQa.GetMult();
    if( dMa < 2 ) return;
    Double_t dMb = vQb.GetMult();
    if( dMb < 2 ) return;
    //fill control histograms
    fCommonHists->FillControlHistograms(anEvent);

    //Normalizing: weight the Q vectors for the subevents
    Double_t dNa = dMa;//fNormalizationType ? dMa: vQa.Mod(); // SP corresponds to true
    Double_t dNb = dMb;//fNormalizationType ? dMb: vQb.Mod(); // SP corresponds to true
    Double_t dWa = 1.;//fNormalizationType ? dMa: 1; // SP corresponds to true
    Double_t dWb = 1.;//NormalizationType ? dMb: 1; // SP corresponds to true

    // in the original SP method, event weights correspond to multiplicity
    // for the V0 this does not necessarily make sense as the eq mult
    // does not scale directly with charged track mult / centrality

    // this task does not support EP style running, as this is 
    // ambiguous at best

    //Scalar product of the two subevents vectors
    Double_t dQaQb = (vQa*vQb);

    //      01    10     11   <===   fTotalQVector
    // Q ?= Qa or Qb or QaQb
    AliFlowVector vQm;
    vQm.Set(0.0,0.0);
    AliFlowVector vQTPC;
    vQTPC.Set(anEvent->GetPsi2(), anEvent->GetPsi3());

    Double_t dNq=0;
    if( (fTotalQvector%2)>0 ) { // 01 or 11
        vQm += vQa;
        dNq += dMa;
    }
    if( fTotalQvector>1 ) { // 10 or 11
        vQm += vQb;
        dNq += dMb;
    }
    Double_t dWq = 1.;//fNormalizationType ? dNq: 1; // SP corresponds to true


    // if we dont normalize the q-vectors to the multiplicity, just use 1 here
    if(!fScaling) {
        dNa = 1.;
        dNb = 1.;
    }

    // this profile stores the scalar product of the two 
    // subevent q-vectors (the denominator or 'resolution' of the
    // measurement)
    // fHistProQaQbNorm->Fill(1., dQaQb/dNa/dNb); 

    //loop over the tracks of the event
    AliFlowTrackSimple*   pTrack = NULL; 
    Int_t iNumberOfTracks = anEvent->NumberOfTracks(); 
    Double_t dMq = 0;
    for (Int_t i=0;i<iNumberOfTracks;i++) {
        // so this is a track loop ...
        pTrack = anEvent->GetTrack(i) ; 
        if (!pTrack) continue;
        Double_t dPhi = pTrack->Phi();
        Double_t dPt  = pTrack->Pt();
        Double_t dEta = pTrack->Eta();

        // make sure it's not a vzero track
        if(TMath::Abs(dEta) > 0.9) continue;

        //calculate vU
        // this is the track q-vecotr (small u)
        TVector2 vU;
        Double_t dUX = TMath::Cos(fHarmonic*dPhi);
        Double_t dUY = TMath::Sin(fHarmonic*dPhi);
        vU.Set(dUX,dUY);

        //      01    10     11   <===   fTotalQVector
        // Q ?= Qa or Qb or QaQb
        // this will be the vector for hte scalar product itself
        vQm.Set(0.0,0.0); //start the loop fresh
        dMq=0;
        if( (fTotalQvector%2)>0 ) { // 01 or 11
            vQm += vQa;
            dMq += dMa;
        }
        if( fTotalQvector>1 ) { // 10 or 11
            vQm += vQb;
            dMq += dMb;
        }

        dNq = dMq;//fNormalizationType ? dMq : vQm.Mod();
        dWq = 1;//fNormalizationType ? dMq : 1;

        // this little guy is the enumerator of the final equation
        Double_t dUQ = vU*vQm;
        // if we dont want scaling, disable it here
        if(!fScaling) dNq = 1;
        //fill the profile histograms
        for(Int_t iPOI=0; iPOI!=2; ++iPOI) {
            if( (iPOI==0)&&(!pTrack->InRPSelection()) )
                continue;
            if( (iPOI==1)&&(!pTrack->InPOISelection(fPOItype)) )
                continue;
            fHistProUQ[iPOI][0]->Fill(dPt ,fCentralityWeight*dUQ/dNq); //Fill (uQ/Nq') with weight (Nq')
            fHistProUQ[iPOI][1]->Fill(dEta,fCentralityWeight*dUQ/dNq); //Fill (uQ/Nq') with weight (Nq')
            //fHistProUQQaQb[iPOI][0]-> Fill(dPt,(dUQ*dUQ/dNq)/(dQaQb/dNa/dNb)); //Fill [Qu/Nq']*[QaQb/NaNb] with weight (Nq')NaNb
        }
    }//loop over tracks
    //Filling QA (for compatibility with previous version)
    if(!fMinimalBook) {
        fHistQaQb->Fill(vQa*vQb/dNa/dNb, fCentralityWeight);

        fHistQaQc->Fill(vQa*vQTPC, fCentralityWeight);
        fHistQbQc->Fill(vQb*vQTPC, fCentralityWeight);

    }


}

//--------------------------------------------------------------------  
void AliFlowAnalysisWithSimpleSP::GetOutputHistograms(TList *outputListHistos){
    //get pointers to all output histograms (called before Finish())
    fHistList = outputListHistos;

    fCommonHists = (AliFlowCommonHist*) fHistList->FindObject("AliFlowCommonHist_SP");
    //  fCommonHistsuQ = (AliFlowCommonHist*) fHistList->FindObject("AliFlowCommonHist_uQ");
    fCommonHistsRes = (AliFlowCommonHistResults*) fHistList->FindObject("AliFlowCommonHistResults_SP");
    fHistProConfig = (TProfile*) fHistList->FindObject("FlowPro_Flags_SP");
    if(!fHistProConfig) printf("Error loading fHistProConfig\n");
    TList *uQ = (TList*) fHistList->FindObject("uQ");
    TList *nua = (TList*) fHistList->FindObject("NUA");
    TList *error = (TList*) fHistList->FindObject("error");
    TList* tQARelated = (TList*) fHistList->FindObject("QA");

    fHistProQaQbNorm = (TProfile*) error->FindObject("FlowPro_QaQbNorm_SP");
    if(!fHistProQaQbNorm) printf("Error loading fHistProQaQbNorm\n");
    fHistProNUAq = (TProfile*) nua->FindObject("FlowPro_NUAq_SP");
    if(!fHistProNUAq) printf("Error loading fHistProNUAq\n");
    fHistSumOfWeights = (TH1D*) error->FindObject("Flow_SumOfWeights_SP");
    if(!fHistSumOfWeights) printf("Error loading fHistSumOfWeights\n");

    TString sPOI[2] = {"RP","POI"};
    TString sEta[2] = {"Pt","eta"};
    TString sWeights[3] = {"uQ","uQuQ","uQQaQb"};
    for(Int_t iPOI=0; iPOI!=2; ++iPOI) for(Int_t iSpace=0; iSpace!=2; ++iSpace) {
        fHistProUQ[iPOI][iSpace] = (TProfile*) uQ->FindObject( Form( "FlowPro_UQ_%s%s_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ) );
        if(!fHistProUQ[iPOI][iSpace]) printf("Error loading fHistProUQ[%d][%d]\n",iPOI,iSpace);
        fHistProNUAu[iPOI][iSpace][0] = (TProfile*) nua->FindObject( Form("FlowPro_NUAu_%s%s_IM_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ) );
        if(!fHistProNUAu[iPOI][iSpace][0]) printf("Error loading fHistProNUAu[%d][%d][0]\n",iPOI,iSpace);
        fHistProNUAu[iPOI][iSpace][1] = (TProfile*) nua->FindObject( Form("FlowPro_NUAu_%s%s_RE_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ) );
        if(!fHistProNUAu[iPOI][iSpace][1]) printf("Error loading fHistProNUAu[%d][%d][1]\n",iPOI,iSpace);
        fHistProUQQaQb[iPOI][iSpace] = (TProfile*) error->FindObject( Form("FlowPro_UQQaQb_%s%s_SP", sEta[iSpace].Data(), sPOI[iPOI].Data() ) );
        for(Int_t i=0; i!=3; ++i){
            fHistSumOfWeightsu[iPOI][iSpace][i] = (TH1D*) error->FindObject( Form("Flow_SumOfWeights_%s%s_%s_SP",sWeights[i].Data(),sPOI[iPOI].Data(),sEta[iSpace].Data()) );
            if(!fHistSumOfWeightsu[iPOI][iSpace][i]) printf("Error loading fHistSumOfWeightsu[%d][%d][%d]\n",iPOI,iSpace,i);
        }
    }
    if(fHistProConfig) {
        fApplyCorrectionForNUA = (Int_t) fHistProConfig->GetBinContent(1);
        fNormalizationType  = (Int_t) fHistProConfig->GetBinContent(2);
        fUsePhiWeights = (Int_t) fHistProConfig->GetBinContent(3);
        fHarmonic = (Int_t) fHistProConfig->GetBinContent(4);
    }
    fHistQaQb = (TH1D*)tQARelated->FindObject("Flow_QaQb_SP");
    fHistQaQc = (TH1D*)tQARelated->FindObject("Flow_QaQc_SP");
    fHistQbQc = (TH1D*)tQARelated->FindObject("Flow_QbQc_SP");


}   



//--------------------------------------------------------------------            
void AliFlowAnalysisWithSimpleSP::Finish(Bool_t A) {
    //calculate flow and fill the AliFlowCommonHistResults
    printf("AliFlowAnalysisWithSimpleSP::Finish()\n");

    // access harmonic:
    fApplyCorrectionForNUA = 0;//(Int_t)(fHistProConfig->GetBinContent(1));
    fNormalizationType = 1;//(Int_t)(fHistProConfig->GetBinContent(2));
    fHarmonic = (Int_t)(fHistProConfig->GetBinContent(4));

    printf("*************************************\n");
    printf("*************************************\n");
    printf("      Integrated flow from  SIMPLE         \n");
    printf("         Scalar Product              \n\n");

    //Calculate reference flow
    //----------------------------------
    //weighted average over (QaQb/NaNb) with weight (NaNb)


    if(!fHistQaQb) {
        printf(" PANIC: run with full booking !");
        return;
    }
    Double_t dEntriesQaQb = fHistQaQb->GetEntries();
    if( dEntriesQaQb < 1 ) {
        printf(" fHistQaQb has less than 1 entry, probably sub-events are misconfigured \n");
        return;
    }
    //fHistQaQb->GetXaxis()->SetRangeUser(-10.,10.);
    //fHistQaQC->GetXaxis()->SetRangeUser(-1000., 1000.);
    //fHistQbQc->GetXaxis()->SetRangeUser(-1000., 1000.);
    
    Double_t dQaQb = fHistQaQb->GetMean();
    Double_t dQaQc = fHistQaQc->GetMean();
    Double_t dQbQc = fHistQbQc->GetMean();

    if(A) {
      // now it's the resolution of A (which is vzero C) 
      dQaQb = (dQaQb*dQaQc)/dQbQc;
    } else {
      // else we select the resolution of B (which is V0A)
      dQaQb = (dQaQb*dQbQc)/dQaQc;
    }

    Double_t dSpreadQaQb = fHistQaQb->GetRMS();
    // this is the `resolution'
    if(dQaQb <= .0 ) {
        printf(" Panic! the average of QaQb <= 0! Probably you need to run on more events !\n");
        printf("  \tusing dummy value 1 to avoid segfault \n");
        dQaQb = 1.;
        dSpreadQaQb = 1.;
    } 
    Double_t dV = TMath::Sqrt(dQaQb);

    printf("ResSub (sqrt of scalar product of sub-event qvectors) = %f\n", dV );
    printf("fTotalQvector %d \n",fTotalQvector);
    // we just take the spread as the uncertainty

    Double_t dStatErrorQaQb = dSpreadQaQb;

    Double_t dVerr = 0.;
    if(dQaQb > 0.) dVerr = (1./(2.*pow(dQaQb,0.5)))*dStatErrorQaQb;
    fCommonHistsRes->FillIntegratedFlow(dV,dVerr);
    printf("v%d(subevents) = %f +- %f\n",fHarmonic,dV,dVerr);

    Int_t iNbins[2];
    iNbins[0] = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
    iNbins[1] = AliFlowCommonConstants::GetMaster()->GetNbinsEta();

    //Calculate differential flow and integrated flow (RP, POI)
    //---------------------------------------------------------
    //v as a function of eta for RP selection
    for(Int_t iRFPorPOI=0; iRFPorPOI != 2; ++iRFPorPOI)
        for(Int_t iPTorETA=0; iPTorETA != 2; ++iPTorETA)
            for(Int_t b=1; b != iNbins[iPTorETA]+1; ++b) {
                Double_t duQpro = fHistProUQ[iRFPorPOI][iPTorETA]->GetBinContent(b);
                Double_t dv2pro = -999.;
                if( TMath::Abs(dV!=0.) && fExternalResolution < 0 ) { dv2pro = duQpro/dV; }
                else if(fExternalResolution > 0) { dv2pro = duQpro / fExternalResolution; }
                Double_t dv2ProErr = fHistProUQ[iRFPorPOI][iPTorETA]->GetBinError(b);
                if(fHistProUQ[iRFPorPOI][iPTorETA]->GetBinEntries(b) > 1) dv2ProErr/=TMath::Sqrt(fHistProUQ[iRFPorPOI][iPTorETA]->GetBinEntries(b)); 
                if(fExternalResErr > 0) { 
                    // do default error propagation for a fraction where
                    Double_t a(0), b(0), t(0);
                    if(dv2ProErr > 0) a = duQpro*duQpro/(dv2ProErr*dv2ProErr);      // squared relative error
                    if(fExternalResErr > 0) b = fExternalResolution*fExternalResolution/(fExternalResErr*fExternalResErr);
                    t = duQpro*fExternalResolution/(dv2ProErr*fExternalResErr);
                    t = dv2pro*dv2pro*(a + b -2.*t);
                    if(t>0) dv2ProErr =  TMath::Sqrt(t);
                }

                if( (iRFPorPOI==0)&&(iPTorETA==0) )
                    fCommonHistsRes->FillDifferentialFlowPtRP(  b, dv2pro, dv2ProErr);   
                if( (iRFPorPOI==0)&&(iPTorETA==1) )
                    fCommonHistsRes->FillDifferentialFlowEtaRP( b, dv2pro, dv2ProErr);   
                if( (iRFPorPOI==1)&&(iPTorETA==0) ) 
                    fCommonHistsRes->FillDifferentialFlowPtPOI( b, dv2pro, dv2ProErr); 
                if( (iRFPorPOI==1)&&(iPTorETA==1) )
                    fCommonHistsRes->FillDifferentialFlowEtaPOI(b, dv2pro, dv2ProErr);
            }

    printf("\n");
    printf("*************************************\n");
    printf("*************************************\n");
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleSP::WriteHistograms(TDirectoryFile *outputFileName) const {
    //store the final results in output .root file
    outputFileName->Add(fHistList);
    outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}

//--------------------------------------------------------------------            
