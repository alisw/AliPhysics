#include <TChain.h>
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TFile.h"
#include "TF1.h"
#include "TAxis.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

//#include "AliESDEvent.h"
//#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliGenHijingEventHeader.h"

#include "AliCentrality.h"
#include "AliAODVertex.h"

//#include "AliAODVZERO.h"

#include "AliForwardBackwardAnalysis.h"


#include <iostream>
using namespace std;

const double ptMin[N_PT_BINS] = { 0.2, 0.2, 0.8 };
const double ptMax[N_PT_BINS] = { 2.0, 5.0, 5.0 };

const double cWidth[FB_N_CENTR_BINNINGS] = { 10, 5, 2 };

// analysis task


ClassImp(AliForwardBackwardAnalysis)


//________________________________________________________________________
AliForwardBackwardAnalysis::AliForwardBackwardAnalysis(const char *name)
    : AliAnalysisTaskSE(name)
    //    , fAOD(0)
    , fOutputList(0)

    , fCentralityEstimator("V0M")
    , fAODTrackCutBit(128)
    , fRunNumber(0)
    , fNumberOfTracks(0)
    , fVertexX(0)
    , fVertexY(0)
    , fVertexZ(0)

    , fHistVertexZ(0)
    , fHistVertexXY(0)
    , fHistCentrality(0)
    , fHistTrackPt(0)
    , fHistTrackEta(0)
    , fHistTrackPhi(0)

    , fHistLog(0)
    //    , fEventTree(0)
    , fHist2D_ESTIMATOR_VS_multTPC(0)
    , fHist2D_ESTIMATOR_VS_multTPC_afterCuts(0)

    , fBorderToCutOutliersLower(0)
    , fBorderToCutOutliersUpper(0)


{
    // Constructor

    DefineOutput(1, TList::Class());
    //    DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
void AliForwardBackwardAnalysis::UserCreateOutputObjects()
{
    // Create histograms  and add them to the respective lists
    // Called once

    fOutputList = new TList();

    fOutputList->SetOwner(kTRUE);

    fHistCentrality = new	TH1F("fHistCentrality", "Centrality of Events; Centrality Percentile (%); N",50,0.0,100.0);
    fHistVertexZ = new TH1F("fHistVertexZ","Z Coord. of primary vertex;v_{z}, cm;n events",80,-20.0,20.0);
    fHistVertexXY = new TH2F("fHistVertexXY","XY Coord. of primary vertex;v_{x}, cm;v_{y}, cm;n events",100,-0.3,0.3,100,-0.3,0.3);

    fHistTrackPt = new TH1F("fHistTrackPt", "PT distribution;p_{T} (GeV/c);dN/dp_{T}", 200, 0., 5);
    fHistTrackEta = new TH1F("fHistTrackEta", "Eta distribution of Tracks;#eta;dN/d#eta", 160, -2.0, 2.0);
    fHistTrackPhi = new TH1F("fHistTrackPhi", "Phi distribution of Tracks;#eta;dN/d#eta", 160, 0, TMath::TwoPi() );

    fHistLog = new TH1F("fHistLog","Log Variables",100, -0.5, 99.5);

    fOutputList->Add(fHistVertexXY);
    fOutputList->Add(fHistVertexZ);
    fOutputList->Add(fHistCentrality);
    fOutputList->Add(fHistTrackPt);
    fOutputList->Add(fHistTrackEta);
    fOutputList->Add(fHistTrackPhi);
    fOutputList->Add(fHistLog);

    fHist2D_ESTIMATOR_VS_multTPC = new TH2D( "hist2D_ESTIMATOR_VS_multTPC", "hist2D_ESTIMATOR_VS_multTPC;estimator;mult in TPC", 4080, -2, 100, 301, -0.5, 3000.5);
    fOutputList->Add(fHist2D_ESTIMATOR_VS_multTPC);
    fHist2D_ESTIMATOR_VS_multTPC_afterCuts = new TH2D( "fHist2D_ESTIMATOR_VS_multTPC_afterCuts", "fHist2D_ESTIMATOR_VS_multTPC_afterCuts;estimator;mult in TPC", 4080, -2, 100, 301, -0.5, 3000.5);
    fOutputList->Add(fHist2D_ESTIMATOR_VS_multTPC_afterCuts);

    // Adding a Tree



    if(0)for( Int_t i=0; i < FB_N_PHI_BINS; i++ )
        for( Int_t j=0; j < FB_N_PHI_BINS; j++ )
        {
            TString strName = Form( "PtPt_phiBin_%d_%d", i, j );
            fHistPtPt[i][j] = new TH2D( strName, "PtPt;F;B", 320, 0.8, 4.0, 320, 0.8, 4.0 );
            fOutputList->Add( fHistPtPt[i][j] );

            strName = Form( "NN_phiBin_%d_%d", i, j );
            fHistNN[i][j] = new TH2D( strName, "NN;F;B", 51, -0.5, 50.5, 51, -0.5, 50.5 );
            fOutputList->Add( fHistNN[i][j] );
        }

    // ### for FB analysis:
    // const int nFB_etaPhi_bins = 7;
    //    for( Int_t pt=0; pt < N_PT_BINS; pt++ )
    //    {
    //        fHistFB_eachEvent[pt] = new TH2D( "hist_FB_routine_each_event", "FB_routine_each_event;;#eta#varphi_{ID}",
    //                                          4, -0.5, 3.5, nFB_etaPhi_bins, -0.5, nFB_etaPhi_bins-0.5 );
    //    }


    // tune FB pairs:
    for( Int_t pt=0; pt < N_PT_BINS; pt++ )
    {
        cout << "#### pT bin " << pt << endl;
        for( Int_t ep=0; ep < FB_N_ETAPHI_BINNINGS; ep++ )
        {
            WinPairInfo *wp = &winPairs[pt][ep];

            // pt bins
            wp->setPtRanges( ptMin[pt], ptMax[pt], ptMin[pt], ptMax[pt] );

            // eta-phi bins
            wp->setPhiRanges( 0, TMath::TwoPi(), 0, TMath::TwoPi() );

            if ( ep < 7 ) // eW 0.4
                wp->setEtaRanges( 0.6-ep*0.1, 0.8-ep*0.1, -0.8+ep*0.1, -0.6+ep*0.1 );
            else if ( ep < 12 ) // eW 0.2
            {
                int ep2 = ep-7;
                wp->setEtaRanges( 0.4-ep2*0.1, 0.8-ep2*0.1, -0.8+ep2*0.1, -0.4+ep2*0.1 );
            }
            else
            {
                int ep2 = ep-12;

                const int nPhiBins = 32;
                if ( ep2 < nPhiBins*nPhiBins )
                    wp->setEtaRanges( 0.4, 0.8, -0.8, -0.4 );
                else
                    wp->setEtaRanges( -0.8, 0.8, -0.8, 0.8 );

                const double phiStep = TMath::TwoPi()/nPhiBins;
                //                for( Int_t phi1=0; phi1 < nPhiBins; phi1++ )
                //                    for( Int_t phi2=0; phi2 < nPhiBins; phi2++ )
                Int_t phi1 = ep2/nPhiBins;
                if ( ep2 >= nPhiBins*nPhiBins )
                    phi1 -= nPhiBins;
                Int_t phi2 = ep2%nPhiBins;

                wp->setPhiRanges( phi1*phiStep, (phi1+1)*phiStep, phi2*phiStep, (phi2+1)*phiStep );
                //if(ep2>1022)
//                {
//                    cout << "phi: " << phi1 << "..." << phi2 << endl;
//                    int a;
//                    cin >> a;
//                }
            }


            if(1)cout << "pt: " << wp->pt[0] << "..." << wp->pt[1] << ", " << wp->pt[2] << "..." << wp->pt[3]
                      << ", eta: " << wp->eta[0] << "..." << wp->eta[1] << ", " << wp->eta[2] << "..." << wp->eta[3]
                      << ", phi: " << wp->phi[0] << "..." << wp->phi[1] << ", " << wp->phi[2] << "..." << wp->phi[3]
                      << endl;
        }
    }

    // tune 3D histos ingredients-etaPhiWins-centralities
    const int nFBlabels = 31; //28; //24; //16;
    TString strFB_binLabels[nFBlabels] =
    {
        "nEvents",
        "nEventsF",
        "nEventsB",
        "nEventsFB",

        "nF",
        "nB",
        "nFnB",
        "nF2",
        "nB2",

        "ptF",
        "ptB",
        "ptFptB",
        "ptF2",
        "ptB2",

        "ptFnB",
        "ptBnF",

        // ### for C when nPairs OUTSIDE sum:
        "pipjF",
        "pipjB",
        "(nF-1)*sum_pF",
        "(nB-1)*sum_pB",
        "nF*(nF-1)",
        "nB*(nB-1)",

        // ### for C when nPairs INSIDE sum: (like in GOOD_Ck_definition_STAR_2005_0504031.pdf)
        "pipjF_avPerEv",
        "pipjB_avPerEv",
        "(nF-1)*sum_pF_avPerEv",
        "(nB-1)*sum_pB_avPerEv",

        // ### for dptdpt
        "piFpjB",
        "nF*sum_pB",
        "nB*sum_pF",

        // to calc mean pT in F and B for all tracks in all events:
        "sumPtAllEvF",
        "sumPtAllEvB",
    };
    for( Int_t pt=0; pt < N_PT_BINS; pt++ )
    {
        for( Int_t cW=0; cW < FB_N_CENTR_BINNINGS; cW++ )
        {
            TString strName = Form( "fHistFB_ingredients_etaPhi_ptBin%d_cW%d", pt, cW );
            TString strTitle = Form( "%s;FB win ID;FB label;centrality", strName.Data() );
            fHistFB_ingredients_etaPhi_centr[pt][cW]
                    = new TH3D( strName, strTitle
                                , FB_N_ETAPHI_BINNINGS, -0.5, FB_N_ETAPHI_BINNINGS-0.5 // eta-phi window pairs
                                , nFBlabels, -0.5, nFBlabels-0.5 // ingredients
                                , 80/cWidth[cW], 0, 80  // centrality
                                );
            for( Int_t i=0; i < nFBlabels; i++ )
                fHistFB_ingredients_etaPhi_centr[pt][cW]->GetYaxis()->SetBinLabel(i+1, strFB_binLabels[i]);

            fOutputList->Add( fHistFB_ingredients_etaPhi_centr[pt][cW] );
        }
    }



    // eta-phi binning:
    int arrNEtaBins[] = {20, };//   4, 2,      4, 2,    2    };
    int arrNPhiBins[] = {72, };//   64, 64,    32, 32,  49   };

    // ### for DptDpt:
    for( Int_t c=0; c < DPT_N_CENTR_BINS; c++ )
    {
        for( Int_t pt=0; pt < N_PT_BINS; pt++ )
        {
            for( Int_t ep=0; ep < DPT_N_ETAPHI_BINNINGS; ep++ )
            {
                TString strName = Form( "fHistDptDpt_p1p2_etaPhiBinningId%d_ptBin%d_c%d", ep, pt, c );
                TString strTitle = Form( "%s;#Delta#eta;#Delta#varphi", strName.Data() );
                fHistDptDpt_p1p2[pt][c][ep]    = new TH2D( strName, strTitle, arrNEtaBins[ep], 0, 1.6, arrNPhiBins[ep], 0, TMath::TwoPi() );


                strName = Form( "fHistDptDpt_p1_etaPhiBinningId%d_ptBin%d_c%d", ep, pt, c );
                strTitle = Form( "%s;#Delta#eta;#Delta#varphi", strName.Data() );
                fHistDptDpt_p1[pt][c][ep]      = new TH2D( strName, strTitle, arrNEtaBins[ep], 0, 1.6, arrNPhiBins[ep], 0, TMath::TwoPi() );

                strName = Form( "fHistDptDpt_p2_etaPhiBinningId%d_ptBin%d_c%d", ep, pt, c );
                strTitle = Form( "%s;#Delta#eta;#Delta#varphi", strName.Data() );
                fHistDptDpt_p2[pt][c][ep]      = new TH2D( strName, strTitle, arrNEtaBins[ep], 0, 1.6, arrNPhiBins[ep], 0, TMath::TwoPi() );

                strName = Form( "fHistDptDpt_n1n2_etaPhiBinningId%d_ptBin%d_c%d", ep, pt, c );
                strTitle = Form( "%s;#Delta#eta;#Delta#varphi", strName.Data() );
                fHistDptDpt_n1n2[pt][c][ep]    = new TH2D( strName, strTitle, arrNEtaBins[ep], 0, 1.6, arrNPhiBins[ep], 0, TMath::TwoPi() );

                fOutputList->Add( fHistDptDpt_p1p2[pt][c][ep] );
                fOutputList->Add( fHistDptDpt_p1[pt][c][ep]   );
                fOutputList->Add( fHistDptDpt_p2[pt][c][ep]   );
                fOutputList->Add( fHistDptDpt_n1n2[pt][c][ep] );

            }
            // additional histos to calculate meanPt
            TString strName = Form( "fHistDptDpt_forMeanPt_ptBin%d_c%d", pt, c );
            fHistDptDpt_meanPt[pt][c] = new TH1D( strName, strName, 1200, 0, 6.);
            fOutputList->Add( fHistDptDpt_meanPt[pt][c] );

        }
    }


    fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","[0]+[1]*exp([2]*x)",0,90);
    fBorderToCutOutliersLower->SetParameters( -100, 1500, -0.04 );
    //    fBorderToCutOutliersLower->DrawClone("same");


    fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","[0]+[1]*exp([2]*x)",0,90);
    fBorderToCutOutliersUpper->SetParameters( 50, 2500, -0.04 );
    //    fBorderToCutOutliersUpper->DrawClone("same");



    // NEW HISTO added to fOutputList here
    PostData(1, fOutputList); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliForwardBackwardAnalysis::UserExec(Option_t *)
{
    // Main loop
    // Called for each event

    const double cBinMin_dptdpt[] = { 0, };// 5, 10, 20,    30, 40, 50, 60,     70, 0 };
    const double cBinMax_dptdpt[] = { 5, };// 10, 20, 30,   40, 50, 60, 70,     80, 2 };

    const int myFilterBit = fAODTrackCutBit;//128; //768

    // Post output data.
    AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD) {
        printf("ERROR: fAOD not available\n");
        fHistLog->Fill(0);
        return;
    }

    fHistLog->Fill(1);

    fRunNumber = fAOD->GetRunNumber();

    // ##### event centrality
    AliCentrality* centrality = fAOD->GetCentrality();
    fCentPercentile[0] = (Float_t)centrality->GetCentralityPercentile("V0M");
    fCentPercentile[1] = (Float_t)centrality->GetCentralityPercentile("TRK");
    // This is either set upto <~44% or is 100. i.e. cent is set to 100 if 44%<cent<100:
    fCentPercentile[2] = (Float_t)centrality->GetCentralityPercentile("ZEMvsZDC");

    if(fCentPercentile[0]<0.0) {
        fHistLog->Fill(3);
        return;
    }
    double centr = (Float_t)centrality->GetCentralityPercentile( fCentralityEstimator.Data() );
    //  cout<<" Got Event of Centrality by V0M = "<<fCentPercentile[0]<<endl;

    // ##### choose a particular trigger
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    if(isSelected==kFALSE) {
        //cout<<" Not selected by trigger"<<endl;
        fHistLog->Fill(13); // rejected by trigger condition
        return;
    }

    // ##### obtain vertex and apply vertex cuts
    const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
    if (!vtx)
    {
        fHistLog->Fill(9); // rejected due to no vertex
        return;
    }
    Double32_t fCov[6];
    vtx->GetCovarianceMatrix(fCov);
    if(vtx->GetNContributors() <= 0) { fHistLog->Fill(9); return;}
    if(fCov[5] == 0) {fHistLog->Fill(9); return;}

    fVertexX= (Float_t)vtx->GetX();
    fVertexY= (Float_t)vtx->GetY();
    fVertexZ= (Float_t)vtx->GetZ();
    //  cout<<" Vertex XYZ are="<<fVertexX<<" "<<fVertexY<<" "<<fVertexZ<<endl;

    if ( fabs(fVertexZ) > 8 )
        return;

    fHistVertexZ->Fill(fVertexZ);
    fHistVertexXY->Fill(fVertexX,fVertexY);

    fHistLog->Fill(50);

    //    if ( centr > 2 )
    //        return;

    //  cout<<"Centrality Percentile is "<<fCentPercentile[0]<<" "<<fCentPercentile[1]<<" "<<fCentPercentile[2]<<endl;
    fHistCentrality->Fill(centr);

    // ##### count tracks and fill multTPC vs centrality hist2D
    Int_t nAllTracks = fAOD->GetNumberOfTracks();
    Int_t multTPC = 0;
    for (Int_t i = 0; i < nAllTracks; i++)
    {
        AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(i);
        if (!track) {
            printf("ERROR: Could not receive track %d\n", i);
            continue;
        }
        if(!track->TestFilterBit(768)) continue;
        if ( fabs(track->Eta()) > 0.8 )
            continue;
        if ( fabs(track->Pt()) < 0.2 || fabs(track->Pt()) > 2.0 )
            continue;
        multTPC++;
    }
    fHist2D_ESTIMATOR_VS_multTPC->Fill(centr, multTPC);
    // !!! CUT ON MULT_TPC vs centr HIST!
    if ( multTPC < fBorderToCutOutliersLower->Eval(centr)
         || multTPC > fBorderToCutOutliersUpper->Eval(centr) )
        return;
    fHist2D_ESTIMATOR_VS_multTPC_afterCuts->Fill(centr, multTPC);



    // ##### for FB 2D analysis:
    int winInfo_n[FB_N_PHI_BINS];
    double winInfo_pt[FB_N_PHI_BINS];
    for ( int phiBin = 0; phiBin < FB_N_PHI_BINS; phiBin++ )
    {
        winInfo_n[phiBin] = 0;
        winInfo_pt[phiBin] = 0;
    }

    // ##### for FB analysis:
    for( Int_t pt=0; pt < N_PT_BINS; pt++ )
        for( Int_t ep=0; ep < FB_N_ETAPHI_BINNINGS; ep++ )
            winPairs[pt][ep].reset();




    double phiStep = TMath::TwoPi()/FB_N_PHI_BINS;

    // Track loop
    //  printf("There are %d tracks in this event\n", fAOD->GetNumberOfTracks());
    Int_t NTracks = 0; // This is the accepted track number
    //    Int_t nAllTracks = (fAOD->GetTracks())->GetEntriesFast();
    for (Int_t i = 0; i < nAllTracks; i++)
    {
        AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(i);
        if (!track) {
            printf("ERROR: Could not receive track %d\n", i);
            continue;
        }

        if( !track->TestFilterBit( myFilterBit ) )
            continue;

        if ( fabs( track->Eta() ) > 1 )
            continue;


        NTracks++;

        double pt  = (Float_t)track->Pt();
        double phi = (Float_t)track->Phi();
        double eta = (Float_t)track->Eta();
        int charge = (Int_t)track->Charge();

        fHistTrackPt->Fill(pt);
        fHistTrackEta->Fill(eta);
        fHistTrackPhi->Fill(phi);

        // ##### for FB 2D analysis:
        if ( pt > 0.8 && pt < 10.0 )
        {
            if ( eta > -0.8 && eta < 0.8 )
            {
                for ( int phiBin = 0; phiBin < FB_N_PHI_BINS; phiBin++ )
                {
                    if ( phi >= phiStep*phiBin && phi < phiStep*(phiBin+1) )
                    {
                        winInfo_n[phiBin]++;
                        winInfo_pt[phiBin] += pt;
                    }
                }
            }
        }


        // ##### for FB analysis:
        for( Int_t ptBin=0; ptBin < N_PT_BINS; ptBin++ )
            for( Int_t ep=0; ep < FB_N_ETAPHI_BINNINGS; ep++ )
                winPairs[ptBin][ep].addTrack( pt, eta, phi, charge, 0 );


        // ###### DptDpt analysis
        if ( fabs(eta) < 0.8
             && ( pt > 0.2 && pt < 5.0 )
             )
        {
            // fill mean pT in each pt- and centr-bin
            for( Int_t cBin=0; cBin < DPT_N_CENTR_BINS; cBin++ )
            {
                if ( centr < cBinMin_dptdpt[cBin] || centr > cBinMax_dptdpt[cBin] )
                    continue;
                for( Int_t ptBin=0; ptBin < N_PT_BINS; ptBin++ )
                {
                    if ( pt > ptMin[ptBin] && pt < ptMax[ptBin] )
                        fHistDptDpt_meanPt[ptBin][cBin]->Fill(pt);
                }
            }

            // inner loop over tracks:
            for (Int_t j = i+1; j < nAllTracks; j++)
            {
                AliAODTrack* track2 = (AliAODTrack*)fAOD->GetTrack(j);
                if (!track2) {
                    printf("ERROR: Could not receive track %d\n", j);
                    continue;
                }
                if( !track2->TestFilterBit( myFilterBit ) )
                    continue;


                double pt2  = (Float_t)track2->Pt();
                double phi2 = (Float_t)track2->Phi();
                double eta2 = (Float_t)track2->Eta();

                if ( fabs(eta2) > 0.8 )
                    continue;
                if ( pt2 < 0.2 || pt2 > 5.0 )
                    continue;

                double dEta = fabs(eta2-eta);
                double dPhi = phi2-phi;
                if ( dPhi < 0 )
                    dPhi += TMath::TwoPi();

                for( Int_t cBin=0; cBin < DPT_N_CENTR_BINS; cBin++ )
                {
                    if ( centr < cBinMin_dptdpt[cBin] || centr > cBinMax_dptdpt[cBin] )
                        continue;
                    for( Int_t ptBin=0; ptBin < N_PT_BINS; ptBin++ )
                    {
                        if ( (pt > ptMin[ptBin] && pt < ptMax[ptBin] )
                             && (pt2 > ptMin[ptBin] && pt2 < ptMax[ptBin] )
                             )
                        {
                            for( Int_t ep=0; ep < DPT_N_ETAPHI_BINNINGS; ep++ )
                            {
                                fHistDptDpt_p1p2[ptBin][cBin][ep]->Fill( dEta, dPhi, pt*pt2 );
                                fHistDptDpt_p1[ptBin][cBin][ep]->Fill( dEta, dPhi, pt );
                                fHistDptDpt_p2[ptBin][cBin][ep]->Fill( dEta, dPhi, pt2 );
                                fHistDptDpt_n1n2[ptBin][cBin][ep]->Fill( dEta, dPhi, 1 );
                            }
                        }
                    } // end of ptBins loop
                } // end of centrBins loop
            } //end of track2 loop
        } // end of if for eta and pt

    } //end of track loop
    fNumberOfTracks = NTracks;

    //final actions wins for FB 2D:
    for ( int phiBin = 0; phiBin < FB_N_PHI_BINS; phiBin++ )
    {
        if ( winInfo_n[phiBin] > 0 )
            winInfo_pt[phiBin] /= winInfo_n[phiBin];
        else
            winInfo_pt[phiBin] = -1;
    }

    //fill PtPt hist:
    if(0)for ( int binF = 0; binF < FB_N_PHI_BINS; binF++ )
        for ( int binB = 0; binB < FB_N_PHI_BINS; binB++ )
        {
            if ( winInfo_pt[binF]>=0 && winInfo_pt[binB]>=0 )
                fHistPtPt[binF][binB]->Fill( winInfo_pt[binF], winInfo_pt[binB] );
            fHistNN[binF][binB]->Fill( winInfo_n[binF], winInfo_n[binB] );
        }

    // ##### final actions for FB, fill hist:
    for( Int_t pt=0; pt < N_PT_BINS; pt++ )
    {
        for( Int_t ep=0; ep < FB_N_ETAPHI_BINNINGS; ep++ )
        {
            WinPairInfo *wp = &winPairs[pt][ep];
            wp->finalActionsForEvent();

            for( Int_t cW=0; cW < FB_N_CENTR_BINNINGS; cW++ )
            {

                //            if ( centr < 10 && wp->nF>0 && wp->nB>0 && pt==0)
                //            {
                //                cout << "AHNUNG: ptBin=" << pt << ", epBin=" << ep << endl;
                //            }

                // fill FB info into hist bins:
                TH3D *hist = fHistFB_ingredients_etaPhi_centr[pt][cW];

                hist->Fill( ep,   "nEvents",    centr, 1 );
                hist->Fill( ep,   "nF",         centr, wp->nF );
                hist->Fill( ep,   "nB",         centr, wp->nB );
                hist->Fill( ep,   "nF2",        centr, wp->nF*wp->nF );
                hist->Fill( ep,   "nB2",        centr, wp->nB*wp->nB );
                hist->Fill( ep,   "nFnB",       centr, wp->nF*wp->nB );

                hist->Fill( ep,   "sumPtAllEvF",        centr, wp->nF*wp->PtF );
                hist->Fill( ep,   "sumPtAllEvB",        centr, wp->nB*wp->PtB );

                if ( wp->nF > 0 )
                {
                    hist->Fill( ep,   "nEventsF",   centr, 1 );
                    hist->Fill( ep,   "ptF",        centr, wp->PtF );
                    hist->Fill( ep,   "ptF2",       centr, wp->PtF*wp->PtF );
                    hist->Fill( ep,   "ptFnB",      centr, wp->PtF*wp->nB );

                    // for C when av over pairs is OUTSIDE sum:
                    hist->Fill( ep,   "pipjF",              centr, wp->pipjF );
                    hist->Fill( ep,   "(nF-1)*sum_pF",      centr, (wp->nF-1)*(wp->nF*wp->PtF) ); // (n-1)*sumPt
                    hist->Fill( ep,   "nF*(nF-1)",          centr, (wp->nF-1)*wp->nF ); // (n-1)*n

                    // for C when av over pairs is INSIDE sum:
                    int nPairsF = (wp->nF-1)*wp->nF;
                    if ( nPairsF > 0 )
                    {
                        hist->Fill( ep,   "pipjF_avPerEv",              centr, wp->pipjF / nPairsF );
                        hist->Fill( ep,   "(nF-1)*sum_pF_avPerEv",      centr, (wp->nF-1)*(wp->nF*wp->PtF) / nPairsF );
                    }
                }
                if ( wp->nB > 0 )
                {
                    hist->Fill( ep,   "nEventsB",   centr, 1 );
                    hist->Fill( ep,   "ptB",        centr, wp->PtB );
                    hist->Fill( ep,   "ptB2",       centr, wp->PtB*wp->PtB );
                    hist->Fill( ep,   "ptBnF",      centr, wp->PtB*wp->nF );

                    // for C when av over pairs is OUTSIDE sum:
                    hist->Fill( ep,   "pipjB",              centr, wp->pipjB );
                    hist->Fill( ep,   "(nB-1)*sum_pB",      centr, (wp->nB-1)*(wp->nB*wp->PtB) ); // (n-1)*sumPt
                    hist->Fill( ep,   "nB*(nB-1)",          centr, (wp->nB-1)*wp->nB ); // (n-1)*n

                    // for C when av over pairs is INSIDE sum: (like in GOOD_Ck_definition_STAR_2005_0504031.pdf)
                    int nPairsB = (wp->nB-1)*wp->nB;
                    if ( nPairsB > 0 )
                    {
                        hist->Fill( ep,   "pipjB_avPerEv",              centr, wp->pipjB / nPairsB );
                        hist->Fill( ep,   "(nB-1)*sum_pB_avPerEv",      centr, (wp->nB-1)*(wp->nB*wp->PtB) / nPairsB );
                    }
                }
                if ( wp->nF > 0 && wp->nB > 0 )
                {
                    hist->Fill( ep,   "nEventsFB",  centr, 1 );
                    hist->Fill( ep,   "ptFptB",     centr, wp->PtF*wp->PtB );

                    // for dptdpt:
                    hist->Fill( ep,   "piFpjB",         centr, wp->piFpjB );
                    hist->Fill( ep,   "nF*sum_pB",      centr, wp->nF*(wp->nB*wp->PtB) );
                    hist->Fill( ep,   "nB*sum_pF",      centr, wp->nB*(wp->nF*wp->PtF) );
                }
            } // end of centr binnings
        } // end of eta-phi bins
    } // end of pt bins

    //  cout<<" Number of AOD tracks are "<<fNumberOfTracks<<endl;

    PostData(1, fOutputList);
}

//________________________________________________________________________
void AliForwardBackwardAnalysis::Terminate(Option_t *)
{
    //    fHist2D_ESTIMATOR_VS_multTPC->DrawClone("colz");

}


