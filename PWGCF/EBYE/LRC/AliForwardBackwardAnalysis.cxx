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

#include "TRandom3.h"

#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliForwardBackwardAnalysis.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include <AliMultiplicity.h>
#include <AliESDv0.h>

#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliStack.h>
#include <AliESDtrackCuts.h>
#include <AliKalmanTrack.h>
//#include <AliPhysicsSelection.h>
#include <AliCentrality.h>
#include <AliESDZDC.h>
#include <AliESDFMD.h>


#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliGenHijingEventHeader.h"

#include "AliCentrality.h"
#include "AliAODVertex.h"

#include "AliAnalysisUtils.h"
#include "AliAODTracklets.h"
#include "AliMultSelection.h"

#include "AliAODVZERO.h"

#include <AliPIDResponse.h>
#include <AliPIDCombined.h>




#include <iostream>
using namespace std;

//const double ptMin[N_PT_BINS] = { 0.2, 0.2, 0.8 };
//const double ptMax[N_PT_BINS] = { 2.0, 5.0, 5.0 };

//const double ptMin[N_PT_BINS] = { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, /*12*/ 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,  /*7*/ } ;//0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  /*7*/};
//const double ptMax[N_PT_BINS] = { 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, /*12*/ 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, /*7*/ } ;//1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, /*7*/};

const double ptMin[N_PT_BINS] = { 0.2, 0.2, 0.8 };
const double ptMax[N_PT_BINS] = { 2.0, 5.0, 5.0 };


//const double cWidth[FB_N_CENTR_BINNINGS] = { 10, 5, 2 };
//const double cWidth[FB_N_CENTR_BINNINGS] = { 5 }; //2
const double cWidth[FB_N_CENTR_BINNINGS] = { 2 };


//const int nArrEtaPhiBins[FB_N_ETAPHI_BINNINGS] = { 7, 5 };
//const int nArrPhiBins[FB_N_ETAPHI_BINNINGS] = { 1, 1 };

//const int nArrEtaPhiBins[FB_N_ETAPHI_BINNINGS] = { 7, 5,   32, 32,   48, 48 }; //,     64, 64, 72, 72 };
//const int nArrPhiBins[FB_N_ETAPHI_BINNINGS] = { 1, 1,   32, 32,   48, 48 }; //,     64, 64, 72, 72 };

const int nArrEtaPhiBins[FB_N_ETAPHI_BINNINGS] = { 7, 5,   32, 32 };
const int nArrPhiBins[FB_N_ETAPHI_BINNINGS] = { 1, 1,   32, 32    };

//const int nArrEtaPhiBins[FB_N_ETAPHI_BINNINGS] = { 7, 5 , 231 };
//const int nArrPhiBins[FB_N_ETAPHI_BINNINGS] = { 1, 1, 1 };

//const int nEtaBins[FB_N_ETAPHI_BINNINGS] = { 7, 5, 1, 1, 1, 1 };

const char *AliForwardBackwardAnalysis::fgkBinMomDesc[AliForwardBackwardAnalysis::kPtBins] = {
    " 0 <= p < 0.5 GeV/c",
    " 0.5 <= p < 0.7 GeV/c",
    " 0.7 <= p < 1.0 GeV/c",
    " 1.0 <= p < 1.5 GeV/c",
    " 1.5 <= p < 2.0 GeV/c",
    " p >= 2.0 GeV/c"
};


ClassImp(AliForwardBackwardAnalysis)


//________________________________________________________________________
AliForwardBackwardAnalysis::AliForwardBackwardAnalysis(const char *name, Bool_t runKine)
    : AliAnalysisTaskSE(name)
    //    , fAOD(0)
    , fOutputList(0)
    , fRunKine(false)

    , fCentralityEstimator("V0M")
    , fAODTrackCutBit(128)
    , fMinCentrality(0)
    , fMaxCentrality(80)
    , fRunNumber(0)
    , fNumberOfTracks(0)

    , fParameterWhichLHCperiod(0)
    , fHistCentrality(0)

    , fHistEventCutStats(0)
    , fHistVertexStats (0)
    , fHistVertexZ(0)
    , fHistVertexXY(0)
    , fHistVzDiffSPDvsTracks(0)
    , fHistVzDiffTPCvsTracks(0)
    , fHist2DmultTOFvsTracks(0)
    , fHist2DmultTOFvsTracksAfterCuts(0)
    , fHist2DmultV0MvsTracksTPCoutBEFOREcutOnLargeMultEsd(0)
    , fHist2DmultV0MvsTracksTPCout(0)
    , fHistVertexNconributors(0)
    , fHistNumberOfPileupVerticesTracks(0)
    , fHistNumberOfPileupVerticesSPD(0)
    , fHistNumberOfVertices(0)


    , fHistTrackPt(0)
    , fHistTrackEta(0)
    , fHistTrackPhi(0)

    , fHistLog(0)
    , fHistAcceptedTracks(0)

    //    , fEventTree(0)
    , fHist2D_ESTIMATOR_VS_multTPC(0)
    , fHist2D_ESTIMATOR_VS_multTPC_afterCuts(0)

    , fBorderToCutOutliersLower(0)
    , fBorderToCutOutliersUpper(0)

    //    , fHist3D_EffFromOutside(0x0)
    , fAmplifiedInefficiencyForHist3D(-1)

    // /*
    // ########### NEW PID 7.11.2016
    , fPIDResponse(0x0)
    , fPIDCombined(0x0)

    , fHistPtPID()
    , fTPCsignalInPtBins()
    , fSignalDeltaTPCTOF()
    , fnSigmasTPCTOF()

    // from ANALYSIS/ANALYSISalice/AliAnalysisTaskPIDCombined.cxx:
    , fProbTPCnSigma()
    , fProbTOFnSigma()
    , fProbTPCTOFnSigmaTPC()
    , fProbTPC()
    , fProbTOF()
    , fProbTPCTOF()
    , fPriors()

    , fProbTPCTOFnSigTPCMom()
    , fProbTPCnSigTPCMom()
    , fProbTOFnSigTOFMom()
    , fPriorsUsed()


    , fDeDx(0x0)
    , fDeDxTuned(0x0)
    , fDeDxAfterElectronCut(0x0)

    , fDeDx_SPECIES_nSigmaCutTPC()
    , fDeDx_SPECIES_nSigmaCutTOF()
    , fDeDx_SPECIES_nSigmaCutTPCTOF()

    , fTOFvsMom(0x0)
    , fTrackLength(0x0)

    , fProbAllDets()
    , fHistPidDirtyMaxProbability(0)
    , fHistPidPureMaxProbability(0)

    , fHistParticlesDistr(0)
    , fHistParticlesDistrDirty(0)

    // */

    , fPIDnSigmaCut(2.)	// default value is 2
    , fPIDforAnalysis(-1)
    // ######### end of NEW PID 7.11.2016



{
    // Constructor
    fRunKine = runKine;

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


    // Event level QA
    const Int_t nEventStatBins = 16;
    fHistEventCutStats = new TH1I("fHistEventCutStats","Event statistics;;N_{events}", nEventStatBins,-0.5,nEventStatBins-0.5);
    TString gEventCutBinNames[nEventStatBins] =
    {
        "Total",
        "IsPileUpEvent",
        //"No trigger",
        "Not kMB",
        "Not kINT7",
        "IsIncompleteDAQ",
        "BadV0Decision",
        "IsPileupSPD",
        "Wrong centrality",
        "No reco vertex",
        "After vertex cuts",
        "vertex ESD Status=0",
        "Bad vertex (VertexerZ)",
        "Few SPD tracklets",
        "HighMult cut",
        "LowMult cut",
        "Analyzed"
    };
    for(Int_t i = 1; i <= nEventStatBins; i++)fHistEventCutStats->GetXaxis()->SetBinLabel(i,gEventCutBinNames[i-1].Data());
    fOutputList->Add(fHistEventCutStats);

    // Vertex reconstruction QA
    const Int_t nVertexStatBins = 12;
    fHistVertexStats = new TH1I("fHistVertexStats","Vertex statistics;;N_{events}", nVertexStatBins,-0.5,nVertexStatBins-0.5);
    TString gVertexBinNames[] =
    {
        "Total",
        "nContr=0 || cov5=0",
        "HasPrimaryVertexSPD",
        "vtxSPD->IsFromVertexer3D",
        "vtxSPD->IsFromVertexerZ",
        "zRes>fMaxResol",
        "|vSPDz-vtx_z|>0.5",
        "comb.cut on abs and n#sigma dist",
        "too many TPC clusters",
        "isPileupSPD",
        "isPileupFromMV",
        "wrong TOF timing",
        "after cuts"
    };
    for(Int_t i = 1; i <= nVertexStatBins; i++)fHistVertexStats->GetXaxis()->SetBinLabel(i,gVertexBinNames[i-1].Data());
    fOutputList->Add(fHistVertexStats);


    // Vertex distributions
    fHistVzDiffSPDvsTracks = new TH1D("fHistVzDiffSPDvsTracks","diff vertex z (vSPD-vert);diff v_{z} (cm);Entries",300,-15.,15.);
    fHistVzDiffTPCvsTracks = new TH1D("fHistVzDiffTPCvsTracks","diff vertex z (vTPC-vert);diff v_{z} (cm);Entries",300,-15.,15.);
    fHist2DmultTOFvsTracks  = new TH2F("fHist2DmultTOFvsTracks","Pileup from out-of-bunch using TOF;mult (bit32);mult (bit32+TOF)",100,0,5000,100,0,2500);
    fHist2DmultTOFvsTracksAfterCuts  = new TH2F("fHist2DmultTOFvsTracksAfterCuts","Pileup from out-of-bunch using TOF;mult (bit32);mult (bit32+TOF)",100,0,5000,100,0,2500);

    fHist2DmultV0MvsTracksTPCoutBEFOREcutOnLargeMultEsd  = new TH2F("fHist2DmultV0MvsTracksTPCoutBEFOREcutOnLargeMultEsd","V0M vs nTracksWithTPCout BEFORE cut on nEsdTracks;nTracksWithTPCout;multV0M"
                                             ,300,0,30000,500,0,50000);
    fHist2DmultV0MvsTracksTPCout  = new TH2F("fHist2DmultV0MvsTracksTPCout","V0M vs nTracksWithTPCout;nTracksWithTPCout;multV0M"
                                             ,300,0,30000,500,0,50000);


    fHistVertexNconributors = new TH1I("fHistVertexNconributors","Primary vertex n contributors;N contributors;Entries",101,-0.5,100.5);
    fHistNumberOfPileupVerticesTracks = new TH1I("fHistNumberOfPileupVerticesTracks","Number of pileup vertices (by tracks);N vertices;Entries",11,-0.5,10.5);
    fHistNumberOfPileupVerticesSPD = new TH1I("fHistNumberOfPileupVerticesSPD","Number of pileup vertices (by SPD);N vertices;Entries",11,-0.5,10.5);
    fHistNumberOfVertices = new TH1I("fHistNumberOfVertices","Number of vertices;N vertices;Entries",51,-0.5,50.5);

    fOutputList->Add(fHistVzDiffSPDvsTracks);
    fOutputList->Add(fHistVzDiffTPCvsTracks);
    fOutputList->Add(fHist2DmultTOFvsTracks);
    fOutputList->Add(fHist2DmultTOFvsTracksAfterCuts);
    fOutputList->Add(fHist2DmultV0MvsTracksTPCoutBEFOREcutOnLargeMultEsd);
    fOutputList->Add(fHist2DmultV0MvsTracksTPCout);
    fOutputList->Add(fHistVertexNconributors);
    fOutputList->Add(fHistNumberOfPileupVerticesTracks);
    fOutputList->Add(fHistNumberOfPileupVerticesSPD);
    fOutputList->Add(fHistNumberOfVertices);

    fHistVertexZ = new TH1F("fHistVertexZ","Z Coord. of primary vertex;v_{z}, cm;n events",100,-20.0,20.0);
    fHistVertexXY = new TH2F("fHistVertexXY","XY Coord. of primary vertex;v_{x}, cm;v_{y}, cm;n events",80,-0.8,0.8,80,-0.8,0.8);

    fOutputList->Add(fHistVertexXY);
    fOutputList->Add(fHistVertexZ);


    fHist2D_ESTIMATOR_VS_multTPC = new TH2D( "hist2D_ESTIMATOR_VS_multTPC", "hist2D_ESTIMATOR_VS_multTPC;estimator;mult in TPC (bit768)", 4080, -2, 100, 301, -0.5, 3000.5);
    fOutputList->Add(fHist2D_ESTIMATOR_VS_multTPC);
    fHist2D_ESTIMATOR_VS_multTPC_afterCuts = new TH2D( "fHist2D_ESTIMATOR_VS_multTPC_afterCuts", "fHist2D_ESTIMATOR_VS_multTPC_afterCuts;estimator;mult in TPC (bit768)", 4080, -2, 100, 301, -0.5, 3000.5);
    fOutputList->Add(fHist2D_ESTIMATOR_VS_multTPC_afterCuts);


    // ##### lower and upper boundaries for multTPC vs centrality plots (to cut outliers!)
    fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","[0]+[1]*exp([2]*x)",0,90);
    fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","[0]+[1]*exp([2]*x)",0,90);


    if ( fParameterWhichLHCperiod == 0 ) // LHC10h:
    {
        fBorderToCutOutliersLower->SetParameters( -120, 1700, -0.04 );
        fBorderToCutOutliersUpper->SetParameters( -30, 2500, -0.036 );
    }
    else if ( fParameterWhichLHCperiod == 1 ) // LHC15o:
    {
        fBorderToCutOutliersLower->SetParameters( -130, 1900, -0.046 );
        fBorderToCutOutliersUpper->SetParameters( 50, 2900, -0.036 );
    }
    fBorderToCutOutliersLower->DrawClone("same");
    fBorderToCutOutliersUpper->DrawClone("same");



    //
    fHistTrackPt = new TH1F("fHistTrackPt", "PT distribution;p_{T} (GeV/c);dN/dp_{T}", 200, 0., 5);
    fHistTrackEta = new TH1F("fHistTrackEta", "Eta distribution of Tracks;#eta;dN/d#eta", 160, -2.0, 2.0);
    fHistTrackPhi = new TH1F("fHistTrackPhi", "Phi distribution of Tracks;#varphi;dN/d#varphi", 360, 0, TMath::TwoPi() );

    fHistLog = new TH1F("fHistLog","Log Variables",100, -0.5, 99.5);

    fHistAcceptedTracks = new TH1D("fHistAcceptedTracks","N_{ch} - accepted tracks for LRC;N_{ch} accepted;Entries"
                                   ,4000,-0.5,4000-0.5);

    fOutputList->Add(fHistCentrality);
    fOutputList->Add(fHistTrackPt);
    fOutputList->Add(fHistTrackEta);
    fOutputList->Add(fHistTrackPhi);
    fOutputList->Add(fHistLog);
    fOutputList->Add(fHistAcceptedTracks);




    // Setting the FB analysis
    if(0)for( Int_t i=0; i < FB_2D_N_PHI_BINS; i++ )
        for( Int_t j=0; j < FB_2D_N_PHI_BINS; j++ )
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
        for( Int_t binning=0; binning < FB_N_ETAPHI_BINNINGS; binning++ )
        {
            int nBins = nArrEtaPhiBins[binning] * nArrPhiBins[binning];
            winPairs[pt][binning] = new WinPairInfo[ nBins ];
            for( Int_t ep=0; ep < nBins; ep++ )
            {
                WinPairInfo *wp = &winPairs[pt][binning][ep];
                // wp->setMaxNumberOfTracks( binning < 2 ? 1200 : 150 ); // for larger acceptance - larger arrays!

                // pt bins
                wp->setPtRanges( ptMin[pt], ptMax[pt], ptMin[pt], ptMax[pt] );

                // eta-phi bins
                wp->setPhiRanges( 0, TMath::TwoPi(), 0, TMath::TwoPi() );

                if ( binning == 0 )
                {
                    wp->setMaxNumberOfTracks( 1200 ); // for larger acceptance - larger arrays!
                    //cout << "binning=" << binning << endl;
                    wp->setEtaRanges( 0.6-ep*0.1, 0.8-ep*0.1, -0.8+ep*0.1, -0.6+ep*0.1 );
                }
                //else if ( ep < 12 ) // eW 0.4
                else if ( binning == 1 )
                {
                    wp->setMaxNumberOfTracks( 1200 ); // for larger acceptance - larger arrays!
                    wp->setEtaRanges( 0.4-ep*0.1, 0.8-ep*0.1, -0.8+ep*0.1, -0.4+ep*0.1 );
                }
                else if ( 0 ) //if ( binning == 2 ) // FOR STRONGLY INTENSIVE REFINED STUDIES!
                {
                    wp->setMaxNumberOfTracks( 600 ); // for larger acceptance - larger arrays!
                    //cout << "lala" << endl;
                    const int nSpecBins = 21;
                    if (ep%nSpecBins==0)
                        cout << "### test " << ep/nSpecBins << endl;
                    double middle = -0.5+0.1*(ep/nSpecBins);
                    double smallStep = (ep%nSpecBins)*0.005;
                    wp->setEtaRanges( middle+smallStep, middle+0.2+smallStep, middle-0.2-smallStep, middle-smallStep );
                }
                else if ( binning > 1 ) // FOR PHI BINS STUDIES!!!
                {
                    wp->setMaxNumberOfTracks( 150 );

                    int nPhiBins = nArrPhiBins[binning];
                    if ( binning == 2 || binning == 4 || binning == 6 )
                        wp->setEtaRanges( 0.4, 0.8, -0.8, -0.4 );
                    else if ( binning == 3 || binning == 5 || binning == 7 )
                        wp->setEtaRanges( -0.8, 0.8, -0.8, 0.8 );

                    Int_t phi1 = ep/nPhiBins;
                    if ( ep >= nPhiBins*nPhiBins )
                        phi1 -= nPhiBins;
                    Int_t phi2 = ep%nPhiBins;

                    const double phiStep = TMath::TwoPi()/nPhiBins;

                    wp->setPhiRanges( phi1*phiStep, (phi1+1)*phiStep, phi2*phiStep, (phi2+1)*phiStep );
                }


                if(1)cout << "pt: " << wp->pt[0] << "..." << wp->pt[1] << ", " << wp->pt[2] << "..." << wp->pt[3]
                          << ", eta: " << wp->eta[0] << "..." << wp->eta[1] << ", " << wp->eta[2] << "..." << wp->eta[3]
                          << ", phi: " << wp->phi[0] << "..." << wp->phi[1] << ", " << wp->phi[2] << "..." << wp->phi[3]
                          << endl;
            } // end of eta-phi win
        } // end of eta-phi binning
    } // end of pt bins

    // tune 3D histos ingredients-etaPhiWins-centralities
    const int nFBlabels = 36; //34; //31; //28; //24; //16;
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

        "nFnBptFptB",
        "nFnBptF",
        "nFnBptB",
        "nFptF",
        "nBptB",


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
            for( Int_t binning=0; binning < FB_N_ETAPHI_BINNINGS; binning++ )
            {
                int nCentrBinsCalculated = (fMaxCentrality-fMinCentrality)/cWidth[cW];
                int nBins = nArrEtaPhiBins[binning];

                // QA mult:
                int nMultBins = ( binning < 2 ? 1000 : 140 );
                TString strName = Form( "fHistMiltF_ptBin%d_etaPhiBinning%d_cW%d", pt, binning, cW );
                TString strTitle = Form( "%s;FB win ID;mult;centrality", strName.Data() );
                fHistMiltF[pt][cW][binning]
                        = new TH3D( strName, strTitle
                                    , nBins, -0.5, nBins-0.5 // eta-phi window pairs
                                    , nMultBins, -0.5, nMultBins-0.5 // ingredients
                                    , nCentrBinsCalculated, fMinCentrality, fMaxCentrality  // centrality
                                    );
                strName = Form( "fHistMiltB_ptBin%d_etaPhiBinning%d_cW%d", pt, binning, cW );
                strTitle = Form( "%s;FB win ID;mult;centrality", strName.Data() );
                fHistMiltB[pt][cW][binning]
                        = new TH3D( strName, strTitle
                                    , nBins, -0.5, nBins-0.5 // eta-phi window pairs
                                    , nMultBins, -0.5, nMultBins-0.5 // ingredients
                                    , nCentrBinsCalculated, fMinCentrality, fMaxCentrality  // centrality
                                    );


                // hist 3D:
                strName = Form( "fHistFB_ingredients_etaPhi_ptBin%d_etaPhiBinning%d_cW%d", pt, binning, cW );
                strTitle = Form( "%s;FB win ID;FB label;centrality", strName.Data() );

                if ( nArrPhiBins[binning] > 1 )//binning >=2 ) //phi bins -> square nBins to have all combinations
                    nBins *= nBins;

                fHistFB_ingredients_etaPhi_centr[pt][cW][binning]
                        = new TH3D( strName, strTitle
                                    , nBins, -0.5, nBins-0.5 // eta-phi window pairs
                                    , nFBlabels, -0.5, nFBlabels-0.5 // ingredients
                                    , nCentrBinsCalculated, fMinCentrality, fMaxCentrality  // centrality
                                    );
                for( Int_t i=0; i < nFBlabels; i++ )
                    fHistFB_ingredients_etaPhi_centr[pt][cW][binning]->GetYaxis()->SetBinLabel(i+1, strFB_binLabels[i]);

                fOutputList->Add( fHistMiltF[pt][cW][binning] );
                fOutputList->Add( fHistMiltB[pt][cW][binning] );
                fOutputList->Add( fHistFB_ingredients_etaPhi_centr[pt][cW][binning] );
            }
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

                if(0)
                {
                    fOutputList->Add( fHistDptDpt_p1p2[pt][c][ep] );
                    fOutputList->Add( fHistDptDpt_p1[pt][c][ep]   );
                    fOutputList->Add( fHistDptDpt_p2[pt][c][ep]   );
                    fOutputList->Add( fHistDptDpt_n1n2[pt][c][ep] );
                }

            }
            // additional histos to calculate meanPt
            TString strName = Form( "fHistDptDpt_forMeanPt_ptBin%d_c%d", pt, c );
            fHistDptDpt_meanPt[pt][c] = new TH1D( strName, strName, 1200, 0, 6.);

            if(0)
                fOutputList->Add( fHistDptDpt_meanPt[pt][c] );

        }
    }


    //##### PID stuff
    //The common PID object can then be retrieved from the input handler:
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    //if (!fPIDResponse) !!!!!!!!!!!!
    //	AliFatal("This Task needs the PID response attached to the inputHandler");
    for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
    {
        fProbAllDets[ispec]=new TH2D(Form("prob%s_mom_AllDets",AliPID::ParticleName(ispec)),
                                     Form("%s probability vs. momentum;momentum;probability",AliPID::ParticleName(ispec)),
                                     100,0.,20.,50,0.,1.);
        if ( 0 ) // !fFlagSuppressAddingSomeHistos )
        {
            fOutputList->Add(fProbAllDets[ispec]);
        }
    }

    fHistPidDirtyMaxProbability = new TH1D("fHistPidDirtyMaxProbability","HistPidMaxProbability;Probability;Entries",100,0,1);
    fOutputList->Add(fHistPidDirtyMaxProbability);

    fHistPidPureMaxProbability = new TH1D("fHistPidPureMaxProbability","HistPidMaxProbability;Probability;Entries",100,0,1);
    fOutputList->Add(fHistPidPureMaxProbability);

    const int kMySpecies = 6; // AliPID EParticleTypes + my "NotDefined"
    TString gBinParticleNames[kMySpecies] = {"Electron","Muon","Pion","Kaon", "Proton","NotDefined"};

    fHistParticlesDistr = new TH1D("fHistParticlesDistr","Particles Distr;Particle;Entries",kMySpecies,-0.5,kMySpecies-0.5);
    for(Int_t i = 1; i <= kMySpecies; i++)fHistParticlesDistr->GetXaxis()->SetBinLabel(i,gBinParticleNames[i-1].Data());
    fOutputList->Add(fHistParticlesDistr);

    fHistParticlesDistrDirty = new TH1D("fHistParticlesDistrDirty","Particles Distr;Particle;Entries",kMySpecies,-0.5,kMySpecies-0.5);
    for(Int_t i = 1; i <= kMySpecies; i++)fHistParticlesDistrDirty->GetXaxis()->SetBinLabel(i,gBinParticleNames[i-1].Data());
    fOutputList->Add(fHistParticlesDistrDirty);


    // ########### NEW PID 7.11.2016
    fDeDx = new TH2D("hDeDx",";p_{TPC};dE/dx (a.u.)",500,0,5,500,0,500);
    fOutputList->Add(fDeDx);
    fDeDxTuned = new TH2D("hDeDxTuned",";p_{TPC};dE/dx (a.u.)",500,0,5,500,0,500);
    fOutputList->Add(fDeDxTuned);

    fDeDxAfterElectronCut = new TH2D("hDeDxAfterElectronCut",";p_{TPC};dE/dx (a.u.)",500,0,5,500,0,500);
    fOutputList->Add(fDeDxAfterElectronCut);

    fTOFvsMom = new TH2D("fTOFvsMom",";p_{TPC};TOF #beta",500,0,5,500,0, 2);//-4000,4000);
    fOutputList->Add(fTOFvsMom);

    fTrackLength = new TH1D("fTrackLength",";length(cm);n tracks",500,0,10);
    fOutputList->Add(fTrackLength);

    //
    for (Int_t ibin=0;ibin<kPtBins;ibin++)
    {
        fTPCsignalInPtBins[ibin] = new TH1D( Form("fTPCsignalInPtBins_pt%s"
                                                  ,fgkBinMomDesc[ibin] ), ";#Delta dE/dx", 500, 0, 500 );
        fOutputList->Add(fTPCsignalInPtBins[ibin]);
    }

    //
    for (Int_t ispec=0; ispec<5; ++ispec)
        for (Int_t ibin=0;ibin<kPtBins;ibin++)
        {
            fSignalDeltaTPCTOF[ibin][ispec] = new TH2D( Form("fSignalDeltaTPCTOF_pt%s_%s"
                                                             ,fgkBinMomDesc[ibin], AliPID::ParticleName(ispec) ), ";#Delta dE/dx;#Delta t_{TOF}", 200, -30, 30, 200, -500, 1500 );
            fOutputList->Add(fSignalDeltaTPCTOF[ibin][ispec]);
        }
    //
    for (Int_t ispec=0; ispec<5; ++ispec)
        for (Int_t ibin=0;ibin<kPtBins;ibin++)
        {
            fnSigmasTPCTOF[ibin][ispec] = new TH2D( Form("fnSigmasTPCTOF_pt%s_%s"
                                                         ,fgkBinMomDesc[ibin], AliPID::ParticleName(ispec) ), ";n#sigma TPC;n#sigma TOF", 200, -10, 10, 200, -10, 10 );
            fOutputList->Add(fnSigmasTPCTOF[ibin][ispec]);
        }


    // fDeDx_SPECIES_nSigmaCut
    for (Int_t ispec=0; ispec<5; ++ispec)
    {
        //TPC
        fDeDx_SPECIES_nSigmaCutTPC[ispec] = new TH2D( Form("fDeDx_%s_in_sigma_cut_TPC",AliPID::ParticleName(ispec) ),
                                                      Form("%s;p_{TPC};dE/dx (a.u.)",AliPID::ParticleName(ispec)),
                                                      500,0,5,500,0,500);
        fOutputList->Add(fDeDx_SPECIES_nSigmaCutTPC[ispec]);

        //TOF
        fDeDx_SPECIES_nSigmaCutTOF[ispec] = new TH2D( Form("fDeDx_%s_in_sigma_cut_TOF",AliPID::ParticleName(ispec) ),
                                                      Form("%s;p_{TPC};dE/dx (a.u.)",AliPID::ParticleName(ispec)),
                                                      500,0,5,500,0,500);
        fOutputList->Add(fDeDx_SPECIES_nSigmaCutTOF[ispec]);

        //TPCTOF
        fDeDx_SPECIES_nSigmaCutTPCTOF[ispec] = new TH2D( Form("fDeDx_%s_in_sigma_cut_TPCTOF",AliPID::ParticleName(ispec) ),
                                                      Form("%s;p_{TPC};dE/dx (a.u.)",AliPID::ParticleName(ispec)),
                                                      500,0,5,500,0,500);
        fOutputList->Add(fDeDx_SPECIES_nSigmaCutTPCTOF[ispec]);
    }

    // pT spectra for species
    for (Int_t ispec=0; ispec<5; ++ispec)
    {
        fHistPtPID[ispec] = new TH1D(Form("fHistPtPID_%s", AliPID::ParticleName(ispec) ), Form("p_{T} distribution %s; p_{T}, GeV/c; N/dp_{T}", AliPID::ParticleName(ispec) ), 800, 0.0, 20.0);
        fOutputList->Add(fHistPtPID[ispec]);
    }

     /*


    // ####### from ANALYSIS/ANALYSISalice/AliAnalysisTaskPIDCombined.cxx:
    // ------- setup PIDCombined
    fPIDCombined = new AliPIDCombined;
    fPIDCombined->SetDefaultTPCPriors();
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);

    // no light nuclei - no need to call it, this is default
    //  fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);

    for (Int_t ispec=0; ispec<5; ++ispec)
    {
        fProbTPC[ispec]=new TH2D(Form("prob%s_mom_TPC",AliPID::ParticleName(ispec)),
                                 Form("%s probability vs. momentum;momentum;probability",AliPID::ParticleName(ispec)),
                                 100,0.,20.,50,0.,1.);
        fOutputList->Add(fProbTPC[ispec]);
        fProbTPCnSigma[ispec]=new TH2D(Form("prob%s_nSigma_TPC",AliPID::ParticleName(ispec)),
                                       Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
                                       20,-5.,5.,50,0.,1.);
        fOutputList->Add(fProbTPCnSigma[ispec]);

        for (Int_t ibin=0;ibin<kPtBins;ibin++) {
            fProbTPCnSigTPCMom[ibin][ispec]=new TH2D(Form("prob%s_nSigma_TPC (%s)",AliPID::ParticleName(ispec),fgkBinMomDesc[ibin]),
                                                     Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
                                                     20,-5.,5.,50,0.,1.);
            fOutputList->Add(fProbTPCnSigTPCMom[ibin][ispec]);
        }


        fProbTOF[ispec]=new TH2D(Form("prob%s_mom_TOF",AliPID::ParticleName(ispec)),
                                 Form("%s probability vs. momentum;momentum;probability",AliPID::ParticleName(ispec)),
                                 100,0.,20.,50,0.,1.);
        fOutputList->Add(fProbTOF[ispec]);
        fProbTOFnSigma[ispec]=new TH2D(Form("prob%s_nSigma_TOF",AliPID::ParticleName(ispec)),
                                       Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
                                       20,-5.,5.,50,0.,1.);
        fOutputList->Add(fProbTOFnSigma[ispec]);
        for (Int_t ibin=0;ibin<kPtBins;ibin++) {
            fProbTOFnSigTOFMom[ibin][ispec]=new TH2D(Form("prob%s_nSigma_TOF (%s)",AliPID::ParticleName(ispec),fgkBinMomDesc[ibin]),
                                                     Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
                                                     20,-5.,5.,50,0.,1.);
            fOutputList->Add(fProbTOFnSigTOFMom[ibin][ispec]);
        }

        fProbTPCTOF[ispec]=new TH2D(Form("prob%s_mom_TPCTOF",AliPID::ParticleName(ispec)),
                                    Form("%s probability vs. momentum;momentum;probability",AliPID::ParticleName(ispec)),
                                    100,0.,20.,50,0.,1.);
        fOutputList->Add(fProbTPCTOF[ispec]);
        fProbTPCTOFnSigmaTPC[ispec]=new TH2D(Form("probTPCTOF%s_nSigma_TPC",AliPID::ParticleName(ispec)),
                                             Form("%s TPCTOF probability vs. n#sigmaTPC;n#sigma;probability",AliPID::ParticleName(ispec)),
                                             20,-5.,5.,50,0.,1.);
        fOutputList->Add(fProbTPCTOFnSigmaTPC[ispec]);
        for (Int_t ibin=0;ibin<kPtBins;ibin++)
        {
            fProbTPCTOFnSigTPCMom[ibin][ispec]=new TH2D(Form("probTPCTOF%s_nSigma_TPC (%s)",AliPID::ParticleName(ispec),fgkBinMomDesc[ibin]),
                                                        Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
                                                        20,-5.,5.,50,0.,1.);
            fOutputList->Add(fProbTPCTOFnSigTPCMom[ibin][ispec]);
        }



        // basic priors
        fPriors[ispec]=new TH1F(Form("%s_priors",AliPID::ParticleName(ispec)),
                                Form("%s priors vs momentum",AliPID::ParticleName(ispec)),
                                100,0.,20.);
        //        fOutputList->Add(fPriors[ispec]);
        //        switch (ispec) {
        //        case AliPID::kElectron:
        //            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.02);
        //            break;
        //        case AliPID::kMuon:
        //            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.02);
        //            break;
        //        case AliPID::kPion:
        //            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.56);
        //            break;
        //        case AliPID::kKaon:
        //            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.20);
        //            break;
        //        case AliPID::kProton:
        //            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.20);
        //            break;
        //        default:
        //            break;
        //        }
        //        fPIDCombined->SetPriorDistribution((AliPID::EParticleType)ispec,fPriors[ispec]);

        // priors used
        fPriorsUsed[ispec] = new TH2D(Form("%s_priorsUsed",AliPID::ParticleName(ispec)),
                                      Form("%s priors vs transverse momentum;p_{t} (GeV/c);priors",AliPID::ParticleName(ispec)),
                                      100,0.,20.,101,0,1.01);
        fOutputList->Add(fPriorsUsed[ispec]);
    }



    // ####### end of NEW PID 7.11.2016
    */


    // NEW HISTO added to fOutputList here
    PostData(1, fOutputList); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliForwardBackwardAnalysis::UserExec(Option_t *)
{
    //    cout << "TEST! " << GetName() << endl;
    // Main loop
    // Called for each event

    const double cBinMin_dptdpt[] = { 0, };// 5, 10, 20,    30, 40, 50, 60,     70, 0 };
    const double cBinMax_dptdpt[] = { 5, };// 10, 20, 30,   40, 50, 60, 70,     80, 2 };

    const int myFilterBit = fAODTrackCutBit;//128; //768

    // Post output data.
    AliVEvent *event = InputEvent();
    AliAODEvent *lAOD = dynamic_cast<AliAODEvent*>(event);

    if (!lAOD) {
        printf("ERROR: lAOD not available\n");
        fHistLog->Fill(0);
        return;
    }

    fHistLog->Fill(1);

    fRunNumber = lAOD->GetRunNumber();

    // ##### event centrality
    AliCentrality* centrality = lAOD->GetCentrality();
    double centr = (Float_t)centrality->GetCentralityPercentile(fCentralityEstimator);

    if (fParameterWhichLHCperiod==1) // LHC15o
    {
        AliMultSelection *MultSelection = (AliMultSelection * ) lAOD->FindListObject("MultSelection");

        if( !MultSelection) {
            //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
            AliWarning("AliMultSelection object not found!");
        }
        else{
            centr = MultSelection->GetMultiplicityPercentile(fCentralityEstimator);
            //            cout << centr << endl;
            //            cout << " !!! " << MultSelection->GetMultiplicityPercentile("V0M") << endl;
            // !!!: Supported Estimators: V0M, CL0, CL1
        }
    }


    //fCentPercentile[1] = (Float_t)centrality->GetCentralityPercentile("TRK");
    // This is either set upto <~44% or is 100. i.e. cent is set to 100 if 44%<cent<100:
    //fCentPercentile[2] = (Float_t)centrality->GetCentralityPercentile("ZEMvsZDC");



    if(centr<0.0) {
        fHistLog->Fill(3);
        return;
    }
    //float centr = fCentPercentile[0];
    //  cout<<" Got Event of Centrality by V0M = "<<fCentPercentile[0]<<endl;

    if ( centr < fMinCentrality || centr > fMaxCentrality )
        return;

    // ##### choose a particular trigger
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    if(isSelected==kFALSE) {
        //cout<<" Not selected by trigger"<<endl;
        fHistEventCutStats->Fill("Not kMB", 1);
        //return;
    }

    if( !((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7 )
        fHistEventCutStats->Fill("Not kINT7", 1);




    // #### Vertex stuff

    // NOTE: GetPrimaryVertex returns:
    // the vertex reconstructed from tracks when available
    // if the track vertex is not available it returns the SPD vertex
    // if also the SPD vertex is not available it returns the vertex from TPC tracks

    fHistVertexStats->Fill( "Total", 1 );
    const AliVVertex *vertex = lAOD->GetPrimaryVertex();

    Double32_t lCov[6];
    vertex->GetCovarianceMatrix(lCov);

    //check nContributors and z*z>0
    //cut on reco vertex params
    if( !(( vertex->GetNContributors() > 0 ) && ( lCov[5] != 0 )) )
    {
        fHistVertexStats->Fill( "nContr=0 || cov5=0", 1 );
        PostData(1, fOutputList);
        return;
    }

    fHistVertexNconributors->Fill( vertex->GetNContributors() ); // # of tracklets/tracks used for the estimate

    // from https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
    // Selection on SPD vertex type
    const AliVVertex* vtxSPD = lAOD->GetPrimaryVertexSPD();
    if(vtxSPD->IsFromVertexer3D()){
        // vertex with 3D (x,y,z) reconstruction from SPD tracklets
        fHistVertexStats->Fill( "vtxSPD->IsFromVertexer3D", 1 );
        //PostData(1, fOutputList);
        //        return;
    }
    if(vtxSPD->IsFromVertexerZ()){
        // vertex with Z reconstruction from SPD tracklets and x,y from mean vertex
        fHistVertexStats->Fill( "vtxSPD->IsFromVertexerZ", 1 );
        //PostData(1, fOutputList);
        //        return;
    }

    // Special selections for SPD vertex
    // Reject SPD vertices for which only z coordinate is reconstructed and it is determined with poor resolution.
    // Suggested value: fMaxResol=0.25
    const double fMaxResol = 0.25;
    Double_t cov[6]={0};
    vtxSPD->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if ( vtxSPD->IsFromVertexerZ() && (zRes>fMaxResol)) {
        // vertex Z reconstructed with very poor resolution, suggestion is to reject these events
        fHistVertexStats->Fill( "zRes>fMaxResol", 1 );
        PostData(1, fOutputList);
        return;
    }

    // For ESD events, one can apply an additional cleanup cut on the vertex dispersion
    // This is the size of the window opened to define the trackletsused int he determination of the vertex
    // Suggested value: fMaxDipersion=0.03
    //        AliESDVertex* esdVtxSPD = esdEvent->GetPrimaryVertexSPD();
    //        if ( esdVtxSPD->IsFromVertexerZ() && esdVtxSPD->GetDispersion()>fMaxDipersion ) {
    //            // vertex Z reconstructed with very poor resolution, suggestion is to reject these events
    //        }


    // Cut on absolute distance between track and SPD vertices (used in Pb-Pb 2011):
    // Rejects events with badly reconstructed track vertex, typically events with true collision at very large z (~20cm from the center of the detector) with misreconstructed track vetex (details)
    Float_t lDiff_vZ_SPD_tracks = vtxSPD->GetZ() - vertex->GetZ();
    fHistVzDiffSPDvsTracks->Fill( lDiff_vZ_SPD_tracks );
    if( fabs( lDiff_vZ_SPD_tracks ) > 0.5 ) {
        // reject, bad reconstructed track vertex
        fHistVertexStats->Fill( "|vSPDz-vtx_z|>0.5", 1 );
        PostData(1, fOutputList);
        return;
    }

    // IA: check also TPC vertex
    const AliVVertex* vtxTPC = lAOD->GetPrimaryVertexTPC();
    Float_t lDiff_vZ_TPC_tracks = vtxTPC->GetZ() - vertex->GetZ();
    fHistVzDiffTPCvsTracks->Fill( lDiff_vZ_TPC_tracks );

    // Combined cut on absolute and nsigma distance between track and SPD vertices (to be used for Pb-Pb 2015):
    // Rejects events with badly reconstructed track vertex
    // Not needed in pp and p-Pb where the "multivertex" reconstruction of interaction vertices from tracks is used
    if(1)
    {
        double covTrc[6],covSPD[6];
        vertex->GetCovarianceMatrix(covTrc);
        vtxSPD->GetCovarianceMatrix(covSPD);
        double dz = vertex->GetZ()-vtxSPD->GetZ();
        double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
        double errTrc = TMath::Sqrt(covTrc[5]);
        double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
        if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20){
            // reject, bad reconstructed track vertex
            fHistVertexStats->Fill( "comb.cut on abs and n#sigma dist", 1 );
            PostData(1, fOutputList);
            return;
        }
    }

    // ##### Pileup removal (pp and p-Pb)
    if(0)
    {
        // Selections based on multiple interaction vertices (pp and p-Pb)
        // Pileup tagging with SPD multiple vertices
        // Tags events with more than one reconstructed SPD vertex
        // Sensitive to same-bunch pileup and to out-of-bunch pileup for collisions occurring within the SPD readout time window (300 ns)
        // Suggested cut values:
        // fSPDZDiffCut=0.8 (vertices with z separation lower than 8 mm are not considered for tagging)
        // fSPDContributorsCut=3 for low multiplicity pp events (higher efficiency, but also higher contamination from false positives), fSPDContributorsCut=5 at high multiplicity
        Bool_t lEvSel_isPileupSPD = lAOD->IsPileupFromSPD(3,0.8,3.,2.,5.);
        if ( lEvSel_isPileupSPD )
        {
            fHistVertexStats->Fill( "isPileupSPD", 1 );
            PostData(1, fOutputList);
            return;
        }
    }

    // Pileup tagging with track multivertexer (pp and p-Pb)
    // default configuration:
    if(0)
    {
        //        double minContributors=5; double minChi2=5.; double minWeiZDiff=15; double checkPlpFromDifferentBC=kFALSE;
        //        AliAnalysisUtils utils;
        //        utils.SetMinPlpContribMV(minContributors);
        //        utils.SetMaxPlpChi2MV(minChi2);
        //        utils.SetMinWDistMV(minWeiZDiff);
        //        utils.SetCheckPlpFromDifferentBCMV(checkPlpFromDifferentBC);
        //        Bool_t lEvSel_isPileupFromMV = utils.IsPileUpMV(lAOD);
        //        if ( lEvSel_isPileupFromMV )
        //        {
        //            fHistVertexStats->Fill( "isPileupFromMV", 1 );
        //            PostData(1, fOutputList);
        //            return;
        //        }
    }


    //fill QA V0mult vs TPC out BEFORE CUTS on large multEsd (!)
    if (1)
    {
        const Int_t nTracks = lAOD->GetNumberOfTracks();
        Int_t multTrkTPCout=0;
        for (Int_t it = 0; it < nTracks; it++)
        {
            AliAODTrack* aodTrk = (AliAODTrack*)lAOD->GetTrack(it);
            if (!aodTrk) continue;
            if (aodTrk->GetFlags()&AliESDtrack::kTPCout)
                multTrkTPCout++;
        }
        const AliAODVZERO* vzrData = lAOD->GetVZEROData();
        fHist2DmultV0MvsTracksTPCoutBEFOREcutOnLargeMultEsd->Fill( multTrkTPCout, vzrData->GetMTotV0A() + vzrData->GetMTotV0C() );
    }

    // ##### Additional cleanup cuts for LHC15o (https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp)
    // Events with large number of TPC clusters
    if (1)
    {
        // Events with very large number of TPC clusters (~6M, i.e. factor ~3 of single central event), i.e. with strong TPC pile-up.
        // Worse track quality, larger chi2/nclu, significant fraction of tracks rejected by the cut chi2/clu<4.
        // Fall outside the correlation band between number of tracks passing FilterBit32 selection and V0 centrality
        // In ESDs they can be removed via a cut on the number of TPC clusters AliESDEvent::GetNumberOfTPCClusters()
        // In AODs a selection based on the total number of tracks AliAODHeader::GetNumberOfESDTracks() and on the number of TPC only track FIlterBit128
        // Cut still to be fine tuned for very low multiplicities
        Int_t multEsd = ((AliAODHeader*)lAOD->GetHeader())->GetNumberOfESDTracks();
        const Int_t nTracks = lAOD->GetNumberOfTracks();
        Int_t multTPC=0;
        for (Int_t it1 = 0; it1 < nTracks; it1++) {
            AliAODTrack* aodTrk1 = (AliAODTrack*)lAOD->GetTrack(it1);
            if (!aodTrk1) continue;
            if (aodTrk1->TestFilterBit(128)) multTPC++;
        }

        //        if(multEsd -3.38*multTPC<15000){
        if(multEsd -3.38*multTPC<700) { // Pileup cut 15o from Alexandru Dobrin 29 May 2017
            // keep the event
        }else{
            // reject the event
            fHistVertexStats->Fill( "too many TPC clusters", 1 );
            PostData(1, fOutputList);
            return;
        }
    }

    // Out of bunch pileup removal using TOF information
    if (1)
    {
        // Count tracks passing FiltBit32 selections (TPC+ITS+primary selections) and those passing FiltBit32 + TOF Bunch Crossing ID=0 request (tracks originating from the collision that fired the trigger)
        // Apply a 4sigma cut on the correlation
        const Int_t nTracks = lAOD->GetNumberOfTracks();
        Int_t multTrk=0;
        Int_t multTrkTOF=0;
        for (Int_t it1 = 0; it1 < nTracks; it1++) {
            AliAODTrack* aodTrk1 = (AliAODTrack*)lAOD->GetTrack(it1);
            if (!aodTrk1) continue;
            if (aodTrk1->TestFilterBit(32)){
                multTrk++;
                if ( TMath::Abs(aodTrk1->GetTOFsignalDz()) <= 10 && aodTrk1->GetTOFsignal() >= 12000 && aodTrk1->GetTOFsignal() <= 25000) multTrkTOF++;
            }
        }
        fHist2DmultTOFvsTracks->Fill( multTrk, multTrkTOF );
        //        lEvSel_multTrk_bit32 = multTrk;
        //        lEvSel_multTrk_bit32_TOF = multTrkTOF;

        if( ( fParameterWhichLHCperiod==1 // LHC15o
              && ( multTrkTOF < 200./1400*(multTrk-100)
                   || multTrkTOF < 900./2000*(multTrk-400) ) )
                || ( fParameterWhichLHCperiod==0 // LHC10h
                     && ( multTrkTOF < 200./1400*(multTrk-100)) )
                )
        {
            fHistVertexStats->Fill( "wrong TOF timing", 1 );
            PostData(1, fOutputList);
            return;
        }

        fHist2DmultTOFvsTracksAfterCuts->Fill( multTrk, multTrkTOF );
    } // end of out-of-bunch pileup

    //fill QA V0mult vs TPC out
    if (1)
    {
        const Int_t nTracks = lAOD->GetNumberOfTracks();
        Int_t multTrkTPCout=0;
        for (Int_t it = 0; it < nTracks; it++)
        {
            AliAODTrack* aodTrk = (AliAODTrack*)lAOD->GetTrack(it);
            if (!aodTrk) continue;
            if (aodTrk->GetFlags()&AliESDtrack::kTPCout)
                multTrkTPCout++;
        }
        const AliAODVZERO* vzrData = lAOD->GetVZEROData();
        fHist2DmultV0MvsTracksTPCout->Fill( multTrkTPCout, vzrData->GetMTotV0A() + vzrData->GetMTotV0C() );
    }


    fHistVertexStats->Fill( "after QA cuts", 1 );

    Float_t lVertexZ = vertex->GetZ();
    if ( fabs(lVertexZ) > 8 )
    {
        PostData(1, fOutputList);
        return;
    }
    fHistVertexStats->Fill( "after cuts on range", 1 );
    fHistEventCutStats->Fill("After vertex cuts",1);

    Float_t lVertexX = vertex->GetX();
    Float_t lVertexY = vertex->GetY();
    fHistVertexXY->Fill( lVertexX, lVertexY );
    fHistVertexZ->Fill( lVertexZ );

    // #### end of vertex stuff





    fHistLog->Fill(50);

    //    if ( centr > 2 )
    //        return;

    //  cout<<"Centrality Percentile is "<<fCentPercentile[0]<<" "<<fCentPercentile[1]<<" "<<fCentPercentile[2]<<endl;
    fHistCentrality->Fill(centr);

    // ##### count tracks and fill multTPC vs centrality hist2D
    Int_t nAllTracks = lAOD->GetNumberOfTracks();
    Int_t multTPC = 0;
    for (Int_t i = 0; i < nAllTracks; i++)
    {
        AliAODTrack* track = (AliAODTrack*)lAOD->GetTrack(i);
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
    // ##### end of multTPC vs centrality



    // ##### for FB 2D analysis:
    int winInfo_n[FB_2D_N_PHI_BINS];
    double winInfo_pt[FB_2D_N_PHI_BINS];
    for ( int phiBin = 0; phiBin < FB_2D_N_PHI_BINS; phiBin++ )
    {
        winInfo_n[phiBin] = 0;
        winInfo_pt[phiBin] = 0;
    }

    // ##### for FB analysis:
    for( Int_t ptBin=0; ptBin < N_PT_BINS; ptBin++ )
        for( Int_t binning=0; binning < FB_N_ETAPHI_BINNINGS; binning++ )
            for( Int_t ep=0; ep < nArrEtaPhiBins[binning]*nArrPhiBins[binning]; ep++ )
                winPairs[ptBin][binning][ep].reset();



    double phiStep = TMath::TwoPi()/FB_2D_N_PHI_BINS;

    // ##### Track loop
    //  printf("There are %d tracks in this event\n", lAOD->GetNumberOfTracks());
    Int_t NTracks = 0; // This is the accepted track number
    //    Int_t nAllTracks = (lAOD->GetTracks())->GetEntriesFast();
    for (Int_t iTrack = 0; iTrack < nAllTracks; iTrack++)
    {
        AliAODTrack* track = (AliAODTrack*)lAOD->GetTrack(iTrack);
        if (!track) {
            printf("ERROR: Could not receive track %d\n", iTrack);
            continue;
        }

        if( !track->TestFilterBit( myFilterBit ) )
            continue;

        if ( fabs( track->Eta() ) > 1 )
            continue;

        double pt  = (Float_t)track->Pt();
        double phi = (Float_t)track->Phi();
        double eta = (Float_t)track->Eta();
        int charge = (Int_t)track->Charge();


        // ##### GET PID!!!

        bool flagPIDdecision = PIDdecisionLogic( event, iTrack );
        if ( fPIDforAnalysis != -1 && !flagPIDdecision )
            continue;

        // ##### end of PID


        // try ineff by hand!
        if ( 1 && !AcceptTrackAfterIneffMap( eta, pt, centr ) )
            continue;

        NTracks++;


        fHistTrackPt->Fill(pt);
        fHistTrackEta->Fill(eta);
        fHistTrackPhi->Fill(phi);





        // ##### for FB 2D analysis:
        if ( pt > 0.8 && pt < 10.0 )
        {
            if ( eta > -0.8 && eta < 0.8 )
            {
                for ( int phiBin = 0; phiBin < FB_2D_N_PHI_BINS; phiBin++ )
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
            for( Int_t binning=0; binning < FB_N_ETAPHI_BINNINGS; binning++ )
                for( Int_t ep=0; ep < nArrEtaPhiBins[binning]*nArrPhiBins[binning]; ep++ )
                    winPairs[ptBin][binning][ep].addTrack( pt, eta, phi, charge, 0 );


        // ###### DptDpt analysis
        if(0)if ( fabs(eta) < 0.8
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
            for (Int_t jTrack = iTrack+1; jTrack < nAllTracks; jTrack++)
            {
                AliAODTrack* track2 = (AliAODTrack*)lAOD->GetTrack(jTrack);
                if (!track2) {
                    printf("ERROR: Could not receive track %d\n", jTrack);
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
    fHistAcceptedTracks->Fill( fNumberOfTracks );

    //final actions wins for FB 2D:
    for ( int phiBin = 0; phiBin < FB_2D_N_PHI_BINS; phiBin++ )
    {
        if ( winInfo_n[phiBin] > 0 )
            winInfo_pt[phiBin] /= winInfo_n[phiBin];
        else
            winInfo_pt[phiBin] = -1;
    }

    //fill PtPt hist:
    if(0)for ( int binF = 0; binF < FB_2D_N_PHI_BINS; binF++ )
        for ( int binB = 0; binB < FB_2D_N_PHI_BINS; binB++ )
        {
            if ( winInfo_pt[binF]>=0 && winInfo_pt[binB]>=0 )
                fHistPtPt[binF][binB]->Fill( winInfo_pt[binF], winInfo_pt[binB] );
            fHistNN[binF][binB]->Fill( winInfo_n[binF], winInfo_n[binB] );
        }

    // ##### final actions for FB, fill hist:
    for( Int_t ptBin=0; ptBin < N_PT_BINS; ptBin++ )
        for( Int_t binning=0; binning < FB_N_ETAPHI_BINNINGS; binning++ )
            for( Int_t ep=0; ep < nArrEtaPhiBins[binning]*nArrPhiBins[binning]; ep++ )
            {
                WinPairInfo *wp = &winPairs[ptBin][binning][ep];
                wp->finalActionsForEvent();

                //corresponding deltaEta-deltaPhi bin in 3D-histos:

                int binToFill = ep;
                //                if ( binning >= 2 ) // i.e. phi wins
                //                {
                //                    int nPhiBins = nArrPhiBins[binning];
                //                    Int_t phi1 = ep/nPhiBins;
                //                    if ( ep >= nPhiBins*nPhiBins )
                //                        phi1 -= nPhiBins;
                //                    Int_t phi2 = ep%nPhiBins;

                //                    binToFill = ( phi1 > phi2 ? nPhiBins-fabs(phi2-phi1) : fabs(phi2-phi1) );
                //                    if(0)cout << "phi1, phi2: " << phi1 << ", " << phi2 << ", binToFill: " << binToFill << endl;
                //                }

                for( Int_t cW=0; cW < FB_N_CENTR_BINNINGS; cW++ )
                {

                    //            if ( centr < 10 && wp->nF>0 && wp->nB>0 && pt==0)
                    //            {
                    //                cout << "AHNUNG: ptBin=" << pt << ", epBin=" << ep << endl;
                    //            }

                    // fill QA mult:
                    if ( ep < nArrEtaPhiBins[binning] )
                    {
                        fHistMiltF[ptBin][cW][binning]->Fill( binToFill,   wp->nF,    centr );
                        fHistMiltB[ptBin][cW][binning]->Fill( binToFill,   wp->nB,    centr );
                    }


                    // fill FB info into hist bins:
                    TH3D *hist = fHistFB_ingredients_etaPhi_centr[ptBin][cW][binning];

                    hist->Fill( binToFill,   "nEvents",    centr, 1 );
                    hist->Fill( binToFill,   "nF",         centr, wp->nF );
                    hist->Fill( binToFill,   "nB",         centr, wp->nB );
                    hist->Fill( binToFill,   "nF2",        centr, wp->nF*wp->nF );
                    hist->Fill( binToFill,   "nB2",        centr, wp->nB*wp->nB );
                    hist->Fill( binToFill,   "nFnB",       centr, wp->nF*wp->nB );

                    hist->Fill( binToFill,   "sumPtAllEvF",        centr, wp->nF*wp->PtF );
                    hist->Fill( binToFill,   "sumPtAllEvB",        centr, wp->nB*wp->PtB );

                    if ( wp->nF > 0 )
                    {
                        hist->Fill( binToFill,   "nEventsF",   centr, 1 );
                        hist->Fill( binToFill,   "ptF",        centr, wp->PtF );
                        hist->Fill( binToFill,   "ptF2",       centr, wp->PtF*wp->PtF );
                        hist->Fill( binToFill,   "ptFnB",      centr, wp->PtF*wp->nB );
                        hist->Fill( binToFill,   "nFnBptF",     centr, wp->nF*wp->nB*wp->PtF );
                        hist->Fill( binToFill,   "nFptF",     centr, wp->nF*wp->PtF );

                        // for C when av over pairs is OUTSIDE sum:
                        hist->Fill( binToFill,   "pipjF",              centr, wp->pipjF );
                        hist->Fill( binToFill,   "(nF-1)*sum_pF",      centr, (wp->nF-1)*(wp->nF*wp->PtF) ); // (n-1)*sumPt
                        hist->Fill( binToFill,   "nF*(nF-1)",          centr, (wp->nF-1)*wp->nF ); // (n-1)*n

                        // for C when av over pairs is INSIDE sum:
                        int nPairsF = (wp->nF-1)*wp->nF;
                        if ( nPairsF > 0 )
                        {
                            hist->Fill( binToFill,   "pipjF_avPerEv",              centr, wp->pipjF / nPairsF );
                            hist->Fill( binToFill,   "(nF-1)*sum_pF_avPerEv",      centr, (wp->nF-1)*(wp->nF*wp->PtF) / nPairsF );
                        }
                    }
                    if ( wp->nB > 0 )
                    {
                        hist->Fill( binToFill,   "nEventsB",   centr, 1 );
                        hist->Fill( binToFill,   "ptB",        centr, wp->PtB );
                        hist->Fill( binToFill,   "ptB2",       centr, wp->PtB*wp->PtB );
                        hist->Fill( binToFill,   "ptBnF",      centr, wp->PtB*wp->nF );
                        hist->Fill( binToFill,   "nFnBptB",     centr, wp->nF*wp->nB*wp->PtB );
                        hist->Fill( binToFill,   "nBptB",     centr, wp->nB*wp->PtB );

                        // for C when av over pairs is OUTSIDE sum:
                        hist->Fill( binToFill,   "pipjB",              centr, wp->pipjB );
                        hist->Fill( binToFill,   "(nB-1)*sum_pB",      centr, (wp->nB-1)*(wp->nB*wp->PtB) ); // (n-1)*sumPt
                        hist->Fill( binToFill,   "nB*(nB-1)",          centr, (wp->nB-1)*wp->nB ); // (n-1)*n

                        // for C when av over pairs is INSIDE sum: (like in GOOD_Ck_definition_STAR_2005_0504031.pdf)
                        int nPairsB = (wp->nB-1)*wp->nB;
                        if ( nPairsB > 0 )
                        {
                            hist->Fill( binToFill,   "pipjB_avPerEv",              centr, wp->pipjB / nPairsB );
                            hist->Fill( binToFill,   "(nB-1)*sum_pB_avPerEv",      centr, (wp->nB-1)*(wp->nB*wp->PtB) / nPairsB );
                        }
                    }
                    if ( wp->nF > 0 && wp->nB > 0 )
                    {
                        hist->Fill( binToFill,   "nEventsFB",  centr, 1 );
                        hist->Fill( binToFill,   "ptFptB",     centr, wp->PtF*wp->PtB );
                        hist->Fill( binToFill,   "nFnBptFptB",     centr, wp->nF*wp->nB*wp->PtF*wp->PtB );

                        // for dptdpt:
                        hist->Fill( binToFill,   "piFpjB",         centr, wp->piFpjB );
                        hist->Fill( binToFill,   "nF*sum_pB",      centr, wp->nF*(wp->nB*wp->PtB) );
                        hist->Fill( binToFill,   "nB*sum_pF",      centr, wp->nB*(wp->nF*wp->PtF) );
                    }
                } // end of centr binnings
            } // end of eta-phi bins

    //  cout<<" Number of AOD tracks are "<<fNumberOfTracks<<endl;

    PostData(1, fOutputList);
}

//________________________________________________________________________
void AliForwardBackwardAnalysis::Terminate(Option_t *)
{
    //    fHist2D_ESTIMATOR_VS_multTPC->DrawClone("colz");

}




bool AliForwardBackwardAnalysis::AcceptTrackAfterIneffMap( double &eta, double &pt, double &centr )
{
    // 10.03.2016: TMP!!! ??? apply ineff from 3D hist:  (same as for FB tree below! -> tmp?..)
    if ( fAmplifiedInefficiencyForHist3D >= 0 ) //&& ( fHist3D_EffFromOutside || fHist1D_EffFromOutside ) )
    {
        if ( centr > 90 || fabs(eta) > 0.8 )
            return false;

        // 10.06.2016: to set centr bins BY HAND when no V0M signal (no reco)
        //        if ( fFlagSetCentrBinsByHandForAMPTkine )
        //        {
        //            if ( fNumberOfKineTracksInTPC_eta08_pt02 > 1700 )
        //                centr = 2.5;
        //            else if ( fNumberOfKineTracksInTPC_eta08_pt02 > 1400 )
        //                centr = 7.5;
        //            else if ( fNumberOfKineTracksInTPC_eta08_pt02 > 900 )
        //                centr = 15;
        //            else if ( fNumberOfKineTracksInTPC_eta08_pt02 > 600 )
        //                centr = 25;
        //            else if ( fNumberOfKineTracksInTPC_eta08_pt02 > 250 )
        //                centr = 40;
        //            else
        //                centr = 65;
        //        }

        //        cout << " test 1 inside... " << endl;
        int globalBinEff = -1;
        float effFromHist = -1;

        if ( fHist3D_EffFromOutside.GetNbinsX()>0 ) // we have root-file with eff
        {
            globalBinEff = fHist3D_EffFromOutside.FindBin( eta, pt, centr );
            effFromHist = fHist3D_EffFromOutside.GetBinContent(globalBinEff);
        }

        // !!! root file contains eff only for pt 0.2-2.0 (bit768).
        // !!! Complement it up to 10 GeV/c by hand!
        if ( pt > 2.0 ) // !!!
            effFromHist = 0.83;



        if ( effFromHist < 0.01 )
            return false;

        //        if (effFromHist < 0.1 || effFromHist > 0.97 )
        //            cout << " globalBinEff " << globalBinEff << " " <<  effFromHist << endl;

        //        cout << " test 3 inside... " << endl;


        if ( gRandom->Uniform() < fAmplifiedInefficiencyForHist3D*(1-effFromHist) )
            return false;
    }
    //    cout << "returning true" << endl;
    return true;
}





Int_t AliForwardBackwardAnalysis::GetMomBin(Float_t mom)
{
    //
    // Given momentum return histogram to be filled
    //
    if (mom>0. && mom < 0.5) return 0;
    if (mom>=0.5 && mom < 0.7) return 1;
    if (mom>=0.7 && mom < 1.0) return 2;
    if (mom>=1.0 && mom < 1.5) return 3;
    if (mom>=1.5 && mom < 2.0) return 4;
    return kPtBins-1;
}




void AliForwardBackwardAnalysis::GetPID(AliVEvent *event, Int_t iTrack
                                        , float *nSigmas_TPC, float *nSigmas_TOF
                                        , bool &flagIsElectron )
{
    // ###### New PID
    Int_t lMostProbablePIDPure = AliPID::kSPECIES; //-1;
    Int_t lMostProbablePIDdirty = AliPID::kSPECIES; //-1;
    Double_t lMaxProb = 0;

    Double_t probTPC[AliPID::kSPECIES]={0.};
    Double_t probTOF[AliPID::kSPECIES]={0.};
    Double_t probTPCTOF[AliPID::kSPECIES]={0.};

    Double_t length = -999., beta =-999, tofTime = -999., tof = -999., t0 = -999.;
    Double_t c = TMath::C()*1.E-9;// m/ns


    if ( !fRunKine && fPIDResponse )
    {
        AliVTrack *trackV = (AliVTrack*)event->GetTrack(iTrack);
        double eta = trackV->Eta();
        double pt = trackV->Pt();

        Double_t momTPC = trackV->GetTPCmomentum();
        Double_t tpcSignal = trackV->GetTPCsignal();

        Int_t ibin=GetMomBin(momTPC);

        //Double_t signalTOF = track->GetTOFsignal();

        // EDetector detCode = fPIDResponse->AliPIDResponse::kTPC;
        //AliPIDResponse::EDetPidStatus pidStatus = fPIDResponse->CheckPIDStatus( AliPIDResponse::kTPC, track );


        for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
        {
            nSigmas_TPC[ispec] = -1000;
            nSigmas_TOF[ispec] = -1000;

            probTPC[ispec] = -1000;
            probTOF[ispec] = -1000;
            probTPCTOF[ispec] = -1000;
        }

        // dE/dx
        fDeDx->Fill(momTPC,tpcSignal);
        fDeDxTuned->Fill(momTPC,trackV->GetTPCsignalTunedOnData());


        // n sigma TPC 1D
        fTPCsignalInPtBins[ibin]->Fill( tpcSignal );

        // Signal Deltas
        UInt_t detUsedTPC = 0;
        for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
        {
            Double_t signalDeltaTPC, signalDeltaTOF;
            detUsedTPC = fPIDResponse->GetSignalDelta( AliPIDResponse::kTPC, trackV, (AliPID::EParticleType)ispec, signalDeltaTPC );//, Bool_t ratio=kFALSE) const;

            if ( detUsedTPC )
            {
                // nSigmaTPC
                Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(trackV,(AliPID::EParticleType)ispec);
                nSigmas_TPC[ispec] = nSigmaTPC;
                if ( fabs(nSigmas_TPC[ispec]) < fPIDnSigmaCut )
                    fDeDx_SPECIES_nSigmaCutTPC[ispec]->Fill(momTPC,tpcSignal);

                // ### if TOF available
                UInt_t detUsedTOF = fPIDResponse->GetSignalDelta( AliPIDResponse::kTOF, trackV, (AliPID::EParticleType)ispec, signalDeltaTOF );//, Bool_t ratio=kFALSE) const;
                if ( detUsedTOF )
                {
                    fSignalDeltaTPCTOF[ibin][ispec]->Fill( signalDeltaTPC, signalDeltaTOF );

                    // nSigmaTOF
                    Double_t nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(trackV,(AliPID::EParticleType)ispec);
                    nSigmas_TOF[ispec] = nSigmaTOF;
                    if ( fabs(nSigmas_TOF[ispec]) < fPIDnSigmaCut )
                        fDeDx_SPECIES_nSigmaCutTOF[ispec]->Fill(momTPC,tpcSignal);

                    if ( sqrt(nSigmas_TPC[ispec]*nSigmas_TPC[ispec] + nSigmas_TOF[ispec]*nSigmas_TOF[ispec]) < fPIDnSigmaCut )
                        fDeDx_SPECIES_nSigmaCutTPCTOF[ispec]->Fill(momTPC,tpcSignal);

                    // n sigma TPC vs TOF
                    fnSigmasTPCTOF[ibin][ispec]->Fill( nSigmaTPC, nSigmaTOF );
                }
            }
        }


        // !!! ##### check electrons: (from M.Janik AnNote pp: https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/majanik/2015-Sep-09-analysis_note-UntriggeredAngularWithPID_PaperProposal_v1.pdf)
        if(detUsedTPC)
        {
            if ( fabs(nSigmas_TPC[0]) < 3 && fabs(nSigmas_TPC[2]) > 3 && fabs(nSigmas_TPC[3]) > 3 && fabs(nSigmas_TPC[4]) > 3 )
                flagIsElectron = true;
            else
                fDeDxAfterElectronCut->Fill(momTPC,tpcSignal);
//            cout << "flagIsElectron=" << flagIsElectron << endl;
        }




        // TOF signal vs p (for QA)
        double tmpSignalDeltaTOF;
        UInt_t detUsedTOF = fPIDResponse->GetSignalDelta( AliPIDResponse::kTOF, trackV, AliPID::kPion, tmpSignalDeltaTOF );
        if( 1 && detUsedTOF) // from PWG/FLOW/Tasks/AliAnalysisTaskPIDconfig.cxx
        {

            tofTime = trackV->GetTOFsignal();//in ps
            t0 = fPIDResponse->GetTOFResponse().GetStartTime(trackV->P());
            length = trackV->GetIntegratedLength();

            tof -= t0;
            tof = tofTime*1E-3; // ns
            //cout<<"tof = "<<tof<<endl;
            //cout<<"length = "<<length<<endl;
            if (tof > 0  && length > 0)
            {
                length = length*0.01; // in meters
                tof = tof*c;
                beta = length/tof;


                if ( fabs(eta) < 0.9 && pt > 0.2 ) // && length > 3.7 ) // TOF: The internal radius of the cylinder is 3.70 m (https://edms.cern.ch/ui/file/460192/1/ALICE-DOC-2004-002.pdf)
                {
                    fTOFvsMom->Fill( trackV->P()*trackV->Charge(), beta ); //track->GetTOFsignal() );
                    fTrackLength->Fill( length );
                }
            }
            // from Performance paper:
            // TOF detector is a large area array of Multigap Resistive Plate Chambers (MRPC),
            // positioned at 370–399 cm from the beam axis and covering the full azimuth and |η| < 0.9
        }

        /*

        // ####### PID July 9 2014 - verified on 7 November 2016
        // ##### from Pid Combined Task - $ALICE_ROOT/ANALYSIS/AliAnalyTaskPIDCombined
        //Double_t mom=trackV->GetTPCmomentum();
        //            Int_t ibin=GetMomBin(mom);

        fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
        UInt_t detUsed = fPIDCombined->ComputeProbabilities(trackV, fPIDResponse, probTPC);

        if (detUsed != 0) // TPC is available
        {
            for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
            {
                Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(trackV,(AliPID::EParticleType)ispec);
                fProbTPC[ispec]->Fill(momTPC,probTPC[ispec]);
                fProbTPCnSigma[ispec]->Fill(nSigmaTPC,probTPC[ispec]);
                fProbTPCnSigTPCMom[ibin][ispec]->Fill(nSigmaTPC,probTPC[ispec]);


                if ( ( fabs(nSigmaTPC) < fPIDnSigmaCut ) && ( probTPC[ispec] > lMaxProb ) ) //take prob from TPC
                {
                    lMaxProb = probTPC[ispec];
                    lMostProbablePIDPure = ispec;
                }
            }

            // ### compute priors for TPC+TOF, even if we ask just TOF for PID
            fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
            detUsed = fPIDCombined->ComputeProbabilities(trackV, fPIDResponse, probTOF);
            Double_t priors[5];     // check priors used for TOF
            fPIDCombined->GetPriors(trackV,priors,fPIDResponse,detUsed);
            for(Int_t ispec=0;ispec<5;ispec++) fPriorsUsed[ispec]->Fill(TMath::Abs(trackV->Pt()),priors[ispec]);

            if (detUsed != 0)
            {  // TOF is available
                for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
                {
                    Double_t nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(trackV,(AliPID::EParticleType)ispec);
                    fProbTOF[ispec]->Fill(momTPC,probTOF[ispec]);
                    fProbTOFnSigma[ispec]->Fill(nSigmaTOF,probTOF[ispec]);
                    fProbTOFnSigTOFMom[ibin][ispec]->Fill(nSigmaTOF,probTOF[ispec]);

                }
            } // end of TOF is available

            fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
            detUsed = fPIDCombined->ComputeProbabilities(trackV, fPIDResponse, probTPCTOF);
            if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask() )
            {
                for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
                {
                    Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(trackV,(AliPID::EParticleType)ispec);
                    fProbTPCTOF[ispec]->Fill(momTPC,probTPCTOF[ispec]);
                    fProbTPCTOFnSigmaTPC[ispec]->Fill(nSigmaTPC,probTPCTOF[ispec]);
                    fProbTPCTOFnSigTPCMom[ibin][ispec]->Fill(nSigmaTPC,probTPCTOF[ispec]);

                    bool isInNsigma = (fabs(nSigmaTPC) < fPIDnSigmaCut);
                    // try "elliptic" nSigmaCut if pT>0.5 GeV/c:
                    if ( momTPC > 0.5 )
                    {
                        Double_t nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(trackV,(AliPID::EParticleType)ispec);
                        Double_t ell_nSigma = sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF);
                        //cout << ">>> ell_nSigma = " << ell_nSigma << endl;
                        isInNsigma = (ell_nSigma < fPIDnSigmaCut);
                    }

                    if ( isInNsigma && (probTPCTOF[ispec] > lMaxProb) ) //update prob from TPC+TOF
                    {
                        lMaxProb = probTPCTOF[ispec];
                        lMostProbablePIDPure = ispec;
                    }
                }
            }

            // QA printing final:
            if(0)
            {
                printf( "momTPC=%.2f, lMaxProb = %.2f \n",  momTPC, lMaxProb );
                for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
                {
                    printf( "%d: sigmaTPC: %.1f, sigmaTOF: %.1f, prob TPC, TOF, TPCTOF: %.2f %.2f --> %.2f\n"
                            , ispec, nSigmas_TPC[ispec], nSigmas_TOF[ispec], probTPC[ispec], probTOF[ispec], probTPCTOF[ispec] );
                }
            }

            // !!! ##### check electrons: (from M.Janik AnNote pp: https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/majanik/2015-Sep-09-analysis_note-UntriggeredAngularWithPID_PaperProposal_v1.pdf)
            if ( fabs(nSigmas_TPC[0]) < 3 && fabs(nSigmas_TPC[2]) > 3 && fabs(nSigmas_TPC[3]) > 3 && fabs(nSigmas_TPC[4]) > 3 )
            {
                if(0)cout << "!!! PID electron warning: nSigma_TPC_electron = " << nSigmas_TPC[0]
                          << " lMostProbablePIDPure = " << lMostProbablePIDPure
                          << " lMaxProb = " << lMaxProb << endl;
                lMaxProb = probTPC[0]; // assign prob of electron
                lMostProbablePIDPure = AliPID::kElectron; // assign id of electron
            }

            //fill for watching at max prob-s in pid species arrays
            fHistPidPureMaxProbability->Fill( lMaxProb );
            if ( lMaxProb < 0.5 ) //here is MY CUT on max probability! IMPORTANT TO REMEMBER!!!
                lMostProbablePIDPure = AliPID::kSPECIES; //my "NotDefined" bin
            fHistParticlesDistr->Fill( lMostProbablePIDPure );


            // fill pT distributions for each PID:
            if ( lMostProbablePIDPure < AliPID::kSPECIES )
                fHistPtPID[lMostProbablePIDPure]->Fill( pt );
        }

        // ####### end of PID July 9 2014

        const bool useDirtyPID = true;
        if ( useDirtyPID )
        {
            //Int_t lMostProbablePIDdirty = -1;
            Double_t lMaxProbDirty = 0;
            for ( Int_t ispec=0; ispec < (int)(AliPID::kSPECIES); ++ispec )
            {
                fProbAllDets[ispec]->Fill( pt, probTPCTOF[ispec] );
                if ( probTPCTOF[ispec] > lMaxProbDirty )
                {
                    lMaxProbDirty = probTPCTOF[ispec];
                    lMostProbablePIDdirty = ispec;	// define most probable particle!!!
                }
            }
            //fill for watching at max prob-s in pid species arrays
            fHistPidDirtyMaxProbability->Fill( lMaxProbDirty );
            fHistParticlesDistrDirty->Fill( lMostProbablePIDdirty );
        }

         */

    }
    // end of new PID filling
}




bool AliForwardBackwardAnalysis::PIDdecisionLogic(AliVEvent *event, Int_t iTrack
                                                  //, float *nSigmas_TPC, float *nSigmas_TOF
                                                  )
{
    // return false;
    float nSigmas_TPC[AliPID::kSPECIES] = {-1000};
    float nSigmas_TOF[AliPID::kSPECIES] = {-1000};
    bool flagIsElectron = false;
    GetPID( event, iTrack, nSigmas_TPC, nSigmas_TOF, flagIsElectron );

    //   return false;
    //   /*

    AliVTrack *trackV = (AliVTrack*)event->GetTrack(iTrack);
    Double_t momTPC = trackV->GetTPCmomentum();
    if(0)
    {
        cout << "track pT: momTPC= " << momTPC << endl;
        cout << "nSigmaTPC: " << nSigmas_TPC[0] << " " << nSigmas_TPC[1] << " " << nSigmas_TPC[2] << " "
             <<  nSigmas_TPC[3] << " " << nSigmas_TPC[4] << " " << endl;
        cout << "nSigmas_TOF: " << nSigmas_TOF[0] << " " << nSigmas_TOF[1] << " " << nSigmas_TOF[2] << " "
             <<  nSigmas_TOF[3] << " " << nSigmas_TOF[4] << " " << endl;
    }

    // from https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/rbertens/2017-May-29-analysis_note-vn_run2_note.pdf:
    // Only tracks that pass the filter bit 32 requirements [11] were selected within the pseudo-rapidity range
    // |η| < 0.8 and 0.2 < pT < 50. GeV/c. Particle identification (PID) is done using the n-σ approach from
    // the AliPIDResponse class by combining information from TPC and TOF sqrt(nSigmaTPC^2+nSigmaTOF^2)<3
    // for pT >= 0.5 GeV/c, while n-σT PC is used for pT < 0.5 GeV/c. The AliPIDHelper class was employed
    // to remove tracks from the overlapping n-σ bands of different particle species and double counting. P

    // ##### selection by PID:
    if ( fPIDforAnalysis != -1 )
    {
        // reject electrons!
        if ( flagIsElectron )
            return false;

        int idMin = 1000;
        double minSigma = 1000;
        double minSigmaPrev = 1000;

        for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
        {
            double nSigma;
            if ( (momTPC < 0.5) // TPC-only for pT<0.5
                 || ( fabs(nSigmas_TOF[ispec]) > 999 )  // TOF is unavailable?
                 )
                nSigma = fabs(nSigmas_TPC[ispec]); // use TPC nSigma
            else // TPC+TOF
                nSigma = sqrt( nSigmas_TPC[ispec]*nSigmas_TPC[ispec] + nSigmas_TOF[ispec]*nSigmas_TOF[ispec] );

            if ( ispec == 1 ) // omit muons
                continue;

            if ( nSigma < minSigma )
            {
                idMin = ispec;
                minSigmaPrev = minSigma;
                minSigma = nSigma;
            }
        }

        if(0)cout << "track momTPC= " << momTPC << ", idMin: "
                  << idMin << ", minSigma=" << minSigma << ", minSigmaPrev=" << minSigmaPrev
                  << ", nSigmas_TPC[ispec]=" << nSigmas_TPC[idMin] << ", nSigmas_TOF[ispec]=" << nSigmas_TOF[idMin] << endl;

        // PID decision:
        if ( idMin == fPIDforAnalysis
             && minSigma < fPIDnSigmaCut
             && minSigmaPrev > 3 // !! reject if some other species within 3 sigma!
             )
            return true;
    }

    return false;
    // */
}






