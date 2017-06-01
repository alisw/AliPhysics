/**************************************************************************
 * Authors: Andrey Ivanov, Igor Altsybeev, St.Petersburg State University *
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
// Analysis task for Long Range Correlation (LRC) analysis using TPC data
// This includes a TList of AliLRCBase objects that are processing LRC analysis
// for a given Eta and Phi windows

// Email: Igor.Altsybeev@cern.ch

#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliAnalysisTaskLRC.h"
#include <AliLRCBase.h>
//#include <AliLRCProcess.h>
#include <AliVEvent.h>
#include <AliMCEvent.h>
#include <AliESDEvent.h>     
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"

#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"

#include <AliPID.h>
#include <AliPIDCombined.h>
#include <AliPIDResponse.h>
#include <AliMultiplicity.h>
#include <AliESDv0.h>

#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliStack.h>
#include <AliESDtrackCuts.h>
#include <AliKalmanTrack.h>
#include <AliPhysicsSelection.h>
#include <AliCentrality.h>
#include <AliESDZDC.h>
#include <AliESDFMD.h>

#include "AliEventplane.h"


#include <TParticle.h>
#include <TH1.h>
#include <TH2.h>
#include "TH1I.h"
#include "TRandom3.h"
#include <TComplex.h>

#include <TF1.h>

//#include "AliSimpleEvent.h"

//flag use my Event Tree
//#define fSetIncludeEventTreeInOutput 1


#include "TStopwatch.h"

using std::endl;
using std::cout;

ClassImp(AliAnalysisTaskLRC)

//________________________________________________________________________
AliAnalysisTaskLRC::AliAnalysisTaskLRC( const char *name, Bool_t runKine)
    : AliAnalysisTaskSE(name)
    ,fNumberOfPhiSectors(1)
    ,fFlagSuppressAddingSomeHistos(kTRUE)
    ,fAnalysisLevel("ESD")
    ,fEsdTrackCuts(0)
    ,fAODtrackCutBit(128)
    ,fArrTrackCuts(0x0)
    ,fHistCutsNamesBins(0)
    //    ,fNumberOfCutsToRemember(0)
    ,fSwitchToListingCuts(0)
    
    //    ,fAnalysisType(en_AnalysisType_ESD)
    ,fMinNumberOfSPDtracklets(0)
    ,fMaxPtLimit(100.0)
    ,fMinPtLimit(0.0)
    ,fKineLowPtCut(0.0)
    ,fMinAcceptedTracksCut(0)
    ,fMaxAcceptedTracksCut(0)
    ,fCheckForkVtx(kTRUE)
    ,fCheckForVtxPosition(kFALSE)
    ,fVxMax(0)
    ,fVyMax(0)
    ,fVzMax(0)
    ,fLRCproc(0)
    ,fOutList(0)
    ,fRunKine(0)
    ,fAnalyseMCESD(0)
    ,fShowEventStats(kFALSE)
    ,fShowPerTrackStats(kFALSE)
    ,fEtaRegionForTests(0.8)
    ,fMultCutInEtaRegion(0)
    ,fHistEventCutStats(0)
    ,fHistTrackCutStats(0)
    ,fHistAODTrackStats(0)
    ,fHistVx(0)
    ,fHistVy(0)
    ,fHistVz(0)
    ,fHistVxMCrecoDiff(0)
    ,fHistVyMCrecoDiff(0)
    ,fHistVzMCrecoDiff(0)
    ,fHistVertexNconributors(0)
    ,fHistNumberOfPileupVerticesTracks(0)
    ,fHistNumberOfPileupVerticesSPD(0)
    ,fHistEventPlane(0)
    ,fHistPt(0)
    ,fHistEta(0)
    ,fHistEtaAODpure(0)
    ,fHistPhi(0)
    ,fHistEtaPhi(0)
    ,fHistPhiLRCrotationsCheck(0)
    ,fHistPhiArtificialProfilingCheck(0)
    ,fHistPhiArtificialProfilingCheckWrtEvPlane(0)
    ,fHistPhiArtificialEvPlane(0)
    ,fHistEtaVsZvCoverage(0)
    ,fHistEtaVsZvCoverageAccepted(0)
    ,fHistMultBeforeCuts(0)
    ,fHistAcceptedMult(0)
    ,fHistAcceptedTracks(0)
    ,fHistMultiplicityInEtaRegion(0)
    ,fHistMultiplicityInEtaRegionAfterPtCuts(0)
    ,fHist2DMultiplicityMCESDInEtaRegion(0)
    ,fHistAcceptedTracksAfterPtCuts(0)
    ,fHistAcceptedTPCtracks(0)
    ,fHistClustersTPC(0)
    ,fHistClustersTPCafterCuts(0)
    ,fHistCrossedRowsTPC(0)
    ,fHistCrossedRowsTPCafterCuts(0)
    ,fHistClustersITS(0)
    ,fHistTrackletsITS(0)
    ,fHist2DClustersTPCvsPt(0)
    ,fHist2DClustersTPCvsEta(0)
    ,fHist2DAcceptedTracksPtvsEta(0)
    ,fHistMClabels(0)
    ,fHistRejectedTracksCharge(0)
    ,fHistTracksCharge(0)
    ,fHistPtPlus(0)
    ,fHistPtMinus(0)
    ,fHistNetChargeVsPt(0)
    ,fHistChargePlusVsPtTmp(0)
    ,fHistChargeMinusVsPtTmp(0)
    ,fHist2DNetChargeVsPt(0)
    ,fHist2DNetChargeVsPtCorrectedOnEventMean(0)
    ,fHist2DNetChargeVsPtCorrectedOnEventMeanNormOnNch(0)
    ,fHistProbabilitiesPID(0)
    ,fHistESDtrackMass(0)
    ,fHistProbabilityPion(0)
    ,fHistProbabilityKaon(0)
    ,fHistProbabilityProton(0)
    ,fHistParticlesDistr(0)
    ,fHistParticlesDistrBeforeCuts(0)
    //    ,fNumberOfSectors(1)
    ,fHistCentralityPercentile(0)
    ,fHistCentralityClass10(0)
    ,fHistCentralityClass5(0)
    //,fHistZDCenergy(0x0)
    ,fHistZDCparticipants(0)
    ,fHistV0multiplicity(0)
    ,fHistV0Amultiplicity(0)
    ,fHistV0Cmultiplicity(0)
    ,fHist2DV0ACmultiplicity(0)
    ,fHist2DTracksAcceptedVsV0multiplicity(0)
    //,fHistV0spectra(0)
    
    ,fHistV0cells    (0)
    ,fHistV0Acells   (0)
    ,fHistV0Ccells   (0)
    ,fHist2DV0ACcells(0)
    
    ,fThresholdOnV0mult(0)
    
    
    ,fMinCentralityClass(-0.01)
    ,fMaxCentralityClass(98)
    ,fIsIonsAnalysis(kFALSE)
    ,fEtInsteadOfPt(kFALSE)
    ,fUsePhiShufflingByHand(kFALSE)
    ,fUseToyEvents(kFALSE)
    ,fNumberOfToyEvents(1)
    ,fTmpCounter(0)
    ,fPIDResponse(0x0)
    ,fPIDCombined(0x0)
    
    ,fHistPidMaxProbability(0)
    ,fHistPidPureMaxProbability(0)
    ,fStrPIDforFwd("any")
    ,fStrPIDforBwd("any")

    ,fPIDsensingFlag(kFALSE)
    ,fArtificialInefficiency(-1.)
    ,fHistNumberOfDroppedByHandTracks(0)
    ,fRand(0)
    
    
    ,fPhiArtificialGapBegin(0)
    ,fPhiArtificialGapEnd(0)
    
    ,fFlagWatchZDC(0)
    ,fFlagWatchV0 (0)
    ,fFlagWatchFMD(0)
    
    ,fAnalysisTimer(0)
    
    //MC qa and studies
    //    ,fHistMCvertexRdeltaFromParent(0)
    //    ,fHistMCparentsStat(0)
    //,fHistMCparentsEta(0)
    //,fHistMCchildsEta(0)
    //,fHistMCdeltaEtaChildParent(0)
    //,fProbabilitiesPID(0)
    
    
    //    ,fSimpleEvent(0)
    //    ,fNsimpleEvents(0)
    //    ,fEventTree(0)
    //    ,fSetIncludeEventTreeInOutput(0)
{
    fAnalysisTimer = new TStopwatch;
    //Init
    fRunKine = runKine;
    for ( int i = 0; i < 5; i++ )
    {
        fHistZDCenergy[i] = 0x0;
    }
    for ( int i = 0; i < 4; i++ )
    {
        fHistV0AmultiplicityRing[i] = 0x0;
        fHistV0CmultiplicityRing[i] = 0x0;
        fHist2DV0ACmultiplicityRing[i] = 0x0;
        fHist2DTracksAcceptedVsV0AmultiplicityRing[i] = 0x0;
        fHist2DTracksAcceptedVsV0CmultiplicityRing[i] = 0x0;
    }
    //fProbabilitiesPID = new Double_t[AliPID::kSPECIES];
    // Output slot #1 writes into a TList container for common task data and QA
    DefineOutput(1, TList::Class());
    
    //Defining output slots for each LRC processor (required to avoid TList of TLists on merging)
    for(Int_t i=0; i < fLRCproc.GetEntries(); i++)
    {
        DefineOutput(Proc(i)->GetOutputSlotNumber(),TList::Class());
    }
}


// ---------------------------------------  Setters ------------------

void AliAnalysisTaskLRC::SetMaxPtLimit(Double_t MaxPtLimit)
{
    //Sets  Max Pt filter
    fMaxPtLimit = MaxPtLimit;
}
void AliAnalysisTaskLRC::SetMinPtLimit(Double_t MinPtLimit)
{
    //Sets  Min Pt filter
    fMinPtLimit = MinPtLimit;
}
void AliAnalysisTaskLRC::SetKineLowPtCut(Double_t kineLowPtCut)
{
    //Sets  Min Pt filter for Kine tracks
    fKineLowPtCut = kineLowPtCut;
}
AliLRCBase*  AliAnalysisTaskLRC::Proc(Int_t index)
{
    // Get Processor i
    return (dynamic_cast<AliLRCBase*> (fLRCproc.At(index)));
}

void AliAnalysisTaskLRC::SetMinNumberOfSPDtracklets( Int_t MinSPDtracklets )
{
    //Sets  Min SPD tracklets number
    fMinNumberOfSPDtracklets = MinSPDtracklets;
}

//________________________________________________________________________
void AliAnalysisTaskLRC::UserCreateOutputObjects()
{
    // --------- Output list
    fOutList = new TList();
    
    //### added 13.12.2011
    fOutList->SetOwner();  // IMPORTANT!
    
    fAnalysisTimer->Start();
    
    
    //##### PID stuff
    // ------- setup PIDCombined
    fPIDCombined=new AliPIDCombined;
    fPIDCombined->SetDefaultTPCPriors();
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
    
    // no light nuclei
    fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
    //Int_t lNumberOfPidSpecies = (Int_t)AliPID::kSPECIES;
    for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
    {
        fProbTPCTOF[ispec]=new TH2D(Form("prob%s_mom_TPCTOF",AliPID::ParticleName(ispec)),
                                    Form("%s probability vs. momentum;momentum;probability",AliPID::ParticleName(ispec)),
                                    100,0.,20.,50,0.,1.);
        fProbAllDets[ispec]=new TH2D(Form("prob%s_mom_AllDets",AliPID::ParticleName(ispec)),
                                     Form("%s probability vs. momentum;momentum;probability",AliPID::ParticleName(ispec)),
                                     100,0.,20.,50,0.,1.);
        if ( !fFlagSuppressAddingSomeHistos )
        {
            fOutList->Add(fProbTPCTOF[ispec]);
            fOutList->Add(fProbAllDets[ispec]);
        }
        
        
        // basic priors
        fPriors[ispec] = new TH1F(Form("%s_priors",AliPID::ParticleName(ispec)),
                                  Form("%s priors vs momentum",AliPID::ParticleName(ispec)),
                                  100,0.,20.);
        if ( !fFlagSuppressAddingSomeHistos )
            fOutList->Add(fPriors[ispec]);
        switch (ispec) {
        case AliPID::kElectron:
            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.02);
            break;
        case AliPID::kMuon:
            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.02);
            break;
        case AliPID::kPion:
            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.56);
            break;
        case AliPID::kKaon:
            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.20);
            break;
        case AliPID::kProton:
            for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.20);
            break;
        default:
            break;
        }
        fPIDCombined->SetPriorDistribution((AliPID::EParticleType)ispec,fPriors[ispec]);
        
        // priors used
        fPriorsUsed[ispec] = new TH2D(Form("%s_priorsUsed",AliPID::ParticleName(ispec)),
                                      Form("%s priors vs transverse momentum;p_{t} (GeV/c);priors",AliPID::ParticleName(ispec)),
                                      100,0.,20.,101,0,1.01);
        if ( !fFlagSuppressAddingSomeHistos )
            fOutList->Add(fPriorsUsed[ispec]);
    }
    
    fHistPidMaxProbability = new TH1D("fHistPidMaxProbability","HistPidMaxProbability;Probability;Entries",100,0,1);
    fOutList->Add(fHistPidMaxProbability);
    
    fHistPidPureMaxProbability = new TH1D("fHistPidPureMaxProbability","HistPidMaxProbability;Probability;Entries",100,0,1);
    fOutList->Add(fHistPidPureMaxProbability);
    
    // ##### end PID part
    Printf("UserCreateOutputObjects.........");
    //Disabling "replacing existing TH2D ...etc" warning
    
    Bool_t lTH1oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    
    Bool_t lTH2oldStatus = TH2::AddDirectoryStatus();
    TH2::AddDirectory(kFALSE);
    
    
    
    //LRC processors init
    Int_t lLrcNum = fLRCproc.GetEntries();
    for(Int_t i = 0; i < lLrcNum; i++)
    {
        AliLRCBase *lrcBase = dynamic_cast<AliLRCBase*> (fLRCproc.At(i));
        if(lrcBase)
            lrcBase->InitDataMembers();
        else continue;
        //remember pointer to be used in the analysis
        //fLRCprocArrayPointers[i] = lrcBase;
    }
    
    
    //fOutList->Add( fEsdTrackCuts );
    
    
    //create here array with bin names
    //15.03.2013 - list with cuts to apply on tracks and remember decisions
    if ( fSwitchToListingCuts )
    {
        int lNumberOfCutsToRemember = fArrTrackCuts.GetEntries();
        fHistCutsNamesBins = new TH1I("fHistCutsNamesBins","Tracks in different cut conditions;;N_{tracks}", lNumberOfCutsToRemember,-0.5,lNumberOfCutsToRemember-0.5);
        for(Int_t i = 1; i <= lNumberOfCutsToRemember; i++)
        {
            //            TString cutsName = (TString*)fArrCutsNames[i-1];
            //            fHistCutsNamesBins->GetXaxis()->SetBinLabel(i,((TString*)fArrCutsNames[i-1])->Data());
            fHistCutsNamesBins->GetXaxis()->SetBinLabel(i,fArrCutsNames[i-1].Data());
            
        }
        fOutList->Add(fHistCutsNamesBins);
    }
    
    
    //Event level QA
    const Int_t nEventStatBins = 10;
    fHistEventCutStats = new TH1I("fHistEventCutStats","Event statistics;;N_{events}", nEventStatBins,-0.5,nEventStatBins-0.5);
    TString gEventCutBinNames[nEventStatBins] = {"Total","No trigger","Wrong centrality","No reco vertex","Bad vertex params", "Bad vertex position","Few SPD tracklets","HighMult cut","LowMult cut","Analyzed"};
    for(Int_t i = 1; i <= nEventStatBins; i++)fHistEventCutStats->GetXaxis()->SetBinLabel(i,gEventCutBinNames[i-1].Data());
    fOutList->Add(fHistEventCutStats);
    
    //Track level QA
    const Int_t nTracksStatBins = 6;
    fHistTrackCutStats = new TH1I("fHistTrackCutStats","Track statistics;;N_{tracks}", nTracksStatBins,-0.5,nTracksStatBins-0.5);
    TString gTrackCutBinNames[nTracksStatBins] = {"Total","QA track cuts","HighPt cut","LowPt cut","Good", "Not phys primary"};
    for(Int_t i = 1; i <= nTracksStatBins; i++)fHistTrackCutStats->GetXaxis()->SetBinLabel(i,gTrackCutBinNames[i-1].Data());
    fOutList->Add(fHistTrackCutStats);
    
    //AOD track cut bits stats
    const int nAODtrackStats = 16;
    fHistAODTrackStats = new TH1I("fHistAODTrackStats","AOD tracks statistics;TrackFilterBit;N_{tracks}",nAODtrackStats,-0.5,nAODtrackStats-0.5);
    if ( fAnalysisLevel == "AOD" )
        fOutList->Add(fHistAODTrackStats);
    
    
    //Vertex distributions
    fHistVx = new TH1D("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm);Entries",200,-0.4,0.4);
    fOutList->Add(fHistVx);
    fHistVy = new TH1D("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm);Entries",200,-0.4,0.4);
    fOutList->Add(fHistVy);
    fHistVz = new TH1D("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Entries",200,-20.,20.);
    fOutList->Add(fHistVz);
    
    if ( fRunKine )
    {
        fHistVxMCrecoDiff = new TH1D("fHistVxMCrecoDiff","Primary vertex mc-reco diff - x coordinate;V_{x} (cm);Entries",200,-0.2,0.2);
        fOutList->Add(fHistVxMCrecoDiff);
        fHistVyMCrecoDiff = new TH1D("fHistVyMCrecoDiff","Primary vertex mc-reco diff - y coordinate;V_{y} (cm);Entries",200,-0.2,0.2);
        fOutList->Add(fHistVyMCrecoDiff);
        fHistVzMCrecoDiff = new TH1D("fHistVMCrecoDiffz","Primary vertex mc-reco diff - z coordinate;V_{z} (cm);Entries",200,-10.,10.);
        fOutList->Add(fHistVzMCrecoDiff);
    }
    
    fHistVertexNconributors = new TH1I("fHistVertexNconributors","Primary vertex n contributors;N contributors;Entries",101,-0.5,100.5);
    fOutList->Add(fHistVertexNconributors);
    
    fHistNumberOfPileupVerticesTracks = new TH1I("fHistNumberOfPileupVerticesTracks","Number of pilup verteces (by tracks);N verteces;Entries",11,-0.5,10.5);
    fOutList->Add(fHistNumberOfPileupVerticesTracks);
    
    fHistNumberOfPileupVerticesSPD = new TH1I("fHistNumberOfPileupVerticesSPD","Number of pilup verteces (by SPD);N verteces;Entries",11,-0.5,10.5);
    fOutList->Add(fHistNumberOfPileupVerticesSPD);
    
    
    if ( fIsIonsAnalysis )
    {
        //Event plane
        fHistEventPlane = new TH2F("fHistEventPlane",";#Psi, rad.;Centrality percentile;Counts",80,0,TMath::Pi()/*0,360.*/,101,-0.5,100.5);
        fOutList->Add(fHistEventPlane);
    }
    //pt, eta, phi checkplots
    fHistPt = new TH1F("fHistPt", "p_{T} distribution", 120, 0.0, 6.0);
    fHistPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistPt->GetYaxis()->SetTitle("dN/dp_{T} (c/GeV)");
    fHistPt->SetMarkerStyle(kFullCircle);
    fOutList->Add(fHistPt);
    
    fHistEta = new TH1F("fHistEta", "#eta distribution", 200, -4, 4);
    fHistEta->GetXaxis()->SetTitle("#eta");
    fHistEta->GetYaxis()->SetTitle("dN/#eta");
    fHistEta->SetMarkerStyle(kFullCircle);
    fOutList->Add(fHistEta);
    
    fHistEtaAODpure = new TH1F("fHistEtaAODpure", "#eta distribution AOD before cuts;#eta;dN/d#eta", 200, -4, 4);
    fOutList->Add(fHistEtaAODpure);
    
    fHistPhi = new TH1F("fHistPhi", "#phi distribution", 200, 0, 2*TMath::Pi() );
    fHistPhi->GetXaxis()->SetTitle("#phi");
    fHistPhi->GetYaxis()->SetTitle("dN/#phi");
    fHistPhi->SetMarkerStyle(kFullCircle);
    fOutList->Add(fHistPhi);
    
    fHistEtaPhi = new TH2D("fHistEtaPhi","N tracks in (#eta, #phi);#eta;#phi",25,-0.8,0.8,25,0,2*TMath::Pi());
    fOutList->Add(fHistEtaPhi);
    
    fHistPhiLRCrotationsCheck = new TH1F("fHistPhiLRCrotationsCheck", "#phi distribution in LRC rotations;#phi;dN/#phi", 200, -4*TMath::Pi(), 4*TMath::Pi() );
    fOutList->Add(fHistPhiLRCrotationsCheck);
    
    //some histos for tracks manipulations by hand
    fHistPhiArtificialProfilingCheck = new TH1F("fHistPhiArtificialProfilingCheck", "#phi distribution tracks distr by hand;#phi;dN/#phi", 200, -4*TMath::Pi(), 4*TMath::Pi() );
    fHistPhiArtificialProfilingCheckWrtEvPlane = new TH1F("fHistPhiArtificialProfilingCheckWrtEvPlane", "#phi distribution by hand wrt event plane;#phi;dN/#phi", 200, -4*TMath::Pi(), 4*TMath::Pi() );
    fHistPhiArtificialEvPlane = new TH1F("fHistPhiArtificialEvPlane", "#phi distribution of event plane;#phi;dN/#phi", 200, -4*TMath::Pi(), 4*TMath::Pi() );
    if ( fUsePhiShufflingByHand )
    {
        fOutList->Add(fHistPhiArtificialProfilingCheck);
        fOutList->Add(fHistPhiArtificialProfilingCheckWrtEvPlane);
        fOutList->Add(fHistPhiArtificialEvPlane);
    }
    
    //Eta vs Zv
    fHistEtaVsZvCoverage = new TH2D("fHistEtaVsZvCoverage","TPC tracks #eta vs Zv;V_{z} (cm);#eta",100,-20,20,50,-2,2);
    fHistEtaVsZvCoverageAccepted = new TH2D("fHistEtaVsZvCoverageAccepted","Accepted TPC tracks #eta vs Zv;V_{z} (cm);#eta",100,-20,20,50,-2,2);
    fOutList->Add(fHistEtaVsZvCoverage);
    fOutList->Add(fHistEtaVsZvCoverageAccepted);
    
    int lMaxAcceptedTracksInHist = fIsIonsAnalysis ? 4001 : 101;
    fHistMultBeforeCuts = new TH1D("fHistMultBeforeCuts","N_{ch} - tracks;N_{ch} tracks;Entries"
                                   ,lMaxAcceptedTracksInHist,-0.5,lMaxAcceptedTracksInHist-0.5);
    fOutList->Add(fHistMultBeforeCuts);
    
    fHistAcceptedMult = new TH1D("fHistAcceptedMult","N_{ch} - accepted tracks;N_{ch} accepted;Entries"
                                 ,lMaxAcceptedTracksInHist,-0.5,lMaxAcceptedTracksInHist-0.5);
    fOutList->Add(fHistAcceptedMult);
    
    fHistAcceptedTracks = new TH1D("fHistAcceptedTracks","N_{ch} - accepted tracks for LRC;N_{ch} accepted;Entries"
                                   ,lMaxAcceptedTracksInHist,-0.5,lMaxAcceptedTracksInHist-0.5);
    fOutList->Add(fHistAcceptedTracks);
    
    fHistMultiplicityInEtaRegion = new TH1D("fHistMultiplicityInEtaRegion","N_{ch} in |#eta| region;N_{ch};Entries"
                                            ,lMaxAcceptedTracksInHist,-0.5,lMaxAcceptedTracksInHist-0.5);
    fOutList->Add(fHistMultiplicityInEtaRegion);
    
    fHistMultiplicityInEtaRegionAfterPtCuts = new TH1D("fHistMultiplicityInEtaRegionAfterPtCuts","N_{ch} in |#eta| region after p_{T} cuts;N_{ch};Entries"
                                                       ,lMaxAcceptedTracksInHist,-0.5,lMaxAcceptedTracksInHist-0.5);
    fOutList->Add(fHistMultiplicityInEtaRegionAfterPtCuts);
    
    if ( fAnalyseMCESD )
    {
        fHist2DMultiplicityMCESDInEtaRegion = new TH2D("fHist2DMultiplicityMCESDInEtaRegion","N_{ch} in |#eta| region;N_{ch} MC;N_{ch} ESD;Entries"
                                                       ,101,-0.5,100.5,101,-0.5,100.5);
        fOutList->Add(fHist2DMultiplicityMCESDInEtaRegion);
    }
    
    
    fHistAcceptedTracksAfterPtCuts = new TH1D("fHistAcceptedTracksAfterPtCuts","N_{ch} - accepted tracks for LRC after Pt cuts;N_{ch} accepted;Entries"
                                              ,lMaxAcceptedTracksInHist,-0.5,lMaxAcceptedTracksInHist-0.5);
    fOutList->Add(fHistAcceptedTracksAfterPtCuts);
    
    fHistAcceptedTPCtracks = new TH1D("fHistAcceptedTPCtracks","N_{ch} - accepted tracks with TPC inner param;N_{ch} accepted;Entries"
                                      ,lMaxAcceptedTracksInHist,-0.5,lMaxAcceptedTracksInHist-0.5);
    fOutList->Add(fHistAcceptedTPCtracks);
    
    fHistClustersTPC = new TH1D("fHistClustersTPC","N Clusters TPC;N_{TPC clusters};Entries",161,-0.5,160.5);
    fOutList->Add(fHistClustersTPC);
    
    fHistClustersTPCafterCuts = new TH1D("fHistClustersTPCafterCuts","N Clusters TPC after cuts;N_TPC_clusters_after_cuts;Entries",161,-0.5,160.5);
    fOutList->Add(fHistClustersTPCafterCuts);
    
    
    fHistCrossedRowsTPC = new TH1D("fHistCrossedRowsTPC","N Crossed Rows TPC;N_{TPC CrossedRows};Entries",161,-0.5,160.5);
    fOutList->Add(fHistCrossedRowsTPC);
    
    fHistCrossedRowsTPCafterCuts = new TH1D("fHistCrossedRowsTPCafterCuts","N CrossedRows TPC after cuts;N_TPC_clusters_after_cuts;Entries",161,-0.5,160.5);
    fOutList->Add(fHistCrossedRowsTPCafterCuts);
    
    
    
    fHistTrackletsITS = new TH1D("fHistTrackletsITS","N Tracklets ITS;N_ITS_tracklets;Entries",101,-0.5,100.5);
    fOutList->Add(fHistTrackletsITS);
    
    
    fHistClustersITS = new TH1D("fHistClustersITS","N Clusters ITS;N_ITS_clusters;Entries",11,-0.5,10.5);
    fOutList->Add(fHistClustersITS);
    
    
    fHist2DClustersTPCvsPt = new TH2D("fHist2DClustersTPCvsPt","Num TPC clusters vs Pt;P_t;N_clusters",50,0,2,161,-0.5,160.5);
    fOutList->Add(fHist2DClustersTPCvsPt);
    
    fHist2DClustersTPCvsEta = new TH2D("fHist2DClustersTPCvsEta","Num TPC clusters vs Eta;#eta;N_clusters",50,-1.5,1.5,161,-0.5,160.5);
    fOutList->Add(fHist2DClustersTPCvsEta);
    
    fHist2DAcceptedTracksPtvsEta = new TH2D("fHist2DAcceptedTracksPtvsEta","Accepted Tracks Pt vs Eta;pt;#eta",50,0,2,50,-1.5,1.5);
    fOutList->Add(fHist2DAcceptedTracksPtvsEta);
    
    fHistProbabilitiesPID = new TH1D("fHistProbabilitiesPID","PID Probabilities;pid;Entries",10,-0.5,9.5 );//AliPID::kSPECIES,-0.5,AliPID::kSPECIES-0.5);
    fOutList->Add(fHistProbabilitiesPID);
    
    fHistESDtrackMass = new TH1D("fHistESDtrackMass","ESD Mass;Mass;Entries",100,0,1);
    fOutList->Add(fHistESDtrackMass);
    
    fHistProbabilityPion = new TH1D("fHistProbabilityPion","Probability Pion;Probability;Entries",100,0,1);
    fOutList->Add(fHistProbabilityPion);
    
    fHistProbabilityKaon = new TH1D("fHistProbabilityKaon","Probability Kaon;Probability;Entries",100,0,1);
    fOutList->Add(fHistProbabilityKaon);
    
    fHistProbabilityProton = new TH1D("fHistProbabilityProton","Probability Proton;Probability;Entries",100,0,1);
    fOutList->Add(fHistProbabilityProton);
    
    
    
    fHistParticlesDistrBeforeCuts = new TH1D("fHistParticlesDistrBeforeCuts","Particles Distr;Particle;Entries",5,-0.5,4.5);
    TString gBinParticleNames[5] = {/*"Other",*/"Electron","Muon","Pion","Kaon", "Proton"};
    for(Int_t i = 1; i <= 5; i++)fHistParticlesDistrBeforeCuts->GetXaxis()->SetBinLabel(i,gBinParticleNames[i-1].Data());
    fOutList->Add(fHistParticlesDistrBeforeCuts);
    
    fHistParticlesDistr = new TH1D("fHistParticlesDistr","Particles Distr;Particle;Entries",5,-0.5,4.5);
    //TString gBinParticleNames[6] = {"Undefined","Electron","Muon","Pion","Kaon", "Proton"};
    for(Int_t i = 1; i <= 5; i++)fHistParticlesDistr->GetXaxis()->SetBinLabel(i,gBinParticleNames[i-1].Data());
    fOutList->Add(fHistParticlesDistr);
    
    fHistCentralityPercentile = new TH1D("fHistCentralityPercentile","Centrality Percentile;Centrality;Entries",101,-1.01,100.01);
    if ( fIsIonsAnalysis ) fOutList->Add(fHistCentralityPercentile);
    
    fHistCentralityClass10 = new TH1D("fHistCentralityClass10","Centrality Class 10;Centrality;Entries",10,-0.5,9.5);
    if ( fIsIonsAnalysis ) fOutList->Add(fHistCentralityClass10);
    
    fHistCentralityClass5 = new TH1D("fHistCentralityClass5","Centrality Class 5;Centrality;Entries",20,-0.5,19.5);
    if ( fIsIonsAnalysis ) fOutList->Add(fHistCentralityClass5);
    
    
    fRand = new TRandom3;
    //ZDC additional study
    
    fMultForZDCstudy[0] = 0;
    fMultForZDCstudy[1] = 200;
    fMultForZDCstudy[2] = 220;
    fMultForZDCstudy[3] = 240;
    fMultForZDCstudy[4] = 250;
    
    TString strZDCenergyName;
    TString strZDCenergyTitle;
    for ( int i = 0; i < 5; i++ )
    {
        strZDCenergyName = Form( "fHistZDCenergy%d", fMultForZDCstudy[i] );
        strZDCenergyTitle = Form(  "ZDC Energy, multThreshold %d;Energy;Entries", fMultForZDCstudy[i] );
        fHistZDCenergy[i] = new TH1D(strZDCenergyName, strZDCenergyTitle, 200, 0, 10000);
        if ( fFlagWatchZDC && fIsIonsAnalysis )
            fOutList->Add(fHistZDCenergy[i]);
    }
    
    
    fHistZDCparticipants = new TH1D("fHistZDCparticipants","ZDC Participants;Participants;Entries",800,0,800);
    if ( fFlagWatchZDC && fIsIonsAnalysis )
        fOutList->Add(fHistZDCparticipants);
    
    if ( fFlagWatchV0 )//fIsIonsAnalysis )
    {
        const int nBinsV0multForNonIons = 700;
        const int nBinsV0multForIons = 220;
        const int nMaxV0multForIons = 22000;
        fHistV0multiplicity = new TH1D("fHistV0multiplicity","V0 Multiplicity;Multiplicity;Entries",
                                       fIsIonsAnalysis ? nBinsV0multForIons : nBinsV0multForNonIons, 0, fIsIonsAnalysis ? nMaxV0multForIons : nBinsV0multForNonIons);
        fHistV0Amultiplicity = new TH1D("fHistV0Amultiplicity","V0-A Multiplicity;Multiplicity;Entries",
                                        fIsIonsAnalysis ? nBinsV0multForIons : nBinsV0multForNonIons, 0, fIsIonsAnalysis ? nMaxV0multForIons : nBinsV0multForNonIons);
        fHistV0Cmultiplicity = new TH1D("fHistV0Cmultiplicity","V0-C Multiplicity;Multiplicity;Entries",
                                        fIsIonsAnalysis ? nBinsV0multForIons : nBinsV0multForNonIons, 0, fIsIonsAnalysis ? nMaxV0multForIons : nBinsV0multForNonIons);
        fHist2DV0ACmultiplicity = new TH2D("fHist2DV0ACmultiplicity","V0 A-C Multiplicity;Multiplicity A;Multiplicity C;Entries"
                                           ,50, 0, fIsIonsAnalysis ? nMaxV0multForIons/2 : nBinsV0multForNonIons/2
                                                                     ,50, 0, fIsIonsAnalysis ? nMaxV0multForIons/2 : nBinsV0multForNonIons/2);
        
        fHist2DTracksAcceptedVsV0multiplicity = new TH2D("fHist2DTracksAcceptedVsV0multiplicity","V0 Multiplicity vs Accepted;N Accepted tracks;Multiplicity V0;Entries"
                                                         ,70, 0, fIsIonsAnalysis ? 3000 : 70
                                                                                   ,50, 0, fIsIonsAnalysis ? nMaxV0multForIons/2 : 100 );
        
        //number of fired cells
        const int nV0cells = 32;
        //fHistV0cells = new TH1D("fHistV0cells","V0 cells;N cells;Entries", nV0cells+1, -0.5, nV0cells + 0.5 );
        fHistV0Acells = new TH1D("fHistV0Acells","V0-A cells;N cells;Entries", nV0cells+1, -0.5, nV0cells + 0.5 );
        fHistV0Ccells = new TH1D("fHistV0Ccells","V0-C cells;N cells;Entries", nV0cells+1, -0.5, nV0cells + 0.5 );
        fHist2DV0ACcells = new TH2D("fHist2DV0ACcells","V0 A-C cells;N cells A;N cells C;Entries", nV0cells+1, -0.5, nV0cells + 0.5, nV0cells+1, -0.5, nV0cells + 0.5 );
        
        
        fOutList->Add(fHistV0multiplicity);
        fOutList->Add(fHistV0Amultiplicity);
        fOutList->Add(fHistV0Cmultiplicity);
        fOutList->Add(fHist2DV0ACmultiplicity);
        fOutList->Add(fHist2DTracksAcceptedVsV0multiplicity);
        //cells
        //fOutList->Add(fHistV0cells    );
        fOutList->Add(fHistV0Acells   );
        fOutList->Add(fHistV0Ccells   );
        fOutList->Add(fHist2DV0ACcells);
        
        //mult in rings
        //        char strV0ringName[200];
        //        char strV0ringTitle[200];
        TString strV0ringName;
        TString strV0ringTitle;
        for ( int i = 0; i < 4; i++ )
        {
            strV0ringName =  Form( "fHistV0AmultiplicityRing%d", i );
            strV0ringTitle = Form( strV0ringTitle, "V0-A Multiplicity Ring %d;Multiplicity;Entries", i );
            fHistV0AmultiplicityRing[i] = new TH1D( strV0ringName, strV0ringTitle,
                                                    fIsIonsAnalysis ? nBinsV0multForIons : nBinsV0multForNonIons, 0, fIsIonsAnalysis ? nMaxV0multForIons : nBinsV0multForNonIons);
            strV0ringName =  Form(  "fHistV0CmultiplicityRing%d", i );
            strV0ringTitle = Form( strV0ringTitle, "V0-C Multiplicity Ring %d;Multiplicity;Entries", i );
            fHistV0CmultiplicityRing[i] = new TH1D( strV0ringName, strV0ringTitle,
                                                    fIsIonsAnalysis ? nBinsV0multForIons : nBinsV0multForNonIons, 0, fIsIonsAnalysis ? nMaxV0multForIons : nBinsV0multForNonIons);
            strV0ringName =  Form( "fHist2DV0ACmultiplicityRing%d", i );
            strV0ringTitle = Form(  "V0-AC Multiplicity Ring %d;Multiplicity A;Multiplicity C;Entries", i );
            fHist2DV0ACmultiplicityRing[i] = new TH2D( strV0ringName, strV0ringTitle
                                                       , 100, 0, fIsIonsAnalysis ? nMaxV0multForIons/2 : 100
                                                                                   , 100, 0, fIsIonsAnalysis ? nMaxV0multForIons/2 : 100 );
            fOutList->Add(fHistV0AmultiplicityRing[i]);
            fOutList->Add(fHistV0CmultiplicityRing[i]);
            fOutList->Add(fHist2DV0ACmultiplicityRing[i]);
            
            //mult in barrel vs V0 rings
            strV0ringName =  Form( "fHist2DTracksAcceptedVsV0AmultiplicityRing%d", i );
            strV0ringTitle = Form(  "Accepted tracks vs V0-A Multiplicity in Ring %d;N Accepted tracks;Multiplicity V0A;Entries", i );
            fHist2DTracksAcceptedVsV0AmultiplicityRing[i] = new TH2D( strV0ringName, strV0ringTitle
                                                                      , 100, 0, fIsIonsAnalysis ? nMaxV0multForIons/2 : 100
                                                                                                  , 100, 0, fIsIonsAnalysis ? nMaxV0multForIons/2 : 100 );
            strV0ringName =  Form( "fHist2DTracksAcceptedVsV0CmultiplicityRing%d", i );
            strV0ringTitle = Form(  "Accepted tracks vs V0-C Multiplicity in Ring %d;N Accepted tracks;Multiplicity V0C;Entries", i );
            fHist2DTracksAcceptedVsV0CmultiplicityRing[i] = new TH2D( strV0ringName, strV0ringTitle
                                                                      , 100, 0, fIsIonsAnalysis ? nMaxV0multForIons/2 : 100
                                                                                                  , 100, 0, fIsIonsAnalysis ? nMaxV0multForIons/2 : 100 );
            
            fOutList->Add(fHist2DTracksAcceptedVsV0AmultiplicityRing[i]);
            fOutList->Add(fHist2DTracksAcceptedVsV0CmultiplicityRing[i]);
        }
        
    }
    
    
    
    
    //    fHistV0spectra = new TH1D("fHistV0spectra","V0 spectra;Mass, GeV;Entries",500,0,5000);
    //    fOutList->Add(fHistV0spectra);
    
    
    fHistMClabels = new TH1D("fHistMClabels","MC label;label;Entries",102,-100.5,100.5);
    fOutList->Add(fHistMClabels);
    
    fHistRejectedTracksCharge = new TH1D("fHistRejectedTracksCharge","Rejected tracks charge;charge;Entries",3,-1.5,1.5);
    fOutList->Add(fHistRejectedTracksCharge);
    
    fHistTracksCharge = new TH1D("fHistTracksCharge","Accepted tracks charge;charge;Entries",3,-1.5,1.5);
    fOutList->Add(fHistTracksCharge);

    // ##### net charge vs pt study
    const int kPtNetChargePtBins = 200;
    const double kPtNetChargePtMin = 0.1;
    const double kPtNetChargePtMax = 2.1;

    //pt plus
    fHistPtPlus = new TH1D("fHistPtPlus","p_{T} +;p_{T};dN/dpT"
                                  , kPtNetChargePtBins,  kPtNetChargePtMin, kPtNetChargePtMax
                                  );
    fOutList->Add(fHistPtPlus);

    //pt minus
    fHistPtMinus = new TH1D("fHistPtMinus","p_{T} -;p_{T};dN/dpT"
                                  , kPtNetChargePtBins,  kPtNetChargePtMin, kPtNetChargePtMax
                                  );
    fOutList->Add(fHistPtMinus);

    //net charge vs pT
    fHistNetChargeVsPt = new TH1D("fHistNetChargeVsPt","charge vs p_{T};p_{T};Q"
                                  , kPtNetChargePtBins,  kPtNetChargePtMin, kPtNetChargePtMax
                                  );
    fOutList->Add(fHistNetChargeVsPt);

    //plus
    fHistChargePlusVsPtTmp = new TH1D( "fHistChargePlusVsPtTmp","charge plus vs p_{T};p_{T};q plus"
                                       , kPtNetChargePtBins,  kPtNetChargePtMin, kPtNetChargePtMax
                                       );
    //    fOutList->Add(fHistChargePlusVsPtTmp);

    //minus
    fHistChargeMinusVsPtTmp = new TH1D( "fHistChargeMinusVsPtTmp","charge minus vs p_{T};p_{T};q minus"
                                        , kPtNetChargePtBins,  kPtNetChargePtMin, kPtNetChargePtMax
                                        );
    //    fOutList->Add(fHistChargeMinusVsPtTmp);

    //NetChargeVsPt
    fHist2DNetChargeVsPt = new TH2D( "fHist2DNetChargeVsPt","Net charge vs p_{T};p_{T};Q"
                                     , kPtNetChargePtBins,  kPtNetChargePtMin, kPtNetChargePtMax
                                     , 40,  -20, 20
                                     );
    fOutList->Add(fHist2DNetChargeVsPt);

    //NetChargeVsPt CorrectedOnEventMean
    fHist2DNetChargeVsPtCorrectedOnEventMean = new TH2D( "fHist2DNetChargeVsPtCorrectedOnEventMean","Net charge corrected on mean e-by-e vs p_{T};p_{T};Q"
                                                         , kPtNetChargePtBins,  kPtNetChargePtMin, kPtNetChargePtMax
                                                         , 40,  -20, 20
                                                         );
    fOutList->Add(fHist2DNetChargeVsPtCorrectedOnEventMean);

    //NetChargeVsPt CorrectedOnEventMean normalized on Nch
    fHist2DNetChargeVsPtCorrectedOnEventMeanNormOnNch = new TH2D( "fHist2DNetChargeVsPtCorrectedOnEventMeanNormOnNch","Net charge vs p_{T} corrected on mean e-by-e normalized on nCh;p_{T};Q"
                                                                  , kPtNetChargePtBins,  kPtNetChargePtMin, kPtNetChargePtMax
                                                                  , 40,  -20, 20
                                                                  );
    fOutList->Add(fHist2DNetChargeVsPtCorrectedOnEventMeanNormOnNch);



    if ( fArtificialInefficiency >= 0 ) //i.e. have this kind of analysis
    {
        fHistNumberOfDroppedByHandTracks = new TH2D("fHistNumberOfDroppedByHandTracks","Accepted tracks vs Dropped artificially;N_{ch} accepted;N_{ch} dropped;Entries", 71, -0.5, 70.5,   71, -0.5, 70.5);
        fOutList->Add(fHistNumberOfDroppedByHandTracks);
    }
    
    
    // Returning TH1,TH2 AddDirectory status
    TH1::AddDirectory(lTH1oldStatus);
    TH2::AddDirectory(lTH2oldStatus);
    
    //fHistProbabilitiesPID->Fill(1,fStrPIDforFwd);
    //fHistProbabilitiesPID->Fill(2,fStrPIDforBwd);
    Printf("UserCreateOutputObjects done.");
    
    // NEW HISTO added to fOutput here
    PostData(1, fOutList); // Post data for ALL output slots >0 here, to get at least an empty histogram
    
    for( Int_t i = 0; i < fLRCproc.GetEntries(); i++)
    {
        PostData( Proc(i)->GetOutputSlotNumber(),Proc(i)->CreateOutput() );
    }
    
    
    //if ( fSetIncludeEventTreeInOutput )
    //    PostData(2, fEventTree);
    //int a;
    //cin >> a;
    
    
}

//________________________________________________________________________
//void AliAnalysisTaskLRC::UserExec(Option_t *)
//{

//    //cout << "TEST in Task" << endl;
//    // ###### tuning phi

//    //rotate phi when N_phi_sectors > 1 and call UserExecLoop:
//    double phiStep = 2 * TMath::Pi();
//    if ( fNumberOfSectors > 1 )
//        phiStep /= fNumberOfSectors;

//    double lPhiRotatedExtra = 0; //additional phi rotation of all tracks
//    for ( Int_t sectorId = 0; sectorId < fNumberOfSectors; sectorId++ )
//    {
//        //cout << "loop " << sectorId << endl;
//        UserExecLoop( lPhiRotatedExtra );
//        lPhiRotatedExtra += phiStep; //rotate track
//    }


//}

//________________________________________________________________________
void AliAnalysisTaskLRC::UserExec(Option_t *)   //UserExecLoop( Double_t phiAdditional )//Option_t *)
{
    // ########### if use toy events
    if ( fUseToyEvents )
    {
        //generate all events here and fill LRC processors
        UseToyEvents();

        //just return
        return;
    }
    
    // Main loop
    // Called for each event
    //printf( "starting UserExec...\n" );
    //Pid setting
    if ( fPIDsensingFlag )
    {
        SetParticleTypeToProcessors( 0 /*fwd*/, fStrPIDforFwd );
        SetParticleTypeToProcessors( 1 /*backward*/, fStrPIDforBwd );
        fPIDsensingFlag = kFALSE;
    }
    
    //printf( "starting event...\n" );
    //The common PID object can then be retrieved from the input handler:
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    //if (!fPIDResponse) !!!!!!!!!!!!
    //	AliFatal("This Task needs the PID response attached to the inputHandler");
    
    
    AliVEvent *event = InputEvent();
    AliESDEvent *lESD = 0x0;
    AliAODEvent *lAOD = 0x0;
    
    AliStack *stack = 0x0;
    AliMCEvent *eventMC = 0x0;
    
    
    if ( fAnalysisLevel == "ESD" ) //fAnalysisType == en_AnalysisType_AOD )
    {
        lESD = dynamic_cast<AliESDEvent*> (event) ;
        //printf( "cast to ESD is Ok\n." );
    }
    else if ( fAnalysisLevel == "AOD" ) //fAnalysisType == en_AnalysisType_AOD )
    {
        //cout << "test AOD analysis... " << endl;
        lAOD = dynamic_cast<AliAODEvent*>(event); // from TaskSE
        //Int_t bMagneticFieldSign = (lAOD->GetMagneticField() > 0) ? 1 : -1;
        //printf( "Number of AOD tracks = %d, magnetic field = %f\n", lAOD->GetNumberOfTracks(), lAOD->GetMagneticField() );
        //return;
        
    }
    
    if( fRunKine ) //|| fAnalyseMCESD )
    {
        eventMC = MCEvent();
        stack = eventMC->Stack();
        //Printf("Number of primaries: %d",stack->GetNprimary());
    }
    
    if (!event)
    {
        //Printf("ERROR: Could not retrieve event");
        return;
    }
    //cout << "test start Event!" << endl;
    
    if( (lESD ) && ( !fEsdTrackCuts ) )
    {
        AliDebug(AliLog::kError, "No ESD track cuts avalible");
        return;
    }
    
    //make combined flag for esd selection during MC kine analysis
    Bool_t lRunKine = ( fRunKine && !fAnalyseMCESD );
    
    // Processing event selection
    fHistEventCutStats->Fill( "Total", 1 );
    
    
    
    //cout << "test: nOfTracks = " << event->GetNumberOfTracks() << endl;
    
    //Trigger
    Bool_t lTrigger = kTRUE;
    if( lESD && !lRunKine ) //just for tests; usually it is done before UserExecLoop in PhysSel task
    {
        if(fShowEventStats)
            Printf("Trigger classes: %s:", lESD->GetFiredTriggerClasses().Data());
        
        lTrigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
        
        //8 TeV runs: TEST triggers! (29.12.2012)
        //        Bool_t lTrigger1 = (lESD->IsTriggerClassFired("CINT7-B-NOPF-ALLNOTRD")) ? 1 : 0;
        //        Bool_t lTrigger2 = (lESD->IsTriggerClassFired("CINT7WU-B-NOPF-ALL")) ? 1 : 0;
        //        lTrigger = lTrigger1 && lTrigger2;
        
        if( !lTrigger )
        {
            if( fShowEventStats )
                Printf("Rejected!");
            fHistEventCutStats->Fill( "No trigger", 1 );
            PostData(1, fOutList);
            return;
        }
    }
    //Centrality      //https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies
    Double_t lCentralityPercentile 	= 0.;
    Int_t lCentralityClass10 		= 0;
    Int_t lCentralityClass5 		= 0;
    Float_t lReactionPlane       = -1.;
    
    
    if( !lRunKine && fIsIonsAnalysis )
    {
        AliCentrality *centrality = 0x0;
        
        if ( lESD )
            centrality = lESD->GetCentrality();
        else if ( lAOD )
            centrality = lAOD->GetCentrality();
        
        if ( centrality )
        {
            lCentralityPercentile 	= centrality->GetCentralityPercentile("V0M"); 	// returns the centrality percentile, a float from 0 to 100 (or to the trigger efficiency)
            lCentralityClass10 		= centrality->GetCentralityClass10("V0M");//"FMD"); 		// returns centrality class for 10% width (a integer from 0 to 10)
            lCentralityClass5 		= centrality->GetCentralityClass5("V0M");//"TKL"); 		// returns centrality class for 5% width (a integer from 0 to 20)
            Bool_t lCentralityInClass = centrality->IsEventInCentralityClass(fMinCentralityClass,fMaxCentralityClass,"V0M"); // returns kTRUE if the centrality of the event is between a and b, otherwise kFALSE
            //cout << lCentralityPercentile << " "
            //<< fMinCentralityClass << " "
            //<< fMaxCentralityClass << " "
            //<< lCentralityInClass << endl;
            
            if ( !lCentralityInClass )
            {
                //cout << "outside of centrality class!" << endl;
                fHistEventCutStats->Fill("Wrong centrality", 1);
                PostData(1, fOutList);
                return;
            }
            
            fHistCentralityPercentile->Fill( lCentralityPercentile );
            fHistCentralityClass10->Fill( lCentralityClass10 );
            fHistCentralityClass5->Fill( lCentralityClass5 );
            
            // get the reaction plane
            lReactionPlane = GetEventPlane( event );
            fHistEventPlane->Fill( lReactionPlane, lCentralityPercentile );
        }
    }
    
    //number of verteces in ESD (pileup)
    if ( lESD )
    {
        int lNumberOfPileUpVerteces = 0;
        TClonesArray *lPileupVertecesTracks = lESD->GetPileupVerticesTracks();
        lNumberOfPileUpVerteces = lPileupVertecesTracks->GetSize();
        fHistNumberOfPileupVerticesTracks->Fill( lNumberOfPileUpVerteces );
        
        TClonesArray *lPileupVertecesSPD = lESD->GetPileupVerticesSPD();
        lNumberOfPileUpVerteces = lPileupVertecesSPD->GetSize();
        fHistNumberOfPileupVerticesSPD->Fill( lNumberOfPileUpVerteces );
    }
    else if ( lAOD ) //number of verteces in AOD (pileup)
    {
        int lNumberOfPileUpVerteces = 0;
        lNumberOfPileUpVerteces = lAOD->GetNumberOfPileupVerticesTracks();
        fHistNumberOfPileupVerticesTracks->Fill( lNumberOfPileUpVerteces );
        lNumberOfPileUpVerteces = lAOD->GetNumberOfPileupVerticesSPD();
        fHistNumberOfPileupVerticesSPD->Fill( lNumberOfPileUpVerteces );
    }
    
    
    // Vertex present
    //    const AliESDVertex *vertex = lESD->GetPrimaryVertex();
    Double_t lVertexX = -1000;
    Double_t lVertexY = -1000;
    Double_t lVertexZ = -1000;
    Double_t lMCVertexX = -1000;
    Double_t lMCVertexY = -1000;
    Double_t lMCVertexZ = -1000;
    
    Bool_t lReconstructedVertexPresent = kFALSE;
    Bool_t lVertexAcceptable = kFALSE;
    
    if ( fAnalysisLevel == "ESD" ||  fAnalysisLevel == "AOD" )
    {
        const AliVVertex *vertex = event->GetPrimaryVertex();
        fHistVertexNconributors->Fill( vertex->GetNContributors() );
        
        Double32_t lCov[6];
        vertex->GetCovarianceMatrix(lCov);
        
        //check nContributors and z*z>0
        lReconstructedVertexPresent = ( ( vertex->GetNContributors() > 0 ) && ( lCov[5] != 0 ) );
        //        if ( (vertex) && (vertex->GetNContributors() > 0 ) && ( vertex->GetZRes() == 0 ) )
        //        {
        //            int aa;
        //            cout << "stop with nContributors = " << vertex->GetNContributors() << endl;
        //            cin >> aa;
        //        }
        lVertexX = vertex->GetX();
        lVertexY = vertex->GetY();
        lVertexZ = vertex->GetZ();
    }
    
    if ( fRunKine ) //take MC vertex params and rewrite position
    {
        AliGenEventHeader *header = eventMC->GenEventHeader();
        if(!header) return;
        
        TArrayF gVertexArray;
        header->PrimaryVertex(gVertexArray);
        lMCVertexX = gVertexArray.At(0);
        lMCVertexY = gVertexArray.At(1);
        lMCVertexZ = gVertexArray.At(2);
        
        //fill mc-reco diff
        if ( lReconstructedVertexPresent )
        {
            fHistVxMCrecoDiff->Fill( lMCVertexX - lVertexX );
            fHistVyMCrecoDiff->Fill( lMCVertexY - lVertexY );
            fHistVzMCrecoDiff->Fill( lMCVertexZ - lVertexZ );
        }
        
    }
    
    //cut on reco vertex params
    if ( !lRunKine )
    {
        if( ( !lReconstructedVertexPresent ) && fCheckForkVtx )
        {
            if( fShowEventStats )
                Printf("Bad vertex params");
            fHistEventCutStats->Fill("Bad vertex params",1);
            PostData(1, fOutList);
            return;
        }
    }
    else
    {
        //rewrite vert x,y,z into mc values
        lVertexX = lMCVertexX;
        lVertexY = lMCVertexY;
        lVertexZ = lMCVertexZ;
    }
    
    
    // Vertex in range
    lVertexAcceptable = (TMath::Abs(lVertexX) < fVxMax) && (TMath::Abs(lVertexY) < fVyMax);
    if( lVertexAcceptable )
    {
        if( fVzMax > 0 )   //   fVzMax < 0 -> select Zv outside selected range
            lVertexAcceptable = (TMath::Abs(lVertexZ) < fVzMax);
        else
            lVertexAcceptable = (TMath::Abs(lVertexZ) > -fVzMax);
    }
    
    if( (!lVertexAcceptable) && fCheckForVtxPosition && !lRunKine )
    {
        if(fShowEventStats)
            Printf("Vertex out of range");
        fHistEventCutStats->Fill("Bad vertex position",1);
        PostData(1, fOutList);
        return;
    }
    
    //fill bin if no reco vertex in MC analysis
    if ( lRunKine && !lReconstructedVertexPresent )
        fHistEventCutStats->Fill( "No reco vertex", 1 );
    
    fHistVx->Fill(lVertexX);
    fHistVy->Fill(lVertexY);
    fHistVz->Fill(lVertexZ);
    //end with vertex
    
    
    
    
    //cut on number of SPD tracklets (25.03.2012)
    if( lESD && !lRunKine )
    {
        //How to get tracklets
        const AliMultiplicity *tracklets = lESD->GetMultiplicity();
        Int_t multSPD = tracklets->GetNumberOfTracklets();
        //Int_t nITStracklets = kalmanTrack->GetNumberOfTracklets();
        fHistTrackletsITS->Fill( multSPD );
        if ( multSPD < fMinNumberOfSPDtracklets )
        {
            if(fShowEventStats)
                Printf("Too few SPD tracklets");
            fHistEventCutStats->Fill( "Few SPD tracklets", 1 );
            PostData(1, fOutList);
            return;
        }
    }
    
    
    //tmp eta ranges (for ZDC respond study)
    double etaBmin = -0.8;
    double etaBmax = -0.6;
    double etaFmin = 0.6;
    double etaFmax = 0.8;
    int countTracksEtaB = 0;
    int countTracksEtaF = 0;
    
    //n accepted tracks
    int lNchTrigger = 0;
    int lTPCtracks  = 0; // added 25.03.2012
    fHistMultBeforeCuts->Fill( event->GetNumberOfTracks() );
    
    // Pre event loop
    if( !lRunKine )
    {
        for ( Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++)
        {
            if ( fAnalysisLevel == "ESD" )
            {
                AliESDtrack *track = lESD->GetTrack(iTracks);
                if ( fEsdTrackCuts->AcceptTrack(track) )
                {
                    lNchTrigger++;
                    //todo: try to implement pt cuts here?....
                    //How to get TPC standalone tracks
                    AliExternalTrackParam *tpcTrack
                            = (AliExternalTrackParam *)track->GetTPCInnerParam();
                    if ( tpcTrack )
                        lTPCtracks++;
                }
                else
                {
                    fHistRejectedTracksCharge->Fill( track->Charge() );
                }
            }
            else if ( fAnalysisLevel == "AOD" )
            {
                AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
                if( aodTrack->TestFilterBit( fAODtrackCutBit ) )
                    lNchTrigger++;
            }
        }
    }
    fHistAcceptedTPCtracks->Fill( lTPCtracks );
    
    
    
    
    //if(lRunKine)lNchTrigger=eventMC->GetNumberOfPrimaries();
    
    // Nch bins cut
    if( ( lNchTrigger > fMaxAcceptedTracksCut ) && ( fMaxAcceptedTracksCut != 0 ) )
    {
        fHistEventCutStats->Fill( "HighMult cut", 1 );
        PostData(1, fOutList);
        return;
    }
    
    if( lNchTrigger < fMinAcceptedTracksCut )
    {
        fHistEventCutStats->Fill( "LowMult cut", 1 );
        PostData(1, fOutList);
        return;
    }
    fHistAcceptedMult->Fill(lNchTrigger);
    fHistEventCutStats->Fill("Analyzed",1);
    
    if ( lESD && fFlagWatchV0 ) // cut on V0 multiplicity "radius" in 2D-hist for both A and C sides
    {
        const AliESDVZERO* vzrData = lESD->GetVZEROData(); //aod the same
        
        //        const int lThresholdMultV0RingId = 3; //which ring is considered
        //        Double_t lThisEventV0MultSumRing3 = vzrData->GetMRingV0A(2) + vzrData->GetMRingV0C(2);
        //        Double_t lThisEventV0MultSumRing4 = vzrData->GetMRingV0A(3) + vzrData->GetMRingV0C(3);
        //        Double_t lThisEventV0MultSum = lThisEventV0MultSumRing3 + lThisEventV0MultSumRing4;
        double sumV0Amult = vzrData->GetMTotV0A();
        double sumV0Cmult = vzrData->GetMTotV0C();

        Double_t lThisEventV0MultSum = sumV0Amult + sumV0Cmult;

        //sqrt ( pow ( vzrData->GetMRingV0A(lThresholdMultV0RingId), 2 ) + pow ( vzrData->GetMRingV0C(lThresholdMultV0RingId), 2 ) );
        if ( lThisEventV0MultSum < fThresholdOnV0mult
             || sumV0Amult < fThresholdOnV0mult / 5 //additionally, to "exclude" diffraction
             || sumV0Cmult < fThresholdOnV0mult / 5 //additionally, to "exclude" diffraction
             )
        {
            PostData(1, fOutList);
            return;
        }
    }
    
    //    if ( fAnalysisLevel == "AOD" )
    //        return;
    

    //reset tmp histos
    fHistChargePlusVsPtTmp->Reset();
    fHistChargeMinusVsPtTmp->Reset();

    //Track selection counters
    int lNaccept=0;
    int lNacceptEtaInRegion=0;
    int lNacceptAfterPtCuts=0;
    int lNacceptEtaInRegionAfterPtCuts=0;
    //    int lNacceptPlusEtaInRegion=0;
    //    int lNacceptMinusEtaInRegion=0;
    int lPtOver=0;
    int lPtUnder=0;
    int lNoCharge=0;
    int lCutfail=0;
    int lDecay=0;
    //AliLRCBase::LRCparticleType lParticleType = AliLRCBase::kLRCany;
    
    int lNumberOfDropByHandTracks = 0; //number of artificially dropped tracks
    
    //Double_t probTPC[AliPID::kSPECIES]={0.};
    Double_t probTOF[AliPID::kSPECIES]={0.};
    Double_t probTPCTOF[AliPID::kSPECIES]={0.};
    
    Int_t nTracksInEvent = ( !lRunKine )
            ? event->GetNumberOfTracks()
            : eventMC->GetNumberOfPrimaries();//GetNumberOfTracks();
    
    //        if ( fSetIncludeEventTreeInOutput )
    //        {
    //            fSimpleEvent->SetHeader( fNsimpleEvents, -1, -1, lCentralityPercentile, lReactionPlane );
    //            fSimpleEvent->SetVertexPos( vertex->GetX(), vertex->GetY(), vertex->GetZ() );
    //            fNsimpleEvents++;
    //            fSimpleEvent->StartEventFilling();
    //        }
    //if ( !lRunKine )
    int lNumberOfAcceptedTracksForLRC = 0;
    // Track loop -----------------------------------------------------------------
    for (Int_t iTracks = 0; iTracks < nTracksInEvent/*event->GetNumberOfTracks()*/; iTracks++)
    {
        
        //Track variables
        double lPt = 0;   // Temp Pt
        double lEta = -100;	  // Temp ETA
        double lPhi = -100;    // Temp Phi
        double lMass = 0;    // Temp Mass
        double lEt = 0;
        Short_t lCharge = 0;
        
        AliVParticle* track = ( !lRunKine )
                ? event->GetTrack(iTracks)
                : eventMC->GetTrack(iTracks);
        
        
        if ( !track ) {
            Printf("ERROR: Could not receive track %d", iTracks);
            continue;
        }
        
        fHistTrackCutStats->Fill("Total",1);
        fHistEtaVsZvCoverage->Fill(lVertexZ,track->Eta());
        
        lPt = track->Pt();
        //if ( iTracks == 0 ) cout << "pt of 1st track is " << lPt << endl;
        lEta = track->Eta();
        lPhi = track->Phi();
        
        if ( !lRunKine && fAnalysisLevel == "AOD" )
        {
            AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
            if ( aodTrack->TestFilterBit( fAODtrackCutBit ) )
                fHistEtaAODpure->Fill( lEta );
        }
        
        //1.10.2012 - "inefficiency" in phi-acceptance "by hand"
        if ( lPhi > fPhiArtificialGapBegin &&
             lPhi < fPhiArtificialGapEnd )
            continue; // drop track in this "ineffective" area
        
        
        //lPhi += phiAdditional; // add here extra phi angle! (when applying tracks rotation by hand)
        
        
        //        if ( lPhi > 2 * TMath::Pi() )
        //            lPhi -= 2 * TMath::Pi();
        //cout << "track pt = " <<  lPt << endl;
        lCharge = track->Charge();
        
        Int_t lTriggerMask = -1; //for manual event tree extraction
        
        // ESD or AOD track cuts
        if( !lRunKine )
        {
            if ( fAnalysisLevel == "ESD" )
            {
                AliESDtrack* lESDtrack = lESD->GetTrack(iTracks);
                if( fShowPerTrackStats )
                    Printf("ESD Track N%d , Eta=%f, Pt=%f , TPC Nclusters =%d Sigma=%f ",  iTracks , lEta , lPt, lESDtrack-> GetTPCNcls(),fEsdTrackCuts->GetSigmaToVertex( lESDtrack) );
                //cluster histograms (for watching)
                fHistClustersTPC->Fill( lESDtrack->GetTPCNcls() );
                fHistCrossedRowsTPC->Fill( lESDtrack->GetTPCCrossedRows() );
                Int_t nITSclusters = lESDtrack->GetITSclusters(0);
                fHistClustersITS->Fill( nITSclusters );//kalmanTrack->GetNumberOfClusters()  );
                //Printf( " after \n" );
                
                fHist2DClustersTPCvsPt->Fill( lPt, lESDtrack->GetTPCNcls() );
                fHist2DClustersTPCvsEta->Fill( lEta, lESDtrack->GetTPCNcls() );
                
                //now check track cuts: either look at fEsdTrackCuts or take cuts from array and look
                if ( !fSwitchToListingCuts )
                {
                    if( !fEsdTrackCuts->AcceptTrack(lESDtrack) )
                    {
                        lCutfail++;
                        if( fShowPerTrackStats )
                            Printf("Rejected by cuts");
                        fHistTrackCutStats->Fill( "QA track cuts", 1 );
                        continue;
                    }
                    else
                    {
                        if(fShowPerTrackStats)
                            Printf("OK");
                    }
                }
                else
                {
                    lTriggerMask = 0;
                    //                    cout << "cuts: " << endl;
                    for( Int_t cutsId = 0; cutsId < fArrTrackCuts.GetEntries()/*fNumberOfCutsToRemember*/; cutsId++ )
                    {
                        //                        if ( !fArrTrackCuts[cutsId] )
                        //                            continue;
                        int cutsPassed = ( ((AliESDtrackCuts*)fArrTrackCuts[cutsId])->AcceptTrack( lESDtrack ) );
                        if ( cutsPassed )
                        {
                            fHistCutsNamesBins->Fill( cutsId );
                            Int_t cutMask = 1;
                            for( int iPow = 0; iPow < cutsId; iPow++ )
                                cutMask *= 2;
                            lTriggerMask ^=  cutMask;//(cutsId+1);
                        }
                        //                        cout << cutsPassed;
                    }
                }
                
                //fHist2DRejectedTracksPtvsEta->Fill( lPt, lEta );
                fHistClustersTPCafterCuts->Fill( lESDtrack->GetTPCNcls() );
                fHistCrossedRowsTPCafterCuts->Fill( lESDtrack->GetTPCCrossedRows() );
            }
            else if ( fAnalysisLevel == "AOD" )
            {
                AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
                if (!aodTrack)
                {
                    AliError(Form("Could not receive track %d", iTracks));
                    continue;
                }
                
                for(Int_t iTrackBit = 0; iTrackBit < 16; iTrackBit++){
                    fHistAODTrackStats->Fill(iTrackBit,aodTrack->TestFilterBit(1<<iTrackBit));
                    //                    cout << aodTrack->TestFilterBit(1<<iTrackBit) << " ";
                }
                //                cout << endl;
                if( !aodTrack->TestFilterBit( fAODtrackCutBit ) )
                {
                    //                    cout << "fAODtrackCutBit=" << fAODtrackCutBit << endl;
                    fHistTrackCutStats->Fill( "QA track cuts", 1 );
                    continue;
                }
                //                int a;
                //                cin >> a;
            }
            
            fHistTracksCharge->Fill(lCharge);

        }
        //        lEta = fRand->Uniform(-1,1); //!!!!!!!!!!!!!!!!!!!!!!! TEST!!!!!!
        // end of ESD or AOD track cuts
        
        
        //MC truth
        //Int_t label = -1000;
        TParticle * part = 0x0;
        if( lRunKine )   // Dropping undetectable particles in Kine
        {
            //            if ( fAnalyseMCESD ) //run MC+ESD
            //            {
            //                label = track->GetLabel();//TMath::Abs(track->GetLabel()); //abs!!!!
            //                fHistMClabels->Fill( label );
            //                if ( label < 0 ) //by Peter Hristov - it's ghost tracks
            //                    continue;
            //                part = stack->Particle(label);
            //            }
            //            else
            {
                AliMCParticle* trackMC = dynamic_cast<AliMCParticle *>(eventMC->GetTrack(iTracks));
                if(!trackMC)
                    continue;
                part = trackMC->Particle();
                if( !part ) continue;
            }
            
            Int_t lNoD = part->GetNDaughters();
            if(fShowPerTrackStats)
            {
                printf("Track %d, Mother %d, track Pt=%f, MC Pt=%f, PDG=%d  Nof Dothers=%d  ETA=%f   ## ",  iTracks , part->GetFirstMother() ,lPt,      part->Pt() ,part->GetPdgCode(),lNoD,lEta);
                printf("%s", part->GetName());
                printf(" ");
                printf("%s", part->GetPDG()->ParticleClass());
            }
            
            if( !stack->IsPhysicalPrimary(iTracks) )//label))
            {
                fHistTrackCutStats->Fill( "Not phys primary", 1 );
                continue;
            }
            
            //charge in TParticlePDG is in units of |e|/3
            fHistTracksCharge->Fill( part->GetPDG()->Charge() / 3 ); // 0-charge bin fills only for MC truth events
            //cout << " charge = " << part->GetPDG()->Charge();
            if ( part->GetPDG()->Charge() == 0 )
            {
                if(fShowPerTrackStats)
                    Printf(" ChargeReject");
                lNoCharge++;
                continue;
            }
            
            //cut on low-pt kine particles
            if ( lPt < fKineLowPtCut )
                continue;
        }    //End of  Kine particle filter
        
        //now decide that we have a good track
        lNaccept++;
        
        if( fabs(lEta) < fEtaRegionForTests ) //look at eta region
        {
            lNacceptEtaInRegion++;
            //net charge fillings
            fHistNetChargeVsPt->Fill( lPt, lCharge );
            if ( lCharge > 0 )
            {
                fHistChargePlusVsPtTmp->Fill( lPt );
                fHistPtPlus->Fill( lPt );
            }
            else if ( lCharge < 0 )
            {
                fHistChargeMinusVsPtTmp->Fill( lPt );
                fHistPtMinus->Fill( lPt );
            }
        }
        
        if( lPt > fMaxPtLimit )
        {
            lPtOver++;
            fHistTrackCutStats->Fill(2);
            continue;
        } // Dropping tracks with hi Pt
        
        if( lPt < fMinPtLimit )
        {
            lPtUnder++;
            fHistTrackCutStats->Fill(3);
            continue;
        } // Dropping tracks with low Pt
        lNacceptAfterPtCuts++;
        fHistTrackCutStats->Fill(4);
        fHistPt->Fill( lPt );
        fHistEta->Fill( lEta );
        fHistPhi->Fill( lPhi );
        
        if( fabs(lEta) < fEtaRegionForTests ) //look at eta region
            lNacceptEtaInRegionAfterPtCuts++;
        
        //fill counters for ZDC-comparing with a 'shelve'
        if ( lEta > etaBmin && lEta < etaBmax )
            countTracksEtaB++;
        if ( lEta > etaFmin && lEta < etaFmax )
            countTracksEtaF++;
        
        if ( fabs( lEta ) < 0.8 ) //eta-phi for tracks in "ALICE barrel acceptance"
            fHistEtaPhi->Fill(lEta, lPhi);
        
        
        // ###### New PID
        Int_t lMostProbablePIDPure = -1;
        Int_t lMostProbablePIDdirty = -1;
        if ( (lESD || lAOD) && !lRunKine && fPIDResponse )
        {
            AliVTrack *trackV = (AliVTrack*)event->GetTrack(iTracks);
            //##### from Pid Combined Task
            // compute priors for TPC+TOF, even if we ask just TOF for PID
            fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
            UInt_t detUsed = fPIDCombined->ComputeProbabilities(trackV, fPIDResponse, probTOF);
            Double_t priors[5]; 	// check priors used for TOF
            fPIDCombined->GetPriors(trackV,priors,fPIDResponse,detUsed);
            for(Int_t ispec=0;ispec<5;ispec++) fPriorsUsed[ispec]->Fill(TMath::Abs(trackV->Pt()),priors[ispec]);
            
            fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
            detUsed = fPIDCombined->ComputeProbabilities(trackV, fPIDResponse, probTPCTOF);
            fHistProbabilityPion->Fill( probTPCTOF[AliPID::kPion] );////probDensity[iParticle] );
            fHistProbabilityKaon->Fill( probTPCTOF[AliPID::kKaon] );
            fHistProbabilityProton->Fill( probTPCTOF[AliPID::kProton] );////probDensity[iParticle] );
            //cout << "chto-to s detId: " << detUsed << endl;
//            for ( Int_t ispec=0; ispec < (int)(AliPID::kSPECIES); ++ispec ) {
//                //cout << probTPCTOF[ispec] << "  " << endl;
//            }
            //cout << endl;
            if ( detUsed == (UInt_t)fPIDCombined->GetDetectorMask() )
            {
                Double_t lMaxProb = 0;
                for ( Int_t ispec=0; ispec < (int)(AliPID::kSPECIES); ++ispec )
                {
                    //Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(trackV,(AliPID::EParticleType)ispec);
                    fProbTPCTOF[ispec]->Fill( lPt, probTPCTOF[ispec] );
                    //fProbTPCTOFnSigmaTPC[ispec]->Fill(nSigmaTPC,probTPCTOF[ispec]);
                    //fProbTPCTOFnSigTPCMom[ibin][ispec]->Fill(nSigmaTPC,probTPCTOF[ispec]);
                    
                    if ( probTPCTOF[ispec] > lMaxProb )
                    {
                        lMaxProb = probTPCTOF[ispec];
                        lMostProbablePIDPure = ispec;
                    }
                }
                //fill for watching at max prob-s in pid species arrays
                fHistPidPureMaxProbability->Fill( lMaxProb );
            }
            
            const bool useDirtyPID = true;
            if ( useDirtyPID )
            {
                //Int_t lMostProbablePIDdirty = -1;
                
                Double_t lMaxProb = 0;
                for ( Int_t ispec=0; ispec < (int)(AliPID::kSPECIES); ++ispec )
                {
                    fProbAllDets[ispec]->Fill( lPt, probTPCTOF[ispec] );
                    if ( probTPCTOF[ispec] > lMaxProb )
                    {
                        lMaxProb = probTPCTOF[ispec];
                        lMostProbablePIDdirty = ispec;	// define most probable particle!!!
                    }
                }
                //fill for watching at max prob-s in pid species arrays
                fHistPidMaxProbability->Fill( lMaxProb );
                fHistParticlesDistr->Fill( lMostProbablePIDdirty );
            }
            
            
        }
        // end of new PID filling
        
        fHist2DAcceptedTracksPtvsEta->Fill( lPt, lEta );
        
        //get particle mass and Et
        if ( lMostProbablePIDPure >= 0 )
        {
            lMass = AliPID::ParticleMass( lMostProbablePIDPure );
            fHistESDtrackMass->Fill( lMass );
            lEt = sqrt( lPt*lPt + lMass*lMass );
            //cout << lEt << " " << lMass << " " << lPt << endl;
        }
        
        //eta vs vertex z coverage
        fHistEtaVsZvCoverageAccepted->Fill(lVertexZ,track->Eta());
        
        //artificial inefficiency (27.09.12) - check, what if we drop some fraction of tracks randomly
        Bool_t lDropDecision = ( ( 1.0 - fArtificialInefficiency ) < fRand->Uniform(0,1) );
        
        if ( lDropDecision )
        {
            lNumberOfDropByHandTracks++;
            continue;
        }
        
        if ( lNumberOfAcceptedTracksForLRC >= kMaxParticlesNumber )
        {
            AliDebug(AliLog::kError, "lNumberOfAcceptedTracksForLRC is too large...");
            break;
        }
        
        //fill arrays with track data for LRC
        fArrayTracksPt[lNumberOfAcceptedTracksForLRC]       = (fEtInsteadOfPt ? lEt : lPt);
        fArrayTracksEta[lNumberOfAcceptedTracksForLRC]      = lEta;
        fArrayTracksPhi[lNumberOfAcceptedTracksForLRC]      = lPhi;
        fArrayTracksCharge[lNumberOfAcceptedTracksForLRC]   = lCharge;
        fArrayTracksPID[lNumberOfAcceptedTracksForLRC]      = lMostProbablePIDPure;
        
        lNumberOfAcceptedTracksForLRC++;
        
    } //end of track loop
    
    //net charge vs pt
//    if ( lNacceptEtaInRegion > 20 )
//    {
//        fHistChargePlusVsPtTmp->DrawClone();
//        fHistChargeMinusVsPtTmp->SetLineColor( kBlue );
//        fHistChargeMinusVsPtTmp->DrawClone("same");
//    }


    //final actions for net-charge
    double integralPlus     = fHistChargePlusVsPtTmp->ComputeIntegral();
    double integralMinus    = fHistChargeMinusVsPtTmp->ComputeIntegral();
    double xAxisLowEdge     = fHistChargePlusVsPtTmp->GetBinLowEdge(1);
    double xAxisHighEdge    = fHistChargePlusVsPtTmp->GetBinLowEdge(fHistChargePlusVsPtTmp->GetNbinsX()) + fHistChargePlusVsPtTmp->GetBinWidth(fHistChargePlusVsPtTmp->GetNbinsX());
    double ptRange =  xAxisHighEdge - xAxisLowEdge;
//    cout << "integralPlus=" << integralPlus << ", integralMinus=" << integralMinus << ", ptRange=" << ptRange << endl;

    for ( int bin = 0; bin < fHistChargePlusVsPtTmp->GetNbinsX(); bin++)
    {
        double binCenter = fHistChargePlusVsPtTmp->GetBinCenter(bin+1);
        double binContentPlus = fHistChargePlusVsPtTmp->GetBinContent(bin+1);
        double binContentMinus = fHistChargeMinusVsPtTmp->GetBinContent(bin+1);



        double lNetCharge = binContentPlus-binContentMinus;
        fHist2DNetChargeVsPt->Fill( binCenter, lNetCharge );
        if ( ptRange>0 )
            fHist2DNetChargeVsPtCorrectedOnEventMean->Fill( binCenter, lNetCharge - (integralPlus-integralMinus)/ptRange );
        //        fHist2DNetChargeVsPtCorrectedOnEventMeanNormOnNch->Fill( );

    }

    //############ MC ESD QA
    if ( fRunKine && fAnalyseMCESD ) //run MC+ESD
    {
        Double_t lEta;
        Double_t lPt;
//        Double_t lY;
        Int_t lMCtracksAcceptedEtaInRegionAfterPtCuts = 0;
        for (Int_t iMCtracks = 0; iMCtracks < eventMC->GetNumberOfPrimaries(); iMCtracks++)
        {
            AliMCParticle* track = dynamic_cast<AliMCParticle *>(eventMC->GetTrack(iMCtracks));
            if (!track) {
                Printf("ERROR: Could not receive particle %d", iMCtracks);
                continue;
            }
            
            //exclude non stable particles
            if(!eventMC->IsPhysicalPrimary(iMCtracks)) continue;
            
            if ( track->Particle()->GetPDG()->Charge() == 0 )
                continue;
            
            lEta    = track->Eta();
            lPt     = track->Pt();
//            lY      = track->Y();
            
            if( fabs(lEta) > fEtaRegionForTests ) //look at eta region
                continue;
            
            if( lPt > fMaxPtLimit || lPt < fMinPtLimit )
                continue;
            
            
            lMCtracksAcceptedEtaInRegionAfterPtCuts++;
        }
        fHist2DMultiplicityMCESDInEtaRegion->Fill( lMCtracksAcceptedEtaInRegionAfterPtCuts, lNacceptEtaInRegionAfterPtCuts );
    }
    
    
    
    
    //profiling tracks in phi BY HAND if requested
    //    if ( fUsePhiShufflingByHand )
    //    {
    //        double lFictionEventPlanePhi = fRand->Uniform(0,2*TMath::Pi());
    //        fHistPhiArtificialEvPlane->Fill( lFictionEventPlanePhi );
    //        TF1 lFictionEventPlaneFunc( "fictionEventPlane", "0.75 + 0.25*TMath::Cos(x)", 0, 2*TMath::Pi() );
    //        for ( Int_t trackId = 0; trackId < lNumberOfAcceptedTracksForLRC; trackId++ )
    //        {
    //            //check phi wrt event plane
    //            Double_t lProfiledPhiWrtEvPlane = lFictionEventPlaneFunc.GetRandom();
    //            //FixAngleInTwoPi( lProfiledPhiWrtEvPlane );
    //            fHistPhiArtificialProfilingCheckWrtEvPlane->Fill( lProfiledPhiWrtEvPlane );
    
    //            //  double lOppositePhi = ( ( fRand->Uniform(0,1) > 0.5 ) ? -TMath::Pi() : 0 );
    //            //  double lRandomConeAngleEta = fRand->Gaus(0,0.5);
    //            //  double lRandomConeAnglePhi = fRand->Gaus(0,TMath::PiOver2());
    
    //            //double lFictionPhi = fRand->Uniform(0,2*TMath::Pi());
    //            //double lPhiSign = ( ( fRand->Uniform(0,1) > 0.5 ) ? TMath::Pi() : 0 );
    
    //            //add event plane phi angle
    //            Double_t lProfiledPhi = lProfiledPhiWrtEvPlane + lFictionEventPlanePhi;
    //            FixAngleInTwoPi( lProfiledPhi );
    //            fHistPhiArtificialProfilingCheck->Fill( lProfiledPhi );
    
    //            fArrayTracksPhi[trackId] = lProfiledPhi; //lPhi + lRandomConeAnglePhi; //lPhi;
    
    //            //        fArrayTracksPhi[lNumberOfAcceptedTracksForLRC]      = lPhi - lRandomConeAnglePhi; //lPhi;
    //            //        if ( fArrayTracksPhi[lNumberOfAcceptedTracksForLRC] > 2 * TMath::Pi() )
    //            //            fArrayTracksPhi[lNumberOfAcceptedTracksForLRC] -= 2 * TMath::Pi();
    //        }
    //    }
    
    
    // ######### fill some QA plots
    fHistAcceptedTracks->Fill(lNaccept);
    fHistMultiplicityInEtaRegion->Fill( lNacceptEtaInRegion );
    fHistMultiplicityInEtaRegionAfterPtCuts->Fill( lNacceptEtaInRegionAfterPtCuts );
    fHistAcceptedTracksAfterPtCuts->Fill(lNacceptAfterPtCuts);
    if ( fArtificialInefficiency >= 0 ) //i.e. have ineff analysis ON
        fHistNumberOfDroppedByHandTracks->Fill( lNacceptAfterPtCuts, lNumberOfDropByHandTracks );
    
    // ########### ProfilePhiByHand
    if ( fUsePhiShufflingByHand )
        ProfilePhiByHand( lNumberOfAcceptedTracksForLRC);
    
    // ########### LRC processors filling
    Bool_t lFlagAcceptEventForLRC = ( lNacceptEtaInRegion >= fMultCutInEtaRegion );
    if ( lFlagAcceptEventForLRC )
        FillLRCProcessors( lNumberOfAcceptedTracksForLRC, lCentralityPercentile );
    
    // ##### ZDC, VZERO, etc
    if ( fFlagWatchZDC && fIsIonsAnalysis )
    {
        AliESDZDC *zdcData = lESD->GetESDZDC();
        //fHistZDCenergy->Fill( zdcData->GetZDCEMEnergy(0) );
        fHistZDCparticipants->Fill( zdcData->GetZDCParticipants ());
        for ( int i = 0; i < 5; i++ )
        {
            //if multiplicity in eta windows above threshold - fill ZDC energy hist
            if ( countTracksEtaB >= fMultForZDCstudy[i] && countTracksEtaF >= fMultForZDCstudy[i] )
                fHistZDCenergy[i]->Fill( zdcData->GetZDCEMEnergy(0) );
        }
    }
    if ( fFlagWatchV0 )
    {
        const AliESDVZERO* vzrData = lESD->GetVZEROData(); //aod the same

        
        //V0 multiplicity
        float lV0mult[64];
        //float lV0Amult[64];
        //float lV0Cmult[64];
        double sumV0mult = 0;
        double sumV0Amult = 0;
        double sumV0Cmult = 0;
        
        for (int i = 0; i < 64; i++ )
        {
            lV0mult[i] = (float)(vzrData->GetMultiplicity(i));
            sumV0mult += lV0mult[i];
        }
        fHistV0multiplicity->Fill( sumV0mult );
        
        //            for (int i=0; i<32; i++)
        //            {
        //                lV0Cmult[i] = (float)(vzrData->GetMultiplicityV0C(i));
        //                sumV0Cmult += lV0Cmult[i];
        //                lV0Amult[i] = (float)(vzrData->GetMultiplicityV0A(i));
        //                sumV0Amult += lV0Amult[i];
        //            }
        sumV0Amult = vzrData->GetMTotV0A();
        sumV0Cmult = vzrData->GetMTotV0C();
        fHistV0Amultiplicity->Fill( sumV0Amult );
        fHistV0Cmultiplicity->Fill( sumV0Cmult );
        fHist2DV0ACmultiplicity->Fill( sumV0Amult, sumV0Cmult );
        
        fHist2DTracksAcceptedVsV0multiplicity->Fill( lNaccept, sumV0Amult );
        
        //cells
        //fHistV0cells->Fill(  );
        fHistV0Acells->Fill( vzrData->GetNbPMV0A() );
        fHistV0Ccells->Fill( vzrData->GetNbPMV0C() );
        fHist2DV0ACcells->Fill( vzrData->GetNbPMV0A(), vzrData->GetNbPMV0C() );
        
        //rings
        for ( int i = 0; i < 4; i++ )
        {
            int lMultRingV0A = vzrData->GetMRingV0A(i);
            int lMultRingV0C = vzrData->GetMRingV0C(i);
            fHistV0AmultiplicityRing[i]->Fill( lMultRingV0A );
            fHistV0CmultiplicityRing[i]->Fill( lMultRingV0C );
            fHist2DV0ACmultiplicityRing[i]->Fill( lMultRingV0A, lMultRingV0C );
            
            fHist2DTracksAcceptedVsV0AmultiplicityRing[i]->Fill( lNaccept, lMultRingV0A );
            fHist2DTracksAcceptedVsV0CmultiplicityRing[i]->Fill( lNaccept, lMultRingV0C );
        }
        
    }
    if ( fFlagWatchFMD )
    {
        // multiplicity
        //const AliESDFMD *lFMDdata = lESD->GetFMDData();
        
        
    }
    
    //        //look at V0s
    //        AliESDv0 *lV0 = 0x0;
    //        //cout << ">>>> showing v0's:" << endl;
    //        for ( int iV0 = 0; iV0 < lESD->GetNumberOfV0s(); iV0++ )
    //        {
    //            lV0 = lESD->GetV0( iV0 );
    
    //            int lV0pdgCode = lV0->GetPdgCode();
    //            double lMassV0 = lV0->GetEffMass();//ChangeMassHypothesis( lV0pdgCode );  //GetEffMass();
    //            //if ( abs(lV0pdgCode) > 10 )//== 3122 )
    //            {
    //                //cout << ">>>>>>>>>>>>>>>>>>>>>>>>> " << lV0pdgCode << ": mass = " << lMassV0 << endl;
    //                fHistV0spectra->Fill(lV0pdgCode);//lMassV0);
    //            }
    
    //            if ( lMassV0 > 1. )// abs(lPdgCode) == 3122 || abs(lPdgCode) == 421 ) //abs(lPdgCode) % 1000 == 122 )
    //            {
    //                //cout << ">>>>>>>>>>>>>>>>>>>>>>>>> " << lV0pdgCode << ": mass = " << lMassV0 << endl;
    //            }
    //        }
    
    
    
    
    
    
    //Debuging output of track filter
    if( fShowPerTrackStats )
        Printf("NofTracks= %d , accepted %d , LowPt %d, HighPt %d, LowCharge %d,  ESD cut fail %d , Decay Filer %d", event->GetNumberOfTracks(), lNaccept, lPtUnder, lPtOver ,lNoCharge , lCutfail  , lDecay);
    
    // Post output data
    PostData(1, fOutList);
    
    for( Int_t i = 0; i < fLRCproc.GetEntries(); i++)
    {
        PostData( Proc(i)->GetOutputSlotNumber(),Proc(i)->CreateOutput() );
    }
    
}


void AliAnalysisTaskLRC::SetParticleTypeToProcessors( int windowId, TString strPid ) //char* strPid )//AliLRCBase::LRCparticleType whichParticleToFill )
{
    //dummy actions
    windowId+=1;
    if(0) printf("%s", strPid.Data());

    //Set pid types for LRC processors (windowId = 0 - fwd window, =1 - backward window)
    //    Int_t lLrcNum=fLRCproc.GetEntries();
    
    //    for(Int_t i=0; i < lLrcNum; i++)
    //    {
    //        AliLRCBase *lrcProc = dynamic_cast<AliLRCBase*> (fLRCproc.At(i));
    //        if(lrcProc)
    //        {
    //            //tmp lines!
    //            windowId++;
    //            if (0)
    //                printf("%s",strPid);
    //            //commented for this commit
    //            //            if ( windowId == 0 )
    //            //                lrcProc->SetParticleType( "fwd", strPid );
    //            //            else if ( windowId == 1 )
    //            //                lrcProc->SetParticleType( "bkwd", strPid );
    
    //            //fTmpCounter++;
    //        }
    //        else continue;
    //        //fTmpCounter++;
    //    }
    
}

void AliAnalysisTaskLRC::SetParticleTypeForTask( TString strF, TString strB ) //char* strF, char* strB )
{
    // Set PID sensitivity to 'true' and write pids for windows
    fPIDsensingFlag = kTRUE;
    fStrPIDforFwd = strF;
    fStrPIDforBwd = strB;

    //    sprintf ( fStrPIDforFwd, "%s", strF );
    //    sprintf ( fStrPIDforBwd, "%s", strB );
    //cout << "we just set fwd win to " << fStrPIDforFwd << " and bwd win to " << fStrPIDforBwd << endl;
}



//________________________________________________________________________
void AliAnalysisTaskLRC::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    //fOutList = dynamic_cast<TList*> (GetOutputData(0));
    
    //lESDTrackCuts->DrawHistograms();
    
    
    fAnalysisTimer->Stop();
    Double_t rtime = fAnalysisTimer->RealTime();
    Double_t ctime = fAnalysisTimer->CpuTime();
    
    printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
    
}


void AliAnalysisTaskLRC::AddLRCProcess(AliLRCBase *newProc)
{
    // Add new AliLRCBase (Main LRC processor per ETA window) to the processing list
    // Used to add new ETA window to AnalysisTask
    if(!newProc)
    {
        Printf("ERROR:No AliLRCBase object -  NULL pointer!");
        return;
    }
    
    fLRCproc.Add(newProc);
    newProc->SetOutputSlotNumber(fLRCproc.GetEntries() + 1 );//( fSetIncludeEventTreeInOutput ? 2 : 1 ) );
    DefineOutput(newProc->GetOutputSlotNumber(),TList::Class());
    return ;
}


Double_t AliAnalysisTaskLRC::GetEventPlane(AliVEvent *event)
{
    // Get the event plane
    
    if ( !event )
        return -1;
    //TString gAnalysisLevel = fBalance->GetAnalysisLevel();
    
    Float_t lVZEROEventPlane    = -10.;
    Float_t lReactionPlane      = -10.;
    Double_t qxTot = 0.0, qyTot = 0.0;
    
    //MC: from reaction plane
    if( fRunKine )//gAnalysisLevel == "MC"){
    {
        AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(dynamic_cast<AliMCEvent*>(event)->GenEventHeader());
        if (headerH)
        {
            lReactionPlane = headerH->ReactionPlaneAngle();
            //lReactionPlane *= TMath::RadToDeg();
        }
    }//MC
    
    // AOD,ESD,ESDMC: from VZERO Event Plane
    else
    {
        AliEventplane *ep = event->GetEventplane();
        if( ep )
        {
            lVZEROEventPlane = ep->CalculateVZEROEventPlane( event, 10, 2, qxTot, qyTot);
            if( lVZEROEventPlane < 0. )
                lVZEROEventPlane += TMath::Pi();
            //lReactionPlane = lVZEROEventPlane*TMath::RadToDeg();
            lReactionPlane = lVZEROEventPlane;
            //            cout << "lReactionPlane: " << lReactionPlane << endl;
            //            int aa;
            //            cin >> aa;
        }
    }//AOD,ESD,ESDMC
    
    return lReactionPlane;
}


void AliAnalysisTaskLRC::AddTrackCutForBits(AliESDtrackCuts*  const cuts, TString cutsName )
{
    //TString *strPtrCutsName = new TString(cutsName.Data());
//    cout << "ADDING CUTS " << cutsName << endl;
    //    fArrTrackCuts[fNumberOfCutsToRemember] = cuts;
    int lNumberOfCutsToRemember = fArrTrackCuts.GetEntries();
    
    fArrCutsNames[lNumberOfCutsToRemember] = cutsName;
    fArrTrackCuts.Add( cuts );// = cuts;
    //    fArrCutsNames.Add( strPtrCutsName );
    //    fNumberOfCutsToRemember++;
}

void AliAnalysisTaskLRC::UseToyEvents()
{
    //generate all events here and fill LRC processors
    for ( Int_t toyEventId = 0; toyEventId < fNumberOfToyEvents; toyEventId++ )
    {
        if ( toyEventId % 10000 == 0 )
            printf("Processing %d event...\n", toyEventId );

        //        cout << "toyEventId=" << toyEventId << endl;
        int nToyTracks =  0;
        const double kExpParam = 0.42;
        const double kEtaGausSigma = 0.3;
        const double kPhiGausSigma = TMath::PiOver4();

        //try to do "jets" with some particles
        const int nJets = 2;

        //jets eta angles
        double lJets_Eta[nJets];
        lJets_Eta[0] = fRand->Uniform(0, 0.8);
        lJets_Eta[1] = -lJets_Eta[0];

        //jets phi angles
        double lJets_Phi[nJets];
        lJets_Phi[0] = fRand->Uniform(0, 2*TMath::Pi());
        lJets_Phi[1] = lJets_Phi[0] + TMath::Pi();
        FixAngleInTwoPi( lJets_Phi[1] );


        //jets loop
        for ( Int_t jetId = 0; jetId < nJets; jetId++ )
        {
            int nPartInJet = fRand->Gaus(4,1);
            //jet tracks params
            for ( Int_t toyTrackId = 0; toyTrackId < nPartInJet; toyTrackId++ )
            {
                fArrayTracksPt[nToyTracks] = fRand->Exp( kExpParam );
                fArrayTracksEta[nToyTracks] = fRand->Gaus( lJets_Eta[jetId], kEtaGausSigma );
                fArrayTracksPhi[nToyTracks] = fRand->Gaus( lJets_Phi[jetId], kPhiGausSigma );
                FixAngleInTwoPi( fArrayTracksPhi[nToyTracks] );

                nToyTracks++;
            }
        }

        FillLRCProcessors( nToyTracks, 0 );

        //        for ( Int_t toyTrackId = 0; toyTrackId < nToyTracks; toyTrackId++ )
        //        {
        //            // ########### ProfilePhiByHand
        //            ProfilePhiByHand( nToyTracks);
        //        }

    }
}


void AliAnalysisTaskLRC::FillLRCProcessors( int numberOfAcceptedTracksForLRC, Double_t eventCentrality )
{
    //pass signal to LRC-based analysis to start event
    Int_t lLrcNum = fLRCproc.GetEntries(); // Number of processors attached
    
    //prepare phi rotation step
    Double_t phiStep = 2 * TMath::Pi();
    if ( fNumberOfPhiSectors > 1 )
        phiStep /= fNumberOfPhiSectors;
    Double_t lPhiRotatedExtra = 0; //additional phi rotation of all tracks in case
    
    
    //start LRC filling
    for ( Int_t sectorId = 0; sectorId < fNumberOfPhiSectors; sectorId++ )
    {
        if ( lLrcNum > kMaxParticlesNumber )
        {
            AliDebug(AliLog::kError, "lLrcNum is too large... Mistake???");
            break;
        }
        
        //pass signal to LRC-based analysis to start event
        for( Int_t lrcProcessorId = 0; lrcProcessorId < lLrcNum; lrcProcessorId++ )
        {
            AliLRCBase *lrcBase = dynamic_cast<AliLRCBase*> (fLRCproc.At(lrcProcessorId));
            if (!lrcBase)
            {
                AliDebug(AliLog::kError, "can't cast lrcBase to AliLRCBase!" );
                break;
            }
            lrcBase->StartEvent();
            //pass the centrality
            lrcBase->SetEventCentrality( eventCentrality );
            //        }
            
            //pass track data to LRC-based analysis
            for ( Int_t trackId = 0; trackId < numberOfAcceptedTracksForLRC; trackId++ )
            {
                if ( fEtInsteadOfPt && (fArrayTracksPt[trackId]  == 0 ) ) //in Et-mode we didn't get Et! exit
                {
                    AliDebug(AliLog::kError, "pT=0 Mistake???" );
                    break;
                }               //if ( fEtInsteadOfPt && lMostProbablePIDPure == -1 ) //don't have pure PID
                //	continue;
                
                //rotate track phi
                Double_t lPhi = fArrayTracksPhi[trackId] + lPhiRotatedExtra;
                FixAngleInTwoPi( lPhi );
                
                fHistPhiLRCrotationsCheck->Fill( lPhi );
                
                //            for(Int_t i = 0; i < lLrcNum; i++ )
                //            {
                lrcBase->AddTrackPtEta(
                            fArrayTracksPt[trackId]
                            , fArrayTracksEta[trackId]
                            , lPhi
                            , fArrayTracksCharge[trackId]
                            , fArrayTracksPID[trackId]
                            );
                //            }
                
            }
            
            //take event only if at least 1 track in this event fulfill the requirements! //21.11.11
            Bool_t lDontTakeEventDecision =  kFALSE; //lNaccept > 0 ? kFALSE : kTRUE;
            //pass signal to LRC-based analysis to finish the event
            //        for( Int_t i = 0; i < lLrcNum; i++ )
            //        {
            lrcBase->FinishEvent( lDontTakeEventDecision );
        }
        
        lPhiRotatedExtra += phiStep; //increase phi step to rotate tracks
    }
}

void AliAnalysisTaskLRC::ProfilePhiByHand( int numberOfAcceptedTracksForLRC )
{
    //profiling tracks in phi BY HAND if requested
    
    double lFictionEventPlanePhi = fRand->Uniform(0,2*TMath::Pi());
    fHistPhiArtificialEvPlane->Fill( lFictionEventPlanePhi );
    TF1 lFictionEventPlaneFunc( "fictionEventPlane", "0.75 + 0.25*TMath::Cos(x)", 0, 2*TMath::Pi() );
    for ( Int_t trackId = 0; trackId < numberOfAcceptedTracksForLRC; trackId++ )
    {
        //check phi wrt event plane
        Double_t lProfiledPhiWrtEvPlane = lFictionEventPlaneFunc.GetRandom();
        //FixAngleInTwoPi( lProfiledPhiWrtEvPlane );
        fHistPhiArtificialProfilingCheckWrtEvPlane->Fill( lProfiledPhiWrtEvPlane );
        
        //  double lOppositePhi = ( ( fRand->Uniform(0,1) > 0.5 ) ? -TMath::Pi() : 0 );
        //  double lRandomConeAngleEta = fRand->Gaus(0,0.5);
        //  double lRandomConeAnglePhi = fRand->Gaus(0,TMath::PiOver2());
        
        //double lFictionPhi = fRand->Uniform(0,2*TMath::Pi());
        //double lPhiSign = ( ( fRand->Uniform(0,1) > 0.5 ) ? TMath::Pi() : 0 );
        
        //add event plane phi angle
        Double_t lProfiledPhi = lProfiledPhiWrtEvPlane + lFictionEventPlanePhi;
        FixAngleInTwoPi( lProfiledPhi );
        fHistPhiArtificialProfilingCheck->Fill( lProfiledPhi );
        
        fArrayTracksPhi[trackId] = lProfiledPhi; //lPhi + lRandomConeAnglePhi; //lPhi;
        
        //        fArrayTracksPhi[lNumberOfAcceptedTracksForLRC]      = lPhi - lRandomConeAnglePhi; //lPhi;
        //        if ( fArrayTracksPhi[lNumberOfAcceptedTracksForLRC] > 2 * TMath::Pi() )
        //            fArrayTracksPhi[lNumberOfAcceptedTracksForLRC] -= 2 * TMath::Pi();
    }
}
