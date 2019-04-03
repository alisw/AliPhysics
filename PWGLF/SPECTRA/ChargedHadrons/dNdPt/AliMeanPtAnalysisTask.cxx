#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliPhysicsSelection.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"

#include "AliMCEvent.h"

#include "AliMeanPtAnalysisTask.h"

/// \cond CLASSIMP
ClassImp(AliMeanPtAnalysisTask);
/// \endcond


//________________________________________________________________________
AliMeanPtAnalysisTask::AliMeanPtAnalysisTask(const char* name) : AliAnalysisTaskSE(name),
  //General member variables
  fPRECISION(1e-6),
  fOutputList(0),
  fEvent(0),
  fMCEvent(0),
  fESDtrackCuts(0),
  fUtils(0),
  //Toggles
  fIsESD(kTRUE),
  fIsMC(kFALSE),
  fIs2013pA(kFALSE),
  fIs2015data(kFALSE),
  fTPCRefit(kFALSE),
  fITSRefit(kFALSE),
  fAcceptKinks(kTRUE),
  fRequiresClusterITS(kTRUE),
  fDCAToVertex2D(kFALSE),
  fSigmaToVertex(kFALSE),
  fUseGeomCut(kFALSE),
  fIncludeCrosscheckHistos(kFALSE),
  // Cut Parameters
  fTriggerMask(AliVEvent::kINT7),
  fMinEta(-10),
  fMaxEta(10),
  fMinPt(0),
  fMaxPt(999),
  fSigmaMeanXYZv(),
  fMeanXYZv(),
  fZvtx(10),
  fMinNCrossedRowsTPC(0),
  fMinRatioCrossedRowsOverFindableClustersTPC(0),
  fMaxFractionSharedClustersTPC(0),
  fMaxChi2PerTPCCluster(0),
  fMaxChi2PerITSCluster(0),
  fMaxDCAzITSTPC(0),
  fDCAToVertexXYPtDep("0"),
  fDCAToVertexXY(0),
  fMaxChi2TPCConstrained(0),
  fMinActiveLength(0),
  fDeadZoneWidth(2),
  fCutGeoNcrNclLenght(130),
  fCutGeoNcrNclGeom1Pt(1.5),
  fCutGeoNcrNclFractionNcr(0.85),
  fCutGeoNcrNclFractionNcl(0.7),
  //Arrays for Binning
  fBinsMult(0),
  fBinsCent(0),
  fBinsPt(0),
  fBinsEta(0),
  fBinsZv(0),
  fBinsPtReso(0),
  fBins1Pt(0),
  fBinsSigma1Pt(0),
  //Event-Histograms
  fEventCount(0),
  fHistMCTrackParticle(0),
  fHistEvent(0),
  fHistMCResponseMat(0),
  fHistMCResponseMatTracks(0),
  //Track-Histograms
  fHistTrack(0),
  fHistRelPtResoFromCov(0),
  fHistMCRecTrack(0),
  fHistMCGenPrimTrack(0),
  fHistMCRecPrimTrack(0),
  fHistMCRecSecTrack(0),
  fHistMCMultPtGenerated(0),
  fHistMCTrackMultGen(0),
  fHistMCPtRes(0),
  fHistMCRelPtReso(0),
  fHistMCEtaRes(0),
  fHistMCMultRes(0),
  fHistMCParticle(0)
{
  // Set default binning
  Double_t binsMultDefault[2] = {0., 10000.};
  Double_t binsCentDefault[2] = {0., 100.};
  Double_t binsPtDefault[49] = {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0};
  Double_t binsEtaDefault[19] = {-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
  //  Double_t binsPtDefault[69] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};
  //  Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};

  // binning for relative pT resolution
  const Int_t nBinsPtReso = 300;
  Double_t binsPtReso[nBinsPtReso+1];
  SetFixedBinEdges(binsPtReso, 0., 0.3, nBinsPtReso);
  SetBinsPtReso(nBinsPtReso, binsPtReso);

  // binning for 1/pt
  const Int_t nBins1Pt = 200;
  Double_t bins1Pt[nBins1Pt+1];
  SetFixedBinEdges(bins1Pt, 0., 10., nBins1Pt);
  SetBins1Pt(nBins1Pt, bins1Pt);

  // binning for sigma 1/pt
  const Int_t nBinsSigma1Pt = 200;
  Double_t binsSigma1Pt[nBinsSigma1Pt+1];
  SetFixedBinEdges(binsSigma1Pt, 0., 0.1, nBinsSigma1Pt);
  SetBinsSigma1Pt(nBinsSigma1Pt, binsSigma1Pt);


  SetBinsMult(1,binsMultDefault);
  SetBinsCent(1,binsCentDefault);
  SetBinsPt(48,binsPtDefault);
  SetBinsEta(18,binsEtaDefault);
  SetBinsZv(12,binsZvDefault);
  SetMeanXYZv(0.0,0.0,0.0);
  SetSigmaMeanXYZv(1.0,1.0,10.0);
  SetZvtx(10.);
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliMeanPtAnalysisTask::UserCreateOutputObjects(){
  // Create histograms here (function is called once)
  OpenFile(1,"recreate");
  fOutputList = new TList();
  fOutputList -> SetOwner();


  // Control histogram to check the effect of event cuts
  fEventCount = new TH1F("fEventCount","Number of events after cuts",4,0.5,4.5);
  fEventCount->GetYaxis()->SetTitle("#it{N}_{events}");
  fEventCount->GetXaxis()->SetBinLabel(1, "all");
  fEventCount->GetXaxis()->SetBinLabel(2, "triggered");
  fEventCount->GetXaxis()->SetBinLabel(3, "Vertex ok, no Pileup");
  fEventCount->GetXaxis()->SetBinLabel(4, "without #it{N}_{ch} overflow");

  /// Event histogram Nacc:cent
  Int_t nBinsEvent[2]={fBinsMult->GetSize()-1, fBinsCent->GetSize()-1};
  Double_t minEvent[2]={fBinsMult->GetAt(0), fBinsCent->GetAt(0)};
  Double_t maxEvent[2]={fBinsMult->GetAt(fBinsMult->GetSize()-1), fBinsCent->GetAt(fBinsCent->GetSize()-1)};

  fHistEvent = new THnF("fHistEvent", "Histogram for Events mult vs cent",2,nBinsEvent,minEvent,maxEvent);
  fHistEvent -> SetBinEdges(0,fBinsMult->GetArray());
  fHistEvent -> SetBinEdges(1,fBinsCent->GetArray());
  fHistEvent->GetAxis(0)->SetTitle("#it{N}_{acc}");
  fHistEvent->GetAxis(1)->SetTitle("Centrality (%)");
  fHistEvent -> Sumw2();

  /// Track histogram pt:eta:zV:multcent
  Int_t nBinsTrack[4]={fBinsPt->GetSize()-1,fBinsEta->GetSize()-1,fBinsMult->GetSize()-1,fBinsCent->GetSize()-1};
  Double_t minTrack[4]={fBinsPt->GetAt(0),fBinsEta->GetAt(0),fBinsMult->GetAt(0),fBinsCent->GetAt(0)};
  Double_t maxTrack[4]={fBinsPt->GetAt(fBinsPt->GetSize()-1),fBinsEta->GetAt(fBinsEta->GetSize()-1),fBinsMult->GetAt(fBinsMult->GetSize()-1),fBinsCent->GetAt(fBinsCent->GetSize()-1)};

  fHistTrack = new THnF("fHistTrack", "Histogram for Tracks",4,nBinsTrack,minTrack,maxTrack);
  fHistTrack -> SetBinEdges(0,fBinsPt->GetArray());
  fHistTrack -> SetBinEdges(1,fBinsEta->GetArray());
  fHistTrack -> SetBinEdges(2,fBinsMult->GetArray());
  fHistTrack -> SetBinEdges(3,fBinsCent->GetArray());
  fHistTrack->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistTrack->GetAxis(1)->SetTitle("#eta");
  fHistTrack->GetAxis(2)->SetTitle("#it{N}_{acc}");
  fHistTrack->GetAxis(3)->SetTitle("Centrality (%)");
  fHistTrack -> Sumw2();


  /// relative pT resolution from covariance matrix (global tracks) as a function of pt and centrality
  Int_t nBinsRelPtReso[3]  = {fBinsPtReso->GetSize()-1,fBinsPt->GetSize()-1, fBinsCent->GetSize()-1};
  Double_t minRelPtReso[3] = {fBinsPtReso->GetAt(0),fBinsPt->GetAt(0), fBinsCent->GetAt(0)};
  Double_t maxRelPtReso[3] = {fBinsPtReso->GetAt(fBinsPtReso->GetSize()-1), fBinsPt->GetAt(fBinsPt->GetSize()-1), fBinsCent->GetAt(fBinsCent->GetSize()-1)};

  fHistRelPtResoFromCov = new THnF("fHistRelPtResoFromCov", "Relative pT resolution from covariance matrix", 3, nBinsRelPtReso, minRelPtReso, maxRelPtReso);
  fHistRelPtResoFromCov -> SetBinEdges(0,fBinsPtReso->GetArray());
  fHistRelPtResoFromCov -> SetBinEdges(1,fBinsPt->GetArray());
  fHistRelPtResoFromCov -> SetBinEdges(2,fBinsCent->GetArray());
  fHistRelPtResoFromCov ->GetAxis(0)->SetTitle("#sigma(#it{p}_{T}) / #it{p}_{T}");
  fHistRelPtResoFromCov ->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistRelPtResoFromCov ->GetAxis(2)->SetTitle("Centrality (%)");
  fHistRelPtResoFromCov -> Sumw2();


  if(fIsMC){

    // Control histogram showing reconstructed Tracks vs. reconstructed Particles
    fHistMCTrackParticle = new TH1F("fHistMCTrackParticle","Reconstructed MC Tracks and Particles",2,0.5,2.5);
    fHistMCTrackParticle->GetYaxis()->SetTitle("#it{N}");
    fHistMCTrackParticle->GetXaxis()->SetBinLabel(1, "#it{N}_{acc}");
    fHistMCTrackParticle->GetXaxis()->SetBinLabel(2, "#it{N}_{rec}");

    // Response Matrix Histogram
    Int_t nBinsRespMat[2]={fBinsMult->GetSize()-1,fBinsMult->GetSize()-1};
    Double_t minMultRespMat[2]={fBinsMult->GetAt(0),fBinsMult->GetAt(0)};
    Double_t maxMultRespMat[2]={fBinsMult->GetAt(fBinsMult->GetSize()-1),fBinsMult->GetAt(fBinsMult->GetSize()-1)};

    fHistMCResponseMat = new THnSparseF("fHistMCResponseMat","Response Matrix",2,nBinsRespMat, minMultRespMat, maxMultRespMat);
    fHistMCResponseMat->SetBinEdges(0,fBinsMult->GetArray());
    fHistMCResponseMat->SetBinEdges(1,fBinsMult->GetArray());
    fHistMCResponseMat->GetAxis(0)->SetTitle("#it{N}_{acc}");
    fHistMCResponseMat->GetAxis(1)->SetTitle("#it{N}_{ch}");
    fHistMCResponseMat->Sumw2();

    /// Track histogram with multiplicity correlation pt:multacc:multgen
    Int_t nBinsMultTrack[3]={fBinsMult->GetSize()-1,fBinsMult->GetSize()-1, fBinsPt->GetSize()-1};
    Double_t minMultTrack[3]={fBinsMult->GetAt(0),fBinsMult->GetAt(0), fBinsPt->GetAt(0)};
    Double_t maxMultTrack[3]={fBinsMult->GetAt(fBinsMult->GetSize()-1),fBinsMult->GetAt(fBinsMult->GetSize()-1), fBinsPt->GetAt(fBinsPt->GetSize()-1)};

    fHistMCResponseMatTracks = new THnSparseF("fHistMCResponseMatTracks","Response Matrix for primary particles",3,nBinsMultTrack, minMultTrack, maxMultTrack);
    fHistMCResponseMatTracks->SetBinEdges(0,fBinsMult->GetArray());
    fHistMCResponseMatTracks->SetBinEdges(1,fBinsMult->GetArray());
    fHistMCResponseMatTracks->SetBinEdges(2,fBinsPt->GetArray());
    fHistMCResponseMatTracks->GetAxis(0)->SetTitle("#it{N}_{acc}");
    fHistMCResponseMatTracks->GetAxis(1)->SetTitle("#it{N}_{ch}");
    fHistMCResponseMatTracks->GetAxis(2)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistMCResponseMatTracks->Sumw2();


    /// pT Resolution Histogram
    Int_t nBinsMCPtPt[2]={fBinsPt->GetSize()-1, fBinsPt->GetSize()-1};
    Double_t minBinsMCPtPt[2]={fBinsPt->GetAt(0),fBinsPt->GetAt(0)};
    Double_t maxBinsMCPtPt[2]={fBinsPt->GetAt(fBinsPt->GetSize()-1), fBinsPt->GetAt(fBinsPt->GetSize()-1)};

    fHistMCPtRes = new THnF("fHistMCPtRes","#it{p}_{T} resolution", 2, nBinsMCPtPt,minBinsMCPtPt,maxBinsMCPtPt);
    fHistMCPtRes -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCPtRes -> SetBinEdges(1,fBinsPt->GetArray());
    fHistMCPtRes->GetAxis(0)->SetTitle("#it{p}_{T}^{track} (GeV/#it{c})");
    fHistMCPtRes->GetAxis(1)->SetTitle("#it{p}_{T}^{particle} (GeV/#it{c})");


    /// relative pT resolution from MC as a function of pt generated, reconstructed and centrality
    Int_t nBinsMCRelPtReso[4]={fBinsPtReso->GetSize()-1,fBinsPt->GetSize()-1, fBinsPt->GetSize()-1, fBinsCent->GetSize()-1};
    Double_t minMCRelPtReso[4]={fBinsPtReso->GetAt(0),fBinsPt->GetAt(0), fBinsPt->GetAt(0), fBinsCent->GetAt(0)};
    Double_t maxMCRelPtReso[4]={fBinsPtReso->GetAt(fBinsPtReso->GetSize()-1),fBinsPt->GetAt(fBinsPt->GetSize()-1), fBinsPt->GetAt(fBinsPt->GetSize()-1), fBinsCent->GetAt(fBinsCent->GetSize()-1)};

    fHistMCRelPtReso = new THnF("fHistMCRelPtReso", "Relative pT resolution MC", 4, nBinsMCRelPtReso, minMCRelPtReso, maxMCRelPtReso);
    fHistMCRelPtReso -> SetBinEdges(0,fBinsPtReso->GetArray());
    fHistMCRelPtReso -> SetBinEdges(1,fBinsPt->GetArray());
    fHistMCRelPtReso -> SetBinEdges(2,fBinsPt->GetArray());
    fHistMCRelPtReso -> SetBinEdges(3,fBinsCent->GetArray());
    fHistMCRelPtReso ->GetAxis(0)->SetTitle("#Delta(#it{p}_{T}) / #it{p}_{T}");
    fHistMCRelPtReso ->GetAxis(1)->SetTitle("#it{p}^{MC, gen}_{T} (GeV/#it{c})");
    fHistMCRelPtReso ->GetAxis(2)->SetTitle("#it{p}^{MC, rec}_{T} (GeV/#it{c})");
    fHistMCRelPtReso ->GetAxis(3)->SetTitle("Centrality (%)");
    fHistMCRelPtReso -> Sumw2();



    /// Eta Resolution Histogram
    Int_t nBinsMCEtaEta[2]={fBinsEta->GetSize()-1, fBinsEta->GetSize()-1};
    Double_t minBinsMCEtaEta[2]={fBinsEta->GetAt(0), fBinsEta->GetAt(0)};
    Double_t maxBinsMCEtaEta[2]={fBinsEta->GetAt(fBinsEta->GetSize()-1), fBinsEta->GetAt(fBinsEta->GetSize()-1)};

    fHistMCEtaRes = new THnF("fHistMCEtaRes","#eta resolution", 2, nBinsMCEtaEta,minBinsMCEtaEta,maxBinsMCEtaEta);
    fHistMCEtaRes -> SetBinEdges(0,fBinsEta->GetArray());
    fHistMCEtaRes -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCEtaRes->GetAxis(0)->SetTitle("#eta^{track}");
    fHistMCEtaRes->GetAxis(1)->SetTitle("#eta^{particle}");

    // Binning for the correction histograms
    Int_t nBinsMCTrack[3] = {fBinsPt->GetSize()-1,fBinsEta->GetSize()-1, fBinsCent->GetSize()-1};
    Double_t minMCTrack[3] = {fBinsPt->GetAt(0),fBinsEta->GetAt(0), fBinsCent->GetAt(0)};
    Double_t maxMCTrack[3] = {fBinsPt->GetAt(fBinsPt->GetSize()-1),fBinsEta->GetAt(fBinsEta->GetSize()-1), fBinsCent->GetAt(fBinsCent->GetSize()-1)};

    fHistMCRecTrack = new THnF("fHistMCRecTrack", "Histogram for reconstructed MC Tracks",3,nBinsMCTrack,minMCTrack,maxMCTrack);
    fHistMCRecTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecTrack -> SetBinEdges(2,fBinsCent->GetArray());
    fHistMCRecTrack->GetAxis(0)->SetTitle("#it{p}_{T}^{MC} (GeV/#it{c})");
    fHistMCRecTrack->GetAxis(1)->SetTitle("#eta^{MC}");
    fHistMCRecTrack->GetAxis(2)->SetTitle("Centrality (%)");
    fHistMCRecTrack -> Sumw2();

    fHistMCGenPrimTrack = new THnF("fHistMCGenPrimTrack", "Histogram for generated MC Tracks",3,nBinsMCTrack,minMCTrack,maxMCTrack);
    fHistMCGenPrimTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCGenPrimTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCGenPrimTrack -> SetBinEdges(2,fBinsCent->GetArray());
    fHistMCGenPrimTrack->GetAxis(0)->SetTitle("#it{p}_{T}^{MC} (GeV/#it{c})");
    fHistMCGenPrimTrack->GetAxis(1)->SetTitle("#eta^{MC}");
    fHistMCGenPrimTrack->GetAxis(2)->SetTitle("Centrality (%)");
    fHistMCGenPrimTrack -> Sumw2();

    fHistMCRecPrimTrack = new THnF("fHistMCRecPrimTrack", "Histogram for reconstructed primary MC Tracks",3,nBinsMCTrack,minMCTrack,maxMCTrack);
    fHistMCRecPrimTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecPrimTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecPrimTrack -> SetBinEdges(2,fBinsCent->GetArray());
    fHistMCRecPrimTrack->GetAxis(0)->SetTitle("#it{p}_{T}^{MC} (GeV/#it{c})");
    fHistMCRecPrimTrack->GetAxis(1)->SetTitle("#eta^{MC}");
    fHistMCRecPrimTrack->GetAxis(2)->SetTitle("Centrality (%)");
    fHistMCRecPrimTrack -> Sumw2();

    fHistMCRecSecTrack = new THnF("fHistMCRecSecTrack", "Histogram for reconstructed secondary MC Tracks",3,nBinsMCTrack,minMCTrack,maxMCTrack);
    fHistMCRecSecTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecSecTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecSecTrack -> SetBinEdges(2,fBinsCent->GetArray());
    fHistMCRecSecTrack->GetAxis(0)->SetTitle("#it{p}_{T}^{MC} (GeV/#it{c})");
    fHistMCRecSecTrack->GetAxis(1)->SetTitle("#eta^{MC}");
    fHistMCRecSecTrack->GetAxis(2)->SetTitle("Centrality (%)");
    fHistMCRecSecTrack -> Sumw2();


    // MC truth pt vs. Nch (reference for closure test)
    Int_t nBinsMultPt[3]={fBinsMult->GetSize()-1,fBinsPt->GetSize()-1, fBinsCent->GetSize()-1};
    Double_t minMultPt[3]={fBinsMult->GetAt(0),fBinsPt->GetAt(0), fBinsCent->GetAt(0)};
    Double_t maxMultPt[3]={fBinsMult->GetAt(fBinsMult->GetSize()-1),fBinsPt->GetAt(fBinsPt->GetSize()-1), fBinsCent->GetAt(fBinsCent->GetSize()-1)};

    fHistMCMultPtGenerated = new THnF("fHistMCMultPtGenerated","Histogram for generator comparison and closure test",3,nBinsMultPt, minMultPt, maxMultPt);
    fHistMCMultPtGenerated->SetBinEdges(0,fBinsMult->GetArray());
    fHistMCMultPtGenerated->SetBinEdges(1,fBinsPt->GetArray());
    fHistMCMultPtGenerated->SetBinEdges(2,fBinsCent->GetArray());
    fHistMCMultPtGenerated->GetAxis(0)->SetTitle("#it{N}_{ch}");
    fHistMCMultPtGenerated->GetAxis(1)->SetTitle("#it{p}_{T}^{MC}");
    fHistMCMultPtGenerated->GetAxis(2)->SetTitle("Centrality (%)");
    fHistMCMultPtGenerated->Sumw2();


    if (fIncludeCrosscheckHistos){

        // Histogram to crosscheck pt resolution effect on resulting <pt> vs Nch
        // (use this instead of fHistTrack as input for unfolding)
        // could also be useful as estimate for systematic effect of pt resolution
        fHistMCParticle = new THnF("fHistMCParticle", "Histogram for reconstructed MC particles",4,nBinsTrack,minTrack,maxTrack);
        fHistMCParticle -> SetBinEdges(0,fBinsPt->GetArray());
        fHistMCParticle -> SetBinEdges(1,fBinsEta->GetArray());
        fHistMCParticle -> SetBinEdges(2,fBinsMult->GetArray());
        fHistMCParticle -> SetBinEdges(3,fBinsCent->GetArray());
        fHistMCParticle->GetAxis(0)->SetTitle("#it{p}_{T}^{MC} (GeV/#it{c})");
        fHistMCParticle->GetAxis(1)->SetTitle("#eta^{MC}");
        fHistMCParticle->GetAxis(2)->SetTitle("#it{N}_{acc}");
        fHistMCParticle->GetAxis(3)->SetTitle("Centrality (%)");
        fHistMCParticle -> Sumw2();


        // Histogram to illustrate <pt>(Nacc) for fixed Nch and  <pt>(Nch) for fixed Nacc
        fHistMCTrackMultGen = new THnSparseF("fHistMCTrackMultGen", "True Tracks as function of measured and true Mult", 3, nBinsMultTrack, minMultTrack, maxMultTrack);
        fHistMCTrackMultGen -> SetBinEdges(0,fBinsMult->GetArray());
        fHistMCTrackMultGen -> SetBinEdges(1,fBinsMult->GetArray());
        fHistMCTrackMultGen -> SetBinEdges(2,fBinsPt->GetArray());
        fHistMCTrackMultGen->GetAxis(0)->SetTitle("#it{N}_{acc}");
        fHistMCTrackMultGen->GetAxis(1)->SetTitle("#it{N}_{ch}");
        fHistMCTrackMultGen->GetAxis(2)->SetTitle("#it{p}_{T}^{MC}");
        fHistMCTrackMultGen -> Sumw2();

        // Correlation between reconstructed tracks and reconstructed particles
        fHistMCMultRes = new THnSparseF("fHistMCMultRes","Effect of momentum resolution on reconstructed particle multiplicities",2,nBinsRespMat, minMultRespMat, maxMultRespMat);
        fHistMCMultRes->SetBinEdges(0,fBinsMult->GetArray());
        fHistMCMultRes->SetBinEdges(1,fBinsMult->GetArray());
        fHistMCMultRes->GetAxis(0)->SetTitle("#it{N}_{acc}");
        fHistMCMultRes->GetAxis(1)->SetTitle("#it{N}_{part}");
      }
  }

  fOutputList->Add(fEventCount);

  fOutputList->Add(fHistEvent);
  fOutputList->Add(fHistTrack);
  fOutputList->Add(fHistRelPtResoFromCov);

  if(fIsMC){
    fOutputList->Add(fHistMCTrackParticle);
    fOutputList->Add(fHistMCPtRes);
    fOutputList->Add(fHistMCRelPtReso);
    fOutputList->Add(fHistMCEtaRes);

    fOutputList->Add(fHistMCRecTrack);
    fOutputList->Add(fHistMCGenPrimTrack);
    fOutputList->Add(fHistMCRecPrimTrack);
    fOutputList->Add(fHistMCRecSecTrack);

    fOutputList->Add(fHistMCResponseMat);
    fOutputList->Add(fHistMCResponseMatTracks);
    fOutputList->Add(fHistMCMultPtGenerated);

    if(fIncludeCrosscheckHistos){
      fOutputList->Add(fHistMCTrackMultGen);
      fOutputList->Add(fHistMCMultRes);
      fOutputList->Add(fHistMCParticle);
    }
  }

  PostData(1, fOutputList);

  if(fIsESD) InitESDTrackCuts();
  if((fIs2013pA || fIs2015data) && !fUtils){fUtils = new AliAnalysisUtils();}

}

/// Destructor
AliMeanPtAnalysisTask::~AliMeanPtAnalysisTask(){
  if(fUtils){delete fUtils; fUtils = NULL;}
  if(fESDtrackCuts){delete fESDtrackCuts; fESDtrackCuts = NULL;}
}

///________________________________________________________________________
void AliMeanPtAnalysisTask::UserExec(Option_t *){ // Main loop (called for each event)

  /// ====================== Initialize variables ===============================

  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler){Printf("ERROR: Could not receive inputHandler"); return;}

  AliPhysicsSelection* physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
  if(!physicsSelection) {Printf("ERROR: Could not receive physicsSelection"); return;}

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fEvent) {Printf("ERROR: fEvent not available\n"); return;}

  if(fIsMC){
    fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!fMCEvent) {Printf("ERROR: fMCEvent not available\n"); return;}
  }

  Bool_t isEventTriggered = inputHandler->IsEventSelected() & GetTriggerMask();

  Double_t multAccTracks = 0;   	 /// N_acc
  Double_t multGenPart = 0;  		   /// N_ch
  Double_t multRecPart = 0;		     /// N_rec


  /// ==================== Fill Histogramms ======================================

  fEventCount->Fill(1); // Generated Events

  if(!isEventTriggered) return;

  fEventCount->Fill(2); // Triggered Events

  if(!IsEventVertexOK(fEvent)) return;
  if (fIs2013pA){	if(!IsEventAccepted2013pA(fEvent)) return;	}
  if (fIs2015data){	if(!IsEventAccepted2015data(fEvent)) return;	}  // Requiring IncompleteDAQ, SPD background and Pileup cuts

  fEventCount->Fill(3); // Events after Vertex Cut and Pileup rejection

  Double_t centrality = 50;
  if((fBinsCent->GetSize()-1) > 1) centrality = GetCentrality(fEvent);

  /// ------------------ Count Multiplicities --------------------------------------

  // True Multiplicity Nch:
  if(fIsMC){
    for(Int_t iGenPart = 1; iGenPart < fMCEvent->GetNumberOfTracks(); iGenPart++) {
      AliMCParticle* mcGenParticle  = (AliMCParticle*)fMCEvent->GetTrack(iGenPart);
      if(!mcGenParticle) {Printf("ERROR: mcGenParticle  not available\n"); continue;}
      if(!IsParticleInKinematicRange(mcGenParticle)) continue;
      if(IsChargedPrimary(iGenPart)) multGenPart++;
    }
  }

  // Measured Multiplicity Nacc:
  AliVTrack* track = NULL;
  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++){
    track = fEvent->GetVTrack(iTrack);
    if (!track){Printf("ERROR: Could not receive track %d\n", iTrack); continue;}
    if(!IsTrackInKinematicRange(track)) continue;
    if(!IsTrackAcceptedQuality(track)) continue;
    multAccTracks++;
  }

  if(fIsMC){
    // Response Matrix
    Double_t responseMatrixTuple[2] = {multAccTracks, multGenPart};
    fHistMCResponseMat->Fill(responseMatrixTuple);

    if(multGenPart >= fBinsMult->GetAt(fBinsMult->GetSize()-1)) return;
    fEventCount->Fill(4); // Events after excluding Nch overflow
  }

  /// ------------------ Event Histogram ---------------------------------------

  Double_t eventValues[2] = {multAccTracks, centrality};
  fHistEvent->Fill(eventValues);

  ///--------------- Loop over measured Tracks ---------------------------------

  track = NULL;
  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++){
    track = fEvent->GetVTrack(iTrack);
    if(!track) {Printf("ERROR: Could not receive track %d\n", iTrack); continue;}

    if(!IsTrackInKinematicRange(track)) continue;
    if(!IsTrackAcceptedQuality(track)) continue;

    Double_t trackValues[4] = {track->Pt(), track->Eta(), multAccTracks, centrality};
    fHistTrack->Fill(trackValues);


    ///TODO this only works for ESD tracks!!
    Double_t ptResoValues[3] = {1./TMath::Abs(dynamic_cast<AliESDtrack*>(track)->GetSigned1Pt())*TMath::Sqrt(dynamic_cast<AliESDtrack*>(track)->GetSigma1Pt2()), track->Pt(), centrality};
    fHistRelPtResoFromCov->Fill(ptResoValues);

    /// Find original particle in MC-Stack
    if(fIsMC){
      Int_t mcLabel = TMath::Abs(track->GetLabel()); // negative label means bad quality track
      AliMCParticle* mcParticle  = (AliMCParticle*)fMCEvent->GetTrack(mcLabel);
      if(!mcParticle) {Printf("ERROR: mcParticle not available\n"); continue;}

      Double_t ptResValues[2] = {track->Pt(), mcParticle->Pt()};
      fHistMCPtRes->Fill(ptResValues);
      Double_t etaResValues[2] = {track->Eta(), mcParticle->Eta()};
      fHistMCEtaRes->Fill(etaResValues);

      if(!IsParticleInKinematicRange(mcParticle)) continue;
      multRecPart++;

      Double_t relPtReso[4] = {TMath::Abs(track->Pt() - mcParticle->Pt())/mcParticle->Pt(), mcParticle->Pt(), track->Pt(), centrality};
      fHistMCRelPtReso->Fill(relPtReso);

      if(fIncludeCrosscheckHistos){
        // Analogon to fHistTrack but with MC truth information
        Double_t particleValues[4] = {mcParticle->Pt(), mcParticle->Eta(), multAccTracks, centrality};
        fHistMCParticle->Fill(particleValues);
      }

      // Histograms for efficiency x acceptance correction
      Double_t mcRecTrackValue[3] = {mcParticle->Pt(), mcParticle->Eta(), centrality};
      fHistMCRecTrack->Fill(mcRecTrackValue);

      if(IsChargedPrimary(mcLabel))
      {
        Double_t mcPrimTrackValue[3] = {mcParticle->Pt(), mcParticle->Eta(), centrality};
        fHistMCRecPrimTrack->Fill(mcPrimTrackValue);

        Double_t responseMatrixTracksTuple[3] = {multAccTracks, multGenPart, mcParticle->Pt()};
        fHistMCResponseMatTracks->Fill(responseMatrixTracksTuple);

      }else{
        Double_t mcSecTrackValue[3] = {mcParticle->Pt(), mcParticle->Eta(), centrality};
        fHistMCRecSecTrack->Fill(mcSecTrackValue);
      }
    }
  }

  if(fIsMC){
    // Control Histogram to check how many tracks come from a particle outside of kinematic range
    fHistMCTrackParticle->Fill(1, multAccTracks);
    fHistMCTrackParticle->Fill(2, multRecPart);

    if(fIncludeCrosscheckHistos){
      Double_t multResoTuple[2] = {multAccTracks, multRecPart};
      fHistMCMultRes->Fill(multResoTuple);
    }
  }


  ///------------------- Loop over Generated Tracks (True MC)------------------------------
  if (fIsMC){

    for(Int_t iGenPart = 1; iGenPart < fMCEvent->GetNumberOfTracks(); iGenPart++) {
      AliMCParticle* mcGenParticle  = (AliMCParticle*)fMCEvent->GetTrack(iGenPart);
      if(!mcGenParticle) {Printf("ERROR: mcGenParticle  not available\n"); continue;}

      if(!IsParticleInKinematicRange(mcGenParticle)) continue;

      if(IsChargedPrimary(iGenPart)){

        Double_t mcGenPrimTrackValue[3] = {mcGenParticle->Pt(), mcGenParticle->Eta(), centrality};
        fHistMCGenPrimTrack->Fill(mcGenPrimTrackValue);

      	// Reference histogram for MC closure test and generator comparison
    	  Double_t mcMultPtGenerated[3] = {multGenPart, mcGenParticle->Pt(), centrality};
    	  fHistMCMultPtGenerated->Fill(mcMultPtGenerated);

        if(fIncludeCrosscheckHistos){
      	  Double_t trackValuesMult[3] = {multAccTracks, multGenPart, mcGenParticle->Pt()};
      	  fHistMCTrackMultGen->Fill(trackValuesMult);
        }
      }
    }
  }
  PostData(1, fOutputList);
}

//____________________________________________________________________________________________
void AliMeanPtAnalysisTask::Terminate(Option_t*)
{

}

Bool_t AliMeanPtAnalysisTask::IsChargedPrimary(Int_t mcLabel)
{
  if(!fMCEvent->IsPhysicalPrimary(mcLabel)) return kFALSE;
  AliMCParticle* mcParticle  = (AliMCParticle*)fMCEvent->GetTrack(mcLabel);
  if(!mcParticle) {Printf("ERROR: mcGenParticle  not available\n"); return kFALSE;}
  if(!(TMath::Abs(mcParticle->Charge()) > 0.01)) return kFALSE;
  return kTRUE;
}


/// Function implementing Track Acceptance cuts for tracks.
Bool_t AliMeanPtAnalysisTask::IsTrackInKinematicRange(AliVTrack* track)
{
  if(!track) return kFALSE;

  Double_t eta = track->Eta();
  Double_t pt  = track->Pt();

  if(eta <= fMinEta + fPRECISION)  return kFALSE;
  if(eta >= fMaxEta - fPRECISION)  return kFALSE;
  if(pt  <= fMinPt  + fPRECISION)  return kFALSE;
  if(pt  >= fMaxPt  - fPRECISION)  return kFALSE;
  return kTRUE;
}


/// Function implementing Track Acceptance cuts for MC particles.
Bool_t AliMeanPtAnalysisTask::IsParticleInKinematicRange(AliMCParticle* mcParticle)
{
  if(!mcParticle) return kFALSE;

  Double_t eta = mcParticle->Eta();
  Double_t pt  = mcParticle->Pt();

  if(eta <= fMinEta + fPRECISION)  return kFALSE;
  if(eta >= fMaxEta - fPRECISION)  return kFALSE;
  if(pt  <= fMinPt  + fPRECISION)  return kFALSE;
  if(pt  >= fMaxPt  - fPRECISION)  return kFALSE;
  return kTRUE;
}

/// Function implementing Track Quality cuts.
Bool_t AliMeanPtAnalysisTask::IsTrackAcceptedQuality(AliVTrack* track){
  if(fIsESD){
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (track);
    if(!fESDtrackCuts->AcceptTrack(esdTrack)) return kFALSE;
  }
  return kTRUE;
}

/// Function for event selection and pileup rejection for 2013 p-A
Bool_t AliMeanPtAnalysisTask::IsEventAccepted2013pA(AliVEvent* event)
{
  if (fUtils->IsFirstEventInChunk(event)) { return kFALSE;  }
  if (!fUtils->IsVertexSelected2013pA(event)) { return kFALSE;  }
  if (fUtils->IsPileUpEvent(event)) { return kFALSE;  }
  return kTRUE;
}

/// Function for Event cuts and pileup rejection for 2015 data: pp and Pb-Pb
Bool_t AliMeanPtAnalysisTask::IsEventAccepted2015data(AliVEvent* event)
{
  AliESDEvent* ESDevent = dynamic_cast<AliESDEvent*>(event);
  // if (fUtils->IsFirstEventInChunk(event)) { return kFALSE;  }
  // if (!fUtils->IsVertexSelected2013pA(event)) { return kFALSE;  }
  if (ESDevent->IsIncompleteDAQ()) { return kFALSE; }
  if (fUtils->IsSPDClusterVsTrackletBG(event)) { return kFALSE; }
  if (ESDevent->IsPileupFromSPD(5,0.8)) {return kFALSE; }
  return kTRUE;
}

/// Function for Event Acceptance cuts.
Bool_t AliMeanPtAnalysisTask::IsEventVertexOK(AliVEvent* event)
{
  if(!event) return kFALSE;
  if(TMath::Abs(event->GetPrimaryVertex()->GetZ()) > fZvtx) return kFALSE;
  if(!IsVertexOK(event)) return kFALSE;
  return kTRUE;
}


/// Function to get event centrality based on V0 measurement
Double_t AliMeanPtAnalysisTask::GetCentrality(AliVEvent* event)
{
  Double_t centrality = -1;
  AliMultSelection* multSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
  if(!multSelection){AliInfo("ERROR: No MultSelection found!"); return 999;}
  centrality = multSelection->GetMultiplicityPercentile("V0M"/*, lEmbedEventSelection = kFALSE*/);
  if(centrality > 100) {AliInfo("ERROR: Centrality determination does not work proprely!"); return 999;}
  return centrality;
}

/// Function to initialize the ESD track cuts object
void AliMeanPtAnalysisTask::InitESDTrackCuts(){

  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
  if(!fESDtrackCuts) {printf("ERROR: fESDtrackCuts not available\n"); return;}

  fESDtrackCuts->SetRequireTPCRefit(fTPCRefit);
  fESDtrackCuts->SetRequireITSRefit(fITSRefit);
  fESDtrackCuts->SetAcceptKinkDaughters(fAcceptKinks);
  if(fMinNCrossedRowsTPC > 0) fESDtrackCuts->SetMinNCrossedRowsTPC(fMinNCrossedRowsTPC);
  if(fMinRatioCrossedRowsOverFindableClustersTPC > 0) fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(fMinRatioCrossedRowsOverFindableClustersTPC);
  if(fMaxFractionSharedClustersTPC > 0) fESDtrackCuts->SetMaxFractionSharedTPCClusters(fMaxFractionSharedClustersTPC);
  if(fMaxChi2PerTPCCluster > 0) fESDtrackCuts->SetMaxChi2PerClusterTPC(fMaxChi2PerTPCCluster);
  if(fRequiresClusterITS) fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  if(!fRequiresClusterITS)fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
  if(fMaxChi2PerITSCluster > 0) fESDtrackCuts->SetMaxChi2PerClusterITS(fMaxChi2PerITSCluster);
  if(fDCAToVertex2D > 0) fESDtrackCuts->SetDCAToVertex2D(fDCAToVertex2D);
  if(fSigmaToVertex > 0) fESDtrackCuts->SetRequireSigmaToVertex(fSigmaToVertex);
  if(fMaxDCAzITSTPC > 0) fESDtrackCuts->SetMaxDCAToVertexZ(fMaxDCAzITSTPC);
  if(fDCAToVertexXY > 0) fESDtrackCuts->SetMaxDCAToVertexXY(fDCAToVertexXY);
  if(fDCAToVertexXYPtDep)          fESDtrackCuts->SetMaxDCAToVertexXYPtDep(fDCAToVertexXYPtDep);
  if(fMaxChi2TPCConstrained > 0) fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(fMaxChi2TPCConstrained);
  if(fMinActiveLength> 0) fESDtrackCuts->SetMinLengthActiveVolumeTPC(fMinActiveLength);
  if(fUseGeomCut) fESDtrackCuts->SetCutGeoNcrNcl(fDeadZoneWidth,fCutGeoNcrNclLenght,fCutGeoNcrNclGeom1Pt,fCutGeoNcrNclFractionNcr,fCutGeoNcrNclFractionNcl);

}

///Function to check vertex quality
Bool_t AliMeanPtAnalysisTask::IsVertexOK(AliVEvent* event){
  if(fIsESD){
    Double_t requiredZResolution = 1000;
    AliESDEvent* ESDevent = dynamic_cast<AliESDEvent*>(event);
    const AliESDVertex* esdVertex = ESDevent->GetPrimaryVertexTracks();
    if(!esdVertex){Printf("ERROR: vertex not available\n"); return kFALSE;}
    if(esdVertex->GetNContributors() < 1) {
      // SPD vertex
      esdVertex = ESDevent->GetPrimaryVertexSPD();
    }
    //     AliESDVertex *esdVertex = dynamic_cast<AliESDVertex*> (vertex);
    if(!esdVertex->GetStatus()){return kFALSE;}
    Double_t zRes = esdVertex->GetZRes();
    if (zRes > requiredZResolution) return kFALSE;

    const AliESDVertex* vertexSPD = ESDevent->GetPrimaryVertexSPD();
    // always check for SPD vertex
    const AliESDVertex* trkVertex = ESDevent->GetPrimaryVertexTracks();
    if(!vertexSPD) return kFALSE;
    if(!vertexSPD->GetStatus()) return kFALSE;
    if(!trkVertex->GetStatus()) return kFALSE;
    if (vertexSPD->IsFromVertexerZ() && !(vertexSPD->GetDispersion() < 0.04 && vertexSPD->GetZRes() < 0.25)) return kFALSE;
    //if (vertexSPD->IsFromVertexerZ() && vertexSPD->GetDispersion() > 0.04) return kFALSE; /// vertexSPD->GetDispersion() > 0.02 to 0.04
    if ((TMath::Abs(vertexSPD->GetZ() - trkVertex->GetZ()) > 0.5)) return kFALSE;
  }else{
    // TODO: AOD code goes here
  }
  return kTRUE;
}

void AliMeanPtAnalysisTask::SetFixedBinEdges(Double_t* array, Double_t lowerEdge, Double_t upperEdge, Int_t nBins){
  for(Int_t i = 0; i <= nBins; i++){
    array[i] = lowerEdge + i*(upperEdge - lowerEdge)/nBins;
  }
}
