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
#include "AlidNdPtEventCuts.h"

#include "AliPhysicsSelection.h"
#include "AliVMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"

#include "AliMCEvent.h"
#include "AliStack.h"

#include "AlidNdPtUnifiedAnalysisTask.h"


/// \cond CLASSIMP
ClassImp(AlidNdPtUnifiedAnalysisTask);
/// \endcond



//________________________________________________________________________
AlidNdPtUnifiedAnalysisTask::AlidNdPtUnifiedAnalysisTask(const char *name) : AliAnalysisTaskSE(name),
  //General member variables
  fOutputList(0),
  fEvent(0),
  fMCEvent(0),
  fMCStack(0),
  fFunTrkEff(0),
  fEventCuts(0),
  fESDtrackCuts(0),
  fUtils(0),
  //Toggles
  fIsESD(kTRUE),
  fUseMultiplicity(kTRUE),
  fIsMC(kFALSE),
  fEventTriggerRequired(kTRUE),
  fIs2013pA(kFALSE),
  fIs2015data(kFALSE),
  fUseTOFBunchCrossing(kFALSE),
  fTPCRefit(kFALSE),
  fITSRefit(kFALSE),
  fAcceptKinks(kTRUE),
  fRequiresClusterITS(kTRUE),
  fDCAToVertex2D(kFALSE),
  fSigmaToVertex(kFALSE),
  fUseGeomCut(kFALSE),
  fUseCountedMult(kFALSE),
  //Event-Histograms
  fHistEvent(0),
  fHistMultEvent(0),
  fHistMCGenEvent(0),
  fHistMCRecEvent(0),
  fHistMultMCGenEvent(0),		//rename fHistMCGenMultEvent
  fHistMultMCRecEvent(0),		//rename fHistMCRecMultEvent
  fHistMCGenINEL0Event(0),
  fHistMCRecINEL0Event(0),
  fHistMultCorrelation(0),		//rename fHistMCMultCorrelationEvent
  fHistMCResponseMat(0),
  fHistMCTrigEvent(0),
  fHistMCTrigINEL0Event(0),
  //Track-Histograms
  fHistTrack(0),
  fHistTrackCharge(0),		//rename fHistChargeTrack?
  fHistMCRecTrack(0),
  fHistMCRecTrackParticle(0),		//rename fHistMCRecPIDTrack
  fHistMCGenPrimTrack(0),
  fHistMCRecPrimTrack(0),
  fHistMCGenPrimTrackParticle(0),	//rename fHistMCGenPrimPIDTrack
  fHistMCRecPrimTrackParticle(0),	//rename fHistMCRecPrimPIDTrack
  fHistMCRecSecTrack(0),
  fHistMCRecSecTrackParticle(0),	//rename fHistMCRecSecPIDTrack
  //fHistMCGenTrackINEL0(0),
  // Cut Parameters
  fTriggerMask(AliVEvent::kMB),
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
  fBinsMultCent(0),
  fBinsPt(0),
  fBinsEta(0),
  fBinsZv(0),
  fHistMCMultPt(0),
  fUseCentralityCut(kFALSE),
  fLowerCentralityBound(0.),
  fUpperCentralityBound(0.),
  fIncludeSigmas(kTRUE),
  fHistV0Amp(0)
{
  // Set default binning
  Double_t binsMultCentDefault[2] = {0,10000};
  Double_t binsPtDefault[69] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};
  Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
  Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
  SetBinsPt(68,binsPtDefault);
  SetBinsEta(30,binsEtaDefault);
  SetBinsMultCent(1,binsMultCentDefault);
  SetBinsZv(12,binsZvDefault);
  SetMeanXYZv(0.0,0.0,0.0);
  SetSigmaMeanXYZv(1.0,1.0,10.0);
  SetZvtx(10.);
  SetEventTriggerRequired(kTRUE);
  TF1 *constant = new TF1("constant","1",0,100);
  SetTrkEffParametrisation(constant);
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AlidNdPtUnifiedAnalysisTask::UserCreateOutputObjects(){
  // Create histograms here (function is called once)
  OpenFile(1,"recreate");
  fOutputList = new TList();
  fOutputList -> SetOwner();

  /// Standard track histogram pt:eta:zV:multcent
  Int_t nBinsTrack[4]={fBinsPt->GetSize()-1,fBinsEta->GetSize()-1,fBinsZv->GetSize()-1,fBinsMultCent->GetSize()-1};
  Double_t minTrack[4]={fBinsPt->GetAt(0),fBinsEta->GetAt(0),fBinsZv->GetAt(0),fBinsMultCent->GetAt(0)};
  Double_t maxTrack[4]={fBinsPt->GetAt(fBinsPt->GetSize()-1),fBinsEta->GetAt(fBinsEta->GetSize()-1),fBinsZv->GetAt(fBinsZv->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1)};

  /// Control histogram charged patricles pt:eta:multcent:charge
  Int_t nBinsTrackCharge[4]={fBinsPt->GetSize()-1,fBinsEta->GetSize()-1,fBinsMultCent->GetSize()-1,2};
  Double_t minTrackCharge[4]={fBinsPt->GetAt(0),fBinsEta->GetAt(0),fBinsMultCent->GetAt(0),-2.0};
  Double_t maxTrackCharge[4]={fBinsPt->GetAt(fBinsPt->GetSize()-1),fBinsEta->GetAt(fBinsEta->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1),2.0};

  /// Control histogram charged patricles pt:eta:multcent:charge]
  const Int_t iNumberOfParticles = 11;
  Int_t nBinsTrackParticle[4]={fBinsPt->GetSize()-1,fBinsEta->GetSize()-1,fBinsMultCent->GetSize()-1,iNumberOfParticles};
  Double_t minTrackParticle[4]={fBinsPt->GetAt(0),fBinsEta->GetAt(0),fBinsMultCent->GetAt(0),0};
  Double_t maxTrackParticle[4]={fBinsPt->GetAt(fBinsPt->GetSize()-1),fBinsEta->GetAt(fBinsEta->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1),iNumberOfParticles};
  TString binNameTrackParticle[iNumberOfParticles]={"Pion","Kaon","Proton","SigmaPlus", "SigmaMinus", "OmegaMinus", "XiMinus", "Electron", "Muon", "Rest", "Lambda"};

  /// Standard event histogram zV:multcent
  Int_t nBinsEvent[2]={fBinsZv->GetSize()-1,fBinsMultCent->GetSize()-1};
  Double_t minEvent[2]={fBinsZv->GetAt(0),fBinsMultCent->GetAt(0)};
  Double_t maxEvent[2]={fBinsZv->GetAt(fBinsZv->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1)};

  /// Event multiplicity investigation histograms multcent:multacc:multcorr
  Int_t nBinsMultEvent[3]={fBinsMultCent->GetSize()-1,fBinsMultCent->GetSize()-1,fBinsMultCent->GetSize()-1};
  Double_t minMultEvent[3]={fBinsMultCent->GetAt(0),fBinsMultCent->GetAt(0),fBinsMultCent->GetAt(0)};
  Double_t maxMultEvent[3]={fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1)};

    
  /// Event multiplicity investigation histograms multcent:multacc
  Int_t nBinsMultEventCorrelation[2]={fBinsMultCent->GetSize()-1,fBinsMultCent->GetSize()-1};
  Double_t minMultEventCorrelation[2]={fBinsMultCent->GetAt(0),fBinsMultCent->GetAt(0)};
  Double_t maxMultEventCorrelation[2]={fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1)};
  

  /// Closure Test histogram (pt vs. Nch)
  Int_t nBinsMCMultPt[2]={fBinsMultCent->GetSize()-1,fBinsPt->GetSize()-1};
  Double_t minMCMultPt[2]={fBinsMultCent->GetAt(0),fBinsPt->GetAt(0)};
  Double_t maxMCMultPt[2]={fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1),fBinsPt->GetAt(fBinsPt->GetSize()-1)};
  
  
  fHistTrack = new THnF("fHistTrack", "Histogram for Tracks",4,nBinsTrack,minTrack,maxTrack);
  fHistTrack -> SetBinEdges(0,fBinsPt->GetArray());
  fHistTrack -> SetBinEdges(1,fBinsEta->GetArray());
  fHistTrack -> SetBinEdges(2,fBinsZv->GetArray());
  fHistTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
  fHistTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fHistTrack->GetAxis(1)->SetTitle("#eta");
  fHistTrack->GetAxis(2)->SetTitle("Zv (cm)");
  fHistTrack->GetAxis(3)->SetTitle("multiplicity (multCuts)");
  fHistTrack -> Sumw2();

  fHistTrackCharge = new THnF("fHistTrackCharge", "Control histogram for charged tracks",4,nBinsTrackCharge,minTrackCharge,maxTrackCharge);
  fHistTrackCharge -> SetBinEdges(0,fBinsPt->GetArray());
  fHistTrackCharge -> SetBinEdges(1,fBinsEta->GetArray());
  fHistTrackCharge -> SetBinEdges(2,fBinsMultCent->GetArray());
  fHistTrackCharge->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fHistTrackCharge->GetAxis(1)->SetTitle("#eta");
  fHistTrackCharge->GetAxis(2)->SetTitle("multiplicity (multCuts)");
  fHistTrackCharge->GetAxis(3)->SetTitle("track charge");
  fHistTrackCharge -> Sumw2();

  fHistEvent = new THnF("fHistEvent", "Histogram for Events",2,nBinsEvent,minEvent,maxEvent);
  fHistEvent -> SetBinEdges(0,fBinsZv->GetArray());
  fHistEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
  fHistEvent->GetAxis(0)->SetTitle("Zv (cm)");
  fHistEvent->GetAxis(1)->SetTitle("multiplicity (multCuts)");
  fHistEvent -> Sumw2();

  fHistMultEvent = new THnF("fHistMultEvent", "Histogram for MultEvents",3,nBinsMultEvent,minMultEvent,maxMultEvent);
  fHistMultEvent -> SetBinEdges(0,fBinsMultCent->GetArray());
  fHistMultEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
  fHistMultEvent -> SetBinEdges(2,fBinsMultCent->GetArray());
  fHistMultEvent->GetAxis(0)->SetTitle("Event multipicity");
  fHistMultEvent->GetAxis(1)->SetTitle("Accepted track multiplicity");
  fHistMultEvent->GetAxis(2)->SetTitle("Corrected track multiplicity");
  fHistMultEvent -> Sumw2();

// temporary histogram  
  fHistV0Amp = new TH1D("V0Amp", "V0Amp",2000,0,200);
  fHistV0Amp -> Sumw2();

  
  if(fIsMC){

    fHistMCGenPrimTrack = new THnF("fHistMCGenPrimTrack", "Histogram for generated MC Tracks",4,nBinsTrack,minTrack,maxTrack);
    fHistMCGenPrimTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCGenPrimTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCGenPrimTrack -> SetBinEdges(2,fBinsZv->GetArray());
    fHistMCGenPrimTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
    fHistMCGenPrimTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCGenPrimTrack->GetAxis(1)->SetTitle("#eta");
    fHistMCGenPrimTrack->GetAxis(2)->SetTitle("Zv (cm)");
    fHistMCGenPrimTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
    fHistMCGenPrimTrack -> Sumw2();

    fHistMCRecTrack = new THnF("fHistMCRecTrack", "Histogram for reconstructed MC Tracks",4,nBinsTrack,minTrack,maxTrack);
    fHistMCRecTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecTrack -> SetBinEdges(2,fBinsZv->GetArray());
    fHistMCRecTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
    fHistMCRecTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCRecTrack->GetAxis(1)->SetTitle("#eta");
    fHistMCRecTrack->GetAxis(2)->SetTitle("Zv (cm)");
    fHistMCRecTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
    fHistMCRecTrack -> Sumw2();

    fHistMCRecPrimTrack = new THnF("fHistMCRecPrimTrack", "Histogram for reconstructed primary MC Tracks",4,nBinsTrack,minTrack,maxTrack);
    fHistMCRecPrimTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecPrimTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecPrimTrack -> SetBinEdges(2,fBinsZv->GetArray());
    fHistMCRecPrimTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
    fHistMCRecPrimTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCRecPrimTrack->GetAxis(1)->SetTitle("#eta");
    fHistMCRecPrimTrack->GetAxis(2)->SetTitle("Zv (cm)");
    fHistMCRecPrimTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
    fHistMCRecPrimTrack -> Sumw2();

    fHistMCRecSecTrack = new THnF("fHistMCRecSecTrack", "Histogram for reconstructed secondary MC Tracks",4,nBinsTrack,minTrack,maxTrack);
    fHistMCRecSecTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecSecTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecSecTrack -> SetBinEdges(2,fBinsZv->GetArray());
    fHistMCRecSecTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
    fHistMCRecSecTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCRecSecTrack->GetAxis(1)->SetTitle("#eta");
    fHistMCRecSecTrack->GetAxis(2)->SetTitle("Zv (cm)");
    fHistMCRecSecTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
    fHistMCRecSecTrack -> Sumw2();

    fHistMCGenPrimTrackParticle = new THnF("fHistMCGenPrimTrackParticle", "Control histogram for charged tracks",4,nBinsTrackParticle,minTrackParticle,maxTrackParticle);
    fHistMCGenPrimTrackParticle -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCGenPrimTrackParticle -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCGenPrimTrackParticle -> SetBinEdges(2,fBinsMultCent->GetArray());
    fHistMCGenPrimTrackParticle->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCGenPrimTrackParticle->GetAxis(1)->SetTitle("#eta");
    fHistMCGenPrimTrackParticle->GetAxis(2)->SetTitle("multiplicity (multCuts)");
    fHistMCGenPrimTrackParticle->GetAxis(3)->SetTitle("Particle type");
    fHistMCGenPrimTrackParticle -> Sumw2();

    fHistMCRecTrackParticle = new THnF("fHistMCRecTrackParticle", "Control histogram for charged tracks",4,nBinsTrackParticle,minTrackParticle,maxTrackParticle);
    fHistMCRecTrackParticle -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecTrackParticle -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecTrackParticle -> SetBinEdges(2,fBinsMultCent->GetArray());
    fHistMCRecTrackParticle->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCRecTrackParticle->GetAxis(1)->SetTitle("#eta");
    fHistMCRecTrackParticle->GetAxis(2)->SetTitle("multiplicity (multCuts)");
    fHistMCRecTrackParticle->GetAxis(3)->SetTitle("Particle type");
    fHistMCRecTrackParticle -> Sumw2();

    fHistMCRecPrimTrackParticle = new THnF("fHistMCRecPrimTrackParticle", "Control histogram for charged tracks",4,nBinsTrackParticle,minTrackParticle,maxTrackParticle);
    fHistMCRecPrimTrackParticle -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecPrimTrackParticle -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecPrimTrackParticle -> SetBinEdges(2,fBinsMultCent->GetArray());
    fHistMCRecPrimTrackParticle->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCRecPrimTrackParticle->GetAxis(1)->SetTitle("#eta");
    fHistMCRecPrimTrackParticle->GetAxis(2)->SetTitle("multiplicity (multCuts)");
    fHistMCRecPrimTrackParticle->GetAxis(3)->SetTitle("Particle type");
    fHistMCRecPrimTrackParticle -> Sumw2();

    fHistMCRecSecTrackParticle = new THnF("fHistMCRecSecTrackParticle", "Control histogram for charged tracks",4,nBinsTrackParticle,minTrackParticle,maxTrackParticle);
    fHistMCRecSecTrackParticle -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecSecTrackParticle -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecSecTrackParticle -> SetBinEdges(2,fBinsMultCent->GetArray());
    fHistMCRecSecTrackParticle->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCRecSecTrackParticle->GetAxis(1)->SetTitle("#eta");
    fHistMCRecSecTrackParticle->GetAxis(2)->SetTitle("multiplicity (multCuts)");
    fHistMCRecSecTrackParticle->GetAxis(3)->SetTitle("Particle type");
    fHistMCRecSecTrackParticle -> Sumw2();

    for( Int_t ii  = 1; ii<= fHistMCRecSecTrackParticle->GetAxis(3)->GetNbins(); ii++){
      fHistMCGenPrimTrackParticle->GetAxis(3)->SetBinLabel(ii,binNameTrackParticle[ii-1].Data());
      fHistMCRecTrackParticle->GetAxis(3)->SetBinLabel(ii,binNameTrackParticle[ii-1].Data());
      fHistMCRecPrimTrackParticle->GetAxis(3)->SetBinLabel(ii,binNameTrackParticle[ii-1].Data());
      fHistMCRecSecTrackParticle->GetAxis(3)->SetBinLabel(ii,binNameTrackParticle[ii-1].Data());
    }

    fHistMCGenEvent = new THnF("fHistMCGenEvent", "Histogram for generated MC Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCGenEvent -> SetBinEdges(0,fBinsZv->GetArray());
    fHistMCGenEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCGenEvent->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCGenEvent->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCGenEvent -> Sumw2();

    fHistMultMCGenEvent = new THnF("fHistMultMCGenEvent", "Mult histogram for rec MC events",3,nBinsMultEvent,minMultEvent,maxMultEvent);
    fHistMultMCGenEvent -> SetBinEdges(0,fBinsMultCent->GetArray());
    fHistMultMCGenEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMultMCGenEvent -> SetBinEdges(2,fBinsMultCent->GetArray());
    fHistMultMCGenEvent->GetAxis(0)->SetTitle("Event multipicity");
    fHistMultMCGenEvent->GetAxis(1)->SetTitle("Accepted track multiplicity");
    fHistMultMCGenEvent->GetAxis(2)->SetTitle("Corrected track multiplicity");
    fHistMultMCGenEvent -> Sumw2();

    fHistMCTrigEvent = new THnF("fHistMCTrigEvent", "Histogram for triggered MC Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCTrigEvent -> SetBinEdges(0,fBinsZv->GetArray());
    fHistMCTrigEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCTrigEvent->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCTrigEvent->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCTrigEvent -> Sumw2();

    fHistMCRecEvent = new THnF("fHistMCRecEvent", "Histogram for reconstructed MC  Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCRecEvent -> SetBinEdges(0,fBinsZv->GetArray());
    fHistMCRecEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCRecEvent->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCRecEvent->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCRecEvent -> Sumw2();

    fHistMultMCRecEvent = new THnF("fHistMultMCRecEvent", "Mult histogram for gen MC events",3,nBinsMultEvent,minMultEvent,maxMultEvent);
    fHistMultMCRecEvent -> SetBinEdges(0,fBinsMultCent->GetArray());
    fHistMultMCRecEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMultMCRecEvent -> SetBinEdges(2,fBinsMultCent->GetArray());
    fHistMultMCRecEvent->GetAxis(0)->SetTitle("Event multipicity");
    fHistMultMCRecEvent->GetAxis(1)->SetTitle("Accepted track multiplicity");
    fHistMultMCRecEvent->GetAxis(2)->SetTitle("Corrected track multiplicity");
    fHistMultMCRecEvent -> Sumw2();

    fHistMultCorrelation = new THnF("fHistMultCorrelation", "Histogram for MultEventCorrelations",2,nBinsMultEventCorrelation,minMultEventCorrelation,maxMultEventCorrelation);
    fHistMultCorrelation -> SetBinEdges(0,fBinsMultCent->GetArray());
    fHistMultCorrelation -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMultCorrelation -> SetBinEdges(2,fBinsMultCent->GetArray());
    fHistMultCorrelation->GetAxis(0)->SetTitle("Gen. Corrected track multiplicity");
    fHistMultCorrelation->GetAxis(1)->SetTitle("Rec. Corrected track multiplicity");
    fHistMultCorrelation -> Sumw2();


    fHistMCGenINEL0Event = new THnF("fHistMCGenINEL0Event", "Histogram for generated INEL>0 MC Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCGenINEL0Event -> SetBinEdges(0,fBinsZv->GetArray());
    fHistMCGenINEL0Event -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCGenINEL0Event->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCGenINEL0Event->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCGenINEL0Event -> Sumw2();

    fHistMCTrigINEL0Event = new THnF("fHistMCTrigINEL0Event","Histogram for triggered INEL>0 MC Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCTrigINEL0Event->SetBinEdges(0,fBinsZv->GetArray());
    fHistMCTrigINEL0Event->SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCTrigINEL0Event->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCTrigINEL0Event->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCTrigINEL0Event->Sumw2();

    fHistMCRecINEL0Event = new THnF("fHistMCRecINEL0Event","Histogram for reconstructed INEL>0 MC Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCRecINEL0Event->SetBinEdges(0,fBinsZv->GetArray());
    fHistMCRecINEL0Event->SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCRecINEL0Event->GetAxis(0)->SetTitle("mcZv (cm)");
    fHistMCRecINEL0Event->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCRecINEL0Event->Sumw2();

           
    fHistMCResponseMat = new THnF("fHistMCResponseMat","Histogram for MC Response Matrix N_{ch} vs. N_{acc}",2,nBinsMultEventCorrelation, minMultEventCorrelation, maxMultEventCorrelation);
    fHistMCResponseMat->SetBinEdges(0,fBinsMultCent->GetArray());
    fHistMCResponseMat->SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCResponseMat->GetAxis(0)->SetTitle("reconstructed track multiplicity N_{acc}");
    fHistMCResponseMat->GetAxis(1)->SetTitle("generated particle multiplicity N_{ch}");
    fHistMCResponseMat->Sumw2();
    
    fHistMCMultPt = new THnF("fHistMCMultPt","Histogram for MC Closure test #it{p}_T vs. N_{acc}",2,nBinsMCMultPt, minMCMultPt, maxMCMultPt);
    fHistMCMultPt->SetBinEdges(0,fBinsMultCent->GetArray());
    fHistMCMultPt->SetBinEdges(1,fBinsPt->GetArray());
    fHistMCMultPt->GetAxis(0)->SetTitle("N_{ch}");
    fHistMCMultPt->GetAxis(1)->SetTitle("#it{p}_T");
    fHistMCMultPt->Sumw2();


    /*
       fHistMCGenTrackINEL0 = new THnF("fHistMCGenTrackINEL0","Histogram for generated tracks for INEL>0 MC Events",4,nBinsTrack,minTrack,maxTrack);
       fHistMCGenTrackINEL0->SetBinEdges(0,fBinsPt->GetArray());
       fHistMCGenTrackINEL0->SetBinEdges(1,fBinsEta->GetArray());
       fHistMCGenTrackINEL0->SetBinEdges(2,fBinsZv->GetArray());
       fHistMCGenTrackINEL0->SetBinEdges(3,fBinsMultCent->GetArray());
       fHistMCGenTrackINEL0->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
       fHistMCGenTrackINEL0->GetAxis(1)->SetTitle("#eta");
       fHistMCGenTrackINEL0->GetAxis(2)->SetTitle("Zv (cm)");
       fHistMCGenTrackINEL0->GetAxis(3)->SetTitle("true multiplicity (MC)");
       fHistMCGenTrackINEL0 -> Sumw2();
       */
  }

  fOutputList->Add(fHistTrack);
  fOutputList->Add(fHistEvent);
  fOutputList->Add(fHistMultEvent);
  fOutputList->Add(fHistTrackCharge);

  fOutputList->Add(fHistV0Amp);


  if(fIsMC){
    fOutputList->Add(fHistMCGenPrimTrack);
    fOutputList->Add(fHistMCRecTrack);
    fOutputList->Add(fHistMCRecPrimTrack);
    fOutputList->Add(fHistMCRecSecTrack);

    fOutputList->Add(fHistMCGenPrimTrackParticle);
    fOutputList->Add(fHistMCRecTrackParticle);
    fOutputList->Add(fHistMCRecPrimTrackParticle);
    fOutputList->Add(fHistMCRecSecTrackParticle);

    fOutputList->Add(fHistMCGenEvent);
    fOutputList->Add(fHistMultMCGenEvent);
    fOutputList->Add(fHistMCTrigEvent);
    fOutputList->Add(fHistMCRecEvent);
    fOutputList->Add(fHistMultMCRecEvent);
    fOutputList->Add(fHistMultCorrelation);
    fOutputList->Add(fHistMCGenINEL0Event);
    fOutputList->Add(fHistMCTrigINEL0Event);
    fOutputList->Add(fHistMCRecINEL0Event);
    fOutputList->Add(fHistMCResponseMat);
    fOutputList->Add(fHistMCMultPt);
    
    
    //     fOutputList->Add(fHistMCGenTrackINEL0);
  }
  PostData(1, fOutputList);

  /// Create and initialize Analysis Objects here instead of in UserExec() to save resources
  //InitdNdPtEventCuts(); (does nothing atm so I just leave it out for the moment)
  if(fIsESD) InitESDTrackCuts();
  if((fIs2013pA || fIs2015data) && !fUtils){fUtils = new AliAnalysisUtils();}

}

/// Destructor
AlidNdPtUnifiedAnalysisTask::~AlidNdPtUnifiedAnalysisTask(){
  if(fUtils){delete fUtils; fUtils=0;}
  if(fESDtrackCuts){delete fESDtrackCuts; fESDtrackCuts=0;}
  //if(fEventCuts){delete fEventCuts; fEventCuts=0;}
}

///________________________________________________________________________
void AlidNdPtUnifiedAnalysisTask::UserExec(Option_t *){ // Main loop (called for each event)


  /// ====================== Initialize variables ===============================

  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler){ Printf("ERROR: Could not receive inputHandler"); return; }

  AliPhysicsSelection *physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
  if(!physicsSelection) {Printf("ERROR: Could not receive physicsSelection"); return;}
  //AliTriggerAnalysis* triggerAnalysis = NULL;

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fEvent) {printf("ERROR: fEvent not available\n"); return;}

  Bool_t isEventINEL0 = kFALSE;
  if(fIsMC){
    fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!fMCEvent) {printf("ERROR: fMCEvent not available\n"); return;}

    fMCStack = fMCEvent->Stack();
    if (!fMCStack) {printf("ERROR: fMCStack not available\n"); return;}

    // Check if MC event is in inel0 class (1 charged particle in abs(eta)<1.0, pt>0)
    isEventINEL0 = IsMCEventINEL0(fMCEvent,0,1.0);
  }

  Bool_t isEventTriggered = inputHandler->IsEventSelected() & GetTriggerMask();
  Double_t multEvent = GetEventMultCent(fEvent);
  Double_t zVertEvent = fEvent->GetPrimaryVertex()->GetZ();
  Double_t eventValues[2] = {zVertEvent, multEvent};

  AliVVZERO * vZeroHandler = fEvent->GetVZEROData();
  if (!vZeroHandler) {printf("ERROR: vZeroHandler not available\n"); return;}
  Double_t v0Mult = (Double_t) vZeroHandler->GetMTotV0A();
  
  // take a look at the range of V0 amplitude
  fHistV0Amp->Fill(v0Mult);

  Double_t multAccTracks = 0;   	/// N_acc (of tracks!!)
  Double_t multAccCorrTracks = 0;	/// N_acc (of tracks!!) corrected with trk efficiency

  Double_t multGenPart = 0;  		/// N_ch
  Double_t multRecPart = 0;		/// N_acc (of tracks that can be assigned to a real particle)
  Double_t multRecCorrPart = 0;		/// N_acc corrected with trk efficiency

  // for PbPb only use specified centrality interval because particle composition correction depends on centrality
  if(fUseCentralityCut && !IsSelectedCentrality()){ return; }

  /// ==================== Fill Histogramms ======================================

  /// -------------------- Generated Events --------------------------------------

  if(fIsMC){
    fHistMCGenEvent->Fill(eventValues);
    if(isEventINEL0) fHistMCGenINEL0Event->Fill(eventValues);
  }

  /// \li Event Trigger Cut
  if(!isEventTriggered) return;

  /// -------------------- Triggered Events ---------------------------------------

  if(fIsMC){
    fHistMCTrigEvent->Fill(eventValues);
    if(isEventINEL0) fHistMCTrigINEL0Event->Fill(eventValues);
  }

  /// \li Event Acceptance Cuts
  if(!IsEventAcceptedGeometrics(fEvent)) return;
  if(!IsEventAcceptedQuality(fEvent)) return;
  if (fIs2013pA){	if(!IsEventAccepted2013pA(fEvent)) return;	}
  if (fIs2015data){	if(!IsEventAccepted2015data(fEvent)) return;	}


  /// ------------------ Reconstructed Events --------------------------------------

  if(fIsMC){
    fHistMCRecEvent->Fill(eventValues);
    if(isEventINEL0) fHistMCRecINEL0Event->Fill(eventValues);
  }

  ///--------------- Loop over measured Tracks ---------------------------------

  /// Multiplicity counting loop
  AliVTrack *track = NULL;
  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++){
    track = fEvent->GetVTrack(iTrack);
    if (!track){ printf("ERROR: Could not receive track %d\n", iTrack); continue; }

    /// \li Track Acceptance Cuts
    if(!IsTrackAcceptedKinematics(track)) continue;
    if(!IsTrackAcceptedQuality(track)) continue;

    if(fUseTOFBunchCrossing && fIs2015data){
      if(TMath::Abs(track->GetTOFsignalDz())>10) continue;
      if((track->GetTOFsignal())<12000) continue;
      if((track->GetTOFsignal())>25000) continue;
    }


    /// \li Count Track Multiplicity
    multAccTracks++;

    //correct multiplicity via pt-dependant tracking efficiency function
    Double_t dTrackingEff = fFunTrkEff->Eval(track->Pt());
    if (dTrackingEff != 0) multAccCorrTracks = multAccCorrTracks + 1/dTrackingEff;
    else multAccCorrTracks = multAccCorrTracks + 1;

  }
  
  // from now on use actual multipicity in all histograms
  if (fUseCountedMult) {multEvent = multAccTracks; eventValues[1] = multAccTracks;}

  fHistEvent->Fill(eventValues);


  track = NULL;
  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++){
    track = fEvent->GetVTrack(iTrack);
    if (!track){ printf("ERROR: Could not receive track %d\n", iTrack); continue; }

    /// \li Track Acceptance Cuts
    if(!IsTrackAcceptedKinematics(track)) continue;
    if(!IsTrackAcceptedQuality(track)) continue;

    if(fUseTOFBunchCrossing && fIs2015data){
      if(TMath::Abs(track->GetTOFsignalDz())>10) continue;
      if((track->GetTOFsignal())<12000) continue;
      if((track->GetTOFsignal())>25000) continue;
    }


    /// \li Fill Track Histograms

    Double_t trackValues[4] = {track->Pt(), track->Eta(), zVertEvent, multEvent};
    fHistTrack->Fill(trackValues);

    Double_t trackChargeValues[4] = {track->Pt(), track->Eta(), multEvent, ((Double_t) track->Charge())/3.};
    fHistTrackCharge->Fill(trackChargeValues);


    /// \li Find original particle in MC-Stack
    if(fIsMC){
      Int_t mcLabel = TMath::Abs(track->GetLabel());	///TODO doesn't work without abs for some reason?? source of errors!!!
      TParticle *mcParticle = fMCStack->Particle(mcLabel);
      if(!mcParticle) {printf("ERROR: mcParticle not available\n"); continue;}

      //is original particle also fulfilling the required conditions (not only its reconstructed track)?
      if(!IsTrackAcceptedKinematics(mcParticle)) continue;

      Double_t mcRecTrackValue[4] = {mcParticle->Pt(), mcParticle->Eta(), zVertEvent, multEvent};
      fHistMCRecTrack->Fill(mcRecTrackValue);

      Double_t mcRecTrackParticleValue[4] = {mcParticle->Pt(), mcParticle->Eta(), multEvent,((Double_t) IdentifyMCParticle(mcLabel)) - 0.5};
      fHistMCRecTrackParticle->Fill(mcRecTrackParticleValue);

      if(IsChargedPrimary(mcLabel))
      {
        Double_t mcPrimTrackValue[4] = {mcParticle->Pt(), mcParticle->Eta(), zVertEvent, multEvent};
        fHistMCRecPrimTrack->Fill(mcPrimTrackValue);
	
        multRecPart++;

        Double_t dTrackingEff = fFunTrkEff->Eval(mcParticle->Pt());
        if (dTrackingEff != 0) multRecCorrPart = multRecCorrPart + 1/dTrackingEff;
        else multRecCorrPart = multRecCorrPart + 1;
      }else{//works because tracks are always charged
        Double_t mcSecTrackValue[4] = {mcParticle->Pt(), mcParticle->Eta(), zVertEvent, multEvent};
        fHistMCRecSecTrack->Fill(mcSecTrackValue);

        Double_t mcRecSecTrackParticleValue[4] = {mcParticle->Pt(), mcParticle->Eta(), multEvent,((Double_t) IdentifyMCParticle(mcLabel)) - 0.5};
        fHistMCRecSecTrackParticle->Fill(mcRecSecTrackParticleValue);
      }
      if(IsChargedPrimaryOrLambda(mcLabel))
      {
        Double_t mcRecPrimTrackParticleValue[4] = {mcParticle->Pt(), mcParticle->Eta(), multEvent, ((Double_t) IdentifyMCParticle(mcLabel)) - 0.5};
        fHistMCRecPrimTrackParticle->Fill(mcRecPrimTrackParticleValue);
      }
    }
  }// end of Track-loop

  //Fill histogram with overall,accepted and corrected multiplicity
  Double_t eventMultValues[3] = {multEvent, multAccTracks, multAccCorrTracks};
  fHistMultEvent->Fill(eventMultValues);

  if (fIsMC){
    Double_t eventRecMultValues[3] = {multEvent, multRecPart, multRecCorrPart};
    fHistMultMCRecEvent->Fill(eventRecMultValues);
  }

  ///------------------- Loop over Generated Tracks (True MC)------------------------------
  if (fIsMC){

    /// count generated multipicity
    for (Int_t iParticle = 0; iParticle < fMCStack->GetNtrack(); iParticle++){
      TParticle *mcGenParticle = fMCStack->Particle(iParticle);
      if(!mcGenParticle) {printf("ERROR: mcGenParticle  not available\n"); continue;}

      /// \li Acceptance cuts for generated particles
      // lower pt cut is disabled for mpt analysis! (Nch should be counted down to pt=0)!
      if(!IsTrackAcceptedKinematics(mcGenParticle, kTRUE)) continue;

      if(IsChargedPrimary(iParticle)){
	
        multGenPart++;

      }
    }

    
    for (Int_t iParticle = 0; iParticle < fMCStack->GetNtrack(); iParticle++){
      TParticle *mcGenParticle = fMCStack->Particle(iParticle);
      if(!mcGenParticle) {printf("ERROR: mcGenParticle  not available\n"); continue;}

      /// \li Acceptance cuts for generated particles
      // lower pt cut is disabled for mpt analysis! (Nch should be counted down to pt=0)!
      if(!IsTrackAcceptedKinematics(mcGenParticle, kTRUE)) continue;

      if(IsChargedPrimary(iParticle)){

        Double_t mcGenPrimTrackValue[4] = {mcGenParticle->Pt(), mcGenParticle->Eta(), zVertEvent, multEvent};
        fHistMCGenPrimTrack->Fill(mcGenPrimTrackValue);

	// multPt hist for MC closure test
	if(multEvent > 0.1){
	  Double_t mcMultPt[2] = {multGenPart, mcGenParticle->Pt()};
	  fHistMCMultPt->Fill(mcMultPt);
	}

      }
      if(IsChargedPrimaryOrLambda(iParticle))
      {
        Double_t mcGenPrimTrackParticleValue[4] = {mcGenParticle->Pt(), mcGenParticle->Eta(), multEvent,((Double_t) IdentifyMCParticle(iParticle)) - 0.5};
        fHistMCGenPrimTrackParticle->Fill(mcGenPrimTrackParticleValue);
      }
    }

    Double_t eventGenMultValues[3] = {multEvent, multGenPart, multGenPart};//?
    fHistMultMCGenEvent->Fill(eventGenMultValues);

    Double_t eventMultCorrelValues[2] = {multGenPart, multRecCorrPart};
    fHistMultCorrelation->Fill(eventMultCorrelValues);

    Double_t responseMatrixTuple[2] = {multEvent, multGenPart};
    fHistMCResponseMat->Fill(responseMatrixTuple);
    
    

  }
    
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AlidNdPtUnifiedAnalysisTask::Terminate(Option_t *)
{

}

Bool_t AlidNdPtUnifiedAnalysisTask::IsChargedPrimary(Int_t stackIndex){
  if (fMCStack->IsPhysicalPrimary(stackIndex) && (TMath::Abs(fMCStack->Particle(stackIndex)->GetPDG()->Charge()) > 0.01)){
   
    // In case charged particles are defined as charged particles without the sigmas exclude them
    Int_t pID = IdentifyMCParticle(stackIndex);
    if(!fIncludeSigmas && (pID == kSigmaPlus || pID == kSigmaMinus)) return kFALSE;
    
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AlidNdPtUnifiedAnalysisTask::IsChargedPrimaryOrLambda(Int_t stackIndex){
  if (fMCStack->IsPhysicalPrimary(stackIndex) && (TMath::Abs(fMCStack->Particle(stackIndex)->GetPDG()->Charge()) > 0.01 || TMath::Abs(fMCStack->Particle(stackIndex)->GetPdgCode())==3122)) return kTRUE;
  return kFALSE;
}

/// Track Acceptance cuts, for tracks.
///
/// \param AliVTrack Input track
///
/// \return Is track accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedKinematics(AliVTrack *track)
{
  if(!track) return kFALSE;

  Float_t eta = track->Eta();	///TODO why float?
  Float_t pt = track->Pt();

  if(eta < fMinEta) return kFALSE;
  if(eta > fMaxEta) return kFALSE;
  if(pt < fMinPt) return kFALSE;
  if(pt > fMaxPt) return kFALSE;

  return kTRUE;
}

/// Track Acceptance cuts, for MC particles.
///
/// \param TParticle Input particle
///
/// \return Is particle accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedKinematics(TParticle *mcTrack, Bool_t useLowerPtCut)
{
  if(!mcTrack) return kFALSE;

  Float_t eta = mcTrack->Eta();		///TODO why float?
  Float_t pt = mcTrack->Pt();

  if(eta < fMinEta) return kFALSE;
  if(eta > fMaxEta) return kFALSE;
  if((pt < fMinPt) && useLowerPtCut) return kFALSE;
  if(pt > fMaxPt) return kFALSE;
  return kTRUE;
}

/// Track Quality cuts.
///
/// \param AliVTrack Input track
///
/// \return Is track accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedQuality(AliVTrack *track){
  if(fIsESD){
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*> (track);
    if(!fESDtrackCuts->AcceptTrack(esdTrack)) return kFALSE;
  }
  return kTRUE;
}

/// selection and pileup rejection for 2013 p-A
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAccepted2013pA(AliVEvent *event)
{
  if (fUtils->IsFirstEventInChunk(event)) { return kFALSE;  }
  if (!fUtils->IsVertexSelected2013pA(event)) { return kFALSE;  }
  if (fUtils->IsPileUpEvent(event)) { return kFALSE;  }
  return kTRUE;
}

/// Event cuts and pileup rejection for 2015 data: pp and Pb-Pb
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAccepted2015data(AliVEvent *event)
{
  AliESDEvent* ESDevent = dynamic_cast<AliESDEvent*>(event);
  // if (fUtils->IsFirstEventInChunk(event)) { return kFALSE;  }
  // if (!fUtils->IsVertexSelected2013pA(event)) { return kFALSE;  }
  if (ESDevent->IsIncompleteDAQ()) { return kFALSE; }
  if (fUtils->IsSPDClusterVsTrackletBG(event)) { return kFALSE; }
  if (ESDevent->IsPileupFromSPD(5,0.8)) {return kFALSE; }
  return kTRUE;
}

/// Event Acceptance cuts.
///
/// \param AliVEvent Input event
///
/// \return Is event accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAcceptedGeometrics(AliVEvent *event)
{
  // when fEventCuts is used (i.e. something happens in InitdNdPtEventCuts), we could use  something like if(!fEventCuts->AcceptEvent)
  // maybe even without this* function...
  if(TMath::Abs(event->GetPrimaryVertex()->GetZ())>fZvtx) return kFALSE;
  return kTRUE;
}

/// Event Quality cuts.
///
/// \param AliVEvent Input event
///
/// \return Is event accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAcceptedQuality(AliVEvent *event)
{
  //   AliVVertex *vertex = event->GetPrimaryVertexTracks();
  if(!event) return kFALSE;
  if(!IsVertexOK(event)) return kFALSE;
  return kTRUE;
}

/// Function to either fill Multiplicity or centrality
///
/// \param AliVEvent event to be analised
///
/// \return Double_t with centrality or multiplicity
Double_t AlidNdPtUnifiedAnalysisTask::GetEventMultCent(AliVEvent *event)
{
  if(fUseMultiplicity)
  {
    AliVMultiplicity* multiplicity = event->GetMultiplicity();
    if(!multiplicity) {printf("ERROR: multiplicity not available\n"); return 999;}
    Int_t mult = multiplicity->GetNumberOfTracklets();
    return mult;
  }
  else
  {
    Float_t centralityF = -1;
    AliMultSelection *MultSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
    if ( MultSelection ){
      centralityF = MultSelection->GetMultiplicityPercentile("V0M"/*, lEmbedEventSelection = kFALSE*/);
      if(centralityF>100) return 999;
      return centralityF;
    }else{
      AliInfo("Didn't find MultSelection!");
    }
  }
}

/// Function obtain multipicity percentile
Bool_t AlidNdPtUnifiedAnalysisTask::IsSelectedCentrality(){

    // centrality determination
    Float_t centralityF = -1.;
    AliMultSelection *MultSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");

    if (MultSelection){
      centralityF = MultSelection->GetMultiplicityPercentile("V0M");	
      if(centralityF >= fLowerCentralityBound  && centralityF <= fUpperCentralityBound) return kTRUE; 
    }
    else{Printf("ERROR: Could not receive mult selection"); AliInfo("Didn't find MultSelection!");}

    return kFALSE;
}

/// Function to initialize the ESD track cuts
void AlidNdPtUnifiedAnalysisTask::InitESDTrackCuts(){

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

/// Function to initialize the dNdPtEventCuts (CURRENTLY not used...)
void AlidNdPtUnifiedAnalysisTask::InitdNdPtEventCuts(){

  fEventCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  if(!fEventCuts) {printf("ERROR: fEventCuts not available\n"); return;}

  fEventCuts->SetMeanXYZv(fMeanXYZv[0],fMeanXYZv[1],fMeanXYZv[2]);
  fEventCuts->SetSigmaMeanXYZv(fSigmaMeanXYZv[0],fSigmaMeanXYZv[1],fSigmaMeanXYZv[2]);
  fEventCuts->SetTriggerRequired(fEventTriggerRequired);
}

///Function to check vertex quality
Bool_t AlidNdPtUnifiedAnalysisTask::IsVertexOK(AliVEvent *event){
  if(fIsESD){
    Float_t requiredZResolution = 1000;
    AliESDEvent* ESDevent = dynamic_cast<AliESDEvent*>(event);
    const AliESDVertex *esdVertex = ESDevent->GetPrimaryVertexTracks();
    if(!esdVertex){printf("ERROR: vertex not available\n"); return kFALSE;}
    if(esdVertex->GetNContributors()<1) {
      // SPD vertex
      esdVertex = ESDevent->GetPrimaryVertexSPD();
    }
    //     AliESDVertex *esdVertex = dynamic_cast<AliESDVertex*> (vertex);
    if(!esdVertex->GetStatus()){return kFALSE;}
    Double_t zRes = esdVertex->GetZRes();
    if (zRes > requiredZResolution) return kFALSE;

    const AliESDVertex *vertexSPD = ESDevent->GetPrimaryVertexSPD();
    // always check for SPD vertex
    if(!vertexSPD) return kFALSE;
    if(!vertexSPD->GetStatus()) return kFALSE;
    if (vertexSPD->IsFromVertexerZ() && vertexSPD->GetDispersion() > 0.04) return kFALSE; /// vertexSPD->GetDispersion() > 0.02 to 0.04
  }else{
    //AOD code goes here
  }
  return kTRUE;
}

/// Function to determine if MC Event is in INEL>0 Class
/// select INEL>0 events with at least
/// one prompt (MC primary) particle in acceptance
/// pT>0, |eta|<1.0 for normalization
Bool_t AlidNdPtUnifiedAnalysisTask::IsMCEventINEL0(AliMCEvent* mcEvent, Double_t ptmin, Double_t etarange){
  if(!mcEvent) return kFALSE;
  // AliMCEvent* fmcEvent = static_cast<AliMCEvent*>(mcEvent);
  AliStack* stack = mcEvent->Stack();
  if(!stack) return kFALSE;

  Int_t count = 0;
  for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc)
  {
    TParticle* particle = stack->Particle(iMc);
    if (!particle) continue;

    // only charged particles
    if(!particle->GetPDG()) continue;
    Double_t charge = particle->GetPDG()->Charge()/3.;
    if(charge == 0) continue;

    // physical primary
    Bool_t prim = stack->IsPhysicalPrimary(iMc);
    if(!prim) continue;

    if(particle->Pt() < ptmin) continue;
    if(TMath::Abs(particle->Eta()) > etarange) continue;

    count++;
  }

  if(count > 0) return kTRUE;
  else return kFALSE;
}

///TODO Not implemented yet
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedGeometricalCut(AliVTrack *tr, Double_t bMagZ){
  return kTRUE;
}

/// Function to return Particle ID for Histograms
Int_t AlidNdPtUnifiedAnalysisTask::IdentifyMCParticle(Int_t mcLabel){

  Int_t ipdg = TMath::Abs(fMCStack->Particle(mcLabel)->GetPdgCode()); // Abs() because antiparticles are negaitve...
  if(ipdg==211)  return AlidNdPtUnifiedAnalysisTask::kPion;
  if(ipdg==321)  return AlidNdPtUnifiedAnalysisTask::kKaon;
  if(ipdg==2212) return AlidNdPtUnifiedAnalysisTask::kProtons;
  if(ipdg==3222) return AlidNdPtUnifiedAnalysisTask::kSigmaPlus;
  if(ipdg==3112) return AlidNdPtUnifiedAnalysisTask::kSigmaMinus;
  if(ipdg==3334) return AlidNdPtUnifiedAnalysisTask::kOmegaMinus;
  if(ipdg==3312) return AlidNdPtUnifiedAnalysisTask::kXiMinus;
  if(ipdg==11) return AlidNdPtUnifiedAnalysisTask::kElectron;
  if(ipdg==13) return AlidNdPtUnifiedAnalysisTask::kMuon;
  if(ipdg==3122) return AlidNdPtUnifiedAnalysisTask::kLambda;
  return AlidNdPtUnifiedAnalysisTask::kRest;
}
