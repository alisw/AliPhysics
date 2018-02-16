  // $Id$
  //
  // Track analysis base task.
  //
  // Authors: D.Lodato



#include <TClonesArray.h>
#include <TChain.h>
#include <TList.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THnSparse.h>
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliClusterContainer.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVCaloCells.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliTrackContainer.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliEmcalParticle.h"
#include "AliParticleContainer.h"
#include "AliPicoTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEMCALRecoUtils.h"
#include "AliLog.h"
#include "TF1.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "AliGenPythiaEventHeader.h"


#include "AliAnalysisTaskEMCALTrackIsolation.h"

  /// \cond CLASSIMP
ClassImp(AliAnalysisTaskEMCALTrackIsolation);
  /// \endcond

using std::cout;
using std::endl;
  //________________________________________________________________________
AliAnalysisTaskEMCALTrackIsolation::AliAnalysisTaskEMCALTrackIsolation() :
AliAnalysisTaskEmcal("AliAnalysisTaskEMCALTrackIsolation",kTRUE),
fAOD(0),
fVevent(0),
fNTracks(0),
fAODMCParticles(0),
fmcHeader(0),
fTracksAna(0),
fStack(0),
fEMCALRecoUtils(new AliEMCALRecoUtils),
fIsoConeRadius(0.4),
fptIsoMethod(0),
fptIsoThreshold(2),
fIsoMethod(0),
fdetacut(0.025),
fdphicut(0.03),
fQA(0),
fUEMethod(0),
fIsMC(0),
fTPC4Iso(0),
fNDimensions(0),
fMCDimensions(0),
fMCQAdim(0),
fisLTAnalysis(kFALSE),
fRejectionEventWithoutTracks(kTRUE),
fAnalysispPb(kFALSE),
fTriggerLevel1(0),
fTest1(0),
fTest2(0),
fMCtruth(0),
fBinsPt(),
fBinsPtiso(),
fBinsPtue(),
fBinsEta(),
fBinsPhi(),
fBinsLabel(),
fBinsPDG(),
fBinsMomPDG(),
fBinsTrackPDG(),
fBinsDx(),
fBinsDz(),
fBinsDecay(),
fTrackMult(0),
fEtaPhiTrack(0),
fPT(0),
fP(0),
fVz(0),
fEvents(0),
fPtaftTime(0),
fPtaftFC(0),
fTrackTime(0),
fPtIsoTrack(0),
fPhiBandUETracks(0),
fEtaBandUETracks(0),
fPerpConesUETracks(0),
fTPCWithoutIsoConeB2BbandUE(0),
fPtIsolatedNTracks(0),
fTestIndex(0),
fTestIndexPt(0),
fTestLocalIndexPt(0),
fTrackMultvsPt(0),
fHistXsection(0),
fHistTrials(0),
fPtTracksVSpTTR(0),
fPtTracksVSpTTR_MC(0),
fPhiTracksVStrackPt(0),
fEtaTracksVStrackPt(0),
fOutputTHnS(0),
fOutMCTruth(0),
fOutTrackMC(0)
{
  SetMakeGeneralHistograms(kTRUE);
}

  //________________________________________________________________________
AliAnalysisTaskEMCALTrackIsolation::AliAnalysisTaskEMCALTrackIsolation(const char *name, Bool_t histo) :
AliAnalysisTaskEmcal(name, histo),
fAOD(0),
fVevent(0),
fNTracks(0),
fAODMCParticles(0),
fmcHeader(0),
fTracksAna(0),
fStack(0),
fEMCALRecoUtils(new AliEMCALRecoUtils),
fIsoConeRadius(0.4),
fptIsoMethod(0),
fptIsoThreshold(2),
fIsoMethod(0),
fdetacut(0.025),
fdphicut(0.03),
fQA(0),
fUEMethod(0),
fIsMC(0),
fTPC4Iso(0),
fNDimensions(0),
fMCDimensions(0),
fMCQAdim(0),
fisLTAnalysis(kFALSE),
fRejectionEventWithoutTracks(kTRUE),
fAnalysispPb(kFALSE),
fTriggerLevel1(0),
fTest1(0),
fTest2(0),
fMCtruth(0),
fBinsPt(),
fBinsPtiso(),
fBinsPtue(),
fBinsEta(),
fBinsPhi(),
fBinsLabel(),
fBinsPDG(),
fBinsMomPDG(),
fBinsTrackPDG(),
fBinsDx(),
fBinsDz(),
fBinsDecay(),
fTrackMult(0),
fEtaPhiTrack(0),
fPT(0),
fP(0),
fVz(0),
fEvents(0),
fPtaftTime(0),
fPtaftFC(0),
fTrackTime(0),
fPtIsoTrack(0),
fPhiBandUETracks(0),
fEtaBandUETracks(0),
fPerpConesUETracks(0),
fTPCWithoutIsoConeB2BbandUE(0),
fPtIsolatedNTracks(0),
fTestIndex(0),
fTestIndexPt(0),
fTestLocalIndexPt(0),
fTrackMultvsPt(0),
fHistXsection(0),
fHistTrials(0),
fPtTracksVSpTTR(0),
fPtTracksVSpTTR_MC(0),
fPhiTracksVStrackPt(0),
fEtaTracksVStrackPt(0),
fOutputTHnS(0),
fOutMCTruth(0),
fOutTrackMC(0)
{
  SetMakeGeneralHistograms(kTRUE);
}
  //________________________________________________________________________
AliAnalysisTaskEMCALTrackIsolation::~AliAnalysisTaskEMCALTrackIsolation(){
    // Destructor
}


  //________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::UserCreateOutputObjects(){
    // Create ouput histograms and THnSparse and TTree
  
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
    //printf("Up here all good");
  
  fOutput = new AliEmcalList(); // RH: Leak? fOutput already exists in base class
  fOutput->SetOwner();
  
  TString sIsoMethod="\0",sUEMethod="\0",sBoundaries="\0";
  TString sTitle;
  
  if(fTPC4Iso)
    sBoundaries = "TPC Acceptance";
  else
    sBoundaries = "EMCAL Acceptance";
  
  if(fUEMethod==0)
    sUEMethod = "PhiBand";
  else if(fUEMethod==1)
    sUEMethod = "EtaBand";
  else if(fUEMethod==2)
    sUEMethod = "PerpCones";
  else if(fUEMethod==3)
    sUEMethod = "FullTPC";
  
  Int_t binPT = fBinsPt.size()-1;
  Int_t binPTiso = fBinsPtiso.size()-1;
  Int_t binPTUE = fBinsPtue.size()-1;
  Int_t binetatr= fBinsEta.size()-1;
  Int_t binphitr = fBinsPhi.size()-1;
  Int_t binlabel = fBinsLabel.size()-1;
  
  Int_t binMCPDG = fBinsPDG.size()-1;
  Int_t binMCMotherPDG = fBinsMomPDG.size()-1;
  Int_t binMCTrackPDG = fBinsTrackPDG.size()-1;
  Int_t bindx = fBinsDx.size()-1;
  Int_t bindz = fBinsDz.size()-1 ;
  Int_t binDecayType = fBinsDecay.size()-1;
    //bincells=20;
  
    //   Initialization of all the Common THistos for the Three different outputs
  fEvents = new TH1D("hEvents_TR","Events",100,0.,100.);
  fEvents->Sumw2();
  fOutput->Add(fEvents);
  
  fVz = new TH1D("hVz_TR","Vertex Z distribution",100,-50.,50.);
  fVz->Sumw2();
  fOutput->Add(fVz);
  
  fPtaftTime = new TH1D("hPtaftTime_TR","p_{T} distribution for tracks after tracks time cut",200,0.,100.);
  fPtaftTime->Sumw2();
  fOutput->Add(fPtaftTime);
  
  fPtaftFC = new TH1D("hPtaftFC_TR","p_{T} distribution for tracks after fiducial cut",200,0.,100.);
  fPtaftFC->Sumw2();
  fOutput->Add(fPtaftFC);
  
  fPtTracksVSpTTR = new TH2F ("hPtTracksVSpTTR","Charged Particle spectrum vs pT Candidate",70,0.,70.,200,0.,20.);
  fPtTracksVSpTTR->Sumw2();
  fOutput->Add(fPtTracksVSpTTR);

    //Include QA plots to the OutputList //DEFINE BETTER THE BINNING AND THE AXES LIMITS
  fTrackMult = new TH1D ("hTrackMult","Tracks multiplicity Distribution",100,0.,100.);
  fTrackMult->Sumw2();
  fOutput->Add(fTrackMult);
  
  fTrackTime = new TH1D("hTrackTime_TR","Time distribution for tracks",2000,1e4,3e4);
  fTrackTime->Sumw2();
  fOutput->Add(fTrackTime);
  
  fEtaPhiTrack = new TH2D ("hEtaPhiTrackActivity","",250,-0.8,0.8, 250, 1.2, 3.4);
  fEtaPhiTrack->Sumw2();
  fOutput->Add(fEtaPhiTrack);
  
  fPT = new TH1D("hPt_TR","P_{T} distribution for tracks",100,0.,100.);
  fPT->Sumw2();
  fOutput->Add(fPT);
  
  fP = new TH1D("hE_TR","E distribution for tracks",200,0.,100.);
  fP->Sumw2();
  fOutput->Add(fP);
  
  fTestIndexPt= new TH2D("hTestIndexPt","Test index vs transverse momentum for tracks",200,0.,100.,100,0.,100.);
  fTestIndexPt->SetXTitle("tracks pT");
  fTestIndexPt->SetYTitle("index");
  fTestIndexPt->Sumw2();
  fOutput->Add(fTestIndexPt);
  
  Int_t bins[] = {binPT, binPTiso, binPTUE, binetatr, binphitr};
  
  fNDimensions = sizeof(bins)/sizeof(Int_t);
  const Int_t ndims =   fNDimensions;
  
  sTitle = Form("Direct Photons: p_{T} , P_{T} Iso tracks in %s acceptance, P_{T} UE tracks in %s, #eta_{track} distr,#phi_{track} distr; p_{T} (GeV/c); P_{T}^{isoTracks} (GeV/c) ; P_{T}^{UE%s} (GeV/c); #eta_{tr}; #phi_{tr}", sBoundaries.Data(), sUEMethod.Data(), sUEMethod.Data());
  
  fOutputTHnS =  new THnSparseF("fHnOutput",sTitle.Data(), ndims, bins);
  fOutputTHnS->SetBinEdges(0,fBinsPt.data());
  fOutputTHnS->SetBinEdges(1,fBinsPtiso.data());
  fOutputTHnS->SetBinEdges(2,fBinsPtue.data());
  fOutputTHnS->SetBinEdges(3,fBinsEta.data());
  fOutputTHnS->SetBinEdges(4,fBinsPhi.data());
  fOutputTHnS->Sumw2();
  
  fOutput->Add(fOutputTHnS);
  
  if(fIsMC){
    Int_t binsMC[] = {binPT, binPTiso, binPTUE, binMCPDG ,binetatr,binphitr,binlabel};
    Int_t binsSMC[] = {binPT, binMCTrackPDG, binMCMotherPDG, binPT, bindx, bindz, binPTiso,binPTiso};
    
    fMCDimensions = sizeof(binsMC)/sizeof(Int_t);
    const Int_t ndimsMC = fMCDimensions;
    
    fOutMCTruth = new THnSparseF ("fOutMCTruth","E_{#gamma}, E_{T}^{iso cone}, E_{T}^{UE}, MomPDG, Eta, Phi, Label; E_{T}^{#gamma} (GeV/c); p_{T}^{Iso}(GeV/c);E_{T} ^{UE} (GeV/c); PDG; #eta; #phi; Label",ndimsMC,binsMC);
    fOutMCTruth->SetBinEdges(0,fBinsPt.data());
    fOutMCTruth->SetBinEdges(1,fBinsPtiso.data());
    fOutMCTruth->SetBinEdges(2,fBinsPtue.data());
    fOutMCTruth->SetBinEdges(3,fBinsPDG.data());
    fOutMCTruth->SetBinEdges(4,fBinsEta.data());
    fOutMCTruth->SetBinEdges(5,fBinsPhi.data());
    fOutMCTruth->SetBinEdges(6,fBinsLabel.data());
    fOutMCTruth->Sumw2();
    fOutput->Add(fOutMCTruth);
    
    fMCQAdim = sizeof(binsSMC)/sizeof(Int_t);
    const Int_t ndimsMCQA = fMCQAdim;
    
    Double_t xminbismix[] = {0., -3000, -400,  0.,-1., -1., -10.,-10.};
    Double_t xmaxbismix[] = {70., 3000,  400, 70., 1.,  1., 100.,100.};
    fOutTrackMC = new THnSparseF ("fOutTracktMC", "E_{T}^{track}, PDG, MOM PDG, E_{T}^{true}, #Deltax, #Deltaz, E_{T}^{iso},Label;E_{T}^{reco} (GeV/c); PDG Code; Mothers' PDG Code; E_{T}^{MCtrue} (GeV/c); #Delta#phi; #Delta#eta; E_{T}^{iso,reco} (Gev/c);E_{T}^{iso,truth} (Gev/c)",ndimsMCQA,binsSMC);
    fOutTrackMC->SetBinEdges(0,fBinsPt.data());
    fOutTrackMC->SetBinEdges(1,fBinsTrackPDG.data());
    fOutTrackMC->SetBinEdges(2,fBinsMomPDG.data());
    fOutTrackMC->SetBinEdges(3,fBinsPt.data());
    fOutTrackMC->SetBinEdges(4,fBinsDx.data());
    fOutTrackMC->SetBinEdges(5,fBinsDz.data());
    fOutTrackMC->SetBinEdges(6,fBinsPtiso.data());
    fOutTrackMC->SetBinEdges(7,fBinsPtiso.data());
    fOutTrackMC->Sumw2();
    fOutput->Add(fOutTrackMC);
  }
  fPtIsoTrack = new TH2D("hPtIsoTrack_TR"," #Sigma p_{T}^{iso cone} in iso cone distribution for Neutral tracks with Tracks",200,0.,100.,200,0.,100.);
  fPtIsoTrack->SetYTitle("#Sigma p_{T}^{iso cone} (GeV/c)");
  fPtIsoTrack->SetXTitle("p_{T}^{track}");
  fPtIsoTrack->Sumw2();
  fOutput->Add(fPtIsoTrack);
  
  fPhiBandUETracks = new TH2D(Form("hPhiBandUE_TPC"),Form("UE Estimation with Phi Band TPC "),200,0.,100.,250,0.,100.);
  fPhiBandUETracks->SetXTitle("E_{T}");
  fPhiBandUETracks->SetYTitle("#Sigma P_{T}^{UE}");
  fPhiBandUETracks->Sumw2();
  fOutput->Add(fPhiBandUETracks);
  
  fEtaBandUETracks = new TH2D(Form("hEtaBandUE_TPC"),Form("UE Estimation with Eta Band and TPC"),200,0.,100.,250,0.,100.);
  fEtaBandUETracks->SetXTitle("E_{T}");
  fEtaBandUETracks->SetYTitle("#Sigma P_{T}^{UE}");
  fEtaBandUETracks->Sumw2();
  fOutput->Add(fEtaBandUETracks);
  
  fPerpConesUETracks = new TH2D("hConesUE","UE Estimation with Perpendicular Cones in TPC",200,0.,100.,250,0.,100.);
  fPerpConesUETracks->SetXTitle("E_{T}");
  fPerpConesUETracks->SetYTitle("#Sigma P_{T}^{UE}");
  fPerpConesUETracks->Sumw2();
  fOutput->Add(fPerpConesUETracks);
  
  fTPCWithoutIsoConeB2BbandUE = new TH2D("hFullTPCUE","UE Estimation with almost Full TPC",200,0.,100.,250,0.,100.);
  fTPCWithoutIsoConeB2BbandUE->SetXTitle("E_{T}");
  fTPCWithoutIsoConeB2BbandUE->SetYTitle("#Sigma E_{T}^{UE}");
  fTPCWithoutIsoConeB2BbandUE->Sumw2();
  fOutput->Add(fTPCWithoutIsoConeB2BbandUE);
  
  fPtIsolatedNTracks = new TH1D("hEtIsolatedNTracks","p_{T} distribution for isolated tracks; #Sigma p_{T}^{iso cone}<Pthres",200,0.,100.);
  fPtIsolatedNTracks->SetXTitle("p_{T}");
  fPtIsolatedNTracks->Sumw2();
  fOutput->Add(fPtIsolatedNTracks);
    //
    // TO CHECK IF THERE ARE OTHER HISTOGRAMS TO BE INSERTED HERE
    //
  
  fTestIndex= new TH2D("hTestIndex","Test index for tracks",100,0.,100.,100,0.,100.);
  fTestIndex->SetXTitle("index");
  fTestIndex->SetYTitle("local index");
  fTestIndex->Sumw2();
  fOutput->Add(fTestIndex);
  
  fTestLocalIndexPt= new TH2D("hTestLocalIndexE","Test local index vs tranverse momentum for tracks",200,0.,100.,100,0.,100.);
  fTestLocalIndexPt->SetXTitle("tracks pT");
  fTestLocalIndexPt->SetYTitle("local index");
  fTestLocalIndexPt->Sumw2();
  fOutput->Add(fTestLocalIndexPt);
  
  fTrackMultvsPt = new TH2D("hTrackMultvsPt","Track Multiplicity vs  p_{T} track for clusters",100,0.,100.,200,0.,100.);
  fTrackMultvsPt->Sumw2();
  fOutput->Add(fTrackMultvsPt);
  
  
  if(fIsMC){
      //CREATE THE TH2 specific for the MC Analysis Maybe to be added in the THNSparse, or cloning two or three axes and add the specific MC Truth info
    fHistXsection = new TH1F("fHistXsection", "fHistXsection", 1, 0, 1);
    fHistXsection->GetXaxis()->SetBinLabel(1,"<sigma>");
    fOutput->Add(fHistXsection);
    
    fHistTrials = new TH1F("fHistTrials", "fHistTrials", 1, 0, 1);
    fHistTrials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
    fOutput->Add(fHistTrials);
    
    fPtTracksVSpTTR_MC = new TH2F ("hChargedptSpecVSpT_MC","Charged Particle spectrum vs pT Candidate",70,0.,70.,200,0.,50.);
    fPtTracksVSpTTR_MC->Sumw2();
    fOutput->Add(fPtTracksVSpTTR_MC);
  }
  PostData(1, fOutput);
    //     //   return;
}

  //________________________________________________________________________
Double_t* AliAnalysisTaskEMCALTrackIsolation::GenerateFixedBinArray(Int_t n, Double_t min, Double_t max) const
{
    // Generate the bin array for the ThnSparse
  
  Double_t *bins = new Double_t[n+1];
  
  Double_t binWidth = (max-min)/n;
  bins[0] = min;
  for (Int_t i = 1; i <= n; i++) {
    bins[i] = bins[i-1]+binWidth;
  }
  
  return bins;
}

  //________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::ExecOnce()
{
    //tracks for Isolation and analysis
  AliTrackContainer *tracksANA = GetTrackContainer("filterTracksAna");
    //  //Printf("name of the second track container: %s", tracksANA->GetClassName().Data());
  
  if (!tracksANA) {
    AliError(Form("%s: This task needs a 2particle container!", GetName()));
    return;
  }
  if(fIsoMethod==1){
    AliClusterContainer *clusters = GetClusterContainer(0);
    if (!clusters) {
      AliError(Form("%s: This task needs a cluster container!", GetName()));
      return;
    }
  }
    
    //Init the EMCAL Framework
  AliAnalysisTaskEmcal::ExecOnce();
  if (!fLocalInitialized) {
    
    AliError(Form("AliAnalysisTask not initialized"));
    return;
  }
}

  //______________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALTrackIsolation::SelectCandidate(AliAODTrack *aodToi)
{
  if(!(aodToi->IsHybridGlobalConstrainedGlobal())) return kFALSE;
  
    //Printf("Pt %f",aodToi->Pt());
  if((aodToi->Pt())<0.2)
    return kFALSE;
  
    //Printf("With ITSrefit? %s",((aodToi->GetStatus() & AliAODTrack::kITSrefit)==AliAODTrack::kITSrefit)?"yes":"no");
  if((aodToi->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit)
    return kFALSE;
  
  
  Double_t frac = 0;
  Float_t ncls  = Float_t(aodToi->GetTPCncls ());
  Float_t nclsS = Float_t(aodToi->GetTPCnclsS());
  if (  ncls> 0 )  frac =  nclsS / ncls ;
  
  if(frac > 0.4) return kFALSE;
  
  Int_t index=aodToi->GetID();
  aodToi->P();
  Double_t toiTOF = aodToi->GetTOFsignal();
  fTrackTime->Fill(toiTOF);
  
  if(fQA)  FillQAHistograms(aodToi);
  
  if(!fIsMC){
    if(toiTOF< 1.2e4 || toiTOF > 2.5e4)
      return kFALSE;
  }
  fPtaftTime->Fill(aodToi->Pt());
  
  if(!CheckBoundaries(aodToi))
    return kFALSE;
  
    //Printf("ok2!");
  fPtaftFC->Fill(aodToi->Pt());
  
  if(fQA){
    fTestIndexPt->Fill(aodToi->Pt(),index);
  }
    //Printf("ok3!");
  if(aodToi->Pt()<5.) return kFALSE;
  
  fEtaPhiTrack->Fill(aodToi->Eta(),aodToi->Phi());
  
  return kTRUE;
}
  //______________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALTrackIsolation::Run()
{
    // Run the analysis.
    //fTest1+=1;
    //vertex cuts
  AliTrackContainer *tracksANA = GetTrackContainer("filterTracksAna");
  if(!tracksANA){
      //Printf("Cannot find the tracks for Isolation");
    return kFALSE;
  }
  
    //  //Printf("FilterType of the tracks for Analysis: %d \t(should be %d)", tracksANA->GetTrackFilterType(),AliEmcalTrackSelection::kHybridTracks);
  
    //    AliError(Form("\n\n\n\nGO CHECK the Settings!!!! Is Isolation calculated with filteredTracks?\n\n\n\n"));
  if(tracksANA->GetTrackFilterType()!= AliEmcalTrackSelection::kHybridTracks){
    AliWarning(Form("Isolation NOT calculated with HybridTracks"));
    AliWarning(Form("You better be sure of what you are doing"));
  }
  
  
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  
    //Printf("Vertex Z coordinate for NF: %lf", fVertex[2]);
  
  if (fVertex[2]>10. || fVertex[2]<-10.) return kFALSE;
    //  AliError(Form("La task tourne bien"));
  
  Int_t nbTracksEvent;
  nbTracksEvent =InputEvent()->GetNumberOfTracks();
    //Printf("nbTracksEvent %d",nbTracksEvent);
  if(fRejectionEventWithoutTracks && nbTracksEvent==0)
    return kFALSE;
  
    // Fill events number histogram
  fEvents->Fill(0);
  
  Double_t minPhi= 0., maxPhi= 2*TMath::Pi(), minEta = -0.87, maxEta= 0.87;
  
  if(fIsMC){
    Float_t fXSection=0;
    Int_t fTrials=0;
    if(fPythiaHeader){
      fXSection = fPythiaHeader->GetXsection();
      fTrials = fPythiaHeader->Trials();
      if(fTrials>0){
        fHistXsection->Fill("<sigma>",fXSection);
        fHistTrials->Fill("#sum{ntrials}",fTrials);
      }
    }
      // AliError(Form("EMCAL L1 trigger for MC simulation anchored to LHC13 data"));
      //    if(!MCSimTrigger(fVevent,fTriggerLevel1) && fAnalysispPb) return kFALSE;
  }
  
    //Fill Vertex Z histogram
  fVz->Fill(fVertex[2]);
  
    // delete output USEFUL LATER FOR CONTAINER CREATION !!
    //fOuttracks->Delete();
  Int_t index=0;
  
  if(fIsMC){
    AliAODMCHeader *mcHeader;
    
    fAODMCParticles = static_cast <TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
    
    fmcHeader = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));
      //    //Printf("%d",fmcHeader->GetEventType());
    if (!fIsMC)
      return kFALSE;
      //AliInfo(Form("It's a MC analysis %e",fAODMCParticles));
    if(!fStack && !fAODMCParticles){
      AliError("no MC stack saved\n");
      return kFALSE;
    }
    
      //cout<<"there's a List of particles"<<endl;
      //DO THIS ALSO FOR ESDs
    
    if(fAODMCParticles->GetEntries() < 1){
      AliError("number of MC particles insufficient");
      return kFALSE;
    }
      //    //Printf("%d",fMCtruth);
      //    //Printf("Passe analyze MC");
    
      //    if(fMCtruth || (fmcHeader->GetEventType()==14 || fmcHeader->GetEventType()==29)){
      //      //Printf("Analysing mc");
    AnalyzeMC();
      //    }
  }
  
  if(fisLTAnalysis){
    AliVTrack *toi = (tracksANA->GetLeadingTrack());
    if(!toi) {
      AliError("No Leading track found");
      return kFALSE;
    }
    AliAODTrack *aodToi=static_cast<AliAODTrack*>(toi);
    
    Bool_t isSelected=SelectCandidate(aodToi);
    
    if(isSelected){
      for (auto it : tracksANA->accepted()){
        AliAODTrack *tr = static_cast<AliAODTrack*>(it);
        if(!tr)
          continue;
        
        fPtTracksVSpTTR->Fill(aodToi->Pt(),tr->Pt());
          //      cout<<"ok4.2"<<endl;
          //      fPhiTracksVStrackPt->Fill(aodToi->Pt(),tr->Phi());
          //      cout<<"ok4.3"<<endl;
          //      fEtaTracksVStrackPt->Fill(aodToi->Pt(),tr->Eta());
      }
        //Printf("ok5!");
      Int_t index=aodToi->GetID();
      
      FillGeneralHistograms(aodToi,index);
    }//is Selected Track
  }//Leading Track Analysis
  else{
    AliAODTrack *aodToi;
    for (auto it : tracksANA->accepted()){
      AliVTrack *toi = static_cast<AliVTrack*>(it);
      if(!toi) {
        AliError("No tracks found");
        return kFALSE;
      }
        //
      aodToi = static_cast<AliAODTrack*>(toi);
        //Printf("hybrid? %s", aodToi->IsHybridGlobalConstrainedGlobal()?"yes":"no");
        //Printf("ok4!");
        //Printf("Inside Run: Passing to FillGeneralHistograms for tracks with ID: %d, Pt: %.3f, Eta: %.3f, Phi: %.3f",index,aodToi->Pt(),aodToi->Eta(),aodToi->Phi());
      Bool_t isSelected=SelectCandidate(aodToi);
      
      if(isSelected){
        for (auto it : tracksANA->accepted()){
          AliAODTrack *tr = static_cast<AliAODTrack*>(it);
          if(!tr)
            continue;
          
          cout<<"ok4.1"<<endl;
          fPtTracksVSpTTR->Fill(aodToi->Pt(),tr->Pt());
            //      cout<<"ok4.2"<<endl;
            //      fPhiTracksVStrackPt->Fill(aodToi->Pt(),tr->Phi());
            //      cout<<"ok4.3"<<endl;
            //      fEtaTracksVStrackPt->Fill(aodToi->Pt(),tr->Eta());
        }
          //Printf("ok5!");
        Int_t index=aodToi->GetID();
        FillGeneralHistograms(aodToi,index);
      }//is Selected Track
    }//For Loop on Accepted Tracks
  }//NO Leading Track Analysis
  return kTRUE;
}

  //_________________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::FillQAHistograms(AliAODTrack *toi){
  
  fPT->Fill(toi->Pt());
  fP->Fill(toi->P());
  
}

  //___________________________________________________________________________________
  //Bool_t AliAnalysisTaskEMCALTrackIsolation::MCSimTrigger(AliVEvent *eventIn, Int_t triggerLevel){
  //
  //  if(triggerLevel<1 || triggerLevel>2){
  //    AliError(Form("No trigger level have been choosed for the MC analysis"));
  //    return kFALSE;
  //  }
  //  Double_t threshold=0;
  //  Double_t spread=0;
  //
  //    //Int_t runNumber = InputEvent()->GetRunNumber();
  //
  //  Int_t runNumber = eventIn->GetRunNumber();
  //    //AliError(Form("Le numéro de run est le %d",runNumber));
  //
  //  if(!runNumber) return kFALSE;
  //
  //
  //  if(runNumber < 195180 || runNumber > 197469) return kFALSE;
  //
  //    //AliError(Form("est bien après le numéro de run %d",runNumber));
  //
  //    // TString fired = InputEvent() ->GetFiredTriggerClasses();
  //
  //    //AliError(Form("Trigger utilisés dans les evts %s",fired));
  //
  //  if(triggerLevel==1){
  //    threshold = 11.5;
  //    spread = 0.5;
  //  }
  //
  //  if(triggerLevel==2){
  //    threshold = 7.2;
  //    spread = 0.3;
  //  }
  //
  //  if (spread != 0.){
  //    TF1* triggerSmearing =  new TF1("triggerSmearing","[0]*exp(-0.5*((x-[1])/[2])**2)",0,15);
  //    triggerSmearing->SetParameter(0, 1/(spread*TMath::Sqrt(TMath::Pi()*2)));
  //    triggerSmearing->SetParameter(1, threshold);
  //    triggerSmearing->SetParameter(2, spread);
  //    threshold = triggerSmearing->GetRandom();
  //    delete triggerSmearing;
  //  }
  //
  //    //  AliError(Form("passe bien la définition de la fonction du trigger"));
  //  AliTrackContainer *tracks = GetTrackContainer(0);
  //  Int_t localIndex=0;
  //
  //  for(auto it : tracks->accepted()){
  //    AliVTrack* toi = static_cast<AliVTrack*>(it);
  //      //  AliError(Form("Recupère bien les tracks"));
  //
  //    if(!toi) continue;
  //      //    AliError(Form("Recupère bien le tracks"));
  //    if(toi->E() > threshold){
  //        //	AliError(Form("Un tracks passe le critère en énergie"));
  //      return kTRUE;
  //    }
  //  }
  //  return kFALSE;
  //}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::EtIsoClusOnlyEtaBand(AliAODTrack *t, Double_t &ptIso, Double_t &etaBandtrack, Int_t index){
    // Underlying events study with tracks in eta band in EMCAL acceptance
    //Printf("Inside EtaBand");
  
    //INSERT BOUNDARIES ACCORDING TO THE FLAG WE WANT TO USE
  Double_t sumpTConeClust=0., sumpTEtaBandClust=0.;
  Double_t minPhi= 0., maxPhi= 2*TMath::Pi(), minEta = -0.87, maxEta= 0.87;
  AliTrackContainer* tracks = GetTrackContainer(0);

  if(!fTPC4Iso){
    minEta = -0.67;
    maxEta = 0.67;
    minPhi = 1.4;
    maxPhi = 3.15;
  }
  AliClusterContainer *clustAna = GetClusterContainer(0);

    //Printf("tracks ANA : %d" ,tracksAna);
  
  if(!clustAna){
    AliError(Form("Could not retrieve clusters !"));
    return;
  }
  
    //  //Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  Double_t phiClust,etaClust,radius,eTcluster,distCT;
  Int_t nbMObj;
  AliAODCaloCluster *aodClust=0x0;
  TLorentzVector nClust;
  for(auto it : clustAna->accepted()){
    AliAODCaloCluster* coi = static_cast<AliAODCaloCluster*>(it);
    
    coi->GetMomentum(nClust,fVertex);

    phiClust =nClust.Phi();
    etaClust= nClust.Eta();
    eTcluster=0;
    
    Double_t clustTOF = coi->GetTOF()*1e9;
    if(!fIsMC)
      if(clustTOF<-30. || clustTOF>30.) continue;
      //taking off CPV cluster
    Bool_t matched=0;
    nbMObj = coi -> GetNTracksMatched();
    if(nbMObj> 0){
      AliVTrack *mt =0;
      UInt_t rejectionReason = 0;

      for(Int_t i=0;i<nbMObj;i++){
        mt= static_cast<AliVTrack*>(coi->GetTrackMatched(i));
        if (!tracks->AcceptParticle(mt, rejectionReason)) mt = 0;
        
        if(!mt) continue;
        Double_t deltaEta,deltaPhi;
        Double_t veta = mt->GetTrackEtaOnEMCal();
        Double_t vphi = mt->GetTrackPhiOnEMCal();
        deltaEta= veta-etaClust;
        deltaPhi= vphi-phiClust;
        if(!matched && (TMath::Abs(deltaEta)<fdetacut && TMath::Abs(deltaPhi)<fdphicut)){
          matched=kTRUE;
        }
        else
          continue;
      }
    }
    if(matched) continue;
    
    if(nClust.Pt()<0.3) continue;
      //computing distance cluster from *candidate* track.
    radius=TMath::Sqrt(TMath::Power(t->Eta()-etaClust,2)+TMath::Power(t->Phi()-phiClust,2));
    if(radius<fIsoConeRadius && radius!=0.) //In-Cone
      sumpTConeClust += nClust.Pt();
    else
      if(TMath::Abs(etaClust - t->Eta()) < fIsoConeRadius) //Eta band
        sumpTEtaBandClust += nClust.Pt();
  }
  ptIso = sumpTConeClust;
  etaBandtrack = sumpTEtaBandClust;
    //Printf("Returning to FillGeneralHistograms");
}


  //__________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::PtIsoTrackPhiBand(AliAODTrack *t, Double_t &ptIso, Double_t &phiBandtrack, Int_t index){
    // Underlying events study with tracks in phi band in EMCAL acceptance
    //Printf("Inside PhiBand");
  
    //INSERT BOUNDARIES ACCORDING TO THE FLAG WE WANT TO USE
  Double_t sumpTConeCharged=0., sumpTPhiBandTrack=0.;
  Double_t minPhi= 0., maxPhi= 2*TMath::Pi(), minEta = -0.87, maxEta= 0.87;
  
  if(!fTPC4Iso){
    minEta = -0.67;
    maxEta = 0.67;
    minPhi = 1.4;
    maxPhi = 3.15;
  }
  
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  
    //  //Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  tracksAna->ResetCurrentID();
  Int_t localIndex=0;
  
  Double_t phiTrack,etaTrack,radius;
  Int_t iTracksCone = 0;
  AliAODTrack *aodEtrack=0x0;
  while((aodEtrack = static_cast<AliAODTrack*>(tracksAna->GetNextAcceptParticle()))){
    if(!aodEtrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    localIndex=aodEtrack->GetID();
    if(localIndex==index)
      continue;
    
    if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())) continue;
    
    if((aodEtrack->Pt())<0.2)
      continue;
    
    if((aodEtrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    Double_t frac = 0;
    Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
    Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
    if (  ncls> 0 )  frac =  nclsS / ncls ;
    
    if(frac > 0.4) continue;
    
      //CHECK IF TRACK IS IN BOUNDARIES
    phiTrack = aodEtrack->Phi();
    etaTrack = aodEtrack->Eta();
      // define the radius between the leading tracks and the considered tracks
      //redefine phi/t->Eta() from the tracks we passed to the function
    if( (phiTrack < maxPhi) && (phiTrack > minPhi) && (etaTrack < maxEta) && (etaTrack > minEta)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - t->Phi(),2)+TMath::Power(etaTrack - t->Eta(),2));
      
      if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study
        
          // actually phi band here --- ADD Boundaries conditions
        if(TMath::Abs(etaTrack - t->Eta()) < fIsoConeRadius)
          sumpTPhiBandTrack += aodEtrack->Pt();
      }
      else if(radius<fIsoConeRadius && radius!=0.){
        sumpTConeCharged += aodEtrack->Pt(); // should not double count if the track matching is already done
        if(fQA){
          fTestIndex->Fill(index,localIndex);
          fTestLocalIndexPt->Fill(aodEtrack->Pt(),localIndex);
        }
        iTracksCone++;
      }
    }
  }
  fTrackMultvsPt->Fill(iTracksCone,t->Pt());
  ptIso = sumpTConeCharged;
  phiBandtrack = sumpTPhiBandTrack;
}


  //__________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::PtIsoTrackEtaBand(AliAODTrack *t, Double_t &ptIso, Double_t &etaBandtrack, Int_t index){
    // Underlying events study with tracks in eta band in EMCAL acceptance
    //Printf("Inside EtaBand");
  
    //INSERT BOUNDARIES ACCORDING TO THE FLAG WE WANT TO USE
  Double_t sumpTConeCharged=0., sumpTEtaBandTrack=0.;
  Double_t minPhi= 0., maxPhi= 2*TMath::Pi(), minEta = -0.87, maxEta= 0.87;
  
  if(!fTPC4Iso){
    minEta = -0.67;
    maxEta = 0.67;
    minPhi = 1.4;
    maxPhi = 3.15;
  }
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  Int_t localIndex=0;
  
    //Printf("tracks ANA : %d" ,tracksAna);
  
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  
    //  //Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  tracksAna->ResetCurrentID();
  Double_t phiTrack,etaTrack,radius;
  Int_t iTracksCone = 0;
  AliAODTrack *aodEtrack=0x0;
  while((aodEtrack = static_cast<AliAODTrack*>(tracksAna->GetNextAcceptParticle()))){
    if(!aodEtrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    localIndex=aodEtrack->GetID();
    if(localIndex==index)
      continue;
    
      //Printf("\t\thybrid? %s", aodEtrack->IsHybridGlobalConstrainedGlobal()?"yes":"no");
    if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())) continue;
    
    if((aodEtrack->Pt())<0.2)
      continue;
    
      //Printf("\t\tWith ITSrefit? %s",((aodEtrack->GetStatus() & AliAODTrack::kITSrefit)==AliAODTrack::kITSrefit)?"yes":"no");
    if((aodEtrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    Double_t frac = 0;
    Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
    Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
    if (  ncls> 0 )  frac =  nclsS / ncls ;
    
    if(frac > 0.4) continue;
    
    
      //CHECK IF TRACK IS IN BOUNDARIES
    phiTrack = aodEtrack->Phi();
    etaTrack = aodEtrack->Eta();
      //Printf("\t\ttracks with ID: %d, Pt: %.3f, Eta: %.3f, Phi: %.3f",localIndex,aodEtrack->Pt(),etaTrack,phiTrack);
      //Printf("allowed range :Phi in (%f, %f),Eta in (%f, %f)",minPhi,maxPhi,minEta,maxEta);
      // define the radius between the leading tracks and the considered tracks
      //redefine phi/t->Eta() from the tracks we passed to the function
    if( (phiTrack < maxPhi) && (phiTrack > minPhi) && (etaTrack < maxEta) && (etaTrack > minEta)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - t->Phi(),2)+TMath::Power(etaTrack - t->Eta(),2));
        //Printf("\t\tdistance to candidate: %lf with allowed %f",radius,fIsoConeRadius);
      if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study
        
          // actually phi band here --- ADD Boundaries conditions
        if(TMath::Abs(phiTrack - t->Phi()) < fIsoConeRadius)
          sumpTEtaBandTrack += aodEtrack->Pt();
      }
      else if(radius<fIsoConeRadius && radius!=0.){
          //Printf("\t\tInside Cone");
        sumpTConeCharged += aodEtrack->Pt(); // should not double count if the track matching is already done
        if(fQA){
          fTestIndex->Fill(index,localIndex);
          fTestLocalIndexPt->Fill(aodEtrack->Pt(),localIndex);
        }
        iTracksCone++;
      }
    }
  }
    //Printf("filling some histo %f...%f",iTracksCone,t->Pt());
  fTrackMultvsPt->Fill(iTracksCone,t->Pt());
  ptIso = sumpTConeCharged;
  etaBandtrack = sumpTEtaBandTrack;
    //Printf("Returning to FillGeneralHistograms");
}


  //__________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::PtIsoTrackOrthCones(AliAODTrack *t, Double_t &ptIso, Double_t &cones, Int_t index){
    // Underlying events study with tracks in orthogonal cones in TPC
    //Printf("Inside PerpCones");
  
  Double_t sumpTConeCharged=0., sumpTPerpConeTrack=0.;
  
  Double_t etaCandidate = t->Eta();
  Double_t phiCandidate = t->Phi();
  Double_t phiCone1 = phiCandidate - TMath::PiOver2();
  phiCone1=fmod(phiCone1,TMath::TwoPi());
  Double_t phiCone2 = phiCandidate + TMath::PiOver2();
  phiCone2=fmod(phiCone2,TMath::TwoPi());
  
  if (phiCone1 < 0.) phiCone1 += 2*TMath::Pi();
  
  AliTrackContainer *tracks = GetTrackContainer(0);
  Int_t localIndex;
  
  if(!tracks){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  
    //  //Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  tracks->ResetCurrentID();
  AliAODTrack *aodEtrack=0x0;
  Double_t phiTrack,etaTrack,dist2track,dist2Cone1,dist2Cone2;
  Int_t iTracksCone = 0;
  
  while((aodEtrack = static_cast<AliAODTrack*>(tracks->GetNextAcceptParticle()))){
    if(!aodEtrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    localIndex=aodEtrack->GetID();
    if(localIndex==index)
      continue;
    
    if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())) continue;
    
    if((aodEtrack->Pt())<0.2)
      continue;
    
    if((aodEtrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    Double_t frac = 0;
    Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
    Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
    if (  ncls> 0 )  frac =  nclsS / ncls ;
    
    if(frac > 0.4) continue;
    
    phiTrack = aodEtrack->Phi();
    etaTrack = aodEtrack->Eta();
    dist2track = TMath::Sqrt(TMath::Power(etaTrack-etaCandidate, 2)+TMath::Power(phiTrack-phiCandidate, 2));
    
    if (dist2track>fIsoConeRadius){
      dist2Cone1 = TMath::Sqrt(TMath::Power(etaTrack-etaCandidate, 2)+TMath::Power(phiTrack-phiCone1, 2));
      dist2Cone2 = TMath::Sqrt(TMath::Power(etaTrack-etaCandidate, 2)+TMath::Power(phiTrack-phiCone2, 2));
      
        //Is the Track Inside one of the two Cones ->Add to UE
      if((dist2Cone1 < fIsoConeRadius) || (dist2Cone2 < fIsoConeRadius))
        sumpTPerpConeTrack += aodEtrack->Pt();//Distances from the centres of the two Orthogonal Cones
    }
    else if (dist2track<fIsoConeRadius && dist2track!=0){//tracks outside the IsoCone
      sumpTConeCharged += aodEtrack->Pt(); // tracks are inside the isolation cone
      if(fQA){
        fTestIndex->Fill(index,localIndex);
        fTestLocalIndexPt->Fill(aodEtrack->Pt(),localIndex);
      }
      iTracksCone++;
    }
  }
  fTrackMultvsPt->Fill(iTracksCone,t->Pt());
  
  ptIso = sumpTConeCharged;
  cones = sumpTPerpConeTrack;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::PtIsoTrackFullTPC(AliAODTrack *t, Double_t &ptIso, Double_t &full, Int_t index){
    // Underlying events study with tracks in full TPC except back to back bands
    //Printf("Inside FullTPC");
  Double_t sumpTConeCharged=0., sumpTTPCexceptB2B=0.;
  
  AliTrackContainer *tracks = GetTrackContainer(0);
  
  if(!tracks){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  Int_t localIndex;
    //  //Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  tracks->ResetCurrentID();
  Double_t phiTrack,etaTrack,radius, dphiUp, dphiDown;
  Int_t iTracksCone = 0;
  AliAODTrack *aodEtrack=0x0;
  while((aodEtrack = static_cast<AliAODTrack*>(tracks->GetNextAcceptParticle()))){
    if(!aodEtrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    localIndex=aodEtrack->GetID();
    if(localIndex==index)
      continue;
    
    if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())) continue;
    
    if((aodEtrack->Pt())<0.2)
      continue;
    
    if((aodEtrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    Double_t frac = 0;
    Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
    Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
    if (  ncls> 0 )  frac =  nclsS / ncls ;
    
    if(frac > 0.4) continue;
    
    phiTrack = aodEtrack->Phi();
    etaTrack = aodEtrack->Eta();
      //redefine phi/t->Eta() from the tracks we passed to the function
    radius = TMath::Sqrt(TMath::Power(phiTrack-t->Phi(),2)+TMath::Power(etaTrack-t->Eta(),2)); // define the radius between the leading tracks and the considered tracks
    
    if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study UE
      dphiUp = t->Phi() + TMath::Pi() - fIsoConeRadius;
      dphiDown = t->Phi() + TMath::Pi() + fIsoConeRadius;
        // TPC except B2B
      if(phiTrack < dphiDown && phiTrack> dphiUp) sumpTTPCexceptB2B += aodEtrack->Pt();
    }
    else if(radius<fIsoConeRadius && radius!=0){
      sumpTConeCharged += aodEtrack->Pt(); // should not double count if the track matching is already done
      if(fQA){
        fTestIndex->Fill(index,localIndex);
        fTestLocalIndexPt->Fill(aodEtrack->Pt(),localIndex);
      }
      iTracksCone++;
    }
  }
  fTrackMultvsPt->Fill(iTracksCone,t->Pt());
  ptIso = sumpTConeCharged;
  full = sumpTTPCexceptB2B;
}

  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALTrackIsolation::CheckBoundaries(AliAODTrack *t){
    // Check if the cone around the considered tracks is in EMCAL acceptance
    ////Printf("Inside CheckBoundaries\n");
  
  Double_t minPhiBound= 0. , minEtaBound= 0., maxPhiBound= 0., maxEtaBound= 0.;
  Bool_t isINBoundaries;
  
  if(fTPC4Iso) // to avoid to have the isolation partially outside the TPC acceptance in eta
  {
    minPhiBound = 1.4;
    
    if(!fAnalysispPb)
      maxPhiBound = 2.740; // normally 110° but shorter cut to avoid EMCAL border
    else
      maxPhiBound = 3.15; // normally 110° but shorter cut to avoid EMCAL border
    
    minEtaBound = -0.87+fIsoConeRadius; // ""
    maxEtaBound = 0.87-fIsoConeRadius; // ""
  }
  else{
    minEtaBound = -0.67+fIsoConeRadius; // ""
    maxEtaBound = 0.67-fIsoConeRadius; // ""
    
    if(!fAnalysispPb){
      minPhiBound = 1.798;
      maxPhiBound = 2.740; // normally 110° but shorter cut to avoid EMCAL border
    }
    else{
      minPhiBound = 1.4+fIsoConeRadius;
      maxPhiBound = 3.15-fIsoConeRadius; // normally 110° but shorter cut to avoid EMCAL border
    }
  }
  
  if(t->Eta() > maxEtaBound || t->Eta() < minEtaBound || t->Phi() > maxPhiBound || t->Phi() <minPhiBound)
    isINBoundaries=kFALSE;
  else
    isINBoundaries=kTRUE;
  
  return isINBoundaries;
}

  //_________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::LookforParticle(Int_t tracklabel, Double_t ptTRK, Double_t phiTRK, Double_t etaTRK, Double_t isolation){
  
  cout<<"\n\n\n\n\n\n\nInside Look4Particle \n For tracks "<<tracklabel<<"\t\t"<<ptTRK<<"\t\t"<<etaTRK<<"\t\t"<<phiTRK<<"\t\t"<<isolation<<"\n\n\n\n"<<endl;
  
  Int_t trackPDG,momidx;
  Float_t pTTrue , dPhi,dEta, isolationTrue=0., phiTrue, etaTrue;
  
  if (!fIsMC){
    AliWarning("not a montecarlo run!!!!!!");
    return;
  } //AliInfo(Form("It's a MC analysis %e",fAODMCParticles));
  if(!fStack && !fAODMCParticles){
    AliWarning("No Particle Stack !!!!!");
      //cout<<"No Particle Stack !!!!!"<<endl;
    return;
  }
  if(fAODMCParticles->GetEntries() < 1){
    AliWarning("number of tracks insufficient!!!!");
    return;
  }
  
  Int_t ndimsMCmix = fMCQAdim;
  Double_t outputvalueMCmix[ndimsMCmix];
    //    //cout<<"dimensions of the array: "<<ndimsMCmix<<endl;
  
  Int_t npart=fAODMCParticles->GetEntries();
    //cout<<"Number of particles in the event: "<<npart<<endl;
  
  AliAODMCParticle *particle2Check, *momP2Check;
    //    //DEVELOP THE CODE TO LOOK FOR TRACKS IN MC
    //  if(tracklabel<0)
    //    return;
  
  particle2Check=static_cast<AliAODMCParticle*>(fAODMCParticles->At(TMath::Abs(tracklabel)));
  
  trackPDG=particle2Check->GetPdgCode();
  etaTrue=particle2Check->Eta();
  phiTrue=particle2Check->Phi();
  pTTrue=particle2Check->Pt();
  
  dPhi = phiTRK - phiTrue;
  dEta = etaTRK - etaTrue;
  
  momidx=particle2Check->GetMother();
  if(momidx<0 || momidx>npart)
    return;
  
  momP2Check=static_cast<AliAODMCParticle*>(fAODMCParticles->At(momidx));
  AliAODMCParticle *part;
  
  for(int ipart=0;ipart<npart;ipart++){
    if(ipart==TMath::Abs(tracklabel))
      continue;
    
    part=static_cast<AliAODMCParticle*>(fAODMCParticles->At(ipart));
    
    if(part->GetStatus()>10 || (!part->IsPrimary()) || (!part->IsPhysicalPrimary()))
      continue;
    
    if(fIsoMethod==0 && part->Charge()==0)
      continue;
    if(fIsoMethod==1 && part->Charge()!=0)
      continue;

    Float_t phiPart = part->Phi();
    Float_t etaPart = part->Eta();
    
    Float_t dist2Candidate= TMath::Sqrt(((phiPart-phiTrue)*(phiPart-phiTrue))+((etaPart-etaTrue)*(etaPart-etaTrue)));
    
    if((dist2Candidate < fIsoConeRadius) && dist2Candidate!=0)
      isolationTrue += part->Pt();
    else
      continue;
  }
    //
  outputvalueMCmix[0] = ptTRK;
  outputvalueMCmix[1] = trackPDG;
  outputvalueMCmix[2] = momP2Check->GetPdgCode();
  outputvalueMCmix[3] = pTTrue;
  outputvalueMCmix[4] = dPhi;
  outputvalueMCmix[5] = dEta;
  outputvalueMCmix[6] = isolation;
  outputvalueMCmix[7] = isolationTrue;
  
  fOutTrackMC->Fill(outputvalueMCmix);
  
  return;
  
  
}
  //__________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::IsolationAndUEinEMCAL(AliAODTrack *toi, Double_t& isolation,Double_t& ue,Double_t pTThreshold, Int_t index){
  
    //Printf("Inside IsolationAndUEinEMCAL");
    //EMCAL (Only tracks for the Isolation since Isolated Tracks analysis)
  Double_t isoConeArea = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  
  Double_t etaBandAreaTr = 1.74*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandAreaTr = 2.*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  
  toi->GetP();
  Double_t pTTOI=toi->Pt();
  
  switch (fUEMethod)
  {
    case 0: //phi band
      PtIsoTrackPhiBand(toi, isolation, ue, index);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / phiBandAreaTr);
      fPhiBandUETracks->Fill(toi->Pt() , ue);
      ue = ue * (isoConeArea / phiBandAreaTr);
        // fill histograms for isolation
        //      if(fWho==2) fPtvsM02vsSum->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction
      
        //      if(fWho==2) {
        //        fPtvsM02vsSumUE->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      fPtIsoTrack->Fill(toi->Pt(), isolation);
        //      }
      
      if(isolation<pTThreshold)
      {
          //        FillInvMassHistograms(kTRUE, vecTOI, index,isolation);
          //        if(fWho==2) {
          //          fPtvsM02iso->Fill(vecTOI.Pt(),toi->GetM02());
        fPtIsolatedNTracks->Fill(toi->Pt());
          //        }
          //        fPtisoT=vecTOI.Pt();
          //        fM02isoT=m02TOI;
        
          //        if(fM02mincut < m02TOI && m02TOI < fM02maxcut)
          //        {
          //          if(fWho==2) fEtIsolatedTracks->Fill(eTTOI);
          //          fEtisolatedT=eTTOI;
          //          fPtisolatedT=vecTOI.Pt();
          //        }
      }
        //      else
        //      {
        //          //        if(isolation>3.)   FillInvMassHistograms(kFALSE, m02TOI, vecTOI, index,isolation);
        //          //        if(fWho==2) fPtvsM02noiso->Fill(vecTOI.Pt(),toi->GetM02());
        //        fPtNotIsolatedNtracks=vecTOI.Pt();
        //          //        fM02noisoT=m02TOI;
        //      }
      break;
    case 1: //eta band
      if(fIsoMethod==1)
        EtIsoClusOnlyEtaBand(toi, isolation, ue, index);
      else
        PtIsoTrackEtaBand(toi, isolation, ue, index);
      
      fEtaBandUETracks->Fill(toi->Pt() , ue);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / etaBandAreaTr);
        //      if(fWho==2)
      ue = ue * (isoConeArea / etaBandAreaTr);
        // fill histograms for isolation
        //      if(fWho==2) fPtvsM02vsSum->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction
      
        //      if(fWho==2) {
        //        fPtvsM02vsSumUE->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      fPtIsoTrack->Fill(toi->Pt(), isolation);
        //      }
      if(isolation<pTThreshold)
      {
          //        FillInvMassHistograms(kTRUE, vecTOI, index,isolation);
          //        if(fWho==2) {
          //          fPtvsM02iso->Fill(vecTOI.Pt(),toi->GetM02());
        fPtIsolatedNTracks->Fill(toi->Pt());
          //        }
          //        fPtisoT=vecTOI.Pt();
          //        fM02isoT=m02TOI;
        
          //        if(fM02mincut < m02TOI && m02TOI < fM02maxcut)
          //        {
          //          if(fWho==2) fEtIsolatedTracks->Fill(eTTOI);
          //          fEtisolatedT=eTTOI;
          //          fPtisolatedT=vecTOI.Pt();
          //        }
      }
        //      else
        //      {
        //          //        if(isolation>3.)   FillInvMassHistograms(kFALSE, m02TOI, vecTOI, index,isolation);
        //          //        if(fWho==2) fPtvsM02noiso->Fill(vecTOI.Pt(),toi->GetM02());
        //        fPtnoisoT=vecTOI.Pt();
        //          //        fM02noisoT=m02TOI;
        //      }
      break;
  }
  
}



  //__________________________________________________________________________
void AliAnalysisTaskEMCALTrackIsolation::IsolationAndUEinTPC(AliAODTrack *toi, Double_t& isolation,Double_t& ue,Double_t pTThreshold, Int_t index){
  
    //EMCAL + TPC (Only tracks for the Isolation since IsoCone Goes Out of EMCAL)
  
  Double_t isoConeArea = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  
  Double_t etaBandAreaTr = 1.74*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandAreaTr = 2.*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  Double_t perpConesArea = 2.*isoConeArea;
  Double_t fullTPCArea = (1.8*2.*TMath::Pi())-(fIsoConeRadius*2*1.8)-isoConeArea;
  
  toi->GetP();
  Double_t pTTOI=toi->Pt();
  
  switch (fUEMethod)
  {
    case 0: //phi band
      PtIsoTrackPhiBand(toi, isolation, ue, index);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / phiBandAreaTr);
      fPhiBandUETracks->Fill(toi->Pt() , ue);
      ue = ue * (isoConeArea / phiBandAreaTr);
        // fill histograms for isolation
        //      if(fWho==2) fPtvsM02vsSum->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction
      
        //      if(fWho==2) {
        //        fPtvsM02vsSumUE->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      fPtIsoTrack->Fill(toi->Pt(), isolation);
        //      }
      
      if(isolation<pTThreshold)
      {
          //        FillInvMassHistograms(kTRUE, vecTOI, index,isolation);
          //        if(fWho==2) {
          //          fPtvsM02iso->Fill(vecTOI.Pt(),toi->GetM02());
        fPtIsolatedNTracks->Fill(toi->Pt());
          //        }
          //        fPtisoT=vecTOI.Pt();
          //        fM02isoT=m02TOI;
        
          //        if(fM02mincut < m02TOI && m02TOI < fM02maxcut)
          //        {
          //          if(fWho==2) fEtIsolatedTracks->Fill(eTTOI);
          //          fEtisolatedT=eTTOI;
          //          fPtisolatedT=vecTOI.Pt();
          //        }
      }
        //      else
        //      {
        //          //        if(isolation>3.)   FillInvMassHistograms(kFALSE, m02TOI, vecTOI, index,isolation);
        //          //        if(fWho==2) fPtvsM02noiso->Fill(vecTOI.Pt(),toi->GetM02());
        //        fPtNotIsolatedNtracks=vecTOI.Pt();
        //          //        fM02noisoT=m02TOI;
        //      }
      break;
    case 1: //eta band
      PtIsoTrackEtaBand(toi, isolation, ue, index);
      fEtaBandUETracks->Fill(toi->Pt() , ue);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / etaBandAreaTr);
        //      if(fWho==2)
      ue = ue * (isoConeArea / etaBandAreaTr);
        // fill histograms for isolation
        //      if(fWho==2) fPtvsM02vsSum->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction
      
        //      if(fWho==2) {
        //        fPtvsM02vsSumUE->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      fPtIsoTrack->Fill(toi->Pt(), isolation);
        //      }
      if(isolation<pTThreshold)
      {
          //        FillInvMassHistograms(kTRUE, vecTOI, index,isolation);
          //        if(fWho==2) {
          //          fPtvsM02iso->Fill(vecTOI.Pt(),toi->GetM02());
        fPtIsolatedNTracks->Fill(toi->Pt());
          //        }
          //        fPtisoT=vecTOI.Pt();
          //        fM02isoT=m02TOI;
        
          //        if(fM02mincut < m02TOI && m02TOI < fM02maxcut)
          //        {
          //          if(fWho==2) fEtIsolatedTracks->Fill(eTTOI);
          //          fEtisolatedT=eTTOI;
          //          fPtisolatedT=vecTOI.Pt();
          //        }
      }
        //      else
        //      {
        //          //        if(isolation>3.)   FillInvMassHistograms(kFALSE, m02TOI, vecTOI, index,isolation);
        //          //        if(fWho==2) fPtvsM02noiso->Fill(vecTOI.Pt(),toi->GetM02());
        //        fPtnoisoT=vecTOI.Pt();
        //          //        fM02noisoT=m02TOI;
        //      }
      break;
    case 2: //Cones
      PtIsoTrackOrthCones(toi, isolation, ue, index);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / perpConesArea);
        //      if(fWho==2)
      fPerpConesUETracks ->Fill(toi->Pt() , ue);
      ue = ue * (isoConeArea / perpConesArea);
        // fill histograms for isolation
        //      if(fWho==2) fPtvsM02vsSum->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction
      
        //      if(fWho==2) {
        //        fPtvsM02vsSumUE->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      fPtIsoTrack->Fill(toi->Pt(), isolation);
        //      }
      if(isolation<pTThreshold)
      {
          //        FillInvMassHistograms(kTRUE, m02TOI, vecTOI, index, isolation);
          //        if(fWho==2) {
          //          fPtvsM02iso->Fill(vecTOI.Pt(),toi->GetM02());
        fPtIsolatedNTracks->Fill(toi->Pt());
          //        }
          //        fPtisoT=vecTOI.Pt();
          //        fM02isoT=m02TOI;
        
          //        if(fM02mincut < m02TOI && m02TOI < fM02maxcut)
          //        {
          //          if(fWho==2) fEtIsolatedTracks->Fill(eTTOI);
          //          fEtisolatedT=eTTOI;
          //          fPtisolatedT=vecTOI.Pt();
          //        }
      }
        //      else
        //      {
        //          //        if(isolation>3.)  FillInvMassHistograms(kFALSE, m02TOI, vecTOI, index, isolation);
        //          //        if(fWho==2) fPtvsM02noiso->Fill(vecTOI.Pt(),toi->GetM02());
        //          //        fPtnoisoT=vecTOI.Pt();
        //          //        fM02noisoT=m02TOI;
        //      }
      break;
    case 3: //Full TPC
      PtIsoTrackFullTPC(toi, isolation, ue, index);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / fullTPCArea);
        //      if(fWho==2)
      fTPCWithoutIsoConeB2BbandUE->Fill(toi->Pt() , ue);
      ue = ue * (isoConeArea / fullTPCArea);
        // fill histograms for isolation
      
        //      if(fWho==2) fPtvsM02vsSum->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction
      
        //      if(fWho==2) {
        //        fPtvsM02vsSumUE->Fill(vecTOI.Pt(),toi->GetM02(),isolation);
      fPtIsoTrack->Fill(toi->Pt(), isolation);
        //      } //	fTracksConeEtaPt->Fill(isolation, vecTOI.Eta(), vecTOI.Pt());
        //        //	fTracksConeEtaM02->Fill(isolation, vecTOI.Eta(), toi->GetM02());
      if(isolation<pTThreshold)
      {
          //        FillInvMassHistograms(kTRUE, m02TOI, vecTOI, index, isolation);
          //        if(fWho==2) {
          //          fPtvsM02iso->Fill(vecTOI.Pt(),toi->GetM02());
        fPtIsolatedNTracks->Fill(toi->Pt());
          //        }
          //        fPtisoT=vecTOI.Pt();
          //        fM02isoT=m02TOI;
          //
          //        if(fM02mincut < m02TOI && m02TOI < fM02maxcut)
          //        {
          //          if(fWho==2) fEtIsolatedTracks->Fill(eTTOI);
          //          fEtisolatedT=eTTOI;
          //          fPtisolatedT=vecTOI.Pt();
          //        }
      }
        //      else
        //      {
        //        if(isolation>3.) FillInvMassHistograms(kFALSE, m02TOI, vecTOI, index, isolation);
        //        if(fWho==2) fPtvsM02noiso->Fill(vecTOI.Pt(),toi->GetM02());
        //        fPtnoisoT=vecTOI.Pt();
        //        fM02noisoT=m02TOI;
        //      }
      break;
  }
  
}


  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALTrackIsolation::FillGeneralHistograms(AliAODTrack *toi, Int_t index){
    //AliInfo("Inside FillGeneralHistograms\n");
    //Printf("Inside FillGeneralHistograms\n");
    // Fill the histograms for underlying events and isolation studies
    // AliError(Form("Arrive bien dans fill general histograms"));
  
    // I would like to remove this part and fill the tracks multiplicity histogram in FillQAHistograms, is that ok for thnSparses? (especially cause here the histogram is filled several times per event)
    // AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
    //Printf("Name of the tracks used for Isolation: %s",(tracksAna->GetClassName()).Data());
  
  tracksAna->ResetCurrentID();
  
  const Int_t nTracks = tracksAna->GetNAcceptedTracks();
    //Printf("Ntracks for the event with this tracks: %d", nTracks);
  if(fQA)
    fTrackMult->Fill(nTracks);
  
  Double_t pTTOI = 0.;
  
    //Definition of the Array for Davide's Output
  const Int_t ndims =   fNDimensions;
  Double_t outputValues[ndims];
  
  pTTOI = toi->Pt();
  
    //Printf("tracks ID %d ",toi->GetID());
  
    // ******** Isolation and UE calculation with different methods *********
  
  Double_t pTThreshold = 2.;
  
  switch(fptIsoMethod)
  {
    case 0:  // SumEt<EtThr
      pTThreshold = fptIsoThreshold;
      break;
      
    case 1:  // SumEt<%Ephoton
      pTThreshold = fptIsoThreshold * pTTOI;
      break;
      
    case 2: // Etmax<pTThreshold
      pTThreshold = fptIsoThreshold;
  }
  
  Double_t isolation=0, ue=0;
  
  if(!fTPC4Iso)
    IsolationAndUEinEMCAL(toi,isolation,ue,pTThreshold,index);
  else
    IsolationAndUEinTPC(toi,isolation,ue,pTThreshold,index);
  
    //Printf("after Isolation and UE in %s",fTPC4Iso?"TPC":"EMCAL");
  if(fIsMC)
    LookforParticle(TMath::Abs(toi->GetLabel()),toi->Pt(),toi->Phi(),toi->Eta(),isolation);
  
    //Printf("after Look4Particle");
  /*  Here we should call something to know the number of tracks...
   Soon I'll put in this version the "old way"; please let me know if
   any of you could do the same with the JET framework*/
  
  
  outputValues[0] = pTTOI;
  outputValues[1] = isolation;
  outputValues[2] = ue;
  outputValues[3] = toi->Eta();
  outputValues[4] = toi->Phi();
  fOutputTHnS -> Fill(outputValues);
  
  return kTRUE;
}

  //_________________________________________________________________________

void AliAnalysisTaskEMCALTrackIsolation::AddParticleToUEMC(Double_t& sumUE,AliAODMCParticle* mcpp, Double_t eta,Double_t phi){
  
  Double_t etap=mcpp->Eta();
  Double_t phip=mcpp->Phi();
  
  if(!fTPC4Iso){
    if(TMath::Abs(etap)>=0.7 || (phip<=1.4 || phip>= TMath::Pi()))
      return;
    else{
      switch(fUEMethod){
        case 0: //Phi band
          if(TMath::Abs(eta-etap)<fIsoConeRadius) //to be changed with fIsoConeRadius
            sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
          else
            return;
          
          break;
        case 1: //Eta band
          if(TMath::Abs(phi-phip)<fIsoConeRadius)
            sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
          else
            return;
          
          break;
      }
    }
  }
  else{
    if(TMath::Abs(etap)>=1.0)
      return;
    else{
      switch(fUEMethod){
        case 0: //Phi band
        {if(TMath::Abs(eta-etap)<fIsoConeRadius) //to be changed with fIsoConeRadius
          sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
        else
          return;
          break;
        }
        case 1: //Eta band
        {  if(TMath::Abs(phi-phip)<fIsoConeRadius)
          sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
        else
          return;
          
          break;
        }
        case 2: //Orthogonal Cones
        { double etacone1= eta;
          double etacone2= eta;
          double phicone1= phi - TMath::PiOver2();
          double phicone2= phi + TMath::PiOver2();
          
          if (phicone1 < 0.) phicone1 += 2*TMath::Pi();
          
          if(TMath::Sqrt(TMath::Power(etap-etacone1,2)+TMath::Power(phip-phicone1,2))< fIsoConeRadius ||
             TMath::Sqrt(TMath::Power(etap-etacone2,2)+TMath::Power(phip-phicone2,2))< fIsoConeRadius) //to be changed with fIsoConeRadius
          {sumUE += mcpp->Pt();}
          else
            return;
          
          break;
        }
        case 3: //Full TPC
        {    //                  Double_t phiup= phi +TMath::Pi()+fIsoConeRadius;
             //                  Double_t phidown= phi +TMath::Pi()-fIsoConeRadius;
             //
             //                  if(phip < phidown || phip > phiup ) //TO BE CHECKED
             //                    continue;
          break;
        }
      }
    }
  }
}
  //_________________________________________________________________________

void AliAnalysisTaskEMCALTrackIsolation::CalculateUEDensityMC(Double_t& sumUE){
  
  Double_t isoConeArea = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  Double_t etaBandAreaTr,phiBandAreaTr,perpConesArea,fullTPCArea;
  if(fIsoMethod==0){
    etaBandAreaTr = (1.8 - 2.*fIsoConeRadius)*2.*fIsoConeRadius;
    phiBandAreaTr = (2.*TMath::Pi() - 4.*fIsoConeRadius)*2.*fIsoConeRadius;//there is a 2 more because of the JET CONE B2B
    perpConesArea = 2.*isoConeArea;
    fullTPCArea = 1.8*2.*TMath::Pi()-fIsoConeRadius*(TMath::Pi()*fIsoConeRadius + 3.6);
  }
  else{
    etaBandAreaTr = (1.4 - 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  }
  switch (fUEMethod){
    case 0:
      sumUE = sumUE * (isoConeArea / phiBandAreaTr);
      break;
    case 1:
      sumUE = sumUE * (isoConeArea / etaBandAreaTr);
      break;
    case 2:
      sumUE = sumUE * (isoConeArea / perpConesArea);
      break;
    case 3:
      sumUE = sumUE * (isoConeArea / fullTPCArea);
      break;
  }
}

  //_________________________________________________________________________

void AliAnalysisTaskEMCALTrackIsolation::AnalyzeMC(){
  
    //printf("New Event...Analysing Stack!");
  if (!fIsMC)
    return;
    //AliInfo(Form("It's a MC analysis %e",fAODMCParticles));
  if(!fStack && !fAODMCParticles)
  {cout<<"no stack saved\n"; return;}
  
    //cout<<"there's a List of particles"<<endl;
    //DO THIS ALSO FOR ESDs
  
  Double_t pT, sumpTiso, sumUE,phi, eta, distance, phip, etap, mcfirstPt;
  
  if(fAODMCParticles->GetEntries() < 1){
    AliError("number of tracks insufficient");
    return;
  }
  int nDimMC = fMCDimensions;
  Double_t outputValuesMC[nDimMC];
  
  Int_t nTracks = fAODMCParticles->GetEntriesFast();
  Int_t nFSParticles = 0;
  AliAODMCParticle *multTracks;
  
  for(int a=0; a<nTracks; a++){
    
    multTracks = static_cast<AliAODMCParticle*>(fAODMCParticles->At(a));
    
    if(multTracks->IsPrimary() && multTracks->IsPhysicalPrimary() && multTracks->GetStatus()<10){
      if(TMath::Abs(multTracks->Eta())<=0.9 && multTracks->Charge()!= 0)
        nFSParticles++;
      else
        continue;
    }//implement final state particle condition
    else
      continue;
  }
    //AliInfo(Form("number of particles in the array %d",nTracks));
  AliAODMCParticle *mcpart, *mom, *mcpp,*mcsearch, *mcfirst, *mcfirstmom,*matchingtrack, *mum;
  
    // Bool_t prompt=kFALSE;
  Double_t mcpT, maxpT;
  Int_t pdg, mompdg, tracklabel;
  Double_t mcFirstEta=0., mcFirstPhi=0.;
  Double_t charge;
  
    // AliAODMCParticle *mcfirst = static_cast<AliAODMCParticle*>(fAODMCParticles->At(0));
    //AliAODMCParticle *mcp, *mcpmaxE, *mcpp, *mom;
    //  if(!fisLCAnalysis){
    //Loop on the event
  if(fisLTAnalysis){
    maxpT=0.;
    int indexmaxpT=0;
      //getting the index of the particle with the maximum energy.
    for(int iTr=0;iTr<nTracks;iTr++){
      mcsearch= static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTr));
      
      if(!mcsearch) continue;
      
      if(mcsearch->GetStatus()>10) continue;
      if(mcsearch->Charge()==0) continue;
      if(!mcsearch->IsPrimary()) continue;
      
      if(TMath::Abs(mcsearch->Eta())>0.67-fIsoConeRadius) continue;
      if(!fTPC4Iso){
        if((TMath::Abs(mcsearch->Eta())>0.67-fIsoConeRadius ) || (mcsearch->Phi() < 1.398 + fIsoConeRadius || mcsearch->Phi()>(TMath::Pi()-fIsoConeRadius-0.03)))
          continue;
      }
      else {
        if((TMath::Abs(mcsearch->Eta())>0.87-fIsoConeRadius ) || (mcsearch->Phi() < 1.398 || mcsearch->Phi()>(TMath::Pi()-0.03)))
          continue;
      }
      
      mcfirstPt= mcsearch->Pt();
      if(mcfirstPt>maxpT){
        maxpT=mcfirstPt;
        indexmaxpT=iTr;
      }
      else continue;
    }
    mcfirst= static_cast<AliAODMCParticle*>(fAODMCParticles->At(indexmaxpT));
    mcfirstPt=mcfirst->Pt();
    
    int momidx= mcfirst->GetMother();
    if(momidx>0){
      mcfirstmom =  static_cast<AliAODMCParticle*>(fAODMCParticles->At(momidx));
      mompdg= TMath::Abs(mcfirstmom->GetPdgCode());
    }
    else
      mompdg=mcfirst->GetPdgCode();
    
    mcFirstEta = mcfirst->Eta();
    mcFirstPhi = mcfirst->Phi();
    
    phip=0., etap=0.;
    sumpTiso=0,sumUE=0;
    
    for(Int_t iTrack=1;iTrack<nTracks ;iTrack++){
      if(iTrack==indexmaxpT) continue;
      
      mcpp= static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
      if(!mcpp)
        continue;
      
      if(fIsoMethod==0 && mcpp->Charge()==0) continue;
      if(fIsoMethod==1 && mcpp->Charge()!=0) continue;
      
      phip = mcpp->Phi();
      etap = mcpp->Eta();
      
      if(mcpp->GetStatus()>10) continue;
      if(fIsoMethod==0 && !mcpp->IsPrimary())continue;
      if(fIsoMethod==1 && !mcpp->IsPhysicalPrimary()) continue;
      
      distance=0.;
      distance= TMath::Sqrt((mcFirstPhi- phip)*(mcFirstPhi- phip) + (mcFirstEta- etap)*(mcFirstEta- etap));
      
      if(distance<=fIsoConeRadius){
        sumpTiso += mcpp->Pt();
      }
      else{
        AddParticleToUEMC(sumUE,mcpp,mcFirstEta,mcFirstPhi);
      }
    }
      //  cout<<"\n\nTotal Energy inside the Isolation Cone : "<<sumEiso<<endl;
    CalculateUEDensityMC(sumUE);
      //cout<<"Total UE Energy : "<<sumUE<<" calculated with method "<<fUEMethod<<endl;
    outputValuesMC[0] = mcfirstPt;
    outputValuesMC[1] = sumpTiso;
    outputValuesMC[2] = sumUE;
    outputValuesMC[3] = mompdg;
    outputValuesMC[4] = mcFirstEta;
    outputValuesMC[5] = mcFirstPhi;
    outputValuesMC[6] = mcfirst->GetLabel();
      // EtaPhiMCPhoton
      // EtMC
      // EtIsoCone
      // EtMother
      // UE Et
      // Mother PDG
      //fill some histograms or a THnSparse or a TTree.
      //	AliError(Form("Fill something in Analize MC"));
    fOutMCTruth -> Fill(outputValuesMC);

      //Fill the Output TTree for MC Truth
  }
  else{
    
    for(int iTr=0;iTr<nTracks;iTr++){
      
      mcpT=0.;pT =0; phi=0.; eta=0.;
      
      mcpart = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTr));
      
      if(mcpart->GetStatus()>10) {continue;}
      if(!mcpart->IsPrimary()) {continue;}
      if(!mcpart->IsPhysicalPrimary()) {continue;}
      
      charge = mcpart->Charge();
      if(charge==0)
        continue;
      
      eta = mcpart->Eta();
      phi = mcpart->Phi();
      
      if(!fTPC4Iso){
        if((TMath::Abs(eta)>0.67-fIsoConeRadius ) || (phi < 1.398 + fIsoConeRadius || phi>(TMath::Pi()-fIsoConeRadius-0.03)))
          continue;
      }
      else{
        if((TMath::Abs(eta)>0.87-fIsoConeRadius ) || (phi < 1.398 || phi>(TMath::Pi()-0.03)))
          continue;
      }
        //printf("\nParticle Position %d  and Label: %d  PDG: %d  Pt: %f  Eta: %f  Phi: %f",iTr, mcpart->GetLabel(),pdg,mcpart->Pt(), eta, phi);
      
      tracklabel = iTr;
      int momidx = mcpart->GetMother();
      
      if(momidx>0){
        mom = static_cast<AliAODMCParticle*>(fAODMCParticles->At(momidx));
        mompdg= TMath::Abs(mom->GetPdgCode());
      }
      else
        mompdg=mcpart->GetPdgCode();
      
        //printf("With Mother at %d with label %d which is a %d",momidx, mom->GetLabel(), mompdg);
      
      pT= mcpart->Pt(); //transform to transverse Energy
      
      bool foundmatch=kFALSE;
        //This loop excludes tracks which share the same direction with other particles in their isolation cone
      for(int m=0;m<nTracks && foundmatch==kFALSE;m++){
          //not the same track
        if(m==iTr) continue;
        
        matchingtrack = static_cast<AliAODMCParticle*>(fAODMCParticles->At(m));
        
        if(! matchingtrack->IsPrimary()) continue;
        if(! matchingtrack->IsPhysicalPrimary()) continue;
        if(matchingtrack->GetStatus()> 10 ) continue;
        
        Double_t etamatching = matchingtrack->Eta();
        Double_t phimatching = matchingtrack->Phi();
        
        if(TMath::Abs(eta-etamatching)<=fdetacut && TMath::Abs(phi-phimatching)<=fdphicut)
          foundmatch=kTRUE;
      }
      if(foundmatch) continue;
      
      distance=0.;
      phip=0., etap=0.;
      sumpTiso=0.,sumUE=0.;
      
      for(int iTrack=0;iTrack<nTracks;iTrack++){
        
        if(iTrack==tracklabel)
          continue;
        
        mcpp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
        
        if(!mcpp) {continue;}
        
        if(fIsoMethod==0 && mcpp->Charge()==0) continue;
        if(fIsoMethod==1 && mcpp->Charge()!=0) continue;
        
        if(mcpp->GetStatus()<10)
          fPtTracksVSpTTR_MC->Fill(pT,mcpp->Pt());
        else
          continue;
        
        int mumidx=mcpp->GetMother();
        if (mumidx<0 || mumidx>nTracks) continue;
        
        mum = static_cast<AliAODMCParticle*>(fAODMCParticles->At(mumidx));
        if(mumidx == tracklabel || mum->GetPdgCode()==22) continue;
        
        phip = mcpp->Phi();
        etap = mcpp->Eta();
        
          //Depending on which Isolation method and UE method is considered.
        distance= TMath::Sqrt((phi-phip)*(phi-phip) + (eta-etap)*(eta-etap));
        
        if(distance <= fIsoConeRadius){
            //cout<<iTrack<<"\t"<<photonlabel<<endl;
            //mcpp->Print();
          sumpTiso += mcpp->Pt();
        }
        else{
          AddParticleToUEMC(sumUE,mcpp, eta, phi);}
      }
      CalculateUEDensityMC(sumUE);
      
        //printf("Storing Particle: Label %d  PDG: %d  Eta: %f  Phi: %f",mcpart->GetLabel(),pdg,eta,phi);
        //printf("With Mother at %d with label %d which is a %d",momidx, mom->GetLabel(), mompdg);
      outputValuesMC[0] = pT;
      outputValuesMC[1] = sumpTiso;
      outputValuesMC[2] = sumUE;
      outputValuesMC[3] = mompdg;
      outputValuesMC[4] = eta;
      outputValuesMC[5] = phi;
      outputValuesMC[6] = mcpart->GetLabel();
        // EtaPhiMCPhoton
        // EtMC
        // EtIsoCone
        // EtMother
        // UE Et
        // Mother PDG
        //fill some histograms or a THnSparse or a TTree.
        //	AliError(Form("Fill something in Analize MC"));
      fOutMCTruth -> Fill(outputValuesMC);
    }
  }
  
  return;
}

