#include "TList.h"
#include "TTree.h"
#include "TH3F.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliTRDTriggerAnalysis.h"
#include "AliAnalysisTaskFragmentationTriggered.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAODHeader.h"
#include "AliMultiplicity.h"
#include <iostream>
#include <string>
#include <vector>
#include "TVector2.h"
#include "TFormula.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"

//MC stuff
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliMCEvent.h"
#include "AliStack.h"
//#include "AliAODMCHeader.h"
//#include "AliAODMCParticle.h"

//random generator
#include "TRandom3.h"


using std::cout;
using std::endl;

AliAnalysisTaskFragmentationTriggered::AliAnalysisTaskFragmentationTriggered(const char *name) :
  //general objects
  AliAnalysisTaskSE(name),
  fLeading(0),
  kmbtrigger(kTRUE),
  khjttrigger(kFALSE),
  //MC Objects
  fMCRC(0x0), //random cone
  fMCjetstracks(0x0),
  fMCjets(0x0),
  fMCtracks_pcp(0x0),
  fMCtracks_pcm(0x0),
  //reco objects
  fRC(0x0),
  fJets(0x0),
  fJetsTracks(0x0),
  fTracks_pcp(0x0),
  fTracks_pcm(0x0),
  fTracks_pcp_LJCM(0x0),
  fTracks_LJCM(0x0),
  fCuts(0x0),
  fCuts_noSPD(0x0),
  fCuts_noITS(0x0),
  fAnaUtils(0x0),
  fMCTruthParticles(0x0),
  fRandomGen(0x0),
  fOutputList(0x0),   
  //MC Truth histos:
  fHistMCTruthPtHard(0x0),
  fHistMCTruthPtHardTrials(0x0),
  fHistMCTruthNTracks(0x0),
  fHistMCTruthTrackPt(0x0),
  fHistMCTruthTrackPhi(0x0),
  fHistMCTruthTrackEta(0x0),
  fMCTruthRecTrackPtGen(0x0),
  fMCTruthRecTrackSecPtGen(0x0),
  fMCTruthRecTrackSecPtRec(0x0),
  fHistMCTruthTrackPtRatio(0x0),
  fHistMCTruthTrackPtRatiowSPD(0x0),
  fHistMCTruthTrackPtRationoSPD(0x0),

  
  //MC Truth jets:
  fMCTruthnJets(0x0),
  fMCTruthAntiktHistPt(0x0),
  fMCTruthAntiktHistAreavsPt(0x0),
  fMCTruthAntiktHistAreaoverCirclevsPt(0x0),
  fMCTruthAntiktHistxipt(0x0),
  fMCTruthAntiktHistxipt_pcp(0x0),
  fMCTruthAntiktHistxipt_pcm(0x0),
  fMCTruthAntiktHistxipt_pc(0x0),

  //leading jets
  fMCTruthAntiktHistPtLeading(0x0),
  fMCTruthAntiktHistEtaLeading(0x0),
  fMCTruthAntiktHistPhiLeading(0x0),
  fMCTruthAntiktHistxiptLeading(0x0),
  fMCTruthAntiktHistxiptLeading_pcp(0x0),

  //Random cone histos
  fMCTruthRCntry(0x0),
  fMCTruthAntiktHistPhiRC(0x0),
  fMCTruthAntiktHistEtaRC(0x0),
  fMCTruthRCnTracksOverlap(0x0),
  fMCTruthHistRCPt(0x0),
  fMCTruthHistRCnTracks(0x0),
  fMCTruthHistRCWithOverlapPt(0x0),
  fMCTruthHistRCWithOverlapnTracks(0x0),
  fMCTruthHistRCTrackEta(0x0),
  fMCTruthHistRCTrackPhi(0x0),
  fMCTruthHistRCTrackEtaDiscarded(0x0),
  fMCTruthHistRCTrackPhiDiscarded(0x0),
  fMCTruthHistRCxipt(0x0),
  fMCTruthHistRCxiptLeading(0x0),
  fMCTruthHistRCxiptLeadingDiscarded(0x0),
 
  //rho histos
  fMCTruthHistRho(0x0),
  fMCTruthHistRhovsLJPt(0x0),
  fMCTruthHistRhoSubPt(0x0),
  fMCTruthHistRhoSubAreaOverCirclevsPt(0x0),
  fMCTruthHistRhoSubPhiTracksJetvsPt(0x0),
  fMCTruthHistRhoSubEtaTracksJetvsPt(0x0),
  fMCTruthHistRhoSubPhiTracksPCvsPt(0x0),
  fMCTruthHistRhoSubEtaTracksPCvsPt(0x0),
  fMCTruthHistRhoSubXiPt(0x0),
  fMCTruthHistRhoSubXiPCPt(0x0),

  //  track histos:
  fHistSumw2Debug(0x0),
  fHistEventDebug(0x0),
  fHistCuts(0x0),
  fHistTrackCuts(0x0),
  fHistTrackCutsDep(0x0),
  fHistPt(0x0),
  fHistnTracks(0x0),
  fHistEta(0x0),
  fHistPhi(0x0),
  fHistEtaPhi(0x0),
  fHistPhiHybridUnconstrained(0x0),
  fHistPhiHybridConstrainedwITS(0x0),
  fHistPhiHybridConstrainedwoITS(0x0),
  fHistITSNcls(0x0),
  fHistITShitsinlayer(0x0),
  fHistSPDNclsvstracklets(0x0),
  fHistITSchi2overhits(0x0),
  fHistMulthybridvsglobal(0x0),
  fHistTPCNcls(0x0),
  fHistTPCchi2(0x0),
  fHistTPCchi2overtpchits(0x0),
  // fHistVertex(0x0),
  fHistImpactParameterglobalxy(0x0),
  fHistImpactParameterglobalz(0x0),
  fHistImpactParameternoSPDxy(0x0),
  fHistImpactParameternoSPDz(0x0),
  fHistImpactParameternoITSxy(0x0),
  fHistImpactParameternoITSz(0x0),
  fHistVertexXY(0x0),
  fHistVertexZ(0x0),
  fHistContribtoVtx(0x0),
  
  //antikt histos
  fHistnJets(0x0),
  fAntiktHistPt(0x0),
  fAntiktHistPtwTRD(0x0),
  fAntiktHistPtwoTRD(0x0),
  fAntiktHistEta(0x0),
  fAntiktHistEtavsPt(0x0),
  fAntiktHistPhi(0x0),
  fAntiktHistEtaPhi(0x0),
  fAntiktHistAreavsPt(0x0),
  fAntiktHistAreaoverCirclevsPt(0x0),
  fAntiktHistnTracks(0x0), 
  fAntiktHistz(0x0),
  fAntiktHistzpt(0x0),
  fAntiktHistxi(0x0),
  fAntiktHistxipt(0x0),
  fAntiktHistxipt_pcp(0x0),  
  fAntiktHistxipt_pcm(0x0),
  fAntiktHistxipt_pc(0x0),
  fAntiktHistxipt_RC(0x0),
  fAntiktHistxiptwTRD(0x0),
  fAntiktHistxiptwoTRD(0x0),

  // fHistnJetsLeading(0x0),
  
  //random cone histos
  fHistRCntry(0x0),
  fHistPhiRC(0x0),
  fHistEtaRC(0x0),
  fHistRCnTracksOverlap(0x0),
  fHistRCPt(0x0),
  fHistRCnTracks(0x0),
  fHistRCnTracksvsLJPt(0x0),
  fHistRCPtvsLJPt(0x0),
 
  //Leading jet histos:
  fAntiktHistPtLeading(0x0),
  fAntiktHistPtLeadingvsHybridMult(0x0),
  fAntiktHistEtaLeading(0x0),
  fAntiktHistPhiLeading(0x0),
  fAntiktHistEtaPhiLeading(0x0),
  fAntiktHistAreaLeading(0x0),
  fAntiktHistAreaoverCircleLeading(0x0),
  fAntiktHistAreaLeadingvsPt(0x0),
  fAntiktHistAreaoverCircleLeadingvsPt(0x0),
  fAntiktHistAreaoverCircleLeadingvsNtracks(0x0),
  fAntiktHistnTracksLeading(0x0),
  fAntiktHistnTracksLeadingvsPt(0x0),
  fAntiktHistTrackPtLeading(0x0),
  fAntiktHistTrackPhiLeading(0x0),
  fAntiktHistTrackEtaLeading(0x0),
  fAntiktHistzLeading(0x0),
  fAntiktHistzptLeading(0x0),
  fAntiktHistxiLeading(0x0),
  fAntiktHistxiptLeading(0x0),
  fAntiktHistxitimesptptLeading(0x0),

  //Perpendicular Cone histos:
  //Plus histos
  fLJPCPTrackPt(0x0),
  fLJPCPTrackEta(0x0),
  fLJPCPTrackPhi(0x0),
  fLJPCPTrackEtaPhi(0x0),
  fLJPCPnTracks(0x0),
  fLJPCPnTracksPhi(0x0),
  fLJPCPPt(0x0),
  fLJPCPPtvsLJPt(0x0),
  fLJPCPEta(0x0),
  fLJPCPPhi(0x0),
  fLJPCPzpt(0x0),
  fLJPCPxipt(0x0),

  //Minus histos:

  fLJPCMTrackPt(0x0),
  fLJPCMTrackEta(0x0),  
  fLJPCMTrackPhi(0x0), 
  fLJPCMTrackEtaPhi(0x0), 
  fLJPCMnTracks(0x0),
  fLJPCMnTracksPhi(0x0),
  fLJPCMPt(0x0),
  fLJPCMPtvsLJPt(0x0),
  fLJPCMEta(0x0),
  fLJPCMPhi(0x0),
  fLJPCMzpt(0x0),
  fLJPCMxipt(0x0),

  //histos using both cones:
  fLJPCPMMPt(0x0),
  fLJPCPt(0x0),
  fLJPCPtvsLJPt(0x0),

  //Perp cone subtracted Histos

  fAntiktHistPtLeading_sub_pcp(0x0),
  fAntiktHistxiptLeading_sub_pcp(0x0),
  fAntiktHistTrackPtLeading_sub_pcp(0x0),

  fAntiktHistxiptLeading_pcp_sub_pcp(0x0),
  fAntiktHistTrackPtLeading_pcp_sub_pcp(0x0),
  
  //Cone around leading jet method
  fAntiktHistxiptLeading_pcp_LJCM(0x0),
  fAntiktHistPtLeading_LJCM(0x0),
  fAntiktHistxiptLeading_LJCM(0x0),
  fAntiktHistTrackPtLeading_LJCM(0x0),

  //leading track cut
  fAntiktHistxiptLeading_LTC(0x0),
  fAntiktHistxiptLeading_LTC_pcp(0x0),
  fAntiktHistPtLeading_LTC(0x0),

  //UE pT subtraction
  fHistRho(0x0),
  fHistRhovsLJPt(0x0),
  fHistRhoSubPt(0x0),
  fHistRhoSubAreaOverCirclevsPt(0x0),
  fHistRhoSubPhiTracksJetvsPt(0x0),
  fHistRhoSubEtaTracksJetvsPt(0x0),
  fHistRhoSubPhiTracksPCvsPt(0x0),
  fHistRhoSubEtaTracksPCvsPt(0x0),
  fHistRhoSubXiPt(0x0),
  fHistRhoSubXiPCPt(0x0),
  
  //RC RHO Histos
  fHistRhoRC(0x0),
  fHistRhoRCvsLJPt(0x0),
  fHistRhoRCSubPt(0x0),
  fHistRhoRCSubAreaOverCirclevsPt(0x0),
  fHistRhoRCSubPhiTracksJetvsPt(0x0),
  fHistRhoRCSubEtaTracksJetvsPt(0x0),
  fHistRhoRCSubPhiTracksPCvsPt(0x0),
  fHistRhoRCSubEtaTracksPCvsPt(0x0),
  fHistRhoRCSubXiPt(0x0),
  fHistRhoRCSubXiPCPt(0x0),

  //Delta pT histos
  fHistDeltaPt(0x0),
  fHistDeltaPtvsLJPt(0x0)



  // default constructor (with optional name)
 {                                                                                                                                                                                       
  DefineOutput(1, TList::Class() );
 }

AliAnalysisTaskFragmentationTriggered::~AliAnalysisTaskFragmentationTriggered()
{
  //destructor
  delete fMCRC;
  delete fMCjetstracks;
  delete fMCjets;
  delete fMCtracks_pcp;
  delete fMCtracks_pcm;
  delete fRC;
  delete fJets;
  delete fJetsTracks;
  delete fTracks_pcp;
  delete fTracks_pcm;
  delete fTracks_pcp_LJCM;
  delete fTracks_LJCM;
  delete fCuts;
  delete fCuts_noSPD;
  delete fCuts_noITS;
  delete fAnaUtils;
  delete fMCTruthParticles;
  delete fRandomGen;
 
   }

void AliAnalysisTaskFragmentationTriggered::UserCreateOutputObjects()
{
  //call initialize function
  Initialize();

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  fHistSumw2Debug =  new TH1F("Sumw2Debug", "Error should be 0.1", 5, -0.5, 4.5);
  TAxis *SDxaxis = fHistSumw2Debug->GetXaxis();

  SDxaxis->SetBinLabel(1, "100 entries scaled by 1/100");
  fOutputList->Add(fHistSumw2Debug);


  // track histos:
  fHistCuts = new TH1F("Cuts", "Number of Events passing Cut", 10, -0.5, 9.5);
  TAxis *xaxis = fHistCuts->GetXaxis();

  xaxis->SetBinLabel(1, "NoCut");  
  xaxis->SetBinLabel(2, "After first Trigger selection");
  xaxis->SetBinLabel(3, "After second Trigger selection");
  xaxis->SetBinLabel(4, "AnaUtils in-bunch pileup");
  xaxis->SetBinLabel(5, "SPD tracklet vs cluster PU");
  xaxis->SetBinLabel(6, "out of bunch PU");
  xaxis->SetBinLabel(7, "Vertex exists");
  xaxis->SetBinLabel(8, "> 3 contrib tracks to vtx");
  xaxis->SetBinLabel(9, "z<=10cm");
  
  fOutputList->Add(fHistCuts);  
  
  fHistEventDebug = new TH1F("EventDebug", "Number of existing Events", 5, -0.5, 4.5);
  TAxis *EDxaxis = fHistEventDebug->GetXaxis();

  EDxaxis->SetBinLabel(1, "MCEvent");
  EDxaxis->SetBinLabel(2, "ESDEvent");
  EDxaxis->SetBinLabel(3, "AODEvent");
  fOutputList->Add(fHistEventDebug);



  fHistTrackCuts = new TH1F("NTracksCuts", "Number of Tracks after Cuts", 10, -0.5, 9.5 );
  TAxis * trackcutsx = fHistTrackCuts->GetXaxis();

  trackcutsx->SetBinLabel(1, "NoCut");
  trackcutsx->SetBinLabel(2, "EtaCut"); 
  trackcutsx->SetBinLabel(3, "pT_range");
  trackcutsx->SetBinLabel(4, "TPCNCluster");
  trackcutsx->SetBinLabel(5, "DCA-xy");
  trackcutsx->SetBinLabel(6, "DCA-z");
  trackcutsx->SetBinLabel(7, "TPC-chi2overhits");
  trackcutsx->SetBinLabel(8, "ITS-chi2overhits");


  trackcutsx->SetBinLabel(10, "AllCuts");
  fOutputList->Add(fHistTrackCuts);


  fHistTrackCutsDep = new TH1F("NTracksCutsDep", "Number of Tracks after Cuts (dependent)", 10, -0.5, 9.5 );
  TAxis * trackcutsdepx = fHistTrackCutsDep->GetXaxis();

  trackcutsdepx->SetBinLabel(1, "NoCut");
  trackcutsdepx->SetBinLabel(2, "EtaCut");
  trackcutsdepx->SetBinLabel(3, "pT_range");
  trackcutsdepx->SetBinLabel(4, "TPCNCluster");
  trackcutsdepx->SetBinLabel(5, "DCA-xy");
  trackcutsdepx->SetBinLabel(6, "DCA-z");
  trackcutsdepx->SetBinLabel(7, "TPC-chi2overhits");
  trackcutsdepx->SetBinLabel(8, "ITS-chi2overhits");
  
  trackcutsdepx->SetBinLabel(10, "AllCuts");

  fOutputList->Add(fHistTrackCutsDep);

  //particle level
  fHistMCTruthPtHard = new TH1F("MCTruthPtHard", "Generated pT hard; p_{T} (GeV/#it{c});Counts", 350, -0.5, 349.5 );
  fOutputList->Add(fHistMCTruthPtHard);

  fHistMCTruthPtHardTrials = new TH1F("MCTruthPtHardTrials", "Generated pT hard weighted with nTrials;  p_{T} (GeV/#it{c});Counts", 350, -0.5, 349.5);
  fOutputList->Add(fHistMCTruthPtHardTrials);

  fHistMCTruthNTracks = new TH1F("MCTruthNTracks","Number of tracks in event from MC Truth; Number of tracks; Counts", 200, -0.5, 199.5);
  fOutputList->Add(fHistMCTruthNTracks);

  fHistMCTruthTrackPt = new TH1F("MCTruthTrackPt", "Transverse Momentum of Tracks;p_{T} (GeV/#it{c});Counts", 400, 0., 200.);
  fOutputList->Add(fHistMCTruthTrackPt);
   
  fHistMCTruthTrackPhi = new TH1F("MCTruthTrackPhi","Phi of Tracks;#phi;Counts",100,0, 2*TMath::Pi() );
  fOutputList->Add(fHistMCTruthTrackPhi);

  fHistMCTruthTrackEta = new TH1F("MCTruthTrackEta","Eta of Tracks;#eta;Counts",100,-5,5);
  fOutputList->Add(fHistMCTruthTrackEta);
    
  fMCTruthRecTrackPtGen = new TH1F("MCTruthRecTrackPtGen", "Generated Transverse Momentum of reco Tracks;p_{T} (GeV/#it{c});Counts", 400, 0., 200.);
  fOutputList->Add(fMCTruthRecTrackPtGen);

  fMCTruthRecTrackSecPtGen = new TH1F("MCTruthRecTrackSecPtGen", "Generated Transverse Momentum of secondaries;p_{T} (GeV/#it{c});Counts", 400, 0., 200.);
  fOutputList->Add(fMCTruthRecTrackSecPtGen);

  fMCTruthRecTrackSecPtRec = new TH1F("MCTruthRecTrackSecPtRec", "Reco Transverse Momentum of secondaries;p_{T} (GeV/#it{c});Counts", 400, 0., 200.);
  fOutputList->Add(fMCTruthRecTrackSecPtRec);

  fHistMCTruthTrackPtRatio = new TH2F("MCTruthTrackPtRatio", "Ratio of reconstructed track pt and generated pt", 200, 0, 2, 200, 0, 100);
  fOutputList->Add(fHistMCTruthTrackPtRatio);

  fHistMCTruthTrackPtRatiowSPD = new TH2F("MCTruthTrackPtRatiowSPD", "reconstructed track pt / generated pt", 200, 0, 2, 200, 0, 100);
  fOutputList->Add(fHistMCTruthTrackPtRatiowSPD);

  fHistMCTruthTrackPtRationoSPD = new TH2F("MCTruthTrackPtRationoSPD", "reconstructed track pt / generated pt", 200, 0, 2, 200, 0, 100);
  fOutputList->Add(fHistMCTruthTrackPtRationoSPD);

  //MCTruth Jets
  fMCTruthnJets = new TH1F("MCTruthnJets", "Number of jets above threshold and in acceptance", 20, -0.5, 19.5);
  fOutputList->Add(fMCTruthnJets);

  fMCTruthAntiktHistPt = new TH1F("MCTruthAntiktHistPt", "Transverse Momentum of Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
			   400, 0., 200.);
  fOutputList->Add(fMCTruthAntiktHistPt);

  fMCTruthAntiktHistAreavsPt = new TH2F("MCTruthAntiktHistAreavsPt", "Jetarea vs Pt; Area; p_{T} (GeV/#it{c})", 100, 0, 1, 40, 0, 200);
  fOutputList->Add(fMCTruthAntiktHistAreavsPt);

  fMCTruthAntiktHistAreaoverCirclevsPt = new TH2F("MCTruthAntiktHistAreaoverCirclevsPt", "Jetarea divided by Pi*RR; Area/Circle; p_{T} (GeV/#it{c})", 100, 0.5, 1.5, 40, 0, 200);
  fOutputList->Add(fMCTruthAntiktHistAreaoverCirclevsPt);

  fMCTruthAntiktHistxipt = new TH2F("MCTruthAntiktHistxipt", "#xi vs p_{T} of charged Tracks in Jets; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fMCTruthAntiktHistxipt);

  fMCTruthAntiktHistxipt_pcp =  new TH2F("MCTruthAntiktHistxipt_pcp", "#xi vs p_{T} of Tracks in perp cone plus; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fMCTruthAntiktHistxipt_pcp);

  fMCTruthAntiktHistxipt_pcm =  new TH2F("MCTruthAntiktHistxipt_pcm", "#xi vs p_{T} of Tracks in perp cone minus; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fMCTruthAntiktHistxipt_pcm);

  fMCTruthAntiktHistxipt_pc =  new TH2F("MCTruthAntiktHistxipt_pc", "#xi vs p_{T} of Tracks in perp cone plus; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fMCTruthAntiktHistxipt_pc);

  fMCTruthAntiktHistPtLeading = new TH1F("MCTruthAntiktHistPtLeading", "Transverse Momentum of leading Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
				  400, 0., 200.);
  fOutputList->Add(fMCTruthAntiktHistPtLeading);

  fMCTruthAntiktHistEtaLeading = new TH1F("MCTruthAntiktHistEtaLeading","Eta of leading Jets; #eta^{Jet,ch}; Counts",100,-1, 1);
  fOutputList->Add(fMCTruthAntiktHistEtaLeading);

  fMCTruthAntiktHistPhiLeading = new TH1F("MCTruthAntiktHistPhiLeading","Phi of leading Jets; #phi^{Jet,ch}; Counts", 100, 0, 2*TMath::Pi());
  fOutputList->Add(fMCTruthAntiktHistPhiLeading);

  fMCTruthAntiktHistxiptLeading =  new TH2F("MCTruthAntiktHistxiptLeading", "#xi vs p_{T} of charged Tracks in leading Jets; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fMCTruthAntiktHistxiptLeading);

  fMCTruthAntiktHistxiptLeading_pcp = new TH2F("MCTruthAntiktHistxiptLeading_pcp", "#xi vs p_{T} of charged Tracks in pcp; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fMCTruthAntiktHistxiptLeading_pcp);

  fMCTruthRCntry = new TH1F("MCTruthRCntry", "Number of thrown random cones until no overlap with signal jets", 50, -0.5, 49.5);
  fOutputList->Add(fMCTruthRCntry);

  fMCTruthAntiktHistPhiRC = new TH1F("MCTruthAntiktHistPhiRC","Phi of Jets; #phi^{Jet,ch}; Counts", 100, 0, 2*TMath::Pi());
  fOutputList->Add(fMCTruthAntiktHistPhiRC);

  fMCTruthAntiktHistEtaRC = new TH1F("MCTruthAntiktHistEtaRC","Eta of Jets; #eta^{Jet,ch}; Counts",100,-0.5,0.5);
  fOutputList->Add(fMCTruthAntiktHistEtaRC);

  fMCTruthRCnTracksOverlap = new TH1F("MCTruthRCnTracksOverlap", "Number of overlapping tracks from signal jets with RC; Number of Tracks; Counts", 20, -0.5, 19.5);
  fOutputList->Add(fMCTruthRCnTracksOverlap);

  fMCTruthHistRCPt = new TH1F("MCTruthHistRCPt", "Pt in random cone;p_{T} (GeV/#it{c});Counts", 400, 0, 20);
  fOutputList->Add(fMCTruthHistRCPt);

  fMCTruthHistRCnTracks = new TH1F("MCTruthHistRCnTracks", "Number of tracks in random cone; Number of tracks; Counts", 20, -0.5, 19.5);
  fOutputList->Add(fMCTruthHistRCnTracks);

  fMCTruthHistRCWithOverlapPt = new TH1F("MCTruthHistRCWithOverlapPt", "Pt in random cones that have overlap with signal jet; p_{T} (GeV/#it{c});Counts", 400, 0, 20);
  fOutputList->Add(fMCTruthHistRCWithOverlapPt);

  fMCTruthHistRCWithOverlapnTracks = new TH1F("MCTruthHistRCWithOverlapnTracks", "Number of tracks in random cone; Number of tracks; Counts", 20, -0.5, 19.5);
  fOutputList->Add(fMCTruthHistRCWithOverlapnTracks);

  fMCTruthHistRCTrackEta = new TH1F("MCTruthHistRCEta","Eta of tracks in random cones; #eta^{Jet,ch}; Counts",100,-1, 1);
  fOutputList->Add(fMCTruthHistRCTrackEta);

  fMCTruthHistRCTrackPhi = new TH1F("MCTruthHistRCTrackPhi","Phi of tracks in random cones; #phi^{Jet,ch}; Counts", 100, 0, 2*TMath::Pi());
  fOutputList->Add(fMCTruthHistRCTrackPhi);

  fMCTruthHistRCTrackEtaDiscarded = new TH1F("MCTruthHistRCEtaDiscarded","Eta of tracks in discarded random cones; #eta^{Jet,ch}; Counts",100,-1, 1);
  fOutputList->Add(fMCTruthHistRCTrackEtaDiscarded);

  fMCTruthHistRCTrackPhiDiscarded = new TH1F("MCTruthHistRCTrackPhiDiscarded","Phi of tracks in discarded random cones; #phi^{Jet,ch}; Counts", 100, 0, 2*TMath::Pi());
  fOutputList->Add(fMCTruthHistRCTrackPhiDiscarded);

  fMCTruthHistRCxipt = new TH2F("MCTruthHistRCxipt", "#xi vs p_{T} of charged Tracks in random Cones using inclusive jet pt; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0,   10, 40 , 0, 200 );
  fOutputList->Add(fMCTruthHistRCxipt);

  fMCTruthHistRCxiptLeading = new TH2F("MCTruthHistRCxiptLeading", "#xi vs p_{T} of charged Tracks in random Cones using leading jet pt; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fMCTruthHistRCxiptLeading);

  fMCTruthHistRCxiptLeadingDiscarded = new TH2F("MCTruthHistRCxiptLeadingDiscarded", "#xi vs p_{T} of charged Tracks in discarded RCs using leading jet pt; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fMCTruthHistRCxiptLeadingDiscarded);

  //rho plots
  fMCTruthHistRho = new TH1F("MCTruthHistRho", "Rho for lj perp cones ;#rho (GeV/#it{c});Counts", 1000, 0., 50.);
  fOutputList->Add(fMCTruthHistRho);

  fMCTruthHistRhovsLJPt = new TH2F("MCTruthHistRhovsLJPt", "Rho for lj perp cones vs leading jet pT; #rho (GeV/#it{c}); p_{T} (GeV/#it{c})", 1000, 0., 50., 40, 0, 200);
  fOutputList->Add(fMCTruthHistRhovsLJPt);

  fMCTruthHistRhoSubPt = new TH1F("MCTruthHistRhoSubPt", "Rho subtracted transverse Momentum of Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
  400, 0., 200.);
  fOutputList->Add(fMCTruthHistRhoSubPt);

  fMCTruthHistRhoSubAreaOverCirclevsPt = new TH2F("MCTruthHistRhoSubAreaOverCirclevsPt", "Jetarea divided by Pi*RR; Area/Circle; p_{T} (GeV/#it{c})", 100, 0.5, 1.5, 40, 0, 200);
  fOutputList->Add(fMCTruthHistRhoSubAreaOverCirclevsPt);

  fMCTruthHistRhoSubPhiTracksJetvsPt = new TH2F("MCTruthHistRhoSubPhiTracksJetvsPt", "Phi of Tracks in Jets vs Jet pt; #Phi^{Jet,ch}; Counts",100,0,TMath::Pi ()*2, 40, 0, 200);
  fOutputList->Add(fMCTruthHistRhoSubPhiTracksJetvsPt);

  fMCTruthHistRhoSubEtaTracksJetvsPt = new TH2F("MCTruthHistRhoSubEtaTracksJetvsPt", "Eta of Tracks in Jets vs Jet pt; #eta^{Jet,ch}; Counts",100,-1,1, 40, 0, 200);
  fOutputList->Add(fMCTruthHistRhoSubEtaTracksJetvsPt);

  fMCTruthHistRhoSubPhiTracksPCvsPt = new TH2F("MCTruthHistRhoSubPhiTracksPCvsPt", "Phi of Tracks in PCs vs Jet pt; #Phi^{Jet,ch}; Counts",100,0,TMath::Pi()*2, 40, 0, 200);
  fOutputList->Add(fMCTruthHistRhoSubPhiTracksPCvsPt);

  fMCTruthHistRhoSubEtaTracksPCvsPt = new TH2F("MCTruthHistRhoSubEtaTracksPCvsPt", "Eta of Tracks in PCs vs Jet pt; #eta^{Jet,ch}; Counts",100,-1,1, 40, 0, 200);
  fOutputList->Add(fMCTruthHistRhoSubEtaTracksPCvsPt);

  fMCTruthHistRhoSubXiPt = new TH2F("MCTruthHistRhoSubXiPt", "#xi vs p_{T}^{sub} of charged Tracks in Jets; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 ,   0, 200 );
  fOutputList->Add(fMCTruthHistRhoSubXiPt);

  fMCTruthHistRhoSubXiPCPt = new TH2F("MCTruthHistRhoSubXiPCPt", "#xi vs p_{T}^{sub} of charged Tracks in perp. cones; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50,   0, 10, 40 ,  0, 200 );
  fOutputList->Add(fMCTruthHistRhoSubXiPCPt);



  //detector level
  fHistPt = new TH1F("TrackPt", "Transverse Momentum of Tracks;p_{T} (GeV/#it{c});Counts", 400, 0., 200.);
  fOutputList->Add(fHistPt);

  fHistnTracks = new TH1F("HistnTracks","Number of global hybrid tracks in |#eta|<0.9 in event; Number of tracks; Counts", 200, -0.5, 199.5);
  fOutputList->Add(fHistnTracks);    

  fHistEta = new TH1F("TrackEta","Eta of Tracks;#eta;Counts",100,-1,1);
  fOutputList->Add(fHistEta);

  fHistPhi = new TH1F("TrackPhi","Phi of Tracks;#phi;Counts",100,0, 2*TMath::Pi() );
  fOutputList->Add(fHistPhi);

  fHistEtaPhi = new TH2F("TrackEtaPhi","Eta vs. Phi of Tracks;#eta; #phi", 100, -1, 1, 100, 0, 2*TMath::Pi());
  fOutputList->Add(fHistEtaPhi); 

  fHistPhiHybridUnconstrained = new TH1F("TrackPhiHybridUnconstrained","Phi of unconstrained hybrid Tracks;#phi;Counts",100,0, 2*TMath::Pi() );
  fOutputList->Add(fHistPhiHybridUnconstrained);

  fHistPhiHybridConstrainedwITS = new TH1F("TrackPhiHybridConstrainedwITS","Phi of constrained hybrid Tracks with ITS refit;#phi;Counts",100,0, 2*TMath::Pi() );
  fOutputList->Add(fHistPhiHybridConstrainedwITS);

  fHistPhiHybridConstrainedwoITS = new TH1F("TrackPhiHybridConstrainedwoITS","Phi of constrained hybrid Tracks without ITS refit (TPC only tracks);#phi;Counts",100,0, 2*TMath::Pi() );
  fOutputList->Add(fHistPhiHybridConstrainedwoITS);

  fHistITSNcls = new TH1F("TrackITSNHits","Number of ITS Hits of Tracks; Number of Hits; Counts", 10, -0.5, 9.5);
  fOutputList->Add(fHistITSNcls);

  fHistITShitsinlayer = new TH1F("TrackITShitsinlayer","Number of Tracks that have hit in layer no. ; Number of Layer; Counts", 10, -0.5, 9.5);
  fOutputList->Add(fHistITShitsinlayer);

  fHistSPDNclsvstracklets = new TH2F("SPDnclsvstracklets","Number of SPD cluster vs tracklets; Number of Tracklets; Number of Cluster", 100, 0, 200, 200, 0, 1000);
  fOutputList->Add(fHistSPDNclsvstracklets);

  fHistITSchi2overhits = new TH1F("TrackITSchi2overHits", "#chi^{2} / Number of ITS Hits; #chi^{2} / ITS Hits; Counts", 100 , 0, 50);
  fOutputList->Add(fHistITSchi2overhits);

  fHistMulthybridvsglobal = new TH2F("MultHybridvsGlobal","Hybrid multiplicity vs global multiplicity; Global Multiplicity |#eta|<0.9; Hybrid Multiplicity |#eta|<0.9", 100, 0, 200, 100, 0, 200);
  fOutputList->Add(fHistMulthybridvsglobal);

  fHistTPCNcls = new TH1F("TrackTPCNHits","Number of TPC Hits; Number of Hits; Counts", 200, -0.5, 199.5);
  fOutputList->Add(fHistTPCNcls);

  fHistTPCchi2 = new TH1F("TrackTPCchi2","#chi^{2} of TPC Tracks; #chi^{2}; Counts",100,0,100);
  fOutputList->Add(fHistTPCchi2);
 
  fHistTPCchi2overtpchits = new TH1F("TrackTPCchi2overTPCHits", "#chi^{2} / Number of TPC Hits; #chi^{2} / TPC Hits; Counts",100 , 0, 15);
  fOutputList->Add(fHistTPCchi2overtpchits);

  //fHistVertex = new TH3F("VertexXYZ","Position of Vertex; X/cm; Y/cm; Z/cm", 50, -0.5, 0.5, 50 , -0.5, 0.5, 100 , -20,  20);
  //fOutputList->Add(fHistVertex);
 
  fHistImpactParameterglobalxy=new TH1F("ImpactParameterofglobaltracksxy", "IP of tracks with SPD and ITS refit in XY-Plane; IP in XY-plane; Counts", 200, -5, 5 );
  fOutputList->Add(fHistImpactParameterglobalxy);

  fHistImpactParameterglobalz=new TH1F("ImpactParameterofglobaltracksz", "IP of tracks with SPD and ITS refit in Z Direction; IP in z-direction; Counts", 200, -10, 10);
  fOutputList->Add(fHistImpactParameterglobalz);

  fHistImpactParameternoSPDxy=new TH1F("ImpactParameterofnoSPDtracksxy", "IP of tracks wo SPD and w ITS refit in XY-Plane; IP in XY-plane; Counts", 200, -5, 5 );
  fOutputList->Add(fHistImpactParameternoSPDxy);

  fHistImpactParameternoSPDz=new TH1F("ImpactParameterofnoSPDtracksz", "IP of tracks wo SPD and w ITS refit in Z Direction; IP in z-direction; Counts", 200, -10, 10);
  fOutputList->Add(fHistImpactParameternoSPDz);

  fHistImpactParameternoITSxy=new TH1F("ImpactParameterofnoITStracksxy", "IP of tracks wo ITS refit in XY-Plane; IP in XY-plane; Counts", 200, -5, 5 );
  fOutputList->Add(fHistImpactParameternoITSxy);

  fHistImpactParameternoITSz=new TH1F("ImpactParameterofnoITStracksz", "IP of tracks wo ITS refit in Z Direction; IP in z-direction; Counts", 200, -10, 10);
  fOutputList->Add(fHistImpactParameternoITSz);
  
  fHistVertexXY = new TH2F("VertexXY", "Position of Vertex in XY Plane; X/cm; Y/cm", 200, -0.5, 0.5, 200, -0.5, 0.5);
  fOutputList->Add(fHistVertexXY);
 
  fHistVertexZ = new TH1F("VertexZ","Z Position of Vertex;Z Position/cm;Counts", 100, -20, 20);
  fOutputList->Add(fHistVertexZ); 

  fHistContribtoVtx = new TH1F("ContributerstoVertex","Number of Contributing Tracks to Primary Vertex; Number of Contributing Tracks; Counts", 100, 0, 200);
  fOutputList->Add(fHistContribtoVtx);
 
   //jet histos:

  fHistnJets = new TH1F("NumberofJetsperEvent", "Number of Jets per Event;Number of Jets; Counts", 50, -0.5 , 49.5);
  fOutputList->Add(fHistnJets);
  
  fAntiktHistPt = new TH1F("AntiktHistPt", "Transverse Momentum of Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
		     400, 0., 200.);
  fOutputList->Add(fAntiktHistPt);
    
  fAntiktHistPtwTRD = new TH1F("AntiktHistPtwTRD", "Transverse Momentum of Jets in TRD #phi range ;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
			   400, 0., 200.);
  fOutputList->Add(fAntiktHistPtwTRD);

  fAntiktHistPtwoTRD = new TH1F("AntiktHistPtwoTRD", "Transverse Momentum of Jets in no TRD #phi range ;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
			   400, 0., 200.);
  fOutputList->Add(fAntiktHistPtwoTRD);

  fAntiktHistEta = new TH1F("AntiktHistEta","Eta of Jets; #eta^{Jet,ch}; Counts",100,-0.5,0.5);
  fOutputList->Add(fAntiktHistEta);

  fAntiktHistEtavsPt = new TH2F("AntiktHistEtavsPt", "Eta of Jets; #eta^{Jet,ch}; Counts",100,-0.5,0.5, 40, 0, 200);
  fOutputList->Add(fAntiktHistEtavsPt);

  fAntiktHistPhi = new TH1F("AntiktHistPhi","Phi of Jets; #phi^{Jet,ch}; Counts", 100, 0, 2*TMath::Pi());
  fOutputList->Add(fAntiktHistPhi);

  fAntiktHistEtaPhi = new TH2F("AntiktHistEtaPhi","Eta vs Phi of Jets; #eta^{Jet,ch}; #phi^{Jet,ch}", 50, -0.5, 0.5, 50, 0, 2*TMath::Pi());
  fOutputList->Add(fAntiktHistEtaPhi); 

  fAntiktHistAreavsPt = new TH2F("AntiktHistAreavsPt", "Jetarea vs Pt; Area; p_{T} (GeV/#it{c})", 100, 0, 1, 40, 0, 200);
  fOutputList->Add(fAntiktHistAreavsPt);

  fAntiktHistAreaoverCirclevsPt = new TH2F("AntiktHistAreaoverCirclevsPt", "Jetarea divided by Pi*RR; Area/Circle; p_{T} (GeV/#it{c})", 100, 0.5, 1.5, 40, 0, 200);
  fOutputList->Add(fAntiktHistAreaoverCirclevsPt);

  fAntiktHistnTracks = new TH1F("AntiktHistnTracks","Number of Tracks in Jet; Number of Tracks; Counts", 100,-0.5, 99.5);
  fOutputList->Add(fAntiktHistnTracks);

  fAntiktHistz = new TH1F("AntiktHistz","z of charged Tracks in Jets; z=p_{T}^{particle}/p_{T}^{jet}; Counts", 50, 0, 1);
  fOutputList->Add(fAntiktHistz);

  fAntiktHistzpt = new TH2F("AntiktHistzpt","z vs p_{T} of charged Tracks in Jets; z=p_{T}^{particle}/p_{T}^{jet}; p_{T} (GeV/#it{c})", 50, 0, 1, 40, 0, 200);
  fOutputList->Add(fAntiktHistzpt);

  fAntiktHistxi = new TH1F("AntiktHistxi","#xi of charged Tracks in Jets; #xi=ln(1/z); Counts", 50, 0, 10);
  fOutputList->Add(fAntiktHistxi);

  fAntiktHistxipt = new TH2F("Antikthistxipt", "#xi vs p_{T} of charged Tracks in Jets; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxipt);
 
  fAntiktHistxipt_pcp = new TH2F("Antikthistxipt_pcp", "#xi vs jet p_{T} of charged Tracks in pcp; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxipt_pcp);

  fAntiktHistxipt_pcm = new TH2F("Antikthistxipt_pcm", "#xi vs jet p_{T} of charged Tracks in pcm; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxipt_pcm);

  fAntiktHistxipt_pc = new TH2F("Antikthistxipt_pc", "#xi vs jet p_{T} of charged Tracks in pc; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxipt_pc);

  fAntiktHistxipt_RC = new TH2F("Antikthistxipt_RC", "#xi vs jet p_{T} of charged Tracks in RC; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxipt_RC);

  fAntiktHistxiptwTRD = new TH2F("AntikthistxiptwTRD", "#xi vs p_{T} of charged Tracks in Jets in #phi areas with TRD; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxiptwTRD);

  fAntiktHistxiptwoTRD = new TH2F("AntikthistxiptwoTRD", "#xi vs p_{T} of charged Tracks in Jets in #phi areas without TRD; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxiptwoTRD);

  //===============================
  //random cone histos:
  //=============================

  fHistRCntry = new TH1F("HistRCntry", "Number of thrown random cones until no overlap with signal jets", 50, -0.5, 49.5);
  fOutputList->Add(fHistRCntry);

  fHistPhiRC = new TH1F("HistPhiRC","Phi of random cones; #phi; Counts", 100, 0, 2*TMath::Pi());
  fOutputList->Add(fHistPhiRC);

  fHistEtaRC = new TH1F("HistEtaRC","Eta of random cones; #eta; Counts",100,-0.5,0.5);
  fOutputList->Add(fHistEtaRC);

  fHistRCnTracksOverlap = new TH1F("HistRCnTracksOverlap", "Number of overlapping tracks from signal jets with RC; Number of Tracks; Counts", 20, -0.5, 19.5);
  fOutputList->Add(fHistRCnTracksOverlap);

  fHistRCPt = new TH1F("HistRCPt", "Pt in random cone;p_{T} (GeV/#it{c});Counts", 400, 0, 20);
  fOutputList->Add(fHistRCPt);

  fHistRCnTracks = new TH1F("HistRCnTracks", "Number of tracks in random cone; Number of tracks; Counts", 20, -0.5, 19.5);
  fOutputList->Add(fHistRCnTracks);
  
  fHistRCnTracksvsLJPt = new TH2F("HistRCnTracksvsLJPt", "Number of tracks in random conevs leading jet pT; Number of tracks; Counts", 20, -0.5, 19.5,  40, 0, 200);
  fOutputList->Add(fHistRCnTracksvsLJPt);

  fHistRCPtvsLJPt = new TH2F("HistRCPtvsLJPt", "Random Cone p_{T} vs leading jet p_{T};  p_{T} (GeV/#it{c}); p_{T} (GeV/#it{c})",  400, 0, 20, 40, 0, 200);
  fOutputList->Add(fHistRCPtvsLJPt);


  //========================================
  //leading jet histos:
  //=======================================

  fAntiktHistPtLeading = new TH1F("AntiktHistPtLeading", "Transverse Momentum of leading Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
			   400, 0., 200.);
  fOutputList->Add(fAntiktHistPtLeading);

  fAntiktHistPtLeadingvsHybridMult = new TH2F("AntiktHistPtLeadingvsHybridMult", "Transverse Momentum of leading Jets vs hybrid Multiplicity; p_{T}^{Jet,ch} (GeV/#it{c}); Hybrid Multiplicity", 200, 0, 200, 200, -0.5, 199.5);
  fOutputList->Add(fAntiktHistPtLeadingvsHybridMult);

  fAntiktHistEtaLeading = new TH1F("AntiktHistEtaLeading","Eta of leading Jets; #eta^{Jet,ch}; Counts",100,-0.5,0.5);
  fOutputList->Add(fAntiktHistEtaLeading);

  fAntiktHistPhiLeading = new TH1F("AntiktHistPhiLeading","Phi of leading Jets; #phi^{Jet,ch}; Counts", 100, 0, 2*TMath::Pi());
  fOutputList->Add(fAntiktHistPhiLeading);

  fAntiktHistEtaPhiLeading = new TH2F("AntiktHistEtaPhiLeading","Eta vs Phi of leading Jets; #eta^{Jet,ch}; #phi^{Jet,ch}", 50, -0.5, 0.5, 50, 0, 2*TMath::Pi());
  fOutputList->Add(fAntiktHistEtaPhiLeading);

  fAntiktHistAreaLeading = new TH1F("AntiktHistAreaLeading", "Area of leading Jet", 100, 0, 1);
  fOutputList->Add(fAntiktHistAreaLeading);

  fAntiktHistAreaoverCircleLeading = new TH1F("AntiktHistAreaoverCircleLeading", "Area of leading Jet divided by Area of Circle with Jet Radius", 100, 0.5, 1.5);
  fOutputList->Add(fAntiktHistAreaoverCircleLeading);

  fAntiktHistAreaLeadingvsPt = new TH2F("AntiktHistAreaLeadingvsPt", "Area of leading Jet vs Pt; Area; p_{T} (GeV/#it{c})", 100, 0, 1, 40, 0, 200);
  fOutputList->Add(fAntiktHistAreaLeadingvsPt);

  fAntiktHistAreaoverCircleLeadingvsPt = new TH2F("AntiktHistAreaoverCircleLeadingvsPt", "Area of leading Jet divided by Area of Circle with Jet Radius; Area/Circle; p_{T} (GeV/#it{c})", 100, 0.5, 1.5, 40, 0, 200);
  fOutputList->Add(fAntiktHistAreaoverCircleLeadingvsPt);

  fAntiktHistAreaoverCircleLeadingvsNtracks = new TH2F("AntiktHistAreaoverCircleLeadingvsNtracks", "Area of leading Jet divided by Area of Circle with Jet Radius; Area/Circle; N_{tracks}", 100, 0.5, 1.5, 20, -0.5, 19.5);
  fOutputList->Add(fAntiktHistAreaoverCircleLeadingvsNtracks);

  fAntiktHistnTracksLeading = new TH1F("AntiktHistnTracksLeading","Number of Tracks in leading Jet; Number of Tracks; Counts", 100,-0.5, 99.5);
  fOutputList->Add(fAntiktHistnTracksLeading);

  fAntiktHistnTracksLeadingvsPt = new TH2F("AntiktHistnTracksLeadingvsPt","Number of Tracks in leading Jet vs Pt; Number of Tracks; p_{T} (GeV/#it{c})", 100,-0.5, 99.5, 40, 0, 200);
  fOutputList->Add(fAntiktHistnTracksLeadingvsPt);

  fAntiktHistTrackPtLeading = new TH1F("AntiktHistTrackPtLeading", "Transverse Momentum of tracks in leading Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
				  400, 0., 200.);
  fOutputList->Add(fAntiktHistTrackPtLeading);

  fAntiktHistTrackPhiLeading = new TH1F("AntiktHistTrackPhiLeading","Phi of tracks in leading Jets; #phi^{Jet,ch}; Counts", 100, 0, 2*TMath::Pi());
  fOutputList->Add(fAntiktHistTrackPhiLeading);

  fAntiktHistTrackEtaLeading = new TH1F("AntiktHistTrackEtaLeading","Eta of tracks in leading Jets; #eta^{Jet,ch}; Counts", 100, -1, 1 );
  fOutputList->Add(fAntiktHistTrackEtaLeading);

  fAntiktHistzLeading = new TH1F("AntiktHistzLeading","z of charged Tracks in leading Jets; z=p_{T}^{particle}/p_{T}^{jet}; Counts", 50, 0, 1);
  fOutputList->Add(fAntiktHistzLeading);

  fAntiktHistzptLeading = new TH2F("AntiktHistzptLeading","z vs p_{T} of charged Tracks in leading Jets; z=p_{T}^{particle}/p_{T}^{jet}; p_{T} (GeV/#it{c})", 50, 0, 1, 40, 0, 200);
  fOutputList->Add(fAntiktHistzptLeading);

  fAntiktHistxiLeading = new TH1F("AntiktHistxiLeading","#xi of charged Tracks in leading Jets; #xi=ln(1/z); Counts", 50, 0, 10);
  fOutputList->Add(fAntiktHistxiLeading);

  fAntiktHistxiptLeading = new TH2F("AntikthistxiptLeading", "#xi vs p_{T} of charged Tracks in leading Jets; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxiptLeading);

  fAntiktHistxitimesptptLeading = new TH2F("AntiktHistxitimesptptLeading", "p_{T}^{track}*#xi vs p_{T}^{jet} of charged Tracks in leading Jets; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxitimesptptLeading);

  //================================================
  //perpendicular cones for leading jet
  //================================================

  //Histos for plus cone:
  //properties of tracks:
  fLJPCPTrackPt= new TH1F("LJPCPTrackPt","p_{T} of charged Tracks in perpendicular cone 'plus' wrt leading Jet; p_{T} (GeV/#it{c}); Counts", 400, 0, 20);
  fOutputList->Add(fLJPCPTrackPt);

  fLJPCPTrackEta = new TH1F("LJPCPTrackEta","Eta of Tracks in perpendicular cone 'plus' wrt leading Jet;#eta;Counts",100,-1,1);
  fOutputList->Add(fLJPCPTrackEta);

  fLJPCPTrackPhi = new TH1F("LJPCPTrackPhi","Phi of Tracks in perpendicular cone 'plus' wrt leading Jet ;#phi;Counts",100,0, 2*TMath::Pi() );
  fOutputList->Add(fLJPCPTrackPhi);

  fLJPCPTrackEtaPhi = new TH2F("LJPCPTrackEtaPhi","Eta vs Phi of Tracks in perpendicular cone 'plus' wrt leading Jet; #eta^{Jet,ch}; #phi^{Jet,ch}", 100, -1, 1, 100, 0, 2*TMath::Pi());
  fOutputList->Add(fLJPCPTrackEtaPhi);

  fLJPCPnTracks = new TH1F("LJPCPnTracks", "Number of charged Tracks in perpendicular cone 'plus' wrt leading Jet; Number of Tracks; Counts",100, -0.5, 99.5 );
  fOutputList->Add(fLJPCPnTracks);

  fLJPCPnTracksPhi = new TH2F("LJPCPnTracksPhi", "Number of charged Tracks in perpendicular cone 'plus' wrt leading Jet; Number of Tracks; #Phi",100, -0.5, 99.5, 100, -TMath::Pi()/2, 2.5*TMath::Pi() );
  fOutputList->Add(fLJPCPnTracksPhi);

  //properties of cone:
  fLJPCPPt = new TH1F("LJPCPPt", "p_{T} of perpendicular cone 'plus' wrt leading Jet; p_{T} (GeV/#it{c}); Counts", 400, 0, 20);
  fOutputList->Add(fLJPCPPt);

  fLJPCPPtvsLJPt = new TH2F("LJPCPPtvsLJPt", "p_{T} of perpendicular cone 'plus' wrt leading Jet vs leading jet pt ; p_{T} (GeV/#it{c}); p_{T} (GeV/#it{c}) ", 400, 0, 20, 40, 0, 200);
  fOutputList->Add(fLJPCPPtvsLJPt);

  fLJPCPEta = new TH1F("LJPCPEta","Eta of perpendicular cone 'plus' wrt leading Jet ;#eta;Counts",100,-1,1);
  fOutputList->Add(fLJPCPEta);

  fLJPCPPhi = new TH1F("LJPCPPhi","Phi of perpendicular cone 'plus' wrt leading Jet ;#phi;Counts",100, -TMath::Pi()/2, 2.5*TMath::Pi() );
  fOutputList->Add(fLJPCPPhi);
 
 //fragmentation functions:
  fLJPCPzpt = new TH2F("LJPCPzpt","z vs p_{T} of charged Tracks in perpendicular cone 'plus' wrt leading Jet ; z=p_{T}^{particle}/p_{T}^{jet}; p_{T} (GeV/#it{c})", 50, 0, 1, 40, 0, 200);
  fOutputList->Add(fLJPCPzpt);

  fLJPCPxipt = new TH2F("LJPCPxipt", "#xi vs p_{T} of charged Tracks in perpendicular cone 'plus' wrt leading Jet ; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fLJPCPxipt);


 
  //histos for minus cone:
  //properties of tracks:                                                                                                                                                                                                                    
  fLJPCMTrackPt= new TH1F("LJPCMTrackPt","p_{T} of charged Tracks in perpendicular cone 'minus' wrt leading Jet; p_{T} (GeV/#it{c}); Counts", 400, 0, 20);
  fOutputList->Add(fLJPCMTrackPt);

  fLJPCMTrackEta = new TH1F("LJPCMTrackEta","Eta of Tracks in perpendicular cone 'minus' wrt leading Jet;#eta;Counts",100,-1,1);
  fOutputList->Add(fLJPCMTrackEta);

  fLJPCMTrackPhi = new TH1F("fLJPCMTrackPhi","Phi of Tracks in perpendicular cone 'minus' wrt leading Jet ;#phi;Counts",100,0, 2*TMath::Pi() );
  fOutputList->Add(fLJPCMTrackPhi);

  fLJPCMTrackEtaPhi = new TH2F("LJPCMTrackEtaPhi","Eta vs Phi of Tracks in perpendicular cone 'minus' wrt leading Jet; #eta^{Jet,ch}; #phi^{Jet,ch}", 100, -1, 1, 100, 0, 2*TMath::Pi());
  fOutputList->Add(fLJPCMTrackEtaPhi);

  fLJPCMnTracks = new TH1F("LJPCMnTracks", "Number of charged Tracks in perpendicular cone 'minus' wrt leading Jet; Number of Tracks; Counts",100, -0.5, 99.5 );
  fOutputList->Add(fLJPCMnTracks);

  fLJPCMnTracksPhi = new TH2F("LJPCMnTracksPhi", "Number of charged Tracks in perpendicular cone 'minus' wrt leading Jet; Number of Tracks; #Phi",100, -0.5, 99.5, 100, -TMath::Pi()/2, 2.5*TMath::Pi() );
  fOutputList->Add(fLJPCMnTracksPhi);

  //properties of cone:                                                                                                                                                                                                                      
  fLJPCMPt = new TH1F("LJPCMPt", "p_{T} of perpendicular cone 'minus' wrt leading Jet; p_{T} (GeV/#it{c}); Counts", 400, 0, 20);
  fOutputList->Add(fLJPCMPt);

  fLJPCMPtvsLJPt = new TH2F("LJPCMPtvsLJPt", "p_{T} of perpendicular cone 'minus' wrt leading Jet vs leading jet pt ; p_{T} (GeV/#it{c}); p_{T} (GeV/#it{c}) ", 400, 0, 20, 40, 0, 200);
  fOutputList->Add(fLJPCMPtvsLJPt);

  fLJPCMEta = new TH1F("LJPCMEta","Eta of perpendicular cone 'minus' wrt leading Jet ;#eta;Counts",100,-1,1);
  fOutputList->Add(fLJPCMEta);

  fLJPCMPhi = new TH1F("LJPCMPhi","Phi of perpendicular cone 'minus' wrt leading Jet ;#phi;Counts",100,-TMath::Pi()/2, 2.5*TMath::Pi() );
  fOutputList->Add(fLJPCMPhi);

  //fragmentation functions:                                                                                                                                                                                                                  
  fLJPCMzpt = new TH2F("LJPCMzpt","z vs p_{T} of charged Tracks in perpendicular cone 'minus' wrt leading Jet ; z=p_{T}^{particle}/p_{T}^{jet}; p_{T} (GeV/#it{c})", 50, 0, 1, 40, 0, 200);
  fOutputList->Add(fLJPCMzpt);

  fLJPCMxipt = new TH2F("LJPCMxipt", "#xi vs p_{T} of charged Tracks in perpendicular cone 'minus' wrt leading Jet ; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fLJPCMxipt);
 
   //mean of both cone pts
  fLJPCPMMPt = new TH1F("LJPCPMMPt", "mean of p_{T} of perpendicular cone 'minus' and 'plus'  wrt leading Jet; p_{T} (GeV/#it{c}); Counts", 400, 0, 20);
  fOutputList->Add(fLJPCPMMPt);

  //Plots using both cones:
  fLJPCPt = new TH1F("LJPCPt", "p_{T} of perpendicular cone wrt leading Jet including both cones; p_{T} (GeV/#it{c}); Counts", 400, 0, 20);
  fOutputList->Add(fLJPCPt);

  fLJPCPtvsLJPt = new TH2F("LJPCPtvsLJPt", "p_{T} of perpendicular cone vs leading jet pt ; p_{T} (GeV/#it{c}); p_{T} (GeV/#it{c}) ", 400, 0, 20, 40, 0, 200);
  fOutputList->Add(fLJPCPtvsLJPt);


  //perp cone subtracted plots

  fAntiktHistPtLeading_sub_pcp = new TH1F("AntiktHistPtLeading_sub_pcp", "Transverse Momentum of leading Jets subtracted by PCP ;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
				  400, 0., 200.);
  fOutputList->Add(fAntiktHistPtLeading_sub_pcp);

  fAntiktHistxiptLeading_sub_pcp = new TH2F("AntikthistxiptLeading_sub_pcp", "#xi vs p_{T} of charged Tracks in leading Jets subtracted by pcp pt; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxiptLeading_sub_pcp);

  fAntiktHistTrackPtLeading_sub_pcp = new TH1F("AntiktHistTrackPtLeading_sub_pcp", "Transverse Momentum of tracks in leading Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
				       400, 0., 200.);
  fOutputList->Add(fAntiktHistTrackPtLeading_sub_pcp);

  fAntiktHistxiptLeading_pcp_sub_pcp = new TH2F("AntikthistxiptLeading_pcp_sub_pcp", "#xi vs p_{T} of charged Tracks in leading Jets subtracted by pcp pt; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxiptLeading_pcp_sub_pcp);

  fAntiktHistTrackPtLeading_pcp_sub_pcp = new TH1F("AntiktHistTrackPtLeading_pcp_sub_pcp", "Transverse Momentum of tracks in leading Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
					       400, 0., 200.);
  fOutputList->Add(fAntiktHistTrackPtLeading_pcp_sub_pcp);



  //leading jet cone method histos

  fAntiktHistxiptLeading_pcp_LJCM = new TH2F("AntikthistxiptLeading_pcp_LJCM", "#xi vs p_{T} of charged Tracks in perp Cone for jet pt  with Cone method; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxiptLeading_pcp_LJCM);

  fAntiktHistPtLeading_LJCM = new TH1F("AntiktHistPtLeading_LJCM", "Transverse Momentum of leading Jets with Cone method ;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
					  400, 0., 200.);
  fOutputList->Add(fAntiktHistPtLeading_LJCM);

  fAntiktHistxiptLeading_LJCM = new TH2F("AntikthistxiptLeading_LJCM", "#xi vs p_{T} of charged Tracks in leading Jets with cone method; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxiptLeading_LJCM);

  fAntiktHistTrackPtLeading_LJCM = new TH1F("AntiktHistTrackPtLeading_LJCM", "Transverse Momentum of tracks in leading Jets with cone method;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
					       400, 0., 200.);
  fOutputList->Add(fAntiktHistTrackPtLeading_LJCM);

  //leading track cut

  fAntiktHistxiptLeading_LTC = new TH2F("AntikthistxiptLeading_LTC", "#xi vs p_{T} of charged Tracks in leading Jets with leading track cut; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxiptLeading_LTC);

  fAntiktHistxiptLeading_LTC_pcp = new TH2F("AntikthistxiptLeading_LTC_pcp", "#xi vs p_{T} of charged Tracks in perp cone plus wrt leading Jet with leading track cut; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fAntiktHistxiptLeading_LTC_pcp); 

  fAntiktHistPtLeading_LTC = new TH1F("AntiktHistPtLeading_LTC", "Transverse Momentum of leading Jets with leading track cut ;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
				       400, 0., 200.);
  fOutputList->Add(fAntiktHistPtLeading_LTC);

  //plots for UE pT subtraction of perp cones for inclusive jets:
  fHistRho = new TH1F("HistRho", "Rho for lj perp cones ;#rho (GeV/#it{c});Counts", 1000, 0., 50.);
  fOutputList->Add(fHistRho);

  fHistRhovsLJPt = new TH2F("HistRhovsLJPt", "Rho for lj perp cones vs leading jet pT; #rho (GeV/#it{c}); p_{T} (GeV/#it{c})", 1000, 0., 50., 40, 0, 200);
  fOutputList->Add(fHistRhovsLJPt);

  fHistRhoSubPt = new TH1F("HistRhoSubPt", "Rho subtracted transverse Momentum of Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
			   400, 0., 200.);
  fOutputList->Add(fHistRhoSubPt);

  fHistRhoSubAreaOverCirclevsPt = new TH2F("HistRhoSubAreaOverCirclevsPt", "Jetarea divided by Pi*RR; Area/Circle; p_{T} (GeV/#it{c})", 100, 0.5, 1.5, 40, 0, 200);
  fOutputList->Add(fHistRhoSubAreaOverCirclevsPt);

  fHistRhoSubPhiTracksJetvsPt = new TH2F("HistRhoSubPhiTracksJetvsPt", "Phi of Tracks in Jets vs Jet pt; #Phi^{Jet,ch}; Counts",100,0,TMath::Pi()*2, 40, 0, 200);
  fOutputList->Add(fHistRhoSubPhiTracksJetvsPt);

  fHistRhoSubEtaTracksJetvsPt = new TH2F("HistRhoSubEtaTracksJetvsPt", "Eta of Tracks in Jets vs Jet pt; #eta^{Jet,ch}; Counts",100,-1,1, 40, 0, 200);
  fOutputList->Add(fHistRhoSubEtaTracksJetvsPt);

  fHistRhoSubPhiTracksPCvsPt = new TH2F("HistRhoSubPhiTracksPCvsPt", "Phi of Tracks in PCs vs Jet pt; #Phi^{Jet,ch}; Counts",100,0,TMath::Pi()*2, 40, 0, 200);
  fOutputList->Add(fHistRhoSubPhiTracksPCvsPt);

  fHistRhoSubEtaTracksPCvsPt = new TH2F("HistRhoSubEtaTracksPCvsPt", "Eta of Tracks in PCs vs Jet pt; #eta^{Jet,ch}; Counts",100,-1,1, 40, 0, 200);
  fOutputList->Add(fHistRhoSubEtaTracksPCvsPt);

  fHistRhoSubXiPt = new TH2F("HistRhoSubXiPt", "#xi vs p_{T}^{sub} of charged Tracks in Jets; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fHistRhoSubXiPt);

  fHistRhoSubXiPCPt = new TH2F("HistRhoSubXiPCPt", "#xi vs p_{T}^{sub} of charged Tracks in perp. cones; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 ,  0, 200 );
  fOutputList->Add(fHistRhoSubXiPCPt);


  //plots for UE pT subtraction of RC:                                                                                                                                
  fHistRhoRC = new TH1F("HistRhoRC", "Rho for RC ;#rho (GeV/#it{c});Counts", 1000, 0., 50.);
  fOutputList->Add(fHistRhoRC);

  fHistRhoRCvsLJPt = new TH2F("HistRhoRCvsLJPt", "Rho for RCs vs leading jet pT; #rho (GeV/#it{c}); p_{T} (GeV/#it{c})", 1000, 0., 50., 40, 0, 200);
  fOutputList->Add(fHistRhoRCvsLJPt);

  fHistRhoRCSubPt = new TH1F("HistRhoRCSubPt", "Rho RC subtracted transverse Momentum of Jets;p_{T}^{Jet,ch} (GeV/#it{c});Counts",
                           400, 0., 200.);
  fOutputList->Add(fHistRhoRCSubPt);

  fHistRhoRCSubAreaOverCirclevsPt = new TH2F("HistRhoRCSubAreaOverCirclevsPt", "Jetarea divided by Pi*RR; Area/Circle; p_{T} (GeV/#it{c})", 100, 0.5, 1.5, 40, 0, 200);
  fOutputList->Add(fHistRhoRCSubAreaOverCirclevsPt);
  
  fHistRhoRCSubPhiTracksJetvsPt = new TH2F("HistRhoRCSubPhiTracksJetvsPt", "Phi of Tracks in Jets vs Jet pt; #Phi^{Jet,ch}; Counts",100,0,TMath::Pi()*2, 40, 0, 200);
  fOutputList->Add(fHistRhoRCSubPhiTracksJetvsPt);

  fHistRhoRCSubEtaTracksJetvsPt = new TH2F("HistRhoRCSubEtaTracksJetvsPt", "Eta of Tracks in Jets vs Jet pt; #eta^{Jet,ch}; Counts",100,-1,1, 40, 0, 200);
  fOutputList->Add(fHistRhoRCSubEtaTracksJetvsPt);

  fHistRhoRCSubPhiTracksPCvsPt = new TH2F("HistRhoRCSubPhiTracksPCvsPt", "Phi of Tracks in PCs vs Jet pt; #Phi^{Jet,ch}; Counts",100,0,TMath::Pi()*2, 40, 0, 200);
  fOutputList->Add(fHistRhoRCSubPhiTracksPCvsPt);

  fHistRhoRCSubEtaTracksPCvsPt = new TH2F("HistRhoRCSubEtaTracksPCvsPt", "Eta of Tracks in PCs vs Jet pt; #eta^{Jet,ch}; Counts",100,-1,1, 40, 0, 200);
  fOutputList->Add(fHistRhoRCSubEtaTracksPCvsPt);

  fHistRhoRCSubXiPt = new TH2F("HistRhoRCSubXiPt", "#xi vs p_{T}^{sub} of charged Tracks in Jets; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 , 0, 200 );
  fOutputList->Add(fHistRhoRCSubXiPt);

  fHistRhoRCSubXiPCPt = new TH2F("HistRhoRCSubXiPCPt", "#xi vs p_{T}^{sub} of charged Tracks in perp. cones; #xi=ln(1/z); p_{T} (GeV/#it{c})", 50, 0, 10, 40 ,  0, 200 );
  fOutputList->Add(fHistRhoRCSubXiPCPt);
  

  //delta pT plots
  fHistDeltaPt = new TH1F("HistDeltaPt", "Delta pT ;#delta p_{T} (GeV/#it{c});Counts", 1000, -50, 50.);
  fOutputList->Add(fHistDeltaPt);

  fHistDeltaPtvsLJPt = new TH2F("HisDeltaPtvsLJPt", "Delta pT vs leading jet pT; #delta p_{T} (GeV/#it{c}); p_{T} (GeV/#it{c})", 1000, -50, 50., 40, 0, 200);
  fOutputList->Add(fHistDeltaPtvsLJPt);



  //now: loop over entries in fOutputList and create Sumw2 for all histograms (important for adding pt hard bins)



  for (Int_t iCount = 0; iCount < fOutputList->GetEntries(); iCount++ ) {

    if (fOutputList->At(iCount)->InheritsFrom("TH2") ){
      TH2F *hist =(TH2F*)fOutputList->At(iCount);

      //std::cout << "this object inherits from th2" << endl;

      hist->Sumw2(); 
    }

    else {
      TH1F *hist2 =(TH1F*)fOutputList->At(iCount);

      //std::cout << "this is a TH1" << endl;
      hist2->Sumw2();
    }
  }//close loop


  //fill sumw2 debug histo
  for (Int_t i = 0; i < 100; i++){
  fHistSumw2Debug->Fill(1);
  }

  fHistSumw2Debug->Scale(0.01);




  PostData(1, fOutputList);
}



void AliAnalysisTaskFragmentationTriggered::Initialize()
{

  //monte carlo truth objects:
  fMCRC = new TObjArray();
  fMCjetstracks = new TObjArray();
  fMCjets = new TObjArray();
  fMCtracks_pcp = new TObjArray();
  fMCtracks_pcm = new TObjArray();

  //reconstructed objects:
  fRC = new TObjArray();
  fJets = new TObjArray();
  fJetsTracks = new TObjArray();
  fTracks_pcp = new TObjArray();
  fTracks_pcm = new TObjArray();
  fTracks_pcp_LJCM = new TObjArray();
  fTracks_LJCM = new TObjArray();

  //Define set of trackcuts: 
  fCuts = new AliESDtrackCuts();
  fCuts_noSPD=new AliESDtrackCuts();
  fCuts_noITS=new AliESDtrackCuts();

  //other stuff
  fAnaUtils = new AliAnalysisUtils();
  fMCTruthParticles = new TObjArray();
  fRandomGen = new TRandom3();

  //set random seed for random generator:
  fRandomGen->SetSeed(0);

  //GOOD GLOBAL TRACKS WITH SPD 
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  fCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
  fCuts->SetMinNClustersTPC(70);
  fCuts->SetMaxChi2PerClusterTPC(4);
  fCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1                                                                                            
  fCuts->SetAcceptKinkDaughters(kFALSE);
  fCuts->SetRequireTPCRefit(kTRUE);
  fCuts->SetMaxFractionSharedTPCClusters(0.4);
  // ITS                                                                                                                                                                    
  fCuts->SetRequireITSRefit(kTRUE);
  //accept secondaries                                                                                                                                                      
  fCuts->SetMaxDCAToVertexXY(2.4);
  fCuts->SetMaxDCAToVertexZ(3.2);
  fCuts->SetDCAToVertex2D(kTRUE);
  //reject fakes                                                                                                                                                            
  fCuts->SetMaxChi2PerClusterITS(36);
  fCuts->SetMaxChi2TPCConstrainedGlobal(36);

  fCuts->SetRequireSigmaToVertex(kFALSE);

  fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

  fCuts->SetEtaRange(-.9, .9);
  fCuts->SetPtRange(0.15, 1E+15);

  

  //TRACKS W/O SPD
  fCuts_noSPD->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
  fCuts_noSPD->SetMinNClustersTPC(70);
  fCuts_noSPD->SetMaxChi2PerClusterTPC(4);
  fCuts_noSPD->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1                                                                                                                                                             
  fCuts_noSPD->SetAcceptKinkDaughters(kFALSE);
  fCuts_noSPD->SetRequireTPCRefit(kTRUE);
  fCuts_noSPD->SetMaxFractionSharedTPCClusters(0.4);
  // ITS                                                                                                                                                                                                                                     
  fCuts_noSPD->SetRequireITSRefit(kTRUE);
  //accept secondaries                                                                                                                                                                                                                       
  fCuts_noSPD->SetMaxDCAToVertexXY(2.4);
  fCuts_noSPD->SetMaxDCAToVertexZ(3.2);
  fCuts_noSPD->SetDCAToVertex2D(kTRUE);
  //reject fakes                                                                                                                                                                                                                             
  fCuts_noSPD->SetMaxChi2PerClusterITS(36);
  fCuts_noSPD->SetMaxChi2TPCConstrainedGlobal(36);

  fCuts_noSPD->SetRequireSigmaToVertex(kFALSE);

  fCuts_noSPD->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);

  fCuts_noSPD->SetEtaRange(-.9, .9);
  fCuts_noSPD->SetPtRange(0.15, 1E+15);


  //TRACKS W/O ITS
  //usable for TPC multiplicity estimate?

  fCuts_noITS->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
  fCuts_noITS->SetMinNClustersTPC(70);
  fCuts_noITS->SetMaxChi2PerClusterTPC(4);
  fCuts_noITS->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1                                                                                                                                                     
                                                                                                                                                                                                                                             
  fCuts_noITS->SetAcceptKinkDaughters(kFALSE);
  fCuts_noITS->SetRequireTPCRefit(kTRUE);
  fCuts_noITS->SetMaxFractionSharedTPCClusters(0.4);
  // ITS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  fCuts_noITS->SetRequireITSRefit(kFALSE);
  //accept secondaries                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  fCuts_noITS->SetMaxDCAToVertexXY(2.4);
  fCuts_noITS->SetMaxDCAToVertexZ(3.2);
  fCuts_noITS->SetDCAToVertex2D(kTRUE);
  //reject fakes                                                                                                                                                                                                                            
                                                                                                                                                                                                                                            
  fCuts_noITS->SetMaxChi2PerClusterITS(36);
  fCuts_noITS->SetMaxChi2TPCConstrainedGlobal(36);

  fCuts_noITS->SetRequireSigmaToVertex(kFALSE);

  //fCuts_noSPD->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);

  fCuts_noITS->SetEtaRange(-.9, .9);
  fCuts_noITS->SetPtRange(0.15, 1E+15);
   

}



void AliAnalysisTaskFragmentationTriggered::UserExec(Option_t* /* option */)
{
  
  //  std::cout << "begin of userexec" << endl;
  //=========================== 
  //NO USEREXEC!!
  //  return;
  //NO USEREXEC!!
  //===========================
 
  if (fLeading>0) {std::cout << "fLeading was not set to 0 - Clearing it now..." << endl; 
    ClearLeading();
  }
   fHistCuts->Fill(0); // No Cut
  //fHistCuts->Fill("NoCut");

                                                                                                     
  //=======================================
  //Event Selection
  //=======================================
   
 AliInputEventHandler *inputHandler =					\
    (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();  
  

 AliMCEvent *mcEvent = MCEvent();
 AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*> (InputEvent());                                                                                                                                                           
 AliAODEvent *aodEvent = AODEvent();                                                                                                                                                                                                
 //event debug
 if (mcEvent) fHistEventDebug->Fill(0);
 if (esdEvent) fHistEventDebug->Fill(1);
 if (aodEvent) fHistEventDebug->Fill(2);


 // const AliESDVertex *vtx_esd = esdEvent->GetPrimaryVertex();                                                                                                                                                                                 
 /*
 if (!mcEvent) std::cout << "No MC Event found for this event" << endl;
 else std::cout << "Wow. MCEvent was found. Very MC. Much Simulation" << endl;

 if (esdEvent) std::cout << "esdEvent available" << endl;
 else std::cout << "esdEvent not available" << endl;

 if (aodEvent) std::cout << "aodEvent available" << endl;
 else std::cout << "aodEvent not available" << endl;

 //END EXECUTION HERE
 return;
 */

 // if (aodEvent) {

const AliAODVertex *vtx = aodEvent->GetPrimaryVertex();  

// if(inputHandler->IsEventSelected()) std::cout << "Event is selected" << endl;
// else std::cout << "event is not selected" << endl;
 
//TString triggerclasses = esdEvent->GetFiredTriggerClasses();

// std::cout << "Before trigger selection" << endl;

 if(!mcEvent){

 if (kmbtrigger){
   if (!( inputHandler->IsEventSelected() & AliVEvent::kINT7) )  
     return;
   // std::cout <<"The used trigger classes after kAnyInt selection are.." << triggerclasses << endl;
  fHistCuts->Fill(1); // after trigger selection 
 }
 
                                                                                      
                                                                                                                                                                                                                                           
  

 if (khjttrigger){
  
   if (!(inputHandler->IsEventSelected() & AliVEvent::kTRD) )
     return;
   fHistCuts->Fill(1);
   //   std::cout<<"The used trigger classes after kTRD selection are.."<< triggerclasses << endl;
  //select TRD Jet trigger:
  AliTRDTriggerAnalysis trdSelection;
  //trdSlection.CalcTriggers(vEvent);
  trdSelection.CalcTriggers(aodEvent);

  if (!( trdSelection.HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHJT)) )
    return;
  fHistCuts->Fill(2);
 }
 // std::cout << "after trigger selection, before pileup rejection" << endl;
 //=======================================
 //PILE UP REJECTION
 //========================================  
 
 //in bunch pileup
 if(aodEvent){
 if (fAnaUtils->IsPileUpEvent(aodEvent))
   return;}
 fHistCuts->Fill(3);

 //unclear if the following statement works for aods..
 //probably not because spd tracklets are not included in aods that are produced with PWGJE GSI filter

 //ATTENTION: SPD Cluster vs Tracklet BG works only after vAN-26-08....
 if(esdEvent){
 if (fAnaUtils->IsSPDClusterVsTrackletBG(esdEvent))
   return;}
  fHistCuts->Fill(4);
    
 }//close event selection for real data.. 

 AliAODHeader *aodhead =(AliAODHeader*) aodEvent->GetHeader();
 

 fHistCuts->Fill(5);//out of bunch pileup without vertex matching

    if( vtx->GetNContributors()<2 )
  return;
  fHistCuts->Fill(6); // Vertex exists 
   // fHistCuts->Fill("Vertex exists");
   
   fHistVertexZ->Fill(vtx->GetZ());
  
   fHistContribtoVtx->Fill(vtx->GetNContributors());

 
   if (!(vtx->GetNContributors() >= 3) )  
      return;
   fHistCuts->Fill(7); // at least 3 contributing tracks to the vertex
    // fHistCuts->Fill("Vertex Cuts");

   if (! (TMath::Abs(vtx->GetZ()) <=10 ))
   return;
  fHistCuts->Fill(8); // vertex should lie around |z|<=10cm  


 //fHistVertex->Fill(vtx->GetX(), vtx->GetY(), vtx->GetZ());
  fHistVertexXY->Fill(vtx->GetX(), vtx->GetY() );
  
  
  //=====================================================================================================
  //end of event selection
  //=====================================================================================================
  //


  //================
  //MCEVENT
  //================

  Double_t ptHard = 0.; //parton energy bins -> energy of particle                                                                                                               
  Double_t nTrials = 1; // trials for MC trigger weight for real data   


  if (mcEvent){
    //Get the generator headers  
      AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
      AliGenPythiaEventHeader*  pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader); 
     
      if(pythiaGenHeader){

	nTrials = pythiaGenHeader->Trials();
	ptHard  = pythiaGenHeader->GetPtHard();

	fHistMCTruthPtHard->Fill(ptHard);
	fHistMCTruthPtHardTrials->Fill(ptHard,nTrials);

      }

      // fh1Trials->Fill("#sum{ntrials}",fAvgTrials);
  



    AliStack *stackMC = mcEvent->Stack();
     if (!stackMC) {
       std::cout << "ERROR: MC Stack is not available!" << endl;
       }

     // std::cout << "The MC Event contains " << mcEvent->GetNumberOfTracks() << " tracks" << endl;
     // std::cout << "The MC Stack contains " << stackMC->GetNtrack() << " tracks" << endl;
     // std::cout << "The MC Event contains " << mcEvent->GetNumberOfPrimaries() << "primaries " << endl;
     // std::cout << "The MC Stack contains " << stackMC->GetNprimary() <<  "primaries" << endl;

    // Int_t NTracksMC = stackMC->GetEntriesFast();

    //std::cout << "The MC stack has " << NTracksMC << " entries" << endl;
    // fHistMCTruthNTracks->Fill(mcEvent->GetNumberOfTracks());
    Int_t MCTruthNTracks = 0;

    //Loop over all generated Particles
      for(int iPart = 1; iPart < (mcEvent->GetNumberOfTracks()); iPart++) {
   AliMCParticle *mcPart  = (AliMCParticle*)mcEvent->GetTrack(iPart);
   //ATTENTION: THIS LOOPS OVER ALL GENERATED PARTICLES. I HAVE TO MAKE SURE I USE ONLY THOSE THAT ARE OF INTEREST FOR ME...   
   //MC particles are all final state particles including particles from material interaction
   
   /*
   if (mcEvent->IsPhysicalPrimary(iPart)){
   std::cout << "The MC Particle " << iPart << " is a physical primary " <<  endl;
   std::cout << "It has the charge: " << mcPart->Charge() << endl;
   std::cout << "It has the mass " << mcPart->M() << endl;
   }
   */
   // std::cout << "The charge of this mc particle is: " << mcPart->Charge() << endl; 

   //mcPart->GetLabel();
   //for testing purposes fill arrays here..
   //all possible generated particles....
  

   //Use only primary charged tracks in eta acceptance with pt larger than 150MeV
   //(secondaries from weak decays and material interaction are neglected)

   if (TMath::Abs(mcPart->Eta()) < 0.9 && mcEvent->IsPhysicalPrimary(iPart) && TMath::Abs(mcPart->Charge())>0 && mcPart->Pt()>0.15){

   fMCTruthParticles->Add(mcPart);
   fHistMCTruthTrackEta->Fill(mcPart->Eta());
   fHistMCTruthTrackPt->Fill(mcPart->Pt());
   fHistMCTruthTrackPhi->Fill(mcPart->Phi());
   MCTruthNTracks++;
   // std::cout << "The Label of this MCTruth track is: " << mcPart->GetLabel() << endl;

   // std::cout << "The number of track references for this MCParticle is: " << mcPart->GetNumberOfTrackReferences() << endl;
     }//close track cuts

    
      }//close track loop

      fHistMCTruthNTracks->Fill(MCTruthNTracks);
      // std::cout << "The number of track references for this MC Event is: " << mcEvent->GetNumberOfTrackReferences() << endl;

      //========================================
      //INCLUSIVE JETS
      //=======================================

      //MC TRUTH JETS:
      Double_t RMC=0.4; //radius of jet                                                       
      Double_t jet_acceptanceMC = 0.9 - RMC;
      Double_t MCjetptleading =0;
      Double_t MCjet_ptthreshold = 5;
   
  TClonesArray *MCjetArray =
    dynamic_cast<TClonesArray*> (aodEvent->FindListObject("clustersAODMC2_ANTIKT04_B0_Filter00048_Cut00150_Skip00"));

  Int_t NJetsMC = MCjetArray->GetEntriesFast();

  for(Int_t ijetmc = 0; ijetmc < NJetsMC; ijetmc++){

    //NEED TO introduce the same cuts as for reconstructed jets...
    
    AliAODJet *mcjet = (AliAODJet*) MCjetArray->At(ijetmc);
    //eta acceptance cut
    if(TMath::Abs(mcjet->Eta()) > jet_acceptanceMC) continue;

    //std::cout << "The pt of this jet is: " << mcjet->Pt() << endl;
   
    //set the leading jet:
    SetLeading(mcjet->Pt());
    //std::cout << "The pt of the leading jet is atm is: " << MCjetptleading << endl;

 //plot all the jet properties
 fMCTruthAntiktHistPt->Fill(mcjet->Pt());
 fMCTruthAntiktHistAreavsPt->Fill(mcjet->EffectiveAreaCharged(), mcjet->Pt());
 fMCTruthAntiktHistAreaoverCirclevsPt->Fill(mcjet->EffectiveAreaCharged()/(TMath::Pi()*RMC*RMC), mcjet->Pt());


 //loop over tracks in jet to plot FFs
 
//could introduce jet pt cut here
//and fill all the labels from tracks in those jets in an array
//then for the cones: loop over this array and check if same label is found --> reject cone
//also: do some cross checks: no of rejected cones over total cones etc...

//INTRODUCE LOWER PT THRESHOLD HERE...

 if (mcjet->Pt() < MCjet_ptthreshold) continue;

 //add all jets in jet acceptance and above threshold to array:
 fMCjets->Add(mcjet);


 TRefArray* MCtrackarrayjet = mcjet->GetRefTracks();
 const Int_t MCntracksjet = MCtrackarrayjet->GetEntriesFast();

 //std::cout << "The number of tracks in this MC jet is: " << MCntracksjet << endl; 

 //loop over tracks associated with the jet                                                                                                                                                                                          
 for (Int_t itrack=0; itrack < MCntracksjet; itrack++){
   AliAODTrack *MCtrackjet = (AliAODTrack*) mcjet->GetTrack(itrack);

   // std::cout << "the pt of this track in the mc jet is: " << MCtrackjet->Pt() << endl;

   //write tracks in jets above some threshold to array
   fMCjetstracks->Add(MCtrackjet);
   //maybe also write number of tracks in jet into array to have a connection
    Double_t z = MCtrackjet->Pt()/mcjet->Pt();
    Double_t xi = TMath::Log(1/z);

   // fAntiktHistz->Fill(z);
   //fAntiktHistzpt->Fill(z, jet->Pt() );
   //fAntiktHistxi->Fill(xi);
   fMCTruthAntiktHistxipt->Fill(xi, mcjet->Pt() );

   //need perp cone for inclusive jets? 



 }//close loop over tracks in jet



 //in jet loop:
 //now: loop over all hybrid tracks in event and find tracks that lie in perp cones, cone around jet axis..



  }//close loop over mc jets


  MCjetptleading = GetLeading(); 
  ClearLeading(); // IMPORTANT: RESETS THE COUNTER

  fMCTruthAntiktHistPtLeading->Fill(MCjetptleading);

  //std:: cout << "The leading jet pt is: " << MCjetptleading << endl;

  //now: loop over jets to find leading jet and plot all the interesting properties.... 


  






  //================================
  //RANDOM CONE
  //================================
  //One random cone per jet event. Independent of all jets.

  //flag event if jet above lower threshold is present
  //if flag = true find tracks in random cone

  if (MCjetptleading > MCjet_ptthreshold){

    //throw random cones until there is no overlap with signal jets...
    //then use this cone to plot all the interesting properties
    //idea: use do{bla} while{ntracksoverlap>0} loop

    //Better rename these variables to sth related with MCRC..

    
    Double_t MCRCPhi=0;
    Double_t MCRCEta=0;
    Double_t MCRCdeta = -9;
    Double_t MCRCdphi = -9;
    Double_t MCRCpt = 0;
    Int_t MCRCntracks = 0;
    Int_t MCRCntracksoverlap = 0;
    Int_t MCRCntry = 0; //how many random cones until no overlap
    Int_t MCRCntrybreak = 50;
    Bool_t breakMCRC = kFALSE;

    //Find tracks in random cone and check overlap with signal jets. do until no overlap.

    do {

      //set number of overlapping tracks to zero
      MCRCntracksoverlap=0;

      //clear random cone properties: (so that only for no overlap these variables are used)
      fMCRC->Clear();
      MCRCpt=0;
      MCRCntracks=0;
      MCRCntry++;

      //need abort criterion: if there are too much signal jets (high multiplicity events), the loop will be endless...                     
      if (MCRCntry > MCRCntrybreak) {
	//break do while loop                                                                                                               
	breakMCRC = kTRUE;
	MCRCntry = 0; // set to value visible in histogram and unique                                                                         
	break;
      }


  //get random distributed numbers in (0, 2pi)                                                                                                                                                                                                
  MCRCPhi = fRandomGen->Rndm()*2*TMath::Pi();

  //get random distributed numbers in ( -jetacceptance, jetacceptance)                     
  //stretch interval by 2*jetacceptance first, then subtract jetacceptance               
                                                                                                                                                                     
  MCRCEta = fRandomGen->Rndm()*2*jet_acceptanceMC - jet_acceptanceMC;

 
  for (Int_t iMCPart=0; iMCPart < fMCTruthParticles->GetEntriesFast(); iMCPart++ ) {
    //use MC Particles from array that have cuts                                                                                                                                                                                              
    AliMCParticle *mcPart  = (AliMCParticle*) fMCTruthParticles->At(iMCPart);
    //random cone with veto for cone, that overlaps with jet for wich perp cones are built                                                                                                                                                   
    //Area of random cones 
    MCRCdeta = mcPart->Eta() - MCRCEta;
    MCRCdphi = mcPart->Phi() - MCRCPhi;

    MCRCdeta = TVector2::Phi_mpi_pi(MCRCdeta);
    MCRCdphi = TVector2::Phi_mpi_pi(MCRCdphi);



    if (MCRCdeta*MCRCdeta + MCRCdphi*MCRCdphi < RMC*RMC){
      //Add tracks to Random Cone Array                                                                                                                                                                                                     

      //check overlap with jets: for every MC particle in RC loop over tracks in signal jets
      for (Int_t iMCtrackjet=0; iMCtrackjet < fMCjetstracks->GetEntriesFast(); iMCtrackjet++){  
	AliMCParticle *mctrackinsignaljet = (AliMCParticle*) fMCjetstracks->At(iMCtrackjet);
	if (TMath::Abs(mctrackinsignaljet->GetLabel()) == TMath::Abs(mcPart->GetLabel())) MCRCntracksoverlap++;

      }//close loop over all tracks in signal jets

      //if running this in a loop the values should be set to zero in every cicle
      fMCRC->Add(mcPart);
      MCRCpt+=mcPart->Pt();
      MCRCntracks++;
  
     }//close if track in RC
  }//close loop over mc particles in array

    }//close "do"
    while(MCRCntracksoverlap>0);


  //Some RC properties:
    if (!breakMCRC) {
  fMCTruthRCntry->Fill(MCRCntry);
  fMCTruthRCnTracksOverlap->Fill(MCRCntracksoverlap);
  fMCTruthAntiktHistPhiRC->Fill(MCRCPhi);                                                                                                                                        
  fMCTruthAntiktHistEtaRC->Fill(MCRCEta);  
  fMCTruthHistRCPt->Fill(MCRCpt);
  fMCTruthHistRCnTracks->Fill(MCRCntracks);
    }


  //CHeck no of overlapping tracks:

  if (MCRCntracksoverlap==0){
    //Plot all the RC properties
    //    fMCTruthHistRCPt->Fill(MCRCpt);
    //fMCTruthHistRCnTracks->Fill(MCRCntracks);

    //loop over particles in random cone to plot fragmentation variables. Different possibilities: inclusive cones or leading cones    

    //here: Using leading jet pt

     for (Int_t itrackRC=0; itrackRC < fMCRC->GetEntriesFast(); itrackRC++){

       AliMCParticle *RCPart  = (AliMCParticle*) fMCRC->At(itrackRC);

       //General Properties: eta-phi....
       fMCTruthHistRCTrackEta->Fill(RCPart->Eta());
       fMCTruthHistRCTrackPhi->Fill(RCPart->Phi());

       //Contribution to FF of leading jet: 
       Double_t z = RCPart->Pt()/MCjetptleading;
       Double_t xi = TMath::Log(1/z);

       fMCTruthHistRCxiptLeading->Fill(xi, MCjetptleading);
    
       //plot number of jets above threshold and in acceptance
       fMCTruthnJets->Fill(fMCjets->GetEntriesFast());

       //Contribution to FF of inclusive jets
       for(Int_t iMCjets=0; iMCjets<  fMCjets->GetEntriesFast(); iMCjets++ ) {

	 AliAODJet *mcjet = (AliAODJet*) fMCjets->At(iMCjets);

	 Double_t z_incl = RCPart->Pt()/mcjet->Pt();
	 Double_t xi_incl = TMath::Log(1/z_incl);

	 //FIll histograms..
	 fMCTruthHistRCxipt->Fill(xi_incl, mcjet->Pt());

       }//close loop over inclusive jets above threshold and in acceptance

    }//close track loop
  
  }//close if no overlap

  //if overlap...
  //should not be the case ... use these histograms as check if sth went wrong...
  else {
    //Plot RC properties for discarded cones
    fMCTruthHistRCWithOverlapPt->Fill(MCRCpt);
    fMCTruthHistRCWithOverlapnTracks->Fill(MCRCntracks);

    //here: Using leading jet pt                                                                                                                                                   

    for (Int_t itrackRC=0; itrackRC < fMCRC->GetEntriesFast(); itrackRC++){

      AliMCParticle *RCPart  = (AliMCParticle*) fMCRC->At(itrackRC);
      //General Properties:

      fMCTruthHistRCTrackEtaDiscarded->Fill(RCPart->Eta());
      fMCTruthHistRCTrackPhiDiscarded->Fill(RCPart->Phi());


      //Contribution to FF of leading jet:                                                                                                                                        
      Double_t z = RCPart->Pt()/MCjetptleading;
      Double_t xi = TMath::Log(1/z);

      fMCTruthHistRCxiptLeadingDiscarded->Fill(xi, MCjetptleading);
    }//close track loop
  }//close if overlap

  //    std::cout << "The pt in the random cone is: " << RCpt << endl;
  //  std::cout << "The no of trakcs in RC is: " << RCntracks << endl;


  Double_t rhoMC=-1;

  //========================
  //LEADING JETS
  //=======================  

  //loop over all jets and find leading jet

  for(Int_t ijetmc = 0; ijetmc < NJetsMC; ijetmc++){

   //NEED TO introduce the same cuts as for reconstructed jets...                                                                                                              
    AliAODJet *mcjet = (AliAODJet*) MCjetArray->At(ijetmc);
    //eta acceptance cut                                                                                                                                                         
    if(TMath::Abs(mcjet->Eta()) > jet_acceptanceMC) continue;
    //jet pt threshold
    if (mcjet->Pt() < MCjet_ptthreshold) continue;

    //if jet pt = leading jet pt then its the leading jet...

    if (TMath::Abs(MCjetptleading - mcjet->Pt()) < 0.001) {
      //plot all the interesting leading jet properties...

      fMCTruthAntiktHistEtaLeading->Fill(mcjet->Eta());
      fMCTruthAntiktHistPhiLeading->Fill(mcjet->Phi());


      //FFs etc...

      //loop over tracks associated with the leading jet   
      TRefArray* MCtrackarrayjet = mcjet->GetRefTracks();
      const Int_t MCntracksjet = MCtrackarrayjet->GetEntriesFast();
 
      for (Int_t itrack=0; itrack < MCntracksjet; itrack++){
	AliAODTrack *MCtrackjet = (AliAODTrack*) mcjet->GetTrack(itrack);

	// std::cout << "the pt of this track in the mc jet is: " << MCtrackjet->Pt() << endl;                                                                                                                                                                                                              
	                                                                                                      
	Double_t z = MCtrackjet->Pt()/mcjet->Pt();
	Double_t xi = TMath::Log(1/z);

                                                                                                                                                                    
	fMCTruthAntiktHistxiptLeading->Fill(xi, mcjet->Pt() );                                                                                                                                               

      }//close loop over tracks in jet   

      //define perp cones... 
      //loop over all tracks in event and find tracks in perp cones using the position of the leading jet...

      Double_t detaPC = 0;
      Double_t dphiPCP = 0;
      Double_t dphiPCM = 0;
      Double_t ptPCP = 0;
      Int_t ntracksPCP = 0;
      Double_t ptPCM = 0;
      Int_t ntracksPCM = 0;


      for (Int_t iMCPart=0; iMCPart < fMCTruthParticles->GetEntriesFast(); iMCPart++ ) {
	//use MC Particles from array that have cuts                            
                                                                                                                                                                                 
	AliMCParticle *mcPart  = (AliMCParticle*) fMCTruthParticles->At(iMCPart);
	                                                                                                                                                                          
	//Area of perp cones:                                                                                                                                        
	detaPC = mcPart->Eta() - mcjet->Eta();
	dphiPCP = mcPart->Phi() - mcjet->Phi() - TMath::Pi()/2;
	dphiPCM = mcPart->Phi() - mcjet->Phi() + TMath::Pi()/2;

	detaPC = TVector2::Phi_mpi_pi(detaPC);
	dphiPCP = TVector2::Phi_mpi_pi(dphiPCP);
	dphiPCM = TVector2::Phi_mpi_pi(dphiPCM);

	//perp cone plus:
	if (detaPC*detaPC + dphiPCP*dphiPCP < RMC*RMC) {
	  //add track to array
	  fMCtracks_pcp->Add(mcPart);
	  ntracksPCP++;
	  ptPCP+=mcPart->Pt();
	}

	//perp cone minus:
	if (detaPC*detaPC + dphiPCM*dphiPCM < RMC*RMC){
	  //add track to array
	  fMCtracks_pcm->Add(mcPart);
	  ntracksPCM++;
	  ptPCM+=mcPart->Pt();
        }


      }//close loop over all mc particles


      //============================================================                             
      //Event by event UE subtraction                                                                                                          
      //============================================================                                                                                                       

      //Better: For inclusive jets: Define RHO and then subtract RHO from inclusive jets in loop                                                                           
      rhoMC = (ptPCP+ptPCM)/(2*TMath::Pi()*RMC*RMC);

      fMCTruthHistRho->Fill(rhoMC);
      fMCTruthHistRhovsLJPt->Fill(rhoMC, MCjetptleading);




      //now loop over particles in perp cone array and plot ff's
      //PCP using leading jet pt:

      for (Int_t itrackPCP=0; itrackPCP < fMCtracks_pcp->GetEntriesFast(); itrackPCP++) {
        AliMCParticle *mcPart  = (AliMCParticle*) fMCtracks_pcp->At(itrackPCP);

	Double_t z = mcPart->Pt()/mcjet->Pt();
        Double_t xi = TMath::Log(1/z);


        fMCTruthAntiktHistxiptLeading_pcp->Fill(xi, mcjet->Pt() );       

      }// close loop over tracks in perp cone plus

      //need also loop over tracks in pcm?

    } // close if leading jet


    //======================================       
    //EVENT BY EVENT SUBTRACTION                                                                                                               
    //======================================                                                                                                                                   
    //(we are still in the loop over inclusive jets)                                                                                                                           
    //now: subtract RHO from jet pt and plot observables of interest with subtracted jet pt 
    Double_t ptsubMC = mcjet->Pt() - rhoMC*mcjet->EffectiveAreaCharged();
    //Now: plot all jet observables with ptsub                                                                                                                                 

    fMCTruthHistRhoSubPt->Fill(ptsubMC);
    //eta, phi, area                                                                                                                                                           
    fMCTruthHistRhoSubAreaOverCirclevsPt->Fill(mcjet->EffectiveAreaCharged()/(TMath::Pi()*RMC*RMC), ptsubMC);

    //need loop over tracks in signal jet?
    TRefArray* mctrackarrayjet = mcjet->GetRefTracks();
    const Int_t mcntracksjet = mctrackarrayjet->GetEntriesFast();
    //xi jet                                                                                                                                
    for (Int_t itrack=0; itrack < mcntracksjet; itrack++){
      AliAODTrack *mctrackjet = (AliAODTrack*) mcjet->GetTrack(itrack);

      //eta phi of tracks in jets                                                                                                           
      fMCTruthHistRhoSubPhiTracksJetvsPt->Fill(mctrackjet->Phi(), ptsubMC);
      fMCTruthHistRhoSubEtaTracksJetvsPt->Fill(mctrackjet->Eta(), ptsubMC);


      Double_t z = mctrackjet->Pt()/ptsubMC;
      Double_t xi = TMath::Log(1/z);

      fMCTruthHistRhoSubXiPt->Fill(xi, ptsubMC );

    }//close loop over tracks in jet      



    //now: in jet loop: for a given jet: loop over tracks in perp cones and plot ff's using the inclusive jet pt...


    for (Int_t itrackPCP=0; itrackPCP < fMCtracks_pcp->GetEntriesFast(); itrackPCP++) {
      AliMCParticle *mcPart  = (AliMCParticle*) fMCtracks_pcp->At(itrackPCP);
      //unsubtracted:
      Double_t z = mcPart->Pt()/mcjet->Pt();
      Double_t xi = TMath::Log(1/z);

      fMCTruthAntiktHistxipt_pcp->Fill(xi, mcjet->Pt() );
      fMCTruthAntiktHistxipt_pc->Fill(xi, mcjet->Pt() );

      //subtracted:
      Double_t zsub = mcPart->Pt()/ptsubMC;
      Double_t xisub = TMath::Log(1/zsub);

      fMCTruthHistRhoSubXiPCPt->Fill(xisub, ptsubMC);

      //eta phi of tracks in PC                                                                                                              
      fMCTruthHistRhoSubPhiTracksPCvsPt->Fill(mcPart->Phi(), ptsubMC);
      fMCTruthHistRhoSubEtaTracksPCvsPt->Fill(mcPart->Eta(), ptsubMC);





    }//close loop over pcp tracks

      
      //need also perp cone minus...

    for (Int_t itrackPCM=0; itrackPCM < fMCtracks_pcm->GetEntriesFast(); itrackPCM++) {
      AliMCParticle *mcPart  = (AliMCParticle*) fMCtracks_pcm->At(itrackPCM);
      //unsubtracted:
      Double_t z = mcPart->Pt()/mcjet->Pt();
      Double_t xi = TMath::Log(1/z);

      fMCTruthAntiktHistxipt_pcm->Fill(xi, mcjet->Pt() );
      fMCTruthAntiktHistxipt_pc->Fill(xi, mcjet->Pt() );

      //subtracted:
      Double_t zsub = mcPart->Pt()/ptsubMC;
      Double_t xisub = TMath::Log(1/zsub);

      fMCTruthHistRhoSubXiPCPt->Fill(xisub, ptsubMC);

      //eta phi of tracks in PC                                                                                                             
      fMCTruthHistRhoSubPhiTracksPCvsPt->Fill(mcPart->Phi(), ptsubMC);
      fMCTruthHistRhoSubEtaTracksPCvsPt->Fill(mcPart->Eta(), ptsubMC);


    }//close loop over pcp tracks 



          }// close jet loop

     }//close if leading jet above threshold 
 
  }//close MC Event




  //====================================
  //ESDEVENT 
  //===================================


  if (esdEvent) {

    Float_t xy_global=-19;
    Float_t z_global=-19;
    Float_t xy_noSPD=-19;
    Float_t z_noSPD=-19;
    Float_t xy_noITS=-19;
    Float_t z_noITS=-19;
    Int_t tpcmult=0;

      // loop over tracks, ...
     const Int_t nTracks = esdEvent->GetNumberOfTracks();
   
     //std::cout << "The esdEvent contains " << nTracks << " tracks" << endl;
     //   const  AliMultiplicity *multipl= esdEvent->GetMultiplicity();
     // std::cout <<"no of tracklets from esdevent is: " <<  multipl->GetNumberOfTracklets() << endl;

   //As used in AnalysisUtils:                                                                                                                                                                                                              
   Int_t nClustersLayer0 = esdEvent->GetNumberOfITSClusters(0);
   Int_t nClustersLayer1 = esdEvent->GetNumberOfITSClusters(1);
   Int_t nTracklets      = esdEvent->GetMultiplicity()->GetNumberOfTracklets();
   Int_t sumCluster = nClustersLayer0 + nClustersLayer1;

   //std::cout << "(ESD) as in anautils: "<< nClustersLayer0<< " , " << nClustersLayer1 <<" , " << nTracklets << endl;

   fHistSPDNclsvstracklets->Fill(nTracklets, sumCluster);

    
   for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
     
      AliESDtrack *trk = esdEvent->GetTrack(iTrack);
      fHistTrackCuts->Fill(0);
      fHistTrackCutsDep->Fill(0);       
    
  
     if (fCuts->AcceptTrack(trk) ) { 
      
       trk->GetImpactParameters(xy_global,z_global);
       fHistImpactParameterglobalxy->Fill(xy_global);
       fHistImpactParameterglobalz->Fill(z_global);
      
         
     }//close track cuts

     if (fCuts_noSPD->AcceptTrack(trk) ) {
       trk->GetImpactParameters(xy_noSPD,z_noSPD);
       fHistImpactParameternoSPDxy->Fill(xy_noSPD);
       fHistImpactParameternoSPDz->Fill(z_noSPD);
      

}

     if (fCuts_noITS->AcceptTrack(trk) && ( !(trk->GetStatus() & AliESDtrack::kITSrefit))  ) {
       trk->GetImpactParameters(xy_noITS,z_noITS);
       fHistImpactParameternoITSxy->Fill(xy_noITS);
       fHistImpactParameternoITSz->Fill(z_noITS);    
      
 }

     //TPC Multiplicity:
     if (fCuts_noITS->AcceptTrack(trk)) ++tpcmult;

   
    }//close track loop

   // std::cout << "The tpc multiplicity is: " << tpcmult << endl;

  } // close if esdevent
  //========================================================================================


  //=====================================================
  //AODEVENT
  //=====================================================


  //number of hybrid tracks in event
  Double_t ntrackshybrid = 0;
  //  Int_t nSPDhits_trk = 0;
  Int_t nSPDhits_l0 = 0;
  Int_t nSPDhits_l1 = 0;
  Int_t nSPDhits_evt=0;

  // get aod data and loop over jets
  if (aodEvent) {
    // cout << "aodEvent exists" << endl;
    //look out for Radius of jet: ANTIKTR_
   
 //  TClonesArray *jetArray =
    //      dynamic_cast<TClonesArray*> (aodEvent->FindListObject("clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00"));


    //     TClonesArray *jetArray =
    // dynamic_cast<TClonesArray*> (aodEvent->FindListObject("clustersAOD_ANTIKT02_B0_Filter00048_Cut00150_Skip00"));


    //   TClonesArray *jetArray =
    // dynamic_cast<TClonesArray*> (aodEvent->FindListObject("clustersAOD_ANTIKT03_B0_Filter00048_Cut00150_Skip00"));

     
    TClonesArray *jetArray =
   dynamic_cast<TClonesArray*> (aodEvent->FindListObject("clustersAOD_ANTIKT04_B0_Filter00048_Cut00150_Skip00"));
    
    
   //convert clusteraod.... to a pointer to a tclonesarray
  
    /*        if (!aodEvent) {
     cout<< "aodEvent doesnt exist" << endl; 
     }  

     if (!jetArray){
     cout << "jetArray doesnt exist" << endl;
     aodEvent->Print();
     }
    */
    //As used in AnalysisUtils:    
    // Int_t nClustersLayer0 = aodEvent->GetNumberOfITSClusters(0);
    // Int_t nClustersLayer1 = aodEvent->GetNumberOfITSClusters(1);
    // Int_t nTracklets      = aodEvent->GetMultiplicity()->GetNumberOfTracklets();

    // std::cout << "as in anautils: "<< nClustersLayer0<< " , " << nClustersLayer1 <<" , " << nTracklets << endl;


    // AliAODTracklets *tracklets = aodEvent->GetTracklets();

    const Int_t nTracks = aodEvent->GetNumberOfTracks();


    //    Int_t aodiddimension = aodEvent->GetNumberOfTracks();
    // std::vector<int> track_IDs(aodiddimension);

    // std::cout << nTracks << endl;

    for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
      AliAODTrack *trk = (AliAODTrack*) aodEvent->GetTrack(iTrack);

    // Int_t ntracklets = tracklets->GetNumberOfTracklets();
    //    std::cout << "the number of spd tracklets is: " << ntracklets << endl;

    //std::cout << "haspointonitslayer(0) returns:  " << trk->HasPointOnITSLayer(0) << endl;
    // std::cout << "haspointonitslayer(1) returns:  " << trk->HasPointOnITSLayer(1) << endl;
    //    nSPDhits_trk=trk->HasPointOnITSLayer(0)+trk->HasPointOnITSLayer(1);    
    //std::cout << "the sum of spd hits is: " << nSPDhits << endl;
    if(trk->HasPointOnITSLayer(0)){
     ++nSPDhits_l0;
     fHistITShitsinlayer->Fill(1);}
   
    if(trk->HasPointOnITSLayer(1)){
     ++nSPDhits_l1;
     fHistITShitsinlayer->Fill(2);}

    if(trk->HasPointOnITSLayer(2)) fHistITShitsinlayer->Fill(3);
    if(trk->HasPointOnITSLayer(3)) fHistITShitsinlayer->Fill(4);
    if(trk->HasPointOnITSLayer(4)) fHistITShitsinlayer->Fill(5);
    if(trk->HasPointOnITSLayer(5)) fHistITShitsinlayer->Fill(6);

    fHistTrackCuts->Fill(0);                                                                                                                                                                                                               
    fHistTrackCutsDep->Fill(0);                                                                                                                                                                                                            
   
                                                                                                                       
    if (TMath::Abs(trk->Eta())<=0.9 ) {fHistTrackCuts->Fill(1);}                                                                                                                                                                           
    if (trk->Pt() > 0.15 && trk->Pt() < 1E+15) {fHistTrackCuts->Fill(2);}                                                                                                                                                                  
    if (trk->GetTPCNcls() > 70 ) {fHistTrackCuts->Fill(3);}    

    if (TMath::Abs(trk->Eta())<=0.9) {fHistTrackCutsDep->Fill(1);}                                                                                                                                                                         
    if (TMath::Abs(trk->Eta())<=0.9 && trk->Pt() > 0.15 && trk->Pt() < 1E+15) {fHistTrackCutsDep->Fill(2);}                                                                                                                                
    if (TMath::Abs(trk->Eta())<=0.9 && trk->Pt() > 0.15 && trk->Pt() < 1E+15 && trk->GetTPCNcls() > 70) {fHistTrackCutsDep->Fill(3);}                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                             
    if(!trk) {cout <<" ERROR - pointer to track is NULL " << endl; }                                                                                                                                                                       
    if(!fCuts) {cout << "ERROR - pointer to trackcuts is NULL" << endl;}                                                                                                                                                                                                                                                                                                                                           
      
    if (trk->IsHybridGlobalConstrainedGlobal() && (trk->GetStatus() & AliAODTrack::kITSrefit)) {
      // use hybrid track cuts                                                                                                                                                                                                              

      if (mcEvent){    

	AliMCParticle *mctruthrectrack =(AliMCParticle*) mcEvent->GetTrack(TMath::Abs(trk->GetLabel() ));

	//	std::cout << "reco label: " << TMath::Abs(trk->GetLabel()) << "truth label: " << TMath::Abs(mctruthrectrack->Label()) << endl;
	//is primary of interest?
	//	if (){
	if (TMath::Abs(mctruthrectrack->Eta()) < 0.9 && mcEvent->IsPhysicalPrimary(TMath::Abs(trk->GetLabel() )) && TMath::Abs(mctruthrectrack->Charge())>0 && mctruthrectrack->Pt()>0.15){
	fMCTruthRecTrackPtGen->Fill(mctruthrectrack->Pt());

	}
	  //secondaries:
	else if(!mcEvent->IsPhysicalPrimary(TMath::Abs(trk->GetLabel() ))){
	fMCTruthRecTrackSecPtGen->Fill(mctruthrectrack->Pt());
	fMCTruthRecTrackSecPtRec->Fill(trk->Pt());
	}


       //Get MCTruth track with same label
      // fMCTruthLabel->GetEntriesFast();
      for(Int_t ilabelentry=0; ilabelentry < fMCTruthParticles->GetEntriesFast(); ilabelentry++){


	AliMCParticle *mcparticle = (AliMCParticle*) fMCTruthParticles->At(ilabelentry);
	//ATTENTION: LABEL NOs in AOD or ESD can be negative while in MC they cant
	if(TMath::Abs(trk->GetLabel()) == TMath::Abs(mcparticle->Label()) ){

	  fHistMCTruthTrackPtRatio->Fill(trk->Pt()/mcparticle->Pt(), mcparticle->Pt());	 
	  //plot also for the different classes of tracks: with and without SPD hit

	  if (trk->HasPointOnITSLayer(0) || trk->HasPointOnITSLayer(1)) {
	    //plot for tracks with spd hit
	    fHistMCTruthTrackPtRatiowSPD->Fill(trk->Pt()/mcparticle->Pt(), mcparticle->Pt());

	  }//close if spd hit
          else {
	    //plot for tracks without spd hit
	    fHistMCTruthTrackPtRationoSPD->Fill(trk->Pt()/mcparticle->Pt(), mcparticle->Pt());

         }	 

         
	}//close if same label
      }//close loop over mctruthlabelarray
      }//close if mcevent

      if (TMath::Abs(trk->Eta())<0.9) {++ntrackshybrid;}
      fHistTrackCuts->Fill(9);
      fHistTrackCutsDep->Fill(9);
     fHistPt->Fill(trk->Pt());                                                                                                                                                                                                              
    
      if (trk->Pt()>3){
       fHistEta->Fill(trk->Eta());                                                                                                                                                      fHistPhi->Fill(trk->Phi());                                                                                                                                                       fHistEtaPhi->Fill(trk->Eta(), trk->Phi());                                                                                                                                           }  

     

      if(!(trk->IsGlobalConstrained())) {fHistPhiHybridUnconstrained->Fill(trk->Phi());
	                  
      } //close track not globally constrained

      if(trk->IsGlobalConstrained()){
	if(trk->GetStatus() & AliAODTrack::kITSrefit) {fHistPhiHybridConstrainedwITS->Fill(trk->Phi());
	 
	}//close if its refitted

        else if ( !(trk->GetStatus() & AliAODTrack::kITSrefit) ){fHistPhiHybridConstrainedwoITS->Fill(trk->Phi());
	 
	} //close no its refit tracks
	}//close is global constrained

   //Plot no of ITS hits                                                                                                                                                                                                                  
     fHistITSNcls->Fill(trk->GetITSNcls());                                                                                                                                                                                                   
   //plot no of TPC hits                                                                                                                                                                                                                  
    fHistTPCNcls->Fill(trk->GetTPCNcls());                                                                                                                                                                                                   //fHistTPCchi2->Fill(trk->GetTPCchi2());                                                                                                                                                                                                 
   // fHistTPCchi2overtpchits->Fill(trk->GetTPCchi2() / trk->GetTPCNcls() );                                                                                                                                                      
   //fHistImpactParameterxy->Fill(dcat);                                                                                                                                                                                                    
   // fHistImpactParameterz->Fill(dcaz);       

    }//close hybrid track cuts
    
    }//close track loop

    fHistnTracks->Fill(ntrackshybrid);
    nSPDhits_evt = nSPDhits_l0 + nSPDhits_l1;    
    
    //already defined upstairs:
    // AliAODHeader *aodhead = aodEvent->GetHeader();
    Int_t multipl_aodhead = aodhead->GetRefMultiplicityComb08();
    //ntrackshybrid: no of hybrid tracks in |eta|<0.9
    fHistMulthybridvsglobal->Fill(multipl_aodhead, ntrackshybrid);
    //  std::cout << "multiplicity from aodheader is: " << multipl_aodhead << endl;

    // std::cout << "first layer: " << nSPDhits_l0 << "second layer: " << nSPDhits_l1 << endl;
    //fHistSPDNclsvstracklets->Fill(aodEvent->GetTracklets(), nSPDhits_evt);

    
    //    AliAODTracklets *tracklets = aodEvent->GetTracklets();
    // Int_t ntracklets = tracklets->GetNumberOfTracklets();
    // std::cout << "the number of spd tracklets is: " << ntracklets << endl;

    //use different R for cone method for some crosscheck
            Double_t R_LJCM = 0.4;
    //       Double_t R=0.2; //radius of jet       
    //      Double_t R=0.3; //radius of jet   
            Double_t R=0.4; //radius of jet    
	    Double_t jet_acceptance = 0.9 - R;
	    Double_t jet_ptthreshold = 5;

    Int_t njetsevent=0;

     if (jetArray) {
       // cout << "jetArray exists" << endl;
      const Int_t nJets = jetArray->GetEntriesFast();
      
      // fHistnJets->Fill(nJets);
    
      Double_t ptleading=0; 
      // loop over jets in event
        for (Int_t iJet = 0; iJet < nJets; ++iJet) {
	AliAODJet *jet = (AliAODJet*) jetArray->At(iJet);

        //eta acceptance cut:
        if (TMath::Abs(jet->Eta()) > jet_acceptance ) continue;
	
	//Pt Spectrum (only eta acceptance cut)                                                                                                                                   
        fAntiktHistPt->Fill(jet->Pt());

	//pt cut
        if (jet->Pt() < jet_ptthreshold) continue;
             ++njetsevent;
	     fJets->Add(jet);
	

        //regions without trd:                                                                                                                                                                                                              
	if ((jet->Phi()>4*TMath::Pi()/9 && jet->Phi()< 6*TMath::Pi()/9) || (jet->Phi()>12*TMath::Pi()/9 && jet->Phi()< 15*TMath::Pi()/9)  ) {fAntiktHistPtwoTRD->Fill(jet->Pt());}
	//regions with trd:                                                                                                                                                                                                                 
	else {fAntiktHistPtwTRD->Fill(jet->Pt());}

        
        fAntiktHistPhi->Fill(jet->Phi());
        fAntiktHistEta->Fill(jet->Eta());
	fAntiktHistEtavsPt->Fill(jet->Eta(), jet->Pt());
        fAntiktHistEtaPhi->Fill(jet->Eta(), jet->Phi());
	//inclusive jet area:
	fAntiktHistAreavsPt->Fill(jet->EffectiveAreaCharged(), jet->Pt());
	fAntiktHistAreaoverCirclevsPt->Fill(jet->EffectiveAreaCharged()/(TMath::Pi()*R*R), jet->Pt());


        //get jet with largest pt in event:
        if (jet->Pt()>ptleading ){
	  ptleading = jet->Pt();
	  //for some events (with no jet in event) comparison is not true and thus ptleading wont be equal to jet->pt... 
          }
    
        //now: get track references from jet and get track properties to plot inclusive ff etc..
   
         TRefArray* trackarrayjet = jet->GetRefTracks();
	 const Int_t ntracksjet = trackarrayjet->GetEntriesFast();
       
	   fAntiktHistnTracks->Fill(ntracksjet);	 

         //loop over tracks associated with the jet
         for (Int_t itrack=0; itrack < ntracksjet; itrack++){
	   AliAODTrack *trackjet = (AliAODTrack*) jet->GetTrack(itrack);
	 
	   //fill tracks in signal jet into array:
	   fJetsTracks->Add(trackjet);
 
         Double_t z = trackjet->Pt()/jet->Pt();
         Double_t xi = TMath::Log(1/z);
   
         fAntiktHistz->Fill(z);
         fAntiktHistzpt->Fill(z, jet->Pt() );
         fAntiktHistxi->Fill(xi); 
         fAntiktHistxipt->Fill(xi, jet->Pt() );        
	 //plot fragmentation functions for different regions in phi (with and without TRD)
	 //regions without trd:
	 if ((jet->Phi()>4*TMath::Pi()/9 && jet->Phi()< 6*TMath::Pi()/9) || (jet->Phi()>12*TMath::Pi()/9 && jet->Phi()< 15*TMath::Pi()/9)  ) {fAntiktHistxiptwoTRD->Fill(xi, jet->Pt());}
	 //regions with trd:
	 else {fAntiktHistxiptwTRD->Fill(xi, jet->Pt());}


          }//close loop over tracks in jet

       } // close jet loop


	//	std::cout << "number of jets in array: " << fJets->GetEntriesFast() << " number of tracks in signal jets: " << fJetsTracks->GetEntriesFast() << endl;

	fHistnJets->Fill(njetsevent);
	fAntiktHistPtLeading->Fill(ptleading);
	fAntiktHistPtLeadingvsHybridMult->Fill(ptleading, ntrackshybrid);	      


	//================================================
	//Random Cone
	//================================================

        Double_t rhoRC = -1; //initialise                                                                                                                                          
	//only introduce random cone if a jet above a certain threshold is present in the event
	if (ptleading > jet_ptthreshold) {

	Double_t RCPhi=0;
	Double_t RCEta=0;
	Double_t RCdeta = -9;
	Double_t RCdphi = -9;
	Double_t RCpt = 0;
	Int_t RCntracks = 0;
	Int_t RCntracksoverlap = 0;
	Int_t RCntry = 0; //how many random cones until no overlap 
	Int_t RCntrybreak = 50;
	Bool_t breakRC = kFALSE;

	do {

	  //set number of overlapping tracks to zero        
	  RCntracksoverlap=0;

	  //clear random cone properties: (so that only for no overlap these variables are used) 
	  fRC->Clear();
	  RCpt=0;
	  RCntracks=0;
	  RCntry++;

	  //need abort criterion: if there are too much signal jets (high multiplicity events), the loop will be endless...
	  if (RCntry > RCntrybreak) {
	    //break do while loop
	    breakRC = kTRUE;
	    RCntry = 0; // set to value visible in histogram and unique
	    break;
         }

	  //get random distributed numbers in (0, 2pi)  
	  RCPhi = fRandomGen->Rndm()*2*TMath::Pi();

	  //get random distributed numbers in ( -jetacceptance, jetacceptance) 
	  //stretch interval by 2*jetacceptance first, then subtract jetacceptance                                                                                                         
	  RCEta = fRandomGen->Rndm()*2*jet_acceptance - jet_acceptance;

	  //loop over all tracks in event, find tracks in random cones and check overlap to signal jets
                                                                                                                 
              for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
	      AliAODTrack *trk = (AliAODTrack*) aodEvent->GetTrack(iTrack);
	      //use hybrid tracks!!!!!!                                                                                                                                              
	      if (!(trk->IsHybridGlobalConstrainedGlobal() && (trk->GetStatus() & AliAODTrack::kITSrefit) )) continue;
                                                  
                                                                                                                                 
	    //Area of random cones                                                                                                                                                         
	    RCdeta = trk->Eta() - RCEta;
	    RCdphi = trk->Phi() - RCPhi;

	    RCdeta = TVector2::Phi_mpi_pi(RCdeta);
	    RCdphi = TVector2::Phi_mpi_pi(RCdphi);

	    if (RCdeta*RCdeta + RCdphi*RCdphi < R*R){
	    
	      //fill RC array, pt, ntracks..
            
	      fRC->Add(trk);
	      RCpt+=trk->Pt();
	      RCntracks++;

	      //std::cout << "This track is in the random cone" << endl;
	      //number of trys:    
              //RCntry++; why here?? not correct..


               //loop over all tracks in signal jets and check if same track as in random cone..
	      for (Int_t itrackjet=0; itrackjet < fJetsTracks->GetEntriesFast(); itrackjet++){
		AliAODTrack *trackinsignaljet = (AliAODTrack*) fJetsTracks->At(itrackjet);
		//		std::cout << "the label of the track in signal jet is: " << trackinsignaljet->GetLabel() << endl;
		//	std::cout << "label of track in random cone: " << trk->GetLabel() << " label of track in signal jet: " << trackinsignaljet->GetLabel() << endl;

		//if (TMath::Abs(trackinsignaljet->GetLabel()) == TMath::Abs(trk->GetLabel())) RCntracksoverlap++;

		//GETLABEL WORKS ONLY FOR MC!!
		if (TMath::Abs(trackinsignaljet->Pt() - trk->Pt()) < 0.0001) RCntracksoverlap++;	  
     
	      }//close loop over tracks in signal jets

	    }//close if track in RC


	  }//close track loop

	      // std::cout << "number of tracks in this RC: " << RCntracks << " Overlap tracks: " << RCntracksoverlap << endl;

	}//close do
	while (RCntracksoverlap>0);

	//plot all the RC properties..
	//but only if the loop was not broken because of too much trys
	
	fHistRCntry->Fill(RCntry);
	if (!breakRC){
        fHistPhiRC->Fill(RCPhi);
	fHistEtaRC->Fill(RCEta);
	fHistRCnTracksOverlap->Fill(RCntracksoverlap);
	fHistRCPt->Fill(RCpt);
	fHistRCnTracks->Fill(RCntracks);
	fHistRCnTracksvsLJPt->Fill(RCntracks, ptleading);
	//random cone pt vs leading jet pt for UE
	fHistRCPtvsLJPt->Fill(RCpt, ptleading);
	
	//define rho for RC:
	rhoRC = RCpt/(TMath::Pi()*R*R);

	fHistRhoRC->Fill(rhoRC);
        fHistRhoRCvsLJPt->Fill(rhoRC, ptleading);


        
        }//close if RC loop was not breaked
	}//close if jet above threshold
	
	//====================================================================          
        //second jet loop for leading jets:
	//===================================================================
	//initialize rho
	Double_t rho = -1;                                                                                                                                                                                                                                                                              

	for (Int_t iJet = 0; iJet < fJets->GetEntriesFast(); ++iJet) {
	  //use jets from array of signal jets.. (dont need to introduce cuts again...)
        AliAODJet *jet = (AliAODJet*) fJets->At(iJet);
	        
	//	if (iJet==0) std::cout << "new event " << endl;
	//	std::cout << "jet pt: " << jet->Pt() << endl;

        if(TMath::Abs(jet->Pt() - ptleading) < 0.001){
	  //	  AliAODJet *leadingjet = (AliAODJet*) jetArray->At(iJet);
	  //at this point 2 % of leading jets are lost... (Events where no Jet exists in Jet acceptance)
	  //  fAntiktHistPtLeading->Fill(ptleading);
	  if (jet->Pt() > jet_ptthreshold){
          fAntiktHistPhiLeading->Fill(jet->Phi());                                                                                                                                          fAntiktHistEtaLeading->Fill(jet->Eta());                                                                                                                                          fAntiktHistEtaPhiLeading->Fill(jet->Eta(), jet->Phi()); 
	  }

	  TRefArray* trackarrayjet = jet->GetRefTracks();
          const Int_t ntracksjet = trackarrayjet->GetEntriesFast();

	  //if (ntracksjet==7){
	  fAntiktHistAreaLeading->Fill(jet->EffectiveAreaCharged());
          fAntiktHistAreaoverCircleLeading->Fill(jet->EffectiveAreaCharged()/(TMath::Pi()*R*R));
	  fAntiktHistAreaLeadingvsPt->Fill(jet->EffectiveAreaCharged(), jet->Pt());
          fAntiktHistAreaoverCircleLeadingvsPt->Fill(jet->EffectiveAreaCharged()/(TMath::Pi()*R*R), jet->Pt());
	  if(jet->Pt()>jet_ptthreshold){
	    fAntiktHistAreaoverCircleLeadingvsNtracks->Fill(jet->EffectiveAreaCharged()/(TMath::Pi()*R*R), ntracksjet);
	  }
           //}

	  //now: get track references from jet and get track properties to plot ff etc..                                                                                                                                                       

	  //	  TRefArray* trackarrayjet = jet->GetRefTracks();
	  // const Int_t ntracksjet = trackarrayjet->GetEntriesFast();
	 
	  //ATTENTION: THIS IS JUST FOR INFORMATION.... TAKE OUT BEFORE NEXT TRAIN
	  // if (jet->EffectiveAreaCharged()/(TMath::Pi()*R*R) > 1.03 && jet->EffectiveAreaCharged()/(TMath::Pi()*R*R) < 1.05 ){
          fAntiktHistnTracksLeading->Fill(ntracksjet);
	  fAntiktHistnTracksLeadingvsPt->Fill(ntracksjet, jet->Pt());
	  

	  Double_t leadingtrackpt=0;
	  Bool_t acceptjetleadingtrack=kFALSE;
	  //loop over tracks associated with the jet                                                                                                                                                                                          
	  for (Int_t itrack=0; itrack < ntracksjet; itrack++){
	    AliAODTrack *trackjet = (AliAODTrack*) jet->GetTrack(itrack);



	    Double_t z = trackjet->Pt()/jet->Pt();
	    Double_t xi = TMath::Log(1/z);

	    fAntiktHistzLeading->Fill(z);
	    fAntiktHistzptLeading->Fill(z, jet->Pt() );
	    fAntiktHistxiLeading->Fill(xi);
	    fAntiktHistxiptLeading->Fill(xi, jet->Pt() );
	    fAntiktHistxitimesptptLeading->Fill(trackjet->Pt()*xi, jet->Pt());

	    if (jet->Pt()>jet_ptthreshold){
            fAntiktHistTrackPtLeading->Fill(trackjet->Pt());
            fAntiktHistTrackPhiLeading->Fill(trackjet->Phi());
            fAntiktHistTrackEtaLeading->Fill(trackjet->Eta());
	    }
	    //determine pt of leading track:
	    if(trackjet->Pt()>leadingtrackpt) leadingtrackpt = trackjet->Pt();

	  }//close loop over tracks in jet       

	  //flag jets with leading track > 5GeV as good jets:
	  if(leadingtrackpt>5) acceptjetleadingtrack=kTRUE;


	  //  }//CLOSE IF IN JET AREA...

	  //==========================================================================
	  //Estimate Background with Perpendicular Cones wrt Leading Jet
	  //==========================================================================

      
	  //Could use the more elegant way...



	  //Need two perpendicular cones for each jet: cone 'plus' and 'minus'
      Double_t ljpcppt=0; // leading jet perp cone plus pt
      Double_t ljpcmpt=0; //leading jet perpendicular cone minus pt                                                                                                                                                                         
      Int_t ljpcpnt=0; //leading jet perpendicular cone plus number of tracks
      Int_t ljpcmnt=0; // lj perp cone minus number of tracks 
     
      //everything that follows is obsolete?
      /*
      Double_t PhiPCP=-1;
      Double_t PhiPCM=-1;
      Double_t PhiPCPA=-1; //Additional cone for cases 4 and 5 
      Double_t PhiPCMA=-1; //Additional cone for cases 4 and 5
      bool case4 = kFALSE; //initialize case 4 and 5 to be false
      bool case5 = kFALSE;
      
      //=========================================================================================================================================
      //Case 1: Jet in the medium phi range         
      if (jet->Phi() > (TMath::Pi()/2 + R) && jet->Phi() < (3*TMath::Pi()/2 - R) ){
	//std::cout << "case 1 is true for this jet"<< endl;

      // define phi of perp cones                                                                                                                                                                                                    
      PhiPCP = jet->Phi() + TMath::Pi()/2;
      PhiPCM = jet->Phi() - TMath::Pi()/2;

     }//close if case 1
      //=========================================================================================================================================
                     
      
      
      //Case 2: Jet in lowest phi range -> pc minus cone has to be shifted in other direction                                                                                                                                                

       else if (jet->Phi() > 0 && jet->Phi() < (TMath::Pi()/2 - R) ){

        // define phi of perp cones                                                                                                                                                                                                                                   
	//std::cout << "case 2 is true for this jet"<< endl;
                                                                                                                                                                                                                 
         PhiPCP = jet->Phi() + TMath::Pi()/2;
         PhiPCM = jet->Phi() + 3*TMath::Pi()/2;

      }//close if case 2
      //==========================================================================================================================================
      

      
      
      //case 3: Jet in highest phi range -> pc plus has to be shifted in other direction 
      
      else if (jet->Phi() > (3*TMath::Pi()/2 + R)  && jet->Phi() < 2*TMath::Pi()  ){

        // define phi of perp cones    
       // std::cout << "case 3 is true for this jet"<< endl;
                                                                                                                                                                                                                
	PhiPCP = jet->Phi() - 3*TMath::Pi()/2;
	PhiPCM = jet->Phi() - TMath::Pi()/2;

      }//close if case 3               
      //==========================================================================================================================================
      
      //case 4: Jet in lower intermediate phi range so that pc minus cone is on the lower edge

      else if (jet->Phi() > (TMath::Pi()/2 - R)  && jet->Phi() < (TMath::Pi()/2 + R)  ){

	case4=kTRUE;
        // define phi of perp cones                             
                                                                                                                                                                          
	//std::cout << "case 4 is true for this jet"<< endl;
                                                                                                                                                                                                                                           
        PhiPCP = jet->Phi() + TMath::Pi()/2;
	PhiPCM = jet->Phi() - TMath::Pi()/2;
        PhiPCMA = jet->Phi() + 3*TMath::Pi()/2;

      }//close if case 4  
      //==========================================================================================================================================

      
      //case 5: Jet in upper intermediate phi range so that pc plus cone is on the upper edge                                                                                                                                               

      else if (jet->Phi() > (3*TMath::Pi()/2 - R)  && jet->Phi() < (3*TMath::Pi()/2 + R)  ){

	case5 = kTRUE;
        // define phi of perp cones                              
	//    std::cout << "case 5 is true for this jet"<< endl;

        PhiPCP = jet->Phi() + TMath::Pi()/2;
        PhiPCPA = jet->Phi() - 3*TMath::Pi()/2;
        PhiPCM = jet->Phi() - TMath::Pi()/2;                                                                                                                                                                                             

      }//close if case 5
      //==========================================================================================================================================
      //
      
      //only needed if not all cases in 0, 2pi are implemented:
     else {// std::cout << "ERROR - no case is true for this jet and the phi of the jet is "<< jet->Phi() << endl; 
        continue;
      }
      //Plot properties of perpendicular cones:
      fLJPCPPhi->Fill(PhiPCP);
      if (case5){ fLJPCPPhi->Fill(PhiPCPA);}
      fLJPCMPhi->Fill(PhiPCM);
      if (case4) {fLJPCMPhi->Fill(PhiPCMA);}
      fLJPCPEta->Fill(jet->Eta());
      fLJPCMEta->Fill(jet->Eta());

      */

      Double_t deta=-9;
      Double_t dphim=-9;
      Double_t dphip=-9;
      bool accepttrackpcm=kFALSE;
      bool accepttrackpcp=kFALSE;
     
      //  bool accepttrackmmm=kFALSE;
      //bool accepttrackmmp=kFALSE;
      ///bool accepttrackmmma=kFALSE;
      //bool accepttrackmmpa=kFALSE;
      //Int_t pcp_count = 0;

      Double_t ljptCM = 0;
      Double_t dphiLJCM=-9;
      Double_t pcppt=0;	


            //now: loop over all tracks in event and find tracks in perpendicular background cones wrt leadingjet
	   for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
	     AliAODTrack *trk = (AliAODTrack*) aodEvent->GetTrack(iTrack);
	    //use hybrid tracks!!!!!!
	    if (!(trk->IsHybridGlobalConstrainedGlobal() && (trk->GetStatus() & AliAODTrack::kITSrefit) )) continue;
	    //elegant method:
	    deta = trk->Eta()-jet->Eta();
            dphim = trk->Phi() - jet->Phi() + TMath::Pi()/2;
	    dphip = trk->Phi() - jet->Phi() - TMath::Pi()/2;
	    dphiLJCM = trk->Phi() - jet->Phi();

	    dphim = TVector2::Phi_mpi_pi(dphim);
            dphip = TVector2::Phi_mpi_pi(dphip);
	    dphiLJCM = TVector2::Phi_mpi_pi(dphiLJCM);

	    if (deta*deta + dphiLJCM*dphiLJCM < R_LJCM*R_LJCM){
	      fTracks_LJCM->Add(trk);
	      ljptCM+=trk->Pt();
	    }

	    if (deta*deta + dphim*dphim < R*R) {
	     accepttrackpcm = kTRUE;
	     //should add track to perp cone minus array..
	     fTracks_pcm->Add(trk);
                 }
	    else accepttrackpcm=kFALSE;

	    // std::cout << "The value of accepttrackjm is " << accepttrackjm << endl;

	    if (deta*deta + dphip*dphip < R_LJCM*R_LJCM) {
	      fTracks_pcp_LJCM->Add(trk);
             }

	    if (deta*deta + dphip*dphip < R*R) {
	     accepttrackpcp = kTRUE;
	     //     ++pcp_count;
	     //Add track in pcp to pcp array:
	     fTracks_pcp->Add(trk); 
	     pcppt+=trk->Pt();
                 } 
	    else accepttrackpcp=kFALSE;

	    //Better: loop over tracks in conearray and plot interesting properties...


	    //My very elegant method...:

              //Find tracks in perp cone plus
	    if (accepttrackpcp) {
		//fill all the histograms...
		fLJPCPTrackPt->Fill(trk->Pt());
                fLJPCPTrackEta->Fill(trk->Eta());
                fLJPCPTrackPhi->Fill(trk->Phi());
		fLJPCPTrackEtaPhi->Fill(trk->Eta(), trk->Phi());
                fLJPCPzpt->Fill(trk->Pt()/jet->Pt(),jet->Pt());
                fLJPCPxipt->Fill(TMath::Log(jet->Pt()/trk->Pt()), jet->Pt());

		ljpcppt+=trk->Pt();
                ++ljpcpnt;

		//	accepttrackmmp=kTRUE;
              }//close perp cone plus
	    // else accepttrackmmp=kFALSE;
	    /*
	      //Find tracks in additional perp cone plus for case 5:
	      if (case5){
               if ( ( (jet->Eta()-trk->Eta())*(jet->Eta()-trk->Eta()) + (PhiPCPA - trk->Phi())*(PhiPCPA - trk->Phi()) ) < R*R ) {
                //fill all the histograms...                                                                                                                                                                                                 
                fLJPCPTrackPt->Fill(trk->Pt());
                fLJPCPTrackEta->Fill(trk->Eta());
                fLJPCPTrackPhi->Fill(trk->Phi());
		fLJPCPTrackEtaPhi->Fill(trk->Eta(), trk->Phi());
		fLJPCPzpt->Fill(trk->Pt()/jet->Pt(),jet->Pt()); 
                fLJPCPxipt->Fill(TMath::Log(jet->Pt()/trk->Pt()), jet->Pt());


                ljpcppt+=trk->Pt();
                ++ljpcpnt;
		accepttrackmmpa=kTRUE;
              }//close additional perp cone plus      
	       else accepttrackmmpa=kFALSE;
	      }//close if case 5
	    */
             //Find tracks in perp cone minus

	      if ( accepttrackpcm ) {

		//fill all the histograms...                                                                                                                                                                                                 
                fLJPCMTrackPt->Fill(trk->Pt());
                fLJPCMTrackEta->Fill(trk->Eta());
                fLJPCMTrackPhi->Fill(trk->Phi());
		fLJPCMTrackEtaPhi->Fill(trk->Eta(), trk->Phi());
		fLJPCMzpt->Fill(trk->Pt()/jet->Pt(),jet->Pt()); 
                fLJPCMxipt->Fill(TMath::Log(jet->Pt()/trk->Pt()), jet->Pt());


                ljpcmpt+=trk->Pt();
                ++ljpcmnt;
		//accepttrackmmm=kTRUE;
	      }//close perp cone minus
	      // else accepttrackmmm=kFALSE;
	      /*
	      //Find tracks in additional perp cone minus for case 4:
	      if (case4){
              if ( ( (jet->Eta()-trk->Eta())*(jet->Eta()-trk->Eta()) + (PhiPCMA - trk->Phi())*(PhiPCMA - trk->Phi()) ) < R*R ) {

                //fill all the histograms...                                                                                                                                                                                                                                             
                fLJPCMTrackPt->Fill(trk->Pt());
                fLJPCMTrackEta->Fill(trk->Eta());
                fLJPCMTrackPhi->Fill(trk->Phi());
		fLJPCMTrackEtaPhi->Fill(trk->Eta(), trk->Phi());
                fLJPCMzpt->Fill(trk->Pt()/jet->Pt(),jet->Pt()); 
                fLJPCMxipt->Fill(TMath::Log(jet->Pt()/trk->Pt()), jet->Pt());


                ljpcmpt+=trk->Pt();
                ++ljpcmnt;
		accepttrackmmma=kTRUE;
              }//close additional perp cone minus 
	      else accepttrackmmma=kFALSE;       
	      }//close if case 4

	      */
	      //now:check if same tracks are accepted by both methods:
	      //check if both tracks are accepted by both methods
	      // if (accepttrackjm != (accepttrackmmm || accepttrackmmma)) {std:: cout << "ERROR - Track is not accepted by both methods "<< endl;}
	      //else {std::cout << "SUCCES - Track is accepted or declined by both methods " << endl;}

	      // if (accepttrackjp != (accepttrackmmp || accepttrackmmpa)) {std:: cout << "ERROR - Track is not accepted by both methods "<< endl;}
	      // else {std::cout << "SUCCES - Track is accepted or declined by both methods " << endl;}

	  
	   }//close track loop

              //Fill pt and number of tracks in perp cones
	      
              fLJPCPPt->Fill(ljpcppt);
	      fLJPCPnTracks->Fill(ljpcpnt);
	      // fLJPCPnTracksPhi->Fill(ljpcpnt, PhiPCP);
	      //      if (case5){ fLJPCPnTracksPhi->Fill(ljpcpnt,PhiPCPA);} makes no sense because the number of tracks is already the sum of normal and additional cone
	      
	      // fLJPCMnTracksPhi->Fill(ljpcmnt, PhiPCM);
	      // if (case4) {fLJPCMnTracksPhi->Fill(ljpcmnt, PhiPCMA);}
              fLJPCMPt->Fill(ljpcmpt);
              fLJPCMnTracks->Fill(ljpcmnt);

	      //Fill histos perp cone pt vs leading jet pt
	      fLJPCMPtvsLJPt->Fill(ljpcmpt, jet->Pt());
	      fLJPCPPtvsLJPt->Fill(ljpcppt, jet->Pt());
	      

	      //Use both cones. Up to now no veto, but one could think about...
	      fLJPCPt->Fill(ljpcmpt);
              fLJPCPt->Fill(ljpcppt);
	      fLJPCPtvsLJPt->Fill(ljpcmpt, jet->Pt());
	      fLJPCPtvsLJPt->Fill(ljpcppt, jet->Pt());

	      //=====================================================================================================================================
	      //============================================================
              //Event by event UE subtraction
	      //============================================================

	      //Better: For inclusive jets: Define RHO and then subtract RHO from inclusive jets in loop 
	      rho = (ljpcppt+ljpcmpt)/(2*TMath::Pi()*R*R);  

	      fHistRho->Fill(rho);                                                                      
	      fHistRhovsLJPt->Fill(rho, ptleading); 	  


	      //First attempt: subtract mean pt from both cones from leading jet pt   
                //try: use mean value of PCP and PCM
	      Double_t ptleading_sub_pcp= jet->Pt() - (ljpcppt+ljpcmpt)*jet->EffectiveAreaCharged()/(2*TMath::Pi()*R*R); 

	      //Average UE pt subtraction method
	      //     Double_t ptleading_sub_pcp = 0;
	      fLJPCPMMPt->Fill((ljpcppt+ljpcmpt)/2);

	      //p-Pb:
	      /*
            a  if (jet->Pt()>5 && jet->Pt()<10) ptleading_sub_pcp = jet->Pt() - 1.4;
	      else if (jet->Pt()>10 && jet->Pt()<15) ptleading_sub_pcp = jet->Pt() - 1.66;
	      else if (jet->Pt()>15 && jet->Pt()<20) ptleading_sub_pcp = jet->Pt() - 1.62;
	      else if (jet->Pt()>20) ptleading_sub_pcp = jet->Pt() - 1.55;
	      */
		      //pp
	      // if (jet->Pt()>5 && jet->Pt()<10) ptleading_sub_pcp = jet->Pt() - 0.72;
	      // else if(jet->Pt()>10) ptleading_sub_pcp = jet->Pt() - 0.75;

	      //CROSSCHECK:
              //Double_t ptleading_sub_pcp= jet->Pt();
              fAntiktHistPtLeading_sub_pcp->Fill(ptleading_sub_pcp); 
             //now: loop again over all tracks in leading jet to plot FF with subtracted jet pt

	      for (Int_t itrack=0; itrack < ntracksjet; itrack++){
		AliAODTrack *trackjet = (AliAODTrack*) jet->GetTrack(itrack);



		Double_t z = trackjet->Pt()/ptleading_sub_pcp;
		Double_t xi = TMath::Log(1/z);

		//fAntiktHistzLeading->Fill(z);
		//fAntiktHistzptLeading->Fill(z, ptleading_sub_pcp );
		
                //fAntiktHistxiLeading_sub_pcp->Fill(xi);
		fAntiktHistxiptLeading_sub_pcp->Fill(xi, ptleading_sub_pcp );
		//fAntiktHistxitimesptptLeading_sub_pcp->Fill(trackjet->Pt()*xi, ptleading_sub_pcp);

		if (ptleading_sub_pcp>5){
		  fAntiktHistTrackPtLeading_sub_pcp->Fill(trackjet->Pt());
		 }
	      }//close loop over tracks in jet   

	     //now: loop over all tracks in perp cone 'plus' to plot FF with subtracted jet pt
	      //ATTENTION: JUST USE LJCM NOW AS A CROSSCHECK

	     Int_t ntracks_pcp_LJCM = fTracks_pcp_LJCM->GetEntriesFast();
	     //	     if (ntracks_pcp != ljpcpnt){	    
	     // std::cout << "The number of entries in pcp array and number of tracks in pcp don't match" << endl;
	     // }

	     for(Int_t itrack_pcp_LJCM=0; itrack_pcp_LJCM < ntracks_pcp_LJCM; itrack_pcp_LJCM++){

	       AliAODTrack *track_pcp = (AliAODTrack*) fTracks_pcp_LJCM->At(itrack_pcp_LJCM);
	       
               Double_t z = track_pcp->Pt()/ptleading_sub_pcp;
	       Double_t xi = TMath::Log(1/z);

                                                                                                                                                    
	        fAntiktHistxiptLeading_pcp_sub_pcp->Fill(xi, ptleading_sub_pcp );
	       //fAntiktHistxitimesptptLeading_sub_pcp->Fill(trackjet->Pt()*xi, ptleading_sub_pcp);                                                                                                                                         

	        fAntiktHistTrackPtLeading_pcp_sub_pcp->Fill(track_pcp->Pt());


		//perp cone ff's for leading jet cone method:	       
		Double_t z_LJCM = track_pcp->Pt()/ljptCM;
		Double_t xi_LJCM = TMath::Log(1/z_LJCM);
		fAntiktHistxiptLeading_pcp_LJCM->Fill(xi_LJCM, ljptCM );



	     }//close loop over tracks in perp cone plus

	     //=========================================
	     //CONE AROUNG LEADING JET AXIS METHOD
	     //=============================================
	     Int_t ntracks_LJCM=  fTracks_LJCM->GetEntriesFast();
	     fAntiktHistPtLeading_LJCM->Fill(ljptCM);

	     for (Int_t itrack_LJCM=0; itrack_LJCM< ntracks_LJCM; itrack_LJCM++){

       AliAODTrack *track_LJCM = (AliAODTrack*) fTracks_LJCM->At(itrack_LJCM);

       Double_t z = track_LJCM->Pt()/ljptCM;
       Double_t xi = TMath::Log(1/z);


      fAntiktHistxiptLeading_LJCM->Fill(xi, ljptCM );
 //fAntiktHistxitimesptptLeading_sub_pcp->Fill(trackjet->Pt()*xi, ptleading_sub_pcp);                                                                                                                                         

      fAntiktHistTrackPtLeading_LJCM->Fill(track_LJCM->Pt());

	     }//close loop over tracks in leading jet with cone method

     //===================================
     //LEADING TRACK CUT (>5GeV)
     //=====================================

	     if(acceptjetleadingtrack){

	       fAntiktHistPtLeading_LTC->Fill(jet->Pt());
	    
	     for (Int_t itrack=0; itrack < ntracksjet; itrack++){
	       AliAODTrack *trackjet = (AliAODTrack*) jet->GetTrack(itrack);

	    
	       Double_t z = trackjet->Pt()/jet->Pt();
	       Double_t xi = TMath::Log(1/z);

	       //fAntiktHistzLeading->Fill(z);
	       //fAntiktHistzptLeading->Fill(z, jet->Pt() );
	       //fAntiktHistxiLeading->Fill(xi);
	       fAntiktHistxiptLeading_LTC->Fill(xi, jet->Pt() );
	     }//close loop over tracks in leading jet

	     Int_t ntracks_pcp = fTracks_pcp->GetEntriesFast();

	     for(Int_t itrack_pcp=0; itrack_pcp < ntracks_pcp; itrack_pcp++){

               AliAODTrack *track_pcp = (AliAODTrack*) fTracks_pcp->At(itrack_pcp);

               Double_t z = track_pcp->Pt()/jet->Pt();
               Double_t xi = TMath::Log(1/z);


	       fAntiktHistxiptLeading_LTC_pcp->Fill(xi, jet->Pt() );
                                                                                                                                                     
                                                                                                                                                  
	        
	          }//close loop over tracks in perp cone plus   
	      	     
	     }//close if accepted jet (leading track cut)
	   
           }//close if leadingjet       
	}//close second jet loop


	//need new jet loop!
        for (Int_t iJet = 0; iJet < fJets->GetEntriesFast(); ++iJet) {
          //use jets from array of signal jets.. (dont need to introduce cuts again...)                                                                                            
	  AliAODJet *jet = (AliAODJet*) fJets->At(iJet);


	//here: in loop over all signal jets..
	//now loop over tracks in both perp cones to plot inclusive jet ff contribution
	  /*
	//plus cone:
	for(Int_t itrack_pcp=0; itrack_pcp < fTracks_pcp->GetEntriesFast(); itrack_pcp++){

	  AliAODTrack *track_pcp = (AliAODTrack*) fTracks_pcp->At(itrack_pcp);

	  Double_t z = track_pcp->Pt()/jet->Pt();
	  Double_t xi = TMath::Log(1/z);


	  fAntiktHistxipt_pcp->Fill(xi, jet->Pt() );
	  fAntiktHistxipt_pc->Fill(xi, jet->Pt() );

	}//close loop over tracks in perp cone plus

	//minus cone:                                                                                                                              
        for(Int_t itrack_pcm=0; itrack_pcm < fTracks_pcm->GetEntriesFast(); itrack_pcm++){

          AliAODTrack *track_pcm = (AliAODTrack*) fTracks_pcm->At(itrack_pcm);

          Double_t z = track_pcm->Pt()/jet->Pt();
          Double_t xi = TMath::Log(1/z);


          fAntiktHistxipt_pcm->Fill(xi, jet->Pt() );
	  fAntiktHistxipt_pc->Fill(xi, jet->Pt() );
        }//close loop over tracks in perp cone minus
	  */
	//======================================
	//EVENT BY EVENT SUBTRACTION
	//======================================
	//(we are still in the loop over inclusive jets)	

        //now: subtract RHO from jet pt and plot observables of interest with subtracted jet pt

	Double_t ptsub = jet->Pt() - rho*jet->EffectiveAreaCharged();
	Double_t ptsubRC = jet->Pt() - rhoRC*jet->EffectiveAreaCharged();

	//delta pT
	Double_t deltapt = (rhoRC - rho)*TMath::Pi()*R*R;
 
	fHistDeltaPt->Fill(deltapt);
	fHistDeltaPtvsLJPt->Fill(deltapt, ptleading);


        //Now: plot all jet observables with ptsub

	fHistRhoSubPt->Fill(ptsub);
	fHistRhoRCSubPt->Fill(ptsubRC);

	//eta, phi, area
	fHistRhoSubAreaOverCirclevsPt->Fill(jet->EffectiveAreaCharged()/(TMath::Pi()*R*R), ptsub);
        fHistRhoRCSubAreaOverCirclevsPt->Fill(jet->EffectiveAreaCharged()/(TMath::Pi()*R*R), ptsubRC);



	TRefArray* trackarrayjet = jet->GetRefTracks();
	const Int_t ntracksjet = trackarrayjet->GetEntriesFast();
	//xi jet 
	for (Int_t itrack=0; itrack < ntracksjet; itrack++){
	  AliAODTrack *trackjet = (AliAODTrack*) jet->GetTrack(itrack);
	  //PC RHO subtracted 
	  //eta phi of tracks in jets
          fHistRhoSubPhiTracksJetvsPt->Fill(trackjet->Phi(), ptsub);
	  fHistRhoSubEtaTracksJetvsPt->Fill(trackjet->Eta(), ptsub);


	  Double_t z = trackjet->Pt()/ptsub;
	  Double_t xi = TMath::Log(1/z);
                                                                                                
	  fHistRhoSubXiPt->Fill(xi, ptsub );

	  //RC RHO subtracted


	  fHistRhoRCSubPhiTracksJetvsPt->Fill(trackjet->Phi(), ptsubRC);
          fHistRhoRCSubEtaTracksJetvsPt->Fill(trackjet->Eta(), ptsubRC);


          Double_t zRC = trackjet->Pt()/ptsubRC;
          Double_t xiRC = TMath::Log(1/zRC);

          fHistRhoRCSubXiPt->Fill(xiRC, ptsubRC );

	}//close loop over tracks in jet 
 
	//xi perp cone

	//plus cone:                                                                                                                           
        for(Int_t itrack_pcp=0; itrack_pcp < fTracks_pcp->GetEntriesFast(); itrack_pcp++){

          AliAODTrack *track_pcp = (AliAODTrack*) fTracks_pcp->At(itrack_pcp);
	  //unsubtracted:
          Double_t z = track_pcp->Pt()/jet->Pt();
          Double_t xi = TMath::Log(1/z);
	  fAntiktHistxipt_pcp->Fill(xi, jet->Pt() );
          fAntiktHistxipt_pc->Fill(xi, jet->Pt() );


	  //subtracted:
	  Double_t zsub = track_pcp->Pt()/ptsub;
          Double_t xisub = TMath::Log(1/zsub);

          fHistRhoSubXiPCPt->Fill(xisub, ptsub);
	 
         //eta phi of tracks in PC
	  fHistRhoSubPhiTracksPCvsPt->Fill(track_pcp->Phi(), ptsub);
          fHistRhoSubEtaTracksPCvsPt->Fill(track_pcp->Eta(), ptsub);



        }//close loop over tracks in perp cone plus                                                               

        //minus cone:                                                                                                  
                                                                                                                                               
        for(Int_t itrack_pcm=0; itrack_pcm < fTracks_pcm->GetEntriesFast(); itrack_pcm++){

          AliAODTrack *track_pcm = (AliAODTrack*) fTracks_pcm->At(itrack_pcm);
	  //unsubtracted:
          Double_t z = track_pcm->Pt()/jet->Pt();
          Double_t xi = TMath::Log(1/z);

          fAntiktHistxipt_pcm->Fill(xi, jet->Pt() );
          fAntiktHistxipt_pc->Fill(xi, jet->Pt() );
	  //subtracted:
	  Double_t zsub = track_pcm->Pt()/ptsub;
          Double_t xisub = TMath::Log(1/zsub);

          fHistRhoSubXiPCPt->Fill(xisub, ptsub);
	  
          //eta phi of tracks in PC
	  fHistRhoSubPhiTracksPCvsPt->Fill(track_pcm->Phi(), ptsub);
          fHistRhoSubEtaTracksPCvsPt->Fill(track_pcm->Eta(), ptsub);


        }//close loop over tracks in perp cone minus    


	//Random Cone:

	                                                                                                                                                        

        for(Int_t itrack_RC=0; itrack_RC < fRC->GetEntriesFast(); itrack_RC++){

          AliAODTrack *track_RC = (AliAODTrack*) fRC->At(itrack_RC);
          //unsubtracted:                                                                                                                                                          
          Double_t z = track_RC->Pt()/jet->Pt();
          Double_t xi = TMath::Log(1/z);

          fAntiktHistxipt_RC->Fill(xi, jet->Pt() );
          //subtracted:                                                                                                                                                            
          Double_t zsubRC = track_RC->Pt()/ptsubRC;
          Double_t xisubRC = TMath::Log(1/zsubRC);

          fHistRhoRCSubXiPCPt->Fill(xisubRC, ptsubRC);

          //eta phi of tracks in PC                                                                                                                                                
          fHistRhoRCSubPhiTracksPCvsPt->Fill(track_RC->Phi(), ptsubRC);
          fHistRhoRCSubEtaTracksPCvsPt->Fill(track_RC->Eta(), ptsubRC);


        }//close loop over tracks in random cone                                                                                                              






	}//close third jet loop

	//clear event arrays:
	fRC->Clear();
	fJets->Clear();
	fJetsTracks->Clear();
	fTracks_pcp->Clear();
	fTracks_pcm->Clear();
	fTracks_pcp_LJCM->Clear();
	fTracks_LJCM->Clear();

	//MC Truth arrays
	fMCTruthParticles->Clear();
	fMCRC->Clear();
	fMCtracks_pcp->Clear();
	fMCtracks_pcm->Clear();
	fMCjets->Clear();
        fMCjetstracks->Clear();

     } // close if jetarray
  }// close if aodevent


  //  std:: cout << track_IDs[1] << endl;

  PostData(1, fOutputList);
      }//. close userexec

void AliAnalysisTaskFragmentationTriggered::Terminate(const Option_t * /* option */)
{
  return;
}
