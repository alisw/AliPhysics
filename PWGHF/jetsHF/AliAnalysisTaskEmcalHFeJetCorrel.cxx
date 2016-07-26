//////////////////////////////////////////////////////
//													//
// Heavy Flavour Correlations and Jet Analysis		//
//													//
// Authors: Diogenes D. Gimenez, Andrea Rossi	//
//////////////////////////////////////////////////////


#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include "TRandom.h"
#include "TRandom3.h"

#include "THnSparse.h"

#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskEmcalHFeJetCorrel.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"


#include "AliAODEvent.h"

#include "AliAODRecoDecay.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"

#include "AliAODMCHeader.h"

#include "AliAnalysisTaskEmcalHFeJetCorrel.h"
#include "AliRDHFJetsCuts.h"
#include "AliAnalysisVertexingHF.h"

#include "AliPIDResponse.h"
#include "AliVParticle.h"
#include "AliVCluster.h"
#include <TArrayI.h>



ClassImp(AliAnalysisTaskEmcalHFeJetCorrel)

using std::cout;
using std::endl;

//________________________________________________________________________
AliAnalysisTaskEmcalHFeJetCorrel::AliAnalysisTaskEmcalHFeJetCorrel() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalHFeJetCorrel", kTRUE),
//======================================================================
//Containers
//fOutput(0),
  fJetsCont(0),
  fMCJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0),
//======================================================================

//======================================================================
//Flags
  fReadMC(0),
  kQA(kFALSE),
  kCheckClusterMatching(kFALSE),
  fCheckVzero(kFALSE),
  kAnalysis(kFALSE),
  kGeneralSpectra(kFALSE),
  kDetector(kTRUE),
//======================================================================

//======================================================================
//Variables
  fdebug(-1),
  fNentries(0x0),
  fNRejected(0x0),
  fcandEleTPC(0x0),
  fLastEleTPC(0),
  fCutsHFjets(0x0),
  fCutsElectron(0x0),
  fpidResp(0x0),
  ffilterbit(AliAODTrack::kTrkTPCOnly),
  fPhiTrackClusterDistance(0.05),
  fEtaTrackClusterDistance(0.05),
  fsigmaTPCmin(-10.),
  fsigmaTPCmax(3.5),
  fminpt(.1),
  fEtaMin(-0.9),
  fEtaMax(0.9),
  fMinEoverP(0.),
  fMaxEoverP(3.),
  fMassPhotonicCut(0.5),
  fminNcell(2),
  fmaxM20(0.6),
  fmaxM02(0.6),
  fMCParticlesName(0x0),
  fMCJetsBranch(0x0),
//======================================================================

//======================================================================
//Histograms
  fhTrackRejection(0x0),
  fhEleRejection(0x0),
  fSparseHFSpectrum(0x0),
  fSparseHFTriggers(0x0),
  fSparseJet(0x0),
  fSparseMCJet(0x0),
  fSparseQATracks(0x0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistPtDEtaDPhiTrackClus(0),
  fhPhotonicEle(0x0),
  fDetector(0x0),
  fRandomCones(0x0)
//======================================================================

{
  // Default constructor.
  SetMakeGeneralHistograms(kFALSE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalHFeJetCorrel::AliAnalysisTaskEmcalHFeJetCorrel(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),

//======================================================================
//Containers
//fOutput(0),
  fJetsCont(0),
  fMCJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0),
//======================================================================

//======================================================================
//Flags
  fReadMC(0),
  kQA(kFALSE),
  kCheckClusterMatching(kFALSE),
  fCheckVzero(kFALSE),
  kAnalysis(kFALSE),
  kGeneralSpectra(kFALSE),
  kDetector(kTRUE),
//======================================================================

//======================================================================
//Variables
  fdebug(-1),
  fNentries(0x0),
  fNRejected(0x0),
  fcandEleTPC(0x0),
  fLastEleTPC(0),
  fCutsHFjets(0x0),
  fCutsElectron(0x0),
  fpidResp(0x0),
  ffilterbit(AliAODTrack::kTrkTPCOnly),
  fPhiTrackClusterDistance(0.05),
  fEtaTrackClusterDistance(0.05),
  fsigmaTPCmin(-10.),
  fsigmaTPCmax(3.5),
  fminpt(.1),
  fEtaMin(-0.9),
  fEtaMax(0.9),
  fMinEoverP(0.),
  fMaxEoverP(3.),
  fMassPhotonicCut(0.5),
  fminNcell(2),
  fmaxM20(0.6),
  fmaxM02(0.6),
  fMCParticlesName(0x0),
  fMCJetsBranch(0x0),
//======================================================================

//======================================================================
//Histograms
  fhTrackRejection(0x0),
  fhEleRejection(0x0),
  fSparseHFSpectrum(0x0),
  fSparseHFTriggers(0x0),
  fSparseJet(0x0),
  fSparseMCJet(0x0),
  fSparseQATracks(0x0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistPtDEtaDPhiTrackClus(0),
  fhPhotonicEle(0x0),
  fDetector(0x0),
  fRandomCones(0x0)
//======================================================================

{
  // Standard constructor.
  SetMakeGeneralHistograms(kFALSE);
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalHFeJetCorrel::~AliAnalysisTaskEmcalHFeJetCorrel()
{
  // Destructor.
//if(fOutput) delete fOutput;

  if(fJetsCont) delete fJetsCont;
  if(fMCJetsCont) delete fMCJetsCont;
  if(fTracksCont) delete fTracksCont;
  if(fCaloClustersCont) delete fCaloClustersCont;
  if(fNentries) delete fNentries;
  if(fNRejected) delete fNRejected;
  if(fcandEleTPC) delete fcandEleTPC;

  if(fhTrackRejection) delete fhTrackRejection;
  if(fhEleRejection) delete fhEleRejection;
  if(fCutsHFjets) delete fCutsHFjets;
  if(fSparseHFSpectrum) delete fSparseHFSpectrum;
  if(fSparseHFTriggers) delete fSparseHFTriggers;
  if(fSparseJet) delete fSparseJet;
  if(fSparseMCJet) delete fSparseMCJet;
  if(fSparseQATracks) delete fSparseQATracks;
  if(fHistPtDEtaDPhiClusTrack) delete fHistPtDEtaDPhiClusTrack;
  if(fHistPtDEtaDPhiTrackClus) delete fHistPtDEtaDPhiTrackClus;
  if(fhPhotonicEle) delete fhPhotonicEle;
  if(fCutsElectron) delete fCutsElectron;
  if(fDetector) delete fDetector;
  if(fRandomCones) delete fRandomCones;
  if(fpidResp) delete fpidResp;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::UserCreateOutputObjects()
{
  // Create user output.
//Printf("STARTING CREATION OF OUTPUT OBJECTS");
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fJetsCont                 = GetJetContainer(0);
    if(fJetsCont)
  { //get particles and clusters connected to jets
    fTracksCont       = fJetsCont->GetParticleContainer();
    fCaloClustersCont = fJetsCont->GetClusterContainer();
  }
  else
  {//no jets, just analysis tracks and clusters
    fTracksCont       = GetParticleContainer(0);
    fCaloClustersCont = GetClusterContainer(0);
  }
  fTracksCont->SetClassName("AliPicoTrack");
  fCaloClustersCont->SetClassName("AliVCluster");

//  fOutput = new TList();
//  fOutput->SetOwner(kTRUE);
//=====================================================================================================================
//Andrea's Outputs

//=====================================================================================================================
//Event Rejection Histogram
  fNentries=new TH1F("nentriesChFr", "Analyzed sample properties", 9,-0.5,8.5);
  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nEvSel");
  fNentries->GetXaxis()->SetBinLabel(3,"nEvPile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(4,"nEvGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"nEvRejVtxZ");
  fNentries->GetXaxis()->SetBinLabel(6,"nTracksEv");
  fNentries->GetXaxis()->SetBinLabel(7,"nJetsCand");
  fNentries->GetXaxis()->SetBinLabel(8,"nJetsTagged");
  fNentries->GetXaxis()->SetBinLabel(9,"nUnexpError");
  fOutput->Add(fNentries);
//=====================================================================================================================

//=====================================================================================================================
//Event Selection Rejection Histogram
fNRejected=new TH1F("SelectionRejection", "Analyzed Rejection Reasons", 9,-0.5,8.5);
fNRejected->GetXaxis()->SetBinLabel(1,"NoCentrality");
fNRejected->GetXaxis()->SetBinLabel(2,"PileUp");
fNRejected->GetXaxis()->SetBinLabel(3,"Centrality1");
fNRejected->GetXaxis()->SetBinLabel(4,"Centrality2");
fNRejected->GetXaxis()->SetBinLabel(5,"Centrality3");
fNRejected->GetXaxis()->SetBinLabel(6,"TriggerClass");
fNRejected->GetXaxis()->SetBinLabel(7,"Vertex");
fNRejected->GetXaxis()->SetBinLabel(8,"PhysicsSelection");
fOutput->Add(fNRejected);
//=====================================================================================================================

//=====================================================================================================================
//MCElectron Rejection Histogram
if(fReadMC)
{
  fhEleRejection=new TH3F("fhEleRejection","fhEleRejection; step; pt; nsigma;",20,-0.5,19.5,102,-1.,50.,12,-3,3);
  fhEleRejection->GetXaxis()->SetBinLabel(1,"   ");
  fhEleRejection->GetXaxis()->SetBinLabel(2,"All");
  fhEleRejection->GetXaxis()->SetBinLabel(3,"No Clst");
  fhEleRejection->GetXaxis()->SetBinLabel(4,"Filt Bit");
  fhEleRejection->GetXaxis()->SetBinLabel(5,"NoEMCal");
  fhEleRejection->GetXaxis()->SetBinLabel(6,"TrkClstDistAndNcell");
  fhEleRejection->GetXaxis()->SetBinLabel(7,"M20");
  fhEleRejection->GetXaxis()->SetBinLabel(8,"Min pt");
  fhEleRejection->GetXaxis()->SetBinLabel(9,"IsDaugSel");
  fhEleRejection->GetXaxis()->SetBinLabel(10,"nsigmaTPC");
  fhEleRejection->GetXaxis()->SetBinLabel(11,"E/p");
    fhEleRejection->GetXaxis()->SetBinLabel(12,"Is Photonic");
  fhEleRejection->GetXaxis()->SetBinLabel(13,"SelTrack");
  fhEleRejection->GetXaxis()->SetBinLabel(14,"IsVzeroDaugh");
  fOutput->Add(fhEleRejection);
}
//=====================================================================================================================

//=====================================================================================================================
//Track Rejection Histogram
  fhTrackRejection=new TH3F("fhTrackRejection","fhTrackRejection; step; pt; nsigma;",20,-0.5,19.5,102,-1.,50.,12,-3,3);
  fhTrackRejection->GetXaxis()->SetBinLabel(1,"   ");
  fhTrackRejection->GetXaxis()->SetBinLabel(2,"All");
  fhTrackRejection->GetXaxis()->SetBinLabel(3,"No Clst");
  fhTrackRejection->GetXaxis()->SetBinLabel(4,"Filt Bit");
  fhTrackRejection->GetXaxis()->SetBinLabel(5,"NoEMCal");
  fhTrackRejection->GetXaxis()->SetBinLabel(6,"TrkClstDistAndNcell");
  fhTrackRejection->GetXaxis()->SetBinLabel(7,"M20");
  fhTrackRejection->GetXaxis()->SetBinLabel(8,"Min pt");
  fhTrackRejection->GetXaxis()->SetBinLabel(9,"IsDaugSel");
  fhTrackRejection->GetXaxis()->SetBinLabel(10,"nsigmaTPC");
  fhTrackRejection->GetXaxis()->SetBinLabel(11,"E/p");
  
  fhTrackRejection->GetXaxis()->SetBinLabel(12,"Is Photonic");
  fhTrackRejection->GetXaxis()->SetBinLabel(13,"SelTrack");
  fhTrackRejection->GetXaxis()->SetBinLabel(14,"IsVzeroDaugh");
  fOutput->Add(fhTrackRejection);
//=====================================================================================================================


//=====================================================================================================================
if(kAnalysis)
{  //study of ELE-jets:pt ele, eta ele, phi ele, nsigma, pt jet, deltaeta, deltaphi, flag istrackinjet
  // MC histo (ele id from pdg code)
if(fdebug>1) cout<<"DDG Outputs 0"<<endl;

  const int NTracks = 16;
    //                                1,  2,   3,  4,   5,   6,   7,  8,  9, 10,11,  12,       13,14, 15,16
  Int_t nbinsTrack[NTracks]		=  {50 ,150,  90, 90,  54,   8,   6,   6,  9, 16, 5, 300,       20,40, 40, 3};
  Double_t binlowTrack[NTracks]	=  { 0., 0 ,  70, 70,-10.,-0.5, 0.0, 0.0,0.0,0.0,-1,-100,-(TMath::Pi()*0.5),-2,  0,-1};
  Double_t binupTrack[NTracks]	=  {50., 3., 160,160, 3.5, 7.5,0.06,0.06,0.9,1.6, 4,200.,(1.5*TMath::Pi()) , 2,  4, 2};
    //                              pT,E/p,TPCcl,PID,NSig, ITS, Phi, Eta,M20,M02,NLM,pT,       Phi,Eta,R,InJet
  TString titles = "fSparseHFSpectrum; p_{T}; E/p; N_{clusters}^{TPC}; N_{PID}^{TPC}; N_{#sigma TPC}; N_{clusters}^{ITS}; #Delta #varphi ^{Match}; #Delta #eta ^{Match}; M20; M02; NLM; p_{T}^{jet}; #Delta #varphi; #Delta #eta; R_{#varphi #eta};InJet;";
  // data histo (ele id from detector id, filled in MC as well)
  fSparseHFSpectrum = new THnSparseF("fSparseHFSpectrum", titles.Data(),NTracks, nbinsTrack, binlowTrack, binupTrack);
  fOutput->Add(fSparseHFSpectrum);

    const int NTriggers = 11;
    //                                        1,  2,   3,  4,   5,   6,   7,   8,  9, 10,11
    Int_t nbinsTriggers[NTriggers]		=  {50 ,150,  90, 90,  54,   8,   6,   6,  9, 16, 5};
    Double_t binlowTriggers[NTriggers]	=  { 0., 0 ,  70, 70,-10.,-0.5, 0.0, 0.0,0.0,0.0,-1};
    Double_t binupTriggers[NTriggers]	=  {50., 3., 160,160, 3.5, 7.5,0.06,0.06,0.9,1.6, 4};
    //                                      pT,E/p,TPCcl,PID,NSig, ITS, Phi, Eta,M20,M02,NLM
    TString titlesTriggers = "fSparseHFTriggers; p_{T}; E/p; N_{clusters}^{TPC}; N_{PID}^{TPC}; N_{#sigma TPC}; N_{clusters}^{ITS}; #Delta #varphi ^{Match}; #Delta #eta ^{Match}; M20; M02; NLM;";
    // data histo (ele id from detector id, filled in MC as well)
    fSparseHFTriggers = new THnSparseF("fSparseHFTriggers", titlesTriggers.Data(),NTriggers, nbinsTriggers, binlowTriggers, binupTriggers);
    fOutput->Add(fSparseHFTriggers);
}
if(fReadMC)
{
}
if(fdebug>1) cout<<"DDG Outputs 1"<<endl;
//=====================================================================================================================

//=====================================================================================================================
//Diogenes's Histograms
//=====================================================================================================================

if(fdebug>1) cout<<"DDG Outputs 2"<<endl;

  const int NCones = 6;
  Int_t nbinsCones[NCones]		= {20, 50,  50, 300, 40,   50};
  Double_t binlowCones[NCones]	= { 0, -1,   0,-100, -1,    0};
  Double_t binupCones[NCones]	= {20,  1,2*TMath::Pi(),200.,  1,2.*TMath::Pi()};
  fRandomCones = new THnSparseF("fRandomCones", "fRandomCones; Constituents; #eta _{Random}; #varphi _{Random}; p_{T};#eta; #varphi;", NCones, nbinsCones, binlowCones, binupCones);
  fOutput->Add(fRandomCones);

//=====================================================================================================================

//=====================================================================================================================
//Invariant Mass Histogram
if(fdebug>1) cout<<"DDG Outputs 3"<<endl;
  const int NPhotonic = 22;
//                                        1,  2,   3,  4,   5,   6,   7,  8,  9, 10,11,  12,       13,14, 15,16,  17,  18,  19, 20, 21,  22
Int_t nbinsPhotonic[NPhotonic]		=  {50 ,150,  90, 90,  54,   8,   6,   6,  9, 16, 5, 300,       20,40, 40, 3,  50,  90,  54, 80,  6,   3};
Double_t binlowPhotonic[NPhotonic]	=  { 0., 0 ,  70, 70,-10.,-0.5, 0.0, 0.0,0.0,0.0,-1,-100,-(TMath::Pi()*0.5),-2,  0,-1,   0,  70, -10,  0,0.0,-1.5};
Double_t binupPhotonic[NPhotonic]	=  {50., 3., 160,160, 3.5, 7.5,0.06,0.06,0.9,1.6, 4,200.,(1.5*TMath::Pi()) , 2,  4, 2,  50, 160, 3.5,0.8,0.6, 1.5};
//                                       pT,E/p,TPCcl,PID,NSig,ITS, Phi, Eta,M20,M02,NLM,pT,       Phi,Eta,R,InJet,pT,TPCcl,Nsig,m_,theta,Sign";
 TString titlesPhotonic = "fhPhotonicEle; p_{T}; E/p; N_{clusters}^{TPC}; N_{PID}^{TPC}; N_{#sigma TPC}; N_{clusters}^{ITS}; #Delta #varphi ^{Match}; #Delta #eta ^{Match}; M20; M02; NLM; p_{T}^{jet}; #Delta #varphi; #Delta #eta; R_{#varphi #eta};InJet; p_{T}^{partner}; N_{clusters}^{TPC, part}; N_{#sigma TPC}^{e, part};m_{inv};#theta_{open};Sign;";
    
 fhPhotonicEle		= new THnSparseF("fhPhotonicEle",titlesPhotonic.Data(),NPhotonic,nbinsPhotonic, binlowPhotonic, binupPhotonic);
  fOutput->Add(fhPhotonicEle);
//=====================================================================================================================

//=====================================================================================================================
//Container - Track QA
if(kQA)
{
const int NQA = 7;
Int_t nbinsQATracks[NQA]	 	= {200 ,200 ,100, 20  , 20, 90, 80};
Double_t binlowQATracks[NQA] 	= {0.  ,0.  ,0. , 0.0 ,-1., 0.,-10};
Double_t binupQATracks[NQA]	= {100.,100.,100, 2*TMath::Pi(), 1., 9., 10};
fSparseQATracks = new THnSparseF("fSparseQATracks", "fSparseQATracks; p_{T}; p; E; #varphi; #eta; E/p; n_{#sigma TPC};",NQA, nbinsQATracks, binlowQATracks, binupQATracks);
fOutput->Add(fSparseQATracks);

//study of JET:pt, eta, phi
const int NJets = 4;
Int_t nbinsJet[NJets]		= { 300, 300,40 ,80	};
Double_t binlowJet[NJets]	= {-100,-100,-2.,0	};
Double_t binupJet[NJets]	= { 200,200., 2.,2*TMath::Pi()};
fSparseJet = new THnSparseF("fSparseJet", "fSparseJet; p_{T}; p_{T}^{Corr};#eta;#varphi;",NJets, nbinsJet, binlowJet, binupJet);
fOutput->Add(fSparseJet);
}

if(kGeneralSpectra)
{
if(fdebug>1) cout<<"DDG Outputs 4"<<endl;

  if(fReadMC)
  {
  }  


}
//======================================================================

//======================================================================
//Cluster Matching Histograms
if(kCheckClusterMatching) 
{
  fHistPtDEtaDPhiClusTrack = new TH3F("fHistPtDEtaDPhiClusTrack","fHistPtDEtaDPhiClusTrack;#it{p}_{T}^{clus};#Delta#eta;#Delta#varphi",100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiClusTrack);

  fHistPtDEtaDPhiTrackClus = new TH3F("fHistPtDEtaDPhiTrackClus","fHistPtDEtaDPhiTrackClus;#it{p}_{T}^{clus};#Delta#eta;#Delta#varphi",100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiTrackClus);
}
//======================================================================


if(kDetector)
{
  const int NDetectors = 4;
  Int_t nbinsDetectors[NDetectors]	 	= {160 , 7  , 2  , 2  };
  Double_t binlowDetectors[NDetectors] 	= { 0  ,-0.5,-0.5,-0.5};
  Double_t binupDetectors[NDetectors]	= {160., 6.5,1.5 , 1.5};
  fDetector = new THnSparseF("fDetector","fDetector; TPCNcls; ITSNcls; ITS1; ITS2;", NDetectors, nbinsDetectors, binlowDetectors, binupDetectors);
  fOutput->Add(fDetector);
}
//======================================================================

PostData(1, fOutput); // Post data for ALL output slots > 0 here.
//cout<<"EXITING CREATING OUTPUT OBJECTS")<<endl;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHFeJetCorrel::FillHistograms()
{
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::CheckClusTrackMatching()
{
  if(!fTracksCont || !fCaloClustersCont) return;

  Double_t deta = 999;
  Double_t dphi = 999;

  //Get closest cluster to track
  fTracksCont->ResetCurrentID();
  AliPicoTrack* PicoTrack = static_cast<AliPicoTrack*>(fTracksCont->GetNextAcceptParticle());
  while(PicoTrack)
  {
	AliVTrack* track = PicoTrack->GetTrack();

    //Get matched cluster
    Int_t emc1 = track->GetEMCALcluster();
    if(emc1>=0)
    {

      AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
      if(clusMatch)
      {
		AliPicoTrack::GetEtaPhiDiff(track, clusMatch, dphi, deta);

		fHistPtDEtaDPhiTrackClus->Fill(track->Pt(),deta,dphi);
      }
    }
    PicoTrack = static_cast<AliPicoTrack*>(fTracksCont->GetNextAcceptParticle());
  }

  //Get closest track to cluster
  fCaloClustersCont->ResetCurrentID();
  AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(); 
  while(cluster)
  {
    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);
    //Get matched track
    AliVTrack *mt = NULL;      
    AliAODCaloCluster *acl = dynamic_cast<AliAODCaloCluster*>(cluster);
    if(acl || (acl->GetNTracksMatched()>0) ) mt = static_cast<AliVTrack*>(acl->GetTrackMatched(0));
    if(mt)
    {
      AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
      fHistPtDEtaDPhiClusTrack->Fill(nPart.Pt(),deta,dphi);
    }
    cluster = fCaloClustersCont->GetNextAcceptCluster();
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHFeJetCorrel::Run()
{
// Execute analysis for current event:
// heavy flavor candidates association to MC truth

fdebug = -2;
//=========================================================================================

if(fdebug>1) cout<<"DDG Run 0"<<endl;

//========================================================================================
//Getting the Event  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!aod)
  {
    Printf("ERROR: AOD not available");
    PostData(1,fOutput);
    return kFALSE;
  }
  fNentries->Fill(0);
    
    //It must be done in the CUTS object for EVENTS
    //TString triggerClasses = InputEvent()->GetFiredTriggerClasses();
    //if( (! fReadMC) && (! triggerClasses.Contains( "EG1" )) ) return kFALSE; 
    //DONE IN THE CUT OBJECT

if(fdebug>1) cout<<"DDG Run 1"<<endl;

//Copying Deepa way to get the event
    AliVEvent *fVevent = dynamic_cast<AliVEvent*> (InputEvent());
  if(!fVevent)
  {
    Printf("ERROR: VEvent not available");
    PostData(1,fOutput);
    return kFALSE;
  }

//=========================================================================================

if(fdebug>1) cout<<"DDG Run 2"<<endl;

//tracks rejection
fNentries->Fill(5);

if(!fJetsCont)
{
	AliWarning("AliAnalysisTaskSEHFjets::Run: Jets Container not found!");
	PostData(1,fOutput);
	return kFALSE;
}
fNentries->Fill(6);
//======================================================================

if(fdebug>1) cout<<"DDG Run 3"<<endl;
  if(!fpidResp) SetupPIDresponse();
if(fdebug>1) cout<<"DDG Run 4"<<endl;
//======================================================================



//=========================================================================================

bool IsSelected = fCutsHFjets->IsEventSelected(fVevent);
int kRejected = fCutsHFjets->GetWhyRejection();
if(!IsSelected) fNRejected->Fill(kRejected);

  if(!IsSelected) //SelectedDDG )
  {
    AliWarning("Rejecting event ");
    PostData(1,fOutput);
    return kFALSE;
  }
if(fdebug>1) cout<<"DDG Run 5"<<endl;

  fNentries->Fill(1);

//=========================================================================================

//PileUp rejection
  fNentries->Fill(2);
  if(!fcandEleTPC)fcandEleTPC=new TArrayI(200);
  else fcandEleTPC->Reset(-1);

  fLastEleTPC=0;
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0)
  { 
    fNentries->Fill(3);
  }
  else
  {
	PostData(1,fOutput);
	return kFALSE;  
  }
if(fdebug>1) cout<<"DDG Run 6"<<endl;


//=========================================================================================

//======================================================================
  // MC information: Did not change anything inside here
  TClonesArray *arrayMC=0x0;
  AliAODMCHeader *aodmcHeader=0x0;
  Double_t vtxTrue[3];

  if(fReadMC)
  {
    // load MC particles
    
	//arrayMC = (TClonesArray*)aod->GetList()->FindObject(fMCParticlesName);
    arrayMC = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC)
    {
      Printf("AliAnalysisTaskSEHFjets::UserExec: MC particles branch not found!\n");
      PostData(1,fOutput);
      return kFALSE;
    }
    //arrayMC->SetClass("AliAODMCParticle");
    // load MC header
    aodmcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!aodmcHeader)
    {
      Printf("AliAnalysisTaskSEHFjets::UserExec: MC header branch not found!\n");
      PostData(1,fOutput);
      return kFALSE;
    }
    // MC primary vertex
    aodmcHeader->GetVertex(vtxTrue);
 
      //LOAD MC JETS
      /*
      fMCJetsCont = (AliJetContainer*)aod->GetList()->FindObject(fMCJetsBranch);
      if (! fMCJetsCont)
      {
          Printf("AliAnalysisTaskSEHFjets::UserExec: MC JETS branch not found!\n");          
          PostData(1,fOutput);
          return kFALSE;
      }
      */
  }
//======================================================================


//=========================================================================================

  //AliVVertex *vprimary = (AliVVertex*)aod->GetPrimaryVertex();
  Double_t primvtx[3],primcov[6];
  vtx1->GetXYZ(primvtx);
  vtx1->GetCovarianceMatrix(primcov);
  Int_t nPrimContr=vtx1->GetNContributors();
  Double_t chi2=vtx1->GetChi2();
  
  AliESDVertex* v1 = new AliESDVertex(primvtx,primcov,chi2,nPrimContr);


//______________________________________________________________________
//Preparing the Objetcs for the Analysis

//======================================================================
//Find Electrons Using TPC for Invariant Mass Study
FindTPCElectrons(3.5,-3.5);
//======================================================================
if(fdebug>1) cout<<"DDG Run 7"<<endl;

//======================================================================
// JET STUDY
JetStudy();
RandomCones(fJetsCont->GetJetRadius(), fJetsCont->GetRhoVal());
//if(fReadMC) MCJetStudy();
//======================================================================
//______________________________________________________________________

if(fdebug>1) cout<<"DDG Run 7.0"<<endl;

//======================================================================


if(fdebug>1) cout<<"DDG Run 8"<<endl;

//=========================================================================================
//MAIN LOOP
for(int itr = 0; itr < fTracksCont->GetNAcceptedParticles(); itr++)
{
	AliPicoTrack *PicoTrack = static_cast<AliPicoTrack*>(fTracksCont->GetAcceptParticle(itr));
	if(!PicoTrack) continue;

	AliVTrack* track = PicoTrack->GetTrack();
	if(!track) continue;
	
if(fdebug>1) cout<<"DDG Run 9"<<endl;

    if(!track)
    {
      printf("ERROR: Could not receive track %d\n", itr);
      continue;
    }

if(fdebug>1) cout<<"DDG Run 10"<<endl;
//======================================================================
//Monte Carlo Studies
    bool isMCele = kFALSE;
    double IsHeavy = -1;
    if(fReadMC)
    {
	  int Label = track->GetLabel();
	  if(fdebug>2) cout<<"DDG Run 10.1"<<endl;
      if( (Label>=0)&&(Label < arrayMC->GetEntriesFast()) )
      {
		  AliAODMCParticle* pele = (AliAODMCParticle*)arrayMC->At(Label);
		  if(fdebug>2) cout<<"DDG Run 10.2"<<endl;
		  if(pele)
		  {
			if(fdebug>2) cout<<"DDG Run 10.2.2"<<endl;

			if(fdebug>2) cout<<"DDG Run 10.3"<<endl;
			int PDG = pele->GetPdgCode();
			if(fdebug>2) cout<<"DDG Run 10.4"<<endl;
			bool isPhysPrim= pele->IsPhysicalPrimary();

			if( (TMath::Abs(PDG)==11)&&(isPhysPrim) )
			{
				isMCele=kTRUE;
			
				int MotherLabel = -1;
				AliAODMCParticle* OriginParton = GetMCPartonOrigin(arrayMC,pele, MotherLabel);
				if(OriginParton)
				{
					IsHeavy = OriginParton->GetPdgCode();
				}
				if( (abs(IsHeavy) > 5)||(abs(IsHeavy) < 1) ) IsHeavy = -1;
				if ( abs(IsHeavy) == 21) IsHeavy = 0; //gluons
					
			}
		    if(fdebug>2) cout<<"DDG Run 10.5"<<endl;

		  }
		  if(fdebug>2) cout<<"DDG Run 10.6"<<endl;
	  }
    }
    //if(isMCele)fhEleRejection->Fill(0,-1,nsigma);
    
    Double_t p=track->P();
    Double_t pt=track->Pt();


//======================================================================
//ELECTRONS SELECTION STARTS HERE

	//======================================================================
	//Check whether it is an electron candidate: TPC CUT ONLY FOR THE MOMENT: TOF to BE INCLUDED, 
	//do not reject here hadrons to allow QA checks of (E/p,nsigmaTPC,pt)
    Double_t nsigma=fpidResp->NumberOfSigmasTPC(track, AliPID::kElectron);
    int TPCPID = track->GetTPCsignalN();

    AliAODTrack* atrack = dynamic_cast<AliAODTrack*>(track);
    if(!atrack){continue;}
    int TPCNcls = atrack->GetTPCNcls();
    int ITSNcls = atrack->GetITSNcls();

    if(fdebug>1) cout<<"DDG Run 11"<<endl;
	//All Tracks
    fhTrackRejection->Fill(1,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(1,pt,nsigma);
	if(fdebug>1) cout<<"DDG Run 11.1"<<endl;


	//CHECK WHETHER THERE IS A EMCAL CLUSTER
    Int_t nClsId = track->GetEMCALcluster();

    if(fdebug>1) cout<<"DDG Run 11.2"<<endl;
    if(nClsId <0) {continue;}
    
if(fdebug>1) cout<<"DDG Run 11: With Clusters!"<<endl;
    fhTrackRejection->Fill(2,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(2,pt,nsigma);      
//======================================================================

//======================================================================
    // CHECK FILTER BIT and track ID
    if( (!(atrack->TestFilterBit(ffilterbit))) || (track->GetID()<0) ) {continue;}

    fhTrackRejection->Fill(3,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(3,pt,nsigma);
//======================================================================

//==============================================================================
//Cluster Cuts
	AliVCluster *cluster = fCaloClustersCont->GetCluster(nClsId);

    if(!cluster) {continue;}
    if(!cluster->IsEMCAL()) {continue;}

	fhTrackRejection->Fill(4,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(4,pt,nsigma);
//======================================================================

    double ClusterEnergy = cluster->E();
    double EoP = ClusterEnergy/p;

if(fdebug>1) cout<<"DDG Run 12"<<endl;

//LOOSE CUTS
//======================================================================
    if(TMath::Abs(cluster->GetTrackDx())>fPhiTrackClusterDistance
       || TMath::Abs(cluster->GetTrackDz())>fEtaTrackClusterDistance
       || cluster->GetNCells()<fminNcell) {continue;}

	fhTrackRejection->Fill(5,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(5,pt,nsigma);
//======================================================================

if(fdebug>1) cout<<"DDG Run 12.2"<<endl;

//======================================================================
    double M20 = cluster->GetM20();
    double M02 = cluster->GetM02();
    int NLM = -1;

if(fdebug>1) cout<<"DDG Run 12.3"<<endl;

    AliVCaloCells* cells = InputEvent()->GetEMCALCells();
if(fdebug>1) cout<<"DDG Run 12.4"<<endl;
//    AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
//if(fdebug>1) cout<<"DDG Run 12.5"<<endl;
//    if(cells && fCaloUtils) NLM = fCaloUtils->GetNumberOfLocalMaxima(cluster, cells);
if(fdebug>1) cout<<"DDG Run 12.6"<<endl;

    
    if( (M02>fmaxM02)||(M20>fmaxM20)||(NLM > 3) ) {continue;} 

    double PhiMatch = - 1;
    double EtaMatch = -1;
    AliPicoTrack::GetEtaPhiDiff(track, cluster, PhiMatch, EtaMatch);
    
    fhTrackRejection->Fill(6,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(6,pt,nsigma);
//==============================================================================


//==============================================================================
//Track Cuts: ETA CUT NEEDED    
    if( (pt<fminpt)||(track->Eta() > fEtaMax)||(track->Eta() < fEtaMin) ) {continue;}

    fhTrackRejection->Fill(7,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(7,pt,nsigma);
    
if(fdebug>1) cout<<"DDG Run 13"<<endl;
//==============================================================================

//==============================================================================
if(kQA)
{
	double poiv[7]={pt, p, ClusterEnergy,track->Phi(),track->Eta(), EoP, nsigma};
	fSparseQATracks->Fill(poiv);
}

//======================================================================
// !!!!!!!!!!!!!!!!!!!!!!!!!! PROBLEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//The Eletron Cuts are created here!
//Something is not working properly
fCutsElectron = new AliRDHFJetsCuts("fCutsElectron", "fCutsElectron");
//======================================================================
//======================================================================
//Verify Daughter
    if( !fCutsElectron->IsDaughterSelected(atrack,v1,fCutsElectron->GetTrackCuts())){continue;}

	fhTrackRejection->Fill(8,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(8,pt,nsigma);
if(fdebug>1) cout<<"DDG Run 14"<<endl;
//=======================================================================


//======================================================================
//Rejecting tracks outside the TPC region for electrons + Pions
	if( (nsigma<fsigmaTPCmin)||(nsigma>fsigmaTPCmax) ) continue;
	fhTrackRejection->Fill(9,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(9,pt,nsigma);
//======================================================================

//======================================================================
//Rejecting the tracks with E/p outside the electron + pions interval
	if(EoP<fMinEoverP||EoP>fMaxEoverP) continue;
    fhTrackRejection->Fill(10,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(10,pt,nsigma);
//======================================================================
if(fdebug>1) cout<<"DDG Run 14.1"<<endl;



if(fdebug>1) cout<<"DDG Run 14.5"<<endl;

	fhTrackRejection->Fill(11,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(11,pt,nsigma);
    
    // semi inclusive Electron candidate selected
    fhTrackRejection->Fill(12,pt,nsigma);
    if(isMCele)fhEleRejection->Fill(12,pt,nsigma);

    if(fCheckVzero)
    {
		// check if the electron is a vzero candidate
		for(Int_t vz=0;vz<=aod->GetNumberOfV0s();vz++)
		{
			AliAODv0 *v0=aod->GetV0(vz);
			for(Int_t vzd=0;vzd<v0->GetNDaughters();vzd++)
			{
			  AliAODTrack *aodvz=(AliAODTrack*)v0->GetDaughter(vzd);
			  if(aodvz->GetID()>0)
			  {
				if(aodvz->GetID()==track->GetID())
				{
				  fhTrackRejection->Fill(13,pt,nsigma);
				  if(isMCele)fhEleRejection->Fill(13,pt,nsigma);
				}
			  }
			}
		 }
    }


//=======================================================================
//Correlate electrons with jets
if(kAnalysis)
{
  double PointHFTrigger[16] = {track->Pt(),EoP,(double)TPCNcls,(double)TPCPID,nsigma,(double)ITSNcls,PhiMatch, EtaMatch,M20, M02, (double)NLM};
    fSparseHFTriggers->Fill(PointHFTrigger);
    
    for(int ithJet = 0; ithJet < fJetsCont->GetNAcceptedJets(); ithJet++)
    {
        AliEmcalJet* Emcaljet = fJetsCont->GetAcceptJet(ithJet); 
        if(!Emcaljet) continue;
        
        Double_t phiele = track->Phi();
        Double_t phijet= Emcaljet->Phi();
        Double_t deltaphi = PhiInterval(phijet - phiele);
        
        double etaele = track->Eta();
        double etajet = Emcaljet->Eta();
        double deltaeta = etajet - etaele;
        
        double R = RPhiEta(phijet, etajet, phiele, etaele);
        double InJet = -0.5;//= IsInsideJet(track, Jet);
        
        double PointHF[16] = {track->Pt(),EoP,(double)TPCNcls,(double)TPCPID,nsigma,(double)ITSNcls,PhiMatch, EtaMatch,M20, M02, (double)NLM,Emcaljet->PtSub(), deltaphi, deltaeta, R, InJet};
        fSparseHFSpectrum->Fill(PointHF);
        
        //======================================================================
        //Photonic Selection
        Bool_t isPhotonic=kFALSE;
        //Electron mass
        double Me = 0.0005109989;
        double Me2 = Me*Me;
        //NOW CHECK IF IT IS FROM CONVERSIONS
        for(Int_t kph=0;kph<fLastEleTPC;kph++)
        {
            Int_t trelenum=fcandEleTPC->At(kph);
            if(trelenum<0)continue;
            if(trelenum==itr)continue;
            
            //Using the VEvent
            AliPicoTrack *PicoTracktrEle = static_cast<AliPicoTrack*>(fTracksCont->GetAcceptParticle(trelenum));
            if(!PicoTracktrEle) continue;
            AliVTrack* trEle = PicoTracktrEle->GetTrack();
            if(!trEle) continue;
            
            if(fdebug>3) cout<<"DDG Run 14.2"<<endl;
            if (!trEle)
            {
                printf("ERROR: Could not receive track %d for Photonic calculation\n", trelenum);
                continue;
            }
            Double_t pxyz[3];
            track->PxPyPz(pxyz);
            Double_t pxyzB[3];
            trEle->PxPyPz(pxyzB);
            
            Double_t enA=TMath::Sqrt(Me2 + pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2]);
            Double_t enB=TMath::Sqrt(Me2 + pxyzB[0]*pxyzB[0]+pxyzB[1]*pxyzB[1]+pxyzB[2]*pxyzB[2]);
            //Double_t mass=TMath::Sqrt((enA+enB)*(enA+enB)-(pxyz[0]+pxyzB[0])*(pxyz[0]+pxyzB[0])-(pxyz[1]+pxyzB[1])*(pxyz[1]+pxyzB[1])-(pxyz[2]+pxyzB[2])*(pxyz[2]+pxyzB[2]));
            Double_t mass=TMath::Sqrt(2.*(Me2 + enA*enB-pxyz[0]*pxyzB[0]-pxyz[1]*pxyzB[1]-pxyz[2]*pxyzB[2])); // Using mass=sqrt(m2+m2+2 p p); NOT converted to keV
            
            //====================================================================================================
            //PARTNER VARIABLES
            double PtPartner = trEle->Pt();
            int TPCNclsPartner = trEle->GetTPCNcls();
            double nsigmaPartner = fpidResp->NumberOfSigmasTPC(trEle, AliPID::kElectron);
            double ThetaOpen = 0;
            double dd = TMath::Sqrt( (track->Phi() - trEle->Phi())*(track->Phi() - trEle->Phi()) + (track->Eta() - trEle->Eta())*(track->Eta() - trEle->Eta()) );
            if(dd > 0.0) ThetaOpen = 0.0;
            if(dd > 0.1) ThetaOpen = 0.1;
            if(dd > 0.2) ThetaOpen = 0.2;
            if(dd > 0.3) ThetaOpen = 0.3;
            if(dd > 0.4) ThetaOpen = 0.4;
            if(dd > 0.5) ThetaOpen = 0.5;
            double Sign = 0;
            if(trEle->Charge()*track->Charge() >  0.00001) Sign = +1;
            if(trEle->Charge()*track->Charge() < -0.00001) Sign = -1;
            //====================================================================================================
            
            
            if( (mass<fMassPhotonicCut) )
            {
                if(fdebug>3) cout<<"DDG Run 14.4"<<endl;
                double Point[22] = {track->Pt(),EoP,(double)TPCNcls,(double)TPCPID,nsigma,(double)ITSNcls,PhiMatch, EtaMatch,M20, M02, (double)NLM,Emcaljet->PtSub(), deltaphi, deltaeta, R, InJet, PtPartner, (double)TPCNclsPartner, nsigmaPartner, mass, ThetaOpen, Sign};
                fhPhotonicEle->Fill(Point);

            }
            
        }
        //=======================================================================

    }
}
if(fdebug>1) cout<<"DDG Run 15"<<endl;
//=======================================================================



//======================================================================
//Detector Minimum Requirements
if(kDetector)
{
	double Poid[4] = {(double)atrack->GetTPCNcls(), (double)atrack->GetITSNcls(), (double)atrack->HasPointOnITSLayer(0), (double)atrack->HasPointOnITSLayer(1)};
	fDetector->Fill(Poid);
}
//======================================================================


}
if(fdebug>1) cout<<"DDG Run 16"<<endl;
  delete v1;

//Check Cluster Track Matching
if(kCheckClusterMatching) CheckClusTrackMatching();

PostData(1, fOutput);


  return kFALSE;
  //return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::Terminate(Option_t *) 
{
	//cout<<"DDG10"<<endl;
  // Called once at the end of the analysis.
}

//---------------------------------------------------------------
AliAODMCParticle* AliAnalysisTaskEmcalHFeJetCorrel::IsMCJet(TClonesArray *arrayMC,const AliEmcalJet *jet, Double_t &contribution){// assignment of parton ID to jet
  // method by L. Feldkamp
  std::vector< int >           idx;
  std::vector< int >           idx2;
  std::vector< double >     weight;

  int counter =0;
  //int num = jet->GetRefTracks()->GetEntries();
  int num = jet->GetNumberOfTracks();
  for(int k=0;k<num;++k){

    AliVParticle* track = fTracksCont->GetAcceptParticle(jet->TrackAt(k));
//Verifing wether this method gets the RIGHT particles
    //TestHisto(jet, fTracksCont);
    
    //DDG upperbond for the label
    if( (track->GetLabel() >=0)&&(track->GetLabel() < arrayMC->GetEntriesFast()) )
    {
      AliAODMCParticle* part =  (AliAODMCParticle*)  arrayMC->At(track->GetLabel());
      if(!part)continue;

      int label =0 ;
      AliAODMCParticle* motherParton=GetMCPartonOrigin(arrayMC,part, label);
      if (!motherParton) ;//Printf("no mother");
      else {
	counter++;
	idx.push_back(label);                       //! Label  of Mother
	idx2.push_back(label);        
	weight.push_back(track->Pt());  //! Weight : P_t trak /  P_t jet ... the latter used at the end
      }
    }///END LOOP OVER REFERENCED TRACKS   
  }
  //! Remove duplicate indices for counting
  sort( idx2.begin(), idx2.end() );
  idx2.erase( unique( idx2.begin(), idx2.end() ), idx2.end() );
  if (idx2.size() == 0) return 0x0;
  Double_t* arrayOfWeights = new Double_t [(int)idx2.size()];
  for(unsigned int ii=0;ii<idx2.size();ii++) arrayOfWeights[ii]=0;

  for (unsigned int idxloop =0 ;idxloop<idx2.size();idxloop++){
    for (unsigned int z=0; z< idx.size() ; ++z){
      int     a = idx.at(z);
      double w = weight.at(z);
      if(a == idx2.at(idxloop))    arrayOfWeights[idxloop] += w;
    }
  }
  
  int winner = -1;
  double c=-1.;
  for (unsigned int idxloop =0 ;idxloop<idx2.size();idxloop++){
    if(c < arrayOfWeights[idxloop]){
      winner =idxloop; 
      c=arrayOfWeights[idxloop];
    }
  }
  
  AliAODMCParticle *parton = 0x0;
  if( (idx.at(winner) > 0)&&(idx.at(winner) < arrayMC->GetEntriesFast()) )
  {
	parton=(AliAODMCParticle*)arrayMC->At(idx.at(winner));
	contribution = arrayOfWeights[winner]/jet->PtSub();
  }
  delete[] arrayOfWeights;

  return parton;  
}

//---------------------------------------------------------------
AliAODMCParticle *AliAnalysisTaskEmcalHFeJetCorrel::GetMCPartonOrigin(TClonesArray* &arrayMC,AliAODMCParticle* &p, Int_t &idx)
{  //Follows chain of track mothers until q/g or idx = -1	
  AliAODMCParticle *p2=0x0;
  Int_t mlabel = TMath::Abs(p->GetMother()) ; 
//  Double_t pz=0.;
  while( (mlabel > 1)&&(mlabel < arrayMC->GetEntriesFast()) )
  {
    p2 = (AliAODMCParticle*)arrayMC->At(mlabel);
//    pz=TMath::Abs(p2->Pz());
    //printf("Mother label %d, pdg %d, pz %f\n",mlabel,p2->PdgCode(),pz);
    int PDG = abs(p2->PdgCode());
    if( (PDG == 21) || ( (PDG > 0) && (PDG <6) ) )
      {
		idx = mlabel; 
		return p2;
      }
    mlabel = TMath::Abs(p2->GetMother()); 
  }
  idx=-1;
  return p2;

} 
//======================================================================

//____________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::FillCorrelationHisto(const AliVParticle *part, const double Ep, const Int_t idlab, const AliEmcalJet *jet,Double_t nsigma, THnSparseF* &h){
// Fill histo with ele properties and deltaphi
  
  int num = jet->GetNumberOfTracks();
  //Printf("numtrackjet %d/n",num);
  Double_t flag=0.5;
  Int_t trlab;
  
  for(int w=0;w<num;++w)
  {
    AliVParticle* trjet = fTracksCont->GetAcceptParticle(jet->TrackAt(w));
	//Verifing whether this method gets the RIGHT particles
    //TestHisto(jet, fTracksCont);

      trlab=trjet->GetLabel();
      if(idlab==trlab)
      {
		//if ele is a track in the jet, flag=1.5
		flag=1.5;
		break;
      }
    
  }
  Double_t phiele = part->Phi();
  Double_t phijet= jet->Phi();
  Double_t deltaphi = PhiInterval(phijet - phiele);

  double etaele = part->Eta();
  double etajet = jet->Eta();
  double deltaeta = etajet - etaele;

  double R = RPhiEta(phijet, etajet, phiele, etaele);

  Double_t pois[8]={part->Pt(), Ep, nsigma,jet->PtSub(), deltaphi, deltaeta, R, flag};
  h->Fill(pois);
}
//======================================================================

//____________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::FillMCCorrelationHisto(const AliVParticle *part, const double Ep, const Int_t idlab, const AliEmcalJet *jet,Double_t nsigma, THnSparseF* &h, double IsHeavy, double JetHeavy, double Contribution){
// Fill histo with ele properties and deltaphi
  
  int num = jet->GetNumberOfTracks();
  //Printf("numtrackjet %d/n",num);
  Double_t flag=0.5;
  Int_t trlab;
  
  for(int w=0;w<num;++w)
  {
    AliVParticle* trjet = fTracksCont->GetAcceptParticle(jet->TrackAt(w));
	//Verifing whether this method gets the RIGHT particles
    //TestHisto(jet, fTracksCont);

      trlab=trjet->GetLabel();
      if(idlab==trlab)
      {
		//if ele is a track in the jet, flag=1.5
		flag=1.5;
		break;
      }
    
  }
  Double_t phiele = part->Phi();
  Double_t phijet= jet->Phi();
  Double_t deltaphi = PhiInterval(phijet - phiele);

  double etaele = part->Eta();
  double etajet = jet->Eta();
  double deltaeta = etajet - etaele;

  double R = RPhiEta(phijet, etajet, phiele, etaele);

  Double_t pois[11]={part->Pt(), Ep, nsigma,jet->PtSub(), deltaphi, deltaeta, R, flag, IsHeavy, JetHeavy, Contribution};
  h->Fill(pois);


}
//======================================================================


//____________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::FillJetHisto(const AliEmcalJet *jet){
// Fill histo with jet properties

  Double_t pais[4]={jet->Pt(), jet->PtSub(), jet->Eta(), jet->Phi()};

  fSparseJet->Fill(pais);

}
//======================================================================

//____________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::FillMCJetHisto(const AliEmcalJet *jet){
    // Fill histo with jet properties
    
    Double_t pais[3]={jet->Pt(), jet->Eta(), jet->Phi()};
    
    //fSparseMCJet->Fill(pais);
    
}
//======================================================================


//____________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::FillMassHisto(const AliVParticle *part, const AliEmcalJet *jet, double pgamma, double mass, THnSparseF* &h){
// Fill histo with ele properties and deltaphi
  Double_t phiele = part->Phi();
  Double_t phijet= jet->Phi();
  Double_t deltaphi = PhiInterval(phijet - phiele);

  double etaele = part->Eta();
  double etajet = jet->Eta();
  double deltaeta = etajet - etaele;
  
  double R = RPhiEta(phijet, etajet, phiele, etaele);
            
    
//    double Point[22] = {part->Pt(),nEid,TPCNcls,TPCPID,nsigma,ITScls,PhiMatch, EtaMatch,M20, M02, NLM,jet->PtSub(), deltaphi, deltaeta, R, InJet, Partner->Pt(), TPCNclsPartner, nsigmaPartner, mass, ThetaOpen, Sign};

    Double_t pois[7]={part->Pt(), pgamma, mass, deltaphi, deltaeta, R, jet->PtSub()};
  
  h->Fill(pois);
}
//======================================================================

//____________________________________________________________
Int_t AliAnalysisTaskEmcalHFeJetCorrel::TagJetMC(){// METHOD NOT IMPLEMENTED YET
  return -1;
}
//======================================================================

//____________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::SetupPIDresponse(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler(); 
  fpidResp=inputHandler->GetPIDResponse();
  if(!fpidResp)AliFatal("No PID response could be set");
}
//======================================================================


//____________________________________________________________
/*void AliAnalysisTaskEmcalHFeJetCorrel::TestHisto(const AliEmcalJet *jet,AliParticleContainer* TracksCont)
{// Fill histo with constituents properties

for(int k  =0 ; k < jet->GetNumberOfTracks(); k++)
{
	AliVParticle* track = TracksCont->GetParticle(jet->TrackAt(k));
	if(!track) continue;
	double DeltaPtJet = track->Pt() - jet->PtSub();
	double DeltaEtaJet = track->Eta() - jet->Eta();
	double DeltaPhiJet = track->Phi() - jet->Phi();
	Double_t pais[3]={DeltaPtJet,DeltaEtaJet,DeltaPhiJet};

	fSparseConstituents->Fill(pais);
}
}
*/
//======================================================================


//______________________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::TestContainerEvent(AliParticleContainer* &Cont, AliClusterContainer* &CaloClustersCont, TClonesArray* &Tracks_tender, TClonesArray* &CaloClusters_tender, AliVEvent* &Event, THnSparseF* &fTHN, THnSparseF* &fTHN2, THnSparseF* &fTHN3)
{
double ContainerEvent=1;

for(int itrx = 0; itrx < fTracksCont->GetNAcceptedParticles(); itrx++)
{
	AliPicoTrack *PicoTrack = static_cast<AliPicoTrack*>(Cont->GetAcceptParticle(itrx));
	AliVTrack* ContTrack = PicoTrack->GetTrack();
	if(!ContTrack) continue;

	double Label = ContTrack->GetLabel();
	double Pt = ContTrack->Pt();
	double Eta = ContTrack->Eta();
	double Phi = ContTrack->Phi();
	int ID = ContTrack->GetEMCALcluster();
	double Energy = -1;
	if(ID >= 0)
	{
		AliVCluster *cluster = CaloClustersCont->GetCluster(ID);
		if(cluster && cluster->IsEMCAL()) Energy = cluster->E();
	}
	Double_t pais[6]={Energy,Pt,Eta,Phi, Label, ContainerEvent};
	fTHN->Fill(pais);
}

ContainerEvent=2;
for(Int_t itr = 0; itr < Event->GetNumberOfTracks(); itr++)
{
	AliVTrack* VTrack = static_cast<AliVTrack*>(Event->GetTrack(itr));

	if(!VTrack) continue;
	
	double Label2 = VTrack->GetLabel();
	double Pt2 = VTrack->Pt();
	double Eta2 = VTrack->Eta();
	double Phi2 = VTrack->Phi();
	double Energy2 = -1;
	int ID2 = VTrack->GetEMCALcluster();
	if(ID2 >= 0)
	{
		AliVCluster *cluster2 = Event->GetCaloCluster(ID2);
		if(cluster2 && cluster2->IsEMCAL()) Energy2 = cluster2->E();
	}
	
	Double_t pais2[6]={Energy2,Pt2,Eta2,Phi2, Label2, ContainerEvent};
	fTHN2->Fill(pais2);
}

ContainerEvent=3;
for(int itr = 0; itr < Tracks_tender->GetEntries(); itr++)
{
	AliVTrack* VTrack = static_cast<AliVTrack*>(Tracks_tender->At(itr));

	if(!VTrack) continue;
	
	double Label3 = VTrack->GetLabel();
	double Pt3 = VTrack->Pt();
	double Eta3 = VTrack->Eta();
	double Phi3 = VTrack->Phi();
	double Energy3 = -1;
	int ID3 = VTrack->GetEMCALcluster();
	if(ID3 >= 0)
	{
		AliVCluster *cluster3 = static_cast<AliVCluster*>(CaloClusters_tender->At(ID3));
		if(cluster3 && cluster3->IsEMCAL()) Energy3 = cluster3->E();
	}
	
	Double_t pais3[6]={Energy3,Pt3,Eta3,Phi3, Label3, ContainerEvent};
	fTHN3->Fill(pais3);
}

}
//======================================================================

//____________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::ClusterHisto(AliClusterContainer* &Cont, THnSparseF* &fTHN)
{

	for(int ithClus = 0; ithClus < Cont->GetNAcceptedClusters(); ithClus++)
	{
		AliVCluster *cluster = Cont->GetAcceptCluster(ithClus);
		if(!cluster) continue;

		double Energy = cluster->E();
		double EMCal = cluster->IsEMCAL();
		double NMatch = cluster->GetNTracksMatched();	
		Double_t pais[3]={NMatch, Energy, EMCal};
		fTHN->Fill(pais);
	}
}
//======================================================================

//____________________________________________________________
/*void FillRejection(TString Label, double pT, double NSigma, TH3F *fRejection)
{
	//How to set the label? Maybe using another function to save them if not set yet
	kNRejection++;
	(fRejection->GetAxis(kNREjection))->SetLabel();
	fRejection->Fill(kNRejection,pT,NSigma)
}
*/

//______________________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::FindTPCElectrons(double TPCMax, double TPCMin)
{//FIND TPC ELECTRONS HERE, OUTSIDE THE MAIN LOOP

for(int itrTPC = 0; itrTPC < fTracksCont->GetNAcceptedParticles(); itrTPC++)
{
	AliPicoTrack *PicoTrackTPC = static_cast<AliPicoTrack*>(fTracksCont->GetAcceptParticle(itrTPC));
	if(!PicoTrackTPC) continue;

	AliVTrack* trackTPC = PicoTrackTPC->GetTrack();
	if(!trackTPC) continue;

// Check whether it is an electron candidate: TPC CUT ONLY FOR THE MOMENT: TOF to BE INCLUDED, 
    Double_t nsigmaTPC=fpidResp->NumberOfSigmasTPC(trackTPC, AliPID::kElectron);
    if( (nsigmaTPC < TPCMax)&&(nsigmaTPC >TPCMin) )
    {
		fcandEleTPC->AddAt(itrTPC,fLastEleTPC);
		fLastEleTPC++;
	}
    if(fLastEleTPC > fcandEleTPC->GetSize()-1)fcandEleTPC->Set(fcandEleTPC->GetSize()+50);

}
}
//======================================================================


//______________________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::JetStudy()
{
  for(int ithJet = 0; ithJet < fJetsCont->GetNAcceptedJets(); ithJet++)
  {
	  AliEmcalJet *EmCalJet = fJetsCont->GetAcceptJet(ithJet);
	  if(!EmCalJet) continue;

	  //Sets PtSub: pT after background subtraction
	  EmCalJet->SetPtSub( EmCalJet->Pt() - fJetsCont->GetRhoVal() * EmCalJet->Area() );
      if(fdebug > 1) cout<<"Jet pT = "<<EmCalJet->PtSub()<<endl;
	  if(kGeneralSpectra) FillJetHisto(EmCalJet);
  }
}
//======================================================================

/*
//______________________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::MCJetStudy()
{
    if(fdebug>1) cout<<"DDG MC Jets Container 0"<<endl;
    if(fdebug>1) cout<<"DDG MC Jets Number = "<<fMCJetsCont->GetNAcceptedJets()<<endl;
    for(int ithJet = 0; ithJet < fMCJetsCont->GetNAcceptedJets(); ithJet++)
    {
        if(fdebug>1) cout<<"DDG MC Jets Container 1.1"<<endl;
        AliEmcalJet *EmCalJet = fMCJetsCont->GetAcceptJet(ithJet);
        if(!EmCalJet) continue;
        
        FillMCJetHisto(EmCalJet);
    }
    if(fdebug>1) cout<<"DDG MC Jets Container 2"<<endl;
}
//======================================================================
*/

//______________________________________________________________________
void AliAnalysisTaskEmcalHFeJetCorrel::RandomCones(double RParameter, double Rho, int TMax)
{
if(fdebug>1) cout<<"DDG RandomCones 0"<<endl;

TRandom *Random = new TRandom3();
Random->SetSeed(time(0));

for(int T = 0; T < TMax; T++)
{
int Constituents = 0;
//while(Constituents == 0)
{
double RandomPhi = Random->Uniform(0, 2*TMath::Pi());

double RandomEta = Random->Uniform(-1+RParameter,1-RParameter);

double ConepT = 0;
double ConePx = 0;
double ConePy = 0;
double ConePz = 0;
double ConePhi = 0;
double ConeEta = 0;
double ConepT2 = 0;


if(fdebug>1) cout<<"DDG RandomCones 1"<<endl;
	for(int itrx = 0; itrx < fTracksCont->GetNAcceptedParticles(); itrx++)
	{
		if(fdebug>1) cout<<"DDG RandomCones 1.1"<<endl;
        AliPicoTrack *PicoTrack = static_cast<AliPicoTrack*>(fTracksCont->GetAcceptParticle(itrx));
		if(!PicoTrack) continue;
		AliVTrack* ContTrack = PicoTrack->GetTrack();
		if(!ContTrack) continue;

		double Label = ContTrack->GetLabel();
		double Pt = ContTrack->Pt();
		double Eta = ContTrack->Eta();
		double Phi = ContTrack->Phi();
		double Px = ContTrack->Px();
		double Py = ContTrack->Py();
		double Pz = ContTrack->Pz();
		
		if(fdebug>1) cout<<"DDG RandomCones 1.2"<<endl;
        double deltaphi = PhiInterval(RandomPhi - Phi);
        if(deltaphi > TMath::Pi()) deltaphi = 2.0*TMath::Pi() - deltaphi;
        if(deltaphi < 0) deltaphi = -1.0*deltaphi;
		
        double DeltaR = deltaphi*deltaphi + (RandomEta - Eta)*(RandomEta - Eta);
		
		if(DeltaR < RParameter*RParameter)
		{
			ConepT = ConepT + Pt;
			ConePx = ConePx + Px;
			ConePy = ConePy + Py;
			ConePz = ConePz + Pz;
			ConePhi = ConePhi + Phi*Pt;
			ConeEta = ConeEta + Eta*Pt;
            Constituents++;
		}
		if(fdebug>1) cout<<"DDG RandomCones 1.3"<<endl;
	}

if(fdebug>1) cout<<"DDG RandomCones 2"<<endl;

//if(Constituents>0)
{
    ConepT = ConepT - Rho*TMath::Pi()*RParameter*RParameter;
    ConepT2 = TMath::Sqrt(ConePx*ConePx + ConePy*ConePy) - Rho*TMath::Pi()*RParameter*RParameter;
    ConePhi = ConePhi/ConepT;
    ConePhi = PhiInterval(ConePhi,0,2*TMath::Pi());
    ConeEta = ConeEta/ConepT;
    if(fdebug>1) cout<<"DDG RandomCones 3"<<endl;
    Double_t point[6]={(double)Constituents, RandomEta, RandomPhi, ConepT, ConeEta, ConePhi};
    fRandomCones->Fill(point);
}
    
}
}
}
//======================================================================

//______________________________________________________________________
double AliAnalysisTaskEmcalHFeJetCorrel::PhiInterval(double Phi, double Low, double Up)
{
    while( (Phi < Low)||(Phi > Up) )
    {
        if(Phi< Low){Phi= Phi + 2.*TMath::Pi();}
        if(Phi > Up) {Phi= Phi - 2.*TMath::Pi();}
    }
    
    return Phi;
}
//======================================================================

//______________________________________________________________________
double AliAnalysisTaskEmcalHFeJetCorrel::RPhiEta(double PhiRef, double EtaRef, double Phi, double Eta)
{
    double DeltaPhi = PhiInterval(PhiRef - Phi + TMath::Pi(), 0, 2*TMath::Pi());
    if(DeltaPhi > TMath::Pi()) DeltaPhi =  2*TMath::Pi() - DeltaPhi;
    double R = TMath::Sqrt(DeltaPhi*DeltaPhi + (EtaRef + Eta)*(EtaRef + Eta));
    return R;
}
//======================================================================
