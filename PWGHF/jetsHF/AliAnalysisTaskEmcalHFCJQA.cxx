//////////////////////////////////////////////////////
//													//
// QA Task for Analysing pp2012 MC anchored to Data	//
//													//
// Authors: Andrea Rossi and Diógenes D. Gimenez	//
// August, 2015										//
//////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

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

#include "AliAnalysisTaskEmcalHFCJQA.h"

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

#include "AliAnalysisTaskEmcalHFCJQA.h"
#include "AliRDHFJetsCuts.h"
#include "AliAnalysisVertexingHF.h"

#include "AliPIDResponse.h"
#include "AliVParticle.h"
#include "AliVCluster.h"
#include <TArrayI.h>


ClassImp(AliAnalysisTaskEmcalHFCJQA)

using std::cout;
using std::endl;
//________________________________________________________________________
AliAnalysisTaskEmcalHFCJQA::AliAnalysisTaskEmcalHFCJQA() : 
AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalHFCJQA", kTRUE),

//======================================================================
//Containers
fJetsCont(0),
fTracksCont(0),
fCaloClustersCont(0),
//======================================================================

//======================================================================
//Flags
fReadMC(0),
fDebug(0),
//======================================================================

//======================================================================
//Variables
fpidResp(0x0),
fCuts(0),
ffilterbit(AliAODTrack::kTrkTPCOnly),
fKeepTrackNegID(kTRUE),
//======================================================================

//======================================================================
//Histograms
fhEventCounter(0),
fhTriggerCounter(0),
fhTriggerMaskCounter(0),
fhTriggerBitCounter(0),
fClustersEnergydistribution(0),
fEventsThreshold(0),
fSparseRecoJets(0),
fhSparseFilterMask(0),
fhSparseFilterMaskPico(0),
fhSparseFilterMaskTrackAcc(0),
fhSparseFilterMaskTrackAccPico(0),
fhSparseFilterMaskImpPar(0),
fhSparseFilterMaskImpParPico(0),
fhImpParResolITSsel(0),
fhImpParResolITSselGoodTracks(0),
fhnSigmaTPCTOFEle(0),
fhnSigmaTPCTOFPion(0),
fhnSigmaTPCTOFKaon(0),
fhnSigmaTPCTOFProton(0),
fhSparseEoverPeleTPC(0),
fhSparseShowShapeEleTPC(0),
fhTrackEMCal(0),
fhptSpectrum(0),
fhptSpectrumPico(0),
fCheckPico(kTRUE)
//======================================================================

{
	// Default constructor.
	SetMakeGeneralHistograms(kFALSE);
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalHFCJQA::AliAnalysisTaskEmcalHFCJQA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),

//======================================================================
//Containers
fJetsCont(0),
fTracksCont(0),
fCaloClustersCont(0),
//======================================================================

//======================================================================
//Flags
fReadMC(0),
fDebug(0),
//======================================================================

//======================================================================
//Variables
fpidResp(0x0),
fCuts(0),
ffilterbit(AliAODTrack::kTrkTPCOnly),
fKeepTrackNegID(kTRUE),
//======================================================================

//======================================================================
//Histograms
fhEventCounter(0),
fhTriggerCounter(0),
fhTriggerMaskCounter(0),
fhTriggerBitCounter(0),
fClustersEnergydistribution(0),
fEventsThreshold(0),
fSparseRecoJets(0),
fhSparseFilterMask(0),
fhSparseFilterMaskPico(0),
fhSparseFilterMaskTrackAcc(0),
fhSparseFilterMaskTrackAccPico(0),
fhSparseFilterMaskImpPar(0),
fhSparseFilterMaskImpParPico(0),
fhImpParResolITSsel(0),
fhImpParResolITSselGoodTracks(0),
fhnSigmaTPCTOFEle(0),
fhnSigmaTPCTOFPion(0),
fhnSigmaTPCTOFKaon(0),
fhnSigmaTPCTOFProton(0),
fhSparseEoverPeleTPC(0),
fhSparseShowShapeEleTPC(0),
  fhTrackEMCal(0),
  fhptSpectrum(0),
fhptSpectrumPico(0),
fCheckPico(kTRUE)
//======================================================================

{
	// Standard constructor.
	SetMakeGeneralHistograms(kFALSE);
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalHFCJQA::~AliAnalysisTaskEmcalHFCJQA()
{
  // Destructor.

  if(fJetsCont)						delete fJetsCont;
  if(fTracksCont)					delete fTracksCont;
  if(fCaloClustersCont)				delete fCaloClustersCont;
  if(fpidResp)						delete fpidResp;
  if(fCuts)							delete fCuts;
  if(fhEventCounter)				delete fhEventCounter;
  if(fhTriggerMaskCounter)			delete fhTriggerMaskCounter;
  if(fhTriggerBitCounter)			delete fhTriggerBitCounter;
  if(fhTriggerCounter)				delete fhTriggerCounter;
  if(fClustersEnergydistribution)   delete fClustersEnergydistribution;
  if(fEventsThreshold)				delete fEventsThreshold;
  if(fSparseRecoJets)				delete fSparseRecoJets;
  if(fhSparseFilterMask)			delete fhSparseFilterMask;
  if(fhSparseFilterMaskPico)		delete fhSparseFilterMaskPico;
  if(fhSparseFilterMaskTrackAcc)	delete fhSparseFilterMaskTrackAcc;
  if(fhSparseFilterMaskTrackAccPico)delete fhSparseFilterMaskTrackAccPico;
  if(fhSparseFilterMaskImpPar)		delete fhSparseFilterMaskImpPar;
  if(fhSparseFilterMaskImpParPico)	delete fhSparseFilterMaskImpParPico;
  if(fhImpParResolITSsel)			delete fhImpParResolITSsel;
  if(fhImpParResolITSselGoodTracks)	delete fhImpParResolITSselGoodTracks;
  if(fhnSigmaTPCTOFEle) 			delete fhnSigmaTPCTOFEle;
  if(fhnSigmaTPCTOFPion)			delete fhnSigmaTPCTOFPion;
  if(fhnSigmaTPCTOFKaon)			delete fhnSigmaTPCTOFKaon;
  if(fhnSigmaTPCTOFProton)			delete fhnSigmaTPCTOFProton;
  if(fhSparseEoverPeleTPC)			delete fhSparseEoverPeleTPC;
  if(fhSparseShowShapeEleTPC)		delete fhSparseShowShapeEleTPC;
  if(fhTrackEMCal)					delete fhTrackEMCal;
  if(fhptSpectrum)					delete fhptSpectrum;
  if(fhptSpectrumPico)              delete fhptSpectrumPico;

}

//________________________________________________________________________
void AliAnalysisTaskEmcalHFCJQA::UserCreateOutputObjects()
{
  // Create user output.
//Printf("STARTING CREATION OF OUTPUT OBJECTS");
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

//======================================================================
//Get the Containers
  fJetsCont           = GetJetContainer(0);
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
//======================================================================

//--------------------------
double pi = TMath::Pi();
//--------------------------

//== Histograms ==//

//  ########### DEFINE THE EVENT COUNTER ############
fhEventCounter=new TH1F("fhEventCounter","Counter of event selected",20,-0.5,19.5);
fhEventCounter->GetXaxis()->SetBinLabel(1,"Events analyzed");
fhEventCounter->GetXaxis()->SetBinLabel(2,"Event selected");
fhEventCounter->GetXaxis()->SetBinLabel(3,"Jet array present");
fhEventCounter->GetXaxis()->SetBinLabel(4,"Vtx Track Ncont");
fOutput->Add(fhEventCounter);
//======================================================================

//  ########### DEFINE THE EVENT TRIGGER ############
    fhTriggerCounter=new TH1F("fhTriggerCounter","Counter of Trigger Events",10,-0.5,9.5);
    fhTriggerCounter->GetXaxis()->SetBinLabel(1,"ALL");
    fhTriggerCounter->GetXaxis()->SetBinLabel(2,"EMC7");
    fhTriggerCounter->GetXaxis()->SetBinLabel(3,"EG1");
    fhTriggerCounter->GetXaxis()->SetBinLabel(4,"EG2");
    fhTriggerCounter->GetXaxis()->SetBinLabel(5,"EJ1");
    fhTriggerCounter->GetXaxis()->SetBinLabel(6,"EJ2");
    fhTriggerCounter->GetXaxis()->SetBinLabel(7,"EGA");
    fhTriggerCounter->GetXaxis()->SetBinLabel(8,"EJE");
    fOutput->Add(fhTriggerCounter);
//======================================================================
    
    //  ########### DEFINE THE EVENT TRIGGER MASK ############
    fhTriggerMaskCounter=new TH1F("fhTriggerMaskCounter","Counter of Trigger Masks Events",10,-0.5,9.5);
    fhTriggerMaskCounter->GetXaxis()->SetBinLabel(1,"ALL");
    fhTriggerMaskCounter->GetXaxis()->SetBinLabel(2,"kMC7");
    fhTriggerMaskCounter->GetXaxis()->SetBinLabel(3,"kEMCEGA");
    fhTriggerMaskCounter->GetXaxis()->SetBinLabel(4,"kEMCEJE");
    fhTriggerMaskCounter->GetXaxis()->SetBinLabel(5,"MC7+EGA");
    fhTriggerMaskCounter->GetXaxis()->SetBinLabel(6,"MC7+EJE");
    fhTriggerMaskCounter->GetXaxis()->SetBinLabel(7,"EGA+EJE");
    
    fOutput->Add(fhTriggerMaskCounter);
    
    //  ########### DEFINE THE EVENT TRIGGER MASK ############
    fhTriggerBitCounter=new TH1F("fhTriggerBitCounter","Counter of Trigger Bit Events",13,-0.5,12.5);
    
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(1,"ALL");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(2,"L0");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(3,"EGA1");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(4,"EGA2");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(5,"EJE1");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(6,"EJE2");
    
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(7,"EGA12");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(8,"EJE12");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(9,"EGA1EJE1");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(10,"EGA1EJE2");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(11,"EGA2EJE1");
    fhTriggerBitCounter->GetXaxis()->SetBinLabel(12,"EGA2EJE2");
    
    fOutput->Add(fhTriggerBitCounter);

    //  ########### DEFINE THE TRIGGER MASK ENERGY DISTRIBUTION ############    
    fClustersEnergydistribution = new TH2F("fClustersEnergydistribution","Counter of Trigger Masks Energy Distributions; Bit; Energy",12,-0.5,11.5,300,0,150);
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(1,"ALL");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(2,"L0");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(3,"EGA1");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(4,"EGA2");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(5,"EJE1");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(6,"EJE2");
    
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(7,"EGA12");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(8,"EJE12");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(9,"EGA1EJE1");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(10,"EGA1EJE2");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(11,"EGA2EJE1");
    fClustersEnergydistribution->GetXaxis()->SetBinLabel(12,"EGA2EJE2");

    fOutput->Add(fClustersEnergydistribution);
//======================================================================


//======================================================================
fEventsThreshold=new TH1F("fhTriggerThreshold","Counter of Trigger Thrseholds",10,-0.5,9.5);
fEventsThreshold->GetXaxis()->SetBinLabel(1,"ALL");
fEventsThreshold->GetXaxis()->SetBinLabel(2,"2GeV");
fEventsThreshold->GetXaxis()->SetBinLabel(3,"4GeV");
fEventsThreshold->GetXaxis()->SetBinLabel(4,"5GeV");
fEventsThreshold->GetXaxis()->SetBinLabel(5,"6GeV");
fEventsThreshold->GetXaxis()->SetBinLabel(6,"8GeV");
fEventsThreshold->GetXaxis()->SetBinLabel(7,"10GeV");
    
fOutput->Add(fEventsThreshold);
//======================================================================


Int_t nbinsRecoJets[8]			= { 50, 20,20,   20,   5,   5,10, 60};
Double_t binlowRecoJets[8]		= { 5.,-1.,pi, 0.99,  0.,-0.5, 0, 0.};
Double_t binupRecoJets[8]		= {55., 1.,pi,20.99,4.99, 4.5,2.,60.};
fSparseRecoJets = new THnSparseF("fSparseRecoJets","fSparseRecoJets;jetpt;eta;phi;ntrks;nEle;parton;partContr;ptPart;",8,nbinsRecoJets,binlowRecoJets,binupRecoJets);
fOutput->Add(fSparseRecoJets);

// Num axes: 10 filter bits + ID + TPCrefit,ITSrefit,bothTPCandITSrefit + kAny,kFirst,kBoth+ 20Nclust TPC+ 20 NTPC crossed padRows + 20 DCA
Int_t nbinsFilterMask[16]		={   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   3,   3, 20, 20,  70};
Double_t binlowFilterMask[16]	={-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,  0,  0,-3.5};
Double_t binupFilterMask[16]	={ 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2.5, 2.5,160,160, 3.5};
fhSparseFilterMask = new THnSparseF("fhSparseFilterMask","fhSparseFilterMask;fb0;fb1;fb2;fb3;fb4;fb5;fb6;fb7;fb8;fb9;trackID;refitting;SPD;NTPCclust;NTPCcrossRows;DCA",16,nbinsFilterMask,binlowFilterMask,binupFilterMask);
fOutput->Add(fhSparseFilterMask);

fhSparseFilterMaskPico = new THnSparseF("fhSparseFilterMaskPico","fhSparseFilterMaskPico;fb0;fb1;fb2;fb3;fb4;fb5;fb6;fb7;fb8;fb9;trackID;refitting;SPD;NTPCclust;NTPCcrossRows;DCA",16,nbinsFilterMask,binlowFilterMask,binupFilterMask);
fOutput->Add(fhSparseFilterMaskPico);


// Num axes: 10 filter bits + ID*isSelected 5 + kAny,kFirst,kBoth + phi+ eta +pt 
Int_t nbinsFilterMaskTrackAcc[15]		={   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   5,   3,  36,  30, 30};
Double_t binlowFilterMaskTrackAcc[15]	={-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-2.5,-0.5,  0.,-1.5, 0.};
Double_t binupFilterMaskTrackAcc[15]	={ 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2.5, 2.5,2*pi, 1.5,15.};
fhSparseFilterMaskTrackAcc=new THnSparseF("fhSparseFilterMaskTrackAcc","fhSparseFilterMaskTrackAcc;fb0;fb1;fb2;fb3;fb4;fb5;fb6;fb7;fb8;fb9;trackIDandsel;SPD;#phi (rad);#eta;p_{T} (GeV/c)",15,nbinsFilterMaskTrackAcc,binlowFilterMaskTrackAcc,binupFilterMaskTrackAcc);
fOutput->Add(fhSparseFilterMaskTrackAcc); 
fhSparseFilterMaskTrackAccPico=new THnSparseF("fhSparseFilterMaskTrackAccPico","fhSparseFilterMaskTrackAccPico;fb0;fb1;fb2;fb3;fb4;fb5;fb6;fb7;fb8;fb9;trackIDandsel;SPD;#phi (rad);#eta;p_{T} (GeV/c)",15,nbinsFilterMaskTrackAcc,binlowFilterMaskTrackAcc,binupFilterMaskTrackAcc);
fOutput->Add(fhSparseFilterMaskTrackAccPico);  

// Num axes: ID*isSelected 5 + kAny,kFirst,kBoth + phi+ eta +pt + imp par
Int_t nbinsFilterMaskImpPar[6]		={   5,   3,  36,  30, 30,  200};
Double_t binlowFilterMaskImpPar[6]	={-2.5,-0.5,  0.,-1.5, 0.,-300.};
Double_t binupFilterMaskImpPar[6]	={ 2.5, 2.5,2*pi, 1.5,15., 300.};
fhSparseFilterMaskImpPar=new THnSparseF("fhSparseFilterMaskImpPar","fhSparseFilterMaskImpPar;fb0;fb1;fb2;fb3;fb4;fb5;fb6;fb7;fb8;fb9;trackIDandsel;SPD;#phi (rad);#eta;p_{T} (GeV/c; imp par (#mum)",6,nbinsFilterMaskImpPar,binlowFilterMaskImpPar,binupFilterMaskImpPar);
fOutput->Add(fhSparseFilterMaskImpPar);  

fhSparseFilterMaskImpParPico=new THnSparseF("fhSparseFilterMaskImpParPico","fhSparseFilterMaskImpParPico;fb0;fb1;fb2;fb3;fb4;fb5;fb6;fb7;fb8;fb9;trackIDandsel;SPD;#phi (rad);#eta;p_{T} (GeV/c; imp par (#mum)",6,nbinsFilterMaskImpPar,binlowFilterMaskImpPar,binupFilterMaskImpPar);
fOutput->Add(fhSparseFilterMaskImpParPico);  

fhImpParResolITSsel=new TH3F("fhImpParResolITSsel","fhImpParResolITSsel;p_{T} (GeV/c);imp. par (#mum);ITS clust",50.,0.,10.,200,-800.,800.,38,-0.5,37.5);
// the convention for ITS clust axis: number between 0 and 48, first digit (units) = number of clust in ITS, second digits (x10): SPD status: 0 = none, 1=kFirst, 2=kSecond; 3= kBoth 
fOutput->Add(fhImpParResolITSsel);

fhImpParResolITSselGoodTracks=new TH3F("fhImpParResolITSselGoodTracks","fhImpParResolITSselGoodTracks;p_{T} (GeV/c);imp. par (#mum);ITS clust",50.,0.,10.,200,-800.,800.,38,-0.5,37.5);
// the convention for ITS clust axis: number between 0 and 48, first digit (units) = number of clust in ITS, second digits (x10): SPD status: 0 = none, 1=kFirst, 2=kSecond; 3= kBoth 
fOutput->Add(fhImpParResolITSselGoodTracks);

// PID PLOTS
fhnSigmaTPCTOFEle	= new TH3F("fhnSigmaTPCTOFEle","fhnSigmaTPCTOFEle;p (GeV/c);n#sigma_{TPC}^{ele};n#sigma_{TOF}^{ele}",120,0.,30,64,-8.,8.,64,-8.,8.);
fhnSigmaTPCTOFPion	= new TH3F("fhnSigmaTPCTOFPion","fhnSigmaTPCTOFPion;p (GeV/c);n#sigma_{TPC}^{ele};n#sigma_{TOF}^{ele}",120,0.,30,64,-8.,8.,64,-8.,8.);
fhnSigmaTPCTOFKaon	= new TH3F("fhnSigmaTPCTOFKaon","fhnSigmaTPCTOFKaon;p (GeV/c);n#sigma_{TPC}^{ele};n#sigma_{TOF}^{ele}",120,0.,30,64,-8.,8.,64,-8.,8.);
fhnSigmaTPCTOFProton= new TH3F("fhnSigmaTPCTOFProton","fhnSigmaTPCTOFProton;p (GeV/c);n#sigma_{TPC}^{ele};n#sigma_{TOF}^{ele}",120,0.,30,64,-8.,8.,64,-8.,8.);

fOutput->Add(fhnSigmaTPCTOFEle);
fOutput->Add(fhnSigmaTPCTOFPion);
fOutput->Add(fhnSigmaTPCTOFKaon);
fOutput->Add(fhnSigmaTPCTOFProton);


//=====================================================================================================================
/// EMCAL PLOTS

//study of NsigmaTPC, e/p, p
Int_t nbinsEoP[4]		= { 202, 400, 100,400};
Double_t binlowEoP[4]	= { -1.,-20.,-20., -1};
Double_t binupEoP[4]	= {100.,  20, 20.,  9};
fhSparseEoverPeleTPC 	= new THnSparseF("fhSparseEoP", "fhSparseEoP; p;nsigmatpc;nsigmaElePIDresp;E/p;",4, nbinsEoP, binlowEoP, binupEoP);
fOutput->Add(fhSparseEoverPeleTPC);
//fSparseEoverPallHadr = new THnSparseF("fSparseEoP", "fSparseEoP; p;nsigmatpc;E/p;",3, nbinsEoP, binlowEoP, binupEoP);

Int_t nbinsEmShSh[7]		= { 35, 120,100,50,50,50,50};
Double_t binlowEmShSh[7]	= { 5.,-20., -1, 0, 0, 0, 0};
Double_t binupEmShSh[7]		= {40.,  10,  9, 1, 1, 2,50};
fhSparseShowShapeEleTPC 	= new THnSparseF("fhSparseShowShapeEleTPC", "fhSparseShowShapeEleTPC; pt;nsigmatpc;E/p;M02;M20;disp;Nclust",7, nbinsEmShSh, binlowEmShSh, binupEmShSh);
fOutput->Add(fhSparseShowShapeEleTPC);  
  
Int_t nbinsTrEM[6]		={124, 20, 124,25,20,20};
Double_t binlowTrEM[6]	={ -1,-1.,-1.0,0.,0.,0.};
Double_t binupTrEM[6]	={ 31, 1., 31.,5.,1.,1.};

fhTrackEMCal=new THnSparseF("fhTrackEMCal","fhTrackEMCal;p(GeV/c);eta;clusterE(GeV/c);E/p;trkDistX;trkDistZ",6,nbinsTrEM,binlowTrEM,binupTrEM);
fOutput->Add(fhTrackEMCal);
//=====================================================================================================================


//=====================================================================================================================
//MCElectron Rejection Histogram
if(fReadMC){}
//=====================================================================================================================

fhptSpectrum=new TH1F("fhptSpectrum", "fhptSpectrum;p_{T};", 100, 0, 50);
fOutput->Add(fhptSpectrum);

 fhptSpectrumPico=new TH1F("fhptSpectrumPico", "fhptSpectrum Pico;p_{T};", 100, 0, 50);
fOutput->Add(fhptSpectrumPico);

PostData(1, fOutput); // Post data for ALL output slots > 0 here.
//cout<<"EXITING CREATING OUTPUT OBJECTS")<<endl;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHFCJQA::FillHistograms()
{
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHFCJQA::CheckClusTrackMatching()
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
    }
    cluster = fCaloClustersCont->GetNextAcceptCluster();
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHFCJQA::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHFCJQA::Run()
{

PrintDebug(0,"Begin");

//========================================================================================
//Getting the Event  
AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
if(!aod)
{
	Printf("ERROR: AOD not available");
	PostData(1,fOutput);
	return kFALSE;
}

AliVEvent *fVevent = dynamic_cast<AliVEvent*> (InputEvent());
if(!fVevent)
{
	Printf("ERROR: VEvent not available");
	PostData(1,fOutput);
	return kFALSE;
}
fhEventCounter->Fill(0);

//TRIGGERS
fhTriggerCounter->Fill(0);
TriggersHistogram(InputEvent()->GetFiredTriggerClasses());

fhTriggerMaskCounter->Fill(0);
TriggersMaskHistogram(InputEvent()->GetTriggerMask());    
TriggersBitHistogram(aod, fReadMC);


//===========================
//========================================================================================

PrintDebug(1,"Event");

//======================================================================
if(!fpidResp) SetupPIDresponse();
//======================================================================

PrintDebug(2,"PIDResponse");

//======================================================================
//Event Selection
bool isSelected = fCuts->IsEventSelected(aod);
if(!isSelected)
{
	AliWarning("Rejecting event ");
	PostData(1,fOutput);
	return kFALSE;
}
fhEventCounter->Fill(1);
//======================================================================

PrintDebug(3,"Event Selection");


//======================================================================
// AOD primary vertex
AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
TString primTitle = vtx1->GetTitle();
if(!(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0)){
  PostData(1,fOutput);
  return kFALSE;  
 }
//======================================================================

PrintDebug(4,"Primary Vertex");

//======================================================================
//Load MC information
TClonesArray *arrayMC=0x0;
AliAODMCHeader *aodmcHeader=0x0;
Double_t vtxTrue[3];

if(fReadMC){
  arrayMC = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!arrayMC)
    {
      Printf("AliAnalysisTaskSEHFjets::UserExec: MC particles branch not found!\n");
      PostData(1,fOutput);
      return kFALSE;
    }
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
 }
//======================================================================

PrintDebug(5,"MC Loading");

//======================================================================
//ESD Vertex
Double_t pos[3],cov[6];
vtx1->GetXYZ(pos);
vtx1->GetCovarianceMatrix(cov);
const AliESDVertex vESD(pos,cov,100.,100);
Double_t magfield=aod->GetMagneticField();
//======================================================================

PrintDebug(6,"ESD Vertex");

// AOD TRACKS, NOT PICO LOOP

for(Int_t itraod=0;itraod<aod->GetNumberOfTracks();itraod++){
  
  AliAODTrack *atrack=(AliAODTrack*)aod->GetTrack(itraod);
  if(!atrack){
    PrintDebug(7,"AOD track loop","null pointer to track");
  }
  //======================================================================
  //Getting the tracks' information
  Double_t p=atrack->P();
  Double_t pt=atrack->Pt();
  Double_t eta=atrack->Eta();
  Double_t pxyz[3];
  atrack->PxPyPz(pxyz);

  //======================================================================

//======================================================================
  //======================================================================
  //PID
  //Double_t nsigma=fpidResp->NumberOfSigmasTPC(atrack, AliPID::kElectron);
  //  PrintDebug(9,"Tracks Loop","PID");

  //======================================================================
  // CHECK FILTER BIT and track ID
  //  AliAODTrack* atrack = dynamic_cast<AliAODTrack*>(track);
  //  if( (!atrack) || (!(atrack->TestFilterBit(ffilterbit))) || (track->GetID()<0) ) continue;
  //  PrintDebug(9,"Tracks Loop","FilterBit");
  //======================================================================

  //======================================================================
  // CHECK FILTER MAPS
  if(!FillTrackHistosAndSelectTrack(atrack,&vESD,magfield,kFALSE))continue;
  fhptSpectrum->Fill(pt);
  //======================================================================
    
  PrintDebug(7,"AOD Tracks Loop","FilterMaps");
  if(!fCheckPico){
    //======================================================================
    // START PID: TPC
    Double_t nsigmaEleTPC=fpidResp->NumberOfSigmasTPC(atrack, AliPID::kElectron);
    Double_t nsigmaPionTPC=fpidResp->NumberOfSigmasTPC(atrack, AliPID::kPion);
    Double_t nsigmaKaonTPC=fpidResp->NumberOfSigmasTPC(atrack, AliPID::kKaon);
    Double_t nsigmaProtonTPC=fpidResp->NumberOfSigmasTPC(atrack, AliPID::kProton);
    
    // TOF
    Double_t nsigmaEleTOF=fpidResp->NumberOfSigmasTOF(atrack, AliPID::kElectron);
    Double_t nsigmaPionTOF=fpidResp->NumberOfSigmasTOF(atrack, AliPID::kPion);
    Double_t nsigmaKaonTOF=fpidResp->NumberOfSigmasTOF(atrack, AliPID::kKaon);
    Double_t nsigmaProtonTOF=fpidResp->NumberOfSigmasTOF(atrack, AliPID::kProton);
    
    
    fhnSigmaTPCTOFEle->Fill(p,nsigmaEleTPC,nsigmaEleTOF);
    fhnSigmaTPCTOFPion->Fill(p,nsigmaPionTPC,nsigmaPionTOF);
    fhnSigmaTPCTOFKaon->Fill(p,nsigmaKaonTPC,nsigmaKaonTOF);
    fhnSigmaTPCTOFProton->Fill(p,nsigmaProtonTPC,nsigmaProtonTOF);
  }
 }


//======================================================================
//Check EMCAL Clusters Energy "trigger"
EnergyTriggers();


//======================================================================
//JetContainer
if(!fJetsCont)
{
	AliWarning("AliAnalysisTaskSEHFjets::Run: Jets Container not found!");
	PostData(1,fOutput);
	return kFALSE;
}
fhEventCounter->Fill(2);
//======================================================================

PrintDebug(8,"Jets Container");

//======================================================================
//TracksContainer
if(!fTracksCont)
  {
    AliWarning("AliAnalysisTaskSEHFjets::Run: Tracks Container not found!");
    PostData(1,fOutput);
    return kFALSE;
  }

fhEventCounter->Fill(3);
//======================================================================

PrintDebug(9,"Tracks Container");


//=========================================================================================
//PICO TRACK LOOP
 if(fCheckPico){
   int NTracks = fTracksCont->GetNAcceptedParticles();
   for(int itr = 0; itr < NTracks; itr++)
     {
       AliPicoTrack *PicoTrack = static_cast<AliPicoTrack*>(fTracksCont->GetAcceptParticle(itr));
       if(!PicoTrack) continue;
       
       AliVTrack* track = PicoTrack->GetTrack();      
       if(!track)
	 {
	   printf("ERROR: Could not receive track %d\n", itr);
	   continue;
	 }
       
       PrintDebug(10,"Tracks Loop");
       //======================================================================
       //Monte Carlo Studies
       bool isMCele = kFALSE;
       if(fReadMC)
	 {
	   int label = track->GetLabel();
	   
	   if( (label>=0)&&(label < arrayMC->GetEntriesFast()) )
	     {
	       AliAODMCParticle* pele = (AliAODMCParticle*)arrayMC->At(label);
	    if(pele)
	      {
		int pdg = pele->GetPdgCode();
		bool isPhysPrim= pele->IsPhysicalPrimary();
		
	      }
	    PrintDebug(9,"Tracks Loop","MC");
	     }
	 }
       //======================================================================
       
       //======================================================================
       //Getting the tracks' information
       Double_t p=track->P();
       Double_t pt=track->Pt();
       Double_t eta=track->Eta();
       Double_t pxyz[3];
       track->PxPyPz(pxyz);
       //======================================================================
       
       
       //======================================================================      
       //======================================================================
       // CHECK FILTER BIT and track ID       
//        if( (!atrack) || (!(atrack->TestFilterBit(ffilterbit))) || (track->GetID()<0) ) continue;
//        PrintDebug(9,"Tracks Loop","FilterBit");
       //======================================================================
       AliAODTrack* atrack = dynamic_cast<AliAODTrack*>(track);
       //======================================================================
       // CHECK FILTER MAPS
       if(!FillTrackHistosAndSelectTrack(atrack,&vESD,magfield,kTRUE))continue;
       //======================================================================

       fhptSpectrumPico->Fill(pt);       
       PrintDebug(9,"Tracks Loop","FilterMaps");
       
       //PID
       //======================================================================
       // START PID: TPC
       Double_t nsigmaEleTPC=fpidResp->NumberOfSigmasTPC(atrack, AliPID::kElectron);
       Double_t nsigmaPionTPC=fpidResp->NumberOfSigmasTPC(atrack, AliPID::kPion);
       Double_t nsigmaKaonTPC=fpidResp->NumberOfSigmasTPC(atrack, AliPID::kKaon);
       Double_t nsigmaProtonTPC=fpidResp->NumberOfSigmasTPC(atrack, AliPID::kProton);
       
       // TOF
       Double_t nsigmaEleTOF=fpidResp->NumberOfSigmasTOF(atrack, AliPID::kElectron);
       Double_t nsigmaPionTOF=fpidResp->NumberOfSigmasTOF(atrack, AliPID::kPion);
       Double_t nsigmaKaonTOF=fpidResp->NumberOfSigmasTOF(atrack, AliPID::kKaon);
       Double_t nsigmaProtonTOF=fpidResp->NumberOfSigmasTOF(atrack, AliPID::kProton);
       
       
       fhnSigmaTPCTOFEle->Fill(p,nsigmaEleTPC,nsigmaEleTOF);
       fhnSigmaTPCTOFPion->Fill(p,nsigmaPionTPC,nsigmaPionTOF);
       fhnSigmaTPCTOFKaon->Fill(p,nsigmaKaonTPC,nsigmaKaonTOF);
       fhnSigmaTPCTOFProton->Fill(p,nsigmaProtonTPC,nsigmaProtonTOF);
       //======================================================================
            
       //======================================================================
       //EMCal Clusters
       Int_t nClsId = track->GetEMCALcluster();
       
       if(nClsId <0) continue;
       
       PrintDebug(9,"Tracks Loop","WithClusters");
       
    //======================================================================
       
    //==============================================================================
    //Cluster Cuts
       AliVCluster *cluster = fCaloClustersCont->GetCluster(nClsId);
       if(!cluster) continue;
       if(!cluster->IsEMCAL()) continue;
       PrintDebug(9,"Tracks Loop","EMCal Cluster");
       //======================================================================
       
       Double_t clsE = cluster->E();
       Double_t nEoverP = clsE/p;
       Double_t eOverPpidResp;
       Double_t showerShape[4];
       Double_t nsigmaEleEMCal=fpidResp->NumberOfSigmasEMCAL(atrack,AliPID::kElectron,eOverPpidResp,showerShape);
       Double_t poix[4]={p, nsigmaEleTPC, nsigmaEleEMCal, nEoverP};
       fhSparseEoverPeleTPC->Fill(poix);
       
       Double_t emcTrackDx=cluster->GetTrackDx();
       Double_t emcTrackDz=cluster->GetTrackDz();
       Double_t pointET[6]={p,eta,clsE,nEoverP,emcTrackDx,emcTrackDz};
       fhTrackEMCal->Fill(pointET);    
       
       Double_t pointEmShSh[7]={atrack->Pt(), static_cast<Double_t>(nsigmaEleTPC), static_cast<Double_t>(nEoverP),cluster->GetM02(),cluster->GetM20(),cluster->GetDispersion(), static_cast<Double_t>(cluster->GetNCells())};
       
       fhSparseShowShapeEleTPC->Fill(pointEmShSh);
    //======================================================================
    
    
       
     }
 }

//=========================================================================================
//JET LOOP
int NJets = fJetsCont->GetNAcceptedJets();
for(int ithJet = 0; ithJet < NJets; ithJet++)
  {
    PrintDebug(10,"Jets Loop");
    AliEmcalJet *emCalJet = fJetsCont->GetAcceptJet(ithJet);
    if(!emCalJet) continue;
    
    Double_t contribution=0,ptpart=-1;
    Int_t partonnat=0;
    
    if(fReadMC)
      {
		AliAODMCParticle *parton=IsMCJet(arrayMC,emCalJet,contribution);
		if(parton)
		{
			Int_t pdg=TMath::Abs(parton->PdgCode());
			//printf("pdg parton: %d \n",pdg);
			if(pdg==21)partonnat=1;
			else if(pdg<4)partonnat=2;
			else if(pdg==4)partonnat=3;
			else if(pdg==5)partonnat=4;
			ptpart=parton->Pt();
		}

    }
    
	//Sets PtSub: pT after background subtraction
	emCalJet->SetPtSub( emCalJet->Pt() - fJetsCont->GetRhoVal() * emCalJet->Area() );

	FillJetRecoHisto(emCalJet,partonnat,contribution,ptpart);
}
//=========================================================================================

PostData(1, fOutput);

PrintDebug(11,"End");
return kFALSE;
//return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHFCJQA::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

//---------------------------------------------------------------
AliAODMCParticle* AliAnalysisTaskEmcalHFCJQA::IsMCJet(TClonesArray *arrayMC,const AliEmcalJet *jet, Double_t &contribution){// assignment of parton ID to jet
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
  if(arrayOfWeights)    delete arrayOfWeights;

  return parton;  
}

//---------------------------------------------------------------
AliAODMCParticle *AliAnalysisTaskEmcalHFCJQA::GetMCPartonOrigin(TClonesArray* &arrayMC,AliAODMCParticle* &p, Int_t &idx)
{  //Follows chain of track mothers until q/g or idx = -1	
  AliAODMCParticle *p2=0x0;
  Int_t mlabel = TMath::Abs(p->GetMother()) ; 
//  Double_t pz=0.;
  while( (mlabel > 1)&&(mlabel < arrayMC->GetEntriesFast()) )
  {
    p2 = (AliAODMCParticle*)arrayMC->At(mlabel);
//    pz=TMath::Abs(p2->Pz());
    //printf("Mother label %d, pdg %d, pz %f\n",mlabel,p2->PdgCode(),pz);
    int pdg = abs(p2->PdgCode());
    if( (pdg == 21) || ( (pdg > 0) && (pdg <6) ) )
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
Int_t AliAnalysisTaskEmcalHFCJQA::TagJetMC(){// METHOD NOT IMPLEMENTED YET
  return -1;
}
//======================================================================

//______________________________________________________________________
void AliAnalysisTaskEmcalHFCJQA::FillJetRecoHisto(const AliEmcalJet *jet,Int_t partonnat,Double_t contribution,Double_t ptpart)
{
//FIll sparse with reco jets properties
  Double_t point[8]={jet->Pt(),jet->Eta(),jet->Phi()-TMath::Pi(), static_cast<Double_t>(jet->GetNumberOfTracks()),0, static_cast<Double_t>(partonnat),contribution,ptpart};
  fSparseRecoJets->Fill(point);
}
//======================================================================

//____________________________________________________________
void AliAnalysisTaskEmcalHFCJQA::SetupPIDresponse(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler(); 
  fpidResp=inputHandler->GetPIDResponse();
  if(!fpidResp)AliFatal("No PID response could be set");
}
//======================================================================

//______________________________________________________________________
Bool_t AliAnalysisTaskEmcalHFCJQA::FillTrackHistosAndSelectTrack(AliAODTrack *aodtr, const AliESDVertex *primary, const Double_t magfield,Bool_t isPico)
{  
  Bool_t retval=kTRUE;
  // THnSparse for filter bits
  Double_t point[16]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,999.,-999.,-999.,-999.,-999.,-999.};
  Double_t pointAcc[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,999.,-999.,-999.,-999.,-999.};
  Double_t pointImpPar[6]={999.,-999.,-999.,-999.,-999.,-999.};

  for(Int_t j=0;j<10;j++){
    if(aodtr->TestFilterBit(TMath::Power(2,j))){
      point[j]=1;      
      pointAcc[j]=1;
    }
  }
  // check ID
  Int_t trid=aodtr->GetID();
  if(aodtr->GetID()>0)point[10]=1.;
  else if(aodtr->GetID()>0)point[10]=0.;
  if(aodtr->GetID()==0)point[10]=-1.;

  Float_t iparxy,iparz;
  Float_t pt=aodtr->Pt();
  Int_t clustITS=aodtr->GetITSNcls();// for histo
  AliESDtrack esdtrack(aodtr);
  // needed to calculate impact parameters
  

  // check refit status
  Int_t refit=-1;
  UInt_t status = esdtrack.GetStatus();
  if(status&AliESDtrack::kTPCrefit)refit+=1;
  if(status&AliESDtrack::kITSrefit)refit+=2;
  point[11]=refit;
  // CHECK SPD
  point[12]=-1;
  if(aodtr->HasPointOnITSLayer(0)){
    clustITS+=10;
    point[12]+=1;
  }
  if(aodtr->HasPointOnITSLayer(1)){
    point[12]+=2;
    clustITS+=20;
  }
  point[13]=aodtr->GetTPCNcls();
  point[14]=aodtr->GetTPCNCrossedRows();
  esdtrack.RelateToVertex(primary,magfield,4.);// CHECK THIS : I put 4.. usually we set it to 3 
  esdtrack.GetImpactParameters(iparxy,iparz);
  point[15]=iparxy;



  
  if(!isPico)fhSparseFilterMask->Fill(point);
  else fhSparseFilterMaskPico->Fill(point);

  if(!aodtr->TestFilterBit(ffilterbit))retval =kFALSE;
  if(retval) PrintDebug(9,"FilterMaps","FilterBit",1);

  if(aodtr->GetID()<0&&!fKeepTrackNegID)retval = kFALSE;
  if(retval){
    if(isPico)fhImpParResolITSsel->Fill(pt,iparxy*10000.,clustITS); 
    else if(!fCheckPico)fhImpParResolITSsel->Fill(pt,iparxy*10000.,clustITS);
  }
  if(retval) PrintDebug(9,"FilterMaps","NegativeTracks",1);
  
  AliESDtrackCuts *cuts=fCuts->GetTrackCuts();
  if(!cuts->IsSelected(&esdtrack)) retval = kFALSE;
  if(retval) PrintDebug(9,"FilterMaps","TrackCuts",1);
  
  if(fCuts->GetUseKinkRejection()){
    AliAODVertex *maybeKink=aodtr->GetProdVertex();
    if(maybeKink->GetType()==AliAODVertex::kKink) retval=kFALSE;
  }
  if(retval) PrintDebug(9,"FilterMaps","KinkRejection",1);

  if(retval){
    pointAcc[10]=1*trid;
  }
  else {
    if(trid!=0)
    pointAcc[10]=2*trid;
    else pointAcc[10]=-3;
  }

  pointAcc[11]=point[12];
  pointAcc[12]=aodtr->Phi();
  if(pointAcc[12]<0.)pointAcc[12]+=2.*TMath::Pi();    
  pointAcc[13]=aodtr->Eta();
  pointAcc[14]=pt;
  if(!isPico)fhSparseFilterMaskTrackAcc->Fill(pointAcc);
  else fhSparseFilterMaskTrackAccPico->Fill(pointAcc);
  pointImpPar[0]=pointAcc[10];
  pointImpPar[1]=pointAcc[11];
  pointImpPar[2]=pointAcc[12];
  pointImpPar[3]=pointAcc[13];
  pointImpPar[4]=pointAcc[14];
  pointImpPar[5]=iparxy;
  if(!isPico)fhSparseFilterMaskImpPar->Fill(pointImpPar);
  else fhSparseFilterMaskImpParPico->Fill(pointImpPar);

  if(retval)  {
    if(isPico)fhImpParResolITSselGoodTracks->Fill(pt,iparxy*10000.,clustITS);  
    else if(!fCheckPico)fhImpParResolITSselGoodTracks->Fill(pt,iparxy*10000.,clustITS);  
  }
  
  return retval;
  
}
//======================================================================

//======================================================================
void AliAnalysisTaskEmcalHFCJQA::PrintDebug(int N, TString Section, TString Sub, int LEVEL)
{
	if(fDebug>LEVEL) cout<<"Debug "<<N<<"-"<<Section<<"-"<<Sub<<endl;
}
//======================================================================

//======================================================================
void AliAnalysisTaskEmcalHFCJQA::TriggersHistogram(TString TriggerClass)
{
  //    printf(TriggerClass.Data());
    if( TriggerClass.Contains( "EMC7" ) && ! TriggerClass.Contains( "EMC7W" ) ) fhTriggerCounter->Fill(1);

    if( TriggerClass.Contains( "EG1" ) ) fhTriggerCounter->Fill(2);
    else if( TriggerClass.Contains( "EG2" ) ) fhTriggerCounter->Fill(3);
    else if( TriggerClass.Contains( "EJ1" ) ) fhTriggerCounter->Fill(4);
    else if( TriggerClass.Contains( "EJ2" ) ) fhTriggerCounter->Fill(5);
    else if( TriggerClass.Contains( "EGA" ) ) fhTriggerCounter->Fill(6);
    else if( TriggerClass.Contains( "EJE" ) ) fhTriggerCounter->Fill(7);
}
//======================================================================

//======================================================================
void AliAnalysisTaskEmcalHFCJQA::TriggersMaskHistogram(int kMask)
{
    //printf(TriggerClass.Data());
    if(kMask == AliVEvent::kEMC7) fhTriggerMaskCounter->Fill(1);
    if(kMask == AliVEvent::kEMCEGA) fhTriggerMaskCounter->Fill(2);
    if(kMask == AliVEvent::kEMCEJE) fhTriggerMaskCounter->Fill(3);
    if( (kMask == AliVEvent::kEMC7)&&(kMask == AliVEvent::kEMCEGA) ) fhTriggerMaskCounter->Fill(4);
    if( (kMask == AliVEvent::kEMC7)&&(kMask == AliVEvent::kEMCEJE) ) fhTriggerMaskCounter->Fill(5);
    if( (kMask == AliVEvent::kEMCEGA)&&(kMask == AliVEvent::kEMCEJE) ) fhTriggerMaskCounter->Fill(6);
}
//======================================================================

//======================================================================
void AliAnalysisTaskEmcalHFCJQA::TriggersBitHistogram(AliAODEvent* aod, bool ReadMC)
{
    fhTriggerBitCounter->Fill(0);
    AliVCaloTrigger& trg = *(aod->GetCaloTrigger("EMCAL"));
    
    bool isL0;
    bool EGA1 = kFALSE;
    bool EGA2 = kFALSE;
    bool EJE1 = kFALSE;
    bool EJE2 = kFALSE;
    
    trg.Reset();
    while (trg.Next())
    {
        Int_t nTimes = 0;
        trg.GetNL0Times(nTimes);
        if (nTimes) isL0 = kTRUE;
        
        Int_t col, row;
        trg.GetPosition(col, row);
        if (col > -1 && row > -1)
        {            
            Int_t bit = 0;
            
            trg.GetTriggerBits(bit);
            if(ReadMC)
            {
                if ((bit >> 1) & 0x1) EGA1 = kTRUE;
                if ((bit >> 2) & 0x1) EGA2 = kTRUE;
                if ((bit >> 3) & 0x1) EJE1 = kTRUE;
                if ((bit >> 4) & 0x1) EJE2 = kTRUE;
            }
            else
            {
                if ((bit >> 6) & 0x1) EGA1 = kTRUE;
                if ((bit >> 7) & 0x1) EGA2 = kTRUE;
                if ((bit >> 8) & 0x1) EJE1 = kTRUE;
                if ((bit >> 9) & 0x1) EJE2 = kTRUE;
            }
            
        }
    }
    
    //Fill the histogram
    if(isL0) fhTriggerBitCounter->Fill(1);// L0
    
    if(EGA1) fhTriggerBitCounter->Fill(2);// EGA1
    if(EGA2) fhTriggerBitCounter->Fill(3);// EGA2
    if(EJE1) fhTriggerBitCounter->Fill(4);// EJE1
    if(EJE2) fhTriggerBitCounter->Fill(5);// EJE2
    
    if(EGA1 && EGA2) fhTriggerBitCounter->Fill(6);
    if(EJE1 && EJE2) fhTriggerBitCounter->Fill(7);
    
    if(EGA1 && EJE1) fhTriggerBitCounter->Fill(8);
    if(EGA1 && EJE2) fhTriggerBitCounter->Fill(9);
    if(EGA2 && EJE1) fhTriggerBitCounter->Fill(10);
    if(EGA2 && EJE2) fhTriggerBitCounter->Fill(11);
    
    //Fill the Energy Distribution
    ClustersEnergyDistribution(isL0,EGA1,EGA2,EJE1,EJE2);
    
    
}
//======================================================================

//======================================================================
//Testing the Clusters Energy Trigger per Event
void AliAnalysisTaskEmcalHFCJQA::EnergyTriggers()
{
    //Energy thresholds to be tested
    double Threshold[6] = {2,4,5,6,8,10};//2, 4,5,6,8,10
    bool kOne = kFALSE;
    bool kTwo = kFALSE;
    bool kThree = kFALSE;
    bool kFour = kFALSE;
    bool kFive = kFALSE;
    bool kSix = kFALSE;

for(int itr = 0; itr < fCaloClustersCont->GetNAcceptedClusters(); itr++)
{
    AliVCluster *cluster = fCaloClustersCont->GetAcceptCluster(itr);
    if(! cluster) continue;
    if(! cluster->IsEMCAL()) continue;
    double E = cluster->E();
    if(E > Threshold[0]) kOne = kTRUE;
    if(E > Threshold[1]) kTwo = kTRUE;
    if(E > Threshold[2]) kThree = kTRUE;
    if(E > Threshold[3]) kFour = kTRUE;
    if(E > Threshold[4]) kFive = kTRUE;
    if(E > Threshold[5]) kSix = kTRUE;
}
fEventsThreshold->Fill(0);
if(kOne)fEventsThreshold->Fill(1);
if(kTwo)fEventsThreshold->Fill(2);
if(kThree)fEventsThreshold->Fill(3);
if(kFour)fEventsThreshold->Fill(4);
if(kFive)fEventsThreshold->Fill(5);
if(kSix)fEventsThreshold->Fill(6);
}
//======================================================================
//======================================================================
//Testing the Clusters Energy Trigger per Event
void AliAnalysisTaskEmcalHFCJQA::ClustersEnergyDistribution(bool isL0, bool EGA1,bool EGA2,bool EJE1,bool EJE2)
{
    for(int itr = 0; itr < fCaloClustersCont->GetNAcceptedClusters(); itr++)
    {
        AliVCluster *cluster = fCaloClustersCont->GetAcceptCluster(itr);
        if(! cluster) continue;
        if(! cluster->IsEMCAL()) continue;
        double E = cluster->E();

        fClustersEnergydistribution->Fill(0.0,E);
        
        if(isL0) fClustersEnergydistribution->Fill(1.,E);
        
        if(EGA1) fClustersEnergydistribution->Fill(2.,E);
        if(EGA2) fClustersEnergydistribution->Fill(3.,E);
        if(EJE1) fClustersEnergydistribution->Fill(4.,E);
        if(EJE2) fClustersEnergydistribution->Fill(5.,E);
        
        if(EGA1 && EGA2) fClustersEnergydistribution->Fill(6.,E);
        if(EJE1 && EJE2) fClustersEnergydistribution->Fill(7.,E);
        
        if(EGA1 && EJE1) fClustersEnergydistribution->Fill(8.,E);
        if(EGA1 && EJE2) fClustersEnergydistribution->Fill(9.,E);
        if(EGA2 && EJE1) fClustersEnergydistribution->Fill(10.,E);
        if(EGA2 && EJE2) fClustersEnergydistribution->Fill(11.,E);
    }

}
//======================================================================
