/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnaysisTaskCaloHFEpp
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskCaloHFEpp.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h" 
#include "AliPID.h" 
#include "AliKFParticle.h"
#include "AliESDtrackCuts.h" 
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"


//using std::cout;
//using std::endl;

class AliAnalysisTaskCaloHFEpp;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskCaloHFEpp) // classimp: necessary for root

AliAnalysisTaskCaloHFEpp::AliAnalysisTaskCaloHFEpp() : AliAnalysisTaskSE(), 
				fAOD(0), 
				fOutputList(0), 
				fVevent(0), 
				fTracks_tender(0),
				fCaloClusters_tender(0),
				fMultSelection(0), 
				fpidResponse(0), 
				fHist_trackPt(0),
				fHistMatchPt(0),
				fHistSelectPt(0),
				fHist_ClustE(0),
				fHist_SelectClustE(0),
				fHistMatchE(0),
				fHistEta_track(0),
				fHistPhi_track(0),
				fHistEta_EMcal(0),
				fHistPhi_EMcal(0),
				fHistScatter_EMcal(0),
				fHist_VertexZ(0),
				fHist_Centrality(0),
				fHist_Mult(0),
				fNevents(0),
				fTrigMulti(0),
				fdEdx(0),
				fTPCnsig(0),
				fHistNsigEop(0),
				fM02(0),
				fM20(0),
				fM02_2(0),
				fM20_2(0),
				fTPCNcls(0),
				fITSNcls(0),
				fEopPt_ele_loose(0),
				fEopPt_ele_tight(0),
				fEopPt_had(0),
				fInvmassLS(0),
				fInvmassULS(0),
				fEtadiff(0),
				fPhidiff(0),
				fHistoNCells(0),
				fEop_electron(0),
				fEop_hadron(0),
				fInv_pT_ULS(0),
				fInv_pT_LS(0),
				//==== Tender flag ====
				fUseTender(kTRUE),
				//==== Trigger or Calorimeter flag ====
				fEMCEG1(kFALSE),
				fEMCEG2(kFALSE),
				fFlagClsTypeEMC(kTRUE),
				fFlagClsTypeDCAL(kFALSE),
				//==== MC output ===
				fMCcheckMother(0),
				fMCarray(0),
				fMCparticle(0),
				fMCheader(0),
				fHistPhoReco0(0),
				fHistPhoReco1(0),
				fHistPhoReco2(0),
				fCheckEtaMC(0),
				fHistMCorgPi0(0),
				fHistMCorgEta(0),
				fHistMCorgD(0),
				fHistMCorgB(0),
				NembMCpi0(0),
				NembMCeta(0),
				NpureMCproc(0),
				fHistPhoPi0(0), 
				fHistPhoPi1(0),
				fHistPhoEta0(0),
				fHistPhoEta1(0),
				fPi010(0),
				fEta010(0),
				fHistPt_HFE_MC_D(0),
				fHistPt_HFE_MC_B(0),
				//fHist_eff_pretrack(0),
				//fHist_eff_posttrack(0),
				fHist_eff_HFE(0),
				fHist_eff_match(0),
				fHist_eff_TPC(0)

{
				// default constructor, don't allocate memory here!
				// this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCaloHFEpp::AliAnalysisTaskCaloHFEpp(const char* name) : AliAnalysisTaskSE(name),
				fAOD(0), 
				fOutputList(0), 
				fVevent(0), 
				fTracks_tender(0),
				fCaloClusters_tender(0),
				fMultSelection(0), 
				fpidResponse(0), 
				fHist_trackPt(0),
				fHistMatchPt(0),
				fHistSelectPt(0),
				fHist_ClustE(0),
				fHist_SelectClustE(0),
				fHistMatchE(0),
				fHistEta_track(0),
				fHistPhi_track(0),
				fHistEta_EMcal(0),
				fHistPhi_EMcal(0),
				fHistScatter_EMcal(0),
				fHist_VertexZ(0),
				fHist_Centrality(0),
				fHist_Mult(0),
				fNevents(0),
				fTrigMulti(0),
				fdEdx(0),
				fTPCnsig(0),
				fHistNsigEop(0),
				fM02(0),
				fM20(0),
				fM02_2(0),
				fM20_2(0),
				fTPCNcls(0),
				fITSNcls(0),
				fEopPt_ele_loose(0),
				fEopPt_ele_tight(0),
				fEopPt_had(0),
				fInvmassLS(0),
				fInvmassULS(0),
				fEtadiff(0),
				fPhidiff(0),
				fHistoNCells(0),
				fEop_electron(0),
				fEop_hadron(0),
				fInv_pT_ULS(0),
				fInv_pT_LS(0),
				//==== Tender flag ====
				fUseTender(kTRUE),
				//==== Trigger or Calorimeter flag ====
				fEMCEG1(kFALSE),
				fEMCEG2(kFALSE),
				fFlagClsTypeEMC(kTRUE),
				fFlagClsTypeDCAL(kFALSE),
				//==== MC output ===
				fMCcheckMother(0),
				fMCarray(0),
				fMCparticle(0),
				fMCheader(0),
				fHistPhoReco0(0),
				fHistPhoReco1(0),
				fHistPhoReco2(0),
				fCheckEtaMC(0),
				fHistMCorgPi0(0),
				fHistMCorgEta(0),
				fHistMCorgD(0),
				fHistMCorgB(0),
				NembMCpi0(0),
				NembMCeta(0),
				NpureMCproc(0),
				fHistPhoPi0(0), 
				fHistPhoPi1(0),
				fHistPhoEta0(0),
				fHistPhoEta1(0),
				fPi010(0),
				fEta010(0),
				fHistPt_HFE_MC_D(0),
				fHistPt_HFE_MC_B(0),
				//fHist_eff_pretrack(0),
				//fHist_eff_postHFE(0),
				fHist_eff_HFE(0),
				fHist_eff_match(0),
				fHist_eff_TPC(0)
{
				// constructor
				DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
				// this chain is created by the analysis manager, so no need to worry about it, 
				// it does its work automatically
				DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
				// you can add more output objects by calling DefineOutput(2, classname::Class())
				// if you add more output objects, make sure to call PostData for all of them, and to
				// make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskCaloHFEpp::~AliAnalysisTaskCaloHFEpp()
{
				// destructor
				if(fOutputList) {
								delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
				}
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::UserCreateOutputObjects()
{
				// create output objects
				//
				// this function is called ONCE at the start of your analysis (RUNTIME)
				// here you ceate the histograms that you want to use 
				//
				// the histograms are in this case added to a tlist, this list is in the end saved
				// to an output file
				//
				fOutputList = new TList();          // this is a list which will contain all of your histograms
				// at the end of the analysis, the contents of this list are written
				// to the output file
				fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
				// if requested (dont worry about this now)

				// example of a histogram
				fHist_trackPt = new TH1F("fHist_trackPt", "EMCAL cluster pt distributiont; pt(GeV/c); counts", 1000, 0, 100);       // create your histogra
				fHistMatchPt = new TH1F("fHistMatchPt", "EMCAL matched cluster pt distributiont; pt(GeV/c); counts", 1000, 0, 100);       // create your histogra
				fHistSelectPt = new TH1F("fHistSelectPt", "EMCAL Slected cluster pt distributiont; pt(GeV/c); counts", 1000, 0, 100);       // create your histogra
				fHist_ClustE = new TH1F("fHist_ClustE", "fHist_ClustE; EMCAL cluster energy distribution; counts", 1000, 0, 100);      
				fHist_SelectClustE = new TH1F("fHist_SelectClustE", "fHistE; EMCAL cluster energy distribution before selection; counts", 1000, 0, 100);      
				fHistMatchE = new TH1F("fHistMatchE", "fHistMatchE; EMCAL Matched cluster energy distribution; counts", 1000, 0, 100);      
				fHistEta_track = new TH1F("fHistEta_track", "Track #eta distribution; pt(GeV/c); counts", 200, -4, 4);    
				fHistPhi_track = new TH1F("fHistPhi_track", "Track #phi distribution; #phi; counts", 200, 0, 10);    
				fHistEta_EMcal = new TH1F("fHistEta_EMcal", "EMCAL selected cluster #eta distribution; #eta; counts", 200, -4, 4);    
				fHistPhi_EMcal = new TH1F("fHistPhi_EMcal", "EMCAL selected cluster #phi distribution; #phi; counts", 200, 0, 10);    
				fHist_VertexZ = new TH1F("fHist_VertexZ", "Z Vertex position; Vtx_{z}; counts", 200, -25, 25);     
				fHist_Centrality = new TH1F("fHist_Centrality", "Centrality", 100, 0, 100);
				fNevents = new TH1F("fNevents","No of events",3,-0.5,2.5);                                                                                                                 
				fTPCNcls = new TH1F("fTPCNcls","No of TPC clusters; N^{TPC}_{cls}; counts",100,0.0,200.);           
				fITSNcls = new TH1F("fITSNcls","No of ITS clusters; N^{ITS}_{cls}; counts",100,0.0,20.); 
				fInvmassLS = new TH1F("fInvmassLS", "Invmass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
				fInvmassULS = new TH1F("fInvmassULS", "Invmass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
				fEtadiff = new TH1F("fEtadiff", "Distance of EMCAL to its closest track(Eta)", 60,-0.3,0.3);
				fPhidiff = new TH1F("fPhidiff", "Distance of EMCAL to its closest track(Phi)", 60,-0.3,0.3);
				fEop_electron = new TH1F("fEop_electron","E/p distribution(-1<n^{TPC}_{#sigma}<3); E/p; counts",60,0,3.0);
				fEop_hadron = new TH1F("fEop_hadron","E/p distribution(n^{TPC}_{#sigma}<-3.5); E/p; counts",60,0,3.0);
				fMCcheckMother = new TH1F("fMCcheckMother", "Mother MC PDG", 1000,-0.5,999.5);
				fHistPhoReco0 = new TH1D("fHistPhoReco0", "total pho in sample; p_{T}(GeV/c)", 60,0,60);
				fHistPhoReco1 = new TH1D("fHistPhoReco1", "reco pho in sample; p_{T}(GeV/c)", 60,0,60);
				fHistPhoReco2 = new TH1D("fHistPhoReco2", "non-reco pho in sample; p_{T}(GeV/c)", 60,0,60);
				fHistPhoPi0 = new TH1D("fHistPhoPi0", "total pi0 in sample; p_{T}(GeV/c)", 60,0,60);
				fHistPhoPi0->Sumw2(); 
				fHistPhoPi1 = new TH1D("fHistPhoPi1", "reco pi0 in sample; p_{T}(GeV/c)", 60,0,60);
				fHistPhoPi1->Sumw2(); 
				fHistPhoEta0 = new TH1D("fHistPhoEta0", "total Eta in sample; p_{T}(GeV/c)", 60,0,60);
				fHistPhoEta0->Sumw2(); 
				fHistPhoEta1 = new TH1D("fHistPhoEta1", "reco Eta in sample; p_{T}(GeV/c)", 60,0,60);
				fHistPhoEta1->Sumw2(); 
				fCheckEtaMC = new TH1F("fCheckEtaMC","check Eta range cut in MC",160,-0.8,0.8);
				fHistMCorgD = new TH1F("fHistMCorgD","MC org D",60,0,60);
				fHistMCorgB = new TH1F("fHistMCorgB","MC org B",60,0,60);
				fHistPt_HFE_MC_D  = new TH1F("fHistPt_HFE_MC_D","HFE from D MC",60,0,60);
				fHistPt_HFE_MC_B  = new TH1F("fHistPt_HFE_MC_B","HFE fron B MC",60,0,60);
				//fHist_eff_pretrack   = new TH1F("fHist_eff_pretrack","efficiency :: before track cut",60,0,60);
				//fHist_eff_posttrack   = new TH1F("fHist_eff_posttrack","efficiency :: afger track cut",60,0,60);
				fHist_eff_HFE     = new TH1F("fHist_eff_HFE","efficiency :: HFE",60,0,60);
				fHist_eff_match   = new TH1F("fHist_eff_match","efficiency :: matched cluster",60,0,60);
				fHist_eff_TPC     = new TH1F("fHist_eff_TPC","efficiency :: TPC cut",60,0,60);





				/////////////////
				//pi0 weight			
				/////////////////
				fPi010 = new TF1("fPi010","[0]*x/pow([1]+x/[2]+x*x/[3],[4])");
			  fPi010->SetParameters(3.68528e+01,5.43694e-02,1.99270e+00,5.33945e+00,3.08814e+00);

				/////////////////
				// Eta weight
				/////////////////
				fEta010 = new TF1("fEta010","[0]*x/pow([1]+x/[2]+x*x/[3],[4])");
				fEta010->SetParameters(1.50102e+01,2.08498e-01,2.95617e+00,5.05032e+00,2.95377e+00);    

				fHist_Mult = new TH2F("fMult","Track multiplicity",100,0,100,20000,0,20000);
				fHistScatter_EMcal = new TH2F("fHistScatter_EMcal", "EMCAL cluster scatter plot; #eta; #phi", 200,0.,6.,200, -1., 1.);       // create your histogra
				fHistScatter_EMcal_aftMatch = new TH2F("fHistScatter_EMcal_aftMatch", "EMCAL cluster scatter plot after track matching; #eta; #phi", 40,-1.0,1.0,200, 0., 6.);       // create your histogra
				fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",500,0,50,500,0,160);
				fTPCnsig = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
				fHistNsigEop = new TH2F ("fHistNsigEop", "E/p vs TPC nsig; E/p; #sigme_{TPC-dE/dX}",60, 0.0, 3.0, 200, -10,10);   
				fM02 = new TH2F ("fM02","M02 vs pt distribution; pt(GeV/c); M02",500,0,50,400,0,2);
				fM20 = new TH2F ("fM20","M20 vs pt distribution; pt(GeV/c); M20",500,0,50,400,0,2);
				fM02_2 = new TH2F ("fM02_2","M02 vs pt distribution (-1<nSigma<3 & 0.9<E/p<1.3); pt(GeV/c); M02",500,0,50,400,0,2);
				fM20_2 = new TH2F ("fM20_2","M20 vs pt distribution (-1<nSigma<3 & 0.9<E/p<1.3); pt(GeV/c); M20",500,0,50,400,0,2);
				fEopPt_ele_loose = new TH2F ("fEopPt_ele_loose","pt vs E/p distribution (-3<nSigma<3); pt(GeV/c); E/p",500,0,50,60,0,3.0);
				fEopPt_ele_tight = new TH2F ("fEopPt_ele_tight","pt vs E/p distribution (-1<nSigma<3); pt(GeV/c); E/p",500,0,50,60,0,3.0);
				fEopPt_had = new TH2F ("fEopPt_had","pt vs E/p distribution (nSigma<-3.5); pt(GeV/c); E/p",500,0,50,60,0,3.0);
				fHistoNCells = new TH2F("fHistoNCells", "No of EMCAL cells in a cluster; Cluster E; N^{EMC}_{cells}",500,0,50,30,0,30);
				fInv_pT_ULS = new TH2F("fInv_pT_ULS", "Invariant mass vs p_{T} distribution(ULS) ; pt(GeV/c) ; mass(GeV/c^2)",500,0,50,1000,0,1.0);
				fInv_pT_LS = new TH2F("fInv_pT_LS", "Invariant mass vs p_{T} distribution(LS) ; pt(GeV/c) ; mass(GeV/c^2)",500,0,50,1000,0,1.0);
				fHistMCorgPi0 = new TH2F("fHistMCorgPi0","MC org Pi0",2,-0.5,1.5,100,0,50);
				fHistMCorgEta = new TH2F("fHistMCorgEta","MC org Eta",2,-0.5,1.5,100,0,50);
				fTrigMulti = new TH2F("fTrigMulti","Multiplicity distribution for different triggers; Trigger type; multiplicity",11,-1,10,2000,0,2000);

				fOutputList->Add(fHist_trackPt);          // don't forget to add it to the list! the list will be written to file, so if you want
				fOutputList->Add(fHistMatchPt);          
				fOutputList->Add(fHistSelectPt);          
				fOutputList->Add(fHist_ClustE);          
				fOutputList->Add(fHist_SelectClustE);          
				fOutputList->Add(fHistMatchE);          
				fOutputList->Add(fHistEta_track);         
				fOutputList->Add(fHistPhi_track);         
				fOutputList->Add(fHistEta_EMcal);         
				fOutputList->Add(fHistPhi_EMcal);         
				fOutputList->Add(fHistScatter_EMcal);     
				fOutputList->Add(fHistScatter_EMcal_aftMatch);     
				fOutputList->Add(fHist_VertexZ);          
				fOutputList->Add(fHist_Centrality);       
				fOutputList->Add(fHist_Mult);           
				fOutputList->Add(fTrigMulti);
				fOutputList->Add(fNevents);
				fOutputList->Add(fdEdx);
				fOutputList->Add(fTPCnsig);
				fOutputList->Add(fHistNsigEop);
				fOutputList->Add(fM02);
				fOutputList->Add(fM20);
				fOutputList->Add(fM02_2);
				fOutputList->Add(fM20_2);
				fOutputList->Add(fTPCNcls);
				fOutputList->Add(fITSNcls);
				fOutputList->Add(fEopPt_ele_loose);
				fOutputList->Add(fEopPt_ele_tight);
				fOutputList->Add(fEopPt_had);
				fOutputList->Add(fInvmassLS);
				fOutputList->Add(fInvmassULS);
				fOutputList->Add(fEtadiff);
				fOutputList->Add(fPhidiff);
				fOutputList->Add(fHistoNCells);
				fOutputList->Add(fEop_electron);
				fOutputList->Add(fEop_hadron);
				fOutputList->Add(fInv_pT_ULS);
				fOutputList->Add(fInv_pT_LS);
				//==== MC output ====
				fOutputList->Add(fMCcheckMother);
				fOutputList->Add(fHistPhoReco0);
				fOutputList->Add(fHistPhoReco1);
				fOutputList->Add(fHistPhoReco2);
				fOutputList->Add(fCheckEtaMC);
				fOutputList->Add(fHistMCorgPi0);
				fOutputList->Add(fHistMCorgEta);
				fOutputList->Add(fHistMCorgD);
				fOutputList->Add(fHistMCorgB);
				fOutputList->Add(fHistPhoPi0);
				fOutputList->Add(fHistPhoPi1);
				fOutputList->Add(fHistPhoEta0);
				fOutputList->Add(fHistPhoEta1);
				fOutputList->Add(fHistPt_HFE_MC_D);
				fOutputList->Add(fHistPt_HFE_MC_B);
				//fOutputList->Add(fHist_eff_pretrack);
				//fOutputList->Add(fHist_eff_posttrack);
				fOutputList->Add(fHist_eff_HFE); 
				fOutputList->Add(fHist_eff_match); 
				fOutputList->Add(fHist_eff_TPC); 




				PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
				// fOutputList object. the manager will in the end take care of writing your output to file
				// so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::UserExec(Option_t *)
{

				UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
				////////////////////
				//cuts initialised//
				////////////////////
				AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
				esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
				esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
				esdTrackCutsH->SetDCAToVertex2D(kTRUE);
				esdTrackCutsH->SetMinNClustersTPC(80);
				esdTrackCutsH->SetMinNClustersITS(3);
				esdTrackCutsH->SetRequireTPCRefit(kTRUE);
				esdTrackCutsH->SetRequireITSRefit(kTRUE);
				esdTrackCutsH->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
				esdTrackCutsH->SetAcceptKinkDaughters(kFALSE);
				esdTrackCutsH->SetMaxChi2PerClusterITS(6); //test.....


				// user exec
				// this function is called once for each event
				// the manager will take care of reading the events from file, and with the static function InputEvent() you 
				// have access to the current event. 
				// once you return from the UserExec function, the manager will retrieve the next event from the chain
				fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
				// there's another event format (ESD) which works in a similar wya
				// but is more cpu/memory unfriendly. for now, we'll stick with aod's
				if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
				// example part: i'll show how to loop over the tracks in an event 
				// and extract some information from them which we'll store in a histogram

				fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
				fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));

				fVevent = dynamic_cast<AliVEvent*>(InputEvent());
				if (!fVevent) {
								printf("ERROR: fVEvent not available\n");
								return;
				}   

       
				//////////////////////////////
				//Get Tender?  
				//////////////////////////////
				Bool_t fFlagEMCalCorrection = kTRUE;
				if(fFlagEMCalCorrection){
				TString fTenderClusterName("caloClusters"); //default name
				TString fTenderTrackName("tracks"); //default name
				fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderTrackName));
				fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderClusterName));
				}

				//////////////////////////////
				//PID initialised   
				//////////////////////////////
				fpidResponse = fInputHandler->GetPIDResponse();

				//////////////////////////////
				//Vertex 
				//////////////////////////////
				fNevents->Fill(0);//all enent
				const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
				Double_t NcontV = pVtx->GetNContributors();
				if(NcontV<2)return;
				fNevents->Fill(1); //events with 2 tracks

				Double_t Xvertex = pVtx->GetX();
				Double_t Yvertex = pVtx->GetY();
				Double_t Zvertex = pVtx->GetZ();
				fHist_VertexZ->Fill(Zvertex);                     // plot the pt value of the track in a histogram



				/////////////////
				//trigger check//
				/////////////////
				TString firedTrigger;
				TString TriggerEG1("EG1");
				TString TriggerEG2("EG2");
				fVevent->GetFiredTriggerClasses();
				if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();

				Bool_t EG1tr = kFALSE;
				Bool_t EG2tr = kFALSE;

				if(firedTrigger.Contains(TriggerEG1))EG1tr = kTRUE;
				if(firedTrigger.Contains(TriggerEG2))EG2tr = kTRUE;

				if(fEMCEG1){if(!firedTrigger.Contains(TriggerEG1))return;}
				if(fEMCEG2){if(!firedTrigger.Contains(TriggerEG2))return;}

				Int_t trigger = -1;
				if (fAOD){
								AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
								if(!header) AliFatal("Not a standard AOD");
								Double_t multiplicity = header->GetRefMultiplicity();

								fTrigMulti->Fill(-0.5, multiplicity);
								if(evSelMask & AliVEvent::kAny) fTrigMulti->Fill(0.5, multiplicity);
								if(evSelMask & AliVEvent::kMB) fTrigMulti->Fill(1.5, multiplicity);
								if(evSelMask & AliVEvent::kINT7) fTrigMulti->Fill(2.5, multiplicity);
								if(evSelMask & AliVEvent::kINT8) fTrigMulti->Fill(3.5, multiplicity);
								if(evSelMask & AliVEvent::kEMC1) fTrigMulti->Fill(4.5, multiplicity);
								if(evSelMask & AliVEvent::kEMC7) fTrigMulti->Fill(5.5, multiplicity);
								if(evSelMask & AliVEvent::kEMC8) fTrigMulti->Fill(6.5, multiplicity);
								if(evSelMask & AliVEvent::kEMCEJE) fTrigMulti->Fill(7.5, multiplicity);
								if(evSelMask & AliVEvent::kEMCEGA) fTrigMulti->Fill(8.5, multiplicity);
								if(evSelMask & AliVEvent::kEMCEGA & EG2tr) fTrigMulti->Fill(9.5, multiplicity);

								if(evSelMask & AliVEvent::kMB) trigger =0;
								if(evSelMask & AliVEvent::kINT7) trigger =1;
								if(evSelMask & AliVEvent::kINT8) trigger =2;
								if(evSelMask & AliVEvent::kEMC1) trigger =3;
								if(evSelMask & AliVEvent::kEMC7) trigger =4;
								if(evSelMask & AliVEvent::kEMC8) trigger =5;
								if(evSelMask & AliVEvent::kEMCEJE) trigger =6;
								if(evSelMask & AliVEvent::kEMCEGA) trigger =7;
				}


				//////////////////////////////
				//Centarality
				//////////////////////////////
				Double_t centrality = -1;
				AliCentrality *fCentrality = (AliCentrality*)fAOD->GetCentrality();
				//centrality = fCentrality->GetCentralityPercentile("V0M");
				//cout << "Centrality == " << fCentrality->GetCentralityPercentile("V0M") << endl;
				if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection"); 
				if(!fMultSelection) {
								//If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
								//AliWarning("AliMultSelection object not found!");
								centrality = fCentrality->GetCentralityPercentile("V0M");
				}else{
								//lPercentile = fMultSelection->GetMultiplicityPercentile("V0M");
								centrality = fMultSelection->GetMultiplicityPercentile("V0M"); 
				}


				//////////////////////////////
				// Event selection
				//////////////////////////////
				if(TMath::Abs(Zvertex)>10.0)return;
				fNevents->Fill(2); //Zvertex < 10
				fHist_Centrality -> Fill(centrality);

				if(fMCarray)CheckMCgen(fMCheader);

				//////////////////////////////
				// EMCal cluster loop
				//////////////////////////////
				Int_t Nclust = -999;
				if(!fFlagEMCalCorrection)Nclust =  fVevent->GetNumberOfCaloClusters();
				if(fFlagEMCalCorrection)Nclust =  fCaloClusters_tender->GetEntries();

				Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;;

				for(Int_t icl=0; icl<Nclust; icl++)
				{
								AliVCluster *clust = 0x0;     
								if(!fFlagEMCalCorrection)clust = (AliVCluster*)fVevent->GetCaloCluster(icl); // address cluster matched to track
								if(fFlagEMCalCorrection)clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(icl)); // address cluster matched to track

								fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;

								if(clust && clust->IsEMCAL())
								{
												Double_t clustE = clust->E();
												fHist_SelectClustE -> Fill(clustE);

												Float_t clustpos[3] = {0.};
												clust->GetPosition(clustpos);

												TVector3 pos(clustpos);

												Double_t Phi =  pos.Phi();
												if(Phi <0){Phi += 2*TMath::Pi();}
												//fHistEta_EMcal->Fill(pos.Eta());                     // plot the pt value of the track in a histogram
												//fHistPhi_EMcal->Fill(pos.Phi());                     // plot the pt value of the track in a histogram

												if(Phi > 1.39 && Phi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187     
												if(Phi > 4.53 && Phi < 5.708) fClsTypeDCAL = kTRUE; //DCAL  : 260 < phi < 327

												//----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
												if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
																if(!fClsTypeEMC) continue; //selecting only EMCAL clusters

												if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
																if(!fClsTypeDCAL) continue; //selecting only DCAL clusters


												fHistScatter_EMcal->Fill(Phi,pos.Eta());                     // plot the pt value of the track in a histogram
												fHistoNCells -> Fill(clustE, clust->GetNCells());
												fHist_ClustE->Fill(clustE);                     // plot the pt value of the track in a histogram
								}
				}

				//////////////////////////////
				// Track loop
				//////////////////////////////
				//Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
				//Int_t iTracks(fVevent->GetNumberOfTracks());           // see how many tracks there are in the event
				Int_t iTracks = -999;
				if(!fFlagEMCalCorrection)iTracks = fVevent->GetNumberOfTracks();           // see how many tracks there are in the event
				if(fFlagEMCalCorrection)iTracks = fTracks_tender->GetEntries();           // see how many tracks there are in the event
				//fHist_Mult->Fill(centrality,nTracks);

				for(Int_t i(0); i < iTracks; i++) {                 // loop ove rall these tracks
								Double_t fTPCnSigma = -999, dEdx = -999, TrkP = -999, TrkPt = -999; 

								AliAODTrack* track;         // get a track (type AliAODTrack) from the event
								if(fFlagEMCalCorrection){
												AliVParticle* Vtrack = 0x0;
												Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(i));
												track = dynamic_cast<AliAODTrack*>(Vtrack);
								}
								if(!fFlagEMCalCorrection)track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));   // get a track (type AliAODTrack) from the event

								if(!track) continue;                            // if we failed, skip this track
								fHist_trackPt->Fill(track->Pt());                     // plot the pt value of the track in a histogram
								fHistEta_track->Fill(track->Eta());                     // plot the pt value of the track in a histogram
								fHistPhi_track->Fill(track->Phi());                     // plot the pt value of the track in a histogram
								dEdx = track->GetTPCsignal();
								TrkP = track->P();
								TrkPt = track->Pt();
								fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron); 
								fTPCNcls->Fill(track->GetTPCNcls());
								fITSNcls->Fill(track->GetITSNcls());
								fdEdx->Fill(TrkP,dEdx);
								fTPCnsig->Fill(TrkP,fTPCnSigma);
								//printf( "TPCnSigma ::  %f \n" ,fTPCnSigma); 

								Int_t EMCalIndex = -1;
								EMCalIndex = track->GetEMCALcluster();  // get index of EMCal cluster which matched to track

								//fHist_eff_pretrack->Fill(TrkPt);

								/////////////////////////
								// track cut
								/////////////////////////
								if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
								if(track->GetTPCNcls() < 80) continue; //TPC cluster cut
								if(track->GetITSNcls() < 3) continue;  //ITS cluster cut
								if(!(track -> HasPointOnITSLayer(0) || track -> HasPointOnITSLayer(1))) continue;
								if(track->Eta()>0.6 || track->Eta()<-0.6) continue; //Eta cut

								Double_t DCA[2] = {-999.,-999.}, covar[3]; //DCA cut
								if(track -> PropagateToDCA(pVtx,fVevent -> GetMagneticField(),20.,DCA,covar))
								{
												if(TMath::Abs(DCA[0]) > 2.4 || TMath::Abs(DCA[1]) > 3.2) continue;
								}

								//fHist_eff_posttrack->Fill(TrkPt);

								///////////////////////
								// Get MC information//
								///////////////////////
								Int_t ilabel = TMath::Abs(track->GetLabel());
								Int_t pdg = -999;
								Double_t pid_ele = 0.0;
								Double_t pTmom = -1.0;
								Int_t pidM = -1;
								Int_t ilabelM = -1;
								Bool_t iEmbPi0 = kFALSE; 
								Bool_t iEmbEta = kFALSE;
								Bool_t pid_eleD = kFALSE;
								Bool_t pid_eleB = kFALSE;
								Bool_t pid_eleP = kFALSE;

								if(ilabel>0 && fMCarray)
								{
												fMCparticle = (AliAODMCParticle*) fMCarray->At(ilabel);
												pdg = fMCparticle->GetPdgCode();
												if(TMath::Abs(pdg)==11)pid_ele = 1.0;
												if(pid_ele==1.0)FindMother(fMCparticle, ilabelM, pidM, pTmom);

												pid_eleD = IsDdecay(pidM);
												pid_eleB = IsBdecay(pidM);
												pid_eleP = IsPdecay(pidM);

												if(pidM==111)
												{
																if(ilabelM>=NembMCpi0 && ilabelM<NembMCeta)iEmbPi0 = kTRUE;
																if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;
												}
												if(pidM==221)
												{
																if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;
												}

												if(pidM==22) // from pi0 & eta
												{
																AliAODMCParticle* fMCparticleM = (AliAODMCParticle*) fMCarray->At(ilabelM);
																FindMother(fMCparticleM, ilabelM, pidM, pTmom);

																if(pidM==111)
																{
																				if(ilabelM>=NembMCpi0 && ilabelM<NembMCeta)iEmbPi0 = kTRUE;
																				if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;
																}
																if(pidM==221)
																{
																				if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;
																}
												}
												fMCcheckMother->Fill(abs(pidM));
								}

								if(pidM==443)continue; // remove enhanced J/psi in MC !
								if(pidM==-99)continue; // remove e from no mother !

								if(pid_eleB || pid_eleD) {
												fHist_eff_HFE->Fill(TrkPt);
												if(fTPCnSigma<3 && fTPCnSigma>-1) fHist_eff_TPC->Fill(TrkPt);
								}				


								//if(pid_ele==1.0)cout << "pidM = " << pidM << " ; " << pid_eleP << endl;  
								//if(pid_eleD)fHistDCAde->Fill(track->Pt(),DCAxy);
								//if(pid_eleB)fHistDCAbe->Fill(track->Pt(),DCAxy);
								//if(pid_eleP)fHistDCApe->Fill(track->Pt(),DCAxy);


								//////////////////////////////////////
								//calculate weight of photon for MC  
								//////////////////////////////////////
								Double_t WeightPho = -1.0;

								if(iEmbPi0)
								{
												WeightPho = fPi010->Eval(pTmom);
								}
								if(iEmbEta)
								{
												WeightPho = fEta010->Eval(pTmom);
								}


								AliVCluster *clustMatch=0x0;
								//cout << "EMCalIndex = " << EMCalIndex << endl;
								//if(EMCalIndex>=0)continue;
								if(EMCalIndex<0)continue;
								if(!fFlagEMCalCorrection)clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex); // address cluster matched to track
								if(fFlagEMCalCorrection) clustMatch = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex));
								fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
								if(clustMatch && clustMatch->IsEMCAL())
								{
												fHistMatchPt->Fill(TrkPt);

												///////get position of clustMatch/////////
												Float_t clustMatchpos[3] = {0.};
												clustMatch->GetPosition(clustMatchpos);
												TVector3 cpos(clustMatchpos);
												Double_t Matcheta = cpos.Eta();
												Double_t Matchphi = cpos.Phi();

												///////calculate phi and eta difference between a track and a cluster//////////
												Double_t phidiff = -999;
												Double_t etadiff = -999;
												etadiff = track->GetTrackEtaOnEMCal()-Matcheta;
												phidiff = TVector2::Phi_mpi_pi(track->GetTrackPhiOnEMCal()-Matchphi);
												fEtadiff->Fill(etadiff); 
												fPhidiff->Fill(phidiff); 
												Matchphi=TVector2::Phi_mpi_pi(Matchphi);

												if(TMath::Abs(etadiff)>0.05 || TMath::Abs(phidiff)>0.05) continue;
												if(Matchphi>1.39 && Matchphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187     
												if(Matchphi>4.53 && Matchphi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327

												//----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
												if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
																if(!fClsTypeEMC) continue; //selecting only EMCAL clusters

												if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
																if(!fClsTypeDCAL) continue; //selecting only DCAL clusters

												fHistSelectPt->Fill(TrkPt);
												if(TrkPt>1.0){
												fHistEta_EMcal->Fill(track->Eta()); 
												fHistPhi_EMcal->Fill(track->Phi());
												}


												if(fTPCnSigma<3 && fTPCnSigma>-1){
																if(pid_eleB || pid_eleD) fHist_eff_match->Fill(TrkPt);
												}

												fHistScatter_EMcal_aftMatch->Fill(track->Eta(),track->Phi());	

												Double_t eop = -1.0;
												Double_t clE = clustMatch->E();
												Double_t m20 = clustMatch->GetM20();
												Double_t m02 = clustMatch->GetM02();
												fHistMatchE -> Fill(clE);
												if(TrkP>0)eop= clE/TrkP;
												if(TrkPt>3){
																fHistNsigEop -> Fill(eop,fTPCnSigma);
																if(fTPCnSigma<3 && fTPCnSigma>-1){
																				fEop_electron -> Fill(eop);
																}
																if(fTPCnSigma<-3.5){
																				fEop_hadron -> Fill(eop);
																}
												}

												fM02->Fill(TrkPt,m02);
												fM20->Fill(TrkPt,m20);

												Bool_t fFlagNonHFE=kFALSE; 
												//if(fTPCnSigma<6 && fTPCnSigma>-6 && eop < 1.2&& eop > 0.8 && m20>0.02 && m20<0.25){ //for MC
												if(fTPCnSigma<3 && fTPCnSigma>-1 && eop < 1.2&& eop > 0.8 && m20>0.02 && m20<0.25){
																fM02_2->Fill(TrkPt,m02);
																fM20_2->Fill(TrkPt,m20);

																///////-----Identify Non-HFE////////////////////////////
																SelectPhotonicElectron(iTracks,track,fFlagNonHFE,pidM,TrkPt);
																if(pid_eleP)
																{
																				fHistPhoReco0->Fill(track->Pt()); // reco pho
																				if(iEmbPi0)fHistPhoPi0->Fill(track->Pt(),WeightPho); // reco pho
																				if(iEmbEta)fHistPhoEta0->Fill(track->Pt(),WeightPho); // reco pho

																				if(fFlagNonHFE)
																				{
																								fHistPhoReco1->Fill(track->Pt()); // reco pho
																								if(iEmbPi0)fHistPhoPi1->Fill(track->Pt(),WeightPho); // reco pho
																								if(iEmbEta)fHistPhoEta1->Fill(track->Pt(),WeightPho); // reco pho
																				}
																				else
																				{
																								fHistPhoReco2->Fill(track->Pt()); // org pho
																				}
																}
												}



												if(fTPCnSigma<3 && fTPCnSigma>-3 && m20>0.02 && m20<0.25){
																fEopPt_ele_loose -> Fill(TrkPt,eop);
												}
												if(fTPCnSigma<3 && fTPCnSigma>-1 && m20>0.02 && m20<0.25){
																fEopPt_ele_tight -> Fill(TrkPt,eop);
																if(eop < 1.2&& eop > 0.8){
																				if(pid_eleB) fHistPt_HFE_MC_B -> Fill(track->Pt());
																				if(pid_eleD) fHistPt_HFE_MC_D -> Fill(track->Pt());
																}
																//if(ilabelM<NpureMC){fHistPt_HFE_PYTHIA -> Fill(track->Pt());}
																//else {fHistPt_HFE_emb -> Fill(track->Pt());}
												}
												if(fTPCnSigma<-3.5 && m20>0.02 && m20<0.25){
																fEopPt_had -> Fill(TrkPt,eop);
												}
								}

								}                                                   // continue until all the tracks are processed
				PostData(1, fOutputList);                           // stream the results the analysis of this event to
				// the output manager which will take care of writing
				// it to a file
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::Terminate(Option_t *)
{
				// terminate
				// called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC, Double_t TrkPt)
{
				////// ////////////////////////////////////
				//////Non-HFE - Invariant mass method//////
				///////////////////////////////////////////

				AliESDtrackCuts* esdTrackCutsAsso = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
				esdTrackCutsAsso->SetAcceptKinkDaughters(kFALSE);
				esdTrackCutsAsso->SetRequireTPCRefit(kTRUE);
				esdTrackCutsAsso->SetRequireITSRefit(kTRUE);
				esdTrackCutsAsso->SetEtaRange(-0.9,0.9);
				esdTrackCutsAsso->SetMaxChi2PerClusterTPC(4);
				esdTrackCutsAsso->SetMinNClustersTPC(70);
				esdTrackCutsAsso->SetMaxDCAToVertexZ(3.2);
				esdTrackCutsAsso->SetMaxDCAToVertexXY(2.4);
				esdTrackCutsAsso->SetDCAToVertex2D(kTRUE);

				Bool_t flagPhotonicElec = kFALSE;
				Bool_t fFlagEMCalCorrection = kTRUE;


				Int_t ntracks = -999;
				if(!fFlagEMCalCorrection)ntracks = fVevent->GetNumberOfTracks();
				if(fFlagEMCalCorrection) ntracks = fTracks_tender->GetEntries();

				for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
								AliVParticle* VAssotrack = 0x0;
								if(!fFlagEMCalCorrection) VAssotrack  = fVevent->GetTrack(jtrack);
								if(fFlagEMCalCorrection) VAssotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list

								if (!VAssotrack) {
												printf("ERROR: Could not receive track %d\n", jtrack);
												continue;
								}

								AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
								AliESDtrack *eAssotrack = dynamic_cast<AliESDtrack*>(VAssotrack);
								AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);

								//------reject same track
								if(jtrack==itrack) continue;
								if(aAssotrack->Px()==track->Px() && aAssotrack->Py()==track->Py() && aAssotrack->Pz()==track->Pz())continue;

								Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
								Double_t ptAsso=-999., nsigma=-999.0, mass=-999., width = -999;
								Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;

								nsigma = fpidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
								ptAsso = Assotrack->Pt();
								Int_t chargeAsso = Assotrack->Charge();
								Int_t charge = track->Charge();
								if(charge>0) fPDGe1 = -11;
								if(chargeAsso>0) fPDGe2 = -11;
								if(charge == chargeAsso) fFlagLS = kTRUE;
								if(charge != chargeAsso) fFlagULS = kTRUE;


								//------track cuts applied
								if(fAOD) {
												if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
												//if(aAssotrack->GetTPCNcls() < 70) continue;
												if(aAssotrack->GetTPCNcls() < 80) continue;
												if(aAssotrack->GetITSNcls() < 3 ) continue;
												if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
								}
								else{
												if(!esdTrackCutsAsso->AcceptTrack(eAssotrack)) continue;
								}

								//-------loose cut on partner electron
								if(ptAsso <0.2) continue;
								if(aAssotrack->Eta()<-0.9 || aAssotrack->Eta()>0.9) continue;
								if(nsigma < -3 || nsigma > 3) continue;

								//-------define KFParticle to get mass
								AliKFParticle::SetField(fVevent->GetMagneticField());
								AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
								AliKFParticle ge2 = AliKFParticle(*Assotrack, fPDGe2);
								AliKFParticle recg(ge1, ge2);

								if(recg.GetNDF()<1) continue;
								Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
								if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

								//-------Get mass
								Int_t MassCorrect;
								MassCorrect = recg.GetMass(mass,width);

								if(fFlagLS){
												if(mass < 0.002)cout <<"Px="<<aAssotrack->Px() <<" Py="<<aAssotrack->Py()<<" Pz="<<aAssotrack->Pz()<<endl;
												if(mass < 0.002)cout <<"Px="<<track->Px() <<" Py="<<track->Py()<<" Pz="<<track->Pz()<<endl;
												if(track->Pt()>1) {fInvmassLS->Fill(mass);
																           fInv_pT_LS->Fill(TrkPt,mass);}
								}
								if(fFlagULS){
												if(track->Pt()>1) {fInvmassULS->Fill(mass);
																           fInv_pT_ULS->Fill(TrkPt,mass);}
								}

								if(iMC>0)
								{
												Int_t iMCbin = -999;
												if(iMC == 111) 
												{
																iMCbin = 1;
												}
												else if(iMC == 221)
												{
																iMCbin = 2;
												}
												else
												{
																iMCbin = -999;
												}  

												//if(fFlagULS && track->Pt()>1.5 && iMCbin!=-999)fInvmassULS_MCtrue->Fill(iMCbin,mass);  
								}

								if(mass<0.1 && fFlagULS && !flagPhotonicElec)
												flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised) 
				}
				fFlagPhotonicElec = flagPhotonicElec;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCaloHFEpp::IsDdecay(int mpid)
{
 int abmpid = TMath::Abs(mpid);
 if(abmpid==411 || abmpid==421 || abmpid==413 || abmpid==423 || abmpid==431 || abmpid==433)
   {
    return kTRUE;
   }
 else
   {
    return kFALSE;
   } 
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskCaloHFEpp::IsBdecay(int mpid)
{
 int abmpid = TMath::Abs(mpid);
 if(abmpid==511 || abmpid==521 || abmpid==513 || abmpid==523 || abmpid==531 || abmpid==533)
   {
    return kTRUE;
   }
 else
   {
    return kFALSE;
   } 
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskCaloHFEpp::IsPdecay(int mpid)
{
 int abmpid = TMath::Abs(mpid);
 if(abmpid==22 || abmpid==111 || abmpid==221)
   {
		//fMCcheckMother->Fill(abs(abmpid));
    return kTRUE;
   }
 else
   {
    return kFALSE;
   } 
}


//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom)
{

 if(part->GetMother()>-1)
   {
    label = part->GetMother();
    AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(label);
    pid = partM->GetPdgCode();
		ptmom = partM->Pt();
   }
 else
   {
    pid = -99;
   } 
   //cout << "Find Mother : label = " << label << " ; pid" << pid << endl;
}

//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::CheckMCgen(AliAODMCHeader* fMCheader)
{
 TList *lh=fMCheader->GetCocktailHeaders();
 Int_t NpureMC = 0;
 NpureMCproc = 0;
 NembMCpi0 = 0;
 NembMCeta = 0;
 TString MCgen;
 TString embpi0("pi");
 TString embeta("eta");

 if(lh)
    {     
     //cout << "<------- lh = " << lh << " ; NproAll = "<<  lh->GetEntries() << endl; 

     for(int igene=0; igene<lh->GetEntries(); igene++)
        {
         AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
         if(gh)
           {
            cout << "<------- imc = "<< gh->GetName() << endl;     
            MCgen =  gh->GetName();     
            cout << "<-------- Ncont = " << gh->NProduced() << endl;
            if(igene==0)NpureMC = gh->NProduced();  // generate by PYTHIA or HIJING
           
            //if(MCgen.Contains(embpi0))cout << MCgen << endl;
            //if(MCgen.Contains(embeta))cout << MCgen << endl;
            if(MCgen.Contains(embpi0))NembMCpi0 = NpureMCproc;
            if(MCgen.Contains(embeta))NembMCeta = NpureMCproc;

            NpureMCproc += gh->NProduced();  // generate by PYTHIA or HIJING
           }
        }
    }

 //cout << "NpureMC =" << NpureMC << endl;
 //cout << "NembMCpi0 =" << NembMCpi0 << endl;
 //cout << "NembMCeta =" << NembMCeta << endl;
 //cout << "NpureMCproc =" << NpureMCproc << endl;

 //for(int imc=0; imc<fMCarray->GetEntries(); imc++)
 for(int imc=0; imc<NpureMCproc; imc++)
     {
      //cout << "imc = " << imc << endl;
      Bool_t iEnhance = kFALSE;
      if(imc>=NpureMC)iEnhance = kTRUE;
      Int_t iHijing = 1;  // select particles from Hijing or PYTHIA

      //if(imc==NpureMC)cout << "========================" << endl;  

      fMCparticle = (AliAODMCParticle*) fMCarray->At(imc);
      Int_t pdgGen = TMath::Abs(fMCparticle->GetPdgCode());
      Double_t pdgEta = fMCparticle->Eta(); 
      Double_t pTtrue = fMCparticle->Pt(); 
			if(imc>=NembMCeta && imc<NpureMCproc){
			cout<<"!!!!!!!!!!!!!!PGD = "<<pdgGen<<" imc = "<<imc<<endl;
			}

			//cout << "fetarange = " << fetarange << endl;
			if(TMath::Abs(pdgEta)>0.6)continue;

			fCheckEtaMC->Fill(pdgEta);

      Int_t pdgMom = -99;
      Int_t labelMom = -1;
      Double_t pTmom = -1.0;
      //cout << "check Mother" << endl;
      FindMother(fMCparticle,labelMom,pdgMom,pTmom);
      if(pdgMom==-99 && iEnhance)iHijing = 0;  // particles from enhance
      if(pdgMom>0 && iEnhance)iHijing = -1;  // particles from enhance but feeddown
      //if(pdgGen==111)cout << "pdg = " << pdgGen << " ; enhance = " << iEnhance << " ; HIJIJG = " << iHijing << " ; mother = " << pdgMom  << endl;

      if(iHijing>-1)
        {
         if(pdgGen==111)fHistMCorgPi0->Fill(iHijing,pTtrue);
         if(pdgGen==221)fHistMCorgEta->Fill(iHijing,pTtrue);
        }

      if(TMath::Abs(pdgGen)!=11)continue;
      if(pTtrue<2.0)continue;


      //if(iHijing ==0)
      //if(pdgMom>0)
      if(pdgMom!=0)
       {
         AliAODMCParticle* fMCparticleMom = (AliAODMCParticle*) fMCarray->At(labelMom);
         //if(pdgMom==411 || pdgMom==421 || pdgMom==413 || pdgMom==423 || pdgMom==431 || pdgMom==433)
				 if(IsDdecay(pdgMom))
            {
             fHistMCorgD->Fill(fMCparticle->Pt());
             //cout << "orgD : " << pdgMom << " ; " << pdgGen << endl;
            }
         //if(pdgMom==511 || pdgMom==521 || pdgMom==513 || pdgMom==523 || pdgMom==531 || pdgMom==533)
				 if(IsBdecay(pdgMom))
           {
            fHistMCorgB->Fill(fMCparticle->Pt());
            //cout << "orgB : " << pdgMom << " ; " << pdgGen << endl;
           }
        }

     }

 return;
}
