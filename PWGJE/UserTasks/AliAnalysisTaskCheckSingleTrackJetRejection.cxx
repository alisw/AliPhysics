//#include <string.h>
//#include <TStyle.h>
#include <list>
#include <string>

#include "TTree.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliInputEventHandler.h"
#include "AliESDtrack.h"
#include "AliAODVertex.h"
#include "AliAODCluster.h"

#include <TROOT.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TRefArray.h>

#include "TDatabasePDG.h"
#include "AliAnalysisManager.h"
#include "AliJetFinder.h"
#include "AliJetHeader.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliUA1JetHeaderV1.h"
#include "AliSISConeJetHeader.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliVParticle.h"
#include "AliAODJet.h"
#include "AliAODJetEventBackground.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAODHeader.h"
#include "AliAODMCHeader.h"
#include "AliJetKineReaderHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliGenDPMjetEventHeader.h"



#include "AliAnalysisTaskCheckSingleTrackJetRejection.h"
#include "AliAnalysisTaskPhiCorrelations.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliPWG4HighPtQAMC.h"




ClassImp(AliAnalysisTaskCheckSingleTrackJetRejection)

				//________________________________________________________________________
				AliAnalysisTaskCheckSingleTrackJetRejection::AliAnalysisTaskCheckSingleTrackJetRejection(): 
								AliAnalysisTaskSE(),
								fUseAODInput(kFALSE),
								fFillAOD(kFALSE),
								fNonStdFile(""),
								fJetBranch("jets"),
								fAODIn(0x0),
								fAODOut(0x0),
								fAODExtension(0x0),
								JFAlg("ANTIKT"),
								Radius(0.4),
								Filtermask(256),
								BackM(0),
								TrackPtcut(0.15),
								SkipCone(0),
								IsMC(kTRUE),
								fHistList(0x0),
								fxsec(0.),
								ftrial(1.),
								fJetRecEtaWindow(0.5),
								fMinJetPt(0.15),
								fH1Xsec(0x0),
								fH1Trials(0x0),
								fH1Events(0x0),

								fH2jetMCAKT04_Jetpt_maxpt(0x0),
								fH2jetAKT04_Jetpt_maxpt  (0x0)

{
				for(int j=0;j<6;j++){
								fH1jetMCAKT04_pt         [j]=0;
								fH1jetAKT04_pt           [j]=0;

								fH2jetMCAKT04_Eratio       [j]=0;
								fH1jetMCAKT04_match        [j]=0;
				}

				// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskCheckSingleTrackJetRejection::AliAnalysisTaskCheckSingleTrackJetRejection(const char *name): 
				AliAnalysisTaskSE(name),
				fUseAODInput(kFALSE),
				fFillAOD(kFALSE),
				fNonStdFile(""),
				fJetBranch("jets"),
				fAODIn(0x0), 
				fAODOut(0x0), 
				fAODExtension(0x0),
				JFAlg("ANTIKT"),
				Radius(0.4),
				Filtermask(256),
				BackM(0),
				TrackPtcut(0.15),
				SkipCone(0),
				IsMC(kTRUE),
				fHistList(0x0),
				fxsec(0.),
				ftrial(1.),
				fJetRecEtaWindow(0.5),
				fMinJetPt(0.15),
				fH1Xsec(0x0),
				fH1Trials(0x0),
				fH1Events(0x0),
				fH2jetMCAKT04_Jetpt_maxpt(0x0),
				fH2jetAKT04_Jetpt_maxpt  (0x0)

{

				for(int j=0;j<6;j++){
								fH1jetMCAKT04_pt         [j]=0;
								fH1jetAKT04_pt           [j]=0;

								fH2jetMCAKT04_Eratio       [j]=0;
								fH1jetMCAKT04_match        [j]=0;
				}

				// Constructor
				// Define input and output slots here
				// Input slot #0 works with a TChain
				DefineInput(0, TChain::Class());
				// Output slot #0 id reserved by the base class for AOD
				// Output slot #1 writes into a TH1 container
				DefineOutput(1, TList::Class());

}

//________________________________________________________________________
void AliAnalysisTaskCheckSingleTrackJetRejection::UserCreateOutputObjects()
{
				// Create histograms
				// Called once

				if (!fHistList){ fHistList = new TList();fHistList->SetOwner(kTRUE); cout<<"TList is created for output "<<endl;}

				Bool_t oldStatus = TH1::AddDirectoryStatus();
				TH1::AddDirectory(kFALSE);

				char *histname;
				if(IsMC){
								fH1Xsec           = new TProfile("Xsec","Xsec",1,0,1);
								fH1Trials         = new TH1F    ("Trials","Trials",1,0,1);
								histname = Form("fH2jetMCAKT04_Jetpt_maxpt");
								fH2jetMCAKT04_Jetpt_maxpt = new TH2F(histname,histname,400,0,400,400,0,400);
								histname = Form("fH2jetAKT04_Jetpt_maxpt");
								fH2jetAKT04_Jetpt_maxpt = new TH2F(histname,histname,400,0,400,400,0,400);
								for(int j=0;j<6;j++){
												histname = Form("fH1jetMCAKT04_pt%d",j);
												fH1jetMCAKT04_pt[j] = new TH1F(histname,histname,200,0,200);
												histname = Form("fH2jetAKT04_pt%d",j);
												fH1jetAKT04_pt[j] = new TH1F(histname,histname,200,0,200);

												histname = Form("fH2jetMCAKT04_Eratio%d",j);
												fH2jetMCAKT04_Eratio[j] = new TH2F(histname,histname,200,0,200,100,0,30);
												histname = Form("fH1jetMCAKT04_match%d",j);
												fH1jetMCAKT04_match[j] = new TH1F(histname,histname,200,0,200);
								}
								fHistList->Add(fH1Xsec);
								fHistList->Add(fH1Trials);
								fHistList->Add(fH2jetMCAKT04_Jetpt_maxpt);
								fHistList->Add(fH2jetAKT04_Jetpt_maxpt);
								for(int j=0;j<6;j++){
												fHistList->Add(fH1jetAKT04_pt[j]);
												fHistList->Add(fH1jetMCAKT04_pt[j]);
												fHistList->Add(fH2jetMCAKT04_Eratio[j]);
												fHistList->Add(fH1jetMCAKT04_match[j]);
								}
				}
				else    {
								fH1Events         = new TH1F("Events","Number of Events",1,0,1);
								histname = Form("fH2jetAKT04_Jetpt_maxpt");
								fH2jetAKT04_Jetpt_maxpt = new TH2F(histname,histname,400,0,400,400,0,400);
								for(int j=0;j<6;j++){
												histname = Form("fH2jetAKT04_pt%d",j);
												fH1jetAKT04_pt[j] = new TH1F(histname,histname,200,0,200);
								}
								fHistList->Add(fH1Events);
								fHistList->Add(fH2jetAKT04_Jetpt_maxpt);
								for(int j=0;j<6;j++){
												fHistList->Add(fH1jetAKT04_pt[j]);
								}
				}


				// =========== Switch on Sumw2 for all histos ===========
				for (Int_t i=0; i<fHistList->GetEntries(); ++i) 
				{
								TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
								if (h1)
								{
												h1->Sumw2();
												continue;
								}
								THnSparse *hn = dynamic_cast<THnSparse*>(fHistList->At(i));
								if(hn)hn->Sumw2();
				}
				TH1::AddDirectory(oldStatus);


				PostData(1,fHistList);

}

//----------------------------------------------------------------------                                                 
void AliAnalysisTaskCheckSingleTrackJetRejection::Init()
{
				// Initialization                                                                                                    
				if (fDebug) printf("AnalysisTaskCheckSingleTrackJetRejection::Init() \n");

}

Bool_t AliAnalysisTaskCheckSingleTrackJetRejection::Notify()
{
				fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
				fAODOut = AODEvent();
				if(fNonStdFile.Length()!=0){
								AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
								fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
								if(fAODExtension){
												if(fDebug>1)Printf("AODExtension found for %s ",fNonStdFile.Data());
								}
				}

				TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
				fxsec=0.;
				ftrial=1.;

				if(tree){
								TFile *curfile = tree->GetCurrentFile();
								if(!curfile){
												Error("Notify","No current file");
												return kFALSE;
								}

								if(IsMC){
												PythiaInfoFromFile(curfile->GetName(),fxsec,ftrial);
												//cout<<" Xsec "<<fxsec<<" trial "<<ftrial<<endl;
												fH1Xsec  ->Fill(0.,fxsec);
												fH1Trials->Fill(0.,ftrial);
								}
								else{
												Float_t totalEvent;
												totalEvent = GetTotalEvents(curfile->GetName());
												fH1Events->Fill(0.,totalEvent);
								}

				}

				printf("Reading File %s ",fInputHandler->GetTree()->GetCurrentFile()->GetName());
				return kTRUE;
}

void AliAnalysisTaskCheckSingleTrackJetRejection::FinishTaskOutput()
{
}



//________________________________________________________________________
void AliAnalysisTaskCheckSingleTrackJetRejection::UserExec(Option_t *) 
{


				// Main loop (called each event)
				// Execute analysis for current event

				// start jet analysis

				Double_t Jet_n  [20];
				Double_t Jet_pt [20][10000];
				Double_t Jet_eta[20][10000];
				Double_t Jet_phi[20][10000];

				for(int i=0;i<20;i++){
								Jet_n[i]=0;
								for(int j=0;j<10000;j++){
												Jet_pt[i][j]=0.;
												Jet_phi[i][j]=999.;
												Jet_eta[i][j]=999.;
								}
				}

				fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
				if (!fAODIn) {
								Printf("ERROR: fAODIn not available");
								return;
				}

				TString cAdd = "";
				cAdd += Form("%02d_",(int)((Radius+0.01)*10.));
				cAdd += Form("B%d",(int)BackM);
				cAdd += Form("_Filter%05d",Filtermask);
				cAdd += Form("_Cut%05d",(int)(1000.*TrackPtcut));
				cAdd += Form("_Skip%02d",SkipCone);
				TString Branchname_gen,Branchname_rec;
				Branchname_gen = Form("clustersMCKINE2_%s%s",JFAlg.Data(),cAdd.Data());
				Branchname_rec = Form("clustersAOD_%s%s",JFAlg.Data(),cAdd.Data());


				bool fFIND[6][1000];
				double maxpt[2][1000];for(int i=0;i<1000;i++){maxpt[0][i]=0;maxpt[1][i]=0;}
				int nearrecID[1000];  for(int i=0;i<1000;i++){nearrecID[i]=99999;}
				AliAODJet* jetsAOD;
				for(int algorithm=0;algorithm<2;algorithm++){
								//for LHC11a1  LHC11a2 official
								if((!IsMC&&algorithm==0))continue;

								if(algorithm==0)fJetBranch   = Branchname_gen.Data();
								if(algorithm==1)fJetBranch   = Branchname_rec.Data();

								TClonesArray* jets = dynamic_cast <TClonesArray*> (fAODIn->FindListObject(fJetBranch.Data()));
								if(!jets)continue;
								Int_t nj = jets->GetEntriesFast();
								if (fDebug) printf("There are %5d jets in the event \n", nj);


								Jet_n[algorithm] = nj;
								for(int njet =0;njet<nj;njet++){
												jetsAOD = (AliAODJet*) (jets->At(njet));
												Jet_pt   [algorithm][njet] = jetsAOD->Pt();
												Jet_phi  [algorithm][njet] = jetsAOD->Phi();  
												Jet_eta  [algorithm][njet] = jetsAOD->Eta();
												double eta_cut_Jet=0.5;
												TRefArray *reftracks = jetsAOD->GetRefTracks();
												int ntracks=reftracks->GetEntriesFast();
												if(TMath::Abs(Jet_eta[algorithm][njet])<eta_cut_Jet){
																//------------calc max pt in Jet----------------
																double trackpt;
																double sumtrackpt=0;//test
																for(int ntr=0;ntr<ntracks;ntr++){// calc. max pt of track which is in Jet
																				AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(reftracks->At(ntr));
																				if(AODtrack){
																								if(AODtrack->TestFilterMask(256)){
																												trackpt = AODtrack->Pt();
																												sumtrackpt += trackpt;
																												if(trackpt>maxpt[algorithm][njet]){
																																maxpt[algorithm][njet] = trackpt;
																												}
																								}
																				}
																}// track Loop

																if(algorithm==0){
																				fH1jetMCAKT04_pt[0]->Fill(Jet_pt[algorithm][njet]);
																				if(maxpt[algorithm][njet]<40)fH1jetMCAKT04_pt[1]->Fill(Jet_pt[algorithm][njet]);
																				if(ntracks>1)     fH1jetMCAKT04_pt[2]->Fill(Jet_pt[algorithm][njet]);
																				if(ntracks>2)     fH1jetMCAKT04_pt[3]->Fill(Jet_pt[algorithm][njet]);
																				if(ntracks>3)     fH1jetMCAKT04_pt[4]->Fill(Jet_pt[algorithm][njet]);
																				if((maxpt[algorithm][njet]<40)&&(ntracks>1))fH1jetMCAKT04_pt[5]->Fill(Jet_pt[algorithm][njet]);
																				fH2jetMCAKT04_Jetpt_maxpt->Fill(maxpt[algorithm][njet],Jet_pt[algorithm][njet]);
																}
																if(algorithm==1){
																				fH1jetAKT04_pt[0]->Fill(Jet_pt[algorithm][njet]);
																				if(maxpt[algorithm][njet]<40)fH1jetAKT04_pt[1]->Fill(Jet_pt[algorithm][njet]);
																				if(ntracks>1)     fH1jetAKT04_pt[2]->Fill(Jet_pt[algorithm][njet]);
																				if(ntracks>2)     fH1jetAKT04_pt[3]->Fill(Jet_pt[algorithm][njet]);
																				if(ntracks>3)     fH1jetAKT04_pt[4]->Fill(Jet_pt[algorithm][njet]);
																				if((maxpt[algorithm][njet]<40)&&(ntracks>1))fH1jetAKT04_pt[5]->Fill(Jet_pt[algorithm][njet]);
																				fH2jetAKT04_Jetpt_maxpt->Fill(maxpt[algorithm][njet],Jet_pt[algorithm][njet]);
																}
												}
								}
								if(!(IsMC&&algorithm==1))continue;
								for(int njetMC =0;njetMC<Jet_n[0];njetMC++){
												double eta_cut_Jet=0.5;
												if(TMath::Abs(Jet_eta[0][njetMC])<eta_cut_Jet){
																//Find muched jet pare=====================================
																for(int cut=0;cut<6;cut++){
																				double min_R=10.;
																				for(int njetAOD=0;njetAOD<Jet_n[1];njetAOD++){
																				  fJetBranch   = "clustersAOD_ANTIKT04_B0_Filter00256_Cut00150_Skip00";
																				  jets = dynamic_cast <TClonesArray*> (fAODIn->FindListObject(fJetBranch.Data()));
																				  if(!jets)continue;
																				  jetsAOD = (AliAODJet*) (jets->At(njetAOD));
																				  TRefArray *reftracks = jetsAOD->GetRefTracks();
																				  int ntracks=reftracks->GetEntriesFast();
																								if(cut==1){if(maxpt[1][njetAOD]>=40.)continue;}
																								if(cut==2){if(ntracks==1)continue;}
																								if(cut==3){if(ntracks<=2)continue;}
																								if(cut==4){if(ntracks<=3)continue;}
																								if(cut==5){if(maxpt[1][njetAOD]>=40.)continue;if(ntracks==1)continue;}
																								double DelR = DeltaR(Jet_phi[0][njetMC],Jet_phi[1][njetAOD],Jet_eta[0][njetMC],Jet_eta[1][njetAOD]);
																								if(DelR<min_R){
																												nearrecID[njetMC]=njetAOD;
																												min_R=DelR;
																								}
																				}
																				if(min_R<0.4){
																								min_R=10.;
																								int neargenID=99999;
																								for(int njet =0;njet<Jet_n[0];njet++){
																												double DelR = DeltaR(Jet_phi[0][njet],Jet_phi[1][nearrecID[njetMC]],Jet_eta[0][njet],Jet_eta[1][nearrecID[njetMC]]);
																												if(DelR<min_R){
																																neargenID=njet;
																																min_R=DelR;
																												}
																								}
																								if((min_R<0.4)&&(neargenID==njetMC))fFIND[cut][njetMC]=true;
																				}
																}
																//======================================================
												}
								}
				}//algorithm
				if(IsMC){
								for(int njetMC =0;njetMC<Jet_n[0];njetMC++){
												double eta_cut_Jet=0.5;
												if(TMath::Abs(Jet_eta[0][njetMC])<eta_cut_Jet){
																for(int cut=0;cut<6;cut++){
																				if(fFIND[cut][njetMC]==true){
																								fH1jetMCAKT04_match[cut]->Fill(Jet_pt[0][njetMC]);
																								fH2jetMCAKT04_Eratio[cut]->Fill(Jet_pt[0][njetMC],Jet_pt[1][nearrecID[njetMC]]/Jet_pt[0][njetMC]);
																				}
																}
												}//eta window
								}//njet loop
				}

				// End of jet analysis ---------
				// Post output data.
				PostData(1, fHistList);
				return;
}      

//________________________________________________________________________
void AliAnalysisTaskCheckSingleTrackJetRejection::Terminate(Option_t *) 
{
				// Terminate analysis
				//
				if (fDebug) printf("AnalysisTaskPt: Terminate() \n");

}


Bool_t  AliAnalysisTaskCheckSingleTrackJetRejection::JetSelected(AliAODJet *jet){
				Bool_t selected = false;

				if(!jet)return selected;

				if(fabs(jet->Eta())<fJetRecEtaWindow&&jet->Pt()>fMinJetPt){
								selected = kTRUE;
				}
				return selected;

}


Double_t AliAnalysisTaskCheckSingleTrackJetRejection::DeltaPhi(Double_t phi1,Double_t phi2){
				Float_t pi=TMath::Pi();
				Double_t dphi = phi1-phi2;
				if     (dphi<(-1./2*pi))dphi = dphi +2*pi;
				else if(dphi>( 3./2*pi))dphi = dphi -2*pi;
				return dphi;
}
Double_t AliAnalysisTaskCheckSingleTrackJetRejection::DeltaR(Double_t phi1,Double_t phi2,Double_t eta1,Double_t eta2){
				Float_t pi=TMath::Pi();
				Double_t dphi = DeltaPhi(phi1,phi2);
				if     (dphi<(-1./2*pi))dphi = dphi +2*pi;
				else if(dphi>( 3./2*pi))dphi = dphi -2*pi;
				Double_t Deta = eta1 - eta2;
				Double_t deltaR = sqrt(pow(dphi,2)+pow(Deta,2));
				return deltaR;
}

Bool_t AliAnalysisTaskCheckSingleTrackJetRejection::PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials){
				//
				// get the cross section and the trails either from pyxsec.root or from pysec_hists.root
				// This is to called in Notify and should provide the path to the AOD/ESD file
				// Copied from AliAnalysisTaskJetSpectrum2
				//

				TString file(currFile);
				fXsec = 0;
				fTrials = 1;

				if(file.Contains("root_archive.zip#")){
								Ssiz_t pos1 = file.Index("root_archive",12,TString::kExact);
								Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
								file.Replace(pos+1,20,"");
				}
				else {
								// not an archive take the basename....
								file.ReplaceAll(gSystem->BaseName(file.Data()),"");
				}
				//  Printf("%s",file.Data());


				TFile *filexsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really 
				//test the existance of a file in a archive so we have to lvie with open error message from root
								if(!filexsec){
												// next trial fetch the histgram file
												filexsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
												if(!filexsec){
																// not a severe condition but inciate that we have no information
																return kFALSE;
												}
												else{
																// find the tlist we want to be independtent of the name so use the Tkey
																TKey* key = (TKey*)filexsec->GetListOfKeys()->At(0);
																if(!key){
																				filexsec->Close();
																				return kFALSE;
																}
																TList *list = dynamic_cast<TList*>(key->ReadObj());
																if(!list){
																				filexsec->Close();
																				return kFALSE;
																}
																fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
																fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
																filexsec->Close();
												}
								} // no tree pyxsec.root
								else {
												TTree *xtree = (TTree*)filexsec->Get("Xsection");
												if(!xtree){
																filexsec->Close();
																return kFALSE;
												}
												UInt_t   ntrials  = 0;
												Double_t  xsection  = 0;
												xtree->SetBranchAddress("xsection",&xsection);
												xtree->SetBranchAddress("ntrials",&ntrials);
												xtree->GetEntry(0);
												fTrials = ntrials;
												fXsec = xsection;
												filexsec->Close();
								}
				return kTRUE;
}


Float_t AliAnalysisTaskCheckSingleTrackJetRejection::GetTotalEvents(const char* currFile){
				Float_t totalevent;
				TString file_es(currFile);
				if(file_es.Contains("root_archive.zip#")){
								Ssiz_t pos1 = file_es.Index("root_archive",12,TString::kExact);
								Ssiz_t pos = file_es.Index("#",1,pos1,TString::kExact);
								file_es.Replace(pos+1,20,"");
				}
				else {
								// not an archive take the basename....
								file_es.ReplaceAll(gSystem->BaseName(file_es.Data()),"");
				}

				TString cAdd = "";
				cAdd += Form("%02d_",(int)((Radius+0.01)*10.));
				cAdd += Form("B%d",(int)BackM);
				cAdd += Form("_Filter%05d",Filtermask);
				cAdd += Form("_Cut%05d",(int)(1000.*TrackPtcut));
				cAdd += Form("_Skip%02d",SkipCone);
				TString Dirname,Listname;
				Dirname  = Form("PWG4_cluster_AOD__%s%s",JFAlg.Data(),cAdd.Data());
				Listname = Form("pwg4cluster_AOD__%s%s",JFAlg.Data(),cAdd.Data());

				TFile *feventstat = TFile::Open(Form("%s%s",file_es.Data(),"JetTasksOutput.root"));
				//gROOT->Cd("PWG4_cluster_AOD__KT06_B0_Filter00256_Cut00150_Skip00");
				//TList *templist     = (TList*)gROOT->FindObject("pwg4cluster_AOD__KT06_B0_Filter00256_Cut00150_Skip00");
				gROOT->Cd(Dirname.Data());
				TList *templist     = (TList*)gROOT->FindObject(Listname.Data());
				TH1F* temphist = (TH1F*)templist->FindObject("fh1Trials");
				totalevent = temphist->Integral();
				//cout<<temphist->Integral()<<endl;
				delete feventstat;
				return totalevent;

}


