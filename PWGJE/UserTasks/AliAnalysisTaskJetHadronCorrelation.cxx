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
#include "AliAODJet.h"
#include "AliAODJetEventBackground.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAODHeader.h"
#include "AliAODMCHeader.h"
//#include "AliGenPythiaEventHeader.h"
#include "AliJetKineReaderHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliGenDPMjetEventHeader.h"



#include "AliAnalysisTaskJetHadronCorrelation.h"
#include "AliAnalysisTaskPhiCorrelations.h"
//#include "AliAnalysisHelperJetTasks.h"
#include "AliPWG4HighPtQAMC.h"




ClassImp(AliAnalysisTaskJetHadronCorrelation)

				//________________________________________________________________________
				AliAnalysisTaskJetHadronCorrelation::AliAnalysisTaskJetHadronCorrelation():
								AliAnalysisTaskSE(),
								fUseAODInput(kFALSE),
								fFillAOD(kFALSE),
								fJetBranch("jets"),
								fNonStdFile(""),
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
                fxsec(0.),
                ftrial(1.),
                fJetRecEtaWindow(0.5),       // eta window for rec jets
                fMinJetPt(10), 
								fHistList(0x0), // Output list
								fIfiles(0),
								fH1Events(0x0),
								fH1Xsec(0x0),
								fH1Trials(0x0),
								fH1JetAKT04_pt                (0x0),
								fH1leadJetAKT04_pt            (0x0),
								fH1leadJetAKT04_pt_dijet      (0x0),
								fH1subJetAKT04_pt_dijet       (0x0),
								fH2JetsJetAKT04_dphi          (0x0),
								fH2JetsJetAKT04_deta          (0x0),
								fH2JetsJetAKT04_Aj            (0x0),
								fH2JetsJetAKT04_pt            (0x0),
								fH1JetMCAKT04_pt              (0x0),
								fH1leadJetMCAKT04_pt          (0x0),
								fH1leadJetMCAKT04_pt_dijet    (0x0),
								fH1subJetMCAKT04_pt_dijet     (0x0),
								fH2JetsJetMCAKT04_dphi        (0x0),
								fH2JetsJetMCAKT04_deta        (0x0),
								fH2JetsJetMCAKT04_Aj          (0x0),
								fH2JetsJetMCAKT04_pt          (0x0)


{


				if(IsMC){
								for(int j=0;j<5;j++){
												fH1AKT04_ndiJ_ediv[j]=0;
												fH1leadJetMCAKT04_dphiResolution          [j]=0;
												fH1subJetMCAKT04_dphiResolution           [j]=0;
												for(int k=0;k<5;k++){
																fH1JetHadronAKT04_dphi_ediv               [j][k]=0;
																fH1JetHadronAKT04_dphi_tptweight_ediv     [j][k]=0;
																fH1JetHadronAKT04_dphi_tJptweight_ediv    [j][k]=0;
												}
								}
				}else{
								for(int j=0;j<5;j++){
												fH1AKT04_ndiJ_ediv[j]=0;
												for(int k=0;k<5;k++){
																fH1JetHadronAKT04_dphi_ediv             [j][k]=0;
																fH1JetHadronAKT04_dphi_tptweight_ediv   [j][k]=0;
																fH1JetHadronAKT04_dphi_tJptweight_ediv  [j][k]=0;
												}
								}
				}


				// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskJetHadronCorrelation::AliAnalysisTaskJetHadronCorrelation(const char *name):
				AliAnalysisTaskSE(name),
				fUseAODInput(kFALSE),
				fFillAOD(kFALSE),
				fJetBranch("jets"),
				fNonStdFile(""),
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
				fxsec(0.),
				ftrial(1.),
				fJetRecEtaWindow(0.5),       // eta window for rec jets
				fMinJetPt(10), 
				fHistList(0x0), // Output list
				fIfiles(0),

				fH1Events(0x0),
				fH1Xsec(0x0),
				fH1Trials(0x0),
				fH1JetAKT04_pt                (0x0),
				fH1leadJetAKT04_pt            (0x0),
				fH1leadJetAKT04_pt_dijet      (0x0),
				fH1subJetAKT04_pt_dijet       (0x0),
				fH2JetsJetAKT04_dphi          (0x0),
				fH2JetsJetAKT04_deta          (0x0),
				fH2JetsJetAKT04_Aj            (0x0),
				fH2JetsJetAKT04_pt            (0x0),
				fH1JetMCAKT04_pt              (0x0),
				fH1leadJetMCAKT04_pt          (0x0),
				fH1leadJetMCAKT04_pt_dijet    (0x0),
				fH1subJetMCAKT04_pt_dijet     (0x0),
				fH2JetsJetMCAKT04_dphi        (0x0),
				fH2JetsJetMCAKT04_deta        (0x0),
				fH2JetsJetMCAKT04_Aj          (0x0),
				fH2JetsJetMCAKT04_pt          (0x0)


{
				if(IsMC){
								for(int j=0;j<5;j++){
												fH1AKT04_ndiJ_ediv[j]=0;
												fH1leadJetMCAKT04_dphiResolution          [j]=0;
												fH1subJetMCAKT04_dphiResolution           [j]=0;
												for(int k=0;k<5;k++){
																fH1JetHadronAKT04_dphi_ediv               [j][k]=0;
																fH1JetHadronAKT04_dphi_tptweight_ediv     [j][k]=0;
																fH1JetHadronAKT04_dphi_tJptweight_ediv    [j][k]=0;
												}
								}
				}else{
								for(int j=0;j<5;j++){
												fH1AKT04_ndiJ_ediv[j]=0;
												for(int k=0;k<5;k++){
																fH1JetHadronAKT04_dphi_ediv             [j][k]=0;
																fH1JetHadronAKT04_dphi_tptweight_ediv   [j][k]=0;
																fH1JetHadronAKT04_dphi_tJptweight_ediv  [j][k]=0;
												}
								}
				}

				// Default constructor

				// Constructor

				// Define input and output slots here
				// Input slot #0 works with a TChain
				DefineInput(0, TChain::Class());
				// Output slot #0 id reserved by the base class for AOD
				// Output slot #1 writes into a TH1 container
				DefineOutput(1, TList::Class());

}

//________________________________________________________________________
void AliAnalysisTaskJetHadronCorrelation::UserCreateOutputObjects()
{
				// Create histograms
				// Called once


				fHistList = new TList();fHistList->SetOwner(kTRUE); cout<<"TList is created for output "<<endl;
				//if (!fHistList){ fHistList = new TList();fHistList->SetOwner(kTRUE); cout<<"TList is created for output "<<endl;}

				Bool_t oldStatus = TH1::AddDirectoryStatus();
				TH1::AddDirectory(kFALSE);

				Float_t pi=TMath::Pi();
				//gStyle->SetPalette(1);


				char *histname;
				if(IsMC){
								fH1Xsec                         = new TProfile("Xsec"               ,"Xsec"                   ,1,0,1);
								fH1Trials                       = new TH1F    ("Trials"             ,"Trials"                 ,1,0,1);
								fH1JetMCAKT04_pt                = new TH1F("JetMCAKT04_pt"          ,"JetMCAKT04_pt"          ,400,0,400);
								fH1leadJetMCAKT04_pt            = new TH1F("leadJetMCAKT04_pt"      ,"leadJetMCAKT04_pt"      ,400,0,400);
								fH1leadJetMCAKT04_pt_dijet      = new TH1F("leadJetMCAKT04_pt_dijet","leadJetMCAKT04_pt_dijet",400,0,400);
								fH1subJetMCAKT04_pt_dijet       = new TH1F("subJetMCAKT04_pt_dijet" ,"subJetMCAKT04_pt_dijet" ,400,0,400);
								fH1JetAKT04_pt                  = new TH1F("JetAKT04_pt"            ,"JetAKT04_pt"            ,400,0,400);
								fH1leadJetAKT04_pt              = new TH1F("leadJetAKT04_pt"        ,"leadJetAKT04_pt"        ,400,0,400);
								fH1leadJetAKT04_pt_dijet        = new TH1F("leadJetAKT04_pt_dijet"  ,"leadJetAKT04_pt_dijet"  ,400,0,400);
								fH1subJetAKT04_pt_dijet         = new TH1F("subJetAKT04_pt_dijet"   ,"subJetAKT04_pt_dijet"   ,400,0,400);
								histname = Form("JetsJetMCAKT04_dphi");
								fH2JetsJetMCAKT04_dphi          = new TH2F(histname,histname,200,0,400,100,-2*pi,2*pi);
								histname = Form("JetsJetMCAKT04_deta");
								fH2JetsJetMCAKT04_deta          = new TH2F(histname,histname,200,0,400,100,-1.5,1.5);
								histname = Form("JetsJetMCAKT04_Aj");
								fH2JetsJetMCAKT04_Aj            = new TH2F(histname,histname,200,0,400,100,0,1.2);
								histname = Form("JetsJetMCAKT04_pt");
								fH2JetsJetMCAKT04_pt            = new TH2F(histname,histname,200,0,400,200,0,400);
								histname = Form("JetsJetAKT04_dphi");
								fH2JetsJetAKT04_dphi            = new TH2F(histname,histname,200,0,400,100,-2*pi,2*pi);
								histname = Form("JetsJetAKT04_deta");
								fH2JetsJetAKT04_deta            = new TH2F(histname,histname,200,0,400,100,-1.5,1.5);
								histname = Form("JetsJetAKT04_Aj");
								fH2JetsJetAKT04_Aj              = new TH2F(histname,histname,200,0,400,100,0,1.2);
								histname = Form("JetsJetAKT04_pt");
								fH2JetsJetAKT04_pt              = new TH2F(histname,histname,200,0,400,200,0,400);
								for(int j=0;j<5;j++){
												histname = Form("AKT04_ndiJ_ediv%d",j);
												fH1AKT04_ndiJ_ediv[j]= new TH1F(histname,histname,1,1,2);
												histname = Form("leadJetMCAKT04_dphiResolution%d",j);
												fH1leadJetMCAKT04_dphiResolution[j] = new TH1F(histname,histname,200,-2*pi,2*pi);
												histname = Form("subJetMCAKT04_dphiResolution%d",j);
												fH1subJetMCAKT04_dphiResolution[j] = new TH1F(histname,histname,200,-2*pi,2*pi);
												for(int k=0;k<5;k++){
																histname = Form("JetHadronAKT04_dphi_ediv%d%d",j,k);
																fH1JetHadronAKT04_dphi_ediv             [j][k]= new TH1F(histname,histname,200,-2*pi,2*pi);
																histname = Form("JetHadronAKT04_dphi_tptweight_ediv%d%d",j,k);
																fH1JetHadronAKT04_dphi_tptweight_ediv   [j][k]= new TH1F(histname,histname,200,-2*pi,2*pi);
																histname = Form("JetHadronAKT04_dphi_tJptweight_ediv%d%d",j,k);
																fH1JetHadronAKT04_dphi_tJptweight_ediv  [j][k]= new TH1F(histname,histname,200,-2*pi,2*pi);
												}
								}
								fHistList->Add(fH1Xsec);
								fHistList->Add(fH1Trials);
								fHistList->Add(fH1JetMCAKT04_pt          );
								fHistList->Add(fH1leadJetMCAKT04_pt      );
								fHistList->Add(fH1leadJetMCAKT04_pt_dijet);
								fHistList->Add(fH1subJetMCAKT04_pt_dijet );
								fHistList->Add(fH1JetAKT04_pt            );
								fHistList->Add(fH1leadJetAKT04_pt        );
								fHistList->Add(fH1leadJetAKT04_pt_dijet  );
								fHistList->Add(fH1subJetAKT04_pt_dijet   );
								fHistList->Add(fH2JetsJetMCAKT04_dphi);
								fHistList->Add(fH2JetsJetMCAKT04_deta);
								fHistList->Add(fH2JetsJetMCAKT04_Aj  );
								fHistList->Add(fH2JetsJetMCAKT04_pt  );
								fHistList->Add(fH2JetsJetAKT04_dphi  );
								fHistList->Add(fH2JetsJetAKT04_deta  );
								fHistList->Add(fH2JetsJetAKT04_Aj    );
								fHistList->Add(fH2JetsJetAKT04_pt    );
								for(int j=0;j<5;j++){
												fHistList->Add(fH1AKT04_ndiJ_ediv    [j]);
												fHistList->Add(fH1leadJetMCAKT04_dphiResolution          [j]);
												fHistList->Add(fH1subJetMCAKT04_dphiResolution           [j]);
												for(int k=0;k<5;k++){
																fHistList->Add(fH1JetHadronAKT04_dphi_ediv               [j][k]);
																fHistList->Add(fH1JetHadronAKT04_dphi_tptweight_ediv     [j][k]);
																fHistList->Add(fH1JetHadronAKT04_dphi_tJptweight_ediv    [j][k]);
												}
								}
				}
				else{
								fH1Events                     = new TH1F("Events"                  ,"Events"                  ,1,0,1);
								fH1JetAKT04_pt                = new TH1F("JetAKT04_pt"             ,"JetAKT04_pt"             ,400,0,400);
								fH1leadJetAKT04_pt            = new TH1F("leadJetAKT04_pt"         ,"leadJetAKT04_pt"         ,400,0,400);
								fH1leadJetAKT04_pt_dijet      = new TH1F("leadJetAKT04_pt_dijet"   ,"leadJetAKT04_pt_dijet"   ,400,0,400);
								fH1subJetAKT04_pt_dijet       = new TH1F("subJetAKT04_pt_dijet"    ,"subJetAKT04_pt_dijet"    ,400,0,400);
								histname = Form("JetsJetAKT04_dphi");
								fH2JetsJetAKT04_dphi          = new TH2F(histname,histname,200,0,400,100,-2*pi,2*pi);
								histname = Form("JetsJetAKT04_deta");
								fH2JetsJetAKT04_deta          = new TH2F(histname,histname,200,0,400,100,-1.5,1.5);
								histname = Form("JetsJetAKT04_Aj");
								fH2JetsJetAKT04_Aj            = new TH2F(histname,histname,200,0,400,100,0,1.2);
								histname = Form("JetsJetAKT04_pt");
								fH2JetsJetAKT04_pt            = new TH2F(histname,histname,200,0,400,200,0,400);
								for(int j=0;j<5;j++){
												histname = Form("AKT04_ndiJ_ediv%d",j);
												fH1AKT04_ndiJ_ediv[j]= new TH1F(histname,histname,1,1,2);
												for(int k=0;k<5;k++){
																histname = Form("JetHadronAKT04_dphi_ediv%d%d",j,k);
																fH1JetHadronAKT04_dphi_ediv             [j][k]= new TH1F(histname,histname,200,-2*pi,2*pi);
																histname = Form("JetHadronAKT04_dphi_tptweight_ediv%d%d",j,k);
																fH1JetHadronAKT04_dphi_tptweight_ediv   [j][k]= new TH1F(histname,histname,200,-2*pi,2*pi);
																histname = Form("JetHadronAKT04_dphi_tJptweight_ediv%d%d",j,k);
																fH1JetHadronAKT04_dphi_tJptweight_ediv  [j][k]= new TH1F(histname,histname,200,-2*pi,2*pi);
												}
								}
								fHistList->Add(fH1Events);
								fHistList->Add(fH1JetAKT04_pt          );
								fHistList->Add(fH1leadJetAKT04_pt      );
								fHistList->Add(fH1leadJetAKT04_pt_dijet);
								fHistList->Add(fH1subJetAKT04_pt_dijet );
								fHistList->Add(fH2JetsJetAKT04_dphi  );
								fHistList->Add(fH2JetsJetAKT04_deta  );
								fHistList->Add(fH2JetsJetAKT04_Aj    );
								fHistList->Add(fH2JetsJetAKT04_pt    );
								for(int j=0;j<5;j++){
												fHistList->Add(fH1AKT04_ndiJ_ediv    [j]);
												for(int k=0;k<5;k++){
																fHistList->Add(fH1JetHadronAKT04_dphi_ediv             [j][k]);
																fHistList->Add(fH1JetHadronAKT04_dphi_tptweight_ediv   [j][k]);
																fHistList->Add(fH1JetHadronAKT04_dphi_tJptweight_ediv  [j][k]);
												}
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
void AliAnalysisTaskJetHadronCorrelation::Init()
{
				// Initialization                                                                                                    
				if (fDebug) printf("AnalysisTaskJetHadronCorrelation::Init() \n");

}

Bool_t AliAnalysisTaskJetHadronCorrelation::Notify()
{


				fIfiles++;
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
				fxsec=0;
				ftrial=1;

				if(tree){
								TFile *curfile = tree->GetCurrentFile();
								if(!curfile){
												Error("Notify","No current file");
												return kFALSE;
								}


								if(IsMC){
												AliPWG4HighPtQAMC::PythiaInfoFromFile(curfile->GetName(),fxsec,ftrial);
												//cout<<" Xsec "<<fxsec<<" trial "<<ftrial<<endl;
												fH1Xsec  ->Fill(0.,fxsec);
												fH1Trials->Fill(0.,ftrial);
								}else{
												Float_t totalEvent;
												totalEvent = GetTotalEvents(curfile->GetName());
												fH1Events->Fill(0.,totalEvent);
								}

				}

				printf("Reading File %s ",fInputHandler->GetTree()->GetCurrentFile()->GetName());
				return kTRUE;
}
void AliAnalysisTaskJetHadronCorrelation::FinishTaskOutput()
{
}



//________________________________________________________________________
void AliAnalysisTaskJetHadronCorrelation::UserExec(Option_t *) 
{


				// Main loop (called each event)
				// Execute analysis for current event

				// start jet analysis

				Double_t Jet_n  [20];
				Double_t Jet_pt [20][1000];
				Double_t Jet_eta[20][1000];
				Double_t Jet_phi[20][1000];
				Double_t subJet_n  [20];
				Double_t subJet_pt [20][1000];
				Double_t subJet_eta[20][1000];
				Double_t subJet_phi[20][1000];
				Double_t Track_n  ;
				Double_t Track_pt [1000];
				Double_t Track_eta[1000];
				Double_t Track_phi[1000];

				Track_n=0;
				for(int i=0;i<20;i++){
								Jet_n[i]=0;
								subJet_n[i]=0;
								for(int j=0;j<1000;j++){
												Jet_pt[i][j]=0.;
												Jet_phi[i][j]=999.;
												Jet_eta[i][j]=999.;
												subJet_pt[i][j]=0.;
												subJet_phi[i][j]=999.;
												subJet_eta[i][j]=999.;
												Track_pt [j]=0.;
												Track_phi[j]=999.;
												Track_eta[j]=999.;
								}
				}

				fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
				if (!fAODIn) {
								Printf("ERROR: fAODIn not available");
								return;
				}

				//////-----------------------------------------------------------------------------------


				int nLJetAOD=999; double ptLJetAOD=0;double phiLJetAOD=999.;double etaLJetAOD=999.;int nsLJetAOD=900;double ptsLJetAOD=0;double phisLJetAOD=900.;double etasLJetAOD=900.;
				int nLJetMC2=999; double ptLJetMC2=0;double phiLJetMC2=999.;double etaLJetMC2=999.;int nsLJetMC2=900;double ptsLJetMC2=0;double phisLJetMC2=900.;double etasLJetMC2=900.;
				int nLJetMC =999; double ptLJetMC =0;double phiLJetMC =999.;double etaLJetMC =999.;int nsLJetMC =900;double ptsLJetMC =0;double phisLJetMC =900.;double etasLJetMC =900.;
				bool findLJetAOD=false;bool findsLJetAOD=false;bool findsLJetAOD_temp=false;
				bool findLJetMC2=false;bool findsLJetMC2=false;bool findsLJetMC2_temp=false;
				bool findLJetMC =false;bool findsLJetMC =false;bool findsLJetMC_temp =false;
				bool findsLJet=false;
				int nLJet = 999;


				TString cAdd = "";
				cAdd += Form("%02d_",(int)((Radius+0.01)*10.));
				cAdd += Form("B%d",(int)BackM);
				cAdd += Form("_Filter%05d",Filtermask);
				cAdd += Form("_Cut%05d",(int)(1000.*TrackPtcut));
				cAdd += Form("_Skip%02d",SkipCone);
				TString Branchname_gen,Branchname_gen2,Branchname_rec;
				Branchname_gen  = Form("clustersMCKINE_%s%s",JFAlg.Data(),cAdd.Data());
				Branchname_gen2 = Form("clustersMCKINE2_%s%s",JFAlg.Data(),cAdd.Data());
				Branchname_rec  = Form("clustersAOD_%s%s",JFAlg.Data(),cAdd.Data());


				for(int algorithm=0;algorithm<3;algorithm++){
								//for LHC11a1  LHC11a2
								//if(algorithm==0)fJetBranch   = "clustersAOD_ANTIKT04_B0_Filter00256_Cut00150_Skip00";
								//if(algorithm==1)fJetBranch   = "clustersMCKINE2_ANTIKT04_B0_Filter00256_Cut00150_Skip00";
								//if(algorithm==2)fJetBranch   = "clustersMCKINE_ANTIKT04_B0_Filter00256_Cut00150_Skip00";
								if(algorithm==0)fJetBranch   = Branchname_rec.Data();
								if(algorithm==1)fJetBranch   = Branchname_gen2.Data();
								if(algorithm==2)fJetBranch   = Branchname_gen.Data();

								if((!IsMC&&(algorithm==1||algorithm==2)))continue;

								TClonesArray* jets = dynamic_cast <TClonesArray*> (fAODIn->FindListObject(fJetBranch.Data()));
								if(!jets)continue;
								Int_t nj = jets->GetEntriesFast();
								if (fDebug) printf("There are %5d jets in the event \n", nj);
								AliAODJet* jetsAOD;
								Jet_n[algorithm] = nj;
								int nLjet_in05_pthdiv[20][20];
								int nsLjet_in05_pthdiv[20][20];
								for(int i=0;i<20;i++){
												for(int j=0;j<20;j++){
																nLjet_in05_pthdiv [i][j]=0;
																nsLjet_in05_pthdiv[i][j]=0;
												}
								}
								//Find Leading Jet -------------------------------------------------------
								for(int njet =0;njet<nj;njet++){
												jetsAOD = (AliAODJet*) (jets->At(njet));
												Jet_pt   [algorithm][njet] = jetsAOD->Pt();
												Jet_phi  [algorithm][njet] = jetsAOD->Phi();  
												Jet_eta  [algorithm][njet] = jetsAOD->Eta();
												//TRefArray *reftracks = jetsAOD->GetRefTracks();
												double eta_cut_Jet=0.5;
												if((TMath::Abs(Jet_eta[algorithm][njet])<eta_cut_Jet)&&(Jet_pt[algorithm][njet]>10.)){
																if(algorithm==0){
																				 fH1JetAKT04_pt->Fill(Jet_pt[algorithm][njet]);  
																				if(Jet_pt[algorithm][njet]>ptLJetAOD){
																								findLJetAOD=true;
																								nLJetAOD=njet;ptLJetAOD=Jet_pt[algorithm][njet];phiLJetAOD=Jet_phi[algorithm][njet];etaLJetAOD=Jet_eta[algorithm][njet];
																				}
																}
																if(algorithm==1){
																				 fH1JetMCAKT04_pt->Fill(Jet_pt[algorithm][njet]);  
																				if(Jet_pt[algorithm][njet]>ptLJetMC2){
																								findLJetMC2=true;
																								nLJetMC2=njet;ptLJetMC2=Jet_pt[algorithm][njet];phiLJetMC2=Jet_phi[algorithm][njet];etaLJetMC2=Jet_eta[algorithm][njet];
																				}
																}
																if(algorithm==2){
																				if(Jet_pt[algorithm][njet]>ptLJetMC){
																								findLJetMC=true;
																								nLJetMC=njet;ptLJetMC=Jet_pt[algorithm][njet];phiLJetMC=Jet_phi[algorithm][njet];etaLJetMC=Jet_eta[algorithm][njet];
																				}
																}
												}
								}//njet loop
								//Leading Jet -----------------------------------------------------------

								if(algorithm==0){nLJet=nLJetAOD;fH1leadJetAKT04_pt->Fill(Jet_pt[algorithm][nLJet]);}
								if(algorithm==1){nLJet=nLJetMC2;fH1leadJetMCAKT04_pt->Fill(Jet_pt[algorithm][nLJet]);}
								if(algorithm==2){nLJet=nLJetMC2;}

								if(nj<2)continue;
								//Find Sub leading Jet ==================================================
								for(int njet=0;njet<nj;njet++){
												if(njet==nLJet)continue;
												jetsAOD = (AliAODJet *)jets->At(njet);
												subJet_pt [algorithm][njet] = jetsAOD->Pt();
												subJet_phi[algorithm][njet] = jetsAOD->Phi();
												subJet_eta[algorithm][njet] = jetsAOD->Eta();
												double eta_cut_Jet=0.5;
												if((TMath::Abs(subJet_eta[algorithm][njet])<eta_cut_Jet) && (subJet_pt[algorithm][njet]>10.)){
																if(subJet_pt[algorithm][njet]>ptsLJetAOD&&algorithm==0){
																				findsLJetAOD_temp=true;
																				nsLJetAOD=njet;ptsLJetAOD=Jet_pt[algorithm][njet];phisLJetAOD=Jet_phi[algorithm][njet];etasLJetAOD=Jet_eta[algorithm][njet];
																}
																if(subJet_pt[algorithm][njet]>ptsLJetMC2 &&algorithm==1){
																				findsLJetMC2_temp=true;
																				nsLJetMC2=njet;ptsLJetMC2=Jet_pt[algorithm][njet];phisLJetMC2=Jet_phi[algorithm][njet];etasLJetMC2=Jet_eta[algorithm][njet];
																}
																if(subJet_pt[algorithm][njet]>ptsLJetMC &&algorithm==2){
																				findsLJetMC_temp=true;
																				nsLJetMC =njet;ptsLJetMC =Jet_pt[algorithm][njet];phisLJetMC =Jet_phi[algorithm][njet];etasLJetMC =Jet_eta[algorithm][njet];
																}
												}
								}
								////Sub leading Jet ======================================================

								double Leading_pt=0.;double Leading_phi=999.;double Leading_eta=999.;double sLeading_pt=0.;double sLeading_phi=999.;double sLeading_eta=999.;
								if(algorithm==0){Leading_pt=ptLJetAOD;Leading_phi=phiLJetAOD;Leading_eta=etaLJetAOD;sLeading_pt=ptsLJetAOD;sLeading_phi=phisLJetAOD;sLeading_eta=etasLJetAOD;}
								if(algorithm==1){Leading_pt=ptLJetMC2;Leading_phi=phiLJetMC2;Leading_eta=etaLJetMC2;sLeading_pt=ptsLJetMC2;sLeading_phi=phisLJetMC2;sLeading_eta=etasLJetMC2;}
								if(algorithm==2){Leading_pt=ptLJetMC ;Leading_phi=phiLJetMC ;Leading_eta=etaLJetMC ;sLeading_pt=ptsLJetMC ;sLeading_phi=phisLJetMC ;sLeading_eta=etasLJetMC ;}

								////Di-Jet event trigger +++++++++++++++++++++++++++++++++++++++++++++++
								double DPhi = DeltaPhi(Leading_phi,sLeading_phi);
								double DEta = Leading_eta-sLeading_eta;
								if(algorithm==0){
												fH2JetsJetAKT04_dphi->Fill(Leading_pt,DPhi);
												fH2JetsJetAKT04_deta->Fill(Leading_pt,DEta);
								}
								if(algorithm==1){
												fH2JetsJetMCAKT04_dphi->Fill(Leading_pt,DPhi);
												fH2JetsJetMCAKT04_deta->Fill(Leading_pt,DEta);
								}
								if((TMath::Cos(DPhi)<-0.5)&&(Leading_pt>10.)&&(sLeading_pt>10.)){
												if(algorithm==0)findsLJetAOD=true;                                                         
												if(algorithm==1)findsLJetMC=true;                                                         

												double Aj = (Leading_pt-sLeading_pt)/(Leading_pt+sLeading_pt);
												if(algorithm==0){
																fH1leadJetAKT04_pt_dijet->Fill(Leading_pt);
																fH1subJetAKT04_pt_dijet ->Fill(sLeading_pt);
																fH2JetsJetAKT04_Aj->Fill(Leading_pt,Aj);
																fH2JetsJetAKT04_pt->Fill(Leading_pt,sLeading_pt);
												}
												if(algorithm==1){
																fH1leadJetMCAKT04_pt_dijet->Fill(Leading_pt);
																fH1subJetMCAKT04_pt_dijet ->Fill(sLeading_pt);
																fH2JetsJetMCAKT04_Aj->Fill(Leading_pt,Aj);
																fH2JetsJetMCAKT04_pt->Fill(Leading_pt,sLeading_pt);
												}
								}
								if(algorithm==0)findsLJet=findsLJetAOD;
								if(algorithm==1)findsLJet=findsLJetMC2;
								if(algorithm==2)findsLJet=findsLJetMC;
								////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

								if(algorithm>=2)continue;

								if((findsLJet)&&(Leading_pt>10.)&&(sLeading_pt>10.)){
												for(int eb=0;eb<5;eb++){//count number of Di-Jet in pt bin
																if(TMath::Abs(Leading_pt -20.*(eb+1))<10.){
																				if(algorithm==0)fH1AKT04_ndiJ_ediv[eb]->Fill(1);
																}
												}
												fJetBranch = "tracks";
												TClonesArray* tracks = dynamic_cast <TClonesArray*> (fAODIn->FindListObject(fJetBranch.Data()));
								
												if(!tracks)continue;
												Int_t nt = tracks->GetEntriesFast();
												AliAODTrack* trackAOD;
												Track_n = nt;
												for(int njet=0;njet<Jet_n[algorithm];njet++){
																if(njet!=nLJet)continue;
																double eta_cut_Jet=0.5;
																if(TMath::Abs(Jet_eta[algorithm][njet])<eta_cut_Jet){
																				for(int eb=0;eb<5;eb++){
																								if(TMath::Abs(Jet_pt[algorithm][njet] -20.*(eb+1))<10.){
																												for(int ntrack =0;ntrack<nt;ntrack++){
																																trackAOD = (AliAODTrack*) (tracks->At(ntrack));
																																if(trackAOD->TestFilterMask(Filtermask)){
																																				Track_pt   [ntrack]      = trackAOD->Pt();
																																				Track_phi  [ntrack]      = trackAOD->Phi();
																																				Track_eta  [ntrack]      = trackAOD->Eta();
																																				double DelPhi = DeltaPhi(Jet_phi[algorithm][njet],Track_phi[ntrack]);
																																				if(TMath::Abs(Track_eta[ntrack])<0.9){
																																								for(int teb=0;teb<5;teb++){
																																												if(teb==0){if(!(Track_pt[ntrack]>0.15))continue;}
																																												if(teb==1){if(!((Track_pt[ntrack]<1.5)&&(Track_pt[ntrack]>0.15)))continue;}
																																												if(teb==2){if(!((Track_pt[ntrack]<3.0)&&(Track_pt[ntrack]>1.5)))continue;}
																																												if(teb==3){if(!((Track_pt[ntrack]<4.5)&&(Track_pt[ntrack]>3.0)))continue;}
																																												if(teb==4){if(!(Track_pt[ntrack]>4.5))continue;}
																																												if(algorithm==0){
																																															fH1JetHadronAKT04_dphi_ediv                [eb][teb]->Fill(DelPhi); 
																																															fH1JetHadronAKT04_dphi_tptweight_ediv      [eb][teb]->Fill(DelPhi,Track_pt[ntrack]);
																																															fH1JetHadronAKT04_dphi_tJptweight_ediv     [eb][teb]->Fill(DelPhi,Track_pt[ntrack]/Jet_pt[algorithm][njet]);
																																												}
																																								}
																																				}
																																}
																												}//Track Loop
																								}
																				}
																}//eta cut
												}//Jet Loop
								}// Di-Jet
				}// algorithm LOOP
				if(IsMC){
								for(int eb=0;eb<5;eb++){
												double DPhi;
												if(TMath::Abs(ptLJetAOD -20.*(eb+1))<10.){
																DPhi = DeltaPhi(phiLJetMC,phiLJetAOD);
																fH1leadJetMCAKT04_dphiResolution[eb]->Fill(DPhi);
																DPhi = DeltaPhi(phisLJetMC,phisLJetAOD);
																fH1subJetMCAKT04_dphiResolution[eb]->Fill(DPhi);
												}
								}
				}





				PostData(1, fHistList);
				return;
}      

//________________________________________________________________________
void AliAnalysisTaskJetHadronCorrelation::Terminate(Option_t *) 
{
				// Terminate analysis
				if (fDebug) printf("AnalysisTaskPt: Terminate() \n");

}


Int_t  AliAnalysisTaskJetHadronCorrelation::GetListOfJets(TList *list,TClonesArray* jarray,Int_t type){

				if(fDebug>2)Printf("%s:%d Selecting jets with cuts %d",(char*)__FILE__,__LINE__,type);
				Int_t iCount = 0;

				if(!jarray){
								Printf("%s:%d no Jet array",(char*)__FILE__,__LINE__);
								return iCount;
				}


				for(int ij=0;ij<jarray->GetEntries();ij++){
								AliAODJet* jet = (AliAODJet*)jarray->At(ij);
								if(!jet)continue;
								if(type==0){
												// No cut at all, main purpose here is sorting      
												list->Add(jet);
												iCount++;
								}
								else if(type == 1){
												// eta cut
												if(JetSelected(jet)){
																list->Add(jet);
																iCount++;
												}
								}
				}

				list->Sort();
				return iCount;

}



Bool_t  AliAnalysisTaskJetHadronCorrelation::JetSelected(AliAODJet *jet){
				Bool_t selected = false;

				if(!jet)return selected;

				if(fabs(jet->Eta())<fJetRecEtaWindow&&jet->Pt()>fMinJetPt){
								selected = kTRUE;
				}
				return selected;

}


Double_t AliAnalysisTaskJetHadronCorrelation::DeltaPhi(Double_t phi1,Double_t phi2){
				Float_t pi=TMath::Pi();
				Double_t dphi = phi1-phi2;
				if     (dphi<(-1./2*pi))dphi = dphi +2*pi;
				else if(dphi>( 3./2*pi))dphi = dphi -2*pi;
				return dphi;
}

Float_t AliAnalysisTaskJetHadronCorrelation::GetTotalEvents(const char* currFile){
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


