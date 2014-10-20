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
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
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

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetHadronCorrelation)

				//________________________________________________________________________
				AliAnalysisTaskJetHadronCorrelation::AliAnalysisTaskJetHadronCorrelation():
								AliAnalysisTaskSE(),
								fUseAODInput(kFALSE),
								fJetBranch("jets"),
								fNonStdFile(""),
								fAODIn(0x0),
								fAODOut(0x0),
								fAODExtension(0x0),
								JFAlg("ANTIKT"),         
								Radius(0.4),
								Filtermask(272),
								BackM(0),
								TrackPtcut(0.15),
								SkipCone(0),
								IsMC(kTRUE),
								JetEScale(1.),
								TrackEScale(1.),
								fxsec(0.),
								ftrial(1.),
								fHistList(0x0), // Output list
  fIfiles(0),
								fH1Events(0x0),
								fH1Xsec(0x0),
								fH1Trials(0x0),
								fH1Track_pt              (0x0),
								fH1Track_phi             (0x0),
								fH1Track_eta             (0x0),
								fH1MCTrack_pt            (0x0),
								fH1MCTrack_phi           (0x0),
								fH1MCTrack_eta           (0x0),
								fH1MCPrimTrack_pt        (0x0),
								fH1MCPrimTrack_phi       (0x0),
								fH1MCPrimTrack_eta       (0x0),
								fH1Jet_pt                (0x0),
								fH1Jet_phi               (0x0),
								fH1Jet_eta               (0x0),
								fH1leadJet_pt            (0x0),
								fH1leadJet_pt_dijet      (0x0),
								fH1subJet_pt_dijet       (0x0),
								fH1JetMC_pt              (0x0),
								fH1leadJetMC_pt          (0x0),
								fH1leadJetMC_pt_dijet    (0x0),
								fH1subJetMC_pt_dijet     (0x0),
								fH2JetsJet_dphi          (0x0),
								fH2JetsJet_deta          (0x0),
								fH2JetsJet_Aj            (0x0),
								fH2JetsJet_pt            (0x0),
								fH2JetsJetMC_dphi        (0x0),
								fH2JetsJetMC_deta        (0x0),
								fH2JetsJetMC_Aj          (0x0),
								fH2JetsJetMC_pt          (0x0),
								fH2Mult_Mtrack           (0x0),
								fH2Mult_Mlead            (0x0),
								fH2Mult_Mjet             (0x0),
								fH2Mult_Njet             (0x0),
								fH2Mult_Aj               (0x0),
								fH2Mlead_Aj              (0x0),
								fH2Jet_pt_Mlead          (0x0),
								fH2Jet_pt_Munder         (0x0),
								fH2leadJetMCptResolution (0x0),
								fH2TrackMCptResolution   (0x0),
								fH2TrackMCptEfficiency   (0x0),
								fH2AjCorrelation_MCRec   (0x0),
								fH2MleadCorrelation_MCRec(0x0)
{
				for(int j=0;j<5;j++){
								fH1ndiJ_ediv                             [j]=0;
								fH1Aj                                    [j]=0;
								fH1Mlead                                 [j]=0;
								fH1leadJetMC_dphiResolution              [j]=0;
								fH1subJetMC_dphiResolution               [j]=0;
								fH1leadJetMC_Efficiency                  [j]=0;
								fH1subJetMC_Efficiency                   [j]=0;
								for(int k=0;k<5;k++){
												fH1JetHadron_dphi_ediv           [j][k]=0;
												fH1JetHadron_dphi_tptweight_ediv [j][k]=0;
												fH1JetHadron_dphi_tJptweight_ediv[j][k]=0;
												fH1JetHadronMC_dphi_ediv           [j][k]=0;
												fH1JetHadronMC_dphi_tptweight_ediv [j][k]=0;
												fH1JetHadronMC_dphi_tJptweight_ediv[j][k]=0;
												fH1JetHadronMCPrim_dphi_ediv           [j][k]=0;
												fH1JetHadronMCPrim_dphi_tptweight_ediv [j][k]=0;
												fH1JetHadronMCPrim_dphi_tJptweight_ediv[j][k]=0;
								}
				}
				for(int j=0;j<3;j++){
								fH1ndiJ_2040Mlead                               [j]=0;
								fH1ndiJ_2040Aj                                  [j]=0;
								for(int k=0;k<5;k++){
												fH1JetHadron_dphi_tptweight2040_Mleaddep[j][k]=0;
												fH1JetHadron_dphi_tptweight2040_Ajdep   [j][k]=0;
												fH1JetHadronMC_dphi_tptweight2040_Mleaddep[j][k]=0;
												fH1JetHadronMC_dphi_tptweight2040_Ajdep   [j][k]=0;
												fH1JetHadronMCPrim_dphi_tptweight2040_Mleaddep[j][k]=0;
												fH1JetHadronMCPrim_dphi_tptweight2040_Ajdep   [j][k]=0;
								}
				}
				// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskJetHadronCorrelation::AliAnalysisTaskJetHadronCorrelation(const char *name):
				AliAnalysisTaskSE(name),
				fUseAODInput(kFALSE),
				fJetBranch("jets"),
				fNonStdFile(""),
				fAODIn(0x0), 
				fAODOut(0x0), 
				fAODExtension(0x0),
				JFAlg("ANTIKT"),         
				Radius(0.4),
				Filtermask(272),
				BackM(0),
				TrackPtcut(0.15),
				SkipCone(0),
				IsMC(kTRUE),
				JetEScale(1.),
				TrackEScale(1.),
				fxsec(0.),
				ftrial(1.),
				fHistList(0x0), // Output list
				fIfiles(0),
				fH1Events(0x0),
				fH1Xsec(0x0),
				fH1Trials(0x0),
				fH1Track_pt              (0x0),
				fH1Track_phi             (0x0),
				fH1Track_eta             (0x0),
				fH1MCTrack_pt            (0x0),
				fH1MCTrack_phi           (0x0),
				fH1MCTrack_eta           (0x0),
				fH1MCPrimTrack_pt        (0x0),
				fH1MCPrimTrack_phi       (0x0),
				fH1MCPrimTrack_eta       (0x0),
				fH1Jet_pt                (0x0),
				fH1Jet_phi               (0x0),
				fH1Jet_eta               (0x0),
				fH1leadJet_pt            (0x0),
				fH1leadJet_pt_dijet      (0x0),
				fH1subJet_pt_dijet       (0x0),
				fH1JetMC_pt              (0x0),
				fH1leadJetMC_pt          (0x0),
				fH1leadJetMC_pt_dijet    (0x0),
				fH1subJetMC_pt_dijet     (0x0),
				fH2JetsJet_dphi          (0x0),
				fH2JetsJet_deta          (0x0),
				fH2JetsJet_Aj            (0x0),
				fH2JetsJet_pt            (0x0),
				fH2JetsJetMC_dphi        (0x0),
				fH2JetsJetMC_deta        (0x0),
				fH2JetsJetMC_Aj          (0x0),
				fH2JetsJetMC_pt          (0x0),
				fH2Mult_Mtrack           (0x0),
				fH2Mult_Mlead            (0x0),
				fH2Mult_Mjet             (0x0),
				fH2Mult_Njet             (0x0),
				fH2Mult_Aj               (0x0),
				fH2Mlead_Aj              (0x0),
				fH2Jet_pt_Mlead          (0x0),
				fH2Jet_pt_Munder         (0x0),
				fH2leadJetMCptResolution (0x0),
				fH2TrackMCptResolution   (0x0),
				fH2TrackMCptEfficiency   (0x0),
				fH2AjCorrelation_MCRec   (0x0),
				fH2MleadCorrelation_MCRec(0x0)
{

				for(int j=0;j<5;j++){
								fH1ndiJ_ediv                             [j]=0;
								fH1Aj                                    [j]=0;
								fH1Mlead                                 [j]=0;
								fH1leadJetMC_dphiResolution              [j]=0;
								fH1subJetMC_dphiResolution               [j]=0;
								fH1leadJetMC_Efficiency                  [j]=0;
								fH1subJetMC_Efficiency                   [j]=0;
								for(int k=0;k<5;k++){
												fH1JetHadron_dphi_ediv           [j][k]=0;
												fH1JetHadron_dphi_tptweight_ediv [j][k]=0;
												fH1JetHadron_dphi_tJptweight_ediv[j][k]=0;
												fH1JetHadronMC_dphi_ediv           [j][k]=0;
												fH1JetHadronMC_dphi_tptweight_ediv [j][k]=0;
												fH1JetHadronMC_dphi_tJptweight_ediv[j][k]=0;
												fH1JetHadronMCPrim_dphi_ediv           [j][k]=0;
												fH1JetHadronMCPrim_dphi_tptweight_ediv [j][k]=0;
												fH1JetHadronMCPrim_dphi_tJptweight_ediv[j][k]=0;
								}
				}
				for(int j=0;j<3;j++){
								fH1ndiJ_2040Mlead                               [j]=0;
								fH1ndiJ_2040Aj                                  [j]=0;
								for(int k=0;k<5;k++){
												fH1JetHadron_dphi_tptweight2040_Mleaddep[j][k]=0;
												fH1JetHadron_dphi_tptweight2040_Ajdep   [j][k]=0;
												fH1JetHadronMC_dphi_tptweight2040_Mleaddep[j][k]=0;
												fH1JetHadronMC_dphi_tptweight2040_Ajdep   [j][k]=0;
												fH1JetHadronMCPrim_dphi_tptweight2040_Mleaddep[j][k]=0;
												fH1JetHadronMCPrim_dphi_tptweight2040_Ajdep   [j][k]=0;
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


				char *histname;

				fH1Events                  = new TH1F    ("Events"        ,"Events"            ,1,0,1);
				fH1Xsec                    = new TProfile("Xsec"          ,"Xsec"              ,1,0,1);
				fH1Trials                  = new TH1F    ("Trials"        ,"Trials"            ,1,0,1);

				fH1Track_pt                = new TH1F("Track_pt"          ,"Track_pt"          ,200,0,200);
				fH1Track_phi               = new TH1F("Track_phi"         ,"Track_phi"         ,100,0,2*pi);
				fH1Track_eta               = new TH1F("Track_eta"         ,"Track_eta"         ,100,-1.,1);
				fH1MCTrack_pt              = new TH1F("MCTrack_pt"        ,"MCTrack_pt"        ,200,0,200);
				fH1MCTrack_phi             = new TH1F("MCTrack_phi"       ,"MCTrack_phi"       ,100,0,2*pi);
				fH1MCTrack_eta             = new TH1F("MCTrack_eta"       ,"MCTrack_eta"       ,100,-1.,1);
				fH1MCPrimTrack_pt          = new TH1F("MCPrimTrack_pt"    ,"MCPrimTrack_pt"    ,200,0,200);
				fH1MCPrimTrack_phi         = new TH1F("MCPrimTrack_phi"   ,"MCPrimTrack_phi"   ,100,0,2*pi);
				fH1MCPrimTrack_eta         = new TH1F("MCPrimTrack_eta"   ,"MCPrimTrack_eta"   ,100,-1.,1);
				fH1Jet_pt                  = new TH1F("Jet_pt"            ,"Jet_pt"            ,200,0,200);
				fH1Jet_phi                 = new TH1F("Jet_phi"           ,"Jet_pt"            ,100,0,2*pi);
				fH1Jet_eta                 = new TH1F("Jet_eta"           ,"Jet_pt"            ,100,-1.,1);
				fH1leadJet_pt              = new TH1F("leadJet_pt"        ,"leadJet_pt"        ,200,0,200);
				fH1leadJet_pt_dijet        = new TH1F("leadJet_pt_dijet"  ,"leadJet_pt_dijet"  ,200,0,200);
				fH1subJet_pt_dijet         = new TH1F("subJet_pt_dijet"   ,"subJet_pt_dijet"   ,200,0,200);
				fH1JetMC_pt                = new TH1F("JetMC_pt"          ,"JetMC_pt"          ,200,0,200);
				fH1leadJetMC_pt            = new TH1F("leadJetMC_pt"      ,"leadJetMC_pt"      ,200,0,200);
				fH1leadJetMC_pt_dijet      = new TH1F("leadJetMC_pt_dijet","leadJetMC_pt_dijet",200,0,200);
				fH1subJetMC_pt_dijet       = new TH1F("subJetMC_pt_dijet" ,"subJetMC_pt_dijet" ,200,0,200);
				fH2JetsJet_dphi            = new TH2F("JetsJet_dphi"      ,"JetsJet_dphi"      ,200,0,200,100,-2*pi,2*pi);
				fH2JetsJet_deta            = new TH2F("JetsJet_deta"      ,"JetsJet_deta"      ,200,0,200,100,-1.5,1.5);
				fH2JetsJet_Aj              = new TH2F("JetsJet_Aj"        ,"JetsJet_Aj"        ,200,0,200,100,0,1.2);
				fH2JetsJet_pt              = new TH2F("JetsJet_pt"        ,"JetsJet_pt"        ,200,0,200,200,0,200);
				fH2JetsJetMC_dphi          = new TH2F("JetsJetMC_dphi"    ,"JetsJetMC_dphi"    ,200,0,200,100,-2*pi,2*pi);
				fH2JetsJetMC_deta          = new TH2F("JetsJetMC_deta"    ,"JetsJetMC_deta"    ,200,0,200,100,-1.5,1.5);
				fH2JetsJetMC_Aj            = new TH2F("JetsJetMC_Aj"      ,"JetsJetMC_Aj"      ,200,0,200,100,0,1.2);
				fH2JetsJetMC_pt            = new TH2F("JetsJetMC_pt"      ,"JetsJetMC_pt"      ,200,0,200,200,0,200);
				fH2Mult_Mtrack             = new TH2F("Mult_Mtrack"       ,"Mult_Mtrack"       ,50,0,250,50,0,250); 
				fH2Mult_Mlead              = new TH2F("Mult_Mlead"        ,"Mult_Mlead"        ,50,0,250,25,0,25);   
				fH2Mult_Mjet               = new TH2F("Mult_Mjet"         ,"Mult_Mjet"         ,50,0,250,50,0,100);  
				fH2Mult_Njet               = new TH2F("Mult_Njet"         ,"Mult_Njet"         ,50,0,250,50,0,50);   
				fH2Mult_Aj                 = new TH2F("Mult_Aj"           ,"Mult_Aj"           ,50,0,250,25,0,1.2);  
				fH2Mlead_Aj                = new TH2F("Mlead_Aj"          ,"Mlead_Aj"          ,25,0,25,25,0,1.2);   
				fH2Jet_pt_Mlead            = new TH2F("Jet_pt_Mlead"      ,"Jet_pt_Mlead"      ,50,0,200,25,0,25);   
				fH2Jet_pt_Munder           = new TH2F("Jet_pt_Munder"     ,"Jet_pt_Munder"     ,50,0,200,25,0,5);    
				fH2leadJetMCptResolution   = new TH2F("leadJetMCptResolution" ,"leadJetMCptResolution" ,200,0,200,200,0,200);    
				fH2TrackMCptResolution     = new TH2F("TrackMCptResolution"   ,"TrackMCptResolution"   ,200,0,200,200,0,200);    
				fH2TrackMCptEfficiency     = new TH2F("TrackMCptEfficiency"   ,"TrackMCptEfficiency"   ,200,0,200,100,0,1.2);    
				fH2AjCorrelation_MCRec     = new TH2F("AjCorrelation_MCRec"   ,"AjCorrelation_MCRec"   ,60,0,1.2,60,0,1.2);    
				fH2MleadCorrelation_MCRec  = new TH2F("MleadCorrelation_MCRec","MleadCorrelation_MCRec",60,0,60,60,0,60);    

				for(int j=0;j<5;j++){
								histname = Form("ndiJ_ediv%d",j);
								fH1ndiJ_ediv[j]= new TH1F(histname,histname,1,1,2);
								histname        = Form("Aj%d",j);                    
								fH1Aj[j]    = new TH1F(histname,histname,50,0,1.2);  
								histname        = Form("Mlead%d",j);                 
								fH1Mlead[j] = new TH1F(histname,histname,50,0,50);   
								histname = Form("leadJetMC_dphiResolution%d",j);
								fH1leadJetMC_dphiResolution[j] = new TH1F(histname,histname,200,-2*pi,2*pi);
								histname = Form("subJetMC_dphiResolution%d",j);
								fH1subJetMC_dphiResolution[j] = new TH1F(histname,histname,200,-2*pi,2*pi);
								histname = Form("leadJetMC_Efficiency%d",j);
								fH1leadJetMC_Efficiency[j] = new TH1F(histname,histname,100,0,1.2);
								histname = Form("subJetMC_Efficiency%d",j);
								fH1subJetMC_Efficiency[j] = new TH1F(histname,histname,100,0,1.2);
								for(int k=0;k<5;k++){
												histname = Form("JetHadron_dphi_ediv%d%d",j,k);
												fH1JetHadron_dphi_ediv             [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadron_dphi_tptweight_ediv%d%d",j,k);
												fH1JetHadron_dphi_tptweight_ediv   [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadron_dphi_tJptweight_ediv%d%d",j,k);
												fH1JetHadron_dphi_tJptweight_ediv  [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMC_dphi_ediv%d%d",j,k);
												fH1JetHadronMC_dphi_ediv             [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMC_dphi_tptweight_ediv%d%d",j,k);
												fH1JetHadronMC_dphi_tptweight_ediv   [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMC_dphi_tJptweight_ediv%d%d",j,k);
												fH1JetHadronMC_dphi_tJptweight_ediv  [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMCPrim_dphi_ediv%d%d",j,k);
												fH1JetHadronMCPrim_dphi_ediv             [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMCPrim_dphi_tptweight_ediv%d%d",j,k);
												fH1JetHadronMCPrim_dphi_tptweight_ediv   [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMCPrim_dphi_tJptweight_ediv%d%d",j,k);
												fH1JetHadronMCPrim_dphi_tJptweight_ediv  [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
								}
				}
				for(int j=0;j<3;j++){
								histname = Form("ndiJ_2040Mlead%d",j);
								fH1ndiJ_2040Mlead[j]= new TH1F(histname,histname,1,1,2);
								histname = Form("ndiJ_2040Aj%d",j);
								fH1ndiJ_2040Aj[j]= new TH1F(histname,histname,1,1,2);
								for(int k=0;k<5;k++){
												histname = Form("JetHadron_dphi_tptweight2040_Mleaddep%d%d",j,k);
												fH1JetHadron_dphi_tptweight2040_Mleaddep   [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadron_dphi_tptweight2040_Ajdep%d%d",j,k);
												fH1JetHadron_dphi_tptweight2040_Ajdep      [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMC_dphi_tptweight2040_Mleaddep%d%d",j,k);
												fH1JetHadronMC_dphi_tptweight2040_Mleaddep   [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMC_dphi_tptweight2040_Ajdep%d%d",j,k);
												fH1JetHadronMC_dphi_tptweight2040_Ajdep      [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMCPrim_dphi_tptweight2040_Mleaddep%d%d",j,k);
												fH1JetHadronMCPrim_dphi_tptweight2040_Mleaddep   [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
												histname = Form("JetHadronMCPrim_dphi_tptweight2040_Ajdep%d%d",j,k);
												fH1JetHadronMCPrim_dphi_tptweight2040_Ajdep      [j][k]= new TH1F(histname,histname,200,-1./2.*pi,3./2.*pi);
								}
				}


				if(IsMC){
								fHistList->Add(fH1Events            );
								fHistList->Add(fH1Xsec              );
								fHistList->Add(fH1Trials            );
								fHistList->Add(fH1Track_pt          );
								fHistList->Add(fH1Track_phi         );
								fHistList->Add(fH1Track_eta         );
								fHistList->Add(fH1MCTrack_pt        );
								fHistList->Add(fH1MCTrack_phi       );
								fHistList->Add(fH1MCTrack_eta       );
								fHistList->Add(fH1MCPrimTrack_pt    );
								fHistList->Add(fH1MCPrimTrack_phi   );
								fHistList->Add(fH1MCPrimTrack_eta   );
								fHistList->Add(fH1Jet_pt            );
								fHistList->Add(fH1Jet_phi           );
								fHistList->Add(fH1Jet_eta           );
								fHistList->Add(fH1leadJet_pt        );
								fHistList->Add(fH1leadJet_pt_dijet  );
								fHistList->Add(fH1subJet_pt_dijet   );
								fHistList->Add(fH1JetMC_pt          );
								fHistList->Add(fH1leadJetMC_pt      );
								fHistList->Add(fH1leadJetMC_pt_dijet);
								fHistList->Add(fH1subJetMC_pt_dijet );
								fHistList->Add(fH2JetsJet_dphi      );
								fHistList->Add(fH2JetsJet_deta      );
								fHistList->Add(fH2JetsJet_Aj        );
								fHistList->Add(fH2JetsJet_pt        );
								fHistList->Add(fH2JetsJetMC_dphi    );
								fHistList->Add(fH2JetsJetMC_deta    );
								fHistList->Add(fH2JetsJetMC_Aj      );
								fHistList->Add(fH2JetsJetMC_pt      );
								fHistList->Add(fH2Mult_Mtrack       );
								fHistList->Add(fH2Mult_Mlead        ); 
								fHistList->Add(fH2Mult_Mjet         );  
								fHistList->Add(fH2Mult_Njet         );  
								fHistList->Add(fH2Mult_Aj           );    
								fHistList->Add(fH2Mlead_Aj          );   
								fHistList->Add(fH2Jet_pt_Mlead      );   
								fHistList->Add(fH2Jet_pt_Munder     );  
								fHistList->Add(fH2leadJetMCptResolution );
								fHistList->Add(fH2TrackMCptResolution   );
								fHistList->Add(fH2TrackMCptEfficiency   );
								fHistList->Add(fH2AjCorrelation_MCRec   );
								fHistList->Add(fH2MleadCorrelation_MCRec);
								for(int j=0;j<5;j++){
												fHistList->Add(fH1ndiJ_ediv                        [j]);
												fHistList->Add(fH1Aj                               [j]);
												fHistList->Add(fH1Mlead                            [j]);
												fHistList->Add(fH1leadJetMC_dphiResolution         [j]);
												fHistList->Add(fH1subJetMC_dphiResolution          [j]);
												fHistList->Add(fH1leadJetMC_Efficiency             [j]);
												fHistList->Add(fH1subJetMC_Efficiency              [j]);
												for(int k=0;k<5;k++){
																fHistList->Add(fH1JetHadron_dphi_ediv               [j][k]);
																fHistList->Add(fH1JetHadron_dphi_tptweight_ediv     [j][k]);
																fHistList->Add(fH1JetHadron_dphi_tJptweight_ediv    [j][k]);
																fHistList->Add(fH1JetHadronMC_dphi_ediv               [j][k]);
																fHistList->Add(fH1JetHadronMC_dphi_tptweight_ediv     [j][k]);
																fHistList->Add(fH1JetHadronMC_dphi_tJptweight_ediv    [j][k]);
																fHistList->Add(fH1JetHadronMCPrim_dphi_ediv               [j][k]);
																fHistList->Add(fH1JetHadronMCPrim_dphi_tptweight_ediv     [j][k]);
																fHistList->Add(fH1JetHadronMCPrim_dphi_tJptweight_ediv    [j][k]);
												}
								}
								for(int j=0;j<3;j++){
												fHistList->Add(fH1ndiJ_2040Mlead    [j]);
												fHistList->Add(fH1ndiJ_2040Aj       [j]);
												for(int k=0;k<5;k++){
																fHistList->Add(fH1JetHadron_dphi_tptweight2040_Mleaddep     [j][k]);
																fHistList->Add(fH1JetHadron_dphi_tptweight2040_Ajdep        [j][k]);
																fHistList->Add(fH1JetHadronMC_dphi_tptweight2040_Mleaddep     [j][k]);
																fHistList->Add(fH1JetHadronMC_dphi_tptweight2040_Ajdep        [j][k]);
																fHistList->Add(fH1JetHadronMCPrim_dphi_tptweight2040_Mleaddep     [j][k]);
																fHistList->Add(fH1JetHadronMCPrim_dphi_tptweight2040_Ajdep        [j][k]);
												}
								}
				}
				else{
								fHistList->Add(fH1Events            );
								fHistList->Add(fH1Track_pt          );
								fHistList->Add(fH1Track_phi         );
								fHistList->Add(fH1Track_eta         );
								fHistList->Add(fH1Jet_pt            );
								fHistList->Add(fH1Jet_phi           );
								fHistList->Add(fH1Jet_eta           );
								fHistList->Add(fH1leadJet_pt        );
								fHistList->Add(fH1leadJet_pt_dijet  );
								fHistList->Add(fH1subJet_pt_dijet   );
								fHistList->Add(fH2JetsJet_dphi      );
								fHistList->Add(fH2JetsJet_deta      );
								fHistList->Add(fH2JetsJet_Aj        );
								fHistList->Add(fH2JetsJet_pt        );
								fHistList->Add(fH2Mult_Mtrack       );
								fHistList->Add(fH2Mult_Mlead        ); 
								fHistList->Add(fH2Mult_Mjet         );  
								fHistList->Add(fH2Mult_Njet         );  
								fHistList->Add(fH2Mult_Aj           );    
								fHistList->Add(fH2Mlead_Aj          );   
								fHistList->Add(fH2Jet_pt_Mlead      );   
								fHistList->Add(fH2Jet_pt_Munder     );  
								for(int j=0;j<5;j++){
												fHistList->Add(fH1ndiJ_ediv                        [j]);
												fHistList->Add(fH1Aj                               [j]);
												fHistList->Add(fH1Mlead                            [j]);
												for(int k=0;k<5;k++){
																fHistList->Add(fH1JetHadron_dphi_ediv               [j][k]);
																fHistList->Add(fH1JetHadron_dphi_tptweight_ediv     [j][k]);
																fHistList->Add(fH1JetHadron_dphi_tJptweight_ediv    [j][k]);
												}
								}
								for(int j=0;j<3;j++){
												fHistList->Add(fH1ndiJ_2040Mlead    [j]);
												fHistList->Add(fH1ndiJ_2040Aj       [j]);
												for(int k=0;k<5;k++){
																fHistList->Add(fH1JetHadron_dphi_tptweight2040_Mleaddep     [j][k]);
																fHistList->Add(fH1JetHadron_dphi_tptweight2040_Ajdep        [j][k]);
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
												fH1Xsec  ->Fill(0.,fxsec);
												fH1Trials->Fill(0.,ftrial);
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

				AliAODEvent *fAOD;
				TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
				if( handler && handler->InheritsFrom("AliAODInputHandler") ) {
								fAOD  =  ((AliAODInputHandler*)handler)->GetEvent();
								if(fUseAODInput) fAODIn = fAOD;
								if (fDebug > 1)  Printf("%s:%d AOD event from input", (char*)__FILE__,__LINE__);
				}
				else {
								handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
								if( handler && handler->InheritsFrom("AliAODHandler") ) {
												fAOD = ((AliAODHandler*)handler)->GetAOD();
												fAODIn = fAOD;
												if (fDebug > 1)  Printf("%s:%d AOD event from output", (char*)__FILE__,__LINE__);
								}
				}

				if(!fAODIn && !fUseAODInput){ // case we have AOD in input & output and want jets from output
								TObject* outHandler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
								if( outHandler && outHandler->InheritsFrom("AliAODHandler") ) {
												fAODIn = ((AliAODHandler*)outHandler)->GetAOD();
												if (fDebug > 1)  Printf("%s:%d jets from output AOD", (char*)__FILE__,__LINE__);
								}
				}
				if (!fAODIn) {
								Printf("ERROR %s : fAODIn not available",(char*)__FILE__);
								return;
				}

				AliInputEventHandler* inputHandler = (AliInputEventHandler*)
								((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
				if(!(inputHandler->IsEventSelected() & AliVEvent::kMB)){
								if (fDebug > 1 ) Printf(" Trigger Selection: event REJECTED ... ");
								return;
				}
				fH1Events->Fill(0);

				AliAODHeader* aliH = dynamic_cast <AliAODHeader*> (fAODIn->GetHeader());
				if(!aliH){
								Printf("ERROR: AliAODHeader not available");
								return;
				}
				double Mult = aliH->GetRefMultiplicity();

				// start jet analysis  -------------------------Init.
				Float_t pi=TMath::Pi();

				Double_t Jet_pt   [20][5000];
				Double_t Jet_phi  [20][5000];
				Double_t Jet_eta  [20][5000];
				Double_t Jet_area [20][5000];
				Double_t subJet_pt [20][5000];
				Double_t subJet_eta[20][5000];
				Double_t Track_n  ;
				Double_t Track_pt [5000];
				Double_t Track_eta[5000];
				Double_t Track_phi[5000];
				Double_t MCTrack_n  ;
				Double_t MCTrack_pt [5000];
				Double_t MCTrack_eta[5000];
				Double_t MCTrack_phi[5000];

				Track_n=0;MCTrack_n=0;
				for(int i=0;i<20;i++){
								for(int j=0;j<1000;j++){
												Jet_pt[i][j]=0.;
												Jet_phi[i][j]=999.;
												Jet_eta[i][j]=999.;
												Jet_area[i][j]=999.;
												subJet_pt[i][j]=0.;
												subJet_eta[i][j]=999.;
												Track_pt [j]=0.;
												Track_phi[j]=999.;
												Track_eta[j]=999.;
												MCTrack_pt [j]=0.;
												MCTrack_phi[j]=999.;
												MCTrack_eta[j]=999.;
								}
				}

				int nLJetAOD=999; double ptLJetAOD=0;double phiLJetAOD=999.;double etaLJetAOD=999.;double ptsLJetAOD=0;double phisLJetAOD=900.;double etasLJetAOD=900.;
				int nLJetMC2=999; double ptLJetMC2=0;double phiLJetMC2=999.;double etaLJetMC2=999.;double ptsLJetMC2=0;double phisLJetMC2=900.;double etasLJetMC2=900.;
				int nLJetMC =999; double ptLJetMC =0;double phiLJetMC =999.;double etaLJetMC =999.;double ptsLJetMC =0;double phisLJetMC =900.;double etasLJetMC =900.;
				bool findLJetAOD=false;
				bool findLJetMC2=false;
				bool findDiJet=false,findDiJetMC=false;
				int nLJet = 999;
				int Mjet_tot =0;
				int Njet_tot =0;

				double Aj=99.,AjMC=99.;
				double Mlead=99.,MleadMC=99.;
				int    Munder=99.;

				////--------------------------------------------------------------------Init.

				TString cAdd = "";
				TString Branchname_gen,Branchname_gen2,Branchname_rec;
				if((JFAlg=="ANTIKT")||(JFAlg=="KT")){
								cAdd += Form("%02d_",(int)((Radius+0.01)*10.));
								cAdd += Form("B%d",(int)BackM);
								cAdd += Form("_Filter%05d",Filtermask);
								cAdd += Form("_Cut%05d",(int)(1000.*TrackPtcut));
								cAdd += Form("_Skip%02d",SkipCone);
								Branchname_gen  = Form("clustersAODMC_%s%s",JFAlg.Data(),cAdd.Data());
								Branchname_gen2 = Form("clustersAODMC2_%s%s",JFAlg.Data(),cAdd.Data());
								Branchname_rec  = Form("clustersAOD_%s%s",JFAlg.Data(),cAdd.Data());
				}
				else{
								cAdd += Form("%02d_",(int)((Radius+0.01)*10.));
								cAdd += Form("B%d",(int)BackM);
								cAdd += Form("_Filter%05d",Filtermask);
								cAdd += Form("_Cut%05d",(int)(1000.*TrackPtcut));
								Branchname_gen  = Form("jetsAODMC_%s%s",JFAlg.Data(),cAdd.Data());
								Branchname_gen2 = Form("jetsAODMC2_%s%s",JFAlg.Data(),cAdd.Data());
								Branchname_rec  = Form("jetsAOD_%s%s",JFAlg.Data(),cAdd.Data());
				}



				//count number of tracks@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
				//Reconstructed Track
				TClonesArray* tracks = dynamic_cast <TClonesArray*> (fAODIn->GetTracks());
				if(!tracks){
								if (fDebug > 1)  Printf("%s:%d could not get AODtracks", (char*)__FILE__,__LINE__);
								return;
				}

				Bool_t TrackEff[5000];
				for(int i=0;i<5000;i++){
								TrackEff[i]=false;
				}
				Int_t nt = fAODIn->GetNumberOfTracks();
				AliAODTrack* trackAOD=NULL;
				for(int ntrack =0;ntrack<nt;ntrack++){
								trackAOD = (AliAODTrack*) (tracks->At(ntrack));
								Bool_t bgoodT=false;
								if(Filtermask!=272){if(trackAOD->TestFilterMask(Filtermask))bgoodT=true;}
								else               {if(trackAOD->IsHybridGlobalConstrainedGlobal())bgoodT=true;} //for hybrid Track cuts
								if(!bgoodT)continue;
								if(TMath::Abs(trackAOD->Eta())<0.9){
												Track_n++;
												fH1Track_pt ->Fill(trackAOD->Pt()*TrackEScale);
												fH1Track_phi->Fill(trackAOD->Phi());
												fH1Track_eta->Fill(trackAOD->Eta());
												//cout<<"Scale "<<TrackEScale<<"  org pt "<<trackAOD->Pt()<< " scaled pt "<< trackAOD->Pt()*TrackEScale <<endl;
												if(IsMC){
																// track pt resplution-------------------
																Int_t MCID = TMath::Abs(trackAOD->GetLabel());
																TClonesArray* mctracks = dynamic_cast <TClonesArray*> (fAODIn->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
																if(!mctracks){
																				if (fDebug > 1)  Printf("%s:%d could not get AODMCtracks", (char*)__FILE__,__LINE__);
																				continue;
																}
																AliAODMCParticle *trackMCAOD = (AliAODMCParticle*) mctracks->At(MCID);
																fH2TrackMCptResolution->Fill(trackMCAOD->Pt(),trackAOD->Pt());
																TrackEff[MCID]=true;
																// --------------------------------------
												}
								}
				}
				if(IsMC){
								//MC Track
								TClonesArray* mctracks = dynamic_cast <TClonesArray*> (fAODIn->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
								if(!mctracks){
												if (fDebug > 1)  Printf("%s:%d could not get AODMCtracks", (char*)__FILE__,__LINE__);
												return;
								}
								Int_t ntmc = mctracks->GetEntriesFast();
								AliAODMCParticle* trackMCAOD;
								int lastprim=0;
								for(int ntrack =0;ntrack<ntmc;ntrack++){
												trackMCAOD = (AliAODMCParticle*) (mctracks->At(ntrack));
												if((trackMCAOD->IsPhysicalPrimary())==1)lastprim=ntrack;
								}
								for(int ntrack =0;ntrack<ntmc;ntrack++){
												trackMCAOD = (AliAODMCParticle*) (mctracks->At(ntrack));
												if((trackMCAOD->GetPdgCode()>10)&&((trackMCAOD->GetMother())>1)&&(ntrack>lastprim)&&(trackMCAOD->Charge())){// for Decay particles
																if(TMath::Abs(trackMCAOD->Eta())<0.9){
																				fH1MCTrack_pt ->Fill(trackMCAOD->Pt());
																				fH1MCTrack_phi->Fill(trackMCAOD->Phi());
																				fH1MCTrack_eta->Fill(trackMCAOD->Eta());
																				if(TrackEff[ntrack])fH2TrackMCptEfficiency->Fill(trackMCAOD->Pt(),1);
																				else                fH2TrackMCptEfficiency->Fill(trackMCAOD->Pt(),0);
																}
												}
												if((trackMCAOD->IsPhysicalPrimary())&&(trackMCAOD->Charge())){// for Physical particles
																if(TMath::Abs(trackMCAOD->Eta())<0.9){
																				MCTrack_n++;
																				fH1MCTrack_pt ->Fill(trackMCAOD->Pt());
																				fH1MCTrack_phi->Fill(trackMCAOD->Phi());
																				fH1MCTrack_eta->Fill(trackMCAOD->Eta());
																				fH1MCPrimTrack_pt ->Fill(trackMCAOD->Pt());
																				fH1MCPrimTrack_phi->Fill(trackMCAOD->Phi());
																				fH1MCPrimTrack_eta->Fill(trackMCAOD->Eta());
																				if(TrackEff[ntrack])fH2TrackMCptEfficiency->Fill(trackMCAOD->Pt(),1);
																				else                fH2TrackMCptEfficiency->Fill(trackMCAOD->Pt(),0);
																}
												}
								}
				}
				//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  count number of tracks




				for(int algorithm=0;algorithm<3;algorithm++){

								if(algorithm==0)fJetBranch   = Branchname_rec.Data();
								if(algorithm==1)fJetBranch   = Branchname_gen2.Data();
								if(algorithm==2)fJetBranch   = Branchname_gen.Data();

								if((!IsMC&&(algorithm==1||algorithm==2)))continue;

								TClonesArray* jets = dynamic_cast <TClonesArray*> (fAODIn->FindListObject(fJetBranch.Data()));
								if(!jets){
												printf(" Tere are no Branch named %s \n",fJetBranch.Data());
												continue;
								}
								Int_t nj = jets->GetEntriesFast();
								if (fDebug) printf("There are %5d jets in the event \n", nj);
								AliAODJet* jetsAOD;
								//Find Leading Jet -------------------------------------------------------
								for(int njet =0;njet<nj;njet++){
												jetsAOD = (AliAODJet*) (jets->At(njet));
												Jet_pt   [algorithm][njet] = jetsAOD->Pt()*JetEScale;
												Jet_phi  [algorithm][njet] = jetsAOD->Phi();  
												Jet_eta  [algorithm][njet] = jetsAOD->Eta();
												Jet_area [algorithm][njet] = jetsAOD->EffectiveAreaCharged();


												TRefArray *reftracks = jetsAOD->GetRefTracks();
												if(algorithm==0){if(Jet_pt[algorithm][njet]>1.)Mjet_tot +=  reftracks->GetEntriesFast();Njet_tot++;}
												double eta_cut_Jet=0.5;
												if((TMath::Abs(Jet_eta[algorithm][njet])<eta_cut_Jet)&&(Jet_pt[algorithm][njet]>10.)){
																if(algorithm==0){
																				fH1Jet_pt ->Fill(Jet_pt [algorithm][njet]);  
																				fH1Jet_phi->Fill(Jet_phi[algorithm][njet]);  
																				fH1Jet_eta->Fill(Jet_eta[algorithm][njet]);  
																				if(Jet_pt[algorithm][njet]>ptLJetAOD){
																								findLJetAOD=true;
																								nLJetAOD=njet;ptLJetAOD=Jet_pt[algorithm][njet];phiLJetAOD=Jet_phi[algorithm][njet];etaLJetAOD=Jet_eta[algorithm][njet];
																				}
																}
																if(algorithm==1){
																				fH1JetMC_pt->Fill(Jet_pt[algorithm][njet]);  
																				if(Jet_pt[algorithm][njet]>ptLJetMC2){
																								findLJetMC2=true;
																								nLJetMC2=njet;ptLJetMC2=Jet_pt[algorithm][njet];phiLJetMC2=Jet_phi[algorithm][njet];etaLJetMC2=Jet_eta[algorithm][njet];
																				}
																}
																if(algorithm==2){
																				if(Jet_pt[algorithm][njet]>ptLJetMC){
																								nLJetMC=njet;ptLJetMC=Jet_pt[algorithm][njet];phiLJetMC=Jet_phi[algorithm][njet];etaLJetMC=Jet_eta[algorithm][njet];
																				}
																}
												}
								}//njet loop
								if(algorithm==0){nLJet=nLJetAOD;fH1leadJet_pt  ->Fill(Jet_pt[algorithm][nLJet]);}
								if(algorithm==1){nLJet=nLJetMC2;fH1leadJetMC_pt->Fill(Jet_pt[algorithm][nLJet]);}
								if(algorithm==2){nLJet=nLJetMC;}
								if(findLJetAOD&&(algorithm==0)){
												jetsAOD = (AliAODJet*) (jets->At(nLJet));
												TRefArray *reftracks = jetsAOD->GetRefTracks();
												Mlead = reftracks->GetEntriesFast();
								}
								if(findLJetMC2&&(algorithm==1)){
												jetsAOD = (AliAODJet*) (jets->At(nLJetMC2));
												TRefArray *reftracks = jetsAOD->GetRefTracks();
												MleadMC = reftracks->GetEntriesFast();
								}
								//----------------------------------------------------------- Leading Jet
								if(nj<2)continue;
								//Find Sub leading Jet ==================================================
								for(int njet=0;njet<nj;njet++){
												if(njet==nLJet)continue;
												jetsAOD = (AliAODJet *)jets->At(njet);
												subJet_pt [algorithm][njet] = jetsAOD->Pt()*JetEScale;
												subJet_eta[algorithm][njet] = jetsAOD->Eta();
												double eta_cut_Jet=0.5;
												if((TMath::Abs(subJet_eta[algorithm][njet])<eta_cut_Jet) && (subJet_pt[algorithm][njet]>10.)){
																if(subJet_pt[algorithm][njet]>ptsLJetAOD&&algorithm==0){
																				ptsLJetAOD=Jet_pt[algorithm][njet];phisLJetAOD=Jet_phi[algorithm][njet];etasLJetAOD=Jet_eta[algorithm][njet];
																}
																if(subJet_pt[algorithm][njet]>ptsLJetMC2 &&algorithm==1){
																				ptsLJetMC2=Jet_pt[algorithm][njet];phisLJetMC2=Jet_phi[algorithm][njet];etasLJetMC2=Jet_eta[algorithm][njet];
																}
																if(subJet_pt[algorithm][njet]>ptsLJetMC &&algorithm==2){
																				ptsLJetMC =Jet_pt[algorithm][njet];phisLJetMC =Jet_phi[algorithm][njet];etasLJetMC =Jet_eta[algorithm][njet];
																}
												}
								}
								//====================================================== Sub leading Jet 

								double Leading_pt=0.;double Leading_phi=999.;double Leading_eta=999.;double sLeading_pt=0.;double sLeading_phi=999.;double sLeading_eta=999.;
								if(algorithm==0){Leading_pt=ptLJetAOD;Leading_phi=phiLJetAOD;Leading_eta=etaLJetAOD;sLeading_pt=ptsLJetAOD;sLeading_phi=phisLJetAOD;sLeading_eta=etasLJetAOD;}
								if(algorithm==1){Leading_pt=ptLJetMC2;Leading_phi=phiLJetMC2;Leading_eta=etaLJetMC2;sLeading_pt=ptsLJetMC2;sLeading_phi=phisLJetMC2;sLeading_eta=etasLJetMC2;}
								if(algorithm==2){Leading_pt=ptLJetMC ;Leading_phi=phiLJetMC ;Leading_eta=etaLJetMC ;sLeading_pt=ptsLJetMC ;sLeading_phi=phisLJetMC ;sLeading_eta=etasLJetMC ;}

								////Di-Jet event trigger +++++++++++++++++++++++++++++++++++++++++++++++
								double DPhi = DeltaPhi(Leading_phi,sLeading_phi);
								double DEta = Leading_eta-sLeading_eta;
								if(algorithm==0){
												fH2JetsJet_dphi->Fill(Leading_pt,DPhi);
												fH2JetsJet_deta->Fill(Leading_pt,DEta);
								}
								if(algorithm==1){
												fH2JetsJetMC_dphi->Fill(Leading_pt,DPhi);
												fH2JetsJetMC_deta->Fill(Leading_pt,DEta);
								}
								if((TMath::Cos(DPhi)<-0.5)&&(Leading_pt>10.)&&(sLeading_pt>10.)){
												if(algorithm==0)Aj   = (Leading_pt-sLeading_pt)/(Leading_pt+sLeading_pt);
												if(algorithm==1)AjMC = (Leading_pt-sLeading_pt)/(Leading_pt+sLeading_pt);
												if(algorithm==0){
																fH1subJet_pt_dijet ->Fill(sLeading_pt);
																fH1leadJet_pt_dijet->Fill(Leading_pt);
																fH2JetsJet_Aj      ->Fill(Leading_pt,Aj);
																fH2JetsJet_pt      ->Fill(Leading_pt,sLeading_pt);
																fH2Mult_Aj         ->Fill(Mult,Aj); 
																fH2Mlead_Aj        ->Fill(Mlead,Aj);
																for(int eb=0;eb<5;eb++){
																				if(TMath::Abs(Leading_pt -10.-20.*(eb))<10.){
																								fH1Aj[eb]   ->Fill(Aj);
																				}
																}
												}
												if(algorithm==1){
																fH1leadJetMC_pt_dijet->Fill(Leading_pt);
																fH1subJetMC_pt_dijet ->Fill(sLeading_pt);
																fH2JetsJetMC_Aj      ->Fill(Leading_pt,AjMC);
																fH2JetsJetMC_pt      ->Fill(Leading_pt,sLeading_pt);
																findDiJetMC=true;
												}
												findDiJet=true;

								}
								////+++++++++++++++++++++++++++++++++++++++++++++++ Di-Jet event trigger 

								if(algorithm!=0)continue;// for only data & reconstructed Jets


								//Jet-Hadron Correlation###############################################
								if((findDiJet)&&(Leading_pt>10.)&&(sLeading_pt>10.)){
												double eta_cut_Jet=0.5;
												if(TMath::Abs(Leading_eta)<eta_cut_Jet){
																for(int eb=0;eb<5;eb++){
																				if(TMath::Abs(Leading_pt -10.-20.*(eb))<10.){
																								fH1ndiJ_ediv[eb]->Fill(1);
																								if(eb==1){
																												if((0<Mlead)&&Mlead<7)              {fH1ndiJ_2040Mlead[0]->Fill(1);}
																												else if((7<=Mlead)&&(Mlead<10))     {fH1ndiJ_2040Mlead[1]->Fill(1);}
																												else                                {fH1ndiJ_2040Mlead[2]->Fill(1);}
																												if((0<Aj)&&(Aj<0.19))               {fH1ndiJ_2040Aj   [0]->Fill(1);}
																												else if((0.19<=Aj)&&(Aj<0.38))      {fH1ndiJ_2040Aj   [1]->Fill(1);}
																												else                                {fH1ndiJ_2040Aj   [2]->Fill(1);}
																								}
																								fH1Mlead[eb]->Fill(Mlead);
																								for(int ntrack =0;ntrack<nt;ntrack++){
																												trackAOD = (AliAODTrack*) (fAODIn->GetTrack(ntrack));
																												Bool_t bgoodT=false;
																												if(Filtermask!=272){if(trackAOD->TestFilterMask(Filtermask))bgoodT=true;}
																												else{               if(trackAOD->IsHybridGlobalConstrainedGlobal())bgoodT=true;} //for hybrid Track cuts
																												if(!bgoodT)continue;
																												Track_pt   [ntrack]      = (trackAOD->Pt()*TrackEScale);
																												Track_phi  [ntrack]      = trackAOD->Phi();
																												Track_eta  [ntrack]      = trackAOD->Eta();

																												//cout<<"Scale "<<TrackEScale<<"  org pt "<<trackAOD->Pt()<< " scaled pt "<< trackAOD->Pt()*TrackEScale <<endl;

																												double DelPhi = DeltaPhi(Leading_phi,Track_phi[ntrack]);
																												if(TMath::Abs(Track_eta[ntrack])<0.9){
																																if((TMath::Abs(DelPhi-pi/2.)<pi/8.)||((DelPhi+pi/2.)<pi/8.)||(TMath::Abs(DelPhi-3./2.*pi)<pi/8.))Munder++;
																																for(int teb=0;teb<5;teb++){
																																				if(teb==0){if(!( Track_pt[ntrack]>0.15))continue;}
																																				if(teb==1){if(!((Track_pt[ntrack]<1.5)&&(Track_pt[ntrack]>0.15)))continue;}
																																				if(teb==2){if(!((Track_pt[ntrack]<3.0)&&(Track_pt[ntrack]>1.5)))continue;}
																																				if(teb==3){if(!((Track_pt[ntrack]<4.5)&&(Track_pt[ntrack]>3.0)))continue;}
																																				if(teb==4){if(!( Track_pt[ntrack]>4.5))continue;}
																																				fH1JetHadron_dphi_ediv                [eb][teb]->Fill(DelPhi); 
																																				fH1JetHadron_dphi_tptweight_ediv      [eb][teb]->Fill(DelPhi,Track_pt[ntrack]);
																																				fH1JetHadron_dphi_tJptweight_ediv     [eb][teb]->Fill(DelPhi,Track_pt[ntrack]/Leading_pt);
																																				if(eb==1){
																																								if((0<Mlead)&&Mlead<7)         {fH1JetHadron_dphi_tptweight2040_Mleaddep[0][teb]->Fill(DelPhi,Track_pt[ntrack]);}
																																								else if((7<=Mlead)&&(Mlead<10)){fH1JetHadron_dphi_tptweight2040_Mleaddep[1][teb]->Fill(DelPhi,Track_pt[ntrack]);}
																																								else                           {fH1JetHadron_dphi_tptweight2040_Mleaddep[2][teb]->Fill(DelPhi,Track_pt[ntrack]);}
																																								if((0<Aj)&&(Aj<0.19))          {fH1JetHadron_dphi_tptweight2040_Ajdep   [0][teb]->Fill(DelPhi,Track_pt[ntrack]);}
																																								else if((0.19<=Aj)&&(Aj<0.38)) {fH1JetHadron_dphi_tptweight2040_Ajdep   [1][teb]->Fill(DelPhi,Track_pt[ntrack]);}
																																								else                           {fH1JetHadron_dphi_tptweight2040_Ajdep   [2][teb]->Fill(DelPhi,Track_pt[ntrack]);}
																																				}
																																}
																												}
																								}//Track Loop
																								if(IsMC){
																												//MC Track
																												TClonesArray* mctracks = dynamic_cast <TClonesArray*> (fAODIn->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
																												if(!mctracks){
																																if (fDebug > 1)  Printf("%s:%d could not get AODMCtracks", (char*)__FILE__,__LINE__);
																																continue;
																												}
																												Int_t ntmc = mctracks->GetEntriesFast();
																												AliAODMCParticle* trackMCAOD;
																												int lastprim=0;
																												for(int ntrack =0;ntrack<ntmc;ntrack++){
																																trackMCAOD = (AliAODMCParticle*) (mctracks->At(ntrack));
																																if((trackMCAOD->IsPhysicalPrimary())==1)lastprim=ntrack;
																												}
																												for(int ntrack =0;ntrack<ntmc;ntrack++){
																																trackMCAOD = (AliAODMCParticle*) (mctracks->At(ntrack));
																																if((trackMCAOD->GetPdgCode()>10)&&((trackMCAOD->GetMother())>1)&&(ntrack>lastprim)&&(trackMCAOD->Charge())){// for Decay particles
																																				MCTrack_pt   [ntrack]      = trackMCAOD->Pt();
																																				MCTrack_phi  [ntrack]      = trackMCAOD->Phi();
																																				MCTrack_eta  [ntrack]      = trackMCAOD->Eta();
																																				double DelPhi = DeltaPhi(Leading_phi,MCTrack_phi[ntrack]);
																																				if(TMath::Abs(MCTrack_eta[ntrack])<0.9){
																																								for(int teb=0;teb<5;teb++){
																																												if(teb==0){if(!( MCTrack_pt[ntrack]>0.15))continue;}
																																												if(teb==1){if(!((MCTrack_pt[ntrack]<1.5)&&(MCTrack_pt[ntrack]>0.15)))continue;}
																																												if(teb==2){if(!((MCTrack_pt[ntrack]<3.0)&&(MCTrack_pt[ntrack]>1.5)))continue;}
																																												if(teb==3){if(!((MCTrack_pt[ntrack]<4.5)&&(MCTrack_pt[ntrack]>3.0)))continue;}
																																												if(teb==4){if(!( MCTrack_pt[ntrack]>4.5))continue;}
																																												fH1JetHadronMC_dphi_ediv                [eb][teb]->Fill(DelPhi); 
																																												fH1JetHadronMC_dphi_tptweight_ediv      [eb][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);
																																												fH1JetHadronMC_dphi_tJptweight_ediv     [eb][teb]->Fill(DelPhi,MCTrack_pt[ntrack]/Leading_pt);
																																												if(eb==1){
																																																if((0<Mlead)&&Mlead<7)         {fH1JetHadronMC_dphi_tptweight2040_Mleaddep[0][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else if((7<=Mlead)&&(Mlead<10)){fH1JetHadronMC_dphi_tptweight2040_Mleaddep[1][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else                           {fH1JetHadronMC_dphi_tptweight2040_Mleaddep[2][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																if((0<Aj)&&(Aj<0.19))          {fH1JetHadronMC_dphi_tptweight2040_Ajdep   [0][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else if((0.19<=Aj)&&(Aj<0.38)) {fH1JetHadronMC_dphi_tptweight2040_Ajdep   [1][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else                           {fH1JetHadronMC_dphi_tptweight2040_Ajdep   [2][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																												}
																																								}
																																				}
																																}
																																if((trackMCAOD->IsPhysicalPrimary())&&(trackMCAOD->Charge())){// for Physical particles
																																				MCTrack_pt   [ntrack]      = trackMCAOD->Pt();
																																				MCTrack_phi  [ntrack]      = trackMCAOD->Phi();
																																				MCTrack_eta  [ntrack]      = trackMCAOD->Eta();
																																				double DelPhi = DeltaPhi(Leading_phi,MCTrack_phi[ntrack]);
																																				if(TMath::Abs(MCTrack_eta[ntrack])<0.9){
																																								for(int teb=0;teb<5;teb++){
																																												if(teb==0){if(!( MCTrack_pt[ntrack]>0.15))continue;}
																																												if(teb==1){if(!((MCTrack_pt[ntrack]<1.5)&&(MCTrack_pt[ntrack]>0.15)))continue;}
																																												if(teb==2){if(!((MCTrack_pt[ntrack]<3.0)&&(MCTrack_pt[ntrack]>1.5)))continue;}
																																												if(teb==3){if(!((MCTrack_pt[ntrack]<4.5)&&(MCTrack_pt[ntrack]>3.0)))continue;}
																																												if(teb==4){if(!( MCTrack_pt[ntrack]>4.5))continue;}
																																												fH1JetHadronMC_dphi_ediv                [eb][teb]->Fill(DelPhi); 
																																												fH1JetHadronMC_dphi_tptweight_ediv      [eb][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);
																																												fH1JetHadronMC_dphi_tJptweight_ediv     [eb][teb]->Fill(DelPhi,MCTrack_pt[ntrack]/Leading_pt);
																																												fH1JetHadronMCPrim_dphi_ediv                [eb][teb]->Fill(DelPhi); 
																																												fH1JetHadronMCPrim_dphi_tptweight_ediv      [eb][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);
																																												fH1JetHadronMCPrim_dphi_tJptweight_ediv     [eb][teb]->Fill(DelPhi,MCTrack_pt[ntrack]/Leading_pt);
																																												if(eb==1){
																																																if((0<Mlead)&&Mlead<7)         {fH1JetHadronMC_dphi_tptweight2040_Mleaddep[0][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else if((7<=Mlead)&&(Mlead<10)){fH1JetHadronMC_dphi_tptweight2040_Mleaddep[1][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else                           {fH1JetHadronMC_dphi_tptweight2040_Mleaddep[2][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																if((0<Aj)&&(Aj<0.19))          {fH1JetHadronMC_dphi_tptweight2040_Ajdep   [0][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else if((0.19<=Aj)&&(Aj<0.38)) {fH1JetHadronMC_dphi_tptweight2040_Ajdep   [1][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else                           {fH1JetHadronMC_dphi_tptweight2040_Ajdep   [2][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}

																																																if((0<Mlead)&&Mlead<7)         {fH1JetHadronMCPrim_dphi_tptweight2040_Mleaddep[0][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else if((7<=Mlead)&&(Mlead<10)){fH1JetHadronMCPrim_dphi_tptweight2040_Mleaddep[1][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else                           {fH1JetHadronMCPrim_dphi_tptweight2040_Mleaddep[2][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																if((0<Aj)&&(Aj<0.19))          {fH1JetHadronMCPrim_dphi_tptweight2040_Ajdep   [0][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else if((0.19<=Aj)&&(Aj<0.38)) {fH1JetHadronMCPrim_dphi_tptweight2040_Ajdep   [1][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																																else                           {fH1JetHadronMCPrim_dphi_tptweight2040_Ajdep   [2][teb]->Fill(DelPhi,MCTrack_pt[ntrack]);}
																																												}
																																								}
																																				}
																																}
																												}
																								}
																				}
																}// Momentum Loop Jet
																fH2Jet_pt_Munder   ->Fill(Leading_pt,(double)Munder/(1.8*pi/2.)*Jet_area[0][nLJet]);
																fH2Jet_pt_Mlead    ->Fill(Leading_pt,Mlead);
												}//eta cut
								}// Di-Jet
								//############################################### Jet-Hadron Correlation
				}// algorithm LOOP
				if(IsMC){
								for(int eb=0;eb<5;eb++){
												double DPhi,DEta;
												if(TMath::Abs(ptLJetAOD -10.-20.*(eb))<10.){
																DPhi = DeltaPhi(phiLJetMC,phiLJetAOD);
																DEta = etaLJetMC-etaLJetAOD;
																fH1leadJetMC_dphiResolution[eb]->Fill(DPhi);
																if(sqrt(pow(DPhi,2)+pow(DEta,2))<0.4)fH1leadJetMC_Efficiency[eb]->Fill(1);
																else                                 fH1leadJetMC_Efficiency[eb]->Fill(0);
																DPhi = DeltaPhi(phisLJetMC,phisLJetAOD);
																DEta = etasLJetMC-etasLJetAOD;
																fH1subJetMC_dphiResolution[eb]->Fill(DPhi);
																if(sqrt(pow(DPhi,2)+pow(DEta,2))<0.4)fH1subJetMC_Efficiency[eb]->Fill(1);
																else                                 fH1subJetMC_Efficiency[eb]->Fill(0);
																DPhi = DeltaPhi(phiLJetMC2,phiLJetAOD);
																DEta = etaLJetMC2-etaLJetAOD;

																if(sqrt(pow(DPhi,2)+pow(DEta,2))<0.4)fH2leadJetMCptResolution->Fill(ptLJetMC2,ptLJetAOD);
																if(findDiJetMC)fH2AjCorrelation_MCRec   ->Fill(AjMC,Aj);
																if(findDiJetMC)fH2MleadCorrelation_MCRec->Fill(MleadMC,Mlead);
												}
								}
								fH2Mult_Mtrack->Fill(Mult,Track_n); 
								fH2Mult_Mjet  ->Fill(Mult,Mjet_tot);
								fH2Mult_Njet  ->Fill(Mult,Njet_tot);
								if(findLJetAOD)fH2Mult_Mlead ->Fill(Mult,Mlead);   
				}
				else{
								fH2Mult_Mtrack->Fill(Mult,Track_n); 
								fH2Mult_Mjet  ->Fill(Mult,Mjet_tot);
								fH2Mult_Njet  ->Fill(Mult,Njet_tot);
								if(findLJetAOD)fH2Mult_Mlead ->Fill(Mult,Mlead);   
				}

				PostData(1, fHistList);
				return;
}      

//________________________________________________________________________
void AliAnalysisTaskJetHadronCorrelation::Terminate(Option_t *){
				// Terminate analysis
				if (fDebug) printf("AnalysisTaskPt: Terminate() \n");
}

Double_t AliAnalysisTaskJetHadronCorrelation::DeltaPhi(Double_t phi1,Double_t phi2){
				Float_t pi=TMath::Pi();
				Double_t dphi = phi1-phi2;
				if     (dphi<(-1./2*pi))dphi = dphi +2*pi;
				else if(dphi>( 3./2*pi))dphi = dphi -2*pi;
				return dphi;
}
