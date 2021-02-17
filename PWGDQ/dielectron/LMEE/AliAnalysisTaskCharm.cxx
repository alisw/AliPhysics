/*************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//    Class to produce charm cocktail                  (analysis task)   //
//                                        (description in .h file)       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// c++ includes
#include <iostream>
// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TChain.h"
#include "TVector2.h"
#include "TGraph.h"
#include <TMath.h>
//#include <TDatabasePDG.h>
#include <TClonesArray.h>
#include "TParticle.h"
#include "TLorentzVector.h"
// AliRoot includes
#include "AliAnalysisTaskCharm.h" // this task
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliVParticle.h"
#include "AliInputEventHandler.h"

ClassImp(AliAnalysisTaskCharm)

//______________________________| Default Constructor
AliAnalysisTaskCharm::AliAnalysisTaskCharm():
AliAnalysisTaskSE(),
AliCocktailSmearing(),
fMcEvent(0x0),
fMcHandler(0x0),
fNbEvents(0),
fProcessType(0),
fPtCutHigh(1e6),
fPtCutLow(0.2),
fEtamin(-0.8),
fEtamax(0.8),
fScaleByRAA(kFALSE),
fScaleByCNM(kFALSE),
fgraphCNM(0),
fTakeptOfDCNM(kFALSE),
fApplywm(kFALSE),
fEventWeight(kFALSE),
fSelectoneccbar(kFALSE),
fSelectcleanhistory(kFALSE),
fNbEvent(-1),
fNbEventCounter(0),
hNEvents(0),
hNEventsW(0),
hQuarkMethod1(0),
hQuarkMethod2(0),
hCharm(0),
hDpm(0),
hD0(0),
hDs(0),
hDstar(0),
hLambdac(0),
hDeltaPhi_D0(0),
hDeltaPhi_D0_LHCb(0),
hDeltaPhi_D0_LHCb_e(0),
hRapElectronPt200(0),
hRapElectronPt500(0),
hPtCharmQuarkMethod1(0),
hPtCharmQuarkMethod1CNM(0),
hRapCharmQuarkMethod1(0),
hRapCharmQuarkMethod2(0),
hRapCharmQuarkMethod1Pair(0),
hRapCharmQuarkMethod2Pair(0),
hRapCharmQuark2DMethod1(0),
hRapCharmQuark2DMethod2(0),
hPhiCharmQuark1DMethod1(0),
hPhiCharmQuark1DMethod2(0),
hPtCharmQuark2DMethod1(0),
hPtCharmQuark2DMethod2(0),
hRapCquarkDmeson(0),
hPtCquarkDmeson(0),
hRapDmesonElectron(0),
hPtDmesonElectron(0),
hRapDMom2e(0),
hPtEtaElectron(0),
hPtEtaElectronE(0),
hPtEtaElectronP(0),
hEtaPositronElectron(0),
hEtaPositronElectronDelta(0),
hPhiPositronElectronDelta(0),
hPhiPositronElectronDeltaonehighDmeson(0),
hPhiPositronElectronDeltabothhighDmeson(0),
hPte_eta08(0),
hPteP_eta08(0),
hPteM_eta08(0),
hPte_y08(0),
hPteP_y08(0),
hPteM_y08(0),
hMee_ULS_simulated(0),
hMee_LS_simulated(0),
hMeePtee_ULS_eta05(0),
hMeePtee_LS_eta05(0),
hMeePtee_ULS_eta035(0),
hMeePtee_LS_eta035(0),
hPhiee_ULS_eta08_pt200(0),
hPhiee_LS_eta08_pt200(0),
hPhiee_ULS_eta08_pt400(0),
hPhiee_LS_eta08_pt400(0),
hMeePtee_ULS_eta035_phenixacc(0),
hMeePtee_LS_eta035_phenixacc(0),
hPhiee_ULS_eta08_highoneDmeson_pt200(0),
hPhiee_LS_eta08_highoneDmeson_pt200(0),
hPhiee_ULS_eta08_highbothDmeson_pt200(0),
hPhiee_LS_eta08_highbothDmeson_pt200(0),
hPhiee_ULS_eta08_highoneDmeson_pt400(0),
hPhiee_LS_eta08_highoneDmeson_pt400(0),
hPhiee_ULS_eta08_highbothDmeson_pt400(0),
hPhiee_LS_eta08_highbothDmeson_pt400(0),
hMee_ULS_eta08(0),
hMee_LS_eta08(0),
hMee_ULS_eta08_pt200(0),
hMee_LS_eta08_pt200(0),
hMeePtee_ULS_eta08(0),
hMeePtee_LS_eta08(0),
hMeePtee_ULS_eta08_pt200(0),
hMeePtee_LS_eta08_pt200(0),
hMeePtee_ULS_eta_pt(0),
hMeePtee_LS_eta_pt(0),
hMeePtee_ULS_eta08_pt300(0),
hMeePtee_LS_eta08_pt300(0),
hMeePtee_ULS_eta08_pt400(0),
hMeePtee_LS_eta08_pt400(0),
hMeePtee_ULS_eta08_pt200_opAngle50(0),
hMeePtee_LS_eta08_pt200_opAngle50(0),
hMeePtee_ULS_eta08_pt300_opAngle50(0),
hMeePtee_LS_eta08_pt300_opAngle50(0),
hMeePtee_ULS_eta08_pt400_opAngle50(0),
hMeePtee_LS_eta08_pt400_opAngle50(0),
hMeeOpAngle_ULS_eta08_pt200(0),
hMeeOpAngle_LS_eta08_pt200(0),
hMotherPt_ULS_eta08_pt200(0),
hMotherPt_LS_eta08_pt200(0),
hMotherPt_ULS_eta08_pt400(0),
hMotherPt_LS_eta08_pt400(0),
hweightcnmD(0),
hweightcnmHFE(0),
hweightcnmee(0),
fOutputList(0)
{

}

//______________________________| Specific Constructor
AliAnalysisTaskCharm::AliAnalysisTaskCharm(const Char_t* name) :
AliAnalysisTaskSE(name),
AliCocktailSmearing(),
fMcEvent(0x0),
fMcHandler(0x0),
fNbEvents(0),
fProcessType(0),
fPtCutHigh(1e6),
fPtCutLow(0.2),
fEtamin(-0.8),
fEtamax(0.8),
fScaleByRAA(kFALSE),
fScaleByCNM(kFALSE),
fgraphCNM(0),
fTakeptOfDCNM(kFALSE),
fApplywm(kFALSE),
fEventWeight(kFALSE),
fSelectoneccbar(kFALSE),
fSelectcleanhistory(kFALSE),
fNbEvent(-1),
fNbEventCounter(0),
hNEvents(0),
hNEventsW(0),
hQuarkMethod1(0),
hQuarkMethod2(0),
hCharm(0),
hDpm(0),
hD0(0),
hDs(0),
hDstar(0),
hLambdac(0),
hDeltaPhi_D0(0),
hDeltaPhi_D0_LHCb(0),
hDeltaPhi_D0_LHCb_e(0),
hRapElectronPt200(0),
hRapElectronPt500(0),
hPtCharmQuarkMethod1(0),
hPtCharmQuarkMethod1CNM(0),
hRapCharmQuarkMethod1(0),
hRapCharmQuarkMethod2(0),
hRapCharmQuarkMethod1Pair(0),
hRapCharmQuarkMethod2Pair(0),
hRapCharmQuark2DMethod1(0),
hRapCharmQuark2DMethod2(0),
hPhiCharmQuark1DMethod1(0),
hPhiCharmQuark1DMethod2(0),
hPtCharmQuark2DMethod1(0),
hPtCharmQuark2DMethod2(0),
hRapCquarkDmeson(0),
hPtCquarkDmeson(0),
hRapDmesonElectron(0),
hPtDmesonElectron(0),
hRapDMom2e(0),
hPtEtaElectron(0),
hPtEtaElectronE(0),
hPtEtaElectronP(0),
hEtaPositronElectron(0),
hEtaPositronElectronDelta(0),
hPhiPositronElectronDelta(0),
hPhiPositronElectronDeltaonehighDmeson(0),
hPhiPositronElectronDeltabothhighDmeson(0),
hPte_eta08(0),
hPteP_eta08(0),
hPteM_eta08(0),
hPte_y08(0),
hPteP_y08(0),
hPteM_y08(0),
hMee_ULS_simulated(0),
hMee_LS_simulated(0),
hMeePtee_ULS_eta05(0),
hMeePtee_LS_eta05(0),
hMeePtee_ULS_eta035(0),
hMeePtee_LS_eta035(0),
hPhiee_ULS_eta08_pt200(0),
hPhiee_LS_eta08_pt200(0),
hPhiee_ULS_eta08_pt400(0),
hPhiee_LS_eta08_pt400(0),
hMeePtee_ULS_eta035_phenixacc(0),
hMeePtee_LS_eta035_phenixacc(0),
hPhiee_ULS_eta08_highoneDmeson_pt200(0),
hPhiee_LS_eta08_highoneDmeson_pt200(0),
hPhiee_ULS_eta08_highbothDmeson_pt200(0),
hPhiee_LS_eta08_highbothDmeson_pt200(0),
hPhiee_ULS_eta08_highoneDmeson_pt400(0),
hPhiee_LS_eta08_highoneDmeson_pt400(0),
hPhiee_ULS_eta08_highbothDmeson_pt400(0),
hPhiee_LS_eta08_highbothDmeson_pt400(0),
hMee_ULS_eta08(0),
hMee_LS_eta08(0),
hMee_ULS_eta08_pt200(0),
hMee_LS_eta08_pt200(0),
hMeePtee_ULS_eta08(0),
hMeePtee_LS_eta08(0),
hMeePtee_ULS_eta08_pt200(0),
hMeePtee_LS_eta08_pt200(0),
hMeePtee_ULS_eta_pt(0),
hMeePtee_LS_eta_pt(0),
hMeePtee_ULS_eta08_pt300(0),
hMeePtee_LS_eta08_pt300(0),
hMeePtee_ULS_eta08_pt400(0),
hMeePtee_LS_eta08_pt400(0),
hMeePtee_ULS_eta08_pt200_opAngle50(0),
hMeePtee_LS_eta08_pt200_opAngle50(0),
hMeePtee_ULS_eta08_pt300_opAngle50(0),
hMeePtee_LS_eta08_pt300_opAngle50(0),
hMeePtee_ULS_eta08_pt400_opAngle50(0),
hMeePtee_LS_eta08_pt400_opAngle50(0),
hMeeOpAngle_ULS_eta08_pt200(0),
hMeeOpAngle_LS_eta08_pt200(0),
hMotherPt_ULS_eta08_pt200(0),
hMotherPt_LS_eta08_pt200(0),
hMotherPt_ULS_eta08_pt400(0),
hMotherPt_LS_eta08_pt400(0),
hweightcnmD(0),
hweightcnmHFE(0),
hweightcnmee(0),
fOutputList(0)
{
  Info("AliAnalysisTaskCharm","Calling Constructor");
  // Output slot #1 writes into a TList container (nevents histogram)
  DefineInput(0, TChain::Class());
  DefineOutput(1,TList::Class());
}

//______________________________| Destructor
AliAnalysisTaskCharm::~AliAnalysisTaskCharm()
{
  // Destructor
  Info("~AliAnalysisTaskCharm","Calling Destructor");
  if (fOutputList) delete fOutputList;
}

//______________________________| User Output
void AliAnalysisTaskCharm::UserCreateOutputObjects()
{
  Info("AliAnalysisTaskCharm","CreateOutputObjects of task %s", GetName());
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  printf("\n\n========================================\n  Configuration of task: \n========================================\n\n");
  printf("  process type:        %d\n", fProcessType);
  printf("  select only event with one ccbar pair:    %s\n", fSelectoneccbar?"YES":"NO");
  printf("  select event with clean history:    %s\n", fSelectcleanhistory?"YES":"NO");
  printf("  high-pt cut:         %f\n", fPtCutHigh);
  printf("  low-pt cut:         %f\n", fPtCutLow);
  printf("  etamin %f and etamax %f\n", fEtamin, fEtamax);
  printf("  use CNM scaling:    %s\n", fScaleByCNM?"YES":"NO");
  printf("  Take pt of D meson:    %s\n", fTakeptOfDCNM?"YES":"NO");
  printf("  use R_AA scaling:    %s\n", fScaleByRAA?"YES":"NO");
  printf("  use Decay weight:    %s\n", fApplywm?"YES":"NO");
  printf("  use Event weight:    %s\n", fEventWeight?"YES":"NO");
 
  std::cout << std::endl;

  AliCocktailSmearing::Print();

  CreateHistos();
  fNbEvents = 0;

  PostData(1, fOutputList);
}

//______________________________| Init
void AliAnalysisTaskCharm::Init()
{
  if(fDebug > 1) printf("AliAnalysisTaskCharm::Init() \n");
  fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}


//______________________________| User Exec
void AliAnalysisTaskCharm::UserExec(Option_t *)
{


  //
  //Nb event
  //
  
  fNbEventCounter++;
  if(fNbEvent>0 && fNbEventCounter<fNbEvent) return;

  if(fDebug > 1) printf("AliAnalysisTaskCharm::UserExec \n");
  Init();

  if(fMcHandler){
    fMcEvent = fMcHandler->MCEvent();
  }else{
    if(fDebug > 1) printf("AliAnalysisTaskCharm::Handler() fMcHandler=NULL\n");
    return;
  }
  if(!fMcEvent){
    if(fDebug > 1) printf("AliAnalysisTaskCharm::UserExec()   fMcEvent=NULL \n");
    return;
  }


  int nparticles = fMcEvent->GetNumberOfTracks();
 
  double opAngle50   =0.05; // 50 mrad

  // one can intialize it out of the event loop (if one reset this in the event loop after each event) :
  bool IsDmeson_totaly    = kFALSE;
  bool IsDmeson_y1        = kFALSE;
  bool IsCharm_totaly     = kFALSE;
  bool IsCharm_y1         = kFALSE;
  bool IsElectron_totaly  = kFALSE;
  bool IsElectron_y1      = kFALSE;
  
  //
  //bool foundanticharm = kFALSE;
  bool foundbeauty = kFALSE;
  Bool_t oneccbarpair = kFALSE;
  //

  //// Debug at quark level: two methods
  // Search c quarks which fragment
  Int_t Method1ncquarkintheeventp = 0;
  Int_t Method1ncquarkintheeventn = 0;
  Double_t Method1yquarkp = -1000.;
  Double_t Method1yquarkn = -1000.;
  Double_t Method1phiquarkp = -1000.;
  Double_t Method1phiquarkn = -1000.;
  Double_t Method1ptquarkp = -1000.;
  Double_t Method1ptquarkn = -1000.;
  AliVParticle *Method1partquarkp = 0x0;
  AliVParticle *Method1partquarkn = 0x0;
  // Search first two c quarks in the tree
  Int_t Method2ncquarkintheeventp = 0;
  Int_t Method2ncquarkintheeventn = 0;
  Double_t Method2yquarkp = -1000.;
  Double_t Method2yquarkn = -1000.;
  Double_t Method2phiquarkp = -1000.;
  Double_t Method2phiquarkn = -1000.;
  Double_t Method2ptquarkp = -1000.;
  Double_t Method2ptquarkn = -1000.;
  AliVParticle *Method2partquarkp = 0x0;
  AliVParticle *Method2partquarkn = 0x0;
  // Two correlated electrons
  Int_t highptDmesons = 0;
  Int_t Numberofelectron = 0;
  Int_t Numberofpositron = 0;
  Double_t etaelectron = -1000.;
  Double_t etapositron = -1000.;
  Double_t phielectron = -1000.;
  Double_t phipositron = -1000.;

 
  //printf("how many particles %d\n",nparticles);

  //// First loop over particles
  //============================ single loop =================================

  Float_t crosssection = 1.;
  // Find process + cross section
  AliGenEventHeader* header = fMcEvent->GenEventHeader();
  if(header) {
    //printf("Find the header %s\n",header->ClassName());
    AliGenPythiaEventHeader *pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(header);
    if(pythiaHeader){
      //printf("Pythia header found\n");
      Int_t process = pythiaHeader->ProcessType();
      //printf("Process type is %d and select %d\n",process,fProcessType);
      if((fProcessType>0) && (process!=fProcessType)) return;
      crosssection = pythiaHeader->GetXsection();
      //printf("Get cross section %f\n",crosssection);
    }
  }
  else printf("Did not find the header\n");

  Double_t eventw = 1.;
  Double_t eventcnm = 1.;
  if(fEventWeight) eventw = crosssection;

  //
  // First loop to check the number of heavy quarks in the events
  // To select events with only one ccbar and/or events with clean history
  //
  Bool_t elecwithhistory = kTRUE;
  for(int iparticle=0; iparticle<nparticles;iparticle++){

    AliVParticle * p = fMcEvent->GetTrack(iparticle);
    if (!p) continue;
    int pdg = TMath::Abs( p->PdgCode() );


    // Found charm quarks
    if(pdg==4){
      int k1 = p->GetDaughterFirst();
      int k2 = p->GetDaughterLast();
      // Look for charm quarks which fragments
      for(int d=k1; d <= k2; d++) {
        if(d>0){
          AliVParticle *decay = fMcEvent->GetTrack(d);
	  if(!decay) continue;
          int pdgdecay = decay->PdgCode();
          if ( int(TMath::Abs(pdgdecay)/100.) == 4 || int(TMath::Abs(pdgdecay)/1000.) == 4 ) {
            if(p->PdgCode()==4) {
              Method1ncquarkintheeventp ++;
	      Method1ptquarkp = p->Pt();
	      //printf("pdg decay %d (%d) for 4 quark %d\n",pdgdecay,d,iparticle);
            }
            if(p->PdgCode()==-4) {
              Method1ncquarkintheeventn ++;
	      Method1ptquarkn = p->Pt();
	      //printf("pdg decay %d (%d) for -4 quark %d\n",pdgdecay,d,iparticle);
            }
	  }
	}
      }
    }
    

    // Found electrons with history
    if ((pdg == 11) && (!p->IsPrimary())) { 
      
	int num = p->GetMother();
	AliVParticle *mom = fMcEvent->GetTrack ( num );
	if(!mom) continue;
	int ppid = TMath::Abs( mom->PdgCode() );

	// fill the single electron histogram, if mother is charm and doughter is electron
	bool is_charm = kFALSE;
	bool is_beauty2charm = kFALSE;
	bool is_cquark = kFALSE;
	//Int_t icquark = -1;
	if (((abs(ppid)>=400) && (abs(ppid)<=439)) || ((abs(ppid)>=4000) && (abs(ppid)<=4399))) is_charm = kTRUE;
	// Check if the D meson comes from B
	if(is_charm == kTRUE) {
	  AliVParticle *gm_i = 0x0;
	  Int_t indexx_i = mom->GetMother(); // First mother of D
	  while (indexx_i > 0) { // recursive loop to check if it comes from beauty.
	    gm_i = fMcEvent->GetTrack ( indexx_i );
	    if(gm_i) {
	      int pid_gm_i = TMath::Abs( gm_i->PdgCode() );// pdg of the mother
	      if (((pid_gm_i>=500) && (pid_gm_i<=549)) || ((pid_gm_i>=5000) && (pid_gm_i<=5499)) || (pid_gm_i == 5)) {
		is_beauty2charm = kTRUE;
	      }
	      if(pid_gm_i==4){
		is_cquark = kTRUE;
		//icquark = indexx_i;
	      }
	    indexx_i = gm_i->GetMother(); //
	    }
	    else indexx_i = -1;
	  }
	}
	if(is_charm && !is_beauty2charm) {
	  //printf("Found one HFE %d from open heavy-flavour %d with c quark %d\n",iparticle,ppid,icquark);
	  if(!is_cquark) elecwithhistory = kFALSE;
	}
    }
  }

  //printf("Number of charm quarks found with heavy-flavour hadrons as daughter %d and %d\n",Method1ncquarkintheeventp,Method1ncquarkintheeventn);
  // Reject events without one ccbar pair
  if((Method1ncquarkintheeventp==1) && (Method1ncquarkintheeventn==1)) oneccbarpair = kTRUE;
  if(!oneccbarpair && fSelectoneccbar) {
    PostData(1, fOutputList);
    return;
  }

  //Reject events without clear history
  if(!elecwithhistory && fSelectcleanhistory) {
    PostData(1, fOutputList);
    return;
  }

  // Calculate the CNM event weight if select only one ccbar pair
  if(fScaleByCNM && fSelectoneccbar && !fSelectcleanhistory){
    eventcnm = 0.5*(scale_CNM(Method1ptquarkp)+scale_CNM(Method1ptquarkn));
    //printf("CNM mean weight %f\n",eventcnm);
  }

  // Reset the variable for the next loop
  Method1ncquarkintheeventp = 0;
  Method1ncquarkintheeventn = 0;
  Method1ptquarkp = -1000.;
  Method1ptquarkn = -1000.;

  ///////////////////////////////////////////////////////////////////////////////
  // Second loop for the quark level for the normalization
  // In general no CNM weight used there to be consistent in the normalization
  // One histo to cross check consistency: hPtCharmQuarkMethodCNM
  //////////////////////////////////////////////////////////////////////////////

  for(int iparticle=0; iparticle<nparticles;iparticle++){

    AliVParticle * p = fMcEvent->GetTrack(iparticle);
    if (!p) continue;
    int pdg = TMath::Abs( p->PdgCode() );

    double px  = p->Px();
    double py  = p->Py();
    //double pz  = p->Pz();
    //double E   = p->Energy();
    //double eta = p->Eta();
    double rap = p->Y();
    double pt  = sqrt(px*px+py*py);


    //-------------------------- if Dmeson  ---------------------//
    if(pdg == 411  || pdg == 421  || pdg == 431 || pdg ==4122 ) {
      IsDmeson_totaly=kTRUE;
      // in |y|<1
      if(fabs(rap)<1.) IsDmeson_y1=kTRUE;
    }
    //-----------------------------------------------------------//

    //------------------ if charm quark -----------------//
    if (pdg == 4)  IsCharm_totaly=kTRUE;                 //
    // in |y|<1
    if (pdg == 4 && fabs(rap)<1.)  IsCharm_y1=kTRUE;//
    //---------------------------------------------------//

    // Found beauty quarks
    if(pdg==5) foundbeauty=kTRUE;

    // Found the charm quarks
    if(pdg==4){
      // Method 2
      if((p->PdgCode()==4) && (Method2ncquarkintheeventp==0)) {
        Method2ncquarkintheeventp ++;
        Method2yquarkp = rap;
	Method2phiquarkp = p->Phi();
	Method2ptquarkp = p->Pt();
	Method2partquarkp = p;
        hRapCharmQuarkMethod2->Fill(rap,eventw);
	hQuarkMethod2->Fill(rap,p->Pt(),eventw);
      }
      if((p->PdgCode()==-4) && (Method2ncquarkintheeventn==0)) {
        Method2ncquarkintheeventn ++;
        Method2yquarkn = rap;
	Method2phiquarkn = p->Phi();
	Method2ptquarkn = p->Pt();
	Method2partquarkn = p;
        hRapCharmQuarkMethod2->Fill(rap,eventw);
	hQuarkMethod2->Fill(rap,p->Pt(),eventw);
      }
      // Method 1 now
      //printf("Found a c quark %d with rapidity %f and pseudorapidity %f\n",p->PdgCode(),rap,eta);
      int k1 = p->GetDaughterFirst();
      int k2 = p->GetDaughterLast();
      //printf("k1 %d and k2 %d\n",k1,k2);
      bool foundhadrons = kFALSE;
      for(int d=k1; d <= k2; d++) {
        if(d>0){
          AliVParticle *decay = fMcEvent->GetTrack(d);
	  if(!decay) continue;
          int pdgdecay = decay->PdgCode();
          //printf("Mesons %d and pdg %d, index %d\n",d,pdgdecay,d);
          if ( int(TMath::Abs(pdgdecay)/100.) == 4 || int(TMath::Abs(pdgdecay)/1000.) == 4 ) {
            foundhadrons = kTRUE;
            if(p->PdgCode()==4) {
              Method1ncquarkintheeventp ++;
              Method1yquarkp = rap;
	      Method1phiquarkp = p->Phi();
	      Method1ptquarkp = p->Pt();
	      Method1partquarkp= p;
            }
            if(p->PdgCode()==-4) {
              Method1ncquarkintheeventn ++;
              Method1yquarkn = rap;
	      Method1phiquarkn = p->Phi();
	      Method1ptquarkn = p->Pt();
	      Method1partquarkn = p;
            }
            //printf("Found a quark with rapidity %f and pdg %d\n",rap,p->PdgCode());
            double rapD = decay->Y();
	    double ptD = decay->Pt();
            int kk1 = decay->GetDaughterFirst();
            int kk2 = decay->GetDaughterLast();
            //cout << "first doughter decay: " << kk1 << " last doughter decay: " << kk2 << endl;
            hRapCquarkDmeson->Fill(rap,rapD,eventw);
	    if(TMath::Abs(rap)<1.) hPtCquarkDmeson->Fill(pt,ptD,eventw);
            // Daughter of D mesons or D baryons
            for(int dd=kk1; dd <= kk2; dd++) {
              if(dd>0){
                AliVParticle *ddecay = fMcEvent->GetTrack(dd);
		if(!ddecay) continue;
                int pdgddecay = ddecay->PdgCode();
                //printf("Leptons %d and pdg %d, IsPrimary %d\n",dd,pdgddecay,ddecay->IsPrimary());
                if ( int(TMath::Abs(pdgddecay)) == 11 ) {
		  //printf("Found the decay electron, is it primary %d\n",(Int_t) ddecay->IsPrimary());
                  double rape = ddecay->Y();
		  double pte = ddecay->Pt();
                  hRapDmesonElectron->Fill(rapD,rape,eventw);
		  if(TMath::Abs(rap)<1.) hPtDmesonElectron->Fill(ptD,pte,eventw);
                }
		if(pdgddecay==11) {
		  Numberofpositron ++;
		  if(TMath::Abs(decay->Pt()) > 3.) highptDmesons++;
		  if(Numberofpositron==1) {
		    etapositron = ddecay->Eta();
		    phipositron = ddecay->Phi();
		  }
		}
		if(pdgddecay==-11) {
		  Numberofelectron ++;
		  if(TMath::Abs(decay->Pt()) > 3.) highptDmesons++;
		  if(Numberofelectron==1) {
		    etaelectron = ddecay->Eta();
		    phielectron = ddecay->Phi();
		  }
		}
		// Grand daughter of D mesons or D baryons (in case of D*(2007)->D0 for example)
		if ( int(TMath::Abs(pdgddecay)/100.) == 4 || int(TMath::Abs(pdgddecay)/1000.) == 4 ) {
		  int kkk1 = ddecay->GetDaughterFirst();
		  int kkk2 = ddecay->GetDaughterLast();
		  for(int ddd=kkk1; ddd <= kkk2; ddd++) {
		    if(ddd>0){
		      AliVParticle *dddecay = fMcEvent->GetTrack(ddd);
		      if(!dddecay) continue;
		      int pdgdddecay = dddecay->PdgCode();
		      //printf("New Leptons %d and pdg %d, IsPrimary %d\n",ddd,pdgdddecay,dddecay->IsPrimary());
		      if(pdgdddecay==11) {
			//printf("Found the decay electron, is it primary %d\n",(Int_t) ddecay->IsPrimary());
			Numberofpositron ++;
			if(TMath::Abs(ddecay->Pt()) > 3.) highptDmesons++;
			if(Numberofpositron==1) {
			  etapositron = dddecay->Eta();
			  phipositron = dddecay->Phi();
			}
		      }
		      if(pdgdddecay==-11) {
			//printf("Found the decay electron, is it primary %d\n",(Int_t) ddecay->IsPrimary());
			Numberofelectron ++;
			if(TMath::Abs(ddecay->Pt()) > 3.) highptDmesons++;
			// Take the first pairs found
			if(Numberofelectron==1) {
			  etaelectron = dddecay->Eta();
			  phielectron = dddecay->Phi();
			}
		      }
		    }
		  }
                }
              }
            }
          }
        }
      }
      //printf("Found hadrons %d\n",foundhadrons);
      if(foundhadrons)  {
	hQuarkMethod1->Fill(rap,p->Pt(),eventw);
	hRapCharmQuarkMethod1  ->Fill(rap,eventw);
	hPtCharmQuarkMethod1->Fill(pt,eventw);
	if(fScaleByCNM){
	  Double_t wcnm = 1.;
	  if(fSelectoneccbar && !fSelectcleanhistory){
	    wcnm = eventcnm;}
	  else {
	    wcnm = scale_CNM(pt);
	  }
	  hPtCharmQuarkMethod1CNM->Fill(pt,wcnm*eventw);
	}
      }
    }
  }

     
  /////////////////////////////////////////
  // Second loop for single D and e spectra
  //////////////////////////////////////////
  
  for(int iparticle=0; iparticle<nparticles;iparticle++){

    AliVParticle * p = fMcEvent->GetTrack(iparticle);
    if (!p) continue;
    int pdg = TMath::Abs( p->PdgCode() );

    double px  = p->Px();
    double py  = p->Py();
    //double pz  = p->Pz();
    //double E   = p->Energy();
    double eta = p->Eta();
    double rap = p->Y();
    double pt  = sqrt(px*px+py*py);


    //
    // D mesons spectra
    //

    // D from B mesons
    
    Bool_t is_Dmeson = kFALSE;
    Bool_t is_DfromBmeson = kFALSE;
    Double_t pt_cquark = -1.;
    Double_t pt_Dmeson = -1.;
    //Int_t i_c_quark = -1;
    if (((abs(pdg)>=400) && (abs(pdg)<=439)) || ((abs(pdg)>=4000) && (abs(pdg)<=4399))) {
      is_Dmeson = kTRUE;
      pt_Dmeson = p->Pt();
    }
    // Check if the D meson comes from B
    if(is_Dmeson == kTRUE) {
      //printf("Found one D mesons %d with pdg %d\n",iparticle,pdg);
      AliVParticle *gm_i = 0x0;
      Int_t indexx_i = p->GetMother(); // First mother of D
      while (indexx_i > 0) { // recursive loop to check if it comes from beauty.
	gm_i = fMcEvent->GetTrack ( indexx_i );
	if(gm_i) {
	  int pid_gm_i = TMath::Abs( gm_i->PdgCode() );// pdg of the mother
	  //printf("First mother %d with pdg %d\n",indexx_i,pid_gm_i);
	  if (((pid_gm_i>=500) && (pid_gm_i<=549)) || ((pid_gm_i>=5000) && (pid_gm_i<=5499)) || (pid_gm_i == 5)) {
	    is_DfromBmeson = kTRUE;
	  }
	  // Take the last open-charmed hadrons in the chain
	  if (((pid_gm_i>=400) && (pid_gm_i<=439)) || ((pid_gm_i>=4000) && (pid_gm_i<=4399))) {
	    pt_Dmeson = gm_i->Pt();
	  }
	  if((pid_gm_i==4) && (pt_cquark<0.)){
	    //i_c_quark = indexx_i;
	    pt_cquark = gm_i->Pt();
	  }
	  indexx_i = gm_i->GetMother(); //
	} else indexx_i = -1;
      }
    }
    
    // weight
    Double_t wcnm = 1.;
    if(fScaleByCNM && is_Dmeson && !is_DfromBmeson){
      if(fSelectoneccbar && !fSelectcleanhistory){
	wcnm = eventcnm;
      } 
      else if(!fTakeptOfDCNM) {
	if(pt_cquark > 0.){
	  wcnm = scale_CNM(pt_cquark);
	}
	else {
	  wcnm = -1.;
	}
      }
      else {
	if(pt_Dmeson > 0.){
	  wcnm = scale_CNM(pt_Dmeson);
	}
	else {
	  wcnm = -1.;
	}
      }
      if(!fTakeptOfDCNM) hweightcnmD->Fill(wcnm,pt_cquark);
      else hweightcnmD->Fill(wcnm,pt_Dmeson);
    }

    // **** default ****
    
    // fill the single histograms
    // mesons w/ charm quarks
    if(is_Dmeson && !is_DfromBmeson)         hCharm        ->Fill(rap,pt,eventw*wcnm);
    //D^{pm}
    if(pdg == 411 && !is_DfromBmeson)        hDpm          ->Fill(rap,pt,eventw*wcnm);
    //D0
    if(pdg == 421 && !is_DfromBmeson)        hD0           ->Fill(rap,pt,eventw*wcnm);
    //Ds
    if(pdg == 431 && !is_DfromBmeson)        hDs           ->Fill(rap,pt,eventw*wcnm);
    //Dstar
    if(pdg == 413 && !is_DfromBmeson)        hDstar        ->Fill(rap,pt,eventw*wcnm);
    //Lamdac
    if(pdg == 4122 && !is_DfromBmeson)       hLambdac      ->Fill(rap,pt,eventw*wcnm);


    /////********


      
    //------------->>>> electrons are daughters
    //
    if (!( pdg == 11 )) continue;
    if ( p->IsPrimary() ) continue;

      
    int num = p->GetMother();
    AliVParticle *mom = fMcEvent->GetTrack ( num );
    if(!mom) continue;
    int ppid = TMath::Abs( mom->PdgCode() );

    //double pxMom  = mom->Px();
    //double pyMom  = mom->Py();
    double rapMom = mom->Y();
    //double ptMom  = sqrt(pxMom*pxMom+pyMom*pyMom);

    // fill the single electron histogram, if mother is charm and doughter is electron
    bool is_charm = kFALSE;
    bool is_beauty2charm = kFALSE;
    pt_cquark = -1;
    //i_c_quark = -1;
    pt_Dmeson = -1;
    if (((abs(ppid)>=400) && (abs(ppid)<=439)) || ((abs(ppid)>=4000) && (abs(ppid)<=4399))) {
      is_charm = kTRUE;
      pt_Dmeson = mom->Pt();
    }
    // Check if the D meson comes from B
    if(is_charm == kTRUE) {
      AliVParticle *gm_i = 0x0;
      Int_t indexx_i = mom->GetMother(); // First mother of D
      while (indexx_i > 0) { // recursive loop to check if it comes from beauty.
	gm_i = fMcEvent->GetTrack ( indexx_i );
	if(gm_i) {
	  int pid_gm_i = TMath::Abs( gm_i->PdgCode() );// pdg of the mother
	  if (((pid_gm_i>=500) && (pid_gm_i<=549)) || ((pid_gm_i>=5000) && (pid_gm_i<=5499)) || (pid_gm_i == 5)) {
	    is_beauty2charm = kTRUE;
	  }
	  if (((pid_gm_i>=400) && (pid_gm_i<=439)) || ((pid_gm_i>=4000) && (pid_gm_i<=4399))) {
	    pt_Dmeson = gm_i->Pt();
	  }
	  if((pid_gm_i==4) && (pt_cquark < 0.)){
	    //i_c_quark = indexx_i;
	    pt_cquark = gm_i->Pt();
	  }
	  indexx_i = gm_i->GetMother(); //
	} else indexx_i = -1;
      }
    }


    // Weight with Branching ratios
    Double_t wm = 1.;
    if(fApplywm) {
      if(ppid==411) wm = 0.1607;
      else if(ppid==421) wm = 0.0649;
      else if(ppid==431) wm = 0.065; 
      else if(ppid==4122) wm = 0.036;
      else wm = 0.;
    }
    wm = wm*eventw;
    
    
    // weight CNM
    if(fScaleByCNM && is_charm && !is_beauty2charm){
      Double_t wcnme = 1.;
      if(fSelectoneccbar && !fSelectcleanhistory){
	wcnme = eventcnm;
	// For check
	//if(pt_cquark > 0.){
	  //printf("HFE: mean weight %f and single weight %f\n",wcnme,scale_CNM(pt_cquark));
	//}
      } else if(!fTakeptOfDCNM) {
	if(pt_cquark>0){
	  wcnme = scale_CNM(pt_cquark);
	}
	else {
	  wcnme = -1.;
	  //printf("HFE %d with c quark %d\n",iparticle,i_c_quark);
	}
      } else {
	if(pt_Dmeson>0){
	  wcnme = scale_CNM(pt_Dmeson);
	}
	else {
	  wcnme = -1.;
	  //printf("HFE %d with c quark %d\n",iparticle,i_c_quark);
	}
      }
      //printf("HFE: wcnme %f\n",wcnme);
      wm = wm*wcnme;
      if(!fTakeptOfDCNM) hweightcnmHFE->Fill(wcnme,pt_cquark);
      else hweightcnmHFE->Fill(wcnme,pt_Dmeson);
    }
   
    if(is_charm && !is_beauty2charm){
     
      hRapDMom2e->Fill(rapMom,rap,wm);
      hPtEtaElectron->Fill(pt, eta,wm); // generated pt vs eta of electrons from charm
      if(p->PdgCode() ==   11) hPtEtaElectronE->Fill(pt, eta,wm); // generated pt vs eta of electrons from charm
      if(p->PdgCode() ==   -11) hPtEtaElectronP->Fill(pt, eta,wm); // generated pt vs eta of electrons from charm

      //--rapidity and pseudorapidity distributions of electrons for pt> 0.5 GeV/c ---//
      if(pt>0.2) hRapElectronPt200->Fill(rap,wm);
      if(pt>0.5) hRapElectronPt500->Fill(rap,wm);
      //-- single leg pt distributions of electrons ----------------------------------//
      if(fabs(eta)<0.8)  hPte_eta08 ->Fill(pt,wm); // all e  in |eta|<0.8
      if(p->PdgCode() ==   11 && fabs(eta)<0.8)  hPteP_eta08->Fill(pt,wm); // all e+ in |eta|<0.8
      if(p->PdgCode() ==  -11 && fabs(eta)<0.8)  hPteM_eta08->Fill(pt,wm); // all e- in |eta|<0.8
      if(fabs(rap)<0.8)  hPte_y08   ->Fill(pt,wm); // all e  in |y|<0.8
      if(p->PdgCode() ==  11  && fabs(rap)<0.8)  hPteP_y08  ->Fill(pt,wm); // all e+ in |y|<0.8
      if(p->PdgCode() == -11  && fabs(rap)<0.8)  hPteM_y08  ->Fill(pt,wm); // all e+ in |y|<0.8

      //-------------------- if electron from charm --------------//
      IsElectron_totaly=kTRUE;//
      // if electron from charm, in |y|<1:
      if(fabs(rap)<1.) IsElectron_y1=kTRUE;    //
      //----------------------------------------------------------//
    }

  }//end of single loop

  //
  //================================== D meson correlation =================================
  // For LHCb paper at 7 TeV
  // 
  for(int i=1; i<nparticles;i++) {
    AliVParticle * p_i = fMcEvent->GetTrack(i);
    if (!p_i) continue;
   
    int pid_i = p_i->PdgCode();
    if( ! ((TMath::Abs( pid_i ) == 421) || (TMath::Abs( pid_i ) == 411) || (TMath::Abs( pid_i ) == 431)) ) continue;

    // Check not from beauty
    int num_i = p_i->GetMother();
    AliVParticle *mom_i = fMcEvent->GetTrack ( num_i );
    if(!mom_i) continue;
    int ppid_i = TMath::Abs( mom_i->PdgCode());
    if ( int(TMath::Abs(ppid_i)/100.) == 5 || int(TMath::Abs(ppid_i)/1000.) == 5 ) continue;

    double px_i  = p_i->Px();
    double py_i  = p_i->Py();
    //double eta_i = p_i->Eta();
    double rap_i = p_i->Y();
    double pt_i  = sqrt(px_i*px_i+py_i*py_i);
    //double phi_i = TVector2::Phi_mpi_pi(p_i->Phi());
    //double phicheck_i = p_i->Phi();

   

    if((pt_i < 3) || (pt_i > 12)) continue;

    int k1_i = p_i->GetDaughterFirst();
    int k2_i = p_i->GetDaughterLast();
    Int_t founde_i = 0;
    TLorentzVector ppoe_i;
    for(int d=k1_i; d <= k2_i; d++) {
      if(d>0){
	AliVParticle *decay = fMcEvent->GetTrack(d);
	if(!decay) continue;
	int pdgdecay = decay->PdgCode();
	if(TMath::Abs(pdgdecay)==11) {
	  decay->Momentum(ppoe_i);
	  founde_i = pdgdecay; 
	}
      }
    }

    for(int j=i+1; j<nparticles;j++) {
      AliVParticle * p_j = fMcEvent->GetTrack(j);
      if(!p_j) continue;
      
      int pid_j = p_j->PdgCode();
      //if ( p_j->IsPrimary() ) continue;
      if( ! ((TMath::Abs( pid_j ) == 421) || (TMath::Abs( pid_j ) == 411) || (TMath::Abs( pid_j ) == 431)) ) continue;

      // Check not from beauty
      int num_j = p_j->GetMother();
      AliVParticle *mom_j = fMcEvent->GetTrack ( num_j );
      if(!mom_j) continue;
      int ppid_j = TMath::Abs( mom_j->PdgCode());
      if ( int(TMath::Abs(ppid_j)/100.) == 5 || int(TMath::Abs(ppid_j)/1000.) == 5 ) continue;
      
      double px_j  = p_j->Px();
      double py_j  = p_j->Py();
      //double eta_j = p_j->Eta();
      double rap_j = p_j->Y();
      double pt_j  = sqrt(px_j*px_j+py_j*py_j);
      //double phi_j = TVector2::Phi_mpi_pi(p_j->Phi());
      //double phicheck_j = p_j->Phi();
      

      if((pt_j < 3) || (pt_j > 12)) continue;

      TLorentzVector ppo_i,ppo_j;
      p_i->Momentum(ppo_i);
      p_j->Momentum(ppo_j);

      Double_t deltaphie = TMath::Abs(ppo_i.DeltaPhi(ppo_j));

      int k1_j = p_j->GetDaughterFirst();
      int k2_j = p_j->GetDaughterLast();
      Int_t founde_j = 0;
      TLorentzVector ppoe_j;
      for(int d=k1_j; d <= k2_j; d++) {
	if(d>0){
	  AliVParticle *decay = fMcEvent->GetTrack(d);
	  if(!decay) continue;
	  int pdgdecay = decay->PdgCode();
	  if(TMath::Abs(pdgdecay)==11) {
	    decay->Momentum(ppoe_j);
	    founde_j = pdgdecay; 
	  }
	}
      }

      //Double_t deltaphie = TMath::Abs(TVector2::Phi_mpi_pi(phi_i-phi_j));
      //Double_t deltaphiecheck = TMath::Abs(TVector2::Phi_mpi_pi(phicheck_i-phicheck_j));
      //printf("deltaphi %f and %f and %f\n",deltaphie,deltaphiecheck,deltaphiv);

      // ULS
      if(pid_i * pid_j <0.){
	  hDeltaPhi_D0->Fill(deltaphie);
	  if((rap_i > 2) && (rap_i < 4) && (rap_j > 2) && (rap_j < 4)) {
	    hDeltaPhi_D0_LHCb->Fill(deltaphie);
	  }
	}

      // ULS decay
      if(founde_i * founde_j <0.){
	if((rap_i > 2) && (rap_i < 4) && (rap_j > 2) && (rap_j < 4)) {
	  Double_t deltaphiee = TMath::Abs(ppoe_i.DeltaPhi(ppoe_j));
	    hDeltaPhi_D0_LHCb_e->Fill(deltaphiee);
	}
      }

    }
  }

  //
  //================================== pair loop =================================
  // For dielectron final spectra 
  //

  for(int i=1; i<nparticles;i++) {
    AliVParticle * p_i = fMcEvent->GetTrack(i);
    if (!p_i) continue;
    if ( p_i->IsPrimary() ) continue;

    int pid_i = p_i->PdgCode();
    if( ! (TMath::Abs( pid_i ) == 11) ) continue;


    for(int j=i+1; j<nparticles;j++) {
      AliVParticle * p_j = fMcEvent->GetTrack(j);
      if (!p_j) continue;
      if ( p_j->IsPrimary() ) continue;

      int pid_j = p_j->PdgCode();
      if( ! (TMath::Abs( pid_j ) == 11) ) continue;

      //-----------------------------------------
      // variable without bremsstrahlung+momentum resolution
      TLorentzVector ppo_i,ppo_j;
      p_i->Momentum(ppo_i);
      p_j->Momentum(ppo_j);

      double etao_i  = ppo_i.Eta();
      double pto_i   = ppo_i.Pt();
      double phio_i  = TVector2::Phi_0_2pi(ppo_i.Phi());

      double etao_j  = ppo_j.Eta();
      double pto_j   = ppo_j.Pt();
      double phio_j  = TVector2::Phi_0_2pi(ppo_j.Phi());

      double masso    = (ppo_i + ppo_j).M();
      //double pto_pair = (ppo_i + ppo_j).Pt();
      //double opAngleo = ppo_i.Angle(ppo_j.Vect());
      //double deltaphio = TMath::Abs(ppo_i.DeltaPhi(ppo_j));

      //-----------------------------------------
      // apply bremsstrahlung+momentum resolution
      //-----------------------------------------
      Short_t ch_i = 1;
      if(pid_i>0) ch_i = -1; // 11 is electron, -11 is positron
      Short_t ch_j = 1;
      if(pid_j>0) ch_j = -1; // 11 is electron, -11 is positron
      TLorentzVector pp_i,pp_j;
      if(!fOton) {
	pp_i = Smear(p_i);
	pp_j = Smear(p_j);
	// smearing of p is only linear, so we need to smear the opening angle explicitly in addition.
	// in order not to overdo the corrections, the vectors are just rotated in this procedure.
	SmearOpeningAngle(pp_i, pp_j);
      } else {
	pp_i = ApplySmearingOton(ppo_i,ch_i);
	pp_j = ApplySmearingOton(ppo_j,ch_j);
      }
      

      double eta_i  = pp_i.Eta();
      double pt_i   = pp_i.Pt();
      //double phi_i  = pp_i.Phi();

      double eta_j  = pp_j.Eta();
      double pt_j   = pp_j.Pt();
      //double phi_j  = pp_j.Phi();

      double mass     = (pp_i + pp_j).M();
      double pt_pair  = (pp_i + pp_j).Pt();
      double opAngle  = pp_i.Angle(pp_j.Vect());
      double deltaphi = TMath::Abs(pp_i.DeltaPhi(pp_j));
      //-----------------------------------------


      int num_i = p_i->GetMother();
      AliVParticle *mom_i = fMcEvent->GetTrack ( num_i );
      if(!mom_i) continue;
      int ppid_i = TMath::Abs( mom_i->PdgCode() );
      double ppx_i  = mom_i->Px();
      double ppy_i  = mom_i->Py();
      double ppt_i  = sqrt(ppx_i*ppx_i+ppy_i*ppy_i);
      

      int num_j = p_j->GetMother();
      AliVParticle *mom_j = fMcEvent->GetTrack ( num_j );
      if(!mom_j) continue;
      int ppid_j = TMath::Abs( mom_j->PdgCode() );
      double ppx_j  = mom_j->Px();
      double ppy_j  = mom_j->Py();
      double ppt_j  = sqrt(ppx_j*ppx_j+ppy_j*ppy_j);
      
       // i electron
      bool i_is_charm = kFALSE;
      bool i_is_beauty2charm = kFALSE;
      Double_t i_pt_c = -1.;
      Double_t i_pt_D = -1.;
      if (((abs(ppid_i)>=400) && (abs(ppid_i)<=439)) || ((abs(ppid_i)>=4000) && (abs(ppid_i)<=4399))) {
	i_is_charm = kTRUE;
	i_pt_D = mom_i->Pt();
      }
      // Check if the D meson comes from B
      if(i_is_charm == kTRUE) {
	AliVParticle *gm_i = 0x0;
	Int_t indexx_i = mom_i->GetMother(); // First mother of D
	while (indexx_i > 0) { // recursive loop to check if it comes from beauty.
	  gm_i = fMcEvent->GetTrack ( indexx_i );
	  if(gm_i) {
	    int pid_gm_i = TMath::Abs( gm_i->PdgCode() );// pdg of the mother
	    if (((pid_gm_i>=500) && (pid_gm_i<=549)) || ((pid_gm_i>=5000) && (pid_gm_i<=5499)) || (pid_gm_i == 5)) {
	      i_is_beauty2charm = kTRUE;
	    }
	    if (((pid_gm_i>=400) && (pid_gm_i<=439)) || ((pid_gm_i>=4000) && (pid_gm_i<=4399))) {
	      i_pt_D = gm_i->Pt();
	    }
	    if((pid_gm_i==4) && (i_pt_c < 0.)) {
	      i_pt_c = gm_i->Pt();
	    }
	    indexx_i = gm_i->GetMother(); //
	  } else indexx_i = -1;
	}
	if(i_is_beauty2charm == kTRUE) i_is_charm = kFALSE; // it is a charm that comes from beauty.
      }
      // j electron
      bool j_is_charm = kFALSE;
      bool j_is_beauty2charm = kFALSE;
      Double_t j_pt_c = -1.;
      Double_t j_pt_D = -1.;
      if (((abs(ppid_j)>=400) && (abs(ppid_j)<=439)) || ((abs(ppid_j)>=4000) && (abs(ppid_j)<=4399))) {
	j_is_charm = kTRUE;
	j_pt_D = mom_j->Pt();
      }
      // Check if the D meson comes from B
      if(j_is_charm == kTRUE) {
	AliVParticle *gm_j = 0x0;
	Int_t indexx_j = mom_j->GetMother(); // First mother of D
	while (indexx_j > 0) { // recursive loop to check if it comes from beauty.
	  gm_j = fMcEvent->GetTrack ( indexx_j );
	  if(gm_j){
	    int pid_gm_j = TMath::Abs( gm_j->PdgCode() );// pdg of the mother
	    if (((pid_gm_j>=500) && (pid_gm_j<=549)) || ((pid_gm_j>=5000) && (pid_gm_j<=5499)) || (pid_gm_j == 5)) {
	      j_is_beauty2charm = kTRUE;
	    }
	    if (((pid_gm_j>=400) && (pid_gm_j<=439)) || ((pid_gm_j>=4000) && (pid_gm_j<=4399))) {
	      j_pt_D = gm_j->Pt();
	    }
	    if((pid_gm_j==4) && (j_pt_c < 0.)){
	      j_pt_c = gm_j->Pt();
	    }
	    indexx_j = gm_j->GetMother(); //
	  } else indexx_j = -1;
	}
	if(j_is_beauty2charm == kTRUE) j_is_charm = kFALSE; // it is a charm that comes from beauty.
      }
      
      //
      Double_t wm_i = 1.;
      Double_t wm_j = 1.;
      Double_t wm = 1.;
      if(fApplywm) {
	if(ppid_i==411) wm_i = 0.1607;
	else if(ppid_i==421) wm_i = 0.0649;
	else if(ppid_i==431) wm_i = 0.065; 
	else if(ppid_i==4122) wm_i = 0.036;
	else wm_i = 0.;
	if(ppid_j==411) wm_j = 0.1607;
	else if(ppid_j==421) wm_j = 0.0649;
	else if(ppid_j==431) wm_j = 0.065; 
	else if(ppid_j==4122) wm_j = 0.036;
	else wm_j = 0.;
	wm = wm_i * wm_j;
      }
      wm = wm*eventw;


      //------------------


      // track pt cuts
      double ptweight2 = pt_cut200(pt_i) * pt_cut200(pt_j) * pt_cutHigh(pt_i) * pt_cutHigh(pt_j); // pT>0.2
      double ptweight3 = pt_cut300(pt_i) * pt_cut300(pt_j) * pt_cutHigh(pt_i) * pt_cutHigh(pt_j); // pT>0.3
      double ptweight4 = pt_cut400(pt_i) * pt_cut400(pt_j) * pt_cutHigh(pt_i) * pt_cutHigh(pt_j); // pT>0.4
      double ptweight5 = pt_cutLow(pt_i) * pt_cutLow(pt_j) * pt_cutHigh(pt_i) * pt_cutHigh(pt_j); // variable

       // R_AA CNM quark
      if (fScaleByCNM && i_is_charm && j_is_charm) {
	Double_t cnmw = 1.;
	if(fSelectoneccbar && !fSelectcleanhistory){
	  cnmw = eventcnm;
	  // For check
	  //if((i_pt_c > 0.) && (j_pt_c>0.)){
	    //printf("ee: mean weight %f and single weight %f\n",cnmw,0.5*(scale_CNM(i_pt_c)+scale_CNM(j_pt_c)));
	  //}
	} else if(!fTakeptOfDCNM) {
	  if((i_pt_c > 0.) && (j_pt_c>0.)) {
	    cnmw = 0.5*(scale_CNM(i_pt_c)+scale_CNM(j_pt_c));
	  }
	  else {
	    cnmw = -1.;
	    //printf("ee: Event with negative weight %d and pt1 %f pt2 %f\n",fNbEventCounter,i_pt_c,j_pt_c);
	  }
	} else {
	  if((i_pt_D > 0.) && (j_pt_D>0.)) {
	    cnmw = 0.5*(scale_CNM(i_pt_D)+scale_CNM(j_pt_D));
	  }
	  else {
	    cnmw = -1.;
	    //printf("ee: Event with negative weight %d and pt1 %f pt2 %f\n",fNbEventCounter,i_pt_c,j_pt_c);
	  }
	}
	//printf("ee: cnmw %f\n",cnmw);
	if(!fTakeptOfDCNM) {
	  hweightcnmee->Fill(cnmw,i_pt_c);
	  hweightcnmee->Fill(cnmw,j_pt_c);
	} else {
	  hweightcnmee->Fill(cnmw,i_pt_D);
	  hweightcnmee->Fill(cnmw,j_pt_D);
	}
	ptweight2 *= cnmw;
	ptweight3 *= cnmw;
	ptweight4 *= cnmw;
	ptweight5 *= cnmw;
      }
      
      
      // HFE R_AA scaling
      if (fScaleByRAA) {
        ptweight2 *= scale_RAA(pt_i) * scale_RAA(pt_j);
        ptweight3 *= scale_RAA(pt_i) * scale_RAA(pt_j);
        ptweight4 *= scale_RAA(pt_i) * scale_RAA(pt_j);
	ptweight5 *= scale_RAA(pt_i) * scale_RAA(pt_j);
      }
      // if not apply w or event w, then 1 for wm
      ptweight2 *= wm;
      ptweight3 *= wm;
      ptweight4 *= wm;
      ptweight5 *= wm;
      

      //--------------------------------------- ULS pairs -------------------------------------------//

      if(( pid_i * pid_j < 0.) ) { // +-pairs


        if(i_is_charm && j_is_charm) {

          hMee_ULS_simulated->Fill(mass,wm);

          //______________ histograms in phenix acceptance
          if(fabs(etao_i)<0.5  && fabs(etao_j)<0.5)    {
            hMeePtee_ULS_eta05->Fill(masso,pt_pair,wm);
          } //|eta|<0.5

          if(fabs(etao_i)<0.35  && fabs(etao_j)<0.35)  {
            hMeePtee_ULS_eta035->Fill(masso,pt_pair,wm);
            if(Inphenixacc(phio_i,pto_i,pid_i) && Inphenixacc(phio_j,pto_j,pid_j)) {
	      hMeePtee_ULS_eta035_phenixacc->Fill(masso,pt_pair,wm);// dphi cut, pt>0.2
	    }
          }//|eta|<0.35
          //_________________________________________________________

	  if((fEtamin < eta_i) && (eta_i < fEtamax)  && (fEtamin < eta_i) && (eta_i < fEtamax))  {
            hMeePtee_ULS_eta_pt->Fill(mass,pt_pair,ptweight5); // pt>fptmin
	  }
	  

          if(fabs(eta_i)<0.8  && fabs(eta_j)<0.8)  {
            hMee_ULS_eta08->Fill(mass,wm);
            hMee_ULS_eta08_pt200->Fill(mass,ptweight2); // pt>0.2
            // 2D:
            hMeePtee_ULS_eta08->Fill(mass,pt_pair,wm);
            hMeePtee_ULS_eta08_pt200->Fill(mass,pt_pair,ptweight2); // pt>0.2
            hMeePtee_ULS_eta08_pt300->Fill(mass,pt_pair,ptweight3); // pt>0.3
            hMeePtee_ULS_eta08_pt400->Fill(mass,pt_pair,ptweight4); // pt>0.4
            hMeeOpAngle_ULS_eta08_pt200->Fill(mass,opAngle,ptweight2); // opening angle
            if (opAngle > opAngle50) { // cut on smeared opening angle
              hMeePtee_ULS_eta08_pt200_opAngle50->Fill(mass,pt_pair,ptweight2); // pt>0.2
              hMeePtee_ULS_eta08_pt300_opAngle50->Fill(mass,pt_pair,ptweight3); // pt>0.3
              hMeePtee_ULS_eta08_pt400_opAngle50->Fill(mass,pt_pair,ptweight4); // pt>0.4
            }
	    // 3D
	    hPhiee_ULS_eta08_pt200->Fill(deltaphi,mass,pt_pair,ptweight2);
	    if((ppt_i>3) || (ppt_j>3))  hPhiee_ULS_eta08_highoneDmeson_pt200->Fill(deltaphi,mass,pt_pair,ptweight2);
	    if((ppt_i>3) && (ppt_j>3))  hPhiee_ULS_eta08_highbothDmeson_pt200->Fill(deltaphi,mass,pt_pair,ptweight2);
	    hPhiee_ULS_eta08_pt400->Fill(deltaphi,mass,pt_pair,ptweight4);
	    if((ppt_i>3) || (ppt_j>3))  hPhiee_ULS_eta08_highoneDmeson_pt400->Fill(deltaphi,mass,pt_pair,ptweight4);
	    if((ppt_i>3) && (ppt_j>3))  hPhiee_ULS_eta08_highbothDmeson_pt400->Fill(deltaphi,mass,pt_pair,ptweight4);
	    // D meson pt
	    hMotherPt_ULS_eta08_pt200->Fill(ppt_i,mass,pt_pair,ptweight2);
	    hMotherPt_ULS_eta08_pt200->Fill(ppt_j,mass,pt_pair,ptweight2);
	    hMotherPt_ULS_eta08_pt400->Fill(ppt_i,mass,pt_pair,ptweight4);
	    hMotherPt_ULS_eta08_pt400->Fill(ppt_j,mass,pt_pair,ptweight4);
          } //  |eta|<0.8

        }//charm
      }//unlike-sign


      //--------------------------------------- LS pairs -------------------------------------------//

      if(( pid_i * pid_j > 0.))  { // ++ and -- pairs


        if(i_is_charm && j_is_charm) {

          hMee_LS_simulated->Fill(mass,wm);

          /*cout << "================ LS ================" << endl;
           cout << "pid_i * pid_j: " << pid_i * pid_j << endl;
           cout << "mother_i: " << ppid_i << " " <<  "daughter_i: " << pid_i << endl;
           cout << "mother_j: " << ppid_j << " " <<  "daughter_j: " << pid_j << endl;
           */


          //______________ new histograms in phenix acceptance
          if(fabs(etao_i)<0.5  && fabs(etao_j)<0.5)    {
            hMeePtee_LS_eta05->Fill(masso,pt_pair,wm);
          } //|eta|<0.5

          if(fabs(etao_i)<0.35  && fabs(etao_j)<0.35)  {
            hMeePtee_LS_eta035->Fill(masso,pt_pair,wm);
	    if(Inphenixacc(phio_i,pto_i,pid_i) && Inphenixacc(phio_j,pto_j,pid_j)) {
	      hMeePtee_LS_eta035_phenixacc->Fill(masso,pt_pair,wm);// dphi cut, pt>0.2
	    }
          }//|eta|<0.35
          //_________________________________________________________
	  
	  if((fEtamin < eta_i) && (eta_i < fEtamax)  && (fEtamin < eta_i) && (eta_i < fEtamax))  {
            hMeePtee_LS_eta_pt->Fill(mass,pt_pair,ptweight5); // pt>fptmin
	  }

          if(fabs(eta_i)<0.8  && fabs(eta_j)<0.8)  {
            hMee_LS_eta08->Fill(mass,wm);
            hMee_LS_eta08_pt200->Fill(mass,ptweight2); // pt>0.2
            // 2D:
            hMeePtee_LS_eta08->Fill(mass,pt_pair,wm);
            hMeePtee_LS_eta08_pt200->Fill(mass,pt_pair,ptweight2); // pt>0.2
            hMeePtee_LS_eta08_pt300->Fill(mass,pt_pair,ptweight3); // pt>0.3
            hMeePtee_LS_eta08_pt400->Fill(mass,pt_pair,ptweight4); // pt>0.4
            hMeeOpAngle_LS_eta08_pt200->Fill(mass,opAngle,ptweight2); // opening angle
            if (opAngle > opAngle50) { // cut on smeared opening angle
              hMeePtee_LS_eta08_pt200_opAngle50->Fill(mass,pt_pair,ptweight2); // pt>0.2
              hMeePtee_LS_eta08_pt300_opAngle50->Fill(mass,pt_pair,ptweight3); // pt>0.3
              hMeePtee_LS_eta08_pt400_opAngle50->Fill(mass,pt_pair,ptweight4); // pt>0.4
            }
	    // 3D
	    hPhiee_LS_eta08_pt200->Fill(deltaphi,mass,pt_pair,ptweight2);
	    if((ppt_i>3) || (ppt_j>3))  hPhiee_LS_eta08_highoneDmeson_pt200->Fill(deltaphi,mass,pt_pair,ptweight2);
	    if((ppt_i>3) && (ppt_j>3))  hPhiee_LS_eta08_highbothDmeson_pt200->Fill(deltaphi,mass,pt_pair,ptweight2);
	    hPhiee_LS_eta08_pt400->Fill(deltaphi,mass,pt_pair,ptweight4);
	    if((ppt_i>3) || (ppt_j>3))  hPhiee_LS_eta08_highoneDmeson_pt400->Fill(deltaphi,mass,pt_pair,ptweight4);
	    if((ppt_i>3) && (ppt_j>3))  hPhiee_LS_eta08_highbothDmeson_pt400->Fill(deltaphi,mass,pt_pair,ptweight4);
	     // D meson pt
	    hMotherPt_LS_eta08_pt200->Fill(ppt_i,mass,pt_pair,ptweight2);
	    hMotherPt_LS_eta08_pt200->Fill(ppt_j,mass,pt_pair,ptweight2);
	    hMotherPt_LS_eta08_pt400->Fill(ppt_i,mass,pt_pair,ptweight4);
	    hMotherPt_LS_eta08_pt400->Fill(ppt_j,mass,pt_pair,ptweight4);
          } // |eta|<0.8
        }//charm
      }//like-sign


    }//end of pair loop for j^th particle

  }//end of pair loop for i^th particle

  hNEventsW->Fill(crosssection); // Fill with cross section
  hNEvents->Fill(0.5); // all events w/o any rapidity cut
  if(IsDmeson_totaly==kTRUE)     hNEvents->Fill(1.5); // all Dmeson events in whole y
  if(IsDmeson_y1==kTRUE)         hNEvents->Fill(2.5); // all Dmeson events in |y|<1
  if(IsCharm_totaly==kTRUE)      hNEvents->Fill(3.5); // all charm events in whole y
  if(IsCharm_y1==kTRUE)          hNEvents->Fill(4.5); // all charm events in |y|<1
  if(IsElectron_totaly==kTRUE)   hNEvents->Fill(5.5); // all electrons events in whole y
  if(IsElectron_y1==kTRUE)       hNEvents->Fill(6.5); // all electrons in |y|<1 from charm
  if(foundbeauty==kTRUE)         hNEvents->Fill(7.5); // beauty quarks in the event
  fNbEvents++;


  if((Method1ncquarkintheeventp==1) && (Method1ncquarkintheeventn==1)) {
    hRapCharmQuark2DMethod1->Fill(Method1yquarkp,Method1yquarkn,eventw);
    if((TMath::Abs(Method1yquarkp)<1.) && (TMath::Abs(Method1yquarkn)<1.)) hPtCharmQuark2DMethod1->Fill(Method1ptquarkp,Method1ptquarkn,eventw);
    Double_t deltaphiquark = TVector2::Phi_mpi_pi(Method1phiquarkp-Method1phiquarkn);
    hPhiCharmQuark1DMethod1->Fill(deltaphiquark,eventw);
    // Rapidity of the pair
    if(Method1partquarkp && Method1partquarkn) {
      TLorentzVector q_n,q_p;
      Method1partquarkp->Momentum(q_n);
      Method1partquarkn->Momentum(q_p);
      Double_t rapiditypairquark    = (q_n + q_p).Rapidity();
      hRapCharmQuarkMethod1Pair  ->Fill(rapiditypairquark,eventw);
    }
  }
  if((Method2ncquarkintheeventp==1) && (Method2ncquarkintheeventn==1)) {
    hRapCharmQuark2DMethod2->Fill(Method2yquarkp,Method2yquarkn,eventw);
    if((TMath::Abs(Method2yquarkp)<1.) && (TMath::Abs(Method2yquarkn)<1.)) hPtCharmQuark2DMethod2->Fill(Method2ptquarkp,Method2ptquarkn,eventw);
    Double_t deltaphiquark = TVector2::Phi_mpi_pi(Method2phiquarkp-Method2phiquarkn);
    hPhiCharmQuark1DMethod2->Fill(deltaphiquark,eventw);
    // Rapidity of the pair
    if(Method2partquarkp && Method2partquarkn) {
      TLorentzVector q_n,q_p;
      Method2partquarkp->Momentum(q_n);
      Method2partquarkn->Momentum(q_p);
      Double_t rapiditypairquark    = (q_n + q_p).Rapidity();
      hRapCharmQuarkMethod2Pair  ->Fill(rapiditypairquark,eventw);
    }
  }
  if((Numberofelectron==1) && (Numberofpositron==1)) {
    hEtaPositronElectron->Fill(etapositron,etaelectron,eventw);
    Double_t deltaeta = etapositron-etaelectron;
    Double_t meaneta = (etapositron+etaelectron)/2.;
    Double_t deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(phipositron-phielectron));
    hEtaPositronElectronDelta->Fill(meaneta,deltaeta,eventw);
    hPhiPositronElectronDelta->Fill(deltaphi,eventw);
    if(highptDmesons == 1) hPhiPositronElectronDeltaonehighDmeson->Fill(deltaphi,eventw);
    if(highptDmesons == 2) hPhiPositronElectronDeltabothhighDmeson->Fill(deltaphi,eventw);
  }
  //printf("Number of electron found %d\n",Numberofelectron);
  //printf("Number of positron found %d\n",Numberofpositron);
  //printf("Number of high D mesons %d\n",highptDmesons);
  
  PostData(1, fOutputList);

  return;
}


//______________________________| Terminate
void AliAnalysisTaskCharm::Terminate(Option_t*)
{
  Info("Terminate","Start and end of Method");

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    Printf("ERROR: Output list not available");
    return;
  }
  return;
}

//______________________________| Definitions of histo etc.
void AliAnalysisTaskCharm::CreateHistos(){
  //
  // Init the histos
  //

  hNEvents    = new TH1F("hNEvents","",8,0.,8.);
  hNEventsW   = new TH1F("hNEventsW","",1000,0.,100.);

 
  Int_t   nBinRap   =  200;
  Float_t RapMin    = -10.;
  Float_t RapMax    =  10.;
  Int_t   nBinPtHad =  720;
  Float_t PtHadMin  =   0.;
  Float_t PtHadMax  =  36.;
  Int_t   nBinPtEle =  200;
  Float_t PtEleMin  =   0.;
  Float_t PtEleMax  =  10.;
  
  // Histograms for mesons including charm quark
  hQuarkMethod1  = new TH2F("hQuarkMethod1"  ,"rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hQuarkMethod2  = new TH2F("hQuarkMethod2"  ,"rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hCharm      = new TH2F("hCharm"  ,"rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hDpm        = new TH2F("hDpm"    ,"rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hD0         = new TH2F("hD0"     ,"rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hDs         = new TH2F("hDs"     ,"rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hDstar      = new TH2F("hDstar"  ,"rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hLambdac    = new TH2F("hLambdac","rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  //
  hDeltaPhi_D0     = new TH1F("hDeltaPhi_D0","deltaphi",20,0,TMath::Pi());
  hDeltaPhi_D0_LHCb  = new TH1F("hDeltaPhi_D0_LHCb","deltaphi",20,0,TMath::Pi());
  hDeltaPhi_D0_LHCb_e  = new TH1F("hDeltaPhi_D0_LHCb_e","deltaphi",20,0,TMath::Pi());


  // Histograms for rapidity distributions
  hRapElectronPt200         = new TH1F("hRapElectronPt200"       ,"y_{e},  pt>0.2 GeV/c"   ,100,-7.2,7.2); // for comparison to FONLL
  hRapElectronPt500         = new TH1F("hRapElectronPt500"       ,"y_{e},  pt>0.5 GeV/c"   ,100,-7.2,7.2); // for comparison to FONLL
  //
  hPtCharmQuarkMethod1     = new TH1F("hPtCharmQuarkMethod1"   ,"p_{T,c}"  ,nBinPtHad,PtHadMin,30.);
  hPtCharmQuarkMethod1CNM     = new TH1F("hPtCharmQuarkMethod1CNM"   ,"p_{T,c}" ,nBinPtHad,PtHadMin,30);
  hRapCharmQuarkMethod1     = new TH1F("hRapCharmQuarkMethod1"   ,"y_{c}"                  ,nBinRap,RapMin,RapMax);
  hRapCharmQuarkMethod2     = new TH1F("hRapCharmQuarkMethod2"   ,"y_{c}"                  ,nBinRap,RapMin,RapMax);
  hRapCharmQuarkMethod1Pair     = new TH1F("hRapCharmQuarkMethod1Pair"   ,"y_{c}"                  ,nBinRap,RapMin,RapMax);
  hRapCharmQuarkMethod2Pair     = new TH1F("hRapCharmQuarkMethod2Pair"   ,"y_{c}"                  ,nBinRap,RapMin,RapMax);
  hRapCharmQuark2DMethod1   = new TH2F("hRapCharmQuark2DMethod1" ,"y_{c} vs y_{cbar}"      ,nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hRapCharmQuark2DMethod2   = new TH2F("hRapCharmQuark2DMethod2" ,"y_{c} vs y_{cbar}"      ,nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hPhiCharmQuark1DMethod1   = new TH1F("hPhiCharmQuark1DMethod1" ,"#phi_{c} vs #phi_{cbar}"      ,200,-TMath::Pi(),TMath::Pi());
  hPhiCharmQuark1DMethod2   = new TH1F("hPhiCharmQuark1DMethod2" ,"#phi_{c} vs #phi_{cbar}"      ,200,-TMath::Pi(),TMath::Pi());
  hPtCharmQuark2DMethod1   = new TH2F("hPtCharmQuark2DMethod1" ,"p_{T,c} vs p_{T,cbar}"  ,nBinPtHad,PtHadMin,PtHadMax,nBinPtHad,PtHadMin,PtHadMax);
  hPtCharmQuark2DMethod2   = new TH2F("hPtCharmQuark2DMethod2" ,"p_{T,c} vs p_{T,cbar}"  ,nBinPtHad,PtHadMin,PtHadMax,nBinPtHad,PtHadMin,PtHadMax);
  //
  hRapCquarkDmeson  = new TH2F("hRapCquarkDmeson " ,"y_{c} vs y_{D}",nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hPtCquarkDmeson  = new TH2F("hPtCquarkDmeson " ,"p_{T,c} vs p_{T,D}",100,0.,20.,100,0.,20.);
  hRapDmesonElectron = new TH2F("hRapDmesonElectron ","y_{D} vs y_{e}",nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hPtDmesonElectron = new TH2F("hPtDmesonElectron ","p_{T,D} vs p_{T,e}",100,0.,20.,100,0.,20.);
  // hRapDmesonElectron is empty, so, wanted to fill it in another way
  hRapDMom2e         = new TH2F("hRapDMom2e " ,"y_{Dmom} vs y_{e}",nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hPtEtaElectron     = new TH2F("hPtEtaElectron","gen. electrons from charm;p_{T,e};#eta_{e}",nBinPtEle,PtEleMin,PtEleMax,nBinRap,RapMin,RapMax);
  hPtEtaElectronE     = new TH2F("hPtEtaElectronE","gen. electrons from charm;p_{T,e};#eta_{e}",nBinPtEle,PtEleMin,PtEleMax,nBinRap,RapMin,RapMax);
  hPtEtaElectronP     = new TH2F("hPtEtaElectronP","gen. electrons from charm;p_{T,e};#eta_{e}",nBinPtEle,PtEleMin,PtEleMax,nBinRap,RapMin,RapMax);
  hEtaPositronElectron  = new TH2F("hEtaPositronElectron" ,"#eta_{e^{+}} vs #eta_{e^{-}};#eta_{e^{+}};#eta_{e^{-}}",nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hEtaPositronElectronDelta  = new TH2F("hEtaPositronElectronDelta" ,"#Delta#eta_{e} vs Mean_{#eta_{e}};(#eta_{e^{+}}+#eta_{e^{-}})/2.;#eta_{e^{+}}-#eta_{e^{-}}",nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hPhiPositronElectronDelta  = new TH1F("hPhiPositronElectron" ,"#phi_{e^{+}}-#phi_{e^{-}} vs pair pt",100,0.,TMath::Pi());
  hPhiPositronElectronDeltaonehighDmeson  = new TH1F("hPhiPositronElectrononehighDmeson" ,"#phi_{e^{+}}-#phi_{e^{-}} vs pair pt",100,0.,TMath::Pi());
  hPhiPositronElectronDeltabothhighDmeson  = new TH1F("hPhiPositronElectronbothhighDmeson" ,"#phi_{e^{+}}-#phi_{e^{-}} vs pair pt",100,0.,TMath::Pi());

  Int_t   nBinMass   = 800;
  Float_t MassMin    = 0.;
  Float_t MassMax    = 8.;
  Int_t   nBinPairPt = 400;
  Float_t PairPtMin  = 0.;
  Float_t PairPtMax  = 10.;

  Int_t   nBinMassl   = 35;
  Float_t MassMinl    = 0.;
  Float_t MassMaxl    = 3.5;
  Int_t   nBinPairPtl = 80;
  Float_t PairPtMinl  = 0.;
  Float_t PairPtMaxl  = 8.;

  // Histograms for Pt spectra : c-->e , cBar->e
  hPte_eta08       = new TH1F("hPte_eta08" ,"e,   |eta|<0.8;p_{T}" ,nBinPtEle,PtEleMin,PtEleMax);
  hPteP_eta08      = new TH1F("hPteP_eta08","e+,  |eta|<0.8;p_{T}" ,nBinPtEle,PtEleMin,PtEleMax);
  hPteM_eta08      = new TH1F("hPteM_eta08","e-,  |eta|<0.8;p_{T}" ,nBinPtEle,PtEleMin,PtEleMax);
  hPte_y08         = new TH1F("hPte_y08 "  ,"e,   |y|<0.8;p_{T}"   ,nBinPtEle,PtEleMin,PtEleMax);
  hPteP_y08        = new TH1F("hPteP_y08 " ,"e+,  |y|<0.8;p_{T}"   ,nBinPtEle,PtEleMin,PtEleMax);
  hPteM_y08        = new TH1F("hPteM_y08 " ,"e-,  |y|<0.8;p_{T}"   ,nBinPtEle,PtEleMin,PtEleMax);

  // Histograms for invariant mass spectra (ULS,LS): c-->e , cBar->e
  hMee_ULS_simulated    = new TH1F("hMee_ULS_simulated"  ,"e+e-;m_{ee}"                                 ,nBinMass,MassMin,MassMax);
  hMee_LS_simulated     = new TH1F("hMee_LS_simulated"   ,"e+e+ & e-e-;m_{ee}"                          ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta08        = new TH1F("hMee_ULS_eta08"      ,"e+e-,        |eta|<0.8;m_{ee}"               ,nBinMass,MassMin,MassMax);
  hMee_LS_eta08         = new TH1F("hMee_LS_eta08"       ,"e+e+ & e-e-, |eta|<0.8;m_{ee}"               ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta08_pt200  = new TH1F("hMee_ULS_eta08_pt200","e+e-,        |eta|<0.8 , pt>0.2 GeV/c;m_{ee}",nBinMass,MassMin,MassMax);
  hMee_LS_eta08_pt200   = new TH1F("hMee_LS_eta08_pt200" ,"e+e+ & e-e-, |eta|<0.8 , pt>0.2 GeV/c;m_{ee}",nBinMass,MassMin,MassMax);
  // new histograms for phenix acceptance
  hMeePtee_ULS_eta05            = new TH2F("hMeePtee_ULS_eta05"            ,"e+e-,        |eta|<0.5;m_{ee}"                ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta05             = new TH2F("hMeePtee_LS_eta05"             ,"e+e+ & e-e-, |eta|<0.5;m_{ee}"                ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta035           = new TH2F("hMeePtee_ULS_eta035"           ,"e+e-,        |eta|<0.35;m_{ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta035            = new TH2F("hMeePtee_LS_eta035"            ,"e+e+ & e-e-, |eta|<0.35;m_{ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hPhiee_ULS_eta08_pt200         = new TH3F("hPhiee_ULS_eta08_pt200"           ,"",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_LS_eta08_pt200          = new TH3F("hPhiee_LS_eta08_pt200"            ,"",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_ULS_eta08_pt400         = new TH3F("hPhiee_ULS_eta08_pt400"           ,"",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_LS_eta08_pt400          = new TH3F("hPhiee_LS_eta08_pt400"            ,"",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hMeePtee_ULS_eta035_phenixacc = new TH2F("hMeePtee_ULS_eta035_phenixacc" ,"e+e-,        |eta|<0.35 , pt>0.2 GeV/c;m_{ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta035_phenixacc  = new TH2F("hMeePtee_LS_eta035_phenixacc"  ,"e+e- & e-e-, |eta|<0.35 , pt>0.2 GeV/c;m_{ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hPhiee_ULS_eta08_highoneDmeson_pt200  = new TH3F("hPhiee_ULS_eta08_highoneDmeson_pt200","",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_LS_eta08_highoneDmeson_pt200   = new TH3F("hPhiee_LS_eta08_highoneDmeson_pt200" ,"",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_ULS_eta08_highbothDmeson_pt200  = new TH3F("hPhiee_ULS_eta08_highbothDmeson_pt200","",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_LS_eta08_highbothDmeson_pt200   = new TH3F("hPhiee_LS_eta08_highbothDmeson_pt200" ,"",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_ULS_eta08_highoneDmeson_pt400  = new TH3F("hPhiee_ULS_eta08_highoneDmeson_pt400","",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_LS_eta08_highoneDmeson_pt400   = new TH3F("hPhiee_LS_eta08_highoneDmeson_pt400" ,"",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_ULS_eta08_highbothDmeson_pt400  = new TH3F("hPhiee_ULS_eta08_highbothDmeson_pt400","",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hPhiee_LS_eta08_highbothDmeson_pt400   = new TH3F("hPhiee_LS_eta08_highbothDmeson_pt400" ,"",100,0.,TMath::Pi(),nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  //------

  // Histograms for invariant mass vs pair pt (ULS,LS): c-->e , cBar->e
  hMeePtee_ULS_eta08       = new TH2F("hMeePtee_ULS_eta08"      ,"e+e-,        |eta|<0.8;m_{ee};p_{T,ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08        = new TH2F("hMeePtee_LS_eta08"       ,"e+e- & e-e-, |eta|<0.8;m_{ee};p_{T,ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt200 = new TH2F("hMeePtee_ULS_eta08_pt200","e+e-,        |eta|<0.8 , pt>0.2 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt200  = new TH2F("hMeePtee_LS_eta08_pt200" ,"e+e- & e-e-, |eta|<0.8 , pt>0.2 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta_pt = new TH2F("hMeePtee_ULS_eta_pt","e+e- ;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta_pt  = new TH2F("hMeePtee_LS_eta_pt" ,"e+e- & e-e- ;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt300 = new TH2F("hMeePtee_ULS_eta08_pt300","e+e-,        |eta|<0.8 , pt>0.3 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt300  = new TH2F("hMeePtee_LS_eta08_pt300" ,"e+e- & e-e-, |eta|<0.8 , pt>0.3 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt400 = new TH2F("hMeePtee_ULS_eta08_pt400","e+e-,        |eta|<0.8 , pt>0.4 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt400  = new TH2F("hMeePtee_LS_eta08_pt400" ,"e+e- & e-e-, |eta|<0.8 , pt>0.4 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt200_opAngle50 = new TH2F("hMeePtee_ULS_eta08_pt200_opAngle50","e+e-,        |eta|<0.8 , pt>0.2 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt200_opAngle50  = new TH2F("hMeePtee_LS_eta08_pt200_opAngle50" ,"e+e- & e-e-, |eta|<0.8 , pt>0.2 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt300_opAngle50 = new TH2F("hMeePtee_ULS_eta08_pt300_opAngle50","e+e-,        |eta|<0.8 , pt>0.3 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt300_opAngle50  = new TH2F("hMeePtee_LS_eta08_pt300_opAngle50" ,"e+e- & e-e-, |eta|<0.8 , pt>0.3 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt400_opAngle50 = new TH2F("hMeePtee_ULS_eta08_pt400_opAngle50","e+e-,        |eta|<0.8 , pt>0.4 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt400_opAngle50  = new TH2F("hMeePtee_LS_eta08_pt400_opAngle50" ,"e+e- & e-e-, |eta|<0.8 , pt>0.4 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  // opening angle
  hMeeOpAngle_ULS_eta08_pt200 = new TH2F("hMeeOpAngle_ULS_eta08_pt200","e+e-,        |eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax,128,0,3.2);
  hMeeOpAngle_LS_eta08_pt200  = new TH2F("hMeeOpAngle_LS_eta08_pt200" ,"e+e-,        |eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax,128,0,3.2);

  // Mean pt of mother D mesons
  hMotherPt_ULS_eta08_pt200          = new TH3F("hMotherPt_ULS_eta08_pt200","",nBinPairPtl,PairPtMinl,PairPtMaxl,nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hMotherPt_LS_eta08_pt200          = new TH3F("hMotherPt_LS_eta08_pt200","",nBinPairPtl,PairPtMinl,PairPtMaxl,nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hMotherPt_ULS_eta08_pt400          = new TH3F("hMotherPt_ULS_eta08_pt400","",nBinPairPtl,PairPtMinl,PairPtMaxl,nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);
  hMotherPt_LS_eta08_pt400          = new TH3F("hMotherPt_LS_eta08_pt400","",nBinPairPtl,PairPtMinl,PairPtMaxl,nBinMassl,MassMinl,MassMaxl,nBinPairPtl,PairPtMinl,PairPtMaxl);

  // Monitor CNM weighting procedure
  hweightcnmD      = new TH2F("hweightcnmD"  ,"w vs pt;w;p_{T}"   ,50,-2.,2.,nBinPtHad,PtHadMin,PtHadMax); 
  hweightcnmHFE    = new TH2F("hweightcnmHFE"  ,"w vs pt;w;p_{T}"   ,50,-2.,2.,nBinPtHad,PtHadMin,PtHadMax); 
  hweightcnmee      = new TH2F("hweightcnmee"  ,"w vs pt;w;p_{T}"   ,50,-2.,2.,nBinPtHad,PtHadMin,PtHadMax); 
     

  hNEventsW->Sumw2();
  hNEvents->Sumw2();
  hQuarkMethod1->Sumw2();
  hQuarkMethod2->Sumw2();
  hCharm->Sumw2();
  hDpm->Sumw2();
  hD0->Sumw2();
  hDs->Sumw2();
  hDstar->Sumw2();
  hLambdac->Sumw2();
  hDeltaPhi_D0->Sumw2();
  hDeltaPhi_D0_LHCb->Sumw2();
  hDeltaPhi_D0_LHCb_e->Sumw2();
  hRapElectronPt200->Sumw2();
  hRapElectronPt500->Sumw2();
  hPtCharmQuarkMethod1->Sumw2();
  hPtCharmQuarkMethod1CNM->Sumw2();
  hRapCharmQuarkMethod1->Sumw2();
  hRapCharmQuarkMethod2->Sumw2();
  hRapCharmQuarkMethod1Pair->Sumw2();
  hRapCharmQuarkMethod2Pair->Sumw2();
  hRapCharmQuark2DMethod1->Sumw2();
  hRapCharmQuark2DMethod2->Sumw2();
  hPhiCharmQuark1DMethod1->Sumw2();
  hPhiCharmQuark1DMethod2->Sumw2();
  hPtCharmQuark2DMethod1->Sumw2();
  hPtCharmQuark2DMethod2->Sumw2();
  hRapCquarkDmeson->Sumw2();
  hPtCquarkDmeson->Sumw2();
  hRapDmesonElectron->Sumw2();
  hPtDmesonElectron->Sumw2();
  hRapDMom2e->Sumw2();
  hPtEtaElectron->Sumw2();
  hPtEtaElectronE->Sumw2();
  hPtEtaElectronP->Sumw2();
  hEtaPositronElectron->Sumw2();
  hEtaPositronElectronDelta->Sumw2();
  hPte_eta08->Sumw2();
  hPteP_eta08->Sumw2();
  hPteM_eta08->Sumw2();
  hPte_y08->Sumw2();
  hPteP_y08->Sumw2();
  hPteM_y08->Sumw2();
  hPhiPositronElectronDelta->Sumw2();
  hPhiPositronElectronDeltaonehighDmeson->Sumw2();
  hPhiPositronElectronDeltabothhighDmeson->Sumw2();
  hMee_ULS_simulated->Sumw2();
  hMee_LS_simulated->Sumw2();
  hMeePtee_ULS_eta05->Sumw2();
  hMeePtee_LS_eta05->Sumw2();
  hMeePtee_ULS_eta035->Sumw2();
  hMeePtee_LS_eta035->Sumw2();
  hPhiee_ULS_eta08_pt200->Sumw2();
  hPhiee_LS_eta08_pt200->Sumw2();
  hPhiee_ULS_eta08_pt400->Sumw2();
  hPhiee_LS_eta08_pt400->Sumw2();
  hMeePtee_ULS_eta035_phenixacc->Sumw2();
  hMeePtee_LS_eta035_phenixacc->Sumw2();
  hPhiee_ULS_eta08_highoneDmeson_pt200->Sumw2();
  hPhiee_LS_eta08_highoneDmeson_pt200->Sumw2();
  hPhiee_ULS_eta08_highbothDmeson_pt200->Sumw2();
  hPhiee_LS_eta08_highbothDmeson_pt200->Sumw2();
  hPhiee_ULS_eta08_highoneDmeson_pt400->Sumw2();
  hPhiee_LS_eta08_highoneDmeson_pt400->Sumw2();
  hPhiee_ULS_eta08_highbothDmeson_pt400->Sumw2();
  hPhiee_LS_eta08_highbothDmeson_pt400->Sumw2();
  hMee_ULS_eta08->Sumw2();
  hMee_LS_eta08->Sumw2();
  hMee_ULS_eta08_pt200->Sumw2();
  hMee_LS_eta08_pt200->Sumw2();
  hMeePtee_ULS_eta08->Sumw2();
  hMeePtee_LS_eta08->Sumw2();
  hMeePtee_ULS_eta08_pt200->Sumw2();
  hMeePtee_LS_eta08_pt200->Sumw2();
  hMeePtee_ULS_eta_pt->Sumw2();
  hMeePtee_LS_eta_pt->Sumw2();
  hMeePtee_ULS_eta08_pt300->Sumw2();
  hMeePtee_LS_eta08_pt300->Sumw2();
  hMeePtee_ULS_eta08_pt400->Sumw2();
  hMeePtee_LS_eta08_pt400->Sumw2();
  hMeePtee_ULS_eta08_pt200_opAngle50->Sumw2();
  hMeePtee_LS_eta08_pt200_opAngle50->Sumw2();
  hMeePtee_ULS_eta08_pt300_opAngle50->Sumw2();
  hMeePtee_LS_eta08_pt300_opAngle50->Sumw2();
  hMeePtee_ULS_eta08_pt400_opAngle50->Sumw2();
  hMeePtee_LS_eta08_pt400_opAngle50->Sumw2();
  hMeeOpAngle_ULS_eta08_pt200->Sumw2();
  hMeeOpAngle_LS_eta08_pt200->Sumw2();
  hMotherPt_ULS_eta08_pt200->Sumw2();
  hMotherPt_LS_eta08_pt200->Sumw2();
  hMotherPt_ULS_eta08_pt400->Sumw2();
  hMotherPt_LS_eta08_pt400->Sumw2();
  hweightcnmD->Sumw2();
  hweightcnmHFE->Sumw2();
  hweightcnmee->Sumw2();


  fOutputList->Add(hNEventsW);
  fOutputList->Add(hNEvents);
  fOutputList->Add(hQuarkMethod1);
  fOutputList->Add(hQuarkMethod2);
  fOutputList->Add(hCharm);
  fOutputList->Add(hDpm);
  fOutputList->Add(hD0);
  fOutputList->Add(hDs);
  fOutputList->Add(hDstar);
  fOutputList->Add(hLambdac);
  fOutputList->Add(hDeltaPhi_D0);
  fOutputList->Add(hDeltaPhi_D0_LHCb);
  fOutputList->Add(hDeltaPhi_D0_LHCb_e);
  fOutputList->Add(hRapElectronPt200);
  fOutputList->Add(hRapElectronPt500);
  fOutputList->Add(hPtCharmQuarkMethod1);
  fOutputList->Add(hPtCharmQuarkMethod1CNM);
  fOutputList->Add(hRapCharmQuarkMethod1);
  fOutputList->Add(hRapCharmQuarkMethod2);
  fOutputList->Add(hRapCharmQuarkMethod1Pair);
  fOutputList->Add(hRapCharmQuarkMethod2Pair);
  fOutputList->Add(hRapCharmQuark2DMethod1);
  fOutputList->Add(hRapCharmQuark2DMethod2);
  fOutputList->Add(hPhiCharmQuark1DMethod1);
  fOutputList->Add(hPhiCharmQuark1DMethod2);
  fOutputList->Add(hPtCharmQuark2DMethod1);
  fOutputList->Add(hPtCharmQuark2DMethod2);
  fOutputList->Add(hRapCquarkDmeson);
  fOutputList->Add(hPtCquarkDmeson);
  fOutputList->Add(hRapDmesonElectron);
  fOutputList->Add(hPtDmesonElectron);
  fOutputList->Add(hRapDMom2e);
  fOutputList->Add(hPtEtaElectron);
  fOutputList->Add(hPtEtaElectronE);
  fOutputList->Add(hPtEtaElectronP);
  fOutputList->Add(hEtaPositronElectron);
  fOutputList->Add(hEtaPositronElectronDelta);
  fOutputList->Add(hPhiPositronElectronDelta);
  fOutputList->Add(hPhiPositronElectronDeltaonehighDmeson);
  fOutputList->Add(hPhiPositronElectronDeltabothhighDmeson);
  fOutputList->Add(hPte_eta08);
  fOutputList->Add(hPteP_eta08);
  fOutputList->Add(hPteM_eta08);
  fOutputList->Add(hPte_y08);
  fOutputList->Add(hPteP_y08);
  fOutputList->Add(hPteM_y08);
  fOutputList->Add(hMee_ULS_simulated);
  fOutputList->Add(hMee_LS_simulated);
  fOutputList->Add(hMeePtee_ULS_eta05);
  fOutputList->Add(hMeePtee_LS_eta05);
  fOutputList->Add(hMeePtee_ULS_eta035);
  fOutputList->Add(hMeePtee_LS_eta035);
  fOutputList->Add(hPhiee_ULS_eta08_pt200);
  fOutputList->Add(hPhiee_LS_eta08_pt200);
  fOutputList->Add(hPhiee_ULS_eta08_pt400);
  fOutputList->Add(hPhiee_LS_eta08_pt400);
  fOutputList->Add(hMeePtee_ULS_eta035_phenixacc);
  fOutputList->Add(hMeePtee_LS_eta035_phenixacc);
  fOutputList->Add(hPhiee_ULS_eta08_highoneDmeson_pt200);
  fOutputList->Add(hPhiee_LS_eta08_highoneDmeson_pt200);
  fOutputList->Add(hPhiee_ULS_eta08_highbothDmeson_pt200);
  fOutputList->Add(hPhiee_LS_eta08_highbothDmeson_pt200);
  fOutputList->Add(hPhiee_ULS_eta08_highoneDmeson_pt400);
  fOutputList->Add(hPhiee_LS_eta08_highoneDmeson_pt400);
  fOutputList->Add(hPhiee_ULS_eta08_highbothDmeson_pt400);
  fOutputList->Add(hPhiee_LS_eta08_highbothDmeson_pt400);
  fOutputList->Add(hMee_ULS_eta08);
  fOutputList->Add(hMee_LS_eta08);
  fOutputList->Add(hMee_ULS_eta08_pt200);
  fOutputList->Add(hMee_LS_eta08_pt200);
  fOutputList->Add(hMeePtee_ULS_eta08);
  fOutputList->Add(hMeePtee_LS_eta08);
  fOutputList->Add(hMeePtee_ULS_eta08_pt200);
  fOutputList->Add(hMeePtee_LS_eta08_pt200);
  fOutputList->Add(hMeePtee_ULS_eta_pt);
  fOutputList->Add(hMeePtee_LS_eta_pt);
  fOutputList->Add(hMeePtee_ULS_eta08_pt400);
  fOutputList->Add(hMeePtee_LS_eta08_pt400);
  fOutputList->Add(hMotherPt_ULS_eta08_pt200);
  fOutputList->Add(hMotherPt_LS_eta08_pt200);
  fOutputList->Add(hMotherPt_ULS_eta08_pt400);
  fOutputList->Add(hMotherPt_LS_eta08_pt400);
  fOutputList->Add(hweightcnmD);
  fOutputList->Add(hweightcnmHFE);
  fOutputList->Add(hweightcnmee);
 

}


Double_t AliAnalysisTaskCharm::pt_cut200(Double_t pT) {
  //Double_t weight=1.0; // maybe a bit quicker without memory allocation.
  if (pT<0.2) return 0.0;
  return 1.0;
}

Double_t AliAnalysisTaskCharm::pt_cut300(Double_t pT) {
  if (pT<0.3) return 0.0;
  return 1.0;
}

Double_t AliAnalysisTaskCharm::pt_cut400(Double_t pT) {
  if (pT<0.4) return 0.0;
  return 1.0;
}

Double_t AliAnalysisTaskCharm::pt_cutHigh(Double_t pT) {
  if (pT>fPtCutHigh) return 0.0;
  return 1.0;
}

Double_t AliAnalysisTaskCharm::pt_cutLow(Double_t pT) {
  if (pT<fPtCutLow) return 0.0;
  return 1.0;
}

Double_t AliAnalysisTaskCharm::scale_RAA(Double_t pT) {
  //return 8.46749e-01 + (-1.30301e-01) * pT; // fit 0-10%
  //return 9.03546e-01 + (-1.09215e-01) * pT; // fit 20-40%
  return 0.885416 + (-0.114357) * pT; // fits for 10-20% + 20-40% combined
}

Double_t AliAnalysisTaskCharm::scale_CNM(Double_t pT) {
  if(!fgraphCNM) return 0.;
  return fgraphCNM->Eval(pT); // 
}

Bool_t AliAnalysisTaskCharm::Inphenixacc(Double_t phi,Double_t pt, Int_t pdg){
  //
  // Check the phenix acceptance in phi and pt cut at 0.2
  //
  Double_t kdc = 0.206;
  Double_t krich= 0.309;
  // Acceptance defined in [0-2Pi] and not [-pi/2,3pi/2] like for PHENIX
  Double_t phicuta1 = -3./16.*TMath::Pi() + 0.5*TMath::Pi();
  Double_t phicuta2 = 5./16.*TMath::Pi() + 0.5*TMath::Pi();
  Double_t phicutb1 = 11./16.*TMath::Pi() + 0.5*TMath::Pi();
  Double_t phicutb2 = 19./16.*TMath::Pi() + 0.5*TMath::Pi();

  Double_t charge = 1.;
  if(pdg>0) charge = -1.;

  if(pt<0.2) return kFALSE;

  // Check that phi is in [0-2Pi]
  Double_t phia = TVector2::Phi_0_2pi(phi + charge*kdc/pt);
  Double_t phib = TVector2::Phi_0_2pi(phi + charge*krich/pt);

  // Two arms
  Bool_t firstarm = kTRUE;
  if((phia<phicuta1) || (phia>phicuta2)) firstarm = kFALSE;
  Bool_t secondarm = kTRUE;
  if((phib<phicutb1) || (phib>phicutb2)) secondarm = kFALSE;

  // return
  if(firstarm || secondarm) return kTRUE;
  else return kFALSE;

}
