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
//using namespace std;
// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TChain.h"
#include "TVector2.h"
#include <TMath.h>
//#include <TDatabasePDG.h>
#include "TParticle.h"
#include "TLorentzVector.h"
// AliRoot includes
#include "AliAnalysisTaskBeauty.h" // this task
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliVParticle.h"
#include "AliInputEventHandler.h"

ClassImp(AliAnalysisTaskBeauty)

//______________________________| Default Constructor
AliAnalysisTaskBeauty::AliAnalysisTaskBeauty():
AliAnalysisTaskSE(),
AliCocktailSmearing(),
fMcEvent(0x0),
fMcHandler(0x0),
fNbEvents(0),
fProcessType(0),
fEventWeight(0),
fPtCutHigh(1e6),
fScaleByRAA(kFALSE),
hNEvents(0),
hNEventsW(0),
hQuarkMethod1(0),                     
hQuarkMethod2(0),                     
hBeauty(0),
hBpm(0),
hB0(0),
hBs(0),
hLambdab(0),
hRapElectron_be_Pt200(0),
hRapElectron_be_Pt500(0),
hRapElectron_bce_Pt200(0),
hRapElectron_bce_Pt500(0),
hRapBeautyQuarkMethod1(0),
hRapBeautyQuarkMethod2(0),
hRapBeautyQuark2DMethod1(0),
hRapBeautyQuark2DMethod2(0),
hRapBMom2e(0),
hRapBMom2c2e(0),
hPtEtaElectron_be(0),
hPtEtaElectron_bce(0),
hPte_eta08_be(0),
hPteP_eta08_be(0),
hPteM_eta08_be(0),
hPte_y08_be(0),
hPteP_y08_be(0),
hPteM_y08_be(0),
hPte_eta08_bce(0),
hPteP_eta08_bce(0),
hPteM_eta08_bce(0),
hPte_y08_bce(0),
hPteP_y08_bce(0),
hPteM_y08_bce(0),
hMee_ULS_simulated(0),
hMee_LS_simulated(0),
hMee_ULS_eta05(0),
hMee_LS_eta05(0),
hMee_ULS_eta08(0),
hMee_LS_eta08(0),
hMee_ULS_eta035(0),
hMee_LS_eta035(0),
hMee_ULS_eta08_pt200(0),
hMee_LS_eta08_pt200(0),
hMee_ULS_eta08_pt400(0),
hMee_LS_eta08_pt400(0),
hMee_ULS_eta035_phenixacc(0),
hMee_LS_eta035_phenixacc(0),
hMeePtee_ULS_eta08(0),
hMeePtee_LS_eta08(0),
hMeePtee_ULS_eta08_pt200(0),
hMeePtee_LS_eta08_pt200(0),
hMeePtee_ULS_eta08_pt400(0),
hMeePtee_LS_eta08_pt400(0),
hMeePtee_ULS_eta08_pt200_opAngle50(0),
hMeePtee_LS_eta08_pt200_opAngle50(0),
hMeePtee_ULS_eta08_pt300_opAngle50(0),
hMeePtee_LS_eta08_pt300_opAngle50(0),
hMeePtee_ULS_eta08_pt400_opAngle50(0),
hMeePtee_LS_eta08_pt400_opAngle50(0),
hMee_ULS_simulated_be(0),
hMee_LS_simulated_be(0),
hMee_ULS_eta05_be(0),
hMee_LS_eta05_be(0),
hMee_ULS_eta08_be(0),
hMee_LS_eta08_be(0),
hMee_ULS_eta035_be(0),
hMee_LS_eta035_be(0),
hMee_ULS_eta08_pt200_be(0),
hMee_LS_eta08_pt200_be(0),
hMee_ULS_eta08_pt400_be(0),
hMee_LS_eta08_pt400_be(0),
hMee_ULS_eta035_phenixacc_be(0),
hMee_LS_eta035_phenixacc_be(0),
hMeePtee_ULS_eta08_be(0),
hMeePtee_LS_eta08_be(0),
hMeePtee_ULS_eta08_pt200_be(0),
hMeePtee_LS_eta08_pt200_be(0),
hMeePtee_ULS_eta08_pt400_be(0),
hMeePtee_LS_eta08_pt400_be(0),
hMee_ULS_simulated_bce(0),
hMee_LS_simulated_bce(0),
hMee_ULS_eta05_bce(0),
hMee_LS_eta05_bce(0),
hMee_ULS_eta08_bce(0),
hMee_LS_eta08_bce(0),
hMee_ULS_eta035_bce(0),
hMee_LS_eta035_bce(0),
hMee_ULS_eta08_pt200_bce(0),
hMee_LS_eta08_pt200_bce(0),
hMee_ULS_eta08_pt400_bce(0),
hMee_LS_eta08_pt400_bce(0),
hMee_ULS_eta035_phenixacc_bce(0),
hMee_LS_eta035_phenixacc_bce(0),
hMeePtee_ULS_eta08_bce(0),
hMeePtee_LS_eta08_bce(0),
hMeePtee_ULS_eta08_pt200_bce(0),
hMeePtee_LS_eta08_pt200_bce(0),
hMeePtee_ULS_eta08_pt400_bce(0),
hMeePtee_LS_eta08_pt400_bce(0),
hMeeOpAngle_ULS_eta08_pt200(0),
hMeeOpAngle_LS_eta08_pt200(0),
fOutputList(0)
{

}

//______________________________| Specific Constructor
AliAnalysisTaskBeauty::AliAnalysisTaskBeauty(const Char_t* name) :
AliAnalysisTaskSE(name),
AliCocktailSmearing(),
fMcEvent(0x0),
fMcHandler(0x0),
fNbEvents(0),
fProcessType(0),
fEventWeight(0),
fPtCutHigh(1e6),
fScaleByRAA(kFALSE),
hNEvents(0),
hNEventsW(0),
hQuarkMethod1(0),                     
hQuarkMethod2(0),
hBeauty(0),
hBpm(0),
hB0(0),
hBs(0),
hLambdab(0),
hRapElectron_be_Pt200(0),
hRapElectron_be_Pt500(0),
hRapElectron_bce_Pt200(0),
hRapElectron_bce_Pt500(0),
hRapBeautyQuarkMethod1(0),
hRapBeautyQuarkMethod2(0),
hRapBeautyQuark2DMethod1(0),
hRapBeautyQuark2DMethod2(0),
hRapBMom2e(0),
hRapBMom2c2e(0),
hPtEtaElectron_be(0),
hPtEtaElectron_bce(0),
hPte_eta08_be(0),
hPteP_eta08_be(0),
hPteM_eta08_be(0),
hPte_y08_be(0),
hPteP_y08_be(0),
hPteM_y08_be(0),
hPte_eta08_bce(0),
hPteP_eta08_bce(0),
hPteM_eta08_bce(0),
hPte_y08_bce(0),
hPteP_y08_bce(0),
hPteM_y08_bce(0),
hMee_ULS_simulated(0),
hMee_LS_simulated(0),
hMee_ULS_eta05(0),
hMee_LS_eta05(0),
hMee_ULS_eta08(0),
hMee_LS_eta08(0),
hMee_ULS_eta035(0),
hMee_LS_eta035(0),
hMee_ULS_eta08_pt200(0),
hMee_LS_eta08_pt200(0),
hMee_ULS_eta08_pt400(0),
hMee_LS_eta08_pt400(0),
hMee_ULS_eta035_phenixacc(0),
hMee_LS_eta035_phenixacc(0),
hMeePtee_ULS_eta08(0),
hMeePtee_LS_eta08(0),
hMeePtee_ULS_eta08_pt200(0),
hMeePtee_LS_eta08_pt200(0),
hMeePtee_ULS_eta08_pt400(0),
hMeePtee_LS_eta08_pt400(0),
hMeePtee_ULS_eta08_pt200_opAngle50(0),
hMeePtee_LS_eta08_pt200_opAngle50(0),
hMeePtee_ULS_eta08_pt300_opAngle50(0),
hMeePtee_LS_eta08_pt300_opAngle50(0),
hMeePtee_ULS_eta08_pt400_opAngle50(0),
hMeePtee_LS_eta08_pt400_opAngle50(0),
hMee_ULS_simulated_be(0),
hMee_LS_simulated_be(0),
hMee_ULS_eta05_be(0),
hMee_LS_eta05_be(0),
hMee_ULS_eta08_be(0),
hMee_LS_eta08_be(0),
hMee_ULS_eta035_be(0),
hMee_LS_eta035_be(0),
hMee_ULS_eta08_pt200_be(0),
hMee_LS_eta08_pt200_be(0),
hMee_ULS_eta08_pt400_be(0),
hMee_LS_eta08_pt400_be(0),
hMee_ULS_eta035_phenixacc_be(0),
hMee_LS_eta035_phenixacc_be(0),
hMeePtee_ULS_eta08_be(0),
hMeePtee_LS_eta08_be(0),
hMeePtee_ULS_eta08_pt200_be(0),
hMeePtee_LS_eta08_pt200_be(0),
hMeePtee_ULS_eta08_pt400_be(0),
hMeePtee_LS_eta08_pt400_be(0),
hMee_ULS_simulated_bce(0),
hMee_LS_simulated_bce(0),
hMee_ULS_eta05_bce(0),
hMee_LS_eta05_bce(0),
hMee_ULS_eta08_bce(0),
hMee_LS_eta08_bce(0),
hMee_ULS_eta035_bce(0),
hMee_LS_eta035_bce(0),
hMee_ULS_eta08_pt200_bce(0),
hMee_LS_eta08_pt200_bce(0),
hMee_ULS_eta08_pt400_bce(0),
hMee_LS_eta08_pt400_bce(0),
hMee_ULS_eta035_phenixacc_bce(0),
hMee_LS_eta035_phenixacc_bce(0),
hMeePtee_ULS_eta08_bce(0),
hMeePtee_LS_eta08_bce(0),
hMeePtee_ULS_eta08_pt200_bce(0),
hMeePtee_LS_eta08_pt200_bce(0),
hMeePtee_ULS_eta08_pt400_bce(0),
hMeePtee_LS_eta08_pt400_bce(0),
hMeeOpAngle_ULS_eta08_pt200(0),
hMeeOpAngle_LS_eta08_pt200(0),
fOutputList(0)
{
  Info("AliAnalysisTaskBeauty","Calling Constructor");
  // Output slot #1 writes into a TList container (nevents histogram)
  DefineInput(0, TChain::Class());
  DefineOutput(1,TList::Class());
}

//______________________________| Destructor
AliAnalysisTaskBeauty::~AliAnalysisTaskBeauty()
{
  // Destructor
  Info("~AliAnalysisTaskBeauty","Calling Destructor");
  if (fOutputList) delete fOutputList;
}

//______________________________| User Output
void AliAnalysisTaskBeauty::UserCreateOutputObjects()
{
  Info("AliAnalysisTaskBeauty","CreateOutputObjects of task %s", GetName());
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  printf("\n\n========================================\n  Configuration of task: \n========================================\n\n");
  printf("  process type:        %d\n", fProcessType);
  printf("  high-pt cut:         %f\n", fPtCutHigh);
  printf("  use R_AA scaling:    %s\n", fScaleByRAA?"YES":"NO");
  printf("  use Event weight:    %s\n", fEventWeight?"YES":"NO");
  std::cout << std::endl;

  AliCocktailSmearing::Print();

  CreateHistos();
  fNbEvents = 0;

  PostData(1, fOutputList);
}

//______________________________| Init
void AliAnalysisTaskBeauty::Init()
{
  if(fDebug > 1) printf("AliAnalysisTaskBeauty::Init() \n");
  fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}


//______________________________| User Exec
void AliAnalysisTaskBeauty::UserExec(Option_t *)
{

  if(fDebug > 1) printf("AliAnalysisTaskBeauty::UserExec \n");
  Init();

  if(fMcHandler){
    fMcEvent = fMcHandler->MCEvent();
  }else{
    if(fDebug > 1) printf("AliAnalysisTaskBeauty::Handler() fMcHandler=NULL\n");
    return;
  }
  if(!fMcEvent){
    if(fDebug > 1) printf("AliAnalysisTaskBeauty::UserExec()   fMcEvent=NULL \n");
    return;
  }


  int nparticles = fMcEvent->GetNumberOfTracks();

  double EtaCut1    =1.;
  double EtaCut     =0.8;
  //double EtaCutTight=0.5;
  //double EtaCutPhenix=0.35;

  double opAngle50   =0.05; // 50 mrad

  // one can intialize it out of the event loop (if one reset this in the event loop after each event) :
  bool IsBmeson      = kFALSE; //kTRUE if all Bmeson events, in |y|<1
  bool IsBeauty      = kFALSE; //kTRUE if all beauty events, in |y|<1
  bool IsElectronB   = kFALSE; //kTRUE if all events for electrons coming from beauty quark, in |y|<1
  bool IsElectronD   = kFALSE;
  //
  // these booleans are per track so they have to be resetted in the single and pair loops!
  bool is_beauty         = kFALSE;
  bool is_charm2e        = kFALSE;
  bool is_beauty2charm   = kFALSE;
  bool i_is_beauty       = kFALSE;
  bool j_is_beauty       = kFALSE;
  bool i_is_charm2e      = kFALSE;
  bool j_is_charm2e      = kFALSE;
  bool i_is_beauty2charm = kFALSE;
  bool j_is_beauty2charm = kFALSE;
  bool IsCharmQ = kFALSE;

  //// Debug at quark level: two methods
  // Search b quarks which fragment
  Int_t Method1nbquarkintheeventp = 0;
  Int_t Method1nbquarkintheeventn = 0;
  Double_t Method1yquarkp = -1000.;
  Double_t Method1yquarkn = -1000.;
  // Search first two b quarks in the tree
  Int_t Method2nbquarkintheeventp = 0;
  Int_t Method2nbquarkintheeventn = 0;
  Double_t Method2yquarkp = -1000.;
  Double_t Method2yquarkn = -1000.;

  //============================ single loop =================================

  Float_t crosssection = 1.;
  // Find process
  AliGenEventHeader* header = fMcEvent->GenEventHeader();
  if(header) {
    //printf("Find the header %s\n",header->ClassName());
    AliGenPythiaEventHeader *pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(header);
    if(pythiaHeader){
      //printf("Pythia header found\n");
      Int_t process = pythiaHeader->ProcessType();
      //printf("Process type is %d\n",process);
      if((fProcessType) && (process!=fProcessType)) return;
      crosssection = pythiaHeader->GetXsection();
      //printf("Get cross section %f\n",crosssection);
    }
  }
  //else printf("Did not find the header\n");

  Double_t eventw = 1.;
  if(fEventWeight) eventw = crosssection;

  //printf("Event number %d\n",fNbEvents);
  

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

    //-------------------- if Bmeson in |y|<1 -------------------//
    if(pdg == 511  || pdg == 521  || pdg == 531 || pdg ==5122 ) {
      if(fabs(rap)<EtaCut1) IsBmeson=kTRUE;
    }
    //-----------------------------------------------------------//

    //-------------- if beauty quark in |y|<1 -----------//
    if (pdg == 5 && fabs(rap)<EtaCut1)  IsBeauty=kTRUE;
    //---------------------------------------------------//
    if(pdg==4) IsCharmQ=kTRUE;
    
    // Found the quarks
    if(pdg==5){
      // Method 2
      if((p->PdgCode()==5) && (Method2nbquarkintheeventp==0)) {
        Method2nbquarkintheeventp ++;
        Method2yquarkp = rap;
        hRapBeautyQuarkMethod2->Fill(rap,eventw);
	hQuarkMethod2->Fill(rap,p->Pt(),eventw);
      }
      if((p->PdgCode()==-5) && (Method2nbquarkintheeventn==0)) {
        Method2nbquarkintheeventn ++;
        Method2yquarkn = rap;
        hRapBeautyQuarkMethod2->Fill(rap,eventw);
	hQuarkMethod2->Fill(rap,p->Pt(),eventw);
      }
      // Method 1 now
      //printf("Found a b quark %d with rapidity %f and pseudorapidity %f and status %d\n",p->PdgCode(),rap,eta,p->GetStatusCode());
      int k1 = p->GetDaughterFirst();
      int k2 = p->GetDaughterLast();
      //printf("k1 %d and k2 %d\n",k1,k2);
      bool foundhadrons = kFALSE;
      for(int d=k1; d <= k2; d++) {
        if(d>0){
          AliVParticle *decay = fMcEvent->GetTrack(d);
          int pdgdecay = decay->PdgCode();
          //printf("Mesons %d and dpg %d, index %d\n",d,pdgdecay,d);
          if ( int(TMath::Abs(pdgdecay)/100.) == 5 || int(TMath::Abs(pdgdecay)/1000.) == 5 ) {
            foundhadrons = kTRUE;
            if(p->PdgCode()==5) {
              Method1nbquarkintheeventp ++;
              Method1yquarkp = rap;
            }
            if(p->PdgCode()==-5) {
              Method1nbquarkintheeventn ++;
              Method1yquarkn = rap;
            }
            //printf("Found a quark with rapidity %f and pdg %d\n",rap,p->PdgCode());
          }
        }
      }
      //printf("Found hadrons %d\n",foundhadrons);
      if(foundhadrons) {
	hRapBeautyQuarkMethod1->Fill(rap,eventw);
	hQuarkMethod1->Fill(rap,p->Pt(),eventw);
      }
    }
    //printf("Pdg code of the particle %d and index %d\n",p->PdgCode(),iparticle);

    // fill the single histograms
    // mesons w/ beauty quarks
    if(pdg == 511 ||
       pdg == 521 ||
       pdg == 531 ||
       pdg ==5122 ||
       pdg ==5132 ||
       pdg ==5232 ||
       pdg ==5332)                           hBeauty       ->Fill(rap,pt,eventw);
    //B^{pm}
    if(pdg == 511)                           hBpm          ->Fill(rap,pt,eventw);
    //B0
    if(pdg == 521)                           hB0           ->Fill(rap,pt,eventw);
    //Bs
    if(pdg == 531)                           hBs           ->Fill(rap,pt,eventw);
    //Lamdab
    if(pdg == 5122)                          hLambdab      ->Fill(rap,pt,eventw);


    //------------->>>> electrons are daughters
    //


    if (!( pdg == 11 )) continue;

    if ( p->IsPrimary() ) continue;

    //_____________________________ b->e _______________________________
    int num = p->GetMother();
    AliVParticle *mom = fMcEvent->GetTrack ( num );
    int ppid = TMath::Abs( mom->PdgCode() );


    //double pxMom  = mom->Px();
    //double pyMom  = mom->Py();
    double rapMom = mom->Y();
    //double ptMom  = sqrt(pxMom*pxMom+pyMom*pyMom);


    // fill the single electron histogram, if mother is beuaty and doughter is electron, for |y|<1
    is_beauty = kFALSE;
    //if (ppid == 511 || ppid == 521 || ppid == 531 || ppid == 5122) is_beauty = kTRUE; //Bpm,B0,Bs,Lambdab
    if (((ppid>=500) && (ppid<=549)) || ((ppid>=5000) && (ppid<=5499))) is_beauty = kTRUE; //open beauty mesons and baryons
    /*if (is_beauty) {
     cout << "ppid_beauty2e :" << ppid << endl;
     } */

    if(is_beauty == kTRUE){

      //printf("Found one electron with mother pdg %d and index %d\n",mom->PdgCode(),num);
      // cout << "ppid_beauty2e :" << ppid << endl;

      hRapBMom2e->Fill(rapMom,rap,eventw);
      hPtEtaElectron_be->Fill(pt, eta,eventw); // generated pt vs eta of electrons from beauty
      
      //--rapidity and pseudorapidity distributions of electrons for pt> 0.5 GeV/c ---//
      if(pdg == 11 && pt>0.5) hRapElectron_be_Pt500->Fill(rap,eventw);
      if(pdg == 11 && pt>0.2) hRapElectron_be_Pt200->Fill(rap,eventw);
      //
      if(pdg == 11               && fabs(eta)<EtaCut)   hPte_eta08_be ->Fill(pt,eventw); // all e  in |eta|<0.8
      if(p->PdgCode() ==   11 && fabs(eta)<EtaCut)   hPteP_eta08_be->Fill(pt,eventw); // all e+ in |eta|<0.8
      if(p->PdgCode() ==  -11 && fabs(eta)<EtaCut)   hPteM_eta08_be->Fill(pt,eventw); // all e- in |eta|<0.8

      if(pdg == 11               && fabs(rap)<EtaCut)   hPte_y08_be ->Fill(pt,eventw); // all e  in |y|<0.8
      if(p->PdgCode() ==   11 && fabs(rap)<EtaCut)   hPteP_y08_be->Fill(pt,eventw); // all e+ in |y|<0.8
      if(p->PdgCode() ==  -11 && fabs(rap)<EtaCut)   hPteM_y08_be->Fill(pt,eventw); // all e- in |y|<0.8


      //--------------------- if electron from beauty, in |y|<1 ---------------//
      //if(pdg == 11               && fabs(eta)<EtaCut1)     IsElectronB=kTRUE;  //
      if(pdg == 11               && fabs(rap)<EtaCut1)     IsElectronB=kTRUE;  //
      //-----------------------------------------------------------------------//

    } // b->e

    //_____________________________ b->c->e _______________________________
    //
    is_charm2e = kFALSE;
    is_beauty2charm = kFALSE;

    int num_ce = p->GetMother(); // c is the first mother of e
    AliVParticle *mom_ce = fMcEvent->GetTrack ( num_ce );
    int ppid_ce = TMath::Abs( mom_ce->PdgCode() ); // c
    //if ( ppid_ce == 411 || ppid_ce == 421 || ppid_ce == 431 || ppid_ce ==4122 ) is_charm2e = kTRUE;
    if (((ppid_ce>=400) && (ppid_ce<=439)) || ((ppid_ce>=4000) && (ppid_ce<=4399))) is_charm2e = kTRUE;
    // Check if the D meson comes from a B meson
    if(is_charm2e == kTRUE) {
      AliVParticle *gm = 0x0;
      Int_t indexx = mom_ce->GetMother(); // First mother of D
      //printf("Found one electron from charm mother pdg %d and index %d\n",mom_ce->PdgCode(),num_ce);
      while (indexx > 0) { // recursive loop to check if it comes from beauty.
        gm = fMcEvent->GetTrack ( indexx );
        int pid_gm = TMath::Abs( gm->PdgCode() );// pdg of the mother
        //printf("Next mother is %d with index %d\n",gm->PdgCode(),indexx);
        if (((pid_gm>=500) && (pid_gm<=549)) || ((pid_gm>=5000) && (pid_gm<=5499)) || (pid_gm == 5)) {
          is_beauty2charm = kTRUE;
          break;
        }
        indexx = gm->GetMother(); //
      }
      if(is_beauty2charm == kFALSE) {
        is_charm2e = kFALSE; // it is not a charm that comes from beauty.
      }
    }


    if(is_charm2e == kTRUE) {

      hRapBMom2c2e->Fill(rapMom,rap,eventw);
      hPtEtaElectron_bce->Fill(pt, eta,eventw); // generated pt vs eta of electrons from beauty to charm

      //--rapidity and pseudorapidity distributions of electrons for pt> 0.5 GeV/c ---//
      if(pdg ==  11 && pt>0.5) hRapElectron_bce_Pt500->Fill(rap,eventw);
      if(pdg ==  11 && pt>0.2) hRapElectron_bce_Pt200->Fill(rap,eventw);
      //--single leg pt distributions of electrons--------------------------------------//
      if(pdg ==  11 && fabs(eta)<EtaCut)                     hPte_eta08_bce ->Fill(pt,eventw); // all e  in |eta|<0.8
      if(p->PdgCode()       ==  11 && fabs(eta)<EtaCut)   hPteP_eta08_bce->Fill(pt,eventw); // all e+ in |eta|<0.8
      if(p->PdgCode()       == -11 && fabs(eta)<EtaCut)   hPteM_eta08_bce->Fill(pt,eventw); // all e- in |eta|<0.8
      if(pdg ==  11 && fabs(rap)<EtaCut)                     hPte_y08_bce ->Fill(pt,eventw);   // all e  in |y|<0.8
      if(p->PdgCode()       ==  11 && fabs(rap)<EtaCut)   hPteP_y08_bce->Fill(pt,eventw);   // all e+ in |y|<0.8
      if(p->PdgCode()       == -11 && fabs(rap)<EtaCut)   hPteM_y08_bce->Fill(pt,eventw);   // all e- in |y|<0.8

      //--------------------- if electron from charm, with |y|<1 --------------------//
      if(fabs(p->PdgCode()) == 11  &&
         fabs(rap)<EtaCut1)               IsElectronD=kTRUE;
      //-----------------------------------------------------------------------------//

    } // b->c->e

  }//end of single loop



  //================================== pair loop =================================
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
      // variable without Bremsstrahlung+momentum resolution
      TLorentzVector ppo_i,ppo_j;
      p_i->Momentum(ppo_i);
      p_j->Momentum(ppo_j);

      double etao_i  = ppo_i.Eta();
      double pto_i   = ppo_i.Pt();
      double phio_i  = ppo_i.Phi();

      double etao_j  = ppo_j.Eta();
      double pto_j   = ppo_j.Pt();
      double phio_j  = ppo_j.Phi();

      //double masso    = (ppo_i + ppo_j).M();
      //double pto_pair = (ppo_i + ppo_j).Pt();
      //double opAngleo = ppo_i.Angle(ppo_j.Vect());

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

      // Variable after smearing
      double eta_i  = pp_i.Eta();
      double pt_i   = pp_i.Pt();
      //double phi_i  = pp_i.Phi();

      double eta_j  = pp_j.Eta();
      double pt_j   = pp_j.Pt();
      //double phi_j  = pp_j.Phi();

      double mass     = (pp_i + pp_j).M();
      double pt_pair  = (pp_i + pp_j).Pt();
      double opAngle  = pp_i.Angle(pp_j.Vect());
      //-----------------------------------------


      //_____________________________ b->e _______________________________
      //

      int num_i = p_i->GetMother();
      AliVParticle *mom_i = fMcEvent->GetTrack ( num_i );
      int ppid_i = TMath::Abs( mom_i->PdgCode() );
      //if ( !(ppid_i == 511) || !(ppid_i == 521) || !(ppid_i == 531) || !(ppid_i == 5122) ) continue;
      //
      int num_j = p_j->GetMother();
      AliVParticle *mom_j = fMcEvent->GetTrack ( num_j );
      int ppid_j = TMath::Abs( mom_j->PdgCode() );
      //if ( !(ppid_j == 511) || !(ppid_j == 521) || !(ppid_j == 531) || !(ppid_j == 5122) ) continue;
      //
      // if mother is from Bpm, B0, Bs, Lambdab
      i_is_beauty = kFALSE;
      j_is_beauty = kFALSE;
      //if (ppid_i == 511 || ppid_i == 521 || ppid_i == 531 || ppid_i == 5122)  i_is_beauty = kTRUE;
      //if (ppid_j == 511 || ppid_j == 521 || ppid_j == 531 || ppid_j == 5122)  j_is_beauty = kTRUE;
      if (((ppid_i>=500) && (ppid_i<=549)) || ((ppid_i>=5000) && (ppid_i<=5499))) i_is_beauty = kTRUE;
      if (((ppid_j>=500) && (ppid_j<=549)) || ((ppid_j>=5000) && (ppid_j<=5499))) j_is_beauty = kTRUE;

      //__________________________________________________________________

      //_____________________________ b->c->e _______________________________
      //

      // c->e
      int num_ce_i = p_i->GetMother(); // c is the first mother of e
      AliVParticle *mom_ce_i = fMcEvent->GetTrack ( num_ce_i );
      int ppid_ce_i = TMath::Abs( mom_ce_i->PdgCode() ); // c
      //
      int num_ce_j = p_j->GetMother(); // c is the first mother of e
      AliVParticle *mom_ce_j = fMcEvent->GetTrack ( num_ce_j );
      int ppid_ce_j = TMath::Abs( mom_ce_j->PdgCode() ); // c
      //
      i_is_charm2e = kFALSE;
      j_is_charm2e = kFALSE;
      i_is_beauty2charm = kFALSE;
      j_is_beauty2charm = kFALSE;
      //if (ppid_ce_i == 411 ||  ppid_ce_i == 421 || ppid_ce_i == 431 || ppid_ce_i ==4122 )  i_is_charm2e = kTRUE;
      //if (ppid_ce_j == 411 ||  ppid_ce_j == 421 || ppid_ce_j == 431 || ppid_ce_j ==4122 )  j_is_charm2e = kTRUE;
      if (((ppid_ce_i>=400) && (ppid_ce_i<=439)) || ((ppid_ce_i>=4000) && (ppid_ce_i<=4399))) i_is_charm2e = kTRUE;
      if (((ppid_ce_j>=400) && (ppid_ce_j<=439)) || ((ppid_ce_j>=4000) && (ppid_ce_j<=4399))) j_is_charm2e = kTRUE;

      //int num_bce_j = mom_ce_j->GetMother(); // b is the first mother of c
      //AliVParticle *mom_bce_j = fMcEvent->GetTrack ( num_bce_j );
      //int ppid_bce_j = TMath::Abs( mom_bce_j->PdgCode() );// b
      // Check if the D meson comes from a B
      if(i_is_charm2e == kTRUE) {
        AliVParticle *gm_i = 0x0;
        Int_t indexx_i = mom_ce_i->GetMother(); // First mother of D
        while (indexx_i > 0) { // recursive loop to check if it comes from beauty.
          gm_i = fMcEvent->GetTrack ( indexx_i );
          int pid_gm_i = TMath::Abs( gm_i->PdgCode() );// pdg of the mother
          if (((pid_gm_i>=500) && (pid_gm_i<=549)) || ((pid_gm_i>=5000) && (pid_gm_i<=5499)) || (pid_gm_i == 5)) {
            i_is_beauty2charm = kTRUE;
            break;
          }
          indexx_i = gm_i->GetMother(); //
        }
        if(i_is_beauty2charm == kFALSE) i_is_charm2e = kFALSE; // it is not a charm that comes from beauty.
      }
      if(j_is_charm2e == kTRUE) {
        AliVParticle *gm_j = 0x0;
        Int_t indexx_j = mom_ce_j->GetMother(); // First mother of D
        while (indexx_j > 0) { // recursive loop to check if it comes from beauty.
          gm_j = fMcEvent->GetTrack ( indexx_j );
          int pid_gm_j = TMath::Abs( gm_j->PdgCode() );// pdg of the mother
          if (((pid_gm_j>=500) && (pid_gm_j<=549)) || ((pid_gm_j>=5000) && (pid_gm_j<=5499)) || (pid_gm_j == 5)) {
            j_is_beauty2charm = kTRUE;
            break;
          }
          indexx_j = gm_j->GetMother(); //
        }
        if(j_is_beauty2charm == kFALSE) j_is_charm2e = kFALSE; // it is not a charm that comes from beauty.
      }


      // track pt cuts
      double ptweight2 = pt_cut200(pt_i) * pt_cut200(pt_j) * pt_cutHigh(pt_i) * pt_cutHigh(pt_j); // pT>0.2
      double ptweight3 = pt_cut300(pt_i) * pt_cut300(pt_j) * pt_cutHigh(pt_i) * pt_cutHigh(pt_j); // pT>0.3
      double ptweight4 = pt_cut400(pt_i) * pt_cut400(pt_j) * pt_cutHigh(pt_i) * pt_cutHigh(pt_j); // pT>0.4

      // HFE R_AA scaling
      if (fScaleByRAA) {
        ptweight2 *= scale_RAA(pt_i) * scale_RAA(pt_j);
        ptweight3 *= scale_RAA(pt_i) * scale_RAA(pt_j);
        ptweight4 *= scale_RAA(pt_i) * scale_RAA(pt_j);
      }

      ptweight2 = ptweight2*eventw;
      ptweight3 = ptweight3*eventw;
      ptweight4 = ptweight4*eventw;

      //--------------------------------------- ULS pairs -------------------------------------------//

      if ( pid_i*pid_j < 0. ){

        if( ((i_is_beauty == kTRUE) || (i_is_charm2e == kTRUE))   &&  ((j_is_beauty == kTRUE) || (j_is_charm2e == kTRUE)))
        { // b->e & b->e, b->c->e & b->c->e, b->e & b->c->e, b->c->e & b->e
          // all acceptance
          hMee_ULS_simulated->Fill(mass,eventw);
          // Cut on MC variables for PHENIX acceptance
          if(fabs(etao_i)<0.5 && fabs(etao_j)<0.5) {
            hMee_ULS_eta05->Fill(mass,eventw);
          }
          if(fabs(etao_i)<0.35 && fabs(etao_j)<0.35) {
            hMee_ULS_eta035->Fill(mass,eventw);
            if(Inphenixacc(phio_i,pto_i,pid_i) && Inphenixacc(phio_j,pto_j,pid_j)) hMee_ULS_eta035_phenixacc->Fill(mass,eventw);
          }
          // Cut on smeared variables for ALICE acceptance
          if(fabs(eta_i)<EtaCut  && fabs(eta_j)<EtaCut)  {
            hMee_ULS_eta08->Fill(mass,eventw);
            hMee_ULS_eta08_pt200->Fill(mass,ptweight2);  // pt>0.2
            hMee_ULS_eta08_pt400->Fill(mass,ptweight4); // pt>0.4
            // 2D:
            hMeePtee_ULS_eta08->Fill(mass,pt_pair,eventw);
            hMeePtee_ULS_eta08_pt200->Fill(mass,pt_pair,ptweight2);  // pt>0.2
            hMeePtee_ULS_eta08_pt400->Fill(mass,pt_pair,ptweight4); // pt>0.4
            hMeeOpAngle_ULS_eta08_pt200->Fill(mass,opAngle,ptweight2); // opening angle
            if (opAngle > opAngle50) { // cut on smeared opening angle
              hMeePtee_ULS_eta08_pt200_opAngle50->Fill(mass,pt_pair,ptweight2); // pt>0.2
              hMeePtee_ULS_eta08_pt300_opAngle50->Fill(mass,pt_pair,ptweight3); // pt>0.3
              hMeePtee_ULS_eta08_pt400_opAngle50->Fill(mass,pt_pair,ptweight4); // pt>0.4
            }

          } //  |eta|<0.8
        } // all combinations for b->e and b->c->e


        if((i_is_beauty == kTRUE) && (j_is_beauty == kTRUE)) { // b->e & b->e
          // all acceptance
          hMee_ULS_simulated_be->Fill(mass,eventw);
          // Cut on MC variables for PHENIX acceptance
          if(fabs(etao_i)<0.5 && fabs(etao_j)<0.5) {
            hMee_ULS_eta05_be->Fill(mass,eventw);
          }
          if(fabs(etao_i)<0.35 && fabs(etao_j)<0.35) {
            hMee_ULS_eta035_be->Fill(mass,eventw);
            if(Inphenixacc(phio_i,pto_i,pid_i) && Inphenixacc(phio_j,pto_j,pid_j)) hMee_ULS_eta035_phenixacc_be->Fill(mass,eventw);
          }
          // Cut on smeared variables for ALICE acceptance
          if(fabs(eta_i)<EtaCut  && fabs(eta_j)<EtaCut)  {
            hMee_ULS_eta08_be->Fill(mass,eventw);
            hMee_ULS_eta08_pt200_be->Fill(mass,ptweight2);  // pt>0.2
            hMee_ULS_eta08_pt400_be->Fill(mass,ptweight4); // pt>0.4
            hMeePtee_ULS_eta08_be->Fill(mass,pt_pair,eventw);
            hMeePtee_ULS_eta08_pt200_be->Fill(mass,pt_pair,ptweight2);  // pt>0.2
            hMeePtee_ULS_eta08_pt400_be->Fill(mass,pt_pair,ptweight4); // pt>0.4
          } //  |eta|<0.8
        }// b->e & b->e


        if((i_is_charm2e == kTRUE) && (j_is_charm2e == kTRUE)){ // b->c->e & b->c->e
          // all acceptance
          hMee_ULS_simulated_bce->Fill(mass,eventw);
          // Cut on MC variables for PHENIX acceptance
          if(fabs(etao_i)<0.5 && fabs(etao_j)<0.5) {
            hMee_ULS_eta05_bce->Fill(mass,eventw);
          }
          if(fabs(etao_i)<0.35 && fabs(etao_j)<0.35) {
            hMee_ULS_eta035_bce->Fill(mass,eventw);
            if(Inphenixacc(phio_i,pto_i,pid_i) && Inphenixacc(phio_j,pto_j,pid_j)) hMee_ULS_eta035_phenixacc_bce->Fill(mass,eventw);
          }
          // Cut on smeared variables for ALICE acceptance
          if(fabs(eta_i)<EtaCut  && fabs(eta_j)<EtaCut)  {
            hMee_ULS_eta08_bce->Fill(mass,eventw);
            hMee_ULS_eta08_pt200_bce->Fill(mass,ptweight2);  // pt>0.2
            hMee_ULS_eta08_pt400_bce->Fill(mass,ptweight4); // pt>0.4
            hMeePtee_ULS_eta08_bce->Fill(mass,pt_pair,eventw);
            hMeePtee_ULS_eta08_pt200_bce->Fill(mass,pt_pair,ptweight2);  // pt>0.2
            hMeePtee_ULS_eta08_pt400_bce->Fill(mass,pt_pair,ptweight4); // pt>0.4
          } //  |eta|<0.8
        } // b->c->e & b->c->e

      }//unlike-sign



      //--------------------------------------- LS pairs -------------------------------------------//


      if ( pid_i*pid_j > 0. ){

        if( ((i_is_beauty == kTRUE) || (i_is_charm2e == kTRUE))   &&  ((j_is_beauty == kTRUE) || (j_is_charm2e == kTRUE))) {
          // all acceptance
          hMee_LS_simulated->Fill(mass,eventw);
          // Cut on MC variables for PHENIX acceptance
          if(fabs(etao_i)<0.5 && fabs(etao_j)<0.5) {
            hMee_LS_eta05->Fill(mass,eventw);
          }
          if(fabs(etao_i)<0.35 && fabs(etao_j)<0.35) {
            hMee_LS_eta035->Fill(mass,eventw);
            if(Inphenixacc(phio_i,pto_i,pid_i) && Inphenixacc(phio_j,pto_j,pid_j)) hMee_LS_eta035_phenixacc->Fill(mass,eventw);
          }
          // Cut on smeared variables for ALICE acceptance
          if(fabs(eta_i)<EtaCut  && fabs(eta_j)<EtaCut)  {
            hMee_LS_eta08->Fill(mass,eventw);
            hMee_LS_eta08_pt200->Fill(mass,ptweight2); // pt>0.2
            hMee_LS_eta08_pt400->Fill(mass,ptweight4); // pt>0.4
            // 2D:
            hMeePtee_LS_eta08->Fill(mass,pt_pair,eventw);
            hMeePtee_LS_eta08_pt200->Fill(mass,pt_pair,ptweight2); // pt>0.2
            hMeePtee_LS_eta08_pt400->Fill(mass,pt_pair,ptweight4); // pt>0.4
            hMeeOpAngle_LS_eta08_pt200->Fill(mass,opAngle,ptweight2); // opening angle
            if (opAngle > opAngle50) { // cut on smeared opening angle
              hMeePtee_LS_eta08_pt200_opAngle50->Fill(mass,pt_pair,ptweight2); // pt>0.2
              hMeePtee_LS_eta08_pt300_opAngle50->Fill(mass,pt_pair,ptweight3); // pt>0.3
              hMeePtee_LS_eta08_pt400_opAngle50->Fill(mass,pt_pair,ptweight4); // pt>0.4
            }
          } // |eta|<0.8
        } // all combinations for b->e and b->c->e


        if((i_is_beauty == kTRUE) && (j_is_beauty == kTRUE)) { // b->e & b->e
          // all acceptance
          hMee_LS_simulated_be->Fill(mass,eventw);
          // Cut on MC variables for PHENIX acceptance
          if(fabs(etao_i)<0.5 && fabs(etao_j)<0.5) {
            hMee_LS_eta05_be->Fill(mass,eventw);
          }
          if(fabs(etao_i)<0.35 && fabs(etao_j)<0.35) {
            hMee_LS_eta035_be->Fill(mass,eventw);
            if(Inphenixacc(phio_i,pto_i,pid_i) && Inphenixacc(phio_j,pto_j,pid_j)) hMee_LS_eta035_phenixacc_be->Fill(mass,eventw);
          }
          // Cut on smeared variables for ALICE acceptance
          if(fabs(eta_i)<EtaCut  && fabs(eta_j)<EtaCut)  {
            hMee_LS_eta08_be->Fill(mass,eventw);
            hMee_LS_eta08_pt200_be->Fill(mass,ptweight2);  // pt>0.2
            hMee_LS_eta08_pt400_be->Fill(mass,ptweight4); // pt>0.4
            hMeePtee_LS_eta08_be->Fill(mass,pt_pair,eventw);
            hMeePtee_LS_eta08_pt200_be->Fill(mass,pt_pair,ptweight2);  // pt>0.2
            hMeePtee_LS_eta08_pt400_be->Fill(mass,pt_pair,ptweight4); // pt>0.4
          } // |eta|<0.8
        }//beauty->e


        if((i_is_charm2e == kTRUE) && (j_is_charm2e == kTRUE)){
          // all acceptance
          hMee_LS_simulated_bce->Fill(mass,eventw);
          // Cut on MC variables for PHENIX acceptance
          if(fabs(etao_i)<0.5 && fabs(etao_j)<0.5) {
            hMee_LS_eta05_bce->Fill(mass,eventw);
          }
          if(fabs(etao_i)<0.35 && fabs(etao_j)<0.35) {
            hMee_LS_eta035_bce->Fill(mass,eventw);
            if(Inphenixacc(phio_i,pto_i,pid_i) && Inphenixacc(phio_j,pto_j,pid_j)) hMee_LS_eta035_phenixacc_bce->Fill(mass,eventw);
          }
          // Cut on smeared variables for ALICE acceptance
          if(fabs(eta_i)<EtaCut  && fabs(eta_j)<EtaCut)  {
            hMee_LS_eta08_bce->Fill(mass,eventw);
            hMee_LS_eta08_pt200_bce->Fill(mass,ptweight2);  // pt>0.2
            hMee_LS_eta08_pt400_bce->Fill(mass,ptweight4); // pt>0.4
            hMeePtee_LS_eta08_bce->Fill(mass,pt_pair,eventw);
            hMeePtee_LS_eta08_pt200_bce->Fill(mass,pt_pair,ptweight2);  // pt>0.2
            hMeePtee_LS_eta08_pt400_bce->Fill(mass,pt_pair,ptweight4); // pt>0.4
          } // |eta|<0.8
        } // b->c->e
      }// LS pairs


    } //end of pair loop for j^th particle

  }//end of pair loop for i^th particle


  hNEventsW->Fill(crosssection); // Fill with cross section
  hNEvents->Fill(0.5); // all events w/o any rapidity cut
  if(IsBmeson==kTRUE)        hNEvents->Fill(1.5); // all Bmeson events,in |y|<1
  if(IsBeauty==kTRUE)        hNEvents->Fill(2.5); // all beauty events,in |y|<1
  if(IsElectronB==kTRUE)     hNEvents->Fill(3.5); // all electrons in |y|<1 from beauty
  if(IsElectronD==kTRUE)     hNEvents->Fill(4.5); // all electrons in |y|<1 from beauty
  if(IsCharmQ==kTRUE)        hNEvents->Fill(5.5); // all electrons in |y|<1 from beauty
  fNbEvents++;

  //printf("Method 1: number of quarks positive %d and number of quarks negative %d\n",Method1nbquarkintheeventp,Method1nbquarkintheeventn);
  //printf("Method 1: y of quarks positive %f and y of quarks negative %f\n",Method1yquarkp,Method1yquarkn);
  //printf("Method 2: number of quarks positive %d and number of quarks negative %d\n",Method2nbquarkintheeventp,Method2nbquarkintheeventn);
  //printf("Method 2: y of quarks positive %f and y of quarks negative %f\n",Method2yquarkp,Method2yquarkn);
  // Correlation between the two quarks (when we found two....)
  if((Method1nbquarkintheeventp==1) && (Method1nbquarkintheeventn==1)) {
    hRapBeautyQuark2DMethod1->Fill(Method1yquarkp,Method1yquarkn);
  }
  if((Method2nbquarkintheeventp==1) && (Method2nbquarkintheeventn==1)) {
    hRapBeautyQuark2DMethod2->Fill(Method2yquarkp,Method2yquarkn);
  }


  PostData(1, fOutputList);

  return;

}
//______________________________| Terminate
void AliAnalysisTaskBeauty::Terminate(Option_t*)
{
  Info("Terminate","Start and end of Method");


  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    Printf("ERROR: Output list not available");
    return;
  }

  return;
}
//______________________________| Definations of histo etc.
void AliAnalysisTaskBeauty::CreateHistos(){
  //
  // Init the histos
  //

  // Histogram with 5 bins: for # of events
  hNEvents    = new TH1F("hNEvents","",6,0.,6.);
  hNEventsW    = new TH1F("hNEventsW","",1000,0.,1.);

  Int_t   nBinRap   =  200;
  Float_t RapMin    = -10.;
  Float_t RapMax    =  10.;
  Int_t   nBinPtHad =  720;
  Float_t PtHadMin  =   0.;
  Float_t PtHadMax  =  36.;
  Int_t   nBinPtEle =  200;
  Float_t PtEleMin  =   0.;
  Float_t PtEleMax  =  10.;

  // Histograms for mesons including beauty quark
  hQuarkMethod1  = new TH2F("hQuarkMethod1"  ,"rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hQuarkMethod2  = new TH2F("hQuarkMethod2"  ,"rap vs pt;y;p_{T}"   ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hBeauty     = new TH2F("hBeauty" ,"rap vs pt"    ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hBpm        = new TH2F("hBpm"    ,"rap vs pt"    ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hB0         = new TH2F("hB0"     ,"rap vs pt"    ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hBs         = new TH2F("hBs"     ,"rap vs pt"    ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
  hLambdab    = new TH2F("hLambdab","rap vs pt"    ,nBinRap,RapMin,RapMax,nBinPtHad,PtHadMin,PtHadMax);
 


  // Histograms for rapidity distributions
  // rap and eta do not make any difference for electrons, nevertheless, wanted to have both cases in hRapElectron and hEtaElectron
  hRapElectron_be_Pt500     = new TH1F("hRapElectron_be_Pt500"   ,"y_{e},  pt>0.5 GeV/c"   ,100,-7.2,7.2);// for comparison to FONLL
  hRapElectron_bce_Pt500    = new TH1F("hRapElectron_bce_Pt500"  ,"y_{e},  pt>0.5 GeV/c"   ,100,-7.2,7.2);// for comparison to FONLL
  hRapElectron_be_Pt200     = new TH1F("hRapElectron_be_Pt200"   ,"y_{e},  pt>0.2 GeV/c"   ,100,-7.2,7.2);// for comparison to FONLL
  hRapElectron_bce_Pt200    = new TH1F("hRapElectron_bce_Pt200"  ,"y_{e},  pt>0.2 GeV/c"   ,100,-7.2,7.2);// for comparison to FONLL
  hRapBeautyQuarkMethod1    = new TH1F("hRapBeautyQuarkMethod1"  ,"y_{b}"                  ,nBinRap,RapMin,RapMax);
  hRapBeautyQuarkMethod2    = new TH1F("hRapBeautyQuarkMethod2"  ,"y_{b}"                  ,nBinRap,RapMin,RapMax);
  hRapBeautyQuark2DMethod1  = new TH2F("hRapBeautyQuark2DMethod1","y_{b} vs y_{bbar}"      ,nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hRapBeautyQuark2DMethod2  = new TH2F("hRapBeautyQuark2DMethod2","y_{b} vs y_{bbar}"      ,nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  // hRapBmesonElectron is empty, so, wanted to fill it in another way
  hRapBMom2e         = new TH2F("hRapBMom2e"   ,"y_{Bmom} vs y_{e}",nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hRapBMom2c2e       = new TH2F("hRapBMom2c2e" ,"y_{Bmom} vs y_{e}",nBinRap,RapMin,RapMax,nBinRap,RapMin,RapMax);
  hPtEtaElectron_be  = new TH2F("hPtEtaElectron_be" ,"gen. electrons from b->e   ;p_{T,e};#eta_{e}",nBinPtEle,PtEleMin,PtEleMax,nBinRap,RapMin,RapMax);
  hPtEtaElectron_bce = new TH2F("hPtEtaElectron_bce","gen. electrons from b->c->e;p_{T,e};#eta_{e}",nBinPtEle,PtEleMin,PtEleMax,nBinRap,RapMin,RapMax);

  Int_t   nBinMass   = 800;
  Float_t MassMin    = 0.;
  Float_t MassMax    = 8.;
  Int_t   nBinPairPt = 400;
  Float_t PairPtMin  = 0.;
  Float_t PairPtMax  = 10.;


  // Histograms for Pt spectra, b-->e , bBar->e
  hPte_eta08_be       = new TH1F("hPte_eta08_be " ,"e,   |eta|<0.8" ,nBinPtEle,PtEleMin,PtEleMax);
  hPteP_eta08_be      = new TH1F("hPteP_eta08_be ","e+,  |eta|<0.8" ,nBinPtEle,PtEleMin,PtEleMax);
  hPteM_eta08_be      = new TH1F("hPteM_eta08_be ","e-,  |eta|<0.8" ,nBinPtEle,PtEleMin,PtEleMax);
  hPte_y08_be         = new TH1F("hPte_y08_be"    ,"e,   |y|<0.8"   ,nBinPtEle,PtEleMin,PtEleMax);
  hPteP_y08_be        = new TH1F("hPteP_y08_be"   ,"e+,  |y|<0.8"   ,nBinPtEle,PtEleMin,PtEleMax);
  hPteM_y08_be        = new TH1F("hPteM_y08_be"   ,"e-,  |y|<0.8"   ,nBinPtEle,PtEleMin,PtEleMax);
  // Histograms for Pt spectra, b->c->e , bBar->cBar->e
  hPte_eta08_bce       = new TH1F("hPte_eta08_bce " ,"e,   |eta|<0.8" ,nBinPtEle,PtEleMin,PtEleMax);
  hPteP_eta08_bce      = new TH1F("hPteP_eta08_bce ","e+,  |eta|<0.8" ,nBinPtEle,PtEleMin,PtEleMax);
  hPteM_eta08_bce      = new TH1F("hPteM_eta08_bce ","e-,  |eta|<0.8" ,nBinPtEle,PtEleMin,PtEleMax);
  hPte_y08_bce         = new TH1F("hPte_y08_bce"    ,"e,   |y|<0.8"   ,nBinPtEle,PtEleMin,PtEleMax);
  hPteP_y08_bce        = new TH1F("hPteP_y08_bce"   ,"e+,  |y|<0.8"   ,nBinPtEle,PtEleMin,PtEleMax);
  hPteM_y08_bce        = new TH1F("hPteM_y08_bce"   ,"e-,  |y|<0.8"   ,nBinPtEle,PtEleMin,PtEleMax);

  // Histograms for invariant mass spectra (ULS,LS), for all combinations of  b->e and b->c->e
  // titles (eta, phenix) are messed up...
  hMee_ULS_simulated    = new TH1F("hMee_ULS_simulated"  ,"e+e-"                                ,nBinMass,MassMin,MassMax);
  hMee_LS_simulated     = new TH1F("hMee_LS_simulated"   ,"e+e+ & e-e-"                         ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta05        = new TH1F("hMee_ULS_eta05"      ,"e+e-,       |eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_LS_eta05         = new TH1F("hMee_LS_eta05"       ,"e+e+ & e-e-,|eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta08        = new TH1F("hMee_ULS_eta08"      ,"e+e-,       |eta|<0.8"               ,nBinMass,MassMin,MassMax);
  hMee_LS_eta08         = new TH1F("hMee_LS_eta08"       ,"e+e- & e-e-,|eta|<0.8"               ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta035       = new TH1F("hMee_ULS_eta035"     ,"e+e-,       |eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_LS_eta035        = new TH1F("hMee_LS_eta035"      ,"e+e+ & e-e-,|eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta08_pt200  = new TH1F("hMee_ULS_eta08_pt200","e+e-,       |eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax);
  hMee_LS_eta08_pt200   = new TH1F("hMee_LS_eta08_pt200" ,"e+e- & e-e-,|eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax);
  hMee_ULS_eta08_pt400  = new TH1F("hMee_ULS_eta08_pt400","e+e-,       |eta|<0.8 , pt>0.4 GeV/c",nBinMass,MassMin,MassMax);
  hMee_LS_eta08_pt400   = new TH1F("hMee_LS_eta08_pt400" ,"e+e- & e-e-,|eta|<0.8 , pt>0.4 GeV/c",nBinMass,MassMin,MassMax);
  hMee_ULS_eta035_phenixacc = new TH1F("hMee_ULS_eta035_phenixacc"     ,"e+e-,       |eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_LS_eta035_phenixacc  = new TH1F("hMee_LS_eta035_phenixacc"      ,"e+e+ & e-e-,|eta|<1"                 ,nBinMass,MassMin,MassMax);
  //invariant mass vs pair pt
  hMeePtee_ULS_eta08       = new TH2F("hMeePtee_ULS_eta08"      ,"e+e-,       |eta|<0.8;m_{ee};p_{T,ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08        = new TH2F("hMeePtee_LS_eta08"       ,"e+e+ & e-e-,|eta|<0.8;m_{ee};p_{T,ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt200 = new TH2F("hMeePtee_ULS_eta08_pt200","e+e-,       |eta|<0.8 , pt>0.2 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt200  = new TH2F("hMeePtee_LS_eta08_pt200" ,"e+e+ & e-e-,|eta|<0.8 , pt>0.2 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt400 = new TH2F("hMeePtee_ULS_eta08_pt400","e+e-,       |eta|<0.8 , pt>0.4 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt400  = new TH2F("hMeePtee_LS_eta08_pt400" ,"e+e+ & e-e-,|eta|<0.8 , pt>0.4 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt200_opAngle50 = new TH2F("hMeePtee_ULS_eta08_pt200_opAngle50","e+e-,        |eta|<0.8 , pt>0.2 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt200_opAngle50  = new TH2F("hMeePtee_LS_eta08_pt200_opAngle50" ,"e+e- & e-e-, |eta|<0.8 , pt>0.2 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt300_opAngle50 = new TH2F("hMeePtee_ULS_eta08_pt300_opAngle50","e+e-,        |eta|<0.8 , pt>0.3 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt300_opAngle50  = new TH2F("hMeePtee_LS_eta08_pt300_opAngle50" ,"e+e- & e-e-, |eta|<0.8 , pt>0.3 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt400_opAngle50 = new TH2F("hMeePtee_ULS_eta08_pt400_opAngle50","e+e-,        |eta|<0.8 , pt>0.4 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt400_opAngle50  = new TH2F("hMeePtee_LS_eta08_pt400_opAngle50" ,"e+e- & e-e-, |eta|<0.8 , pt>0.4 GeV/c , #theta_{ee}>50 mrad;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);

  // Histograms for invariant mass spectra (ULS,LS),  b-->e , bBar->e
  // titles (eta, phenix) are messed up...
  hMee_ULS_simulated_be = new TH1F("hMee_ULS_simulated_be","e+e-"                                ,nBinMass,MassMin,MassMax);
  hMee_LS_simulated_be  = new TH1F("hMee_LS_simulated_be" ,"e+e+ & e-e-"                         ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta05_be      = new TH1F("hMee_ULS_eta05_be"     ,"e+e-,       |eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_LS_eta05_be       = new TH1F("hMee_LS_eta05_be"      ,"e+e+ & e-e-,|eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta08_be     = new TH1F("hMee_ULS_eta08_be"    ,"e+e-,       |eta|<0.8"               ,nBinMass,MassMin,MassMax);
  hMee_LS_eta08_be      = new TH1F("hMee_LS_eta08_be"     ,"e+e- & e-e-,|eta|<0.8"               ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta035_be      = new TH1F("hMee_ULS_eta035_be"     ,"e+e-,       |eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_LS_eta035_be       = new TH1F("hMee_LS_eta035_be"      ,"e+e+ & e-e-,|eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta08_pt200_be = new TH1F("hMee_ULS_eta08_pt200_be","e+e-,       |eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax);
  hMee_LS_eta08_pt200_be  = new TH1F("hMee_LS_eta08_pt200_be" ,"e+e- & e-e-,|eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax);
  hMee_ULS_eta08_pt400_be = new TH1F("hMee_ULS_eta08_pt400_be","e+e-,       |eta|<0.8 , pt>0.4 GeV/c",nBinMass,MassMin,MassMax);
  hMee_LS_eta08_pt400_be  = new TH1F("hMee_LS_eta08_pt400_be" ,"e+e- & e-e-,|eta|<0.8 , pt>0.4 GeV/c",nBinMass,MassMin,MassMax);
  hMee_ULS_eta035_phenixacc_be = new TH1F("hMee_ULS_eta035_phenixacc_be"     ,"e+e-,       |eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_LS_eta035_phenixacc_be  = new TH1F("hMee_LS_eta035_phenixacc_be"      ,"e+e+ & e-e-,|eta|<1"                 ,nBinMass,MassMin,MassMax);
  // Histograms for invariant mass vs pair pt (ULS,LS), b-->e , bBar->e
  hMeePtee_ULS_eta08_be     = new TH2F("hMeePtee_ULS_eta08_be"    ,"e+e-,       |eta|<0.8;m_{ee};p_{T,ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_be      = new TH2F("hMeePtee_LS_eta08_be"     ,"e+e+ & e-e-,|eta|<0.8;m_{ee};p_{T,ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt200_be = new TH2F("hMeePtee_ULS_eta08_pt200_be","e+e-,       |eta|<0.8 , pt>0.2 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt200_be  = new TH2F("hMeePtee_LS_eta08_pt200_be" ,"e+e+ & e-e-,|eta|<0.8 , pt>0.2 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt400_be = new TH2F("hMeePtee_ULS_eta08_pt400_be","e+e-,       |eta|<0.8 , pt>0.4 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt400_be  = new TH2F("hMeePtee_LS_eta08_pt400_be" ,"e+e+ & e-e-,|eta|<0.8 , pt>0.4 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);

  // Histograms for invariant mass spectra (ULS,LS),  b->c->e , bBar->c->e
  // titles (eta, phenix) are messed up...
  hMee_ULS_simulated_bce = new TH1F("hMee_ULS_simulated_bce","e+e-"                                ,nBinMass,MassMin,MassMax);
  hMee_LS_simulated_bce  = new TH1F("hMee_LS_simulated_bce" ,"e+e+ & e-e-"                         ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta05_bce      = new TH1F("hMee_ULS_eta05_bce"     ,"e+e-,       |eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_LS_eta05_bce       = new TH1F("hMee_LS_eta05_bce"      ,"e+e+ & e-e-,|eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta08_bce     = new TH1F("hMee_ULS_eta08_bce"    ,"e+e-,       |eta|<0.8"               ,nBinMass,MassMin,MassMax);
  hMee_LS_eta08_bce      = new TH1F("hMee_LS_eta08_bce"     ,"e+e- & e-e-,|eta|<0.8"               ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta035_bce      = new TH1F("hMee_ULS_eta035_bce"     ,"e+e-,       |eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_LS_eta035_bce       = new TH1F("hMee_LS_eta035_bce"      ,"e+e+ & e-e-,|eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_ULS_eta08_pt200_bce = new TH1F("hMee_ULS_eta08_pt200_bce","e+e-,       |eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax);
  hMee_LS_eta08_pt200_bce  = new TH1F("hMee_LS_eta08_pt200_bce" ,"e+e- & e-e-,|eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax);
  hMee_ULS_eta08_pt400_bce = new TH1F("hMee_ULS_eta08_pt400_bce","e+e-,       |eta|<0.8 , pt>0.4 GeV/c",nBinMass,MassMin,MassMax);
  hMee_LS_eta08_pt400_bce  = new TH1F("hMee_LS_eta08_pt400_bce" ,"e+e- & e-e-,|eta|<0.8 , pt>0.4 GeV/c",nBinMass,MassMin,MassMax);
  hMee_ULS_eta035_phenixacc_bce = new TH1F("hMee_ULS_eta035_phenixacc_bce"     ,"e+e-,       |eta|<1"                 ,nBinMass,MassMin,MassMax);
  hMee_LS_eta035_phenixacc_bce = new TH1F("hMee_LS_eta035_phenixacc_bce"      ,"e+e+ & e-e-,|eta|<1"                 ,nBinMass,MassMin,MassMax);
  // Histograms for invariant mass vs pair pt (ULS,LS), b->c->e , bBar->c->e
  hMeePtee_ULS_eta08_bce     = new TH2F("hMeePtee_ULS_eta08_bce"    ,"e+e-,       |eta|<0.8;m_{ee};p_{T,ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_bce      = new TH2F("hMeePtee_LS_eta08_bce"     ,"e+e+ & e-e-,|eta|<0.8;m_{ee};p_{T,ee}"               ,nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt200_bce = new TH2F("hMeePtee_ULS_eta08_pt200_bce","e+e-,       |eta|<0.8 , pt>0.2 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt200_bce  = new TH2F("hMeePtee_LS_eta08_pt200_bce" ,"e+e+ & e-e-,|eta|<0.8 , pt>0.2 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_ULS_eta08_pt400_bce = new TH2F("hMeePtee_ULS_eta08_pt400_bce","e+e-,       |eta|<0.8 , pt>0.4 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  hMeePtee_LS_eta08_pt400_bce  = new TH2F("hMeePtee_LS_eta08_pt400_bce" ,"e+e+ & e-e-,|eta|<0.8 , pt>0.4 GeV/c;m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  // opening angle
  hMeeOpAngle_ULS_eta08_pt200 = new TH2F("hMeeOpAngle_ULS_eta08_pt200","e+e-,        |eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax,128,0,3.2);
  hMeeOpAngle_LS_eta08_pt200  = new TH2F("hMeeOpAngle_LS_eta08_pt200" ,"e+e-,        |eta|<0.8 , pt>0.2 GeV/c",nBinMass,MassMin,MassMax,128,0,3.2);


  hNEvents->Sumw2();
  hNEventsW->Sumw2();
  hQuarkMethod1->Sumw2();
  hQuarkMethod2->Sumw2();
  hBeauty->Sumw2();
  hBpm->Sumw2();
  hB0->Sumw2();
  hBs->Sumw2();
  hLambdab->Sumw2();
  hRapElectron_be_Pt200->Sumw2();
  hRapElectron_be_Pt500->Sumw2();
  hRapElectron_bce_Pt200->Sumw2();
  hRapElectron_bce_Pt500->Sumw2();
  hRapBeautyQuarkMethod1->Sumw2();
  hRapBeautyQuarkMethod2->Sumw2();
  hRapBeautyQuark2DMethod1->Sumw2();
  hRapBeautyQuark2DMethod2->Sumw2();
  hRapBMom2e->Sumw2();
  hRapBMom2c2e->Sumw2();
  hPtEtaElectron_be->Sumw2();
  hPtEtaElectron_bce->Sumw2();
  hPte_eta08_be->Sumw2();
  hPteP_eta08_be->Sumw2();
  hPteM_eta08_be->Sumw2();
  hPte_y08_be->Sumw2();
  hPteP_y08_be->Sumw2();
  hPteM_y08_be->Sumw2();
  hPte_eta08_bce->Sumw2();
  hPteP_eta08_bce->Sumw2();
  hPteM_eta08_bce->Sumw2();
  hPte_y08_bce->Sumw2();
  hPteP_y08_bce->Sumw2();
  hPteM_y08_bce->Sumw2();
  hMee_ULS_simulated->Sumw2();
  hMee_LS_simulated->Sumw2();
  hMee_ULS_eta05->Sumw2();
  hMee_LS_eta05->Sumw2();
  hMee_ULS_eta08->Sumw2();
  hMee_LS_eta08->Sumw2();
  hMee_ULS_eta035->Sumw2();
  hMee_LS_eta035->Sumw2();
  hMee_ULS_eta08_pt200->Sumw2();
  hMee_LS_eta08_pt200->Sumw2();
  hMee_ULS_eta08_pt400->Sumw2();
  hMee_LS_eta08_pt400->Sumw2();
  hMee_ULS_eta035_phenixacc->Sumw2();
  hMee_LS_eta035_phenixacc->Sumw2();
  hMeePtee_ULS_eta08->Sumw2();
  hMeePtee_LS_eta08->Sumw2();
  hMeePtee_ULS_eta08_pt200->Sumw2();
  hMeePtee_LS_eta08_pt200->Sumw2();
  hMeePtee_ULS_eta08_pt400->Sumw2();
  hMeePtee_LS_eta08_pt400->Sumw2();
  hMeePtee_ULS_eta08_pt200_opAngle50->Sumw2();
  hMeePtee_LS_eta08_pt200_opAngle50->Sumw2();
  hMeePtee_ULS_eta08_pt300_opAngle50->Sumw2();
  hMeePtee_LS_eta08_pt300_opAngle50->Sumw2();
  hMeePtee_ULS_eta08_pt400_opAngle50->Sumw2();
  hMeePtee_LS_eta08_pt400_opAngle50->Sumw2();
  hMee_ULS_simulated_be->Sumw2();
  hMee_LS_simulated_be->Sumw2();
  hMee_ULS_eta05_be->Sumw2();
  hMee_LS_eta05_be->Sumw2();
  hMee_ULS_eta08_be->Sumw2();
  hMee_LS_eta08_be->Sumw2();
  hMee_ULS_eta035_be->Sumw2();
  hMee_LS_eta035_be->Sumw2();
  hMee_ULS_eta08_pt200_be->Sumw2();
  hMee_LS_eta08_pt200_be->Sumw2();
  hMee_ULS_eta08_pt400_be->Sumw2();
  hMee_LS_eta08_pt400_be->Sumw2();
  hMee_ULS_eta035_phenixacc_be->Sumw2();
  hMee_LS_eta035_phenixacc_be->Sumw2();
  hMeePtee_ULS_eta08_be->Sumw2();
  hMeePtee_LS_eta08_be->Sumw2();
  hMeePtee_ULS_eta08_pt200_be->Sumw2();
  hMeePtee_LS_eta08_pt200_be->Sumw2();
  hMeePtee_ULS_eta08_pt400_be->Sumw2();
  hMeePtee_LS_eta08_pt400_be->Sumw2();
  hMee_ULS_simulated_bce->Sumw2();
  hMee_LS_simulated_bce->Sumw2();
  hMee_ULS_eta05_bce->Sumw2();
  hMee_LS_eta05_bce->Sumw2();
  hMee_ULS_eta08_bce->Sumw2();
  hMee_LS_eta08_bce->Sumw2();
  hMee_ULS_eta035_bce->Sumw2();
  hMee_LS_eta035_bce->Sumw2();
  hMee_ULS_eta08_pt200_bce->Sumw2();
  hMee_LS_eta08_pt200_bce->Sumw2();
  hMee_ULS_eta08_pt400_bce->Sumw2();
  hMee_LS_eta08_pt400_bce->Sumw2();
  hMee_ULS_eta035_phenixacc_bce->Sumw2();
  hMee_LS_eta035_phenixacc_bce->Sumw2();
  hMeePtee_ULS_eta08_bce->Sumw2();
  hMeePtee_LS_eta08_bce->Sumw2();
  hMeePtee_ULS_eta08_pt200_bce->Sumw2();
  hMeePtee_LS_eta08_pt200_bce->Sumw2();
  hMeePtee_ULS_eta08_pt400_bce->Sumw2();
  hMeePtee_LS_eta08_pt400_bce->Sumw2();
  hMeeOpAngle_ULS_eta08_pt200->Sumw2();
  hMeeOpAngle_LS_eta08_pt200->Sumw2();


  fOutputList->Add(hNEvents);
  fOutputList->Add(hNEventsW);
  fOutputList->Add(hQuarkMethod1);
  fOutputList->Add(hQuarkMethod2);
  fOutputList->Add(hBeauty);
  fOutputList->Add(hBpm);
  fOutputList->Add(hB0);
  fOutputList->Add(hBs);
  fOutputList->Add(hLambdab);
  fOutputList->Add(hRapElectron_be_Pt200);
  fOutputList->Add(hRapElectron_be_Pt500);
  fOutputList->Add(hRapElectron_bce_Pt200);
  fOutputList->Add(hRapElectron_bce_Pt500);
  fOutputList->Add(hRapBeautyQuarkMethod1);
  fOutputList->Add(hRapBeautyQuarkMethod2);
  fOutputList->Add(hRapBeautyQuark2DMethod1);
  fOutputList->Add(hRapBeautyQuark2DMethod2);
  fOutputList->Add(hRapBMom2e);
  fOutputList->Add(hRapBMom2c2e);
  fOutputList->Add(hPtEtaElectron_be);
  fOutputList->Add(hPtEtaElectron_bce);
  fOutputList->Add(hPte_eta08_be);
  fOutputList->Add(hPteP_eta08_be);
  fOutputList->Add(hPteM_eta08_be);
  fOutputList->Add(hPte_y08_be);
  fOutputList->Add(hPteP_y08_be);
  fOutputList->Add(hPteM_y08_be);
  fOutputList->Add(hPte_eta08_bce);
  fOutputList->Add(hPteP_eta08_bce);
  fOutputList->Add(hPteM_eta08_bce);
  fOutputList->Add(hPte_y08_bce);
  fOutputList->Add(hPteP_y08_bce);
  fOutputList->Add(hPteM_y08_bce);
  fOutputList->Add(hMee_ULS_simulated);
  fOutputList->Add(hMee_LS_simulated);
  fOutputList->Add(hMee_ULS_eta05);
  fOutputList->Add(hMee_LS_eta05);
  fOutputList->Add(hMee_ULS_eta08);
  fOutputList->Add(hMee_LS_eta08);
  fOutputList->Add(hMee_ULS_eta035);
  fOutputList->Add(hMee_LS_eta035);
  fOutputList->Add(hMee_ULS_eta08_pt200);
  fOutputList->Add(hMee_LS_eta08_pt200);
  fOutputList->Add(hMee_ULS_eta08_pt400);
  fOutputList->Add(hMee_LS_eta08_pt400);
  fOutputList->Add(hMee_ULS_eta035_phenixacc);
  fOutputList->Add(hMee_LS_eta035_phenixacc);
  fOutputList->Add(hMeePtee_ULS_eta08);
  fOutputList->Add(hMeePtee_LS_eta08);
  fOutputList->Add(hMeePtee_ULS_eta08_pt200);
  fOutputList->Add(hMeePtee_LS_eta08_pt200);
  fOutputList->Add(hMeePtee_ULS_eta08_pt400);
  fOutputList->Add(hMeePtee_LS_eta08_pt400);
  fOutputList->Add(hMeePtee_ULS_eta08_pt200_opAngle50);
  fOutputList->Add(hMeePtee_LS_eta08_pt200_opAngle50);
  fOutputList->Add(hMeePtee_ULS_eta08_pt300_opAngle50);
  fOutputList->Add(hMeePtee_LS_eta08_pt300_opAngle50);
  fOutputList->Add(hMeePtee_ULS_eta08_pt400_opAngle50);
  fOutputList->Add(hMeePtee_LS_eta08_pt400_opAngle50);
  fOutputList->Add(hMee_ULS_simulated_be);
  fOutputList->Add(hMee_LS_simulated_be);
  fOutputList->Add(hMee_ULS_eta05_be);
  fOutputList->Add(hMee_LS_eta05_be);
  fOutputList->Add(hMee_ULS_eta08_be);
  fOutputList->Add(hMee_LS_eta08_be);
  fOutputList->Add(hMee_ULS_eta035_be);
  fOutputList->Add(hMee_LS_eta035_be);
  fOutputList->Add(hMee_ULS_eta08_pt200_be);
  fOutputList->Add(hMee_LS_eta08_pt200_be);
  fOutputList->Add(hMee_ULS_eta08_pt400_be);
  fOutputList->Add(hMee_LS_eta08_pt400_be);
  fOutputList->Add(hMee_ULS_eta035_phenixacc_be);
  fOutputList->Add(hMee_LS_eta035_phenixacc_be);
  fOutputList->Add(hMeePtee_ULS_eta08_be);
  fOutputList->Add(hMeePtee_LS_eta08_be);
  fOutputList->Add(hMeePtee_ULS_eta08_pt200_be);
  fOutputList->Add(hMeePtee_LS_eta08_pt200_be);
  fOutputList->Add(hMeePtee_ULS_eta08_pt400_be);
  fOutputList->Add(hMeePtee_LS_eta08_pt400_be);
  fOutputList->Add(hMee_ULS_simulated_bce);
  fOutputList->Add(hMee_LS_simulated_bce);
  fOutputList->Add(hMee_ULS_eta05_bce);
  fOutputList->Add(hMee_LS_eta05_bce);
  fOutputList->Add(hMee_ULS_eta08_bce);
  fOutputList->Add(hMee_LS_eta08_bce);
  fOutputList->Add(hMee_ULS_eta035_bce);
  fOutputList->Add(hMee_LS_eta035_bce);
  fOutputList->Add(hMee_ULS_eta08_pt200_bce);
  fOutputList->Add(hMee_LS_eta08_pt200_bce);
  fOutputList->Add(hMee_ULS_eta08_pt400_bce);
  fOutputList->Add(hMee_LS_eta08_pt400_bce);
  fOutputList->Add(hMee_ULS_eta035_phenixacc_bce);
  fOutputList->Add(hMee_LS_eta035_phenixacc_bce);
  fOutputList->Add(hMeePtee_ULS_eta08_bce);
  fOutputList->Add(hMeePtee_LS_eta08_bce);
  fOutputList->Add(hMeePtee_ULS_eta08_pt200_bce);
  fOutputList->Add(hMeePtee_LS_eta08_pt200_bce);
  fOutputList->Add(hMeePtee_ULS_eta08_pt400_bce);
  fOutputList->Add(hMeePtee_LS_eta08_pt400_bce);
  fOutputList->Add(hMeeOpAngle_ULS_eta08_pt200);
  fOutputList->Add(hMeeOpAngle_LS_eta08_pt200);

}


Double_t AliAnalysisTaskBeauty::pt_cut200(Double_t pT) {
  //Double_t weight=1.0; // maybe a bit quicker without memory allocation.
  if (pT<0.2) return 0.0;
  return 1.0;
}

Double_t AliAnalysisTaskBeauty::pt_cut300(Double_t pT) {
  if (pT<0.3) return 0.0;
  return 1.0;
}

Double_t AliAnalysisTaskBeauty::pt_cut400(Double_t pT) {
  if (pT<0.4) return 0.0;
  return 1.0;
}

Double_t AliAnalysisTaskBeauty::pt_cutHigh(Double_t pT) {
  if (pT>fPtCutHigh) return 0.0;
  return 1.0;
}

Double_t AliAnalysisTaskBeauty::scale_RAA(Double_t pT) {
  //return 8.46749e-01 + (-1.30301e-01) * pT; // fit 0-10%
  //return 9.03546e-01 + (-1.09215e-01) * pT; // fit 20-40%
  return 0.885416 + (-0.114357) * pT; // fits for 10-20% + 20-40% combined
}


Bool_t AliAnalysisTaskBeauty::Inphenixacc(Double_t phi,Double_t pt, Int_t pdg){
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
