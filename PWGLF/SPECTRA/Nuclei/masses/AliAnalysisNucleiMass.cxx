#include "AliAnalysisNucleiMass.h"

// ROOT includes
#include <TMath.h>
#include "TChain.h"

// AliRoot includes
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliCentrality.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TProfile.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisManager.h"
#include "TFile.h"

ClassImp(AliAnalysisNucleiMass)

//_____________________________________________________________________________
AliAnalysisNucleiMass::AliAnalysisNucleiMass():
  AliAnalysisTaskSE(),
  fMC(kFALSE),
  FilterBit(16),
  NminTPCcluster(0),
  DCAzCUT(100),
  DCAxyCUT(0.1),
  kTPCcut(kTRUE),
  kTPC(0),
  kTOF(0),
  iBconf(0),
  isSignalCheck(kTRUE),
//NsigmaTPCCut(2.0),
//MomType(1),
  fAOD(NULL),
  fESD(NULL),
  fEvent(NULL),
//  fPIDResponse(NULL),
//  fmism(NULL),
  hmism(NULL),
  fchDist(NULL),
  hChDist(NULL)
/*fBetaTofVSp(NULL),
  fCentrality(NULL),
  hNevent(NULL),
  hNeventSelected(NULL),
  hTOFSignalPion(NULL),
  hEtaDistribution(NULL),
  hZvertex(NULL)*/
{
   // Default constructor (should not be used)
  fList1[0]=new TList();
  fList1[0]->SetName("results");
  
  fList1[1]=new TList();
  fList1[1]->SetName("results2");
}
//______________________________________________________________________________
AliAnalysisNucleiMass::AliAnalysisNucleiMass(const char *name):
  AliAnalysisTaskSE(name),
  fMC(kFALSE),
  FilterBit(16),
  NminTPCcluster(0),
  DCAzCUT(100),
  DCAxyCUT(0.1),
  kTPCcut(kTRUE),
  kTPC(0),
  kTOF(0),
  iBconf(0),
  isSignalCheck(kTRUE),
  //  NsigmaTPCCut(2.0),
  //MomType(1),
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  //fPIDResponse(NULL),
  //fmism(NULL),
  hmism(NULL),
  fchDist(NULL),
  hChDist(NULL)
  /*fBetaTofVSp(NULL),
  fCentrality(NULL),
  hNevent(NULL),
  hNeventSelected(NULL),
  hTOFSignalPion(NULL),
  hZvertex(NULL)*/
{
  fList1[0]=new TList();
  DefineOutput(1, TList::Class());
  fList1[0]->SetName("results");
  
  fList1[1]=new TList();
  DefineOutput(2, TList::Class());
  fList1[1]->SetName("results2");
}
//_____________________________________________________________________________
AliAnalysisNucleiMass::~AliAnalysisNucleiMass()
{
  if(fList1[0]) delete fList1[0];
  if(fList1[1]) delete fList1[1];
}
//______________________________________________________________________________
void AliAnalysisNucleiMass::UserCreateOutputObjects()
{
  
  fmism = new TFile("$ALICE_ROOT/TOF/data/TOFmismatchDistr.root");
  hmism = (TH1F *)fmism->Get("TOFmismDistr");

  fchDist = new TFile("$ALICE_ROOT/TOF/data/TOFchannelDist.root");
  hChDist = (TH1D *)fchDist->Get("hTOFchanDist");

  for(Int_t iB=0;iB<2;iB++) {

    hNevent[iB] = new TH1F("hNevent_Analyzed","Centrality(analyzed)",20,0,100);
    
    hNeventSelected[iB] = new TH1F("hNevent_Selected","Centrality(selected)",20,0,100);

    hZvertex[iB] = new TH1F("hZvertex","Vertex distribution of selected events; z vertex (cm)",240,-30,30);

    hEtaDistribution[iB][0] = new TH1F("hEtaDistribution_BeforeTRDcut","Eta distribution of the tracks_BeforeTRDcut(if there is); |#eta|",11,-0.1,1.0);
    hEtaDistribution[iB][1] = new TH1F("hEtaDistribution_TrackAnalyzed","Eta distribution of the tracks_TrackAnalyzed; |#eta|",11,-0.1,1.0);

    hTOFSignalPion[iB] = new TH1F("hTOFSignalPion","TOF signal 0.9<p_{T}<1.0; t-t_{exp}^{#pi} (ps)",1500,-1500,1500);

    hNminTPCcl[iB] = new TH1F("hNminTPCcl","hNminTPCcl",300,0,300);
   
    hPhi[iB][0] = new TH1F("hPhi_NoTRDCut","hPhi_NoTRDCut;#phi (rad.)",90,0,6.3);//each TRD supermodule is divided for 5 (DeltaPhi(TRD)=0.35 theoretical)
    hPhi[iB][1] = new TH1F("hPhi_kTRDin","hPhi_kTRDin;#phi (rad.)",90,0,6.3); 
    hPhi[iB][2] = new TH1F("hPhi_kTRDout","hPhi_kTRDout;#phi (rad.)",90,0,6.3); 
    hPhi[iB][3] = new TH1F("hPhi_kTRDin&out","hPhi_kTRDin&out;#phi (rad.)",90,0,6.3); 
    hPhi[iB][4] = new TH1F("hPhi_NoTRD","hPhi_NoTRD;#phi (rad.)",90,0,6.3); 
    hPhi[iB][5] = new TH1F("hPhi_TrackAnalyzed","hPhi_TrackAnalyzed;#phi (rad.)",90,0,6.3);

    fEtaPhi[iB][0] = new TH2F("fEtaPhi_NoTRDCut","fEtaPhi_NoTRDCut;|#eta|;#phi (rad.)",10,0.0,1.0,90,0,6.3);
    fEtaPhi[iB][1] = new TH2F("fEtaPhi_kTRDin","fEtaPhi_kTRDin;|#eta|;#phi (rad.)",10,0.0,1.0,90,0,6.3); 
    fEtaPhi[iB][2] = new TH2F("fEtaPhi_kTRDout","fEtaPhi_kTRDout;|#eta|;#phi (rad.)",10,0.0,1.0,90,0,6.3); 
    fEtaPhi[iB][3] = new TH2F("fEtaPhi_kTRDin&out","fEtaPhi_kTRDin&out;|#eta|;#phi (rad.)",10,0.0,1.0,90,0,6.3); 
    fEtaPhi[iB][4] = new TH2F("fEtaPhi_NoTRD","fEtaPhi_NoTRD;|#eta|;#phi (rad.)",10,0.0,1.0,90,0,6.3); 
    fEtaPhi[iB][5] = new TH2F("fEtaPhi_TrackAnalyzed","fEtaPhi_TrackAnalyzed;|#eta|;#phi (rad.)",10,0.0,1.0,90,0,6.3);

    char namePart[9][30];
    char namePart_par_TPC[9][40];
    char namePart_title_TPC[9][120];
    
    char namePart_par_TOF[9][40];
    char namePart_title_TOF[9][120];
    
    char namePart_par_ProfileTPC[9][40];
    char namePart_title_ProfileTPC[9][80];
    
    char namePart_par_ProfileTOF[9][40];
    char namePart_title_ProfileTOF[9][80];
        
    snprintf(namePart[0],20,"e");
    snprintf(namePart[1],20,"#mu");
    snprintf(namePart[2],20,"#pi");
    snprintf(namePart[3],20,"K");
    snprintf(namePart[4],20,"p");
    snprintf(namePart[5],20,"d");
    snprintf(namePart[6],20,"t");
    snprintf(namePart[7],20,"He3");
    snprintf(namePart[8],20,"He4");
    
    for(Int_t i=0;i<9;i++) {
      snprintf(namePart_par_TPC[i],40,"NsigmaTPC_%s",namePart[i]);
      snprintf(namePart_title_TPC[i],120,"NsigmaTPC_%s;p_{T} (GeV/c);n_{#sigma_{TPC}}^{%s}",namePart[i],namePart[i]);
      
      snprintf(namePart_par_TOF[i],40,"NsigmaTOF_%s",namePart[i]);
      snprintf(namePart_title_TOF[i],120,"NsigmaTOF_%s;p_{T} (GeV/c);n_{#sigma_{TOF}}^{%s}",namePart[i],namePart[i]);
      
      snprintf(namePart_par_ProfileTPC[i],40,"hDeDxExp_%s",namePart[i]);
      snprintf(namePart_title_ProfileTPC[i],80,"hDeDxExp_%s;p (GeV/c);dE/dx_{TPC} (a.u.)",namePart[i]);
      
      snprintf(namePart_par_ProfileTOF[i],40,"hBetaVsP_Exp_%s",namePart[i]);
      snprintf(namePart_title_ProfileTOF[i],80,"hBetaVsP_Exp%s;p (GeV/c); #beta_{TOF}",namePart[i]);
    }
    
    char namePart_par_TPCvsP_kTOFtrue[18][80];
    char namePart_title_TPCvsP_kTOFtrue[18][120];
    
    char name[18][30];
    
    snprintf(name[0],20,"e^{+}");
    snprintf(name[1],20,"#mu^{+}");
    snprintf(name[2],20,"#pi^{+}");
    snprintf(name[3],20,"K^{+}");
    snprintf(name[4],20,"p");
    snprintf(name[5],20,"d");
    snprintf(name[6],20,"t");
    snprintf(name[7],20,"He3");
    snprintf(name[8],20,"He4");
    
    snprintf(name[9],20,"e^{-}");
    snprintf(name[10],20,"#mu^{-}");
    snprintf(name[11],20,"#pi^{-}");
    snprintf(name[12],20,"K^{-}");
    snprintf(name[13],20,"#bar{p}");
    snprintf(name[14],20,"#bar{d}");
    snprintf(name[15],20,"#bar{t}");
    snprintf(name[16],20,"#bar{He3}");
    snprintf(name[17],20,"#bar{He4}");
    
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(namePart_par_TPCvsP_kTOFtrue[iS],40,"NsigmaTPCvsP_kTOFout&&kTIME_%s",name[iS]);
      snprintf(namePart_title_TPCvsP_kTOFtrue[iS],120,"NsigmaTPCvsP_kTOFout&&kTIME_%s;p (GeV/c);n_{#sigma_{TPC}}^{%s}",name[iS],name[iS]);
     }
 
    char name_par_MvsP[18][60];
    char name_title_MvsP[18][150];
    
    for(Int_t i=0;i<18;i++) {
      snprintf(name_par_MvsP[i],60,"M2vsP_%s",name[i]);
      snprintf(name_title_MvsP[i],150,"M_{TOF}^{2}_%s_2#sigma_TPCcut;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4});p/|Z| (GeV/c)",name[i]);
    }
    
    char name_par_MvsP_DCAxyCut[18][60];
    char name_title_MvsP_DCAxyCut[18][150];
    
    for(Int_t i=0;i<18;i++) {
      snprintf(name_par_MvsP_DCAxyCut[i],60,"M2vsP_DCAxyCut_%s",name[i]);
      snprintf(name_title_MvsP_DCAxyCut[i],150,"M_{TOF}^{2}_%s_2#sigma_TPCcut_DCAxyCut;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4});p/|Z| (GeV/c)",name[i]);
    }
    
    if(isSignalCheck){
      
      fdEdxVSp[iB][0] = new TH2F("fdEdxVSp","dE/dx vs p; p/|Z| (GeV/c); dE/dx_{TPC} (a.u.)",500,0,5,2000,0,1000);
      fdEdxVSp[iB][1] = new TH2F("fdEdxVSp_pos","dE/dx vs p positive charge; p/|Z| (GeV/c); dE/dx_{TPC} (a.u.)",500,0,5,2000,0,1000);
      fdEdxVSp[iB][2] = new TH2F("fdEdxVSp_neg","dE/dx vs p negative charge; p/|Z| (GeV/c); dE/dx_{TPC} (a.u.)",500,0,5,2000,0,1000);
      
      fBetaTofVSp[iB] = new TH2F("fBetaTofVSp","#beta_{TOF} vs p; p(GeV/c); #beta_{TOF}",1000,0,5,1300,0.4,1.05);
      
      fM2vsP_NoTpcCut[iB][0] = new TH2F("fM2vsP_NoTpcCut","M_{TOF}^{2} vs p; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
      fM2vsP_NoTpcCut[iB][1] = new TH2F("fM2vsP_NoTpcCut_Positive","M_{TOF}^{2} vs p Pos Part; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
      fM2vsP_NoTpcCut[iB][2] = new TH2F("fM2vsP_NoTpcCut_Negative","M_{TOF}^{2} vs p Neg Part; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
      
      fM2vsP_NoTpcCut_DCAxyCut[iB][0] = new TH2F("fM2vsP_NoTpcCut_DCAxycut","M_{TOF}^{2} vs p with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
      fM2vsP_NoTpcCut_DCAxyCut[iB][1] = new TH2F("fM2vsP_NoTpcCut_Positive_DCAxycut","M_{TOF}^{2} vs p Pos Part with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
      fM2vsP_NoTpcCut_DCAxyCut[iB][2] = new TH2F("fM2vsP_NoTpcCut_Negative_DCAxycut","M_{TOF}^{2} vs p Neg Part with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
      
      fM2vsZ[iB][0] = new TH2F("fM2vsZ","M_{TOF}^{2} vs Z^{2} Integrated p_{T};Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][1] = new TH2F("fM2vsZ_0.3pT0.5","M_{TOF}^{2} vs Z^{2} 0.3<pT<0.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][2] = new TH2F("fM2vsZ_0.5pT1.0","M_{TOF}^{2} vs Z^{2} 0.5<pT<1.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][3] = new TH2F("fM2vsZ_1.0pT1.5","M_{TOF}^{2} vs Z^{2} 1.0<pT<1.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][4] = new TH2F("fM2vsZ_1.5pT2.0","M_{TOF}^{2} vs Z^{2} 1.5<pT<2.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][5] = new TH2F("fM2vsZ_2.0pT2.5","M_{TOF}^{2} vs Z^{2} 2.0<pT<2.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][6] = new TH2F("fM2vsZ_2.5pT3.0","M_{TOF}^{2} vs Z^{2} 2.5<pT<3.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][7] = new TH2F("fM2vsZ_3.0pT3.5","M_{TOF}^{2} vs Z^{2} 3.0<pT<3.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][8] = new TH2F("fM2vsZ_3.5pT4.0","M_{TOF}^{2} vs Z^{2} 3.5<pT<4.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][9] = new TH2F("fM2vsZ_4.0pT4.5","M_{TOF}^{2} vs Z^{2} 4.0<pT<4.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][10] = new TH2F("fM2vsZ_0.0pT1.0","M_{TOF}^{2} vs Z^{2} 2.0<pT<2.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][11] = new TH2F("fM2vsZ_1.0pT2.0","M_{TOF}^{2} vs Z^{2} 2.5<pT<3.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][12] = new TH2F("fM2vsZ_2.0pT3.0","M_{TOF}^{2} vs Z^{2} 3.0<pT<3.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][13] = new TH2F("fM2vsZ_3.0pT4.0","M_{TOF}^{2} vs Z^{2} 3.5<pT<4.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZ[iB][14] = new TH2F("fM2vsZ_4.0pT4.5","M_{TOF}^{2} vs Z^{2} 4.0<pT<4.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);

      fM2vsZwithTPC[iB][0] = new TH2F("fM2vsZwithTPC","M_{TOF}^{2} vs Z^{2} Integrated p_{T} withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][1] = new TH2F("fM2vsZwithTPC_0.3pT0.5","M_{TOF}^{2} vs Z^{2} 0.3<pT<0.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][2] = new TH2F("fM2vsZwithTPC_0.5pT1.0","M_{TOF}^{2} vs Z^{2} 0.5<pT<1.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][3] = new TH2F("fM2vsZwithTPC_1.0pT1.5","M_{TOF}^{2} vs Z^{2} 1.0<pT<1.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][4] = new TH2F("fM2vsZwithTPC_1.5pT2.0","M_{TOF}^{2} vs Z^{2} 1.5<pT<2.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][5] = new TH2F("fM2vsZwithTPC_2.0pT2.5","M_{TOF}^{2} vs Z^{2} 2.0<pT<2.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][6] = new TH2F("fM2vsZwithTPC_2.5pT3.0","M_{TOF}^{2} vs Z^{2} 2.5<pT<3.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][7] = new TH2F("fM2vsZwithTPC_3.0pT3.5","M_{TOF}^{2} vs Z^{2} 3.0<pT<3.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][8] = new TH2F("fM2vsZwithTPC_3.5pT4.0","M_{TOF}^{2} vs Z^{2} 3.5<pT<4.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][9] = new TH2F("fM2vsZwithTPC_4.0pT4.5","M_{TOF}^{2} vs Z^{2} 4.0<pT<4.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][10] = new TH2F("fM2vsZwithTPC_0.0pT1.0","M_{TOF}^{2} vs Z^{2} 2.0<pT<2.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][11] = new TH2F("fM2vsZwithTPC_1.0pT2.0","M_{TOF}^{2} vs Z^{2} 2.5<pT<3.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][12] = new TH2F("fM2vsZwithTPC_2.0pT3.0","M_{TOF}^{2} vs Z^{2} 3.0<pT<3.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][13] = new TH2F("fM2vsZwithTPC_3.0pT4.0","M_{TOF}^{2} vs Z^{2} 3.5<pT<4.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      fM2vsZwithTPC[iB][14] = new TH2F("fM2vsZwithTPC_4.0pT4.5","M_{TOF}^{2} vs Z^{2} 4.0<pT<4.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
      
      for(Int_t i=0;i<9;i++) {
	fNsigmaTPC[iB][i] = new TH2F(namePart_par_TPC[i],namePart_title_TPC[i],125,0,5,100,-5,5);
	fNsigmaTPC[iB][i]->GetYaxis()->CenterTitle();
	fNsigmaTOF[iB][i] = new TH2F(namePart_par_TOF[i],namePart_title_TOF[i],125,0,5,100,-5,5);
	fNsigmaTOF[iB][i]->GetYaxis()->CenterTitle();
	hDeDxExp[iB][i] = new TProfile(namePart_par_ProfileTPC[i],namePart_title_ProfileTPC[i],500,0,5,0,1000,"");
	hBetaExp[iB][i] = new TProfile(namePart_par_ProfileTOF[i],namePart_title_ProfileTOF[i],400,0,5,0.4,1.05,"");
      }
      
      for(Int_t iS=0;iS<18;iS++) {
	fNsigmaTPCvsP_kTOFtrue[iB][iS] = new TH2F(namePart_par_TPCvsP_kTOFtrue[iS],namePart_title_TPCvsP_kTOFtrue[iS],125,0,5,100,-5,5);
	fNsigmaTPCvsP_kTOFtrue[iB][iS]->GetYaxis()->CenterTitle();
      }
      
      for (Int_t i=0;i<18;i++) fM2vsP[iB][i] = new TH2F(name_par_MvsP[i],name_title_MvsP[i],8000,0,10,200,0,10);
      
      for (Int_t i=0;i<18;i++) fM2vsP_DCAxyCut[iB][i] = new TH2F(name_par_MvsP_DCAxyCut[i],name_title_MvsP_DCAxyCut[i],8000,0,10,200,0,10);
      
    }
    
    else{//IsSignalCheck is kFALSE
      
      fdEdxVSp[iB][0] = new TH2F("fdEdxVSp","dE/dx vs p; p/|Z| (GeV/c); dE/dx_{TPC} (a.u.)",1,0,5,1,0,1000);
      fdEdxVSp[iB][1] = new TH2F("fdEdxVSp_pos","dE/dx vs p positive charge; p/|Z| (GeV/c); dE/dx_{TPC} (a.u.)",1,0,5,1,0,1000);
      fdEdxVSp[iB][2] = new TH2F("fdEdxVSp_neg","dE/dx vs p negative charge; p/|Z| (GeV/c); dE/dx_{TPC} (a.u.)",1,0,5,1,0,1000);
    
      fBetaTofVSp[iB] = new TH2F("fBetaTofVSp","#beta_{TOF} vs p; p(GeV/c); #beta_{TOF}",1,0,5,1,0.4,1.05);
      
      fM2vsP_NoTpcCut[iB][0] = new TH2F("fM2vsP_NoTpcCut","M_{TOF}^{2} vs p; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,8);//1250,...,80
      fM2vsP_NoTpcCut[iB][1] = new TH2F("fM2vsP_NoTpcCut_Positive","M_{TOF}^{2} vs p Pos Part; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,8);
      fM2vsP_NoTpcCut[iB][2] = new TH2F("fM2vsP_NoTpcCut_Negative","M_{TOF}^{2} vs p Neg Part; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,8);
      
      fM2vsP_NoTpcCut_DCAxyCut[iB][0] = new TH2F("fM2vsP_NoTpcCut_DCAxycut","M_{TOF}^{2} vs p with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,8);
      fM2vsP_NoTpcCut_DCAxyCut[iB][1] = new TH2F("fM2vsP_NoTpcCut_Positive_DCAxycut","M_{TOF}^{2} vs p Pos Part with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,8);
      fM2vsP_NoTpcCut_DCAxyCut[iB][2] = new TH2F("fM2vsP_NoTpcCut_Negative_DCAxycut","M_{TOF}^{2} vs p Neg Part with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,8);

      /*fM2vsP_NoTpcCut[iB][0] = new TH2F("fM2vsP_NoTpcCut","M_{TOF}^{2} vs p; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,10);
      fM2vsP_NoTpcCut[iB][1] = new TH2F("fM2vsP_NoTpcCut_Positive","M_{TOF}^{2} vs p Pos Part; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,10);
      fM2vsP_NoTpcCut[iB][2] = new TH2F("fM2vsP_NoTpcCut_Negative","M_{TOF}^{2} vs p Neg Part; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,10);
      
      fM2vsP_NoTpcCut_DCAxyCut[iB][0] = new TH2F("fM2vsP_NoTpcCut_DCAxycut","M_{TOF}^{2} vs p with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,10);
      fM2vsP_NoTpcCut_DCAxyCut[iB][1] = new TH2F("fM2vsP_NoTpcCut_Positive_DCAxycut","M_{TOF}^{2} vs p Pos Part with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,10);
      fM2vsP_NoTpcCut_DCAxyCut[iB][2] = new TH2F("fM2vsP_NoTpcCut_Negative_DCAxycut","M_{TOF}^{2} vs p Neg Part with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",1,0,10,1,0,10);*/
      
      fM2vsZ[iB][0] = new TH2F("fM2vsZ","M_{TOF}^{2} vs Z^{2} Integrated p_{T};Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][1] = new TH2F("fM2vsZ_0.3pT0.5","M_{TOF}^{2} vs Z^{2} 0.3<pT<0.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][2] = new TH2F("fM2vsZ_0.5pT1.0","M_{TOF}^{2} vs Z^{2} 0.5<pT<1.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][3] = new TH2F("fM2vsZ_1.0pT1.5","M_{TOF}^{2} vs Z^{2} 1.0<pT<1.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][4] = new TH2F("fM2vsZ_1.5pT2.0","M_{TOF}^{2} vs Z^{2} 1.5<pT<2.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][5] = new TH2F("fM2vsZ_2.0pT2.5","M_{TOF}^{2} vs Z^{2} 2.0<pT<2.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][6] = new TH2F("fM2vsZ_2.5pT3.0","M_{TOF}^{2} vs Z^{2} 2.5<pT<3.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][7] = new TH2F("fM2vsZ_3.0pT3.5","M_{TOF}^{2} vs Z^{2} 3.0<pT<3.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][8] = new TH2F("fM2vsZ_3.5pT4.0","M_{TOF}^{2} vs Z^{2} 3.5<pT<4.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][9] = new TH2F("fM2vsZ_4.0pT4.5","M_{TOF}^{2} vs Z^{2} 4.0<pT<4.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][10] = new TH2F("fM2vsZ_0.0pT1.0","M_{TOF}^{2} vs Z^{2} 2.0<pT<2.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][11] = new TH2F("fM2vsZ_1.0pT2.0","M_{TOF}^{2} vs Z^{2} 2.5<pT<3.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][12] = new TH2F("fM2vsZ_2.0pT3.0","M_{TOF}^{2} vs Z^{2} 3.0<pT<3.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][13] = new TH2F("fM2vsZ_3.0pT4.0","M_{TOF}^{2} vs Z^{2} 3.5<pT<4.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZ[iB][14] = new TH2F("fM2vsZ_4.0pT4.5","M_{TOF}^{2} vs Z^{2} 4.0<pT<4.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);

      fM2vsZwithTPC[iB][0] = new TH2F("fM2vsZwithTPC","M_{TOF}^{2} vs Z^{2} Integrated p_{T} withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][1] = new TH2F("fM2vsZwithTPC_0.3pT0.5","M_{TOF}^{2} vs Z^{2} 0.3<pT<0.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][2] = new TH2F("fM2vsZwithTPC_0.5pT1.0","M_{TOF}^{2} vs Z^{2} 0.5<pT<1.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][3] = new TH2F("fM2vsZwithTPC_1.0pT1.5","M_{TOF}^{2} vs Z^{2} 1.0<pT<1.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][4] = new TH2F("fM2vsZwithTPC_1.5pT2.0","M_{TOF}^{2} vs Z^{2} 1.5<pT<2.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][5] = new TH2F("fM2vsZwithTPC_2.0pT2.5","M_{TOF}^{2} vs Z^{2} 2.0<pT<2.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][6] = new TH2F("fM2vsZwithTPC_2.5pT3.0","M_{TOF}^{2} vs Z^{2} 2.5<pT<3.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][7] = new TH2F("fM2vsZwithTPC_3.0pT3.5","M_{TOF}^{2} vs Z^{2} 3.0<pT<3.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][8] = new TH2F("fM2vsZwithTPC_3.5pT4.0","M_{TOF}^{2} vs Z^{2} 3.5<pT<4.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][9] = new TH2F("fM2vsZwithTPC_4.0pT4.5","M_{TOF}^{2} vs Z^{2} 4.0<pT<4.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][10] = new TH2F("fM2vsZwithTPC_0.0pT1.0","M_{TOF}^{2} vs Z^{2} 2.0<pT<2.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][11] = new TH2F("fM2vsZwithTPC_1.0pT2.0","M_{TOF}^{2} vs Z^{2} 2.5<pT<3.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][12] = new TH2F("fM2vsZwithTPC_2.0pT3.0","M_{TOF}^{2} vs Z^{2} 3.0<pT<3.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][13] = new TH2F("fM2vsZwithTPC_3.0pT4.0","M_{TOF}^{2} vs Z^{2} 3.5<pT<4.0 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      fM2vsZwithTPC[iB][14] = new TH2F("fM2vsZwithTPC_4.0pT4.5","M_{TOF}^{2} vs Z^{2} 4.0<pT<4.5 withTPCcut;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",1,-4,4,1,0,10);
      
      for(Int_t i=0;i<9;i++) {
	fNsigmaTPC[iB][i] = new TH2F(namePart_par_TPC[i],namePart_title_TPC[i],1,0,5,1,-5,5);
	fNsigmaTPC[iB][i]->GetYaxis()->CenterTitle();
	fNsigmaTOF[iB][i] = new TH2F(namePart_par_TOF[i],namePart_title_TOF[i],1,0,5,1,-5,5);
	fNsigmaTOF[iB][i]->GetYaxis()->CenterTitle();
	hDeDxExp[iB][i] = new TProfile(namePart_par_ProfileTPC[i],namePart_title_ProfileTPC[i],1,0,5,0,1000,"");
	hBetaExp[iB][i] = new TProfile(namePart_par_ProfileTOF[i],namePart_title_ProfileTOF[i],1,0,5,0.4,1.05,"");
      }
      
      for(Int_t iS=0;iS<18;iS++) {
	fNsigmaTPCvsP_kTOFtrue[iB][iS] = new TH2F(namePart_par_TPCvsP_kTOFtrue[iS],namePart_title_TPCvsP_kTOFtrue[iS],1,0,5,1,-5,5);
	fNsigmaTPCvsP_kTOFtrue[iB][iS]->GetYaxis()->CenterTitle();
      }
      
      /*for (Int_t i=0;i<18;i++) fM2vsP[iB][i] = new TH2F(name_par_MvsP[i],name_title_MvsP[i],1,0,10,1,0,10);
     
	for (Int_t i=0;i<18;i++) fM2vsP_DCAxyCut[iB][i] = new TH2F(name_par_MvsP_DCAxyCut[i],name_title_MvsP_DCAxyCut[i],1,0,10,1,0,10);*/
      
      for (Int_t i=0;i<18;i++) fM2vsP[iB][i] = new TH2F(name_par_MvsP[i],name_title_MvsP[i],1000,0,6,60,0,6);//1250,0,10,80,0,8
      
      for (Int_t i=0;i<18;i++) fM2vsP_DCAxyCut[iB][i] = new TH2F(name_par_MvsP_DCAxyCut[i],name_title_MvsP_DCAxyCut[i],1000,0,6,60,0,6);//1250,0,10,80,0,8

    }

    Char_t namefEtaSpecies[18][300];
    Char_t titlefEtaSpecies[18][300];
    
    for(Int_t iS=0;iS<18;iS++) {
      sprintf(namefEtaSpecies[iS],"fEtaSpecies_kTOF_%s",name[iS]);
      sprintf(titlefEtaSpecies[iS],"fEtaSpecies_kTOF_%s;|#eta|;p_{T} GeV/c",name[iS]);
    }
    
    for(Int_t iS=0;iS<18;iS++) {
      fEtaSpecies[iB][iS] = new TH2F(namefEtaSpecies[iS],titlefEtaSpecies[iS],10,0,1,200,0,10);
    }

    Char_t namefPhiSpecies[18][300];
    Char_t titlefPhiSpecies[18][300];
    
    for(Int_t iS=0;iS<18;iS++) {
      sprintf(namefPhiSpecies[iS],"fPhiSpecies_kTOF_%s",name[iS]);
      sprintf(titlefPhiSpecies[iS],"fPhiSpecies_kTOF_%s;#phi (rad.);p_{T} GeV/c",name[iS]);
    }
    
    for(Int_t iS=0;iS<18;iS++) {
      fPhiSpecies[iB][iS] = new TH2F(namefPhiSpecies[iS],titlefPhiSpecies[iS],90,0,6.3,200,0,10);
    }

    Float_t binPt[nbin+1];
    for(Int_t i=0;i<nbin+1;i++) {
      binPt[i]=0.4+0.1*i;
    }
        
    Char_t par_name_nbin[nbin][30];
    
    for(Int_t j=0;j<nbin;j++) {
      snprintf(par_name_nbin[j],30,"%.1f<Pt<%.1f",binPt[j],binPt[j+1]);
    }
    
    Char_t par_name_nbin_pbin[nbin][30];
    
    for(Int_t j=0;j<nbin;j++) {
      snprintf(par_name_nbin_pbin[j],30,"%.1f<P<%.1f",binPt[j],binPt[j+1]);
    }
    
    Char_t par_name_nbin_pTpcbin[nbin][30];
    
    for(Int_t j=0;j<nbin;j++) {
      snprintf(par_name_nbin_pTpcbin[j],30,"%.1f<PTpc<%.1f",binPt[j],binPt[j+1]);
    }
    
    Char_t nameDCAxy[3][18][nbin][120];
    Char_t titleDCAxy[3][18][nbin][120];
    
    Char_t nameDCAz[3][18][nbin][120];
    Char_t titleDCAz[3][18][nbin][120];
    
    Char_t nameM2CutDCAxy[3][18][nbin][120];
    Char_t titleM2CutDCAxy[3][18][nbin][120];
    
    Char_t nameM2CutGroundDCAxy[3][18][nbin][120];
    Char_t titleM2CutGroundDCAxy[3][18][nbin][120];
        
    Char_t nameM2BkgMism[3][nbin][120];
    Char_t titleM2BkgMism[3][nbin][120];

    for(Int_t iS=0;iS<18;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	snprintf(nameDCAxy[0][iS][j],120,"hDCAxy_%s_%s",name[iS],par_name_nbin[j]);
	snprintf(titleDCAxy[0][iS][j],120,"hDCAxy_%s_%s;DCA_{xy} (cm)",name[iS],par_name_nbin[j]);
	
	snprintf(nameDCAz[0][iS][j],120,"hDCAz_%s_%s",name[iS],par_name_nbin[j]);
	snprintf(titleDCAz[0][iS][j],120,"hDCAz_%s_%s;DCA_{z} (cm)",name[iS],par_name_nbin[j]);
	
	snprintf(nameM2CutDCAxy[0][iS][j],120,"hM2_CutDCAxy_%s_%s",name[iS],par_name_nbin[j]);
	snprintf(titleM2CutDCAxy[0][iS][j],120,"hM2_CutDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin[j]);
	
	snprintf(nameM2CutGroundDCAxy[0][iS][j],120,"hM2_GroundCatDCAxy_%s_%s",name[iS],par_name_nbin[j]);
	snprintf(titleM2CutGroundDCAxy[0][iS][j],120,"hM2_GroundCatDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin[j]);


	snprintf(nameDCAxy[1][iS][j],120,"hDCAxy_pbin_%s_%s",name[iS],par_name_nbin_pbin[j]);
	snprintf(titleDCAxy[1][iS][j],120,"hDCAxy_pbin_%s_%s;DCA_{xy} (cm)",name[iS],par_name_nbin_pbin[j]);
	
	snprintf(nameDCAz[1][iS][j],120,"hDCAz_pbin_%s_%s",name[iS],par_name_nbin_pbin[j]);
	snprintf(titleDCAz[1][iS][j],120,"hDCAz_pbin_%s_%s;DCA_{z} (cm)",name[iS],par_name_nbin_pbin[j]);
	
	snprintf(nameM2CutDCAxy[1][iS][j],120,"hM2_pbin_CutDCAxy_%s_%s",name[iS],par_name_nbin_pbin[j]);
	snprintf(titleM2CutDCAxy[1][iS][j],120,"hM2_pbin_CutDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin_pbin[j]);
	
	snprintf(nameM2CutGroundDCAxy[1][iS][j],120,"hM2_pbin_GroundCatDCAxy_%s_%s",name[iS],par_name_nbin_pbin[j]);
	snprintf(titleM2CutGroundDCAxy[1][iS][j],120,"hM2_pbin_GroundCatDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin_pbin[j]);
      

	snprintf(nameDCAxy[2][iS][j],120,"hDCAxy_pTpcbin_%s_%s",name[iS],par_name_nbin_pTpcbin[j]);
	snprintf(titleDCAxy[2][iS][j],120,"hDCAxy_pbin_%s_%s;DCA_{xy} (cm)",name[iS],par_name_nbin_pTpcbin[j]);
	
	snprintf(nameDCAz[2][iS][j],120,"hDCAz_pTpcbin_%s_%s",name[iS],par_name_nbin_pTpcbin[j]);
	snprintf(titleDCAz[2][iS][j],120,"hDCAz_pTpcbin_%s_%s;DCA_{z} (cm)",name[iS],par_name_nbin_pTpcbin[j]);
	
	snprintf(nameM2CutDCAxy[2][iS][j],120,"hM2_pTpcbin_CutDCAxy_%s_%s",name[iS],par_name_nbin_pTpcbin[j]);
	snprintf(titleM2CutDCAxy[2][iS][j],120,"hM2_pTpcbin_CutDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin_pTpcbin[j]);
	
	snprintf(nameM2CutGroundDCAxy[2][iS][j],120,"hM2_pTpcbin_GroundCatDCAxy_%s_%s",name[iS],par_name_nbin_pTpcbin[j]);
	snprintf(titleM2CutGroundDCAxy[2][iS][j],120,"hM2_pTpcbin_GroundCatDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin_pTpcbin[j]);
      }
    }
    
    for(Int_t j=0;j<nbin;j++) {
      snprintf(nameM2BkgMism[0][j],120,"hM2_BkgMism_%s",par_name_nbin[j]);
      snprintf(titleM2BkgMism[0][j],120,"hM2_BkgMism_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",par_name_nbin[j]);

      snprintf(nameM2BkgMism[1][j],120,"hM2_pbin_BkgMism_%s",par_name_nbin_pbin[j]);
      snprintf(titleM2BkgMism[1][j],120,"hM2_pbin_BkgMism_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",par_name_nbin_pbin[j]);

      snprintf(nameM2BkgMism[2][j],120,"hM2_pTpcbin_BkgMism_%s",par_name_nbin_pTpcbin[j]);
      snprintf(titleM2BkgMism[2][j],120,"hM2_pTpcbin_BkgMism_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",par_name_nbin_pTpcbin[j]);
    }

    for(Int_t iS=0;iS<18;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	hDCAxy[iB][iS][j] = new TH1D(nameDCAxy[0][iS][j],titleDCAxy[0][iS][j],875,-3.5,3.5);//125 bins
	hDCAxy[iB][iS][j]->GetXaxis()->CenterTitle();
	hDCAz[iB][iS][j] = new TH1D(nameDCAz[0][iS][j],titleDCAz[0][iS][j],875,-3.5,3.5);//125 bins
	hDCAz[iB][iS][j]->GetXaxis()->CenterTitle();
	
	hDCAxy_pbin[iB][iS][j] = new TH1D(nameDCAxy[1][iS][j],titleDCAxy[1][iS][j],875,-3.5,3.5);//125 bins
	hDCAxy_pbin[iB][iS][j]->GetXaxis()->CenterTitle();
	hDCAz_pbin[iB][iS][j] = new TH1D(nameDCAz[1][iS][j],titleDCAz[1][iS][j],875,-3.5,3.5);//125 bins
	hDCAz_pbin[iB][iS][j]->GetXaxis()->CenterTitle();
	
	hDCAxy_pTpcbin[iB][iS][j] = new TH1D(nameDCAxy[2][iS][j],titleDCAxy[2][iS][j],875,-3.5,3.5);//125 bins
	hDCAxy_pTpcbin[iB][iS][j]->GetXaxis()->CenterTitle();
	hDCAz_pTpcbin[iB][iS][j] = new TH1D(nameDCAz[2][iS][j],titleDCAz[2][iS][j],875,-3.5,3.5);//125 bins
	hDCAz_pTpcbin[iB][iS][j]->GetXaxis()->CenterTitle();
      }
    }

    for(Int_t iBinMom=0;iBinMom<3;iBinMom++) {
      for(Int_t j=0;j<nbin;j++) {
	hM2BkgMism[iB][iBinMom][j]=new TH1D(nameM2BkgMism[iBinMom][j],titleM2BkgMism[iBinMom][j],500,0,6);//125 bins
	hM2BkgMism[iB][iBinMom][j]->GetXaxis()->CenterTitle();
      }
    }
    
    const Int_t BinM2pT[9]={1,1,600,250,500,500,1000,400,600};
    const Float_t RangeM2min[9]={0.0,0.0,-0.1,0.0,0.0,0.0,0.0,0.0,0.0};
    const Float_t RangeM2max[9]={1.0,1.0,0.5,2.0,4.0,6.0,12.0,4.0,6.0};
    
    

    for(Int_t iSp=0;iSp<9;iSp++) {
      for(Int_t j=0;j<nbin;j++) {
	hM2CutDCAxy[iB][iSp][j] = new TH1D(nameM2CutDCAxy[0][iSp][j],titleM2CutDCAxy[0][iSp][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	hM2CutGroundDCAxy[iB][iSp][j] = new TH1D(nameM2CutGroundDCAxy[0][iSp][j],titleM2CutGroundDCAxy[0][iSp][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	hM2CutDCAxy[iB][iSp+9][j] = new TH1D(nameM2CutDCAxy[0][iSp+9][j],titleM2CutDCAxy[0][iSp+9][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	hM2CutGroundDCAxy[iB][iSp+9][j] = new TH1D(nameM2CutGroundDCAxy[0][iSp+9][j],titleM2CutGroundDCAxy[0][iSp+9][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	
	hM2CutDCAxy_pbin[iB][iSp][j] = new TH1D(nameM2CutDCAxy[1][iSp][j],titleM2CutDCAxy[1][iSp][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	hM2CutGroundDCAxy_pbin[iB][iSp][j] = new TH1D(nameM2CutGroundDCAxy[1][iSp][j],titleM2CutGroundDCAxy[1][iSp][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	hM2CutDCAxy_pbin[iB][iSp+9][j] = new TH1D(nameM2CutDCAxy[1][iSp+9][j],titleM2CutDCAxy[1][iSp+9][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	hM2CutGroundDCAxy_pbin[iB][iSp+9][j] = new TH1D(nameM2CutGroundDCAxy[1][iSp+9][j],titleM2CutGroundDCAxy[1][iSp+9][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	
	hM2CutDCAxy_pTpcbin[iB][iSp][j] = new TH1D(nameM2CutDCAxy[2][iSp][j],titleM2CutDCAxy[2][iSp][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	hM2CutGroundDCAxy_pTpcbin[iB][iSp][j] = new TH1D(nameM2CutGroundDCAxy[2][iSp][j],titleM2CutGroundDCAxy[2][iSp][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	hM2CutDCAxy_pTpcbin[iB][iSp+9][j] = new TH1D(nameM2CutDCAxy[2][iSp+9][j],titleM2CutDCAxy[2][iSp+9][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
	hM2CutGroundDCAxy_pTpcbin[iB][iSp+9][j] = new TH1D(nameM2CutGroundDCAxy[2][iSp+9][j],titleM2CutGroundDCAxy[2][iSp+9][j],BinM2pT[iSp],RangeM2min[iSp],RangeM2max[iSp]);
      }
    }
    
    fList1[iB]->Add(hNeventSelected[iB]);
    fList1[iB]->Add(hNevent[iB]);
    fList1[iB]->Add(hZvertex[iB]);
    fList1[iB]->Add(hNminTPCcl[iB]);
    fList1[iB]->Add(hTOFSignalPion[iB]);
    for(Int_t i=0;i<2;i++)fList1[iB]->Add(hEtaDistribution[iB][i]);
    for(Int_t i=0;i<6;i++) fList1[iB]->Add(hPhi[iB][i]);
    for(Int_t i=0;i<6;i++) fList1[iB]->Add(fEtaPhi[iB][i]);
    for(Int_t iS=0;iS<18;iS++) fList1[iB]->Add(fEtaSpecies[iB][iS]);
    for(Int_t iS=0;iS<18;iS++) fList1[iB]->Add(fPhiSpecies[iB][iS]);

    for(Int_t iS=0;iS<18;iS++) fList1[iB]->Add(fNsigmaTPCvsP_kTOFtrue[iB][iS]);
    
    for(Int_t i=0;i<3;i++) fList1[iB]->Add(fdEdxVSp[iB][i]);
    for(Int_t i=0;i<9;i++) fList1[iB]->Add(hDeDxExp[iB][i]);
    fList1[iB]->Add(fBetaTofVSp[iB]);
    for(Int_t i=0;i<9;i++) fList1[iB]->Add(hBetaExp[iB][i]);
    for(Int_t i=0;i<9;i++) fList1[iB]->Add(fNsigmaTPC[iB][i]);
    for(Int_t i=0;i<9;i++) fList1[iB]->Add(fNsigmaTOF[iB][i]);
    
    for(Int_t i=0;i<15;i++) fList1[iB]->Add(fM2vsZ[iB][i]);
    for(Int_t i=0;i<15;i++) fList1[iB]->Add(fM2vsZwithTPC[iB][i]);
    
    for(Int_t i=0;i<3;i++) fList1[iB]->Add(fM2vsP_NoTpcCut[iB][i]);
    for(Int_t i=3;i<6;i++) {//for(Int_t i=0;i<18;i++)
      fList1[iB]->Add(fM2vsP[iB][i]);
      fList1[iB]->Add(fM2vsP[iB][i+9]);//via-^
    }

    for(Int_t i=0;i<3;i++) fList1[iB]->Add(fM2vsP_NoTpcCut_DCAxyCut[iB][i]);
    for(Int_t i=3;i<6;i++) {//for(Int_t i=0;i<18;i++)
      fList1[iB]->Add(fM2vsP_DCAxyCut[iB][i]);
      fList1[iB]->Add(fM2vsP_DCAxyCut[iB][i+9]);//via-^
    }    

    if(MomType & 1) {
      for(Int_t iSp=3;iSp<6;iSp++) {//for(Int_t iSp=2;iSp<9;iSp++)
	for(Int_t j=0;j<nbin;j++) {
	  fList1[iB]->Add(hDCAxy[iB][iSp][j]);
	  fList1[iB]->Add(hDCAz[iB][iSp][j]);
	  fList1[iB]->Add(hM2CutDCAxy[iB][iSp][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy[iB][iSp][j]);
	  fList1[iB]->Add(hDCAxy[iB][iSp+9][j]);
	  fList1[iB]->Add(hDCAz[iB][iSp+9][j]);
	  fList1[iB]->Add(hM2CutDCAxy[iB][iSp+9][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy[iB][iSp+9][j]);
	}
      }
      for(Int_t j=0;j<nbin;j++) {//3he
	
	fList1[iB]->Add(hDCAxy[iB][7][j]);
	fList1[iB]->Add(hDCAz[iB][7][j]);
	fList1[iB]->Add(hM2CutDCAxy[iB][7][j]);
	fList1[iB]->Add(hM2CutGroundDCAxy[iB][7][j]);
	fList1[iB]->Add(hDCAxy[iB][7+9][j]);
	fList1[iB]->Add(hDCAz[iB][7+9][j]);
	fList1[iB]->Add(hM2CutDCAxy[iB][7+9][j]);
	fList1[iB]->Add(hM2CutGroundDCAxy[iB][7+9][j]);

      }
    }
    if(MomType & 2) {
      for(Int_t iSp=3;iSp<6;iSp++) {//for(Int_t iSp=2;iSp<9;iSp++)
	for(Int_t j=0;j<nbin;j++) {
	  fList1[iB]->Add(hDCAxy_pbin[iB][iSp][j]);
	  fList1[iB]->Add(hDCAz_pbin[iB][iSp][j]);
	  fList1[iB]->Add(hM2CutDCAxy_pbin[iB][iSp][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy_pbin[iB][iSp][j]);
	  fList1[iB]->Add(hDCAxy_pbin[iB][iSp+9][j]);
	  fList1[iB]->Add(hDCAz_pbin[iB][iSp+9][j]);
	  fList1[iB]->Add(hM2CutDCAxy_pbin[iB][iSp+9][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy_pbin[iB][iSp+9][j]);
	}
      }
      for(Int_t j=0;j<nbin;j++) {//3he
	  fList1[iB]->Add(hDCAxy_pbin[iB][7][j]);
	  fList1[iB]->Add(hDCAz_pbin[iB][7][j]);
	  fList1[iB]->Add(hM2CutDCAxy_pbin[iB][7][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy_pbin[iB][7][j]);
	  fList1[iB]->Add(hDCAxy_pbin[iB][7+9][j]);
	  fList1[iB]->Add(hDCAz_pbin[iB][7+9][j]);
	  fList1[iB]->Add(hM2CutDCAxy_pbin[iB][7+9][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy_pbin[iB][7+9][j]);
      }
    }
    if(MomType & 4) {
      for(Int_t iSp=3;iSp<6;iSp++) {//for(Int_t iSp=2;iSp<9;iSp++)
	for(Int_t j=0;j<nbin;j++) {
	  fList1[iB]->Add(hDCAxy_pTpcbin[iB][iSp][j]);
	  fList1[iB]->Add(hDCAz_pTpcbin[iB][iSp][j]);
	  fList1[iB]->Add(hM2CutDCAxy_pTpcbin[iB][iSp][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy_pTpcbin[iB][iSp][j]);
	  fList1[iB]->Add(hDCAxy_pTpcbin[iB][iSp+9][j]);
	  fList1[iB]->Add(hDCAz_pTpcbin[iB][iSp+9][j]);
	  fList1[iB]->Add(hM2CutDCAxy_pTpcbin[iB][iSp+9][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy_pTpcbin[iB][iSp+9][j]);
	}
      }
      for(Int_t j=0;j<nbin;j++) {//3he
	  fList1[iB]->Add(hDCAxy_pTpcbin[iB][7][j]);
	  fList1[iB]->Add(hDCAz_pTpcbin[iB][7][j]);
	  fList1[iB]->Add(hM2CutDCAxy_pTpcbin[iB][7][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy_pTpcbin[iB][7][j]);
	  fList1[iB]->Add(hDCAxy_pTpcbin[iB][7+9][j]);
	  fList1[iB]->Add(hDCAz_pTpcbin[iB][7+9][j]);
	  fList1[iB]->Add(hM2CutDCAxy_pTpcbin[iB][7+9][j]);
	  fList1[iB]->Add(hM2CutGroundDCAxy_pTpcbin[iB][7+9][j]);
      }
    }
    
    if(MomType & 1) {
      for(Int_t j=0;j<nbin;j++) {
	fList1[iB]->Add(hM2BkgMism[iB][0][j]);
      }
    }
    if(MomType & 2) {
      for(Int_t j=0;j<nbin;j++) {
	fList1[iB]->Add(hM2BkgMism[iB][1][j]);
      }
    }
    if(MomType & 4) {
      for(Int_t j=0;j<nbin;j++) {
	fList1[iB]->Add(hM2BkgMism[iB][2][j]);
      }
    }

    // Post output data.
    PostData(1, fList1[0]);
    PostData(2, fList1[1]);
  }//close the iB loop
}
//______________________________________________________________________________
void AliAnalysisNucleiMass::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fAOD && !fESD){
    Printf("%s:%d AODEvent and ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    return;
  }
  
  if(fESD) fEvent = fESD;
  else fEvent = fAOD;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse(); // data member di tipo "const AliPIDResponse *fPIDResponse;"

  //Centrality
  Float_t v0Centr  = -10.;
  //Float_t trkCentr  = -10.;
  AliCentrality *centrality = fEvent->GetCentrality();
  if (centrality){
    v0Centr  = centrality->GetCentralityPercentile("V0M"); // VZERO
    //trkCentr = centrality->GetCentralityPercentile("TRK"); // TPC
  }
  
  Double_t fBfield=fEvent->GetMagneticField();
  if(fBfield<0.0) iBconf=0;//B--
  else iBconf=1;//B++

  hNeventSelected[iBconf]->Fill(v0Centr);//selected events

  //const AliAODVertex* vtxEVENT = (AliAODVertex*) fEvent->GetPrimaryVertex();
  
  const AliVVertex* vtxEVENT = fEvent->GetPrimaryVertex();

  Float_t zvtx = 10000.0;

  if(vtxEVENT->GetNContributors()>0)
    zvtx = vtxEVENT->GetZ();
  
   hZvertex[iBconf]->Fill(zvtx);
    
  if(TMath::Abs(zvtx) < 10.0){ // consistency cut on centrality selection AND d(zPrimaryVertez;NominalPointInteraction)<10cm
    
    //Bool_t isTrack=1;
    
    Int_t nTracks = fEvent->GetNumberOfTracks();
    
    if(v0Centr>=fCentrality[0] && v0Centr<=fCentrality[1]) {//window cut centrality open
      
      hNevent[iBconf]->Fill(v0Centr);//analyzed events
                  
      for(Int_t iT = 0; iT < nTracks; iT++) { // loop on the tracks
	AliVTrack* track = (AliVTrack *) fEvent->GetTrack(iT);
	
	if (!track){
	  continue;
	}
	
	Bool_t trkFlag = 0;
	trkFlag = ((AliAODTrack *) track)->TestFilterBit(FilterBit);
	//TestFilterBit(16) -- Standard Cuts with very loose DCA: GetStandardITSTPCTrackCuts2011(kFALSE) && SetMaxDCAToVertexXY(2.4) && SetMaxDCAToVertexZ(3.2) && SetDCaToVertex2D(kTRUE)
	//TestFilterBit(32) (STARDARD) -- Standard Cuts with very tight DCA cut ( 7sigma^primaries: 7*(0.0015+0.0050/pt^1.1) ) : GetStandardITSTPCTrackCuts2011(). 
		
       	Int_t NTpcCls=track->GetTPCNcls();
	if(NTpcCls>NminTPCcluster) kTPC=kTRUE;
	else kTPC=kFALSE;
	
	Float_t etaAbs = TMath::Abs(track->Eta());
		
	//if(etaAbs<EtaLimit[0] && etaAbs>EtaLimit[1]) continue;
      
	if(etaAbs<EtaLimit[0]) continue;
	if(etaAbs>EtaLimit[1]) continue;

	if ((track->Pt() < 0.2) || !trkFlag || !kTPC){
	  continue;
	}	
	
	Float_t phi= track->Phi();
	hPhi[iBconf][0]->Fill(phi);
	fEtaPhi[iBconf][0]->Fill(etaAbs,phi);

	Int_t iTRDtemp=1;
	if(kTRDana) {//TRD analysis
	  if((track->GetStatus() & AliVTrack::kTRDin) && (track->GetStatus() & AliVTrack::kTRDout)) {
	    iTRDtemp=4;//YES TRD
	  }
	  else if (!(track->GetStatus() & AliVTrack::kTRDin) && !(track->GetStatus() & AliVTrack::kTRDout)){
	    iTRDtemp=2;//NO TRD
	  }
	}
	else {//NO TRD analysis
	  iTRDtemp=1;
	}
	
	if(track->GetStatus() & AliVTrack::kTRDin) {
	  hPhi[iBconf][1]->Fill(phi);
	  fEtaPhi[iBconf][1]->Fill(etaAbs,phi);
	}
	
	if(track->GetStatus() & AliVTrack::kTRDout) {
	  hPhi[iBconf][2]->Fill(phi);
	  fEtaPhi[iBconf][2]->Fill(etaAbs,phi);
	}
	if((track->GetStatus() & AliVTrack::kTRDin) && (track->GetStatus() & AliVTrack::kTRDout)) {
	    //YES TRD
	    hPhi[iBconf][3]->Fill(phi);
	    fEtaPhi[iBconf][3]->Fill(etaAbs,phi);
	}
	else if (!(track->GetStatus() & AliVTrack::kTRDin) && !(track->GetStatus() & AliVTrack::kTRDout)){
	    //NO TRD
	    hPhi[iBconf][4]->Fill(phi);
	    fEtaPhi[iBconf][4]->Fill(etaAbs,phi);
	}	

	hEtaDistribution[iBconf][0]->Fill(etaAbs);

	if(!(iTRDtemp & iTRD)) {
	  continue;
	}
	
	hNminTPCcl[iBconf]->Fill(NTpcCls);

	Double_t b[2] = {-99., -99.};
	Double_t bCov[3] = {-99., -99., -99.};
	if (!track->PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 100., b, bCov))
	  continue;
	
	//Float_t etaAbs = TMath::Abs(track->Eta());
	Float_t charge = (Float_t)track->Charge();
	Float_t p = track->P();
	Float_t pt = track->Pt();
	Float_t dedx = track->GetTPCsignal();
	Float_t tof  = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
	Float_t pTPC = track->GetTPCmomentum();
	Float_t beta = 0.0;
	Float_t M2 = 1000.0;
	Float_t M = 1000.0;
	Float_t Z2 = 1000.0;
	Float_t DCAxy = b[0];
	Float_t DCAz = b[1];
	
	hEtaDistribution[iBconf][1]->Fill(etaAbs);

	hPhi[iBconf][5]->Fill(phi);
	fEtaPhi[iBconf][5]->Fill(etaAbs,phi);

	if(TMath::Abs(DCAz)>DCAzCUT)//CUT ON DCAz
	  continue;
	
	Bool_t kTpcPure;
	kTpcPure = track->GetTPCsignal()>10;
	if(kTpcPure==kFALSE) continue;
	
	kTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);

       	Float_t nsigmaTPC[9];
	Float_t nsigmaTOF[9];
	
	Double_t expdedx[9];
	
	for(Int_t iS=0;iS < 9;iS++){ //TPC expected signal
	  expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
	}
	
	for(Int_t iS=0;iS < 9;iS++){
	  nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
	  fNsigmaTPC[iBconf][iS]->Fill(pt,nsigmaTPC[iS]);
	  hDeDxExp[iBconf][iS]->Fill(pTPC,expdedx[iS]);
	}
	fdEdxVSp[iBconf][0]->Fill(pTPC,dedx);
	if(charge>0) fdEdxVSp[iBconf][1]->Fill(pTPC,dedx);
	else fdEdxVSp[iBconf][2]->Fill(pTPC,dedx);
	
	Float_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.877837,2.817402,1.408701,1.877837};
	
	Double_t exptimes[9]; // TOF expected times
	track->GetIntegratedTimes(exptimes);
	exptimes[5] = exptimes[0] / p * massOverZ[5] * TMath::Sqrt(1+p*p/massOverZ[5]/massOverZ[5]);
	exptimes[6] = exptimes[0] / p * massOverZ[6] * TMath::Sqrt(1+p*p/massOverZ[6]/massOverZ[6]);
	exptimes[7] = exptimes[0] / p * massOverZ[7] * TMath::Sqrt(1+p*p/massOverZ[7]/massOverZ[7]);
	exptimes[8] = exptimes[0] / p * massOverZ[8] * TMath::Sqrt(1+p*p/massOverZ[8]/massOverZ[8]);
	
	beta=exptimes[0];//expected times of the electron (it will be diveded for the T.o.f.)
	
	if(kTOF) {
	  for(Int_t iS=0;iS < 9;iS++){
	    nsigmaTOF[iS] = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType) iS);
	    fNsigmaTOF[iBconf][iS]->Fill(pt,nsigmaTOF[iS]);
	    hBetaExp[iBconf][iS]->Fill(p,beta/exptimes[iS]);
	  }
	  if(pt>0.9 && pt<1.0) hTOFSignalPion[iBconf]->Fill(tof-exptimes[2]);
	  beta=beta/tof;
	  fBetaTofVSp[iBconf]->Fill(p,beta);
	}
	
	Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,#mu,#pi,K,p,d,t,3He,4He
	Int_t FlagPid = 0;

	Float_t binPt[nbin+1];
	for(Int_t i=0;i<nbin+1;i++) {
	  binPt[i]=0.4+i*0.1;
	}
	
	//M2 background distribution from mismatch (START):

	//Hit channel in the TOF: 
	Int_t channel = (Int_t)(4334.09-4758.36*etaAbs-1989.71*etaAbs*etaAbs+1957.62*etaAbs*etaAbs*etaAbs);
	
	// get distance
	channel = channel % 8736;
	Float_t distIP = hChDist->GetBinContent(channel);
	
	// generate random time
	Float_t timeRandom = hmism->GetRandom() + distIP*3.35655419905265973e+00; 
	Float_t betaRandom = 1.0;
	Float_t M2Random = 1000.0;
	
	if(kTOF) {
	  betaRandom=exptimes[0];
	  if(timeRandom!=0.0)betaRandom=betaRandom/timeRandom;
	  M2Random = (p*p*(1-betaRandom*betaRandom))/(betaRandom*betaRandom);
	  if(MomType & 1) {
	    for(Int_t j=0;j<nbin;j++) {
	      if(pt>binPt[j] && pt<binPt[j+1]) {
		hM2BkgMism[iBconf][0][j]->Fill(M2Random);
		break;
	      }
	    }
	  }
	  if(MomType & 2) {
	    for(Int_t j=0;j<nbin;j++) {
	      if(p>binPt[j] && p<binPt[j+1]) {
		hM2BkgMism[iBconf][1][j]->Fill(M2Random);
		break;
	      }
	      }
	  }
	  if(MomType & 4) {
	    for(Int_t j=0;j<nbin;j++) {
	      if(pTPC>binPt[j] && pTPC<binPt[j+1]) {
		hM2BkgMism[iBconf][2][j]->Fill(M2Random);
		break;
	      }
	    }	  
	  }
	}		
	//M2 background distribution from mismatch (FINISH).

	Float_t binCutPt[10] = {0.3,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5};

	Float_t binCutLargePt[6] = {0.0,1.0,2.0,3.0,4.0,5.0};

	if(kTOF) {
	  
	  M2 = (p*p*(1-beta*beta))/(beta*beta);
	  
	  fM2vsP_NoTpcCut[iBconf][0]->Fill(M2,p);
	  if(TMath::Abs(DCAxy)<DCAxyCUT) fM2vsP_NoTpcCut_DCAxyCut[iBconf][0]->Fill(M2,p);
	  
	  if(M2>0.0) {
	    M=TMath::Sqrt(M2);
	    Z2 = TMath::Power(dedx/fPIDResponse->GetTPCResponse().GetExpectedSignal(pTPC*massOverZ[4]/M, AliPID::kProton),0.862);
	    fM2vsZ[iBconf][0]->Fill(charge*TMath::Sqrt(Z2),M2);
	    
	    for(Int_t i=0;i<9;i++) {
	      if(pt>binCutPt[i] && pt<binCutPt[i+1]){
		fM2vsZ[iBconf][i+1]->Fill(charge*TMath::Sqrt(Z2),M2);
		break;
	      }
	    }
	    for(Int_t i=0;i<5;i++) {
	      if(pt>binCutLargePt[i] && pt<binCutLargePt[i+1]){
		fM2vsZ[iBconf][10+i]->Fill(charge*TMath::Sqrt(Z2),M2);
		break;		
	      }	      
	    }
	  }
	  
	  if(charge>0) {
	    fM2vsP_NoTpcCut[iBconf][1]->Fill(M2,p);
	    if(TMath::Abs(DCAxy)<DCAxyCUT) fM2vsP_NoTpcCut_DCAxyCut[iBconf][1]->Fill(M2,p);
	    for(Int_t iS=0;iS<9;iS++){
	      fNsigmaTPCvsP_kTOFtrue[iBconf][iS]->Fill(p,nsigmaTPC[iS]);
	    }
	  }
	  else {//else charge<0
	    fM2vsP_NoTpcCut[iBconf][2]->Fill(M2,p);
	    if(TMath::Abs(DCAxy)<DCAxyCUT) fM2vsP_NoTpcCut_DCAxyCut[iBconf][2]->Fill(M2,p);
	    
	    for(Int_t iS=0;iS < 9;iS++){
	      fNsigmaTPCvsP_kTOFtrue[iBconf][iS+9]->Fill(p,nsigmaTPC[iS]);
	    }
	  }
	  
	  for(Int_t iS=0;iS<9;iS++) {
	    if(TMath::Abs(nsigmaTPC[iS])<NsigmaTPCCut) {
	      FlagPid += ((Int_t)TMath::Power(2,iS));
	    }
	  }
	 	 
	  if(M2>0.0) {
	    for(Int_t iS=0;iS<9;iS++) {
	      if(FlagPid & stdFlagPid[iS]) {
		fM2vsZwithTPC[iBconf][0]->Fill(charge*TMath::Sqrt(Z2),M2);
		for(Int_t i=0;i<9;i++) {
		  if(pt>binCutPt[i] && pt<binCutPt[i+1]) {
		    fM2vsZwithTPC[iBconf][i+1]->Fill(charge*TMath::Sqrt(Z2),M2);
		    break;
		  }
		}
		for(Int_t i=0;i<5;i++) {
		  if(pt>binCutLargePt[i] && pt<binCutLargePt[i+1]){
		    fM2vsZwithTPC[iBconf][10+i]->Fill(charge*TMath::Sqrt(Z2),M2);
		    break;		
		  }	      
		}
	      }
	    }
	  }	  
	  
	  for(Int_t iS=0;iS<9;iS++) {
	    if(FlagPid & stdFlagPid[iS] || !kTPCcut) {
	      if(charge>0) {
		fEtaSpecies[iBconf][iS]->Fill(etaAbs,pt);
		fPhiSpecies[iBconf][iS]->Fill(phi,pt);
		fM2vsP[iBconf][iS]->Fill(M2,p);
		if(TMath::Abs(DCAxy)<DCAxyCUT) {
		  fM2vsP_DCAxyCut[iBconf][iS]->Fill(M2,p);
		}
		if(MomType & 1) {
		  for(Int_t j=0;j<nbin;j++) {
		    if(pt>binPt[j] && pt<binPt[j+1]) {
		      hDCAxy[iBconf][iS][j]->Fill(DCAxy);
		      hDCAxy[iBconf][iS][j]->Fill(-DCAxy);
		      hDCAz[iBconf][iS][j]->Fill(DCAz);
		      hDCAz[iBconf][iS][j]->Fill(-DCAz);
		      if(TMath::Abs(DCAxy)<DCAxyCUT) {
			hM2CutDCAxy[iBconf][iS][j]->Fill(M2);
		      }
		      if(TMath::Abs(DCAxy+0.5)<DCAxyCUT) hM2CutGroundDCAxy[iBconf][iS][j]->Fill(M2);
		      break;
		    }
		  }
		}
		if(MomType & 2) {
		  for(Int_t j=0;j<nbin;j++) {
		    if(p>binPt[j] && p<binPt[j+1]) {
		      hDCAxy_pbin[iBconf][iS][j]->Fill(DCAxy);
		      hDCAxy_pbin[iBconf][iS][j]->Fill(-DCAxy);
		      hDCAz_pbin[iBconf][iS][j]->Fill(DCAz);
		      hDCAz_pbin[iBconf][iS][j]->Fill(-DCAz);
		      if(TMath::Abs(DCAxy)<DCAxyCUT) {
			hM2CutDCAxy_pbin[iBconf][iS][j]->Fill(M2);
		      }
		      if(TMath::Abs(DCAxy+0.5)<DCAxyCUT) hM2CutGroundDCAxy_pbin[iBconf][iS][j]->Fill(M2);
		      break;
		    }
		  }
		}
		if(MomType & 4) {
		  for(Int_t j=0;j<nbin;j++) {
		    if(pTPC>binPt[j] && pTPC<binPt[j+1]) {
		      hDCAxy_pTpcbin[iBconf][iS][j]->Fill(DCAxy);
		      hDCAxy_pTpcbin[iBconf][iS][j]->Fill(-DCAxy);
		      hDCAz_pTpcbin[iBconf][iS][j]->Fill(DCAz);
		      hDCAz_pTpcbin[iBconf][iS][j]->Fill(-DCAz);
		      if(TMath::Abs(DCAxy)<DCAxyCUT) {
			hM2CutDCAxy_pTpcbin[iBconf][iS][j]->Fill(M2);
		      }
		      if(TMath::Abs(DCAxy+0.5)<DCAxyCUT) hM2CutGroundDCAxy_pTpcbin[iBconf][iS][j]->Fill(M2);
		      break;
		    }
		  }
		}
	      }
	      else {//if(charge<0)
		fM2vsP[iBconf][iS+9]->Fill(M2,p);
		fEtaSpecies[iBconf][iS+9]->Fill(etaAbs,pt);
		fPhiSpecies[iBconf][iS+9]->Fill(phi,pt);
		if(TMath::Abs(DCAxy)<DCAxyCUT) {
		  fM2vsP_DCAxyCut[iBconf][iS+9]->Fill(M2,p);
		}
		if(MomType & 1) {
		  for(Int_t j=0;j<nbin;j++) {
		    if(pt>binPt[j] && pt<binPt[j+1]) {
		      hDCAxy[iBconf][iS+9][j]->Fill(DCAxy);
		      hDCAxy[iBconf][iS+9][j]->Fill(-DCAxy);
		      hDCAz[iBconf][iS+9][j]->Fill(DCAz);
		      hDCAz[iBconf][iS+9][j]->Fill(-DCAz);
		      if(TMath::Abs(DCAxy)<DCAxyCUT) {
			hM2CutDCAxy[iBconf][iS+9][j]->Fill(M2);
		      }
		      if(TMath::Abs(DCAxy+0.5)<DCAxyCUT) hM2CutGroundDCAxy[iBconf][iS+9][j]->Fill(M2);
		      break;
		    }
		  }
		}
		if(MomType & 2) {
		  for(Int_t j=0;j<nbin;j++) {
		    if(p>binPt[j] && p<binPt[j+1]) {
		      hDCAxy_pbin[iBconf][iS+9][j]->Fill(DCAxy);
		      hDCAxy_pbin[iBconf][iS+9][j]->Fill(-DCAxy);
		      hDCAz_pbin[iBconf][iS+9][j]->Fill(DCAz);
		      hDCAz_pbin[iBconf][iS+9][j]->Fill(-DCAz);
		      if(TMath::Abs(DCAxy)<DCAxyCUT) {
			hM2CutDCAxy_pbin[iBconf][iS+9][j]->Fill(M2);
		      }
		      if(TMath::Abs(DCAxy+0.5)<DCAxyCUT) hM2CutGroundDCAxy_pbin[iBconf][iS+9][j]->Fill(M2);
		      break;
		    }
		  }
		}
		if(MomType & 4) {
		  for(Int_t j=0;j<nbin;j++) {
		    if(pTPC>binPt[j] && pTPC<binPt[j+1]) {
		      hDCAxy_pTpcbin[iBconf][iS+9][j]->Fill(DCAxy);
		      hDCAxy_pTpcbin[iBconf][iS+9][j]->Fill(-DCAxy);
		      hDCAz_pTpcbin[iBconf][iS+9][j]->Fill(DCAz);
		      hDCAz_pTpcbin[iBconf][iS+9][j]->Fill(-DCAz);
		      if(TMath::Abs(DCAxy)<DCAxyCUT) {
			hM2CutDCAxy_pTpcbin[iBconf][iS+9][j]->Fill(M2);
		      }
		      if(TMath::Abs(DCAxy+0.5)<DCAxyCUT) hM2CutGroundDCAxy_pTpcbin[iBconf][iS+9][j]->Fill(M2);
		      break;
		    }
		  }
		}
	      }
	    }
	  }
	}// close (KTOF request)
	
	
      } // end track loop
      
    } //window cut centrality close
    
  }  
  
}

//_____________________________________________________________________________
void AliAnalysisNucleiMass::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
