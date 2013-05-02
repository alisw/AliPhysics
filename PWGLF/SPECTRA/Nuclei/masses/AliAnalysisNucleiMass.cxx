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
  fAOD(NULL),
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL)
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
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL)
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
  
  for(Int_t iB=0;iB<2;iB++) {

    hNevent[iB] = new TH1F("hNevent_Analyzed","Centrality(analyzed)",20,0,100);
    
    hNeventSelected[iB] = new TH1F("hNevent_Selected","Centrality(selected)",20,0,100);

    hZvertex[iB] = new TH1F("hZvertex","Vertex distribution of selected events; z vertex (cm)",240,-30,30);
    
    hTOFSignalPion[iB] = new TH1F("hTOFSignalPion","TOF signal 0.9<p_{T}<1.0; t-t_{exp}^{#pi} (ps)",1500,-1500,1500);
    
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
    
    for(Int_t i=0;i<9;i++) {
      fNsigmaTPC[iB][i] = new TH2F(namePart_par_TPC[i],namePart_title_TPC[i],250,0,5,200,-5,5);
      fNsigmaTPC[iB][i]->GetYaxis()->CenterTitle();
      fNsigmaTOF[iB][i] = new TH2F(namePart_par_TOF[i],namePart_title_TOF[i],250,0,5,200,-5,5);
      fNsigmaTOF[iB][i]->GetYaxis()->CenterTitle();
      hDeDxExp[iB][i] = new TProfile(namePart_par_ProfileTPC[i],namePart_title_ProfileTPC[i],500,0,5,0,1000,"");
      hBetaExp[iB][i] = new TProfile(namePart_par_ProfileTOF[i],namePart_title_ProfileTOF[i],400,0,5,0.4,1.05,"");
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
    
    for(Int_t iS=0;iS<18;iS++) {
      fNsigmaTPCvsP_kTOFtrue[iB][iS] = new TH2F(namePart_par_TPCvsP_kTOFtrue[iS],namePart_title_TPCvsP_kTOFtrue[iS],250,0,5,200,-5,5);
      fNsigmaTPCvsP_kTOFtrue[iB][iS]->GetYaxis()->CenterTitle();
    }
    
    char name_par_MvsP[18][60];
    char name_title_MvsP[18][150];
    
    for(Int_t i=0;i<18;i++) {
      snprintf(name_par_MvsP[i],60,"M2vsP_%s",name[i]);
      snprintf(name_title_MvsP[i],150,"M_{TOF}^{2}_%s_2#sigma_TPCcut;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4});p/|Z| (GeV/c)",name[i]);
    }
    
    for (Int_t i=0;i<18;i++) fM2vsP[iB][i] = new TH2F(name_par_MvsP[i],name_title_MvsP[i],8000,0,10,200,0,10);
    
    char name_par_MvsP_DCAxyCut[18][60];
    char name_title_MvsP_DCAxyCut[18][150];
    
    for(Int_t i=0;i<18;i++) {
      snprintf(name_par_MvsP_DCAxyCut[i],60,"M2vsP_DCAxyCut_%s",name[i]);
      snprintf(name_title_MvsP_DCAxyCut[i],150,"M_{TOF}^{2}_%s_2#sigma_TPCcut_DCAxyCut;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4});p/|Z| (GeV/c)",name[i]);
    }
    
    for (Int_t i=0;i<18;i++) fM2vsP_DCAxyCut[iB][i] = new TH2F(name_par_MvsP_DCAxyCut[i],name_title_MvsP_DCAxyCut[i],8000,0,10,200,0,10);
    
    
    Char_t par_name_nbin[nbin][30];
    
    snprintf(par_name_nbin[0],30,"0.4<Pt<0.5");
    snprintf(par_name_nbin[1],30,"0.5<Pt<0.6");
    snprintf(par_name_nbin[2],30,"0.6<Pt<0.7");
    snprintf(par_name_nbin[3],30,"0.7<Pt<0.8");
    snprintf(par_name_nbin[4],30,"0.8<Pt<0.9");
    snprintf(par_name_nbin[5],30,"0.9<Pt<1.0");
    snprintf(par_name_nbin[6],30,"1.0<Pt<1.1");
    snprintf(par_name_nbin[7],30,"1.1<Pt<1.2");
    snprintf(par_name_nbin[8],30,"1.2<Pt<1.3");
    snprintf(par_name_nbin[9],30,"1.3<Pt<1.4");
    snprintf(par_name_nbin[10],30,"1.4<Pt<1.5");
    snprintf(par_name_nbin[11],30,"1.5<Pt<1.6");
    snprintf(par_name_nbin[12],30,"1.6<Pt<1.7");
    snprintf(par_name_nbin[13],30,"1.7<Pt<1.8");
    snprintf(par_name_nbin[14],30,"1.8<Pt<1.9");
    snprintf(par_name_nbin[15],30,"1.9<Pt<2.0");
    snprintf(par_name_nbin[16],30,"2.0<Pt<2.1");
    snprintf(par_name_nbin[17],30,"2.1<Pt<2.2");
    snprintf(par_name_nbin[18],30,"2.2<Pt<2.3");
    snprintf(par_name_nbin[19],30,"2.3<Pt<2.4");
    snprintf(par_name_nbin[20],30,"2.4<Pt<2.5");
    snprintf(par_name_nbin[21],30,"2.5<Pt<2.6");
    snprintf(par_name_nbin[22],30,"2.6<Pt<2.7");
    snprintf(par_name_nbin[23],30,"2.7<Pt<2.8");
    snprintf(par_name_nbin[24],30,"2.8<Pt<2.9");
    snprintf(par_name_nbin[25],30,"2.9<Pt<3.0");
    snprintf(par_name_nbin[26],30,"3.0<Pt<3.1");
    snprintf(par_name_nbin[27],30,"3.1<Pt<3.2");
    snprintf(par_name_nbin[28],30,"3.2<Pt<3.3");
    snprintf(par_name_nbin[29],30,"3.3<Pt<3.4");
    snprintf(par_name_nbin[30],30,"3.4<Pt<3.5");
    snprintf(par_name_nbin[31],30,"3.5<Pt<3.6");
    snprintf(par_name_nbin[32],30,"3.6<Pt<3.7");
    snprintf(par_name_nbin[33],30,"3.7<Pt<3.8");
    snprintf(par_name_nbin[34],30,"3.8<Pt<3.9");
    snprintf(par_name_nbin[35],30,"3.9<Pt<4.0");
    snprintf(par_name_nbin[36],30,"4.0<Pt<4.1");
    snprintf(par_name_nbin[37],30,"4.1<Pt<4.2");
    snprintf(par_name_nbin[38],30,"4.2<Pt<4.3");
    snprintf(par_name_nbin[39],30,"4.3<Pt<4.4");
    snprintf(par_name_nbin[40],30,"4.4<Pt<4.5");
    snprintf(par_name_nbin[41],30,"4.5<Pt<4.6");
    snprintf(par_name_nbin[42],30,"4.6<Pt<4.7");
    snprintf(par_name_nbin[43],30,"4.7<Pt<4.8");
    snprintf(par_name_nbin[44],30,"4.8<Pt<4.9");
    snprintf(par_name_nbin[45],30,"4.9<Pt<5.0");
    
    
    Char_t nameDCAxy[18][nbin][120];
    Char_t titleDCAxy[18][nbin][120];
    
    Char_t nameDCAz[18][nbin][120];
    Char_t titleDCAz[18][nbin][120];
    
    Char_t nameM2CutDCAxy[18][nbin][120];
    Char_t titleM2CutDCAxy[18][nbin][120];
    
    Char_t nameM2CutGroundDCAxy[18][nbin][120];
    Char_t titleM2CutGroundDCAxy[18][nbin][120];
    
    
    for(Int_t iS=0;iS<18;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	snprintf(nameDCAxy[iS][j],120,"hDCAxy_%s_%s",name[iS],par_name_nbin[j]);
	snprintf(titleDCAxy[iS][j],120,"hDCAxy_%s_%s;DCA_{xy} (cm)",name[iS],par_name_nbin[j]);
	
	snprintf(nameDCAz[iS][j],120,"hDCAz_%s_%s",name[iS],par_name_nbin[j]);
	snprintf(titleDCAz[iS][j],120,"hDCAz_%s_%s;DCA_{z} (cm)",name[iS],par_name_nbin[j]);
	
	snprintf(nameM2CutDCAxy[iS][j],120,"hM2_CutDCAxy_%s_%s",name[iS],par_name_nbin[j]);
	snprintf(titleM2CutDCAxy[iS][j],120,"hM2_CutDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin[j]);
	
	snprintf(nameM2CutGroundDCAxy[iS][j],120,"hM2_GroundCatDCAxy_%s_%s",name[iS],par_name_nbin[j]);
	snprintf(titleM2CutGroundDCAxy[iS][j],120,"hM2_GroundCatDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin[j]);
      }
    }
    
    for(Int_t iS=0;iS<18;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	hDCAxy[iB][iS][j] = new TH1D(nameDCAxy[iS][j],titleDCAxy[iS][j],875,-3.5,3.5);//125 bins
	hDCAxy[iB][iS][j]->GetXaxis()->CenterTitle();
	hDCAz[iB][iS][j] = new TH1D(nameDCAz[iS][j],titleDCAz[iS][j],875,-3.5,3.5);//125 bins
	hDCAz[iB][iS][j]->GetXaxis()->CenterTitle();
      }
    }
    
    //for e,#mu,#pi and antiparticle (e and #mu will not be drawn)
    //the binning is chosen for #pi distribution:
    for(Int_t iSp=0;iSp<3;iSp++) {
      for(Int_t j=0;j<nbin;j++) {
	hM2CutDCAxy[iB][iSp][j] = new TH1D(nameM2CutDCAxy[iSp][j],titleM2CutDCAxy[iSp][j],600,-0.1,0.5);
	hM2CutGroundDCAxy[iB][iSp][j] = new TH1D(nameM2CutGroundDCAxy[iSp][j],titleM2CutGroundDCAxy[iSp][j],600,-0.1,0.5);
	hM2CutDCAxy[iB][iSp+9][j] = new TH1D(nameM2CutDCAxy[iSp+9][j],titleM2CutDCAxy[iSp+9][j],600,-0.1,0.5);
	hM2CutGroundDCAxy[iB][iSp+9][j] = new TH1D(nameM2CutGroundDCAxy[iSp+9][j],titleM2CutGroundDCAxy[iSp+9][j],600,-0.1,0.5);
      }
    }
    
    for(Int_t j=0;j<nbin;j++) {
      hM2CutDCAxy[iB][3][j] = new TH1D(nameM2CutDCAxy[3][j],titleM2CutDCAxy[3][j],400,0,1);
      hM2CutGroundDCAxy[iB][3][j] = new TH1D(nameM2CutGroundDCAxy[3][j],titleM2CutGroundDCAxy[3][j],400,0,1);
      hM2CutDCAxy[iB][3+9][j] = new TH1D(nameM2CutDCAxy[3+9][j],titleM2CutDCAxy[3+9][j],400,0,1);
      hM2CutGroundDCAxy[iB][3+9][j] = new TH1D(nameM2CutGroundDCAxy[3+9][j],titleM2CutGroundDCAxy[3+9][j],400,0,1);
    }
    
    for(Int_t j=0;j<nbin;j++) {
      hM2CutDCAxy[iB][4][j] = new TH1D(nameM2CutDCAxy[4][j],titleM2CutDCAxy[4][j],500,0,4);
      hM2CutGroundDCAxy[iB][4][j] = new TH1D(nameM2CutGroundDCAxy[4][j],titleM2CutGroundDCAxy[4][j],500,0,4);
      hM2CutDCAxy[iB][4+9][j] = new TH1D(nameM2CutDCAxy[4+9][j],titleM2CutDCAxy[4+9][j],500,0,4);
      hM2CutGroundDCAxy[iB][4+9][j] = new TH1D(nameM2CutGroundDCAxy[4+9][j],titleM2CutGroundDCAxy[4+9][j],500,0,4);
    }
    
    for(Int_t j=0;j<nbin;j++) {
      hM2CutDCAxy[iB][5][j] = new TH1D(nameM2CutDCAxy[5][j],titleM2CutDCAxy[5][j],500,0,6);
      hM2CutGroundDCAxy[iB][5][j] = new TH1D(nameM2CutGroundDCAxy[5][j],titleM2CutGroundDCAxy[5][j],500,0,6);
      hM2CutDCAxy[iB][5+9][j] = new TH1D(nameM2CutDCAxy[5+9][j],titleM2CutDCAxy[5+9][j],500,0,6);
      hM2CutGroundDCAxy[iB][5+9][j] = new TH1D(nameM2CutGroundDCAxy[5+9][j],titleM2CutGroundDCAxy[5+9][j],500,0,6);
    }
    
    for(Int_t j=0;j<nbin;j++) {
      hM2CutDCAxy[iB][6][j] = new TH1D(nameM2CutDCAxy[6][j],titleM2CutDCAxy[6][j],1000,0,12);
      hM2CutGroundDCAxy[iB][6][j] = new TH1D(nameM2CutGroundDCAxy[6][j],titleM2CutGroundDCAxy[6][j],1000,0,12);
      hM2CutDCAxy[iB][6+9][j] = new TH1D(nameM2CutDCAxy[6+9][j],titleM2CutDCAxy[6+9][j],1000,0,12);
      hM2CutGroundDCAxy[iB][6+9][j] = new TH1D(nameM2CutGroundDCAxy[6+9][j],titleM2CutGroundDCAxy[6+9][j],1000,0,12);
    }
    
    for(Int_t j=0;j<nbin;j++) {
      hM2CutDCAxy[iB][7][j] = new TH1D(nameM2CutDCAxy[7][j],titleM2CutDCAxy[7][j],200,0,4);
      hM2CutGroundDCAxy[iB][7][j] = new TH1D(nameM2CutGroundDCAxy[7][j],titleM2CutGroundDCAxy[7][j],200,0,4);
      hM2CutDCAxy[iB][7+9][j] = new TH1D(nameM2CutDCAxy[7+9][j],titleM2CutDCAxy[7+9][j],200,0,4);
      hM2CutGroundDCAxy[iB][7+9][j] = new TH1D(nameM2CutGroundDCAxy[7+9][j],titleM2CutGroundDCAxy[7+9][j],200,0,4);
    }
    
    for(Int_t j=0;j<nbin;j++) {
      hM2CutDCAxy[iB][8][j] = new TH1D(nameM2CutDCAxy[8][j],titleM2CutDCAxy[8][j],600,0,6);
      hM2CutGroundDCAxy[iB][8][j] = new TH1D(nameM2CutGroundDCAxy[8][j],titleM2CutGroundDCAxy[8][j],600,0,6);
      hM2CutDCAxy[iB][8+9][j] = new TH1D(nameM2CutDCAxy[8+9][j],titleM2CutDCAxy[8+9][j],600,0,6);
      hM2CutGroundDCAxy[iB][8+9][j] = new TH1D(nameM2CutGroundDCAxy[8+9][j],titleM2CutGroundDCAxy[8+9][j],600,0,6);
    }
    
    fList1[iB]->Add(hNeventSelected[iB]);
    fList1[iB]->Add(hNevent[iB]);
    fList1[iB]->Add(hZvertex[iB]);
    fList1[iB]->Add(hTOFSignalPion[iB]);
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
    for(Int_t i=0;i<18;i++) fList1[iB]->Add(fM2vsP[iB][i]);
    
    for(Int_t i=0;i<3;i++) fList1[iB]->Add(fM2vsP_NoTpcCut_DCAxyCut[iB][i]);
    for(Int_t i=0;i<18;i++) fList1[iB]->Add(fM2vsP_DCAxyCut[iB][i]);
    
    /*
      for(Int_t j=0;j<nbin;j++) {//electron
      fList1[iB]->Add(hDCAxy[iB][0][j]);
      fList1[iB]->Add(hDCAz[iB][0][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][0][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][0][j]);
      fList1[iB]->Add(hDCAxy[iB][9][j]);
      fList1[iB]->Add(hDCAz[iB][9][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][9][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][9][j]);
      }
      
      for(Int_t j=0;j<nbin;j++) {//muon
      fList1[iB]->Add(hDCAxy[iB][1][j]);
      fList1[iB]->Add(hDCAz[iB][1][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][1][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][1][j]);
      fList1[iB]->Add(hDCAxy[iB][10][j]);
      fList1[iB]->Add(hDCAz[iB][10][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][10][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][10][j]);
      }
    */
    
    for(Int_t j=0;j<nbin;j++) {//pion
      fList1[iB]->Add(hDCAxy[iB][2][j]);
      fList1[iB]->Add(hDCAz[iB][2][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][2][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][2][j]);
      fList1[iB]->Add(hDCAxy[iB][11][j]);
      fList1[iB]->Add(hDCAz[iB][11][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][11][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][11][j]);
    }

    for(Int_t j=0;j<nbin;j++) {//kaon
      fList1[iB]->Add(hDCAxy[iB][3][j]);
      fList1[iB]->Add(hDCAz[iB][3][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][3][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][3][j]);
      fList1[iB]->Add(hDCAxy[iB][12][j]);
      fList1[iB]->Add(hDCAz[iB][12][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][12][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][12][j]);
    }
    
    for(Int_t j=0;j<nbin;j++) {//proton
      fList1[iB]->Add(hDCAxy[iB][4][j]);
      fList1[iB]->Add(hDCAz[iB][4][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][4][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][4][j]);
      fList1[iB]->Add(hDCAxy[iB][13][j]);
      fList1[iB]->Add(hDCAz[iB][13][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][13][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][13][j]);
    }
    
    for(Int_t j=0;j<nbin;j++) {//deuteron
      fList1[iB]->Add(hDCAxy[iB][5][j]);
      fList1[iB]->Add(hDCAz[iB][5][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][5][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][5][j]);
      fList1[iB]->Add(hDCAxy[iB][14][j]);
      fList1[iB]->Add(hDCAz[iB][14][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][14][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][14][j]);
    }
    
    /*
      for(Int_t j=0;j<nbin;j++) {//triton
      fList1[iB]->Add(hDCAxy[iB][6][j]);
      fList1[iB]->Add(hDCAz[iB][6][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][6][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][6][j]);
      fList1[iB]->Add(hDCAxy[iB][15][j]);
      fList1[iB]->Add(hDCAz[iB][15][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][15][j]);
    fList1[iB]->Add(hM2CutGroundDCAxy[iB][15][j]);
    }
    */
    
    for(Int_t j=0;j<nbin;j++) {//He3
      fList1[iB]->Add(hDCAxy[iB][7][j]);
      fList1[iB]->Add(hDCAz[iB][7][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][7][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][7][j]);
      fList1[iB]->Add(hDCAxy[iB][16][j]);
      fList1[iB]->Add(hDCAz[iB][16][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][16][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][16][j]);
    }
    
    /*
      for(Int_t j=0;j<nbin;j++) {//4He
      fList1[iB]->Add(hDCAxy[iB][8][j]);
      fList1[iB]->Add(hDCAz[iB][8][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][8][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][8][j]);
      fList1[iB]->Add(hDCAxy[iB][17][j]);
      fList1[iB]->Add(hDCAz[iB][17][j]);
      fList1[iB]->Add(hM2CutDCAxy[iB][17][j]);
      fList1[iB]->Add(hM2CutGroundDCAxy[iB][17][j]);
      }
    */
    
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
  
  const AliAODVertex* vtxAOD = fAOD->GetPrimaryVertex();
  
  Float_t zvtx = 10000.0;

  if(vtxAOD->GetNContributors()>0)
    zvtx = vtxAOD->GetZ();
  
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

	if ((TMath::Abs(track->Eta()) > 0.8) || (track->Pt() < 0.2) || !trkFlag || !kTPC){
	  continue;
	}
	
	Double_t b[2] = {-99., -99.};
	Double_t bCov[3] = {-99., -99., -99.};
	if (!track->PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 100., b, bCov))
	  continue;
	
//	Float_t Eta = TMath::Abs(track->Eta());
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
	    if(TMath::Abs(nsigmaTPC[iS])<2.0) {
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
		fM2vsP[iBconf][iS]->Fill(M2,p);
		if(TMath::Abs(DCAxy)<DCAxyCUT) {
		  fM2vsP_DCAxyCut[iBconf][iS]->Fill(M2,p);
		}
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
	      else {//if(charge<0)
		fM2vsP[iBconf][iS+9]->Fill(M2,p);
		if(TMath::Abs(DCAxy)<DCAxyCUT) {
		  fM2vsP_DCAxyCut[iBconf][iS+9]->Fill(M2,p);
		}
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
