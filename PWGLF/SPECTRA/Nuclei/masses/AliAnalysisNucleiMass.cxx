
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
  fAOD(NULL),
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL),
  fList1(new TList()),
  hNevent(NULL),
  hZvertex(NULL),    
  fBetaTofVSp(NULL),
  hTOFSignalPion(NULL)
{
   // Default constructor (should not be used)
  fList1->SetName("results");
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
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL),
  fList1(new TList()),
  hNevent(NULL),
  hZvertex(NULL),    
  fBetaTofVSp(NULL),
  hTOFSignalPion(NULL)
{
  DefineOutput(1, TList::Class());
  fList1->SetName("results");
}
//_____________________________________________________________________________
AliAnalysisNucleiMass::~AliAnalysisNucleiMass()
{
  if(fList1) delete fList1;
}
//______________________________________________________________________________
void AliAnalysisNucleiMass::UserCreateOutputObjects()
{
  
  hNevent = new TH1F("hNevent","Centrality",20,0,100);

  hZvertex = new TH1F("hZvertex","Vertex distribution of accepted events; z vertex (cm)",120,-15,15);

  hTOFSignalPion = new TH1F("hTOFSignalPion","TOF signal 0.9<p_{T}<1.0; t-t_{exp}^{#pi} (ps)",1500,-1500,1500);
  
  fdEdxVSp[0] = new TH2F("fdEdxVSp","dE/dx vs p; p/|Z| (GeV/c); dE/dx_{TPC} (a.u.)",500,0,5,2000,0,1000);
  fdEdxVSp[1] = new TH2F("fdEdxVSp_pos","dE/dx vs p positive charge; p/|Z| (GeV/c); dE/dx_{TPC} (a.u.)",500,0,5,2000,0,1000);
  fdEdxVSp[2] = new TH2F("fdEdxVSp_neg","dE/dx vs p negative charge; p/|Z| (GeV/c); dE/dx_{TPC} (a.u.)",500,0,5,2000,0,1000);
  
  fBetaTofVSp = new TH2F("fBetaTofVSp","#beta_{TOF} vs p; p(GeV/c); #beta_{TOF}",1000,0,5,1300,0.4,1.05);
  
  fM2vsP_NoTpcCut[0] = new TH2F("fM2vsP_NoTpcCut","M_{TOF}^{2} vs p; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
  fM2vsP_NoTpcCut[1] = new TH2F("fM2vsP_NoTpcCut_Positive","M_{TOF}^{2} vs p Pos Part; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
  fM2vsP_NoTpcCut[2] = new TH2F("fM2vsP_NoTpcCut_Negative","M_{TOF}^{2} vs p Neg Part; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
  
  fM2vsP_NoTpcCut_DCAxyCut[0] = new TH2F("fM2vsP_NoTpcCut_DCAxycut","M_{TOF}^{2} vs p with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
  fM2vsP_NoTpcCut_DCAxyCut[1] = new TH2F("fM2vsP_NoTpcCut_Positive_DCAxycut","M_{TOF}^{2} vs p Pos Part with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);
  fM2vsP_NoTpcCut_DCAxyCut[2] = new TH2F("fM2vsP_NoTpcCut_Negative_DCAxycut","M_{TOF}^{2} vs p Neg Part with DCAxy cut; M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4}); p/|Z| (GeV/c)",8000,0,10,200,0,10);

  fM2vsZ[0] = new TH2F("fM2vsZ","M_{TOF}^{2} vs Z^{2} Integrated p_{T};Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  fM2vsZ[1] = new TH2F("fM2vsZ_0.3pT0.5","M_{TOF}^{2} vs Z^{2} 0.3<pT<0.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  fM2vsZ[2] = new TH2F("fM2vsZ_0.5pT1.0","M_{TOF}^{2} vs Z^{2} 0.5<pT<1.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  fM2vsZ[3] = new TH2F("fM2vsZ_1.0pT1.5","M_{TOF}^{2} vs Z^{2} 1.0<pT<1.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  fM2vsZ[4] = new TH2F("fM2vsZ_1.5pT2.0","M_{TOF}^{2} vs Z^{2} 1.5<pT<2.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  fM2vsZ[5] = new TH2F("fM2vsZ_2.0pT2.5","M_{TOF}^{2} vs Z^{2} 2.0<pT<2.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  fM2vsZ[6] = new TH2F("fM2vsZ_2.5pT3.0","M_{TOF}^{2} vs Z^{2} 2.5<pT<3.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  fM2vsZ[7] = new TH2F("fM2vsZ_3.0pT3.5","M_{TOF}^{2} vs Z^{2} 3.0<pT<3.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  fM2vsZ[8] = new TH2F("fM2vsZ_3.5pT4.0","M_{TOF}^{2} vs Z^{2} 3.5<pT<4.0;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  fM2vsZ[9] = new TH2F("fM2vsZ_4.0pT4.5","M_{TOF}^{2} vs Z^{2} 4.0<pT<4.5;Z;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4})",4000,-4,4,1000,0,10);
  
  
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
    fNsigmaTPC[i] = new TH2F(namePart_par_TPC[i],namePart_title_TPC[i],250,0,5,200,-5,5);
    fNsigmaTPC[i]->GetYaxis()->CenterTitle();
    fNsigmaTOF[i] = new TH2F(namePart_par_TOF[i],namePart_title_TOF[i],250,0,5,200,-5,5);
    fNsigmaTOF[i]->GetYaxis()->CenterTitle();
    hDeDxExp[i] = new TProfile(namePart_par_ProfileTPC[i],namePart_title_ProfileTPC[i],500,0,5,0,1000,"");
    hBetaExp[i] = new TProfile(namePart_par_ProfileTOF[i],namePart_title_ProfileTOF[i],400,0,5,0.4,1.05,"");
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
    snprintf(namePart_title_TPCvsP_kTOFtrue[iS],150,"NsigmaTPCvsP_kTOFout&&kTIME_%s;p (GeV/c);n_{#sigma_{TPC}}^{%s}",name[iS],name[iS]);
  }
  
  for(Int_t iS=0;iS<18;iS++) {
    fNsigmaTPCvsP_kTOFtrue[iS] = new TH2F(namePart_par_TPCvsP_kTOFtrue[iS],namePart_title_TPCvsP_kTOFtrue[iS],250,0,5,200,-5,5);
    fNsigmaTPCvsP_kTOFtrue[iS]->GetYaxis()->CenterTitle();
  }

  char name_par_MvsP[18][60];
  char name_title_MvsP[18][150];
  
  for(Int_t i=0;i<18;i++) {
    snprintf(name_par_MvsP[i],60,"M2vsP_%s",name[i]);
    snprintf(name_title_MvsP[i],150,"M_{TOF}^{2}_%s_2#sigma_TPCcut;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4});p/|Z| (GeV/c)",name[i]);
  }

  for (Int_t i=0;i<18;i++) fM2vsP[i] = new TH2F(name_par_MvsP[i],name_title_MvsP[i],8000,0,10,200,0,10);
  
  char name_par_MvsP_DCAxyCut[18][60];
  char name_title_MvsP_DCAxyCut[18][150];
  
  for(Int_t i=0;i<18;i++) {
    snprintf(name_par_MvsP_DCAxyCut[i],60,"M2vsP_DCAxyCut_%s",name[i]);
    snprintf(name_title_MvsP_DCAxyCut[i],150,"M_{TOF}^{2}_%s_2#sigma_TPCcut_DCAxyCut;M_{TOF}^{2}/Z^{2} (GeV^{2}/c^{4});p/|Z| (GeV/c)",name[i]);
  }

  for (Int_t i=0;i<18;i++) fM2vsP_DCAxyCut[i] = new TH2F(name_par_MvsP_DCAxyCut[i],name_title_MvsP_DCAxyCut[i],8000,0,10,200,0,10);


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

  Char_t nameM2CutDCAxy[18][nbin][120];
  Char_t titleM2CutDCAxy[18][nbin][120];
 
  Char_t nameM2CutGroundDCAxy[18][nbin][120];
  Char_t titleM2CutGroundDCAxy[18][nbin][120];

  
  for(Int_t iS=0;iS<18;iS++) {
    for(Int_t j=0;j<nbin;j++) {
      snprintf(nameDCAxy[iS][j],120,"hDCAxy_%s_%s",name[iS],par_name_nbin[j]);
      snprintf(titleDCAxy[iS][j],120,"hDCAxy_%s_%s;DCA_{xy} (cm)",name[iS],par_name_nbin[j]);
    
      snprintf(nameM2CutDCAxy[iS][j],120,"hM2_CutDCAxy_%s_%s",name[iS],par_name_nbin[j]);
      snprintf(titleM2CutDCAxy[iS][j],120,"hM2_CutDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin[j]);
      
      snprintf(nameM2CutGroundDCAxy[iS][j],120,"hM2_GroundCatDCAxy_%s_%s",name[iS],par_name_nbin[j]);
      snprintf(titleM2CutGroundDCAxy[iS][j],120,"hM2_GroundCatDCAxy_%s_%s;M^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],par_name_nbin[j]);
    }
  }
 
  for(Int_t iS=0;iS<18;iS++) {
    for(Int_t j=0;j<nbin;j++) {
      hDCAxy[iS][j] = new TH1D(nameDCAxy[iS][j],titleDCAxy[iS][j],1000,-2.5,2.5);//125 bins
      hDCAxy[iS][j]->GetXaxis()->CenterTitle();
      //hM2CutDCAxy[iS][j] = new TH1D(nameM2CutDCAxy[iS][j],titleM2CutDCAxy[iS][j],4000,0,10);
      //hM2CutGroundDCAxy[iS][j] = new TH1D(nameM2CutGroundDCAxy[iS][j],titleM2CutGroundDCAxy[iS][j],4000,0,10);
    }
  }
  
  //for e,#mu,#pi and antiparticle (e and #mu will not be drawn)
  //the binning is chosen for #pi distribution:
  for(Int_t iSp=0;iSp<3;iSp++) {
    for(Int_t j=0;j<nbin;j++) {
      hM2CutDCAxy[iSp][j] = new TH1D(nameM2CutDCAxy[iSp][j],titleM2CutDCAxy[iSp][j],600,-0.1,0.5);
      hM2CutGroundDCAxy[iSp][j] = new TH1D(nameM2CutGroundDCAxy[iSp][j],titleM2CutGroundDCAxy[iSp][j],600,-0.1,0.5);
      hM2CutDCAxy[iSp+9][j] = new TH1D(nameM2CutDCAxy[iSp+9][j],titleM2CutDCAxy[iSp+9][j],600,-0.1,0.5);
      hM2CutGroundDCAxy[iSp+9][j] = new TH1D(nameM2CutGroundDCAxy[iSp+9][j],titleM2CutGroundDCAxy[iSp+9][j],600,-0.1,0.5);
    }
  }

  for(Int_t j=0;j<nbin;j++) {
    hM2CutDCAxy[3][j] = new TH1D(nameM2CutDCAxy[3][j],titleM2CutDCAxy[3][j],400,0,1);
    hM2CutGroundDCAxy[3][j] = new TH1D(nameM2CutGroundDCAxy[3][j],titleM2CutGroundDCAxy[3][j],400,0,1);
    hM2CutDCAxy[3+9][j] = new TH1D(nameM2CutDCAxy[3+9][j],titleM2CutDCAxy[3+9][j],400,0,1);
    hM2CutGroundDCAxy[3+9][j] = new TH1D(nameM2CutGroundDCAxy[3+9][j],titleM2CutGroundDCAxy[3+9][j],400,0,1);
  }

  for(Int_t j=0;j<nbin;j++) {
    hM2CutDCAxy[4][j] = new TH1D(nameM2CutDCAxy[4][j],titleM2CutDCAxy[4][j],500,0,4);
    hM2CutGroundDCAxy[4][j] = new TH1D(nameM2CutGroundDCAxy[4][j],titleM2CutGroundDCAxy[4][j],500,0,4);
    hM2CutDCAxy[4+9][j] = new TH1D(nameM2CutDCAxy[4+9][j],titleM2CutDCAxy[4+9][j],500,0,4);
    hM2CutGroundDCAxy[4+9][j] = new TH1D(nameM2CutGroundDCAxy[4+9][j],titleM2CutGroundDCAxy[4+9][j],500,0,4);
  }

  for(Int_t j=0;j<nbin;j++) {
    hM2CutDCAxy[5][j] = new TH1D(nameM2CutDCAxy[5][j],titleM2CutDCAxy[5][j],500,0,6);
    hM2CutGroundDCAxy[5][j] = new TH1D(nameM2CutGroundDCAxy[5][j],titleM2CutGroundDCAxy[5][j],500,0,6);
    hM2CutDCAxy[5+9][j] = new TH1D(nameM2CutDCAxy[5+9][j],titleM2CutDCAxy[5+9][j],500,0,6);
    hM2CutGroundDCAxy[5+9][j] = new TH1D(nameM2CutGroundDCAxy[5+9][j],titleM2CutGroundDCAxy[5+9][j],500,0,6);
  }

  for(Int_t j=0;j<nbin;j++) {
    hM2CutDCAxy[6][j] = new TH1D(nameM2CutDCAxy[6][j],titleM2CutDCAxy[6][j],1000,0,12);
    hM2CutGroundDCAxy[6][j] = new TH1D(nameM2CutGroundDCAxy[6][j],titleM2CutGroundDCAxy[6][j],1000,0,12);
    hM2CutDCAxy[6+9][j] = new TH1D(nameM2CutDCAxy[6+9][j],titleM2CutDCAxy[6+9][j],1000,0,12);
    hM2CutGroundDCAxy[6+9][j] = new TH1D(nameM2CutGroundDCAxy[6+9][j],titleM2CutGroundDCAxy[6+9][j],1000,0,12);
  }

  for(Int_t j=0;j<nbin;j++) {
    hM2CutDCAxy[7][j] = new TH1D(nameM2CutDCAxy[7][j],titleM2CutDCAxy[7][j],200,0,4);
    hM2CutGroundDCAxy[7][j] = new TH1D(nameM2CutGroundDCAxy[7][j],titleM2CutGroundDCAxy[7][j],200,0,4);
    hM2CutDCAxy[7+9][j] = new TH1D(nameM2CutDCAxy[7+9][j],titleM2CutDCAxy[7+9][j],200,0,4);
    hM2CutGroundDCAxy[7+9][j] = new TH1D(nameM2CutGroundDCAxy[7+9][j],titleM2CutGroundDCAxy[7+9][j],200,0,4);
  }

  for(Int_t j=0;j<nbin;j++) {
    hM2CutDCAxy[8][j] = new TH1D(nameM2CutDCAxy[8][j],titleM2CutDCAxy[8][j],600,0,6);
    hM2CutGroundDCAxy[8][j] = new TH1D(nameM2CutGroundDCAxy[8][j],titleM2CutGroundDCAxy[8][j],600,0,6);
    hM2CutDCAxy[8+9][j] = new TH1D(nameM2CutDCAxy[8+9][j],titleM2CutDCAxy[8+9][j],600,0,6);
    hM2CutGroundDCAxy[8+9][j] = new TH1D(nameM2CutGroundDCAxy[8+9][j],titleM2CutGroundDCAxy[8+9][j],600,0,6);
  }

 
  fList1->Add(hNevent);
  fList1->Add(hZvertex);
  fList1->Add(hTOFSignalPion);
  for(Int_t iS=0;iS<18;iS++) fList1->Add(fNsigmaTPCvsP_kTOFtrue[iS]);
   
  for(Int_t i=0;i<3;i++) fList1->Add(fdEdxVSp[i]);
  for(Int_t i=0;i<9;i++) fList1->Add(hDeDxExp[i]);
  fList1->Add(fBetaTofVSp);
  for(Int_t i=0;i<9;i++) fList1->Add(hBetaExp[i]);
  for(Int_t i=0;i<9;i++) fList1->Add(fNsigmaTPC[i]);
  for(Int_t i=0;i<9;i++) fList1->Add(fNsigmaTOF[i]);

  for(Int_t i=0;i<10;i++) fList1->Add(fM2vsZ[i]);
  
  for(Int_t i=0;i<3;i++) fList1->Add(fM2vsP_NoTpcCut[i]);
  for(Int_t i=0;i<18;i++) fList1->Add(fM2vsP[i]);
  
  for(Int_t i=0;i<3;i++) fList1->Add(fM2vsP_NoTpcCut_DCAxyCut[i]);
  for(Int_t i=0;i<18;i++) fList1->Add(fM2vsP_DCAxyCut[i]);

  /*
  for(Int_t j=0;j<nbin;j++) {//electron
    fList1->Add(hDCAxy[0][j]);
    fList1->Add(hM2CutDCAxy[0][j]);
    fList1->Add(hM2CutGroundDCAxy[0][j]);
    fList1->Add(hDCAxy[9][j]);
    fList1->Add(hM2CutDCAxy[9][j]);
    fList1->Add(hM2CutGroundDCAxy[9][j]);
  }

  for(Int_t j=0;j<nbin;j++) {//muon
    fList1->Add(hDCAxy[1][j]);
    fList1->Add(hM2CutDCAxy[1][j]);
    fList1->Add(hM2CutGroundDCAxy[1][j]);
    fList1->Add(hDCAxy[10][j]);
    fList1->Add(hM2CutDCAxy[10][j]);
    fList1->Add(hM2CutGroundDCAxy[10][j]);
  }
  */

  for(Int_t j=0;j<nbin;j++) {//pion
    fList1->Add(hDCAxy[2][j]);
    fList1->Add(hM2CutDCAxy[2][j]);
    fList1->Add(hM2CutGroundDCAxy[2][j]);
    fList1->Add(hDCAxy[11][j]);
    fList1->Add(hM2CutDCAxy[11][j]);
    fList1->Add(hM2CutGroundDCAxy[11][j]);
  }

  for(Int_t j=0;j<nbin;j++) {//kaon
    fList1->Add(hDCAxy[3][j]);
    fList1->Add(hM2CutDCAxy[3][j]);
    fList1->Add(hM2CutGroundDCAxy[3][j]);
    fList1->Add(hDCAxy[12][j]);
    fList1->Add(hM2CutDCAxy[12][j]);
    fList1->Add(hM2CutGroundDCAxy[12][j]);
  }

  for(Int_t j=0;j<nbin;j++) {//proton
    fList1->Add(hDCAxy[4][j]);
    fList1->Add(hM2CutDCAxy[4][j]);
    fList1->Add(hM2CutGroundDCAxy[4][j]);
    fList1->Add(hDCAxy[13][j]);
    fList1->Add(hM2CutDCAxy[13][j]);
    fList1->Add(hM2CutGroundDCAxy[13][j]);
  }

  for(Int_t j=0;j<nbin;j++) {//deuteron
    fList1->Add(hDCAxy[5][j]);
    fList1->Add(hM2CutDCAxy[5][j]);
    fList1->Add(hM2CutGroundDCAxy[5][j]);
    fList1->Add(hDCAxy[14][j]);
    fList1->Add(hM2CutDCAxy[14][j]);
    fList1->Add(hM2CutGroundDCAxy[14][j]);
  }

  /*
  for(Int_t j=0;j<nbin;j++) {//triton
    fList1->Add(hDCAxy[6][j]);
    fList1->Add(hM2CutDCAxy[6][j]);
    fList1->Add(hM2CutGroundDCAxy[6][j]);
    fList1->Add(hDCAxy[15][j]);
    fList1->Add(hM2CutDCAxy[15][j]);
    fList1->Add(hM2CutGroundDCAxy[15][j]);
  }
  */

  for(Int_t j=0;j<nbin;j++) {//He3
    fList1->Add(hDCAxy[7][j]);
    fList1->Add(hM2CutDCAxy[7][j]);
    fList1->Add(hM2CutGroundDCAxy[7][j]);
    fList1->Add(hDCAxy[16][j]);
    fList1->Add(hM2CutDCAxy[16][j]);
    fList1->Add(hM2CutGroundDCAxy[16][j]);
  }
  
  /*
  for(Int_t j=0;j<nbin;j++) {//4He
    fList1->Add(hDCAxy[8][j]);
    fList1->Add(hM2CutDCAxy[8][j]);
    fList1->Add(hM2CutGroundDCAxy[8][j]);
    fList1->Add(hDCAxy[17][j]);
    fList1->Add(hM2CutDCAxy[17][j]);
    fList1->Add(hM2CutGroundDCAxy[17][j]);
  }
  */

  //fTrackFilter = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  
  // Post output data.
  PostData(1, fList1);
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
    this->Dump();
    return;
  }
  
  if(fESD) fEvent = fESD;
  else fEvent = fAOD;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse(); // data member di tipo "const AliPIDResponse *fPIDResponse;"

  //Centrality
  Float_t v0Centr  = -10.;
  Float_t trkCentr  = -10.;
  AliCentrality *centrality = fEvent->GetCentrality();
  if (centrality){
    v0Centr  = centrality->GetCentralityPercentile("V0M"); // VZERO
    trkCentr = centrality->GetCentralityPercentile("TRK"); // TPC
  }

  const AliAODVertex* vtxAOD = fAOD->GetPrimaryVertex();
  
  Float_t zvtx = 10000.0;

  if(vtxAOD->GetNContributors()>0)
    zvtx = vtxAOD->GetZ();
  
  if(TMath::Abs(v0Centr - trkCentr) < 5.0 && TMath::Abs(zvtx) < 10.0){ // consistency cut on centrality selection AND d(zPrimaryVertez;NominalPointInteraction)<10cm
    
    //Bool_t isTrack=1;
    
    Int_t nTracks = fEvent->GetNumberOfTracks();
    
    if(v0Centr>=fCentrality[0] && v0Centr<=fCentrality[1]) {//window cut centrality open
      
      hZvertex->Fill(zvtx);
      hNevent->Fill(v0Centr);
      
      for(Int_t iT = 0; iT < nTracks; iT++) { // loop on the tracks
	AliVTrack* track = (AliVTrack *) fEvent->GetTrack(iT);
	
	if (!track){
	  track->Delete();
	  //isTrack=0;
	  continue;
	}
	
	Bool_t trkFlag = 0;
	
	if(fAOD)
	  trkFlag = ((AliAODTrack *) track)->TestFilterBit(FilterBit);
	
	  //TestFilterBit(16) -- Standard Cuts with very loose DCA: GetStandardITSTPCTrackCuts2011(kFALSE) && SetMaxDCAToVertexXY(2.4) && SetMaxDCAToVertexZ(3.2) && SetDCaToVertex2D(kTRUE)
	
	  //TestFilterBit(32) (STARDARD) -- Standard Cuts with very tight DCA cut ( 7sigma^primaries: 7*(0.0015+0.0050/pt^1.1) ) : GetStandardITSTPCTrackCuts2011(). 

	else{
	  //trkFlag = fTrackFilter->IsSelected(((AliESDtrack *) track));
	}
	
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
	Float_t DCAxy = b[1];
	Float_t DCAz = b[0];

	if(TMath::Abs(DCAz)>DCAzCUT)//CUT ON DCAz
	  continue;
	
	Bool_t kTpcPure;
	kTpcPure = track->GetTPCsignal()>10;
	
	kTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);

       	Float_t nsigmaTPC[9];
	Float_t nsigmaTOF[9];
	
	Double_t expdedx[9];
//	Double_t dedxRes[9];
	
	for(Int_t iS=0;iS < 9;iS++){ //TPC expected signal
	  //expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(pTPC, (AliPID::EParticleType) iS);// Deprecated function: Temporary solution to measure also He (The other method has a small bag) . This must be replace from:
	  expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
//	  dedxRes[iS] = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
	  //if(iS>6) dedxRes[iS] = expdedx[iS]*0.08; //Temporary solution to measure also He. This line must be removed 
	}
	
	if(kTpcPure) {
	  for(Int_t iS=0;iS < 9;iS++){
	    nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
	    //if(iS>6) nsigmaTPC[iS] = (dedx-expdedx[iS])/dedxRes[iS]; // Temporary solution to measure also He. This line must be removed.
	    fNsigmaTPC[iS]->Fill(pt,nsigmaTPC[iS]);
	    hDeDxExp[iS]->Fill(pTPC,expdedx[iS]);
	  }
	  fdEdxVSp[0]->Fill(pTPC,dedx);
	  if(charge>0) fdEdxVSp[1]->Fill(pTPC,dedx);
	  else fdEdxVSp[2]->Fill(pTPC,dedx);
	}
	
	Float_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.877837,2.817402,1.408701,1.877837};
	
	Double_t exptimes[9]; // TOF expected times
	track->GetIntegratedTimes(exptimes);
	exptimes[5] = exptimes[0] / p * massOverZ[5] * TMath::Sqrt(1+p*p/massOverZ[5]/massOverZ[5]);
	exptimes[6] = exptimes[0] / p * massOverZ[6] * TMath::Sqrt(1+p*p/massOverZ[6]/massOverZ[6]);
	exptimes[7] = exptimes[0] / p * massOverZ[7] * TMath::Sqrt(1+p*p/massOverZ[7]/massOverZ[7]);
	exptimes[8] = exptimes[0] / p * massOverZ[8] * TMath::Sqrt(1+p*p/massOverZ[8]/massOverZ[8]);
	
//	Double_t tofRes[9] = {1000,1000,1000,1000,1000,1000,1000,1000,1000};
//	if(kTOF){
//	  for(Int_t iS=0;iS < 9;iS++) tofRes[iS] = fPIDResponse->GetTOFResponse().GetExpectedSigma(p, exptimes[iS], massOverZ[iS]);
//	}
	
	beta=exptimes[0];
	
	if(kTOF) {
	  for(Int_t iS=0;iS < 9;iS++){
	    nsigmaTOF[iS] = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType) iS);
	    fNsigmaTOF[iS]->Fill(pt,nsigmaTOF[iS]);
	    hBetaExp[iS]->Fill(p,beta/exptimes[iS]);
	  }
	  if(pt>0.9 && pt<1.0) hTOFSignalPion->Fill(tof-exptimes[2]);
	  beta=beta/tof;
	  fBetaTofVSp->Fill(p,beta);
	}
	
	Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,#mu,#pi,K,p,d,t,3He,4He
	Int_t FlagPid = 0;
	Float_t binPt[nbin+1];

	for(Int_t i=0;i<nbin+1;i++) {
	  binPt[i]=0.4+i*0.1;
	}
	
	Float_t binCutPt[10] = {0.3,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5};

	if(kTOF && kTpcPure) {
	  
	  M2 = (p*p*(1-beta*beta))/(beta*beta);
	  
	  fM2vsP_NoTpcCut[0]->Fill(M2,p);
	  if(TMath::Abs(DCAxy)<DCAxyCUT) fM2vsP_NoTpcCut_DCAxyCut[0]->Fill(M2,p);
	  
	  if(M2>0.0) {
	    M=TMath::Sqrt(M2);
	    Z2 = TMath::Power(dedx/fPIDResponse->GetTPCResponse().GetExpectedSignal(pTPC*massOverZ[4]/M, AliPID::kProton),0.862);
	    fM2vsZ[0]->Fill(charge*TMath::Sqrt(Z2),M2);
	    
	    for(Int_t i=0;i<9;i++) {
	      if(pt>binCutPt[i] && pt<binCutPt[i+1]){
		fM2vsZ[i+1]->Fill(charge*TMath::Sqrt(Z2),M2);
		break;
	      }
	    }
	  }
	  
	  if(charge>0) {
	    fM2vsP_NoTpcCut[1]->Fill(M2,p);
	    if(TMath::Abs(DCAxy)<DCAxyCUT) fM2vsP_NoTpcCut_DCAxyCut[1]->Fill(M2,p);
	    for(Int_t iS=0;iS < 9;iS++){
	      fNsigmaTPCvsP_kTOFtrue[iS]->Fill(p,nsigmaTPC[iS]);
	    }
	  }
	  else {//else charge<0
	    fM2vsP_NoTpcCut[2]->Fill(M2,p);
	    if(TMath::Abs(DCAxy)<DCAxyCUT) fM2vsP_NoTpcCut_DCAxyCut[2]->Fill(M2,p);
	    
	    for(Int_t iS=0;iS < 9;iS++){
	      fNsigmaTPCvsP_kTOFtrue[iS+9]->Fill(p,nsigmaTPC[iS]);
	    }
	  }
	  
	  for(Int_t iS=0;iS<9;iS++) {
	    if(TMath::Abs(nsigmaTPC[iS])<2) {
	      FlagPid += ((Int_t)TMath::Power(2,iS));
	    }
	  }
	    
	  for(Int_t iS=0;iS<9;iS++) {
	    if(FlagPid & stdFlagPid[iS] || !kTPCcut) {
	      if(charge>0) {
		fM2vsP[iS]->Fill(M2,p);
		for(Int_t j=0;j<nbin;j++) {
		  if(pt>binPt[j] && pt<binPt[j+1]) {
		    hDCAxy[iS][j]->Fill(DCAxy);
		    hDCAxy[iS][j]->Fill(-DCAxy);
		    if(TMath::Abs(DCAxy)<DCAxyCUT) {
		      hM2CutDCAxy[iS][j]->Fill(M2);
		      fM2vsP_DCAxyCut[iS]->Fill(M2,p);
		    }
		    if(TMath::Abs(DCAxy+0.5)<DCAxyCUT) hM2CutGroundDCAxy[iS][j]->Fill(M2);
		    break;
		  }
		}
	      }
	      else {//if(charge<0)
		fM2vsP[iS+9]->Fill(M2,p);
		for(Int_t j=0;j<nbin;j++) {
		  if(pt>binPt[j] && pt<binPt[j+1]) {
		    hDCAxy[iS+9][j]->Fill(DCAxy);
		    hDCAxy[iS+9][j]->Fill(-DCAxy);
		    if(TMath::Abs(DCAxy)<DCAxyCUT) {
		      hM2CutDCAxy[iS+9][j]->Fill(M2);
		      fM2vsP_DCAxyCut[iS]->Fill(M2,p);
		    }
		    if(TMath::Abs(DCAxy+0.5)<DCAxyCUT) hM2CutGroundDCAxy[iS+9][j]->Fill(M2);
		    break;
		  }
		}
	      }
	    }
	  }
	}// close (KTOF && kTpcPure request)
	
	
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
