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
#include "TF2.h"
#include "TGraph.h"
#include "TProfile.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisManager.h"
#include "TFile.h"

ClassImp(AliAnalysisNucleiMass)

//_____________________________________________________________________________
AliAnalysisNucleiMass::AliAnalysisNucleiMass():
  AliAnalysisTaskSE(),
  Centrality(),                         
  FilterBit(16),                                
  EtaLimit(),                           
  DCAxyCut(0.1),                               
  DCAzCut(1000.),                                
  NsigmaTpcCut(2.0),                           
  NminTpcCluster(0),                           
  iTrdCut(0),
  kSignalCheck(1),
  iMtof(1),
  kPvtxCorr(1),
  iBconf(0),                                  
  kTOF(0),               
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL)
{
  Centrality[0]=0.0;
  Centrality[1]=100.0;

  EtaLimit[0]=-99.0;
  EtaLimit[1]=99.0;

  fList[0]=new TList();
  fList[0]->SetName("results");
  
  fList[1]=new TList();
  fList[1]->SetName("results2");
}
//______________________________________________________________________________
AliAnalysisNucleiMass::AliAnalysisNucleiMass(const char *name):
  AliAnalysisTaskSE(name),
  Centrality(),                         
  FilterBit(16),                                
  EtaLimit(),                           
  DCAxyCut(0.1),                               
  DCAzCut(1000.),                                
  NsigmaTpcCut(2.0),                           
  NminTpcCluster(0),
  iTrdCut(0),
  kSignalCheck(1),
  iMtof(1),
  kPvtxCorr(1),
  iBconf(0),                                  
  kTOF(0),               
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL)
{

  Centrality[0]=0.0;
  Centrality[1]=100.0;

  EtaLimit[0]=-99.0;
  EtaLimit[1]=99.0;

  fList[0]=new TList();
  DefineOutput(1, TList::Class());
  fList[0]->SetName("results");
  
  fList[1]=new TList();
  DefineOutput(2, TList::Class());
  fList[1]->SetName("results2");
}
//_____________________________________________________________________________
AliAnalysisNucleiMass::~AliAnalysisNucleiMass()
{
  if(fList[0]) delete fList[0];
  if(fList[1]) delete fList[1];
}
//______________________________________________________________________________
void AliAnalysisNucleiMass::UserCreateOutputObjects()
{
  Char_t namePart[nPart][30];
  snprintf(namePart[0],30,"e");
  snprintf(namePart[1],30,"#mu");
  snprintf(namePart[2],30,"#pi");
  snprintf(namePart[3],30,"K");
  snprintf(namePart[4],30,"p");
  snprintf(namePart[5],30,"d");
  snprintf(namePart[6],30,"t");
  snprintf(namePart[7],30,"He3");
  snprintf(namePart[8],30,"He4");
  
  Char_t name[nSpec][30];
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
  
  Double_t binP[nbin+1];
  for(Int_t i=0;i<nbin+1;i++) {
    binP[i]=0.4+0.1*i;
  }
  
  Char_t name_nbin[nbin][200];
  for(Int_t j=0;j<nbin;j++) {
    snprintf(name_nbin[j],200,"%.1f<P<%.1f",binP[j],binP[j+1]);
  }
  
  for(Int_t iB=0;iB<nBconf;iB++) {
    
    htemp[iB] = new TH1F("htemp","htemp (avoid the problem with the empty list...);B field",20,-10,10);

    hCentrality[iB][0] = new TH1F("hCentrality_Selected","Centrality (selected events);centrality(%)",20,0,100);
    hCentrality[iB][1] = new TH1F("hCentrality_Analyzed","Centrality (analyzed events);centrality (%)",20,0,100);
    
    hZvertex[iB][0] = new TH1F("hZvertex_Selected","Vertex distribution of selected events;z vertex (cm)",240,-30,30);
    hZvertex[iB][1] = new TH1F("hZvertex_Analyzed","Vertex distribution of analyzed events;z vertex (cm)",240,-30,30);
    
    hEta[iB] = new TH1F("hEta_Analyzed","|#eta| distribution after the track cuts;|#eta|",200,-1.0,1.0);
    
    hPhi[iB] = new TH1F("hPhi_Analyzed","#phi distribution after the track cuts;#phi (rad.)",90,0,6.3);//Each TRD supermodule is divided for 5 (DeltaPhi(TRD)=0.35 theoretical)
    
    Int_t hbins[2];
    if(kSignalCheck!=0) {hbins[0]=1; hbins[1]=1;}//{hbins[0]=100; hbins[1]=90;} to reduce RAM consuming (toram)
    else {hbins[0]=1; hbins[1]=1;}
    fEtaPhi[iB] = new TH2F("fEtaPhi_Analyzed","#eta vs. #phi after the track cuts;|#eta|;#phi (rad.)",hbins[0],-1.0,1.0,hbins[1],0,6.3);

    hNTpcCluster[iB] = new TH1F("hNTpcCluster","Number of the TPC clusters after the track cuts;n_{cl}^{TPC}",300,0,300);

    hNTrdSlices[iB] = new TH1F("hNTrdSlices","Number of the TRD slices after the track cuts;n_{slices}^{TRD}",40,0,40);

    if(kSignalCheck==1) {hbins[0]=500; hbins[1]=2000;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=1; hbins[1]=1;}//{hbins[0]=100; hbins[1]=500;} toram
    fdEdxVSp[iB][0] = new TH2F("fdEdxVSp_pos","dE/dx vs p (positive charge); p/|z| (GeV/c); dE/dx_{TPC} (a.u.)",hbins[0],0,5,hbins[1],0,1000);
    fdEdxVSp[iB][1] = new TH2F("fdEdxVSp_neg","dE/dx vs p (negative charge); p/|z| (GeV/c); dE/dx_{TPC} (a.u.)",hbins[0],0,5,hbins[1],0,1000);

    Char_t name_hDeDxExp[nPart][200];
    Char_t title_hDeDxExp[nPart][200];
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_hDeDxExp[i],200,"hDeDxExp_%s",namePart[i]);
      snprintf(title_hDeDxExp[i],200,"Expected dE/dx of %s in the TPC;p/|z| (GeV/c);dE/dx_{TPC} (a.u.)",namePart[i]);
      hDeDxExp[iB][i] = new TProfile(name_hDeDxExp[i],title_hDeDxExp[i],1,0,5,0,1,"");//,500,0,5,0,1000,""); toram
    }

    Char_t name_fNsigmaTpc[nPart][200];
    Char_t title_fNsigmaTpc[nPart][200];
    if(kSignalCheck==1) {hbins[0]=1; hbins[1]=1;}//{hbins[0]=100; hbins[1]=100;} toram
    else {hbins[0]=1; hbins[1]=1;}
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_fNsigmaTpc[i],200,"NsigmaTpc_%s",namePart[i]);
      snprintf(title_fNsigmaTpc[i],200,"NsigmaTpc_%s;p_{TPC}/|z| (GeV/c);n_{#sigma_{TPC}}^{%s}",namePart[i],namePart[i]);
      fNsigmaTpc[iB][i] = new TH2F(name_fNsigmaTpc[i],title_fNsigmaTpc[i],hbins[0],0,5,hbins[1],-5,5);
    }
    
    if(kSignalCheck>0) {hbins[0]=1; hbins[1]=1;}//{hbins[0]=100; hbins[1]=100;} toram
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    Char_t name_fNsigmaTpc_kTOF[nSpec][200];
    Char_t title_fNsigmaTpc_kTOF[nSpec][200];
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fNsigmaTpc_kTOF[i],200,"NsigmaTpc_%s_kTOF",name[i]);
      snprintf(title_fNsigmaTpc_kTOF[i],200,"NsigmaTpc_kTOF_%s in DCAxyCut;p/|z| (GeV/c);n_{#sigma_{TPC}}^{%s}",name[i],name[i]);
      fNsigmaTpc_kTOF[iB][i] = new TH2F(name_fNsigmaTpc_kTOF[i],title_fNsigmaTpc_kTOF[i],hbins[0],0,5,hbins[1],-5,5);
    }

    if(kSignalCheck==1) {hbins[0]=1000; hbins[1]=1300;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=1; hbins[1]=1;}//{hbins[0]=100; hbins[1]=260;}
    fBetaTofVSp[iB][0] = new TH2F("fBetaTofVSp_pos","#beta_{TOF} vs p/|z| (positive charge);p(GeV/c);#beta_{TOF}",hbins[0],0,5,hbins[1],0.4,1.05);
    fBetaTofVSp[iB][1] = new TH2F("fBetaTofVSp_neg","#beta_{TOF} vs p/|z| (negative charge);p(GeV/c);#beta_{TOF}",hbins[0],0,5,hbins[1],0.4,1.05);
    
    Char_t name_hBetaExp[nPart][200];
    Char_t title_hBetaExp[nPart][200];
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_hBetaExp[i],200,"hBetaTofVsP_Exp_%s",namePart[i]);
      snprintf(title_hBetaExp[i],200,"Expected #beta_{TOF} vs p/|z| of %s;p/|z| (GeV/c); #beta_{TOF}",namePart[i]);
      hBetaExp[iB][i] = new TProfile(name_hBetaExp[i],title_hBetaExp[i],1,0,5,0.4,1.05,"");//,400,0,5,0.4,1.05,""); toram
    }
    
    Char_t name_fNsigmaTof[nPart][200];
    Char_t title_fNsigmaTof[nPart][200];    
    if(kSignalCheck==1) {hbins[0]=100; hbins[1]=100;}
    else {hbins[0]=1; hbins[1]=1;}
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_fNsigmaTof[i],200,"NsigmaTof_%s",namePart[i]);
      snprintf(title_fNsigmaTof[i],200,"NsigmaTof_%s;p_{T}/|z| (GeV/c);n_{#sigma_{TOF}}^{%s}",namePart[i],namePart[i]);
      fNsigmaTof[iB][i] = new TH2F(name_fNsigmaTof[i],title_fNsigmaTof[i],hbins[0],0,5,hbins[1],-5,5);
    }

    Char_t name_fNsigmaTof_DcaCut[nSpec][200];
    Char_t title_fNsigmaTof_DcaCut[nSpec][200];    
    if(kSignalCheck==1) {hbins[0]=100; hbins[1]=100;}
    else {hbins[0]=1; hbins[1]=1;}
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fNsigmaTof_DcaCut[i],200,"NsigmaTof_DcaCut_%s",name[i]);
      snprintf(title_fNsigmaTof_DcaCut[i],200,"NsigmaTof_%s with DCAxyCut;p_{T}/|z| (GeV/c);n_{#sigma_{TOF}}^{%s}",name[i],name[i]);
      fNsigmaTof_DcaCut[iB][i] = new TH2F(name_fNsigmaTof_DcaCut[i],title_fNsigmaTof_DcaCut[i],hbins[0],0,5,hbins[1],-5,5);
    }

    if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
    else {hbins[0]=1; hbins[1]=1;}
    fM2vsP_NoTpcCut[iB][0][0] = new TH2F("fM2vsP_NoTpcCut_pos","m^{2}/z^{2}_{TOF} vs p/|z| (positive charge);m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,5);
    fM2vsP_NoTpcCut[iB][0][1] = new TH2F("fM2vsP_NoTpcCut_neg","m^{2}/z^{2}_{TOF} vs p/|z| (negative charge);m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,5);

    if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=1; hbins[1]=1;}// {hbins[0]=1000; hbins[1]=100;} toram
    fM2vsP_NoTpcCut[iB][1][0] = new TH2F("fM2vsP_NoTpcCut_DCAxyCut_pos","m^{2}/z^{2}_{TOF} vs p/|z| (positive charge) with DCAxy cut;m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,5);
    fM2vsP_NoTpcCut[iB][1][1] = new TH2F("fM2vsP_NoTpcCut_DCAxyCut_neg","m^{2}/z^{2}_{TOF} vs p/|z| (negative charge) with DCAxy cut;m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,5);

    Char_t name_fM2vsP[2][18][300]; 
    Char_t title_fM2vsP[2][18][300]; 

    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fM2vsP[0][i],300,"fM2vsPc_%s",name[i]);
      snprintf(title_fM2vsP[0][i],300,"m^{2}/z^{2}_{TOF} vs p/|z| of %s with a NsigmaTpcCut (pReco->pTrue for nuclei);m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p/|z| (GeV/c)",name[i]);
      
      snprintf(name_fM2vsP[1][i],300,"fM2vsPc_%s_DCAxyCut",name[i]);
      snprintf(title_fM2vsP[1][i],300,"m^{2}/z^{2}_{TOF} vs p/|z| of %s with a NsigmaTpcCut and with the DCAxy cut (pReco->pTrue for nuclei);m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p/|z| (GeV/c)",name[i]);

      if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
      else {hbins[0]=1; hbins[1]=1;}
      fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],hbins[0],0,10,hbins[1],0,5);
      
      if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
      else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
      else if(kSignalCheck==2) {hbins[0]=1; hbins[1]=1;}//{hbins[0]=1000; hbins[1]=100;} toram
      fM2vsP[iB][1][i] = new TH2F(name_fM2vsP[1][i],title_fM2vsP[1][i],hbins[0],0,10,hbins[1],0,5);
    }
    
    if(kSignalCheck==1) {hbins[0]=4000; hbins[1]=1000;}
    else {hbins[0]=1; hbins[1]=1;}
    fM2vsZ[iB][0] = new TH2F("fM2vsZ","m^{2}/z^{2}_{TOF} vs z_{TPC} Integrated p_{T};z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][1] = new TH2F("fM2vsZ_0.5pT1.0","m^{2}/z^{2}_{TOF} vs z_{TPC} 0.5<pT<1.0;z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][2] = new TH2F("fM2vsZ_1.0pT1.5","m^{2}/z^{2}_{TOF} vs z_{TPC} 1.0<pT<1.5;z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][3] = new TH2F("fM2vsZ_1.5pT2.0","m^{2}/z^{2}_{TOF} vs z_{TPC} 1.5<pT<2.0;z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][4] = new TH2F("fM2vsZ_2.0pT2.5","m^{2}/z^{2}_{TOF} vs z_{TPC} 2.0<pT<2.5;z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][5] = new TH2F("fM2vsZ_2.5pT3.0","m^{2}/z^{2}_{TOF} vs z_{TPC} 2.5<pT<3.0;z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][6] = new TH2F("fM2vsZ_3.0pT3.5","m^{2}/z^{2}_{TOF} vs z_{TPC} 3.0<pT<3.5;z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][7] = new TH2F("fM2vsZ_3.5pT4.0","m^{2}/z^{2}_{TOF} vs z_{TPC} 3.5<pT<4.0;z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][8] = new TH2F("fM2vsZ_4.0pT4.5","m^{2}/z^{2}_{TOF} vs z_{TPC} 4.0<pT<4.5;z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][9] = new TH2F("fM2vsZ_4.5pT5.0","m^{2}/z^{2}_{TOF} vs z_{TPC} 2.0<pT<2.5;z_{TPC};m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
        
    Char_t name_hDCAxy[18][nbin][200];
    Char_t title_hDCAxy[18][nbin][200];
    Char_t name_hDCAz[18][nbin][200];
    Char_t title_hDCAz[18][nbin][200];
    for(Int_t iS=0;iS<nSpec;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	snprintf(name_hDCAxy[iS][j],200,"hDCAxy_%s_%s",name[iS],name_nbin[j]);
	snprintf(title_hDCAxy[iS][j],200,"hDCAxy_%s_%s;DCA_{xy} (cm)",name[iS],name_nbin[j]);
	hDCAxy[iB][iS][j] = new TH1D(name_hDCAxy[iS][j],title_hDCAxy[iS][j],875,-3.5,3.5);

	snprintf(name_hDCAz[iS][j],200,"hDCAz_%s_%s",name[iS],name_nbin[j]);
	snprintf(title_hDCAz[iS][j],200,"hDCAz_%s_%s;DCA_{z} (cm)",name[iS],name_nbin[j]);
	hDCAz[iB][iS][j] = new TH1D(name_hDCAz[iS][j],title_hDCAz[iS][j],875,-3.5,3.5);
      }
    }
  
    Char_t name_hM2CutDCAxy[18][nbin][200];
    Char_t title_hM2CutDCAxy[18][nbin][200];
    Char_t name_hM2CutGroundDCAxy[18][nbin][200];
    Char_t title_hM2CutGroundDCAxy[18][nbin][200];
    for(Int_t iS=0;iS<nSpec;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	snprintf(name_hM2CutDCAxy[iS][j],200,"hM2_CutDCAxy_%s_%s",name[iS],name_nbin[j]);
	snprintf(title_hM2CutDCAxy[iS][j],200,"m^{2}/z^{2} Tof distribution of %s in DCAxy cut and in %s;m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],name_nbin[j]);
	snprintf(name_hM2CutGroundDCAxy[iS][j],200,"hM2_GroundCatDCAxy_%s_%s",name[iS],name_nbin[j]);
	snprintf(title_hM2CutGroundDCAxy[iS][j],200,"m^{2}/z^{2} Tof distribution of %s in the bkg. of DCAxy and in %s;m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],name_nbin[j]);
      }
    }

    const Int_t BinM2pT[nPart]={1,1,600,250,500,500,1000,400,600};
    const Double_t RangeM2min[nPart]={0.0,0.0,-0.1,0.0,0.0,0.0,0.0,0.0,0.0};
    const Double_t RangeM2max[nPart]={1.0,1.0,0.5,2.0,4.0,6.0,12.0,4.0,6.0};

    for(Int_t iS=0;iS<nPart;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	
	hM2CutDCAxy[iB][iS][j] = new TH1D(name_hM2CutDCAxy[iS][j],title_hM2CutDCAxy[iS][j],BinM2pT[iS],RangeM2min[iS],RangeM2max[iS]);
	hM2CutGroundDCAxy[iB][iS][j] = new TH1D(name_hM2CutGroundDCAxy[iS][j],title_hM2CutGroundDCAxy[iS][j],BinM2pT[iS],RangeM2min[iS],RangeM2max[iS]);
      
	hM2CutDCAxy[iB][iS+nPart][j] = new TH1D(name_hM2CutDCAxy[iS+nPart][j],title_hM2CutDCAxy[iS+nPart][j],BinM2pT[iS],RangeM2min[iS],RangeM2max[iS]);
	hM2CutGroundDCAxy[iB][iS+nPart][j] = new TH1D(name_hM2CutGroundDCAxy[iS+nPart][j],title_hM2CutGroundDCAxy[iS+nPart][j],BinM2pT[iS],RangeM2min[iS],RangeM2max[iS]);
      }
    }

    Char_t name_fPmeanVsBetaGamma[18][200];
    Char_t title_fPmeanVsBetaGamma[18][200];
    
    hbins[0]=200; hbins[1]=200;
    for(Int_t iS=0;iS<nSpec;iS++) {
      snprintf(name_fPmeanVsBetaGamma[iS],200,"fPmeanVsPvtx_%s",name[iS]);
      snprintf(title_fPmeanVsBetaGamma[iS],200,"<p>/p_{vtx} vs #beta#gamma of %s (in DCAxyCut);p_{vtx}/m_{%s};<p>_{%s}/p_{vtx}",name[iS],name[iS],name[iS]);
      fPmeanVsBetaGamma[iB][iS]=new TH2F(name_fPmeanVsBetaGamma[iS],title_fPmeanVsBetaGamma[iS],hbins[0],0,10,hbins[1],0.8,1.2);
    }	
    
    Char_t name_prPmeanVsBetaGamma[18][200];
    Char_t title_prPmeanVsBetaGamma[18][200];
    
    for(Int_t iS=0;iS<nSpec;iS++) {
      snprintf(name_prPmeanVsBetaGamma[iS],200,"prPmeanVsPvtx_%s",name[iS]);
      snprintf(title_prPmeanVsBetaGamma[iS],200,"<p>/p_{vtx} vs #beta#gamma of %s (in DCAxyCut);p_{vtx}/m_{%s};<p>_{%s}/p_{vtx}",name[iS],name[iS],name[iS]);
      prPmeanVsBetaGamma[iB][iS]=new TProfile(name_prPmeanVsBetaGamma[iS],title_prPmeanVsBetaGamma[iS],hbins[0],0,10,0.8,1.2,"");
    }	
    
    //for (bar)d
    fPvtxTrueVsReco[0]=new TF2("fcorr_d","([0]*TMath::Power(x,[1])+[2])*(TMath::Power((TMath::Exp([3]*x)+[4]),[5]*TMath::Power(y,[6])));|#eta|;p_{true}/p_{reco}",0.0001,100,0,1);//for (bar)d
    fPvtxTrueVsReco[0]->SetParameter(0,0.031263);
    fPvtxTrueVsReco[0]->SetParameter(1,-3.276770);
    fPvtxTrueVsReco[0]->SetParameter(2,1.000113);
    fPvtxTrueVsReco[0]->SetParameter(3,-5.195875);
    fPvtxTrueVsReco[0]->SetParameter(4,1.000674);
    fPvtxTrueVsReco[0]->SetParameter(5,2.870503);
    fPvtxTrueVsReco[0]->SetParameter(6,3.777729);
    
    fPvtxTrueVsReco[0]->SetNpx(fPvtxTrueVsReco[0]->GetNpx()*10);
        
    //for (bar)He3
    fPvtxTrueVsReco[1]=new TF2("fcorr_He","([0]*TMath::Power(x,[1])+[2])*(TMath::Power((TMath::Exp([3]*x)+[4]),[5]*TMath::Power(y,[6])));|#eta|;p_{true}/p_{reco}",0.0001,100,0,1);//for (bar)He3
    fPvtxTrueVsReco[1]->SetParameter(0,0.037986);
    fPvtxTrueVsReco[1]->SetParameter(1,-2.707620);
    fPvtxTrueVsReco[1]->SetParameter(2,1.000742);
    fPvtxTrueVsReco[1]->SetParameter(3,-4.934743);
    fPvtxTrueVsReco[1]->SetParameter(4,1.001640);
    fPvtxTrueVsReco[1]->SetParameter(5,2.744372);
    fPvtxTrueVsReco[1]->SetParameter(6,3.528561);

    fPvtxTrueVsReco[1]->SetNpx(fPvtxTrueVsReco[1]->GetNpx()*10);

    prPvtxTrueVsReco[iB][0]=new TProfile("prPvtxTrueVsReco_d","p_{true} vs p_{reco} of d and dbar;p_{reco} (GeV/c); p_{true}/p_{reco} (d)",200,0,10);
    prPvtxTrueVsReco[iB][1]=new TProfile("prPvtxTrueVsReco_He3","p_{true} vs p_{reco} of He3 and He3bar;p_{reco} (GeV/c);p_{true}/p_{reco} (He3)",200,0,10);

    Char_t nameTemp[10][200];
    snprintf(nameTemp[0],200,"#pi^{+}");
    snprintf(nameTemp[1],200,"K^{+}");
    snprintf(nameTemp[2],200,"p");
    snprintf(nameTemp[3],200,"d");
    snprintf(nameTemp[4],200,"He3");
    snprintf(nameTemp[5],200,"#pi^{-}");
    snprintf(nameTemp[6],200,"K^{-}");
    snprintf(nameTemp[7],200,"#bar{p}");
    snprintf(nameTemp[8],200,"#bar{d}");
    snprintf(nameTemp[9],200,"#bar{He3}");
    
    Char_t name_fPmeanVsBGcorr[10][200];
    for(Int_t i=0;i<10;i++) {
      snprintf(name_fPmeanVsBGcorr[i],200,"fPmeanVsBGcorr_%s",nameTemp[i]);
      fPmeanVsBGcorr[0][i]=new TF1(name_fPmeanVsBGcorr[i],"[2]-[0]*TMath::Power(x,[1]);p_{vtx}/m;<p>/p",0.0001,100);
      fPmeanVsBGcorr[1][i]=new TF1(name_fPmeanVsBGcorr[i],"[2]-[0]*TMath::Power(x,[1]);p_{vtx}/m;<p>/p",0.0001,100);
      //fPmeanVsBGcorr[i]->SetParameters(pars_fPmeanVsBGcorr[i]);
      //fPmeanVsBGcorr[i]->SetNpx(fPmeanVsBGcorr[i]->GetNpx()*10);
    }
    SetPmeanCorrections();
       
    Char_t name_prPmeanVsBGcorr[10][200];
    Char_t title_prPmeanVsBGcorr[10][200];
   
    hbins[0]=200;
    for(Int_t iS=0;iS<10;iS++) {
      snprintf(name_prPmeanVsBGcorr[iS],200,"prPmeanVsBGcorr_%s",nameTemp[iS]);
      snprintf(title_prPmeanVsBGcorr[iS],200,"<p>/p_{vtx} vs #beta#gamma of %s as parameterized in input TF1 (in DCAxyCut);p_{vtx}/m_{%s};<p>_{%s}/p_{vtx}",nameTemp[iS],nameTemp[iS],nameTemp[iS]);
      prPmeanVsBGcorr[iB][iS]=new TProfile(name_prPmeanVsBGcorr[iS],title_prPmeanVsBGcorr[iS],hbins[0],0,10,0.8,1.2,"");
    }	

    fList[iB]->Add(htemp[iB]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(hCentrality[iB][i]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(hZvertex[iB][i]);
    fList[iB]->Add(hEta[iB]);
    fList[iB]->Add(hPhi[iB]);
    //fList[iB]->Add(fEtaPhi[iB]);
    fList[iB]->Add(hNTpcCluster[iB]);
    fList[iB]->Add(hNTrdSlices[iB]);
    //for(Int_t i=0;i<2;i++) fList[iB]->Add(fdEdxVSp[iB][i]);
    //for(Int_t i=0;i<nPart;i++) fList[iB]->Add(hDeDxExp[iB][i]);
    //for(Int_t i=0;i<nPart;i++) fList[iB]->Add(fNsigmaTpc[iB][i]);
    for(Int_t i=0;i<nPart;i++) {
      if(kSignalCheck!=1) 
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      //fList[iB]->Add(fNsigmaTpc_kTOF[iB][i]);
      //fList[iB]->Add(fNsigmaTpc_kTOF[iB][i+nPart]);
    }
    //for(Int_t i=0;i<2;i++) fList[iB]->Add(fBetaTofVSp[iB][i]);
    //for(Int_t i=0;i<nPart;i++) fList[iB]->Add(hBetaExp[iB][i]);
    //for(Int_t i=0;i<nPart;i++) fList[iB]->Add(fNsigmaTof[iB][i]);
    for(Int_t i=0;i<nPart;i++) {
      if(kSignalCheck!=1) 
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      //fList[iB]->Add(fNsigmaTof_DcaCut[iB][i]);
      //fList[iB]->Add(fNsigmaTof_DcaCut[iB][i+nPart]);
    }
    //for(Int_t i=0;i<2;i++) fList[iB]->Add(fM2vsP_NoTpcCut[iB][0][i]);
    //for(Int_t i=0;i<2;i++) fList[iB]->Add(fM2vsP_NoTpcCut[iB][1][i]);
    for(Int_t i=0;i<nPart;i++) {
      if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      //fList[iB]->Add(fM2vsP[iB][0][i]);
      //fList[iB]->Add(fM2vsP[iB][0][i+nPart]);
    }
    for(Int_t i=0;i<nPart;i++){
      if(kSignalCheck!=1) 
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      fList[iB]->Add(fM2vsP[iB][1][i]);
      fList[iB]->Add(fM2vsP[iB][1][i+nPart]);
    }
    for(Int_t i=0;i<2;i++){
      //fList[iB]->Add(fPvtxTrueVsReco[i]);
      fList[iB]->Add(prPvtxTrueVsReco[iB][i]);
    }
    if(iMtof!=1) {
      for(Int_t i=0;i<nPart;i++){
	if(i<2) continue;//e,mu excluded
	fList[iB]->Add(fPmeanVsBetaGamma[iB][i]);
	fList[iB]->Add(prPmeanVsBetaGamma[iB][i]);
	fList[iB]->Add(fPmeanVsBetaGamma[iB][i+nPart]);
	fList[iB]->Add(prPmeanVsBetaGamma[iB][i+nPart]);
      }
    }
    if(iMtof>2) {
      //for(Int_t i=0;i<10;i++)fList[iB]->Add(fPmeanVsBGcorr[i]);
      for(Int_t i=0;i<10;i++)fList[iB]->Add(prPmeanVsBGcorr[iB][i]);
    }
    //for(Int_t i=0;i<10;i++) fList[iB]->Add(fM2vsZ[iB][i]);
    for(Int_t i=0;i<nPart;i++){
      if(kSignalCheck!=1) 
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      for(Int_t j=0;j<nbin;j++){
	fList[iB]->Add(hDCAxy[iB][i][j]);
	fList[iB]->Add(hDCAz[iB][i][j]);
	fList[iB]->Add(hM2CutDCAxy[iB][i][j]);
	fList[iB]->Add(hM2CutGroundDCAxy[iB][i][j]);
	fList[iB]->Add(hDCAxy[iB][i+nPart][j]);
	fList[iB]->Add(hDCAz[iB][i+nPart][j]);
	fList[iB]->Add(hM2CutDCAxy[iB][i+nPart][j]);
	fList[iB]->Add(hM2CutGroundDCAxy[iB][i+nPart][j]);
      }
    }
    
    // Post output data.
    PostData(1, fList[0]);
    PostData(2, fList[1]);
        
  }//end iB loop
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
  fPIDResponse=inputHandler->GetPIDResponse();
  
  //--------------------------Magnetic field polarity--------------------
  Double_t fBfield=fEvent->GetMagneticField();
  if(fBfield<0.0) iBconf=0;//B--
  else iBconf=1;//B++
  for(Int_t i=0;i<nBconf;i++) htemp[i]->Fill(fBfield);
    
  //--------------------------Centrality--------------------------------
  Double_t v0Centr  = -10.;
  AliCentrality *centrality = fEvent->GetCentrality();
  if (centrality){
    v0Centr=centrality->GetCentralityPercentile("V0M"); // VZERO
  }
  hCentrality[iBconf][0]->Fill(v0Centr);

  //-------------------------zVertex determination of event----------------
  Double_t zvtx = 9999.9;
  const AliVVertex* vtxEVENT = fEvent->GetPrimaryVertex();
  if(vtxEVENT->GetNContributors()>0) zvtx = vtxEVENT->GetZ();
  
  hZvertex[iBconf][0]->Fill(zvtx);
  
  //---------------------------EVENT CUTS-----------------------------
  if(TMath::Abs(zvtx) < 10.0 && v0Centr>Centrality[0] && v0Centr<Centrality[1]){

    hCentrality[iBconf][1]->Fill(v0Centr);
    hZvertex[iBconf][1]->Fill(zvtx);
    
    Int_t nTracks = fEvent->GetNumberOfTracks();
    
    //----------------------loop on the TRACKS-----------------------------
    for(Int_t iT = 0; iT < nTracks; iT++) { 
      AliVTrack* track = (AliVTrack *) fEvent->GetTrack(iT);
      
      if (!track){
	continue;
      }
      
     //For the geometrical cuts
      Double_t eta = track->Eta();
      
      Bool_t trkFlag = 0;
      trkFlag = ((AliAODTrack *) track)->TestFilterBit(FilterBit);
      //TestFilterBit(16) -- Standard Cuts with very loose DCA: GetStandardITSTPCTrackCuts2011(kFALSE) && SetMaxDCAToVertexXY(2.4) && SetMaxDCAToVertexZ(3.2) && SetDCaToVertex2D(kTRUE)
      //TestFilterBit(32) (STARDARD) -- Standard Cuts with very tight DCA cut ( 7sigma^primaries: 7*(0.0015+0.0050/pt^1.1) ) : GetStandardITSTPCTrackCuts2011(). 
      
      //Cut on the Minumum Number of the TPC clusters
      Bool_t isMinTpcCluster=kFALSE;
      Int_t nTpcCluster=0;
      nTpcCluster=track->GetTPCNcls();
      if(nTpcCluster>NminTpcCluster) isMinTpcCluster=kTRUE;

      //-------------------------------------start TRACK CUTS----------------------------------
      if ((track->Pt() < 0.2) || (eta<EtaLimit[0]) || (eta>EtaLimit[1]) || !trkFlag || !isMinTpcCluster)
	continue;
      
      //Vertex determination
      Double_t b[2] = {-99., -99.};
      Double_t bCov[3] = {-99., -99., -99.};
      if (!track->PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 100., b, bCov))
	continue;
      
      Double_t DCAxy = b[0];
      Double_t DCAz = b[1];
      
      //Cut on the DCAz
      Bool_t isDCAzCut=kFALSE;
      if(DCAz<DCAzCut) isDCAzCut=kTRUE;

      if(!isDCAzCut)
	continue;
      
      //For the Tpc purity cut
      Double_t dedx = track->GetTPCsignal();
      if(dedx<10) continue;

      Int_t nTrdSlices = track->GetNumberOfTRDslices();
      if(nTrdSlices<2 && iTrdCut==1) continue; 
      if(nTrdSlices>0 && iTrdCut==2) continue;
      
      //-------------------------------------end TRACK CUTS----------------------------------

      //-------------------------------------Track info--------------------------------------
      Double_t phi= track->Phi();
      
      hEta[iBconf]->Fill(eta);
      hPhi[iBconf]->Fill(phi);
      fEtaPhi[iBconf]->Fill(eta,phi);
      hNTpcCluster[iBconf]->Fill(nTpcCluster);
      hNTrdSlices[iBconf]->Fill(nTrdSlices);
	
      Double_t charge = (Double_t)track->Charge();
      Double_t p = track->P();
      Double_t pt = track->Pt();
      Double_t tof  = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
      Double_t pTPC = track->GetTPCmomentum();
      Double_t beta = 0.0;
      Double_t M2 = 999.9;
      Double_t Z2 = 999.9;
           
      kTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
      
      //-----------------------------TPC info------------------------------
      Double_t nsigmaTPC[nPart];
      Double_t expdedx[nPart];
      
      Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,#mu,#pi,K,p,d,t,3He,4He
      Int_t FlagPid = 0;
      
      for(Int_t iS=0;iS<9;iS++){
	nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
	//TPC identification:
	if(TMath::Abs(nsigmaTPC[iS])<NsigmaTpcCut) {
	  FlagPid += ((Int_t)TMath::Power(2,iS));
	}
      }
      //Correction of the momentum to the vertex for (anti)nuclei
      Double_t pC[9];
      for(Int_t iS=0;iS<9;iS++)	pC[iS]=p;
      this->MomVertexCorrection(p,pC,eta,FlagPid);

      //More TPC info:
      for(Int_t iS=0;iS<9;iS++){
	expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
	hDeDxExp[iBconf][iS]->Fill(pTPC,expdedx[iS]);
	nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
	fNsigmaTpc[iBconf][iS]->Fill(pTPC,nsigmaTPC[iS]);
	if(charge>0) {//positive particle
	  if(kTOF && (TMath::Abs(DCAxy)<DCAxyCut)) fNsigmaTpc_kTOF[iBconf][iS]->Fill(p,nsigmaTPC[iS]);
	}
	else {//negative particle
	  if(kTOF && (TMath::Abs(DCAxy)<DCAxyCut)) fNsigmaTpc_kTOF[iBconf][iS+nPart]->Fill(p,nsigmaTPC[iS]);
	}
	/*
	  if(TMath::Abs(nsigmaTPC[iS])<NsigmaTpcCut) {
	  FlagPid += ((Int_t)TMath::Power(2,iS));
	  }*/
      }
          
      if(charge>0) fdEdxVSp[iBconf][0]->Fill(pTPC,dedx);
      else fdEdxVSp[iBconf][1]->Fill(pTPC,dedx);

      //-----------------------------TOF info------------------------------
      
      Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};

      //----------------------------------------kTOF available-----------------------------
           
      if(kTOF) {
	Double_t exptimes[9];
	track->GetIntegratedTimes(exptimes);
	//Integrated times of the Nuclei:
	for(Int_t iN=5;iN<9;iN++) {
	  exptimes[iN] = exptimes[4]*exptimes[4]*(massOverZ[iN]*massOverZ[iN]/p/p+1)/(massOverZ[4]*massOverZ[4]/p/p+1);
	  exptimes[iN] = TMath::Sqrt(exptimes[iN]);
	}  
	
	beta=exptimes[0];
	beta=beta/tof;//beta = L/tof/c = t_e/tof
	
	Int_t FlagPidTof = 0;
	Double_t NsigmaTofCut = 2.0;
	
	Double_t nsigmaTOF[9];
	for(Int_t iS=0;iS<9;iS++){
	  nsigmaTOF[iS] = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType) iS);
	  fNsigmaTof[iBconf][iS]->Fill(pt,nsigmaTOF[iS]);
	  if(charge>0) {
	    hBetaExp[iBconf][iS]->Fill(p,exptimes[0]/exptimes[iS]);
	    if(TMath::Abs(DCAxy)<DCAxyCut) fNsigmaTof_DcaCut[iBconf][iS]->Fill(pt,nsigmaTOF[iS]);
	  }
	  else {
	    hBetaExp[iBconf][iS+nPart]->Fill(p,exptimes[0]/exptimes[iS]);
	    if(TMath::Abs(DCAxy)<DCAxyCut) fNsigmaTof_DcaCut[iBconf][iS+nPart]->Fill(pt,nsigmaTOF[iS]);
	  }

	  //TOF identification:
	  if(TMath::Abs(nsigmaTOF[iS])<NsigmaTofCut) {
	    FlagPidTof += ((Int_t)TMath::Power(2,iS));
	  }
	}
	
	if(charge>0) fBetaTofVSp[iBconf][0]->Fill(p,beta);
	else fBetaTofVSp[iBconf][1]->Fill(p,beta);
		
	this->GetMassFromPvertex(beta,p,M2);
	this->GetZTpc(dedx,pTPC,M2,Z2);
	
	Double_t Mass2[9];
	//-----------------------------M2 as a function of momentum to the primary vertex if iMtof==1---------------------------------
	if(iMtof==1) this->GetMassFromPvertexCorrected(beta,pC,Mass2);

	if(iMtof>1) this->GetPmeanVsBetaGamma(exptimes,pC,FlagPid,FlagPidTof,charge,DCAxy);
	
	//-----------------------------M2 as a function of expected times---------------------------------
	if(iMtof==2) this->GetMassFromExpTimes(beta,exptimes,Mass2);
	 
	//-----------------------------M2 as a function of mean momentum calculated from expected time and extrapolated to the (anti)nuclei---------------------------------
	if(iMtof>2) this->GetMassFromMeanMom(beta,exptimes,pC,charge,Mass2,FlagPid,FlagPidTof,DCAxy);

	//-------------------------------Squared Mass TH2 distributions-----------------------
	if(charge>0) {
	  //without TPC
	  fM2vsP_NoTpcCut[iBconf][0][0]->Fill(M2,p);
	  if(TMath::Abs(DCAxy)<DCAxyCut) fM2vsP_NoTpcCut[iBconf][1][0]->Fill(M2,p);
	  //with TPC
	  for(Int_t iS=0;iS<9;iS++) {
	    M2=999.9;
	    M2=Mass2[iS];
	    //-----------------
	    if(FlagPid & stdFlagPid[iS]) {
	      fM2vsP[iBconf][0][iS]->Fill(M2,pC[iS]);
	      if(TMath::Abs(DCAxy)<DCAxyCut) fM2vsP[iBconf][1][iS]->Fill(M2,pC[iS]);
	    }
	  }
	}
	else {//charge<0
	  //without TPC
	  fM2vsP_NoTpcCut[iBconf][0][1]->Fill(M2,p);
	  if(TMath::Abs(DCAxy)<DCAxyCut) fM2vsP_NoTpcCut[iBconf][1][1]->Fill(M2,p);
	   //with TPC
	  for(Int_t iS=0;iS<9;iS++) {
	    M2=999.9;
	    M2=Mass2[iS];
	    //-----------------
	    if(FlagPid & stdFlagPid[iS]) {
	      fM2vsP[iBconf][0][iS+nPart]->Fill(M2,pC[iS]);
	      if(TMath::Abs(DCAxy)<DCAxyCut) fM2vsP[iBconf][1][iS+nPart]->Fill(M2,pC[iS]);
	    }
	  }
	}
	
	//------------------------------start DCA and Squared Mass TH1 distributions-------------------------
	Double_t binP[nbin+1];
	for(Int_t i=0;i<nbin+1;i++) {
	  binP[i]=0.4+i*0.1;
	}

	if(charge>0) {
	  for(Int_t iS=0;iS<9;iS++) {
	    M2=999.9;
	    M2=Mass2[iS];
	    
	    if(FlagPid & stdFlagPid[iS]) {
	      for(Int_t j=0;j<nbin;j++) {
		if(pC[iS]>binP[j] && pC[iS]<binP[j+1]) {
		  hDCAxy[iBconf][iS][j]->Fill(DCAxy);
		  hDCAxy[iBconf][iS][j]->Fill(-DCAxy);
		  hDCAz[iBconf][iS][j]->Fill(DCAz);
		  hDCAz[iBconf][iS][j]->Fill(-DCAz);
		  if(TMath::Abs(DCAxy)<DCAxyCut) {
		    hM2CutDCAxy[iBconf][iS][j]->Fill(M2);
		  }
		  if(TMath::Abs(DCAxy+0.5)<DCAxyCut) {
		    hM2CutGroundDCAxy[iBconf][iS][j]->Fill(M2);
		  }
		  break;
		}
	      }//end loop on the p bins (j)
	    }
	  }//end loop on the particle species (iS)
	}
	else {//charge<0
	  for(Int_t iS=0;iS<9;iS++) {
	    M2=999.9;
	    M2=Mass2[iS];
	    	    
	    if(FlagPid & stdFlagPid[iS]) {
	      for(Int_t j=0;j<nbin;j++) {
		if(pC[iS]>binP[j] && pC[iS]<binP[j+1]) {
		  hDCAxy[iBconf][iS+nPart][j]->Fill(DCAxy);
		  hDCAxy[iBconf][iS+nPart][j]->Fill(-DCAxy);
		  hDCAz[iBconf][iS+nPart][j]->Fill(DCAz);
		  hDCAz[iBconf][iS+nPart][j]->Fill(-DCAz);
		  if(TMath::Abs(DCAxy)<DCAxyCut) {
		    hM2CutDCAxy[iBconf][iS+nPart][j]->Fill(M2);
		  }
		  if(TMath::Abs(DCAxy+0.5)<DCAxyCut) {
		    hM2CutGroundDCAxy[iBconf][iS+nPart][j]->Fill(M2);
		  }
		  break;
		}
	      }//end loop on the p bins (j)
	    }
	  }//end loop on the particle species (iS)
	}
		
	//-------------------------------------------------M2/Z2 vs Z-------------------------
	

	Double_t binCutPt[10] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
	Double_t Z=999.9;
	if(Z2>0) Z=TMath::Sqrt(Z2);
	
	fM2vsZ[iBconf][0]->Fill(charge*TMath::Sqrt(Z2),M2);
	for(Int_t i=1;i<10;i++) {
	  if(pt>binCutPt[i-1] && pt<binCutPt[i]){
	    fM2vsZ[iBconf][i]->Fill(charge*Z,M2);
	    break;
	  }
	}
	
      }//end kTOF available
    }//end track loop
  }//end loop on the events
}

//_____________________________________________________________________________
void AliAnalysisNucleiMass::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
//_____________________________________________________________________________
void AliAnalysisNucleiMass::MomVertexCorrection(Double_t p, Double_t *pC, Double_t eta, Int_t FlagPid){

  Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,#mu,#pi,K,p,d,t,3He,4He

  for(Int_t iS=0;iS<9;iS++) {
    if(FlagPid & stdFlagPid[iS]) {
      if(iS==5) {
	if(kPvtxCorr==1) pC[iS]=pC[iS]*fPvtxTrueVsReco[0]->Eval(pC[iS],TMath::Abs(eta));//for (bar)d
	prPvtxTrueVsReco[iBconf][0]->Fill(p,pC[iS]/p);
      }
      else if(iS==7) {
	if(kPvtxCorr==1) pC[iS]=pC[iS]*fPvtxTrueVsReco[1]->Eval(pC[iS],TMath::Abs(eta));//for (bar)He3
	prPvtxTrueVsReco[iBconf][1]->Fill(p,pC[iS]/p);
      }
    }
  }
  
  return;
  
}
//_____________________________________________________________________________
void AliAnalysisNucleiMass::GetMassFromPvertex(Double_t beta, Double_t p, Double_t &M2) {
  
  M2 = p*p*(1-beta*beta)/(beta*beta);

  return;
  
}
//_________________________________________________________________________________________________________________________
void AliAnalysisNucleiMass::GetZTpc(Double_t dedx, Double_t pTPC, Double_t M2, Double_t &Z2) {

  //z^2_tpc = dedx^{Tpc} / dedx^{exp,Tof}_{z=1}
  
  Z2=999.9;
  
  Double_t M=999.9;
  Double_t pTPC_pr=999.9;//rescaling of the pTPC for the proton
  Double_t expdedx_Tof=999.9;
  
  if(M2>0) {
    M=TMath::Sqrt(M2);
    pTPC_pr=pTPC*0.938272/M;
    expdedx_Tof=fPIDResponse->GetTPCResponse().GetExpectedSignal(pTPC_pr,AliPID::kProton);
    if((dedx/expdedx_Tof)<0) return;
    Z2=TMath::Power(dedx/expdedx_Tof,0.862);
  }
  
  return;
}
//_________________________________________________________________________________________________________________________
void AliAnalysisNucleiMass::GetMassFromPvertexCorrected(Double_t beta, Double_t *pC, Double_t *Mass2) {
  
  for(Int_t iS=0;iS<9;iS++) Mass2[iS] = pC[iS]*pC[iS]*(1-beta*beta)/(beta*beta);

  return;
}  
//____________________________________________________________________________________________________________
void AliAnalysisNucleiMass::GetMassFromExpTimes(Double_t beta, Double_t *IntTimes, Double_t *Mass2) {
 
  // m = p_exp/beta/gamma where p_exp = mPDG*beta_exp*gamma_exp; beta_exp = L/t_exp/c = t_e/t_exp ; beta=L/tof/c = t_e/tof
  // In this way m_tof = mPDG only if tof=t_exp
  
  Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};
    
  Double_t beta2Exp[9];
  Double_t p2Exp[9];
  
  //Double_t pExp[9];
  
  for(Int_t iS=0;iS<9;iS++) {
    beta2Exp[iS]=IntTimes[0]/IntTimes[iS];//beta = L/tof*c = t_e/tof
    beta2Exp[iS]=beta2Exp[iS]*beta2Exp[iS];
    if((1-beta2Exp[iS])==0) {
      Mass2[iS]=999.9;
      continue;
    }
    p2Exp[iS]=massOverZ[iS]*massOverZ[iS]*beta2Exp[iS]/(1-beta2Exp[iS]);
    
    //--------------------for MC corrections
    if(p2Exp[iS]<0) {
      Mass2[iS]=999.9;
      continue;
    }
    //pExp[iS]=TMath::Sqrt(p2Exp[iS]);
    
    //------------
    Mass2[iS]=p2Exp[iS]*(1-beta*beta)/(beta*beta);
  }//end loop on the particle species
  
  return;
}
//____________________________________________________________________________________________________________
void AliAnalysisNucleiMass::GetPmeanVsBetaGamma(Double_t *IntTimes, Double_t *pVtx, Int_t FlagPid, Int_t FlagPidTof, Double_t charge, Double_t DCAxy) {
 
  // m = p_exp/beta/gamma where p_exp = mPDG*beta_exp*gamma_exp; beta_exp = L/t_exp/c = t_e/t_exp ; beta=L/tof/c = t_e/tof
  // In this way m_tof = mPDG only if tof=t_exp
  
  Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};
    
  Double_t beta2Exp[9];
  Double_t p2Exp[9];
  
  Double_t pExp[9];
  
  Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,#mu,#pi,K,p,d,t,3He,4He

  for(Int_t iS=0;iS<9;iS++) {
    beta2Exp[iS]=IntTimes[0]/IntTimes[iS];//beta = L/tof*c = t_e/tof
    beta2Exp[iS]=beta2Exp[iS]*beta2Exp[iS];
    if((1-beta2Exp[iS])==0) {
      continue;
    }
    p2Exp[iS]=massOverZ[iS]*massOverZ[iS]*beta2Exp[iS]/(1-beta2Exp[iS]);
    
    if(p2Exp[iS]<0) {
      continue;
    }
    pExp[iS]=TMath::Sqrt(p2Exp[iS]);
       
    if((FlagPid & stdFlagPid[iS]) && (FlagPidTof & stdFlagPid[iS])) {
      if(TMath::Abs(DCAxy)>DCAxyCut) continue;
      if(charge>0){
	fPmeanVsBetaGamma[iBconf][iS]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
	prPmeanVsBetaGamma[iBconf][iS]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
      }
      else {
	fPmeanVsBetaGamma[iBconf][iS+nPart]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
	prPmeanVsBetaGamma[iBconf][iS+nPart]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
      }
    }
  }//end loop on the particle species
  
  return;
  
}
//____________________________________________________________________________________________________________
void AliAnalysisNucleiMass::GetMassFromMeanMom(Double_t beta, Double_t *IntTimes, Double_t *pVtx, Double_t charge, Double_t *Mass2, Int_t FlagPid, Int_t FlagPidTof, Double_t DCAxy) {//Double_t *Mass2, Int_t iCorr
 
  // m = p_exp/beta/gamma where p_exp = mPDG*beta_exp*gamma_exp; beta_exp = L/t_exp/c = t_e/t_exp ; beta=L/tof/c = t_e/tof
  // In this way m_tof = mPDG only if tof=t_exp
  
  Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};
    
  Double_t beta2Exp[9];
  Double_t p2Exp[9];
  
  Double_t pExp[9];
  
  Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,#mu,#pi,K,p,d,t,3He,4He

  for(Int_t iS=0;iS<9;iS++) {
    if(iS==2 || iS==3 || iS==4 || iS==5 || iS==7) {
      if(charge>0) {
	if(iS!=7) p2Exp[iS]=pVtx[iS]*fPmeanVsBGcorr[iBconf][iS-2]->Eval(pVtx[iS]/massOverZ[iS]);
	else p2Exp[iS]=pVtx[iS]*fPmeanVsBGcorr[iBconf][iS-3]->Eval(pVtx[iS]/massOverZ[iS]);
      }
      else if(charge<0) {
	if(iS!=7) p2Exp[iS]=pVtx[iS]*fPmeanVsBGcorr[iBconf][iS+3]->Eval(pVtx[iS]/massOverZ[iS]);
	else p2Exp[iS]=pVtx[iS]*fPmeanVsBGcorr[iBconf][iS+2]->Eval(pVtx[iS]/massOverZ[iS]);
      }
      p2Exp[iS]*=p2Exp[iS];
    }
    else {
      beta2Exp[iS]=IntTimes[0]/IntTimes[iS];//beta = L/tof*c = t_e/tof
      beta2Exp[iS]=beta2Exp[iS]*beta2Exp[iS];
      if((1-beta2Exp[iS])==0) {
	Mass2[iS]=999.9;
	continue;
      }
      p2Exp[iS]=massOverZ[iS]*massOverZ[iS]*beta2Exp[iS]/(1-beta2Exp[iS]);
    }
    //--------------------for MC corrections
    if(p2Exp[iS]<0) {
      Mass2[iS]=999.9;
      continue;
    }
    pExp[iS]=TMath::Sqrt(p2Exp[iS]);
    
    //------------
    Mass2[iS]=p2Exp[iS]*(1-beta*beta)/(beta*beta);
    
    //-----------
    if(TMath::Abs(DCAxy)>DCAxyCut) continue;
    if(iS==2 || iS==3 || iS==4 || iS==5 || iS==7) {
      if((FlagPid & stdFlagPid[iS]) && (FlagPidTof & stdFlagPid[iS])) {
	if(charge>0) {
	  if(iS!=7) prPmeanVsBGcorr[iBconf][iS-2]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
	  else prPmeanVsBGcorr[iBconf][iS-3]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
	}
	else if(charge<0) {
	  if(iS!=7) prPmeanVsBGcorr[iBconf][iS+3]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
	  else prPmeanVsBGcorr[iBconf][iS+2]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
	}
      }
    }

  }//end loop on the particle species
  
  return;
  
}
//________________________________________________________________________________________
void AliAnalysisNucleiMass::SetPmeanCorrections(){
  
  //iMtof==8 -> different particle and antiparticle parameterization 
  
  Double_t pars_fPmeanVsBGcorr[nBconf][10][3];
  //particle

  Double_t etaMin=0.0;
  Double_t etaMax=0.8;

  if(EtaLimit[0]<0.0 || EtaLimit[1]<0.0) {
    etaMin=TMath::Abs(EtaLimit[1]);
    etaMax=TMath::Abs(EtaLimit[0]);
  }
  else {
    etaMin=TMath::Abs(EtaLimit[0]);
    etaMax=TMath::Abs(EtaLimit[1]);
  }

  if(etaMin>-0.00001 && etaMax<0.10001) {
    //printf("EtaLimit[0]== %f and EtaLimit[1]== %fAAA\n",EtaLimit[0],EtaLimit[1]);
    for(Int_t i=0;i<5;i++) {
      if(i==0) {//pi
	pars_fPmeanVsBGcorr[0][i][0]=4.16853e-02; pars_fPmeanVsBGcorr[0][i][1]=-7.67091e-01; pars_fPmeanVsBGcorr[0][i][2]=9.98035e-01;//B--
	pars_fPmeanVsBGcorr[1][i][0]=5.51380e-02; pars_fPmeanVsBGcorr[1][i][1]=-7.58112e-01; pars_fPmeanVsBGcorr[1][i][2]=1.00360e+00;//B++
      }
      else if(i==1) {//K
	pars_fPmeanVsBGcorr[0][i][0]=2.73697e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.43042e+00; pars_fPmeanVsBGcorr[0][i][2]=9.93148e-01;
	pars_fPmeanVsBGcorr[1][i][0]=3.19397e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.08037e+00; pars_fPmeanVsBGcorr[1][i][2]=9.98016e-01;
      }
      else if(i==2) {//p
	pars_fPmeanVsBGcorr[0][i][0]=1.35721e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.80958e+00; pars_fPmeanVsBGcorr[0][i][2]=9.93925e-01;
	pars_fPmeanVsBGcorr[1][i][0]=1.63564e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.55914e+00; pars_fPmeanVsBGcorr[1][i][2]=9.98106e-01;
      }
      else if(i==3) {//d
	pars_fPmeanVsBGcorr[0][i][0]=0.009609; pars_fPmeanVsBGcorr[0][i][1]=-2.534810; pars_fPmeanVsBGcorr[0][i][2]=0.993507;
	pars_fPmeanVsBGcorr[1][i][0]=0.011580; pars_fPmeanVsBGcorr[1][i][1]=-2.308857; pars_fPmeanVsBGcorr[1][i][2]=0.998126;
      }
      else if(i==4) {//He3
	pars_fPmeanVsBGcorr[0][i][0]=0.026420; pars_fPmeanVsBGcorr[0][i][1]=-2.253066; pars_fPmeanVsBGcorr[0][i][2]=0.993507;
	pars_fPmeanVsBGcorr[1][i][0]=0.031840; pars_fPmeanVsBGcorr[1][i][1]=-2.052228; pars_fPmeanVsBGcorr[1][i][2]=0.998126;
      }
    }
  }
  else if(etaMin>0.09999 && etaMax<0.20001) {
    //printf("EtaLimit[0]== %f and EtaLimit[1]== %fBBB\n",EtaLimit[0],EtaLimit[1]);
    for(Int_t i=0;i<5;i++) {
      if(i==0) {//pi
	pars_fPmeanVsBGcorr[0][i][0]=4.98872e-02; pars_fPmeanVsBGcorr[0][i][1]=-3.56884e-01; pars_fPmeanVsBGcorr[0][i][2]=1.01356e+00;
	pars_fPmeanVsBGcorr[1][i][0]=6.11287e-02; pars_fPmeanVsBGcorr[1][i][1]=-3.65072e-01; pars_fPmeanVsBGcorr[1][i][2]=1.02074e+00;
      }
      else if(i==1) {//K
	pars_fPmeanVsBGcorr[0][i][0]=2.85027e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.04376e+00; pars_fPmeanVsBGcorr[0][i][2]=9.94804e-01;
	pars_fPmeanVsBGcorr[1][i][0]=3.30937e-02; pars_fPmeanVsBGcorr[1][i][1]=-1.72959e+00; pars_fPmeanVsBGcorr[1][i][2]=9.99966e-01;
      }
      else if(i==2) {//p
	pars_fPmeanVsBGcorr[0][i][0]=1.38640e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.71621e+00; pars_fPmeanVsBGcorr[0][i][2]=9.94151e-01;
	pars_fPmeanVsBGcorr[1][i][0]=1.74869e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.38269e+00; pars_fPmeanVsBGcorr[1][i][2]=9.98776e-01;
      }
      else if(i==3) {//d
	pars_fPmeanVsBGcorr[0][i][0]=0.009816; pars_fPmeanVsBGcorr[0][i][1]=-2.450567; pars_fPmeanVsBGcorr[0][i][2]=0.994465;
	pars_fPmeanVsBGcorr[1][i][0]=0.012381; pars_fPmeanVsBGcorr[1][i][1]=-2.149671; pars_fPmeanVsBGcorr[1][i][2]=0.999302;
      }
      else if(i==4) {//He3
	pars_fPmeanVsBGcorr[0][i][0]=0.026988; pars_fPmeanVsBGcorr[0][i][1]=-2.178186; pars_fPmeanVsBGcorr[0][i][2]=0.994465;
	pars_fPmeanVsBGcorr[1][i][0]=0.034041; pars_fPmeanVsBGcorr[1][i][1]=-1.910736; pars_fPmeanVsBGcorr[1][i][2]=0.999302;
      }
    }
  }
  else if(etaMin>0.19999 && etaMax<0.30001) {
    //printf("EtaLimit[0]== %f and EtaLimit[1]== %fCCC\n",EtaLimit[0],EtaLimit[1]);
    for(Int_t i=0;i<5;i++) {
      if(i==0) {//pi
	pars_fPmeanVsBGcorr[0][i][0]=4.71844e-02; pars_fPmeanVsBGcorr[0][i][1]=-6.24048e-01; pars_fPmeanVsBGcorr[0][i][2]=1.00525e+00;
	pars_fPmeanVsBGcorr[1][i][0]=5.45281e-02; pars_fPmeanVsBGcorr[1][i][1]=-5.87331e-01; pars_fPmeanVsBGcorr[1][i][2]=1.01029e+00;
      }
      else if(i==1) {//K
	pars_fPmeanVsBGcorr[0][i][0]=2.92060e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.15537e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97130e-01;
	pars_fPmeanVsBGcorr[1][i][0]=3.24550e-02; pars_fPmeanVsBGcorr[1][i][1]=-1.97289e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00059e+00;
      }
      else if(i==2) {//p
	pars_fPmeanVsBGcorr[0][i][0]=1.33594e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.86707e+00; pars_fPmeanVsBGcorr[0][i][2]=9.96053e-01;
	pars_fPmeanVsBGcorr[1][i][0]=1.57187e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.62957e+00; pars_fPmeanVsBGcorr[1][i][2]=9.99431e-01;
      }
      else if(i==3) {//d
	pars_fPmeanVsBGcorr[0][i][0]=0.009458; pars_fPmeanVsBGcorr[0][i][1]=-2.586677; pars_fPmeanVsBGcorr[0][i][2]=0.996592;
	pars_fPmeanVsBGcorr[1][i][0]=0.011129; pars_fPmeanVsBGcorr[1][i][1]=-2.372404; pars_fPmeanVsBGcorr[1][i][2]=1.000024;
      }
      else if(i==4) {//He3
	pars_fPmeanVsBGcorr[0][i][0]=0.026006; pars_fPmeanVsBGcorr[0][i][1]=-2.299168; pars_fPmeanVsBGcorr[0][i][2]=0.996592;
	pars_fPmeanVsBGcorr[1][i][0]=0.030599; pars_fPmeanVsBGcorr[1][i][1]=-2.108711; pars_fPmeanVsBGcorr[1][i][2]=1.000024;
      }
    }
  }
  else if(etaMin>0.29999 && etaMax<0.40001) {
    //printf("EtaLimit[0]== %f and EtaLimit[1]== %fDDD\n",EtaLimit[0],EtaLimit[1]);
    for(Int_t i=0;i<5;i++) {
      if(i==0) {//pi
	pars_fPmeanVsBGcorr[0][i][0]=5.25262e-02; pars_fPmeanVsBGcorr[0][i][1]=-3.04325e-01; pars_fPmeanVsBGcorr[0][i][2]=1.02056e+00;
	pars_fPmeanVsBGcorr[1][i][0]=5.70585e-02; pars_fPmeanVsBGcorr[1][i][1]=-5.95375e-01; pars_fPmeanVsBGcorr[1][i][2]=1.01130e+00;
      }
      else if(i==1) {//K
	pars_fPmeanVsBGcorr[0][i][0]=2.96035e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.17931e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97539e-01;
	pars_fPmeanVsBGcorr[1][i][0]=3.35067e-02; pars_fPmeanVsBGcorr[1][i][1]=-1.99656e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00128e+00;
      }
      else if(i==2) {//p
	pars_fPmeanVsBGcorr[0][i][0]=1.44529e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.77844e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97130e-01;
	pars_fPmeanVsBGcorr[1][i][0]=1.68180e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.56489e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00070e+00;
      }
      else if(i==3) {//d
	pars_fPmeanVsBGcorr[0][i][0]=0.010233; pars_fPmeanVsBGcorr[0][i][1]=-2.506714; pars_fPmeanVsBGcorr[0][i][2]=0.997341;
	pars_fPmeanVsBGcorr[1][i][0]=0.011907; pars_fPmeanVsBGcorr[1][i][1]=-2.314052; pars_fPmeanVsBGcorr[1][i][2]=1.001048;
      }
      else if(i==4) {//He3
	pars_fPmeanVsBGcorr[0][i][0]=0.028135; pars_fPmeanVsBGcorr[0][i][1]=-2.228093; pars_fPmeanVsBGcorr[0][i][2]=0.997341;
	pars_fPmeanVsBGcorr[1][i][0]=0.032739; pars_fPmeanVsBGcorr[1][i][1]=-2.056845; pars_fPmeanVsBGcorr[1][i][2]=1.001048;
      }
    }
  }
  else if(etaMin>0.39999 && etaMax<0.50001) {
    //printf("EtaLimit[0]== %f and EtaLimit[1]== %fEEE\n",EtaLimit[0],EtaLimit[1]);
    for(Int_t i=0;i<5;i++) {
      if(i==0) {//pi
	pars_fPmeanVsBGcorr[0][i][0]=5.72833e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.51868e-01; pars_fPmeanVsBGcorr[0][i][2]=1.02665e+00;
	pars_fPmeanVsBGcorr[1][i][0]=6.59446e-02; pars_fPmeanVsBGcorr[1][i][1]=-9.09587e-01; pars_fPmeanVsBGcorr[1][i][2]=1.00472e+00;
      }
      else if(i==1) {//K
	pars_fPmeanVsBGcorr[0][i][0]=3.00754e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.18175e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97758e-01;
	pars_fPmeanVsBGcorr[1][i][0]=3.36764e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.08206e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00094e+00;
      }
      else if(i==2) {//p
	pars_fPmeanVsBGcorr[0][i][0]=1.54832e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.70549e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97921e-01;
	pars_fPmeanVsBGcorr[1][i][0]=1.75353e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.52898e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00121e+00;
      }
      else if(i==3) {//d
	pars_fPmeanVsBGcorr[0][i][0]=0.010962; pars_fPmeanVsBGcorr[0][i][1]=-2.440895; pars_fPmeanVsBGcorr[0][i][2]=0.997846;
	pars_fPmeanVsBGcorr[1][i][0]=0.012415; pars_fPmeanVsBGcorr[1][i][1]=-2.281648; pars_fPmeanVsBGcorr[1][i][2]=1.001189;
      }
      else if(i==4) {//He3
	pars_fPmeanVsBGcorr[0][i][0]=0.030140; pars_fPmeanVsBGcorr[0][i][1]=-2.169590; pars_fPmeanVsBGcorr[0][i][2]=0.997846;
	pars_fPmeanVsBGcorr[1][i][0]=0.034135; pars_fPmeanVsBGcorr[1][i][1]=-2.028043; pars_fPmeanVsBGcorr[1][i][2]=1.001189;
      }
    }
  }
  else if(etaMin>0.49999 && etaMax<0.60001) {
    //printf("EtaLimit[0]== %f and EtaLimit[1]== %fFFF\n",EtaLimit[0],EtaLimit[1]);
    for(Int_t i=0;i<5;i++) {
      if(i==0) {//pi
	pars_fPmeanVsBGcorr[0][i][0]=5.29436e-02; pars_fPmeanVsBGcorr[0][i][1]=-5.04070e-01; pars_fPmeanVsBGcorr[0][i][2]=1.00951e+00;
	pars_fPmeanVsBGcorr[1][i][0]=1.04356e-01; pars_fPmeanVsBGcorr[1][i][1]=-1.19297e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00197e+00;
      }
      else if(i==1) {//K
	pars_fPmeanVsBGcorr[0][i][0]=3.36237e-02; pars_fPmeanVsBGcorr[0][i][1]=-1.89739e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97921e-01;
	pars_fPmeanVsBGcorr[1][i][0]=3.76386e-02; pars_fPmeanVsBGcorr[1][i][1]=-1.89484e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00097e+00;
      }
      else if(i==2) {//p
	pars_fPmeanVsBGcorr[0][i][0]=1.93889e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.38744e+00; pars_fPmeanVsBGcorr[0][i][2]=9.98551e-01;
	pars_fPmeanVsBGcorr[1][i][0]=2.12666e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.29606e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00174e+00;
      }
      else if(i==3) {//d
	pars_fPmeanVsBGcorr[0][i][0]=0.013727; pars_fPmeanVsBGcorr[0][i][1]=-2.153951; pars_fPmeanVsBGcorr[0][i][2]=0.998275;
	pars_fPmeanVsBGcorr[1][i][0]=0.015057; pars_fPmeanVsBGcorr[1][i][1]=-2.071511; pars_fPmeanVsBGcorr[1][i][2]=1.001396;
      }
      else if(i==4) {//He3
	pars_fPmeanVsBGcorr[0][i][0]=0.037743; pars_fPmeanVsBGcorr[0][i][1]=-1.914539; pars_fPmeanVsBGcorr[0][i][2]=0.998275;
	pars_fPmeanVsBGcorr[1][i][0]=0.041398; pars_fPmeanVsBGcorr[1][i][1]=-1.841262; pars_fPmeanVsBGcorr[1][i][2]=1.001396;
      }
    }
  }
  else if(etaMin>0.59999 && etaMax<0.70001) {
    //printf("EtaLimit[0]== %f and EtaLimit[1]== %fGGG\n",EtaLimit[0],EtaLimit[1]);
    for(Int_t i=0;i<5;i++) {
      if(i==0) {//pi
	pars_fPmeanVsBGcorr[0][i][0]=7.18462e-02; pars_fPmeanVsBGcorr[0][i][1]=-1.15676e+00; pars_fPmeanVsBGcorr[0][i][2]=9.99111e-01;
	pars_fPmeanVsBGcorr[1][i][0]=8.52428e-02; pars_fPmeanVsBGcorr[1][i][1]=-9.11048e-01; pars_fPmeanVsBGcorr[1][i][2]=1.00777e+00;
      }
      else if(i==1) {//K
	pars_fPmeanVsBGcorr[0][i][0]=3.32328e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.30015e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97683e-01;
	pars_fPmeanVsBGcorr[1][i][0]=3.88555e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.18109e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00130e+00;
      }
      else if(i==2) {//p
	pars_fPmeanVsBGcorr[0][i][0]=1.71488e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.73349e+00; pars_fPmeanVsBGcorr[0][i][2]=9.98517e-01;
	pars_fPmeanVsBGcorr[1][i][0]=1.89105e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.65564e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00159e+00;
      }
      else if(i==3) {//d
	pars_fPmeanVsBGcorr[0][i][0]=0.012141; pars_fPmeanVsBGcorr[0][i][1]=-2.466162; pars_fPmeanVsBGcorr[0][i][2]=0.998125;
	pars_fPmeanVsBGcorr[1][i][0]=0.013389; pars_fPmeanVsBGcorr[1][i][1]=-2.395927; pars_fPmeanVsBGcorr[1][i][2]=1.001633;
      }
      else if(i==4) {//He3
	pars_fPmeanVsBGcorr[0][i][0]=0.033383; pars_fPmeanVsBGcorr[0][i][1]=-2.192048; pars_fPmeanVsBGcorr[0][i][2]=0.998125;
	pars_fPmeanVsBGcorr[1][i][0]=0.036812; pars_fPmeanVsBGcorr[1][i][1]=-2.129620; pars_fPmeanVsBGcorr[1][i][2]=1.001633;
      }
    }
  }
  else if(etaMin>0.69999 && etaMax<0.80001) {
    //printf("EtaLimit[0]== %f and EtaLimit[1]== %fHHH\n",EtaLimit[0],EtaLimit[1]);
    for(Int_t i=0;i<5;i++) {
      if(i==0) {//pi
	pars_fPmeanVsBGcorr[0][i][0]=9.56419e-02; pars_fPmeanVsBGcorr[0][i][1]=-1.31941e+00; pars_fPmeanVsBGcorr[0][i][2]=9.98375e-01;
	pars_fPmeanVsBGcorr[1][i][0]=8.30340e-02; pars_fPmeanVsBGcorr[1][i][1]=-4.46775e-01; pars_fPmeanVsBGcorr[1][i][2]=1.02721e+00;
      }
      else if(i==1) {//K
	pars_fPmeanVsBGcorr[0][i][0]=3.55532e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.25782e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97746e-01;
	pars_fPmeanVsBGcorr[1][i][0]=4.26998e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.10431e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00185e+00;
      }
      else if(i==2) {//p
	pars_fPmeanVsBGcorr[0][i][0]=1.87103e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.65814e+00; pars_fPmeanVsBGcorr[0][i][2]=9.98847e-01;
	pars_fPmeanVsBGcorr[1][i][0]=2.07010e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.60124e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00210e+00;
      }
      else if(i==3) {//d
	pars_fPmeanVsBGcorr[0][i][0]=0.013247; pars_fPmeanVsBGcorr[0][i][1]=-2.398177; pars_fPmeanVsBGcorr[0][i][2]=0.998269;
	pars_fPmeanVsBGcorr[1][i][0]=0.014656; pars_fPmeanVsBGcorr[1][i][1]=-2.346847; pars_fPmeanVsBGcorr[1][i][2]=1.002033;
      }
      else if(i==4) {//He3
	pars_fPmeanVsBGcorr[0][i][0]=0.036422; pars_fPmeanVsBGcorr[0][i][1]=-2.131620; pars_fPmeanVsBGcorr[0][i][2]=0.998269;
	pars_fPmeanVsBGcorr[1][i][0]=0.040298; pars_fPmeanVsBGcorr[1][i][1]=-2.085995; pars_fPmeanVsBGcorr[1][i][2]=1.002033;
      }
    }
  }
  else {//for all eta
    //printf("EtaLimit[0]== %f and EtaLimit[1]== %fIII\n",EtaLimit[0],EtaLimit[1]);
    for(Int_t i=0;i<5;i++) {
      if(i==0) {//pi
	pars_fPmeanVsBGcorr[0][i][0]=4.89956e-02; pars_fPmeanVsBGcorr[0][i][1]=-6.46308e-01; pars_fPmeanVsBGcorr[0][i][2]=1.00462e+00;
	pars_fPmeanVsBGcorr[1][i][0]=6.36672e-02; pars_fPmeanVsBGcorr[1][i][1]=-6.10966e-01; pars_fPmeanVsBGcorr[1][i][2]=1.01188e+00;
      }
      else if(i==1) {//K
	pars_fPmeanVsBGcorr[0][i][0]=3.06216e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.10247e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97142e-01;
	pars_fPmeanVsBGcorr[1][i][0]=3.48865e-02; pars_fPmeanVsBGcorr[1][i][1]=-1.89213e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00123e+00;
      }
      else if(i==2) {//p
	pars_fPmeanVsBGcorr[0][i][0]=1.58652e-02; pars_fPmeanVsBGcorr[0][i][1]=-2.64898e+00; pars_fPmeanVsBGcorr[0][i][2]=9.97176e-01;
	pars_fPmeanVsBGcorr[1][i][0]=1.83264e-02; pars_fPmeanVsBGcorr[1][i][1]=-2.45858e+00; pars_fPmeanVsBGcorr[1][i][2]=1.00079e+00;
      }
      else if(i==3) {//d
	pars_fPmeanVsBGcorr[0][i][0]=0.011233; pars_fPmeanVsBGcorr[0][i][1]=-2.389911; pars_fPmeanVsBGcorr[0][i][2]=0.997176;//0.997210
	pars_fPmeanVsBGcorr[1][i][0]=0.012975; pars_fPmeanVsBGcorr[1][i][1]=-2.218137; pars_fPmeanVsBGcorr[1][i][2]=1.001091;
      }
      else if(i==4) {//He3
	pars_fPmeanVsBGcorr[0][i][0]=0.030884; pars_fPmeanVsBGcorr[0][i][1]=-2.124273; pars_fPmeanVsBGcorr[0][i][2]=0.997176;//0.997210
	pars_fPmeanVsBGcorr[1][i][0]=0.035675; pars_fPmeanVsBGcorr[1][i][1]=-1.971591; pars_fPmeanVsBGcorr[1][i][2]=1.001091;
      }
    }
  }

  /*
  for(Int_t iB=0;iB<nBconf;iB++) {
    for(Int_t i=0;i<5;i++) {
      pars_fPmeanVsBGcorr[iB][i][0]=0.02; pars_fPmeanVsBGcorr[iB][i][1]=-2.0; pars_fPmeanVsBGcorr[iB][i][2]=1.0;
    }
    }*/

  for(Int_t iB=0;iB<nBconf;iB++) {
    for(Int_t i=5;i<10;i++) {
      pars_fPmeanVsBGcorr[iB][i][0]=pars_fPmeanVsBGcorr[iB][i-5][0];
      pars_fPmeanVsBGcorr[iB][i][1]=pars_fPmeanVsBGcorr[iB][i-5][1];
      pars_fPmeanVsBGcorr[iB][i][2]=pars_fPmeanVsBGcorr[iB][i-5][2];
    }
  }

  for(Int_t iB=0;iB<nBconf;iB++) {
    for(Int_t i=0;i<10;i++) {
      fPmeanVsBGcorr[iB][i]->SetParameters(pars_fPmeanVsBGcorr[iB][i]);
      fPmeanVsBGcorr[iB][i]->SetNpx(fPmeanVsBGcorr[iB][i]->GetNpx()*10);
    }
  }
    
  return;
  
}

