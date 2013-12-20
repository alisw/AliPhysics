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

ClassImp(AliAnalysisNucleiMass)//...

//_____________________________________________________________________________
AliAnalysisNucleiMass::AliAnalysisNucleiMass():
  AliAnalysisTaskSE(),
//Centrality(NULL),                         
  FilterBit(16),                                
//  EtaLimit(NULL),                           
  DCAxyCut(0.1),                               
  DCAzCut(1000.),                                
  NsigmaTpcCut(2.0),                           
  NminTpcCluster(0),                           
  iTrdCut(0),
  kSignalCheck(1),
  iMtof(1),
  iBconf(0),                                  
  kTOF(0),               
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL)
{
  fList[0]=new TList();
  fList[0]->SetName("results");
  
  fList[1]=new TList();
  fList[1]->SetName("results2");
}
//______________________________________________________________________________
AliAnalysisNucleiMass::AliAnalysisNucleiMass(const char *name):
  AliAnalysisTaskSE(name),
  //Centrality(NULL),                         
  FilterBit(16),                                
  //EtaLimit(NULL),                           
  DCAxyCut(0.1),                               
  DCAzCut(1000.),                                
  NsigmaTpcCut(2.0),                           
  NminTpcCluster(0),
  iTrdCut(0),
  kSignalCheck(1),
  iMtof(1),
  iBconf(0),                                  
  kTOF(0),               
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fPIDResponse(NULL)
{
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
  
  Double_t binPt[nbin+1];
  for(Int_t i=0;i<nbin+1;i++) {
    binPt[i]=0.4+0.1*i;
  }
  
  Char_t name_nbin[nbin][200];
  for(Int_t j=0;j<nbin;j++) {
    snprintf(name_nbin[j],200,"%.1f<Pt<%.1f",binPt[j],binPt[j+1]);
  }
  
  for(Int_t iB=0;iB<nBconf;iB++) {
    
    htemp[iB] = new TH1F("htemp","htemp (avoid the problem with the empty list...);B field",20,-10,10);

    hCentrality[iB][0] = new TH1F("hCentrality_Selected","Centrality (selected events);centrality(%)",20,0,100);
    hCentrality[iB][1] = new TH1F("hCentrality_Analyzed","Centrality (analyzed events);centrality (%)",20,0,100);
    
    hZvertex[iB][0] = new TH1F("hZvertex_Selected","Vertex distribution of selected events;z vertex (cm)",240,-30,30);
    hZvertex[iB][1] = new TH1F("hZvertex_Analyzed","Vertex distribution of analyzed events;z vertex (cm)",240,-30,30);
    
    hEta[iB] = new TH1F("hEta_Analyzed","|#eta| distribution after the track cuts;|#eta|",100,0.0,1.0);
    
    hPhi[iB] = new TH1F("hPhi_Analyzed","#phi distribution after the track cuts;#phi (rad.)",90,0,6.3);//Each TRD supermodule is divided for 5 (DeltaPhi(TRD)=0.35 theoretical)
    
    Int_t hbins[2];
    if(kSignalCheck>1) {hbins[0]=100; hbins[1]=90;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    fEtaPhi[iB] = new TH2F("fEtaPhi_Analyzed","|#eta| vs. #phi after the track cuts;|#eta|;#phi (rad.)",hbins[0],0.0,1.0,hbins[1],0,6.3);

    hNTpcCluster[iB] = new TH1F("hNTpcCluster","Number of the TPC clusters after the track cuts;n_{cl}^{TPC}",300,0,300);

    hNTrdSlices[iB] = new TH1F("hNTrdSlices","Number of the TRD slices after the track cuts;n_{slices}^{TRD}",40,0,40);

    if(kSignalCheck==1) {hbins[0]=500; hbins[1]=2000;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=100; hbins[1]=500;}
    fdEdxVSp[iB][0] = new TH2F("fdEdxVSp_pos","dE/dx vs p (positive charge); p/|z| (GeV/c); dE/dx_{TPC} (a.u.)",hbins[0],0,5,hbins[1],0,1000);
    fdEdxVSp[iB][1] = new TH2F("fdEdxVSp_neg","dE/dx vs p (negative charge); p/|z| (GeV/c); dE/dx_{TPC} (a.u.)",hbins[0],0,5,hbins[1],0,1000);

    Char_t name_hDeDxExp[nPart][200];
    Char_t title_hDeDxExp[nPart][200];
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_hDeDxExp[i],200,"hDeDxExp_%s",namePart[i]);
      snprintf(title_hDeDxExp[i],200,"Expected dE/dx of %s in the TPC;p/|z| (GeV/c);dE/dx_{TPC} (a.u.)",namePart[i]);
      hDeDxExp[iB][i] = new TProfile(name_hDeDxExp[i],title_hDeDxExp[i],500,0,5,0,1000,"");
    }

    Char_t name_fNsigmaTpc[nPart][200];
    Char_t title_fNsigmaTpc[nPart][200];
    if(kSignalCheck==1) {hbins[0]=100; hbins[1]=100;}
    else {hbins[0]=1; hbins[1]=1;}
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_fNsigmaTpc[i],200,"NsigmaTpc_%s",namePart[i]);
      snprintf(title_fNsigmaTpc[i],200,"NsigmaTpc_%s;p_{TPC}/|z| (GeV/c);n_{#sigma_{TPC}}^{%s}",namePart[i],namePart[i]);
      fNsigmaTpc[iB][i] = new TH2F(name_fNsigmaTpc[i],title_fNsigmaTpc[i],hbins[0],0,5,hbins[1],-5,5);
    }
    
    if(kSignalCheck>1) {hbins[0]=100; hbins[1]=100;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    Char_t name_fNsigmaTpc_kTOF[nSpec][200];
    Char_t title_fNsigmaTpc_kTOF[nSpec][200];
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fNsigmaTpc_kTOF[i],200,"NsigmaTpc_%s_kTOF",name[i]);
      snprintf(title_fNsigmaTpc_kTOF[i],200,"NsigmaTpc_kTOF_%s in DCAxyCut;p_{T}/|z| (GeV/c);n_{#sigma_{TPC}}^{%s}",name[i],name[i]);
      fNsigmaTpc_kTOF[iB][i] = new TH2F(name_fNsigmaTpc_kTOF[i],title_fNsigmaTpc_kTOF[i],hbins[0],0,5,hbins[1],-5,5);
    }

    if(kSignalCheck==1) {hbins[0]=1000; hbins[1]=1300;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=100; hbins[1]=260;}
    fBetaTofVSp[iB][0] = new TH2F("fBetaTofVSp_pos","#beta_{TOF} vs p/|z| (positive charge);p(GeV/c);#beta_{TOF}",hbins[0],0,5,hbins[1],0.4,1.05);
    fBetaTofVSp[iB][1] = new TH2F("fBetaTofVSp_neg","#beta_{TOF} vs p/|z| (negative charge);p(GeV/c);#beta_{TOF}",hbins[0],0,5,hbins[1],0.4,1.05);
    
    Char_t name_hBetaExp[nPart][200];
    Char_t title_hBetaExp[nPart][200];
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_hBetaExp[i],200,"hBetaTofVsP_Exp_%s",namePart[i]);
      snprintf(title_hBetaExp[i],200,"Expected #beta_{TOF} vs p/|z| of %s;p/|z| (GeV/c); #beta_{TOF}",namePart[i]);
      hBetaExp[iB][i] = new TProfile(name_hBetaExp[i],title_hBetaExp[i],400,0,5,0.4,1.05,"");
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
    if(kSignalCheck>1) {hbins[0]=100; hbins[1]=100;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fNsigmaTof_DcaCut[i],200,"NsigmaTof_DcaCut_%s",name[i]);
      snprintf(title_fNsigmaTof_DcaCut[i],200,"NsigmaTof_%s with DCAxyCut;p_{T}/|z| (GeV/c);n_{#sigma_{TOF}}^{%s}",name[i],name[i]);
      fNsigmaTof_DcaCut[iB][i] = new TH2F(name_fNsigmaTof_DcaCut[i],title_fNsigmaTof_DcaCut[i],hbins[0],0,5,hbins[1],-5,5);
    }

    if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
    else {hbins[0]=1; hbins[1]=1;}
    fM2vsPt_NoTpcCut[iB][0][0] = new TH2F("fM2vsPt_NoTpcCut_pos","m^{2}/z^{2}_{TOF} vs p_{T}/|z| (positive charge);m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p_{T}/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,5);
    fM2vsPt_NoTpcCut[iB][0][1] = new TH2F("fM2vsPt_NoTpcCut_neg","m^{2}/z^{2}_{TOF} vs p_{T}/|z| (negative charge);m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p_{T}/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,5);

    if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=1000; hbins[1]=100;}
    fM2vsPt_NoTpcCut[iB][1][0] = new TH2F("fM2vsPt_NoTpcCut_DCAxyCut_pos","m^{2}/z^{2}_{TOF} vs p_{T}/|z| (positive charge) with DCAxy cut;m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p_{T}/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,5);
    fM2vsPt_NoTpcCut[iB][1][1] = new TH2F("fM2vsPt_NoTpcCut_DCAxyCut_neg","m^{2}/z^{2}_{TOF} vs p_{T}/|z| (negative charge) with DCAxy cut;m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p_{T}/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,5);

    Char_t name_fM2vsPt[2][18][300]; 
    Char_t title_fM2vsPt[2][18][300]; 

    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fM2vsPt[0][i],300,"fM2vsPt_%s",name[i]);
      snprintf(title_fM2vsPt[0][i],300,"m^{2}/z^{2}_{TOF} vs p_{T}/|z| of %s with a NsigmaTpcCut;m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p_{T}/|z| (GeV/c)",name[i]);
      
      snprintf(name_fM2vsPt[1][i],300,"fM2vsPt_%s_DCAxyCut",name[i]);
      snprintf(title_fM2vsPt[1][i],300,"m^{2}/z^{2}_{TOF} vs p_{T}/|z| of %s with a NsigmaTpcCut and with the DCAxy cut;m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4});p_{T}/|z| (GeV/c)",name[i]);

      if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
      else {hbins[0]=1; hbins[1]=1;}
      fM2vsPt[iB][0][i] = new TH2F(name_fM2vsPt[0][i],title_fM2vsPt[0][i],hbins[0],0,10,hbins[1],0,5);
      
      if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
      else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
      else if(kSignalCheck==2) {hbins[0]=1000; hbins[1]=100;}
      fM2vsPt[iB][1][i] = new TH2F(name_fM2vsPt[1][i],title_fM2vsPt[1][i],hbins[0],0,10,hbins[1],0,5);
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

    //Parameterizations:
    fPmeanVsPexp[0]=new TF1("fPmeanVsPexp_p","[2]-[0]*TMath::Exp(-(TMath::Max(x,[3])*[1]))",0,20);
    fPmeanVsPexp[1]=new TF1("fPmeanVsPexp_d","[2]-[0]*TMath::Exp(-(TMath::Max(x,[3])*[1]))",0,20);
    fPmeanVsPexp[2]=new TF1("fPmeanVsPexp_He3","[2]-[0]*TMath::Exp(-(TMath::Max(x,[3])*[1]))",0,20);
           
    Double_t fpars_p[4]={5.14500484596484148e-03,9.74729863202270397e-01,0.0,1.00607413672776569e+00};
    Double_t fpars_d[4]={3.16023942908439243e-02,1.24005027514358490e+00,-1.50000000000000003e-03,1.40607413672776560e+00};
    Double_t fpars_He3[4]={2.73329079591698026e-02,1.53005942367188852e+00,-4.10231310888738848e-03,1.20607413672776564e+00};
    fPmeanVsPexp[0]->SetParameters(fpars_p);
    fPmeanVsPexp[1]->SetParameters(fpars_d);
    fPmeanVsPexp[2]->SetParameters(fpars_He3);

    /*Char_t title_Xaxis[3][200];
    Char_t title_Yaxis[3][200];
    snprintf(title_Xaxis[0],200,"p(t_{exp}^{%s})",namePart[4]);
    snprintf(title_Yaxis[0],200,"p(t_{TOF})-p(t_{exp}^{%s})/p(t_{exp}^{%s})",namePart[4],namePart[4]);
    snprintf(title_Xaxis[1],200,"p(t_{exp}^{%s})",namePart[5]);
    snprintf(title_Yaxis[1],200,"p(t_{TOF})-p(t_{exp}^{%s})/p(t_{exp}^{%s})",namePart[5],namePart[5]);
    snprintf(title_Xaxis[2],200,"p(t_{exp}^{%s})",namePart[7]);
    snprintf(title_Yaxis[2],200,"p(t_{TOF})-p(t_{exp}^{%s})/p(t_{exp}^{%s})",namePart[7],namePart[7]);
    for(Int_t i=0;i<3;i++){
      fPmeanVsPexp[i]->GetXaxis()->SetTitle(title_Xaxis[i]);
      fPmeanVsPexp[i]->GetYaxis()->SetTitle(title_Yaxis[i]);
      fPmeanVsPexp[i]->SetTitle("Parameterization calculated with Monte Carlo (LHC13d15)");
      }*/
    //end parameterizations

    Char_t name_fPmeanVsBetaGamma[18][200];
    Char_t title_fPmeanVsBetaGamma[18][200];
    
    hbins[0]=200; hbins[1]=200;
    for(Int_t iS=0;iS<nSpec;iS++) {
      snprintf(name_fPmeanVsBetaGamma[iS],200,"fPmeanVsPvtx_%s",name[iS]);
      snprintf(title_fPmeanVsBetaGamma[iS],200,"<p>/p_{vtx} vs #beta#gamma of %s;p_{vtx}/m_{%s};<p>_{%s}/p_{vtx}",name[iS],name[iS],name[iS]);
      fPmeanVsBetaGamma[iB][iS]=new TH2F(name_fPmeanVsBetaGamma[iS],title_fPmeanVsBetaGamma[iS],hbins[0],0,10,hbins[1],0.8,1.2);
    }	
    
    Char_t name_prPmeanVsBetaGamma[18][200];
    Char_t title_prPmeanVsBetaGamma[18][200];
    
    for(Int_t iS=0;iS<nSpec;iS++) {
      snprintf(name_prPmeanVsBetaGamma[iS],200,"prPmeanVsPvtx_%s",name[iS]);
      snprintf(title_prPmeanVsBetaGamma[iS],200,"<p>/p_{vtx} vs #beta#gamma of %s;p_{vtx}/m_{%s};<p>_{%s}/p_{vtx}",name[iS],name[iS],name[iS]);
      prPmeanVsBetaGamma[iB][iS]=new TProfile(name_prPmeanVsBetaGamma[iS],title_prPmeanVsBetaGamma[iS],hbins[0],0,10,0.8,1.2,"");
    }	

    fList[iB]->Add(htemp[iB]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(hCentrality[iB][i]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(hZvertex[iB][i]);
    fList[iB]->Add(hEta[iB]);
    fList[iB]->Add(hPhi[iB]);
    fList[iB]->Add(fEtaPhi[iB]);
    fList[iB]->Add(hNTpcCluster[iB]);
    fList[iB]->Add(hNTrdSlices[iB]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fdEdxVSp[iB][i]);
    for(Int_t i=0;i<nPart;i++) fList[iB]->Add(hDeDxExp[iB][i]);
    for(Int_t i=0;i<nPart;i++) fList[iB]->Add(fNsigmaTpc[iB][i]);
    for(Int_t i=0;i<nPart;i++) {
      if(kSignalCheck!=1) 
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      fList[iB]->Add(fNsigmaTpc_kTOF[iB][i]);
      fList[iB]->Add(fNsigmaTpc_kTOF[iB][i+nPart]);
    }
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fBetaTofVSp[iB][i]);
    for(Int_t i=0;i<nPart;i++) fList[iB]->Add(hBetaExp[iB][i]);
    for(Int_t i=0;i<nPart;i++) fList[iB]->Add(fNsigmaTof[iB][i]);
    for(Int_t i=0;i<nPart;i++) {
      if(kSignalCheck!=1) 
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      fList[iB]->Add(fNsigmaTof_DcaCut[iB][i]);
      fList[iB]->Add(fNsigmaTof_DcaCut[iB][i+nPart]);
    }
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fM2vsPt_NoTpcCut[iB][0][i]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fM2vsPt_NoTpcCut[iB][1][i]);
    for(Int_t i=0;i<nPart;i++) {
      if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      fList[iB]->Add(fM2vsPt[iB][0][i]);
      fList[iB]->Add(fM2vsPt[iB][0][i+nPart]);
    }
    for(Int_t i=0;i<nPart;i++){
      if(kSignalCheck!=1) 
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      fList[iB]->Add(fM2vsPt[iB][1][i]);
      fList[iB]->Add(fM2vsPt[iB][1][i+nPart]);
    }
    if(iMtof!=1) {
      for(Int_t i=0;i<nPart;i++){
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
	fList[iB]->Add(fPmeanVsBetaGamma[iB][i]);
	fList[iB]->Add(prPmeanVsBetaGamma[iB][i]);
	fList[iB]->Add(fPmeanVsBetaGamma[iB][i+nPart]);
	fList[iB]->Add(prPmeanVsBetaGamma[iB][i+nPart]);
      }
    }
    if(iMtof==8) for(Int_t i=0;i<3;i++) fList[iB]->Add(fPmeanVsPexp[i]);
    else if(iMtof==4) for(Int_t i=1;i<3;i++) fList[iB]->Add(fPmeanVsPexp[i]);
    for(Int_t i=0;i<10;i++) fList[iB]->Add(fM2vsZ[iB][i]);
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
      Double_t etaAbs = TMath::Abs(track->Eta());
      
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
      if ((track->Pt() < 0.2) || (etaAbs<EtaLimit[0]) || (etaAbs>EtaLimit[1]) || !trkFlag || !isMinTpcCluster)
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

      //Track info:
      Double_t phi= track->Phi();
      
      hEta[iBconf]->Fill(etaAbs);
      hPhi[iBconf]->Fill(phi);
      fEtaPhi[iBconf]->Fill(etaAbs,phi);
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
	expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
	hDeDxExp[iBconf][iS]->Fill(pTPC,expdedx[iS]);
	nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
	fNsigmaTpc[iBconf][iS]->Fill(pTPC,nsigmaTPC[iS]);
	if(charge>0) {//positive particle
	  if(kTOF && (TMath::Abs(DCAxy)<DCAxyCut)) fNsigmaTpc_kTOF[iBconf][iS]->Fill(pt,nsigmaTPC[iS]);
	}
	else {//negative particle
	  if(kTOF && (TMath::Abs(DCAxy)<DCAxyCut)) fNsigmaTpc_kTOF[iBconf][iS+nPart]->Fill(pt,nsigmaTPC[iS]);
	}

	//TPC identification:
	if(TMath::Abs(nsigmaTPC[iS])<NsigmaTpcCut) {
	  FlagPid += ((Int_t)TMath::Power(2,iS));
	}
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
	}
	if(charge>0) fBetaTofVSp[iBconf][0]->Fill(p,beta);
	else fBetaTofVSp[iBconf][1]->Fill(p,beta);
		
	this->GetMassFromPvertex(beta,p,M2);
	this->GetZTpc(dedx,pTPC,M2,Z2);
		
	//-----------------------------M2 as a function of expected times, if iMtof>1---------------------------------
	Double_t Mass2[9];
	if(iMtof>1) this->GetMassFromExpTimes(beta,exptimes,Mass2,iMtof,p,FlagPid,charge);
		
	//-------------------------------Squared Mass TH2 distributions-----------------------
	if(charge>0) {
	  //without TPC
	  fM2vsPt_NoTpcCut[iBconf][0][0]->Fill(M2,pt);
	  if(TMath::Abs(DCAxy)<DCAxyCut) fM2vsPt_NoTpcCut[iBconf][1][0]->Fill(M2,pt);
	  //with TPC
	  for(Int_t iS=0;iS<9;iS++) {
	    //-----------------------------M2 as a function of expected times, if iMtof>1---------------------------------
	    if(iMtof>1) {
	      M2=999.9;
	      M2=Mass2[iS];
	    }
	    //-----------------
	    if(FlagPid & stdFlagPid[iS]) {
	      fM2vsPt[iBconf][0][iS]->Fill(M2,pt);
	      if(TMath::Abs(DCAxy)<DCAxyCut) fM2vsPt[iBconf][1][iS]->Fill(M2,pt);
	    }
	  }
	}
	else {//charge<0
	  //without TPC
	  fM2vsPt_NoTpcCut[iBconf][0][1]->Fill(M2,pt);
	  if(TMath::Abs(DCAxy)<DCAxyCut) fM2vsPt_NoTpcCut[iBconf][1][1]->Fill(M2,pt);
	   //with TPC
	  for(Int_t iS=0;iS<9;iS++) {
	    //-----------------------------M2 as a function of expected times, if iMtof>1---------------------------------
	    if(iMtof>1) {
	      M2=999.9;
	      M2=Mass2[iS];
	    }
	    //-----------------
	    if(FlagPid & stdFlagPid[iS]) {
	      fM2vsPt[iBconf][0][iS+nPart]->Fill(M2,pt);
	      if(TMath::Abs(DCAxy)<DCAxyCut) fM2vsPt[iBconf][1][iS+nPart]->Fill(M2,pt);
	    }
	  }
	}
	
	//------------------------------start DCA and Squared Mass TH1 distributions-------------------------
	Double_t binPt[nbin+1];
	for(Int_t i=0;i<nbin+1;i++) {
	  binPt[i]=0.4+i*0.1;
	}

	if(charge>0) {
	  for(Int_t iS=0;iS<9;iS++) {
	    //-----------------------------M2 as a function of expected times, if iMtof>1---------------------------------
	    if(iMtof>1) {
	      M2=999.9;
	      M2=Mass2[iS];
	    }
	    //-----------------
	    if(FlagPid & stdFlagPid[iS]) {
	      for(Int_t j=0;j<nbin;j++) {
		if(pt>binPt[j] && pt<binPt[j+1]) {
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
	      }//end loop on the pT bins (j)
	    }
	  }//end loop on the particle species (iS)
	}
	else {//charge<0
	  for(Int_t iS=0;iS<9;iS++) {
	    //-----------------------------M2 as a function of expected times, if iMtof>1---------------------------------
	    if(iMtof>1) {
	      M2=999.9;
	      M2=Mass2[iS];
	    }
	    //-----------------
	    if(FlagPid & stdFlagPid[iS]) {
	      for(Int_t j=0;j<nbin;j++) {
		if(pt>binPt[j] && pt<binPt[j+1]) {
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
	      }//end loop on the pT bins (j)
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
void AliAnalysisNucleiMass::GetMassFromPvertex(Double_t beta, Double_t p, Double_t &M2) {
  
  M2 = p*p*(1-beta*beta)/(beta*beta);

  return;
  
}
//____________________________________________________________________________________________________________
void AliAnalysisNucleiMass::GetMassFromExpTimes(Double_t beta, Double_t *IntTimes, Double_t *Mass2, Int_t iCorr, Double_t pVtx, Int_t FlagPid, Double_t charge) {
 
  // m = p_exp/beta/gamma where p_exp = mPDG*beta_exp*gamma_exp; beta_exp = L/t_exp/c = t_e/t_exp ; beta=L/tof/c = t_e/tof
  // In this way m_tof = mPDG only if tof=t_exp
  
  Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};
    
  Double_t beta2Exp[9];
  Double_t p2Exp[9];
  
  Double_t pExp[9];
  Double_t CorrFactor=0.0;

  Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,#mu,#pi,K,p,d,t,3He,4He

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
    pExp[iS]=TMath::Sqrt(p2Exp[iS]);
    
    CorrFactor=0.0;
    if(iCorr & 12) {//iCorr==4 || iCorr==8
      if(iCorr==8 && iS==4) CorrFactor=fPmeanVsPexp[0]->Eval(pExp[iS]);
      
      if(iS==5) CorrFactor=fPmeanVsPexp[1]->Eval(pExp[iS]);
      else if(iS==7) CorrFactor=fPmeanVsPexp[2]->Eval(pExp[iS]);
      CorrFactor=pExp[iS]*CorrFactor;
      pExp[iS]=pExp[iS]+CorrFactor;//CorrFactor is negative so pExp(Corrected)<pExp
    }
    p2Exp[iS]=pExp[iS]*pExp[iS];
    //------------
    Mass2[iS]=p2Exp[iS]*(1-beta*beta)/(beta*beta);
  
    //------------
    if(FlagPid & stdFlagPid[iS]) {
      if(charge>0){
	fPmeanVsBetaGamma[iBconf][iS]->Fill(pVtx/massOverZ[iS],pExp[iS]/pVtx);
	prPmeanVsBetaGamma[iBconf][iS]->Fill(pVtx/massOverZ[iS],pExp[iS]/pVtx);
      }
      else {
	fPmeanVsBetaGamma[iBconf][iS+nPart]->Fill(pVtx/massOverZ[iS],pExp[iS]/pVtx);
	prPmeanVsBetaGamma[iBconf][iS+nPart]->Fill(pVtx/massOverZ[iS],pExp[iS]/pVtx);
      }
    }
  }//end loop on the particle species
  
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
