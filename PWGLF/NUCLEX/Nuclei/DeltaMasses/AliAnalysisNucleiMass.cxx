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
//iTriggerSel(-99),
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
  //iTriggerSel(-99),
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

    //htriggerbits[iB] = new TH1I("htriggerbits","htriggerbits; bits",10,-5,5);
    htriggerbits[iB][0] = new TH1I("htriggerbits_0","trigger mask; bits",45,-5,40);
    htriggerbits[iB][1] = new TH1I("htriggerbits_1","trigger bits (exclusive); bits",45,-5,40);
    
    hCentrality[iB][0] = new TH1F("hCentrality_Selected","Centrality (selected events);centrality(%)",20,0,100);//20,0,100
    hCentrality[iB][1] = new TH1F("hCentrality_Analyzed","Centrality (analyzed events);centrality (%)",20,0,100);//20,0,100
    
    hZvertex[iB][0] = new TH1F("hZvertex_Selected","Vertex distribution of selected events;z vertex (cm)",240,-30,30);
    hZvertex[iB][1] = new TH1F("hZvertex_Analyzed","Vertex distribution of analyzed events;z vertex (cm)",240,-30,30);
    
    hEta[iB] = new TH1F("hEta_Analyzed","|#eta| distribution after the track cuts;#eta",200,-1.0,1.0);
    
    hPhi[iB] = new TH1F("hPhi_Analyzed","#phi distribution after the track cuts;#phi (rad.)",90,0,6.3);//Each TRD supermodule is divided for 5 (DeltaPhi(TRD)=0.35 theoretical)
    
    Int_t hbins[2];
    if(kSignalCheck!=0) {hbins[0]=1; hbins[1]=1;}//{hbins[0]=100; hbins[1]=90;} to reduce RAM consuming (toram)
    else {hbins[0]=1; hbins[1]=1;}
    fEtaPhi[iB] = new TH2F("fEtaPhi_Analyzed","#eta vs. #phi after the track cuts;#eta;#phi (rad.)",hbins[0],-1.0,1.0,hbins[1],0,6.3);

    hNTpcCluster[iB] = new TH1F("hNTpcCluster","Number of the TPC clusters after the track cuts;n_{cl}^{TPC}",300,0,300);

    hNTrdSlices[iB] = new TH1F("hNTrdSlices","Number of the TRD slices after the track cuts;n_{slices}^{TRD}",40,0,40);

    if(kSignalCheck==1) {hbins[0]=500; hbins[1]=2000;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=200; hbins[1]=1000;}//{hbins[0]=100; hbins[1]=500;} toram
    fdEdxVSp[iB][0] = new TH2F("fdEdxVSp_pos","dE/dx vs p (positive charge); p/|z| (GeV/c); dE/dx_{TPC} (a.u.)",hbins[0],0,5,hbins[1],0,2000);
    fdEdxVSp[iB][1] = new TH2F("fdEdxVSp_neg","dE/dx vs p (negative charge); p/|z| (GeV/c); dE/dx_{TPC} (a.u.)",hbins[0],0,5,hbins[1],0,2000);

    Char_t name_hDeDxExp[nPart][200];
    Char_t title_hDeDxExp[nPart][200];
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_hDeDxExp[i],200,"hDeDxExp_%s",namePart[i]);
      snprintf(title_hDeDxExp[i],200,"Expected dE/dx of %s in the TPC;p/|z| (GeV/c);dE/dx_{TPC} (a.u.)",namePart[i]);
      hDeDxExp[iB][i] = new TProfile(name_hDeDxExp[i],title_hDeDxExp[i],200,0,5,0,2000,"");//,500,0,5,0,1000,""); toram
    }

    Char_t name_fNsigmaTpc[nSpec][200];
    Char_t title_fNsigmaTpc[nSpec][200];
    if(kSignalCheck==1) {hbins[0]=100; hbins[1]=100;}//{hbins[0]=100; hbins[1]=100;} toram
    else {hbins[0]=100; hbins[1]=200;}//temp!
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fNsigmaTpc[i],200,"NsigmaTpc_%s",name[i]);
      snprintf(title_fNsigmaTpc[i],200,"NsigmaTpc_%s;p_{TPC}/|z| (GeV/c);n_{#sigma_{TPC}}^{%s}",name[i],name[i]);
      fNsigmaTpc[iB][i] = new TH2F(name_fNsigmaTpc[i],title_fNsigmaTpc[i],hbins[0],0,5,hbins[1],-10,10);
    }
    
    if(kSignalCheck>0) {hbins[0]=1; hbins[1]=1;}//{hbins[0]=100; hbins[1]=100;} toram
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    Char_t name_fNsigmaTpc_kTOF[nSpec][200];
    Char_t title_fNsigmaTpc_kTOF[nSpec][200];
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fNsigmaTpc_kTOF[i],200,"NsigmaTpc_%s_kTOF",name[i]);
      snprintf(title_fNsigmaTpc_kTOF[i],200,"NsigmaTpc_kTOF_%s;p/|z| (GeV/c);n_{#sigma_{TPC}}^{%s}",name[i],name[i]);
      fNsigmaTpc_kTOF[iB][i] = new TH2F(name_fNsigmaTpc_kTOF[i],title_fNsigmaTpc_kTOF[i],hbins[0],0,5,hbins[1],-5,5);
    }

    if(kSignalCheck==1) {hbins[0]=1000; hbins[1]=1300;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=200; hbins[1]=260;}//{hbins[0]=100; hbins[1]=260;}
    fBetaTofVSp[iB][0] = new TH2F("fBetaTofVSp_pos","#beta_{TOF} vs p/|z| (positive charge);p(GeV/c);#beta_{TOF}",hbins[0],0,10,hbins[1],0.4,1.05);
    fBetaTofVSp[iB][1] = new TH2F("fBetaTofVSp_neg","#beta_{TOF} vs p/|z| (negative charge);p(GeV/c);#beta_{TOF}",hbins[0],0,10,hbins[1],0.4,1.05);
    
    Char_t name_hBetaExp[nPart][200];
    Char_t title_hBetaExp[nPart][200];
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_hBetaExp[i],200,"hBetaTofVsP_Exp_%s",namePart[i]);
      snprintf(title_hBetaExp[i],200,"Expected #beta_{TOF} vs p/|z| of %s;p/|z| (GeV/c); #beta_{TOF}",namePart[i]);
      hBetaExp[iB][i] = new TProfile(name_hBetaExp[i],title_hBetaExp[i],200,0,10,0.4,1.05,"");//,400,0,5,0.4,1.05,""); toram
    }
    
    Char_t name_fNsigmaTof[nPart][200];
    Char_t title_fNsigmaTof[nPart][200];    
    if(kSignalCheck==1) {hbins[0]=100; hbins[1]=100;}
    else {hbins[0]=100; hbins[1]=200;}
    for(Int_t i=0;i<nPart;i++) {
      snprintf(name_fNsigmaTof[i],200,"NsigmaTof_%s",namePart[i]);
      snprintf(title_fNsigmaTof[i],200,"NsigmaTof_%s;p_{T}/|z| (GeV/c);n_{#sigma_{TOF}}^{%s}",namePart[i],namePart[i]);
      fNsigmaTof[iB][i] = new TH2F(name_fNsigmaTof[i],title_fNsigmaTof[i],hbins[0],0,5,hbins[1],-10,10);
    }

    if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=800; hbins[1]=200;}// {hbins[0]=1000; hbins[1]=100;} toram
    fM2vsP_NoTpcCut[iB][0][0] = new TH2F("fM2vsP_NoTpcCut_pos","(m/z)^{2}_{TOF} vs p/|z| (positive charge);(m/z)^{2}_{TOF} (GeV^{2}/c^{4});p/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,10);
    fM2vsP_NoTpcCut[iB][0][1] = new TH2F("fM2vsP_NoTpcCut_neg","(m/z)^{2}_{TOF} vs p/|z| (negative charge);(m/z)^{2}_{TOF} (GeV^{2}/c^{4});p/|z| (GeV/c)",hbins[0],0,10,hbins[1],0,10);

    Char_t name_fM2vsP[1][18][300]; 
    Char_t title_fM2vsP[1][18][300]; 

    if(kSignalCheck==1) {hbins[0]=8000; hbins[1]=100;}
    else {hbins[0]=1; hbins[1]=1;}
    for(Int_t i=0;i<nSpec;i++) {
      snprintf(name_fM2vsP[0][i],300,"fM2vsPc_%s",name[i]);
      snprintf(title_fM2vsP[0][i],300,"(m/z)^{2}_{TOF} vs p/|z| of %s with a NsigmaTpcCut (pReco->pTrue for nuclei (if kMomVtxCorr==1));(m/z)^{2}_{TOF} (GeV^{2}/c^{4});p/|z| (GeV/c)",name[i]);
      if(i==0 || i==0+9) fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],100,0.0,0.2,100,0,5);
      if(i==1 || i==1+9) fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],100,0.0,0.2,100,0,5);
      if(i==2 || i==2+9) fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],200,0.0,0.4,100,0,5);
      if(i==3 || i==3+9) fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],100,0.1,0.5,100,0,5);
      if(i==4 || i==4+9) fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],100,0.5,1.5,100,0,5);
      if(i==5 || i==5+9) fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],200,2.5,5,100,0,5);
      if(i==6 || i==6+9) fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],100,7,9,100,0,5);
      if(i==7 || i==7+9) fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],100,1,3,100,0,5);
      if(i==8 || i==8+9) fM2vsP[iB][0][i] = new TH2F(name_fM2vsP[0][i],title_fM2vsP[0][i],100,2.5,5,100,0,5);
    }
    
    if(kSignalCheck==1) {hbins[0]=4000; hbins[1]=1000;}
    else if(kSignalCheck==0) {hbins[0]=1; hbins[1]=1;}
    else if(kSignalCheck==2) {hbins[0]=4000; hbins[1]=1000;}//{hbins[0]=1000 oppure 500; hbins[1]=100;} toram
    fM2vsZ[iB][0] = new TH2F("fM2vsZ_0","(m/z)^{2}_{TOF} vs z_{TPC} over all p/z (dEdx>80, DCAxy<0.1cm, DCAz<0.1cm);z_{TPC};(m/z)^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
    fM2vsZ[iB][1] = new TH2F("fM2vsZ_1","(m/z)^{2}_{TOF} vs z_{TPC} over all p/z (dEdx>80, DCAxy<0.1cm, DCAz<3.2cm);z_{TPC};(m/z)^{2}_{TOF} (GeV^{2}/c^{4})",hbins[0],-4,4,hbins[1],0,10);
        
    Char_t name_h2DCAap[18][200];
    Char_t title_h2DCAap[18][200];
    
    for(Int_t iS=0;iS<nSpec;iS++) {
      snprintf(name_h2DCAap[iS],200,"h2DCAap_%s",name[iS]);
      snprintf(title_h2DCAap[iS],200,"h2DCA_%s in for p/z<1.5GeV;DCA_{xy} (cm);DCA_{z} (cm)",name[iS]);
      if(iS==5 || iS==7 || iS==5+9 || iS==7+9) h2DCAap[iB][iS] = new TH2F(name_h2DCAap[iS],title_h2DCAap[iS],1750,-3.5,3.5,1750,-3.5,3.5);//1750,-3.5,3.5,1750,-3.5,3.5
      else h2DCAap[iB][iS] = new TH2F(name_h2DCAap[iS],title_h2DCAap[iS],1750,-3.5,3.5,1750,-3.5,3.5);//1750,-3.5,3.5,1750,-3.5,3.5
    }
        
    Char_t name_hDCAxy[18][nbin][200];
    Char_t title_hDCAxy[18][nbin][200];
    Char_t name_hDCAz[18][nbin][200];
    Char_t title_hDCAz[18][nbin][200];
    
    //Char_t name_h2DCA[18][nbin][200];
    //Char_t title_h2DCA[18][nbin][200];
    
    for(Int_t iS=0;iS<nSpec;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	snprintf(name_hDCAxy[iS][j],200,"hDCAxy_%s_%s",name[iS],name_nbin[j]);
	snprintf(title_hDCAxy[iS][j],200,"hDCAxy_%s_%s in DCAzCut;DCA_{xy} (cm)",name[iS],name_nbin[j]);
	if(iS==5 || iS==7 || iS==5+9 || iS==7+9) hDCAxy[iB][iS][j] = new TH1D(name_hDCAxy[iS][j],title_hDCAxy[iS][j],875,-3.5,3.5);
	else hDCAxy[iB][iS][j] = new TH1D(name_hDCAxy[iS][j],title_hDCAxy[iS][j],1,-3.5,3.5);

	snprintf(name_hDCAz[iS][j],200,"hDCAz_%s_%s",name[iS],name_nbin[j]);
	snprintf(title_hDCAz[iS][j],200,"hDCAz_%s_%s in DCAxyCut;DCA_{z} (cm)",name[iS],name_nbin[j]);
	if(iS==5 || iS==7 || iS==5+9 || iS==7+9) hDCAz[iB][iS][j] = new TH1D(name_hDCAz[iS][j],title_hDCAz[iS][j],875,-3.5,3.5);
	else hDCAz[iB][iS][j] = new TH1D(name_hDCAz[iS][j],title_hDCAz[iS][j],1,-3.5,3.5);
      
	//snprintf(name_h2DCA[iS][j],200,"h2DCA_%s_%s",name[iS],name_nbin[j]);
	//snprintf(title_h2DCA[iS][j],200,"h2DCA_%s_%s;DCA_{xy} (cm);DCA_{z} (cm)",name[iS],name_nbin[j]);
	//if(iS==2 || iS==5 || iS==7 || iS==2+9 || iS==5+9 || iS==7+9) h2DCA[iB][iS][j] = new TH2F(name_h2DCA[iS][j],title_h2DCA[iS][j],1,-4,4,1,-4,4);//,160,-4,4,160,-4,4);
	//else h2DCA[iB][iS][j] = new TH2F(name_h2DCA[iS][j],title_h2DCA[iS][j],1,-4,4,1,-4,4);
      }
    }
      
    Char_t name_hM2CutDCAxy[18][nbin][200];
    Char_t title_hM2CutDCAxy[18][nbin][200];
    for(Int_t iS=0;iS<nSpec;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	snprintf(name_hM2CutDCAxy[iS][j],200,"hM2_CutDCAxy_%s_%s",name[iS],name_nbin[j]);
	snprintf(title_hM2CutDCAxy[iS][j],200,"m^{2}/z^{2} Tof distribution of %s and in %s;m^{2}/z^{2}_{TOF} (GeV^{2}/c^{4})",name[iS],name_nbin[j]);
      }
    }

    const Int_t BinM2pT[nPart]={1,1,1,1000,500,500,1,400,1};//1,1,600,250,500,500,1000,400,600//1,1,1,250,500,500,1,400,1//1,1,1,1000,500,500,1,400,1
    const Double_t RangeM2min[nPart]={0.0,0.0,-0.1,0.0,0.0,0.0,0.0,0.0,0.0};
    const Double_t RangeM2max[nPart]={1.0,1.0,0.5,2.0,4.0,6.0,12.0,4.0,6.0};

    for(Int_t iS=0;iS<nPart;iS++) {
      for(Int_t j=0;j<nbin;j++) {
	hM2CutDCAxy[iB][iS][j] = new TH1D(name_hM2CutDCAxy[iS][j],title_hM2CutDCAxy[iS][j],BinM2pT[iS],RangeM2min[iS],RangeM2max[iS]);
	hM2CutDCAxy[iB][iS+nPart][j] = new TH1D(name_hM2CutDCAxy[iS+nPart][j],title_hM2CutDCAxy[iS+nPart][j],BinM2pT[iS],RangeM2min[iS],RangeM2max[iS]);
      }
    }

    Char_t name_fPmeanVsBetaGamma[18][200];
    Char_t title_fPmeanVsBetaGamma[18][200];
    
    if(iMtof==2) {hbins[0]=1; hbins[1]=1;}//if(iMtof==2) {hbins[0]=200; hbins[1]=200;}
    else {hbins[0]=1; hbins[1]=1;}
    for(Int_t iS=0;iS<nSpec;iS++) {
      snprintf(name_fPmeanVsBetaGamma[iS],200,"fPmeanVsPvtx_%s",name[iS]);
      snprintf(title_fPmeanVsBetaGamma[iS],200,"<p>/p_{vtx} vs #beta#gamma of %s;p_{vtx}/m_{%s};<p>_{%s}/p_{vtx}",name[iS],name[iS],name[iS]);
      fPmeanVsBetaGamma[iB][iS]=new TH2F(name_fPmeanVsBetaGamma[iS],title_fPmeanVsBetaGamma[iS],hbins[0],0,10,hbins[1],0.8,1.2);
    }	
    
    Char_t name_prPmeanVsBetaGamma[18][200];
    Char_t title_prPmeanVsBetaGamma[18][200];
    
    if(iMtof==2) {hbins[0]=1; hbins[1]=1;}//if(iMtof==2) {hbins[0]=200; hbins[1]=200;}
    else {hbins[0]=1; hbins[1]=1;}
    for(Int_t iS=0;iS<nSpec;iS++) {
      snprintf(name_prPmeanVsBetaGamma[iS],200,"prPmeanVsPvtx_%s",name[iS]);
      snprintf(title_prPmeanVsBetaGamma[iS],200,"<p>/p_{vtx} vs #beta#gamma of %s;p_{vtx}/m_{%s};<p>_{%s}/p_{vtx}",name[iS],name[iS],name[iS]);
      prPmeanVsBetaGamma[iB][iS]=new TProfile(name_prPmeanVsBetaGamma[iS],title_prPmeanVsBetaGamma[iS],hbins[0],0,10,0.8,1.2,"");
    }	
    
    SetPvtxCorrections();

    prPvtxTrueVsReco[iB][0]=new TProfile("prPvtxTrueVsReco_d","p_{true} vs p_{reco} of d and dbar;p_{reco} (GeV/c); p_{true}/p_{reco} (d)",1,0,10);//,100,0,10
    prPvtxTrueVsReco[iB][1]=new TProfile("prPvtxTrueVsReco_t","p_{true} vs p_{reco} of t and tbar;p_{reco} (GeV/c);p_{true}/p_{reco} (t)",1,0,10);//,100,0,10
    prPvtxTrueVsReco[iB][2]=new TProfile("prPvtxTrueVsReco_He3","p_{true} vs p_{reco} of He3 and He3bar;p_{reco} (GeV/c);p_{true}/p_{reco} (He3)",1,0,10);//,100,0,10
    prPvtxTrueVsReco[iB][3]=new TProfile("prPvtxTrueVsReco_He4","p_{true} vs p_{reco} of He4 and He4bar;p_{reco} (GeV/c);p_{true}/p_{reco} (He4)",1,0,10);//,100,0,10

    SetPmeanCorrections();
       
    Char_t nameTemp[14][200];
    snprintf(nameTemp[0],200,"#pi^{+}");
    snprintf(nameTemp[1],200,"K^{+}");
    snprintf(nameTemp[2],200,"p");
    snprintf(nameTemp[3],200,"d");
    snprintf(nameTemp[4],200,"t");
    snprintf(nameTemp[5],200,"He3");
    snprintf(nameTemp[6],200,"He4");
    snprintf(nameTemp[7],200,"#pi^{-}");
    snprintf(nameTemp[8],200,"K^{-}");
    snprintf(nameTemp[9],200,"#bar{p}");
    snprintf(nameTemp[10],200,"#bar{d}");
    snprintf(nameTemp[11],200,"#bar{t}");
    snprintf(nameTemp[12],200,"#bar{He3}");
    snprintf(nameTemp[13],200,"#bar{He4}");
    Char_t name_prPmeanVsBGcorr[14][200];
    Char_t title_prPmeanVsBGcorr[14][200];
   
    hbins[0]=200;
    for(Int_t iS=0;iS<14;iS++) {
      snprintf(name_prPmeanVsBGcorr[iS],200,"prPmeanVsBGcorr_%s",nameTemp[iS]);
      snprintf(title_prPmeanVsBGcorr[iS],200,"<p>/p_{vtx} vs #beta#gamma of %s as parameterized in input TF1;p_{vtx}/m_{%s};<p>_{%s}/p_{vtx}",nameTemp[iS],nameTemp[iS],nameTemp[iS]);
      prPmeanVsBGcorr[iB][iS]=new TProfile(name_prPmeanVsBGcorr[iS],title_prPmeanVsBGcorr[iS],hbins[0],0,20,0.8,1.2,"");
    }	

    fList[iB]->Add(htemp[iB]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(htriggerbits[iB][i]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(hCentrality[iB][i]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(hZvertex[iB][i]);
    fList[iB]->Add(hEta[iB]);
    fList[iB]->Add(hPhi[iB]);
    //fList[iB]->Add(fEtaPhi[iB]);
    fList[iB]->Add(hNTpcCluster[iB]);
    fList[iB]->Add(hNTrdSlices[iB]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fdEdxVSp[iB][i]);
    for(Int_t i=0;i<nPart;i++) fList[iB]->Add(hDeDxExp[iB][i]);
    for(Int_t i=0;i<nSpec;i++) fList[iB]->Add(fNsigmaTpc[iB][i]);
    for(Int_t i=0;i<nPart;i++) {
      if(kSignalCheck!=1) 
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      //fList[iB]->Add(fNsigmaTpc_kTOF[iB][i]);
      //fList[iB]->Add(fNsigmaTpc_kTOF[iB][i+nPart]);
    }
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fBetaTofVSp[iB][i]);
    for(Int_t i=0;i<nPart;i++) fList[iB]->Add(hBetaExp[iB][i]);
    for(Int_t i=0;i<nPart;i++) fList[iB]->Add(fNsigmaTof[iB][i]);
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fM2vsP_NoTpcCut[iB][0][i]);
    for(Int_t i=0;i<nPart;i++) {
      //if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      fList[iB]->Add(fM2vsP[iB][0][i]);
      fList[iB]->Add(fM2vsP[iB][0][i+nPart]);
    }
  
    for(Int_t i=0;i<2;i++){
      //fList[iB]->Add(fPvtxTrueVsReco[i]);
      //fList[iB]->Add(prPvtxTrueVsReco[iB][i]);
    }
    if(iMtof==2) {
      for(Int_t i=0;i<nPart;i++){
	if(i<2) continue;//e,mu excluded
	//fList[iB]->Add(fPmeanVsBetaGamma[iB][i]);
	//fList[iB]->Add(prPmeanVsBetaGamma[iB][i]);
	//fList[iB]->Add(fPmeanVsBetaGamma[iB][i+nPart]);
	//fList[iB]->Add(prPmeanVsBetaGamma[iB][i+nPart]);
      }
    }
    if(iMtof>2) {
      //for(Int_t i=0;i<14;i++)fList[iB]->Add(fPmeanVsBGcorr[i]);
      //for(Int_t i=0;i<14;i++)fList[iB]->Add(prPmeanVsBGcorr[iB][i]);
    }
    for(Int_t i=0;i<nPart;i++) {
      if(i<5 || i==6 || i==8) continue;//e,mu,pi,K,p,t,he4 excluded//i<5 || i==6 || i==8
      fList[iB]->Add(h2DCAap[iB][i]);
      fList[iB]->Add(h2DCAap[iB][i+nPart]);
    }
    /*
      for(Int_t i=0;i<nPart;i++) {
      if(i<5 || i==6 || i==8) continue;//e,mu,pi,K,p,t,he4 excluded//i<5 || i==6 || i==8
      for(Int_t j=0;j<nbin;j++){
      fList[iB]->Add(h2DCA[iB][i][j]);
      fList[iB]->Add(h2DCA[iB][i+nPart][j]);
      }
      }
    */
    for(Int_t i=0;i<2;i++) fList[iB]->Add(fM2vsZ[iB][i]);
    for(Int_t i=0;i<nPart;i++){
      if(kSignalCheck!=1) 
	if(i<3 || i==6 || i==8) continue;//e,mu,pi,t,he4 excluded
      for(Int_t j=0;j<nbin;j++){
	fList[iB]->Add(hDCAxy[iB][i][j]);
	fList[iB]->Add(hDCAz[iB][i][j]);
	fList[iB]->Add(hM2CutDCAxy[iB][i][j]);
	fList[iB]->Add(hDCAxy[iB][i+nPart][j]);
	fList[iB]->Add(hDCAz[iB][i+nPart][j]);
	fList[iB]->Add(hM2CutDCAxy[iB][i+nPart][j]);
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

    //TRIGGER SELECTION
    Int_t iTrigger=-2;

    if(inputHandler->IsEventSelected() & AliVEvent::kMB) iTrigger = 0;
    if(inputHandler->IsEventSelected() & AliVEvent::kCentral) iTrigger = 16;
    if(inputHandler->IsEventSelected() & AliVEvent::kSemiCentral) iTrigger = 17;
    //if((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & AliVEvent::kAny) iTrigger = 35;
    
    if(iTriggerSel!=-99) {//if a dedicated trigger is required
      if(iTrigger!=iTriggerSel) return;
    }
    
    for(Int_t i=0;i<32;i++) {
      Int_t bit=(1<<i);
      if(inputHandler->IsEventSelected() & bit) htriggerbits[iBconf][0]->Fill(i);
    }
    if(inputHandler->IsEventSelected() & AliVEvent::kAny) htriggerbits[iBconf][0]->Fill(35);
    if(inputHandler->IsEventSelected() & AliVEvent::kAnyINT) htriggerbits[iBconf][0]->Fill(36);
    
    htriggerbits[iBconf][1]->Fill(iTrigger);
    
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

      //-------------------------------------start TRACK CUTS (I): for (II) see below--------
      if ((track->Pt() < 0.2) || (eta<EtaLimit[0]) || (eta>EtaLimit[1]) || !trkFlag || !isMinTpcCluster)
	continue;
      
      //For the Tpc purity cut
      Double_t dedx = track->GetTPCsignal();
      if(dedx<10) continue;

      Int_t nTrdSlices = track->GetNumberOfTRDslices();
      if(nTrdSlices<2 && iTrdCut==1) continue; 
      if(nTrdSlices>0 && iTrdCut==2) continue;
      
      //-------------------------------------end TRACK CUTS (I)----------------------------------

      //-------------------------------------Track info--------------------------------------
      Double_t phi= track->Phi();
      Double_t charge = (Double_t)track->Charge();
      Double_t p = track->P();
      Double_t pt = track->Pt();
      Double_t tof  = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
      Double_t pTPC = track->GetTPCmomentum();
      Double_t beta = 0.0;
      Double_t M2 = 999.9;
      Double_t Z2 = 999.9;
      
       //Vertex determination
      Double_t b[2] = {-99., -99.};
      Double_t bCov[3] = {-99., -99., -99.};
      if (!track->PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 100., b, bCov))
	continue;
      
      Double_t DCAxy = b[0];
      Double_t DCAz = b[1];
     
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
      if(kPvtxCorr) this->MomVertexCorrection(p,pC,eta,FlagPid);
      
      this->FillDCAdist(DCAxy,DCAz,charge,FlagPid,stdFlagPid,pC);
      
      //-------------------------------------start TRACK CUTS (II)-------------------------------------
      //Cut on the DCAxy
      Bool_t isDCAxyCut=kFALSE;
      if(TMath::Abs(DCAxy)<DCAxyCut) isDCAxyCut=kTRUE;
      
      //Cut on the DCAz
      Bool_t isDCAzCut=kFALSE;
      if(TMath::Abs(DCAz)<DCAzCut) isDCAzCut=kTRUE;
      
      if (!isDCAxyCut || !isDCAzCut)
	continue;
          
      //-------------------------------------end TRACK CUTS (II)----------------------------------
      
      hEta[iBconf]->Fill(eta);
      hPhi[iBconf]->Fill(phi);
      fEtaPhi[iBconf]->Fill(eta,phi);
      hNTpcCluster[iBconf]->Fill(nTpcCluster);
      hNTrdSlices[iBconf]->Fill(nTrdSlices);
      
      //More TPC info:
      for(Int_t iS=0;iS<9;iS++){
	expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
	hDeDxExp[iBconf][iS]->Fill(pTPC,expdedx[iS]);
	nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
	//fNsigmaTpc[iBconf][iS]->Fill(pTPC,nsigmaTPC[iS]);
	if(charge>0) {//positive particle
	  fNsigmaTpc[iBconf][iS]->Fill(pTPC,nsigmaTPC[iS]);
	  if(kTOF) fNsigmaTpc_kTOF[iBconf][iS]->Fill(p,nsigmaTPC[iS]);
	}
	else {//negative particle
	  fNsigmaTpc[iBconf][iS+nPart]->Fill(pTPC,nsigmaTPC[iS]);
	  if(kTOF) fNsigmaTpc_kTOF[iBconf][iS+nPart]->Fill(p,nsigmaTPC[iS]);
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
	  }
	  else {
	    hBetaExp[iBconf][iS+nPart]->Fill(p,exptimes[0]/exptimes[iS]);
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

	if(iMtof==2) this->GetPmeanVsBetaGamma(exptimes,pC,FlagPid,FlagPidTof,charge);
	
	//-----------------------------M2 as a function of expected times---------------------------------
	if(iMtof==2) this->GetMassFromExpTimes(beta,exptimes,Mass2);
	 
	//-----------------------------M2 as a function of mean momentum calculated from expected time and extrapolated to the (anti)nuclei---------------------------------
	if(iMtof>2) this->GetMassFromMeanMom(beta,exptimes,pC,eta,charge,Mass2,FlagPid,FlagPidTof);

	//-------------------------------Squared Mass TH2 distributions-----------------------
	if(charge>0) {
	  //without TPC
	  fM2vsP_NoTpcCut[iBconf][0][0]->Fill(M2,p);
	  //with TPC
	  for(Int_t iS=0;iS<9;iS++) {
	    M2=999.9;
	    M2=Mass2[iS];
	    //-----------------
	    if(FlagPid & stdFlagPid[iS]) {
	      fM2vsP[iBconf][0][iS]->Fill(M2,pC[iS]);
	    }
	  }
	}
	else {//charge<0
	  //without TPC
	  fM2vsP_NoTpcCut[iBconf][0][1]->Fill(M2,p);
	  //with TPC
	  for(Int_t iS=0;iS<9;iS++) {
	    M2=999.9;
	    M2=Mass2[iS];
	    //-----------------
	    if(FlagPid & stdFlagPid[iS]) {
	      fM2vsP[iBconf][0][iS+nPart]->Fill(M2,pC[iS]);
	    }
	  }
	}
	
	//------------------------------start Squared Mass TH1 distributions-------------------------
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
		  hM2CutDCAxy[iBconf][iS][j]->Fill(M2);
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
		  hM2CutDCAxy[iBconf][iS+nPart][j]->Fill(M2);
		  break;
		}
	      }//end loop on the p bins (j)
	    }
	  }//end loop on the particle species (iS)
	}
		
	//-------------------------------------------------M2/Z2 vs Z-------------------------
	

	// Double_t binCutPt[10] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
	Double_t Z=999.9;
	if(Z2>0) Z=TMath::Sqrt(Z2);
	
	if(dedx>80) {
	  if(TMath::Abs(DCAz)<0.1) {
	    fM2vsZ[iBconf][0]->Fill(charge*Z,M2);
	  }
	  fM2vsZ[iBconf][1]->Fill(charge*Z,M2);
	}//end dedx>80 requirement

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
      else if(iS==6) {
	if(kPvtxCorr==1) pC[iS]=pC[iS]*fPvtxTrueVsReco[1]->Eval(pC[iS],TMath::Abs(eta));//for (bar)t
	prPvtxTrueVsReco[iBconf][1]->Fill(p,pC[iS]/p);
      }
      else if(iS==7) {
	if(kPvtxCorr==1) pC[iS]=pC[iS]*fPvtxTrueVsReco[2]->Eval(pC[iS],TMath::Abs(eta));//for (bar)He3
	prPvtxTrueVsReco[iBconf][2]->Fill(p,pC[iS]/p);
      }
      else if(iS==8) {
	if(kPvtxCorr==1) pC[iS]=pC[iS]*fPvtxTrueVsReco[3]->Eval(pC[iS],TMath::Abs(eta));//for (bar)He3
	prPvtxTrueVsReco[iBconf][3]->Fill(p,pC[iS]/p);
      }
    }
  }
  
  return;
  
}
//__________________________________________________________________________________________________
void AliAnalysisNucleiMass::FillDCAdist(Double_t DCAxy, Double_t DCAz, Double_t charge, Int_t FlagPid, Int_t stdFlagPid[9], Double_t *pC){

  Double_t binP[nbin+1];
  for(Int_t i=0;i<nbin+1;i++) {
    binP[i]=0.4+i*0.1;
  }

  if(charge>0) {
    for(Int_t iS=0;iS<9;iS++) {
      if(FlagPid & stdFlagPid[iS]) {
	if(pC[iS]<1.5) {
	  h2DCAap[iBconf][iS]->Fill(DCAxy,DCAz);
	  h2DCAap[iBconf][iS]->Fill(-DCAxy,-DCAz);
	}
	for(Int_t j=0;j<nbin;j++) {
	  if(pC[iS]>binP[j] && pC[iS]<binP[j+1]) {
	    if(TMath::Abs(DCAz)<DCAzCut) {
	      hDCAxy[iBconf][iS][j]->Fill(DCAxy);
	      hDCAxy[iBconf][iS][j]->Fill(-DCAxy);
	    }
	    if(TMath::Abs(DCAxy)<DCAxyCut) {
	      hDCAz[iBconf][iS][j]->Fill(DCAz);
	      hDCAz[iBconf][iS][j]->Fill(-DCAz);
	    }
	    //h2DCA[iBconf][iS][j]->Fill(DCAxy,DCAz);
	    //h2DCA[iBconf][iS][j]->Fill(-DCAxy,-DCAz);
	    break;
	  }
	}//end loop on the p bins (j)
      }
    }//end loop on the particle species (iS)
  }
  else {//charge<0
    for(Int_t iS=0;iS<9;iS++) {
      if(FlagPid & stdFlagPid[iS]) {
	if(pC[iS]<1.5) {
	  h2DCAap[iBconf][iS+nPart]->Fill(DCAxy,DCAz);
	  h2DCAap[iBconf][iS+nPart]->Fill(-DCAxy,-DCAz);
	}
	for(Int_t j=0;j<nbin;j++) {
	  if(pC[iS]>binP[j] && pC[iS]<binP[j+1]) {
	    if(TMath::Abs(DCAz)<DCAzCut) {
	      hDCAxy[iBconf][iS+nPart][j]->Fill(DCAxy);
	      hDCAxy[iBconf][iS+nPart][j]->Fill(-DCAxy);
	    }
	    if(TMath::Abs(DCAxy)<DCAxyCut) {
	      hDCAz[iBconf][iS+nPart][j]->Fill(DCAz);
	      hDCAz[iBconf][iS+nPart][j]->Fill(-DCAz);
	    }
	    //h2DCA[iBconf][iS+nPart][j]->Fill(DCAxy,DCAz);
	    //h2DCA[iBconf][iS+nPart][j]->Fill(-DCAxy,-DCAz);
	    break;
	  }
	}//end loop on the p bins (j)
      }
    }//end loop on the particle species (iS)
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
void AliAnalysisNucleiMass::GetPmeanVsBetaGamma(Double_t *IntTimes, Double_t *pVtx, Int_t FlagPid, Int_t FlagPidTof, Double_t charge) {
 
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
void AliAnalysisNucleiMass::GetMassFromMeanMom(Double_t beta, Double_t *IntTimes, Double_t *pVtx, Double_t eta, Double_t charge, Double_t *Mass2, Int_t FlagPid, Int_t FlagPidTof) {//Double_t *Mass2, Int_t iCorr
 
  // m = p_exp/beta/gamma where p_exp = mPDG*beta_exp*gamma_exp; beta_exp = L/t_exp/c = t_e/t_exp ; beta=L/tof/c = t_e/tof
  // In this way m_tof = mPDG only if tof=t_exp
  
  Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};
  
  Double_t beta2Exp[9];
  Double_t p2Exp[9];
  
  Double_t pExp[9];
  
  Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,#mu,#pi,K,p,d,t,3He,4He
  
  for(Int_t iS=0;iS<9;iS++) {
    if(iS>1) {
      p2Exp[iS]=pVtx[iS]*fPmeanVsBGcorr[iS-2]->Eval(pVtx[iS]/massOverZ[iS],TMath::Abs(eta));
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
    
    if(p2Exp[iS]<0) {
      Mass2[iS]=999.9;
      continue;
    }
    pExp[iS]=TMath::Sqrt(p2Exp[iS]);
    
    //------------
    Mass2[iS]=p2Exp[iS]*(1-beta*beta)/(beta*beta);
    
    //-----------
    
    if(iS>1) {
      if((FlagPid & stdFlagPid[iS]) && (FlagPidTof & stdFlagPid[iS])) {
	if(charge>0) {
	  prPmeanVsBGcorr[iBconf][iS-2]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
	}
	else if(charge<0) {
	  prPmeanVsBGcorr[iBconf][iS-2+7]->Fill(pVtx[iS]/massOverZ[iS],pExp[iS]/pVtx[iS]);
	}
      }
    }
  }//end loop on the particle species
    
  return;
  
}
//________________________________________________________________________________________
void AliAnalysisNucleiMass::SetPvtxCorrections(){
  //for (bar)d
  fPvtxTrueVsReco[0]=new TF2("fcorr_d","([0]*TMath::Power(x,[1])+[2])*(TMath::Power((TMath::Exp([3]*x)+[4]),[5]*TMath::Power(y,[6])));p_{reco};|#eta|;p_{true}/p_{reco}",0.0001,100,0,1);//for (bar)d
  fPvtxTrueVsReco[0]->SetParameter(0,0.031263);
  fPvtxTrueVsReco[0]->SetParameter(1,-3.276770);
  fPvtxTrueVsReco[0]->SetParameter(2,1.000113);
  fPvtxTrueVsReco[0]->SetParameter(3,-5.195875);
  fPvtxTrueVsReco[0]->SetParameter(4,1.000674);
  fPvtxTrueVsReco[0]->SetParameter(5,2.870503);
  fPvtxTrueVsReco[0]->SetParameter(6,3.777729);
  
  //for (bar)t
  fPvtxTrueVsReco[1]=new TF2("fcorr_t","([0]*TMath::Power(x,[1])+[2])+[3]*y;p_{reco};|#eta|;p_{true}/p_{reco}",0.0001,100,0,1);//for (bar)He3
  fPvtxTrueVsReco[1]->SetParameter(0,8.79761e-02);
  fPvtxTrueVsReco[1]->SetParameter(1,-3.23189e+00);
  fPvtxTrueVsReco[1]->SetParameter(2,9.99578e-01);
  fPvtxTrueVsReco[1]->SetParameter(3,0.0);
  
  //for (bar)He3
  fPvtxTrueVsReco[2]=new TF2("fcorr_He","([0]*TMath::Power(x,[1])+[2])*(TMath::Power((TMath::Exp([3]*x)+[4]),[5]*TMath::Power(y,[6])));p_{reco};|#eta|;p_{true}/p_{reco}",0.0001,100,0,1);//for (bar)He3
  fPvtxTrueVsReco[2]->SetParameter(0,0.037986);
  fPvtxTrueVsReco[2]->SetParameter(1,-2.707620);
  fPvtxTrueVsReco[2]->SetParameter(2,1.000742);
  fPvtxTrueVsReco[2]->SetParameter(3,-4.934743);
  fPvtxTrueVsReco[2]->SetParameter(4,1.001640);
  fPvtxTrueVsReco[2]->SetParameter(5,2.744372);
  fPvtxTrueVsReco[2]->SetParameter(6,3.528561);
  
  //for (bar)He4
  fPvtxTrueVsReco[3]=new TF2("fcorr_He4","([0]*TMath::Power(x,[1])+[2])+[3]*y;p_{reco};|#eta|;p_{true}/p_{reco}",0.0001,100,0,1);//for (bar)He3
  fPvtxTrueVsReco[3]->SetParameter(0,7.08785e-02);
  fPvtxTrueVsReco[3]->SetParameter(1,-2.87201e+00);
  fPvtxTrueVsReco[3]->SetParameter(2,1.00070e+00);
  fPvtxTrueVsReco[3]->SetParameter(3,0.0);
  
  for(Int_t i=0;i<4;i++) {
    fPvtxTrueVsReco[i]->SetNpx(fPvtxTrueVsReco[i]->GetNpx()*10.0);
  }
}
//________________________________________________________________________________________
void AliAnalysisNucleiMass::SetPmeanCorrections(){
  
  Char_t nameTemp[14][200];
  snprintf(nameTemp[0],200,"#pi^{+}");
  snprintf(nameTemp[1],200,"K^{+}");
  snprintf(nameTemp[2],200,"p");
  snprintf(nameTemp[3],200,"d");
  snprintf(nameTemp[4],200,"t");
  snprintf(nameTemp[5],200,"He3");
  snprintf(nameTemp[6],200,"He4");
  snprintf(nameTemp[7],200,"#pi^{-}");
  snprintf(nameTemp[8],200,"K^{-}");
  snprintf(nameTemp[9],200,"#bar{p}");
  snprintf(nameTemp[10],200,"#bar{d}");
  snprintf(nameTemp[11],200,"#bar{t}");
  snprintf(nameTemp[12],200,"#bar{He3}");
  snprintf(nameTemp[13],200,"#bar{He4}");
  
  Char_t name_fPmeanVsBGcorr[14][200];
  for(Int_t i=0;i<14;i++) {
    snprintf(name_fPmeanVsBGcorr[i],200,"fPmeanVsBGcorr_%s",nameTemp[i]);
  }

  //Pions
  fPmeanVsBGcorr[0]=new TF2(name_fPmeanVsBGcorr[0],"(x>[5])*([2]-[0]*TMath::Power(x,[1]))*([3]+[4]*y*y)+(x<=[5])*[6]",0.0001,100,0,0.8);
  fPmeanVsBGcorr[0]->SetParameter(0,-0.179607);
  fPmeanVsBGcorr[0]->SetParameter(1,-0.384809);
  fPmeanVsBGcorr[0]->SetParameter(2,0.885534);
  fPmeanVsBGcorr[0]->SetParameter(3,0.992710);
  fPmeanVsBGcorr[0]->SetParameter(4,0.011390);
  fPmeanVsBGcorr[0]->SetParameter(5,3.231000);
  fPmeanVsBGcorr[0]->SetParameter(6,0.999900);
 
  //Kaons
  fPmeanVsBGcorr[1]=new TF2(name_fPmeanVsBGcorr[1],"(x>[8])*([2]-[0]*TMath::Power(x,[1]))*TMath::Power([3]+[4]*TMath::Exp([5]*x),[6]+[7]*y*y)+(x<=[8])*[9]",0.0001,20,0,0.8);
  fPmeanVsBGcorr[1]->SetParameter(0,0.033500);
  fPmeanVsBGcorr[1]->SetParameter(1,-2.461673);
  fPmeanVsBGcorr[1]->SetParameter(2,0.996501);
  fPmeanVsBGcorr[1]->SetParameter(3,1.000000);
  fPmeanVsBGcorr[1]->SetParameter(4,0.089715);
  fPmeanVsBGcorr[1]->SetParameter(5,-2.473531);
  fPmeanVsBGcorr[1]->SetParameter(6,1.000000);
  fPmeanVsBGcorr[1]->SetParameter(7,-1.562500);
  fPmeanVsBGcorr[1]->SetParameter(8,0.253000);
  fPmeanVsBGcorr[1]->SetParameter(9,0.009387);
  
  //Protons
  fPmeanVsBGcorr[2]=new TF2(name_fPmeanVsBGcorr[2],"(x>[8])*([2]-[0]*TMath::Power(x,[1]))*TMath::Power([3]+[4]*TMath::Exp([5]*x),[6]+[7]*y*y)+(x<=[8])*[9]",0.0001,20,0,0.8);
  fPmeanVsBGcorr[2]->SetParameter(0,0.015081);
  fPmeanVsBGcorr[2]->SetParameter(1,-2.927557);
  fPmeanVsBGcorr[2]->SetParameter(2,0.997904);
  fPmeanVsBGcorr[2]->SetParameter(3,1.000000);
  fPmeanVsBGcorr[2]->SetParameter(4,0.102697);
  fPmeanVsBGcorr[2]->SetParameter(5,-3.399528);
  fPmeanVsBGcorr[2]->SetParameter(6,1.000000);
  fPmeanVsBGcorr[2]->SetParameter(7,-1.562500);
  fPmeanVsBGcorr[2]->SetParameter(8,0.239000);
  fPmeanVsBGcorr[2]->SetParameter(9,0.002054);

  //Deuterons
  fPmeanVsBGcorr[3]=new TF2(name_fPmeanVsBGcorr[3],"(x>[8])*([2]-[0]*TMath::Power(x,[1]))*TMath::Power([3]+[4]*TMath::Exp([5]*x),[6]+[7]*y*y)+(x<=[8])*[9]",0.0001,20,0,0.8);
  fPmeanVsBGcorr[3]->SetParameter(0,0.008672);
  fPmeanVsBGcorr[3]->SetParameter(1,-2.712343);
  fPmeanVsBGcorr[3]->SetParameter(2,0.997639);
  fPmeanVsBGcorr[3]->SetParameter(3,1.000000);
  fPmeanVsBGcorr[3]->SetParameter(4,0.039627);
  fPmeanVsBGcorr[3]->SetParameter(5,-2.768122);
  fPmeanVsBGcorr[3]->SetParameter(6,1.000000);
  fPmeanVsBGcorr[3]->SetParameter(7,-1.562500);
  fPmeanVsBGcorr[3]->SetParameter(8,0.174000);
  fPmeanVsBGcorr[3]->SetParameter(9,0.002189);

  //Triton
  fPmeanVsBGcorr[4]=new TF2(name_fPmeanVsBGcorr[4],"(x>[4])*([2]-[0]*TMath::Power(x,[1])+[3]*y)+(x<=[4])*[5]",0.0001,20,0,0.8);
  fPmeanVsBGcorr[4]->SetParameter(0,6.79641e-03);
  fPmeanVsBGcorr[4]->SetParameter(1,-1.92801e+00);
  fPmeanVsBGcorr[4]->SetParameter(2,1.000000);
  fPmeanVsBGcorr[4]->SetParameter(3,0.0);
  fPmeanVsBGcorr[4]->SetParameter(4,0.076);
  fPmeanVsBGcorr[4]->SetParameter(5,2.25779e-02);

  //Helium-3
  fPmeanVsBGcorr[5]=new TF2(name_fPmeanVsBGcorr[5],"(x>[8])*([2]-[0]*TMath::Power(x,[1]))*TMath::Power([3]+[4]*TMath::Exp([5]*x),[6]+[7]*y*y)+(x<=[8])*[9]",0.0001,20,0,0.8);
  fPmeanVsBGcorr[5]->SetParameter(0,0.024339);
  fPmeanVsBGcorr[5]->SetParameter(1,-2.922613);
  fPmeanVsBGcorr[5]->SetParameter(2,0.993761);
  fPmeanVsBGcorr[5]->SetParameter(3,1.000000);
  fPmeanVsBGcorr[5]->SetParameter(4,1.087549);
  fPmeanVsBGcorr[5]->SetParameter(5,-6.216154);
  fPmeanVsBGcorr[5]->SetParameter(6,1.000000);
  fPmeanVsBGcorr[5]->SetParameter(7,-1.562500);
  fPmeanVsBGcorr[5]->SetParameter(8,0.282000);
  fPmeanVsBGcorr[5]->SetParameter(9,0.009711);

  //Helium-4
  fPmeanVsBGcorr[6]=new TF2(name_fPmeanVsBGcorr[6],"(x>[4])*([2]-[0]*TMath::Power(x,[1])+[3]*y)+(x<=[4])*[5]",0.0001,20,0,0.8);
  fPmeanVsBGcorr[6]->SetParameter(0,2.34185e-02);
  fPmeanVsBGcorr[6]->SetParameter(1,-2.31200e+00);
  fPmeanVsBGcorr[6]->SetParameter(2,1.000000);
  fPmeanVsBGcorr[6]->SetParameter(3,0.0);
  fPmeanVsBGcorr[6]->SetParameter(4,0.198);
  fPmeanVsBGcorr[6]->SetParameter(5,9.9226e-03);

  for(Int_t i=7;i<14;i++) {
    fPmeanVsBGcorr[i]=(TF2 *)fPmeanVsBGcorr[i-7]->Clone();
    fPmeanVsBGcorr[i]->SetName(name_fPmeanVsBGcorr[i]);
  }
    
  for(Int_t i=0;i<14;i++) {
    fPmeanVsBGcorr[i]->SetNpx(fPmeanVsBGcorr[i]->GetNpx()*100.0);
  }

  return;
  
}

