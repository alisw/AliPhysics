#if !defined(__CINT__) || defined(__MAKECINT__)
#include<Riostream.h>
#include<TROOT.h>
#include<TArrayI.h>
#include<TBranch.h>
#include<TCanvas.h>
#include<TClassTable.h>
#include<TClonesArray.h>
#include<TFile.h>
#include<TStyle.h>
#include<TH1.h>
#include<TH2.h>
#include<TLatex.h>
#include <TInterpreter.h>
#include <TGeoManager.h>
#include<TObject.h>
#include<TObjArray.h>
#include<TTree.h>
#include<TNtuple.h>
#include<TParticle.h>
#include "AliStack.h"
#include "AliRun.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliITS.h"
#include "AliITSgeomTGeo.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPoint.h"
#include "AliITSdigit.h"
#include "AliITSdigitSSD.h"
#include "AliITShit.h"
#include "AliITSmodule.h" 
#include "AliITSsegmentation.h"
#include "AliITSsegmentationSPD.h" 
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliHeader.h"
#endif

/*  $Id$    */

// macro to display the coordindates (local+global) and the energy deposit
// of ITS hits 


void PlotITSHits() {

 
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    printf("Setting a local default storage and run number 0\n");
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man->SetRun(0);
  }
  else {
    printf("Using deafult storage \n");
  }
  // retrives geometry 
  if(!gGeoManager){
    AliGeomManager::LoadGeometry("geometry.root");
  }
  AliGeomManager::ApplyAlignObjsFromCDB("ITS");

  Int_t totmod=AliITSgeomTGeo::GetNModules();
  Int_t modmin=AliITSgeomTGeo::GetModuleIndex(1,1,1);
  Int_t modmax=AliITSgeomTGeo::GetNModules()-1;

  Float_t xlim[6]={4.5,7.5,16.,26.,40.,45.};
  Float_t zlim[6]={15.,15.,22.,30.,45.,55.};

  TH1F* hlayer=new TH1F("hlayer","",6,-0.5,5.5);
  TH1F** hmod=new TH1F*[6];
  TH1F** hxl=new TH1F*[6];
  TH1F** hzl=new TH1F*[6];
  TH1F** hxg=new TH1F*[6];
  TH1F** hyg=new TH1F*[6];
  TH1F** hzg=new TH1F*[6];
  TH1F** hr=new TH1F*[6];
  TH1F** hphi=new TH1F*[6];
  TH2F** hzphi=new TH2F*[6];
  TH1F** hq=new TH1F*[6];

  Char_t name[10];
  for(Int_t iLay=0;iLay<6;iLay++){
    sprintf(name,"hmod%d",iLay+1);
    hmod[iLay]=new TH1F(name,"",totmod,modmin-0.5,modmax+0.5);
    hmod[iLay]->GetXaxis()->SetTitle("Module");
    hmod[iLay]->GetXaxis()->CenterTitle();
    sprintf(name,"hxloc%d",iLay+1);
    hxl[iLay]=new TH1F(name,"",100,-4.,4.);
    hxl[iLay]->GetXaxis()->SetTitle("Xloc (cm)");
    hxl[iLay]->GetXaxis()->CenterTitle();
    sprintf(name,"hzloc%d",iLay+1);
    hzl[iLay]=new TH1F(name,"",100,-4.,4.);
    hzl[iLay]->GetXaxis()->SetTitle("Zloc (cm)");
    hzl[iLay]->GetXaxis()->CenterTitle();
    sprintf(name,"hxgl%d",iLay+1);
    hxg[iLay]=new TH1F(name,"",100,-xlim[iLay],xlim[iLay]);
    hxg[iLay]->GetXaxis()->SetTitle("Xglob (cm)");
    hxg[iLay]->GetXaxis()->CenterTitle();
    sprintf(name,"hygl%d",iLay+1);
    hyg[iLay]=new TH1F(name,"",100,-xlim[iLay],xlim[iLay]);
    hyg[iLay]->GetXaxis()->SetTitle("Yglob (cm)");
    hyg[iLay]->GetXaxis()->CenterTitle();
    sprintf(name,"hzgl%d",iLay+1);
    hzg[iLay]=new TH1F(name,"",100,-zlim[iLay],zlim[iLay]);
    hzg[iLay]->GetXaxis()->SetTitle("Zglob (cm)");
    hzg[iLay]->GetXaxis()->CenterTitle();
    sprintf(name,"hr%d",iLay+1);
    hr[iLay]=new TH1F(name,"",100,0.,50.);
    hr[iLay]->GetXaxis()->SetTitle("r (cm)");
    hr[iLay]->GetXaxis()->CenterTitle();
    sprintf(name,"hphi%d",iLay+1);
    hphi[iLay]=new TH1F(name,"",100,-TMath::Pi(),TMath::Pi());    
    hphi[iLay]->GetXaxis()->SetTitle("#varphi (rad)");
    hphi[iLay]->GetXaxis()->CenterTitle();
    sprintf(name,"hq%d",iLay+1);
    hq[iLay]=new TH1F(name,"",100,0.,300.);    
    hq[iLay]->GetXaxis()->SetTitle("Charge (keV)");
    hq[iLay]->GetXaxis()->CenterTitle();
    sprintf(name,"hzphi%d",iLay+1);
    hzphi[iLay]=new TH2F(name,Form("Layer %d",iLay+1),50,-TMath::Pi(),TMath::Pi(),50,-zlim[iLay],zlim[iLay]);
    hzphi[iLay]->GetXaxis()->SetTitle("#varphi (rad)");
    hzphi[iLay]->GetYaxis()->SetTitle("Zglob (cm)");
    hzphi[iLay]->SetStats(0);
  }


  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0){
    cerr<<"Can not open session RL=NULL"<< endl;
    return;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval){
    cerr<<"LoadgAlice returned error"<<endl;
    return;
  }
  gAlice=rl->GetAliRun();
  
  retval = rl->LoadHeader();
  if (retval){
    cerr<<"LoadHeader returned error"<<endl;
    return;
  }
  
  retval = rl->LoadKinematics();
  if (retval){
    cerr<<"LoadKinematics returned error"<<endl;
    return;
  }
  
  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  if(!ITSloader){
    cerr<<"ITS loader not found"<<endl;
    return;
  }

  ITSloader->LoadHits("read");
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  ITS->SetTreeAddress();
  Int_t totev=rl->GetNumberOfEvents();

  for(Int_t iev=0; iev<totev; iev++){
    rl->GetEvent(iev);

    // HITS
    TTree *TH = ITSloader->TreeH();
    printf("Event %d  Tracks %d\n",iev,(Int_t)TH->GetEntries());
    
    // ITS
    Int_t nmodules;
    ITS->InitModules(-1,nmodules);
    ITS->FillModules(0,0,nmodules," "," ");  
   
    
    for (Int_t mod=modmin; mod<=modmax; mod++){
      Int_t lay,lad,det;
      AliITSgeomTGeo::GetModuleId(mod,lay,lad,det);  
      lay--;

      // Hits
      AliITSmodule *modu = ITS->GetModule(mod);
      TObjArray *arrHits = modu->GetHits();
      Int_t nhits = arrHits->GetEntriesFast();
      for (Int_t iHit=0;iHit<nhits;iHit++) {
	AliITShit *hit = (AliITShit*) arrHits->At(iHit);
	Int_t iMod=hit->GetModule();
	hlayer->Fill(lay);
	Double_t xl,yl,zl,xl0,yl0,zl0;
	Double_t xg,yg,zg,tof,xg0,yg0,zg0,tof0;
	Float_t hitloc[3],hitglo[3];
	hit->GetPositionL(xl,yl,zl,tof);
	hit->GetPositionL0(xl0,yl0,zl0,tof0);
	hit->GetPositionG(xg,yg,zg,tof);
	hit->GetPositionG0(xg0,yg0,zg0,tof0);
	Double_t hitlen=TMath::Abs(yl-yl0);
	if(hitlen<0.005) continue; // remove hits "shorter" than 50 um 
	if(lay > 1 && hitlen<0.025) continue; // remove hits "shorter" than 250 um in SDD,SSD
	hitloc[0]=0.5*(xl+xl0);
	hitloc[1]=0.5*(yl+yl0);
	hitloc[2]=0.5*(zl+zl0);
	hitglo[0]=0.5*(xg+xg0);
	hitglo[1]=0.5*(yg+yg0);
	hitglo[2]=0.5*(zg+zg0);
	Float_t edep=hit->GetIonization()*1000000;
	Float_t rad=TMath::Sqrt(hitglo[0]*hitglo[0]+hitglo[1]*hitglo[1]); 
	Float_t phi=TMath::ATan2(hitglo[1],hitglo[0]);
	hmod[lay]->Fill(iMod);
	hzl[lay]->Fill(hitloc[2]);
	hxl[lay]->Fill(hitloc[0]);
	hzg[lay]->Fill(hitglo[2]);
	hyg[lay]->Fill(hitglo[1]);
	hxg[lay]->Fill(hitglo[0]);
	hr[lay]->Fill(rad);
	hphi[lay]->Fill(phi);
	hq[lay]->Fill(edep);
	hzphi[lay]->Fill(phi,hitglo[2]);
      }
    }
  }

  gStyle->SetOptStat(10);
  gStyle->SetPadBottomMargin(0.14);


  TCanvas **c=new TCanvas*[6];
  Char_t ctit[30];
  for(Int_t iLay=0;iLay<6;iLay++){
    sprintf(name,"can%d",iLay+1);
    sprintf(ctit,"Layer %d",iLay+1);
    c[iLay]=new TCanvas(name,ctit,1200,900);
    c[iLay]->Divide(3,3,0.001,0.001);
    c[iLay]->cd(1);
    hmod[iLay]->Draw();
    c[iLay]->cd(2);
    hxl[iLay]->Draw();
    c[iLay]->cd(3);
    hzl[iLay]->Draw();
    c[iLay]->cd(4);
    hxg[iLay]->Draw();
    c[iLay]->cd(5);
    hyg[iLay]->Draw();
    c[iLay]->cd(6);
    hzg[iLay]->Draw();
    c[iLay]->cd(7);
    hr[iLay]->Draw();
    c[iLay]->cd(8);
    hphi[iLay]->Draw();    
    c[iLay]->cd(9);
    hq[iLay]->Draw();    
  }

  gStyle->SetPalette(1);
  TLatex* tstat=new TLatex();
  tstat->SetNDC();
  TCanvas* cspd=new TCanvas("cspd","SPD",1000,600);
  cspd->Divide(2,1);
  cspd->cd(1);
  hzphi[0]->Draw("colz");
  tstat->DrawLatex(0.6,0.92,Form("# Clusters = %d",int(hzphi[0]->GetEntries())));
  cspd->cd(2);
  hzphi[1]->Draw("colz");
  tstat->DrawLatex(0.6,0.92,Form("# Clusters = %d",int(hzphi[1]->GetEntries())));

  TCanvas* csdd=new TCanvas("csdd","SDD",1000,600);
  csdd->Divide(2,1);
  csdd->cd(1);
  hzphi[2]->Draw("colz");
  tstat->DrawLatex(0.6,0.92,Form("# Clusters = %d",int(hzphi[2]->GetEntries())));  
  csdd->cd(2);
  hzphi[3]->Draw("colz");
  tstat->DrawLatex(0.6,0.92,Form("# Clusters = %d",int(hzphi[3]->GetEntries())));  

  TCanvas* cssd=new TCanvas("cssd","SSD",1000,600);
  cssd->Divide(2,1);
  cssd->cd(1);
  hzphi[4]->Draw("colz");
  tstat->DrawLatex(0.6,0.92,Form("# Clusters = %d",int(hzphi[4]->GetEntries())));  
  cssd->cd(2);
  hzphi[5]->Draw("colz");
  tstat->DrawLatex(0.6,0.92,Form("# Clusters = %d",int(hzphi[5]->GetEntries())));  

}
