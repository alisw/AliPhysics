#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TClassTable.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGeoManager.h>
#include <TH1.h>
#include <TInterpreter.h>
#include "AliGeomManager.h"
#include "AliHeader.h"
#include "AliITS.h"
#include "AliITSDetTypeRec.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliRun.h"
#endif

Int_t ShowITSRecPoints(Int_t nevfordisp=0){
  ///////////////////////////////////////////////////////////////////////
  // Macro to check clusters in the 6 ITS layers                       //
  // Provides:                                                         //
  //  6 canvases with 9 plots each (1 canvas for each layer)           //
  //  3 canvases with cluster XY coordinates for the first 3 events    //
  ///////////////////////////////////////////////////////////////////////

  if (gClassTable->GetID("AliRun") < 0) {
    gInterpreter->ExecuteMacro("loadlibs.C");
  }
 
  // retrives geometry 
  if(!gGeoManager){
    AliGeomManager::LoadGeometry("geometry.root");
  }

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0){
    cerr<<"Can not open session RL=NULL"<< endl;
    return -1;
  }

  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  if(!ITSloader){
    cerr<<"ITS loader not found"<<endl;
    return -1;
  }
  ITSloader->LoadRecPoints("read");

  Float_t cluglo[3]={0.,0.,0.}; 
  AliITSDetTypeRec* detTypeRec = new AliITSDetTypeRec();

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
  }

  TGraph *gptsXY=new TGraph(0);
  TGraph *gptsRZ=new TGraph(0);

  Int_t totev=rl->GetNumberOfEvents();
  printf("Total Number of events = %d\n",totev);

  for(Int_t iev=0;iev<totev;iev++){
    rl->GetEvent(iev);
    TTree *TR = ITSloader->TreeR();
    TClonesArray *ITSrec  = detTypeRec->RecPoints();
    TBranch *branch = 0;
    if(TR && ITSrec){
      branch = ITSloader->TreeR()->GetBranch("ITSRecPoints");
      if(branch)branch->SetAddress(&ITSrec);
    }
    if(iev%100==0) printf("Event #%d\n",iev);


    Int_t ipt=0;
    for (Int_t mod=modmin; mod<=modmax; mod++){
      detTypeRec->ResetRecPoints();
      branch->GetEvent(mod);
      Int_t nrecp = ITSrec->GetEntries();
      if(nrecp>0){
	for(Int_t irec=0;irec<nrecp;irec++) {
	  AliITSRecPoint *recp = (AliITSRecPoint*)ITSrec->At(irec);
	  Int_t lay=recp->GetLayer();
	  hlayer->Fill(lay);
	  recp->GetGlobalXYZ(cluglo);
	  Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
	  Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
	  if(iev==nevfordisp){
	    gptsXY->SetPoint(ipt,cluglo[0],cluglo[1]);
	    if(cluglo[1]>0) gptsRZ->SetPoint(ipt,cluglo[2],rad);
	    else gptsRZ->SetPoint(ipt,cluglo[2],-rad);
	    ipt++;
	  }
	  hmod[lay]->Fill(mod);
	  hzl[lay]->Fill(recp->GetDetLocalZ());
	  hxl[lay]->Fill(recp->GetDetLocalX());
	  hzg[lay]->Fill(cluglo[2]);
	  hyg[lay]->Fill(cluglo[1]);
	  hxg[lay]->Fill(cluglo[0]);
	  hr[lay]->Fill(rad);
	  hphi[lay]->Fill(phi);
	  hq[lay]->Fill(recp->GetQ());
	}
      }
    }
  }
  
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

  TCanvas *cev0;
  sprintf(ctit,"Event %d XY",nevfordisp);
  cev0=new TCanvas("cev0",ctit,600,600);
  if(gptsXY->GetN()>0){
    gptsXY->SetMarkerStyle(7);
    gptsXY->SetTitle(0);
    gptsXY->Draw("AP");
    gptsXY->GetXaxis()->SetTitle("Xglob");
    gptsXY->GetYaxis()->SetTitle("Yglob");
   }

  TCanvas *cev1;
  sprintf(ctit,"Event %d Zr",nevfordisp);
  cev1=new TCanvas("cev1",ctit,600,600);
  if(gptsRZ->GetN()>0){
    gptsRZ->SetMarkerStyle(7);
    gptsRZ->SetTitle(0);
    gptsRZ->Draw("AP");
    gptsRZ->GetXaxis()->SetTitle("Zglob");
    gptsRZ->GetYaxis()->SetTitle("Radius");
  }

  return 0;
}


