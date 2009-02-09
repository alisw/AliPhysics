#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TClassTable.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGeoManager.h>
#include <TH1.h>
#include <TInterpreter.h>
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliHeader.h"
#include "AliITS.h"
#include "AliITSDetTypeRec.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliRun.h"
#endif

Int_t ShowITSHitsRecPoints(Bool_t align=kFALSE,
			   TString alignfile="ITSfullv11Misalignment.root")
{
  ///////////////////////////////////////////////////////////////////////
  // Macro to check clusters and hits in the 6 ITS layers              //
  ///////////////////////////////////////////////////////////////////////

  if (gClassTable->GetID("AliRun") < 0) {
    gInterpreter->ExecuteMacro("loadlibs.C");
  }
  else { 
    if(gAlice){
      delete AliRunLoader::Instance();
      delete gAlice;
      gAlice=0;
    }
  }
  // Set OCDB if needed
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
  if(align) {
    TFile f(alignfile.Data());
    TClonesArray* ar = (TClonesArray*)f.Get("ITSAlignObjs");
    AliGeomManager::ApplyAlignObjsToGeom(*ar);
    f.Close();
  }

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0){
    cerr<<"Can not open session RL=NULL"<< endl;
    return -1;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval){
    cerr<<"LoadgAlice returned error"<<endl;
    return -1;
  }
  gAlice=rl->GetAliRun();

  retval = rl->LoadHeader();
  if (retval){
    cerr<<"LoadHeader returned error"<<endl;
    return -1;
  }


  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  if(!ITSloader){
    cerr<<"ITS loader not found"<<endl;
    return -1;
  }
  ITSloader->LoadRecPoints("read");
  ITSloader->LoadHits("read");


  Float_t cluglo[3]={0.,0.,0.}; 
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  ITS->SetTreeAddress();
  AliITSgeom *geom = ITS->GetITSgeom();
  AliITSDetTypeRec* detTypeRec = new AliITSDetTypeRec();
  detTypeRec->SetITSgeom(ITSloader->GetITSgeom());
  detTypeRec->SetDefaults();

  Int_t modmin=geom->GetStartDet(0);
  Int_t modmax=geom->GetLastDet(2);
  Int_t totmod=modmax-modmin;
  Float_t xlim[6]={4.5,7.5,16.,26.,40.,45.};
  Float_t zlim[6]={15.,15.,22.,30.,45.,55.};

  TH1F* hlayer=new TH1F("hlayer","",6,0.5,6.5);
  TH1F** hmod=new TH1F*[6];
  TH1F** hxl=new TH1F*[6];
  TH1F** hzl=new TH1F*[6];
  TH1F** hxg=new TH1F*[6];
  TH1F** hyg=new TH1F*[6];
  TH2F** hxyg=new TH2F*[6];
  TH2F** hxygHits=new TH2F*[6];
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
    sprintf(name,"hxygl%d",iLay+1);
    hxyg[iLay]=new TH2F(name,"",1000,-xlim[iLay],xlim[iLay],1000,-xlim[iLay],xlim[iLay]);
    hxyg[iLay]->GetXaxis()->SetTitle("Xglob (cm)");
    hxyg[iLay]->GetYaxis()->SetTitle("Yglob (cm)");
    sprintf(name,"hxygHitsl%d",iLay+1);
    hxygHits[iLay]=new TH2F(name,"",1000,-xlim[iLay],xlim[iLay],1000,-xlim[iLay],xlim[iLay]);
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

  TGraph **gpts=new TGraph*[3];
  TGraph2D **gpts3d=new TGraph2D*[3];
  TGraph **gRP=new TGraph*[6];
  TGraph **gHend=new TGraph*[6];
  TGraph **gHstart=new TGraph*[6];
  for(Int_t i=0;i<3;i++){
    gpts[i]=new TGraph(0);
    gpts3d[i]=new TGraph2D(0);
  }
  for(Int_t i=0;i<6;i++){
    gRP[i]=new TGraph(0);
    gRP[i]->GetXaxis()->SetTitle("global x [cm]");
    gRP[i]->GetYaxis()->SetTitle("global y [cm]");
    gHend[i]=new TGraph(0);
    gHstart[i]=new TGraph(0);
  }
  Int_t totev=rl->GetNumberOfEvents();
  printf("Total Number of events = %d\n",totev);

  Int_t iRP[6]={0},iH[6]={0};

  for(Int_t iev=0;iev<totev;iev++){
    rl->GetEvent(iev);
    TTree *TR = ITSloader->TreeR();
    TClonesArray *ITSrec  = detTypeRec->RecPoints();
    TBranch *branch = 0;
    if(TR && ITSrec){
      branch = ITSloader->TreeR()->GetBranch("ITSRecPoints");
      if(branch)branch->SetAddress(&ITSrec);
    }
    TTree *TH = ITSloader->TreeH();
    TClonesArray *hits=new TClonesArray("AliITShit",10000);
    TH->SetBranchAddress("ITS",&hits);

    Int_t nparticles = rl->GetHeader()->GetNtrack();
    cout<<"Event #"<<iev<<"   #Particles="<<nparticles<<endl;


    Int_t ipt=0;
    for(Int_t subd=0;subd<3;subd++){

      Int_t first = geom->GetStartDet(subd);
      Int_t last = geom->GetLastDet(subd);

      for (Int_t mod=first; mod<=last; mod++){
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
	    if(iev<3){
 	      gpts[iev]->SetPoint(ipt,cluglo[0],cluglo[1]);
 	      gpts3d[iev]->SetPoint(ipt,cluglo[0],cluglo[1],cluglo[2]);
	      ipt++;
	    }
	    gRP[lay]->SetPoint(iRP[lay],cluglo[0],cluglo[1]); iRP[lay]++; 
	    hmod[lay]->Fill(mod);
	    hzl[lay]->Fill(recp->GetDetLocalZ());
	    hxl[lay]->Fill(recp->GetDetLocalX());
	    hzg[lay]->Fill(cluglo[2]);
	    hyg[lay]->Fill(cluglo[1]);
	    hxyg[lay]->Fill(cluglo[0],cluglo[1]);
	    hxg[lay]->Fill(cluglo[0]);
	    hr[lay]->Fill(rad);
	    hphi[lay]->Fill(phi);
	    hq[lay]->Fill(recp->GetQ());
	  }
	}
      }
    }


    cout<<" Now read hits "<<endl;

    Int_t nentrHits=(Int_t)TH->GetEntries();
    for (Int_t i=0; i<nentrHits; i++) {
      TH->GetEvent(i);
      Int_t nhit=hits->GetEntriesFast();
      for (Int_t ih=0; ih<nhit; ih++) {
	AliITShit *h=(AliITShit*)hits->UncheckedAt(ih);
	if(h->StatusExiting()) {
	  Double_t xl,yl,zl,tl,xl0,yl0,zl0,tl0;
	  h->GetPositionL(xl,yl,zl,tl);
	  h->GetPositionL0(xl0,yl0,zl0,tl0);
	  //if(TMath::Abs(yl-yl0)<0.0290) continue;
	  hxygHits[h->GetLayer()-1]->Fill(h->GetXG(),h->GetYG());
	  gHend[h->GetLayer()-1]->SetPoint(iH[h->GetLayer()-1],h->GetXG(),h->GetYG()); 
	  Double_t x0,y0,z0,t0;
	  h->GetPositionG0(x0,y0,z0,t0);
	  gHstart[h->GetLayer()-1]->SetPoint(iH[h->GetLayer()-1],x0,y0); 
	  iH[h->GetLayer()-1]++;
	  //printf("layer %d hit length xy %f hit length yl %f\n",h->GetLayer()-1,TMath::Sqrt((x0-h->GetXG())*(x0-h->GetXG())+(y0-h->GetYG())*(y0-h->GetYG())),TMath::Abs(yl-yl0));
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
    hxyg[iLay]->Draw();    
  }

  TCanvas *cev0;
  cev0=new TCanvas("cev0","Event 0",600,600);
  gpts[0]->SetMarkerStyle(7);
  gpts[0]->SetTitle(0);
  gpts[0]->Draw("AP");

  TCanvas *cev1;
  cev1=new TCanvas("cev1","Event 1",600,600);
  gpts[1]->SetMarkerStyle(7);
  gpts[1]->SetTitle(0);
  gpts[1]->Draw("AP");

  TCanvas *cev2;
  cev2=new TCanvas("cev2","Event 2",600,600);
  gpts[2]->SetMarkerStyle(7);
  gpts[2]->SetTitle(0);
  gpts[2]->Draw("AP");

  TCanvas *chr = new TCanvas("chr","chr");
  chr->Divide(3,2);
  for(Int_t i=0;i<6;i++) {
    chr->cd(i+1);
    gHend[i]->SetMarkerStyle(7);
    gHend[i]->SetMarkerColor(1);
    gHend[i]->Draw("A,P");
    gHstart[i]->SetMarkerStyle(7);
    gHstart[i]->SetMarkerColor(4);
    gHstart[i]->Draw("P");
    gRP[i]->SetMarkerStyle(7);
    gRP[i]->SetMarkerColor(2);
    gRP[i]->Draw("P");
  }


  return 0;
}


Int_t ShowITSHitsRecPointsNtuple(Bool_t align=kFALSE,
				 TString alignfile="ITSfullv11Misalignment.root")
{
  ///////////////////////////////////////////////////////////////////////
  // Macro to check clusters and hits in the 6 ITS layers              //
  // Creates also ntuple                                               //
  ///////////////////////////////////////////////////////////////////////

  if (gClassTable->GetID("AliRun") < 0) {
    gInterpreter->ExecuteMacro("loadlibs.C");
  }
  else { 
    if(gAlice){
      delete AliRunLoader::Instance();
      delete gAlice;
      gAlice=0;
    }
  }
  // Set OCDB if needed
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
  if(align) {
    TFile f(alignfile.Data());
    TClonesArray* ar = (TClonesArray*)f.Get("ITSAlignObjs");
    AliGeomManager::ApplyAlignObjsToGeom(*ar);
    f.Close();
  }

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0){
    cerr<<"Can not open session RL=NULL"<< endl;
    return -1;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval){
    cerr<<"LoadgAlice returned error"<<endl;
    return -1;
  }
  gAlice=rl->GetAliRun();

  retval = rl->LoadHeader();
  if (retval){
    cerr<<"LoadHeader returned error"<<endl;
    return -1;
  }


  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  if(!ITSloader){
    cerr<<"ITS loader not found"<<endl;
    return -1;
  }
  ITSloader->LoadRecPoints("read");
  ITSloader->LoadHits("read");


  Float_t cluglo[3]={0.,0.,0.}; 
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  ITS->SetTreeAddress();
  AliITSgeom *geom = ITS->GetITSgeom();
  AliITSDetTypeRec* detTypeRec = new AliITSDetTypeRec();
  detTypeRec->SetITSgeom(ITSloader->GetITSgeom());
  detTypeRec->SetDefaults();

  Int_t modmin=geom->GetStartDet(0);
  Int_t modmax=geom->GetLastDet(2);
  Int_t totmod=modmax-modmin;
  Float_t xlim[6]={4.5,7.5,16.,26.,40.,45.};
  Float_t zlim[6]={15.,15.,22.,30.,45.,55.};

  TH1F* hlayer=new TH1F("hlayer","",6,0.5,6.5);
  TH1F** hmod=new TH1F*[6];
  TH1F** hxl=new TH1F*[6];
  TH1F** hzl=new TH1F*[6];
  TH1F** hxg=new TH1F*[6];
  TH1F** hyg=new TH1F*[6];
  TH2F** hxyg=new TH2F*[6];
  TH2F** hxygHits=new TH2F*[6];
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
    sprintf(name,"hxygl%d",iLay+1);
    hxyg[iLay]=new TH2F(name,"",1000,-xlim[iLay],xlim[iLay],1000,-xlim[iLay],xlim[iLay]);
    hxyg[iLay]->GetXaxis()->SetTitle("Xglob (cm)");
    hxyg[iLay]->GetYaxis()->SetTitle("Yglob (cm)");
    sprintf(name,"hxygHitsl%d",iLay+1);
    hxygHits[iLay]=new TH2F(name,"",1000,-xlim[iLay],xlim[iLay],1000,-xlim[iLay],xlim[iLay]);
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

  TGraph **gpts=new TGraph*[3];
  TGraph2D **gpts3d=new TGraph2D*[3];
  TGraph **gRP=new TGraph*[6];
  TGraph **gHend=new TGraph*[6];
  TGraph **gHstart=new TGraph*[6];
  for(Int_t i=0;i<3;i++){
    gpts[i]=new TGraph(0);
    gpts3d[i]=new TGraph2D(0);
  }
  for(Int_t i=0;i<6;i++){
    gRP[i]=new TGraph(0);
    gRP[i]->GetXaxis()->SetTitle("global x [cm]");
    gRP[i]->GetYaxis()->SetTitle("global y [cm]");
    gHend[i]=new TGraph(0);
    gHstart[i]=new TGraph(0);
  }
  Int_t totev=rl->GetNumberOfEvents();
  printf("Total Number of events = %d\n",totev);

  Int_t iRP[6]={0},iH[6]={0};

  TNtuple *nt = new TNtuple("nt","ntuple","lay:mod:deltaxl:deltax:deltay:deltaz:phi:z:xl:hlength");

  for(Int_t iev=0;iev<totev;iev++){
    rl->GetEvent(iev);
    TTree *TR = ITSloader->TreeR();
    TClonesArray *ITSrec  = detTypeRec->RecPoints();
    TBranch *branch = 0;
    if(TR && ITSrec){
      branch = ITSloader->TreeR()->GetBranch("ITSRecPoints");
      if(branch)branch->SetAddress(&ITSrec);
    }
    TTree *TH = ITSloader->TreeH();
    TClonesArray *hits=new TClonesArray("AliITShit",10000);
    TH->SetBranchAddress("ITS",&hits);

    Int_t nparticles = rl->GetHeader()->GetNtrack();
    cout<<"Event #"<<iev<<"   #Particles="<<nparticles<<endl;


    Int_t ipt=0;
    for(Int_t subd=0;subd<3;subd++){

      Int_t first = geom->GetStartDet(subd);
      Int_t last = geom->GetLastDet(subd);

      for (Int_t mod=first; mod<=last; mod++){
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
	    if(cluglo[2]>0) {gRP[lay]->SetPoint(iRP[lay],cluglo[0],cluglo[1]); iRP[lay]++;} 
	    hmod[lay]->Fill(mod);
	    hzl[lay]->Fill(recp->GetDetLocalZ());
	    hxl[lay]->Fill(recp->GetDetLocalX());
	    hzg[lay]->Fill(cluglo[2]);
	    hyg[lay]->Fill(cluglo[1]);
	    hxyg[lay]->Fill(cluglo[0],cluglo[1]);
	    hxg[lay]->Fill(cluglo[0]);
	    hr[lay]->Fill(rad);
	    hphi[lay]->Fill(phi);
	    hq[lay]->Fill(recp->GetQ());

	    Double_t hlength=0.;
	    for(Int_t jhits=0;jhits<TH->GetEntries();jhits++) {
	      TH->GetEvent(jhits);
	      Int_t nhit=hits->GetEntriesFast();
	      for (Int_t ih=0; ih<nhit; ih++) {
		AliITShit *h=(AliITShit*)hits->UncheckedAt(ih);
		if(h->GetTrack()!=recp->GetLabel(0)) continue;
		if(h->GetLayer()-1!=lay) continue;
		if(h->GetModule()!=mod) continue;
		//hxygHits[h->GetLayer()-1]->Fill(h->GetXG(),h->GetYG());
		//gHend[h->GetLayer()-1]->SetPoint(iH[h->GetLayer()-1],h->GetXG(),h->GetYG()); 
		Double_t x0,y0,z0,t0;
		Double_t x,y,z,t;
		h->GetPositionG0(x0,y0,z0,t0);
		h->GetPositionG(x,y,z,t);
		//gHstart[h->GetLayer()-1]->SetPoint(iH[h->GetLayer()-1],x0,y0); 
		//iH[h->GetLayer()-1]++;
		Double_t xl,yl,zl,tl,xl0,yl0,zl0,tl0;
		h->GetPositionL(xl,yl,zl,tl);
		h->GetPositionL0(xl0,yl0,zl0,tl0);
		hlength += TMath::Abs(yl-yl0);
		if(!h->StatusExiting()) continue;
		nt->Fill(lay,
			 mod,
			 1.e4*(0.5*(xl+xl0)-recp->GetDetLocalX()),
			 1.e4*(0.5*(x+x0)-cluglo[0]),1.e4*(0.5*(y+y0)-cluglo[1]),
			 1.e4*(0.5*(z+z0)-cluglo[2]),
			 phi,
			 cluglo[2],
			 recp->GetDetLocalX(),
			 hlength);
		hlength = 0.;
	      }
	    }


	  }
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
    hxyg[iLay]->Draw();    
  }

  TCanvas *cev0;
  cev0=new TCanvas("cev0","Event 0",600,600);
  gpts[0]->SetMarkerStyle(7);
  gpts[0]->SetTitle(0);
  gpts[0]->Draw("AP");

  TCanvas *cev1;
  cev1=new TCanvas("cev1","Event 1",600,600);
  gpts[1]->SetMarkerStyle(7);
  gpts[1]->SetTitle(0);
  gpts[1]->Draw("AP");

  TCanvas *cev2;
  cev2=new TCanvas("cev2","Event 2",600,600);
  gpts[2]->SetMarkerStyle(7);
  gpts[2]->SetTitle(0);
  gpts[2]->Draw("AP");

  TFile *out= new TFile("ntHitsRecPoints.root","recreate");
  nt->Write();
  out->Close();

  return 0;
}
