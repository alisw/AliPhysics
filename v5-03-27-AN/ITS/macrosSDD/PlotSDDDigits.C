#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TClassTable.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGeoManager.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TInterpreter.h>
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliHeader.h"
#include "AliITS.h"
#include "AliITSdigit.h"
#include "AliITSDetTypeRec.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliRun.h"
#endif

// Macro to display the SDD digits for 1 ladder
// Origin: F. Prino,   prino@to.infn.it


Int_t PlotSDDDigits(Int_t lay, Int_t lad, Int_t iev=0){
    if (gClassTable->GetID("AliRun") < 0) {
    gInterpreter->ExecuteMacro("loadlibs.C");
  }


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

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0){
    cerr<<"Can not open session RL=NULL"<< endl;
    return -1;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval){
    cerr<<"AliITSGeoPlot.C : LoadgAlice returned error"<<endl;
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
  AliITSgeom *geom = ITSloader->GetITSgeom();

  ITSloader->LoadDigits("read");
  rl->GetEvent(iev);
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  ITS->SetTreeAddress();

  TTree *TD = ITSloader->TreeD();
  AliITSDetTypeRec* detTypeRec = new AliITSDetTypeRec();
  detTypeRec->SetITSgeom(geom);
  detTypeRec->SetDefaults();
  TClonesArray *ITSdigits  = ITS->DigitsAddress(1);
  detTypeRec->ResetDigits();
  Int_t maxmod=6;
  if(lay==4) maxmod=8;
  TH2F **hdig0=new TH2F*[8];
  TH2F **hdig1=new TH2F*[8];
  Char_t histtit[100], histname[10];
  Int_t nbytes;
  for(Int_t i=0; i<maxmod; i++){
    Int_t nmod=AliITSgeomTGeo::GetModuleIndex(lay,lad,i+1);
    nbytes = TD->GetEvent(nmod);
    Int_t ndigits = ITSdigits->GetEntries();
    printf("Module %d  ndigits=%d\n",nmod,ndigits);
    sprintf(histname,"hdig%ds0",nmod);
    sprintf(histtit,"Digit mod %d side 0",nmod);
    hdig0[i]=new TH2F(histname,histtit,256,-0.5,255.5,256,-0.5,255.5);
    sprintf(histname,"hdig%ds1",nmod);
    sprintf(histtit,"Digit mod %d side 1",nmod);
    hdig1[i]=new TH2F(histname,histtit,256,-0.5,255.5,256,-0.5,255.5);
    hdig0[i]->SetStats(0);
    hdig1[i]->SetStats(0);
    for (Int_t idig=0; idig<ndigits; idig++) {
      AliITSdigit *dig=(AliITSdigit*)ITSdigits->UncheckedAt(idig);
      Int_t iz=dig->GetCoord1();  // cell number z
      Int_t ix=dig->GetCoord2();  // cell number x
      Int_t sig=dig->GetSignal();
      if(iz<256){
	hdig0[i]->SetBinContent(ix+1,iz+1,sig);
      }else{
	iz-=256;
	hdig1[i]->SetBinContent(ix+1,iz+1,sig);
      }
    }
  }
  gStyle->SetPalette(1);
  TCanvas *c1=new TCanvas("c1","",1200,900);
  c1->Divide(4,maxmod/2);
  for(Int_t i=0;i<maxmod;i++){
    c1->cd(1+2*i);
    hdig0[i]->Draw("colz");
    hdig0[i]->GetXaxis()->SetTitle("Time bin");
    hdig0[i]->GetYaxis()->SetTitle("Anode");
    c1->cd(2+2*i);
    hdig1[i]->Draw("colz");
    hdig1[i]->GetXaxis()->SetTitle("Time bin");
    hdig1[i]->GetYaxis()->SetTitle("Anode");
  }
  return 0;
}
