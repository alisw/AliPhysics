#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TGrid.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TLine.h>
#include <TView.h>
#include <TMath.h>
#include <TGeoManager.h>
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliITSAlignMille2Module.h"
#endif

void Draw3D(AliITSOnlineCalibrationSPDhandler *h);

void ShowSPDConfiguration(Int_t runNb=0, const char *ocdblocation="local://$ALICE_ROOT/OCDB", bool grid=kFALSE,bool threed=kFALSE){


  gStyle->SetOptStat(0);

  if(grid){
    TGrid::Connect("alien://");
    if(!gGrid){
      printf("no grid connection is available, exiting.\n");
      return;
    }
  }
  AliITSOnlineCalibrationSPDhandler *h = new AliITSOnlineCalibrationSPDhandler();
  h->ReadDeadFromDB(runNb,ocdblocation);
  AliCDBManager::Instance();
  AliCDBManager::Instance()->SetRun(runNb);
  AliCDBManager::Instance()->SetDefaultStorage(ocdblocation);
  AliGeomManager::LoadGeometry();

  if(threed) {
    Draw3D(h);  
    return;
  }

  Double_t nact[2]={0.,0.};

  TCanvas *c = new TCanvas("c","Active Modules ",500,700);
  c->Divide(1,2);
  c->cd(1);
  TH2D *hPhiZInner = new TH2D("hPhiZInner","Inner layer Active Modules ",200,-20,20,3,0,2*TMath::Pi());
  hPhiZInner->SetXTitle("Z (cm)");
  hPhiZInner->SetYTitle("#varphi (rad)");
  hPhiZInner->Draw();

  for(Int_t i=0; i<80; i++){
    if((h->GetNrBad(i))<1) {
      TGeoHMatrix matrix;
      int vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(i);
      AliITSAlignMille2Module::SensVolMatrix(vid,&matrix);
      Double_t local0[3],local1[3],local2[3],local3[3]; // local position of the four angles
      local0[0]=-0.6375; local0[1]=0; local0[2]= 3.48; 
      local1[0]=-0.6375; local1[1]=0; local1[2]= -3.48; 
      local2[0]=0.6375; local2[1]=0; local2[2]= -3.48; 
      local2[0]=0.6375; local3[1]=0; local3[2]= 3.48; 
      Double_t global0[3],global1[3],global2[3],global3[3];
      matrix.LocalToMaster(local0,global0);
      matrix.LocalToMaster(local1,global1);
      matrix.LocalToMaster(local2,global2);
      matrix.LocalToMaster(local3,global3);
      Double_t phiUp = atan2(global0[1],global0[0]);
      if(phiUp<0) phiUp+=2*TMath::Pi();
      Double_t phiDown = atan2(global2[1],global2[0]);
      if(phiDown<0) phiDown+=2*TMath::Pi();
      TLine *lhor1 = new TLine(global0[2],phiDown,global1[2],phiDown); lhor1->Draw("same");
      lhor1->SetLineColor(kBlue);
      lhor1->SetLineWidth(3); 
      TLine *lver1 = new TLine(global1[2],phiDown,global2[2],phiUp); lver1->Draw("same");
      lver1->SetLineColor(kBlue);
      lver1->SetLineWidth(3); 
      TLine *lhor2 = new TLine(global2[2],phiUp,global3[2],phiUp); lhor2->Draw("same");
      lhor2->SetLineColor(kBlue);
      lhor2->SetLineWidth(3); 
      TLine *lver2 = new TLine(global3[2],phiUp,global0[2],phiDown); lver2->Draw("same");
      lver2->SetLineColor(kBlue);
      lver2->SetLineWidth(3); 
      nact[0]++;
    } 
  }
  c->cd(2);

  TH2D *hPhiZOuter = new TH2D("hPhiZOuter","Outer layer Active Modules ",200,-20,20,3,0,2*TMath::Pi());
  hPhiZOuter->SetXTitle("Z (cm)");
  hPhiZOuter->SetYTitle("#varphi (rad)");
  hPhiZOuter->Draw();
  for(Int_t i=80; i<240; i++){
    TGeoHMatrix matrix;
    int vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(i);
    AliITSAlignMille2Module::SensVolMatrix(vid,&matrix);
    Double_t local[4][3]; // local position of the four angles
    local[0][0]=0.6375; local[0][1]=0; local[0][2]= -3.48;
    local[1][0]=0.6375; local[1][1]=0; local[1][2]= 3.48;
    local[2][0]=-0.6375; local[2][1]=0; local[2][2]= 3.48;
    local[2][0]=-0.6375; local[3][1]=0; local[3][2]= -3.48;
    Double_t global[4][3];
    for(Int_t j=0; j<4; j++){
      matrix.LocalToMaster(local[j],global[j]);
    }
    Double_t phiUp = atan2(global[0][1],global[0][0]);
    if(phiUp<0) phiUp+=2*TMath::Pi();
    Double_t phiDown = atan2(global[2][1],global[2][0]);
    if(phiDown<0) phiDown+=2*TMath::Pi();
    if(i>235) if(phiUp < 0.1) phiUp = TMath::Pi()*2;
    //  printf("module %i  -   phiDown %f   phiUp %f \n",i,phiDown,phiUp); 
    if((h->GetNrBad(i))<1) {
      TLine *lhor1 = new TLine(global[0][2],phiUp,global[1][2],phiUp); lhor1->Draw("same");
      lhor1->SetLineColor(kBlue);
      lhor1->SetLineWidth(2); 
      TLine *lver1 = new TLine(global[1][2],phiUp,global[2][2],phiDown); lver1->Draw("same");
      lver1->SetLineColor(kBlue);
      lver1->SetLineWidth(2); 
      TLine *lhor2 = new TLine(global[2][2],phiDown,global[3][2],phiDown); lhor2->Draw("same");
      lhor2->SetLineColor(kBlue);
      lhor2->SetLineWidth(2); 
      TLine *lver2 = new TLine(global[3][2],phiDown,global[0][2],phiUp); lver2->Draw("same");
      lver2->SetLineColor(kBlue);
      lver2->SetLineWidth(2); 
      nact[1]++; 
    } 
  }
  printf("  \n   Number of Active SPD modules (->Total)  : inner %3.0f (80) %f   -  outer %3.0f (160)  %f \n",nact[0],nact[0]/80.,nact[1],nact[1]/160.);
  c->SaveAs(Form("active%i.png",runNb));
}


void Draw3D(AliITSOnlineCalibrationSPDhandler *h){

  TGeoHMatrix m2t[240];
  for(Int_t imod=0; imod<240; imod++){
    int vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(imod);
    AliITSAlignMille2Module::SensVolMatrix(vid,&m2t[imod]);
  }

  delete gGeoManager;

  new TGeoManager("SPD","active");

  TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
  TGeoMedium *none = new TGeoMedium("Vacuum",0,vacuum);
  TGeoVolume *top = gGeoManager->MakeBox("TOP",none,500,500,500);
  gGeoManager->SetTopVolume(top);

  TGeoVolume *ladder = gGeoManager->MakeBox("ladder",none,0.6375,0.001/2,3.48);

  Int_t nActive[2]={0,0};
  for(Int_t imod=0; imod<240; imod++){
    TGeoRotation *rot  = new TGeoRotation();
    rot->SetMatrix(m2t[imod].GetRotationMatrix());
    TGeoCombiTrans *matrix = new TGeoCombiTrans(m2t[imod].GetTranslation()[0],m2t[imod].GetTranslation()[1],m2t[imod].GetTranslation()[2],rot);
    if((40960-h->GetNrBad(imod))>0) {
      top->AddNode(ladder,imod,matrix);
      if(imod<80) nActive[0]++;
      else nActive[1]++;
    }
  }

  printf("  \n\n   Number of Active SPD modules (->Total)  : inner %i (80) outer %i (160) \n\n\n",nActive[0],nActive[1]);
  gGeoManager->CloseGeometry();
  top->Draw("ogl");
  gPad->GetView()->ShowAxis();

}

