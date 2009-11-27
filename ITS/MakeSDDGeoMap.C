#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TClassTable.h>
#include <TInterpreter.h>
#include <TGeoManager.h>
#include <TNtuple.h>
#include "AliRun.h"
#include "AliGeomManager.h"
#include "AliITSgeomTGeo.h"
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#endif

void MakeSDDGeoMap(){
///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Macro to Create SDD geometrical map starting from geometry.root file  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

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
  // retrives geometry 
  if(!gGeoManager){
    AliGeomManager::LoadGeometry("geometry.root");
  }
  
  TNtuple *ntsddgeo=new TNtuple("ntsddgeo","SDD module positions","iMod:lay:lad:det:x:y:z:r:theta:phi");
  Float_t xnt[10];
  FILE *outfil;
  outfil=fopen("SDDgeomap.data","w");
  
  Int_t first = AliITSgeomTGeo::GetModuleIndex(3,1,1);
  Int_t last = AliITSgeomTGeo::GetModuleIndex(5,1,1)-1;
  Double_t pos[3];
  Int_t lay,lad,det;
  //  TPaveText **text=new TPaveText*[36];
  TLatex **text=new TLatex*[36];
  TLatex **ltext=new TLatex*[36];
  TLine **linxy=new TLine*[36];
  Char_t modtxt[3];
  Char_t ladtxt[8];
  Int_t it=0;
  Float_t D=3.5; //SDD half length along drift
  for (Int_t iMod=first; iMod<=last; iMod++){
    fprintf(outfil,"=========================================================\n");
    AliITSgeomTGeo::GetModuleId(iMod,lay,lad,det);
    AliITSgeomTGeo::GetTranslation(iMod,pos);
    Float_t rad=TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    Float_t theta=90.-TMath::ATan2(pos[2],rad)*TMath::RadToDeg();
    Float_t phi=TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
    fprintf(outfil,"ModuleId=%d   --- Layer %d  Ladder %d  Det %d\n",iMod,lay,lad,det);
    fprintf(outfil,"Global coordinates: x=%f \t y=%f \t z=%f\n",pos[0],pos[1],pos[2]);
    fprintf(outfil,"                    r=%f \t theta=%f \t phi=%f\n",rad,theta,phi);
    xnt[0]=(Float_t)iMod;
    xnt[1]=(Float_t)lay;    
    xnt[2]=(Float_t)lad;    
    xnt[3]=(Float_t)det;
    xnt[4]=pos[0];
    xnt[5]=pos[1];
    xnt[6]=pos[2];
    xnt[7]=rad;
    xnt[8]=theta;
    xnt[9]=phi;
    ntsddgeo->Fill(xnt);
    sprintf(modtxt,"%d",iMod);
    sprintf(ladtxt,"L%d",lad);
    if(det==1){ 
      text[it]=new TLatex(1.1*pos[0],1.1*pos[1],modtxt);
      text[it]->SetTextAlign(22);
      text[it]->SetTextColor(2);
      text[it]->SetTextSize(0.03);      
      ltext[it]=new TLatex(.9*pos[0],.9*pos[1],ladtxt);
      ltext[it]->SetTextAlign(22);
      ltext[it]->SetTextColor(4);
      ltext[it]->SetTextSize(0.03);
      Float_t deltax=D*TMath::Sin(phi*TMath::DegToRad());
      Float_t deltay=D*TMath::Cos(phi*TMath::DegToRad());
      linxy[it]=new TLine(pos[0]-deltax,pos[1]+deltay,pos[0]+deltax,pos[1]-deltay);
      linxy[it]->SetLineWidth(2);
      it++;
    }
  }
  fclose(outfil);
  ntsddgeo->SetMarkerStyle(7);
  TCanvas *c1;
  c1=new TCanvas("c1","",800,800);
  ntsddgeo->Draw("y:x");
  for(Int_t i=0;i<36;i++){
    linxy[i]->Draw();
    text[i]->Draw();
    ltext[i]->Draw();
  }
  

  TCanvas *c2=new TCanvas("c2","",1200,800);
  c2->Divide(3,2);
  c2->cd(1);
  ntsddgeo->SetMarkerColor(1);
  ntsddgeo->Draw("x:iMod");
  ntsddgeo->SetMarkerColor(2);
  ntsddgeo->Draw("x:iMod","lay==4","same");
  c2->cd(2);
  ntsddgeo->SetMarkerColor(1);
  ntsddgeo->Draw("y:iMod");
  ntsddgeo->SetMarkerColor(2);
  ntsddgeo->Draw("y:iMod","lay==4","same");
  c2->cd(3);
  ntsddgeo->SetMarkerColor(1);
  ntsddgeo->Draw("z:iMod");
  ntsddgeo->SetMarkerColor(2);
  ntsddgeo->Draw("z:iMod","lay==4","same");
  c2->cd(4);
  ntsddgeo->SetMarkerColor(1);
  ntsddgeo->Draw("r:iMod");
  ntsddgeo->SetMarkerColor(2);
  ntsddgeo->Draw("r:iMod","lay==4","same");
  TLatex *t3=new TLatex(0.68,0.45,"Layer 3");
  t3->SetTextSize(0.05);
  t3->SetNDC();
  t3->Draw();
  TLatex *t4=new TLatex(0.68,0.37,"Layer 4");
  t4->SetTextSize(0.05);
  t4->SetNDC();
  t4->SetTextColor(2);
  t4->Draw();
  c2->cd(5);
  ntsddgeo->SetMarkerColor(1);
  ntsddgeo->Draw("phi:iMod");
  ntsddgeo->SetMarkerColor(2);
  ntsddgeo->Draw("phi:iMod","lay==4","same");
  t3->Draw();
  t4->Draw();
  c2->cd(6);
  ntsddgeo->SetMarkerColor(1);
  ntsddgeo->Draw("theta:iMod");
  ntsddgeo->SetMarkerColor(2);
  ntsddgeo->Draw("theta:iMod","lay==4","same");
  t3->Draw();
  t4->Draw();
  
  TFile *f=new TFile("SDDgeomap.root","recreate");
  f->cd();
  ntsddgeo->Write();
  f->Close();
  
}


