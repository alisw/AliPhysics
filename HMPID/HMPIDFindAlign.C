#include <TSystem.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TText.h>
#include <TF1.h>
#include <AliHMPIDParam.h>
#include <AliGeomManager.h>
#include <AliCDBEntry.h>
#include <TClonesArray.h>

TVector3 AlignAxis(Int_t ch,Int_t axis);
Double_t FindChCenterD(Int_t ch);
void FindAligners(Int_t ch, Double_t xShift,Double_t yShift, Double_t dist, AliAlignObjMatrix *matrix);
TVector3 FindNewCenter(Int_t ch, Double_t xShLoc,Double_t yShLoc, Double_t zShLoc);
void help();

//--------------------------------------------------
void HMPIDFindAlign(Int_t cham=-1)
{
  Printf("");
  Printf("******************************************");
  Printf("* You can type help() for instructions...*");
  Printf("******************************************");
  Printf("");
  
  gSystem->Sleep(3000);
  
  AliAlignObjMatrix *matrix[7];
  
  AliGeomManager::LoadGeometry("geometry.root");

  TFile *file = TFile::Open("real_geometry.root");
  AliCDBEntry *entry = (AliCDBEntry*)file->Get("AliCDBEntry");
  TClonesArray *array = (TClonesArray*)entry->GetObject(); 
  for(Int_t h=0; h<7; h++) matrix[h] = (AliAlignObjMatrix*)array->At(h);
  AliGeomManager::ApplyAlignObjsToGeom(*array);      
  AliHMPIDParam::Instance(); //just to print at the beginning...

  Int_t ichMin = 0;
  Int_t ichMax = 6;         
  if(cham>=0&&cham<=6) { ichMin = cham; ichMax = cham;}
  
  for(Int_t ich=ichMin;ich<=ichMax;ich++) {
    
    Printf("");
    Printf("*********SUMMARY FOR CHAMBER %i **********",ich);
    Printf ("Distance from IP of chamber %i = %f (cm)",ich,FindChCenterD(ich));
    
    TVector3 xAl = AlignAxis(ich,0);
    TVector3 yAl = AlignAxis(ich,1);
    
    Double_t p1 = 1/(xAl.Z()*xAl.Z());
    Double_t p2 = 1/(yAl.Z()*yAl.Z());
    
    Double_t distMean = (p1*xAl.Y()+p2*yAl.Y())/(p1+p2);
      
    FindAligners(ich,xAl.X(),yAl.X(),distMean,matrix[ich]); 
  }
}      
//--------------------------------------------------
TVector3 AlignAxis(Int_t ch,Int_t axis)
{
  
// Axis = 0  X Local Axis;
// Axis = 1  Y Local Axis;
  
  
  TH1F* delta[6];
  TH2F *map2D = new TH2F(Form("map2DCham%i",ch),Form("map2DCham%i",ch),520,0,130,520,0,130);
  map2D->SetStats(kFALSE);
  
  TFile* file = TFile::Open("histos.root","old"); 
  TTree *pNtuple = (TTree*)file->Get("Tree")); 
      
  for(Int_t j=0;j<6;j++) {
    if(axis == 0) {
      delta[j] = new TH1F(Form("delta%i",j),Form("delta%i",j),100,-25.,25.);
    } else {
      delta[j] = new TH1F(Form("delta%i",j),Form("delta%i",j),100,-25.,25.);
    }
  }
  Double_t mean[6];
  Double_t sigma[6];
  Double_t siz[6];
  Double_t bin[6];
  
  TCanvas *c = new TCanvas();

  Int_t index = 0;  
  for(Int_t j=0;j<6;j++) {
    Double_t xMin = j*20;
    Double_t xMax = xMin+20;
    if(axis == 0) {
      pNtuple->Draw(Form("Xpc-X>>delta%i",j),Form("Chamber==%i&&Charge>100&&X>%f&&X<%f&&NumTPCclust>70",ch,xMin,xMax),"0");
    } else {
      pNtuple->Draw(Form("Ypc-Y>>delta%i",j),Form("Chamber==%i&&Charge>100&&Y>%f&&Y<%f&&NumTPCclust>70",ch,xMin,xMax),"0");
    }
    if(delta[j]->GetEntries()<50)  continue;

    Double_t dIntFit = 5.;
    siz[index] = 0.5*(xMax-xMin);
    bin[index]= 0.5*(xMax+xMin);
    Double_t xMaxPos = delta[j]->GetBinCenter(delta[j]->FindFirstBinAbove(delta[j]->GetMaximum()-1));
    delta[j]->Fit("gaus","Q","",xMaxPos-dIntFit,xMaxPos+dIntFit);
    mean[index]  = (delta[j]->GetFunction("gaus"))->GetParameter(1);
    sigma[index] = (delta[j]->GetFunction("gaus"))->GetParError(1);
    index++;
  }
  c->Clear();
  c->Divide(3,3);
  

  for(Int_t j=0;j<6;j++) {
    c->cd(j+1);delta[j]->Draw();
  }
      
  c->cd(8);
  
  TH1F *h1 = new TH1F("h1","h1",6,0,120);
      
  h1->SetTitle(Form("LRS for chamber %i",ch));
  if(axis == 0) {
    h1->GetXaxis()->SetTitle("Xmip (cm)");
    h1->GetYaxis()->SetTitle("deltaX (cm)");
  } else {  
    h1->GetXaxis()->SetTitle("Ymip (cm)");
    h1->GetYaxis()->SetTitle("DeltaY (cm)");
  }
  h1->SetMinimum(-5.);
  h1->SetMaximum( 5.);
  
  for(Int_t jBin=0; jBin<6; jBin++) {h1->SetBinContent(jBin+1,mean[jBin]); h1->SetBinError(jBin+1,sigma[jBin]);}
  //pGr->SetMinimum(TMath::MinElement(index,mean)-1.);
  //pGr->SetMaximum(TMath::MaxElement(index,mean)+1.);
  h1->SetMarkerStyle(20);
  h1->Draw("E");
  h1->Fit("pol1","Q");
//  Double_t  shift      = (h1->GetFunction("pol1"))->GetParameter(0);
  Double_t  shift      = (h1->GetFunction("pol1"))->Eval(60.);
  Double_t  errshift   = (h1->GetFunction("pol1"))->GetParError(0);
  Double_t  coefang    = (h1->GetFunction("pol1"))->GetParameter(1);
  Double_t  errcoefang = (h1->GetFunction("pol1"))->GetParError(1);
  Double_t  prob       = (h1->GetFunction("pol1"))->GetProb();
  Double_t dist = FindChCenterD(ch);
  
  Double_t rfShift    = -coefang/(1.+coefang)*dist;
  Double_t errrfShift = TMath::Abs(dist/((1.+coefang)*(1+coefang))*errcoefang);
  if (axis == 0) {
    Printf(" x shift = %5.3f +/- %5.3f dist shift = %5.3f +/- %5.3f prob chi2 %5.1f",shift,errshift,rfShift,errrfShift,prob*100);
  } else {
    Printf(" y shift = %5.3f +/- %5.3f dist shift = %5.3f +/- %5.3f prob chi2 %5.1f",shift,errshift,rfShift,errrfShift,prob*100);
  }

  c->cd(7);
  
  TPaveText *pPa = new TPaveText();
  pPa->DrawPave(0.2,0.2,0.8,0.8);
  TText *t = new TText();
  t->SetTextSize(0.08);
  t->DrawText(0.21,0.70,Form(" shift %5.3f +/- %5.3f",shift,errshift));
  t->DrawText(0.21,0.50,Form(" dist shift %5.3f +/- %5.3f",rfShift,errrfShift));
  t->DrawText(0.21,0.30,Form(" dist from IP %5.3f ",dist));
  
  
  c->cd(9); 
  map2D->SetTitle("2D MAP");
  pNtuple->Draw("Y:X>>map2D",Form("Chamber==%i&&Charge>100&&NumTPCclust>70",ch),"0");  
  map2D->SetTitle("2D MAP");

  
//  c1->Clear();
  if (axis == 0) {
    c->SaveAs(Form("ShiftXChamber%i.png",ch));
  } else {
    c->SaveAs(Form("ShiftYChamber%i.png",ch));
  }
    
  TVector3 align(shift,rfShift,errrfShift);
  return align;
}
//------------------------------------------------------------------
Double_t FindChCenterD(Int_t ch)
{
  
  AliHMPIDParam *pParam = AliHMPIDParam::Instance();
  TVector3 xyz = pParam->Lors2Mars(ch,0.5*pParam->SizeAllX(),0.5*pParam->SizeAllY());
  return xyz.Mag();
}
//------------------------------------------------------------------
void FindAligners(Int_t ch, Double_t xShift,Double_t yShift, Double_t dist, AliAlignObjMatrix *matrix) 
{
  TVector3 oldC = FindNewCenter(ch,0,0,0);
    
  Printf("");
  Printf("Old X Y Z of the chamber n. %i",ch);
  Printf("XC = %f",oldC.X());
  Printf("YC = %f",oldC.Y());
  Printf("ZC = %f",oldC.Z());
  Printf("Dist. from IP = %f",oldC.Mag());
  Printf("");
  

  Printf("Local X shift = %f",xShift);
  Printf("Local Y shift = %f",yShift);
  Printf("RPhi    shift = %f",dist);
    
  TVector3 newC = FindNewCenter(ch,xShift,yShift,dist);
  
  Printf("");
  Printf("");
  Printf("(next line are the current values in the Align Object for chamber %i)",ich);
  Printf("");
  matrix->Print();
  
  Double_t tr[] = {0.,0.,0.};
  
  matrix->GetTranslation(tr);
  
  Printf("");
  Printf("New X Y Z of the chamber n. %i",ch);
  Printf("");
  Printf("XC = %f",newC.X());
  Printf("YC = %f",newC.Y());
  Printf("ZC = %f",newC.Z());
  Printf("Dist. from IP = %f",newC.Mag());
  Printf("");
  Printf("shift with respect to the old chamber center");
  Printf("");
  Printf("XC shift = %5.2f YC shift = %5.2f ZC shift = %5.2f  for chamber %i",newC.X()-oldC.X(),newC.Y()-oldC.Y(),newC.Z()-oldC.Z(),ch);
  Printf("");
  Printf("next line are the values to be inserted in the Align Object");
  Printf("");
  Printf("XC shift = %5.2f YC shift = %5.2f ZC shift = %5.2f  for chamber %i",newC.X()-oldC.X() + tr[0],newC.Y()-oldC.Y() + tr[1],newC.Z()-oldC.Z() + tr[2],ch);
  
}  
//------------------------------------------------------------------
TVector3 FindNewCenter(Int_t ch, Double_t xShLoc,Double_t yShLoc, Double_t zShLoc) 
{
//First shift in (LocX,LocY) and go in MRS
//Then shift radially...
  AliHMPIDParam *pParam = AliHMPIDParam::Instance();
  TVector3 xyz = pParam->Lors2Mars(ch,0.5*pParam->SizeAllX()+xShLoc,0.5*pParam->SizeAllY()+yShLoc,-1);
  TVector3 xyzSh(1,1,1);
  xyzSh.SetTheta(xyz.Theta());
  xyzSh.SetPhi(xyz.Phi());
  xyzSh.SetMag(xyz.Mag()+zShLoc);
  return xyzSh;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void help()
{
  Printf("Instructions to run this macros:");
  Printf("    1 - In the same dir where the macro runs it must exist the file geometry.root");
  Printf("    2 - assign symbolic links as it follows:");
  Printf("       a -  file containing the alignment objects in OCDB/HMPID/Align/Data linked to real_geometry.root");
  Printf("       b -  file with the NTupla of ESD Data (produced from AliHMPIDAnalysisTask.C) linked as histos.root");
  Printf("    ");
  Printf("   The macro have to be run in the following way:");
  Printf("    ");
  Printf("    > aliroot");
  Printf("    ");
  Printf("    AliRoot> .L HMPIDFindAlign.C");
  Printf("    AliRoot> HMPIDFindAlign(); summary.txt");
  Printf("    AliRoot> .q");
  Printf("    ");
  Printf("The results of the alignment procedures is in the text file summary.txt");
  Printf("In addition 14 Canvases (7 chamb. x X and Y) are produced to show the residual distrib. in Local Coord. and other 7 Canvases with the 2D map on for chamber");
  Printf("if some plots are funny, please try to run for a signle chamber modifying some parameters in the macro:");
  Printf("until you find a reasonable result.");
  Printf("    ");
  Printf("AliRoot> HMPIDFindAlign(chamber)");
  Printf("    ");
  Printf("Ciao!."); 
}  
