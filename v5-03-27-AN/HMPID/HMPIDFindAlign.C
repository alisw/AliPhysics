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
void FindAligners(Int_t ch, Double_t xShift,Double_t yShift, Double_t dist);
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
  
  AliGeomManager::LoadGeometry("geometry.root");

  TFile *file = TFile::Open("real_geometry.root");
  if(file) {
    AliCDBEntry *entry = (AliCDBEntry*)file->Get("AliCDBEntry");
    TClonesArray *array = (TClonesArray*)entry->GetObject();
    AliGeomManager::ApplyAlignObjsToGeom(*array);
}      
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
  
    FindAligners(ich,xAl.X(),yAl.X(),distMean); 
  }
}      
//--------------------------------------------------
TVector3 AlignAxis(Int_t ch,Int_t axis)
{
  
// Axis = 0  X Local Axis;
// Axis = 1  Y Local Axis;
  
  
  TH1F* delta[6];
  
  TFile* file = TFile::Open("histos.root","old"); 
  TNtuple *pNtuple;
  if(!(TNtuple*)file->Get("nTupla")) {
// new version because Giacomo uses now a Task to strio the data...
  TList *hmpoutput = (TList*)(file->FindObjectAny("hmpoutput"));
  pNtuple = (TNtuple*)hmpoutput->At(0);
} else {
  pNtuple = (TNtuple*)file->Get("nTupla");
}
  
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
      pNtuple->Draw(Form("Xtrk-Xmip>>delta%i",j),Form("ch==%i&&mipCharge>50&&Xmip>%f&&Xmip<%f",ch,xMin,xMax),"0");
    } else {
      pNtuple->Draw(Form("Ytrk-Ymip>>delta%i",j),Form("ch==%i&&mipCharge>50&&Ymip>%f&&Ymip<%f",ch,xMin,xMax),"0");
    }
    if(delta[j]->GetEntries()<50)  continue;

    Double_t dIntFit = 5.;
    siz[index] = 0.5*(xMax-xMin);
    bin[index]= 0.5*(xMax+xMin);
    Double_t xMaxPos = delta[j]->GetBinCenter(delta[j]->FindFirstBinAbove(delta[j]->GetMaximum()-1));
    delta[j]->Fit("gaus","Q0","",xMaxPos-dIntFit,xMaxPos+dIntFit);
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
  
  TGraphErrors *pGr = new TGraphErrors(index,bin,mean,siz,sigma);
  
  pGr->SetTitle(Form("LRS for chamber %i",ch));
  if(axis == 0) {
    pGr->GetXaxis()->SetTitle("Xmip (cm)");
    pGr->GetYaxis()->SetTitle("deltaX (cm)");
  } else {  
    pGr->GetXaxis()->SetTitle("Ymip (cm)");
    pGr->GetYaxis()->SetTitle("DeltaY (cm)");
  }
  pGr->SetMinimum(-5.);
  pGr->SetMaximum( 5.);
  //pGr->SetMinimum(TMath::MinElement(index,mean)-1.);
  //pGr->SetMaximum(TMath::MaxElement(index,mean)+1.);
  pGr->Draw("ALP");
  pGr->Fit("pol1","Q");  
  Double_t  shift      = (pGr->GetFunction("pol1"))->GetParameter(0);
  Double_t  errshift   = (pGr->GetFunction("pol1"))->GetParError(0);
  Double_t  coefang    = (pGr->GetFunction("pol1"))->GetParameter(1);
  Double_t  errcoefang = (pGr->GetFunction("pol1"))->GetParError(1);
  Double_t  prob       = (pGr->GetFunction("pol1"))->GetProb();
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
  
  if (axis == 0) {
    c->SaveAs(Form("ShiftXChamber%i.gif",ch));
  } else {
    c->SaveAs(Form("ShiftYChamber%i.gif",ch));
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
void FindAligners(Int_t ch, Double_t xShift,Double_t yShift, Double_t dist) 
{
  TVector3 old = FindNewCenter(ch,0,0,0);
    
  Printf("");
  Printf("Old X Y Z of the chamber n. %i",ch);
  Printf("XC = %f",old.X());
  Printf("YC = %f",old.Y());
  Printf("ZC = %f",old.Z());
  Printf("Dist. from IP = %f",old.Mag());
  Printf("");
  

  Printf("Local X shift = %f",xShift);
  Printf("Local Y shift = %f",yShift);
  Printf("RPhi    shift = %f",dist);
    
  TVector3 shift = FindNewCenter(ch,xShift,yShift,dist);

  Printf("");
  Printf("New X Y Z of the chamber n. %i",ch);
  Printf("");
  Printf("XC = %f YC = %f ZC = %f",shift.X(),shift.Y(),shift.Z());
  Printf("Dist. from IP = %f",shift.Mag());
  Printf("");
  Printf("(next line are the values to be inserted in the Align Object)");
  Printf("");
  Printf("XC shift = %5.2f YC shift = %5.2f ZC shift = %5.2f  for chamber %i",shift.X()-old.X(),shift.Y()-old.Y(),shift.Z()-old.Z(),ch);
  Printf("");
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
  Printf("       a -  file from OCDB/HMPID/Align/Data linked to real_geometry.root");
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
  Printf("In addition 14 Canvases (7 chamb. x X and Y) are produced to show the residual distrib. in Local Coord.");
  Printf("if some plots are funny, please try to run for a signle chamber modifying some parameters in the macro:");
  Printf("until you find a reasonable result.");
  Printf("    ");
  Printf("AliRoot> HMPIDFindAlign(chamber)");
  Printf("    ");
  Printf("Ciao!."); 
}  
