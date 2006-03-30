//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************
#include "AliRICHParam.h" //class header
#include "AliESD.h"
#include <TCanvas.h>      //TestXXX() 
#include <TLatex.h>
#include <THStack.h>
#include <TLegend.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TPolyLine.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TRotation.h>
#include <AliCDBManager.h> //CdbRead()
#include <AliCDBStorage.h> //CdbRead()
#include <AliCDBEntry.h>   //CdbRead()
#include <AliRunLoader.h>  //Stack()
#include <AliStack.h>      //Stack()
#include <TParticle.h>     //Stack()    
#include "AliRICHHelix.h"  //TestTrans()

ClassImp(AliRICHParam)
AliRICHParam * AliRICHParam::fgInstance             =0x0;     //singleton pointer               
Bool_t         AliRICHParam::fgIsWireSag            =kTRUE;   //take ware sagita into account?
Bool_t         AliRICHParam::fgIsResolveClusters    =kTRUE;   //do cluster resolving?
Bool_t         AliRICHParam::fgIsFeedback           =kTRUE;   //generate feedback photons?
Bool_t         AliRICHParam::fgIsTestBeam           =kFALSE;  //special test beam configuration

Int_t    AliRICHParam::fgHV[kNsectors]        ={2050,2050,2050,2050,2050,2050};
Int_t    AliRICHParam::fgNsigmaTh             =4;
Float_t  AliRICHParam::fgSigmaThMean          =1.132; //QDC 
Float_t  AliRICHParam::fgSigmaThSpread        =0.035; //     
Double_t AliRICHParam::fgErrChrom[4][330];                       //
Double_t AliRICHParam::fgErrGeom[4][330];                        //
Double_t AliRICHParam::fgErrLoc[4][330];                         //Chromatic, Geometric and Localization array to parametrize SigmaCerenkov
Double_t AliRICHParam::fgMass[5]              ={0.00051,0.10566,0.13957,0.49360,0.93828};  


Double_t AliRICHParam::fEckovMin=5.5e-9; //GeV
Double_t AliRICHParam::fEckovMax=8.5e-9; //GeV
TF1 AliRICHParam::fgAbsC6F14("RabsC4F14","6512.39*(x<=7.75e-9)+(x>7.75e-9)*0.039/(-0.166+0.063e9*x-8.01e7*x^2+3.39e5*x^3)" ,fEckovMin,fEckovMax); 
TF1 AliRICHParam::fgAbsSiO2 ("RabsSiO2" ,"333"                                                                             ,fEckovMin,fEckovMax); 
TF1 AliRICHParam::fgAbsCH4  ("RabsCH4"  ,"6512.39*(x<=7.75e-9)+(x>7.75e-9)*0.039/(-0.166+0.063e9*x-8.01e7*x^2+3.39e5*x^3)" ,fEckovMin,fEckovMax);
TF1 AliRICHParam::fgAbsAir  ("RabsAir"  ,"500"                                                                             ,fEckovMin,fEckovMax);  //len ???

TF1 AliRICHParam::fgIdxAir  ("RidxAir"  ,"1+1e-8*(8342.13 + 2406030/(130-(1.23984e-9/x)^2)+15597/(38.9-(1.23984e-9/x)^2))" ,fEckovMin,fEckovMax);  //???
TF1 AliRICHParam::fgIdxSiO2 ("RidxSiO2" ,"sqrt(1+46.411/(10.666*10.666-x*x*1e18)+228.71/(18.125*18.125-x*x*1e18))"         ,fEckovMin,fEckovMax);  //TDR p.35
TF1 AliRICHParam::fgIdxCH4  ("RidxCH4"  ,"1+0.12489e-6/(2.62e-4 - (1239.84e-9/x)^-2)"                                      ,fEckovMin,fEckovMax);  //Olav preprint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliRICHParam::AliRICHParam():TNamed("RichParam","default version") 
{
// Here all the intitializition is taken place when AliRICHParam::Instance() is invoked for the first time.
// In particulare, matrices to be used for LORS<->MARS trasnformations are initialized from TGeo structure.    
// Note that TGeoManager should be already initialized from geometry.root file  
  for(Int_t iCh=0;iCh<kNchambers;iCh++) fMatrix[iCh]=(TGeoHMatrix*)gGeoManager->GetVolume("ALIC")->GetNode(Form("RICH_%i",iCh+1))->GetMatrix();
  CdbRead(0,0);
  fgInstance=this; 
}//ctor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Float_t AliRICHParam::AbsCH4(Float_t eV)
{
// Evaluate the absorbtion lenght of CH4 for a photon of energy eV in electron-volts
  const Float_t kLoschmidt=2.686763e19;                                      // LOSCHMIDT NUMBER IN CM-3
  const Float_t kPressure=750.0;          //mm of Hg
  const Float_t kTemperature=283.0;       //K (10 grad C)                               
  const Float_t kPn=kPressure/760.;
  const Float_t kTn=kTemperature/273.16;
  const Float_t kC0=-1.655279e-1;
  const Float_t kC1= 6.307392e-2;
  const Float_t kC2=-8.011441e-3;
  const Float_t kC3= 3.392126e-4;
    		
  Float_t crossSection=0;                        
  if (eV<7.75) 
    crossSection=0.06e-22;
  else                 //------ METHANE CROSS SECTION cm-2 ASTROPH. J. 214, L47 (1978)                                               
    crossSection=(kC0+kC1*eV+kC2*eV*eV+kC3*eV*eV*eV)*1.e-18;
    
    Float_t density=kLoschmidt*kPn/kTn; //CH4 molecular concentration (cm^-3)
    return 1.0/(density*crossSection);
}//AbsoCH4()
//__________________________________________________________________________________________________sss
void AliRICHParam::CdbRead(Int_t run,Int_t version)
{
// This methode read all the calibration information and initialise corresponding fields for requested run number
// Arguments: run - run number for which to retrieve calibration
//            version- version number   
//   Returns: none      

  AliCDBEntry *pEntry=AliCDBManager::Instance()->Get("RICH/RICHConfig/RefIdxC6F14",run,0,version); //try to get from common local storage  
  if(pEntry){
    fIdxC6F14=(TF2*)pEntry->GetObject(); delete pEntry;
  }else{
    AliWarning("No valid calibarion, the hardcoded will be used!");
    fIdxC6F14=new TF2("RidxC4F14","sqrt(1+0.554*(1239.84e-9/x)^2/((1239.84e-9/x)^2-5796)-0.0005*(y-20))",5.5e-9,8.5e-9,0,50); //DiMauro mail
    fIdxC6F14->SetUniqueID(20);//T=20 deg C
  }
}//CdbRead()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHParam::Print(Option_t* opt) const
{
// print some usefull (hopefully) info on some internal guts of RICH parametrisation 
  Printf("Pads in chamber (%3i,%3i) in sector (%2i,%2i) pad size (%4.2f,%4.2f)",NpadsX(),NpadsY(),NpadsXsec(),NpadsYsec(),PadSizeX(),PadSizeY());
  Printf("Resolve clusters %i sagita %i",IsResolveClusters(),IsWireSag()); 
  
  for(Int_t i=0;i<kNchambers;i++) fMatrix[i]->Print(opt);
}//Print()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHParam::TestSeg()
{
// Provides a set of pictures to test segementation currently in use.    
// Arguments: none
//   Returns: none    
  new TCanvas("pads","PC segmentation - pads display",700,600);
  gPad->Range(-5,-5,PcSizeX()+5,PcSizeY()+15);
  TVector p(2);   TVector2 c;    TVector2 b;   //current: pad, pad center, pad boundary
// list of corners:
  Double_t x0=0,x1=SecSizeX(),x2=SecSizeX()+DeadZone()                                                      ,x3=PcSizeX();  
  Double_t y0=0,y1=SecSizeY(),y2=SecSizeY()+DeadZone(),y3=2*SecSizeY()+DeadZone(),y4=PcSizeY()-SecSizeY(),y5=PcSizeY();  
  DrawSectors();
//header 
  TLatex t;
  t.SetTextSize(0.02); t.SetTextColor(kBlack); t.SetTextAlign(11);
  t.DrawLatex(0,PcSizeY()+10,Form("IP in front of this page. pad size %.2fx%.2fcm   dead zone %.2fcm",PadSizeX(),PadSizeY(),DeadZone()));
  t.DrawLatex(0,PcSizeY()+ 5,Form("Pc  %.2fx%.2f cm %ix%i pads               Sec %.2fx%.2f cm %ix%i pads",
                                         PcSizeX()     , PcSizeY()     , NpadsX()    , NpadsY()                                 ,
                                         SecSizeX() , SecSizeY() , NpadsXsec() , NpadsYsec()                              ));
//sectors  
  t.SetTextSize(0.015); t.SetTextColor(kRed); t.SetTextAlign(22);
  c=Pad2Loc( 40, 24); t.DrawText(c.X(),c.Y(),Form("sec 1 (%.2f,%.2f)",c.X(),c.Y()  ));  
  c=Pad2Loc( 40, 75); t.DrawText(c.X(),c.Y(),Form("sec 3 (%.2f,%.2f)",c.X(),c.Y()  ));  
  c=Pad2Loc( 40,121); t.DrawText(c.X(),c.Y(),Form("sec 5 (%.2f,%.2f)",c.X(),c.Y()  ));  
  c=Pad2Loc(120, 24); t.DrawText(c.X(),c.Y(),Form("sec 2 (%.2f,%.2f)",c.X(),c.Y()  ));  
  c=Pad2Loc(120, 75); t.DrawText(c.X(),c.Y(),Form("sec 4 (%.2f,%.2f)",c.X(),c.Y()  ));  
  c=Pad2Loc(120,121); t.DrawText(c.X(),c.Y(),Form("sec 6 (%.2f,%.2f)",c.X(),c.Y()  ));  
//coners  
  t.SetTextSize(0.015); t.SetTextColor(kBlue);

  b.Set(x0,y0);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(11);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x0,y1);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(13);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x0,y2);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(11);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x0,y3);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(13);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x0,y4);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(11);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x0,y5);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(13);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  
  b.Set(x1,y0);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(31);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x1,y1);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(33);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x1,y2);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(31);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x1,y3);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(33);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x1,y4);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(31);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x1,y5);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(33);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  
  b.Set(x2,y0);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(11);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x2,y1);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(13);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x2,y2);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(11);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x2,y3);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(13);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x2,y4);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(11);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x2,y5);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(13);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  
  b.Set(x3,y0);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(31);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x3,y1);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(33);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x3,y2);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(31);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x3,y3);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(33);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x3,y4);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(31);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
  b.Set(x3,y5);p=Loc2Pad(b);c=Pad2Loc(p);t.SetTextAlign(33);t.DrawText(c.X(),c.Y(),Form("(%.2f,%.2f)-(%.0f,%.0f)-(%.2f,%.2f)",b.X(),b.Y(),p(0),p(1),c.X(),c.Y()));
}//TestSeg()
//__________________________________________________________________________________________________
void AliRICHParam::TestResp()
{
// Provides a set of plot to check the response parametrisation currently in use.  
// Arguments: none
//   Returns: none    
  TCanvas *pC=new TCanvas("c","Amplification test",900,800);
  pC->Divide(1,2);
  
  
  const Int_t kNpoints=8;
  THStack *pStackPhot=new THStack("StackPhot","photons");
  THStack *pStackMip =new THStack("StackMip","mips");
  TLegend *pLeg=new TLegend(0.6,0.2,0.9,0.5,"legend");    
  TH1F *apHphot[kNpoints];
  TH1F *apHmip[kNpoints];
  
  Double_t starty=0;
  Double_t deltay=AliRICHParam::SecSizeY()/kNpoints;
  
  for(int i=0;i<kNpoints;i++){
    apHphot[i]=new TH1F(Form("hphot%i",i),"Qdc for Photon;QDC;Counts",500,0,500); apHphot[i]->SetLineColor(i);pStackPhot->Add(apHphot[i]);                 
    apHmip[i] =new TH1F(Form("hmip%i",i),"Qdc for Mip;QDC;Counts",4000,0,4000);   apHmip[i]->SetLineColor(i);pStackMip->Add(apHmip[i]);                 
    
    pLeg->AddEntry(apHphot[i],Form("@(10,%5.2f->%5.2f)",starty+i*deltay,starty+i*deltay-SecSizeY()/2));
  }
        
  
  TVector2 x2(0,0);  
  for(Int_t i=0;i<10000;i++){//events loop
    for(int j=0;j<kNpoints;j++){
      x2.Set(10,starty+j*deltay);
      apHphot[j]->Fill(TotQdc(x2,0));
      apHmip[j]->Fill(TotQdc(x2,gRandom->Landau(600,150)*1e-9));
    }
  }
  
  pC->cd(1);  pStackMip->Draw("nostack");
  pC->cd(2);  pStackPhot->Draw("nostack"); pLeg->Draw();
}//TestResp()
//__________________________________________________________________________________________________
void AliRICHParam::TestTrans()
{
// Tests transformation methods
// Arguments: none
//   Returns: none    
  
  AliRICHParam *pParam=AliRICHParam::Instance();
  Int_t iNpointsX=50,iNpointsY=50;  
  new TCanvas("trasform","Test LORS-MARS transform"); TLatex t; t.SetTextSize(0.02);
  
  TView *pView=new TView(1);  pView->SetRange(-400,-400,-400,400,400,400);
  DrawAxis();  
  for(Int_t iCham=1;iCham<=7;iCham++){//chamber loop
    AliRICHHelix helix(2.5,Norm(iCham).Theta()*TMath::RadToDeg(),Norm(iCham).Phi()*TMath::RadToDeg());
    helix.RichIntersect(AliRICHParam::Instance());
    TPolyMarker3D *pChamber=new TPolyMarker3D(iNpointsX*iNpointsY);
    Int_t i=0;
    for(Double_t x=0;x<PcSizeX();x+=PcSizeX()/iNpointsX)
      for(Double_t y=0;y<PcSizeY();y+=PcSizeY()/iNpointsY){//step loop
        TVector3 v3=pParam->Lors2Mars(iCham,x,y,kPc); TVector2 v2=pParam->Mars2Lors(iCham,v3,kPc);//LORS->MARS->LORS
        Double_t dx=v2.X()-x , dy=v2.Y()-y; 
        if(dx>0.000001 || dy>0.000001) Printf("Problem in MARS<->LORS transformations dx=%f dy=%f!!!",dx,dy);
        pChamber->SetPoint(i++,v3.X(),v3.Y(),v3.Z());//Pc plane point in MARS
      }//step loop
    pChamber->SetMarkerSize(1);
    pChamber->SetMarkerColor(iCham);
    pChamber->Draw();  
    helix.Draw();
    t.SetNDC();t.SetTextColor(iCham); t.DrawText(0.1,iCham*0.1,Form("Chamber %i",iCham));    
  }//chambers loop   
}//TestTrans()
//__________________________________________________________________________________________________
void AliRICHParam::DrawAxis()
{
// This utility methode draws axis on geometry scene  
// Arguments: none
//   Returns: none    
  Double_t x[6]={0,0,0,300,0,0};  Double_t y[6]={0,0,0,0,300,0};  Double_t z[6]={0,0,0,0,0,300};  
  TPolyLine3D *pXaxis=new TPolyLine3D(2,x);pXaxis->SetLineColor(kRed);   pXaxis->Draw();
  TPolyLine3D *pYaxis=new TPolyLine3D(2,y);pYaxis->SetLineColor(kGreen); pYaxis->Draw();
  TPolyLine3D *pZaxis=new TPolyLine3D(2,z);pZaxis->SetLineColor(kBlue);  pZaxis->Draw();  
}
//__________________________________________________________________________________________________
void AliRICHParam::DrawSectors() 
{ 
// Utility methode draws RICH chamber sectors on event display.
// Arguments: none
//   Returns: none      
  Double_t xLeft[5]  = {0,0,SecSizeX(),SecSizeX(),0};
  Double_t xRight[5] = {SecSizeX()+DeadZone(),SecSizeX()+DeadZone(),PcSizeX(),PcSizeX(),SecSizeX()+DeadZone()};
  
  Double_t yDown[5]   = {0,SecSizeY(),SecSizeY(),0,0};
  Double_t yCenter[5] = {  SecSizeY()+DeadZone(),2*SecSizeY()+DeadZone(),2*SecSizeY()+DeadZone(),
                           SecSizeY()+DeadZone(),SecSizeY()+DeadZone()};  
  Double_t yUp[5]     = {2*SecSizeY()+2*DeadZone(),PcSizeY(),PcSizeY(),2*SecSizeY()+2*DeadZone(),2*SecSizeY()+2*DeadZone()};
  
  TPolyLine *sec1 = new TPolyLine(5,xLeft ,yDown);    sec1->SetLineColor(21);  sec1->Draw();
  TPolyLine *sec2 = new TPolyLine(5,xRight,yDown);    sec2->SetLineColor(21);  sec2->Draw();
  TPolyLine *sec3 = new TPolyLine(5,xLeft ,yCenter);  sec3->SetLineColor(21);  sec3->Draw();
  TPolyLine *sec4 = new TPolyLine(5,xRight,yCenter);  sec4->SetLineColor(21);  sec4->Draw();
  TPolyLine *sec5 = new TPolyLine(5,xLeft, yUp);      sec5->SetLineColor(21);  sec5->Draw();
  TPolyLine *sec6 = new TPolyLine(5,xRight,yUp);      sec6->SetLineColor(21);  sec6->Draw();
}//DrawSectors()
//__________________________________________________________________________________________________
void AliRICHParam::ReadErrFiles()
{
// Read the three files corresponding to Chrom,Geom and Loc They are parameters of a polynomial of 6th order...  ????????? go to CDB?
// Arguments: none
//   Returns: none    
  
  static Bool_t count = kFALSE;
  
  Float_t c0,c1,c2,c3,c;
  Float_t g0,g1,g2,g3,g;
  Float_t l0,l1,l2,l3,l;
  
  FILE *pChromErr, *pGeomErr, *pLocErr;  

  if(!count) {
     AliInfoGeneral("ReadErrFiles","reading RICH error parameters...");
     pChromErr = fopen(Form("%s/RICH/RICHConfig/SigmaChromErr.txt",gSystem->Getenv("ALICE_ROOT")),"r");
     pGeomErr  = fopen(Form("%s/RICH/RICHConfig/SigmaGeomErr.txt",gSystem->Getenv("ALICE_ROOT")),"r");
     pLocErr   = fopen(Form("%s/RICH/RICHConfig/SigmaLocErr.txt",gSystem->Getenv("ALICE_ROOT")),"r");
     if(!pChromErr||!pGeomErr||!pLocErr) {AliErrorGeneral("ReadErrFiles"," RICH ERROR READING Parameter FILES: can't open files!!! ");return;}
     for(Int_t i=0;i<330;i++) {
       fscanf(pChromErr,"%f%f%f%f%f\n",&c0,&c1,&c2,&c3,&c);
       fscanf(pGeomErr,"%f%f%f%f%f\n",&g0,&g1,&g2,&g3,&g);
       fscanf(pLocErr,"%f%f%f%f%f\n",&l0,&l1,&l2,&l3,&l);
       fgErrChrom[0][i] = c0;
       fgErrChrom[1][i] = c1;
       fgErrChrom[2][i] = c2;
       fgErrChrom[3][i] = c3;	
       fgErrGeom[0][i] = g0;
       fgErrGeom[1][i] = g1;
       fgErrGeom[2][i] = g2;
       fgErrGeom[3][i] = g3;	
       fgErrLoc[0][i] = l0;
       fgErrLoc[1][i] = l1;
       fgErrLoc[2][i] = l2;
       fgErrLoc[3][i] = l3;	
     }
     AliInfoGeneral("ReadErrFiles","DONE successfully!");
     fclose(pChromErr);
     fclose(pGeomErr);
     fclose(pLocErr);
  }
  count = kTRUE;
}//ReadErrFiles()
//__________________________________________________________________________________________________
TVector3 AliRICHParam::SigmaSinglePhoton(Int_t partID, Double_t mom, Double_t theta, Double_t phi)

{
// Find sigma for single photon. It returns the thrree different errors. If you want
// to have the error---> TVector3.Mag()
// partID = 0,1,2,3,4 ---> e,mu,pi,k,p in agreement with AliPID
  TVector3 v(-999,-999,-999);
  Double_t pmom;

  ReadErrFiles();
  Double_t mass = fgMass[partID];
  Double_t massRef = fgMass[4]; // all the files are calculated for protons...so mass ref is proton mass
  pmom = mom*massRef/mass; // normalized momentum respect to proton...
  if(pmom>PmodMax()) pmom = PmodMax();
  Double_t oneOverRefIndex = 1/IdxC6F14(EckovMean());
  Double_t pmin = mass*oneOverRefIndex/TMath::Sqrt(1-oneOverRefIndex*oneOverRefIndex);
  if(pmom<pmin) return v;
  Double_t Theta = theta*TMath::RadToDeg();
  if(phi<0) phi+=TMath::TwoPi();
  Double_t Phi = phi*TMath::RadToDeg();
  v.SetX(Interpolate(fgErrChrom,pmom,Theta,Phi));
  v.SetY(Interpolate(fgErrGeom,pmom,Theta,Phi));
  v.SetZ(Interpolate(fgErrLoc,pmom,Theta,Phi));
//  v*=1.5; // take into account bigger errors due to multiplicity...to change in future

  return v;
}//SigmaSinglePhoton
//__________________________________________________________________________________________________
TVector3 AliRICHParam::SigmaSinglePhoton(Double_t thetaCer, Double_t theta, Double_t phi)

{
// Find sigma for single photon. It returns the thrree different errors. If you want
// to have the error---> TVector3.Mag()
// partID = 0,1,2,3,4 ---> e,mu,pi,k,p in agreement with AliPID
  TVector3 v(-999,-999,-999);
  Double_t pmom;

  ReadErrFiles();
  Double_t massRef = fgMass[4]; // all the files are calculated for protons...so mass ref is proton mass
  Double_t beta=1./(IdxC6F14(EckovMean())*TMath::Cos(thetaCer));
  if(beta>=1) {
    pmom=6.5; // above physical limi the error is calculated at the saturation...
  } else {
    Double_t gamma=1./TMath::Sqrt(1-beta*beta);
    pmom = beta*gamma*massRef; // normalized momentum respect to proton...
  }
  if(pmom>PmodMax()) pmom = PmodMax();
  Double_t oneOverRefIndex = 1/IdxC6F14(EckovMean());
  Double_t pmin = massRef*oneOverRefIndex/TMath::Sqrt(1-oneOverRefIndex*oneOverRefIndex);
  if(pmom<pmin) return v;
  Double_t Theta = theta*TMath::RadToDeg();
  if(phi<0) phi+=TMath::TwoPi();
  Double_t Phi = phi*TMath::RadToDeg();
  v.SetX(Interpolate(fgErrChrom,pmom,Theta,Phi));
  v.SetY(Interpolate(fgErrGeom,pmom,Theta,Phi));
  v.SetZ(Interpolate(fgErrLoc,pmom,Theta,Phi));
//  v*=1.5; // take into account bigger errors due to multiplicity...to change in future

  return v;
}//SigmaSinglePhoton
//__________________________________________________________________________________________________
Double_t AliRICHParam::Interpolate(Double_t par[4][330], Double_t x, Double_t y, Double_t phi)
  
{
  static Double_t amin = 1.15; static Double_t astep  = 0.2;
  static Double_t bmin = 0; static Double_t bstep  = 1;
  
  Double_t Phi = (phi - 180)/300.;
  
  Double_t Sigma[30][11];
  
  for(Int_t j=0;j<11;j++) { for(Int_t i=0;i<30;i++) {
    Sigma[i][j] = par[0][j+11*i] + par[1][j+11*i]*Phi*Phi + par[2][j+11*i]*TMath::Power(Phi,4) + par[3][j+11*i]*TMath::Power(Phi,6);
    }
  }

  Int_t i=0;Int_t j=0;
  
  i = (Int_t)((x-amin)/astep);
  j = (Int_t)((y-bmin)/bstep);
  Double_t ai = amin+i*astep;
  Double_t ai1 = ai+astep;
  Double_t bj = bmin+j*bstep;
  Double_t bj1 = bj+bstep;
  Double_t t = (x-ai)/(ai1-ai);
  Double_t gj = (1-t)*Sigma[i][j]+t*Sigma[i+1][j];
  Double_t gj1 = (1-t)*Sigma[i][j+1]+t*Sigma[i+1][j+1];
  Double_t u = (y-bj)/(bj1-bj);
  return (1-u)*gj+u*gj1;
}//Interpolate
//__________________________________________________________________________________________________
TVector3 AliRICHParam::ForwardTracing(TVector3 entranceTrackPoint, TVector3 vectorTrack, Double_t thetaC, Double_t phiC)
{
// Trace a single Ckov photon from a given emission point up to photocathode taking into account ref indexes of materials it travereses
  TVector3 vBad(-999,-999,-999);
  TVector3 nPlane(0,0,1);
  Double_t planeZposition = 0.5*RadThick();
  TVector3 planePoint(0,0,0.5*RadThick()); //this is plane parallel to window which contains emission point 
  TVector3 emissionPoint = PlaneIntersect(vectorTrack,entranceTrackPoint,nPlane,planePoint);
  Double_t thetaout,phiout;
  AnglesInDRS(vectorTrack.Theta(),vectorTrack.Phi(),thetaC,phiC,thetaout,phiout);
  TVector3 vectorPhotonInC6F14;  
  vectorPhotonInC6F14.SetMagThetaPhi(1,thetaout,phiout);
  planeZposition=RadThick();
  planePoint.SetXYZ(0,0,planeZposition);
  TVector3 entranceToSiO2Point = PlaneIntersect(vectorPhotonInC6F14,emissionPoint,nPlane,planePoint);

  Double_t photonEn = EckovMean();
  Double_t angleInSiO2 = SnellAngle(IdxC6F14(EckovMean()),IdxSiO2(EckovMean()),vectorPhotonInC6F14.Theta());if(angleInSiO2<0) return vBad;
  TVector3 vectorPhotonInSiO2;
  vectorPhotonInSiO2.SetMagThetaPhi(1,angleInSiO2,phiout);
//  planeZposition+=AliRICHParam::SiO2Thickness();
  planeZposition+=WinThick();
  planePoint.SetXYZ(0,0,planeZposition);
  TVector3 entranceToCH4 = PlaneIntersect(vectorPhotonInSiO2,entranceToSiO2Point,nPlane,planePoint);
//  entranceToCH4.Dump();

  //  Double_t angleInCH4 = SnellAngle(AliRICHParam::IndOfRefSiO2(6.755),AliRICHParam::IndOfRefCH4,angleInSiO2);
  Double_t angleInCH4 = SnellAngle(IdxSiO2(photonEn),IdxCH4(photonEn),vectorPhotonInSiO2.Theta());if(angleInCH4<0) return vBad;
  TVector3 vectorPhotonInCH4;
  vectorPhotonInCH4.SetMagThetaPhi(1,angleInCH4,phiout);
//  planeZposition+=AliRICHParam::GapProx();
  planeZposition+=Pc2Win();
  planePoint.SetXYZ(0,0,planeZposition);
  TVector3 impactToPC = PlaneIntersect(vectorPhotonInCH4,entranceToCH4,nPlane,planePoint);
//  impactToPC.Dump();
  return impactToPC;
}//FowardTracing
//__________________________________________________________________________________________________
TVector3 AliRICHParam::PlaneIntersect(const TVector3 &lineDir,const TVector3 &linePoint,const TVector3 &planeNorm,const TVector3 &planePoint)
{
// Finds an intersection point between a line and plane.
// Arguments:  lineDir,linePoint    - vector along the line and any point of the line
//             planeNorm,planePoint - vector normal to the plane and any point of the plane
//   Returns:  point of intersection if any
  if(planeNorm*lineDir==0) return   TVector3(-999,-999,-999);
  TVector3 diff=planePoint-linePoint;
  Double_t sint=(planeNorm*diff)/(planeNorm*lineDir);
  return linePoint+sint*lineDir;
}//PlaneIntersect
//__________________________________________________________________________________________________ 
Double_t AliRICHParam::SnellAngle(Float_t n1, Float_t n2, Float_t theta1)
{
// Compute the angle of refraction out of Snell law
// Arguments: n1 - ref idx of first substance  
//            n2 - ref idx of second substance  
//            n1 - photon impact angle in the first substance i.e. angle between the photon direction and vector normal to the surface (radians)
//   Returns: photon refraction angle, i.e. angle in the second substance (radians)
  Double_t sinref=(n1/n2)*TMath::Sin(theta1);
  if(sinref>1.)    return -999;  
  else             return TMath::ASin(sinref);
}//SnellAngle
//__________________________________________________________________________________________________
void AliRICHParam::AnglesInDRS(Double_t trackTheta,Double_t trackPhi,Double_t thetaCerenkov,Double_t phiCerenkov,Double_t &tout,Double_t &pout)
{
// Setup the rotation matrix of the track...

  TRotation mtheta;
  TRotation mphi;
  TRotation minv;
  TRotation mrot;
  
  mtheta.RotateY(trackTheta);
  mphi.RotateZ(trackPhi);
  
  mrot = mphi * mtheta;
    //  minv = mrot.Inverse();

  TVector3 photonInRadiator(1,1,1);

  photonInRadiator.SetTheta(thetaCerenkov);
  photonInRadiator.SetPhi(phiCerenkov);
  photonInRadiator = mrot * photonInRadiator;
  tout=photonInRadiator.Theta();
  pout=photonInRadiator.Phi();
}//AnglesInDRS
/*
void DrawRing()
{

  //  Float_t xGraph[1000],yGraph[1000];

  Float_t type;
  Float_t MassOfParticle;
  Float_t beta;
  Float_t nfreon;

  Float_t ThetaCerenkov;

  Float_t Xtoentr = GetEntranceX();
  Float_t Ytoentr = GetEntranceY();

  Float_t pmod = GetTrackMomentum();
  Float_t TrackTheta = GetTrackTheta();
  Float_t TrackPhi = GetTrackPhi();

  SetPhotonEnergy(AliRICHParam::MeanCkovEnergy());
  SetFreonRefractiveIndex();

  SetEmissionPoint(RadiatorWidth/2.);

  ThetaCerenkov = GetThetaCerenkov();
  FindBetaFromTheta(ThetaCerenkov);
  nfreon = GetFreonRefractiveIndex();
  
  Int_t nPoints = 100;

  Int_t nPointsToDraw = 0;
  for(Int_t i=0;i<nPoints;i++)
    {
      Float_t phpad = 2*TMath::Pi()*i/nPoints;
      SetThetaPhotonInTRS(thetacer);
      SetPhiPhotonInTRS(phpad);
      FindPhotonAnglesInDRS();
      Float_t Radius = FromEmissionToCathode();
      if (Radius == 999.) continue;
      xGraph[nPointsToDraw] = GetXPointOnCathode() + GetShiftX();
      yGraph[nPointsToDraw] = GetYPointOnCathode() + GetShiftY();
      nPointsToDraw++;
    }
  gra = new TGraph(nPointsToDraw,xGraph,yGraph);
  gra->Draw("AC"); 
}
//__________________________________________________________________________________________________
*/
void AliRICHParam::TestHit2SDigs(Double_t x,Double_t y,Double_t e,Bool_t isNew)
{
//Test  hit->sdigits procedures
//Arguments: isNew - if true use new (abs pad) procedure else use old one (TVector)
//  Returns: none
  TClonesArray *pSDigLst=new TClonesArray("AliRICHDigit");
  Int_t iQtot=-1;
  if(isNew){
    iQtot=Hit2SDigs(10101,e,pSDigLst);        //new technique
  }else{
    iQtot=Hit2SDigs(TVector2(x,y),e,pSDigLst);//old technique
  }
  pSDigLst->Print();
  Double_t dQsum=0;
  for(Int_t i=0;i<pSDigLst->GetEntriesFast();i++)
    dQsum+=((AliRICHDigit*)pSDigLst->At(i))->Qdc();
  Printf("Qtot=%i Qsum=%.2f ",iQtot,dQsum);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHParam::Stack(Int_t evt,Int_t tid)
{
// Prints some usefull info from stack
// Arguments: evt - event number. if not -1 print info only for that event
//            tid - track id. if not -1 then print it and all it's mothers if any   
//   Returns: mother tid of the given tid if any
  AliRunLoader *pAL=AliRunLoader::Open(); 
  if(pAL->LoadHeader()) return -1;
  if(pAL->LoadKinematics()) return -1;
  
  Int_t mtid=-1;
  Int_t iNevt=pAL->GetNumberOfEvents();     Printf("This session contains %i event(s)",iNevt);
  
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){//events loop
    if(evt!=-1 && evt!=iEvt) continue; //in case one needs to print the requested event, ignore all others
    pAL->GetEvent(iEvt);    
    AliStack *pStack=pAL->Stack();  
    if(tid==-1){                        //print all tids for this event
      for(Int_t i=0;i<pStack->GetNtrack();i++) pStack->Particle(i)->Print();
      Printf("totally %i tracks including %i primaries for event %i out of %i event(s)",pStack->GetNtrack(),pStack->GetNprimary(),iEvt,iNevt);
    }else{                              //print only this tid and it;s mothers
      if(tid<0 || tid>pStack->GetNtrack()) {Printf("Wrong tid, valid tid range for event %i is 0-%i",iEvt,pStack->GetNtrack());break;}
      TParticle *pTrack=pStack->Particle(tid); mtid=pTrack->GetFirstMother();
      TString str=pTrack->GetName();
      while((tid=pTrack->GetFirstMother()) >= 0){
        pTrack=pStack->Particle(tid);
        str+=" from ";str+=pTrack->GetName();
      } 
      Printf("%s",str.Data());       
    }//if(tid==-1)      
  }//events loop
  pAL->UnloadHeader();  pAL->UnloadKinematics();
  return mtid;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHParam::StackCount(Int_t pid,Int_t evt)
{
// Counts total number of particles of given sort (including secondary) for a given event
  AliRunLoader *pAL=AliRunLoader::Open(); 
  pAL->GetEvent(evt);    
  if(pAL->LoadHeader()) return 0;
  if(pAL->LoadKinematics()) return 0;
  AliStack *pStack=pAL->Stack();
  
  Int_t iCnt=0;
  for(Int_t i=0;i<pStack->GetNtrack();i++) if(pStack->Particle(i)->GetPdgCode()==pid) iCnt++;
  
  pAL->UnloadHeader();  pAL->UnloadKinematics();
  return iCnt;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
