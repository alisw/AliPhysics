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
#include "AliRICHParam.h"
#include "AliRICHChamber.h"
#include "AliRICHDisplFast.h"
#include <TCanvas.h>
#include <TLatex.h>
#include <THStack.h>
#include <TLegend.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

ClassImp(AliRICHParam)
Bool_t   AliRICHParam::fgIsWireSag            =kTRUE;   //take ware sagita into account?
Bool_t   AliRICHParam::fgIsResolveClusters    =kTRUE;   //do cluster resolving?
Bool_t   AliRICHParam::fgIsRadioSrc           =kFALSE;  //put radioactive source instead of radiators?
Bool_t   AliRICHParam::fgIsAerogel            =kFALSE;  //special aerogel configuration
Bool_t   AliRICHParam::fgIsTestBeam           =kFALSE;  //special test beam configuration

Int_t    AliRICHParam::fgHV[kNsectors]        ={2050,2050,2050,2050,2050,2050};
Int_t    AliRICHParam::fgNsigmaTh             =4;
Float_t  AliRICHParam::fgSigmaThMean          =1.132; //QDC 
Float_t  AliRICHParam::fgSigmaThSpread        =0.035; //     

//__________________________________________________________________________________________________
void AliRICHParam::Print(Option_t*)
{
  AliInfo(Form("Pads in chamber (%3i,%3i) in sector (%2i,%2i)",NpadsX(),NpadsY(),NpadsXsec(),NpadsYsec()));
  AliInfo(Form("Resolve clusters %i sagita %i Radio source %i Aerogel %i TestBeam %i",
                IsResolveClusters(),IsWireSag(),IsRadioSrc(),IsAerogel(),IsTestBeam())); 
  fpChambers->Print();
}//Print()
//__________________________________________________________________________________________________
void AliRICHParam::CreateChambers()
{
//Create all RICH Chambers on each call. Previous chambers deleted.
  if(fpChambers) delete fpChambers;
  if(fgIsTestBeam){ 
    fpChambers=new TObjArray(1);//test beam configuration 1 chamber
    fpChambers->AddAt(new AliRICHChamber(0),0);  
  }else{ 
    fpChambers=new TObjArray(kNchambers);//normal configuration 7 chambers
    for(int iChamberN=0;iChamberN<kNchambers;iChamberN++)  fpChambers->AddAt(new AliRICHChamber(iChamberN+1),iChamberN);  
  }
  fpChambers->SetOwner();
}//CreateChambers()
//__________________________________________________________________________________________________
Float_t AliRICHParam::AbsCH4(Float_t eV)
{
//Evaluate the absorbtion lenght of CH4 for a photon of energy eV in electron-volts
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
//__________________________________________________________________________________________________
void AliRICHParam::TestSeg()
{
  
  new TCanvas("name","PC segmentation");
  gPad->Range(-20,-20,200,150);
  AliRICHDisplFast::DrawSectors();
  
  TLatex t; t.SetTextSize(0.02);
  t.DrawText(0,140,"View from interaction point");
  t.DrawLatex(PcSizeX()+10,120,Form("Pc  %6.2fx%6.2fcm %3ix%3ipads",PcSizeX()    ,PcSizeY(),    NpadsX()   ,NpadsY()));
  t.DrawLatex(PcSizeX()+10,115,Form("Sec %6.2fx%5.2fcm %3ix%2ipads",SectorSizeX(),SectorSizeY(),NpadsXsec(),NpadsYsec()));
  t.DrawLatex(PcSizeX()+10,110,Form("Pad %6.2fx%4.2fcm DeadZone %6.2fcm",PadSizeX()   ,PadSizeY(),DeadZone()));
  
  TVector2 v2;
  t.SetTextAlign(12);
  v2=AliRICHParam::Pad2Loc( 1,24);  t.DrawText(v2.X(),v2.Y(),"sec 1");  
  v2=AliRICHParam::Pad2Loc(81,24);  t.DrawText(v2.X(),v2.Y(),"sec 2");  
  v2=AliRICHParam::Pad2Loc( 1,70);  t.DrawText(v2.X(),v2.Y(),"sec 3");  
  v2=AliRICHParam::Pad2Loc(81,70);  t.DrawText(v2.X(),v2.Y(),"sec 4");  
  v2=AliRICHParam::Pad2Loc( 1,120); t.DrawText(v2.X(),v2.Y(),"sec 5");  
  v2=AliRICHParam::Pad2Loc(81,120); t.DrawText(v2.X(),v2.Y(),"sec 6");  
  
//  TGaxis *pAx=new TGaxis(0,0,140,  0,0,140,510,"-="); pAx->SetTitle("x, cm"); pAx->SetTextSize(0.05); pAx->Draw();
//  TGaxis *pAy=new TGaxis(0,0,  0,140,0,140,510,"-="); pAy->SetTitle("y, cm"); pAy->SetTextSize(0.05); pAy->Draw();
   
  
  t.SetTextColor(kBlue);  

  Double_t margin=5;
  Double_t x0=0; Double_t x1=SectorSizeX(); Double_t x2=SectorSizeX()+DeadZone(); Double_t x3=PcSizeX();
  Double_t y0=0; Double_t y1=SectorSizeY(); Double_t y2=SectorSizeY()+DeadZone(); 
  Double_t y3=2*SectorSizeY()+DeadZone(); Double_t y4=PcSizeY()-SectorSizeY();
  Double_t y5=PcSizeY();
   
//write pad numbers along x  
  t.SetTextAlign(13); t.DrawText(x0,y0,"1");  
  t.SetTextAlign(33); t.DrawText(x1,y0,"80");
  t.SetTextAlign(13); t.DrawText(x2,y0,"81");   
  t.SetTextAlign(33); t.DrawText(x3,y0,"160");
//write pad numbers along y    
  t.SetTextAlign(31); t.DrawText(x0,y0,"1");  
  t.SetTextAlign(33); t.DrawText(x0,y1,"48"); 
  t.SetTextAlign(31); t.DrawText(x0,y2,"49");
  t.SetTextAlign(33); t.DrawText(x0,y3,"96"); 
  t.SetTextAlign(31); t.DrawText(x0,y4,"97");
  t.SetTextAlign(33); t.DrawText(x0,y5,"144");        
  
   
//positions along x   
  t.SetTextColor(kGreen);
  t.SetTextAlign(11);t.DrawText(x0,y0-margin,Form("%5.2f",x0));   
  t.SetTextAlign(31);t.DrawText(x1,y0-margin,Form("%5.2f",x1));   
  t.SetTextAlign(11);t.DrawText(x2,y0-margin,Form("%5.2f",x2));   
  t.SetTextAlign(31);t.DrawText(x3,y0-margin,Form("%5.2f",x3));   
//positions along y
  t.SetTextAlign(31);t.DrawText(x0-margin,y0,Form("%5.2f",y0));   
  t.SetTextAlign(33);t.DrawText(x0-margin,y1,Form("%5.2f",y1));   
  t.SetTextAlign(31);t.DrawText(x0-margin,y2,Form("%5.2f",y2));   
  t.SetTextAlign(33);t.DrawText(x0-margin,y3,Form("%5.2f",y3));   
  t.SetTextAlign(31);t.DrawText(x0-margin,y4,Form("%5.2f",y4));   
  t.SetTextAlign(33);t.DrawText(x0-margin,y5,Form("%5.2f",y5));   
//coners      
  t.SetTextColor(kRed);
  t.SetTextSize(0.01);
  TVector pad(2);
//sector 1  
  v2=Pad2Loc(Loc2Pad(TVector2(x0,y0)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x0,y1)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x1,y1)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x1,y0)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
//secr 3
  v2=Pad2Loc(Loc2Pad(TVector2(x0,y2)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x0,y3)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x1,y3)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x1,y2)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
//secr 5   
  v2=Pad2Loc(Loc2Pad(TVector2(x0,y4)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x0,y5)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x1,y5)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x1,y4)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));

  v2=Pad2Loc(Loc2Pad(TVector2(x2,y4)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x2,y5)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x3,y5)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x3,y4)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));

  v2=Pad2Loc(Loc2Pad(TVector2(x2,y2)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x2,y3)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x3,y3)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x3,y2)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));

  v2=Pad2Loc(Loc2Pad(TVector2(x2,y0)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x2,y1)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x3,y1)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=Pad2Loc(Loc2Pad(TVector2(x3,y0)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
}//TestSeg()
//__________________________________________________________________________________________________
void AliRICHParam::TestResp()
{
//test the response set of methodes  
  TCanvas *pC=new TCanvas("c","Amplification test",900,800);
  pC->Divide(1,2);
  
  
  const Int_t nPoints=8;
  THStack *pStackPhot=new THStack("StackPhot","photons");
  THStack *pStackMip =new THStack("StackMip","mips");
  TLegend *pLeg=new TLegend(0.6,0.2,0.9,0.5,"legend");    
  TH1F *apHphot[nPoints];
  TH1F *apHmip[nPoints];
  
  Double_t starty=0;
  Double_t deltay=AliRICHParam::SectorSizeY()/nPoints;
  
  for(int i=0;i<nPoints;i++){
    apHphot[i]=new TH1F(Form("hphot%i",i),"Qdc for Photon;QDC;Counts",500,0,500); apHphot[i]->SetLineColor(i);pStackPhot->Add(apHphot[i]);                 
    apHmip[i] =new TH1F(Form("hmip%i",i),"Qdc for Mip;QDC;Counts",4000,0,4000);   apHmip[i]->SetLineColor(i);pStackMip->Add(apHmip[i]);                 
    
    pLeg->AddEntry(apHphot[i],Form("@(10,%5.2f->%5.2f)",starty+i*deltay,starty+i*deltay-SectorSizeY()/2));
  }
        
  
  TVector2 x2(0,0);  
  for(Int_t i=0;i<10000;i++){//events loop
    for(int j=0;j<nPoints;j++){
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
//test the set of transformation methods
  new TCanvas("trasform","Test LRS-MRS transform");
  TLatex t; t.SetTextSize(0.02);
  
  TView *pView=new TView(1);
  pView->SetRange(-600,-600,-600,600,600,600);
  DrawAxis();  
//Draw PC for all chambers by trasfering Pc plane using Pc2Mrs methode  
  Int_t iNpointsX=50,iNpointsY=50;  
  for(Int_t iChamberN=1;iChamberN<=7;iChamberN++){//chamber loop
    TPolyMarker3D *pChamber=new TPolyMarker3D(iNpointsX*iNpointsY);
    Int_t i=0;
    for(Double_t x=0;x<PcSizeX();x+=PcSizeX()/iNpointsX)
      for(Double_t y=0;y<PcSizeY();y+=PcSizeY()/iNpointsY){//step loop
        TVector3 v3=C(iChamberN)->Pc2Mrs(TVector2(x,y));//from regular grid of local PC points to MRS presentation
        pChamber->SetPoint(i++,v3.X(),v3.Y(),v3.Z());//Pc plane poing in MRS
      }//step loop
    pChamber->SetMarkerSize(1);
    pChamber->SetMarkerColor(iChamberN);
    pChamber->Draw();  
    t.SetNDC();t.SetTextColor(iChamberN); t.DrawText(0.1,iChamberN*0.1,Form("Chamber %i",iChamberN));    
  }//chamber loop   
//  gPad->GetView()->RotateView(94,45);
}//TestTrans()
//__________________________________________________________________________________________________
void AliRICHParam::DrawAxis()
{
  Double_t X[6]={0,0,0,300,0,0};  Double_t Y[6]={0,0,0,0,300,0};  Double_t Z[6]={0,0,0,0,0,300};  
  TPolyLine3D *pXaxis=new TPolyLine3D(2,X);pXaxis->SetLineColor(kRed);   pXaxis->Draw();
  TPolyLine3D *pYaxis=new TPolyLine3D(2,Y);pYaxis->SetLineColor(kGreen); pYaxis->Draw();
  TPolyLine3D *pZaxis=new TPolyLine3D(2,Z);pZaxis->SetLineColor(kBlue);  pZaxis->Draw();  
}
