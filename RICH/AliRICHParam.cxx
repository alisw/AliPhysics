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
#include <TPolyLine.h>

ClassImp(AliRICHParam)
Bool_t   AliRICHParam::fgIsWireSag            =kTRUE;   //take ware sagita into account?
Bool_t   AliRICHParam::fgIsResolveClusters    =kTRUE;   //do cluster resolving?
Bool_t   AliRICHParam::fgIsFeedback           =kTRUE;   //generate feedback photons?
Bool_t   AliRICHParam::fgIsRadioSrc           =kFALSE;  //put radioactive source instead of radiators?
Bool_t   AliRICHParam::fgIsAerogel            =kFALSE;  //special aerogel configuration
Bool_t   AliRICHParam::fgIsTestBeam           =kFALSE;  //special test beam configuration

Int_t    AliRICHParam::fgHV[kNsectors]        ={2050,2050,2050,2050,2050,2050};
Int_t    AliRICHParam::fgNsigmaTh             =4;
Float_t  AliRICHParam::fgSigmaThMean          =1.132; //QDC 
Float_t  AliRICHParam::fgSigmaThSpread        =0.035; //     

//__________________________________________________________________________________________________
void AliRICHParam::Print(Option_t*) const
{
//print some usefull (hopefully) info on some internal guts of RICH parametrisation 
  AliInfo(Form("Pads in chamber (%3i,%3i) in sector (%2i,%2i) pad size (%4.2f,%4.2f)",NpadsX(),NpadsY(),NpadsXsec(),NpadsYsec(),PadSizeX(),PadSizeY()));
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
//Provides a set of pictures to test segementation currently in use.    
  new TCanvas("pads","PC segmentation - pads display",700,600);
  gPad->Range(-5,-5,PcSizeX()+5,PcSizeY()+15);
  TVector p(2);   TVector2 c;    TVector2 b;   //current: pad, pad center, pad boundary
// list of corners:
  Double_t x0=0,x1=SectorSizeX(),x2=SectorSizeX()+DeadZone(),                                                                 x3=PcSizeX();  
  Double_t y0=0,y1=SectorSizeY(),y2=SectorSizeY()+DeadZone(),y3=2*SectorSizeY()+DeadZone(),y4=PcSizeY()-SectorSizeY(),y5=PcSizeY();  
  DrawSectors();
//header 
  TLatex t;
  t.SetTextSize(0.02); t.SetTextColor(kBlack); t.SetTextAlign(11);
  t.DrawLatex(0,PcSizeY()+10,Form("IP in front of this page. pad size %.2fx%.2fcm   dead zone %.2fcm",PadSizeX(),PadSizeY(),DeadZone()));
  t.DrawLatex(0,PcSizeY()+ 5,Form("Pc  %.2fx%.2f cm %ix%i pads               Sec %.2fx%.2f cm %ix%i pads",
                                         PcSizeX()     , PcSizeY()     , NpadsX()    , NpadsY()                                 ,
                                         SectorSizeX() , SectorSizeY() , NpadsXsec() , NpadsYsec()                              ));
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
//Provides a set of plot to check the response parametrisation currently in use.  
  TCanvas *pC=new TCanvas("c","Amplification test",900,800);
  pC->Divide(1,2);
  
  
  const Int_t kNpoints=8;
  THStack *pStackPhot=new THStack("StackPhot","photons");
  THStack *pStackMip =new THStack("StackMip","mips");
  TLegend *pLeg=new TLegend(0.6,0.2,0.9,0.5,"legend");    
  TH1F *apHphot[kNpoints];
  TH1F *apHmip[kNpoints];
  
  Double_t starty=0;
  Double_t deltay=AliRICHParam::SectorSizeY()/kNpoints;
  
  for(int i=0;i<kNpoints;i++){
    apHphot[i]=new TH1F(Form("hphot%i",i),"Qdc for Photon;QDC;Counts",500,0,500); apHphot[i]->SetLineColor(i);pStackPhot->Add(apHphot[i]);                 
    apHmip[i] =new TH1F(Form("hmip%i",i),"Qdc for Mip;QDC;Counts",4000,0,4000);   apHmip[i]->SetLineColor(i);pStackMip->Add(apHmip[i]);                 
    
    pLeg->AddEntry(apHphot[i],Form("@(10,%5.2f->%5.2f)",starty+i*deltay,starty+i*deltay-SectorSizeY()/2));
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
//Provides a set of plots to test transformation methods
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
//Utility: draws axis on geometry scene  
  Double_t x[6]={0,0,0,300,0,0};  Double_t y[6]={0,0,0,0,300,0};  Double_t z[6]={0,0,0,0,0,300};  
  TPolyLine3D *pXaxis=new TPolyLine3D(2,x);pXaxis->SetLineColor(kRed);   pXaxis->Draw();
  TPolyLine3D *pYaxis=new TPolyLine3D(2,y);pYaxis->SetLineColor(kGreen); pYaxis->Draw();
  TPolyLine3D *pZaxis=new TPolyLine3D(2,z);pZaxis->SetLineColor(kBlue);  pZaxis->Draw();  
}
//__________________________________________________________________________________________________
void AliRICHParam::DrawSectors() 
{ 
//Utility: draws RICH chamber sectors on event display.
  Double_t xLeft[5]  = {0,0,SectorSizeX(),SectorSizeX(),0};
  Double_t xRight[5] = {SectorSizeX()+DeadZone(),SectorSizeX()+DeadZone(),PcSizeX(),PcSizeX(),SectorSizeX()+DeadZone()};
  
  Double_t yDown[5]   = {0,SectorSizeY(),SectorSizeY(),0,0};
  Double_t yCenter[5] = {  SectorSizeY()+DeadZone(),2*SectorSizeY()+DeadZone(),2*SectorSizeY()+DeadZone(),
                           SectorSizeY()+DeadZone(),SectorSizeY()+DeadZone()};  
  Double_t yUp[5]     = {2*SectorSizeY()+2*DeadZone(),PcSizeY(),PcSizeY(),2*SectorSizeY()+2*DeadZone(),2*SectorSizeY()+2*DeadZone()};
  
  TPolyLine *sec1 = new TPolyLine(5,xLeft ,yDown);    sec1->SetLineColor(21);  sec1->Draw();
  TPolyLine *sec2 = new TPolyLine(5,xRight,yDown);    sec2->SetLineColor(21);  sec2->Draw();
  TPolyLine *sec3 = new TPolyLine(5,xLeft ,yCenter);  sec3->SetLineColor(21);  sec3->Draw();
  TPolyLine *sec4 = new TPolyLine(5,xRight,yCenter);  sec4->SetLineColor(21);  sec4->Draw();
  TPolyLine *sec5 = new TPolyLine(5,xLeft, yUp);      sec5->SetLineColor(21);  sec5->Draw();
  TPolyLine *sec6 = new TPolyLine(5,xRight,yUp);      sec6->SetLineColor(21);  sec6->Draw();
}//DrawSectors()
