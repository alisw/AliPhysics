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
#include "AliESD.h"
#include "AliRICHChamber.h"
#include <TCanvas.h>
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
Double_t AliRICHParam::fgErrChrom[4][330];                       //
Double_t AliRICHParam::fgErrGeom[4][330];                        //
Double_t AliRICHParam::fgErrLoc[4][330];                         //Chromatic, Geometric and Localization array to parametrize SigmaCerenkov
  

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
//__________________________________________________________________________________________________
void AliRICHParam::ReadErrFiles()
{
// Read the three files corresponding to Chrom,Geom and Loc
// They are parameters of a polynomial of 6th order...
  
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
  
  TVector3 v(-999,-999,-999);
  Double_t pmom;

  ReadErrFiles();
  Double_t mass = AliPID::ParticleMass(partID);
  Double_t massRef = AliPID::ParticleMass(AliPID::kProton); // all the files are calculated for protons...so mass ref is proton mass
  pmom = mom*massRef/mass; // normalized momentum respect to proton...
  if(pmom>6.5) pmom = 6.5;
  Double_t oneOverRefIndex = 1/RefIdxC6F14(6.755);
  Double_t pmin = mass*oneOverRefIndex/TMath::Sqrt(1-oneOverRefIndex*oneOverRefIndex);
  if(pmom<pmin) return v;
  v.SetX(Interpolate(fgErrChrom,pmom,theta,phi));
  v.SetY(Interpolate(fgErrGeom,pmom,theta,phi));
  v.SetZ(Interpolate(fgErrLoc,pmom,theta,phi));

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
//
  TVector3 vBad(-999,-999,-999);
  TVector3 nPlane(0,0,1);
  Double_t planeZposition = 0.5*Zfreon();
  TVector3 planePoint(0,0,planeZposition);
  TVector3 emissionPoint = PlaneIntersect(vectorTrack,entranceTrackPoint,nPlane,planePoint);
//  emissionPoint.Dump();
  Double_t thetaout,phiout;
  AnglesInDRS(vectorTrack.Theta(),vectorTrack.Phi(),thetaC,phiC,thetaout,phiout);
//  cout << "thetaout "<<thetaout << " phiout " << phiout << endl;
  TVector3 vectorPhotonInC6F14;  
  vectorPhotonInC6F14.SetMagThetaPhi(1,thetaout,phiout);
//  vectorPhotonInC6F14.Dump();
//  planeZposition=AliRICHParam::C6F14Thickness();
  planeZposition=Zfreon();
  planePoint.SetXYZ(0,0,planeZposition);
  TVector3 entranceToSiO2Point = PlaneIntersect(vectorPhotonInC6F14,emissionPoint,nPlane,planePoint);
//  entranceToSiO2Point.Dump();

  Double_t photonEn = MeanCkovEnergy();
  Double_t angleInSiO2 = SnellAngle(RefIdxC6F14(photonEn),RefIdxSiO2(photonEn),vectorPhotonInC6F14.Theta());if(angleInSiO2<0) return vBad;
  TVector3 vectorPhotonInSiO2;
  vectorPhotonInSiO2.SetMagThetaPhi(1,angleInSiO2,phiout);
//  planeZposition+=AliRICHParam::SiO2Thickness();
  planeZposition+=Zwin();
  planePoint.SetXYZ(0,0,planeZposition);
  TVector3 entranceToCH4 = PlaneIntersect(vectorPhotonInSiO2,entranceToSiO2Point,nPlane,planePoint);
//  entranceToCH4.Dump();

  //  Double_t angleInCH4 = SnellAngle(AliRICHParam::IndOfRefSiO2(6.755),AliRICHParam::IndOfRefCH4,angleInSiO2);
  Double_t angleInCH4 = SnellAngle(RefIdxSiO2(photonEn),RefIdxCH4(photonEn),vectorPhotonInSiO2.Theta());if(angleInCH4<0) return vBad;
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
TVector3 AliRICHParam::PlaneIntersect(TVector3 vstart,TVector3 p0,TVector3 n,TVector3 v0)
{
//
  TVector3 parallel(-999,-999,-999);
  // vstart = given vector
  // p0 = origin of the given vector
  // n = normal to a given plane
  // v0 = point of the given plane
//  cout << " n*vstart = " << n*vstart << endl;
  if(n*vstart==0) return parallel;
  TVector3 diff=v0-p0;
  Double_t sint=(n*diff)/(n*vstart);
  return p0+sint*vstart;
}//PlaneIntersect
//__________________________________________________________________________________________________ 
Double_t AliRICHParam::SnellAngle(Float_t n1, Float_t n2, Float_t theta1)
{
// Snell law
// Compute the Snell angle

  Double_t sinrefractangle;
  Double_t refractangle;

  sinrefractangle = (n1/n2)*sin(theta1);

  if(sinrefractangle>1.) {
    //    cout << " PROBLEMS IN SNELL ANGLE !!!!! " << endl;
    refractangle = -999.;
    return refractangle;
  }

  refractangle = asin(sinrefractangle);
  return refractangle;
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
//__________________________________________________________________________________________________

//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
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
