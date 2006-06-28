// **************************************************************************
// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// *                                                                        *
// * Author: The ALICE Off-line Project.                                    *
// * Contributors are mentioned in the code where appropriate.              *
// *                                                                        *
// * Permission to use, copy, modify and distribute this software and its   *
// * documentation strictly for non-commercial purposes is hereby granted   *
// * without fee, provided that the above copyright notice appears in all   *
// * copies and that both the copyright notice and this permission notice   *
// * appear in the supporting documentation. The authors make no claims     *
// * about the suitability of this software for any purpose. It is          *
// * provided "as is" without express or implied warranty.                  *
// **************************************************************************


#include "AliRICHv1.h"     //class header
#include "AliRICHParam.h"
#include <TParticle.h>     //Hits2SDigits()
#include <TRandom.h> 
#include <TVirtualMC.h>    //StepManager() for gMC
#include <TPDGCode.h>      //StepHistory() 
#include <AliStack.h>      //StepManager(),Hits2SDigits()
#include <AliLoader.h>        //Hits2SDigits()
#include <AliRunLoader.h>     //Hits2SDigits()
#include <AliConst.h>
#include <AliPDG.h>
#include <AliMC.h>            //StepManager()      
#include <AliRawDataHeader.h> //Digits2Raw()
#include <AliDAQ.h>           //Digits2Raw()
#include <AliRun.h>           //CreateMaterials()    
#include <AliMagF.h>          //CreateMaterials()
#include <TGeoManager.h>      //CreateGeometry()
#include <TMultiGraph.h>      //Optics() 
#include <TGraph.h>           //Optics() 
#include <TLegend.h>          //Optics() 
#include <TCanvas.h>          //Optics() 
#include <TF2.h>              //CreateMaterials()
#include <AliCDBManager.h>    //CreateMaterials()
#include <AliCDBEntry.h>      //CreateMaterials()
 
ClassImp(AliRICHv1)    

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::CreateMaterials()
{
// Definition of available RICH materials  
// Arguments: none
//   Returns: none    
  AliDebug(1,"Start v1 RICH.");
  
  const Int_t kNbins=30;       //number of photon energy points
  
  Float_t aAbsC6F14[kNbins]={//New values from A.DiMauro 28.10.03 total 30
    32701.4219, 17996.1141, 10039.7281, 1799.1230, 1799.1231, 1799.1231, 1241.4091, 179.0987, 179.0986, 179.0987,
      179.0987,   118.9800,    39.5058,   23.7244,   11.1283,    7.1573,    3.6249,   2.1236,   0.7362,   0.5348,
        0.3387,     0.3074,     0.3050,    0.0001,    0.0001,    0.0001,    0.0001,   0.0001,   0.0001,   0.0001};    
    
  Float_t aAbsSiO2[kNbins]={//New values from A.DiMauro 28.10.03 total 30
       34.4338, 30.5424, 30.2584, 31.4928, 31.7868, 17.8397, 9.3410, 6.4492, 6.1128, 5.8128,
        5.5589,  5.2877,  5.0162,  4.7999,  4.5734,  4.2135, 3.7471, 2.6033, 1.5223, 0.9658,
        0.4242,  0.2500,  0.1426,  0.0863,  0.0793,  0.0724, 0.0655, 0.0587, 0.0001, 0.0001};
    
  Float_t aQeCsI[kNbins] = {//New values from A.DiMauro 28.10.03 total 31 the last one cut to provide 30
                            0.0002, 0.0006, 0.0007, 0.0010, 0.0049, 0.0073, 0.0104, 0.0519, 0.0936, 0.1299,
                            0.1560, 0.1768, 0.1872, 0.1976, 0.2142, 0.2288, 0.2434, 0.2599, 0.2673, 0.2808,
                            0.2859, 0.2954, 0.3016, 0.3120, 0.3172, 0.3224, 0.3266, 0.3328, 0.3359, 0.3390}; //0.3431};
                              
  Float_t aQeCsIold[kNbins]={//previous values 26 in total added 0.0001 to the 30
    0.0002, 0.0006, 0.0007, 0.0050, 0.0075, 0.0101, 0.0243, 0.0405, 0.0689, 0.1053, 
    0.1215, 0.1417, 0.1579, 0.1620, 0.1661, 0.1677, 0.1743, 0.1768, 0.1793, 0.1826,
    0.1859, 0.1876, 0.1892, 0.1909, 0.2075, 0.2158, 0.0001, 0.0001, 0.0001, 0.0001 };      
                            
//         radiator            window             gas               metal 
  Float_t aIdxC6F14[kNbins] , aIdxSiO2[kNbins] , aIdxCH4[kNbins] , aIdx0[kNbins]   , aIdx1[kNbins] ; 
  Float_t                                        aAbsCH4[kNbins] , aAbsMet[kNbins]                 ;
  Float_t                                                                            aQe1[kNbins]  ;    //QE for all but PC
  Float_t aEckov[kNbins];  //Ckov energy in GeV
                            
  
  Double_t eCkovMin=5.5e-9,eCkovMax=8.5e-9; //in GeV
  TF2 idxC6F14("RidxC4F14","sqrt(1+0.554*(1239.84e-9/x)^2/((1239.84e-9/x)^2-5796)-0.0005*(y-20))"   ,eCkovMin,eCkovMax,0,50); //DiMauro mail temp 0-50 degrees C
  TF1 idxSiO2( "RidxSiO2" ,"sqrt(1+46.411/(10.666*10.666-x*x*1e18)+228.71/(18.125*18.125-x*x*1e18))",eCkovMin,eCkovMax);  //TDR p.35
  
  
  for(Int_t i=0;i<kNbins;i++){
    aEckov     [i] =eCkovMin+0.1e-9*i;//Ckov energy in GeV
    
    aIdxC6F14  [i] =idxC6F14.Eval(aEckov[i],20);      //Simulation for 20 degress C       
    aIdxSiO2   [i] =idxSiO2 .Eval(aEckov[i]);
    aIdxCH4    [i] =AliRICHParam::IdxCH4   (aEckov[i]); aAbsCH4 [i]  =AliRICHParam::AbsCH4     (aEckov[i]); 
    aAbsMet    [i] =0.0001;                              //metal has absorption probability
    aIdx0      [i] =0;                                   //metal ref idx must be 0 in order to reflect photon
    aIdx1      [i] =1;                                   //metal ref idx must be 1 in order to apply photon to QE conversion 
    aQe1       [i] =1;                                   //QE for all other materials except for PC must be 1.
  }
  
//data from PDG booklet 2002     density [gr/cm^3] rad len [cm] abs len [cm]    
  Float_t   aAir[4]={12,14,16,36}    ,   zAir[4]={6,7,8,18} ,   wAir[4]={0.000124,0.755267,0.231781,0.012827} , dAir=0.00120479; Int_t nAir=4;//mixture 0.9999999
  Float_t aC6F14[2]={ 12.01 , 18.99} , zC6F14[2]={ 6 , 9}   , wC6F14[2]={6 , 14} , dC6F14=1.68    ; Int_t nC6F14=-2;
  Float_t  aSiO2[2]={ 28.09 , 15.99} ,  zSiO2[2]={14 , 8}   ,  wSiO2[2]={1 ,  2} ,  dSiO2=2.64    ; Int_t  nSiO2=-2; 
  Float_t   aCH4[2]={ 12.01 ,  1.01} ,   zCH4[2]={ 6 , 1}   ,   wCH4[2]={1 ,  4} ,   dCH4=7.17e-4 ; Int_t   nCH4=-2; 
  Float_t   aCsI[2]={132.90 ,126.90} ,   zCsI[2]={55 ,53}   ,   wCsI[2]={1 ,  1} ,   dCsI=0.1     ; Int_t   nCsI=-2; 
  Float_t     aRoha= 12.01 ,               zRoha=  6 ,                              dRoha=  0.10 , radRoha= 18.80 , absRoha=  86.3/dRoha; //special material- quazi carbon
  Float_t       aCu= 63.55 ,                 zCu= 29 ,                                dCu=  8.96 ,   radCu=  1.43 ,   absCu= 134.9/dCu  ;
  Float_t        aW=183.84 ,                  zW= 74 ,                                 dW= 19.30 ,    radW=  0.35 ,    absW= 185.0/dW   ;
  Float_t       aAl= 26.98 ,                 zAl= 13 ,                                dAl=  2.70 ,   radAl=  8.90 ,   absAl= 106.4/dAl  ;
    
  Int_t   matId=0;                           //tmp material id number
  Int_t   unsens =  0, sens=1;               //sensitive or unsensitive medium
  Int_t   itgfld = gAlice->Field()->Integ(); //type of field intergration 0 no field -1 user in guswim 1 Runge Kutta 2 helix 3 const field along z
  Float_t maxfld = gAlice->Field()->Max();   //max field value
  Float_t tmaxfd = -10.0;                    //max deflection angle due to magnetic field in one step
  Float_t deemax = - 0.2;                    //max fractional energy loss in one step   
  Float_t stemax = - 0.1;                    //mas step allowed [cm]
  Float_t epsil  =   0.001;                  //abs tracking precision [cm]   
  Float_t stmin  = - 0.001;                  //min step size [cm] in continius process transport, negative value: choose it automatically
  AliMixture(++matId,"Air"  ,aAir  ,zAir  ,dAir  ,nAir  ,wAir  ); AliMedium(kAir  ,"Air"  ,matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  AliMixture(++matId,"C6F14",aC6F14,zC6F14,dC6F14,nC6F14,wC6F14); AliMedium(kC6F14,"C6F14",matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);      
  AliMixture(++matId,"SiO2" ,aSiO2 ,zSiO2 ,dSiO2 ,nSiO2 ,wSiO2 ); AliMedium(kSiO2 ,"SiO2" ,matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);    
  AliMixture(++matId,"CH4"  ,aCH4  ,zCH4  ,dCH4  ,nCH4  ,wCH4  ); AliMedium(kCH4  ,"CH4"  ,matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);  
  AliMixture(++matId,"CsI"  ,aCsI  ,zCsI  ,dCsI  ,nCsI  ,wCsI  ); AliMedium(kCsI  ,"CsI"  ,matId,   sens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);//sensitive
  
  AliMaterial(++matId,"Roha",aRoha,zRoha,dRoha,radRoha,absRoha);  AliMedium(kRoha,"Roha", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  AliMaterial(++matId,"Cu"  ,aCu  ,zCu  ,dCu  ,radCu  ,absCu  );  AliMedium(kCu  ,"Cu"  , matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  AliMaterial(++matId,"W"   ,aW   ,zW   ,dW   ,radW   ,absW   );  AliMedium(kW   ,"W"   , matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  AliMaterial(++matId,"Al"  ,aAl  ,zAl  ,dAl  ,radAl  ,absAl  );  AliMedium(kAl  ,"Al"  , matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  gMC->SetCerenkov((*fIdtmed)[kC6F14]    , kNbins, aEckov, aAbsC6F14  , aQe1  , aIdxC6F14 );    
  gMC->SetCerenkov((*fIdtmed)[kSiO2]     , kNbins, aEckov, aAbsSiO2   , aQe1  , aIdxSiO2  );    
  gMC->SetCerenkov((*fIdtmed)[kCH4]      , kNbins, aEckov, aAbsCH4    , aQe1  , aIdxCH4   );    
  gMC->SetCerenkov((*fIdtmed)[kCu]       , kNbins, aEckov, aAbsMet    , aQe1  , aIdx0     );    
  gMC->SetCerenkov((*fIdtmed)[kW]        , kNbins, aEckov, aAbsMet    , aQe1  , aIdx0     ); //n=0 means reflect photons       
  gMC->SetCerenkov((*fIdtmed)[kCsI]      , kNbins, aEckov, aAbsMet    , aQeCsI, aIdx1     ); //n=1 means convert photons    
  gMC->SetCerenkov((*fIdtmed)[kAl]       , kNbins, aEckov, aAbsMet    , aQe1  , aIdx0     );    
  
  AliDebug(1,"Stop v1 RICH.");
  TString ttl=GetTitle();
  if(!ttl.Contains("ShowOptics")) return; //do not plot optical curves
  {
  const Double_t kWidth=0.25,kHeight=0.2;  
  const Int_t kC6F14Marker=24 , kC6F14Color=kRed;  
  const Int_t kCH4Marker  =25 , kCH4Color  =kGreen;  
  const Int_t kSiO2M      =26 , kSiO2Color =kBlue;  
  const Int_t kCsIMarker  = 2 , kCsIColor  =kMagenta;  
  
//Ref index           
  TGraph *pIdxC6F14=new TGraph(kNbins,aEckov,aIdxC6F14);pIdxC6F14->SetMarkerStyle(kC6F14Marker);pIdxC6F14->SetMarkerColor(kC6F14Color);
  TGraph *pIdxSiO2 =new TGraph(kNbins,aEckov,aIdxSiO2); pIdxSiO2 ->SetMarkerStyle(kSiO2M)      ;pIdxSiO2 ->SetMarkerColor(kSiO2Color);  
  TGraph *pIdxCH4  =new TGraph(kNbins,aEckov,aIdxCH4);  pIdxCH4  ->SetMarkerStyle(kCH4Marker)  ;pIdxCH4  ->SetMarkerColor(kCH4Color);
  TMultiGraph *pIdxMG=new TMultiGraph("refidx","Ref index;E_{#check{C}} [GeV]");  TLegend *pIdxLe=new TLegend(0.5,0.21,0.5+kWidth,0.21+kHeight);
  pIdxMG->Add(pIdxC6F14); pIdxLe->AddEntry(pIdxC6F14,"C6F14"  ,"p");            
  pIdxMG->Add(pIdxSiO2) ; pIdxLe->AddEntry(pIdxSiO2 ,"SiO2"   ,"p");          
  pIdxMG->Add(pIdxCH4)  ; pIdxLe->AddEntry(pIdxCH4  ,"CH4"    ,"p");          
//Absorbtion
  TGraph *pAbsC6F14=new TGraph(kNbins,aEckov,aAbsC6F14);pAbsC6F14->SetMarkerStyle(kC6F14Marker); pAbsC6F14->SetMarkerColor(kC6F14Color);
  TGraph *pAbsSiO2 =new TGraph(kNbins,aEckov,aAbsSiO2) ;pAbsSiO2 ->SetMarkerStyle(kSiO2M)      ; pAbsSiO2 ->SetMarkerColor(kSiO2Color);
  TGraph *pAbsCH4  =new TGraph(kNbins,aEckov,aAbsCH4)  ;pAbsCH4  ->SetMarkerStyle(kCH4Marker)  ; pAbsCH4  ->SetMarkerColor(kCH4Color);
  
  TMultiGraph *pAbsMG=new TMultiGraph("abs","Absorption [cm];E_{#check{C}} [GeV]");  TLegend *pAbsLe=new TLegend(0.2,0.15,0.2+kWidth,0.15+kHeight);
  pAbsMG->Add(pAbsC6F14);      pAbsLe->AddEntry(pAbsC6F14,  "C6F14"    ,"p"); 
  pAbsMG->Add(pAbsSiO2);       pAbsLe->AddEntry(pAbsSiO2 ,  "SiO2"     ,"p"); 
  pAbsMG->Add(pAbsCH4);        pAbsLe->AddEntry(pAbsCH4  ,  "CH4"      ,"p"); 
//QE new and old
  TGraph *pQeCsI   =new TGraph(kNbins,aEckov,aQeCsI);    pQeCsI   ->SetMarkerStyle(kCsIMarker);  pQeCsI   ->SetMarkerColor(kCsIColor);
  TGraph *pQeCsIold=new TGraph(kNbins,aEckov,aQeCsIold); pQeCsIold->SetMarkerStyle(kC6F14Marker);pQeCsIold->SetMarkerColor(kC6F14Color);    
  TMultiGraph *pCompMG=new TMultiGraph("qe","QE;E_{#check{C}} [GeV]");  TLegend *pCompLe=new TLegend(0.2,0.6,0.2+kWidth,0.6+kHeight);
  pCompMG->Add(pQeCsI);       pCompLe->AddEntry(pQeCsI,    "QE new 30.10.03", "p");  
  pCompMG->Add(pQeCsIold);    pCompLe->AddEntry(pQeCsIold, "QE old 01.01.02", "p");
//transmission  
  Float_t aTrC6F14[kNbins],aTrSiO2[kNbins],aTrCH4[kNbins],aTrTotal[kNbins];
  for(Int_t i=0;i<kNbins;i++){//calculate probability for photon to survive during transversing a volume of material with absorption length  
    aTrC6F14[i] =TMath::Exp(-AliRICHParam::RadThick() / (aAbsC6F14[i]+0.0001)); //radiator
    aTrSiO2[i]  =TMath::Exp(-AliRICHParam::WinThick() / (aAbsSiO2[i] +0.0001)); //window
    aTrCH4[i]   =TMath::Exp(-AliRICHParam::Pc2Win()   / (aAbsCH4[i]  +0.0001)); //from window to PC   
    aTrTotal[i] =aTrC6F14[i]*aTrSiO2[i]*aTrCH4[i]*aQeCsI[i];
  }
  TGraph *pTrC6F14=new TGraph(kNbins,aEckov,aTrC6F14)  ;pTrC6F14->SetMarkerStyle(kC6F14Marker);pTrC6F14->SetMarkerColor(kC6F14Color);  
  TGraph *pTrSiO2 =new TGraph(kNbins,aEckov,aTrSiO2)   ;pTrSiO2 ->SetMarkerStyle(kSiO2M)      ;pTrSiO2 ->SetMarkerColor(kSiO2Color);  
  TGraph *pTrCH4  =new TGraph(kNbins,aEckov,aTrCH4)    ;pTrCH4  ->SetMarkerStyle(kCH4Marker)  ;pTrCH4  ->SetMarkerColor(kCH4Color);   
  TGraph *pTrTotal=new TGraph(kNbins,aEckov,aTrTotal)  ;pTrTotal->SetMarkerStyle(30)          ;pTrTotal->SetMarkerColor(kYellow);  
  TMultiGraph *pTrMG=new TMultiGraph("trans","Transmission;E_{#check{C}} [GeV]");  TLegend *pTrLe=new TLegend(0.2,0.4,0.2+kWidth,0.4+kHeight);
  pTrMG->Add(pQeCsI);       pTrLe->AddEntry(pQeCsI,   "CsI QE", "p");  
  pTrMG->Add(pTrC6F14);     pTrLe->AddEntry(pTrC6F14, "C6F14" , "p");  
  pTrMG->Add(pTrSiO2);      pTrLe->AddEntry(pTrSiO2,  "SiO2"  , "p");          
  pTrMG->Add(pTrCH4);       pTrLe->AddEntry(pTrCH4,   "CH4"   , "p");          
  pTrMG->Add(pTrTotal);     pTrLe->AddEntry(pTrTotal, "total" , "p");          
  
  TCanvas *pC=new TCanvas("c1","RICH optics to check",1100,900);  pC->Divide(2,2);           
  pC->cd(1);                    pIdxMG ->Draw("AP");  pIdxLe ->Draw();      //ref idx       
  pC->cd(2);  gPad->SetLogy();  pAbsMG ->Draw("AP");  pAbsLe ->Draw();      //absorption
  pC->cd(3);                    pCompMG->Draw("AP");  pCompLe->Draw();      //QE      
  pC->cd(4);                    pTrMG  ->Draw("AP");  pTrLe  ->Draw();      //transmission
  }
}//void AliRICH::CreateMaterials()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::CreateGeometry()
{
//Creates detailed geometry simulation (currently GEANT volumes tree)         
  AliDebug(1,"Start main.");
  if(!gMC->IsRootGeometrySupported()) return;     
  
  Double_t cm=1,mm=0.1*cm,mkm=0.001*mm,dx,dy,dz;//default is cm
  
  TGeoVolume *pRich=gGeoManager->MakeBox("RICH",gGeoManager->GetMedium("RICH_CH4"),dx=(6*mm+1681*mm+6*mm)/2,  //main RICH volume
                                                                                   dy=(6*mm+1466*mm+6*mm)/2,
                                                                                   dz=(80*mm+40*mm)*2/2);     //x,y taken from 2033P1  z from p84 TDR  
  const Double_t kAngHor=19.5; //  horizontal angle between chambers  19.5 grad
  const Double_t kAngVer=20;   //  vertical angle between chambers    20   grad     
  const Double_t kAngCom=30;   //  common RICH rotation with respect to x axis  30   grad     
  const Double_t trans[3]={490,0,0}; //center of the chamber is on window-gap surface
  for(Int_t iCh=1;iCh<=7;iCh++){//place 7 chambers
    TGeoHMatrix *pMatrix=new TGeoHMatrix;
    pMatrix->RotateY(90);           //rotate around y since initial position is in XY plane -> now in YZ plane
    pMatrix->SetTranslation(trans); //now plane in YZ is shifted along x 
    switch(iCh){
      case 1:                pMatrix->RotateY(kAngHor);  pMatrix->RotateZ(-kAngVer);  break; //right and down 
      case 2:                                            pMatrix->RotateZ(-kAngVer);  break; //down              
      case 3:                pMatrix->RotateY(kAngHor);                               break; //right 
      case 4:                                                                         break; //no rotation
      case 5:                pMatrix->RotateY(-kAngHor);                              break; //left   
      case 6:                                            pMatrix->RotateZ(kAngVer);   break; //up
      case 7:                pMatrix->RotateY(-kAngHor); pMatrix->RotateZ(kAngVer);   break; //left and up 
    }
    pMatrix->RotateZ(kAngCom);     //apply common rotation  in XY plane    
    gGeoManager->GetVolume("ALIC")->AddNode(pRich,iCh,pMatrix);
  }

  Float_t par[3];
  Int_t matrixIdReturn=0; //matrix id returned by AliMatrix
//Pad Panel frame  6 sectors
  par[0]=648*mm/2;par[1]=  411*mm/2;par[2]=40  *mm/2;gMC->Gsvolu("Rppf"     ,"BOX ",(*fIdtmed)[kAl]  ,par,3);//PPF 2001P2 inner size of the slab by 1mm more
  par[0]=181*mm/2;par[1]=89.25*mm/2;par[2]=38.3*mm/2;gMC->Gsvolu("RppfLarge","BOX ",(*fIdtmed)[kAir] ,par,3);//large whole
  par[0]=114*mm/2;par[1]=89.25*mm/2;par[2]=38.3*mm/2;gMC->Gsvolu("RppfSmall","BOX ",(*fIdtmed)[kAir] ,par,3);//small whole
  par[0]=644*mm/2;par[1]=  407*mm/2;par[2]= 1.7*mm/2;gMC->Gsvolu("Rpc"      ,"BOX ",(*fIdtmed)[kCsI] ,par,3);//by 0.2 mm more then actual size (PCB 2006P1)
  
  gMC->Gspos("Rppf",1,"RICH",    -335*mm,      -433*mm,  8*cm+20*mm,  0,"ONLY");//F1 2040P1 z p.84 TDR
  gMC->Gspos("Rppf",2,"RICH",    +335*mm,      -433*mm,  8*cm+20*mm,  0,"ONLY");
  gMC->Gspos("Rppf",3,"RICH",    -335*mm,         0*mm,  8*cm+20*mm,  0,"ONLY");
  gMC->Gspos("Rppf",4,"RICH",    +335*mm,         0*mm,  8*cm+20*mm,  0,"ONLY");
  gMC->Gspos("Rppf",5,"RICH",    -335*mm,      +433*mm,  8*cm+20*mm,  0,"ONLY");
  gMC->Gspos("Rppf",6,"RICH",    +335*mm,      +433*mm,  8*cm+20*mm,  0,"ONLY");  
    gMC->Gspos("Rpc"      ,1,"Rppf",       0*mm,         0*mm,   -19.15*mm,  0,"ONLY");//PPF 2001P2 
    gMC->Gspos("RppfLarge",1,"Rppf",  -224.5*mm,  -151.875*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfLarge",2,"Rppf",  -224.5*mm,  - 50.625*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfLarge",3,"Rppf",  -224.5*mm,  + 50.625*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfLarge",4,"Rppf",  -224.5*mm,  +151.875*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfSmall",1,"Rppf",  - 65.0*mm,  -151.875*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfSmall",2,"Rppf",  - 65.0*mm,  - 50.625*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfSmall",3,"Rppf",  - 65.0*mm,  + 50.625*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfSmall",4,"Rppf",  - 65.0*mm,  +151.875*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfSmall",5,"Rppf",  + 65.0*mm,  -151.875*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfSmall",6,"Rppf",  + 65.0*mm,  - 50.625*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfSmall",7,"Rppf",  + 65.0*mm,  + 50.625*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfSmall",8,"Rppf",  + 65.0*mm,  +151.875*mm,     0.85*mm,  0,"ONLY"); 
    gMC->Gspos("RppfLarge",5,"Rppf",  +224.5*mm,  -151.875*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfLarge",6,"Rppf",  +224.5*mm,  - 50.625*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfLarge",7,"Rppf",  +224.5*mm,  + 50.625*mm,     0.85*mm,  0,"ONLY");
    gMC->Gspos("RppfLarge",8,"Rppf",  +224.5*mm,  +151.875*mm,     0.85*mm,  0,"ONLY");
//Gap - anod wires 6 copies to RICH
  par[0]=648*mm/2;par[1]=  411*mm/2 ;par[2]=4.45*mm/2;gMC->Gsvolu("Rgap","BOX ",(*fIdtmed)[kCH4] ,par,3);//xy as PPF 2001P2 z WP 2099P1
  par[0]=  0*mm  ;par[1]=  20*mkm/2 ;par[2]= 648*mm/2;gMC->Gsvolu("Rano","TUBE",(*fIdtmed)[kW]   ,par,3);//WP 2099P1 z = gap x PPF 2001P2
  AliMatrix(matrixIdReturn,180,0, 90,90, 90,0); //wires along x
  
  gMC->Gspos("Rgap",1,"RICH",    -335*mm,      -433*mm,8*cm-2.225*mm, 0,"ONLY"); //F1 2040P1 z WP 2099P1
  gMC->Gspos("Rgap",2,"RICH",    +335*mm,      -433*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("Rgap",3,"RICH",    -335*mm,         0*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("Rgap",4,"RICH",    +335*mm,         0*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("Rgap",5,"RICH",    -335*mm,      +433*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("Rgap",6,"RICH",    +335*mm,      +433*mm,8*cm-2.225*mm, 0,"ONLY"); 
  for(int i=1;i<=96;i++)
    gMC->Gspos("Rano",i,"Rgap",     0*mm, -411/2*mm+i*4*mm, 0.185*mm, matrixIdReturn,"ONLY"); //WP 2099P1  
//Defines radiators geometry  
  par[0]=1330*mm/2 ;par[1]= 413*mm/2  ;par[2]=  24*mm/2;  gMC->Gsvolu("Rrad"      ,"BOX ",(*fIdtmed)[kC6F14]     ,par,3); // Rad 2011P1
  par[0]=1330*mm/2 ;par[1]= 413*mm/2  ;par[2]=   4*mm/2;  gMC->Gsvolu("RradFront" ,"BOX ",(*fIdtmed)[kRoha]      ,par,3); //front 
  par[0]=1330*mm/2 ;par[1]= 413*mm/2  ;par[2]=   5*mm/2;  gMC->Gsvolu("RradWin"   ,"BOX ",(*fIdtmed)[kSiO2]      ,par,3); //window
  par[0]=1330*mm/2 ;par[1]=   5*mm/2  ;par[2]=  15*mm/2;  gMC->Gsvolu("RradLong"  ,"BOX ",(*fIdtmed)[kRoha]      ,par,3); //long side  
  par[0]=  10*mm/2 ;par[1]= 403*mm/2  ;par[2]=  15*mm/2;  gMC->Gsvolu("RradShort" ,"BOX ",(*fIdtmed)[kRoha]      ,par,3); //short side 
  par[0]=   0      ;par[1]=  10*mm/2  ;par[2]=  15*mm/2;  gMC->Gsvolu("RradSpacer","TUBE",(*fIdtmed)[kSiO2]      ,par,3); //spacer        
    
  gMC->Gspos("Rrad",1,"RICH",   0*mm,-434*mm,   -12*mm,  0,"ONLY"); //3 radiators to RICH
  gMC->Gspos("Rrad",2,"RICH",   0*mm,   0*mm,   -12*mm,  0,"ONLY"); 
  gMC->Gspos("Rrad",3,"RICH",   0*mm,+434*mm,   -12*mm,  0,"ONLY"); 
    gMC->Gspos("RradFront",1,"Rrad",   0*mm,   0*mm, -10.0*mm,  0,"ONLY"); //front cover 
    gMC->Gspos("RradWin"  ,1,"Rrad",   0*mm,   0*mm,   9.5*mm,  0,"ONLY"); //quartz window (back cover)
    gMC->Gspos("RradLong" ,1,"Rrad",   0*mm,-204*mm,  -0.5*mm,  0,"ONLY"); //long side
    gMC->Gspos("RradLong" ,2,"Rrad",   0*mm,+204*mm,  -0.5*mm,  0,"ONLY"); //long side
    gMC->Gspos("RradShort",1,"Rrad",-660*mm,   0*mm,  -0.5*mm,  0,"ONLY"); //short side
    gMC->Gspos("RradShort",2,"Rrad",+660*mm,   0*mm,  -0.5*mm,  0,"ONLY"); //short side 
    for(int i=0;i<3;i++)
      for(int j=0;j<10;j++)
        gMC->Gspos("RradSpacer",10*i+j,"Rrad",-1330*mm/2+116*mm+j*122*mm,(i-1)*105*mm,-0.5*mm,0,"ONLY");//spacers
//Defines SandBox geometry
  par[0]=1419*mm/2 ;par[1]=1378*mm/2;par[2]=50.5*mm/2; gMC->Gsvolu("Rsb"     ,"BOX ",(*fIdtmed)[kAir]  ,par,3);  //2072P1   
  par[0]=1419*mm/2 ;par[1]=1378*mm/2;par[2]= 0.5*mm/2; gMC->Gsvolu("RsbCover","BOX ",(*fIdtmed)[kAl]   ,par,3);  //cover
  par[0]=1359*mm/2 ;par[1]=1318*mm/2;par[2]=49.5*mm/2; gMC->Gsvolu("RsbComb" ,"BOX ",(*fIdtmed)[kRoha] ,par,3);  //honeycomb structure 
  
  gMC->Gspos("Rsb",1,"RICH",   0*mm, 0*mm, -73.75*mm, 0,"ONLY"); //p.84 TDR sandbox to rich
    gMC->Gspos("RsbComb" ,1,"Rsb", 0*mm, 0*mm,      0*mm, 0,"ONLY"); //2072P1 honeycomv to sandbox
    gMC->Gspos("RsbCover",1,"Rsb", 0*mm, 0*mm,    +25*mm, 0,"ONLY"); //cover to sandbox
    gMC->Gspos("RsbCover",2,"Rsb", 0*mm, 0*mm,    -25*mm, 0,"ONLY"); //cover to sandbox
  AliDebug(1,"Stop v1. HMPID option");  
}//CreateGeometry()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::Init()
{
// This methode defines ID for sensitive volumes, i.e. such geometry volumes for which there are if(gMC->CurrentVolID()==XXX) statements in StepManager()
// Arguments: none
//   Returns: none      
  AliDebug(1,"Start v1 HMPID.");    
  fIdRad     = gMC->VolId("Rrad");
  fIdWin     = gMC->VolId("RradWin");
  fIdPc      = gMC->VolId("Rpc");
  fIdAmpGap  = gMC->VolId("Rgap");
  fIdProxGap = gMC->VolId("Rgap");
  AliDebug(1,"Stop v1 HMPID.");    
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliRICHv1::IsLostByFresnel()
{
// Calculate probability for the photon to be lost by Fresnel reflection.
  TLorentzVector p4;
  Double_t mom[3],localMom[3];
  gMC->TrackMomentum(p4);   mom[0]=p4(1);   mom[1]=p4(2);   mom[2]=p4(3);
  localMom[0]=0; localMom[1]=0; localMom[2]=0;
  gMC->Gmtod(mom,localMom,2);
  Double_t localTc    = localMom[0]*localMom[0]+localMom[2]*localMom[2];
  Double_t localTheta = TMath::ATan2(TMath::Sqrt(localTc),localMom[1]);
  Double_t cotheta = TMath::Abs(TMath::Cos(localTheta));
  if(gMC->GetRandom()->Rndm() < Fresnel(p4.E()*1e9,cotheta,1)){
    AliDebug(1,"Photon lost");
    return kTRUE;
  }else
    return kFALSE;
}//IsLostByFresnel()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::GenFee(Int_t iChamber,Float_t eloss)
{
// Generate FeedBack photons for the current particle. To be invoked from StepManager().
// eloss=0 means photon so only pulse height distribution is to be analysed. This one is done in AliRICHParam::TotQdc()
  
  TLorentzVector x4;
  gMC->TrackPosition(x4);  
  TVector2 x2=AliRICHParam::Instance()->Mars2Lors(iChamber,x4.Vect(),AliRICHParam::kPc);//hit position on photocathode plane
  TVector2 xspe=x2;
  Int_t sector=AliRICHParam::Loc2Sec(xspe);  if(sector==-1) return; //hit in dead zone, nothing to produce
  Int_t iTotQdc=AliRICHParam::TotQdc(x2,eloss);
  Int_t iNphotons=gMC->GetRandom()->Poisson(AliRICHParam::AlphaFeedback(iChamber,sector)*iTotQdc);    
  AliDebug(1,Form("N photons=%i",iNphotons));
  Int_t j;
  Float_t cthf, phif, enfp = 0, sthf, e1[3], e2[3], e3[3], vmod, uswop,dir[3], phi,pol[3], mom[4];
//Generate photons
  for(Int_t i=0;i<iNphotons;i++){//feedbacks loop
    Double_t ranf[2];
    gMC->GetRandom()->RndmArray(2,ranf);    //Sample direction
    cthf=ranf[0]*2-1.0;
    if(cthf<0) continue;
    sthf = TMath::Sqrt((1 - cthf) * (1 + cthf));
    phif = ranf[1] * 2 * TMath::Pi();
    
    if(Double_t randomNumber=gMC->GetRandom()->Rndm()<=0.57)
      enfp = 7.5e-9;
    else if(randomNumber<=0.7)
      enfp = 6.4e-9;
    else
      enfp = 7.9e-9;
    

    dir[0] = sthf * TMath::Sin(phif);    dir[1] = cthf;    dir[2] = sthf * TMath::Cos(phif);
    gMC->Gdtom(dir, mom, 2);
    mom[0]*=enfp;    mom[1]*=enfp;    mom[2]*=enfp;
    mom[3] = TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
    
    // Polarisation
    e1[0]=      0;    e1[1]=-dir[2];    e1[2]= dir[1];
    e2[0]=-dir[1];    e2[1]= dir[0];    e2[2]=      0;
    e3[0]= dir[1];    e3[1]=      0;    e3[2]=-dir[0];
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e1[j];
      e1[j]=e3[j];
      e3[j]=uswop;
    }
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e2[j];
      e2[j]=e3[j];
      e3[j]=uswop;
    }
    
    vmod=0;  for(j=0;j<3;j++) vmod+=e1[j]*e1[j];  vmod=TMath::Sqrt(1/vmod);  for(j=0;j<3;j++) e1[j]*=vmod;    
    vmod=0;  for(j=0;j<3;j++) vmod+=e2[j]*e2[j];  vmod=TMath::Sqrt(1/vmod);  for(j=0;j<3;j++) e2[j]*=vmod;
    
    phi = gMC->GetRandom()->Rndm()* 2 * TMath::Pi();
    for(j=0;j<3;j++) pol[j]=e1[j]*TMath::Sin(phi)+e2[j]*TMath::Cos(phi);
    gMC->Gdtom(pol, pol, 2);
    Int_t outputNtracksStored;    
    gAlice->GetMCApp()->PushTrack(1,                             //transport
                     gAlice->GetMCApp()->GetCurrentTrackNumber(),//parent track 
                     kFeedback,                                  //PID
		     mom[0],mom[1],mom[2],mom[3],                //track momentum  
                     x4.X(),x4.Y(),x4.Z(),x4.T(),                //track origin 
                     pol[0],pol[1],pol[2],                       //polarization
		     kPFeedBackPhoton,                           //process ID   
                     outputNtracksStored,                        //on return how many new photons stored on stack
                     1.0);                                       //weight
  }//feedbacks loop
  AliDebug(1,"Stop.");
}//GenerateFeedbacks()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::Hits2SDigits()
{
// Create a list of sdigits corresponding to list of hits. Every hit generates one or more sdigits.
  AliDebug(1,"Start.");
  for(Int_t iEventN=0;iEventN<GetLoader()->GetRunLoader()->GetAliRun()->GetEventsPerRun();iEventN++){//events loop
    GetLoader()->GetRunLoader()->GetEvent(iEventN);//get next event
  
    if(!GetLoader()->TreeH()) GetLoader()->LoadHits();    GetLoader()->GetRunLoader()->LoadHeader(); 
    if(!GetLoader()->GetRunLoader()->TreeK())             GetLoader()->GetRunLoader()->LoadKinematics();//from
    if(!GetLoader()->TreeS()) GetLoader()->MakeTree("S"); MakeBranch("S");//to
          
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      GetLoader()->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){//hits loop 
        AliRICHHit *pHit=(AliRICHHit*)Hits()->At(iHitN);//get current hit                
        TVector2 x2 = AliRICHParam::Instance()->Mars2Lors(pHit->C(),0.5*(pHit->InX3()+pHit->OutX3()),AliRICHParam::kAnod);//hit position in the anod plane
        Int_t iTotQdc=AliRICHParam::TotQdc(x2,pHit->Eloss());//total charge produced by hit, 0 if hit in dead zone
        if(iTotQdc==0) continue;
        //
        //need to quantize the anod....
        TVector padHit=AliRICHParam::Loc2Pad(x2);
        TVector2 padHitXY=AliRICHParam::Pad2Loc(padHit);
        TVector2 anod;
        if((x2.Y()-padHitXY.Y())>0) anod.Set(x2.X(),padHitXY.Y()+AliRICHParam::AnodPitch()/2);
        else anod.Set(x2.X(),padHitXY.Y()-AliRICHParam::AnodPitch()/2);
        //end to quantize anod
        //
        TVector area=AliRICHParam::Loc2Area(anod);//determine affected pads, dead zones analysed inside
        AliDebug(1,Form("hitanod(%6.2f,%6.2f)->area(%3.0f,%3.0f)-(%3.0f,%3.0f) QDC=%4i",anod.X(),anod.Y(),area[0],area[1],area[2],area[3],iTotQdc));
        TVector pad(2);
        for(pad[1]=area[1];pad[1]<=area[3];pad[1]++)//affected pads loop
          for(pad[0]=area[0];pad[0]<=area[2];pad[0]++){                    
            Double_t padQdc=iTotQdc*AliRICHParam::FracQdc(anod,pad);
            AliDebug(1,Form("current pad(%3.0f,%3.0f) with QDC  =%6.2f",pad[0],pad[1],padQdc));
            if(padQdc>0.1) SDigAdd(pHit->C(),pad,padQdc,GetLoader()->GetRunLoader()->Stack()->Particle(pHit->GetTrack())->GetPdgCode(),pHit->GetTrack());
          }//affected pads loop 
      }//hits loop
    }//prims loop
    GetLoader()->TreeS()->Fill();
    GetLoader()->WriteSDigits("OVERWRITE");
    SDigReset();
  }//events loop  
  GetLoader()->UnloadHits(); GetLoader()->GetRunLoader()->UnloadHeader(); GetLoader()->GetRunLoader()->UnloadKinematics();
  GetLoader()->UnloadSDigits();  
  AliDebug(1,"Stop.");
}//Hits2SDigits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::Digits2Raw()
{
//Creates raw data files in DDL format. Invoked by AliSimulation
//loop over events is done outside in AliSimulation
//Arguments: none
//  Returns: none    
  AliDebug(1,"Start.");
  GetLoader()->LoadDigits();
  GetLoader()->TreeD()->GetEntry(0);
  
  ofstream file[AliRICHDigit::kNddls];   //output streams
  Int_t    cnt[AliRICHDigit::kNddls];        //data words counters for DDLs
  AliRawDataHeader header; //empty DDL header
  UInt_t w32=0;            //32 bits data word 
  
  for(Int_t i=0;i<AliRICHDigit::kNddls;i++){//open all 14 DDL in parallel
    file[i].open(AliDAQ::DdlFileName("RICH",i));
    file[i].write((char*)&header,sizeof(header));                //write dummy header as place holder, actual will be written later when total size of DDL is known
    cnt[i]=0; //reset counters
  }
  
  for(Int_t iCh=0;iCh<fNcham;iCh++){ //digits are stored on chamber by chamber basis   
    TClonesArray *pDigs=(TClonesArray*)fDig->UncheckedAt(iCh);
    for(Int_t iDig=0;iDig<pDigs->GetEntriesFast();iDig++){//digits loop for a given chamber
      AliRICHDigit *pDig=(AliRICHDigit*)pDigs->At(iDig);
      Int_t ddl=pDig->Dig2Raw(w32);  //ddl is 0..13 
      file[ddl].write((char*)&w32,sizeof(w32));  cnt[ddl]++;//write formated digit to the propriate file (as decided in Dig2Raw) and increment corresponding counter
    }//digits loop for a given chamber
  }//chambers loop    
  for(Int_t i=0;i<AliRICHDigit::kNddls;i++){
    header.fSize=sizeof(header)+cnt[i]*sizeof(w32); //now calculate total number of bytes for each DDL file  
    header.SetAttribute(0); 
    file[i].seekp(0); file[i].write((char*)&header,sizeof(header));//rewrite DDL header with fSize field properly set
    file[i].close();                                               //close DDL file
  }
  GetLoader()->UnloadDigits();
  AliDebug(1,"Stop.");      
}//Digits2Raw()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Float_t AliRICHv1::Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)
{
// Correction for Fresnel   ???????????
// Arguments:   ene - photon energy [GeV],
//              PDOTI=COS(INC.ANG.), PDOTR=COS(POL.PLANE ROT.ANG.)
//   Returns:  
    Float_t en[36] = {5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,
		      6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,
		      7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5};
    Float_t csin[36] = {2.14,2.21,2.33,2.48,2.76,2.97,2.99,2.59,2.81,3.05,
			2.86,2.53,2.55,2.66,2.79,2.96,3.18,3.05,2.84,2.81,2.38,2.11,
			2.01,2.13,2.39,2.73,3.08,3.15,2.95,2.73,2.56,2.41,2.12,1.95,
			1.72,1.53};
    Float_t csik[36] = {0.,0.,0.,0.,0.,0.196,0.408,0.208,0.118,0.49,0.784,0.543,
	 		0.424,0.404,0.371,0.514,0.922,1.102,1.139,1.376,1.461,1.253,0.878,
			0.69,0.612,0.649,0.824,1.347,1.571,1.678,1.763,1.857,1.824,1.824,
			1.714,1.498};
    Float_t xe=ene;
    Int_t  j=Int_t(xe*10)-49;
    Float_t cn=csin[j]+((csin[j+1]-csin[j])/0.1)*(xe-en[j]);
    Float_t ck=csik[j]+((csik[j+1]-csik[j])/0.1)*(xe-en[j]);

    //FORMULAE FROM HANDBOOK OF OPTICS, 33.23 OR
    //W.R. HUNTER, J.O.S.A. 54 (1964),15 , J.O.S.A. 55(1965),1197

    Float_t sinin=TMath::Sqrt(1-pdoti*pdoti);
    Float_t tanin=sinin/pdoti;

    Float_t c1=cn*cn-ck*ck-sinin*sinin;
    Float_t c2=4*cn*cn*ck*ck;
    Float_t aO=TMath::Sqrt(0.5*(TMath::Sqrt(c1*c1+c2)+c1));
    Float_t b2=0.5*(TMath::Sqrt(c1*c1+c2)-c1);
    
    Float_t rs=((aO-pdoti)*(aO-pdoti)+b2)/((aO+pdoti)*(aO+pdoti)+b2);
    Float_t rp=rs*((aO-sinin*tanin)*(aO-sinin*tanin)+b2)/((aO+sinin*tanin)*(aO+sinin*tanin)+b2);
    

    //CORRECTION FACTOR FOR SURFACE ROUGHNESS
    //B.J. STAGG  APPLIED OPTICS, 30(1991),4113

    Float_t sigraf=18.;
    Float_t lamb=1240/ene;
    Float_t fresn;
 
    Float_t  rO=TMath::Exp(-(4*TMath::Pi()*pdoti*sigraf/lamb)*(4*TMath::Pi()*pdoti*sigraf/lamb));

    if(pola)
    {
	Float_t pdotr=0.8;                                 //DEGREE OF POLARIZATION : 1->P , -1->S
	fresn=0.5*(rp*(1+pdotr)+rs*(1-pdotr));
    }
    else
	fresn=0.5*(rp+rs);
      
    fresn = fresn*rO;
    return fresn;
}//Fresnel()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::Print(Option_t *option)const
{
// Debug printout
  TObject::Print(option);
  Printf("Total number of MIP reached radiator        %9i",fCounters(kMipEnterRad));
  Printf("Total number of Ckov created                %9i",fCounters(kCkovNew));
  Printf("number of Ckov created in radiator          %9i",fCounters(kCkovNewRad));
  Printf("number of Ckov created in window            %9i",fCounters(kCkovNewWin));
  Printf("number of Ckov created in proximity gap     %9i",fCounters(kCkovNewProxGap));
  Printf("number of Ckov created in amplification gap %9i",fCounters(kCkovNewAmpGap));
  Printf("number of Ckov reached PC                   %9i",fCounters(kCkovEnterPc));
  Printf("number of photelectrons                     %9i",fCounters(kPhotoEle));
}//void AliRICH::Print(Option_t *option)const
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::StepCount()
{
// Count number of ckovs created  
  Int_t copy;
  if(gMC->TrackCharge()        &&gMC->CurrentVolID(copy)==fIdRad    &&gMC->IsTrackEntering()                ) fCounters(kMipEnterRad)++;  
  if(gMC->TrackPid()==kCerenkov                                     &&gMC->IsNewTrack()                     ) fCounters(kCkovNew)++;      
  if(gMC->TrackPid()==kCerenkov&&gMC->CurrentVolID(copy)==fIdRad    &&gMC->IsNewTrack()                     ) fCounters(kCkovNewRad)++;   
  if(gMC->TrackPid()==kCerenkov&&gMC->CurrentVolID(copy)==fIdWin    &&gMC->IsNewTrack()                     ) fCounters(kCkovNewWin)++;   
  if(gMC->TrackPid()==kCerenkov&&gMC->CurrentVolID(copy)==fIdProxGap&&gMC->IsNewTrack()                     ) fCounters(kCkovNewProxGap)++;
  if(gMC->TrackPid()==kCerenkov&&gMC->CurrentVolID(copy)==fIdAmpGap &&gMC->IsNewTrack()                     ) fCounters(kCkovNewAmpGap)++;
  if(gMC->TrackPid()==kCerenkov&&gMC->CurrentVolID(copy)==fIdPc     &&gMC->IsTrackEntering()                ) fCounters(kCkovEnterPc)++;  
  if(gMC->TrackPid()==kCerenkov&&gMC->CurrentVolID(copy)==fIdPc     &&gMC->IsTrackEntering() &&gMC->Edep()>0) fCounters(kPhotoEle)++;       
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::StepHistory()
{
// This methode is invoked from StepManager() in order to print out 
  static Int_t iStepN;
  const char *sParticle;
  switch(gMC->TrackPid()){
    case kProton:      sParticle="PROTON"    ;break;
    case kNeutron:     sParticle="neutron"   ;break;
    case kGamma:       sParticle="gamma"     ;break;
    case kCerenkov:    sParticle="CKOV"    ;break;
    case kPi0:         sParticle="Pi0"       ;break;  
    case kPiPlus:      sParticle="Pi+"       ;break;  
    case kPiMinus:     sParticle="Pi-"       ;break;  
    case kElectron:    sParticle="electron"  ;break;  
    default:           sParticle="not known" ;break;
  }

  TString flag="fanny combination";
  if(gMC->IsTrackAlive())
      if(gMC->IsTrackEntering())      flag="enters to";
      else if(gMC->IsTrackExiting())  flag="exits from";
      else if(gMC->IsTrackInside())   flag="inside";
  else
      if(gMC->IsTrackStop())          flag="stoped in";        
  
  Int_t vid=0,copy=0;
  TString path=gMC->CurrentVolName(); path.Prepend("-");path.Prepend(gMC->CurrentVolOffName(1));//current volume and his mother are always there
  vid=gMC->CurrentVolOffID(2,copy);  if(vid) {path.Prepend("-");path.Prepend(gMC->VolName(vid));}
  vid=gMC->CurrentVolOffID(3,copy);  if(vid) {path.Prepend("-");path.Prepend(gMC->VolName(vid));}
  
  Printf("Step %i: %s (%i) %s %s m=%.6f GeV q=%.1f dEdX=%.4f",iStepN,sParticle,gMC->TrackPid(),flag.Data(),path.Data(),gMC->TrackMass(),gMC->TrackCharge(),gMC->Edep()*1e9);
  
  Printf("Step %i: tid=%i flags alive=%i disap=%i enter=%i exit=%i inside=%i out=%i stop=%i new=%i",
                            iStepN, gAlice->GetMCApp()->GetCurrentTrackNumber(),
                            gMC->IsTrackAlive(), gMC->IsTrackDisappeared(),gMC->IsTrackEntering(), gMC->IsTrackExiting(),
                            gMC->IsTrackInside(),gMC->IsTrackOut(),        gMC->IsTrackStop(),     gMC->IsNewTrack());
  
  Float_t a,z,den,rad,abs; a=z=den=rad=abs=-1;
  Int_t mid=gMC->CurrentMaterial(a,z,den,rad,abs);
  Printf("Step %i: id=%i a=%7.2f z=%7.2f den=%9.4f rad=%9.2f abs=%9.2f\n\n",iStepN,mid,a,z,den,rad,abs);
  iStepN++;
}//StepHistory()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHv1::StepManager()
{
// Full Step Manager.
// Arguments: none
//   Returns: none           
//  StepHistory(); return; //uncomment to print tracks history
//  StepCount(); return;     //uncomment to count photons
  
  Int_t          copy; //volume copy aka node
  static Int_t   iCham;
//Treat photons    
  static TLorentzVector cerX4;
  if((gMC->TrackPid()==kCerenkov||gMC->TrackPid()==kFeedback)&&gMC->CurrentVolID(copy)==fIdPc){//photon in PC
    if(gMC->Edep()>0){//photon survided QE test i.e. produces electron
      if(IsLostByFresnel()){ gMC->StopTrack(); return;} //photon lost due to fersnel reflection        
      gMC->TrackPosition(cerX4); gMC->CurrentVolOffID(2,iCham);//RICH-Rppf-Rpc
      HitAdd(iCham,gMC->GetStack()->GetCurrentTrackNumber(),gMC->TrackPid(),cerX4.Vect(),cerX4.Vect());//HIT for PHOTON in conditions CF+CSI+DE
      GenFee(iCham);
    }//photon in PC and DE >0 
  }//photon in PC
  
//Treat charged particles  
  static Float_t eloss;
  static TLorentzVector mipInX4,mipOutX4;
  if(gMC->TrackCharge() && gMC->CurrentVolID(copy)==fIdAmpGap){//MIP in amplification gap
    gMC->CurrentVolOffID(1,iCham);//RICH-Rgap
    if(gMC->IsTrackEntering()||gMC->IsNewTrack()) {//MIP in GAP entering or newly created
      eloss=0;                                                           
      gMC->TrackPosition(mipInX4);
    }else if(gMC->IsTrackExiting()||gMC->IsTrackStop()||gMC->IsTrackDisappeared()){//MIP in GAP exiting or disappeared
      eloss+=gMC->Edep();//take into account last step dEdX
      gMC->TrackPosition(mipOutX4);  
      HitAdd(iCham,gMC->GetStack()->GetCurrentTrackNumber(),gMC->TrackPid(),mipInX4.Vect(),mipOutX4.Vect(),eloss);//HIT for MIP: MIP in GAP Exiting
      GenFee(iCham,eloss);//MIP+GAP+Exit
    }else//MIP in GAP going inside
      eloss   += gMC->Edep();
  }//MIP in GAP
}//StepManager()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
