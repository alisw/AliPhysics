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

#include "AliRICH.h"
#include "AliRICHParam.h"
#include "AliRICHChamber.h"
#include "AliRICHClusterFinder.h"
#include <TArrayF.h>
#include <TGeometry.h>
#include <TBRIK.h>
#include <TTUBE.h>
#include <TFile.h>
#include <TNode.h> 
#include <TObjArray.h>
#include <TParticle.h>
#include <AliStack.h>
#include <AliMagF.h>
#include <AliRun.h>
#include <AliRunDigitizer.h>
#include <AliMC.h>
#include <AliESD.h>
#include <TVirtualMC.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <AliLog.h>
 
ClassImp(AliRICHhit)
//__________________________________________________________________________________________________
void AliRICHhit::Print(Option_t*)const
{
  AliInfo(Form("Ch=%1i, TID=%6i, eloss=%9.3f eV, in-out dist=%9.4f, OUT(%7.2f,%7.2f,%7.2f)"
      ,fChamber,fTrack,fEloss*1e9,Length(),fOutX3.X(),fOutX3.Y(),fOutX3.Z()));
}
//__________________________________________________________________________________________________
ClassImp(AliRICHdigit)
//__________________________________________________________________________________________________
void AliRICHdigit::Print(Option_t*)const
{
  AliInfo(Form("cfm=%9i, cs=%2i, x=%3i, y=%3i, q=%8.3f, TID1=%5i, TID2=%5i, TID3=%5i",
                  fCFM,fChamber,fPadX,fPadY,fQdc,fTracks[0],fTracks[1],fTracks[2]));
}
//__________________________________________________________________________________________________
ClassImp(AliRICHcluster)
//__________________________________________________________________________________________________
void AliRICHcluster::Print(Option_t*)const
{
  char *status=0;
  switch(fStatus){
    case      kRaw: status="raw"     ;break;
    case kResolved: status="resolved";break;
    case    kEmpty: status="empty"   ;break;
  }
  if(fDigits)    
    ::Info("cluster","cfm=%10i, cs=%2i, SiMa=%6i, Shape=%5i, x=%7.3f, y=%7.3f, Q=%6i, %s with %i digits",
                             fCFM,fChamber,fSize,fShape,fX,fY,fQdc,status,fDigits->GetEntriesFast());
  else
    AliInfo(Form("cfm=%10i, cs=%2i, SiMa=%6i, Shape=%5i, x=%7.3f, y=%7.3f, Q=%6i, %s with %i digits",
                             fCFM,fChamber,fSize,fShape,fX,fY,fQdc,status,0));
    
}
//__________________________________________________________________________________________________
ClassImp(AliRICH)    
//__________________________________________________________________________________________________
// RICH manager class   
//BEGIN_HTML
/*
  <img src="gif/alirich.gif">
*/
//END_HTML
//__________________________________________________________________________________________________
AliRICH::AliRICH():AliDetector(),fpParam(0),  fSdigits(0),fNsdigits(0),fDigitsNew(0),fClusters(0) 
{
//Default ctor should not contain any new operators
//AliDetector ctor deals with Hits and Digits  
  for(int i=0;i<kNchambers;i++) fNdigitsNew[i]  =0;
  for(int i=0;i<kNchambers;i++) fNclusters[i]=0;
//  fCounters.ResizeTo(20); fCounters.Zero();
}//AliRICH::AliRICH()
//__________________________________________________________________________________________________
AliRICH::AliRICH(const char *name, const char *title)
        :AliDetector(name,title),fpParam(new AliRICHParam),fSdigits(0),fNsdigits(0),fDigitsNew(0),fClusters(0)
{
//Named ctor
  AliDebug(1,"Start.");
//AliDetector ctor deals with Hits and Digits (reset them to 0, does not create them)
  CreateHits();          gAlice->GetMCApp()->AddHitList(fHits);
  fCounters.ResizeTo(20); fCounters.Zero();
  AliDebug(1,"Stop.");
}//AliRICH::AliRICH(const char *name, const char *title)
//__________________________________________________________________________________________________
AliRICH::~AliRICH()
{
//dtor
  AliDebug(1,"Start.");

  if(fpParam)    delete fpParam;
  
  if(fHits)      delete fHits;
  if(fSdigits)   delete fSdigits;
  if(fDigits)    delete fDigits;
  if(fDigitsNew) {fDigitsNew->Delete();   delete fDigitsNew;}
  if(fClusters)  {fClusters->Delete();    delete fClusters;}
  AliDebug(1,"Stop.");    
}//AliRICH::~AliRICH()
//__________________________________________________________________________________________________
void AliRICH::Hits2SDigits()
{
// Create a list of sdigits corresponding to list of hits. Every hit generates one or more sdigits.
  AliDebug(1,"Start.");
  for(Int_t iEventN=0;iEventN<GetLoader()->GetRunLoader()->GetAliRun()->GetEventsPerRun();iEventN++){//events loop
    GetLoader()->GetRunLoader()->GetEvent(iEventN);//get next event
  
    if(!GetLoader()->TreeH()) GetLoader()->LoadHits();    GetLoader()->GetRunLoader()->LoadHeader(); 
                                                          GetLoader()->GetRunLoader()->LoadKinematics();//from
    if(!GetLoader()->TreeS()) GetLoader()->MakeTree("S"); MakeBranch("S");//to
          
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      GetLoader()->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){//hits loop 
        AliRICHhit *pHit=(AliRICHhit*)Hits()->At(iHitN);//get current hit                
        TVector2 x2 = C(pHit->C())->Mrs2Pc(pHit->OutX3());//hit position in photocathode plane
        Int_t iTotQdc=P()->TotQdc(x2,pHit->Eloss());//total charge produced by hit, 0 if hit in dead zone
        if(iTotQdc==0) continue;
        //
        //need to quantize the anod....
        TVector padHit=AliRICHParam::Loc2Pad(x2);
        TVector2 padHitXY=AliRICHParam::Pad2Loc(padHit);
        TVector2 anod;
        if((x2.Y()-padHitXY.Y())>0) anod.Set(x2.X(),padHitXY.Y()+AliRICHParam::PitchAnod()/2);
        else anod.Set(x2.X(),padHitXY.Y()-AliRICHParam::PitchAnod()/2);
        //end to quantize anod
        //
        TVector area=P()->Loc2Area(anod);//determine affected pads, dead zones analysed inside
        AliDebug(1,Form("hitanod(%6.2f,%6.2f)->area(%3.0f,%3.0f)-(%3.0f,%3.0f) QDC=%4i",anod.X(),anod.Y(),area[0],area[1],area[2],area[3],iTotQdc));
        TVector pad(2);
        for(pad[1]=area[1];pad[1]<=area[3];pad[1]++)//affected pads loop
          for(pad[0]=area[0];pad[0]<=area[2];pad[0]++){                    
            Double_t padQdc=iTotQdc*P()->FracQdc(anod,pad);
            AliDebug(1,Form("current pad(%3.0f,%3.0f) with QDC  =%6.2f",pad[0],pad[1],padQdc));
            if(padQdc>0.1) AddSDigit(pHit->C(),pad,padQdc,GetLoader()->GetRunLoader()->Stack()->Particle(pHit->GetTrack())->GetPdgCode(),pHit->GetTrack());
          }//affected pads loop 
      }//hits loop
    }//prims loop
    GetLoader()->TreeS()->Fill();
    GetLoader()->WriteSDigits("OVERWRITE");
    ResetSDigits();
  }//events loop  
  GetLoader()->UnloadHits(); GetLoader()->GetRunLoader()->UnloadHeader(); GetLoader()->GetRunLoader()->UnloadKinematics();
  GetLoader()->UnloadSDigits();  
  AliDebug(1,"Stop.");
}//Hits2SDigits()
//__________________________________________________________________________________________________
void AliRICH::BuildGeometry() 
{
//Builds a TNode geometry for event display
  AliInfo("Start.");
  
  TNode *node, *subnode, *top;
  top=gAlice->GetGeometry()->GetNode("alice");

  Float_t widx=P()->SectorSizeX();
  Float_t leny=P()->SectorSizeY();
  Float_t dz = P()->Zfreon()+P()->Zwin()+P()->Pc2Win();
  Float_t dead = P()->DeadZone();

  new TBRIK("RICH","RICH","void",widx+dead/2,leny+leny/2+dead,dz+0.1); //RICH chamber
  new TBRIK("RPC" ,"RPC" ,"void",widx/2,leny/2,0.01);                  //RICH sector 

  for(int i=1;i<=kNchambers;i++){
    top->cd();
    node = new TNode(Form("RICH%i",i),Form("RICH%i",i),"RICH",C(i)->Center().X(),C(i)->Center().Y(),C(i)->Center().Z(),C(i)->RotMatrixName());
    node->SetLineColor(kRed);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","RPC",-widx/2-dead/2,-leny-dead/2,dz,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","RPC", widx/2+dead/2,-leny-dead/2,dz,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","RPC",-widx/2-dead/2,           0,dz,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","RPC", widx/2+dead/2,           0,dz,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","RPC",-widx/2-dead/2, leny+dead/2,dz,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","RPC", widx/2+dead/2, leny+dead/2,dz,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);
  }

  AliDebug(1,"Stop.");    
}//void AliRICH::BuildGeometry()

//______________________________________________________________________________
void AliRICH::CreateMaterials()
{
// Definition of available RICH materials  
        
  Int_t   material=0; //tmp material id number
  Float_t a=0,z=0,den=0,radl=0,absl=0; //tmp material parameters
  Float_t tmaxfd=-10.0, deemax=-0.2, stemax=-0.1,epsil=0.001, stmin=-0.001; 
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
    
  Float_t aAir[4]={12.,14.,16.,36.};  Float_t zAir[4]={6.,7.,8.,18.}; Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};//total 0.9999999
  AliMixture(++material, "RichAir",aAir,zAir,den=0.00120479,4,wAir);                                          //1 (Air) 0.01% C 75% N  23% O 1% Ar
  AliMedium(kAir, "RichAir",material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++material, "RichRohacell", a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);                   //2 Rohacell 51 C-equiv radl rad cover
  AliMedium(kRoha, "RichRohacell", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aQuartz[2]={28.09,16.0};  Float_t  zQuartz[2]={14.00, 8.0};  Float_t  wQuartz[2]={1,2};
  AliMixture(++material, "RichSiO2",aQuartz,zQuartz,den=2.64,-2, wQuartz);                                    //3 Quarz (SiO2) -trasparent rad window
  AliMedium(kSiO2, "RichSiO2",material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aFreon[2]={12,19};  Float_t  zFreon[2]={6,9};  Float_t wmatFreon[2]={6,14};
  AliMixture(++material, "RichC6F14",aFreon,zFreon,den=1.68,-2,wmatFreon);                                    //4 Freon (C6F14) 
  AliMedium(kC6F14, "RichC6F14",material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t aMethane[2]={12.01,1}; Float_t zMethane[2]={6,1}; Float_t wMethane[2]={1,4};
  AliMixture (++material, "RichCH4", aMethane, zMethane, den=7.17e-4,-2, wMethane);                        //5,9 methane (CH4) normal and for Gap    
  AliMedium(kCH4, "RichCH4"   , material, 1, isxfld, sxmgmx, tmaxfd, stemax,  deemax, epsil,  stmin);  
  AliMixture (++material, "RichCH4", aMethane, zMethane, den=7.17e-4,-2, wMethane);                        //5,9 methane (CH4) normal and for Gap    
  AliMedium(kGap, "RichCH4gap$", material, 1, isxfld, sxmgmx, tmaxfd, 0.1   , -deemax, epsil, -stmin);
    
  AliMaterial(++material, "RichCsI",      a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);                   //6 CsI-radl equivalent
  AliMedium(kCsI, "RichCsI$", material, 1, isxfld, sxmgmx,tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++material, "GridCu",    a=63.54,z=29.0,den=8.96,    radl=1.43,   absl=0);                   //7 anode grid (Cu) 
  AliMedium(kGridCu, "GRID$", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture (++material, "OpSiO2",aQuartz, zQuartz, den=2.64, -2, wQuartz);                             //8 Quarz (SiO2) - opaque
  AliMedium(kOpSiO2, "QUARTZO$",material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++material, "ALUM",     a=26.98,z=13.0,den=2.699,     radl=8.9,    absl=0);                 //10 aluminium sheet (Al)
  AliMedium(kAl, "ALUMINUM$", material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t aGlass[5]={12.01,28.09,16,10.8,23}; Float_t zGlass[5]={6,14,8,5,11};  Float_t wGlass[5]={0.5,0.105,0.355,0.03,0.01};
  AliMixture(++material,"GLASS",aGlass, zGlass, den=1.74, 5, wGlass);                                    //11 Glass 50%-C 10.5%-Si 35.5%-O 3%-B 1%-Na
  AliMedium(kGlass, "GLASS", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++material, "COPPER$",  a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);                   //12 Cu
  AliMedium(kCu, "PCB_COPPER", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  
  AliMaterial(++material, "W$",  a=183.84,z=74.0,den=19.3,    radl=0.35,    absl=185.0/den);              //13 W - anod wires
  AliMedium(kW, "W", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  
  if(P()->IsRadioSrc()){
    AliInfo("Special radioactive source materials");
    AliMaterial(++material, "STEEL$",  a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);                  //14 Steel
    AliMedium(kSteel, "STEEL", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
    AliMaterial(++material, "PERPEX$",  a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);                 //15 Perpex
    AliMedium(kPerpex, "PERPEX", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    
    AliMaterial(++material, "STRONZIUM$",  a=87.62,z=38.0,den=2.54,    radl=4.24,    absl=0);             //16 Sr90
    AliMedium(kSr90, "STRONZIUM", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  }
  
//Optical properties:
#include "Opticals.h"
  gMC->SetCerenkov((*fIdtmed)[kAir]      , kNbins, aPckov, aAbsCH4,    aQeAll, aIdxCH4);       //1 Air
  gMC->SetCerenkov((*fIdtmed)[kRoha]     , kNbins, aPckov, aAbsCH4,    aQeAll, aIdxCH4);       //2 Honeycomb  
  gMC->SetCerenkov((*fIdtmed)[kSiO2]     , kNbins, aPckov, aAbsSiO2,   aQeAll, aIdxSiO2);      //3 Quartz SiO2 
  gMC->SetCerenkov((*fIdtmed)[kC6F14]    , kNbins, aPckov, aAbsC6F14,  aQeAll, aIdxC6F14);     //4 Freon C6F14
  gMC->SetCerenkov((*fIdtmed)[kCH4]      , kNbins, aPckov, aAbsCH4,    aQeAll, aIdxCH4);       //5 Methane CH4 
  gMC->SetCerenkov((*fIdtmed)[kCsI]      , kNbins, aPckov, aAbsCsI,    aQeCsI, aIdxCH4);       //6 CsI
  gMC->SetCerenkov((*fIdtmed)[kGridCu]   , kNbins, aPckov, aAbsGrid,   aQeAll, aIdxGrid);      //7 grid Cu
  gMC->SetCerenkov((*fIdtmed)[kOpSiO2]   , kNbins, aPckov, aAbsOpSiO2, aQeAll, aIdxOpSiO2);    //8 Opaque quartz SiO2
  gMC->SetCerenkov((*fIdtmed)[kGap]      , kNbins, aPckov, aAbsCH4,    aQeAll, aIdxCH4);       //9 Special methane gap
  gMC->SetCerenkov((*fIdtmed)[kAl]       , kNbins, aPckov, aAbsGrid,   aQeAll, aIdxGrid);      //10 Aluminium
  gMC->SetCerenkov((*fIdtmed)[kGlass]    , kNbins, aPckov, aAbsOpSiO2, aQeAll, aIdxOpSiO2);    //11 Glass    
}//void AliRICH::CreateMaterials()
//__________________________________________________________________________________________________
Float_t AliRICH::Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)const
{

    //ENE(EV), PDOTI=COS(INC.ANG.), PDOTR=COS(POL.PLANE ROT.ANG.)
    
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
    return(fresn);
}//Fresnel()
//__________________________________________________________________________________________________
Float_t AliRICH::AbsoCH4(Float_t x)const
{
//Evaluate the absorbtion lenght of CH4
  Float_t sch4[9] = {.12,.16,.23,.38,.86,2.8,7.9,28.,80.};              //MB X 10^22
  Float_t em[9] = {8.1,8.158,8.212,8.267,8.322,8.378,8.435,8.493,8.55};
  const Float_t kLoschmidt=2.686763e19;                                      // LOSCHMIDT NUMBER IN CM-3
  const Float_t kPressure=750.,kTemperature=283.;                                      
  const Float_t kPn=kPressure/760.;
  const Float_t kTn=kTemperature/273.16;
  const Float_t kC0=-1.655279e-1;
  const Float_t kC1=6.307392e-2;
  const Float_t kC2=-8.011441e-3;
  const Float_t kC3=3.392126e-4;
    		
  Float_t crossSection=0;                        
  if (x<7.75) 
    crossSection=.06e-22;
  else if(x>=7.75 && x<=8.1){                 //------ METHANE CROSS SECTION cm-2 ASTROPH. J. 214, L47 (1978)                                               
	crossSection=(kC0+kC1*x+kC2*x*x+kC3*x*x*x)*1.e-18;
  }else if (x> 8.1){
    Int_t j=0;
    while (x<=em[j] || x>=em[j+1]){
      j++;
      Float_t a=(sch4[j+1]-sch4[j])/(em[j+1]-em[j]);
      crossSection=(sch4[j]+a*(x-em[j]))*1e-22;
    }
  }//if
    
    Float_t density=kLoschmidt*kPn/kTn; //CH4 molecular density 1/cm-3
    return 1./(density*crossSection);
}//AbsoCH4()
//__________________________________________________________________________________________________
void AliRICH::MakeBranch(Option_t* option)
{
//Create Tree branches for the RICH.
  AliDebug(1,Form("Start with option= %s.",option));
    
  const Int_t kBufferSize = 4000;
      
  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");
  const char *cS = strstr(option,"S");

  if(cH&&TreeH()){//H
    CreateHits();      //branch will be created in AliDetector::MakeBranch
  }//H     
  AliDetector::MakeBranch(option);//this is after cH because we need to guarantee that fHits array is created
      
  if(cS&&fLoader->TreeS()){//S  
    CreateSDigits();   MakeBranchInTree(fLoader->TreeS(),"RICH",&fSdigits,kBufferSize,0) ;
  }//S
   
  if(cD&&fLoader->TreeD()){//D
    CreateDigits();
    for(Int_t i=0;i<kNchambers;i++){ 
      MakeBranchInTree(fLoader->TreeD(),Form("%s%d",GetName(),i+1),&((*fDigitsNew)[i]),kBufferSize,0);
    }
  }//D
  
  if(cR&&fLoader->TreeR()){//R
    CreateClusters();
    for(Int_t i=0;i<kNchambers;i++)
      MakeBranchInTree(fLoader->TreeR(),Form("%sClusters%d",GetName(),i+1), &((*fClusters)[i]), kBufferSize, 0);    
  }//R
  AliDebug(1,"Stop.");   
}//void AliRICH::MakeBranch(Option_t* option)
//__________________________________________________________________________________________________
void AliRICH::SetTreeAddress()
{
//Set branch address for the Hits and Digits Tree.
  AliDebug(1,"Start.");
      
  TBranch *branch;
    
  if(fLoader->TreeH()){//H
    AliDebug(1,"tree H is requested.");
    CreateHits();//branch map will be in AliDetector::SetTreeAddress    
  }//H
  AliDetector::SetTreeAddress();//this is after TreeH because we need to guarantee that fHits array is created

  if(fLoader->TreeS()){//S
    AliDebug(1,"tree S is requested.");
    branch=fLoader->TreeS()->GetBranch(GetName());        if(branch){CreateSDigits();   branch->SetAddress(&fSdigits);}
  }//S
    
  if(fLoader->TreeD()){//D    
    AliDebug(1,"tree D is requested.");
    for(int i=0;i<kNchambers;i++){      
      branch=fLoader->TreeD()->GetBranch(Form("%s%d",GetName(),i+1)); 
      if(branch){CreateDigits(); branch->SetAddress(&((*fDigitsNew)[i]));}
    }
  }//D
    
  if(fLoader->TreeR()){//R
    AliDebug(1,"tree R is requested.");
    for(int i=0;i<kNchambers;i++){         
      branch=fLoader->TreeR()->GetBranch(Form("%sClusters%d" ,GetName(),i+1));
      if(branch){CreateClusters(); branch->SetAddress(&((*fClusters)[i]));}
    }
  }//R
  AliDebug(1,"Stop.");
}//void AliRICH::SetTreeAddress()
//__________________________________________________________________________________________________
void AliRICH::Print(Option_t *option)const
{
//Debug printout
  TObject::Print(option);
  P()->Print();
  fCounters.Print();
}//void AliRICH::Print(Option_t *option)const
//__________________________________________________________________________________________________
void AliRICH::CreateGeometry()
{
//Creates detailed geometry simulation (currently GEANT volumes tree)         
  AliDebug(1,"Start main.");
  Double_t cm=1,mm=0.1*cm,mkm=0.001*cm;
  Float_t par[3];
  Int_t id=0;
       
//place chambers into mother volume ALIC
  par[0]=(6*mm+1681*mm+6*mm)/2;par[1]=(6*mm+1466*mm+6*mm)/2;par[2]=(80*mm+40*mm)*2/2;gMC->Gsvolu("RICH","BOX ",(*fIdtmed)[kCH4],par,3);//2033P1  z p84 TDR
  if(P()->IsRadioSrc()){
    AliInfo("Special test beam geometry");
    gMC->Gspos("RICH",1,"ALIC",0,0,0,0, "ONLY");  //single RICH chamber without rotation test beam geometry
  }else
    for(int i=1;i<=kNchambers;i++){ //full RICH with 7 chambers
      AliMatrix(id,C(i)->ThetaXd(),C(i)->PhiXd(),  C(i)->ThetaYd(),C(i)->PhiYd(),  C(i)->ThetaZd(),C(i)->PhiZd());
      gMC->Gspos("RICH",i,"ALIC",C(i)->Center().X(),C(i)->Center().Y(),C(i)->Center().Z(),id, "ONLY");
    }
//Pad Panel frame  6 sectors
  par[0]=648*mm/2;par[1]=  411*mm/2;par[2]=40  *mm/2;gMC->Gsvolu("RPPF","BOX ",(*fIdtmed)[kAl]  ,par,3);//PPF 2001P2 inner size of the slab by 1mm more
  par[0]=181*mm/2;par[1]=89.25*mm/2;par[2]=38.3*mm/2;gMC->Gsvolu("PPFL","BOX ",(*fIdtmed)[kAir] ,par,3);//large whole
  par[0]=114*mm/2;par[1]=89.25*mm/2;par[2]=38.3*mm/2;gMC->Gsvolu("PPFS","BOX ",(*fIdtmed)[kAir] ,par,3);//small whole
  par[0]=644*mm/2;par[1]=  407*mm/2;par[2]= 1.7*mm/2;gMC->Gsvolu("RPC ","BOX ",(*fIdtmed)[kCsI] ,par,3);//by 0.2 mm more then actual size (PCB 2006P1)
  
  gMC->Gspos("RPPF",1,"RICH",    -335*mm,      -433*mm,  8*cm+20*mm,  0,"ONLY");//F1 2040P1 z p.84 TDR
  gMC->Gspos("RPPF",2,"RICH",    +335*mm,      -433*mm,  8*cm+20*mm,  0,"ONLY");
  gMC->Gspos("RPPF",3,"RICH",    -335*mm,         0*mm,  8*cm+20*mm,  0,"ONLY");
  gMC->Gspos("RPPF",4,"RICH",    +335*mm,         0*mm,  8*cm+20*mm,  0,"ONLY");
  gMC->Gspos("RPPF",5,"RICH",    -335*mm,      +433*mm,  8*cm+20*mm,  0,"ONLY");
  gMC->Gspos("RPPF",6,"RICH",    +335*mm,      +433*mm,  8*cm+20*mm,  0,"ONLY");  
  gMC->Gspos("RPC ",1,"RPPF",       0*mm,         0*mm,   -19.15*mm,  0,"ONLY");//PPF 2001P2 
  gMC->Gspos("PPFL",1,"RPPF",  -224.5*mm,  -151.875*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFL",2,"RPPF",  -224.5*mm,  - 50.625*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFL",3,"RPPF",  -224.5*mm,  + 50.625*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFL",4,"RPPF",  -224.5*mm,  +151.875*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFS",1,"RPPF",  - 65.0*mm,  -151.875*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFS",2,"RPPF",  - 65.0*mm,  - 50.625*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFS",3,"RPPF",  - 65.0*mm,  + 50.625*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFS",4,"RPPF",  - 65.0*mm,  +151.875*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFS",5,"RPPF",  + 65.0*mm,  -151.875*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFS",6,"RPPF",  + 65.0*mm,  - 50.625*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFS",7,"RPPF",  + 65.0*mm,  + 50.625*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFS",8,"RPPF",  + 65.0*mm,  +151.875*mm,     0.85*mm,  0,"ONLY"); 
  gMC->Gspos("PPFL",5,"RPPF",  +224.5*mm,  -151.875*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFL",6,"RPPF",  +224.5*mm,  - 50.625*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFL",7,"RPPF",  +224.5*mm,  + 50.625*mm,     0.85*mm,  0,"ONLY");
  gMC->Gspos("PPFL",8,"RPPF",  +224.5*mm,  +151.875*mm,     0.85*mm,  0,"ONLY");
//Gap - anod wires 6 copies to RICH
  par[0]=648*mm/2;par[1]=  411*mm/2 ;par[2]=4.45*mm/2;gMC->Gsvolu("RGAP ","BOX ",(*fIdtmed)[kCH4] ,par,3);//xy as PPF 2001P2 z WP 2099P1
  par[0]=  0*mm  ;par[1]=  20*mkm/2 ;par[2]= 648*mm/2;gMC->Gsvolu("RANO","TUBE",(*fIdtmed)[kW]   ,par,3);//WP 2099P1 z = gap x PPF 2001P2
  AliMatrix(id,180,0, 90,90, 90,0); //wires along x
  
  gMC->Gspos("RGAP",1,"RICH",    -335*mm,      -433*mm,8*cm-2.225*mm, 0,"ONLY"); //F1 2040P1 z WP 2099P1
  gMC->Gspos("RGAP",2,"RICH",    +335*mm,      -433*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("RGAP",3,"RICH",    -335*mm,         0*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("RGAP",4,"RICH",    +335*mm,         0*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("RGAP",5,"RICH",    -335*mm,      +433*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("RGAP",6,"RICH",    +335*mm,      +433*mm,8*cm-2.225*mm, 0,"ONLY"); 
  for(int i=1;i<=96;i++)
    gMC->Gspos("RANO",i,"RGAP",     0*mm, -411/2*mm+i*4*mm, 0.185*mm, id,"ONLY"); //WP 2099P1  
//Radiator 3 modules
  par[0]=1330*mm/2 ;par[1]= 413*mm/2  ;par[2]=  24*mm/2;  gMC->Gsvolu("RRAD","BOX ",(*fIdtmed)[kC6F14]     ,par,3); // Rad 2011P1
  par[0]=1330*mm/2 ;par[1]= 413*mm/2  ;par[2]=   4*mm/2;  gMC->Gsvolu("RRFR","BOX ",(*fIdtmed)[kRoha]      ,par,3); //front 
  par[0]=1330*mm/2 ;par[1]= 413*mm/2  ;par[2]=   5*mm/2;  gMC->Gsvolu("RRWI","BOX ",(*fIdtmed)[kSiO2]      ,par,3); //window
  par[0]=1330*mm/2 ;par[1]=   5*mm/2  ;par[2]=  15*mm/2;  gMC->Gsvolu("RRLO","BOX ",(*fIdtmed)[kRoha]      ,par,3); //long side  
  par[0]=  10*mm/2 ;par[1]= 403*mm/2  ;par[2]=  15*mm/2;  gMC->Gsvolu("RRSH","BOX ",(*fIdtmed)[kRoha]      ,par,3); //short side 
  par[0]=   0      ;par[1]=  10*mm/2  ;par[2]=  15*mm/2;  gMC->Gsvolu("RRSP","TUBE",(*fIdtmed)[kSiO2]      ,par,3); //spacer        
    
  gMC->Gspos("RRAD",1,"RICH",   0*mm,-434*mm,   -12*mm,  0,"ONLY"); //radiator to RICH
  gMC->Gspos("RRAD",2,"RICH",   0*mm,   0*mm,   -12*mm,  0,"ONLY"); 
  gMC->Gspos("RRAD",3,"RICH",   0*mm,+434*mm,   -12*mm,  0,"ONLY"); 
    gMC->Gspos("RRFR",1,"RRAD",   0*mm,   0*mm, -10.0*mm,  0,"ONLY"); //front cover 
    gMC->Gspos("RRWI",1,"RRAD",   0*mm,   0*mm,   9.5*mm,  0,"ONLY"); //quartz window (back cover)
    gMC->Gspos("RRLO",1,"RRAD",   0*mm,-204*mm,  -0.5*mm,  0,"ONLY"); //long side
    gMC->Gspos("RRLO",2,"RRAD",   0*mm,+204*mm,  -0.5*mm,  0,"ONLY"); //long side
    gMC->Gspos("RRSH",1,"RRAD",-660*mm,   0*mm,  -0.5*mm,  0,"ONLY"); //short side
    gMC->Gspos("RRSH",2,"RRAD",+660*mm,   0*mm,  -0.5*mm,  0,"ONLY"); //short side 
    for(int i=0;i<3;i++)
      for(int j=0;j<10;j++)
        gMC->Gspos("RRSP",10*i+j,"RRAD",-1330*mm/2+116*mm+j*122*mm,(i-1)*105*mm,-0.5*mm,0,"ONLY");

//Sandbox  
  par[0]=1419*mm/2 ;par[1]=1378*mm/2;par[2]=50.5*mm/2; gMC->Gsvolu("RSNB","BOX ",(*fIdtmed)[kAir]  ,par,3);  //2072P1   
  par[0]=1419*mm/2 ;par[1]=1378*mm/2;par[2]= 0.5*mm/2; gMC->Gsvolu("RSCO","BOX ",(*fIdtmed)[kAl]   ,par,3);  //cover
  par[0]=1359*mm/2 ;par[1]=1318*mm/2;par[2]=49.5*mm/2; gMC->Gsvolu("RSHO","BOX ",(*fIdtmed)[kRoha] ,par,3); //honeycomb structure 
  
  gMC->Gspos("RSNB",1,"RICH",   0*mm, 0*mm, -73.75*mm, 0,"ONLY"); //p.84 TDR sandbox to rich
    gMC->Gspos("RSHO",1,"RSNB", 0*mm, 0*mm,      0*mm, 0,"ONLY"); //2072P1 honeycomv to sandbox
    gMC->Gspos("RSCO",1,"RSNB", 0*mm, 0*mm,    +25*mm, 0,"ONLY"); //cover to sandbox
    gMC->Gspos("RSCO",2,"RSNB", 0*mm, 0*mm,    -25*mm, 0,"ONLY"); //cover to sandbox
             
  AliDebug(1,"Stop main.");  
}//void AliRICH::CreateGeometry()
//__________________________________________________________________________________________________
void AliRICH::ControlPlots()
{ 
// Creates a set of hists to control the results of simulation. Hists are in file $HOME/RCP.root
     
  TH1F *pHxD=0,*pHyD=0,*pNumClusH1=0,
                   *pQdcH1=0,       *pSizeH1=0,
                   *pPureMipQdcH1=0,*pPureMipSizeH1=0,
                   *pPureCerQdcH1=0,*pPureCerSizeH1=0,
                   *pPureFeeQdcH1=0,*pPureFeeSizeH1=0,
                   *pMipQdcH1=0,    *pPhotQdcH1=0;  
  TH2F *pMapH2=0,*pPureMipMapH2=0,*pPureCerMapH2=0,*pPureFeeMapH2=0;
  
  Bool_t isDig =!GetLoader()->LoadDigits();
  Bool_t isClus=!GetLoader()->LoadRecPoints();

  if(!isDig && !isClus){AliError("No digits and clusters! Nothing to do.");return;}
  
  TStopwatch sw;TDatime time;
    
  TFile *pFile = new TFile("$(HOME)/RCP.root","RECREATE");   
  
  if(isDig){
    AliInfo("Digits available");
    pHxD=new TH1F("HitDigitDiffX","Hit-Digits diff X all chambers;diff [cm]",100,-10,10); 
    pHyD=new TH1F("HitDigitDiffY","Hit-Digits diff Y all chambers;diff [cm]",100,-10,10); 
  }//isDig
  
  if(isClus){ 
    cout<<"Clusters available\n";
    pNumClusH1=new TH1F("NumClusPerEvent","Number of clusters per event;number",50,0,49);
    
    pQdcH1        =new TH1F("ClusQdc",   "Cluster Charge all chambers;q [QDC]",R()->P()->MaxQdc(),0,R()->P()->MaxQdc());
    pSizeH1       =new TH1F("ClusSize",  "Cluster size all chambers;size [number of pads in cluster]",100,0,100);
    pMapH2        =new TH2F("ClusMap",   "Cluster map;x [cm];y [cm]",1000,0,R()->P()->PcSizeX(),1000,0,R()->P()->PcSizeY());
  
    pMipQdcH1     =new TH1F("QdcMip"      ,"MIP Cluster Charge all chambers;q [QDC]",R()->P()->MaxQdc(),0,R()->P()->MaxQdc());
    pPhotQdcH1    =new TH1F("QdcPhot"     ,"Cer+Fee Cluster Charge all chambers;q [QDC]",R()->P()->MaxQdc(),0,R()->P()->MaxQdc());
        
    pPureMipQdcH1 =new TH1F("QdcPureMip"  ,"MIP only Cluster Charge all chambers;q [QDC]",R()->P()->MaxQdc(),0,R()->P()->MaxQdc());
    pPureMipSizeH1=new TH1F("SizePureMip" ,"MIP only Cluster size all chambers;size [number of pads in cluster]",100,0,100);
    pPureMipMapH2 =new TH2F("MapPureMip"  ,"MIP only Cluster map;x [cm];y [cm]",1000,0,R()->P()->PcSizeX(),1000,0,R()->P()->PcSizeY());
  
    pPureCerQdcH1 =new TH1F("QdcPureCer"  ,"Cerenkov only Cluster Charge all chambers;q [QDC]",R()->P()->MaxQdc(),0,R()->P()->MaxQdc());
    pPureCerSizeH1=new TH1F("SizePureCer" ,"Cernekov only Cluster size all chambers;size [number of pads in cluster]",100,0,100);
    pPureCerMapH2 =new TH2F("MapPureCer"  ,"Cerenkov only Cluster map;x [cm];y [cm]",1000,0,R()->P()->PcSizeX(),1000,0,R()->P()->PcSizeY());
    
    pPureFeeQdcH1 =new TH1F("QdcPureFee"  ,"Feedback only Cluster Charge all chambers;q [QDC]",R()->P()->MaxQdc(),0,R()->P()->MaxQdc());
    pPureFeeSizeH1=new TH1F("SizePureFee" ,"Feedback only Cluster size all chambers;size [number of pads in cluster]",100,0,100);
    pPureFeeMapH2 =new TH2F("MapPureFee"  ,"Feedback only Cluster map;x [cm];y [cm]",1000,0,R()->P()->PcSizeX(),1000,0,R()->P()->PcSizeY());
  }//isClus
  
  for(Int_t iEvtN=0;iEvtN < GetLoader()->GetRunLoader()->GetAliRun()->GetEventsPerRun();iEvtN++){//events loop
    GetLoader()->GetRunLoader()->GetEvent(iEvtN);    //gets current event
    
    if(!GetLoader()->TreeH()) GetLoader()->LoadHits();
    for(Int_t iPrimN=0;iPrimN < GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
       GetLoader()->TreeH()->GetEntry(iPrimN);      
    }
    
    if(isClus) GetLoader()->TreeR()->GetEntry(0);
    if(isDig)  GetLoader()->TreeD()->GetEntry(0);  
    
    for(Int_t iChamN=1;iChamN<=7;iChamN++){//chambers loop
      if(isClus){
        Int_t iNclusCham=Clusters(iChamN)->GetEntries(); if(iNclusCham) pNumClusH1->Fill(iNclusCham);//number of clusters per event
        for(Int_t iClusN=0;iClusN<iNclusCham;iClusN++){//clusters loop
          AliRICHcluster *pClus=(AliRICHcluster*)Clusters(iChamN)->At(iClusN);
                                       pQdcH1        ->Fill(pClus->Q());   
                                       pSizeH1       ->Fill(pClus->Size());  
                                       pMapH2        ->Fill(pClus->X(),pClus->Y()); //common
                                       
           if(pClus->IsSingleMip())     {pPureMipQdcH1 ->Fill(pClus->Q());
                                       pPureMipSizeH1->Fill(pClus->Size());
                                       pPureMipMapH2 ->Fill(pClus->X(),pClus->Y());}//Pure Mips
                                       
           if(pClus->IsSingleCerenkov()){pPureCerQdcH1 ->Fill(pClus->Q());
                                       pPureCerSizeH1->Fill(pClus->Size());
                                       pPureCerMapH2 ->Fill(pClus->X(),pClus->Y());}//Pure Cerenkovs
                                       
           if(pClus->IsSingleFeedback()){pPureFeeQdcH1 ->Fill(pClus->Q());
                                       pPureFeeSizeH1->Fill(pClus->Size());
                                       pPureFeeMapH2 ->Fill(pClus->X(),pClus->Y());}//Pure Feedbacks
           
           if(pClus->IsMip()) {pMipQdcH1 ->Fill(pClus->Q());} //MIP+ other contributions
           if(!pClus->IsPureMip())     pPhotQdcH1->Fill(pClus->Q());  //not MIP
        }//clusters loop
      }//isClus
      if(isDig){
        for(Int_t iDigN=0;iDigN<R()->Digits(iChamN)->GetEntries();iDigN++){//digits loop
          AliRICHdigit *pDig=(AliRICHdigit*)R()->Digits(iChamN)->At(iDigN);
          AliRICHhit   *pHit=Hit(pDig->GetTrack(0));//get first hit of this digit
          TVector2 hitV2=R()->C(iChamN)->Mrs2Pc(pHit->OutX3()); TVector2 digV2=R()->P()->Pad2Loc(pDig->Pad());//center of pad for digit
          pHxD->Fill(hitV2.X()-digV2.X()); pHyD->Fill(hitV2.Y()-digV2.Y());
        }//digits loop
      }//isDig
    }//chambers loop
    Info("ControlPlots","Event %i processed.",iEvtN);
  }//events loop 
  
  if(isDig)  GetLoader()->UnloadDigits();
  if(isClus) GetLoader()->UnloadRecPoints();
  
  pFile->Write(); delete pFile;
  sw.Print();time.Print();
}//ControlPlots()
//__________________________________________________________________________________________________
AliRICHhit* AliRICH::Hit(Int_t tid)
{
//defines which hit provided by given tid for the currently loaded event
  R()->GetLoader()->LoadHits();
  for(Int_t iPrimN=0;iPrimN<R()->GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop      
    R()->GetLoader()->TreeH()->GetEntry(iPrimN);
    for(Int_t iHitN=0;iHitN<R()->Hits()->GetEntries();iHitN++){
      AliRICHhit *pHit=(AliRICHhit*)R()->Hits()->At(iHitN);
      if(tid==pHit->Track()) {R()->GetLoader()->UnloadHits();return pHit;}
    }//hits
  }//prims loop
  R()->GetLoader()->UnloadHits();
  return 0;
}
//__________________________________________________________________________________________________
void AliRICH::PrintHits(Int_t iEvtN)
{
//Prints a list of RICH hits for a given event. Default is event number 0.
  AliInfo(Form("List of RICH hits for event %i",iEvtN));
  R()->GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(R()->GetLoader()->LoadHits()) return;
  
  Int_t iTotalHits=0;
  for(Int_t iPrimN=0;iPrimN<R()->GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
    R()->GetLoader()->TreeH()->GetEntry(iPrimN);      
    R()->Hits()->Print();
    iTotalHits+=R()->Hits()->GetEntries();
  }
  R()->GetLoader()->UnloadHits();
  AliInfo(Form("totally %i hits",iTotalHits));
}
//__________________________________________________________________________________________________
void AliRICH::PrintSDigits(Int_t iEvtN)
{
//prints a list of RICH sdigits  for a given event
  Info("PrintSDigits","List of RICH sdigits for event %i",iEvtN);
  R()->GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(R()->GetLoader()->LoadSDigits()) return;
  
  R()->GetLoader()->TreeS()->GetEntry(0);
  R()->SDigits()->Print();
  R()->GetLoader()->UnloadSDigits();
  Info("PrintSDigits","totally %i sdigits",R()->SDigits()->GetEntries());
}
//__________________________________________________________________________________________________
void AliRICH::PrintDigits(Int_t iEvtN)
{
//prints a list of RICH digits  for a given event
  Info("PrintDigits","List of RICH digits for event %i",iEvtN);
  R()->GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(R()->GetLoader()->LoadDigits()) return;
  
  Int_t iTotalDigits=0;
  R()->GetLoader()->TreeD()->GetEntry(0);
  for(Int_t iChamber=1;iChamber<=kNchambers;iChamber++){
    R()->Digits(iChamber)->Print();
    iTotalDigits+=R()->Digits(iChamber)->GetEntries();
  }
  R()->GetLoader()->UnloadDigits();
  Info("PrintDigits","totally %i Digits",iTotalDigits);
}
//__________________________________________________________________________________________________
void AliRICH::PrintClusters(Int_t iEvtN)
{
//prints a list of RICH clusters  for a given event
  Info("PrintClusters","List of RICH clusters for event %i",iEvtN);
  R()->GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(R()->GetLoader()->LoadRecPoints()) return;
  
  Int_t iTotalClusters=0;
  R()->GetLoader()->TreeR()->GetEntry(0);
  for(Int_t iChamber=1;iChamber<=kNchambers;iChamber++){
    R()->Clusters(iChamber)->Print();
    iTotalClusters+=R()->Clusters(iChamber)->GetEntries();
  }
  R()->GetLoader()->UnloadRecPoints();
  Info("PrintClusters","totally %i clusters",iTotalClusters);
}
//__________________________________________________________________________________________________
void AliRICH::PrintTracks(Int_t iEvtN)
{
//prints a list of tracks (including secondaru) for a given event
  Info("PrintTracks","List of all tracks for event %i",iEvtN);
  R()->GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(R()->GetLoader()->GetRunLoader()->LoadHeader()) return;
  if(R()->GetLoader()->GetRunLoader()->LoadKinematics()) return;
  Int_t iTracksCounter=0;
  Info("PrintTracks","totally %i tracks",iTracksCounter);
}
