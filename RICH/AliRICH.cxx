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
 
ClassImp(AliRICHhit)
//__________________________________________________________________________________________________
void AliRICHhit::Print(Option_t*)const
{
  ::Info("hit","Ch=%1i, TID=%6i, eloss=%9.3f eV, in-out dist=%9.4f, OUT(%7.2f,%7.2f,%7.2f)"
      ,fChamber,fTrack,fEloss*1e9,Length(),fOutX3.X(),fOutX3.Y(),fOutX3.Z());
}
//__________________________________________________________________________________________________
ClassImp(AliRICHdigit)
//__________________________________________________________________________________________________
void AliRICHdigit::Print(Option_t*)const
{
  ::Info("digit","cfm=%9i, cs=%2i, x=%3i, y=%3i, q=%8.3f, TID1=%5i, TID2=%5i, TID3=%5i",
                  fCFM,fChamber,fPadX,fPadY,fQdc,fTracks[0],fTracks[1],fTracks[2]);
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
    ::Info("cluster","cfm=%10i, cs=%2i, SiMa=%6i, Shape=%5i, x=%7.3f, y=%7.3f, Q=%6i, %s with %i digits",
                             fCFM,fChamber,fSize,fShape,fX,fY,fQdc,status,0);
    
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
  if(GetDebug())Info("named ctor","Start.");
//AliDetector ctor deals with Hits and Digits (reset them to 0, does not create them)
  CreateHits();          gAlice->GetMCApp()->AddHitList(fHits);
  fCounters.ResizeTo(20); fCounters.Zero();
  if(GetDebug())Info("named ctor","Stop.");
}//AliRICH::AliRICH(const char *name, const char *title)
//__________________________________________________________________________________________________
AliRICH::~AliRICH()
{
//dtor
  if(GetDebug()) Info("dtor","Start.");

  if(fpParam)    delete fpParam;
  
  if(fHits)      delete fHits;
  if(fSdigits)   delete fSdigits;
  if(fDigits)    delete fDigits;
  if(fDigitsNew) {fDigitsNew->Delete();   delete fDigitsNew;}
  if(fClusters)  {fClusters->Delete();    delete fClusters;}
  if(GetDebug()) Info("dtor","Stop.");    
}//AliRICH::~AliRICH()
//__________________________________________________________________________________________________
void AliRICH::Hits2SDigits()
{
// Create a list of sdigits corresponding to list of hits. Every hit generates one or more sdigits.
  if(GetDebug()) Info("Hit2SDigits","Start.");
  for(Int_t iEventN=0;iEventN<GetLoader()->GetRunLoader()->GetAliRun()->GetEventsPerRun();iEventN++){//events loop
    GetLoader()->GetRunLoader()->GetEvent(iEventN);//get next event
  
    if(!GetLoader()->TreeH()) GetLoader()->LoadHits();    GetLoader()->GetRunLoader()->LoadHeader(); 
                                                          GetLoader()->GetRunLoader()->LoadKinematics();//from
    if(!GetLoader()->TreeS()) GetLoader()->MakeTree("S"); MakeBranch("S");//to
          
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      GetLoader()->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){//hits loop 
        AliRICHhit *pHit=(AliRICHhit*)Hits()->At(iHitN);//get current hit                
        TVector2 x2 = C(pHit->C())->Glob2Loc(pHit->OutX3());//hit position in the chamber local system
        Int_t iTotQdc=P()->TotQdc(x2,pHit->Eloss());//total charge produced by hit, 0 if hit in dead zone
        if(iTotQdc==0) continue;
        TVector area=P()->Loc2Area(x2);//determine affected pads, dead zones analysed inside
        if(GetDebug()) Info("Hits2SDigits","hit(%6.2f,%6.2f)->area(%3.0f,%3.0f)-(%3.0f,%3.0f) QDC=%4i",x2.X(),x2.Y(),area[0],area[1],area[2],area[3],iTotQdc);
        TVector pad(2);
        for(pad[1]=area[1];pad[1]<=area[3];pad[1]++)//affected pads loop
          for(pad[0]=area[0];pad[0]<=area[2];pad[0]++){                    
            Double_t padQdc=iTotQdc*P()->FracQdc(x2,pad);
            if(GetDebug()) Info("Hits2SDigits","current pad(%3.0f,%3.0f) with QDC  =%6.2f",pad[0],pad[1],padQdc);
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
  if(GetDebug()) Info("Hit2SDigits","Stop.");
}//Hits2SDigits()
//__________________________________________________________________________________________________
void AliRICH::BuildGeometry() 
{
//Builds a TNode geometry for event display
  if(GetDebug())Info("BuildGeometry","Start.");
  
  TNode *node, *subnode, *top;
  top=gAlice->GetGeometry()->GetNode("alice");
  
  new TBRIK("S_RICH","S_RICH","void",71.09999,11.5,73.15);

  Float_t wid=P()->SectorSizeX();
  Float_t len=P()->SectorSizeY();
  new TBRIK("PHOTO","PHOTO","void",wid/2,0.1,len/2);
  
  for(int i=1;i<=kNchambers;i++){
    top->cd();
    node = new TNode(Form("RICH%i",i),Form("RICH%i",i),"S_RICH",C(i)->X(),C(i)->Y(),C(i)->Z(),C(i)->RotMatrixName());
    node->SetLineColor(kRed);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",wid+P()->DeadZone(),5,len/2+P()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,len/2+P()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-wid-P()->DeadZone(),5,len/2+P()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",wid+P()->DeadZone(),5,-len/2-P()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-len/2 -P()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-wid-P()->DeadZone(),5,-len/2 - P()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);
  }  
  if(GetDebug())Info("BuildGeometry","Stop.");    
}//void AliRICH::BuildGeometry()

//______________________________________________________________________________
void AliRICH::CreateMaterials()
{
// Definition of available RICH materials  
#include "Opticals.h"
        
  Float_t a=0,z=0,den=0,radl=0,absl=0;
  Float_t tmaxfd=-10.0, deemax=-0.2, stemax=-0.1,epsil=0.001, stmin=-0.001; 
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  Int_t material;
    
  AliMaterial(material=1, "Air     $",a=14.61,z=7.3, den=0.001205,radl=30420.0,absl=67500);//(Air)
  AliMedium(1, "DEFAULT MEDIUM AIR$",material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial( 6, "HON",      a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);    //(C)-equivalent radl
  AliMedium(2, "HONEYCOMB$", 6, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(16, "CSI",      a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);    //CsI-radl equivalent
  AliMedium(kCSI, "CSI$", 16, 1, isxfld, sxmgmx,tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(11, "GRI",      a=63.54,z=29.0,den=8.96,    radl=1.43,   absl=0);    //anode grid (Cu) 
  AliMedium(7, "GRID$", 11, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(50, "ALUM",     a=26.98,z=13.0,den=2.699,     radl=8.9,    absl=0);    //aluminium sheet (Al)
  AliMedium(10, "ALUMINUM$", 50, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(material=31, "COPPER$",  a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);    //(Cu)
  AliMedium(12, "PCB_COPPER", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aQuartz[2]={28.09,16.0};  Float_t  zQuartz[2]={14.00, 8.0};  Float_t  wmatQuartz[2]={1,2};
  AliMixture (20, "QUA",aQuartz,zQuartz,den=2.64,-2, wmatQuartz);//Quarz (SiO2) - trasnparent 
  AliMedium(3, "QUARTZ$", 20, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture (21, "QUAO",aQuartz, zQuartz, den=2.64, -2, wmatQuartz);//Quarz (SiO2) - opaque
  AliMedium(8, "QUARTZO$", 21, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aFreon[2]={12,19};  Float_t  zFreon[2]={6,9};  Float_t wmatFreon[2]={6,14};
  AliMixture (material=30, "C6F14",aFreon,zFreon,den=1.68,-2,wmatFreon);//Freon (C6F14) 
  AliMedium(4, "C6F14$",material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t aMethane[2]={12.01,1}; Float_t zMethane[2]={6,1}; Float_t wmatMethane[2]={1,4};
  AliMixture (material=40, "CH4", aMethane, zMethane, den=7.17e-4,-2, wmatMethane);//methane (CH4)     
  AliMedium(kCH4, "CH4$"   , material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);  
  AliMedium(kGAP, "CH4GAP$", material, 1, isxfld, sxmgmx,tmaxfd, 0.1, -deemax, epsil, -stmin);
  
  if(P()->IsRadioSrc()){
    AliMaterial(material=45, "STEEL$",  a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);    //Steel
    AliMedium(kSteel, "STEEL", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
    AliMaterial(material=46, "PERPEX$",  a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);    //Perpex
    AliMedium(kPerpex, "PERPEX", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    
    AliMaterial(material=47, "STRONZIUM$",  a=87.62,z=38.0,den=2.54,    radl=4.24,    absl=0); //Sr90
    AliMedium(kSr90, "STRONZIUM", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  }
  
  Float_t aGlass[5]={12.01, 28.09, 16.,   10.8,  23.};
  Float_t zGlass[5]={ 6.,   14.,    8.,    5.,   11.};
  Float_t wGlass[5]={ 0.5,  0.105, 0.355, 0.03,  0.01};
  AliMixture(material=32, "GLASS",aGlass, zGlass, den=1.74, 5, wGlass);//Glass 50%C+10.5%Si+35.5%O+3% + 1%
  AliMedium(11, "GLASS", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
            
  Int_t *idtmed = fIdtmed->GetArray()-999;
  gMC->SetCerenkov(idtmed[1000], kNbins, aPckov, aAbsCH4,    aQeAll, aIdxCH4);
  gMC->SetCerenkov(idtmed[1001], kNbins, aPckov, aAbsCH4,    aQeAll, aIdxCH4);
  gMC->SetCerenkov(idtmed[1002], kNbins, aPckov, aAbsSiO2,   aQeAll, aIdxSiO2);
  gMC->SetCerenkov(idtmed[1003], kNbins, aPckov, aAbsC6F14,  aQeAll, aIdxC6F14);
  gMC->SetCerenkov(idtmed[1004], kNbins, aPckov, aAbsCH4,    aQeAll, aIdxCH4);
  gMC->SetCerenkov(idtmed[1005], kNbins, aPckov, aAbsCsI,    aQeCsI, aIdxCH4);
  gMC->SetCerenkov(idtmed[1006], kNbins, aPckov, aAbsGrid,   aQeAll, aIdxGrid);
  gMC->SetCerenkov(idtmed[1007], kNbins, aPckov, aAbsOpSiO2, aQeAll, aIdxOpSiO2);
  gMC->SetCerenkov(idtmed[1008], kNbins, aPckov, aAbsCH4,    aQeAll, aIdxCH4);
  gMC->SetCerenkov(idtmed[1009], kNbins, aPckov, aAbsGrid,   aQeAll, aIdxGrid);
  gMC->SetCerenkov(idtmed[1010], kNbins, aPckov, aAbsOpSiO2, aQeAll, aIdxOpSiO2);
    
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
  if(GetDebug())Info("MakeBranch","Start with option= %s.",option);
    
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
  if(GetDebug())Info("MakeBranch","Stop.");   
}//void AliRICH::MakeBranch(Option_t* option)
//__________________________________________________________________________________________________
void AliRICH::SetTreeAddress()
{
//Set branch address for the Hits and Digits Tree.
  if(GetDebug())Info("SetTreeAddress","Start.");
      
  TBranch *branch;
    
  if(fLoader->TreeH()){//H
    if(GetDebug())Info("SetTreeAddress","tree H is requested.");
    CreateHits();//branch map will be in AliDetector::SetTreeAddress    
  }//H
  AliDetector::SetTreeAddress();//this is after TreeH because we need to guarantee that fHits array is created

  if(fLoader->TreeS()){//S
    if(GetDebug())Info("SetTreeAddress","tree S is requested.");
    branch=fLoader->TreeS()->GetBranch(GetName());        if(branch){CreateSDigits();   branch->SetAddress(&fSdigits);}
  }//S
    
  if(fLoader->TreeD()){//D    
    if(GetDebug())Info("SetTreeAddress","tree D is requested.");
    for(int i=0;i<kNchambers;i++){      
      branch=fLoader->TreeD()->GetBranch(Form("%s%d",GetName(),i+1)); 
      if(branch){CreateDigits(); branch->SetAddress(&((*fDigitsNew)[i]));}
    }
  }//D
    
  if(fLoader->TreeR()){//R
    if(GetDebug())Info("SetTreeAddress","tree R is requested.");
    for(int i=0;i<kNchambers;i++){         
      branch=fLoader->TreeR()->GetBranch(Form("%sClusters%d" ,GetName(),i+1));
      if(branch){CreateClusters(); branch->SetAddress(&((*fClusters)[i]));}
    }
  }//R
  if(GetDebug())Info("SetTreeAddress","Stop.");
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
  if(GetDebug())Info("CreateGeometry","Start main.");
//Opaque quartz thickness
  Float_t oquaThickness = .5;
//CsI dimensions
  Float_t pcX=P()->PcSizeX();
  Float_t pcY=P()->PcSizeY();
  
  Int_t *idtmed = fIdtmed->GetArray()-999;
    
  Int_t i;
  Float_t zs;
  Int_t idrotm[1099];
  Float_t par[3];
  par[0]=68.8;par[1]=13   ;par[2]=70.86;  gMC->Gsvolu("RICH","BOX ",(*fIdtmed)[kAl], par,3);//External aluminium box 
  par[0]=66.3;par[1]=13   ;par[2]=68.35;  gMC->Gsvolu("SRIC","BOX ",(*fIdtmed)[kAir],par,3);//Air 
  par[0]=66.3;par[1]=0.025;par[2]=68.35;  gMC->Gsvolu("ALUM","BOX ",(*fIdtmed)[kAl], par,3);//Aluminium sheet 
//Air 2 (cutting the lower part of the box)
//  par[0]=1.25;    par[1] = 3;    par[2] = 70.86;   gMC->Gsvolu("AIR2", "BOX ", idtmed[1000], par, 3);
//Air 3 (cutting the lower part of the box)
//  par[0]=66.3;    par[1] = 3;  par[2] = 1.2505;    gMC->Gsvolu("AIR3", "BOX ", idtmed[1000], par, 3);
//Honeycomb 
  par[0]=66.3;par[1]=0.188;  par[2] = 68.35;       gMC->Gsvolu("HONE", "BOX ", idtmed[1001], par, 3);

  //par[0] = 66.5; par[1] = .025; par[2] = 63.1;
//Quartz 
  par[0]=P()->QuartzWidth()/2;par[1]=P()->QuartzThickness()/2;par[2]=P()->QuartzLength()/2;
  gMC->Gsvolu("QUAR", "BOX ", idtmed[1002], par, 3);
//Spacers (cylinders) 
  par[0]=0.;par[1]=.5;par[2]=P()->FreonThickness()/2;  gMC->Gsvolu("SPAC", "TUBE", idtmed[1002], par, 3);    
//Feet (freon slabs supports)
  par[0] = .7;  par[1] = .3;  par[2] = 1.9;        gMC->Gsvolu("FOOT", "BOX", idtmed[1009], par, 3);
//Opaque quartz 
  par[0]=P()->QuartzWidth()/2;par[1]= .2;par[2]=P()->QuartzLength()/2;
  gMC->Gsvolu("OQUA", "BOX ", idtmed[1007], par, 3);
//Frame of opaque quartz
  par[0]=P()->OuterFreonWidth()/2;par[1]=P()->FreonThickness()/2;par[2]=P()->OuterFreonLength()/2; 
  gMC->Gsvolu("OQF1", "BOX ", idtmed[1007], par, 3);
  par[0]=P()->InnerFreonWidth()/2;par[1]=P()->FreonThickness()/2;par[2]=P()->InnerFreonLength()/2; 
  gMC->Gsvolu("OQF2", "BOX ", idtmed[1007], par, 3);
//Freon 
  par[0]=P()->OuterFreonWidth()/2 - oquaThickness;
  par[1]=P()->FreonThickness()/2;
  par[2]=P()->OuterFreonLength()/2 - 2*oquaThickness; 
  gMC->Gsvolu("FRE1", "BOX ", idtmed[1003], par, 3);

  par[0]=P()->InnerFreonWidth()/2 - oquaThickness;
  par[1]=P()->FreonThickness()/2;
  par[2]=P()->InnerFreonLength()/2 - 2*oquaThickness; 
  gMC->Gsvolu("FRE2", "BOX ", idtmed[1003], par, 3);    
//Methane 
  par[0]=pcX/2;par[1]=P()->GapThickness()/2;par[2]=pcY/2;         gMC->Gsvolu("META","BOX ",idtmed[1004], par, 3);
//Methane gap 
  par[0]=pcX/2;par[1]=P()->GapAmp()/2;par[2]=pcY/2;gMC->Gsvolu("GAP ","BOX ",(*fIdtmed)[kGAP],par,3);
//CsI PC
  par[0]=pcX/2;par[1]=.25;par[2]=pcY/2;  gMC->Gsvolu("CSI ", "BOX ", (*fIdtmed)[kCSI], par, 3);
//Anode grid 
  par[0] = 0.;par[1] = .001;par[2] = 20.;  gMC->Gsvolu("GRID", "TUBE", idtmed[1006], par, 3);

//Wire supports
//Bar of metal
  par[0]=pcX/2;par[1]=1.05;par[2]=1.05;  gMC->Gsvolu("WSMe", "BOX ", idtmed[1009], par, 3);
//Ceramic pick up (base)
  par[0]=pcX/2;par[1]= .25;par[2]=1.05;  gMC->Gsvolu("WSG1", "BOX ", idtmed[1010], par, 3);
//Ceramic pick up (head)
  par[0] = pcX/2;par[1] = .1;par[2] = .1;  gMC->Gsvolu("WSG2", "BOX ", idtmed[1010], par, 3);

//Aluminium supports for methane and CsI
//Short bar
  par[0]=pcX/2;par[1]=P()->GapThickness()/2 + .25; par[2] = (68.35 - pcY/2)/2;
  gMC->Gsvolu("SMSH", "BOX", idtmed[1009], par, 3);
//Long bar
  par[0]=(66.3 - pcX/2)/2;par[1]=P()->GapThickness()/2+.25;par[2]=pcY/2+68.35-pcY/2;
  gMC->Gsvolu("SMLG", "BOX", idtmed[1009], par, 3);
    
//Aluminium supports for freon
//Short bar
  par[0] = P()->QuartzWidth()/2; par[1] = .3; par[2] = (68.35 - P()->QuartzLength()/2)/2;
  gMC->Gsvolu("SFSH", "BOX", idtmed[1009], par, 3);    
//Long bar
  par[0] = (66.3 - P()->QuartzWidth()/2)/2; par[1] = .3;
  par[2] = P()->QuartzLength()/2 + 68.35 - P()->QuartzLength()/2;
  gMC->Gsvolu("SFLG", "BOX", idtmed[1009], par, 3);    
//PCB backplane
  par[0] = pcX/2;par[1] = .25; par[2] = pcY/4 -.5025;  gMC->Gsvolu("PCB ", "BOX", idtmed[1011], par, 3);

//Backplane supports
//Aluminium slab
  par[0] = 33.15;par[1] = 2;par[2] = 21.65;  gMC->Gsvolu("BACK", "BOX", idtmed[1009], par, 3);    
//Big hole
  par[0] = 9.05; par[1] = 2; par[2] = 4.4625;  gMC->Gsvolu("BKHL", "BOX", idtmed[1000], par, 3);
//Small hole
  par[0] = 5.7;par[1] = 2;par[2] = 4.4625;  gMC->Gsvolu("BKHS", "BOX", idtmed[1000], par, 3);
//Place holes inside backplane support
  gMC->Gspos("BKHS", 1, "BACK", .8 + 5.7,0., .6 + 4.4625, 0, "ONLY");
  gMC->Gspos("BKHS", 2, "BACK", -.8 - 5.7,0., .6 + 4.4625, 0, "ONLY");
  gMC->Gspos("BKHS", 3, "BACK", .8 + 5.7,0., -.6 - 4.4625, 0, "ONLY");
  gMC->Gspos("BKHS", 4, "BACK", -.8 - 5.7,0., -.6 - 4.4625, 0, "ONLY");
  gMC->Gspos("BKHS", 5, "BACK", .8 + 5.7,0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
  gMC->Gspos("BKHS", 6, "BACK", -.8 - 5.7,0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
  gMC->Gspos("BKHS", 7, "BACK", .8 + 5.7,0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
  gMC->Gspos("BKHS", 8, "BACK", -.8 - 5.7,0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
  gMC->Gspos("BKHL", 1, "BACK", .8 + 11.4 + 1.6 + 9.05, 0., .6 + 4.4625, 0, "ONLY");
  gMC->Gspos("BKHL", 2, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., .6 + 4.4625, 0, "ONLY");
  gMC->Gspos("BKHL", 3, "BACK", .8 + 11.4 + 1.6 + 9.05, 0., -.6 - 4.4625, 0, "ONLY");
  gMC->Gspos("BKHL", 4, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., -.6 - 4.4625, 0, "ONLY");
  gMC->Gspos("BKHL", 5, "BACK", .8 + 11.4+ 1.6 + 9.05, 0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
  gMC->Gspos("BKHL", 6, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
  gMC->Gspos("BKHL", 7, "BACK", .8 + 11.4 + 1.6 + 9.05, 0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
  gMC->Gspos("BKHL", 8, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
//Place material inside RICH 
  gMC->Gspos("SRIC", 1, "RICH", 0.,0., 0., 0, "ONLY");
//  gMC->Gspos("AIR2", 1, "RICH", 66.3 + 1.2505, 1.276-P()->GapThickness()/2-P()->QuartzThickness()-P()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
//  gMC->Gspos("AIR2", 2, "RICH", -66.3 - 1.2505,1.276-P()->GapThickness()/2-P()->QuartzThickness()-P()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
//  gMC->Gspos("AIR3", 1, "RICH", 0., 1.276-P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, -68.35 - 1.25, 0, "ONLY");
//  gMC->Gspos("AIR3", 2, "RICH", 0., 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35,  68.35 + 1.25, 0, "ONLY");
  gMC->Gspos("ALUM", 1, "SRIC", 0., 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .6 - .05 - .376 -.025, 0., 0, "ONLY");
  gMC->Gspos("HONE", 1, "SRIC", 0., 1.276- P()->GapThickness()/2  - P()->QuartzThickness() - P()->FreonThickness()- .4 - .6 - .05 - .188, 0., 0, "ONLY");
  gMC->Gspos("ALUM", 2, "SRIC", 0., 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .6 - .025, 0., 0, "ONLY");
  gMC->Gspos("FOOT", 1, "SRIC", 64.95, 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .3, 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 2, "SRIC", 21.65, 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .3 , 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 3, "SRIC", -21.65, 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .3, 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 4, "SRIC", -64.95, 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .3, 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 5, "SRIC", 64.95, 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 6, "SRIC", 21.65, 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 7, "SRIC", -21.65, 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 8, "SRIC", -64.95, 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("OQUA", 1, "SRIC", 0., 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .2, 0., 0, "ONLY");
// Methane supports
  gMC->Gspos("SMLG", 1, "SRIC", pcX/2 + (66.3 - pcX/2)/2, 1.276 + .25, 0., 0, "ONLY");
  gMC->Gspos("SMLG", 2, "SRIC", - pcX/2 - (66.3 - pcX/2)/2, 1.276 + .25, 0., 0, "ONLY");
  gMC->Gspos("SMSH", 1, "SRIC", 0., 1.276 + .25, pcY/2 + (68.35 - pcY/2)/2, 0, "ONLY");
  gMC->Gspos("SMSH", 2, "SRIC", 0., 1.276 + .25, - pcY/2 - (68.35 - pcY/2)/2, 0, "ONLY");
//Freon supports
  Float_t suppY = 1.276 - P()->GapThickness()/2- P()->QuartzThickness() -P()->FreonThickness() - .2 + .3; //y position of freon supports
  gMC->Gspos("SFLG", 1, "SRIC", P()->QuartzWidth()/2 + (66.3 - P()->QuartzWidth()/2)/2, suppY, 0., 0, "ONLY");
  gMC->Gspos("SFLG", 2, "SRIC", - P()->QuartzWidth()/2 - (66.3 - P()->QuartzWidth()/2)/2, suppY, 0., 0, "ONLY");
  gMC->Gspos("SFSH", 1, "SRIC", 0., suppY, P()->QuartzLength()/2 + (68.35 - P()->QuartzLength()/2)/2, 0, "ONLY");
  gMC->Gspos("SFSH", 2, "SRIC", 0., suppY, - P()->QuartzLength()/2 - (68.35 - P()->QuartzLength()/2)/2, 0, "ONLY");
  AliMatrix(idrotm[1019], 0., 0., 90., 0., 90., 90.);
//Place spacers
  Int_t nspacers = 30;
  for (i = 0; i < nspacers/3; i++) {
    zs = -11.6/2 + (TMath::Abs(nspacers/6) - i) * 12.2;
    gMC->Gspos("SPAC", i, "FRE1", 10.5, 0., zs, idrotm[1019], "ONLY");  //Original settings 
  }
  for (i = nspacers/3; i < (nspacers*2)/3; i++) {
    zs = -11.6/2 + (nspacers/3 + TMath::Abs(nspacers/6) - i) * 12.2;
    gMC->Gspos("SPAC", i, "FRE1", 0, 0., zs, idrotm[1019], "ONLY");  //Original settings 
  }
  for (i = (nspacers*2)/3; i < nspacers; ++i) {
    zs = -11.6/2 + ((nspacers*2)/3 + TMath::Abs(nspacers/6) - i) * 12.2;
    gMC->Gspos("SPAC", i, "FRE1", -10.5, 0., zs, idrotm[1019], "ONLY"); //Original settings  
  }
  for (i = 0; i < nspacers/3; i++) {
    zs = -11.6/2 + (TMath::Abs(nspacers/6) - i) * 12.2;
    gMC->Gspos("SPAC", i, "FRE2", 10.5, 0., zs, idrotm[1019], "ONLY");  //Original settings 
  }
  for (i = nspacers/3; i < (nspacers*2)/3; i++) {
    zs = -11.6/2 + (nspacers/3 + TMath::Abs(nspacers/6) - i) * 12.2;
    gMC->Gspos("SPAC", i, "FRE2", 0, 0., zs, idrotm[1019], "ONLY");  //Original settings 
  }
  for (i = (nspacers*2)/3; i < nspacers; ++i) {
    zs = -11.6/2 + ((nspacers*2)/3 + TMath::Abs(nspacers/6) - i) * 12.2;
    gMC->Gspos("SPAC", i, "FRE2", -10.5, 0., zs, idrotm[1019], "ONLY"); //Original settings  
  }
  gMC->Gspos("FRE1", 1, "OQF1", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("FRE2", 1, "OQF2", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("OQF1", 1, "SRIC", P()->OuterFreonWidth()/2 + P()->InnerFreonWidth()/2 + 2, 1.276 - P()->GapThickness()/2- P()->QuartzThickness() -P()->FreonThickness()/2, 0., 0, "ONLY"); //Original settings (31.3)
  gMC->Gspos("OQF2", 2, "SRIC", 0., 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()/2, 0., 0, "ONLY");          //Original settings 
  gMC->Gspos("OQF1", 3, "SRIC", - (P()->OuterFreonWidth()/2 + P()->InnerFreonWidth()/2) - 2, 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()/2, 0., 0, "ONLY");       //Original settings (-31.3)
  gMC->Gspos("QUAR", 1, "SRIC", 0., 1.276 - P()->GapThickness()/2 - P()->QuartzThickness()/2, 0., 0, "ONLY");
  gMC->Gspos("GAP ", 1, "META", 0., P()->GapThickness()/2 - P()->GapAmp()/2 - 0.0001, 0., 0, "ONLY");
  gMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
  gMC->Gspos("CSI ", 1, "SRIC", 0., 1.276 + P()->GapThickness()/2 + .25, 0., 0, "ONLY");
//Wire support placing
  gMC->Gspos("WSG2", 1, "GAP ", 0., P()->GapAmp()/2 - .1, 0., 0, "ONLY");
  gMC->Gspos("WSG1", 1, "CSI ", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("WSMe", 1, "SRIC ", 0., 1.276 + P()->GapThickness()/2 + .5 + 1.05, 0., 0, "ONLY");
//Backplane placing
  gMC->Gspos("BACK", 1, "SRIC ", -33.15, 1.276 + P()->GapThickness()/2 + .5 + 2.1 + 2, 43.3, 0, "ONLY");
  gMC->Gspos("BACK", 2, "SRIC ", 33.15, 1.276 + P()->GapThickness()/2 + .5 + 2.1 + 2 , 43.3, 0, "ONLY");
  gMC->Gspos("BACK", 3, "SRIC ", -33.15, 1.276 + P()->GapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
  gMC->Gspos("BACK", 4, "SRIC ", 33.15, 1.276 + P()->GapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
  gMC->Gspos("BACK", 5, "SRIC ", 33.15, 1.276 + P()->GapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");
  gMC->Gspos("BACK", 6, "SRIC ", -33.15, 1.276 + P()->GapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");
//PCB placing
  gMC->Gspos("PCB ", 1, "SRIC ", 0.,  1.276 + P()->GapThickness()/2 + .5 + 1.05, pcX/4 + .5025 + 2.5, 0, "ONLY");
  gMC->Gspos("PCB ", 2, "SRIC ", 0.,  1.276 + P()->GapThickness()/2 + .5 + 1.05, -pcX/4 - .5025 - 2.5, 0, "ONLY");
  
//place chambers into mother volume ALIC
  for(int i=1;i<=kNchambers;i++){
    AliMatrix(idrotm[1000+i],C(i)->ThetaXd(),C(i)->PhiXd(),
                             C(i)->ThetaYd(),C(i)->PhiYd(),
                             C(i)->ThetaZd(),C(i)->PhiZd());
    gMC->Gspos("RICH",i,"ALIC",C(i)->X(),C(i)->Y(),C(i)->Z(),idrotm[1000+i], "ONLY");
  }

         
  if(GetDebug())Info("CreateGeometry","Stop main.");  
}//void AliRICH::CreateGeometry()
//__________________________________________________________________________________________________
void AliRICH::Reconstruct()const
{
  AliRICHClusterFinder finder(const_cast<AliRICH*>(this));
  finder.Exec();
}
//__________________________________________________________________________________________________
void AliRICH::ControlPlots()
{ 
//Creates a set of hists to control the results of simulation
     
  TH1F *pHxD=0,*pHyD=0,*pNumClusH1=0,
                   *pQdcH1=0,       *pSizeH1=0,
                   *pPureMipQdcH1=0,*pPureMipSizeH1=0,
                   *pPureCerQdcH1=0,*pPureCerSizeH1=0,
                   *pPureFeeQdcH1=0,*pPureFeeSizeH1=0,
                   *pMipQdcH1=0,    *pPhotQdcH1=0;  
  TH2F *pMapH2=0,*pPureMipMapH2=0,*pPureCerMapH2=0,*pPureFeeMapH2=0;
  
  Bool_t isDig =!GetLoader()->LoadDigits();
  Bool_t isClus=!GetLoader()->LoadRecPoints();

  if(!isDig && !isClus){Error("ControlPlots","No digits and clusters! Nothing to do.");return;}
  
  TStopwatch sw;TDatime time;
    
  TFile *pFile = new TFile("$(HOME)/RCP.root","RECREATE");   
  
  if(isDig){
    cout<<"Digits available\n";
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
          AliRICHhit   *pHit=Hit(pDig->GetTrack(0));
          TVector2 hitV2=R()->C(iChamN)->Glob2Loc(pHit->OutX3()); TVector2 digV2=R()->P()->Pad2Loc(pDig->Pad());
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
  Info("PrintHits","List of RICH hits for event %i",iEvtN);
  R()->GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(R()->GetLoader()->LoadHits()) return;
  
  Int_t iTotalHits=0;
  for(Int_t iPrimN=0;iPrimN<R()->GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
    R()->GetLoader()->TreeH()->GetEntry(iPrimN);      
    R()->Hits()->Print();
    iTotalHits+=R()->Hits()->GetEntries();
  }
  R()->GetLoader()->UnloadHits();
  Info("PrintHits","totally %i hits",iTotalHits);
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
void AliRICH::FillESD(AliESD *pESD)const
{
//This methode fills AliESDtrack with information from RICH  
  Info("FillESD","Start with %i tracks",pESD->GetNumberOfTracks());
  const Double_t masses[5]={0.000511,0.105658,0.139567,0.493677,0.93828};//electron,muon,pion,kaon,proton
  const Double_t refIndex = 1.29052;

  Double_t thetaExp = 0.7;
  Double_t thetaTh[5];
  Double_t sinThetaThNorm;
  Double_t sigmaThetaTh[5];
  Double_t height[5];
  Double_t totalHeight=0;
  
  for(Int_t iTrackN=0;iTrackN<pESD->GetNumberOfTracks();iTrackN++){//ESD tracks loop
    AliESDtrack *pTrack = pESD->GetTrack(iTrackN);
    
    pTrack->Print("");
//    TVector2 x2=P()->HelixCross(pTrack);//returns cross point of track with RICH PC in LRS
    Double_t pmod = pTrack->GetP();
    
    for(Int_t iPart=4;iPart>=0;iPart--){
      Double_t cosThetaTh = TMath::Sqrt(masses[iPart]*masses[iPart]+pmod*pmod)/(refIndex*pmod);
      if(cosThetaTh>=1) {pTrack->ResetRICH(); break;}
      thetaTh[iPart] = TMath::ACos(cosThetaTh);
      sinThetaThNorm = TMath::Sin(thetaTh[iPart])/TMath::Sqrt(1-1/(refIndex*refIndex));
      sigmaThetaTh[iPart] = (0.014*(1/sinThetaThNorm-1) + 0.0043)*1.25;
      height[iPart] = TMath::Gaus(thetaExp,thetaTh[iPart],sigmaThetaTh[iPart]);
      totalHeight +=height[iPart];
    }
    
    pTrack->SetRICHsignal(thetaExp);
    Double_t richPID[5];
    for(Int_t iPart=0;iPart<5;iPart++) richPID[iPart] = height[iPart]/totalHeight;    
    pTrack->SetRICHpid(richPID); 
  }//ESD tracks loop
}
