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
#include <TVirtualMC.h>
 
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
  ::Info("digit","csxy=%6i, cfm=%9i, c=%2i, x=%3i, y=%3i, q=%8.3f, TID1=%5i, TID2=%5i, TID3=%5i",
                  Id(),fChFbMip,fChamber,fPadX,fPadY,fQdc,fTracks[0],fTracks[1],fTracks[2]);
}
//__________________________________________________________________________________________________
ClassImp(AliRICHcluster)
//__________________________________________________________________________________________________
void AliRICHcluster::Print(Option_t*)const
{
  ::Info("cluster","CombiPid=%10i, c=%2i, size=%6i, dim=%5i, x=%7.3f, y=%7.3f, Q=%6i, st=%i",
           fCombiPid,fChamber,fSize,fDimXY,fX,fY,fQdc,fStatus);
}
//__________________________________________________________________________________________________
ClassImp(AliRICHreco)
//__________________________________________________________________________________________________
void AliRICHreco::Print(Option_t*)const
{
  ::Info("reco","ThetaCherenkov=%9.6f, Nphotons=%4i, TID=%9i",fThetaCherenkov,fNphotons,fTid);
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
AliRICH::AliRICH()
        :AliDetector() 
{
//Default ctor should not contain any new operators
  fpParam     =0;
  fChambers   =0;   
//AliDetector ctor deals with Hits and Digits  
  fSdigits    =0; fNsdigits   =0;
  fDigitsNew  =0; for(int i=0;i<kNCH;i++) fNdigitsNew[i]  =0;
  fClusters   =0; for(int i=0;i<kNCH;i++) fNclusters[i]=0;
  fRecos      =0; fNrecos     =0;
}//AliRICH::AliRICH()
//__________________________________________________________________________________________________
AliRICH::AliRICH(const char *name, const char *title)
        :AliDetector(name,title)
{
//Named ctor
  if(GetDebug())Info("named ctor","Start.");
  fpParam     =   new AliRICHParam;
  fChambers = 0;  CreateChambers();
//AliDetector ctor deals with Hits and Digits (reset them to 0, does not create them)
  fHits=       0;     CreateHits();          gAlice->GetMCApp()->AddHitList(fHits);
  fSdigits=    0;
  fDigitsNew=  0;
  fClusters=   0;
  fRecos      =0;
  if(GetDebug())Info("named ctor","Stop.");
}//AliRICH::AliRICH(const char *name, const char *title)
//__________________________________________________________________________________________________
AliRICH::~AliRICH()
{
//dtor
  if(GetDebug()) Info("dtor","Start.");

  if(fpParam)    delete fpParam;
  if(fChambers)  delete fChambers;
  
  if(fHits)      delete fHits;
  if(fSdigits)   delete fSdigits;
  if(fDigits)    delete fDigits;
  if(fDigitsNew) {fDigitsNew->Delete();   delete fDigitsNew;}
  if(fClusters)  {fClusters->Delete();    delete fClusters;}
  if(fRecos)     delete fRecos;
  if(GetDebug()) Info("dtor","Stop.");    
}//AliRICH::~AliRICH()
//__________________________________________________________________________________________________
void AliRICH::Hits2SDigits()
{
// Create a list of sdigits corresponding to list of hits. Every hit generates one or more sdigits.
//   
  if(GetDebug()) Info("Hit2SDigits","Start.");

  AliLoader * richLoader = GetLoader();
  AliRunLoader * runLoader = GetLoader()->GetRunLoader();

  for(Int_t iEventN=0;iEventN<GetLoader()->GetRunLoader()->GetAliRun()->GetEventsPerRun();iEventN++){//events loop
    runLoader->GetEvent(iEventN);
  
    if (!richLoader->TreeH()) richLoader->LoadHits();
    if (!runLoader->TreeE()) runLoader->LoadHeader(); 
    if (!runLoader->TreeK()) runLoader->LoadKinematics();//from
    if (!richLoader->TreeS()) richLoader->MakeTree("S"); MakeBranch("S");//to
          
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      richLoader->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){//hits loop 
        AliRICHhit *pHit=(AliRICHhit*)Hits()->At(iHitN);                
        TVector2 x2 = P()->ShiftToWirePos(C(pHit->C())->Glob2Loc(pHit->OutX3()));                
        Int_t iTotQdc=P()->TotQdc(x2,pHit->Eloss());
        if(iTotQdc==0) continue;
        Int_t iPadXmin,iPadXmax,iPadYmin,iPadYmax;
        P()->Loc2Area(x2,iPadXmin,iPadYmin,iPadXmax,iPadYmax);//determine affected pads
        if(GetDebug()) Info("Hits2SDigits","left-down=(%i,%i) right-up=(%i,%i)",iPadXmin,iPadYmin,iPadXmax,iPadYmax);
        for(Int_t iPadY=iPadYmin;iPadY<=iPadYmax;iPadY++)//affected pads loop
          for(Int_t iPadX=iPadXmin;iPadX<=iPadXmax;iPadX++){
            Double_t padQdc=iTotQdc*P()->FracQdc(x2,iPadX,iPadY);
            if(padQdc>0.1) AddSDigit(pHit->C(),iPadX,iPadY,padQdc,
              runLoader->Stack()->Particle(pHit->GetTrack())->GetPdgCode(),pHit->GetTrack());
          }//affected pads loop 
      }//hits loop
    }//prims loop
    richLoader->TreeS()->Fill();
    richLoader->WriteSDigits("OVERWRITE");
    ResetSDigits();
  }//events loop  
  richLoader->UnloadHits();
  runLoader->UnloadHeader();
  runLoader->UnloadKinematics();
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
  
  for(int i=1;i<=kNCH;i++){
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
    //
    // *** DEFINITION OF AVAILABLE RICH MATERIALS *** 
  
#include "Opticals.h"
        
  Float_t a=0,z=0,den=0,radl=0,absl=0;
  Float_t tmaxfd=-10.0, deemax=-0.2, stemax=-0.1,epsil=0.001, stmin=-0.001; 
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
    
  AliMaterial( 1, "Air     $",a=14.61,z=7.3, den=0.001205,radl=30420.0,absl=67500);//(Air)
  AliMedium(1, "DEFAULT MEDIUM AIR$", 1, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial( 6, "HON",      a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);    //(C)-equivalent radl
  AliMedium(2, "HONEYCOMB$", 6, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(16, "CSI",      a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);    //CsI-radl equivalent
  AliMedium(kCSI, "CSI$", 16, 1, isxfld, sxmgmx,tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(11, "GRI",      a=63.54,z=29.0,den=8.96,    radl=1.43,   absl=0);    //anode grid (Cu) 
  AliMedium(7, "GRID$", 11, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(50, "ALUM",     a=26.98,z=13.0,den=2.7,     radl=8.9,    absl=0);    //aluminium sheet (Al)
  AliMedium(10, "ALUMINUM$", 50, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(31, "COPPER$",  a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);    //(Cu)
  AliMedium(12, "PCB_COPPER", 31, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aQuartz[2]={28.09,16.0};  Float_t  zQuartz[2]={14.00, 8.0};  Float_t  wmatQuartz[2]={1,2};
  AliMixture (20, "QUA",aQuartz,zQuartz,den=2.64,-2, wmatQuartz);//Quarz (SiO2) - trasnparent 
  AliMedium(3, "QUARTZ$", 20, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture (21, "QUAO",aQuartz, zQuartz, den=2.64, -2, wmatQuartz);//Quarz (SiO2) - opaque
  AliMedium(8, "QUARTZO$", 21, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aFreon[2]={12,19};  Float_t  zFreon[2]={6,9};  Float_t wmatFreon[2]={6,14};
  AliMixture (30, "FRE",aFreon,zFreon,den=1.7,-2,wmatFreon);//Freon (C6F14) 
  AliMedium(4, "FREON$", 30, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t aMethane[2]={12.01,1}; Float_t zMethane[2]={6,1}; Float_t wmatMethane[2]={1,4};
  AliMixture (40, "MET", aMethane, zMethane, den=7.17e-4,-2, wmatMethane);//methane (CH4)     
  AliMedium(5, "METHANE$", 40, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture (41, "METG", aMethane, zMethane, den=7.17e-4, -2, wmatMethane);
  AliMedium(kGAP, "GAP$", 41, 1, isxfld, sxmgmx,tmaxfd, 0.1, -deemax, epsil, -stmin);
  
  Float_t aGlass[5]={12.01, 28.09, 16.,   10.8,  23.};
  Float_t zGlass[5]={ 6.,   14.,    8.,    5.,   11.};
  Float_t wGlass[5]={ 0.5,  0.105, 0.355, 0.03,  0.01};
  AliMixture (32, "GLASS",aGlass, zGlass, den=1.74, 5, wGlass);//Glass 50%C+10.5%Si+35.5%O+3% + 1%
  AliMedium(11, "GLASS", 32, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
            
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
    for(Int_t i=0;i<kNCH;i++){ 
      MakeBranchInTree(fLoader->TreeD(),Form("%s%d",GetName(),i+1),&((*fDigitsNew)[i]),kBufferSize,0);
    }
  }//D
  
  if(cR&&fLoader->TreeR()){//R
    CreateClusters();
    for(Int_t i=0;i<kNCH;i++)
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
    for(int i=0;i<kNCH;i++){      
      branch=fLoader->TreeD()->GetBranch(Form("%s%d",GetName(),i+1)); 
      if(branch){CreateDigits(); branch->SetAddress(&((*fDigitsNew)[i]));}
    }
  }//D
    
  if(fLoader->TreeR()){//R
    if(GetDebug())Info("SetTreeAddress","tree R is requested.");
    for(int i=0;i<kNCH;i++){         
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
  P()->Dump();
  fChambers->Print(option);  
}//void AliRICH::Print(Option_t *option)const
//__________________________________________________________________________________________________
void AliRICH::CreateGeometry()
{
//Creates detailed geometry simulation (currently GEANT volumes tree)         
  if(GetDebug())Info("CreateGeometry","Start.");
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
    
//External aluminium box 
  par[0]=68.8;par[1]=13;par[2]=70.86;  gMC->Gsvolu("RICH", "BOX ", idtmed[1009], par, 3);
//Air 
  par[0]=66.3;   par[1] = 13; par[2] = 68.35;      gMC->Gsvolu("SRIC", "BOX ", idtmed[1000], par, 3); 
//Air 2 (cutting the lower part of the box)
  par[0]=1.25;    par[1] = 3;    par[2] = 70.86;   gMC->Gsvolu("AIR2", "BOX ", idtmed[1000], par, 3);
//Air 3 (cutting the lower part of the box)
  par[0]=66.3;    par[1] = 3;  par[2] = 1.2505;    gMC->Gsvolu("AIR3", "BOX ", idtmed[1000], par, 3);
//Honeycomb 
  par[0]=66.3;par[1]=0.188;  par[2] = 68.35;       gMC->Gsvolu("HONE", "BOX ", idtmed[1001], par, 3);
//Aluminium sheet 
  par[0]=66.3;par[1]=0.025;par[2]=68.35;           gMC->Gsvolu("ALUM", "BOX ", idtmed[1009], par, 3);
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
  par[0]=pcX/2;par[1]=P()->ProximityGap()/2;par[2]=pcY/2;gMC->Gsvolu("GAP ","BOX ",(*fIdtmed)[kGAP],par,3);
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
  gMC->Gspos("AIR2", 1, "RICH", 66.3 + 1.2505, 1.276-P()->GapThickness()/2-P()->QuartzThickness()-P()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
  gMC->Gspos("AIR2", 2, "RICH", -66.3 - 1.2505,1.276-P()->GapThickness()/2-P()->QuartzThickness()-P()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
  gMC->Gspos("AIR3", 1, "RICH", 0., 1.276-P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, -68.35 - 1.25, 0, "ONLY");
  gMC->Gspos("AIR3", 2, "RICH", 0., 1.276 - P()->GapThickness()/2 - P()->QuartzThickness() - P()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35,  68.35 + 1.25, 0, "ONLY");
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
  gMC->Gspos("GAP ", 1, "META", 0., P()->GapThickness()/2 - P()->ProximityGap()/2 - 0.0001, 0., 0, "ONLY");
  gMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
  gMC->Gspos("CSI ", 1, "SRIC", 0., 1.276 + P()->GapThickness()/2 + .25, 0., 0, "ONLY");
//Wire support placing
  gMC->Gspos("WSG2", 1, "GAP ", 0., P()->ProximityGap()/2 - .1, 0., 0, "ONLY");
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
  for(int i=1;i<=kNCH;i++){
    AliMatrix(idrotm[1000+i],C(i)->ThetaXd(),C(i)->PhiXd(),
                             C(i)->ThetaYd(),C(i)->PhiYd(),
                             C(i)->ThetaZd(),C(i)->PhiZd());
    gMC->Gspos("RICH",i,"ALIC",C(i)->X(),C(i)->Y(),C(i)->Z(),idrotm[1000+i], "ONLY");
  }

  if(GetDebug())Info("CreateGeometry","Stop.");  
}//void AliRICH::CreateGeometry()
//__________________________________________________________________________________________________
void AliRICH::CreateChambers()
{
//create all RICH Chambers on each call. Previous chambers deleted
  if(fChambers) delete fChambers;
  if(GetDebug())Info("CreateChambers","Creating RICH chambers.");
  fChambers=new TObjArray(kNCH);
  fChambers->SetOwner();
  for(int i=0;i<kNCH;i++)  fChambers->AddAt(new AliRICHChamber(i+1,P()),i);  
}//void AliRICH::CreateChambers()
//__________________________________________________________________________________________________
void AliRICH::GenerateFeedbacks(Int_t iChamber,Float_t eloss)
{
// Generate FeedBack photons 
// eloss=0 means photon so only pulse height distribution is to be analysed. This one is done in AliRICHParam::TotQdc()
  
  TLorentzVector x4;
  gMC->TrackPosition(x4);  
  TVector2 x2=C(iChamber)->Glob2Loc(x4);
  Int_t sector=P()->Sector(x2);  if(sector==kBad) return; //hit in dead zone nothing to produce
  Int_t iTotQdc=P()->TotQdc(x2,eloss);
  Int_t iNphotons=gMC->GetRandom()->Poisson(P()->AlphaFeedback(sector)*iTotQdc);    
  if(GetDebug())Info("GenerateFeedbacks","N photons=%i",iNphotons);
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
    gAlice->GetMCApp()->PushTrack(1,                 //do not transport
                     gAlice->GetMCApp()->GetCurrentTrackNumber(),//parent track 
                     kFeedback,                      //PID
		     mom[0],mom[1],mom[2],mom[3],    //track momentum  
                     x4.X(),x4.Y(),x4.Z(),x4.T(),    //track origin 
                     pol[0],pol[1],pol[2],           //polarization
		     kPFeedBackPhoton,
                     outputNtracksStored,
                     1.0);    
  }//feedbacks loop
}//GenerateFeedbacks()

