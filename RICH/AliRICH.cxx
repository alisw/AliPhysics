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
#include "AliRICHTracker.h"
#include "AliRICHHelix.h"
//#include <TArrayF.h>
#include <TGeometry.h>
#include <TBRIK.h>
#include <TTUBE.h>
#include <TFile.h>
#include <TNode.h> 
#include <TObjArray.h>
#include <TParticle.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliGenHijingEventHeader.h>
#include <AliMagF.h>
#include <AliRun.h>
#include <AliRunDigitizer.h>
#include <AliMC.h>
#include <AliESD.h>
#include <TVirtualMC.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TBenchmark.h>
#include <AliLog.h>
#include <TNtupleD.h>
#include <AliTracker.h>
#include <AliRawDataHeader.h>
#include <TLatex.h> //Display()
#include <TCanvas.h> //Display()
#include <TGraph.h> //Display()
#include <TStyle.h> //Display()
#include <TMarker.h> //Display()
 
ClassImp(AliRICH)    
//__________________________________________________________________________________________________
// RICH manager class   
//BEGIN_HTML
/*
  <img src="gif/alirich.gif">
*/
//END_HTML
//__________________________________________________________________________________________________
AliRICH::AliRICH():AliDetector(),fParam(0),  fSdigits(0),fNsdigits(0),fDigs(0),fClus(0) 
{
//Default ctor should not contain any new operators
//AliDetector ctor deals with Hits and Digits  
  for(int i=0;i<kNchambers;i++) fNdigs[i]  =0;
  for(int i=0;i<kNchambers;i++) fNclus[i]=0;
//  fCounters.ResizeTo(20); fCounters.Zero();
}//AliRICH::AliRICH()
//__________________________________________________________________________________________________
AliRICH::AliRICH(const char *name, const char *title)
        :AliDetector(name,title),fParam(new AliRICHParam),fSdigits(0),fNsdigits(0),fDigs(0),fClus(0)
{
//Named ctor
  AliDebug(1,"Start.");
//AliDetector ctor deals with Hits and Digits (reset them to 0, does not create them)
  HitsCreate();          gAlice->GetMCApp()->AddHitList(fHits);
  fCounters.ResizeTo(20); fCounters.Zero();
  AliDebug(1,"Stop.");
}//AliRICH::AliRICH(const char *name, const char *title)
//__________________________________________________________________________________________________
AliRICH::~AliRICH()
{
//dtor
  AliDebug(1,"Start.");

  if(fParam)     delete fParam;
  
  if(fHits)      delete fHits;
  if(fSdigits)   delete fSdigits;
  if(fDigits)    delete fDigits;
  if(fDigs)      {fDigs->Delete();   delete fDigs;}
  if(fClus)      {fClus->Delete();   delete fClus;}
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
    if(!GetLoader()->GetRunLoader()->TreeK())             GetLoader()->GetRunLoader()->LoadKinematics();//from
    if(!GetLoader()->TreeS()) GetLoader()->MakeTree("S"); MakeBranch("S");//to
          
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      GetLoader()->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){//hits loop 
        AliRICHHit *pHit=(AliRICHHit*)Hits()->At(iHitN);//get current hit                
        TVector2 x2 = C(pHit->C())->Mrs2Anod(0.5*(pHit->InX3()+pHit->OutX3()));//hit position in the anod plane
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
            if(padQdc>0.1) SDigitAdd(pHit->C(),pad,padQdc,GetLoader()->GetRunLoader()->Stack()->Particle(pHit->GetTrack())->GetPdgCode(),pHit->GetTrack());
          }//affected pads loop 
      }//hits loop
    }//prims loop
    GetLoader()->TreeS()->Fill();
    GetLoader()->WriteSDigits("OVERWRITE");
    SDigitsReset();
  }//events loop  
  GetLoader()->UnloadHits(); GetLoader()->GetRunLoader()->UnloadHeader(); GetLoader()->GetRunLoader()->UnloadKinematics();
  GetLoader()->UnloadSDigits();  
  AliDebug(1,"Stop.");
}//Hits2SDigits()
//__________________________________________________________________________________________________
void AliRICH::Digits2Raw()
{
//Loops over all digits and creates raw data files in DDL format. GetEvent() is done outside (AliSimulation)
//RICH has 2 DDL per chamber, even number for right part(2-4-6) odd number for left part(1-3-5) 
//RICH has no any propriate header just uses the common one
//Each PC is divided by 8 rows counted from 1 to 8 from top to bottom for left PCs(1-3-5) and from bottom to top for right PCc(2-4-6)     (denoted  rrrrr 5 bits 32 values)
//Each raw is composed from 10 DILOGIC chips counted from left to right from 1 to 10                                                      (denoted   dddd 4 bits 16 values)
//Each DILOGIC chip serves 48 channels counted from 0 to 47                                                                               (denoted aaaaaa 6 bits 64 values)
//So RICH info word is  32 bits word with structure:   0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq (five 0 five r six a twelve q) with QDC    (denoted q...q 12 bits 4096 values)
  AliDebug(1,"Start.");
  GetLoader()->LoadDigits();
  GetLoader()->TreeD()->GetEntry(0);
  
  Int_t kRichOffset=0x700; //currently one DDL per 3 sectors
  
  ofstream file135,file246;//output streams 2 DDL per chamber
  AliRawDataHeader header;//empty DDL miniheader
  UInt_t word32=1;        //32 bits data word 
  
  for(Int_t iChN=1;iChN<=kNchambers;iChN++){ //2 DDL per chamber open both in parallel   
    file135.open(Form("RICH_%4i.ddl",kRichOffset+2*iChN-1));   //left part of chamber; sectors 1-3-5 odd DDL number
    file246.open(Form("RICH_%4i.ddl",kRichOffset+2*iChN));     //right part of chamber; sectors 2-4-6 even DDL number
//common DDL header defined by standart, now dummy as the total number of bytes is not yet known    
    file135.write((char*)&header,sizeof(header)); //dummy header just place holder
    file246.write((char*)&header,sizeof(header)); //actual will be written later
    
    Int_t counter135=0,counter246=0;//counts total number of records per DDL 
    
    for(Int_t iDigN=0;iDigN<Digits(iChN)->GetEntriesFast();iDigN++){//digits loop for a given chamber
      AliRICHDigit *pDig=(AliRICHDigit*)Digits(iChN)->At(iDigN);
      word32=UInt_t (pDig->Q()+0x400*pDig->X()+0x4000*pDig->Y());  //arbitrary structure
      switch(pDig->S()){//readout by vertical sectors: 1,3,5-left DDL 2,4,6-right DDL
        case 1: case 3: case 5: file135.write((char*)&word32,sizeof(word32)); counter135++; break;
        case 2: case 4: case 6: file246.write((char*)&word32,sizeof(word32)); counter246++; break;
      }//switch sector  
    }//digits loop for a given chamber
//now count total byte number for each DDL file and rewrite actual header    
    header.fSize=sizeof(header)+counter135*sizeof(word32);   header.SetAttribute(0); file135.seekp(0); file135.write((char*)&header,sizeof(header));
    header.fSize=sizeof(header)+counter246*sizeof(word32);   header.SetAttribute(0); file246.seekp(0); file246.write((char*)&header,sizeof(header));
    file135.close(); file246.close();
  }//chambers loop  
  GetLoader()->UnloadDigits();
  AliDebug(1,"Stop.");      
}//Digits2Raw()
//__________________________________________________________________________________________________
void AliRICH::BuildGeometry() 
{
//Builds a TNode geometry for event display
  AliInfo("Start.");
  
  TNode *node, *subnode, *top;
  top=gAlice->GetGeometry()->GetNode("alice");

  Float_t widx =P()->SectorSizeX();
  Float_t leny =P()->SectorSizeY();
  Float_t dz   =P()->Zfreon()+P()->Zwin()+P()->Pc2Win();
  Float_t dead =P()->DeadZone();

  new TBRIK("RICH","RICH","void",widx+dead/2,leny+leny/2+dead,dz+0.1); //RICH chamber
  new TBRIK("RPC" ,"RPC" ,"void",widx/2,leny/2,0.01);                  //RICH sector 

  for(int i=1;i<=P()->Nchambers();i++){
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
//__________________________________________________________________________________________________
void AliRICH::CreateMaterials()
{
// Definition of available RICH materials  
        
  Int_t   material=0; //tmp material id number
  Float_t a=0,z=0,den=0,radl=0,absl=0; //tmp material parameters
  
  Float_t tmaxfd=-10.0, deemax=-0.2, stemax=-0.1,epsil=0.001, stmin=-0.001; 
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
    
  Float_t aAir[4]={12,14,16,36};  Float_t zAir[4]={6,7,8,18}; Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};//total 0.9999999
  AliMixture(++material, "RichAir",aAir,zAir,den=0.00120479,4,wAir);                                          //1 (Air) 0.01% C 75% N  23% O 1% Ar
  AliMedium(kAir, "RichAir",material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++material, "RichAerogel",aAir,zAir,den=P()->DenGel(),4,wAir);                     //Aerogel represented by Air
  AliMedium(kGel, "RichAerogel",material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++material, "RichAerogelReflector",aAir,zAir,den=P()->DenGel(),4,wAir);           //Aerogel reflector represented by Air
  AliMedium(kReflector, "RichAerogelReflector",material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++material, "RichRohacell", a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);                   //2 Rohacell 51 C-equiv radl rad cover
  AliMedium(kRoha, "RichRohacell", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aQuartz[2]={28.09,16.0};  Float_t  zQuartz[2]={14.00, 8.0};  Float_t  wQuartz[2]={1,2};
  AliMixture(++material, "RichSiO2",aQuartz,zQuartz,den=2.64,-2, wQuartz);                                    //3 Quarz (SiO2) -trasparent rad window
  AliMedium(kSiO2, "RichSiO2",material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aFreon[2]={12,19};  Float_t  zFreon[2]={6,9};  Float_t wmatFreon[2]={6,14};                        // C12-6 F19-9   
  AliMixture(++material, "RichC6F14",aFreon,zFreon,den=1.68,-2,wmatFreon);                                    //4 Freon (C6F14) 
  AliMedium(kC6F14, "RichC6F14",material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t aMethane[2]={12.01,1}; Float_t zMethane[2]={6,1}; Float_t wMethane[2]={1,4};
  AliMixture (++material, "RichCH4", aMethane, zMethane, den=7.17e-4,-2, wMethane);                        //5,9 methane (CH4) normal and for Gap    
  AliMedium(kCH4, "RichCH4"   , material, 1, isxfld, sxmgmx, tmaxfd, stemax,  deemax, epsil,  stmin);  
  AliMixture (++material, "RichCH4gap", aMethane, zMethane, den=7.17e-4,-2, wMethane);                      //5,9 methane (CH4) normal and for Gap    
  AliMedium(kGap, "RichCH4gap", material, 1, isxfld, sxmgmx, tmaxfd, 0.1   , -deemax, epsil, -stmin);
    
  AliMaterial(++material, "RichCsI",      a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);                   //6 CsI-radl equivalent
  AliMedium(kCsI, "RichCsI", material, 1, isxfld, sxmgmx,tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++material, "RichGridCu",    a=63.54,z=29.0,den=8.96,    radl=1.43,   absl=0);                   //7 anode grid (Cu) 
  AliMedium(kGridCu, "RichGridCu", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    
  AliMaterial(++material, "RichPcbCu",     a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);                   //12 Cu
  AliMedium(kCu, "RichPcbCu", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture (++material, "RichOpSiO2",aQuartz, zQuartz, den=2.64, -2, wQuartz);                             //8 Quarz (SiO2) - opaque
  AliMedium(kOpSiO2, "RichOpSiO2",material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++material, "RichAl",     a=26.98,z=13.0,den=2.699,     radl=8.9,    absl=0);                 //10 aluminium sheet (Al)
  AliMedium(kAl, "RichAl", material, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t aGlass[5]={12.01,28.09,16,10.8,23}; Float_t zGlass[5]={6,14,8,5,11};  Float_t wGlass[5]={0.5,0.105,0.355,0.03,0.01};
  AliMixture(++material,"RichGlass",aGlass, zGlass, den=1.74, 5, wGlass);                                    //11 Glass 50%-C 10.5%-Si 35.5%-O 3%-B 1%-Na
  AliMedium(kGlass, "RichGlass", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  den=19.3;
  AliMaterial(++material, "RichW",  a=183.84,z=74.0,den,    radl=0.35,    absl=185.0/den);              //13 W - anod wires
  AliMedium(kW, "RichW", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  
  if(P()->IsRadioSrc()){
    AliInfo("Special radioactive source materials");
    den=7.87;
    AliMaterial(++material, "RichSteel",  a=55.845,z=26.0,den,    radl=1.76,    absl=131.9/den);        //14 Steel (Fe)
    AliMedium(kSteel, "RichSteel", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
    AliMaterial(++material, "RichPerpex",  a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);                 //15 Perpex
    AliMedium(kPerpex, "RichPerpex", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    
    AliMaterial(++material, "RichSr90",  a=87.62,z=38.0,den=2.54,    radl=4.24,    absl=0);                  //16 Sr90
    AliMedium(kSr90, "RichSr90", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    
    Float_t aMylar[5]={12.01,1,16}; Float_t zMylar[5]={6,1,8};  Float_t wMylar[5]={5,4,5};                  //17 Mylar C5H4O5
    AliMixture(++material,"RichMylar",aMylar, zMylar, den=1.39, -3, wMylar); 
    AliMedium(kMylar, "RichMylar", material, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  }
  
//Optical properties:
#include "Opticals.h"
  gMC->SetCerenkov((*fIdtmed)[kAir]      , kNbins, aPckov, aAbsCH4    , aQeAll, aIdxCH4);       //1 Air
  gMC->SetCerenkov((*fIdtmed)[kRoha]     , kNbins, aPckov, aAbsCH4    , aQeAll, aIdxCH4);       //2 Honeycomb  
  gMC->SetCerenkov((*fIdtmed)[kSiO2]     , kNbins, aPckov, aAbsSiO2   , aQeAll, aIdxSiO2);      //3 Quartz SiO2 
  gMC->SetCerenkov((*fIdtmed)[kC6F14]    , kNbins, aPckov, aAbsC6F14  , aQeAll, aIdxC6F14);     //4 Freon C6F14
  gMC->SetCerenkov((*fIdtmed)[kCH4]      , kNbins, aPckov, aAbsCH4    , aQeAll, aIdxCH4);       //5 Methane CH4 
  gMC->SetCerenkov((*fIdtmed)[kCsI]      , kNbins, aPckov, aAbsCsI    , aQeCsI, aIdxCH4);       //6 CsI
  gMC->SetCerenkov((*fIdtmed)[kGridCu]   , kNbins, aPckov, aAbsGrid   , aQeAll, aIdxMetal);     //7 grid Cu
  gMC->SetCerenkov((*fIdtmed)[kOpSiO2]   , kNbins, aPckov, aAbsOpSiO2 , aQeAll, aIdxMetal);     //8 Opaque quartz SiO2
  gMC->SetCerenkov((*fIdtmed)[kGap]      , kNbins, aPckov, aAbsCH4    , aQeAll, aIdxCH4);       //9 Special methane gap
  gMC->SetCerenkov((*fIdtmed)[kAl]       , kNbins, aPckov, aAbsGrid   , aQeAll, aIdxMetal);     //10 Aluminium
  gMC->SetCerenkov((*fIdtmed)[kGlass]    , kNbins, aPckov, aAbsOpSiO2 , aQeAll, aIdxMetal);     //11 Glass    
  gMC->SetCerenkov((*fIdtmed)[kGel]      , kNbins, aPckov, aAbsGel    , aQeAll, aIdxGel);       //12 Aerogel
  gMC->SetCerenkov((*fIdtmed)[kReflector], kNbins, aPckov, aAbsRef    , aQeAll, aIdxMetal);     //13 Aerogel reflector
}//void AliRICH::CreateMaterials()
//__________________________________________________________________________________________________
Float_t AliRICH::Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)
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
    HitsCreate();      //branch will be created in AliDetector::MakeBranch
  }//H     
  AliDetector::MakeBranch(option);//this is after cH because we need to guarantee that fHits array is created
      
  if(cS&&fLoader->TreeS()){//S  
    SDigitsCreate();   MakeBranchInTree(fLoader->TreeS(),"RICH",&fSdigits,kBufferSize,0) ;
  }//S
   
  if(cD&&fLoader->TreeD()){//D
    DigitsCreate();
    for(Int_t i=0;i<kNchambers;i++){ 
      MakeBranchInTree(fLoader->TreeD(),Form("%s%d",GetName(),i+1),&((*fDigs)[i]),kBufferSize,0);
    }
  }//D
  
  if(cR&&fLoader->TreeR()){//R
    ClustersCreate();
    for(Int_t i=0;i<kNchambers;i++)
      MakeBranchInTree(fLoader->TreeR(),Form("%sClusters%d",GetName(),i+1), &((*fClus)[i]), kBufferSize, 0);    
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
    HitsCreate();//branch map will be in AliDetector::SetTreeAddress    
  }//H
  AliDetector::SetTreeAddress();//this is after TreeH because we need to guarantee that fHits array is created

  if(fLoader->TreeS()){//S
    AliDebug(1,"tree S is requested.");
    branch=fLoader->TreeS()->GetBranch(GetName());        if(branch){SDigitsCreate();   branch->SetAddress(&fSdigits);}
  }//S
    
  if(fLoader->TreeD()){//D    
    AliDebug(1,"tree D is requested.");
    for(int i=0;i<kNchambers;i++){      
      branch=fLoader->TreeD()->GetBranch(Form("%s%d",GetName(),i+1)); 
      if(branch){DigitsCreate(); branch->SetAddress(&((*fDigs)[i]));}
    }
  }//D
    
  if(fLoader->TreeR()){//R
    AliDebug(1,"tree R is requested.");
    for(int i=0;i<kNchambers;i++){         
      branch=fLoader->TreeR()->GetBranch(Form("%sClusters%d" ,GetName(),i+1));
      if(branch){ClustersCreate(); branch->SetAddress(&((*fClus)[i]));}
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
void AliRICH::ControlPlots()
{ 
//Creates a set of QA hists to control the results of simulation. Hists are in file $HOME/RCP.root     
  TH1F             *pElecP=0 ,*pMuonP=0 ,*pPionP=0 ,*pKaonP=0 ,*pProtP=0,  //stack particles
                   *pHxD=0,*pHyD=0,*pHxSd=0,*pHySd=0,      //diff hit position - digit sdigit position 
                   *pNumClusH1=0,
                   *pQdcH1=0,       *pSizeH1=0,
                   *pPureMipQdcH1=0,*pPureMipSizeH1=0,
                   *pPureCerQdcH1=0,*pPureCerSizeH1=0,
                   *pPureFeeQdcH1=0,*pPureFeeSizeH1=0,
                   *pMipQdcH1=0,    *pPhotQdcH1=0;  
  TH2F *pMapH2=0,*pPureMipMapH2=0,*pPureCerMapH2=0,*pPureFeeMapH2=0;
  TH1F *pelecRadius=0,*pprotRadius=0,*pprotbarRadius=0;
//load all information  
                 GetLoader()->GetRunLoader()->LoadHeader();  
                 GetLoader()->GetRunLoader()->LoadKinematics();  
                 GetLoader()->LoadHits();  
  Bool_t isSdig=0;//!GetLoader()->LoadSDigits();
  Bool_t isDig =0;//!GetLoader()->LoadDigits();
  Bool_t isClus=!GetLoader()->LoadRecPoints();

  gBenchmark->Start("ControlPlots");
    
  TFile *pFile = new TFile("$(HOME)/RCP.root","RECREATE");   
  
  pElecP=new TH1F("Pelec","Electrons made hit in RICH;p [GeV]",1000,-30,30); 
  pMuonP=new TH1F("Pmuon","Muons made hit in RICH;p [GeV]",1000,-30,30); 
  pPionP=new TH1F("Ppion","Pions made hit in RICH;p [GeV]",1000,-30,30); 
  pKaonP=new TH1F("Pkaon","Kaon made hit in RICH;p [GeV]",1000,-30,30); 
  pProtP=new TH1F("Pprot","Protons made hit in RICH;p [GeV]",1000,-30,30); 
  pelecRadius=new TH1F("elecRadius","elec",600,0.,600.);  
  pprotRadius=new TH1F("protRadius","elec",600,0.,600.);  
  pprotbarRadius=new TH1F("protbarRadius","elec",600,0.,600.);
    
  if(isSdig){
    AliInfo("SDigits available");
    pHxSd=new TH1F("DiffHitSDigitX","Hit-SDigit diff X all chambers;diff [cm]",300,-10,10); 
    pHySd=new TH1F("DiffHitSDigitY","Hit-SDigit diff Y all chambers;diff [cm]",300,-10,10); 
  }//isSdig
  
  if(isDig){
    AliInfo("Digits available");
    pHxD=new TH1F("DiffHitDigitX","Hit-Digit diff X all chambers;diff [cm]",300,-10,10); 
    pHyD=new TH1F("DiffHitDigitY","Hit-Digit diff Y all chambers;diff [cm]",300,-10,10); 
  }//isDig
  
  if(isClus){ 
    AliInfo("Clusters available");
    pNumClusH1=new TH1F("NumClusPerEvent","Number of clusters per event;number",50,0,49);
    
    pQdcH1        =new TH1F("ClusQdc",   "Cluster Charge all chambers;q [QDC]",P()->MaxQdc(),0,P()->MaxQdc());
    pSizeH1       =new TH1F("ClusSize",  "Cluster size all chambers;size [number of pads in cluster]",100,0,100);
    pMapH2        =new TH2F("ClusMap",   "Cluster map;x [cm];y [cm]",1000,0,P()->PcSizeX(),1000,0,P()->PcSizeY());
  
    pMipQdcH1     =new TH1F("QdcMip"      ,"MIP Cluster Charge all chambers;q [QDC]",P()->MaxQdc(),0,P()->MaxQdc());
    pPhotQdcH1    =new TH1F("QdcPhot"     ,"Cer+Fee Cluster Charge all chambers;q [QDC]",P()->MaxQdc(),0,P()->MaxQdc());
        
    pPureMipQdcH1 =new TH1F("QdcPureMip"  ,"MIP only Cluster Charge all chambers;q [QDC]",P()->MaxQdc(),0,P()->MaxQdc());
    pPureMipSizeH1=new TH1F("SizePureMip" ,"MIP only Cluster size all chambers;size [number of pads in cluster]",100,0,100);
    pPureMipMapH2 =new TH2F("MapPureMip"  ,"MIP only Cluster map;x [cm];y [cm]",1000,0,P()->PcSizeX(),1000,0,P()->PcSizeY());
  
    pPureCerQdcH1 =new TH1F("QdcPureCer"  ,"Cerenkov only Cluster Charge all chambers;q [QDC]",P()->MaxQdc(),0,P()->MaxQdc());
    pPureCerSizeH1=new TH1F("SizePureCer" ,"Cernekov only Cluster size all chambers;size [number of pads in cluster]",100,0,100);
    pPureCerMapH2 =new TH2F("MapPureCer"  ,"Cerenkov only Cluster map;x [cm];y [cm]",1000,0,P()->PcSizeX(),1000,0,P()->PcSizeY());
    
    pPureFeeQdcH1 =new TH1F("QdcPureFee"  ,"Feedback only Cluster Charge all chambers;q [QDC]",P()->MaxQdc(),0,P()->MaxQdc());
    pPureFeeSizeH1=new TH1F("SizePureFee" ,"Feedback only Cluster size all chambers;size [number of pads in cluster]",100,0,100);
    pPureFeeMapH2 =new TH2F("MapPureFee"  ,"Feedback only Cluster map;x [cm];y [cm]",1000,0,P()->PcSizeX(),1000,0,P()->PcSizeY());

  }//isClus
//end of hists booking  
  for(Int_t iEvtN=0;iEvtN < GetLoader()->GetRunLoader()->GetAliRun()->GetEventsPerRun();iEvtN++){//events loop
    GetLoader()->GetRunLoader()->GetEvent(iEvtN);    //get current event
    
    for(Int_t iPrimN=0;iPrimN < GetLoader()->TreeH()->GetEntries();iPrimN++){//hit tree loop
      GetLoader()->TreeH()->GetEntry(iPrimN);      
      for(Int_t j=0;j<Hits()->GetEntries();j++){//hits loop
        AliRICHHit *pHit = (AliRICHHit*)Hits()->At(j);
        TParticle *pParticle = GetLoader()->GetRunLoader()->Stack()->Particle(pHit->GetTrack());//get particle produced this hit
        Double_t dRadius = TMath::Sqrt(pParticle->Vx()*pParticle->Vx()+pParticle->Vy()*pParticle->Vy()+pParticle->Vz()*pParticle->Vz());
        switch(pParticle->GetPdgCode()){
          case kPositron : pElecP->Fill( pParticle->P());pelecRadius->Fill(dRadius); break;
          case kElectron : pElecP->Fill(-pParticle->P());pelecRadius->Fill(dRadius); break;
          
          case kMuonPlus : pMuonP->Fill( pParticle->P()); break;
          case kMuonMinus: pMuonP->Fill(-pParticle->P()); break;
                    
          case kPiPlus   : pPionP->Fill( pParticle->P()); break;
          case kPiMinus  : pPionP->Fill(-pParticle->P()); break;
          
          case kKPlus    : pKaonP->Fill( pParticle->P()); break;
          case kKMinus   : pKaonP->Fill(-pParticle->P()); break;
          
          case kProton   : pProtP->Fill( pParticle->P()); pprotRadius->Fill(dRadius); break;
          case kProtonBar: pProtP->Fill(-pParticle->P()); pprotbarRadius->Fill(dRadius); break;
              
        }//switch PdgCode
            
      }//hits loop
    }//hit tree loop
    
    if(isSdig){
      GetLoader()->TreeS()->GetEntry(0);  
      for(Int_t iSdigN=0;iSdigN<SDigits()->GetEntries();iSdigN++){//sdigits loop 
        AliRICHDigit *pSdig=(AliRICHDigit*)SDigits()->At(iSdigN); //get current sdigit pointer  
        AliRICHHit   *pHit=Hit(pSdig->GetTrack(0));               //get hit of this sdigit (always one)
        TVector2 hit2 =C(pHit->C())->Mrs2Pc(pHit->OutX3());       //this hit position  in local system
        TVector2 sdig2=P()->Pad2Loc(pSdig->Pad());                //center of pad for this sdigit
        pHxSd->Fill(hit2.X()-sdig2.X());        
        pHySd->Fill(hit2.Y()-sdig2.Y());      
      }//sdigits loop
    }//if(isSdig)
        
    if(isDig)  GetLoader()->TreeD()->GetEntry(0);  
    if(isClus) GetLoader()->TreeR()->GetEntry(0);
    
    for(Int_t iChamN=1;iChamN<=7;iChamN++){//chambers loop
      if(isDig){
        for(Int_t iDigN=0;iDigN<Digits(iChamN)->GetEntries();iDigN++){//digits loop
          AliRICHDigit *pDig=(AliRICHDigit*)Digits(iChamN)->At(iDigN);
          AliRICHHit   *pHit=Hit(pDig->GetTrack(0));           //get first hit of this digit
          TVector2 hitV2=C(iChamN)->Mrs2Pc(pHit->OutX3()); 
          TVector2 digV2=P()->Pad2Loc(pDig->Pad());            //center of pad for this digit
          pHxD->Fill(hitV2.X()-digV2.X());        pHyD->Fill(hitV2.Y()-digV2.Y());
        }//digits loop
      }//isDig
      if(isClus){
        Int_t iNclusCham=Clusters(iChamN)->GetEntries(); if(iNclusCham) pNumClusH1->Fill(iNclusCham);//number of clusters per event
        for(Int_t iClusN=0;iClusN<iNclusCham;iClusN++){//clusters loop
          AliRICHCluster *pClus=(AliRICHCluster*)Clusters(iChamN)->At(iClusN);
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
    }//chambers loop
    Info("ControlPlots","Event %i processed.",iEvtN);
  }//events loop 
             GetLoader()->UnloadHits();
  if(isSdig) GetLoader()->UnloadSDigits();
  if(isDig)  GetLoader()->UnloadDigits();
  if(isClus) GetLoader()->UnloadRecPoints();
  
  GetLoader()->GetRunLoader()->UnloadHeader();  
  GetLoader()->GetRunLoader()->UnloadKinematics();  
  
  pFile->Write(); delete pFile;
  
  gBenchmark->Show("ControlPlots");
}//ControlPlots()
//__________________________________________________________________________________________________
AliRICHHit* AliRICH::Hit(Int_t tid)const
{
//defines which hit provided by given tid for the currently loaded event
  GetLoader()->LoadHits();
  for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop      
    GetLoader()->TreeH()->GetEntry(iPrimN);
    for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){
      AliRICHHit *pHit=(AliRICHHit*)Hits()->At(iHitN);
      if(tid==pHit->Track()) {GetLoader()->UnloadHits();return pHit;}
    }//hits
  }//prims loop
  GetLoader()->UnloadHits();
  return 0;
}
//__________________________________________________________________________________________________
void AliRICH::HitsPrint(Int_t iEvtN)const
{
//Prints a list of RICH hits for a given event. Default is event number 0.
  if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
  AliInfo(Form("List of RICH hits for event %i",iEvtN));
  if(GetLoader()->LoadHits()) return;
  
  Int_t iTotalHits=0;
  for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
    GetLoader()->TreeH()->GetEntry(iPrimN);      
    Hits()->Print();
    iTotalHits+=Hits()->GetEntries();
  }
  GetLoader()->UnloadHits();
  AliInfo(Form("totally %i hits",iTotalHits));
}
//__________________________________________________________________________________________________
void AliRICH::SDigitsPrint(Int_t iEvtN)const
{
//prints a list of RICH sdigits  for a given event
  if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
  Info("PrintSDigits","List of RICH sdigits for event %i",iEvtN);
  if(GetLoader()->LoadSDigits()) return;
  
  GetLoader()->TreeS()->GetEntry(0);
  SDigits()->Print();
  GetLoader()->UnloadSDigits();
  Info("PrintSDigits","totally %i sdigits",SDigits()->GetEntries());
}
//__________________________________________________________________________________________________
void AliRICH::DigitsPrint(Int_t iEvtN)const
{
//prints a list of RICH digits  for a given event
  if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
  Info("PrintDigits","List of RICH digits for event %i",iEvtN);
  if(GetLoader()->LoadDigits()) return;
  
  Int_t iTotalDigits=0;
  GetLoader()->TreeD()->GetEntry(0);
  for(Int_t iChamber=1;iChamber<=kNchambers;iChamber++){
    Digits(iChamber)->Print();
    iTotalDigits+=Digits(iChamber)->GetEntries();
  }
  GetLoader()->UnloadDigits();
  Info("PrintDigits","totally %i Digits",iTotalDigits);
}
//__________________________________________________________________________________________________
void AliRICH::OccupancyPrint(Int_t iEvtNreq)const
{
//prints occupancy for each chamber in a given event
  Int_t iEvtNmin,iEvtNmax;
  if(iEvtNreq==-1){
    iEvtNmin=0;
    iEvtNmax=gAlice->GetEventsPerRun();
  } else { 
    iEvtNmin=iEvtNreq;iEvtNmax=iEvtNreq+1;
  }
    
  if(GetLoader()->GetRunLoader()->LoadHeader()) return;    
  if(GetLoader()->GetRunLoader()->LoadKinematics()) return;    
  
//  Info("Occupancy","for event %i",iEvtN);
  if(GetLoader()->LoadHits()) return;
  if(GetLoader()->LoadDigits()) return;

  Int_t totPadsPerChamber = AliRICHParam::NpadsX()*AliRICHParam::NpadsY();  

  Int_t nDigCh[kNchambers]={0,0,0,0,0,0,0};  
  Int_t iChHits[kNchambers]={0,0,0,0,0,0,0};
  Int_t nPrim[kNchambers]={0,0,0,0,0,0,0};
  Int_t nSec[kNchambers]={0,0,0,0,0,0,0};
  
  for(Int_t iEvtN=iEvtNmin;iEvtN<iEvtNmax;iEvtN++){
    if(iEvtN%10==0) AliInfo(Form("events processed %i",iEvtN));
    if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
    AliStack *pStack = GetLoader()->GetRunLoader()->Stack();
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      GetLoader()->TreeH()->GetEntry(iPrimN);      
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){
        AliRICHHit *pHit = (AliRICHHit *)Hits()->At(iHitN);
        if(pHit->Eloss()>0){
          iChHits[pHit->C()-1]++;
          if(pStack->Particle(pHit->GetTrack())->Rho()<0.01) nPrim[pHit->C()-1]++;else nSec[pHit->C()-1]++;
        }
      }
    }
    GetLoader()->TreeD()->GetEntry(0);
    for(Int_t iChamber=1;iChamber<=kNchambers;iChamber++) nDigCh[iChamber-1]= Digits(iChamber)->GetEntries();
  }  
  for(Int_t iChamber=1;iChamber<=kNchambers;iChamber++){
    Double_t occupancy = (Double_t)nDigCh[iChamber-1]/(Double_t)totPadsPerChamber;
    Info("Occupancy","for chamber %i = %4.2f %% and charged prim tracks %i and sec. tracks %i with total %i",
        iChamber,occupancy*100.,nPrim[iChamber-1],nSec[iChamber-1],iChHits[iChamber-1]);
  }    
  GetLoader()->UnloadHits();
  GetLoader()->UnloadDigits();
  GetLoader()->GetRunLoader()->UnloadHeader();    
  GetLoader()->GetRunLoader()->UnloadKinematics();    
}
//__________________________________________________________________________________________________
void AliRICH::ClustersPrint(Int_t iEvtN)const
{
//prints a list of RICH clusters  for a given event
  AliInfo(Form("List of RICH clusters for event %i",iEvtN));
  GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(GetLoader()->LoadRecPoints()) return;
  
  Int_t iTotalClusters=0;
  GetLoader()->TreeR()->GetEntry(0);
  for(Int_t iChamber=1;iChamber<=kNchambers;iChamber++){
    Clusters(iChamber)->Print();
    iTotalClusters+=Clusters(iChamber)->GetEntries();
  }
  GetLoader()->UnloadRecPoints();
  AliInfo(Form("totally %i clusters for event %i",iTotalClusters,iEvtN));
}
//__________________________________________________________________________________________________
void AliRICH::PrintTracks(Int_t iEvtN)
{
//prints a list of tracks (including secondary) for a given event
  AliInfo(Form("List of all tracks for event %i",iEvtN));
  GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(GetLoader()->GetRunLoader()->LoadHeader()) return;
  if(GetLoader()->GetRunLoader()->LoadKinematics()) return;
  AliStack *pStack=GetLoader()->GetRunLoader()->Stack();
  
  for(Int_t i=0;i<pStack->GetNtrack();i++) pStack->Particle(i)->Print();
  
  AliInfo(Form("totally %i tracks including %i primaries for event %i",pStack->GetNtrack(),pStack->GetNprimary(),iEvtN));
  GetLoader()->GetRunLoader()->UnloadHeader();
  GetLoader()->GetRunLoader()->UnloadKinematics();
}
//__________________________________________________________________________________________________
void AliRICH::GeomPadPanelFrame()const
{
//Pad Panel frame  6 sectors
  Double_t cm=1,mm=0.1*cm;//default is cm
  Float_t par[3];
  
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
}//GeomPadPanelFrame()
//__________________________________________________________________________________________________
void AliRICH::GeomAmpGap()const
{
//Gap - anod wires 6 copies to RICH
  Double_t cm=1,mm=0.1*cm,mkm=0.001*mm;//default is cm
  Int_t matrixIdReturn=0; //matrix id returned by AliMatrix
  Float_t par[3];


  par[0]=648*mm/2;par[1]=  411*mm/2 ;par[2]=4.45*mm/2;gMC->Gsvolu("RGAP","BOX ",(*fIdtmed)[kCH4] ,par,3);//xy as PPF 2001P2 z WP 2099P1
  par[0]=  0*mm  ;par[1]=  20*mkm/2 ;par[2]= 648*mm/2;gMC->Gsvolu("RANO","TUBE",(*fIdtmed)[kW]   ,par,3);//WP 2099P1 z = gap x PPF 2001P2
  AliMatrix(matrixIdReturn,180,0, 90,90, 90,0); //wires along x
  
  gMC->Gspos("RGAP",1,"RICH",    -335*mm,      -433*mm,8*cm-2.225*mm, 0,"ONLY"); //F1 2040P1 z WP 2099P1
  gMC->Gspos("RGAP",2,"RICH",    +335*mm,      -433*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("RGAP",3,"RICH",    -335*mm,         0*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("RGAP",4,"RICH",    +335*mm,         0*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("RGAP",5,"RICH",    -335*mm,      +433*mm,8*cm-2.225*mm, 0,"ONLY"); 
  gMC->Gspos("RGAP",6,"RICH",    +335*mm,      +433*mm,8*cm-2.225*mm, 0,"ONLY"); 
  for(int i=1;i<=96;i++)
    gMC->Gspos("RANO",i,"RGAP",     0*mm, -411/2*mm+i*4*mm, 0.185*mm, matrixIdReturn,"ONLY"); //WP 2099P1  
}//GeomAmpGap()
//__________________________________________________________________________________________________
void AliRICH::GeomRadiators()const
{
//Defines radiators geometry  
  Double_t mm=0.1;//default is cm
  Float_t par[3];
  par[0]=1330*mm/2 ;par[1]= 413*mm/2  ;par[2]=  24*mm/2;  gMC->Gsvolu("RRAD","BOX ",(*fIdtmed)[kC6F14]     ,par,3); // Rad 2011P1
  par[0]=1330*mm/2 ;par[1]= 413*mm/2  ;par[2]=   4*mm/2;  gMC->Gsvolu("RRFR","BOX ",(*fIdtmed)[kRoha]      ,par,3); //front 
  par[0]=1330*mm/2 ;par[1]= 413*mm/2  ;par[2]=   5*mm/2;  gMC->Gsvolu("RRWI","BOX ",(*fIdtmed)[kSiO2]      ,par,3); //window
  par[0]=1330*mm/2 ;par[1]=   5*mm/2  ;par[2]=  15*mm/2;  gMC->Gsvolu("RRLO","BOX ",(*fIdtmed)[kRoha]      ,par,3); //long side  
  par[0]=  10*mm/2 ;par[1]= 403*mm/2  ;par[2]=  15*mm/2;  gMC->Gsvolu("RRSH","BOX ",(*fIdtmed)[kRoha]      ,par,3); //short side 
  par[0]=   0      ;par[1]=  10*mm/2  ;par[2]=  15*mm/2;  gMC->Gsvolu("RRSP","TUBE",(*fIdtmed)[kSiO2]      ,par,3); //spacer        
    
  gMC->Gspos("RRAD",1,"RICH",   0*mm,-434*mm,   -12*mm,  0,"ONLY"); //3 radiators to RICH
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
        gMC->Gspos("RRSP",10*i+j,"RRAD",-1330*mm/2+116*mm+j*122*mm,(i-1)*105*mm,-0.5*mm,0,"ONLY");//spacers
}//GeomRadiators()
//__________________________________________________________________________________________________
void AliRICH::GeomSandBox()const
{
//Defines SandBox geometry
  Double_t mm=0.1;//default is cm
  Float_t par[3];
  par[0]=1419*mm/2 ;par[1]=1378*mm/2;par[2]=50.5*mm/2; gMC->Gsvolu("RSNB","BOX ",(*fIdtmed)[kAir]  ,par,3);  //2072P1   
  par[0]=1419*mm/2 ;par[1]=1378*mm/2;par[2]= 0.5*mm/2; gMC->Gsvolu("RSCO","BOX ",(*fIdtmed)[kAl]   ,par,3);  //cover
  par[0]=1359*mm/2 ;par[1]=1318*mm/2;par[2]=49.5*mm/2; gMC->Gsvolu("RSHO","BOX ",(*fIdtmed)[kRoha] ,par,3); //honeycomb structure 
  
  gMC->Gspos("RSNB",1,"RICH",   0*mm, 0*mm, -73.75*mm, 0,"ONLY"); //p.84 TDR sandbox to rich
    gMC->Gspos("RSHO",1,"RSNB", 0*mm, 0*mm,      0*mm, 0,"ONLY"); //2072P1 honeycomv to sandbox
    gMC->Gspos("RSCO",1,"RSNB", 0*mm, 0*mm,    +25*mm, 0,"ONLY"); //cover to sandbox
    gMC->Gspos("RSCO",2,"RSNB", 0*mm, 0*mm,    -25*mm, 0,"ONLY"); //cover to sandbox
}//GeomSandBox()
//__________________________________________________________________________________________________
void AliRICH::GeomRadioSrc()const
{
// Defines geometry for radioactive source  
  Double_t cm=1,mm=0.1*cm,mkm=0.001*cm;
  Float_t par[3];
  
  par[0]=0 ;par[1]= 70*mm/2  ;par[2]=  30*mm/2;      gMC->Gsvolu("RSRC","TUBE",(*fIdtmed)[kCH4]    ,par,3); //top src container
    par[0]=0 ;par[1]= 38*mm/2  ;par[2]=  21.8*mm/2;  gMC->Gsvolu("RSAG","TUBE",(*fIdtmed)[kAl]     ,par,3); //Al glass
      par[0]=0 ;par[1]= 34*mm/2  ;par[2]=  20*mm/2;  gMC->Gsvolu("RSPP","TUBE",(*fIdtmed)[kPerpex] ,par,3); //perpex plug
        par[0]=0 ;par[1]= 5*mm/2  ;par[2]=  15*mm/2; gMC->Gsvolu("RSSC","TUBE",(*fIdtmed)[kSteel]  ,par,3); //steel screw in center of perpex
        par[0]=0 ;par[1]= 2*mm/2  ;par[2]=  10*mm/2; gMC->Gsvolu("RSSS","TUBE",(*fIdtmed)[kSteel]  ,par,3); //Steel screw to support Sr90 
          par[0]=0 ;par[1]= 1*mm/2  ;par[2]= 1*mm/2; gMC->Gsvolu("RSSR","TUBE",(*fIdtmed)[kSr90]   ,par,3); //Sr90 source
        par[0]=0 ;par[1]= 4*mm/2  ;par[2]= 10*mm/2;  gMC->Gsvolu("RSWP","TUBE",(*fIdtmed)[kAir]    ,par,3); //Air hole in perpex plug      
      par[0]=0 ;par[1]= 5*mm/2  ;par[2]= 1.8*mm/2;   gMC->Gsvolu("RSWA","TUBE",(*fIdtmed)[kAir]    ,par,3); //Air hole in Al glass bottom
    par[0]=0 ;par[1]= 30*mm/2  ;par[2]= 50*mkm/2;    gMC->Gsvolu("RSMF","TUBE",(*fIdtmed)[kMylar]  ,par,3); //Mylar foil                
    
  gMC->Gspos("RSRC",1,"RICH",       30*cm,        0,     1*cm, 0,"ONLY"); //source to RICH
    gMC->Gspos("RSMF",1,"RSRC",         0,        0,21.8*mm/2+50*mkm/2, 0,"ONLY");//mylar foil to top src volume
    gMC->Gspos("RSAG",1,"RSRC",         0,        0,        0, 0,"ONLY");//Al glass to fake Src volume 
      gMC->Gspos("RSWA",1,"RSAG",    6*mm,        0,   -10*mm, 0,"ONLY");//air whole in al glass bottom
      gMC->Gspos("RSPP",1,"RSAG",       0,        0,   0.9*mm, 0,"ONLY");//perpex plug to Al glass
        gMC->Gspos("RSWP",1,"RSPP",  6*mm,        0,    -5*mm, 0,"ONLY");//air whole in perpex plug
        gMC->Gspos("RSSC",1,"RSPP",     0,        0,   2.5*mm, 0,"ONLY");//steel screw in center of perpex plug
        gMC->Gspos("RSSS",1,"RSPP",  6*mm,        0,     5*mm, 0,"ONLY");//steel screw to support Sr90  in perpex plug
          gMC->Gspos("RSSR",1,"RSSS",   0,        0,  -4.5*mm, 0,"ONLY");//Sr90  in support steel screw
}//GeomSr90()
//__________________________________________________________________________________________________
void AliRICH::GeomAerogel()const
{
//Creates detailed geometry for aerogel study.
  AliDebug(1,"Start.");
  Double_t cm=1;
  Float_t par[3]; //tmp array for volume dimentions
       
  par[0]=10.1*cm/2;par[1]=10.1*cm/2;par[2]=10.1*cm/2;
  gMC->Gsvolu("RREF","BOX ",(*fIdtmed)[kReflector],par,3);//reflector box
  gMC->Gspos("RREF",1,"RICH",0,0,0,0, "ONLY");            //put it to RICH volume
  
  par[0]=10*cm/2;par[1]=10*cm/2;par[2]=10*cm/2;
  gMC->Gsvolu("RGEL","BOX ",(*fIdtmed)[kGel],par,3);//10x10x10 cm^3 cubic of aerogel
  gMC->Gspos("RGEL",1,"RREF",0,0,0,0,"ONLY");//put gel cell to reflector
  AliDebug(1,"Stop.");  
}//GeomAerogel()
//__________________________________________________________________________________________________
void AliRICH::CreateGeometry()
{
//Creates detailed geometry simulation (currently GEANT volumes tree)         
  AliDebug(1,"Start main.");
  Double_t mm=0.1;//default is cm
  Float_t par[3];
  Int_t matrixIdReturn=0; //matrix id returned by AliMatrix
       
//place chambers into mother volume ALIC
  par[0]=(6*mm+1681*mm+6*mm)/2;par[1]=(6*mm+1466*mm+6*mm)/2;par[2]=(80*mm+40*mm)*2/2;
  gMC->Gsvolu("RICH","BOX ",(*fIdtmed)[kCH4],par,3);//2033P1  z p84 TDR
  for(int i=1;i<=P()->Nchambers();i++){ //test configuration with single chamber is taken into account automaticaly in AliRICHParam
    AliMatrix(matrixIdReturn,
                   C(i)->ThetaXd(),C(i)->PhiXd(),  
                   C(i)->ThetaYd(),C(i)->PhiYd(),  
                   C(i)->ThetaZd(),C(i)->PhiZd());
    gMC->Gspos("RICH",i,"ALIC",C(i)->Center().X(),
                               C(i)->Center().Y(),
                               C(i)->Center().Z(),matrixIdReturn, "ONLY");
  }
  
  if(P()->IsAerogel()) 
    GeomAerogel();
  else{
    GeomPadPanelFrame();
    GeomAmpGap();
    if(P()->IsRadioSrc())    GeomRadioSrc(); else GeomRadiators(); 
    GeomSandBox();           
  }
  AliDebug(1,"Stop main.");  
}//CreateGeometry()
//__________________________________________________________________________________________________
void AliRICH::CheckPR()const
{
//Pattern recognition with stack particles
  TFile *pFile = new TFile("$(HOME)/RPR.root","RECREATE","RICH Pattern Recognition");
  TNtupleD *hn = new TNtupleD("hn","ntuple","Pmod:Charge:TrackTheta:TrackPhi:TrackX:TrackY:MinX:MinY:ChargeMIP:ThetaCerenkov:NPhotons:MipIndex:Chamber:Particle");
//  printf("\n\n");
//  printf("Pattern Recognition done for event %5i",0);
  AliMagF * magf = gAlice->Field();
  AliTracker::SetFieldMap(magf);
  for(Int_t iEvtN=0;iEvtN<GetLoader()->GetRunLoader()->GetNumberOfEvents();iEvtN++) {
    GetLoader()->GetRunLoader()->GetEvent(iEvtN);
    AliRICHTracker *tr = new AliRICHTracker();
    tr->RecWithStack(hn);
    Info("CheckPR","Pattern Recognition done for event %i \b",iEvtN);
//    printf("\b\b\b\b\b%5i",iEvtN+1);
  }
  printf("\n\n");
  pFile->Write();pFile->Close();
}
//__________________________________________________________________________________________________
void AliRICH::DisplayEvent(Int_t iEvtNmin,Int_t iEvtNmax)const
{
  TH2F *pDigitsH2[8];

  Bool_t isDigits  =!GetLoader()->LoadDigits();
  if(!isDigits){Error("ShoEvent","No digits. Nothing to display.");return;}
  
  TCanvas *canvas = new TCanvas("RICHDisplay","RICH Display",0,0,1226,900);     
  gStyle->SetPalette(1);

  
  for(Int_t iChamber=1;iChamber<=7;iChamber++) {
    pDigitsH2[iChamber] = new TH2F(Form("pDigitsH2_%i",iChamber),Form("Chamber %i",iChamber),165,0,P()->PcSizeX(),144,0,P()->PcSizeY());
    pDigitsH2[iChamber]->SetMarkerColor(kGreen); 
    pDigitsH2[iChamber]->SetMarkerStyle(29); 
    pDigitsH2[iChamber]->SetMarkerSize(0.4);
    pDigitsH2[iChamber]->SetStats(kFALSE);
    pDigitsH2[iChamber]->SetMaximum(300);
  }
  
  if(iEvtNmax>gAlice->GetEventsPerRun()||iEvtNmax==0) iEvtNmax=gAlice->GetEventsPerRun()-1;

  TLatex t;  t.SetTextSize(0.1);
  for(Int_t iEventN=iEvtNmin;iEventN<=iEvtNmax;iEventN++) {//events loop
    canvas->Divide(3,3);    
    canvas->cd(1);
    t.DrawText(0.2,0.4,Form("Event Number %i",iEventN));        

    GetLoader()->GetRunLoader()->GetEvent(iEventN); //get event
    GetLoader()->TreeD()->GetEntry(0);              //get list of digits 
    for(Int_t iChamber=1;iChamber<=7;iChamber++) {//chambers loop
      pDigitsH2[iChamber]->Reset();    
      for(Int_t j=0;j<Digits(iChamber)->GetEntries();j++) {//digits loop
        AliRICHDigit *pDig = (AliRICHDigit*)Digits(iChamber)->At(j);
        TVector2 x2=AliRICHParam::Pad2Loc(pDig->Pad());
        pDigitsH2[iChamber]->Fill(x2.X(),x2.Y(),pDig->Q());
      }//digits loop
      if(iChamber==1) canvas->cd(7);
      if(iChamber==2) canvas->cd(8);
      if(iChamber==3) canvas->cd(4);
      if(iChamber==4) canvas->cd(5);
      if(iChamber==5) canvas->cd(6);
      if(iChamber==6) canvas->cd(2);
      if(iChamber==7) canvas->cd(3);
      pDigitsH2[iChamber]->Draw("col");
      ReadESD(iEventN,iChamber);
      AliRICHParam::DrawSectors();
    }//chambers loop
    canvas->Update();
    canvas->Modified();

    if(iEventN<iEvtNmax) {gPad->WaitPrimitive();canvas->Clear();}
  }//events loop
}//ShowEvent()
//__________________________________________________________________________________________________
void AliRICH::Display()const
{
//Provides fast event display
//For RICH only, full display is .x Display.C    
  Bool_t isHits    =!GetLoader()->LoadHits();
  Bool_t isDigits  =!GetLoader()->LoadDigits();
  Bool_t isClusters=!GetLoader()->LoadRecPoints();
  
  if(!isHits && !isDigits && !isClusters){Error("Exec","No hits digits and clusters. Nothing to display.");return;}
  
  TCanvas *pCanvas = new TCanvas("Display","RICH Display",0,0,600,600);
  
  TH2F *pHitsH2=0,*pDigitsH2=0,*pClustersH2=0;
  
  if(isHits)     pHitsH2     = new TH2F("pHitsH2"  ,  "Event Display;x,cm;y,cm",165,0,AliRICHParam::PcSizeX(),
                                                                                144,0,AliRICHParam::PcSizeY());
  if(pHitsH2)    pHitsH2->SetStats(kFALSE);
  
  if(isDigits)   pDigitsH2   = new TH2F("pDigitsH2"  ,"Event Display",165,0,AliRICHParam::PcSizeX(),
                                                                      144,0,AliRICHParam::PcSizeY());
  if(isClusters) pClustersH2 = new TH2F("pClustersH2","Event Display",165,0,AliRICHParam::PcSizeX(),
                                                                      144,0,AliRICHParam::PcSizeY());
  
  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events Loop
    GetLoader()->GetRunLoader()->GetEvent(iEventN);  
//display all the staff on chamber by chamber basis           
    for(Int_t iChamber=1;iChamber<=7;iChamber++){//chambers loop       
      if(isHits)     pHitsH2    ->Reset();     
      if(isDigits)   pDigitsH2  ->Reset();     
      if(isClusters) pClustersH2->Reset();
//deals with hits
      for(Int_t i=0;i<GetLoader()->TreeH()->GetEntries();i++){//TreeH loop
        GetLoader()->TreeH()->GetEntry(i);
        for(Int_t j=0;j<Hits()->GetEntries();j++){//hits loop
          AliRICHHit *pHit = (AliRICHHit*)Hits()->At(j);
          if(pHit->C()==iChamber){
            TVector3 hitGlobX3= pHit->OutX3();
            TVector2 hitLocX2 = C(iChamber)->Mrs2Pc(hitGlobX3);
            pHitsH2->Fill(hitLocX2.X(),hitLocX2.Y(),200);
          }//if
        }//hits loop         
      }//TreeH loop
      pHitsH2->SetTitle(Form("event %i chamber %2i",iEventN,iChamber));
      pHitsH2->SetMarkerColor(kRed); pHitsH2->SetMarkerStyle(29); pHitsH2->SetMarkerSize(0.4);
      ReadESD(iEventN,iChamber);
      pHitsH2->Draw();
      ReadESD(iEventN,iChamber);
      AliRICHParam::DrawSectors();
      TLatex l; l.SetNDC(); l.SetTextSize(0.02);
      if(!isHits)     {l.SetTextColor(kRed)  ;l.DrawLatex(0.1,0.01,"No Hits"    );}
      if(!isDigits)   {l.SetTextColor(kGreen);l.DrawLatex(0.4,0.01,"No DIGITS"  );}
      if(!isClusters) {l.SetTextColor(kBlue) ;l.DrawLatex(0.8,0.01,"No CLUSTERS");}
      pCanvas->Update();        pCanvas->Modified();       gPad->WaitPrimitive();
//deals with digits      
      if(isDigits){
        GetLoader()->TreeD()->GetEntry(0);
        for(Int_t j=0;j<Digits(iChamber)->GetEntries();j++){//digits loop
          AliRICHDigit *pDig = (AliRICHDigit*)Digits(iChamber)->At(j);
	  TVector2 x2=AliRICHParam::Pad2Loc(pDig->Pad());
	  pDigitsH2->Fill(x2.X(),x2.Y(),100);
        }//digits loop
        pDigitsH2->SetMarkerColor(kGreen); pDigitsH2->SetMarkerStyle(29); pDigitsH2->SetMarkerSize(0.4);
        pDigitsH2->Draw("same");
        pCanvas->Update();        pCanvas->Modified();       gPad->WaitPrimitive();
      }//if(isDigits)      
//deals with clusters      
      if(isClusters){
        GetLoader()->TreeR()->GetEntry(0);
        for(Int_t j=0;j<Clusters(iChamber)->GetEntries();j++){//clusters loop
          AliRICHCluster *pClus = (AliRICHCluster*)Clusters(iChamber)->At(j);
          pClustersH2->Fill(pClus->X(),pClus->Y(),50);
        }//clusters loop
        pClustersH2->SetMarkerColor(kBlue); pClustersH2->SetMarkerStyle(29);  pClustersH2->SetMarkerSize(0.4);
        pClustersH2->Draw("same");
        pCanvas->Update();        pCanvas->Modified();       gPad->WaitPrimitive();
      }//if(isClusters)
    }//chambers loop
  }//events Loop
  
  delete pCanvas;
  GetLoader()->UnloadHits();
  if(isDigits)   GetLoader()->UnloadDigits();
  if(isClusters) GetLoader()->UnloadRecPoints();
}//Display()
//__________________________________________________________________________________________________
Int_t AliRICH::Nparticles(Int_t iPartID,Int_t iEvtN,AliRunLoader *pRL)
{
//counts total number of particles of given type (including secondary) for a given event
  pRL->GetEvent(iEvtN);    
  if(pRL->LoadHeader()) return 0;
  if(pRL->LoadKinematics()) return 0;
  AliStack *pStack=pRL->Stack();
  
  Int_t iCounter=0;
  for(Int_t i=0;i<pStack->GetNtrack();i++){
    if(pStack->Particle(i)->GetPdgCode()==iPartID) iCounter++;
  }
  
  pRL->UnloadHeader();
  pRL->UnloadKinematics();
  return iCounter;
}
//__________________________________________________________________________________________________
void AliRICH::ReadESD(Int_t iEventN, Int_t iChamber)const
{
//
  AliInfo("Start.");
  TFile *pFile=TFile::Open("AliESDs.root","read");
  if(!pFile || !pFile->IsOpen()) {AliInfo("ESD file not open.");return;}      //open AliESDs.root                                                                    
  TTree *pTree = (TTree*) pFile->Get("esdTree");
  if(!pTree){AliInfo("ESD not found.");return;}                               //get ESD tree
  
  AliInfo("ESD found. Try to draw ring");
                                                                 
  AliESD *pESD=new AliESD;  pTree->SetBranchAddress("ESD", &pESD);
  
  pTree->GetEvent(iEventN);
  
  Double_t b = pESD->GetMagneticField()/10.;
  
  Int_t iNtracks=pESD->GetNumberOfTracks();    
  
  for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){//ESD tracks loop
    AliESDtrack *pTrack = pESD->GetTrack(iTrackN);// get next reconstructed track
    Int_t charge = (Int_t)(-TMath::Sign(1.,pTrack->GetSign()*b));
    AliRICHHelix helix(pTrack->X3(),pTrack->P3(),charge,b);
    Int_t iChamberOnRICH=helix.RichIntersect(P());        
    if(iChamberOnRICH==iChamber) {
//
      TMarker *trackImpact = new TMarker(helix.PosPc().X(),helix.PosPc().Y(),kStar);
      trackImpact->SetMarkerColor(kRed);
      trackImpact->Draw();
//
      Int_t iChamberRecon = pTrack->GetRICHcluster()/100000;
      if(iChamberRecon==iChamber) {
        Double_t thetaCer = pTrack->GetRICHsignal();
        if(thetaCer<0) continue;
        TVector3 entrance(helix.PosRad().X(),helix.PosRad().Y(),0);
        Double_t thetaTrack,phiTrack;
        pTrack->GetRICHthetaPhi(thetaTrack,phiTrack);
        TVector3 vectorTrack;
        vectorTrack.SetMagThetaPhi(pTrack->GetP(),thetaTrack,phiTrack);
        AliInfo(Form("Draw ring started for track %i on chamber %i",iTrackN,iChamber));
        AliInfo(Form("ThetaCer %f TrackTheta %f TrackPhi %f Momentum %f",thetaCer,thetaTrack,phiTrack,pTrack->GetP()));
        Double_t dx,dy;
        pTrack->GetRICHdxdy(dx,dy);
        AliInfo(Form("dx %f dy %f ",dx,dy));
        DrawRing(entrance,vectorTrack,thetaCer);
      }
    }
  }
  delete pESD;  pFile->Close();//close AliESDs.root
  AliInfo("Stop.");
}
//__________________________________________________________________________________________________
void AliRICH::DrawRing(TVector3 entrance,TVector3 vectorTrack,Double_t thetaCer)const
{
  Double_t xGraph[100],yGraph[100];
  Int_t nPointsToDraw = 0;
  for(Int_t i=0;i<100;i++) {
    Double_t phiCer = 2*TMath::Pi()*i/100;
    TVector3 pos = AliRICHParam::ForwardTracing(entrance,vectorTrack,thetaCer,phiCer);
    if(pos.X()==-999) continue;
    xGraph[nPointsToDraw] = pos.X();yGraph[nPointsToDraw] = pos.Y();nPointsToDraw++;
  }
//  AliInfo(Form("Npoints per ring %i",nPointsToDraw));
  TGraph *gra = new TGraph(nPointsToDraw,xGraph,yGraph);
  gra->Draw("C");  
}
//__________________________________________________________________________________________________
void AliRICH::SummaryOfEvent(Int_t iEvtN) const
{
//prints a summary for a given event
  AliInfo(Form("Summary of event %i",iEvtN));
  GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(GetLoader()->GetRunLoader()->LoadHeader()) return;
  if(GetLoader()->GetRunLoader()->LoadKinematics()) return;
  AliStack *pStack=GetLoader()->GetRunLoader()->Stack();
  
  AliGenEventHeader* pGenHeader =  gAlice->GetHeader()->GenEventHeader();
  if(pGenHeader->InheritsFrom("AliGenHijingEventHeader")) {
    AliInfo(Form(" Hijing event with impact parameter b = %.2f (fm)",((AliGenHijingEventHeader*) pGenHeader)->ImpactParameter()));
  }
  Int_t nChargedPrimaries=0;
  for(Int_t i=0;i<pStack->GetNtrack();i++) {
    TParticle *pParticle = pStack->Particle(i);
    if(pParticle->IsPrimary()&&pParticle->GetPDG()->Charge()!=0) nChargedPrimaries++;
    }
  AliInfo(Form("Total number of         primaries %i",pStack->GetNprimary()));
  AliInfo(Form("Total number of charged primaries %i",nChargedPrimaries));
  AliInfo(Form("Total n. of tracks in stack(+sec) %i",pStack->GetNtrack()));
  GetLoader()->GetRunLoader()->UnloadHeader();
  GetLoader()->GetRunLoader()->UnloadKinematics();
}
//__________________________________________________________________________________________________
void AliRICH::HitsQA(Double_t cut,Double_t cutele,Double_t cutR)
{
// Provides a set of control plots intended primarily for charged particle flux analisys
// Arguments: cut (GeV)    - cut on momentum of any charged particles but electrons, 
//            cetele (GeV) - the same for electrons-positrons
//            cutR (cm)    - cut on production vertex radius (cylindrical system)        
  gBenchmark->Start("HitsAna");
  
  Double_t cutPantiproton    =cut;
  Double_t cutPkaonminus     =cut;
  Double_t cutPpionminus     =cut;
  Double_t cutPmuonminus     =cut;
  Double_t cutPpositron      =cutele;
                    
  Double_t cutPelectron      =cutele;
  Double_t cutPmuonplus      =cut;
  Double_t cutPpionplus      =cut;
  Double_t cutPkaonplus      =cut;
  Double_t cutPproton        =cut;
                       

  TH2F *pEleHitRZ    =new TH2F("EleHitRZ"    ,Form("e^{+} e^{-} hit %s;z[cm];R[cm]" ,GetName())     , 400,-300,300 ,400,-500,500);   //R-z plot 0cm<R<550cm -300cm<z<300cm  
  TH2F *pEleHitRP    =new TH2F("EleHitRP"    ,Form("e^{+} e^{-} hit %s;p[GeV];R[cm]",GetName())     ,1000,-1  ,1   ,400,   0,550);   //R-p plot 0cm<R<550cm -1GeV<p<1GeV 
  TH1F *pEleAllP     =new TH1F("EleAllP"     ,     "e^{+} e^{-} all;p[GeV]"                         ,1000,-1  ,1                );  
  TH1F *pEleHitP     =new TH1F("EleHitP"     ,Form("e^{+} e^{-} hit %s;p[GeV]"      ,GetName())     ,1000,-1  ,1                );   
  TH1F *pMuoHitP     =new TH1F("MuoHitP"     ,Form("#mu^{-} #mu^{+} hit %s;p[GeV]"  ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pPioHitP     =new TH1F("PioHitP"     ,Form("#pi^{-} #pi^{+} hit %s;p[GeV]"  ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pKaoHitP     =new TH1F("KaoHitP"     ,Form("K^{-} K^{+} hit %s;p[GeV]"      ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pProHitP     =new TH1F("ProHitP"     ,Form("p^{-} p^{+} hit %s;p[GeV]"      ,GetName())     ,1000,-4  ,4                ); 
  TH2F *pFlux        =new TH2F("flux"        ,Form("%s flux with Rvertex<%.1fcm"    ,GetName(),cutR),10  ,-5  ,5   , 10,0    ,10); //special text hist
  TH2F *pVertex      =new TH2F("vertex"      ,Form("%s 2D vertex of RICH hit;x;y"   ,GetName())     ,120 ,0   ,600 ,120,0    ,600); //special text hist
  TH1F *pRho         =new TH1F("rho"         ,Form("%s r of RICH hit"               ,GetName())     ,600 ,0   ,600); //special text hist
  pFlux->SetStats(0);
  pFlux->GetXaxis()->SetBinLabel(1 ,Form("p^{-}>%.3fGeV/c"   ,cutPantiproton));        
  pFlux->GetXaxis()->SetBinLabel(2 ,Form("K^{-}>%.3fGeV/c"   ,cutPkaonminus ));        
  pFlux->GetXaxis()->SetBinLabel(3 ,Form("#pi^{-}>%.3fGeV/c" ,cutPpionminus ));      
  pFlux->GetXaxis()->SetBinLabel(4 ,Form("#mu^{-}>%.3fGeV/c" ,cutPmuonminus ));      
  pFlux->GetXaxis()->SetBinLabel(5 ,Form("e^{+}>%.3fGeV/c"   ,cutPpositron  ));        
  
  pFlux->GetXaxis()->SetBinLabel(6 ,Form("e^{-}>%.3fGeV/c"   ,cutPelectron  ));        
  pFlux->GetXaxis()->SetBinLabel(7 ,Form("#mu^{+}>%.3fGeV/c" ,cutPmuonplus  ));      
  pFlux->GetXaxis()->SetBinLabel(8 ,Form("#pi^{+}>%.3fGeV/c" ,cutPpionplus  ));      
  pFlux->GetXaxis()->SetBinLabel(9 ,Form("K^{+}>%.3fGeV/c"   ,cutPkaonplus  ));        
  pFlux->GetXaxis()->SetBinLabel(10,Form("p^{+}>%.3fGeV/c"   ,cutPproton    ));        
  
  pFlux->GetYaxis()->SetBinLabel(1,"sum");  
  pFlux->GetYaxis()->SetBinLabel(2,"ch1");  
  pFlux->GetYaxis()->SetBinLabel(3,"ch2");  
  pFlux->GetYaxis()->SetBinLabel(4,"ch3");  
  pFlux->GetYaxis()->SetBinLabel(5,"ch4");  
  pFlux->GetYaxis()->SetBinLabel(6,"ch5");  
  pFlux->GetYaxis()->SetBinLabel(7,"ch6");  
  pFlux->GetYaxis()->SetBinLabel(8,"ch7");  
  pFlux->GetYaxis()->SetBinLabel(9,"prim"); 
  pFlux->GetYaxis()->SetBinLabel(10,"tot");  
  
//end of hists definition
   
  Int_t iNevents=fLoader->GetRunLoader()->GetAliRun()->GetEventsPerRun(),iCntPrimParts=0,iCntTotParts=0;
//load all needed trees   
  fLoader->LoadHits(); 
  fLoader->GetRunLoader()->LoadHeader(); 
  fLoader->GetRunLoader()->LoadKinematics();  
  
  for(Int_t iEvtN=0;iEvtN < iNevents;iEvtN++){//events loop
    fLoader->GetRunLoader()->GetEvent(iEvtN);
    AliInfo(Form(" %i event processes",fLoader->GetRunLoader()->GetEventNumber()));
    AliStack *pStack= fLoader->GetRunLoader()->Stack(); 
    
    for(Int_t iParticleN=0;iParticleN<pStack->GetNtrack();iParticleN++){//stack loop
      TParticle *pPart=pStack->Particle(iParticleN);

      if(iParticleN%10000==0) AliInfo(Form(" %i particles read",iParticleN));
    
      switch(pPart->GetPdgCode()){
        case kProtonBar: pFlux->Fill(-4.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-4.5,8); break;
        case kKMinus:    pFlux->Fill(-3.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-3.5,8); break;
        case kPiMinus:   pFlux->Fill(-2.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-2.5,8); break;
        case kMuonMinus: pFlux->Fill(-1.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-1.5,8); break;
        case kPositron:  pFlux->Fill(-0.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-0.5,8); pEleAllP->Fill(-pPart->P()); break;
      
        case kElectron:  pFlux->Fill( 0.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 0.5,8); pEleAllP->Fill( pPart->P()); break;      
        case kMuonPlus:  pFlux->Fill( 1.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 1.5,8); break;      
        case kPiPlus:    pFlux->Fill( 2.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 2.5,8); break;      
        case kKPlus:     pFlux->Fill( 3.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 3.5,8); break;      
        case kProton:    pFlux->Fill( 4.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 4.5,8); break;            
      }//switch
    }//stack loop
//now hits analiser        
    for(Int_t iEntryN=0;iEntryN < fLoader->TreeH()->GetEntries();iEntryN++){//TreeH loop
      fLoader->TreeH()->GetEntry(iEntryN);                                  //get current entry (prim)                
      for(Int_t iHitN=0;iHitN < Hits()->GetEntries();iHitN++){//hits loop
        AliRICHHit *pHit = (AliRICHHit*)Hits()->At(iHitN);            //get current hit
        TParticle  *pPart=pStack->Particle(pHit->GetTrack());      //get stack particle which produced the current hit
        
        if(pPart->GetPDG()->Charge()!=0&&pPart->Rho()>0.1) pVertex->Fill(pPart->Vx(),pPart->Vy()); //safe margin for sec.
        if(pPart->GetPDG()->Charge()!=0) pRho->Fill(pPart->Rho()); //safe margin for sec.
        if(pPart->R()>cutR) continue;                                   //cut on production radius (cylindrical system) 
      
        switch(pPart->GetPdgCode()){
          case kProtonBar: if(pPart->P()>cutPantiproton) {pProHitP->Fill(-pPart->P()); pFlux->Fill(-4.5,pHit->C());}break;
          case kKMinus   : if(pPart->P()>cutPkaonminus)  {pKaoHitP->Fill(-pPart->P()); pFlux->Fill(-3.5,pHit->C());}break;
          case kPiMinus  : if(pPart->P()>cutPpionminus)  {pPioHitP->Fill(-pPart->P()); pFlux->Fill(-2.5,pHit->C());}break;
          case kMuonMinus: if(pPart->P()>cutPmuonminus)  {pMuoHitP->Fill(-pPart->P()); pFlux->Fill(-1.5,pHit->C());}break;        
          case kPositron : if(pPart->P()>cutPpositron)   {pEleHitP->Fill(-pPart->P()); pFlux->Fill(-0.5,pHit->C()); 
               pEleHitRP->Fill(-pPart->P(),pPart->R());  pEleHitRZ->Fill(pPart->Vz(),pPart->R()); }break;
          
          case kElectron : if(pPart->P()>cutPelectron)   {pEleHitP->Fill( pPart->P()); pFlux->Fill( 0.5,pHit->C()); 
               pEleHitRP->Fill( pPart->P(),pPart->R());  pEleHitRZ->Fill(pPart->Vz(),pPart->R()); }break;
          case kMuonPlus : if(pPart->P()>cutPmuonplus)   {pMuoHitP->Fill( pPart->P()); pFlux->Fill( 1.5,pHit->C());}break;                     
          case kPiPlus   : if(pPart->P()>cutPpionplus)   {pPioHitP->Fill( pPart->P()); pFlux->Fill( 2.5,pHit->C());}break;           
          case kKPlus    : if(pPart->P()>cutPkaonplus)   {pKaoHitP->Fill( pPart->P()); pFlux->Fill( 3.5,pHit->C());}break;           
          case kProton   : if(pPart->P()>cutPproton)     {pProHitP->Fill( pPart->P()); pFlux->Fill( 4.5,pHit->C());}break;
        }
      }//hits loop      
    }//TreeH loop
    iCntPrimParts +=pStack->GetNprimary();
    iCntTotParts  +=pStack->GetNtrack();
  }//events loop                        
//unload all loaded staff  
  fLoader->UnloadHits();  
  fLoader->GetRunLoader()->UnloadHeader(); 
  fLoader->GetRunLoader()->UnloadKinematics();  
//Calculater some sums
  Stat_t sum=0;
//sum row, sum over rows  
  for(Int_t i=1;i<=pFlux->GetNbinsX();i++){
    sum=0; for(Int_t j=2;j<=8;j++)    sum+=pFlux->GetBinContent(i,j);    
    pFlux->SetBinContent(i,1,sum);
  }
    
//display everything  
  new TCanvas("canvas1",Form("Events %i Nprims=%i Nparticles=%i",iNevents,iCntPrimParts,iCntTotParts),1000,900); pFlux->Draw("text");  gPad->SetGrid();  
//total prims and particles
  TLatex latex; latex.SetTextSize(0.02);
  sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,10);    latex.DrawLatex(5.1,9.5,Form("%.0f",sum));
  sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,9);     latex.DrawLatex(5.1,8.5,Form("%.0f",sum));
  for(Int_t iChN=1;iChN<=kNchambers;iChN++) {
    sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,iChN+1);latex.DrawLatex(5.1,iChN+0.5,Form("%.0f",sum));
  }  
  sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,1);    latex.DrawLatex(5.1,0.5,Form("%.0f",sum));
  
  new TCanvas("cEleAllP"   ,"e" ,200,100); pEleAllP->Draw();
  new TCanvas("cEleHitRP"  ,"e" ,200,100); pEleHitRP->Draw();
  new TCanvas("cEleHitRZ"  ,"e" ,200,100); pEleHitRZ->Draw();
  new TCanvas("cEleHitP"   ,"e" ,200,100); pEleHitP->Draw();
  new TCanvas("cMuoHitP"   ,"mu",200,100); pMuoHitP->Draw();
  new TCanvas("cPioHitP"   ,"pi",200,100); pPioHitP->Draw();
  new TCanvas("cKaoHitP"   ,"K" ,200,100); pKaoHitP->Draw();
  new TCanvas("cProHitP"   ,"p" ,200,100); pProHitP->Draw();
  new TCanvas("cVertex"    ,"2d vertex" ,200,100); pVertex->Draw();
  new TCanvas("cRho"    ,"Rho of sec" ,200,100); pRho->Draw();
  
  gBenchmark->Show("HitsPlots");
}//HitsPlots()
