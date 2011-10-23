/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber                                                  //
//  This class contains the basic functions for the Time Projection Chamber  //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTPCClass.gif">
*/
//End_Html
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////

//

#include <Riostream.h>
#include <stdlib.h>

#include <TF2.h>
#include <TFile.h>  
#include <TGeoGlobalMagField.h>
#include <TInterpreter.h>
#include <TMath.h>
#include <TMatrixF.h>
#include <TObjectTable.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>     
#include <TTree.h>
#include <TVector.h>
#include <TVirtualMC.h>
#include <TParameter.h>

#include "AliDigits.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliSimDigits.h"
#include "AliTPC.h"
#include "AliTPC.h"
#include "AliTPCDigitsArray.h"
#include "AliTPCLoader.h"
#include "AliTPCPRF2D.h"
#include "AliTPCParamSR.h"
#include "AliTPCRF1D.h"
#include "AliTPCTrackHitsV2.h"
#include "AliTrackReference.h"
#include "AliMC.h"
#include "AliStack.h"
#include "AliTPCDigitizer.h"
#include "AliTPCBuffer.h"
#include "AliTPCDDLRawData.h"
#include "AliLog.h"
#include "AliTPCcalibDB.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "AliTPCExB.h"
#include "AliRawReader.h"
#include "AliTPCRawStreamV3.h"
#include "TTreeStream.h"

ClassImp(AliTPC) 
//_____________________________________________________________________________
  AliTPC::AliTPC():AliDetector(),
		   fDefaults(0),
		   fSens(0),
		   fNsectors(0),
		   fDigitsArray(0),
		   fTPCParam(0),
		   fTrackHits(0),
		   fHitType(0),
		   fDigitsSwitch(0),
		   fSide(0),
                   fPrimaryIonisation(0),
		   fNoiseDepth(0),
		   fNoiseTable(0),
		   fCurrentNoise(0),
		   fActiveSectors(0),
                   fGainFactor(1.),
                   fDebugStreamer(0),
                   fLHCclockPhaseSw(0)

{
  //
  // Default constructor
  //
  fIshunt   = 0;
  for(Int_t i=0;i<4;i++) fCurrentIndex[i]=0;
 
  //  fTrackHitsOld = 0;   
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  fHitType = 4; // ROOT containers
#else
  fHitType = 2; //default CONTAINERS - based on ROOT structure
#endif 
}
 
//_____________________________________________________________________________
AliTPC::AliTPC(const char *name, const char *title)
  : AliDetector(name,title),
                   fDefaults(0),
		   fSens(0),
		   fNsectors(0),
		   fDigitsArray(0),
		   fTPCParam(0),
		   fTrackHits(0),
		   fHitType(0),
		   fDigitsSwitch(0),
		   fSide(0),
                   fPrimaryIonisation(0),
		   fNoiseDepth(0),
		   fNoiseTable(0),
		   fCurrentNoise(0),
                   fActiveSectors(0),
                   fGainFactor(1.),
                   fDebugStreamer(0),
                   fLHCclockPhaseSw(0)
                  
{
  //
  // Standard constructor
  //

  //
  // Initialise arrays of hits and digits 
  fHits     = new TClonesArray("AliTPChit",  176);
  gAlice->GetMCApp()->AddHitList(fHits); 
  //
  fTrackHits = new AliTPCTrackHitsV2;  
  fTrackHits->SetHitPrecision(0.002);
  fTrackHits->SetStepPrecision(0.003);  
  fTrackHits->SetMaxDistance(100);

  //fTrackHitsOld = new AliTPCTrackHits;  //MI - 13.09.2000
  //fTrackHitsOld->SetHitPrecision(0.002);
  //fTrackHitsOld->SetStepPrecision(0.003);  
  //fTrackHitsOld->SetMaxDistance(100); 


#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  fHitType = 4; // ROOT containers
#else
  fHitType = 2;
#endif

  for(Int_t i=0;i<4;i++) fCurrentIndex[i]=0;

  //
  fIshunt     =  0;
  //
  // Initialise color attributes
  //PH SetMarkerColor(kYellow);

  //
  //  Set TPC parameters
  //


  if (!strcmp(title,"Default")) {       
    //fTPCParam = new AliTPCParamSR;
    fTPCParam = AliTPCcalibDB::Instance()->GetParameters();
  } else {
    AliWarning("In Config.C you must set non-default parameters.");
    fTPCParam=0;
  }

}
void AliTPC::CreateDebugStremer(){
  //
  // Create Debug streamer to check simulation
  // 
  fDebugStreamer = new TTreeSRedirector("TPCSimdebug.root");
}
//_____________________________________________________________________________
AliTPC::~AliTPC()
{
  //
  // TPC destructor
  //

  fIshunt   = 0;
  delete fHits;
  delete fDigits;
  //delete fTPCParam;
  delete fTrackHits; //MI 15.09.2000
  //  delete fTrackHitsOld; //MI 10.12.2001
  
  fDigitsArray = 0x0;
  delete [] fNoiseTable;
  delete [] fActiveSectors;
  if (fDebugStreamer) delete fDebugStreamer;
}

//_____________________________________________________________________________
void AliTPC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a hit to the list
  //
  if (fHitType&1){
    TClonesArray &lhits = *fHits;
    new(lhits[fNhits++]) AliTPChit(fIshunt,track,vol,hits);
  }
  if (fHitType>1)
    AddHit2(track,vol,hits);
}

//_____________________________________________________________________________
void AliTPC::CreateMaterials()
{
  //-----------------------------------------------
  // Create Materials for for TPC simulations
  //-----------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

   Int_t iSXFLD=((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  Float_t sXMGMX=((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();

  Float_t amat[5]; // atomic numbers
  Float_t zmat[5]; // z
  Float_t wmat[5]; // proportions

  Float_t density;
 


  //***************** Gases *************************

 
  //--------------------------------------------------------------
  // gases - air and CO2
  //--------------------------------------------------------------

  // CO2

  amat[0]=12.011;
  amat[1]=15.9994;

  zmat[0]=6.;
  zmat[1]=8.;

  wmat[0]=0.2729;
  wmat[1]=0.7271;

  density=0.001754609;


  AliMixture(10,"CO2",amat,zmat,density,2,wmat);
  //
  // Air
  //
  amat[0]=15.9994;
  amat[1]=14.007;
  //
  zmat[0]=8.;
  zmat[1]=7.;
  //
  wmat[0]=0.233;
  wmat[1]=0.767;
  //
  density=0.001205;

  AliMixture(11,"Air",amat,zmat,density,2,wmat);
  
  //----------------------------------------------------------------
  // drift gases 
  //----------------------------------------------------------------

  //
  // Drift gases 1 - nonsensitive, 2 - sensitive
  // Ne-CO2 (90-10) (volume) values at 20deg and 1 atm.
  // rho(Ne) = 0.839 g/cm^3, rho(CO2) = 1.842 g/cm^3
  

  amat[0]= 20.18;
  amat[1]=12.011;
  amat[2]=15.9994;
  // amat[3]=14.007;

  zmat[0]= 10.; 
  zmat[1]=6.;
  zmat[2]=8.;
  // zmat[3]=7.;

  //wmat[0]=0.756992632;
  wmat[0]=0.8038965;
  //wmat[1]=0.056235789;
  wmat[1]= 0.053519;
  //wmat[2]=0.128469474;
  wmat[2]= 0.1425743;
  // wmat[3]=0.058395789;
 
  density=0.0009393;

  AliMixture(12,"Ne-CO2-1",amat,zmat,density,3,wmat);
  AliMixture(13,"Ne-CO2-2",amat,zmat,density,3,wmat);
  AliMixture(35,"Ne-CO2-3",amat,zmat,density,3,wmat);
  //----------------------------------------------------------------------
  //               solid materials
  //----------------------------------------------------------------------


  // Kevlar C14H22O2N2

  amat[0] = 12.011;
  amat[1] = 1.;
  amat[2] = 15.999;
  amat[3] = 14.006;

  zmat[0] = 6.;
  zmat[1] = 1.;
  zmat[2] = 8.;
  zmat[3] = 7.;

  wmat[0] = 14.;
  wmat[1] = 22.;
  wmat[2] = 2.;
  wmat[3] = 2.;

  density = 1.45;

  AliMixture(14,"Kevlar",amat,zmat,density,-4,wmat);  

  // NOMEX

  amat[0] = 12.011;
  amat[1] = 1.;
  amat[2] = 15.999;
  amat[3] = 14.006;

  zmat[0] = 6.;
  zmat[1] = 1.;
  zmat[2] = 8.;
  zmat[3] = 7.;

  wmat[0] = 14.;
  wmat[1] = 22.;
  wmat[2] = 2.;
  wmat[3] = 2.;

  density = 0.029;
 
  AliMixture(15,"NOMEX",amat,zmat,density,-4,wmat);

  // Makrolon C16H18O3

  amat[0] = 12.011;
  amat[1] = 1.;
  amat[2] = 15.999;

  zmat[0] = 6.;
  zmat[1] = 1.;
  zmat[2] = 8.;

  wmat[0] = 16.;
  wmat[1] = 18.;
  wmat[2] = 3.;
  
  density = 1.2;

  AliMixture(16,"Makrolon",amat,zmat,density,-3,wmat);

  // Tedlar C2H3F

  amat[0] = 12.011;
  amat[1] = 1.;
  amat[2] = 18.998;

  zmat[0] = 6.;
  zmat[1] = 1.;
  zmat[2] = 9.;

  wmat[0] = 2.;
  wmat[1] = 3.; 
  wmat[2] = 1.;

  density = 1.71;

  AliMixture(17, "Tedlar",amat,zmat,density,-3,wmat);  
  
  // Mylar C5H4O2

  amat[0]=12.011;
  amat[1]=1.;
  amat[2]=15.9994;

  zmat[0]=6.;
  zmat[1]=1.;
  zmat[2]=8.;

  wmat[0]=5.;
  wmat[1]=4.;
  wmat[2]=2.; 

  density = 1.39;
  
  AliMixture(18, "Mylar",amat,zmat,density,-3,wmat); 
  // material for "prepregs"
  // Epoxy - C14 H20 O3
  // Quartz SiO2
  // Carbon C
  // prepreg1 60% C-fiber, 40% epoxy (vol)
  amat[0]=12.011;
  amat[1]=1.;
  amat[2]=15.994;

  zmat[0]=6.;
  zmat[1]=1.;
  zmat[2]=8.;

  wmat[0]=0.923;
  wmat[1]=0.023;
  wmat[2]=0.054;

  density=1.859;

  AliMixture(19, "Prepreg1",amat,zmat,density,3,wmat);

  //prepreg2 60% glass-fiber, 40% epoxy

  amat[0]=12.01;
  amat[1]=1.;
  amat[2]=15.994;
  amat[3]=28.086;

  zmat[0]=6.;
  zmat[1]=1.;
  zmat[2]=8.;
  zmat[3]=14.;

  wmat[0]=0.194;
  wmat[1]=0.023;
  wmat[2]=0.443;
  wmat[3]=0.34;

  density=1.82;

  AliMixture(20, "Prepreg2",amat,zmat,density,4,wmat);

  //prepreg3 50% glass-fiber, 50% epoxy

  amat[0]=12.01;
  amat[1]=1.;
  amat[2]=15.994;
  amat[3]=28.086;

  zmat[0]=6.;
  zmat[1]=1.;
  zmat[2]=8.;
  zmat[3]=14.;

  wmat[0]=0.257;
  wmat[1]=0.03;
  wmat[2]=0.412;
  wmat[3]=0.3;

  density=1.725;

  AliMixture(21, "Prepreg3",amat,zmat,density,4,wmat);

  // G10 60% SiO2 40% epoxy

  amat[0]=12.01;
  amat[1]=1.;
  amat[2]=15.994;
  amat[3]=28.086;

  zmat[0]=6.;
  zmat[1]=1.;
  zmat[2]=8.;
  zmat[3]=14.;

  wmat[0]=0.194;
  wmat[1]=0.023;
  wmat[2]=0.443;
  wmat[3]=0.340;

  density=1.7;

  AliMixture(22, "G10",amat,zmat,density,4,wmat);
 
  // Al

  amat[0] = 26.98;
  zmat[0] = 13.;

  density = 2.7;

  AliMaterial(23,"Al",amat[0],zmat[0],density,999.,999.);

  // Si (for electronics

  amat[0] = 28.086;
  zmat[0] = 14.;

  density = 2.33;

  AliMaterial(24,"Si",amat[0],zmat[0],density,999.,999.);

  // Cu

  amat[0] = 63.546;
  zmat[0] = 29.;

  density = 8.96;

  AliMaterial(25,"Cu",amat[0],zmat[0],density,999.,999.);

  // brass

  amat[0] = 63.546;
  zmat[0] = 29.;
  //
  amat[1]= 65.409;
  zmat[1]= 30.;
  //
  wmat[0]= 0.6;
  wmat[1]= 0.4;

  //
  density = 8.23;
  
 
  //
  AliMixture(33,"Brass",amat,zmat,density,2,wmat);
  
  // Epoxy - C14 H20 O3
 
  amat[0]=12.011;
  amat[1]=1.;
  amat[2]=15.9994;

  zmat[0]=6.;
  zmat[1]=1.;
  zmat[2]=8.;

  wmat[0]=14.;
  wmat[1]=20.;
  wmat[2]=3.;

  density=1.25;

  AliMixture(26,"Epoxy",amat,zmat,density,-3,wmat);
  //
  // epoxy film - 90% epoxy, 10% glass fiber 
  //
  amat[0]=12.01;
  amat[1]=1.;
  amat[2]=15.994;
  amat[3]=28.086;

  zmat[0]=6.;
  zmat[1]=1.;
  zmat[2]=8.;
  zmat[3]=14.;

  wmat[0]=0.596;
  wmat[1]=0.071;
  wmat[2]=0.257;
  wmat[3]=0.076;


  density=1.345;

  AliMixture(34, "Epoxy-film",amat,zmat,density,4,wmat);

  // Plexiglas  C5H8O2

  amat[0]=12.011;
  amat[1]=1.;
  amat[2]=15.9994;

  zmat[0]=6.;
  zmat[1]=1.;
  zmat[2]=8.;

  wmat[0]=5.;
  wmat[1]=8.;
  wmat[2]=2.;

  density=1.18;

  AliMixture(27,"Plexiglas",amat,zmat,density,-3,wmat);

  // Carbon

  amat[0]=12.011;
  zmat[0]=6.;
  density= 2.265;

  AliMaterial(28,"C",amat[0],zmat[0],density,999.,999.);

  // Fe (steel for the inner heat screen)
 
  amat[0]=55.845;

  zmat[0]=26.;

  density=7.87;

  AliMaterial(29,"Fe",amat[0],zmat[0],density,999.,999.);
  //
  // Peek - (C6H4-O-OC6H4-O-C6H4-CO)n
  amat[0]=12.011;
  amat[1]=1.;
  amat[2]=15.9994;

  zmat[0]=6.;
  zmat[1]=1.;
  zmat[2]=8.;

  wmat[0]=19.;
  wmat[1]=12.;
  wmat[2]=3.;
  //
  density=1.3;
  //
  AliMixture(30,"Peek",amat,zmat,density,-3,wmat);  
  //
  //  Ceramics - Al2O3
  //
  amat[0] = 26.98;
  amat[1]= 15.9994;

  zmat[0] = 13.;
  zmat[1]=8.;
 
  wmat[0]=2.;
  wmat[1]=3.;
 
  density = 3.97;

  AliMixture(31,"Alumina",amat,zmat,density,-2,wmat);   

  //
  // liquids
  //

  // water

  amat[0]=1.;
  amat[1]=15.9994;

  zmat[0]=1.;
  zmat[1]=8.;

  wmat[0]=2.;
  wmat[1]=1.;

  density=1.;

  AliMixture(32,"Water",amat,zmat,density,-2,wmat);  

 
  //----------------------------------------------------------
  // tracking media for gases
  //----------------------------------------------------------

  AliMedium(0, "Air", 11, 0, iSXFLD, sXMGMX, 10., 999., .1, .01, .1);
  AliMedium(1, "Ne-CO2-1", 12, 0, iSXFLD, sXMGMX, 10., 999.,.1,.001, .001);
  AliMedium(2, "Ne-CO2-2", 13, 1, iSXFLD, sXMGMX, 10., 999.,.1,.001, .001);
  AliMedium(3,"CO2",10,0, iSXFLD, sXMGMX, 10., 999.,.1, .001, .001); 
  AliMedium(20, "Ne-CO2-3", 35, 1, iSXFLD, sXMGMX, 10., 999.,.1,.001, .001);
  //-----------------------------------------------------------  
  // tracking media for solids
  //-----------------------------------------------------------
  
  AliMedium(4,"Al",23,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(5,"Kevlar",14,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(6,"Nomex",15,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(7,"Makrolon",16,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(8,"Mylar",18,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(9,"Tedlar",17,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  //
  AliMedium(10,"Prepreg1",19,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(11,"Prepreg2",20,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(12,"Prepreg3",21,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(13,"Epoxy",26,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);

  AliMedium(14,"Cu",25,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(15,"Si",24,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(16,"G10",22,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(17,"Plexiglas",27,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(18,"Steel",29,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001); 
  AliMedium(19,"Peek",30,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(21,"Alumina",31,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);    
  AliMedium(22,"Water",32,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(23,"Brass",33,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(24,"Epoxyfm",34,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);  
}

void AliTPC::GenerNoise(Int_t tablesize)
{
  //
  //Generate table with noise
  //
  if (fTPCParam==0) {
    // error message
    fNoiseDepth=0;
    return;
  }
  if (fNoiseTable)  delete[] fNoiseTable;
  fNoiseTable = new Float_t[tablesize];
  fNoiseDepth = tablesize; 
  fCurrentNoise =0; //!index of the noise in  the noise table 
  
  Float_t norm = fTPCParam->GetNoise()*fTPCParam->GetNoiseNormFac();
  for (Int_t i=0;i<tablesize;i++) fNoiseTable[i]= gRandom->Gaus(0,norm);      
}

Float_t AliTPC::GetNoise()
{
  // get noise from table
  //  if ((fCurrentNoise%10)==0) 
  //  fCurrentNoise= gRandom->Rndm()*fNoiseDepth;
  if (fCurrentNoise>=fNoiseDepth) fCurrentNoise=0;
  return fNoiseTable[fCurrentNoise++];
  //gRandom->Gaus(0, fTPCParam->GetNoise()*fTPCParam->GetNoiseNormFac()); 
}


Bool_t  AliTPC::IsSectorActive(Int_t sec) const
{
  //
  // check if the sector is active
  if (!fActiveSectors) return kTRUE;
  else return fActiveSectors[sec]; 
}

void    AliTPC::SetActiveSectors(Int_t * sectors, Int_t n)
{
  // activate interesting sectors
  SetTreeAddress();//just for security
  if (!fActiveSectors) fActiveSectors = new Bool_t[fTPCParam->GetNSector()];
  for (Int_t i=0;i<fTPCParam->GetNSector();i++) fActiveSectors[i]=kFALSE;
  for (Int_t i=0;i<n;i++) 
    if ((sectors[i]>=0) && sectors[i]<fTPCParam->GetNSector())  fActiveSectors[sectors[i]]=kTRUE;
    
}

void    AliTPC::SetActiveSectors(Int_t flag)
{
  //
  // activate sectors which were hitted by tracks 
  //loop over tracks
  SetTreeAddress();//just for security
  if (fHitType==0) return;  // if Clones hit - not short volume ID information
  if (!fActiveSectors) fActiveSectors = new Bool_t[fTPCParam->GetNSector()];
  if (flag) {
    for (Int_t i=0;i<fTPCParam->GetNSector();i++) fActiveSectors[i]=kTRUE;
    return;
  }
  for (Int_t i=0;i<fTPCParam->GetNSector();i++) fActiveSectors[i]=kFALSE;
  //TBranch * branch=0;
  if (fLoader->TreeH() == 0x0)
   {
     AliFatal("Can not find TreeH in folder");
     return;
   }
  //if (fHitType>1) branch = fLoader->TreeH()->GetBranch("TPC2");
  if (fHitType>1) fLoader->TreeH()->GetBranch("TPC2");
  //else branch = fLoader->TreeH()->GetBranch("TPC");
  else fLoader->TreeH()->GetBranch("TPC");
  Stat_t ntracks = fLoader->TreeH()->GetEntries();
  // loop over all hits
  AliDebug(1,Form("Got %d tracks", (Int_t) ntracks));
  
  for(Int_t track=0;track<ntracks;track++) {
    ResetHits();
    //
    if (fTrackHits && fHitType&4) {
      TBranch * br1 = fLoader->TreeH()->GetBranch("fVolumes");
      TBranch * br2 = fLoader->TreeH()->GetBranch("fNVolumes");
      br1->GetEvent(track);
      br2->GetEvent(track);
      Int_t *volumes = fTrackHits->GetVolumes();
      for (Int_t j=0;j<fTrackHits->GetNVolumes(); j++) {
	if (volumes[j]>-1 && volumes[j]<fTPCParam->GetNSector()) {
	  fActiveSectors[volumes[j]]=kTRUE;
	}
	else {
	    AliError(Form("Volume %d -> sector number %d is outside (0..%d)",
			  j,
			  volumes[j],
			  fTPCParam->GetNSector()));
	}
      }
    }
    
    //
//     if (fTrackHitsOld && fHitType&2) {
//       TBranch * br = fLoader->TreeH()->GetBranch("fTrackHitsInfo");
//       br->GetEvent(track);
//       AliObjectArray * ar = fTrackHitsOld->fTrackHitsInfo;
//       for (UInt_t j=0;j<ar->GetSize();j++){
// 	fActiveSectors[((AliTrackHitsInfo*)ar->At(j))->fVolumeID] =kTRUE;
//       } 
//     }    
  }
}  




//_____________________________________________________________________________
void AliTPC::Digits2Raw()
{
// convert digits of the current event to raw data

  static const Int_t kThreshold = 0;

  fLoader->LoadDigits();
  TTree* digits = fLoader->TreeD();
  if (!digits) {
    AliError("No digits tree");
    return;
  }

  //
  AliSimDigits digarr;
  AliSimDigits* digrow = &digarr;
  digits->GetBranch("Segment")->SetAddress(&digrow);

  const char* fileName = "AliTPCDDL.dat";
  AliTPCBuffer* buffer  = new AliTPCBuffer(fileName);
  //Verbose level
  // 0: Silent
  // 1: cout messages
  // 2: txt files with digits 
  //BE CAREFUL, verbose level 2 MUST be used only for debugging and
  //it is highly suggested to use this mode only for debugging digits files
  //reasonably small, because otherwise the size of the txt files can reach
  //quickly several MB wasting time and disk space.
  buffer->SetVerbose(0);

  Int_t nEntries = Int_t(digits->GetEntries());
  Int_t previousSector = -1;
  Int_t subSector = 0;
  for (Int_t i = 0; i < nEntries; i++) {
    digits->GetEntry(i);
    Int_t sector, row;
    fTPCParam->AdjustSectorRow(digarr.GetID(), sector, row);
    if(previousSector != sector) {
      subSector = 0;
      previousSector = sector;
    }

    if (sector < 36) { //inner sector [0;35]
      if (row != 30) {
	//the whole row is written into the output file
	buffer->WriteRowBinary(kThreshold, digrow, 0, 0, 0, 
			       sector, subSector, row);
      } else {
	//only the pads in the range [37;48] are written into the output file
	buffer->WriteRowBinary(kThreshold, digrow, 37, 48, 1, 
			       sector, subSector, row);
	subSector = 1;
	//only the pads outside the range [37;48] are written into the output file
	buffer->WriteRowBinary(kThreshold, digrow, 37, 48, 2, 
			       sector, subSector, row);
      }//end else

    } else { //outer sector [36;71]
      if (row == 54) subSector = 2;
      if ((row != 27) && (row != 76)) {
	buffer->WriteRowBinary(kThreshold, digrow, 0, 0, 0,
			       sector, subSector, row);
      } else if (row == 27) {
	//only the pads outside the range [43;46] are written into the output file
	buffer->WriteRowBinary(kThreshold, digrow, 43, 46, 2,
				 sector, subSector, row);
	subSector = 1;
	//only the pads in the range [43;46] are written into the output file
	buffer->WriteRowBinary(kThreshold, digrow, 43, 46, 1,
				 sector, subSector, row);
      } else if (row == 76) {
	//only the pads outside the range [33;88] are written into the output file
	buffer->WriteRowBinary(kThreshold, digrow, 33, 88, 2,
			       sector, subSector, row);
	subSector = 3;
	//only the pads in the range [33;88] are written into the output file
	buffer->WriteRowBinary(kThreshold, digrow, 33, 88, 1,
				 sector, subSector, row);
      }
    }//end else
  }//end for

  delete buffer;
  fLoader->UnloadDigits();

  AliTPCDDLRawData rawWriter;
  rawWriter.SetVerbose(0);

  rawWriter.RawData(fileName);
  gSystem->Unlink(fileName);

}


//_____________________________________________________________________________
Bool_t AliTPC::Raw2SDigits(AliRawReader* rawReader){
  // Converts the TPC raw data into summable digits
  // The method is used for merging simulated and
  // real data events
  if (fLoader->TreeS() == 0x0 ) {
    fLoader->MakeTree("S");
  }

  if(fDefaults == 0) SetDefaults();  // check if the parameters are set

  //setup TPCDigitsArray 
  if(GetDigitsArray()) delete GetDigitsArray();

  AliTPCDigitsArray *arr = new AliTPCDigitsArray; 
  arr->SetClass("AliSimDigits");
  arr->Setup(fTPCParam);
  arr->MakeTree(fLoader->TreeS());

  SetDigitsArray(arr);

  // set zero suppression to "0"
  fTPCParam->SetZeroSup(0);

  // Loop over sectors
  const Int_t kmaxTime = fTPCParam->GetMaxTBin();
  const Int_t kNIS = fTPCParam->GetNInnerSector();
  const Int_t kNOS = fTPCParam->GetNOuterSector();
  const Int_t kNS = kNIS + kNOS;

  // Setup storage
  AliTPCROC * roc = AliTPCROC::Instance();
  Int_t nRowsMax = roc->GetNRows(roc->GetNSector()-1);
  Int_t nPadsMax = roc->GetNPads(roc->GetNSector()-1,nRowsMax-1);
  Short_t** allBins = new Short_t*[nRowsMax];
  for (Int_t iRow = 0; iRow < nRowsMax; iRow++) {
    Int_t maxBin = kmaxTime*nPadsMax;
    allBins[iRow] = new Short_t[maxBin];
    memset(allBins[iRow],0,sizeof(Short_t)*maxBin);
  }

  for(Int_t iSector = 0; iSector < kNS; iSector++) {
    
    Int_t nRows = fTPCParam->GetNRow(iSector);
    Int_t nDDLs = 0, indexDDL = 0;
    if (iSector < kNIS) {
      nDDLs = 2;
      indexDDL = iSector * 2;
    }
    else {
      nDDLs = 4;
      indexDDL = (iSector-kNIS) * 4 + kNIS * 2;
    }

    // Load the raw data for corresponding DDLs
    rawReader->Reset();

    AliTPCAltroMapping** mapping =AliTPCcalibDB::Instance()->GetMapping();
    AliTPCRawStreamV3 input(rawReader,(AliAltroMapping**)mapping);
    rawReader->Select("TPC",indexDDL,indexDDL+nDDLs-1);

    // Clean storage
    for (Int_t iRow = 0; iRow < nRowsMax; iRow++) {
      Int_t maxBin = kmaxTime*nPadsMax;
      memset(allBins[iRow],0,sizeof(Short_t)*maxBin);
    }

    // Begin loop over altro data
    while (input.NextDDL()) {

      if (input.GetSector() != iSector)
	AliFatal(Form("Sector index mismatch ! Expected (%d), but got (%d) !",iSector,input.GetSector()));

      //loop over pads
      while ( input.NextChannel() ) {

        Int_t iRow = input.GetRow();
        if (iRow < 0 || iRow >= nRows)
          AliFatal(Form("Pad-row index (%d) outside the range (%d -> %d) !",
                        iRow, 0, nRows -1));
        Int_t iPad = input.GetPad();

        Int_t maxPad = fTPCParam->GetNPads(iSector,iRow);

        if (iPad < 0 || iPad >= maxPad)
          AliFatal(Form("Pad index (%d) outside the range (%d -> %d) !",
                        iPad, 0, maxPad -1));

        //loop over bunches
        while ( input.NextBunch() ){
          Int_t  startTbin    = (Int_t)input.GetStartTimeBin();
          Int_t  bunchlength  = (Int_t)input.GetBunchLength();
          const UShort_t *sig = input.GetSignals();
          for (Int_t iTime = 0; iTime<bunchlength; iTime++){
            Int_t iTimeBin=startTbin-iTime;
            if ( iTimeBin < 0 || iTimeBin >= kmaxTime) {
              continue;
              //AliFatal(Form("Timebin index (%d) outside the range (%d -> %d) !",
              //               iTimeBin, 0, kmaxTime -1));
            }

            Int_t maxBin = kmaxTime*maxPad;
            if (((iPad*kmaxTime+iTimeBin) >= maxBin) ||
                ((iPad*kmaxTime+iTimeBin) < 0))
              AliFatal(Form("Index outside the allowed range"
                            " Sector=%d Row=%d Pad=%d Timebin=%d"
                            " (Max.index=%d)",iSector,iRow,iPad,iTimeBin,maxBin));
            allBins[iRow][iPad*kmaxTime+iTimeBin] = sig[iTime];
          }
        }
      } // End loop over altro data
    }

    // Now fill the digits array
    if (fDigitsArray->GetTree()==0) {
      AliFatal("Tree not set in fDigitsArray");
    }

    for (Int_t iRow = 0; iRow < nRows; iRow++) {
      AliDigits * dig = fDigitsArray->CreateRow(iSector,iRow);
      Int_t maxPad = fTPCParam->GetNPads(iSector,iRow);
      for(Int_t iPad = 0; iPad < maxPad; iPad++) {
	for(Int_t iTimeBin = 0; iTimeBin < kmaxTime; iTimeBin++) {
	  Short_t q = allBins[iRow][iPad*kmaxTime + iTimeBin];
	  if (q <= 0) continue;
	  q *= 16;
	  dig->SetDigitFast((Short_t)q,iTimeBin,iPad);
	}
      }
      fDigitsArray->StoreRow(iSector,iRow);
      Int_t ndig = dig->GetDigitSize(); 
	
      AliDebug(10,
	       Form("*** Sector, row, compressed digits %d %d %d ***\n",
		    iSector,iRow,ndig));        
 	
      fDigitsArray->ClearRow(iSector,iRow);  

    } // end of the sector digitization
  }
  // get LHC clock phase from the digits tree

  TParameter<float> *ph; 
  Float_t phase;
  TTree *digtree = fLoader->TreeD();
  //
  if(digtree){ // if TreeD exists
    ph = (TParameter<float>*)digtree->GetUserInfo()->FindObject("lhcphase0");
    phase = ph->GetVal();
  }
  else{ //TreeD does not exist
    phase = 0.; 
  }
    //
    // store lhc clock phase in S-digits tree
    //
    fLoader->TreeS()->GetUserInfo()->Add(new TParameter<float>("lhcphase0",phase));
   //
   fLoader->WriteSDigits("OVERWRITE");

  if(GetDigitsArray()) delete GetDigitsArray();
  SetDigitsArray(0x0);

  // cleanup storage
  for (Int_t iRow = 0; iRow < nRowsMax; iRow++)
    delete [] allBins[iRow];
  delete [] allBins;

  return kTRUE;
}

//______________________________________________________________________
AliDigitizer* AliTPC::CreateDigitizer(AliDigitizationInput* digInput) const
{
  return new AliTPCDigitizer(digInput);
}
//__
void AliTPC::SDigits2Digits2(Int_t /*eventnumber*/)  
{
  //create digits from summable digits
  GenerNoise(500000); //create teble with noise

  //conect tree with sSDigits
  TTree *t = fLoader->TreeS();

  if (t == 0x0) {
    fLoader->LoadSDigits("READ");
    t = fLoader->TreeS();
    if (t == 0x0) {
      AliError("Can not get input TreeS");
      return;
    }
  }
  
  if (fLoader->TreeD() == 0x0) fLoader->MakeTree("D");
  
  AliSimDigits digarr, *dummy=&digarr;
  TBranch* sdb = t->GetBranch("Segment");
  if (sdb == 0x0) {
    AliError("Can not find branch with segments in TreeS.");
    return;
  }  

  sdb->SetAddress(&dummy);
      
  Stat_t nentries = t->GetEntries();

  // set zero suppression

  fTPCParam->SetZeroSup(2);

  // get zero suppression

  Int_t zerosup = fTPCParam->GetZeroSup();

  //make tree with digits 
  
  AliTPCDigitsArray *arr = new AliTPCDigitsArray; 
  arr->SetClass("AliSimDigits");
  arr->Setup(fTPCParam);
  arr->MakeTree(fLoader->TreeD());
  
  AliTPCParam * par = fTPCParam;

  //Loop over segments of the TPC

  for (Int_t n=0; n<nentries; n++) {
    t->GetEvent(n);
    Int_t sec, row;
    if (!par->AdjustSectorRow(digarr.GetID(),sec,row)) {
      AliWarning(Form("Invalid segment ID ! %d",digarr.GetID()));
      continue;
    }
    if (!IsSectorActive(sec)) continue;
    
    AliSimDigits * digrow =(AliSimDigits*) arr->CreateRow(sec,row);
    Int_t nrows = digrow->GetNRows();
    Int_t ncols = digrow->GetNCols();

    digrow->ExpandBuffer();
    digarr.ExpandBuffer();
    digrow->ExpandTrackBuffer();
    digarr.ExpandTrackBuffer();

    
    Short_t * pamp0 = digarr.GetDigits();
    Int_t   * ptracks0 = digarr.GetTracks();
    Short_t * pamp1 = digrow->GetDigits();
    Int_t   * ptracks1 = digrow->GetTracks();
    Int_t  nelems =nrows*ncols;
    Int_t saturation = fTPCParam->GetADCSat() - 1;
    //use internal structure of the AliDigits - for speed reason
    //if you cahnge implementation
    //of the Alidigits - it must be rewriten -
    for (Int_t i= 0; i<nelems; i++){
      Float_t q = TMath::Nint(Float_t(*pamp0)/16.+GetNoise());
      if (q>zerosup){
	if (q>saturation) q=saturation;      
	*pamp1=(Short_t)q;

	ptracks1[0]=ptracks0[0];	
	ptracks1[nelems]=ptracks0[nelems];
	ptracks1[2*nelems]=ptracks0[2*nelems];
      }
      pamp0++;
      pamp1++;
      ptracks0++;
      ptracks1++;	 
    }

    arr->StoreRow(sec,row);
    arr->ClearRow(sec,row);   
  }  

    
  //write results
  fLoader->WriteDigits("OVERWRITE");
   
  delete arr;
}
//__________________________________________________________________
void AliTPC::SetDefaults(){
  //
  // setting the defaults
  //
   
  // Set response functions

  //
  AliRunLoader* rl = (AliRunLoader*)fLoader->GetEventFolder()->FindObject(AliRunLoader::GetRunLoaderName());
  rl->CdGAFile();
  //AliTPCParamSR *param=(AliTPCParamSR*)gDirectory->Get("75x40_100x60");
  //gDirectory->Get("75x40_100x60");
  AliTPCParamSR *param = (AliTPCParamSR*)AliTPCcalibDB::Instance()->GetParameters();
  if(!param){
    AliFatal("No TPC parameters found");
    return;
  }
  if (!param->IsGeoRead()){
      //
      // read transformation matrices for gGeoManager
      //
      param->ReadGeoMatrices();
    }



  AliTPCPRF2D    * prfinner   = new AliTPCPRF2D;
  AliTPCPRF2D    * prfouter1   = new AliTPCPRF2D;
  AliTPCPRF2D    * prfouter2   = new AliTPCPRF2D;  

  
  //AliTPCRF1D     * rf    = new AliTPCRF1D(kTRUE);
  //rf->SetGauss(param->GetZSigma(),param->GetZWidth(),1.);
  //rf->SetOffset(3*param->GetZSigma());
  //rf->Update();
  //
  // Use gamma 4
  //
  char  strgamma4[1000];
  //sprintf(strgamma4,"AliTPCRF1D::Gamma4((x-0.135+%f)*%f,55,160)",3*param->GetZSigma(), 1000000000*param->GetTSample()/param->GetZWidth());
  
  snprintf(strgamma4,1000,"AliTPCRF1D::Gamma4((x-0.135+%f)*%f,55,160)",3*param->GetZSigma(), 1000000000*param->GetTSample()/param->GetZWidth());
  TF1 * fgamma4 = new TF1("fgamma4",strgamma4, -1,1);
  AliTPCRF1D     * rf    = new AliTPCRF1D(kTRUE,1000);
  rf->SetParam(fgamma4,param->GetZWidth(), 1,0.2);
  rf->SetOffset(3*param->GetZSigma()); 
  rf->Update();

  TDirectory *savedir=gDirectory;
  TFile *f=TFile::Open("$ALICE_ROOT/TPC/AliTPCprf2d.root");
  if (!f->IsOpen()) 
    AliFatal("Can't open $ALICE_ROOT/TPC/AliTPCprf2d.root !");

  TString s;
  prfinner->Read("prf_07504_Gati_056068_d02");
  //PH Set different names
  s=prfinner->GetGRF()->GetName();
  s+="in";
  prfinner->GetGRF()->SetName(s.Data());

  prfouter1->Read("prf_10006_Gati_047051_d03");
  s=prfouter1->GetGRF()->GetName();
  s+="out1";
  prfouter1->GetGRF()->SetName(s.Data());

  prfouter2->Read("prf_15006_Gati_047051_d03");  
  s=prfouter2->GetGRF()->GetName();
  s+="out2";
  prfouter2->GetGRF()->SetName(s.Data());

  f->Close();
  savedir->cd();

  param->SetInnerPRF(prfinner);
  param->SetOuter1PRF(prfouter1); 
  param->SetOuter2PRF(prfouter2);
  param->SetTimeRF(rf);

  // set fTPCParam

  SetParam(param);


  fDefaults = 1;

}
//__________________________________________________________________  
void AliTPC::Hits2Digits()  
{
  //
  // creates digits from hits
  //
  if (!fTPCParam->IsGeoRead()){
    //
    // read transformation matrices for gGeoManager
    //
    fTPCParam->ReadGeoMatrices();
  }

  fLoader->LoadHits("read");
  fLoader->LoadDigits("recreate");
  AliRunLoader* runLoader = fLoader->GetRunLoader(); 

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    //PH    runLoader->GetEvent(iEvent);
    Hits2Digits(iEvent);
  }

  fLoader->UnloadHits();
  fLoader->UnloadDigits();
} 
//__________________________________________________________________  
void AliTPC::Hits2Digits(Int_t eventnumber)  
{ 
 //----------------------------------------------------
 // Loop over all sectors for a single event
 //----------------------------------------------------
  AliRunLoader* rl = (AliRunLoader*)fLoader->GetEventFolder()->FindObject(AliRunLoader::GetRunLoaderName());
  rl->GetEvent(eventnumber);
  SetActiveSectors();   
  if (fLoader->TreeH() == 0x0) {
    if(fLoader->LoadHits()) {
      AliError("Can not load hits.");
    }
  }
  SetTreeAddress();
  
  if (fLoader->TreeD() == 0x0 ) {
    fLoader->MakeTree("D");
    if (fLoader->TreeD() == 0x0 ) {
      AliError("Can not get TreeD");
      return;
    }
  }

  if(fDefaults == 0) SetDefaults();  // check if the parameters are set
  GenerNoise(500000); //create teble with noise

  //setup TPCDigitsArray 

  if(GetDigitsArray()) delete GetDigitsArray();
  
  AliTPCDigitsArray *arr = new AliTPCDigitsArray; 
  arr->SetClass("AliSimDigits");
  arr->Setup(fTPCParam);

  arr->MakeTree(fLoader->TreeD());
  SetDigitsArray(arr);

  fDigitsSwitch=0; // standard digits
  // here LHC clock phase
  Float_t lhcph = 0.;
  switch (fLHCclockPhaseSw){
  case 0: 
    // no phase
    lhcph=0.;
    break;
  case 1:
    // random phase
    lhcph = (Int_t)(gRandom->Rndm()/0.25);    
    break;
  case 2:
    lhcph=0.;
    // not implemented yet
    break;
  }
  // adding phase to the TreeD user info 
  fLoader->TreeD()->GetUserInfo()->Add(new TParameter<float>("lhcphase0",lhcph));
  //
  for(Int_t isec=0;isec<fTPCParam->GetNSector();isec++) 
    if (IsSectorActive(isec)) {
      AliDebug(1,Form("Hits2Digits: Sector %d is active.",isec));
      Hits2DigitsSector(isec);
    }
    else {
      AliDebug(1,Form("Hits2Digits: Sector %d is NOT active.",isec));
    }
  
  fLoader->WriteDigits("OVERWRITE"); 
  
//this line prevents the crash in the similar one
//on the beginning of this method
//destructor attempts to reset the tree, which is deleted by the loader
//need to be redesign
  if(GetDigitsArray()) delete GetDigitsArray();
  SetDigitsArray(0x0);
  
}

//__________________________________________________________________
void AliTPC::Hits2SDigits2(Int_t eventnumber)  
{ 

  //-----------------------------------------------------------
  //   summable digits - 16 bit "ADC", no noise, no saturation
  //-----------------------------------------------------------

  //----------------------------------------------------
  // Loop over all sectors for a single event
  //----------------------------------------------------

  AliRunLoader* rl = fLoader->GetRunLoader();

  rl->GetEvent(eventnumber);
  if (fLoader->TreeH() == 0x0) {
    if(fLoader->LoadHits()) {
      AliError("Can not load hits.");
      return;
    }
  }
  SetTreeAddress();


  if (fLoader->TreeS() == 0x0 ) {
    fLoader->MakeTree("S");
  }
  
  if(fDefaults == 0) SetDefaults();
  
  GenerNoise(500000); //create table with noise
  //setup TPCDigitsArray 

  if(GetDigitsArray()) delete GetDigitsArray();

  
  AliTPCDigitsArray *arr = new AliTPCDigitsArray; 
  arr->SetClass("AliSimDigits");
  arr->Setup(fTPCParam);
  arr->MakeTree(fLoader->TreeS());

  SetDigitsArray(arr);

  fDigitsSwitch=1; // summable digits
  
    // set zero suppression to "0"
  // here LHC clock phase
  Float_t lhcph = 0.;
  switch (fLHCclockPhaseSw){
  case 0: 
    // no phase
    lhcph=0.;
    break;
  case 1:
    // random phase
    lhcph = (Int_t)(gRandom->Rndm()/0.25);    
    break;
  case 2:
    lhcph=0.;
    // not implemented yet
    break;
  }
  // adding phase to the TreeS user info 
  
  fLoader->TreeS()->GetUserInfo()->Add(new TParameter<float>("lhcphase0",lhcph));

  fTPCParam->SetZeroSup(0);

  for(Int_t isec=0;isec<fTPCParam->GetNSector();isec++) 
    if (IsSectorActive(isec)) {
      Hits2DigitsSector(isec);
    }

  fLoader->WriteSDigits("OVERWRITE");

//this line prevents the crash in the similar one
//on the beginning of this method
//destructor attempts to reset the tree, which is deleted by the loader
//need to be redesign
  if(GetDigitsArray()) delete GetDigitsArray();
  SetDigitsArray(0x0);
}
//__________________________________________________________________

void AliTPC::Hits2SDigits()  
{ 

  //-----------------------------------------------------------
  //   summable digits - 16 bit "ADC", no noise, no saturation
  //-----------------------------------------------------------

  if (!fTPCParam->IsGeoRead()){
    //
    // read transformation matrices for gGeoManager
    //
    fTPCParam->ReadGeoMatrices();
  }
  
  fLoader->LoadHits("read");
  fLoader->LoadSDigits("recreate");
  AliRunLoader* runLoader = fLoader->GetRunLoader(); 

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);
    SetTreeAddress();
    SetActiveSectors();
    Hits2SDigits2(iEvent);
  }
  
  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
  if (fDebugStreamer) {
    delete fDebugStreamer;
    fDebugStreamer=0;
  }    
}
//_____________________________________________________________________________

void AliTPC::Hits2DigitsSector(Int_t isec)
{
  //-------------------------------------------------------------------
  // TPC conversion from hits to digits.
  //------------------------------------------------------------------- 

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  //-------------------------------------------------------
  //  Get the access to the track hits
  //-------------------------------------------------------

  // check if the parameters are set - important if one calls this method
  // directly, not from the Hits2Digits

  if(fDefaults == 0) SetDefaults();

  TTree *tH = fLoader->TreeH(); // pointer to the hits tree
  if (tH == 0x0) {
    AliFatal("Can not find TreeH in folder");
    return;
  }

  Stat_t ntracks = tH->GetEntries();

    Int_t nrows =fTPCParam->GetNRow(isec);

    TObjArray **row=new TObjArray* [nrows+2]; // 2 extra rows for cross talk
    for(Int_t j=0;j<nrows+2;j++) row[j]=0;
    
    MakeSector(isec,nrows,tH,ntracks,row);

    //--------------------------------------------------------
    //   Digitize this sector, row by row
    //   row[i] is the pointer to the TObjArray of TVectors,
    //   each one containing electrons accepted on this
    //   row, assigned into tracks
    //--------------------------------------------------------

    Int_t i;

    if (fDigitsArray->GetTree()==0) {
      AliFatal("Tree not set in fDigitsArray");
    }

    for (i=0;i<nrows;i++){
      
      AliDigits * dig = fDigitsArray->CreateRow(isec,i); 

      DigitizeRow(i,isec,row);

      fDigitsArray->StoreRow(isec,i);

      Int_t ndig = dig->GetDigitSize(); 
	
      AliDebug(10,
	       Form("*** Sector, row, compressed digits %d %d %d ***\n",
		    isec,i,ndig));        
 	
      fDigitsArray->ClearRow(isec,i);  

   
    } // end of the sector digitization

    for(i=0;i<nrows+2;i++){
      row[i]->Delete();  
      delete row[i];   
    }
      
    delete [] row; // delete the array of pointers to TObjArray-s
        

} // end of Hits2DigitsSector


//_____________________________________________________________________________
void AliTPC::DigitizeRow(Int_t irow,Int_t isec,TObjArray **rows)
{
  //-----------------------------------------------------------
  // Single row digitization, coupling from the neighbouring
  // rows taken into account
  //-----------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  // Modified: Marian Ivanov GSI Darmstadt, m.ivanov@gsi.de
  //-----------------------------------------------------------------
 
  Float_t zerosup = fTPCParam->GetZeroSup();
  AliTPCCalPad * gainTPC = AliTPCcalibDB::Instance()->GetDedxGainFactor(); 
  AliTPCCalPad * noiseTPC = AliTPCcalibDB::Instance()->GetPadNoise(); 
  AliTPCCalROC * gainROC = gainTPC->GetCalROC(isec);  // pad gains per given sector
  AliTPCCalROC * noiseROC = noiseTPC->GetCalROC(isec);  // noise per given sector


  fCurrentIndex[1]= isec;
  

  Int_t nofPads = fTPCParam->GetNPads(isec,irow);
  Int_t nofTbins = fTPCParam->GetMaxTBin();
  Int_t indexRange[4];
  //
  //  Integrated signal for this row
  //  and a single track signal
  //    

  TMatrixF *m1 = new TMatrixF(0,nofPads,0,nofTbins); // integrated
  TMatrixF *m2 = new TMatrixF(0,nofPads,0,nofTbins); // single
  //
  TMatrixF &total  = *m1;

  //  Array of pointers to the label-signal list

  Int_t nofDigits = nofPads*nofTbins; // number of digits for this row
  Float_t  **pList = new Float_t* [nofDigits]; 

  Int_t lp;
  Int_t i1;   
  for(lp=0;lp<nofDigits;lp++)pList[lp]=0; // set all pointers to NULL
  //
  //calculate signal 
  //
  Int_t row1=irow;
  Int_t row2=irow+2; 
  for (Int_t row= row1;row<=row2;row++){
    Int_t nTracks= rows[row]->GetEntries();
    for (i1=0;i1<nTracks;i1++){
      fCurrentIndex[2]= row;
      fCurrentIndex[3]=irow+1;
      if (row==irow+1){
	m2->Zero();  // clear single track signal matrix
	Float_t trackLabel = GetSignal(rows[row],i1,m2,m1,indexRange); 
	GetList(trackLabel,nofPads,m2,indexRange,pList);
      }
      else   GetSignal(rows[row],i1,0,m1,indexRange);
    }
  }
         
  Int_t tracks[3];

  AliDigits *dig = fDigitsArray->GetRow(isec,irow);
  Int_t gi=-1;
  Float_t fzerosup = zerosup+0.5;
  for(Int_t it=0;it<nofTbins;it++){
    for(Int_t ip=0;ip<nofPads;ip++){
      gi++;
      Float_t q=total(ip,it);      
      if(fDigitsSwitch == 0){	
	Float_t gain = gainROC->GetValue(irow,ip);  // get gain for given - pad-row pad	
	Float_t noisePad = noiseROC->GetValue(irow,ip);	
	//
	q*=gain;
	q+=GetNoise()*noisePad;
        if(q <=fzerosup) continue; // do not fill zeros
        q = TMath::Nint(q);
        if(q >= fTPCParam->GetADCSat()) q = fTPCParam->GetADCSat() - 1;  // saturation

      }

      else {
	if(q <= 0.) continue; // do not fill zeros
	if(q>2000.) q=2000.;
	q *= 16.;
	q = TMath::Nint(q);
      }

      //
      //  "real" signal or electronic noise (list = -1)?
      //    

      for(Int_t j1=0;j1<3;j1++){
	tracks[j1] = (pList[gi]) ?(Int_t)(*(pList[gi]+j1)) : -2;
      }

//Begin_Html
/*
  <A NAME="AliDigits"></A>
  using of AliDigits object
*/
//End_Html
      dig->SetDigitFast((Short_t)q,it,ip);
      if (fDigitsArray->IsSimulated()) {
	((AliSimDigits*)dig)->SetTrackIDFast(tracks[0],it,ip,0);
	((AliSimDigits*)dig)->SetTrackIDFast(tracks[1],it,ip,1);
	((AliSimDigits*)dig)->SetTrackIDFast(tracks[2],it,ip,2);
      }
    
    } // end of loop over time buckets
  }  // end of lop over pads 
  //
  // test
  //
  //

  // glitch filters if normal simulated digits
  //
  if(!fDigitsSwitch) ((AliSimDigits*)dig)->GlitchFilter();
  //
  //  This row has been digitized, delete nonused stuff
  //

  for(lp=0;lp<nofDigits;lp++){
    if(pList[lp]) delete [] pList[lp];
  }
  
  delete [] pList;

  delete m1;
  delete m2;

} // end of DigitizeRow

//_____________________________________________________________________________

Float_t AliTPC::GetSignal(TObjArray *p1, Int_t ntr, 
             TMatrixF *m1, TMatrixF *m2,Int_t *indexRange)
{

  //---------------------------------------------------------------
  //  Calculates 2-D signal (pad,time) for a single track,
  //  returns a pointer to the signal matrix and the track label 
  //  No digitization is performed at this level!!!
  //---------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  // Modified: Marian Ivanov 
  //-----------------------------------------------------------------

  TVector *tv;

  tv = (TVector*)p1->At(ntr); // pointer to a track
  TVector &v = *tv;
  
  Float_t label = v(0);
  Int_t centralPad = (fTPCParam->GetNPads(fCurrentIndex[1],fCurrentIndex[3]-1))/2;

  Int_t nElectrons = (tv->GetNrows()-1)/5;
  indexRange[0]=9999; // min pad
  indexRange[1]=-1; // max pad
  indexRange[2]=9999; //min time
  indexRange[3]=-1; // max time

  TMatrixF &signal = *m1;
  TMatrixF &total = *m2;
  //
  // Get LHC clock phase
  //
  TParameter<float> *ph;
  if(fDigitsSwitch){// s-digits
    ph = (TParameter<float>*)fLoader->TreeS()->GetUserInfo()->FindObject("lhcphase0");  
  }
  else{ // normal digits
    ph = (TParameter<float>*)fLoader->TreeD()->GetUserInfo()->FindObject("lhcphase0");
  } 
  //  Loop over all electrons
  //
  for(Int_t nel=0; nel<nElectrons; nel++){
    Int_t idx=nel*5;
    Float_t aval =  v(idx+4);
    Float_t eltoadcfac=aval*fTPCParam->GetTotalNormFac(); 
    Float_t xyz[4]={v(idx+1),v(idx+2),v(idx+3),v(idx+5)};
    Int_t n = ((AliTPCParamSR*)fTPCParam)->CalcResponseFast(xyz,fCurrentIndex,
							    fCurrentIndex[3],ph->GetVal());

    Int_t *index = fTPCParam->GetResBin(0);  
    Float_t *weight = & (fTPCParam->GetResWeight(0));

    if (n>0) for (Int_t i =0; i<n; i++){       
      Int_t pad=index[1]+centralPad;  //in digit coordinates central pad has coordinate 0

      if (pad>=0){
	Int_t time=index[2];	 
	Float_t qweight = *(weight)*eltoadcfac;
	
	if (m1!=0) signal(pad,time)+=qweight;
	total(pad,time)+=qweight;
	if (indexRange[0]>pad) indexRange[0]=pad;
	if (indexRange[1]<pad) indexRange[1]=pad;
	if (indexRange[2]>time) indexRange[2]=time;
	if (indexRange[3]<time) indexRange[3]=time;
	
	index+=3;
	weight++;	

      }	 
    }
  } // end of loop over electrons
  
  return label; // returns track label when finished
}

//_____________________________________________________________________________
void AliTPC::GetList(Float_t label,Int_t np,TMatrixF *m,
                     Int_t *indexRange, Float_t **pList)
{
  //----------------------------------------------------------------------
  //  Updates the list of tracks contributing to digits for a given row
  //----------------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  TMatrixF &signal = *m;

  // lop over nonzero digits

  for(Int_t it=indexRange[2];it<indexRange[3]+1;it++){
    for(Int_t ip=indexRange[0];ip<indexRange[1]+1;ip++){


      // accept only the contribution larger than 500 electrons (1/2 s_noise)

      if(signal(ip,it)<0.5) continue; 

      Int_t globalIndex = it*np+ip; // globalIndex starts from 0!
        
      if(!pList[globalIndex]){
        
	// 
	// Create new list (6 elements - 3 signals and 3 labels),
	//

	pList[globalIndex] = new Float_t [6];

	// set list to -1 
	
	*pList[globalIndex] = -1.;
	*(pList[globalIndex]+1) = -1.;
	*(pList[globalIndex]+2) = -1.;
	*(pList[globalIndex]+3) = -1.;
	*(pList[globalIndex]+4) = -1.;
	*(pList[globalIndex]+5) = -1.;

	*pList[globalIndex] = label;
	*(pList[globalIndex]+3) = signal(ip,it);
      }
      else {

	// check the signal magnitude

	Float_t highest = *(pList[globalIndex]+3);
	Float_t middle = *(pList[globalIndex]+4);
	Float_t lowest = *(pList[globalIndex]+5);
	
	//
	//  compare the new signal with already existing list
	//
	
	if(signal(ip,it)<lowest) continue; // neglect this track

	//

	if (signal(ip,it)>highest){
	  *(pList[globalIndex]+5) = middle;
	  *(pList[globalIndex]+4) = highest;
	  *(pList[globalIndex]+3) = signal(ip,it);
	  
	  *(pList[globalIndex]+2) = *(pList[globalIndex]+1);
	  *(pList[globalIndex]+1) = *pList[globalIndex];
	  *pList[globalIndex] = label;
	}
	else if (signal(ip,it)>middle){
	  *(pList[globalIndex]+5) = middle;
	  *(pList[globalIndex]+4) = signal(ip,it);
	  
	  *(pList[globalIndex]+2) = *(pList[globalIndex]+1);
	  *(pList[globalIndex]+1) = label;
	}
	else{
	  *(pList[globalIndex]+5) = signal(ip,it);
	  *(pList[globalIndex]+2) = label;
	}
      }
      
    } // end of loop over pads
  } // end of loop over time bins

}//end of GetList
//___________________________________________________________________
void AliTPC::MakeSector(Int_t isec,Int_t nrows,TTree *TH,
                        Stat_t ntracks,TObjArray **row)
{

  //-----------------------------------------------------------------
  // Prepares the sector digitization, creates the vectors of
  // tracks for each row of this sector. The track vector
  // contains the track label and the position of electrons.
  //-----------------------------------------------------------------

  // 
  // The trasport of the electrons through TPC drift volume
  //    Drift (drift velocity + velocity map(not yet implemented)))
  //    Application of the random processes (diffusion, gas gain)
  //    Systematic effects (ExB effect in drift volume + ROCs)  
  //
  // Algorithm:
  // Loop over primary electrons:
  //    Creation of the secondary electrons
  //    Loop over electrons (primary+ secondaries)
  //        Global coordinate frame:
  //          1. Skip electrons if attached  
  //          2. ExB effect in drift volume
  //             a.) Simulation   calib->GetExB()->CorrectInverse(dxyz0,dxyz1);
  //             b.) Reconstruction -  calib->GetExB()->CorrectInverse(dxyz0,dxyz1);
  //          3. Generation of gas gain (Random - Exponential distribution) 
  //          4. TransportElectron function (diffusion)
  //
  //        5. Conversion to the local coordinate frame  pad-row, pad, timebin
  //        6. Apply Time0 shift - AliTPCCalPad class 
  //            a.) Plus sign in simulation
  //            b.) Minus sign in reconstruction 
  // end of loop          
  //
  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  // Origin: Marian Ivanov,  marian.ivanov@cern.ch
  //-----------------------------------------------------------------
  AliTPCcalibDB* const calib=AliTPCcalibDB::Instance();
  if (gAlice){ // Set correctly the magnetic field in the ExB calculation
    if (!calib->GetExB()){
      AliMagF * field = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField()); 
      if (field) {
	calib->SetExBField(field);
      }
    }
  }

  Float_t gasgain = fTPCParam->GetGasGain();
  gasgain = gasgain/fGainFactor;
  Int_t i;
  Float_t xyz[5]; 

  AliTPChit *tpcHit; // pointer to a sigle TPC hit    
  //MI change
  TBranch * branch=0;
  if (fHitType>1) branch = TH->GetBranch("TPC2");
  else branch = TH->GetBranch("TPC");

 
  //----------------------------------------------
  // Create TObjArray-s, one for each row,
  // each TObjArray will store the TVectors
  // of electrons, one TVectors per each track.
  //---------------------------------------------- 
    
  Int_t *nofElectrons = new Int_t [nrows+2]; // electron counter for each row
  TVector **tracks = new TVector* [nrows+2]; //pointers to the track vectors

  for(i=0; i<nrows+2; i++){
    row[i] = new TObjArray;
    nofElectrons[i]=0;
    tracks[i]=0;
  }

 

  //--------------------------------------------------------------------
  //  Loop over tracks, the "track" contains the full history
  //--------------------------------------------------------------------
  
  Int_t previousTrack,currentTrack;
  previousTrack = -1; // nothing to store so far!

  for(Int_t track=0;track<ntracks;track++){
    Bool_t isInSector=kTRUE;
    ResetHits();
    isInSector = TrackInVolume(isec,track);
    if (!isInSector) continue;
    //MI change
    branch->GetEntry(track); // get next track
    
    //M.I. changes

    tpcHit = (AliTPChit*)FirstHit(-1);

    //--------------------------------------------------------------
    //  Loop over hits
    //--------------------------------------------------------------


    while(tpcHit){
      
      Int_t sector=tpcHit->fSector; // sector number
      if(sector != isec){
	tpcHit = (AliTPChit*) NextHit();
	continue; 
      }

      // Remove hits which arrive before the TPC opening gate signal
      if(((fTPCParam->GetZLength(isec)-TMath::Abs(tpcHit->Z()))
	  /fTPCParam->GetDriftV()+tpcHit->Time())<fTPCParam->GetGateDelay()) {
	tpcHit = (AliTPChit*) NextHit();
	continue;
      }

      currentTrack = tpcHit->Track(); // track number

      if(currentTrack != previousTrack){
                          
	// store already filled fTrack
              
	for(i=0;i<nrows+2;i++){
	  if(previousTrack != -1){
	    if(nofElectrons[i]>0){
	      TVector &v = *tracks[i];
	      v(0) = previousTrack;
	      tracks[i]->ResizeTo(5*nofElectrons[i]+1); // shrink if necessary
	      row[i]->Add(tracks[i]);                     
	    }
	    else {
	      delete tracks[i]; // delete empty TVector
	      tracks[i]=0;
	    }
	  }

	  nofElectrons[i]=0;
	  tracks[i] = new TVector(601); // TVectors for the next fTrack

	} // end of loop over rows
	       
	previousTrack=currentTrack; // update track label 
      }
	   
      Int_t qI = (Int_t) (tpcHit->fQ); // energy loss (number of electrons)

      //---------------------------------------------------
      //  Calculate the electron attachment probability
      //---------------------------------------------------


      Float_t time = 1.e6*(fTPCParam->GetZLength(isec)-TMath::Abs(tpcHit->Z()))
	/fTPCParam->GetDriftV(); 
      // in microseconds!	
      Float_t attProb = fTPCParam->GetAttCoef()*
	fTPCParam->GetOxyCont()*time; //  fraction! 
   
      //-----------------------------------------------
      //  Loop over electrons
      //-----------------------------------------------
      Int_t index[3];
      index[1]=isec;
      for(Int_t nel=0;nel<qI;nel++){
	// skip if electron lost due to the attachment
	if((gRandom->Rndm(0)) < attProb) continue; // electron lost!
	
	//
	// ExB effect
	//
	Double_t dxyz0[3],dxyz1[3];
	dxyz0[0]=tpcHit->X();
	dxyz0[1]=tpcHit->Y();
	dxyz0[2]=tpcHit->Z(); 	
	if (calib->GetExB()){
	  calib->GetExB()->CorrectInverse(dxyz0,dxyz1);
	}else{
	  AliError("Not valid ExB calibration");
	  dxyz1[0]=tpcHit->X();
	  dxyz1[1]=tpcHit->Y();
	  dxyz1[2]=tpcHit->Z(); 	
	}
	xyz[0]=dxyz1[0];
	xyz[1]=dxyz1[1];
	xyz[2]=dxyz1[2]; 	
	//
	//
	//
	// protection for the nonphysical avalanche size (10**6 maximum)
	//  
	Double_t rn=TMath::Max(gRandom->Rndm(0),1.93e-22);
	xyz[3]= (Float_t) (-gasgain*TMath::Log(rn)); 
	index[0]=1;
	  
	TransportElectron(xyz,index);    
	Int_t rowNumber;
	Int_t padrow = fTPCParam->GetPadRow(xyz,index); 
	//
	// Add Time0 correction due unisochronity
	// xyz[0] - pad row coordinate 
	// xyz[1] - pad coordinate
	// xyz[2] - is in now time bin coordinate system
	Float_t correction =0;
	if (calib->GetPadTime0()){
	  if (!calib->GetPadTime0()->GetCalROC(isec)) continue;	  
	  Int_t npads = fTPCParam->GetNPads(isec,padrow);
	  //	  Int_t pad  = TMath::Nint(xyz[1]+fTPCParam->GetNPads(isec,TMath::Nint(xyz[0]))*0.5);
	  // pad numbering from -npads/2 .. npads/2-1
	  Int_t pad  = TMath::Nint(xyz[1]+npads/2);
	  if (pad<0) pad=0;
	  if (pad>=npads) pad=npads-1;
	  correction = calib->GetPadTime0()->GetCalROC(isec)->GetValue(padrow,pad);
	  //	  printf("%d\t%d\t%d\t%f\n",isec,padrow,pad,correction);
	  if (fDebugStreamer){
	    (*fDebugStreamer)<<"Time0"<<
	      "isec="<<isec<<
	      "padrow="<<padrow<<
	      "pad="<<pad<<
	      "x0="<<xyz[0]<<
	      "x1="<<xyz[1]<<
	      "x2="<<xyz[2]<<
	      "hit.="<<tpcHit<<
	      "cor="<<correction<<
	      "\n";
	  }
	}
	xyz[2]+=correction;
	xyz[2]+=fTPCParam->GetNTBinsL1();    // adding Level 1 time bin offset
	//
	// Electron track time (for pileup simulation)
	xyz[2]+=tpcHit->Time()/fTPCParam->GetTSample(); // adding time of flight
	xyz[4] =0;

	//
	// row 0 - cross talk from the innermost row
	// row fNRow+1 cross talk from the outermost row
	rowNumber = index[2]+1; 
	//transform position to local digit coordinates
	//relative to nearest pad row 
	if ((rowNumber<0)||rowNumber>fTPCParam->GetNRow(isec)+1) continue;
	/*	Float_t x1,y1;
	if (isec <fTPCParam->GetNInnerSector()) {
	  x1 = xyz[1]*fTPCParam->GetInnerPadPitchWidth();
	  y1 = fTPCParam->GetYInner(rowNumber);
	}
	else{
	  x1=xyz[1]*fTPCParam->GetOuterPadPitchWidth();
	  y1 = fTPCParam->GetYOuter(rowNumber);
	}
	// gain inefficiency at the wires edges - linear
	x1=TMath::Abs(x1);
	y1-=1.;
	if(x1>y1) xyz[3]*=TMath::Max(1.e-6,(y1-x1+1.));	*/
	
	nofElectrons[rowNumber]++;	  
	//----------------------------------
	// Expand vector if necessary
	//----------------------------------
	if(nofElectrons[rowNumber]>120){
	  Int_t range = tracks[rowNumber]->GetNrows();
	  if((nofElectrons[rowNumber])>(range-1)/5){
	    
	    tracks[rowNumber]->ResizeTo(range+500); // Add 100 electrons
	  }
	}
	
	TVector &v = *tracks[rowNumber];
	Int_t idx = 5*nofElectrons[rowNumber]-4;
	Real_t * position = &(((TVector&)v)(idx)); //make code faster
	memcpy(position,xyz,5*sizeof(Float_t));
	
      } // end of loop over electrons

      tpcHit = (AliTPChit*)NextHit();
      
    } // end of loop over hits
  } // end of loop over tracks

    //
    //   store remaining track (the last one) if not empty
    //
  
  for(i=0;i<nrows+2;i++){
    if(nofElectrons[i]>0){
      TVector &v = *tracks[i];
      v(0) = previousTrack;
      tracks[i]->ResizeTo(5*nofElectrons[i]+1); // shrink if necessary
      row[i]->Add(tracks[i]);  
    }
    else{
      delete tracks[i];
      tracks[i]=0;
    }  
  }  
  
  delete [] tracks;
  delete [] nofElectrons;

} // end of MakeSector


//_____________________________________________________________________________
void AliTPC::Init()
{
  //
  // Initialise TPC detector after definition of geometry
  //
  AliDebug(1,"*********************************************");
}

//_____________________________________________________________________________
void AliTPC::ResetDigits()
{
  //
  // Reset number of digits and the digits array for this detector
  //
  fNdigits   = 0;
  if (fDigits)   fDigits->Clear();
}



//_____________________________________________________________________________
void AliTPC::SetSens(Int_t sens)
{

  //-------------------------------------------------------------
  // Activates/deactivates the sensitive strips at the center of
  // the pad row -- this is for the space-point resolution calculations
  //-------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  fSens = sens;
}

 
void AliTPC::SetSide(Float_t side=0.)
{
  // choice of the TPC side

  fSide = side;
 
}
//_____________________________________________________________________________

void AliTPC::TransportElectron(Float_t *xyz, Int_t *index)
{
  //
  // electron transport taking into account:
  // 1. diffusion, 
  // 2.ExB at the wires
  // 3. nonisochronity
  //
  // xyz and index must be already transformed to system 1
  //

  fTPCParam->Transform1to2(xyz,index);  // mis-alignment applied in this step
  
  //add diffusion
  Float_t driftl=xyz[2];
  if(driftl<0.01) driftl=0.01;
  driftl=TMath::Sqrt(driftl);
  Float_t sigT = driftl*(fTPCParam->GetDiffT());
  Float_t sigL = driftl*(fTPCParam->GetDiffL());
  xyz[0]=gRandom->Gaus(xyz[0],sigT);
  xyz[1]=gRandom->Gaus(xyz[1],sigT);
  xyz[2]=gRandom->Gaus(xyz[2],sigL);

  // ExB
  
  if (fTPCParam->GetMWPCReadout()==kTRUE){
    Float_t dx = fTPCParam->Transform2to2NearestWire(xyz,index);
    xyz[1]+=dx*(fTPCParam->GetOmegaTau());
  }
  //add nonisochronity (not implemented yet) 
 
  
}
  
ClassImp(AliTPChit)
  //______________________________________________________________________
  AliTPChit::AliTPChit()
            :AliHit(),
	     fSector(0),
	     fPadRow(0),
	     fQ(0),
	     fTime(0)
{
  //
  // default
  //

}
//_____________________________________________________________________________
AliTPChit::AliTPChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits)
          :AliHit(shunt,track),
	     fSector(0),
	     fPadRow(0),
	     fQ(0),
	     fTime(0)
{
  //
  // Creates a TPC hit object
  //
  fSector     = vol[0];
  fPadRow     = vol[1];
  fX          = hits[0];
  fY          = hits[1];
  fZ          = hits[2];
  fQ          = hits[3];
  fTime       = hits[4];
}
 
//________________________________________________________________________
// Additional code because of the AliTPCTrackHitsV2

void AliTPC::MakeBranch(Option_t *option)
{
  //
  // Create a new branch in the current Root Tree
  // The branch of fHits is automatically split
  // MI change 14.09.2000
  AliDebug(1,"");
  if (fHitType<2) return;
  char branchname[10];
  //sprintf(branchname,"%s2",GetName()); 
  snprintf(branchname,10,"%s2",GetName()); 
  //
  // Get the pointer to the header
  const char *cH = strstr(option,"H");
  //
  if (fTrackHits   && fLoader->TreeH() && cH && fHitType&4) {
    AliDebug(1,"Making branch for Type 4 Hits");
    fLoader->TreeH()->Branch(branchname,"AliTPCTrackHitsV2",&fTrackHits,fBufferSize,99);
  }

//   if (fTrackHitsOld   && fLoader->TreeH() && cH && fHitType&2) {    
//     AliDebug(1,"Making branch for Type 2 Hits");
//     AliObjectBranch * branch = new AliObjectBranch(branchname,"AliTPCTrackHits",&fTrackHitsOld, 
//                                                    fLoader->TreeH(),fBufferSize,99);
//     fLoader->TreeH()->GetListOfBranches()->Add(branch);
//   }	
}

void AliTPC::SetTreeAddress()
{
  //Sets tree address for hits  
  if (fHitType<=1) {
    if (fHits == 0x0 ) fHits = new TClonesArray("AliTPChit", 176);//skowron 20.06.03
    AliDetector::SetTreeAddress();
  }
  if (fHitType>1) SetTreeAddress2();
}

void AliTPC::SetTreeAddress2()
{
  //
  // Set branch address for the TrackHits Tree
  // 
  AliDebug(1,"");
  
  TBranch *branch;
  char branchname[20];
  //sprintf(branchname,"%s2",GetName());
  snprintf(branchname,20,"%s2",GetName());
  //
  // Branch address for hit tree
  TTree *treeH = fLoader->TreeH();
  if ((treeH)&&(fHitType&4)) {
    branch = treeH->GetBranch(branchname);
    if (branch) {
      branch->SetAddress(&fTrackHits);
      AliDebug(1,"fHitType&4 Setting");
    }
    else 
      AliDebug(1,"fHitType&4 Failed (can not find branch)");
    
  }
 //  if ((treeH)&&(fHitType&2)) {
//     branch = treeH->GetBranch(branchname);
//     if (branch) {
//       branch->SetAddress(&fTrackHitsOld);
//       AliDebug(1,"fHitType&2 Setting");
//     }
//     else
//       AliDebug(1,"fHitType&2 Failed (can not find branch)");
//   }
}

void AliTPC::FinishPrimary()
{
  if (fTrackHits &&fHitType&4)      fTrackHits->FlushHitStack();  
  //  if (fTrackHitsOld && fHitType&2)  fTrackHitsOld->FlushHitStack();  
}


void AliTPC::AddHit2(Int_t track, Int_t *vol, Float_t *hits)
{ 
  //
  // add hit to the list

  Int_t rtrack;
  if (fIshunt) {
    int primary = gAlice->GetMCApp()->GetPrimary(track);
    gAlice->GetMCApp()->Particle(primary)->SetBit(kKeepBit);
    rtrack=primary;
  } else {
    rtrack=track;
    gAlice->GetMCApp()->FlagTrack(track);
  }  
  if (fTrackHits && fHitType&4) 
    fTrackHits->AddHitKartez(vol[0],rtrack, hits[0],
                             hits[1],hits[2],(Int_t)hits[3],hits[4]);
 //  if (fTrackHitsOld &&fHitType&2 ) 
//     fTrackHitsOld->AddHitKartez(vol[0],rtrack, hits[0],
//                                 hits[1],hits[2],(Int_t)hits[3]);
  
}

void AliTPC::ResetHits()
{ 
  if (fHitType&1) AliDetector::ResetHits();
  if (fHitType>1) ResetHits2();
}

void AliTPC::ResetHits2()
{
  //
  //reset hits
  if (fTrackHits && fHitType&4) fTrackHits->Clear();
  // if (fTrackHitsOld && fHitType&2) fTrackHitsOld->Clear();

}   

AliHit* AliTPC::FirstHit(Int_t track)
{
  if (fHitType>1) return FirstHit2(track);
  return AliDetector::FirstHit(track);
}
AliHit* AliTPC::NextHit()
{
  //
  // gets next hit
  //
  if (fHitType>1) return NextHit2();
  
  return AliDetector::NextHit();
}

AliHit* AliTPC::FirstHit2(Int_t track)
{
  //
  // Initialise the hit iterator
  // Return the address of the first hit for track
  // If track>=0 the track is read from disk
  // while if track<0 the first hit of the current
  // track is returned
  // 
  if(track>=0) {
    gAlice->GetMCApp()->ResetHits();
    fLoader->TreeH()->GetEvent(track);
  }
  //
  if (fTrackHits && fHitType&4) {
    fTrackHits->First();
    return fTrackHits->GetHit();
  }
 //  if (fTrackHitsOld && fHitType&2) {
//     fTrackHitsOld->First();
//     return fTrackHitsOld->GetHit();
//   }

  else return 0;
}

AliHit* AliTPC::NextHit2()
{
  //
  //Return the next hit for the current track


//   if (fTrackHitsOld && fHitType&2) {
//     fTrackHitsOld->Next();
//     return fTrackHitsOld->GetHit();
//   }
  if (fTrackHits) {
    fTrackHits->Next();
    return fTrackHits->GetHit();
  }
  else 
    return 0;
}

void AliTPC::RemapTrackHitIDs(Int_t *map)
{
  //
  // remapping
  //
  if (!fTrackHits) return;
  
//   if (fTrackHitsOld && fHitType&2){
//     AliObjectArray * arr = fTrackHitsOld->fTrackHitsInfo;
//     for (UInt_t i=0;i<arr->GetSize();i++){
//       AliTrackHitsInfo * info = (AliTrackHitsInfo *)(arr->At(i));
//       info->fTrackID = map[info->fTrackID];
//     }
//   }
//  if (fTrackHitsOld && fHitType&4){
  if (fTrackHits && fHitType&4){
    TClonesArray * arr = fTrackHits->GetArray();;
    for (Int_t i=0;i<arr->GetEntriesFast();i++){
      AliTrackHitsParamV2 * info = (AliTrackHitsParamV2 *)(arr->At(i));
      info->SetTrackID(map[info->GetTrackID()]);
    }
  }
}

Bool_t   AliTPC::TrackInVolume(Int_t id,Int_t track)
{
  //return bool information - is track in given volume
  //load only part of the track information 
  //return true if current track is in volume
  //
  //  return kTRUE;
 //  if (fTrackHitsOld && fHitType&2) {
//     TBranch * br = fLoader->TreeH()->GetBranch("fTrackHitsInfo");
//     br->GetEvent(track);
//     AliObjectArray * ar = fTrackHitsOld->fTrackHitsInfo;
//     for (UInt_t j=0;j<ar->GetSize();j++){
//       if (  ((AliTrackHitsInfo*)ar->At(j))->fVolumeID==id) return kTRUE;
//     } 
//   }

  if (fTrackHits && fHitType&4) {
    TBranch * br1 = fLoader->TreeH()->GetBranch("fVolumes");
    TBranch * br2 = fLoader->TreeH()->GetBranch("fNVolumes");    
    br2->GetEvent(track);
    br1->GetEvent(track);    
    Int_t *volumes = fTrackHits->GetVolumes();
    Int_t nvolumes = fTrackHits->GetNVolumes();
    if (!volumes && nvolumes>0) {
      AliWarning(Form("Problematic track\t%d\t%d",track,nvolumes));
      return kFALSE;
    }
    for (Int_t j=0;j<nvolumes; j++)
      if (volumes[j]==id) return kTRUE;    
  }

  if (fHitType&1) {
    TBranch * br = fLoader->TreeH()->GetBranch("fSector");
    br->GetEvent(track);
    for (Int_t j=0;j<fHits->GetEntriesFast();j++){
      if (  ((AliTPChit*)fHits->At(j))->fSector==id) return kTRUE;
    } 
  }
  return kFALSE;  

}


AliLoader* AliTPC::MakeLoader(const char* topfoldername)
{
  //Makes TPC loader
  fLoader = new AliTPCLoader(GetName(),topfoldername);
  return fLoader;
}

////////////////////////////////////////////////////////////////////////
AliTPCParam* AliTPC::LoadTPCParam(TFile *file) {
//
// load TPC paarmeters from a given file or create new if the object
// is not found there
// 12/05/2003 This method should be moved to the AliTPCLoader
// and one has to decide where to store the TPC parameters
// M.Kowalski
  char paramName[50];
  //sprintf(paramName,"75x40_100x60_150x60");
  snprintf(paramName,50,"75x40_100x60_150x60");
  AliTPCParam *paramTPC=(AliTPCParam*)file->Get(paramName);
  if (paramTPC) {
    AliDebugClass(1,Form("TPC parameters %s found.",paramName));
  } else {
    AliWarningClass("TPC parameters not found. Create new (they may be incorrect)");
    //paramTPC = new AliTPCParamSR;
    paramTPC = AliTPCcalibDB::Instance()->GetParameters();
    if (!paramTPC->IsGeoRead()){
      //
      // read transformation matrices for gGeoManager
      //
      paramTPC->ReadGeoMatrices();
    }
  
  }
  return paramTPC;

// the older version of parameters can be accessed with this code.
// In some cases, we have old parameters saved in the file but 
// digits were created with new parameters, it can be distinguish 
// by the name of TPC TreeD. The code here is just for the case 
// we would need to compare with old data, uncomment it if needed.
//
//  char paramName[50];
//  sprintf(paramName,"75x40_100x60");
//  AliTPCParam *paramTPC=(AliTPCParam*)in->Get(paramName);
//  if (paramTPC) {
//    cout<<"TPC parameters "<<paramName<<" found."<<endl;
//  } else {
//    sprintf(paramName,"75x40_100x60_150x60");
//    paramTPC=(AliTPCParam*)in->Get(paramName);
//    if (paramTPC) {
//	cout<<"TPC parameters "<<paramName<<" found."<<endl;
//    } else {
//	cerr<<"TPC parameters not found. Create new (they may be incorrect)."
//	    <<endl;    
//	paramTPC = new AliTPCParamSR;
//    }
//  }
//  return paramTPC;

}


