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

#include <TFile.h>  
#include <TGeometry.h>
#include <TInterpreter.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TNode.h>
#include <TObjectTable.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TSystem.h>     
#include <TTUBS.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <TString.h>
#include <TF2.h>
#include <TStopwatch.h>

#include "AliArrayBranch.h"
#include "AliDigits.h"
#include "AliMagF.h"
#include "AliPoints.h"
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
#include "AliTPCTrackHits.h"
#include "AliTPCTrackHitsV2.h"
#include "AliTrackReference.h"
#include "AliMC.h"
#include "AliTPCDigitizer.h"
#include "AliTPCBuffer.h"
#include "AliTPCDDLRawData.h"


ClassImp(AliTPC) 

//_____________________________________________________________________________
// helper class for fast matrix and vector manipulation - no range checking
// origin - Marian Ivanov

class AliTPCFastMatrix : public TMatrix {
public :
  AliTPCFastMatrix(Int_t rowlwb, Int_t rowupb, Int_t collwb, Int_t colupb);
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  Float_t & UncheckedAt(Int_t rown, Int_t coln) const  {return fElements[(rown-fRowLwb)*fNcols+(coln-fColLwb)];} //fast acces
  Float_t   UncheckedAtFast(Int_t rown, Int_t coln) const  {return fElements[(rown-fRowLwb)*fNcols+(coln-fColLwb)];} //fast acces
#else
  Float_t & UncheckedAt(Int_t rown, Int_t coln) const  {return  (fIndex[coln])[rown];} //fast acces   
  Float_t   UncheckedAtFast(Int_t rown, Int_t coln) const  {return  (fIndex[coln])[rown];} //fast acces   
#endif
};

AliTPCFastMatrix::AliTPCFastMatrix(Int_t rowlwb, Int_t rowupb, Int_t collwb, Int_t colupb):
  TMatrix(rowlwb, rowupb,collwb,colupb)
   {
   };
//_____________________________________________________________________________
AliTPC::AliTPC()
{
  //
  // Default constructor
  //
  fIshunt   = 0;
  fHits     = 0;
  fDigits   = 0;
  fNsectors = 0;
  fDigitsArray = 0;
  fDefaults = 0;
  fTrackHits = 0; 
  fTrackHitsOld = 0;   
  fHitType = 2; //default CONTAINERS - based on ROOT structure 
  fTPCParam = 0;    
  fNoiseTable = 0;
  fActiveSectors =0;

}
 
//_____________________________________________________________________________
AliTPC::AliTPC(const char *name, const char *title)
      : AliDetector(name,title)
{
  //
  // Standard constructor
  //

  //
  // Initialise arrays of hits and digits 
  fHits     = new TClonesArray("AliTPChit",  176);
  gAlice->GetMCApp()->AddHitList(fHits); 
  fDigitsArray = 0;
  fDefaults = 0;
  //
  fTrackHits = new AliTPCTrackHitsV2;  
  fTrackHits->SetHitPrecision(0.002);
  fTrackHits->SetStepPrecision(0.003);  
  fTrackHits->SetMaxDistance(100);

  fTrackHitsOld = new AliTPCTrackHits;  //MI - 13.09.2000
  fTrackHitsOld->SetHitPrecision(0.002);
  fTrackHitsOld->SetStepPrecision(0.003);  
  fTrackHitsOld->SetMaxDistance(100); 

  fNoiseTable =0;

  fHitType = 2;
  fActiveSectors = 0;
  //
  // Initialise counters
  fNsectors = 0;

  //
  fIshunt     =  0;
  //
  // Initialise color attributes
  SetMarkerColor(kYellow);

  //
  //  Set TPC parameters
  //


  if (!strcmp(title,"Default")) {       
    fTPCParam = new AliTPCParamSR;
  } else {
    cerr<<"AliTPC warning: in Config.C you must set non-default parameters\n";
    fTPCParam=0;
  }

}

//_____________________________________________________________________________
AliTPC::AliTPC(const AliTPC& t):AliDetector(t){
  //
  // dummy copy constructor
  //
}
AliTPC::~AliTPC()
{
  //
  // TPC destructor
  //

  fIshunt   = 0;
  delete fHits;
  delete fDigits;
  delete fTPCParam;
  delete fTrackHits; //MI 15.09.2000
  delete fTrackHitsOld; //MI 10.12.2001
  if (fNoiseTable) delete [] fNoiseTable;

}

//_____________________________________________________________________________
void AliTPC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a hit to the list
  //
  //  TClonesArray &lhits = *fHits;
  //  new(lhits[fNhits++]) AliTPChit(fIshunt,track,vol,hits);
  if (fHitType&1){
    TClonesArray &lhits = *fHits;
    new(lhits[fNhits++]) AliTPChit(fIshunt,track,vol,hits);
  }
  if (fHitType>1)
   AddHit2(track,vol,hits);
}

//_____________________________________________________________________________
void AliTPC::BuildGeometry()
{

  //
  // Build TPC ROOT TNode geometry for the event display
  //
  TNode *nNode, *nTop;
  TTUBS *tubs;
  Int_t i;
  const int kColorTPC=19;
  char name[5], title[25];
  const Double_t kDegrad=TMath::Pi()/180;
  const Double_t kRaddeg=180./TMath::Pi();


  Float_t innerOpenAngle = fTPCParam->GetInnerAngle();
  Float_t outerOpenAngle = fTPCParam->GetOuterAngle();

  Float_t innerAngleShift = fTPCParam->GetInnerAngleShift();
  Float_t outerAngleShift = fTPCParam->GetOuterAngleShift();

  Int_t nLo = fTPCParam->GetNInnerSector()/2;
  Int_t nHi = fTPCParam->GetNOuterSector()/2;  

  const Double_t kloAng = (Double_t)TMath::Nint(innerOpenAngle*kRaddeg);
  const Double_t khiAng = (Double_t)TMath::Nint(outerOpenAngle*kRaddeg);
  const Double_t kloAngSh = (Double_t)TMath::Nint(innerAngleShift*kRaddeg);
  const Double_t khiAngSh = (Double_t)TMath::Nint(outerAngleShift*kRaddeg);  


  const Double_t kloCorr = 1/TMath::Cos(0.5*kloAng*kDegrad);
  const Double_t khiCorr = 1/TMath::Cos(0.5*khiAng*kDegrad);

  Double_t rl,ru;
  

  //
  // Get ALICE top node
  //

  nTop=gAlice->GetGeometry()->GetNode("alice");

  //  inner sectors

  rl = fTPCParam->GetInnerRadiusLow();
  ru = fTPCParam->GetInnerRadiusUp();
 

  for(i=0;i<nLo;i++) {
    sprintf(name,"LS%2.2d",i);
    name[4]='\0';
    sprintf(title,"TPC low sector %3d",i);
    title[24]='\0';
    
    tubs = new TTUBS(name,title,"void",rl*kloCorr,ru*kloCorr,250.,
                     kloAng*(i-0.5)+kloAngSh,kloAng*(i+0.5)+kloAngSh);
    tubs->SetNumberOfDivisions(1);
    nTop->cd();
    nNode = new TNode(name,title,name,0,0,0,"");
    nNode->SetLineColor(kColorTPC);
    fNodes->Add(nNode);
  }

  // Outer sectors

  rl = fTPCParam->GetOuterRadiusLow();
  ru = fTPCParam->GetOuterRadiusUp();

  for(i=0;i<nHi;i++) {
    sprintf(name,"US%2.2d",i);
    name[4]='\0';
    sprintf(title,"TPC upper sector %d",i);
    title[24]='\0';
    tubs = new TTUBS(name,title,"void",rl*khiCorr,ru*khiCorr,250,
                     khiAng*(i-0.5)+khiAngSh,khiAng*(i+0.5)+khiAngSh);
    tubs->SetNumberOfDivisions(1);
    nTop->cd();
    nNode = new TNode(name,title,name,0,0,0,"");
    nNode->SetLineColor(kColorTPC);
    fNodes->Add(nNode);
  }

}    

//_____________________________________________________________________________
Int_t AliTPC::DistancetoPrimitive(Int_t , Int_t ) const
{
  //
  // Calculate distance from TPC to mouse on the display
  // Dummy procedure
  //
  return 9999;
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

  Int_t iSXFLD=gAlice->Field()->Integ();
  Float_t sXMGMX=gAlice->Field()->Max();

  Float_t amat[5]; // atomic numbers
  Float_t zmat[5]; // z
  Float_t wmat[5]; // proportions

  Float_t density;
  Float_t apure[2];


  //***************** Gases *************************
  
  //-------------------------------------------------
  // pure gases
  //-------------------------------------------------

  // Neon


  amat[0]= 20.18;
  zmat[0]= 10.;  
  density = 0.0009;
 
  apure[0]=amat[0];

  AliMaterial(20,"Ne",amat[0],zmat[0],density,999.,999.);

  // Argon

  amat[0]= 39.948;
  zmat[0]= 18.;  
  density = 0.001782;  

  apure[1]=amat[0];

  AliMaterial(21,"Ar",amat[0],zmat[0],density,999.,999.);
 

  //--------------------------------------------------------------
  // gases - compounds
  //--------------------------------------------------------------

  Float_t amol[3];

  // CO2

  amat[0]=12.011;
  amat[1]=15.9994;

  zmat[0]=6.;
  zmat[1]=8.;

  wmat[0]=1.;
  wmat[1]=2.;

  density=0.001977;

  amol[0] = amat[0]*wmat[0]+amat[1]*wmat[1];

  AliMixture(10,"CO2",amat,zmat,density,-2,wmat);
  
  // CF4

  amat[0]=12.011;
  amat[1]=18.998;

  zmat[0]=6.;
  zmat[1]=9.;
 
  wmat[0]=1.;
  wmat[1]=4.;
 
  density=0.003034;

  amol[1] = amat[0]*wmat[0]+amat[1]*wmat[1];

  AliMixture(11,"CF4",amat,zmat,density,-2,wmat); 


  // CH4

  amat[0]=12.011;
  amat[1]=1.;

  zmat[0]=6.;
  zmat[1]=1.;

  wmat[0]=1.;
  wmat[1]=4.;

  density=0.000717;

  amol[2] = amat[0]*wmat[0]+amat[1]*wmat[1];

  AliMixture(12,"CH4",amat,zmat,density,-2,wmat);

  //----------------------------------------------------------------
  // gases - mixtures, ID >= 20 pure gases, <= 10 ID < 20 -compounds
  //----------------------------------------------------------------

  char namate[21]=""; 
  density = 0.;
  Float_t am=0;
  Int_t nc;
  Float_t rho,absl,x0,buf[1];
  Int_t nbuf;
  Float_t a,z;

  for(nc = 0;nc<fNoComp;nc++)
    {
    
      // retrive material constants
      
      gMC->Gfmate((*fIdmate)[fMixtComp[nc]],namate,a,z,rho,x0,absl,buf,nbuf);

      amat[nc] = a;
      zmat[nc] = z;

      Int_t nnc = (fMixtComp[nc]>=20) ? fMixtComp[nc]%20 : fMixtComp[nc]%10;
 
      am += fMixtProp[nc]*((fMixtComp[nc]>=20) ? apure[nnc] : amol[nnc]); 
      density += fMixtProp[nc]*rho;  // density of the mixture
      
    }

  // mixture proportions by weight!

  for(nc = 0;nc<fNoComp;nc++)
    {

      Int_t nnc = (fMixtComp[nc]>=20) ? fMixtComp[nc]%20 : fMixtComp[nc]%10;

      wmat[nc] = fMixtProp[nc]*((fMixtComp[nc]>=20) ? 
                 apure[nnc] : amol[nnc])/am;

    } 

  // Drift gases 1 - nonsensitive, 2 - sensitive

  AliMixture(31,"Drift gas 1",amat,zmat,density,fNoComp,wmat);
  AliMixture(32,"Drift gas 2",amat,zmat,density,fNoComp,wmat);


  // Air

  amat[0] = 14.61;
  zmat[0] = 7.3;
  density = 0.001205;

  AliMaterial(24,"Air",amat[0],zmat[0],density,999.,999.); 


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

  AliMixture(34,"Kevlar",amat,zmat,density,-4,wmat);  

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

  density = 0.03;

  
  AliMixture(35,"NOMEX",amat,zmat,density,-4,wmat);

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

  AliMixture(36,"Makrolon",amat,zmat,density,-3,wmat);
  
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
  
  AliMixture(37, "Mylar",amat,zmat,density,-3,wmat); 

  // SiO2 - used later for the glass fiber

  amat[0]=28.086;
  amat[1]=15.9994;

  zmat[0]=14.;
  zmat[1]=8.;

  wmat[0]=1.;
  wmat[1]=2.;


  AliMixture(38,"SiO2",amat,zmat,2.2,-2,wmat); //SiO2 - quartz (rho=2.2)

  // Al

  amat[0] = 26.98;
  zmat[0] = 13.;

  density = 2.7;

  AliMaterial(40,"Al",amat[0],zmat[0],density,999.,999.);

  // Si

  amat[0] = 28.086;
  zmat[0] = 14.;

  density = 2.33;

  AliMaterial(41,"Si",amat[0],zmat[0],density,999.,999.);

  // Cu

  amat[0] = 63.546;
  zmat[0] = 29.;

  density = 8.96;

  AliMaterial(42,"Cu",amat[0],zmat[0],density,999.,999.);

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

  AliMixture(43, "Tedlar",amat,zmat,density,-3,wmat);  


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

  AliMixture(44,"Plexiglas",amat,zmat,density,-3,wmat);

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

  AliMixture(45,"Epoxy",amat,zmat,density,-3,wmat);

  // Carbon

  amat[0]=12.011;
  zmat[0]=6.;
  density= 2.265;

  AliMaterial(46,"C",amat[0],zmat[0],density,999.,999.);

  // get epoxy

  gMC->Gfmate((*fIdmate)[45],namate,amat[1],zmat[1],rho,x0,absl,buf,nbuf);

  // Carbon fiber

  wmat[0]=0.644; // by weight!
  wmat[1]=0.356;

  density=0.5*(1.25+2.265);

  AliMixture(47,"Cfiber",amat,zmat,density,2,wmat);

  // get SiO2

  gMC->Gfmate((*fIdmate)[38],namate,amat[0],zmat[0],rho,x0,absl,buf,nbuf); 

  wmat[0]=0.725; // by weight!
  wmat[1]=0.275;

  density=1.7;

  AliMixture(39,"G10",amat,zmat,density,2,wmat);

 


  //----------------------------------------------------------
  // tracking media for gases
  //----------------------------------------------------------

  AliMedium(0, "Air", 24, 0, iSXFLD, sXMGMX, 10., 999., .1, .01, .1);
  AliMedium(1, "Drift gas 1", 31, 0, iSXFLD, sXMGMX, 10., 999.,.1,.001, .001);
  AliMedium(2, "Drift gas 2", 32, 1, iSXFLD, sXMGMX, 10., 999.,.1,.001, .001);
  AliMedium(3,"CO2",10,0, iSXFLD, sXMGMX, 10., 999.,.1, .001, .001); 

  //-----------------------------------------------------------  
  // tracking media for solids
  //-----------------------------------------------------------
  
  AliMedium(4,"Al",40,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(5,"Kevlar",34,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(6,"Nomex",35,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(7,"Makrolon",36,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(8,"Mylar",37,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(9,"Tedlar",43,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(10,"Cu",42,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(11,"Si",41,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(12,"G10",39,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(13,"Plexiglas",44,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
  AliMedium(14,"Epoxy",45,0, iSXFLD, sXMGMX, 10., 999., .1, .0005, .001);
  AliMedium(15,"Cfiber",47,0, iSXFLD, sXMGMX, 10., 999., .1, .001, .001);
    
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
  if (fActiveSectors) delete [] fActiveSectors;
  fActiveSectors = new Bool_t[fTPCParam->GetNSector()];
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
  if (fActiveSectors) delete [] fActiveSectors;
  fActiveSectors = new Bool_t[fTPCParam->GetNSector()];
  if (flag) {
    for (Int_t i=0;i<fTPCParam->GetNSector();i++) fActiveSectors[i]=kTRUE;
    return;
  }
  for (Int_t i=0;i<fTPCParam->GetNSector();i++) fActiveSectors[i]=kFALSE;
  TBranch * branch=0;
  if (TreeH() == 0x0)
   {
     Fatal("SetActiveSectors","Can not find TreeH in folder");
     return;
   }
  if (fHitType>1) branch = TreeH()->GetBranch("TPC2");
  else branch = TreeH()->GetBranch("TPC");
  Stat_t ntracks = TreeH()->GetEntries();
  // loop over all hits
  if (GetDebug()) cout<<"\nAliTPC::SetActiveSectors():  Got "<<ntracks<<" tracks\n";
  
  for(Int_t track=0;track<ntracks;track++)
   {
    ResetHits();
    //
    if (fTrackHits && fHitType&4) {
      TBranch * br1 = TreeH()->GetBranch("fVolumes");
      TBranch * br2 = TreeH()->GetBranch("fNVolumes");
      br1->GetEvent(track);
      br2->GetEvent(track);
      Int_t *volumes = fTrackHits->GetVolumes();
      for (Int_t j=0;j<fTrackHits->GetNVolumes(); j++)
	fActiveSectors[volumes[j]]=kTRUE;
    }
    
    //
    if (fTrackHitsOld && fHitType&2) {
      TBranch * br = TreeH()->GetBranch("fTrackHitsInfo");
      br->GetEvent(track);
      AliObjectArray * ar = fTrackHitsOld->fTrackHitsInfo;
      for (UInt_t j=0;j<ar->GetSize();j++){
	fActiveSectors[((AliTrackHitsInfo*)ar->At(j))->fVolumeID] =kTRUE;
      } 
    }    
  }
}  




//_____________________________________________________________________________
void AliTPC::Digits2Raw()
{
// convert digits of the current event to raw data

  static const Int_t kThreshold = 0;
  static const Bool_t kCompress = kTRUE;

  fLoader->LoadDigits();
  TTree* digits = fLoader->TreeD();
  if (!digits) {
    Error("Digits2Raw", "no digits tree");
    return;
  }

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

  if (kCompress) {
    Info("Digits2Raw", "compressing raw data");
    rawWriter.RawDataCompDecompress(kTRUE);
    gSystem->Unlink("Statistics");
  }
}



//______________________________________________________________________
AliDigitizer* AliTPC::CreateDigitizer(AliRunDigitizer* manager) const
{
  return new AliTPCDigitizer(manager);
}
//__
void AliTPC::SDigits2Digits2(Int_t /*eventnumber*/)  
{
  //create digits from summable digits
  GenerNoise(500000); //create teble with noise

  //conect tree with sSDigits
  TTree *t = fLoader->TreeS();

  if (t == 0x0) 
   {
     fLoader->LoadSDigits("READ");
     t = fLoader->TreeS();
     if (t == 0x0)
      {
        Error("SDigits2Digits2","Can not get input TreeS");
        return;
      }
   }
  
  if (fLoader->TreeD() == 0x0) fLoader->MakeTree("D");
  
  AliSimDigits digarr, *dummy=&digarr;
  TBranch* sdb = t->GetBranch("Segment");
  if (sdb == 0x0)
   {
     Error("SDigits2Digits2","Can not find branch with segments in TreeS.");
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
      cerr<<"AliTPC warning: invalid segment ID ! "<<digarr.GetID()<<endl;
      continue;
    }
    if (!IsSectorActive(sec)) 
     {
//       cout<<n<<" NOT Active \n";
       continue;
     }
    else
     {
//       cout<<n<<" Active \n";
     }
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
    Int_t saturation = fTPCParam->GetADCSat();
    //use internal structure of the AliDigits - for speed reason
    //if you cahnge implementation
    //of the Alidigits - it must be rewriten -
    for (Int_t i= 0; i<nelems; i++){
      //      Float_t q = *pamp0;
      //q/=16.;  //conversion faktor
      //Float_t noise= GetNoise(); 
      //q+=noise;      
      //q= TMath::Nint(q);
      Float_t q = TMath::Nint(Float_t(*pamp0)/16.+GetNoise());
      if (q>zerosup){
	if (q>saturation) q=saturation;      
	*pamp1=(Short_t)q;
	//if (ptracks0[0]==0)
	//  ptracks1[0]=1;
	//else
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
    // cerr<<sec<<"\t"<<row<<"\n";   
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
   
  //   cerr<<"Setting default parameters...\n";

  // Set response functions

  //
  AliRunLoader* rl = (AliRunLoader*)fLoader->GetEventFolder()->FindObject(AliRunLoader::GetRunLoaderName());
  rl->CdGAFile();
  AliTPCParamSR *param=(AliTPCParamSR*)gDirectory->Get("75x40_100x60");
  if(param){
    printf("You are using 2 pad-length geom hits with 3 pad-lenght geom digits...\n");
    delete param;
    param = new AliTPCParamSR();
  }
  else {
    param=(AliTPCParamSR*)gDirectory->Get("75x40_100x60_150x60");
  }
  if(!param){
    printf("No TPC parameters found\n");
    exit(4);
  }


  AliTPCPRF2D    * prfinner   = new AliTPCPRF2D;
  AliTPCPRF2D    * prfouter1   = new AliTPCPRF2D;
  AliTPCPRF2D    * prfouter2   = new AliTPCPRF2D;  
  AliTPCRF1D     * rf    = new AliTPCRF1D(kTRUE);
  rf->SetGauss(param->GetZSigma(),param->GetZWidth(),1.);
  rf->SetOffset(3*param->GetZSigma());
  rf->Update();
  
  TDirectory *savedir=gDirectory;
  TFile *f=TFile::Open("$ALICE_ROOT/TPC/AliTPCprf2d.root");
  if (!f->IsOpen()) { 
    cerr<<"Can't open $ALICE_ROOT/TPC/AliTPCprf2d.root !\n" ;
     exit(3);
  }

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

  fLoader->LoadHits("read");
  fLoader->LoadDigits("recreate");
  AliRunLoader* runLoader = fLoader->GetRunLoader(); 

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);
    SetActiveSectors();   
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
  if (fLoader->TreeH() == 0x0)
   {
     if(fLoader->LoadHits())
      {
        Error("Hits2Digits","Can not load hits.");
      }
   }
  SetTreeAddress();
  
  if (fLoader->TreeD() == 0x0 ) 
   {
     fLoader->MakeTree("D");
     if (fLoader->TreeD() == 0x0 ) 
      {
       Error("Hits2Digits","Can not get TreeD");
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

  //  cerr<<"Digitizing TPC -- normal digits...\n";

 for(Int_t isec=0;isec<fTPCParam->GetNSector();isec++) 
  if (IsSectorActive(isec)) 
   {
    if (fDebug) Info("Hits2Digits","Sector %d is active.",isec);
    Hits2DigitsSector(isec);
   }
  else
   {
    if (fDebug) Info("Hits2Digits","Sector %d is NOT active.",isec);
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
//  AliRunLoader* rl = (AliRunLoader*)fLoader->GetEventFolder()->FindObject(AliRunLoader::fgkRunLoaderName);

  AliRunLoader* rl = fLoader->GetRunLoader();

  rl->GetEvent(eventnumber);
  if (fLoader->TreeH() == 0x0)
   {
     if(fLoader->LoadHits())
      {
        Error("Hits2Digits","Can not load hits.");
        return;
      }
   }
  SetTreeAddress();


  if (fLoader->TreeS() == 0x0 ) 
   {
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

  //  cerr<<"Digitizing TPC -- summable digits...\n"; 

  fDigitsSwitch=1; // summable digits
  
    // set zero suppression to "0"

  fTPCParam->SetZeroSup(0);

 for(Int_t isec=0;isec<fTPCParam->GetNSector();isec++) 
  if (IsSectorActive(isec)) 
   {
//    cout<<"Sector "<<isec<<" is active\n";
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

  TTree *tH = TreeH(); // pointer to the hits tree
  if (tH == 0x0)
   {
     Fatal("Hits2DigitsSector","Can not find TreeH in folder");
     return;
   }

  Stat_t ntracks = tH->GetEntries();

  if( ntracks > 0){

  //------------------------------------------- 
  //  Only if there are any tracks...
  //-------------------------------------------

    TObjArray **row;
    
    //printf("*** Processing sector number %d ***\n",isec);

      Int_t nrows =fTPCParam->GetNRow(isec);

      row= new TObjArray* [nrows+2]; // 2 extra rows for cross talk
    
      MakeSector(isec,nrows,tH,ntracks,row);

      //--------------------------------------------------------
      //   Digitize this sector, row by row
      //   row[i] is the pointer to the TObjArray of AliTPCFastVectors,
      //   each one containing electrons accepted on this
      //   row, assigned into tracks
      //--------------------------------------------------------

      Int_t i;

      if (fDigitsArray->GetTree()==0) 
       {
         Fatal("Hits2DigitsSector","Tree not set in fDigitsArray");
       }

      for (i=0;i<nrows;i++){

	AliDigits * dig = fDigitsArray->CreateRow(isec,i); 

	DigitizeRow(i,isec,row);

	fDigitsArray->StoreRow(isec,i);

	Int_t ndig = dig->GetDigitSize(); 
	
	if (gDebug > 10) 
	printf("*** Sector, row, compressed digits %d %d %d ***\n",isec,i,ndig);        
 	
        fDigitsArray->ClearRow(isec,i);  

   
       } // end of the sector digitization

      for(i=0;i<nrows+2;i++){
        row[i]->Delete();  
        delete row[i];   
      }
      
       delete [] row; // delete the array of pointers to TObjArray-s
        
  } // ntracks >0

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
  //  Int_t nrows =fTPCParam->GetNRow(isec);
  fCurrentIndex[1]= isec;
  

  Int_t nofPads = fTPCParam->GetNPads(isec,irow);
  Int_t nofTbins = fTPCParam->GetMaxTBin();
  Int_t indexRange[4];
  //
  //  Integrated signal for this row
  //  and a single track signal
  //    

  AliTPCFastMatrix *m1 = new AliTPCFastMatrix(0,nofPads,0,nofTbins); // integrated
  AliTPCFastMatrix *m2 = new AliTPCFastMatrix(0,nofPads,0,nofTbins); // single
  //
  AliTPCFastMatrix &total  = *m1;

  //  Array of pointers to the label-signal list

  Int_t nofDigits = nofPads*nofTbins; // number of digits for this row
  Float_t  **pList = new Float_t* [nofDigits]; 

  Int_t lp;
  Int_t i1;   
  for(lp=0;lp<nofDigits;lp++)pList[lp]=0; // set all pointers to NULL
  //
  //calculate signal 
  //
  //Int_t row1 = TMath::Max(irow-fTPCParam->GetNCrossRows(),0);
  //Int_t row2 = TMath::Min(irow+fTPCParam->GetNCrossRows(),nrows-1);
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
    Float_t *pq = &(total.UncheckedAt(0,it));
    for(Int_t ip=0;ip<nofPads;ip++){
      gi++;
      Float_t q=*pq;      
      pq++;
      if(fDigitsSwitch == 0){
	q+=GetNoise();
        if(q <=fzerosup) continue; // do not fill zeros
        q = TMath::Nint(q);
        if(q > fTPCParam->GetADCSat()) q = fTPCParam->GetADCSat();  // saturation

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
      if (fDigitsArray->IsSimulated())
	{
	 ((AliSimDigits*)dig)->SetTrackIDFast(tracks[0],it,ip,0);
	 ((AliSimDigits*)dig)->SetTrackIDFast(tracks[1],it,ip,1);
	 ((AliSimDigits*)dig)->SetTrackIDFast(tracks[2],it,ip,2);
	}
     
    
    } // end of loop over time buckets
  }  // end of lop over pads 

  //
  //  This row has been digitized, delete nonused stuff
  //

  for(lp=0;lp<nofDigits;lp++){
    if(pList[lp]) delete [] pList[lp];
  }
  
  delete [] pList;

  delete m1;
  delete m2;
  //  delete m3;

} // end of DigitizeRow

//_____________________________________________________________________________

Float_t AliTPC::GetSignal(TObjArray *p1, Int_t ntr, 
             AliTPCFastMatrix *m1, AliTPCFastMatrix *m2,Int_t *indexRange)
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

  AliTPCFastVector *tv;

  tv = (AliTPCFastVector*)p1->At(ntr); // pointer to a track
  AliTPCFastVector &v = *tv;
  
  Float_t label = v(0);
  Int_t centralPad = (fTPCParam->GetNPads(fCurrentIndex[1],fCurrentIndex[3]-1)-1)/2;

  Int_t nElectrons = (tv->GetNrows()-1)/4;
  indexRange[0]=9999; // min pad
  indexRange[1]=-1; // max pad
  indexRange[2]=9999; //min time
  indexRange[3]=-1; // max time

  AliTPCFastMatrix &signal = *m1;
  AliTPCFastMatrix &total = *m2;
  //
  //  Loop over all electrons
  //
  for(Int_t nel=0; nel<nElectrons; nel++){
    Int_t idx=nel*4;
    Float_t aval =  v(idx+4);
    Float_t eltoadcfac=aval*fTPCParam->GetTotalNormFac(); 
    Float_t xyz[3]={v(idx+1),v(idx+2),v(idx+3)};
    Int_t n = ((AliTPCParamSR*)fTPCParam)->CalcResponseFast(xyz,fCurrentIndex,fCurrentIndex[3]);

    Int_t *index = fTPCParam->GetResBin(0);  
    Float_t *weight = & (fTPCParam->GetResWeight(0));

    if (n>0) for (Int_t i =0; i<n; i++){       
       Int_t pad=index[1]+centralPad;  //in digit coordinates central pad has coordinate 0

         if (pad>=0){
	 Int_t time=index[2];	 
         Float_t qweight = *(weight)*eltoadcfac;
	 
	 if (m1!=0) signal.UncheckedAt(pad,time)+=qweight;
         total.UncheckedAt(pad,time)+=qweight;
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
void AliTPC::GetList(Float_t label,Int_t np,AliTPCFastMatrix *m,
                     Int_t *indexRange, Float_t **pList)
{
  //----------------------------------------------------------------------
  //  Updates the list of tracks contributing to digits for a given row
  //----------------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  AliTPCFastMatrix &signal = *m;

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
        else{

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

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  Float_t gasgain = fTPCParam->GetGasGain();
  Int_t i;
  Float_t xyz[4]; 

  AliTPChit *tpcHit; // pointer to a sigle TPC hit    
  //MI change
  TBranch * branch=0;
  if (fHitType>1) branch = TH->GetBranch("TPC2");
  else branch = TH->GetBranch("TPC");

 
  //----------------------------------------------
  // Create TObjArray-s, one for each row,
  // each TObjArray will store the AliTPCFastVectors
  // of electrons, one AliTPCFastVectors per each track.
  //---------------------------------------------- 
    
  Int_t *nofElectrons = new Int_t [nrows+2]; // electron counter for each row
  AliTPCFastVector **tracks = new AliTPCFastVector* [nrows+2]; //pointers to the track vectors

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

	currentTrack = tpcHit->Track(); // track number


        if(currentTrack != previousTrack){
                          
           // store already filled fTrack
              
	   for(i=0;i<nrows+2;i++){
             if(previousTrack != -1){
	       if(nofElectrons[i]>0){
                 AliTPCFastVector &v = *tracks[i];
		 v(0) = previousTrack;
                 tracks[i]->ResizeTo(4*nofElectrons[i]+1); // shrink if necessary
	         row[i]->Add(tracks[i]);                     
	       }
               else{
                 delete tracks[i]; // delete empty AliTPCFastVector
                 tracks[i]=0;
	       }
	     }

             nofElectrons[i]=0;
             tracks[i] = new AliTPCFastVector(481); // AliTPCFastVectors for the next fTrack

	   } // end of loop over rows
	       
           previousTrack=currentTrack; // update track label 
	}
	   
	Int_t qI = (Int_t) (tpcHit->fQ); // energy loss (number of electrons)

       //---------------------------------------------------
       //  Calculate the electron attachment probability
       //---------------------------------------------------


        Float_t time = 1.e6*(fTPCParam->GetZLength()-TMath::Abs(tpcHit->Z()))
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
	  xyz[0]=tpcHit->X();
	  xyz[1]=tpcHit->Y();
	  xyz[2]=tpcHit->Z();	
	  //
	  // protection for the nonphysical avalanche size (10**6 maximum)
	  //  
          Double_t rn=TMath::Max(gRandom->Rndm(0),1.93e-22);
	  xyz[3]= (Float_t) (-gasgain*TMath::Log(rn)); 
	  index[0]=1;
	  
	  TransportElectron(xyz,index);    
	  Int_t rowNumber;
	  fTPCParam->GetPadRow(xyz,index); 
	  // row 0 - cross talk from the innermost row
	  // row fNRow+1 cross talk from the outermost row
	  rowNumber = index[2]+1; 
	  //transform position to local digit coordinates
	  //relative to nearest pad row 
	  if ((rowNumber<0)||rowNumber>fTPCParam->GetNRow(isec)+1) continue;
          Float_t x1,y1;
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
          if(x1>y1) xyz[3]*=TMath::Max(1.e-6,(y1-x1+1.));	
       
	  nofElectrons[rowNumber]++;	  
	  //----------------------------------
	  // Expand vector if necessary
	  //----------------------------------
	  if(nofElectrons[rowNumber]>120){
	    Int_t range = tracks[rowNumber]->GetNrows();
	    if((nofElectrons[rowNumber])>(range-1)/4){
        
	      tracks[rowNumber]->ResizeTo(range+400); // Add 100 electrons
	    }
	  }
	  
          AliTPCFastVector &v = *tracks[rowNumber];
	  Int_t idx = 4*nofElectrons[rowNumber]-3;
	  Real_t * position = &(((AliTPCFastVector&)v).UncheckedAt(idx)); //make code faster
	  memcpy(position,xyz,4*sizeof(Float_t));
 
	} // end of loop over electrons

        tpcHit = (AliTPChit*)NextHit();
        
      } // end of loop over hits
    } // end of loop over tracks

    //
    //   store remaining track (the last one) if not empty
    //

     for(i=0;i<nrows+2;i++){
       if(nofElectrons[i]>0){
          AliTPCFastVector &v = *tracks[i];
	  v(0) = previousTrack;
          tracks[i]->ResizeTo(4*nofElectrons[i]+1); // shrink if necessary
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
  Int_t i;
  //
  if(fDebug) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" TPC_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

//_____________________________________________________________________________
void AliTPC::MakeBranch(Option_t* option)
{
  //
  // Create Tree branches for the TPC.
  //
  if(GetDebug()) Info("MakeBranch","");
  Int_t buffersize = 4000;
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  
  const char *h = strstr(option,"H");

  if ( h && (fHitType<=1) && (fHits == 0x0)) fHits = new TClonesArray("AliTPChit", 176);//skowron 20.06.03
  
  AliDetector::MakeBranch(option);

  const char *d = strstr(option,"D");
 
  if (fDigits   && fLoader->TreeD() && d) 
   {
      MakeBranchInTree(gAlice->TreeD(), branchname, &fDigits, buffersize, 0);
   }	

  if (fHitType>1) MakeBranch2(option,0); // MI change 14.09.2000
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
void AliTPC::SetSecAL(Int_t sec)
{
  //---------------------------------------------------
  // Activate/deactivate selection for lower sectors
  //---------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------
  fSecAL = sec;
}

//_____________________________________________________________________________
void AliTPC::SetSecAU(Int_t sec)
{
  //----------------------------------------------------
  // Activate/deactivate selection for upper sectors
  //---------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------
  fSecAU = sec;
}

//_____________________________________________________________________________
void AliTPC::SetSecLows(Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5, Int_t s6)
{
  //----------------------------------------
  // Select active lower sectors
  //----------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  fSecLows[0] = s1;
  fSecLows[1] = s2;
  fSecLows[2] = s3;
  fSecLows[3] = s4;
  fSecLows[4] = s5;
  fSecLows[5] = s6;
}

//_____________________________________________________________________________
void AliTPC::SetSecUps(Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5, Int_t s6,
                       Int_t s7, Int_t s8 ,Int_t s9 ,Int_t s10, 
                       Int_t s11 , Int_t s12)
{
  //--------------------------------
  // Select active upper sectors
  //--------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  fSecUps[0] = s1;
  fSecUps[1] = s2;
  fSecUps[2] = s3;
  fSecUps[3] = s4;
  fSecUps[4] = s5;
  fSecUps[5] = s6;
  fSecUps[6] = s7;
  fSecUps[7] = s8;
  fSecUps[8] = s9;
  fSecUps[9] = s10;
  fSecUps[10] = s11;
  fSecUps[11] = s12;
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
//____________________________________________________________________________
void AliTPC::SetGasMixt(Int_t nc,Int_t c1,Int_t c2,Int_t c3,Float_t p1,
                           Float_t p2,Float_t p3)
{

  // gax mixture definition

 fNoComp = nc;
 
 fMixtComp[0]=c1;
 fMixtComp[1]=c2;
 fMixtComp[2]=c3;

 fMixtProp[0]=p1;
 fMixtProp[1]=p2;
 fMixtProp[2]=p3; 
 
 
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

  fTPCParam->Transform1to2(xyz,index);
  
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
 
//_____________________________________________________________________________
AliTPChit::AliTPChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
AliHit(shunt,track)
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
}
 
//________________________________________________________________________
// Additional code because of the AliTPCTrackHitsV2

void AliTPC::MakeBranch2(Option_t *option,const char */*file*/)
{
  //
  // Create a new branch in the current Root Tree
  // The branch of fHits is automatically split
  // MI change 14.09.2000
  if(GetDebug()) Info("MakeBranch2","");
  if (fHitType<2) return;
  char branchname[10];
  sprintf(branchname,"%s2",GetName());  
  //
  // Get the pointer to the header
  const char *cH = strstr(option,"H");
  //
  if (fTrackHits   && TreeH() && cH && fHitType&4) 
   {
    if(GetDebug()) Info("MakeBranch2","Making branch for Type 4 Hits");
    TreeH()->Branch(branchname,"AliTPCTrackHitsV2",&fTrackHits,fBufferSize,99);
   }	

  if (fTrackHitsOld   && TreeH() && cH && fHitType&2) 
   {    
    if(GetDebug()) Info("MakeBranch2","Making branch for Type 2 Hits");
    AliObjectBranch * branch = new AliObjectBranch(branchname,"AliTPCTrackHits",&fTrackHitsOld, 
                                                   TreeH(),fBufferSize,99);
    TreeH()->GetListOfBranches()->Add(branch);
   }	
}

void AliTPC::SetTreeAddress()
{
//Sets tree address for hits  
  if (fHitType<=1)
   {
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
  if(GetDebug()) Info("SetTreeAddress2","");
  
  TBranch *branch;
  char branchname[20];
  sprintf(branchname,"%s2",GetName());
  //
  // Branch address for hit tree
  TTree *treeH = TreeH();
  if ((treeH)&&(fHitType&4)) {
    branch = treeH->GetBranch(branchname);
    if (branch) 
     {
       branch->SetAddress(&fTrackHits);
       if (GetDebug()) Info("SetTreeAddress2","fHitType&4 Setting");
     }
    else 
    if (GetDebug()) Info("SetTreeAddress2","fHitType&4 Failed (can not find branch)");
    
  }
  if ((treeH)&&(fHitType&2)) {
    branch = treeH->GetBranch(branchname);
    if (branch) 
     {
       branch->SetAddress(&fTrackHitsOld);
       if (GetDebug()) Info("SetTreeAddress2","fHitType&2 Setting");
     }
    else if (GetDebug()) 
      Info("SetTreeAddress2","fHitType&2 Failed (can not find branch)");
  }
  //set address to TREETR

  TTree *treeTR = TreeTR();
  if (treeTR && fTrackReferences) {
    branch = treeTR->GetBranch(GetName());
    if (branch) branch->SetAddress(&fTrackReferences);
  }

}

void AliTPC::FinishPrimary()
{
  if (fTrackHits &&fHitType&4)      fTrackHits->FlushHitStack();  
  if (fTrackHitsOld && fHitType&2)  fTrackHitsOld->FlushHitStack();  
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
  //AliTPChit *hit = (AliTPChit*)fHits->UncheckedAt(fNhits-1);
  //if (hit->fTrack!=rtrack)
  //  cout<<"bad track number\n";
  if (fTrackHits && fHitType&4) 
    fTrackHits->AddHitKartez(vol[0],rtrack, hits[0],
                             hits[1],hits[2],(Int_t)hits[3]);
  if (fTrackHitsOld &&fHitType&2 ) 
    fTrackHitsOld->AddHitKartez(vol[0],rtrack, hits[0],
                                hits[1],hits[2],(Int_t)hits[3]);
  
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
  if (fTrackHitsOld && fHitType&2) fTrackHitsOld->Clear();

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
    gAlice->ResetHits();
    TreeH()->GetEvent(track);
  }
  //
  if (fTrackHits && fHitType&4) {
    fTrackHits->First();
    return fTrackHits->GetHit();
  }
  if (fTrackHitsOld && fHitType&2) {
    fTrackHitsOld->First();
    return fTrackHitsOld->GetHit();
  }

  else return 0;
}

AliHit* AliTPC::NextHit2()
{
  //
  //Return the next hit for the current track


  if (fTrackHitsOld && fHitType&2) {
    fTrackHitsOld->Next();
    return fTrackHitsOld->GetHit();
  }
  if (fTrackHits) {
    fTrackHits->Next();
    return fTrackHits->GetHit();
  }
  else 
    return 0;
}

void AliTPC::LoadPoints(Int_t)
{
  //
  Int_t a = 0;
  /*  if(fHitType==1) return AliDetector::LoadPoints(a);
  LoadPoints2(a);
  */
  if(fHitType==1) AliDetector::LoadPoints(a);
  else LoadPoints2(a);
   
  // LoadPoints3(a);

}


void AliTPC::RemapTrackHitIDs(Int_t *map)
{
  //
  // remapping
  //
  if (!fTrackHits) return;
  
  if (fTrackHitsOld && fHitType&2){
    AliObjectArray * arr = fTrackHitsOld->fTrackHitsInfo;
    for (UInt_t i=0;i<arr->GetSize();i++){
      AliTrackHitsInfo * info = (AliTrackHitsInfo *)(arr->At(i));
      info->fTrackID = map[info->fTrackID];
    }
  }
  if (fTrackHitsOld && fHitType&4){
    TClonesArray * arr = fTrackHits->GetArray();;
    for (Int_t i=0;i<arr->GetEntriesFast();i++){
      AliTrackHitsParamV2 * info = (AliTrackHitsParamV2 *)(arr->At(i));
      info->fTrackID = map[info->fTrackID];
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
  if (fTrackHitsOld && fHitType&2) {
    TBranch * br = TreeH()->GetBranch("fTrackHitsInfo");
    br->GetEvent(track);
    AliObjectArray * ar = fTrackHitsOld->fTrackHitsInfo;
    for (UInt_t j=0;j<ar->GetSize();j++){
      if (  ((AliTrackHitsInfo*)ar->At(j))->fVolumeID==id) return kTRUE;
    } 
  }

  if (fTrackHits && fHitType&4) {
    TBranch * br1 = TreeH()->GetBranch("fVolumes");
    TBranch * br2 = TreeH()->GetBranch("fNVolumes");    
    br2->GetEvent(track);
    br1->GetEvent(track);    
    Int_t *volumes = fTrackHits->GetVolumes();
    Int_t nvolumes = fTrackHits->GetNVolumes();
    if (!volumes && nvolumes>0) {
      printf("Problematic track\t%d\t%d",track,nvolumes);
      return kFALSE;
    }
    for (Int_t j=0;j<nvolumes; j++)
      if (volumes[j]==id) return kTRUE;    
  }

  if (fHitType&1) {
    TBranch * br = TreeH()->GetBranch("fSector");
    br->GetEvent(track);
    for (Int_t j=0;j<fHits->GetEntriesFast();j++){
      if (  ((AliTPChit*)fHits->At(j))->fSector==id) return kTRUE;
    } 
  }
  return kFALSE;  

}

//_____________________________________________________________________________
void AliTPC::LoadPoints2(Int_t)
{
  //
  // Store x, y, z of all hits in memory
  //
  if (fTrackHits == 0 && fTrackHitsOld==0) return;
  //
  Int_t nhits =0;
  if (fHitType&4) nhits = fTrackHits->GetEntriesFast();
  if (fHitType&2) nhits = fTrackHitsOld->GetEntriesFast();
  
  if (nhits == 0) return;
  Int_t tracks = gAlice->GetMCApp()->GetNtrack();
  if (fPoints == 0) fPoints = new TObjArray(tracks);
  AliHit *ahit;
  //
  Int_t *ntrk=new Int_t[tracks];
  Int_t *limi=new Int_t[tracks];
  Float_t **coor=new Float_t*[tracks];
  for(Int_t i=0;i<tracks;i++) {
    ntrk[i]=0;
    coor[i]=0;
    limi[i]=0;
  }
  //
  AliPoints *points = 0;
  Float_t *fp=0;
  Int_t trk;
  Int_t chunk=nhits/4+1;
  //
  // Loop over all the hits and store their position
  //
  ahit = FirstHit2(-1);
  while (ahit){
    trk=ahit->GetTrack();
    if(ntrk[trk]==limi[trk]) {
      //
      // Initialise a new track
      fp=new Float_t[3*(limi[trk]+chunk)];
      if(coor[trk]) {
	memcpy(fp,coor[trk],sizeof(Float_t)*3*limi[trk]);
	delete [] coor[trk];
      }
      limi[trk]+=chunk;
      coor[trk] = fp;
    } else {
      fp = coor[trk];
    }
    fp[3*ntrk[trk]  ] = ahit->X();
    fp[3*ntrk[trk]+1] = ahit->Y();
    fp[3*ntrk[trk]+2] = ahit->Z();
    ntrk[trk]++;
    ahit = NextHit2();
  }



  //
  for(trk=0; trk<tracks; ++trk) {
    if(ntrk[trk]) {
      points = new AliPoints();
      points->SetMarkerColor(GetMarkerColor());
      points->SetMarkerSize(GetMarkerSize());
      points->SetDetector(this);
      points->SetParticle(trk);
      points->SetPolyMarker(ntrk[trk],coor[trk],GetMarkerStyle());
      fPoints->AddAt(points,trk);
      delete [] coor[trk];
      coor[trk]=0;
    }
  }
  delete [] coor;
  delete [] ntrk;
  delete [] limi;
}


//_____________________________________________________________________________
void AliTPC::LoadPoints3(Int_t)
{
  //
  // Store x, y, z of all hits in memory
  // - only intersection point with pad row
  if (fTrackHits == 0) return;
  //
  Int_t nhits = fTrackHits->GetEntriesFast();
  if (nhits == 0) return;
  Int_t tracks = gAlice->GetMCApp()->GetNtrack();
  if (fPoints == 0) fPoints = new TObjArray(2*tracks);
  fPoints->Expand(2*tracks);
  AliHit *ahit;
  //
  Int_t *ntrk=new Int_t[tracks];
  Int_t *limi=new Int_t[tracks];
  Float_t **coor=new Float_t*[tracks];
  for(Int_t i=0;i<tracks;i++) {
    ntrk[i]=0;
    coor[i]=0;
    limi[i]=0;
  }
  //
  AliPoints *points = 0;
  Float_t *fp=0;
  Int_t trk;
  Int_t chunk=nhits/4+1;
  //
  // Loop over all the hits and store their position
  //
  ahit = FirstHit2(-1);
  //for (Int_t hit=0;hit<nhits;hit++) {

  Int_t lastrow = -1;
  while (ahit){
    //    ahit = (AliHit*)fHits->UncheckedAt(hit);
    trk=ahit->GetTrack(); 
    Float_t  x[3]={ahit->X(),ahit->Y(),ahit->Z()};
    Int_t    index[3]={1,((AliTPChit*)ahit)->fSector,0};
    Int_t    currentrow = fTPCParam->GetPadRow(x,index) ;
    if (currentrow!=lastrow){
      lastrow = currentrow;
      //later calculate intersection point           
      if(ntrk[trk]==limi[trk]) {
	//
	// Initialise a new track
	fp=new Float_t[3*(limi[trk]+chunk)];
	if(coor[trk]) {
	  memcpy(fp,coor[trk],sizeof(Float_t)*3*limi[trk]);
	  delete [] coor[trk];
	}
	limi[trk]+=chunk;
	coor[trk] = fp;
      } else {
	fp = coor[trk];
      }
      fp[3*ntrk[trk]  ] = ahit->X();
      fp[3*ntrk[trk]+1] = ahit->Y();
      fp[3*ntrk[trk]+2] = ahit->Z();
      ntrk[trk]++;
    }
    ahit = NextHit2();
  }
  
  //
  for(trk=0; trk<tracks; ++trk) {
    if(ntrk[trk]) {
      points = new AliPoints();
      points->SetMarkerColor(GetMarkerColor()+1);
      points->SetMarkerStyle(5);
      points->SetMarkerSize(0.2);
      points->SetDetector(this);
      points->SetParticle(trk);
      //      points->SetPolyMarker(ntrk[trk],coor[trk],GetMarkerStyle()20);
      points->SetPolyMarker(ntrk[trk],coor[trk],30);
      fPoints->AddAt(points,tracks+trk);
      delete [] coor[trk];
      coor[trk]=0;
    }
  }
  delete [] coor;
  delete [] ntrk;
  delete [] limi;
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
  sprintf(paramName,"75x40_100x60_150x60");
  AliTPCParam *paramTPC=(AliTPCParam*)file->Get(paramName);
  if (paramTPC) {
    //    cout<<"TPC parameters "<<paramName<<" found."<<endl;
  } else {
    cerr<<"TPC parameters not found. Create new (they may be incorrect)."
	<<endl;    
    paramTPC = new AliTPCParamSR;
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


