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

/*
$Log$
Revision 1.26  2001/01/10 07:59:43  kowal2
Corrections to load points from the noncompressed hits.

Revision 1.25  2000/11/02 07:25:31  kowal2
Changes due to the new hit structure.
Memory leak removed.

Revision 1.24  2000/10/05 16:06:09  kowal2
Forward declarations. Changes due to a new class AliComplexCluster.

Revision 1.23  2000/10/02 21:28:18  fca
Removal of useless dependecies via forward declarations

Revision 1.22  2000/07/10 20:57:39  hristov
Update of TPC code and macros by M.Kowalski

Revision 1.19.2.4  2000/06/26 07:39:42  kowal2
Changes to obey the coding rules

Revision 1.19.2.3  2000/06/25 08:38:41  kowal2
Splitted from AliTPCtracking

Revision 1.19.2.2  2000/06/16 12:59:28  kowal2
Changed parameter settings

Revision 1.19.2.1  2000/06/09 07:15:07  kowal2

Defaults loaded automatically (hard-wired)
Optional parameters can be set via macro called in the constructor

Revision 1.19  2000/04/18 19:00:59  fca
Small bug fixes to TPC files

Revision 1.18  2000/04/17 09:37:33  kowal2
removed obsolete AliTPCDigitsDisplay.C

Revision 1.17.2.2  2000/04/10 08:15:12  kowal2

New, experimental data structure from M. Ivanov
New tracking algorithm
Different pad geometry for different sectors
Digitization rewritten

Revision 1.17.2.1  2000/04/10 07:56:53  kowal2
Not used anymore - removed

Revision 1.17  2000/01/19 17:17:30  fca
Introducing a list of lists of hits -- more hits allowed for detector now

Revision 1.16  1999/11/05 09:29:23  fca
Accept only signals > 0

Revision 1.15  1999/10/08 06:26:53  fca
Removed ClustersIndex - not used anymore

Revision 1.14  1999/09/29 09:24:33  fca
Introduction of the Copyright and cvs Log

*/

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

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TMatrix.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBS.h>
#include <TObjectTable.h>
#include "TParticle.h"
#include "AliTPC.h"
#include <TFile.h>       
#include "AliRun.h"
#include <iostream.h>
#include <fstream.h>
#include "AliMC.h"
#include "AliMagF.h"


#include "AliTPCParamSR.h"
#include "AliTPCPRF2D.h"
#include "AliTPCRF1D.h"
#include "AliDigits.h"
#include "AliSimDigits.h"
#include "AliTPCTrackHits.h"
#include "AliPoints.h"
#include "AliArrayBranch.h"


#include "AliTPCDigitsArray.h"
#include "AliComplexCluster.h"
#include "AliClusters.h"
#include "AliTPCClustersRow.h"
#include "AliTPCClustersArray.h"

#include "AliTPCcluster.h"
#include "AliTPCclusterer.h"
#include "AliTPCtracker.h"

#include <TInterpreter.h>
#include <TTree.h>



ClassImp(AliTPC) 

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
  //MI changes
  fDigitsArray = 0;
  fClustersArray = 0;
  fTPCParam=0;
  fTrackHits = 0;  
  fHitType = 2;  
  fTPCParam = 0; 
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
  gAlice->AddHitList(fHits);
  //MI change  
  fDigitsArray = 0;
  fClustersArray= 0;
  //
  fTrackHits = new AliTPCTrackHits;  //MI - 13.09.2000
  fTrackHits->SetHitPrecision(0.002);
  fTrackHits->SetStepPrecision(0.003);  
  fTrackHits->SetMaxDistance(100); 

  fHitType = 2;
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
  if (fHitType&2)
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
Int_t AliTPC::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Calculate distance from TPC to mouse on the display
  // Dummy procedure
  //
  return 9999;
}

void AliTPC::Clusters2Tracks(TFile *of) {
  //-----------------------------------------------------------------
  // This is a track finder.
  //-----------------------------------------------------------------
  AliTPCtracker::Clusters2Tracks(fTPCParam,of);
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

  char namate[21]; 
  density = 0.;
  Float_t am=0;
  Int_t nc;
  Float_t rho,absl,X0,buf[1];
  Int_t nbuf;
  Float_t a,z;

  for(nc = 0;nc<fNoComp;nc++)
    {
    
      // retrive material constants
      
      gMC->Gfmate((*fIdmate)[fMixtComp[nc]],namate,a,z,rho,X0,absl,buf,nbuf);

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

  // G10 60% SiO2 + 40% epoxy, I use A and Z for SiO2

  amat[0]=28.086;
  amat[1]=15.9994;

  zmat[0]=14.;
  zmat[1]=8.;

  wmat[0]=1.;
  wmat[1]=2.;

  density = 1.7;

  AliMixture(38,"SiO2",amat,zmat,2.2,-2,wmat); //SiO2 - quartz
 
  gMC->Gfmate((*fIdmate)[38],namate,amat[0],zmat[0],rho,X0,absl,buf,nbuf); 

  AliMaterial(39,"G10",amat[0],zmat[0],density,999.,999.); 

  // Al

  amat[0] = 26.98;
  zmat[0] = 13.;

  density = 2.7;

  AliMaterial(40,"Al",amat[0],zmat[0],density,999.,999.);

  // Si

  amat[0] = 28.086;
  zmat[0] = 14.,

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
    
}


void AliTPC::Digits2Clusters(TFile *of)
{
  //-----------------------------------------------------------------
  // This is a simple cluster finder.
  //-----------------------------------------------------------------
  AliTPCclusterer::Digits2Clusters(fTPCParam,of);
}

extern Double_t SigmaY2(Double_t, Double_t, Double_t);
extern Double_t SigmaZ2(Double_t, Double_t);
//_____________________________________________________________________________
void AliTPC::Hits2Clusters(TFile *of)
{
  //--------------------------------------------------------
  // TPC simple cluster generator from hits
  // obtained from the TPC Fast Simulator
  // The point errors are taken from the parametrization
  //--------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------
  // Adopted to Marian's cluster data structure by I.Belikov, CERN,
  // Jouri.Belikov@cern.ch
  //----------------------------------------------------------------
  
  /////////////////////////////////////////////////////////////////////////////
  //
  //---------------------------------------------------------------------
  //   ALICE TPC Cluster Parameters
  //--------------------------------------------------------------------
       
  

  // Cluster width in rphi
  const Float_t kACrphi=0.18322;
  const Float_t kBCrphi=0.59551e-3;
  const Float_t kCCrphi=0.60952e-1;
  // Cluster width in z
  const Float_t kACz=0.19081;
  const Float_t kBCz=0.55938e-3;
  const Float_t kCCz=0.30428;

  TDirectory *savedir=gDirectory; 

  if (!of->IsOpen()) {
     cerr<<"AliTPC::Hits2Clusters(): output file not open !\n";
     return;
  }

   if(fTPCParam == 0){
     printf("AliTPCParam MUST be created firstly\n");
     return;
   }

  Float_t sigmaRphi,sigmaZ,clRphi,clZ;
  //
  TParticle *particle; // pointer to a given particle
  AliTPChit *tpcHit; // pointer to a sigle TPC hit
  TClonesArray *particles; //pointer to the particle list
  Int_t sector;
  Int_t ipart;
  Float_t xyz[5];
  Float_t pl,pt,tanth,rpad,ratio;
  Float_t cph,sph;
  
  //---------------------------------------------------------------
  //  Get the access to the tracks 
  //---------------------------------------------------------------
  
  TTree *tH = gAlice->TreeH();
  Stat_t ntracks = tH->GetEntries();
  particles=gAlice->Particles();

  //Switch to the output file
  of->cd();

  fTPCParam->Write(fTPCParam->GetTitle());
  AliTPCClustersArray carray;
  carray.Setup(fTPCParam);
  carray.SetClusterType("AliTPCcluster");
  carray.MakeTree();

  Int_t nclusters=0; //cluster counter
  
  //------------------------------------------------------------
  // Loop over all sectors (72 sectors for 20 deg
  // segmentation for both lower and upper sectors)
  // Sectors 0-35 are lower sectors, 0-17 z>0, 17-35 z<0
  // Sectors 36-71 are upper sectors, 36-53 z>0, 54-71 z<0
  //
  // First cluster for sector 0 starts at "0"
  //------------------------------------------------------------
   
  for(Int_t isec=0;isec<fTPCParam->GetNSector();isec++){
    //MI change
    fTPCParam->AdjustCosSin(isec,cph,sph);
    
    //------------------------------------------------------------
    // Loop over tracks
    //------------------------------------------------------------
    
    for(Int_t track=0;track<ntracks;track++){
      ResetHits();
      tH->GetEvent(track);
      //
      //  Get number of the TPC hits
      //
      // nhits=fHits->GetEntriesFast();
      //
     
       tpcHit = (AliTPChit*)FirstHit(-1);

      // Loop over hits
      //
       //   for(Int_t hit=0;hit<nhits;hit++){
       //tpcHit=(AliTPChit*)fHits->UncheckedAt(hit);

       while(tpcHit){
 
	 if (tpcHit->fQ == 0.) {
           tpcHit = (AliTPChit*) NextHit();
           continue; //information about track (I.Belikov)
	 }
	sector=tpcHit->fSector; // sector number


	//	if(sector != isec) continue; //terminate iteration

       if(sector != isec){
	 tpcHit = (AliTPChit*) NextHit();
	 continue; 
       }
	ipart=tpcHit->Track();
	particle=(TParticle*)particles->UncheckedAt(ipart);
	pl=particle->Pz();
	pt=particle->Pt();
	if(pt < 1.e-9) pt=1.e-9;
	tanth=pl/pt;
	tanth = TMath::Abs(tanth);
	rpad=TMath::Sqrt(tpcHit->X()*tpcHit->X() + tpcHit->Y()*tpcHit->Y());
	ratio=0.001*rpad/pt; // pt must be in MeV/c - historical reason

	//   space-point resolutions
	
	sigmaRphi=SigmaY2(rpad,tanth,pt);
	sigmaZ   =SigmaZ2(rpad,tanth   );
	
	//   cluster widths
	
	clRphi=kACrphi-kBCrphi*rpad*tanth+kCCrphi*ratio*ratio;
	clZ=kACz-kBCz*rpad*tanth+kCCz*tanth*tanth;
	
	// temporary protection
	
	if(sigmaRphi < 0.) sigmaRphi=0.4e-3;
	if(sigmaZ < 0.) sigmaZ=0.4e-3;
	if(clRphi < 0.) clRphi=2.5e-3;
	if(clZ < 0.) clZ=2.5e-5;
	
	//
	
	//
	// smearing --> rotate to the 1 (13) or to the 25 (49) sector,
	// then the inaccuracy in a X-Y plane is only along Y (pad row)!
	//
        Float_t xprim= tpcHit->X()*cph + tpcHit->Y()*sph;
	Float_t yprim=-tpcHit->X()*sph + tpcHit->Y()*cph;
	xyz[0]=gRandom->Gaus(yprim,TMath::Sqrt(sigmaRphi));   // y
          Float_t alpha=(isec < fTPCParam->GetNInnerSector()) ?
	  fTPCParam->GetInnerAngle() : fTPCParam->GetOuterAngle();
          Float_t ymax=xprim*TMath::Tan(0.5*alpha);
          if (TMath::Abs(xyz[0])>ymax) xyz[0]=yprim; 
	xyz[1]=gRandom->Gaus(tpcHit->Z(),TMath::Sqrt(sigmaZ)); // z
          if (TMath::Abs(xyz[1])>fTPCParam->GetZLength()) xyz[1]=tpcHit->Z(); 
	xyz[2]=sigmaRphi;                                     // fSigmaY2
	xyz[3]=sigmaZ;                                        // fSigmaZ2
	xyz[4]=tpcHit->fQ;                                    // q

        AliTPCClustersRow *clrow=carray.GetRow(sector,tpcHit->fPadRow);
        if (!clrow) clrow=carray.CreateRow(sector,tpcHit->fPadRow);	

        Int_t tracks[3]={tpcHit->Track(), -1, -1};
	AliTPCcluster cluster(tracks,xyz);

        clrow->InsertCluster(&cluster); nclusters++;

        tpcHit = (AliTPChit*)NextHit();
        

      } // end of loop over hits

    }   // end of loop over tracks

    Int_t nrows=fTPCParam->GetNRow(isec);
    for (Int_t irow=0; irow<nrows; irow++) {
        AliTPCClustersRow *clrow=carray.GetRow(isec,irow);
        if (!clrow) continue;
        carray.StoreRow(isec,irow);
        carray.ClearRow(isec,irow);
    }

  } // end of loop over sectors  

  cerr<<"Number of made clusters : "<<nclusters<<"                        \n";

  carray.GetTree()->Write();

  savedir->cd(); //switch back to the input file
  
} // end of function

//_________________________________________________________________
void AliTPC::Hits2ExactClustersSector(Int_t isec)
{
  //--------------------------------------------------------
  //calculate exact cross point of track and given pad row
  //resulting values are expressed in "digit" coordinata
  //--------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marian Ivanov  GSI Darmstadt, m.ivanov@gsi.de
  //-----------------------------------------------------------------
  //
  if (fClustersArray==0){    
    return;
  }
  //
  TParticle *particle; // pointer to a given particle
  AliTPChit *tpcHit; // pointer to a sigle TPC hit
  TClonesArray *particles; //pointer to the particle list
  Int_t sector,nhits;
  Int_t ipart;
  const Int_t kcmaxhits=30000;
  TVector * xxxx = new TVector(kcmaxhits*4);
  TVector & xxx = *xxxx;
  Int_t maxhits = kcmaxhits;
  //construct array for each padrow
  for (Int_t i=0; i<fTPCParam->GetNRow(isec);i++) 
    fClustersArray->CreateRow(isec,i);
  
  //---------------------------------------------------------------
  //  Get the access to the tracks 
  //---------------------------------------------------------------
  
  TTree *tH = gAlice->TreeH();
  Stat_t ntracks = tH->GetEntries();
  particles=gAlice->Particles();
  Int_t npart = particles->GetEntriesFast();
    
  //------------------------------------------------------------
  // Loop over tracks
  //------------------------------------------------------------
  
  for(Int_t track=0;track<ntracks;track++){
    ResetHits();
    tH->GetEvent(track);
    //
    //  Get number of the TPC hits and a pointer
    //  to the particles
    //
    nhits=fHits->GetEntriesFast();
    //
    // Loop over hits
    //
    Int_t currentIndex=0;
    Int_t lastrow=-1;  //last writen row
    for(Int_t hit=0;hit<nhits;hit++){
      tpcHit=(AliTPChit*)fHits->UncheckedAt(hit);
      if (tpcHit==0) continue;
      sector=tpcHit->fSector; // sector number
      if(sector != isec) continue; 
      ipart=tpcHit->Track();
      if (ipart<npart) particle=(TParticle*)particles->UncheckedAt(ipart);
      
      //find row number

      Float_t  x[3]={tpcHit->X(),tpcHit->Y(),tpcHit->Z()};
      Int_t    index[3]={1,isec,0};
      Int_t    currentrow = fTPCParam->GetPadRow(x,index) ;	
      if (currentrow<0) continue;
      if (lastrow<0) lastrow=currentrow;
      if (currentrow==lastrow){
	if ( currentIndex>=maxhits){
	  maxhits+=kcmaxhits;
	  xxx.ResizeTo(4*maxhits);
	}     
	xxx(currentIndex*4)=x[0];
	xxx(currentIndex*4+1)=x[1];
	xxx(currentIndex*4+2)=x[2];	
	xxx(currentIndex*4+3)=tpcHit->fQ;
	currentIndex++;	
      }
      else 
	if (currentIndex>2){
	  Float_t sumx=0;
	  Float_t sumx2=0;
	  Float_t sumx3=0;
	  Float_t sumx4=0;
	  Float_t sumy=0;
	  Float_t sumxy=0;
	  Float_t sumx2y=0;
	  Float_t sumz=0;
	  Float_t sumxz=0;
	  Float_t sumx2z=0;
	  Float_t sumq=0;
	  for (Int_t index=0;index<currentIndex;index++){
	    Float_t x,x2,x3,x4;
	    x=x2=x3=x4=xxx(index*4);
	    x2*=x;
	    x3*=x2;
	    x4*=x3;
	    sumx+=x;
	    sumx2+=x2;
	    sumx3+=x3;
	    sumx4+=x4;
	    sumy+=xxx(index*4+1);
	    sumxy+=xxx(index*4+1)*x;
	    sumx2y+=xxx(index*4+1)*x2;
	    sumz+=xxx(index*4+2);
	    sumxz+=xxx(index*4+2)*x;
	    sumx2z+=xxx(index*4+2)*x2;	 
	    sumq+=xxx(index*4+3);
	  }
	  Float_t centralPad = (fTPCParam->GetNPads(isec,lastrow)-1)/2;
	  Float_t det=currentIndex*(sumx2*sumx4-sumx3*sumx3)-sumx*(sumx*sumx4-sumx2*sumx3)+
	    sumx2*(sumx*sumx3-sumx2*sumx2);
	  
	  Float_t detay=sumy*(sumx2*sumx4-sumx3*sumx3)-sumx*(sumxy*sumx4-sumx2y*sumx3)+
	    sumx2*(sumxy*sumx3-sumx2y*sumx2);
	  Float_t detaz=sumz*(sumx2*sumx4-sumx3*sumx3)-sumx*(sumxz*sumx4-sumx2z*sumx3)+
	    sumx2*(sumxz*sumx3-sumx2z*sumx2);
	  
	  Float_t detby=currentIndex*(sumxy*sumx4-sumx2y*sumx3)-sumy*(sumx*sumx4-sumx2*sumx3)+
	    sumx2*(sumx*sumx2y-sumx2*sumxy);
	  Float_t detbz=currentIndex*(sumxz*sumx4-sumx2z*sumx3)-sumz*(sumx*sumx4-sumx2*sumx3)+
	    sumx2*(sumx*sumx2z-sumx2*sumxz);
	  
	  Float_t y=detay/det+centralPad;
	  Float_t z=detaz/det;	
	  Float_t by=detby/det; //y angle
	  Float_t bz=detbz/det; //z angle
	  sumy/=Float_t(currentIndex);
	  sumz/=Float_t(currentIndex);
	  AliComplexCluster cl;
	  cl.fX=z;
	  cl.fY=y;
	  cl.fQ=sumq;
	  cl.fSigmaX2=bz;
	  cl.fSigmaY2=by;
	  cl.fTracks[0]=ipart;
	  
	  AliTPCClustersRow * row = (fClustersArray->GetRow(isec,lastrow));
	  if (row!=0) row->InsertCluster(&cl);
	  currentIndex=0;
	  lastrow=currentrow;
	} //end of calculating cluster for given row
	
	
	
    } // end of loop over hits
  }   // end of loop over tracks 
  //write padrows to tree 
  for (Int_t ii=0; ii<fTPCParam->GetNRow(isec);ii++) {
    fClustersArray->StoreRow(isec,ii);    
    fClustersArray->ClearRow(isec,ii);        
  }
  xxxx->Delete();
 
}

//__________________________________________________________________  
void AliTPC::Hits2Digits()  
{ 
 //----------------------------------------------------
 // Loop over all sectors
 //----------------------------------------------------

  if(fTPCParam == 0){
    printf("AliTPCParam MUST be created firstly\n");
    return;
  } 

 for(Int_t isec=0;isec<fTPCParam->GetNSector();isec++) Hits2DigitsSector(isec);

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


  TTree *tH = gAlice->TreeH(); // pointer to the hits tree
  Stat_t ntracks = tH->GetEntries();

  if( ntracks > 0){

  //------------------------------------------- 
  //  Only if there are any tracks...
  //-------------------------------------------

    TObjArray **row;
    
      printf("*** Processing sector number %d ***\n",isec);

      Int_t nrows =fTPCParam->GetNRow(isec);

      row= new TObjArray* [nrows];
    
      MakeSector(isec,nrows,tH,ntracks,row);

      //--------------------------------------------------------
      //   Digitize this sector, row by row
      //   row[i] is the pointer to the TObjArray of TVectors,
      //   each one containing electrons accepted on this
      //   row, assigned into tracks
      //--------------------------------------------------------

      Int_t i;

      if (fDigitsArray->GetTree()==0) fDigitsArray->MakeTree();

      for (i=0;i<nrows;i++){

	AliDigits * dig = fDigitsArray->CreateRow(isec,i); 

	DigitizeRow(i,isec,row);

	fDigitsArray->StoreRow(isec,i);

	Int_t ndig = dig->GetDigitSize(); 
 
	printf("*** Sector, row, compressed digits %d %d %d ***\n",isec,i,ndig);
	
        fDigitsArray->ClearRow(isec,i);  

   
       } // end of the sector digitization

      for(i=0;i<nrows;i++){
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
  Int_t nrows =fTPCParam->GetNRow(isec);
  fCurrentIndex[1]= isec;
  

  Int_t nofPads = fTPCParam->GetNPads(isec,irow);
  Int_t nofTbins = fTPCParam->GetMaxTBin();
  Int_t indexRange[4];
  //
  //  Integrated signal for this row
  //  and a single track signal
  //    
  TMatrix *m1   = new TMatrix(0,nofPads,0,nofTbins); // integrated
  TMatrix *m2   = new TMatrix(0,nofPads,0,nofTbins); // single
  //
  TMatrix &total  = *m1;

  //  Array of pointers to the label-signal list

  Int_t nofDigits = nofPads*nofTbins; // number of digits for this row
  Float_t  **pList = new Float_t* [nofDigits]; 

  Int_t lp;
  Int_t i1;   
  for(lp=0;lp<nofDigits;lp++)pList[lp]=0; // set all pointers to NULL
  //
  //calculate signal 
  //
  Int_t row1 = TMath::Max(irow-fTPCParam->GetNCrossRows(),0);
  Int_t row2 = TMath::Min(irow+fTPCParam->GetNCrossRows(),nrows-1);
  for (Int_t row= row1;row<=row2;row++){
    Int_t nTracks= rows[row]->GetEntries();
    for (i1=0;i1<nTracks;i1++){
      fCurrentIndex[2]= row;
      fCurrentIndex[3]=irow;
      if (row==irow){
	m2->Zero();  // clear single track signal matrix
	Float_t trackLabel = GetSignal(rows[row],i1,m2,m1,indexRange); 
	GetList(trackLabel,nofPads,m2,indexRange,pList);
      }
      else   GetSignal(rows[row],i1,0,m1,indexRange);
    }
  }
         
  Int_t tracks[3];

  AliDigits *dig = fDigitsArray->GetRow(isec,irow);
  for(Int_t ip=0;ip<nofPads;ip++){
    for(Int_t it=0;it<nofTbins;it++){

      Float_t q = total(ip,it);

      Int_t gi =it*nofPads+ip; // global index

      q = gRandom->Gaus(q,fTPCParam->GetNoise()*fTPCParam->GetNoiseNormFac()); 

      q = (Int_t)q;

      if(q <=zerosup) continue; // do not fill zeros
      if(q > fTPCParam->GetADCSat()) q = fTPCParam->GetADCSat();  // saturation

      //
      //  "real" signal or electronic noise (list = -1)?
      //    

      for(Int_t j1=0;j1<3;j1++){
        tracks[j1] = (pList[gi]) ?(Int_t)(*(pList[gi]+j1)) : -1;
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

Float_t AliTPC::GetSignal(TObjArray *p1, Int_t ntr, TMatrix *m1, TMatrix *m2,
                          Int_t *indexRange)
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
  Int_t centralPad = (fTPCParam->GetNPads(fCurrentIndex[1],fCurrentIndex[3])-1)/2;

  Int_t nElectrons = (tv->GetNrows()-1)/4;
  indexRange[0]=9999; // min pad
  indexRange[1]=-1; // max pad
  indexRange[2]=9999; //min time
  indexRange[3]=-1; // max time

  //  Float_t IneffFactor = 0.5; // inefficiency in the gain close to the edge, as above

  TMatrix &signal = *m1;
  TMatrix &total = *m2;
  //
  //  Loop over all electrons
  //
  for(Int_t nel=0; nel<nElectrons; nel++){
    Int_t idx=nel*4;
    Float_t aval =  v(idx+4);
    Float_t eltoadcfac=aval*fTPCParam->GetTotalNormFac(); 
    Float_t xyz[3]={v(idx+1),v(idx+2),v(idx+3)};
    Int_t n = fTPCParam->CalcResponse(xyz,fCurrentIndex,fCurrentIndex[3]);
    
    if (n>0) for (Int_t i =0; i<n; i++){
       Int_t *index = fTPCParam->GetResBin(i);        
       Int_t pad=index[1]+centralPad;  //in digit coordinates central pad has coordinate 0
       if ( ( pad<(fTPCParam->GetNPads(fCurrentIndex[1],fCurrentIndex[3]))) && (pad>0)) {
	 Int_t time=index[2];	 
	 Float_t weight = fTPCParam->GetResWeight(i); //we normalise response to ADC channel
	 weight *= eltoadcfac;
	 
	 if (m1!=0) signal(pad,time)+=weight; 
	 total(pad,time)+=weight;
	 indexRange[0]=TMath::Min(indexRange[0],pad);
	 indexRange[1]=TMath::Max(indexRange[1],pad);
	 indexRange[2]=TMath::Min(indexRange[2],time);
	 indexRange[3]=TMath::Max(indexRange[3],time); 
       }	 
    }
  } // end of loop over electrons
  
  return label; // returns track label when finished
}

//_____________________________________________________________________________
void AliTPC::GetList(Float_t label,Int_t np,TMatrix *m,Int_t *indexRange,
                     Float_t **pList)
{
  //----------------------------------------------------------------------
  //  Updates the list of tracks contributing to digits for a given row
  //----------------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  TMatrix &signal = *m;

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
  if (fHitType&2) branch = TH->GetBranch("TPC2");
  else branch = TH->GetBranch("TPC");

 
  //----------------------------------------------
  // Create TObjArray-s, one for each row,
  // each TObjArray will store the TVectors
  // of electrons, one TVector per each track.
  //---------------------------------------------- 
    
  Int_t *nofElectrons = new Int_t [nrows]; // electron counter for each row
  TVector **tracks = new TVector* [nrows]; //pointers to the track vectors
  for(i=0; i<nrows; i++){
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

    if (fHitType&2) {
      isInSector=kFALSE;
      TBranch * br = TH->GetBranch("fTrackHitsInfo");
      br->GetEvent(track);
      AliObjectArray * ar = fTrackHits->fTrackHitsInfo;
      for (UInt_t j=0;j<ar->GetSize();j++){
	if (  ((AliTrackHitsInfo*)ar->At(j))->fVolumeID==isec) isInSector=kTRUE;
      }
    }
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
      //      if(sector != isec) continue; 
      if(sector != isec){
	tpcHit = (AliTPChit*) NextHit();
	continue; 
      }

	currentTrack = tpcHit->Track(); // track number


        if(currentTrack != previousTrack){
                          
           // store already filled fTrack
              
	   for(i=0;i<nrows;i++){
             if(previousTrack != -1){
	       if(nofElectrons[i]>0){
	         TVector &v = *tracks[i];
		 v(0) = previousTrack;
                 tracks[i]->ResizeTo(4*nofElectrons[i]+1); // shrink if necessary
	         row[i]->Add(tracks[i]);                     
	       }
               else{
                 delete tracks[i]; // delete empty TVector
                 tracks[i]=0;
	       }
	     }

             nofElectrons[i]=0;
             tracks[i] = new TVector(481); // TVectors for the next fTrack

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
	  xyz[3]= (Float_t) (-gasgain*TMath::Log(gRandom->Rndm()));
	  index[0]=1;
	  
	  TransportElectron(xyz,index); //MI change -august	  
	  Int_t rowNumber;
	  fTPCParam->GetPadRow(xyz,index); //MI change august
	  rowNumber = index[2];
	  //transform position to local digit coordinates
	  //relative to nearest pad row 
	  if ((rowNumber<0)||rowNumber>=fTPCParam->GetNRow(isec)) continue;	  
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
	  
	  TVector &v = *tracks[rowNumber];
	  Int_t idx = 4*nofElectrons[rowNumber]-3;

	  v(idx)=  xyz[0];   // X - pad row coordinate
	  v(idx+1)=xyz[1];   // Y - pad coordinate (along the pad-row)
          v(idx+2)=xyz[2];   // Z - time bin coordinate
	  v(idx+3)=xyz[3];   // avalanche size  
	} // end of loop over electrons

        tpcHit = (AliTPChit*)NextHit();
        
      } // end of loop over hits
    } // end of loop over tracks

    //
    //   store remaining track (the last one) if not empty
    //

     for(i=0;i<nrows;i++){
       if(nofElectrons[i]>0){
          TVector &v = *tracks[i];
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
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" TPC_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

//_____________________________________________________________________________
void AliTPC::MakeBranch(Option_t* option)
{
  //
  // Create Tree branches for the TPC.
  //
  Int_t buffersize = 4000;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);

  char *d = strstr(option,"D");

  if (fDigits   && gAlice->TreeD() && d) {
    gAlice->TreeD()->Branch(branchname,&fDigits, buffersize);
    printf("Making Branch %s for digits\n",branchname);
  }	

  if (fHitType&2) MakeBranch2(option); // MI change 14.09.2000
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
    Float_t x1=xyz[0];
    fTPCParam->Transform2to2NearestWire(xyz,index);
    Float_t dx=xyz[0]-x1;
    xyz[1]+=dx*(fTPCParam->GetOmegaTau());
  }
  //add nonisochronity (not implemented yet)
  
}
//_____________________________________________________________________________
void AliTPC::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliTPC.
  //
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliDetector::Streamer(R__b);
      R__b >> fTPCParam;
      if (R__v < 2) return;
      R__b >> fNsectors;
      R__b >> fHitType;

   } else {
      R__b.WriteVersion(AliTPC::IsA());
      AliDetector::Streamer(R__b);
      R__b << fTPCParam;
      R__b << fNsectors;
      R__b << fHitType; 
   }
}
  
ClassImp(AliTPCdigit)
 
//_____________________________________________________________________________
AliTPCdigit::AliTPCdigit(Int_t *tracks, Int_t *digits):
  AliDigit(tracks)
{
  //
  // Creates a TPC digit object
  //
  fSector     = digits[0];
  fPadRow     = digits[1];
  fPad        = digits[2];
  fTime       = digits[3];
  fSignal     = digits[4];
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
// Additional code because of the AliTPCTrackHits

void AliTPC::MakeBranch2(Option_t *option)
{
  //
  // Create a new branch in the current Root Tree
  // The branch of fHits is automatically split
  // MI change 14.09.2000
  if (fHitType&2==0) return;
  char branchname[10];
  sprintf(branchname,"%s2",GetName());  
  //
  // Get the pointer to the header
  char *cH = strstr(option,"H");
  //
  if (fTrackHits   && gAlice->TreeH() && cH) {    
    //    gAlice->TreeH()->Branch(branchname,&fTrackHits, fBufferSize);
    AliObjectBranch * branch = new AliObjectBranch(branchname,"AliTPCTrackHits",&fTrackHits, 
						   gAlice->TreeH(),fBufferSize,1);
    gAlice->TreeH()->GetListOfBranches()->Add(branch);
    printf("* AliDetector::MakeBranch * Making Branch %s for trackhits\n",branchname);
  }	
}

void AliTPC::SetTreeAddress()
{
  if (fHitType&1) AliDetector::SetTreeAddress();
  if (fHitType&2) SetTreeAddress2();
}

void AliTPC::SetTreeAddress2()
{
  //
  // Set branch address for the TrackHits Tree
  // 
  TBranch *branch;
  char branchname[20];
  sprintf(branchname,"%s2",GetName());
  //
  // Branch address for hit tree
  TTree *treeH = gAlice->TreeH();
  if (treeH) {
    branch = treeH->GetBranch(branchname);
    if (branch) branch->SetAddress(&fTrackHits);
  }
}

void AliTPC::FinishPrimary()
{
  if (fTrackHits) fTrackHits->FlushHitStack();  
}


void AliTPC::AddHit2(Int_t track, Int_t *vol, Float_t *hits)
{ 
  //
  // add hit to the list  
  TClonesArray &particles = *(gAlice->Particles());
  Int_t rtrack;
  if (fIshunt) {
    int primary = gAlice->GetPrimary(track);
    ((TParticle *)particles[primary])->SetBit(kKeepBit);
    rtrack=primary;
  } else {
    rtrack=track;
    gAlice->FlagTrack(track);
  }  
  //AliTPChit *hit = (AliTPChit*)fHits->UncheckedAt(fNhits-1);
  //if (hit->fTrack!=rtrack)
  //  cout<<"bad track number\n";
  if (fTrackHits) 
    fTrackHits->AddHitKartez(vol[0],rtrack, hits[0],
			     hits[1],hits[2],(Int_t)hits[3]);
}

void AliTPC::ResetHits()
{
  if (fHitType&1) AliDetector::ResetHits();
  if (fHitType&2) ResetHits2();
}

void AliTPC::ResetHits2()
{
  //
  //reset hits
  if (fTrackHits) fTrackHits->Clear();
}   

AliHit* AliTPC::FirstHit(Int_t track)
{
  if (fHitType&2) return FirstHit2(track);
  return AliDetector::FirstHit(track);
}
AliHit* AliTPC::NextHit()
{
  if (fHitType&2) return NextHit2();
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
    gAlice->TreeH()->GetEvent(track);
  }
  //
  if (fTrackHits) {
    fTrackHits->First();
    return fTrackHits->GetHit();
  }
  else return 0;
}

AliHit* AliTPC::NextHit2()
{
  //
  //Return the next hit for the current track

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
  if (!fTrackHits) return;
  AliObjectArray * arr = fTrackHits->fTrackHitsInfo;
  for (UInt_t i=0;i<arr->GetSize();i++){
    AliTrackHitsInfo * info = (AliTrackHitsInfo *)(arr->At(i));
    info->fTrackID = map[info->fTrackID];
  }
  
}


//_____________________________________________________________________________
void AliTPC::LoadPoints2(Int_t)
{
  //
  // Store x, y, z of all hits in memory
  //
  if (fTrackHits == 0) return;
  //
  Int_t nhits = fTrackHits->GetEntriesFast();
  if (nhits == 0) return;
  Int_t tracks = gAlice->GetNtrack();
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
  //for (Int_t hit=0;hit<nhits;hit++) {
  while (ahit){
    //    ahit = (AliHit*)fHits->UncheckedAt(hit);
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
  Int_t tracks = gAlice->GetNtrack();
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



void AliTPC::FindTrackHitsIntersection(TClonesArray * arr)
{

  //
  //fill clones array with intersection of current point with the
  //middle of the row
  Int_t sector;
  Int_t ipart;
  
  const Int_t kcmaxhits=30000;
  TVector * xxxx = new TVector(kcmaxhits*4);
  TVector & xxx = *xxxx;
  Int_t maxhits = kcmaxhits;
      
  //
  AliTPChit * tpcHit=0;
  tpcHit = (AliTPChit*)FirstHit2(-1);
  Int_t currentIndex=0;
  Int_t lastrow=-1;  //last writen row

  while (tpcHit){
    if (tpcHit==0) continue;
    sector=tpcHit->fSector; // sector number
    ipart=tpcHit->Track();
    
    //find row number
    
    Float_t  x[3]={tpcHit->X(),tpcHit->Y(),tpcHit->Z()};
    Int_t    index[3]={1,sector,0};
    Int_t    currentrow = fTPCParam->GetPadRow(x,index) ;	
    if (currentrow<0) continue;
    if (lastrow<0) lastrow=currentrow;
    if (currentrow==lastrow){
      if ( currentIndex>=maxhits){
	maxhits+=kcmaxhits;
	xxx.ResizeTo(4*maxhits);
      }     
      xxx(currentIndex*4)=x[0];
      xxx(currentIndex*4+1)=x[1];
      xxx(currentIndex*4+2)=x[2];	
      xxx(currentIndex*4+3)=tpcHit->fQ;
      currentIndex++;	
    }
    else 
      if (currentIndex>2){
	Float_t sumx=0;
	Float_t sumx2=0;
	Float_t sumx3=0;
	Float_t sumx4=0;
	Float_t sumy=0;
	Float_t sumxy=0;
	Float_t sumx2y=0;
	Float_t sumz=0;
	Float_t sumxz=0;
	Float_t sumx2z=0;
	Float_t sumq=0;
	for (Int_t index=0;index<currentIndex;index++){
	  Float_t x,x2,x3,x4;
	  x=x2=x3=x4=xxx(index*4);
	  x2*=x;
	  x3*=x2;
	  x4*=x3;
	  sumx+=x;
	  sumx2+=x2;
	  sumx3+=x3;
	  sumx4+=x4;
	  sumy+=xxx(index*4+1);
	  sumxy+=xxx(index*4+1)*x;
	  sumx2y+=xxx(index*4+1)*x2;
	  sumz+=xxx(index*4+2);
	  sumxz+=xxx(index*4+2)*x;
	  sumx2z+=xxx(index*4+2)*x2;	 
	  sumq+=xxx(index*4+3);
	}
	Float_t centralPad = (fTPCParam->GetNPads(sector,lastrow)-1)/2;
	Float_t det=currentIndex*(sumx2*sumx4-sumx3*sumx3)-sumx*(sumx*sumx4-sumx2*sumx3)+
	  sumx2*(sumx*sumx3-sumx2*sumx2);
	
	Float_t detay=sumy*(sumx2*sumx4-sumx3*sumx3)-sumx*(sumxy*sumx4-sumx2y*sumx3)+
	  sumx2*(sumxy*sumx3-sumx2y*sumx2);
	Float_t detaz=sumz*(sumx2*sumx4-sumx3*sumx3)-sumx*(sumxz*sumx4-sumx2z*sumx3)+
	  sumx2*(sumxz*sumx3-sumx2z*sumx2);
	
	Float_t detby=currentIndex*(sumxy*sumx4-sumx2y*sumx3)-sumy*(sumx*sumx4-sumx2*sumx3)+
	  sumx2*(sumx*sumx2y-sumx2*sumxy);
	Float_t detbz=currentIndex*(sumxz*sumx4-sumx2z*sumx3)-sumz*(sumx*sumx4-sumx2*sumx3)+
	  sumx2*(sumx*sumx2z-sumx2*sumxz);
	
	Float_t y=detay/det+centralPad;
	Float_t z=detaz/det;	
	Float_t by=detby/det; //y angle
	Float_t bz=detbz/det; //z angle
	sumy/=Float_t(currentIndex);
	sumz/=Float_t(currentIndex);
	
	AliComplexCluster cl;
	cl.fX=z;
	cl.fY=y;
	cl.fQ=sumq;
	cl.fSigmaX2=bz;
	cl.fSigmaY2=by;
	cl.fTracks[0]=ipart;
	
	AliTPCClustersRow * row = (fClustersArray->GetRow(sector,lastrow));
	if (row!=0) row->InsertCluster(&cl);
	currentIndex=0;
	lastrow=currentrow;
      } //end of calculating cluster for given row
        	
  } // end of loop over hits
  xxxx->Delete();
 


}
