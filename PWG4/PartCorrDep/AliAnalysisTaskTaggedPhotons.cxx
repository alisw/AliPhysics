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

/* $Id: */

//_________________________________________________________________________
// Analysis for Tagged Photons
// Prepares all necessary histograms for calculation of 
// the yield of pi0 decay photon in calorimeter:
// Marks photons which makes pi0 with some other and
// fill invariant mass distributions for estimate background below pi0
// peak so that estimate proportion of fake pairs. 
// Fills as well controll MC histograms with clasification of the photon origin 
// and check of the proportion of truly tagged photons.
// 
//
//*-- Dmitry Blau 
//////////////////////////////////////////////////////////////////////////////

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TROOT.h>

#include "AliAnalysisTaskTaggedPhotons.h" 
#include "AliAnalysisManager.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h" 
#include "AliAODEvent.h" 
#include "AliESDCaloCluster.h" 
#include "AliAODPWG4Particle.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "TGeoManager.h"
#include "AliMCAnalysisUtils.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeoUtils.h"



//______________________________________________________________________________
AliAnalysisTaskTaggedPhotons::AliAnalysisTaskTaggedPhotons() : AliAnalysisTaskSE(),
  fPHOSgeom(0x0),
  fEMCALgeom(0x0),
  fStack(0x0),fPHOS(1),
  fPhotonId(1.0),fMinEnergyCut(0.4),
  fPi0MeanP0(0.136),fPi0MeanP1(0.0),fPi0MeanP2(0.0),fPi0MeanP3(0.0),
  fPi0SigmaP0(0.004508),fPi0SigmaP1(0.005497),fPi0SigmaP2(0.00000006),
  fZmax(0.),fZmin(0.),fPhimax(0.),fPhimin(0.),
  fOutputList(0x0),fEventList(0x0),
  fhRecAllArea1(0x0),fhRecAllArea2(0x0),fhRecAllArea3(0x0),
  fhRecPhoton(0x0),fhRecOther(0x0),fhRecPhotPi0(0x0),fhRecPhotEta(0x0),
  fhRecPhotOmega(0x0),fhRecPhotEtapr(0x0),fhRecPhotConv(0x0),fhRecPhotHadron(0x0),
  fhRecPhotDirect(0x0),fhRecPhotOther(0x0),
  fhDecWMCPartner(0x0),fhDecWMissedPartnerNotPhoton(0x0),fhDecWMissedPartnerAll(0x0),
  fhDecWMissedPartnerEmin(0x0),fhDecWMissedPartnerConv(0x0),fhDecWMissedPartnerStack(0x0),
  fhDecWMissedPartnerGeom0(0x0),
  fhDecWMissedPartnerGeom1(0x0),fhDecWMissedPartnerGeom2(0x0),fhDecWMissedPartnerGeom3(0x0),
  fhPartnerMCReg(0x0),fhPartnerMissedEmin(0x0),fhPartnerMissedConv(0x0),fhPartnerMissedGeo(0x0),
  fhTaggedAll(0x0),fhTaggedArea1(0x0),fhTaggedArea2(0x0),fhTaggedArea3(0x0),fhTaggedMult(0x0),
  fhTaggedMCTrue(0x0),fhMCMissedTagging(0x0),fhMCFakeTagged(0x0),
  fhMCMissedTaggingMass(0x0),
  fhConversionRadius(0x0),fhInteractionRadius(0x0),fhEvents(0x0)
{
  //Default constructor
  for (Int_t i = 0; i < 4; i++)
    {
      fhRecAll[i] = 0;
      fhRecPhotonPID[i] = 0;
      fhRecOtherPID[i] = 0;
      fhTaggedPID[i] = 0;
      fhInvMassReal[i] = 0;
      fhInvMassMixed[i] = 0;
    }

}
//______________________________________________________________________________
AliAnalysisTaskTaggedPhotons::AliAnalysisTaskTaggedPhotons(const char *name) : 
  AliAnalysisTaskSE(name),
  fPHOSgeom(0x0),
  fEMCALgeom(0x0),
  fStack(0x0),fPHOS(1),
  fPhotonId(1.0),fMinEnergyCut(0.4),
  fPi0MeanP0(0.136),fPi0MeanP1(0.0),fPi0MeanP2(0.0),fPi0MeanP3(0.0),
  fPi0SigmaP0(0.004508),fPi0SigmaP1(0.005497),fPi0SigmaP2(0.00000006),
  fZmax(0.),fZmin(0.),fPhimax(0.),fPhimin(0.),
  fOutputList(0x0),fEventList(0x0),
  fhRecAllArea1(0x0),fhRecAllArea2(0x0),fhRecAllArea3(0x0),
  fhRecPhoton(0x0),fhRecOther(0x0),fhRecPhotPi0(0x0),fhRecPhotEta(0x0),
  fhRecPhotOmega(0x0),fhRecPhotEtapr(0x0),fhRecPhotConv(0x0),fhRecPhotHadron(0x0),
  fhRecPhotDirect(0x0),fhRecPhotOther(0x0),
  fhDecWMCPartner(0x0),fhDecWMissedPartnerNotPhoton(0x0),fhDecWMissedPartnerAll(0x0),
  fhDecWMissedPartnerEmin(0x0),fhDecWMissedPartnerConv(0x0),fhDecWMissedPartnerStack(0x0),
  fhDecWMissedPartnerGeom0(0x0),
  fhDecWMissedPartnerGeom1(0x0),fhDecWMissedPartnerGeom2(0x0),fhDecWMissedPartnerGeom3(0x0),
  fhPartnerMCReg(0x0),fhPartnerMissedEmin(0x0),fhPartnerMissedConv(0x0),fhPartnerMissedGeo(0x0),
  fhTaggedAll(0x0),fhTaggedArea1(0x0),fhTaggedArea2(0x0),fhTaggedArea3(0x0),fhTaggedMult(0x0),
  fhTaggedMCTrue(0x0),fhMCMissedTagging(0x0),fhMCFakeTagged(0x0),
  fhMCMissedTaggingMass(0x0),
  fhConversionRadius(0x0),fhInteractionRadius(0x0),fhEvents(0x0)
{
  // Constructor.
  for (Int_t i = 0; i < 4; i++)
    {
      fhRecAll[i] = 0;
      fhRecPhotonPID[i] = 0;
      fhRecOtherPID[i] = 0;
      fhTaggedPID[i] = 0;
      fhInvMassReal[i] = 0;
      fhInvMassMixed[i] = 0;
    }

  // Output slots 
  DefineOutput(1,  TList::Class()) ; 
}

//____________________________________________________________________________
AliAnalysisTaskTaggedPhotons::AliAnalysisTaskTaggedPhotons(const AliAnalysisTaskTaggedPhotons& ap) :
  AliAnalysisTaskSE(ap.GetName()),  
  fPHOSgeom(0x0),
  fEMCALgeom(0x0),
  fStack(0x0),fPHOS(ap.fPHOS),
  fPhotonId(ap.fPhotonId),fMinEnergyCut(ap.fMinEnergyCut),
  fPi0MeanP0(ap.fPi0MeanP0),fPi0MeanP1(ap.fPi0MeanP1),fPi0MeanP2(ap.fPi0MeanP2),fPi0MeanP3(ap.fPi0MeanP3),
  fPi0SigmaP0(ap.fPi0SigmaP0),fPi0SigmaP1(ap.fPi0SigmaP1),fPi0SigmaP2(ap.fPi0SigmaP2),
  fZmax(ap.fZmax),fZmin(ap.fZmin),fPhimax(ap.fPhimax),fPhimin(ap.fPhimin),
  fOutputList(0x0),fEventList(0x0),
  fhRecAllArea1(0x0),fhRecAllArea2(0x0),fhRecAllArea3(0x0),
  fhRecPhoton(0x0),fhRecOther(0x0),fhRecPhotPi0(0x0),fhRecPhotEta(0x0),
  fhRecPhotOmega(0x0),fhRecPhotEtapr(0x0),fhRecPhotConv(0x0),fhRecPhotHadron(0x0),
  fhRecPhotDirect(0x0),fhRecPhotOther(0x0),
  fhDecWMCPartner(0x0),fhDecWMissedPartnerNotPhoton(0x0),fhDecWMissedPartnerAll(0x0),
  fhDecWMissedPartnerEmin(0x0),fhDecWMissedPartnerConv(0x0),fhDecWMissedPartnerStack(0x0),
  fhDecWMissedPartnerGeom0(0x0),
  fhDecWMissedPartnerGeom1(0x0),fhDecWMissedPartnerGeom2(0x0),fhDecWMissedPartnerGeom3(0x0),
  fhPartnerMCReg(0x0),fhPartnerMissedEmin(0x0),fhPartnerMissedConv(0x0),
  fhPartnerMissedGeo(0x0),
  fhTaggedAll(0x0),fhTaggedArea1(0x0),fhTaggedArea2(0x0),fhTaggedArea3(0x0),fhTaggedMult(0x0),
  fhTaggedMCTrue(0x0),fhMCMissedTagging(0x0),fhMCFakeTagged(0x0),
  fhMCMissedTaggingMass(0x0),
  fhConversionRadius(0x0),fhInteractionRadius(0x0),fhEvents(0x0)
{
  // cpy ctor
  for (Int_t i = 0; i < 4; i++)
    {
      fhRecAll[i] = 0;
      fhRecPhotonPID[i] = 0;
      fhRecOtherPID[i] = 0;
      fhTaggedPID[i] = 0;
      fhInvMassReal[i] = 0;
      fhInvMassMixed[i] = 0;
    }

}

//_____________________________________________________________________________
AliAnalysisTaskTaggedPhotons& AliAnalysisTaskTaggedPhotons::operator = (const AliAnalysisTaskTaggedPhotons& ap)
{
// assignment operator

  this->~AliAnalysisTaskTaggedPhotons();
  new(this) AliAnalysisTaskTaggedPhotons(ap);
  return *this;
}

//______________________________________________________________________________
AliAnalysisTaskTaggedPhotons::~AliAnalysisTaskTaggedPhotons()
{
  // dtor
  if(fOutputList) {
    fOutputList->Clear() ; 
    delete fOutputList ;
  }
}


//________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::UserCreateOutputObjects()
{ 


  //Load geometry
  //if file "geometry.root" exists, load it
  //otherwise use misaligned matrixes stored in ESD

  // Create the outputs containers

  OpenFile(1) ; 
  const Int_t nPtBins=52 ;
  Double_t ptBins[nPtBins+1] ;
  for(Int_t i=0;i<=20;i++)ptBins[i]=0.1*i ; //0-2 GeV:  0.1 GeV/bin
  for(Int_t i=21;i<=30;i++)ptBins[i]=2.+0.2*(i-20) ; //2-4 GeV:  0.2 GeV/bin                   
  for(Int_t i=31;i<=38;i++)ptBins[i]=4.+0.5*(i-30) ; //4-8 GeV:  0.5 GeV/bin                   
  for(Int_t i=39;i<=42;i++)ptBins[i]=8.+1.0*(i-38) ; //8-12 GeV: 1. GeV/bin                   
  for(Int_t i=43;i<=52;i++)ptBins[i]=12.+2.0*(i-42) ; //12-30 GeV: 2. GeV/bin                   


  //Reconstructed spectra
    fhRecAll[0]   = new TH1D("fhRecAll0","Spectrum of all reconstructed particles, no PID", nPtBins, ptBins ) ;
    fhRecAll[1]   = new TH1D("fhRecAll1","Spectrum of all reconstructed particles, PID=1", nPtBins, ptBins ) ;
    fhRecAll[2]   = new TH1D("fhRecAll2","Spectrum of all reconstructed particles, PID=2", nPtBins, ptBins ) ;
    fhRecAll[3]   = new TH1D("fhRecAll3","Spectrum of all reconstructed particles, PID=3", nPtBins, ptBins ) ;

    fhRecAllArea1 = new TH1D("fhRecAllArea1","Spectrum of rec particles in Fid. Area 1", nPtBins, ptBins ) ;
    fhRecAllArea2 = new TH1D("fhRecAllArea2","Spectrum of rec particles in Fid. Area 2", nPtBins, ptBins ) ;
    fhRecAllArea3 = new TH1D("fhRecAllArea3","Spectrum of rec particles in Fid. Area 3", nPtBins, ptBins ) ;

    //Sort registered particles spectra according MC information
    fhRecPhoton   = new TH1D("fhRecPhoton","Spectrum of rec. with primary==22 and no PID criteria", nPtBins, ptBins ) ;
    fhRecOther    = new TH1D("fhRecOther"," Spectrum of rec. with primary!=22 and no PID criteria", nPtBins, ptBins ) ;
    fhRecPhotonPID[0]= new TH1D("fhRecPhotonPID0","Spectrum of rec. with primary==22 and PID=0", nPtBins, ptBins ) ;
    fhRecPhotonPID[1]= new TH1D("fhRecPhotonPID1","Spectrum of rec. with primary==22 and PID=1", nPtBins, ptBins ) ;
    fhRecPhotonPID[2]= new TH1D("fhRecPhotonPID2","Spectrum of rec. with primary==22 and PID=2", nPtBins, ptBins ) ;
    fhRecPhotonPID[3]= new TH1D("fhRecPhotonPID3","Spectrum of rec. with primary==22 and PID=3", nPtBins, ptBins ) ;
    fhRecOtherPID[0]= new TH1D("fhRecOtherPID0","Spectrum of rec. with primary!=22 and PID=0", nPtBins, ptBins ) ;
    fhRecOtherPID[1]= new TH1D("fhRecOtherPID1","Spectrum of rec. with primary!=22 and PID=1", nPtBins, ptBins ) ;
    fhRecOtherPID[2]= new TH1D("fhRecOtherPID2","Spectrum of rec. with primary!=22 and PID=2", nPtBins, ptBins ) ;
    fhRecOtherPID[3]= new TH1D("fhRecOtherPID3","Spectrum of rec. with primary!=22 and PID=3", nPtBins, ptBins ) ;
    fhRecPhotPi0    = new TH1D("fhRecPhotPi0","Spectrum of rec. photons from pi0 decays", nPtBins, ptBins ) ;
    fhRecPhotEta    = new TH1D("fhRecPhotEta","Spectrum of rec. photons from eta decays", nPtBins, ptBins ) ;
    fhRecPhotOmega  = new TH1D("fhRecPhotOmega","Spectrum of rec. photons from omega decays", nPtBins, ptBins ) ;
    fhRecPhotEtapr  = new TH1D("fhRecPhotEtapr","Spectrum of rec. photons from eta prime decays", nPtBins, ptBins ) ;
    fhRecPhotConv   = new TH1D("fhRecPhotConv"," Spectrum of rec. photons from conversion", nPtBins, ptBins ) ;
    fhRecPhotHadron = new TH1D("fhRecPhotHadron","Spectrum of rec. photons from hadron-matter interactions", nPtBins, ptBins ) ;
    fhRecPhotDirect = new TH1D("fhRecPhotDirect","Spectrum of rec. photons direct or no primary", nPtBins, ptBins ) ;
    fhRecPhotOther  = new TH1D("fhRecPhotOther","Spectrum of rec. photons from other hadron decays", nPtBins, ptBins ) ;

    //MC tagging: reasons of partner loss etc.
    fhDecWMCPartner        = new TH1D("fhDecWMCPartner","pi0 decay photon which partner should be registered according to MC", nPtBins, ptBins ) ;
    fhDecWMissedPartnerNotPhoton = new TH1D("fhDecWMissedPartnerNotPhoton","Decay photon with missed non-photon partner", nPtBins, ptBins ) ;
    fhDecWMissedPartnerAll  = new TH1D("fhDecWMissedPartnerAll","Decay photons with partner missed due to some reason", nPtBins, ptBins ) ;
    fhDecWMissedPartnerEmin = new TH1D("fhDecWMissedPartnerEmin","Decay photons with partner missed due to low energy", nPtBins, ptBins ) ;
    fhDecWMissedPartnerConv = new TH1D("fhDecWMissedPartnerConv","Decay photons with partner missed due to conversion", nPtBins, ptBins ) ;
    fhDecWMissedPartnerStack= new TH1D("fhDecWMissedPartnerStack","Decay photons with partner not in Stack", nPtBins, ptBins ) ;
    fhDecWMissedPartnerGeom0 = new TH1D("fhDecWMissedPartnerGeom0","Decay photons with partner missed due geometry", nPtBins, ptBins ) ;
    fhDecWMissedPartnerGeom1 = new TH1D("fhDecWMissedPartnerGeom1","Decay photons with partner missed due geometry Fid. area. 1", nPtBins, ptBins ) ;
    fhDecWMissedPartnerGeom2 = new TH1D("fhDecWMissedPartnerGeom2","Decay photons with partner missed due geometry Fid. area. 2", nPtBins, ptBins ) ;
    fhDecWMissedPartnerGeom3 = new TH1D("fhDecWMissedPartnerGeom3","Decay photons with partner missed due geometry Fid. area. 3", nPtBins, ptBins ) ;

    //MC tagging: Decay partners spectra
    fhPartnerMCReg      = new TH1D("fhPartnerMCReg","Spectrum of decay partners which should be registered (MC)", nPtBins, ptBins ) ;
    fhPartnerMissedEmin = new TH1D("fhPartnerMissedEmin","Spectrum of decay partners lost due to Emin cut", nPtBins, ptBins ) ;
    fhPartnerMissedConv = new TH1D("fhPartnerMissedConv","Spectrum of decay partners lost due to conversion", nPtBins, ptBins ) ;
    fhPartnerMissedGeo  = new TH1D("fhPartnerMissedGeo","Spectrum of decay partners lost due to acceptance", nPtBins, ptBins ) ;

    //Tagging
    fhTaggedAll    = new TH1D("fhTaggedAll","Spectrum of all tagged photons", nPtBins, ptBins ) ;
    fhTaggedArea1  = new TH1D("fhTaggedArea1","Spectrum of all tagged photons Fid. area1", nPtBins, ptBins ) ;
    fhTaggedArea2  = new TH1D("fhTaggedArea2","Spectrum of all tagged photons Fid. area2", nPtBins, ptBins ) ;
    fhTaggedArea3  = new TH1D("fhTaggedArea3","Spectrum of all tagged photons Fid. area3", nPtBins, ptBins ) ;
    fhTaggedPID[0] = new TH1D("fhTaggedPID0","Spectrum of tagged photons for PID=0", nPtBins, ptBins ) ;
    fhTaggedPID[1] = new TH1D("fhTaggedPID1","Spectrum of tagged photons for PID=1", nPtBins, ptBins ) ;
    fhTaggedPID[2] = new TH1D("fhTaggedPID2","Spectrum of tagged photons for PID=2", nPtBins, ptBins ) ;
    fhTaggedPID[3] = new TH1D("fhTaggedPID3","Spectrum of tagged photons for PID=3", nPtBins, ptBins ) ;
    fhTaggedMult   = new TH1D("fhTaggedMult","Spectrum of multiply tagged photons", nPtBins, ptBins ) ;

    //Tagging: use MC information if available
    fhTaggedMCTrue    = new TH1D("fhTaggedMCTrue","Spectrum of correctly tagged pi0 decay photons", nPtBins, ptBins ) ;
    fhMCMissedTagging = new TH1D("fhMCMissedTagging","Spectrum of pi0 decay photons missed tagging due to wrong pair mass", nPtBins, ptBins ) ;
    fhMCFakeTagged    = new TH1D("fhMCFakeTagged","Spectrum of photons wrongly tagged according to MC", nPtBins, ptBins ) ;

    //Invariant mass distributions for fake corrections
    const Int_t nmass=500 ;
    Double_t masses[nmass+1] ;
    for(Int_t i=0;i<=nmass;i++)masses[i]=0.002*i ;
    fhInvMassReal[0] = new TH2D("hInvMassRealPID0","Two-photon inv. mass vs first photon pt",nmass,masses,nPtBins, ptBins ) ;
    fhInvMassReal[1] = new TH2D("hInvMassRealPID1","Two-photon inv. mass vs first photon pt",nmass,masses,nPtBins, ptBins ) ;
    fhInvMassReal[2] = new TH2D("hInvMassRealPID2","Two-photon inv. mass vs first photon pt",nmass,masses,nPtBins, ptBins ) ;
    fhInvMassReal[3] = new TH2D("hInvMassRealPID3","Two-photon inv. mass vs first photon pt",nmass,masses,nPtBins, ptBins ) ;
    fhInvMassMixed[0]  = new TH2D("hInvMassMixedPID0","Two-photon inv. mass vs first photon pt",nmass,masses,nPtBins, ptBins ) ;
    fhInvMassMixed[1]  = new TH2D("hInvMassMixedPID1","Two-photon inv. mass vs first photon pt",nmass,masses,nPtBins, ptBins ) ;
    fhInvMassMixed[2]  = new TH2D("hInvMassMixedPID2","Two-photon inv. mass vs first photon pt",nmass,masses,nPtBins, ptBins ) ;
    fhInvMassMixed[3]  = new TH2D("hInvMassMixedPID3","Two-photon inv. mass vs first photon pt",nmass,masses,nPtBins, ptBins ) ;
    fhMCMissedTaggingMass= new TH2D("hMCMissedTaggingMass","Inv mass of pairs missed tagging",nmass,masses,nPtBins, ptBins ) ;

    //Conversion and annihilation radius distributions
    fhConversionRadius  = new TH1D("fhConversionRadius","Radii of photon production (conversion)",100,0.,500.) ;
    fhInteractionRadius = new TH1D("fhInteractionRadius","Radii of photon production (hadron interaction)",100,0.,500.);

    fhEvents               = new TH1D("hEvents", "Number of Events processed", 1, 0, 1);

// create output container
  
  fOutputList = new TList() ;
  fEventList = new TList() ; 
  fOutputList->SetName(GetName()) ; 

  fOutputList->AddAt(fhRecAll[0],                            0) ;
  fOutputList->AddAt(fhRecAll[1],                            1) ;
  fOutputList->AddAt(fhRecAll[2],                            2) ;
  fOutputList->AddAt(fhRecAll[3],                            3) ;
  fOutputList->AddAt(fhRecAllArea1,                          4) ;
  fOutputList->AddAt(fhRecAllArea2,                          5) ;
  fOutputList->AddAt(fhRecAllArea3,                          6) ;

  fOutputList->AddAt(fhRecPhoton,                            7) ;
  fOutputList->AddAt(fhRecOther,                             8) ;
  fOutputList->AddAt(fhRecPhotonPID[0],                      9) ;
  fOutputList->AddAt(fhRecPhotonPID[1],                     10) ;
  fOutputList->AddAt(fhRecPhotonPID[2],                     11) ;
  fOutputList->AddAt(fhRecPhotonPID[3],                     12) ;
  fOutputList->AddAt(fhRecOtherPID[0],                      13) ;
  fOutputList->AddAt(fhRecOtherPID[1],                      14) ;
  fOutputList->AddAt(fhRecOtherPID[2],                      15) ;
  fOutputList->AddAt(fhRecOtherPID[3],                      16) ;
  fOutputList->AddAt(fhRecPhotPi0,                          17) ;
  fOutputList->AddAt(fhRecPhotEta,                          18) ;
  fOutputList->AddAt(fhRecPhotOmega,                        19) ;
  fOutputList->AddAt(fhRecPhotEtapr,                        20) ;
  fOutputList->AddAt(fhRecPhotConv,                         21) ;
  fOutputList->AddAt(fhRecPhotHadron,                       22) ;
  fOutputList->AddAt(fhRecPhotDirect,                       23) ;
  fOutputList->AddAt(fhRecPhotOther,                        24) ;

  fOutputList->AddAt(fhDecWMCPartner,                       25) ;
  fOutputList->AddAt(fhDecWMissedPartnerNotPhoton,          26) ;
  fOutputList->AddAt(fhDecWMissedPartnerAll,                27) ;
  fOutputList->AddAt(fhDecWMissedPartnerEmin,               28) ;
  fOutputList->AddAt(fhDecWMissedPartnerConv,               29) ;
  fOutputList->AddAt(fhDecWMissedPartnerStack,              30) ;
  fOutputList->AddAt(fhDecWMissedPartnerGeom0,              31) ;
  fOutputList->AddAt(fhDecWMissedPartnerGeom1,              32) ;
  fOutputList->AddAt(fhDecWMissedPartnerGeom2,              33) ;
  fOutputList->AddAt(fhDecWMissedPartnerGeom3,              34) ;

  fOutputList->AddAt(fhPartnerMCReg,                        35) ;
  fOutputList->AddAt(fhPartnerMissedEmin,                   36) ;
  fOutputList->AddAt(fhPartnerMissedConv,                   37) ;
  fOutputList->AddAt(fhPartnerMissedGeo,                    38) ;

  fOutputList->AddAt(fhTaggedAll,                           39) ;
  fOutputList->AddAt(fhTaggedArea1,                         40) ;
  fOutputList->AddAt(fhTaggedArea2,                         41) ;
  fOutputList->AddAt(fhTaggedArea3,                         42) ;
  fOutputList->AddAt(fhTaggedPID[0],                        43) ;
  fOutputList->AddAt(fhTaggedPID[1],                        44) ;
  fOutputList->AddAt(fhTaggedPID[2],                        45) ;
  fOutputList->AddAt(fhTaggedPID[3],                        46) ;
  fOutputList->AddAt(fhTaggedMult,                          47) ;

  fOutputList->AddAt(fhTaggedMCTrue,                        48) ;
  fOutputList->AddAt(fhMCMissedTagging,                     49) ;
  fOutputList->AddAt(fhMCFakeTagged,                        50) ;

  fOutputList->AddAt(fhInvMassReal[0],                      51) ;
  fOutputList->AddAt(fhInvMassReal[1],                      52) ;
  fOutputList->AddAt(fhInvMassReal[2],                      53) ;
  fOutputList->AddAt(fhInvMassReal[3],                      54) ;
  fOutputList->AddAt(fhInvMassMixed[0],                     55) ;
  fOutputList->AddAt(fhInvMassMixed[1],                     56) ;
  fOutputList->AddAt(fhInvMassMixed[2],                     57) ;
  fOutputList->AddAt(fhInvMassMixed[3],                     58) ;
  fOutputList->AddAt(fhMCMissedTaggingMass,                 59) ;

  fOutputList->AddAt(fhConversionRadius,                    60) ;
  fOutputList->AddAt(fhInteractionRadius,                   61) ;

  fOutputList->AddAt(fhEvents,                              62) ;

}

//______________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::UserExec(Option_t *) 
{
  //Fill all histograms


  // We prepared 4 PID histograms, choose which PID combinations to use
  // see AliAODPWG4Particle for more choises
  //    0=all,2=Disp,4=Ch,6=Ch+Disp
  const Int_t pidMap[4]={0,2,4,6} ;

  // Processing of one event
  if(fDebug>1)
    AliInfo(Form("\n\n Processing event # %lld",  Entry())) ; 
  AliESDEvent* esd = (AliESDEvent*)InputEvent();


  //read geometry if not read yet
  if((fPHOS && fPHOSgeom==0) || (!fPHOS && fEMCALgeom==0))
    InitGeometry() ;
  if(fPHOS && fPHOSgeom!=0){
  }
  if(!fPHOS && fEMCALgeom!=0){ //EMCAL geometry initialized, but still need to set matrixes
    AliAODEvent * aod = 0x0 ;
    if(!esd)
      aod=dynamic_cast<AliAODEvent*>(InputEvent()) ;
    if(!esd && !aod)
      AliFatal("Can not read geometry matrixes from ESD/AOD: NO ESD") ;
    for(Int_t mod=0; mod < (fEMCALgeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
      if(esd){
        const TGeoHMatrix* m=esd->GetEMCALMatrix(mod) ;
        fEMCALgeom->SetMisalMatrix(m, mod) ;
      }
      else{
        const TGeoHMatrix* m=aod->GetHeader()->GetEMCALMatrix(mod) ;
        fEMCALgeom->SetMisalMatrix(m, mod) ;
      }
    }
  }
  //MC stack init
  fStack=0x0 ;
  if(AliAnalysisManager::GetAnalysisManager()){
    AliMCEventHandler* mctruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if(mctruth)
      fStack = mctruth->MCEvent()->Stack();
  }

  if(!fStack && gDebug>1)
    AliInfo("No stack! \n");


//  Here one can choose trigger. But according to framework it should be chosen 
//  by external filters???
//  TString trigClasses = esd->GetFiredTriggerClasses();
//  if (!fStack && !trigClasses.Contains("CINT1B-ABCE-NOPF-ALL") ){
//    Printf("Skip event with trigger class \"%s\"",trigClasses.Data());
//    return;
//  }

  fhEvents->Fill(0.);

  //************************  PHOS/EMCAL *************************************
  TRefArray * caloClustersArr  = new TRefArray();  
  if(fPHOS)
    esd->GetPHOSClusters(caloClustersArr);
  else
    esd->GetEMCALClusters(caloClustersArr);
  const Int_t kNumberOfClusters = caloClustersArr->GetEntries() ;  

  TClonesArray * fCaloPhotonsArr   = new TClonesArray("AliAODPWG4Particle",kNumberOfClusters);
  Int_t inList = 0; //counter of caloClusters

  Int_t cluster ; 

  // loop over Clusters
  for(cluster = 0 ; cluster < kNumberOfClusters ; cluster++) {
    AliESDCaloCluster * caloCluster = static_cast<AliESDCaloCluster*>(caloClustersArr->At(cluster)) ;
  
    if((fPHOS && !caloCluster->IsPHOS()) ||
       (!fPHOS && caloCluster->IsPHOS()))
      continue ; 

    Double_t v[3] ; //vertex ;
    esd->GetVertex()->GetXYZ(v) ;

    TLorentzVector momentum ;
    caloCluster->GetMomentum(momentum, v);
    new ((*fCaloPhotonsArr)[inList]) AliAODPWG4Particle(momentum.Px(),momentum.Py(),momentum.Pz(),caloCluster->E() );
    AliAODPWG4Particle *p = static_cast<AliAODPWG4Particle*>(fCaloPhotonsArr->At(inList));
    inList++;

    p->SetCaloLabel(cluster,-1); //This and partner cluster
    p->SetDistToBad(Int_t(caloCluster->GetDistanceToBadChannel()));

    p->SetTag(AliMCAnalysisUtils::kMCUnknown);
    p->SetTagged(kFALSE);   //Reconstructed pairs found
    p->SetLabel(caloCluster->GetLabel());
    Float_t pos[3] ;
    caloCluster->GetPosition(pos) ;
    p->SetFiducialArea(GetFiducialArea(pos)) ;

    //PID criteria
    p->SetDispBit(TestDisp(caloCluster->GetM02(),caloCluster->GetM20(),caloCluster->E())) ;
    p->SetTOFBit(TestTOF(caloCluster->GetTOF(),caloCluster->E())) ;
    p->SetChargedBit(TestCharged(caloCluster->GetEmcCpvDistance(),caloCluster->E())) ;

    for(Int_t iPID=0; iPID<4; iPID++){
      if(p->IsPIDOK(pidMap[iPID],22))
        fhRecAll[iPID]->Fill( p->Pt() ) ; //All recontructed particles
    }
    Int_t iFidArea = p->GetFiducialArea(); 
    if(iFidArea>0){
      fhRecAllArea1->Fill(p->Pt() ) ; 
      if(iFidArea>1){
        fhRecAllArea2->Fill(p->Pt() ) ; 
        if(iFidArea>2){
          fhRecAllArea3->Fill(p->Pt() ) ; 
        }
      }
    }

    if(fStack){
      TParticle * prim = fStack->Particle(caloCluster->GetLabel()) ;
      if(fDebug>2)
        printf("Pdgcode = %d\n",prim->GetPdgCode());

      if(prim->GetPdgCode()!=22){ //not photon
//<--DP        p->SetPhoton(kFALSE);
        fhRecOther->Fill(p->Pt()); //not photons real spectra
        for(Int_t iPID=0; iPID<4; iPID++){
          if(p->IsPIDOK(pidMap[iPID],22))
            fhRecOtherPID[iPID]->Fill(p->Pt());
        }
      }
      else{ //Primary photon (as in MC)
//<--DP        p->SetPhoton(kTRUE);
        fhRecPhoton->Fill(p->Pt()); //Reconstructed with primary photon
        for(Int_t iPID=0; iPID<4; iPID++){
          if(p->IsPIDOK(pidMap[iPID],22))
            fhRecPhotonPID[iPID]->Fill(p->Pt());
        }
        Int_t pi0i=prim->GetFirstMother();
        Int_t grandpaPDG=-1 ;
        TParticle * pi0p = 0;
        if(pi0i>=0){
          pi0p=fStack->Particle(pi0i);
          grandpaPDG=pi0p->GetPdgCode() ;
        }
        switch(grandpaPDG){
        case 111: //Pi0 decay
                //Primary decay photon (as in MC)
                fhRecPhotPi0->Fill(p->Pt());
                break ;
        case  11:
        case -11: //electron/positron conversion
//<--DP                p->SetConverted(1);
                fhRecPhotConv->Fill(p->Pt());  //Reconstructed with photon from conversion primary
                fhConversionRadius->Fill(prim->R());
                break ;
        case -2212:
        case -2112: //antineutron & antiproton conversion
                fhRecPhotHadron->Fill(p->Pt());  //Reconstructed with photon from antibaryon annihilation
                fhInteractionRadius->Fill(prim->R());
                break ;
          
        case 221: //eta decay
                fhRecPhotEta->Fill(p->Pt());
                break ;  
          
        case 223: //omega meson decay
                fhRecPhotOmega->Fill(p->Pt());
                break ;
            
        case 331: //eta' decay
                fhRecPhotEtapr->Fill(p->Pt());
                break ;
              
        case -1: //direct photon or no primary
                fhRecPhotDirect->Fill(p->Pt());
                break ;
              
        default:  
                fhRecPhotOther->Fill(p->Pt());
                break ;
        }  

        //Now classify pi0 decay photon
        if(grandpaPDG==111){
//<--DP          p->Pi0Decay(kTRUE); //Mark this photon as primary decayed
//<--DP          p->Pi0Id(pi0i); //remember id of the parent pi0

          //Now check if second (partner) photon from pi0 decay hits PHOS or not
          //i.e. both photons can be tagged or it's the systematic error
//<--DP          p->SetPartnerPt(0.); 
          Int_t indexdecay1,indexdecay2;

          indexdecay1=pi0p->GetFirstDaughter();
          indexdecay2=pi0p->GetLastDaughter();
          Int_t indexdecay=-1;
          if(fDebug>2)
                printf("checking pi0 decay...index1=%d, index2=%d, index_pi0=%d, index_ph_prim=%d\n", indexdecay1,indexdecay2,pi0i,caloCluster->GetLabel());
                
          if(indexdecay1!=caloCluster->GetLabel()) 
            indexdecay=indexdecay1;
          if(indexdecay2!=caloCluster->GetLabel()) 
            indexdecay=indexdecay2;
          if(indexdecay==-1){
            if(fDebug>2){
              printf("Probably the other photon is not in the stack!\n");
              printf("Number of daughters: %d\n",pi0p->GetNDaughters());
            }
            fhDecWMissedPartnerStack->Fill(p->Pt()) ;
          }  
          else{
            TParticle *partner = fStack->Particle(indexdecay);
//<--DP            p->SetPartnerPt(partner->Pt());
            if(partner->GetPdgCode()==22){ 
              Bool_t isPartnerLost=kFALSE; //If partner is lost for some reason
              if(partner->GetNDaughters()!=0){ //this photon is converted before it is registered by some detector
                if(fDebug>2)
                  printf("P_Conv, daughters=%d\n",partner->GetNDaughters());
//<--DP                p->SetConvertedPartner(1);
                fhPartnerMissedConv->Fill(partner->Pt());
                fhDecWMissedPartnerConv->Fill(p->Pt()) ;  //Spectrum of tagged with missed partner
                isPartnerLost=kTRUE;
              }
              Bool_t impact = kFALSE ;
              if(fPHOS){
                Int_t modulenum ;
                Double_t ztmp,xtmp ;
                impact=fPHOSgeom->ImpactOnEmc(partner,modulenum,ztmp,xtmp) ;
                if(fDebug>2){
                  printf("Impact on PHOS: module: %d, x tower: %f, z tower: %f\n", modulenum,xtmp,ztmp);
                }
              }
              else{
                impact = fEMCALgeom->Impact(partner) ;
              }
              if(!impact){ //this photon cannot hit PHOS
                if(fDebug>2)
                  printf("P_Geo\n");
                fhPartnerMissedGeo->Fill(partner->Pt());
                fhDecWMissedPartnerGeom0->Fill(p->Pt()) ;  //Spectrum of tagged with missed partner
                if(iFidArea>0){
                  fhDecWMissedPartnerGeom1->Fill(p->Pt()) ;  //Spectrum of tagged with missed partner
                  if(iFidArea>1){
                    fhDecWMissedPartnerGeom2->Fill(p->Pt()) ;  //Spectrum of tagged with missed partner
                    if(iFidArea>2){
                      fhDecWMissedPartnerGeom3->Fill(p->Pt()) ;  //Spectrum of tagged with missed partner
                    }
                  }
                }
                isPartnerLost=kTRUE;
              }
              if(!isPartnerLost && partner->Energy()<fMinEnergyCut){ //energy is not enough to be registered by PHOS
                if(fDebug>2)
                  printf("P_Reg, E=%f\n",partner->Energy());
                fhPartnerMissedEmin->Fill(partner->Pt());  //Spectrum of missed partners
                fhDecWMissedPartnerEmin->Fill(p->Pt()) ;  //Spectrum of tagged with missed partner
                isPartnerLost=kTRUE;
              }
              if(!isPartnerLost){
//                p->SetMCTagged(1); //set this photon as primary tagged
                fhDecWMCPartner->Fill(p->Pt());
                fhPartnerMCReg->Fill(partner->Pt());
                if(fDebug>2){
                  printf("both photons are inside PHOS. Energy: %f, Pt of pair photon: %f, E of pair photon: %f, Px: %f Py: %f Pz: %f, num of daughters: %d \n", caloCluster->E(),partner->Pt(),partner->Energy(),partner->Px(),partner->Py(),partner->Pz(),partner->GetNDaughters());
                }
              }
              else{
                fhDecWMissedPartnerAll->Fill(p->Pt());
              }
            }//Partner - photon
            else{//partner not photon
              fhDecWMissedPartnerNotPhoton->Fill(p->Pt());                
            }
          }
        }
      }
    }
  } //PHOS/EMCAL clusters
    
  if(fDebug>1)   
    printf("number of clusters: %d\n",inList);

  //Invariant Mass analysis
  for(Int_t phosPhoton1 = 0 ; phosPhoton1 < inList ; phosPhoton1++) {
    AliAODPWG4Particle * p1 = static_cast<AliAODPWG4Particle*>(fCaloPhotonsArr->At(phosPhoton1));        
    for(Int_t phosPhoton2 = 0 ; phosPhoton2 < inList ; phosPhoton2++) {
      if(phosPhoton1==phosPhoton2)
        continue ;
      AliAODPWG4Particle * p2 = static_cast<AliAODPWG4Particle*>(fCaloPhotonsArr->At(phosPhoton2));

      Double_t invMass = p1->GetPairMass(p2);
      for(Int_t iPID=0; iPID<4; iPID++){
        if(p1->IsPIDOK(pidMap[iPID],22))
          fhInvMassReal[iPID]->Fill(invMass,p1->Pt());
        if(p2->IsPIDOK(pidMap[iPID],22))
          fhInvMassReal[iPID]->Fill(invMass,p2->Pt());
      }
      if(fDebug>2)
          printf("Pair i=%d,j=%d, M=%f\n",phosPhoton1,phosPhoton2,invMass);

      Bool_t makePi0=IsInPi0Band(invMass,p1->Pt());

      if(makePi0 && p1->IsTagged()){//Multiple tagging
        fhTaggedMult->Fill(p1->Pt());
      }  
      if(makePi0 && !p1->IsTagged()){//Each photon should enter histogram once, even if tagged several times
        fhTaggedAll->Fill(p1->Pt());
        Int_t iFidArea = p1->GetFiducialArea(); 
        if(iFidArea>0){
          fhTaggedArea1->Fill(p1->Pt() ) ; 
          if(iFidArea>1){
            fhTaggedArea2->Fill(p1->Pt() ) ; 
            if(iFidArea>2){
              fhTaggedArea3->Fill(p1->Pt() ) ; 
            }
          }
        }

        for(Int_t iPID=0; iPID<4; iPID++){
          if(p1->IsPIDOK(pidMap[iPID],22))
            fhTaggedPID[iPID]->Fill(p1->Pt());
        }

        p1->SetTagged(kTRUE) ;
      }  
      
      //First chesk if this is true pi0 pair
      if(IsSamePi0(p1,p2)){ //Correctly tagged - from the same pi0
//        p1->SetTrueTagged(1);
//        p2->SetTrueTagged(1);
        if(makePi0)//Correctly tagged photons
          fhTaggedMCTrue->Fill(p1->Pt());
        else{ //Decay pair missed tagging      
          fhMCMissedTagging->Fill(p1->Pt());
          fhMCMissedTaggingMass->Fill(invMass,p1->Pt()) ;
          //Clussify why missed tagging (todo)
          //Converted
          //Partner not a photon
          //Tagged not a photon
          //Just wrong inv.mass          
        }  
      }
      else{//Fake tagged - not from the same pi0
        if(makePi0)//Fake pair
          fhMCFakeTagged->Fill(p1->Pt());
      }
    }
  }

  //Fill Mixed InvMass distributions:
  TIter nextEv(fEventList) ;
  while(TClonesArray * event2 = static_cast<TClonesArray*>(nextEv())){
    Int_t nPhotons2 = event2->GetEntriesFast() ;
    for(Int_t i=0; i < inList ; i++){
      AliAODPWG4Particle * p1 = static_cast<AliAODPWG4Particle*>(fCaloPhotonsArr->At(i)) ;
      for(Int_t j=0; j < nPhotons2 ; j++){
        AliAODPWG4Particle * p2 = static_cast<AliAODPWG4Particle*>(event2->At(j)) ;
        Double_t invMass = p1->GetPairMass(p2);
        for(Int_t iPID=0; iPID<4; iPID++){
          if(p1->IsPIDOK(pidMap[iPID],22))
            fhInvMassMixed[iPID]->Fill(invMass,p1->Pt());
          if(p2->IsPIDOK(pidMap[iPID],22))
            fhInvMassMixed[iPID]->Fill(invMass,p2->Pt());
        }
      }
    }
  }

  //Remove old events
  fEventList->AddFirst(fCaloPhotonsArr);
  if(fEventList->GetSize() > 10){
    TClonesArray *tmp = static_cast <TClonesArray*> (fEventList->Last());
    fEventList->Remove(tmp);
    delete tmp;
  }

  delete caloClustersArr;

  PostData(1, fOutputList);
}


//______________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::Init()
{
  // Intialisation of parameters
  AliInfo("Doing initialisation") ; 
  SetPhotonId(0.9) ; 
  SetMinEnergyCut(0.4);
  SetPi0MeanParameters(0.136,0.,0.0,0.0);
//  SetPi0MeanParameters(0.1377,-0.002566,0.001216,-0.0001256);
  SetPi0SigmaParameters(0.004508,0.005497,0.00000006);
}

//______________________________________________________________________________
void AliAnalysisTaskTaggedPhotons::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  if (fDebug > 1) Printf("Terminate()");
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskTaggedPhotons::IsInPi0Band(Double_t m, Double_t pt)const
{
  //Parameterization of the pi0 peak region
  Double_t mpi0mean = fPi0MeanP0 + fPi0MeanP1 * pt + fPi0MeanP2 * pt*pt + fPi0MeanP3 * pt*pt*pt;
  Double_t mpi0sigma = TMath::Sqrt(fPi0SigmaP0 * fPi0SigmaP0 / pt + fPi0SigmaP1 * fPi0SigmaP1 + fPi0SigmaP2 * fPi0SigmaP2 / pt / pt);
 
  return (m>mpi0mean-2*mpi0sigma && m<mpi0mean+2*mpi0sigma) ;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskTaggedPhotons::IsSamePi0(const AliAODPWG4Particle *p1, const AliAODPWG4Particle *p2)const{
  //Looks through parents and finds if there was commont pi0 among ancestors

  if(!fStack)
    return kFALSE ; //can not say anything

  Int_t prim1 = p1->GetLabel();
  while(prim1!=-1){ 
    Int_t prim2 = p2->GetLabel();
    while(prim2!=-1){ 
      if(prim1==prim2){
        if(fStack->Particle(prim1)->GetPdgCode()==111)
          return kTRUE ;
        else
          return kFALSE ;
      }
      prim2=fStack->Particle(prim2)->GetFirstMother() ;
    }
    prim1=fStack->Particle(prim1)->GetFirstMother() ;
  }
  return kFALSE ;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskTaggedPhotons::GetFiducialArea(Float_t * pos)const{
  //calculates in which kind of fiducial area photon hit
  Double_t phi=TMath::ATan2(pos[1],pos[0]) ;
  Double_t z=pos[2] ;
  while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
  while(phi<0.)phi+=TMath::TwoPi() ;
  if(fDebug>2){
    printf("FiducialArea: phi=%f, z=%f \n",phi,z) ;
        printf("   fZmax=%f, fZmin=%f, fPhimax=%f, fPhimin=%f \n",fZmax,fZmin,fPhimax,fPhimin) ;
  }
  if(fPHOS){
    //From active PHOS area remove bands in 10 cm
    const Double_t kphi=TMath::ATan(10./460.) ; //angular band width
    Double_t dzMax=TMath::Ceil((fZmax-z)/10.) ;
    Double_t dzMin=TMath::Ceil((z-fZmin)/10.) ;
    Double_t dphiMax=TMath::Ceil((fPhimax-phi)/kphi);
    Double_t dphiMin=TMath::Ceil((phi-fPhimin)/kphi);
    if(fDebug>2){
     printf("In PHOS \n") ;
    printf("    dzMax=%f, dzMin=%f, dphiMax=%f, dphiMin=%f ret=%d\n",dzMax,dzMin,dphiMax,dphiMin,(Int_t)TMath::Min(TMath::Min(dzMax,dzMin),TMath::Min(dphiMax,dphiMin))) ;
    }
    return (Int_t)TMath::Min(TMath::Min(dzMax,dzMin),TMath::Min(dphiMax,dphiMin)); 
  }
  else{//EMCAL
    //From active EMCAL area remove bands in 20 cm
    const Double_t kphi=TMath::ATan(20./428.) ; //angular band width
    Double_t dzMax=TMath::Ceil((fZmax-z)/20.) ;
    Double_t dzMin=TMath::Ceil((z-fZmin)/20.) ;
    Double_t dphiMax=TMath::Ceil((fPhimax-phi)/kphi);
    Double_t dphiMin=TMath::Ceil((phi-fPhimin)/kphi);
    return (Int_t)TMath::Min(TMath::Min(dzMax,dzMin),TMath::Min(dphiMax,dphiMin)); 
  }
}
//______________________________________________________________________________
Bool_t  AliAnalysisTaskTaggedPhotons::TestDisp(Double_t l0, Double_t l1, Double_t /*e*/)const{
  //test if dispersion corresponds to those of photon
  if(fPHOS){
//    Double_t l0mean=1.38736+0.490405*TMath::Exp(-e*0.286170) ;
//    Double_t l1mean=1.09786-0.323469*TMath::Exp(-e*0.918719) ;
//    Double_t l0sigma=0.159905+0.829831/e-0.158067/e/e ;
//    Double_t l1sigma=0.133170+0.404387/e-0.0426302/e/e ;
//    Double_t c =-0.382233 ; 
//    return ((l0-l0mean)*(l0-l0mean)/l0sigma/l0sigma + (l1-l1mean)*(l1-l1mean)/l1sigma/l1sigma+c*(l0-l0mean)*(l1-l1mean)/l0sigma/l1sigma)<1. ;

if(l1>= 0   && l0>= 0   && l1 < 0.1 && l0 < 0.1) return kFALSE;
if(l1>= 0   && l0 > 0.5 && l1 < 0.1 && l0 < 1.5) return kTRUE;
if(l1>= 0   && l0 > 2.0 && l1 < 0.1 && l0 < 2.7) return kFALSE;
if(l1>= 0   && l0 > 2.7 && l1 < 0.1 && l0 < 4.0) return kFALSE;
if(l1 > 0.1 && l1 < 0.7 && l0 > 0.7 && l0 < 2.1) return kTRUE;
if(l1 > 0.1 && l1 < 0.3 && l0 > 3.0 && l0 < 5.0) return kFALSE;
if(l1 > 0.3 && l1 < 0.7 && l0 > 2.5 && l0 < 4.0) return kFALSE;
if(l1 > 0.7 && l1 < 1.3 && l0 > 1.0 && l0 < 1.6) return kTRUE;
if(l1 > 0.7 && l1 < 1.3 && l0 > 1.6 && l0 < 3.5) return kTRUE; 
if(l1 > 1.3 && l1 < 3.5 && l0 > 1.3 && l0 < 3.5) return kTRUE; 
   return kFALSE ;

  }
  else{ //EMCAL: not ready yet
   return kTRUE ;

  }

}
//______________________________________________________________________________^M
Bool_t  AliAnalysisTaskTaggedPhotons::TestCharged(Double_t dr,Double_t /*en*/)const{

  if(dr<15.) return kFALSE ;
  return kTRUE ;

}
//______________________________________________________________________________^M
void  AliAnalysisTaskTaggedPhotons::InitGeometry(){
  //Read rotation matrixes from ESD

  if(fDebug>1)printf("Init geometry \n") ;
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent()) ;
  AliAODEvent * aod = 0x0 ;
  if(!esd)
    aod=dynamic_cast<AliAODEvent*>(InputEvent()) ;
  if(!esd && !aod)
    AliFatal("Can not read geometry matrixes from ESD/AOD: NO ESD") ;
  if(fPHOS){//reading PHOS matrixes
    fPHOSgeom = new AliPHOSGeoUtils("IHEP","");
    Bool_t activeMod[5]={0} ;
    for(Int_t mod=0; mod<5; mod++){
      if(esd){
        const TGeoHMatrix* m=esd->GetPHOSMatrix(mod) ;
        fPHOSgeom->SetMisalMatrix(m, mod) ;
        if(m!=0)
         activeMod[mod]=kTRUE ;
      }
      else{
        const TGeoHMatrix* m=aod->GetHeader()->GetPHOSMatrix(mod) ;
        fPHOSgeom->SetMisalMatrix(m, mod) ;
        if(m!=0)
         activeMod[mod]=kTRUE ;
      }
    }

    if(fZmax>fZmin) //already  set manually
      return ;
    fZmax= 999. ;
    fZmin=-999. ;
    fPhimax=-999. ;
    fPhimin= 999. ;
    for(Int_t imod=1; imod<=5; imod++){
      if(!activeMod[imod-1])
 	continue ;
    //Find exact coordinates of PHOS corners
      Int_t relId[4]={imod,0,1,1} ;
      Float_t x,z ;
      fPHOSgeom->RelPosInModule(relId,x,z) ;
      TVector3 pos ;
      fPHOSgeom->Local2Global(imod,x,z,pos) ;
      Double_t phi=pos.Phi() ;
      while(phi<0.)phi+=TMath::TwoPi() ;
      while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
      fPhimin=TMath::Min(fPhimin,float(phi)) ;
      fZmin=TMath::Max(fZmin,float(pos.Z())) ;
      relId[2]=64 ;
      relId[3]=56 ;
      fPHOSgeom->RelPosInModule(relId,x,z) ;
      fPHOSgeom->Local2Global(imod,x,z,pos) ;
      phi=pos.Phi() ;
      while(phi<0.)phi+=TMath::TwoPi() ;
      while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
      fPhimax=TMath::Max(fPhimax,float(phi)) ;
      fZmax=TMath::Min(fZmax,float(pos.Z())) ;
    }
  }
  else{ //EMCAL
    fEMCALgeom = new AliEMCALGeoUtils("EMCAL_COMPLETE");
    for(Int_t mod=0; mod < (fEMCALgeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){ //<---Gustavo, could you check???
      if(esd){
        const TGeoHMatrix* m=esd->GetEMCALMatrix(mod) ;
        fEMCALgeom->SetMisalMatrix(m, mod) ;
      }
      else{
        const TGeoHMatrix* m=aod->GetHeader()->GetEMCALMatrix(mod) ;
        fEMCALgeom->SetMisalMatrix(m, mod) ;
      }
    }
    if(fZmax>fZmin) //already  set manually
      return ;

    fZmax= 999. ;
    fZmin=-999. ;
    fPhimax=-999. ;
    fPhimin= 999. ;
    for(Int_t imod=0; imod<12; imod++){
    //Find exact coordinates of SM corners
      Int_t absId = fEMCALgeom->GetAbsCellIdFromCellIndexes(imod, 0, 0);
      TVector3 pos ;
      //Get the position of this tower.
      fEMCALgeom->RelPosCellInSModule(absId,pos);
      Double_t phi=pos.Phi() ;
      while(phi<0.)phi+=TMath::TwoPi() ;
      while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
      fPhimin=TMath::Min(fPhimin,float(phi)) ;
      fZmin=TMath::Max(fZmin,float(pos.Z())) ;
      absId = fEMCALgeom->GetAbsCellIdFromCellIndexes(imod, 24, 48);
      fEMCALgeom->RelPosCellInSModule(absId,pos);
      phi=pos.Phi() ;
      while(phi<0.)phi+=TMath::TwoPi() ;
      while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
      fPhimax=TMath::Max(fPhimax,float(phi)) ;
      fZmax=TMath::Min(fZmax,float(pos.Z())) ;
    }
  }
}
