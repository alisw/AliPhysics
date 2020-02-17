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

#include <TMath.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include "AliConst.h"
#include "AliDecayerPythia.h"
#include "AliGenPythiaPlus.h"
#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliPythiaBase.h"
#include "AliPythiaRndm.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliRunLoader.h"
#include "AliMC.h"
#include "AliLog.h"
#include "PyquenCommon.h"

ClassImp(AliGenPythiaPlus)

AliGenPythiaPlus::AliGenPythiaPlus():
    AliGenMC(),
    fPythia(0),
    fProcess(kPyCharm),          
    fStrucFunc(kCTEQ5L), 
    fKineBias(0.),
    fTrials(0),
    fTrialsRun(0),
    fQ(0.),
    fX1(0.),
    fX2(0.),
    fEventTime(0.),
    fInteractionRate(0.),
    fTimeWindow(0.),
    fCurSubEvent(0),
    fEventsTime(0),
    fNev(0),
    fFlavorSelect(0),
    fXsection(0.),
    fPtHardMin(0.),
    fPtHardMax(1.e4),
    fYHardMin(-1.e10),
    fYHardMax(1.e10),
    fGinit(1),
    fGfinal(1),
    fHadronisation(1),
    fNpartons(0),
    fReadFromFile(0),
    fQuench(0),
    fPtKick(1.),
    fFullEvent(kTRUE),
    fDecayer(0),
    fDebugEventFirst(-1),
    fDebugEventLast(-1),
    fEtMinJet(0.),      
    fEtMaxJet(1.e4),      
    fEtaMinJet(-20.),     
    fEtaMaxJet(20.),     
    fPhiMinJet(0.),     
    fPhiMaxJet(2.* TMath::Pi()),     
    fJetReconstruction(kCell),
    fEtaMinGamma(-20.),      
    fEtaMaxGamma(20.),      
    fPhiMinGamma(0.),      
    fPhiMaxGamma(2. * TMath::Pi()),      
    fUseYCutHQ(kFALSE),
    fYMinHQ(-20.),
    fYMaxHQ(20.),
    fPycellEtaMax(2.),     
    fPycellNEta(274),       
    fPycellNPhi(432),       
    fPycellThreshold(0.),  
    fPycellEtSeed(4.),     
    fPycellMinEtJet(10.),  
    fPycellMaxRadius(1.), 
    fStackFillOpt(kFlavorSelection),   
    fFeedDownOpt(kTRUE),    
    fFragmentation(kTRUE),
    fSetNuclei(kFALSE),
    fNewMIS(kFALSE),   
    fHFoff(kFALSE),    
    fTriggerParticle(0),
    fTriggerEta(0.9),
    fTriggerMultiplicity(0),
    fTriggerMultiplicityEta(0),
    fTriggerMultiplicityEtaMin(0),
    fTriggerMultiplicityEtaMax(0),
    fTriggerMultiplicityPtMin(0),  
    fCountMode(kCountAll),      
    fHeader(0),  
    fRL(0),      
    fFileName(0),
    fFragPhotonInCalo(kFALSE),
    fHadronInCalo(kFALSE) ,
    fPi0InCalo(kFALSE) ,
    fEtaInCalo(kFALSE) ,
    fPhotonInCalo(kFALSE), // not in use
    fDecayPhotonInCalo(kFALSE),
    fForceNeutralMeson2PhotonDecay(kFALSE),
    fEleInCalo(kFALSE),
    fEleInEMCAL(kFALSE), // not in use
    fCheckBarrel(kFALSE),
    fCheckBarrelCalos(kFALSE),
    fCheckEMCAL(kFALSE),
    fCheckPHOS(kFALSE),
    fCheckPHOSeta(kFALSE),
    fPHOSRotateCandidate(-1),
    fTriggerParticleMinPt(0), 
    fPhotonMinPt(0), // not in use
    fElectronMinPt(0), // not in use
    fPHOSMinPhi(250.),
    fPHOSMaxPhi(320.),
    fPHOSEta(0.13),
    fEMCALMinPhi(80.),
    fEMCALMaxPhi(187.),
    fEMCALEta(0.7),
    fDCALMinPhi(260.),
    fDCALMaxPhi(320.),
    fDCALMinEta(0.22),
    fDCALMaxEta(0.7),
    fDCALMinPhiThird(320.),
    fDCALMaxPhiThird(327.),
    fDCALEtaThird(0.7),
    fItune(-1), 
    fInfo(1) 
{
// Default Constructor
  fEnergyCMS = 5500.;
 
  if (!AliPythiaRndm::GetPythiaRandom()) 
      AliPythiaRndm::SetPythiaRandom(GetRandom());
}

AliGenPythiaPlus::AliGenPythiaPlus(AliPythiaBase* pythia)
    :AliGenMC(-1),
     fPythia(pythia),
     fProcess(kPyCharm),          
     fStrucFunc(kCTEQ5L), 
     fKineBias(0.),
     fTrials(0),
     fTrialsRun(0),
     fQ(0.),
     fX1(0.),
     fX2(0.),
     fEventTime(0.),
     fInteractionRate(0.),
     fTimeWindow(0.),
     fCurSubEvent(0),
     fEventsTime(0),
     fNev(0),
     fFlavorSelect(0),
     fXsection(0.),
     fPtHardMin(0.),
     fPtHardMax(1.e4),
     fYHardMin(-1.e10),
     fYHardMax(1.e10),
     fGinit(kTRUE),
     fGfinal(kTRUE),
     fHadronisation(kTRUE),
     fNpartons(0),
     fReadFromFile(kFALSE),
     fQuench(kFALSE),
     fPtKick(1.),
     fFullEvent(kTRUE),
     fDecayer(0),
     fDebugEventFirst(-1),
     fDebugEventLast(-1),
     fEtMinJet(0.),      
     fEtMaxJet(1.e4),      
     fEtaMinJet(-20.),     
     fEtaMaxJet(20.),     
     fPhiMinJet(0.),     
     fPhiMaxJet(2.* TMath::Pi()),     
     fJetReconstruction(kCell),
     fEtaMinGamma(-20.),      
     fEtaMaxGamma(20.),      
     fPhiMinGamma(0.),      
     fPhiMaxGamma(2. * TMath::Pi()),      
     fUseYCutHQ(kFALSE),
     fYMinHQ(-20.),
     fYMaxHQ(20.),
     fPycellEtaMax(2.),     
     fPycellNEta(274),       
     fPycellNPhi(432),       
     fPycellThreshold(0.),  
     fPycellEtSeed(4.),     
     fPycellMinEtJet(10.),  
     fPycellMaxRadius(1.), 
     fStackFillOpt(kFlavorSelection),   
     fFeedDownOpt(kTRUE),    
     fFragmentation(kTRUE),
     fSetNuclei(kFALSE),
     fNewMIS(kFALSE),   
     fHFoff(kFALSE),    
     fTriggerParticle(0),
     fTriggerEta(0.9),     
     fTriggerMultiplicity(0),
     fTriggerMultiplicityEta(0),
     fTriggerMultiplicityEtaMin(0),
     fTriggerMultiplicityEtaMax(0),
     fTriggerMultiplicityPtMin(0),  
     fCountMode(kCountAll),      
     fHeader(0),  
     fRL(0),      
     fFileName(0),
     fFragPhotonInCalo(kFALSE),
     fHadronInCalo(kFALSE) ,
     fPi0InCalo(kFALSE) ,
     fEtaInCalo(kFALSE) ,
     fPhotonInCalo(kFALSE), // not in use
     fDecayPhotonInCalo(kFALSE),
     fForceNeutralMeson2PhotonDecay(kFALSE),
     fEleInCalo(kFALSE),
     fEleInEMCAL(kFALSE), // not in use
     fCheckBarrel(kFALSE),
     fCheckBarrelCalos(kFALSE),
     fCheckEMCAL(kFALSE),
     fCheckPHOS(kFALSE),
     fCheckPHOSeta(kFALSE),
     fPHOSRotateCandidate(-1),
     fTriggerParticleMinPt(0), 
     fPhotonMinPt(0), // not in use
     fElectronMinPt(0), // not in use
     fPHOSMinPhi(250.),
     fPHOSMaxPhi(320.),
     fPHOSEta(0.13),
     fEMCALMinPhi(80.),
     fEMCALMaxPhi(187.),
     fEMCALEta(0.7),
     fDCALMinPhi(260.),
     fDCALMaxPhi(320.),
     fDCALMinEta(0.22),
     fDCALMaxEta(0.7),
     fDCALMinPhiThird(320.),
     fDCALMaxPhiThird(327.),
     fDCALEtaThird(0.7),
     fItune(-1),
     fInfo(1) 
{
// default charm production at 5. 5 TeV
// semimuonic decay
// structure function GRVHO
//
    fEnergyCMS = 5500.;
    fName = "Pythia";
    fTitle= "Particle Generator using PYTHIA";
    fDecayer = pythia->Decayer();
    SetForceDecay();
    // Set random number generator 
    if (!AliPythiaRndm::GetPythiaRandom()) 
      AliPythiaRndm::SetPythiaRandom(GetRandom());
 }

AliGenPythiaPlus::~AliGenPythiaPlus()
{
// Destructor
  if(fEventsTime) delete fEventsTime;
}

void AliGenPythiaPlus::SetInteractionRate(Float_t rate,Float_t timewindow)
{
// Generate pileup using user specified rate
    fInteractionRate = rate;
    fTimeWindow = timewindow;
    GeneratePileup();
}

void AliGenPythiaPlus::GeneratePileup()
{
// Generate sub events time for pileup
    fEventsTime = 0;
    if(fInteractionRate == 0.) {
      Warning("GeneratePileup","Zero interaction specified. Skipping pileup generation.\n");
      return;
    }

    Int_t npart = NumberParticles();
    if(npart < 0) {
      Warning("GeneratePileup","Negative number of particles. Skipping pileup generation.\n");
      return;
    }

    if(fEventsTime) delete fEventsTime;
    fEventsTime = new TArrayF(npart);
    TArrayF &array = *fEventsTime;
    for(Int_t ipart = 0; ipart < npart; ipart++)
      array[ipart] = 0.;

    Float_t eventtime = 0.;
    while(1)
      {
	eventtime += (AliPythiaRndm::GetPythiaRandom())->Exp(1./fInteractionRate);
	if(eventtime > fTimeWindow) break;
	array.Set(array.GetSize()+1);
	array[array.GetSize()-1] = eventtime;
      }

    eventtime = 0.;
    while(1)
      {
	eventtime -= (AliPythiaRndm::GetPythiaRandom())->Exp(1./fInteractionRate);
	if(TMath::Abs(eventtime) > fTimeWindow) break;
	array.Set(array.GetSize()+1);
	array[array.GetSize()-1] = eventtime;
      }

    SetNumberParticles(fEventsTime->GetSize());
}

void AliGenPythiaPlus::SetPycellParameters(Float_t etamax, Int_t neta, Int_t nphi,
				       Float_t thresh, Float_t etseed, Float_t minet, Float_t r)
{
// Set pycell parameters
    fPycellEtaMax    =  etamax;
    fPycellNEta      =  neta;
    fPycellNPhi      =  nphi;
    fPycellThreshold =  thresh;
    fPycellEtSeed    =  etseed;
    fPycellMinEtJet  =  minet;
    fPycellMaxRadius =  r;
}



void AliGenPythiaPlus::SetEventListRange(Int_t eventFirst, Int_t eventLast)
{
  // Set a range of event numbers, for which a table
  // of generated particle will be printed
  fDebugEventFirst = eventFirst;
  fDebugEventLast  = eventLast;
  if (fDebugEventLast==-1) fDebugEventLast=fDebugEventFirst;
}

void AliGenPythiaPlus::Init()
{
// Initialisation
    
//    SetMC(AliPythia::Instance());
//    fPythia=(AliPythia*) fMCEvGen;
  
    // Coeffs to go from mm / mm to meter / second
    SetGeneratorUnitsForMeterSecond(1.e-3, 1e-3/TMath::C()); 
  
//
    fParentWeight=1./Float_t(fNpart);
//

    
    fPythia->SetPtHardRange(fPtHardMin, fPtHardMax);
    fPythia->SetYHardRange(fYHardMin, fYHardMax);
    
    if (fAProjectile > 0 && fATarget > 0) fPythia->SetNuclei(fAProjectile, fATarget);  
    // Fragmentation?
    if (fFragmentation) {
      fPythia->SetFragmentation(1);
    } else {
      fPythia->SetFragmentation(0);
    }


//  initial state radiation   
    fPythia->SetInitialAndFinalStateRadiation(fGinit, fGfinal);

//  pt - kick
//  10/4/2017
//  this line has been removed since it overwrites the settings of the standard tunes
//  fPythia->SetIntrinsicKt(fPtKick);

    if (fReadFromFile) {
	fRL  =  AliRunLoader::Open(fFileName, "Partons");
	fRL->LoadKinematics();
	fRL->LoadHeader();
    } else {
	fRL = 0x0;
    }
 //
    fPythia->ProcInit(fProcess, fEnergyCMS, fStrucFunc, fItune);
    //  Forward Paramters to the AliPythia object
    fDecayer->SetForceDecay(fForceDecay);    
// Switch off Heavy Flavors on request  
    if (fHFoff) {
	fPythia->SwitchHFOff();
	// Switch off g->QQbar splitting in decay table
	((AliDecayerPythia*) fDecayer)->HeavyFlavourOff();
    }

    fDecayer->Init();


//  Parent and Children Selection
    switch (fProcess) 
    {
    case kPyOldUEQ2ordered:
    case kPyOldUEQ2ordered2:
    case kPyOldPopcorn:
      break;
    case kPyCharm:
    case kPyCharmUnforced:
    case kPyCharmPbPbMNR:
    case kPyCharmpPbMNR:
    case kPyCharmppMNR:
    case kPyCharmppMNRwmi:
    case kPyCharmPWHG:
	fParentSelect[0] =   411;
	fParentSelect[1] =   421;
	fParentSelect[2] =   431;
	fParentSelect[3] =  4122;
	fParentSelect[4] =  4232;
	fParentSelect[5] =  4132;
	fParentSelect[6] =  4332;
	fFlavorSelect    =  4;	
	break;
    case kPyD0PbPbMNR:
    case kPyD0pPbMNR:
    case kPyD0ppMNR:
	fParentSelect[0] =   421;
	fFlavorSelect    =   4;	
	break;
    case kPyDPlusPbPbMNR:
    case kPyDPluspPbMNR:
    case kPyDPlusppMNR:
	fParentSelect[0] =   411;
	fFlavorSelect    =   4;	
	break;
    case kPyDPlusStrangePbPbMNR:
    case kPyDPlusStrangepPbMNR:
    case kPyDPlusStrangeppMNR:
	fParentSelect[0] =   431;
	fFlavorSelect    =   4;	
	break;
    case kPyLambdacppMNR:
	fParentSelect[0] =  4122;
	fFlavorSelect    =   4;	
	break;      
    case kPyBeauty:
    case kPyBeautyJets:
    case kPyBeautyPbPbMNR:
    case kPyBeautypPbMNR:
    case kPyBeautyppMNR:
    case kPyBeautyppMNRwmi:
    case kPyBeautyPWHG:
	fParentSelect[0]=  511;
	fParentSelect[1]=  521;
	fParentSelect[2]=  531;
	fParentSelect[3]= 5122;
	fParentSelect[4]= 5132;
	fParentSelect[5]= 5232;
	fParentSelect[6]= 5332;
	fFlavorSelect   = 5;	
	break;
    case kPyHeavyFlavppMNRwmi:
	fParentSelect[0]=  511;  //settings to selct decay products
	fParentSelect[1]=  521;
	fParentSelect[2]=  531;
	fParentSelect[3]= 5122;
	fParentSelect[4]= 5132;
	fParentSelect[5]= 5232;
	fParentSelect[6]= 5332;
	fFlavorSelect    =  5;
	break;
    case kPyBeautyUnforced:
	fParentSelect[0] =  511;
	fParentSelect[1] =  521;
	fParentSelect[2] =  531;
	fParentSelect[3] = 5122;
	fParentSelect[4] = 5132;
	fParentSelect[5] = 5232;
	fParentSelect[6] = 5332;
	fFlavorSelect    = 5;	
	break;
    case kPyJpsiChi:
    case kPyJpsi:
	fParentSelect[0] = 443;
	break;
    case kPyMbAtlasTuneMC09:
    case kPyMbDefault:
    case kPyMb:
    case kPyMbWithDirectPhoton:
    case kPyMbNonDiffr:
    case kPyMbMSEL1:
    case kPyJets:
    case kPyJetsPWHG:
    case kPyDirectGamma:
    case kPyLhwgMb:	
	break;
    case kPyWPWHG:
    case kPyW:
    case kPyZ:
    case kPyZgamma:
    case kPyMBRSingleDiffraction:
    case kPyMBRDoubleDiffraction:
    case kPyMBRCentralDiffraction:
        break;
    }
//
//
//  JetFinder for Trigger
//
//  Configure detector (EMCAL like)
//
    fPythia->SetPycellParameters(fPycellEtaMax,fPycellNEta, fPycellNPhi,
				 fPycellThreshold, fPycellEtSeed, 
				 fPycellMinEtJet, fPycellMaxRadius);
//
//  This counts the total number of calls to Pyevnt() per run.
    fTrialsRun = 0;
    fQ         = 0.;
    fX1        = 0.;
    fX2        = 0.;    
    fNev       = 0 ;
//    
//
//
    AliGenMC::Init();
//
//
//  
    if (fSetNuclei) {
	fDyBoost = 0;
	Warning("Init","SetNuclei used. Use SetProjectile + SetTarget instead. fDyBoost has been reset to 0\n");
    }
    
    if (fQuench) {
	fPythia->InitQuenching(0., 0.1, 0.6e6, 0, 0.97, 30);
    }

//    fPythia->SetPARJ(200, 0.0);

//    if (fQuench == 3) {
//	// Nestor's change of the splittings
//	fPythia->SetPARJ(200, 0.8);
//	fPythia->SetMSTJ(41, 1);  // QCD radiation only
//	fPythia->SetMSTJ(42, 2);  // angular ordering
//	fPythia->SetMSTJ(44, 2);  // option to run alpha_s
//	fPythia->SetMSTJ(47, 0);  // No correction back to hard scattering element
//	fPythia->SetMSTJ(50, 0);  // No coherence in first branching
//	fPythia->SetPARJ(82, 1.); // Cut off for parton showers
//    }
}

void AliGenPythiaPlus::SetSeed(UInt_t seed)
{
  fPythia->SetSeed(seed);
}


void AliGenPythiaPlus::Generate()
{
// Generate one event
    
    fDecayer->ForceDecay();

    Double_t polar[3]   =   {0,0,0};
    Double_t origin[3]  =   {0,0,0};
    Double_t p[4];
//  converts from mm/c to s
    const Float_t kconv=0.001/TMath::C();
//
    Int_t nt=0;
    Int_t jev=0;
    Int_t j, kf;
    fTrials=0;
    fEventTime = 0.;
    
    

    //  Set collision vertex position 
    if (fVertexSmear == kPerEvent) Vertex();
    
//  event loop    
    while(1)
    {
//
// Produce event
//
//
// Switch hadronisation off
//
//	fPythia->SwitchHadronisationOff();
//
// Either produce new event or read partons from file
//	
	if (!fReadFromFile) {
	    if (!fNewMIS) {
		fPythia->GenerateEvent();
	    } else {
		fPythia->GenerateMIEvent();
	    }
	    fNpartons = fPythia->GetNumberOfParticles();
	} else {
	    printf("Loading Event %d\n",AliRunLoader::Instance()->GetEventNumber());
	    fRL->GetEvent(AliRunLoader::Instance()->GetEventNumber());
	    fPythia->SetNumberOfParticles(0);
	    fPythia->LoadEvent(fRL->Stack(), 0 , 1);
	    fPythia->EditEventList(21);
	}
	
//
//  Run quenching routine 
//
	if (fQuench == 1) {
	    fPythia->Quench();
	} else if (fQuench == 2){
	    fPythia->Pyquen(208., 0, 0.);
	} else if (fQuench == 3) {
	    // Quenching is via multiplicative correction of the splittings
	}
	
//
// Switch hadronisation on
//
//	fPythia->SwitchHadronisationOn();
//
// .. and perform hadronisation
//	printf("Calling hadronisation %d\n", fPythia->GetN());
//	fPythia->HadronizeEvent();	
	fTrials++;
	fPythia->GetParticles(&fParticles);
	Boost();
	if(TMath::Abs(fXingAngleY) > 1.e-10) BeamCrossAngle();
//
//
//
	Int_t i;
	
	fNprimaries = 0;
	Int_t np = fParticles.GetEntriesFast();
	
	if (np == 0) continue;
//
	
//
	Int_t* pParent   = new Int_t[np];
	Int_t* pSelected = new Int_t[np];
	Int_t* trackIt   = new Int_t[np];
	for (i = 0; i < np; i++) {
	  pParent[i]   = -1;
	  pSelected[i] =  0;
	  trackIt[i]   =  0;

	  if(fPythia->Version() == 8){
	    // order parent quarks by flavour
	    TParticle* iparticle = (TParticle *) fParticles.At(i);
	    Int_t pdgPart = TMath::Abs(iparticle->GetPdgCode());
	    if(pdgPart>=100){ // hadrons
	      Int_t kfl=pdgPart;
	      if (kfl > 100000) kfl %= 100000; // resonance
	      if (kfl > 10000)  kfl %= 10000;  // resonance
	      if (kfl > 10) kfl/=100; // meson
	      if (kfl > 10) kfl/=10; // baryon
	      Int_t iMo1 = iparticle->GetFirstMother();
	      Int_t iMo2 = iparticle->GetSecondMother();
	      if(iMo1 >=0 && iMo2 >=0){ // particle with two mothers
		TParticle* mother1 = (TParticle *) fParticles.At(iMo1);
		TParticle* mother2 = (TParticle *) fParticles.At(iMo2);
		Int_t absPdgMo1 = TMath::Abs(mother1->GetPdgCode());
		Int_t absPdgMo2 = TMath::Abs(mother2->GetPdgCode());
		if( (absPdgMo1<=6 || absPdgMo1==21) && (absPdgMo2<=6 || absPdgMo2==21) ){
		  // parton parents
		  if(absPdgMo1!=kfl && absPdgMo2==kfl){
		    // assign as first mother the quark with same flavour as the hadron
		    iparticle->SetFirstMother(iMo2);
		    iparticle->SetLastMother(iMo1);
		  }
		  if(absPdgMo1!=kfl && absPdgMo2!=kfl && absPdgMo1<absPdgMo2){
		    // in case no quark has the flavour of the hadron
		    // assign as first mother the one with higher pdg code
		    iparticle->SetFirstMother(iMo2);
		    iparticle->SetLastMother(iMo1);
		  }
		}
	      }
	    }
	  }
	}

	Int_t nc = 0;        // Total n. of selected particles
	Int_t nParents = 0;  // Selected parents
	Int_t nTkbles = 0;   // Trackable particles
	if (fProcess != kPyMbDefault && 
	    fProcess != kPyMb && 
	    fProcess != kPyMbWithDirectPhoton && 
	    fProcess != kPyJets && 
	    fProcess != kPyDirectGamma &&
	    fProcess != kPyMbNonDiffr  &&
	    fProcess != kPyMbMSEL1     &&
	    fProcess != kPyW && 
	    fProcess != kPyZ &&
      fProcess != kPyZgamma &&
	    fProcess != kPyCharmppMNRwmi && 
	    fProcess != kPyBeautyppMNRwmi &&
	    fProcess != kPyHeavyFlavppMNRwmi &&
      fProcess != kPyWPWHG &&
	    fProcess != kPyJetsPWHG &&
            fProcess != kPyCharmPWHG &&
     fProcess != kPyBeautyPWHG) {
	    
	    for (i = 0; i < np; i++) {
		TParticle* iparticle = (TParticle *) fParticles.At(i);
		Int_t ks = iparticle->GetStatusCode();
		kf = CheckPDGCode(iparticle->GetPdgCode());
// No initial state partons
		if (ks==21) continue;
//
// Heavy Flavor Selection
//
		// quark ?
		kf = TMath::Abs(kf);
		Int_t kfl = kf;
		// Resonance

		if (kfl > 100000) kfl %= 100000;
		if (kfl > 10000)  kfl %= 10000;
		// meson ?
		if  (kfl > 10) kfl/=100;
		// baryon
		if (kfl > 10) kfl/=10;
		Int_t ipa = (fPythia->Version() == 6) ? (iparticle->GetFirstMother() - 1) :(iparticle->GetFirstMother()) ;
		Int_t kfMo = 0;
//
// Establish mother daughter relation between heavy quarks and mesons
//
		if (kf >= fFlavorSelect && kf <= 6) {
		    Int_t idau = (fPythia->Version() == 6) ? (iparticle->GetFirstDaughter() - 1) :(iparticle->GetFirstDaughter());
		    if (idau > -1) {
			TParticle* daughter = (TParticle *) fParticles.At(idau);
			Int_t pdgD = daughter->GetPdgCode();
			if (pdgD == 91 || pdgD == 92) {
			    Int_t jmin = (fPythia->Version() == 6) ? (daughter->GetFirstDaughter() - 1) : (daughter->GetFirstDaughter());
			    Int_t jmax = (fPythia->Version() == 6) ? (daughter->GetLastDaughter() - 1)  : (daughter->GetLastDaughter());

			    for (Int_t jp = jmin; jp <= jmax; jp++)
				((TParticle *) fParticles.At(jp))->SetFirstMother(i+1);
			} // is string or cluster
		    } // has daughter
		} // heavy quark
		

		if (ipa > -1) {
		    TParticle *  mother = (TParticle *) fParticles.At(ipa);
		    kfMo = TMath::Abs(mother->GetPdgCode());
		}
		
		// What to keep in Stack?
		Bool_t flavorOK = kFALSE;
		Bool_t selectOK = kFALSE;
		if (fFeedDownOpt) {
		    if (kfl >= fFlavorSelect) flavorOK = kTRUE;
		} else {
		    if (kfl > fFlavorSelect) {
			nc = -1;
			break;
		    }
		    if (kfl == fFlavorSelect) flavorOK = kTRUE;
		}
		switch (fStackFillOpt) {
		case kHeavyFlavor:
		case kFlavorSelection:
		    selectOK = kTRUE;
		    break;
		case kParentSelection:
		    if (ParentSelected(kf) || kf <= 10) selectOK = kTRUE;
		    break;
		}
		if (flavorOK && selectOK) { 
//
// Heavy flavor hadron or quark
//
// Kinematic seletion on final state heavy flavor mesons
		    if (ParentSelected(kf) && !KinematicSelection(iparticle, 0)) 
		    {
			continue;
		    }
		    pSelected[i] = 1;
		    if (ParentSelected(kf)) ++nParents; // Update parent count
//		    printf("\n particle (HF)  %d %d %d", i, pSelected[i], kf);
		} else {
// Kinematic seletion on decay products
		    if (fCutOnChild && ParentSelected(kfMo) && ChildSelected(kf) 
			&& !KinematicSelection(iparticle, 1)) 
		    {
			continue;
		    }
//
// Decay products 
// Select if mother was selected and is not tracked

		    if (pSelected[ipa] && 
			!trackIt[ipa]  &&     // mother will be  tracked ?
			kfMo !=  5 &&         // mother is b-quark, don't store fragments          
			kfMo !=  4 &&         // mother is c-quark, don't store fragments 
			kf   != 92)           // don't store string
		    {
//
// Semi-stable or de-selected: diselect decay products:
// 
//
			if (pSelected[i] == -1 ||  fDecayer->GetLifetime(kf) > fMaxLifeTime)
			{
			    Int_t ipF = iparticle->GetFirstDaughter();
			    Int_t ipL = iparticle->GetLastDaughter();	
			    if (ipF > 0) for (j = ipF-1; j < ipL; j++) pSelected[j] = -1;
			}
//			printf("\n particle (decay)  %d %d %d", i, pSelected[i], kf);
			pSelected[i] = (pSelected[i] == -1) ? 0 : 1;
		    }
		}
		if (pSelected[i] == -1) pSelected[i] = 0;
		if (!pSelected[i]) continue;
		// Count quarks only if you did not include fragmentation
		if (fFragmentation && kf <= 10) continue;

		nc++;
// Decision on tracking
		trackIt[i] = 0;
//
// Track final state particle
		if (ks == 1) trackIt[i] = 1;
// Track semi-stable particles
		if ((ks == 1) || (fDecayer->GetLifetime(kf) > fMaxLifeTime))  trackIt[i] = 1;
// Track particles selected by process if undecayed. 
		if (fForceDecay == kNoDecay) {
		    if (ParentSelected(kf)) trackIt[i] = 1;
		} else {
		    if (ParentSelected(kf)) trackIt[i] = 0;
		}
		if (trackIt[i] == 1) ++nTkbles; // Update trackable counter
//
//

  	    } // particle selection loop
	    if (nc > 0) {
		for (i = 0; i < np; i++) {
		    if (!pSelected[i]) continue;
		    TParticle *  iparticle = (TParticle *) fParticles.At(i);
		    kf = CheckPDGCode(iparticle->GetPdgCode());
		    Int_t ks = iparticle->GetStatusCode();  
		    p[0] = iparticle->Px();
		    p[1] = iparticle->Py();
		    p[2] = iparticle->Pz();
		    p[3] = iparticle->Energy();
		    
		    origin[0] = fVertex[0]+iparticle->Vx()/10; // [cm]
		    origin[1] = fVertex[1]+iparticle->Vy()/10; // [cm]
		    origin[2] = fVertex[2]+iparticle->Vz()/10; // [cm]
		    
		    Float_t tof   = fTime + kconv*iparticle->T();
		    Int_t ipa = (fPythia->Version() == 6) ? (iparticle->GetFirstMother() - 1) :(iparticle->GetFirstMother()) ;
		    Int_t iparent = (ipa > -1) ? pParent[ipa] : -1;
 
		    PushTrack(fTrackIt*trackIt[i], iparent, kf, 
			      p[0], p[1], p[2], p[3], 
			      origin[0], origin[1], origin[2], tof, 
			      polar[0], polar[1], polar[2],
			      kPPrimary, nt, 1., ks);
		    pParent[i] = nt;
		    KeepTrack(nt);
		    fNprimaries++;
		} //  PushTrack loop
	    }
  	} else {
	    nc = GenerateMB();
	} // mb ?
	
	GetSubEventTime();

	delete[] pParent;
	delete[] pSelected;
	delete[] trackIt;

	if (nc > 0) {
	  switch (fCountMode) {
	  case kCountAll:
	    // printf(" Count all \n");
	    jev += nc;
	    break;
	  case kCountParents:
	    // printf(" Count parents \n");
	    jev += nParents;
	    break;
	  case kCountTrackables:
	    // printf(" Count trackable \n");
	    jev += nTkbles;
	    break;
	  }
	    if (jev >= fNpart || fNpart == -1) {
		fKineBias=Float_t(fNpart)/Float_t(fTrials);
		if (fInfo) fPythia->GetXandQ(fX1, fX2, fQ);
		fTrialsRun += fTrials;
		fNev++;
		MakeHeader();
		break;
	    }
	}
    } // event loop
    SetHighWaterMark(nt);
//  Adjust weight due to kinematic selection
//  AdjustWeights();
//  get cross-section
    fXsection = fPythia->GetXSection();
}

Int_t  AliGenPythiaPlus::GenerateMB()
{
//
// Min Bias selection and other global selections
//
    Int_t i, kf, nt, iparent;
    Int_t nc = 0;
    Double_t p[4];
    Double_t polar[3]   =   {0,0,0};
    Double_t origin[3]  =   {0,0,0};
//  converts from mm/c to s
    const Float_t kconv = 0.001 / 2.999792458e8;
    
    Int_t np = (fHadronisation) ? fParticles.GetEntriesFast() : fNpartons;
    Int_t* pParent = new Int_t[np];
    for (i=0; i< np; i++) pParent[i] = -1;
  
    if (fProcess == kPyJets || fProcess == kPyDirectGamma || fProcess == kPyJetsPWHG || 
        fProcess == kPyCharmPWHG || fProcess == kPyBeautyPWHG ) {
       // 6,7 particles in PYTHIA6

       TParticle* jet1 = (TParticle *) fParticles.At(4);
       TParticle* jet2 = (TParticle *) fParticles.At(5);

       if (!CheckTrigger(jet1, jet2)) {
	  delete [] pParent;
	  return 0;
	}
    }
  
    // Select events with fragmentation photon, decay photon, pi0, eta 
    // or other hadrons going to PHOS or EMCAL or central barrel,
    // implemented primarily for kPyJets, but extended to any kind of process.
    if ( ( fFragPhotonInCalo || fPi0InCalo || 
           fEtaInCalo        || fEleInCalo || 
           fHadronInCalo     || fDecayPhotonInCalo   ) && 
         ( fCheckPHOS || fCheckEMCAL || fCheckBarrel || fCheckBarrelCalos)    ) 
    {
      Bool_t ok = TriggerOnSelectedParticles(np);
      
      if ( !ok ) 
      {
        delete[] pParent;
        return 0;
      }
    }
    
    // Check for minimum multiplicity
    if (fTriggerMultiplicity > 0) {
      Int_t multiplicity = 0;
      for (i = 0; i < np; i++) {
	TParticle *  iparticle = (TParticle *) fParticles.At(i);
	
	Int_t statusCode = iparticle->GetStatusCode();
	
	// Initial state particle
	if (statusCode != 1)
	  continue;
	// eta cut
	if (fTriggerMultiplicityEta > 0 && TMath::Abs(iparticle->Eta()) > fTriggerMultiplicityEta)
	  continue;
	//multiplicity check for a given eta range
	if ((fTriggerMultiplicityEtaMin != fTriggerMultiplicityEtaMax) && 
	    (iparticle->Eta() < fTriggerMultiplicityEtaMin || iparticle->Eta() > fTriggerMultiplicityEtaMax))
	  continue;
	// pt cut
	if (iparticle->Pt() < fTriggerMultiplicityPtMin) 
	    continue;

	TParticlePDG* pdgPart = iparticle->GetPDG();
	if (pdgPart && pdgPart->Charge() == 0)
	  continue;
	
	++multiplicity;
      }

      if (multiplicity < fTriggerMultiplicity) {
	delete [] pParent;
	return 0;
      }
      Printf("Triggered on event with multiplicity of %d >= %d", multiplicity, fTriggerMultiplicity);
    }
    
    if (fTriggerParticle) {
	Bool_t triggered = kFALSE;
	for (i = 0; i < np; i++) {
	    TParticle *  iparticle = (TParticle *) fParticles.At(i);
	    kf = CheckPDGCode(iparticle->GetPdgCode());
	    if (kf != fTriggerParticle) continue;
	    if (iparticle->Pt() == 0.) continue;
	    if (TMath::Abs(iparticle->Eta()) > fTriggerEta) continue;
	    triggered = kTRUE;
	    break;
	}
	if (!triggered) {
	  delete [] pParent;
	  return 0;
	}
    }
	

    // Check if there is a ccbar or bbbar pair with at least one of the two
    // in fYMin < y < fYMax
    if (fProcess == kPyCharmppMNRwmi || fProcess == kPyBeautyppMNRwmi || fProcess == kPyHeavyFlavppMNRwmi ) {
      TParticle *partCheck;
      TParticle *mother;
      Bool_t  theQ=kFALSE,theQbar=kFALSE,inYcut=kFALSE;
      Bool_t  theChild=kFALSE;
      Float_t y;  
      Int_t   pdg,mpdg,mpdgUpperFamily;
      for(i=0; i<np; i++) {
	partCheck = (TParticle*)fParticles.At(i);
	pdg = partCheck->GetPdgCode();  
	Bool_t flavSel=kFALSE;
	if(TMath::Abs(pdg) == fFlavorSelect) flavSel=kTRUE;
	if(fProcess == kPyHeavyFlavppMNRwmi && (TMath::Abs(pdg) == 4 || TMath::Abs(pdg) == 5)) flavSel=kTRUE;
	if(flavSel) { // quark  
	  if(pdg>0) { theQ=kTRUE; } else { theQbar=kTRUE; }

      	if(partCheck->Energy()-TMath::Abs(partCheck->Pz()) > FLT_EPSILON) y = 0.5*TMath::Log((partCheck->Energy()+partCheck->Pz()+1.e-13)/
			     (partCheck->Energy()-partCheck->Pz()+1.e-13));
	  	else y = 9999.;
	  if(fUseYCutHQ && y>fYMinHQ && y<fYMaxHQ) inYcut=kTRUE;
	  if(!fUseYCutHQ && y>fYMin && y<fYMax) inYcut=kTRUE;
	}

	if(fCutOnChild && TMath::Abs(pdg) == fPdgCodeParticleforAcceptanceCut) {
	  Int_t mi = (fPythia->Version() == 6) ? (partCheck->GetFirstMother() - 1) :(partCheck->GetFirstMother()) ;
	  if(mi<0) continue;
	  mother = (TParticle*)fParticles.At(mi);
	  mpdg=TMath::Abs(mother->GetPdgCode());
	  mpdgUpperFamily=(mpdg>1000 ? mpdg+1000 : mpdg+100); // keep e from c from b
	  if ( ParentSelected(mpdg) || 
	      (fFlavorSelect==5 && ParentSelected(mpdgUpperFamily))) {
	    if (KinematicSelection(partCheck,1)) {
	      theChild=kTRUE;
	    }
	  }
	}
      }
      if (!theQ || !theQbar || !inYcut) { // one of the c/b conditions not satisfied
	delete[] pParent;
	return 0;
      }
      if (fCutOnChild && !theChild) { // one of the child conditions not satisfied
	delete[] pParent;
	return 0;	
      }

    }

    //Introducing child cuts in case kPyW, kPyZ, kPyMb, and kPyMbNonDiff
    if ( (
    fProcess == kPyWPWHG ||
    fProcess == kPyW ||
	  fProcess == kPyZ ||
    fProcess == kPyZgamma ||
	  fProcess == kPyMbDefault ||
	  fProcess == kPyMb ||
	  fProcess == kPyMbWithDirectPhoton ||
	  fProcess == kPyMbNonDiffr)  
	 && (fCutOnChild == 1) ) {
      if ( !CheckKinematicsOnChild() ) {
	delete[] pParent;
	return 0;
      }
    }
  

    for (i = 0; i < np; i++) {
	Int_t trackIt = 0;
	TParticle *  iparticle = (TParticle *) fParticles.At(i);
	kf = CheckPDGCode(iparticle->GetPdgCode());
	Int_t ks = iparticle->GetStatusCode();
	Int_t km = iparticle->GetFirstMother();
	if ((ks == 1  && kf!=0 && KinematicSelection(iparticle, 0)) ||
	    (ks != 1) ||
	    ((fProcess == kPyJets || fProcess == kPyJetsPWHG) && ks == 21 && km == 0 && i>1)) {

	    if(fStackFillOpt == kHeavyFlavor && !IsFromHeavyFlavor(i)) continue;
	    nc++;
	    if (ks == 1) trackIt = 1;

	    Int_t ipa = (fPythia->Version() == 6) ? (iparticle->GetFirstMother() - 1) :(iparticle->GetFirstMother()) ;
	    iparent = (ipa > -1) ? pParent[ipa] : -1;
	    if (ipa >= np) fPythia->EventListing();
	    
//
// store track information
	    p[0] = iparticle->Px();
	    p[1] = iparticle->Py();
	    p[2] = iparticle->Pz();
	    p[3] = iparticle->Energy();

	    
	    origin[0] = fVertex[0]+iparticle->Vx()/10; // [cm]
	    origin[1] = fVertex[1]+iparticle->Vy()/10; // [cm]
	    origin[2] = fVertex[2]+iparticle->Vz()/10; // [cm]
	    
	    Float_t tof = fTime + fEventTime + kconv * iparticle->T();

	    PushTrack(fTrackIt*trackIt, iparent, kf, 
		      p[0], p[1], p[2], p[3], 
		      origin[0], origin[1], origin[2], tof, 
		      polar[0], polar[1], polar[2],
		      kPPrimary, nt, 1., ks);
	    fNprimaries++;

	    
	    //
	    // Special Treatment to store color-flow
	    //
	    /*
	    if (ks == 3 || ks == 13 || ks == 14) {
		TParticle* particle = 0;
		if (fStack) {
		    particle = fStack->Particle(nt);
		} else {
  		    particle = AliRunLoader::Instance()->Stack()->Particle(nt);
		}
		particle->SetFirstDaughter(fPythia->GetK(2, i));
		particle->SetLastDaughter(fPythia->GetK(3, i));		
	    }
	    */  
	    KeepTrack(nt);
	    pParent[i] = nt;
	    SetHighWaterMark(nt);
	    
	} // select particle
    } // particle loop 

    delete[] pParent;
    
    return 1;
}


void AliGenPythiaPlus::FinishRun()
{
// Print x-section summary
    fPythia->PrintStatistics();

    if (fNev > 0.) {
	fQ  /= fNev;
	fX1 /= fNev;
	fX2 /= fNev;    
    }
    
    printf("\nTotal number of Pyevnt() calls %d\n", fTrialsRun);
    printf("\nMean Q, x1, x2: %f %f %f\n", fQ, fX1, fX2);
    WriteXsection();
}

void AliGenPythiaPlus::AdjustWeights() const
{
// Adjust the weights after generation of all events
//
    if (gAlice) {
	TParticle *part;
	Int_t ntrack=gAlice->GetMCApp()->GetNtrack();
	for (Int_t i=0; i<ntrack; i++) {
	    part= gAlice->GetMCApp()->Particle(i);
	    part->SetWeight(part->GetWeight()*fKineBias);
	}
    }
}
    
void AliGenPythiaPlus::SetNuclei(Int_t a1, Int_t a2)
{
// Treat protons as inside nuclei with mass numbers a1 and a2  
    fAProjectile = a1;
    fATarget     = a2;
    fSetNuclei   = kTRUE;
}


void AliGenPythiaPlus::MakeHeader()
{
//
// Make header for the simulated event
// 
  if (gAlice) {
    if (gAlice->GetEvNumber()>=fDebugEventFirst &&
	gAlice->GetEvNumber()<=fDebugEventLast) fPythia->EventListing();
  }

// Builds the event header, to be called after each event
    if (fHeader) delete fHeader;
    fHeader = new AliGenPythiaEventHeader("Pythia");
    fHeader->SetTitle(GetTitle());

//
// Event type  
    ((AliGenPythiaEventHeader*) fHeader)->SetProcessType(fPythia->ProcessCode());
//
// Number of trials
    ((AliGenPythiaEventHeader*) fHeader)->SetTrials(fTrials);
// Number of MPI
    ((AliGenPythiaEventHeader*) fHeader)->SetNMPI(fPythia->GetNMPI());
//
// Ncoll (for event superposition)
    ((AliGenPythiaEventHeader*) fHeader)->SetNSuperpositions(fPythia->GetNSuperpositions());
//
// Event Vertex 
    fHeader->SetPrimaryVertex(fVertex);
    fHeader->SetInteractionTime(fTime+fEventTime);    
//
// Number of primaries
    fHeader->SetNProduced(fNprimaries);
//
// Jets that have triggered

    if (fProcess == kPyJets || fProcess == kPyJetsPWHG)
    {
	Int_t ntrig, njet;
	Float_t jets[4][10];
	GetJets(njet, ntrig, jets);

	
	for (Int_t i = 0; i < ntrig; i++) {
	    ((AliGenPythiaEventHeader*) fHeader)->AddJet(jets[0][i], jets[1][i], jets[2][i], 
							jets[3][i]);
	}
    }
//
// Copy relevant information from external header, if present.
//
    Float_t uqJet[4];
    
    if (fRL) {
	AliGenPythiaEventHeader* exHeader = (AliGenPythiaEventHeader*) (fRL->GetHeader()->GenEventHeader());
	for (Int_t i = 0; i < exHeader->NTriggerJets(); i++)
	{
	    printf("Adding Jet %d %d \n", i,  exHeader->NTriggerJets());
	    
	    
	    exHeader->TriggerJet(i, uqJet);
	    ((AliGenPythiaEventHeader*) fHeader)->AddUQJet(uqJet[0], uqJet[1], uqJet[2], uqJet[3]);
	}
    }
//
// Store quenching parameters
//
    if (fQuench){
        Double_t z[4] = {0.};
	Double_t xp = 0.;
	Double_t yp = 0.;
	if (fQuench == 1) {
	    // Pythia::Quench()
	    fPythia->GetQuenchingParameters(xp, yp, z);
	} else {
	    // Pyquen
	    Double_t r1 = PARIMP.rb1;
	    Double_t r2 = PARIMP.rb2;
	    Double_t b  = PARIMP.b1;
	    Double_t r   = 0.5 * TMath::Sqrt(2. * (r1 * r1 + r2 * r2) - b * b);
	    Double_t phi = PARIMP.psib1;
	    xp = r * TMath::Cos(phi);
	    yp = r * TMath::Sin(phi);
	    
	}
	    ((AliGenPythiaEventHeader*) fHeader)->SetXYJet(xp, yp);
	    ((AliGenPythiaEventHeader*) fHeader)->SetZQuench(z);
	}
//
// Store pt^hard 
    ((AliGenPythiaEventHeader*) fHeader)->SetPtHard(fPythia->GetPtHard());
    ((AliGenPythiaEventHeader*) fHeader)->SetXsection(fPythia->GetXSection());

//
//  Pass header
//
    AddHeader(fHeader);
    fHeader = 0x0;
}

Bool_t AliGenPythiaPlus::CheckTrigger(const TParticle* jet1, const TParticle* jet2)
{
// Check the kinematic trigger condition
//
    Double_t eta[2];
    eta[0] = jet1->Eta();
    eta[1] = jet2->Eta();
    Double_t phi[2];
    phi[0] = jet1->Phi();
    phi[1] = jet2->Phi();
    Int_t    pdg[2]; 
    pdg[0] = jet1->GetPdgCode();
    pdg[1] = jet2->GetPdgCode();    
    Bool_t   triggered = kFALSE;

    if (fProcess == kPyJets || fProcess == kPyJetsPWHG) {
	Int_t njets = 0;
	Int_t ntrig = 0;
	Float_t jets[4][10];
//
// Use Pythia clustering on parton level to determine jet axis
//
	GetJets(njets, ntrig, jets);
	
	if (ntrig || fEtMinJet == 0.) triggered = kTRUE;
//
    } else {
	Int_t ij = 0;
	Int_t ig = 1;
	if (pdg[0] == kGamma) {
	    ij = 1;
	    ig = 0;
	}
      
      // Search the physical direct photon
      // and recover its eta and phi
      if ( fProcess == kPyDirectGamma )      
      {
        const TParticle * jets[] = {jet1,jet2};
        Int_t status = jets[ig]->GetStatusCode();
        Int_t index  = jets[ig]->GetDaughter(0);
        Int_t pdg    = jets[ig]->GetPdgCode();
        Int_t np     = fParticles.GetEntriesFast();
        
        //printf("Search physical photon...\n");
        TParticle * photon = 0;
        while ( status!=1 )
        {
          if ( pdg!=kGamma || index < 0 || index >= np)
          {
            AliWarning(Form("Photon with daughter PDG <%d> or negative/large index <%d>/<%d>, skip event",pdg,index,np));
            return kFALSE;
          }
          
          //printf("\t daught %d, status %d, pdg %d, eta %2.2f, phi %2.2f\n",index,status,pdg,eta[ig],phi[ig]);
          photon = (TParticle*) fParticles.At(index);
          if(!photon) return kFALSE;
          
          status = photon->GetStatusCode();
          index  = photon->GetDaughter(0);
          pdg    = photon->GetPdgCode();
          eta[ig]= photon->Eta();
          phi[ig]= photon->Phi();
        }
        //printf("final:   daught %d, status %d, pdg %d, eta %2.2f, phi %2.2f\n",index,status,pdg,eta[ig],phi[ig]);
        //printf("...found\n");
      }
      
	//Check eta range first...
	if ((eta[ij] < fEtaMaxJet   && eta[ij] > fEtaMinJet) &&
	    (eta[ig] < fEtaMaxGamma && eta[ig] > fEtaMinGamma))
	{
	    //Eta is okay, now check phi range
	    if ((phi[ij] < fPhiMaxJet   && phi[ij] > fPhiMinJet) &&
		(phi[ig] < fPhiMaxGamma && phi[ig] > fPhiMinGamma))
	    {
		triggered = kTRUE;
                if ( fCheckBarrelCalos )
                {
                  Float_t phiGJ = phi[ig]*180./TMath::Pi(); //Convert to degrees
                  if ( !IsInBarrelCalorimeters(phiGJ,TMath::Abs(eta[ig])) ) triggered = kFALSE;
                }
	    }
	}
    }
    return triggered;
}



Bool_t AliGenPythiaPlus::CheckKinematicsOnChild(){
//
//Checking Kinematics on Child (status code 1, particle code ?, kin cuts
//
    Bool_t checking = kFALSE;
    Int_t j, kcode, ks, km;
    Int_t nPartAcc = 0; //number of particles in the acceptance range
    Int_t numberOfAcceptedParticles = 1;
    if (fNumberOfAcceptedParticles != 0) { numberOfAcceptedParticles = fNumberOfAcceptedParticles; }
    Int_t npart = fParticles.GetEntriesFast();
    
    for (j = 0; j<npart; j++) {
	TParticle *  jparticle = (TParticle *) fParticles.At(j);
	kcode = TMath::Abs( CheckPDGCode(jparticle->GetPdgCode()) );
	ks = jparticle->GetStatusCode();
	km = jparticle->GetFirstMother(); 
	
	if( (ks == 1)  &&  (kcode == fPdgCodeParticleforAcceptanceCut)  &&  (KinematicSelection(jparticle,1)) ){
	    nPartAcc++;
	}
	if( numberOfAcceptedParticles <= nPartAcc){
	  checking = kTRUE;
	  break;
	}
    }

    return checking;
}

void AliGenPythiaPlus::RecJetsUA1(Int_t& njets, Float_t jets [4][50])
{
//
//  Calls the Pythia jet finding algorithm to find jets in the current event
//
//
//
//  Save jets
//
//  Run Jet Finder
    fPythia->Pycell(njets);
    Int_t i;
    for (i = 0; i < njets; i++) {
	Float_t px, py, pz, e;
	fPythia->GetJet(i, px, py, pz, e);
	jets[0][i] = px;
	jets[1][i] = py;
	jets[2][i] = pz;
	jets[3][i] = e;
    }
}


void  AliGenPythiaPlus::GetJets(Int_t& nJets, Int_t& nJetsTrig, Float_t jets[4][10])
{
//
//  Calls the Pythia clustering algorithm to find jets in the current event
//
    nJets       = 0;
    nJetsTrig   = 0;
    if (fJetReconstruction == kCluster) {
//
//  Configure cluster algorithm
//    
//	fPythia->SetPARU(43, 2.);
//	fPythia->SetMSTU(41, 1);
//
//  Call cluster algorithm
//    
	fPythia->Pyclus(nJets);
//
//  Loading jets from common block
//
    } else {

//
//  Run Jet Finder
	fPythia->Pycell(nJets);
    }

    Int_t i;
    for (i = 0; i < nJets; i++) {
	Float_t px, py, pz, e;
	fPythia->GetJet(i, px, py, pz, e);
	Float_t pt    = TMath::Sqrt(px * px + py * py);
	Float_t phi   = TMath::Pi() + TMath::ATan2(-py, -px);  
	Float_t theta = TMath::ATan2(pt,pz);
	Float_t et    = e * TMath::Sin(theta);
	Float_t eta   = -TMath::Log(TMath::Tan(theta / 2.));
	if (
	    eta > fEtaMinJet && eta < fEtaMaxJet && 
	    phi > fPhiMinJet && phi < fPhiMaxJet &&
	    et  > fEtMinJet  && et  < fEtMaxJet     
	    ) 
	{
	    jets[0][nJetsTrig] = px;
	    jets[1][nJetsTrig] = py;
	    jets[2][nJetsTrig] = pz;
	    jets[3][nJetsTrig] = e;
	    nJetsTrig++;
	} else {
	}
    }
}

void AliGenPythiaPlus::GetSubEventTime()
{
  // Calculates time of the next subevent
  fEventTime = 0.;
  if (fEventsTime) {
    TArrayF &array = *fEventsTime;
    fEventTime = array[fCurSubEvent++];
  }
  //  printf(" Event time: %d %f %p",fCurSubEvent,fEventTime,fEventsTime);
  return;
}

///
/// \return true if particle in central barrel acceptance
///
/// \param eta        particle pseudorapidity,  etamin=-etamax
///
Bool_t AliGenPythiaPlus::IsInBarrel(Float_t eta) const
{
  return ( eta < fTriggerEta ) ;
}

///
/// \return true if particle in EMCAL or DCAL or PHOS acceptance
/// Acceptance slightly larger to be considered.
///
/// \param phi        particle azimuthal angle in degrees
/// \param eta        particle pseudorapidity,  etamin=-etamax
///
Bool_t AliGenPythiaPlus::IsInBarrelCalorimeters(Float_t phi, Float_t eta) 
{
  if      ( IsInEMCAL(phi,eta   ) ) return kTRUE ;
  else if ( IsInDCAL (phi,eta   ) ) return kTRUE ;
  else if ( IsInPHOS (phi,eta,-1) ) return kTRUE ;
  else                              return kFALSE;
}

///
/// \return true if particle in EMCAL acceptance
/// Acceptance slightly larger to be considered.
///
/// \param phi        particle azimuthal angle in degrees
/// \param eta        particle pseudorapidity,  etamin=-etamax
///
Bool_t AliGenPythiaPlus::IsInEMCAL(Float_t phi, Float_t eta) const
{
  if(phi > fEMCALMinPhi  && phi < fEMCALMaxPhi && 
     eta < fEMCALEta  ) 
    return kTRUE;
  else 
    return kFALSE;
}

///
/// \return true if particle in DCAL acceptance
/// Acceptance slightly larger to be considered.
///
/// \param phi        particle azimuthal angle in degrees
/// \param eta        particle pseudorapidity,  etamin=-etamax
///
Bool_t AliGenPythiaPlus::IsInDCAL(Float_t phi, Float_t eta) const
{
  Bool_t fullSM  = kFALSE;
  Bool_t thirdSM = kFALSE;
  
  if(phi > fDCALMinPhi  && phi < fDCALMaxPhi && 
     eta > fDCALMinEta  && eta < fDCALMaxEta   ) fullSM = kTRUE;
  
  if(phi > fDCALMinPhiThird  && phi < fDCALMaxPhiThird && 
     eta < fDCALEtaThird  ) thirdSM = kTRUE;
  
  if ( fullSM || thirdSM )
    return kTRUE;
  else 
    return kFALSE;
}

///
/// \return true if particle in PHOS acceptance
/// Acceptance slightly larger to be considered.
///
/// \param phi        particle azimuthal angle in degrees
/// \param eta        particle pseudorapidity,  etamin=-etamax
/// \param iparticle  particle index, for PHOS event rotation
///
Bool_t AliGenPythiaPlus::IsInPHOS(Float_t phi, Float_t eta, Int_t iparticle)
{
  if(phi > fPHOSMinPhi  && phi < fPHOSMaxPhi && 
     eta < fPHOSEta  ) 
    return kTRUE;
  else 
  {
    if( fCheckPHOSeta && eta < fPHOSEta) fPHOSRotateCandidate = iparticle;
    
    return kFALSE;
  }
}

///
/// Rotate event in phi to enhance events in PHOS acceptance
///
void AliGenPythiaPlus::RotatePhi(Bool_t& okdd)
{  
  if(fPHOSRotateCandidate < 0) return ; 
  
  //calculate the new position random between fPHOSMinPhi and fPHOSMaxPhi 
  Double_t phiPHOSmin = TMath::Pi()*fPHOSMinPhi/180;
  Double_t phiPHOSmax = TMath::Pi()*fPHOSMaxPhi/180;
  Double_t phiPHOS = (AliPythiaRndm::GetPythiaRandom())->Uniform(phiPHOSmin,phiPHOSmax);
  
  //calculate deltaphi
  TParticle* ph = (TParticle *) fParticles.At(fPHOSRotateCandidate);
  Double_t phphi = ph->Phi();
  Double_t deltaphi = phiPHOS - phphi;
  
  //loop for all particles and produce the phi rotation
  Int_t np = (fHadronisation) ? fParticles.GetEntriesFast() : fNpartons;
  Double_t oldphi, newphi;
  Double_t newVx, newVy, r, vZ, time; 
  Double_t newPx, newPy, pt, pz, e;
  for(Int_t i=0; i< np; i++) 
  {
    TParticle* iparticle = (TParticle *) fParticles.At(i);
    oldphi = iparticle->Phi();
    newphi = oldphi + deltaphi;
    if(newphi < 0) newphi = 2*TMath::Pi() + newphi; // correct angle 
    if(newphi > 2*TMath::Pi()) newphi = newphi - 2*TMath::Pi(); // correct angle
    
    r = iparticle->R();
    newVx = r * TMath::Cos(newphi);
    newVy = r * TMath::Sin(newphi);
    vZ   = iparticle->Vz(); // don't transform
    time = iparticle->T(); // don't transform
    
    pt = iparticle->Pt();
    newPx = pt * TMath::Cos(newphi);
    newPy = pt * TMath::Sin(newphi);
    pz = iparticle->Pz(); // don't transform
    e  = iparticle->Energy(); // don't transform
    
    // apply rotation 
    iparticle->SetProductionVertex(newVx, newVy, vZ, time);
    iparticle->SetMomentum(newPx, newPy, pz, e);
    
  } //end particle loop 
  
  // now let's check that we put correctly the candidate photon in PHOS
  Float_t phi = ph->Phi()*180./TMath::Pi(); //Convert to degrees
  Float_t eta =TMath::Abs(ph->Eta());//in calos etamin=-etamax 
  if(IsInPHOS(phi,eta,-1)) 
  okdd = kTRUE;
  
  // reset the value for next event
  fPHOSRotateCandidate = -1;
}

/// Check if this is a heavy flavor decay product
Bool_t AliGenPythiaPlus::IsFromHeavyFlavor(Int_t ipart)
{
  TParticle *  part = (TParticle *) fParticles.At(ipart);
  Int_t mpdg = TMath::Abs(part->GetPdgCode());
  Int_t mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));
  if (mfl >= 4 && mfl < 6) return kTRUE; // HF hadron
  
  // seach if originates from HF haron decay
  Int_t imo = (fPythia->Version() == 6) ? (part->GetFirstMother()-1) : (part->GetFirstMother());
  TParticle* pm = part;
  while (imo >  -1) {
    pm  =  (TParticle*)fParticles.At(imo);
    mpdg = TMath::Abs(pm->GetPdgCode());
    mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));
    if ((mfl > 3) && (mfl <6) && mpdg > 400) return kTRUE;
    imo = (fPythia->Version() == 6) ? (pm->GetFirstMother()-1) : (pm->GetFirstMother());
  }
  return kFALSE;
}

///
/// Check the eta/phi correspond to the detectors acceptance
/// iparticle is the index of the particle to be checked, for PHOS rotation case
///
/// \return true if particle in selected acceptance
/// \param phi        particle azimuthal angle in degrees
/// \param eta        particle pseudorapidity,  etamin=-etamax
/// \param iparticle  particle index, for PHOS event rotation
Bool_t AliGenPythiaPlus::CheckDetectorAcceptance(Float_t phi, Float_t eta, Int_t iparticle)
{
  if     (fCheckPHOS   && IsInPHOS  (phi,eta,iparticle)) return kTRUE;
  else if(fCheckEMCAL  && IsInEMCAL (phi,eta)) return kTRUE;
  else if(fCheckBarrel && IsInBarrel(    eta)) return kTRUE;
  else if(fCheckBarrelCalos && IsInBarrelCalorimeters(phi,eta)) return kTRUE;
  else                                         return kFALSE;
}

///
/// Select events with fragmentation photon, decay photon, pi0, eta or other 
/// hadrons going to PHOS or EMCAL or central barrel,
/// implemented primarily for kPyJets, but extended to any kind of process.
///
/// \param np : number of generated particles
/// \return true if accept the event
///
Bool_t AliGenPythiaPlus::TriggerOnSelectedParticles(Int_t np)
{
  AliDebug(1,Form("** Check: frag photon %d, pi0 %d, eta %d, electron %d, hadron %d, decay %d; "
                  " in PHOS %d, EMCAL %d, Barrel %d; Barrel calos %d, with pT > %2.2f GeV/c**",
                  fFragPhotonInCalo,fPi0InCalo, fEtaInCalo,fEleInCalo,fHadronInCalo,fDecayPhotonInCalo,
                  fCheckPHOS,fCheckEMCAL, fCheckBarrel,fCheckBarrelCalos,fTriggerParticleMinPt)); 
    
  Bool_t ok = kFALSE;
  for (Int_t i=0; i< np; i++) 
  {
    TParticle* iparticle = (TParticle *) fParticles.At(i);
    
    Int_t pdg          = iparticle->GetPdgCode();
    Int_t status       = iparticle->GetStatusCode();
    Int_t imother      = iparticle->GetFirstMother();//- 1;
    
    TParticle* pmother = 0x0;
    Int_t momStatus    = -1;
    Int_t momPdg       = -1;
    if(imother > 0 )
    {  
      pmother = (TParticle *) fParticles.At(imother);
      momStatus    = pmother->GetStatusCode();
      momPdg       = pmother->GetPdgCode();
    }
    
    ok = kFALSE;

    //
    // Check the particle type: hadron (not pi0 or eta), electron, decay photon (from pi0 or eta or any), pi0 or eta
    //
    
    // Hadron
    if (fHadronInCalo && status == 1)
    {
      // avoid quarks, bosons, photons, electrons, muons, neutrinos and eta or pi0 
      // (in case neutral mesons were declared stable)
      if ( TMath::Abs(pdg) > 100 && pdg != 221 && pdg != 111 ) 
      ok = kTRUE;
    }
    
    // Electron
    else if ( fEleInCalo && status == 1 && TMath::Abs(pdg) == 11 )
    {
      ok = kTRUE;
    }
    
    // Fragmentation photon
    else if ( fFragPhotonInCalo && pdg == 22 && status == 1 && TMath::Abs(momPdg) < 23 ) 
    {        
      // mother can be only q/g/strings/photon, no photon from hadron decay
      ok = kTRUE ;  
    }
    
    // Decay photon
    else if ( fDecayPhotonInCalo && !fForceNeutralMeson2PhotonDecay && pdg == 22 && TMath::Abs(momPdg) > 100 ) 
    {              
      if( momStatus == 0)
      {
        //if(iparticle->Pt() > fTriggerParticleMinPt) 
        //printf("Decay photon! pdg %d, status %d, pt %2.2f, mom: pdg %d, pt %2.2f\n",
        //        pdg,status,iparticle->Pt(),momPdg,pmother->Pt());
        ok = kTRUE ;  // photon from hadron decay
        
        // In case only decays from pi0 or eta requested
        if ( fPi0InCalo && momPdg!=111 ) ok = kFALSE;
        if ( fEtaInCalo && momPdg!=221 ) ok = kFALSE;
      }
    }
    // Pi0 or Eta particle
    else if ( fPi0InCalo || fEtaInCalo )
    {
      if ( fDecayPhotonInCalo && !fForceNeutralMeson2PhotonDecay ) continue ;
      
      if (fPi0InCalo && pdg == 111) // pi0 status can be 1 or 11 depending on decay settings
      {
        //if(iparticle->Pt() > fTriggerParticleMinPt) printf("Pi0! pdg %d, status %d, pt %2.2f\n",pdg,status,iparticle->Pt());
        ok = kTRUE;
      }      
      else if (fEtaInCalo && pdg == 221) 
      {
        //if(iparticle->Pt() > fTriggerParticleMinPt) printf("Eta! pdg %d, status %d, pt %2.2f\n",pdg,status,iparticle->Pt());
        ok = kTRUE;
      }
      
    }// pi0 or eta
    
    //
    // Check that the selected particle is in the calorimeter acceptance
    //
    if(ok && iparticle->Pt() > fTriggerParticleMinPt)
    {
      // Just check if the selected particle falls in the acceptance
      if(!fForceNeutralMeson2PhotonDecay )
      {
        //printf("\t Check acceptance! \n");
        Float_t phi = iparticle->Phi()*180./TMath::Pi(); //Convert to degrees
        Float_t eta =TMath::Abs(iparticle->Eta()); //in calos etamin=-etamax	  
        
        if(CheckDetectorAcceptance(phi,eta,i))
        {
          ok =kTRUE;
          
          AliInfo(Form("Selected trigger pdg %d, status %d, pt %2.2f, eta %2.2f, phi %2.2f; mother pdg %d status %d",
                       pdg,status,iparticle->Pt(), eta, phi, momPdg, momStatus));
          
          break;
        }
        else ok = kFALSE;
      }
      // Mesons have several decay modes, select only those decaying into 2 photons
      else if ( fForceNeutralMeson2PhotonDecay && (fPi0InCalo || fEtaInCalo) )
      {
        // In case we want the pi0/eta trigger, 
        // check the decay mode (2 photons)
        
        //printf("\t Force decay 2 gamma\n");          
        
        Int_t ndaughters = iparticle->GetNDaughters();
        if(ndaughters != 2)
        {
          ok=kFALSE;
          continue;
        }
        
        TParticle*d1 = (TParticle *) fParticles.At(iparticle->GetDaughter(0)-1);
        TParticle*d2 = (TParticle *) fParticles.At(iparticle->GetDaughter(1)-1);
        if(!d1 || !d2) 
        {
          ok=kFALSE;
          continue;
        }
        
        //iparticle->Print();
        //d1->Print();
        //d2->Print();
        
        Int_t pdgD1 = d1->GetPdgCode();
        Int_t pdgD2 = d2->GetPdgCode();
        //printf("\t \t 2 daughters pdg = %d - %d\n",pdgD1,pdgD2);
        //printf("mother %d - %d\n",d1->GetFirstMother(),d2->GetFirstMother());
        
        if(pdgD1 != 22  || pdgD2 != 22)
        { 
          ok = kFALSE;
          continue;
        }
        
        //printf("\t accept decay\n");
        
        // Trigger on the meson, not on the daughter
        if ( !fDecayPhotonInCalo )
        {    
          Float_t phi = iparticle->Phi()*180./TMath::Pi(); //Convert to degrees
          Float_t eta =TMath::Abs(iparticle->Eta()); //in calos etamin=-etamax	
          
          if(CheckDetectorAcceptance(phi,eta,i))
          {
            //printf("\t Accept meson pdg %d\n",pdg);
            ok =kTRUE;
            
            AliDebug(1,Form("Selected trigger pdg %d (decay), status %d, pt %2.2f, eta %2.2f, phi %2.2f; mother pdg %d status %d",
                            pdg,status,iparticle->Pt(), eta, phi, momPdg, momStatus));
            
            break;
          } 
          else
          {
            ok=kFALSE;
            continue;
          }
        }
        
        //printf("Check daughters acceptance\n");
        
        // Trigger on the meson daughters
        // Photon 1
        Float_t phi = d1->Phi()*180./TMath::Pi(); //Convert to degrees
        Float_t eta =TMath::Abs(d1->Eta()); //in calos etamin=-etamax	  
        if(d1->Pt() > fTriggerParticleMinPt)
        {
          if(CheckDetectorAcceptance(phi,eta,i))
          {
            ok =kTRUE;
            
            AliDebug(1,Form("Selected trigger pdg %d (decay), status %d, pt %2.2f, eta %2.2f, phi %2.2f; mother pdg %d status %d",
                            pdg,status,iparticle->Pt(), eta, phi, momPdg, momStatus));
            
            break;
          }
          else ok = kFALSE;
        } // pt cut
        else  ok = kFALSE;
        
        //Photon 2
        phi = d2->Phi()*180./TMath::Pi(); //Convert to degrees
        eta =TMath::Abs(d2->Eta()); //in calos etamin=-etamax	
        
        if(d2->Pt() > fTriggerParticleMinPt)
        {
          if(CheckDetectorAcceptance(phi,eta,i))
          {
            ok =kTRUE;

            AliDebug(1,Form("Selected trigger pdg %d (decay), status %d, pt %2.2f, eta %2.2f, phi %2.2f\n",
                            pdg,status,iparticle->Pt(), eta, phi));
            
            break;
          } 
          else ok = kFALSE;         
        } // pt cut
        else ok = kFALSE;
      } // force 2 photon daughters in pi0/eta decays
      else ok = kFALSE;
    } else ok = kFALSE; // check acceptance
  } // primary loop
  
  //
  // If requested, rotate the particles event in phi to enhance/speed PHOS selection
  // A particle passing all trigger conditions except phi position in PHOS, is used as reference
  //
  if(fCheckPHOSeta)
  {
    RotatePhi(ok);
  }
  
  return ok;
}

void AliGenPythiaPlus::WriteXsection(const Char_t *fname) {
  //
  // Write cross section and Ntrials to a tree in a file
  //   Used for pt-hard bin productions
  //
  TFile fout(fname,"recreate");
  TTree tree("Xsection","Pythia cross section");
  // Convert to expected types for backwards compatibility
  Double_t xsec = fXsection;
  UInt_t trials = fTrialsRun;
  tree.Branch("xsection", &xsec, "X/D");
  tree.Branch("ntrials" , &trials , "X/i");
  tree.Fill();
  tree.Write();
  fout.Close();
}

#ifdef never
void AliGenPythiaPlus::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliGenPythia.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliGenerator::Streamer(R__b);
      R__b >> (Int_t&)fProcess;
      R__b >> (Int_t&)fStrucFunc;
      R__b >> (Int_t&)fForceDecay;
      R__b >> fEnergyCMS;
      R__b >> fKineBias;
      R__b >> fTrials;
      fParentSelect.Streamer(R__b);
      fChildSelect.Streamer(R__b);
      R__b >> fXsection;
//      (AliPythia::Instance())->Streamer(R__b);
      R__b >> fPtHardMin;
      R__b >> fPtHardMax;
//      if (fDecayer) fDecayer->Streamer(R__b);
   } else {
      R__b.WriteVersion(AliGenPythiaPlus::IsA());
      AliGenerator::Streamer(R__b);
      R__b << (Int_t)fProcess;
      R__b << (Int_t)fStrucFunc;
      R__b << (Int_t)fForceDecay;
      R__b << fEnergyCMS;
      R__b << fKineBias;
      R__b << fTrials;
      fParentSelect.Streamer(R__b);
      fChildSelect.Streamer(R__b);
      R__b << fXsection;
//      R__b << fPythia;
      R__b << fPtHardMin;
      R__b << fPtHardMax;
      //     fDecayer->Streamer(R__b);
   }
}
#endif



