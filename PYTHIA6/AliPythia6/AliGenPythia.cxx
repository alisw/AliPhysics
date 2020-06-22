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

//
// Generator using the TPythia interface (via AliPythia)
// to generate pp collisions.
// Using SetNuclei() also nuclear modifications to the structure functions
// can be taken into account. This makes, of course, only sense for the
// generation of the products of hard processes (heavy flavor, jets ...)
//
// andreas.morsch@cern.ch
//

#include <TMath.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>
#include "AliConst.h"
#include "AliDecayerPythia.h"
#include "AliGenPythia.h"
#include "AliFastGlauber.h"
#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliPythia.h"
#include "AliPythiaRndm.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliRunLoader.h"
#include "AliMC.h"
#include "AliLog.h"
#include "PyquenCommon.h"

ClassImp(AliGenPythia)


AliGenPythia::AliGenPythia():
    AliGenMC(),
    fProcess(kPyCharm),
    fItune(-1),
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
    fPythia(0),
    fWeightPower(0.),
    fPtHardMin(0.),
    fPtHardMax(1.e4),
    fYHardMin(-1.e10),
    fYHardMax(1.e10),
    fGinit(1),
    fGfinal(1),
    fCRoff(0),
    fHadronisation(1),
    fPatchOmegaDalitz(0),
    fDecayerExodus(0),
    fNpartons(0),
    fReadFromFile(0),
    fReadLHEF(0),
    fQuench(0),
    fQhat(0.),
    fLength(0.),
    fpyquenT(1.),
    fpyquenTau(0.1),
    fpyquenNf(0),
    fpyquenEloss(0),
    fpyquenAngle(0),
    fImpact(0.),
    fPtKick(1.),
    fFullEvent(kTRUE),
    fDecayer(new AliDecayerPythia()),
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
    fUseNuclearPDF(kFALSE),
    fUseLorentzBoost(kTRUE),
    fNewMIS(kFALSE),
    fHFoff(kFALSE),
    fNucPdf(0),
    fTriggerParticle(0),
    fTriggerEta(0.9),
    fTriggerY(999.),
    fTriggerEtaMin(0.9),
    fTriggerMinPt(-1),
    fTriggerMaxPt(1000),
    fTriggerMultiplicity(0),
    fTriggerMultiplicityEta(0),
    fTriggerMultiplicityEtaMin(0),
    fTriggerMultiplicityEtaMax(0),
    fTriggerMultiplicityPtMin(0),
    fCountMode(kCountAll),
    fHeader(0),
    fRL(0),
    fkFileName(0),
    fkNameLHEF(0),
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
    fkTuneForDiff(0),
    fProcDiff(0)
{
// Default Constructor
  fEnergyCMS = 5500.;
  if (!AliPythiaRndm::GetPythiaRandom())
      AliPythiaRndm::SetPythiaRandom(GetRandom());
}

AliGenPythia::AliGenPythia(Int_t npart)
    :AliGenMC(npart),
     fProcess(kPyCharm),
     fItune(-1),
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
     fPythia(0),
     fWeightPower(0.),
     fPtHardMin(0.),
     fPtHardMax(1.e4),
     fYHardMin(-1.e10),
     fYHardMax(1.e10),
     fGinit(kTRUE),
     fGfinal(kTRUE),
     fCRoff(kFALSE),
     fHadronisation(kTRUE),
     fPatchOmegaDalitz(0),
     fDecayerExodus(0),
     fNpartons(0),
     fReadFromFile(kFALSE),
     fReadLHEF(0),
     fQuench(kFALSE),
     fQhat(0.),
     fLength(0.),
     fpyquenT(1.),
     fpyquenTau(0.1),
     fpyquenNf(0),
     fpyquenEloss(0),
     fpyquenAngle(0),
     fImpact(0.),
     fPtKick(1.),
     fFullEvent(kTRUE),
     fDecayer(new AliDecayerPythia()),
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
     fUseNuclearPDF(kFALSE),
     fUseLorentzBoost(kTRUE),
     fNewMIS(kFALSE),
     fHFoff(kFALSE),
     fNucPdf(0),
     fTriggerParticle(0),
     fTriggerEta(0.9),
     fTriggerY(999.),
     fTriggerEtaMin(0.9),
     fTriggerMinPt(-1),
     fTriggerMaxPt(1000),
     fTriggerMultiplicity(0),
     fTriggerMultiplicityEta(0),
     fTriggerMultiplicityEtaMin(0),
     fTriggerMultiplicityEtaMax(0),
     fTriggerMultiplicityPtMin(0),
     fCountMode(kCountAll),
     fHeader(0),
     fRL(0),
     fkFileName(0),
     fkNameLHEF(0),
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
     fkTuneForDiff(0),
     fProcDiff(0)
{
// default charm production at 5. 5 TeV
// semimuonic decay
// structure function GRVHO
//
    fEnergyCMS = 5500.;
    fName = "Pythia";
    fTitle= "Particle Generator using PYTHIA";
    SetForceDecay();
    // Set random number generator
    if (!AliPythiaRndm::GetPythiaRandom())
      AliPythiaRndm::SetPythiaRandom(GetRandom());
 }

AliGenPythia::~AliGenPythia()
{
// Destructor
  if(fEventsTime) delete fEventsTime;
}

void AliGenPythia::SetInteractionRate(Float_t rate,Float_t timewindow)
{
// Generate pileup using user specified rate
    fInteractionRate = rate;
    fTimeWindow = timewindow;
    GeneratePileup();
}

void AliGenPythia::GeneratePileup()
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

void AliGenPythia::SetPycellParameters(Float_t etamax, Int_t neta, Int_t nphi,
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



void AliGenPythia::SetEventListRange(Int_t eventFirst, Int_t eventLast)
{
  // Set a range of event numbers, for which a table
  // of generated particle will be printed
  fDebugEventFirst = eventFirst;
  fDebugEventLast  = eventLast;
  if (fDebugEventLast==-1) fDebugEventLast=fDebugEventFirst;
}

void AliGenPythia::Init()
{
// Initialisation

    // Coeffs to go from mm / mm to meter / second
    SetGeneratorUnitsForMeterSecond(1.e-3, 1e-3/TMath::C());

    SetMC(AliPythia::Instance());
    fPythia=(AliPythia*) fMCEvGen;

//
    fParentWeight=1./Float_t(fNpart);
//
    if (fWeightPower != 0)
      fPythia->SetWeightPower(fWeightPower);
    fPythia->SetCKIN(3,fPtHardMin);
    fPythia->SetCKIN(4,fPtHardMax);
    fPythia->SetCKIN(7,fYHardMin);
    fPythia->SetCKIN(8,fYHardMax);

    if (fProjectile != "p" || fTarget != "p") fPythia->SetCollisionSystem(fProjectile,fTarget);

    if(fUseNuclearPDF)
      fPythia->SetNuclei(fAProjectile, fATarget, fNucPdf);
    // Fragmentation?
    if (fFragmentation) {
      fPythia->SetMSTP(111,1);
    } else {
      fPythia->SetMSTP(111,0);
    }


//  initial state radiation
    fPythia->SetMSTP(61,fGinit);
//  final state radiation
    fPythia->SetMSTP(71,fGfinal);
    //color reconnection strength
    if(fCRoff==1)fPythia->SetMSTP(95,0);
//  pt - kick
    if (fPtKick > 0.) {
	fPythia->SetMSTP(91,1);
	fPythia->SetPARP(91,fPtKick);
	fPythia->SetPARP(93, 4. * fPtKick);
    } else {
	fPythia->SetMSTP(91,0);
    }

	   if (fReadLHEF) fPythia->OpenFortranFile(97, const_cast<char*>(fkNameLHEF));

    if (fReadFromFile) {
	fRL  =  AliRunLoader::Open(fkFileName, "Partons");
	fRL->LoadKinematics();
	fRL->LoadHeader();
    } else {
	fRL = 0x0;
    }
 //
    fPythia->ProcInit(fProcess,fEnergyCMS,fStrucFunc, fItune);
    //  Forward Paramters to the AliPythia object
    fDecayer->SetForceDecay(fForceDecay);
// Switch off Heavy Flavors on request
    if (fHFoff) {
	// Maximum number of quark flavours used in pdf
	fPythia->SetMSTP(58, 3);
	// Maximum number of flavors that can be used in showers
	fPythia->SetMSTJ(45, 3);
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
    case kPyMbDefault:
    case kPyMbAtlasTuneMC09:
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
    fPythia->SetPARU(51, fPycellEtaMax);
    fPythia->SetMSTU(51, fPycellNEta);
    fPythia->SetMSTU(52, fPycellNPhi);
//
//  Configure Jet Finder
//
    fPythia->SetPARU(58,  fPycellThreshold);
    fPythia->SetPARU(52,  fPycellEtSeed);
    fPythia->SetPARU(53,  fPycellMinEtJet);
    fPythia->SetPARU(54,  fPycellMaxRadius);
    fPythia->SetMSTU(54,  2);
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

    // Reset Lorentz boost if demanded
    if(!fUseLorentzBoost) {
      fDyBoost = 0;
      Warning("Init","Demand to discard Lorentz boost.\n");
    }
//
//
//
    if (fSetNuclei) {
      fDyBoost = 0;
      Warning("Init","Deprecated function SetNuclei() used (nPDFs + no boost). Use SetProjectile + SetTarget + SetUseNuclearPDF + SetUseLorentzBoost instead.\n");
    }
    fPythia->SetPARJ(200, 0.0);
    fPythia->SetPARJ(199, 0.0);
    fPythia->SetPARJ(198, 0.0);
    fPythia->SetPARJ(197, 0.0);

    if (fQuench == 1) {
	fPythia->InitQuenching(0., 0.1, 0.6e6, 0);
    }

    if(fQuench ==2){fPythia->SetPyquenParameters(fpyquenT,fpyquenTau,fpyquenNf,fpyquenEloss,fpyquenAngle);}

    if (fQuench == 3) {
	// Nestor's change of the splittings
	fPythia->SetPARJ(200, 0.8);
	fPythia->SetMSTJ(41, 1);  // QCD radiation only
	fPythia->SetMSTJ(42, 2);  // angular ordering
	fPythia->SetMSTJ(44, 2);  // option to run alpha_s
	fPythia->SetMSTJ(47, 0);  // No correction back to hard scattering element
	fPythia->SetMSTJ(50, 0);  // No coherence in first branching
	fPythia->SetPARJ(82, 1.); // Cut off for parton showers
    } else if (fQuench == 4) {
	// Armesto-Cunqueiro-Salgado change of the splittings.
	AliFastGlauber* glauber = AliFastGlauber::Instance();
	glauber->Init(2);
	//read and store transverse almonds corresponding to differnt
	//impact parameters.
	glauber->SetCentralityClass(0.,0.1);
	fPythia->SetPARJ(200, 1.);
	fPythia->SetPARJ(198, fQhat);
	fPythia->SetPARJ(199, fLength);
	fPythia->SetMSTJ(42, 2);  // angular ordering
	fPythia->SetMSTJ(44, 2);  // option to run alpha_s
	fPythia->SetPARJ(82, 1.); // Cut off for parton showers
    }

  if ( AliLog::GetDebugLevel("","AliGenPythia") >= 1 ) {
    fPythia->Pystat(4);
    fPythia->Pystat(5);
  }
}

void AliGenPythia::Generate()
{
// Generate one event
    if (!fPythia) fPythia=(AliPythia*) fMCEvGen;
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
	fPythia->SetMSTJ(1, 0);

	if (fQuench ==4){
	    Double_t bimp;
	    // Quenching comes through medium-modified splitting functions.
	    AliFastGlauber::Instance()->GetRandomBHard(bimp);
	    fPythia->SetPARJ(197, bimp);
	    fImpact = bimp;
	    fPythia->Qpygin0();
	}
//
// Either produce new event or read partons from file
//
	if (!fReadFromFile) {
	    if (!fNewMIS) {
		fPythia->Pyevnt();
	    } else {
		fPythia->Pyevnw();
	    }
	    fNpartons = fPythia->GetN();
	} else {
	    printf("Loading Event %d\n",AliRunLoader::Instance()->GetEventNumber());
	    fRL->GetEvent(AliRunLoader::Instance()->GetEventNumber());
	    fPythia->SetN(0);
	    LoadEvent(fRL->Stack(), 0 , 1);
	    fPythia->Pyedit(21);
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
	if (fHadronisation) {
	    fPythia->SetMSTJ(1, 1);
//
// .. and perform hadronisation
//	printf("Calling hadronisation %d\n", fPythia->GetN());

	    if (fPatchOmegaDalitz) {
	      fPythia->SetMDCY(fPythia->Pycomp(111) ,1, 0);
	      fPythia->Pyexec();
	      fPythia->DalitzDecays();
	      fPythia->SetMDCY(fPythia->Pycomp(111) ,1, 1);
	    }

	    else  if (fDecayerExodus) {

              fPythia->SetMDCY(fPythia->Pycomp(22) ,1, 0);
              fPythia->SetMDCY(fPythia->Pycomp(111) ,1, 0);
              fPythia->SetMDCY(fPythia->Pycomp(221) ,1, 0);
              fPythia->SetMDCY(fPythia->Pycomp(223) ,1, 0);
              fPythia->Pyexec();
              fPythia->EtaprimeDalitz();
              fPythia->SetMDCY(fPythia->Pycomp(223) ,1, 1);
              fPythia->OmegaDalitz();
              fPythia->PhiDalitz();
              fPythia->SetMDCY(fPythia->Pycomp(111) ,1, 1);
              fPythia->PizeroDalitz();
              fPythia->SetMDCY(fPythia->Pycomp(221) ,1, 1);
              fPythia->EtaDalitz();
              fPythia->SetMDCY(fPythia->Pycomp(22) ,1, 1);
              fPythia->RhoDirect();
              fPythia->OmegaDirect();
              fPythia->PhiDirect();
              fPythia->JPsiDirect();

	    }

	    fPythia->Pyexec();
	}
	fTrials++;
	fPythia->ImportParticles(&fParticles,"All");

	if (TMath::Abs(fDyBoost) > 1.e-4) Boost();
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
	}

	Int_t nc = 0;        // Total n. of selected particles
	Int_t nParents = 0;  // Selected parents
	Int_t nTkbles = 0;   // Trackable particles


	if (fProcess != kPyMbDefault &&
	    fProcess != kPyMb &&
	    fProcess != kPyMbAtlasTuneMC09 &&
	    fProcess != kPyMbWithDirectPhoton &&
	    fProcess != kPyJets &&
	    fProcess != kPyDirectGamma &&
	    fProcess != kPyMbNonDiffr  &&
	    fProcess != kPyMbMSEL1     &&
	    fProcess != kPyW &&
	    fProcess != kPyZ &&
	    fProcess != kPyJpsi &&
            fProcess != kPyZgamma &&
	    fProcess != kPyCharmppMNRwmi &&
	    fProcess != kPyBeautyppMNRwmi &&
	    fProcess != kPyHeavyFlavppMNRwmi &&
	    fProcess != kPyBeautyJets &&
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
		Int_t ipa = iparticle->GetFirstMother()-1;
		Int_t kfMo = 0;
//
// Establish mother daughter relation between heavy quarks and mesons
//
		if (kf >= fFlavorSelect && kf <= 6) {
		    Int_t idau = iparticle->GetFirstDaughter() - 1;
		    if (idau > -1) {
			TParticle* daughter = (TParticle *) fParticles.At(idau);
			Int_t pdgD = daughter->GetPdgCode();
			if (pdgD == 91 || pdgD == 92) {
			    Int_t jmin = daughter->GetFirstDaughter() - 1;
			    Int_t jmax = daughter->GetLastDaughter()  - 1;
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
//		if ((ks == 1) || (fDecayer->GetLifetime(kf) > fMaxLifeTime))  trackIt[i] = 1;
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
		for (i = 0; i<np; i++) {
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
		    Int_t ipa     = iparticle->GetFirstMother()-1;
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

  ///---------------------------------------------------------------------------
  // Application of the user trigger
  if(!ApplyUserTrigger()) continue;
  ///---------------------------------------------------------------------------

	GetSubEventTime();

	delete[] pParent;
	delete[] pSelected;
	delete[] trackIt;

	if (nc > 0) {
	  switch (fCountMode) {
	  case kCountAll:
	    //printf(" Count all \n");
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

		fQ  += fPythia->GetVINT(51);
		fX1 += fPythia->GetVINT(41);
		fX2 += fPythia->GetVINT(42);
		fTrialsRun += fTrials;
		fNev++;
		MakeHeader();
		break;
	    }
	}
    } // event loop
    SetHighWaterMark(nt);
//  adjust weight due to kinematic selection
//    AdjustWeights();
//  get cross-section
    fXsection=fPythia->GetPARI(1);
}

Int_t  AliGenPythia::GenerateMB()
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
    const Float_t kconv=0.001/2.999792458e8;


    Int_t np = (fHadronisation) ? fParticles.GetEntriesFast() : fNpartons;



    Int_t* pParent = new Int_t[np];
    for (i=0; i< np; i++) pParent[i] = -1;
    if ((fProcess == kPyJets || fProcess == kPyDirectGamma || fProcess == kPyBeautyJets || fProcess == kPyBeautyppMNRwmi || fProcess == kPyJetsPWHG || fProcess == kPyCharmPWHG || fProcess == kPyBeautyPWHG)) {
 	TParticle* jet1 = (TParticle *) fParticles.At(6);
	TParticle* jet2 = (TParticle *) fParticles.At(7);

	if (!jet1 || ! jet2 || !CheckTrigger(jet1, jet2)) {
	  delete [] pParent;
	  return 0;
	}
    }


  // Select events with fragmentation photon, decay photon, pi0, eta or other hadrons going to PHOS or EMCAL or central barrel,
  // implemented primaryly for kPyJets, but extended to any kind of process.
  if ((fFragPhotonInCalo || fPi0InCalo || fEtaInCalo || fEleInCalo || fHadronInCalo || fDecayPhotonInCalo) &&
      (fCheckPHOS || fCheckEMCAL || fCheckBarrel || fCheckBarrelCalos) ) {
    Bool_t ok = TriggerOnSelectedParticles(np);

    if(!ok) {
      delete[] pParent;
      return 0;
    }
  }

    // Check for diffraction
    if(fkTuneForDiff) {
      if(fItune==320 && ( (TMath::Abs(fEnergyCMS - 900) < 1) || (TMath::Abs(fEnergyCMS - 2760) < 1) || (TMath::Abs(fEnergyCMS - 7000) < 1)) ) {
	if(!CheckDiffraction()) {
	  delete [] pParent;
	  return 0;
	}
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
	    if (TMath::Abs(iparticle->Y()) > fTriggerY) continue;
	    if (fTriggerEtaMin == fTriggerEta) {
	      if (TMath::Abs(iparticle->Eta()) > fTriggerEta) continue;
	    } else {
	      if (iparticle->Eta() > fTriggerEta || iparticle->Eta() < fTriggerEtaMin) continue;
	    }
	    if ( iparticle->Pt() > fTriggerMaxPt || iparticle->Pt() < fTriggerMinPt ) continue;
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

    if (fProcess == kPyCharmppMNRwmi || fProcess == kPyBeautyppMNRwmi || fProcess == kPyHeavyFlavppMNRwmi || fProcess == kPyBeautyJets) {
      TParticle *partCheck;
      TParticle *mother;
      Bool_t  theQ=kFALSE,theQbar=kFALSE,inYcut=kFALSE;
      Bool_t  theChild=kFALSE;
      Float_t y;
      Int_t   pdg, mpdg, mpdgUpperFamily;
      Bool_t  theHFFound=kFALSE;
      if(fCutOnChild==2){ //for Pairs
	theChild=kTRUE;
      }

      for(i = 0; i < np; i++) {
	partCheck = (TParticle*)fParticles.At(i);
	pdg = partCheck->GetPdgCode();
	Bool_t flavSel=kFALSE;
	if(TMath::Abs(pdg) == fFlavorSelect) flavSel=kTRUE;
	if(fProcess == kPyHeavyFlavppMNRwmi && (TMath::Abs(pdg) == 4 || TMath::Abs(pdg) == 5)) flavSel=kTRUE;
	if(flavSel){ // quark
	  if(pdg>0) { theQ=kTRUE; } else { theQbar=kTRUE; }
	  y = 0.5*TMath::Log((partCheck->Energy()+partCheck->Pz()+1.e-13)/
			     (partCheck->Energy()-partCheck->Pz()+1.e-13));
	  if(fUseYCutHQ && y>fYMinHQ && y<fYMaxHQ) inYcut=kTRUE;
	  if(!fUseYCutHQ && y>fYMin && y<fYMax) inYcut=kTRUE;
	}
	if(fCutOnChild==1 && TMath::Abs(pdg) == fPdgCodeParticleforAcceptanceCut) {
	  Int_t mi = partCheck->GetFirstMother() - 1;
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

	/// for pairs
	if(fCutOnChild==2 && TMath::Abs(pdg) == fPdgCodeParticleforAcceptanceCut) {
	  Int_t mi = partCheck->GetFirstMother() - 1;
	  if(mi<0) continue;
	  mother = (TParticle*)fParticles.At(mi);
	  mpdg=TMath::Abs(mother->GetPdgCode());
	  //mpdgUpperFamily=(mpdg>1000 ? mpdg+1000 : mpdg+100); // keep e from c from b
	  //above is not necesary. we would like to keep bb->ee in KinematicSelection.

	  //if ( ParentSelected(mpdg) ||
	  //(fFlavorSelect==5 && ParentSelected(mpdgUpperFamily))) {
	  if(ParentSelected(mpdg)){
	    theHFFound = kTRUE; // found heavy flavors
	    if (!KinematicSelection(partCheck,1)) {
	      theChild=kFALSE;
	    }else{
	      printf(" ********** found leg in acceptance... %d %d %d %f %f\n",
		     partCheck->GetPdgCode(), mpdg, i, partCheck->Y(), partCheck->Pt());
	    }
	  }
	}
      }
      if (!theQ || !theQbar || !inYcut) { // one of the c/b conditions not satisfied
	delete[] pParent;
	return 0;
      }
      if (fCutOnChild==1 && !theChild) { // one of the child conditions not satisfied
	delete[] pParent;
	return 0;
      }
      if (fCutOnChild==2 && !(theHFFound && theChild)) { // both of the child conditions not satisfied
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
	  fProcess == kPyMbAtlasTuneMC09 ||
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
	if (
	    (((ks == 1  && kf!=0 && KinematicSelection(iparticle, 0)) || (ks !=1)) && ((fStackFillOpt != kHeavyFlavor) || IsFromHeavyFlavor(i))) ||
	    ((fProcess == kPyJets || fProcess == kPyBeautyJets || fProcess == kPyJetsPWHG || fProcess == kPyCharmPWHG || fProcess == kPyBeautyPWHG) && ks == 21 && km == 0 && i>1)
	    )
	  {
            nc++;
	    if (ks == 1) trackIt = 1;
	    Int_t ipa = iparticle->GetFirstMother()-1;

	    iparent = (ipa > -1) ? pParent[ipa] : -1;

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
	    KeepTrack(nt);
	    pParent[i] = nt;
	    SetHighWaterMark(nt);

	} // select particle
    } // particle loop

    delete[] pParent;

    return 1;
}


void AliGenPythia::FinishRun()
{
// Print x-section summary
    fPythia->Pystat(1);

    if (fNev > 0.) {
	fQ  /= fNev;
	fX1 /= fNev;
	fX2 /= fNev;
    }

    printf("\nTotal number of Pyevnt() calls %d\n", fTrialsRun);
    printf("\nMean Q, x1, x2: %f %f %f\n", fQ, fX1, fX2);

    WriteXsection();
}

void AliGenPythia::AdjustWeights() const
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

void AliGenPythia::SetNuclei(Int_t a1, Int_t a2, Int_t pdfset)
{
// Treat protons as inside nuclei with mass numbers a1 and a2

    fAProjectile   = a1;
    fATarget       = a2;
    fNucPdf        = pdfset;  // 0 EKS98 9 EPS09LO 19 EPS09NLO
    fUseNuclearPDF = kTRUE;
    fSetNuclei     = kTRUE;
}


void AliGenPythia::MakeHeader()
{
//
// Make header for the simulated event
//
  if (gAlice) {
    if (gAlice->GetEvNumber()>=fDebugEventFirst &&
	gAlice->GetEvNumber()<=fDebugEventLast) fPythia->Pylist(2);
  }

// Builds the event header, to be called after each event
    if (fHeader) delete fHeader;
    fHeader = new AliGenPythiaEventHeader("Pythia");
    fHeader->SetTitle(GetTitle());
//
// Event type
    if(fProcDiff>0){
      //      if(fProcDiff == 92 || fProcDiff == 93) printf("\n\n\n\n\n");
      //      printf("fPythia->GetMSTI(1) = %d   fProcDiff = %d\n",fPythia->GetMSTI(1), fProcDiff);
    ((AliGenPythiaEventHeader*) fHeader)->SetProcessType(fProcDiff);
    }
    else
    ((AliGenPythiaEventHeader*) fHeader)->SetProcessType(fPythia->GetMSTI(1));
//
// Number of trials
    ((AliGenPythiaEventHeader*) fHeader)->SetTrials(fTrials);
//
// Event Vertex
    fHeader->SetPrimaryVertex(fVertex);
    fHeader->SetInteractionTime(fTime+fEventTime);
//
// Number of primaries
    fHeader->SetNProduced(fNprimaries);
//
// Event weight
    fHeader->SetEventWeight(fPythia->GetVINT(97));
//
// Number of MPI
    fHeader->SetNMPI(fPythia->GetNMPI());
//
// Jets that have triggered

    //Need to store jets for b-jet studies too!
    if (fProcess == kPyJets || fProcess == kPyDirectGamma || fProcess == kPyBeautyJets || fProcess == kPyBeautyppMNRwmi || fProcess == kPyJetsPWHG || fProcess == kPyCharmPWHG || fProcess == kPyBeautyPWHG)
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
	} else if (fQuench == 2){
	    // Pyquen
	    Double_t r1 = PARIMP.rb1;
	    Double_t r2 = PARIMP.rb2;
	    Double_t b  = PARIMP.b1;
	    Double_t r   = 0.5 * TMath::Sqrt(2. * (r1 * r1 + r2 * r2) - b * b);
	    Double_t phi = PARIMP.psib1;
	    xp = r * TMath::Cos(phi);
	    yp = r * TMath::Sin(phi);

	} else if (fQuench == 4) {
	    // QPythia
	    Double_t xy[2];
	    Double_t i0i1[2];
	    AliFastGlauber::Instance()->GetSavedXY(xy);
	    AliFastGlauber::Instance()->GetSavedI0I1(i0i1);
	    xp = xy[0];
	    yp = xy[1];
	    ((AliGenPythiaEventHeader*) fHeader)->SetImpactParameter(fImpact);
	}

	    ((AliGenPythiaEventHeader*) fHeader)->SetXYJet(xp, yp);
	    ((AliGenPythiaEventHeader*) fHeader)->SetZQuench(z);
    }
//
// Store pt^hard and cross-section
    ((AliGenPythiaEventHeader*) fHeader)->SetPtHard(fPythia->GetVINT(47));
    ((AliGenPythiaEventHeader*) fHeader)->SetXsection(fPythia->GetPARI(1));

//
// Store Event Weight
    ((AliGenPythiaEventHeader*) fHeader)->SetEventWeight(fPythia->GetPARI(7)*fPythia->GetPARI(10));
    // PARI(7) is 1 or -1, for weighted generation with accept/reject, e.g. POWHEG
    // PARI(10) is a weight associated with reweighted generation, using Pyevwt
//
//  Pass header
//
    AddHeader(fHeader);
    fHeader = 0x0;
}

Bool_t AliGenPythia::CheckTrigger(const TParticle* jet1, const TParticle* jet2)
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

    if (fProcess == kPyJets || fProcess == kPyBeautyJets || fProcess ==  kPyBeautyppMNRwmi || fProcess == kPyJetsPWHG || fProcess == kPyCharmPWHG || fProcess == kPyBeautyPWHG) {
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



Bool_t AliGenPythia::CheckKinematicsOnChild(){
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

void  AliGenPythia::LoadEvent(AliStack* stack, Int_t flag, Int_t reHadr)
{
  //
  // Load event into Pythia Common Block
  //

  Int_t npart = stack -> GetNprimary();
  Int_t n0 = 0;

  if (!flag) {
    (fPythia->GetPyjets())->N = npart;
  } else {
    n0 = (fPythia->GetPyjets())->N;
    (fPythia->GetPyjets())->N = n0 + npart;
  }


  for (Int_t part = 0; part < npart; part++) {
    TParticle *mPart = stack->Particle(part);

    Int_t kf     =  mPart->GetPdgCode();
    Int_t ks     =  mPart->GetStatusCode();
    Int_t idf    =  mPart->GetFirstDaughter();
    Int_t idl    =  mPart->GetLastDaughter();

    if (reHadr) {
	    if (ks == 11 || ks == 12) {
        ks  -= 10;
        idf  = -1;
        idl  = -1;
	    }
    }

    Float_t px = mPart->Px();
    Float_t py = mPart->Py();
    Float_t pz = mPart->Pz();
    Float_t e  = mPart->Energy();
    Float_t m  = mPart->GetCalcMass();


    (fPythia->GetPyjets())->P[0][part+n0] = px;
    (fPythia->GetPyjets())->P[1][part+n0] = py;
    (fPythia->GetPyjets())->P[2][part+n0] = pz;
    (fPythia->GetPyjets())->P[3][part+n0] = e;
    (fPythia->GetPyjets())->P[4][part+n0] = m;

    (fPythia->GetPyjets())->K[1][part+n0] = kf;
    (fPythia->GetPyjets())->K[0][part+n0] = ks;
    (fPythia->GetPyjets())->K[3][part+n0] = idf + 1;
    (fPythia->GetPyjets())->K[4][part+n0] = idl + 1;
    (fPythia->GetPyjets())->K[2][part+n0] = mPart->GetFirstMother() + 1;
  }
}

void  AliGenPythia::LoadEvent(const TObjArray* stack, Int_t flag, Int_t reHadr)
{
  //
  // Load event into Pythia Common Block
  //

  Int_t npart = stack -> GetEntries();
  Int_t n0 = 0;

  if (!flag) {
    (fPythia->GetPyjets())->N = npart;
  } else {
    n0 = (fPythia->GetPyjets())->N;
    (fPythia->GetPyjets())->N = n0 + npart;
  }


  for (Int_t part = 0; part < npart; part++) {
    TParticle *mPart = dynamic_cast<TParticle *>(stack->At(part));
    if (!mPart) continue;

    Int_t kf     =  mPart->GetPdgCode();
    Int_t ks     =  mPart->GetStatusCode();
    Int_t idf    =  mPart->GetFirstDaughter();
    Int_t idl    =  mPart->GetLastDaughter();

    if (reHadr) {
	if (ks == 11 || ks == 12) {
	    ks  -= 10;
	    idf  = -1;
	    idl  = -1;
	}
    }

    Float_t px = mPart->Px();
    Float_t py = mPart->Py();
    Float_t pz = mPart->Pz();
    Float_t e  = mPart->Energy();
    Float_t m  = mPart->GetCalcMass();


    (fPythia->GetPyjets())->P[0][part+n0] = px;
    (fPythia->GetPyjets())->P[1][part+n0] = py;
    (fPythia->GetPyjets())->P[2][part+n0] = pz;
    (fPythia->GetPyjets())->P[3][part+n0] = e;
    (fPythia->GetPyjets())->P[4][part+n0] = m;

    (fPythia->GetPyjets())->K[1][part+n0] = kf;
    (fPythia->GetPyjets())->K[0][part+n0] = ks;
    (fPythia->GetPyjets())->K[3][part+n0] = idf + 1;
    (fPythia->GetPyjets())->K[4][part+n0] = idl + 1;
    (fPythia->GetPyjets())->K[2][part+n0] = mPart->GetFirstMother() + 1;
  }
}


void AliGenPythia::RecJetsUA1(Int_t& njets, Float_t jets [4][50])
{
//
//  Calls the Pythia jet finding algorithm to find jets in the current event
//
//
//
//  Save jets
    Int_t n     = fPythia->GetN();

//
//  Run Jet Finder
    fPythia->Pycell(njets);
    Int_t i;
    for (i = 0; i < njets; i++) {
	Float_t px    = (fPythia->GetPyjets())->P[0][n+i];
	Float_t py    = (fPythia->GetPyjets())->P[1][n+i];
	Float_t pz    = (fPythia->GetPyjets())->P[2][n+i];
	Float_t e     = (fPythia->GetPyjets())->P[3][n+i];

	jets[0][i] = px;
	jets[1][i] = py;
	jets[2][i] = pz;
	jets[3][i] = e;
    }
}



void  AliGenPythia::GetJets(Int_t& nJets, Int_t& nJetsTrig, Float_t jets[4][10])
{
//
//  Calls the Pythia clustering algorithm to find jets in the current event
//
    Int_t n     = fPythia->GetN();
    nJets       = 0;
    nJetsTrig   = 0;
    if (fJetReconstruction == kCluster) {
//
//  Configure cluster algorithm
//
	fPythia->SetPARU(43, 2.);
	fPythia->SetMSTU(41, 1);
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
	Float_t px    = (fPythia->GetPyjets())->P[0][n+i];
	Float_t py    = (fPythia->GetPyjets())->P[1][n+i];
	Float_t pz    = (fPythia->GetPyjets())->P[2][n+i];
	Float_t e     = (fPythia->GetPyjets())->P[3][n+i];
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
//	    printf("\n........-Jet #%d: %10.3f %10.3f %10.3f %10.3f \n", i, pt, et, eta, phi * kRaddeg);
	} else {
//	    printf("\n........-Jet #%d: %10.3f %10.3f %10.3f %10.3f \n", i, pt, et, eta, phi * kRaddeg);
	}
    }
}

void AliGenPythia::GetSubEventTime()
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

Bool_t AliGenPythia::IsInBarrel(Float_t eta) const
{
  // Is particle in Central Barrel acceptance?
  // etamin=-etamax
  if( eta < fTriggerEta  )
    return kTRUE;
  else
    return kFALSE;
}

Bool_t AliGenPythia::IsInBarrelCalorimeters(Float_t phi, Float_t eta)
{
  if      ( IsInEMCAL(phi,eta   ) ) return kTRUE ;
  else if ( IsInDCAL (phi,eta   ) ) return kTRUE ;
  else if ( IsInPHOS (phi,eta,-1) ) return kTRUE ;
  else                              return kFALSE;
}

Bool_t AliGenPythia::IsInEMCAL(Float_t phi, Float_t eta) const
{
  // Is particle in EMCAL acceptance?
  // phi in degrees, etamin=-etamax
  if(phi > fEMCALMinPhi  && phi < fEMCALMaxPhi &&
     eta < fEMCALEta  )
    return kTRUE;
  else
    return kFALSE;
}


Bool_t AliGenPythia::IsInDCAL(Float_t phi, Float_t eta) const
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

Bool_t AliGenPythia::IsInPHOS(Float_t phi, Float_t eta, Int_t iparticle)
{
  // Is particle in PHOS acceptance?
  // Acceptance slightly larger considered.
  // phi in degrees, etamin=-etamax
  // iparticle is the index of the particle to be checked, for PHOS rotation case

  if(phi > fPHOSMinPhi  && phi < fPHOSMaxPhi &&
     eta < fPHOSEta  )
    return kTRUE;
  else
  {
    if( fCheckPHOSeta && eta < fPHOSEta) fPHOSRotateCandidate = iparticle;

    return kFALSE;
  }
}

void AliGenPythia::RotatePhi(Bool_t& okdd)
{
  //Rotate event in phi to enhance events in PHOS acceptance

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
  for(Int_t i=0; i< np; i++) {
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


Bool_t AliGenPythia::CheckDiffraction()
{
  // use this method only with Perugia-0 tune!

  //  printf("AAA\n");

   Int_t np = (fHadronisation) ? fParticles.GetEntriesFast() : fNpartons;

   Int_t iPart1=-1;
   Int_t iPart2=-1;

   Double_t y1 = 1e10;
   Double_t y2 = -1e10;

  const Int_t kNstable=20;
  const Int_t pdgStable[20] = {
    22,             // Photon
    11,             // Electron
    12,             // Electron Neutrino
    13,             // Muon
    14,             // Muon Neutrino
    15,             // Tau
    16,             // Tau Neutrino
    211,            // Pion
    321,            // Kaon
    311,            // K0
    130,            // K0s
    310,            // K0l
    2212,           // Proton
    2112,           // Neutron
    3122,           // Lambda_0
    3112,           // Sigma Minus
    3222,           // Sigma Plus
    3312,           // Xsi Minus
    3322,           // Xsi0
    3334            // Omega
  };

     for (Int_t i = 0; i < np; i++) {
	TParticle *  part = (TParticle *) fParticles.At(i);

	Int_t statusCode = part->GetStatusCode();

	// Initial state particle
	if (statusCode != 1)
	  continue;

	Int_t pdg = TMath::Abs(part->GetPdgCode());
	Bool_t isStable = kFALSE;
	for (Int_t i1 = 0; i1 < kNstable; i1++) {
	  if (pdg == pdgStable[i1]) {
	    isStable = kTRUE;
	    break;
	  }
	}
	if(!isStable)
	  continue;

	Double_t y = part->Y();

	if (y < y1)
	  {
	    y1 = y;
	    iPart1 = i;
	  }
	if (y > y2)
	{
	  y2 = y;
	  iPart2 = i;
	}
     }

     if(iPart1<0 || iPart2<0) return kFALSE;

     y1=TMath::Abs(y1);
     y2=TMath::Abs(y2);

     TParticle *  part1 = (TParticle *) fParticles.At(iPart1);
     TParticle *  part2 = (TParticle *) fParticles.At(iPart2);

     Int_t pdg1 = part1->GetPdgCode();
     Int_t pdg2 = part2->GetPdgCode();


     Int_t iPart = -1;
     if (pdg1 == 2212 && pdg2 == 2212)
       {
	 if(y1 > y2)
	   iPart = iPart1;
	 else if(y1 < y2)
	   iPart = iPart2;
	 else {
	   iPart = iPart1;
	   if((AliPythiaRndm::GetPythiaRandom())->Uniform(0.,1.)>0.5) iPart = iPart2;
	 }
       }
     else if (pdg1 == 2212)
       iPart = iPart1;
     else if (pdg2 == 2212)
       iPart = iPart2;





     Double_t M=-1.;
     if(iPart>0) {
       TParticle *  part = (TParticle *) fParticles.At(iPart);
       Double_t E= part->Energy();
       Double_t P= part->P();
       Double_t M2 = (fEnergyCMS-E-P)*(fEnergyCMS-E+P);
       if(M2<0)  return kFALSE;
       M= TMath::Sqrt(M2);
     }

     Double_t Mmin, Mmax, wSD, wDD, wND;
     if(!GetWeightsDiffraction(M, Mmin, Mmax, wSD, wDD, wND)) return kFALSE;

     if(M>-1 && M<Mmin) return kFALSE;
     if(M>Mmax) M=-1;

     Int_t procType=fPythia->GetMSTI(1);
     Int_t proc0=2;
     if(procType== 94) proc0=1;
     if(procType== 92 || procType== 93) proc0=0;

     Int_t proc=2;
     if(M>0) proc=0;
     else if(proc0==1) proc=1;

     if(proc==1 && (AliPythiaRndm::GetPythiaRandom())->Uniform(0.,1.) > wDD) return kFALSE;
     if(proc==2 && (AliPythiaRndm::GetPythiaRandom())->Uniform(0.,1.) > wND) return kFALSE;


    //     if(proc==1 || proc==2) return kFALSE;

     if(proc!=0) {
       if(proc0!=0) fProcDiff = procType;
       else       fProcDiff = 95;
       return kTRUE;
     }

    if(wSD<0)  AliError("wSD<0 ! \n");

    if((AliPythiaRndm::GetPythiaRandom())->Uniform(0.,1.)> wSD) return kFALSE;

    //    printf("iPart = %d\n", iPart);

    if(iPart==iPart1) fProcDiff=93;
    else if(iPart==iPart2) fProcDiff=92;
    else {
      printf("EROOR:  iPart!=iPart1 && iPart!=iPart2\n");

    }

    return kTRUE;
}



Bool_t AliGenPythia::GetWeightsDiffraction(Double_t M, Double_t &Mmin, Double_t &Mmax,
                                                       Double_t &wSD, Double_t &wDD, Double_t &wND)
{

  // 900 GeV
  if(TMath::Abs(fEnergyCMS-900)<1 ){

const Int_t nbin=400;
Double_t bin[]={
1.080000, 1.577300, 2.074600, 2.571900, 3.069200, 3.566500,
4.063800, 4.561100, 5.058400, 5.555700, 6.053000, 6.550300,
7.047600, 7.544900, 8.042200, 8.539500, 9.036800, 9.534100,
10.031400, 10.528700, 11.026000, 11.523300, 12.020600, 12.517900,
13.015200, 13.512500, 14.009800, 14.507100, 15.004400, 15.501700,
15.999000, 16.496300, 16.993600, 17.490900, 17.988200, 18.485500,
18.982800, 19.480100, 19.977400, 20.474700, 20.972000, 21.469300,
21.966600, 22.463900, 22.961200, 23.458500, 23.955800, 24.453100,
24.950400, 25.447700, 25.945000, 26.442300, 26.939600, 27.436900,
27.934200, 28.431500, 28.928800, 29.426100, 29.923400, 30.420700,
30.918000, 31.415300, 31.912600, 32.409900, 32.907200, 33.404500,
33.901800, 34.399100, 34.896400, 35.393700, 35.891000, 36.388300,
36.885600, 37.382900, 37.880200, 38.377500, 38.874800, 39.372100,
39.869400, 40.366700, 40.864000, 41.361300, 41.858600, 42.355900,
42.853200, 43.350500, 43.847800, 44.345100, 44.842400, 45.339700,
45.837000, 46.334300, 46.831600, 47.328900, 47.826200, 48.323500,
48.820800, 49.318100, 49.815400, 50.312700, 50.810000, 51.307300,
51.804600, 52.301900, 52.799200, 53.296500, 53.793800, 54.291100,
54.788400, 55.285700, 55.783000, 56.280300, 56.777600, 57.274900,
57.772200, 58.269500, 58.766800, 59.264100, 59.761400, 60.258700,
60.756000, 61.253300, 61.750600, 62.247900, 62.745200, 63.242500,
63.739800, 64.237100, 64.734400, 65.231700, 65.729000, 66.226300,
66.723600, 67.220900, 67.718200, 68.215500, 68.712800, 69.210100,
69.707400, 70.204700, 70.702000, 71.199300, 71.696600, 72.193900,
72.691200, 73.188500, 73.685800, 74.183100, 74.680400, 75.177700,
75.675000, 76.172300, 76.669600, 77.166900, 77.664200, 78.161500,
78.658800, 79.156100, 79.653400, 80.150700, 80.648000, 81.145300,
81.642600, 82.139900, 82.637200, 83.134500, 83.631800, 84.129100,
84.626400, 85.123700, 85.621000, 86.118300, 86.615600, 87.112900,
87.610200, 88.107500, 88.604800, 89.102100, 89.599400, 90.096700,
90.594000, 91.091300, 91.588600, 92.085900, 92.583200, 93.080500,
93.577800, 94.075100, 94.572400, 95.069700, 95.567000, 96.064300,
96.561600, 97.058900, 97.556200, 98.053500, 98.550800, 99.048100,
99.545400, 100.042700, 100.540000, 101.037300, 101.534600, 102.031900,
102.529200, 103.026500, 103.523800, 104.021100, 104.518400, 105.015700,
105.513000, 106.010300, 106.507600, 107.004900, 107.502200, 107.999500,
108.496800, 108.994100, 109.491400, 109.988700, 110.486000, 110.983300,
111.480600, 111.977900, 112.475200, 112.972500, 113.469800, 113.967100,
114.464400, 114.961700, 115.459000, 115.956300, 116.453600, 116.950900,
117.448200, 117.945500, 118.442800, 118.940100, 119.437400, 119.934700,
120.432000, 120.929300, 121.426600, 121.923900, 122.421200, 122.918500,
123.415800, 123.913100, 124.410400, 124.907700, 125.405000, 125.902300,
126.399600, 126.896900, 127.394200, 127.891500, 128.388800, 128.886100,
129.383400, 129.880700, 130.378000, 130.875300, 131.372600, 131.869900,
132.367200, 132.864500, 133.361800, 133.859100, 134.356400, 134.853700,
135.351000, 135.848300, 136.345600, 136.842900, 137.340200, 137.837500,
138.334800, 138.832100, 139.329400, 139.826700, 140.324000, 140.821300,
141.318600, 141.815900, 142.313200, 142.810500, 143.307800, 143.805100,
144.302400, 144.799700, 145.297000, 145.794300, 146.291600, 146.788900,
147.286200, 147.783500, 148.280800, 148.778100, 149.275400, 149.772700,
150.270000, 150.767300, 151.264600, 151.761900, 152.259200, 152.756500,
153.253800, 153.751100, 154.248400, 154.745700, 155.243000, 155.740300,
156.237600, 156.734900, 157.232200, 157.729500, 158.226800, 158.724100,
159.221400, 159.718700, 160.216000, 160.713300, 161.210600, 161.707900,
162.205200, 162.702500, 163.199800, 163.697100, 164.194400, 164.691700,
165.189000, 165.686300, 166.183600, 166.680900, 167.178200, 167.675500,
168.172800, 168.670100, 169.167400, 169.664700, 170.162000, 170.659300,
171.156600, 171.653900, 172.151200, 172.648500, 173.145800, 173.643100,
174.140400, 174.637700, 175.135000, 175.632300, 176.129600, 176.626900,
177.124200, 177.621500, 178.118800, 178.616100, 179.113400, 179.610700,
180.108000, 180.605300, 181.102600, 181.599900, 182.097200, 182.594500,
183.091800, 183.589100, 184.086400, 184.583700, 185.081000, 185.578300,
186.075600, 186.572900, 187.070200, 187.567500, 188.064800, 188.562100,
189.059400, 189.556700, 190.054000, 190.551300, 191.048600, 191.545900,
192.043200, 192.540500, 193.037800, 193.535100, 194.032400, 194.529700,
195.027000, 195.524300, 196.021600, 196.518900, 197.016200, 197.513500,
198.010800, 198.508100, 199.005400, 199.502700, 200.000000};
Double_t w[]={
1.000000, 0.643457, 0.645609, 0.648347, 0.604563, 0.605002,
0.602819, 0.611473, 0.576412, 0.562354, 0.550216, 0.529285,
0.534558, 0.534364, 0.530358, 0.518475, 0.489253, 0.469754,
0.469825, 0.450513, 0.455849, 0.435312, 0.437210, 0.456686,
0.413577, 0.427093, 0.426894, 0.418834, 0.409475, 0.388483,
0.412208, 0.388912, 0.389611, 0.382572, 0.389220, 0.370964,
0.380463, 0.370873, 0.363701, 0.369363, 0.357361, 0.365759,
0.348566, 0.337062, 0.348190, 0.332330, 0.359001, 0.335836,
0.339154, 0.335599, 0.336035, 0.335204, 0.353440, 0.337836,
0.333874, 0.307120, 0.294963, 0.324978, 0.313359, 0.317688,
0.323758, 0.319304, 0.335317, 0.301765, 0.317257, 0.356331,
0.323119, 0.297732, 0.303188, 0.326102, 0.316467, 0.294728,
0.308135, 0.288857, 0.325692, 0.312493, 0.291100, 0.325921,
0.313317, 0.295980, 0.308481, 0.328380, 0.313081, 0.296763,
0.295431, 0.317325, 0.320462, 0.286918, 0.316035, 0.335208,
0.283193, 0.333945, 0.292534, 0.294164, 0.330680, 0.296992,
0.285509, 0.317260, 0.298269, 0.311299, 0.312129, 0.286822,
0.287442, 0.319139, 0.283314, 0.318454, 0.297727, 0.301597,
0.282483, 0.294792, 0.305569, 0.290957, 0.297817, 0.282908,
0.272401, 0.305584, 0.300220, 0.297020, 0.298781, 0.278008,
0.277727, 0.323777, 0.287419, 0.342074, 0.287259, 0.303658,
0.302668, 0.279622, 0.280586, 0.313630, 0.276068, 0.257051,
0.309996, 0.265534, 0.297138, 0.281738, 0.294610, 0.292882,
0.286860, 0.312686, 0.293244, 0.293744, 0.271375, 0.278734,
0.280308, 0.304739, 0.287907, 0.285261, 0.311180, 0.313476,
0.289660, 0.289282, 0.319505, 0.271285, 0.272008, 0.289245,
0.281038, 0.285284, 0.295836, 0.281416, 0.283501, 0.295417,
0.304372, 0.297764, 0.291378, 0.321530, 0.315604, 0.329507,
0.282609, 0.275576, 0.283721, 0.311714, 0.283760, 0.273188,
0.312193, 0.264347, 0.281532, 0.301226, 0.281718, 0.336408,
0.283157, 0.332010, 0.289974, 0.290256, 0.301569, 0.332228,
0.288282, 0.326339, 0.313653, 0.300361, 0.289470, 0.264830,
0.298659, 0.272359, 0.278878, 0.306001, 0.328168, 0.294991,
0.327737, 0.278056, 0.302435, 0.284183, 0.279270, 0.307279,
0.307917, 0.315196, 0.283803, 0.313333, 0.315730, 0.304818,
0.307171, 0.295223, 0.333741, 0.346911, 0.310143, 0.336686,
0.275459, 0.334781, 0.295405, 0.275816, 0.301678, 0.327242,
0.320717, 0.309230, 0.292145, 0.294489, 0.305088, 0.300969,
0.277438, 0.326159, 0.297065, 0.301177, 0.303843, 0.275382,
0.304019, 0.284166, 0.289610, 0.331611, 0.317131, 0.310880,
0.360456, 0.294052, 0.342694, 0.327166, 0.336797, 0.298040,
0.295767, 0.260053, 0.325544, 0.335310, 0.320182, 0.301072,
0.313117, 0.283407, 0.299206, 0.293525, 0.305067, 0.255978,
0.327055, 0.316382, 0.317700, 0.278993, 0.283120, 0.314000,
0.274396, 0.291208, 0.348813, 0.319603, 0.313076, 0.289155,
0.343988, 0.311426, 0.322896, 0.328726, 0.337062, 0.389093,
0.284122, 0.312184, 0.304008, 0.319170, 0.320778, 0.288390,
0.337272, 0.356273, 0.343310, 0.312209, 0.330709, 0.297977,
0.346146, 0.369162, 0.324385, 0.339831, 0.337037, 0.318739,
0.343157, 0.277720, 0.368407, 0.321330, 0.338997, 0.314220,
0.328861, 0.321824, 0.328013, 0.356925, 0.359144, 0.296314,
0.345415, 0.396711, 0.347032, 0.294928, 0.343799, 0.322331,
0.328656, 0.326098, 0.337338, 0.337038, 0.300179, 0.351391,
0.324337, 0.330896, 0.302842, 0.310522, 0.337052, 0.359989,
0.383250, 0.359355, 0.315382, 0.333113, 0.342598, 0.355348,
0.320751, 0.320475, 0.351762, 0.351475, 0.338358, 0.326153,
0.302507, 0.340048, 0.318685, 0.381646, 0.339320, 0.299453,
0.426599, 0.393515, 0.353929, 0.328435, 0.413976, 0.292558,
0.379340, 0.358344, 0.409259, 0.313821, 0.336675, 0.324521,
0.408382, 0.346273, 0.312939, 0.362453, 0.343152, 0.330577,
0.332831, 0.353299, 0.347745, 0.334818, 0.332234, 0.385585,
0.395483, 0.395316, 0.326972, 0.349434, 0.317604, 0.328980,
0.375056, 0.317290, 0.357083, 0.346165, 0.310444, 0.356873,
0.359523, 0.365308, 0.365122, 0.383685, 0.370975, 0.396928,
0.407654, 0.307755, 0.323033, 0.350580, 0.345231, 0.342462,
0.400000, 0.318309, 0.403570, 0.322856, 0.383053, 0.422252,
0.386112, 0.364314, 0.434375, 0.334629};
wDD = 0.379611;
wND = 0.496961;
wSD = -1;

    Mmin = bin[0];
    Mmax = bin[nbin];
    if(M<Mmin || M>Mmax) return kTRUE;

    Int_t ibin=nbin-1;
    for(Int_t i=1; i<=nbin; i++)
      if(M<=bin[i]) {
	ibin=i-1;
	//	printf("Mi> %f && Mi< %f\n", bin[i-1], bin[i]);
	break;
      }
    wSD=w[ibin];
    return kTRUE;
  }
 else if(TMath::Abs(fEnergyCMS-2760)<1 ){

const Int_t nbin=400;
Double_t bin[]={
1.080000, 1.577300, 2.074600, 2.571900, 3.069200, 3.566500,
4.063800, 4.561100, 5.058400, 5.555700, 6.053000, 6.550300,
7.047600, 7.544900, 8.042200, 8.539500, 9.036800, 9.534100,
10.031400, 10.528700, 11.026000, 11.523300, 12.020600, 12.517900,
13.015200, 13.512500, 14.009800, 14.507100, 15.004400, 15.501700,
15.999000, 16.496300, 16.993600, 17.490900, 17.988200, 18.485500,
18.982800, 19.480100, 19.977400, 20.474700, 20.972000, 21.469300,
21.966600, 22.463900, 22.961200, 23.458500, 23.955800, 24.453100,
24.950400, 25.447700, 25.945000, 26.442300, 26.939600, 27.436900,
27.934200, 28.431500, 28.928800, 29.426100, 29.923400, 30.420700,
30.918000, 31.415300, 31.912600, 32.409900, 32.907200, 33.404500,
33.901800, 34.399100, 34.896400, 35.393700, 35.891000, 36.388300,
36.885600, 37.382900, 37.880200, 38.377500, 38.874800, 39.372100,
39.869400, 40.366700, 40.864000, 41.361300, 41.858600, 42.355900,
42.853200, 43.350500, 43.847800, 44.345100, 44.842400, 45.339700,
45.837000, 46.334300, 46.831600, 47.328900, 47.826200, 48.323500,
48.820800, 49.318100, 49.815400, 50.312700, 50.810000, 51.307300,
51.804600, 52.301900, 52.799200, 53.296500, 53.793800, 54.291100,
54.788400, 55.285700, 55.783000, 56.280300, 56.777600, 57.274900,
57.772200, 58.269500, 58.766800, 59.264100, 59.761400, 60.258700,
60.756000, 61.253300, 61.750600, 62.247900, 62.745200, 63.242500,
63.739800, 64.237100, 64.734400, 65.231700, 65.729000, 66.226300,
66.723600, 67.220900, 67.718200, 68.215500, 68.712800, 69.210100,
69.707400, 70.204700, 70.702000, 71.199300, 71.696600, 72.193900,
72.691200, 73.188500, 73.685800, 74.183100, 74.680400, 75.177700,
75.675000, 76.172300, 76.669600, 77.166900, 77.664200, 78.161500,
78.658800, 79.156100, 79.653400, 80.150700, 80.648000, 81.145300,
81.642600, 82.139900, 82.637200, 83.134500, 83.631800, 84.129100,
84.626400, 85.123700, 85.621000, 86.118300, 86.615600, 87.112900,
87.610200, 88.107500, 88.604800, 89.102100, 89.599400, 90.096700,
90.594000, 91.091300, 91.588600, 92.085900, 92.583200, 93.080500,
93.577800, 94.075100, 94.572400, 95.069700, 95.567000, 96.064300,
96.561600, 97.058900, 97.556200, 98.053500, 98.550800, 99.048100,
99.545400, 100.042700, 100.540000, 101.037300, 101.534600, 102.031900,
102.529200, 103.026500, 103.523800, 104.021100, 104.518400, 105.015700,
105.513000, 106.010300, 106.507600, 107.004900, 107.502200, 107.999500,
108.496800, 108.994100, 109.491400, 109.988700, 110.486000, 110.983300,
111.480600, 111.977900, 112.475200, 112.972500, 113.469800, 113.967100,
114.464400, 114.961700, 115.459000, 115.956300, 116.453600, 116.950900,
117.448200, 117.945500, 118.442800, 118.940100, 119.437400, 119.934700,
120.432000, 120.929300, 121.426600, 121.923900, 122.421200, 122.918500,
123.415800, 123.913100, 124.410400, 124.907700, 125.405000, 125.902300,
126.399600, 126.896900, 127.394200, 127.891500, 128.388800, 128.886100,
129.383400, 129.880700, 130.378000, 130.875300, 131.372600, 131.869900,
132.367200, 132.864500, 133.361800, 133.859100, 134.356400, 134.853700,
135.351000, 135.848300, 136.345600, 136.842900, 137.340200, 137.837500,
138.334800, 138.832100, 139.329400, 139.826700, 140.324000, 140.821300,
141.318600, 141.815900, 142.313200, 142.810500, 143.307800, 143.805100,
144.302400, 144.799700, 145.297000, 145.794300, 146.291600, 146.788900,
147.286200, 147.783500, 148.280800, 148.778100, 149.275400, 149.772700,
150.270000, 150.767300, 151.264600, 151.761900, 152.259200, 152.756500,
153.253800, 153.751100, 154.248400, 154.745700, 155.243000, 155.740300,
156.237600, 156.734900, 157.232200, 157.729500, 158.226800, 158.724100,
159.221400, 159.718700, 160.216000, 160.713300, 161.210600, 161.707900,
162.205200, 162.702500, 163.199800, 163.697100, 164.194400, 164.691700,
165.189000, 165.686300, 166.183600, 166.680900, 167.178200, 167.675500,
168.172800, 168.670100, 169.167400, 169.664700, 170.162000, 170.659300,
171.156600, 171.653900, 172.151200, 172.648500, 173.145800, 173.643100,
174.140400, 174.637700, 175.135000, 175.632300, 176.129600, 176.626900,
177.124200, 177.621500, 178.118800, 178.616100, 179.113400, 179.610700,
180.108000, 180.605300, 181.102600, 181.599900, 182.097200, 182.594500,
183.091800, 183.589100, 184.086400, 184.583700, 185.081000, 185.578300,
186.075600, 186.572900, 187.070200, 187.567500, 188.064800, 188.562100,
189.059400, 189.556700, 190.054000, 190.551300, 191.048600, 191.545900,
192.043200, 192.540500, 193.037800, 193.535100, 194.032400, 194.529700,
195.027000, 195.524300, 196.021600, 196.518900, 197.016200, 197.513500,
198.010800, 198.508100, 199.005400, 199.502700, 200.000000};
Double_t w[]={
1.000000, 0.692593, 0.673384, 0.666273, 0.657285, 0.637723,
0.625881, 0.643590, 0.606100, 0.589007, 0.567824, 0.578705,
0.538530, 0.517937, 0.528278, 0.515918, 0.539461, 0.466186,
0.489869, 0.468402, 0.465017, 0.453336, 0.460769, 0.474638,
0.456347, 0.434471, 0.427478, 0.435435, 0.410934, 0.366431,
0.382240, 0.379513, 0.394249, 0.386837, 0.353103, 0.382138,
0.377497, 0.389479, 0.378736, 0.347933, 0.354605, 0.352077,
0.324443, 0.358792, 0.339968, 0.359052, 0.330734, 0.318291,
0.333703, 0.358644, 0.335819, 0.332213, 0.309051, 0.309975,
0.331626, 0.304407, 0.309819, 0.312097, 0.312462, 0.320411,
0.280401, 0.302311, 0.315863, 0.281479, 0.310003, 0.296911,
0.313676, 0.281071, 0.294163, 0.306500, 0.283462, 0.274867,
0.307149, 0.270555, 0.282264, 0.287373, 0.307849, 0.278675,
0.286990, 0.278269, 0.300105, 0.286799, 0.265674, 0.275140,
0.285702, 0.257352, 0.267714, 0.248204, 0.252220, 0.255678,
0.282946, 0.268257, 0.282375, 0.262675, 0.275564, 0.248345,
0.236259, 0.291914, 0.259936, 0.241338, 0.267389, 0.285044,
0.289419, 0.253594, 0.284568, 0.231840, 0.260008, 0.268527,
0.275363, 0.224115, 0.281260, 0.257913, 0.295152, 0.264399,
0.232287, 0.282533, 0.223431, 0.255756, 0.244471, 0.221695,
0.272450, 0.284244, 0.253682, 0.270717, 0.275403, 0.240323,
0.245081, 0.241859, 0.216340, 0.244789, 0.220291, 0.238478,
0.224691, 0.244058, 0.266117, 0.271478, 0.242012, 0.267321,
0.248494, 0.253343, 0.255606, 0.235458, 0.241079, 0.233223,
0.226813, 0.259224, 0.234239, 0.258606, 0.210892, 0.238186,
0.243271, 0.222678, 0.213437, 0.273939, 0.247966, 0.232548,
0.263438, 0.222089, 0.272111, 0.248818, 0.244295, 0.238368,
0.236908, 0.248776, 0.232604, 0.231194, 0.227117, 0.231152,
0.282140, 0.229778, 0.232631, 0.261794, 0.216633, 0.253471,
0.242157, 0.227406, 0.269335, 0.230547, 0.210618, 0.285872,
0.248776, 0.229875, 0.242728, 0.227388, 0.220567, 0.222062,
0.235950, 0.224087, 0.228895, 0.208287, 0.235999, 0.208696,
0.230367, 0.267667, 0.220484, 0.233402, 0.233815, 0.250455,
0.253120, 0.219556, 0.230980, 0.236661, 0.222395, 0.226111,
0.198315, 0.210555, 0.202813, 0.208594, 0.235976, 0.221490,
0.243059, 0.204901, 0.216987, 0.229039, 0.231466, 0.221975,
0.231220, 0.253638, 0.250448, 0.260291, 0.328345, 0.205739,
0.222014, 0.251513, 0.279427, 0.270506, 0.248409, 0.222472,
0.291632, 0.227796, 0.248769, 0.276896, 0.214742, 0.200139,
0.230693, 0.226031, 0.268900, 0.185160, 0.245353, 0.205843,
0.231155, 0.219122, 0.214811, 0.199763, 0.274179, 0.217598,
0.274988, 0.237244, 0.211820, 0.225459, 0.252799, 0.235948,
0.224986, 0.245385, 0.237770, 0.213373, 0.229737, 0.215487,
0.234453, 0.249684, 0.239435, 0.250422, 0.257194, 0.214762,
0.212266, 0.228988, 0.253798, 0.201607, 0.239946, 0.205245,
0.231670, 0.212774, 0.206768, 0.231563, 0.189388, 0.227926,
0.227203, 0.237754, 0.221628, 0.211138, 0.203322, 0.200985,
0.231780, 0.220294, 0.232686, 0.234243, 0.218264, 0.255870,
0.223213, 0.238670, 0.213713, 0.213064, 0.246700, 0.233446,
0.221503, 0.206767, 0.200722, 0.226179, 0.237425, 0.239229,
0.238611, 0.240419, 0.247806, 0.215923, 0.205298, 0.232778,
0.272312, 0.226773, 0.258103, 0.223287, 0.269404, 0.203398,
0.223782, 0.204213, 0.229664, 0.234040, 0.228419, 0.203936,
0.263686, 0.199141, 0.236127, 0.214058, 0.204611, 0.224324,
0.292140, 0.190735, 0.235157, 0.213018, 0.257085, 0.190554,
0.203197, 0.213044, 0.237023, 0.214243, 0.193562, 0.262403,
0.206256, 0.221396, 0.233588, 0.256611, 0.249731, 0.226683,
0.199330, 0.251026, 0.222596, 0.201941, 0.186374, 0.221038,
0.196555, 0.222560, 0.299419, 0.231979, 0.242924, 0.198310,
0.217628, 0.235458, 0.278595, 0.218624, 0.277305, 0.239109,
0.205600, 0.253715, 0.221173, 0.218195, 0.277647, 0.241974,
0.268748, 0.268128, 0.292258, 0.249420, 0.191034, 0.219506,
0.216502, 0.250677, 0.193386, 0.201310, 0.259464, 0.255351,
0.269628, 0.221063, 0.294079, 0.196726, 0.233634, 0.221870,
0.216236, 0.197259, 0.247433, 0.272765, 0.294079, 0.236336,
0.206396, 0.238524, 0.247846, 0.269519, 0.237141, 0.230611,
0.201712, 0.242225, 0.255565, 0.258738};
wDD = 0.512813;
wND = 0.518820;
wSD = -1;

    Mmin = bin[0];
    Mmax = bin[nbin];
    if(M<Mmin || M>Mmax) return kTRUE;

    Int_t ibin=nbin-1;
    for(Int_t i=1; i<=nbin; i++)
      if(M<=bin[i]) {
	ibin=i-1;
	//	printf("Mi> %f && Mi< %f\n", bin[i-1], bin[i]);
	break;
      }
    wSD=w[ibin];
    return kTRUE;
  }


  else if(TMath::Abs(fEnergyCMS-7000)<1 ){
const Int_t nbin=400;
Double_t bin[]={
1.080000, 1.577300, 2.074600, 2.571900, 3.069200, 3.566500,
4.063800, 4.561100, 5.058400, 5.555700, 6.053000, 6.550300,
7.047600, 7.544900, 8.042200, 8.539500, 9.036800, 9.534100,
10.031400, 10.528700, 11.026000, 11.523300, 12.020600, 12.517900,
13.015200, 13.512500, 14.009800, 14.507100, 15.004400, 15.501700,
15.999000, 16.496300, 16.993600, 17.490900, 17.988200, 18.485500,
18.982800, 19.480100, 19.977400, 20.474700, 20.972000, 21.469300,
21.966600, 22.463900, 22.961200, 23.458500, 23.955800, 24.453100,
24.950400, 25.447700, 25.945000, 26.442300, 26.939600, 27.436900,
27.934200, 28.431500, 28.928800, 29.426100, 29.923400, 30.420700,
30.918000, 31.415300, 31.912600, 32.409900, 32.907200, 33.404500,
33.901800, 34.399100, 34.896400, 35.393700, 35.891000, 36.388300,
36.885600, 37.382900, 37.880200, 38.377500, 38.874800, 39.372100,
39.869400, 40.366700, 40.864000, 41.361300, 41.858600, 42.355900,
42.853200, 43.350500, 43.847800, 44.345100, 44.842400, 45.339700,
45.837000, 46.334300, 46.831600, 47.328900, 47.826200, 48.323500,
48.820800, 49.318100, 49.815400, 50.312700, 50.810000, 51.307300,
51.804600, 52.301900, 52.799200, 53.296500, 53.793800, 54.291100,
54.788400, 55.285700, 55.783000, 56.280300, 56.777600, 57.274900,
57.772200, 58.269500, 58.766800, 59.264100, 59.761400, 60.258700,
60.756000, 61.253300, 61.750600, 62.247900, 62.745200, 63.242500,
63.739800, 64.237100, 64.734400, 65.231700, 65.729000, 66.226300,
66.723600, 67.220900, 67.718200, 68.215500, 68.712800, 69.210100,
69.707400, 70.204700, 70.702000, 71.199300, 71.696600, 72.193900,
72.691200, 73.188500, 73.685800, 74.183100, 74.680400, 75.177700,
75.675000, 76.172300, 76.669600, 77.166900, 77.664200, 78.161500,
78.658800, 79.156100, 79.653400, 80.150700, 80.648000, 81.145300,
81.642600, 82.139900, 82.637200, 83.134500, 83.631800, 84.129100,
84.626400, 85.123700, 85.621000, 86.118300, 86.615600, 87.112900,
87.610200, 88.107500, 88.604800, 89.102100, 89.599400, 90.096700,
90.594000, 91.091300, 91.588600, 92.085900, 92.583200, 93.080500,
93.577800, 94.075100, 94.572400, 95.069700, 95.567000, 96.064300,
96.561600, 97.058900, 97.556200, 98.053500, 98.550800, 99.048100,
99.545400, 100.042700, 100.540000, 101.037300, 101.534600, 102.031900,
102.529200, 103.026500, 103.523800, 104.021100, 104.518400, 105.015700,
105.513000, 106.010300, 106.507600, 107.004900, 107.502200, 107.999500,
108.496800, 108.994100, 109.491400, 109.988700, 110.486000, 110.983300,
111.480600, 111.977900, 112.475200, 112.972500, 113.469800, 113.967100,
114.464400, 114.961700, 115.459000, 115.956300, 116.453600, 116.950900,
117.448200, 117.945500, 118.442800, 118.940100, 119.437400, 119.934700,
120.432000, 120.929300, 121.426600, 121.923900, 122.421200, 122.918500,
123.415800, 123.913100, 124.410400, 124.907700, 125.405000, 125.902300,
126.399600, 126.896900, 127.394200, 127.891500, 128.388800, 128.886100,
129.383400, 129.880700, 130.378000, 130.875300, 131.372600, 131.869900,
132.367200, 132.864500, 133.361800, 133.859100, 134.356400, 134.853700,
135.351000, 135.848300, 136.345600, 136.842900, 137.340200, 137.837500,
138.334800, 138.832100, 139.329400, 139.826700, 140.324000, 140.821300,
141.318600, 141.815900, 142.313200, 142.810500, 143.307800, 143.805100,
144.302400, 144.799700, 145.297000, 145.794300, 146.291600, 146.788900,
147.286200, 147.783500, 148.280800, 148.778100, 149.275400, 149.772700,
150.270000, 150.767300, 151.264600, 151.761900, 152.259200, 152.756500,
153.253800, 153.751100, 154.248400, 154.745700, 155.243000, 155.740300,
156.237600, 156.734900, 157.232200, 157.729500, 158.226800, 158.724100,
159.221400, 159.718700, 160.216000, 160.713300, 161.210600, 161.707900,
162.205200, 162.702500, 163.199800, 163.697100, 164.194400, 164.691700,
165.189000, 165.686300, 166.183600, 166.680900, 167.178200, 167.675500,
168.172800, 168.670100, 169.167400, 169.664700, 170.162000, 170.659300,
171.156600, 171.653900, 172.151200, 172.648500, 173.145800, 173.643100,
174.140400, 174.637700, 175.135000, 175.632300, 176.129600, 176.626900,
177.124200, 177.621500, 178.118800, 178.616100, 179.113400, 179.610700,
180.108000, 180.605300, 181.102600, 181.599900, 182.097200, 182.594500,
183.091800, 183.589100, 184.086400, 184.583700, 185.081000, 185.578300,
186.075600, 186.572900, 187.070200, 187.567500, 188.064800, 188.562100,
189.059400, 189.556700, 190.054000, 190.551300, 191.048600, 191.545900,
192.043200, 192.540500, 193.037800, 193.535100, 194.032400, 194.529700,
195.027000, 195.524300, 196.021600, 196.518900, 197.016200, 197.513500,
198.010800, 198.508100, 199.005400, 199.502700, 200.000000};
Double_t w[]={
1.000000, 0.798640, 0.770197, 0.717990, 0.699434, 0.692093,
0.620940, 0.655294, 0.636496, 0.633483, 0.600518, 0.594548,
0.550759, 0.550413, 0.528926, 0.525877, 0.506701, 0.504117,
0.486973, 0.492930, 0.461833, 0.488695, 0.438874, 0.460274,
0.451551, 0.454424, 0.449270, 0.387571, 0.427556, 0.428740,
0.401682, 0.402260, 0.408961, 0.395417, 0.387426, 0.378602,
0.357826, 0.359125, 0.348494, 0.363710, 0.356117, 0.340363,
0.337637, 0.396084, 0.341777, 0.340551, 0.348838, 0.344014,
0.340468, 0.317654, 0.355584, 0.326023, 0.373416, 0.312298,
0.326522, 0.290540, 0.335557, 0.318689, 0.327544, 0.319501,
0.331754, 0.312728, 0.282263, 0.274937, 0.303867, 0.307820,
0.289344, 0.268934, 0.288908, 0.290018, 0.291369, 0.295242,
0.289067, 0.277685, 0.267957, 0.267559, 0.320229, 0.265060,
0.305931, 0.305352, 0.262064, 0.281879, 0.287780, 0.270033,
0.270814, 0.276667, 0.271531, 0.275283, 0.258189, 0.287969,
0.251247, 0.301527, 0.267230, 0.245860, 0.293125, 0.253421,
0.272396, 0.243637, 0.236206, 0.278452, 0.246544, 0.263165,
0.230484, 0.231102, 0.258281, 0.244707, 0.270111, 0.248295,
0.246942, 0.245592, 0.272766, 0.254243, 0.199647, 0.262590,
0.226710, 0.243836, 0.214153, 0.233535, 0.235207, 0.234245,
0.247698, 0.248379, 0.241463, 0.265230, 0.223242, 0.236191,
0.252700, 0.231865, 0.228352, 0.229101, 0.237385, 0.246485,
0.254706, 0.245565, 0.224841, 0.257301, 0.240968, 0.202379,
0.236943, 0.241683, 0.220809, 0.219014, 0.213015, 0.244470,
0.221042, 0.198996, 0.225295, 0.264466, 0.200600, 0.228143,
0.250138, 0.216784, 0.268002, 0.275734, 0.218144, 0.229866,
0.235443, 0.208909, 0.215067, 0.189743, 0.216741, 0.242686,
0.200234, 0.218882, 0.245991, 0.222815, 0.219576, 0.209773,
0.205247, 0.203855, 0.202534, 0.192536, 0.223447, 0.225810,
0.220583, 0.248421, 0.223424, 0.206033, 0.203130, 0.221225,
0.223763, 0.239216, 0.252311, 0.206893, 0.228461, 0.233591,
0.201512, 0.242382, 0.215147, 0.232578, 0.207116, 0.239213,
0.196215, 0.184367, 0.235135, 0.189768, 0.274084, 0.206267,
0.203283, 0.248828, 0.250214, 0.217465, 0.232080, 0.215150,
0.207936, 0.176789, 0.191338, 0.188655, 0.196181, 0.223691,
0.254257, 0.216874, 0.213536, 0.221399, 0.192024, 0.212534,
0.218169, 0.226635, 0.201191, 0.212700, 0.211634, 0.232353,
0.223636, 0.208605, 0.249132, 0.183681, 0.221842, 0.187136,
0.203772, 0.249575, 0.217713, 0.205193, 0.193941, 0.219707,
0.264226, 0.182105, 0.207183, 0.220845, 0.236571, 0.182390,
0.212914, 0.186266, 0.195361, 0.217665, 0.204527, 0.188804,
0.222832, 0.191193, 0.201073, 0.185616, 0.195011, 0.200828,
0.181434, 0.233707, 0.202925, 0.211992, 0.173100, 0.205258,
0.182695, 0.216520, 0.218422, 0.209490, 0.211257, 0.215801,
0.220652, 0.211409, 0.195731, 0.194957, 0.198888, 0.234926,
0.221377, 0.229822, 0.176700, 0.172322, 0.212265, 0.206133,
0.170355, 0.253305, 0.198688, 0.240043, 0.225384, 0.174729,
0.195820, 0.200093, 0.196912, 0.212308, 0.200490, 0.240415,
0.159744, 0.173686, 0.223853, 0.213604, 0.193779, 0.179328,
0.180873, 0.237481, 0.179640, 0.235857, 0.202847, 0.167690,
0.245093, 0.215504, 0.205848, 0.184408, 0.201626, 0.209651,
0.236306, 0.217803, 0.188534, 0.187861, 0.161663, 0.221718,
0.152044, 0.190412, 0.173505, 0.208995, 0.174575, 0.180734,
0.247704, 0.203388, 0.194184, 0.169679, 0.182703, 0.213402,
0.191808, 0.178604, 0.178116, 0.198452, 0.174687, 0.169809,
0.222851, 0.156811, 0.229170, 0.210359, 0.178557, 0.248570,
0.208536, 0.192571, 0.178912, 0.224505, 0.170822, 0.205492,
0.330973, 0.160924, 0.203724, 0.213255, 0.205827, 0.187162,
0.181252, 0.191723, 0.238106, 0.182398, 0.202358, 0.212066,
0.201255, 0.168159, 0.185219, 0.229176, 0.158222, 0.164092,
0.215405, 0.200724, 0.234811, 0.184488, 0.213112, 0.198577,
0.219622, 0.160671, 0.179349, 0.206681, 0.206091, 0.165618,
0.165203, 0.174442, 0.179287, 0.187318, 0.163481, 0.217440,
0.160381, 0.177025, 0.204385, 0.163676, 0.210733, 0.186519,
0.230701, 0.191764, 0.185119, 0.190770, 0.242987, 0.186812,
0.202906, 0.161935, 0.182426, 0.197922, 0.181475, 0.155903,
0.175006, 0.223482, 0.202706, 0.218108};
wDD = 0.207705;
wND = 0.289628;
wSD = -1;

    Mmin = bin[0];
    Mmax = bin[nbin];

    if(M<Mmin || M>Mmax) return kTRUE;

    Int_t ibin=nbin-1;
    for(Int_t i=1; i<=nbin; i++)
      if(M<=bin[i]) {
	ibin=i-1;
	//	printf("Mi> %f && Mi< %f\n", bin[i-1], bin[i]);
	break;
      }
    wSD=w[ibin];
    return kTRUE;
  }

  return kFALSE;
}

Bool_t AliGenPythia::IsFromHeavyFlavor(Int_t ipart)
{
// Check if this is a heavy flavor decay product
  TParticle *  part = (TParticle *) fParticles.At(ipart);
  Int_t mpdg = TMath::Abs(part->GetPdgCode());
  Int_t mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));
  //
  // Light hadron
  if (mfl >= 4 && mfl < 6) return kTRUE;

  Int_t imo = part->GetFirstMother()-1;
  TParticle* pm = part;
  //
  // Heavy flavor hadron produced by generator
  while (imo >  -1) {
    pm  =  (TParticle*)fParticles.At(imo);
    mpdg = TMath::Abs(pm->GetPdgCode());
    mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));
    if ((mfl > 3) && (mfl <6) && mpdg > 400) return kTRUE;
    imo = pm->GetFirstMother()-1;
  }
  return kFALSE;
}

Bool_t AliGenPythia::CheckDetectorAcceptance(Float_t phi, Float_t eta, Int_t iparticle)
{
  // check the eta/phi correspond to the detectors acceptance
  // iparticle is the index of the particle to be checked, for PHOS rotation case
  if     (fCheckPHOS   && IsInPHOS  (phi,eta,iparticle)) return kTRUE;
  else if(fCheckEMCAL  && IsInEMCAL (phi,eta)) return kTRUE;
  else if(fCheckBarrel && IsInBarrel(    eta)) return kTRUE;
  else if(fCheckBarrelCalos && IsInBarrelCalorimeters(phi,eta)) return kTRUE;
  else                                         return kFALSE;
}

Bool_t AliGenPythia::TriggerOnSelectedParticles(Int_t np)
{
  // Select events with fragmentation photon, decay photon, pi0, eta or other hadrons going to PHOS or EMCAL or central barrel,
  // implemented primaryly for kPyJets, but extended to any kind of process.

  //printf("Check: frag photon %d, pi0 %d, eta %d, electron %d, hadron %d, decay %d, PHOS %d, EMCAL %d, Barrel %d Barrel calos %d\n",
  //       fFragPhotonInCalo,fPi0InCalo, fEtaInCalo,fEleInCalo,fHadronInCalo,fDecayPhotonInCalo,fCheckPHOS,fCheckEMCAL, fCheckBarrel, fCheckBarrelCalos);

  Bool_t ok = kFALSE;
  for (Int_t i=0; i< np; i++) {

    TParticle* iparticle = (TParticle *) fParticles.At(i);

    Int_t pdg          = iparticle->GetPdgCode();
    Int_t status       = iparticle->GetStatusCode();
    Int_t imother      = iparticle->GetFirstMother() - 1;

    TParticle* pmother = 0x0;
    Int_t momStatus    = -1;
    Int_t momPdg       = -1;
    if(imother > 0 ){
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
      if(TMath::Abs(pdg) > 23 && pdg !=221 && pdg != 111) // avoid photons, electrons, muons, neutrinos and eta or pi0
        // (in case neutral mesons were declared stable)
        ok = kTRUE;
    }
    //Electron
    else if (fEleInCalo && status == 1 && TMath::Abs(pdg) == 11)
    {
        ok = kTRUE;
    }
    //Fragmentation photon
    else if (fFragPhotonInCalo && pdg == 22 && status == 1)
    {
      if(momStatus != 11) ok = kTRUE ;  // No photon from hadron decay
    }
    // Decay photon
    else if (fDecayPhotonInCalo && !fForceNeutralMeson2PhotonDecay && pdg == 22) // pi0 status can be 1 or 11 depending on decay settings, work only for 11
    {
      if( momStatus == 11)
      {
        //if(iparticle->Pt() > fTriggerParticleMinPt) printf("Decay photon! pdg %d, status %d, pt %2.2f, mom: pdg %d, pt %2.2f\n",
        //                                                   pdg,status,iparticle->Pt(),momPdg,pmother->Pt());
        ok = kTRUE ;  // photon from hadron decay

        //In case only decays from pi0 or eta requested
        if(fPi0InCalo && momPdg!=111) ok = kFALSE;
        if(fEtaInCalo && momPdg!=221) ok = kFALSE;
      }

    }
    // Pi0 or Eta particle
    else if ((fPi0InCalo || fEtaInCalo))
    {
      if(fDecayPhotonInCalo && !fForceNeutralMeson2PhotonDecay ) continue ;

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
      //Just check if the selected particle falls in the acceptance
      if(!fForceNeutralMeson2PhotonDecay )
      {
        //printf("\t Check acceptance! \n");
        Float_t phi = iparticle->Phi()*180./TMath::Pi(); //Convert to degrees
        Float_t eta =TMath::Abs(iparticle->Eta()); //in calos etamin=-etamax

        if(CheckDetectorAcceptance(phi,eta,i))
        {
          ok =kTRUE;
          AliDebug(1,Form("Selected trigger pdg %d, status %d, pt %2.2f, eta %2.2f, phi %2.2f\n",pdg,status,iparticle->Pt(), eta, phi));
          //printf("\t Accept \n");
          break;
        }
        else ok = kFALSE;
      }
      //Mesons have several decay modes, select only those decaying into 2 photons
      else if(fForceNeutralMeson2PhotonDecay && (fPi0InCalo || fEtaInCalo))
      {
        // In case we want the pi0/eta trigger,
        // check the decay mode (2 photons)

        //printf("\t Force decay 2 gamma\n");

        Int_t ndaughters = iparticle->GetNDaughters();
        if(ndaughters != 2){
          ok=kFALSE;
          continue;
        }

        TParticle*d1 = (TParticle *) fParticles.At(iparticle->GetDaughter(0)-1);
        TParticle*d2 = (TParticle *) fParticles.At(iparticle->GetDaughter(1)-1);
        if(!d1 || !d2) {
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

        if(pdgD1 != 22  || pdgD2 != 22){
          ok = kFALSE;
          continue;
        }

        //printf("\t accept decay\n");

        //Trigger on the meson, not on the daughter
        if(!fDecayPhotonInCalo){

          Float_t phi = iparticle->Phi()*180./TMath::Pi(); //Convert to degrees
          Float_t eta =TMath::Abs(iparticle->Eta()); //in calos etamin=-etamax

          if(CheckDetectorAcceptance(phi,eta,i))
          {
            //printf("\t Accept meson pdg %d\n",pdg);
            ok =kTRUE;
            AliDebug(1,Form("Selected trigger pdg %d (decay), status %d, pt %2.2f, eta %2.2f, phi %2.2f\n",pdg,status,iparticle->Pt(), eta, phi));
            break;
          } else {
            ok=kFALSE;
            continue;
          }
        }

        //printf("Check daughters acceptance\n");

        //Trigger on the meson daughters
        //Photon 1
        Float_t phi = d1->Phi()*180./TMath::Pi(); //Convert to degrees
        Float_t eta =TMath::Abs(d1->Eta()); //in calos etamin=-etamax
        if(d1->Pt() > fTriggerParticleMinPt)
        {
          //printf("\t Check acceptance photon 1! \n");
          if(CheckDetectorAcceptance(phi,eta,i))
          {
            //printf("\t Accept Photon 1\n");
            ok =kTRUE;
            AliDebug(1,Form("Selected trigger pdg %d (decay), status %d, pt %2.2f, eta %2.2f, phi %2.2f\n",pdg,status,iparticle->Pt(), eta, phi));
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
          //printf("\t Check acceptance photon 2! \n");
          if(CheckDetectorAcceptance(phi,eta,i))
          {
            //printf("\t Accept Photon 2\n");
            ok =kTRUE;
            AliDebug(1,Form("Selected trigger pdg %d (decay), status %d, pt %2.2f, eta %2.2f, phi %2.2f\n",pdg,status,iparticle->Pt(), eta, phi));
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

void AliGenPythia::SetSeed(UInt_t seed)
{
  GetRandom()->SetSeed(seed);
}

void AliGenPythia::WriteXsection(const Char_t *fname) {
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
