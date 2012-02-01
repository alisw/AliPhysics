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
//
//The design of this class is based on the AliFMDParameters class. Its purpose
//is to hold parameters for the analysis such as background correction and 
//fit functions.
//
//Author: Hans Hjersing Dalsgaard, NBI, hans.dalsgaard@cern.ch
//

// #include "AliFMDDebug.h"		   // ALILOG_H
#include "AliFMDAnaParameters.h"	   // ALIFMDPARAMETERS_H
//#include <AliCDBManager.h>         // ALICDBMANAGER_H
//#include <AliCDBEntry.h>           // ALICDBMANAGER_H
//#include "AliFMDRing.h"
#include <AliLog.h>
#include <Riostream.h>
#include <sstream>
#include <TSystem.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliMultiplicity.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliFMDAnaCalibBackgroundCorrection.h"
#include "AliFMDAnaCalibEnergyDistribution.h"
#include "AliFMDAnaCalibEventSelectionEfficiency.h"
#include "AliFMDAnaCalibSharingEfficiency.h"
#include "TFile.h"

//====================================================================
ClassImp(AliFMDAnaParameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//const char* AliFMDAnaParameters::fgkBackgroundCorrection  = "FMD/Correction/Background";
//const char* AliFMDAnaParameters::fgkEnergyDists = "FMD/Correction/EnergyDistribution";
const char* AliFMDAnaParameters::fgkBackgroundID         = "background";
const char* AliFMDAnaParameters::fgkEnergyDistributionID = "energydistributions";
const char* AliFMDAnaParameters::fgkEventSelectionEffID  = "eventselectionefficiency";
const char* AliFMDAnaParameters::fgkSharingEffID         = "sharingefficiency";
//____________________________________________________________________
AliFMDAnaParameters* AliFMDAnaParameters::fgInstance = 0;

//____________________________________________________________________

AliFMDAnaParameters* 
AliFMDAnaParameters::Instance() 
{
  // Get static instance 
  if (!fgInstance) fgInstance = new AliFMDAnaParameters;
  return fgInstance;
}

//____________________________________________________________________
AliFMDAnaParameters::AliFMDAnaParameters() :
  fIsInit(kFALSE),
  fBackground(0),
  fEnergyDistribution(0),
  fEventSelectionEfficiency(0),
  fSharingEfficiency(0),
  fCorner1(4.2231, 26.6638),
  fCorner2(1.8357, 27.9500),
  fEnergyPath("$ALICE_ROOT/PWG2/FORWARD/corrections/EnergyDistribution"),
  fBackgroundPath("$ALICE_ROOT/PWG2/FORWARD/corrections/Background"),
  fEventSelectionEffPath("$ALICE_ROOT/PWG2/FORWARD/corrections/EventSelectionEfficiency"),
  fSharingEffPath("$ALICE_ROOT/PWG2/FORWARD/corrections/SharingEfficiency"),
  fProcessPrimary(kFALSE),
  fProcessHits(kFALSE),
  fTrigger(kMB1),
  fEnergy(k900),
  fMagField(k5G),
  fSpecies(kPP),
  fPhysicsSelection(0),
  fRealData(kFALSE),
  fSPDlowLimit(0),
  fSPDhighLimit(999999999),
  fCentralSelection(kFALSE),
  fSharingObjectPresent(kTRUE),
  fNumberOfEtaBinsToCut(1),
  fEtaLowBinLimits(),
  fEtaHighBinLimits(),
  fTriggerInel(kFALSE),
  fTriggerNSD(kFALSE),
  fTriggerEmpty(kFALSE),
  fUseBuiltInNSD(kFALSE),
  fInelGtZero(kFALSE),
  fRunDndeta(kTRUE),
  fRunBFCorrelation(kTRUE),
  fRunMultiplicity(kTRUE)
{
  // Default constructor 
  //  fPhysicsSelection = new AliPhysicsSelection;
  //  fPhysicsSelection->SetAnalyzeMC(kTRUE); //For the background correction. This is reset in Init if relevant
  // fPhysicsSelection->SetUseBXNumbers(kFALSE);
  
  
  //  AliBackgroundSelection* backgroundSelection = new AliBackgroundSelection("bg","bg");
  //  backgroundSelection->Init();
  //  fPhysicsSelection->AddBackgroundIdentification(backgroundSelection);
  //fPhysicsSelection->Initialize(104792);
  // Do not use this - it is only for IO 
  fgInstance = this;
  
}
//____________________________________________________________________
const char* AliFMDAnaParameters::GetPath(const char* species) const
{
  //Get path of object
  static TString*  path = new TString();
 
  if(species == fgkBackgroundID)
    path->Form("%s/%s_%d_%d_%d_%d_%d_%d.root",
	      fBackgroundPath.Data(),
	      fgkBackgroundID,
	      fEnergy,
	      fTrigger,
	      fMagField,
	      fSpecies,
	      fInelGtZero,
	      0);
  if(species == fgkEnergyDistributionID)
    path->Form("%s/%s_%d_%d_%d_%d_%d_%d.root",
	      fEnergyPath.Data(),
	      fgkEnergyDistributionID,
	      fEnergy,
	      fTrigger,
	      fMagField,
	      fSpecies,
	      fRealData,
	      0);
  if(species == fgkEventSelectionEffID)
    path->Form("%s/%s_%d_%d_%d_%d_%d_%d.root",
	      fEventSelectionEffPath.Data(),
	      fgkEventSelectionEffID,
	      fEnergy,
	      fTrigger,
	      fMagField,
	      fSpecies,
	      fInelGtZero,
	      0);
  if(species == fgkSharingEffID)
    path->Form("%s/%s_%d_%d_%d_%d_%d_%d.root",
	      fSharingEffPath.Data(),
	      fgkSharingEffID,
	      fEnergy,
	      fTrigger,
	      fMagField,
	      fSpecies,
	      0,
	      0);
  return path->Data();
}
//____________________________________________________________________
void AliFMDAnaParameters::Init(Bool_t forceReInit, UInt_t what)
{
  // Initialize the parameters manager.  We need to get stuff from files here.
  
  /*  AliPhysicsSelection* test = (AliPhysicsSelection*)((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetEventSelection();
  
    if(fPhysicsSelection) {
    if(!fRealData) {
      fPhysicsSelection->SetAnalyzeMC(kTRUE);
    }
    else fPhysicsSelection->SetAnalyzeMC(kFALSE);
  }
  */
  if (forceReInit) fIsInit = kFALSE;
  if (fIsInit) return;
  if (what & kBackgroundCorrection)       InitBackground();
  if (what & kEnergyDistributions)        InitEnergyDists();
  if (what & kEventSelectionEfficiency)   InitEventSelectionEff();
  if (what & kSharingEfficiency)          InitSharingEff();
  
  fIsInit = kTRUE;
  
  if(fBackground)
    FindEtaLimits();
  
 
    
}
//____________________________________________________________________

void AliFMDAnaParameters::InitBackground() {
  //Init background correction objects.
  
  TFile* fin = TFile::Open(GetPath(fgkBackgroundID));
  
  if (!fin) return;
  
  fBackground = dynamic_cast<AliFMDAnaCalibBackgroundCorrection*>(fin->Get(fgkBackgroundID));
  if (!fBackground) AliFatal("Invalid background object from CDB");
  
  

}

//____________________________________________________________________

void AliFMDAnaParameters::InitEnergyDists() {
  //Init energy distributions
    
  TFile* fin = TFile::Open(GetPath(fgkEnergyDistributionID));
  
  if (!fin) return;
  
  fEnergyDistribution = dynamic_cast<AliFMDAnaCalibEnergyDistribution*>(fin->Get(fgkEnergyDistributionID));
  
  if (!fEnergyDistribution) AliFatal("Invalid background object from CDB");
  
}

//____________________________________________________________________

void AliFMDAnaParameters::InitEventSelectionEff() {
  //Init event selection objects

  TFile* fin = TFile::Open(GetPath(fgkEventSelectionEffID));
			    
  if (!fin) return;
  
  fEventSelectionEfficiency = dynamic_cast<AliFMDAnaCalibEventSelectionEfficiency*>(fin->Get(fgkEventSelectionEffID));
  if (!fEventSelectionEfficiency) AliFatal("Invalid background object from CDB");
  
}

//____________________________________________________________________

void AliFMDAnaParameters::InitSharingEff() {
  //Initialize the sharing efficiency
  fSharingObjectPresent = kTRUE;
  TFile* fin = TFile::Open(GetPath(fgkSharingEffID));
			    
  if (!fin) {
    fSharingObjectPresent = kFALSE;
    return; 
  }
  
  fSharingEfficiency = dynamic_cast<AliFMDAnaCalibSharingEfficiency*>(fin->Get(fgkSharingEffID));
  if (!fSharingEfficiency) {
    fSharingObjectPresent = kFALSE;
    return; 
  }
  
}
//____________________________________________________________________
void AliFMDAnaParameters::FindEtaLimits() {
  //Find eta limits for analysis
  fEtaLowBinLimits.SetBins(4,0,4,2,0,2,GetNvtxBins(),0,GetNvtxBins());
  fEtaHighBinLimits.SetBins(4,0,4,2,0,2,GetNvtxBins(),0,GetNvtxBins());
  for(Int_t det=0; det<=3;det++) {
    Int_t nRings = (det<=1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t ringChar = (ir == 0 ? 'I' : 'O');
      for(Int_t v =0; v<GetNvtxBins(); v++) {
	fEtaLowBinLimits.SetBinContent(det,ir,v,GetFirstEtaBinFromMap(v, det, ringChar));
	fEtaHighBinLimits.SetBinContent(det,ir,v,GetLastEtaBinFromMap(v, det, ringChar));
	//std::cout<<det<<"   "<<ringChar<<"   "<<fEtaLowBinLimits.GetBinContent(det,ir,v)<<"   "<<fEtaHighBinLimits.GetBinContent(det,ir,v)<<std::endl;
      }
    }
  }

}
//____________________________________________________________________
void AliFMDAnaParameters::SetEnergy(Float_t cmsNNGeV) 
{
  //Set energy
  if (TMath::Abs(cmsNNGeV - 900.)   < 10)  fEnergy = k900;
  if (TMath::Abs(cmsNNGeV - 2400.)  < 10)  fEnergy = k2400;
  if (TMath::Abs(cmsNNGeV - 2750.)  < 10)  fEnergy = k2750;
  if (TMath::Abs(cmsNNGeV - 5500.)  < 40)  fEnergy = k5500;
  if (TMath::Abs(cmsNNGeV - 7000.)  < 10)  fEnergy = k7000;
  if (TMath::Abs(cmsNNGeV - 10000.) < 10)  fEnergy = k10000;
  if (TMath::Abs(cmsNNGeV - 14000.) < 10)  fEnergy = k14000;
}
//____________________________________________________________________
void AliFMDAnaParameters::SetMagField(Float_t bkG) 
{
  //Set magnetic field
  if (TMath::Abs(bkG - 5.) < 1 ) fMagField = k5G;
  if (TMath::Abs(bkG + 5.) < 1 ) fMagField = k5Gnegative;
  if (TMath::Abs(bkG) < 1)       fMagField = k0G;
}
//____________________________________________________________________
void AliFMDAnaParameters::SetCollisionSystem(const TString& sys) 
{
  //Set the collision system
  TString s(sys);
  s.ToLower();
  if      (s.Contains("p-p")   || s.Contains("pp"))   fSpecies = kPP; 
  else if (s.Contains("pb-pb") || s.Contains("pbpb")) fSpecies = kPbPb;
  else if (s.Contains("a-a")   || s.Contains("aa"))   fSpecies = kPbPb;
}
  
//____________________________________________________________________
void AliFMDAnaParameters::SetParametersFromESD(AliESDEvent* esd) 
{
  //Set the parameters from the ESD header information

  SetCollisionSystem(esd->GetBeamType());

  Float_t energy   = esd->GetBeamEnergy();
  // Correct to center of mass per nucleon (cmsNN) - LHC gives it as 
  // cmsNN * Z / 2 
  if (fSpecies == kPbPb) energy = energy / 208 * 82;
  SetEnergy(2*energy); 
  SetMagField(esd->GetMagneticField());

  // Float_t magfield = esd->GetCurrentL3(); 
  // if (TMath::Abs(magfield - 30000.) < 10 ) fMagField = k5G;
  // if (TMath::Abs(magfield + 30000.) < 10 ) fMagField = k5Gnegative;
  // if (TMath::Abs(magfield) < 10 )          fMagField = k0G;
  
  
  
  Init(kTRUE);
  
}

//____________________________________________________________________

void AliFMDAnaParameters::PrintStatus(Bool_t showpaths) const
{
  //Print current status
  TString energystring;
  switch(fEnergy) {
  case k900:
    energystring.Form("900 GeV");   break;
  case k2400:
    energystring.Form("2400 GeV"); break;
  case k2750:
    energystring.Form("2750 GeV"); break;
  case k5500:
    energystring.Form("5500 GeV"); break;
  case k7000:
    energystring.Form("7000 GeV");  break;
  case k10000:
    energystring.Form("10000 GeV"); break;
  case k14000:
    energystring.Form("14000 GeV"); break;
  default:
    energystring.Form("invalid energy"); break;
  }
  TString triggerstring;
  switch(fTrigger) {
  case kMB1:
    triggerstring.Form("Minimum bias 1");   break;
  case kMB2:
    triggerstring.Form("Minimum bias 2");   break;
  case kSPDFASTOR:
    triggerstring.Form("SPD FAST OR");   break;
  case kNOCTP:
    triggerstring.Form("NO TRIGGER TEST");   break;
  default:
    energystring.Form("invalid trigger"); break;
  }
  TString magstring;
  switch(fMagField) {
  case k5G:
    magstring.Form("+5 kGaus");   break;
  case k0G:
    magstring.Form("0 kGaus");   break;
  case k5Gnegative:
    magstring.Form("-5 kGaus");   break;
  default:
    magstring.Form("invalid mag field %d", fMagField); break;
  }
  TString collsystemstring;
  switch(fSpecies) {
  case kPP:
    collsystemstring.Form("p-p");   break;
  case kPbPb:
    collsystemstring.Form("Pb-Pb");   break;
  default:
    collsystemstring.Form("invalid collision system");   break;
  }
  
  TString datastring;
  
  if(fRealData) 
    datastring.Form("Nature");
  else 
    datastring.Form("MC"); 
  
  TString inelString;
  if(fInelGtZero) inelString = "INEL > 0";
  else inelString = "INEL";
  
  std::cout<<"Energy       = "<<energystring.Data()<<std::endl;
  std::cout<<"Trigger      = "<<triggerstring.Data()<<std::endl;
  std::cout<<"Mag Field    = "<<magstring.Data()<<std::endl;
  std::cout<<"Coll System  = "<<collsystemstring.Data()<<std::endl;
  std::cout<<"Data origin  = "<<datastring.Data()<<std::endl;
  std::cout<<"Basic trigger: "<<inelString.Data()<<std::endl;

  if (showpaths) {
    TString bg = GetPath(fgkBackgroundID);
    TString es = GetPath(fgkEventSelectionEffID);
    TString ed = GetPath(fgkEnergyDistributionID);
    TString me = GetPath(fgkSharingEffID);
    std::cout << "2nd maps:     " << bg << "\n"
	      << "Event sel.:   " << es << "\n"
	      << "Energy dist.: " << ed << "\n"
	      << "Merge eff.:   " << me << std::endl;
  }
  
}

//____________________________________________________________________
Float_t AliFMDAnaParameters::GetVtxCutZ() {
  //Get the z vtx cut in analysis
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return -1;
  }
  
  return fBackground->GetVtxCutZ();
}

//____________________________________________________________________
Int_t AliFMDAnaParameters::GetNvtxBins() {
  //Get number of vtx bins
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return -1;
  }
  
  return fBackground->GetNvtxBins();
}
//____________________________________________________________________
TH1F* AliFMDAnaParameters::GetEnergyDistribution(Int_t det, Char_t ring, Float_t eta) {
  
  return fEnergyDistribution->GetEnergyDistribution(det, ring, eta);
}
//____________________________________________________________________
TH1F* AliFMDAnaParameters::GetEmptyEnergyDistribution(Int_t det, Char_t ring) {
  
  return fEnergyDistribution->GetEmptyEnergyDistribution(det, ring);
}
//____________________________________________________________________
TH1F* AliFMDAnaParameters::GetRingEnergyDistribution(Int_t det, Char_t ring) {
  
  return fEnergyDistribution->GetRingEnergyDistribution(det, ring);
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetSigma(Int_t det, Char_t ring, Float_t eta) {
  //Get sigma of Landau fits
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist       = GetEnergyDistribution(det,ring, eta);
  TF1*  fitFunc           = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) {
    //AliWarning(Form("No function for FMD%d%c, eta %f",det,ring,eta));
    return 1024;
  }
  Float_t sigma           = fitFunc->GetParameter(2);
  return sigma;
}


//____________________________________________________________________
Float_t AliFMDAnaParameters::GetMPV(Int_t det, Char_t ring, Float_t eta) {
  //Get MPV of landau fits
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  //AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  
  if(!fitFunc) {
    AliWarning(Form("No function for FMD%d%c, eta %f (%d)",det,ring,eta,
		    GetEtaBin(eta)));
    return 1024;
  }
  
  Float_t mpv           = fitFunc->GetParameter(1);
  return mpv;
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetConstant(Int_t det, Char_t ring, Float_t eta) {
  //Get constant parameter of Landau fits
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) {
    AliWarning(Form("No function for FMD%d%c, eta %f",det,ring,eta));
    return 0;
  }
    
  Float_t mpv           = fitFunc->GetParameter(0);
  return mpv;
}
//____________________________________________________________________
Float_t 
AliFMDAnaParameters::Get2MIPWeight(Int_t det, Char_t ring, Float_t eta) 
{
  //Get 2 MIP weights of convoluted Landau fits
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) return 0;
  Float_t twoMIPweight    = fitFunc->GetParameter(3);
  
  
  
  if(twoMIPweight < 1e-05)
    twoMIPweight = 0;
  
  return twoMIPweight;
}
//____________________________________________________________________
Float_t 
AliFMDAnaParameters::Get3MIPWeight(Int_t det, Char_t ring, Float_t eta) 
{
  //Get 3 MIP weights of convoluted Landau fits
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) return 0;
  Float_t threeMIPweight    = fitFunc->GetParameter(4);
  
  if(threeMIPweight < 1e-05)
    threeMIPweight = 0;
  
  Float_t twoMIPweight    = fitFunc->GetParameter(3);
  
  if(twoMIPweight < 1e-05)
    threeMIPweight = 0;
    
  return threeMIPweight;
}
//____________________________________________________________________
Int_t AliFMDAnaParameters::GetNetaBins() 
{
  return GetBackgroundCorrection(1,'I',5)->GetNbinsX();  
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetEtaMin() 
{
  return GetBackgroundCorrection(1,'I',5)->GetXaxis()->GetXmin();
} 
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetEtaMax() 
{
  return GetBackgroundCorrection(1,'I',5)->GetXaxis()->GetXmax();
}
//____________________________________________________________________
Int_t AliFMDAnaParameters::GetEtaBin(Float_t eta) 
{  
  TAxis testaxis(GetNetaBins(),GetEtaMin(),GetEtaMax());
  Int_t binnumber = testaxis.FindBin(eta) ;
  
  return binnumber;
}
//____________________________________________________________________

TH2F* AliFMDAnaParameters::GetBackgroundCorrection(Int_t det, 
						   Char_t ring, 
						   Int_t vtxbin) {
  //Get background correction histogram
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  
  
  if(vtxbin > fBackground->GetNvtxBins()) {
    AliWarning(Form("No background object for vertex bin %d", vtxbin));
    return 0;
  } 
  
  return fBackground->GetBgCorrection(det,ring,vtxbin);
}
//____________________________________________________________________
TH2F* AliFMDAnaParameters::GetBackgroundCorrectionNSD(Int_t det, 
						      Char_t ring, 
						      Int_t vtxbin) {
  //Get background correction histogram for NSD event class
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  if(vtxbin > fBackground->GetNvtxBins()) {
    AliWarning(Form("No background object for vertex bin %d", vtxbin));
    return 0;
  } 
  
  if(fBackground->GetNSDBgCorrection(det,ring,vtxbin))
    return fBackground->GetNSDBgCorrection(det,ring,vtxbin);
  else 
    AliWarning("No NSD background map. You get usual one. "
	       "Difference is probably negligible");
  
  return fBackground->GetBgCorrection(det,ring,vtxbin); 
}
//____________________________________________________________________

TH1F* AliFMDAnaParameters::GetDoubleHitCorrection(Int_t det, 
						  Char_t ring) {
  //Get correction for several hits in strips for p+p
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  return fBackground->GetDoubleHitCorrection(det,ring);
}
//_____________________________________________________________________
TH1F* AliFMDAnaParameters::GetSPDDeadCorrection(Int_t vtxbin) {
						
  //Get correction for several hits in strips for p+p
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  return fBackground->GetSPDDeadCorrection(vtxbin);
}
//_____________________________________________________________________
TH1F* AliFMDAnaParameters::GetFMDDeadCorrection(Int_t vtxbin) {
						
  //Get correction for several hits in strips for p+p
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  return fBackground->GetFMDDeadCorrection(vtxbin);
}
//_____________________________________________________________________
Float_t AliFMDAnaParameters::GetEventSelectionEfficiency(Int_t vtxbin) {
  //Get event selection efficiency object
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  return fEventSelectionEfficiency->GetCorrection(vtxbin);

}
//_____________________________________________________________________
Float_t AliFMDAnaParameters::GetVtxSelectionEffFromMC() {
  //Get the vtx selection from MC calculation
   if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
   }
   return fEventSelectionEfficiency->GetVtxToTriggerRatio();


}
//_____________________________________________________________________
TH2F* AliFMDAnaParameters::GetEventSelectionEfficiency(TString trig, Int_t vtxbin, Char_t ring) {
  //Get event selection efficiency object
  
  //TString test = trig;
  if(!trig.Contains("NSD") && !trig.Contains("INEL")) {
    AliWarning("Event selection efficiency only available for INEL and NSD");
    return 0;
  }
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  return fEventSelectionEfficiency->GetCorrection(trig,vtxbin,ring);

}
//_____________________________________________________________________
TH1F* AliFMDAnaParameters::GetSharingEfficiency(Int_t det, Char_t ring, Int_t vtxbin) {
  //Get sharing efficiency object
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  return fSharingEfficiency->GetSharingEff(det,ring,vtxbin);

}
//_____________________________________________________________________
TH1F* AliFMDAnaParameters::GetSharingEfficiencyTrVtx(Int_t det, Char_t ring, Int_t vtxbin) {
  //Get sharing efficiency object TrVtx
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  return fSharingEfficiency->GetSharingEffTrVtx(det,ring,vtxbin);

}
//_____________________________________________________________________
Float_t AliFMDAnaParameters::GetMaxR(Char_t ring) const {
  //Get max R of ring
  Float_t radius = 0;
  if(ring == 'I')
    radius = 17.2;
  else if(ring == 'O')
    radius = 28.0;
  else
    AliWarning("Unknown ring - must be I or O!");
  
  return radius;
}
//_____________________________________________________________________
Float_t AliFMDAnaParameters::GetMinR(Char_t ring) const{
  //Get min R of ring
  Float_t radius = 0;
  if(ring == 'I')
    radius = 4.5213;
  else if(ring == 'O')
    radius = 15.4;
  else
    AliWarning("Unknown ring - must be I or O!");
  
  return radius;

}
//_____________________________________________________________________
void AliFMDAnaParameters::SetCorners(Char_t ring) {
  //Set corners (taken from nominal geometry)
  if(ring == 'I') {
    fCorner1.Set(4.9895, 15.3560);
    fCorner2.Set(1.8007, 17.2000);
  }
  else {
    fCorner1.Set(4.2231, 26.6638);
    fCorner2.Set(1.8357, 27.9500);
  }
  
}
//_____________________________________________________________________
Float_t AliFMDAnaParameters::GetPhiFromSector(UShort_t det, Char_t ring, UShort_t sec) const
{
  //Get phi from sector
  Int_t nsec = (ring == 'I' ? 20 : 40);
  Float_t basephi = 0;
  if(det == 1) 
    basephi = 1.72787594; 
  if(det == 2 && ring == 'I')
    basephi = 0.15707963;
  if(det == 2 && ring == 'O')
    basephi = 0.078539818;
  if(det == 3 && ring == 'I')
    basephi = 2.984513044;
  if(det == 3 && ring == 'O')
    basephi = 3.06305289;
  
  Float_t step = 2*TMath::Pi() / nsec;
  Float_t phi = 0;
  if(det == 3)
    phi = basephi - sec*step;
  else
    phi = basephi + sec*step;
  
  if(phi < 0) 
    phi = phi +2*TMath::Pi();
  if(phi > 2*TMath::Pi() )
    phi = phi - 2*TMath::Pi();
  
  return phi;
}
//_____________________________________________________________________
Float_t AliFMDAnaParameters::GetEtaFromStrip(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip, Float_t zvtx) const
{
  //Calculate eta from strip with vertex (redundant with AliESDFMD::Eta)
  Float_t   rad       = GetMaxR(ring)-GetMinR(ring);
  Float_t   nStrips   = (ring == 'I' ? 512 : 256);
  Float_t   segment   = rad / nStrips;
  Float_t   r         = GetMinR(ring) + segment*strip;
  Float_t   z         = 0;
  Int_t hybrid = sec / 2;
  
  if(det == 1) {
    if(!(hybrid%2)) z = 320.266; else z = 319.766;
  }
  if(det == 2 && ring == 'I' ) {
    if(!(hybrid%2)) z = 83.666; else z = 83.166;
  }
  if(det == 2 && ring == 'O' ) {
    if(!(hybrid%2)) z = 74.966; else z = 75.466;
  }
  if(det == 3 && ring == 'I' ) {
    if(!(hybrid%2)) z = -63.066; else z = -62.566;
  }
  if(det == 3 && ring == 'O' ) {
    if(!(hybrid%2)) z = -74.966; else z = -75.466;
  }
  
  //std::cout<<det<<"   "<<ring<<"   "<<sec<<"   "<<hybrid<<"    "<<z<<std::endl;
  
  // Float_t   r     = TMath::Sqrt(TMath::Power(x,2)+TMath::Power(y,2));
  Float_t   theta = TMath::ATan2(r,z-zvtx);
  Float_t   eta   = -1*TMath::Log(TMath::Tan(0.5*theta));
  
  return eta;
}

//_____________________________________________________________________

Bool_t AliFMDAnaParameters::GetVertex(const AliESDEvent* esd, Double_t* vertexXYZ) 
{
  //Get the vertex from the ESD
  const AliESDVertex* vertex = esd->GetPrimaryVertexSPD();
  
  if (!vertex) return kFALSE;

  vertex->GetXYZ(vertexXYZ);

  //if(vertexXYZ[0] == 0 || vertexXYZ[1] == 0 )
  //  return kFALSE;
 
  if(vertex->GetNContributors() <= 0)
    return kFALSE;
  
  if(vertex->GetZRes() > 0.1 ) 
    return kFALSE;
  
  return vertex->GetStatus();
  
}
//____________________________________________________________________
void AliFMDAnaParameters::SetTriggerStatus(const AliESDEvent *esd) {

  //ULong64_t triggerMask = esd->GetTriggerMask();
  
  AliPhysicsSelection* centralPhysicsSelection = (AliPhysicsSelection*)((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetEventSelection();
  
  if(!centralPhysicsSelection && !fPhysicsSelection) {
    
    std::cout<<"Creating AliPhysicsSelection object due to absence of central object"<<std::endl;
    fPhysicsSelection = new AliPhysicsSelection;
    fPhysicsSelection->SetAnalyzeMC(!fRealData);
    // fPhysicsSelection->SetUseBXNumbers(kFALSE);
    
    
    AliBackgroundSelection* backgroundSelection = new AliBackgroundSelection("bg","bg");
    backgroundSelection->Init();
    fPhysicsSelection->AddBackgroundIdentification(backgroundSelection);
    
  }
  
  TString triggers = esd->GetFiredTriggerClasses();
  
  AliTriggerAnalysis tAna;
  
  fTriggerInel  = kFALSE;
  fTriggerNSD   = kFALSE;
  fTriggerEmpty = kFALSE;
  
  UInt_t inel = kFALSE;
  
  if(centralPhysicsSelection) {
    inel= ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  }
  else
    inel = fPhysicsSelection->IsCollisionCandidate(esd);
  
  
  
  if(fInelGtZero && inel) {
    const AliMultiplicity* spdmult = esd->GetMultiplicity();
    Int_t nCentralTracklets = 0;
    Int_t j = 0;
    while( nCentralTracklets < 1 && j< spdmult->GetNumberOfTracklets() ) {
      if(TMath::Abs(spdmult->GetEta(j)) < 1) nCentralTracklets++;
      j++;
    }
    if(nCentralTracklets < 1) inel = kFALSE;      
  }
  
  //std::cout<<fTriggerInel<<std::endl;
  if(inel) {
    fTriggerInel = kTRUE;
  }
  
  Bool_t nsd = kFALSE;
  if(fUseBuiltInNSD) { 
    if ((tAna.FMDTrigger(esd, AliTriggerAnalysis::kASide) || tAna.V0Trigger(esd, AliTriggerAnalysis::kASide, kFALSE) == AliTriggerAnalysis::kV0BB) && (tAna.FMDTrigger(esd, AliTriggerAnalysis::kCSide) || tAna.V0Trigger(esd, AliTriggerAnalysis::kCSide, kFALSE) == AliTriggerAnalysis::kV0BB))
      nsd = kTRUE;
  }
  else nsd = tAna.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kNSD1);
  
  if(fTriggerInel && nsd) {
    fTriggerNSD = kTRUE;
  }
  if(triggers.Contains("CBEAMB-ABCE-NOPF-ALL")) {
    fTriggerEmpty = kTRUE;
  }
  
  
  /*switch (fTrigger) {
  case kMB1: {
    if( fPhysicsSelection->IsCollisionCandidate(esd)) {
      fTriggerInel = kTRUE;
    }
    break;
    
  }
  case kMB2: { 
    if (triggerMask & spdFO && ((triggerMask & v0left) || (triggerMask & v0right)))
      return kTRUE;
    break;
  }
  case kSPDFASTOR: {
    if (triggerMask & spdFO)
      return kTRUE;
    break;
  }
  case kNOCTP: {
    return kTRUE;
    break;
  }
  case kEMPTY: {
  if(triggers.Contains("CBEAMB-ABCE-NOPF-ALL")) 
  return kTRUE;
    break;
  }
  case kNSD: {
    if(fPhysicsSelection->IsCollisionCandidate(esd) && tAna.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kNSD1))
      return kTRUE;
    break;
  }
    
  }//switch
  */
  
}
/*
//____________________________________________________________________
Bool_t AliFMDAnaParameters::IsEventTriggered(const AliESDEvent *esd, Trigger trig) {
  //Did we have trig trigger ?
  Trigger old = fTrigger;
  fTrigger = trig;
  Bool_t retval = IsEventTriggered(esd);
  fTrigger = old;
  return retval;

}
*/
//____________________________________________________________________
Bool_t AliFMDAnaParameters::IsEventTriggered(Trigger trigger) {
  // check if the event was triggered
  
  if (fCentralSelection) return kTRUE;
  switch (trigger) {
  
  case kMB1:
    return fTriggerInel;
    break;
  case kNSD:
    return fTriggerNSD;
    break;
  case kEMPTY:
    return fTriggerEmpty;
    break;
  case kNOCTP:
    return kTRUE;
    break;
  default:
    AliWarning("Trigger not implemented!!!");
    break;
    
    
  }
  return kFALSE;
  
}
    
//____________________________________________________________________
Float_t 
AliFMDAnaParameters::GetStripLength(Char_t ring, UShort_t strip)  
{
  //Get length of a strip
  
  Float_t rad        = GetMaxR(ring)-GetMinR(ring);
  Float_t   nStrips   = (ring == 'I' ? 512 : 256);
  Float_t segment    = rad / nStrips;
  
  //TVector2* corner1  = fmdring.GetVertex(2);  
  // TVector2* corner2  = fmdring.GetVertex(3);
  
  SetCorners(ring);
  /*
  std::cout<<GetMaxR(ring)<<"   "<<fmdring.GetMaxR()<<std::endl;
  std::cout<<GetMinR(ring)<<"   "<<fmdring.GetMinR()<<std::endl;
  std::cout<<corner1->X()<<"   "<<fCorner1.X()<<std::endl;
  std::cout<<corner2->X()<<"   "<<fCorner2.X()<<std::endl;
  std::cout<<corner1->Y()<<"   "<<fCorner1.Y()<<std::endl;
  std::cout<<corner2->Y()<<"   "<<fCorner2.Y()<<std::endl;*/
  Float_t slope      = (fCorner1.Y() - fCorner2.Y()) / (fCorner1.X() - fCorner2.X());
  Float_t constant   = (fCorner2.Y()*fCorner1.X()-(fCorner2.X()*fCorner1.Y())) / (fCorner1.X() - fCorner2.X());
  Float_t radius     = GetMinR(ring) + strip*segment;
  
  Float_t d          = TMath::Power(TMath::Abs(radius*slope),2) + TMath::Power(radius,2) - TMath::Power(constant,2);
  
  Float_t arclength  = GetBaseStripLength(ring,strip);
  if(d>0) {
    
    Float_t x        = (-1*TMath::Sqrt(d) -slope*constant) / (1+TMath::Power(slope,2));
    Float_t y        = slope*x + constant;
    Float_t theta    = TMath::ATan2(x,y);
    
    if(x < fCorner1.X() && y > fCorner1.Y()) {
      arclength = radius*theta;                        //One sector since theta is by definition half-hybrid
      
    }
    
  }
  
  return arclength;
  
  
}
//____________________________________________________________________
Float_t 
AliFMDAnaParameters::GetBaseStripLength(Char_t ring, UShort_t strip) const  
{  
  //Get length of strip assuming that corners are not cut away
  Float_t rad             = GetMaxR(ring)-GetMinR(ring);
  Float_t nStrips         = (ring == 'I' ? 512 : 256);
  Float_t nSec            = (ring == 'I' ? 20 : 40);
  Float_t segment         = rad / nStrips;
  Float_t basearc         = 2*TMath::Pi() / (0.5*nSec); // One hybrid: 36 degrees inner, 18 outer
  Float_t radius          = GetMinR(ring) + strip*segment;
  Float_t basearclength   = 0.5*basearc * radius;                // One sector   
  
  return basearclength;
}
//____________________________________________________________________
Int_t    AliFMDAnaParameters::GetFirstEtaBinFromMap(Int_t vtxbin, Int_t det, Char_t ring)
{
  //Get the first eta bin from the bg map
  TH2F* hBg = GetBackgroundCorrection(det,ring,vtxbin);
  
  if(det == 0) return hBg->GetXaxis()->FindBin(-1.95);
  
  Int_t firstbin = -1;
  Int_t nNonZeroFirst = 0;
  
  for(Int_t i=1;i<=hBg->GetNbinsX();i++) {
    if(nNonZeroFirst == fNumberOfEtaBinsToCut && firstbin==-1) firstbin = i;
    
    for(Int_t j=1;j<=hBg->GetNbinsY();j++) {
      
      Float_t value = hBg->GetBinContent(i,j);
      
      if(value > 0.001 && nNonZeroFirst<fNumberOfEtaBinsToCut)
	{nNonZeroFirst++; break;}
      
      
    }
  }
  
  return firstbin;

}
//____________________________________________________________________
Int_t    AliFMDAnaParameters::GetLastEtaBinFromMap(Int_t vtxbin, Int_t det, Char_t ring)
{
  //Get the last eta bin from the map
  TH2F* hBg = GetBackgroundCorrection(det,ring,vtxbin);
  Int_t lastbin=-1;
  Int_t nNonZeroLast = 0;
  
  if(det == 0) return hBg->GetXaxis()->FindBin(1.95);
  
  for(Int_t i=hBg->GetNbinsX();i>0;i--) {
    if(nNonZeroLast == fNumberOfEtaBinsToCut && lastbin==-1) lastbin = i;
    
    for(Int_t j=1;j<=hBg->GetNbinsY();j++) {
      
      Float_t value = hBg->GetBinContent(i,j);
      
      if(value > 0.001 && nNonZeroLast<fNumberOfEtaBinsToCut)
	{nNonZeroLast++; break; }
      
      
    }
  }
  
  return lastbin;
}

//____________________________________________________________________
Int_t    AliFMDAnaParameters::GetFirstEtaBinToInclude(Int_t vtxbin, Int_t det, Char_t ring)
{
  Int_t ringNumber = (ring == 'I' ? 0 : 1);
  return (Int_t)fEtaLowBinLimits.GetBinContent(det,ringNumber,vtxbin);

}

//____________________________________________________________________
Int_t    AliFMDAnaParameters::GetLastEtaBinToInclude(Int_t vtxbin, Int_t det, Char_t ring)
{
  Int_t ringNumber = (ring == 'I' ? 0 : 1);
  return (Int_t)fEtaHighBinLimits.GetBinContent(det,ringNumber,vtxbin);
  
}
//____________________________________________________________________
//
// EOF
//
