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

#include "AliFMDDebug.h"		   // ALILOG_H
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

//====================================================================
ClassImp(AliFMDAnaParameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

const char* AliFMDAnaParameters::fgkBackgroundCorrection  = "FMD/Correction/Background";
const char* AliFMDAnaParameters::fgkEnergyDists = "FMD/Correction/EnergyDistribution";
const char* AliFMDAnaParameters::fkBackgroundID = "background";
const char* AliFMDAnaParameters::fkEnergyDistributionID = "energydistributions";
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
  fCorner1(4.2231, 26.6638),
  fCorner2(1.8357, 27.9500),
  fEnergyPath("$ALICE_ROOT/FMD/Correction/EnergyDistribution/energydistributions.root"),
  fBackgroundPath("$ALICE_ROOT/FMD/Correction/Background/background.root")
{
  
  
  //fVerticies.Add(new TVector2(4.2231, 26.6638));
  // fVerticies.Add(new TVector2(1.8357, 27.9500));
  // Default constructor 
}
//____________________________________________________________________
void AliFMDAnaParameters::Init(Bool_t forceReInit, UInt_t what)
{
  // Initialize the parameters manager.  We need to get stuff from the
  // CDB here. 
  if (forceReInit) fIsInit = kFALSE;
  if (fIsInit) return;
  if (what & kBackgroundCorrection)       InitBackground();
  if (what & kEnergyDistributions)        InitEnergyDists();
  
  
  fIsInit = kTRUE;
}
//____________________________________________________________________

void AliFMDAnaParameters::InitBackground() {
  
  //AliCDBEntry*   background = GetEntry(fgkBackgroundCorrection);
  TFile* fin = TFile::Open(fBackgroundPath.Data());
  
  if (!fin) return;
  
  fBackground = dynamic_cast<AliFMDAnaCalibBackgroundCorrection*>(fin->Get(fkBackgroundID));
  if (!fBackground) AliFatal("Invalid background object from CDB");
  
}
//____________________________________________________________________

void AliFMDAnaParameters::InitEnergyDists() {
  
  TFile* fin = TFile::Open(fEnergyPath.Data());
  //AliCDBEntry*   edist = GetEntry(fgkEnergyDists);
  if (!fin) return;
  
  fEnergyDistribution = dynamic_cast<AliFMDAnaCalibEnergyDistribution*>(fin->Get(fkEnergyDistributionID));
  
  if (!fEnergyDistribution) AliFatal("Invalid background object from CDB");
  
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetVtxCutZ() {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return -1;
  }
  
  return fBackground->GetVtxCutZ();
}

//____________________________________________________________________
Int_t AliFMDAnaParameters::GetNvtxBins() {
  
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
Float_t AliFMDAnaParameters::GetSigma(Int_t det, Char_t ring, Float_t eta) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist       = GetEnergyDistribution(det,ring, eta);
  TF1*  fitFunc           = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) {
    AliWarning(Form("No function for FMD%d%c, eta %f",det,ring,eta));
    return 1024;
  }
  Float_t sigma           = fitFunc->GetParameter(2);
  return sigma;
}


//____________________________________________________________________
Float_t AliFMDAnaParameters::GetMPV(Int_t det, Char_t ring, Float_t eta) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) {
    AliWarning(Form("No function for FMD%d%c, eta %f",det,ring,eta));
    return 1024;
  }
    
  Float_t MPV           = fitFunc->GetParameter(1);
  return MPV;
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::Get2MIPWeight(Int_t det, Char_t ring, Float_t eta) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) return 0;
  Float_t TwoMIPweight    = fitFunc->GetParameter(3);
  
  
  
  if(TwoMIPweight < 1e-05)
    TwoMIPweight = 0;
  
  return TwoMIPweight;
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::Get3MIPWeight(Int_t det, Char_t ring, Float_t eta) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) return 0;
  Float_t ThreeMIPweight    = fitFunc->GetParameter(4);
  
  if(ThreeMIPweight < 1e-05)
    ThreeMIPweight = 0;
  
  Float_t TwoMIPweight    = fitFunc->GetParameter(3);
  
  if(TwoMIPweight < 1e-05)
    ThreeMIPweight = 0;
    
  return ThreeMIPweight;
}
//____________________________________________________________________
Int_t AliFMDAnaParameters::GetNetaBins() {
  return GetBackgroundCorrection(1,'I',0)->GetNbinsX();
  
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetEtaMin() {

  return GetBackgroundCorrection(1,'I',0)->GetXaxis()->GetXmin();
} 
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetEtaMax() {

return GetBackgroundCorrection(1,'I',0)->GetXaxis()->GetXmax();

}
//____________________________________________________________________

TH2F* AliFMDAnaParameters::GetBackgroundCorrection(Int_t det, 
						   Char_t ring, 
						   Int_t vtxbin) {
  
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
//_____________________________________________________________________
Float_t AliFMDAnaParameters::GetMaxR(Char_t ring) {
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
Float_t AliFMDAnaParameters::GetMinR(Char_t ring) {
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
Float_t AliFMDAnaParameters::GetPhiFromSector(UShort_t det, Char_t ring, UShort_t sec)
{
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
Float_t AliFMDAnaParameters::GetEtaFromStrip(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip, Float_t zvtx)
{
  // AliFMDRing fmdring(ring);
  // fmdring.Init();
  Float_t   rad       = GetMaxR(ring)-GetMinR(ring);
  Float_t   Nstrips   = (ring == 'I' ? 512 : 256);
  Float_t   segment   = rad / Nstrips;
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
}/*
//____________________________________________________________________
AliCDBEntry* AliFMDAnaParameters::GetEntry(const char* path, Bool_t fatal) const
{
  // Get an entry from the CDB or via preprocessor 
  AliCDBEntry* entry = 0;
  AliCDBManager* cdb = AliCDBManager::Instance();
  entry              = cdb->Get(path);
  
  if (!entry) { 
    TString msg(Form("No %s found in CDB, perhaps you need to "
		     "use AliFMDCalibFaker?", path));
    if (fatal) { AliFatal(msg.Data()); }
    else       AliLog::Message(AliLog::kWarning, msg.Data(), "FMD", 
			       "AliFMDParameters", "GetEntry", __FILE__, 
			       __LINE__);
    return 0;
  }
  return entry;
}
 */
//____________________________________________________________________
Float_t 
AliFMDAnaParameters::GetStripLength(Char_t ring, UShort_t strip)  
{
  //AliFMDRing fmdring(ring);
  // fmdring.Init();
  
  Float_t rad        = GetMaxR(ring)-GetMinR(ring);
  Float_t   Nstrips   = (ring == 'I' ? 512 : 256);
  Float_t segment    = rad / Nstrips;
  
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
AliFMDAnaParameters::GetBaseStripLength(Char_t ring, UShort_t strip)  
{  
  // AliFMDRing fmdring(ring);
  // fmdring.Init();
  Float_t rad             = GetMaxR(ring)-GetMinR(ring);
  Float_t Nstrips         = (ring == 'I' ? 512 : 256);
  Float_t Nsec            = (ring == 'I' ? 20 : 40);
  Float_t segment         = rad / Nstrips;
  Float_t basearc         = 2*TMath::Pi() / (0.5*Nsec); // One hybrid: 36 degrees inner, 18 outer
  Float_t radius          = GetMinR(ring) + strip*segment;
  Float_t basearclength   = 0.5*basearc * radius;                // One sector   
  
  return basearclength;
}
//____________________________________________________________________
//
// EOF
//
