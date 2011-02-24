#ifndef ALIDIELECTRONVARMANAGER_H
#define ALIDIELECTRONVARMANAGER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//#############################################################
//#                                                           # 
//#         Class AliDielectronVarManager                     #
//#         Class for management of available variables       #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Markus    Köhler,   GSI / M.Koehler@gsi.de              #
//#   Frederick Kramer,   Uni Ffm / Frederick.Kramer@cern.ch  #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################


#include <TNamed.h>

#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliVVertex.h>
#include <AliESDVertex.h>

#include <AliVParticle.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliAODPid.h>
#include <AliKFParticle.h>
#include <AliKFVertex.h>
#include <AliMCParticle.h>
#include <AliAODMCParticle.h>
#include <AliVTrack.h>  // ?

#include <AliExternalTrackParam.h>
#include <AliESDpid.h>
#include <AliCentrality.h>
#include <AliAODpidUtil.h>
#include <AliPID.h>

#include "AliDielectronPair.h"
#include "AliDielectronMC.h"
#include "AliDielectronPID.h"

class AliVEvent;

//________________________________________________________________
class AliDielectronVarManager : public TNamed {
  
public:

  // Particle specific variables
  enum ValueTypes {
    kPx = 0,                 // px
    kPy,                     // py
    kPz,                     // pz
    kPt,                     // transverse momentum
    kP,                      // momentum
    kXv,                     // vertex position in x
    kYv,                     // vertex position in y
    kZv,                     // vertex position in z
    kOneOverPt,              // 1/pt
    kPhi,                    // phi angle
    kTheta,                  // theta angle
    kEta,                    // pseudo-rapidity
    kY,                      // rapidity
    kE,                      // energy
    kM,                      // mass
    kCharge,                 // charge
    kNclsITS,                // number of clusters assigned in the ITS
    kNclsTPC,                // number of clusters assigned in the TPC
    kNclsTPCiter1,           // number of clusters assigned in the TPC after first iteration
    kNFclsTPC,               // number of findable clusters in the TPC
    kNFclsTPCr,              // number of findable clusters in the TPC with more robust definition
    kNFclsTPCrFrac,          // number of found/findable clusters in the TPC with more robust definition
    kTPCsignalN,             // number of points used for dEdx
    kTPCchi2Cl,              // chi2/cl in TPC
    kTrackStatus,            // track status bits
    
    kNclsTRD,                // number of clusters assigned in the TRD
    kTRDntracklets,          // number of TRD tracklets used for tracking/PID TODO: correct getter
    kTRDpidQuality,          // number of TRD tracklets used for PID
    kTRDprobEle,             // TRD electron pid probability
    kTRDprobPio,             // TRD electron pid probability
      
    kImpactParXY,            // Impact parameter in XY plane
    kImpactParZ,             // Impact parameter in Z
    kTrackLength,            // Track length
    kPdgCode,                // PDG code

    kPdgCodeMother,          // PDG code of the mother

    kNumberOfDaughters,      // number of daughters
    kHaveSameMother,         // check that particles have the same mother (MC)
    kIsJpsiPrimary,          // check if the particle is primary (MC)
    kITSsignal,		     // ITS dE/dx signal
    kITSsignalSSD1,	     // SSD1 dE/dx signal
    kITSsignalSSD2,	     // SSD2 dE/dx signal
    kITSsignalSDD1,	     // SDD1 dE/dx signal
    kITSsignalSDD2,	     // SDD2 dE/dx signal
    kITSclusterMap,      // ITS cluster map
    kITSnSigmaEle,           // number of sigmas to the dE/dx electron line in the ITS
    kITSnSigmaPio,           // number of sigmas to the dE/dx pion line in the ITS
    kITSnSigmaMuo,           // number of sigmas to the dE/dx muon line in the ITS
    kITSnSigmaKao,           // number of sigmas to the dE/dx kaon line in the ITS
    kITSnSigmaPro,           // number of sigmas to the dE/dx proton line in the ITS

    kPIn,                    // momentum at inner wall of TPC (if available), used for PID
    kTPCsignal,              // TPC dE/dx signal
      
    kTPCnSigmaEle,           // number of sigmas to the dE/dx electron line in the TPC
    kTPCnSigmaPio,           // number of sigmas to the dE/dx pion line in the TPC
    kTPCnSigmaMuo,           // number of sigmas to the dE/dx muon line in the TPC
    kTPCnSigmaKao,           // number of sigmas to the dE/dx kaon line in the TPC
    kTPCnSigmaPro,           // number of sigmas to the dE/dx proton line in the TPC
      
    kTOFnSigmaEle,           // number of sigmas to the pion line in the TOF
    kTOFnSigmaPio,           // number of sigmas to the pion line in the TOF
    kTOFnSigmaMuo,           // number of sigmas to the muon line in the TOF
    kTOFnSigmaKao,           // number of sigmas to the kaon line in the TOF
    kTOFnSigmaPro,           // number of sigmas to the proton line in the TOF

    kKinkIndex0,             // kink index 0
      
    kParticleMax,             //
    // TODO: kRNClusters ??
  // AliDielectronPair specific variables
    kChi2NDF = kParticleMax, // Chi^2/NDF
    kDecayLength,            // decay length
    kR,                      // distance to the origin
    kOpeningAngle,           // opening angle
    // helicity picture: Z-axis is considered the direction of the mother's 3-momentum vector
    kThetaHE,                // theta in mother's rest frame in the helicity picture 
    kPhiHE,                  // phi in mother's rest frame in the helicity picture
    // Collins-Soper picture: Z-axis is considered the direction of the vectorial difference between 
    // the 3-mom vectors of target and projectile beams
    kThetaCS,                // theta in mother's rest frame in Collins-Soper picture
    kPhiCS,                  // phi in mother's rest frame in Collins-Soper picture
    kLegDist,                // distance of the legs
    kLegDistXY,              // distance of the legs in XY
    kDeltaEta,         // Absolute value of Delta Eta for the legs
    kDeltaPhi,           // Absolute value of Delta Phi for the legs
    kMerr,                   // error of mass calculation
    kDCA,                    // distance of closest approach TODO: not implemented yet
    kPairType,               // type of the pair, like like sign ++ unlikesign ...
    kPairMax,                 //
  // Event specific variables
    kXvPrim=kPairMax,        // prim vertex
    kYvPrim,                 // prim vertex
    kZvPrim,                 // prim vertex
    kXRes,                   // primary vertex x-resolution
    kYRes,                   // primary vertex y-resolution
    kZRes,                   // primary vertex z-resolution
    kNTrk,                   // number of tracks (or tracklets)
    kTracks,                 // ESD tracks
    kCentrality,             // event centrality fraction
    kNevents,                // event counter
    kNMaxValues              //
    // TODO: (for A+A) ZDCEnergy, impact parameter, Iflag??
  };
  

  AliDielectronVarManager();
  AliDielectronVarManager(const char* name, const char* title);
  virtual ~AliDielectronVarManager();
  static void Fill(const TObject* particle, Double_t * const values);

  static void InitESDpid(Int_t type=0);
  static void InitAODpidUtil(Int_t type=0);
  static void SetESDpid(AliESDpid * const pid) {fgESDpid=pid;}
  static AliESDpid* GetESDpid() {return fgESDpid;}
  static AliAODpidUtil* GetAODpidUtil() {return fgAODpidUtil;}
  static void SetEvent(AliVEvent * const ev);
  static Bool_t GetDCA(const AliAODTrack *track, Double_t d0z0[2]);
    
  static const char* GetValueName(Int_t i) { return (i>=0&&i<kNMaxValues)?fgkParticleNames[i]:""; }
private:

  static const char* fgkParticleNames[kNMaxValues];  //variable names

  static void FillVarVParticle(const AliVParticle *particle,         Double_t * const values);
  static void FillVarESDtrack(const AliESDtrack *particle,           Double_t * const values);
  static void FillVarAODTrack(const AliAODTrack *particle,           Double_t * const values);
  static void FillVarMCParticle(const AliMCParticle *particle,       Double_t * const values);
  static void FillVarAODMCParticle(const AliAODMCParticle *particle, Double_t * const values);
  static void FillVarDielectronPair(const AliDielectronPair *pair,   Double_t * const values);
  static void FillVarVEvent(const AliVEvent *event,                  Double_t * const values);
  static void FillVarESDEvent(const AliESDEvent *event,              Double_t * const values);
  static void FillVarAODEvent(const AliAODEvent *event,              Double_t * const values);
  static void FillVarMCEvent(const AliMCEvent *event,                Double_t * const values);
  
  static AliESDpid     *fgESDpid;             // ESD pid object
  static AliAODpidUtil *fgAODpidUtil;         // AOD pid object
  static AliVEvent     *fgEvent;              // current event pointer
  static AliKFVertex   *fgKFVertex;           // kf vertex
  
  AliDielectronVarManager(const AliDielectronVarManager &c);
  AliDielectronVarManager &operator=(const AliDielectronVarManager &c);
  
  ClassDef(AliDielectronVarManager,1);
};


//Inline functions
inline void AliDielectronVarManager::Fill(const TObject* object, Double_t * const values)
{
  //
  // Main function to fill all available variables according to the type of particle
  //

  if      (object->IsA() == AliESDtrack::Class())       FillVarESDtrack(static_cast<const AliESDtrack*>(object), values);
  else if (object->IsA() == AliAODTrack::Class())       FillVarAODTrack(static_cast<const AliAODTrack*>(object), values);
  else if (object->IsA() == AliMCParticle::Class())     FillVarMCParticle(static_cast<const AliMCParticle*>(object), values);
  else if (object->IsA() == AliAODMCParticle::Class())  FillVarAODMCParticle(static_cast<const AliAODMCParticle*>(object), values);
  else if (object->IsA() == AliDielectronPair::Class()) FillVarDielectronPair(static_cast<const AliDielectronPair*>(object), values);

  // Main function to fill all available variables according to the type of event
  
  else if (object->IsA() == AliVEvent::Class())         FillVarVEvent(static_cast<const AliVEvent*>(object), values);
  else if (object->IsA() == AliESDEvent::Class())       FillVarESDEvent(static_cast<const AliESDEvent*>(object), values);
  else if (object->IsA() == AliAODEvent::Class())       FillVarAODEvent(static_cast<const AliAODEvent*>(object), values);
  else if (object->IsA() == AliMCEvent::Class())        FillVarMCEvent(static_cast<const AliMCEvent*>(object), values);
//   else Error("Fill",Form("Type %s is not supported by AliDielectronVarManager!", object->ClassName())); //TODO: implement without object needed
}

inline void AliDielectronVarManager::FillVarVParticle(const AliVParticle *particle, Double_t * const values)
{
  //
  // Fill track information available in AliVParticle into an array
  //
  values[AliDielectronVarManager::kPx]        = particle->Px();
  values[AliDielectronVarManager::kPy]        = particle->Py();
  values[AliDielectronVarManager::kPz]        = particle->Pz();
  values[AliDielectronVarManager::kPt]        = particle->Pt();
  values[AliDielectronVarManager::kP]         = particle->P();

  values[AliDielectronVarManager::kXv]        = particle->Xv();
  values[AliDielectronVarManager::kYv]        = particle->Yv();
  values[AliDielectronVarManager::kZv]        = particle->Zv();

  values[AliDielectronVarManager::kOneOverPt] = particle->OneOverPt();
  values[AliDielectronVarManager::kPhi]       = particle->Phi();
  values[AliDielectronVarManager::kTheta]     = particle->Theta();
  values[AliDielectronVarManager::kEta]       = particle->Eta();
  values[AliDielectronVarManager::kY]         = particle->Y();
  
  values[AliDielectronVarManager::kE]         = particle->E();
  values[AliDielectronVarManager::kM]         = particle->M();
  values[AliDielectronVarManager::kCharge]    = particle->Charge();
  
  values[AliDielectronVarManager::kPdgCode]   = particle->PdgCode();
    
  if ( fgEvent ) AliDielectronVarManager::Fill(fgEvent, values);
}

inline void AliDielectronVarManager::FillVarESDtrack(const AliESDtrack *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);

  Double_t pidProbs[AliPID::kSPECIES];
  // Fill AliESDtrack interface specific information
  Double_t tpcNcls=particle->GetTPCNcls();
  values[AliDielectronVarManager::kNclsITS]       = particle->GetNcls(0); // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kNclsTPC]       = tpcNcls; // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kNclsTPCiter1]  = particle->GetTPCNclsIter1(); // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kNFclsTPC]      = particle->GetTPCNclsF();
  values[AliDielectronVarManager::kNFclsTPCr]     = particle->GetTPCClusterInfo(2,1);
  values[AliDielectronVarManager::kNFclsTPCrFrac] = particle->GetTPCClusterInfo(2);
  values[AliDielectronVarManager::kTPCsignalN]    = particle->GetTPCsignalN();
  values[AliDielectronVarManager::kNclsTRD]       = particle->GetNcls(2); // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kTRDntracklets] = particle->GetTRDntracklets(); // TODO: GetTRDtracklets/GetTRDntracklets?
  values[AliDielectronVarManager::kTRDpidQuality] = particle->GetTRDpidQuality();
  values[AliDielectronVarManager::kTrackStatus]   = (Double_t)particle->GetStatus();
  
  
  values[AliDielectronVarManager::kTPCchi2Cl] = -1;
  if (tpcNcls>0) values[AliDielectronVarManager::kTPCchi2Cl] = particle->GetTPCchi2() / tpcNcls;
  //TRD pidProbs
  particle->GetTRDpid(pidProbs);
  values[AliDielectronVarManager::kTRDprobEle]    = pidProbs[AliPID::kElectron];
  values[AliDielectronVarManager::kTRDprobPio]    = pidProbs[AliPID::kPion];

  values[AliDielectronVarManager::kKinkIndex0]    = particle->GetKinkIndex(0);
  
  Float_t impactParXY, impactParZ;
  particle->GetImpactParameters(impactParXY, impactParZ);
  values[AliDielectronVarManager::kImpactParXY]   = impactParXY;
  values[AliDielectronVarManager::kImpactParZ]    = impactParZ;


  values[AliDielectronVarManager::kPdgCode]=0;
  values[AliDielectronVarManager::kPdgCodeMother]=0;

  values[AliDielectronVarManager::kNumberOfDaughters]=-999;
  
  AliDielectronMC *mc=AliDielectronMC::Instance();
  
  if (mc->HasMC()){
    if (mc->GetMCTrack(particle))
      values[AliDielectronVarManager::kPdgCode]=
      mc->GetMCTrack(particle)->PdgCode();
    
    Int_t pdgMother=mc->GetMotherPDG(particle);
    if (pdgMother!=-999)
      values[AliDielectronVarManager::kPdgCodeMother]=pdgMother;

    values[AliDielectronVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
  } //if(mc->HasMC())
  


  values[AliDielectronVarManager::kITSsignal]   =   particle->GetITSsignal();
  
  Double_t itsdEdx[4];
  particle->GetITSdEdxSamples(itsdEdx);

  values[AliDielectronVarManager::kITSsignalSSD1]   =   itsdEdx[0];
  values[AliDielectronVarManager::kITSsignalSSD2]   =   itsdEdx[1];
  values[AliDielectronVarManager::kITSsignalSDD1]   =   itsdEdx[2];
  values[AliDielectronVarManager::kITSsignalSDD2]   =   itsdEdx[3];
  values[AliDielectronVarManager::kITSclusterMap]   =   particle->GetITSClusterMap();
  
  values[AliDielectronVarManager::kTrackLength]   = particle->GetIntegratedLength();
  //dEdx information
  Double_t mom = particle->GetP();
  const AliExternalTrackParam *in=particle->GetInnerParam();
  if (in) mom = in->GetP();
  values[AliDielectronVarManager::kPIn]=mom;
  values[AliDielectronVarManager::kTPCsignal]=particle->GetTPCsignal();
  // nsigma to Electron band
  // TODO: for the moment we set the bethe bloch parameters manually
  //       this should be changed in future!
  
  values[AliDielectronVarManager::kTPCnSigmaEle]=fgESDpid->NumberOfSigmasTPC(particle,AliPID::kElectron)-AliDielectronPID::GetCorrVal();
  values[AliDielectronVarManager::kTPCnSigmaPio]=fgESDpid->NumberOfSigmasTPC(particle,AliPID::kPion);
  values[AliDielectronVarManager::kTPCnSigmaMuo]=fgESDpid->NumberOfSigmasTPC(particle,AliPID::kMuon);
  values[AliDielectronVarManager::kTPCnSigmaKao]=fgESDpid->NumberOfSigmasTPC(particle,AliPID::kKaon);
  values[AliDielectronVarManager::kTPCnSigmaPro]=fgESDpid->NumberOfSigmasTPC(particle,AliPID::kProton);

  values[AliDielectronVarManager::kITSnSigmaEle]=fgESDpid->NumberOfSigmasITS(particle,AliPID::kElectron);
  values[AliDielectronVarManager::kITSnSigmaPio]=fgESDpid->NumberOfSigmasITS(particle,AliPID::kPion);
  values[AliDielectronVarManager::kITSnSigmaMuo]=fgESDpid->NumberOfSigmasITS(particle,AliPID::kMuon);
  values[AliDielectronVarManager::kITSnSigmaKao]=fgESDpid->NumberOfSigmasITS(particle,AliPID::kKaon);
  values[AliDielectronVarManager::kITSnSigmaPro]=fgESDpid->NumberOfSigmasITS(particle,AliPID::kProton);

  Double_t t0=fgESDpid->GetTOFResponse().GetTimeZero();
  values[AliDielectronVarManager::kTOFnSigmaEle]=fgESDpid->NumberOfSigmasTOF(particle,AliPID::kElectron,t0);
  values[AliDielectronVarManager::kTOFnSigmaPio]=fgESDpid->NumberOfSigmasTOF(particle,AliPID::kPion,t0);
  values[AliDielectronVarManager::kTOFnSigmaMuo]=fgESDpid->NumberOfSigmasTOF(particle,AliPID::kMuon,t0);
  values[AliDielectronVarManager::kTOFnSigmaKao]=fgESDpid->NumberOfSigmasTOF(particle,AliPID::kKaon,t0);
  values[AliDielectronVarManager::kTOFnSigmaPro]=fgESDpid->NumberOfSigmasTOF(particle,AliPID::kProton,t0);
}

inline void AliDielectronVarManager::FillVarAODTrack(const AliAODTrack *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);
  
  // Reset AliESDtrack interface specific information
  values[AliDielectronVarManager::kNclsITS]       = 0;
  values[AliDielectronVarManager::kNclsTPC]       = particle->GetTPCNcls();
  values[AliDielectronVarManager::kNclsTPCiter1]  = particle->GetTPCNcls(); // not really available in AOD
  values[AliDielectronVarManager::kNFclsTPC]      = 0;
  values[AliDielectronVarManager::kNFclsTPCr]     = 0;
  values[AliDielectronVarManager::kNFclsTPCrFrac] = 0;
  values[AliDielectronVarManager::kNclsTRD]       = 0;
  values[AliDielectronVarManager::kTRDntracklets] = 0;
  values[AliDielectronVarManager::kTRDpidQuality] = 0;
  
  values[AliDielectronVarManager::kTPCchi2Cl] = -1;
  values[AliDielectronVarManager::kTrackStatus]   = (Double_t)particle->GetStatus();
  
  //TRD pidProbs
  //TODO: set correctly
  values[AliDielectronVarManager::kTRDprobEle]    = 0;
  values[AliDielectronVarManager::kTRDprobPio]    = 0;
  
  //TODO: This is only an approximation!!!
  values[AliDielectronVarManager::kTPCsignalN]    = 0;
  
  // Fill AliAODTrack interface information
  //
  Double_t d0z0[2];
  GetDCA(particle, d0z0);
  values[AliDielectronVarManager::kImpactParXY]   = d0z0[0];
  values[AliDielectronVarManager::kImpactParZ]    = d0z0[1];

  values[AliDielectronVarManager::kPIn]=0;
  values[AliDielectronVarManager::kTPCsignal]=0;

  values[AliDielectronVarManager::kTPCnSigmaEle]=0;
  values[AliDielectronVarManager::kTPCnSigmaPio]=0;
  values[AliDielectronVarManager::kTPCnSigmaMuo]=0;
  values[AliDielectronVarManager::kTPCnSigmaKao]=0;
  values[AliDielectronVarManager::kTPCnSigmaPro]=0;

  values[AliDielectronVarManager::kITSclusterMap]   =   particle->GetITSClusterMap();
  
  AliAODPid *pid=particle->GetDetPid();
  if (pid){
    Double_t mom =pid->GetTPCmomentum();
    values[AliDielectronVarManager::kTPCsignalN] = pid->GetTPCsignalN();
    Double_t tpcNsigmaEle=fgAODpidUtil->NumberOfSigmasTPC(particle,AliPID::kElectron);
    Double_t tpcNsigmaPio=fgAODpidUtil->NumberOfSigmasTPC(particle,AliPID::kPion);
    Double_t tpcNsigmaMuo=fgAODpidUtil->NumberOfSigmasTPC(particle,AliPID::kMuon);
    Double_t tpcNsigmaKao=fgAODpidUtil->NumberOfSigmasTPC(particle,AliPID::kKaon);
    Double_t tpcNsigmaPro=fgAODpidUtil->NumberOfSigmasTPC(particle,AliPID::kProton);
    
    values[AliDielectronVarManager::kPIn]=mom;
    values[AliDielectronVarManager::kTPCsignal]=pid->GetTPCsignal();

    values[AliDielectronVarManager::kTPCnSigmaEle]=tpcNsigmaEle;
    values[AliDielectronVarManager::kTPCnSigmaPio]=tpcNsigmaPio;
    values[AliDielectronVarManager::kTPCnSigmaMuo]=tpcNsigmaMuo;
    values[AliDielectronVarManager::kTPCnSigmaKao]=tpcNsigmaKao;
    values[AliDielectronVarManager::kTPCnSigmaPro]=tpcNsigmaPro;
    
    values[AliDielectronVarManager::kTRDntracklets] = 0;
    values[AliDielectronVarManager::kTRDpidQuality] = 0;
    
  }
}

inline void AliDielectronVarManager::FillVarMCParticle(const AliMCParticle *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  values[AliDielectronVarManager::kNclsITS]       = 0;
  values[AliDielectronVarManager::kNclsTPC]       = 0;
  values[AliDielectronVarManager::kNclsTPCiter1]  = 0; 
  values[AliDielectronVarManager::kNFclsTPC]      = 0;
  values[AliDielectronVarManager::kNFclsTPCr]     = 0;
  values[AliDielectronVarManager::kNFclsTPCrFrac] = 0;
  values[AliDielectronVarManager::kNclsTRD]       = 0;
  values[AliDielectronVarManager::kTRDntracklets] = 0;
  values[AliDielectronVarManager::kTRDpidQuality] = 0;
  values[AliDielectronVarManager::kTPCchi2Cl]     = 0;
  values[AliDielectronVarManager::kTrackStatus]   = 0;
  values[AliDielectronVarManager::kTRDprobEle]    = 0;
  values[AliDielectronVarManager::kTRDprobPio]    = 0;
  values[AliDielectronVarManager::kTPCsignalN]    = 0;
  values[AliDielectronVarManager::kImpactParXY]   = 0;
  values[AliDielectronVarManager::kImpactParZ]    = 0;
  values[AliDielectronVarManager::kPIn]           = 0;
  values[AliDielectronVarManager::kTPCsignal]     = 0;
  values[AliDielectronVarManager::kTPCnSigmaEle]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPio]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaMuo]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaKao]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPro]  = 0;
  values[AliDielectronVarManager::kITSclusterMap] = 0;
  
  values[AliDielectronVarManager::kPdgCode]       = 0;
  values[AliDielectronVarManager::kPdgCodeMother] = 0;
  
  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);
  
  AliDielectronMC *mc=AliDielectronMC::Instance();

  // Fill AliMCParticle interface specific information
  values[AliDielectronVarManager::kPdgCode] = particle->PdgCode();

  AliMCParticle *mother = mc->GetMCTrackMother(particle);
  if (mother) values[AliDielectronVarManager::kPdgCodeMother] = mother->PdgCode();

  values[AliDielectronVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(particle);
  values[AliDielectronVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
}

inline void AliDielectronVarManager::FillVarAODMCParticle(const AliAODMCParticle *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  values[AliDielectronVarManager::kNclsITS]       = 0;
  values[AliDielectronVarManager::kNclsTPC]       = 0;
  values[AliDielectronVarManager::kNclsTPCiter1]  = 0;
  values[AliDielectronVarManager::kNFclsTPC]      = 0;
  values[AliDielectronVarManager::kNclsTRD]       = 0;
  values[AliDielectronVarManager::kTRDntracklets] = 0;
  values[AliDielectronVarManager::kTRDpidQuality] = 0;
  values[AliDielectronVarManager::kTPCchi2Cl]     = 0;
  values[AliDielectronVarManager::kTrackStatus]   = 0;
  values[AliDielectronVarManager::kTRDprobEle]    = 0;
  values[AliDielectronVarManager::kTRDprobPio]    = 0;
  values[AliDielectronVarManager::kTPCsignalN]    = 0;
  values[AliDielectronVarManager::kImpactParXY]   = 0;
  values[AliDielectronVarManager::kImpactParZ]    = 0;
  values[AliDielectronVarManager::kPIn]           = 0;
  values[AliDielectronVarManager::kTPCsignal]     = 0;
  values[AliDielectronVarManager::kTPCnSigmaEle]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPio]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaMuo]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaKao]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPro]  = 0;
  values[AliDielectronVarManager::kITSclusterMap] = 0;
  
  values[AliDielectronVarManager::kPdgCode]       = 0;
  values[AliDielectronVarManager::kPdgCodeMother] = 0;
  
  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);
  
  AliDielectronMC *mc=AliDielectronMC::Instance();


  // Fill AliAODMCParticle interface specific information
  values[AliDielectronVarManager::kPdgCode] = particle->PdgCode();

  AliVParticle *mother = mc->GetMCTrackMother(particle);
  if (mother) values[AliDielectronVarManager::kPdgCodeMother] = mother->PdgCode();
  values[AliDielectronVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(particle);

  values[AliDielectronVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
}

inline void AliDielectronVarManager::FillVarDielectronPair(const AliDielectronPair *pair, Double_t * const values)
{
  //
  // Fill pair information available for histogramming into an array
  //
  
  values[AliDielectronVarManager::kPdgCode]=0;
  values[AliDielectronVarManager::kPdgCodeMother]=0;
  
  // Fill common AliVParticle interface information
  FillVarVParticle(pair, values);

  // Fill AliDielectronPair specific information
  const AliKFParticle &kfPair = pair->GetKFParticle();

  Double_t thetaHE=0;
  Double_t phiHE=0;
  Double_t thetaCS=0;
  Double_t phiCS=0;

  pair->GetThetaPhiCM(thetaHE,phiHE,thetaCS,phiCS);
    
  values[AliDielectronVarManager::kChi2NDF]      = kfPair.GetChi2()/kfPair.GetNDF();
  values[AliDielectronVarManager::kDecayLength]  = kfPair.GetDecayLength();
  values[AliDielectronVarManager::kR]            = kfPair.GetR();
  values[AliDielectronVarManager::kOpeningAngle] = pair->OpeningAngle();
  values[AliDielectronVarManager::kThetaHE]      = thetaHE;
  values[AliDielectronVarManager::kPhiHE]        = phiHE;
  values[AliDielectronVarManager::kThetaCS]      = thetaCS;
  values[AliDielectronVarManager::kPhiCS]        = phiCS;
  values[AliDielectronVarManager::kLegDist]      = pair->DistanceDaughters();
  values[AliDielectronVarManager::kLegDistXY]    = pair->DistanceDaughtersXY();
  values[AliDielectronVarManager::kDeltaEta]     = pair->DeltaEta();
  values[AliDielectronVarManager::kDeltaPhi]     = pair->DeltaPhi();
  values[AliDielectronVarManager::kMerr]         = kfPair.GetErrMass()>1e-30&&kfPair.GetMass()>1e-30?kfPair.GetErrMass()/kfPair.GetMass():1000000;
  values[AliDielectronVarManager::kPairType]     = pair->GetType();


  
  AliDielectronMC *mc=AliDielectronMC::Instance();
  
  if (mc->HasMC()){
  Bool_t samemother =  mc->HaveSameMother(pair);
  values[AliDielectronVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(pair);
  values[AliDielectronVarManager::kHaveSameMother] = samemother ;
  }//if (mc->HasMC())


}


inline void AliDielectronVarManager::FillVarVEvent(const AliVEvent *event, Double_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //
  const AliVVertex *primVtx = event->GetPrimaryVertex();
  
  values[AliDielectronVarManager::kXvPrim]       = 0;
  values[AliDielectronVarManager::kYvPrim]       = 0;
  values[AliDielectronVarManager::kZvPrim]       = 0;
  values[AliDielectronVarManager::kChi2NDF]      = 0;
  
  values[AliDielectronVarManager::kNTrk]         = 0;
  values[AliDielectronVarManager::kNevents]      = 0; //always fill bin 0;

  if (!primVtx) return;
  
  values[AliDielectronVarManager::kXvPrim]       = primVtx->GetX();
  values[AliDielectronVarManager::kYvPrim]       = primVtx->GetY();
  values[AliDielectronVarManager::kZvPrim]       = primVtx->GetZ();
  values[AliDielectronVarManager::kChi2NDF]      = primVtx->GetChi2perNDF();

  values[AliDielectronVarManager::kNTrk]         = event->GetNumberOfTracks();
}

inline void AliDielectronVarManager::FillVarESDEvent(const AliESDEvent *event, Double_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  // 
  
  // Fill common AliVEvent interface information
  FillVarVEvent(event, values);

  Double_t centralityF=-1;
  AliCentrality *esdCentrality = const_cast<AliESDEvent*>(event)->GetCentrality();
  if (esdCentrality) centralityF = esdCentrality->GetCentralityPercentile("V0M");
  
  // Fill AliESDEvent interface specific information
  const AliESDVertex *primVtx = event->GetPrimaryVertex();
  values[AliDielectronVarManager::kXRes]       = primVtx->GetXRes();
  values[AliDielectronVarManager::kYRes]       = primVtx->GetYRes();
  values[AliDielectronVarManager::kZRes]       = primVtx->GetZRes();
  values[AliDielectronVarManager::kCentrality] = centralityF;
}
  
inline void AliDielectronVarManager::FillVarAODEvent(const AliAODEvent *event, Double_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //   

  // Fill common AliVEvent interface information
  FillVarVEvent(event, values);

  // Fill AliAODEvent interface specific information
}
  
inline void AliDielectronVarManager::FillVarMCEvent(const AliMCEvent *event, Double_t * const values)
{ 
  //
  // Fill event information available for histogramming into an array
  //   
        
  // Fill common AliVEvent interface information
  FillVarVEvent(event, values);

  // Fill AliMCEvent interface specific information
} 
  
inline void AliDielectronVarManager::InitESDpid(Int_t type)
{
  //
  // initialize PID parameters
  // type=0 is simulation
  // type=1 is data

  if (!fgESDpid) fgESDpid=new AliESDpid;
  Double_t alephParameters[5];
  // simulation
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;
  fgESDpid->GetTOFResponse().SetTimeResolution(80.);
  
  // data
  if (type==1){    
    alephParameters[0] = 0.0283086/0.97;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;
    fgESDpid->GetTOFResponse().SetTimeResolution(130.);
    fgESDpid->GetTPCResponse().SetMip(50.);
  }

  fgESDpid->GetTPCResponse().SetBetheBlochParameters(
    alephParameters[0],alephParameters[1],alephParameters[2],
    alephParameters[3],alephParameters[4]);
  
  fgESDpid->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
}

inline void AliDielectronVarManager::InitAODpidUtil(Int_t type)
{
  if (!fgAODpidUtil) fgAODpidUtil=new AliAODpidUtil;
  Double_t alephParameters[5];
  // simulation
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;
  fgAODpidUtil->GetTOFResponse().SetTimeResolution(80.);
  
  // data
  if (type==1){
    alephParameters[0] = 0.0283086/0.97;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;
    fgAODpidUtil->GetTOFResponse().SetTimeResolution(130.);
    fgAODpidUtil->GetTPCResponse().SetMip(50.);
  }
  
  fgAODpidUtil->GetTPCResponse().SetBetheBlochParameters(
    alephParameters[0],alephParameters[1],alephParameters[2],
    alephParameters[3],alephParameters[4]);
  
  fgAODpidUtil->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
}


inline void AliDielectronVarManager::SetEvent(AliVEvent * const ev)
{
  
  fgEvent = ev;
  if (fgKFVertex) delete fgKFVertex;
  fgKFVertex=0x0;
  if (ev && ev->GetPrimaryVertex()) fgKFVertex=new AliKFVertex(*ev->GetPrimaryVertex());
}


inline Bool_t AliDielectronVarManager::GetDCA(const AliAODTrack *track, Double_t d0z0[2])
{
  if(track->TestBit(AliAODTrack::kIsDCA)){
    d0z0[0]=track->DCA();
    d0z0[1]=track->ZAtDCA();
    return kTRUE;
  }
  
  Double_t covd0z0[3];
  AliAODTrack copy(*track);
  AliAODVertex *vtx =(AliAODVertex*)(fgEvent->GetPrimaryVertex());
  Double_t fBzkG = fgEvent->GetMagneticField(); // z componenent of field in kG
  Bool_t ok = copy.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
  if(!ok){
    d0z0[0]=-999.;
    d0z0[1]=-999.;
  }
  return ok;
}

/*
inline void AliDielectronVarManager::FillValues(const TParticle *particle, Double_t *values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill TParticle interface information
  values[AliDielectronVarManager::kPx]     = particle->Px();
  values[AliDielectronVarManager::kPy]     = particle->Py();
  values[AliDielectronVarManager::kPz]     = particle->Pz();
  values[AliDielectronVarManager::kPt]     = particle->Pt();
  values[AliDielectronVarManager::kP]      = particle->P();

  values[AliDielectronVarManager::kXv]     = particle->Vx();
  values[AliDielectronVarManager::kYv]     = particle->Vy();
  values[AliDielectronVarManager::kZv]     = particle->Vz();

  values[AliDielectronVarManager::kOneOverPt] = 1./particle->Pt();
  values[AliDielectronVarManager::kPhi]    = particle->Phi();
  values[AliDielectronVarManager::kTheta]  = 
  values[AliDielectronVarManager::kEta]    = particle->Eta();
  values[AliDielectronVarManager::kY]      = 

  values[AliDielectronVarManager::kE]      = particle->Energy();
  values[AliDielectronVarManager::kM]      = particle->GetMass();

  values[AliDielectronVarManager::kCharge] = particle->GetPDG()->Charge()/3; // uggly

}*/

#endif

