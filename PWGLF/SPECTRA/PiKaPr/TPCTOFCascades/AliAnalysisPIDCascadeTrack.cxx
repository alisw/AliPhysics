#include "AliAnalysisPIDCascadeTrack.h"
#include "AliAnalysisPIDCascadeEvent.h"
//#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliMCEvent.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "AliPIDResponse.h"
#include "AliLog.h"

ClassImp(AliAnalysisPIDCascadeTrack)

//___________________________________________________________

Float_t AliAnalysisPIDCascadeTrack::fgMatchTrackDeltaX = 10.;
Float_t AliAnalysisPIDCascadeTrack::fgMatchTrackDeltaZ = 10;

//___________________________________________________________

TLorentzVector AliAnalysisPIDCascadeTrack::fgLorentzVector;
AliTPCPIDResponse *AliAnalysisPIDCascadeTrack::fgTPCResponse = NULL;
AliTOFPIDResponse *AliAnalysisPIDCascadeTrack::fgTOFResponse = NULL;
TH2F *AliAnalysisPIDCascadeTrack::hTOFtuned_th[AliPID::kSPECIES] = {
  NULL, NULL, NULL, NULL, NULL
};


Float_t
AliAnalysisPIDCascadeTrack::GetY(Float_t mass) const
{
  
  fgLorentzVector.SetPtEtaPhiM(fPt, fEta, fPhi, mass);
  return fgLorentzVector.Rapidity();
}

//___________________________________________________________

AliAnalysisPIDCascadeTrack::AliAnalysisPIDCascadeTrack() :
  //General
  TObject(),
  fP(),
  fPt(),
  fEta(),
  fPhi(),
  fSign(0.),
  fStatus(0x0),
  fLabel(0),
  fImpactParameter(),
  //TPC
  fTPCdEdx(0.),
  fTPCdEdxN(0),
  //TOF
  fTOFIndex(0),
  fTOFLength(0.),
  fTOFDeltaX(0.),
  fTOFDeltaZ(0.),
  fTOFLabel(),
  //MC
  fMCPrimary(kFALSE),
  fMCPdgCode(0),
  fMCMotherPrimary(kFALSE),
  fMCMotherPdgCode(0),
  fMCMotherLabel(0),
  fMCPrimaryPdgCode(0),
  fMCPrimaryLabel(0),
  fMCTOFMatchPrimary(kFALSE),
  fMCTOFMatchPdgCode(0),
  fMCTOFMatchLevel(-1),
  fMCTOFTime(0.),
  fMCTOFLength(0.),
  fMCSecondaryWeak(kFALSE),
  fMCSecondaryMaterial(kFALSE),
//PID
  nSigmaPionTPC(0.),
  nSigmaKaonTPC(0.),
  nSigmaProtonTPC(0.),
  nSigmaElectronTPC(0.),
  nSigmaPionTOF(0.),
  nSigmaKaonTOF(0.),
  nSigmaProtonTOF(0.),
  nSigmaElectronTOF(0.),
  fTrackCutFlag(0)
{
  /*
   * default constructor
   */

  Double_t bbParam[6] = { /* R+ fit on minimum-bias PbPb (run 138275, pass2) */
    5.22879e+01,
    2.80863e-02,
    2.58364e+01,
    5.15102e-07,
    2.42169e+00,
    5.38930e+00
  };
  
  if (!fgTPCResponse) {
    fgTPCResponse = new AliTPCPIDResponse();
    fgTPCResponse->SetMip(bbParam[0]);
    fgTPCResponse->SetBetheBlochParameters(bbParam[1], bbParam[2], bbParam[3], bbParam[4], bbParam[5]);
  }

  /* reset */
  Reset();
}

//___________________________________________________________
AliAnalysisPIDCascadeTrack::AliAnalysisPIDCascadeTrack(const AliAnalysisPIDCascadeTrack &source) :
  //General
  TObject(source),
  fP(source.fP),
  fPt(source.fPt),
  fEta(source.fEta),
  fPhi(source.fPhi),
  fSign(source.fSign),
  fStatus(source.fStatus),
  fLabel(source.fLabel),
  fImpactParameter(),
  //TPC
  fTPCdEdx(source.fTPCdEdx),
  fTPCdEdxN(source.fTPCdEdxN),
  //TOF
  fTOFIndex(source.fTOFIndex),
  fTOFLength(source.fTOFLength),
  fTOFDeltaX(source.fTOFDeltaX),
  fTOFDeltaZ(source.fTOFDeltaZ),
  fTOFLabel(),
  //MC
  fMCPrimary(source.fMCPrimary),
  fMCPdgCode(source.fMCPdgCode),
  fMCMotherPrimary(source.fMCMotherPrimary),
  fMCMotherPdgCode(source.fMCMotherPdgCode),
  fMCMotherLabel(source.fMCMotherLabel),
  fMCPrimaryPdgCode(source.fMCPrimaryPdgCode),
  fMCPrimaryLabel(source.fMCPrimaryLabel),
  fMCTOFMatchPrimary(source.fMCTOFMatchPrimary),
  fMCTOFMatchPdgCode(source.fMCTOFMatchPdgCode),
  fMCTOFMatchLevel(source.fMCTOFMatchLevel),
  fMCTOFTime(source.fMCTOFTime),
  fMCTOFLength(source.fMCTOFLength),
  fMCSecondaryWeak(source.fMCSecondaryWeak),
  fMCSecondaryMaterial(source.fMCSecondaryMaterial),
//PID
  nSigmaPionTPC(source.nSigmaPionTPC),
  nSigmaKaonTPC(source.nSigmaKaonTPC),
  nSigmaProtonTPC(source.nSigmaProtonTPC),
  nSigmaElectronTPC(source.nSigmaElectronTPC),
  nSigmaPionTOF(source.nSigmaPionTOF),
  nSigmaKaonTOF(source.nSigmaKaonTOF),
  nSigmaProtonTOF(source.nSigmaProtonTOF),
  nSigmaElectronTOF(source.nSigmaElectronTOF),
  fTrackCutFlag(source.fTrackCutFlag)
{
  /*
   * copy constructor
   */

  for (Int_t i = 0; i < 2; i++) fImpactParameter[i] = source.fImpactParameter[i];
  for (Int_t i = 0; i < 3; i++) fTOFLabel[i] = source.fTOFLabel[i];

}

//___________________________________________________________

AliAnalysisPIDCascadeTrack &
AliAnalysisPIDCascadeTrack::operator=(const AliAnalysisPIDCascadeTrack &source)
{
  /*
   * operator=
   */

  if (&source == this) return *this;
  TObject::operator=(source);
  //General
  fP = source.fP;
  fPt = source.fPt;
  fEta = source.fEta;
  fPhi = source.fPhi;
  fSign = source.fSign;
  fStatus = source.fStatus;
  fLabel = source.fLabel;
  for (Int_t i = 0; i < 2; i++) fImpactParameter[i] = source.fImpactParameter[i];
  //TPC
  fTPCdEdx = source.fTPCdEdx;
  fTPCdEdxN = source.fTPCdEdxN;
  //TOF
  fTOFIndex = source.fTOFIndex;
  fTOFLength = source.fTOFLength;
  fTOFDeltaX = source.fTOFDeltaX;
  fTOFDeltaZ = source.fTOFDeltaZ;
  for (Int_t i = 0; i < 3; i++) fTOFLabel[i] = source.fTOFLabel[i];
  //MC
  fMCPrimary = source.fMCPrimary;
  fMCPdgCode = source.fMCPdgCode;
  fMCMotherPrimary = source.fMCMotherPrimary;
  fMCMotherPdgCode = source.fMCMotherPdgCode;
  fMCMotherLabel = source.fMCMotherLabel;
  fMCPrimaryPdgCode = source.fMCPrimaryPdgCode;
  fMCPrimaryLabel = source.fMCPrimaryLabel;
  fMCTOFMatchPrimary = source.fMCTOFMatchPrimary;
  fMCTOFMatchPdgCode = source.fMCTOFMatchPdgCode;
  fMCTOFMatchLevel = source.fMCTOFMatchLevel;
  fMCTOFTime = source.fMCTOFTime;
  fMCTOFLength = source.fMCTOFLength;
  fMCSecondaryWeak = source.fMCSecondaryWeak;
  fMCSecondaryMaterial = source.fMCSecondaryMaterial;
  //PID
  nSigmaPionTPC = source.nSigmaPionTPC;
  nSigmaKaonTPC = source.nSigmaKaonTPC;
  nSigmaProtonTPC = source.nSigmaProtonTPC;
  nSigmaElectronTPC = source.nSigmaElectronTPC;
  nSigmaPionTOF = source.nSigmaPionTOF;
  nSigmaKaonTOF = source.nSigmaKaonTOF;
  nSigmaProtonTOF = source.nSigmaProtonTOF;
  nSigmaElectronTOF = source.nSigmaElectronTOF;
  fTrackCutFlag = source.fTrackCutFlag;

  return *this;
}

//___________________________________________________________

AliAnalysisPIDCascadeTrack::~AliAnalysisPIDCascadeTrack()
{
  /*
   * default destructor
   */

}

void
AliAnalysisPIDCascadeTrack::Reset()
{
  /*
   * reset
   */
  //General
  fP = 0.;
  fPt = 0.;
  fEta = 0.;
  fPhi = 0.;
  fSign = 0.;
  fStatus = 0;
  fLabel = 0;
  for (Int_t i = 0; i < 2; i++) fImpactParameter[i] = 0.;
  //TPC
  fTPCdEdx = 0.;
  fTPCdEdxN = 0;
  //TOF
  fTOFIndex = 0;
  fTOFLength = 0.;
  fTOFDeltaX = 0.;
  fTOFDeltaZ = 0.;
  for (Int_t i = 0; i < 3; i++) fTOFLabel[i] = 0;
  //MC
  fMCPrimary = kFALSE;
  fMCPdgCode = 0;
  fMCMotherPrimary = kFALSE;
  fMCMotherPdgCode = 0;
  fMCMotherLabel = 0;
  fMCPrimaryPdgCode = 0;
  fMCPrimaryLabel = 0;
  fMCTOFMatchPrimary = kFALSE;
  fMCTOFMatchPdgCode = 0;
  fMCTOFMatchLevel = -1;
  fMCTOFTime = 0.;
  fMCTOFLength = 0.;
  fMCSecondaryWeak = 0.;
  fMCSecondaryMaterial = 0.;
  //PID
  nSigmaPionTPC = 0.;
  nSigmaKaonTPC = 0;
  nSigmaProtonTPC = 0;
  nSigmaElectronTPC = 0;
  nSigmaPionTOF = 0;
  nSigmaKaonTOF = 0;
  nSigmaProtonTOF = 0;
  nSigmaElectronTOF = 0;
  fTrackCutFlag = 0;
  
}

//___________________________________________________________

void
AliAnalysisPIDCascadeTrack::Update(AliESDtrack *track, AliMCEvent *mcevent, AliPIDResponse *PIDRes, Int_t TrackCutFlag)
{
  /*
   * update
   */

  fP = track->P();
  fPt = track->Pt();
  fEta = track->Eta();
  fPhi = track->Phi();
  fSign = track->GetSign();
  fStatus = track->GetStatus();
  fLabel = track->GetLabel();
  track->GetImpactParameters(&fImpactParameter[0],&fImpactParameter[1]);
  fTPCdEdx = track->GetTPCsignal();
  fTPCdEdxN = track->GetTPCsignalN();

  nSigmaPionTPC = PIDRes->NumberOfSigmasTPC(track,AliPID::kPion);
  nSigmaKaonTPC = PIDRes->NumberOfSigmasTPC(track,AliPID::kKaon);
  nSigmaProtonTPC = PIDRes->NumberOfSigmasTPC(track,AliPID::kProton);
  nSigmaElectronTPC = PIDRes->NumberOfSigmasTPC(track,AliPID::kElectron);
  nSigmaPionTOF = PIDRes->NumberOfSigmasTOF(track,AliPID::kPion);
  nSigmaKaonTOF = PIDRes->NumberOfSigmasTOF(track,AliPID::kKaon);
  nSigmaProtonTOF = PIDRes->NumberOfSigmasTOF(track,AliPID::kProton);
  nSigmaElectronTOF = PIDRes->NumberOfSigmasTOF(track,AliPID::kElectron);

  fTOFLength = track->GetIntegratedLength();
  fTOFDeltaX = track->GetTOFsignalDx();
  fTOFDeltaZ = track->GetTOFsignalDz();
  track->GetTOFLabel(fTOFLabel);
  fTrackCutFlag = TrackCutFlag;
  /* info from track references */
  fMCTOFTime = 0.;
  fMCTOFLength = 0.;
  if (mcevent && fTOFLabel[0] > 0) {
    TParticle *particle;
    TClonesArray *arrayTR;
    AliTrackReference *trackRef;
    mcevent->GetParticleAndTR(fTOFLabel[0], particle, arrayTR);
    if(arrayTR) {
      for (Int_t itr = 0; itr < arrayTR->GetEntries(); itr++) {
	trackRef = (AliTrackReference *)arrayTR->At(itr);
	if (!trackRef || trackRef->DetectorId() != AliTrackReference::kTOF) continue;
	fMCTOFTime = trackRef->GetTime() * 1.e12; /* s -> ps */
	fMCTOFLength = trackRef->GetLength();
	/* break as soon as we get it */
	break;
      }
    }
  }
  /* info from stack */
  fMCPrimary = kFALSE;
  fMCPdgCode = 0;
  fMCMotherPrimary = kFALSE;
  fMCMotherPdgCode = 0;
  fMCPrimaryPdgCode = 0;
  fMCSecondaryWeak = kFALSE;
  fMCSecondaryMaterial = kFALSE;
  if (mcevent) {
    Int_t index = TMath::Abs(fLabel);
    if (index < 0) {
      printf("index = %d\n", index);
      return;
    }
    TParticle *particle = mcevent->Particle(index);
    fMCPrimary = mcevent->IsPhysicalPrimary(index);
    fMCSecondaryWeak = mcevent->IsSecondaryFromWeakDecay(index);
    fMCSecondaryMaterial = mcevent->IsSecondaryFromMaterial(index);
    fMCPdgCode = particle->GetPdgCode();
    Int_t indexm = particle->GetFirstMother();
    if (indexm < 0) {
      fMCMotherPrimary = kFALSE;
      fMCMotherPdgCode = 0;
    }
    else {
      TParticle *particlem = mcevent->Particle(indexm);
      fMCMotherPrimary = mcevent->IsPhysicalPrimary(indexm);
      fMCMotherPdgCode = particlem->GetPdgCode();
      fMCMotherLabel = indexm;
    }
    if (fMCPrimary) {
      fMCPrimaryPdgCode = fMCPdgCode;
      fMCPrimaryLabel = fLabel;
    } else {
      Bool_t primary = kFALSE;
      while (!primary) {
	if (indexm < 0) {
	  fMCPrimaryPdgCode = 0;
	  break;
	}
	TParticle *particlem = mcevent->Particle(indexm);
	primary = mcevent->IsPhysicalPrimary(indexm);
	if (primary) {
	  fMCPrimaryPdgCode = particlem->GetPdgCode();
	  fMCPrimaryLabel = indexm;
	} else
	  indexm = particlem->GetFirstMother();
      }
    }
    
    /* check TOF match */
    fMCTOFMatchPrimary = kFALSE;
    fMCTOFMatchPdgCode = 0;
    fMCTOFMatchLevel = -1;
    if (fTOFLabel[0] > 0) {
      index = fTOFLabel[0];
      particle = mcevent->Particle(index);
      fMCTOFMatchPrimary = mcevent->IsPhysicalPrimary(index);
      fMCTOFMatchPdgCode = particle->GetPdgCode();
      Int_t tracklabel = TMath::Abs(fLabel);
      Int_t matchlevel = -1;
      for (Int_t ilevel = 0; index > 0; ilevel++) {
	if (index == tracklabel) {
	  matchlevel = ilevel;
	  break;
	}
	index = mcevent->Particle(index)->GetFirstMother();
      }
      fMCTOFMatchLevel = matchlevel;
    }
    
  }
  fTimeZeroSigma = 0.;
  
}

//___________________________________________________________

Bool_t
AliAnalysisPIDCascadeTrack::HasTPCPID() const
{
  /*
   * has TPC PID
   */

  /* check PID signal */
  if (fTPCdEdx <= 0. || fTPCdEdxN == 0) return kFALSE;
  return kTRUE;
}
  
//___________________________________________________________

Bool_t
AliAnalysisPIDCascadeTrack::HasTOFPID(TH1 *henabled) const
{
  /*
   * has TOF PID
   */

  /* check channel enabled */
  if (henabled && henabled->GetBinContent(fTOFIndex + 1) == 0) return kFALSE;
  
  /* check TOF matched track */
  if (!(fStatus & AliESDtrack::kTOFout)||
      !(fStatus & AliESDtrack::kTIME)) return kFALSE;
  /* check integrated length */
  if (fTOFLength < 350.) return kFALSE;
  /* check deltax and deltaz */
  if (TMath::Abs(fTOFDeltaX) > fgMatchTrackDeltaX ||
      TMath::Abs(fTOFDeltaZ) > fgMatchTrackDeltaZ) return kFALSE;
  return kTRUE;
}
  
//___________________________________________________________  
//___________________________________________________________
//___________________________________________________________
//___________________________________________________________

Int_t
AliAnalysisPIDCascadeTrack::GetMCPID() const
{
  /*
   * get MC PID
   */

  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
    if (TMath::Abs(fMCPdgCode) == AliPID::ParticleCode(ipart))
      return ipart;
  return -1;

}

//___________________________________________________________

Int_t
AliAnalysisPIDCascadeTrack::GetMCCharge() const
{
  /*
   * get MC charge
   */
  
  TDatabasePDG *dbpdg = TDatabasePDG::Instance();
  TParticlePDG *pdg = dbpdg->GetParticle(fMCPdgCode);
  if (!pdg) return 0;
  return (Int_t)TMath::Sign(1., pdg->Charge());
};

//___________________________________________________________
//___________________________________________________________

// Bool_t
// AliAnalysisPIDCascadeTrack::AcceptTrack(Bool_t selPrimaries)
// {
//   /*
//    * accept track
//    */

//   /* check charge */
//   if (fSign == 0.) return kFALSE;
//   /* check eta */
//   if (TMath::Abs(fEta) > fgEtaCut) return kFALSE;
//   if (TMath::Abs(fEta) < fgEtaReject) return kFALSE;
//   /* check max DCA to vertex Z */
//   if (TMath::Abs(fImpactParameter[1]) > fgMaxDCAToVertexZCut) return kFALSE;
//   /* check max DCA to vertex XY */
//   if (selPrimaries) {
//     Float_t maxdca = fgMaxDCAToVertexXYPtDepFormula->Eval(fPt);
//     if (TMath::Abs(fImpactParameter[0]) > maxdca) return kFALSE;
//   }
//   /* check max chi2 per cluster TPC */
//   Float_t chi2 = fTPCchi2 / (Float_t)fTPCNcls;
//   if (chi2 > fgMaxChi2PerClusterTPC) return kFALSE;
//   /* check min N clusters TPC */
//   if (fgAcceptTrackClusterCut != 1) {
//     if (fTPCNcls < fgMinNClustersTPC) return kFALSE;
//   }
//   if (fgAcceptTrackClusterCut == 1) {
//     /* check min N crossed rows TPC */
//     if (fTPCNcr < fgMinNCrossedRowsTPC) return kFALSE;
//     /* check min crossed rows over findable clusters ratio */
//     Float_t crratio = 1.;
//     if (fTPCNclsF > 0) crratio = fTPCNcr / fTPCNclsF;
//     if (crratio < fgMinRatioCrossedRowsOverFindableClustersTPC) return kFALSE;
//   }
//   /* check accept track status cut */
//   if (fgAcceptTrackStatusCut && (fStatus & fgAcceptTrackStatusCut) == 0) return kFALSE;
//   /* check reject track status cut */
//   if (fgRejectTrackStatusCut && (fStatus & fgRejectTrackStatusCut) != 0) return kFALSE;
//   /* reject ITS fakes if requested */
//   if (fgRejectITSFakes && fITSFakeFlag) return kFALSE;

//   /* accept track */
//   return kTRUE;
// }

// Bool_t AliAnalysisPIDCascadeTrack::CheckExtraCuts(Float_t minTPCNcr, Float_t maxChi2PerFindableCluster, Float_t maxDCAz) {
//   if(fTPCNcr<minTPCNcr) return kFALSE;
//   if(fTPCNclsF<=0) return kFALSE;
//   if(fTPCchi2/fTPCNclsF > maxChi2PerFindableCluster) return kFALSE;
//   if(TMath::Abs(fImpactParameter[1])>maxDCAz) return kFALSE;
//   return kTRUE;
// };
