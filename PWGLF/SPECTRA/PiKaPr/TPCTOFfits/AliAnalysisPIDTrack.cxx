#include "AliAnalysisPIDTrack.h"
#include "AliAnalysisPIDEvent.h"
#include "AliStack.h"
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

ClassImp(AliAnalysisPIDTrack)

//___________________________________________________________

Float_t AliAnalysisPIDTrack::fgEtaCut = 0.8;
Float_t AliAnalysisPIDTrack::fgEtaReject = 0.;
//TFormula *AliAnalysisPIDTrack::fgMaxDCAToVertexXYPtDepFormula = new TFormula("f1MaxDCAToVertexXYPtDepFormula", "0.0182+0.0350/x^1.01");
TFormula *AliAnalysisPIDTrack::fgMaxDCAToVertexXYPtDepFormula = new TFormula("f1MaxDCAToVertexXYPtDepFormula", "0.0105+0.0350/x^1.1");
UShort_t AliAnalysisPIDTrack::fgMinNClustersTPC = 70;
UShort_t AliAnalysisPIDTrack::fgMinNCrossedRowsTPC = 70;
Float_t AliAnalysisPIDTrack::fgMinRatioCrossedRowsOverFindableClustersTPC = 0.8;
//Int_t AliAnalysisPIDTrack::fgAcceptTrackClusterCut = 0;
Int_t AliAnalysisPIDTrack::fgAcceptTrackClusterCut = 1;
ULong_t AliAnalysisPIDTrack::fgAcceptTrackStatusCut = 0;
ULong_t AliAnalysisPIDTrack::fgRejectTrackStatusCut = 0;
Float_t AliAnalysisPIDTrack::fgMaxDCAToVertexZCut = 2.;
Float_t AliAnalysisPIDTrack::fgMaxChi2PerClusterTPC = 4.;
Bool_t AliAnalysisPIDTrack::fgRejectITSFakes = kFALSE;
Float_t AliAnalysisPIDTrack::fgMatchTrackDeltaX = 10.;
Float_t AliAnalysisPIDTrack::fgMatchTrackDeltaZ = 10.;

//___________________________________________________________

TLorentzVector AliAnalysisPIDTrack::fgLorentzVector;
//AliTOFGeometry AliAnalysisPIDTrack::fgTOFGeometry;
//AliTOFcalibHisto AliAnalysisPIDTrack::fgTOFcalibHisto;
//Bool_t AliAnalysisPIDTrack::fgTOFcalibHistoFlag = kFALSE;
AliTPCPIDResponse *AliAnalysisPIDTrack::fgTPCResponse = NULL;
AliTOFPIDResponse *AliAnalysisPIDTrack::fgTOFResponse = NULL;
TH2F *AliAnalysisPIDTrack::hTOFtuned_th[AliPID::kSPECIES] = {
  NULL, NULL, NULL, NULL, NULL
};

//___________________________________________________________

Float_t
AliAnalysisPIDTrack::GetY(Float_t mass) const
{
  
  fgLorentzVector.SetPtEtaPhiM(fPt, fEta, fPhi, mass);
  return fgLorentzVector.Rapidity();
}

//___________________________________________________________

AliAnalysisPIDTrack::AliAnalysisPIDTrack() :
  TObject(),
  fP(),
  fPt(),
  fEta(),
  fPhi(),
  fSign(0.),
  fStatus(0x0),
  fLabel(0),
  fImpactParameter(),
  //fImpactParameterCov(),
  fTPCmomentum(0.),
  fTPCdEdx(0.),
  fTPCdEdxN(0),
  fTPCNcls(0),
  fTPCNclsF(0),
  fTPCNcr(0.),
  fTOFIndex(0),
  fTOFTime(0.),
  fTOFExpTime(),
  fTOFLength(0.),
  fTOFDeltaX(0.),
  fTOFDeltaZ(0.),
  fTOFLabel(),
  fMCPrimary(kFALSE),
  fMCPdgCode(0),
  fMCMotherPrimary(kFALSE),
  fMCMotherPdgCode(0),
  fMCMotherLabel(0),
  fMCTOFMatchPrimary(kFALSE),
  fMCTOFMatchPdgCode(0),
  fMCTOFMatchLevel(-1),
  fMCTOFTime(0.),
  fMCTOFLength(0.),
  fMCSecondaryWeak(kFALSE),
  fMCSecondaryMaterial(kFALSE),
//fHMPIDmomentum(0.),
//fHMPIDsignal(0.),
  nSigmaPionTPC(0.),
  nSigmaKaonTPC(0.),
  nSigmaProtonTPC(0.),
  nSigmaPionTOF(0.),
  nSigmaKaonTOF(0.),
  nSigmaProtonTOF(0.),
  fTPCchi2(0.),
  fITSFakeFlag(kFALSE),
  fTrackCutFlag(0),
  fEMCalE(-999),
  fEMCalP(-999),
  fHasEMCal(kFALSE),
  fTimeZeroSigma(0.)
{
  /*
   * default constructor
   */

  /* load calib histo */
  //  if (!fgTOFcalibHistoFlag) {
  //    fgTOFcalibHisto.LoadCalibHisto();
  //    fgTOFcalibHistoFlag = kTRUE;
  //  }

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

AliAnalysisPIDTrack::AliAnalysisPIDTrack(const AliAnalysisPIDTrack &source) :
  TObject(source),
  fP(source.fP),
  fPt(source.fPt),
  fEta(source.fEta),
  fPhi(source.fPhi),
  fSign(source.fSign),
  fStatus(source.fStatus),
  fLabel(source.fLabel),
  fImpactParameter(),
  //fImpactParameterCov(),
  fTPCmomentum(source.fTPCmomentum),
  fTPCdEdx(source.fTPCdEdx),
  fTPCdEdxN(source.fTPCdEdxN),
  fTPCNcls(source.fTPCNcls),
  fTPCNclsF(source.fTPCNclsF),
  fTPCNcr(source.fTPCNcr),  
  fTOFIndex(source.fTOFIndex),
  fTOFTime(source.fTOFTime),
  fTOFExpTime(),
  fTOFLength(source.fTOFLength),
  fTOFDeltaX(source.fTOFDeltaX),
  fTOFDeltaZ(source.fTOFDeltaZ),
  fTOFLabel(),
  fMCPrimary(source.fMCPrimary),
  fMCPdgCode(source.fMCPdgCode),
  fMCMotherPrimary(source.fMCMotherPrimary),
  fMCMotherPdgCode(source.fMCMotherPdgCode),
  fMCMotherLabel(source.fMCMotherLabel),
  fMCTOFMatchPrimary(source.fMCTOFMatchPrimary),
  fMCTOFMatchPdgCode(source.fMCTOFMatchPdgCode),
  fMCTOFMatchLevel(source.fMCTOFMatchLevel),
  fMCTOFTime(source.fMCTOFTime),
  fMCTOFLength(source.fMCTOFLength),
  fMCSecondaryWeak(source.fMCSecondaryWeak),
  fMCSecondaryMaterial(source.fMCSecondaryMaterial),
//fHMPIDmomentum(source.fHMPIDmomentum),
//fHMPIDsignal(source.fHMPIDsignal),
  nSigmaPionTPC(source.nSigmaPionTPC),
  nSigmaKaonTPC(source.nSigmaKaonTPC),
  nSigmaProtonTPC(source.nSigmaProtonTPC),
  nSigmaPionTOF(source.nSigmaPionTOF),
  nSigmaKaonTOF(source.nSigmaKaonTOF),
  nSigmaProtonTOF(source.nSigmaProtonTOF),
  fTPCchi2(source.fTPCchi2),
  fITSFakeFlag(source.fITSFakeFlag),
  fTrackCutFlag(source.fTrackCutFlag),
  fEMCalE(source.fEMCalE),
  fEMCalP(source.fEMCalP),
  fHasEMCal(source.fHasEMCal),
  fTimeZeroSigma(source.fTimeZeroSigma)
{
  /*
   * copy constructor
   */

  for (Int_t i = 0; i < 2; i++) fImpactParameter[i] = source.fImpactParameter[i];
  //for (Int_t i = 0; i < 3; i++) fImpactParameterCov[i] = source.fImpactParameterCov[i];
  for (Int_t i = 0; i < 5; i++) fTOFExpTime[i] = source.fTOFExpTime[i];
  for (Int_t i = 0; i < 3; i++) fTOFLabel[i] = source.fTOFLabel[i];

}

//___________________________________________________________

AliAnalysisPIDTrack &
AliAnalysisPIDTrack::operator=(const AliAnalysisPIDTrack &source)
{
  /*
   * operator=
   */

  if (&source == this) return *this;
  TObject::operator=(source);

  fP = source.fP;
  fPt = source.fPt;
  fEta = source.fEta;
  fPhi = source.fPhi;
  fSign = source.fSign;
  fStatus = source.fStatus;
  fLabel = source.fLabel;
  for (Int_t i = 0; i < 2; i++) fImpactParameter[i] = source.fImpactParameter[i];
  //for (Int_t i = 0; i < 3; i++) fImpactParameterCov[i] = source.fImpactParameterCov[i];
  fTPCmomentum = source.fTPCmomentum;
  fTPCdEdx = source.fTPCdEdx;
  fTPCdEdxN = source.fTPCdEdxN;
  fTPCNcls = source.fTPCNcls;
  fTPCNclsF = source.fTPCNclsF;
  fTPCNcr = source.fTPCNcr;
  fTOFIndex = source.fTOFIndex;
  fTOFTime = source.fTOFTime;
  for (Int_t i = 0; i < 5; i++) fTOFExpTime[i] = source.fTOFExpTime[i];
  fTOFLength = source.fTOFLength;
  fTOFDeltaX = source.fTOFDeltaX;
  fTOFDeltaZ = source.fTOFDeltaZ;
  for (Int_t i = 0; i < 3; i++) fTOFLabel[i] = source.fTOFLabel[i];
  fMCPrimary = source.fMCPrimary;
  fMCPdgCode = source.fMCPdgCode;
  fMCMotherPrimary = source.fMCMotherPrimary;
  fMCMotherPdgCode = source.fMCMotherPdgCode;
  fMCMotherLabel = source.fMCMotherLabel;
  fMCTOFMatchPrimary = source.fMCTOFMatchPrimary;
  fMCTOFMatchPdgCode = source.fMCTOFMatchPdgCode;
  fMCTOFMatchLevel = source.fMCTOFMatchLevel;
  fMCTOFTime = source.fMCTOFTime;
  fMCTOFLength = source.fMCTOFLength;
  fMCSecondaryWeak = source.fMCSecondaryWeak;
  fMCSecondaryMaterial = source.fMCSecondaryMaterial;
  //fHMPIDmomentum = source.fHMPIDmomentum;
  //fHMPIDsignal = source.fHMPIDsignal;
  nSigmaPionTPC = source.nSigmaPionTPC;
  nSigmaKaonTPC = source.nSigmaKaonTPC;
  nSigmaProtonTPC = source.nSigmaProtonTPC;
  nSigmaPionTOF = source.nSigmaPionTOF;
  nSigmaKaonTOF = source.nSigmaKaonTOF;
  nSigmaProtonTOF = source.nSigmaProtonTOF;
  fTPCchi2 = source.fTPCchi2;
  fITSFakeFlag = source.fITSFakeFlag;
  fTrackCutFlag = source.fTrackCutFlag;
  fEMCalE = source.fEMCalE;
  fEMCalP = source.fEMCalP;
  fHasEMCal = source.fHasEMCal;
  fTimeZeroSigma = source.fTimeZeroSigma;

  return *this;
}

//___________________________________________________________

AliAnalysisPIDTrack::~AliAnalysisPIDTrack()
{
  /*
   * default destructor
   */

}

//___________________________________________________________

void
AliAnalysisPIDTrack::Reset()
{
  /*
   * reset
   */

  fP = 0.;
  fPt = 0.;
  fEta = 0.;
  fPhi = 0.;
  fSign = 0.;
  fStatus = 0;
  fLabel = 0;
  for (Int_t i = 0; i < 2; i++) fImpactParameter[i] = 0.;
  //for (Int_t i = 0; i < 3; i++) fImpactParameterCov[i] = 0.;
  fTPCmomentum = 0.;
  fTPCdEdx = 0.;
  fTPCdEdxN = 0;
  fTPCNcls = 0;
  fTPCNclsF = 0;
  fTPCNcr = 0.;
  fTOFIndex = 0;
  fTOFTime = 0.;
  for (Int_t i = 0; i < 5; i++) fTOFExpTime[i] = 0.;
  fTOFLength = 0.;
  fTOFDeltaX = 0.;
  fTOFDeltaZ = 0.;
  for (Int_t i = 0; i < 3; i++) fTOFLabel[i] = 0;
  fMCPrimary = kFALSE;
  fMCPdgCode = 0;
  fMCMotherPrimary = kFALSE;
  fMCMotherPdgCode = 0;
  fMCMotherLabel = 0;
  fMCTOFMatchPrimary = kFALSE;
  fMCTOFMatchPdgCode = 0;
  fMCTOFMatchLevel = -1;
  fMCTOFTime = 0.;
  fMCTOFLength = 0.;
  fMCSecondaryWeak = 0.;
  fMCSecondaryMaterial = 0.;
  //fHMPIDmomentum = 0.;
  //fHMPIDsignal = 0.;
  nSigmaPionTPC = 0.;
  nSigmaKaonTPC = 0;
  nSigmaProtonTPC = 0;
  nSigmaPionTOF = 0;
  nSigmaKaonTOF = 0;
  nSigmaProtonTOF = 0;
  fTPCchi2 = 0.;
  fITSFakeFlag = kFALSE;
  fTrackCutFlag = 0;
  fTimeZeroSigma = 0.;
  
}

//___________________________________________________________

void
AliAnalysisPIDTrack::Update(AliESDtrack *track, AliStack *stack, AliMCEvent *mcevent, AliPIDResponse *PIDRes, Int_t TrackCutFlag)
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
  fTPCmomentum = track->GetInnerParam() ? track->GetInnerParam()->P() : 0.;
  fTPCdEdx = track->GetTPCsignal();
  fTPCdEdxN = track->GetTPCsignalN();
  fTPCNcls = track->GetTPCNcls();
  fTPCNclsF = track->GetTPCNclsF();
  fTPCNcr = track->GetTPCClusterInfo(2,1);
  fTPCchi2 = track->GetTPCchi2();
  fITSFakeFlag = track->GetITSFakeFlag();
  fTOFIndex = track->GetTOFCalChannel();
  fTOFTime = track->GetTOFsignal();
  nSigmaPionTPC = PIDRes->NumberOfSigmasTPC(track,AliPID::kPion);
  nSigmaKaonTPC = PIDRes->NumberOfSigmasTPC(track,AliPID::kKaon);
  nSigmaProtonTPC = PIDRes->NumberOfSigmasTPC(track,AliPID::kProton);
  nSigmaPionTOF = PIDRes->NumberOfSigmasTOF(track,AliPID::kPion);
  nSigmaKaonTOF = PIDRes->NumberOfSigmasTOF(track,AliPID::kKaon);
  nSigmaProtonTOF = PIDRes->NumberOfSigmasTOF(track,AliPID::kProton);
  Double_t timei[5];
  track->GetIntegratedTimes(timei);
  for (Int_t i = 0; i < 5; i++) fTOFExpTime[i] = timei[i];
  fTOFLength = track->GetIntegratedLength();
  fTOFDeltaX = track->GetTOFsignalDx();
  fTOFDeltaZ = track->GetTOFsignalDz();
  track->GetTOFLabel(fTOFLabel);
  fTrackCutFlag = TrackCutFlag;
  //fHMPIDmomentum = track->GetOuterHmpParam() ? track->GetOuterHmpParam()->P() : 0.;
  /* HMPID signal with cuts */
  /*Float_t xPc=0., yPc=0., xMip=0., yMip=0., thetaTrk=0., phiTrk=0.;
  Int_t nPhot=0, qMip=0;
  track->GetHMPIDtrk(xPc,yPc,thetaTrk,phiTrk);
  track->GetHMPIDmip(xMip,yMip,qMip,nPhot);
  Float_t dist = TMath::Sqrt((xPc-xMip)*(xPc-xMip) + (yPc-yMip)*(yPc-yMip));    
  if (dist < 0.7 && nPhot < 30 && qMip > 100)
    fHMPIDsignal = track->GetHMPIDsignal();
  else
    fHMPIDsignal = 0.;*/
  /* info from track references */
  fMCTOFTime = 0.;
  fMCTOFLength = 0.;
  if (mcevent && fTOFLabel[0] > 0) {
    TParticle *particle;
    TClonesArray *arrayTR;
    AliTrackReference *trackRef;
    mcevent->GetParticleAndTR(fTOFLabel[0], particle, arrayTR);
    if(!arrayTR) {
      AliFatal("AddAnalysisTaskTPCTOFPID needs AliTrackReference!");
      return;
    };
    for (Int_t itr = 0; itr < arrayTR->GetEntries(); itr++) {
      trackRef = (AliTrackReference *)arrayTR->At(itr);
      if (!trackRef || trackRef->DetectorId() != AliTrackReference::kTOF) continue;
      fMCTOFTime = trackRef->GetTime() * 1.e12; /* s -> ps */
      fMCTOFLength = trackRef->GetLength();
      /* break as soon as we get it */
      break;
    }
  }
  /* info from stack */
  fMCPrimary = kFALSE;
  fMCPdgCode = 0;
  fMCMotherPrimary = kFALSE;
  fMCMotherPdgCode = 0;
  fMCSecondaryWeak = kFALSE;
  fMCSecondaryMaterial = kFALSE;
  if (stack) {
    Int_t index = TMath::Abs(fLabel);
    if (index < 0) {
      printf("index = %d\n", index);
      return;
    }
    TParticle *particle = stack->Particle(index);
    fMCPrimary = stack->IsPhysicalPrimary(index);
    fMCSecondaryWeak = stack->IsSecondaryFromWeakDecay(index);
    fMCSecondaryMaterial = stack->IsSecondaryFromMaterial(index);
    fMCPdgCode = particle->GetPdgCode();
    Int_t indexm = particle->GetFirstMother();
    if (indexm < 0) {
      fMCMotherPrimary = kFALSE;
      fMCMotherPdgCode = 0;
    }
    else {
      TParticle *particlem = stack->Particle(indexm);
      fMCMotherPrimary = stack->IsPhysicalPrimary(indexm);
      fMCMotherPdgCode = particlem->GetPdgCode();
      fMCMotherLabel = indexm;
    }

    /* check TOF match */
    fMCTOFMatchPrimary = kFALSE;
    fMCTOFMatchPdgCode = 0;
    fMCTOFMatchLevel = -1;
    if (fTOFLabel[0] > 0) {
      index = fTOFLabel[0];
      particle = stack->Particle(index);
      fMCTOFMatchPrimary = stack->IsPhysicalPrimary(index);
      fMCTOFMatchPdgCode = particle->GetPdgCode();
      Int_t tracklabel = TMath::Abs(fLabel);
      Int_t matchlevel = -1;
      for (Int_t ilevel = 0; index > 0; ilevel++) {
	if (index == tracklabel) {
	  matchlevel = ilevel;
	  break;
	}
	index = stack->Particle(index)->GetFirstMother();
      }
      fMCTOFMatchLevel = matchlevel;
    }
    
  }
  fTimeZeroSigma = 0.;
  
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::HasTPCPID() const
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
AliAnalysisPIDTrack::HasTOFPID(TH1 *henabled) const
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

Bool_t
AliAnalysisPIDTrack::IsTPCPID(Int_t ipart, Float_t cutTPC) const
{
  /*
   * is TPC PID
   */

  Float_t nSigma[AliPID::kSPECIES];
  if (!MakeTPCPID(nSigma)) return kFALSE;
  if (TMath::Abs(nSigma[ipart]) > cutTPC) return kFALSE;
  return kTRUE;
}
  
//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::MakeTPCPID(Float_t *nSigma) const
{
  /*
   * make TPC PID
   */
  
  /* check TPC PID */
  if (!HasTPCPID()) return kFALSE;
  /* loop over particles */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    Double_t bethe = fgTPCResponse->GetExpectedSignal(fTPCmomentum, (AliPID::EParticleType)ipart);
    Double_t diff = fTPCdEdx - bethe;
    Double_t sigma = fgTPCResponse->GetExpectedSigma(fTPCmomentum, fTPCdEdxN, (AliPID::EParticleType)ipart);
    nSigma[ipart] = diff / sigma;
  }

  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::MakeTOFPID(Float_t *nSigma) const
{
  /*
   * make TOF PID
   */
  
  /* check TOF PID */
  if (!HasTOFPID()) return kFALSE;

  /* loop over particles */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    Double_t timez = fTOFTime - fTOFExpTime[ipart] - fgTOFResponse->GetStartTime(fP);
    Double_t sigma = fgTOFResponse->GetExpectedSigma(fP, fTOFExpTime[ipart], AliPID::ParticleMass(ipart));
    nSigma[ipart] = timez / sigma;
  }

  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::IsMismatch(const Float_t *nSigmaTPC, const Float_t *nSigmaTOF, Float_t cutTPC, Float_t cutTOF) const
{
  /*
   * is mismatch
   */

  /* search for valid and compatible PID */
  Bool_t valid = kFALSE;
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    if (TMath::Abs(nSigmaTPC[ipart]) < cutTPC) {
      valid = kTRUE;
      if (TMath::Abs(nSigmaTOF[ipart]) < cutTOF) return kFALSE;
    }
  }
  /* check valid TPC PID */
  if (!valid) return kFALSE;
  /* check beta-gamma compatible */
  if (IsBetaGammaCompatible(cutTPC, cutTOF)) return kFALSE;
  /* mismatch */
  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::IsMismatch(Float_t cutTPC, Float_t cutTOF) const
{
  /*
   * is mismatch
   */

  Float_t nSigmaTPC[AliPID::kSPECIES];
  Float_t nSigmaTOF[AliPID::kSPECIES];
  Bool_t hasTPCPID = MakeTPCPID(nSigmaTPC);
  Bool_t hasTOFPID = MakeTOFPID(nSigmaTOF);
  /* check both TPC and TOF PID */
  if (!hasTPCPID || !hasTOFPID) return kFALSE;
  return IsMismatch(nSigmaTPC, nSigmaTOF, cutTPC, cutTOF);
}

//___________________________________________________________

void
AliAnalysisPIDTrack::RemoveTimeZero(const AliAnalysisPIDEvent *analysisEvent)
{
  /*
   * remove time-zero
   */

  if (!analysisEvent) return;
  fTOFTime -= analysisEvent->GetTimeZeroSafe(fP);
  fTimeZeroSigma = analysisEvent->GetTimeZeroSafeSigma(fP);
  //  fTOFTime -= analysisEvent->GetTimeZeroBest(fP);
  //  fTimeZeroSigma = analysisEvent->GetTimeZeroBestSigma(fP);
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::IsTPCDeuton() const
{
  /*
   * is TPC deuton
   */

  /* check TPC PID */
  if (!HasTPCPID()) return kFALSE;
  Double_t mass = 1.8756; /* GeV */
  Double_t bethe = fgTPCResponse->Bethe(fTPCmomentum / mass);
  Double_t diff = fTPCdEdx - bethe;
  Double_t sigma = bethe * 0.07;
  if (TMath::Abs(diff / sigma) < 5.) return kTRUE;
  return kFALSE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::IsTPCTriton() const
{
  /*
   * is TPC triton
   */

  /* check TPC PID */
  if (!HasTPCPID()) return kFALSE;
  Double_t mass = 2. * 0.939565560 + AliPID::ParticleMass(4); /* GeV */
  Double_t bethe = fgTPCResponse->Bethe(fTPCmomentum / mass);
  Double_t diff = fTPCdEdx - bethe;
  Double_t sigma = bethe * 0.07;
  if (TMath::Abs(diff / sigma) < 5.) return kTRUE;
  return kFALSE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::IsTPCHeavy() const
{
  /*
   * is TPC heavy
   */
  
  /* check TPC heavier than proton */
  if (!HasTPCPID()) return kFALSE;
  Double_t mass = 1.8756; /* GeV */
  Double_t bethe = fgTPCResponse->Bethe(fTPCmomentum / mass);
  Double_t diff = fTPCdEdx - bethe;
  Double_t sigma = bethe * 0.07;
  if (diff / sigma > -5.) return kTRUE;
  return kFALSE;
}

//___________________________________________________________

Float_t
AliAnalysisPIDTrack::GetTOFBetaSigma() const
{
  /*
   * get TOF beta sigma
   */

  Float_t tofReso = TMath::Sqrt(80. * 80. + fTimeZeroSigma * fTimeZeroSigma);
  return GetTOFBeta() * tofReso / fTOFTime;
}

//___________________________________________________________

Float_t
AliAnalysisPIDTrack::GetTPCdEdxTh(Float_t betagamma) const
{
  /*
   * get TPD dEdx th
   */

  return fgTPCResponse->Bethe(betagamma);
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::IsBetaGammaCompatible(Float_t cutTPC, Float_t cutTOF) const
{
  /*
   * is beta-gamma compatible
   */

  Double_t epsilon = 1.e-10;
  /* check TPC PID */
  if (!HasTPCPID()) return kTRUE;
  /* get beta and sigma */
  Double_t beta = GetTOFBeta();
  Double_t betasigma = GetTOFBetaSigma();
  /* get and check beta min/max */
  Double_t betamin = beta - cutTOF * betasigma;
  Double_t betamax = beta + cutTOF * betasigma;
  if (betamin <= 0.) betamin = epsilon;
  if (betamin >= 1.) betamin = 1. - epsilon;
  if (betamax <= 0.) betamax = epsilon;
  if (betamax >= 1.) betamax = 1. - epsilon;
  /* check TOF beta-gamma min/max compatible with TPC dEdx */
  Double_t betagammamin = betamin / TMath::Sqrt(1. - betamin * betamin);
  Double_t betagammamax = betamax / TMath::Sqrt(1. - betamax * betamax);
  Double_t bethemin = fgTPCResponse->Bethe(betagammamin);
  Double_t bethemax = fgTPCResponse->Bethe(betagammamax);
  Double_t diffmin = fTPCdEdx - bethemin;
  Double_t diffmax = fTPCdEdx - bethemax;
  Double_t sigmamin = bethemin * 0.07;
  Double_t sigmamax = bethemax * 0.07;
  Double_t nsigmamin = diffmin / sigmamin;
  Double_t nsigmamax = diffmax / sigmamax;
  if (TMath::Abs(nsigmamin) < cutTPC) return kTRUE;
  if (TMath::Abs(nsigmamax) < cutTPC) return kTRUE;
  return kFALSE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::IsHeavyAndCompatible(const Float_t *nSigmaTOF, Float_t cutTPC, Float_t cutTOF) const
{
  /*
   * is heavy and compatible
   */

  /* check TOF time larger than proton expected time */
  if (nSigmaTOF[AliPID::kProton] < cutTOF) return kFALSE;
  return IsBetaGammaCompatible(cutTPC, cutTOF);
}

//___________________________________________________________


//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::LoadTuningExpTimeTh(const Char_t *filename)
{
  /*
   * load tuning exp time th
   */

  TFile *file = TFile::Open(filename);
  if (!file || !file->IsOpen()) return kFALSE;
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
    hTOFtuned_th[ipart] = (TH2F *)file->Get(Form("hTOFtuned_th_%s", AliPID::ParticleName(ipart)));
  return kTRUE;
  
}

//___________________________________________________________

Float_t
AliAnalysisPIDTrack::GetTuningExpTimeTh(Int_t ipart) const
{
  /*
   * get tuning exp time th
   */

  if (!hTOFtuned_th[ipart]) return 0.;
  if (fP > 1.5) return 0.;
  Int_t bin = hTOFtuned_th[ipart]->FindBin(fP, fEta);
  Float_t value = hTOFtuned_th[ipart]->GetBinContent(bin);
  return value;
}

//___________________________________________________________

/*Int_t
AliAnalysisPIDTrack::GetTOFVolumeIndex(Int_t i)
{


  if (i < 0 || i > 4) return -1;
  Int_t det[5];
  fgTOFGeometry.GetVolumeIndices(fTOFIndex, det);
  return det[i];

}*/

//___________________________________________________________

Int_t
AliAnalysisPIDTrack::GetMCPID() const
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
AliAnalysisPIDTrack::GetMCCharge() const
{
  /*
   * get MC charge
   */
  
  TDatabasePDG *dbpdg = TDatabasePDG::Instance();
  TParticlePDG *pdg = dbpdg->GetParticle(fMCPdgCode);
  if (!pdg) return 0;
  return (Int_t)TMath::Sign(1., pdg->Charge());
}

//___________________________________________________________

Float_t
AliAnalysisPIDTrack::GetTOFExpTimeSigma(Int_t i) const
{
  /*
   * get TOF integrated times sigma
   */

  Double_t par[5][3] = {
    {1.84217e-04, 2.50350e-04, -2.06718e+00},
    {1.84217e-04, 2.50350e-04, -2.06718e+00},
    {1.84217e-04, 2.50350e-04, -2.06718e+00},
    {1.30006e-04, 1.65197e-03, -1.21072e+00},
    {1.89219e-03, 2.06852e-03, -1.85165e+00}
  };
  
  return 1.5 * (par[i][0] + par[i][1] * TMath::Power(fP, par[i][2])) * fTOFExpTime[i];

}

//___________________________________________________________

Float_t
AliAnalysisPIDTrack::GetTOFExpTimeCorrection(Int_t i, Int_t chargeCorr) const
{
  /*
   * get TOF integrated times correction
   */

  /* for both charges */
  Double_t par[AliPID::kSPECIES][5] = {
    {-40., 0., 0., 0., 1.},
    {-4.26411e+00, -8.17232e+01, 2.84390e+02, -1.10976e-02, -3.76187e+00},
    {-4.26411e+00, -8.17232e+01, 2.84390e+02, -1.10976e-02, -3.76187e+00},
    {3.74297e+00, -4.90544e+03, 0.00000e+00, 0.00000e+00, -8.79427e+00},
    {3.36701e+00, 1.60189e+04, -3.43919e+04, 8.27698e-02, -6.49010e+00}
  };
  /* split charges */
  Double_t parpos[AliPID::kSPECIES][2] = {
    {0., 0.},
    {1.97815e+00, -1.76842e-01},
    {1.97815e+00, -1.76842e-01},
    {1.66051e+01, -8.98377e-01},
    {2.09552e+01, -3.18122e-01}
  };
  Double_t parneg[AliPID::kSPECIES][2] = {
    {0., 0.},
    {-5.91210e+00, -1.04768e+00},
    {-5.91210e+00, -1.04768e+00},
    {-1.58687e+01, -7.93304e-01},
    {-2.90316e+01, -4.05392e-01}
  };
  
  Double_t x = fP;
  Double_t corr = par[i][0]+(par[i][1]*x+par[i][2]*x*x+par[i][3]*x*x*x)*TMath::Exp(par[i][4]*x);
  Double_t chcorr;
  if (fSign > 0.)
    chcorr = parpos[i][0]*TMath::Exp(parpos[i][1]*x);
  else if (fSign < 0.)
    chcorr = parneg[i][0]*TMath::Exp(parneg[i][1]*x);
  else
    return 0.;
  
  return corr + chargeCorr * chcorr;
}

//___________________________________________________________

void
AliAnalysisPIDTrack::UpdateTOFResponse(AliAnalysisPIDEvent *analysisEvent)
{
  /*
   * update TOF response
   */

  fgTOFResponse->SetT0event(analysisEvent->GetTimeZeroTOF());
  fgTOFResponse->SetT0resolution(analysisEvent->GetTimeZeroTOFSigma());
}

//___________________________________________________________

Float_t
AliAnalysisPIDTrack::GetTOFExpectedSigma(Int_t i) const
{
  /*
   * get TOF expected sigma
   */

  return fgTOFResponse->GetExpectedSigma(fP, fTOFExpTime[i], AliPID::ParticleMass(i));

}

//___________________________________________________________

void
AliAnalysisPIDTrack::ApplyTOFExpectedTimeCorrection(Int_t chargeCorr)
{
  /*
   * apply TOF expected time correction
   */

  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
    fTOFExpTime[ipart] += GetTOFExpTimeCorrection(ipart, chargeCorr);
  
}

//___________________________________________________________

Bool_t
AliAnalysisPIDTrack::AcceptTrack(Bool_t selPrimaries)
{
  /*
   * accept track
   */

  /* check charge */
  if (fSign == 0.) return kFALSE;
  /* check eta */
  if (TMath::Abs(fEta) > fgEtaCut) return kFALSE;
  if (TMath::Abs(fEta) < fgEtaReject) return kFALSE;
  /* check max DCA to vertex Z */
  if (TMath::Abs(fImpactParameter[1]) > fgMaxDCAToVertexZCut) return kFALSE;
  /* check max DCA to vertex XY */
  if (selPrimaries) {
    Float_t maxdca = fgMaxDCAToVertexXYPtDepFormula->Eval(fPt);
    if (TMath::Abs(fImpactParameter[0]) > maxdca) return kFALSE;
  }
  /* check max chi2 per cluster TPC */
  Float_t chi2 = fTPCchi2 / (Float_t)fTPCNcls;
  if (chi2 > fgMaxChi2PerClusterTPC) return kFALSE;
  /* check min N clusters TPC */
  if (fgAcceptTrackClusterCut != 1) {
    if (fTPCNcls < fgMinNClustersTPC) return kFALSE;
  }
  if (fgAcceptTrackClusterCut == 1) {
    /* check min N crossed rows TPC */
    if (fTPCNcr < fgMinNCrossedRowsTPC) return kFALSE;
    /* check min crossed rows over findable clusters ratio */
    Float_t crratio = 1.;
    if (fTPCNclsF > 0) crratio = fTPCNcr / fTPCNclsF;
    if (crratio < fgMinRatioCrossedRowsOverFindableClustersTPC) return kFALSE;
  }
  /* check accept track status cut */
  if (fgAcceptTrackStatusCut && (fStatus & fgAcceptTrackStatusCut) == 0) return kFALSE;
  /* check reject track status cut */
  if (fgRejectTrackStatusCut && (fStatus & fgRejectTrackStatusCut) != 0) return kFALSE;
  /* reject ITS fakes if requested */
  if (fgRejectITSFakes && fITSFakeFlag) return kFALSE;

  /* accept track */
  return kTRUE;
}
