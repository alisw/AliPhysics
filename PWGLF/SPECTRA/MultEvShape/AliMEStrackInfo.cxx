//////////////////////////////////////////////////////////////////////////////
//
// A much longer description of track info
//
//////////////////////////////////////////////////////////////////////////////

#include <AliLog.h>
#include <AliESDtrack.h>
#include <AliPIDResponse.h>
#include <AliPIDCombined.h>
#include <AliStack.h>
#include <AliMCParticle.h>

#include "AliMESbaseTask.h"
#include "AliMEStrackInfo.h"

ClassImp(AliMEStrackInfo)
ClassImp(AliMEStrackInfo::AliMESpid)
ClassImp(AliMEStrackInfo::AliMESfilterParam)

//______________________________________________________________
AliMEStrackInfo::AliMEStrackInfo()
  : AliVParticle()
  ,fTrackId(-1)
  ,fOrigin(0)
  ,fTOFbc(0)
  ,fFilterId(0)
  ,fDetStat(0)
  ,fPt(0.)
  ,fP(0.)
  ,fPz(0.)
  ,fEta(-100.)
  ,fPhi(-100.)
  ,fY(0.)
  ,fdEdx(0.)
  ,fBeta(0.)
  ,fPID()
  ,fFilterParam(NULL)
{
  //
  // Constructor
  //
  memset(fPosition, 0, 3*sizeof(Double_t));
  memset(fDCA, 0, 2*sizeof(Double_t));
}
//______________________________________________________________
AliMEStrackInfo::AliMEStrackInfo(const AliMEStrackInfo &t)
: AliVParticle((const AliVParticle&)t)
,fTrackId(t.fTrackId)
,fOrigin(t.fOrigin)
,fTOFbc(t.fTOFbc)
,fFilterId(t.fFilterId)
,fDetStat(t.fDetStat)
,fPt(t.fPt)
,fP(t.fP)
,fPz(t.fPz)
,fEta(t.fEta)
,fPhi(t.fPhi)
,fY(t.fY)
,fdEdx(t.fdEdx)
,fBeta(t.fBeta)
,fPID()
,fFilterParam(NULL)
{
	//
	// Copy Constructor
	//
	memcpy(fPID.fNsigma, t.fPID.fNsigma, kNdet*AliPID::kSPECIES*sizeof(Int_t));
	fPID.fTOFmisProb = t.fPID.fTOFmisProb;
	memcpy(fPID.fProb, t.fPID.fProb, kNdet*AliPID::kSPECIES*sizeof(Double_t));
	memcpy(fPID.fRaw, t.fPID.fRaw, kNdet*sizeof(Double_t));

	if(t.fFilterParam){
		fFilterParam = new AliMESfilterParam();
		fFilterParam->fNcl = t.fFilterParam->fNcl;
		fFilterParam->fChi2Cl = t.fFilterParam->fChi2Cl;
	}
	memcpy(fDCA, t.fDCA, 3*sizeof(Double_t));
}

//______________________________________________________________
AliMEStrackInfo::AliMEStrackInfo(AliESDtrack *t, AliPIDResponse *rpid, AliPIDCombined *pidComb)
  : AliVParticle()
  ,fTrackId(TMath::Abs(t->GetLabel()))
  ,fOrigin(0)
  // ,fTOFbc(t->GetTOFBunchCrossing())
  ,fTOFbc(0)
  ,fFilterId(0)
  ,fDetStat(0)
  ,fPt(t->GetSignedPt())
//   ,fPt(t->Charge()*t->Pt())
  ,fP(t->P())
  ,fPz(t->Pz())
  ,fEta(t->Eta())
  ,fPhi(t->Phi())
  ,fY(0.)
  ,fdEdx(t->GetTPCsignal())
  ,fBeta(-1.)
  ,fPID()
  ,fFilterParam(NULL)
{
  //
  // Constructor from reconstructed track
  // to be checked and further implemented
  //
  if(fTOFbc==AliVTrack::kTOFBCNA) fTOFbc=0; // reset TOF bc
  // set status
  if(t->GetStatus()&AliESDtrack::kTRDin) SetDetStat(kTRD);
  // fill Position
  t->GetXYZ(fPosition);
  // combined PID
  Double_t bayesProb[AliPID::kSPECIES];

  // set PID ITS
  pidComb->SetDetectorMask(AliPIDResponse::kDetITS|AliPIDResponse::kDetTPC);
  pidComb->ComputeProbabilities(t, rpid, bayesProb);
  for(Int_t is(0); is<AliPID::kSPECIES; is++) fPID.fNsigma[kITS][is] = rpid->NumberOfSigmasITS(t, AliPID::EParticleType(is));
  memcpy(fPID.fProb[kITS], bayesProb, AliPID::kSPECIES*sizeof(Double_t));
  fPID.fRaw[kITS] = t->GetITSsignal();

  // set PID TPC
  pidComb->SetDetectorMask(AliPIDResponse::kDetTPC);
  pidComb->ComputeProbabilities(t, rpid, bayesProb);
  for(Int_t is(0); is<AliPID::kSPECIES; is++) fPID.fNsigma[kITS][is] = rpid->NumberOfSigmasTPC(t, AliPID::EParticleType(is));
  memcpy(fPID.fProb[kTPC], bayesProb, AliPID::kSPECIES*sizeof(Double_t));
  fPID.fRaw[kTPC] = t->GetTPCsignal();

  // set PID TOF
  // pidComb->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);
  pidComb->SetDetectorMask(AliPIDResponse::kDetITS|AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);
  pidComb->ComputeProbabilities(t, rpid, bayesProb);
  for(Int_t is(0); is<AliPID::kSPECIES; is++) fPID.fNsigma[kITS][is] = rpid->NumberOfSigmasTOF(t, AliPID::EParticleType(is));
  memcpy(fPID.fProb[kTOF], bayesProb, AliPID::kSPECIES*sizeof(Double_t));
  fPID.fRaw[kTOF] = (t->GetIntegratedLength()/(t->GetTOFsignal()*TMath::C()))*10e9;

  // set origin // TODO
  SetOrigin(kPrimary);
  // -> analyze V0

  // set DCAxy and DCAz
//   Float_t dca[2];  // 0 = xy; 1 = z
  Float_t bCov[3];
  t->GetImpactParameters(fDCA, bCov);

  // set the PID QA info

  // beta
  Double_t tof = t->GetTOFsignal();
  if(tof > 0.) fBeta = ( (t->GetIntegratedLength()) / (tof * TMath::C()) ) *10e9;
  else fBeta = -1.;

  // fill in debug mode also the filter param
  if(AliMESbaseTask::DebugUsers()){
    fFilterParam = new AliMESfilterParam;
    fFilterParam->fNcl = t->GetNcls(1);
    fFilterParam->fChi2Cl = t->GetTPCchi2();
  }
}

//______________________________________________________________
AliMEStrackInfo::AliMEStrackInfo(AliMCParticle *t, AliStack *mc)
  : AliVParticle()
  // ,fTrackId(-1)
  ,fTrackId(t->GetLabel())
  ,fOrigin(0)
  ,fTOFbc(0.)
  ,fFilterId(0)
  ,fDetStat(0)
  ,fPt((t->Charge()>0?1:-1)*t->Pt())
  ,fP(t->P())
  ,fPz(t->Pz())
  ,fEta(t->Eta())
  ,fPhi(t->Phi())
  ,fY(t->Y())
  ,fdEdx(-1.)
  ,fBeta(-1.)
  ,fPID()
  ,fFilterParam(NULL)
{
  //
  // Constructor from reconstructed track
  // to be checked and further implemented

  memset(fPosition, 0, 3*sizeof(Double_t));
  memset(fDCA, 0, 2*sizeof(Double_t));

  Int_t pdgCode = t->PdgCode();
  memset(fPID.fProb, 0, kNdet*AliPID::kSPECIES*sizeof(Double_t));
  switch(TMath::Abs(pdgCode)){
  case 11:   fPID.fProb[kITS][0] = 1.; break;
  case 13:   fPID.fProb[kITS][1] = 1.; break;
  case 211:  fPID.fProb[kITS][2] = 1.; break;
  case 321:  fPID.fProb[kITS][3] = 1.; break;
  case 2212: fPID.fProb[kITS][4] = 1.; break;
  }

  // set origin
  TParticle* particle = t->Particle(); //mc->Particle(TMath::Abs(label));
  if (mc->IsPhysicalPrimary(TMath::Abs(t->Label()))) SetOrigin(kPrimary);
  else{
    Int_t uniqueID=particle->GetUniqueID();
    if(uniqueID==kPDecay) SetOrigin(kSecondary);
    else if(uniqueID==kPHadronic) SetOrigin(kMaterial);
  }

//   AliInfo(Form("MC charge: %i\n", TMath::Sign(1,(t->Charge()))));
}

//______________________________________________________________
AliMEStrackInfo::~AliMEStrackInfo()
{
  // Destructor
  if(fFilterParam) delete fFilterParam;
}

//______________________________________________________________
Bool_t AliMEStrackInfo::PxPyPz(Double_t*) const
{
  // implemented from AliVParticle
  AliInfo("Not used. Use Px() Py() Pz() instead.");
  return kTRUE;
}


//______________________________________________________________
Double_t AliMEStrackInfo::M() const
{
  // implemented from AliVParticle
  AliInfo("TODO - to be implemented.");
  return 0.;
}

//______________________________________________________________
Int_t AliMEStrackInfo::PdgCode() const
{
  // implemented from AliVParticle
  AliInfo("TODO - to be implemented.");
  return 0;
}

//______________________________________________________________
const Double_t* AliMEStrackInfo::PID() const
{
  // implemented from AliVParticle
  AliInfo("Not used. Use GetPID() instead.");
  return NULL;
}

//______________________________________________________________
AliMEStrackInfo::AliMESpid::AliMESpid()
  : TObject()
  , fTOFmisProb(0.)
{
  //
  // Constructor
  //
  memset(fNsigma, 0, kNdet*AliPID::kSPECIES*sizeof(Int_t));
  memset(fProb, 0, kNdet*AliPID::kSPECIES*sizeof(Double_t));
  memset(fRaw, 0, kNdet*sizeof(Double_t));
}

//______________________________________________________________
AliMEStrackInfo::AliMESfilterParam::AliMESfilterParam()
  : TObject()
  ,fNcl(0)
  ,fChi2Cl(0.)
{
  //
  // Constructor
  //
}
