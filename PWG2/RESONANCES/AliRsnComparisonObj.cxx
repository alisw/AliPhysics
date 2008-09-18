//
// *** Class AliRsnComparisonObj ***
//
//  TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnMCInfo.h"

#include "AliRsnComparisonObj.h"

ClassImp(AliRsnComparisonObj)

//_____________________________________________________________________________
AliRsnComparisonObj::AliRsnComparisonObj(const char *name)
    : TNamed(name,name),fCurrentComparisonType(kParticleInfo),fCurrentESDPID(kEsd),
    fESDstatus(0),fITSClusters(0),fTPCClusters(0),fTRDClusters(0),fPIDDivValue(0.)
{
//=========================================================
// Default constructor
//=========================================================
  TH1::SetDefaultSumw2();

  for (Int_t i=0;i<kLastFormat;i++)
    for (Int_t j=0;j<kLastHistoType;j++)
      for (Int_t k=0;k<AliRsnPID::kSpecies+1;k++)
        fHistosPID[i][j][k] = 0;

  for (Int_t i=0;i<kLastParameterType;i++)
    for (Int_t j=0;j<2;j++)
      for (Int_t k=0;k<AliRsnPID::kSpeciesAll;k++)
        fHistosPartInfo[i][j][k] = 0;

  fPriorProbs[0]=0.02;
  fPriorProbs[1]=0.02;
  fPriorProbs[2]=0.83;
  fPriorProbs[3]=0.07;
  fPriorProbs[4]=0.06;
  fITSClusters = -1;
  fTPCClusters = -1;
  fTRDClusters = -1;
  fPIDDivValue = 0.0;

}
//_____________________________________________________________________________
AliRsnComparisonObj::~AliRsnComparisonObj()
{
//=========================================================
// Destructor
//=========================================================
}

//________________________________________________________________________
TList * AliRsnComparisonObj::GeneratePIDHistogramList(TString prefix)
{
//=========================================================
// Generate histograms
//=========================================================

  TList *list = new TList();
  list->SetName(GetName());

  TString varname("p");
  for (Int_t i=0;i<kLastFormat;i++)
    for (Int_t j=0;j<kLastHistoType;j++)
      for (Int_t k=0;k<AliRsnPID::kSpecies+1;k++)
      {
        if ((i == kMC) && (j>kIndent)) continue;

        fHistosPID[i][j][k]= new TH1D(Form("%s_%s_%s_%s_%s_%s",prefix.Data(),GetName(),varname.Data(),AliRsnPID::ParticleName((AliRsnPID::EType) k),GetFormatName((EFormat) i).Data(),GetHistoTypeName((EHistoType) j).Data()),Form("%s %s %s %s %s",GetName(),varname.Data(),AliRsnPID::ParticleNameLatex((AliRsnPID::EType) k),GetFormatName((EFormat) i).Data(),GetHistoTypeName((EHistoType) j).Data()),1000,0,5);
        list->Add(fHistosPID[i][j][k]);
      }

  return list;
}

TList * AliRsnComparisonObj::GenerateParticleInfoHistogramList(TString prefix)
{
  TList *list = new TList();
  list->SetName(GetName());
  TString charge;
  if (!prefix.IsNull())
    prefix+="_";

  Int_t bins[] = {1000,1000,1000,1000};
  Double_t min[] = {0.0,0.0,-10.0,-10.0};
  Double_t max[] = {5.0,5.0,10.0,10.0};

  for (Int_t i=0;i<kLastParameterType;i++)
    for (Int_t j=0;j<2;j++)
      for (Int_t k=0;k<AliRsnPID::kSpeciesAll;k++)
      {
//         if (k == AliRsnPID::kUnknown) continue;
//         if ((j>0) && (k>AliRsnPID::kSpecies-1)) continue;
        charge = (j>0) ? "-":"+";
        AliInfo(Form("%s%s %d%d",AliRsnPID::ParticleName((AliRsnPID::EType) k),charge.Data(),j,k));
//         if (k>AliRsnPID::kSpecies-1) charge="";
        fHistosPartInfo[i][j][k]= new TH1D(Form("%s%s%s_%s",prefix.Data(),AliRsnPID::ParticleName((AliRsnPID::EType) k),charge.Data(),GetParameterName((EParameterType) i).Data()),Form("%s%s %s",AliRsnPID::ParticleNameLatex((AliRsnPID::EType) k),charge.Data(),GetParameterName((EParameterType) i).Data()),bins[i],min[i],max[i]);
        list->Add(fHistosPartInfo[i][j][k]);
      }
  return list;
}

//________________________________________________________________________
void AliRsnComparisonObj::FillPIDHistograms(AliRsnDaughter * daughter)
{
//=========================================================
// Fill Histograms with AliRsnDaughter info
//=========================================================

  if (!(daughter->CheckFlag(fESDstatus))) return;

  daughter->AssignRealisticPID();

  Double_t prob;
  Double_t p = daughter->P();
  Int_t indexType = (Int_t) daughter->PIDType(prob);
  Int_t pdg = daughter->GetMCInfo()->PDG();
  Int_t indexTypeFromPDG = (Int_t) AliRsnPID::InternalType(pdg);

  fHistosPID[kRSN][kIndent][indexType]->Fill(p);
  fHistosPID[kRSN][kTrue][indexTypeFromPDG]->Fill(p);

  if (indexType == indexTypeFromPDG)
    fHistosPID[kRSN][kGood][indexType]->Fill(p);
  else
    fHistosPID[kRSN][kFake][indexType]->Fill(p);
}

//________________________________________________________________________
void AliRsnComparisonObj::FillPIDHistograms(AliESDtrack * track,AliMCEvent *mc)
{
//=========================================================
// Fill Histograms with AliRsnDaughter info
//=========================================================

  if (fITSClusters>0)
    if (track->GetITSclusters(0) < fITSClusters) return;

  if (fTPCClusters>0)
    if (track->GetTPCclusters(0) < fTPCClusters) return;

  if (fTRDClusters>0)
    if (track->GetTRDclusters(0) < fTRDClusters) return;

  ULong_t status = track->GetStatus();
  if (!((fESDstatus & status) == fESDstatus)) return;

  Int_t label = track->GetLabel();
  AliMCParticle *mctrack = mc->GetTrack(TMath::Abs(label));

  Double_t p = track->P();
  Int_t pdg =0;
  Int_t indexTypeFromPDG = AliRsnPID::kUnknown;
  if (mctrack)
  {
    pdg = mctrack->Particle()->GetPdgCode();
    indexTypeFromPDG = (Int_t) AliRsnPID::InternalType(pdg);
  }
  Double_t pid[5];
  GetESDPID(track,pid,p);


  if (fCurrentESDPID == kITS_TPC_TOF_SP)
    if (p>0.7)
      if (!track->IsOn(AliESDtrack::kTOFpid)) return;

  Int_t i;
  Double_t sum = 0.0, prob[AliRsnPID::kSpecies];
  for (i = 0; i < AliRsnPID::kSpecies; i++)
  {
    prob[i] = fPriorProbs[i] * pid[i];
    sum += prob[i];
  }

  if (sum <= (Double_t) 0.0)
  {
    AliError(Form("Sum of weights = %f <= 0", sum));
    return;
  }


  // normalize
  for (i = 0; i < AliRsnPID::kSpecies; i++)
  {
    prob[i] /= sum;
  }

  Int_t indexType = 0;
  Double_t pmax = prob[0];
  for (i = 1; i < AliRsnPID::kSpecies; i++)
  {
    if (prob[i] > pmax)
    {
      indexType = i;
      pmax = prob[i];
    }
  }

  fHistosPID[kESD][kIndent][indexType]->Fill(p);
  fHistosPID[kESD][kTrue][indexTypeFromPDG]->Fill(p);

  if (indexType == indexTypeFromPDG)
    fHistosPID[kESD][kGood][indexType]->Fill(p);
  else
    fHistosPID[kESD][kFake][indexType]->Fill(p);

}

//________________________________________________________________________
void AliRsnComparisonObj::FillPIDHistograms(AliMCParticle * mctrack)
{
//=========================================================
// Fill Histograms with mctrack info
//=========================================================

  if (!mctrack) { AliError("mctrack is null"); return;}

  Double_t p = mctrack->P();
  Int_t pdg = mctrack->Particle()->GetPdgCode();
  Int_t indexTypeFromPDG = (Int_t) AliRsnPID::InternalType(pdg);
  if (indexTypeFromPDG < AliRsnPID::kSpecies)
    fHistosPID[kMC][kIndent][indexTypeFromPDG]->Fill(p);
}

void AliRsnComparisonObj::FillHistograms(AliMCParticle * mctrack)
{
  Int_t charge;
  if (!mctrack) { AliError("mctrack is null"); return;}

  Int_t pdg = mctrack->Particle()->GetPdgCode();
//   if (pdg>0) charge=0;
  charge = (pdg>0) ? 0 : 1;
//   if (TMath::Abs(pdg) == 333  ) AliInfo(Form("%d",charge));
  if (pdg == 11) charge = 1;
  if (pdg == -11) charge = 0;

  Int_t indexTypeFromPDG = (Int_t) AliRsnPID::InternalType(pdg);
  for (Int_t i=0;i<kLastParameterType;i++)
    fHistosPartInfo[i][charge][indexTypeFromPDG]->Fill(GetParameterNameValue((EParameterType)i,mctrack));
}

//________________________________________________________________________
TString AliRsnComparisonObj::GetFormatName(EFormat type)
{
//=========================================================
// returns Format Name
//=========================================================
  switch (type)
  {
    case kRSN :
      return "RSN";
      break;
    case kESD :
      return "ESD";
      break;
    case kMC :
      return "MC";
      break;
    default:
      AliWarning("Unrecognized value of EOutputType argument");
      break;
  }
  return "";
}

//________________________________________________________________________
TString AliRsnComparisonObj::GetHistoTypeName(EHistoType type)
{
//=========================================================
// returns Histo Type Name
//=========================================================
  switch (type)
  {
    case kIndent :
      return "Ident";
    case kGood :
      return "Good";
    case kFake :
      return "Fake";
    case kTrue :
      return "True";
    default:
      AliWarning("Unrecognized value of EHistoType argument");
      break;
  }
  return "XXX";
}

void AliRsnComparisonObj::GetESDPID(AliESDtrack *track,Double_t *pid,Double_t p)
{
  Double_t ctmp[AliRsnPID::kSpecies];
  switch (fCurrentESDPID)
  {
    case kEsd :
      track->GetESDpid(pid);
      break;
    case kITS :
      track->GetITSpid(pid);
      break;
    case kTPC :
      track->GetTPCpid(pid);
      break;
    case kTOF :
      track->GetTOFpid(pid);
      break;
    case kITS_TPC :
      track->GetITSpid(pid);
      track->GetTPCpid(ctmp);
      for (Int_t i=0;i<5;i++) pid[i]*=ctmp[i];
      break;
    case kITS_TOF :
      track->GetITSpid(pid);
      track->GetTOFpid(ctmp);
      for (Int_t i=0;i<AliRsnPID::kSpecies;i++) pid[i]*=ctmp[i];
      break;
    case kTPC_TOF :
      track->GetTPCpid(pid);
      track->GetTOFpid(ctmp);
      for (Int_t i=0;i<AliRsnPID::kSpecies;i++) pid[i]*=ctmp[i];
      break;
    case kITS_TPC_TOF :
      track->GetITSpid(pid);
      track->GetTPCpid(ctmp);
      for (Int_t i=0;i<AliRsnPID::kSpecies;i++) pid[i]*=ctmp[i];
      track->GetTOFpid(ctmp);
      for (Int_t i=0;i<AliRsnPID::kSpecies;i++) pid[i]*=ctmp[i];
      break;
    case kITS_TPC_TOF_SP :

      if (p<fPIDDivValue)
      {
        track->GetITSpid(pid);
        track->GetTPCpid(ctmp);
        for (Int_t i=0;i<AliRsnPID::kSpecies;i++) pid[i]*=ctmp[i];
      }
      else
      {
        track->GetTOFpid(pid);
      }
      break;
    default:
      AliWarning("Unrecognized value of EPIDType argument");
      for (Int_t i=0;i<AliRsnPID::kSpecies;i++) pid[i]=1.0;
      break;
  }
}

void AliRsnComparisonObj::SetESDTrackQualityCuts(const Int_t & its, const Int_t & tpc, const Int_t & trd)
{
  fITSClusters = its;
  fTPCClusters = tpc;
  fTRDClusters = trd;
}

void AliRsnComparisonObj::SetESDstatus(const ULong_t status)
{
  fESDstatus = status;
}

void AliRsnComparisonObj::SetCurrentESDPID(const EPIDType & type, const Double_t & divValue)
{
  fCurrentESDPID = type;
  fPIDDivValue = divValue;
}

TString AliRsnComparisonObj::GetParameterName(EParameterType type)
{
  switch (type)
  {
    case kP :
      return "p";
    case kPt :
      return "pt";
    case kEta :
      return "eta";
    case kY :
      return "y";
    default:
      break;
  }
  return "?";
}

Double_t AliRsnComparisonObj::GetParameterNameValue(EParameterType type,AliMCParticle * mctrack)
{
  switch (type)
  {
    case kP :
      return mctrack->P();
    case kPt :
      return mctrack->Pt();
    case kEta :
      return mctrack->Eta();
    case kY :
      return mctrack->Y();
    default:
      break;
  }
  return -10000.0;
}

