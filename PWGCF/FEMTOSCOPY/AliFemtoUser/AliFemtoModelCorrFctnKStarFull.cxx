///
/// \file AliFemtoModelCorrFctnKStarFull.cxx
/// \authors Jesse Buxton and Andrew Kubera
///

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelCorrFctnKStarFull);
  /// \endcond
#endif

#include "AliFemtoModelCorrFctnKStarFull.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelAllHiddenInfo.h"
#include "AliFemtoPair.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include <vector>
#include <string>

static const char* true_types[] =
{ "Other"     // 0
, "Primary"   // 1
, "Xi"        // 2
, "Xi0"       // 3
, "Sigma0"    // 4
, "Omega"     // 5
, "Lamba_{c}" // 6
, "Sigma_{c}" // 7
, NULL
};

static const int xi_bin = 2
        , xi0_bin = 3
        , sigma0_bin = 4
        , omega_bin = 5
        , lambda_c_bin = 6
        , sigma_c_bin = 7
        ;

static const int xi_code = 3312
        , xi0_code = 3322
        , sigma0_code = 3212
        , omega_code = 3334
        , lambda_c_code = 4122
        , sigma0_c_code = 3214
        , sigma_c_code = 3224
        ;

const size_t true_type_count = sizeof(true_types) / sizeof(true_types[0]);
const double verySmall = std::numeric_limits < double >::min();

//_________________________________________
AliFemtoModelCorrFctnKStarFull::AliFemtoModelCorrFctnKStarFull() :
  AliFemtoModelCorrFctn(),
  fTitle(""),
  fNbinsKStar(200),
  fKStarLow(0.),
  fKStarHigh(1.),
  fRemoveMisidentified(false),
  fExpectedTrack1Code(3122),
  fExpectedTrack2Code(310),
  fBuildBaseClassHistograms(true),

  fBuildUnitWeights(false),
  fNumTrueUnitWeights(NULL),
  fNumTrueIdealUnitWeights(NULL),

  fBuildParentInfo(false),
  fNumTrueIdealwParentInfo(NULL),
  fDenIdealwParentInfo(NULL),

  fBuildTrueVsRec(true),
  fKTrueVsKRecSame(NULL),
  fKTrueVsKRecMixed(NULL),

  fBuildRotated(false),
  fKTrueVsKRecRotSame(NULL),
  fKTrueVsKRecRotMixed(NULL)
{
  //default constructor
  SetBuildUnitWeights(fBuildUnitWeights);
  SetBuildParentInfo(fBuildParentInfo);
  SetBuildTrueVsRec(fBuildTrueVsRec);
  SetBuildRotated(fBuildRotated);

  //I do not use this histogram, so I might as well delete it
  delete fQgenQrec;
  fQgenQrec = NULL;

}

AliFemtoModelCorrFctnKStarFull::AliFemtoModelCorrFctnKStarFull(const char *title, int aNbins, double aKStarLo, double aKStarHi) :
  AliFemtoModelCorrFctn(title,aNbins,aKStarLo,aKStarHi),
  fTitle(title),
  fNbinsKStar(aNbins),
  fKStarLow(aKStarLo),
  fKStarHigh(aKStarHi),
  fRemoveMisidentified(false),
  fExpectedTrack1Code(3122),
  fExpectedTrack2Code(310),
  fBuildBaseClassHistograms(true),

  fBuildUnitWeights(false),
  fNumTrueUnitWeights(NULL),
  fNumTrueIdealUnitWeights(NULL),

  fBuildParentInfo(false),
  fNumTrueIdealwParentInfo(NULL),
  fDenIdealwParentInfo(NULL),

  fBuildTrueVsRec(true),
  fKTrueVsKRecSame(NULL),
  fKTrueVsKRecMixed(NULL),

  fBuildRotated(false),
  fKTrueVsKRecRotSame(NULL),
  fKTrueVsKRecRotMixed(NULL)
{
  //normal constructor
  SetBuildUnitWeights(fBuildUnitWeights);
  SetBuildParentInfo(fBuildParentInfo);
  SetBuildTrueVsRec(fBuildTrueVsRec);
  SetBuildRotated(fBuildRotated);

  //I do not use this histogram, so I might as well delete it
  delete fQgenQrec;
  fQgenQrec = NULL;
}

//_________________________________________
AliFemtoModelCorrFctnKStarFull::~AliFemtoModelCorrFctnKStarFull()
{
  //destructor
  delete fNumTrueUnitWeights;
  delete fNumTrueIdealUnitWeights;
  delete fNumTrueIdealwParentInfo;
  delete fDenIdealwParentInfo;
  delete fKTrueVsKRecSame;
  delete fKTrueVsKRecMixed;
  delete fKTrueVsKRecRotSame;
  delete fKTrueVsKRecRotMixed;
}

//_________________________________________
AliFemtoModelCorrFctnKStarFull::AliFemtoModelCorrFctnKStarFull(const AliFemtoModelCorrFctnKStarFull& aCorrFctn) :
  AliFemtoModelCorrFctn(aCorrFctn),
  fNbinsKStar(aCorrFctn.fNbinsKStar),
  fKStarLow(aCorrFctn.fKStarLow),
  fKStarHigh(aCorrFctn.fKStarHigh),
  fRemoveMisidentified(aCorrFctn.fRemoveMisidentified),
  fExpectedTrack1Code(aCorrFctn.fExpectedTrack1Code),
  fExpectedTrack2Code(aCorrFctn.fExpectedTrack2Code),
  fBuildBaseClassHistograms(aCorrFctn.fBuildBaseClassHistograms),
  fBuildUnitWeights(aCorrFctn.fBuildUnitWeights),
  fBuildParentInfo(aCorrFctn.fBuildParentInfo),
  fBuildTrueVsRec(aCorrFctn.fBuildTrueVsRec),
  fBuildRotated(aCorrFctn.fBuildRotated)

{
  //copy constructor
  // copy constructor
  if(aCorrFctn.fTitle)
  {
    fTitle = new char[strlen(aCorrFctn.fTitle)+1];
    strcpy(const_cast<char*>(fTitle), aCorrFctn.fTitle);
  }
  else fTitle = 0;

  if (aCorrFctn.fNumTrueUnitWeights)
    fNumTrueUnitWeights = new TH1F(*(aCorrFctn.fNumTrueUnitWeights));
  else fNumTrueUnitWeights = 0;

  if (aCorrFctn.fNumTrueIdealUnitWeights)
    fNumTrueIdealUnitWeights = new TH1F(*(aCorrFctn.fNumTrueIdealUnitWeights));
  else fNumTrueIdealUnitWeights = 0;

  if (aCorrFctn.fNumTrueIdealwParentInfo)
    fNumTrueIdealwParentInfo = new TH3F(*(aCorrFctn.fNumTrueIdealwParentInfo));
  else fNumTrueIdealwParentInfo = 0;

  if (aCorrFctn.fDenIdealwParentInfo)
    fDenIdealwParentInfo = new TH3F(*(aCorrFctn.fDenIdealwParentInfo));
  else fDenIdealwParentInfo = 0;

  if (aCorrFctn.fKTrueVsKRecSame)
    fKTrueVsKRecSame = new TH2F(*(aCorrFctn.fKTrueVsKRecSame));
  else fKTrueVsKRecSame = 0;

  if (aCorrFctn.fKTrueVsKRecMixed)
    fKTrueVsKRecMixed = new TH2F(*(aCorrFctn.fKTrueVsKRecMixed));
  else fKTrueVsKRecMixed = 0;

  if (aCorrFctn.fKTrueVsKRecRotSame)
    fKTrueVsKRecRotSame = new TH2F(*(aCorrFctn.fKTrueVsKRecRotSame));
  else fKTrueVsKRecRotSame = 0;

  if (aCorrFctn.fKTrueVsKRecRotMixed)
    fKTrueVsKRecRotMixed = new TH2F(*(aCorrFctn.fKTrueVsKRecRotMixed));
  else fKTrueVsKRecRotMixed = 0;

}

//_________________________________________
AliFemtoModelCorrFctnKStarFull& AliFemtoModelCorrFctnKStarFull::operator=(const AliFemtoModelCorrFctnKStarFull& aCorrFctn)
{
  //assignment operator
  if (this == &aCorrFctn) return *this;

  AliFemtoModelCorrFctn::operator=(aCorrFctn);

  fNbinsKStar = aCorrFctn.fNbinsKStar;
  fKStarLow = aCorrFctn.fKStarLow;
  fKStarHigh = aCorrFctn.fKStarHigh;

  fRemoveMisidentified = aCorrFctn.fRemoveMisidentified;
  fExpectedTrack1Code = aCorrFctn.fExpectedTrack1Code;
  fExpectedTrack2Code = aCorrFctn.fExpectedTrack2Code;
  fBuildBaseClassHistograms = aCorrFctn.fBuildBaseClassHistograms;

  fBuildUnitWeights = aCorrFctn.fBuildUnitWeights;
  fBuildParentInfo = aCorrFctn.fBuildParentInfo;
  fBuildTrueVsRec = aCorrFctn.fBuildTrueVsRec;
  fBuildRotated = aCorrFctn.fBuildRotated;

  if(aCorrFctn.fTitle)
  {
    fTitle = new char[strlen(aCorrFctn.fTitle)+1];
    strcpy(const_cast<char*>(fTitle), aCorrFctn.fTitle);
  }
  else fTitle = 0;

  if (aCorrFctn.fNumTrueUnitWeights)
    fNumTrueUnitWeights = new TH1F(*(aCorrFctn.fNumTrueUnitWeights));
  else fNumTrueUnitWeights = 0;

  if (aCorrFctn.fNumTrueIdealUnitWeights)
    fNumTrueIdealUnitWeights = new TH1F(*(aCorrFctn.fNumTrueIdealUnitWeights));
  else fNumTrueIdealUnitWeights = 0;

  if (aCorrFctn.fNumTrueIdealwParentInfo)
    fNumTrueIdealwParentInfo = new TH3F(*(aCorrFctn.fNumTrueIdealwParentInfo));
  else fNumTrueIdealwParentInfo = 0;

  if (aCorrFctn.fDenIdealwParentInfo)
    fDenIdealwParentInfo = new TH3F(*(aCorrFctn.fDenIdealwParentInfo));
  else fDenIdealwParentInfo = 0;

  if (aCorrFctn.fKTrueVsKRecSame)
    fKTrueVsKRecSame = new TH2F(*(aCorrFctn.fKTrueVsKRecSame));
  else fKTrueVsKRecSame = 0;

  if (aCorrFctn.fKTrueVsKRecMixed)
    fKTrueVsKRecMixed = new TH2F(*(aCorrFctn.fKTrueVsKRecMixed));
  else fKTrueVsKRecMixed = 0;

  if (aCorrFctn.fKTrueVsKRecRotSame)
    fKTrueVsKRecRotSame = new TH2F(*(aCorrFctn.fKTrueVsKRecRotSame));
  else fKTrueVsKRecRotSame = 0;

  if (aCorrFctn.fKTrueVsKRecRotMixed)
    fKTrueVsKRecRotMixed = new TH2F(*(aCorrFctn.fKTrueVsKRecRotMixed));
  else fKTrueVsKRecRotMixed = 0;

  return *this;
}

//_________________________________________
void AliFemtoModelCorrFctnKStarFull::AddRealPair(AliFemtoPair* aPair)
{
  double tWeight;
  if(fManager) tWeight = fManager->GetWeight(aPair);
  else tWeight = 1.0;

  double tKStarTrue = GetKStarTrue(aPair);
  double tKStarRec = aPair->KStar();

//  if(tWeight != 0 && tKStarTrue != 0 && tKStarTrue > -999)  //make sure we have a good particle
  if(tWeight > 0. && tKStarTrue > 0. && tKStarRec > 0.)  //make sure we have a good particle
  {
    if(fBuildBaseClassHistograms)
    {
      fNumeratorTrue->Fill(tKStarRec,tWeight);
      fNumeratorTrueIdeal->Fill(tKStarTrue,tWeight);
    }

    if(fBuildUnitWeights)
    {
      fNumTrueUnitWeights->Fill(tKStarRec,1.0);
      fNumTrueIdealUnitWeights->Fill(tKStarTrue,1.0);
    }

    if(fBuildTrueVsRec) fKTrueVsKRecSame->Fill(tKStarTrue,tKStarRec);

    if(fBuildRotated)
    {
      double tX = 0.5*(tKStarRec + tKStarTrue);
      double tY = (tKStarRec - tKStarTrue)/sqrt(2);
      fKTrueVsKRecRotSame->Fill(tX,tY);
    }

    if(fBuildParentInfo)
    {
      int tMotherBin1 = GetMotherBin((AliFemtoModelHiddenInfo*)aPair->Track1()->GetHiddenInfo());
      int tMotherBin2 = GetMotherBin((AliFemtoModelHiddenInfo*)aPair->Track2()->GetHiddenInfo());
      fNumTrueIdealwParentInfo->Fill(tKStarTrue,tMotherBin1,tMotherBin2);
    }
  }
}

//_________________________________________
void AliFemtoModelCorrFctnKStarFull::AddMixedPair(AliFemtoPair* aPair)
{
  double tWeight;
  if(fManager) tWeight = fManager->GetWeight(aPair);
  else tWeight = 1.0;

  double tKStarTrue = GetKStarTrue(aPair);
  double tKStarRec = aPair->KStar();

//  if(tWeight > 0 && tKStarTrue > 0 && tKStarTrue > -999)  //make sure we have a good particle
  if(tWeight > 0. && tKStarTrue > 0. && tKStarRec > 0.)  //make sure we have a good particle
  {
    if(fBuildBaseClassHistograms)
    {
      fNumeratorFake->Fill(tKStarRec,tWeight);
      fDenominator->Fill(tKStarRec, 1.0);

      fNumeratorFakeIdeal->Fill(tKStarTrue,tWeight);
      fDenominatorIdeal->Fill(tKStarTrue, 1.0);
    }

    if(fBuildTrueVsRec) fKTrueVsKRecMixed->Fill(tKStarTrue,tKStarRec);

    if(fBuildRotated)
    {
      double tX = 0.5*(tKStarRec + tKStarTrue);
      double tY = (tKStarRec - tKStarTrue)/sqrt(2);
      fKTrueVsKRecRotMixed->Fill(tX,tY);
    }

    if(fBuildParentInfo)
    {
      int tMotherBin1 = GetMotherBin((AliFemtoModelHiddenInfo*)aPair->Track1()->GetHiddenInfo());
      int tMotherBin2 = GetMotherBin((AliFemtoModelHiddenInfo*)aPair->Track2()->GetHiddenInfo());
      fDenIdealwParentInfo->Fill(tKStarTrue,tMotherBin1,tMotherBin2);
    }
  }
}

//_______________________
void AliFemtoModelCorrFctnKStarFull::Write()
{
  // Write out data histos
  if(fBuildBaseClassHistograms)
  {
    fNumeratorTrue->Write();
    fNumeratorFake->Write();
    fDenominator->Write();

    fNumeratorTrueIdeal->Write();
    fNumeratorFakeIdeal->Write();
    fDenominatorIdeal->Write();
  }

  if(fBuildUnitWeights)
  {
    fNumTrueUnitWeights->Write();
    fNumTrueIdealUnitWeights->Write();
  }
  if(fBuildParentInfo)
  {
    fNumTrueIdealwParentInfo->Write();
    fDenIdealwParentInfo->Write();
  }

  if(fBuildTrueVsRec)
  {
    fKTrueVsKRecSame->Write();
    fKTrueVsKRecMixed->Write();
  }
  if(fBuildRotated)
  {
    fKTrueVsKRecRotSame->Write();
    fKTrueVsKRecRotMixed->Write();
  }
}

//_________________________
TList* AliFemtoModelCorrFctnKStarFull::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  if(fBuildBaseClassHistograms)
  {
    tOutputList->Add(fNumeratorTrue);
    tOutputList->Add(fNumeratorFake);
    tOutputList->Add(fDenominator);

    tOutputList->Add(fNumeratorTrueIdeal);
    tOutputList->Add(fNumeratorFakeIdeal);
    tOutputList->Add(fDenominatorIdeal);
  }

  if(fBuildUnitWeights)
  {
    tOutputList->Add(fNumTrueUnitWeights);
    tOutputList->Add(fNumTrueIdealUnitWeights);
  }
  if(fBuildParentInfo)
  {
    tOutputList->Add(fNumTrueIdealwParentInfo);
    tOutputList->Add(fDenIdealwParentInfo);
  }

  if(fBuildTrueVsRec)
  {
    tOutputList->Add(fKTrueVsKRecSame);
    tOutputList->Add(fKTrueVsKRecMixed);
  }
  if(fBuildRotated)
  {
    tOutputList->Add(fKTrueVsKRecRotSame);
    tOutputList->Add(fKTrueVsKRecRotMixed);
  }

  return tOutputList;
}


//_________________________________________
void AliFemtoModelCorrFctnKStarFull::SetExpectedPDGCodes(const Int_t track1_pdg_code,
                                                     const Int_t track2_pdg_code)
{
  fExpectedTrack1Code = track1_pdg_code;
  fExpectedTrack2Code = track2_pdg_code;
}

//_________________________________________
bool AliFemtoModelCorrFctnKStarFull::PairContainsExpectedTypes(const AliFemtoPair *pair)
{
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelAllHiddenInfo*>(pair->Track1()->HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelAllHiddenInfo*>(pair->Track2()->HiddenInfo());

  return PairContainsExpectedTypes(info1, info2);
}

//_________________________________________
bool AliFemtoModelCorrFctnKStarFull::PairContainsExpectedTypes(
  const AliFemtoModelHiddenInfo *info1,
  const AliFemtoModelHiddenInfo *info2)
{
  const bool track_1_is_expected = fExpectedTrack1Code == info1->GetPDGPid(),
             track_2_is_expected = fExpectedTrack2Code == info2->GetPDGPid();

  return track_1_is_expected && track_2_is_expected;
}

//_________________________________________
double AliFemtoModelCorrFctnKStarFull::CalcKStar(const AliFemtoLorentzVector &p1, const AliFemtoLorentzVector &p2)
{
  const double p_inv = (p1 + p2).m2(),
               q_inv = (p1 - p2).m2(),
           mass_diff = p1.m2() - p2.m2();

  const double tQ = ::pow(mass_diff, 2) / p_inv - q_inv;
  return ::sqrt(tQ) / 2.0;
}

//_________________________________________
double AliFemtoModelCorrFctnKStarFull::GetKStarTrue(AliFemtoPair* aPair)
{
  // Generate a simple femtoscopic weight coming simple Lednicky equation
  // The weight is generated using the TrueMomentum in the HiddenInfo

  AliFemtoParticle *tPart1 = (AliFemtoParticle*)aPair->Track1();
  AliFemtoParticle *tPart2 = (AliFemtoParticle*)aPair->Track2();
  if(tPart1 == NULL || tPart2 == NULL) return -9999.0;

  AliFemtoModelHiddenInfo *tPart1HiddenInfo = (AliFemtoModelHiddenInfo*)tPart1->GetHiddenInfo();
  AliFemtoModelHiddenInfo *tPart2HiddenInfo = (AliFemtoModelHiddenInfo*)tPart2->GetHiddenInfo();
  if(tPart1HiddenInfo == NULL || tPart2HiddenInfo == NULL) return -9999.0;

  if(fRemoveMisidentified)
  {
    bool tPass = PairContainsExpectedTypes(tPart1HiddenInfo,tPart2HiddenInfo);
    if(!tPass) return -9999.0;
  }

  AliFemtoThreeVector *tP1Vec = tPart1HiddenInfo->GetTrueMomentum();
  AliFemtoThreeVector *tP2Vec = tPart2HiddenInfo->GetTrueMomentum();
  if(tP1Vec->Mag2() < verySmall || tP2Vec->Mag2() < verySmall) return -9999.0; //p1 or p2 has zero momentum!

  double mass1 = tPart1HiddenInfo->GetMass();
  double mass2 = tPart2HiddenInfo->GetMass();

  double E1 = std::sqrt(mass1*mass1 + tP1Vec->Mag2());
  double E2 = std::sqrt(mass2*mass2 + tP2Vec->Mag2());

  return CalcKStar(AliFemtoLorentzVector(E1, *tP1Vec), AliFemtoLorentzVector(E2, *tP2Vec));
}


//_________________________________________
int AliFemtoModelCorrFctnKStarFull::GetMotherBin(const AliFemtoModelHiddenInfo *info)
{
  switch (abs(info->GetMotherPdgCode())) {
    default: return 0;
    case 0: return 1;
    case xi_code: return xi_bin;
    case xi0_code: return xi0_bin;
    case sigma0_code: return sigma0_bin;
    case omega_code: return omega_bin;
    case lambda_c_code: return lambda_c_bin;
    case sigma0_c_code:
    case sigma_c_code: return sigma_c_bin;
  }
}

//_________________________________________
void AliFemtoModelCorrFctnKStarFull::SetBuildBaseClassHistograms(bool aBuild)
{
  fBuildBaseClassHistograms = aBuild;
}

//_________________________________________
void AliFemtoModelCorrFctnKStarFull::SetBuildUnitWeights(bool aBuild)
{
  fBuildUnitWeights = aBuild;

  if(fBuildUnitWeights && (!fNumTrueUnitWeights || !fNumTrueIdealUnitWeights))
  {
    fNumTrueUnitWeights = new TH1F(
      TString::Format("NumTrueUnitWeights%s", fTitle),
      "KStar - True Numerator w/ Unit Weights; k*(GeV);",
      fNbinsKStar, fKStarLow, fKStarHigh);
    fNumTrueUnitWeights->Sumw2();

    fNumTrueIdealUnitWeights = new TH1F(
      TString::Format("NumTrueIdealUnitWeights%s", fTitle),
      "KStar - True Ideal Numerator w/ Unit Weights; k*(GeV);",
      fNbinsKStar, fKStarLow, fKStarHigh);
    fNumTrueIdealUnitWeights->Sumw2();
  }
}

//_________________________________________
void AliFemtoModelCorrFctnKStarFull::SetBuildParentInfo(bool aBuild)
{
  fBuildParentInfo = aBuild;

  if(fBuildParentInfo && (!fNumTrueIdealwParentInfo || !fDenIdealwParentInfo))
  {
    fNumTrueIdealwParentInfo = new TH3F(
      TString::Format("NumTrueIdealwParentInfo%s", fTitle),
      "KStar - True Ideal Numerator binned in parent PIDs; k*(GeV);",
      fNbinsKStar, fKStarLow, fKStarHigh,
      true_type_count, -0.5, true_type_count - 0.5,
      true_type_count, -0.5, true_type_count - 0.5);
    fNumTrueIdealwParentInfo->Sumw2();

    fDenIdealwParentInfo = new TH3F(
      TString::Format("DenIdealwParentInfo%s", fTitle),
      "KStar - Ideal Denominator binned in parent PIDs; k*(GeV);",
      fNbinsKStar, fKStarLow, fKStarHigh,
      true_type_count, -0.5, true_type_count - 0.5,
      true_type_count, -0.5, true_type_count - 0.5);
    fDenIdealwParentInfo->Sumw2();

    for (size_t i = 0; i < true_type_count; ++i) {
      fNumTrueIdealwParentInfo->GetYaxis()->SetBinLabel(i+1, true_types[i]);
      fDenIdealwParentInfo->GetYaxis()->SetBinLabel(i+1, true_types[i]);
    }
  }
}

//_________________________________________
void AliFemtoModelCorrFctnKStarFull::SetBuildTrueVsRec(bool aBuild)
{
  fBuildTrueVsRec = aBuild;

  if(fBuildTrueVsRec && (!fKTrueVsKRecSame || !fKTrueVsKRecMixed))
  {
    fKTrueVsKRecSame = new TH2F(
      TString::Format("KTrueVsKRecSame%s", fTitle), 
      "k*_{true} vs k*_{rec} for Same Event Pairs; k*(GeV/c);",
       fNbinsKStar, fKStarLow, fKStarHigh,
       fNbinsKStar, fKStarLow, fKStarHigh);
    fKTrueVsKRecSame->Sumw2();

    fKTrueVsKRecMixed = new TH2F(
      TString::Format("KTrueVsKRecMixed%s", fTitle), 
      "k*_{true} vs k*_{rec} for Mixed Event Pairs; k*(GeV/c);",
       fNbinsKStar, fKStarLow, fKStarHigh,
       fNbinsKStar, fKStarLow, fKStarHigh);
    fKTrueVsKRecMixed->Sumw2();
  }
}


//_________________________________________
void AliFemtoModelCorrFctnKStarFull::SetBuildRotated(bool aBuild)
{
  fBuildRotated = aBuild;

  if(fBuildRotated && (!fKTrueVsKRecRotSame || !fKTrueVsKRecRotMixed))
  {

    double tYRange = fKStarHigh-fKStarLow;
    double tRotYMin = -0.5*tYRange;
    double tRotYMax = 0.5*tYRange;

    fKTrueVsKRecRotSame = new TH2F(
      TString::Format("KTrueVsKRecRotSame%s", fTitle), 
      "(k*_{true} + k*_{rec})/2 vs (k*_{true} - k*_{rec})/sqrt(2) for Same Event Pairs; k*(GeV/c);",
       fNbinsKStar, fKStarLow, fKStarHigh,
       fNbinsKStar, tRotYMin, tRotYMax);
    fKTrueVsKRecRotSame->Sumw2();

    fKTrueVsKRecRotMixed = new TH2F(
      TString::Format("KTrueVsKRecRotMixed%s", fTitle), 
      "(k*_{true} + k*_{rec})/2 vs (k*_{true} - k*_{rec})/sqrt(2) for Mixed Event Pairs; k*(GeV/c);",
       fNbinsKStar, fKStarLow, fKStarHigh,
       fNbinsKStar, tRotYMin, tRotYMax);
    fKTrueVsKRecRotMixed->Sumw2();
  }
}


