#include "AliFemtoV0TrackCutNSigmaFilter.h"
#include "AliESDtrack.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoV0TrackCutNSigmaFilter);
  /// \endcond
#endif


AliFemtoV0TrackCutNSigmaFilter::AliFemtoV0TrackCutNSigmaFilter() :
  AliFemtoV0TrackCut(),
  fUseCustomPionNSigmaFilter(false),
  fUseCustomKaonNSigmaFilter(false),
  fUseCustomProtonNSigmaFilter(false),

  fPionNSigmaFilter(0),
  fKaonNSigmaFilter(0),
  fProtonNSigmaFilter(0),

  fBuildMinvHisto(false),
  fMinvHisto(0)
{
  // Default constructor
}
//------------------------------
AliFemtoV0TrackCutNSigmaFilter::~AliFemtoV0TrackCutNSigmaFilter()
{
  /* noop */
}
//------------------------------
bool AliFemtoV0TrackCutNSigmaFilter::Pass(const AliFemtoV0* aV0)
{
  // test the particle and return
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria

  Float_t pt = aV0->PtV0();
  Float_t eta = aV0->EtaV0();

  //kinematic cuts
  if (TMath::Abs(eta) > fEta) return false;  //put in kinematic cuts by hand
  if (pt < fPtMin || fPtMax < pt) return false;
  if (TMath::Abs(aV0->EtaPos()) > fMaxEtaDaughters) return false;
  if (TMath::Abs(aV0->EtaNeg()) > fMaxEtaDaughters) return false;

  if (aV0->PtPos() < fPtMinPosDaughter) return false;
  if (aV0->PtNeg() < fPtMinNegDaughter) return false;
  if (aV0->PtPos() > fPtMaxPosDaughter) return false;
  if (aV0->PtNeg() > fPtMaxNegDaughter) return false;

  //V0 from kinematics information
  if (fParticleType == kLambdaMC ) {
    if (!(aV0->MassLambda() > fInvMassLambdaMin && aV0->MassLambda() < fInvMassLambdaMax) || !(aV0->PosNSigmaTPCP() == 0)) {
      return false;
    } else {
      return true;
    }
  } else if (fParticleType == kAntiLambdaMC) {
    if (!(aV0->MassLambda() > fInvMassLambdaMin && aV0->MassLambda() < fInvMassLambdaMax) || !(aV0->NegNSigmaTPCP() == 0)) {
      return false;
    } else {
      return true;
    }
  } else if (fParticleType == kAll) {
    if (!(aV0->MassK0Short() > fInvMassK0sMin && aV0->MassK0Short() < fInvMassK0sMax) || !(aV0->NegNSigmaTPCP() == 0))
      return false;
    else {
      return true;
    }
  }

  //quality cuts
  if (aV0->OnFlyStatusV0() != fOnFlyStatus) return false;
  if (aV0->StatusNeg() == 999 || aV0->StatusPos() == 999) return false;
  if (aV0->TPCNclsPos() < fTPCNclsDaughters) return false;
  if (aV0->TPCNclsNeg() < fTPCNclsDaughters) return false;
  if (aV0->NdofPos() > fNdofDaughters) return false;
  if (aV0->NdofNeg() > fNdofDaughters) return false;
  if (!(aV0->StatusNeg() & fStatusDaughters)) return false;
  if (!(aV0->StatusPos() & fStatusDaughters)) return false;


  //DCA between daughter particles
  if (TMath::Abs(aV0->DcaV0Daughters()) > fMaxDcaV0Daughters)
    return false;

  //DCA of daughters to primary vertex
  if (TMath::Abs(aV0->DcaPosToPrimVertex()) < fMinDcaDaughterPosToVert || TMath::Abs(aV0->DcaNegToPrimVertex()) < fMinDcaDaughterNegToVert)
    return false;

  //DCA V0 to prim vertex
  if (TMath::Abs(aV0->DcaV0ToPrimVertex()) > fMaxDcaV0 || TMath::Abs(aV0->DcaV0ToPrimVertex()) < fMinDcaV0)
    return false;

  //becomes obsolete - wrong name of the data memeber and the corresnponding methods (by default is set to fMaxCosPointingAngle = 0.0)
  //cos pointing angle
  if (aV0->CosPointingAngle() < fMaxCosPointingAngle)
    return false;

  //this is the correct name of the data member and the corresponding methods (we are accepting cos(pointing angle bigger than certain minimum)
  //cos pointing angle
  if (aV0->CosPointingAngle() < fMinCosPointingAngle)
    return false;

  //decay length
  if (aV0->DecayLengthV0() > fMaxDecayLength)
    return false;


  if (fParticleType == kAll)
    return true;


  bool pid_check = false;
  // Looking for lambdas = proton + pim
  if (fParticleType == kLambda) {
    if (IsProtonNSigma(aV0->PtPos(), aV0->PosNSigmaTPCP(), aV0->PosNSigmaTOFP())) //proton
      if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) { //pion-
        pid_check = true;
        if(fBuildMinvHisto) fMinvHisto->Fill(aV0->MassLambda());
        //invariant mass lambda
        if (aV0->MassLambda() < fInvMassLambdaMin || aV0->MassLambda() > fInvMassLambdaMax)
          return false;
      }

  }//Looking for antilambdas =  antiproton + pip
  else if (fParticleType == kAntiLambda) {
    if (IsProtonNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP())) //antiproton
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) { //pion+
        pid_check = true;
        if(fBuildMinvHisto) fMinvHisto->Fill(aV0->MassAntiLambda());
        //invariant mass antilambda
        if (aV0->MassAntiLambda() < fInvMassLambdaMin || aV0->MassAntiLambda() > fInvMassLambdaMax)
          return false;
      }
  }//Looking for K0s = pip + pim
  else if (fParticleType == kK0s) {
    if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) //pion-
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) { //pion+
        pid_check = true;
        if(fBuildMinvHisto) fMinvHisto->Fill(aV0->MassK0Short());
        //invariant mass K0s
        if (aV0->MassK0Short() < fInvMassK0sMin || aV0->MassK0Short() > fInvMassK0sMax)
          return false;
      }
  }

  if (!pid_check) return false;

  return true;
}
//------------------------------
AliFemtoString AliFemtoV0TrackCutNSigmaFilter::Report()
{
  // Prepare report from the execution
  AliFemtoString returnThis = AliFemtoV0TrackCut::Report();

  string tStemp;
  char tCtemp[100];
  snprintf(tCtemp, 100, "Usings custom Pion NSigma Filter:\t%i\n", fUseCustomPionNSigmaFilter);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Usings custom Kaon NSigma Filter:\t%i\n", fUseCustomKaonNSigmaFilter);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Usings custom Proton NSigma Filter:\t%i\n", fUseCustomProtonNSigmaFilter);
  tStemp += tCtemp;

  returnThis += tStemp;

  return returnThis;
}
TList *AliFemtoV0TrackCutNSigmaFilter::ListSettings()
{
  TList *tListSettings = AliFemtoV0TrackCut::ListSettings();

  //Add any additional information to tLIstSettings

  return tListSettings;
}


bool AliFemtoV0TrackCutNSigmaFilter::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if(fUseCustomKaonNSigmaFilter) {return fKaonNSigmaFilter->Pass(mom,nsigmaTPCK,nsigmaTOFK);}
  else {return AliFemtoV0TrackCut::IsKaonNSigma(mom,nsigmaTPCK,nsigmaTOFK);}  
}

bool AliFemtoV0TrackCutNSigmaFilter::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if(fUseCustomPionNSigmaFilter) {return fPionNSigmaFilter->Pass(mom,nsigmaTPCPi,nsigmaTOFPi);}
  else {return AliFemtoV0TrackCut::IsPionNSigma(mom,nsigmaTPCPi,nsigmaTOFPi);}
}

bool AliFemtoV0TrackCutNSigmaFilter::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
  if(fUseCustomProtonNSigmaFilter) {return fProtonNSigmaFilter->Pass(mom,nsigmaTPCP,nsigmaTOFP);}
  else {return AliFemtoV0TrackCut::IsProtonNSigma(mom,nsigmaTPCP,nsigmaTOFP);}
}




//!!!!!----- 14/12/2015 ---------------------------------------------------------------------------
void AliFemtoV0TrackCutNSigmaFilter::CreateCustomNSigmaFilter(DaughterParticleType aDaughterType)
{
  if(aDaughterType == kPion)
  {
    fUseCustomPionNSigmaFilter = true;
    fPionNSigmaFilter = new AliFemtoNSigmaFilter();
  }
  else if(aDaughterType == kKaon)
  {
    fUseCustomKaonNSigmaFilter = true;
    fKaonNSigmaFilter = new AliFemtoNSigmaFilter();
  }
  else if(aDaughterType == kProton)
  {
    fUseCustomProtonNSigmaFilter = true;
    fProtonNSigmaFilter = new AliFemtoNSigmaFilter();
  }
  else {/*cerr << "Invalid DaughterParticleType selection in CreateCustomNSigmaFilter.  No custom filter will be initialized!!!!!" << endl;*/}

} 
void AliFemtoV0TrackCutNSigmaFilter::CreateCustomPionNSigmaFilter() {CreateCustomNSigmaFilter(kPion);}
void AliFemtoV0TrackCutNSigmaFilter::CreateCustomKaonNSigmaFilter() {CreateCustomNSigmaFilter(kKaon);}
void AliFemtoV0TrackCutNSigmaFilter::CreateCustomProtonNSigmaFilter() {CreateCustomNSigmaFilter(kProton);}


void AliFemtoV0TrackCutNSigmaFilter::AddTPCAndTOFNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  if(aDaughterType == kPion) {fPionNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
  else if(aDaughterType == kKaon) {fKaonNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
  else if(aDaughterType == kProton) {fProtonNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
  else {/*cerr << "Invalid DaughterParticleType selection in AddTPCAndTOFNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}
void AliFemtoV0TrackCutNSigmaFilter::AddPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
  {AddTPCAndTOFNSigmaCut(kPion,aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
void AliFemtoV0TrackCutNSigmaFilter::AddKaonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
  {AddTPCAndTOFNSigmaCut(kKaon,aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}
void AliFemtoV0TrackCutNSigmaFilter::AddProtonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
  {AddTPCAndTOFNSigmaCut(kProton,aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);}



void AliFemtoV0TrackCutNSigmaFilter::AddTPCNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  if(aDaughterType == kPion) {fPionNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);}
  else if(aDaughterType == kKaon) {fKaonNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);}
  else if(aDaughterType == kProton) {fProtonNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);}
  else {/*cerr << "Invalid DaughterParticleType selection in AddTPCNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}
void AliFemtoV0TrackCutNSigmaFilter::AddPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
  {AddTPCNSigmaCut(kPion,aMomMin,aMomMax,aNSigmaValueTPC);}
void AliFemtoV0TrackCutNSigmaFilter::AddKaonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
  {AddTPCNSigmaCut(kKaon,aMomMin,aMomMax,aNSigmaValueTPC);}
void AliFemtoV0TrackCutNSigmaFilter::AddProtonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
  {AddTPCNSigmaCut(kProton,aMomMin,aMomMax,aNSigmaValueTPC);}



void AliFemtoV0TrackCutNSigmaFilter::AddTOFNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  if(aDaughterType == kPion) {fPionNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);}
  else if(aDaughterType == kKaon) {fKaonNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);}
  else if(aDaughterType == kProton) {fProtonNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);}
  else {/*cerr << "Invalid DaughterParticleType selection in AddTOFNSigmaCut.  No cut will be added to the filter!!!!!" << endl;*/}
}
void AliFemtoV0TrackCutNSigmaFilter::AddPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
  {AddTOFNSigmaCut(kPion,aMomMin,aMomMax,aNSigmaValueTOF);}
void AliFemtoV0TrackCutNSigmaFilter::AddKaonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
  {AddTOFNSigmaCut(kKaon,aMomMin,aMomMax,aNSigmaValueTOF);}
void AliFemtoV0TrackCutNSigmaFilter::AddProtonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
  {AddTOFNSigmaCut(kProton,aMomMin,aMomMax,aNSigmaValueTOF);}


void AliFemtoV0TrackCutNSigmaFilter::SetMinvHisto(const char* title, const int& nbins, const float& aInvMassMin, const float& aInvMassMax)
{
  fBuildMinvHisto = true;
  fMinvHisto = new TH1D(title,"MinvHistogram",nbins,aInvMassMin,aInvMassMax);
  fMinvHisto->Sumw2();
}

TH1D* AliFemtoV0TrackCutNSigmaFilter::GetMinvHisto() {return fMinvHisto;}


