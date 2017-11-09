#include "AliFemtoV0TrackCut.h"
#include "AliESDtrack.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoV0TrackCut);
  /// \endcond
#endif

//------------------------------
AliFemtoV0TrackCut::AliFemtoV0TrackCut():
  fInvMassLambdaMin(0.0),
  fInvMassLambdaMax(99.0),
  fInvMassK0sMin(0.0),
  fInvMassK0sMax(99.0),

  fMinDcaDaughterPosToVert(0.0),
  fMinDcaDaughterNegToVert(0.0),
  fMaxDcaV0Daughters(99.0),
  fMaxDcaV0(9999.0),
  fMinDcaV0(0.0),
  fMaxDecayLength(9999.0),
  fMaxCosPointingAngle(0.0),
  fMinCosPointingAngle(0.0),
  fParticleType(99.0),
  fEta(0.8),
  fPtMin(0.0),
  fPtMax(100.0),
  fOnFlyStatus(kFALSE),
  fMaxEtaDaughters(100.0),
  fTPCNclsDaughters(0),
  fNdofDaughters(10),
  fStatusDaughters(0),
  fPtMinPosDaughter(0.0),
  fPtMaxPosDaughter(99.0),
  fPtMinNegDaughter(0.0),
  fPtMaxNegDaughter(99.0),

  
  fMinAvgSepDaughters(0),


  fNsigmaPosDaughterTPC(3.0),
  fNsigmaNegDaughterTPC(3.0),
  fNsigmaPosDaughterTOF(3.0),
  fNsigmaNegDaughterTOF(3.0),

  fRadiusV0Min(0.0),
  fRadiusV0Max(99999.0),
  fRequireTOFPion(false),
  fRequireTOFProton(true),
  
  fBuildPurityAidV0(false),
  fMinvPurityAidHistoV0(0),

  fUseLooseInvMassCut(false),
  fLooseInvMassMin(0),
  fLooseInvMassMax(99),

  fRemoveMisidentified(true),
  fUseSimpleMisIDCut(true),
  fBuildMisIDHistograms(false),
  fInvMassRejectK0sMin(0.0),
  fInvMassRejectK0sMax(0.0),
  fInvMassRejectLambdaMin(0.0),
  fInvMassRejectLambdaMax(0.0),
  fInvMassRejectAntiLambdaMin(0.0),
  fInvMassRejectAntiLambdaMax(0.0),
  fK0sMassOfMisIDV0(0),
  fLambdaMassOfMisIDV0(0),
  fAntiLambdaMassOfMisIDV0(0),
  fIgnoreOnFlyStatus(false)

{
  // Default constructor
}

//------------------------------
AliFemtoV0TrackCut::~AliFemtoV0TrackCut()
{
  /* noop */
}

//------------------------------
AliFemtoV0TrackCut::AliFemtoV0TrackCut(const AliFemtoV0TrackCut& aCut) :
  AliFemtoParticleCut(aCut),

  fInvMassLambdaMin(aCut.fInvMassLambdaMin),
  fInvMassLambdaMax(aCut.fInvMassLambdaMax),
  fInvMassK0sMin(aCut.fInvMassK0sMin),
  fInvMassK0sMax(aCut.fInvMassK0sMax),

  fMinDcaDaughterPosToVert(aCut.fMinDcaDaughterPosToVert),
  fMinDcaDaughterNegToVert(aCut.fMinDcaDaughterNegToVert),
  fMaxDcaV0Daughters(aCut.fMaxDcaV0Daughters),
  fMaxDcaV0(aCut.fMaxDcaV0),
  fMinDcaV0(aCut.fMinDcaV0),
  fMaxDecayLength(aCut.fMaxDecayLength),
  fMinTransverseDistancePrimSecVtx(aCut.fMinTransverseDistancePrimSecVtx),
  
  fMaxCosPointingAngle(aCut.fMaxCosPointingAngle),
  fMinCosPointingAngle(aCut.fMinCosPointingAngle),
  fParticleType(aCut.fParticleType),
  fEta(aCut.fEta),
  fPtMin(aCut.fPtMin),
  fPtMax(aCut.fPtMax),
  fOnFlyStatus(aCut.fOnFlyStatus),

  fMaxEtaDaughters(aCut.fMaxEtaDaughters),
  fTPCNclsDaughters(aCut.fTPCNclsDaughters),
  fNdofDaughters(aCut.fNdofDaughters),
  fStatusDaughters(aCut.fStatusDaughters),
  fPtMinPosDaughter(aCut.fPtMinPosDaughter),
  fPtMaxPosDaughter(aCut.fPtMaxPosDaughter),
  fPtMinNegDaughter(aCut.fPtMinNegDaughter),
  fPtMaxNegDaughter(aCut.fPtMaxNegDaughter),
  fMinAvgSepDaughters(aCut.fMinAvgSepDaughters),

  fNsigmaPosDaughterTPC(aCut.fNsigmaPosDaughterTPC),
  fNsigmaNegDaughterTPC(aCut.fNsigmaNegDaughterTPC),
  fNsigmaPosDaughterTOF(aCut.fNsigmaPosDaughterTOF),
  fNsigmaNegDaughterTOF(aCut.fNsigmaNegDaughterTOF),

  fRadiusV0Min(aCut.fRadiusV0Min),
  fRadiusV0Max(aCut.fRadiusV0Max),

  fRequireTOFPion(aCut.fRequireTOFPion),
  fRequireTOFProton(aCut.fRequireTOFProton),

  fBuildPurityAidV0(aCut.fBuildPurityAidV0),

  fUseLooseInvMassCut(aCut.fUseLooseInvMassCut),
  fLooseInvMassMin(aCut.fLooseInvMassMin),
  fLooseInvMassMax(aCut.fLooseInvMassMax),

  fRemoveMisidentified(aCut.fRemoveMisidentified),
  fUseSimpleMisIDCut(aCut.fUseSimpleMisIDCut),
  fBuildMisIDHistograms(aCut.fBuildMisIDHistograms),
  fInvMassRejectK0sMin(aCut.fInvMassRejectK0sMin),
  fInvMassRejectK0sMax(aCut.fInvMassRejectK0sMax),
  fInvMassRejectLambdaMin(aCut.fInvMassRejectLambdaMin),
  fInvMassRejectLambdaMax(aCut.fInvMassRejectLambdaMax),
  fInvMassRejectAntiLambdaMin(aCut.fInvMassRejectAntiLambdaMin),
  fInvMassRejectAntiLambdaMax(aCut.fInvMassRejectAntiLambdaMax),
  fIgnoreOnFlyStatus(aCut.fIgnoreOnFlyStatus)
{
  //copy constructor
  if(aCut.fMinvPurityAidHistoV0) fMinvPurityAidHistoV0 = new TH1D(*aCut.fMinvPurityAidHistoV0);
    else fMinvPurityAidHistoV0 = 0;
  if(aCut.fK0sMassOfMisIDV0) fK0sMassOfMisIDV0 = new TH1D(*aCut.fK0sMassOfMisIDV0);
    else fK0sMassOfMisIDV0 = 0;
  if(aCut.fLambdaMassOfMisIDV0) fLambdaMassOfMisIDV0 = new TH1D(*aCut.fLambdaMassOfMisIDV0);
    else fLambdaMassOfMisIDV0 = 0;
  if(aCut.fAntiLambdaMassOfMisIDV0) fAntiLambdaMassOfMisIDV0 = new TH1D(*aCut.fAntiLambdaMassOfMisIDV0);
    else fAntiLambdaMassOfMisIDV0 = 0;
}

//------------------------------
AliFemtoV0TrackCut& AliFemtoV0TrackCut::operator=(const AliFemtoV0TrackCut& aCut)
{
  //assignment operator
  if (this == &aCut) return *this;

  AliFemtoParticleCut::operator=(aCut);

  fInvMassLambdaMin = aCut.fInvMassLambdaMin;
  fInvMassLambdaMax = aCut.fInvMassLambdaMax;
  fInvMassK0sMin = aCut.fInvMassK0sMin;
  fInvMassK0sMax = aCut.fInvMassK0sMax;

  fMinDcaDaughterPosToVert = aCut.fMinDcaDaughterPosToVert;
  fMinDcaDaughterNegToVert = aCut.fMinDcaDaughterNegToVert;
  fMaxDcaV0Daughters = aCut.fMaxDcaV0Daughters;
  fMaxDcaV0 = aCut.fMaxDcaV0;
  fMinDcaV0 = aCut.fMinDcaV0;
  fMaxDecayLength = aCut.fMaxDecayLength;
  fMinTransverseDistancePrimSecVtx = aCut.fMinTransverseDistancePrimSecVtx;
  
  fMaxCosPointingAngle = aCut.fMaxCosPointingAngle;
  fMinCosPointingAngle = aCut.fMinCosPointingAngle;
  fParticleType = aCut.fParticleType;
  fEta = aCut.fEta;
  fPtMin = aCut.fPtMin;
  fPtMax = aCut.fPtMax;
  fOnFlyStatus = aCut.fOnFlyStatus;

  fMaxEtaDaughters = aCut.fMaxEtaDaughters;
  fTPCNclsDaughters = aCut.fTPCNclsDaughters;
  fNdofDaughters = aCut.fNdofDaughters;
  fStatusDaughters = aCut.fStatusDaughters;
  fPtMinPosDaughter = aCut.fPtMinPosDaughter;
  fPtMaxPosDaughter = aCut.fPtMaxPosDaughter;
  fPtMinNegDaughter = aCut.fPtMinNegDaughter;
  fPtMaxNegDaughter = aCut.fPtMaxNegDaughter;
  fMinAvgSepDaughters = aCut.fMinAvgSepDaughters;

  fNsigmaPosDaughterTPC = aCut.fNsigmaPosDaughterTPC;
  fNsigmaNegDaughterTPC = aCut.fNsigmaNegDaughterTPC;
  fNsigmaPosDaughterTOF = aCut.fNsigmaPosDaughterTOF;
  fNsigmaNegDaughterTOF = aCut.fNsigmaNegDaughterTOF;

  fBuildPurityAidV0 = aCut.fBuildPurityAidV0;

  fUseLooseInvMassCut = aCut.fUseLooseInvMassCut;
  fLooseInvMassMin = aCut.fLooseInvMassMin;
  fLooseInvMassMax = aCut.fLooseInvMassMax;

  fRadiusV0Min = aCut.fRadiusV0Min;
  fRadiusV0Max = aCut.fRadiusV0Max;

  fRequireTOFPion = aCut.fRequireTOFPion;
  fRequireTOFProton = aCut.fRequireTOFProton;

  if(aCut.fMinvPurityAidHistoV0) fMinvPurityAidHistoV0 = new TH1D(*aCut.fMinvPurityAidHistoV0);
    else fMinvPurityAidHistoV0 = 0;

  fRemoveMisidentified = aCut.fRemoveMisidentified;
  fUseSimpleMisIDCut = aCut.fUseSimpleMisIDCut;
  fBuildMisIDHistograms = aCut.fBuildMisIDHistograms;
  fInvMassRejectK0sMin = aCut.fInvMassRejectK0sMin;
  fInvMassRejectK0sMax = aCut.fInvMassRejectK0sMax;
  fInvMassRejectLambdaMin = aCut.fInvMassRejectLambdaMin;
  fInvMassRejectLambdaMax = aCut.fInvMassRejectLambdaMax;
  fInvMassRejectAntiLambdaMin = aCut.fInvMassRejectAntiLambdaMin;
  fInvMassRejectAntiLambdaMax = aCut.fInvMassRejectAntiLambdaMax;

  if(aCut.fK0sMassOfMisIDV0) fK0sMassOfMisIDV0 = new TH1D(*aCut.fK0sMassOfMisIDV0);
    else fK0sMassOfMisIDV0 = 0;
  if(aCut.fLambdaMassOfMisIDV0) fLambdaMassOfMisIDV0 = new TH1D(*aCut.fLambdaMassOfMisIDV0);
    else fLambdaMassOfMisIDV0 = 0;
  if(aCut.fAntiLambdaMassOfMisIDV0) fAntiLambdaMassOfMisIDV0 = new TH1D(*aCut.fAntiLambdaMassOfMisIDV0);
    else fAntiLambdaMassOfMisIDV0 = 0;

  fIgnoreOnFlyStatus = aCut.fIgnoreOnFlyStatus;

  return *this;
}


//------------------------------
AliFemtoV0TrackCut* AliFemtoV0TrackCut::Clone()
{
  return(new AliFemtoV0TrackCut(*this));
}


//------------------------------
bool AliFemtoV0TrackCut::Pass(const AliFemtoV0* aV0)
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
  } else if (fParticleType == kK0sMC) {
    if (!(aV0->MassK0Short() > fInvMassK0sMin && aV0->MassK0Short() < fInvMassK0sMax) || !(aV0->NegNSigmaTPCP() == 0))
      return false;
    else {
      return true;
    }
  }



  //quality cuts
  if(!fIgnoreOnFlyStatus) {if (aV0->OnFlyStatusV0() != fOnFlyStatus) return false;}
  if (aV0->StatusNeg() == 999 || aV0->StatusPos() == 999) return false;
  if (aV0->TPCNclsPos() < fTPCNclsDaughters) return false;
  if (aV0->TPCNclsNeg() < fTPCNclsDaughters) return false;
  if (aV0->NdofPos() > fNdofDaughters) return false;
  if (aV0->NdofNeg() > fNdofDaughters) return false;
  if (!(aV0->StatusNeg() & fStatusDaughters)) return false;
  if (!(aV0->StatusPos() & fStatusDaughters)) return false;

  //fiducial volume radius
  if(aV0->RadiusV0()<fRadiusV0Min || aV0->RadiusV0()>fRadiusV0Max)
    return false;

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

  //transverse decay vertex
  //TO BE IMPLEMENTED (maybe in future)

  if (fParticleType == kAll)
    return true;


  bool pid_check = false;
  // Looking for lambdas = proton + pim
  if (fParticleType == kLambda) {
    if (IsProtonNSigma(aV0->PtPos(), aV0->PosNSigmaTPCP(), aV0->PosNSigmaTOFP(),fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF,fRequireTOFProton)) //proton
    {
      if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi(),fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF,fRequireTOFPion)) //pion
      {
        pid_check = true;
        if(fUseLooseInvMassCut)
        {
          if(aV0->MassLambda()<fLooseInvMassMin || aV0->MassLambda()>fLooseInvMassMax) return false;
        }
        if(fRemoveMisidentified)
        {
          if(IsMisIDK0s(aV0))
          {
            if(fBuildMisIDHistograms)
            {
              fLambdaMassOfMisIDV0->Fill(aV0->MassLambda());
              fK0sMassOfMisIDV0->Fill(aV0->MassK0Short());
            }
            return false;
          }
        }
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aV0->MassLambda());
        //invariant mass lambda
        if (aV0->MassLambda() < fInvMassLambdaMin || aV0->MassLambda() > fInvMassLambdaMax) return false;    
      }
    }
  }

  //Looking for antilambdas =  antiproton + pip
  else if (fParticleType == kAntiLambda) {
    if (IsProtonNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP(),fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF,fRequireTOFProton)) //proton
    {
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi(),fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF,fRequireTOFPion)) //pion
      {
        pid_check = true;
        if(fUseLooseInvMassCut)
        {
          if(aV0->MassAntiLambda()<fLooseInvMassMin || aV0->MassAntiLambda()>fLooseInvMassMax) return false;
        }
        if(fRemoveMisidentified)
        {
          if(IsMisIDK0s(aV0))
          {
            if(fBuildMisIDHistograms)
            {
              fAntiLambdaMassOfMisIDV0->Fill(aV0->MassAntiLambda());
              fK0sMassOfMisIDV0->Fill(aV0->MassK0Short());
            }
            return false;
          }
        }
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aV0->MassAntiLambda());
        //invariant mass antilambda
        if (aV0->MassAntiLambda() < fInvMassLambdaMin || aV0->MassAntiLambda() > fInvMassLambdaMax) return false;
      }
    }
  }

  //Looking for K0s = pip + pim
  else if (fParticleType == kK0s) {
    if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi(),fNsigmaNegDaughterTPC,fNsigmaNegDaughterTOF,fRequireTOFPion)) //pion-
    {
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi(),fNsigmaPosDaughterTPC,fNsigmaPosDaughterTOF,fRequireTOFPion)) //pion+
      {
        pid_check = true;
        if(fUseLooseInvMassCut)
        {
          if(aV0->MassK0Short()<fLooseInvMassMin || aV0->MassK0Short()>fLooseInvMassMax) return false;
        }
        if(fRemoveMisidentified)
        {
          if(IsMisIDLambda(aV0))
          {
            if(fBuildMisIDHistograms)
            {
              fK0sMassOfMisIDV0->Fill(aV0->MassK0Short());
              fLambdaMassOfMisIDV0->Fill(aV0->MassLambda());
            }
            return false;
          }
          else if(IsMisIDAntiLambda(aV0))
          {
            if(fBuildMisIDHistograms)
            {
              fK0sMassOfMisIDV0->Fill(aV0->MassK0Short());
              fAntiLambdaMassOfMisIDV0->Fill(aV0->MassAntiLambda());
            }
            return false;
          }
        }
        if(fBuildPurityAidV0) fMinvPurityAidHistoV0->Fill(aV0->MassK0Short());
        //invariant mass K0s
        if (aV0->MassK0Short() < fInvMassK0sMin || aV0->MassK0Short() > fInvMassK0sMax) return false;
      }
    }
  }

  if (!pid_check) return false;
  return true;
}
//------------------------------
AliFemtoString AliFemtoV0TrackCut::Report()
{
  // Prepare report from the execution
  string tStemp;
  char tCtemp[100];
  snprintf(tCtemp, 100, "Minimum of Invariant Mass assuming Lambda:\t%lf\n", fInvMassLambdaMin);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Maximum of Invariant Mass assuming Lambda:\t%lf\n", fInvMassLambdaMax);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Minimum DCA of positive daughter to primary vertex:\t%lf\n", fMinDcaDaughterPosToVert);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Minimum DCA of negative daughter to primary vertex:\t%lf\n", fMinDcaDaughterNegToVert);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Max DCA of daughters:\t%lf\n", fMaxDcaV0Daughters);
  tStemp += tCtemp;

  AliFemtoString returnThis = tStemp;
  return returnThis;
}
TList *AliFemtoV0TrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoV0TrackCut.InvMassLambdaMin=%lf", fInvMassLambdaMin);
  tListSetttings->AddLast(new TObjString(buf));
  return tListSetttings;
}

void AliFemtoV0TrackCut::SetMinDaughtersToPrimVertex(double minPos, double minNeg)
{
  fMinDcaDaughterPosToVert = minPos;
  fMinDcaDaughterNegToVert = minNeg;
}

void AliFemtoV0TrackCut::SetMaxDcaV0Daughters(double max)
{
  fMaxDcaV0Daughters = max;
}

void AliFemtoV0TrackCut::SetMaxV0DecayLength(double max)
{
  fMaxDecayLength = max;
}

void AliFemtoV0TrackCut::SetMinTransverseDistancePrimSecVtx(double max)
{
  fMinTransverseDistancePrimSecVtx = max;
}

void AliFemtoV0TrackCut::SetMaxDcaV0(double max)
{
  fMaxDcaV0 = max;
}

void AliFemtoV0TrackCut::SetMinDcaV0(double min)
{
  fMinDcaV0 = min;
}

void AliFemtoV0TrackCut::SetMaxCosPointingAngle(double max) //obsolete
{
  fMaxCosPointingAngle = max;
}

void AliFemtoV0TrackCut::SetMinCosPointingAngle(double min) //correct
{
  fMinCosPointingAngle = min;
}

void AliFemtoV0TrackCut::SetParticleType(short x)
{
  fParticleType = x;
}

void AliFemtoV0TrackCut::SetEta(double x)
{
  fEta = x;
}

void AliFemtoV0TrackCut::SetPt(double min, double max)
{
  fPtMin = min;
  fPtMax = max;
}

void AliFemtoV0TrackCut::SetEtaDaughters(float x)
{
  fMaxEtaDaughters = x;
}
void AliFemtoV0TrackCut::SetTPCnclsDaughters(int x)
{
  fTPCNclsDaughters = x;
}
void AliFemtoV0TrackCut::SetNdofDaughters(int x)
{
  fNdofDaughters = x;
}
void AliFemtoV0TrackCut::SetStatusDaughters(unsigned long x)
{
  fStatusDaughters = x;
}
void AliFemtoV0TrackCut::SetPtPosDaughter(float min, float max)
{
  fPtMinPosDaughter = min;
  fPtMaxPosDaughter = max;
}

void AliFemtoV0TrackCut::SetPtNegDaughter(float min, float max)
{
  fPtMinNegDaughter = min;
  fPtMaxNegDaughter = max;
}

void AliFemtoV0TrackCut::SetOnFlyStatus(bool x)
{
  fOnFlyStatus = x;
}

void AliFemtoV0TrackCut::SetInvariantMassLambda(double min, double max)
{
  fInvMassLambdaMin = min;
  fInvMassLambdaMax = max;
}

void AliFemtoV0TrackCut::SetInvariantMassK0s(double min, double max)
{
  fInvMassK0sMin = min;
  fInvMassK0sMax = max;
}

void AliFemtoV0TrackCut::SetMinAvgSeparation(double minSep)
{
  fMinAvgSepDaughters = minSep;
}

void AliFemtoV0TrackCut::SetNsigmaPosDaughter(double sigma)
{
  fNsigmaPosDaughterTPC = sigma;
  fNsigmaPosDaughterTOF = sigma;
}

void AliFemtoV0TrackCut::SetNsigmaNegDaughter(double sigma)
{
  fNsigmaNegDaughterTPC = sigma;
  fNsigmaNegDaughterTOF = sigma;
}

void AliFemtoV0TrackCut::SetNsigmaPosDaughter(double sigmaTPC, double sigmaTOF)
{
  fNsigmaPosDaughterTPC = sigmaTPC;
  fNsigmaPosDaughterTOF = sigmaTOF;
}

void AliFemtoV0TrackCut::SetNsigmaNegDaughter(double sigmaTPC, double sigmaTOF)
{
  fNsigmaNegDaughterTPC = sigmaTPC;
  fNsigmaNegDaughterTOF = sigmaTOF;
}

void AliFemtoV0TrackCut::SetRadiusV0Min(double radmin)
{
  fRadiusV0Min = radmin;
}

void AliFemtoV0TrackCut::SetRadiusV0Max(double radmax)
{
  fRadiusV0Max = radmax;
}

void AliFemtoV0TrackCut::SetRequireTOFPion(bool req)
{
  fRequireTOFPion = req;
}

void AliFemtoV0TrackCut::SetRequireTOFProton(bool req)
{
  fRequireTOFProton = req;
}


//---------------------PID n Sigma ---------------------------------//
bool AliFemtoV0TrackCut::IsKaonTPCdEdxNSigma(float mom, float nsigmaK)
{
  if (mom < 0.35 && TMath::Abs(nsigmaK) < 5.0) return true;
  if (mom >= 0.35 && mom < 0.5 && TMath::Abs(nsigmaK) < 3.0) return true;
  if (mom >= 0.5 && mom < 0.7 && TMath::Abs(nsigmaK) < 2.0) return true;

  return false;
}


bool AliFemtoV0TrackCut::IsKaonTOFNSigma(float mom, float nsigmaK)
{
  if (mom >= 1.5 && TMath::Abs(nsigmaK) < 2.0) return true;
  return false;
}

bool AliFemtoV0TrackCut::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if (mom < 0.4) {
    if (nsigmaTOFK < -999.) {
      if (TMath::Abs(nsigmaTPCK) < 1.0) return true;
    } else if (TMath::Abs(nsigmaTOFK) < 3.0 && TMath::Abs(nsigmaTPCK) < 3.0) return true;
  } else if (mom >= 0.4 && mom <= 0.6) {
    if (nsigmaTOFK < -999.) {
      if (TMath::Abs(nsigmaTPCK) < 2.0) return true;
    } else if (TMath::Abs(nsigmaTOFK) < 3.0 && TMath::Abs(nsigmaTPCK) < 3.0) return true;
  } else if (nsigmaTOFK < -999.) {
    return false;
  } else if (TMath::Abs(nsigmaTOFK) < 3.0 && TMath::Abs(nsigmaTPCK) < 3.0) return true;

  return false;
}



bool AliFemtoV0TrackCut::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, double nsigmacutTPC, double nsigmacutTOF, bool requireTOF)
{
  
  if (mom < 0.5) {
    if (TMath::Abs(nsigmaTPCPi) < nsigmacutTPC) return true;
  } else {

    if(requireTOF)
      {
    
	if (nsigmaTOFPi < -999.) {
	  if (TMath::Abs(nsigmaTPCPi) < nsigmacutTPC) return true;
	} else {
	  if (TMath::Abs(nsigmaTPCPi) < nsigmacutTPC && TMath::Abs(nsigmaTOFPi) < nsigmacutTOF) return true;
	}

      }
    else
      {
	if (TMath::Abs(nsigmaTPCPi) < nsigmacutTPC) return true;
      }
    
  }

  return false;
}

bool AliFemtoV0TrackCut::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, double nsigmacutTPC, double nsigmacutTOF, bool requireTOF)
{
  if (mom < 0.8) {
    if (TMath::Abs(nsigmaTPCP) < nsigmacutTPC) return true;
  } else {

    if(requireTOF)
      {
    
	if (nsigmaTOFP < -999.) {
	  if (TMath::Abs(nsigmaTPCP) < nsigmacutTPC) return true;
	} else {
	  if (TMath::Abs(nsigmaTPCP) < nsigmacutTPC && TMath::Abs(nsigmaTOFP) < nsigmacutTOF) return true;
	}

      }
    else
      {
	if (TMath::Abs(nsigmaTPCP) < nsigmacutTPC) return true;
      }
    
  }

  return false;
}


void AliFemtoV0TrackCut::SetMinvPurityAidHistoV0(const char* name, const char* title, const int& nbins, const float& aInvMassMin, const float& aInvMassMax)
{
  fBuildPurityAidV0 = true;
  fMinvPurityAidHistoV0 = new TH1D(name,title,nbins,aInvMassMin,aInvMassMax);
  fMinvPurityAidHistoV0->Sumw2();
}

void AliFemtoV0TrackCut::SetLooseInvMassCut(bool aUseCut, double aInvMassMin, double aInvMassMax)
{
  fUseLooseInvMassCut = aUseCut;
  fLooseInvMassMin = aInvMassMin;
  fLooseInvMassMax = aInvMassMax;
}


void AliFemtoV0TrackCut::SetRemoveMisidentified(bool aRemove)
{
  fRemoveMisidentified = aRemove;
}

void AliFemtoV0TrackCut::SetUseSimpleMisIDCut(bool aUse)
{
  fUseSimpleMisIDCut = aUse;
}

void AliFemtoV0TrackCut::SetInvariantMassRejectK0s(double min, double max, bool aRemoveMisidentified)
{
  fInvMassRejectK0sMin = min;
  fInvMassRejectK0sMax = max;

  SetRemoveMisidentified(aRemoveMisidentified);
}

void AliFemtoV0TrackCut::SetInvariantMassRejectLambda(double min, double max, bool aRemoveMisidentified)
{
  fInvMassRejectLambdaMin = min;
  fInvMassRejectLambdaMax = max;

  SetRemoveMisidentified(aRemoveMisidentified);
}

void AliFemtoV0TrackCut::SetInvariantMassRejectAntiLambda(double min, double max, bool aRemoveMisidentified)
{
  fInvMassRejectAntiLambdaMin = min;
  fInvMassRejectAntiLambdaMax = max;

  SetRemoveMisidentified(aRemoveMisidentified);
}

void AliFemtoV0TrackCut::SetInvMassReject(AliFemtoV0Type aV0Type, double aInvMassMin, double aInvMassMax, bool aRemoveMisidentified)
{
  switch(aV0Type)
  {
    case kK0s:
      SetInvariantMassRejectK0s(aInvMassMin,aInvMassMax,aRemoveMisidentified);
      break;

    case kLambda:
      SetInvariantMassRejectLambda(aInvMassMin,aInvMassMax,aRemoveMisidentified);
      break;

    case kAntiLambda:
      SetInvariantMassRejectAntiLambda(aInvMassMin,aInvMassMax,aRemoveMisidentified);
      break;
  }

}

void AliFemtoV0TrackCut::SetBuildMisIDHistograms(bool aBuild)
{
  fBuildMisIDHistograms = aBuild;
  if(fBuildMisIDHistograms) SetDefaultMisIDHistos();
}

void AliFemtoV0TrackCut::SetMisIDHisto(AliFemtoV0Type aMisIDV0Type, const int& nbins, const float& aInvMassMin, const float& aInvMassMax)
{
  TString tTitle1, tTitle2, tTitle3, tTitle;
  tTitle2 = TString("OfMisID");  

  switch(fParticleType)
  {
    case kLambda:
      tTitle3 = TString("Lambda");
      break;
    case kAntiLambda:
      tTitle3 = TString("AntiLambda");
      break;
    case kK0s:
      tTitle3 = TString("K0s");
      break;
  }

  switch(aMisIDV0Type)
  {
    case kLambda:
      tTitle1 = TString("LambdaMass");
      tTitle = tTitle1+tTitle2+tTitle3;
      fLambdaMassOfMisIDV0 = new TH1D(tTitle, "Mass ass. Lambda hyp. of MisID V0", nbins, aInvMassMin, aInvMassMax);
      fLambdaMassOfMisIDV0->Sumw2();
      break;
    case kAntiLambda:
      tTitle1 = TString("AntiLambdaMass");
      tTitle = tTitle1+tTitle2+tTitle3;
      fAntiLambdaMassOfMisIDV0 = new TH1D(tTitle, "Mass ass. AntiLambda hyp. of MisID V0", nbins, aInvMassMin, aInvMassMax);
      fAntiLambdaMassOfMisIDV0->Sumw2();
      break;
    case kK0s:
      tTitle1 = TString("K0sMass");
      tTitle = tTitle1+tTitle2+tTitle3;
      fK0sMassOfMisIDV0 = new TH1D(tTitle, "Mass ass. K0s hyp. of MisID V0", nbins, aInvMassMin, aInvMassMax);
      fK0sMassOfMisIDV0->Sumw2();
      break;
  }

}

void AliFemtoV0TrackCut::SetDefaultMisIDHistos()
{
  double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
  switch(fParticleType)
  {
    case kLambda:
      SetMisIDHisto(kLambda,100,tLambdaMass-0.035,tLambdaMass+0.035);
      SetMisIDHisto(kK0s,100,tK0ShortMass-0.070,tK0ShortMass+0.070);
      break;
    case kAntiLambda:
      SetMisIDHisto(kAntiLambda,100,tLambdaMass-0.035,tLambdaMass+0.035);
      SetMisIDHisto(kK0s,100,tK0ShortMass-0.070,tK0ShortMass+0.070);
      break;
    case kK0s:
      SetMisIDHisto(kK0s,100,tK0ShortMass-0.070,tK0ShortMass+0.070);
      SetMisIDHisto(kLambda,100,tLambdaMass-0.035,tLambdaMass+0.035);
      SetMisIDHisto(kAntiLambda,100,tLambdaMass-0.035,tLambdaMass+0.035);
      break;
  }
}

TObjArray *AliFemtoV0TrackCut::GetMisIDHistos()
{
  TObjArray *tReturnArray = new TObjArray();
  TString tArrayName = "MisID";

  switch(fParticleType)
  {
    case kLambda:
      tArrayName += TString("Lambda");
      tReturnArray->SetName(tArrayName);
      tReturnArray->Add(fLambdaMassOfMisIDV0);
      tReturnArray->Add(fK0sMassOfMisIDV0);
      break;
    case kAntiLambda:
      tArrayName += TString("AntiLambda");
      tReturnArray->SetName(tArrayName);
      tReturnArray->Add(fAntiLambdaMassOfMisIDV0);
      tReturnArray->Add(fK0sMassOfMisIDV0);
      break;
    case kK0s:
      tArrayName += TString("K0s");
      tReturnArray->SetName(tArrayName);
      tReturnArray->Add(fLambdaMassOfMisIDV0);
      tReturnArray->Add(fAntiLambdaMassOfMisIDV0);
      tReturnArray->Add(fK0sMassOfMisIDV0);
      break;
  }

  return tReturnArray;
}

bool AliFemtoV0TrackCut::IsMisIDK0s(const AliFemtoV0* aV0)
{
  if(fUseSimpleMisIDCut)
  {
    if(aV0->MassK0Short()>fInvMassRejectK0sMin && aV0->MassK0Short()<fInvMassRejectK0sMax) return true;
  }

  else
  {
    if(aV0->MassK0Short()>fInvMassRejectK0sMin && aV0->MassK0Short()<fInvMassRejectK0sMax)
    {
      if(IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) //pi+
      {
        if(IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) //pi-
	{
	  return true;
	}
      }
    }

  }
  return false;
}

bool AliFemtoV0TrackCut::IsMisIDLambda(const AliFemtoV0* aV0)
{
  if(fUseSimpleMisIDCut)
  {
    if(aV0->MassLambda()>fInvMassRejectLambdaMin && aV0->MassLambda()<fInvMassRejectLambdaMax) return true;
  }
  else
  {
    if(aV0->MassLambda()>fInvMassRejectLambdaMin && aV0->MassLambda()<fInvMassRejectLambdaMax)
    {
      if(IsProtonNSigma(aV0->PtPos(), aV0->PosNSigmaTPCP(), aV0->PosNSigmaTOFP())) //proton+
      {
        if(IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) //pi-
	{
	  return true;
	}
      }
    }
  }
  return false;
}

bool AliFemtoV0TrackCut::IsMisIDAntiLambda(const AliFemtoV0* aV0)
{
  if(fUseSimpleMisIDCut)
  {
    if(aV0->MassAntiLambda()>fInvMassRejectAntiLambdaMin && aV0->MassAntiLambda()<fInvMassRejectAntiLambdaMax) return true;
  }

  else
  {
    if(aV0->MassAntiLambda()>fInvMassRejectAntiLambdaMin && aV0->MassAntiLambda()<fInvMassRejectAntiLambdaMax)
    {
      if(IsProtonNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP())) //antiproton-
      {
        if(IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) //pi+
	{
	  return true;
	}
      }
    }

  }
  return false;
}

TList *AliFemtoV0TrackCut::GetOutputList()
{
  TList *tOutputList = AliFemtoCutMonitorHandler::GetOutputList();  //add all of the typical objects

  if(fBuildPurityAidV0) tOutputList->Add(fMinvPurityAidHistoV0);
  if(fBuildMisIDHistograms) tOutputList->Add(GetMisIDHistos());

  return tOutputList;
}

