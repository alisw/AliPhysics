#include "AliFemtoXiTrackCut.h"
#include "AliESDtrack.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoXiTrackCut);
  /// \endcond
#endif


//------------------------------
AliFemtoXiTrackCut::AliFemtoXiTrackCut():
  AliFemtoV0TrackCut()
  , fMaxEtaXi(100)
  , fMinPtXi(0)
  , fMaxPtXi(100)
  , fChargeXi(1.0)
  , fMaxEtaBac(100)
  , fMinPtBac(0)
  , fMaxPtBac(100)
  , fTPCNclsBac(0)
  , fNdofBac(100)
  , fStatusBac(1)
  , fMaxDcaXi(1000)
  , fMinDcaXiBac(0)
  , fMaxDcaXiDaughters(1000)
  , fMinCosPointingAngleXi(0)
  , fMinCosPointingAngleV0toXi(0)
  , fMaxDecayLengthXi(100.0)
  , fInvMassXiMin(0)
  , fInvMassXiMax(1000)
  , fInvMassOmegaMin(0)
  , fInvMassOmegaMax(1000) 
  , fInvMassRejectXiMin(0)
  , fInvMassRejectXiMax(1000)
  , fInvMassRejectOmegaMin(0)
  , fInvMassRejectOmegaMax(1000)
  , fParticleTypeXi(kXiMinus)
  , fRadiusXiMin(0.)
  , fRadiusXiMax(99999.0)
  , fBuildPurityAidXi(false)
  , fMinvPurityAidHistoXi(nullptr)
  , fSidebandAnalysis(false)
  , fInvMassRange1XiMin(0)
  , fInvMassRange1XiMax(1000)
  , fInvMassRange2XiMin(0)
  , fInvMassRange2XiMax(1000)
  , fInvMassRange1OmegaMin(0)
  , fInvMassRange1OmegaMax(1000)
  , fInvMassRange2OmegaMin(0)
  , fInvMassRange2OmegaMax(1000)
{
  // Default constructor
}

//------------------------------
AliFemtoXiTrackCut::~AliFemtoXiTrackCut()
{
  /* noop */
}

//------------------------------
AliFemtoXiTrackCut::AliFemtoXiTrackCut(const AliFemtoXiTrackCut& aCut) :
  AliFemtoV0TrackCut(aCut)
  , fMaxEtaXi(aCut.fMaxEtaXi)
  , fMinPtXi(aCut.fMinPtXi)
  , fMaxPtXi(aCut.fMaxPtXi)
  , fChargeXi(aCut.fChargeXi)
  , fMaxEtaBac(aCut.fMaxEtaBac)
  , fMinPtBac(aCut.fMinPtBac)
  , fMaxPtBac(aCut.fMaxPtBac)
  , fTPCNclsBac(aCut.fTPCNclsBac)
  , fNdofBac(aCut.fNdofBac)
  , fStatusBac(aCut.fStatusBac)
  , fMaxDcaXi(aCut.fMaxDcaXi)
  , fMinDcaXiBac(aCut.fMinDcaXiBac)
  , fMaxDcaXiDaughters(aCut.fMaxDcaXiDaughters)
  , fMinCosPointingAngleXi(aCut.fMinCosPointingAngleXi)
  , fMinCosPointingAngleV0toXi(aCut.fMinCosPointingAngleV0toXi)
  , fMaxDecayLengthXi(aCut.fMaxDecayLengthXi)
  , fInvMassXiMin(aCut.fInvMassXiMin)
  , fInvMassXiMax(aCut.fInvMassXiMax)
  , fInvMassOmegaMin(aCut.fInvMassOmegaMin)
  , fInvMassOmegaMax(aCut.fInvMassOmegaMax)
  , fInvMassRejectXiMin(aCut.fInvMassRejectXiMin)
  , fInvMassRejectXiMax(aCut.fInvMassRejectXiMax)
  , fInvMassRejectOmegaMin(aCut.fInvMassRejectOmegaMin)
  , fInvMassRejectOmegaMax(aCut.fInvMassRejectOmegaMax)
  , fParticleTypeXi(aCut.fParticleTypeXi)
  , fRadiusXiMin(aCut.fRadiusXiMin)
  , fRadiusXiMax(aCut.fRadiusXiMax)
  , fBuildPurityAidXi(aCut.fBuildPurityAidXi)
  , fMinvPurityAidHistoXi(nullptr)
  , fSidebandAnalysis(aCut.fSidebandAnalysis)
  , fInvMassRange1XiMin(aCut.fInvMassRange1XiMin)
  , fInvMassRange1XiMax(aCut.fInvMassRange1XiMax)
  , fInvMassRange2XiMin(aCut.fInvMassRange2XiMin)
  , fInvMassRange2XiMax(aCut.fInvMassRange2XiMax)
  , fInvMassRange1OmegaMin(aCut.fInvMassRange1OmegaMin)
  , fInvMassRange1OmegaMax(aCut.fInvMassRange1OmegaMax)
  , fInvMassRange2OmegaMin(aCut.fInvMassRange2OmegaMin)
  , fInvMassRange2OmegaMax(aCut.fInvMassRange2OmegaMax)
{
  //copy constructor
  if (aCut.fMinvPurityAidHistoXi) fMinvPurityAidHistoXi = new TH1D(*aCut.fMinvPurityAidHistoXi);
}

//------------------------------
AliFemtoXiTrackCut& AliFemtoXiTrackCut::operator=(const AliFemtoXiTrackCut& aCut)
{
  //assignment operator
  if (this == &aCut) return *this;

  AliFemtoV0TrackCut::operator=(aCut);

  fMaxEtaXi = aCut.fMaxEtaXi;
  fMinPtXi = aCut.fMinPtXi;
  fMaxPtXi = aCut.fMaxPtXi;
  fChargeXi = aCut.fChargeXi;
  fMaxEtaBac = aCut.fMaxEtaBac;
  fMinPtBac = aCut.fMinPtBac;
  fMaxPtBac = aCut.fMaxPtBac;
  fTPCNclsBac = aCut.fTPCNclsBac;
  fNdofBac = aCut.fNdofBac;
  fStatusBac = aCut.fStatusBac;
  fMaxDcaXi = aCut.fMaxDcaXi;
  fMinDcaXiBac = aCut.fMinDcaXiBac;
  fMaxDcaXiDaughters = aCut.fMaxDcaXiDaughters;
  fMinCosPointingAngleXi = aCut.fMinCosPointingAngleXi;
  fMinCosPointingAngleV0toXi = aCut.fMinCosPointingAngleV0toXi;
  fMaxDecayLengthXi = aCut.fMaxDecayLengthXi;
  fInvMassXiMin = aCut.fInvMassXiMin;
  fInvMassXiMax = aCut.fInvMassXiMax;
  fInvMassOmegaMin = aCut.fInvMassOmegaMin;
  fInvMassOmegaMax = aCut.fInvMassOmegaMax;
  fInvMassRejectXiMin = aCut.fInvMassRejectXiMin;
  fInvMassRejectXiMax = aCut.fInvMassRejectXiMax;
  fInvMassRejectOmegaMin = aCut.fInvMassRejectOmegaMin;
  fInvMassRejectOmegaMax = aCut.fInvMassRejectOmegaMax;
  fParticleTypeXi = aCut.fParticleTypeXi;
  fRadiusXiMin = aCut.fRadiusXiMin;
  fRadiusXiMax = aCut.fRadiusXiMax;
  fBuildPurityAidXi = aCut.fBuildPurityAidXi;
  fSidebandAnalysis = aCut.fSidebandAnalysis;
  fInvMassRange1XiMin = aCut.fInvMassRange1XiMin;
  fInvMassRange1XiMax = aCut.fInvMassRange1XiMax;
  fInvMassRange2XiMin = aCut.fInvMassRange2XiMin;
  fInvMassRange2XiMax = aCut.fInvMassRange2XiMax;
  fInvMassRange1OmegaMin = aCut.fInvMassRange1OmegaMin;
  fInvMassRange1OmegaMax = aCut.fInvMassRange1OmegaMax;
  fInvMassRange2OmegaMin = aCut.fInvMassRange2OmegaMin;
  fInvMassRange2OmegaMax = aCut.fInvMassRange2OmegaMax;

  if (aCut.fMinvPurityAidHistoXi != nullptr) {
    if (fMinvPurityAidHistoXi == nullptr) {
      fMinvPurityAidHistoXi = new TH1D(*aCut.fMinvPurityAidHistoXi);
    } else {
      *fMinvPurityAidHistoXi = *aCut.fMinvPurityAidHistoXi;
    }
  } else if (fMinvPurityAidHistoXi) {
    delete fMinvPurityAidHistoXi;
    fMinvPurityAidHistoXi = nullptr;
  }

  return *this;
}

//------------------------------
AliFemtoXiTrackCut* AliFemtoXiTrackCut::Clone()
{
  return(new AliFemtoXiTrackCut(*this));
}



bool AliFemtoXiTrackCut::Pass(const AliFemtoXi* aXi)
{
  // test the particle and return
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria
  Float_t pt = aXi->PtXi();
  Float_t eta = aXi->EtaXi();
  if(aXi->ChargeXi()==0)
    return false;


  //ParticleType selection
  //If fParticleTypeXi=kAll, any charged candidate will pass
  if( ((fParticleTypeXi == kXiPlus) || (fParticleTypeXi == kOmegaPlus)) && aXi->ChargeXi() == -1)
    return false;

  if( ((fParticleTypeXi == kXiMinus) || (fParticleTypeXi == kOmegaMinus)) && aXi->ChargeXi() == 1)
    return false;

  //kinematic cuts
  if(TMath::Abs(eta) > fMaxEtaXi) return false;  //put in kinematic cuts by hand
  if(pt < fMinPtXi) return false;
  if(pt > fMaxPtXi) return false;
  if(TMath::Abs(aXi->EtaBac()) > fMaxEtaBac) return false;
  if(!fNanoAODAnalysis){
    if(aXi->PtBac()< fMinPtBac) return false;
    if(aXi->PtBac()> fMaxPtBac) return false;
  }


    //Xi from kinematics information
    if (fParticleTypeXi == kXiMinusMC || fParticleTypeXi == kXiPlusMC) {
      if(!(aXi->MassXi()>fInvMassXiMin && aXi->MassXi()<fInvMassXiMax) || !(aXi->BacNSigmaTPCPi()==0))
	return false;
      else
	{
	  return true;
	}
    }

    //Omega from kinematics information
    if (fParticleTypeXi == kOmegaMinusMC || fParticleTypeXi == kOmegaPlusMC) {
      if(!(aXi->MassOmega()>fInvMassOmegaMin && aXi->MassOmega()<fInvMassOmegaMax) || !(aXi->BacNSigmaTPCK()==0))
	return false;
      else
	{
	  return true;
	}
    }
    
    //quality cuts

    if(aXi->TPCNclsBac()<fTPCNclsBac) return false;
    if(aXi->NdofBac()>fNdofBac) return false;
    if(!fNanoAODAnalysis){
      if(!(aXi->StatusBac()&fStatusBac)) return false;
      if(aXi->StatusBac() == 999) return false;
    }
    //DCA Xi to prim vertex
    if(TMath::Abs(aXi->DcaXiToPrimVertex())>fMaxDcaXi)
      return false;

    //DCA Xi bachelor to prim vertex
    if(TMath::Abs(aXi->DcaBacToPrimVertex())<fMinDcaXiBac)
      return false;

    //DCA Xi daughters
    if(TMath::Abs(aXi->DcaXiDaughters())>fMaxDcaXiDaughters)
      return false;

    //cos pointing angle
    if(aXi->CosPointingAngleXi()<fMinCosPointingAngleXi)
      return false;

    //cos pointing angle of V0 to Xi
    if(aXi->CosPointingAngleV0toXi()<fMinCosPointingAngleV0toXi)
      return false;

    //decay length
    if(aXi->DecayLengthXi()>fMaxDecayLengthXi)
      return false;

    //fiducial volume radius
    if(aXi->RadiusXi()<fRadiusXiMin || aXi->RadiusXi()>fRadiusXiMax)
      return false;

  if(fParticleTypeXi == kAll)
    return true;

  bool pid_check=false;
  // Looking for Xi
  if (fParticleTypeXi == kXiMinus || fParticleTypeXi == kXiPlus) {
    if (IsPionNSigmaBac(aXi->PtBac(), aXi->BacNSigmaTPCPi(), aXi->BacNSigmaTOFPi())) //pion
	{
	  pid_check=true;
	}
  }
  else if (fParticleTypeXi == kOmegaMinus || fParticleTypeXi == kOmegaPlus) {   //looking for Omega
    if (IsKaonNSigmaBac(aXi->PtBac(), aXi->BacNSigmaTPCK(), aXi->BacNSigmaTOFK())) //kaon
	{
	  pid_check=true;
	}
  }


  if (!pid_check) return false;

  if(!AliFemtoV0TrackCut::Pass(aXi))
    return false;


  if(fBuildPurityAidXi) {fMinvPurityAidHistoXi->Fill(aXi->MassXi());}

  //invariant mass Xi
  if((fParticleTypeXi == kXiPlus || fParticleTypeXi == kXiMinus) && !fSidebandAnalysis)
    {
      if(aXi->MassXi()<fInvMassXiMin || aXi->MassXi()>fInvMassXiMax)
	{
	  return false;
	}
    }
  else if((fParticleTypeXi == kXiPlus || fParticleTypeXi == kXiMinus) && fSidebandAnalysis)
    {
      if( (aXi->MassXi()<fInvMassRange1XiMin || aXi->MassXi()>fInvMassRange2XiMax) || (aXi->MassXi()>fInvMassRange1XiMax && aXi->MassXi()<fInvMassRange2XiMin) )
	{
	  return false;
	}
    }
  
  //invariant mass Omega
  if((fParticleTypeXi == kOmegaPlus || fParticleTypeXi == kOmegaMinus) && !fSidebandAnalysis)
    {
      if(aXi->MassOmega()<fInvMassOmegaMin || aXi->MassOmega()>fInvMassOmegaMax)
	{
	  return false;
	}
    }
  else if((fParticleTypeXi == kOmegaPlus || fParticleTypeXi == kOmegaMinus) && fSidebandAnalysis)
    {
      if( (aXi->MassOmega()<fInvMassRange1OmegaMin || aXi->MassOmega()>fInvMassRange2OmegaMax) || (aXi->MassOmega()>fInvMassRange1OmegaMax && aXi->MassOmega()<fInvMassRange2OmegaMin) )
	{
	  return false;
	}
    }

  //removing particles in the given Minv window (to reject omegas in Xi sample)
  if(fParticleTypeXi == kXiPlus || fParticleTypeXi == kXiMinus)
    {
      if(aXi->MassOmega()>fInvMassRejectOmegaMin && aXi->MassOmega()<fInvMassRejectOmegaMax)
	{
	  return false;
	}
    }


  //removing particles in the given Minv window (to reject Xis in Omega sample)
  if(fParticleTypeXi == kOmegaPlus || fParticleTypeXi == kOmegaMinus)
    {
      if(aXi->MassXi()>fInvMassRejectXiMin && aXi->MassXi()<fInvMassRejectXiMax)
	{
	  return false;
	}
    }
  

  return true;
}
//------------------------------
AliFemtoString AliFemtoXiTrackCut::Report()
{
  // Prepare report from the execution
  TString report;
  return AliFemtoString(report.Data());
}

TList *AliFemtoXiTrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();

  return tListSetttings;
}


bool AliFemtoXiTrackCut::IsPionNSigmaBac(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if(TMath::Abs(nsigmaTPCPi)<3.0) return true;


  return false;
}

bool AliFemtoXiTrackCut::IsKaonNSigmaBac(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if(TMath::Abs(nsigmaTPCK)<3.0) return true;


  return false;
}



void AliFemtoXiTrackCut::SetEtaXi(double x){
  fMaxEtaXi = x;
}

void AliFemtoXiTrackCut::SetPtXi(double min, double max){
  fMinPtXi = min;
  fMaxPtXi = max;
}

void AliFemtoXiTrackCut::SetChargeXi(int x){
  fChargeXi = x;
}

void AliFemtoXiTrackCut::SetMaxDecayLengthXi(double x){
  fMaxDecayLengthXi = x;
}

void AliFemtoXiTrackCut::SetEtaBac(double x){
  fMaxEtaBac = x;
}

void AliFemtoXiTrackCut::SetPtBac(double min, double max){
  fMinPtBac = min;
  fMaxPtBac = max;
}


void AliFemtoXiTrackCut::SetTPCnclsBac(int x){
  fTPCNclsBac = x;
}

void AliFemtoXiTrackCut::SetNdofBac(double x){
  fNdofBac = x;
}


void AliFemtoXiTrackCut::SetStatusBac(unsigned long x) {
  fStatusBac = x;
}

void AliFemtoXiTrackCut::SetInvariantMassXi(double min, double max)
{
  fInvMassXiMin = min;
  fInvMassXiMax = max;

}

void AliFemtoXiTrackCut::SetInvariantMassOmega(double min, double max)
{
  fInvMassOmegaMin = min;
  fInvMassOmegaMax = max;

}

void AliFemtoXiTrackCut::SetInvariantMassRejectXi(double min, double max)
{
  fInvMassRejectXiMin = min;
  fInvMassRejectXiMax = max;

}

void AliFemtoXiTrackCut::SetInvariantMassRejectOmega(double min, double max)
{
  fInvMassRejectOmegaMin = min;
  fInvMassRejectOmegaMax = max;

}

void AliFemtoXiTrackCut::SetMaxDcaXi(double x)
{
  fMaxDcaXi = x;
}

void AliFemtoXiTrackCut::SetMinDcaXiBac(double x)
{
  fMinDcaXiBac = x;
}

void AliFemtoXiTrackCut::SetMaxDcaXiDaughters(double x)
{
  fMaxDcaXiDaughters = x;
}

void AliFemtoXiTrackCut::SetMinCosPointingAngleXi(double min)
{
  fMinCosPointingAngleXi = min;
}

void AliFemtoXiTrackCut::SetMinCosPointingAngleV0toXi(double min)
{
  fMinCosPointingAngleV0toXi = min;
}


void AliFemtoXiTrackCut::SetParticleTypeXi(short x)
{
  fParticleTypeXi = x;
}

void AliFemtoXiTrackCut::SetMinvPurityAidHistoXi(const char* name, const char* title, const int& nbins, const float& aInvMassMin, const float& aInvMassMax)
{
  fBuildPurityAidXi = true;
  fMinvPurityAidHistoXi = new TH1D(name,title,nbins,aInvMassMin,aInvMassMax);
  fMinvPurityAidHistoXi->Sumw2();
}


void AliFemtoXiTrackCut::SetSidebandAnalysis(bool sideband)
{
  fSidebandAnalysis = sideband;
}


void AliFemtoXiTrackCut::SetInvariantMassXiSideband(double min1, double max1, double min2, double max2)
{
  fInvMassRange1XiMin = min1;
  fInvMassRange1XiMax = max1;
  fInvMassRange2XiMin = min2;
  fInvMassRange2XiMax = max2;
}

void AliFemtoXiTrackCut::SetInvariantMassOmegaSideband(double min1, double max1, double min2, double max2)
{
  fInvMassRange1OmegaMin = min1;
  fInvMassRange1OmegaMax = max1;
  fInvMassRange2OmegaMin = min2;
  fInvMassRange2OmegaMax = max2;
}


TList *AliFemtoXiTrackCut::GetOutputList()
{
  TList *tOutputList = AliFemtoCutMonitorHandler::GetOutputList();  //add all of the typical objects

  if(fBuildPurityAidXi) tOutputList->Add(fMinvPurityAidHistoXi);
  if(fBuildPurityAidV0) tOutputList->Add(fMinvPurityAidHistoV0);

  return tOutputList;
}
