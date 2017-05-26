#include "AliFemtoXiTrackCut.h"
#include "AliESDtrack.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoXiTrackCut);
  /// \endcond
#endif


//------------------------------
AliFemtoXiTrackCut::AliFemtoXiTrackCut() : AliFemtoV0TrackCut(), fMaxEtaXi(100), fMinPtXi(0), fMaxPtXi(100), fChargeXi(1.0), fMaxEtaBac(100), fMinPtBac(0), fMaxPtBac(100), fTPCNclsBac(0), fNdofBac(100), fStatusBac(0), fMaxDcaXi(0), fMinDcaXiBac(0), fMaxDcaXiDaughters(0), fMinCosPointingAngleXi(0), fMinCosPointingAngleV0toXi(0), fMaxDecayLengthXi(100.0), fInvMassXiMin(0), fInvMassXiMax(1000), fParticleTypeXi(kXiMinus), fRadiusXiMin(0.), fRadiusXiMax(99999.0), fBuildPurityAidXi(false), fMinvPurityAidHistoXi(0)
{
  // Default constructor
}

//------------------------------
AliFemtoXiTrackCut::~AliFemtoXiTrackCut(){
  /* noop */
}

//------------------------------
AliFemtoXiTrackCut::AliFemtoXiTrackCut(const AliFemtoXiTrackCut& aCut) :
  AliFemtoV0TrackCut(aCut),
  fMaxEtaXi(aCut.fMaxEtaXi),fMinPtXi(aCut.fMinPtXi),fMaxPtXi(aCut.fMaxPtXi),fChargeXi(aCut.fChargeXi),
  fMaxEtaBac(aCut.fMaxEtaBac),fMinPtBac(aCut.fMinPtBac),fMaxPtBac(aCut.fMaxPtBac),fTPCNclsBac(aCut.fTPCNclsBac),
  fNdofBac(aCut.fNdofBac),fStatusBac(aCut.fStatusBac),fMaxDcaXi(aCut.fMaxDcaXi),fMinDcaXiBac(aCut.fMinDcaXiBac),
  fMaxDcaXiDaughters(aCut.fMaxDcaXiDaughters),fMinCosPointingAngleXi(aCut.fMinCosPointingAngleXi), fMinCosPointingAngleV0toXi(aCut.fMinCosPointingAngleV0toXi),
  fMaxDecayLengthXi(aCut.fMaxDecayLengthXi),fInvMassXiMin(aCut.fInvMassXiMin),fInvMassXiMax(aCut.fInvMassXiMax),
  fParticleTypeXi(aCut.fParticleTypeXi), fRadiusXiMin(aCut.fRadiusXiMin), fRadiusXiMax(aCut.fRadiusXiMax), fBuildPurityAidXi(aCut.fBuildPurityAidXi)
{
  //copy constructor
  if(aCut.fMinvPurityAidHistoXi) fMinvPurityAidHistoXi = new TH1D(*aCut.fMinvPurityAidHistoXi);
  else fMinvPurityAidHistoXi = 0;
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
  fParticleTypeXi = aCut.fParticleTypeXi;
  fRadiusXiMin = aCut.fRadiusXiMin;
  fRadiusXiMax = aCut.fRadiusXiMax;
  fBuildPurityAidXi = aCut.fBuildPurityAidXi;

  if(aCut.fMinvPurityAidHistoXi) fMinvPurityAidHistoXi = new TH1D(*aCut.fMinvPurityAidHistoXi);
  else fMinvPurityAidHistoXi = 0;

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
  if(fParticleTypeXi == kXiPlus && aXi->ChargeXi() == -1) 
    return false;

  if(fParticleTypeXi == kXiMinus && aXi->ChargeXi() == 1) 
    return false;

  //kinematic cuts
  if(TMath::Abs(eta) > fMaxEtaXi) return false;  //put in kinematic cuts by hand
  if(pt < fMinPtXi) return false;
  if(pt > fMaxPtXi) return false;
  if(TMath::Abs(aXi->EtaBac()) > fMaxEtaBac) return false;
  if(aXi->PtBac()< fMinPtBac) return false; 
  if(aXi->PtBac()> fMaxPtBac) return false; 

  

    //Xi from kinematics information
    if (fParticleTypeXi == kXiMinusMC || fParticleTypeXi == kXiPlusMC) {
      if(!(aXi->MassXi()>fInvMassXiMin && aXi->MassXi()<fInvMassXiMax) || !(aXi->BacNSigmaTPCPi()==0))
	return false; 
      else
	{
	  return true;  
	} 
    }


    //quality cuts
    if(aXi->StatusBac() == 999) return false;
    if(aXi->TPCNclsBac()<fTPCNclsBac) return false;
    if(aXi->NdofBac()>fNdofBac) return false;
    if(!(aXi->StatusBac()&fStatusBac)) return false;


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

  if (!pid_check) return false;
  
  if(!AliFemtoV0TrackCut::Pass(aXi))
    return false;

  if(fBuildPurityAidXi) {fMinvPurityAidHistoXi->Fill(aXi->MassXi());}

   //invariant mass Xi
  if(aXi->MassXi()<fInvMassXiMin || aXi->MassXi()>fInvMassXiMax)
     {
       return false;
     }
  
  
  return true;
}
//------------------------------
AliFemtoString AliFemtoXiTrackCut::Report()
{
  // Prepare report from the execution
  string tStemp;
  char tCtemp[100];


  AliFemtoString returnThis = tStemp;
  return returnThis;
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

TList *AliFemtoXiTrackCut::GetOutputList()
{
  TList *tOutputList = AliFemtoCutMonitorHandler::GetOutputList();  //add all of the typical objects

  if(fBuildPurityAidXi) tOutputList->Add(fMinvPurityAidHistoXi);
  if(fBuildPurityAidV0) tOutputList->Add(fMinvPurityAidHistoV0);

  return tOutputList;
}

