#include "AliFemtoV0TrackCut.h"
#include "AliESDtrack.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoV0TrackCut)
#endif



AliFemtoV0TrackCut::AliFemtoV0TrackCut() :
fInvMassLambdaMin(0),fInvMassLambdaMax(99),fMinDcaDaughtersToVert(0),fMaxDcaV0Daughters(99),fMaxDcaV0(99), fMaxCosPointingAngle(0), fParticleType(99), fEta(0.8), fPtMin(0), fPtMax(100), fOnFlyStatus(kFALSE), fMaxEtaDaughters(100), fTPCNclsDaughters(0), fNdofDaughters(10), fStatusDaughters(0), fPtMinDaughters(0), fPtMaxDaughters(99)
{
  // Default constructor
 }
 //------------------------------
AliFemtoV0TrackCut::~AliFemtoV0TrackCut(){
  /* noop */
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
    if(TMath::Abs(eta) > fEta) return false;  //put in kinematic cuts by hand
    if(pt < fPtMin) return false;
    if(pt > fPtMax) return false;
    if(TMath::Abs(aV0->EtaPos()) > fMaxEtaDaughters) return false;
    if(TMath::Abs(aV0->EtaNeg()) > fMaxEtaDaughters) return false;

    if(aV0->PtPos()< fPtMinDaughters) return false;
    if(aV0->PtNeg()< fPtMinDaughters) return false;
    if(aV0->PtPos()> fPtMaxDaughters) return false;
    if(aV0->PtNeg()> fPtMaxDaughters) return false;

    //V0 from kinematics information
    if (fParticleType == kLambdaMC ) {
      if(!(aV0->MassLambda()>fInvMassLambdaMin && aV0->MassLambda()<fInvMassLambdaMax) || !(aV0->PosNSigmaTPCP()==0))
	return false; 
      else
	{
	  return true;  
	} 
    }
    else if (fParticleType == kAntiLambdaMC) {
      if(!(aV0->MassLambda()>fInvMassLambdaMin &&aV0->MassLambda()<fInvMassLambdaMax) || !(aV0->NegNSigmaTPCP()==0))
	return false; 
      else
	{
	  return true;   
	}
    }

    //quality cuts
    if(fOnFlyStatus) if(aV0->OnFlyStatusV0()) return false;
    if(aV0->StatusNeg() == 999 || aV0->StatusPos() == 999) return false;
    if(aV0->TPCNclsPos()<fTPCNclsDaughters) return false;
    if(aV0->TPCNclsNeg()<fTPCNclsDaughters) return false;
    if(aV0->NdofPos()>fNdofDaughters) return false;
    if(aV0->NdofNeg()>fNdofDaughters) return false;
    if(!(aV0->StatusNeg()&fStatusDaughters)) return false;
    if(!(aV0->StatusPos()&fStatusDaughters)) return false;
    
  //DCA between daughter particles
    if(TMath::Abs(aV0->DcaV0Daughters())>fMaxDcaV0Daughters)
    return false;

   //DCA of daughters to primary vertex
    if(TMath::Abs(aV0->DcaPosToPrimVertex())<fMinDcaDaughtersToVert || TMath::Abs(aV0->DcaNegToPrimVertex())<fMinDcaDaughtersToVert)
    return false;

  //DCA V0 to prim vertex
    if(TMath::Abs(aV0->DcaV0ToPrimVertex())>fMaxDcaV0)
    return false;

  //cos pointing angle
  if(aV0->CosPointingAngle()<fMaxCosPointingAngle)
    return false;

  if(fParticleType == kAll)
    return true;

  bool pid_check=false;
  // Looking for lambdas = proton + pim
  if (fParticleType == 0 ) {
    if (IsProtonNSigma(aV0->PtotPos(), aV0->PosNSigmaTPCP(), aV0->PosNSigmaTOFP())) //proton
      if (IsPionNSigma(aV0->PtotNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) //pion
	{
	  pid_check=true;
	  //invariant mass lambda
	  if(aV0->MassLambda()<fInvMassLambdaMin || aV0->MassLambda()>fInvMassLambdaMax)
	    return false;
	}

  }//Looking for antilambdas =  antiproton + pip
  else if (fParticleType == 1) {
    if (IsProtonNSigma(aV0->PtotNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP())) //proton
      if (IsPionNSigma(aV0->PtotPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) //pion
	{
	  pid_check=true;
	  //invariant mass antilambda
	  if(aV0->MassAntiLambda()<fInvMassLambdaMin || aV0->MassAntiLambda()>fInvMassLambdaMax)
	    return false;
	}
  }


  if(!pid_check) return false;

  

  return true;
    
    
}
//------------------------------
AliFemtoString AliFemtoV0TrackCut::Report()
{
  // Prepare report from the execution
  string tStemp;
  char tCtemp[100];
  snprintf(tCtemp , 100, "Minimum of Invariant Mass assuming Lambda:\t%lf\n",fInvMassLambdaMin);
  tStemp+=tCtemp;
  snprintf(tCtemp , 100, "Maximum of Invariant Mass assuming Lambda:\t%lf\n",fInvMassLambdaMax);
  tStemp+=tCtemp;
  snprintf(tCtemp , 100, "Minimum DCA of daughters to primary vertex:\t%lf\n",fMinDcaDaughtersToVert);
  tStemp+=tCtemp;
  snprintf(tCtemp , 100, "Max DCA of daughters:\t%lf\n",fMaxDcaV0Daughters);
  tStemp+=tCtemp;



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

void AliFemtoV0TrackCut::SetMinDaughtersToPrimVertex(double min)
{
  fMinDcaDaughtersToVert=min;
}

void AliFemtoV0TrackCut::SetMaxDcaV0Daughters(double max)
{
  fMaxDcaV0Daughters=max;
};

void AliFemtoV0TrackCut::SetMaxDcaV0(double max)
{
  fMaxDcaV0=max;
};

void AliFemtoV0TrackCut::SetMaxCosPointingAngle(double max)
{
  fMaxCosPointingAngle = max;
}

void AliFemtoV0TrackCut::SetParticleType(short x)
{
  fParticleType = x;
}

void AliFemtoV0TrackCut::SetEta(double x){
  fEta = x;
}

void AliFemtoV0TrackCut::SetPt(double min, double max){
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
  fNdofDaughters=x;
}
void AliFemtoV0TrackCut::SetStatusDaughters(unsigned long x)
{
  fStatusDaughters=x;
}
void AliFemtoV0TrackCut::SetPtDaughters(float min,float max)
{
  fPtMinDaughters = min;
  fPtMaxDaughters = max;
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

//---------------------PID n Sigma ---------------------------------//
bool AliFemtoV0TrackCut::IsKaonTPCdEdxNSigma(float mom, float nsigmaK)
{
  cout<<" AliFemtoV0TrackCut::IsKaonTPCdEdxNSigma "<<mom<<" "<<nsigmaK<<endl;


  if(mom<0.35 && TMath::Abs(nsigmaK)<5.0)return true;
  if(mom>=0.35 && mom<0.5 && TMath::Abs(nsigmaK)<3.0)return true; 
  if(mom>=0.5 && mom<0.7 && TMath::Abs(nsigmaK)<2.0)return true;

  return false;
}


bool AliFemtoV0TrackCut::IsKaonTOFNSigma(float mom, float nsigmaK)
{
  cout<<" AliFemtoV0TrackCut::IsKaonTPCdEdxNSigma "<<mom<<" "<<nsigmaK<<endl;
  if(mom>=1.5 && TMath::Abs(nsigmaK)<2.0)return true; 
  return false;
}

bool AliFemtoV0TrackCut::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{

      if(TMath::Abs(nsigmaTOFK)<3.0 && mom<=1.5 && TMath::Abs(nsigmaTPCK)<3.0)return true;
      if(TMath::Abs(nsigmaTOFK)<2.0 && mom>1.5 && TMath::Abs(nsigmaTPCK)<3.0)return true;

      //no TOF signal
 
      if(nsigmaTOFK<=-9999){
	//cout <<"/////////////// AliFemtoV0TrackCut::IsKaonNSigma  NO TOF SIGNAL ////////////" <<endl;
	if(TMath::Abs(nsigmaTPCK)<2.0) return true;
	/*if(mom<0.4 && TMath::Abs(nsigmaTPCK)<1.0)return true;
	if(mom>=0.4 && mom<0.5 && TMath::Abs(nsigmaTPCK)<2.0)return true;
	if(mom>=0.5 && mom<0.6 && TMath::Abs(nsigmaTPCK)<2.0)return true;*/
      }
 
  return false;
}



bool AliFemtoV0TrackCut::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  mom = 1; //because of warning in the compilation
  nsigmaTOFPi = 1;
  // cout<<" AliFemtoV0TrackCut::IsKaonNSigma "<<mom<<" tpc "<<nsigmaTPCK<<" tof "<<nsigmaTOFK<<endl;

  //TOF signal
  // if(TMath::Abs(nsigmaTOFPi)<2.0 && mom<=1.5 && TMath::Abs(nsigmaTPCPi)<2.0)return true;
  //  if(TMath::Abs(nsigmaTOFPi)<2.0 && mom>1.5 && TMath::Abs(nsigmaTPCPi)<2.0)return true;


      //no TOF signal
  //  if(nsigmaTOFPi<=-9999){
  //	if(TMath::Abs(nsigmaTPCPi)<2.0) return true;
	/*if(mom<0.35 && TMath::Abs(nsigmaTPCPi)<2.0)return true;
	if(mom>=0.35 && mom<0.5 && TMath::Abs(nsigmaTPCPi)<2.0)return true;
	if(mom>=0.5 && TMath::Abs(nsigmaTPCPi)<2.0)return true;*/
  //      }

  if( TMath::Abs(nsigmaTPCPi)<3.0) return true;

  return false;
}


bool AliFemtoV0TrackCut::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
  mom = 1; //because of warning in the compilation
  nsigmaTOFP = 1;
  // cout<<" AliFemtoV0TrackCut::IsKaonNSigma "<<mom<<" tpc "<<nsigmaTPCK<<" tof "<<nsigmaTOFK<<endl;

      //TOF signal
  //  if(TMath::Abs(nsigmaTOFP)<2.0 && mom<=1.5 && TMath::Abs(nsigmaTPCP)<2.0)return true;
  //  if(TMath::Abs(nsigmaTOFP)<2.0 && mom>1.5 && TMath::Abs(nsigmaTPCP)<2.0)return true;

      //no TOF signal
  //  if(nsigmaTOFP<=-9999){
  //	if(TMath::Abs(nsigmaTPCP)<2.0) return true;
	/*if(mom<0.5 && TMath::Abs(nsigmaTPCP)<2.0)return true;
	  if(mom>=0.5 && mom<0.8 && TMath::Abs(nsigmaTPCP)<2.0)return true;*/
  //  }

  if( TMath::Abs(nsigmaTPCP)<3.0) return true;

  return false;
}
