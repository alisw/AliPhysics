#include "AliFemtoWRzTrackCut.h"


AliFemtoWRzTrackCut::AliFemtoWRzTrackCut():
  AliFemtoESDTrackCut()
  {
     fdEdxcut = true ;
     fNsigmaRejection = 3.0;
     fNSigmaMass = -1;
  }

AliFemtoWRzTrackCut::AliFemtoWRzTrackCut(const AliFemtoWRzTrackCut &aCut) : AliFemtoESDTrackCut(aCut)
{
    //copy constructor 
}

AliFemtoWRzTrackCut::~AliFemtoWRzTrackCut()
{
  // Destructor
}


AliFemtoWRzTrackCut& AliFemtoWRzTrackCut::operator=(const AliFemtoWRzTrackCut& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoESDTrackCut::operator=(aCut);
 
  return *this;
}

bool AliFemtoWRzTrackCut::Pass( const AliFemtoTrack* track)
{
  // almost the same function as in AliFemtoESDTrackCut. Differences in *N Sigma Method* and *Countur Method*

  // test the particle and return
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria

  if (fStatus && (track->Flags() & fStatus) != fStatus) {
    return false;
  }
  if (fRemoveKinks && (track->KinkIndex(0) || track->KinkIndex(1) || track->KinkIndex(2))) {
    return false;
  }
  if (fRemoveITSFake && track->ITSncls() < 0) {
    return false;
  }
  if (fminTPCclsF > track->TPCnclsF()) {
    return false;
  }
  if (fminTPCncls > track->TPCncls()) {
    return false;
  }
  if (fminITScls > track->ITSncls()) {
    return false;
  }

  if (fMaxImpactXY < TMath::Abs(track->ImpactD())) {
    return false;
  }
  if (fMinImpactXY > TMath::Abs(track->ImpactD())) {
    return false;
  }
  if (fMaxImpactZ < TMath::Abs(track->ImpactZ())) {
    return false;
  }
  if (fMaxSigmaToVertex < track->SigmaToVertex()) {
    return false;
  }

  if (track->ITSncls() > 0 && (track->ITSchi2() / track->ITSncls()) > fMaxITSchiNdof) {
    return false;
  }

  if (track->TPCchi2perNDF() > fMaxTPCchiNdof) {
    return false;
  }

  // ITS cluster requirenments
  for (Int_t i = 0; i < 3; i++) {
    if (!CheckITSClusterRequirement(fCutClusterRequirementITS[i], track->HasPointOnITSLayer(i * 2), track->HasPointOnITSLayer(i*2+1))) {
      return false;
    }
  }

  if (fLabel) {
    if (track->Label() < 0) {
      fNTracksFailed++;
      return false;
    }
  }
  if (fCharge != 0 && (track->Charge() != fCharge)) {
    fNTracksFailed++;
    return false;
  }


  Bool_t tTPCPidIn = (track->Flags() & AliFemtoTrack::kTPCpid) > 0;
  Bool_t tITSPidIn = (track->Flags() & AliFemtoTrack::kITSpid) > 0;
  Bool_t tTOFPidIn = (track->Flags() & AliFemtoTrack::kTOFpid) > 0;

  const double momentum = track->P().Mag();

  if (fMinPforTOFpid > 0
      && fMinPforTOFpid < momentum && momentum < fMaxPforTOFpid
      && !tTOFPidIn) {
    fNTracksFailed++;
    return false;
  }

  if (fMinPforTPCpid > 0
      && fMinPforTPCpid < momentum && momentum < fMaxPforTPCpid
      && !tTPCPidIn) {
    fNTracksFailed++;
    return false;
  }

  if (fMinPforITSpid > 0
      && fMinPforITSpid < momentum && momentum < fMaxPforITSpid
      && !tITSPidIn) {
    fNTracksFailed++;
    return false;
  }


  float tEnergy = ::sqrt(track->P().Mag2() + fMass * fMass);
  float tRapidity = 0;
  if (tEnergy-track->P().z() != 0 && (tEnergy + track->P().z()) / (tEnergy-track->P().z()) > 0)
    tRapidity = 0.5 * ::log((tEnergy + track->P().z())/(tEnergy-track->P().z()));
  float tPt = track->P().Perp();
  float tEta = track->P().PseudoRapidity();

  if (fMaxImpactXYPtOff < 999.0) {
    if ((fMaxImpactXYPtOff + fMaxImpactXYPtNrm*TMath::Power(tPt, fMaxImpactXYPtPow)) < TMath::Abs(track->ImpactD())) {
      fNTracksFailed++;
      return false;
    }
  }

  if ((tRapidity < fRapidity[0]) || (tRapidity > fRapidity[1])) {
    fNTracksFailed++;
    return false;
  }
  if ((tEta < fEta[0]) || (tEta > fEta[1])) {
    fNTracksFailed++;
    return false;
  }
  if ((tPt < fPt[0]) || (tPt > fPt[1])) {
    fNTracksFailed++;
    return false;
  }


  if ((track->PidProbElectron() < fPidProbElectron[0]) || (track->PidProbElectron() > fPidProbElectron[1])) {
    fNTracksFailed++;
    return false;
  }
  if ((track->PidProbPion() < fPidProbPion[0]) || (track->PidProbPion() > fPidProbPion[1])) {
    fNTracksFailed++;
    return false;
  }
  if ((track->PidProbKaon() < fPidProbKaon[0]) || (track->PidProbKaon() > fPidProbKaon[1])) {
    fNTracksFailed++;
    return false;
  }
  if ((track->PidProbProton() < fPidProbProton[0]) || (track->PidProbProton() > fPidProbProton[1])) {
    fNTracksFailed++;
    return false;
  }
  if ((track->PidProbMuon() < fPidProbMuon[0]) || (track->PidProbMuon() > fPidProbMuon[1])) {
    fNTracksFailed++;
    return false;
  }


  //****N Sigma Method -- electron rejection****
  if (fElectronRejection)
    if (!IsElectron(track->NSigmaTPCE(),track->NSigmaTPCPi(),track->NSigmaTPCK(), track->NSigmaTPCP()))
      return false;


  if (fMostProbable) {
    int imost=0;

    float ipidmax = 0.0;
    //****N Sigma Method****
    if (fPIDMethod==0) {
      // Looking for pions
      if (fMostProbable == 2) {
        if (IsPionNSigma(track->P().Mag(), track->NSigmaTPCPi(), track->NSigmaTOFPi())) {
          imost = 2;
        }
      }
      else if (fMostProbable == 3) {
        if (IsKaonNSigma(track->P().Mag(), track->NSigmaTPCK(), track->NSigmaTOFK())){
          imost = 3;
        }
      }
      else if (fMostProbable == 4) { // proton nsigma-PID required contour adjusting (in LHC10h)
        if (IsProtonNSigma(track->P().Mag(), track->NSigmaTPCP(), track->NSigmaTOFP()) ) {
          imost = 4;
        }
      }
      if (fMostProbable == 10) {//cut on Nsigma in pT not p
        if (IsPionNSigma(track->Pt(), track->NSigmaTPCPi(), track->NSigmaTOFPi()))
          imost = 10;
      }
      else if (fMostProbable == 11) {//cut on Nsigma in pT not p
        if (IsKaonNSigma(track->Pt(), track->NSigmaTPCK(), track->NSigmaTOFK())){
          imost = 11;
        }
      }
      else if (fMostProbable == 12) { //cut on Nsigma in pT not p
        if ( IsProtonNSigma(track->Pt(), track->NSigmaTPCP(), track->NSigmaTOFP()) )
          imost = 12;
      }
      else if (fMostProbable == 13) { //cut on Nsigma deuteron
        if (IsDeuteronNSigma(track->P().Mag(),track->MassTOF(), fNsigmaMass, track->NSigmaTPCD(), track->NSigmaTOFD()) )
          imost = 13;
        if ( fdEdxcut && (track->P().Mag() < 3) && !IsDeuteronTPCdEdx(track->P().Mag(), track->TPCsignal(), 4) )
          imost = 0;
      }
      else if (fMostProbable == 14) { //cut on Nsigma, EXCLUSIVE PID -- deuteron 
        if ( IsDeuteronNSigma(track->P().Mag(),track->MassTOF(), fNsigmaMass, track->NSigmaTPCD(), track->NSigmaTOFD()) && !IsPionNSigmaRejection(track->P().Mag(),track->NSigmaTPCPi(), track->NSigmaTOFPi()) && !IsKaonNSigmaRejection(track->P().Mag(),track->NSigmaTPCK(), track->NSigmaTOFK()) && !IsProtonNSigmaRejection(track->P().Mag(),track->NSigmaTPCP(), track->NSigmaTOFP()) && !IsElectronNSigmaRejection(track->P().Mag(),track->NSigmaTPCE()) )
          imost = 14;
        if ( fdEdxcut && (track->P().Mag() < 3) && !IsDeuteronTPCdEdx(track->P().Mag(), track->TPCsignal(), 4) )
          imost = 0;
      }
    }

    //****Contour Method****
    if (fPIDMethod==1) {
       if (fMostProbable == 13 || fMostProbable == 14) {
        // Using the TPC dEdx to reject non-deuterons
         if (track->P().Mag() < 2) {
           if ( !(IsDeuteronTPCdEdx(track->P().Mag(), track->TPCsignal(),4)) )
             imost = 0;
           else 
             imost = 13;  
         }
       }

    }

    if (imost != fMostProbable) {
      return false;
    }
  }
  fNTracksPassed++ ;
  return true;


}

bool AliFemtoWRzTrackCut::IsDeuteronTPCdEdx(float mom, float dEdx, float maxmom)
{

  double a1 = -250.0,  b1 = 400.0;
  double a2 = -135.0,  b2 = 270.0;
  double a3 = -80,   b3 = 190.0;
  double a4 = 0.0,   b4 = 20.0;

  double a5 = 125.0,   b5 = -100.0; 

  if (mom < 1.1) {
    if (dEdx < a1*mom+b1) return false;
  }
  else if (mom < 1.4) {
    if (dEdx < a2*mom+b2) return false;
  }
  else if (mom < 2) {
    if (dEdx < a3*mom+b3) return false;
  }
  else if (mom >= 2) {
    if (dEdx < a4*mom+b4) return false;
  }

  if (mom > maxmom) return false;
  
  if (fNsigmaTPConly && fNsigmaMass<0) {
    // for selection with only the TPC detector
    // cutting out the final part of the signal (~1.5 GeV) 
    // that is dominated by missidentified particles 
    // this setting will be probably removed  
    if (dEdx < a5*mom+b5) return false;
  }

  return true;

}

bool AliFemtoWRzTrackCut::IsDeuteronNSigma(float mom, float massTOFPDG,float sigmaMass, float nsigmaTPCD, float nsigmaTOFD)
{
  double massPDGD=1.8756;
  if (fNsigmaTPCTOF) {
    //Identyfication with only TPC for mom<1.4 and TPC&TOF for mom>1.4
    if (mom > 1.4){
      if ((TMath::Abs(nsigmaTPCD) < fNsigma) && (TMath::Abs(nsigmaTOFD) < fNsigma))
        return true;
    }
    else{
      if (TMath::Abs(nsigmaTPCD) < fNsigma)
        return true;
    } 
  }
  else if (fNsigmaTPConly){

      if (TMath::Abs(nsigmaTPCD) < fNsigma){
           if(sigmaMass==-1){
              return true;
           } 
           else{
              //strict nsigma selection. removing visible contamination with TPCnsigma distribution
              //this setting is for tests (for now)
              double line1 = -18.41*mom*mom+56.37*mom-43.59;
              double line2 = 0.235*mom-0.79;
              if(nsigmaTPCD>0)
                 return true; 
              else if (mom<=1.5 && nsigmaTPCD>line1)
                   return true; 
              else if(mom>1.5 && nsigmaTPCD>line2)
                   return true; 
           }
      }

  }
  else{// p dependent mass cut 
       // The default setting of sigmaMass (= -1) should provide similar 
       // resluts to the case with fNsigmaTPCTOF=true and fNsigma=2

      double l1, l2;
      if(sigmaMass > 2){//for sideband analysis -- left band (4 sigmas)
         l1 = 2.897 - 0.558*mom + 0.177*mom*mom - 0.026*mom*mom*mom;  
         l2 = 3.77 - 0.487*mom + 0.126*mom*mom - 0.014*mom*mom*mom;
      }
      else if(sigmaMass > 1){//for sideband analysis -- right band (4 sigmas)
         l1 = 4.5 - 0.42*mom + 0.1*mom*mom;
         l2 = 6.602 -0.981*mom + 0.303*mom*mom;
      }
      else{//signal (dist under the mass peak -- 2 sigmas)
         l1 = 4.002 - 0.627*mom + 0.184*mom*mom - 0.02*mom*mom*mom;
         l2 = 4.35 - 0.399*mom + 0.09*mom*mom;
      }

      if ((TMath::Abs(nsigmaTPCD) < fNsigma) && (massTOFPDG > l1) && (massTOFPDG < l2))
	 return true;

  }
  return false;
}


//rejection methods
bool AliFemtoWRzTrackCut::IsPionNSigmaRejection(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if(fNsigmaTPConly){
    if(TMath::Abs(nsigmaTPCPi) < fNsigmaRejection)
      return true;
  }
  else{
    if(TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < fNsigmaRejection)
      return true;
  }

  return false;
}


bool AliFemtoWRzTrackCut::IsKaonNSigmaRejection(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if(fNsigmaTPConly){
    if(TMath::Abs(nsigmaTPCK) < fNsigmaRejection)
      return true;
  }
  else{
    if(TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < fNsigmaRejection)
      return true;
  }

  return false;
}

bool AliFemtoWRzTrackCut::IsProtonNSigmaRejection(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
  if(fNsigmaTPConly){
    if(TMath::Abs(nsigmaTPCP) < fNsigmaRejection)
      return true;
  }
  else{
    if(TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigmaRejection)
      return true;
  }
 
  return false;
}

bool AliFemtoWRzTrackCut::IsElectronNSigmaRejection(float mom, float nsigmaTPCE)
{
  if(TMath::Abs(nsigmaTPCE) < fNsigmaRejection)
    return true;
 

  return false;
}




