//==============================================================\\
// dowang track cut for p-d/t/He3 analysis                      \\
// deuteron part refer to AliFemtoWRzTrackCut.h                 \\
//==============================================================\\

#include "AliFemtoTrackCutPdtHe3.h"
AliFemtoTrackCutPdtHe3::AliFemtoTrackCutPdtHe3():
AliFemtoESDTrackCut()
{
    fNsigmaP    = 3.;
    fNsigmaD    = 3.;
    fNsigmaT    = 3.;
    fNsigmaHe3  = 3.;
    fNsigmaRejection = 3.;
    
    SwitchMom_p = 0.5;
    SwitchMom_d = 0.8;
    SwitchMom_t = 1.2;
    SwitchMom_He3 = 999;
}

AliFemtoTrackCutPdtHe3::AliFemtoTrackCutPdtHe3(const AliFemtoTrackCutPdtHe3 &aCut) : 
AliFemtoESDTrackCut(aCut)
{
    //copy constructor 
    fNsigmaP = aCut.fNsigmaP;
    fNsigmaD = aCut.fNsigmaD;
    fNsigmaT = aCut.fNsigmaT;
    fNsigmaHe3 = aCut.fNsigmaHe3;
    fNsigmaRejection = aCut.fNsigmaRejection;

    SwitchMom_p = aCut.SwitchMom_p;
    SwitchMom_d = aCut.SwitchMom_d;
    SwitchMom_t = aCut.SwitchMom_t;
    SwitchMom_He3 = aCut.SwitchMom_He3;
}

AliFemtoTrackCutPdtHe3::~AliFemtoTrackCutPdtHe3()
{
  // Destructor
}


AliFemtoTrackCutPdtHe3& AliFemtoTrackCutPdtHe3::operator=(const AliFemtoTrackCutPdtHe3& aCut)
{
    // assignment operator
    if (this == &aCut)
        return *this;

    AliFemtoESDTrackCut::operator=(aCut);

    fNsigmaP = aCut.fNsigmaP;
    fNsigmaD = aCut.fNsigmaD;
    fNsigmaT = aCut.fNsigmaD;
    fNsigmaHe3 = aCut.fNsigmaHe3;

    fNsigmaRejection = aCut.fNsigmaRejection;
    SwitchMom_p = aCut.SwitchMom_p;
    SwitchMom_d = aCut.SwitchMom_d;
    SwitchMom_t = aCut.SwitchMom_t;
    SwitchMom_He3 = aCut.SwitchMom_He3;

    return *this;
}

bool AliFemtoTrackCutPdtHe3::Pass(const AliFemtoTrack* track)
{
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
        //****N Sigma Method****
        if (fPIDMethod==0) {
            if (fMostProbable == 4) { // proton nsigma-PID required contour adjusting (in LHC10h)
                if (IsProtonNSigma(track->P().Mag(), track->NSigmaTPCP(), track->NSigmaTOFP()) ) {
                    imost = 4;
                }
            }
            else if (fMostProbable == 13){   //cut on Nsigma deuteron
                if (IsDeuteronNSigma(track->P().Mag(), track->MassTOF(), fNsigmaMass, track->NSigmaTPCD(), track->NSigmaTOFD())
                    && !IsElectronNSigmaRejection(track->P().Mag(),track->NSigmaTPCE())
                    && !IsPionNSigmaRejection(track->P().Mag(),track->NSigmaTPCPi(), track->NSigmaTOFPi())
                    && !IsKaonNSigmaRejection(track->P().Mag(),track->NSigmaTPCK(), track->NSigmaTOFK())
                    && !IsProtonNSigmaRejection(track->P().Mag(),track->NSigmaTPCP(), track->NSigmaTOFP())){
                        imost = 13;
                }
            }
            else if (fMostProbable == 14){   //cut on Nsigma triton
                if (IsTritonNSigma(track->P().Mag(), track->MassTOF(), fNsigmaMass, track->NSigmaTPCT(), track->NSigmaTOFT())
                    && !IsElectronNSigmaRejection(track->P().Mag(),track->NSigmaTPCE())
                    && !IsPionNSigmaRejection(track->P().Mag(),track->NSigmaTPCPi(), track->NSigmaTOFPi())
                    && !IsKaonNSigmaRejection(track->P().Mag(),track->NSigmaTPCK(), track->NSigmaTOFK())
                    && !IsProtonNSigmaRejection(track->P().Mag(),track->NSigmaTPCP(), track->NSigmaTOFP())){
                    imost = 14;
                }
                   
            }
            else if (fMostProbable == 15){
                if (IsHe3NSigma(track->P().Mag(), track->MassTOF(), fNsigmaMass, track->NSigmaTPCH(), track->NSigmaTOFH())
                    && !IsElectronNSigmaRejection(track->P().Mag(),track->NSigmaTPCE())
                    && !IsPionNSigmaRejection(track->P().Mag(),track->NSigmaTPCPi(), track->NSigmaTOFPi())
                    && !IsKaonNSigmaRejection(track->P().Mag(),track->NSigmaTPCK(), track->NSigmaTOFK())
                    && !IsProtonNSigmaRejection(track->P().Mag(),track->NSigmaTPCP(), track->NSigmaTOFP())){
                        imost = 15;
                }

            }
            if (imost != fMostProbable) {
                return false;
            }
            
        }
    fNTracksPassed++ ;
   
    return true;
    
    }
    
}
bool AliFemtoTrackCutPdtHe3::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
    if (mom > SwitchMom_p) {
	        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigma)
	            return true;	
	}
    else {
        if (TMath::Abs(nsigmaTPCP) < fNsigma)
            return true;
    }
    return false;

}
bool AliFemtoTrackCutPdtHe3::IsDeuteronNSigma(float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCD, float nsigmaTOFD)
{
    //double massPDGD=1.8756;
    if (fNsigmaTPCTOF) {
        //Identyfication with only TPC for mom<1.4 and TPC&TOF for mom>1.4
        if (mom > SwitchMom_d){
            if (TMath::Hypot( nsigmaTPCD, nsigmaTOFD ) < fNsigma)
            //if ((TMath::Abs(nsigmaTPCD) < fNsigma) && (TMath::Abs(nsigmaTOFD) < fNsigma))
                return true;
        }
        else{
            if (TMath::Abs(nsigmaTPCD) < fNsigma)
                return true;
        }
        
    }
    else{
        if (TMath::Abs(nsigmaTPCD) < fNsigma)
            return true;
    } 
    return false;
}
bool AliFemtoTrackCutPdtHe3::IsTritonNSigma(float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCT, float nsigmaTOFT)
{
    //double massPDGD=2.8089;
    if (fNsigmaTPCTOF) {
        //Identyfication with only TPC for mom<1.4 and TPC&TOF for mom>1.4
        if (mom > SwitchMom_t){
            if ((TMath::Abs(nsigmaTPCT) < fNsigma) && (TMath::Abs(nsigmaTOFT) < fNsigma))
                return true;
        }
        else{
            if (TMath::Abs(nsigmaTPCT) < fNsigma)
                return true;
        }
        
    }
    else{
        if (TMath::Abs(nsigmaTPCT) < fNsigma)
            return true;
    } 
    return false;
}
bool AliFemtoTrackCutPdtHe3::IsHe3NSigma(float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCHe3, float nsigmaTOFHe3)
{
    //double massPDGD=2.8089;
    if (fNsigmaTPCTOF) {
        if (mom > SwitchMom_He3){
            if ((TMath::Abs(nsigmaTPCHe3) < fNsigma) && (TMath::Abs(nsigmaTOFHe3) < fNsigma))
                return true;
        }
        else{
            if (TMath::Abs(nsigmaTPCHe3) < fNsigma)
                return true;
        }
    }
    else{
        if (TMath::Abs(nsigmaTPCHe3) < fNsigma)
            return true;
    }
    return false;
  
}

//rejection methods
bool AliFemtoTrackCutPdtHe3::IsElectronNSigmaRejection(float mom, float nsigmaTPCE)
{
  if(TMath::Abs(nsigmaTPCE) < fNsigmaRejection)
    return true;
 

  return false;
}

bool AliFemtoTrackCutPdtHe3::IsPionNSigmaRejection(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
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


bool AliFemtoTrackCutPdtHe3::IsKaonNSigmaRejection(float mom, float nsigmaTPCK, float nsigmaTOFK)
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

bool AliFemtoTrackCutPdtHe3::IsProtonNSigmaRejection(float mom, float nsigmaTPCP, float nsigmaTOFP)
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



