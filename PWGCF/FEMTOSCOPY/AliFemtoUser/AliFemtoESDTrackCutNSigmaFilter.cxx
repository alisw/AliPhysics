////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFemtoESDTrackCutNSigmaFilter:                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include "AliFemtoESDTrackCutNSigmaFilter.h"

#ifdef __ROOT__ 
ClassImp(AliFemtoESDTrackCutNSigmaFilter)
#endif




//________________________________________________________________________________________________________________
AliFemtoESDTrackCutNSigmaFilter::AliFemtoESDTrackCutNSigmaFilter() :
  AliFemtoESDTrackCut()
  , fUseCustomPionNSigmaFilter(false)
  , fUseCustomKaonNSigmaFilter(false)
  , fUseCustomProtonNSigmaFilter(false)
  , fUseCustomElectronNSigmaFilter(false)

  , fPionNSigmaFilter(NULL)
  , fKaonNSigmaFilter(NULL)
  , fProtonNSigmaFilter(NULL)
  , fElectronNSigmaFilter(NULL)
  , fPionRejection(false)
{

}

//________________________________________________________________________________________________________________
AliFemtoESDTrackCutNSigmaFilter::AliFemtoESDTrackCutNSigmaFilter(const AliFemtoESDTrackCutNSigmaFilter& aCut):
  AliFemtoESDTrackCut(aCut),
  fUseCustomPionNSigmaFilter(aCut.fUseCustomPionNSigmaFilter),
  fUseCustomKaonNSigmaFilter(aCut.fUseCustomKaonNSigmaFilter),
  fUseCustomProtonNSigmaFilter(aCut.fUseCustomProtonNSigmaFilter),
  fUseCustomElectronNSigmaFilter(aCut.fUseCustomElectronNSigmaFilter),
  fPionRejection(aCut.fPionRejection)
{
  if(aCut.fPionNSigmaFilter) fPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fPionNSigmaFilter);
  if(aCut.fKaonNSigmaFilter) fKaonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fKaonNSigmaFilter);
  if(aCut.fProtonNSigmaFilter) fProtonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fProtonNSigmaFilter);
  if(aCut.fElectronNSigmaFilter) fElectronNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fElectronNSigmaFilter);
}


//________________________________________________________________________________________________________________
AliFemtoESDTrackCutNSigmaFilter& AliFemtoESDTrackCutNSigmaFilter::operator=(const AliFemtoESDTrackCutNSigmaFilter& aCut)
{
  if(this == &aCut) {return *this;}

  AliFemtoESDTrackCut::operator=(aCut);

  fUseCustomPionNSigmaFilter = aCut.fUseCustomPionNSigmaFilter;
  fUseCustomKaonNSigmaFilter = aCut.fUseCustomKaonNSigmaFilter;
  fUseCustomProtonNSigmaFilter = aCut.fUseCustomProtonNSigmaFilter;
  fUseCustomElectronNSigmaFilter = aCut.fUseCustomElectronNSigmaFilter;
  fPionRejection = aCut.fPionRejection;

  if(aCut.fPionNSigmaFilter) fPionNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fPionNSigmaFilter);
  if(aCut.fKaonNSigmaFilter) fKaonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fKaonNSigmaFilter);
  if(aCut.fProtonNSigmaFilter) fProtonNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fProtonNSigmaFilter);
  if(aCut.fElectronNSigmaFilter) fElectronNSigmaFilter = new AliFemtoNSigmaFilter(*aCut.fElectronNSigmaFilter);

  return *this;
}

//________________________________________________________________________________________________________________
AliFemtoESDTrackCutNSigmaFilter* AliFemtoESDTrackCutNSigmaFilter::Clone()
{
  return(new AliFemtoESDTrackCutNSigmaFilter(*this));
}

//________________________________________________________________________________________________________________
AliFemtoESDTrackCutNSigmaFilter::~AliFemtoESDTrackCutNSigmaFilter()
{
  cout << "AliFemtoESDTrackCutNSigmaFilter object is being deleted!!!" << endl;
}

//________________________________________________________________________________________________________________
bool AliFemtoESDTrackCutNSigmaFilter::Pass(const AliFemtoTrack* track)
{
  //cout<<"AliFemtoESDTrackCut::Pass"<<endl;

  // test the particle and return
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria
  float tMost[5];

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

  if (track->TPCncls() > 0 && (track->TPCchi2() / track->TPCncls()) > fMaxTPCchiNdof) {
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




  //   cout << "Track has pids: "
  //        << track->PidProbElectron() << " "
  //        << track->PidProbMuon() << " "
  //        << track->PidProbPion() << " "
  //        << track->PidProbKaon() << " "
  //        << track->PidProbProton() << " "
  //        << track->PidProbElectron()+track->PidProbMuon()+track->PidProbPion()+track->PidProbKaon()+track->PidProbProton() << endl;


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
  if(fElectronRejection)
  {
    if(fUseCustomElectronNSigmaFilter)
    {
//      if(IsElectronNSigma(track->P().Mag(), track->NSigmaTPCE(), track->NSigmaTOFE())) return false;
      if(IsProbableElectron(track)) return false;
    }
    else
    {
      if(!AliFemtoESDTrackCut::IsElectron(track->NSigmaTPCE(),track->NSigmaTPCPi(),track->NSigmaTPCK(), track->NSigmaTPCP())) return false;
    }
  }

  if(fPionRejection)
  {
    if(IsPionNSigma(track->P().Mag(), track->NSigmaTPCPi(), track->NSigmaTOFPi())) return false;
  }

  if (fMostProbable) {

    int imost=0;
    tMost[0] = track->PidProbElectron()*PidFractionElectron(track->P().Mag());
    tMost[1] = 0.0;
    tMost[2] = track->PidProbPion()*PidFractionPion(track->P().Mag());
    tMost[3] = track->PidProbKaon()*PidFractionKaon(track->P().Mag());
    tMost[4] = track->PidProbProton()*PidFractionProton(track->P().Mag());
    float ipidmax = 0.0;

    //****N Sigma Method****
	if(fPIDMethod==0){
	  // Looking for pions
	  if (fMostProbable == 2) {
	    if (IsPionNSigma(track->P().Mag(), track->NSigmaTPCPi(), track->NSigmaTOFPi()))
	      imost = 2;

	  }
	  else if (fMostProbable == 3) {


	    if (IsKaonNSigma(track->P().Mag(), track->NSigmaTPCK(), track->NSigmaTOFK())){

	      imost = 3;
	    }
	  }
	  else if (fMostProbable == 4) { // proton nsigma-PID required contour adjusting (in LHC10h)
	    if ( IsProtonNSigma(track->P().Mag(), track->NSigmaTPCP(), track->NSigmaTOFP()) // && (TMath::Abs(track->NSigmaTPCP()) < TMath::Abs(track->NSigmaTPCPi())) && (TMath::Abs(track->NSigmaTPCP()) < TMath::Abs(track->NSigmaTPCK())) && (TMath::Abs(track->NSigmaTOFP()) < TMath::Abs(track->NSigmaTOFPi())) && (TMath::Abs(track->NSigmaTOFP()) < TMath::Abs(track->NSigmaTOFK()))
		 // && IsProtonTPCdEdx(track->P().Mag(), track->TPCsignal())
		 )
	      imost = 4;
	  }
	  else if (fMostProbable == 5) { // no-protons
	    if ( !IsProtonNSigma(track->P().Mag(), track->NSigmaTPCP(), track->NSigmaTOFP()) )
	      imost = 5;
	  }
	  else if (fMostProbable == 6) { //pions OR kaons OR protons
	    if (IsPionNSigma(track->P().Mag(), track->NSigmaTPCPi(), track->NSigmaTOFPi()))
	      imost = 6;
	    else if (IsKaonNSigma(track->P().Mag(), track->NSigmaTPCK(), track->NSigmaTOFK()))
	      imost = 6;
	    else if (IsProtonNSigma(track->P().Mag(), track->NSigmaTPCP(), track->NSigmaTOFP()) )
	      imost = 6;
	  }
	  else if (fMostProbable == 7) { // pions OR kaons OR protons OR electrons
	    if (IsPionNSigma(track->P().Mag(), track->NSigmaTPCPi(), track->NSigmaTOFPi()))
	      imost = 7;
	    else if (IsKaonNSigma(track->P().Mag(), track->NSigmaTPCK(), track->NSigmaTOFK()))
	      imost = 7;
	    else if (IsProtonNSigma(track->P().Mag(), track->NSigmaTPCP(), track->NSigmaTOFP()) )
	      imost = 7;
	    else if (TMath::Abs(track->NSigmaTPCE())<3)
	      imost = 7;

	  }
	  else if (fMostProbable == 8) { // TOF matching
	    if(track->NSigmaTOFPi() != -1000 || track->Pt()<0.5){
	      imost = 8;
	    }
	  }
	  else if (fMostProbable == 9) { // Other: no kaons, no pions, no protons
	    if (IsPionNSigma(track->P().Mag(), track->NSigmaTPCPi(), track->NSigmaTOFPi()))
	      imost = -1;
	    else if (IsKaonNSigma(track->P().Mag(), track->NSigmaTPCK(), track->NSigmaTOFK()))
	      imost = -1;
	    else if (IsProtonNSigma(track->P().Mag(), track->NSigmaTPCP(), track->NSigmaTOFP()) )
	      imost = -1;
	    else if(track->NSigmaTOFPi() != -1000 || track->Pt()<0.5){
	      imost = 9;
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


	}



    //****Contour Method****
	if(fPIDMethod==1){
	  for (int ip=0; ip<5; ip++)
	    if (tMost[ip] > ipidmax) { ipidmax = tMost[ip]; imost = ip; };

	  // Looking for pions
	  if (fMostProbable == 2) {
	    if (imost == 2) {
	      // Using the TPC to reject non-pions
	      if (!(IsPionTPCdEdx(track->P().Mag(), track->TPCsignal())))
		imost = 0;
	      if (0) {
		// Using the TOF to reject non-pions
		if (track->P().Mag() < 0.6) {
		  if (tTOFPidIn)
		    if (!IsPionTOFTime(track->P().Mag(), track->TOFpionTime()))
		      imost = 0;
		}
		else {
		  if (tTOFPidIn) {
		    if (!IsPionTOFTime(track->P().Mag(), track->TOFpionTime()))
		      imost = 0;
		  }
		  else {
		    imost = 0;
		  }
		}
	      }
	    }
	  }

	  // Looking for kaons
	  else if (fMostProbable == 3) {
	    //       if (imost == 3) {
	    // Using the TPC to reject non-kaons
	    if (track->P().Mag() < 0.6) {
	      if (!(IsKaonTPCdEdx(track->P().Mag(), track->TPCsignal())))
		imost = 0;
	      else imost = 3;
	      if (1) {
		// Using the TOF to reject non-kaons
		if (tTOFPidIn)
		  if (!IsKaonTOFTime(track->P().Mag(), track->TOFkaonTime()))
		    imost = 0;
	      }
	    }
	    else {
	      if (1) {
		if (tTOFPidIn) {
		  if (!IsKaonTOFTime(track->P().Mag(), track->TOFkaonTime()))
		    imost = 0;
		  else
		    imost = 3;
		}
		else {
		  if (!(IsKaonTPCdEdx(track->P().Mag(), track->TPCsignal())))
		    imost = 0;
		  else
		    imost = 3;
		}
	      }
	    }
	    //       }
	  }

	  // Looking for protons
	  else if (fMostProbable == 4) {
	    //       if (imost == 3) {
	    // Using the TPC to reject non-kaons
	    if (track->P().Mag() < 0.8) {
	      if (!(IsProtonTPCdEdx(track->P().Mag(), track->TPCsignal())))
		imost = 0;
	      else imost = 4;
	      if (0) {
		// Using the TOF to reject non-kaons
		if (tTOFPidIn)
		  if (!IsKaonTOFTime(track->P().Mag(), track->TOFkaonTime()))
		    imost = 0;
	      }
	    }
	    else {
	      if (0) {
		if (tTOFPidIn) {
		  if (!IsKaonTOFTime(track->P().Mag(), track->TOFkaonTime()))
		    imost = 0;
		  else
		    imost = 3;
		}
		else {
		  if (!(IsKaonTPCdEdx(track->P().Mag(), track->TPCsignal())))
		    imost = 0;
		  else
		    imost = 3;
		}
	      }
	    }
	    //       }
	  }
	}
    if (imost != fMostProbable) return false;
  }


  //fan
  //cout<<"****** Go Through the cut ******"<<endl;
  // cout<<fLabel<<" Label="<<track->Label()<<endl;
  // cout<<fCharge<<" Charge="<<track->Charge()<<endl;
  // cout<<fPt[0]<<" < Pt ="<<Pt<<" <"<<fPt[1]<<endl;
  //cout<<fRapidity[0]<<" < Rapidity ="<<tRapidity<<" <"<<fRapidity[1]<<endl;
  //cout<<fPidProbElectron[0]<<" <  e="<<track->PidProbElectron()<<" <"<<fPidProbElectron[1]<<endl;
  //cout<<fPidProbPion[0]<<" <  pi="<<track->PidProbPion()<<" <"<<fPidProbPion[1]<<endl;
  //cout<<fPidProbKaon[0]<<" <  k="<<track->PidProbKaon()<<" <"<<fPidProbKaon[1]<<endl;
  //cout<<fPidProbProton[0]<<" <  p="<<track->PidProbProton()<<" <"<<fPidProbProton[1]<<endl;
  //cout<<fPidProbMuon[0]<<" <  mi="<<track->PidProbMuon()<<" <"<<fPidProbMuon[1]<<endl;
  fNTracksPassed++ ;
  return true;


}

void AliFemtoESDTrackCutNSigmaFilter::CreateCustomNSigmaFilter(ParticleType aType)
{
  switch (aType) {
  case kPion:
    fUseCustomPionNSigmaFilter = true;
    fPionNSigmaFilter = new AliFemtoNSigmaFilter();
    break;

  case kKaon:
    fUseCustomKaonNSigmaFilter = true;
    fKaonNSigmaFilter = new AliFemtoNSigmaFilter();
    break;

  case kProton:
    fUseCustomProtonNSigmaFilter = true;
    fProtonNSigmaFilter = new AliFemtoNSigmaFilter();
    break;

  case kElectron:
    fUseCustomElectronNSigmaFilter = true;
    fElectronNSigmaFilter = new AliFemtoNSigmaFilter();
    break;

  default:
    cerr << "E-AliFemtoESDTrackCutNSigmaFilter::CreateCustomNSigmaFilter: Invalid ParticleType"
            "selection '" << aType << "'.  No custom filter will be initialized!!!!!" << endl;
    break;
  }
}


void AliFemtoESDTrackCutNSigmaFilter::AddTPCAndTOFNSigmaCut(ParticleType aType, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  switch (aType) {
  case kPion:
    fPionNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
    break;

  case kKaon:
    fKaonNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
    break;

  case kProton:
    fProtonNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
    break;

  case kElectron:
    fElectronNSigmaFilter->AddTPCAndTOFCut(aMomMin,aMomMax,aNSigmaValueTPC,aNSigmaValueTOF);
    break;

  default:
    cerr << "E-AliFemtoESDTrackCutNSigmaFilter::AddTPCAndTOFNSigmaCut: Invalid ParticleType"
            "selection '" << aType << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoESDTrackCutNSigmaFilter::AddTPCNSigmaCut(ParticleType aType, double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  switch (aType) {
  case kPion:
    fPionNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);
    break;

  case kKaon:
    fKaonNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);
    break;

  case kProton:
    fProtonNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);
    break;

  case kElectron:
    fElectronNSigmaFilter->AddTPCCut(aMomMin,aMomMax,aNSigmaValueTPC);
    break;

  default:
    cerr << "E-AliFemtoESDTrackCutNSigmaFilter::AddTPCNSigmaCut: Invalid ParticleType"
            "selection '" << aType << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoESDTrackCutNSigmaFilter::AddTOFNSigmaCut(ParticleType aType, double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  switch (aType) {
  case kPion:
    fPionNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);
    break;

  case kKaon:
    fKaonNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);
    break;

  case kProton:
    fProtonNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);
    break;

  case kElectron:
    fElectronNSigmaFilter->AddTOFCut(aMomMin,aMomMax,aNSigmaValueTOF);
    break;

  default:
    cerr << "E-AliFemtoESDTrackCutNSigmaFilter::AddTOFNSigmaCut: Invalid ParticleType"
            "selection '" << aType << "'.  No cut will be initialized!!!!!" << endl;
    break;
  }
}

void AliFemtoESDTrackCutNSigmaFilter::SetOverrideImproperPionNSigmaFilter(ParticleType aType, bool aOverride)
{
  switch (aType) {
  case kPion:
    fPionNSigmaFilter->SetOverrideImproperConfig(aOverride);
    break;

  case kKaon:
    fKaonNSigmaFilter->SetOverrideImproperConfig(aOverride);
    break;

  case kProton:
    fProtonNSigmaFilter->SetOverrideImproperConfig(aOverride);
    break;

  case kElectron:
    fElectronNSigmaFilter->SetOverrideImproperConfig(aOverride);
    break;

  default:
    cerr << "E-AliFemtoESDTrackCutNSigmaFilter::SetOverrideImproperPionNSigmaFilter: Invalid ParticleType"
            "selection '" << aType << "'.  No override call will be made!!!!!" << endl;
    break;
  }
}

//________________________________________________________________________________________________________________
bool AliFemtoESDTrackCutNSigmaFilter::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if (fUseCustomKaonNSigmaFilter) {
    return fKaonNSigmaFilter->Pass(mom, nsigmaTPCK, nsigmaTOFK);
  } else {
    return AliFemtoESDTrackCut::IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
  }
}

//________________________________________________________________________________________________________________
bool AliFemtoESDTrackCutNSigmaFilter::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if (fUseCustomPionNSigmaFilter) {
    return fPionNSigmaFilter->Pass(mom, nsigmaTPCPi, nsigmaTOFPi);
  } else {
    return AliFemtoESDTrackCut::IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
  }
}

//________________________________________________________________________________________________________________
bool AliFemtoESDTrackCutNSigmaFilter::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
  if (fUseCustomProtonNSigmaFilter) {
    return fProtonNSigmaFilter->Pass(mom, nsigmaTPCP, nsigmaTOFP);
  } else {
    return AliFemtoESDTrackCut::IsProtonNSigma(mom, nsigmaTPCP, nsigmaTOFP);
  }
}

//________________________________________________________________________________________________________________
bool AliFemtoESDTrackCutNSigmaFilter::IsElectronNSigma(float mom, float nsigmaTPCE, float nsigmaTOFE)
{
  return fElectronNSigmaFilter->Pass(mom, nsigmaTPCE, nsigmaTOFE);
}

//________________________________________________________________________________________________________________
bool AliFemtoESDTrackCutNSigmaFilter::IsProbableElectron(const AliFemtoTrack* track)
{
  float mom, nsigmaTPCPoi, nsigmaTOFPoi, nsigmaTPCE, nsigmaTOFE; //Poi = particle of interest
  mom = track->P().Mag();
  nsigmaTPCE = track->NSigmaTPCE();
  nsigmaTOFE = track->NSigmaTOFE();

  if(fMostProbable==2) {nsigmaTPCPoi = track->NSigmaTPCPi(); nsigmaTOFPoi = track->NSigmaTOFPi();}
  else if(fMostProbable==3) {nsigmaTPCPoi = track->NSigmaTPCK(); nsigmaTOFPoi = track->NSigmaTOFK();}
  else if(fMostProbable==4) {nsigmaTPCPoi = track->NSigmaTPCP(); nsigmaTOFPoi = track->NSigmaTOFP();}
  else{nsigmaTPCPoi = 0.; nsigmaTOFPoi = 0.;}

  if(!fElectronNSigmaFilter->Pass(mom, nsigmaTPCE, nsigmaTOFE)) return false;
  else if((nsigmaTPCPoi > nsigmaTPCE) && (nsigmaTOFPoi > nsigmaTOFE)) return true;
  else return false;
}

