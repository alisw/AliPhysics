/*
***************************************************************************
*
* $Id$
*
*
***************************************************************************
*
*
*
*
***************************************************************************
*
* $Log$
* Revision 1.3  2007/05/22 09:01:42  akisiel
* Add the possibiloity to save cut settings in the ROOT file
*
* Revision 1.2  2007/05/21 10:38:25  akisiel
* More coding rule conformance
*
* Revision 1.1  2007/05/16 10:25:06  akisiel
* Making the directory structure of AliFemtoUser flat. All files go into one common directory
*
* Revision 1.4  2007/05/03 09:46:10  akisiel
* Fixing Effective C++ warnings
*
* Revision 1.3  2007/04/27 07:25:59  akisiel
* Make revisions needed for compilation from the main AliRoot tree
*
* Revision 1.1.1.1  2007/04/25 15:38:41  panos
* Importing the HBT code dir
*
* Revision 1.4  2007-04-03 16:00:08  mchojnacki
* Changes to iprove memory managing
*
* Revision 1.3  2007/03/13 15:30:03  mchojnacki
* adding reader for simulated data
*
* Revision 1.2  2007/03/08 14:58:03  mchojnacki
* adding some alice stuff
*
* Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
* First version on CVS
*
**************************************************************************/

#include "AliFemtoMJTrackCut.h"
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoMJTrackCut)
#endif


// electron
// 0.13 - 1.8
// 0       7.594129e-02    8.256141e-03
// 1       -5.535827e-01   8.170825e-02
// 2       1.728591e+00    3.104210e-01
// 3       -2.827893e+00   5.827802e-01
// 4       2.503553e+00    5.736207e-01
// 5       -1.125965e+00   2.821170e-01
// 6       2.009036e-01    5.438876e-02

// pion
// 0.13 - 2.0
// 0       1.063457e+00    8.872043e-03
// 1       -4.222208e-01   2.534402e-02
// 2       1.042004e-01    1.503945e-02

// kaon
// 0.18 - 2.0
// 0       -7.289406e-02   1.686074e-03
// 1       4.415666e-01    1.143939e-02
// 2       -2.996790e-01   1.840964e-02
// 3       6.704652e-02    7.783990e-03

// proton
// 0.26 - 2.0
// 0       -3.730200e-02   2.347311e-03
// 1       1.163684e-01    1.319316e-02
// 2       8.354116e-02    1.997948e-02
// 3       -4.608098e-02   8.336400e-03


  AliFemtoMJTrackCut::AliFemtoMJTrackCut() :
    fCharge(0),
    fLabel(0),
    fStatus(0),
    fPIDMethod(knSigma),
  fNsigmaTPCTOF(kFALSE),
  fNsigmaTPConly(kFALSE),
  fNsigma(3.),
  fNsigma2(3.),
  fNsigmaRejection(3.),
  fNsigmaAccept(2.),
    fminTPCclsF(0),
    fminTPCncls(0),
    fminITScls(0),
    fMaxITSchiNdof(1000.0),
    fMaxTPCchiNdof(1000.0),
    fMaxSigmaToVertex(1000.0),
    fNTracksPassed(0),
    fNTracksFailed(0),
    fRemoveKinks(kFALSE),
    fRemoveITSFake(kFALSE),
    fMostProbable(0),
    fMaxImpactXY(1000.0),
  fMinImpactXY(-1000.0),
    fMaxImpactZ(1000.0),
    fMaxImpactXYPtOff(1000.0),
    fMaxImpactXYPtNrm(1000.0),
    fMaxImpactXYPtPow(1000.0),
    fMaxImpactZPtOff(1000.0),
    fMaxImpactZPtNrm(1000.0),
    fMaxImpactZPtPow(1000.0),
    fMinPforTOFpid(0.0),
    fMaxPforTOFpid(10000.0),
    fMinPforTPCpid(0.0),
    fMaxPforTPCpid(10000.0),
    fMinPforITSpid(0.0),
   fMaxPforITSpid(10000.0),
   fElectronRejection(0)
{
  // Default constructor
  fNTracksPassed = fNTracksFailed = 0;
  fCharge = 0;  // takes both charges 0
  fPt[0]=0.0;              fPt[1] = 100.0;//100
  fRapidity[0]=-2;       fRapidity[1]=2;//-2 2
  fEta[0]=-2;       fEta[1]=2;//-2 2
  fPidProbElectron[0]=-1;fPidProbElectron[1]=2;
  fPidProbPion[0]=-1;    fPidProbPion[1]=2;
  fPidProbKaon[0]=-1;fPidProbKaon[1]=2;
  fPidProbProton[0]=-1;fPidProbProton[1]=2;
  fPidProbMuon[0]=-1;fPidProbMuon[1]=2;
  for (Int_t i = 0; i < 3; i++)
    fCutClusterRequirementITS[i] = AliESDtrackCuts::kOff;
  fLabel=false;
  fStatus=0;
  fminTPCclsF=0;
  fminITScls=0;
  fPIDMethod=knSigma;
  fNsigmaTPCTOF=kFALSE;
  fNsigmaTPConly=kFALSE;
  fNsigma=3.;
  fNsigma2=3.;
  fNsigmaRejection=3.;
  fNsigmaAccept=2.;
}
//------------------------------
AliFemtoMJTrackCut::~AliFemtoMJTrackCut(){
  /* noop */
}
//------------------------------
bool AliFemtoMJTrackCut::Pass(const AliFemtoTrack* track)
{
  //cout<<"AliFemtoMJTrackCut::Pass"<<endl;

  // test the particle and return
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria
  float tMost[5];

  //cout<<"AliFemtoESD  cut"<<endl;
  //cout<<fPidProbPion[0]<<" < pi ="<<track->PidProbPion()<<" <"<<fPidProbPion[1]<<endl;
  if (fStatus!=0)
    {
      //cout<<" status "<<track->Label()<<" "<<track->Flags()<<" "<<track->TPCnclsF()<<" "<<track->ITSncls()<<endl;
      if ((track->Flags()&fStatus)!=fStatus)
	{
	  //	  cout<<track->Flags()<<" "<<fStatus<<" no go through status"<<endl;
	  return false;
	}

    }
  if (fRemoveKinks) {
    if ((track->KinkIndex(0)) || (track->KinkIndex(1)) || (track->KinkIndex(2)))
      return false;
  }
  if (fRemoveITSFake) {
    if (track->ITSncls() < 0)
      return false;
  }
  if (fminTPCclsF>track->TPCnclsF())
    {
      //cout<<" No go because TPC Number of ClsF"<<fminTPCclsF<< " "<<track->TPCnclsF()<<endl;
      return false;
    }
  if (fminTPCncls>track->TPCncls())
    {
      //cout<<" No go because TPC Number of ClsF"<<fminTPCclsF<< " "<<track->TPCnclsF()<<endl;
      return false;
    }
  if (fminITScls>track->ITSncls())
    {
      //cout<<" No go because ITS Number of Cls"<<fminITScls<< " "<<track->ITSncls()<<endl;
      return false;
    }

  if (fMaxImpactXY < TMath::Abs(track->ImpactD()))
    return false;

  if (fMinImpactXY > TMath::Abs(track->ImpactD()))
    return false;

  if (fMaxImpactZ < TMath::Abs(track->ImpactZ()))
    return false;

  if (fMaxSigmaToVertex < track->SigmaToVertex()) {
    return false;
  }

  if (track->ITSncls() > 0)
    if ((track->ITSchi2()/track->ITSncls()) > fMaxITSchiNdof) {
      return false;
    }

  if (track->TPCncls() > 0)
    if ((track->TPCchi2()/track->TPCncls()) > fMaxTPCchiNdof) {
      return false;
    }
  //ITS cluster requirenments
  for (Int_t i = 0; i < 3; i++)
    if(!CheckITSClusterRequirement(fCutClusterRequirementITS[i], track->HasPointOnITSLayer(i*2), track->HasPointOnITSLayer(i*2+1)))
      return false;

  if (fLabel)
    {
      //cout<<"labels"<<endl;
      if(track->Label()<0)
	{
	  fNTracksFailed++;
	  //   cout<<"No Go Through the cut"<<endl;
	  //  cout<<fLabel<<" Label="<<track->Label()<<endl;
	  return false;
	}
    }
  if (fCharge!=0)
    {
      //cout<<"AliFemtoESD  cut ch "<<endl;
      //cout<<fCharge<<" Charge="<<track->Charge()<<endl;
      if (track->Charge()!= fCharge)
	{
	  fNTracksFailed++;
	  //  cout<<"No Go Through the cut"<<endl;
	  // cout<<fCharge<<" Charge="<<track->Charge()<<endl;
	  return false;
	}
    }




  Bool_t tTPCPidIn = (track->Flags()&AliFemtoTrack::kTPCpid)>0;
  Bool_t tITSPidIn = (track->Flags()&AliFemtoTrack::kITSpid)>0;
  Bool_t tTOFPidIn = (track->Flags()&AliFemtoTrack::kTOFpid)>0;

  if(fMinPforTOFpid > 0 && track->P().Mag() > fMinPforTOFpid &&
     track->P().Mag() < fMaxPforTOFpid && !tTOFPidIn)
    {
      fNTracksFailed++;
      return false;
    }

  if(fMinPforTPCpid > 0 && track->P().Mag() > fMinPforTPCpid &&
     track->P().Mag() < fMaxPforTPCpid && !tTPCPidIn)
    {
      fNTracksFailed++;
      return false;
    }

  if(fMinPforITSpid > 0 && track->P().Mag() > fMinPforITSpid &&
     track->P().Mag() < fMaxPforITSpid && !tITSPidIn)
    {
      fNTracksFailed++;
      return false;
    }


  float tEnergy = ::sqrt(track->P().Mag2()+fMass*fMass);
  float tRapidity = 0;
  if(tEnergy-track->P().z()!=0 && (tEnergy+track->P().z())/(tEnergy-track->P().z())>0)
    tRapidity = 0.5*::log((tEnergy+track->P().z())/(tEnergy-track->P().z()));
  float tPt = ::sqrt((track->P().x())*(track->P().x())+(track->P().y())*(track->P().y()));
  float tEta = track->P().PseudoRapidity();

  if (fMaxImpactXYPtOff < 999.0) {
    if ((fMaxImpactXYPtOff + fMaxImpactXYPtNrm*TMath::Power(tPt, fMaxImpactXYPtPow)) < TMath::Abs(track->ImpactD())) {
      fNTracksFailed++;
      return false;
    }
  }

  if (fMaxImpactZPtOff < 999.0) {
    if ((fMaxImpactZPtOff + fMaxImpactZPtNrm*TMath::Power(tPt, fMaxImpactZPtPow)) < TMath::Abs(track->ImpactZ())) {
      fNTracksFailed++;
      return false;
    }
  }
  
  if ((tRapidity<fRapidity[0])||(tRapidity>fRapidity[1]))
    {
      fNTracksFailed++;
      //cout<<"No Go Through the cut"<<endl;
      //cout<<fRapidity[0]<<" < Rapidity ="<<tRapidity<<" <"<<fRapidity[1]<<endl;
      return false;
    }
  if ((tEta<fEta[0])||(tEta>fEta[1]))
    {
      fNTracksFailed++;
      //cout<<"No Go Through the cut"<<endl;
      //cout<<fEta[0]<<" < Eta ="<<tEta<<" <"<<fEta[1]<<endl;
      return false;
    }
  if ((tPt<fPt[0])||(tPt>fPt[1]))
    {
      fNTracksFailed++;
      //cout<<"No Go Through the cut"<<endl;
      //cout<<fPt[0]<<" < Pt ="<<Pt<<" <"<<fPt[1]<<endl;
      return false;
    }




  //   cout << "Track has pids: "
  //        << track->PidProbElectron() << " "
  //        << track->PidProbMuon() << " "
  //        << track->PidProbPion() << " "
  //        << track->PidProbKaon() << " "
  //        << track->PidProbProton() << " "
  //        << track->PidProbElectron()+track->PidProbMuon()+track->PidProbPion()+track->PidProbKaon()+track->PidProbProton() << endl;


  if ((track->PidProbElectron()<fPidProbElectron[0])||(track->PidProbElectron()>fPidProbElectron[1]))
    {
      fNTracksFailed++;
      //cout<<"No Go Through the cut"<<endl;
      //cout<<fPidProbElectron[0]<<" < e ="<<track->PidProbElectron()<<" <"<<fPidProbElectron[1]<<endl;
      return false;
    }
  if ((track->PidProbPion()<fPidProbPion[0])||(track->PidProbPion()>fPidProbPion[1]))
    {
      fNTracksFailed++;
      //cout<<"No Go Through the cut"<<endl;
      //cout<<fPidProbPion[0]<<" < pi ="<<track->PidProbPion()<<" <"<<fPidProbPion[1]<<endl;
      return false;
    }
  if ((track->PidProbKaon()<fPidProbKaon[0])||(track->PidProbKaon()>fPidProbKaon[1]))
    {
      fNTracksFailed++;
      //cout<<"No Go Through the cut"<<endl;
      //cout<<fPidProbKaon[0]<<" < k ="<<track->PidProbKaon()<<" <"<<fPidProbKaon[1]<<endl;
      return false;
    }
  if ((track->PidProbProton()<fPidProbProton[0])||(track->PidProbProton()>fPidProbProton[1]))
    {
      fNTracksFailed++;
      //cout<<"No Go Through the cut"<<endl;
      //cout<<fPidProbProton[0]<<" < p  ="<<track->PidProbProton()<<" <"<<fPidProbProton[1]<<endl;
      return false;
    }
  if ((track->PidProbMuon()<fPidProbMuon[0])||(track->PidProbMuon()>fPidProbMuon[1]))
    {
      fNTracksFailed++;
      //cout<<"No Go Through the cut"<<endl;
      //cout<<fPidProbMuon[0]<<" <  mi="<<track->PidProbMuon()<<" <"<<fPidProbMuon[1]<<endl;
      return false;
    }

  //****N Sigma Method -- electron rejection****
  if(fElectronRejection) 
    if(!IsElectron(track->NSigmaTPCE(),track->NSigmaTPCPi(),track->NSigmaTPCK(), track->NSigmaTPCP())) 
      return false;


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
	    if (IsKaonNSigma(track->Pt(), track->NSigmaTPCK(), track->NSigmaTOFK()))
	      imost = 11;
	  }
	  else if (fMostProbable == 12) { //cut on Nsigma in pT not p
	    if ( IsProtonNSigma(track->Pt(), track->NSigmaTPCP(), track->NSigmaTOFP()))
	      imost = 12;
	  }
	  else if (fMostProbable == 13) { //cut on Nsigma in pT not p, EXCLUSIVE PID, pion
	    if ((IsPionNSigma(track->Pt(),track->NSigmaTPCPi(), track->NSigmaTOFPi()) && !IsKaonNSigma(track->Pt(),track->NSigmaTPCK(), track->NSigmaTOFK()) && !IsProtonNSigma(track->Pt(),track->NSigmaTPCP(), track->NSigmaTOFP())))
	      imost = 13;
	  }
	  else if (fMostProbable == 14) { //cut on Nsigma in pT not p, EXCLUSIVE PID, kaon
	    if (!IsPionNSigma(track->Pt(),track->NSigmaTPCPi(), track->NSigmaTOFPi()) && IsKaonNSigma(track->Pt(),track->NSigmaTPCK(), track->NSigmaTOFK()) && !IsProtonNSigma(track->Pt(),track->NSigmaTPCP(), track->NSigmaTOFP()))
	      imost = 14;
	  }
	  else if (fMostProbable == 15) { //cut on Nsigma in pT not p, EXCLUSIVE PID, proton
	    if (!IsPionNSigma(track->Pt(),track->NSigmaTPCPi(), track->NSigmaTOFPi()) && !IsKaonNSigma(track->Pt(), !track->NSigmaTPCK(), track->NSigmaTOFK()) && IsProtonNSigma(track->Pt(),track->NSigmaTPCP(), track->NSigmaTOFP()))
	      imost = 15;
	  }
	  //*************** Without double counting, pions ************************
	  else if (fMostProbable == 16) { //cut on Nsigma in pT not p,  without double counting
	    double nSigmaPIDPi = 0, nSigmaPIDK = 0, nSigmaPIDP = 0;
	    if(track->Pt()<0.5){
	      nSigmaPIDPi = abs(track->NSigmaTPCPi());
	      nSigmaPIDK  = abs(track->NSigmaTPCK());
	      nSigmaPIDP  = abs(track->NSigmaTPCP());
	    }
	    else{
	      nSigmaPIDPi = TMath::Hypot(track->NSigmaTOFPi(), track->NSigmaTPCPi());
	      nSigmaPIDK= TMath::Hypot(track->NSigmaTOFK(), track->NSigmaTPCK());
	      nSigmaPIDP= TMath::Hypot(track->NSigmaTOFP(), track->NSigmaTPCP());
	    }

	    if(nSigmaPIDPi<nSigmaPIDK && nSigmaPIDPi<nSigmaPIDP){
	      bool isPionNsigma = 0;
	      isPionNsigma = (IsPionNSigma(track->Pt(),track->NSigmaTPCPi(), track->NSigmaTOFPi()));
	      if(isPionNsigma) imost=16;
	    }
	  }

	  else if (fMostProbable == 17) { //cut on Nsigma in pT not p,  without double counting
	    double nSigmaPIDPi = 0, nSigmaPIDK = 0, nSigmaPIDP = 0;
	    if(track->Pt()<0.5){
	      nSigmaPIDPi = abs(track->NSigmaTPCPi());
	      nSigmaPIDK  = abs(track->NSigmaTPCK());
	      nSigmaPIDP  = abs(track->NSigmaTPCP());
	    }
	    else{
	      nSigmaPIDPi = TMath::Hypot(track->NSigmaTOFPi(), track->NSigmaTPCPi());
	      nSigmaPIDK= TMath::Hypot(track->NSigmaTOFK(), track->NSigmaTPCK());
	      nSigmaPIDP= TMath::Hypot(track->NSigmaTOFP(), track->NSigmaTPCP());
	    }

	    if(nSigmaPIDK<nSigmaPIDPi && nSigmaPIDK<nSigmaPIDP){
	      bool isKaonNsigma = 0;
	      isKaonNsigma = (IsKaonNSigma(track->Pt(),track->NSigmaTPCK(), track->NSigmaTOFK()));
	      if(isKaonNsigma) imost=17;
	    }
	  }
	  else if (fMostProbable == 18) { //cut on Nsigma in pT not p,  without double counting
	    double nSigmaPIDPi = 0, nSigmaPIDK = 0, nSigmaPIDP = 0;
	    if(track->Pt()<0.5){
	      nSigmaPIDPi = abs(track->NSigmaTPCPi());
	      nSigmaPIDK  = abs(track->NSigmaTPCK());
	      nSigmaPIDP  = abs(track->NSigmaTPCP());
	    }
	    else{
	      nSigmaPIDPi = TMath::Hypot(track->NSigmaTOFPi(), track->NSigmaTPCPi());
	      nSigmaPIDK= TMath::Hypot(track->NSigmaTOFK(), track->NSigmaTPCK());
	      nSigmaPIDP= TMath::Hypot(track->NSigmaTOFP(), track->NSigmaTPCP());
	    }

	    if(nSigmaPIDP<nSigmaPIDPi && nSigmaPIDP<nSigmaPIDK){
	      bool isProtonNsigma  = 0;
	      isProtonNsigma = (IsProtonNSigma(track->Pt(),track->NSigmaTPCP(), track->NSigmaTOFP()));
	      if(isProtonNsigma) imost=18;
	    }
	  }
	  else if (fMostProbable == 19) { //cut on Nsigma in pT not p, EXCLUSIVE PID with different REJECTION, pion
	    if ((IsPionNSigmaAccept(track->Pt(),track->NSigmaTPCPi(), track->NSigmaTOFPi(), track->TOFpionTime()) && !IsKaonNSigmaRejection(track->Pt(),track->NSigmaTPCK(), track->NSigmaTOFK(), track->TOFkaonTime()) && !IsProtonNSigmaRejection(track->Pt(),track->NSigmaTPCP(), track->NSigmaTOFP(), track->TOFprotonTime())))
	      imost = 19;
	  }
	  else if (fMostProbable == 20) { //cut on Nsigma in pT not p, EXCLUSIVE PID with different REJECTION, kaon
	    if (!IsPionNSigmaRejection(track->Pt(),track->NSigmaTPCPi(), track->NSigmaTOFPi(), track->TOFpionTime()) && IsKaonNSigmaAccept(track->Pt(),track->NSigmaTPCK(), track->NSigmaTOFK(), track->TOFkaonTime()) && !IsProtonNSigmaRejection(track->Pt(),track->NSigmaTPCP(), track->NSigmaTOFP(), track->TOFprotonTime()))
	      imost = 20;
	  }
	  else if (fMostProbable == 21) { //cut on Nsigma in pT not p, EXCLUSIVE PID with different REJECTION, proton
	    if (!IsPionNSigmaRejection(track->Pt(),track->NSigmaTPCPi(), track->NSigmaTOFPi(), track->TOFpionTime()) && !IsKaonNSigmaRejection(track->Pt(), !track->NSigmaTPCK(), track->NSigmaTOFK(), track->TOFkaonTime()) && IsProtonNSigmaAccept(track->Pt(),track->NSigmaTPCP(), track->NSigmaTOFP(), track->TOFprotonTime()))
	      imost = 21;
	  }
	  else if(fMostProbable == 30) { //DEUTERONs
	    if(IsDeuteronTPCNSigma(track->Pt(),track->TPCsignal()))
	      imost = 30;	    
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
//------------------------------
AliFemtoString AliFemtoMJTrackCut::Report()
{
  // Prepare report from the execution
  string tStemp;
  char tCtemp[100];
  snprintf(tCtemp , 100, "Particle mass:\t%E\n",this->Mass());
  tStemp=tCtemp;
  snprintf(tCtemp , 100, "Particle charge:\t%d\n",fCharge);
  tStemp+=tCtemp;
  snprintf(tCtemp , 100, "Particle pT:\t%E - %E\n",fPt[0],fPt[1]);
  tStemp+=tCtemp;
  snprintf(tCtemp , 100, "Particle rapidity:\t%E - %E\n",fRapidity[0],fRapidity[1]);
  tStemp+=tCtemp;
  snprintf(tCtemp , 100, "Particle eta:\t%E - %E\n",fEta[0],fEta[1]);
  tStemp+=tCtemp;
  snprintf(tCtemp , 100, "Number of tracks which passed:\t%ld  Number which failed:\t%ld\n",fNTracksPassed,fNTracksFailed);
  tStemp += tCtemp;
  AliFemtoString returnThis = tStemp;
  return returnThis;
}
TList *AliFemtoMJTrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoMJTrackCut.mass=%f", this->Mass());
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoMJTrackCut.charge=%i", fCharge);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobpion.minimum=%f", fPidProbPion[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobpion.maximum=%f", fPidProbPion[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobkaon.minimum=%f", fPidProbKaon[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobkaon.maximum=%f", fPidProbKaon[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobproton.minimum=%f", fPidProbProton[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobproton.maximum=%f", fPidProbProton[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobelectron.minimum=%f", fPidProbElectron[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobelectron.maximum=%f", fPidProbElectron[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobMuon.minimum=%f", fPidProbMuon[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pidprobMuon.maximum=%f", fPidProbMuon[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.minimumtpcclusters=%i", fminTPCclsF);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.minimumitsclusters=%i", fminTPCclsF);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pt.minimum=%f", fPt[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.pt.maximum=%f", fPt[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.rapidity.minimum=%f", fRapidity[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.rapidity.maximum=%f", fRapidity[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.removekinks=%i", fRemoveKinks);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.maxitschindof=%f", fMaxITSchiNdof);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.maxtpcchindof=%f", fMaxTPCchiNdof);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.maxsigmatovertex=%f", fMaxSigmaToVertex);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.maximpactxy=%f", fMaxImpactXY);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoMJTrackCut.maximpactz=%f", fMaxImpactZ);
  tListSetttings->AddLast(new TObjString(buf));
  if (fMostProbable) {
    if (fMostProbable == 2)
      snprintf(buf, 200, "AliFemtoMJTrackCut.mostprobable=%s", "Pion");
    if (fMostProbable == 3)
      snprintf(buf, 200, "AliFemtoMJTrackCut.mostprobable=%s", "Kaon");
    if (fMostProbable == 4)
      snprintf(buf, 200, "AliFemtoMJTrackCut.mostprobable=%s", "Proton");
    tListSetttings->AddLast(new TObjString(buf));
  }
  return tListSetttings;
}
void AliFemtoMJTrackCut::SetRemoveKinks(const bool& flag)
{
  fRemoveKinks = flag;
}

void AliFemtoMJTrackCut::SetRemoveITSFake(const bool& flag)
{
  fRemoveITSFake = flag;
}

// electron
// 0.13 - 1.8
// 0       7.594129e-02    8.256141e-03
// 1       -5.535827e-01   8.170825e-02
// 2       1.728591e+00    3.104210e-01
// 3       -2.827893e+00   5.827802e-01
// 4       2.503553e+00    5.736207e-01
// 5       -1.125965e+00   2.821170e-01
// 6       2.009036e-01    5.438876e-02
float AliFemtoMJTrackCut::PidFractionElectron(float mom) const
{
  // Provide a parameterized fraction of electrons dependent on momentum
  if (mom<0.13)
    return (7.594129e-02
	    -5.535827e-01*0.13
	    +1.728591e+00*0.13*0.13
	    -2.827893e+00*0.13*0.13*0.13
	    +2.503553e+00*0.13*0.13*0.13*0.13
	    -1.125965e+00*0.13*0.13*0.13*0.13*0.13
	    +2.009036e-01*0.13*0.13*0.13*0.13*0.13*0.13);

  if (mom>1.8)
    return (7.594129e-02
	    -5.535827e-01*1.8
	    +1.728591e+00*1.8*1.8
	    -2.827893e+00*1.8*1.8*1.8
	    +2.503553e+00*1.8*1.8*1.8*1.8
	    -1.125965e+00*1.8*1.8*1.8*1.8*1.8
	    +2.009036e-01*1.8*1.8*1.8*1.8*1.8*1.8);
  return (7.594129e-02
	  -5.535827e-01*mom
	  +1.728591e+00*mom*mom
	  -2.827893e+00*mom*mom*mom
	  +2.503553e+00*mom*mom*mom*mom
	  -1.125965e+00*mom*mom*mom*mom*mom
	  +2.009036e-01*mom*mom*mom*mom*mom*mom);
}

// pion
// 0.13 - 2.0
// 0       1.063457e+00    8.872043e-03
// 1       -4.222208e-01   2.534402e-02
// 2       1.042004e-01    1.503945e-02
float AliFemtoMJTrackCut::PidFractionPion(float mom) const
{
  // Provide a parameterized fraction of pions dependent on momentum
  if (mom<0.13)
    return ( 1.063457e+00
	     -4.222208e-01*0.13
	     +1.042004e-01*0.0169);
  if (mom>2.0)
    return ( 1.063457e+00
	     -4.222208e-01*2.0
	     +1.042004e-01*4.0);
  return ( 1.063457e+00
	   -4.222208e-01*mom
	   +1.042004e-01*mom*mom);
}

// kaon
// 0.18 - 2.0
// 0       -7.289406e-02   1.686074e-03
// 1       4.415666e-01    1.143939e-02
// 2       -2.996790e-01   1.840964e-02
// 3       6.704652e-02    7.783990e-03
float AliFemtoMJTrackCut::PidFractionKaon(float mom) const
{
  // Provide a parameterized fraction of kaons dependent on momentum
  if (mom<0.18)
    return (-7.289406e-02
	    +4.415666e-01*0.18
	    -2.996790e-01*0.18*0.18
	    +6.704652e-02*0.18*0.18*0.18);
  if (mom>2.0)
    return (-7.289406e-02
	    +4.415666e-01*2.0
	    -2.996790e-01*2.0*2.0
	    +6.704652e-02*2.0*2.0*2.0);
  return (-7.289406e-02
	  +4.415666e-01*mom
	  -2.996790e-01*mom*mom
	  +6.704652e-02*mom*mom*mom);
}

// proton
// 0.26 - 2.0
// 0       -3.730200e-02   2.347311e-03
// 1       1.163684e-01    1.319316e-02
// 2       8.354116e-02    1.997948e-02
// 3       -4.608098e-02   8.336400e-03
float AliFemtoMJTrackCut::PidFractionProton(float mom) const
{
  // Provide a parameterized fraction of protons dependent on momentum
  if (mom<0.26) return  0.0;
  if (mom>2.0)
    return (-3.730200e-02
	    +1.163684e-01*2.0
	    +8.354116e-02*2.0*2.0
	    -4.608098e-02*2.0*2.0*2.0);
  return (-3.730200e-02
	  +1.163684e-01*mom
	  +8.354116e-02*mom*mom
	  -4.608098e-02*mom*mom*mom);
}

void AliFemtoMJTrackCut::SetMomRangeTOFpidIs(const float& minp, const float& maxp)
{
  fMinPforTOFpid = minp;
  fMaxPforTOFpid = maxp;
}

void AliFemtoMJTrackCut::SetMomRangeTPCpidIs(const float& minp, const float& maxp)
{
  fMinPforTPCpid = minp;
  fMaxPforTPCpid = maxp;
}

void AliFemtoMJTrackCut::SetMomRangeITSpidIs(const float& minp, const float& maxp)
{
  fMinPforITSpid = minp;
  fMaxPforITSpid = maxp;
}

bool AliFemtoMJTrackCut::IsPionTPCdEdx(float mom, float dEdx)
{
  //   double a1 = -95.4545, b1 = 86.5455;
  //   double a2 = 0.0,      b2 = 56.0;
  double a1 = -343.75,  b1 = 168.125;
  double a2 = 0.0,      b2 = 65.0;

  if (mom < 0.32) {
    if (dEdx < a1*mom+b1) return true;
  }
  if (dEdx < a2*mom+b2) return true;

  return false;
}

bool AliFemtoMJTrackCut::IsKaonTPCdEdx(float mom, float dEdx)
{

//   double a1 = -547.0; double b1 =  297.0;
//   double a2 = -125.0; double b2 =  145.0;
//   double a3 = -420.0; double b3 =  357.0;
//   double a4 = -110.0; double b4 =  171.0;
//   double b5 =   72.0;

//   if (mom<0.2) return false;

//   if (mom<0.36) {
//     if (dEdx < a1*mom+b1) return false;
//     if (dEdx > a3*mom+b3) return false;
//   }
//   else if (mom<0.6) {
//     if (dEdx < a2*mom+b2) return false;
//     if (dEdx > a3*mom+b3) return false;
//   }
//   else if (mom<0.9) {
//     if (dEdx > a4*mom+b4) return false;
//     if (dEdx <        b5) return false;
//   }
//   else
//     return false;
//   //   else {
//   //     if (dEdx > b5) return false;
//   //   }

//   return true;

  double a1 = -268.896; double b1 =  198.669;
  double a2 = -49.0012;  double b2 =  88.7214;

  if (mom<0.2) return false;

  if (mom>0.3 && mom<0.5) {
    if (dEdx < a1*mom+b1) return false;
  }
  else  if (mom<1.2) {
    if (dEdx < a2*mom+b2) return false;
  }

  return true;

}

bool AliFemtoMJTrackCut::IsProtonTPCdEdx(float mom, float dEdx)
{
  double a1 = -1800.0; double b1 =  940.0;
  double a2 = -500.0;  double b2 =  420.0;
  double a3 = -216.7;  double b3 =  250.0;

  if (mom<0.2) return false;

  if (mom>0.3 && mom<0.4) {
    if (dEdx < a1*mom+b1) return false;
  }
  else  if (mom<0.6) {
    if (dEdx < a2*mom+b2) return false;
  }
  else  if (mom<0.9) {
    if (dEdx < a3*mom+b3) return false;
  }

  return true;

}

bool AliFemtoMJTrackCut::IsPionTOFTime(float mom, float ttof)
{
  double a1 = -427.0; double b1 =  916.0;
  double a2 =  327.0; double b2 = -888.0;
  if (mom<0.3) return kFALSE;
  if (mom>2.0) return kFALSE;
  if (ttof > a1*mom+b1) return kFALSE;
  if (ttof < a2*mom+b2) return kFALSE;

  return kTRUE;
}

bool AliFemtoMJTrackCut::IsKaonTOFTime(float mom, float ttof)
{
  double a1 =   000.0; double b1 =  -500.0;
  double a2 =   000.0; double b2 =   500.0;
  double a3 =   850.0; double b3 = -1503.0;
  double a4 = -1637.0; double b4 =  3621.0;

  if (mom<0.3) return kFALSE;
  if (mom>2.06) return kFALSE;
  if (mom<1.2) {
    if (ttof > a2*mom+b2) return kFALSE;
    if (ttof < a1*mom+b1) return kFALSE;
  }
  if (mom<1.9) {
    if (ttof > a2*mom+b2) return kFALSE;
    if (ttof < a3*mom+b3) return kFALSE;
  }
  if (mom<2.06) {
    if (ttof > a4*mom+b4) return kFALSE;
    if (ttof < a3*mom+b3) return kFALSE;
  }
  return kTRUE;
}

bool AliFemtoMJTrackCut::IsProtonTOFTime(float mom, float ttof)
{
  double a1 =   000.0; double b1 =  -915.0;
  double a2 =   000.0; double b2 =   600.0;
  double a3 =   572.0; double b3 = -1715.0;

  if (mom<0.3) return kFALSE;
  if (mom>3.0) return kFALSE;
  if (mom<1.4) {
    if (ttof > a2*mom+b2) return kFALSE;
    if (ttof < a1*mom+b1) return kFALSE;
  }
  if (mom<3.0) {
    if (ttof > a2*mom+b2) return kFALSE;
    if (ttof < a3*mom+b3) return kFALSE;
  }
  return kTRUE;
}




bool AliFemtoMJTrackCut::IsKaonTPCdEdxNSigma(float mom, float nsigmaK)
{
//  cout<<" AliFemtoMJTrackCut::IsKaonTPCdEdxNSigma "<<mom<<" "<<nsigmaK<<endl;


  if(mom<0.35 && TMath::Abs(nsigmaK)<5.0)return true;
  if(mom>=0.35 && mom<0.5 && TMath::Abs(nsigmaK)<3.0)return true;
  if(mom>=0.5 && mom<0.7 && TMath::Abs(nsigmaK)<2.0)return true;

  return false;
}

bool AliFemtoMJTrackCut::IsKaonTOFNSigma(float mom, float nsigmaK)
{
//  cout<<" AliFemtoMJTrackCut::IsKaonTPCdEdxNSigma "<<mom<<" "<<nsigmaK<<endl;
  //fan
  //  if(mom<1.5 && TMath::Abs(nsigmaK)<3.0)return true;
  if(mom>=1.5 && TMath::Abs(nsigmaK)<2.0)return true;
  return false;
}

/*
bool AliFemtoMJTrackCut::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{


  if(mom<0.5)
    {
	  if(TMath::Abs(nsigmaTPCK)<2.0)
	   {
	   return true;
	   }
	   else
	   {
	   return false;
	   }
    }


   if(mom>=0.5)
    {
         if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0)
         {
         return true;
         }
         else
         {
         return false;
         }
    }

//   if(mom>1.5 || mom<0.15)return false;


}

*/

//old
bool AliFemtoMJTrackCut::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if (fNsigmaTPCTOF) {
    if (mom > 0.5) {
      //        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
      if(mom < 2.0) {
	if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < fNsigma)
	  return true;
      }
      else
	if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < fNsigma2)
	  return true;	
    }
    else {
      if (TMath::Abs(nsigmaTPCK) < fNsigma)
	return true;
    }
  }
  else {

    if(mom<0.4)
      {
	if(nsigmaTOFK<-999.)
	  {
	    if(TMath::Abs(nsigmaTPCK)<2.0) return true;
	  }
	else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
      }
    else if(mom>=0.4 && mom<=0.6)
      {
	if(nsigmaTOFK<-999.)
	  {
	    if(TMath::Abs(nsigmaTPCK)<2.0) return true;
	  }
	else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
      }
    else if(nsigmaTOFK<-999.)
      {
	return false;
      }
    else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
  }
  return false;
}



bool AliFemtoMJTrackCut::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if (fNsigmaTPCTOF) {
    if (mom > 0.5) {
      //        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
      if(mom < 2.0) {
	if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < fNsigma)
	  return true;
      }
      else
	if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < fNsigma2)
	  return true;	
    }
    else {
      if (TMath::Abs(nsigmaTPCPi) < fNsigma)
	return true;
    }
  }
  else {
    if(mom<0.65)
      {
	if(nsigmaTOFPi<-999.)
	  {
	    if(mom<0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
	    else if(mom<0.5 && mom>=0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
	    else if(mom>=0.5 && TMath::Abs(nsigmaTPCPi)<2.0) return true;
	    else return false;
	  }
	else if(TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
      }
    else if(nsigmaTOFPi<-999.)
      {
	return false;
      }
    else if(mom<1.5 && TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
    else if(mom>=1.5 && TMath::Abs(nsigmaTOFPi)<2.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
  }
  return false;
}


bool AliFemtoMJTrackCut::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{

  if (fNsigmaTPCTOF) {
    if (mom > 0.5) {
//        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
      if(mom < 2.0) {
	if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigma)
	  return true;
      }
      else
	if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigma2)
	  return true;	
	}
    else {
        if (TMath::Abs(nsigmaTPCP) < fNsigma)
            return true;
    }
  }
  else if (fNsigmaTPConly) {
    if (TMath::Abs(nsigmaTPCP) < fNsigma)
      return true;
  }
  else {
    if (mom > 0.8 && mom < 2.5) {
      if ( TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 3.0)
        return true;
    }
    else if (mom > 2.5) {
      if ( TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 2.0)
        return true;
    }
    else {
      if (TMath::Abs(nsigmaTPCP) < 3.0)
        return true;
    }
  }

  return false;
}








//rejection methods
bool AliFemtoMJTrackCut::IsPionNSigmaRejection(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime)
{

    if (mom > 0.5) {
        if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < fNsigmaRejection)
            return true;
        }
    else {
        if (TMath::Abs(nsigmaTPCPi) < fNsigmaRejection)
            return true;
    }

  return false;
}


bool AliFemtoMJTrackCut::IsKaonNSigmaRejection(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
  if (mom > 0.5) {
    //rejection of unwanted contamination
    //if(mom>1 && TOFtime<-400)
    //return false;

    if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < fNsigmaRejection)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCK) < fNsigmaRejection)
      return true;
  }

  return false;
}

bool AliFemtoMJTrackCut::IsProtonNSigmaRejection(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{
  if (mom > 0.5) {
    // if(mom>1.8 && TOFtime<-300)
    //return false;


    if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigmaRejection)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCP) < fNsigmaRejection)
      return true;
  }


  return false;
}


bool AliFemtoMJTrackCut::IsPionNSigmaAccept(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime)
{

    if (mom > 0.5) {
        if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < fNsigmaAccept)
            return true;
        }
    else {
        if (TMath::Abs(nsigmaTPCPi) < fNsigmaAccept)
            return true;
    }

  return false;
}


bool AliFemtoMJTrackCut::IsKaonNSigmaAccept(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
  if (mom > 0.5) {
    //rejection of unwanted contamination
    //if(mom>1 && TOFtime<-400)
    //return false;

    if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < fNsigmaAccept)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCK) < fNsigmaAccept)
      return true;
  }

  return false;
}

bool AliFemtoMJTrackCut::IsProtonNSigmaAccept(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{
  if (mom > 0.5) {
    // if(mom>1.8 && TOFtime<-300)
    // return false;


    if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigmaAccept)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCP) < fNsigmaAccept)
      return true;
  }


  return false;
}










void AliFemtoMJTrackCut::SetPIDMethod(ReadPIDMethodType newMethod)
{
  fPIDMethod = newMethod;
}


void AliFemtoMJTrackCut::SetNsigmaTPCTOF(Bool_t nsigma)
{
  fNsigmaTPCTOF = nsigma;
}

void AliFemtoMJTrackCut::SetNsigmaTPConly(Bool_t nsigma)
{
  fNsigmaTPConly = nsigma;
}

void AliFemtoMJTrackCut::SetNsigma(Double_t nsigma)
{
  fNsigma = nsigma;
  fNsigma2 = nsigma;
  fNsigmaAccept = nsigma;
  fNsigmaRejection = nsigma;
}

void AliFemtoMJTrackCut::SetNsigma2(Double_t nsigma)
{
   fNsigma2 = nsigma;
   fNsigmaRejection = nsigma;
}

void AliFemtoMJTrackCut::SetNsigmaRejection(Double_t nsigma)
{
  fNsigmaRejection = nsigma;
}

void AliFemtoMJTrackCut::SetNsigmaAccept(Double_t nsigma)
{
  fNsigmaAccept = nsigma;
}

void AliFemtoMJTrackCut::SetClusterRequirementITS(AliESDtrackCuts::Detector det, AliESDtrackCuts::ITSClusterRequirement req)
{
  fCutClusterRequirementITS[det] = req;
}

Bool_t AliFemtoMJTrackCut::CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2)
{
  // checks if the cluster requirement is fullfilled (in this case: return kTRUE)

  switch (req)
    {
    case AliESDtrackCuts::kOff:        return kTRUE;
    case AliESDtrackCuts::kNone:       return !clusterL1 && !clusterL2;
    case AliESDtrackCuts::kAny:        return clusterL1 || clusterL2;
    case AliESDtrackCuts::kFirst:      return clusterL1;
    case AliESDtrackCuts::kOnlyFirst:  return clusterL1 && !clusterL2;
    case AliESDtrackCuts::kSecond:     return clusterL2;
    case AliESDtrackCuts::kOnlySecond: return clusterL2 && !clusterL1;
    case AliESDtrackCuts::kBoth:       return clusterL1 && clusterL2;
  }

  return kFALSE;
}

bool AliFemtoMJTrackCut::IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
{
  if(TMath::Abs(nsigmaTPCE)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3)
      return false;
   else
     return true;
}

 bool AliFemtoMJTrackCut::IsDeuteronTPCNSigma(float mom, float nsigmaTPC)
 {   
   if(TMath::Abs(nsigmaTPC)<fNsigma)
     return true;
   else
     return false;
 }
