/***************************************************************************
 *
 * $Id: AliFemtoQATrackCut.cxx 24360 2008-03-10 09:48:27Z akisiel $
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

#include "AliFemtoQATrackCut.h"
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoQATrackCut)
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


AliFemtoQATrackCut::AliFemtoQATrackCut() :
    AliFemtoTrackCut(),
    fCharge(0),
    fLabel(false),
    fStatus(0),
    fminTPCclsF(0),
    fminTPCncls(0),
    fminITScls(0),
    fminTPCchiNdof(0),
    fMaxTPCncls(1000),
    fMaxITSchiNdof(1000.0),
    fMaxTPCchiNdof(1000.0),
    fMaxSigmaToVertex(1000.0),
    fNTracksPassed(0),
    fNTracksFailed(0),
    fRemoveKinks(kFALSE),
    fMostProbable(0),
    fTPCnclsExclusionSwitch(kFALSE),
    fTPCchiNdofExclusionSwitch(kFALSE)
{
  // Default constructor
  fNTracksPassed = fNTracksFailed = 0;
  fCharge = 0;  // takes both charges 0
  fPt[0]=0.0;              fPt[1] = 100.0;//100
  fRapidity[0]=-2;       fRapidity[1]=2;//-2 2
  fPidProbElectron[0]=-1;fPidProbElectron[1]=2;
  fPidProbPion[0]=-1;    fPidProbPion[1]=2;
  fPidProbKaon[0]=-1;fPidProbKaon[1]=2;
  fPidProbProton[0]=-1;fPidProbProton[1]=2;
  fPidProbMuon[0]=-1;fPidProbMuon[1]=2;
  fTPCnclsExclusionSwitch = false;
  fTPCnclsExclusion[0] = 0;
  fTPCnclsExclusion[1] = 1000;
  fTPCchiNdofExclusionSwitch = false;
  fTPCchiNdofExclusion[0] = 0.0;
  fTPCchiNdofExclusion[1] = 1000.0;
}
//------------------------------
AliFemtoQATrackCut::~AliFemtoQATrackCut(){
  /* noop */
}
//------------------------------
bool AliFemtoQATrackCut::Pass(const AliFemtoTrack* track)
{
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
  if (fminTPCclsF>track->TPCnclsF())
    {
      //cout<<" No go because TPC Number of ClsF"<<fminTPCclsF<< " "<<track->TPCnclsF()<<endl;
      return false;
    }

  // TPC number of clusters:
  if (fTPCnclsExclusionSwitch) {
    bool outTPCnclsExclusionZone[2];
      outTPCnclsExclusionZone[0] = false;
      outTPCnclsExclusionZone[1] = false;
    if ( (fminTPCncls > track->TPCncls()) || (fTPCnclsExclusion[0] < track->TPCncls()) ) {
      //cout<<" No go because TPC Number of ClsF"<<fminTPCclsF<< " "<<track->TPCnclsF()<<endl;
      outTPCnclsExclusionZone[0] = true;
    }
    if ( (fMaxTPCncls < track->TPCncls()) || (fTPCnclsExclusion[1] > track->TPCncls()) ) {
      //cout<<" No go because TPC Number of Cls"<<fMaxTPCclsF<< " "<<track->TPCnclsF()<<endl;
      outTPCnclsExclusionZone[1] = true;
    }
    if ( outTPCnclsExclusionZone[0] && outTPCnclsExclusionZone[1] ) { return false; }
  }
  else {
    if (fminTPCncls > track->TPCncls()) {
      //cout<<" No go because TPC Number of ClsF"<<fminTPCclsF<< " "<<track->TPCnclsF()<<endl;
      return false;
    }
    if (fMaxTPCncls < track->TPCncls()) {
      //cout<<" No go because TPC Number of Cls"<<fMaxTPCclsF<< " "<<track->TPCnclsF()<<endl;
      return false;
    }
  }

  if (fminITScls>track->ITSncls())
    {
      //cout<<" No go because ITS Number of Cls"<<fminITScls<< " "<<track->ITSncls()<<endl;
      return false;
    }

  if (fMaxSigmaToVertex < track->SigmaToVertex()) {
    return false;
  }

  if (track->ITSncls() > 0 && (track->ITSchi2()/track->ITSncls()) > fMaxITSchiNdof) {
    return false;
  }

  const double tpc_chi2 = track->TPCchi2perNDF();

  // TPC chiNdof of tracks:
  if (fTPCchiNdofExclusionSwitch && (track->TPCncls() > 0)) {
    bool outTPCchiNdofExclusionZone[2];
      outTPCchiNdofExclusionZone[0] = false;
      outTPCchiNdofExclusionZone[1] = false;
    if ( (fminTPCchiNdof > tpc_chi2) || (fTPCchiNdofExclusion[0] < tpc_chi2) ) {
      //cout<<" No go because TPC Number of ClsF"<<fminTPCclsF<< " "<<track->TPCnclsF()<<endl;
      outTPCchiNdofExclusionZone[0] = true;
    }
    if ( (fMaxTPCchiNdof < tpc_chi2) || (fTPCchiNdofExclusion[1] > tpc_chi2) ) {
      //cout<<" No go because TPC Number of Cls"<<fMaxTPCclsF<< " "<<track->TPCnclsF()<<endl;
      outTPCchiNdofExclusionZone[1] = true;
    }
    if ( outTPCchiNdofExclusionZone[0] && outTPCchiNdofExclusionZone[1] ) {
      return false;
    }
  }
  else {
    if (fminTPCchiNdof > tpc_chi2) {
      //cout<<" No go because TPC Number of ClsF"<<fminTPCclsF<< " "<<track->TPCnclsF()<<endl;
      return false;
    }
    if (fMaxTPCchiNdof < tpc_chi2) {
      //cout<<" No go because TPC Number of Cls"<<fMaxTPCclsF<< " "<<track->TPCnclsF()<<endl;
      return false;
    }
  }

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
  float tEnergy = ::sqrt(track->P().Mag2()+fMass*fMass);
  float tRapidity = 0.5*::log((tEnergy+track->P().z())/(tEnergy-track->P().z()));
  float tPt = ::sqrt((track->P().x())*(track->P().x())+(track->P().y())*(track->P().y()));
  if ((tRapidity<fRapidity[0])||(tRapidity>fRapidity[1]))
    {
      fNTracksFailed++;
      //cout<<"No Go Through the cut"<<endl;
      //cout<<fRapidity[0]<<" < Rapidity ="<<tRapidity<<" <"<<fRapidity[1]<<endl;
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

  if (fMostProbable) {
    tMost[0] = track->PidProbElectron()*PidFractionElectron(track->P().Mag());
    tMost[1] = 0.0;
    tMost[2] = track->PidProbPion()*PidFractionPion(track->P().Mag());
    tMost[3] = track->PidProbKaon()*PidFractionKaon(track->P().Mag());
    tMost[4] = track->PidProbProton()*PidFractionProton(track->P().Mag());
    int imost=0;
    float ipidmax = 0.0;
    for (int ip=0; ip<5; ip++)
      if (tMost[ip] > ipidmax) { ipidmax = tMost[ip]; imost = ip; };
    if (imost != fMostProbable) return false;
  }

  // cout<<"Go Through the cut"<<endl;
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
AliFemtoString AliFemtoQATrackCut::Report()
{
  // Prepare report from the execution
  AliFemtoString report;
  report += Form("Particle mass:\t%E\n",this->Mass());
  report += Form("Particle charge:\t%d\n",fCharge);
  report += Form("Particle pT:\t%E - %E\n",fPt[0],fPt[1]);
  report += Form("Particle rapidity:\t%E - %E\n",fRapidity[0],fRapidity[1]);
  report += Form("Number of tracks which passed:\t%ld  Number which failed:\t%ld\n",fNTracksPassed,fNTracksFailed);
  return report;
}

TList *AliFemtoQATrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  tListSetttings->AddVector(
    new TObjString(Form("AliFemtoQATrackCut.mass=%f", this->Mass())),
    new TObjString(Form("AliFemtoQATrackCut.charge=%i", fCharge)),
    new TObjString(Form("AliFemtoQATrackCut.pidprobpion.minimum=%f", fPidProbPion[0])),
    new TObjString(Form("AliFemtoQATrackCut.pidprobpion.maximum=%f", fPidProbPion[1])),
    new TObjString(Form("AliFemtoQATrackCut.pidprobkaon.minimum=%f", fPidProbKaon[0])),
    new TObjString(Form("AliFemtoQATrackCut.pidprobkaon.maximum=%f", fPidProbKaon[1])),
    new TObjString(Form("AliFemtoQATrackCut.pidprobproton.minimum=%f", fPidProbProton[0])),
    new TObjString(Form("AliFemtoQATrackCut.pidprobproton.maximum=%f", fPidProbProton[1])),
    new TObjString(Form("AliFemtoQATrackCut.pidprobelectron.minimum=%f", fPidProbElectron[0])),
    new TObjString(Form("AliFemtoQATrackCut.pidprobelectron.maximum=%f", fPidProbElectron[1])),
    new TObjString(Form("AliFemtoQATrackCut.pidprobMuon.minimum=%f", fPidProbMuon[0])),
    new TObjString(Form("AliFemtoQATrackCut.pidprobMuon.maximum=%f", fPidProbMuon[1])),
    new TObjString(Form("AliFemtoQATrackCut.minimumtpcclusters=%i", fminTPCclsF)),
    new TObjString(Form("AliFemtoQATrackCut.minimumitsclusters=%i", fminITScls)),
    new TObjString(Form("AliFemtoQATrackCut.pt.minimum=%f", fPt[0])),
    new TObjString(Form("AliFemtoQATrackCut.pt.maximum=%f", fPt[1])),
    new TObjString(Form("AliFemtoQATrackCut.rapidity.minimum=%f", fRapidity[0])),
    new TObjString(Form("AliFemtoQATrackCut.rapidity.maximum=%f", fRapidity[1])),
    new TObjString(Form("AliFemtoQATrackCut.removekinks=%i", fRemoveKinks)),
    new TObjString(Form("AliFemtoQATrackCut.maxitschindof=%f", fMaxITSchiNdof)),
    new TObjString(Form("AliFemtoQATrackCut.maxtpcchindof=%f", fMaxTPCchiNdof)),
    new TObjString(Form("AliFemtoQATrackCut.maxsigmatovertex=%f", fMaxSigmaToVertex)),
    nullptr
  );

  if (fMostProbable) {
    if (fMostProbable == 2)
      tListSetttings->AddLast(new TObjString("AliFemtoQATrackCut.mostprobable=Pion"));
    if (fMostProbable == 3)
      tListSetttings->AddLast(new TObjString("AliFemtoQATrackCut.mostprobable=Kaon"));
    if (fMostProbable == 4)
      tListSetttings->AddLast(new TObjString("AliFemtoQATrackCut.mostprobable=Proton"));
  }
  return tListSetttings;
}
void AliFemtoQATrackCut::SetRemoveKinks(const bool& flag)
{
  fRemoveKinks = flag;
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
float AliFemtoQATrackCut::PidFractionElectron(float mom) const
{
  // Provide a parameterized fraction of electrons dependent on momentum
  if (mom<0.13) return 0.0;
  if (mom>1.8) return 0.0;
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
float AliFemtoQATrackCut::PidFractionPion(float mom) const
{
  // Provide a parameterized fraction of pions dependent on momentum
  if (mom<0.13) return 0.0;
  if (mom>2.0) return 0.0;
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
float AliFemtoQATrackCut::PidFractionKaon(float mom) const
{
  // Provide a parameterized fraction of kaons dependent on momentum
  if (mom<0.18) return 0.0;
  if (mom>2.0) return 0.0;
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
float AliFemtoQATrackCut::PidFractionProton(float mom) const
{
  // Provide a parameterized fraction of protons dependent on momentum
  if (mom<0.26) return  0.0;
  if (mom>2.0) return 0.0;
  return (-3.730200e-02
	  +1.163684e-01*mom
	  +8.354116e-02*mom*mom
	  -4.608098e-02*mom*mom*mom);
}
