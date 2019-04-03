/*
***************************************************************************
*
* $Id: AliFemtoKKTrackCut.cxx 60781 2013-02-08 16:24:10Z akisiel $ 
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

#include "AliFemtoKKTrackCutFull.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoKKTrackCutFull)
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


  AliFemtoKKTrackCutFull::AliFemtoKKTrackCutFull() :
    fCharge(0),
    fLabel(0),
    fStatus(0),
    fPIDMethod(knSigma),
//ml
  fNsigmaTPCle250(3.),
  fNsigmaTPC250_400(3.),
  fNsigmaTPC400_450(3.),
  fNsigmaTPC450_500(3.),
  fNsigmaTPCge500(3.),
  fNsigmaTOF450_500(3.),
  fUseNsigmaTOF450_500(false),
  fNsigmaTOF500_800(3.),
  fNsigmaTOF800_1000(3.),
  fNsigmaTOFge1000(3.),
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
    fMaxImpactZ(1000.0),
    fMaxImpactXYPtOff(1000.0),
    fMaxImpactXYPtNrm(1000.0),
    fMaxImpactXYPtPow(1000.0),
    fMinPforTOFpid(0.0),
    fMaxPforTOFpid(10000.0),
    fMinPforTPCpid(0.0),
    fMaxPforTPCpid(10000.0),
    fMinPforITSpid(0.0),
    fMaxPforITSpid(10000.0)
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

  fNsigmaTPCle250=3.0;
  fNsigmaTPC250_400=3.0;
  fNsigmaTPC400_450=3.0;
  fNsigmaTPC450_500=3.0;
  fNsigmaTPCge500=3.0;

  fNsigmaTOF450_500=3.0;
  fUseNsigmaTOF450_500=false;
  fNsigmaTOF500_800=3.0;
  fNsigmaTOF800_1000=3.0;
  fNsigmaTOFge1000=3.0;



}
//------------------------------
AliFemtoKKTrackCutFull::~AliFemtoKKTrackCutFull(){
  /* noop */
}
//------------------------------
bool AliFemtoKKTrackCutFull::Pass(const AliFemtoTrack* track)
{
  //cout<<"AliFemtoKKTrackCut::Pass"<<endl;

  // test the particle and return 
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria
  float tMost[5];
  
  //cout<<"AliKK  cut"<<endl;
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
      //cout<<"AliKK  cut ch "<<endl;
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
  float tRapidity = 0.5*::log((tEnergy+track->P().z())/(tEnergy-track->P().z()));
  float tPt = ::sqrt((track->P().x())*(track->P().x())+(track->P().y())*(track->P().y()));
  float tEta = track->P().PseudoRapidity();
  
  if (fMaxImpactXYPtOff < 999.0) {
    if ((fMaxImpactXYPtOff + fMaxImpactXYPtNrm*TMath::Power(tPt, fMaxImpactXYPtPow)) < TMath::Abs(track->ImpactD())) {
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
	  else if (fMostProbable == 4) { // proton nsigma-PID required contour adjusting
	    if (IsProtonNSigma(track->P().Mag(), track->NSigmaTPCP(), track->NSigmaTOFP()) // && IsProtonTPCdEdx(track->P().Mag(), track->TPCsignal())
            )
	      imost = 4;
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
AliFemtoString AliFemtoKKTrackCutFull::Report()
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
TList *AliFemtoKKTrackCutFull::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.mass=%f", this->Mass());
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoKKTrackCutFull.charge=%i", fCharge);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobpion.minimum=%f", fPidProbPion[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobpion.maximum=%f", fPidProbPion[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobkaon.minimum=%f", fPidProbKaon[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobkaon.maximum=%f", fPidProbKaon[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobproton.minimum=%f", fPidProbProton[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobproton.maximum=%f", fPidProbProton[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobelectron.minimum=%f", fPidProbElectron[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobelectron.maximum=%f", fPidProbElectron[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobMuon.minimum=%f", fPidProbMuon[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pidprobMuon.maximum=%f", fPidProbMuon[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.minimumtpcclusters=%i", fminTPCclsF);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.minimumitsclusters=%i", fminTPCclsF);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pt.minimum=%f", fPt[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.pt.maximum=%f", fPt[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.rapidity.minimum=%f", fRapidity[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.rapidity.maximum=%f", fRapidity[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.removekinks=%i", fRemoveKinks);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.maxitschindof=%f", fMaxITSchiNdof);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.maxtpcchindof=%f", fMaxTPCchiNdof);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.maxsigmatovertex=%f", fMaxSigmaToVertex);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.maximpactxy=%f", fMaxImpactXY);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoKKTrackCutFull.maximpactz=%f", fMaxImpactZ);
  tListSetttings->AddLast(new TObjString(buf));
  if (fMostProbable) {
    if (fMostProbable == 2)
      snprintf(buf, 200, "AliFemtoKKTrackCutFull.mostprobable=%s", "Pion");
    if (fMostProbable == 3)
      snprintf(buf, 200, "AliFemtoKKTrackCutFull.mostprobable=%s", "Kaon");
    if (fMostProbable == 4)
      snprintf(buf, 200, "AliFemtoKKTrackCutFull.mostprobable=%s", "Proton");
    tListSetttings->AddLast(new TObjString(buf));
  }
  return tListSetttings;
}
void AliFemtoKKTrackCutFull::SetRemoveKinks(const bool& flag)
{
  fRemoveKinks = flag;
}
			    
void AliFemtoKKTrackCutFull::SetRemoveITSFake(const bool& flag)
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
float AliFemtoKKTrackCutFull::PidFractionElectron(float mom) const
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
float AliFemtoKKTrackCutFull::PidFractionPion(float mom) const
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
float AliFemtoKKTrackCutFull::PidFractionKaon(float mom) const
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
float AliFemtoKKTrackCutFull::PidFractionProton(float mom) const
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

void AliFemtoKKTrackCutFull::SetMomRangeTOFpidIs(const float& minp, const float& maxp)
{
  fMinPforTOFpid = minp;
  fMaxPforTOFpid = maxp;
}

void AliFemtoKKTrackCutFull::SetMomRangeTPCpidIs(const float& minp, const float& maxp)
{
  fMinPforTPCpid = minp;
  fMaxPforTPCpid = maxp;
}

void AliFemtoKKTrackCutFull::SetMomRangeITSpidIs(const float& minp, const float& maxp)
{
  fMinPforITSpid = minp;
  fMaxPforITSpid = maxp;
}

bool AliFemtoKKTrackCutFull::IsPionTPCdEdx(float mom, float dEdx)
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

bool AliFemtoKKTrackCutFull::IsKaonTPCdEdx(float mom, float dEdx)
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

bool AliFemtoKKTrackCutFull::IsProtonTPCdEdx(float mom, float dEdx)
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

bool AliFemtoKKTrackCutFull::IsPionTOFTime(float mom, float ttof)
{
  double a1 = -427.0; double b1 =  916.0;
  double a2 =  327.0; double b2 = -888.0;
  if (mom<0.3) return kFALSE;
  if (mom>2.0) return kFALSE;
  if (ttof > a1*mom+b1) return kFALSE;
  if (ttof < a2*mom+b2) return kFALSE;

  return kTRUE;
}

bool AliFemtoKKTrackCutFull::IsKaonTOFTime(float mom, float ttof)
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

bool AliFemtoKKTrackCutFull::IsProtonTOFTime(float mom, float ttof)
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




bool AliFemtoKKTrackCutFull::IsKaonTPCdEdxNSigma(float mom, float nsigmaK)
{
//  cout<<" AliFemtoKKTrackCut::IsKaonTPCdEdxNSigma "<<mom<<" "<<nsigmaK<<endl;


  if(mom<0.35 && TMath::Abs(nsigmaK)<5.0)return true;
  if(mom>=0.35 && mom<0.5 && TMath::Abs(nsigmaK)<3.0)return true; 
  if(mom>=0.5 && mom<0.7 && TMath::Abs(nsigmaK)<2.0)return true;

  return false;
}

bool AliFemtoKKTrackCutFull::IsKaonTOFNSigma(float mom, float nsigmaK)
{
//  cout<<" AliFemtoKKTrackCut::IsKaonTPCdEdxNSigma "<<mom<<" "<<nsigmaK<<endl;
  //fan
  //  if(mom<1.5 && TMath::Abs(nsigmaK)<3.0)return true;
  if(mom>=1.5 && TMath::Abs(nsigmaK)<2.0)return true; 
  return false;
}


bool AliFemtoKKTrackCutFull::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{

  if(mom<0.25)
    {
      if(TMath::Abs(nsigmaTPCK)<fNsigmaTPCle250)
	{ 
	  return true;
	} 
      else 
	{
	  return false;
	}
    }

  if(mom>=0.25 && mom<0.4)
    {
      if(TMath::Abs(nsigmaTPCK)<fNsigmaTPC250_400)
	{ 
	  return true;
	} 
      else 
	{
	  return false;
	}
    }

  if(mom>=0.4 && mom<0.45)
    {
      if(TMath::Abs(nsigmaTPCK)<fNsigmaTPC400_450)
	{ 
	  return true;
	} 
      else 
	{
	  return false;
	}
    }

 if(mom>=0.45 && mom<0.5 && !fUseNsigmaTOF450_500)
    {
      if(TMath::Abs(nsigmaTPCK)<fNsigmaTPC450_500)
	{ 
	  return true;
	} 
      else 
	{
	  return false;
	}
    }

/////////////  
  
 if(fUseNsigmaTOF450_500)
   {
     if(mom>=0.45 && mom<0.5)
       {
	 if(TMath::Abs(nsigmaTOFK)<fNsigmaTOF450_500 && TMath::Abs(nsigmaTPCK)<fNsigmaTPC450_500) 
	   {
	     return true;
	   }
	 else
	   {
	     return false;
	   }
       }
   }
  
  if(mom>=0.5 && mom<0.8)
    {
      if(TMath::Abs(nsigmaTOFK)<fNsigmaTOF500_800 && TMath::Abs(nsigmaTPCK)<fNsigmaTPCge500) 
	{
	  //     cout<<"500-800 "<<fNsigmaTOF500_800<<" "<<nsigmaTOFK<<endl;
	  return true;
	}
      else
	{
	  return false;
	}
    }

  if(mom>=0.8 && mom<1.0)
    {
      if(TMath::Abs(nsigmaTOFK)<fNsigmaTOF800_1000 && TMath::Abs(nsigmaTPCK)<fNsigmaTPCge500) 
	{
	  //          cout<<"800-1000 "<<fNsigmaTOF800_1000<<" "<<nsigmaTOFK<<endl;
	  return true;
	}
      else
	{
	  return false;
	}
    }

 if(mom>=1.0)
    {
      if(TMath::Abs(nsigmaTOFK)<fNsigmaTOFge1000 && TMath::Abs(nsigmaTPCK)<fNsigmaTPCge500) 
	{
	  //      cout<<">1000 "<<fNsigmaTOFge1000<<" "<<nsigmaTOFK<<endl;
	  return true;
	}
      else
	{
	  return false;
	}
    }
    
    
    
//   if(mom>1.5 || mom<0.15)return false;
  return false;

}



/* old
bool AliFemtoKKTrackCut::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{

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

  return false;
}
*/


bool AliFemtoKKTrackCutFull::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
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

  return false;
}


bool AliFemtoKKTrackCutFull::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
    if (mom > 0.8) {
        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
            return true;
	}
    else {
        if (TMath::Abs(nsigmaTPCP) < 3.0)
            return true;
    }

  return false;
}


void AliFemtoKKTrackCutFull::SetPIDMethod(ReadPIDMethodType newMethod)
{
  fPIDMethod = newMethod;
}


void AliFemtoKKTrackCutFull::SetClusterRequirementITS(AliESDtrackCuts::Detector det, AliESDtrackCuts::ITSClusterRequirement req) 
{ 
  fCutClusterRequirementITS[det] = req; 
}

Bool_t AliFemtoKKTrackCutFull::CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2)
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


void AliFemtoKKTrackCutFull::SetNsigmaTPCle250(Double_t nsigma)
{
  fNsigmaTPCle250 = nsigma;
}

void AliFemtoKKTrackCutFull::SetNsigmaTPC250_400(Double_t nsigma)
{
  fNsigmaTPC250_400 = nsigma;
}

void AliFemtoKKTrackCutFull::SetNsigmaTPC400_450(Double_t nsigma)
{
  fNsigmaTPC400_450 = nsigma;
}

void AliFemtoKKTrackCutFull::SetNsigmaTPC450_500(Double_t nsigma)
{
  fNsigmaTPC450_500 = nsigma;
}

void AliFemtoKKTrackCutFull::SetNsigmaTPCge500(Double_t nsigma)
{
  fNsigmaTPCge500 = nsigma;
}

void AliFemtoKKTrackCutFull::SetNsigmaTOF450_500(Double_t nsigma)
{
  fNsigmaTOF450_500 = nsigma;
}

void AliFemtoKKTrackCutFull::SetNsigmaTOF500_800(Double_t nsigma)
{
  fNsigmaTOF500_800 = nsigma;
}


void AliFemtoKKTrackCutFull::SetNsigmaTOF800_1000(Double_t nsigma)
{
  fNsigmaTOF800_1000 = nsigma;
}

void AliFemtoKKTrackCutFull::SetNsigmaTOFge1000(Double_t nsigma)
{
  fNsigmaTOFge1000 = nsigma;
}

void AliFemtoKKTrackCutFull::UseNsigmaTOF450_500(bool usensigma)
{
  fUseNsigmaTOF450_500 = usensigma;
}
