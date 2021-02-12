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
*AliFemtoESDTrackCut
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

#include "AliFemtoESDTrackCutPlusJets25.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
ClassImp(AliFemtoESDTrackCutPlusJets25,1);
  /// \endcond
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


AliFemtoESDTrackCutPlusJets25::AliFemtoESDTrackCutPlusJets25():
    fCharge(0), // takes both charges 0
    fLabel(false),
    fStatus(0),
    fPIDMethod(knSigma),
    fNsigmaTPCTOF(kFALSE),
    fNsigmaTPConly(kFALSE),
    fNsigma(3.),
    fNsigmaMass(-1.0),
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
    fMinPforTOFpid(0.0),
    fMaxPforTOFpid(10000.0),
    fMinPforTPCpid(0.0),
    fMaxPforTPCpid(10000.0),
    fMinPforITSpid(0.0),
    fMaxPforITSpid(10000.0),
    fElectronRejection(0)
{
  // Default constructor
  fPt[0]=0.0;              fPt[1] = 100.0;//100
  fRapidity[0]=-2;       fRapidity[1]=2;//-2 2
  fEta[0]=-2;       fEta[1]=2;//-2 2

  // all Probabilities range from [-1.0, 2.0]
  fPidProbElectron[0] = fPidProbPion[0] = fPidProbKaon[0] = fPidProbProton[0] = fPidProbMuon[0]=-1.0;
  fPidProbElectron[1] = fPidProbPion[1] = fPidProbKaon[1] = fPidProbProton[1] = fPidProbMuon[1]=2.0;

  for (Int_t i = 0; i < 3; i++)
    fCutClusterRequirementITS[i] = AliESDtrackCuts::kOff;

  
}
//------------------------------
AliFemtoESDTrackCutPlusJets25::~AliFemtoESDTrackCutPlusJets25()
{
  /* noop */
}
//------------------------------
bool AliFemtoESDTrackCutPlusJets25::Pass(const AliFemtoTrack* track)
{
  //cout<<"AliFemtoESDTrackCutPlusJets25::Pass"<<endl;

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
        if (IsProtonNSigma(track->P().Mag(), track->NSigmaTPCP(), track->NSigmaTOFP())
         // && (TMath::Abs(track->NSigmaTPCP()) < TMath::Abs(track->NSigmaTPCPi()))
         // && (TMath::Abs(track->NSigmaTPCP()) < TMath::Abs(track->NSigmaTPCK()))
         // && (TMath::Abs(track->NSigmaTOFP()) < TMath::Abs(track->NSigmaTOFPi()))
         // && (TMath::Abs(track->NSigmaTOFP()) < TMath::Abs(track->NSigmaTOFK()))

        // && IsProtonTPCdEdx(track->P().Mag(), track->TPCsignal())
        ) {
          imost = 4;
        }
      }
      else if (fMostProbable == 13) {
        if (IsDeuteronNSigma(track->P().Mag(),track->MassTOF(), fNsigmaMass, track->NSigmaTPCD(), track->NSigmaTOFD()))
          imost = 13;
        if ((track->P().Mag() < 3) &&!(IsDeuteronTPCdEdx(track->P().Mag(), track->TPCsignal())))
          imost = 0;
      }
      else if (fMostProbable == 14) {
          if (IsTritonNSigma(track->P().Mag(), track->NSigmaTPCT(), track->NSigmaTOFT())){
          imost = 14;
        }
      }
      else if (fMostProbable == 15) {
        if ( IsHe3NSigma(track->P().Mag(), track->NSigmaTPCH(), track->NSigmaTOFH())
      )
          imost = 15;
      }
      else if (fMostProbable == 16) {
        if ( IsAlphaNSigma(track->P().Mag(), track->NSigmaTPCA(), track->NSigmaTOFA())
      )
          imost = 16;
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
      else if (fMostProbable == 7) { // pions OR kaons OR protons OR electrons or or or
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
        if (track->NSigmaTOFPi() != -1000 || track->Pt()<0.5){
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
        else if (track->NSigmaTOFPi() != -1000 || track->Pt()<0.5){
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
    if (fPIDMethod==1) {
      float tMost[5];
      tMost[0] = track->PidProbElectron()*PidFractionElectron(track->P().Mag());
      tMost[1] = 0.0;
      tMost[2] = track->PidProbPion()*PidFractionPion(track->P().Mag());
      tMost[3] = track->PidProbKaon()*PidFractionKaon(track->P().Mag());
      tMost[4] = track->PidProbProton()*PidFractionProton(track->P().Mag());

      for (int ip=0; ip<5; ip++) {
        if (tMost[ip] > ipidmax) {
          ipidmax = tMost[ip];
          imost = ip;
        }
      }

      // Looking for pions
      if (fMostProbable == 2) {
        if (imost == 2) {
          // Using the TPC to reject non-pions
          if (!(IsPionTPCdEdx(track->P().Mag(), track->TPCsignal()))) {
            imost = 0;
          }

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
          if (!(IsKaonTPCdEdx(track->P().Mag(), track->TPCsignal()))) {
            imost = 0;
          } else {
            imost = 3;
          }

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
      }





      //--------------
    else if (fMostProbable == 13) {
        //       if (imost == 3) {
        // Using the TPC to reject non-deuterons
        if (track->P().Mag() < 2) {
          if (!(IsDeuteronTPCdEdx(track->P().Mag(), track->TPCsignal()))) {
            imost = 0;
          } else {
            imost = 13;
          }
        }
      }
      //---

    }

    if (imost != fMostProbable) {
      return false;
    }
  }



  ////  ML ///////////////////////////////////////
  
  //Difference between leading particle and other particles

    // SetMonitor(data) should be called before this (Pass) method!

  auto currentMonitor = ((AliFemtoCutMonitorDphiDeta1*)this->PassMonitor(0));

   double ptmax   = 0; 
   double phimax  = 0; 
   double etamax  = 0; 

 /// std::cout<<" before curreentMon fPtmax= "<<ptmax<<" fPhimax= "<<phimax<<" fEtamax= "<<etamax<<std::endl;
 
   if (currentMonitor)
   {
      etamax  = currentMonitor->fEtamax; 
      phimax = currentMonitor->fPhimax;
      ptmax = currentMonitor->fPtmax;
   }

   double tPhi = track->P().Phi();
   double dphi= TMath::Abs(phimax-tPhi);

  // std::cout<<"-------after -- fPtmax= "<<ptmax<<" fPhimax= "<<phimax<<" fEtamax= "<<etamax<<std::endl;
 //  std::cout<<" tPhi= "<<tPhi<< "dphi = "<<dphi<<std::endl;
  // jets phi ~ +/- pi

 //  if(dphi<2.0 || dphi>4.0){
   if(dphi<2.0){
   fNTracksFailed++;
  // std::cout<<" bad ptmax= "<<ptmax<<" dphi "<<dphi<<std::endl;
   return false;
   }

  // std::cout<<" after----if-- fPt<0.1= "<<ptmax<<" dphi "<<dphi<<std::endl;
  /////////////////////////////////////////////////

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
AliFemtoString AliFemtoESDTrackCutPlusJets25::Report()
{
  // Prepare report from the execution
  AliFemtoString report;
  report += Form("Particle mass:\t%E\n",this->Mass());
  report += Form("Particle charge:\t%d\n",fCharge);
  report += Form("Particle pT:\t%E - %E\n",fPt[0],fPt[1]);

  report += Form("Particle rapidity:\t%E - %E\n",fRapidity[0],fRapidity[1]);
  report += Form("Particle eta:\t%E - %E\n",fEta[0],fEta[1]);
  report += Form("Number of tracks which passed:\t%ld  Number which failed:\t%ld\n",fNTracksPassed,fNTracksFailed);

  return report;
}

TList *AliFemtoESDTrackCutPlusJets25::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  tListSetttings->AddVector(
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.mass=%f", this->Mass())),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.charge=%i", fCharge)),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobpion.minimum=%f", fPidProbPion[0])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobpion.maximum=%f", fPidProbPion[1])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobkaon.minimum=%f", fPidProbKaon[0])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobkaon.maximum=%f", fPidProbKaon[1])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobproton.minimum=%f", fPidProbProton[0])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobproton.maximum=%f", fPidProbProton[1])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobelectron.minimum=%f", fPidProbElectron[0])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobelectron.maximum=%f", fPidProbElectron[1])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobMuon.minimum=%f", fPidProbMuon[0])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pidprobMuon.maximum=%f", fPidProbMuon[1])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.minimumtpcclusters=%i", fminTPCclsF)),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.minimumitsclusters=%i", fminTPCclsF)),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pt.minimum=%f", fPt[0])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.pt.maximum=%f", fPt[1])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.rapidity.minimum=%f", fRapidity[0])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.rapidity.maximum=%f", fRapidity[1])),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.removekinks=%i", fRemoveKinks)),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.maxitschindof=%f", fMaxITSchiNdof)),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.maxtpcchindof=%f", fMaxTPCchiNdof)),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.maxsigmatovertex=%f", fMaxSigmaToVertex)),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.maximpactxy=%f", fMaxImpactXY)),
    new TObjString(Form("AliFemtoESDTrackCutPlusJets25.maximpactz=%f", fMaxImpactZ)),
    nullptr
  );

  if (fMostProbable) {
    tListSetttings->Add(
      new TObjString(
        Form("AliFemtoESDTrackCutPlusJets25.mostprobable=%s",
                      ((fMostProbable == 2) ? "Pion"
                     : (fMostProbable == 3) ? "Kaon"
                     : (fMostProbable == 4) ? "Proton"
                     : "??"))));
  }

  return tListSetttings;
}

void AliFemtoESDTrackCutPlusJets25::SetRemoveKinks(const bool& flag)
{
  fRemoveKinks = flag;
}

void AliFemtoESDTrackCutPlusJets25::SetRemoveITSFake(const bool& flag)
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
float AliFemtoESDTrackCutPlusJets25::PidFractionElectron(float mom) const
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
float AliFemtoESDTrackCutPlusJets25::PidFractionPion(float mom) const
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
float AliFemtoESDTrackCutPlusJets25::PidFractionKaon(float mom) const
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
float AliFemtoESDTrackCutPlusJets25::PidFractionProton(float mom) const
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

void AliFemtoESDTrackCutPlusJets25::SetMomRangeTOFpidIs(const float& minp, const float& maxp)
{
  fMinPforTOFpid = minp;
  fMaxPforTOFpid = maxp;
}

void AliFemtoESDTrackCutPlusJets25::SetMomRangeTPCpidIs(const float& minp, const float& maxp)
{
  fMinPforTPCpid = minp;
  fMaxPforTPCpid = maxp;
}

void AliFemtoESDTrackCutPlusJets25::SetMomRangeITSpidIs(const float& minp, const float& maxp)
{
  fMinPforITSpid = minp;
  fMaxPforITSpid = maxp;
}

bool AliFemtoESDTrackCutPlusJets25::IsPionTPCdEdx(float mom, float dEdx)
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

bool AliFemtoESDTrackCutPlusJets25::IsKaonTPCdEdx(float mom, float dEdx)
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

bool AliFemtoESDTrackCutPlusJets25::IsDeuteronTPCdEdx(float mom, float dEdx)
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


  if (!fNsigmaTPCTOF && (fNsigmaMass == -1)) {
    //for selection with only the TPC detector, to remove visible contamination
    if (dEdx < a5*mom+b5) return false;
  }
  else if(fNsigmaTPCTOF && (fNsigmaMass == -1)) {
    //for selection with tpc/tof to finish sample before the region with contamination
    if (mom > 2.2) return false;
  }

  return true;

}

bool AliFemtoESDTrackCutPlusJets25::IsProtonTPCdEdx(float mom, float dEdx)
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

bool AliFemtoESDTrackCutPlusJets25::IsPionTOFTime(float mom, float ttof)
{
  double a1 = -427.0; double b1 =  916.0;
  double a2 =  327.0; double b2 = -888.0;
  if (mom<0.3) return kFALSE;
  if (mom>2.0) return kFALSE;
  if (ttof > a1*mom+b1) return kFALSE;
  if (ttof < a2*mom+b2) return kFALSE;

  return kTRUE;
}

bool AliFemtoESDTrackCutPlusJets25::IsKaonTOFTime(float mom, float ttof)
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

bool AliFemtoESDTrackCutPlusJets25::IsProtonTOFTime(float mom, float ttof)
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



bool AliFemtoESDTrackCutPlusJets25::IsKaonTPCdEdxNSigma(float mom, float nsigmaK)
{
//  cout<<" AliFemtoESDTrackCutPlusJets25::IsKaonTPCdEdxNSigma "<<mom<<" "<<nsigmaK<<endl;


  if (mom<0.35 && TMath::Abs(nsigmaK)<5.0) return true;
  if (mom>=0.35 && mom<0.5 && TMath::Abs(nsigmaK)<3.0) return true;
  if (mom>=0.5 && mom<0.7 && TMath::Abs(nsigmaK)<2.0) return true;

  return false;
}

bool AliFemtoESDTrackCutPlusJets25::IsKaonTOFNSigma(float mom, float nsigmaK)
{
//  cout<<" AliFemtoESDTrackCutPlusJets25::IsKaonTPCdEdxNSigma "<<mom<<" "<<nsigmaK<<endl;
  //fan
  //  if (mom<1.5 && TMath::Abs(nsigmaK)<3.0) return true;
  if (mom>=1.5 && TMath::Abs(nsigmaK)<2.0) return true;
  return false;
}

/*
bool AliFemtoESDTrackCutPlusJets25::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{


  if (mom<0.5)
    {
	  if (TMath::Abs(nsigmaTPCK)<2.0)
	   {
	   return true;
	   }
	   else
	   {
	   return false;
	   }
    }


   if (mom>=0.5)
    {
         if (TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0)
         {
         return true;
         }
         else
         {
         return false;
         }
    }

//   if (mom>1.5 || mom<0.15) return false;


}

*/

//old
bool AliFemtoESDTrackCutPlusJets25::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if (fNsigmaTPCTOF) {
    if (mom > 0.5) {
      //        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
      if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < fNsigma)
        return true;
    }
    else {
      if (TMath::Abs(nsigmaTPCK) < fNsigma)
        return true;
    }
  }
  else {

    if (mom<0.4)
      {
        if (nsigmaTOFK<-999.)
          {
            if (TMath::Abs(nsigmaTPCK)<2.0) return true;
          }
        else if (TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0)
          {
            return true;
          }
      }
    else if (mom>=0.4 && mom<=0.6)
      {
        if (nsigmaTOFK < -999.)
          {
            if (TMath::Abs(nsigmaTPCK)<2.0) return true;
          }
        else if (TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0)
          {
            return true;
          }
      }
    else if (nsigmaTOFK < -999.)
      {
        return false;
      }
    else if (TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
  }
  return false;
}



bool AliFemtoESDTrackCutPlusJets25::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if (fNsigmaTPCTOF) {
    if (mom > 0.5) {
      //        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
      return TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < fNsigma;
    }

    return TMath::Abs(nsigmaTPCPi) < fNsigma;
  }

  if (mom < 0.65) {
    if (nsigmaTOFPi < -999.) {
        if (mom < 0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
        else if (mom<0.5 && mom>=0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
        else if (mom>=0.5 && TMath::Abs(nsigmaTPCPi)<2.0) return true;
        else return false;
    }
    else if (TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<3.0) {
      return true;
    }
  }
  else if (nsigmaTOFPi < -999.) {
    return false;
  }
  else if (mom<1.5 && TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<5.0) {
    return true;
  }
  else if (mom>=1.5 && TMath::Abs(nsigmaTOFPi)<2.0 && TMath::Abs(nsigmaTPCPi)<5.0) {
    return true;
  }

  return false;
}



bool AliFemtoESDTrackCutPlusJets25::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
  if (fNsigmaTPCTOF) {    if (mom > 0.5) {
//        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigma)
            return true;
    } else if (TMath::Abs(nsigmaTPCP) < fNsigma) {
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


/***********************************************************************/

bool AliFemtoESDTrackCutPlusJets25::IsDeuteronNSigma(float mom, float massTOFPDG,float sigmaMass, float nsigmaTPCD, float nsigmaTOFD)
{
  double massPDGD=1.8756;
  if (fNsigmaTPCTOF) {
    if (mom > 1.4) {  //if TOF avaliable: && (nsigmaTOFD != -1000) --> always TOF
      if ((TMath::Abs(nsigmaTPCD) < fNsigma) && (TMath::Abs(nsigmaTOFD) < 2))
        return true;
    }
    else {
      if (TMath::Abs(nsigmaTPCD) < fNsigma)
        return true;
    }
  }
  else{
    if(sigmaMass < 0){
      if (TMath::Abs(nsigmaTPCD) < fNsigma)
	return true;
    }
    else if(sigmaMass > 0){

      //p dependent mass cut 
      double l1, l2;
      if(sigmaMass > 2){
         //left band (4 sigmas)
         l1 = 2.897 - 0.558*mom + 0.177*mom*mom - 0.026*mom*mom*mom;  
         l2 = 3.77 - 0.487*mom + 0.126*mom*mom - 0.014*mom*mom*mom;
      }
      else if(sigmaMass > 1){
         //right band (4 sigmas)
         l1 = 4.5 - 0.42*mom + 0.1*mom*mom;
         l2 = 6.602 -0.981*mom + 0.303*mom*mom;
      }
      else{
         //signal (dist under the mass peak -- 2 sigmas)
         l1 = 4.002 - 0.627*mom + 0.184*mom*mom - 0.02*mom*mom*mom;
         l2 = 4.35 - 0.399*mom + 0.09*mom*mom;
      }

      if ((TMath::Abs(nsigmaTPCD) < fNsigma) && (massTOFPDG > l1) && (massTOFPDG < l2))
	 return true;

    }

  }

  return false;
}


bool AliFemtoESDTrackCutPlusJets25::IsTritonNSigma(float mom, float nsigmaTPCT, float nsigmaTOFT)
{
  if (fNsigmaTPCTOF) {
      return false;
  }
  else {
    if (mom<2 && TMath::Abs(nsigmaTPCT)<fNsigma) return true;
  }
  return false;
}

bool AliFemtoESDTrackCutPlusJets25::IsHe3NSigma(float mom, float nsigmaTPCH, float nsigmaTOFH)
{
  if (fNsigmaTPCTOF) {
      return false;
  }
  else {
    if (mom<3 && TMath::Abs(nsigmaTPCH)<fNsigma) return true;
  }
  return false;
}

bool AliFemtoESDTrackCutPlusJets25::IsAlphaNSigma(float mom, float nsigmaTPCA, float nsigmaTOFA)
{
  if (fNsigmaTPCTOF) {
      return false;
  }
  else {
    if (mom<3 && TMath::Abs(nsigmaTPCA)<fNsigma) return true;
  }
  return false;
}





//
/*********************************************************************/


void AliFemtoESDTrackCutPlusJets25::SetPIDMethod(ReadPIDMethodType newMethod)
{
  fPIDMethod = newMethod;
}

void AliFemtoESDTrackCutPlusJets25::SetNsigmaTPCTOF(Bool_t nsigma)
{
  fNsigmaTPCTOF = nsigma;
}

void AliFemtoESDTrackCutPlusJets25::SetNsigmaTPConly(Bool_t nsigma)
{
  fNsigmaTPConly = nsigma;
}

void AliFemtoESDTrackCutPlusJets25::SetNsigma(Double_t nsigma)
{
  fNsigma = nsigma;
}
void AliFemtoESDTrackCutPlusJets25::SetNsigmaMass(Double_t nsigma)
{
  fNsigmaMass = nsigma;
}


void AliFemtoESDTrackCutPlusJets25::SetClusterRequirementITS(AliESDtrackCuts::Detector det, AliESDtrackCuts::ITSClusterRequirement req)
{
  fCutClusterRequirementITS[det] = req;
}

Bool_t AliFemtoESDTrackCutPlusJets25::CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2)
{
  // checks if the cluster requirement is fullfilled (in this case: return kTRUE)

  switch (req) {
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

bool AliFemtoESDTrackCutPlusJets25::IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
{
  if (TMath::Abs(nsigmaTPCE)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3)
     return false;
  else
     return true;
}

