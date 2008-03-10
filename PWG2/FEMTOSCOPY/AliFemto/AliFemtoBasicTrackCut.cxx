////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoBasicTrackCut - the basic cut for tracks.                          //
// Cuts on particle identification, transverse momentum, rapidity, distance   //
// of closest approach to primary vertex and charge                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoBasicTrackCut.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoBasicTrackCut)
#endif

AliFemtoBasicTrackCut::AliFemtoBasicTrackCut():
  fCharge(0),
  fNTracksPassed(0),
  fNTracksFailed(0)
{
  // Default constructor
  fNTracksPassed = fNTracksFailed = 0;
  fCharge = 1;  // takes both charges 0
  fNSigmaPion[0] = -100.0;   fNSigmaPion[1] = 100.0;
  fNSigmaKaon[0] = -100.0;   fNSigmaKaon[1] = 100.0;
  fNSigmaProton[0] = -100.0; fNSigmaProton[1] = 100.0;
  fNHits[0] = 10;          fNHits[1] = 180;
  fPt[0]=0.0;              fPt[1] = 100.0;//100
  fRapidity[0]=-2;       fRapidity[1]=2;//-2 2
  fDCA[0] = -1.0;           fDCA[1] = 20.0;

}
//------------------------------
//AliFemtoBasicTrackCut::~AliFemtoBasicTrackCut(){
//  /* noop */
//}
//------------------------------
bool AliFemtoBasicTrackCut::Pass(const AliFemtoTrack* track){
  // test the particle and return 
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria

  //  return true ;  // THIS CUT IS A STHBTDUMMY!!

  /*
    cout << endl;
    cout << "#track " << trackCount++;
    cout << " * pion " << (track->NSigmaPion() > fNSigmaPion[0]) && (track->NSigmaPion() < fNSigmaPion[1]);
    cout << " * kaon " << (track->NSigmaKaon() > fNSigmaKaon[0]) && (track->NSigmaKaon() < fNSigmaKaon[1]);
    cout << " * proton " << (track->NSigmaProton() > fNSigmaProton[0]) && (track->NSigmaProton() < fNSigmaProton[1]);
    cout << " * charge " << (track->Charge() == fCharge);
  */
  bool goodPID = 1;  
  /* ----- NOT DOING PID CUTS !!!! ------
  bool goodPID = ((track->NSigmaPion()   > fNSigmaPion[0]) &&
                  (track->NSigmaPion()   < fNSigmaPion[1]) &&
                  (track->NSigmaKaon()   > fNSigmaKaon[0]) &&
                  (track->NSigmaKaon()   < fNSigmaKaon[1]) &&
                  (track->NSigmaProton() > fNSigmaProton[0]) &&
                  (track->NSigmaProton() < fNSigmaProton[1]));
  ----- NOT DOING PID CUTS !!!! ------ */
  if (fCharge !=0){               // if user requests "charge=0" then that means ignore charge
    goodPID = (goodPID&&(track->Charge() == fCharge));
  }
  if (goodPID){
    float tEnergy = ::sqrt(track->P().mag2()+fMass*fMass);
    float tRapidity = 0.5*::log((tEnergy+track->P().z())/
			    (tEnergy-track->P().z()));

    float tPt = ::sqrt((track->P().x())*(track->P().x())+
                    (track->P().y())*(track->P().y()));

    
    /*
      cout << " * DCAxy " << (track->DCAxy()  > fDCA[0]) && (track->DCAxy()  < fDCA[1]);
      cout << " * fDCA[0] " << fDCA[0];
      cout << " * fDCA[1] " << fDCA[1];
      cout << " * track->DCAxy " << track->DCAxy();
      cout << " * NHits " <<  (track->NHits() > fNHits[0]) && (track->NHits() < fNHits[1]); 
      cout << " * tPt " << (tPt > fPt[0]) && (tPt < fPt[1]);
      cout << " * y " << (tRapidity > fRapidity[0]) && (tRapidity < fRapidity[1]);
      cout << endl;
    */

    bool goodTrack=
      (//(track->DCAxy()  > fDCA[0]) &&
     //  (track->DCAxy()  < fDCA[1]) &&
  //     (track->NHits() > fNHits[0]) &&
    //   (track->NHits() < fNHits[1]) &&
       (tPt             > fPt[0]) &&
       (tPt             < fPt[1]) &&
       (tRapidity      > fRapidity[0]) &&
       (tRapidity      < fRapidity[1]));
    //  &&
    //       (track->PidProbPion()>0.5)&&//moje
    //       (track->PidProbMuon()<0.47)&&//moje
    //       (track->Label()>0);//moje

    //    cout << track->DCAxy() << " " << track->NHits() << " " << Pt << " " << tRapidity << " " << tEnergy << endl;

    goodTrack ? fNTracksPassed++ : fNTracksFailed++;
    return (goodTrack);
  }
  else{
    fNTracksFailed++;
    return (goodPID);
  }
}
//------------------------------
AliFemtoString AliFemtoBasicTrackCut::Report(){
  // construct report
  string tStemp;
  char tCtemp[100];
  sprintf(tCtemp,"Particle mass:\t%E\n",this->Mass());
  tStemp=tCtemp;
  sprintf(tCtemp,"Particle charge:\t%d\n",fCharge);
  tStemp=tCtemp;
  sprintf(tCtemp,"Particle Nsigma from pion:\t%E - %E\n",fNSigmaPion[0],fNSigmaPion[1]);
  tStemp+=tCtemp;
  sprintf(tCtemp,"Particle Nsigma from kaon:\t%E - %E\n",fNSigmaKaon[0],fNSigmaKaon[1]);
  tStemp+=tCtemp;
  sprintf(tCtemp,"Particle Nsigma from proton:\t%E - %E\n",fNSigmaProton[0],fNSigmaProton[1]);
  tStemp+=tCtemp;
  sprintf(tCtemp,"Particle #hits:\t%d - %d\n",fNHits[0],fNHits[1]);
  tStemp+=tCtemp;
  sprintf(tCtemp,"Particle pT:\t%E - %E\n",fPt[0],fPt[1]);
  tStemp+=tCtemp;
  sprintf(tCtemp,"Particle rapidity:\t%E - %E\n",fRapidity[0],fRapidity[1]);
  tStemp+=tCtemp;
  sprintf(tCtemp,"Particle DCA:\t%E - %E\n",fDCA[0],fDCA[1]);
  tStemp+=tCtemp;
  sprintf(tCtemp,"Number of tracks which passed:\t%ld  Number which failed:\t%ld\n",fNTracksPassed,fNTracksFailed);
  tStemp += tCtemp;
  AliFemtoString returnThis = tStemp;
  return returnThis;
}

TList *AliFemtoBasicTrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoBasicTrackCut.mass=%f", this->Mass());
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoBasicTrackCut.charge=%i", fCharge);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.nsigmapion.minimum=%f", fNSigmaPion[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.nsigmapion.maximum=%f", fNSigmaPion[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.nsigmakaon.minimum=%f", fNSigmaKaon[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.nsigmakaon.maximum=%f", fNSigmaKaon[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.nsigmaproton.minimum=%f", fNSigmaProton[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.nsigmaproton.maximum=%f", fNSigmaProton[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.nhits.minimum=%i", fNHits[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.nhits.maximum=%i", fNHits[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.pt.minimum=%f", fPt[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.pt.maximum=%f", fPt[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.rapidity.minimum=%f", fRapidity[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.rapidity.maximum=%f", fRapidity[1]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.dca.minimum=%f", fDCA[0]);
  tListSetttings->AddLast(new TObjString(buf));
  snprintf(buf, 200, "AliFemtoBasicTrackCut.dca.maximum=%f", fDCA[1]);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}
