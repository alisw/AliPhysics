/***************************************************************************
 *
 * $Id: 
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *    a simple particle cut that selects on phasespace, #hits, DCA, and PID          
 *
 ***************************************************************************
 *
 * $Log:
 **************************************************************************/

#include "Cut/AliFemtoBasicTrackCut.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoBasicTrackCut)
#endif

AliFemtoBasicTrackCut::AliFemtoBasicTrackCut(){
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
    float TEnergy = ::sqrt(track->P().mag2()+fMass*fMass);
    float TRapidity = 0.5*::log((TEnergy+track->P().z())/
			    (TEnergy-track->P().z()));

    float Pt = ::sqrt((track->P().x())*(track->P().x())+
                    (track->P().y())*(track->P().y()));

    
    /*
      cout << " * DCAxy " << (track->DCAxy()  > fDCA[0]) && (track->DCAxy()  < fDCA[1]);
      cout << " * fDCA[0] " << fDCA[0];
      cout << " * fDCA[1] " << fDCA[1];
      cout << " * track->DCAxy " << track->DCAxy();
      cout << " * NHits " <<  (track->NHits() > fNHits[0]) && (track->NHits() < fNHits[1]); 
      cout << " * Pt " << (Pt > fPt[0]) && (Pt < fPt[1]);
      cout << " * y " << (TRapidity > fRapidity[0]) && (TRapidity < fRapidity[1]);
      cout << endl;
    */

    bool goodTrack=
      (//(track->DCAxy()  > fDCA[0]) &&
     //  (track->DCAxy()  < fDCA[1]) &&
  //     (track->NHits() > fNHits[0]) &&
    //   (track->NHits() < fNHits[1]) &&
       (Pt             > fPt[0]) &&
       (Pt             < fPt[1]) &&
       (TRapidity      > fRapidity[0]) &&
       (TRapidity      < fRapidity[1]))&&
       (track->PidProbPion()>0.5)&&//moje
       (track->PidProbMuon()<0.47)&&//moje
       (track->Label()>0);//moje

    //    cout << track->DCAxy() << " " << track->NHits() << " " << Pt << " " << TRapidity << " " << TEnergy << endl;

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
  string Stemp;
  char Ctemp[100];
  sprintf(Ctemp,"Particle mass:\t%E\n",this->Mass());
  Stemp=Ctemp;
  sprintf(Ctemp,"Particle charge:\t%d\n",fCharge);
  Stemp=Ctemp;
  sprintf(Ctemp,"Particle Nsigma from pion:\t%E - %E\n",fNSigmaPion[0],fNSigmaPion[1]);
  Stemp+=Ctemp;
  sprintf(Ctemp,"Particle Nsigma from kaon:\t%E - %E\n",fNSigmaKaon[0],fNSigmaKaon[1]);
  Stemp+=Ctemp;
  sprintf(Ctemp,"Particle Nsigma from proton:\t%E - %E\n",fNSigmaProton[0],fNSigmaProton[1]);
  Stemp+=Ctemp;
  sprintf(Ctemp,"Particle #hits:\t%d - %d\n",fNHits[0],fNHits[1]);
  Stemp+=Ctemp;
  sprintf(Ctemp,"Particle pT:\t%E - %E\n",fPt[0],fPt[1]);
  Stemp+=Ctemp;
  sprintf(Ctemp,"Particle rapidity:\t%E - %E\n",fRapidity[0],fRapidity[1]);
  Stemp+=Ctemp;
  sprintf(Ctemp,"Particle DCA:\t%E - %E\n",fDCA[0],fDCA[1]);
  Stemp+=Ctemp;
  sprintf(Ctemp,"Number of tracks which passed:\t%ld  Number which failed:\t%ld\n",fNTracksPassed,fNTracksFailed);
  Stemp += Ctemp;
  AliFemtoString returnThis = Stemp;
  return returnThis;
}
