/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoMCTrackCut: A basic track cut that used information from                //
// ALICE MC to accept or reject the track.                                         //
// Enables the selection on charge, transverse momentum, rapidity,                 //
// and PDG of the particle                                                         //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//  	       				                                           //
/////////////////////////////////////////////////////////////////////////////////////



#include "AliFemtoMCTrackCut.h"
#include "AliFemtoModelHiddenInfo.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoMCTrackCut);
  /// \endcond
#endif


AliFemtoMCTrackCut::AliFemtoMCTrackCut() :
    fCharge(0),
    fLabel(false),
    fPDGcode(0),
    fNTracksPassed(0),
    fNTracksFailed(0)
{
  // Default constructor
  fPt[0]=0.0;              fPt[1]=100.0;//100
  fRapidity[0]=-100;       fRapidity[1]=100;//-2 2
  fEta[0]=-2;              fEta[1]=2;//-2 2
}
//------------------------------
AliFemtoMCTrackCut::~AliFemtoMCTrackCut()
{
  /* noop */
}
//------------------------------
bool AliFemtoMCTrackCut::Pass(const AliFemtoTrack* track)
{
  if (fLabel)
    {
      if(track->Label()<0)
        {
          fNTracksFailed++;
          return false;
        }
    }
  if (fCharge!=0)
    {
      if (fCharge==10)
        {
          if(track->Charge()==0){
            fNTracksFailed++;
            return false;
          }
        }
      else if (track->Charge()!= fCharge)
        {
          fNTracksFailed++;
          return false;
        }
    }
  if (fPDGcode!=0)
    {
       if(fPDGcode==11 || fPDGcode==-11 )
         { if(!fMass) fMass=0.000511;
           if (track->PidProbElectron()!=1000)
             {
               fNTracksFailed++;
               return false;
             }
         }
       else if(fPDGcode==13 || fPDGcode==-13)
         {
           if (track->PidProbMuon()!=1000)
             {if(!fMass) fMass=0.105658;
               fNTracksFailed++;
               return false;
             }
         }
       else if(fPDGcode==211 || fPDGcode==-211 )
         {
           if (track->PidProbPion()!=1000)
             {if(!fMass) fMass= 0.1395699;
               fNTracksFailed++;
               return false;
             }
         }
       else if(fPDGcode==2212 || fPDGcode==-2212 )
         { if(!fMass) fMass=0.938272013;
           if (track->PidProbProton()!=1000)
             {
               fNTracksFailed++;
               return false;
             }
         }
       else if(fPDGcode==321 || fPDGcode==-321 )
         { if(!fMass) fMass=0.493677;
           if (track->PidProbKaon()!=1000)
             {
               fNTracksFailed++;
               return false;
             }
         }
       else if(((AliFemtoModelHiddenInfo*)track->GetHiddenInfo())->GetPDGPid()!=fPDGcode)	 // checking particle PDG from hidden info
         {
           fNTracksFailed++;
           return false;
         }
     }
  float tEnergy = ::sqrt(track->P().Mag2()+fMass*fMass);
  //cout<<"MCTrackCut: tEnergy: "<<tEnergy<<endl;
  //cout<<"MCTrackCut: track->P().z(): "<<track->P().z()<<endl;
  //cout<<"MCTrackCut: tEnergy-track->P().z(): "<<tEnergy-track->P().z()<<endl;
  float tRapidity;
  if(tEnergy-track->P().z() == 0 || (tEnergy+track->P().z())/(tEnergy-track->P().z()) == 0)
    {
      fNTracksFailed++;
      return false;
    }
  else
    {
      if((tEnergy+track->P().z())/(tEnergy-track->P().z())>0)
        tRapidity = 0.5*TMath::Log((tEnergy+track->P().z())/(tEnergy-track->P().z()));
      else
        tRapidity = 0;
    }
  float tPt = ::sqrt((track->P().x())*(track->P().x())+(track->P().y())*(track->P().y()));
  float tEta = track->P().PseudoRapidity();
  if ((tRapidity<fRapidity[0])||(tRapidity>fRapidity[1]))
    {
      fNTracksFailed++;
      return false;
    }
  if ((tEta<fEta[0])||(tEta>fEta[1]))
    {
      fNTracksFailed++;
      return false;
    }
  if ((tPt<fPt[0])||(tPt>fPt[1]))
    {
      fNTracksFailed++;
      return false;
    }

  fNTracksPassed++;
  return true;
}
//------------------------------
AliFemtoString AliFemtoMCTrackCut::Report()
{
  // Prepare report from the execution
  AliFemtoString report("AliFemtoMCTrackCut report:\n");
  report += Form("Particle mass:\t%E\n",this->Mass());
  report += Form("Particle charge:\t%d\n",fCharge);
  report += Form("Particle pT:\t%E - %E\n",fPt[0],fPt[1]);
  report += Form("Particle rapidity:\t%E - %E\n",fRapidity[0],fRapidity[1]);
  report += Form("Particle eta:\t%E - %E\n",fEta[0],fEta[1]);
  report += Form("Number of tracks which passed:\t%ld  Number which failed:\t%ld\n",fNTracksPassed,fNTracksFailed);

  return report;
}

TList *AliFemtoMCTrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char *buf = nullptr;

  buf = Form("AliFemtoMCTrackCut.mass=%f", this->Mass());
  tListSetttings->AddLast(new TObjString(buf));

  buf = Form("AliFemtoMCTrackCut.charge=%i", fCharge);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}
