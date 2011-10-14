/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoMCTrackCut: A basic track cut that used information from                //
// ALICE MC to accept or reject the track.                                         //  
// Enables the selection on charge, transverse momentum, rapidity,                 //
// and PDG of the particle							   //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//  	       				                                           //
/////////////////////////////////////////////////////////////////////////////////////



#include "AliFemtoMCTrackCut.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoMCTrackCut)
#endif


AliFemtoMCTrackCut::AliFemtoMCTrackCut() :
    fCharge(0),
    fLabel(0),
    fNTracksPassed(0),
    fNTracksFailed(0)
{
  // Default constructor
  fNTracksPassed = fNTracksFailed = 0;
  fCharge = 0;  // takes both charges 0
  fPt[0]=0.0;              fPt[1] = 100.0;//100
  fPDGcode = 0;
  fRapidity[0]=-2;       fRapidity[1]=2;//-2 2
  fEta[0]=-2;       fEta[1]=2;//-2 2
  fLabel=false;
}
//------------------------------
AliFemtoMCTrackCut::~AliFemtoMCTrackCut(){
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
       if(fPDGcode==13 || fPDGcode==-13)
       {
	 if (track->PidProbMuon()!=1000)
	   {if(!fMass) fMass=0.105658;
	     fNTracksFailed++;
	     return false;
	   }
	   }
       if(fPDGcode==211 || fPDGcode==-211 )
       {
	 if (track->PidProbPion()!=1000)
	   {if(!fMass) fMass= 0.1395699;
	     fNTracksFailed++;
	     return false;
	   }
	   }
       if(fPDGcode==2212 || fPDGcode==-2212 )
       { if(!fMass) fMass=0.938272013;
	 if (track->PidProbProton()!=1000)
	   {
	     fNTracksFailed++;
	     return false;
	   }
	   }
      if(fPDGcode==321 || fPDGcode==-321 )
      { if(!fMass) fMass=0.493677;
	 if (track->PidProbKaon()!=1000)
	   {
	     fNTracksFailed++;
	     return false;
	   }
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
    tRapidity = 0.5*::log((tEnergy+track->P().z())/(tEnergy-track->P().z()));
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

  fNTracksPassed++ ;
  return true;
    
    
}
//------------------------------
AliFemtoString AliFemtoMCTrackCut::Report()
{
  // Prepare report from the execution
  string tStemp;
  char tCtemp[100];
  sprintf(tCtemp,"Particle mass:\t%E\n",this->Mass());
  tStemp=tCtemp;
  sprintf(tCtemp,"Particle charge:\t%d\n",fCharge);
  tStemp+=tCtemp;
  sprintf(tCtemp,"Particle pT:\t%E - %E\n",fPt[0],fPt[1]);
  tStemp+=tCtemp;
  sprintf(tCtemp,"Particle rapidity:\t%E - %E\n",fRapidity[0],fRapidity[1]);
  tStemp+=tCtemp; 
  sprintf(tCtemp,"Particle eta:\t%E - %E\n",fEta[0],fEta[1]);
  tStemp+=tCtemp;
 sprintf(tCtemp,"Number of tracks which passed:\t%ld  Number which failed:\t%ld\n",fNTracksPassed,fNTracksFailed);
  tStemp += tCtemp;
  AliFemtoString returnThis = tStemp;
  return returnThis;
}
TList *AliFemtoMCTrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoMCTrackCut.mass=%f", this->Mass());
  tListSetttings->AddLast(new TObjString(buf));

  snprintf(buf, 200, "AliFemtoMCTrackCut.charge=%i", fCharge);
  tListSetttings->AddLast(new TObjString(buf));
  return tListSetttings;
}

			    
