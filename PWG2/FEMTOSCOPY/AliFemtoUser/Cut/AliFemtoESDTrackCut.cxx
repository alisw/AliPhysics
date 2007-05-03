/***************************************************************************
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

#include "AliFemtoESDTrackCut.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoESDTrackCut)
#endif

AliFemtoESDTrackCut::AliFemtoESDTrackCut() :
  fCharge(0),
  fLabel(0),
  fStatus(0),
  fminTPCclsF(0),
  fminITScls(0),
  fNTracksPassed(0),
  fNTracksFailed(0)
{
    fNTracksPassed = fNTracksFailed = 0;
    fCharge = 0;  // takes both charges 0
    fPt[0]=0.0;              fPt[1] = 100.0;//100
    fRapidity[0]=-2;       fRapidity[1]=2;//-2 2
    fPidProbElectron[0]=-1;fPidProbElectron[1]=2;
    fPidProbPion[0]=-1;    fPidProbPion[1]=2;
    fPidProbKaon[0]=-1;fPidProbKaon[1]=2;
    fPidProbProton[0]=-1;fPidProbProton[1]=2;
    fPidProbMuon[0]=-1;fPidProbMuon[1]=2;
    fLabel=false;
    fStatus=0;
    fminTPCclsF=0;
    fminITScls=0;
    
}
//------------------------------
//AliFemtoESDTrackCut::~AliFemtoESDTrackCut(){
//  /* noop */
//}
//------------------------------
bool AliFemtoESDTrackCut::Pass(const AliFemtoTrack* track)
{
    //cout<<"AliFemtoESD  cut"<<endl;
    //cout<<fPidProbPion[0]<<" < pi ="<<track->PidProbPion()<<" <"<<fPidProbPion[1]<<endl;
    if (fStatus!=0)
    {
	//cout<<" status "<<track->Label()<<" "<<track->Flags()<<" "<<track->TPCnclsF()<<" "<<track->ITSncls()<<endl;
	if ((track->Flags()&fStatus)!=fStatus)
	{
	  //  cout<<track->Flags()<<" "<<fStatus<<" no go through status"<<endl;
	    return false;
	}
	
    }
    if (fminTPCclsF>track->TPCnclsF())
    {
	//cout<<" No go because TPC Number of ClsF"<<fminTPCclsF<< " "<<track->TPCnclsF()<<endl;
	return false;
    }
    if (fminITScls>track->ITSncls())
    {
	//cout<<" No go because ITS Number of Cls"<<fminITScls<< " "<<track->ITSncls()<<endl;
	return false;
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
    float TEnergy = ::sqrt(track->P().mag2()+fMass*fMass);
    float TRapidity = 0.5*::log((TEnergy+track->P().z())/(TEnergy-track->P().z()));
    float Pt = ::sqrt((track->P().x())*(track->P().x())+(track->P().y())*(track->P().y()));
    if ((TRapidity<fRapidity[0])||(TRapidity>fRapidity[1]))
    {
	fNTracksFailed++;
	//cout<<"No Go Through the cut"<<endl;   
	//cout<<fRapidity[0]<<" < Rapidity ="<<TRapidity<<" <"<<fRapidity[1]<<endl;
	return false;
    }
    if ((Pt<fPt[0])||(Pt>fPt[1]))
    {
	fNTracksFailed++;
	//cout<<"No Go Through the cut"<<endl;
	//cout<<fPt[0]<<" < Pt ="<<Pt<<" <"<<fPt[1]<<endl;
	return false;
    }
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
  
   // cout<<"Go Through the cut"<<endl;
   // cout<<fLabel<<" Label="<<track->Label()<<endl;
   // cout<<fCharge<<" Charge="<<track->Charge()<<endl;
    // cout<<fPt[0]<<" < Pt ="<<Pt<<" <"<<fPt[1]<<endl;
    //cout<<fRapidity[0]<<" < Rapidity ="<<TRapidity<<" <"<<fRapidity[1]<<endl;
    //cout<<fPidProbElectron[0]<<" <  e="<<track->PidProbElectron()<<" <"<<fPidProbElectron[1]<<endl;
    //cout<<fPidProbPion[0]<<" <  pi="<<track->PidProbPion()<<" <"<<fPidProbPion[1]<<endl;
    //cout<<fPidProbKaon[0]<<" <  k="<<track->PidProbKaon()<<" <"<<fPidProbKaon[1]<<endl;
    //cout<<fPidProbProton[0]<<" <  p="<<track->PidProbProton()<<" <"<<fPidProbProton[1]<<endl;
    //cout<<fPidProbMuon[0]<<" <  mi="<<track->PidProbMuon()<<" <"<<fPidProbMuon[1]<<endl;
    fNTracksPassed++ ;
    return true;
    
    
}
//------------------------------
AliFemtoString AliFemtoESDTrackCut::Report()
{
    string Stemp;
    char Ctemp[100];
    sprintf(Ctemp,"Particle mass:\t%E\n",this->Mass());
    Stemp=Ctemp;
    sprintf(Ctemp,"Particle charge:\t%d\n",fCharge);
    Stemp+=Ctemp;
    sprintf(Ctemp,"Particle pT:\t%E - %E\n",fPt[0],fPt[1]);
    Stemp+=Ctemp;
    sprintf(Ctemp,"Particle rapidity:\t%E - %E\n",fRapidity[0],fRapidity[1]);
    Stemp+=Ctemp;
    sprintf(Ctemp,"Number of tracks which passed:\t%ld  Number which failed:\t%ld\n",fNTracksPassed,fNTracksFailed);
    Stemp += Ctemp;
    AliFemtoString returnThis = Stemp;
    return returnThis;
}
