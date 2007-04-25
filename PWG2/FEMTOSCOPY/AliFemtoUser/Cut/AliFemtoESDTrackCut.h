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
 * Revision 1.4  2007/03/20 09:37:13  mchojnacki
 * *** empty log message ***
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

#ifndef AliFemtoESDTrackCut_hh
#define AliFemtoESDTrackCut_hh

//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "Base/AliFemtoTrackCut.h"

class AliFemtoESDTrackCut : public AliFemtoTrackCut 
{

    public:
           AliFemtoESDTrackCut();
	   //~AliFemtoESDTrackCut();

	   virtual bool Pass(const AliFemtoTrack*);

	   virtual AliFemtoString Report();

	   void SetPt(const float& lo, const float& hi);
	   void SetRapidity(const float& lo, const float& hi);
	   void SetCharge(const int&);
	   void SetPidProbElectron(const float& lo, const float& hi);
	   void SetPidProbPion(const float& lo, const float& hi);
	   void SetPidProbKaon(const float& lo, const float& hi);
	   void SetPidProbProton(const float& lo, const float& hi);
	   void SetPidProbMuon(const float& lo, const float& hi);
	   void SetLabel(const bool& flag);
	   void SetStatus(const long& );
	   void SetminTPCclsF(const short& );
	   void SetminITScls(const int& );
  
	   private:   // here are the quantities I want to cut on...

	   int               fCharge;
	   float             fPt[2];
	   float             fRapidity[2];
	   float             fPidProbElectron[2]; // new
	   float             fPidProbPion[2]; // new
	   float             fPidProbKaon[2]; // new
	   float             fPidProbProton[2]; // new
	   float             fPidProbMuon[2]; //new 
	   bool              fLabel;//if true label<0 will not pass throught 
	   long              fStatus;//staus flag
	   short             fminTPCclsF;//min number of findable clusters in the TPC
	   int               fminITScls;//min number of clusters assigned in the ITS 
	   long              fNTracksPassed;
	   long              fNTracksFailed;

#ifdef __ROOT__ 
  ClassDef(AliFemtoESDTrackCut, 1)
#endif
};


inline void AliFemtoESDTrackCut::SetPt(const float& lo, const float& hi){fPt[0]=lo; fPt[1]=hi;}
inline void AliFemtoESDTrackCut::SetRapidity(const float& lo,const float& hi){fRapidity[0]=lo; fRapidity[1]=hi;}
inline void AliFemtoESDTrackCut::SetCharge(const int& ch){fCharge = ch;}
inline void AliFemtoESDTrackCut::SetPidProbElectron(const float& lo,const float& hi){fPidProbElectron[0]=lo; fPidProbElectron[1]=hi;}
inline void AliFemtoESDTrackCut::SetPidProbPion(const float& lo,const float& hi){fPidProbPion[0]=lo; fPidProbPion[1]=hi;}
inline void AliFemtoESDTrackCut::SetPidProbKaon(const float& lo,const float& hi){fPidProbKaon[0]=lo; fPidProbKaon[1]=hi;}
inline void AliFemtoESDTrackCut::SetPidProbProton(const float& lo,const float& hi){fPidProbProton[0]=lo; fPidProbProton[1]=hi;}
inline void AliFemtoESDTrackCut::SetPidProbMuon(const float& lo,const float& hi){fPidProbMuon[0]=lo; fPidProbMuon[1]=hi;}
inline void AliFemtoESDTrackCut::SetLabel(const bool& flag){fLabel=flag;}
inline void AliFemtoESDTrackCut::SetStatus(const long& status){fStatus=status;}
inline void AliFemtoESDTrackCut::SetminTPCclsF(const short& minTPCclsF){fminTPCclsF=minTPCclsF;}
inline void AliFemtoESDTrackCut::SetminITScls(const int& minITScls){fminITScls=minITScls;}

#endif

