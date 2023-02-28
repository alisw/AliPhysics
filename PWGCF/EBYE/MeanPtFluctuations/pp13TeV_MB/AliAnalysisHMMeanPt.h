#ifndef AliAnalysisHMMeanPt_H
#define AliAnalysisHMMeanPt_H
class TH1F;
class TH1I;
class TH2;
class TH2F;
class THn;
class TH1D;
class TH2D;
class TList;
class AliAODtrack;
class AliAnalysisUtils;
class TTree;
class AliAODEvent;
class AliVEvent;
class TString;
class TObjArray;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "THnSparse.h"
#include "THn.h"
#include "TTreeStream.h"
#include "TMath.h"
#include "TStatistic.h"
#include "TH1D.h"
#include "TString.h"

class AliAnalysisHMMeanPt : public AliAnalysisTaskSE
{
public:
    AliAnalysisHMMeanPt();
    AliAnalysisHMMeanPt(const char *name);
    virtual                 ~AliAnalysisHMMeanPt();
    
    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    AliEventCuts            fEventCuts;       //! Event cuts
 
//===========

    void                    Seteventvrtx(float vrtxz1, float vrtxz2)	{vzmin	=vrtxz1; vzmax = vrtxz2;};
    void                    SettrackBit(int trackBit)					{ftrackBit =trackBit;};
    void                    Settrackpt(float pt_low, float pt_up)		{ptmin	=pt_low; ptmax	=pt_up;};    
    void                    Settpcrows(int ftpcrows)					{TPCrows = ftpcrows;};    
    void					SetTreeName(TString nam)					{Nameoftree = nam;};

private:

    	bool AcceptTrack(AliAODTrack* aodtrack) const;    
 
     	int                   	ftrackBit; 
     	float                   vzmin;     	
      	float                   vzmax;     	
      	float                   fevtvrtxz;     	
   		float				  	ptmin;
     	float					ptmax;
		int						TPCrows;
		TString					Nameoftree;
     
        
       AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
   		TH1I					*fEventCount;	//!

		THnSparse				*fTHnetaptphi; //!
		
		TH1F					*fHistSignal;	//!
		TH1F					*fHistClustersTPC;	//!
		TH1F					*fHistChi2perNDF;		//!
		TH1F					*fHistTPCFoundFrac;	//!
		TH1F					*fHistTPCcrossrows;		//!

		TTree					*fTreept;			//!
		float					field;			//!
		float 		        	pt[200];             	//!       
		float 					cent;			//!
		int 					charge[200];		//!
		int						nevt;			//!
		int						ns;			//!
		int						nch;			//!
		float					fV0_total;		//!
        int						fRunNumber;		//!
		int						sample;		//!
		int 					Vz;				//!
		//int						TPCrows[500];	//!
		//int						filterbit[500];	//!

		TH2D*					fMulttrkV0;			//!
		TH2D*					fMulttrkV0_meanpt;	//!
		TH2D*					fMulttrk_pair;		//!
		TH2D*					fMulttrk_meanpt;		//!		
						
		TH1F*					fMultV0A;			//!
		TH1F*					fMultV0C;			//!
		TH1F*					fMultV0;			//!
		TH1F*					htrk;				//!
		
		TString 				suffixName;			//!
		
		TH1D                    *fHistMag;  //!
		TH1D                    *fHistVx; //! Vx hist
    	TH1D                    *fHistVy; //! Vy hist
    	TH1D                    *fHistVz; //! Vz hist
    	        
        TH1F                    *fHistEta;      	//!  
        TH1F                    *fHistPt;       	//! 
        TH1F                    *fHistPhi;      	//!          
        TH1F                    *fCent;   
                     
    	TH1D        	    	*hHighMultRuns;		//!
		TH2D 					*hScale;	//!
		TH2D 					*hTracks;	//!
		TH2D 					*hTrackssq;	//!
		TH2D					*hTrackavgpt;		//!
		TH2D					*hTrackpair;		//!
		TH2D					*hTrackavgptsq;		//!
		TH2D					*hTracksumpt;		//!
		TH2D					*hTracksumptiptj;		//!
  
        TH1F                    *fHistPtcal;    			//!
        
    AliAnalysisHMMeanPt(const AliAnalysisHMMeanPt&);
    AliAnalysisHMMeanPt& operator=(const AliAnalysisHMMeanPt&);
    ClassDef(AliAnalysisHMMeanPt, 1);
    
};

#endif


