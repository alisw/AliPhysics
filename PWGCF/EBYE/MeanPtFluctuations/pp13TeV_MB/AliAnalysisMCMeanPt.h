#ifndef AliAnalysisMCMeanPt_H
#define AliAnalysisMCMeanPt_H

class TH1F;
class TH1I;
class TH2;
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
class AliAODMCParticle;
class TObjArray;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliMultSelection.h"
#include "AliMCEvent.h"
#include "THnSparse.h"
#include "THn.h"
#include "TRandom.h"
#include "TTreeStream.h"


class AliAnalysisMCMeanPt : public AliAnalysisTaskSE
{
public:
    AliAnalysisMCMeanPt();
    AliAnalysisMCMeanPt(const char *name);
    virtual                 ~AliAnalysisMCMeanPt();

    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    AliEventCuts 	        fEventCuts;                   // Event cuts

//===added for systematics===========
    void                    SetMaxTPCCluster(int MaxTPCclus){fCutTPCMaxCls	=MaxTPCclus;}
    void                    SetNTPCCluster(int TPCNclus)	{fCutTPCNCls	=TPCNclus;}
    void                    Setzvtxcut(float zvtxcut1, float zvtxcut2)			{fzvtxcutmin		=zvtxcut1; fzvtxcutmax		=zvtxcut2;};
    void                    SettrackBit(int trackBit)		{ftrackBit		=trackBit;};


private:
    	Bool_t AcceptTrack(AliAODTrack* aodtrack) const;

	    int                     fAnalysisType;

//===============
    float                   fzvtxcutmin;
     float                   fzvtxcutmax;
   int                   ftrackBit;
    int                   fCutTPCMaxCls;
    int                   fCutTPCNCls;
//===============

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TClonesArray*			AODMCTrackArray;

		TTree					*fTreeFB768rec;//!
		TTree					*fTreeFB768gen;//!
		float 		        	pt[5000];//!
		float 		        	phi[5000];//!
		float 		        	eta[5000];//!
		float 					cent;//!
		float					field;//!
		int 					charge[5000];//!
		int						TPCrows[500];//!
		int						nevt;//!
		int						ns;//!
		int						isample;//!
        int						fRunNumber;//!
		int 					Vz;//!


   		TH1I					*fEventCount;//!
        TH1F                    *fCentMC;//!
		TH1F					*fHistPtgen;//!
		TH1F					*fHistPtrec;//!

		THnSparse				*fTHnfcentetaptphi;//!
		THnSparse				*fTHnfcentetaptphigen;//!
		THnSparse 				*fTHnfnch;//!
		THnSparse 				*fTHnfnchP;//!

   		TH1D                    *fHistVx; //! Vx hist
    	TH1D                    *fHistVy; //! Vy hist
    	TH1D                    *fHistVz; //! Vz hist
    	TH1D                    *fHistVzcut; //! Vz hist

		TH2D					*hsec_wd;//!
		TH2D					*hsec_mat;//!
		TH2D					*hsec;//!
		TH2D					*hrec;//!
		TH2D           			*hprim;//!
		TH2D					*hgen;//!//!
		TH2D           			*hprimgen;//!

        TH1F                    *fHistEta_rec;//!
        TH1F                    *fHistPt_rec;//!
        TH1F                    *fHistPhi_rec;//!
		TH1F                    *fHistEta_gen;//!
        TH1F                    *fHistPt_gen;//!
        TH1F                    *fHistPhi_gen;//!

		TH2D 					*hScale_gen;//!
		TH2D 					*hTracks_gen;//!
		TH2D 					*hTracksq_gen;//!
		TH2D					*hTrackavgpt_gen;//!
		TH2D					*hTrackpair_gen;//!
		TH2D					*hTrackavgptsq_gen;//!

		TH2D 					*hScale_rec;//!
		TH2D 					*hTracks_rec;//!
		TH2D 					*hTracksq_rec;//!
		TH2D					*hTrackavgpt_rec;//!
		TH2D					*hTrackpair_rec;//!
		TH2D					*hTrackavgptsq_rec;//!

		TH2D					*hTracksumpt_rec;//!
		TH2D					*hTracksumptiptj_rec;//!
		TH2D					*hTracksumpt_gen;//!
		TH2D					*hTracksumptiptj_gen;//!


    AliAnalysisMCMeanPt(const AliAnalysisMCMeanPt&);
    AliAnalysisMCMeanPt& operator=(const AliAnalysisMCMeanPt&);
    ClassDef(AliAnalysisMCMeanPt, 1);
};


#endif
