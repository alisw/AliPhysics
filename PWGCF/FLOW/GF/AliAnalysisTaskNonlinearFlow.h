#ifndef ALIANALYSISTASKSECHARGEDFLOW_H
#define ALIANALYSISTASKSECHARGEDFLOW_H
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include <TComplex.h>

#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TComplex.h>
#include <TBits.h>
#include <TRandom3.h>
// AliRoot includes
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

#ifndef __CINT__
// ROOT includes
// # include <TList.h>
// # include <TH1.h>
// # include <TH2.h>
// # include <TH3.h>
// # include <TProfile.h>
// # include <TComplex.h>
// # include <TBits.h>
# include <TRandom3.h>
// AliRoot includes
// # include "AliESDEvent.h"
// # include "AliAODEvent.h"
// # include "AliVEvent.h"
// # include "AliVTrack.h"
// # include "AliVVertex.h"
// # include "AliAnalysisFilter.h"
// # include "AliESDtrackCuts.h"
#else
#endif
class TList;
class TF1;
class TH1;
class TH2;
class TH3F;
class TH1F;
class TH2F;
class TH3D;
class TProfile;
class TProfile2D;
class TComplex;
class AliESDEvent;
class AliAODEvent;
class AliVEvent;
class AliVTrack;
class AliVVertex;
class AliESDtrackCuts;
class AliAODITSsaTrackCuts;
class AliInputEventHandler;

#include <THnSparse.h>

class PhysicsProfile {
  public:
    // Physics profiles
    TProfile*		fChMCsc4242;								//! same profiles as for real data (or reconstructed part of MC)
    TProfile*		fChMCsc4242Gap0;						//!
    TProfile*		fChMCsc4242Gap2;						//!
    TProfile*		fChMCsc4242_3sub;						//!
    TProfile*		fChMCsc4224_3sub;						//!
    TProfile*		fChMCsc4242_3subGap2;				//!
    TProfile*		fChMCsc4224_3subGap2;				//!
    TProfile*		fChMCsc3232;								//!
    TProfile*		fChMCsc3232Gap0;						//!
    TProfile*		fChMCsc3232Gap2;						//!
    TProfile*		fChMCsc3232_3sub;						//!
    TProfile*		fChMCsc3223_3sub;						//!
    TProfile*		fChMCsc3232_3subGap2;				//!
    TProfile*		fChMCsc3223_3subGap2;				//!

    TProfile*		fMCsc4242[12];							//!
    TProfile*		fMCsc4242Gap0[12];					//!
    TProfile*		fMCsc4242Gap2[12];					//!
    TProfile*		fMCsc4242_3sub[12];					//!
    TProfile*		fMCsc4224_3sub[12];					//!
    TProfile*		fMCsc4242_3subGap2[12];			//!
    TProfile*		fMCsc4224_3subGap2[12];			//!
    TProfile*		fMCsc3232[12];							//!
    TProfile*		fMCsc3232Gap0[12];					//!
    TProfile*		fMCsc3232Gap2[12];					//!
    TProfile*		fMCsc3232_3sub[12];					//!
    TProfile*		fMCsc3223_3sub[12];					//!
    TProfile*		fMCsc3232_3subGap2[12];			//!
    TProfile*		fMCsc3223_3subGap2[12];			//!

    TProfile*		fChMCcn2Ntrks1bin[6];  			//!
    TProfile*		fChMCcn2Gap0Ntrks1bin[6];		//!
    TProfile*		fChMCcn2Gap2Ntrks1bin[6];		//!
    TProfile*		fChMCcn2Gap4Ntrks1bin[6];		//!
    TProfile*   fChMCcn2Gap6Ntrks1bin[6];        //!
    TProfile*		fChMCcn2Gap8Ntrks1bin[6];		//!
    TProfile*		fChMCcn2Gap10Ntrks1bin[6];  	//!
    TProfile*		fChMCcn2Gap14Ntrks1bin[6];	//!
    TProfile*		fChMCcn2Gap16Ntrks1bin[6];	//!
    TProfile*		fChMCcn2Gap18Ntrks1bin[6];	//!

    TProfile*		fChMCcn2_3subLMNtrks1bin[6];			//!
    TProfile*		fChMCcn2_3subRMNtrks1bin[6];			//!
    TProfile*		fChMCcn2_3subGap2LMNtrks1bin[6];	//!
    TProfile*		fChMCcn2_3subGap2RMNtrks1bin[6];	//!

    TProfile*		fChMCcn4Ntrks1bin[6];				//!
    TProfile*		fChMCcn4Gap0Ntrks1bin[6];		//!
    TProfile*		fChMCcn4Gap2Ntrks1bin[6];		//!
    TProfile*		fChMCcn4_3subNtrks1bin[6];	//!
    TProfile*		fChMCcn4_3subGap2Ntrks1bin[6];//!

    TProfile*		fChMCcn6Ntrks1bin[6];				//!
    TProfile*		fChMCcn6Gap0Ntrks1bin[6];		//!

    TProfile*		fChMCcn8Ntrks1bin[6];				//!
    TProfile*		fChMCcn8Gap0Ntrks1bin[6];		//!


    TProfile*		fMCcn2Ntrks1bin[6][12];			//!
    TProfile*		fMCcn2Gap0Ntrks1bin[6][12];	//!
    TProfile*		fMCcn2Gap2Ntrks1bin[6][12];	//!
    TProfile*		fMCcn2Gap4Ntrks1bin[6][12];	//!
    TProfile*   fMCcn2Gap6Ntrks1bin[6][12];    //!
    TProfile*		fMCcn2Gap8Ntrks1bin[6][12];	//!
    TProfile*		fMCcn2Gap10Ntrks1bin[6][12];	//!
    TProfile*		fMCcn2Gap14Ntrks1bin[6][12];//!
    TProfile*		fMCcn2Gap16Ntrks1bin[6][12];//!
    TProfile*		fMCcn2Gap18Ntrks1bin[6][12];//!

    TProfile*		fMCcn2_3subLMNtrks1bin[6][12];			//!
    TProfile*		fMCcn2_3subRMNtrks1bin[6][12];			//!
    TProfile*		fMCcn2_3subGap2LMNtrks1bin[6][12];	//!
    TProfile*		fMCcn2_3subGap2RMNtrks1bin[6][12];	//!

    TProfile*		fMCcn4Ntrks1bin[6][12];			//!
    TProfile*		fMCcn4Gap0Ntrks1bin[6][12];	//!
    TProfile*		fMCcn4Gap2Ntrks1bin[6][12];	//!
    TProfile*		fMCcn4_3subNtrks1bin[6][12];//!
    TProfile*		fMCcn4_3subGap2Ntrks1bin[6][12];//!

    TProfile*		fMCcn6Ntrks1bin[6][12];			//!
    TProfile*		fMCcn6Gap0Ntrks1bin[6][12];	//!

    TProfile*		fMCcn8Ntrks1bin[6][12];			//!
    TProfile*		fMCcn8Gap0Ntrks1bin[6][12];	//!

    // Physics profiles
    TProfile*		fChsc4242;									//! SC(4,2)
    TProfile*		fChsc4242Gap0;							//! SC(4,2) |#Delta#eta| > 0.0
    TProfile*		fChsc4242Gap2;							//! SC(4,2) |#Delta#eta| > 0.2
    TProfile*        fChsc4242Gap4;                            //! SC(4,2) |#Delta#eta| > 0.2
    TProfile*        fChsc4242Gap6;                            //! SC(4,2) |#Delta#eta| > 0.2
    TProfile*        fChsc4242Gap8;                            //! SC(4,2) |#Delta#eta| > 0.2
    TProfile*        fChsc4242Gap10;                            //! SC(4,2) |#Delta#eta| > 0.2
    TProfile*		fChsc4242_3sub;							//! SC(4,2)_A 3subevent method
    TProfile*		fChsc4224_3sub;							//! SC(4,2)_B 3subevent method
    TProfile*		fChsc4242_3subGap2;					//! SC(4,2)_A 3subevent method |#Delta#eta| > 0.2
    TProfile*		fChsc4224_3subGap2;					//! SC(4,2)_B 3subevent method |#Delta#eta| > 0.2
    TProfile*		fChsc3232;									//! SC(3,2)
    TProfile*		fChsc3232Gap0;							//! SC(3,2) |#Delta#eta| > 0.0
    TProfile*		fChsc3232Gap2;							//! SC(3,2) |#Delta#eta| > 0.2
    TProfile*        fChsc3232Gap4;                            //! SC(3,2) |#Delta#eta| > 0.2
    TProfile*        fChsc3232Gap6;                            //! SC(3,2) |#Delta#eta| > 0.2
    TProfile*        fChsc3232Gap8;                            //! SC(3,2) |#Delta#eta| > 0.2
    TProfile*        fChsc3232Gap10;                            //! SC(3,2) |#Delta#eta| > 0.2
    TProfile*		fChsc3232_3sub;							//! SC(3,2)_A 3subevent method
    TProfile*		fChsc3223_3sub;							//! SC(3,2)_B 3subevent method
    TProfile*		fChsc3232_3subGap2;					//! SC(3,2)_A 3subevent method |#Delta#eta| > 0.2
    TProfile*		fChsc3223_3subGap2;					//! SC(3,2)_B 3subevent method |#Delta#eta| > 0.2

    TProfile*		fsc4242[12];								//! the same but for different fBin (sampling)
    TProfile*		fsc4242Gap0[12];						//!
    TProfile*		fsc4242Gap2[12];						//!
    TProfile*        fsc4242Gap4[12];                        //!
    TProfile*        fsc4242Gap6[12];                        //!
    TProfile*        fsc4242Gap8[12];                        //!
    TProfile*        fsc4242Gap10[12];                        //!
    TProfile*		fsc4242_3sub[12];						//!
    TProfile*		fsc4224_3sub[12];						//!
    TProfile*		fsc4242_3subGap2[12];				//!
    TProfile*		fsc4224_3subGap2[12];				//!
    TProfile*		fsc3232[12];								//!
    TProfile*		fsc3232Gap0[12];						//!
    TProfile*		fsc3232Gap2[12];						//!
    TProfile*        fsc3232Gap4[12];                        //!
    TProfile*        fsc3232Gap6[12];                        //!
    TProfile*        fsc3232Gap8[12];                        //!
    TProfile*        fsc3232Gap10[12];                        //!
    TProfile*		fsc3232_3sub[12];						//!
    TProfile*		fsc3223_3sub[12];						//!
    TProfile*		fsc3232_3subGap2[12];				//!
    TProfile*		fsc3223_3subGap2[12];				//!

    // Standard correlation profiles for different harmonics
    TProfile*		fChc422Ntrks1bin; //!
    TProfile*       fChc532Ntrks1bin;//!
    TProfile*        fChc422Gap0ANtrks1bin;    //!
    TProfile*        fChc422Gap0BNtrks1bin;    //!
    TProfile*        fChc532Gap0ANtrks1bin;    //!
    TProfile*        fChc532Gap0BNtrks1bin;    //!
    TProfile*        fChc422Gap2ANtrks1bin;    //!
    TProfile*        fChc422Gap2BNtrks1bin;    //!
    TProfile*        fChc532Gap2ANtrks1bin;    //!
    TProfile*        fChc532Gap2BNtrks1bin;    //!
    TProfile*        fChc422Gap4ANtrks1bin;    //!
    TProfile*        fChc422Gap4BNtrks1bin;    //!
    TProfile*        fChc532Gap4ANtrks1bin;    //!
    TProfile*        fChc532Gap4BNtrks1bin;    //!
    TProfile*        fChc422Gap6ANtrks1bin;    //!
    TProfile*        fChc422Gap6BNtrks1bin;    //!
    TProfile*        fChc532Gap6ANtrks1bin;    //!
    TProfile*        fChc532Gap6BNtrks1bin;    //!
    TProfile*        fChc422Gap8ANtrks1bin;    //!
    TProfile*        fChc422Gap8BNtrks1bin;    //!
    TProfile*        fChc532Gap8ANtrks1bin;    //!
    TProfile*        fChc532Gap8BNtrks1bin;    //!
    TProfile*        fChc422Gap10ANtrks1bin;    //!
    TProfile*        fChc422Gap10BNtrks1bin;    //!
    TProfile*        fChc532Gap10ANtrks1bin;    //!
    TProfile*        fChc532Gap10BNtrks1bin;    //!


    TProfile*		fChcn2Ntrks1bin[6];  				//! <<2>> in unit bins of Ntrks
    TProfile*		fChcn2Gap0Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.0
    TProfile*		fChcn2Gap2Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.2
    TProfile*		fChcn2Gap4Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.4
    TProfile*     fChcn2Gap6Ntrks1bin[6];          //! <<2>> |#Delta#eta| > 0.6
    TProfile*		fChcn2Gap8Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.8
    TProfile*		fChcn2Gap10Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 1.0
    TProfile*		fChcn2Gap14Ntrks1bin[6];  	//! <<2>> |#Delta#eta| > 1.4
    TProfile*		fChcn2Gap16Ntrks1bin[6];  	//! <<2>> |#Delta#eta| > 1.6
    TProfile*		fChcn2Gap18Ntrks1bin[6];  	//! <<2>> |#Delta#eta| > 1.8

    TProfile*		fChcn2_3subLMNtrks1bin[6];  		//! <<2>> left vs. middle subevent
    TProfile*		fChcn2_3subRMNtrks1bin[6];  		//! <<2>> middle vs. right subevent
    TProfile*		fChcn2_3subGap2LMNtrks1bin[6];  //! <<2>> left vs. middle subevent |#Delta#eta| > 0.2
    TProfile*		fChcn2_3subGap2RMNtrks1bin[6];  //! <<2>> middle vs. right subevent |#Delta#eta| > 0.2

    TProfile*		fChcn4Ntrks1bin[6];  				//! <<4>> in unit bins of Ntrks
    TProfile*		fChcn4Gap0Ntrks1bin[6];  		//! <<4>> |#Delta#eta| > 0.0
    TProfile*		fChcn4Gap2Ntrks1bin[6];  		//! <<4>> |#Delta#eta| > 0.2
    TProfile*      fChcn4Gap4Ntrks1bin[6];          //! <<4>> |#Delta#eta| > 0.2
    TProfile*      fChcn4Gap6Ntrks1bin[6];          //! <<4>> |#Delta#eta| > 0.2
    TProfile*      fChcn4Gap8Ntrks1bin[6];          //! <<4>> |#Delta#eta| > 0.2
    TProfile*      fChcn4Gap10Ntrks1bin[6];          //! <<4>> |#Delta#eta| > 0.2
    TProfile*     fChcn4_3subNtrks1bin[6];         //! <<4>> 3subevent method
    TProfile*     fChcn4_3subGap2Ntrks1bin[6];//! <<4>> 3subevent method |#Delta#eta| > 0.2
    TProfile*        fChcn5Ntrks1bin[6];                  //! <<6>> in unit bins of Ntrks
    TProfile*        fChcn5Gap0ANtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.0
    TProfile*        fChcn5Gap2ANtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.2
    TProfile*        fChcn5Gap4ANtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.4
    TProfile*        fChcn5Gap6ANtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.6
    TProfile*        fChcn5Gap8ANtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.8
    TProfile*        fChcn5Gap10ANtrks1bin[6];          //! <<6>> |#Delta#eta| > 1.0
    TProfile*        fChcn5Gap0BNtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.0
    TProfile*        fChcn5Gap2BNtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.2
    TProfile*        fChcn5Gap4BNtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.4
    TProfile*        fChcn5Gap6BNtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.6
    TProfile*        fChcn5Gap8BNtrks1bin[6];          //! <<6>> |#Delta#eta| > 0.8
    TProfile*        fChcn5Gap10BNtrks1bin[6];          //! <<6>> |#Delta#eta| > 1.0
    TProfile*		 fChcn6Ntrks1bin[6];  				//! <<6>> in unit bins of Ntrks
    TProfile*		 fChcn6Gap0Ntrks1bin[6];  		//! <<6>> |#Delta#eta| > 0.0
    TProfile*        fChcn6Gap2Ntrks1bin[6];          //! <<6>> |#Delta#eta| > 0.0
    TProfile*        fChcn6Gap4Ntrks1bin[6];          //! <<6>> |#Delta#eta| > 0.0
    TProfile*        fChcn6Gap6Ntrks1bin[6];          //! <<6>> |#Delta#eta| > 0.0
    TProfile*        fChcn6Gap8Ntrks1bin[6];          //! <<6>> |#Delta#eta| > 0.0
    TProfile*        fChcn6Gap10Ntrks1bin[6];          //! <<6>> |#Delta#eta| > 0.0
    TProfile*		fChcn8Ntrks1bin[6];  				//! <<8>> in unit bins of Ntrks
    TProfile*		fChcn8Gap0Ntrks1bin[6];  		//! <<8>> |#Delta#eta| > 0.0

    // the same profiles but for different fBin (sampling)
    TProfile*		fc422Ntrks1bin[12];					//!
    TProfile*       fc532Ntrks1bin[12];                    //!
    TProfile*        fc422Gap0ANtrks1bin[12];//!
    TProfile*        fc422Gap0BNtrks1bin[12];//!
    TProfile*        fc532Gap0ANtrks1bin[12];//!
    TProfile*        fc532Gap0BNtrks1bin[12];//!
    TProfile*		fc422Gap2ANtrks1bin[12];//!
    TProfile*		fc422Gap2BNtrks1bin[12];//!
    TProfile*		fc532Gap2ANtrks1bin[12];//!
    TProfile*		fc532Gap2BNtrks1bin[12];//!
    TProfile*        fc422Gap4ANtrks1bin[12];//!
    TProfile*        fc422Gap4BNtrks1bin[12];//!
    TProfile*        fc532Gap4ANtrks1bin[12];//!
    TProfile*        fc532Gap4BNtrks1bin[12];//!
    TProfile*        fc422Gap6ANtrks1bin[12];//!
    TProfile*        fc422Gap6BNtrks1bin[12];//!
    TProfile*        fc532Gap6ANtrks1bin[12];//!
    TProfile*        fc532Gap6BNtrks1bin[12];//!
    TProfile*        fc422Gap8ANtrks1bin[12];//!
    TProfile*        fc422Gap8BNtrks1bin[12];//!
    TProfile*        fc532Gap8ANtrks1bin[12];//!
    TProfile*        fc532Gap8BNtrks1bin[12];//!
    TProfile*        fc422Gap10ANtrks1bin[12];//!
    TProfile*        fc422Gap10BNtrks1bin[12];//!
    TProfile*        fc532Gap10ANtrks1bin[12];//!
    TProfile*        fc532Gap10BNtrks1bin[12];//!


    TProfile*		fcn2Ntrks1bin[6][12];  			//!
    TProfile*		fcn2Gap0Ntrks1bin[6][12];		//!
    TProfile*		fcn2Gap2Ntrks1bin[6][12];		//!
    TProfile*		fcn2Gap4Ntrks1bin[6][12];		//!
    TProfile*     fcn2Gap6Ntrks1bin[6][12];        //!
    TProfile*		fcn2Gap8Ntrks1bin[6][12];		//!
    TProfile*		fcn2Gap10Ntrks1bin[6][12];		//!
    TProfile*		fcn2Gap14Ntrks1bin[6][12];	//!
    TProfile*		fcn2Gap16Ntrks1bin[6][12];	//!
    TProfile*		fcn2Gap18Ntrks1bin[6][12];	//!

    TProfile*		fcn2_3subLMNtrks1bin[6][12];  		//!
    TProfile*		fcn2_3subRMNtrks1bin[6][12];  		//!
    TProfile*		fcn2_3subGap2LMNtrks1bin[6][12];  //!
    TProfile*		fcn2_3subGap2RMNtrks1bin[6][12];  //!

    TProfile*		fcn4Ntrks1bin[6][12];				//!
    TProfile*		fcn4Gap0Ntrks1bin[6][12];		//!
    TProfile*		fcn4Gap2Ntrks1bin[6][12];		//!
    TProfile*        fcn4Gap4Ntrks1bin[6][12];        //!
    TProfile*        fcn4Gap6Ntrks1bin[6][12];        //!
    TProfile*        fcn4Gap8Ntrks1bin[6][12];        //!
    TProfile*        fcn4Gap10Ntrks1bin[6][12];        //!
    TProfile*		fcn4_3subNtrks1bin[6][12];	//!
    TProfile*		fcn4_3subGap2Ntrks1bin[6][12];//!
    TProfile*        fcn5Ntrks1bin[6][12];                //!
    TProfile*        fcn5Gap0ANtrks1bin[6][12];        //!
    TProfile*        fcn5Gap2ANtrks1bin[6][12];        //!
    TProfile*        fcn5Gap4ANtrks1bin[6][12];        //!
    TProfile*        fcn5Gap6ANtrks1bin[6][12];        //!
    TProfile*        fcn5Gap8ANtrks1bin[6][12];        //!
    TProfile*        fcn5Gap10ANtrks1bin[6][12];        //!
    TProfile*        fcn5Gap0BNtrks1bin[6][12];        //!
    TProfile*        fcn5Gap2BNtrks1bin[6][12];        //!
    TProfile*        fcn5Gap4BNtrks1bin[6][12];        //!
    TProfile*        fcn5Gap6BNtrks1bin[6][12];        //!
    TProfile*        fcn5Gap8BNtrks1bin[6][12];        //!
    TProfile*        fcn5Gap10BNtrks1bin[6][12];        //!
    TProfile*		fcn6Ntrks1bin[6][12];				//!
    TProfile*		fcn6Gap0Ntrks1bin[6][12];		//!
    TProfile*        fcn6Gap2Ntrks1bin[6][12];        //!
    TProfile*        fcn6Gap4Ntrks1bin[6][12];        //!
    TProfile*        fcn6Gap6Ntrks1bin[6][12];        //!
    TProfile*        fcn6Gap8Ntrks1bin[6][12];        //!
    TProfile*        fcn6Gap10Ntrks1bin[6][12];        //!

    TProfile*		fcn8Ntrks1bin[6][12];				//!
    TProfile*		fcn8Gap0Ntrks1bin[6][12];		//!

};
class AliAnalysisTaskNonlinearFlow : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskNonlinearFlow();
    AliAnalysisTaskNonlinearFlow(const char *name);

    virtual ~AliAnalysisTaskNonlinearFlow();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t* option);
    virtual void   Terminate(Option_t* );

    virtual void	SetFilterbit(Int_t bit){fFilterbit = bit;}
    virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
    virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
    virtual void  SetMinPt(Double_t minPt){fMinPt = minPt;}
    virtual void  SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
    virtual void	SetTPCclusters(Int_t tpcClus){fTPCclusters = tpcClus;}
    virtual void	SetMinITSClusters(Int_t minClus){fMinITSClus = minClus;}
    virtual void	SetMaxChi2(Double_t maxChi){fMaxChi2 = maxChi;}
    virtual void	SetUseDCAzCut(Bool_t usedcaz){fUseDCAzCut = usedcaz;}
    virtual void	SetDCAzCut(Double_t dcaz){fDCAz = dcaz;}
    virtual void	SetUseDCAxyCut(Bool_t usedcaxy){fUseDCAxyCut = usedcaxy;}
    virtual void	SetDCAxyCut(Double_t dcaxy){fDCAxy = dcaxy;}
    virtual void  SetIsSample(Int_t IsSample){fSample = IsSample;}
    virtual void  SetCentFlag(Short_t nCent){fCentFlag = nCent;}
    virtual void	SetTrigger(Int_t trig){fTrigger = trig;}
    virtual void	SetLSFlag(Bool_t LS){fLS = LS;}
    virtual void	SetNUEFlag(Bool_t NUE){fNUE = NUE;}
    virtual void	SetNUA(Bool_t NUA){fNUA = NUA;}
    //..MC
    virtual void	SetPrimariesFlag(Bool_t prim){fPrimaries = prim;}
    virtual void	SetSecondariesFlag(Bool_t sec){fSecondaries = sec;}
    virtual void	SetPionsFlag(Bool_t pion){fPions = pion;}
    virtual void	SetPiKPFlag(Bool_t pikp){fPiKP = pikp;}
    virtual void	SetMCFlag(Bool_t mc){fIsMC = mc;}
    //....
    virtual void	SetPeriod(TString period){fPeriod = period;}

  private:
    AliAnalysisTaskNonlinearFlow(const AliAnalysisTaskNonlinearFlow&);
    AliAnalysisTaskNonlinearFlow& operator=(const AliAnalysisTaskNonlinearFlow&);

    virtual void		AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus);
    Short_t					GetCentrCode(AliVEvent* ev);
    // virtual double	ProcessMCTruth(AliVEvent* aod, float VtxZ, double Nreco, bool fPlus);
    bool 						CheckPrimary(AliVEvent *aod, double label);
    bool						IsGoodPSEvent(AliVEvent *aod);
    bool						IsSPDClusterVsTrackletBG(const AliVEvent* event, bool fillHist);
    bool						IsV0C012vsTklBG         (const AliVEvent* event, bool fillHist);
    bool						IsV0Casym               (const AliVEvent* event, bool fillHist);
    bool						IsV0MOnVsOfPileup       (const AliVEvent* event, bool fillHist);
    bool						IsSPDOnVsOfPileup       (const AliVEvent* event, bool fillHist);
    bool						IsV0PFPileup            (const AliVEvent* event);

    int 						GetRunPart(int run);
    double 					GetWeight(double phi, double eta, double pt, int run, bool fPlus, double vz, double runNumber);
    double 					GetPtWeight(double pt, double eta, float vz, double runNumber);

    TComplex Qvector[20][20];
    TComplex Qvector0M[20][20];
    TComplex Qvector0P[20][20];
    TComplex Qvector2M[20][20];
    TComplex Qvector2P[20][20];
    TComplex Qvector4M[20][20];
    TComplex Qvector4P[20][20];
    TComplex Qvector6M[20][20];
    TComplex Qvector6P[20][20];
    TComplex Qvector8M[20][20];
    TComplex Qvector8P[20][20];
    TComplex Qvector10M[20][20];
    TComplex Qvector10P[20][20];
    TComplex Qvector14M[20][20];
    TComplex Qvector14P[20][20];
    TComplex Qvector16M[20][20];
    TComplex Qvector16P[20][20];
    TComplex Qvector18M[20][20];
    TComplex Qvector18P[20][20];
    TComplex QvectorSubLeft[20][20];
    TComplex QvectorSubMiddle[20][20];
    TComplex QvectorSubRight[20][20];
    TComplex QvectorSubGap2Left[20][20];
    TComplex QvectorSubGap2Middle[20][20];
    TComplex QvectorSubGap2Right[20][20];
    TComplex pvector[20][20];
    TComplex pvectorM[20][20];
    TComplex pvectorP[20][20];
    TComplex qvector[20][20];
    TComplex pvector0M[20][20];
    TComplex pvector0P[20][20];
    TComplex pvector4M[20][20];
    TComplex pvector4P[20][20];
    TComplex pvector8M[20][20];
    TComplex pvector8P[20][20];

    TComplex Q(int n, int p);
    TComplex QGap0M(int n, int p);
    TComplex QGap0P(int n, int p);
    TComplex QGap2M(int n, int p);
    TComplex QGap2P(int n, int p);
    TComplex QGap4M(int n, int p);
    TComplex QGap4P(int n, int p);
    TComplex QGap6M(int n, int p);
    TComplex QGap6P(int n, int p);
    TComplex QGap8M(int n, int p);
    TComplex QGap8P(int n, int p);
    TComplex QGap10M(int n, int p);
    TComplex QGap10P(int n, int p);
    TComplex QGap14M(int n, int p);
    TComplex QGap14P(int n, int p);
    TComplex QGap16M(int n, int p);
    TComplex QGap16P(int n, int p);
    TComplex QGap18M(int n, int p);
    TComplex QGap18P(int n, int p);
    TComplex QsubLeft(int n, int p);
    TComplex QsubMiddle(int n, int p);
    TComplex QsubRight(int n, int p);
    TComplex QsubGap2Left(int n, int p);
    TComplex QsubGap2Middle(int n, int p);
    TComplex QsubGap2Right(int n, int p);
    TComplex p(int n, int p);
    TComplex pGap10M(int n, int p);
    TComplex pGap10P(int n, int p);
    TComplex q(int n, int p);
    TComplex pGap0M(int n, int p);
    TComplex pGap0P(int n, int p);
    TComplex pGap4M(int n, int p);
    TComplex pGap4P(int n, int p);
    TComplex pGap6M(int n, int p);
    TComplex pGap6P(int n, int p);
    TComplex pGap8M(int n, int p);
    TComplex pGap8P(int n, int p);
    void ResetQ(const int nMaxHarm, const int nMaxPow);

    TComplex Two(int n1, int n2);
    TComplex TwoGap0(int n1, int n2);
    TComplex TwoGap2(int n1, int n2);
    TComplex TwoGap4(int n1, int n2);
    TComplex TwoGap6(int n1, int n2);
    TComplex TwoGap8(int n1, int n2);
    TComplex TwoGap10(int n1, int n2);
    TComplex TwoGap14(int n1, int n2);
    TComplex TwoGap16(int n1, int n2);
    TComplex TwoGap18(int n1, int n2);
    TComplex Two_3SubLM(int n1, int n2);
    TComplex Two_3SubRM(int n1, int n2);
    TComplex Two_3SubGap2LM(int n1, int n2);
    TComplex Two_3SubGap2RM(int n1, int n2);
    TComplex TwoGap0M(int n1, int n2);
    TComplex TwoGap2M(int n1, int n2);
    TComplex TwoGap4M(int n1, int n2);
    TComplex TwoGap6M(int n1, int n2);
    TComplex TwoGap8M(int n1, int n2);
    TComplex TwoGap10M(int n1, int n2);
    TComplex TwoGap0P(int n1, int n2);
    TComplex TwoGap2P(int n1, int n2);
    TComplex TwoGap4P(int n1, int n2);
    TComplex TwoGap6P(int n1, int n2);
    TComplex TwoGap8P(int n1, int n2);
    TComplex TwoGap10P(int n1, int n2);

    TComplex Three(int n1, int n2, int n3);
    TComplex ThreeGap0A(int n1, int n2, int n3);
    TComplex ThreeGap0B(int n1, int n2, int n3);
    TComplex ThreeGap2A(int n1, int n2, int n3);
    TComplex ThreeGap2B(int n1, int n2, int n3);
    TComplex ThreeGap4A(int n1, int n2, int n3);
    TComplex ThreeGap4B(int n1, int n2, int n3);
    TComplex ThreeGap6A(int n1, int n2, int n3);
    TComplex ThreeGap6B(int n1, int n2, int n3);
    TComplex ThreeGap8A(int n1, int n2, int n3);
    TComplex ThreeGap8B(int n1, int n2, int n3);
    TComplex ThreeGap10A(int n1, int n2, int n3);
    TComplex ThreeGap10B(int n1, int n2, int n3);

    TComplex ThreeGap0_subM(int n1, int n2, int n3);
    TComplex ThreeGap0_subP(int n1, int n2, int n3);
    TComplex ThreeGap2_subM(int n1, int n2, int n3);
    TComplex ThreeGap2_subP(int n1, int n2, int n3);
    TComplex ThreeGap4_subM(int n1, int n2, int n3);
    TComplex ThreeGap4_subP(int n1, int n2, int n3);
    TComplex ThreeGap6_subM(int n1, int n2, int n3);
    TComplex ThreeGap6_subP(int n1, int n2, int n3);
    TComplex ThreeGap8_subM(int n1, int n2, int n3);
    TComplex ThreeGap8_subP(int n1, int n2, int n3);
    TComplex ThreeGap10_subM(int n1, int n2, int n3);
    TComplex ThreeGap10_subP(int n1, int n2, int n3);
    TComplex Four(int n1, int n2, int n3, int n4);
    TComplex FourGap0M(int n1, int n2, int n3, int n4);
    TComplex FourGap0P(int n1, int n2, int n3, int n4);
    TComplex FourGap0(int n1, int n2, int n3, int n4);
    TComplex FourGap2(int n1, int n2, int n3, int n4);
    TComplex FourGap4(int n1, int n2, int n3, int n4);
    TComplex FourGap6(int n1, int n2, int n3, int n4);
    TComplex FourGap8(int n1, int n2, int n3, int n4);
    TComplex FourGap10(int n1, int n2, int n3, int n4);
    TComplex Four_3SubEvts(int n1, int n2, int n3, int n4);
    TComplex Four_3SubGap2Evts(int n1, int n2, int n3, int n4);
    TComplex Five(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap0A_2(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap0A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap0B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap2A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap2B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap4A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap4B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap6A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap6B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap8A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap8B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap10A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap10B(int n1, int n2, int n3, int n4, int n5);
    TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap0(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap2(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap4(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap6(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap8(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap10(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
    TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
    TComplex EightGap0(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);

    AliEventCuts	fEventCuts;					// Event cuts
    AliAODEvent* fAOD;                //! AOD object
    AliAODITSsaTrackCuts* fitssatrackcuts; //itssatrackcuts object

    // Cuts and options
    Int_t				fFilterbit;					// flag for filter bit
    Double_t		fEtaCut;						// Eta cut used to select particles
    Double_t		fVtxCut;						// Vtx cut on z position in cm
    Double_t		fMinPt;							// Min pt - for histogram limits
    Double_t		fMaxPt;							// Max pt - for histogram limits
    Int_t				fTPCclusters;				// min. TPC clusters
    Int_t				fMinITSClus;				// min ITS clusters, LHC15ijl
    Double_t		fMaxChi2;						// max chi2 per ITS cluster, LHC15ijl
    Bool_t			fUseDCAzCut;				// switch to choose whether I want to use DCAz cut or not (for systematic studies, otherwise it is in FB selection by default)
    Double_t		fDCAz;							// max DCAz, for systematics
    Bool_t			fUseDCAxyCut;				// the same switch as for DCAxy
    Double_t		fDCAxy;							// max DCAxy, for systematics
    Int_t				fSample;						// number of sample
    Short_t			fCentFlag;					// centrality flag
    Int_t				fTrigger;						// flag for trigger
    Bool_t			fLS;								// charge, 1:all, 2:pp,  3: mm
    Bool_t			fNUE;								// flag for NUE correction
    Bool_t			fNUA;								// 0: no NUA correction, 1: NUA correction
    //..MC
    Bool_t			fPrimaries;					// flag for selecting only primary particles
    Bool_t			fSecondaries;				// flag for selecting only secondary particles
    Bool_t			fPions;							// flag to select only pions using PDG information
    Bool_t			fPiKP;							// flag to select only particles with pdg<2211 (pi,K,p,e,mu)
    Bool_t			fIsMC;							// flag to do MonteCarlo part
    //....
    TString			fPeriod;						// period

    // Output objects
    TList*			fListOfObjects;			//! Output list of objects
    TList*			fListOfObjectsMC;		//! Output list of objects for MonteCarlo

    // Cut functions for LHC15o
    TF1*				fMultTOFLowCut;			// cut low for TOF multiplicity outliers
    TF1*				fMultTOFHighCut;		// cut high for TOF multiplicity outliers
    TF1*				fMultCentLowCut;		// cut low for multiplicity centrality outliers

    // NUE
    TFile*			fTrackEfficiency;		//! file with tracking efficiency
    TH3F*				hTrackEfficiency;		//! histogram with tracking efficiency
    TH3F*				hTrackEfficiencyRun;//! histogram with tracking efficiency

    // NUA
    TFile*			fPhiWeight;					//! file with phi weights
    TFile*			fPhiWeightPlus;			//! file with phi weights
    TFile*			fPhiWeightMinus;		//! file with phi weights
    TH3F*				hPhiWeight;					//! 3D weight for all periods except LHC15ijl
    TH3F*				hPhiWeightRun;			//! 3D weight run-by-run for pPb 5TeV LHC16q
    TH1F*				hPhiWeight1D;				//! 1D weight in one MC case (maybe need to redo to 3D weight)
    TH3D*				hPhiWeight_LHC15i_part1;	//! LHC15i, part1 runs
    TH3D*				hPhiWeight_LHC15i_part2;	//! LHC15i, part2 runs
    TH3D*				hPhiWeight_LHC15j_part1;	//! LHC15j, part1 runs
    TH3D*				hPhiWeight_LHC15j_part2;	//! LHC15j, part2 runs
    TH3D*				hPhiWeight_LHC15l_part1;	//! LHC15l, part1 runs
    TH3D*				hPhiWeight_LHC15l_part2;	//! LHC15l, part2 runs
    TH3D*				hPhiWeight_LHC15l_part3;	//! LHC15l, part3 runs

    TH3D*				hPhiWeightPlus_LHC15i_part1;	//! LHC15i, part1 runs
    TH3D*				hPhiWeightPlus_LHC15i_part2;	//! LHC15i, part2 runs
    TH3D*				hPhiWeightPlus_LHC15j_part1;	//! LHC15j, part1 runs
    TH3D*				hPhiWeightPlus_LHC15j_part2;	//! LHC15j, part2 runs
    TH3D*				hPhiWeightPlus_LHC15l_part1;	//! LHC15l, part1 runs
    TH3D*				hPhiWeightPlus_LHC15l_part2;	//! LHC15l, part2 runs
    TH3D*				hPhiWeightPlus_LHC15l_part3;	//! LHC15l, part3 runs

    TH3D*				hPhiWeightMinus_LHC15i_part1;	//! LHC15i, part1 runs
    TH3D*				hPhiWeightMinus_LHC15i_part2;	//! LHC15i, part2 runs
    TH3D*				hPhiWeightMinus_LHC15j_part1;	//! LHC15j, part1 runs
    TH3D*				hPhiWeightMinus_LHC15j_part2;	//! LHC15j, part2 runs
    TH3D*				hPhiWeightMinus_LHC15l_part1;	//! LHC15l, part1 runs
    TH3D*				hPhiWeightMinus_LHC15l_part2;	//! LHC15l, part2 runs
    TH3D*				hPhiWeightMinus_LHC15l_part3;	//! LHC15l, part3 runs

    // Event histograms
    TH1D*				hEventCount;								//! counting events passing given event cuts
    TH1F*				hMult;											//! multiplicity distribution
    TH1F*				hMultfBin[12]; 							//! multiplicity distribution in fBin
    TH1F*				fVtxAfterCuts;							//! Vertex z dist after cuts
    TH1F*				fCentralityDis;							//! distribution of centrality percentile using V0M estimator
    TH1F*				fV0CentralityDis;						//! distribution of V0M/<V0M>
    TH2F*				hMultV0vsNtrksAfterCuts;		//! Number of tracks vs. V0M/<V0M>
    TH2F*				hMultSPDvsNtrksAfterCuts;		//! Number of tracks vs. SPD/<SPD>
    TH2F*				hNtrksVSmultPercentile; 		//! Number of tracks vs. percentile using V0M estimator
    TH2F*				fCentralityV0MCL1;					//! LHC15o: V0M vs. CL1 percentile
    TH2F*				fCentralityV0MCL0;					//! LHC15o: V0M vs. CL0 percentile
    TH2F*				fCentralityCL0CL1;					//! LHC15o: CL0 vs. CL1 percentile
    TH2F*				fMultvsCentr;								//! LHC15o: Number of tracks vs. percentile
    TH2F*				fMult128vsCentr;						//! LHC15o: Number of FB128 tracks vs. percentile
    TH2F*				fMultTPCvsTOF;							//! LHC15o: Number of TPC tracks vs. ToF tracks
    TH2F*				fMultTPCvsESD;							//! LHC15o: Number of TPC tracks vs. ESD tracks

    TH2D*				hSPDClsVsTrk;								//! SPD clusters vs. tracklets without any cuts
    TH2D*				hV0C012vsTkl;								//! V0C mult. in 0,1,2nd ring vs. SPD tracklets without any cuts
    TH2D*				hV0C012vsV0C3;							//! V0C mult. in 0,1,2nd ring vs. V0C mult. in 3rd ring without any cuts
    TH2D*				hV0MOnVsOf;									//! V0M amplitude online vs. offline without any cuts
    TH2D*				hSPDOnVsOf;									//! SPD amplitude online vs. offline without anycuts


    // Track histograms
    TH1F*				fPhiDis1D;									//! phi dis 1D
    TH3F*				fPhiDis;										//! phi dist
    TH1F*				fEtaDis;										//! eta dist
    TH1F*				fEtaBefore;									//! eta dist before track cuts
    TH1F*				fPtDis;											//! pt dist
    TH1F*				fPtBefore;									//! pt dist before track cuts
    TH1F*				hDCAxyBefore; 							//!
    TH1F*				hDCAzBefore; 								//!
    TH1F*				hITSclustersBefore; 				//!
    TH1F*				hChi2Before; 								//!
    TH1F*				hDCAxy; 										//!
    TH1F*				hDCAz; 											//!
    TH1F*				hITSclusters; 							//!
    TH1F*				hChi2; 											//!

    // Correlation histograms
    TH2F*				hNtrksPt0530Pt0230; 				//! for correction of Xaxis of results with 0.5<pt<3.0 to 0.2<pt<3.0
    TH2F*				hNtrksPt0730Pt0230; 				//! for correction of Xaxis of results with 0.7<pt<3.0 to 0.2<pt<3.0
    TH2F*				hNtrksEta09Eta10;						//! for correction of Xaxis of results with |eta|<0.9 to |eta|<1.0
    TH2F*				hNtrksEta08Eta10;						//! for correction of Xaxis of results with |eta|<0.8 to |eta|<1.0
    TH2F*				hNtrksAllNtrksLS; 					//! for correction of Xaxis of LS results to all charged ptcls results
    TH2F*				hNtrksNoGapGap0;						//! for correction of Xaxis of Gap0 to NoGap results
    TH2F*				hNtrksNoGapGap2;						//! for correction of Xaxis of Gap2 to NoGap results
    TH2F*				hNtrksNoGapGap4;						//! for correction of Xaxis of Gap4 to NoGap results
    TH2F*               hNtrksNoGapGap6;                        //! for correction of Xaxis of Gap6 to NoGap results
    TH2F*				hNtrksNoGapGap8;						//! for correction of Xaxis of Gap8 to NoGap results
    TH2F*				hNtrksNoGapGap;							//! for correction of Xaxis of Gap10 to NoGap results
    TH2F*				hNtrksNoGapGap14;						//! for correction of Xaxis of Gap14 to NoGap results
    TH2F*				hNtrksNoGapGap16;						//! for correction of Xaxis of Gap16 to NoGap results
    TH2F*				hNtrksNoGapGap18;						//! for correction of Xaxis of Gap18 to NoGap results
    TH2F*				hNtrksNoGap3sub;						//! for correction of Xaxis of 3sub to NoGap results
    TH2F*				hNtrksNoGap3subGap;					//! for correction of Xaxis of 3sub Gap2 to NoGap results

    // Global variables
    int NtrksAfter = 0;
    int NtrksAfterGap10M = 0;
    int NtrksAfterGap10P = 0;
    int NtrksAfterGap14M = 0;
    int NtrksAfterGap14P = 0;
    int NtrksAfter3subL = 0;
    int NtrksAfter3subM = 0;
    int NtrksAfter3subR = 0;

    PhysicsProfile centProfile;
    PhysicsProfile centProfile_bin[10];
    PhysicsProfile multProfile;
    PhysicsProfile multProfile_bin[10];

    //......
    // MC
    // Event histograms
    TH1F*				hMultMC;										//! Multiplicity distribution
    TH1F*				fPhiDisTruth;								//! Phi distribution
    TH1F*				fEtaDisTruth;								//! Eta distribution
    TH1F*				fPtDisTruth;								//! Pt distribution

    TH3F*				hReco;											//! pt, eta, Vz of reconstructed tracks
    TH3F*				hRecoPion;									//! pt, eta, Vz of reconstructed pions
    TH3F*				hRecoKaon;									//! pt, eta, Vz of reconstructed kaons
    TH3F*				hRecoProton;								//! pt, eta, Vz of reconstructed protons
    TH3F*				hRecoElectron;							//! pt, eta, Vz of reconstructed electrons
    TH3F*				hRecoMuon;									//! pt, eta, Vz of reconstructed muons
    TH3F*				hRecoLSplus;								//! pt, eta, Vz of positive reconstructed tracks
    TH3F*				hRecoLSminus;								//! pt, eta, Vz of negative reconstructed tracks

    TH2F*				hPtRecoNtrks;								//! reconstructed Pt vs. Number of generated particles
    TH2F*				hEtaRecoNtrks;							//! reconstructed Eta vs. Number of generated particles
    TH2F*				hVzRecoNtrks;								//! reconstructed Vz vs. Number of generated particles
    TH2F*				hPtRecoNtrksReco;						//! reconstructed Pt vs. Number of reconstructed tracks
    TH2F*				hEtaRecoNtrksReco;					//! reconstructed Eta vs. Number of reconstructed tracks
    TH2F*				hVzRecoNtrksReco;						//! reconstructed Vz vs. Number of reconstructed tracks

    TH3F*				hTruth;											//! pt, eta, Vz of generated particles
    TH3F*				hTruthPion;									//! pt, eta, Vz of generated pions
    TH3F*				hTruthKaon;									//! pt, eta, Vz of generated kaons
    TH3F*				hTruthProton;								//! pt, eta, Vz of generated protons
    TH3F*				hTruthElectron;							//! pt, eta, Vz of generated electrons
    TH3F*				hTruthMuon;									//! pt, eta, Vz of generated muons
    TH3F*				hTruthLSplus;								//! pt, eta, Vz of positive generated particles
    TH3F*				hTruthLSminus;							//! pt, eta, Vz of negative generated particles

    TH2F*				hPtTruthNtrks;							//! generated Pt vs. Number of generated particles
    TH2F*				hEtaTruthNtrks;							//! generated Eta vs. Number of generated particles
    TH2F*				hVzTruthNtrks;							//! generated Vs vs. Number of generated particles
    TH2F*				hPtTruthNtrksReco;					//! generated Pt vs. Number of reconstructed tracks
    TH2F*				hEtaTruthNtrksReco;					//! generated Eta vs. Number of reconstructed tracks
    TH2F*				hVzTruthNtrksReco;					//! generated Vz vs. Number of reconstructed tracks

    TH1D*				hPrimary;										//! number of primary particles
    TH1D*				hPions;											//! number of pions

    TH2F*				hDCAptMC;										//! DCA vs. Pt
    TH2F*				hDCAptMC_material;					//! DCA vs. Pt of decays from material interactions
    TH2F*				hDCAptMC_weak;							//! DCA vs. Pt of weak decays

    TH2F*				hNtrksRecoNtrksTruth; 			//! Number of generated particles vs. number of reconstructed tracks
    TH2F*				hNtrksRecoCorrNtrksTruth; 	//! Number of generated particles vs. number of reconstructed tracks weighted to match the difference between data and MC seen in LHC15ijl

    TH3F*				hRecoMult[8]; 							//!
    TH3F*				hTruthMult[8]; 							//!

    TH1F*				hTruthMultfBin[12]; 				//!

    TRandom3 rand;
    Int_t bootstrap_value;

    void CalculateProfile(PhysicsProfile& profile, double Ntrks);
    void InitCentralityProfile(PhysicsProfile& profile, TString);
    void InitMultiplicityProfile(PhysicsProfile& profile, TString);

    ClassDef(AliAnalysisTaskNonlinearFlow, 1);    //Analysis task
};

#endif
// Local Variables:
//  mode: C++
// End:
