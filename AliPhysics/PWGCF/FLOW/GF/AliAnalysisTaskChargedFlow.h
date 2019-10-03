#ifndef ALIANALYSISTASKCHARGEDFLOW_H
#define ALIANALYSISTASKCHARGEDFLOW_H

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include <TComplex.h>

// STL includes
#include <iostream>
using std::cout;
using std::endl;

class AliAnalysisTaskChargedFlow : public AliAnalysisTaskSE {
	public:
  	AliAnalysisTaskChargedFlow();
  	AliAnalysisTaskChargedFlow(const char *name);
  	virtual ~AliAnalysisTaskChargedFlow();

  	virtual void   UserCreateOutputObjects();
  	virtual void   UserExec(Option_t *);
  	virtual void   Terminate(const Option_t*); 

  	virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  	virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  	virtual void  SetMinPt(Double_t minPt){fMinPt = minPt;}
  	virtual void  SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
  	virtual void  SetIsSample(Int_t IsSample){fSample = IsSample;}
			
  	virtual void		AnalyzeAOD();  

 protected:
  TComplex Qvector[20][20];
  TComplex QvectorM[20][20];
  TComplex QvectorP[20][20];
  TComplex Qvector0M[20][20];
  TComplex Qvector0P[20][20];
  TComplex Qvector2M[20][20];
  TComplex Qvector2P[20][20];
  TComplex Qvector4M[20][20];
  TComplex Qvector4P[20][20];
  TComplex Qvector8M[20][20];
  TComplex Qvector8P[20][20];
  TComplex Qvector14M[20][20];
  TComplex Qvector14P[20][20];
	TComplex QvectorSubLeft[20][20];
	TComplex QvectorSubMiddle[20][20];
	TComplex QvectorSubRight[20][20];
	
  TComplex Q(int n, int p);
  TComplex QGapM(int n, int p);
  TComplex QGapP(int n, int p);
  TComplex QGap0M(int n, int p);
  TComplex QGap0P(int n, int p);
  TComplex QGap2M(int n, int p);
  TComplex QGap2P(int n, int p);
  TComplex QGap4M(int n, int p);
  TComplex QGap4P(int n, int p);
  TComplex QGap8M(int n, int p);
  TComplex QGap8P(int n, int p);
  TComplex QGap14M(int n, int p);
  TComplex QGap14P(int n, int p);
  TComplex QGap16M(int n, int p);
  TComplex QGap16P(int n, int p);
	TComplex QsubLeft(int n, int p);
	TComplex QsubMiddle(int n, int p);
	TComplex QsubRight(int n, int p);

  TComplex Two(int n1, int n2);
  TComplex TwoGap(int n1, int n2);
  TComplex TwoGap0(int n1, int n2);
  TComplex TwoGap2(int n1, int n2);
  TComplex TwoGap4(int n1, int n2);
  TComplex TwoGap8(int n1, int n2);
  TComplex TwoGap14(int n1, int n2);
  TComplex Two_3SubLM(int n1, int n2);
  TComplex Two_3SubRM(int n1, int n2);
  TComplex Two_3SubLM_perm2(int n1, int n2);
  TComplex Two_3SubLR_perm2(int n1, int n2);
  TComplex Two_3SubRM_perm3(int n1, int n2);
  TComplex Two_3SubRL_perm3(int n1, int n2);

  TComplex Three(int n1, int n2, int n3);
  TComplex ThreeGap0M(int n1, int n2, int n3);
  TComplex ThreeGap0P(int n1, int n2, int n3);
  TComplex Four(int n1, int n2, int n3, int n4);
	TComplex FourGap0(int n1, int n2, int n3, int n4);
	TComplex FourGap0M(int n1, int n2, int n3, int n4);
	TComplex FourGap0P(int n1, int n2, int n3, int n4);
	TComplex Four_3SubEvts(int n1, int n2, int n3, int n4);
	TComplex Four_3SubEvts_perm2(int n1, int n2, int n3, int n4);
	TComplex Four_3SubEvts_perm3(int n1, int n2, int n3, int n4);
	TComplex Five(int n1, int n2, int n3, int n4, int n5);
	TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6);
	TComplex SixGap0(int n1, int n2, int n3, int n4, int n5, int n6);
	TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
	TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
	TComplex EightGap0(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
  
  // Cuts and options
	Double_t		fEtaCut;						// Eta cut used to select particles
	Double_t		fVtxCut;						// Vtx cut on z position in cm
	Double_t		fMinPt;							// Min pt - for histogram limits
	Double_t		fMaxPt;							// Max pt - for histogram limits
	Int_t				fSample;						// number of sample
 
  // Output objects
	TList*			fOutputContainer;			//! Output list of objects

	// Event histograms
  TH2F*				hMult;											//! multiplicity distribution
  TH1F*				hMultTrigger;								//! multiplicity distribution
	TH1F*				fVtxAfterCuts;							//! Vertex z dist after cuts

	// Track histograms
 	TH2F*				fPhiDis;										//! phi dist
	TH2F*				fEtaDis;										//! eta dist
	TH2F*				fEtaBefore;									//! eta dist before track cuts
	TH2F*				fPtDis;											//! pt dist
	TH2F*				fPtBefore;									//! pt dist before track cuts
	TH1F*				hIsPiKp;										//!

	// Physics profiles
	TProfile2D*		fChsc4242;									//! SC(4,2)
	TProfile2D*		fChsc4242Gap0;							//! SC(4,2) |#Delta#eta| > 0.0
	TProfile2D*		fChsc4242_3sub;							//! SC(4,2)_A 3subevent method
	TProfile2D*		fChsc4224_3sub;							//! SC(4,2)_B 3subevent method
	TProfile2D*		fChsc4242_3sub_perm2;				//! SC(4,2)_A 3subevent method
	TProfile2D*		fChsc4224_3sub_perm2;				//! SC(4,2)_B 3subevent method
	TProfile2D*		fChsc4242_3sub_perm3;				//! SC(4,2)_A 3subevent method
	TProfile2D*		fChsc4224_3sub_perm3;				//! SC(4,2)_B 3subevent method

	TProfile2D*		fChsc3232;									//! SC(3,2)
	TProfile2D*		fChsc3232Gap0;							//! SC(3,2) |#Delta#eta| > 0.0
	TProfile2D*		fChsc3232_3sub;							//! SC(3,2)_A 3subevent method
	TProfile2D*		fChsc3223_3sub;							//! SC(3,2)_B 3subevent method
	TProfile2D*		fChsc3232_3sub_perm2;				//! SC(3,2)_A 3subevent method
	TProfile2D*		fChsc3223_3sub_perm2;				//! SC(3,2)_B 3subevent method
	TProfile2D*		fChsc3232_3sub_perm3;				//! SC(3,2)_A 3subevent method
	TProfile2D*		fChsc3223_3sub_perm3;				//! SC(3,2)_B 3subevent method

	TProfile2D*		fsc4242[12];								//! the same but for different fBin (sampling)
	TProfile2D*		fsc4242Gap0[12];						//!
	TProfile2D*		fsc4242_3sub[12];						//!	
	TProfile2D*		fsc4224_3sub[12];						//!
	TProfile2D*		fsc4242_3sub_perm2[12];			//!	
	TProfile2D*		fsc4224_3sub_perm2[12];			//!
	TProfile2D*		fsc4242_3sub_perm3[12];			//!	
	TProfile2D*		fsc4224_3sub_perm3[12];			//!

	TProfile2D*		fsc3232[12];								//!
	TProfile2D*		fsc3232Gap0[12];						//!
	TProfile2D*		fsc3232_3sub[12];						//!
	TProfile2D*		fsc3223_3sub[12];						//!
	TProfile2D*		fsc3232_3sub_perm2[12];			//!
	TProfile2D*		fsc3223_3sub_perm2[12];			//!
	TProfile2D*		fsc3232_3sub_perm3[12];			//!
	TProfile2D*		fsc3223_3sub_perm3[12];			//!

	// Standard correlation profiles for different harmonics
  TProfile2D*		fChcn2Ntrks1bin[6];  				//! <<2>> in unit bins of Ntrks
  TProfile2D*		fChcn2Gap0Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.0
  TProfile2D*		fChcn2Gap2Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.2
  TProfile2D*		fChcn2Gap4Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.4
  TProfile2D*		fChcn2Gap8Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.8
  TProfile2D*		fChcn2GapNtrks1bin[6];  		//! <<2>> |#Delta#eta| > 1.0
  TProfile2D*		fChcn2Gap14Ntrks1bin[6];  	//! <<2>> |#Delta#eta| > 1.4

  TProfile2D*		fChcn2_3subLMNtrks1bin[6];  		//! <<2>> left vs. middle subevent
  TProfile2D*		fChcn2_3subRMNtrks1bin[6];  		//! <<2>> middle vs. right subevent
  TProfile2D*		fChcn2_3subLM_perm2Ntrks1bin[6];//! <<2>> left vs. middle subevent
  TProfile2D*		fChcn2_3subLR_perm2Ntrks1bin[6];//! <<2>> left vs. right subevent
  TProfile2D*		fChcn2_3subRM_perm3Ntrks1bin[6];//! <<2>> right vs. middle subevent
  TProfile2D*		fChcn2_3subRL_perm3Ntrks1bin[6];//! <<2>> right vs. left subevent

  TProfile2D*		fChcn4Ntrks1bin[6];  				//! <<4>> in unit bins of Ntrks
  TProfile2D*		fChcn4Gap0Ntrks1bin[6];  		//! <<4>> |#Delta#eta| > 0.0
	TProfile2D*		fChcn4_3subNtrks1bin[6]; 		//! <<4>> 3subevent method
	TProfile2D*		fChcn4_3sub_perm2Ntrks1bin[6];		//!
	TProfile2D*		fChcn4_3sub_perm3Ntrks1bin[6];		//!

  TProfile2D*		fChcn6Ntrks1bin[6];  				//! <<6>> in unit bins of Ntrks
  TProfile2D*		fChcn6Gap0Ntrks1bin[6];  		//! <<6>> |#Delta#eta| > 0.0

  TProfile2D*		fChcn8Ntrks1bin[6];  				//! <<8>> in unit bins of Ntrks
  TProfile2D*		fChcn8Gap0Ntrks1bin[6];  		//! <<8>> |#Delta#eta| > 0.0

	// the same profiles but for different fBin (sampling)
  TProfile2D*		fcn2Ntrks1bin[6][12];  			//!
  TProfile2D*		fcn2Gap0Ntrks1bin[6][12];		//!
  TProfile2D*		fcn2Gap2Ntrks1bin[6][12];		//!
  TProfile2D*		fcn2Gap4Ntrks1bin[6][12];		//!
  TProfile2D*		fcn2Gap8Ntrks1bin[6][12];		//!
  TProfile2D*		fcn2GapNtrks1bin[6][12];		//!
  TProfile2D*		fcn2Gap14Ntrks1bin[6][12];	//!

  TProfile2D*		fcn2_3subLMNtrks1bin[6][12];  		//!
  TProfile2D*		fcn2_3subRMNtrks1bin[6][12];  		//!
  TProfile2D*		fcn2_3subLM_perm2Ntrks1bin[6][12];//!
  TProfile2D*		fcn2_3subLR_perm2Ntrks1bin[6][12];//!
  TProfile2D*		fcn2_3subRM_perm3Ntrks1bin[6][12];//!
  TProfile2D*		fcn2_3subRL_perm3Ntrks1bin[6][12];//!

  TProfile2D*		fcn4Ntrks1bin[6][12];				//!
  TProfile2D*		fcn4Gap0Ntrks1bin[6][12];		//!
	TProfile2D*		fcn4_3subNtrks1bin[6][12];	//!
	TProfile2D*		fcn4_3sub_perm2Ntrks1bin[6][12];	//!
	TProfile2D*		fcn4_3sub_perm3Ntrks1bin[6][12];	//!

  TProfile2D*		fcn6Ntrks1bin[6][12];				//!
  TProfile2D*		fcn6Gap0Ntrks1bin[6][12];		//!

  TProfile2D*		fcn8Ntrks1bin[6][12];				//!
  TProfile2D*		fcn8Gap0Ntrks1bin[6][12];		//!

 private:
  AliAnalysisTaskChargedFlow(const AliAnalysisTaskChargedFlow&);
  AliAnalysisTaskChargedFlow& operator=(const AliAnalysisTaskChargedFlow&);

  ClassDef(AliAnalysisTaskChargedFlow, 2);    //Analysis task
};

#endif
// Local Variables: 
//  mode: C++
// End:
