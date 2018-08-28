/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskDeuteronAbsorption_H
#define AliAnalysisTaskDeuteronAbsorption_H

#include "AliAnalysisTaskSE.h"

class AliESDEvent;
class AliPIDResponse;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliESDtrackCuts;
// for Monte Carlo:
class AliMCEventHandler;
class AliMCEvent;
class AliStack;
//class AliVParticle;
class AliGenEventHeader;



class AliAnalysisTaskDeuteronAbsorption : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskDeuteronAbsorption();
                                AliAnalysisTaskDeuteronAbsorption(const char *name);
        virtual                 ~AliAnalysisTaskDeuteronAbsorption();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliESDEvent*            fESD;           //! input event
        AliPIDResponse*         fPIDResponse;           //! pid response
	//
        AliESDtrackCuts* fESDtrackCuts; //! input track cuts
	//
        TList*                  fOutputList;    //! output list
	//
	TH1F *fHistZv; //!
	TH2F* fHist2PIDvP; //!
	TH2F* fHist2PIDvka; //!
	TH2F* fHist2PIDvDe; //!
        TH2F* fHist2PIDvTr; //!
	TH2F* fHist2PIDv; //!
	TH2F* fHist2PIDf; //!
	TH2F* fHist2PIDDef; //!
	TH3F* fHist3PIDvDe;  //!
        TH3F* fHist3PIDvTr;  //!
	TH3F* fHist3PIDv;    //!
	TH3F* fHist3PIDf;    //!
	TH3F* fHist3PIDDef;  //!
	TH1F* fHistPhi;      //!
	TH2F* fHistPhi2; //!
        TH2F* fHistPhi2n; //!
        TH2F* fHistPhi2o; //!
        TH2F* fHistPhi2no; //!
        TH2F* fHistmass;//!
	TH2F* fHistmassDe;  //!
        TH2F* fHistmassTr;  //!
	TH2F* fHistmassP;  //!
	TH2F* fHistmassDei; //!
	TH2F* fHistmassDeo; //!
	TH2F* fHistmassPri; //!
        TH2F* fHistmassPro; //!
        TH2F* fHistmassTri; //!
        TH2F* fHistmassTro; //!

        TH2F* fHistMatchAllDeuteronPos; //!
	TH2F* fHistMatchTofDeuteronPos;  //!
	TH2F* fHistMatchAllDeuteronNeg; //!
	TH2F* fHistMatchTofDeuteronNeg; //!

        TH2F* fHistMatchAllDeuteronPoso; //!
        TH2F* fHistMatchTofDeuteronPoso;  //!
        TH2F* fHistMatchAllDeuteronNego; //!
        TH2F* fHistMatchTofDeuteronNego; //!

        TH2F* fHistMatchAllDeuteronPosMC; //!
        TH2F* fHistMatchTofDeuteronPosMC;  //!
        TH2F* fHistMatchAllDeuteronNegMC; //!
        TH2F* fHistMatchTofDeuteronNegMC; //!
        
        TH2F* fHistMatchAllDeuteronPosMCo; //!
        TH2F* fHistMatchAllDeuteronNegMCo; //!
        TH2F* fHistMatchTofDeuteronPosMCo; //!
        TH2F* fHistMatchTofDeuteronNegMCo; //!
          
        TH3F* fHistMatchAllTritonPos; //!
        TH3F* fHistMatchTofTritonPos;  //!
        TH3F* fHistMatchAllTritonNeg; //!
        TH3F* fHistMatchTofTritonNeg; //!

        TH2F* fHistMatchAllTritonPosMC; //!
        TH2F* fHistMatchTofTritonPosMC;  //!
        TH2F* fHistMatchAllTritonNegMC; //!
        TH2F* fHistMatchTofTritonNegMC; //!

       	TH2F* fHistMatchAllProtonPos; //!
	TH2F* fHistMatchTofProtonPos; //!
	TH2F* fHistMatchAllProtonNeg; //!
	TH2F* fHistMatchTofProtonNeg; //!
        
        TH2F* fHistMatchAllProtonPoso; //!
        TH2F* fHistMatchTofProtonPoso; //!
        TH2F* fHistMatchAllProtonNego; //!
        TH2F* fHistMatchTofProtonNego; //!

        
        TH2F* fHistMatchAllProtonPosMC; //!
        TH2F* fHistMatchTofProtonPosMC;  //!
        TH2F* fHistMatchAllProtonNegMC; //!
        TH2F* fHistMatchTofProtonNegMC; //!
         
        TH2F* fHistMatchAllProtonPosMCo; //!
        TH2F* fHistMatchTofProtonPosMCo;  //!
        TH2F* fHistMatchAllProtonNegMCo; //!
        TH2F* fHistMatchTofProtonNegMCo; //!
 
        
        TH1F *      hptRecoDeut;     //!

  	TH1F *      hptRecoAntiDeut;  //!

     	TH1F *	hptMatchDeut; //!

  	TH1F *	hptMatchAntiDeut; //!

 	TH1F *	hptGoodMatchDeut; //!

  	TH1F *	hptGoodMatchAntiDeut; //!

        TH1F *      hptRecoProt;     //!

        TH1F *      hptRecoAntiProt;  //!

        TH1F *  hptMatchProt; //!

        TH1F *  hptMatchAntiProt; //!
        TH1F *  hptGoodMatchProt; //!

        TH1F *  hptGoodMatchAntiProt; //!

	TF1*    fTRDboundariesPos[4];       //! Function with the phi limits of TRD boundaries as a function of pt
        TF1*    fTRDboundariesNeg[4];       //! Function with the phi limits of TRD boundaries as a function of pt

        TH2F * fHistMatchAllDeuteronPosTPCsigma; //!
        TH2F * fHistMatchAllDeuteronPosTPCsigmao; //!
        TH2F * fHistMatchAllDeuteronNegTPCsigma; //!
        TH2F * fHistMatchAllDeuteronNegTPCsigmao; //!



	//
        AliAnalysisTaskDeuteronAbsorption(const AliAnalysisTaskDeuteronAbsorption&); // not implemented
        AliAnalysisTaskDeuteronAbsorption& operator=(const AliAnalysisTaskDeuteronAbsorption&); // not implemented
	//
        ClassDef(AliAnalysisTaskDeuteronAbsorption, 1);
      
};

#endif
