/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskGenMcKnoUe_H
#define AliAnalysisTaskGenMcKnoUe_H

class AliESDtrackCuts;
class AliESDEvent;
class TList;
class TH1D;
class TH2D;
class TH3D;
class TH1I;
class TProfile;
class THnSparse;

#include "AliAnalysisTaskSE.h"


#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliGenEventHeader.h"



class AliAnalysisTaskGenMcKnoUe : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskGenMcKnoUe();
	AliAnalysisTaskGenMcKnoUe(const char *name);
	
    virtual                 ~AliAnalysisTaskGenMcKnoUe();

	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);
    
    Bool_t IsMCEventSelected(TObject* obj);
	void FillHisto(const char* objkey, Double_t x);
    
   void       GetMeanGenUEObservables(std::vector<Double_t> &gen);
    
    void       GetGenLeadingObject();
	void       GetGenUEObservables();
    
  //  bool       HasRecVertex();
    virtual    Double_t DeltaPhi(Double_t phia, Double_t phib,Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
    
	void       SetPtMin(Double_t val)              {fPtMin = val;}   // use differnet ptcuts
	
	//virtual void SetGenerator(TString generator){fGenerator=generator;}

protected:



private:
    
	
    AliMCEvent*  fMC;                                               //! MC Event
    AliInputEventHandler*    fMcHandler;  //!<!
	AliStack*    fMCStack;                                                 //! MC stack
	Double_t fEtaCut;
    Double_t fPtMin;
	TList*                  fOutputList;                                      //! output list in the root file

   // TString     fGenerator;
	
	Double_t fGenLeadPhi; 
	Double_t fGenLeadPt;
	Int_t    fGenLeadIn;
	
	// UE 
	
    TH1I * fHistEvt;         //!<!     QA of event properties
	TH1D * hCounter;
    TH1D * hPtLeadingGenAll;
    TH1D * hPtLeadingTrue;
    
	TH2D * hPtVsPtLeadingTrue[3];
	
    TProfile * pNumDenTrueAll[3];
    TProfile * pSumPtTrueAll[3];

    TH2D * hNumDen[3];
    TH2D * hSumPt[3];
	
    TProfile * pNumDenTrue[3];
    TProfile * pSumPtTrue[3];

	AliAnalysisTaskGenMcKnoUe(const AliAnalysisTaskGenMcKnoUe&);                  // not implemented
	AliAnalysisTaskGenMcKnoUe& operator=(const AliAnalysisTaskGenMcKnoUe&);       // not implemented

	ClassDef(AliAnalysisTaskGenMcKnoUe, 3);
};

#endif
