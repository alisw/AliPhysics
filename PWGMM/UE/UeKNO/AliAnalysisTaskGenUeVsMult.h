/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskGenUeVsMult_H
#define AliAnalysisTaskGenUeVsMult_H

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



class AliAnalysisTaskGenUeVsMult : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskGenUeVsMult();
	AliAnalysisTaskGenUeVsMult(const char *name);
	
    virtual                 ~AliAnalysisTaskGenUeVsMult();

	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);
    
	void       GetMultLeadingObject();
	//void       GetGenUEObservables();
    
    void       GetMeanGenUEObservables(std::vector<Double_t> &gen);

  //  bool       HasRecVertex();
    virtual    Double_t DeltaPhi(Double_t phia, Double_t phib,Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
    
    
    void FillHisto(const char* objkey, Double_t x);
    
    
	void       SetPtMin(Double_t val)              {fPtMin = val;}   // use differnet ptcuts
	
	//virtual void SetGenerator(TString generator){fGenerator=generator;}

protected:



private:
    Int_t   GetMultipliciy(std::vector<Int_t> &mult,std::vector<Int_t> &regionArray);
    Bool_t IsMCEventSelected(TObject* obj);
	Int_t    GetPidCode(Int_t pdgCode) const;
    void GetMultiVsUEObservables(std::vector<Int_t> &mult,std::vector<Int_t> &region);
   
    AliMCEvent*  fMC;                                               //! MC Event
    AliInputEventHandler*    fMcHandler;  //!<!
	
    AliStack*    fMCStack;                                                 //! MC stack
	
    Double_t fGenLeadPhi;
    Double_t fGenLeadPt;
    Int_t    fGenLeadIn;

	TList*                  fOutputList;                                      //! output list in the root file

    Double_t fEtaCut;
    Double_t fPtMin;
    
    //TString     fGenerator;
	
	
	
	
	// UE 
	
    TH1I * fHistEvt;         //!<!     QA of event properties
	TH1D * hCounter;
    TH1D * hPtLeadingGenAll;
    TH1D * hPtLeadingTrue;
    TH1D * fDeltaphiNS;
    TH1D * fDeltaphiAS;
    TH1D * fDeltaphiTS;
	TH1D * fMultiplicyNS;
    TH1D * fMultiplicyAS;
    TH1D * fMultiplicyTS;
    
    //TH1D * fMult2[5];
	

	TH2D * hPtVsPtLeadingTrue[3];
    
	TProfile * pNumDenTrueAll[3];
    TProfile * pSumPtTrueAll[3];
    
    TH2D * hNumDen[3];
    TH2D * hSumPt[3];
    
	TProfile * pNumDenTrue[3];
	
	TProfile * pSumPtTrue[3];
    
    TH1D * fMult[4];
    TH3D * fHistPtLeadingVsNchNS[4];
    TH3D * fHistPtLeadingVsNchAS[4];
    TH3D * fHistPtLeadingVsNchTS[4];
	

	AliAnalysisTaskGenUeVsMult(const AliAnalysisTaskGenUeVsMult&);                  // not implemented
	AliAnalysisTaskGenUeVsMult& operator=(const AliAnalysisTaskGenUeVsMult&);       // not implemented

	ClassDef(AliAnalysisTaskGenUeVsMult, 3);
};

#endif
