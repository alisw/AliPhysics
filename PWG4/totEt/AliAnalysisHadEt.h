//_________________________________________________________________________
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//
// This class is designed for the analysis of the hadronic component of 
// transverse energy.  It is used by AliAnalysisTaskHadEt.
//_________________________________________________________________________
#ifndef ALIANALYSISHADET_H
#define ALIANALYSISHADET_H

#include "TString.h"
#include "AliAnalysisEtCommon.h"

class TH2F;
class TH1F;
class AliVEvent;
class TList;
class AliESDtrackCuts;
class Rtypes;
class TParticle;
class TDatabasePDG;
class AliAnalysisEtCuts;

class AliAnalysisHadEt : public AliAnalysisEtCommon
{
public:
   
  AliAnalysisHadEt();
    virtual ~AliAnalysisHadEt();

    /** Analyse the event! */
    virtual Int_t AnalyseEvent(AliVEvent *event);


    /** Initialise the analysis, must be overloaded. */
    virtual void Init();

    /** Reset event specific values (Et etc.) */
    virtual void ResetEventValues();


    /** Cuts info */
    AliAnalysisEtCuts * GetCuts() const { return fCuts; } 
    virtual void SetCuts(const AliAnalysisEtCuts *cuts) 
    { fCuts = (AliAnalysisEtCuts *) cuts; } 
    
    /** Sum of the total Et for all events */
    Double_t GetSumEt() const { return fSumEt; }

    /** Sum of the total Et within our acceptance for all events */
    Double_t GetSumEtAcc() const { return fSumEtAcc; }

    /** Total Et in the event (without acceptance cuts) */
    Double_t GetTotEt() const { return fTotEt; }

    /** Total Et in the event within the acceptance cuts */
    Double_t GetTotEtAcc() const { return fTotEtAcc; }

   /** Total neutral Et in the event (without acceptance cuts) */
    Double_t GetTotNeutralEt() const { return fTotNeutralEt; }

    /** Total neutral Et in the event within the acceptance cuts */
    Double_t GetTotNeutralEtAcc() const { return fTotNeutralEtAcc; }
    
    /** Total charged Et in the event (without acceptance cuts) */
    Double_t GetTotChargedEt() const { return fTotChargedEt; }

    /** Total charged Et in the event within the acceptance cuts */
    Double_t GetTotChargedEtAcc() const { return fTotChargedEtAcc; }


    void SetHistoList(const TList *mylist){fhistoList = (TList *) mylist;}


protected:   
    
    Double_t fSumEt;/** Sum of the total Et for all events */
    Double_t fSumEtAcc;/** Sum of the total Et within our acceptance for all events */
    Double_t fTotEt;/** Total Et in the event (without acceptance cuts) */
    Double_t fTotEtAcc;/** Total Et in the event within the acceptance cuts */
    
    Double_t fTotNeutralEt;/** Total neutral Et in the event */
    Double_t fTotNeutralEtAcc;/** Total neutral Et in the event within the acceptance cuts */
    Double_t fTotChargedEt;/** Total charged Et in the event */
    Double_t fTotChargedEtAcc;/** Total charged Et in the event within the acceptance cuts */

    Int_t fMultiplicity;/** Multiplicity of particles in the event */
    Int_t fChargedMultiplicity;/** Multiplicity of charged particles in the event */
    Int_t fNeutralMultiplicity; /** Multiplicity of neutral particles in the event */
        
    void CreateEtaPtHisto2D(TString name, TString title);
    void CreateResolutionPtHisto2D(TString name, TString title, TString xtitle, TString ytitle);
    void CreatePtHisto1D(TString name, TString title, TString xtitle, TString ytitle);
    void CreateEtaHisto1D(TString name, TString title);
    void CreateHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh,Int_t ybins,Float_t ylow,Float_t yhigh);
    void CreateHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh);
    void CreateIntHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh);
    void CreateIntHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh,Int_t ybins,Int_t ylow,Int_t yhigh);
    void FillHisto1D(TString histname, Float_t x, Float_t weight);
    void FillHisto2D(TString histname, Float_t x, Float_t y, Float_t weight);
    Bool_t GoodEvent() const {return fGoodEvent;}
    Float_t TrueP(float pTrec) const;

    Float_t Et(TParticle *part, float mass = -1000);
    Float_t Et(Float_t p, Float_t theta, Int_t pid, Short_t charge) const;

    TList *fhistoList;//list of histograms saved out to file
    //static Float_t fgEtaAxis[47];//bins for eta axis of histograms
    static Float_t fgEtaAxis[17];//bins for eta axis of histograms
    static Int_t fgnumOfEtaBins;//number of eta bins
    static Float_t fgPtAxis[117];//bins for pt axis of histograms
    static Int_t fgNumOfPtBins;//number of pt bins
    static Float_t fgResAxis[81];//axis for resolution histograms
    static Int_t fgNumOfResBins;//number of bins for resolution axis
    

    Bool_t fGoodEvent;//boolean to keep track of whether or not this is a good event.

 private:
    //Declare it private to avoid compilation warning
    AliAnalysisHadEt & operator = (const AliAnalysisHadEt & g) ;//cpy assignment
    AliAnalysisHadEt(const AliAnalysisHadEt & g) ; // cpy ctor

    ClassDef(AliAnalysisHadEt, 1);
};

#endif // ALIANALYSISHADET_H
