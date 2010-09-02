//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
#ifndef ALIANALYSISHADET_H
#define ALIANALYSISHADET_H

#include "TString.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "Rtypes.h"
#include "TString.h"
#include "AliESDtrackCuts.h"

class TH2F;
class TH1F;
class AliVEvent;
class TList;

class AliAnalysisHadEt
{
public:
   
  AliAnalysisHadEt();
    virtual ~AliAnalysisHadEt();

    /** Analyse the event! */
    virtual Int_t AnalyseEvent(AliVEvent *event);

    /** Fill the objects you want to output, classes which add new histograms should overload this. */
    virtual void FillOutputList();

    /** Initialise the analysis, must be overloaded. */
    virtual void Init();

    /** 
    * Creates the histograms, must be overloaded if you want to add your own. 
    * Uses the fHistogramNameSuffix to create proper histogram names
    */
    virtual void CreateHistograms();
    
    /** Fills the histograms, must be overloaded if you want to add your own */
    virtual void FillHistograms();

    /** Reset event specific values (Et etc.) */
    virtual void ResetEventValues();
    
    /** Sum of the total Et for all events */
    Double_t GetSumEt() { return fSumEt; }

    /** Sum of the total Et within our acceptance for all events */
    Double_t GetSumEtAcc() { return fSumEtAcc; }

    /** Total Et in the event (without acceptance cuts) */
    Double_t GetTotEt() { return fTotEt; }

    /** Total Et in the event within the acceptance cuts */
    Double_t GetTotEtAcc() { return fTotEtAcc; }

   /** Total neutral Et in the event (without acceptance cuts) */
    Double_t GetTotNeutralEt() { return fTotNeutralEt; }

    /** Total neutral Et in the event within the acceptance cuts */
    Double_t GetTotNeutralEtAcc() { return fTotNeutralEtAcc; }
    
    /** Total charged Et in the event (without acceptance cuts) */
    Double_t GetTotChargedEt() { return fTotChargedEt; }

    /** Total charged Et in the event within the acceptance cuts */
    Double_t GetTotChargedEtAcc() { return fTotChargedEtAcc; }


    void SetHistoList(TList *mylist){histoList = mylist;}

    void SetTPCITSTrackCuts(AliESDtrackCuts *cuts){ esdtrackCutsITSTPC = cuts;}
    void SetTPCOnlyTrackCuts(AliESDtrackCuts *cuts){ esdtrackCutsTPC = cuts;}
    void SetITSTrackCuts(AliESDtrackCuts *cuts){ esdtrackCutsITS = cuts;}

protected:
   
    /** The suffix for the histogram names */
    TString fHistogramNameSuffix;

    /** PDG Database */
    TDatabasePDG *fPdgDB;
    Int_t PiPlusCode;
    Int_t PiMinusCode;
    Int_t KPlusCode;
    Int_t KMinusCode;
    Int_t ProtonCode;
    Int_t AntiProtonCode;
    Int_t LambdaCode;
    Int_t AntiLambdaCode;
    Int_t K0SCode;
    Int_t OmegaCode;
    Int_t AntiOmegaCode;
    Int_t Xi0Code;
    Int_t AntiXi0Code;
    Int_t XiCode;
    Int_t AntiXiCode;
    Int_t SigmaCode;
    Int_t AntiSigmaCode;
    Int_t K0LCode;
    Int_t NeutronCode;
    Int_t AntiNeutronCode;
    Int_t EPlusCode;
    Int_t EMinusCode;
    Float_t PionMass;

    /** Sum of the total Et for all events */
    Double_t fSumEt;

    /** Sum of the total Et within our acceptance for all events */
    Double_t fSumEtAcc;

    /** Total Et in the event (without acceptance cuts) */
    Double_t fTotEt;

    /** Total Et in the event within the acceptance cuts */
    Double_t fTotEtAcc;

    /** Total neutral Et in the event */
    Double_t fTotNeutralEt;

    /** Total neutral Et in the event within the acceptance cuts */
    Double_t fTotNeutralEtAcc;

    /** Total charged Et in the event */
    Double_t fTotChargedEt;

    /** Total charged Et in the event within the acceptance cuts */
    Double_t fTotChargedEtAcc;

    /** Multiplicity of particles in the event */
    Int_t fMultiplicity;
    
    /** Multiplicity of charged particles in the event */
    Int_t fChargedMultiplicity;
    
    /** Multiplicity of neutral particles in the event */
    Int_t fNeutralMultiplicity; 
    
    /** Cut in eta ( normally |eta| < 0.5 */
    Double_t fEtaCut;

    /** Eta cut for our acceptance */
    Double_t fEtaCutAcc;

    /** Vertex cuts */
    Double_t fVertexXCut;
    Double_t fVertexYCut;
    Double_t fVertexZCut;

    /** Impact parameter cuts */
    Double_t fIPxyCut;
    Double_t fIPzCut;


    void CreateEtaPtHisto2D(TString name, TString title);
    void CreateEtaHisto1D(TString name, TString title);
    void CreateHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh,Int_t ybins,Float_t ylow,Float_t yhigh);
    void CreateHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh);
    void CreateIntHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh);
    void CreateIntHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh,Int_t ybins,Int_t ylow,Int_t yhigh);
    void FillHisto1D(TString histname, Float_t x, Float_t weight);
    void FillHisto2D(TString histname, Float_t x, Float_t y, Float_t weight);

    Float_t Et(TParticle *part, float mass = -1000);
    AliESDtrackCuts* esdtrackCutsITSTPC;
    AliESDtrackCuts* esdtrackCutsTPC;
    AliESDtrackCuts* esdtrackCutsITS;

    TList *histoList;
    static Float_t etaAxis[47];
    static Int_t numOfEtaBins;
    static Float_t ptAxis[117];
    static Int_t numOfPtBins;
    

 private:

private:
  //Declare it private to avoid compilation warning
  AliAnalysisHadEt & operator = (const AliAnalysisHadEt & g) ;//cpy assignment
  AliAnalysisHadEt(const AliAnalysisHadEt & g) ; // cpy ctor

    ClassDef(AliAnalysisHadEt, 0);
};

#endif // ALIANALYSISHADET_H
