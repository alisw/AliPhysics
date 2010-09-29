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

class TH2F;
class TH1F;
class AliVEvent;
class TList;
class AliESDtrackCuts;
class Rtypes;
class TParticle;
class TDatabasePDG;
class AliAnalysisEtCuts;

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

    /** Set Particle codes/mass */
    virtual void SetParticleCodes();

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

    void SetTPCITSTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsITSTPC = (AliESDtrackCuts *) cuts;}
    void SetTPCOnlyTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsTPC = (AliESDtrackCuts *) cuts;}
    void SetITSTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsITS = (AliESDtrackCuts *) cuts;}

protected:   
    
    TString fHistogramNameSuffix; /** The suffix for the histogram names */

    AliAnalysisEtCuts *fCuts; // keeper of basic cuts

    /** PDG Database */
    TDatabasePDG *fPdgDB;//data base used for looking up pdg codes
    //these codes are stored as variables because otherwise there were issues using this with the plugin
    Int_t fPiPlusCode;//pdg pi plus code
    Int_t fPiMinusCode;//pdg pi minus code
    Int_t fKPlusCode;// pdg k plus code
    Int_t fKMinusCode;//pdg k minus code
    Int_t fProtonCode;//pdg proton code
    Int_t fAntiProtonCode;//pdg antiproton code
    Int_t fLambdaCode;// pdg lambda code
    Int_t fAntiLambdaCode;//pdg antilambda code
    Int_t fK0SCode;//pdg k0 short code
    Int_t fOmegaCode;//pdg omega code
    Int_t fAntiOmegaCode;//pdg anti-omega code
    Int_t fXi0Code;//pdg xi-0 code
    Int_t fAntiXi0Code;//pdg anti-xi0 code
    Int_t fXiCode;//pdg xi code
    Int_t fAntiXiCode;//pdg anti-xi code
    Int_t fSigmaCode;//pdg sigma code
    Int_t fAntiSigmaCode;//pdg anti-sigma code
    Int_t fK0LCode;//pdg k0 long code
    Int_t fNeutronCode;//pdg neutron code
    Int_t fAntiNeutronCode;//pdg anti-neutron code
    Int_t fEPlusCode;//pdg positron code
    Int_t fEMinusCode;//pdg electron code
    Int_t fGammaCode;//pdg gamma code
    Int_t fPi0Code;//pdg neutral pion code
    Int_t fEtaCode;//pdg eta code
    Int_t fOmega0Code;//pdg eta code
    Float_t fPionMass;//pdg pion mass
    Float_t fKaonMass;//pdg kaon mass
    Float_t fProtonMass;//pdg proton mass
    Float_t fElectronMass;//pdg electron mass

    
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
    void CreateEtaHisto1D(TString name, TString title);
    void CreateHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh,Int_t ybins,Float_t ylow,Float_t yhigh);
    void CreateHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh);
    void CreateIntHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh);
    void CreateIntHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh,Int_t ybins,Int_t ylow,Int_t yhigh);
    void FillHisto1D(TString histname, Float_t x, Float_t weight);
    void FillHisto2D(TString histname, Float_t x, Float_t y, Float_t weight);

    Float_t Et(TParticle *part, float mass = -1000);
    Float_t Et(Float_t p, Float_t theta, Int_t pid, Short_t charge);
    AliESDtrackCuts* fEsdtrackCutsITSTPC;//esd track cuts for ITS+TPC tracks
    AliESDtrackCuts* fEsdtrackCutsTPC;//esd track cuts for TPC tracks (which may also contain ITS hits)
    AliESDtrackCuts* fEsdtrackCutsITS;//esd track cuts for ITS stand alone tracks

    TList *fhistoList;//list of histograms saved out to file
    static Float_t fgEtaAxis[47];//bins for eta axis of histograms
    static Int_t fgnumOfEtaBins;//number of eta bins
    static Float_t fgPtAxis[117];//bins for pt axis of histograms
    static Int_t fgNumOfPtBins;//number of pt bins
    static Float_t fgPtTPCCutOff;//cut off for tracks in TPC
    static Float_t fgPtITSCutOff;//cut off for tracks in ITS
    


 private:
    //Declare it private to avoid compilation warning
    AliAnalysisHadEt & operator = (const AliAnalysisHadEt & g) ;//cpy assignment
    AliAnalysisHadEt(const AliAnalysisHadEt & g) ; // cpy ctor

    ClassDef(AliAnalysisHadEt, 1);
};

#endif // ALIANALYSISHADET_H
