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
    //TDatabasePDG *fPdgDB;//data base used for looking up pdg codes
    //these codes are stored as variables because otherwise there were issues using this with the plugin
    static Int_t fgPiPlusCode;//pdg pi plus code
    static Int_t fgPiMinusCode;//pdg pi minus code
    static Int_t fgKPlusCode;// pdg k plus code
    static Int_t fgKMinusCode;//pdg k minus code
    static Int_t fgProtonCode;//pdg proton code
    static Int_t fgAntiProtonCode;//pdg antiproton code
    static Int_t fgLambdaCode;// pdg lambda code
    static Int_t fgAntiLambdaCode;//pdg antilambda code
    static Int_t fgK0SCode;//pdg k0 short code
    static Int_t fgOmegaCode;//pdg omega code
    static Int_t fgAntiOmegaCode;//pdg anti-omega code
    static Int_t fgXi0Code;//pdg xi-0 code
    static Int_t fgAntiXi0Code;//pdg anti-xi0 code
    static Int_t fgXiCode;//pdg xi code
    static Int_t fgAntiXiCode;//pdg anti-xi code
    static Int_t fgSigmaCode;//pdg sigma code
    static Int_t fgAntiSigmaCode;//pdg anti-sigma code
    static Int_t fgK0LCode;//pdg k0 long code
    static Int_t fgNeutronCode;//pdg neutron code
    static Int_t fgAntiNeutronCode;//pdg anti-neutron code
    static Int_t fgEPlusCode;//pdg positron code
    static Int_t fgEMinusCode;//pdg electron code
    static Int_t fgGammaCode;//pdg gamma code
    static Int_t fgPi0Code;//pdg neutral pion code
    static Int_t fgEtaCode;//pdg eta code
    static Int_t fgOmega0Code;//pdg eta code
    static Float_t fgPionMass;//pdg pion mass
    static Float_t fgKaonMass;//pdg kaon mass
    static Float_t fgProtonMass;//pdg proton mass
    static Float_t fgElectronMass;//pdg electron mass

    
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
    Float_t Et(Float_t p, Float_t theta, Int_t pid, Short_t charge) const;
    AliESDtrackCuts* fEsdtrackCutsITSTPC;//esd track cuts for ITS+TPC tracks
    AliESDtrackCuts* fEsdtrackCutsTPC;//esd track cuts for TPC tracks (which may also contain ITS hits)
    AliESDtrackCuts* fEsdtrackCutsITS;//esd track cuts for ITS stand alone tracks

    TList *fhistoList;//list of histograms saved out to file
    //static Float_t fgEtaAxis[47];//bins for eta axis of histograms
    static Float_t fgEtaAxis[17];//bins for eta axis of histograms
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
