#ifndef ALIANAPI0_H
#define ALIANAPI0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Class to fill two-photon invariant mass histograms 
// to be used to extract pi0 raw yield.
// Input is produced by AliAnaPhoton (or any other analysis producing output AliAODPWG4Particles), 
// it will do nothing if executed alone
//
//-- Author: Dmitri Peressounko (RRC "KI")
//-- Adapted to PartCorr frame by Lamia Benhabib (SUBATECH)
//-- and Gustavo Conesa (INFN-Frascati)

//Root
class TList;
class TH3D ;
class TH2D ;
class TObjString;

//Analysis
#include "AliAnaPartCorrBaseClass.h"
class AliAODEvent ;
class AliESDEvent ;
class AliAODPWG4Particle ;

class AliAnaPi0 : public AliAnaPartCorrBaseClass {
  
 public:   
  AliAnaPi0() ; // default ctor
  virtual ~AliAnaPi0() ;//virtual dtor
 private:
  AliAnaPi0(const AliAnaPi0 & g) ; // cpy ctor
  AliAnaPi0 & operator = (const AliAnaPi0 & api0) ;//cpy assignment
  
 public:

  //-------------------------------
  // General analysis frame methods
  //-------------------------------

  TObjString * GetAnalysisCuts();
  TList      * GetCreateOutputObjects(); 
  void Terminate(TList* outputList);
  void ReadHistograms(TList * outputList); //Fill histograms with histograms in ouput list, needed in Terminate.
  void Print(const Option_t * opt) const;  
  //void MakeAnalysisFillAOD() {;} //Not needed
  void MakeAnalysisFillHistograms();
  //void Init();
  void InitParameters();

  //Calorimeter options
  TString GetCalorimeter()        const { return fCalorimeter; }
  void SetCalorimeter(TString & det)    { fCalorimeter = det ; }
  void SetNumberOfModules(Int_t nmod)   { fNModules    = nmod; }
  
  //-------------------------------
  // EVENT Bin Methods
  //-------------------------------

  virtual Int_t GetEventIndex(AliAODPWG4Particle * part, Double_t * vert)  ;

  //Setters for parameters of event buffers
  void SetNCentrBin(Int_t n=5)    {fNCentrBin=n ;} //number of bins in centrality 
//void SetNZvertBin(Int_t n=5)    {fNZvertBin=n ;} //number of bins for vertex position
//void SetNRPBin(Int_t n=6)       {fNrpBin=n ;}    //number of bins in reaction plain
  void SetNMaxEvMix(Int_t n=20)   {fNmaxMixEv=n ;} //Maximal number of events for mixing
  
  //Switchs for event multiplicity bin option, by default, centrality
  void SwitchOnTrackMultBins()    {fUseTrackMultBins = kTRUE  ; }
  void SwitchOffTrackMultBins()   {fUseTrackMultBins = kFALSE ; }
  
  void SwitchOnPhotonMultBins()   {fUsePhotonMultBins = kTRUE  ; }
  void SwitchOffPhotonMultBins()  {fUsePhotonMultBins = kFALSE ; }
  
  void SwitchOnClusterEBins()     {fUseAverClusterEBins = kTRUE  ; }
  void SwitchOffClusterEBins()    {fUseAverClusterEBins = kFALSE ; }
  
  void SwitchOnCellEBins()        {fUseAverCellEBins = kTRUE  ; }
  void SwitchOffCellEBins()       {fUseAverCellEBins = kFALSE ; }

  void SwitchOnClusterEDenBins()     {fUseAverClusterEDenBins = kTRUE  ; }
  void SwitchOffClusterEDenBins()    {fUseAverClusterEDenBins = kFALSE ; }
  
//  void SwitchOnClusterPairRBins()     {fUseAverClusterPairRBins = kTRUE  ; }
//  void SwitchOffClusterPairRBins()    {fUseAverClusterPairRBins = kFALSE ; }
//  
//  void SwitchOnClusterPairRWeightBins() {fUseAverClusterPairRWeightBins = kTRUE  ; }
//  void SwitchOffClusterPairRWeightBins(){fUseAverClusterPairRWeightBins = kFALSE ; }

//  void SwitchOnEMaxBins()             {fUseEMaxBins = kTRUE  ; }
//  void SwitchOffEMaxBins()            {fUseEMaxBins = kFALSE ; }

  //-------------------------------
	//Opening angle pair selection
  //-------------------------------
  void SwitchOnAngleSelection()      {fUseAngleCut = kTRUE      ; }
  void SwitchOffAngleSelection()     {fUseAngleCut = kFALSE     ; }
  void SwitchOnAngleEDepSelection()  {fUseAngleEDepCut = kTRUE  ; }
  void SwitchOffAngleEDepSelection() {fUseAngleEDepCut = kFALSE ; }
  void SetAngleCut(Float_t a)        {fAngleCut    = a          ; }
  void SetAngleMaxCut(Float_t a)     {fAngleMaxCut = a          ; }

  //-------------------------------
  // Use mixing code of this class
  //-------------------------------
  void SwitchOnOwnMix()              {fDoOwnMix = kTRUE  ; }
  void SwitchOffOwnMix()             {fDoOwnMix = kFALSE ; }

  //------------------------------------------
  //Do analysis only with clusters in same SM
  //------------------------------------------
  void SwitchOnSameSM()              {fSameSM = kTRUE    ; }
  void SwitchOffSameSM()             {fSameSM = kFALSE   ; }
  
  //-------------------------------
  //Histogram filling options off by default
  //-------------------------------
  void SwitchOnInvPtWeight()         {fMakeInvPtPlots = kTRUE  ; }
  void SwitchOffInvPtWeight()        {fMakeInvPtPlots = kFALSE ; }
  
  void SwitchOnFillBadDistHisto()    {fFillBadDistHisto    = kTRUE;}
  void SwitchOffFillBadDistHisto()   {fFillBadDistHisto    = kFALSE;}
  
  //-------------------------------------------
  //Cuts for multiple analysis, off by default
  //-------------------------------------------
  void SwitchOnMultipleCutAnalysis()          {fMultiCutAna    = kTRUE ;}
  void SwitchOffMultipleCutAnalysis()         {fMultiCutAna    = kFALSE;}

  void SetNPtCuts   (Int_t size)              {if(size <= 10)fNPtCuts    = size; }
  void SetNAsymCuts (Int_t size)              {if(size <= 10)fNAsymCuts  = size; }
  void SetNNCellCuts(Int_t size)              {if(size <= 10)fNCellNCuts = size; }
  void SetNPIDBits  (Int_t size)              {if(size <= 10)fNPIDBits   = size; }
  
  void SetPtCutsAt   (Int_t pos,Float_t val)  {if(pos < 10)fPtCuts[pos]    = val;}
  void SetAsymCutsAt (Int_t pos,Float_t val)  {if(pos < 10)fAsymCuts[pos]  = val;}
  void SetNCellCutsAt(Int_t pos,Int_t val)    {if(pos < 10)fCellNCuts[pos] = val;}
  void SetPIDBitsAt  (Int_t pos,Int_t val)    {if(pos < 10)fPIDBits[pos]   = val;}
  
  //MC analysis related methods
  void FillAcceptanceHistograms();
  void FillMCVersusRecDataHistograms(const Int_t    index1,  const Int_t    index2,
                                     const Float_t  pt1,     const Float_t  pt2, 
                                     const Int_t    ncells1, const Int_t    ncells2, 
                                     const Double_t mass,    const Double_t pt,  const Double_t asym,    
                                     const Double_t deta,    const Double_t dphi);
  
  void SwitchOnMultipleCutAnalysisInSimulation()   {fMultiCutAnaSim = kTRUE;}
  void SwitchOffMultipleCutAnalysisInSimulation()  {fMultiCutAnaSim = kFALSE;}

  
  private:
  Bool_t   fDoOwnMix;            // Do combinatorial background not the one provided by the frame
  Int_t    fNCentrBin ;	         // Number of bins in event container for centrality
//Int_t    fNZvertBin ;	         // Number of bins in event container for vertex position
//Int_t    fNrpBin ;	           // Number of bins in event container for reaction plain
  Int_t    fNmaxMixEv ;	         // Maximal number of events stored in buffer for mixing
  TString  fCalorimeter ;        // Select Calorimeter for IM
  Int_t    fNModules ;           // Number of EMCAL/PHOS modules, set as many histogras as modules 
  Bool_t   fUseAngleCut ;        // Select pairs depending on their opening angle
  Bool_t   fUseAngleEDepCut ;    // Select pairs depending on their opening angle
  Float_t  fAngleCut ;           // Select pairs with opening angle larger than a threshold
  Float_t  fAngleMaxCut ;        // Select pairs with opening angle smaller than a threshold
  TList ** fEventsList ;         //![fNCentrBin*GetNZvertBin()*GetNRPBin()] Containers for photons in stored events
  
  //Multiple cuts analysis
  Bool_t   fMultiCutAna;         // Do analysis with several or fixed cut
  Bool_t   fMultiCutAnaSim;      // Do analysis with several or fixed cut, in the simulation related part
  Int_t    fNPtCuts;             // Number of pt cuts
  Float_t  fPtCuts[10];          // Array with different pt cuts
  Int_t    fNAsymCuts;           // Number of assymmetry cuts
  Float_t  fAsymCuts[10];        // Array with different assymetry cuts
  Int_t    fNCellNCuts;          // Number of cuts with number of cells in cluster
  Int_t    fCellNCuts[10];       // Array with different cell number cluster cuts
  Int_t    fNPIDBits ;		       // Number of possible PID bit combinations
  Int_t    fPIDBits[10];         // Array with different PID bits
  
  //Switchs of different analysis options
  Bool_t   fMakeInvPtPlots;      // D plots with inverse pt weight
  Bool_t   fSameSM;              // Select only pairs in same SM;
  Bool_t   fUseTrackMultBins;    // Use track multiplicity and not centrality bins
  Bool_t   fUsePhotonMultBins;   // Use photon multiplicity and not centrality bins
  Bool_t   fUseAverCellEBins;    // Use cell average energy and not centrality bins
  Bool_t   fUseAverClusterEBins; // Use cluster average energy and not centrality bins
  Bool_t   fUseAverClusterEDenBins; // Use cluster average energy density and not centrality bins
//  Bool_t   fUseAverClusterPairRBins; // Use cluster average energy and not centrality bins
//  Bool_t   fUseAverClusterPairRWeightBins; // Use cluster average energy and not centrality bins
//  Bool_t   fUseEMaxBins;         // Use Emax bins
  Bool_t   fFillBadDistHisto;    // Do plots for different distances to bad channels
  
  //Histograms
  
  //Event characterization
  TH1F* fhAverTotECluster;          //! Average number of clusters in SM
  TH1F* fhAverTotECell;             //! Average number of cells    in SM
  TH2F* fhAverTotECellvsCluster;    //! Average number of cells    in SM
  TH1F* fhEDensityCluster;          //! Deposited energy in event per cluster
  TH1F* fhEDensityCell;             //! Deposited energy in event per cell vs cluster
  TH2F* fhEDensityCellvsCluster;    //! Deposited energy in event per cell vs cluster
//  TH1F* fhClusterPairDist;          //! Distance between clusters
//  TH1F* fhClusterPairDistWeight;    //! Distance between clusters weighted by pair energy
//  TH1F* fhAverClusterPairDist;      //! Average distance between cluster pairs
//  TH1F* fhAverClusterPairDistWeight;//! Average distance between cluster pairs weighted with pair energy
//  TH2F* fhAverClusterPairDistvsAverE;        //! Average distance between cluster pairs vs average cluster energy
//  TH2F* fhAverClusterPairDistWeightvsAverE;  //! Average distance between cluster pairs weighted with pair energy vs average cluster energy
//  TH2F* fhAverClusterPairDistvsN;            //! Average distance between cluster pairs vs number of clusters
//  TH2F* fhAverClusterPairDistWeightvsN;      //! Average distance between cluster pairs weighted with pair energy vs number of clusters
//  TH2F* fhMaxEvsClustMult;          //!
//  TH2F* fhMaxEvsClustEDen;          //!

  
  TH2D ** fhReMod ;                 //![fNModules]   REAL  two-photon invariant mass distribution for different calorimeter modules.
  TH2D ** fhReDiffMod ;             //![fNModules+1] REAL  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2D ** fhMiMod ;                 //![fNModules]   MIXED two-photon invariant mass distribution for different calorimeter modules.
  TH2D ** fhMiDiffMod ;             //![fNModules+1] MIXED two-photon invariant mass distribution for different clusters in different calorimeter modules.
  
  // Pairs with at least one cluster tagged as conversion
  TH2D * fhReConv ;                 //! REAL  two-photon invariant mass distribution one of the pair was 2 clusters with small mass 
  TH2D * fhMiConv ;                 //! MIXED two-photon invariant mass distribution one of the pair was 2 clusters with small mass
  TH2D * fhReConv2 ;                //! REAL  two-photon invariant mass distribution both pair photons recombined from 2 clusters with small mass 
  TH2D * fhMiConv2 ;                //! MIXED two-photon invariant mass distribution both pair photons recombined from 2 clusters with small mass

  TH2D ** fhRe1 ;                   //![fNCentrBin*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry 
  TH2D ** fhMi1 ;                   //![fNCentrBin*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry
  TH2D ** fhRe2 ;                   //![fNCentrBin*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry 
  TH2D ** fhMi2 ;                   //![fNCentrBin*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry
  TH2D ** fhRe3 ;                   //![fNCentrBin*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry 
  TH2D ** fhMi3 ;                   //![fNCentrBin*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry
  
  //Histograms weighted by inverse pT
  TH2D ** fhReInvPt1 ;              //![fNCentrBin*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  TH2D ** fhMiInvPt1 ;              //![fNCentrBin*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  TH2D ** fhReInvPt2 ;              //![fNCentrBin*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT 
  TH2D ** fhMiInvPt2 ;              //![fNCentrBin*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  TH2D ** fhReInvPt3 ;              //![fNCentrBin*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  TH2D ** fhMiInvPt3 ;              //![fNCentrBin*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  
  //Multiple cuts: Assymmetry, pt, n cells, PID
  TH2D ** fhRePtNCellAsymCuts ;     //![fNPtCuts*fNAsymCuts*fNCellNCuts] REAL two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry
  TH2D ** fhRePtNCellAsymCutsSM0 ;  //![fNPtCuts*fNAsymCuts*fNCellNCuts] REAL two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry
  TH2D ** fhRePtNCellAsymCutsSM1 ;  //![fNPtCuts*fNAsymCuts*fNCellNCuts] REAL two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry
  TH2D ** fhRePtNCellAsymCutsSM2 ;  //![fNPtCuts*fNAsymCuts*fNCellNCuts] REAL two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry
  TH2D ** fhRePtNCellAsymCutsSM3 ;  //![fNPtCuts*fNAsymCuts*fNCellNCuts] REAL two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry
  TH2D ** fhMiPtNCellAsymCuts ;     //![fNPtCuts*fNAsymCuts*fNCellNCuts] Mixed two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry
  TH2D ** fhRePIDBits ;             //![fNPIDBits]  REAL two-photon invariant mass distribution for different PID bits
  TH3D ** fhRePtMult ;              //![fNAsymCuts] REAL two-photon invariant mass distribution for different track multiplicity and assymetry cuts
  
  // Asymmetry vs pt, in pi0/eta regions
  TH2D *  fhRePtAsym    ;           //! REAL two-photon pt vs asymmetry
  TH2D *  fhRePtAsymPi0 ;           //! REAL two-photon pt vs asymmetry, close to pi0 mass
  TH2D *  fhRePtAsymEta ;           //! REAL two-photon pt vs asymmetry, close to eta mass

  TH3D * fhEvents;                  //! Number of events per centrality, RP, zbin
  
  // Pair opening angle
  TH2D * fhRealOpeningAngle ;       //! Opening angle of pair versus pair energy
  TH2D * fhRealCosOpeningAngle ;    //! Cosinus of opening angle of pair version pair energy
  TH2D * fhMixedOpeningAngle ;      //! Opening angle of pair versus pair energy
  TH2D * fhMixedCosOpeningAngle ;   //! Cosinus of opening angle of pair version pair energy
  
  //MC analysis histograms
  //Pi0 Acceptance
  TH1D * fhPrimPi0Pt ;              //! Spectrum of Primary 
  TH1D * fhPrimPi0AccPt ;           //! Spectrum of primary with accepted daughters 
  TH2D * fhPrimPi0Y ;               //! Rapidity distribution of primary particles  vs pT
  TH2D * fhPrimPi0AccY ;            //! Rapidity distribution of primary with accepted daughters  vs pT
  TH2D * fhPrimPi0Phi ;             //! Azimutal distribution of primary particles  vs pT
  TH2D * fhPrimPi0AccPhi;           //! Azimutal distribution of primary with accepted daughters  vs pT
  TH2D * fhPrimPi0OpeningAngle ;    //! Opening angle of pair versus pair energy, primaries
  TH2D * fhPrimPi0CosOpeningAngle ; //! Cosinus of opening angle of pair version pair energy, primaries
  //Eta acceptance
  TH1D * fhPrimEtaPt ;              //! Spectrum of Primary 
  TH1D * fhPrimEtaAccPt ;           //! Spectrum of primary with accepted daughters 
  TH2D * fhPrimEtaY ;               //! Rapidity distribution of primary particles vs pT
  TH2D * fhPrimEtaAccY ;            //! Rapidity distribution of primary with accepted daughters  vs pT
  TH2D * fhPrimEtaPhi ;             //! Azimutal distribution of primary particles  vs pT
  TH2D * fhPrimEtaAccPhi;           //! Azimutal distribution of primary with accepted daughters	 vs pT
  
  // Primaries origin
  TH2D * fhPrimPi0PtOrigin ;        //! Spectrum of generated pi0 vs mother
  TH2D * fhPrimEtaPtOrigin ;        //! Spectrum of generated eta vs mother
  
  //Pair origin
  //Array of histograms ordered as follows: 0-Photon, 1-electron, 2-pi0, 3-eta, 4-a-proton, 5-a-neutron, 6-stable particles, 
  // 7-other decays, 8-string, 9-final parton, 10-initial parton, intermediate, 11-colliding proton, 12-unrelated
  TH2D * fhMCOrgMass[13];           //! Mass vs pt of real pairs, check common origin of pair
  TH2D * fhMCOrgAsym[13];           //! Asymmetry vs pt of real pairs, check common origin of pair
  TH2D * fhMCOrgDeltaEta[13];       //! Delta Eta vs pt of real pairs, check common origin of pair
  TH2D * fhMCOrgDeltaPhi[13];       //! Delta Phi vs pt of real pairs, check common origin of pair
  
  //Multiple cuts in simulation, origin pi0 or eta
  TH2D ** fhMCPi0MassPtRec;         //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real pi0 pairs, reconstructed mass vs reconstructed pt of original pair  
  TH2D ** fhMCPi0MassPtTrue;        //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real pi0 pairs, reconstructed mass vs generated pt of original pair  
  TH2D ** fhMCPi0PtTruePtRec;       //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real pi0 pairs, reconstructed pt vs generated pt of pair
  TH2D ** fhMCEtaMassPtRec;         //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real eta pairs, reconstructed mass vs reconstructed pt of original pair  
  TH2D ** fhMCEtaMassPtTrue;        //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real eta pairs, reconstructed mass vs generated pt of original pair  
  TH2D ** fhMCEtaPtTruePtRec;       //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real eta pairs, reconstructed pt vs generated pt of pair
  
  TH2D *  fhMCPi0PtOrigin ;         //! Mass of reoconstructed pi0 pairs  in calorimeter vs mother
  TH2D *  fhMCEtaPtOrigin ;         //! Mass of reoconstructed pi0 pairs  in calorimeter vs mother

  ClassDef(AliAnaPi0,15)
} ;


#endif //ALIANAPI0_H



