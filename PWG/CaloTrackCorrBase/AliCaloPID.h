#ifndef ALICALOPID_H
#define ALICALOPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Class for PID selection with calorimeters
// The Output of the main method GetIdentifiedParticleType is a PDG number identifying the cluster, 
// being kPhoton, kElectron, kPi0 ... as defined in the header file
//   - GetIdentifiedParticleType(const AliVCluster * cluster) 
//      Assignes a PID tag to the cluster, right now there is the possibility to : use bayesian weights from reco, 
//      recalculate them (EMCAL) or use other procedures not used in reco.
//      In order to recalculate Bayesian, it is necessary to load the EMCALUtils library
//      and do SwitchOnBayesianRecalculation().
//      To change the PID parameters from Low to High like the ones by default, use the constructor 
//      AliCaloPID(flux)
//      where flux is AliCaloPID::kLow or AliCaloPID::kHigh
//      If it is necessary to change the parameters use the constructor 
//      AliCaloPID(AliEMCALPIDUtils *utils) and set the parameters before.

//   - GetGetIdentifiedParticleTypeFromBayesian(const Double_t * pid, const Float_t energy)
//      Reads the PID weights array of the ESDs and depending on its magnitude identifies the particle, 
//      executed when bayesian is ON by GetIdentifiedParticleType(const AliVCluster * cluster) 
//   - SetPIDBits: Simple PID, depending on the thresholds fLOCut fTOFCut and even the
//     result of the PID bayesian a different PID bit is set. 
//
//
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TObject.h> 
class TString ;
class TLorentzVector ;
#include <TFormula.h>
class TList;
class TH2F ;

//--- AliRoot system ---
class AliVCluster;
class AliVCaloCells;
class AliAODPWG4Particle;
class AliEMCALPIDUtils;
class AliCalorimeterUtils;
class AliVEvent;

class AliCaloPID : public TObject {
	
 public: 
  
  AliCaloPID() ; // ctor
  AliCaloPID(const Int_t particleFlux) ; // ctor, to be used when recalculating bayesian PID
  AliCaloPID(const TNamed * emcalpid) ; // ctor, to be used when recalculating bayesian PID and need different parameters
  virtual ~AliCaloPID() ;//virtual dtor
  	
  enum PidType 
  {
    kPhoton         = 22,
    kPi0            = 111,
    kEta            = 221, 
    kElectron       = 11, 
    kEleCon         =-11, 
    kNeutralHadron  = 2112, 
    kChargedHadron  = 211, 
    kNeutralUnknown = 130, 
    kChargedUnknown = 321
  };
  
  enum TagType {kPi0Decay, kEtaDecay, kOtherDecay, kConversion, kNoTag = -1};
  
  // Main methods
  
  TList *   GetCreateOutputObjects();

  void      InitParameters();
  
  Bool_t    IsInPi0SplitAsymmetryRange(const Float_t energy, const Float_t asy,  const Int_t nlm);

  Bool_t    IsInPi0SplitMassRange     (const Float_t energy, const Float_t mass, const Int_t nlm);
  
  Bool_t    IsInPi0M02Range           (const Float_t energy, const Float_t m02,  const Int_t nlm);
  Bool_t    IsInEtaM02Range           (const Float_t energy, const Float_t m02,  const Int_t nlm);
  Bool_t    IsInConM02Range           (const Float_t energy, const Float_t m02,  const Int_t nlm);
  
  
  Int_t     GetIdentifiedParticleTypeFromBayesWeights(const Bool_t isEMCAL, const Double_t * pid, const Float_t energy) ;

  Int_t     GetIdentifiedParticleTypeFromClusterSplitting(AliVCluster * cluster, AliVCaloCells* cells, 
                                                          AliCalorimeterUtils * caloutils,
                                                          Double_t vertex[3], 
                                                          Int_t & nLocMax, Double_t & mass, Double_t & angle,
                                                          Double_t & e1  , Double_t & e2                     ) ;
  
  Int_t     GetIdentifiedParticleType(const AliVCluster * cluster) ;
  
  TString   GetPIDParametersList();
  
  Bool_t    IsTrackMatched(AliVCluster * cluster, AliCalorimeterUtils* cu, AliVEvent* event) const ;    
  
  void      SetPIDBits(AliVCluster * cluster, AliAODPWG4Particle *aodph, 
                       AliCalorimeterUtils* cu, AliVEvent* event);
  
  void      Print(const Option_t * opt)const;
  
  void      PrintClusterPIDWeights(const Double_t * pid) const;
  
  //Check if cluster photon-like. Uses photon cluster parameterization in real pp data 
  //Returns distance in sigmas. Recommended cut 2.5
  Float_t TestPHOSDispersion(const Double_t pt, const Double_t m20, const Double_t m02) const ; 
  //Checks distance to the closest track. Takes into account 
  //non-perpendicular incidence of tracks.
  Float_t TestPHOSChargedVeto(const Double_t dx,  const Double_t dz, const Double_t ptTrack, 
                              const Int_t chargeTrack, const Double_t mf) const ;
  
  // Setters, getters
  
  void    SetDebug(Int_t deb)                  { fDebug = deb                 ; }
  Int_t   GetDebug()                     const { return fDebug                ; }	

  enum    eventType{kLow,kHigh};
  void    SetLowParticleFlux()                 { fParticleFlux        = kLow  ; }
  void    SetHighParticleFlux()                { fParticleFlux        = kHigh ; }  
  // not really used, only for bayesian recalculation in EMCAL, but could be useful in future
  
  // Bayesian
  
  void    SwitchOnBayesian()                   { fUseBayesianWeights  = kTRUE ; }
  void    SwitchOffBayesian()                  { fUseBayesianWeights  = kFALSE; }
  void    SwitchOnBayesianRecalculation()      { fRecalculateBayesian = kTRUE ; fUseBayesianWeights  = kTRUE ;} // EMCAL
  void    SwitchOffBayesianRecalculation()     { fRecalculateBayesian = kFALSE; } // EMCAL
  
  AliEMCALPIDUtils * GetEMCALPIDUtils() ; 
  
  //Weight getters
  Float_t GetEMCALPhotonWeight()         const { return fEMCALPhotonWeight    ; }
  Float_t GetEMCALPi0Weight()            const { return fEMCALPi0Weight       ; }
  Float_t GetEMCALElectronWeight()       const { return fEMCALElectronWeight  ; }
  Float_t GetEMCALChargeWeight()         const { return fEMCALChargeWeight    ; }
  Float_t GetEMCALNeutralWeight()        const { return fEMCALNeutralWeight   ; }
  Float_t GetPHOSPhotonWeight()          const { return fPHOSPhotonWeight     ; }
  Float_t GetPHOSPi0Weight()             const { return fPHOSPi0Weight        ; }
  Float_t GetPHOSElectronWeight()        const { return fPHOSElectronWeight   ; }
  Float_t GetPHOSChargeWeight()          const { return fPHOSChargeWeight     ; }
  Float_t GetPHOSNeutralWeight()         const { return fPHOSNeutralWeight    ; }
  
  Bool_t  IsPHOSPIDWeightFormulaOn()     const { return fPHOSWeightFormula    ; } 

  TFormula * GetPHOSPhotonWeightFormula()      { 
    if(!fPHOSPhotonWeightFormula) 
      fPHOSPhotonWeightFormula = new TFormula("phos_photon_weight",
                                              fPHOSPhotonWeightFormulaExpression);
    return fPHOSPhotonWeightFormula                                           ; } 
  
  TFormula * GetPHOSPi0WeightFormula()         { 
    if(!fPHOSPi0WeightFormula) 
      fPHOSPi0WeightFormula = new TFormula("phos_pi0_weight",
                                           fPHOSPi0WeightFormulaExpression);
    return fPHOSPi0WeightFormula                                              ; } 
  
  TString  GetPHOSPhotonWeightFormulaExpression() const { return fPHOSPhotonWeightFormulaExpression ; } 
  TString  GetPHOSPi0WeightFormulaExpression()    const { return fPHOSPi0WeightFormulaExpression    ; } 
  
  //Weight setters
  void    SetEMCALPhotonWeight  (Float_t  w)   { fEMCALPhotonWeight   = w     ; }
  void    SetEMCALPi0Weight     (Float_t  w)   { fEMCALPi0Weight      = w     ; }
  void    SetEMCALElectronWeight(Float_t  w)   { fEMCALElectronWeight = w     ; }
  void    SetEMCALChargeWeight  (Float_t  w)   { fEMCALChargeWeight   = w     ; }
  void    SetEMCALNeutralWeight (Float_t  w)   { fEMCALNeutralWeight  = w     ; }
  void    SetPHOSPhotonWeight   (Float_t  w)   { fPHOSPhotonWeight    = w     ; }
  void    SetPHOSPi0Weight      (Float_t  w)   { fPHOSPi0Weight       = w     ; }
  void    SetPHOSElectronWeight (Float_t  w)   { fPHOSElectronWeight  = w     ; }
  void    SetPHOSChargeWeight   (Float_t  w)   { fPHOSChargeWeight    = w     ; }
  void    SetPHOSNeutralWeight  (Float_t  w)   { fPHOSNeutralWeight   = w     ; }
  
  void    UsePHOSPIDWeightFormula   (Bool_t ok )           { fPHOSWeightFormula                 = ok ; } 
  void    SetPHOSPhotonWeightFormulaExpression(TString ph) { fPHOSPhotonWeightFormulaExpression = ph ; } 
  void    SetPHOSPi0WeightFormulaExpression   (TString pi) { fPHOSPi0WeightFormulaExpression    = pi ; }
  
  //PID cuts 
  
  void    SetEMCALLambda0CutMax(Float_t lcut ) { fEMCALL0CutMax     = lcut    ; }
  Float_t GetEMCALLambda0CutMax()        const { return fEMCALL0CutMax        ; }   
  
  void    SetEMCALLambda0CutMin(Float_t lcut ) { fEMCALL0CutMin     = lcut    ; }
  Float_t GetEMCALLambda0CutMin()        const { return fEMCALL0CutMin        ; }   
  
  void    SetEMCALDEtaCut(Float_t dcut )       { fEMCALDEtaCut      = dcut    ; }
  Float_t GetEMCALDEtaCut()              const { return fEMCALDEtaCut         ; }   
  
  void    SetEMCALDPhiCut(Float_t dcut )       { fEMCALDPhiCut      = dcut    ; }
  Float_t GetEMCALDPhiCut()              const { return fEMCALDPhiCut         ; }   
  
  void    SetTOFCut(Float_t tcut )             { fTOFCut            = tcut    ; }
  Float_t GetTOFCut()                    const { return fTOFCut               ; }   
  
  void    SetPHOSRCut(Float_t rcut )           { fPHOSRCut          = rcut    ; }
  Float_t GetPHOSRCut()                  const { return fPHOSRCut             ; }   

  void    SetPHOSDispersionCut(Float_t dcut )  { fPHOSDispersionCut = dcut    ; }
  Float_t GetPHOSDispersionCut()         const { return fPHOSDispersionCut    ; }   
  
  // Cluster splitting analysis
  
  void    SwitchOnClusterSplittingPID()        { fDoClusterSplitting = kTRUE  ; }
  void    SwitchOffClusterplittingPID()        { fDoClusterSplitting = kFALSE ; }
  
  void    SwitchOnSimpleSplitMassCut()         { fUseSimpleMassCut   = kTRUE  ; }
  void    SwitchOffSimpleSplitMassCut()        { fUseSimpleMassCut   = kFALSE ; }

  void    SwitchOnSimpleSplitM02Cut()          { fUseSimpleM02Cut    = kTRUE  ; }
  void    SwitchOffSimpleSplitM02Cut()         { fUseSimpleM02Cut    = kFALSE ; }
  
  void    SwitchOnSplitAsymmetryCut()          { fUseSplitAsyCut     = kTRUE  ; }
  void    SwitchOffSplitAsymmetryCut()         { fUseSplitAsyCut     = kFALSE ; }
  
  void    SetClusterSplittingM02Cut(Float_t min=0, Float_t max=100) 
  { fSplitM02MinCut   = min ; fSplitM02MaxCut  = max ; }
  
  void    SetClusterSplittingMinNCells(Int_t cut)   { fSplitMinNCells = cut   ; }
  
  void    SetSplitEnergyFractionMinimum(Int_t i, Float_t min){ if (i < 3 && i >=0 ) fSplitEFracMin[i]  = min   ; }
  Float_t GetSplitEnergyFractionMinimum(Int_t i) const       { if( i < 3 && i >=0 ) return fSplitEFracMin[i]   ;  else return 0 ; }
  
  Float_t GetPi0MinMass()                const { return fMassPi0Min           ; } // Simple cut case
  Float_t GetEtaMinMass()                const { return fMassEtaMin           ; } // Simple cut case
  Float_t GetPhotonMinMass()             const { return fMassPhoMin           ; }  
  Float_t GetPi0MaxMass()                const { return fMassPi0Max           ; }
  Float_t GetEtaMaxMass()                const { return fMassEtaMax           ; }
  Float_t GetPhotonMaxMass()             const { return fMassPhoMax           ; }
  
  void    SetSplitWidthSigma(Float_t s)        { fSplitWidthSigma        = s  ; }
  void    SetPi0MassWidthSelectionParameters    (Int_t iparam, Float_t param) { if(iparam < 7 ) fMassWidthPi0Param[iparam] = param ; }
  void    SetM02MaximumSelectionParameters      (Int_t inlm, Int_t iparam, Float_t param)
  { if(iparam < 6 && inlm < 2) fM02MaxParam[inlm][iparam] = param ; }
  void    SetM02MinimumSelectionParameters      (Int_t inlm, Int_t iparam, Float_t param)
  { if(iparam < 6 && inlm < 2) fM02MinParam[inlm][iparam] = param ; }
  void    SetAsymmetryMinimumSelectionParameters(Int_t inlm, Int_t iparam, Float_t param)
  { if(iparam < 6 && inlm < 2) fAsyMinParam[inlm][iparam] = param ; }

  void    SetPi0MassRange(Float_t min, Float_t max)    { fMassPi0Min    = min ; fMassPi0Max = max ; } // Simple case
  void    SetEtaMassRange(Float_t min, Float_t max)    { fMassEtaMin    = min ; fMassEtaMax = max ; }
  void    SetPhotonMassRange(Float_t min, Float_t max) { fMassPhoMin    = min ; fMassPhoMax = max ; }
    
private:
  
  Int_t	    fDebug;                             // Debug level
  Int_t     fParticleFlux;                      // Particle flux for setting PID parameters

  // Bayesian
  AliEMCALPIDUtils * fEMCALPIDUtils;            // Pointer to EMCALPID to redo the PID Bayesian calculation
  Bool_t    fUseBayesianWeights;                // Select clusters based on weights calculated in reconstruction
  Bool_t    fRecalculateBayesian;               // Recalculate PID bayesian or use simple PID?

  Float_t   fEMCALPhotonWeight;                 // Bayesian PID weight for photons in EMCAL 
  Float_t   fEMCALPi0Weight;                    // Bayesian PID weight for pi0 in EMCAL 
  Float_t   fEMCALElectronWeight;               // Bayesian PID weight for electrons in EMCAL 
  Float_t   fEMCALChargeWeight;                 // Bayesian PID weight for charged hadrons in EMCAL 
  Float_t   fEMCALNeutralWeight;                // Bayesian PID weight for neutral hadrons in EMCAL 
  Float_t   fPHOSPhotonWeight;                  // Bayesian PID weight for photons in PHOS 
  Float_t   fPHOSPi0Weight;                     // Bayesian PID weight for pi0 in PHOS 
  Float_t   fPHOSElectronWeight;                // Bayesian PID weight for electrons in PHOS 
  Float_t   fPHOSChargeWeight;                  // Bayesian PID weight for charged hadrons in PHOS 
  Float_t   fPHOSNeutralWeight;                 // Bayesian PID weight for neutral hadrons in PHOS 
  
  Bool_t    fPHOSWeightFormula ;                // Use parametrized weight threshold, function of energy
  TFormula *fPHOSPhotonWeightFormula ;          // Formula for photon weight
  TFormula *fPHOSPi0WeightFormula ;             // Formula for pi0 weight
  TString   fPHOSPhotonWeightFormulaExpression; // Photon weight formula in string
  TString   fPHOSPi0WeightFormulaExpression;    // Pi0 weight formula in string

  // PID calculation
  Float_t   fEMCALL0CutMax;                     // Max Cut on shower shape lambda0, used in PID evaluation, only EMCAL
  Float_t   fEMCALL0CutMin;                     // Min Cut on shower shape lambda0, used in PID evaluation, only EMCAL
  Float_t   fEMCALDEtaCut;                      // Track matching cut on Dz
  Float_t   fEMCALDPhiCut;                      // Track matching cut on Dx

  Float_t   fTOFCut;                            // Cut on TOF, used in PID evaluation
  
  Float_t   fPHOSDispersionCut;                 // Shower shape elipse radious cut
  Float_t   fPHOSRCut;                          // Track-Cluster distance cut for track matching in PHOS  
  
  // Cluster splitting mass ranges
  Bool_t    fDoClusterSplitting;                // Cluster splitting analysis
  Bool_t    fUseSimpleMassCut;                  // Use simple min-max pi0 mass cut
  Bool_t    fUseSimpleM02Cut;                   // Use simple min-max M02 cut
  Bool_t    fUseSplitAsyCut ;                   // Remove splitted clusters with too large asymmetry, range defined in AliCaloPID
  Float_t   fSplitM02MaxCut ;                   // Study clusters with l0 smaller than cut
  Float_t   fSplitM02MinCut ;                   // Study clusters with l0 larger than cut  // simple case
  Int_t     fSplitMinNCells ;                   // Study clusters with ncells larger than cut  
  Float_t   fMassEtaMin  ;                      // Min Eta mass
  Float_t   fMassEtaMax  ;                      // Max Eta mass  
  Float_t   fMassPi0Min  ;                      // Min Pi0 mass // simple cut case
  Float_t   fMassPi0Max  ;                      // Min Pi0 mass // simple cut case
  Float_t   fMassPhoMin  ;                      // Min Photon mass
  Float_t   fMassPhoMax  ;                      // Min Photon mass
  Float_t   fMassWidthPi0Param[7] ;             // 3 param for pol2 fit on width, 2 param for mass position NLM=1 and NLM>1 for pi0 selection
  Float_t   fM02MinParam[2][6] ;                // 4 param for pol3 fit on M02 minimum for pi0 selection (maximum for conversions)
  Float_t   fM02MaxParam[2][6] ;                // 4 param for pol3 fit on M02 maximum for pi0 selection
  Float_t   fAsyMinParam[2][6] ;                // 4 param for pol3 fit on asymmetry minimum, for 2 cases, NLM=1 and NLM>=2
  Float_t   fSplitEFracMin[3]  ;                // Do not use clusters with too large energy in cluster compared
                                                // to energy in splitted clusters, depeding on NLM
  Float_t   fSplitWidthSigma;                   // Cut on mass+-width*fSplitWidthSigma


  
  AliCaloPID & operator = (const AliCaloPID & cpid) ; // cpy assignment
  AliCaloPID(              const AliCaloPID & cpid) ; // cpy ctor
  
  ClassDef(AliCaloPID,17)
  
} ;


#endif //ALICALOPID_H



