#ifndef ALIANASHOWERPARAMETER_H
#define ALIANASHOWERPARAMETER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnaShowerParameter.h 27413 2008-07-18 13:28:12Z gconesab $ */

//_________________________________________________________________________
//
// Class cloned from AliAnaPhoton, main aim is shower shape studies
// 
// 
//
//-- Author: Jocelyn Mlynarz (WSU) and Gustavo Conesa (LPSC)

// --- ROOT system ---
class TH3F;
class TH2F ;
class TH1F;
class TString ;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"
//#include "AliStack.h"
//#include "TParticle.h"
class AliStack;
class TParticle;

class TList ;
class AliEMCALGeoUtils;
class AliAnaShowerParameter : public AliAnaPartCorrBaseClass {

public: 

  AliAnaShowerParameter() ; // default ctor
  AliAnaShowerParameter(const AliAnaShowerParameter & g) ; // cpy ctor
  AliAnaShowerParameter & operator = (const AliAnaShowerParameter & g) ;//cpy assignment
  virtual ~AliAnaShowerParameter() ; //virtual dtor
  
  TList *  GetCreateOutputObjects();

  void Init();

  TObjString* GetAnalysisCuts();

  void MakeAnalysisFillAOD()  ;
    
  void MakeAnalysisFillHistograms() ; 
  
  void Print(const Option_t * opt)const;
  
  TString GetCalorimeter()   const {return fCalorimeter ; }
  void SetCalorimeter(TString det)    {fCalorimeter = det ; }

  Bool_t IsTrackMatchRejectionOn()   const {return fRejectTrackMatch ; }
  void SwitchOnTrackMatchRejection()  {fRejectTrackMatch = kTRUE ; }
  void SwitchOffTrackMatchRejection() {fRejectTrackMatch = kFALSE ; }  

  Bool_t IsCheckConversionOn()   const {return fCheckConversion ; }
  void SwitchOnConversionChecker()  {fCheckConversion = kTRUE ; }
  void SwitchOffConversionChecker() {fCheckConversion = kFALSE ; }  
	
  Bool_t AreConvertedPairsInAOD()   const {return fAddConvertedPairsToAOD ; }
  void SwitchOnAdditionConvertedPairsToAOD()  {fAddConvertedPairsToAOD = kTRUE ; }
  void SwitchOffAdditionConvertedPairsToAOD() {fAddConvertedPairsToAOD = kFALSE ; }  
	
  void InitParameters();
 
  void SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) {
    fMinDist = m1;
    fMinDist2 = m2;
    fMinDist3 = m3;
  }
	
  Float_t GetMassCut()    const {return fMassCut ; }
  void SetMassCut(Float_t m)    {fMassCut = m ; }
  void SetNClusterCut(Int_t nc)   {fNumClusters = nc;}
  Int_t GetNClusterCut()   {return fNumClusters;}
  void SetTimeCut(Double_t min, Double_t max) {fTimeCutMin = min; fTimeCutMax = max;}
  Double_t GetTimeCutMin() const {return fTimeCutMin;}
  Double_t GetTimeCutMax() const {return fTimeCutMax;}	
	
  private:
 
  TString fCalorimeter ; // Calorimeter where the gamma is searched;
  Float_t fMinDist ;     // Minimal distance to bad channel to accept cluster
  Float_t fMinDist2;     // Cuts on Minimal distance to study acceptance evaluation
  Float_t fMinDist3;     // One more cut on distance used for acceptance-efficiency study
  Bool_t  fRejectTrackMatch ;      //If PID on, reject clusters which have an associated TPC track
  Bool_t  fCheckConversion;        // Combine pairs of clusters with mass close to 0
  Bool_t  fAddConvertedPairsToAOD; // Put Converted pairs in AOD
  Float_t fMassCut;                // Mass cut for the conversion pairs selection
  Float_t fNCellsCut;
  Double_t fTimeCutMin  ;    // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;    // Remove clusters/cells with time larger than this value, in ns
  AliEMCALGeoUtils * fEMCALGeo;
  Int_t fNumClusters;        //Cut that selects events with a specific number of clusters

  //Histograms   
  TH1F * fhNClusters  ; //! Number of clusters per event.
  TH2F * fhNCellCluster;//! Number of cells per cluster
  TH1F * fhPtCluster   ; //! Number of clusters vs transerse momentum 
  TH2F * fhPhiCluster  ; //! Azimuthal angle of clusters vs transerse momentum 
  TH2F * fhEtaCluster  ; //! Pseudorapidity of clusters vs transerse momentum
  TH1F * fhDeltaPhiClusters  ; //! Delta phi of cluster pairs
  TH1F * fhDeltaEtaClusters  ; //! Delta eta of cluster pairs 
  TH3F * fhLambdaCluster  ; //! Shower parameters of clusters vs transerse momentum
  TH2F * fhDispersionCluster  ; //! Dispersion of the clusters
  TH3F * fhELambdaCluster  ; //! Shower parameters of clusters vs Energy
  TH3F * fhELambdaCellCluster  ; //! Shower parameters of clusters vs Energy

  TH2F * fhNCellPhoton;//! Number of cells per photon
  TH3F * fhLambdaPhoton  ; //! Shower parameters of photons vs transerse momentum
  TH2F * fhDispersionPhoton  ; //! Dispersion of the photons
  TH3F * fhELambdaPhoton  ; //! Shower parameters of photons vs Energy
  TH3F * fhELambdaCellPhoton  ; //! Shower parameters of photons vs Energy

  TH2F * fhNCellPi0;//! Number of cells per neutral pion
  TH3F * fhLambdaPi0  ; //! Shower parameters of neutral pions vs transerse momentum
  TH2F * fhDispersionPi0  ; //! Dispersion of the neutral pions
  TH3F * fhELambdaPi0  ; //! Shower parameters of neutral pions vs Energy
  TH3F * fhELambdaCellPi0  ; //! Shower parameters of neutral pions vs Energy

  TH2F * fhNCellChargedHadron;//! Number of cells per charged hadron
  TH3F * fhLambdaChargedHadron  ; //! Shower parameters of charged hadrons vs transerse momentum
  TH2F * fhDispersionChargedHadron  ; //! Dispersion of the charged hadrons
  TH3F * fhELambdaChargedHadron  ; //! Shower parameters of charged hadrons vs Energy
  TH3F * fhELambdaCellChargedHadron  ; //! Shower parameters of charged hadrons vs Energy

  //MC
  TH1F * fhDeltaE  ; //! MC-Reco E distribution      
  TH1F * fhDeltaPt ; //! MC-Reco pT distribution
  TH1F * fhRatioE  ; //! Reco/MC E distribution      
  TH1F * fhRatioPt ; //! Reco/MC pT distribution
  TH2F * fh2E  ; //! E distribution, Reco vs MC
  TH2F * fh2Pt ; //! pT distribution, Reco vs MC
  TH2F * fhMCPdg; //! Complete list of PDG Codes.
  
  TH1F * fhEMCPhoton;   //! Number of identified gamma 
  TH2F * fhPhiMCPhoton;  //! Phi of identified gamma
  TH2F * fhEtaMCPhoton;  //! eta of identified gamma	
  TH3F * fhLambdaMCPhoton  ; //! Shower parameters of MC photons vs transerse momentum 
  TH3F * fhPhotTrueE;  // MC truth E vs Recons E vs. Lambda of the cluster for MC photons
  TH1F * fhRatioEPhoton  ; //! Reco/MC E distribution for photons      
  
  TH1F * fhEMCPi0;   //! Number of identified  (single shower) Pi0 
  TH2F * fhPhiMCPi0;  //! Phi of identified Pi0
  TH2F * fhEtaMCPi0;  //! eta of identified Pi0	
  TH3F * fhLambdaMCPi0  ; //! Shower parameters of MC Pi0 vs transerse momentum   
  TH3F * fhPi0TrueE;  // MC truth E vs Recons E vs. Lambda of the cluster for MC Pi0s
  TH2F * fhProductionDistance; //! Distance from beam to production of the Pi0
  TH2F * fhProductionRadius; //! R from beam to Pi0
  TH2F * fhDecayAngle;//! Decay angle of the Pi0
  TH1F * fhRatioEPi0  ; //! Reco/MC E distribution for Pi0s

  TH1F * fhEMCPion;   //! Number of identified pions 
  TH2F * fhPhiMCPion;  //! Phi of identified pions
  TH2F * fhEtaMCPion;  //! eta of identified pions	
  TH3F * fhLambdaMCPion; //! Shower parameters of MC pions vs transerse momentum   
  TH3F * fhPionTrueE;  // MC truth E vs Recons E vs. Lambda of the cluster for MC pions
  TH1F * fhRatioEPion  ; //! Reco/MC E distribution for charged pions      

  TH1F * fhEMCProton;   //! Number of identified protons
  TH2F * fhPhiMCProton;  //! Phi of identified protons
  TH2F * fhEtaMCProton;  //! eta of identified protons	
  TH3F * fhLambdaMCProton  ; //! Shower parameters of MC protons vs transerse momentum  
  TH3F * fhProtonTrueE;  // MC truth E vs Recons E vs. Lambda of the cluster for MC protons
  TH1F * fhRatioEProton  ; //! Reco/MC E distribution for protons       

  TH1F * fhEMCAntiProton;   //! Number of identified antiprotons
  TH2F * fhPhiMCAntiProton;  //! Phi of identified antiprotons
  TH2F * fhEtaMCAntiProton;  //! eta of identified antiprotons	
  TH3F * fhLambdaMCAntiProton  ; //! Shower parameters of MC antiprotons vs transerse momentum 
  TH3F * fhAntiProtonTrueE;  // MC truth E vs Recons E vs. Lambda of the cluster for MC antiprotons
  TH1F * fhRatioEAntiProton  ; //! Reco/MC E distribution for antiprotons        

  TH1F * fhEMCNeutron;   //! Number of identified neutrons
  TH2F * fhPhiMCNeutron;  //! Phi of identified neutrons
  TH2F * fhEtaMCNeutron;  //! eta of identified neutrons	
  TH3F * fhLambdaMCNeutron  ; //! Shower parameters of MC neutrons vs transerse momentum
  TH3F * fhNeutronTrueE;  // MC truth E vs Recons E vs. Lambda of the cluster for MC Neutrons
  TH1F * fhRatioENeutron  ; //! Reco/MC E distribution for photons         

  TH1F * fhEMCEta;   //! Number of identified etas
  TH2F * fhPhiMCEta;  //! Phi of identified etas
  TH2F * fhEtaMCEta;  //! eta of identified etas	
  TH3F * fhLambdaMCEta  ; //! Shower parameters of MC etas vs transerse momentum 

  TH1D *  fhPrimPt; //! Pi0 pT spectrum truth.

   ClassDef(AliAnaShowerParameter,1)

} ;
 

#endif//AliAnaShowerParameter_H



