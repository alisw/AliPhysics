#ifndef ALIANACALORIMETERQA_H
#define ALIANACALORIMETERQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Class to check results from simulations or reconstructed real data. 
// Fill few histograms and do some checking plots
//
//-- Author: Gustavo Conesa (INFN-LNF)a

// --- Root system ---
class TH3F;
class TH2F;
class TH1F;

// --- Analysis system --- 
class AliESDCaloCluster;
class AliAODCaloCluster;

#include "AliAnaPartCorrBaseClass.h"
 
class AliAnaCalorimeterQA : public AliAnaPartCorrBaseClass {
  
 public: 
  AliAnaCalorimeterQA() ; // default ctor	
  virtual ~AliAnaCalorimeterQA() {;} //virtual dtor
 private:
  AliAnaCalorimeterQA & operator = (const AliAnaCalorimeterQA & g) ;//cpy assignment
  AliAnaCalorimeterQA(const AliAnaCalorimeterQA & g) ; // cpy ctor
  
public:
  
  void ClusterHistograms(const TLorentzVector mom, const Double_t tof, Float_t *pos, Float_t * showerShape, 
						 const Int_t nCaloCellsPerCluster, const Int_t nModule,
						 const Int_t nTracksMatched, const TObject* track, 
						 const Int_t * labels, const Int_t nLabels);
	
  TList * GetCreateOutputObjects();

  void Init();
  void InitParameters();
  
  void Print(const Option_t * opt) const;
    
  void MakeAnalysisFillHistograms() ; 
  
  void MCHistograms(const TLorentzVector mom, const Int_t pdg);

  TString GetCalorimeter() const {return fCalorimeter ;}
  void SetCalorimeter( TString calo ) {fCalorimeter = calo; }
  TString GetStyleMacro() const {return fStyleMacro ;}
  void SetStyleMacro( TString macro ) {fStyleMacro = macro; }
  
  void SwitchOnPlotsMaking()  {fMakePlots = kTRUE;}
  void SwitchOffPlotsMaking() {fMakePlots = kFALSE;}
	
  void SwitchOnCalorimetersCorrelation()  {fCorrelateCalos = kTRUE;}
  void SwitchOffCalorimetersCorrelation() {fCorrelateCalos = kFALSE;}
  void CorrelateCalorimeters(TRefArray* caloClusters);
	
  void Terminate(TList * outputList);
  void ReadHistograms(TList * outputList); //Fill histograms with histograms in ouput list, needed in Terminate.
	
  void SetNumberOfModules(Int_t nmod) {fNModules = nmod;}

  void SetTimeCut(Double_t min, Double_t max) {fTimeCutMin = min; fTimeCutMax = max;}
  Double_t GetTimeCutMin() const {return fTimeCutMin;}
  Double_t GetTimeCutMax() const {return fTimeCutMax;}

  //Histogram binning setters

  Int_t GetNewRebinForRePlotting(TH1D*histo, const Float_t newXmin, const Float_t newXmax,
				 const Int_t newNbins) const;

  virtual void SetHistoPOverERangeAndNBins(Float_t min, Float_t max, Int_t n) {
	fHistoPOverEBins  = n ;
	fHistoPOverEMax   = max ;
	fHistoPOverEMin   = min ;
  }
	
  Int_t   GetHistoPOverEBins()  const { return fHistoPOverEBins ; }
  Float_t GetHistoPOverEMin()   const { return fHistoPOverEMin ; }
  Float_t GetHistoPOverEMax()   const { return fHistoPOverEMax ; }	
	
	virtual void SetHistodEdxRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistodEdxBins  = n ;
		fHistodEdxMax   = max ;
		fHistodEdxMin   = min ;
	}
	
	Int_t   GetHistodEdxBins()  const { return fHistodEdxBins ; }
	Float_t GetHistodEdxMin()   const { return fHistodEdxMin ; }
	Float_t GetHistodEdxMax()   const { return fHistodEdxMax ; }	

	virtual void SetHistodRRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistodRBins  = n ;
		fHistodRMax   = max ;
		fHistodRMin   = min ;
	}
	
	Int_t   GetHistodRBins()  const { return fHistodRBins ; }
	Float_t GetHistodRMin()   const { return fHistodRMin ; }
	Float_t GetHistodRMax()   const { return fHistodRMax ; }	

	virtual void SetHistoTimeRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistoTimeBins  = n ;
		fHistoTimeMax   = max ;
		fHistoTimeMin   = min ;
	}	
	
	Int_t   GetHistoTimeBins()  const { return fHistoTimeBins ; }
	Float_t GetHistoTimeMin()   const { return fHistoTimeMin ; }
	Float_t GetHistoTimeMax()   const { return fHistoTimeMax ; }	

	virtual void SetHistoNClusterCellRangeAndNBins(Int_t min, Int_t max, Int_t n) {
		fHistoNBins  = n ;
		fHistoNMax   = max ;
		fHistoNMin   = min ;
	}
	
	Int_t GetHistoNClusterCellBins()  const { return fHistoNBins ; }
	Int_t GetHistoNClusterCellMin()   const { return fHistoNMin ; }
	Int_t GetHistoNClusterCellMax()   const { return fHistoNMax ; }	

	virtual void SetHistoRatioRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistoRatioBins  = n ;
		fHistoRatioMax   = max ;
		fHistoRatioMin   = min ;
	}
	
	Int_t   GetHistoRatioBins()  const { return fHistoRatioBins ; }
	Float_t GetHistoRatioMin()   const { return fHistoRatioMin ; }
	Float_t GetHistoRatioMax()   const { return fHistoRatioMax ; }	

	virtual void SetHistoVertexDistRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistoVertexDistBins = n ;
		fHistoVertexDistMax   = max ;
		fHistoVertexDistMin   = min ;
	}
	
	Int_t   GetHistoVertexDistBins()  const { return fHistoVertexDistBins ; }
	Float_t GetHistoVertexDistMin()   const { return fHistoVertexDistMin ; }
	Float_t GetHistoVertexDistMax()   const { return fHistoVertexDistMax ; }	
	
	
	virtual void SetHistoXRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistoXBins  = n ;
		fHistoXMax   = max ;
		fHistoXMin   = min ;
	}
	
	Int_t   GetHistoXBins()  const { return fHistoXBins ; }
	Float_t GetHistoXMin()   const { return fHistoXMin ; }
	Float_t GetHistoXMax()   const { return fHistoXMax ; }	

	
	virtual void SetHistoYRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistoYBins  = n ;
		fHistoYMax   = max ;
		fHistoYMin   = min ;
	}
	
	Int_t   GetHistoYBins()  const { return fHistoYBins ; }
	Float_t GetHistoYMin()   const { return fHistoYMin ; }
	Float_t GetHistoYMax()   const { return fHistoYMax ; }	

	
	virtual void SetHistoZRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistoZBins  = n ;
		fHistoZMax   = max ;
		fHistoZMin   = min ;
	}
	
	Int_t   GetHistoZBins()  const { return fHistoZBins ; }
	Float_t GetHistoZMin()   const { return fHistoZMin ; }
	Float_t GetHistoZMax()   const { return fHistoZMax ; }	
	
	virtual void SetHistoRRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistoRBins  = n ;
		fHistoRMax   = max ;
		fHistoRMin   = min ;
	}
	
	Int_t   GetHistoRBins()  const { return fHistoRBins ; }
	Float_t GetHistoRMin()   const { return fHistoRMin ; }
	Float_t GetHistoRMax()   const { return fHistoRMax ; }	
	
	virtual void SetHistoShowerShapeRangeAndNBins(Float_t min, Float_t max, Int_t n) {
		fHistoSSBins  = n ;
		fHistoSSMax   = max ;
		fHistoSSMin   = min ;
	}
	
	Int_t   GetHistoShowerShapeBins()  const { return fHistoSSBins ; }
	Float_t GetHistoShowerShapeMin()   const { return fHistoSSMin ; }
	Float_t GetHistoShowerShapeMax()   const { return fHistoSSMax ; }	
	
	
 private:
  
  TString  fCalorimeter ;    // Calorimeter selection
  TString  fStyleMacro  ;    // Location of macro for plots style
  Bool_t   fMakePlots   ;    // Print plots
  Bool_t   fCorrelateCalos ; // Correlate PHOS/EMCAL clusters
  Int_t    fNModules    ;    // Number of EMCAL/PHOS modules, set as many histogras as modules 
  Int_t    fNRCU        ;    // Number of EMCAL/PHOS RCU, set as many histogras as RCU 
  Double_t fTimeCutMin  ;    // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;    // Remove clusters/cells with time larger than this value, in ns
	
  //Histograms
  //Histogram Bins
  Int_t   fHistoPOverEBins;        // p/E histogram number of bins
  Float_t fHistoPOverEMax;         // p/E maximum value
  Float_t fHistoPOverEMin;         // p/E minimum value
  Int_t   fHistodEdxBins;          // dEdx histogram number of bins
  Float_t fHistodEdxMax;           // dEdx maximum value
  Float_t fHistodEdxMin;           // dEdx minimum value
  Int_t   fHistodRBins;            // dR histogram number of bins
  Float_t fHistodRMax;             // dR maximum value
  Float_t fHistodRMin;             // dR minimum value
  Int_t   fHistoTimeBins;          // cell time histogram number of bins
  Float_t fHistoTimeMax;           // cell time maximum value
  Float_t fHistoTimeMin;           // cell time minimum value
  Int_t   fHistoNBins;             // number of clusters/cells histogram number of bins
  Int_t   fHistoNMax;              // number maximum value
  Int_t   fHistoNMin;              // number minimum value
  Int_t   fHistoRatioBins;         // ratio histogram number of bins
  Float_t fHistoRatioMax;          // ratio maximum value
  Float_t fHistoRatioMin;          // ratio minimum value
  Int_t   fHistoVertexDistBins;    // vertex distance histogram number of bins
  Float_t fHistoVertexDistMax;     // vertex distance maximum value
  Float_t fHistoVertexDistMin;     // vertex distance minimum value	
  Int_t   fHistoRBins;             // r =sqrt(x^2+y^2+z^2) (cm) position histogram number of bins
  Float_t fHistoRMax;              // r =sqrt(x^2+y^2+z^2) (cm)  maximum value
  Float_t fHistoRMin;              // r =sqrt(x^2+y^2+z^2) (cm)  minimum value	
  Int_t   fHistoXBins;             // x (cm) position histogram number of bins
  Float_t fHistoXMax;              // x (cm) position maximum value
  Float_t fHistoXMin;              // x (cm) position minimum value
  Int_t   fHistoYBins;             // y (cm) position histogram number of bins
  Float_t fHistoYMax;              // y (cm) position maximum value
  Float_t fHistoYMin;              // y (cm) position minimum value
  Int_t   fHistoZBins;             // z (cm) position histogram number of bins
  Float_t fHistoZMax;              // z (cm) position maximum value
  Float_t fHistoZMin;              // z (cm) position minimum value
  Int_t   fHistoSSBins;            // Shower Shape parameter histogram number of bins
  Float_t fHistoSSMax;             // Shower Shape parameter position maximum value
  Float_t fHistoSSMin;             // Shower Shape parameter position minimum value
	
  //CaloClusters 
  TH1F * fhE  ; //! E distribution, Reco
  TH1F * fhPt ; //! pT distribution, Reco
  TH1F * fhPhi; //! phi distribution, Reco 
  TH1F * fhEta; //! eta distribution, Reco 
  TH3F * fhEtaPhiE  ; //! eta vs phi vs E, Reco
  TH1F * fhECharged  ; //! E distribution, Reco, matched with track
  TH1F * fhPtCharged ; //! pT distribution, Reco, matched with track
  TH1F * fhPhiCharged; //! phi distribution, Reco, matched with track 
  TH1F * fhEtaCharged; //! eta distribution, Reco, matched with track 
  TH3F * fhEtaPhiECharged  ; //! eta vs phi vs E, Reco, matched with track 
//  TH1F * fhEChargedNoOut  ; //! E distribution, Reco, matched with track, no outer param
//  TH1F * fhPtChargedNoOut ; //! pT distribution, Reco, matched with track, no outer param
//  TH1F * fhPhiChargedNoOut; //! phi distribution, Reco, matched with track, no outer param
//  TH1F * fhEtaChargedNoOut; //! eta distribution, Reco, matched with track, no outer param
//  TH2F * fhEtaPhiChargedNoOut  ; //! eta vs phi, Reco, matched with track, no outer param
  TH1F * fhDeltaE  ; //! MC-Reco E distribution	
  TH1F * fhDeltaPt ; //! MC-Reco pT distribution
  TH1F * fhDeltaPhi; //! MC-Reco phi distribution
  TH1F * fhDeltaEta; //! MC-Reco eta distribution
  TH1F * fhRatioE  ; //! Reco/MC E distribution	
  TH1F * fhRatioPt ; //! Reco/MC pT distribution
  TH1F * fhRatioPhi; //! Reco/MC phi distribution
  TH1F * fhRatioEta; //! Reco/MC eta distribution
  TH2F * fh2E  ; //! E distribution, Reco vs MC
  TH2F * fh2Pt ; //! pT distribution, Reco vs MC
  TH2F * fh2Phi; //! phi distribution, Reco vs MC
  TH2F * fh2Eta; //! eta distribution, Reco vs MC
  
  TH3F * fhLambda  ;     //! Shower ellipse axis Lambda 0 vs vs Lambda 1 vs E
  TH2F * fhDispersion  ; //! Shower dispersion vs E
	
  TH2F * fhIM; //! cluster pairs invariant mass
  TH2F * fhIMCellCut; //! cluster pairs invariant mass, n cells > 1 per cluster
  TH2F * fhAsym; //! cluster pairs invariant mass	
  
  TH3F * fhNCellsPerCluster;           //! N cells per cluster vs cluster energy vs eta of cluster	
  TH3F * fhNCellsPerClusterMIP;        //! N cells per cluster vs cluster energy vs eta of cluster, finer fixed pT bin for MIP search.
  TH3F * fhNCellsPerClusterMIPCharged; //! N cells per cluster vs cluster energy vs eta of cluster, finer fixed pT bin for MIP search, cluster matched with track.	
  	
  TH1F * fhNClusters; //! Number of clusters
	
  TH2F * fhClusterTimeEnergy;   //! Cluster Time vs Energy 
  TH1F * fhCellTimeSpreadRespectToCellMax; //! Difference of the time of cell with maximum dep energy and the rest of cells
  TH1F * fhCellIdCellLargeTimeSpread;      //! Cells with large time respect to max (diff > 100 ns)
	
  TH2F * fhRNCells ; //! R=sqrt(x^2+y^2+z^2) (cm) cluster distribution vs N cells in cluster
  TH2F * fhXNCells ; //! X (cm) cluster distribution vs N cells in cluster
  TH2F * fhYNCells ; //! Y (cm) cluster distribution vs N cells in cluster
  TH2F * fhZNCells ; //! Z (cm) cluster distribution vs N cells in cluster
	
  TH2F * fhRE ; //! R=sqrt(x^2+y^2+z^2) (cm) cluster distribution vs cluster energy
  TH2F * fhXE ; //! X (cm) cluster distribution vs cluster energy
  TH2F * fhYE ; //! Y (cm) cluster distribution vs cluster energy
  TH2F * fhZE ; //! Z (cm) cluster distribution vs cluster energy
  TH3F * fhXYZ; //! cluster X vs Y vs Z (cm)
	
  TH2F * fhRCellE ; //! R=sqrt(x^2+y^2+z^2) (cm) cell distribution vs cell energy
  TH2F * fhXCellE ; //! X (cm) cell distribution vs cell energy
  TH2F * fhYCellE ; //! Y (cm) cell distribution vs cell energy
  TH2F * fhZCellE ; //! Z (cm) cell distribution vs cell energy
  TH3F * fhXYZCell; //! cell X vs Y vs Z (cm)
		
  TH2F * fhDeltaCellClusterRNCells ; //! R cluster - R cell distribution (cm) vs N cells in cluster
  TH2F * fhDeltaCellClusterXNCells ; //! X cluster - X cell distribution (cm) vs N cells in cluster
  TH2F * fhDeltaCellClusterYNCells ; //! Y cluster - Y cell distribution (cm) vs N cells in cluster
  TH2F * fhDeltaCellClusterZNCells ; //! Z cluster - Z cell distribution (cm) vs N cells in cluster
	
  TH2F * fhDeltaCellClusterRE ; //! R cluster - R cell distribution (cm) vs cluster energy
  TH2F * fhDeltaCellClusterXE ; //! X cluster - X cell distribution (cm) vs cluster energy
  TH2F * fhDeltaCellClusterYE ; //! Y cluster - Y cell distribution (cm) vs cluster energy
  TH2F * fhDeltaCellClusterZE ; //! Z cluster - Z cell distribution (cm) vs cluster energy
	
	
  //Calo Cells
  TH1F * fhNCells;    //! Number of towers/crystals with signal
  TH1F * fhAmplitude; //! Amplitude measured in towers/crystals
  TH2F * fhAmpId;     //! Amplitude measured in towers/crystals vs id of tower.	
  TH3F * fhEtaPhiAmp; //! eta vs phi vs amplitude, cells

  TH1F * fhTime;      //! Time measured in towers/crystals
  TH2F * fhTimeId;    //! Time vs Absolute cell Id
  TH2F * fhTimeAmp;   //! Time vs Amplitude 
//  TH1F * fhT0Time;    //! T0 - EMCAL Time measured in towers/crystals
//  TH2F * fhT0TimeId;  //! T0 - EMCAL Time vs Absolute cell Id
//  TH2F * fhT0TimeAmp; //! T0 - EMCAL Time vs Amplitude 


  //Calorimeters Correlation
  TH2F * fhCaloCorrNClusters; // EMCAL vs PHOS, number of clusters	
  TH2F * fhCaloCorrEClusters; // EMCAL vs PHOS, total measured cluster energy
  TH2F * fhCaloCorrNCells; // EMCAL vs PHOS, number of cells
  TH2F * fhCaloCorrECells; // EMCAL vs PHOS,  total measured cell energy
	
  //Module histograms
  TH1F ** fhEMod  ;               //! E distribution for different module, Reco
  TH1F ** fhNClustersMod ;        //! Number of clusters for different module, Reco
  TH2F ** fhNCellsPerClusterMod ; //! N cells per clusters different module, Reco
  TH1F ** fhNCellsMod ;           //! Number of towers/crystals with signal different module, Reco
  TH2F ** fhGridCellsMod ;        //! Cells ordered in column/row for different module, Reco
  TH2F ** fhGridCellsEMod ;       //! Cells ordered in column/row for different module, weighted with energy, Reco
  TH2F ** fhGridCellsTimeMod ;    //! Cells ordered in column/row for different module, weighted with time, Reco
  TH1F ** fhAmplitudeMod ;        //! Amplitude measured in towers/crystals different module, Reco
  TH1F ** fhAmplitudeModFraction; //! Amplitude measured in towers/crystals different fractions of module, Reco
  TH2F ** fhTimeAmpPerRCU;        //! Time vs Amplitude measured in towers/crystals different RCU
  //TH2F ** fhT0TimeAmpPerRCU;      //! T0 - EMCAL Time vs Amplitude measured in towers/crystals different RCU
  //TH2F ** fhTimeCorrRCU;          //! Correlate time entries in the different RCU, E > 0.3
  TH2F ** fhIMMod;                //! cluster pairs invariant mass, different module,
  TH2F ** fhIMCellCutMod;         //! cluster pairs invariant mass, n cells > 1 per cluster, different module
	
  //MC  
  TH1F *fhGenGamPt  ; // pt of primary gamma
  TH1F *fhGenGamEta ; // eta of primart gamma
  TH1F *fhGenGamPhi ; // phi of primary gamma	
  TH1F *fhGenPi0Pt  ; // pt of primary pi0
  TH1F *fhGenPi0Eta ; // eta of primart pi0
  TH1F *fhGenPi0Phi ; // phi of primary pi0	
  TH1F *fhGenEtaPt  ; // pt of primary eta
  TH1F *fhGenEtaEta ; // eta of primart eta
  TH1F *fhGenEtaPhi ; // phi of primary eta
  TH1F *fhGenOmegaPt  ; // pt of primary omega
  TH1F *fhGenOmegaEta ; // eta of primart omega
  TH1F *fhGenOmegaPhi ; // phi of primary omega
  TH1F *fhGenElePt  ; // pt of primary electron
  TH1F *fhGenEleEta ; // eta of primart electron
  TH1F *fhGenElePhi ; // phi of primary electron
	
  //TH3F * fhEMVxyz    ; // Electromagnetic particle production vertex
  TH2F * fhEMVxyz    ; // Electromagnetic particle production vertex
  TH2F * fhEMR       ; // Electromagnetic distance to vertex vs rec energy  
  //TH3F * fhHaVxyz    ; // Hadron production vertex
  TH2F * fhHaVxyz    ; // Hadron production vertex
  TH2F * fhHaR       ; // Hadron distance to vertex vs rec energy  

  TH2F * fhGamE  ; //! E distribution of generated photons, Reco
  TH2F * fhGamPt ; //! pT distribution of generated photons, Reco
  TH2F * fhGamPhi; //! phi distribution of generated photon, Reco 
  TH2F * fhGamEta; //! eta distribution of generated photons, Reco 
  TH1F * fhGamDeltaE  ; //! MC-Reco E distribution of generated photons	
  TH1F * fhGamDeltaPt ; //! MC-Reco pT distribution of generated photons
  TH1F * fhGamDeltaPhi; //! MC-Reco phi distribution of generated photons
  TH1F * fhGamDeltaEta; //! MC-Reco eta distribution of generated photons
  TH1F * fhGamRatioE  ; //! Reco/MC E distribution of generated photons	
  TH1F * fhGamRatioPt ; //! Reco/MC pT distribution of generated photons
  TH1F * fhGamRatioPhi; //! Reco/MC phi distribution of generated photons
  TH1F * fhGamRatioEta; //! Reco/MC eta distribution of generated photons
  TH2F * fhEleE  ; //! E distribution of generated electrons, Reco
  TH2F * fhElePt ; //! pT distribution of generated electrons, Reco
  TH2F * fhElePhi; //! phi distribution of generated electron, Reco 
  TH2F * fhEleEta; //! eta distribution of generated electrons, Reco 		
  TH2F * fhPi0E  ; //! E distribution of generated pi0, Reco, gamma decay overlapped
  TH2F * fhPi0Pt ; //! pT distribution of generated pi0, Reco, gamma decay overlapped
  TH2F * fhPi0Phi; //! phi distribution of generated pi0, Reco, gamma decay overlapped
  TH2F * fhPi0Eta; //! eta distribution of generated pi0, Reco, gamma decay overlapped
  TH2F * fhNeHadE  ; //! E distribution of generated neutral hadron, Reco
  TH2F * fhNeHadPt ; //! pT distribution of generated neutral hadron, Reco
  TH2F * fhNeHadPhi; //! phi distribution of generated neutral hadron, Reco 
  TH2F * fhNeHadEta; //! eta distribution of generated neutral hadron, Reco 	
  TH2F * fhChHadE  ; //! E distribution of generated charged hadron, Reco
  TH2F * fhChHadPt ; //! pT distribution of generated charged hadron, Reco
  TH2F * fhChHadPhi; //! phi distribution of generated charged hadron, Reco 
  TH2F * fhChHadEta; //! eta distribution of generated charged hadron, Reco 

  TH2F * fhGamECharged  ; //! E distribution of generated photons, Reco, track matched cluster
  TH2F * fhGamPtCharged ; //! pT distribution of generated photons, Reco, track matched cluster
  TH2F * fhGamPhiCharged; //! phi distribution of generated photon, Reco, track matched cluster 
  TH2F * fhGamEtaCharged; //! eta distribution of generated photons, Reco, track matched cluster 
  TH2F * fhEleECharged  ; //! E distribution of generated electrons, Reco, track matched cluster
  TH2F * fhElePtCharged ; //! pT distribution of generated electrons, Reco, track matched cluster
  TH2F * fhElePhiCharged; //! phi distribution of generated electron, Reco, track matched cluster 
  TH2F * fhEleEtaCharged; //! eta distribution of generated electrons, Reco, track matched cluster 		
  TH2F * fhPi0ECharged  ; //! E distribution of generated pi0, Reco, gamma decay overlapped, track matched cluster
  TH2F * fhPi0PtCharged ; //! pT distribution of generated pi0, Reco, gamma decay overlapped, track matched cluster
  TH2F * fhPi0PhiCharged; //! phi distribution of generated pi0, Reco, gamma decay overlapped, track matched cluster
  TH2F * fhPi0EtaCharged; //! eta distribution of generated pi0, Reco, gamma decay overlapped, track matched cluster
  TH2F * fhNeHadECharged  ; //! E distribution of generated neutral hadron, Reco, track matched cluster
  TH2F * fhNeHadPtCharged ; //! pT distribution of generated neutral hadron, Reco, track matched cluster
  TH2F * fhNeHadPhiCharged; //! phi distribution of generated neutral hadron, Reco , track matched cluster
  TH2F * fhNeHadEtaCharged; //! eta distribution of generated neutral hadron, Reco, track matched cluster 	
  TH2F * fhChHadECharged  ; //! E distribution of generated charged hadron, Reco, track matched cluster
  TH2F * fhChHadPtCharged ; //! pT distribution of generated charged hadron, Reco, track matched cluster
  TH2F * fhChHadPhiCharged; //! phi distribution of generated charged hadron, Reco, track matched cluster 
  TH2F * fhChHadEtaCharged; //! eta distribution of generated charged hadron, Reco, track matched cluster 	
	
  TH1F *fhGenGamAccE   ; // E of primary gamma
  TH1F *fhGenGamAccPt  ; // pt of primary gamma
  TH1F *fhGenGamAccEta ; // eta of primart gamma
  TH1F *fhGenGamAccPhi ; // phi of primary gamma	
  TH1F *fhGenPi0AccE   ; // E of primary pi0
  TH1F *fhGenPi0AccPt  ; // pt of primary pi0
  TH1F *fhGenPi0AccEta ; // eta of primart pi0
  TH1F *fhGenPi0AccPhi ; // phi of primary pi0		
	
  //Histograms for track-matching
  TH2F *fh1pOverE;     //! p/E for track-cluster matches
  TH1F *fh1dR;         //! distance between projected track and cluster
  TH2F *fh2EledEdx;    //! dE/dx vs. momentum for electron candidates
  TH2F *fh2MatchdEdx;  //! dE/dx vs. momentum for all matches
	
  TH2F *fhMCEle1pOverE;     //! p/E for track-cluster matches, MC electrons
  TH1F *fhMCEle1dR;         //! distance between projected track and cluster, MC electrons
  TH2F *fhMCEle2MatchdEdx;  //! dE/dx vs. momentum for all matches, MC electrons	
	
  TH2F *fhMCChHad1pOverE;     //! p/E for track-cluster matches, MC charged hadrons
  TH1F *fhMCChHad1dR;         //! distance between projected track and cluster, MC charged hadrons
  TH2F *fhMCChHad2MatchdEdx;  //! dE/dx vs. momentum for all matches, MC charged
	
  TH2F *fhMCNeutral1pOverE;     //! p/E for track-cluster matches, MC neutral
  TH1F *fhMCNeutral1dR;         //! distance between projected track and cluster, MC neutral
  TH2F *fhMCNeutral2MatchdEdx;  //! dE/dx vs. momentum for all matches, MC neutral	
	
  TH2F *fh1pOverER02;           //! p/E for track-cluster matches, dR > 0.2	
  TH2F *fhMCEle1pOverER02;      //! p/E for track-cluster matches, dR > 0.2, MC electrons
  TH2F *fhMCChHad1pOverER02;    //! p/E for track-cluster matches, dR > 0.2, MC charged hadrons
  TH2F *fhMCNeutral1pOverER02;  //! p/E for track-cluster matches, dR > 0.2, MC neutral
	
	ClassDef(AliAnaCalorimeterQA,9)
} ;


#endif //ALIANACALORIMETERQA_H



