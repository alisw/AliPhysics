#ifndef ALIANACALOEXOTICS_H
#define ALIANACALOEXOTICS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaCaloExotics
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Class to study exotic clusters
///
/// Task filling basic histograms to check cluster/cell exoticity.
/// Extracted/based on AliAnaCalorimeterQA
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///_________________________________________________________________________

// --- Root system ---
class TH3F;
class TH2F;
class TH1F;
class TObjString;
class TObjArray;

// --- Analysis system --- 
class AliVCaloCells;
class AliVCaloCluster;
class AliVTrack;

#include "AliAnaCaloTrackCorrBaseClass.h"
 
class AliAnaCaloExotics : public AliAnaCaloTrackCorrBaseClass {
  
public:
    
  AliAnaCaloExotics() ;
    
  /// Virtual destructor. Not implemented.
  virtual ~AliAnaCaloExotics() { ; }
    
  // General methods
  
  TObjString * GetAnalysisCuts();

  TList      * GetCreateOutputObjects();
  
  void         Init();
  
  void         MakeAnalysisFillHistograms() ;
  
  void         Print(const Option_t * opt) const;
    
  // Main methods
  
  void         CellHistograms(AliVCaloCells * cells);
  
  void         ClusterHistograms(const TObjArray * clusters, AliVCaloCells * cells);
  
  // Setters and getters
  
  Float_t      GetCellAmpMin()      const   { return fCellAmpMin          ; } 
  void         SetCellAmpMin(Float_t amp)   { fCellAmpMin = amp           ; }
  
  Float_t      GetEMinForExo()      const   { return fEMinForExo          ; }
  void         SetEMinForExo(Float_t min)   { fEMinForExo = min           ; }
 
  Float_t      GetExoCut()          const   { return fExoCut              ; }
  void         SetExoCut(Float_t exo)       { fExoCut = exo               ; }
 
  Float_t      GetEBinLimit(Int_t i) const  
                           { if ( i < fgkNEBins && i >= 0 ) return fEnergyBins[i] ;
                             else                           return -1     ; }
  void         SetEBinLimit(Int_t i, Float_t en) 
                           { if ( i < fgkNEBins && i >= 0 ) fEnergyBins[i] = en ; }
  
  void         SwitchOnFill1CellHisto()     { fFill1CellHisto    = kTRUE  ; }
  void         SwitchOffFill1CellHisto()    { fFill1CellHisto    = kFALSE ; }
 
  void         SwitchOnFillCellHisto()      { fFillCellHisto     = kTRUE  ; }
  void         SwitchOffFillCellHisto()     { fFillCellHisto     = kFALSE ; }
  
  void         SwitchOnFillMatchingHisto()  { fFillMatchingHisto = kTRUE  ; }
  void         SwitchOffFillMatchingHisto() { fFillMatchingHisto = kFALSE ; }
  
  void         SetConstantTimeShift(Float_t shift) { fConstantTimeShift = shift  ; }
  
 private:
      
 
  // Cell amplitude cut
  
  Float_t  fCellAmpMin;                         ///<  Amplitude Threshold on calorimeter cells 
  
  Float_t  fEMinForExo;                         ///<  Minimum energy cut for integrated histograms over E. Exotics start at about 4 GeV.

  Float_t  fExoCut;                             ///<  Exoticity cut for some histograms 
  
  /// Total number of cluster energy bins histograms
  static const Int_t fgkNEBins = 12;
  
  Float_t  fEnergyBins[fgkNEBins];              ///<  Energy bins for some histograms
  
  Bool_t   fFillCellHisto;                      ///<  Fill histograms single cells
 
  Bool_t   fFill1CellHisto;                     ///<  Fill histograms for 1 cell clusters
  
  Bool_t   fFillMatchingHisto;                  ///<  Fill histograms for track-cluster matching
  
  Float_t  fConstantTimeShift;                  ///<  Apply a 600 ns time shift in case of simulation, shift in ns.
  
  TLorentzVector fClusterMomentum;              //!<! Cluster momentum, temporary container
  
  // Calorimeter Clusters
    
  TH2F *   fhExoticityEClus;                    //!<! Exoticity vs cluster energy
  TH2F *   fhExoticityEMaxCell;                 //!<! Exoticity vs energy of highest energy cell in cluster
  TH2F *   fhExoticityEClusTrackMatch;          //!<! Exoticity vs cluster energy, track matched
  TH2F *   fhExoticity1Cell;                    //!<! Exoticity vs energy for 1 cell clusters

  TH2F *   fhNCellsPerCluster;                  //!<! Cluster energy vs N cells in cluster
  TH2F *   fhNCellsPerClusterTrackMatch;        //!<! Cluster energy vs N cells in cluster, for track-matched clusters 
  TH3F *   fhNCellsPerClusterExo;               //!<! Cluster energy vs N cells in cluster vs Exoticity   
  TH3F *   fhNCellsPerClusterExoTrackMatch;     //!<! Cluster energy vs N cells in cluster vs Exoticity, for track-matched clusters 
  TH3F *   fhNCellsPerClusterM02;               //!<! Cluster energy vs N cells in cluster vs M02   

  TH3F *   fhEtaPhiGridExoEnCut  ;              //!<! column vs row vs exoticity when E > fEMinForExo and n cells > 1
  TH3F *   fhEtaPhiGridEnExoCut  ;              //!<! column vs row vs energy when F+ < 0.97 and n cells > 1
  TH3F *   fhEtaPhiGridEn1Cell;                 //!<! column vs row vs energy for 1 cell clusters 
  
  TH3F *   fhTimeEnergyExo;                     //!<! Cluster Energy vs Time vs Exoticity, n cells > 1
  TH2F *   fhTimeEnergy1Cell;                   //!<! Cluster Energy vs Time vs n cells = 1
  TH3F *   fhTimeDiffClusCellExo;               //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs cluster energy vs exoticity
  TH3F *   fhTimeDiffAmpClusCellExo;            //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs secondary cell energy vs exoticity for E > fEMinForExo
  TH3F *   fhTimeEnergyM02;                     //!<! Cluster Energy vs Time vs M02, n cells > 1
  TH3F *   fhTimeDiffClusCellM02;               //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs cluster energy vs M02

  TH3F *   fhM02EnergyNCell;                    //!<! Cluster M02 vs Energy vs n cells
  TH3F *   fhM02EnergyExo;                      //!<! Cluster M02 vs Energy vs exoticity
  TH3F *   fhM02EnergyExoZoomIn;                //!<! Cluster M02 vs Energy vs exoticity, finer binning in exotic region
  TH3F *   fhM20EnergyExoM02MinCut;             //!<! Cluster M20 vs Energy vs exoticity for M02 > 0.1
  TH3F *   fhM02ExoNCells[fgkNEBins];           //!<! Cluster M02 vs exoticity vs n cells
  
  // Cluster-Track matching
  
  TH3F *   fhTrackMatchedDEtaNegExo;            //!<! Eta distance between - track and cluster vs cluster E vs exoticity, n cells > 1
  TH3F *   fhTrackMatchedDPhiNegExo;            //!<! Phi distance between - track and cluster vs cluster E vs exoticity, n cells > 1
  TH3F *   fhTrackMatchedDEtaDPhiNegExo;        //!<! Eta vs Phi distance between - track and cluster vs exoticity, E cluster > fEMinForExo and n cells > 1
  
  TH3F *   fhTrackMatchedDEtaPosExo;            //!<! Eta distance between + track and cluster vs cluster E, exoticity, n cells > 1
  TH3F *   fhTrackMatchedDPhiPosExo;            //!<! Phi distance between + track and cluster vs cluster E, exoticity, n cells > 1
  TH3F *   fhTrackMatchedDEtaDPhiPosExo;        //!<! Eta vs Phi distance between + track and cluster vs exoticity, E cluster > fEMinForExo and n cells > 1
  
  TH3F *   fhEOverPExo;                         //!<! E/p for track-cluster matches vs exoticity, n cells > 1

  TH2F *   fhTrackMatchedDEtaNeg1Cell;          //!<! Eta distance between - track and cluster vs cluster E, n cells = 1
  TH2F *   fhTrackMatchedDPhiNeg1Cell;          //!<! Phi distance between - track and cluster vs cluster E, n cells = 1
  TH2F *   fhTrackMatchedDEtaDPhiNeg1Cell;      //!<! Eta vs Phi distance between - track and cluster, E cluster > fEMinForExo and n cells = 1
  
  TH2F *   fhTrackMatchedDEtaPos1Cell;         //!<! Eta distance between + track and cluster vs cluster E, n cells = 1
  TH2F *   fhTrackMatchedDPhiPos1Cell;         //!<! Phi distance between + track and cluster vs cluster E, n cells = 1
  TH2F *   fhTrackMatchedDEtaDPhiPos1Cell;     //!<! Eta vs Phi distance between + track and cluster, E cluster > fEMinForExo and n cells = 1
  
  TH2F *   fhEOverP1Cell;                      //!<! E/p for track-cluster matches, n cells = 1
  
  // Calorimeter cells
    
  TH2F *   fhCellExoAmp;                     //!<! Cell amplitude vs exoticity
  TH3F *   fhCellExoAmpTime;                 //!<! Cell amplitude vs time vs exoticity
  TH3F *   fhCellExoGrid ;                   //!<! Cells ordered in column/row vs exoticity when amplitude > fEMinForExo 


  /// Copy constructor not implemented.
  AliAnaCaloExotics & operator = (const AliAnaCaloExotics & qa) ;
    
  /// Assignment operator not implemented.
  AliAnaCaloExotics(              const AliAnaCaloExotics & qa) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaCaloExotics,3) ;
  /// \endcond

} ;

#endif //ALIANACALORIMETERQA_H



