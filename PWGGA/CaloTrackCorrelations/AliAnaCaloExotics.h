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
 
  Float_t      GetHighNCellCut()    const   { return fNCellHighCut        ; }
  void         SetHighNCellCut(Int_t n)     { fNCellHighCut = n           ; }
  
  Double_t     GetTimeCutMin()      const   { return fTimeCutMin          ; }
  Double_t     GetTimeCutMax()      const   { return fTimeCutMax          ; }  
  void         SetTimeCut(Float_t min, Float_t max) { fTimeCutMin = min   ; 
                                                      fTimeCutMax = max   ; }
  
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
  
  Int_t    fNCellHighCut;                       ///<  High N cell with w>0  cut for some histograms 
  
  Float_t  fTimeCutMin  ;                       ///<  Remove clusters with time smaller than this value, in ns
  Float_t  fTimeCutMax  ;                       ///<  Remove clusters with time larger than this value, in ns
  
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
  TH2F *   fhExoticityEClusAllSameTCard;        //!<! Exoticity vs energy, all cells in same T-Card
  TH3F *   fhExoticityEClusPerSM;               //!<! Exoticity vs cluster energy, per SM
  TH2F *   fhExoticityEMaxCell;                 //!<! Exoticity vs energy of highest energy cell in cluster
  TH2F *   fhExoticityEClusTrackMatch;          //!<! Exoticity vs cluster energy, track matched
  TH2F *   fhExoticity1Cell;                    //!<! Exoticity vs energy for 1 cell clusters

  TH2F *   fhNCellsPerCluster;                  //!<! Cluster energy vs N cells in cluster
  TH2F *   fhNCellsPerClusterAllSameTCard;      //!<! Cluster energy vs N cells, all cells in same T-Card
  TH3F *   fhNCellsPerClusterPerSM;             //!<! Cluster energy vs N cells in cluster, per SM
  TH3F *   fhNCellsPerClusterWPerSM;            //!<! Cluster energy vs N cells in cluster, per SM
  TH3F *   fhNCellsPerClusterExo;               //!<! Cluster energy vs N cells in cluster vs Exoticity   
  TH3F *   fhNCellsPerClusterExoPerSM[20];      //!<! Cluster energy vs N cells in cluster vs Exoticity, per SM  
  TH2F *   fhNCellsPerClusterTrackMatch;        //!<! Cluster energy vs N cells in cluster, for track-matched clusters 
  TH3F *   fhNCellsPerClusterExoTrackMatch;     //!<! Cluster energy vs N cells in cluster vs Exoticity, for track-matched clusters 
  TH3F *   fhNCellsPerClusterM02;               //!<! Cluster energy vs N cells in cluster vs M02   

  TH3F *   fhEtaPhiGridExoEnCut  ;              //!<! column vs row vs exoticity when E > fEMinForExo and n cells > 1
  TH3F *   fhEtaPhiGridExoEnCutSameFracCut;     //!<! column vs row vs exoticity when E > fEMinForExo and n cells > 1 and n diff = 0
  TH3F *   fhEtaPhiGridEnExoCut  ;              //!<! column vs row vs energy when F+ > 0.97 and n cells > 1
  TH3F *   fhEtaPhiGridEn1Cell;                 //!<! column vs row vs energy for 1 cell clusters 
  TH3F *   fhEtaPhiGridEnHighNCells;            //!<! column vs row vs energy for n cell >  fNCellCut
  TH3F *   fhEtaPhiGridNCellEnCut;              //!<! column vs row vs n cells for E >  fEMinForExo
  
  TH3F *   fhTimeEnergyExo;                     //!<! Cluster Energy vs Time vs Exoticity, n cells > 1
  TH2F *   fhTimeEnergy1Cell;                   //!<! Cluster Energy vs Time vs n cells = 1
  TH3F *   fhTimeDiffClusCellExo;               //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs cluster energy vs exoticity
  TH3F *   fhTimeDiffAmpClusCellExo;            //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs secondary cell energy vs exoticity for E > fEMinForExo
  TH3F *   fhTimeEnergyM02;                     //!<! Cluster Energy vs Time vs M02, n cells > 1
  TH3F *   fhTimeDiffClusCellM02;               //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs cluster energy vs M02
  TH3F *   fhTimeEnergyNCells;                  //!<! Cluster Energy vs Time vs n cells
  TH3F *   fhTimeEnergyNCellsW;                 //!<! Cluster Energy vs Time vs n cells for w > 0
  TH1F *   fhTimeNCellCut;                      //!<! Cluster Time vs n cells > fNCellCut, larger time range 
  
  TH3F *   fhM02EnergyNCell;                    //!<! Cluster M02 vs Energy vs n cells
  TH2F *   fhM02EnergyAllSameTCard;             //!<! Cluster M02 vs Energy, all cells in same T-Card
  TH3F *   fhM02EnergyExo;                      //!<! Cluster M02 vs Energy vs exoticity
  TH3F *   fhM02EnergyExoZoomIn;                //!<! Cluster M02 vs Energy vs exoticity, finer binning in exotic region
  TH3F *   fhM20EnergyExoM02MinCut;             //!<! Cluster M20 vs Energy vs exoticity for M02 > 0.1
  TH3F *   fhM02ExoNCells[fgkNEBins];           //!<! Cluster M02 vs exoticity vs n cells, different E bins
  
  // Different n cells definitions
  TH2F *   fhNCellsPerClusterW ;                //!<! Cluster E vs n cells with w > 0         
  TH2F *   fhNCellsPerClusterSame;              //!<! Cluster E vs n cells in same T-Card as max E cell
  TH2F *   fhNCellsPerClusterDiff;              //!<! Cluster E vs n cells in different T-Card as max E cell 
  TH2F *   fhNCellsPerClusterSame5;             //!<! Cluster E vs n cells in same T-Card as max E cell, neighbour cells
  TH2F *   fhNCellsPerClusterDiff5;             //!<! Cluster E vs n cells in different T-Card as max E cell, neighbour cells 
  TH2F *   fhNCellsPerClusterSameW ;            //!<! Cluster E vs n cells with w > 0 in same T-Card as max E cell
  TH2F *   fhNCellsPerClusterDiffW ;            //!<! Cluster E vs n cells with w > 0 in same T-Card as max E cell
  TH3F *   fhNCellsPerClusterSameDiff;          //!<! Cluster E vs n cells in same vs different T-Card as max E cell
  TH2F *   fhNCellsPerClusterSameFrac;          //!<! Cluster E vs fraction n cells in same T-Card as max E cell
  TH3F *   fhNCellsPerClusterSameFracExo;       //!<! Cluster E vs fraction n cells in same T-Card as max E cell vs exoticity

  TH2F *   fhExoSame;                           //!<! Cluster E vs 1 - E same TCard / E max
  TH2F *   fhExoDiff;                           //!<! Cluster E vs 1 - E different TCard / E max
  TH2F *   fhExoSame5;                          //!<! Cluster E vs 1 - E same & neighbor TCard / E max
  TH2F *   fhExoDiff5;                          //!<! Cluster E vs 1 - E different & neighbor TCard / E max
  
  TH2F *   fhFracEnDiffSame;                    //!<! Cluster E vs fraction of energy in diff over same T-Card
  TH2F *   fhFracNCellDiffSame;                 //!<! Cluster E vs fraction of n cells in diff over same T-Card
  TH2F *   fhFracEnNCellDiffSame;               //!<! Cluster E vs double fraction of energy/n cells in diff over same T-Card
  
  TH2F *   fhFracEnDiffSameW;                   //!<! Cluster E vs fraction of energy in diff over same T-Card, cells with w > 0
  TH2F *   fhFracNCellDiffSameW;                //!<! Cluster E vs fraction of n cells in diff over same T-Card, cells with w > 0
  TH2F *   fhFracEnNCellDiffSameW;              //!<! Cluster E vs double fraction of energy/n cells in diff over same T-Card, cells with w > 0
  
  TH2F *   fhFracEnDiffSame5;                   //!<! Cluster E vs fraction of energy in diff over same T-Card, neighbour cells
  TH2F *   fhFracNCellDiffSame5;                //!<! Cluster E vs fraction of n cells in diff over same T-Card, neighbour cells
  TH2F *   fhFracEnNCellDiffSame5;              //!<! Cluster E vs double fraction of energy/n cells in diff over same T-Card, neighbour cells

  TH3F *   fhFracEnDiffSameExo;                 //!<! Cluster E vs fraction of energy in diff over same T-Card vs F+
  TH3F *   fhFracNCellDiffSameExo;              //!<! Cluster E vs fraction of n cells in diff over same T-Card vs F+
  TH3F *   fhFracEnNCellDiffSameExo;            //!<! Cluster E vs double fraction of energy/n cells in diff over same T-Card vs F+
  
  TH3F *   fhFracEnDiffSameWExo;                //!<! Cluster E vs fraction of energy in diff over same T-Card, cells with w > 0 vs F+
  TH3F *   fhFracNCellDiffSameWExo;             //!<! Cluster E vs fraction of n cells in diff over same T-Card, cells with w > 0 vs F+
  TH3F *   fhFracEnNCellDiffSameWExo;           //!<! Cluster E vs double fraction of energy/n cells in diff over same T-Card, cells with w > 0 vs F+
  
  TH3F *   fhFracEnDiffSame5Exo;                //!<! Cluster E vs fraction of energy in diff over same T-Card, neighbour cells vs F+
  TH3F *   fhFracNCellDiffSame5Exo;             //!<! Cluster E vs fraction of n cells in diff over same T-Card, neighbour cells vs F+
  TH3F *   fhFracEnNCellDiffSame5Exo;           //!<! Cluster E vs double fraction of energy/n cells in diff over same T-Card, neighbour cells vs F+
  
  TH1F *   fhFracEnDiffSameEnCut;               //!<! Cluster E > 20, fraction of energy in diff over same T-Card
  TH1F *   fhFracNCellDiffSameEnCut;            //!<! Cluster E > 20, fraction of n cells in diff over same T-Card
  TH1F *   fhFracEnNCellDiffSameEnCut;          //!<! Cluster E > 20, double fraction of energy/n cells in diff over same T-Card
  
  TH1F *   fhFracEnDiffSameWEnCut;              //!<! Cluster E > 20, fraction of energy in diff over same T-Card, cells with w > 0
  TH1F *   fhFracNCellDiffSameWEnCut;           //!<! Cluster E > 20, fraction of n cells in diff over same T-Card, cells with w > 0
  TH1F *   fhFracEnNCellDiffSameWEnCut;         //!<! Cluster E > 20, double fraction of energy/n cells in diff over same T-Card, cells with w > 0
  
  TH1F *   fhFracEnDiffSame5EnCut;              //!<! Cluster E > 20, fraction of energy in diff over same T-Card, neighbour cells
  TH1F *   fhFracNCellDiffSame5EnCut;           //!<! Cluster E > 20, fraction of n cells in diff over same T-Card, neighbour cells
  TH1F *   fhFracEnNCellDiffSame5EnCut;         //!<! Cluster E > 20, double fraction of energy/n cells in diff over same T-Card, neighbour cells
  
  TH3F *   fhNCellsSameDiffExo[fgkNEBins];      //!<! N cells in same vs diff vs exoticity, different cluster E bins
  TH3F *   fhEnSameDiffExo    [fgkNEBins];      //!<! Sum of E in same vs diff vs exoticity, different  cluster E bins
  TH3F *   fhCellEnSameExo;                     //!<! Cluster E vs cell E in same T-Card vs exoticity
  TH3F *   fhCellEnDiffExo;                     //!<! Cluster E vs cell E in diff T-Card vs exoticity
  TH3F *   fhCellEnNCellWOpenTime;              //!<! Cluster E vs cell E vs n cells with w > 0, no time cut
  TH3F *   fhCellEnNCellW;                      //!<! Cluster E vs cell E vs n cells with w > 0
  TH3F *   fhCellMaxClusterEnRatioNCellWOpenTime;//!<! Cluster E vs E cell max / E cluster vs n cells with w > 0, no time cut
  TH3F *   fhCellMaxClusterEnRatioNCellW;       //!<!  Cluster E vs E cell max / E cluster vs n cells with w > 0
  TH3F *   fhCellMaxClusterEnRatioExo;          //!<!  Cluster E vs E cell max / E cluster vs n exoticity

  
  // Cluster column-row
  TH3F *   fhClusterColRowExo[2][fgkNEBins];     //!<! Cluster col-row centred in cell max vs exoticity for different cluster E bins
  //TH2F *   fhClusterColRow   [2][fgkNEBins];     //!<! Cluster col-row centred in cell max for different cluster E bins
  //TH3F *   fhClusterColRowExoW[2][fgkNEBins];    //!<! Cluster col-row centred in cell max vs exoticity and w>0 for different cluster E bins 

  TH3F *   fhClusterColRowPerSMHighNCell[fgkNEBins]; //!<! Cluster col-row centred in cell max vs SM for different cluster E bins and n cells > fHighNCellCut

  
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
  
  TH2F *   fhTrackMatchedDEtaPos1Cell;          //!<! Eta distance between + track and cluster vs cluster E, n cells = 1
  TH2F *   fhTrackMatchedDPhiPos1Cell;          //!<! Phi distance between + track and cluster vs cluster E, n cells = 1
  TH2F *   fhTrackMatchedDEtaDPhiPos1Cell;      //!<! Eta vs Phi distance between + track and cluster, E cluster > fEMinForExo and n cells = 1
  
  TH2F *   fhEOverP1Cell;                       //!<! E/p for track-cluster matches, n cells = 1
  
  // Calorimeter cells
    
  TH2F *   fhCellExoAmp;                        //!<! Cell amplitude vs exoticity
  TH3F *   fhCellExoAmpTime;                    //!<! Cell amplitude vs time vs exoticity
  TH3F *   fhCellExoGrid ;                      //!<! Cells ordered in column/row vs exoticity when amplitude > fEMinForExo 


  /// Copy constructor not implemented.
  AliAnaCaloExotics & operator = (const AliAnaCaloExotics & qa) ;
    
  /// Assignment operator not implemented.
  AliAnaCaloExotics(              const AliAnaCaloExotics & qa) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaCaloExotics,5) ;
  /// \endcond

} ;

#endif //ALIANACALORIMETERQA_H



