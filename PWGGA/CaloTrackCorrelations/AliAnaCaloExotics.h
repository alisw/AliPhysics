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
  
  void         StripHistograms(AliVCaloCells *cells);

  void         ClusterHistogramsAfterEventCut(const TObjArray * clusters, AliVCaloCells * cells);
  
  void         TrackHistogramsAfterEventCut();
  
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
  
  void         SetHighEnergyCutSM(Float_t e) { fHighEnergyCutSM = e ; }   
  void         SetHighNCellsCutSM(Int_t   n) { fHighNCellsCutSM = n ; }   

  void         SetLowEnergyCutSM3(Float_t e) { fLowEnergyCutSM3 = e ; }   
  void         SetLowNCellsCutSM3(Int_t   n) { fLowNCellsCutSM3 = n ; }   
  
  void         SetEventMaxNumberOfStrips(Int_t n) { fEventMaxNumberOfStrips = n ; }   
  
  void         SetfHighEnergyCutStrip(Float_t eFull, Float_t eThird)
                { fHighEnergyCutStrip[0] = eFull; fHighEnergyCutStrip[1] = eThird; }
 
  void         SetfHighNCellsCutStrip(Int_t nFull, Int_t nThird)
                { fHighNCellsCutStrip[0] = nFull; fHighNCellsCutStrip[1] = nThird; }
  
  void         SetLowEnergyCutSM3Strip(Float_t e) { fLowEnergyCutSM3Strip = e ; }
  void         SetLowNCellsCutSM3Strip(Int_t   n) { fLowNCellsCutSM3Strip = n ; }
  
  Float_t      GetEBinLimit(Int_t i) const  
                           { if ( i < fgkNEBins && i >= 0 ) return fEnergyBins[i] ;
                             else                           return -1     ; }
  void         SetEBinLimit(Int_t i, Float_t en) 
                           { if ( i < fgkNEBins && i >= 0 ) fEnergyBins[i] = en ; }

  Float_t      GetCellEMinBinLimit(Int_t i) const  
                           { if ( i < fgkNCellEnMinBins && i >= 0 ) return fCellEnMins[i] ;
                             else                                return -1     ; }
  void         SetCellEMinBinLimit(Int_t i, Float_t en) 
                           { if ( i < fgkNCellEnMinBins && i >= 0 ) fCellEnMins[i] = en ; }
  
  Float_t      GetCellEMaxBinLimit()                { return fCellEnMax ; } 
  void         SetCellEMaxBinLimit(Float_t en)      { fCellEnMax = en   ; } 
  
  void         SwitchOnFill1CellHisto()             { fFill1CellHisto            = kTRUE  ; }
  void         SwitchOffFill1CellHisto()            { fFill1CellHisto            = kFALSE ; }
 
  void         SwitchOnFillCellHisto()              { fFillCellHisto             = kTRUE  ; }
  void         SwitchOffFillCellHisto()             { fFillCellHisto             = kFALSE ; }

  void         SwitchOnFillAllCellEventHisto( Int_t opt = 1 ) { fFillAllCellEventParamHisto = opt  ; }
  void         SwitchOffFillAllCellEventHisto()     { fFillAllCellEventParamHisto= 0 ; }
  
  void         SwitchOnFillStripHisto()             { fFillStripHisto            = kTRUE  ; }
  void         SwitchOffFillStripHisto()            { fFillStripHisto            = kFALSE ; }
  
  void         SwitchOnFillMatchingHisto()          { fFillMatchingHisto         = kTRUE  ; }
  void         SwitchOffFillMatchingHisto()         { fFillMatchingHisto         = kFALSE ; }
 
  void         SwitchOnFillFractionSameDiffCells()  { fFillSameDiffFracHisto     = kTRUE  ; }
  void         SwitchOffFillFractionSameDiffCells() { fFillSameDiffFracHisto     = kFALSE ; }
  
  void         SwitchOnFillExoEnMinCutHisto()       { fFillExoEnMinCut           = kTRUE  ; }
  void         SwitchOffFillExoEnMinCutHisto()      { fFillExoEnMinCut           = kFALSE ; }
  
  void         SwitchOnFillAllSameTCardHisto()      { fFillAllCellSameTCardHisto = kTRUE  ; }
  void         SwitchOffFillAllSameTCardHisto()     { fFillAllCellSameTCardHisto = kFALSE ; }
  
  void         SwitchOnFillPerSMHisto()             { fFillPerSMHisto            = kTRUE  ; }
  void         SwitchOffFillPerSMHisto()            { fFillPerSMHisto            = kFALSE ; }
  
  void         SwitchOnClusterColRowHisto()         { fFillClusterColRowHisto    = kTRUE  ; }
  void         SwitchOffFillClusterColRowHisto()    { fFillClusterColRowHisto    = kFALSE ; }
  
  void         SwitchOnFillOpenTimeHisto()          { fFillOpenTimeHisto         = kTRUE  ; }
  void         SwitchOffFillOpenTimeHisto()         { fFillOpenTimeHisto         = kFALSE ; }
 
  void         SwitchOnFillExoticity50nsHisto()     { fFillExo50ns               = kTRUE  ; }
  void         SwitchOffFillExoticity50nsHisto()    { fFillExo50ns               = kFALSE ; }
  //
  void         SwitchOnFillClusterAfterEventCutHisto() { fFillClusterHistoAfterEventCut = kTRUE  ; }
  void         SwitchOffFillClusterAfterEventCutHisto(){ fFillClusterHistoAfterEventCut = kFALSE ; }

  void         SwitchOnFillTrackAfterEventCutHisto()   { fFillTrackHistoAfterEventCut   = kTRUE  ; }
  void         SwitchOffFillTrackAfterEventCutHisto()  { fFillTrackHistoAfterEventCut   = kFALSE ; }
  
  void         SetConstantTimeShift(Float_t shift)  { fConstantTimeShift = shift  ; }
  
 private:
      
 
  // Different analysis cuts
  //
  Float_t  fCellAmpMin;                         ///<  Amplitude Threshold on calorimeter cells 
  
  Float_t  fEMinForExo;                         ///<  Minimum energy cut for integrated histograms over E. Exotics start at about 4 GeV.

  Float_t  fExoCut;                             ///<  Exoticity cut for some histograms 
  
  Int_t    fNCellHighCut;                       ///<  High N cell with w>0  cut for some histograms 
  
  Float_t  fTimeCutMin  ;                       ///<  Remove clusters with time smaller than this value, in ns
  Float_t  fTimeCutMax  ;                       ///<  Remove clusters with time larger than this value, in ns
  
  Float_t  fHighEnergyCutSM;                    ///<  SM is too active if energy above this value
  Int_t    fHighNCellsCutSM;                    ///<  SM is too active if n cells above this value
  Float_t  fLowEnergyCutSM3;                    ///<  SM3 low activity if energy above this value
  Int_t    fLowNCellsCutSM3;                    ///<  SM3 low activity if n cells above this value
  
  Int_t    fEventMaxNumberOfStrips;             ///<  Cut on events on number of too active strips
  Float_t  fHighEnergyCutStrip[2];              ///<  Strip is too active if energy above this value
  Int_t    fHighNCellsCutStrip[2];              ///<  Strip is too active if n cells above this value
  Float_t  fLowEnergyCutSM3Strip;               ///<  SM3 strips low activity if energy above this value
  Int_t    fLowNCellsCutSM3Strip;               ///<  SM3 strips low activity if n cells above this value
  
  /// Total number of cluster energy bins histograms
  static const Int_t fgkNEBins = 14;
  Float_t  fEnergyBins[fgkNEBins];              ///<  Energy bins for some histograms
 
  static const Int_t fgkNCellEnMinBins   = 3;
  static const Int_t fgkNCellEnMinBinsFr = 3;
  Float_t  fCellEnMins[fgkNCellEnMinBins];      ///<  Energy minimum for some histograms for cells
  Float_t  fCellEnMax;                          ///<  Energy maximum for some histograms for cells 
  
  Float_t  fConstantTimeShift;                  ///<  Apply a 600 ns time shift in case of simulation, shift in ns.

  // Internal values containers
  //  
  TLorentzVector fClusterMomentum;              //!<! Cluster momentum, temporary container
  
  Bool_t   fLED20       ;                       ///<  There is at least a cluster with 20 cells with w > 0 in the event; internal
  Bool_t   fLED12       ;                       ///<  There is at least a cluster with 12 cells with w > 0 in the event; internal
  Float_t  fLED20Time   ;                       ///<  Time of LED cluster with 20 cells 
  Float_t  fLED12Time   ;                       ///<  Time of LED cluster with 12 cells
  
  Float_t  fAcceptEvent ;                       ///<  Event passes the SM n cells or energy cuts
  
  Int_t    fEventNStripActive ;                 ///<  Number of strips per event with multiple high energy cells
  Int_t    fEventNStripActiveSM[20];            ///<  Number of strips per event with multiple high energy cells per SM
  
  Float_t  fEnCellsStrip[fgkNCellEnMinBins][20][24]; ///< Number of cells in strip per event, per SM and per cell min En cut
  Int_t    fnCellsStrip [fgkNCellEnMinBins][20][24]; ///< Sum of energy of cells in strip per event, per SM and per cell min En cut
    
  // Histogram/analysis switchs
  //
  Bool_t   fFillCellHisto;                      ///<  Fill histograms single cells
  
  Int_t    fFillAllCellEventParamHisto;         ///<  Fill histograms summing our counting cells per event, 0-not filled, 1-n and sum E, 2 all fractions
   
  Bool_t   fFill1CellHisto;                     ///<  Fill histograms for 1 cell clusters
  
  Bool_t   fFillStripHisto;                     ///<  Fill histograms for strip of cells (2x48)
  
  Bool_t   fFillMatchingHisto;                  ///<  Fill histograms for track-cluster matching
  
  Bool_t   fFillSameDiffFracHisto;              ///<  Fill histograms with fraction of cells or energy cells in same or different T-Card
  
  Bool_t   fFillExoEnMinCut;                    ///<  Fill histograms with exoticity for different EMin cuts
  
  Bool_t   fFillAllCellSameTCardHisto;          ///<  Fill histograms when all cells are in same  T-Card
  
  Bool_t   fFillPerSMHisto;                     ///<  Fill histograms per SM
  
  Bool_t   fFillClusterColRowHisto;             ///<  Fill histograms with cluster coll and row size

  Bool_t   fFillOpenTimeHisto;                  ///<  Fill histograms when time cut not applied (not needed in MC)
  
  Bool_t   fFillExo50ns;                        ///<  Fill histograms  with exoticity with 50 ns cut on difference in time with max
  
  Bool_t   fFillClusterHistoAfterEventCut;      ///<  Fill histograms for clusters after event cuts
  
  Bool_t   fFillTrackHistoAfterEventCut;        ///<  Fill histograms for tracks after event cuts
  
  // Histograms
  //
  TH1F *   fhNClusterPerEventNCellHigh20;        //!<! N clusters with NCells > 20 per event
  TH1F *   fhNClusterPerEventNCellHigh12;        //!<! N clusters with NCells > 12 per event
  
  TH1F *   fhNClusterPerEventExotic;             //!<! N clusters with F+ > 0.97  and E > 5 GeV per event
  TH1F *   fhNClusterPerEventExotic1Cell;        //!<! N clusters with ncell = 1  and E > 5 GeV per event
  TH1F *   fhNClusterPerEventExoticNCell;        //!<! N clusters with F+ > 0.97 and ncell > 1 and E > 5 GeV per event
  
  TH2F *   fh2NClusterPerEventNCellHigh20;        //!<! N clusters total vs with NCells > 20 per event
  TH2F *   fh2NClusterPerEventNCellHigh12;        //!<! N clusters total vs with NCells > 12 per event
  
  TH2F *   fh2NClusterPerEventExotic;             //!<! N clusters total vs with F+ > 0.97  and E > 5 GeV per event
  TH2F *   fh2NClusterPerEventExotic1Cell;        //!<! N clusters total vs with ncell = 1  and E > 5 GeV per event
  TH2F *   fh2NClusterPerEventExoticNCell;        //!<! N clusters total vs with F+ > 0.97 and ncell > 1 and E > 5 GeV per event
  TH2F *   fh2NClusterPerEventExoticAmpMax;       //!<! Highest energy exotic amp max vs n clusters
  
  // Calorimeter Clusters
  //

  // Exoticity
  //
  TH2F *   fhExoticityEClus;                    //!<! Exoticity vs cluster energy
  TH3F *   fhExoticityEClusPerSM;               //!<! Exoticity vs cluster energy, per SM
  TH2F *   fhExoticityEMaxCell;                 //!<! Exoticity vs energy of highest energy cell in cluster
  TH2F *   fhExoticityEClusTrackMatch;          //!<! Exoticity vs cluster energy, track matched
  TH2F *   fhExoticity1Cell;                    //!<! Exoticity vs energy for 1 cell clusters
  
  TH2F *   fhNCellsPerCluster;                  //!<! Cluster energy vs N cells in cluster
  TH2F *   fhNCellsPerClusterW;                 //!<! Cluster energy vs N cells in cluster with w > 0
  TH2F *   fhNCellsPerClusterTimeDiff;          //!<! Cluster energy vs N cells in cluster with t diff < 50 ns
  TH2F *   fhNCellsPerClusterEMaxCell;          //!<! Cell max energy vs N cells in cluster
  TH2F *   fhNCellsPerClusterWEMaxCell;         //!<! Cell max energy vs N cells in cluster with w > 0
  TH2F *   fhNCellsPerClusterOpenTime;          //!<! Cluster energy vs N cells in cluster, no time cut
  TH2F *   fhNCellsPerClusterWOpenTime;         //!<! Cluster energy vs N cells in cluster with w > 0, no time cut
  TH2F *   fhNCellsPerClusterEMaxCellOpenTime;  //!<! Cell max energy vs N cells in cluster, no time cut
  TH2F *   fhNCellsPerClusterWEMaxCellOpenTime; //!<! Cell max energy vs N cells in cluster with w > 0, no time cut

  TH3F *   fhNCellsPerClusterExo;               //!<! Cluster energy vs N cells in cluster vs Exoticity   
  TH3F *   fhNCellsPerClusterPerSM;             //!<! Cluster energy vs N cells in cluster, per SM
  TH3F *   fhNCellsPerClusterWPerSM;            //!<! Cluster energy vs N cells in cluster, per SM
  TH3F *   fhNCellsPerClusterExoPerSM[20];      //!<! Cluster energy vs N cells in cluster vs Exoticity, per SM  
  TH2F *   fhNCellsPerClusterTrackMatch;        //!<! Cluster energy vs N cells in cluster, for track-matched clusters 
  TH3F *   fhNCellsPerClusterExoTrackMatch;     //!<! Cluster energy vs N cells in cluster vs Exoticity, for track-matched clusters 
  TH3F *   fhNCellsPerClusterM02;               //!<! Cluster energy vs N cells in cluster vs M02   

  TH3F *   fhEtaPhiGridExoEnCut  ;              //!<! column vs row vs exoticity when E > fEMinForExo and n cells > 1
  TH3F *   fhEtaPhiGridEnExoCut  ;              //!<! column vs row vs energy when F+ > 0.97 and n cells > 1
  TH3F *   fhEtaPhiGridEn1Cell;                 //!<! column vs row vs energy for 1 cell clusters 
  TH3F *   fhEtaPhiGridEnHighNCells;            //!<! column vs row vs energy for n cell >  fNCellCut
  TH3F *   fhEtaPhiGridNCellEnCut;              //!<! column vs row vs n cells for E >  fEMinForExo
  
  TH3F *   fhTimeEnergyExo;                     //!<! Cluster Energy vs Time vs Exoticity, n cells > 1
  TH2F *   fhTimeEnergy1Cell;                   //!<! Cluster Energy vs Time vs n cells = 1
  TH3F *   fhTimeDiffClusCellExo;               //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs cluster energy vs exoticity
  TH3F *   fhTimeDiffClusCellDiffTCardExo;      //!<! Difference of the time of cell with maximum dep energy and the rest of cells in different T-Card vs cluster energy vs exoticity
  TH3F *   fhTimeDiffClusCellSameTCardExo;      //!<! Difference of the time of cell with maximum dep energy and the rest of cells in same T-Card vs cluster energy vs exoticity
  
  TH2F *   fhTimeDiffClusSameCellColRowDiff[4][4];       //!<! Cluster E vs time max-sec in same T-Card vs distance to max E cell, for exoticity < fExoCut
  TH2F *   fhTimeDiffClusDiffCellColRowDiff[4][4];       //!<! Cluster E vs time max-sec in diff T-Card vs distance to max E cell, for exoticity < fExoCut
  TH2F *   fhTimeDiffClusSameCellColRowDiffExoCut[4][4]; //!<! Cluster E vs time max-sec in same T-Card vs distance to max E cell, for exoticity > fExoCut
  TH2F *   fhTimeDiffClusDiffCellColRowDiffExoCut[4][4]; //!<! Cluster E vs time max-sec in diff T-Card vs distance to max E cell, for exoticity > fExoCut
  
  TH3F *   fhTimeDiffWClusCellExo;              //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs cluster energy vs exoticity, for cells with weight
  TH3F *   fhTimeDiffAmpClusCellExo;            //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs secondary cell energy vs exoticity for E > fEMinForExo
  TH3F *   fhTimeEnergyM02;                     //!<! Cluster Energy vs Time vs M02, n cells > 1
  TH3F *   fhTimeDiffClusCellM02;               //!<! Difference of the time of cell with maximum dep energy and the rest of cells vs cluster energy vs M02
  TH3F *   fhTimeEnergyNCells;                  //!<! Cluster Energy vs Time vs n cells
  TH3F *   fhTimeEnergyNCellsW;                 //!<! Cluster Energy vs Time vs n cells for w > 0
  TH1F *   fhTimeNCellCut;                      //!<! Cluster Time vs n cells > fNCellCut, larger time range 
  
  TH3F *   fhM02EnergyExo;                      //!<! Cluster M02 vs Energy vs exoticity
  TH3F *   fhM20EnergyExoM02MinCut;             //!<! Cluster M20 vs Energy vs exoticity for M02 > 0.1
  TH3F *   fhM02ExoNCells[fgkNEBins];           //!<! Cluster M02 vs exoticity vs n cells, different E bins
  
  // Different n cells definitions
  //
  TH2F *   fhNCellsPerClusterSame;              //!<! Cluster E vs n cells in same T-Card as max E cell
  TH2F *   fhNCellsPerClusterDiff;              //!<! Cluster E vs n cells in different T-Card as max E cell 
  TH2F *   fhNCellsPerClusterSame5;             //!<! Cluster E vs n cells in same T-Card as max E cell, neighbour cells
  TH2F *   fhNCellsPerClusterDiff5;             //!<! Cluster E vs n cells in different T-Card as max E cell, neighbour cells 
  TH2F *   fhNCellsPerClusterSameW ;            //!<! Cluster E vs n cells with w > 0 in same T-Card as max E cell
  TH2F *   fhNCellsPerClusterDiffW ;            //!<! Cluster E vs n cells with w > 0 in diff T-Card as max E cell
  TH2F *   fhNCellsPerClusterSameTimeDiff;      //!<! Cluster E vs n cells with t diff < 50 ns in same T-Card as max E cell
  TH2F *   fhNCellsPerClusterDiffTimeDiff;      //!<! Cluster E vs n cells with t diff < 50 ns in diff T-Card as max E cell
  TH3F *   fhNCellsPerClusterSameDiff;          //!<! Cluster E vs n cells in same vs different T-Card as max E cell
  TH3F *   fhNCellsPerClusterSameDiffW;         //!<! Cluster E vs n cells in same vs different T-Card as max E cell, cells with weight
  TH3F *   fhNCellsPerClusterSameDiffTimeDiff;  //!<! Cluster E vs n cells in same vs different T-Card as max E cell, cells with t diff < 50 ns
  TH2F *   fhNCellsPerClusterSameFrac;          //!<! Cluster E vs fraction n cells in same T-Card as max E cell
  TH3F *   fhNCellsPerClusterSameFracExo;       //!<! Cluster E vs fraction n cells in same T-Card as max E cell vs exoticity
  TH2F *   fhNCellsPerClusterSameFracW;         //!<! Cluster E vs fraction n cells with weight in same T-Card as max E cell
  TH3F *   fhNCellsPerClusterSameFracWExo;      //!<! Cluster E vs fraction n cells with weight in same T-Card as max E cell vs exoticity
  
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
  TH3F *   fhEnNCellsSameDiffExo[fgkNEBins];    //!<! Sum of E in / n cell in same vs diff vs exoticity, different  cluster E bins
  TH3F *   fhCellEnSameExo;                     //!<! Cluster E vs cell E in same T-Card vs exoticity
  TH3F *   fhCellEnDiffExo;                     //!<! Cluster E vs cell E in diff T-Card vs exoticity
  TH2F *   fhCellEnSameColRowDiff[4][4];        //!<! Cluster E vs cell E in same T-Card vs distance to max E cell, for exoticity < fExoCut
  TH2F *   fhCellEnDiffColRowDiff[4][4];        //!<! Cluster E vs cell E in diff T-Card vs distance to max E cell, for exoticity < fExoCut
  TH2F *   fhCellEnSameColRowDiffExoCut[4][4];  //!<! Cluster E vs cell E in same T-Card vs distance to max E cell, for exoticity > fExoCut
  TH2F *   fhCellEnDiffColRowDiffExoCut[4][4];  //!<! Cluster E vs cell E in diff T-Card vs distance to max E cell, for exoticity > fExoCut
  TH3F *   fhCellEnNCellWOpenTime;              //!<! Cluster E vs cell E vs n cells with w > 0, no time cut
  TH3F *   fhCellEnNCellW;                      //!<! Cluster E vs cell E vs n cells with w > 0
  TH3F *   fhCellEnNCellWEMaxOpenTime;          //!<! Cluster cell max E vs cell E vs n cells with w > 0, no time cut
  TH3F *   fhCellEnNCellWEMax;                  //!<! Cluster cell max E vs cell E vs n cells with w > 0
  TH3F *   fhCellTimeDiffNCellWOpenTime;        //!<! Cluster E vs cell time-time max vs n cells with w > 0, no time cut
  TH3F *   fhCellTimeDiffNCellW;                //!<! Cluster E vs cell time-time max vs n cells with w > 0
  TH3F *   fhCellTimeDiffNCellWEMaxOpenTime;    //!<! Cluster cell max E vs cell time-time max vs n cells with w > 0, no time cut
  TH3F *   fhCellTimeDiffNCellWEMax;            //!<! Cluster cell max E vs cell time-time max vs n cells with w > 0
  
  TH2F *   fhCellMaxClusterEnOpenTime;          //!<!  Cluster E (without non lin) vs Cell max E, no time cut 
  TH2F *   fhCellMaxClusterEn;                  //!<!  Cluster E (without non lin) vs Cell max E  
  TH2F *   fhCellMaxClusterEnRatioOpenTime;     //!<!  Cluster E vs E cell max / E cluster, no time cut 
  TH2F *   fhCellMaxClusterEnRatio;             //!<!  Cluster E vs E cell max / E cluster 
  TH3F *   fhCellMaxClusterEnRatioNCellWOpenTime;//!<! Cluster E vs E cell max / E cluster vs n cells with w > 0, no time cut
  TH3F *   fhCellMaxClusterEnRatioNCellW;       //!<!  Cluster E vs E cell max / E cluster vs n cells with w > 0
  TH3F *   fhCellMaxClusterEnRatioExo;          //!<!  Cluster E vs E cell max / E cluster vs n exoticity

  // Cluster column-row
  //
  TH3F *   fhClusterColRowExo[2][fgkNEBins];     //!<! Cluster col-row centred in cell max vs exoticity for different cluster E bins
  //TH2F *   fhClusterColRow   [2][fgkNEBins];     //!<! Cluster col-row centred in cell max for different cluster E bins
  //TH3F *   fhClusterColRowExoW[2][fgkNEBins];    //!<! Cluster col-row centred in cell max vs exoticity and w>0 for different cluster E bins 

  TH3F *   fhClusterColRowPerSMHighNCell[fgkNEBins]; //!<! Cluster col-row centred in cell max vs SM for different cluster E bins and n cells > fHighNCellCut

  // Recalculate exiticity for different EMin thresholds, only when fFillExoEnMinCut=true
  //
  TH2F *   fhExoticityWEClus;                   //!<! Weighted exoticity vs Min E cell threshold vs cluster energy
  TH3F *   fhNCellsPerClusterExoW;              //!<! Cluster energy vs N cells in cluster vs Weighted Exoticity   
  TH3F *   fhTimeEnergyExoW;                    //!<! Cluster Energy vs Time vs Weighted Exoticity, n cells > 1
  TH3F *   fhM02EnergyExoW;                     //!<! Cluster M02 vs Energy vs Weighted exoticity
  TH3F *   fhExoticityECellMinCut;              //!<! Exoticity vs Min E cell threshold vs cluster energy
    
  // Apply cut on number of cells of cluster all in same TCard
  //
  TH3F *   fhNCellsPerClusterMinEnCut;         //!<! Cluster energy vs N cells vs E min in cell  
  TH3F *   fhNCellsPerClusterDiffMinEnCut;     //!<! Cluster energy vs N cells in different T-Card as max cell vs E min in cell 
  TH3F *   fhNCellsPerClusterSameMinEnCut;     //!<! Cluster energy vs N cells in same TCard as max cell vs E min in cell 
  
  TH3F *   fhM02ExoNCellsNotAllSameTCard[fgkNEBins]; //!<! Cluster M02 vs exoticity vs n cells, different E bins for n diff > 0

  TH2F *   fhExoticityEClusAllSameTCard;        //!<! Exoticity vs energy, all cells in same T-Card, n diff = 0
  TH2F *   fhM02EnergyAllSameTCard;             //!<! Cluster M02 vs Energy, all cells in same T-Card, n diff = 0
  TH2F *   fhNCellsPerClusterAllSameTCard;      //!<! Cluster energy vs N cells, all cells in same T-Card, n diff = 0  
  TH3F *   fhEtaPhiGridExoEnCutSameFracCut;     //!<! column vs row vs exoticity when E > fEMinForExo and n cells > 1 and n diff = 0

  TH2F *   fhExoticityEClusAllSameTCardW;       //!<! Exoticity vs energy, all cells in same T-Card, n diff_w =0
  TH2F *   fhM02EnergyAllSameTCardW;            //!<! Cluster M02 vs Energy, all cells in same T-Card, n diff_w =0
  TH2F *   fhNCellsPerClusterAllSameTCardW;     //!<! Cluster energy vs N cells, all cells in same T-Card, n diff_w = 0  
  TH3F *   fhEtaPhiGridExoEnCutSameFracCutW;    //!<! column vs row vs exoticity when E > fEMinForExo and n cells > 1 and n diff_w = 0
 
  TH3F *   fhExoticityEClusAllSameTCardMinEnCut;   //!<! Exoticity vs energy, all cells in same T-Card, n diff =0 for E min cut
  TH3F *   fhM02EnergyAllSameTCardMinEnCut;        //!<! Cluster M02 vs Energy, all cells in same T-Card, n diff =0 for E min cut
  TH3F *   fhNCellsPerClusterAllSameTCardMinEnCut; //!<! Cluster energy vs N cells, all cells in same T-Card, n diff = 0 for E min cut
  
  // Apply 50 ns cut on exoticity
  TH2F *   fhExoticity50nsEClus;                //!<! Exoticity vs cluster energy, tdiff < 50 ns
  TH2F *   fhExoticity50ns1Cell;                //!<! Exoticity vs energy for 1 cell clusters, tdiff < 50 ns
  TH2F *   fhExoticity50nsEClusAllSameTCard;    //!<! Exoticity vs energy, all cells in same T-Card, n diff = 0, tdiff < 50 ns
  TH3F *   fhNCellsPerClusterExo50ns;           //!<! Cluster energy vs N cells in cluster vs Exoticity, tdiff < 50 ns   
  TH3F *   fhM02EnergyExo50ns;                  //!<! Cluster M02 vs Energy vs exoticity, tdiff < 50 ns
  TH3F *   fhM02Exo50nsNCells[fgkNEBins];       //!<! Cluster M02 vs exoticity vs n cells, different E bins, tdiff < 50 ns
  
  // Cluster-Track matching
  //
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
  
  //  
  // Calorimeter cells
    
  TH2F *   fhCellExoAmp;                        //!<! Cell amplitude vs exoticity
  TH3F *   fhCellExoAmpTime;                    //!<! Cell amplitude vs time vs exoticity
  TH2F *   fhCellExo50nsAmp;                    //!<! Cell amplitude vs exoticity, tdiff < 50 ns
  
  TH3F *   fhCellExoGrid ;                      //!<! Cells ordered in column/row vs exoticity when amplitude > fEMinForExo 
  TH3F *   fhCellGridTimeHighNCell20 ;          //!<! Cells ordered in column/row vs cluster time when at least 1 cluster n cellsW > 20
  TH3F *   fhCellGridTimeHighNCell12 ;          //!<! Cells ordered in column/row vs cluster time when at least 1 cluster n cellsW > 12

  // En cuts: 0.5, 1, 2, 5 GeV
  TH3F *   fhCellGridTime         [fgkNCellEnMinBins];        //!<! Cells ordered in column/row vs cell time when E > emin
  TH1F *   fhSumEnCells           [fgkNCellEnMinBins];        //!<! For E cell > emin, sum of cells energy 
  TH1F *   fhNCells               [fgkNCellEnMinBins];        //!<! For E cell > emin, count number of cells 
  TH1F *   fhAverSumEnCells       [fgkNCellEnMinBins];        //!<! For E cell > emin, sum of cells energy / total cells number
  TH1F *   fhFracNCells           [fgkNCellEnMinBinsFr];      //!<! total number of cells with E > emin over 0.5 GeV
  TH1F *   fhFracSumEnCells       [fgkNCellEnMinBinsFr];      //!<! sum of cells with E > emin over 0.5 GeV
  
  TH1F *   fhSumEnCellsNHigh20    [fgkNCellEnMinBins];        //!<! For E cell > emin, sum of cells energy, 1 cluster with n_cell_w>20 
  TH1F *   fhNCellsNHigh20        [fgkNCellEnMinBins];        //!<! For E cell > emin, count number of cells, 1 cluster with n_cell_w>20 
  TH1F *   fhAverSumEnCellsNHigh20[fgkNCellEnMinBins];        //!<! For E cell > emin, sum of cells energy / total cells number, 1 cluster with n_cell_w>20
  TH1F *   fhFracNCellsNHigh20    [fgkNCellEnMinBinsFr];      //!<! total number of cells with E > emin over 0.5 GeV, 1 cluster with n_cell_w>20
  TH1F *   fhFracSumEnCellsNHigh20[fgkNCellEnMinBinsFr];      //!<! sum of cells with E > 1 GeV over emin, 1 cluster with n_cell_w>20
 
  TH1F *   fhSumEnCellsAcceptEvent     [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy, LED rejected 
  TH1F *   fhNCellsAcceptEvent         [fgkNCellEnMinBins];   //!<! For E cell > emin, count number of cells, LED rejected
  TH1F *   fhAverSumEnCellsAcceptEvent [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy / total cells number, LED rejected
  TH1F *   fhFracNCellsAcceptEvent     [fgkNCellEnMinBinsFr]; //!<! total number of cells with E > emin over 0.5 GeV, LED rejected
  TH1F *   fhFracSumEnCellsAcceptEvent [fgkNCellEnMinBinsFr]; //!<! sum of cells with E > 1 GeV over emin, LED rejected
 
  TH1F *   fhSumEnCellsAcceptEventStrip[fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy, LED rejected, based on n strip 
  TH1F *   fhNCellsAcceptEventStrip    [fgkNCellEnMinBins];   //!<! For E cell > emin, count number of cells, LED rejected, based on n strip
  TH1F *   fhSumEnCellsAcceptEventBoth[fgkNCellEnMinBins];    //!<! For E cell > emin, sum of cells energy, LED rejected, based on n strip and SM activity
  TH1F *   fhNCellsAcceptEventBoth    [fgkNCellEnMinBins];    //!<! For E cell > emin, count number of cells, LED rejected, based on n strip and SM activity
  
  TH2F *   fhNSumEnCells                [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy 
  TH2F *   fhNSumEnCellsAcceptEvent     [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy, LED rejected 
  TH2F *   fhNSumEnCellsAcceptEventStrip[fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy, LED rejected, based on n strip 
  TH2F *   fhNSumEnCellsAcceptEventBoth [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy, LED rejected, based on n strip and SM activity 

  // Per SM
  TH2F *   fhSumEnCellsPerSM           [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy, per SM 
  TH2F *   fhNCellsPerSM               [fgkNCellEnMinBins];   //!<! For E cell > emin, count number of cells, per SM 
  TH2F *   fhAverSumEnCellsPerSM       [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy / total cells number, per SM
  TH2F *   fhFracNCellsPerSM           [fgkNCellEnMinBinsFr]; //!<! total number of cells with E > emin over 0.5 GeV, per SM
  TH2F *   fhFracSumEnCellsPerSM       [fgkNCellEnMinBinsFr]; //!<! sum of cells with E > emin over 0.5 GeV, per SM
  
  TH2F *   fhSumEnCellsPerSMNHigh20    [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy, 1 cluster with n_cell_w>20, per SM 
  TH2F *   fhNCellsPerSMNHigh20        [fgkNCellEnMinBins];   //!<! For E cell > emin, count number of cells, 1 cluster with n_cell_w>20, per SM 
  TH2F *   fhAverSumEnCellsPerSMNHigh20[fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy / total cells number, 1 cluster with n_cell_w>20, per SM
  TH2F *   fhFracNCellsPerSMNHigh20    [fgkNCellEnMinBinsFr]; //!<! total number of cells with E > emin over 0.5 GeV, 1 cluster with n_cell_w>20, per SM
  TH2F *   fhFracSumEnCellsPerSMNHigh20[fgkNCellEnMinBinsFr]; //!<! sum of cells with E > 1 GeV over emin, 1 cluster with n_cell_w>20, per SM

  TH2F *   fhSumEnCellsPerSMAcceptEvent     [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy, LED rejected, per SM
  TH2F *   fhNCellsPerSMAcceptEvent         [fgkNCellEnMinBins];   //!<! For E cell > emin, count number of cells, LED rejected, per SM 
  TH2F *   fhAverSumEnCellsPerSMAcceptEvent [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy / total cells number, LED rejected, per SM
  TH2F *   fhFracNCellsPerSMAcceptEvent     [fgkNCellEnMinBinsFr]; //!<! total number of cells with E > emin over 0.5 GeV, LED rejected, per SM
  TH2F *   fhFracSumEnCellsPerSMAcceptEvent [fgkNCellEnMinBinsFr]; //!<! sum of cells with E > 1 GeV over emin, LED rejected, per SM
  
  TH2F *   fhSumEnCellsPerSMAcceptEventStrip[fgkNCellEnMinBins]; //!<! For E cell > emin, sum of cells energy, per SM, LED rejected, based on n strip 
  TH2F *   fhNCellsPerSMAcceptEventStrip    [fgkNCellEnMinBins]; //!<! For E cell > emin, count number of cells, per SM, LED rejected, based on n strip
  TH2F *   fhSumEnCellsPerSMAcceptEventBoth [fgkNCellEnMinBins]; //!<! For E cell > emin, sum of cells energy, per SM, LED rejected, based on n strip and SM activity 
  TH2F *   fhNCellsPerSMAcceptEventBoth     [fgkNCellEnMinBins]; //!<! For E cell > emin, count number of cells, per SM, LED rejected, based on n strip and SM activity 
  
  TH3F *   fhNSumEnCellsPerSM                [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy, per SM 
  TH3F *   fhNSumEnCellsPerSMAcceptEvent     [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy, per SM, LED rejected 
  TH3F *   fhNSumEnCellsPerSMAcceptEventStrip[fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy, per SM, LED rejected, based on n strip 
  TH3F *   fhNSumEnCellsPerSMAcceptEventBoth[fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy, per SM, LED rejected, based on n strip and SM activity 
  
  // Per Strip
  //
  TH1F *   fhSumEnCellsPerStrip        [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy in a strip 
  TH1F *   fhNCellsPerStrip            [fgkNCellEnMinBins];   //!<! For E cell > emin, count number of cells in a strip
  TH1F *   fhSumEnCellsPerStripNHigh20 [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy in a strip, 1 cluster with n_cell_w>20, per SM 
  TH1F *   fhNCellsPerStripNHigh20     [fgkNCellEnMinBins];   //!<! For E cell > emin, count number of cells in a strip, 1 cluster with n_cell_w>20, per SM 
 
  TH2F *   fhSumEnCellsPerStripPerSM   [fgkNCellEnMinBins];   //!<! For E cell > emin, sum of cells energy in a strip, per SM 
  TH2F *   fhNCellsPerStripPerSM       [fgkNCellEnMinBins];   //!<! For E cell > emin, count number of cells in a strip, per SM
  TH2F *   fhSumEnCellsPerStripPerSMNHigh20[fgkNCellEnMinBins]; //!<! For E cell > emin, sum of cells energy in a strip, per SM, 1 cluster with n_cell_w>20, per SM  
  TH2F *   fhNCellsPerStripPerSMNHigh20[fgkNCellEnMinBins];   //!<! For E cell > emin, count number of cells in a strip, per SM, 1 cluster with n_cell_w>20, per SM  
  
  TH1F *   fhSumEnCellsPerStripAcceptEventStrip[fgkNCellEnMinBins]; //!<! For E cell > emin, sum of cells energy in a strip, accept event based on n strip 
  TH1F *   fhNCellsPerStripAcceptEventStrip    [fgkNCellEnMinBins]; //!<! For E cell > emin, count number of cells in a strip, accept event based on n strip
  TH1F *   fhSumEnCellsPerStripAcceptEventBoth [fgkNCellEnMinBins]; //!<! For E cell > emin, sum of cells energy in a strip, accept event based on n strip and SM activity
  TH1F *   fhNCellsPerStripAcceptEventBoth     [fgkNCellEnMinBins]; //!<! For E cell > emin, count number of cells in a strip, accept event based on n strip and SM activity
  TH1F *   fhSumEnCellsPerStripEventAccept     [fgkNCellEnMinBins]; //!<! For E cell > emin, sum of cells energy in a strip, accept event based on SM activity
  TH1F *   fhNCellsPerStripEventAccept         [fgkNCellEnMinBins]; //!<! For E cell > emin, count number of cells in a strip, accept event based on SM activity
  
  TH2F *   fhSumEnCellsPerStripPerSMEventAccept     [fgkNCellEnMinBins]; //!<! For E cell > emin, sum of cells energy in a strip, per SM, accept event based on n strip 
  TH2F *   fhNCellsPerStripPerSMEventAccept         [fgkNCellEnMinBins]; //!<! For E cell > emin, count number of cells in a strip, per SM, accept event based on n strip
  TH2F *   fhSumEnCellsPerStripPerSMAcceptEventStrip[fgkNCellEnMinBins]; //!<! For E cell > emin, sum of cells energy in a strip, per SM, accept event based on n strip 
  TH2F *   fhNCellsPerStripPerSMAcceptEventStrip    [fgkNCellEnMinBins]; //!<! For E cell > emin, count number of cells in a strip, per SM, accept event based on n strip
  TH2F *   fhSumEnCellsPerStripPerSMAcceptEventBoth [fgkNCellEnMinBins]; //!<! For E cell > emin, sum of cells energy in a strip, per SM, accept event based on n strip 
  TH2F *   fhNCellsPerStripPerSMAcceptEventBoth     [fgkNCellEnMinBins]; //!<! For E cell > emin, count number of cells in a strip, per SM, accept event based on n strip
  
  TH2F *   fhNSumEnCellsPerStrip                [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy in a strip
  TH2F *   fhNSumEnCellsPerStripAcceptEvent     [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy in a strip, LED rejected 
  TH2F *   fhNSumEnCellsPerStripAcceptEventStrip[fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy in a strip, LED rejected, based on n strip
  TH2F *   fhNSumEnCellsPerStripAcceptEventBoth[fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy in a strip, LED rejected, based on n strip

  TH3F *   fhNSumEnCellsPerStripPerSM                [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy in a strip
  TH3F *   fhNSumEnCellsPerStripPerSMAcceptEvent     [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy in a strip, per SM, LED rejected 
  TH3F *   fhNSumEnCellsPerStripPerSMAcceptEventStrip[fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy in a strip, per SM, LED rejected, based on n strip
  TH3F *   fhNSumEnCellsPerStripPerSMAcceptEventBoth [fgkNCellEnMinBins];  //!<! For E cell > emin, N vs sum of cells energy in a strip, per SM, LED rejected, based on n strip
  
  //
  
  TH1F *   fhNStripsPerEvent;                                 //!<! Number of active strips per event, low activity SM3
  TH2F *   fhNStripsPerEventPerSM;                            //!<! Number of active strips per event, per SM, low activity SM3
  TH1F *   fhNStripsPerEventAccept;                           //!<! Number of active strips per accepted event, low activity SM3
  TH2F *   fhNStripsPerEventAcceptPerSM;                      //!<! Number of active strips per accepted event, per SM, low activity SM3
  
  TH2F *   fhSM3NCellsSumEnSuspiciousEvents;                  //!<! Control histogram to check SM3 activity when there is high activity event in other SM, after event cut
  TH2F *   fhSumEnCellsPerSMEventSuspicious;                  //!<! For E cell > emin, sum of cells energy, per SM 
  TH2F *   fhNCellsPerSMEventSuspicious;                      //!<! For E cell > emin, count number of cells, per SM 
  TH1F *   fhNStripsPerEventSuspicious;                       //!<! Number of active strips per event, accepted event but suspicious
  TH2F *   fhNStripsPerEventSuspiciousPerSM;                  //!<! Number of active strips per event, accepted event but suspicious
 
  // Control event rejection, 4 cases: no selection, SM LED, Strip LED, both LED
  TH2F *   fhEventCutClusterEnergyTime[4];                    //!<! Cluster energy vs time, for event selection criteria, non exotic
  TH1F *   fhEventCutClusterEnergy[4];                        //!<! Cluster energy, for event selection criteria, non exotic, within time cut range
  TH2F *   fhEventCutClusterEnergyECellMax[4];                //!<! Cluster energy vs max cell energy, for event selection criteria, non exotic, within time cut range
  TH2F *   fhEventCutClusterEtaPhiGrid[4];                    //!<! Highest energy cell location in column-row map, for event selection criteria, non exotic, within time cut range, above fEMinForExo
  
  TH1F *   fhEventCutTrackPt[4];                              //!<! Track pT , for event selection criteria
  TH2F *   fhEventCutTrackEtaPhi[4];                          //!<! Track eta vs phi, pt > 0.5, for event selection criteria
  
  TH1F *   fhEventCutVertexX[4];                              //!<! vertex x  for event selection criteria
  TH1F *   fhEventCutVertexY[4];                              //!<! vertex y  for event selection criteria
  TH1F *   fhEventCutVertexZ[4];                              //!<! vertex z  for event selection criteria

  TH1I *   fhEventCutBunchCrossing[4];                        //!<! Bunch crossing of events for different event selection criteria

  
  /// Copy constructor not implemented.
  AliAnaCaloExotics & operator = (const AliAnaCaloExotics & qa) ;
    
  /// Assignment operator not implemented.
  AliAnaCaloExotics(              const AliAnaCaloExotics & qa) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaCaloExotics,11) ;
  /// \endcond

} ;

#endif //ALIANACALORIMETERQA_H



