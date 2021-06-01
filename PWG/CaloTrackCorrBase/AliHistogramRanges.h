#ifndef ALIHISTOGRAMRANGES_H
#define ALIHISTOGRAMRANGES_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliHistogramRanges
/// \ingroup CaloTrackCorrelationsBase
/// \brief Class containing more common histogram axis types.
///
/// Class containing more common histogram type settings:
/// * Number of bins.
/// * Minimum and maximum histograms axis.
/// The types are the pT, time, energy, phi, eta, etc.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// ROOT
#include <TObject.h>
#include <TArrayD.h>

class AliHistogramRanges : public TObject {
	
public:   
  
  AliHistogramRanges() ;              // default ctor
  virtual ~AliHistogramRanges() { ; } // dtor
  
  void InitParameters() ;
  
  void Print(const Option_t * ) const ;
  
  // Pt, Energy 
  
  Int_t   GetHistoPtBins()              const { return fHistoPtBins          ; }
  Float_t GetHistoPtMin()               const { return fHistoPtMin           ; }
  Float_t GetHistoPtMax()               const { return fHistoPtMax           ; }
  TArrayD GetHistoPtArr()               const { return fHistoPtArr           ; }
  void    SetHistoPtArr(TArrayD &arr)         { arr.Copy(fHistoPtArr)        ; }

  void    SetHistoPtRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoPtBins = n ; fHistoPtMax = max ; fHistoPtMin = min ; }
  
  Int_t   GetHistoEnergyBins()          const { return fHistoPtBins          ; }
  Float_t GetHistoEnergyMin()           const { return fHistoPtMin           ; }
  Float_t GetHistoEnergyMax()           const { return fHistoPtMax           ; }
  TArrayD GetHistoEnergyArr()           const { return fHistoPtArr           ; }
  void    SetHistoEnergyArr(TArrayD &arr)     { arr.Copy(fHistoPtArr)        ; }

  void    SetHistoEnergyRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    SetHistoPtRangeAndNBins(min, max, n) ; }
  
  
  Int_t   GetHistoWidePtBins()          const { return fHistoWidePtBins      ; }
  Float_t GetHistoWidePtMin()           const { return fHistoWidePtMin       ; }
  Float_t GetHistoWidePtMax()           const { return fHistoWidePtMax       ; }  
  TArrayD GetHistoWidePtArr()           const { return fHistoWidePtArr       ; }  
  void    SetHistoWidePtArr(TArrayD &arr)     { arr.Copy(fHistoWidePtArr)    ; }

   void SetHistoWidePtRangeAndNBins      (Float_t min, Float_t max, Int_t n) {
     fHistoWidePtBins = n ; fHistoWidePtMax = max ; fHistoWidePtMin = min ; }
  
  Int_t   GetHistoWideEnBins()          const { return fHistoWidePtBins      ; }
  Float_t GetHistoWideEnMin()           const { return fHistoWidePtMin       ; }
  Float_t GetHistoWideEnMax()           const { return fHistoWidePtMax       ; }  
  TArrayD GetHistoWideEnArr()           const { return fHistoWidePtArr       ; }  
  void    SetHistoWideEnArr(TArrayD &arr)     { arr.Copy(fHistoWidePtArr)    ; }

   void SetHistoWideEnRangeAndNBins      (Float_t min, Float_t max, Int_t n) {
     fHistoWidePtBins = n ; fHistoWidePtMax = max ; fHistoWidePtMin = min ; }
  
  Int_t   GetHistoFinePtBins()          const { return fHistoFinePtBins      ; }
  Float_t GetHistoFinePtMin()           const { return fHistoFinePtMin       ; }
  Float_t GetHistoFinePtMax()           const { return fHistoFinePtMax       ; }  
  TArrayD GetHistoFinePtArr()           const { return fHistoFinePtArr       ; }  
  void    SetHistoFinePtArr(TArrayD &arr)     { arr.Copy(fHistoFinePtArr)    ; }

  void SetHistoFinePtRangeAndNBins      (Float_t min, Float_t max, Int_t n) {
    fHistoFinePtBins = n ; fHistoFinePtMax = max ; fHistoFinePtMin = min ; }
  
  Int_t   GetHistoCellEnBins()          const { return fHistoFinePtBins      ; }
  Float_t GetHistoCellEnMin()           const { return fHistoFinePtMin       ; }
  Float_t GetHistoCellEnMax()           const { return fHistoFinePtMax       ; }  
  TArrayD GetHistoCellEnArr()           const { return fHistoFinePtArr       ; }  
  void    SetHistoCellEnArr(TArrayD &arr)     { arr.Copy(fHistoFinePtArr)    ; }

  void SetHistoCellEnRangeAndNBins      (Float_t min, Float_t max, Int_t n) {
    fHistoFinePtBins = n ; fHistoFinePtMax = max ; fHistoFinePtMin = min ; }
  
  // Azimuthal angle
  
  Int_t   GetHistoPhiBins()             const { return fHistoPhiBins         ; }
  Float_t GetHistoPhiMin()              const { return fHistoPhiMin          ; }
  Float_t GetHistoPhiMax()              const { return fHistoPhiMax          ; }
  TArrayD GetHistoPhiArr()              const { return fHistoPhiArr          ; }
  void    SetHistoPhiArr(TArrayD &arr)        { arr.Copy(fHistoPhiArr)       ; }

  void    SetHistoPhiRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoPhiBins = n ; fHistoPhiMax = max ; fHistoPhiMin = min ; }
  
  // Pseudorapidity-rapidity
   
   Int_t   GetHistoEtaBins()             const { return fHistoEtaBins         ; }
   Float_t GetHistoEtaMin()              const { return fHistoEtaMin          ; }
   Float_t GetHistoEtaMax()              const { return fHistoEtaMax          ; }
   TArrayD GetHistoEtaArr()              const { return fHistoEtaArr          ; }
   void    SetHistoEtaArr(TArrayD &arr)        { arr.Copy(fHistoEtaArr)       ; }

   void    SetHistoEtaRangeAndNBins(Float_t min, Float_t max, Int_t n) {
     fHistoEtaBins = n ; fHistoEtaMax = max ; fHistoEtaMin = min ; }
  
  // Delta azymuthal angle
  
  Int_t   GetHistoDeltaPhiBins()        const { return fHistoDeltaPhiBins    ; }
  Float_t GetHistoDeltaPhiMin()         const { return fHistoDeltaPhiMin     ; }
  Float_t GetHistoDeltaPhiMax()         const { return fHistoDeltaPhiMax     ; }
  TArrayD GetHistoDeltaPhiArr()         const { return fHistoDeltaPhiArr     ; }
  void    SetHistoDeltaPhiArr(TArrayD &arr)   { arr.Copy(fHistoDeltaPhiArr)  ; }

  void    SetHistoDeltaPhiRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoDeltaPhiBins   = n ; fHistoDeltaPhiMax   = max ; fHistoDeltaPhiMin   = min ; }


  // Delta eta angle

  Int_t   GetHistoDeltaEtaBins()        const { return fHistoDeltaEtaBins     ; }
  Float_t GetHistoDeltaEtaMin()         const { return fHistoDeltaEtaMin      ; }
  Float_t GetHistoDeltaEtaMax()         const { return fHistoDeltaEtaMax      ; }
  TArrayD GetHistoDeltaEtaArr()         const { return fHistoDeltaEtaArr      ; }
  void    SetHistoDeltaEtaArr(TArrayD &arr)   { arr.Copy(fHistoDeltaEtaArr)   ; }

  void    SetHistoDeltaEtaRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoDeltaEtaBins  = n ; fHistoDeltaEtaMax   = max ; fHistoDeltaEtaMin   = min ; }
  
  // Mass
	
  Int_t   GetHistoMassBins()            const { return fHistoMassBins        ; }
  Float_t GetHistoMassMin()             const { return fHistoMassMin         ; }
  Float_t GetHistoMassMax()             const { return fHistoMassMax         ; }
  TArrayD GetHistoMassArr()             const { return fHistoMassArr         ; }
  void    SetHistoMassArr(TArrayD &arr)       { arr.Copy(fHistoMassArr)      ; }

  void    SetHistoMassRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoMassBins = n ; fHistoMassMax = max ; fHistoMassMin = min ; }
  
  // Asymetry
	
  Int_t   GetHistoAsymmetryBins()        const { return fHistoAsymBins       ; }
  Float_t GetHistoAsymmetryMin()         const { return fHistoAsymMin        ; }
  Float_t GetHistoAsymmetryMax()         const { return fHistoAsymMax        ; }	
  TArrayD GetHistoAsymmetryArr()         const { return fHistoAsymArr        ; }  
  void    SetHistoAsymmetryArr(TArrayD &arr)   { arr.Copy(fHistoAsymArr)     ; }

  void    SetHistoAsymmetryRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoAsymBins = n ; fHistoAsymMax = max ; fHistoAsymMin = min ; }
  
  // VZero
  
  Int_t   GetHistoV0SignalBins()         const { return fHistoV0SBins        ; }
  Int_t   GetHistoV0SignalMin()          const { return fHistoV0SMin         ; }
  Int_t   GetHistoV0SignalMax()          const { return fHistoV0SMax         ; }
  TArrayD GetHistoV0SignalArr()          const { return fHistoV0SArr         ; }
  void    SetHistoV0SignalArr(TArrayD &arr)    { arr.Copy(fHistoV0SArr)      ; }

  void    SetHistoV0SignalRangeAndNBins(Int_t min, Int_t max, Int_t n) {
    fHistoV0SBins = n ; fHistoV0SMax = max ; fHistoV0SMin = min ; }
  
  Int_t   GetHistoV0MultiplicityBins()   const { return fHistoV0MBins        ; }
  Int_t   GetHistoV0MultiplicityMin()    const { return fHistoV0MMin         ; }
  Int_t   GetHistoV0MultiplicityMax()    const { return fHistoV0MMax         ; }
  TArrayD GetHistoV0MultiplicityArr()    const { return fHistoV0MMax         ; }
  void    SetHistoV0MultiplicityArr(TArrayD &arr) { arr.Copy(fHistoV0MArr)   ; }

  void    SetHistoV0MultiplicityRangeAndNBins(Int_t min, Int_t max, Int_t n) {
    fHistoV0MBins = n ; fHistoV0MMax = max ; fHistoV0MMin = min ; }
  
  // Track multiplicity
  
  Int_t   GetHistoTrackMultiplicityBins() const { return fHistoTrMBins       ; }
  Int_t   GetHistoTrackMultiplicityMin()  const { return fHistoTrMMin        ; }
  Int_t   GetHistoTrackMultiplicityMax()  const { return fHistoTrMMax        ; }
  TArrayD GetHistoTrackMultiplicityArr()  const { return fHistoTrMArr        ; }
  void    SetHistoTrackMultiplicityArr(TArrayD &arr) { arr.Copy(fHistoTrMArr); }

  void    SetHistoTrackMultiplicityRangeAndNBins(Int_t min, Int_t max, Int_t n) {
    fHistoTrMBins = n ; fHistoTrMMax = max ; fHistoTrMMin = min ; }
  
  // dEdx
  
  Int_t   GetHistodEdxBins()             const { return fHistodEdxBins       ; }
  Float_t GetHistodEdxMin()              const { return fHistodEdxMin        ; }
  Float_t GetHistodEdxMax()              const { return fHistodEdxMax        ; }  
  TArrayD GetHistodEdxArr()              const { return fHistodEdxArr        ; }  
  void    SetHistodEdxArr(TArrayD &arr)        { arr.Copy(fHistodEdxArr)     ; }

  void    SetHistodEdxRangeAndNBins        (Float_t min, Float_t max, Int_t n) {
    fHistodEdxBins = n ; fHistodEdxMax = max ; fHistodEdxMin = min ; }
  
  // E over p
  
  Int_t   GetHistoEOverPBins()           const { return fHistoEOverPBins     ; }
  Float_t GetHistoEOverPMin()            const { return fHistoEOverPMin      ; }
  Float_t GetHistoEOverPMax()            const { return fHistoEOverPMax      ; }
  TArrayD GetHistoEOverPArr()            const { return fHistoEOverPArr      ; }
  void    SetHistoEOverPArr(TArrayD &arr)      { arr.Copy(fHistoEOverPArr)   ; }

  void    SetHistoEOverPRangeAndNBins      (Float_t min, Float_t max, Int_t n) {
    fHistoEOverPBins = n ; fHistoEOverPMax = max ; fHistoEOverPMin = min ; }

  Int_t   GetHistoNSigmaBins()           const { return fHistoNSigmaBins     ; }
  Float_t GetHistoNSigmaMin()            const { return fHistoNSigmaMin      ; }
  Float_t GetHistoNSigmaMax()            const { return fHistoNSigmaMax      ; }
  TArrayD GetHistoNSigmaArr()            const { return fHistoNSigmaArr      ; }
  void    SetHistoNSigmaArr(TArrayD &arr)      { arr.Copy(fHistoNSigmaArr)   ; }

  void    SetHistoNSigmaRangeAndNBins      (Float_t min, Float_t max, Int_t n) {
    fHistoNSigmaBins = n ; fHistoNSigmaMax = max ; fHistoNSigmaMin = min ; }
  
  // Number of cells per clusters
  
  Int_t   GetHistoNClusterCellBins()     const { return fHistoNClusCellBins  ; }
  Int_t   GetHistoNClusterCellMin()      const { return fHistoNClusCellMin   ; }
  Int_t   GetHistoNClusterCellMax()      const { return fHistoNClusCellMax   ; }  
  TArrayD GetHistoNClusterCellArr()      const { return fHistoNClusCellArr   ; }  
  void    SetHistoNClusterCellArr(TArrayD &arr){ arr.Copy(fHistoNClusCellArr); }

  void    SetHistoNClusterCellRangeAndNBins(Int_t   min, Int_t   max, Int_t n) {
    fHistoNClusCellBins = n ; fHistoNClusCellMax = max ; fHistoNClusCellMin = min ; }
  
  // Number of clusters
  
  Int_t   GetHistoNClustersBins()        const { return fHistoNClustersBins  ; }
  Int_t   GetHistoNClustersMin()         const { return fHistoNClustersMin   ; }
  Int_t   GetHistoNClustersMax()         const { return fHistoNClustersMax   ; }  
  TArrayD GetHistoNClustersArr()         const { return fHistoNClustersArr   ; }  
  void    SetHistoNClustersArr(TArrayD &arr)   { arr.Copy(fHistoNClustersArr); }

  void    SetHistoNClustersRangeAndNBins   (Int_t   min, Int_t   max, Int_t n) {
    fHistoNClustersBins = n ; fHistoNClustersMax = max ; fHistoNClustersMin = min ; }

  // Number of cells 
  
  Int_t   GetHistoNCellsBins()           const { return fHistoNCellsBins     ; }
  Int_t   GetHistoNCellsMin()            const { return fHistoNCellsMin      ; }
  Int_t   GetHistoNCellsMax()            const { return fHistoNCellsMax      ; }  
  TArrayD GetHistoNCellsArr()            const { return fHistoNCellsArr      ; }  
  void    SetHistoNCellsArr(TArrayD &arr)      { arr.Copy(fHistoNCellsArr)   ; }

  void    SetHistoNCellsRangeAndNBins      (Int_t   min, Int_t   max, Int_t n) {
    fHistoNCellsBins = n ; fHistoNCellsMax = max ; fHistoNCellsMin = min ; }

  // dR
  
  Int_t   GetHistodRBins()               const { return fHistodRBins         ; }
  Float_t GetHistodRMin()                const { return fHistodRMin          ; }
  Float_t GetHistodRMax()                const { return fHistodRMax          ; }  
  TArrayD GetHistodRArr()                const { return fHistodRArr          ; }  
  void    SetHistodRArr(TArrayD &arr)          { arr.Copy(fHistodRArr)       ; }

  void    SetHistodRRangeAndNBins          (Float_t min, Float_t max, Int_t n) {
    fHistodRBins = n ; fHistodRMax = max ; fHistodRMin = min ; }

  // Ratio
  
  Int_t   GetHistoRatioBins()            const { return fHistoRatioBins      ; }
  Float_t GetHistoRatioMin()             const { return fHistoRatioMin       ; }
  Float_t GetHistoRatioMax()             const { return fHistoRatioMax       ; }  
  TArrayD GetHistoRatioArr()             const { return fHistoRatioArr       ; }  
  void    SetHistoRatioArr(TArrayD &arr)       { arr.Copy(fHistoRatioArr)    ; }

  void    SetHistoRatioRangeAndNBins       (Float_t min, Float_t max, Int_t n) {
    fHistoRatioBins = n ; fHistoRatioMax = max ; fHistoRatioMin = min ; }

  // Ratio1 (maximum should be 1)
   
  Int_t   GetHistoRatio1Bins()            const { return fHistoRatio1Bins      ; }
  Float_t GetHistoRatio1Min()             const { return fHistoRatio1Min       ; }
  Float_t GetHistoRatio1Max()             const { return fHistoRatio1Max       ; }  
  TArrayD GetHistoRatio1Arr()             const { return fHistoRatio1Arr       ; }  
  void    SetHistoRatio1Arr(TArrayD &arr)       { arr.Copy(fHistoRatio1Arr)    ; }
  
  void    SetHistoRatio1RangeAndNBins       (Float_t min, Float_t max, Int_t n) {
    fHistoRatio1Bins = n ; fHistoRatio1Max = max ; fHistoRatio1Min = min ; }
  
  // Energy difference
  
  Int_t   GetHistoEDiffBins()            const { return fHistoEDiffBins      ; }
  Float_t GetHistoEDiffMin()             const { return fHistoEDiffMin       ; }
  Float_t GetHistoEDiffMax()             const { return fHistoEDiffMax       ; }  
  TArrayD GetHistoEDiffArr()             const { return fHistoEDiffArr       ; }  
  void    SetHistoEDiffArr(TArrayD &arr)       { arr.Copy(fHistoEDiffArr)    ; }

  void    SetHistoEDiffRangeAndNBins       (Float_t min, Float_t max, Int_t n) {
    fHistoEDiffBins = n ; fHistoEDiffMax = max ; fHistoEDiffMin  = min ; }
  
  // Hump-Backed Plateau
  
  Int_t   GetHistoHBPBins()              const { return fHistoHBPBins        ; }
  Float_t GetHistoHBPMin()               const { return fHistoHBPMin         ; }
  Float_t GetHistoHBPMax()               const { return fHistoHBPMax         ; }
  TArrayD GetHistoHBPArr()               const { return fHistoHBPArr         ; }
  void    SetHistoHBPArr(TArrayD &arr)         { arr.Copy(fHistoHBPArr)      ; }

  void    SetHistoHBPRangeAndNBins       (Float_t min, Float_t max, Int_t n) {
    fHistoHBPBins = n ; fHistoHBPMax = max ; fHistoHBPMin = min ; }

  // Vertex
  
  Int_t   GetHistoVertexDistBins()       const { return fHistoVertexDistBins ; }
  Float_t GetHistoVertexDistMin()        const { return fHistoVertexDistMin  ; }
  Float_t GetHistoVertexDistMax()        const { return fHistoVertexDistMax  ; }  
  TArrayD GetHistoVertexDistArr()        const { return fHistoVertexDistArr  ; }  
  void    SetHistoVertexDistArr(TArrayD &arr)  { arr.Copy(fHistoVertexDistArr); }

  void    SetHistoVertexDistRangeAndNBins  (Float_t min, Float_t max, Int_t n) { 
    fHistoVertexDistBins = n ; fHistoVertexDistMax = max ; fHistoVertexDistMin = min ; }
  
  // R =sqrt(x^2+y^2+z^2) (cm)
  
  Int_t   GetHistoRBins()                const { return fHistoRBins          ; }
  Float_t GetHistoRMin()                 const { return fHistoRMin           ; }
  Float_t GetHistoRMax()                 const { return fHistoRMax           ; }    
  TArrayD GetHistoRArr()                 const { return fHistoRArr           ; }    
  void    SetHistoRArr(TArrayD &arr)           { arr.Copy(fHistoRArr)        ; }

  void    SetHistoRRangeAndNBins           (Float_t min, Float_t max, Int_t n) {
    fHistoRBins = n ; fHistoRMax = max ; fHistoRMin = min ; }
  
  // X position
  
  Int_t   GetHistoXBins()                const { return fHistoXBins          ; }
  Float_t GetHistoXMin()                 const { return fHistoXMin           ; }
  Float_t GetHistoXMax()                 const { return fHistoXMax           ; }
  TArrayD GetHistoXArr()                 const { return fHistoXArr           ; }
  void    SetHistoXArr(TArrayD &arr)           { arr.Copy(fHistoXArr)        ; }

  void    SetHistoXRangeAndNBins           (Float_t min, Float_t max, Int_t n) {
    fHistoXBins = n ; fHistoXMax = max ; fHistoXMin = min ; }

  // Y position
  
  Int_t   GetHistoYBins()                const { return fHistoYBins          ; }
  Float_t GetHistoYMin()                 const { return fHistoYMin           ; }
  Float_t GetHistoYMax()                 const { return fHistoYMax           ; }  
  TArrayD GetHistoYArr()                 const { return fHistoYArr           ; }  
  void    SetHistoYArr(TArrayD &arr)           { arr.Copy(fHistoYArr)        ; }

  void    SetHistoYRangeAndNBins           (Float_t min, Float_t max, Int_t n) {
    fHistoYBins = n ; fHistoYMax = max ; fHistoYMin = min ; }
  
  // Z position
  
  Int_t   GetHistoZBins()                const { return fHistoZBins          ; }
  Float_t GetHistoZMin()                 const { return fHistoZMin           ; }
  Float_t GetHistoZMax()                 const { return fHistoZMax           ; }  
  TArrayD GetHistoZArr()                 const { return fHistoZArr           ; }  
  void    SetHistoZArr(TArrayD &arr)           { arr.Copy(fHistoZArr)        ; }

  void    SetHistoZRangeAndNBins           (Float_t min, Float_t max, Int_t n) {
    fHistoZBins = n ; fHistoZMax = max ; fHistoZMin = min ; }

  // Shower shape parameters
  
  Int_t   GetHistoShowerShapeBins()      const { return fHistoSSBins         ; }
  Float_t GetHistoShowerShapeMin()       const { return fHistoSSMin          ; }
  Float_t GetHistoShowerShapeMax()       const { return fHistoSSMax          ; }  
  TArrayD GetHistoShowerShapeArr()       const { return fHistoSSArr          ; }  
  void    SetHistoShowerShapeArr(TArrayD &arr) { arr.Copy(fHistoSSArr)       ; }

  void    SetHistoShowerShapeRangeAndNBins (Float_t min, Float_t max, Int_t n) {
    fHistoSSBins = n ; fHistoSSMax = max ; fHistoSSMin = min ; }

  // Time
  
  Int_t   GetHistoTimeBins()             const { return fHistoTimeBins       ; }
  Float_t GetHistoTimeMin()              const { return fHistoTimeMin        ; }
  Float_t GetHistoTimeMax()              const { return fHistoTimeMax        ; }  
  TArrayD GetHistoTimeArr()              const { return fHistoTimeArr        ; }  
  void    SetHistoTimeArr(TArrayD &arr)        { arr.Copy(fHistoTimeArr)     ; }

  void    SetHistoTimeRangeAndNBins        (Float_t min, Float_t max, Int_t n) {
    fHistoTimeBins = n ; fHistoTimeMax = max ; fHistoTimeMin = min ; }	  
  
  // Cluster time difference
  
  Int_t   GetHistoDiffTimeBins()         const { return fHistoDiffTimeBins   ; }
	Float_t GetHistoDiffTimeMin()          const { return fHistoDiffTimeMin    ; }
  Float_t GetHistoDiffTimeMax()          const { return fHistoDiffTimeMax    ; }  
  TArrayD GetHistoDiffTimeArr()          const { return fHistoDiffTimeArr    ; }  
  void    SetHistoDiffTimeArr(TArrayD &arr)    { arr.Copy(fHistoDiffTimeArr) ; }

  void    SetHistoDiffTimeRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoDiffTimeBins = n ; fHistoDiffTimeMax = max ; fHistoDiffTimeMin = min   ; }
  
  //    Track matching histogrammes setters and getters
  
  Int_t   GetHistoTrackResidualEtaBins() const { return fHistoTrackResidualEtaBins ; }
	Float_t GetHistoTrackResidualEtaMin()  const { return fHistoTrackResidualEtaMin  ; }
  Float_t GetHistoTrackResidualEtaMax()  const { return fHistoTrackResidualEtaMax  ; }    
  TArrayD GetHistoTrackResidualEtaArr()  const { return fHistoTrackResidualEtaArr  ; }    
  void    SetHistoTrackResidualEtaArr(TArrayD &arr) { arr.Copy(fHistoTrackResidualEtaArr) ; }

  void    SetHistoTrackResidualEtaRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoTrackResidualEtaBins = n ; fHistoTrackResidualEtaMax = max ; fHistoTrackResidualEtaMin = min ; }
  
  Int_t   GetHistoTrackResidualPhiBins() const { return fHistoTrackResidualPhiBins ; }
	Float_t GetHistoTrackResidualPhiMin()  const { return fHistoTrackResidualPhiMin  ; }
  Float_t GetHistoTrackResidualPhiMax()  const { return fHistoTrackResidualPhiMax  ; }  
  TArrayD GetHistoTrackResidualPhiArr()  const { return fHistoTrackResidualPhiArr  ; }  
  void    SetHistoTrackResidualPhiArr(TArrayD &arr) { arr.Copy(fHistoTrackResidualPhiArr) ; }

  void    SetHistoTrackResidualPhiRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoTrackResidualPhiBins = n ; fHistoTrackResidualPhiMax = max ; fHistoTrackResidualPhiMin = min ; }
  
  // Isolation task, sum pt
  
  // Sum in cone
  void    SetHistoPtSumRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNPtSumBins = n ;  fHistoPtSumMax = max ; fHistoPtSumMin = min ; }
  
  Int_t   GetHistoNPtSumBins()           const { return fHistoNPtSumBins     ; }
  Float_t GetHistoPtSumMin()             const { return fHistoPtSumMin       ; }
  Float_t GetHistoPtSumMax()             const { return fHistoPtSumMax       ; }
  TArrayD GetHistoPtSumArr()             const { return fHistoPtSumArr       ; }
  void    SetHistoPtSumArr(TArrayD &arr)       { arr.Copy(fHistoPtSumArr)    ; }

  // Sum in cone after subtraction
  void    SetHistoPtSumSubRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNPtSumSubBins = n ;  fHistoPtSumSubMax = max ; fHistoPtSumSubMin = min ; }
  
  Int_t   GetHistoNPtSumSubBins()           const { return fHistoNPtSumSubBins     ; }
  Float_t GetHistoPtSumSubMin()             const { return fHistoPtSumSubMin       ; }
  Float_t GetHistoPtSumSubMax()             const { return fHistoPtSumSubMax       ; }
  TArrayD GetHistoPtSumSubArr()             const { return fHistoPtSumSubArr       ; }
  void    SetHistoPtSumSubArr(TArrayD &arr)       { arr.Copy(fHistoPtSumSubArr)    ; }
  
  // Pt distribution in cone
  void    SetHistoPtInConeRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNPtInConeBins = n ; fHistoPtInConeMax = max ; fHistoPtInConeMin = min  ; }
  
  Int_t   GetHistoNPtInConeBins()        const { return fHistoNPtInConeBins  ; }
  Float_t GetHistoPtInConeMin()          const { return fHistoPtInConeMin    ; }
  Float_t GetHistoPtInConeMax()          const { return fHistoPtInConeMax    ; }
  TArrayD GetHistoPtInConeArr()          const { return fHistoPtInConeArr    ; }
  void    SetHistoPtInConeArr(TArrayD &arr)    { arr.Copy(fHistoPtInConeArr) ; }

  // Opening angle
  
  void    SetHistoOpeningAngleRangeAndNBins(Float_t min, Float_t max, Int_t n) {
   fHistoOpAngleBins = n ; fHistoOpAngleMax = max ; fHistoOpAngleMin = min  ; }
  
  Int_t   GetHistoNOpeningAngleBins()    const { return fHistoOpAngleBins    ; }
  Float_t GetHistoOpeningAngleMin()      const { return fHistoOpAngleMin     ; }
  Float_t GetHistoOpeningAngleMax()      const { return fHistoOpAngleMax     ; }
  TArrayD GetHistoOpeningAngleArr()      const { return fHistoOpAngleArr     ; }
  void    SetHistoOpeningAngleArr(TArrayD &arr){ arr.Copy(fHistoOpAngleArr)  ; }

  // Centrality
  
  Int_t   GetHistoCentralityBins()       const { return fHistoCenBins        ; }
  Float_t GetHistoCentralityMin()        const { return fHistoCenMin         ; }
  Float_t GetHistoCentralityMax()        const { return fHistoCenMax         ; }
  TArrayD GetHistoCentralityArr()        const { return fHistoCenArr         ; }
  void    SetHistoCentralityArr(TArrayD &arr)  { arr.Copy(fHistoCenArr)      ; }

  void    SetHistoCentralityRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoCenBins = n ; fHistoCenMax = max ; fHistoCenMin = min ; }
  
  // Number of local maxima
  
  Int_t   GetHistoNLMBins()       const { return fHistoNLMBins        ; }
  Int_t   GetHistoNLMMin()        const { return fHistoNLMMin         ; }
  Int_t   GetHistoNLMMax()        const { return fHistoNLMMax         ; }
  TArrayD GetHistoNLMArr()        const { return fHistoNLMArr         ; }
  void    SetHistoNLMArr(TArrayD &arr)  { arr.Copy(fHistoNLMArr)      ; }
  
  void    SetHistoNLMRangeAndNBins(Int_t min, Int_t max, Int_t n) {
    fHistoNLMBins = n ; fHistoNLMMax = max ; fHistoNLMMin = min ; }
  
  // Number of overlaps
  
  Int_t   GetHistoNoverlapBins()       const { return fHistoNoverlapBins    ; }
  Int_t   GetHistoNoverlapMin()        const { return fHistoNoverlapMin     ; }
  Int_t   GetHistoNoverlapMax()        const { return fHistoNoverlapMax     ; }
  TArrayD GetHistoNoverlapArr()        const { return fHistoNoverlapArr     ; }
  void    SetHistoNoverlapArr(TArrayD &arr)  { arr.Copy(fHistoNoverlapArr)  ; }
  
  void    SetHistoNoverlapRangeAndNBins(Int_t min, Int_t max, Int_t n) {
    fHistoNoverlapBins = n ; fHistoNoverlapMax = max ; fHistoNoverlapMin = min ; }

  // Exoticity
  
  Int_t   GetHistoExoticityBins()       const { return fHistoExoticityBins    ; }
  Float_t GetHistoExoticityMin()        const { return fHistoExoticityMin     ; }
  Float_t GetHistoExoticityMax()        const { return fHistoExoticityMax     ; }
  TArrayD GetHistoExoticityArr()        const { return fHistoExoticityArr     ; }
  void    SetHistoExoticityArr(TArrayD &arr)  { arr.Copy(fHistoExoticityArr)  ; }
  
  void    SetHistoExoticityRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoExoticityBins = n ; fHistoExoticityMax = max ; fHistoExoticityMin = min ; }
  
  // Spherocity

  Int_t   GetHistoSpherocityBins()       const { return fHistoSpherocityBins    ; }
  Float_t GetHistoSpherocityMin()        const { return fHistoSpherocityMin     ; }
  Float_t GetHistoSpherocityMax()        const { return fHistoSpherocityMax     ; }
  TArrayD GetHistoSpherocityArr()        const { return fHistoSpherocityArr     ; }
  void    SetHistoSpherocityArr(TArrayD &arr)  { arr.Copy(fHistoSpherocityArr)  ; }

  void    SetHistoSpherocityRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoSpherocityBins = n ; fHistoSpherocityMax = max ; fHistoSpherocityMin = min ; }

private:    
  
  Int_t    fHistoPtBins    ;                  ///< Number of bins in pt axis.
  Float_t  fHistoPtMax     ;                  ///< Maximum value of pt histogram range.
  Float_t  fHistoPtMin     ;                  ///< Minimum value of pt histogram range.
  TArrayD  fHistoPtArr     ;                  ///< Pt histogram lower limit bins.
  Int_t    fHistoPhiBins   ;                  ///< Number of bins in phi axis.
  Float_t  fHistoPhiMax    ;                  ///< Maximum value of phi histogram range.
  Float_t  fHistoPhiMin    ;                  ///< Minimum value of phi histogram range.
  TArrayD  fHistoPhiArr    ;                  ///< Phi histogram lower limit bins.
  Int_t    fHistoEtaBins   ;                  ///< Number of bins in eta axis.
  Float_t  fHistoEtaMax    ;                  ///< Maximum value of eta histogram range.
  Float_t  fHistoEtaMin    ;                  ///< Minimum value of eta histogram range.
  TArrayD  fHistoEtaArr    ;                  ///< eta histogram lower limit bins.
  Int_t    fHistoDeltaPhiBins  ;              ///< Number of bins in delta phi axis.
  Float_t  fHistoDeltaPhiMax   ;              ///< Maximum value of delta phi histogram range.
  Float_t  fHistoDeltaPhiMin   ;              ///< Minimum value of delta phi histogram range.
  TArrayD  fHistoDeltaPhiArr   ;              ///< Delta phi histogram lower limit bins.
  Int_t    fHistoDeltaEtaBins  ;              ///< Number of bins in delta eta axis.
  Float_t  fHistoDeltaEtaMax   ;              ///< Maximum value of delta eta histogram range.
  Float_t  fHistoDeltaEtaMin   ;              ///< Minimum value of delta eta histogram range.
  TArrayD  fHistoDeltaEtaArr   ;              ///< Delta eta histogram lower limit bins.
  Int_t    fHistoMassBins  ;                  ///< Number of bins in mass axis.
  Float_t  fHistoMassMax   ;                  ///< Maximum value of mass histogram range.
  Float_t  fHistoMassMin   ;                  ///< Minimum value of mass histogram range.
  TArrayD  fHistoMassArr   ;                  ///< Mass histogram lower limit bins.
  Int_t    fHistoAsymBins  ;                  ///< Number of bins in asymmetry axis.
  Float_t  fHistoAsymMax   ;                  ///< Maximum value of asymmetry histogram range.
  Float_t  fHistoAsymMin   ;                  ///< Minimum value of asymmetry histogram range.
  TArrayD  fHistoAsymArr   ;                  ///<  Asymmetry histogram lower limit bins.
  Int_t    fHistoV0SBins   ;                  ///< Number of bins in V0 signal axis.
  Int_t    fHistoV0SMax    ;                  ///< Maximum value of V0 signal histogram range.
  Int_t    fHistoV0SMin    ;                  ///< Minimum value of V0 signal histogram range.
  TArrayD  fHistoV0SArr    ;                  ///< V0 signal histogram lower limit bins.
  Int_t    fHistoV0MBins   ;                  ///< Number of bins in V0 multiplicity axis.
  Int_t    fHistoV0MMax    ;                  ///< Maximum value of V0 multiplicity histogram range.
  Int_t    fHistoV0MMin    ;                  ///< Minimum value of V0 multiplicity histogram range.
  TArrayD  fHistoV0MArr    ;                  ///<  V0 multiplicity histogram lower limit bins.
  Int_t    fHistoTrMBins   ;                  ///< Number of bins in V0 multiplicity axis.
  Int_t    fHistoTrMMax    ;                  ///< Maximum value of track multiplicity histogram range.
  Int_t    fHistoTrMMin    ;                  ///< Minimum value of track multiplicity histogram range.
  TArrayD  fHistoTrMArr    ;                  ///< Track multiplicity histogram lower limit bins.
  Int_t    fHistoFinePtBins;                  ///< Fne binning for fhAmpId histogram
  Float_t  fHistoFinePtMax ;                  ///< Maximum pt value for fhAmpId histogram
  Float_t  fHistoFinePtMin ;                  ///< Minimum pt value for fhAmpId histogram
  TArrayD  fHistoFinePtArr ;                  ///< For fhAmpid histogram lower limit bins.
  Int_t    fHistoWidePtBins;                  ///< Wide binning for acceptance dependent histogram
  Float_t  fHistoWidePtMax ;                  ///< Maximum pt value for acceptance dependent histogram
  Float_t  fHistoWidePtMin ;                  ///< Minimum pt value for acceptance dependent histogram
  TArrayD  fHistoWidePtArr ;                  ///< For acceptance dependent histogram lower limit bins.
  Int_t    fHistoEOverPBins;                  ///< E/p histogram number of bins.
  Float_t  fHistoEOverPMax ;                  ///< E/p maximum value.
  Float_t  fHistoEOverPMin ;                  ///< E/p minimum value.
  TArrayD  fHistoEOverPArr ;                  ///<  E/p histogram lower limit bins.
  Int_t    fHistoNSigmaBins;                  ///< TPC nSigma histogram number of bins.
  Float_t  fHistoNSigmaMax ;                  ///< TPC nSigma maximum value.
  Float_t  fHistoNSigmaMin ;                  ///< TPC nSigma minimum value. 
  TArrayD  fHistoNSigmaArr ;                  ///< nSigma  histogram lower limit bins.
  Int_t    fHistodEdxBins  ;                  ///< dEdx histogram number of bins.
  Float_t  fHistodEdxMax   ;                  ///< dEdx maximum value.
  Float_t  fHistodEdxMin   ;                  ///< dEdx minimum value.
  TArrayD  fHistodEdxArr   ;                  ///<  dEdx histogram lower limit bins.
  Int_t    fHistodRBins    ;                  ///< dR histogram number of bins.
  Float_t  fHistodRMax     ;                  ///< dR maximum value.
  Float_t  fHistodRMin     ;                  ///< dR minimum value.
  TArrayD  fHistodRArr     ;                  ///< dR  histogram lower limit bins.
  Int_t    fHistoTimeBins  ;                  ///< Cell time histogram number of bins.
  Float_t  fHistoTimeMax   ;                  ///< Cell time maximum value.
  Float_t  fHistoTimeMin   ;                  ///< Cell time minimum value.
  TArrayD  fHistoTimeArr   ;                  ///< Time  histogram lower limit bins.
  Int_t    fHistoNClusCellBins;               ///< Number of cells per cluster histogram number of bins.
  Int_t    fHistoNClusCellMax ;               ///< Number of cells per cluster maximum value.
  Int_t    fHistoNClusCellMin ;               ///< Number of cells per cluster minimum value.
  TArrayD  fHistoNClusCellArr ;               ///< Number of cells per cluster histogram lower limit bins.
  Int_t    fHistoNCellsBins;                  ///< Number of cells histogram number of bins.
  Int_t    fHistoNCellsMax ;                  ///< Number of cells maximum value.
  Int_t    fHistoNCellsMin ;                  ///< Number of cells minimum value.
  TArrayD  fHistoNCellsArr ;                  ///< Number of cells histogram lower limit bins.
  Int_t    fHistoNClustersBins;               ///< Number of clusters histogram number of bins.
  Int_t    fHistoNClustersMax ;               ///< Number of clusters maximum value.
  Int_t    fHistoNClustersMin ;               ///< Number of clusters minimum value.  
  TArrayD  fHistoNClustersArr ;               ///< Number of clusters histogram lower limit bins.
  Int_t    fHistoRatioBins ;                  ///< Ratio histogram number of bins.
  Float_t  fHistoRatioMax  ;                  ///< Ratio maximum value.
  Float_t  fHistoRatioMin  ;                  ///< Ratio minimum value.
  TArrayD  fHistoRatioArr  ;                  ///< Ratio histogram lower limit bins. Max is 1.
  Int_t    fHistoRatio1Bins ;                 ///< Ratio histogram number of bins. Max is 1.
  Float_t  fHistoRatio1Max  ;                 ///< Ratio maximum value. Max should be 1.
  Float_t  fHistoRatio1Min  ;                 ///< Ratio minimum value. Max is 1.
  TArrayD  fHistoRatio1Arr  ;                 ///< Ratio  histogram lower limit bins. Max is 1.
  Int_t    fHistoEDiffBins ;                  ///< Energy difference histogram number of bins.
  Float_t  fHistoEDiffMax  ;                  ///< Energy difference maximum value.
  Float_t  fHistoEDiffMin  ;                  ///< Energy difference minimum value.
  TArrayD  fHistoEDiffArr  ;                  ///< Energy differencee histogram lower limit bins.
  Int_t    fHistoHBPBins   ;                  ///< Hump-backed plateau histogram number of bins.
  Float_t  fHistoHBPMax    ;                  ///< Hump-backed plateau maximum value.
  Float_t  fHistoHBPMin    ;                  ///< Hump-backed plateau minimum value.
  TArrayD  fHistoHBPArr    ;                  ///< Hump-backerd plateau  histogram lower limit bins.
  Int_t    fHistoVertexDistBins;              ///< Vertex distance histogram number of bins.
  Float_t  fHistoVertexDistMax ;              ///< Vertex distance maximum value.
  Float_t  fHistoVertexDistMin ;              ///< Vertex distance minimum value.	
  TArrayD  fHistoVertexDistArr ;              ///< Vertex distance histogram lower limit bins.
  Int_t    fHistoRBins    ;                   ///< r =sqrt(x^2+y^2+z^2) (cm) position histogram number of bins.
  Float_t  fHistoRMax     ;                   ///< r =sqrt(x^2+y^2+z^2) (cm)  maximum value.
  Float_t  fHistoRMin     ;                   ///< r =sqrt(x^2+y^2+z^2) (cm)  minimum value.	
  TArrayD  fHistoRArr     ;                   ///< r =sqrt(x^2+y^2+z^2) (cm)  histogram lower limit bins.
  Int_t    fHistoXBins    ;                   ///< x (cm) position histogram number of bins.
  Float_t  fHistoXMax     ;                   ///< x (cm) position maximum value.
  Float_t  fHistoXMin     ;                   ///< x (cm) position minimum value.
  TArrayD  fHistoXArr     ;                   ///< x (cm) histogram lower limit bins.
  Int_t    fHistoYBins    ;                   ///< y (cm) position histogram number of bins.
  Float_t  fHistoYMax     ;                   ///< y (cm) position maximum value.
  Float_t  fHistoYMin     ;                   ///< y (cm) position minimum value.
  TArrayD  fHistoYArr     ;                   ///< y (cm) histogram lower limit bins.
  Int_t    fHistoZBins    ;                   ///< z (cm) position histogram number of bins.
  Float_t  fHistoZMax     ;                   ///< z (cm) position maximum value.
  Float_t  fHistoZMin     ;                   ///< z (cm) position minimum value.
  TArrayD  fHistoZArr     ;                   ///< z (cm) histogram lower limit bins.
  Int_t    fHistoSSBins   ;                   ///< Shower Shape parameter histogram number of bins.
  Float_t  fHistoSSMax    ;                   ///< Shower Shape parameter position maximum value.
  Float_t  fHistoSSMin    ;                   ///< Shower Shape parameter position minimum value.
  TArrayD  fHistoSSArr    ;                   ///< Shower Shape histogram lower limit bins.
  Int_t    fHistoDiffTimeBins;                ///< Difference cluster pair time parameter histogram number of bins.
  Float_t  fHistoDiffTimeMax ;                ///< Difference cluster pair time parameter maximum value.
  Float_t  fHistoDiffTimeMin ;                ///< Difference cluster pair time parameter minimum value.  
  TArrayD  fHistoDiffTimeArr ;                ///< Difference cluster pair time histogram lower limit bins.
  Int_t    fHistoTrackResidualEtaBins ;       ///< Number of bins in dEta (cluster-track) axis.
  Float_t  fHistoTrackResidualEtaMax  ;       ///< Maximum value of dEta (cluster-track) histogram range.
  Float_t  fHistoTrackResidualEtaMin  ;       ///< Minimum value of dEta (cluster-track) histogram range.		
  TArrayD  fHistoTrackResidualEtaArr  ;       ///< dEta (cluster-track) histogram lower limit bins.
  Int_t    fHistoTrackResidualPhiBins ;       ///< Number of bins in dPhi axis.
  Float_t  fHistoTrackResidualPhiMax  ;       ///< Maximum value of dPhi (cluster-track) histogram range.
  Float_t  fHistoTrackResidualPhiMin  ;       ///< Minimum value of dPhi (cluster-track) histogram range.
  TArrayD  fHistoTrackResidualPhiArr  ;       ///< dPhi (cluster-track) histogram lower limit bins.
  Int_t    fHistoNPtSumBins;                  ///< Number of bins in Isolation PtSum histograms.
  Float_t  fHistoPtSumMax  ;                  ///< Isolation PtSum maximum in histogram.
  Float_t  fHistoPtSumMin  ;	                ///< Isolation PtSum minimum in histogram.
  TArrayD  fHistoPtSumArr  ;                  ///< Isolation PtSum histogram lower limit bins.
  Int_t    fHistoNPtSumSubBins;               ///< Number of bins in Isolation PtSum histograms.
  Float_t  fHistoPtSumSubMax  ;               ///< Isolation UE subtracted PtSum maximum in histogram.
  Float_t  fHistoPtSumSubMin  ;               ///< Isolation UE subtracted PtSum minimum in histogram.
  TArrayD  fHistoPtSumSubArr  ;               ///< Isolation UE subtracted PtSum histogram lower limit bins.
  Int_t    fHistoNPtInConeBins;               ///< Number of bins in Isolation U subtracted PtInCone histogram.
  Float_t  fHistoPtInConeMax  ;               ///< Isolation PtInCone maximum in histogram.
  Float_t  fHistoPtInConeMin  ;               ///< Isolation PtInCone maximum in histogram.
  TArrayD  fHistoPtInConeArr  ;               ///< Isolation PtInCone histogram lower limit bins.
  Int_t    fHistoOpAngleBins;                 ///< Number of bins in pair opening angle axis.
  Float_t  fHistoOpAngleMax ;                 ///< Maximum value of pair opening angle histogram range.
  Float_t  fHistoOpAngleMin ;                 ///< Minimum value of pair opening angle histogram range.
  TArrayD  fHistoOpAngleArr ;                 ///< Opening angle histogram lower limit bins.
  Int_t    fHistoCenBins;                     ///< Number of bins in centrality axis.
  Float_t  fHistoCenMax ;                     ///< Maximum value of  centrality histogram range.
  Float_t  fHistoCenMin ;                     ///< Minimum value of centrality histogram range.
  TArrayD  fHistoCenArr ;                     ///< Centrality histogram lower limit bins.
  Int_t    fHistoNLMBins;                     ///< Number of bins in number of local maxima axis.
  Int_t    fHistoNLMMax ;                     ///< Maximum value of number of local maxima histogram range.
  Int_t    fHistoNLMMin ;                     ///< Minimum value of number of local maxima histogram range.
  TArrayD  fHistoNLMArr ;                     ///< Number of local maxima histogram lower limit bins.
  Int_t    fHistoNoverlapBins;                ///< Number of bins in number of overlaps axis.
  Int_t    fHistoNoverlapMax ;                ///< Maximum value of number of overlaps histogram range.
  Int_t    fHistoNoverlapMin ;                ///< Minimum value of number of overlaps histogram range.
  TArrayD  fHistoNoverlapArr ;                ///< Number of overlaps histogram lower limit bins.
  Int_t    fHistoExoticityBins;               ///< Number of bins in number of exoticity.
  Float_t  fHistoExoticityMax ;               ///< Maximum value of exoticity histogram range.
  Float_t  fHistoExoticityMin ;               ///< Minimum value of exoticity histogram range.
  TArrayD  fHistoExoticityArr ;               ///< Exoticity histogram lower limit bins.
  Int_t    fHistoSpherocityBins;              ///< Number of bins in number of spherocity.
  Float_t  fHistoSpherocityMax ;              ///< Maximum value of spherocity histogram range.
  Float_t  fHistoSpherocityMin ;              ///< Minimum value of spherocity histogram range.
  TArrayD  fHistoSpherocityArr ;              ///< Spherocity histogram lower limit bins.

  /// Copy constructor not implemented.
  AliHistogramRanges(              const AliHistogramRanges & h) ; 
  
  /// Assignment operator not implemented.
  AliHistogramRanges & operator = (const AliHistogramRanges & h) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliHistogramRanges,11) ;
  /// \endcond

} ;


#endif //ALIHISTOGRAMRANGES_H



