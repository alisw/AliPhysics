#ifndef ALIHISTOGRAMRANGES_H
#define ALIHISTOGRAMRANGES_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Class containing histogram settings:
//    - Number of bins
//    - Min and max of range
//
//-- Author: Gustavo Conesa (LPSC-Grenoble)
//

//ROOT
#include <TObject.h>

class AliHistogramRanges : public TObject {
	
public:   
  
  AliHistogramRanges() ;              // default ctor
  virtual ~AliHistogramRanges() { ; } // dtor
  
  void InitParameters() ;
  
  void Print(const Option_t * ) const ;
  
  //Pt, Energy 
  
  Int_t   GetHistoPtBins()              const { return fHistoPtBins          ; }
  Float_t GetHistoPtMin()               const { return fHistoPtMin           ; }
  Float_t GetHistoPtMax()               const { return fHistoPtMax           ; }
  
  void    SetHistoPtRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoPtBins         = n ; fHistoPtMax         = max ; fHistoPtMin         = min ; }
  
  Int_t   GetHistoEnergyBins()          const { return fHistoPtBins          ; }
  Float_t GetHistoEnergyMin()           const { return fHistoPtMin           ; }
  Float_t GetHistoEnergyMax()           const { return fHistoPtMax           ; }
  
  void    SetHistoEnergyRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    SetHistoPtRangeAndNBins(min, max, n) ; }
  
  Int_t   GetHistoFinePtBins()          const { return fHistoFinePtBins      ; }
  Float_t GetHistoFinePtMin()           const { return fHistoFinePtMin       ; }
  Float_t GetHistoFinePtMax()           const { return fHistoFinePtMax       ; }	
  
  void SetHistoFinePtRangeAndNBins      (Float_t min, Float_t max, Int_t n) {
    fHistoFinePtBins     = n ; fHistoFinePtMax     = max ; fHistoFinePtMin     = min ; }
  
  //Azimuthal angle
  
  Int_t   GetHistoPhiBins()             const { return fHistoPhiBins         ; }
  Float_t GetHistoPhiMin()              const { return fHistoPhiMin          ; }
  Float_t GetHistoPhiMax()              const { return fHistoPhiMax          ; }
  
  void    SetHistoPhiRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoPhiBins        = n ; fHistoPhiMax        = max ; fHistoPhiMin        = min ; }
  
 // Delta azymuthal angle
  Int_t   GetHistoDeltaPhiBins()        const { return fHistoDeltaPhiBins    ; }
  Float_t GetHistoDeltaPhiMin()         const { return fHistoDeltaPhiMin     ; }
  Float_t GetHistoDeltaPhiMax()         const { return fHistoDeltaPhiMax     ; }
  
  void    SetHistoDeltaPhiRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoDeltaPhiBins   = n ; fHistoDeltaPhiMax   = max ; fHistoDeltaPhiMin   = min ; }


// Delta eta angle

  Int_t   GetHistoDeltaEtaBins()        const { return fHistoDeltaEtaBins     ; }
  Float_t GetHistoDeltaEtaMin()         const { return fHistoDeltaEtaMin      ; }
  Float_t GetHistoDeltaEtaMax()         const { return fHistoDeltaEtaMax      ; }
 
  void    SetHistoDeltaEtaRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoDeltaEtaBins  = n ; fHistoDeltaEtaMax   = max ; fHistoDeltaEtaMin   = min ; }


  //Pseudorapidity-rapidity
  
  Int_t   GetHistoEtaBins()             const { return fHistoEtaBins         ; }
  Float_t GetHistoEtaMin()              const { return fHistoEtaMin          ; }
  Float_t GetHistoEtaMax()              const { return fHistoEtaMax          ; }
  
  void    SetHistoEtaRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoEtaBins        = n ; fHistoEtaMax         = max ; fHistoEtaMin       = min ; }
  
  //Mass
	
  Int_t   GetHistoMassBins()            const { return fHistoMassBins        ; }
  Float_t GetHistoMassMin()             const { return fHistoMassMin         ; }
  Float_t GetHistoMassMax()             const { return fHistoMassMax         ; }
	
  void    SetHistoMassRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoMassBins       = n ; fHistoMassMax       = max ; fHistoMassMin       = min ; }
  
  //Asymetry
	
  Int_t   GetHistoAsymmetryBins()        const { return fHistoAsymBins       ; }
  Float_t GetHistoAsymmetryMin()         const { return fHistoAsymMin        ; }
  Float_t GetHistoAsymmetryMax()         const { return fHistoAsymMax        ; }	
  
  void    SetHistoAsymmetryRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoAsymBins       = n ; fHistoAsymMax       = max ; fHistoAsymMin       = min ; }
  
  //VZero
  
  Int_t   GetHistoV0SignalBins()         const { return fHistoV0SBins        ; }
  Int_t   GetHistoV0SignalMin()          const { return fHistoV0SMin         ; }
  Int_t   GetHistoV0SignalMax()          const { return fHistoV0SMax         ; }
	
  void    SetHistoV0SignalRangeAndNBins(Int_t min, Int_t max, Int_t n) {
    fHistoV0SBins        = n ; fHistoV0SMax        = max ; fHistoV0SMin        = min ; }
  
  Int_t   GetHistoV0MultiplicityBins()   const { return fHistoV0MBins        ; }
  Int_t   GetHistoV0MultiplicityMin()    const { return fHistoV0MMin         ; }
  Int_t   GetHistoV0MultiplicityMax()    const { return fHistoV0MMax         ; }
  
  void    SetHistoV0MultiplicityRangeAndNBins(Int_t min, Int_t max, Int_t n) {
    fHistoV0MBins        = n ; fHistoV0MMax        = max ; fHistoV0MMin        = min ; }
  
  // Track multiplicity
  
  Int_t   GetHistoTrackMultiplicityBins() const { return fHistoTrMBins      ; }
  Int_t   GetHistoTrackMultiplicityMin()  const { return fHistoTrMMin       ; }
  Int_t   GetHistoTrackMultiplicityMax()  const { return fHistoTrMMax       ; }
  
  void    SetHistoTrackMultiplicityRangeAndNBins(Int_t min, Int_t max, Int_t n) {
    fHistoTrMBins        = n ; fHistoTrMMax        = max ; fHistoTrMMin        = min ; }
  
  // dEdx
  
  Int_t   GetHistodEdxBins()             const { return fHistodEdxBins       ; }
  Float_t GetHistodEdxMin()              const { return fHistodEdxMin        ; }
  Float_t GetHistodEdxMax()              const { return fHistodEdxMax        ; }	
  
  void    SetHistodEdxRangeAndNBins        (Float_t min, Float_t max, Int_t n) {
    fHistodEdxBins       = n ; fHistodEdxMax       = max ; fHistodEdxMin       = min ; }
  
  // E over p
  
  Int_t   GetHistoPOverEBins()           const { return fHistoPOverEBins     ; }
  Float_t GetHistoPOverEMin()            const { return fHistoPOverEMin      ; }
  Float_t GetHistoPOverEMax()            const { return fHistoPOverEMax      ; }
	
  void    SetHistoPOverERangeAndNBins      (Float_t min, Float_t max, Int_t n) {
    fHistoPOverEBins     = n ; fHistoPOverEMax     = max ; fHistoPOverEMin     = min ; }
  
  // Number of cells per clusters
  
  Int_t   GetHistoNClusterCellBins()     const { return fHistoNClusCellBins  ; }
  Int_t   GetHistoNClusterCellMin()      const { return fHistoNClusCellMin   ; }
  Int_t   GetHistoNClusterCellMax()      const { return fHistoNClusCellMax   ; }	
 
  void    SetHistoNClusterCellRangeAndNBins(Int_t   min, Int_t   max, Int_t n) {
    fHistoNClusCellBins  = n ; fHistoNClusCellMax  = max ; fHistoNClusCellMin  = min ; }
  
  // Number of clusters
  
  Int_t   GetHistoNClustersBins()        const { return fHistoNClustersBins  ; }
  Int_t   GetHistoNClustersMin()         const { return fHistoNClustersMin   ; }
  Int_t   GetHistoNClustersMax()         const { return fHistoNClustersMax   ; }	
  
  void    SetHistoNClustersRangeAndNBins   (Int_t   min, Int_t   max, Int_t n) {
    fHistoNClustersBins  = n ; fHistoNClustersMax  = max ; fHistoNClustersMin  = min ; }

  // Number of cells 
  
  Int_t   GetHistoNCellsBins()           const { return fHistoNCellsBins     ; }
  Int_t   GetHistoNCellsMin()            const { return fHistoNCellsMin      ; }
  Int_t   GetHistoNCellsMax()            const { return fHistoNCellsMax      ; }	
  
  void    SetHistoNCellsRangeAndNBins      (Int_t   min, Int_t   max, Int_t n) {
    fHistoNCellsBins     = n ; fHistoNCellsMax     = max ; fHistoNCellsMin     = min ; }

  // dR
  
  Int_t   GetHistodRBins()               const { return fHistodRBins         ; }
  Float_t GetHistodRMin()                const { return fHistodRMin          ; }
  Float_t GetHistodRMax()                const { return fHistodRMax          ; }	
  
  void    SetHistodRRangeAndNBins          (Float_t min, Float_t max, Int_t n) {
    fHistodRBins         = n ; fHistodRMax         = max ; fHistodRMin         = min ; }

  // Ratio
  
  Int_t   GetHistoRatioBins()            const { return fHistoRatioBins      ; }
  Float_t GetHistoRatioMin()             const { return fHistoRatioMin       ; }
  Float_t GetHistoRatioMax()             const { return fHistoRatioMax       ; }	
  
  void    SetHistoRatioRangeAndNBins       (Float_t min, Float_t max, Int_t n) {
    fHistoRatioBins      = n ; fHistoRatioMax      = max ; fHistoRatioMin      = min ; }

  // Vertex
  
  Int_t   GetHistoVertexDistBins()       const { return fHistoVertexDistBins ; }
  Float_t GetHistoVertexDistMin()        const { return fHistoVertexDistMin  ; }
  Float_t GetHistoVertexDistMax()        const { return fHistoVertexDistMax  ; }	
  
  void    SetHistoVertexDistRangeAndNBins  (Float_t min, Float_t max, Int_t n) { 
    fHistoVertexDistBins = n ; fHistoVertexDistMax = max ; fHistoVertexDistMin = min ; }

  
  // R =sqrt(x^2+y^2+z^2) (cm)
  
  Int_t   GetHistoRBins()                const { return fHistoRBins          ; }
  Float_t GetHistoRMin()                 const { return fHistoRMin           ; }
  Float_t GetHistoRMax()                 const { return fHistoRMax           ; }	  
  
  void    SetHistoRRangeAndNBins           (Float_t min, Float_t max, Int_t n) {
    fHistoRBins         = n ; fHistoRMax           = max ; fHistoRMin          = min ; }
  
  // X position
  
  Int_t   GetHistoXBins()                const { return fHistoXBins          ; }
  Float_t GetHistoXMin()                 const { return fHistoXMin           ; }
  Float_t GetHistoXMax()                 const { return fHistoXMax           ; }
  
  void    SetHistoXRangeAndNBins           (Float_t min, Float_t max, Int_t n) {
    fHistoXBins          = n ; fHistoXMax          = max ; fHistoXMin          = min ; }

  // Y position
  
  Int_t   GetHistoYBins()                const { return fHistoYBins          ; }
  Float_t GetHistoYMin()                 const { return fHistoYMin           ; }
  Float_t GetHistoYMax()                 const { return fHistoYMax           ; }	
  
  void    SetHistoYRangeAndNBins           (Float_t min, Float_t max, Int_t n) {
    fHistoYBins          = n ; fHistoYMax          = max ; fHistoYMin          = min ; }
  
  // Z position
  
  Int_t   GetHistoZBins()                const { return fHistoZBins          ; }
  Float_t GetHistoZMin()                 const { return fHistoZMin           ; }
  Float_t GetHistoZMax()                 const { return fHistoZMax           ; }	
	
  void    SetHistoZRangeAndNBins           (Float_t min, Float_t max, Int_t n) {
    fHistoZBins         = n ; fHistoZMax           = max ; fHistoZMin          = min ; }

  // Shower shape parameters
  
  Int_t   GetHistoShowerShapeBins()      const { return fHistoSSBins         ; }
  Float_t GetHistoShowerShapeMin()       const { return fHistoSSMin          ; }
  Float_t GetHistoShowerShapeMax()       const { return fHistoSSMax          ; }	
  
  void    SetHistoShowerShapeRangeAndNBins (Float_t min, Float_t max, Int_t n) {
    fHistoSSBins        = n ; fHistoSSMax          = max ; fHistoSSMin        = min ; }

  // Time
  
  Int_t   GetHistoTimeBins()             const { return fHistoTimeBins       ; }
  Float_t GetHistoTimeMin()              const { return fHistoTimeMin        ; }
  Float_t GetHistoTimeMax()              const { return fHistoTimeMax        ; }	
  
  void    SetHistoTimeRangeAndNBins        (Float_t min, Float_t max, Int_t n) {
    fHistoTimeBins       = n ; fHistoTimeMax       = max ; fHistoTimeMin       = min ; }	  
  
  // Cluster time difference
  
  Int_t   GetHistoDiffTimeBins()         const { return fHistoDiffTimeBins   ; }
	Float_t GetHistoDiffTimeMin()          const { return fHistoDiffTimeMin    ; }
	Float_t GetHistoDiffTimeMax()          const { return fHistoDiffTimeMax    ; }	
   
  void    SetHistoDiffTimeRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoDiffTimeBins   = n ; fHistoDiffTimeMax   = max ; fHistoDiffTimeMin   = min   ; }
  
  //    Track matching histogrammes setters and getters
  
  Int_t   GetHistoTrackResidualEtaBins() const { return fHistoTrackResidualEtaBins ; }
	Float_t GetHistoTrackResidualEtaMin()  const { return fHistoTrackResidualEtaMin  ; }
	Float_t GetHistoTrackResidualEtaMax()  const { return fHistoTrackResidualEtaMax  ; }	  
  
  void    SetHistoTrackResidualEtaRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoTrackResidualEtaBins = n ; fHistoTrackResidualEtaMax = max ; fHistoTrackResidualEtaMin = min           ; }
  
  Int_t   GetHistoTrackResidualPhiBins() const { return fHistoTrackResidualPhiBins ; }
	Float_t GetHistoTrackResidualPhiMin()  const { return fHistoTrackResidualPhiMin  ; }
	Float_t GetHistoTrackResidualPhiMax()  const { return fHistoTrackResidualPhiMax  ; }	
  
  void    SetHistoTrackResidualPhiRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoTrackResidualPhiBins = n ; fHistoTrackResidualPhiMax = max ; fHistoTrackResidualPhiMin = min           ; }
  
private:    
  
  Int_t    fHistoPtBins   ;                   // Number of bins in pt axis
  Float_t  fHistoPtMax    ;                   // Maximum value of pt histogram range
  Float_t  fHistoPtMin    ;                   // Minimum value of pt histogram range
  Int_t    fHistoPhiBins  ;                   // Number of bins in phi axis
  Float_t  fHistoPhiMax   ;                   // Maximum value of phi histogram range
  Float_t  fHistoPhiMin   ;                   // Minimum value of phi histogram range
  Int_t    fHistoEtaBins  ;                   // Number of bins in eta axis
  Float_t  fHistoEtaMax   ;                   // Maximum value of eta histogram range
  Float_t  fHistoEtaMin   ;                   // Minimum value of eta histogram range
  Int_t    fHistoDeltaPhiBins  ;              // Number of bins in delta phi axis
  Float_t  fHistoDeltaPhiMax   ;              // Maximum value of delta phi histogram range
  Float_t  fHistoDeltaPhiMin   ;              // Minimum value of delta phi histogram range
  Int_t    fHistoDeltaEtaBins  ;              // Number of bins in delta eta axis
  Float_t  fHistoDeltaEtaMax   ;              // Maximum value of delta eta histogram range
  Float_t  fHistoDeltaEtaMin   ;              // Minimum value of delta eta histogram range

  Int_t    fHistoMassBins ;                   // Number of bins in mass axis
  Float_t  fHistoMassMax  ;                   // Maximum value of mass histogram range
  Float_t  fHistoMassMin  ;                   // Minimum value of mass histogram range
  Int_t    fHistoAsymBins ;                   // Number of bins in asymmetry axis
  Float_t  fHistoAsymMax  ;                   // Maximum value of asymmetry histogram range
  Float_t  fHistoAsymMin  ;                   // Minimum value of asymmetry histogram range
  Int_t    fHistoV0SBins  ;                   // Number of bins in V0 signal axis
  Int_t    fHistoV0SMax   ;                   // Maximum value of V0 signal histogram range
  Int_t    fHistoV0SMin   ;                   // Minimum value of V0 signal histogram range
  Int_t    fHistoV0MBins  ;                   // Number of bins in V0 multiplicity axis
  Int_t    fHistoV0MMax   ;                   // Maximum value of V0 multiplicity histogram range
  Int_t    fHistoV0MMin   ;                   // Minimum value of V0 multiplicity histogram range
  Int_t    fHistoTrMBins  ;                   // Number of bins in V0 multiplicity axis
  Int_t    fHistoTrMMax   ;                   // Maximum value of track multiplicity histogram range
  Int_t    fHistoTrMMin   ;                   // Minimum value of track multiplicity histogram range
  Int_t    fHistoFinePtBins;                  // fine binning for fhAmpId histogram
  Float_t  fHistoFinePtMax;                   // maximum pt value for fhAmpId histogram
  Float_t  fHistoFinePtMin;                   // minimum pt value for fhAmpId histogram
  Int_t    fHistoPOverEBins;                  // p/E histogram number of bins
  Float_t  fHistoPOverEMax;                   // p/E maximum value
  Float_t  fHistoPOverEMin;                   // p/E minimum value
  Int_t    fHistodEdxBins;                    // dEdx histogram number of bins
  Float_t  fHistodEdxMax;                     // dEdx maximum value
  Float_t  fHistodEdxMin;                     // dEdx minimum value
  Int_t    fHistodRBins;                      // dR histogram number of bins
  Float_t  fHistodRMax;                       // dR maximum value
  Float_t  fHistodRMin;                       // dR minimum value
  Int_t    fHistoTimeBins;                    // cell time histogram number of bins
  Float_t  fHistoTimeMax;                     // cell time maximum value
  Float_t  fHistoTimeMin;                     // cell time minimum value
  Int_t    fHistoNClusCellBins;               // number of cells per cluster histogram number of bins
  Int_t    fHistoNClusCellMax;                // number of cells per cluster maximum value
  Int_t    fHistoNClusCellMin;                // number of cells per cluster minimum value
  Int_t    fHistoNCellsBins;                  // number of cells histogram number of bins
  Int_t    fHistoNCellsMax;                   // number of cells maximum value
  Int_t    fHistoNCellsMin;                   // number of cells minimum value
  Int_t    fHistoNClustersBins;               // number of clusters histogram number of bins
  Int_t    fHistoNClustersMax;                // number of clusters maximum value
  Int_t    fHistoNClustersMin;                // number of clusters minimum value  
  Int_t    fHistoRatioBins;                   // ratio histogram number of bins
  Float_t  fHistoRatioMax;                    // ratio maximum value
  Float_t  fHistoRatioMin;                    // ratio minimum value
  Int_t    fHistoVertexDistBins;              // vertex distance histogram number of bins
  Float_t  fHistoVertexDistMax;               // vertex distance maximum value
  Float_t  fHistoVertexDistMin;               // vertex distance minimum value	
  Int_t    fHistoRBins;                       // r =sqrt(x^2+y^2+z^2) (cm) position histogram number of bins
  Float_t  fHistoRMax;                        // r =sqrt(x^2+y^2+z^2) (cm)  maximum value
  Float_t  fHistoRMin;                        // r =sqrt(x^2+y^2+z^2) (cm)  minimum value	
  Int_t    fHistoXBins;                       // x (cm) position histogram number of bins
  Float_t  fHistoXMax;                        // x (cm) position maximum value
  Float_t  fHistoXMin;                        // x (cm) position minimum value
  Int_t    fHistoYBins;                       // y (cm) position histogram number of bins
  Float_t  fHistoYMax;                        // y (cm) position maximum value
  Float_t  fHistoYMin;                        // y (cm) position minimum value
  Int_t    fHistoZBins;                       // z (cm) position histogram number of bins
  Float_t  fHistoZMax;                        // z (cm) position maximum value
  Float_t  fHistoZMin;                        // z (cm) position minimum value
  Int_t    fHistoSSBins;                      // Shower Shape parameter histogram number of bins
  Float_t  fHistoSSMax;                       // Shower Shape parameter position maximum value
  Float_t  fHistoSSMin;                       // Shower Shape parameter position minimum value
  Int_t    fHistoDiffTimeBins;                // Difference cluster pair time parameter histogram number of bins
  Float_t  fHistoDiffTimeMax;                 // Difference cluster pair time parameter position maximum value
  Float_t  fHistoDiffTimeMin;                 // Difference cluster pair time parameter position minimum value  
  Int_t    fHistoTrackResidualEtaBins ;       // Number of bins in dEta (cluster-track) axis
  Float_t  fHistoTrackResidualEtaMax  ;       // Maximum value of dEta (cluster-track) histogram range
  Float_t  fHistoTrackResidualEtaMin  ;       // Minimum value of dEta (cluster-track) histogram range		
  Int_t    fHistoTrackResidualPhiBins ;       // Number of bins in dPhi axis
  Float_t  fHistoTrackResidualPhiMax  ;       // Maximum value of dPhi (cluster-track) histogram range
  Float_t  fHistoTrackResidualPhiMin  ;       // Minimum value of dPhi (cluster-track) histogram range
  
  AliHistogramRanges(              const AliHistogramRanges & h) ; // cpy ctor
  AliHistogramRanges & operator = (const AliHistogramRanges & h) ; // cpy assignment
  
  ClassDef(AliHistogramRanges,3)
} ;


#endif //ALIHISTOGRAMRANGES_H



