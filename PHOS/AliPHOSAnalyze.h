#ifndef ALIPHOSANALYZE_H
#define ALIPHOSANALYZE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Algorithm class to analyze PHOSv1 events    
//*-- Author : Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSv1.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSPID.h"
#include "AliPHOSIndexToObject.h"

class AliPHOSAnalyze : public TObject {

public:

  AliPHOSAnalyze() ;              // ctor
  AliPHOSAnalyze(Text_t * name) ; // ctor
  AliPHOSAnalyze(const AliPHOSAnalyze & ana) ; // cpy ctor                   
  virtual ~AliPHOSAnalyze() ;     // dtor

  void ActivePPSD(Int_t Nevents) ;
  void AnalyzeManyEvents(Int_t Nevtents = 100, Int_t Module=0) ;  // analyzes many events   ;
  void ReconstructCPV(Int_t Nevents = 1) ;  // reconstruct points in EMC and CPV
  void AnalyzeResolutions(Int_t Nevtents) ; // analyzes Energy and Position resolutions   ;
  void ReadAndPrintCPV(Int_t Nevents);      // Read & print generated and reconstructed hits in CPV
  void AnalyzeCPV(Int_t Nevents);           // analyzes various CPV characteristics
  void InvariantMass(Int_t Nevents = 100) ; 
  void Reconstruct(Int_t Nevtents = 100,Int_t FirstEvent = 0) ;
  void AnalyzeResolutions(Int_t Nevtents) ;  // analyzes Energy and Position resolutions   ;
  void BookingHistograms() ;                // booking histograms for the ManyEvent analysis ;
  void BookResolutionHistograms() ;         // booking histograms for the Resoluion analysis ;
  void Copy(TObject & obj) ;                // copies an analysis into an other one   
  Float_t CorrectEnergy(Float_t ERecPart) ;   //Corrects reconstracted energy
  Bool_t Init(Int_t evt) ;                  // does various initialisations
  void DisplayKineEvent(Int_t evt = -999) ; // displays the Kine events in ALICE coordinate 
  void DisplayRecParticles() ;              // displays RecParticles in ALICE coordinate  
  void DisplayRecPoints() ;                 // displays RecPoints in module coordinate  
  void DisplayTrackSegments() ;             // displays TrackSegments in module coordinate  
  Bool_t OpenRootFile(Text_t * name) ;      // opens the root file
  void SaveHistograms() ;                   // Save histograms in a root file
  void ResetHistograms() ;                  // 
  AliPHOSAnalyze & operator = (const AliPHOSAnalyze & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }
  void SetDebugLevel(Int_t flag) { fDebugLevel = flag; }
 
 private:
  
  AliPHOSClusterizer * fClu ;         // a clusterizer 
  Int_t fEvt ;                        // the evt number being processed 
  AliPHOSGeometry * fGeom ;           // the PHOS Geometry object
  AliPHOSIndexToObject * fObjGetter ; // provides methods to retrieve objects from their index in a list
  AliPHOSv1 * fPHOS ;                 // the PHOS object from the root file 
  AliPHOSPID * fPID ;                 // a particle identifier
  AliPHOSReconstructioner * fRec ;    // a reconstructioner  
  TFile * fRootFile ;                 // the root file that contains the data
  AliPHOSTrackSegmentMaker * fTrs ;   // a tracksegmentmaker ;

  Int_t fDebugLevel;                  // debug level for analysis

  TH2F * fhEnergyCorrelations ;     //Energy correlations between Eloss in Convertor and PPSD(2)


  TH1F * fhEmcDigit ;               // Histo of digit energies in the Emc 
  TH1F * fhVetoDigit ;              // Histo of digit energies in the Veto 
  TH1F * fhConvertorDigit ;         // Histo of digit energies in the Convertor
  TH1F * fhEmcCluster ;             // Histo of Cluster energies in Emc
  TH1F * fhVetoCluster ;            // Histo of Cluster energies in Veto
  TH1F * fhConvertorCluster ;       // Histo of Cluster energies in Convertor
  TH2F * fhConvertorEmc ;           // 2d Convertor versus Emc energies

  TH2F * fhAllEnergy ;           // Spectrum of detected photons with photon primary
  TH2F * fhPhotEnergy ;        // Total spectrum of detected photons
  TH2F * fhEMEnergy ;        // Spectrum of detected neutral EM with EM primary
  TH2F * fhPPSDEnergy ;

  TH2F * fhAllPosition ;        // Position Resolution of  photons with photon primary
  TH2F * fhPhotPosition ;     // Position Resolution of  photons
  TH2F * fhEMPosition ;     // Position Resolution of neutral EM with EM primary
  TH2F * fhPPSDPosition ;  // Position Resolution of neutral EM

  TH1F * fhPhotonPositionY ;        // Y distribution of detected photons
  TH1F * fhElectronPositionY ;      // Y distribution of detected electrons
  TH1F * fhNeutralHadronPositionY ; // Y distribution of detected neutral hadron
  TH1F * fhNeutralEMPositionY ;     // Y distribution of detected neutral EM
  TH1F * fhChargedHadronPositionY ; // Y distribution of detected charged
  TH1F * fhPhotonHadronPositionY ;  // Y distribution of detected Photon-Hadron

  TH1F * fhPhotReg ;          
  TH1F * fhAllReg ;          
  TH1F * fhNReg ;          
  TH1F * fhNBarReg ;          
  TH1F * fhChargedReg ;          
  TH1F * fhPhotEM ;          
  TH1F * fhAllEM ;          
  TH1F * fhNEM ;          
  TH1F * fhNBarEM ;          
  TH1F * fhChargedEM ;          
  TH1F * fhPhotPPSD ;          
  TH1F * fhAllPPSD ;          
  TH1F * fhNPPSD ;          
  TH1F * fhNBarPPSD ;          
  TH1F * fhChargedPPSD ;          

  TH1F * fhPrimary ;          
  TH1F * fhAllRP ;
  TH1F * fhPPSD ;
  TH1F * fhShape ;
  TH1F * fhVeto ;

  TH1F * fhPhotPhot ;
  TH1F * fhPhotElec ;
  TH1F * fhPhotNeuH ;
  TH1F * fhPhotNuEM ; 
  TH1F * fhPhotChHa ;
  TH1F * fhPhotGaHa ;



ClassDef(AliPHOSAnalyze,1)  // PHOSv1 event analyzis algorithm

};

#endif // AliPHOSANALYZE_H
