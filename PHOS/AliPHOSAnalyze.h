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

  void AnalyzeOneEvent(Int_t evt = -999) ;  // analyzes a single event ;
  void AnalyzeManyEvents(Int_t Nevtents = 100, Int_t Module=0) ;  // analyzes many events   ;
  void AnalyzeResolutions(Int_t Nevtents) ;  // analyzes Energy and Position resolutions   ;
  void BookingHistograms() ;                // booking histograms for the ManyEvent analysis ;
  void BookResolutionHistograms() ;         // booking histograms for the Resoluion analysis ;
  void Copy(TObject & obj) ;                // copies an analysis into an other one   
  Bool_t Init(Int_t evt) ;                  // does various initialisations
  void DisplayKineEvent(Int_t evt = -999) ; // displays the Kine events in ALICE coordinate 
  void DisplayRecParticles() ;              // displays RecParticles in ALICE coordinate  
  void DisplayRecPoints() ;                 // displays RecPoints in module coordinate  
  void DisplayTrackSegments() ;             // displays TrackSegments in module coordinate  
  Bool_t OpenRootFile(Text_t * name) ;      // opens the root file
  void SavingHistograms() ;                 // Save histograms in a root file
  void SaveResolutionHistograms() ;         // Save histograms in a root file

  AliPHOSAnalyze & operator = (const AliPHOSAnalyze & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }
 
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
  TH1F * fhEmcDigit ;               // Histo of digit energies in the Emc 
  TH1F * fhVetoDigit ;              // Histo of digit energies in the Veto 
  TH1F * fhConvertorDigit ;         // Histo of digit energies in the Convertor
  TH1F * fhEmcCluster ;             // Histo of Cluster energies in Emc
  TH1F * fhVetoCluster ;            // Histo of Cluster energies in Veto
  TH1F * fhConvertorCluster ;       // Histo of Cluster energies in Convertor
  TH1F * fhConvertorEmc ;           // 2d Convertor versus Emc energies
  TH2F * fhPhotonEnergy ;           // Spectrum of detected photons
  TH2F * fhElectronEnergy ;         // Spectrum of detected electrons
  TH2F * fhNeutralHadronEnergy ;    // Spectrum of detected neutral hadron
  TH2F * fhNeutralEMEnergy ;        // Spectrum of detected neutral EM
  TH2F * fhChargedHadronEnergy ;    // Spectrum of detected charged
  TH2F * fhPhotonHadronEnergy ;     // Spectrum of detected Photon-Hadron
  TH2F * fhPhotonPosition ;        // Position Resolution of  photons
  TH2F * fhElectronPosition ;      // Position Resolution of electrons
  TH2F * fhNeutralHadronPosition ; // Position Resolution of neutral hadron
  TH2F * fhNeutralEMPosition ;     // Position Resolution of neutral EM
  TH2F * fhChargedHadronPosition ; // Position Resolution of charged
  TH2F * fhPhotonHadronPosition ;  // Position Resolution of Photon-Hadron
  TH1F * fhPhotonPositionY ;        // Y distribution of detected photons
  TH1F * fhElectronPositionY ;      // Y distribution of detected electrons
  TH1F * fhNeutralHadronPositionY ; // Y distribution of detected neutral hadron
  TH1F * fhNeutralEMPositionY ;     // Y distribution of detected neutral EM
  TH1F * fhChargedHadronPositionY ; // Y distribution of detected charged
  TH1F * fhPhotonHadronPositionY ;  // Y distribution of detected Photon-Hadron


ClassDef(AliPHOSAnalyze,1)  // PHOSv1 event analyzis algorithm

};

#endif // AliPHOSANALYZE_H
