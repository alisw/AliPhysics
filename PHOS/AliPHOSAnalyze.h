#ifndef ALIPHOSANALYZE_H
#define ALIPHOSANALYZE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Algorythm class to analyze PHOS events    //
//  Yves Schutz            SUBATECH           //
//                                            //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSv0.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSPID.h"

class AliPHOSAnalyze : public TObject {

public:

  AliPHOSAnalyze() ;              // ctor
  AliPHOSAnalyze(Text_t * name) ; // ctor
  virtual ~AliPHOSAnalyze() ;     // dtor

  void AnalyzeOneEvent(Int_t evt = -999) ;  // analyzes a single event ;
  void AnalyzeManyEvents(Int_t Nevtents = 100, Int_t Module=0) ;  // analyzes many events   ;
  void BookingHistograms() ;                // booking histograms for the ManyEvent analysis ;
  Bool_t Init(Int_t evt) ;                  // does various initialisations
  void DisplayKineEvent(Int_t evt = -999) ; // displays the Kine events in ALICE coordinate 
  void DisplayRecParticles() ;              // displays RecParticles in ALICE coordinate  
  void DisplayRecPoints() ;                 // displays RecPoints in module coordinate  
  void DisplayTrackSegments() ;             // displays TrackSegments in module coordinate  
  Bool_t OpenRootFile(Text_t * name) ;      // opens the root file
  void SavingHistograms() ;                 // Save histograms in a root file
 
private:
  
  AliPHOSClusterizer * fClu ;       // a clusterizer 
  Int_t fEvt ;                      // the evt number being processed 
  AliPHOSGeometry * fGeom;          // the PHOS Geometry object
  AliPHOSv0 * fPHOS ;               // the PHOS object from the root file 
  AliPHOSPID * fPID ;               // a particle identifier
  AliPHOSReconstructioner * fRec ;  // a reconstructioner  
  TFile * fRootFile ;               // the root file that contains the data
  AliPHOSTrackSegmentMaker * fTrs ; // a tracksegmentmaker ;
  TH1F * fhEmcDigit        ;   // Histo of digit energies in the Emc 
  TH1F * fhVetoDigit       ;   // Histo of digit energies in the Veto 
  TH1F * fhConvertorDigit  ;   // Histo of digit energies in the Convertor
  TH1F * fhEmcCluster      ;   // Histo of Cluster energies in Emc
  TH1F * fhVetoCluster     ;   // Histo of Cluster energies in Veto
  TH1F * fhConvertorCluster;   // Histo of Cluster energies in Convertor
  TH2F * fhConvertorEmc    ;   // 2d Convertor versus Emc energies
  TH1F * fhPhotonEnergy    ;   // Spectrum of detected photons
  TH1F * fhElectronEnergy  ;   // Spectrum of detected electrons
  TH1F * fhNeutralHadronEnergy   ;   // Spectrum of detected neutral hadron
  TH1F * fhNeutralEMEnergy   ;   // Spectrum of detected neutral EM
  TH1F * fhChargedHadronEnergy   ;   // Spectrum of detected charged
  TH1F * fhPhotonPositionX ;   // X distribution of detected photons
  TH1F * fhElectronPositionX ; // X distribution of detected electrons
  TH1F * fhNeutralHadronPositionX  ; // X distribution of detected neutral hadron
  TH1F * fhNeutralEMPositionX  ; // X distribution of detected neutral EM
  TH1F * fhChargedHadronPositionX  ; // X distribution of detected charged
  TH1F * fhPhotonPositionY   ; // Y distribution of detected photons
  TH1F * fhElectronPositionY ; // Y distribution of detected electrons
  TH1F * fhNeutralHadronPositionY  ; // Y distribution of detected neutral hadron
  TH1F * fhNeutralEMPositionY  ; // Y distribution of detected neutral EM
  TH1F * fhChargedHadronPositionY  ; // Y distribution of detected charged


ClassDef(AliPHOSAnalyze,1)  // PHOS event analyzis , version 1

};

#endif // AliPHOSANALYZE_H
