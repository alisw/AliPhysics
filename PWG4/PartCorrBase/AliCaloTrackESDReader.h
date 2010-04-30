#ifndef ALICALOTRACKESDREADER_H
#define ALICALOTRACKESDREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Class for reading data (ESDs) in order to do prompt gamma 
//  or other particle identification and correlations
// 
// It is a filtering class, transforms ESD tracks or CaloClusters
// into AOD tracks and calocluters, which are the basic input of the analysis
// classes in this frame.
// It is recommended to use the official filter AliAnalysisTaskESDfilter, and 
// then the reader for AODs AliCaloTrackAODReader.//
//
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system --- 

// --- AliRoot system ---
#include "AliCaloTrackReader.h" 
#include "AliESDEvent.h"

class AliCaloTrackESDReader : public AliCaloTrackReader {
  
  public: 
  
  AliCaloTrackESDReader() ; // ctor
  AliCaloTrackESDReader(const AliCaloTrackESDReader & g) ; // cpy ctor
  //AliCaloTrackESDReader & operator = (const AliCaloTrackESDReader & g) ;//cpy assignment
  virtual ~AliCaloTrackESDReader() {;} //virtual dtor
  
  void FillInputCTS  () ;
  void FillInputEMCAL() ;
  void FillInputPHOS () ;  
  void FillInputEMCALCells() ;
  void FillInputPHOSCells() ;
  
  void GetVertex(Double_t v[3]) const ;
  Double_t GetBField() const;
  void SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ; 
	
  TString GetFiredTriggerClasses() {return ((AliESDEvent*)GetInputEvent())->GetFiredTriggerClasses();}

  // Get calorimeter (Super)module number where this cluster falled
  Int_t GetModuleNumber(AliESDCaloCluster * cluster) const;
	
  ClassDef(AliCaloTrackESDReader,1)
    } ;


#endif //ALICALOTRACKESDREADER_H



