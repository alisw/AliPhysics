#ifndef ALIMUONQATASK_H
#define ALIMUONQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
/// An analysis task to check the MUON data in simulated data
/// This class checks out the ESD tree, providing the matching with
/// the trigger,trigger responses for Low and High Pt cuts
/// (in Single, Unlike Sign and Like Sign) and gives Pt, Y, ITS vertex
/// and multiplicity distributions. All results are in histogram form.
/// The output is a root file and eps files in MUON.tar.gz. 

//*-- Frederic Yermia
//////////////////////////////////////////////////////////////////////////////

#include <TLorentzVector.h>
#include "AliAnalysisTask.h"  

class AliESD ; 
class TH1F ;
class Tree ; 
class AliMUONQATask : public AliAnalysisTask {

public:
  AliMUONQATask(const char *name) ;
  virtual ~AliMUONQATask() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container
  
  TLorentzVector fV1; //! Lorentz Momentum
  Int_t fnTrackTrig ; //! Number of trigger matching tracks
  Int_t ftracktot   ; //! Number of ESD tracks
  Int_t fnevents    ; //! Number of events
  Int_t fSLowpt     ; //! Single trigger response
  Int_t fUSLowpt    ; //! Unlike Sign Low trigger response
  Int_t fUSHighpt   ; //! Unlike Sign High trigger response
  Int_t fLSLowpt    ; //! Like Sign Low trigger response
  Int_t fLSHighpt   ; //! Like Sign High trigger response
  Float_t fmuonMass ; //! Muon mass
  Double_t fthetaX  ; //! Angle of track at vertex in X direction (rad)
  Double_t fthetaY  ; //! Angle of track at vertex in Y direction (rad)
  Double_t fpYZ     ; //! Bending Momentum
  Double_t fPxRec1  ; //! X rec. Momentum
  Double_t fPyRec1  ; //! Y rec. Momentum
  Double_t fPzRec1  ; //! Z rec. Momentum
  Double_t fE1      ; //! Muon Enenrgy
  Int_t fZ1         ; //! Z coordinate (cm)

  // Histograms
  TH1F * fhMUONVertex ; //! ITS Z-Vertex
  TH1F * fhMUONMult   ; //! Track Multilicity
  TH1F * fhPt   ;       //! Track transverse momentum
  TH1F * fhY   ;        //!  Track rapidity
  TH1F * fheffMatchT ;  //! Efficiency of trigger matching
  TH1F * fhSLowpt ;     //! Percent of single response
  TH1F * fhUSLowpt ;    //! Percent of US Low response
  TH1F * fhUSHighpt ;   //! Percent of US High response
  TH1F * fhLSLowpt ;    //! Percent of LS Low response
  TH1F * fhLSHighpt ;   //! Percent of LS High response
  TH1F * fhChi2   ;      //! Track Chi Square by d.o.f.
  TH1F * fhChi2match ;  //! Chi2 of trigger/track matching 

  AliMUONQATask(const AliMUONQATask&); // Not implemented
  AliMUONQATask& operator=(const AliMUONQATask&); // Not implemented

  ClassDef(AliMUONQATask, 0); // a MUON photon analysis task 
};
#endif // ALIMUONQATASK_H
