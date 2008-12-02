#ifndef ALIANAEXAMPLE_H
#define ALIANAEXAMPLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Example class on how to read AODCaloClusters, ESDCaloCells and AODTracks and how 
// fill AODs with PWG4PartCorr analysis frame
// Select the type of detector information that you want to analyze, CTS (tracking), PHOS or EMCAL
// Select the PID custer type of the calorimeters
// Set min momentum of the cluster/tracks
// Fill few histograms
//
//-- Author: Gustavo Conesa (INFN-LNF)

#include "TH2F.h" 
#include "AliAnaPartCorrBaseClass.h"
 
class AliAnaExample : public AliAnaPartCorrBaseClass {
       
  public: 
       
       AliAnaExample() ; // default ctor
       AliAnaExample(const AliAnaExample & g) ; // cpy ctor
       AliAnaExample & operator = (const AliAnaExample & g) ;//cpy assignment
       virtual ~AliAnaExample() {;} //virtual dtor
              
       TList * GetCreateOutputObjects();
       
       void InitParameters();
       
       void Print(const Option_t * opt) const;
       
       void MakeAnalysisFillAOD()  ;
       
       void MakeAnalysisFillHistograms() ; 
 
       Int_t GetPdg() const {return fPdg ;}
       void SetPdg( Int_t pdg ) {fPdg = pdg; }

       TString GetDetector() const {return fDetector ;}
       void SetDetector( TString calo ) {fDetector = calo; }

 private:
       
       Int_t  fPdg ; //identified particle id
       TString fDetector ; //detector selection
       //Histograms
       //CaloClusters 
       TH1F * fhPt; //! pT distribution
       TH1F * fhPhi; //! phi distribution
       TH1F * fhEta; //! eta distribution
       TH2F * fh2Pt; //!pT distribution, reconstructed vs generated
       TH2F * fh2Phi; //! phi distribution, reconstructed vs generated
       TH2F * fh2Eta; //! eta distribution, reconstructed vs generated
       //Calo Cells
       TH1F * fhNCells; //! Number of towers/crystals with signal
       TH1F * fhAmplitude; //! Amplitude measured in towers/crystals
       ClassDef(AliAnaExample,1)
 } ;


#endif //ALIANAEXAMPLE_H



