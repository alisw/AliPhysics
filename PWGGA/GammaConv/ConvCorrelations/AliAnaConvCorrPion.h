/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnaConvCorrPion.h
/// @author Svein Lindal
/// @brief  Class used to find correlations between pions and charged tracks


#ifndef ALIANACONVCORRPION_CXX
#define ALIANACONVCORRPION_CXX

#include "AliAnaConvCorrBase.h"
class TH2D;
//#include "THnSparse.h"

//class AliAODConversionPhoton;
class TClonesArray;

class AliAnaConvCorrPion : public AliAnaConvCorrBase {

public:

  AliAnaConvCorrPion(); 
  AliAnaConvCorrPion(TString name, TString title);
  virtual ~AliAnaConvCorrPion();

  TAxis& GetAxisM() { return fAxisM; }


  //Correlate pions with charged tracks
  //virtual void CorrelateWithHadrons(AliAODConversionPhoton * pion, const TClonesArray * tracks, const Bool_t isolated, const Int_t nSpawn, const Int_t * const spawn );

  void CreateHistograms();

  //void Process(TClonesArray * pions, TClonesArray * photons, TClonesArray * tracks);
  

  void FillTriggerCounters(const AliAODConversionParticle * particle, Int_t leading);
  
 private:

  void InitMassAxis();
  //Get array of track labels of the 4 decay electrons (2gamma * 2 electrons)
  //void GetTrackLabels(const AliAODConversionPhoton * pion, const TClonesArray * photons, Int_t* trackLabels);

  //TH2F * fhPtVsInvMass;

  TH2D * hTriggerPtvsMass[3]; //Histograms containing number of triggers in various bins
  TAxis fAxisM;  //Mass axis

  AliAnaConvCorrPion(const AliAnaConvCorrPion&); // not implemented
  AliAnaConvCorrPion& operator=(const AliAnaConvCorrPion&); // not implemented
  ClassDef(AliAnaConvCorrPion, 2); //

};

#endif
