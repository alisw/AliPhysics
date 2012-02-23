/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnaConvCorrPion.h
/// @author Svein Lindal
/// @brief  Class used to find correlations between pions and charged tracks


#ifndef ALIANACONVCORRPION_H
#define ALIANACONVCORRPION_H

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

void CreateHistograms();

  void FillTriggerCounters(const AliAODConversionParticle * particle, Int_t isolated);
  
 private:

  void InitMassAxis();

  TAxis fAxisM;  //Mass axis

  AliAnaConvCorrPion(const AliAnaConvCorrPion&); // not implemented
  AliAnaConvCorrPion& operator=(const AliAnaConvCorrPion&); // not implemented
  ClassDef(AliAnaConvCorrPion, 3); //

};

#endif
