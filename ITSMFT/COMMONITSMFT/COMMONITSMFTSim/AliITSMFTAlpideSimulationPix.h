#ifndef ALIITSMFTALPIDESIMULATIONPIX_H
#define ALIITSMFTALPIDESIMULATIONPIX_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////
// Simulation class for Alpide upgrade pixels (2016)      //
//                                                        //
// Author: D. Pagano                                      //
// Contact: davide.pagano@cern.ch                         //
////////////////////////////////////////////////////////////


#include "AliITSMFTSimulation.h"

class TLorentzVector;
class AliITSMFTSimuParam;

//-------------------------------------------------------------------

class AliITSMFTAlpideSimulationPix : public AliITSMFTSimulation {
public:
    AliITSMFTAlpideSimulationPix();
    AliITSMFTAlpideSimulationPix(AliITSMFTSimuParam*, AliITSMFTSensMap*);
    virtual   ~AliITSMFTAlpideSimulationPix();

    void      Init();
    void      SDigitiseChip(TClonesArray*);
    void      FinishSDigitiseChip(TObjArray*);
    void      DigitiseChip(TObjArray*);
    void      SetResponseParam(AliITSMFTParamList*);
    Bool_t    AddSDigitsToChip(TSeqCollection*, Int_t);
    void      GenerateCluster();

private:
    void      FrompListToDigits(TObjArray*);
    void      WriteSDigits(TClonesArray*);
    Double_t  ACSFromBetaGamma(Double_t, Double_t) const; // Returns the average cluster size from the betagamma value
    Int_t     CSSampleFromLandau(Double_t, Double_t) const; // Sample the actual cluster size from a Landau distribution
    Double_t  ComputeIncidenceAngle(TLorentzVector) const; // Compute the angle between the particle and the normal to the chip
    Int_t     GetPixelPositionResponse(Int_t, Int_t, Float_t, Float_t, Double_t) const;
    void      CreateDigi(UInt_t, UInt_t, Int_t, Int_t);

protected:
    ClassDef(AliITSMFTAlpideSimulationPix,1)  // Simulation of pixel clusters
};
#endif
