#ifndef ALIITSMFTSIMULATIONPIX_H
#define ALIITSMFTSIMULATIONPIX_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////
// Simulation class for upgrade pixels                    //
////////////////////////////////////////////////////////////
#include <TRandom3.h>
#include <TObjArray.h>
#include "AliITSMFTSimulation.h"
#include "AliITSMFTSegmentationPix.h"
#include "AliITSMFTSegmentationPix.h"
#include "AliITSMFTSimuParam.h"


class TH1F;
class AliITSMFTChip;
class AliITSMFTSDigit;
class AliITSMFTParamList;
class TH2;

//-------------------------------------------------------------------

class AliITSMFTSimulationPix : public AliITSMFTSimulation {
public:
     //
    AliITSMFTSimulationPix();
    AliITSMFTSimulationPix(AliITSMFTSimuParam* sim,AliITSMFTSensMap* map);
    virtual ~AliITSMFTSimulationPix();
    AliITSMFTSimulationPix(const AliITSMFTSimulationPix &source);
    AliITSMFTSimulationPix& operator=(const AliITSMFTSimulationPix &s);

    //
    enum {kCellX1,kCellX2,kCellZ1,kCellZ2,kCellYDepth,kNDtSpread}; // data used for ch. spread integral calc.
    //
    // charge spread functions defined
    enum {kSpreadFunGauss2D                   // single gaussian in 2D, SpreadFunGauss2D
        ,kSpreadFunDoubleGauss2D            // double gaussian in 2D, SpreadFunDoubleGauss2D
        ,kSpreadFunHisto                    // use 2D histo from the object stored in fResponseParam
        ,kNSpreadFuns
    };
    // These are enums for interpretation of the entries in the AliITSMFTParamList*fResponseParam :
    // object to holding the sensor-specific response data
    // fist kParamStart entries of spread fun params are reserved for common parameters
    //___ Noisy pixel type
    enum { kNoisyPixOCDB = 9990, kNoisyPixRnd = 9991 };
    // elements of the SpreadFunGauss2D parameterization (offsetted by kParamStart)
    enum {kG1MeanX=AliITSMFTSimuParam::kParamStart,kG1SigX,kG1MeanZ,kG1SigZ,kNG1Par};
    // elements of the SpreadFunDoubleGauss2D parameterization (offsetted by kParamStart)
    enum {kG2MeanX0=AliITSMFTSimuParam::kParamStart,kG2SigX0,kG2MeanZ0,kG2SigZ0,kG2MeanX1,kG2SigX1,kG2MeanZ1,kG2SigZ1,kG2ScaleG2,kNG2Par};


    void Init();
    //
    void FinishSDigitiseChip(TObjArray *);
    void DigitiseChip(TObjArray *);
    //
    void SDigitiseChip(TClonesArray *);
    void Hits2SDigits();
    void Hits2SDigitsFast();
    void Hits2SDigitsFastDigital();
    void AddNoisyPixels();
    void RemoveDeadPixels();
    //
    Int_t AddRandomNoisePixels(Double_t tof=0);
    Bool_t SetTanLorAngle(Double_t WeightHole=1.0);
    Double_t GetTanLorAngle() const {return fTanLorAng;};
    //
    // For backwards compatibility
    //void SDigitsToDigits(){ FinishSDigitiseChip();}
    //
    Double_t SpreadFunDoubleGauss2D(const Double_t *dtIn);
    Double_t SpreadFunGauss2D(const Double_t *dtIn);
    Double_t SpreadFrom2DHisto(const Double_t *dtIn);
    //
    virtual void SetResponseParam(AliITSMFTParamList* resp);
    //
    Int_t GetReadOutCycle(Int_t row, Int_t col, Double_t hitTime);
    Int_t GetReadOutCycleRollingShutter(Int_t row, Int_t col, Double_t hitTime);
    //
    void CalcDiodeShiftInPixel(Int_t xrow, Int_t zcol, Float_t &x, Float_t &z);
    //
private:
    void WriteSDigits(TClonesArray *);
    void FrompListToDigits(TObjArray *);

    void SpreadCharge2D(Double_t x0,Double_t z0, Double_t dy, Int_t ix0,Int_t iz0,
                        Double_t el, Double_t tof, Int_t tID, Int_t hID);
    void PlaceDigitalPixels(Double_t x0,Double_t z0, Double_t el, Double_t tof, Int_t tID, Int_t hID);

    //
    void SetCoupling(AliITSMFTSDigit* old);     // "New" coupling routine  Tiziano Virgili
    void SetCouplingOld(AliITSMFTSDigit* old);  // "Old" coupling routine  Rocco Caliandro

protected:
    Double_t      fTanLorAng;               //! Tangent of the Lorentz Angle (weighted average for hole and electrons)
    Double_t      fGlobalChargeScale;       // Charge scaling to match Geant and Test beam
    //
    TH2*          fSpread2DHisto;           //! optional 2D histo for charge spread parameterization
    Double_t (AliITSMFTSimulationPix::*fSpreadFun)(const Double_t *dtIn); //! pointer on current spread function
    Int_t    (AliITSMFTSimulationPix::*fROTimeFun)(Int_t row,Int_t col, Double_t hitTime); //! pointer on current R/O time check function

    ClassDef(AliITSMFTSimulationPix,1)  // Simulation of pixel clusters
};
#endif
