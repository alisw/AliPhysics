#ifndef ALIITSSIMULATION_H
#define ALIITSSIMULATION_H
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include <TObject.h>
#include "AliITSDetTypeSim.h"
#include "AliITSpList.h"

class AliITSCalibration;
class AliITSsegmentation;
class AliITSmodule;
class TRandom;
class TClonesArray;

// This is the base class for ITS detector signal simulations. Data members
// include are a pointer to the detectors specific response and segmentation
// classes. See the detector specific implementations for the propper code.

class AliITSsimulation : public TObject {

  public:
    AliITSsimulation(); // Default constructor
    // Standard constructor
    AliITSsimulation(AliITSDetTypeSim *dettyp);
    virtual ~AliITSsimulation(); // destructor
    // copy constructor. See detector specific implementation.
    AliITSsimulation(const AliITSsimulation &source);
    // Assignment opporator. See detector specific implementation.
    virtual AliITSsimulation& operator=(const AliITSsimulation &source);
    // Initialize simulation
    virtual void Init() {};

    // *****************  Hits -> SDigits ******************
    // digitize module using the "slow" detector simulator creating
    // summable digits.
    virtual void SDigitiseModule(AliITSmodule *,Int_t,Int_t){;}

    // ***************** sum of SDigits -> Digits **********
    // Reset module arrays (maps), etc
    virtual void InitSimulationModule(Int_t,Int_t){;}
    // add (sum) a new list of summable digits to module, 
    // add an offset (mask) to the track numbers. Returns kTRUE if there
    // is a "good" signal in this module.
    virtual Bool_t AddSDigitsToModule( TClonesArray *pItemArray, Int_t mask );
    // digitize module using the "slow" detector simulator from
    // the sum of summable digits.
    virtual void FinishSDigitiseModule(){;}

    // **************** Hits -> Digits *********************
    // digitize module using the "slow" detector simulator creating digits.
    virtual void DigitiseModule(AliITSmodule *,Int_t,Int_t) {;}
    // digitizes module using the "fast" detector simulator.
    virtual void CreateFastRecPoints(AliITSmodule *,Int_t,
				     TRandom *,TClonesArray* /*recp*/) {;}
   // Return pointer to Response model
    virtual AliITSCalibration* GetCalibrationModel(Int_t mod = 0){return fDetType->GetCalibrationModel(mod);}
   // set pointer to Response model
    virtual void SetCalibrationModel(Int_t mod, AliITSCalibration *res){fDetType->SetCalibrationModel(mod,res);}
    // Return pointer to Segmentation object
    virtual AliITSsegmentation* GetSegmentationModel(Int_t dt) = 0;
    // set pointer to Segmentation object
    virtual void SetSegmentationModel(Int_t dt,AliITSsegmentation *seg) = 0;
    virtual AliITSpList* GetMap(){return fpList;} // Returns fpList, the map.
    virtual void SetMap(AliITSpList *p){fpList = p;} // Sets fpList, the map.
    virtual void ClearMap(){fpList->ClearMap();} // Clear fpList, map.
    virtual void SetModuleNumber(Int_t mod){fModule=mod;} // Set Module number
    virtual Int_t GetModuleNumber()const {return fModule;}// Gets Module number
    virtual void SetEventNumber(Int_t evnt){fEvent=evnt;} // Set Event number
    virtual Int_t GetEventNumber()const {return fEvent;}// Gets Event number
    // Sets the debug flag for debugging output
    void SetDebug(Int_t level=5){fDebug=level;}
    // Clears the debug flag so no debugging output will be generated
    void SetNoDebug(){fDebug=0;}
    // Returns the debug flag value
    Bool_t GetDebug(Int_t level=1)const {return fDebug>=level;}
    void SetDetType(AliITSDetTypeSim* dettyp) {fDetType=dettyp;}

 protected:
    AliITSDetTypeSim    *fDetType;        //! Access resp and segm via this obj
    AliITSpList         *fpList;          //!
    Int_t                fModule;         //!
    Int_t                fEvent;          //!
    Int_t                fDebug;          //  debug flag

  ClassDef(AliITSsimulation,4)  // Simulation base class 
    
};

#endif
