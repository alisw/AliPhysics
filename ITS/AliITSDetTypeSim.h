#ifndef ALIITSDETTYPESIM_H
#define ALIITSDETTYPESIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id$ 
*/

/////////////////////////////////////////////////////////////////////////
// * This class contains all of the "external" information needed to do//
// * detector specific simulations for the ITS.                        //
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "AliITSLoader.h"
#include "AliITSSimuParam.h"
#include "AliITSFOGeneratorSPD.h"

class TObjArray;
class TClonesArray;
class TTree;
class AliCDBMetaData;
class AliITSdigit;
class AliITSmodule;
class AliITSpListItem;
class AliITSsimulation;
class AliITSsegmentation;
class AliITSresponse;
class AliITSCalibrationSSD;
class AliITSGainSSDv2;
class AliITSBadChannelsSSDv2;
class AliITSNoiseSSDv2;
class AliITSresponseSDD;
class AliITSDDLModuleMapSDD;
class AliITSCalibration;
class AliITSgeom;
class AliITSFOSignalsSPD;
class AliITSTriggerConditions;

class AliITSDetTypeSim : public TObject {
 public:
  
    AliITSDetTypeSim();
    virtual ~AliITSDetTypeSim(); 
    AliITSgeom *GetITSgeom() const {
        if(fLoader)return ((AliITSLoader*)fLoader)->GetITSgeom();
	else return 0;}
    void SetITSgeom(AliITSgeom *geom);
    
    virtual void SetSimulationModel(Int_t dettype,AliITSsimulation *sim);
    virtual AliITSsimulation* GetSimulationModel(Int_t dettype) const;        
    virtual AliITSsimulation* GetSimulationModelByModule(Int_t module) const;
        
    virtual void SetSegmentationModel(Int_t dettype,AliITSsegmentation *seg);
    virtual AliITSsegmentation* GetSegmentationModel(Int_t dettype) const;
    virtual AliITSsegmentation* GetSegmentationModelByModule(Int_t module) const;

    virtual void SetCalibrationModel(Int_t iMod,AliITSCalibration *resp);
    virtual void SetSPDNoisyModel(Int_t iMod, AliITSCalibration *cal);
    virtual void SetSPDSparseDeadModel(Int_t iMod, AliITSCalibration *cal);

    virtual AliITSCalibration* GetCalibrationModel(Int_t iMod) const;
    virtual AliITSCalibration* GetSPDNoisyModel(Int_t iMod) const;
    virtual AliITSTriggerConditions* GetTriggerConditions();

    virtual void SetSimuParam(const AliITSSimuParam* spar){
      if(fSimuPar) delete fSimuPar;
      fSimuPar = new AliITSSimuParam(*spar);
    }
    virtual AliITSSimuParam* GetSimuParam() const {return fSimuPar;}

    virtual AliITSDDLModuleMapSDD* GetDDLModuleMapSDD()const { return fDDLMapSDD;}
    virtual AliITSresponseSDD* GetResponseSDD() const { return fRespSDD;}
    TObjArray* GetCalibrationArray() const {return fCalibration;}
    TObjArray* GetSegmentation() const {return fSegmentation;}
    void ResetCalibrationArray();
    void ResetSegmentation();

    virtual void SetLoader(AliITSLoader* loader);
    AliITSLoader* GetLoader() const {return fLoader;}

    virtual void SetDefaults();
    virtual void SetDefaultSimulation();
    virtual void SetRunNumber(Int_t rn=0){fRunNumber = rn;}
    virtual Int_t GetRunNumber() const {return fRunNumber;}
    virtual void SetTreeAddressS(TTree* treeS, const Char_t* name);
    virtual void SetTreeAddressD(TTree* treeD, const Char_t* name);

    virtual void SetDigits(TObjArray* digits) {fDigits=digits;}
    const TClonesArray* GetSDigits() const { return &fSDigits;}
    TObjArray*    GetDigits() const {return fDigits;}
    Int_t* GetNDigitArray() const {return fNDigits;}
    TClonesArray *DigitsAddress(Int_t id) const {
	return ((TClonesArray*)(*fDigits)[id]);}
    virtual void ResetSDigits(){fNSDigits=0; fSDigits.Clear();}
    virtual void ResetDigits();
    virtual void ResetDigits(Int_t branch);
    virtual void SDigitsToDigits(Option_t *opt, Char_t* name);

    virtual void AddSumDigit(AliITSpListItem &sdig);
    virtual void AddSimDigit(Int_t branch, const AliITSdigit *d);
    virtual void AddSimDigit(Int_t branch,Float_t phys,Int_t* digits,
			     Int_t* tracks,Int_t *hits,Float_t* trkcharges,
			     Int_t sigexpanded=-1000);
    virtual void SetDigitClassName(Int_t i, const Char_t* name) {
	fkDigClassName[i]=name;}
    const Char_t* GetDigitClassName(Int_t i) const {return fkDigClassName[i];}

    virtual void ResetFOSignals() {fFOGenerator.ResetSignals();}
    virtual void ProcessSPDDigitForFastOr(UInt_t module, UInt_t colM, UInt_t rowM);
    virtual void ProcessNoiseForFastOr() {fFOGenerator.ProcessNoise();}
    virtual AliITSFOSignalsSPD* GetFOSignals() {return fFOGenerator.GetFOSignals();}
    virtual void WriteFOSignals();
    virtual Float_t GetAverageGainSDD() const {
      if(fAveGainSDD>0.) return fAveGainSDD;
      else return 1.;
    }


 protected:
    virtual void CreateCalibrationArray(); 
    virtual Bool_t GetCalibration();
    
 private:
    AliITSDetTypeSim(const AliITSDetTypeSim &source);
    AliITSDetTypeSim& operator=(const AliITSDetTypeSim &source);
    void SetDefaultSegmentation(Int_t idet);  // creates def segm.

    //conversion of the old SSD calibration objects tothe new ones
    void ReadOldSSDNoise(const TObjArray *array, 
			 AliITSNoiseSSDv2 *noiseSSD);
    void ReadOldSSDBadChannels(const TObjArray *array, 
			       AliITSBadChannelsSSDv2 *badChannelsSSD);
    void ReadOldSSDGain(const TObjArray *array, 
			AliITSGainSSDv2 *gainSSD);

    static const Int_t fgkNdettypes;          // number of different det. types
    static const Int_t fgkDefaultNModulesSPD; // Total numbers of SPD modules by default
    static const Int_t fgkDefaultNModulesSDD; // Total numbers of SDD modules by default
    static const Int_t fgkDefaultNModulesSSD; // Total numbers of SSD modules by default

    TObjArray    *fSimulation;   //! [NDet]
    TObjArray    *fSegmentation; //! [NDet]
    TObjArray    *fCalibration;  //! [NMod]
    AliITSCalibrationSSD* fSSDCalibration;  //! SSD calibration object
    TObjArray    *fSPDNoisy;     //! [fgkDefaultNModulesSPD]
    TObjArray    *fSPDSparseDead;//! [fgkDefaultNModulesSPD]
    Int_t         fNSDigits;     //! number of SDigits
    TClonesArray  fSDigits;      //! Summable digits
    Int_t*        fNDigits;      //! [NDet] number of Digits for det.
    Int_t         fRunNumber;    //! run number (to access DB)
    TObjArray     *fDigits;      //! [NMod][NDigits]
    AliITSSimuParam *fSimuPar;   //! detector simulation parameters
    AliITSDDLModuleMapSDD *fDDLMapSDD; //! mapping DDL/module -> SDD module number
    AliITSresponseSDD *fRespSDD;  //! SDD response parameters 
    Float_t       fAveGainSDD;    //! Average gain of SDD good anodes
    const Char_t*       fkDigClassName[3]; //! String with digit class name.
    AliITSLoader* fLoader;          //! loader  
    Bool_t        fFirstcall;       //! flag
    AliITSFOGeneratorSPD fFOGenerator; //! Fast-OR generator object
    AliITSTriggerConditions* fTriggerConditions; //! Trigger conditions 
       
    ClassDef(AliITSDetTypeSim,15) // ITS Simulation structure
 
};

#endif
