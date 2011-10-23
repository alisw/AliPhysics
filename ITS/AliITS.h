#ifndef ALIITS_H
#define ALIITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//           Manager class for set: ITS                               //
////////////////////////////////////////////////////////////////////////


#include <TObjArray.h> // used in inline function GetModule.
#include "AliDetector.h"
#include "AliITSTrigger.h"
#include "AliITSDetTypeSim.h"


class TString;
class TTree;
class AliITSpListItem;
class AliITSsimulation;
class AliITSsegmentation;

class AliITSCalibration;
class AliITShit;
class AliITSgeom;
class AliITSdigit;
class AliITSmodule;
class AliDigitizationInput;
class TArrayI;

class AliITS : public AliDetector {

 public:
  enum {kSPD,kSDD,kSSD};
    //================= Standard Classes ===============================
    AliITS();  // Default creator.
    AliITS(const Char_t *title); // standard Creator
    AliITS(const char *name, const char *title); // extended standard Creator
    virtual ~AliITS(); // destructor
    virtual Int_t IsVersion() const {return 1;}

    //===================== Simulation Geometry ========================
    // get geometry version - detailed (major) or coarse (minor)
    virtual Int_t GetMajorVersion() const {return -1;}
    virtual Int_t GetMinorVersion() const {return -1;}
    virtual void  GetGeometryVersion(Int_t &a,Int_t &b) const
	                   {a = GetMajorVersion();b=GetMinorVersion();return;}
    virtual void  SetEUCLID(Bool_t euclid=kTRUE) {fEuclidOut = euclid;}
    virtual Bool_t GetEUCLID()const {return fEuclidOut;}
    //-------------------- Geometry Transformations --------------------

    // ITS geometry functions From Simulation
    AliITSgeom* GetITSgeom() const {return fDetTypeSim->GetITSgeom();}
    void   SetITSgeom(AliITSgeom *geom) {fDetTypeSim->SetITSgeom(geom);}
    // return pointer to the array of modules
    TObjArray *GetModules(){return fITSmodules;}

    AliITSmodule *GetModule(Int_t index){
        return (AliITSmodule*)(fITSmodules->At(index));}
    virtual void SetSimuParam(AliITSSimuParam *sp){
      fSimuParam=sp;
      fDetTypeSim->SetSimuParam(sp);
    }
    AliITSSimuParam* GetSimuParam() const {return fSimuParam;}

    virtual void SetDetTypeSim(AliITSDetTypeSim* dts) {fDetTypeSim=dts;}
    AliITSDetTypeSim* GetDetTypeSim() const {return fDetTypeSim;}
    //================ Necessary general Classes =======================
    virtual void Init();
    virtual AliLoader* MakeLoader(const char* topfoldername);
    virtual void SetDefaults();
    virtual void SetDefaultSimulation();
    virtual void MakeBranch(Option_t *opt=" ");
    virtual void MakeBranchS(const char* fl);
    virtual void MakeBranchD(const char* file);
    virtual void MakeBranchInTreeD(TTree* treeD, const char* file=0);
    virtual void SetTreeAddress();
     //---------- Configuration Methods (per detector type) -------------
    // Determines which ITS subdetectors will be processed. Effects
    // digitization, and Reconstruction only.
    void SetDetectors(Option_t *opt="All"){fOpt = opt;}
    // Returns the list of ITS subdetectors that will be processed.
    Option_t* GetDetectors(){return fOpt;}

    // Set calibration
    virtual void SetCalibrationModel(Int_t dettype, AliITSCalibration *cal){
        fDetTypeSim->SetCalibrationModel(dettype,cal);}
    // Set segmentation for Simulation
    virtual void SetSegmentationModel(Int_t id, AliITSsegmentation *seg){
        fDetTypeSim->SetSegmentationModel(id,seg);}
    // Set simulation 
    virtual void SetSimulationModel(Int_t id, AliITSsimulation *sim){
        fDetTypeSim->SetSimulationModel(id,sim);}
    // Set simulation 
    virtual AliITSsimulation* GetSimulationModel(Int_t id){
	return fDetTypeSim->GetSimulationModel(id);}
    //=================== Hits =========================================
    virtual void StepManager() {} // See Step Manager for specific geometry.
    //------------ sort hits by module for Digitisation ----------------
    virtual void FillModules(Int_t /* evnt */,Int_t bgrev,Int_t /* nmodules */,
			     Option_t *opt, const char *filename); 
    virtual Bool_t InitModules(Int_t size,Int_t &nmodules);  
    virtual void FillModules(TTree *treeH, Int_t mask = 0);
    virtual void ClearModules(){if(fITSmodules) fITSmodules->Delete();}
    virtual void AddHit(Int_t track, Int_t *vol, Float_t *hits);
   // Trigger
    virtual AliTriggerDetector* CreateTriggerDetector() const;

    AliDigitizer* CreateDigitizer(AliDigitizationInput* digInput) const;
    virtual void UpdateInternalGeometry();
    virtual void SDigitsToDigits(Option_t *opt="All");
    virtual void SDigits2Digits(){SDigitsToDigits("All");}
    virtual void Hits2Digits(); 
    virtual void Hits2SDigits();
    virtual void Hits2PreDigits();
    virtual void HitsToDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                          Option_t *option,Option_t *opt,
                          const char *filename);
    virtual void HitsToPreDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                             Option_t *option,Option_t *opt,
                             const char *filename);
    void HitsToSDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                       Option_t *add, Option_t *det, const char *filename)
      {HitsToPreDigits(evNumber,bgrev,size,add,det,filename);}
    
    virtual void ResetDigits();
    virtual void ResetDigits(Int_t branch);
    virtual void AddSumDigit(AliITSpListItem &sdig);
    virtual void AddSimDigit(Int_t branch, AliITSdigit *d);
    virtual void AddSimDigit(Int_t branch,Float_t phys,Int_t* digits,
		     Int_t* tracks,Int_t *hits,Float_t* trkcharges,
		     Int_t sigexpanded=-1000);
    TObjArray* GetDigits()  const {return fDetTypeSim->GetDigits();}
    Int_t* GetNDigitArray() const {return fDetTypeSim->GetNDigitArray();}
    TClonesArray *DigitsAddress(Int_t id) {
	return fDetTypeSim->DigitsAddress(id);}
    //Fast simulation
    virtual void HitsToFastRecPoints(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, const char *filename);
    virtual Int_t Hits2Clusters(TTree *in, TTree *out);
    virtual void CheckLabels(Int_t lab[3]) const;

    //===================== Raw Data IO ================================
    // Write digits into raw data format
    virtual void   Digits2Raw();
    virtual Bool_t Raw2SDigits(AliRawReader*);
    
    //===================== FO signals ================================
    // Write FO signals in UserInfo of SDigits/Digits tree
    void WriteFOSignals();
    void     SetRawID2ClusID(const TArrayI* arr, Int_t iDet) { if (iDet>-1&&iDet<fgkNTYPES) fRawID2ClusID[iDet] = arr;}
    const TArrayI* GetRawID2ClusID(Int_t iDet) const {return (iDet>-1&&iDet<fgkNTYPES) ? fRawID2ClusID[iDet]:0;}

 protected:
    static const Int_t fgkNTYPES=3; //number of detector types
    //================== Data Members ==================================
    AliITSDetTypeSim *fDetTypeSim; //detector type for simulation
    Bool_t        fEuclidOut;  // Flag to write geometry in euclid format
    Option_t     *fOpt;        //! Detector option ="All" unless changed.
    Int_t         fIdN;        // the number of layers
    Int_t        *fIdSens;     //[fIdN] layer identifier
    TString      *fIdName;     //[fIdN] layer identifier
    TObjArray    *fITSmodules; //! Pointer to ITS modules
    Bool_t        fTiming;     // flag to turn on/off timers.
    AliITSSimuParam* fSimuParam; //simulation parameters
    TClonesArray** fModA;      //! Used by Raw2SDigits (one TC per module)
    TClonesArray* fpSDigits;   //! Branch address to build SD from raw data 
    const TArrayI* fRawID2ClusID[fgkNTYPES]; //! optional array for SDigit->Cluster assingment in Raw2SDigit (for embedding)
 private:
    AliITS(const AliITS &source); // copy constructor. Not to be used!
    AliITS& operator=(const AliITS &source); // = operator. Not to be used!
    ClassDef(AliITS,9) // Base class for ITS

};

#endif
