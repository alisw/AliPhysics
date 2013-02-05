#ifndef ALIITSU_H
#define ALIITSU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliITSU.h */

////////////////////////////////////////////////////////////////////////
//           Manager class for set: ITS Upgrade                       //
////////////////////////////////////////////////////////////////////////


#include <TObjArray.h> // used in inline function GetModule.
#include "AliDetector.h"
#include "AliITSTrigger.h"
#include "AliITSDetTypeSim.h"
#include "AliITSUGeomTGeo.h"

class TString;
class TTree;
class AliITSUSDigit;
class AliITSUSimulation;
class AliITSsegmentation;
class AliITSUModule;
class AliITSCalibration;
class AliITSUHit;
class AliITSdigit;
class AliDigitizationInput;
class AliITSUSensMap;
class AliITSUSimuParam;
class AliParamList;

class AliITSU : public AliDetector {

 public:
  //
  // number detector types
  enum {kNDetTypes = AliITSUGeomTGeo::kNDetTypes};
  //
  //================= Standard Classes ===============================
  AliITSU();  // Default creator.
  AliITSU(const char *title, Int_t nlayers); // extended standard Creator
  virtual ~AliITSU(); // destructor
  virtual Int_t IsVersion() const {return 11;}
  
  //===================== Simulation Geometry ========================
  // get geometry version - detailed (major) or coarse (minor)
  virtual Int_t GetMajorVersion() const {return -1;}
  virtual Int_t GetMinorVersion() const {return -1;}
  virtual void  GetGeometryVersion(Int_t &a,Int_t &b) const {a = GetMajorVersion();b=GetMinorVersion();return;}
  virtual void  SetEUCLID(Bool_t euclid=kTRUE) {fEuclidOut = euclid;}
  virtual Bool_t GetEUCLID()const {return fEuclidOut;}
  //-------------------- Geometry Transformations --------------------
  
  // ITS geometry functions From Simulation
  AliITSUGeomTGeo* GetITSGeomTGeo() const {return fGeomTGeo;}
  //RS  AliITSgeom* GetITSgeom() const {return fDetTypeSim->GetITSgeom();}
  //RS  void   SetITSgeom(AliITSgeom *geom) {fDetTypeSim->SetITSgeom(geom);}
  // return pointer to the array of modules
  
  AliITSUModule   * GetModule(Int_t index) {return (AliITSUModule*)fModuleHits->UncheckedAt(index);}
  AliITSUSimuParam* GetSimuParam() const {return fSimuParam;}
    
  //================ Necessary general Classes =======================
  virtual void Init();
  virtual AliLoader* MakeLoader(const char* topfoldername);
  virtual void MakeBranch(Option_t *opt=" ");
  virtual void MakeBranchS(const char* fl);
  virtual void MakeBranchD(const char* file);
  virtual void MakeBranchInTreeD(TTree* treeD, const char* file=0);
  virtual void SetTreeAddress();
  virtual AliITSUSimulation*   GetSimulationModel(Int_t lr)   {return (AliITSUSimulation*)fSimModelLr[lr];}
  virtual AliITSsegmentation*  GetSegmentation(Int_t lr)      {return (AliITSsegmentation*)fSegModelLr[lr];}
  virtual AliParamList*        GetResponseParam(Int_t lr)     {return (AliParamList*)fResponseLr[lr];}
  //=================== Hits =========================================
  virtual void StepManager() {} // See Step Manager for specific geometry.
  //------------ sort hits by module for Digitisation ----------------
  virtual void FillModules(Int_t bgrev, Option_t *opt, const char *filename); 
  virtual void FillModules(TTree *treeH, Int_t mask = 0);
  virtual void ClearModules();
  virtual void AddHit(Int_t track, Int_t *vol, Float_t *hits);
  void         InitSimulation();
  //
  // Trigger
  //  virtual AliTriggerDetector* CreateTriggerDetector() const;

  AliDigitizer* CreateDigitizer(AliDigitizationInput* manager) const;
  virtual void SDigits2Digits();
  virtual void Hits2Digits(); 
  virtual void Hits2SDigits();
  virtual void Hits2Digits(Int_t evNumber,Int_t bgrev,Option_t *option, const char *filename);
  virtual void Hits2SDigits(Int_t evNumber,Int_t bgrev,Option_t *option,const char *filename);
    
  virtual void ResetSDigits()       {if (fSDigits) fSDigits->Clear();}
  virtual void ResetDigits();
  virtual void ResetDigits(Int_t branch);
  virtual void AddSumDigit(AliITSUSDigit &sdig);
  virtual void AddSimDigit(Int_t branch, AliITSdigit *d);
  virtual void AddSimDigit(Int_t branch,Float_t phys,Int_t* digits,Int_t* tracks,Int_t *hits,Float_t* trkcharges,Int_t sigexpanded=-1000);
  TObjArray*   GetDigits()                const {return fDetDigits;}
  TClonesArray *DigitsAddress(Int_t id)  {return fDetDigits ? (TClonesArray*)fDetDigits->At(id) : 0;}
  //Fast simulation
  virtual void  Hits2FastRecPoints(Int_t bgrev,Option_t *opr, const char *filename);
  virtual Int_t Hits2Clusters(TTree *in, TTree *out);
  virtual void  CheckLabels(Int_t lab[3]) const;

  //===================== Raw Data IO ================================
  // Write digits into raw data format
  virtual void   Digits2Raw();
  virtual Bool_t Raw2SDigits(AliRawReader*);
    
  //===================== FO signals ================================
  // Write FO signals in UserInfo of SDigits/Digits tree
  void    WriteFOSignals();

  void    SetRunNumber(Int_t rn=0)        {fRunNumber = rn;}
  Int_t   GetNLayers()              const {return fNLayers;}
  Int_t   GetRunNumber()            const {return fRunNumber;}
  Bool_t  IsSimInitDone()           const {return fSimInitDone;}

 protected:
  void        InitArrays();
  const char* GetDigitClassName(Int_t i) {return Form("AliITSUDigit%s",AliITSUGeomTGeo::GetDetTypeName(i));}
  const char* GetDetTypeName(Int_t i) {return AliITSUGeomTGeo::GetDetTypeName(i);}
  
 protected:
  //================== Data Members ==================================
  Bool_t                fEuclidOut;      // Flag to write geometry in euclid format
  Int_t                 fNLayers;        // the number of layers
  Int_t                *fIdSens;         //[fNLayers] layer identifier
  TString              *fLayerName;      //[fNLayers] layer identifier
  Bool_t                fTiming;         // flag to turn on/off timers.
  AliITSUGeomTGeo*      fGeomTGeo;       //  access to geometry details
  AliITSUSimuParam*     fSimuParam;      //!simulation parameters
  TClonesArray**        fModA;           //! Used by Raw2SDigits (one TC per module)
  TClonesArray*         fpSDigits;       //! Branch address to build SD from raw data 
  TClonesArray*         fSDigits;        //! Branch address to build SD
  TClonesArray*         fDetHits;        //! array of full detector hits
  TObjArray*            fModuleHits;     //! module's hits container in (pointers on the fDetHits)
  TObjArray*            fDetDigits;      //! AliDetector has TClonesArray fDigits, avoid same name
  AliITSUSensMap*       fSensMap;        //! sensor map for digitization
  //
  AliITSUSimulation    **fSimModelLr;     //! simulation objects per layer
  AliITSsegmentation   **fSegModelLr;     //! segmentation objects per layar
  AliParamList         **fResponseLr;     //! response parameters for each layer
  TObjArray            *fCalibration;    //! calibration objects
  Int_t                 fRunNumber;      //! run number
  Bool_t                fSimInitDone;    //! flag initialized simulation

 private:
  AliITSU(const AliITSU &source); // copy constructor. Not to be used!
  AliITSU& operator=(const AliITSU &source); // = operator. Not to be used!
  ClassDef(AliITSU,1) // Base class for ITS
};

#endif
