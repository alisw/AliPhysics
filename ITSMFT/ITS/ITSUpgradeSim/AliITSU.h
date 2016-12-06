#ifndef ALIITSU_H
#define ALIITSU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliITSU.h */

////////////////////////////////////////////////////////////////////////
//           Manager class for set: ITS Upgrade                       //
////////////////////////////////////////////////////////////////////////


#include <TObjArray.h> // used in inline function GetChip.
#include "AliDetector.h"
#include "AliITSUGeomTGeo.h"

class TString;
class TTree;
class AliITSMFTSDigit;
class AliITSMFTSimulation;
class AliITSMFTSegmentationPix;
class AliITSUChip;
class AliITSCalibration;
class AliITSMFTHit;
class AliITSMFTDigitPix;
class AliDigitizationInput;
class AliITSMFTSensMap;
class AliITSMFTSimuParam;
class AliITSMFTParamList;

class AliITSU : public AliDetector {

 public:
  //
  // number detector types
  enum {kNChipTypes = AliITSMFTAux::kNChipTypes};
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
  //RS  AliITSgeom* GetITSgeom() const {return fChipTypeSim->GetITSgeom();}
  //RS  void   SetITSgeom(AliITSgeom *geom) {fChipTypeSim->SetITSgeom(geom);}
  // return pointer to the array of chips
  
  AliITSUChip   * GetChip(Int_t index) {return (AliITSUChip*)fChipHits->UncheckedAt(index);}
  AliITSMFTSimuParam* GetSimuParam() const {return fSimuParam;}
    
  //================ Necessary general Classes =======================
  virtual void Init();
  virtual AliLoader* MakeLoader(const char* topfoldername);
  virtual void MakeBranch(Option_t *opt=" ");
  virtual void MakeBranchS(const char* fl);
  virtual void MakeBranchD(const char* file);
  virtual void MakeBranchInTreeD(TTree* treeD, const char* file=0);
  virtual void SetTreeAddress();
  virtual AliITSMFTSimulation*   GetSimulationModel(Int_t lr)   {return fSimModelLr[lr];}
  virtual AliITSMFTSegmentationPix*  GetSegmentation(Int_t lr)  {return fSegModelLr[lr];}
  virtual AliITSMFTParamList*    GetResponseParam(Int_t lr)     {return fResponseLr[lr];}
  //=================== Hits =========================================
  virtual void StepManager() {} // See Step Manager for specific geometry.
  //------------ sort hits by chip for Digitisation ----------------
  virtual void FillChips(Int_t bgrev, Option_t *opt, const char *filename); 
  virtual void FillChips(TTree *treeH, Int_t mask = 0);
  virtual void ClearChips();
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
  virtual void AddSumDigit(AliITSMFTSDigit &sdig);
  virtual void AddSimDigit(Int_t branch, AliITSMFTDigitPix *d);
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
  const char* GetDigitClassName(Int_t i) {return Form("AliITSMFTDigit%s",AliITSUGeomTGeo::GetChipTypeName(i));}
  const char* GetChipTypeName(Int_t i) {return AliITSUGeomTGeo::GetChipTypeName(i);}
  
 protected:
  //================== Data Members ==================================
  Bool_t                fEuclidOut;      // Flag to write geometry in euclid format
  Int_t                 fNLayers;        // the number of layers
  Int_t                *fIdSens;         //[fNLayers] layer identifier
  TString              *fLayerName;      //[fNLayers] layer identifier
  Bool_t                fTiming;         // flag to turn on/off timers.
  AliITSUGeomTGeo*      fGeomTGeo;       //  access to geometry details
  AliITSMFTSimuParam*     fSimuParam;      //!simulation parameters
  TClonesArray**        fModA;           //! Used by Raw2SDigits (one TC per chip)
  TClonesArray*         fpSDigits;       //! Branch address to build SD from raw data 
  TClonesArray*         fSDigits;        //! Branch address to build SD
  TClonesArray*         fDetHits;        //! array of full detector hits
  TObjArray*            fChipHits;       //! chip's hits container in (pointers on the fDetHits)
  TObjArray*            fDetDigits;      //! AliDetector has TClonesArray fDigits, avoid same name
  AliITSMFTSensMap*       fSensMap;        //! sensor map for digitization
  //
  AliITSMFTSimulation    **fSimModelLr;     //! simulation objects per layer
  AliITSMFTSegmentationPix   **fSegModelLr;     //! segmentation objects per layar
  AliITSMFTParamList     **fResponseLr;     //! response parameters for each layer
  TObjArray            *fCalibration;    //! calibration objects
  Int_t                 fRunNumber;      //! run number
  Bool_t                fSimInitDone;    //! flag initialized simulation
  Bool_t                fUseALPIDESim;   //! use the ALPIDE chip simulation

 private:
  AliITSU(const AliITSU &source); // copy constructor. Not to be used!
  AliITSU& operator=(const AliITSU &source); // = operator. Not to be used!
  ClassDef(AliITSU,2) // Base class for ITS
};

#endif
