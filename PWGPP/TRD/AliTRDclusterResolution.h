#ifndef ALITRDCLUSTERRESOLUTION_H
#define ALITRDCLUSTERRESOLUTION_H

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD cluster error parameterization                  
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TCanvas;
class TObjArray;
class AliTRDclusterResolution : public AliTRDrecoTask
{
public:
  enum EAxisBinning { // bins in z and x direction
    kND  = 1
  };
  enum EResultContainer { // results container type
    kYSys = 0   // cluster2center pad calibration
   ,kYRes   = 1   // resolution on charge dependence
   ,kSigm   = 2   // sigma cluster as func of x and z
   ,kMean   = 3   // shift cluster as func of x and z
   ,kNtasks = 4   // total number of subtasks
  };
  enum ECalibrationParam { // calibration parameters to be used from OCDB
    kVdrift = 0
   ,kT0     = 1
   ,kGain   = 2
  };
  enum ECheckBits {
    kSaveAs     = BIT(21) // save intermediary results
   ,kCalibrated = BIT(22) // load calibration
   ,kGlobal     = BIT(23) // load global position
  };
  AliTRDclusterResolution();
  AliTRDclusterResolution(const char *name);
  virtual ~AliTRDclusterResolution();

  void          UserCreateOutputObjects();
  void          UserExec(Option_t *);
  Int_t         GetDetector() const { return fDet; }
  void          GetPad(Int_t &c, Int_t &r) const { c=fCol, r=fRow; return;}
  inline void   GetDiffCoeff(Float_t &dt, Float_t &dl) const;
  inline Float_t GetExB() const;
  inline Float_t GetVdrift() const;
  inline Float_t GetT0() const;
  inline Float_t GetGain() const;
  Float_t       GetDyRange() const {return fDyRange;}
  Bool_t        GetRefFigure(Int_t ifig);
  Bool_t        HasProcess(EResultContainer bit) const {return TESTBIT(fSubTaskMap, bit);}
  Bool_t        IsCalibrated() const { return TestBit(kCalibrated);}
  Bool_t        HasGlobalPosition() const { return TestBit(kGlobal);}
  Bool_t        IsUsingCalibParam(ECalibrationParam par) const {return TESTBIT(fUseCalib, par);}

  TObjArray*    Histos(); 
  TObjArray*    Results() const {return fResults;}; 

  Bool_t        IsVisual() const {return Bool_t(fCanvas);}
  Bool_t        IsSaveAs() const {return TestBit(kSaveAs);}

  Bool_t        LoadCalibration();
  Bool_t        LoadGlobalChamberPosition();
  Bool_t        PostProcess();
  void          SetCalibrationRegion(Int_t det, Int_t col=-1, Int_t row=-1);
  void          SetVisual();
  void          SetDyRange(Float_t dy) {if(dy>0) fDyRange = dy;}
  void          SetProcess(EResultContainer bit, Bool_t v = kTRUE) {v ? SETBIT(fSubTaskMap, bit) : CLRBIT(fSubTaskMap, bit);}
  void          SetSaveAs(Bool_t v = kTRUE) {SetBit(kSaveAs, v);}
  void          SetUseCalibParam(ECalibrationParam par, Bool_t v = kTRUE) {v ? SETBIT(fUseCalib, par) : CLRBIT(fUseCalib, par);}
  inline void   ResetProcesses();

protected:
  void    ProcessCharge();
  Bool_t  ProcessNormalTracks();
  void    ProcessSigma();
  void    ProcessMean();

private:
  AliTRDclusterResolution(const AliTRDclusterResolution&);  
  AliTRDclusterResolution& operator=(const AliTRDclusterResolution&);

  TCanvas    *fCanvas; //! visualization canvas 
  TObjArray  *fInfo;   //! list of cluster info
  TObjArray  *fResults;// list of result graphs/histos/trees
  UChar_t    fSubTaskMap;  // steer map for subtasks
  UChar_t    fUseCalib;    // steer map for calibration params
  Short_t    fDet;     // detector (-1 for all)
  Char_t     fCol;     // pad column (-1 for all)
  Char_t     fRow;     // pad row (-1 for all)
  Float_t    fExB;     // tg of the Lorentz angle
  Float_t    fDt;      // diffusion coeff. transversal
  Float_t    fDl;      // diffusion coeff. longitudinal
  Float_t    fVdrift;  // mean drift velocity
  Float_t    fT0;      // time 0
  Float_t    fGain;    // gain
  Float_t    fXch;     // anode wire position for chamber
  Float_t    fZch;     // Z position for calibration element
  Float_t    fH;       // tg of tilting angle
  static const Float_t fgkTimeBinLength;// time bin length (invers of sampling frequency)

  // working named variables
  Float_t    fDyRange; // min/max dy
  UChar_t    fLy;      // TRD plane 
  Float_t    fT;       // calibrated time 
  Float_t    fX;       // local drift length 
  Float_t    fY;       // local rphi offset 
  Float_t    fZ;       // local anode wire offset 
  Float_t    fR[4];    // mean/sgm resolution
  Float_t    fP[4];    // mean/sgm pulls
  
  ClassDef(AliTRDclusterResolution, 6)  // cluster resolution
};

//___________________________________________________
inline void AliTRDclusterResolution::GetDiffCoeff(Float_t &dt, Float_t &dl) const
{
  if(!IsCalibrated()) printf(" - W - AliTRDclusterResolution::GetDiffCoeff() : Instance not calibrated.\n");
  dt=fDt; dl=fDl;
  return;
}


//___________________________________________________
inline Float_t AliTRDclusterResolution::GetExB() const
{ 
  if(!IsCalibrated()) printf(" - W - AliTRDclusterResolution::GetExB() : Instance not calibrated.\n");
  return fExB;
}

//___________________________________________________
inline Float_t AliTRDclusterResolution::GetVdrift() const
{ 
  if(!IsCalibrated()) printf(" - W - AliTRDclusterResolution::GetVdrift() : Instance not calibrated.\n");
  return fVdrift;
}

//___________________________________________________
inline Float_t AliTRDclusterResolution::GetT0() const
{
  if(!IsCalibrated()) printf(" - W - AliTRDclusterResolution::GetT0() : Instance not calibrated.\n");
  return fT0;
}

//___________________________________________________
inline Float_t AliTRDclusterResolution::GetGain() const
{
  if(!IsCalibrated()) printf(" - W - AliTRDclusterResolution::GetGain() : Instance not calibrated.\n");
  return fGain;
}

//___________________________________________________
inline void AliTRDclusterResolution::ResetProcesses()
{
  CLRBIT(fSubTaskMap, kYRes);
  CLRBIT(fSubTaskMap, kYSys);
  CLRBIT(fSubTaskMap, kSigm);
  CLRBIT(fSubTaskMap, kMean);
}

#endif

