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
    kND  = 25
  };
  enum EResultContainer { // results container type
    kCenter = 0   // cluster2center pad calibration
   ,kQRes   = 1   // resolution on charge dependence
   ,kSigm   = 2   // sigma cluster as func of x and z
   ,kMean   = 3   // shift cluster as func of x and z
   ,kNtasks = 4   // total number os sub tasks
  };
  enum ECheckBits { // force setting the ExB
    kSaveAs    = BIT(22)
   ,kExB       = BIT(23)
  };
  AliTRDclusterResolution();
  AliTRDclusterResolution(const char *name);
  virtual ~AliTRDclusterResolution();

  void          UserCreateOutputObjects();
  void          UserExec(Option_t *);
  Int_t         GetDetector() const { return fDet; }
  inline Float_t GetExB() const;
  inline Float_t GetVdrift() const;
  Bool_t        GetRefFigure(Int_t ifig);
  Bool_t        HasProcess(EResultContainer bit) const {return TESTBIT(fStatus, bit);}
  Bool_t        HasExB() const { return TestBit(kExB);}

  TObjArray*    Histos(); 
  TObjArray*    Results() const {return fResults;}; 

  Bool_t        IsVisual() const {return Bool_t(fCanvas);}
  Bool_t        IsSaveAs() const {return TestBit(kSaveAs);}

  Bool_t        PostProcess();
  Bool_t        SetExB(Int_t det=-1, Int_t c = 70, Int_t r = 7);
  void          SetVisual();
  void          SetProcess(EResultContainer bit, Bool_t v = kTRUE) {v ? SETBIT(fStatus, bit) : CLRBIT(fStatus, bit);}
  void          SetSaveAs(Bool_t v = kTRUE) {SetBit(kSaveAs, v);}
  inline void   ResetProcesses();

protected:
  void    ProcessCharge();
  void    ProcessCenterPad();
  void    ProcessSigma();
  void    ProcessMean();

private:
  AliTRDclusterResolution(const AliTRDclusterResolution&);  
  AliTRDclusterResolution& operator=(const AliTRDclusterResolution&);

  TCanvas    *fCanvas; //! visualization canvas 
  TObjArray  *fInfo;   //! list of cluster info
  TObjArray  *fResults;// list of result graphs/histos/trees
  UChar_t    fStatus;  // steer parameter of the task
  Short_t    fDet;     // detector (-1 for all)
  Float_t    fExB;     // tg of the Lorentz angle
  Float_t    fVdrift;  // mean drift velocity
  Float_t    fT0;      // time 0
  static const Float_t fgkTimeBinLength;// time bin length (invers of sampling frequency)

  // working named variables
  UChar_t    fLy;      // TRD plane 
  Float_t    fT;       // calibrated time 
  Float_t    fX;       // local drift length 
  Float_t    fY;       // local rphi offset 
  Float_t    fZ;       // local anode wire offset 
  Float_t    fR[4];    // mean/sgm resolution
  Float_t    fP[4];    // mean/sgm pulls
  
  ClassDef(AliTRDclusterResolution, 3)  // cluster resolution
};

//___________________________________________________
inline Float_t AliTRDclusterResolution::GetExB() const
{ 
  if(!HasExB()){
    printf("WARNING :: ExB was not set. Use B=0.\n");
  }
  return fExB;
}

//___________________________________________________
inline Float_t AliTRDclusterResolution::GetVdrift() const
{ 
  if(!HasExB()){
    printf("WARNING :: ExB was not set. Use B=0.\n");
  }
  return fVdrift;
}

//___________________________________________________
inline void AliTRDclusterResolution::ResetProcesses()
{
  CLRBIT(fStatus, kQRes);
  CLRBIT(fStatus, kCenter);
  CLRBIT(fStatus, kSigm);
  CLRBIT(fStatus, kMean);
}

#endif

