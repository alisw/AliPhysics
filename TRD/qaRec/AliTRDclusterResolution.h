#ifndef ALITRDCLUSTERRESOLUTION_H
#define ALITRDCLUSTERRESOLUTION_H


#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TCanvas;
class TObjArray;
class TAxis;
class AliTRDclusterResolution : public AliTRDrecoTask
{
public:
  enum EAxisBinning { // bins in z and x direction
    kNTB = 25
    ,kND = 5
    ,kN  = kND*kNTB
  };
  enum EResultContainers { // results containers
    kQRes   = 0
    ,kCenter= 1
    ,kSigm  = 2
    ,kMean  = 3
  };
  enum ECheckBits { // force setting the ExB
    kExB       = BIT(23)
  };
  enum ESteeringBits { // steering bits for task
    kSaveAs         = 0
    ,kProcCharge    = 1
    ,kProcCenterPad = 2
    ,kProcSigma     = 3
    ,kProcMean      = 4
  };
  AliTRDclusterResolution();
  virtual ~AliTRDclusterResolution();

  void    ConnectInputData(Option_t *);
  void    CreateOutputObjects();
  void    Exec(Option_t *);
  Int_t   GetDetector() const { return fDet; }
  inline Float_t GetExB() const;
  inline Float_t GetVdrift() const;
  Bool_t  GetRefFigure(Int_t ifig);
  Bool_t  HasProcessCharge() const {return TESTBIT(fStatus, kProcCharge);}
  Bool_t  HasProcessCenterPad() const {return TESTBIT(fStatus, kProcCenterPad);}
  Bool_t  HasExB() const { return TestBit(kExB);}
  Bool_t  HasProcessMean() const {return TESTBIT(fStatus, kProcMean);}
  Bool_t  HasProcessSigma() const {return TESTBIT(fStatus, kProcSigma);}

  TObjArray*  Histos(); 

  Bool_t  IsVisual() const {return Bool_t(fCanvas);}
  Bool_t  IsSaveAs() const {return TESTBIT(fStatus, kSaveAs);}

  Bool_t  PostProcess();
  Bool_t  SetExB(Int_t det=-1, Int_t c = 70, Int_t r = 7);
  void    SetVisual();
  void    SetProcessCharge(Bool_t v = kTRUE) {v ? SETBIT(fStatus, kProcCharge) : CLRBIT(fStatus, kProcCharge);}
  void    SetProcessCenterPad(Bool_t v = kTRUE) {v ? SETBIT(fStatus, kProcCenterPad) : CLRBIT(fStatus, kProcCenterPad);}
  void    SetProcessMean(Bool_t v = kTRUE) {v ? SETBIT(fStatus, kProcMean) : CLRBIT(fStatus, kProcMean);}
  void    SetProcessSigma(Bool_t v = kTRUE) {v ? SETBIT(fStatus, kProcSigma) : CLRBIT(fStatus, kProcSigma);}
  void    SetSaveAs(Bool_t v = kTRUE) {v ? SETBIT(fStatus, kSaveAs) : CLRBIT(fStatus, kSaveAs);}

  void    Terminate(Option_t *){};

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
  TObjArray  *fResults;// list of result graphs/histos
  TAxis      *fAt;     //! binning in the x(radial) direction (time)
  TAxis      *fAd;     //! binning in the z direction (drift cell)
  UChar_t    fStatus;  // steer parameter of the task
  Short_t    fDet;     // detector (-1 for all)
  Float_t    fExB;     // tg of the Lorentz angle
  Float_t    fVdrift;  // mean drift velocity

  ClassDef(AliTRDclusterResolution, 1)  // cluster resolution
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
#endif

