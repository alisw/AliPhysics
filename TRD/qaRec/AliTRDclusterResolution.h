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
  enum EResultContainer { // results container type
    kQRes   = 0
   ,kCenter= 1
   ,kSigm  = 2
   ,kMean  = 3
  };
  enum ECheckBits { // force setting the ExB
    kSaveAs    = BIT(22)
   ,kExB       = BIT(23)
  };
  AliTRDclusterResolution(const char *name="ClErrParam");
  virtual ~AliTRDclusterResolution();

  void          ConnectInputData(Option_t *);
  void          CreateOutputObjects();
  void          Exec(Option_t *);
  Int_t         GetDetector() const { return fDet; }
  inline Float_t GetExB() const;
  inline Float_t GetVdrift() const;
  Bool_t        GetRefFigure(Int_t ifig);
  Bool_t        HasProcess(EResultContainer bit) const {return TESTBIT(fStatus, bit);}
  Bool_t        HasExB() const { return TestBit(kExB);}

  TObjArray*    Histos(); 

  Bool_t        IsVisual() const {return Bool_t(fCanvas);}
  Bool_t        IsSaveAs() const {return TestBit(kSaveAs);}

  Bool_t        PostProcess();
  Bool_t        SetExB(Int_t det=-1, Int_t c = 70, Int_t r = 7);
  void          SetVisual();
  void          SetProcess(EResultContainer bit, Bool_t v = kTRUE) {v ? SETBIT(fStatus, bit) : CLRBIT(fStatus, bit);}
  void          SetSaveAs(Bool_t v = kTRUE) {SetBit(kSaveAs, v);}
  inline void   ResetProcess();

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
  static const Float_t fgkTimeBinLength;// time bin length (invers of sampling frequency)

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

//___________________________________________________
inline void AliTRDclusterResolution::ResetProcess()
{
  CLRBIT(fStatus, kQRes);
  CLRBIT(fStatus, kCenter);
  CLRBIT(fStatus, kSigm);
  CLRBIT(fStatus, kMean);
}

#endif

