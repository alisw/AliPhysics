#ifndef ALITRDCLUSTERRESOLUTION_H
#define ALITRDCLUSTERRESOLUTION_H


#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class TAxis;
class AliTRDclusterResolution : public AliTRDrecoTask
{
public:
  enum { // bins in z and x direction
    kNTB = 24
    ,kND = 5
    ,kN  = kND*kNTB
  };
  enum { // results containers
    kQRes      = 0
    ,kYRes     = 1
    ,kSXRes    = 2
    ,kSYRes    = 3
  };
  AliTRDclusterResolution();
  virtual ~AliTRDclusterResolution();

  void    ConnectInputData(Option_t *);
  void    CreateOutputObjects();
  void    Exec(Option_t *);
  Float_t GetExB() const { return fExB;}
  void    GetRefFigure(Int_t ifig);

  TObjArray*  Histos(); 
  
  Bool_t  PostProcess();
  void    SetExB(Float_t exb) {fExB = exb;}
  void    Terminate(Option_t *){};

private:
  AliTRDclusterResolution(const AliTRDclusterResolution&);  
  AliTRDclusterResolution& operator=(const AliTRDclusterResolution&);

  TObjArray  *fInfo;   // list of cluster info
  TObjArray  *fResults;// list of result graphs/histos
  TAxis      *fAt;     // binning in the x(radial) direction (time)
  TAxis      *fAd;     // binning in the z direction (drift cell)
  Float_t    fExB;     // tg of the Lorentz angle

  ClassDef(AliTRDclusterResolution, 0)  // cluster resolution
};

#endif

