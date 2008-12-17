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
  enum { // force setting the ExB
    kExB       = BIT(23)
  };
  AliTRDclusterResolution();
  virtual ~AliTRDclusterResolution();

  void    ConnectInputData(Option_t *);
  void    CreateOutputObjects();
  void    Exec(Option_t *);
  Int_t   GetDetector() const { return fDet; }
  inline Float_t GetExB() const;
  Bool_t  GetRefFigure(Int_t ifig);
  Bool_t  HasExB() const { return TestBit(kExB);}
  TObjArray*  Histos(); 
  
  Bool_t  PostProcess();
  Bool_t  SetExB(Int_t det=-1);
  void    Terminate(Option_t *){};

private:
  AliTRDclusterResolution(const AliTRDclusterResolution&);  
  AliTRDclusterResolution& operator=(const AliTRDclusterResolution&);

  TObjArray  *fInfo;   // list of cluster info
  TObjArray  *fResults;// list of result graphs/histos
  TAxis      *fAt;     // binning in the x(radial) direction (time)
  TAxis      *fAd;     // binning in the z direction (drift cell)
  Float_t    fExB;     // tg of the Lorentz angle
  Short_t    fDet;     // detector (-1 for all)

  ClassDef(AliTRDclusterResolution, 0)  // cluster resolution
};

inline Float_t AliTRDclusterResolution::GetExB() const
{ 
  if(!HasExB()){
    printf("WARNING :: ExB was not set. Use B=0.\n");
  }
  return fExB;
}
#endif

