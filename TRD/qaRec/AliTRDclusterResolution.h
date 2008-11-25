#ifndef ALITRDCLUSTERRESOLUTION_H
#define ALITRDCLUSTERRESOLUTION_H


#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class AliTRDclusterResolution : public AliTRDrecoTask
{
public:
  enum {
    kNTB = 24
    ,kND = 5
    ,kN  = 120
  };
  AliTRDclusterResolution();
  virtual ~AliTRDclusterResolution();

  void    ConnectInputData(Option_t *);
  void    CreateOutputObjects();
  void    GetRefFigure(Int_t ifig);
  TObjArray*  Histos(); 

  void    Exec(Option_t *);
  Bool_t  PostProcess();
  void    Terminate(Option_t *){};
private:
  AliTRDclusterResolution(const AliTRDclusterResolution&);  
  AliTRDclusterResolution& operator=(const AliTRDclusterResolution&);

  TObjArray  *fInfo;

  ClassDef(AliTRDclusterResolution, 0)  // cluster resolution
};

#endif

