#ifndef ITSv5_H
#define ITSv5_H
////////////////////////////////////////////////////////////////////////
//     Manager and hits classes for  ITS version 5
////////////////////////////////////////////////////////////////////////

#include "TString.h"

#include "AliITS.h"
#include "AliITSgeom.h"

class AliITSv5 : public AliITS {
private:
    Int_t fId5N; // The number of layers for geometry version 5
    // The name of the layers as defined in the Geant tree.
    char  **fId5Name;

public:
                         AliITSv5();
			 AliITSv5(const char *name, const char *title);
           virtual       ~AliITSv5() {}
           virtual void  CreateGeometry();
           virtual void  CreateMaterials();
           virtual void  Init();   
    inline virtual Int_t IsVersion() const {return 5;}
           virtual void  StepManager();
  
  ClassDef(AliITSv5,1)  //Hits manager for ITS version 5
};
 
#endif
