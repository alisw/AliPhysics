#ifndef ALIITSSIMULATION_H
#define ALIITSSIMULATION_H

#include <TObject.h>

class AliITSresponse;
class AliITSsegmentation;
class AliITSmodule;

//___________________________________________________

class AliITSsimulation : public TObject {

public:

  AliITSsimulation();
  virtual ~AliITSsimulation() {
    // destructor
  }
  AliITSsimulation(const AliITSsimulation &source); // copy constructor
  AliITSsimulation& operator=(const AliITSsimulation &source); // ass. 

  virtual void DigitiseModule(AliITSmodule *mod,Int_t module,Int_t event) {
    // digitize module
  }

  virtual void CreateFastRecPoints(AliITSmodule *mod) {
    // create fast rec points
  }

protected:

  AliITSresponse      *fResponse;       // response
  AliITSsegmentation  *fSegmentation;   // segmentation

  ClassDef(AliITSsimulation,1)  // Simulation base class 
    
};


#endif
