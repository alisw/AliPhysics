#ifndef ALIITSDETTYPE_H
#define ALIITSDETTYPE_H


#include <TString.h>
#include <TObject.h>

#include "AliITSsegmentation.h"
#include "AliITSresponse.h"


class AliITSClusterFinder;
class AliITSsimulation;

class AliITSDetType:public TObject
{

 public:
    AliITSDetType();
    virtual ~AliITSDetType();
    AliITSDetType(const AliITSDetType &source); // copy constructor
    AliITSDetType& operator=(const AliITSDetType &source); // assign. operator

// Set the defaults
  virtual void   Init() {}

//
  virtual void    SegmentationModel(AliITSsegmentation* thisSegmentation){ 
    // Configure segmentation model
    if(fSegmentation) delete fSegmentation;
    fSegmentation=thisSegmentation;
  }
  //
  virtual void    ResponseModel(AliITSresponse* thisResponse) { 
    // Configure response model
    if(fResponse) delete fResponse;
    fResponse=thisResponse;
  }
  //
  virtual void    SimulationModel(AliITSsimulation *thisSimulation) {
    // Configure simulation model
    fSimulation = thisSimulation;
  }
  //
  virtual void    ReconstructionModel(AliITSClusterFinder *thisReconstruction) {
// Configure reconstruction model
      fReconst = thisReconstruction;
  }
  virtual void    ClassNames(const char *digit, const char *cluster) { 
    // Set class names for digits and clusters
    fDigClassName=digit; fClustClassName=cluster; 
  } 
  
  AliITSsegmentation*      &GetSegmentationModel(){
    //  Get reference to segmentation model
    return fSegmentation;
  }
  AliITSresponse*          &GetResponseModel(){
    //  Get reference to response model
    return fResponse;
  }
  AliITSsimulation*        &GetSimulationModel(){
    //  Get reference to simulation model
    return fSimulation;
  }
  AliITSClusterFinder*     &GetReconstructionModel(){
    //  Get reference to hit reconstruction model
    return fReconst;
  }
  //  
  
  void GetClassNames(char *digit,char *cluster) 
    { 
      // Get class names for digits and rec points
      strcpy(digit,fDigClassName.Data()); strcpy(cluster,fClustClassName.Data()); 
    } 
  
protected:
  
  AliITSClusterFinder  *fReconst;          // cluster finder
  AliITSsimulation     *fSimulation;       // simulation
  AliITSresponse       *fResponse;         // response
  AliITSsegmentation   *fSegmentation;     // segmentation
  
  TString              fDigClassName;      // string
  TString              fClustClassName;    // string
  
  ClassDef(AliITSDetType,1)
    
};

#endif
