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
    ~AliITSDetType(){
      //destructor
    }
    AliITSDetType(const AliITSDetType &source); // copy constructor
    AliITSDetType& operator=(const AliITSDetType &source); // assign. operator

// Set the defaults
  void   Init() {}

//
  void    SegmentationModel(AliITSsegmentation* thisSegmentation){ 
    // Configure segmentation model
    if(fSegmentation) delete fSegmentation;
    fSegmentation=thisSegmentation;
  }
  //
  void    ResponseModel(AliITSresponse* thisResponse) { 
    // Configure response model
    if(fResponse) delete fResponse;
    fResponse=thisResponse;
  }
  //
  void    SimulationModel(AliITSsimulation *thisSimulation) {
    // Configure simulation model
    fSimulation = thisSimulation;
  }
  //
  void    ReconstructionModel(AliITSClusterFinder *thisReconstruction) {
// Configure reconstruction model
      fReconst = thisReconstruction;
  }
  void    ClassNames(char * digit, char * cluster) { 
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
  
  void GetClassNames(char* &digit, char*  &cluster) 
    { 
      // Get class names for digits and rec points
      digit=fDigClassName; cluster=fClustClassName; 
    } 
  
protected:
  
  AliITSClusterFinder  *fReconst;          // cluster finder
  AliITSsimulation     *fSimulation;       // simulation
  AliITSresponse       *fResponse;         // response
  AliITSsegmentation   *fSegmentation;     // segmentation
  
  char*              fDigClassName;      //!
                                         // string
  char*              fClustClassName;    //! 
                                         // string
  
  ClassDef(AliITSDetType,1)
    
    };

#endif
