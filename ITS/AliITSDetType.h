#ifndef ALIITSDETTYPE_H
#define ALIITSDETTYPE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
// This Class owns the classes needed to to detector simulations and
// reconstruction. This includes the detector segmentation classes,
// the detector responce classes, the detector simulatin classes, and
// the detector reconstruction (clustering) classes for all of the ITS
// detectors.
////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TObject.h>

#include "AliITSsegmentation.h"
#include "AliITSresponse.h"

class AliITSClusterFinder;
class AliITSsimulation;

class AliITSDetType:public TObject{

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
    virtual void ReconstructionModel(AliITSClusterFinder *thisReconstruction) {
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
    void GetClassNames(char *digit,char *cluster){
	// Get class names for digits and rec points
	strcpy(digit,fDigClassName.Data());
	strcpy(cluster,fClustClassName.Data()); 
    } 
    // Return the Digit Class name
    TString GetDigitClassName(){ return fDigClassName;}
    // Return the Cluster Class name
    TString GetClusterClassName(){ return fClustClassName;}
  
 protected:
    AliITSresponse       *fResponse;         // response
    AliITSsegmentation   *fSegmentation;     // segmentation
    AliITSsimulation     *fSimulation;       // simulation
    AliITSClusterFinder  *fReconst;          // cluster finder

    TString              fDigClassName;      // string
    TString              fClustClassName;    // string

    ClassDef(AliITSDetType,1) //Detector simulation/reconstruction class holder

};

#endif
