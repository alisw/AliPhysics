#ifndef ALIITSDETTYPESIM_H
#define ALIITSDETTYPESIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id$ 
*/

/**********************************************************************
 * This class contains all of the "external" information needed to do *
 * detector specific simulations for the ITS.                         *
 **********************************************************************/
#include <TObject.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliITSgeom.h"

class AliITSsimulation;
class AliITSsegmentation;
class AliITSresponce;

class AliITSDetTypeSim : public TObject {
  public:
    AliITSDetTypeSim(); // Default constructor
    virtual ~AliITSDetTypeSim(); // Proper Destructor
    virtual AliITSDetTypeSim(const AliITSDetTypeSim &source);

    // AliITSgeom related functions
    // Return the pointer to the AliITSgeom object for simulation
    AliITSgeom *GetITSgeom() const{return fGeom;}
    // Set the pointer to the AliITSgeom object for simulation
    void SetITSgeom(AliITSgeom *geom){fGeom=geom;}

    // AliITSmodule related functions
    // Return Pointer to array of hits by module number
    TObjArray *GetModules(){return fHitModule;}
    // Return Pointer to hits for a specific module number
    AliITSmodule *GetModule(Int_t index){
        return (AliITSmodule*)(fHitModule->At(index));}

    // AliITSsimulation related functions
    // Set the simulation model for a give type of detector
    virtual void SetSimulationModel(Int_t dettype,AliITSsimulation *sim){
        if(fSimulation==0) fSimulation = new TObjArray();
        if(fSimulation->At(dettype)!=0) delete fSimulation->At(dettype);
        fSimulation->AddAt(sim,dettype);}
    // Return the simulation model for a give detector type
    virtual AliITSsimulation* GetSimulationModel(Int_t dettype){
        if(fSimulation==0) return 0; return fSimulation->At(dettype);}
    // Return the simulation model for a given module (checks AliITSgeom)
    virtual AliITSsimulation* GetSimulationModelByModule(Int_t module){
        if(fGeom==0) return 0;
        return GetSimulationModel(fGeom->GetModuleType(module));}

    // AliITSsegmentation related functions
    // Set the segmentation model for a give type of detector
    virtual void SetSegmentationModel(Int_t dettype,AliITSsegmentation *seg){
        if(fSegmentation==0) fSegmentation = new TObjArray();
        if(fSegmentation->At(dettype)!=0) delete fSegmentation->At(dettype);
        fSegmentation->AddAt(seg,dettype);}
    // Return the segmentation model for a give detector type
    virtual AliITSsegmentation* GetSegmentationModel(Int_t dettype){
        if(fSegmentation==0) return 0; return fSegmentation->At(dettype);}
    // Return the segmentation model for a given module (checks AliITSgeom)
    virtual AliITSsegmentation* GetSegmentationModelByModule(Int_t module){
        if(fGeom==0) return 0;
        return GetSegmentationModel(fGeom->GetModuleType(module));}

    // AliITSresponse related functions
    // Set the response model for a give module
    virtual void SetResponseModel(Int_t module,AliITSresponse *resp){
        if(fResponse==0) fResponse = new TObjArray();
        if(fResponse->At(module)!=0) delete fResponse->At(module);
        fResponse->AddAt(resp,module);}
    // Return the response model for a give detector type
    virtual AliITSresponse* GetResponseModel(Int_t module){
        if(fResponse==0) return 0; return fResponse->At(module);}

    // Sorted hit info
    virtual void InitModules(Int_t size,Int_t &nmodules);  
    virtual void FillModules(TTree *treeH, Int_t mask = 0);
    virtual void ClearModules(){if(fHitModule) fHitModule->Delete();};

    // TClonesArray of Hits related functions
    virtual void ResetHits(){fNhits=0;if(fHits!=0) fHits->Clear();}

    // TClonesArray of SDigits related functions
    virtual void ResetSDigits(){fNSDigits=0;if(fSDigits!=0) fSDigits->Clear();}

    // TClonesArray of SDigits related functions
    virtual void ResetDigits(){fNDigits=0;if(fDigits!=0) fDigits->Clear();}

  private:
    AliITSgeom   *fGeom;         //
    TObjArray    *fSimulation;   // [NDet]
    TObjArray    *fSegmentation; // [NDet]
    TObjArray    *fResponse;     // [NMod]
    TObjArray    *fPreProcess;   // [] e.g. Fill fHitModule with hits
    TObjArray    *fPostProcess;  // [] e.g. Wright Raw data
    TObjArray    *fHitModule;    //! [NMod][Nhits]
    Int_t         fNhits;        //! number of hits
    TClonesArray *fHits;         //! Local pointer to hits
    Int_t         fNSDigits;     //! number of SDigits
    TClonesArray *fSDigits;      //! [NMod][NSDigits]
    Int_t         fNDigits;      //! number of Digits
    TClonesArray *fDigits;       //! [NMod][NDigits]
    TString       fHitClassName; // String with Hit class name
    TString       fSDigClassName;// String with SDigit class name.
    TString       fDigClassName; // String with digit class name.

    ClassDef(AliITSDetTypeSim,1) // ITS Simulation structure
}
#endif
