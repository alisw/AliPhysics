#ifndef DISPLAYITSV11_H
#define DISPLAYITSV11_H
/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

//
// Class for the display of ITS // Version 11
// with the new geometry
// A proper description of this class will follow
//

#include <TTask.h>
#include <TGeoTube.h>

class TGeoManager;
class AliITS;
class TGeoVolume;
class TCanvas;

class DisplayITSv11 : public TTask {
  public:
    DisplayITSv11();
    virtual ~DisplayITSv11();
    //
    virtual void Exec(Option_t* opt="");
    virtual void DisplayITS();
    virtual void EngineeringSPDThS();
    virtual void EngineeringSDDCone();
    virtual void EngineeringSDDCylinder();
    virtual void EngineeringSupRB24();
    virtual void EngineeringSupRB26();
    //
    virtual Int_t GetDebugITS() const {return fITSdebug;}
    //
    void SetITSdebugOn(){fITSdebug=1;}
    void SetITSdebugOff(){fITSdebug=0;}
    void SetCircleSegments(Int_t i=80){fNsegments=i;}
    void SetCylindericalCutOn(Double_t phimin=0.0,Double_t phimax=180.){
        fPhimincut = phimin;fPhimaxcut=phimax;fCut=1;}
    void SetCylindericalCutOff(){fCut=0;if(fClip) delete fClip;fClip=0;}
    void SetCylindericalClipVolume(){if(fClip)delete fClip;
         fClip = new TGeoTubeSeg(0.0,fRmax[0],fRmax[2],fPhimincut,fPhimaxcut);}
    void SetPerspectiveOn(){fPerspective=1;}
    void SetPerspectiveOff(){fPerspective=0;}
    void SetSolid(){fSolid=1;}
    void SetWire(){fSolid=0;}
    void SetAxisOn(){fAxis=1;}
    void SetAxisOff(){fAxis=0;}
    void SetDisplayAngles(Double_t lon,Double_t lat,Double_t psi){
        fLongitude=lon;fLatitude=lat;fPsi=psi;}
private:
    DisplayITSv11(const DisplayITSv11 & disp) : TTask(disp) {
      Fatal("copy ctor","Not implemented\n");}
    DisplayITSv11 & operator = (const DisplayITSv11 & ) {
      Fatal("= operator","Not implemented\n"); return *this;}
    void Displaying(TGeoVolume *v,TCanvas *c,Int_t ipad);
    TGeoManager *fmgr;  // pointer to the geometry manager.
    AliITS      *fits;  // pointer to AliITS
    TGeoVolume *fALICE; // Pointer to the true mother volume
    TGeoVolume *fITS;   // Pointer to the ITS mother volume
    TGeoShape  *fClip;  // Clipping Volume.
    //
    Int_t fITSdebug;             // Debug flag
    Int_t fNsegments;            // Number of segment
    Int_t fCut;                  // Cut value
    Int_t fAxis;                 // Axis value
    Int_t fPerspective;          // Perspective option
    Int_t fSolid;                // Solid value
    Double_t fRmin[3];           // Minimum radius value 
    Double_t fRmax[3];           // Maximum radius value
    Double_t fPhimincut;         // Min phi cut
    Double_t fPhimaxcut;         // Max phi cut
    Double_t fLongitude;         // Longitute value
    Double_t fLatitude;          // Latitude value
    Double_t fPsi;               // Psi

    ClassDef(DisplayITSv11,1) // Task to display ITS v11 Geometry
};

#endif
