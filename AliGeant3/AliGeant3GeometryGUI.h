#ifndef ALIGEANT3GEOMETRYGUI_H
#define ALIGEANT3GEOMETRYGUI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "TClonesArray.h"
#include "TGeant3.h"

class AliGuiGeomMain;
class AliDrawVolume;
class TRotMatrix;

class AliGeant3GeometryGUI : public TObject {
 public:
    AliGeant3GeometryGUI(const char* opt = "");
    virtual ~AliGeant3GeometryGUI(){}
   private:
    AliGuiGeomMain *fPanel;      // the main gui panel
    Int_t          fNstack;      // number of volumes
    TClonesArray   *fVolumes;    // array of volumes  
    Int_t          fNMaterials;  // number of materials and media
    TClonesArray   *fMaterials;  // array of materials
    TClonesArray   *fMedia;      // array of materials
    TObjArray      *fRotations;  // Rotation Matrices
 private:
    AliGeant3GeometryGUI(const AliGeant3GeometryGUI&) {}
    AliGeant3GeometryGUI & operator=(const AliGeant3GeometryGUI&) 
    {return *this;}
    
    ClassDef(AliGeant3GeometryGUI,1)  // GUI for Geant3 geometry visualisation
};



#endif
