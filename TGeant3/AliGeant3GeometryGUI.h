#ifndef ALIGEANT3GEOMETRYGUI_H
#define ALIGEANT3GEOMETRYGUI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "TClonesArray.h"
#include "TGeant3.h"

class AliGuiGeomMain;
class AliDrawVolume;

class AliGeant3GeometryGUI : public TObject {
 public:
    AliGeant3GeometryGUI();
    virtual ~AliGeant3GeometryGUI(){}
    
    // Reads the zebra geometry tree and put it into the ListTree
    void  ReadGeometryTree();
    // Read material and media information and put it into ComboBox 
    void  ReadMaterials();
    Float_t Cut(Int_t idmed, Int_t icut);
 private:
    AliGuiGeomMain *fPanel;      // the main gui panel
    Int_t          fNstack;      // number of volumes
    TClonesArray   *fVolumes;    // array of volumes  
    Int_t          fNMaterials;  // number of materials and media
    TClonesArray   *fMaterials;  // array of materials
    TClonesArray   *fMedia;      // array of materials    
// Zebra bank related information	
    Int_t*    fZlq;              // pointer to Zebra bank lq
    Float_t*  fZq;               // pointer to Zebra bank q
    Int_t*    fZiq;              // pointer to Zebra bank iq
    Gclink_t* fGclink;           // pointer to Geant common block 
    Gcnum_t*  fGcnum;            // pointer to Geant common block 

 private:
    virtual AliDrawVolume* Volume(Int_t id)
	{return (AliDrawVolume *) (fVolumes->UncheckedAt(id));}
    // Return number of children for volume idvol
    Int_t NChildren(Int_t idvol);
    // Return child number idc of volume idvol
    Int_t Child(Int_t idvol, Int_t idc);
    // Return medium number for given volume idvol
    Int_t Medium(Int_t idvol);
    // Return material number for given volume idvol
    Int_t Material(Int_t idvol);
    //

  AliGeant3GeometryGUI(const AliGeant3GeometryGUI&) {}
  AliGeant3GeometryGUI & operator=(const AliGeant3GeometryGUI&) 
  {return *this;}
    
    ClassDef(AliGeant3GeometryGUI,1)  // GUI for Geant3 geometry visualisation
};



#endif
