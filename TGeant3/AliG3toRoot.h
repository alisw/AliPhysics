#ifndef ALIG3toRoot_H
#define ALIG3toRoot_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TGeant3.h>

class TClonesArray;
class TObjArray;
class TGPicture;
class TFolder;
class TGeometry;

class AliG3Volume;
class AliNode;


class AliG3toRoot : public TObject 
{
 public:
    AliG3toRoot();
    virtual ~AliG3toRoot(){}
    virtual void G3toRoot();
    void ConvertToRootShapes(TFolder *item=0, AliNode** node=0,
			     Int_t nNodes=1);
    // Setters
    virtual void SetExpandDivisions(Int_t flag = 1)
	{fExpand = flag;}
    // Getters
    TFolder*      GetTopFolder() {return fTopFolder;}
    TClonesArray* GetMaterials() {return fMaterials;}    
    TClonesArray* GetMedia()     {return fMedia;}
 private:
    void ExpandDivisions(AliNode* node=0);
    void ReadGeometryTree();
    void ReadMaterials();    
    void ReadRotations();
    TFolder* AddVolume(TObject * obj, TFolder *parent, const char* name);
    virtual AliG3Volume* Volume(Int_t id);
    Int_t Medium(Int_t idvol);
    Int_t Material(Int_t idvol);
    Float_t Cut(Int_t imed, Int_t icut);
    // Return number of children for volume idvol
    Int_t NChildren(Int_t idvol);
    // Return child number idc of volume idvol
    Int_t Child(Int_t idvol, Int_t idc);
    AliG3toRoot &operator=(const AliG3toRoot &) {return *this;}
 private:
    TClonesArray   *fVolumes;    //! array of volumes  
    TClonesArray   *fMaterials;  //! array of materials
    TClonesArray   *fMedia;      //! array of materials
    TObjArray      *fRotations;  //! Rotation Matrices
    // Zebra bank related information	
    Int_t*     fZlq;              //! pointer to Zebra bank lq
    Float_t*   fZq;               //! pointer to Zebra bank q
    Int_t*     fZiq;              //! pointer to Zebra bank iq
    Gclink_t*  fGclink;           //! pointer to Geant common block 
    Gcnum_t*   fGcnum;            //! pointer to Geant common block
    // List Tree
    TFolder*   fTopFolder;        //! Folder structure containing volumes
    TGeometry* fGeometry;         //  Pointer to geometry
    Int_t      fExpand;           //  Flag for division expansion
    
    ClassDef(AliG3toRoot,1) // Material Object for GUI 
};

#endif








