// $Id$
// Category: interfaces
//
// Author: D. Adamova, I. Hrivnacova
//=============================================================
//
//----------------TG4GeometryGUI.h--------------------------//
//----------------AG4 Geometry Browser----------------------//
//
//===============================================================

#ifndef TG4_GEOMETRY_GUI_H
#define TG4_GEOMETRY_GUI_H

#include <TObject.h>

class TG4GUI;
class G4LogicalVolume;
class TGListTreeItem;
  

class TG4GeometryGUI : public TObject
{
public:
    TG4GeometryGUI();
    virtual ~TG4GeometryGUI(){ ;}
    
    void  ReadGeometryTree();
    void  RegisterLogicalVolume(G4LogicalVolume* lv, TGListTreeItem* itemv); 
 
 private:
    TG4GUI* fPanel;   // the main  panel
  
 private:
  TG4GeometryGUI(const TG4GeometryGUI& gg) {;}
  TG4GeometryGUI& operator=(const TG4GeometryGUI& gg) 
  {return *this;}
    
    ClassDef(TG4GeometryGUI,1)  // GUI for Geant4 geometry  
};

#endif

