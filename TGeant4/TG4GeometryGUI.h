// $Id$
// Category: interfaces
//
// Author: D. Adamova
//==============================================================
//
//----------------TG4GeometryGUI.h--------------------------//
//----------------AG4 Geometry Browser----------------------//
//
//===============================================================

#ifndef TG4_GEOMETRY_GUI_H
#define TG4_GEOMETRY_GUI_H

#include <TObject.h>

class TG4MainFrame;
class G4LogicalVolume;
class TGListTreeItem;
class G4LogicalVolumeStore;  

class TG4GeometryGUI : public TObject
{
public:
    TG4GeometryGUI();
    virtual ~TG4GeometryGUI();
    
    void  ReadGeometryTree();
    void  RegisterLogicalVolume(G4LogicalVolume* lv, TGListTreeItem* itemv);
    void  ReadMaterials() const; 

protected:

    TG4GeometryGUI(const TG4GeometryGUI& gg) ;
    TG4GeometryGUI& operator=(const TG4GeometryGUI& gg) ;

 private:
    TG4MainFrame* fPanel;   // the main  panel
  

    ClassDef(TG4GeometryGUI,1)  // GUI for Geant4 geometry  
};

#endif

