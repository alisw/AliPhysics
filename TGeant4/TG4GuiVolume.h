// $Id$
// Category: interfaces
//
// Author: D. Adamova
//===============================================================
// 
//----------------TG4GuiVolume.h-------------------------------//
//----Creating link for Logical Volume Tree in AG4 Geometry----//
//
//=================================================================

#ifndef TG4_GUI_VOLUME_H
#define TG4_GUI_VOLUME_H
 
#include <TObject.h>
#include <TGListTree.h>
 
class G4LogicalVolume;
 
class TG4GuiVolume : public TObject 
{
public:
    TG4GuiVolume(const char* name, G4LogicalVolume* lvolume);
    virtual ~TG4GuiVolume(){;}

//---> Inlines :

    virtual const char*   GetName() const { return fkName;}
    virtual G4LogicalVolume* GetLogicalVolume() const;
    virtual void  SetItem(TGListTreeItem* item) {fItem = item;}
    virtual TGListTreeItem* GetItem() {return fItem;}
	                                              
private:
    const char*       fkName;    //name of the gui volume                   
    G4LogicalVolume*  fLogicalVolume;    // geant logical volume 
    TGListTreeItem*   fItem; // current item

  TG4GuiVolume(const TG4GuiVolume& gv) {}
  TG4GuiVolume & operator=(const TG4GuiVolume& gv) {return *this;}

    ClassDef(TG4GuiVolume,0)   
};

//
 
#endif
