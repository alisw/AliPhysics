// $Id$
// Category: interfaces
//
// Author: D. Adamova
//================================================================
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

    G4LogicalVolume* GetLogicalVolume() const;
    TGListTreeItem* GetItem() const;
    const char* GetName() const;

//---> Inlines :
    void  SetItem(TGListTreeItem* item) {fItem = item;}

//--------------------------------------------------------------------
protected:
    TG4GuiVolume(const TG4GuiVolume& gv) ;

    // operators
    TG4GuiVolume & operator=(const TG4GuiVolume& gv) ;
    
//---------------------------------------------------------------------    
    
private:      
    G4LogicalVolume*  fLogicalVolume;    // geant logical volume 
    TGListTreeItem*   fItem; // current item

    ClassDef(TG4GuiVolume,0)   
};

//
 
#endif
