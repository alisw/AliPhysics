// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliLVTree
// ---------------------------
// Class provides methods for browsing volumes trees, 
// and setting their visualization attributes.

#ifndef ALI_LV_TREE_H
#define ALI_LV_TREE_H

#include "AliLVTreeMessenger.h"

#include <globals.hh>

class AliLVStructure;

class G4LogicalVolume;
#ifdef G4VIS_USE
class G4Colour;
#endif

class AliLVTree
{
  public:
    // --> protected
    // AliLVTree();
    // AliLVTree(const AliLVTree& right);
    virtual ~AliLVTree();

    // static methods
    static AliLVTree* Instance();

    // methods
    void List(const G4String& lvName) const;
    void List(G4LogicalVolume* lv) const;
    void ListLong(const G4String& lvName) const;
    void ListLong(G4LogicalVolume* lv) const;

#ifdef G4VIS_USE
    void SetLVTreeVisibility(G4LogicalVolume* lv, G4bool visibility) const;
    void SetVolumeVisibility(G4LogicalVolume* lv, G4bool visibility) const;
    void SetLVTreeColour(G4LogicalVolume* lv, const G4String& colName) const;
    void SetVolumeColour(G4LogicalVolume* lv, const G4String& colName) const;     
#endif

  protected:
    AliLVTree(); 
    AliLVTree(const AliLVTree& right);

    // operators
    AliLVTree& operator=(const AliLVTree &right);

  private:
    // methods
    void RegisterLogicalVolume(G4LogicalVolume* lv, const G4String& path, 
                               AliLVStructure& lvStructure) const;
    void Warn(const G4String& where, const G4String& lvName) const;			       
    void Warn(const G4String& where) const;			       

    // static data members
    static AliLVTree* fgInstance;

    // data members
    AliLVTreeMessenger  fMessenger; //messenger     
};

// inline methods

#endif //ALI_LV_TREE_H

