// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliLVStructure
// --------------------
// Class that associates the name tree with logical volumes tree. 
// Used for printing volumes trees.  

#ifndef ALI_LV_STRUCTURE_H
#define ALI_LV_STRUCTURE_H

#include <globals.hh>
//#include <g4rw/tpordvec.h>
#include <g4std/vector>

class G4LogicalVolume;
class AliLVStructure;

typedef G4std::vector<AliLVStructure*>  LVStructuresVector;
typedef G4std::vector<G4LogicalVolume*> LogicalVolumesVector;

class AliLVStructure 
{
  public:
    AliLVStructure(G4String aPath);
    AliLVStructure(const AliLVStructure& right);
    // --> protected 
    // AliLVStructure();
    virtual ~AliLVStructure();

    // operators
    AliLVStructure& operator=(const AliLVStructure& right);
    G4int operator==(const AliLVStructure &right) const;

    // methods
    void AddNewVolume(G4LogicalVolume* lv, const G4String& treeStructure);
    void ListTree() const;
    void ListTreeLong() const;

    // set methods
    void SetVerboseLevel(G4int verbose); 
#ifdef G4VIS_USE
    void SetTreeVisibility(G4bool visibility);       
    void SetTreeColour(const G4String& colName);
#endif             

    // get methods
    G4LogicalVolume* GetVolume(const G4String& name) const;
    G4LogicalVolume* FindVolume(const G4String& name) const;

  protected:
    AliLVStructure(); 

  private:
    // methods
    AliLVStructure* FindSubDirectory(const G4String& subDir) const;
    G4String ExtractDirName(const G4String& path) const;
    void ClearAndDestroy(LVStructuresVector* structures); 

    // data members
    LVStructuresVector    fStructures;     //vector of contained structures
    LogicalVolumesVector  fLogicalVolumes; //vector of contained logical volumes
				           //(parallel to fStructures)
    G4String  fPathName;     //full path name
    G4String  fDirName;      //directory name
    G4int     fVerboseLevel; //verbose level
};

#endif //ALI_LV_STRUCTURE_H

