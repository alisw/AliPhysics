// $Id$
// Category: geometry
//
// Class that associates the name tree with logical volumes tree. 
// Used for printing volumes trees.  

#ifndef ALI_LV_STRUCTURE_H
#define ALI_LV_STRUCTURE_H

#include <G4LogicalVolume.hh>
#include <globals.hh>

#include <g4rw/tpordvec.h>

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
    void AddNewVolume(G4LogicalVolume* lv, G4String treeStructure);
    void ListTree() const;
    void ListTreeLong() const;

    // set methods
    void SetVerboseLevel(G4int verbose); 
#ifdef ALICE_VISUALIZE
    void SetTreeVisibility(G4bool visibility);       
    void SetTreeColour(G4String colName);
#endif             

    // get methods
    G4LogicalVolume* GetVolume(G4String name);
    G4LogicalVolume* FindVolume(G4String name);

  protected:
    AliLVStructure(); 

  private:
    // methods
    AliLVStructure* FindSubDirectory(G4String subDir);
    G4String ExtractDirName(G4String path);

    // data members
    G4RWTPtrOrderedVector<AliLVStructure>   fStructures;                     //.
                                                //vector of
                                                //contained structures
    G4RWTPtrOrderedVector<G4LogicalVolume>  fLogicalVolumes;                 //.
                                                //vector of
                                                //contained logical volumes
						//(parallel to fStructures)
    G4String  fPathName;     //full path name
    G4String  fDirName;      //directory name
    G4int     fVerboseLevel; //verbose level
};

#endif //ALI_LV_STRUCTURE_H

