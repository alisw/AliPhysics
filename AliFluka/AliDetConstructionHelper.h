// $Id$
//
// Author: I. Hrivnacova
//
// Class AliDetConstructionHelper
// ------------------------------
// Helper class that takes care of instantiating AliRoot modules
// and creating detector construction class.

#ifndef ALI_DET_CONSTRUCTION_HELPER_H
#define ALI_DET_CONSTRUCTION_HELPER_H

class G4VUserDetectorConstruction;
class AliFiles;

class AliDetConstructionHelper
{
  public:
    AliDetConstructionHelper();
    virtual ~AliDetConstructionHelper();
    
    // methods
    G4VUserDetectorConstruction* DetConstruction() const;

  private:
    G4VUserDetectorConstruction* fDetConstruction;//detector construction
    AliFiles*                    fFiles;          //file paths  
};

// inline methods

inline G4VUserDetectorConstruction* AliDetConstructionHelper::DetConstruction() const
{ return fDetConstruction; }

#endif //ALI_DETECTOR_CONSTRUCTION_HELPER_H

