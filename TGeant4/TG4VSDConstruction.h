// $Id$
// Category: digits+hits
//
// Author: I.Hrivnacova
//
// Class TG4VSDConstruction
// ------------------------
// Abstract class for construction of user sensitive detectors.
// It has one pure virtual method Construct()
// that has to be be implemented by a derived class.
// Constructed sensitive detectors have to inherit from 
// TG4VSensitiveDetector (see TG4VSensitiveDetector.h description);
// all cloned logical volumes (which a single G3 volume correspond to)
// have to share the same sensitive detector instance.

#ifndef TG4V_SD_CONSTRUCTION_H
#define TG4V_SD_CONSTRUCTION_H

class TG4VSDConstruction
{
  public:
    TG4VSDConstruction();
    virtual ~TG4VSDConstruction();

    // methods
    virtual void Construct() = 0;
};

// inline methods

#endif //TG4V_SD_CONSTRUCTION_H

