// $Id$
//
// Class AliMUONTriggerGeometryBuilder
// -----------------------------------
// MUON Trigger stations geometry construction class.
//
// Author: Philippe Crochette, LPC Clermont-Ferrand

#ifndef ALI_MUON_TRIGGER_GEOMETRY_BUILDER_H
#define ALI_MUON_TRIGGER_GEOMETRY_BUILDER_H

#include "AliMUONVGeometryBuilder.h"

class AliMUON;

class AliMUONTriggerGeometryBuilder : public AliMUONVGeometryBuilder
{
  public:
    AliMUONTriggerGeometryBuilder(AliMUON* muon);
    AliMUONTriggerGeometryBuilder(const AliMUONTriggerGeometryBuilder& rhs);
    AliMUONTriggerGeometryBuilder();
    virtual ~AliMUONTriggerGeometryBuilder();

    // operators  
    AliMUONTriggerGeometryBuilder& operator = (const AliMUONTriggerGeometryBuilder& rhs);
  
    // methods
    virtual void CreateGeometry();
    virtual void SetTransformations();
    virtual void SetSensitiveVolumes();
    
  private:
     AliMUON*  fMUON; // the MUON detector class 
        
  ClassDef(AliMUONTriggerGeometryBuilder,1) // MUON chamber geometry base class
};

#endif //ALI_MUON_TRIGGER_GEOMETRY_BUILDER_H
