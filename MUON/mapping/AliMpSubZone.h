// $Id$
// Category: sector
//
// Class AliMpSubZone
// ------------------
// Class describing a zone segment composed of the 
// line segments with the same motif type.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SUB_ZONE_H
#define ALI_MP_SUB_ZONE_H

#include <TObject.h>

#include "AliMpSectorTypes.h"

class AliMpVMotif;
class AliMpVRowSegment;

class AliMpSubZone : public TObject
{
  public:
    AliMpSubZone(AliMpVMotif* motif);
    AliMpSubZone();
    virtual ~AliMpSubZone();
  
    // methods
    void AddRowSegment(AliMpVRowSegment* rowSegment);
    virtual void Print(const char* /*option*/ = 0) const;

    // access methods
    Int_t              GetNofRowSegments() const;
    AliMpVRowSegment*  GetRowSegment(Int_t i) const;
    AliMpVMotif*       GetMotif() const;

  private:
    // data members
    AliMpVMotif*     fMotif;   // the motif in this subzone
    RowSegmentVector fSegments;// contained row segments
    
  ClassDef(AliMpSubZone,1)  //Zone segment
};

#endif //ALI_MP_SUB_ZONE_H
