/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSubZone.h,v 1.8 2005/09/26 16:12:11 ivana Exp $

/// \ingroup sector
/// \class AliMpSubZone
/// \brief A region in zone composed of the row segments with the same 
/// motif type.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SUB_ZONE_H
#define ALI_MP_SUB_ZONE_H

#include "AliMpContainers.h"

#ifdef WITH_STL
#include <vector>
#endif

#ifdef WITH_ROOT
#include <TList.h>
#endif

#include <TObject.h>

class AliMpVMotif;
class AliMpVRowSegment;

class AliMpSubZone : public TObject
{
  public:
#ifdef WITH_STL
    typedef std::vector<AliMpVRowSegment*>  RowSegmentVector;
#endif
#ifdef WITH_ROOT
    typedef TList  RowSegmentVector;
#endif

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

  protected:
    AliMpSubZone(const AliMpSubZone& right);
    AliMpSubZone&  operator = (const AliMpSubZone& right);

  private:
    // data members
    AliMpVMotif*     fMotif;   // the motif in this subzone
    RowSegmentVector fSegments;// contained row segments
    
  ClassDef(AliMpSubZone,1)  //Zone segment
};

#endif //ALI_MP_SUB_ZONE_H
