/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSubZone.h,v 1.11 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpSubZone
/// \brief A region in zone composed of the row segments with the same 
/// motif type.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SUB_ZONE_H
#define ALI_MP_SUB_ZONE_H

#include <TObject.h>
#include <TList.h>

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
    /// Not implemented
    AliMpSubZone(const AliMpSubZone& right);
    /// Not implemented
    AliMpSubZone&  operator = (const AliMpSubZone& right);

    // data members
    AliMpVMotif*     fMotif;   ///< the motif in this subzone
    TList fSegments;///< contained row segments
    
  ClassDef(AliMpSubZone,1)  // Zone segment
};

#endif //ALI_MP_SUB_ZONE_H
