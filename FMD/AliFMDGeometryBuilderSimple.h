#ifndef ALIFMDGEOMETRYBUILDERSIMPLE_H
#define ALIFMDGEOMETRYBUILDERSIMPLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIFMDGEOMETRYBUILDER
# include <AliFMDGeometryBuilder.h>
#endif
class AliFMDRing;

//____________________________________________________________________
class AliFMDGeometryBuilderSimple : public AliFMDGeometryBuilder
{
public:
  AliFMDGeometryBuilderSimple();
  /** CTOR */
  AliFMDGeometryBuilderSimple(Bool_t detailed);
  virtual ~AliFMDGeometryBuilderSimple() {}
protected:
  /** Make a ring volume 
      @param r Ring geometry 
      @return  Ring volume */
  TGeoVolume* RingGeometry(AliFMDRing* r);
  ClassDef(AliFMDGeometryBuilderSimple,1);
};


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//

