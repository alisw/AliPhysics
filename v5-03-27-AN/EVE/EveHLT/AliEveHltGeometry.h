//-*- Mode: C++ -*-

// $Id$

#ifndef ALIEVEHLTGEOMETRY_H
#define ALIEVEHLTGEOMETRY_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice     
 */

/** @file   AliEveHltGeometry.h
    @author Svein Lindal
    @date
    @brief  Creates the HLT geometry in viewer
*/

class TEveGeoShape;
class TEveGeoTopNode;
class TGeoManager;

class AliEveHltGeometry {

public:
  
  /** default constructor */
  AliEveHltGeometry();

  /** destructor */
  virtual ~AliEveHltGeometry();

  TEveGeoShape * CreateGentleGeometry( Bool_t register_as_global = kTRUE);
  TEveGeoTopNode * CreateEmcalGeometry(TGeoManager * manager);
  TEveGeoShape* geom_gentle_rphi();
  TEveGeoShape* geom_gentle_rhoz();
  TEveGeoShape* geom_gentle_trd();
  TEveGeoShape* geom_gentle_muon(Bool_t updateScene = kTRUE);

private:

  void DrawDeep(TEveGeoShape *gsre);

  ClassDef(AliEveHltGeometry, 0); // Manage connections to HLT data-sources.
};

#endif
