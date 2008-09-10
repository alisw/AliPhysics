#ifndef ALITRDTRACKINGEFFICIENCYCOMBINED_H
#define ALITRDTRACKINGEFFICIENCYCOMBINED_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingEfficiencyCombined.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class AliTRDtrackingEfficiencyCombined : public AliTRDrecoTask{
public:
  AliTRDtrackingEfficiencyCombined();
  virtual ~AliTRDtrackingEfficiencyCombined(){;}
  
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *);
  virtual void Terminate(Option_t *);
    
private:
  AliTRDtrackingEfficiencyCombined(const AliTRDtrackingEfficiencyCombined &);
  AliTRDtrackingEfficiencyCombined& operator=(const AliTRDtrackingEfficiencyCombined &);
    
  ClassDef(AliTRDtrackingEfficiencyCombined, 1); // Combined tracking efficiency
};
		
#endif
