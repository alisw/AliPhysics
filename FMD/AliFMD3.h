// -*- mode: C++ -*-
//
// $Id$
//
#ifndef ALIFMD3_H
#define ALIFMD3_H

#ifndef ALIFMDSUBDETECTOR_H
# include "AliFMDSubDetector.h"
#endif

class AliFMD3 : public AliFMDSubDetector 
{
private:
  Int_t    fVolumeId;
  Double_t fDz;
public:
  AliFMD3();
  virtual ~AliFMD3();
  virtual void   SetupGeometry(Int_t airId, Int_t kaptionId);  
  virtual void   Geometry(const char* mother, Int_t pbRotId, 
			  Int_t idRotId, Double_t z=0);
  ClassDef(AliFMD3,1); // Geometry of FMD3 
};

#endif
//
// EOF
//
