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
public:
  AliFMD3();
  virtual ~AliFMD3();
  virtual void   SetupGeometry(Int_t airId, Int_t kaptionId);  
  virtual void   Geometry(const char* mother, Int_t pbRotId, 
			  Int_t idRotId, Double_t z=0);
protected:
  Int_t    fVolumeId;
  Double_t fDz;
  ClassDef(AliFMD3,1); // Geometry of FMD3 
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
//
// EOF
//
