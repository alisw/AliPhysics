//
// $Id$
//
#ifndef ALIFMD3_H
#define ALIFMD3_H

#ifndef ALIFMDSUBDETECTOR_H
# include "AliFMDSubDetector.h"
#endif
#ifndef ALIFMD3SUPPORT_H
# include "AliFMD3Support.h"
#endif

class AliFMD3 : public AliFMDSubDetector 
{
public:
  AliFMD3();
  virtual ~AliFMD3();
  virtual void   SetupGeometry(Int_t airId, Int_t alId, Int_t cId=0);  
  virtual void   Geometry(const char* mother, Int_t pbRotId, 
			  Int_t idRotId, Double_t z=0);
  virtual void   SimpleGeometry(TList* nodes, TNode* mother, 
				Int_t colour, Double_t zMother);
  virtual void   Gsatt();
protected:
  Int_t          fVolumeId;  // Volume ID
  AliFMD3Support fSupport;   // Support for FMD3 
  ClassDef(AliFMD3,2);       // Geometry of FMD3 
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
