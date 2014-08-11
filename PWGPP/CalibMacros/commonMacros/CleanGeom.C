#include <TGeoManager.h>
#include <TChain.h>
#include "AliAnalysisTask.h"
#include "AliGeomManager.h"

// dummy task to delete geometry in Terminate
class CleanGeom : public AliAnalysisTask {
public:
  CleanGeom(const char *name, const char *title="dummy") : AliAnalysisTask(name,title)  {DefineInput(0, TChain::Class());}
  virtual void Terminate(Option_t * = "") {delete AliGeomManager::GetGeometry();}
  virtual void Exec(Option_t * = "") {}
  ClassDef(CleanGeom,0)
};
ClassImp(CleanGeom)
