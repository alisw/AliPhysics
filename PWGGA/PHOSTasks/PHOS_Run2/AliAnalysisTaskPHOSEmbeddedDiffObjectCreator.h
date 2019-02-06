#ifndef AliAnalysisTaskPHOSEmbeddedDiffObjectCreator_cxx
#define AliAnalysisTaskPHOSEmbeddedDiffObjectCreator_cxx

// Author: Daiki Sekihata (Hiroshima University)
class TF1;
class TObjArray;
class TClonesArray;
class AliPHOSGeometry;
class AliPHOSCalibData;

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPHOSObjectCreator.h"

class AliAnalysisTaskPHOSEmbeddedDiffObjectCreator : public AliAnalysisTaskPHOSObjectCreator {
  public:
    AliAnalysisTaskPHOSEmbeddedDiffObjectCreator(const char *name = "PHOSEmbeddedObjectCreator");
    virtual ~AliAnalysisTaskPHOSEmbeddedDiffObjectCreator(); 

    //AnalysisTaskSE method is overloaded.
    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *option);
    //at least, these 3 functions above are necessary.
    void SetEmbeddedParticle(TString parname) {fParticleName = parname;}

  protected:
    void CalibrateEmbeddedClusters(TClonesArray *array, AliAODCaloCells *cellsEmb);
    Bool_t IsSameCluster(AliVCluster * c1, AliVCluster * c2) const;

  private:
    TString fParticleName;
    AliPHOSCalibData *fPHOSCalibData;

  private:
    AliAnalysisTaskPHOSEmbeddedDiffObjectCreator(const AliAnalysisTaskPHOSEmbeddedDiffObjectCreator&);
    AliAnalysisTaskPHOSEmbeddedDiffObjectCreator& operator=(const AliAnalysisTaskPHOSEmbeddedDiffObjectCreator&);

    ClassDef(AliAnalysisTaskPHOSEmbeddedDiffObjectCreator, 7);
};

#endif

