/// \class AliAnalysisTaskParticleFractions
/// \brief Simple task to retrieve particle fraction with respect to their origin
///
/// \author Maciej Szymanski <maszyman@cern.ch>, Warsaw University of Technology
/// \date May 15, 2015

#ifndef ALIANALYSISTASKPARTICLEFRACTIONS
#define ALIANALYSISTASKPARTICLEFRACTIONS

#include "AliAnalysisTaskSE.h"

class AliESDtrackCuts;
class TList;
class TH2F;
class TH1F;
class AliESDEvent;
class AliMCEvent;
class AliStack;
class TParticle;
class AliESDtrackCuts;
class AliCentrality;
class TObjArray;
class AliEventPoolManager;
class AliInputEventHandler;
class AliESDtrack;
class AliESDVertex;
class AliAODEvent;
class AliAODTrack;
class AliAODMCParticle;
class AliAODHandler;
class AliAODInputHandler;
class AliPIDResponse;

class AliAnalysisTaskParticleFractions :public AliAnalysisTaskSE{
public:

  AliAnalysisTaskParticleFractions() : AliAnalysisTaskSE(), centrality(0), fpidResponse(0), fHistoList(0), fHistEv(0)
  {
    for(Int_t i = 0; i < multbins*parttypes; i++)  {
      fParticleOriginMC[i] = NULL;
      fParticleOriginRec[i] = NULL;
    }
  }

  AliAnalysisTaskParticleFractions(const Char_t *partName); // default constructor
  virtual ~AliAnalysisTaskParticleFractions(); // default destructor
  virtual void UserCreateOutputObjects(); // user create output objects
  virtual void UserExec(Option_t *option); // user exec

  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi);
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP);
  bool IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP);

private:
  AliAnalysisTaskParticleFractions(const AliAnalysisTaskParticleFractions &); // copy constructor
  AliAnalysisTaskParticleFractions &operator=(const AliAnalysisTaskParticleFractions &); // operator=
  AliCentrality *centrality;
  AliPIDResponse *fpidResponse;

  static const int multbins = 1;
  static const int parttypes = 5;

  TList *fHistoList; // histo list
  TH1F *fHistEv;
  TH1F *fParticleOriginMC[multbins*parttypes];
  TH1F *fParticleOriginRec[multbins*parttypes];

  double fV1[3];


  ClassDef(AliAnalysisTaskParticleFractions, 1);

};

#endif /* ALIANALYSISTASKPARTICLEFRACTIONS */
