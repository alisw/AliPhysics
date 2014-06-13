#ifndef ALIANALYSISTASKPARTICLEEFFICIENCY
#define ALIANALYSISTASKPARTICLEEFFICIENCY

#define MULTBINS 1
#define PARTTYPES 4

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

class AliAnalysisTaskParticleEfficiency :public AliAnalysisTaskSE{
 public:
 AliAnalysisTaskParticleEfficiency() : AliAnalysisTaskSE(),centrality(0), fHistoList(0),  fHistEv(0), fpidResponse(0)
    {

      for(Int_t i = 0; i < MULTBINS*PARTTYPES; i++)  {
    
	fGeneratedMCPrimaries[i] = NULL;
	fMCPrimariesThatAreReconstructed[i] = NULL;
	fReconstructedAfterCuts[i] = NULL;
	fReconstructedNotPrimaries[i] = NULL;
	fReconstructedPrimaries[i] = NULL;
	fContamination[i] = NULL;
      }
  
      for ( Int_t i = 0; i < 11; i++) { 
	fHistQA[i] = NULL;
	if(i<3) fHistQA2D[i] = NULL;
      }
    }

  AliAnalysisTaskParticleEfficiency(const Char_t *partName); // default constructor
  virtual ~AliAnalysisTaskParticleEfficiency(); // default destructor
  virtual void UserCreateOutputObjects(); // user create output objects
  virtual void UserExec(Option_t *option); // user exec
  //void Terminate(Option_t *option);
  
 private:
  AliAnalysisTaskParticleEfficiency(const AliAnalysisTaskParticleEfficiency &); // copy constructor
  AliAnalysisTaskParticleEfficiency &operator=(const AliAnalysisTaskParticleEfficiency &); // operator=
  //AliAODEvent *aodEvent;
  AliCentrality *centrality;
  //AliAODTrack *fTpcTracks;
  //AliAODVertex    *vertex;
  //AliAODVertex    *vtxSPD;
  //AliAODMCParticle *MCtrk;
  //AliESDtrackCuts *fTrackCuts; // ESD track cuts
  TList *fHistoList; // histo list
  //TClonesArray *arrayMC;
  TH1F *fHistEv;
  AliPIDResponse *fpidResponse;

  TH1F *fHistQA[11];
  TH2F *fHistQA2D[3];
  TH2F *fHistQAPID[5][PARTTYPES];
  TH1F* fHistEvCuts[MULTBINS];
  //TObjArray *recoParticleArray;
  TH2F *fGeneratedMCPrimaries[MULTBINS*PARTTYPES];
  TH2F *fMCPrimariesThatAreReconstructed[MULTBINS*PARTTYPES];
  TH2F *fReconstructedAfterCuts[MULTBINS*PARTTYPES];
  TH2F *fReconstructedNotPrimaries[MULTBINS*PARTTYPES];
  TH2F *fReconstructedPrimaries[MULTBINS*PARTTYPES];
  TH2F *fContamination[MULTBINS*PARTTYPES];
  TH2F *fMisidentification[MULTBINS*PARTTYPES];

  ClassDef(AliAnalysisTaskParticleEfficiency, 1);

};

#endif /* ALIANALYSISTASKPARTICLEEFFICIENCY */
