#ifndef ALIHBTANALYSIS_H
#define ALIHBTANALYSIS_H
//_________________________________________________________
////////////////////////////////////////////////////////////////////////////
//
// class AliHBTAnalysis
//
// Central Object Of HBTAnalyser: 
// This class performs main looping within HBT Analysis
// User must plug a reader of Type AliHBTReader
// User plugs in coorelation and monitor functions
// as well as monitor functions
//
// HBT Analysis Tool, which is integral part of AliRoot,
// ALICE Off-Line framework:
//
// Piotr.Skowronski@cern.ch
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//
////////////////////////////////////////////////////////////////////////////
//_________________________________________________________

#include <TObject.h>

class AliHBTParticleCut;
class AliHBTCut;
class AliHBTPairCut;
class AliHBTPair;

class AliHBTRun;
class AliHBTEvent;
class AliHBTReader;
class AliHBTOnePairFctn;      
class AliHBTTwoPairFctn;

class AliHBTMonOneParticleFctn;
class AliHBTMonTwoParticleFctn;

class TList;

class AliHBTAnalysis: public TObject
 {
   public:
     AliHBTAnalysis();
     AliHBTAnalysis(const AliHBTAnalysis& in);
     const AliHBTAnalysis& operator=(const AliHBTAnalysis& /*right*/);
     virtual ~AliHBTAnalysis();

     virtual void Process(Option_t* option = "TracksAndParticles");
     
     void SetGlobalPairCut(AliHBTPairCut* cut);
     
     void AddTrackFunction(AliHBTOnePairFctn* f);
     void AddParticleFunction(AliHBTOnePairFctn* f);
     void AddParticleAndTrackFunction(AliHBTTwoPairFctn* f);
     
     void AddParticleMonitorFunction(AliHBTMonOneParticleFctn* f);    //z.ch.
     void AddTrackMonitorFunction(AliHBTMonOneParticleFctn* f);    //z.ch.
     void AddParticleAndTrackMonitorFunction(AliHBTMonTwoParticleFctn* f);//z.ch.

     void AddResolutionFunction(AliHBTTwoPairFctn* f){AddParticleAndTrackFunction(f);}
     
     void SetReader(AliHBTReader* r){fReader = r;}
     
     void WriteFunctions();
     
     void SetBufferSize(Int_t buffsize){fBufferSize=buffsize;}
     void SetOwner(Bool_t owner=kTRUE){fIsOwner=owner;}
     Bool_t IsOwner() const {return fIsOwner;}
     Bool_t IsNonIdentAnalysis();
     void   Init();
     void   ResetFunctions();
     void   SetDisplayInfo(Int_t howoften){fDisplayMixingInfo = howoften;}//defines every each line info about mixing is displayed
   protected:
     
     Bool_t RunCoherencyCheck();
     
     void FilterOut(AliHBTEvent* outpart1, AliHBTEvent* outpart2, AliHBTEvent* inpart,
                    AliHBTEvent* outtrack1, AliHBTEvent* outtrack2, AliHBTEvent* intrack);
     void FilterOut(AliHBTEvent* out1, AliHBTEvent* out2, AliHBTEvent* in);
     void DeleteFunctions();
     
     virtual void ProcessTracks();
     virtual void ProcessParticles();
     virtual void ProcessTracksAndParticles();
     
     virtual void ProcessTracksAndParticlesNonIdentAnal();
     virtual void ProcessParticlesNonIdentAnal();
     virtual void ProcessTracksNonIdentAnal();

     AliHBTReader* fReader;//! Pointer to reader
     
     UInt_t fNTrackFunctions; //! Number of Tracks functions 
     UInt_t fNParticleFunctions; //! Number of particles functions
     UInt_t fNParticleAndTrackFunctions; //! Number of resolution functions
		
     UInt_t fNTrackMonitorFunctions; //! Number of Track Monitor functions 
     UInt_t fNParticleMonitorFunctions; //! Number of Particles Monitor functions 
     UInt_t fNParticleAndTrackMonitorFunctions; //! Number of Resolution Monitor functions 

     AliHBTOnePairFctn**  fTrackFunctions; //!array of pointers to functions that analyze rekonstructed tracks
     AliHBTOnePairFctn**  fParticleFunctions; //!array of pointers to functions that analyze generated particles
     AliHBTTwoPairFctn**  fParticleAndTrackFunctions; //!array of pointers to functions that analyze both 
                                        //reconstructed tracks and generated particles
		//i.e. - resolution analyzers
     AliHBTMonOneParticleFctn**  fParticleMonitorFunctions; //! array of pointers to monitoring functions
     AliHBTMonOneParticleFctn**  fTrackMonitorFunctions; //! which are used for single particle analysis,
     AliHBTMonTwoParticleFctn**  fParticleAndTrackMonitorFunctions;  //! cut monitoring, etc.


     /**********************************************/
     /* Control parameters  */
     /**********************************************/

     AliHBTPairCut *fPairCut;//! Pair cut applied for all mixed particles
      
     Int_t  fBufferSize; //!defines the size of buffer for mixed events; -1==MIX All
     Int_t  fDisplayMixingInfo;//!defines every which particle mixing info is displayed
     Bool_t fIsOwner;//!defines of all functions are supposed to be deleted while by the way of analysis defaulr false

   private:
     static const UInt_t fgkFctnArraySize;//!
     static const UInt_t fgkDefaultMixingInfo;//!
     static const Int_t  fgkDefaultBufferSize;//!

     ClassDef(AliHBTAnalysis,0)
 };
#endif
