#ifndef ALIHBTANALYSIS_H
#define ALIHBTANALYSIS_H
//_________________________________________________________
////////////////////////////////////////////////////////////////////////////
//
// class AliHBTAnalysis
//
// Central Object Of HBTAnalyser: 
// This class performs main looping within HBT Analysis
// User must plug a reader of Type AliReader
// User plugs in coorelation and monitor functions
// as well as monitor functions
//
// HBT Analysis Tool, which is integral part of AliRoot,
// ALICE Off-Line framework:
//
// Piotr.Skowronski@cern.ch
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
//
////////////////////////////////////////////////////////////////////////////
//_________________________________________________________

#include <AliAnalysis.h>

class AliHBTCut;
//class AliHBTPair;

class AliHBTRun;
class AliAOD;
class AliReader;
class AliEventBuffer;
class AliEventCut;

class AliHBTOnePairFctn;      
class AliHBTTwoPairFctn;

class AliHBTMonOneParticleFctn;
class AliHBTMonTwoParticleFctn;

class TList;

class AliHBTAnalysis: public AliAnalysis
 {
   public:
     AliHBTAnalysis();
     AliHBTAnalysis(const AliHBTAnalysis& in);
     AliHBTAnalysis& operator=(const AliHBTAnalysis& /*right*/);
     virtual ~AliHBTAnalysis();

     Int_t  Init();
     Int_t  ProcessEvent(AliAOD* aodrec, AliAOD* aodsim = 0x0);
     Int_t  Finish();
     
     enum   EProcessOption{kReconstructed,kSimulated,kSimulatedAndReconstructed};
     void   SetProcessOption(EProcessOption option){fProcessOption = option;}//Init must be called to make effect
     
     void   Process(Option_t* option = "TracksAndParticles");//Stand alone HBT analysis 
     
     void   SetGlobalPairCut(AliAODPairCut* cut);
     
     void   AddTrackFunction(AliHBTOnePairFctn* f);
     void   AddParticleFunction(AliHBTOnePairFctn* f);
     void   AddParticleAndTrackFunction(AliHBTTwoPairFctn* f);
     
     void   AddParticleMonitorFunction(AliHBTMonOneParticleFctn* f);    //z.ch.
     void   AddTrackMonitorFunction(AliHBTMonOneParticleFctn* f);    //z.ch.
     void   AddParticleAndTrackMonitorFunction(AliHBTMonTwoParticleFctn* f);//z.ch.

     void   AddResolutionFunction(AliHBTTwoPairFctn* f){AddParticleAndTrackFunction(f);}
     
     void   SetReader(AliReader* r){fReader = r;}
     
     void   WriteFunctions();
     void   SetOutputFileName(const char* fname);
          
     void   SetBufferSize(Int_t buffsize){fBufferSize=buffsize;}
     void   SetOwner(Bool_t owner=kTRUE){fIsOwner=owner;}
     Bool_t IsOwner() const {return fIsOwner;}
     Bool_t IsNonIdentAnalysis();
     void   ResetFunctions();
     void   SetDisplayInfo(Int_t howoften){fDisplayMixingInfo = howoften;}//defines every each line info about mixing is displayed
     
     void   SetApparentVertex(Double_t x, Double_t y, Double_t z);//Sets apparent vertex
     
     static void PressAnyKey();//small utility function that helps to make comfortable macros

     
   protected:
     
     /**********************************************/
     /*      E V E N T     P R O C E S S I N G     */
     /**********************************************/
     // NEW AOD schema
     Int_t (AliHBTAnalysis::*fProcEvent)(AliAOD* aodrec, AliAOD* aodsim);//Pointer to the processing method

     virtual Int_t ProcessSim(AliAOD* /*aodrec*/, AliAOD* aodsim);
     virtual Int_t ProcessRec(AliAOD* aodrec, AliAOD* /*aodsim*/);
     virtual Int_t ProcessRecAndSim(AliAOD* aodrec, AliAOD* aodsim);
     
     virtual Int_t ProcessRecAndSimNonId(AliAOD* aodrec, AliAOD* aodsim);
     virtual Int_t ProcessSimNonId(AliAOD* /*aodrec*/, AliAOD* aodsim);
     virtual Int_t ProcessRecNonId(AliAOD* aodrec, AliAOD* /*aodsim*/);
     

     // OLD legacy schema

     void   ProcessTracks();
     void   ProcessParticles();
     void   ProcessTracksAndParticles();
     
     void   ProcessTracksAndParticlesNonIdentAnal();
     void   ProcessParticlesNonIdentAnal();
     void   ProcessTracksNonIdentAnal();

     Bool_t RunCoherencyCheck();
     
     void   FilterOut(AliAOD* outpart1, AliAOD* outpart2, AliAOD* inpart,
                      AliAOD* outtrack1, AliAOD* outtrack2, AliAOD* intrack)const;
     void   FilterOut(AliAOD* out1, AliAOD* out2, AliAOD* in)const;
     void   DeleteFunctions();
     

     AliReader* fReader;//! Pointer to reader
     
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

     AliEventCut*     fBkgEventCut;// We can narrow class of events used in 

     AliEventBuffer*  fPartBuffer;//Sim Particles event buffer
     AliEventBuffer*  fTrackBuffer;//Rec Tracks event buffer

     /**********************************************/
     /* Control parameters                         */
     /**********************************************/
      
     Int_t            fBufferSize; //defines the size of buffer for mixed events; -1==MIX All
     Int_t            fDisplayMixingInfo;//!defines every which particle mixing info is displayed
     Bool_t           fIsOwner;//!defines of all functions are supposed to be deleted while by the way of analysis defaulr false
     
     EProcessOption   fProcessOption;//Option that says waht analysis to do (Rec, Sim or SimAndRec)
     Bool_t           fNoCorrfctns;//Internal flag indicating that no cfs are set by the user (only monitor ones)
     TString*         fOutputFileName;//Fiele name where to dump results, if not specified reults are written to gDirectory
     
     Double_t         fVertexX;//X position of apparent vertex
     Double_t         fVertexY;//Y position of apparent vertex
     Double_t         fVertexZ;//Z position of apparent vertex
     
   private:

     static const UInt_t fgkFctnArraySize;//!
     static const UInt_t fgkDefaultMixingInfo;//!
     static const Int_t  fgkDefaultBufferSize;//!

     ClassDef(AliHBTAnalysis,0)
 };
#endif
