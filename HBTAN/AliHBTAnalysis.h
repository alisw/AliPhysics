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
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//
////////////////////////////////////////////////////////////////////////////
//_________________________________________________________

#include <AliAnalysis.h>
#include "AliAODPairCut.h"
#include "AliAODParticleCut.h"

class AliHBTCut;
//class AliHBTPair;

class AliHBTRun;
class AliAOD;
class AliReader;
class AliEventBuffer;

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
     Int_t ProcessEvent(AliAOD* aodrec, AliAOD* aodsim = 0x0);
     Int_t Finish();
     
     virtual void Process(Option_t* option = "TracksAndParticles");
     
     void SetGlobalPairCut(AliAODPairCut* cut);
     
     void AddTrackFunction(AliHBTOnePairFctn* f);
     void AddParticleFunction(AliHBTOnePairFctn* f);
     void AddParticleAndTrackFunction(AliHBTTwoPairFctn* f);
     
     void AddParticleMonitorFunction(AliHBTMonOneParticleFctn* f);    //z.ch.
     void AddTrackMonitorFunction(AliHBTMonOneParticleFctn* f);    //z.ch.
     void AddParticleAndTrackMonitorFunction(AliHBTMonTwoParticleFctn* f);//z.ch.

     void AddResolutionFunction(AliHBTTwoPairFctn* f){AddParticleAndTrackFunction(f);}
     
     void SetReader(AliReader* r){fReader = r;}
     
     void WriteFunctions();
     
     void SetBufferSize(Int_t buffsize){fBufferSize=buffsize;}
     void SetOwner(Bool_t owner=kTRUE){fIsOwner=owner;}
     Bool_t IsOwner() const {return fIsOwner;}
     Bool_t IsNonIdentAnalysis();
     void   ResetFunctions();
     void   SetDisplayInfo(Int_t howoften){fDisplayMixingInfo = howoften;}//defines every each line info about mixing is displayed
     
     void   SetCutsOnParticles(); // -- aplies only to Process Tracks And Particles
     void   SetCutsOnTracks();// -- aplies only to Process Tracks And Particles
     void   SetCutsOnTracksAndParticles();// Default // -- aplies only to Process Tracks And Particles
     
     static void PressAnyKey();//small utility function that helps to make comfortable macros
   protected:
     
     Bool_t RunCoherencyCheck();
     
     void FilterOut(AliAOD* outpart1, AliAOD* outpart2, AliAOD* inpart,
                    AliAOD* outtrack1, AliAOD* outtrack2, AliAOD* intrack)const;
     void FilterOut(AliAOD* out1, AliAOD* out2, AliAOD* in)const;
     void DeleteFunctions();
     

     /**********************************************/
     /*      E V E N T     P R O C E S S I N G     */
     /**********************************************/
     // NEW AOD schema
     Int_t (AliHBTAnalysis::*fProcEvent)(AliAOD* aodrec, AliAOD* aodsim);//Pointer to the processing method

     virtual Int_t ProcessSim(AliAOD* /*aodrec*/, AliAOD* aodsim);
     virtual Int_t ProcessRec(AliAOD* aodrec, AliAOD* /*aodsim*/);
     virtual Int_t ProcessRecAndSim(AliAOD* aodrec, AliAOD* aodsim);
     
     virtual Int_t ProcessRecAndSimNonId(AliAOD* aodrec, AliAOD* aodsim);
     virtual Int_t ProcessSimNonId(AliAOD* aodrec, AliAOD* /*aodsim*/);
     virtual Int_t ProcessRecNonId(AliAOD* /*aodrec*/, AliAOD* aodsim);
     

     // OLD legacy schema

     virtual void ProcessTracks();
     virtual void ProcessParticles();
     virtual void ProcessTracksAndParticles();
     
     virtual void ProcessTracksAndParticlesNonIdentAnal();
     virtual void ProcessParticlesNonIdentAnal();
     virtual void ProcessTracksNonIdentAnal();


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


     /**********************************************/
     /* Control parameters  */
     /**********************************************/

     AliAODPairCut*   fPairCut;//! Pair cut applied for all mixed particles
      
     Int_t  fBufferSize; //!defines the size of buffer for mixed events; -1==MIX All
     Int_t  fDisplayMixingInfo;//!defines every which particle mixing info is displayed
     Bool_t fIsOwner;//!defines of all functions are supposed to be deleted while by the way of analysis defaulr false
    
     AliEventBuffer* fPartBuffer;//Sim Particles event buffer
     AliEventBuffer* fTrackBuffer;//Rec Tracks event buffer
     
     Bool_t          fNoCorrfctns;//Internal flag indicating that no cfs are set by the user (only monitor ones)
     
   private:
   
     /**********************************************/
     /*                C U T S                     */
     /**********************************************/
   
     Bool_t (AliHBTAnalysis::*fkPass)(AliAODPair* partpair, AliAODPair* trackpair) const;//Pointer to function that performes pair cut
     Bool_t (AliHBTAnalysis::*fkPass1)(AliVAODParticle* partpair, AliVAODParticle* trackpair) const;//Pointer to function that performes cut on first particle
     Bool_t (AliHBTAnalysis::*fkPass2)(AliVAODParticle* partpair, AliVAODParticle* trackpair) const;//Pointer to function that performes cut on second particle
     Bool_t (AliHBTAnalysis::*fkPassPairProp)(AliAODPair* partpair, AliAODPair* trackpair) const;//Pointer to function that performes pair cut
     
     Bool_t PassPartAndTrack (AliAODPair* partpair, AliAODPair* trackpair) const {return (fPairCut->Pass((AliAODPair*)partpair))?kTRUE:fPairCut->Pass((AliAODPair*)trackpair);}
     Bool_t PassPartAndTrack1(AliVAODParticle* part, AliVAODParticle* track) const;
     Bool_t PassPartAndTrack2(AliVAODParticle* part, AliVAODParticle* track) const;
     Bool_t PassPairPropPartAndTrack (AliAODPair* partpair, AliAODPair* trackpair) const {return (fPairCut->PassPairProp((AliAODPair*)partpair))?kTRUE:fPairCut->PassPairProp((AliAODPair*)trackpair);}
     
     Bool_t PassPart (AliAODPair* partpair, AliAODPair* /*trackpair*/) const {return fPairCut->Pass((AliAODPair*)partpair);}
     Bool_t PassPart1(AliVAODParticle* part, AliVAODParticle* /*track*/) const {return fPairCut->GetFirstPartCut()->Pass(part);}
     Bool_t PassPart2(AliVAODParticle* part, AliVAODParticle* /*track*/) const {return fPairCut->GetSecondPartCut()->Pass(part);}
     Bool_t PassPairPropPart (AliAODPair* partpair, AliAODPair* /*trackpair*/) const {return fPairCut->PassPairProp((AliAODPair*)partpair);}
     
     Bool_t PassTrack (AliAODPair* /*partpair*/, AliAODPair* trackpair) const {return fPairCut->Pass((AliAODPair*)trackpair);}
     Bool_t PassTrack1(AliVAODParticle* /*part*/, AliVAODParticle* track) const {return fPairCut->GetFirstPartCut()->Pass(track);}
     Bool_t PassTrack2(AliVAODParticle* /*part*/, AliVAODParticle* track) const {return fPairCut->GetSecondPartCut()->Pass(track);}
     Bool_t PassPairPropTrack (AliAODPair* /*partpair*/, AliAODPair* trackpair) const {return fPairCut->PassPairProp((AliAODPair*)trackpair);}

     static const UInt_t fgkFctnArraySize;//!
     static const UInt_t fgkDefaultMixingInfo;//!
     static const Int_t  fgkDefaultBufferSize;//!

     ClassDef(AliHBTAnalysis,0)
 };
 
inline Bool_t AliHBTAnalysis::PassPartAndTrack1(AliVAODParticle* part,AliVAODParticle* track) const
{
//Checks first particle from both, particle and track pairs
  AliAODParticleCut* pc = fPairCut->GetFirstPartCut();
  return (pc->Pass(part))?kTRUE:pc->Pass(track);
}
/*************************************************************************************/ 
inline Bool_t AliHBTAnalysis::PassPartAndTrack2(AliVAODParticle* part,AliVAODParticle* track) const
{
//Checks second particle from both, particle and track pairs
  AliAODParticleCut* pc = fPairCut->GetSecondPartCut();
  return (pc->Pass(part))?kTRUE:pc->Pass(track);
}
/*************************************************************************************/ 
 
#endif
