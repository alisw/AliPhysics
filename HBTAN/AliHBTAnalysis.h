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
#include "AliHBTPairCut.h"
#include "AliHBTParticleCut.h"

class AliHBTCut;
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
     
     void   SetCutsOnParticles(); // -- aplies only to Process Tracks And Particles
     void   SetCutsOnTracks();// -- aplies only to Process Tracks And Particles
     void   SetCutsOnTracksAndParticles();// Default // -- aplies only to Process Tracks And Particles
     
     static void PressAnyKey();//small utility function that helps to make comfortable macros
   protected:
     
     Bool_t RunCoherencyCheck();
     
     void FilterOut(AliHBTEvent* outpart1, AliHBTEvent* outpart2, AliHBTEvent* inpart,
                    AliHBTEvent* outtrack1, AliHBTEvent* outtrack2, AliHBTEvent* intrack)const;
     void FilterOut(AliHBTEvent* out1, AliHBTEvent* out2, AliHBTEvent* in)const;
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

     AliHBTPairCut*          fPairCut;//! Pair cut applied for all mixed particles
      
     Int_t  fBufferSize; //!defines the size of buffer for mixed events; -1==MIX All
     Int_t  fDisplayMixingInfo;//!defines every which particle mixing info is displayed
     Bool_t fIsOwner;//!defines of all functions are supposed to be deleted while by the way of analysis defaulr false

   private:
     Bool_t (AliHBTAnalysis::*fkPass)(AliHBTPair* partpair, AliHBTPair* trackpair) const;//Pointer to function that performes pair cut
     Bool_t (AliHBTAnalysis::*fkPass1)(AliHBTParticle* partpair, AliHBTParticle* trackpair) const;//Pointer to function that performes cut on first particle
     Bool_t (AliHBTAnalysis::*fkPass2)(AliHBTParticle* partpair, AliHBTParticle* trackpair) const;//Pointer to function that performes cut on second particle
     Bool_t (AliHBTAnalysis::*fkPassPairProp)(AliHBTPair* partpair, AliHBTPair* trackpair) const;//Pointer to function that performes pair cut
     
     Bool_t PassPartAndTrack (AliHBTPair* partpair, AliHBTPair* trackpair) const {return (fPairCut->Pass(partpair))?kTRUE:fPairCut->Pass(trackpair);}
     Bool_t PassPartAndTrack1(AliHBTParticle* part, AliHBTParticle* track) const;
     Bool_t PassPartAndTrack2(AliHBTParticle* part, AliHBTParticle* track) const;
     Bool_t PassPairPropPartAndTrack (AliHBTPair* partpair, AliHBTPair* trackpair) const {return (fPairCut->PassPairProp(partpair))?kTRUE:fPairCut->PassPairProp(trackpair);}
     
     Bool_t PassPart (AliHBTPair* partpair, AliHBTPair* /*trackpair*/) const{return fPairCut->Pass(partpair);}
     Bool_t PassPart1(AliHBTParticle* part, AliHBTParticle* /*track*/) const{return fPairCut->GetFirstPartCut()->Pass(part);}
     Bool_t PassPart2(AliHBTParticle* part, AliHBTParticle* /*track*/) const{return fPairCut->GetSecondPartCut()->Pass(part);}
     Bool_t PassPairPropPart (AliHBTPair* partpair, AliHBTPair* /*trackpair*/) const{return fPairCut->PassPairProp(partpair);}
     
     Bool_t PassTrack (AliHBTPair* /*partpair*/, AliHBTPair* trackpair) const{return fPairCut->Pass(trackpair);}
     Bool_t PassTrack1(AliHBTParticle* /*part*/, AliHBTParticle* track) const{return fPairCut->GetFirstPartCut()->Pass(track);}
     Bool_t PassTrack2(AliHBTParticle* /*part*/, AliHBTParticle* track) const{return fPairCut->GetSecondPartCut()->Pass(track);}
     Bool_t PassPairPropTrack (AliHBTPair* /*partpair*/, AliHBTPair* trackpair) const{return fPairCut->PassPairProp(trackpair);}

     static const UInt_t fgkFctnArraySize;//!
     static const UInt_t fgkDefaultMixingInfo;//!
     static const Int_t  fgkDefaultBufferSize;//!

     ClassDef(AliHBTAnalysis,0)
 };
 
inline Bool_t AliHBTAnalysis::PassPartAndTrack1(AliHBTParticle* part,AliHBTParticle* track) const
{
//Checks first particle from both, particle and track pairs
  AliHBTParticleCut* pc = fPairCut->GetFirstPartCut();
  return (pc->Pass(part))?kTRUE:pc->Pass(track);
}
/*************************************************************************************/ 
inline Bool_t AliHBTAnalysis::PassPartAndTrack2(AliHBTParticle* part,AliHBTParticle* track) const
{
//Checks second particle from both, particle and track pairs
  AliHBTParticleCut* pc = fPairCut->GetSecondPartCut();
  return (pc->Pass(part))?kTRUE:pc->Pass(track);
}
/*************************************************************************************/ 
 
#endif
