#ifndef ALIHBTANALYSIS_H
#define ALIHBTANALYSIS_H

#include <TObject.h>



class AliHBTParticleCut;
class AliHBTCut;
class AliHBTPairCut;
class AliHBTPair;

class AliHBTRun;
class AliHBTReader;
class AliHBTTwoPartFctn;      
class AliHBTFourPartFctn;


class TList;

class AliHBTAnalysis: public TObject
 {
   public:
     AliHBTAnalysis();

     virtual ~AliHBTAnalysis();

     virtual void Process(Option_t* option = "TracksAndParticles");
     

     void SetGlobalPairCut(AliHBTPairCut* cut);
     
     void AddTrackFunction(AliHBTTwoPartFctn*);
     void AddParticleFunction(AliHBTTwoPartFctn*);
     void AddParticleAndTrackFunction(AliHBTFourPartFctn*);
     
     void AddResolutionFunction(AliHBTFourPartFctn* f){AddParticleAndTrackFunction(f);}
     
     void SetReader(AliHBTReader* r){fReader = r;}
     
     void Write();
   protected:
     
     Bool_t RunCoherencyCheck();
     
     
     AliHBTReader* fReader;
     
     virtual void ProcessTracks();
     virtual void ProcessParticles();
     virtual void ProcessTracksAndParticles();
     
     
     AliHBTTwoPartFctn**  fTrackFunctions; //array of pointers to functions that analyze rekonstructed tracks
     AliHBTTwoPartFctn**  fParticleFunctions; //array of pointers to functions that analyze generated particles
     AliHBTFourPartFctn** fParticleAndTrackFunctions; //array of pointers to functions that analyze both 
                                        //reconstructed tracks and generated particles
		//i.e. - resolution analyzers
     UInt_t fNTrackFunctions;
     UInt_t fNParticleFunctions;
     UInt_t fNParticleAndTrackFunctions;
		
     /**********************************************/
     /* Control parameters  */

      AliHBTPairCut *fPairCut;
      
   // AliHBTCut *fParticleCut; 
     /**********************************************/
     
     
   private:
     static const Int_t fgkHbtAnalyzeAll;
     static const UInt_t fgkFctnArraySize;
/*********************************************/   
   public:
     ClassDef(AliHBTAnalysis,0)
 };




#endif
