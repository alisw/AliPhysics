#ifndef AliHBTAnalysisStavinskyMixing_H
#define AliHBTAnalysisStavinskyMixing_H
//_________________________________________________________
////////////////////////////////////////////////////////////////////////////
//
// class AliHBTAnalysisStavinskyMixing
//
// Analysis that does mixing with Stavinsky method i.e. 
// denominator is created out of particles originating from
// the same event, but the second particle within a pair has 
// mirrored momenta.
//
// Piotr.Skowronski@cern.ch
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//
////////////////////////////////////////////////////////////////////////////
//_________________________________________________________

#include <AliHBTAnalysis.h>

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

class AliHBTAnalysisStavinskyMixing: public AliHBTAnalysis
 {
   public:
     AliHBTAnalysisStavinskyMixing(){};
     virtual ~AliHBTAnalysisStavinskyMixing() {}

     
   protected:
     
     /**********************************************/
     /*      E V E N T     P R O C E S S I N G     */
     /**********************************************/
     // NEW AOD schema
     Int_t (AliHBTAnalysisStavinskyMixing::*fProcEvent)(AliAOD* aodrec, AliAOD* aodsim);//Pointer to the processing method

     virtual Int_t ProcessSim(AliAOD* /*aodrec*/, AliAOD* aodsim);
     virtual Int_t ProcessRec(AliAOD* aodrec, AliAOD* /*aodsim*/);
     virtual Int_t ProcessRecAndSim(AliAOD* aodrec, AliAOD* aodsim);
     
     virtual Int_t ProcessRecAndSimNonId(AliAOD* aodrec, AliAOD* aodsim);
     virtual Int_t ProcessSimNonId(AliAOD* /*aodrec*/, AliAOD* aodsim);
     virtual Int_t ProcessRecNonId(AliAOD* aodrec, AliAOD* /*aodsim*/);
     

   private:


     ClassDef(AliHBTAnalysisStavinskyMixing,0)
 };
#endif
