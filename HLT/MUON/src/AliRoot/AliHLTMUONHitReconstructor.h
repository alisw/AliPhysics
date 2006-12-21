///////////////////////////////////////////////
//Author : Indranil Das, SINP, INDIA
//         Sukalyan Chattopadhyay, SINP, INDIA
//         
//
//Email :  indra.das@saha.ac.in
//         sukalyan.chattopadhyay@saha.ac.in 
///////////////////////////////////////////////

#ifndef ALIHLTMUONHITRECONSTRUCTOR_H
#define ALIHLTMUONHITRECONSTRUCTOR_H

#include <TObject.h>
#include <TString.h>
#include <TStopwatch.h>
#include "../PubSub/HLTMUONHitReconstructor.h"

class TTree;
class TFile;
class TClonesArray;

class AliRawReader;


class AliHLTMUONHitReconstructor : public TObject
{

 public:
  AliHLTMUONHitReconstructor(AliRawReader* rawReader);
  virtual ~AliHLTMUONHitReconstructor(void);

  //Needed to initialise this
  Bool_t Init(const char* lutpath, const char* buspatchmappath);

  void SetDCCut(Int_t dcCut) {fkHLTRec.SetDCCut(dcCut);}
  Bool_t WriteDHLTRecHits(Int_t iEvent);

 protected:
   AliHLTMUONHitReconstructor(const AliHLTMUONHitReconstructor& rhs); // copy constructor
   AliHLTMUONHitReconstructor& operator=(const AliHLTMUONHitReconstructor& rhs); // assignment operator
 private:

   HLTMUONHitReconstructor fkHLTRec;               // class used in online hiteconstruction
   AliRawReader* fRawReader;                     // pointer to read rawdata file

   TFile* f1;                                    // file pointer to output root file
   TClonesArray* fRecPoints;                     // TClonesArray of type AliMUONRawCluster
   TString fLutPath;                             // lookupTablePath
   DHLTLut* fLut[8];                             // pointers for lookuptables
   TTree* fDHLTTree;                             // Tree Pointer to dHLT Tree
   TStopwatch fRecTimer;                         // Stopwatch for recpoints loading 

   //private functions to read lookuptables and buspatch mapping. 
   Bool_t ReadLookUpTable(const char* lutpath);
   Bool_t ReadBusPatchToDetElemFile(BusToDetElem& busToDetElem, const char* buspatchmappath);

   ClassDef(AliHLTMUONHitReconstructor, 1)  // dHLT tracker algorithm
};

#endif // ALIHLTMUONHITRECONSTRUCTORV1_H
