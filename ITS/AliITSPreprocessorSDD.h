#ifndef ALIITSPREPROCESSORSDD_H
#define ALIITSPREPROCESSORSDD_H

////////////////////////////////////////////////////
//  Class for the                                 //
//  SDD beam test digit preprocessing             //
//  Origin: E. Crescio - crescio@to.infn.it       //
//                                                //
////////////////////////////////////////////////////


#include "AliPreprocessor.h"


class AliITSPreprocessorSDD : public AliPreprocessor { 
 

 public:
 
  AliITSPreprocessorSDD(const char* detector, AliShuttleInterface* shuttle):
    AliPreprocessor(detector,shuttle){;}
  virtual ~AliITSPreprocessorSDD(){;}



 protected:      
  
  virtual UInt_t Process(TMap* dcsAliasMap);

  static const Int_t fgkNumberOfSDD;       // number of SDD modules 
  static const Int_t fgkNumberOfChannels;  // number of channels per module
  static const char* fgkNameHistoPedestals; // name of pedestal histogram
  static const char* fgkNameHistoNoise;     // name of noise histogram

  ClassDef(AliITSPreprocessorSDD,1)  // Alice ITS-SDD preprocessor.

 };



#endif

    
