//-------------------------------------------------------------------------
//     OADB container for trigger analysis configuration (cut ranges.. ...)
//     Author: Michele Floris, CERN
//-------------------------------------------------------------------------

#include "AliOADBTriggerAnalysis.h"
#include "AliLog.h"
#include "TBrowser.h"
#include "TFolder.h"
#include <iostream>

using namespace std;

ClassImp(AliOADBTriggerAnalysis)

AliOADBTriggerAnalysis::AliOADBTriggerAnalysis() : TNamed("AliOADBTriggerAnalysis", "OADB object storing trigger analysis settings"),   
  fZDCCutRefSumCorr(-65.5),
  fZDCCutRefDeltaCorr(-2.1),
  fZDCCutSigmaSumCorr(6.0),
  fZDCCutSigmaDeltaCorr(1.2)
{
  // default ctor

  
}
AliOADBTriggerAnalysis::AliOADBTriggerAnalysis(char* name) : TNamed(name, "OADB object storing trigger analysis settings"), 
  fZDCCutRefSumCorr(-65.5),
  fZDCCutRefDeltaCorr(-2.1),
  fZDCCutSigmaSumCorr(6.0),
  fZDCCutSigmaDeltaCorr(1.2)
{
  // ctor
  // Init();
}

// void AliOADBTriggerAnalysis::Init() {
//   // initialize pointers
  
  

// }

AliOADBTriggerAnalysis::~AliOADBTriggerAnalysis(){
  // dtor

}
  
void AliOADBTriggerAnalysis::Browse(TBrowser *b)
{
   // Browse this object.
   // If b=0, there is no Browse call TObject::Browse(0) instead.
   //         This means TObject::Inspect() will be invoked indirectly


  static TObjString * strZDCCutRefSumCorr     =0;    
  static TObjString * strZDCCutRefDeltaCorr   =0;  
  static TObjString * strZDCCutSigmaSumCorr   =0;  
  static TObjString * strZDCCutSigmaDeltaCorr =0;

  if(strZDCCutRefSumCorr     ) delete strZDCCutRefSumCorr     ;
  if(strZDCCutRefDeltaCorr   ) delete strZDCCutRefDeltaCorr   ;
  if(strZDCCutSigmaSumCorr   ) delete strZDCCutSigmaSumCorr   ;
  if(strZDCCutSigmaDeltaCorr ) delete strZDCCutSigmaDeltaCorr ;
  
  strZDCCutRefSumCorr     = new TObjString(Form("ZDCCutRefSumCorr     %f", fZDCCutRefSumCorr    )); 
  strZDCCutRefDeltaCorr   = new TObjString(Form("ZDCCutRefDeltaCorr   %f", fZDCCutRefDeltaCorr  )); 
  strZDCCutSigmaSumCorr   = new TObjString(Form("ZDCCutSigmaSumCorr   %f", fZDCCutSigmaSumCorr  )); 
  strZDCCutSigmaDeltaCorr = new TObjString(Form("ZDCCutSigmaDeltaCorr %f", fZDCCutSigmaDeltaCorr)); 



  if (b) {
    // Creates a folder for each beam type containing the list of corresponding bx ids
    b->Add(strZDCCutRefSumCorr    );
    b->Add(strZDCCutRefDeltaCorr  );
    b->Add(strZDCCutSigmaSumCorr  );
    b->Add(strZDCCutSigmaDeltaCorr);
  }     
  else
    TObject::Browse(b);
}

void AliOADBTriggerAnalysis::Print(Option_t* option) const {
  // Print Class contents
  // Option is passed to TMap::Print
  cout << "ZDC configuration" << endl;
  cout << " - ZDCCutRefSumCorr     "<< fZDCCutRefSumCorr     << endl;
  cout << " - ZDCCutRefDeltaCorr   "<< fZDCCutRefDeltaCorr   << endl;
  cout << " - ZDCCutSigmaSumCorr   "<< fZDCCutSigmaSumCorr   << endl;
  cout << " - ZDCCutSigmaDeltaCorr "<< fZDCCutSigmaDeltaCorr << endl;
  cout << option << endl;
  

}
