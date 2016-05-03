//-------------------------------------------------------------------------
//     OADB container for trigger analysis configuration (cut ranges.. ...)
//     Author: Michele Floris, CERN
//-------------------------------------------------------------------------

#include "AliOADBTriggerAnalysis.h"
#include "AliLog.h"
#include "TBrowser.h"
#include "TFolder.h"
#include <iostream>
#include "TObjString.h"

using namespace std;

ClassImp(AliOADBTriggerAnalysis)

AliOADBTriggerAnalysis::AliOADBTriggerAnalysis(TString name) : TNamed(name.Data(), "OADB object storing trigger analysis settings"),
fFMDLowCut(0.2),
fFMDHitCut(0.5),
fZDCCutRefSum(-568.5),
fZDCCutRefDelta(-2.1),
fZDCCutSigmaSum(3.25),
fZDCCutSigmaDelta(2.25),
fZDCCutRefSumCorr(-65.5),
fZDCCutRefDeltaCorr(-2.1),
fZDCCutSigmaSumCorr(6.0),
fZDCCutSigmaDeltaCorr(1.2),
fZDCCutZNATimeCorrMin(0.0),
fZDCCutZNATimeCorrMax(2.0),
fZDCCutZNCTimeCorrMin(0.0),
fZDCCutZNCTimeCorrMax(5.0),
fSPDClsVsTklA(65),
fSPDClsVsTklB(4),
fV0MOnVsOfA(-1200.),
fV0MOnVsOfB(8.),
fSPDOnVsOfA(-20.),
fSPDOnVsOfB(1.0),
fV0CasymA(-25),
fV0CasymB(0.15),
fNBCsPast(4),
fNBCsFuture(7),
fVHMBBAflags(-1),
fVHMBBCflags(-1),
fVHMBGAflags(-1),
fVHMBGCflags(-1),
fVIRBBAflags(10),
fVIRBBCflags(10),
fVIRBGAflags(33),
fVIRBGCflags(33),
fV0MOnThreshold(-1),
fV0MOfThreshold(400),
fSPDGFOThreshold(2),
fTRDptHSE(3.),
fTRDpidHSE(144),
fTRDptHQU(2.),
fTRDpidHQU(164.),
fTRDptHEE(3.),
fTRDpidHEE(144),
fTRDminSectorHEE(6),
fTRDmaxSectorHEE(8),
fTRDptHJT(3.),
fTRDnHJT(3)
{
}

AliOADBTriggerAnalysis::~AliOADBTriggerAnalysis(){
}
  
void AliOADBTriggerAnalysis::Browse(TBrowser *b){
   // Browse this object.
   // If b=0, there is no Browse call TObject::Browse(0) instead.
   //         This means TObject::Inspect() will be invoked indirectly


  static TObjString * strZDCCutRefSumCorr     =0;    
  static TObjString * strZDCCutRefDeltaCorr   =0;  
  static TObjString * strZDCCutSigmaSumCorr   =0;  
  static TObjString * strZDCCutSigmaDeltaCorr =0;
  static TObjString * strZDCCutZNATimeCorrMin =0;
  static TObjString * strZDCCutZNATimeCorrMax =0;
  static TObjString * strZDCCutZNCTimeCorrMin =0;
  static TObjString * strZDCCutZNCTimeCorrMax =0;

  if(strZDCCutRefSumCorr     ) delete strZDCCutRefSumCorr     ;
  if(strZDCCutRefDeltaCorr   ) delete strZDCCutRefDeltaCorr   ;
  if(strZDCCutSigmaSumCorr   ) delete strZDCCutSigmaSumCorr   ;
  if(strZDCCutSigmaDeltaCorr ) delete strZDCCutSigmaDeltaCorr ;
  if(strZDCCutZNATimeCorrMin ) delete strZDCCutZNATimeCorrMin ;
  if(strZDCCutZNATimeCorrMax ) delete strZDCCutZNATimeCorrMax ;
  if(strZDCCutZNCTimeCorrMin ) delete strZDCCutZNCTimeCorrMin ;
  if(strZDCCutZNCTimeCorrMax ) delete strZDCCutZNCTimeCorrMax ;
  
  strZDCCutRefSumCorr     = new TObjString(Form("ZDCCutRefSumCorr     %f", fZDCCutRefSumCorr    )); 
  strZDCCutRefDeltaCorr   = new TObjString(Form("ZDCCutRefDeltaCorr   %f", fZDCCutRefDeltaCorr  )); 
  strZDCCutSigmaSumCorr   = new TObjString(Form("ZDCCutSigmaSumCorr   %f", fZDCCutSigmaSumCorr  )); 
  strZDCCutSigmaDeltaCorr = new TObjString(Form("ZDCCutSigmaDeltaCorr %f", fZDCCutSigmaDeltaCorr)); 
  strZDCCutZNATimeCorrMin = new TObjString(Form("ZDCCutZNATimeCorrMin %f", fZDCCutZNATimeCorrMin));
  strZDCCutZNATimeCorrMax = new TObjString(Form("ZDCCutZNATimeCorrMax %f", fZDCCutZNATimeCorrMax));
  strZDCCutZNCTimeCorrMin = new TObjString(Form("ZDCCutZNCTimeCorrMin %f", fZDCCutZNCTimeCorrMin));
  strZDCCutZNCTimeCorrMax = new TObjString(Form("ZDCCutZNCTimeCorrMax %f", fZDCCutZNCTimeCorrMax));

  if (b) {
    // Creates a folder for each beam type containing the list of corresponding bx ids
    b->Add(strZDCCutRefSumCorr    );
    b->Add(strZDCCutRefDeltaCorr  );
    b->Add(strZDCCutSigmaSumCorr  );
    b->Add(strZDCCutSigmaDeltaCorr);
    b->Add(strZDCCutZNATimeCorrMin);
    b->Add(strZDCCutZNATimeCorrMax);
    b->Add(strZDCCutZNCTimeCorrMin);
    b->Add(strZDCCutZNCTimeCorrMax);
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
  cout << " - ZDCCutZNATimeCorrMin "<< fZDCCutZNATimeCorrMin << endl;
  cout << " - ZDCCutZNATimeCorrMax "<< fZDCCutZNATimeCorrMax << endl;
  cout << " - ZDCCutZNCTimeCorrMin "<< fZDCCutZNCTimeCorrMin << endl;
  cout << " - ZDCCutZNCTimeCorrMax "<< fZDCCutZNCTimeCorrMax << endl;
  cout << option << endl;
}
