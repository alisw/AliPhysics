#ifndef ALICUTHANDLERPCM_H
#define ALICUTHANDLERPCM_H

#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <fstream>
#include <iostream>

class AliCutHandlerPCM{
  public:
    AliCutHandlerPCM();
    AliCutHandlerPCM(Int_t nMax = 10);
    virtual ~AliCutHandlerPCM() {};


    void AddCutPCM(TString eventCut, TString photonCut, TString mesonCut);
    void AddCutPCM(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut);
    void AddCutCalo(TString eventCut, TString clusterCut, TString mesonCut);
    void AddCutPCMCalo(TString eventCut, TString photonCut, TString clusterCut, TString mesonCut);
    void AddCutMergedCalo(TString eventCut, TString clusterCut, TString clusterMergedCut, TString mesonCut);

    TString GetSpecialSettingFromAddConfig (TString additionalTrainConfig, TString configString, TString fileNameMatBudWeights);
    TString GetSpecialFileNameFromString (TString fileNameExternalInputs, TString configString);

    Bool_t AreValid(){return fValidCuts;}
    Int_t GetNCuts();
    TString GetEventCut(Int_t i);
    TString GetPhotonCut(Int_t i);
    TString GetClusterCut(Int_t i);
    TString GetMesonCut(Int_t i);
    TString GetClusterMergedCut(Int_t i);
  protected:
    Int_t fMode;
    Int_t fNCuts;
    Int_t fNMaxCuts;
    Bool_t fValidCuts;
    Bool_t fValidCutsEvent;
    Bool_t fValidCutsPCM;
    Bool_t fValidCutsCalo;
    Bool_t fValidCutsMergedCalo;
    Bool_t fValidCutsMeson;
    TString* fEventCutArray;
    TString* fPhotonCutArray;
    TString* fMesonCutArray;
    TString* fClusterCutArray;
    TString* fMergedClusterCutArray;

  private:
    AliCutHandlerPCM(const AliCutHandlerPCM&);                  // Prevent copy-construction
    AliCutHandlerPCM &operator=(const AliCutHandlerPCM&);       // Prevent assignment

    ClassDef(AliCutHandlerPCM,1);
};

#endif
