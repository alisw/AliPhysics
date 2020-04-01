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
    AliCutHandlerPCM(Int_t );
    virtual ~AliCutHandlerPCM() {};

    void AddCutPCM(TString eventCut, TString photonCut, TString mesonCut);
    void AddCutPCM(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut);
    void AddCutCalo(TString eventCut, TString clusterCut, TString mesonCut);
    void AddCutCaloCalo(TString eventCut, TString clusterCut1, TString clusterCut2, TString mesonCut);
    void AddCutPCMCalo(TString eventCut, TString photonCut, TString clusterCut, TString mesonCut);
    void AddCutMergedCalo(TString eventCut, TString clusterCut, TString clusterMergedCut, TString mesonCut);
    void AddCutPCMDalitz(TString eventCut, TString photonCut, TString mesonCut, TString electronCut);
    void AddCutHeavyMesonPCM(TString eventCut, TString photonCut, TString pionCut, TString ndmCut, TString mesonCut);
    void AddCutHeavyMesonCalo(TString eventCut, TString clusterCut, TString pionCut, TString ndmCut, TString mesonCut);
    void AddCutHeavyMesonPCMCalo(TString eventCut, TString photonCut, TString clusterCut, TString pionCut, TString ndmCut, TString mesonCut);
    void AddCutPCMMaterial(TString eventCut, TString photonCut);
    void AddCutTrackQA(TString eventCut, TString pionCut, TString kaonCut, TString protonCut, TString deuteronCut);
    void AddCutTrackQAPion(TString eventCut, TString pionCut);

    TString GetSpecialSettingFromAddConfig (TString additionalTrainConfig, TString configString, TString fileNameMatBudWeights, TString addTaskName);
    TString GetSpecialFileNameFromString (TString fileNameExternalInputs, TString configString);

    Bool_t AreValid(){return fValidCuts;}
    Int_t GetNCuts();
    TString GetEventCut(Int_t i);
    TString GetPhotonCut(Int_t i);
    TString GetClusterCut(Int_t i);
    TString GetClusterCut2(Int_t i);
    TString GetMesonCut(Int_t i);
    TString GetClusterMergedCut(Int_t i);
    TString GetElectronCut(Int_t i);
    TString GetNDMCut(Int_t i);
    TString GetPionCut(Int_t i);
    TString GetKaonCut(Int_t i);
    TString GetProtonCut(Int_t i);
    TString GetDeuteronCut(Int_t i);


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
    Bool_t fValidCutsElectron;
    Bool_t fValidCutsNDM;
    Bool_t fValidCutsChargedPion;
    Bool_t fValidCutsChargedKaon;
    Bool_t fValidCutsProton;
    Bool_t fValidCutsDeuteron;
    TString* fEventCutArray;
    TString* fPhotonCutArray;
    TString* fMesonCutArray;
    TString* fClusterCutArray;
    TString* fClusterCutArray2;
    TString* fMergedClusterCutArray;
    TString* fElectronCutArray;
    TString* fNeutralDecayMesonCutArray;
    TString* fChargedPionCutArray;
    TString* fChargedKaonCutArray;
    TString* fProtonCutArray;
    TString* fDeuteronCutArray;

  private:
    AliCutHandlerPCM(const AliCutHandlerPCM&);                  // Prevent copy-construction
    AliCutHandlerPCM &operator=(const AliCutHandlerPCM&);       // Prevent assignment

    ClassDef(AliCutHandlerPCM,5);
};

#endif
