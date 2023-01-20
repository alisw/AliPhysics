#ifndef AliMultInput_H
#define AliMultInput_H
#include <TNamed.h>
#include "TProfile.h"
#include <TMap.h>
#include "AliMultVariable.h"
#include "AliOADBMultSelection.h"

class AliMultInput : public TNamed {
  
public:
  AliMultInput();
  AliMultInput(const char * name, const char * title = "MultInput");
  AliMultInput(const AliMultInput& o);
  AliMultInput& operator=(const AliMultInput& o);
  ~AliMultInput();
  
  void     AddVariable ( AliMultVariable *lVar );
  AliMultVariable* GetVariable (const TString& lName) const;
  AliMultVariable* GetVariable (Long_t iIdx) const;
  Long_t GetNVariables         () const { return fNVars; }
  
  void     AddVtxZ ( TProfile *prof );
  TProfile* GetVtxZProfile   (const AliMultVariable *v) const;
  Double_t GetVtxZCorrection (const AliMultVariable *v, Double_t lVtxZ) const;
  Long_t GetNVtxZ         () const { return fNVtxZ; }
  void ClearVtxZ (); //cleanup
  
  void Clear(Option_t* option="");
  void Set(const AliMultInput* other);
  void Print(Option_t* option="") const;
  
  void SetupAutoVtxZCorrection();
  void SetupAutoVtxZCorrection(const AliOADBMultSelection *oadb); 
  
  //Aliases for multiplicity selection criteria
  enum MultSelVar {
    kAmplitude_V0A = 0,
    kAmplitude_V0A1,
    kAmplitude_V0A2,
    kAmplitude_V0A3,
    kAmplitude_V0A4,
    kAmplitude_V0C,
    kAmplitude_V0C1,
    kAmplitude_V0C2,
    kAmplitude_V0C3,
    kAmplitude_V0C4,
    kAmplitude_V0Apartial,
    kAmplitude_V0Cpartial,
    kAmplitude_V0AEq,
    kAmplitude_V0CEq,
    kAmplitude_OnlineV0A,
    kAmplitude_OnlineV0C,
    kAmplitude_V0AADC,
    kAmplitude_V0CADC,
    kFlatenicity_V0,
    kMultiplicity_ADA,
    kMultiplicity_ADC,
    knSPDClusters,
    knSPDClusters0,
    knSPDClusters1,
    knTracklets,
    knTracklets08,
    knTracklets15,
    kRefMultEta5,
    kRefMultEta8,
    kZncEnergy,
    kZpcEnergy,
    kZnaEnergy,
    kZpaEnergy,
    kZem1Energy,
    kZem2Energy,
    kZnaTower,
    kZncTower,
    kZpaTower,
    kZpcTower,
    kZnaFired,
    kZncFired,
    kZpaFired,
    kZpcFired,
    kNTracks,
    kNTracksTPCout,
    kNTracksGlobal2015,
    kNTracksGlobal2015Trigger,
    kNTracksITSsa2010,
    kNTracksINELgtONE,
    kNPartINELgtONE,
    kEvSel_VtxZ,
    kMC_NPart,
    kMC_NColl,
    kMC_NchV0A,
    kMC_NchV0C,
    kMC_NchEta05,
    kMC_NchEta08,
    kMC_NchEta10,
    kMC_NchEta14,
    kMC_b,
    kSpherocityMC,
    kSpherocityTracksMC,
    kNVariables
  };
  
  static const TString VarName[kNVariables];  //! Names of the TTree containers
  static const Bool_t VarIsInteger[kNVariables]; //! Is this an integer or not?
  
private:
  Long_t fNVars;
  TList *fVariableList; //List containing all AliMultVariables
  
  Long_t fNVtxZ;
  TList *fVariableVertexZList; //List containing variable profiles for usage
  
  TMap* fMap; //! Map variable to profile histogram if it exits
  
  ClassDef(AliMultInput, 2)
  //2 - vertex-Z correction
};
#endif
