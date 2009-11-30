#ifndef ALIPROTONSPECTRACORRECTION_H
#define ALIPROTONSPECTRACORRECTION_H

//-------------------------------------------------------------------------
//               Class AliProtonSpectraCorrection
//   This is the class for the absorption corrections used for 
//   the baryon (proton) ratio analysis
//
//    Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"
#include "TH1I.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

class TF1;
class TH2D;
class TH1F;
class TList;

class AliPID;
class AliCFDataGrid;
class AliAODEvent;
class AliAODtrack;
class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;
class AliStack;
class AliESDVertex;
class AliProtonAnalysisBase;
class AliMCEvent;

class AliProtonSpectraCorrection : public TObject {
 public:
  enum {
    kStepGenerated       = 0,
    kStepReconstructible = 1,
    kStepReconstructed   = 2,
    kStepIdentified      = 3,
    kStepSelected        = 4,
    kNSteps = 5
  };

  AliProtonSpectraCorrection();
  virtual ~AliProtonSpectraCorrection();
  
  void SetBaseAnalysis(AliProtonAnalysisBase * const baseAnalysis) {
    fProtonAnalysisBase = baseAnalysis;}
  AliProtonAnalysisBase *GetProtonAnalysisBaseObject() const {
    return fProtonAnalysisBase;}
		
  void InitAnalysisHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY, 
			      Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  void FillCorrectionMaps(AliESDEvent *fESD, 
			  const AliESDVertex *vertex,
			  AliMCEvent *mcEvent);
  void FillCorrectionMaps(AliAODEvent *fAOD);
		
  AliCFContainer *GetProtonContainer() const {
    return fCFManagerProtons->GetParticleContainer();}
  AliCFContainer *GetAntiProtonContainer() const {
    return fCFManagerAntiProtons->GetParticleContainer();}
  
 private:
  AliProtonSpectraCorrection(const AliProtonSpectraCorrection&); // Not implemented
  AliProtonSpectraCorrection& operator=(const AliProtonSpectraCorrection&); // Not implemented
  
  AliProtonAnalysisBase *fProtonAnalysisBase;//base analysis object
  
  Int_t fNBinsY; //number of bins in y or eta
  Double_t fMinY, fMaxY; //min & max value of y or eta
  Int_t fNBinsPt;  //number of bins in pT
  Double_t fMinPt, fMaxPt; //min & max value of pT
  
  //Analysis containers
  AliCFManager   *fCFManagerProtons;      // CF manager protons
  AliCFManager   *fCFManagerAntiProtons;  // CF manager antiprotons
    
  ClassDef(AliProtonSpectraCorrection,1);
};

#endif

