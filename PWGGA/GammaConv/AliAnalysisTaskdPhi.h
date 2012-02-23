////////////////////////////////////////////////
//--------------------------------------------- 
// Class doing conversion gamma dPhi correlations
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#ifndef AliAnalysisTaskdPhi_cxx
#define AliAnalysisTaskdPhi_cxx

#include "AliAnalysisTaskSE.h"

#include <TAxis.h>
#include <TH3I.h>
#include <THnSparse.h>
#include <AliAnalysisFilter.h>
#include <iostream>
#include <AliAnaConvCorrBase.h>
#include <AliLog.h>
#include <AliAnalysisCuts.h>
class AliAnaConvIsolation;
class AliConversionCuts;
class TList;
class TH2I;
//class THnSparseF;

using namespace std;

class AliAnalysisTaskdPhi : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskdPhi(const char *name);
  virtual ~AliAnalysisTaskdPhi();

  virtual void   UserCreateOutputObjects();
  virtual void   SetUpCorrObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  TAxis& GetAxistPt()   { return fAxistPt;   }
  TAxis& GetAxiscPt()   { return fAxiscPt;   }
  TAxis& GetAxisEta()  { return fAxisEta;  }
  TAxis& GetAxisPhi()  { return fAxisPhi;  }
  TAxis& GetAxisZ()    { return fAxisZ;    }
  TAxis& GetAxisCent() { return fAxisCent; }
  TAxis& GetAxisPiMass() { return fAxisPiM; }


  void SetV0Filter(AliConversionCuts * filter) { fV0Filter = filter; }
  AliAnalysisCuts * GetTrackCuts() const { return fTrackCuts; }
  void SetTrackCuts( AliAnalysisCuts * cuts) { if (fTrackCuts) delete fTrackCuts; fTrackCuts = cuts; }
  
  
protected:
  
  TClonesArray * GetConversionGammas(Bool_t isAOD);

private:
  
  THnSparseF * CreateSparse(TString nameString, TString titleString, TList * axesList);
  Int_t GetBin(TAxis &axis, Double_t value);
  THnSparseF * GetMEHistogram(Int_t binz, Int_t binc, TObjArray * array);
  AliAnaConvCorrBase * GetCorrObject(Int_t binz, Int_t binc, TObjArray * array);
  void Process(TObjArray * gammas, TObjArray * tracks, Int_t vertexBin, Int_t centBin);
  void FindDeltaAODBranchName(AliVEvent * event);
  
  TList * fHistograms; //histograms
  TList * fHistoGamma; //gamma histo
  TList * fHistoPion; //pion histo

  AliConversionCuts * fV0Filter; //v0 filter
  AliAnalysisCuts * fTrackCuts;

  TObjArray * fGammas; //gammas
  TObjArray * fPions; //poins

  TObjArray * hMETracks; //mixed event tracks
  TObjArray * hMEPhotons; //photons
  TObjArray * hMEPions; //pions
  TH2I * hMEvents; //event histrogam

  TObjArray * fPhotonCorr; //photon
  TObjArray * fPionCorr; //poin
  AliAnaConvIsolation * fIsoAna; //comment


  TString fDeltaAODBranchName; //comment

  TAxis fAxistPt; //comment
  TAxis fAxiscPt; //comment
  TAxis fAxisEta; //comment
  TAxis fAxisPhi; //comment
  TAxis fAxisCent; //comment
  TAxis fAxisZ; //comment
  TAxis fAxisPiM; //comment
  
  AliAnalysisTaskdPhi(const AliAnalysisTaskdPhi&); // not implemented
  AliAnalysisTaskdPhi& operator=(const AliAnalysisTaskdPhi&); // not implemented
  
  ClassDef(AliAnalysisTaskdPhi, 3); // example of analysis
};

inline THnSparseF * AliAnalysisTaskdPhi::GetMEHistogram(Int_t binz, Int_t binc, TObjArray * array) {
  ///Get Mixed Event histogram
  if(binz < 0 || binz > fAxisZ.GetNbins()) {
	cout << "error out of z axis range: " << binz << endl; 
	return NULL;
  }  
  if(binc < 0 || binc >= fAxisCent.GetNbins()) {
	cout << "error out of centraliy axis range: " << binc << endl; 
	return NULL;
  }  
  
  TObjArray * arrayc = static_cast<TObjArray*>(array->At(binz));
  THnSparseF * histogram = static_cast<THnSparseF*>(arrayc->At(binc));
  return histogram;
}


inline AliAnaConvCorrBase * AliAnalysisTaskdPhi::GetCorrObject(Int_t binz, Int_t binc, TObjArray * array) {
  ///Get correlation object
  if(binc < 0 || binz < 0) {
	  AliError("We have a bad bin!!!");
	  return NULL;
	}

  TObjArray * arrayc = static_cast<TObjArray*>(array->At(binz));
  AliAnaConvCorrBase * corrmaker = static_cast<AliAnaConvCorrBase*>(arrayc->At(binc));
  return corrmaker;

}

inline Int_t AliAnalysisTaskdPhi::GetBin(TAxis & axis, Double_t value) {
  //Return bin - 1 if within range, else return -1
  Int_t bin = axis.FindFixBin(value);
  

  bin = (bin > 0 && bin <= axis.GetNbins()) ? bin -1 : -1;
  return bin;
}

#endif

