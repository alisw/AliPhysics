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
class AliAnaConvIsolation;
//class AliConversionPi0Filter;
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

  AliAnalysisFilter& GetDielV0Filter()      { return fDielV0Filter;      }
  AliAnalysisFilter& GetDielV0TrackFilter() { return fDielV0TrackFilter; }
  AliAnalysisFilter& GetDielTrackFilter()   { return fDielTrackFilter;   }
  AliAnalysisFilter& GetDielPi0Filter()     { return fDielPi0Filter;     }

  TAxis& GetAxistPt()   { return fAxistPt;   }
  TAxis& GetAxiscPt()   { return fAxiscPt;   }
  TAxis& GetAxisEta()  { return fAxisEta;  }
  TAxis& GetAxisPhi()  { return fAxisPhi;  }
  TAxis& GetAxisZ()    { return fAxisZ;    }
  TAxis& GetAxisCent() { return fAxisCent; }
  TAxis& GetAxisPiMass() { return fAxisPiM; }

  // void SetDielV0Filter(AliAnalysisFilter * filter) { fVDielV0Filter = filter; }
  // void SetDielPi0Filter(AliAnalysisFilter * filter) { fDielPi0Filter = filter; }
  // void SetDielV0TrackFilter(AliAnalysisFilter * filter) { fVDielV0TrackFilter = filter; }
  // void SetDielTrackFilter(AliAnalysisFilter * filter) { fTDielrackFilter = filter; }

  void SetV0Filter(AliConversionCuts * filter) { fV0Filter = filter; }
  //void SetPi0Filter(AliConversionPi0Filter * filter) { fPionFilter = filter; }
  
  
  //enum kAxes { kVertexZ, kCentrality, kEta, kPhi, kPt };
  
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

  AliAnalysisFilter  fDielV0TrackFilter; //Track filter
  AliAnalysisFilter  fDielV0Filter; //v0 filter
  AliAnalysisFilter  fDielPi0Filter; //pion filter
  AliAnalysisFilter  fDielTrackFilter; //track filter

  AliConversionCuts * fV0Filter; //v0 filter
  //AliConversionPi0Filter * fPionFilter;

  TObjArray * fGammas; //gammas
  TObjArray * fPions; //poins

  TObjArray * hMETracks; //mixed event tracks
  TObjArray * hMEPhotons; //photons
  TObjArray * hMEPions; //pions
  TH2I * hMEvents; //event histrogam

  TObjArray * fPhotonCorr; //photon
  TObjArray * fPionCorr; //poin
  AliAnaConvIsolation * fIsoAna; //comment

  Int_t fL1; //comment
  Int_t fL2; //comment

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
  
  ClassDef(AliAnalysisTaskdPhi, 2); // example of analysis
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

