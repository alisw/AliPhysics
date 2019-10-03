#ifndef AliAnalysisTaskPHOSTimeCalib_cxx
#define AliAnalysisTaskPHOSTimeCalib_cxx

// class for PHOS calibratoin.
// Author: Daiki Sekihata (Hiroshima University)

class TH1F;
class TH2F;
class TH1I;
class TH2I;
class AliPHOSGeometry;
class AliVEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPHOSTimeCalib : public AliAnalysisTaskSE {
public:

	AliAnalysisTaskPHOSTimeCalib(const char *name = "AliAnalysisTaskPHOSTimeCalib");
	virtual ~AliAnalysisTaskPHOSTimeCalib(); 
	
	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);

protected:
  void CalibrateCellTime(AliVCaloCells *cells, UShort_t BC);
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key

	Int_t WhichSRU(Int_t cellx);
	Int_t WhichBranch(Int_t cellz);
  Int_t WhichDDL(Int_t module, Int_t sru);

protected:
  AliVEvent *fVEvent;
	AliPHOSGeometry *fPHOSGeo;
  THashList *fOutputContainer; //! Output list
	Int_t fRunNumber;

  AliAnalysisTaskPHOSTimeCalib(const AliAnalysisTaskPHOSTimeCalib&); // not implemented
  AliAnalysisTaskPHOSTimeCalib& operator=(const AliAnalysisTaskPHOSTimeCalib&); // not implemented
  
  ClassDef(AliAnalysisTaskPHOSTimeCalib, 0); // example of analysis
};

#endif
