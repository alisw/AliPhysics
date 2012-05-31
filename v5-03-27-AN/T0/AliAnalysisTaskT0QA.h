#ifndef AliAnalysisTaskT0QA_cxx
#define AliAnalysisTaskT0QA_cxx

// task determines mean and sigma of T0 signals  ORA, ORC, ORA-ORC, ORA+ORC/2  
// Authors: FK  

#define NPMT0 24  //number T0 of photomultipliers

class TH1F;
class TObjArray; 
class AliESDEvent;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskT0QA : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskT0QA();
  AliAnalysisTaskT0QA(const char *name);
  virtual ~AliAnalysisTaskT0QA(); 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  TObjArray* GetOffsetHistos() {return fTzeroObject;}
  
 private:
  AliESDEvent *fESD;          //! ESD object
  TObjArray   *fTzeroObject;  // array with CFDi-CFD1 and  CFDi
  TH1F        *fTzeroORA;     //! or A spectrum    
  TH1F        *fTzeroORC;     //! or C spectrum    
  TH1F        *fResolution;   //! or A minus or C spectrum    
  TH1F        *fTzeroORAplusORC; //! ORA+ORC /2 
  int         fRunNumber;
  TH2F        **fTimeVSAmplitude; //! Time vs. Amplitude
  TH2F        *fCFDVSPmtId;   //! CFDi vs pmt id
  TH2F        *fSPDVertexVST0Vertex; //! SPD vertex vs T0 vertex   
  TH2F        *fOrAvsNtracks; //! T0A vs Ntracks
  TH2F        *fOrCvsNtracks; //! T0C vs Ntracks
  TH2F        *fT0vsNtracks; //! T0A vs Ntracks
  
 
  AliAnalysisTaskT0QA(const AliAnalysisTaskT0QA&); // not implemented
  AliAnalysisTaskT0QA& operator=(const AliAnalysisTaskT0QA&); // not implemented
  
  ClassDef(AliAnalysisTaskT0QA, 1); // example of analysis
};

#endif
