#ifndef ALIRUNANALYSIS_H
#define ALIRUNANALYSIS_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliRunAnalysis
//
// Analysis manager
//
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

class AliEventCut;
class TObjArray;
class TFile;

#include <TString.h>
#include <TTask.h>

#include "AliAnalysis.h"

class AliRunAnalysis: public TTask
{
  public: 
    AliRunAnalysis();
    virtual ~AliRunAnalysis();
    
    Int_t Run();
    void  Add(AliAnalysis* a);
    void  ReadKinematics(Bool_t flag){fReadKinematics = flag;}
    
    Int_t GetDebug() {return AliAnalysis::GetDebug();}
    void SetDirs(TObjArray* dirs){fDirs = dirs;} //sets array directories names;
  protected:
    TObjArray*    fAnalysies;//arry with analysies
    TObjArray*    fDirs;//arry with directories to read data from
    
    AliEventCut*  fEventCut;//event cut    
    
    TString       fFileName;//name of the file with ESDs
    Bool_t        fReadKinematics;
    
    TString& GetDirName(Int_t entry);
    TFile* OpenFile(Int_t n);
    
  private:
    void SetName(const char *){}//change SetName to be private
    
    ClassDef(AliRunAnalysis,1)
};

#endif
