
#ifndef AliComparisonSelector_h
#define AliComparisonSelector_h

#include <iostream>
#include <fstream>
using namespace std;
#include <TSelector.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "AliGenInfo.h"
#include "AliRecInfo.h"

class AliESDEvent; 
class AliESD;
class AliESDfriend;
class TH1I;
class AliComparisonDraw;
 
class AliComparisonSelector : public TSelector {
public :
   AliComparisonSelector(TTree *tree=0);
   virtual ~AliComparisonSelector() { /*delete fESD; delete fESDfriend;*/ }
   virtual Int_t   Version() const { return 1; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   ReadEvent(Long64_t entry);
   virtual Int_t   ProcessIn(Long64_t entry);   
  //
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   void            DumpSysInfo(Int_t entry); // dump system info  
  //
  void            Clean();
  //
  //
  //  
protected:
   TTree          *fChain;        //! pointer to the analyzed TTree or TChain
   //
   // System info 
   //
  Int_t           fFileNo;          //!file number
   fstream        *fSysWatch;       //!system watch - Memory and CPU usage 
   fstream        *fFileWatch;      //!file watch   - write the status of the analyzed files
   Int_t              fDebugLevel;  //debug level
  //
  AliMCInfo *fInfoMC;
  AliESDRecInfo *fInfoRC;
  AliComparisonDraw *fComp;

   ClassDef(AliComparisonSelector,1);
};















#endif
