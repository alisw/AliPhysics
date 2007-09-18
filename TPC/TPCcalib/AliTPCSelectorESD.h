
#ifndef AliTPCSelectorESD_h
#define AliTPCSelectorESD_h

#include <iostream>
#include <fstream>
using namespace std;
#include <TSelector.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class AliESDEvent; 
class AliESD;
class AliESDfriend;
class TH1I;

 
class AliTPCSelectorESD : public TSelector {
public :
   AliTPCSelectorESD(TTree *tree=0);
   virtual ~AliTPCSelectorESD() { /*delete fESD; delete fESDfriend;*/ }
   virtual Int_t   Version() const { return 1; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   ReadEvent(Long64_t entry);
   virtual Int_t   ProcessIn(Long64_t entry);   

   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   void            CleanESD();
   void            DumpSysInfo(Int_t entry); // dump system info  
   
protected:
   TTree          *fChain;        //! pointer to the analyzed TTree or TChain
   TTree          *fTreeFriend;   //! pointer to friend tree
  AliESDEvent    *fESDevent;      //! esd event
   AliESD         *fESD;          //! pointer to ESD
   AliESDfriend   *fESDfriend;    //! pointer to friend
  //                USER defined variables
   Int_t           fFileNo;       //! file number
   TH1I          *fNtracks;       //! number of Tracks
   TH1I          *fNtracksFriend; //! number of firend Tracks  
   TH1I          *fNClusters;      //! number of clusters on track
   //
   // System info 
   //
   fstream        *fSysWatch;       // system watch - Memory and CPU usage 
   fstream        *fFileWatch;      // file watch   - write the status of the analyzed files
   Int_t              fDebugLevel;     //debug level

   ClassDef(AliTPCSelectorESD,1);
};





#endif
