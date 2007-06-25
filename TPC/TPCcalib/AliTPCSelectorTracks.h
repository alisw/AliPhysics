
#ifndef AliTPCSelectorTracks_h
#define AliTPCSelectorTracks_h

#include <iostream>
using namespace std;
#include <TSelector.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


class AliESD;
class AliESDfriend;
class TH1I;
class AliTPCcalibTracks;


class AliTPCSelectorTracks : public TSelector {
public :
   AliTPCSelectorTracks(TTree *tree=0);
   virtual ~AliTPCSelectorTracks() { /*delete fESD; delete fESDfriend;*/ }
   virtual Int_t   Version() const { return 1; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
  void            CleanESD();
private:
   TTree          *fChain;        //! pointer to the analyzed TTree or TChain
   TTree          *fTreeFriend;   //! pointer to friend tree
   AliESD         *fESD;          //! pointer to ESD
   AliESDfriend   *fESDfriend;    //! pointer to friend
  //                USER defined variables
   Int_t           fFileNo;       //! file number
   TH1I          *fNtracks;       //! number of Tracks
   TH1I          *fNtracksFriend; //! number of firend Tracks  
   TH1I          *fNClusters;      //! number of clusters on track
  AliTPCcalibTracks *fCalibTracks; //! calib Tracks object
   ClassDef(AliTPCSelectorTracks,1);
};





#endif
