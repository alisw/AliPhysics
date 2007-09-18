
#ifndef AliTPCSelectorTracks_h
#define AliTPCSelectorTracks_h

#include <iostream>
using namespace std;
#include <TSelector.h>

#include <TROOT.h>
#include <TFile.h>


class AliESD;
class AliESDfriend;
class TH1I;
class AliTPCcalibTracks;
class AliTPCcalibTracksGain;


class AliTPCSelectorTracks : public AliTPCSelectorESD {
public :
   AliTPCSelectorTracks(TTree *tree=0);
  virtual ~AliTPCSelectorTracks();
   virtual void    SlaveBegin(TTree *tree);
   virtual Int_t   ProcessIn(Long64_t entry);
   virtual void    Terminate();

private:
  AliTPCcalibTracks *fCalibTracks;         //! calib Tracks object
  AliTPCcalibTracksGain *fCalibTracksGain; //! gain calibration object for tracks
  static const char *fgkOutputFileName;    //! filename of the output root file
//  Int_t              fDebugLevel;        // debug level
  ClassDef(AliTPCSelectorTracks,1);
};





#endif
