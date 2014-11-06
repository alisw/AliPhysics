#ifndef AliTPCCALIBV0_H
#define AliTPCCALIBV0_H


#include <AliTPCcalibBase.h>


class TTreeSRedirector;
class AliTPCROC;
class AliTPCseed;
//class AliESDtrack;
//class AliESDEvent;
class AliVEvent;
class AliVTrack;
class TH3F;
class TH1F;
class TH2F;
class TH1I;
class TDatabasePDG;
class AliKFParticle;
class AliKFVertex;
class AliESDv0;
class TArrayI;
class TTree;
class AliStack;
class  AliExternalTrackParam;
class  AliTPCCorrection;

class AliTPCcalibV0 : public AliTPCcalibBase {
public :

   // List of branches

  AliTPCcalibV0();
  AliTPCcalibV0(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibV0();
  virtual void     Process(AliVEvent *event) {return ProcessESD(event);}
  void FilterV0s(AliVEvent* event);
  Long64_t Merge(TCollection *const li);
  void AddTree(TTree * treeInput);
  void AddTreeHPT(TTree * treeInput);
  static AliExternalTrackParam * RefitTrack(AliTPCseed *seed, AliTPCCorrection * corr, Double_t xstart, Double_t xstop, Double_t mass);
  static void MakeFitTreeTrack(const TObjArray * corrArray, Double_t ptCut, Int_t run);
  static void MakeFitTreeV0(const TObjArray * corrArray, Double_t ptCut, Int_t run);
  //
  //
  //
  void ProcessESD(AliVEvent *event);
  void DumpToTree(AliVEvent *event);
  void DumpToTreeHPT(AliVEvent *event);
  TTree * GetV0Tree(){return fV0Tree;}
  TTree * GetHPTTree(){return fHPTTree;}
  //  
  //  void MakeMC();
  //  void MakeV0s();
  //void ProcessV0(Int_t ftype);
  //void ProcessPI0();
  void BinLogX(TH2F * h);
  //
  //
  //  
  static AliKFParticle * Fit(AliKFVertex & primVtx, AliESDv0 *v0, Int_t PDG1, Int_t PDG2);

protected:
private:

  AliTPCcalibV0(const AliTPCcalibV0&); // Not implemented
  AliTPCcalibV0& operator=(const AliTPCcalibV0&); // Not implemented
  TTree          *fV0Tree;      // tree with full V0 information
  TTree          *fHPTTree;      // tree with high mometa tracks - full calib info
  //
  AliStack       *fStack;        // pointer to kinematic tree        
  AliVEvent      *fEvent;              //! current ED to proccess - NOT OWNER
  TDatabasePDG   *fPdg;              //! particle database
  TObjArray      *fParticles;         // array of selected MC particles
  TObjArray      *fV0s;               // array of V0s
  TObjArray      *fGammas;           // gamma conversion candidates
  //
  //void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliVTrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}  
  //       
  ClassDef(AliTPCcalibV0,3);
};


#endif
