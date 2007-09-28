/* Id: AliSelectorLYZNewMethod.h, v1.0 30/07/2007 kolk Exp */
/* derived from AliSelectorFoF.h, v1.1 01/02/2007 esimili Exp */
/* derived from AliSelector.h,v 1.10 2006/08/15 jgrosseo Exp */

// This selector is only dependent on the ESD library, if you need the whole of AliROOT use AliSelectorRL
#ifndef ALISELECTORLYZNEWMETHOD_H
#define ALISELECTORLYZNEWMETHOD_H

#include <iostream>
using namespace std;

#include <TSelector.h>
#include "AliFlowConstants.h"
#include "AliFlowLYZConstants.h"

class AliFlowEvent;
class AliFlowTrack;
class AliFlowSelection;
class AliFlowMaker;
class AliESD;
class AliESDtrack;

class TFile;
class TTree;
class TObjArray;
class TH1F;
class TProfile;

class AliSelectorLYZNewMethod : public TSelector {


  public:

    AliSelectorLYZNewMethod();
    virtual ~AliSelectorLYZNewMethod();

   // default AliSelector things
    virtual Int_t   Version() const {return 1;}
    virtual void    Begin(TTree*);
    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

  
 protected:
  
    TTree*	  fTree; 	          //! pointer to the TTree containing the events
    AliESD*       fESD;		          //! "ESD" branch in fChain
    Int_t 	  fCountFiles ;	          //! number of processed file
  
  // default AliSelector things
    TTree*  	  GetKinematics();
    void    	  CheckOptions();

 private:
  // flow things
    TFile*        fOutfile; 
    TFile*        fFirstRunFile ;      //! pointer to file from first run
    TFile*        fSecondRunFile ;     //! pointer to file from second run
    AliFlowEvent*             fFlowEvent ;        //! pointer to flow event
    AliFlowTrack*             fFlowTrack;         //! 
    TObjArray*                fFlowTracks;        //! 
    AliFlowSelection*         fFlowSelect;    	  //! pointer to flow selection
    AliFlowMaker*             fFlowMaker;    	  //!

  // enumerators 		          	  
    Int_t         fRunID;	      //! last run ID
    Int_t 	  fEventNumber ;      //! progressive enumeration of ESD events
    Int_t         fNumberOfTracks ;   //! progressive enumeration of ESD tracks
    Int_t         fNumberOfV0s ;      //! progressive enumeration of ESD V0
    Int_t         fNumberOfEvents ;   //! total number of ESD events in file
 
    //histograms
    //input
    TProfile*  h1;    //!
    TProfile*  h2;    //!
    TH1F*      h3;    //!
    TProfile*  p1;    //!
    TProfile*  p2;    //!
    TProfile*  p3;    //!
    TProfile*  p4;    //!
    TProfile*  p5;    //!
    //output
    TProfile*  fHistProFlow;         //!
    TH1F*      fHistQtheta;          //!
    TProfile*  fHistProR0thetaHar2;  //!
    TProfile*  fHistProReDtheta;     //!
    TProfile*  fHistProImDtheta;     //!

  // default AliSelector things
    void 	      DeleteKinematicsFile();
    TFile*            fKineFile;          //! pointer to Kinematics.root if the file was opened

  // to make the code checker happy
    AliSelectorLYZNewMethod(const AliSelectorLYZNewMethod&);
    AliSelectorLYZNewMethod& operator=(const AliSelectorLYZNewMethod&);

  ClassDef(AliSelectorLYZNewMethod,0);
};

#endif
