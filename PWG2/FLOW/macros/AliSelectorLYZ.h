/* Id: AliSelectorLYZ.h, v1.0 30/07/2007 kolk Exp */
/* derived from AliSelectorFoF.h, v1.1 01/02/2007 esimili Exp */
/* derived from AliSelector.h,v 1.10 2006/08/15 jgrosseo Exp */

// This selector is only dependent on the ESD library, if you need the whole of AliROOT use AliSelectorRL
#ifndef ALISELECTORLYZ_H
#define ALISELECTORLYZ_H

#include <iostream>
using namespace std;

#include <TSelector.h>
#include "AliFlowConstants.h"
#include "AliFlowLYZConstants.h"

class AliFlowEvent;
class AliFlowTrack;
class AliFlowSelection;
class AliFlowMaker;
class AliFlowAnalyser; 
class AliFlowWeighter;
class AliFlowLeeYangZerosMaker;
class AliFlowLYZHist1;
class AliFlowLYZHist2;

class AliESD;
class AliESDtrack;
class AliESDv0;

class TFile;
class TTree;
class TString;
class TObjArray;

class AliSelectorLYZ : public TSelector {


  public:

    AliSelectorLYZ();
    virtual ~AliSelectorLYZ();

   // default AliSelector things
    virtual Int_t   Version() const {return 1;}
    virtual void    Begin(TTree*);
    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

  // output file name
    void            SetFlowEventFileName(TString name)          { this->fFlowEventFileName = name ; }
    void            SetFlowLYZAnalysisFileName(TString name)    { this->fFlowLYZAnalysisFileName = name ; }
    TString         GetFlowEventFileName() const		{ return this->fFlowEventFileName ; }
    TString         GetFlowLYZAnalysisFileName()const 		{ return this->fFlowLYZAnalysisFileName ; }
    
  // tasks
    void            SetDoNothing(Bool_t kt = kTRUE)      { this->fDoNothing = kt ; }
    void            SetSaveFlowEvents(Bool_t kt = kTRUE) { this->fSaveFlowEvents = kt ; }
    void            SetOnFlyAnalysis(Bool_t kt = kTRUE)  { this->fOnFlyAnalysis = kt ; }
    void            SetOnFlyWeight(Bool_t kt = kTRUE)    { this->fOnFlyWeight = kt ; }
  //lyz flags
    void            SetFirstRunLYZ(Bool_t kt)            { this->fFirstRunLYZ = kt ;  }
    Bool_t          GetFirstRunLYZ() const               { return this->fFirstRunLYZ ; }
    void            SetUseSumLYZ(Bool_t kt)              { this->fUseSumLYZ = kt ;  }
    Bool_t          GetUseSumLYZ() const                 { return this->fUseSumLYZ ; }



 protected:
  
    TTree*	  fTree; 	          //! pointer to the TTree containing the events
    AliESD*       fESD;		          //! "ESD" branch in fChain
    Int_t 	  fCountFiles ;	          //! number of processed file

  // output file names
    TString 	  fFlowEventFileName;        //! AliFlowEvents file 
    TString 	  fFlowLYZAnalysisFileName;  //! LYZ Analysis Histograms file

  // default AliSelector things
    TTree*  	  GetKinematics();
    void    	  CheckOptions();


 private:

  // flow things
    TFile*                    fFlowfile ;         //! pointer to flow event file
    TFile*                    fFirstRunFile ;     //! pointer to file from first run
    AliFlowEvent*             fFlowEvent ;        //! pointer to flow event
    AliFlowTrack*             fFlowTrack;
    TObjArray*                fFlowTracks;
    AliFlowSelection*         fFlowSelect;    	  //! pointer to flow selection
    AliFlowMaker*             fFlowMaker ;        //! flow evt maker
    AliFlowLeeYangZerosMaker* fFlowLYZ;           //! flow LYZ analysis


  // enumerators 		          	  
    Int_t             fRunID;	          //! last run ID
    Int_t 	      fEventNumber ;      //! progressive enumeration of ESD events
    Int_t             fNumberOfTracks ;   //! progressive enumeration of ESD tracks
    Int_t             fNumberOfV0s ;      //! progressive enumeration of ESD V0
    Int_t             fNumberOfEvents ;   //! total number of ESD events in file
 
  // Flags
    Bool_t            fDoNothing ;        //! flag for a dummy execution 
    Bool_t            fSaveFlowEvents ;   //! flag for saving the flow events 
    Bool_t            fOnFlyAnalysis ;    //! flag for on-fly analysis 
    Bool_t            fOnFlyWeight ;      //! flag for on-fly analysis 
    Bool_t            fFirstRunLYZ ;      //! flag for lyz analysis 
    Bool_t            fUseSumLYZ ;        //! flag for lyz analysis 


  // default AliSelector things
    void 	      DeleteKinematicsFile();
    TFile*            fKineFile;          //! pointer to Kinematics.root if the file was opened

  // to make the code checker happy
    AliSelectorLYZ(const AliSelectorLYZ&);
    AliSelectorLYZ& operator=(const AliSelectorLYZ&);

  ClassDef(AliSelectorLYZ,0);
};

#endif
