/* Id: AliSelectorFoF.h, v1.1 01/02/2007 esimili Exp */
/* derived from AliSelector.h,v 1.10 2006/08/15 jgrosseo Exp */

// This selector is only dependent on the ESD library, if you need the whole of AliROOT use AliSelectorRL
#ifndef ALISELECTORF_H
#define ALISELECTORF_H

#include <TSelector.h>

class AliFlowEvent;
class AliFlowSelection;
class AliFlowMaker;
class AliFlowAnalyser; 
class AliFlowWeighter;

class AliESD;
class AliESDtrack;
class AliESDv0;

class TFile;
class TTree;
class TString;

class AliSelectorFoF : public TSelector {


  public:

    AliSelectorFoF();
    virtual ~AliSelectorFoF();

   // defalt AliSelector things
    virtual Int_t   Version() const {return 1;}
    virtual void    Begin(TTree*);
    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

  // output file name
    void            SetFlowEventFileName(TString name)    { fFlowEventFileName = name ; }
    void            SetFlowWeightFileName(TString name)   { fFlowWgtFileName = name ; }
    void            SetFlowAnalysisFileName(TString name) { fFlowAnalysisFileName = name ; }
    TString         GetFlowEventFileName() 		  { return fFlowEventFileName ; }
    TString         GetFlowWeightFileName()		  { return fFlowWgtFileName ; } 
    TString         GetFlowAnalysisFileName()		  { return fFlowAnalysisFileName ; } 

  // tasks
    void            SetDoNothing(Bool_t kt = kTRUE)      { fDoNothing = kt ; }
    void            SetSaveFlowEvents(Bool_t kt = kTRUE) { fSaveFlowEvents = kt ; }
    void            SetOnFlyAnalysis(Bool_t kt = kTRUE)  { fOnFlyAnalysis = kt ; }
    void            SetOnFlyWeight(Bool_t kt = kTRUE)    { fOnFlyWeight = kt ; }


 protected:
  
    TTree	  *fTree; 	          //! pointer to the TTree containing the events
    AliESD*       fESD;		          //! "ESD" branch in fChain
    Int_t 	  fCountFiles ;	          //! number of processed file

  // output file names
    TString 	  fFlowEventFileName;     //! AliFlowEvents file 
    TString 	  fFlowAnalysisFileName;  //! Analysis Histograms file
    TString 	  fFlowWgtFileName;       //! Weights file

  // defalt AliSelector things
    TTree*  	  GetKinematics();
    void    	  CheckOptions();


 private:

  // flow things
    TFile*            fFlowfile ;         //! pointer to flow event file
    TFile*            fFlowWgtFile ;      //! pointer to flow weight file
    AliFlowEvent*     fFlowEvent ;        //! pointer to flow event
    AliFlowSelection* fFlowSelect;    	  //! pointer to flow selection
    AliFlowMaker*     fFlowMaker ;        //! flow evt maker
    AliFlowWeighter*  fFlowWeighter ;     //! flow phi weights
    AliFlowAnalyser*  fFlowAnal ;         //! flow analysis

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

  // defalt AliSelector things
    void 	      DeleteKinematicsFile();
    TFile*            fKineFile;          //! pointer to Kinematics.root if the file was opened

    AliSelectorFoF(const AliSelectorFoF&);
    AliSelectorFoF& operator=(const AliSelectorFoF&);

  ClassDef(AliSelectorFoF,0);
};

#endif
