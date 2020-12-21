/**  
 *  
 *  \class AliAnalysisTaskSimSpectraLF
 *  
 *  \brief Analysis task to be used for on-the-fly simulations in LEGO trains
 * 
 *  Minimal analysis task meant to be used for on-the-fly simulations in LEGO trains (for generating light flavor particle pt spectra)

 *  \author:	Gyula Bencedi  <Gyula.Bencedi@cern.ch>, WIGNER RCP
 * 
 *  \date: 	June 4, 2017
 * 
*/
#ifndef ALIANALYSISTASKSIMSPECTRALF_H
#define ALIANALYSISTASKSIMSPECTRALF_H
 
// ROOT includes

#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <TArrayD.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjArray.h>

// AliRoot includes

#include <AliAnalysisTaskSE.h>
#include <AliMCEvent.h>
#include <AliStack.h>
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"

class TList;
class TH1F;
class TH1I;
class TParticle;
class AliStack;
class AliVVertex;
class AliVParticle;

class AliAnalysisTaskSimSpectraLF : public AliAnalysisTaskSE { //

 public:
   
  AliAnalysisTaskSimSpectraLF();
  AliAnalysisTaskSimSpectraLF(const char *name);

  virtual ~AliAnalysisTaskSimSpectraLF();

  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  void SetYRange(Float_t y){ fY=y; }

 private:

    void	EventSel(TObject* obj);
    void 	ParticleSel(TObject *obj);
    Int_t	GetPidCode(Int_t pdgCode) const;

 protected:
  
    Bool_t IsMCEventSelected(TObject* obj);
    
    AliMCEvent*              fMcEvent;    //!<! MC event
    AliInputEventHandler*    fMcHandler;  //!<! MCEventHandler
    
    AliStack* 	fStack;

    Float_t	fY;     	///< rapidity cut

    void FillHisto(const char* objkey, Double_t x);

    template <class T> T* InitHisto(const char* hname = "hname", const char* htitle = "htitle", Int_t nxbins = 100, Double_t xmin = 0., Double_t xmax = 20., const char* xtitle = "xtitle", const char* ytitle = "ytitle");
    
    TH1I* fHistEvt;	 	//!<! 	QA of event properties
    TH1I* fHistPart;	 	//!<! 	QA of particle properties
    
    TList* fListOfObjects;	//!<! Output list of objects
    
    
    AliAnalysisTaskSimSpectraLF(const AliAnalysisTaskSimSpectraLF&);            // not implemented
    AliAnalysisTaskSimSpectraLF& operator=(const AliAnalysisTaskSimSpectraLF&); // not implemented

    ClassDef(AliAnalysisTaskSimSpectraLF, 1);	// Analysis task for LF spectra analysis  
};

#endif
