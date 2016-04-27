/**
 * @file   AliBaseMultTask.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Mar 30 13:47:17 2016
 * 
 * @brief  Base task for multiplicity distribution tasks 
 * 
 * 
 */
#ifndef ALIBASEMULTTASK_H
#define ALIBASEMULTTASK_H
#include "AliBaseAODTask.h"

class AliBaseMultTask : public AliBaseAODTask
{
public:

  /** 
   * Base class for eta bins 
   */
  struct Bin : public TNamed
  {
    /**
     * Default Constructor
     */
    Bin();
    /**
     * Constructor
     */
    Bin(Double_t etaLow, Double_t etaHigh);
    /**
     * Copy Constructor
     */ 
    Bin(const Bin&){;}
    /**
     * Assignment Operator
     */
    Bin& operator=(const Bin&) { return *this; }
    /**
     * Destructor
     */
    virtual ~Bin(){}
    /**
     * Define outputs of a single eta bin
     */
    virtual void CreateOutputObjects(TList* cont,  Int_t max);
    /** 
     * Process a single eta bin
     * 
     * @param dndetaForward   Forward observations 
     * @param dndetaCentral   Central observations 
     * @param normForward     Acceptance 
     * @param normCentral     Acceptance 
     * @param mc              Primary "observations"
     * @param ipZ             Interaction point 
     * @param pileup          True if flagged as pile-up
     * @param selectedTrigger Is event selected
     * @param isMCClass       Is event MC NSD 
     * @param isESDClass      Is event real NSD 
     * @param aodevent        Full event 
     * @param minIPz          Least Z coordinate of IP 
     * @param maxIPz          Largest Z coordinate of IP 
     */
    virtual void Process(TH1D*              dndetaForward,
			 TH1D*              dndetaCentral,
			 TH1D*              normForward,
			 TH1D*              normCentral,
			 TH1D*              mc,
			 Double_t           ipZ,
			 Bool_t             pileup,
			 Bool_t             selectedTrigger,
			 Bool_t             isMCClass,
			 Bool_t             isESDClass,
			 const AliAODEvent& aodevent,
			 Double_t           minIPz,
			 Double_t           maxIPz);
    /** 
     * Calculate multiplicity in this bin 
     * 
     * @param dndetaForward Forward observations 
     * @param dndetaCentral Central observations 
     * @param normForward   Acceptance 
     * @param normCentral   Acceptance 
     * @param mc            Primary "observations"
     * @param ipZ           Interaction point 
     * @param statErr       On return, statistical errors 
     * @param sysErr        On return, systematic errors 
     * @param mcMult        On return, true multiplicity (or -1 if not MC)
     * @param mcErr         On return, true stat error (or -1 if not MC) 
     * 
     * @return Event multiplicity in this bin
     */
    virtual Double_t CalcMult(TH1D*              dndetaForward, 
			      TH1D*              dndetaCentral,
			      TH1D*              normForward,   
			      TH1D*              normCentral,
			      TH1D*              mc, 
			      Double_t           ipZ,
			      Double_t&          statErr,
			      Double_t&          sysErr,
			      Double_t&          mcMult, 
			      Double_t&          mcErr);
    /**
     *  Form name of eta bin
     */
    static const Char_t* FormBinName(Double_t etaLow, Double_t etaHigh);
    
    Double_t fEtaLow;                  // low eta limit 
    Double_t fEtaHigh;                 // high eta limit
    TH1D*    fHist;                    //!
    TH1D*    fHistMC;                  //!
    TH2D*    fAcceptance;              //!
    TH2D*    fVtxZvsNdataBins;         //! 
    ClassDef(Bin,1); // Base class for eta bins 
  };

  /**
   * Default constructor for ROOT I/O only
   * 
   */
  AliBaseMultTask();
  /**
   * Default constructor for ROOT I/O only
   * 
   */
  virtual ~AliBaseMultTask() {}
  /** 
   * Add another eta bin to the task
   */
  void AddBin(Double_t etaLow, Double_t etaHigh)
  {
    fBins.Add(MakeBin(etaLow, etaHigh));
  }
  /** 
   * Create a bin.  Must be overloaded in derived class 
   *
   * @param minEta Least pseudorapidity 
   * @param maxEta Largest pseudorapidity 
   *
   * @return newly created bin object 
   */
  virtual Bin* MakeBin(Double_t minEta,Double_t maxEta) = 0;
  /**
   * Create our default bins 
   */
  virtual void DefaultBins();
  /**
   * Create Output Objects
   */
  virtual Bool_t Book();
  /** 
   * Reset before start of evet 
   * 
   * @return Always true
   */
  virtual Bool_t PreEvent() { fIsSelected = false; return true; }
  /**
   * User Exec
   */
  virtual Bool_t Event(AliAODEvent& aod);
  virtual Bool_t Finalize() {return true;}
protected:
  /**
   * NAmed constructor 
   * 
   */
  AliBaseMultTask(const char* name);
  /**
   * Copy constructor 
   * 
   */
  AliBaseMultTask(const AliBaseMultTask&)
    : fBins(),
      fIsSelected(false)
  {}
  /**
   * Assignment operator
   * 
   */
  AliBaseMultTask& operator=(const AliBaseMultTask&) { return *this; }
  /** 
   * Check the event
   * 
   * @param fwd Forwarddata 
   * 
   * @return true on success
   */
  virtual Bool_t CheckEvent(const AliAODForwardMult& fwd);
  /** 
   * Check if event is proper MC class 
   * 
   * @param m Forward object 
   * 
   * @return true if proper MC class 
   */
  virtual Bool_t IsMCClass(AliAODForwardMult* m) const;
  /** 
   * Check if event is proper reconstructed event class 
   * 
   * @param m Forward object 
   * 
   * @return true if proper reconstructed event class 
   */
  virtual Bool_t IsESDClass(AliAODForwardMult* m) const;
  /** 
   * Loop over bins and call Process for each of them 
   * 
   * @param dndetaForward   Forward observations 
   * @param dndetaCentral   Central observations 
   * @param normForward     Acceptance 
   * @param normCentral     Acceptance 
   * @param mc              Primary "observations"
   * @param ipZ             Interaction point 
   * @param pileup          True if flagged as pile-up
   * @param selectedTrigger Is event selected
   * @param isMCClass       Is event MC NSD 
   * @param isESDClass      Is event real NSD 
   * @param aodevent        Full event 
   */
  virtual void Process(TH1D*              dndetaForward,
		       TH1D*              dndetaCentral,
		       TH1D*              normForward,
		       TH1D*              normCentral,
		       TH1D*              mc,
		       Double_t           ipZ,
		       Bool_t             pileup,
		       Bool_t             selectedTrigger,
		       Bool_t             isMCClass,
		       Bool_t             isESDClass,
		       const AliAODEvent& aodevent);
  TList  fBins;
  Bool_t fIsSelected;
  
  ClassDef(AliBaseMultTask,1); // Base class for multiplicity dist tasks,
};
#endif
// Local Variables:
//  mode: C++
// End:

  
