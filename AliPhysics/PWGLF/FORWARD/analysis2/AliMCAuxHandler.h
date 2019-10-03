/**
 * @file   AliMCAuxHandler.h
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 12:04:02 2013
 * 
 * @brief  
 * 
 * 
 */
#ifndef ALIFMDMCHITHANDLER_H
#define ALIFMDMCHITHANDLER_H
#include <AliMCEventHandler.h>
class TFile;
class TTree;

/**
 * This class defines an input handler for simulated data which will
 * connect the FMD Hit tree.  It is intended to be added to the global
 * MC input handler using AliMCEventHandler::AddSubsidiaryHandler
 * 
 */
class AliMCAuxHandler : public AliMCEventHandler
{
public:
  /** 
   * Constructor 
   * 
   * @param name    Name 
   * @param clsName Class name of hits 
   * @param parent  Parent event handler 
   */
  AliMCAuxHandler(const char* name="FMD",
		     const char* clsName="AliFMDHit",
		     AliMCEventHandler* parent=0);
  /** Destructor */
  virtual ~AliMCAuxHandler() {}
  /** 
   * @{
   * @name Interface member functions 
   */
  /** 
   * Intialize 
   * 
   * @param t Not used
   * @param o Not used
   * 
   * @return always true 
   */  
  virtual Bool_t Init(TTree* t,Option_t* o) { return AliMCEventHandler::Init(t,o); }
  /** 
   * Initialize the input
   * 
   * @param opt Options 
   * 
   * @return true on success 
   */
  virtual Bool_t Init(Option_t* opt);
  /**
   * Called at the beginning of an event 
   * 
   * @param entry Entry in tree 
   * 
   * @return true on success
   */
  virtual Bool_t BeginEvent(Long64_t entry);
  /** 
   * Called when the input file is changed 
   * 
   * @return true on success
   */
  virtual Bool_t Notify() { return AliMCEventHandler::Notify(); }
  /** 
   * Called when the input file is changed 
   * 
   * @param path New path 
   *
   * @return true on success
   */
  virtual Bool_t Notify(const char* path);
  /** 
   * Called at the end of an event 
   * 
   * @return true on success
   */
  virtual Bool_t FinishEvent();
  /** 
   * Called at the end of a job 
   * 
   * @return true on success 
   */
  virtual Bool_t Terminate();
  /** 
   * Called at the end of a sub-job
   * 
   * @return true on success
   */  
  virtual Bool_t TerminateIO();
  /** 
   * Reset the I/O
   * 
   */
  virtual void ResetIO();
  /** 
   * Load an event 
   * 
   * @param iev Event number 
   * 
   * @return true on success 
   */
  virtual Bool_t LoadEvent(Int_t iev);
  /** 
   * Set the number of events in the container 
   * 
   * @param nev Number of events 
   */
  virtual void  SetNumberOfEventsInContainer(Int_t nev) {
    this->fNEventsInContainer = nev;}
  /** 
   * Open a file 
   * 
   * @param ev event number 
   * 
   * @return true on success
   */
  virtual Bool_t OpenFile(Int_t ev);
  /* @} */

  /** 
   * Get the parent handler 
   * 
   * @return Parent handler
   */
  AliMCEventHandler* GetParent() { return fParent; }
  /** 
   * Get the tree 
   * 
   * @return The connected hits tree
   */
  virtual TTree*  GetTree() const { return fTree;}  
  /** 
   * Get array of hits 
   * 
   * @return Array of hits
   */
  TClonesArray*  GetArray() const { return fArray; }
  /** 
   * Get the number of entries
   * 
   * @return Entries in array 
   */
  Int_t         GetNEntry() const;
  /** 
   * Get array of track-refs for an entry 
   * 
   * @param entry Entry number 
   * 
   * @return Array 
   */
  TClonesArray* GetEntryArray(Int_t entry);
  /** 
   * Static member function to create and attach this handler 
   * 
   * @param name Name of the handler 
   * @param what What to get 
   * 
   * @return Newly allocated handler or null
   */
  static AliMCAuxHandler* Create(const char* name="FMD",
				    const char* what="Hits");
  /** 
   * Static member function to get the kinamtics array
   * 
   * @param handler   Unput handler
   * @param particle  Particle number 
   * 
   * @return Array of hits 
   */
  static TClonesArray* GetParticleArray(AliMCAuxHandler* handler, 
					Int_t particle);
protected:
  /** 
   * Copy constructor
   * 
   * @param o Object to copy from
   */
  AliMCAuxHandler(const AliMCAuxHandler& o)
    : AliMCEventHandler(), 
      fParent(o.fParent),
      fFile(0), 
      fTree(0), 
      fDir(0),
      fArray(0),
      fNEvents(0), 
      fNEventsPerFile(0),
      fNEventsInContainer(0),
      fEvent(0), 
      fFileNumber(0), 
      fTreeName(""),
      fFileBase("")
  {}
  /** 
   * Assignment operator
   * 
   * @param o Object to assign from 
   * 
   * @return reference to this object 
   */
  AliMCAuxHandler& operator=(const AliMCAuxHandler& o)
  {
    if (&o == this) return *this;
    // AliMCEventHandler::operator=(o);
    fParent = o.fParent;
    fFile   = o.fFile;
    fTree   = o.fTree;
    return *this;
  }
  /** 
   * Get the parent path 
   * 
   * @return Parent path
   */
  TString* GetParentPath() const;
  AliMCEventHandler* fParent; // Parent MC handler 
  TFile*             fFile;                //!
  TTree*             fTree;                //!
  TDirectory*        fDir;                 //!
  TClonesArray*      fArray;               //!
  Int_t              fNEvents;             //!
  Int_t              fNEventsPerFile;      //!
  Int_t              fNEventsInContainer;  //!
  Int_t              fEvent;               //!
  Int_t              fFileNumber;          //!
  TString            fTreeName;            //!
  TString            fFileBase;            //! 
  ClassDef(AliMCAuxHandler,1); // Connect FMD hits tree
};

#endif
// Local Variables:
//  mode: C++
// End:

