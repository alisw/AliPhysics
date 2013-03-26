////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderAOD - the reader class for the Alice AOD                //
// Reads in AOD information and converts it into internal AliFemtoEvent       //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOEVENTREADERAOD_H
#define ALIFEMTOEVENTREADERAOD_H
#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "TChain.h"
#include "TBits.h"
#include "AliAODEvent.h"
#include <list>
//#include "AliPWG2AODTrack.h"
#include "AliAODMCParticle.h"
#include "AliFemtoV0.h"
#include "AliAODpidUtil.h"
#include "AliAODHeader.h"

class AliFemtoEvent;
class AliFemtoTrack;

class AliFemtoEventReaderAOD : public AliFemtoEventReader 
{
 public:
  enum EventMult {kCentrality=0, kGlobalCount=1, kReference=2, kTPCOnlyRef=3, kVZERO=4, kCentralityTRK=5, kCentralityZNA=6, kCentralityCL1=7, kCentralityCND=9};
  typedef enum EventMult EstEventMult;

  AliFemtoEventReaderAOD();
  AliFemtoEventReaderAOD(const AliFemtoEventReaderAOD &aReader);
  virtual ~AliFemtoEventReaderAOD();

  AliFemtoEventReaderAOD& operator=(const AliFemtoEventReaderAOD& aReader);

  virtual AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  void SetInputFile(const char* inputfile);
  void SetFilterBit(UInt_t ibit);
  void SetFilterMask(int ibit);
  void SetReadMC(unsigned char a);
  void SetReadV0(unsigned char a);
  void SetCentralityPreSelection(double min, double max);
  void SetNoCentrality(bool anocent);
  void SetAODpidUtil(AliAODpidUtil *aAODpidUtil);
  void SetAODheader(AliAODHeader *aAODheader);
  void SetMagneticFieldSign(int s);
  void SetEPVZERO(Bool_t);
  void GetGlobalPositionAtGlobalRadiiThroughTPC(AliAODTrack *track, Float_t bfield, Float_t globalPositionsAtRadii[9][3]);
  void SetUseMultiplicity(EstEventMult aType);

 protected:
  virtual void CopyAODtoFemtoEvent(AliFemtoEvent *tEvent);
  virtual void CopyAODtoFemtoTrack( AliAODTrack *tAodTrack,  
				    AliFemtoTrack *tFemtoTrack
				    //				    AliPWG2AODTrack *tPWG2AODTrack
				    ); 
  virtual void CopyAODtoFemtoV0(AliAODv0 *tAODv0, AliFemtoV0 *tFemtoV0);			    
  virtual void CopyPIDtoFemtoTrack( AliAODTrack *tAodTrack, 
				   AliFemtoTrack *tFemtoTrack);

  int            fNumberofEvent;    // number of Events in AOD file
  int            fCurEvent;         // number of current event
  AliAODEvent*   fEvent;            // AOD event
  TBits          fAllTrue;          // Bit set with all true bits
  TBits          fAllFalse;         // Bit set with all false bits
  UInt_t         fFilterBit;        // Bitmap bit for AOD filters
  UInt_t         fFilterMask;       
  //  TClonesArray*  fPWG2AODTracks;    // Link to PWG2 specific AOD information (if it exists)
  
  unsigned char  fReadMC;           // Attempt to read the MC information from the AOD
  unsigned char  fReadV0;           // Read V0 information from the AOD and put it into V0Collection
  unsigned char  fUsePreCent;       // Use centrality pre-selection to speed up analysis
  EstEventMult   fEstEventMult;  // Type of the event multiplicity estimator
  double         fCentRange[2];     // Centrality pre-selection range
  AliAODpidUtil* fAODpidUtil;
  AliAODHeader* fAODheader; 

 private:

  AliAODMCParticle* GetParticleWithLabel(TClonesArray *mcP, Int_t aLabel); 

  string         fInputFile;        // name of input file with AOD filenames
  string         fFileName;         // name of current AOD file
  TChain*        fTree;             // AOD tree
  TFile*         fAodFile;          // AOD file 
  int            fMagFieldSign;     // Magnetic field sign
  Bool_t fisEPVZ; // to get event plane angle from VZERO


#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderAOD, 11)
#endif

};
  
#endif


