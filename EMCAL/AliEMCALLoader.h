#ifndef ALIEMCALLOADER_H
#define ALIEMCALLOADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A singleton that returns various objects 
//  Should be used on the analysis stage to avoid confusing between different
//  branches of reconstruction tree: e.g. reading RecPoints and TS made from 
//  another set of RecPoints.
// 
//  The objects are retrived from folders.  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//    


// --- ROOT system ---
#include "TClonesArray.h"
#include "TFolder.h"  
#include "TTree.h"
class TString ;
class TParticle ;
class TTask ;

// --- Standard library ---

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALTrackSegment.h"
#include "AliEMCALClusterizer.h"
#include "AliEMCALTrackSegmentMaker.h" 
#include "AliEMCALPID.h"
class AliLoader ;
class AliEMCAL ; 
class AliEMCALHit ;
class AliEMCALRecParticle ; 
class AliEMCALGeometry ;
class AliEMCALDigitizer ; 
class AliEMCALSDigitizer ; 
class AliEMCALCalibrationDB ;


//

class AliEMCALLoader : public AliLoader {
  
 public:

  AliEMCALLoader();
  AliEMCALLoader(const AliEMCALLoader & obj):AliLoader(obj){}
  AliEMCALLoader(const Char_t *detname,const Char_t *eventfoldername); 
  
  virtual ~AliEMCALLoader() ; 

  // assignement operator requested by coding convention, but not needed
  const AliEMCALLoader & operator = (const AliEMCALLoader & ) {return *this;}

  Int_t   GetEvent();//extends the method on EMCAL RecPart posting
  Int_t   SetEvent();//extends the method on EMCAL RecPart posting
  
  Bool_t  BranchExists(const TString& recName);
  Int_t   LoadHits(Option_t* opt=""); //reads  from disk and sends them to folder; array as well as tree
  Int_t   LoadSDigits(Option_t* opt="");
  Int_t   LoadDigits(Option_t* opt=""); //reads Digits from disk and sends them to folder; array as well as tree
  Int_t   LoadRecPoints(Option_t* opt=""); //reads RecPoints from disk and sends them to folder; array as well as tree
  Int_t   LoadTracks(Option_t* opt="");  //reads Tracks from disk and sends them to folder; array as well as tree
  Int_t   LoadRecParticles(Option_t* opt="");
 
  void    UnloadRecParticles();
  void    UnloadTracks();
  
  Int_t   PostHits();  //Posts the 
  Int_t   PostSDigits();
  Int_t   PostDigits();
  Int_t   PostRecPoints();
  Int_t   PostTracks();
  Int_t   PostRecParticles();
  
  void    CleanFolders();//cleans all the stuff loaded by this detector + calls AliLoader::Clean

  void    CleanHits();
  void    CleanSDigits();
  void    CleanDigits();
  void    CleanRecPoints();
  void    CleanTracks();
  void    CleanRecParticles();

//up to now it is only here -> no definition about global/incremental tracking/PID
 
//   Int_t   WriteRecParticles(Option_t* opt="");//writes the reconstructed particles
//   Int_t   WritePID(Option_t* opt="");//writes the task for PID to file
//   Bool_t  PostPID  (AliEMCALPID * pid) const {return kTRUE;}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

  TObject** HitsRef(){return GetDetectorDataRef(Hits());}
  TObject** SDigitsRef(){return GetDetectorDataRef(SDigits());}
  TObject** DigitsRef(){return GetDetectorDataRef(Digits());}
  TObject** ECARecPointsRef(){return GetDetectorDataRef(ECARecPoints());}
  TObject** TracksRef(){return GetDetectorDataRef(TrackSegments());}
  TObject** RecParticlesRef(){return GetDetectorDataRef(RecParticles());}

  void   Track(Int_t itrack) ;

  static AliEMCALGeometry* GetEMCALGeometry();
  static AliEMCALLoader* GetEMCALLoader(const  char* eventfoldername);

  //Method to be used when digitizing under AliRunDigitizer, who opens all files etc.
  Int_t  EventNumber()       { return (Int_t) GetRunLoader()->GetEventNumber();}
  Int_t  MaxEvent()          { return (Int_t) GetRunLoader()->TreeE()->GetEntries();}

  const AliEMCAL *         EMCAL();
  const AliEMCALGeometry  *EMCALGeometry() ; 

  /*********************************************/
  /************    TClonesArrays     ***********/
  /*********************************************/
  /****   H i t s  ****/
  TClonesArray*  Hits(void);
  const AliEMCALHit*    Hit(Int_t index);
  void MakeHitsArray();
  /****   S D i g i t s  ****/ 
  TClonesArray*  SDigits();
  const AliEMCALDigit*  SDigit(Int_t index);
  void MakeSDigitsArray();
  /****  D i g i t s  ****/
  TClonesArray*   Digits();
  const AliEMCALDigit *  Digit(Int_t index);
  void MakeDigitsArray();
  /****  R e c P o i n t s  ****/
  TObjArray * ECARecPoints();
  const AliEMCALRecPoint * ECARecPoint(Int_t index) ;
  void MakeRecPointsArray();
  /****   T r a c k S e g m e n t s ****/
  TClonesArray * TrackSegments();
  const AliEMCALTrackSegment * TrackSegment(Int_t index);
  void MakeTrackSegmentsArray();
  /****  R e c P a r t ic l e s   ****/
  TClonesArray * RecParticles() ;
  const AliEMCALRecParticle * RecParticle(Int_t index);
  void MakeRecParticlesArray();

  /*********************************************/
  /************    T A S K S      **************/
  /*********************************************/
  // 
  //  AliEMCALSDigitizer*  EMCALSDigitizer(TString name = AliConfig::GetDefaultEventFolderName());
  //AliEMCALDigitizer*   EMCALDigitizer()  { return  dynamic_cast<AliEMCALDigitizer*>(Digitizer()) ;}

  AliEMCALClusterizer* Clusterizer () const {return dynamic_cast<AliEMCALClusterizer*>(Reconstructioner()) ;}
  Int_t PostClusterizer(TTask* clust) const {return PostReconstructioner(clust);}
  Int_t LoadClusterizer(Option_t * opt="") const {return LoadReconstructioner(opt);}
  Int_t WriteClusterizer(Option_t * opt="") const {return WriteReconstructioner(opt);}

  AliEMCALPID * PID() const {return dynamic_cast<AliEMCALPID*>(PIDTask()) ;}
  Int_t PostPID(TTask* pid) const {return PostPIDTask(pid);}
  Int_t LoadPID(Option_t * opt="") const {return LoadPIDTask(opt);}
  Int_t WritePID(Option_t * opt="") const {return WritePIDTask(opt);}


  AliEMCALTrackSegmentMaker * TrackSegmentMaker () const  { return dynamic_cast<AliEMCALTrackSegmentMaker *>(Tracker()) ;}
  Int_t PostTrackSegmentMaker(TTask* segmaker) const {return PostTracker(segmaker);}
  Int_t LoadTrackSegmentMaker(Option_t * opt="") const {return LoadTracker(opt);}
  Int_t WriteTrackSegmentMaker(Option_t * opt="") const  {return WriteTracker(opt);}

  
  void   SetDebug(Int_t level) {fDebug = level;} // Set debug level
  void   SetBranchTitle(const TString& btitle);
  
  AliEMCALCalibrationDB * CalibrationDB(){return  fcdb; }
  //void ReadCalibrationDB(const char * name, const char * filename);
  

  static TString HitsName() { return fgkHitsName ; }   //Name for TClonesArray with hits from one event
  static TString SDigitsName() { return fgkSDigitsName ;} //Name for TClonesArray 
  static TString DigitsName() { return fgkDigitsName ;} //Name for TClonesArray 
  static TString ECARecPointsName() { return fgkECARecPointsName ;} //Name for TClonesArray y 
  static TString TracksName() { return fgkTracksName ;} //Name for TClonesArray 
  static TString RecParticlesName() { return fgkRecParticlesName ;} //Name for TClonesArray
  static TString ECARecPointsBranchName() { return fgkECARecPointsBranchName ;} //Name for branch
  static TString TrackSegmentsBranchName() { return fgkTrackSegmentsBranchName ;} //Name for branch
  static TString RecParticlesBranchName() { return fgkRecParticlesBranchName ;} //Name for branch

protected:
  TString fBranchTitle;            //Title of the branch
  Bool_t  fRecParticlesLoaded;     //Flag signing if Reconstructed Particles are loaded
  Bool_t  fTracksLoaded;           //Flag signing if Tracks are loaded
  TString fRecParticlesFileOption; //Loading Option for Reconstructed Particles
  AliEMCALCalibrationDB * fcdb ;       //!

private:

  Int_t ReadHits();
  Int_t ReadDigits();
  Int_t ReadSDigits();
  Int_t ReadRecPoints();
  Int_t ReadTracks();
  Int_t ReadRecParticles();

  static const TString fgkHitsName;//Name for TClonesArray with hits from one event
  static const TString fgkSDigitsName;//Name for TClonesArray 
  static const TString fgkDigitsName;//Name for TClonesArray 
  static const TString fgkECARecPointsName;//Name for TClonesArray 
  static const TString fgkTracksName;//Name for TClonesArray 
  static const TString fgkRecParticlesName;//Name for TClonesArray

  static const TString fgkECARecPointsBranchName;//Name for branch
  static const TString fgkTrackSegmentsBranchName;//Name for branch
  static const TString fgkRecParticlesBranchName;//Name for branch
  Int_t  fDebug ;             // Debug level
 

  
  ClassDef(AliEMCALLoader,4)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

/******************************************************************************/
/****************    I N L I N E S     ****************************************/
/******************************************************************************/

inline TClonesArray* AliEMCALLoader::Hits()  
{
 return (TClonesArray*)GetDetectorData(fgkHitsName);
}
/******************************************************************************/

inline const AliEMCALHit* AliEMCALLoader::Hit(Int_t index)  
{
  const TClonesArray* tcarr = Hits();
  if (tcarr)
    return (const AliEMCALHit*) tcarr->At(index);
  return 0x0; 
}
/******************************************************************************/

inline TClonesArray* AliEMCALLoader::SDigits()
{
   return dynamic_cast<TClonesArray*>(GetDetectorData(fgkSDigitsName));
}
/******************************************************************************/

inline const AliEMCALDigit*  AliEMCALLoader::SDigit(Int_t index)
{
  const TClonesArray* tcarr = SDigits();
  if (tcarr)
    return (const AliEMCALDigit*) tcarr->At(index);
  return 0x0; 
}
/******************************************************************************/

inline TClonesArray* AliEMCALLoader::Digits()
{
 return dynamic_cast<TClonesArray*>(GetDetectorData(fgkDigitsName));
}
/******************************************************************************/

inline const AliEMCALDigit*  AliEMCALLoader::Digit(Int_t index)
{
  const TClonesArray* tcarr = Digits();
  if (tcarr)
    return (const AliEMCALDigit*) tcarr->At(index);
  return 0x0; 
}

/******************************************************************************/

inline TObjArray * AliEMCALLoader::ECARecPoints()
{
 return dynamic_cast<TObjArray*>(GetDetectorData(fgkECARecPointsName));
}

/******************************************************************************/

inline const AliEMCALRecPoint * AliEMCALLoader::ECARecPoint(Int_t index)
{
  TObjArray* tcarr = ECARecPoints();
  if (tcarr)
    return dynamic_cast<const AliEMCALRecPoint*>(tcarr->At(index));
  return 0x0; 
}

/******************************************************************************/

inline TClonesArray * AliEMCALLoader::TrackSegments()
{
 return dynamic_cast<TClonesArray*>(GetDetectorData(fgkTracksName));
}
/******************************************************************************/

inline const AliEMCALTrackSegment * AliEMCALLoader::TrackSegment(Int_t index)
{
  const TClonesArray* tcarr = TrackSegments();
  if (tcarr)
    return (const AliEMCALTrackSegment*) tcarr->At(index);
  return 0x0; 
}
/******************************************************************************/

inline TClonesArray * AliEMCALLoader::RecParticles() 
{
 return dynamic_cast<TClonesArray*>(GetDetectorData(fgkRecParticlesName)); 
}
/******************************************************************************/

inline const AliEMCALRecParticle* AliEMCALLoader::RecParticle(Int_t index)
{
  TClonesArray* tcarr = RecParticles();
  if (tcarr)
    return (const AliEMCALRecParticle*) tcarr->At(index);
  return 0x0;  
}

#endif // AliEMCALLOADER_H
