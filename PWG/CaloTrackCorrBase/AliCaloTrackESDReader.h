#ifndef ALICALOTRACKESDREADER_H
#define ALICALOTRACKESDREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloTrackESDReader
/// \ingroup CaloTrackCorrelationsBase
/// \brief Class for event, clusters and tracks filtering and preparation for the ESD analysis.
///
/// Class for accessing/filtering ESD data. Most of the job is done in the mother class
/// here only very specific methods of the ESD format are implemented.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________


class AliESDEvent;

#include "AliCaloTrackReader.h" 

class AliCaloTrackESDReader : public AliCaloTrackReader {
  
public:
  
                   AliCaloTrackESDReader() ; // ctor
  
  virtual         ~AliCaloTrackESDReader() ; // virtual dtor

  Bool_t           CheckForPrimaryVertex() const ;
  
  void             Init();
  
  AliGenEventHeader* GetGenEventHeader()   const ;

  Bool_t           SelectTrack(AliVTrack* track, Double_t* pTrack);
  
  AliESDtrackCuts* GetTrackCuts()                    const { return fESDtrackCuts     ; }
  void             SetTrackCuts(AliESDtrackCuts * cuts) ;
  
  AliESDtrackCuts* GetTrackComplementaryCuts()       const { return fESDtrackComplementaryCuts ; }
  void             SetTrackComplementaryCuts(AliESDtrackCuts * cuts)  ;

  void             SwitchOnConstrainTrackToVertex()        { fConstrainTrack = kTRUE  ; }
  void             SwitchOffConstrainTrackToVertex()       { fConstrainTrack = kFALSE ; }
  
  void             SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ;
	 
private:
  
  Bool_t           fConstrainTrack;            ///< Constrain Track to vertex.
  AliESDtrackCuts* fESDtrackCuts ;             ///< Track cut machinery.
  AliESDtrackCuts* fESDtrackComplementaryCuts; ///< Track cut machinery for complementary cuts for hybrids.
  
  /// Copy constructor not implemented.
  AliCaloTrackESDReader(              const AliCaloTrackESDReader & r) ; 
  
  /// Assignment operator not implemented.
  AliCaloTrackESDReader & operator = (const AliCaloTrackESDReader & r) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliCaloTrackESDReader,2) ;
  /// \endcond 

} ;

#endif //ALICALOTRACKESDREADER_H



