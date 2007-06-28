#ifndef ALIMUONTRIGGERCHAMBEREFF_H
#define ALIMUONTRIGGERCHAMBEREFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup base
/// \class AliMUONTriggerChamberEff
/// \brief trigger chamber efficiency from data
///
//  Author Diego Stocco

#include <TObject.h>
#include <TList.h>

class AliMUONGeometryTransformer;
class AliMUONDigitMaker;
class AliMUONTriggerTrack;
class AliMUONVDigitStore;
class AliMUONVTriggerStore;
class AliMUONVTriggerTrackStore;
class AliMUONVTrackStore;
class TClonesArray;

class AliMUONTriggerChamberEff : public TObject
{
public:
    AliMUONTriggerChamberEff();
    AliMUONTriggerChamberEff(const AliMUONGeometryTransformer* transformer,
			     const AliMUONDigitMaker* digitMaker,
			     Bool_t writeOnESD=kFALSE);
    
    virtual ~AliMUONTriggerChamberEff();

    AliMUONTriggerChamberEff(const AliMUONTriggerChamberEff& other); 
    AliMUONTriggerChamberEff& operator=(const AliMUONTriggerChamberEff& other);
    
    /// Set Reproduce trigger response
    void SetReproduceTrigResponse(Bool_t reproduceTrigRes=kFALSE)
    {fReproduceTrigResponse=reproduceTrigRes;}
    /// Set Print informations on event
    void SetPrintInfo(Bool_t printInfo=kFALSE)
    {fPrintInfo=printInfo;}
    /// Set Debug level
    void SetDebugLevel(Int_t debugLevel)
    {fDebugLevel=debugLevel;}

    void EventChamberEff(const AliMUONVTriggerStore& triggerStore,
			 const AliMUONVTriggerTrackStore& trigTrackStore,
			 const AliMUONVTrackStore& trackStore);
    void WriteEfficiencyMap(const char* outputDir);
    void WriteEfficiencyMapTxt(const char* outputDir);

    
private:
    /*
   /// Not implemented
    AliMUONTriggerChamberEff(const AliMUONTriggerChamberEff& other);
    /// Not implemented
    AliMUONTriggerChamberEff& operator=(const AliMUONTriggerChamberEff& other);
    */
    
    static const Int_t fgkNchambers=4; ///< Number of trigger chambers
    static const Int_t fgkNcathodes=2; ///< Number of cathodes per chamber
    static const Int_t fgkNslats=18;   ///< Number of slats per chamber
    static const Int_t fgkNboards=234; ///< Number of trigger boards per chamber

    Int_t fTrigger34[fgkNchambers][fgkNcathodes]; ///< Array counting # of times chamber was inefficient
    Int_t fTrigger44[fgkNcathodes]; ///< Array counting # of times all chambers were efficient
    Int_t fInefficientSlat[fgkNchambers][fgkNcathodes][fgkNslats]; ///< Array counting # of times slats were inefficient
    Int_t fHitPerSlat[fgkNchambers][fgkNcathodes][fgkNslats]; ///< Array counting # of times slats were efficient
    Int_t fInefficientBoard[fgkNchambers][fgkNcathodes][fgkNboards]; ///< Array counting # of times boards were inefficient
    Int_t fHitPerBoard[fgkNchambers][fgkNcathodes][fgkNboards]; ///< Array counting # of times boards were efficient
    
    const AliMUONGeometryTransformer* fTransformer; //!< geometry transformer
    const AliMUONDigitMaker* fDigitMaker; //!< pointer to digit maker
    Bool_t fReproduceTrigResponse; //!< Reproduce trigger response
    Bool_t fPrintInfo; //!< Print informations on event
    Int_t fWriteOnESD; //!< flag to write on ESD
    Int_t fDebugLevel; //!< Debug level
    const Float_t fkMaxDistance; //!< Maximum distance for reference

    
protected:
    Int_t MatchingPad(AliMUONVDigitStore& digitStore, Int_t &detElemId, Float_t coor[2],
		      Bool_t isMatch[fgkNcathodes], Int_t nboard[fgkNcathodes][4],
		      Float_t zRealMatch[fgkNchambers], Float_t y11);
    Float_t PadMatchTrack(Float_t xPad, Float_t yPad, Float_t dpx, Float_t dpy,
			  Float_t xTrackAtPad, Float_t yTrackAtPad, Int_t chamber);
    void InfoDigit(AliMUONVDigitStore& digitStore);
    void CalculateEfficiency(Int_t trigger44, Int_t trigger34, Float_t &efficiency,
			     Float_t &error, Bool_t failuresAsInput);
    Int_t DetElemIdFromPos(Float_t x, Float_t y, Int_t chamber, Int_t cathode);
    void LocalBoardFromPos(Float_t x, Float_t y, Int_t detElemId,
			   Int_t cathode, Int_t localBoard[4]);
    void ResetArrays();
    Bool_t TriggerDigits(const AliMUONVTriggerStore& triggerStore,
			 AliMUONVDigitStore& digitStore) const;
    Bool_t IsCleanTrack(AliMUONTriggerTrack *triggerTrack,
			const AliMUONVTrackStore& trackStore);
    void SaveInESDFile();

    ClassDef(AliMUONTriggerChamberEff,1) // Trigger chamber efficiency
};
#endif
