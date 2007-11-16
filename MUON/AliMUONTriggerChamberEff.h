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
#include <TArrayI.h>
#include <TArrayF.h>
#include <TH3.h>

class AliMUONGeometryTransformer;
class AliMUONDigitMaker;
class AliMUONTriggerTrack;
class AliMUONVDigitStore;
class AliMUONVTriggerStore;
class AliMUONVTriggerTrackStore;
class AliMUONVTrackStore;

class AliMUONTriggerChamberEff : public TObject
{
public:
    AliMUONTriggerChamberEff();
    AliMUONTriggerChamberEff(const AliMUONGeometryTransformer* transformer,
			     const AliMUONDigitMaker* digitMaker,
			     Bool_t writeOnESD=kFALSE);
    
    virtual ~AliMUONTriggerChamberEff();

    AliMUONTriggerChamberEff(const AliMUONTriggerChamberEff& other); // copy constructor
    AliMUONTriggerChamberEff& operator=(const AliMUONTriggerChamberEff& other); // assignment operator
    
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

protected:
    Int_t MatchingPad(AliMUONVDigitStore& digitStore, Int_t &detElemId, Float_t coor[2],
		      Bool_t isMatch[2], TArrayI nboard[2],
		      TArrayF &zRealMatch, Float_t y11);
    Float_t PadMatchTrack(Float_t xPad, Float_t yPad, Float_t dpx, Float_t dpy,
			  Float_t xTrackAtPad, Float_t yTrackAtPad, Int_t chamber);
    void InfoDigit(AliMUONVDigitStore& digitStore);
    void CalculateEfficiency(Int_t trigger44, Int_t trigger34, Float_t &efficiency,
			     Float_t &error, Bool_t failuresAsInput);
    Int_t DetElemIdFromPos(Float_t x, Float_t y, Int_t chamber, Int_t cathode);
    void LocalBoardFromPos(Float_t x, Float_t y, Int_t detElemId,
			   Int_t cathode, Int_t localBoard[4]);
    void ResetArrays();
    void InitHistos();
    Bool_t TriggerDigits(const AliMUONVTriggerStore& triggerStore,
			 AliMUONVDigitStore& digitStore) const;
    Bool_t IsCleanTrack(AliMUONTriggerTrack *triggerTrack,
			const AliMUONVTrackStore& trackStore);
    void SaveInESDFile();
    void GetEfficiencyHistos(TList &countList, TList &noCountList);

    
private:
    void CheckConstants() const;
    /// Get max number of strips along x
    inline Int_t GetMaxX(Int_t cath){return (cath==0) ? 7 : 112;}
    /// Get max number of strips along x
    inline Int_t GetMaxY(Int_t cath){return (cath==0) ? 64 : 1;}

    const AliMUONGeometryTransformer* fTransformer; //!< geometry transformer
    const AliMUONDigitMaker* fDigitMaker; //!< pointer to digit maker
    Bool_t fReproduceTrigResponse; //!< Reproduce trigger response
    Bool_t fPrintInfo; //!< Print informations on event
    Int_t fWriteOnESD; //!< flag to write on ESD
    Int_t fDebugLevel; //!< Debug level
    const Float_t fkMaxDistance; //!< Maximum distance for reference

    static const Int_t fgkNcathodes=2; ///<Number of cathodes
    static const Int_t fgkNchambers=4; ///<Number of chambers
    static const Int_t fgkNplanes=8;   ///<Number of planes
    static const Int_t fgkNslats=18;   ///<Number of slats
    static const Int_t fgkNlocations=4; ///<Number of locations

    TArrayI fTrigger44; ///< Array counting # of times all chambers were efficient
    TArrayI fTrigger34; ///< Array counting # of times chamber was inefficient
    TArrayI fInefficientSlat[fgkNplanes]; ///< Array counting # of times slats were inefficient
    TArrayI fHitPerSlat[fgkNplanes]; ///< Array counting # of times slats were efficient
    TArrayI fInefficientBoard[fgkNplanes]; ///< Array counting # of times boards were inefficient
    TArrayI fHitPerBoard[fgkNplanes]; ///< Array counting # of times boards were efficient

    TH3F *fPadFired[fgkNcathodes]; ///< Histo counting the fired pads
    
    ClassDef(AliMUONTriggerChamberEff,3) // Trigger chamber efficiency
};
#endif
