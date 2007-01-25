#ifndef ALIMUONTRIGGERCHAMBEREFF_H
#define ALIMUONTRIGGERCHAMBEREFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup base
/// \class AliMUONTriggerChamberEff
/// \brief trigger chamber efficiency from data
/// \author Diego Stocco

#include <TObject.h>
#include <TString.h>

class AliRunLoader;
class AliMUONData;
class AliMUON;
class AliMUONGlobalTrigger;
class TString;

class AliMUONTriggerChamberEff : public TObject
{
public:
    AliMUONTriggerChamberEff(const char* galiceFile, Int_t firstEvent=0, Int_t lastEvent=-1);
    AliMUONTriggerChamberEff(Int_t firstRun, Int_t lastRun, const char* galiceRunDir, Int_t firstEvent=0, Int_t lastEvent=-1);
    virtual ~AliMUONTriggerChamberEff();

    void SetReproduceTrigResponse(Bool_t reproduceTrigRes=kFALSE)
    {fReproduceTrigResponse=reproduceTrigRes;}
    void SetPrintInfo(Bool_t printInfo=kFALSE)
    {fPrintInfo=printInfo;}
    void SetDebugLevel(Int_t debugLevel)
    {fDebugLevel=debugLevel;}

    void PerformTriggerChamberEff(const char* outputDir);
    
protected:
    Bool_t PadMatchTrack(Float_t xPad, Float_t yPad, Float_t dpx, Float_t dpy,
			 Float_t xTrackAtPad, Float_t yTrackAtPad, Int_t chamber);
    Bool_t IsDiffLocalBoard(Int_t currDetElemId, Int_t iy, Int_t detElemIdP1, Int_t iyDigitP1) const;
    void PrintTrigger(AliMUONGlobalTrigger *globalTrig);
    void InfoDigit();
    void CalculateEfficiency(Int_t trigger44, Int_t trigger34, Float_t &efficiency, Float_t &error, Bool_t failuresAsInput);
    Int_t DetElemIdFromPos(Float_t x, Float_t y, Int_t chamber, Int_t cathode);
    void LocalBoardFromPos(Float_t x, Float_t y, Int_t detElemId, Int_t cathode, Int_t localBoard[4]);
    void ResetArrays();
    void WriteOutput(const char* outputDir, Int_t totNumOfTrig[4][2], Int_t atLeast1MuPerEv[4][2]);
    void WriteEfficiencyMap(const char* outputDir);
    
private:
    AliMUONTriggerChamberEff(const AliMUONTriggerChamberEff& other);
    AliMUONTriggerChamberEff& operator=(const AliMUONTriggerChamberEff& other);
    
    void SetGaliceFile(const char* galiceFile);
    void CleanGalice();
    
    Int_t   fFirstEvent; //!< First event to consider
    Int_t   fLastEvent;  //!< Last event to consider
    Int_t   fFirstRun; //!< First run to consider
    Int_t   fLastRun;  //!< Last run to consider
    AliRunLoader* fRunLoader; //!< AliRunLoader pointer
    AliMUONData*  fData; //!< AliMUONData pointer (to access containers)
    Bool_t fReproduceTrigResponse;//!< Reproduce trigger response
    Bool_t fPrintInfo;//!< Print informations on event
    AliMUON *fMUON; //!< AliMUON pointer
    Int_t fDebugLevel; //!< Debug level
    TString fGaliceDir; //!< base directory for many runs.

    static const Int_t fgkNchambers=4; //!< Number of trigger chambers
    static const Int_t fgkNcathodes=2; //!< Number of cathodes per chamber
    static const Int_t fgkNslats=18;   //!< Number of slats per chamber
    static const Int_t fgkNboards=234; //!< Number of trigger boards per chamber

    Int_t fTrigger34[fgkNchambers][fgkNcathodes];//!< Array counting # of times chamber was inefficient
    Int_t fTrigger44[fgkNcathodes];//!< Array counting # of times all chambers were efficient
    Int_t fInefficientSlat[fgkNchambers][fgkNcathodes][fgkNslats];//!< Array counting # of times slats were inefficient
    Int_t fHitPerSlat[fgkNchambers][fgkNcathodes][fgkNslats];//!< Array counting # of times slats were efficient
    Int_t fInefficientBoard[fgkNchambers][fgkNcathodes][fgkNboards];//!< Array counting # of times boards were inefficient
    Int_t fHitPerBoard[fgkNchambers][fgkNcathodes][fgkNboards];//!< Array counting # of times boards were efficient

    ClassDef(AliMUONTriggerChamberEff,0) // Dumper of MUON related data
};
#endif
