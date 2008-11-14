#ifndef AliITSFastOrCalibrationSPD_H
#define AliITSFastOrCalibrationSPD_H

#include <TObject.h>
#include <TBits.h>
#include <AliCDBManager.h>
#include <AliCDBMetaData.h>
#include <AliCDBEntry.h>


class AliITSFastOrCalibrationSPD : public TObject{
 public:
    AliITSFastOrCalibrationSPD(); //default constructor
    virtual ~AliITSFastOrCalibrationSPD(); //destructor

    //setters          
    void   SetFastOrConfiguredChips(UInt_t chipKey) {fFastOrConfiguredChips.SetBitNumber(chipKey);}
    void   ResetFastOrConfiguredChips() {fFastOrConfiguredChips.ResetAllBits();}

    //getters
    TBits  GetFastOrConfiguredChips() const {return fFastOrConfiguredChips;}
    Bool_t TestFastOrConfiguredChips(UInt_t chipKey) const {return fFastOrConfiguredChips.TestBitNumber(chipKey);}

    Bool_t WriteFOConfToDB(Int_t runNrStart, Int_t runNrEnd); 

 private:
    TBits  fFastOrConfiguredChips;     // Map of FastOr configured chips

    ClassDef(AliITSFastOrCalibrationSPD,1) 
};

#endif
