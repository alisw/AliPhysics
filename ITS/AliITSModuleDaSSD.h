#ifndef ALIITSMODULEDASSD_H
#define ALIITSMODULEDASSD_H


#include <iostream>
#include "TObject.h"
#include "AliITSNoiseSSD.h"

#include "AliITSChannelDaSSD.h"

class AliITSModuleDaSSD : public TObject {
  public :
    AliITSModuleDaSSD();
    AliITSModuleDaSSD(const Int_t numberofstrips);
    AliITSModuleDaSSD(const Int_t numberofstrips, const Long_t eventsnumber);
    AliITSModuleDaSSD(const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const UShort_t moduleID);
    AliITSModuleDaSSD(const AliITSModuleDaSSD& module);
    AliITSModuleDaSSD& operator = (const AliITSModuleDaSSD& module);
    virtual ~AliITSModuleDaSSD();
    
    UChar_t      GetDdlId()    const { return fDdlId; }
    UChar_t      GetAD()       const { return fAd; }
    UChar_t      GetADC()      const { return fAdc; }
    UShort_t     GetModuleId() const { return fModuleId; }
    Int_t        GetModuleRorcEquipId()   const { return fEquipId; }
    Int_t        GetModuleRorcEquipType() const { return fEquipType; }
    Int_t        GetNumberOfStrips() const { return fNumberOfStrips; }
    Long_t       GetEventsNumber()   const { return fEventsNumber; }
//    AliITSChannelDaSSD*  GetStrip()  const { return fStrips; }
    AliITSChannelDaSSD*  GetStrip(const Int_t stripnumber)  const 
                                { return (fStrips) ? fStrips[stripnumber] : NULL; }

    AliITSNoiseSSD* GetCalibrationSSDModule() const;
    
    Bool_t  SetEventsNumber(const Long_t eventsnumber);
    Bool_t  SetModuleIdData (const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const UShort_t moduleID);
    void    SetModuleFEEId (const UChar_t ddlID, const UChar_t ad, const UChar_t adc);
    void    SetModuleRorcId (const Int_t equipid, const Int_t equiptype);
    void    SetModuleId (const UShort_t moduleID) { fModuleId = moduleID; }
    void    SetStrip(AliITSChannelDaSSD* strip, const Int_t strID) { if ((fStrips) && (strID <= fNumberOfStrips)) fStrips[strID] = strip; }
    void    DeleteSignal() {if (fStrips) for (Int_t i = 0; i < fNumberOfStrips; i++) 
                                            if (fStrips[i]) fStrips[i]->DeleteSignal(); fEventsNumber = 0; }
    static Int_t GetStripsPerModuleConst() { return  fgkStripsPerModule;  }
    static Int_t GetPNStripsPerModule()    { return  fgkPNStripsPerModule;}
    static Int_t GetStripsPerChip()        { return  fgkStripsPerChip;    }

  protected :
    static const Int_t   fgkStripsPerModule   = 1536;
    static const Int_t   fgkPNStripsPerModule = 768;
    static const Int_t   fgkStripsPerChip     = 128;
    static const UChar_t fgkMaxAdNumber       = 9;
    static const UChar_t fgkMaxAdcNumber      = 13;

    Int_t          fEquipId;        // required to access to rorc
    Int_t          fEquipType;
    UChar_t        fDdlId;          // index of DDL, ITS SSD: 33-48
    UChar_t        fAd;             // index of AD module     0-9
    UChar_t        fAdc;            // index of ADC module    0-5, 8-13
    UShort_t       fModuleId;       // Module number          0-1697
    
    Int_t                 fNumberOfStrips;
    AliITSChannelDaSSD  **fStrips;         //[fNumberOfStrips]  Array of *AliITSChannelDaSSD

    Long_t            fEventsNumber;       // number of events for fsignal memory allocation

  private:
    Bool_t ForbiddenAdcNumber (const UChar_t adcn) const { return ((adcn == 6) || (adcn == 7)); }
 
    ClassDef(AliITSModuleDaSSD, 1) 
 
};

#endif
