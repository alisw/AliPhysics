#ifndef ALIITSHANDLEDASSD_H
#define ALIITSHANDLEDASSD_H

#include <math.h>
#include <sstream>
#include <string>
#include "TObject.h"
#include "TFile.h"
#include "TObjArray.h"
#include "AliRawReaderDate.h"

#include "AliITSChannelDaSSD.h"
#include "AliITSModuleDaSSD.h"

#ifndef PHYSICS_EVENT
#define PHYSICS_EVENT 7
#endif

class AliITSHandleDaSSD : public TObject {
  public :
    AliITSHandleDaSSD();
    explicit AliITSHandleDaSSD(const Int_t numberofmodules);
    AliITSHandleDaSSD(const AliITSHandleDaSSD& ssdadldc);
    AliITSHandleDaSSD& operator = (const AliITSHandleDaSSD& ssdadldc);
    virtual ~AliITSHandleDaSSD();

    Int_t        GetNumberOfModules() const { return fNumberOfModules; }

    AliITSModuleDaSSD* GetModule (const UChar_t ddlID, const UChar_t ad, const UChar_t adc) const;
    AliITSModuleDaSSD* GetModule (const Int_t index) const 
                               {if ((fModules) && (index < fNumberOfModules)) return fModules[index]; else return NULL;}
    AliITSChannelDaSSD* GetStrip (const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const UShort_t stripID) const;
    TObjArray*  GetCalibrationSSDLDC()  const;
    Bool_t      SaveCalibrationSSDLDC(std::string& dafname) const;
    
    Bool_t  SetNumberOfModules (const Int_t numberofmodules);
    Bool_t  SetModule(AliITSModuleDaSSD *const module, const Int_t index); 
    Bool_t  ReadCalibrationDataFile (const char* fileName, const Long_t eventsnumber);

    virtual Bool_t  CalculatePedestal();
    virtual Bool_t  CalculateNoise();
    virtual Bool_t  CalculateNoiseCM();

    void    DeleteSignal() { if (fModules) for (Int_t i = 0; i < fNumberOfModules; i++) if (fModules[i]) fModules[i]->DeleteSignal();}

    static Int_t GetNumberOfSSDModulesConst() { return fgkNumberOfSSDModules; }

  protected :
    static const Int_t fgkNumberOfSSDModules = 1698;

    Int_t                fNumberOfModules;       // number of AliITSModuleDaSSD to allocate
    AliITSModuleDaSSD  **fModules;               //[fNumberOfModules]  array of all SSD 1698 Modules (2608128 strips)
    UInt_t               fLdcId;                 //  LDC number, read from header
    UInt_t               fRunId;                 //  Run number, read from header

  private :
    Bool_t   CalculateCM(const Int_t modind, const Int_t stripind, Float_t* const cm);
    Bool_t   RelocateModules();
    Bool_t   SignalOutOfRange (const Short_t signal) const { return (signal >= AliITSChannelDaSSD::GetOverflowConst()); }

    ClassDef(AliITSHandleDaSSD, 1)

};

#endif
