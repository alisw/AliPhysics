#ifndef ALIITSHANDLEDASSD_H
#define ALIITSHANDLEDASSD_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*                                                                        */
/* $Id$ */

#include "TObject.h"
#include "AliITSModuleDaSSD.h"

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides ITS SSD data handling
/// used by DA. 
///
///////////////////////////////////////////////////////////////////////////////

class TObjArray;

class AliITSHandleDaSSD : public TObject {
  public :
    AliITSHandleDaSSD();
    explicit AliITSHandleDaSSD(Char_t *rdfname);
    AliITSHandleDaSSD(const AliITSHandleDaSSD& ssdadldc);
    AliITSHandleDaSSD& operator = (const AliITSHandleDaSSD& ssdadldc);
    virtual ~AliITSHandleDaSSD();

    virtual Bool_t Init(Char_t *rdfname, const Char_t *configfname = NULL);
    virtual Bool_t ReadConfigurationFile(const Char_t *configfname = NULL) const;
    Bool_t  SetRawFileName (Char_t *rdfname) {return Init(rdfname); }

    Int_t              GetNumberOfModules() const { return fNumberOfModules; }
    UInt_t             GetLdcId() const { return fLdcId; }
    UInt_t             GetRunId() const { return fRunId; }
    AliITSModuleDaSSD* GetModule (const UChar_t ddlID, const UChar_t ad, const UChar_t adc) const;
    AliITSModuleDaSSD* GetModule (const Int_t index) const 
                               {if ((fModules) && (index < fNumberOfModules)) return fModules[index]; else return NULL;}
    Int_t GetModuleIndex (const UChar_t ddlID, const UChar_t ad, const UChar_t adc) const;
    AliITSChannelDaSSD* GetStrip (const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const UShort_t stripID) const;
    TObjArray*  GetCalibrationSSDLDC()  const;
    Bool_t      SaveCalibrationSSDLDC(Char_t*& dafname) const;
    
    void    SetModIndProcessed(Int_t mi) {fModIndProcessed = mi;}
    void    SetModIndRead (Int_t mr)  {fModIndRead = mr;}
    Bool_t  SetNumberOfModules (const Int_t numberofmodules);
    Bool_t  SetModule(AliITSModuleDaSSD *const module, const Int_t index); 
    Bool_t  ReadCalibrationDataFile (char* fileName, const Long_t eventsnumber);
    Int_t   ReadModuleRawData (const Int_t modulesnumber);  

    virtual Bool_t  CalculatePedestal(AliITSModuleDaSSD *const module);
    virtual Bool_t  CalculateNoise(AliITSModuleDaSSD *const module, const Bool_t CorrectCM = kFALSE);
    virtual Bool_t  CalculateNoiseCM(AliITSModuleDaSSD *const module);
    virtual Bool_t  CalculateCM(AliITSModuleDaSSD *const module);
    virtual Bool_t  ProcessRawData(const Int_t nmread = fgkNumberOfSSDModulesPerDdl);
    virtual Bool_t  RelocateModules();
    virtual Bool_t  AllocateSimulatedModules(const Int_t copymodind = 0);
     
    virtual void    Reset();
    virtual Short_t RetrieveModuleId(const UChar_t ddlID, const UChar_t ad, const UChar_t adc) const;
    Bool_t  DumpModInfo(const Float_t meannosethreshold) const;
    Bool_t  PrintModCalibrationData(const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const Char_t *fname = NULL) const;
    void    DumpInitData(const Char_t *str = "") const;
    void    DeleteSignalAll() { if (fModules) for (Int_t i = 0; i < fNumberOfModules; i++) if (fModules[i]) fModules[i]->DeleteSignal();}
    void    DeleteSignal() { if (fModules) for (Int_t i = fModIndProcessed; i < fModIndRead; i++) if (fModules[i]) fModules[i]->DeleteSignal();}

    static Int_t GetNumberOfSSDModulesConst() { return fgkNumberOfSSDModules; }

  protected :
    static const Int_t    fgkNumberOfSSDModules ;        // Number of SSD modules in ITS
    static const Int_t    fgkNumberOfSSDModulesPerDdl;   // Number of SSD modules in ITS
    static const Float_t  fgkPedestalThresholdFactor;    // Defalt value for fPedestalThresholdFactor 
    static const Float_t  fgkCmThresholdFactor;          // Defalt value for fCmThresholdFactor 

    Char_t              *fRawDataFileName;       // Name of the file with raw data
    Int_t                fNumberOfModules;       // number of AliITSModuleDaSSD to allocate
    AliITSModuleDaSSD  **fModules;               //[fNumberOfModules] array of pointer on AliITSModuleDaSSD objects (1698 SSD  Modules)
    Int_t                fModIndProcessed;       //! index of the last module in fModules array with processed data
    Int_t                fModIndRead;            //! index of the last module in fModules array with adc data present (read)
    Long_t               fNumberOfEvents;        // Number of physics or calibration events in raw data file fRawDataFileName
    
    UInt_t               fLdcId;                 //  LDC number, read from header
    UInt_t               fRunId;                 //  Run number, read from header

    Float_t     fPedestalThresholdFactor;        // configuration parameter: ThresholdFactor for pedestal calculation 
    Float_t     fCmThresholdFactor;              // configuration parameter: ThresholdFactor for CM calculation 

  private :
    Bool_t   SignalOutOfRange (const Short_t signal) const { return (signal >= AliITSChannelDaSSD::GetOverflowConst()); }

    ClassDef(AliITSHandleDaSSD, 2)

};

#endif
