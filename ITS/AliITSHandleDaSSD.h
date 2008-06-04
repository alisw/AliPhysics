#ifndef ALIITSHANDLEDASSD_H
#define ALIITSHANDLEDASSD_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*                                                                        */
/* $Id$ */

#include <string>
#include "TObject.h"
#include "AliITSModuleDaSSD.h"

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides ITS SSD data handling
/// used by DA. 
//  Author: Oleksandr Borysov
//  Date: 19/05/2008
///////////////////////////////////////////////////////////////////////////////

using namespace std;

class TObjArray;

class AliITSHandleDaSSD : public TObject {
  public :
    AliITSHandleDaSSD();
    explicit AliITSHandleDaSSD(Char_t *rdfname);
    AliITSHandleDaSSD(const AliITSHandleDaSSD& ssdadldc);
    AliITSHandleDaSSD& operator = (const AliITSHandleDaSSD& ssdadldc);
    virtual ~AliITSHandleDaSSD();

    virtual Bool_t Init(Char_t *rdfname);
    Bool_t  SetRawFileName (Char_t *rdfname) {return Init(rdfname); }

    void    SetZsDefaul(Int_t zs)        { fZsDefault = zs;       }
    void    SetOffsetDefault(Int_t offs) { fOffsetDefault = offs; }
    void    SetZsFactor(Float_t zsf)     { fZsFactor = zsf;       }
    void    SetPedestalThresholdFactor(Float_t pthf) { fPedestalThresholdFactor = pthf; }
    void    SetCmThresholdFactor(Float_t cmthf)      { fCmThresholdFactor = cmthf;      }

    Int_t   GetZsDefault() const     { return fZsDefault;     }
    Int_t   GetOffsetDefault() const { return fOffsetDefault; }
    Float_t GetZsFactor() const      { return fZsFactor;      }
    Float_t GetPedestalThresholdFactor() const { return fPedestalThresholdFactor; }
    Float_t GetCmThresholdFactor() const       { return fCmThresholdFactor;       }
    
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
    virtual Bool_t  ReadStaticBadChannelsMap(const Char_t *filename = NULL);  
    virtual Bool_t  ReadDDLModuleMap(const Char_t *filename = NULL);  
    Int_t   ReadCalibrationDataFile (char* fileName, const Long_t eventsnumber);
    Int_t   ReadModuleRawData (const Int_t modulesnumber);  

    virtual Bool_t  CalculatePedestal(AliITSModuleDaSSD *const module);
    virtual Bool_t  CalculateNoise(AliITSModuleDaSSD *const module);
    virtual Bool_t  CalculateNoiseCM(AliITSModuleDaSSD *const module);
    virtual Bool_t  CalculateCM(AliITSModuleDaSSD *const module);
    virtual Bool_t  AddFeromCm(AliITSModuleDaSSD *const module);
    virtual Bool_t  ProcessRawData(const Int_t nmread = fgkNumberOfSSDModulesPerDdl);
    virtual Bool_t  RelocateModules();
    virtual Bool_t  AllocateSimulatedModules(const Int_t copymodind = 0);

    Bool_t  AdDataPresent(const Int_t ddl, const Int_t ad) const;
    Int_t   DdlToEquipmentId (Int_t ddl) const { return (512 + ddl); }
    Int_t   ChannelIsBad(const UChar_t ddl, const UChar_t ad, const UChar_t adc, const Int_t strn) const;
    Bool_t  SaveEqSlotCalibrationData(const Int_t ddl, const Int_t ad, const Char_t *fname) const;
    ULong_t OffsetValue(const AliITSChannelDaSSD *strip, const UChar_t ddl = 0, const UChar_t ad = 0, 
                                 const UChar_t adc = 0, const Int_t strn = -1) const;
    ULong_t OffsetValue(const UChar_t ddl, const UChar_t ad, const UChar_t adc, const Int_t strn) const;
    ULong_t ZsThreshold(AliITSChannelDaSSD *strip) const;
    ULong_t ZsThreshold(const UChar_t ddl, const UChar_t ad, const UChar_t adc, const Int_t strn) const;
    
    virtual void    Reset();
    virtual Short_t RetrieveModuleId(const UChar_t ddlID, const UChar_t ad, const UChar_t adc) const;
    Bool_t  DumpModInfo(const Float_t meannosethreshold) const;
    Bool_t  PrintModCalibrationData(const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const Char_t *fname = NULL) const;
    void    DumpInitData(const Char_t *str = "") const;
    void    DeleteSignalAll() { if (fModules) for (Int_t i = 0; i < fNumberOfModules; i++) if (fModules[i]) fModules[i]->DeleteSignal();}
    void    DeleteSignal() { if (fModules) for (Int_t i = fModIndProcessed; i < fModIndRead; i++) if (fModules[i]) fModules[i]->DeleteSignal();}
    void    DeleteCMAll() { if (fModules) for (Int_t i = 0; i < fNumberOfModules; i++) if (fModules[i]) fModules[i]->DeleteCM();}
    void    DeleteCM() { if (fModules) for (Int_t i = fModIndProcessed; i < fModIndRead; i++) if (fModules[i]) fModules[i]->DeleteCM();}
    void    DeleteCMFerom() { if (fModules) for (Int_t i = fModIndProcessed; i < fModIndRead; i++) if (fModules[i]) fModules[i]->DeleteCMFerom ();}

    static Int_t GetNumberOfSSDModulesConst() { return fgkNumberOfSSDModules; }

  protected :

    static const Int_t    fgkNumberOfSSDModules ;        // Number of SSD modules in ITS
    static const Int_t    fgkNumberOfSSDModulesPerDdl;   // Number of SSD modules in DDL
    static const Int_t    fgkNumberOfSSDModulesPerSlot;  // Number of SSD modules in Slot
    static const Int_t    fgkNumberOfSSDDDLs;            // Number of DDLs in SSD
    static const Int_t    fgkNumberOfSSDSlotsPerDDL;     // Number of SSD slots per DDL
    static const Float_t  fgkPedestalThresholdFactor;    // Defalt value for fPedestalThresholdFactor 
    static const Float_t  fgkCmThresholdFactor;          // Defalt value for fCmThresholdFactor 
   
    static const UInt_t   fgkZsBitMask ;           // Bit mask for FEROM ZS
    static const UInt_t   fgkOffSetBitMask;        // Bit mask for FEROM Offset correction
    static const UInt_t   fgkBadChannelMask;       // Mask to suppress the channel from the bad channel list
    static const Int_t    fgkAdcPerDBlock;         // FEROM configuration file constant
     
    Char_t              *fRawDataFileName;       // Name of the file with raw data
    Int_t                fNumberOfModules;       // number of AliITSModuleDaSSD to allocate
    AliITSModuleDaSSD  **fModules;               //[fNumberOfModules] array of pointer on AliITSModuleDaSSD objects (1698 SSD  Modules)
    Int_t                fModIndProcessed;       //! index of the last module in fModules array with processed data
    Int_t                fModIndRead;            //! index of the last module in fModules array with adc data present (read)
    Int_t               *fModIndex;              //! index array for fModules
    Long_t               fNumberOfEvents;        // Number of physics or calibration events in raw data file fRawDataFileName

    TObjArray           *fStaticBadChannelsMap;  // Static bad channels map read from the file
    Int_t               *fDDLModuleMap;          //! DDL map  
    
    UInt_t               fLdcId;                 //  LDC number, read from header
    UInt_t               fRunId;                 //  Run number, read from header

    Float_t         fPedestalThresholdFactor;        // configuration parameter: ThresholdFactor for pedestal calculation 
    Float_t         fCmThresholdFactor;              // configuration parameter: ThresholdFactor for CM calculation 
    Int_t           fZsDefault;                      // default value for ZS threshold
    Int_t           fOffsetDefault;                  // default value for offset correction
    Float_t         fZsFactor;                       // zs factor 3.0

  protected :
    Bool_t   SignalOutOfRange (const Short_t signal) const { return (signal >= AliITSChannelDaSSD::GetOverflowConst()); }
    string   ConvBase(const unsigned long value, const long base) const;

    ClassDef(AliITSHandleDaSSD, 5)

};

#endif

