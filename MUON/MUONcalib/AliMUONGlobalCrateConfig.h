/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 

/// \ingroup calib
/// \class AliMUONGlobalCrateConfig
/// \brief The class defines the configuration of global crate
///
/// \author Ch. Finck, Subatech Nantes

#ifndef ALI_MUON_GLOBAL_CRATE_CONFIG_H
#define ALI_MUON_GLOBAL_CRATE_CONFIG_H

#include <TNamed.h>
#include <TString.h>

class AliMUONGlobalCrateConfig : public  TNamed {

  public:
    AliMUONGlobalCrateConfig();
    virtual ~AliMUONGlobalCrateConfig();
    
    // methods
    Int_t ReadData(const TString& fileName = "");

    // global crate enable
          /// set global crate enbale
    void SetGlobalCrateEnable(UInt_t enable) {fGlobalCrateEnable = enable;}
         /// Get global crate enbale
    UInt_t GetGlobalCrateEnable() const {return fGlobalCrateEnable;}
    
    // Jtag
            /// Get Jtag board VME address
    ULong_t GetJtagVmeAddr() const {return fJtagVmeAddr;}
           /// Set Jtag board VME address   
    void    SetJtagVmeAddr(ULong_t addr) {fJtagVmeAddr = addr;}
    
            /// Get Jtag board Clock Divider
    UInt_t  GetJtagClockDiv() const {return fJtagClockDiv;}
            /// Set Jtag board Clock Divider
    void    SetJtagClockDiv(UInt_t clk) {fJtagClockDiv = clk;}
    
            /// Get Jtag board Rx Phase
    UInt_t  GetJtagRxPhase() const {return fJtagRxPhase;}
            /// Set Jtag board Rx Phase
    void    SetJtagRxPhase(UInt_t rx) {fJtagRxPhase = rx;}

            /// Get Jtag board Read out Delay
    UInt_t  GetJtagRdDelay() const {return fJtagRdDelay;}
            /// Set Jtag board Read out Delay
    void    SetJtagRdDelay(UInt_t rd) {fJtagRdDelay = rd;}
    
            /// Get Jtag enabled lines
    Bool_t GetEnableJtag(Int_t index) const;
           /// Set Jtag enable word
    void   SetEnableJtag(UChar_t en) {fEnableJtag = en;} 
           /// Get Jtag enable word
    UChar_t GetEnableJtag() const {return fEnableJtag;}  
           /// Set First Darc enable word
    void   SetEnableFirstDarc(UChar_t en) {fEnableFirstDarc = en;} 
           /// Get First Darc enable word
    UChar_t GetEnableFirstDarc() const {return fEnableFirstDarc;}  
           /// Get First Darc enable lines
    Bool_t GetEnableFirstDarc(Int_t index) const;  
           /// Set Second Darc enable word
    void   SetEnableSecondDarc(UChar_t en) {fEnableSecondDarc = en;} 
           /// Get Second Darc enable word
    UChar_t GetEnableSecondDarc() const {return fEnableSecondDarc;}  
           /// Get Second Darc enable lines
    Bool_t GetEnableSecondDarc(Int_t index) const;  
    
           /// Get Jtag Crate names
    TString GetJtagCrateName(Int_t jtagLine, Int_t index) const;
           /// Set Jtag Crate names
    void    SetJtagCrateName(Int_t index, TString name);
    
    // first Darc Board
            /// Get First Darc board VME address
    ULong_t GetFirstDarcVmeAddr() const        {return fFirstDarcVmeAddr;}
            /// Get First Darc board VME address
    void    SetFirstDarcVmeAddr(ULong_t addr)  {fFirstDarcVmeAddr = addr;}
    
            /// Get type for First Darc board   
    Int_t   GetFirstDarcType() const           {return fFirstDarcType;}
            /// Get type for First Darc board   
    void    SetFirstDarcType(Int_t type)       {fFirstDarcType = type;}
    
            /// Get disable word for First Darc board   
    UChar_t GetFirstDarcDisable() const         {return fFirstDarcDisable;} 
            /// Get disable per regional crate for First Darc board   
    Bool_t  GetFirstDarcDisable(Int_t iCrate) const   {return !((fFirstDarcDisable >> iCrate) & 0x1);} 
            /// Set disable word for First Darc board   
    void    SetFirstDarcDisable(UChar_t en)     {fFirstDarcDisable = en;} 
    
            /// Get L0 Delay for First Darc board   
    UInt_t  GetFirstDarcL0Delay() const        {return fFirstDarcL0Delay;}
            /// Set L0 Delay for First Darc board   
    void    SetFirstDarcL0Delay(UInt_t delay)  {fFirstDarcL0Delay = delay;}
    
            /// Get L1 Time Out for First Darc board   
    UInt_t  GetFirstDarcL1TimeOut() const      {return fFirstDarcL1TimeOut;}
            /// Set L1 Time Out for First Darc board   
    void    SetFirstDarcL1TimeOut(UInt_t time) {fFirstDarcL1TimeOut = time;}
 
            /// Get global L0  delay for First Darc board   
    UInt_t  GetFirstDarcGlobalL0() const      {return fFirstDarcGlobalL0;}
           /// set global L0  delay for First Darc board   
    void    SetFirstDarcGlobalL0(UInt_t time) {fFirstDarcGlobalL0 = time;}
           
            /// Get configuration for First Darc board   
    UInt_t  GetFirstDarcConfig() const      {return fFirstDarcConfig;}
           /// set configuration for First Darc board   
    void    SetFirstDarcConfig(UInt_t conf) {fFirstDarcConfig = conf;}   
    
           /// Get First Darc Crate names
    TString GetFirstDarcCrateName(Int_t index) const;
           /// Set First Darc Crate names
    void    SetFirstDarcCrateName(Int_t index, TString name);
    
    // second Darc Board
            /// Get Second Darc board VME address
    ULong_t GetSecondDarcVmeAddr() const        {return fSecondDarcVmeAddr;}
            /// Set Second Darc board VME address
    void    SetSecondDarcVmeAddr(ULong_t addr)  {fSecondDarcVmeAddr = addr;}
    
            /// Get type for Second Darc board
    Int_t   GetSecondDarcType() const           {return fSecondDarcType;}
            /// Set type for Second Darc board   
    void    SetSecondDarcType(Int_t type)       {fSecondDarcType = type;}
    
            /// Get disable word for Second Darc board   
    UChar_t GetSecondDarcDisable() const         {return fSecondDarcDisable;} 
            /// Get disable per regional crate for Second Darc board   
    Bool_t  GetSecondDarcDisable(Int_t iCrate) const  {return !((fSecondDarcDisable >> iCrate) & 0x1);} 
            /// Set disable word for Second Darc board   
    void    SetSecondDarcDisable(UChar_t en)     {fSecondDarcDisable = en;} 
    
            /// Get L0 Delay for Second Darc board   
    UInt_t  GetSecondDarcL0Delay() const        {return fSecondDarcL0Delay;}
            /// Set L0 Delay for Second Darc board   
    void    SetSecondDarcL0Delay(UInt_t delay)  {fSecondDarcL0Delay = delay;}
            /// Get L1 Time Out for Second Darc board   
    UInt_t  GetSecondDarcL1TimeOut() const      {return fSecondDarcL1TimeOut;}
            /// Set L1 Time Out for Second Darc board   
    void    SetSecondDarcL1TimeOut(UInt_t time) {fSecondDarcL1TimeOut = time;}
    
            /// Get global L0  delay for Second Darc board   
    UInt_t  GetSecondDarcGlobalL0() const      {return fSecondDarcGlobalL0;}
           /// set global L0  delay for Second Darc board   
    void    SetSecondDarcGlobalL0(UInt_t time) {fSecondDarcGlobalL0 = time;}
      
           /// Get configuration for Second Darc board   
    UInt_t  GetSecondDarcConfig() const      {return fSecondDarcConfig;}
           /// set configuration for Second Darc board   
    void    SetSecondDarcConfig(UInt_t conf) {fSecondDarcConfig = conf;}
    
           /// Get Second Darc Crate names
    TString GetSecondDarcCrateName(Int_t index) const;
           /// Set Second Darc Crate names
    void    SetSecondDarcCrateName(Int_t index, TString name);
        
    // global board
            /// Get Global board VME address
    ULong_t GetGlobalVmeAddr() const        {return fGlobalVmeAddr;}
            /// Set Global board VME address
    void    SetGlobalVmeAddr(ULong_t addr)  {fGlobalVmeAddr = addr;}
    
            /// Get register for Global
    UInt_t GetGlobalRegister(Int_t index) const;
           /// Set register for Global
    void   SetGlobalRegister(Int_t index, UInt_t reg);
           /// Get register word for Global
    UInt_t* GetGlobalRegister() {return fGlobalRegisters;}
           /// Set mask for the global input
      void   SetGlobalMask(Int_t index, UInt_t mask);
           /// Get mask for the global input
    UInt_t GetGlobalMask(Int_t index) const;
           /// Indicates if global masks are active on global inputs
    Bool_t GetMasksOn() const;    

    // fet board
            /// Get FET board VME address
    ULong_t GetFetVmeAddr()  const       {return fFetVmeAddr;}
            /// Set FET board VME address
    void    SetFetVmeAddr(ULong_t addr)  {fFetVmeAddr = addr;}
    
            /// Get register for FET
    UInt_t GetFetRegister(Int_t index) const;
            /// Set register for FET
    void   SetFetRegister(Int_t index, UInt_t reg);
            /// Set register word for FET
    UInt_t* GetFetRegister() {return fFetRegisters;}
    
    //static members
            /// Get Jtag Name identifier
    const Char_t* GetJtagName()       const  {return fgkJtagName;}
            /// Get First Darc Name identifier
    const Char_t* GetFirstDarcName()  const  {return fgkFirstDarcName;}
            /// Get Second Darc Name identifier
    const Char_t* GetSecondDarcName() const  {return fgkSecondDarcName;}
            /// Get Global Name identifier   
    const Char_t* GetGlobalName()     const  {return fgkGlobalName;}
            /// Get Global Name identifier   
    const Char_t* GetFetName()        const  {return fgkFetName;}
    
            /// Get number of registers for Global
    Int_t   GetGlobalNofRegisters() const {return fgkGlobalNofRegisters;}
            /// Get number of registers for FET
    Int_t   GetFetNofRegisters()    const {return fgkFetNofRegisters;}
            /// Get number of JTag lines
    Int_t   GetJtagNofLines()       const {return fgkJtagNofLines;}
            /// Get number of Darc Crate lines
    Int_t   GetDarcNofLines()       const {return fgkDarcNofLines;}
    
  private:
    /// Not implemented
    AliMUONGlobalCrateConfig(const AliMUONGlobalCrateConfig& rhs);
    /// Not implemented
    AliMUONGlobalCrateConfig& operator=(const AliMUONGlobalCrateConfig& rhs);

    // data members
    UInt_t       fGlobalCrateEnable;   ///< Global Crate Enable
    ULong_t      fJtagVmeAddr;         ///< JTag VME address
    UInt_t       fJtagClockDiv;        ///< Clock Divider number for JTag
    UInt_t       fJtagRxPhase;         ///< Rx phase number for JTag 
    UInt_t       fJtagRdDelay;         ///< Read delay  for JTag 
    UChar_t      fEnableJtag;          ///< Enable mask for JTag lines
    TString      fJtagCrateName[16];   ///< Crate name for the Jtag lines
    TString      fFirstDarcCrateName[8];   ///< Crate name for the First Darc lines
    TString      fSecondDarcCrateName[8];   ///< Crate name for the Second Darc lines
                                       
    ULong_t      fFirstDarcVmeAddr;    ///< First Darc Board VME Address
    Int_t        fFirstDarcType;       ///< Type of the first Darc Board                             
    UChar_t      fFirstDarcDisable;     ///< disable the readout of the 8 crates connected to this board
    UInt_t       fFirstDarcL0Delay;    ///< L0 delay for this board
    UInt_t       fFirstDarcL1TimeOut;  ///< L1 time out for this board
    UInt_t       fFirstDarcGlobalL0 ;  ///< L0 global l0 delay this board
    UInt_t       fFirstDarcConfig ;    ///< Trigger configuration this board
 
    ULong_t      fSecondDarcVmeAddr;   ///< Second Darc Board VME Address 
    Int_t        fSecondDarcType;      ///< Type of the first Darc Board                             
    UChar_t      fSecondDarcDisable;   ///< disable the readout of the 8 crates connected to this board
    UInt_t       fSecondDarcL0Delay;   ///< L0 delay for this board
    UInt_t       fSecondDarcL1TimeOut; ///< L1 time out for this board
    UInt_t       fSecondDarcGlobalL0;  ///<  Global L0 delay for this board
    UInt_t       fSecondDarcConfig ;   ///< Trigger configuration this board
 
    ULong_t      fGlobalVmeAddr;       ///< Global Board VME Address 
    UInt_t       fGlobalRegisters[13]; ///< Global registers

    ULong_t      fFetVmeAddr;          ///< Fet Board VME Address 
    UInt_t       fFetRegisters[7];     ///< Fet registers                                  
                                       
    UChar_t      fEnableFirstDarc;     ///< Enable mask for First Darc lines
    UChar_t      fEnableSecondDarc;    ///< Enable mask for Second Darc lines

    static const Char_t* fgkJtagName;       ///< JTag Board name                         
    static const Char_t* fgkFirstDarcName;  ///< First DARC board name                         
    static const Char_t* fgkSecondDarcName; ///< Second DARC board name
    static const Char_t* fgkGlobalName;     ///< Global Board name
    static const Char_t* fgkFetName;        ///< FET Board name                     
               
    static const Int_t fgkGlobalNofRegisters;  ///< Number of registers for Global Board                      
    static const Int_t fgkFetNofRegisters;     ///< Number of registers for Fet
    static const Int_t fgkJtagNofLines;        ///< Number of lines for Jtag
    static const Int_t fgkDarcNofLines;        ///< Number of lines for Darc Crate

  ClassDef(AliMUONGlobalCrateConfig,4)  
};

#endif 














