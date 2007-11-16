/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: $ 

/// \ingroup management
/// \class AliMpGlobalCrate
/// \brief The class defines the properties of global crate
///
/// \author Ch. Finck, Subatech Nantes

#ifndef ALI_MP_GLOBAL_CRATE_H
#define ALI_MP_GLOBAL_CRATE_H

#include <TNamed.h>
#include <TString.h>

class AliMpGlobalCrate : public  TNamed {

  public:
    AliMpGlobalCrate();
    AliMpGlobalCrate(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpGlobalCrate();
    

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
    const Int_t   GetGlobalNofRegisters() const {return fgkGlobalNofRegisters;}
            /// Get number of registers for FET
    const Int_t   GetFetNofRegisters()    const {return fgkFetNofRegisters;}
            /// Get number of JTag lines
    const Int_t   GetJtagNofLines()       const {return fgkJtagNofLines;}
    
  private:
    /// Not implemented
    AliMpGlobalCrate(const AliMpGlobalCrate& rhs);
    /// Not implemented
    AliMpGlobalCrate& operator=(const AliMpGlobalCrate& rhs);

    // data members
    ULong_t      fJtagVmeAddr;         ///< JTag VME address
    UInt_t       fJtagClockDiv;        ///< Clock Divider number for JTag
    UInt_t       fJtagRxPhase;         ///< Rx phase number for JTag 
    UInt_t       fJtagRdDelay;         ///< Read delay  for JTag 
    UChar_t      fEnableJtag;          ///< Enable mask for JTag lines
    TString      fJtagCrateName[16];   ///< Crate name for the Jtag lines
                                       
    ULong_t      fFirstDarcVmeAddr;    ///< First Darc Board VME Address
    Int_t        fFirstDarcType;       ///< Type of the first Darc Board                             
    UChar_t      fFirstDarcDisable;     ///< disable the readout of the 8 crates connected to this board
    UInt_t       fFirstDarcL0Delay;    ///< L0 delay for this board
    UInt_t       fFirstDarcL1TimeOut;  ///< L1 time out for this board

    ULong_t      fSecondDarcVmeAddr;   ///< Second Darc Board VME Address 
    Int_t        fSecondDarcType;      ///< Type of the first Darc Board                             
    UChar_t      fSecondDarcDisable;   ///< disable the readout of the 8 crates connected to this board
    UInt_t       fSecondDarcL0Delay;   ///< L0 delay for this board
    UInt_t       fSecondDarcL1TimeOut; ///< L1 time out for this board
                                       
    ULong_t      fGlobalVmeAddr;       ///< Global Board VME Address 
    UInt_t       fGlobalRegisters[13]; ///< Global registers

    ULong_t      fFetVmeAddr;          ///< Fet Board VME Address 
    UInt_t       fFetRegisters[7];     ///< Fet registers                                  
                                       
    static const Char_t* fgkJtagName;       ///< JTag Board name                         
    static const Char_t* fgkFirstDarcName;  ///< First DARC board name                         
    static const Char_t* fgkSecondDarcName; ///< Second DARC board name
    static const Char_t* fgkGlobalName;     ///< Global Board name
    static const Char_t* fgkFetName;        ///< FET Board name                     
               
    static const Int_t fgkGlobalNofRegisters;  ///< Number of registers for Global Board                      
    static const Int_t fgkFetNofRegisters;     ///< Number of registers for Fet
    static const Int_t fgkJtagNofLines;        ///< Number of lines for Jtag

  ClassDef(AliMpGlobalCrate,1)  // The class collectiong electronics properties of DDL
};

#endif //ALI_MP_GLOBAL__CRATE_H














