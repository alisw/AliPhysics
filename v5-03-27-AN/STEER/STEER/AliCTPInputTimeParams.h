#ifndef ALICTPINPUTTIMEPARAMS_H
#define ALICTPINPUTTIMEPARAMS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

class TObject;

class AliCTPInputTimeParams : public TObject {

public:
                         AliCTPInputTimeParams();
                         AliCTPInputTimeParams( TString& name, UInt_t& level, UInt_t delay, TString edge, UInt_t deltamin, UInt_t deltamax );   

              virtual   ~AliCTPInputTimeParams() {}
                         AliCTPInputTimeParams( const AliCTPInputTimeParams &ctptime );
         AliCTPInputTimeParams&   operator=(const AliCTPInputTimeParams& ctptime);
              
      // Getters
              TString    GetInputName()      const { return fName; }
               UInt_t    GetLevel()     const { return fLevel; }       
               UInt_t    GetDelay() const { return fDelay; }  
              TString    GetEdge()    const { return fEdge; }      
               Int_t    GetDeltaMin() const { return fDeltaMin; }  
               Int_t    GetDeltaMax() const { return fDeltaMax; }  
     // Setters
                 void    SetCTPInputTimeParams( TString name, UInt_t level, 
                            UInt_t delay, TString edge, UInt_t deltamin, UInt_t deltamax );
              
         virtual void    Print( const Option_t* opt ="" ) const;
                               
              
protected:
	     TString	fName;  	
              UInt_t    fLevel; 
              UInt_t    fDelay;
	     TString	fEdge;
	      Int_t	fDeltaMin;
	      Int_t	fDeltaMax;
                         
private:                         

   ClassDef( AliCTPInputTimeParams, 3 )  
};                                                                         


#endif
