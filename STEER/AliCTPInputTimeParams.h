#ifndef ALITIMESTAMP_H
#define ALITIMESTAMP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliCTPInputTimeParams.h 22322 2007-11-22 11:43:14Z cvetan $ */

///////////////////////////////////////////////////////////////////////////////
//


class TObject;

class AliCTPInputTimeParams : public TObject {

public:
                         AliCTPInputTimeParams();
                         AliCTPInputTimeParams( TString& name, UInt_t& level, UInt_t delay, TString edge );   

              virtual   ~AliCTPInputTimeParams() {}
                         AliCTPInputTimeParams( const AliCTPInputTimeParams &ctptime );
         AliCTPInputTimeParams&   operator=(const AliCTPInputTimeParams& ctptime);
              
      // Getters
              TString    GetInputName()      const { return fName; }
               UInt_t    GetLevel()     const { return fLevel; }       
               UInt_t    GetDelay() const { return fDelay; }  
              TString    GetEdge()    const { return fEdge; }      

     // Setters
                 void    SetCTPInputTimeParams( TString name, UInt_t level, 
                                       UInt_t delay, TString edge );
              
         virtual void    Print( const Option_t* opt ="" ) const;
                               
              
protected:
	     TString	fName;  	
              UInt_t    fLevel; 
              UInt_t    fDelay;
	     TString	fEdge;
                         
private:                         

   ClassDef( AliCTPInputTimeParams, 1 )  
};                                                                         


#endif
