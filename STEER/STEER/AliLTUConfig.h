#ifndef ALILTUCONFIG_H
#define ALILTUCONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Class that describes a detector LTU configuration                        //
//                                                                           //
//                                                                           //
//  Presently we store a subset of the LTU parameters:                       //
//  FineDelay1 3126 # picosec                                                //
//  FineDelay2 20459 # picosec                                               //
//  BC_DELAY_ADD 18 # ns                                                     //
//                                                                           //
//  cvetan.cheshkov@cern.ch 3/9/2010                                         //
///////////////////////////////////////////////////////////////////////////////


#include <TObject.h>

class AliLTUConfig : public TObject {

public:
                          AliLTUConfig(): TObject(),
                            fDetectorId(-1),
			    fFineDelay1(0.),
			    fFineDelay2(0.),
			    fBCDelayAdd(0.)
			    {}
                          AliLTUConfig(UChar_t detectorId, Float_t fineDelay1, Float_t fineDelay2, Float_t bcDelayAdd): TObject(),
                            fDetectorId(detectorId),
			    fFineDelay1(fineDelay1),
			    fFineDelay2(fineDelay2),
			    fBCDelayAdd(bcDelayAdd)
			    {}
                          AliLTUConfig(AliLTUConfig & ltu): TObject(ltu),
                            fDetectorId(ltu.fDetectorId),
			    fFineDelay1(ltu.fFineDelay1),
			    fFineDelay2(ltu.fFineDelay2),
			    fBCDelayAdd(ltu.fBCDelayAdd)
			    {}
               virtual   ~AliLTUConfig() {}

  //  Getters
	        Char_t    GetDetectorId() const { return fDetectorId; }
	    const char*   GetDetectorName() const;
	       Float_t    GetFineDelay1() const { return fFineDelay1; }
	       Float_t    GetFineDelay2() const { return fFineDelay2; }
	       Float_t    GetBCDelayAdd() const { return fBCDelayAdd; }

           virtual void   Print( const Option_t* opt ="" ) const;

private:
	   AliLTUConfig & operator=(const AliLTUConfig & );

                Char_t    fDetectorId; // Detector ID, see AliDAQ class for more details
               Float_t    fFineDelay1; // Fine delay in ns
	       Float_t    fFineDelay2; // Fine delay in ns
	       Float_t    fBCDelayAdd; // BC_DELAY_ADD in ns 

   ClassDef( AliLTUConfig, 1 )  // LTU Configuration class
};


#endif
