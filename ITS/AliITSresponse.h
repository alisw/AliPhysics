#ifndef ALIITSRESPONSE_H
#define ALIITSRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

#include "AliITSsegmentation.h"

class TF1;
class TString;
class AliITSgeom;

//----------------------------------------------
//
// ITS response virtual base class
//
class AliITSresponse : public TObject {
 public:
    // Default Constructor
    AliITSresponse();
    // Standard Constructor
    AliITSresponse(Double_t Thickness);
    // Destructor.
    virtual ~AliITSresponse() {}
    //
    // Configuration methods
    //

    // Set Electronics
    virtual void    SetElectronics(Int_t) {}
    // Get Electronics
    virtual Int_t Electronics() const {return 0;}

    // Set maximum Adc-count value
    virtual void    SetMaxAdc(Float_t) {}
    // Get maximum Adc-count value
    virtual Float_t MaxAdc() const {return 0.;}

    // Set maximum Adc-top value
    virtual void    SetDynamicRange(Float_t) {}
    // Get maximum Adc-top value
    virtual Float_t DynamicRange() const {return 0.0;}

    // Set Charge Loss Linear Coefficient
    virtual void    SetChargeLoss(Float_t) {}
    // Get Charge Loss Linear Coefficient
    virtual Float_t ChargeLoss() const {return 0.0;}

    // Set GeVcharge value
    virtual void SetGeVToCharge(Float_t gc=2.778E+8){fGeVcharge = gc;}
    // Returns the value fGeVcharge
    virtual Float_t GetGeVToCharge() const {return fGeVcharge;}
    // Converts deposited energy to number of electrons liberated
    virtual Float_t GeVToCharge(Float_t gev) const {return gev*fGeVcharge;}

    // Diffusion coefficient
    virtual void    SetDiffCoeff(Float_t, Float_t) {}
    // Get diffusion coefficients
    virtual void    DiffCoeff(Float_t &,Float_t &) {}

    // Temperature in [degree K]
    virtual void    SetTemperature(Float_t t=300.0) {fT = t;}
    // Get temperature [degree K]
    virtual Float_t Temperature() const {return fT;}
    // Type of data - real or simulated
    virtual void    SetDataType(const char *) {}
    // Set the impurity concentrations in [#/cm^3]
    virtual void SetImpurity(Double_t n=0.0){fN = n;}
    // Returns the impurity consentration in [#/cm^3]
    virtual Double_t Impurity() const {return fN;}
    // Sets the applied ratio distance/voltage [cm/volt]
    virtual void SetDistanceOverVoltage(Double_t d,Double_t v){fdv = d/v;}
    // Sets the applied ration distance/voltage [cm/volt]. Default value
    // is 300E-4cm/80 volts = 0.000375 cm/volts
    virtual void SetDistanceOverVoltage(Double_t dv=0.000375){fdv = dv;}
    // Returns the ration distance/voltage
    virtual Double_t DistanceOverVoltage() const {return fdv;}
    // Get data type
    virtual const char  *DataType() const {return "";}
 
    // Set parameters options: "same" or read from "file" or "SetInvalid" or...
    virtual void   SetParamOptions(const char*,const char*) {}
    // Set noise parameters 
    virtual void   SetNoiseParam(Float_t, Float_t) {}
    // Number of parameters to be set
    virtual  void   SetNDetParam(Int_t) {}
    // Set detector parameters: gain, coupling ...
    virtual  void   SetDetParam(Float_t *) {}

    // Parameters options
    virtual void   ParamOptions(char *,char*) {}
    virtual Int_t  NDetParam() const {return 0;}
    virtual void   GetDetParam(Float_t *) {}
    virtual void   GetNoiseParam(Float_t&, Float_t&) {}

    // Zero-suppression option - could be 1D, 2D or non-ZeroSuppressed
    virtual void   SetZeroSupp(const char*) {}
    // Get zero-suppression option
    virtual const char *ZeroSuppOption() const {return "";}
     // Set thresholds
    virtual void   SetThresholds(Float_t, Float_t) {}
    virtual void   Thresholds(Float_t &, Float_t &) {}
    // Set min val
    virtual void   SetMinVal(Int_t) {};
    virtual Int_t  MinVal() const {return 0;};

    // Set filenames
    virtual void SetFilenames(const char *,const char *,const char *) {}
    // Filenames
    virtual void   Filenames(char*,char*,char*) {}

    virtual Float_t DriftSpeed() const {return 0.;}
    virtual Bool_t  OutputOption() const {return kFALSE;}
    virtual Bool_t  Do10to8() const {return kTRUE;}
    virtual void    GiveCompressParam(Int_t *) {}

    //
    // Detector type response methods
    // Set number of sigmas over which cluster disintegration is performed
    virtual void    SetNSigmaIntegration(Float_t) {}
    // Get number of sigmas over which cluster disintegration is performed
    virtual Float_t NSigmaIntegration() const {return 0.;}
    // Set number of bins for the gaussian lookup table
    virtual void    SetNLookUp(Int_t) {}
    // Get number of bins for the gaussian lookup table
    virtual Int_t GausNLookUp() const {return 0;}
    // Get scaling factor for bin i-th from the gaussian lookup table
    virtual Float_t GausLookUp(Int_t) const {return 0.;}
    // Set sigmas of the charge spread function
    virtual void    SetSigmaSpread(Float_t, Float_t) {}
    // Get sigmas for the charge spread
    virtual void    SigmaSpread(Float_t &,Float_t &) {}

    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t) const {return 0.;}
    // Charge disintegration
    virtual Float_t IntXZ(AliITSsegmentation *) const {return 0.;}
    // Electron mobility in Si. [cm^2/(Volt Sec)]. T in degree K, N in #/cm^3
    virtual Double_t MobilityElectronSiEmp() const ;
    // Hole mobility in Si. [cm^2/(Volt Sec)]  T in degree K, N in #/cm^3
    virtual Double_t MobilityHoleSiEmp() const ;
    // Einstein relation for Diffusion Coefficient of Electrons. [cm^2/sec]
    //  T in degree K, N in #/cm^3
    virtual Double_t DiffusionCoefficientElectron() const ;
    // Einstein relation for Diffusion Coefficient of Holes. [cm^2/sec]
    //  T in [degree K], N in [#/cm^3]
    virtual Double_t DiffusionCoefficientHole() const ;
    // Electron <speed> under an applied electric field E=Volts/cm. [cm/sec]
    // d distance-thickness in [cm], v in [volts], T in [degree K],
    // N in [#/cm^3]
    virtual Double_t SpeedElectron() const ;
    // Holes <speed> under an applied electric field E=Volts/cm. [cm/sec]
    // d distance-thickness in [cm], v in [volts], T in [degree K],
    // N in [#/cm^3]
    virtual Double_t SpeedHole() const ;
    // Returns the Gaussian sigma == <x^2+z^2> [cm^2] due to the defusion of
    // electrons or holes through a distance l [cm] caused by an applied
    // voltage v [volt] through a distance d [cm] in any material at a
    // temperature T [degree K].
    virtual Double_t SigmaDiffusion3D(Double_t) const ;
    // Returns the Gaussian sigma == <x^2 +y^2+z^2> [cm^2] due to the
    // defusion of electrons or holes through a distance l [cm] caused by an
    // applied voltage v [volt] through a distance d [cm] in any material at a
    // temperature T [degree K].
    virtual Double_t SigmaDiffusion2D(Double_t) const ;
    // Returns the Gaussian sigma == <x^2+z^2> [cm^2] due to the defusion of
    // electrons or holes through a distance l [cm] caused by an applied
    // voltage v [volt] through a distance d [cm] in any material at a
    // temperature T [degree K].
    virtual Double_t SigmaDiffusion1D(Double_t) const ;
    // Prints out the content of this class in ASCII format.
    virtual void Print(ostream *os);
    // Reads in the content of this class in the format of Print
    virtual void Read(istream *is);
 private:
    Double_t fdv;  // The parameter d/v where d is the disance over which the
                   // the potential v is applied d/v [cm/volts]
    Double_t fN;   // the impurity consentration of the material in #/cm^3
    Double_t fT;   // The temperature of the Si in Degree K.
    Double_t fGeVcharge; // Energy to ionize (free an electron).

    ClassDef(AliITSresponse,2) // Detector type response virtual base class 
};
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITSresponse &source);
istream& operator>>(istream &os,AliITSresponse &source);
#endif
