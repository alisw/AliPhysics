#ifndef ALIITSRESPONSE_H
#define ALIITSRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TString.h>

class AliITSsegmentation;
class TF1;
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

    // Set GeVcharge value
    virtual void SetGeVToCharge(Double_t gc=2.778E+8){fGeVcharge = gc;}
    // Returns the value fGeVcharge
    virtual Double_t GetGeVToCharge() const {return fGeVcharge;}
    // Converts deposited energy to number of electrons liberated
    virtual Double_t GeVToCharge(Double_t gev) const {return gev*fGeVcharge;}

    // Diffusion coefficient
    virtual void    SetDiffCoeff(Double_t, Double_t) = 0;
    // Get diffusion coefficients
    virtual void    DiffCoeff(Double_t &,Double_t &) const = 0;

    // Temperature in [degree K]
    virtual void    SetTemperature(Double_t t=300.0) {fT = t;}
    // Get temperature [degree K]
    virtual Double_t Temperature() const {return fT;}
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
    virtual const char  *DataType() const {return fDataType.Data();}
    // Type of data - real or simulated
    virtual void    SetDataType(const char *data="simulated") {fDataType=data;}
    // Set parameters options: "same" or read from "file" or "SetInvalid" or...
    virtual void   SetParamOptions(const char*,const char*) = 0;
    // Set noise parameters 
    virtual void   SetNoiseParam(Double_t, Double_t) = 0;
    // Number of parameters to be set
    virtual  void   SetNDetParam(Int_t) = 0;
    // Set detector parameters: gain, coupling ...
    virtual  void   SetDetParam(Double_t *) = 0;

    // Parameters options
    virtual void   ParamOptions(char *,char*) const = 0;
    virtual Int_t  NDetParam() const = 0;
    virtual void   GetDetParam(Double_t *) const = 0;
    virtual void   GetNoiseParam(Double_t&, Double_t&) const = 0;

    // Zero-suppression option - could be 1D, 2D or non-ZeroSuppressed
    virtual void   SetZeroSupp(const char*) = 0;
    // Get zero-suppression option
    virtual const char *ZeroSuppOption() const = 0;
     // Set thresholds
    virtual void   SetThresholds(Double_t, Double_t) = 0;
    virtual void   Thresholds(Double_t &, Double_t &) const = 0;

    // Set filenames
    virtual void SetFilenames(const char *f1="",const char *f2="",
                              const char *f3=""){
	// Set filenames - input, output, parameters ....
	fFileName1=f1; fFileName2=f2; fFileName3=f3;}
    // Filenames
    virtual void   Filenames(char* input,char* baseline,char* param) {
        strcpy(input,fFileName1.Data());  strcpy(baseline,fFileName2.Data());  
        strcpy(param,fFileName3.Data());}

    virtual Double_t DriftSpeed() const {return SpeedElectron();};
    // set output option
    virtual void    SetOutputOption(Bool_t write=kFALSE) {fWrite = write;}
	
    virtual Bool_t  OutputOption() const {return fWrite;}
    virtual Bool_t  Do10to8() const {return kTRUE;}
    virtual void    GiveCompressParam(Int_t *) const =0;
    //
    // Detector type response methods
    // Set number of sigmas over which cluster disintegration is performed
    virtual void    SetNSigmaIntegration(Double_t) = 0;
    // Get number of sigmas over which cluster disintegration is performed
    virtual Double_t NSigmaIntegration() const = 0;
    // Set number of bins for the gaussian lookup table
    virtual void    SetNLookUp(Int_t) = 0;
    // Get number of bins for the gaussian lookup table
    virtual Int_t GausNLookUp() const {return 0;}
    // Get scaling factor for bin i-th from the gaussian lookup table
    virtual Double_t GausLookUp(Int_t) const {return 0.;}
    // Set sigmas of the charge spread function
    virtual void    SetSigmaSpread(Double_t, Double_t) = 0;
    // Get sigmas for the charge spread
    virtual void    SigmaSpread(Double_t &,Double_t &) const = 0;
    // Pulse height from scored quantity (eloss)
    virtual Double_t IntPH(Double_t) const {return 0.;}
    // Charge disintegration
    virtual Double_t IntXZ(AliITSsegmentation *) const {return 0.;}
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
    virtual Double_t SigmaDiffusion3D(Double_t  l) const ;
    // Returns the Gaussian sigma == <x^2 +y^2+z^2> [cm^2] due to the
    // defusion of electrons or holes through a distance l [cm] caused by an
    // applied voltage v [volt] through a distance d [cm] in any material at a
    // temperature T [degree K].
    virtual Double_t SigmaDiffusion2D(Double_t l) const ;
    // Returns the Gaussian sigma == <x^2+z^2> [cm^2] due to the defusion of
    // electrons or holes through a distance l [cm] caused by an applied
    // voltage v [volt] through a distance d [cm] in any material at a
    // temperature T [degree K].
    virtual Double_t SigmaDiffusion1D(Double_t l) const ;
    // Prints out the content of this class in ASCII format.
    virtual void Print(ostream *os) const;
    // Reads in the content of this class in the format of Print
    virtual void Read(istream *is);
 protected:
    void NotImplemented(const char *method) const {if(gDebug>0)
         Warning(method,"This method is not implemented for this sub-class");}

    TString  fDataType;   // data type - real or simulated
 private:
    Double_t fdv;  // The parameter d/v where d is the disance over which the
                   // the potential v is applied d/v [cm/volts]
    Double_t fN;   // the impurity consentration of the material in #/cm^3
    Double_t fT;   // The temperature of the Si in Degree K.
    Double_t fGeVcharge; // Energy to ionize (free an electron).
    TString  fFileName1;        // input keys : run, module #
    TString  fFileName2;        // baseline & noise val or output code
                                // signal or monitored bgr.
    TString  fFileName3;        // param values or output coded signal
    Bool_t     fWrite;          // Write option for the compression algorithms

    ClassDef(AliITSresponse,2) // Detector type response virtual base class 
};
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITSresponse &source);
istream& operator>>(istream &os,AliITSresponse &source);
#endif
