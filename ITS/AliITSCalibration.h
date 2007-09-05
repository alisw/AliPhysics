#ifndef ALIITSCALIBRATION_H
#define ALIITSCALIBRATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////
// Base ITS calibration class               //
//////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>
#include "AliLog.h"
#include "AliITSresponse.h"
#include "AliITSMapSDD.h"

class AliITSsegmentation;
class TF1;
class AliITSgeom;


/////////////////////////////////////////////
//                                         //
// ITS calibration virtual base class      //
/////////////////////////////////////////////
class AliITSCalibration : public TObject {
 public:
    // Default Constructor
    AliITSCalibration();
    // Standard Constructor
    AliITSCalibration(Double_t Thickness);

    // Destructor.
    virtual ~AliITSCalibration() {;}
    //
    // Configuration methods
    //
    // fGeVcharge is set by default 3.6e-9 GeV See for ex. PDG 2004.
    virtual void SetGeVToCharge(Double_t gc=3.6e-9){fGeVcharge = gc;}
    // Returns the value fGeVcharge
    virtual Double_t GetGeVToCharge() const {return fGeVcharge;}
    // Converts deposited energy to number of electrons liberated
    virtual Double_t GeVToCharge(Double_t gev) const {return gev/fGeVcharge;}
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
    //  virtual void   SetParamOptions(const char*,const char*) = 0;
    // Set noise parameters 
    virtual void   SetNoiseParam(Double_t, Double_t) = 0;
    // Number of parameters to be set
    virtual  void   SetNDetParam(Int_t) = 0;
    // Set detector parameters: gain, coupling ...
    virtual  void   SetDetParam(Double_t *) = 0;

    // Parameters options
    //  virtual void   ParamOptions(char *,char*) const = 0;
    virtual Int_t  NDetParam() const = 0;
    virtual void   GetDetParam(Double_t *) const = 0;
    virtual void   GetNoiseParam(Double_t&, Double_t&) const = 0;
    virtual void   SetThresholds(Double_t, Double_t) = 0;
    virtual void   Thresholds(Double_t &, Double_t &) const = 0;
    virtual void   SetMapA(Int_t, AliITSMapSDD*) {AliError("This method must be implemented in a derived class");}
    virtual void   SetMapT(Int_t, AliITSMapSDD*) {AliError("This method must be implemented in a derived class");}
    virtual Double_t DriftSpeed() const {return SpeedElectron();};
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
    virtual Double_t SigmaDiffusion3D(Double_t  l) const;
    // Returns the Gaussian sigma == <x^2 +y^2+z^2> [cm^2] due to the
    // defusion of electrons or holes through a distance l [cm] caused by an
    // applied voltage v [volt] through a distance d [cm] in any material at a
    // temperature T [degree K].
    virtual Double_t SigmaDiffusion2D(Double_t l) const;
    // Returns the Gaussian sigma == <x^2+z^2> [cm^2] due to the defusion of
    // electrons or holes through a distance l [cm] caused by an applied
    // voltage v [volt] through a distance d [cm] in any material at a
    // temperature T [degree K].
    virtual Double_t SigmaDiffusion1D(Double_t l) const;
    // Compute the thickness of the depleted region in a Si detector, version A
    virtual Double_t DepletedRegionThicknessA(Double_t dopCons,
                                              Double_t voltage,
                                              Double_t elecCharge,
                                              Double_t voltBuiltIn=0.5)const;
    // Compute the thickness of the depleted region in a Si detector, version B
    virtual Double_t DepletedRegionThicknessB(Double_t resist,Double_t voltage,
                                              Double_t mobility,
                                              Double_t voltBuiltIn=0.5,
                                              Double_t dielConst=1.E-12)const;
    // Computes the temperature dependance of the reverse bias current
    virtual Double_t ReverseBiasCurrent(Double_t temp,Double_t revBiasCurT1,
                                    Double_t tempT1,Double_t energy=1.2)const;
    // Prints out the content of this class in ASCII format.
    virtual void Print(ostream *os) const;
    // Reads in the content of this class in the format of Print
    virtual void Read(istream *is);
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual Int_t Read(const char *name) {return TObject::Read(name);}

    void SetResponse(AliITSresponse* response) {fResponse=response;}
    AliITSresponse* GetResponse() const {return fResponse;}

    virtual void SetDiffCoeff(Float_t p1, Float_t p2) {fResponse->SetDiffCoeff(p1,p2);}
    virtual void DiffCoeff(Float_t &diff,Float_t &diff1) const {fResponse->DiffCoeff(diff,diff1);} 
    virtual void SetFilenames(const char *f1="",const char *f2="",
				 const char *f3="") {fResponse->SetFilenames(f1,f2,f3);}    
    virtual void Filenames(char* input,char* baseline,char* param) {fResponse->Filenames(input,baseline,param);}
    virtual void SetOutputOption(Bool_t write=kFALSE) {fResponse->SetOutputOption(write);}
    virtual Bool_t OutputOption() const {return fResponse->OutputOption();} 

 protected:
    AliITSCalibration(const AliITSCalibration &ob); // copy constructor
    AliITSCalibration& operator=(const AliITSCalibration& source); // ass.
    void NotImplemented(const char *method) const {if(gDebug>0)
         Warning(method,"This method is not implemented for this sub-class");}

    TString  fDataType;   // data type - real or simulated
    
    Double_t fdv;  // The parameter d/v where d is the disance over which the
                   // the potential v is applied d/v [cm/volts]
    Double_t fN;   // the impurity consentration of the material in #/cm^3
    Float_t fT;   // The temperature of the Si in Degree K.
    Double_t fGeVcharge; // Energy to ionize (free an electron) in GeV
    AliITSresponse* fResponse; //! ptr to base response obj. It is not
                 // deleted here but in AliITSDetTypeSim and AliITSDetTypeRec

    ClassDef(AliITSCalibration,1) // Detector type response virtual base class 
};
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITSCalibration &source);
istream& operator>>(istream &os,AliITSCalibration &source);
#endif
