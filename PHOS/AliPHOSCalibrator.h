#ifndef ALIPHOSCALIBRATOR_H
#define ALIPHOSCALIBRATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.10  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */


//_________________________________________________________________________
//  Class for performing calibration in PHOS     
//                  
//*-- Author: D.Peressounko (RRC KI & SUBATECH)


// --- ROOT system ---
#include "TTask.h"
#include "TObjArray.h"
#include "TH1F.h"  

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSConTableDB ;
class AliPHOSDigit ;

class AliPHOSCalibrator: public TTask{

public:
  AliPHOSCalibrator() ;          // ctor
  AliPHOSCalibrator(const char* run, const char * title = "Default") ;
  AliPHOSCalibrator(const AliPHOSCalibrator & ctor);
 
  virtual ~AliPHOSCalibrator() ; // dtor

  void AddRun(const char * filename) ; //Add one more file to handle

  virtual void Exec(Option_t * option) ; //Steering method 

  void ScanPedestals(Option_t * option = "append") ;
  void CalculatePedestals(void) ; //calulates pedestals
  void ScanGains(Option_t * opt = "append") ; //calculates gains
  void CalculateGains(void) ; //calculates gains

  void PlotPedestal(Int_t channel) ; //plots distribution of pedestals for given channel
  void PlotPedestals(void) ;
  void PlotGain(Int_t channel) ; //Plot histo with gains for a channel
  void PlotGains() ;             //Plot all gains

  virtual void Print(const Option_t * = "") const ;

  TH1F * PedestalHisto(Int_t channel)
    {return dynamic_cast<TH1F* >(fPedHistos->At(channel)) ;} ;
  TH1F * GainHisto(Int_t channel)
    {return dynamic_cast<TH1F* >(fGainHistos->At(channel)) ;} ;

  
  TH1F * Pedestals(void){return fhPedestals ;}
  TH1F * Gains(void){return fhGains ;}

  //Clean resulting histograms
  void Reset(void){if(fhPedestals)fhPedestals->Reset("ICE");
                   if(fhGains)fhGains->Reset("ICE");}

  //Set energy of beam used to calibrate
  void SetBeamEnergy(Float_t e){fBeamEnergy = e;}

  void SetConTableDB(const char * filename, const char * title = "Default") ;
       //Connection table to convert RawId to AbsId

  void SetNChan(UShort_t nch = 100)
    {fNChan = nch ; }         //Sets number of channels in pedestal histos

  void SetNGainBins(Int_t nbin = 100)
    {fNGainBins = nbin ;}  //Set number of bins in gain histograms

  void SetGainMax(Float_t hmax = 0.01)
    {fGainMax = hmax ;}    //Set range of gain histograms

  void ReadFromASCII(const char * filename) ; //Read gains and pedestals from ascii file

  void WritePedestals(const char * version="v1") ;

  void ReadPedestals(const char * version="v1") ;
		      
  void WriteGains(const char * version="v1") ;

  void ReadGains(const char * version="v1") ;

  AliPHOSCalibrator & operator = (const AliPHOSCalibrator & /*rvalue*/){
    Fatal("operator =","assigment operator is not implemented") ;
    return *this ;
 }



private:
  void Init() ;
 
private:
  TList  * fRunList ;          //list of runs to be handled
  TObjArray * fPedHistos ;     //Resulting histograms of pedestals
  TObjArray * fGainHistos;     //Results of Calibration 

  TH1F *   fhPedestals ;      //Mean values of pedestals for different channels
  TH1F *   fhPedestalsWid ;   //Widths of pedestal distributions for different channels
  TH1F *   fhGains ;          //Final Gains from fitting procedure
  TH1F *   fhGainsWid ;       //Width of final gains from fit

  AliPHOSConTableDB * fctdb ; //!Connection map
  TString  fConTableDB ;      //Name of ConTableDB
  TString  fConTableDBFile ;  //File where ConTableDB is stored

  Float_t  fBeamEnergy ;      //Calibration beam energy
  Float_t  fGainAcceptCorr;   //Maximal deviation from mean Gain (factor)
  Float_t  fAcceptCorr ;      //Maximal deviation of Pedestal from mean for good channel

  Float_t  fGainMax ;         //Range used in Gain histos
  Int_t    fNGainBins ;       //Number of bins in Gain histos

  Int_t    fNch ;             //Number of channels to calibrate
  UShort_t fNChan ;           //Number of bins in pedestal histos

  ClassDef(AliPHOSCalibrator,1)  // description 

};

#endif // AliPHOSCALIBRATOR_H
