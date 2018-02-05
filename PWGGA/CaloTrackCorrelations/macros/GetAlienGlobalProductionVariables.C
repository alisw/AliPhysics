/// \file GetAlienGlobalProductionVariables.C
/// \ingroup CaloTrackCorrMacros
/// \brief Get year, collision type, mc/data type and period from alien global variables
///
/// Get year, collision type, mc/data type and period from alien global variables
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)

// Set includes for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TString.h>
#include <TSystem.h>

#endif

///
/// Main method 
///
/// The options that can be passed to the macro are:
/// \param simulation : bool, true MC, false data
/// \param col        : string with pp, pPb, PbPb return. If already set nothing to recover
/// \param period     : string with period return
/// \param year       : int, year, check in case of data set by hand for MC
/// \param print      : bool to print recovered parameters
///
void GetAlienGlobalProductionVariables(Bool_t simulation, 
                                       TString & col, TString & period, Int_t & year, 
                                       Bool_t print = kFALSE)
{
  TString colType  = gSystem->Getenv("ALIEN_JDL_LPMINTERACTIONTYPE");
  TString prodTag  = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTAG");
  TString prodType = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTYPE");
    
  if(col=="") // Check the alien environment 
  {
    if      (colType.Contains( "PbPb")) col = "PbPb"; 
    else if (colType.Contains( "XeXe")) col = "PbPb"; 
    else if (colType.Contains( "AA"  )) col = "PbPb"; 
    else if (colType.Contains( "pA"  )) col = "pPb"; 
    else if (colType.Contains( "Ap"  )) col = "pPb";     
    else if (colType.Contains( "pPb" )) col = "pPb"; 
    else if (colType.Contains( "Pbp" )) col = "pPb"; 
    else if (colType.Contains( "pp"  )) col = "pp" ; 
    
    // Check if production is MC or data, of data recover period name
    if   ( prodType.Contains("MC") ) simulation = kTRUE;
    else                             simulation = kFALSE;
    
    if   ( !simulation && period!="" ) period = prodTag;
    
    // print check on global settings once
    if(print) 
      printf("GetAlienGlobalProductionVariables() - Get the data features from global parameters: "
             "collision <%s> (<%s>), period <%s>, tag <%s>, type <%s>, MC bool <%d> \n",
             colType.Data(),col.Data(),
             period.Data(),prodType.Data(),prodTag.Data(),simulation);
  }
  
  if ( year < 2009 && !simulation )
  {
    if     (period.Contains("16")) year = 2016;
    else if(period.Contains("15")) year = 2015;
    else if(period.Contains("13")) year = 2013;
    else if(period.Contains("12")) year = 2012;
    else if(period.Contains("11")) year = 2011;
    else if(period.Contains("10")) year = 2010;
  } 
  
}
