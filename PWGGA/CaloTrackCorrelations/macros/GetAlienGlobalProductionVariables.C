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
void GetAlienGlobalProductionVariables(Bool_t & simulation, 
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
    if     (period.Contains("18")) year = 2018;
    else if(period.Contains("17")) year = 2017;
    else if(period.Contains("16")) year = 2016;
    else if(period.Contains("15")) year = 2015;
    else if(period.Contains("13")) year = 2013;
    else if(period.Contains("12")) year = 2012;
    else if(period.Contains("11")) year = 2011;
    else if(period.Contains("10")) year = 2010;
  } 
  
  // Check MC production tag name to match with data year and production name
  if ( simulation && period=="" )
  {
    // 2011 MC productions
    if      ( prodType.Contains("14j") )
    {
      year = 2010;
      period = "LHC10";
    } 
    
    // 2011 MC productions
    else if ( prodType.Contains("14ka1") || prodType.Contains("14k1b") || // 7 TeV jet-jet+gamma
              prodType.Contains("12a15g")|| prodType.Contains("13e4")  || // 7 TeV gamma-jet   
              prodType.Contains("13e5")  || // jet-jet pi0   
              prodType.Contains("13d1")  || // MB+pi0
              prodType.Contains("12f2a") || prodType.Contains("12f2b") || // 7 TeV jet-jet+pi0, jet-jet+hadron
              prodType.Contains("12a15f")|| prodType.Contains("12a15a")|| // 7 TeV and 2.76 jet-jet
              prodType.Contains("12a17") || prodType.Contains("14a1")   ) // Pb-Pb LHC11h
    {
      year = 2011;
      period = "LHC11";
    }
    
    // 2012 MC productions
    else if ( prodType.Contains("13d6")  || prodType.Contains("13d7")  ||
              prodType.Contains("13d8")  || prodType.Contains("13d9")  || //MB
              prodType.Contains("14i")   ||
              prodType.Contains("15h1")  || // 8 TeV min bias PYTHIA8
              prodType.Contains("15h2")  || // 8 TeV min bias PHOJET
              prodType.Contains("16c2")  || // 8 TeV jet-jet   
              prodType.Contains("17g5")   ) // 8 TeV jet-jet+gamma, gamma-jet
    {
      year = 2012;
      period = "LHC12";
    }   
    
    // 2013 MC productions
    else if ( prodType.Contains("LHC13") || prodType.Contains("LHC14") ||
              prodType.Contains("15b1")  ||
              prodType.Contains("15a3")  || // pp 2.76 GJ, JJ
              prodType.Contains("15a1")  ||
              prodType.Contains("16c3")  || // JJ+G, GJ Pythia6
              prodType.Contains("17g6")  || // JJ+G, GJ, + BKG Pythia8
              prodType.Contains("16k1a") ||
              prodType.Contains("15g")    )
    {
      year = 2013;
      period = "LHC13";
    }
   else if ( prodType.Contains("16k")    || prodType.Contains("16h")    ||
             prodType.Contains("17d5")   || prodType.Contains("17d6")   || // general purpose and other 
             prodType.Contains("17d7")   || prodType.Contains("17d8")   || // general purpose and other
             prodType.Contains("17e2")   || prodType.Contains("17h5")   || // PbPb
             prodType.Contains("LHC17i4")|| prodType.Contains("LHC17l1")||
             prodType.Contains("18a5a")  || 
             prodType.Contains("18a7")   || 
             prodType.Contains("18b10")  || 
             prodType.Contains("18b11")   ) // jet-jet PbPb
   {
     year = 2015;
     period = "LHC15";
   }
  
   // 2016 MC productions
   else if ( prodType.Contains("17f")    || prodType.Contains("17e")  ||// general purpose and other
             prodType.Contains("17d")    || prodType.Contains("17l2") ||  
             prodType.Contains("17l6")   || prodType.Contains("17l7") || 
             prodType.Contains("17h6")   || // pPb
             prodType.Contains("17h2")   || // pp 13
             prodType.Contains("17h10b") || 
             prodType.Contains("17h4")   || prodType.Contains("17h8") || 
             prodType.Contains("17h9")   || prodType.Contains("18b3_")||
             prodType.Contains("17i3")    ) // pp 13 TeV jet-jet+gamma, gamma-jet
   {
     year = 2016;
     period = "LHC16";
   }
    
   // 2017 MC productions
   else if ( prodType.Contains("17h")    || prodType.Contains("17l")  || 
             prodType.Contains("17k")    ||
             prodType.Contains("18b8")   || // Jet-Jet 5 TeV
             prodType.Contains("18b10")  || // Gamma-Jet 5 TeV
             prodType.Contains("18c12")  || prodType.Contains("18c13")|| 
             prodType.Contains("18a1")   || prodType.Contains("18a3") || 
             prodType.Contains("18a4")   || prodType.Contains("18a8") || 
             prodType.Contains("18a9")    )
   {
     year = 2017;
     period = "LHC17";
   }
    else
    {
      year = 2018;
      period = "LHC18";
    }
  } // Prod MC names
  
}
