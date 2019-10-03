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
  if ( col != "" && year >= 2009 )
  {
    printf("GetAlienGlobalProductionVariables() - Use values set via configuration:"
           " collision <%s>, period <%s>, year %d, MC %d\n",
           col.Data(), period.Data(), year, simulation);
    return;
  }
  
  TString colType  = gSystem->Getenv("ALIEN_JDL_LPMINTERACTIONTYPE");
  TString prodTag  = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTAG");
  TString prodType = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTYPE");
       
  // In case of environmental name for child meta data
  if ( colType == "" )
  {
    printf("GetAlienGlobalProductionVariables() - Default environment not found? check childs\n");
    for(Int_t ichild = 0; ichild < 100; ichild++)
    {
      colType  = gSystem->Getenv(Form("ALIEN_JDL_child_%d_LPMINTERACTIONTYPE",ichild));
      prodTag  = gSystem->Getenv(Form("ALIEN_JDL_child_%d_LPMPRODUCTIONTAG"  ,ichild));
      prodType = gSystem->Getenv(Form("ALIEN_JDL_child_%d_LPMPRODUCTIONTYPE" ,ichild));
      
      printf("\t child %d col <%s>, tag <%s>, type <%s>\n",
             ichild,colType.Data(),prodTag.Data(),prodType.Data());
      
      if ( colType != "" ) break;
    }
  }
  
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
  
  if   ( !simulation && period == "" ) period = prodTag;
  
  // print check on global settings once
  if ( print ) 
    printf("GetAlienGlobalProductionVariables() - Get the data features from global parameters: "
           "collision <%s> (<%s>), period <%s>, tag <%s>, type <%s>, MC bool <%d> \n",
           colType.Data(),col.Data(),period.Data(),prodTag.Data(),prodType.Data(),simulation);
  
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
    else                           year = 2018;
    
    if ( print ) 
      printf("GetAlienGlobalProductionVariables() -  Data year <%d> \n", year);
  } 
  
  // Check MC production tag name to match with data year and production name
  if ( simulation )
  {
    // 2011 MC productions
    if      ( prodTag.Contains("14j") )
    {
      year   = 2010;
      period = "LHC10";
    } 
    
    // 2011 MC productions
    else if ( prodTag.Contains("14ka1") || prodTag.Contains("14k1b") || // 7 TeV jet-jet+gamma
              prodTag.Contains("12a15g")|| prodTag.Contains("13e4")  || // 7 TeV gamma-jet   
              prodTag.Contains("13e5")  || // jet-jet pi0   
              prodTag.Contains("13d1")  || // MB+pi0
              prodTag.Contains("12f2a") || prodTag.Contains("12f2b") || // 7 TeV jet-jet+pi0, jet-jet+hadron
              prodTag.Contains("12a15f")|| prodTag.Contains("12a15a")|| // 7 TeV and 2.76 jet-jet
              prodTag.Contains("12a17") || prodTag.Contains("14a1")   ) // Pb-Pb LHC11h
    {
      year   = 2011;
      period = "LHC11";
    }
    
    // 2012 MC productions
    else if ( prodTag.Contains("13d6")  || prodTag.Contains("13d7")  ||
              prodTag.Contains("13d8")  || prodTag.Contains("13d9")  || //MB
              prodTag.Contains("14i")   ||
              prodTag.Contains("15h1")  || // 8 TeV min bias PYTHIA8
              prodTag.Contains("15h2")  || // 8 TeV min bias PHOJET
              prodTag.Contains("16c2")  || // 8 TeV jet-jet   
              prodTag.Contains("17g5")   ) // 8 TeV jet-jet+gamma, gamma-jet
    {
      year   = 2012;
      period = "LHC12";
    }   
    
    // 2013 MC productions
    else if ( prodTag.Contains("LHC13") || prodTag.Contains("LHC14") ||
              prodTag.Contains("15b1")  ||
              prodTag.Contains("15a3")  || // pp 2.76 GJ, JJ
              prodTag.Contains("15a1")  ||
              prodTag.Contains("16c3")  || // JJ+G, GJ Pythia6
              prodTag.Contains("17g6")  || // JJ+G, GJ, + BKG Pythia8
              prodTag.Contains("16k1a") ||
              prodTag.Contains("15g")    )
    {
      year   = 2013;
      period = "LHC13";
    }
   else if ( prodTag.Contains("16k")    || prodTag.Contains("16h")    || // LHC15o,n
             prodTag.Contains("16g")    || // LHC15o
             prodTag.Contains("17d5")   || prodTag.Contains("17d6")   || // general purpose and other 
             prodTag.Contains("17d7")   || prodTag.Contains("17d8")   || // general purpose and other
             prodTag.Contains("17e2")   || prodTag.Contains("17h5")   || // PbPb
             prodTag.Contains("LHC17i4")|| prodTag.Contains("LHC17l1")||
             prodTag.Contains("18a5a")  || 
             prodTag.Contains("18a7")   || 
             prodTag.Contains("18b10")  || 
             prodTag.Contains("18b11")   ) // jet-jet PbPb
   {
     year   = 2015;
     period = "LHC15";
   }
  
   // 2016 MC productions
   else if ( prodTag.Contains("17f")    || prodTag.Contains("17e")  ||// general purpose and other
             prodTag.Contains("17d")    || prodTag.Contains("17l2") ||  
             prodTag.Contains("17l6")   || prodTag.Contains("17l7") || 
             prodTag.Contains("17h6")   || // pPb
             prodTag.Contains("17h2")   || // pp 13
             prodTag.Contains("17h10b") || 
             prodTag.Contains("17h4")   || prodTag.Contains("17h8") || 
             prodTag.Contains("17h9")   || prodTag.Contains("18b3_")||
             prodTag.Contains("17i3")    ) // pp 13 TeV jet-jet+gamma, gamma-jet
   {
     year   = 2016;
     period = "LHC16";
   }
    
   // 2017 MC productions
   else if ( prodTag.Contains("17h")    || prodTag.Contains("17l")  || 
             prodTag.Contains("17k")    ||
             prodTag.Contains("18b8")   || // Jet-Jet 5 TeV
             prodTag.Contains("18b10")  || // Gamma-Jet 5 TeV
             prodTag.Contains("18c12")  || prodTag.Contains("18c13")|| 
             prodTag.Contains("18a1")   || prodTag.Contains("18a3") || 
             prodTag.Contains("18a4")   || prodTag.Contains("18a8") || 
             prodTag.Contains("18a9")    )
    {
      year   = 2017;
      period = "LHC17";
    }
    else
    {
      year   = 2018;
      period = "LHC18";
    }
    
    if ( print ) 
      printf("GetAlienGlobalProductionVariables() -  Simulation period <%s>, year <%d> \n",
             period.Data(),year);
  } // Prod MC names
  
}
