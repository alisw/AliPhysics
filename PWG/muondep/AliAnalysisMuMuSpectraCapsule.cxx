/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///
/// Mother class to deal with Spectra after the fitting procedure.
///
/// author: Benjamin Audurier (Subatech)
///

#include "AliAnalysisMuMuSpectraCapsule.h"

ClassImp(AliAnalysisMuMuSpectraCapsule)


#include "AliLog.h"
#include "TObject.h"
#include <TString.h>
#include <iostream>
#include <string>
#include "AliAnalysisMuMuResult.h"
#include "AliAnalysisMuMuJpsiResult.h"
#include "AliAnalysisMuMuBinning.h"
#include "TObjArray.h"
#include "TString.h"
#include <fstream>
#include <string>

using std::cout;
using std::endl;
using std::ifstream;


//_____________________________________________________________________________
AliAnalysisMuMuSpectraCapsule::AliAnalysisMuMuSpectraCapsule() : TObject()
{
  /// Default ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectraCapsule::~AliAnalysisMuMuSpectraCapsule()
{
  // dtor
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuSpectraCapsule::SetConstantFromExternFile(const char* file, Double_t* constantArray, const TString* spectraName)
{
  /// Set member constants depending on centrality bin from an ewternfile.
  /// If values are empty and can be obtained from a graph provided by the AliAnalysisMuMu Framework, this value is set by default
  /// For the PP capsule ine could be :
  ///
  /// #centrality Low  High  lumi.    lumi (stat)  lumi (syst. %)  Trigg  Trigg (local board)   Traj.err.(%)  MC Input (%)  Matching(%)  AccEff  dAccEff  NofJpsi Stat.Jpsi SystJpsi
  /// PP          PP   PP    0.         0.0         0.0            0.0    0.0                    04           03            01           0.0     0.0      0.0     0.0        0.0
  ///
  /// Note that for the PP case, the centrality limits are irrelevant
  ///
  /// For PbPb capsule :
  /// #centrality Low  High  <npart>    d<npart>    TAA           dTAA      sys.AP(%)           Traj.err.(%)   Trigg.err.(%) Matching(%)   AccEff   dAccEff  NofJpsi Stat.Jpsi SystJpsi
  /// V0M         00   10    359        31.2        23.4          0.351     2.00                04             03             01            0.1297   0.00040  105159  1693      488

    // Reset on fConstant
    for (int i = 0; i < 13; ++i) constantArray[i]=0.;
    Bool_t ok= kFALSE;
    AliDebug(1,Form("Reading from file %s",file));

    //________Open file
    ifstream infile(file,std::ios::in);
    TString line;
    TObjArray* lineArray;

    if (infile){
        AliDebug(1, " ==== opening file ==== ");
        // Loop until end of file is reached
        while(infile.eof()!=kTRUE){

            //read the line
            line.ReadLine(infile,kFALSE);
            if (line.BeginsWith("#"))continue;
            AliDebug(1,Form(" Read line : %s",line.Data()));

            // Put the line in a TObjArray
            lineArray = line.Tokenize(" ");


            // Select the good interval. Since interval is written in <binAsString>, just need them to match
            TString centrality   =  static_cast<TObjString*>(lineArray->At(0))->String().Data();
            TString intervalLow  =  TString::Format("%.2f",static_cast<TObjString*>(lineArray->At(1))->String().Atof());
            TString intervalHigh =  TString::Format("%.2f",static_cast<TObjString*>(lineArray->At(2))->String().Atof());
            AliDebug(1,Form(" --__--__-- interval low = %s",intervalLow.Data()));
            AliDebug(1,Form(" --__--__-- interval high = %s",intervalHigh.Data()));
            if (intervalLow.EqualTo("0.00")) intervalLow ="00.00";

            // Select the good interval for PbPb case. Since interval is written in <binAsString>, just need them to match
            if(spectraName->Contains(Form("%s",centrality.Data()))&& spectraName->Contains(Form("%s_%s",intervalLow.Data(),intervalHigh.Data())) && spectraName->Contains(Form("%s_%s",centrality.Data(),intervalLow.Data()))){
                AliDebug(1,Form(" spectraName = %s",spectraName->Data()));
                AliDebug(1,Form(" -- line selected -- "));
                ok = kTRUE;
                break;
            }
            // PP case
            else if(centrality.Contains("PP")){
                AliDebug(1,Form(" spectraName = %s",spectraName->Data()));
                AliDebug(1,Form(" -- line selected -- "));
                ok = kTRUE;
                break;
            }
            else continue;
        }
        infile.close();
        AliDebug(1, " ==== closing file ==== ");

        // Store the value
        for (int i =0 ; i<13 ; i++) {
            constantArray[i]= static_cast<TObjString*>(lineArray->At(i+3))->String().Atof();
        }
        return ok;
    }
    else return ok;
}


//_____________________________________________________________________________
void AliAnalysisMuMuSpectraCapsule::PrintNofWhat(const char* what) const
{
  /// Print whar number for each results on terminal.


  //Check point
  if(!GetSpectra() || strcmp(what,"")==1 )
    {
      AliError("No Spectra or no arguments given !");
      return ;
    }

  // Pointers to handle results and subresults and binning
  AliAnalysisMuMuResult    * result;
  AliAnalysisMuMuJpsiResult* subresult;
  AliAnalysisMuMuResult    * sr;
  AliAnalysisMuMuBinning   ::Range* r;

  // Array to store bins for the while loop
  TObjArray * bins=GetSpectra()->Binning()->CreateBinObjArray();// (intrinseque 'new')
  if (!bins)
  {
    AliError(Form("Cannot find bins"));
    return;
  }

  //Counters and Iterator for bin
  Int_t nofResult = 0;
  TIter nextBin(bins);
  nextBin.Reset();

  // Loop on bins
  //==============================================================================
  while ((r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin())))
  {
    // Make bin a MuMuResult
    result = GetSpectra()->GetResultForBin(*r);
    if (!result)
    {
      AliError(Form("Cannot find result "));
      return;
    }
    AliDebug(1, Form("result(%s) = %p ",result->GetName(),result));

    Int_t nofSubResult = 0; // Counter for subresult
    TIter nextSubResult(result->SubResults());// Iterator for subresults
    nextSubResult.Reset();

    //Some variables
    TString  binAsString(r->AsString());// Usefull for the coming loop

    cout << Form(" -_-_-_-_- %s_%s -_-_-_-_- ",binAsString.Data(),GetSpectraName().Data()) << endl;
    // Loop on subresults
    //==============================================================================
    while ((sr = static_cast<AliAnalysisMuMuResult*>(nextSubResult())))
    {
      // Get our final result
      subresult = static_cast<AliAnalysisMuMuJpsiResult*>(result->SubResult(Form("%s",sr->GetName())));
      if (!subresult)
      {
        AliError(Form("Cannot find subresult "));
        return;
      }
      AliDebug(1,Form("subresult(%s) = %p",sr->GetName(),subresult));

      //Get quantities
      Double_t NofJPsiSub      = subresult->GetValue(what);
      Double_t NofJPsiErrorStat = subresult->GetErrorStat(what);

      //Output messages
      cout << Form(" -------- ") << endl;
      cout << Form(" -- subresult %s :  %.0f +/- %.0f ",sr->GetName(),NofJPsiSub,NofJPsiErrorStat) << endl;
      nofSubResult++;

    }
    cout << Form(" -------- ") << endl;
    cout << Form(" ------ Mean :  %.1f +/- %.1f (%.1f %%) +/- %.1f (%.1f %%) ------ ",
      result->GetValue(what),result->GetErrorStat(what),100*result->GetErrorStat(what)/result->GetValue(what),result->GetRMS(what),100*result->GetRMS(what)/result->GetValue(what)) << endl;
    cout << "" << endl;
    nofResult++;
  }
}
