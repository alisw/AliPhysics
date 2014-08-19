/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliCaloRawAnalyzerComparison.h"

#include "AliCaloRawAnalyzerCrude.h"
#include "AliCaloRawAnalyzerKStandard.h"
#include "AliCaloRawAnalyzerFastFit.h"
#include "AliCaloRawAnalyzerNN.h"
#include "AliCaloRawAnalyzerPeakFinder.h"

#include "AliCaloFitResults.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"

#include <iostream>

using namespace std;

AliCaloRawAnalyzerComparison::AliCaloRawAnalyzerComparison() : fMod(0),
							       fMonCol1(1),
							       fMonRow1(15),
							       fMonCol2(1),
							       fMonRow2(16)
{
  // ctor
  
  fReferenceAnalyzer                       = new  AliCaloRawAnalyzerKStandard();

  AliCaloRawAnalyzerCrude *crude           = new AliCaloRawAnalyzerCrude();
  AliCaloRawAnalyzerKStandard *lms         = new AliCaloRawAnalyzerKStandard();
  AliCaloRawAnalyzerFastFit *fastfit       = new AliCaloRawAnalyzerFastFit();
  AliCaloRawAnalyzerNN *neuralnet          = new AliCaloRawAnalyzerNN();
  AliCaloRawAnalyzerPeakFinder *peakfinder = new AliCaloRawAnalyzerPeakFinder();

  fReferenceAnalyzer->SetIsZeroSuppressed();
  peakfinder        ->SetIsZeroSuppressed();
  lms               ->SetIsZeroSuppressed();
  fastfit           ->SetIsZeroSuppressed();
  neuralnet         ->SetIsZeroSuppressed();
  crude             ->SetIsZeroSuppressed();
  
  fRawAnalyzers.push_back( (AliCaloRawAnalyzer*)crude);
  fRawAnalyzers.push_back( (AliCaloRawAnalyzer*)lms);
  fRawAnalyzers.push_back( (AliCaloRawAnalyzer*)fastfit);
  fRawAnalyzers.push_back( (AliCaloRawAnalyzer*)neuralnet);
  
  fRawAnalyzers.push_back( (AliCaloRawAnalyzer*)peakfinder);
  InitHistograms( fRawAnalyzers, fReferenceAnalyzer );
} 

void 
AliCaloRawAnalyzerComparison::Evaluate( const vector<AliCaloBunchInfo> &bunchvector, const UInt_t altrocfg1,  
					const UInt_t altrocfg2, const int event, const int col, const int row )
{ // evaluate methods
  //  cout << __FILE__ << __LINE__ << endl;
  
  AliCaloFitResults ref = fReferenceAnalyzer->Evaluate( bunchvector, altrocfg1, altrocfg2 );
  //  AliCaloFitResults ref;
  
  for(int i=0; i < fRawAnalyzers.size(); i++ )
  {
    AliCaloFitResults an = fRawAnalyzers.at(i)->Evaluate( bunchvector, altrocfg1,  altrocfg2 );
    
    //     cout << __FILE__ << __LINE__ << ":" << fRawAnalyzers.at(i)->GetAlgoAbbr() << "  col = " << col << " row " << row << "  amp= " <<  an.GetAmp() << endl ;
    
    //   fAmplitudeVsEvent[i]->Fill( event, an.GetAmp() );
    //   fTofVsEvent[i]->Fill( event, an.GetTof() );
    
    fRefAmpVsAnalyzers[i]->Fill(ref.GetAmp(), an.GetAmp() );
    fRefTofVsAnalyzers[i]->Fill(ref.GetTof(), an.GetTof() );
    fAmpDiff[i]->Fill(ref.GetAmp() - an.GetAmp() );
    fTofDiff[i]->Fill(ref.GetTof() - an.GetTof() );
    
    if(col  == fMonCol1 && fMonRow1 == row )
    {
      fAmplitudeVsEvent[i]->Fill( event, an.GetAmp() );
      fTofVsEvent[i]->Fill( event, an.GetTof() );
      fMon1[i] = an;
    }
    
    if(col  == fMonCol2 &&  row == fMonRow2  )
    {
      fMon2[i] = an;
    }
    
    if(row >= 0 && col >= 0 )
    {
      fAmpHistograms[i][col][row]->Fill( an.GetAmp() );
    }
  }
}  

void 
AliCaloRawAnalyzerComparison::EventChanged()
{ // new event
  for(int i=0; i < fRawAnalyzers.size(); i++ )
  {
    /*
     if( ( fMon1[i].GetAmp() > 50 &&   fMon2[i].GetAmp() > 50 ) &&  (  fMon1[i].GetAmp() < 1023  &&   fMon2[i].GetAmp() < 1023 )  )
     {
     cout << __FILE__ << __LINE__ << fRawAnalyzers.at(i)->GetAlgoAbbr() << "\tamp1=" <<  fMon1[i].GetAmp()  << "\tamp2 = " <<  fMon2[i].GetAmp() << endl;
     }
     */
    
    if( fMon1[i].GetAmp()  > 50  &&  fMon2[i].GetAmp() > 50  )
    {
      //  if(  fMon1[i].GetTof() != fMon2[i].GetTof())
      {
        //    fTofResDifferential[i]->Fill( ( fMon1[i].GetTof() - fMon2[i].GetTof())/TMath::Sqrt(2) );
        fTofResDifferential[i]->Fill( ( fMon1[i].GetTof() - fMon2[i].GetTof()) );
        fTofResAbsolute[i]->Fill(  fMon1[i].GetTof() );
      }
    }
  }
}


void 
AliCaloRawAnalyzerComparison::WriteHistograms()
{ // write histograms
  TFile *f = new TFile("comparison2.root", "recreate");
  
  /*
   for(int col=0; col < NZCOLSSMOD; col ++ )
   {
   for(int row=0; row < NXROWSSMOD; row ++ )
   {
   fAmpHistograms[col][row]->Write();
   }
   }
   */
  
  for(int i=0; i < fRawAnalyzers.size(); i++ )
  {
    for(int col=0; col < NZCOLSSMOD; col ++ )
    {
      for(int row=0; row < NXROWSSMOD; row ++ )
	    {
	      fAmpHistograms[i][col][row]->Write();
	    }
    }
    
    fAmplitudeVsEvent[i]->Write();
    fTofVsEvent[i]->Write();
    fRefAmpVsAnalyzers[i]->Write();
    fRefTofVsAnalyzers[i]->Write();
    fAmpDiff[i]->Write();
    fTofDiff[i]->Write();
    fTofResDifferential[i]->Write();
    fTofResAbsolute[i]->Write();
    
  }
  
  f->Close();
}


void
AliCaloRawAnalyzerComparison::InitHistograms( vector <AliCaloRawAnalyzer*> analyzers, AliCaloRawAnalyzer* ref )
{ // init histograms
  char tmpname[256];
  
  /*
   for(int col=0; col < NZCOLSSMOD; col ++ )
   {
   for(int row=0; row < NXROWSSMOD; row ++ )
   {
   sprintf(tmpname, "z(col)%d_x(row)%d_amplitude;Counts;Amplitude/ADC counts", col, row);
   fAmpHistograms[col][row] = new TH1D(tmpname, tmpname, 1024, 0, 1023 );
   }
   }
   */
  
  // TH1D *fAmpHistograms[NZCOLSSMOD][NXROWSSMOD];
  
  for(int i=0; i < analyzers.size(); i++ )
  {
    for(int col=0; col < NZCOLSSMOD; col ++ )
    {
      for(int row=0; row < NXROWSSMOD; row ++ )
	    {
	      sprintf(tmpname, "z(col)%d_x(row)%d_amplitude_%s;Amplitude/ADC counts;Counts", col, row,  analyzers.at(i)->GetAlgoAbbr() );
	      fAmpHistograms[i][col][row] = new TH1D(tmpname, tmpname, 1024, 0, 1023 );
	    }
    }
    
    sprintf(tmpname, "%s_amplitude_vs_event_row%d_col%d;Event;Amplitude/ADC counts", analyzers.at(i)->GetAlgoAbbr(), fMonRow1, fMonCol1 );
    fAmplitudeVsEvent[i]=new TH2D(tmpname, tmpname, 8000, 0, 7999, 1024, 0, 1023 );
    
    sprintf(tmpname, "%s_tof_vs_event_row%d_col%d;Event;Amplitude/ADC counts", analyzers.at(i)->GetAlgoAbbr(), fMonRow1, fMonCol1 );
    
    //     fTofVsEvent[i]=new TH2D(tmpname, tmpname, 8000, 0, 7999, 2000, 3000, 5000);
    fTofVsEvent[i]=new TH2D(tmpname, tmpname, 8000, 0, 7999, 2000, 1000, 4000);
    
    
    sprintf(tmpname, "%s_vs_%s_ampltude; Amplitude_{%s}/ADC counts; Amplitude_{%s}/ADC counts", ref->GetAlgoAbbr(),
            analyzers.at(i)->GetAlgoAbbr(), ref->GetAlgoAbbr(), analyzers.at(i)->GetAlgoAbbr() );
    fRefAmpVsAnalyzers[i] = new TH2D(tmpname, tmpname, 1024, 0, 1023, 1024, 0, 1023);
    sprintf(tmpname, "%s_vs_%s_tof; tof_{%s}/ns; tof_{%s}/ns", ref->GetAlgoAbbr(), analyzers.at(i)->GetAlgoAbbr(), ref->GetAlgoAbbr(), analyzers.at(i)->GetAlgoAbbr() );
    fRefTofVsAnalyzers[i] = new TH2D(tmpname, tmpname, 500, 2000, 4999, 500, 2000, 4999 );
    sprintf( tmpname, "%s-%s amplitude;counts;A_{%s} - A_{%s}",   ref->GetAlgoAbbr(),
            analyzers.at(i)->GetAlgoAbbr(), ref->GetAlgoAbbr(), analyzers.at(i)->GetAlgoAbbr());
    fAmpDiff[i] = new TH1D(tmpname, tmpname, 100, -10, 10 );
    sprintf( tmpname, "%s-%s tof;counts;A_{%s} - A_{%s}",   ref->GetAlgoAbbr(),
            analyzers.at(i)->GetAlgoAbbr(), ref->GetAlgoAbbr(), analyzers.at(i)->GetAlgoAbbr());
    fTofDiff[i] = new TH1D(tmpname, tmpname, 1000, -5000, 5000 );
    sprintf( tmpname, "%s Differential tof resolution (%d, %d) vs (%d, %d);#sigma_{tof}^{%s}/ns;Counts",  analyzers.at(i)->GetAlgoAbbr(), fMonCol1, fMonRow1, fMonCol2, fMonRow2,
            analyzers.at(i)->GetAlgoAbbr() );
    fTofResDifferential[i] = new TH1D(tmpname, tmpname, 1000, -250, 250 );
    sprintf( tmpname, "%s Absolute tof distribution (%d, %d);#sigma_{tof}^{%s}/ns;Counts",  analyzers.at(i)->GetAlgoAbbr(), fMonCol1, fMonRow1, analyzers.at(i)->GetAlgoAbbr() );
    
    fTofResAbsolute[i] = new TH1D(tmpname, tmpname, 2000 , -1000, 7000 );
    
    
  }
}

