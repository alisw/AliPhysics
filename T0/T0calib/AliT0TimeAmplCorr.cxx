#ifndef ALIT0TIMEAMPLCORR_CXX
#define ALIT0TIMEAMPLCORR_CXX


/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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


/**************************************************************************
 * Time - Amplitude correction functions for T0 detector
 * Create correction graphs for PMTs by ESD data trees
 *
 * Creation of correction graphs consist of two independent parts:
 * 1)creating Time-Amplitude plots;
 * 2)making projection and fitting it to creating correction graph for each plot.
 *
 * Plots generated with selections by ESD files with C0TVX, CINT7-B,
 * SPDvertex in range +-fVertexRange, CFD fMeanCFD+-CFDsigma*3.0.
 * Plots are normed my mean CFD taken in 1mip QTC range.
 *
 * EXAMPLE of generating PLOTS:
 *
 *   AliT0TimeAmplCorr *T0ampltimeCorrection = new AliT0TimeAmplCorr(); //creating class
 *   T0ampltimeCorrection->SetVerbose(6); //set printout on max
 *   T0ampltimeCorrection->SetTreeType("ESD"); //set type of tree data
 *   T0ampltimeCorrection->SetQTCType("OLD"); //set QTC type
 *   T0ampltimeCorrection->SetRangeVertexSPD(1); //set vertex range selection
 *
 *   T0ampltimeCorrection->SetNumRUN(245407); //set num run
 *   ESDfile = Form("/tzero/alla/alice/ESDtree/2015/%imuon_calo/", numRUN);
 *   ESDfile += Form("alice_cern.ch_user_a_alla_treeLHC15o_output2445064muon_calo_pass1_000_%03d_AnalysisResults.root",part);
 *   if( !T0ampltimeCorrection->ConnectSourceFile("T0DATA/ESD/ESDtree245407.root", "Alla_histograms/chist") ) return;
 *   if( !T0ampltimeCorrection->MakeQTCCFDplots() ) return;
 *   T0ampltimeCorrection->SaveResultsInFile( "comment" );
 *
 * When QTC-CFD projection created or merged, correction function may be created with MakeATCorrectionForPlots()\n
 * see function description for details\n
 *
 * !!!Plots creation from RAW is not ready yet!!!
 *
 * 03.02.2016
 * Dmitry.Finogeev@cern.ch
 *
 ****************************************************************************/

#include "AliT0TimeAmplCorr.h"
#include <Riostream.h>

using std::cout;
using std::endl;

AliT0TimeAmplCorr::AliT0TimeAmplCorr()
{
    fNumRUN = 0;
    fIsNewQTC = 0;
    fIsPlotsLoaded = 0;

    fMonHists = 0;
    fAssay = 0;
    fQTCCFDplot = 0;
    fProjection = 0;
    fGraphs = 0;
    fChiTrends = 0;

    fTimer = new TStopwatch();
    fTreeChain = new TChain();
    fTreeFile = 0;
    fTree = 0;

    ResetMeanValues();

    SetVerbose();

    SetRangeVertexSPD();
    SetRangeCFD();
    SetRangeQT1();
    SetRangeTVDC();
    SetRangeNewQTC();
    SetQTCType();

    SetTreeType();


    fCFDPlotMin = -100;
    fCFDPlotMax = 100;

}

/*******************************************************************************************************************/
void AliT0TimeAmplCorr::Copy(const AliT0TimeAmplCorr &source)
{
    //clone class data exept results arrays
    //does not checked

    ResetArray( &fMonHists, "delete" );
    ResetArray( &fAssay, "delete" );
    ResetArray( &fQTCCFDplot, "delete" );
    ResetArray( &fProjection, "delete" );
    ResetArray( &fGraphs, "delete" );
    ResetArray( &fChiTrends, "delete" );

    if(fTimer) delete fTimer;
    if(fTreeChain) delete fTreeChain;

    if(fTreeFile)
    {
        if(fTreeFile->IsOpen())fTreeFile->Close();
        delete fTreeFile;
    }
    if(fTree) delete fTree;


    fNumRUN = source.fNumRUN;
    fIsNewQTC = source.fIsNewQTC;
    fIsRAWtree = source.fIsRAWtree;
    fVertexRange = source.fVertexRange;
    fVerbose = source.fVerbose;
    frangeCFD = source.frangeCFD;
    frangeQT1 = source.frangeQT1;
    frangeTVDC = source.frangeTVDC;
    frangeNewQTCmin = source.frangeNewQTCmin;
    frangeNewQTCmax = source.frangeNewQTCmax;
    fIsNewQTC = source.fIsNewQTC;
    fQTCPlotMin = source.fQTCPlotMin;
    fQTCPlotMax = source.fQTCPlotMax;
    fCFDPlotMin = source.fCFDPlotMin;
    fCFDPlotMax = source.fCFDPlotMax;

    for(Int_t pmt=0; pmt<NPMT0; pmt++)
    {
        fNewQTCa[pmt] = source.fNewQTCa[pmt];
        fNewQTCb[pmt] = source.fNewQTCb[pmt];
    }

    fIsPlotsLoaded = 0;

    fMonHists = 0;
    fAssay = 0;
    fQTCCFDplot = 0;
    fProjection = 0;
    fGraphs = 0;
    fChiTrends = 0;

    fTimer = new TStopwatch();

    fTreeChain = new TChain(); //clone this in future
    fTreeFile = 0;
    fTree = 0;

    ResetMeanValues();

}

/*******************************************************************************************************************/

/*******************************************************************************************************************/
/*!
 * \brief Set RUN number
 * \param numrun Run number
 * \return 1 if successfully else 0
 */
Int_t AliT0TimeAmplCorr::SetNumRUN(const Int_t numrun)
{
    if(numrun <= 100000)
    {
      //    if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::SetNumRUN::ERROR: invalid run number: %i",numrun)<<endl;
        fNumRUN = 0;
        return 0;
    }

    fNumRUN = numrun;
    return 1;
}

/*******************************************************************************************************************/
/*!
 * \brief Set tree type
 * \param option "RAW", "ESD"
 * \return 1 if successfully else 0
 */
Int_t AliT0TimeAmplCorr::SetTreeType(TString option)
{
    if( (option.Contains("RAW") + option.Contains("ESD")) != 1 )
         {if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::SetTreeType::ERROR: Wrong tree type parameter: %s (must be \"RAW\" or \"ESD\"); tree type was NOT changed",option.Data())<<endl; return 0;}

    fIsRAWtree = option.Contains("RAW");
    return 1;
}

/*******************************************************************************************************************/
/*!
 * \brief Set QTC type
 * \param QTCtype "NEW", "OLD"
 * \return 1 if successfully else 0
 */
Int_t AliT0TimeAmplCorr::SetQTCType(TString QTCtype)
{
    if( (QTCtype.Contains("OLD") + QTCtype.Contains("NEW")) != 1 )
      {if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::SetQTCType::ERROR: Wrong QTC type parameter: %s (must be \"OLD\" or \"NEW\"); QTC type was NOT changed",QTCtype.Data())<<endl; return 0;}

    fIsNewQTC = QTCtype.Contains("NEW");

    fQTCPlotMin = 0;
    fQTCPlotMax = (fIsNewQTC)?2000:10000;

    return 1;
}

/*******************************************************************************************************************/
AliT0TimeAmplCorr::~AliT0TimeAmplCorr()
{
    ResetArray( &fMonHists, "delete" );
    ResetArray( &fAssay, "delete" );
    ResetArray( &fQTCCFDplot, "delete" );
    ResetArray( &fProjection, "delete" );
    ResetArray( &fGraphs, "delete" );
    ResetArray( &fChiTrends, "delete" );

    if(fTimer) delete fTimer;
    if(fTreeChain) delete fTreeChain;

    if(fTreeFile)
    {
        if(fTreeFile->IsOpen())fTreeFile->Close();
        delete fTreeFile;
    }
    if(fTree) delete fTree;
}

/*******************************************************************************************************************/
void AliT0TimeAmplCorr::ResetArray(TObjArray **array, TString option)
{
    Int_t ioption;

    //if option is incorrect
    if(option.Contains("delete")) ioption = 0;
    else if(option.Contains("recreate")) ioption = 1;
    else {if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::ResetArray::ERROR: array \"%s\" was NOT correctly deleted, wrong option \"%s\"",(*array)->GetName(), option.Data())<<endl; return;}

    //if array is empty
    if(!(*array))
    {
        if(ioption == 1){*array = new TObjArray(); return;}
        if(ioption == 0)return;
    }

    for(Int_t element = 0; element <= (*array)->GetLast(); element++)
        if((*array)->At(element)) delete (*array)->At(element);

    (*array)->Clear();

    delete *array;

    if(ioption == 1){*array = new TObjArray(); return;}
    if(ioption == 0){*array = 0; return;}

}

/*******************************************************************************************************************/
/*!
 * \brief Take TTree "treepathname", from file "filename" if tree is in folder
 * treepathname is the path to tree, if tree is in root, "treepathname" is tne tree name
 * tree type must be set with SetTreeType()
 * \param filename File name
 * \param treename TTree name saved in file
 * \return 1 if successfully else 0
 */
Int_t AliT0TimeAmplCorr::ConnectSourceFile(TString filename, TString treepathname)
{
    if(treepathname == "")
    {
        if(fIsRAWtree)treepathname = "t0tree";
        else treepathname = "T0tree";
    }

    //RAW data
    if(fIsRAWtree)
    {
        if( !IsTTreeWritedCorrectly(filename, treepathname) )
        {
            if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::ConnectSourceFile::ERROR: tree \"%s\" from file \"%s\" can NOT be read",treepathname.Data(), filename.Data())<<endl;
            return 0;
        }

        if( fTreeChain->AddFile( filename.Data(), 0, treepathname.Data() ) )
        {
            if(fVerbose >= 2) cout <<Form("Tree \"%s\" from file \"%s\" was added as RAW data",treepathname.Data(), filename.Data())<<endl;
            return 1;
        }else{
            if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::ConnectSourceFile::ERROR: tree \"%s\" from file \"%s\" was NOT added",treepathname.Data(), filename.Data())<<endl;
            return 0;
        }
    }


    //ESD data
    if(!fIsRAWtree)
    {
        if(fTreeFile)
        {
            if(fTreeFile->IsOpen())fTreeFile->Close();
            delete fTreeFile;
            fTree = 0;
        }

        TObject* currObjTree = 0;
        TFile *fTreeFile = new TFile(filename.Data(),  "READONLY");

        if(!(fTreeFile->IsOpen()))
        {
            if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::ConnectSourceFile::ERROR: File \"%s\" was NOT found", filename.Data())<<endl;
            return 0;
        }

        TList *l_list = (TList*) gDirectory->Get(treepathname.Data());

        if(l_list == 0)
        {
            if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::ConnectSourceFile::ERROR: Directory \"%s\" in file \"%s\" was NOT found",treepathname.Data(), filename.Data())<<endl;
            return 0;
        }

        currObjTree = l_list->First();
        if(currObjTree) fTree = (TTree*)currObjTree;


        if( (fTreeFile->IsOpen())&&(fTree) )
        {
            if(fVerbose >= 2) cout <<Form("Tree \"%s\" from file \"%s\" was loaded as ESD data",treepathname.Data(), filename.Data())<<endl;
            fTreeChain->Print();
            return 1;
        }else{
            if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::ConnectSourceFile::ERROR: tree \"%s\" from file \"%s\" was NOT loaded",treepathname.Data(), filename.Data())<<endl;
            return 0;
        }

    }


}

/*******************************************************************************************************************/
/*!
 * \brief Take from file coefficients for new QTC data
 * RUN number must be set with SetNumRUN(), by default or in case of failure set default values
 * set QTC type as "NEW", so does not need to call SetQTCType("NEW")
 * \param filename
 * \return 1 if successfully else 0
 */
Int_t AliT0TimeAmplCorr::SetNewQTCcoefficients(TString filename)
{

    if((filename != "")&&(fNumRUN != 0))
    {
        TFile *file_coef = new TFile(filename.Data(),  "READONLY");
        TTree *tree_coef = (TTree*)file_coef->Get("coefs");

        if(tree_coef)
        {
            //setting branches
            Int_t currRUN = 0;
            Float_t coeffA[NPMT0] = {1.}, coeffB[NPMT0] = {0.};
            tree_coef->SetBranchAddress("num_run", &currRUN);
            for(Int_t pmt = 0; pmt < NPMT0; pmt++)
            {
                tree_coef->SetBranchAddress(Form("coef_a_PMT%i",pmt+1), &coeffA[pmt]);
                tree_coef->SetBranchAddress(Form("coef_b_PMT%i",pmt+1), &coeffB[pmt]);
            }

            //searching most convenient run
            Int_t collRUN = tree_coef->GetEntries();
            Int_t runDiff = 0, minRUNdiff = -111, runChosen = 0, enteryChoosen = 0;
            Int_t enteryForMeanValues = 0;
            Int_t nValidRUNS = 0, runMax = -111, runMin = 111;
            for(Int_t entry = 0; entry < collRUN; entry++)
            {
                tree_coef->GetEntry(entry);
                if(fNumRUN == 0){ enteryForMeanValues = entry; continue; }
                nValidRUNS++;
                if((fNumRUN < runMin)||(runMin == 111))runMin = fNumRUN;
                if((runMax < fNumRUN)||(runMax == -111))runMax = fNumRUN;

                runDiff = TMath::Abs(currRUN - fNumRUN);
                if((runDiff < minRUNdiff)||(minRUNdiff < 0))
                {minRUNdiff = runDiff; runChosen = currRUN; enteryChoosen = entry;}
            }

            if((runMax - runMin) < minRUNdiff) enteryChoosen = enteryForMeanValues;
            tree_coef->GetEntry(enteryChoosen);

            for(Int_t pmt = 0; pmt <= NPMT0; pmt++)
            {
                fNewQTCa[pmt] = coeffA[pmt];
                fNewQTCb[pmt] = coeffB[pmt];
            }

            if(fVerbose >= 2)
            {
                cout <<endl<<Form( "Coefficients for new QTC loaded from file\"%s\"",filename.Data())<<endl;
                cout <<endl<<Form( "RUN %i chosen, current RUN %i",runChosen, fNumRUN)<<endl;

                for(Int_t pmt=0; pmt<24; pmt++) cout <<Form("PMT_%i a: %.03f, b: %.03f",pmt, fNewQTCa[pmt], fNewQTCb[pmt])<<endl;
            }

            file_coef->Close();
            delete file_coef;

            if(fVerbose >= 2) cout <<Form("Plots will be drawn for new QTC")<<endl;
            fIsNewQTC = 1;

            return 1;

        }
        else if(fVerbose >= 2) cout <<Form("AliT0TimeAmplCorr::SetNewQTCcoefficients::ERROR: TTree \"coefs\" was not loaded from file \"%s\"",filename.Data())<<endl;


    }

    //defaults values from Furs's mail:
    Double_t a_mean[24] = {0.842277, 0.9588, 0.957138, 0.985494, 0.983393, 0.960508,
                           0.929485, 0.941351, 0.987975, 0.994604, 0.95595, 0.964502, 0.964731,
                           0.97356, 0.996221, 0.981255, 0.943558, 1.00606, 0.934549, 0.969384,
                           0.97097, 0.989021, 0.95537, 0.987018};

    Double_t b_mean[24] = {35.431, 14.6214, 3.83161, 48.624, 18.9751, 31.1847,
                           40.1204, 24.1932, 20.8494, 17.2447, 23.0073, 20.1561, 22.3656, 27.0031,
                           19.4473, 18.5315, 20.9247, 15.5292, 42.4446, -3.91989, 11.0756, 16.9568,
                           25.8498, 21.2587};

    if(fVerbose >= 2)cout <<endl<<Form( "Coefficients for new QTC was set by default:")<<endl;

    for(Int_t pmt=0; pmt<24; pmt++)
    {
        fNewQTCa[pmt] = a_mean[pmt];
        fNewQTCb[pmt] = b_mean[pmt];
        if(fVerbose >= 2)cout <<Form("PMT_%i a: %.03f, b: %.03f",pmt, fNewQTCa[pmt], fNewQTCb[pmt])<<endl;
    }


    if(fVerbose >= 2) cout <<Form("Plots will be drawn for new QTC")<<endl;
    SetQTCType("NEW");

    return 0;


}

/*******************************************************************************************************************/
/*!
 * \brief Take QTCCFD plots from root file and merge them to already loaded
 * in case plots have a name "plot_for_pmt1", where 1 is PMT number 1-24,
 * use PlotNameMask = "plot_for_pmt%i"
 * \param filename
 * \return 1 if successfully else 0
 */
Int_t AliT0TimeAmplCorr::MergeQTCCFDplots(TString filename, TString PlotNameMask)
{
    if(!fQTCCFDplot) fQTCCFDplot = new TObjArray();

    TFile *file = new TFile(filename.Data(),"READONLY");
    if(file->IsOpen()) if(fVerbose >= 2) cout <<Form("File \"%s\" was opened for QTCCFD plots collection", filename.Data())<<endl;
    if(!file->IsOpen()){ if(fVerbose >= 2) cout <<Form("File \"%s\" was NOT opened\n", filename.Data())<<endl; return 0;}


    TString pmtname;
    TString sourcePlotName, PlotName;
    TObject *currFilePlot = 0, *currArrayPlot = 0;
    Int_t plotsAdded = 0;

    for(Int_t PMT = 0; PMT < 24; PMT++)
    {
        pmtname = Form("%s_%02d", (PMT<=11)?"C":"A", (PMT<=11)?PMT+1:PMT-11);
        PlotName = Form("%s_QTCCFD", pmtname.Data());
        sourcePlotName = (PlotNameMask == "")? PlotName : Form(PlotNameMask.Data(), PMT);

        currFilePlot = file->FindObjectAny( sourcePlotName.Data() );
        currArrayPlot =  fQTCCFDplot->FindObject( PlotName.Data() );
        if(fVerbose >= 3) cout <<Form("Plot \"%s\" was %sfound", sourcePlotName.Data(), (currFilePlot)?"":"NOT " )<<endl;

        if(!currFilePlot) continue;

        ( (TH1*)currFilePlot)->SetName(PlotName.Data() );

        if(currArrayPlot == 0)
        {
            TObject * FilePlotClone = file->CloneObject(currFilePlot, 0);
            fQTCCFDplot->Add( FilePlotClone );
            if(fVerbose >= 2) cout <<Form("Plot \"%s\" was added; total Entries: %.0f", FilePlotClone->GetName(), ((TH2*)FilePlotClone)->GetEntries() )<<endl;

        }
        else
        {
            TList *mergelist = new TList;
            mergelist->Add(currFilePlot);
            ((TH2*)currArrayPlot)->Merge(mergelist);

            if(mergelist) delete mergelist;
            if(fVerbose >= 2) cout <<Form("Plot \"%s\" with events %.0f was merged; total Entries: %.0f", currFilePlot->GetName(), ((TH2*)currFilePlot)->GetEntries(), ((TH2*)currArrayPlot)->GetEntries() )<<endl;
        }


        plotsAdded++;
    }

    file->Close();
    if(file) delete file;

    if(fVerbose >= 2) cout <<Form("Was read %i QTCCFD plots from file %s", plotsAdded, filename.Data())<<endl;

    if(fVerbose >= 5) fQTCCFDplot->Print();

    fIsPlotsLoaded = 1;
    return plotsAdded;
}

/*******************************************************************************************************************/
/*!
 * \brief Save all Results In File
 * file name consist of parameters markers and comment
 * return constructed file name
 * rewrites existing information in file
 * \param file name comment
 */
TString AliT0TimeAmplCorr::SaveResultsInFile(TString filecomment)
{
    TString filename = Form("ATC_%d_%s_QTC%s_%s_%s.root", fNumRUN, (fIsRAWtree)?"RAW":"ESD",
                            (fIsNewQTC)?"new":"old", (fIsPlotsLoaded)?"corr":"plots", filecomment.Data());

    if(fVerbose >= 2) cout <<Form("\n\nWriting results in file \"%s\":", filename.Data())<<endl;

    TFile *file = new TFile(filename.Data(), "RECREATE");

    if(file->IsOpen())
    {
        if(fVerbose >= 2) cout <<Form("\n\nWriting results in file \"%s\" ... ", filename.Data())<<endl;
    }
    else
    {
        if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::SaveResultsInFile::ERROR: File \"%s\" was NOT created",filename.Data())<<endl;
        return filename;
    }

    file->Close();
    if(file) delete file;

    SaveArrayHistsInFile(fMonHists, filename, "monHists");
    SaveArrayHistsInFile(fQTCCFDplot, filename, "QTCCFDplots");
    SaveArrayHistsInFile(fProjection, filename, "projections");
    SaveArrayHistsInFile(fGraphs, filename, "resultGraphs");
    SaveArrayHistsInFile(fChiTrends, filename, "CHIgraphs");
    SaveArrayHistsInFile(fAssay, filename, "HistAssay");

    return filename;
}

/*******************************************************************************************************************/
/*!
 * \brief Take TTree "treename" from files
 * file name is configured by mask format as "/tree/treeRUN%i/tree_file_%i.root"\n
 * mask must include two "%i" for numrun and file index respectively\n
 * file index varies from ffirst to flast\n
 * \param filemask example: "/tree/treeRUN%i/tree_file_%i.root"
 * \param numrun
 * \param treename
 * \param ffirst
 * \param flast
 * \return 1 if successfully else 0
 */
Int_t AliT0TimeAmplCorr::ConnectSourceFilesByMask(TString filemask, const Int_t numrun, TString treename, Int_t ffirst, Int_t flast)
{

    if(SetNumRUN(numrun) < 1) return 0;
    if(ffirst< 0)ffirst = 0;
    if(flast < ffirst)flast = ffirst;

    //may be filemask check

    Bool_t fileadded;
    Int_t nfilescollect = 0;
    TString currFileName;
    for(Int_t nfile=ffirst; nfile <= flast; nfile++)
    {
        currFileName = Form(filemask.Data(),fNumRUN,nfile);

        fileadded = ConnectSourceFile(currFileName, treename);

        if(fileadded) nfilescollect++;
    }

    return nfilescollect;
}

/*******************************************************************************************************************/
/*!
 * \brief Take TTree "treename" from files in directory.
 * File name must contains oblnamepart
 * \param path
 * \param treename
 * \param oblnamepart
 * \return 1 if successfully else 0
 */
Int_t AliT0TimeAmplCorr::ConnectSourceFilesFromDir(TString path, TString treename, TString oblnamepart)
{

    if(path == "") path = gSystem->pwd();
    if(!path.EndsWith("/")) path += "/";

    Int_t nfilescollect = 0;
    TSystemFile *currFile;
    TString currFileName;
    TSystemDirectory SDirectory(path, path);
    TIter nextfile(SDirectory.GetListOfFiles());
    Bool_t fileadded = 0;

    while ((currFile=(TSystemFile*)nextfile()))
    {
        currFileName = currFile->GetName();
        if (!currFile->IsDirectory() && currFileName.EndsWith("root") && currFileName.Contains(oblnamepart.Data()))

            fileadded = ConnectSourceFile(path + currFileName, treename);

        if(fileadded) nfilescollect++;
    }

    return nfilescollect;
}

/*******************************************************************************************************************/
/*!
 * \brief MergeQTCCFDplotsFromDir merge all plot files from dirrectory
 * if(PlotNameMask != "") QTCname = Form(PlotNameMask, PMT+1)
 * \param filename
 * \param PlotNameMask
 * \return
 */
Int_t AliT0TimeAmplCorr::MergeQTCCFDplotsFromDir(TString path)
{
    if(path == "") path = gSystem->pwd();
    if(!path.EndsWith("/")) path += "/";

    Int_t nfilescollect = 0;
    TSystemFile *currFile;
    TString currFileName;
    TSystemDirectory SDirectory(path, path);
    TIter nextfile(SDirectory.GetListOfFiles());
    Bool_t fileadded = 0;

    while ((currFile=(TSystemFile*)nextfile()))
    {
        currFileName = currFile->GetName();
        if (!currFile->IsDirectory() && currFileName.EndsWith("root"))

            fileadded = ( 0 < MergeQTCCFDplots(path + currFileName) );

        if(fileadded) nfilescollect++;
    }

    return nfilescollect;

}

/*******************************************************************************************************************/
Int_t AliT0TimeAmplCorr::SaveArrayHistsInFile(TObjArray *array, TString filename, TString directoryname)
{
    if( !array )
    {
        if(fVerbose >= 1) cout <<Form("WARING: Histogram array \"%s\" does NOT exist", directoryname.Data())<<endl;
        return 0;
    }

    TFile *file = new TFile(filename.Data(), "UPDATE");

    if(file->IsOpen())
    {
        if(fVerbose >= 2) cout <<Form("Writing \"%s\" in file \"%s\" ... ", directoryname.Data(), filename.Data())<<endl;
    }
    else
    {
        if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::SaveArrayHistsInFile::ERROR: File \"%s\" was NOT opened for UPDATE",filename.Data())<<endl;
        return 0;
    }

    if(array->GetLast() < 0)
    {
        if(fVerbose >= 1) cout <<Form("WARING: Histogram array \"%s\" is empty", directoryname.Data())<<endl;
        return 0;
    }

    file->rmdir(directoryname.Data()); //delete old data to rewrite
    TDirectory *CurrDir = file->mkdir(directoryname.Data()); CurrDir->cd();

    for(Int_t n=0; n<array->GetEntries(); n++)
        if(array->At(n)) array->At(n)->Write();

    file->Close();
    delete file;

    if(fVerbose >= 2) cout <<"done"<<endl;

    return 1;
}

/*******************************************************************************************************************/
Int_t AliT0TimeAmplCorr::IsTTreeWritedCorrectly(TString filename, TString treename)
{
    TFile *file = new TFile(filename.Data(), "READONLY");
    if(!file->IsOpen()) return 0;

    TTree *tree = (TTree*)file->Get(treename.Data());
    //tree->Print();

    TString branchNames[7]={"_CFD_0","_QT0_0","_QT1_0","_QTC00_0","_QTC01_0","_QTC10_0","_QTC11_0"};
    TBranch *currBranch;
    TString treePMTname, currBranchName;
    Int_t basketsLoaded;

    for (Int_t pmt=0;pmt<24;pmt++)
    {
        treePMTname = Form("%s_%d", (pmt<=11)?"C":"A", (pmt<=11)?pmt+1:pmt-11);
        for(Int_t branch=0;branch<7;branch++)
        {
            currBranchName = Form("T0_%s%s",treePMTname.Data(), branchNames[branch].Data());
            currBranch = tree->GetBranch( currBranchName.Data() );
            basketsLoaded = currBranch->LoadBaskets();
            if(basketsLoaded <= 0)
            {
                if(fVerbose >= 1)cout <<Form("WARING: baskets was NOT loaded from tree \"%s\" in file \"%s\" for branch \"%s\": %i",
                                             treePMTname.Data(), filename.Data(), currBranchName.Data(), basketsLoaded )<<endl;
                file->Close(); delete file; return 0;
            }
        }
    }

    file->Close(); delete file; return 1;
}

/*******************************************************************************************************************/
void AliT0TimeAmplCorr::ArrayBubbleSortInc(Double_t *const array, const Int_t kElements)
{
    Double_t swap = 0; //for array elements swap
    Int_t swappedTimes = 0; //for break

    while(0 < swappedTimes)
    {
        swappedTimes = 0;
        for(Int_t currelement=0; currelement < (kElements-1); currelement++)
            if(array[currelement+1] < array[currelement])
            {
                swap = array[currelement];
                array[currelement] = array[currelement+1];
                array[currelement+1] = swap;
                swappedTimes++;
            } //for() if(){ }
    } //while
}

/*******************************************************************************************************************/
Int_t AliT0TimeAmplCorr::FitHistogrammInFirstMaxPeak(TH1*const historgamm, Double_t &mean, Double_t &sigma,const Double_t peaksRangeMin,
                                                     const Double_t peaksRangeMax, const Double_t peakRange, const Double_t peakRMSRange)
{
    const Double_t kPeakLeftMinValueRatio = 0.5;

    historgamm->SetAxisRange(peaksRangeMin, peaksRangeMax); //cutting saturation peak
    Double_t peakMax = historgamm->GetBinCenter( historgamm->GetMaximumBin() );
    Double_t peakMaxValue = historgamm->GetBinContent(historgamm->GetMaximumBin());
    Double_t peakLeftEdge = historgamm->GetBinCenter(historgamm->FindFirstBinAbove(peakMaxValue*kPeakLeftMinValueRatio)-2);
    if(fVerbose >= 4) cout <<Form("FitHistogrammInFirstMaxPeak: RANGE[%.2f, %.2f]; Max %.2f (%.2f); Peak left edge %.2f",
                                  peaksRangeMin, peaksRangeMax, peakMax, peakMaxValue, peakLeftEdge)<<endl;

    historgamm->SetAxisRange(peakLeftEdge, peakMax+peakRange); //cutting saturation peak
    Double_t peakRMS = historgamm->GetRMS();

    if(fVerbose >= 4) cout <<Form("FitHistogrammInFirstMaxPeak: Check RANGE[%.2f, %.2f]; RMS %.2f (%.2f); FIT[%.2f, %.2f]",
                                  peakLeftEdge, peakMax+peakRange, peakRMS, peakRMSRange, peakLeftEdge, peakMax+peakRMS*peakRMSRange)<<endl;


    Int_t result = FitHistogramm(historgamm, mean, sigma, peakLeftEdge, peakMax+peakRMS*peakRMSRange);

    historgamm->SetAxisRange(peaksRangeMin, peaksRangeMax); //cutting saturation peak

    return result;
}

/*******************************************************************************************************************/
Int_t AliT0TimeAmplCorr::FitHistogrammInRMSRange(TH1 *const histogramm, Double_t &mean, Double_t &sigma, const Double_t RMSrange)
{
    if( !histogramm ){if(fVerbose >= 1)
            cout <<Form("AliT0TimeAmplCorr::FitHistogrammInRMSRange::ERROR: histogramm pointer is NULL")<<endl; return 0;}

    Double_t histMean = histogramm->GetMean();
    Double_t histRMS = histogramm->GetRMS();

    return FitHistogramm(histogramm, mean, sigma, histMean-RMSrange*histRMS, histMean+RMSrange*histRMS);
}

/*******************************************************************************************************************/
Int_t AliT0TimeAmplCorr::FitHistogramm(TH1 *const histogramm, Double_t &mean, Double_t &sigma, Double_t rangeXmin, Double_t rangeXmax)
{
    static const Int_t kMinFitEvents = 1;       //minimum events in range for fitting
    static const Int_t kEventsNoFit = 50;

    TString fitoptions[] = {"QR", "QRL", "QRI", "END"};

    if( !histogramm ){if(fVerbose >= 1)
            cout <<Form("AliT0TimeAmplCorr::FitHistogramm::ERROR: histogramm pointer is NULL")<<endl; return 0;}


    if(rangeXmax <= rangeXmin)
    {
        rangeXmin = histogramm->GetXaxis()->GetXmin();
        rangeXmax = histogramm->GetXaxis()->GetXmax();
    }

    histogramm->SetAxisRange(rangeXmin, rangeXmax);

    //if too low events for fitting returning 0
    Int_t intBinMin = histogramm->FindBin(rangeXmin);
    Int_t intBinMax = histogramm->FindBin(rangeXmax);
    Int_t integral =  histogramm->Integral(intBinMin, intBinMax);

    if( integral < kMinFitEvents){if(fVerbose >= 1)
            cout <<Form("AliT0TimeAmplCorr::FitHistogramm::ERROR: Too low events in \"%s\" histogram: %i",histogramm->GetName(), integral)<<endl; return 0;}

    Double_t histMean = histogramm->GetMean();
    Double_t histRMS = histogramm->GetRMS();


    //if low events, only RMS and mean;
    if( integral <= kEventsNoFit )
    {
        mean = histMean;
        sigma = histRMS;
        return 1;
    }

    TF1 *currfitfunction = 0;
    TF1 *bestfitfunction = 0;
    Int_t bestOption = -111;

    Double_t bestChisquare = -111., currChisquare = 0;

    for(Int_t opt = 0; !fitoptions[opt].Contains("END"); opt++){

        currfitfunction = new TF1("fitFunction", "gaus", rangeXmin, rangeXmax);

        histogramm->Fit(currfitfunction, fitoptions[opt], "", rangeXmin, rangeXmax);
        currChisquare = currfitfunction->GetChisquare();

        if(fVerbose >= 5)cout<<Form("current fiting option: [%.2f, %.2f] \"%s\", chi-square: %.1f", rangeXmin, rangeXmax,fitoptions[opt].Data(), currChisquare)<<endl;
        if( ( (currChisquare < bestChisquare)&&(currChisquare > 0.0001) )||( (bestChisquare < 0.)&&(currChisquare > 0.0001)) )
        {
            bestChisquare = currChisquare;
            if(bestfitfunction) delete bestfitfunction;
            bestfitfunction = currfitfunction;
            bestOption = opt;
        }
        else
        {
            if(currfitfunction) delete currfitfunction;
        }
    }


    if( 0. < bestChisquare)
    {
        histogramm->Fit(bestfitfunction, fitoptions[bestOption], "", rangeXmin, rangeXmax);
        if(fVerbose >= 5)cout<<Form("best fit: %.2f, %i", bestChisquare, bestOption)<<endl;

        mean = bestfitfunction->GetParameter(1);
        sigma = bestfitfunction->GetParameter(2);
    }
    else
    {
        if(fVerbose >= 5)cout<<Form("fitting failure for all options, mean and RMS returned")<<endl;
        mean = histMean;
        sigma = histRMS;
    }

    return 1;

}

/*******************************************************************************************************************/
Int_t AliT0TimeAmplCorr::FitHistogrammInMaxPeak(TH1 *const histogramm, Double_t &mean, Double_t &sigma)
{
    static const Int_t kMinFitEvents = 1;       //minimum events in range for fitting
    static const Int_t kEventsNoFit = 500;
    static const Int_t kMinPeakRMS = 10; //ch

    TString fitoptions[] = {"QR", "QRL", "QRI", "END"};

    Double_t histMean = histogramm->GetMean();
    Double_t histRMS = histogramm->GetRMS();
    if(histRMS < kMinPeakRMS) histRMS = kMinPeakRMS;
    Double_t histMaximumX = histogramm->GetBinCenter( histogramm->GetMaximumBin() );

    Double_t rangeMin = histMaximumX - histRMS*1.;
    Double_t rangeMax = histMaximumX + histRMS*1.;

    histogramm->SetAxisRange(rangeMin, rangeMax);

    histMean = histogramm->GetMean();
    histRMS = histogramm->GetRMS();
    //if(histRMS < minPeakRMS) histRMS = minPeakRMS;
    histMaximumX = histogramm->GetBinCenter( histogramm->GetMaximumBin() );


    //if too low events for fitting returning 0
    Int_t intBinMin = histogramm->FindBin(histMaximumX-histRMS*1);
    Int_t intBinMax = histogramm->FindBin(histMaximumX+histRMS*1);
    Int_t integral =  histogramm->Integral(intBinMin, intBinMax);

    if( integral < kMinFitEvents){if(fVerbose >= 1)
            cout <<Form("AliT0TimeAmplCorr::FitHistogrammInMaxPeak::ERROR: Too low events in \"%s\" histogram: %i",histogramm->GetName(), integral)<<endl; return 0;}

    //if low events, only RMS and mean;
    if( integral <= kEventsNoFit )
    {
        mean = histMean;
        sigma = histRMS;
        return 1;
    }

    //TF1 *currfitfunction = new TF1("fitFunction", "gaus", RangeMin, RangeMax);
    TF1 *currfitfunction = 0;
    TF1 *bestfitfunction = 0;
    Int_t bestOption = -111;
    Double_t bestHistFitRange = -111.;

    Double_t bestChisquare = -111., currChisquare = 0;

    for(Double_t histFitRange = 1.; histFitRange <= 2.; histFitRange += 0.2)
        for(Int_t opt = 0; !fitoptions[opt].Contains("END"); opt++){
            rangeMin = histMaximumX - histFitRange*histRMS;
            rangeMax = histMaximumX + histFitRange*histRMS;

            currfitfunction = new TF1("fitFunction", "gaus", rangeMin, rangeMax);

            histogramm->Fit(currfitfunction, fitoptions[opt], "", rangeMin, rangeMax);
            currChisquare = currfitfunction->GetChisquare() / histFitRange;

            if(fVerbose >= 6)cout<<Form("current fiting option: \"%s\" range: %.2f, chi-square: %.1f",fitoptions[opt].Data(), histFitRange, currChisquare)<<endl;
            if( ( (currChisquare < bestChisquare)&&(currChisquare > 0.0001) )||( (bestChisquare < 0.)&&(currChisquare > 0.0001)) )
            {
                bestChisquare = currChisquare;
                if(bestfitfunction) delete bestfitfunction;
                bestfitfunction = currfitfunction;
                bestHistFitRange = histFitRange;
                bestOption = opt;
            }
            else
            {
                if(currfitfunction) delete currfitfunction;
            }
        }

    Double_t fitMean;
    Double_t fitSigma;

    if( 0. < bestChisquare)
    {
        histogramm->Fit(bestfitfunction, fitoptions[bestOption], "", histMean-histRMS*bestHistFitRange, histMean+histRMS*bestHistFitRange);
        if(fVerbose >= 6)cout<<Form("best fit: %.2f, %i, %.1f", bestChisquare, bestOption, bestHistFitRange)<<endl;

        fitMean = bestfitfunction->GetParameter(1);
        fitSigma = bestfitfunction->GetParameter(2);
    }
    else
    {
        if(fVerbose >= 6)cout<<Form("fitting failure for all options")<<endl;
        fitMean = 0.;
        fitSigma = 0.;
    }



    if(bestfitfunction) delete bestfitfunction;

    histogramm->SetAxisRange(histMaximumX-histRMS, histMaximumX+histRMS);
    histMean = histogramm->GetMean();
    histRMS = histogramm->GetRMS();

    Double_t testReg = histRMS;

    //integral by mean
    intBinMin = histogramm->FindBin(histMean-testReg);
    intBinMax = histogramm->FindBin(histMean+testReg);
    Int_t integral_mean =  histogramm->Integral(intBinMin, intBinMax);

    //integral by fit
    intBinMin = histogramm->FindBin(fitMean-testReg);
    intBinMax = histogramm->FindBin(fitMean+testReg);
    Int_t integral_fit =  histogramm->Integral(intBinMin, intBinMax);

    if(fVerbose >= 6)cout<<Form("[fit_mean %.2f, sigma %.2f] [max %.2f mean %.2f, RMS %.2f +- %.2f] Integral by fit: %i, Integral by mean: %i",
                                fitMean, fitSigma, histMaximumX, histMean, histRMS,histogramm->GetRMSError() ,integral_fit, integral_mean)<<endl;


    //if( TMath::Abs(histMaximumX - histMean) < TMath::Abs(histMaximumX - FitMean) )
    if(integral_fit < integral_mean)
    {
        mean = histMean;
        sigma = histRMS;
    }
    else
    {
        mean = fitMean;
        sigma = fitSigma;
    }

    return 1;
}

/*******************************************************************************************************************/
Int_t AliT0TimeAmplCorr::SubtractGraphFromHist(TH1 *const histogramm, const TGraph *graph)
{
    if(!graph || !histogramm) return 0;

    if(fVerbose >= 2)cout <<Form("Subtracting graph \"%s\" from histogram \"%s\" ... ", graph->GetName(), histogramm->GetName())<<endl;

    static const Double_t kMaxFuncValue = 10000.;

    Double_t grafMin = graph->GetXaxis()->GetXmin();
    Double_t grafMax = graph->GetXaxis()->GetXmax();
    Double_t currBinCenter, currBinContent, currBinError, currFuncVal;

    for(Int_t bin = 0; bin <= histogramm->GetNbinsX(); bin++)
    {

        currBinCenter = histogramm->GetBinCenter(bin);
        if( (grafMax <= currBinCenter)||(currBinCenter <= grafMin) ) continue;

        currFuncVal = graph->Eval(currBinCenter);
        if( (kMaxFuncValue <= currFuncVal)||(currFuncVal <= -kMaxFuncValue) ) continue;

        currBinContent = histogramm->GetBinContent(bin);
        if(TMath::Abs(currBinContent) < 0.1) continue;

        currBinError = histogramm->GetBinError(bin);

        currBinContent -= currFuncVal;

        currBinError = currBinError;

        histogramm->SetBinContent(bin, currBinContent);
        histogramm->SetBinError(bin, currBinError);

    }

    return 1;
}


/*******************************************************************************************************************/
TH2* AliT0TimeAmplCorr::CreateCFDQTCplotFromESD(const Int_t pmt)
{


    // for monitor histogramms
    static const Double_t kAmplAbsMin = 0.;
    static const Double_t kAmplAbsMax = 10000.;
    static const Double_t kTimeAbsMin = 10000.;
    static const Double_t kTimeAbsMax = 11000;
    static const Double_t kBinPerUnitTime = 5;
    Double_t binPerUnitAmpl = (fIsNewQTC)?2.:10.; //for monitor hist


    // Ranges values
    static const Double_t kTimeSigmaRange = 3.;   //plot is taken only for main CFD peak mean+-Sigma*TimeSigmaRange
    static const Double_t kAmplMIPSigmaRange = .7;   //plot is taken only for main CFD peak mean+-Sigma*TimeSigmaRange

    // CFD and QTC plots:
    Double_t timeMean, timeSigma, amplMIPmean, amplMIPsigma, timeCFDmipMEAN;//, AmplMean, AmplSigma;

    // for QTC-CFD plot
    Double_t amplMin, amplMax, timeMin, timeMax, amplMIPmin, amplMIPmax;
    static const Double_t kBinPerUnitTimeMip = 5.;

    TString pmtname = Form("%s_%d", (pmt<=11)?"C":"A", (pmt<=11)?pmt+1:pmt-11);  //pmt name for histogram name

    //sting variables for draw function
    TString selectionTRGeq;
    TString vertexeq, selectionVertex;
    TString timeCFDeq, selectionCFDeq;
    TString amplitudeQTCeq, selectionQTCeq, selectionQTCmipEq;

    TString C0TVX_B_selEq = Form("(trigger.GetString().Contains(\"C0TVX-B\"))");
    C0TVX_B_selEq = Form("T0trigger&(1<<0)");
    TString CINT7_B_selEq = Form("(trigger.GetString().Contains(\"CINT7-B\"))");
    TString C0TSC_B_selEq = Form("(trigger.GetString().Contains(\"C0TSC-B\"))");
    C0TSC_B_selEq = "(1)";

    selectionTRGeq = Form("(%s && %s && %s)", C0TVX_B_selEq.Data(), CINT7_B_selEq.Data(), C0TSC_B_selEq.Data());

    vertexeq = "vertexSPD";
    selectionVertex = Form("(%i < %s)&&(%s < %i)", -fVertexRange, vertexeq.Data(), vertexeq.Data(), fVertexRange);

    timeCFDeq = Form("time%i", pmt+1);
    selectionCFDeq = Form("(%f < %s)", 0., timeCFDeq.Data());

    amplitudeQTCeq = (fIsNewQTC)?Form("amp%i_new", pmt+1):Form("amp%i", pmt+1);
    selectionQTCeq = Form("(0 < %s)", amplitudeQTCeq.Data());

    TH1 *ptrTH1curr = 0;   //pointers to created hist

    //determine mean cfd selection by CFD and vertex
    fMonHists->Add( ptrTH1curr = new TH1F( Form("%s_CFD", pmtname.Data()), "title", ceil( (kTimeAbsMax-kTimeAbsMin)/kBinPerUnitTime ),kTimeAbsMin,kTimeAbsMax) );
    fTree->Draw( Form("%s >> %s_CFD", timeCFDeq.Data(), pmtname.Data()),
                 Form("%s && %s && %s", selectionCFDeq.Data(), selectionVertex.Data(), selectionTRGeq.Data()) );
    //if(!FitHistogrammInMaxPeak(ptrTH1curr, TimeMean, TimeSigma))
    //     { if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::CreateCFDQTCplotFromESD::ERROR: \"%s\" is incorrect",ptrTH1curr->GetName())<<endl; return NULL; }
    timeMean = ptrTH1curr->GetMean();
    timeSigma = ptrTH1curr->GetRMS();
    timeMin = timeMean - timeSigma*kTimeSigmaRange;
    timeMax = timeMean + timeSigma*kTimeSigmaRange;

    //ptrTH1curr->GetFunction("fitFunction")->SetRange(TimeMin, TimeMax); //to see selected range on func at hist
    ptrTH1curr->SetTitle( Form("CFD plot, mean %.1f, RMS %.1f, [%.1f, %.1f]",timeMean, timeSigma, timeMin, timeMax) );
    if(fVerbose >= 4)cout <<Form("CFD: mean %.1f, sigma %.1f; [%.1f, %.1f]",timeMean, timeSigma, timeMin, timeMax)<<endl;
    selectionCFDeq += Form("&&(%f < %s)&&(%s < %f)",timeMin, timeCFDeq.Data(), timeCFDeq.Data(), timeMax);

    //QTC plot to determine amplitude range and 1mip QTC position selection by CFD, QTC, vertex
    fMonHists->Add( ptrTH1curr = new TH1F( Form("%s_QTC", pmtname.Data()), "title", ceil( (kAmplAbsMax-kAmplAbsMin)/binPerUnitAmpl ),kAmplAbsMin,kAmplAbsMax) );
    fTree->Draw( Form("%s >> %s_QTC",amplitudeQTCeq.Data(), pmtname.Data()),
                 Form("%s && %s && %s && %s",selectionQTCeq.Data(), selectionCFDeq.Data(), selectionVertex.Data(), selectionTRGeq.Data())  );
    if(!FitHistogrammInFirstMaxPeak(ptrTH1curr, amplMIPmean, amplMIPsigma, 10., 500., 100, 1.5))
    { if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::CreateCFDQTCplotFromESD::ERROR: \"%s\" is incorrect",ptrTH1curr->GetName())<<endl; return 0; }

    amplMIPmin = amplMIPmean - amplMIPsigma*kAmplMIPSigmaRange;
    amplMIPmax = amplMIPmean + amplMIPsigma*kAmplMIPSigmaRange;

    /***
    Double_t QTColdMipTemp[24] = {207.9,  184.72, 143.33, 220.60, 217.01, 256.31,
                                  169.01, 127.83, 262.28, 117.52, 261.42, 190.55,
                                  151.75, 146.28, 171.20, 289.69, 231.71, 201.78,
                                  212.87, 230.55, 181.12, 154.03, 152.29, 146.69};

    Double_t QTCnewMipTemp[24] = {148, 159, 126, 194, 159, 207,
                                  144, 136, 196, 97,  207, 130,
                                  182, 171, 187, 257, 205, 174,
                                  207, 152, 178, 155, 168, 107};



    AmplMIPmean = (fIsNewQTC)?QTCnewMipTemp[PMT]:QTColdMipTemp[PMT];
    //AmplMIPmean = ptrTH1curr->GetBinCenter( ptrTH1curr->GetMaximumBin() );
    AmplMIPmin = AmplMIPmean - 10.;
    AmplMIPmax = AmplMIPmean + 10;
    /***/





    if(ptrTH1curr->GetFunction("fitFunction"))
        ptrTH1curr->GetFunction("fitFunction")->SetRange(amplMIPmax, amplMIPmin); //to see selected range on func at hist

    selectionQTCmipEq = Form("(%f < %s)&&(%s < %f)", amplMIPmin, amplitudeQTCeq.Data(), amplitudeQTCeq.Data(), amplMIPmax);
    amplMin = ptrTH1curr->FindFirstBinAbove();
    amplMax = ptrTH1curr->FindLastBinAbove();



    ptrTH1curr->SetTitle( Form("QTC plot 1mip for selection [%.1f, %.1f]", amplMIPmin, amplMIPmax) );
    if(fVerbose >= 4)cout <<Form("QTC 1mip: mean %.1f, sigma %.1f; [%.1f, %.1f]",amplMIPmean, amplMIPsigma, amplMIPmin, amplMIPmax)<<endl;
    selectionQTCeq += Form("&&(%f < %s)&&(%s < %f)",amplMin, amplitudeQTCeq.Data(), amplitudeQTCeq.Data(), amplMax);

    //determine mean cfd 1MIP selection by CFD and vertex
    fMonHists->Add( ptrTH1curr = new TH1F( Form("%s_CFD_mip", pmtname.Data()), "title", ceil( (timeMax-timeMin)/kBinPerUnitTimeMip ),timeMin,timeMax) );
    fTree->Draw( Form("%s >> %s_CFD_mip", timeCFDeq.Data(), pmtname.Data()),
                 Form("%s && %s && %s && %s", selectionCFDeq.Data(), selectionVertex.Data(), selectionQTCmipEq.Data(), selectionTRGeq.Data()) );
    if(!FitHistogramm(ptrTH1curr, timeCFDmipMEAN, timeSigma))
    { if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::CreateCFDQTCplotFromESD::ERROR: \"%s\" is incorrect",ptrTH1curr->GetName())<<endl; return 0; }

    //TimeCFDmipMEAN = ptrTH1curr->GetMean();
    ptrTH1curr->SetTitle( Form("CFD 1mip plot, mean: %.1f", timeCFDmipMEAN) );
    if(fVerbose >= 4)cout <<Form("CFD 1mip selection: mean %.1f",timeCFDmipMEAN)<<endl;


    Int_t BinsAmplPlot = fQTCPlotMax - fQTCPlotMin;
    Int_t BinsTimePlot = fCFDPlotMax - fCFDPlotMin;

    //TH2 *CFDQTCplot = new TH2F( Form("CFDQTCplot_%i",PMT),"QTCCFD plot",10000,0,10000,200,-100,100);
    TH2 *CFDQTCplot = new TH2F( Form("CFDQTCplot_%i",pmt),"QTCCFD plot",BinsAmplPlot, fQTCPlotMin, fQTCPlotMax,
                                BinsTimePlot, fCFDPlotMin, fCFDPlotMax);

    fTree->Draw( Form("%s - %f:%s + %d >> CFDQTCplot_%i", timeCFDeq.Data(), timeCFDmipMEAN, amplitudeQTCeq.Data(), fPedestalQTC[pmt], pmt ),

                 //Form("%s && %s && %s && %s",SelectionQTCeq.Data(), SelectionCFDeq.Data(), SelectionVertex.Data(), SelectionTRGeq.Data())  );
                 Form("%s && %s && %s",selectionCFDeq.Data(), selectionVertex.Data(), selectionTRGeq.Data())  );


    return CFDQTCplot;

    return 0;
}

/*******************************************************************************************************************/
Int_t AliT0TimeAmplCorr::GetRAWmeanValuesFromOCDB()
{
/*
    if(fVerbose >= 2)cout <<Form("Getting mean values from OCDB for RUN %i ... ", fNumRUN)<<endl;

    if( fNumRUN <= 0 ){ if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::GetRAWmeanValuesFromOCDB::ERROR: RUN number is incorrect")<<endl; return 0; }

    TGrid::Connect("alien://");
    AliCDBManager * CDBManager = AliCDBManager::Instance();
    CDBManager->SetDefaultStorage("raw://");

    CDBManager->SetRun(fNumRUN);
    AliT0Parameters *fParam = AliT0Parameters::Instance();
    fParam->Init();

    AliCDBEntry *CDBEntry = AliCDBManager::Instance()->Get("T0/Calib/TimeDelay");
    AliT0CalibTimeEq *T0CalTimeEq = (AliT0CalibTimeEq*)CDBEntry->GetObject();

    fMeanTVDC = fParam->GetMeanTVDC();
    if(fVerbose >= 4)cout <<Form("Mean TVDC from OCDB: %i", fMeanTVDC)<<endl;

    TString pmtname;
    for (int pmt=0; pmt<24; pmt++){
        pmtname = Form("%s_%02d", (pmt<=11)?"C":"A", (pmt<=11)?pmt+1:pmt-11);  //pmt name for histogram name

        //fMeanCFD[PMT] = fParam->GetCFD(PMT);
        fMeanCFD[pmt] = T0CalTimeEq->GetCFDvalue(pmt,0);
        if(fVerbose >= 4)cout <<Form("Mean CFD %s from OCDB: %i",pmtname.Data(), fMeanCFD[pmt])<<endl;

        fMeanQT1[pmt] = fParam->GetQT1(pmt);
        if(fVerbose >= 4)cout <<Form("Mean QT1 %s from OCDB: %i",pmtname.Data(), fMeanQT1[pmt])<<endl;

        fPedestalQTC[pmt] = T0CalTimeEq->GetCFDvalue(pmt,3);
        if(fVerbose >= 4)cout <<Form("Pedestal QTC %s from OCDB: %i",pmtname.Data(), fPedestalQTC[pmt])<<endl;
    }

    return 1;
*/
}

/*******************************************************************************************************************/
Int_t AliT0TimeAmplCorr::GetRAWmeanValuesByFIT(const Int_t statistic)
{

    //procedure parameters:
    static const Int_t kHistMin = 0;            //minimum for monitor histograms
    static const Int_t kHistMax = 30000;        //maximum for monitor histograms
    static const Int_t kCHInBin = 5;           //bin wight
    static const Int_t kNumBinsInHist = ceil((kHistMax-kHistMin)/kCHInBin); //bins in histograms
    static const Double_t kMipQTCrange = 0.7;


    //temp variables
    TH1 *ptrTH1curr;                                        //pointers to created hist
    Double_t histmean, histsigma;
    Int_t fitMeanCFD, fitMeanQT1, bestMean;
    Double_t fitMeanQTCmip, fitSigmaQTCmip;
    TString TREEpmtname, pmtname;                           //pmt name in form "A_9" or "C_12"
    TString timeCFDeq, amplitudeQT1eq, amplitudeQTCeq;      //variables equation for for TTree->Draw():
    TString selectionCFDeq, selectionQT1eq, selectionQTCmipeq; //selection equation for TTree->Draw():

    //drawing tvdc histogram *****************************************
    if(fVerbose >= 2)cout <<Form("Fitting for mean values TVDC histograms ... ")<<endl;
    TString tvdcVertexEq = "T0_VERTEX_0";
    fMonHists->Add( ptrTH1curr = new TH1I("tvdc", "tvdc", kNumBinsInHist,kHistMin,kHistMax) );
    fTreeChain->Draw( Form("%s>>tvdc", tvdcVertexEq.Data()), Form("0 < %s", tvdcVertexEq.Data()), "", statistic );

    //values range
    if( !FitHistogrammInMaxPeak(ptrTH1curr, histmean, histsigma) ){ if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::GetRAWmeanValuesByFIT::ERROR: TVDC is incorrect")<<endl; return 0; }
    Int_t fitMeanTVDC = floor(histmean);
    if(fVerbose >= 4)cout <<Form("TVDC: mean %i; [%i, %i]",fitMeanTVDC, fitMeanTVDC - frangeTVDC, fitMeanTVDC + frangeTVDC);

    if(fMeanTVDC != 0)
    {
        if(fVerbose >= 4)cout <<Form(" from OCDB: %i, (OCDB-fit) %i", fMeanTVDC, fMeanTVDC-fitMeanTVDC);
        Int_t EntIntegralForFit  = ptrTH1curr->Integral(ptrTH1curr->FindBin(fitMeanTVDC - frangeTVDC), ptrTH1curr->FindBin(fitMeanTVDC + frangeTVDC));
        Int_t EntIntegralForOCDB = ptrTH1curr->Integral(ptrTH1curr->FindBin(fMeanTVDC - frangeTVDC), ptrTH1curr->FindBin(fMeanTVDC + frangeTVDC));
        bestMean = (EntIntegralForOCDB < EntIntegralForFit)? fitMeanTVDC:fMeanTVDC;
        cout <<Form(", chosen %i",bestMean);
    } else bestMean = fitMeanTVDC;

    if(fVerbose >= 4)cout<<endl;
    ptrTH1curr->SetTitle( Form("TVDC: Mean %i, OCDB: %i, chosen %i", fitMeanTVDC, fMeanTVDC, bestMean) );
    fMeanTVDC = bestMean;
    if(ptrTH1curr->GetFunction("fitFunction") != 0)
        ptrTH1curr->GetFunction("fitFunction")->SetRange(fMeanTVDC - frangeTVDC, fMeanTVDC + frangeTVDC); //to see selected range on func at hist
    //tvdc selection mean+-frangeVertex
    TString SelectionTVDCeq(  Form("((%i <= %s)&&(%s <= %i))", fMeanTVDC - frangeTVDC, tvdcVertexEq.Data(), tvdcVertexEq.Data(), fMeanTVDC + frangeTVDC)  );


    for(Int_t PMT=0; PMT<NPMT0; PMT++)
    {
        TREEpmtname = Form("%s_%d", (PMT<=11)?"C":"A", (PMT<=11)?PMT+1:PMT-11); //pmt name for branch
        pmtname = Form("%s_%d", (PMT<=11)?"C":"A", (PMT<=11)?PMT+1:PMT-11);  //pmt name for histogram name

        if(fVerbose >= 2)cout <<Form("Fitting for mean values %s PMT histograms ... ",pmtname.Data())<<endl;

        //leafs names
        timeCFDeq =  Form("(T0_%s_CFD_0)",TREEpmtname.Data());

        /*
        (QT00-QT01) -> (QT10-QT11)
        (QT00-QT01) = (QT10-QT11)*a + b - right

        X   X
        ^   ^
        |   |
        |   signal
        |
        No QTC
        */

        if(fIsNewQTC)
        {
            TString QT00eq = Form("T0_%s_QTC00_0",TREEpmtname.Data());
            TString QT01eq = Form("T0_%s_QTC01_0",TREEpmtname.Data());
            TString QT10eq = Form("T0_%s_QTC10_0",TREEpmtname.Data());
            TString QT11eq = Form("T0_%s_QTC11_0",TREEpmtname.Data());

            TString QT0sel = Form("( (%d < %s)&&(%s < %d) )", frangeNewQTCmin, QT01eq.Data(), QT01eq.Data(), frangeNewQTCmax);
            TString QT1sel = Form("( (%d < %s)&&(%s < %d) )", frangeNewQTCmin, QT11eq.Data(), QT11eq.Data(), frangeNewQTCmax);

            TString QTCnew0eq = Form("(  %s * ((%s - %s)*%f + %f)  )",QT0sel.Data(), QT00eq.Data(), QT01eq.Data(), 1., 0.);
            TString QTCnew1eq = Form("(  %s * ((%s - %s)*%f + %f)  )",QT1sel.Data(), QT10eq.Data(), QT11eq.Data(), fNewQTCa[PMT], fNewQTCb[PMT]);

            amplitudeQT1eq = Form("((%s * %s) + (%s * %s))",QT0sel.Data(), QT01eq.Data(), QT1sel.Data(), QT11eq.Data());
            amplitudeQTCeq = Form("(%s + %s)", QTCnew0eq.Data(), QTCnew1eq.Data());

        }
        else
        {
            amplitudeQT1eq =  Form("(T0_%s_QT1_0)",TREEpmtname.Data());
            amplitudeQTCeq =  Form("(T0_%s_QT0_0 - T0_%s_QT1_0)",TREEpmtname.Data(), TREEpmtname.Data());
        }

        //selection equation:
        selectionCFDeq = Form("(0 < %s)", timeCFDeq.Data());      //0<CFD initial selection, after fitting, additional selection will be added
        selectionQT1eq = Form("(0 < %s)", amplitudeQT1eq.Data()); //0<QT1
        selectionQTCmipeq = Form("(0 < %s)", amplitudeQTCeq.Data()); //0<QTC



        //constructing projections for fMonHists array:

        //CFD: tvdc selection *****************************************
        fMonHists->Add(  ptrTH1curr = new TH1I(Form("%s_CFD",pmtname.Data()),
                                               Form("CFD_0 PMT_%s;channels;events",pmtname.Data()),kNumBinsInHist,kHistMin,kHistMax)  );
        fTreeChain->Draw( Form("%s >> %s_CFD", timeCFDeq.Data(), pmtname.Data()),
                          Form("%s&&%s",selectionCFDeq.Data(),SelectionTVDCeq.Data()), "", statistic  );

        //values range
        if( !FitHistogrammInMaxPeak(ptrTH1curr, histmean, histsigma) ){ if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::GetRAWmeanValuesByFIT::ERROR: %s_CFD is incorrect",pmtname.Data())<<endl; return 0; }
        fitMeanCFD = floor(histmean);


        frangeCFD = 50;
        if(fVerbose >= 4)cout <<Form("CFD: mean %d; range %d; [%d, %d]",fitMeanCFD, frangeCFD, (fitMeanCFD - frangeCFD), (fitMeanCFD + frangeCFD));

        if(fMeanCFD[PMT] != 0)
        {
            if(fVerbose >= 4)cout <<Form(" from OCDB: %i, (OCDB-fit) %i", fMeanCFD[PMT], fMeanCFD[PMT]-fitMeanCFD);
            Int_t entIntegralForFit  = ptrTH1curr->Integral(ptrTH1curr->FindBin(fitMeanCFD - frangeCFD), ptrTH1curr->FindBin(fitMeanCFD + frangeCFD));
            Int_t entIntegralForOCDB = ptrTH1curr->Integral(ptrTH1curr->FindBin(fMeanCFD[PMT] - frangeCFD), ptrTH1curr->FindBin(fMeanCFD[PMT] + frangeCFD));
            bestMean = (entIntegralForOCDB < entIntegralForFit)? fitMeanCFD:fMeanCFD[PMT];
            cout <<Form(", chosen %i",bestMean);
        } else bestMean = fitMeanCFD;

        if(fVerbose >= 4)cout <<""<<endl;
        ptrTH1curr->SetTitle( Form("CFD_0 PMT_%s: Mean %i, OCDB: %i, chosen %i; channels; events",pmtname.Data(), fitMeanCFD, fMeanCFD[PMT], bestMean) );
        fMeanCFD[PMT] = bestMean;
        if(ptrTH1curr->GetFunction("fitFunction") != 0)
            ptrTH1curr->GetFunction("fitFunction")->SetRange(fMeanCFD[PMT] - frangeCFD, fMeanCFD[PMT] + frangeCFD); //to see selected range on function at hist
        //after CFD histogram fitted, additional range cut for CFD: mean +- frangeCFD
        selectionCFDeq += Form("&&((%i <= %s)&&(%s <= %i))", fMeanCFD[PMT] - frangeCFD, timeCFDeq.Data(), timeCFDeq.Data(), fMeanCFD[PMT] + frangeCFD);


        //QT1: CFD  && tvdc selection *****************************************
        fMonHists->Add( ptrTH1curr = new TH1I(Form("%s_QT1",pmtname.Data()),
                                              Form("QT1_0 PMT_%s;channels;events",pmtname.Data()),kNumBinsInHist,kHistMin,kHistMax) );
        fTreeChain->Draw( Form("%s >> %s_QT1", amplitudeQT1eq.Data(), pmtname.Data()),
                          Form("%s&&%s&&%s", selectionQT1eq.Data(),selectionCFDeq.Data(), SelectionTVDCeq.Data()), "", statistic  );

        //values range
        if( !FitHistogrammInMaxPeak(ptrTH1curr, histmean, histsigma) ){ if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::GetRAWmeanValuesByFIT::ERROR: %s_QTC1 %s is incorrect",pmtname.Data(),fIsNewQTC?"new":"old")<<endl; return 0; }
        fitMeanQT1 = floor(histmean);
        if(fVerbose >= 4)cout <<Form("QT1: mean %i; [%i, %i]",fitMeanQT1, fitMeanQT1 - frangeQT1, fitMeanQT1 + frangeQT1);

        if(fMeanQT1[PMT] != 0)
        {
            if(fVerbose >= 4)cout <<Form(" from OCDB: %i, (OCDB-fit) %i", fMeanQT1[PMT], fMeanQT1[PMT]-fitMeanQT1);
            Int_t entIntegralForFit  = ptrTH1curr->Integral(ptrTH1curr->FindBin(fitMeanQT1 - frangeQT1), ptrTH1curr->FindBin(fitMeanQT1 + frangeQT1));
            Int_t entIntegralForOCDB = ptrTH1curr->Integral(ptrTH1curr->FindBin(fMeanQT1[PMT] - frangeQT1), ptrTH1curr->FindBin(fMeanQT1[PMT] + frangeQT1));
            bestMean = (entIntegralForOCDB < entIntegralForFit)? fitMeanQT1:fMeanQT1[PMT];
            cout <<Form(", chosen %i",bestMean);
        } else bestMean = fitMeanQT1;

        if(fVerbose >= 4)cout <<""<<endl;
        ptrTH1curr->SetTitle( Form("QT1_0 PMT_%s: Mean %i, OCDB: %i, chosen %i; channels; events",pmtname.Data(), fitMeanQT1, fMeanQT1[PMT], bestMean) );
        fMeanQT1[PMT] = bestMean;
        if(ptrTH1curr->GetFunction("fitFunction") != 0)
            ptrTH1curr->GetFunction("fitFunction")->SetRange(fMeanQT1[PMT] - frangeQT1, fMeanQT1[PMT] + frangeQT1); //to see selected range on function at histogram

        //after QT1 histogram fitted, additional range cut for QT1 mean +- frangeQT1
        selectionQT1eq += Form("&&((%i <= %s)&&(%s <= %i))", fMeanQT1[PMT] - frangeQT1, amplitudeQT1eq.Data(), amplitudeQT1eq.Data(), fMeanQT1[PMT] + frangeQT1 );


        //QTC: tvdc && QTC && CFD && QT1 selection 1mip for CFD mean ***************
        fMonHists->Add( ptrTH1curr = new TH1I(Form("%s_QTC",pmtname.Data()),
                                              Form("QTC_0 PMT_%s;channels;events",pmtname.Data()),kNumBinsInHist,kHistMin,kHistMax) );
        fTreeChain->Draw( Form("%s >> %s_QTC", amplitudeQTCeq.Data(), pmtname.Data()),
                          Form("%s&&%s&&%s&&%s", selectionQTCmipeq.Data(), selectionQT1eq.Data(),selectionCFDeq.Data(), SelectionTVDCeq.Data()), "", statistic  );

        //values range
        if( !FitHistogrammInMaxPeak(ptrTH1curr, fitMeanQTCmip, fitSigmaQTCmip) ){ if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::GetRAWmeanValuesByFIT::ERROR: %s_QTC %s is incorrect",pmtname.Data(),fIsNewQTC?"new":"old")<<endl; return 0; }
        if(fVerbose >= 4)cout <<Form("QTC 1mip: mean %.1f, sigma %.1f",fitMeanQTCmip, fitSigmaQTCmip)<<endl;

        ptrTH1curr->SetTitle( Form("QTC_0 PMT_%s for 1mip CFD mean: %.1f, sigma %.1f;channels;events",pmtname.Data(),fitMeanQTCmip, fitSigmaQTCmip ) );
        if(ptrTH1curr->GetFunction("fitFunction") != 0)
            ptrTH1curr->GetFunction("fitFunction")->SetRange(fitMeanQTCmip - fitSigmaQTCmip*kMipQTCrange, fitMeanQTCmip + fitSigmaQTCmip*kMipQTCrange); //to see selected range on func at hist

        //after QT1 histogram fitted, additional range cut for QT1 mean +- frangeQT1
        selectionQTCmipeq = Form("((%.1f <= %s)&&(%s <= %.1f))", fitMeanQTCmip - fitSigmaQTCmip*kMipQTCrange, amplitudeQTCeq.Data(), amplitudeQTCeq.Data(), fitMeanQTCmip + fitSigmaQTCmip*kMipQTCrange );

        //CFD: tvdc && QTC && CFD && QT1 selection && 1mip QTC***************
        fMonHists->Add(  ptrTH1curr = new TH1I(Form("%s_CFDmip",pmtname.Data()),
                                               Form("CFD_0 PMT_%s;channels;events",pmtname.Data()),kNumBinsInHist,kHistMin,kHistMax)  );
        fTreeChain->Draw( Form("%s >> %s_CFDmip", timeCFDeq.Data(), pmtname.Data()),
                          Form("%s&&%s&&%s&&%s", selectionQTCmipeq.Data(), selectionQT1eq.Data(),selectionCFDeq.Data(), SelectionTVDCeq.Data()), "", statistic  );

        //values range
        if( !FitHistogrammInMaxPeak(ptrTH1curr, histmean, histsigma) ){ if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::GetRAWmeanValuesByFIT::ERROR: %s_CFD is incorrect",pmtname.Data())<<endl; return 0; }
        fMeanCFDmip[PMT] = floor(histmean);

        if(fVerbose >= 4)cout <<Form("CFD: 1mip mean %i", fMeanCFDmip[PMT])<<endl;

        ptrTH1curr->SetTitle( Form("CFD_0, 1 mip QTC, PMT_%s: Mean %i; channels; events",pmtname.Data(), fMeanCFDmip[PMT]) );

    }//PMT

    return 1;
}

/*******************************************************************************************************************/
TGraph *AliT0TimeAmplCorr::CreateCorrectionGraph(TH1*const sourceprojection, const Int_t kMaxfuncpower, TString fitoption, TString fragoption, const Int_t kNslices)
{

    if(fVerbose >= 2)cout <<Form("Creating correction graph for \"%s\" ... ", sourceprojection->GetName())<<endl;

    //procedure parameters:
    static const Float_t kPointsInUnit = 1.;
    static const Int_t kMinBinsInSlice = 6;
    fitoption += "QNR"; //fit option

    Int_t numTotalEv = sourceprojection->GetEntries();
    if(numTotalEv < kMinBinsInSlice){ if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::CreateCorrectionGraph::ERROR: Too low events in \"%s\" hist: %i",sourceprojection->GetName(), numTotalEv)<<endl; return 0;}

    //histogram range
    Int_t histNBins = sourceprojection->GetNbinsX();
    Float_t histMinX = sourceprojection->GetXaxis()->GetXmin();
    Float_t histMaxX = sourceprojection->GetXaxis()->GetXmax();
    Float_t graphXmin = sourceprojection->GetBinCenter(1);
    Float_t graphXmax = sourceprojection->GetBinCenter(histNBins);

    if(fVerbose >= 4)cout <<Form("PROJECTION RANGE: %.1f, %.1f", histMinX, histMaxX)<<endl;

    //result values
    TGraph *resultGraph = 0;
    Int_t maxPointsInGraph = ceil( (graphXmax-graphXmin + 1)*kPointsInUnit ) + 2;
    Int_t numPointsInGraph = 0;
    Double_t *graphXPoint = new Double_t[maxPointsInGraph];
    Double_t *graphYPoint = new Double_t[maxPointsInGraph];
    Double_t graphYlastPoint = 0;

    //Creating graph without fitting if option "nofit" ********************************************************
    if(fragoption.Contains("nofit")) //searching peaks by maximums and minimums
    {
        if(fVerbose >= 3)cout <<Form("Creating graph without fitting ... ")<<endl;

        Double_t currX, currY;

        for(Int_t currBin = 1; currBin <= histNBins; currBin++)
        {
            currY = sourceprojection->GetBinContent(currBin);

            if(currBin == 1)
                currX = sourceprojection->GetBinLowEdge(currBin);

            else if(currBin == histNBins)
                currX = sourceprojection->GetXaxis()->GetBinUpEdge(currBin);

            else
                currX = sourceprojection->GetBinCenter(currBin);


            if(numPointsInGraph >= maxPointsInGraph) break;

            //if((GraphXmin <= currX)&&(currX <= GraphXmax))
            {
                graphXPoint[numPointsInGraph] = currX;
                graphYPoint[numPointsInGraph] = currY;
                numPointsInGraph++;
            }
        }

        resultGraph = new TGraph(numPointsInGraph, graphXPoint, graphYPoint);

        delete[] graphXPoint;
        delete[] graphYPoint;

        return resultGraph;
    }

    //Searching peaks ******************************************************************************************

    Double_t *peaks;
    Int_t totalPeaks = 0;

    //mean values for one bin
    Double_t totErrIntegral = 0; for(Int_t currBin=1; currBin<=histNBins; currBin++) totErrIntegral += sourceprojection->GetBinError(currBin);
    Double_t meanErrIntegralInSlice = totErrIntegral / (Double_t)kNslices;//mean error for one slice
    Double_t meanWidthSlice = (histMaxX-histMinX) / (Double_t)kNslices;
    Int_t meanBinsInSlice = floor( (Float_t)histNBins / (Double_t)kNslices);
    if(meanBinsInSlice < kMinBinsInSlice)meanBinsInSlice = kMinBinsInSlice;



    if(fragoption.Contains("p")) //searching peaks by maximums and minimums
    {
        if(fVerbose >= 3)cout <<Form("Searching peaks ... ")<<endl;

        TSpectrum *spectrum_max = new TSpectrum();
        //searching maximums and minimums on projection:
        TH1 *tempProjection = (TH1*)sourceprojection->Clone("tempprojection");
        spectrum_max->Search(tempProjection, 0.1, "", 0.01);
        Int_t nmaximums = spectrum_max->GetNPeaks();
        //Double_t *maximums = spectrum_max->GetPositionX();

        tempProjection->Scale(-1);
        TSpectrum *spectrum_min = new TSpectrum();
        spectrum_min->Search(tempProjection, 0.1, "", 0.01);
        Int_t nminimums = spectrum_min->GetNPeaks();
        //Double_t *minimums = spectrum_min->GetPositionX();

        //array of histogram peak, spaces between which will be divided in slices
        Int_t nPeaks = nminimums + nmaximums + 2;
        peaks = new Double_t[nPeaks];
        Double_t *unsortedPeaks = new Double_t[nPeaks];
        Int_t totalUnsortPeaks = 0;

	// don't use intermediate pointer to avoid root5 / root6 incompatibility
        for(Int_t max=0;max<nmaximums;max++)unsortedPeaks[totalUnsortPeaks++] = spectrum_max->GetPositionX()[max];//maximums[max];
        for(Int_t min=0;min<nminimums;min++)unsortedPeaks[totalUnsortPeaks++] = spectrum_min->GetPositionX()[min];//minimums[min];

        delete spectrum_max;
        delete spectrum_min;
        delete tempProjection;

        //formation of the correct sequence of points (summarize arrays and sort them):
        unsortedPeaks[totalUnsortPeaks++] = histMaxX;

        ArrayBubbleSortInc(unsortedPeaks, totalUnsortPeaks);



        //copying sorting array, if bins between peaks less then MinBinInSlice, skip next peak
        peaks[totalPeaks++] = histMinX;
        Int_t fBin, lBin;

        for(Int_t i=0;i<totalUnsortPeaks-1;i++)
        {
            if(TMath::Abs(peaks[totalPeaks-1] - unsortedPeaks[i]) < 0.1) continue; //current peak is the same

            fBin = sourceprojection->FindBin(  peaks[totalPeaks-1]  );
            lBin = sourceprojection->FindBin(  unsortedPeaks[i]  );


            if(histNBins-lBin < kMinBinsInSlice) break; //there are not enough space for last slice

            if(lBin-fBin < kMinBinsInSlice) continue; //current slice is too small

            peaks[totalPeaks++] = unsortedPeaks[i];
        }

        peaks[totalPeaks++] = histMaxX;

        delete[] unsortedPeaks;
    }
    else //case no peaks mode
    {
        totalPeaks = 2;
        peaks = new Double_t[totalPeaks];
        peaks[0] = histMinX;
        peaks[1] = histMaxX;
    }

    if(fVerbose >= 4)
    {
        cout <<Form("PEAKS %i (bin): ", totalPeaks);
        for(Int_t peak=0; peak <= totalPeaks; peak++)
            cout <<Form("%.1f(%i) ",peaks[peak],sourceprojection->FindBin(peaks[peak]));
        cout <<""<<endl;
    }

    //Dividing peaks regions in slices *************************************************************************

    //array of low edges of slices for fitting
    Double_t *slices = new Double_t[kNslices];
    slices[0] = histMinX;
    Int_t numSlices = 1;


    //temp variables
    Double_t errIntegralInSlice, errIntegralLeft; //sum of errors
    Double_t widthSlice, widthSliceLeft;
    Int_t binsInSlice, binsInSliceLeft;
    Int_t firstBin, lastBin, numBins;
    Float_t firstX, lastX;

    // TotalPeaks = 5
    //   0    1        2      3      TotalPeaks-1
    //   O    O        O      O       O   peaks
    //   ******************************   bins
    //   [   ][sli|sli][range][       ]<-------- last bin in last range = (bin-0)
    //       ^                ^
    //       |                |
    //       |             last range begins
    // last bin in         on TotalPeaks-2
    // range = (bin-1)

    if(fVerbose >= 3)cout <<Form("Configuration slices region ... ")<<endl;

    if(fVerbose >= 6) cout <<Form("Mean: err %.2f, wdh %.2f, bin %i; minBin %i", meanErrIntegralInSlice, meanWidthSlice, meanBinsInSlice, kMinBinsInSlice)<<endl;


    for(Int_t peak=0; peak<totalPeaks-1; peak++) //Peaks[peak] = low edge of range
    {
        //number of slices is more than maxnslices
        if(numSlices == kNslices-1) break;

        //range edges
        firstBin = sourceprojection->FindBin(slices[numSlices-1]);
        lastBin = sourceprojection->FindBin(peaks[peak+1]) - 1;
        if(peak == totalPeaks-2)lastBin++; //last bin in last range
        numBins = lastBin-firstBin;
        firstX = sourceprojection->GetXaxis()->GetBinLowEdge(firstBin);
        lastX = sourceprojection->GetXaxis()->GetBinUpEdge(lastBin);

        //values left
        errIntegralLeft = 0;
        for(Int_t currBin=firstBin; currBin<=lastBin; currBin++)
            errIntegralLeft += sourceprojection->GetBinError(currBin);
        widthSliceLeft = lastX-firstX;
        binsInSliceLeft = numBins;

        if(fVerbose >= 6) cout <<Form("Peak %i; Left: err %.2f, wdh %.2f, bin %i",peak, errIntegralLeft, widthSliceLeft, binsInSliceLeft)<<endl;

        //if peak range too small
        if(
                ( (errIntegralLeft < meanErrIntegralInSlice) &&(fragoption.Contains("e")) )||
                ( (widthSliceLeft < meanWidthSlice)          &&(fragoption.Contains("r")) )||
                ( (binsInSliceLeft < meanBinsInSlice)        &&(fragoption.Contains("b")) )||
                (binsInSliceLeft < kMinBinsInSlice)
                ) continue;

        //dividing each range between peaks by mean values
        errIntegralInSlice=widthSlice=binsInSlice=0;
        for(Int_t currBin=firstBin+1; currBin<=lastBin; currBin++)
        {

            widthSlice += sourceprojection->GetBinWidth(currBin);
            widthSliceLeft -= sourceprojection->GetBinWidth(currBin);
            binsInSlice += 1;
            binsInSliceLeft -= 1;
            errIntegralInSlice += sourceprojection->GetBinError(currBin);
            errIntegralLeft -= sourceprojection->GetBinError(currBin);

            if(fVerbose >= 5) cout <<Form("Bin %i, nSlices %i; Left: err %.2f, wdh %.2f, bin %i;",
                                          currBin, numSlices, errIntegralLeft, widthSliceLeft, binsInSliceLeft)<<endl;


            //space for only one slice in array
            if(numSlices+1 == kNslices-1)
            {
                Int_t LastBinInHist = sourceprojection->FindBin(peaks[totalPeaks-1]);
                slices[numSlices++] = sourceprojection->GetXaxis()->GetBinUpEdge(LastBinInHist);
                break;
            }

            //last slice thrift
            if(
                    ( (errIntegralLeft < meanErrIntegralInSlice) &&(fragoption.Contains("e")) )||
                    ( (widthSliceLeft < meanWidthSlice)          &&(fragoption.Contains("r")) )||
                    ( (binsInSliceLeft < meanBinsInSlice)        &&(fragoption.Contains("b")) )||
                    (currBin == lastBin)||(binsInSliceLeft <= kMinBinsInSlice)
                    )
            {
                slices[numSlices++] = lastX;
                //printf("Err: %.1f %.1f, Width: %.1f %.1f, Bin: %i %i\n", ErrIntegralLeft,MeanErrIntegralInSlice,
                //      WidthSliceLeft, MeanWidthSlice, numBins,MeanBinsInSlice);
                break;
            }
            //slice saturation condition
            if(
                    (
                        ( (errIntegralInSlice >= meanErrIntegralInSlice) &&(fragoption.Contains("e")) )||
                        ( (widthSlice >= meanWidthSlice)                 &&(fragoption.Contains("r")) )||
                        ( (binsInSlice >= meanBinsInSlice)               &&(fragoption.Contains("b")) )
                        )
                    &&(kMinBinsInSlice <= binsInSlice)
                    )
            {
                slices[numSlices++] = sourceprojection->GetXaxis()->GetBinLowEdge(currBin);
                errIntegralInSlice = 0;
                widthSlice = 0;
                binsInSlice = 0;
            }

        } // for currBin

    } //for peak

    delete[] peaks;

    if(fVerbose >= 4)
    {
        cout <<Form("SLICES(bin) [%d] : %.1f; ", numSlices, slices[0]);
        for(Int_t currslice=1; currslice < numSlices; currslice++)
            cout <<Form("%.1f(%i); ",slices[currslice],
                        sourceprojection->FindBin(slices[currslice]) - sourceprojection->FindBin(slices[currslice-1]));
        cout <<""<<endl;
    }


    //fitting each slice with different pow polynomial, and choose the best **************************************

    //temp variables
    Double_t X1, X2;
    Double_t currChi, bestChi;
    Double_t pointY;

    TF1 *currFunction = 0, *bestFunction = 0;

    if(fVerbose >= 3) cout <<Form("Searching best polynomials ... ")<<endl;

    Int_t currVerbose = fVerbose;
    for(Int_t currslice=0; currslice < numSlices-1; currslice++)
    {
        X1 = slices[currslice];
        X2 = slices[currslice+1];
        //cout <<Form("%.1f, %.1f",X1,X2)<<endl;

        bestChi = -111;

        for(Int_t currPower=0; currPower<=kMaxfuncpower; currPower++)
        {

            //cout <<Form("CURR POL: %s",Form("(x>=%f && x<%f)*pol%d",X1,X2,currPower))<<endl;
            currFunction = new TF1(Form("currfitFunc_%.1f_%.1f",X1,X2), Form("(x>=%f && x<%f)*pol%d",X1,X2,currPower), X1, X2);
            if(!currFunction)
            {
                if(currVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::CreateCorrectionGraph::ERROR: fit function was NOT created")<<endl;
                continue;
            }
            sourceprojection->Fit(currFunction, fitoption,"",X1,X2);
            //cout <<Form("OK")<<endl;

            currChi = currFunction->GetChisquare();
            if(0.00001 < currChi)
                if((currChi < bestChi)||(bestChi < -100))
                {
                    if(bestFunction) delete bestFunction;
                    bestFunction = currFunction;
                    bestChi = currChi;
                }
                else
                {
                    if(currFunction) delete currFunction;
                }

        }//currPower

        if(currVerbose >= 4) cout <<Form("[%.1f, %.1f] BEST POW: %i, chi%f", X1, X2, bestFunction->GetNpar()-1, bestChi)<<endl;

        //calculatin points for graph

        if(currslice == 0)//expand graph begin
        {
            pointY = bestFunction->Eval(graphXmin);
            graphXPoint[numPointsInGraph] = fQTCPlotMin;
            graphYPoint[numPointsInGraph] = pointY;
            numPointsInGraph++;

        }

        for(Double_t PointX = X1; PointX < X2; PointX += 1./kPointsInUnit)
        {
            pointY = bestFunction->Eval(PointX);

            if(numPointsInGraph >= maxPointsInGraph-1)
            {
                if(currVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::CreateCorrectionGraph::ERROR: wrong number points in graph: max %i; point in graph %i; X %.2f; Max X in Graph %.2f",
                                                 maxPointsInGraph, numPointsInGraph, PointX, graphXmax)<<endl;
                break;
            }

            if((graphXmin <= PointX)&&(PointX <= graphXmax))
            {
                graphXPoint[numPointsInGraph] = PointX;
                graphYPoint[numPointsInGraph] = pointY;
                graphYlastPoint = pointY;
                numPointsInGraph++;
            }
        }

    }  //currFunction
    graphXPoint[numPointsInGraph] = fQTCPlotMax;
    graphYPoint[numPointsInGraph] = graphYlastPoint;
    numPointsInGraph++;

    delete[] slices;

    if(bestFunction) delete bestFunction;
    resultGraph = new TGraph(numPointsInGraph, graphXPoint, graphYPoint);
    resultGraph->GetXaxis()->SetLimits(graphXmin, graphXmax );
    //resultGraph->GetXaxis()->SetLimits(0, 10000 );
    delete[] graphXPoint;
    delete[] graphYPoint;

    return resultGraph;

}

/*******************************************************************************************************************/
TH1 *AliT0TimeAmplCorr::CreatePojection(TH2 *const sourcehist, TString rebinoption, const Float_t kMinbinevX,
                                        const Float_t kMinbinevedge, const Int_t kRebinX, const Int_t kRebinY)
{
    static const Int_t kMinBinEvOnEdge = 5;
    static const Float_t kProjectionErrorK = 0.25;
    static const Int_t kProjectionRebin = 5;
    static const Int_t kProjectionFit = 0;


    if(fVerbose >= 2)cout <<Form("Creating Projection for \"%s\" ... ", sourcehist->GetName())<<endl;

    Int_t numTotalEv = sourcehist->GetEntries();
    if(numTotalEv < 1000)
    {
        cout <<Form("AliT0TimeAmplCorr::CreatePojection::ERROR: Too low events in \"%s\" hist: %i",sourcehist->GetName(), numTotalEv)<<endl;
        return 0;
    }

    TH2 *rebinedhist = sourcehist->Rebin2D(kRebinX, kRebinY, Form("%s_rebin",sourcehist->GetName()) );

    if(rebinoption.Contains("s"))rebinedhist->Smooth(3);

    Int_t numBinsY = rebinedhist->GetNbinsY(); //number of bin along CFD
    Int_t numBinsX = rebinedhist->GetNbinsX(); //number of bin along QTC

    //average events in bin, does not calculated empty bins
    Int_t AverageEvInBinX = 0, bincalculated = 0, integral; //events in bin (QTC axis) on the averages
    for(Int_t xbin = 1; xbin <= numBinsX; xbin++)
    {
        integral = rebinedhist->Integral(xbin,xbin,1,numBinsY);
        if(0 < integral){ AverageEvInBinX += integral; bincalculated++;}
    }
    AverageEvInBinX = ceil((Float_t)AverageEvInBinX / (Float_t)bincalculated);


    if(fVerbose >= 4)cout <<Form("Average events in bin: %i", AverageEvInBinX)<<endl;

    //fist and last bin with minbinevedge*AverageEvInBinX
    Int_t firstBinX = -1, lastBinX = -1;
    for(Int_t xbin = 1; xbin <= numBinsX; xbin++)
    {
        if(firstBinX < 0)
            if( (ceil(kMinbinevedge*AverageEvInBinX) <= rebinedhist->Integral(xbin,xbin,1,numBinsY)) &&
                    ( kMinBinEvOnEdge <= rebinedhist->Integral(xbin,xbin,1,numBinsY)) ) firstBinX = xbin;

        if(lastBinX < 0)
            if( (ceil(kMinbinevedge*AverageEvInBinX) <= rebinedhist->Integral(numBinsX-xbin+1,numBinsX-xbin+1,1,numBinsY)) &&
                    ( kMinBinEvOnEdge <= rebinedhist->Integral(numBinsX-xbin+1,numBinsX-xbin+1,1,numBinsY))) lastBinX = numBinsX-xbin+1;

        if( (lastBinX > 0) && (firstBinX > 0) ) break;
    }

    if(fVerbose >= 4)cout <<Form("FIRST BIN: %i, %.1f (%.1f ev); LAST BIN: %i, %.1f (%.1f ev)",
                                 firstBinX, rebinedhist->GetXaxis()->GetBinCenter(firstBinX), rebinedhist->Integral(firstBinX,firstBinX,1,numBinsY),
                                 lastBinX,  rebinedhist->GetXaxis()->GetBinCenter( lastBinX), rebinedhist->Integral(lastBinX, lastBinX, 1,numBinsY))<<endl;
    if(rebinoption.Contains("r"))if(fVerbose >= 4)cout <<Form("Bin before rebinning: %i", lastBinX-firstBinX)<<endl;

    //finding bins: in each bin minimum minQTCCFDbinEv*QTCCFDnumEvInBinXAverage events
    Float_t *newBinning = new Float_t[lastBinX-firstBinX+2];
    Int_t eventsINbin = 0, numberOfBins = 0;

    if(fVerbose >= 5)cout <<Form("BINNING: ");
    newBinning[0] = rebinedhist->GetXaxis()->GetBinLowEdge(firstBinX);
    for(Int_t xbin = firstBinX+1; xbin <= lastBinX; xbin++)
    {
        eventsINbin += rebinedhist->Integral(xbin,xbin,1,numBinsY);

        if( (eventsINbin >= AverageEvInBinX*kMinbinevX) || (xbin == lastBinX) || (!rebinoption.Contains("r")))
        {
            numberOfBins++;
            newBinning[numberOfBins] = rebinedhist->GetXaxis()->GetBinUpEdge(xbin);
            if(fVerbose >= 5)cout <<Form("%.1f ", newBinning[numberOfBins]);
            eventsINbin = 0;
        }
    } //for sbin
    if(fVerbose >= 4)cout <<""<<endl;
    if(rebinoption.Contains("r"))if(fVerbose >= 4)cout <<Form("Bin after rebining: %i", numberOfBins)<<endl;


    //filling 1D histograms with mean and errors by gauss fit
    TH1 *ptrProjection = new TH1F("Projection","title",numberOfBins,0,1);
    ptrProjection->GetXaxis()->Set(numberOfBins, newBinning);
    ptrProjection->SetOption("E1");

    Int_t firstbin, lastbin;
    Double_t mean, sigma;
    TH1 *ptrProjectionY;
    for(Int_t bin = 0; bin < numberOfBins; bin++)
    {
        mean = sigma = 0.;
        firstbin = rebinedhist->GetXaxis()->FindBin( newBinning[bin] );
        lastbin = rebinedhist->GetXaxis()->FindBin( newBinning[bin+1] );

        ptrProjectionY = rebinedhist->ProjectionY(Form("QTCCFDslice_bin_%i_%i",firstbin,lastbin),firstbin,lastbin);
        if(ptrProjectionY)
        {
            ptrProjectionY->Rebin(kProjectionRebin);
            mean = ptrProjectionY->GetMean();
            sigma = ptrProjectionY->GetRMS();
            if(kProjectionFit) FitHistogrammInRMSRange(ptrProjectionY,mean,sigma, 1.);
            delete ptrProjectionY;
        } else {  cout <<Form("WARING: ProjectionY \"%s\" was NOT created", Form("QTCCFDslice_bin_%i_%i",firstbin,lastbin) )<<endl; }

        //sigma /= Projection->GetBinWidth(bin+1);
        sigma *= kProjectionErrorK;
        ptrProjection->SetBinContent(bin+1,mean);
        ptrProjection->SetBinError(bin+1,sigma);
    }

    delete rebinedhist;
    delete[] newBinning;
    return ptrProjection;
}


/*******************************************************************************************************************/
TH2 *AliT0TimeAmplCorr::CreateCFDQTCplotFromRAW(const Int_t pmt)
{
    //procedure parameters:
    static const Int_t kQTChistMin = 0;            //minimum for QTC axis (for CFD axis range is fMeanCFD[PMT] - frangeCFD)
    static const Int_t kQTChistMax = 20000;        //maximum for QTC axis



    TString TREEpmtname = Form("%s_%d", (pmt<=11)?"C":"A", (pmt<=11)?pmt+1:pmt-11);

    //leafs names
    TString tvdcVertexEq = "T0_VERTEX_0"; //leaf name for vertex
    TString timeCFDeq =  Form("(T0_%s_CFD_0)",TREEpmtname.Data());
    TString amplitudeQT1eq, amplitudeQTCeq;

    /*
      if (qtnew01[pmt]>18000 && qtnew01[pmt]<19000 ) qtcnew = qtnew00[pmt]-qtnew01[pmt];

      if (qtnew11[pmt]>18000 && qtnew11[pmt]<19000 ) qtcnew = a[pmt]*(qtnew10[pmt]-qtnew11[pmt]) + b[pmt];

      if (qtcnew>0) {

    (QT00-QT01) -> (QT10-QT11)
    (QT00-QT01)*a + b = (QT10-QT11)
    */

    if(fIsNewQTC)
    {
        TString QT00eq = Form("(T0_%s_QTC00_0)",TREEpmtname.Data());
        TString QT01eq = Form("(T0_%s_QTC01_0)",TREEpmtname.Data());
        TString QT10eq = Form("(T0_%s_QTC10_0)",TREEpmtname.Data());
        TString QT11eq = Form("(T0_%s_QTC11_0)",TREEpmtname.Data());

        TString QT0sel = Form("( (%d < %s)&&(%s < %d) )", frangeNewQTCmin, QT01eq.Data(), QT01eq.Data(), frangeNewQTCmax);
        TString QT1sel = Form("( (%d < %s)&&(%s < %d) )", frangeNewQTCmin, QT11eq.Data(), QT11eq.Data(), frangeNewQTCmax);

        TString QTCnew0eq = Form("(  %s * ((%s - %s)*%f + %f)  )",QT0sel.Data(), QT00eq.Data(), QT01eq.Data(), fNewQTCa[pmt], fNewQTCb[pmt]);
        TString QTCnew1eq = Form("(  %s * ((%s - %s)*%f + %f)  )",QT1sel.Data(), QT10eq.Data(), QT11eq.Data(), fNewQTCa[pmt], fNewQTCb[pmt]);

        amplitudeQT1eq = Form("((%s * %s) + (%s * %s))",QT0sel.Data(), QT01eq.Data(), QT1sel.Data(), QT11eq.Data());
        amplitudeQTCeq = Form("(%s + %s)",QTCnew0eq.Data(), QTCnew1eq.Data());

    }
    else
    {
        amplitudeQT1eq =  Form("(T0_%s_QT1_0)",TREEpmtname.Data());
        amplitudeQTCeq =  Form("(T0_%s_QT0_0 - T0_%s_QT1_0)",TREEpmtname.Data(), TREEpmtname.Data());
    }


    //values ranges
    Int_t cfdMin = fMeanCFD[pmt] - frangeCFD;
    Int_t cfdMax = fMeanCFD[pmt] + frangeCFD;

    //selection equation:
    TString selectionTVDCeq = Form("((%i <= %s)&&(%s <= %i))", fMeanTVDC - frangeTVDC, tvdcVertexEq.Data(),
                                   tvdcVertexEq.Data(), fMeanTVDC + frangeTVDC );
    TString selectionCFDeq = Form("((%i <= %s)&&(%s <= %i))", cfdMin, timeCFDeq.Data(), timeCFDeq.Data(), cfdMax );
    TString selectionQT1eq = Form("((%i <= %s)&&(%s <= %i))", fMeanQT1[pmt] - frangeQT1, amplitudeQT1eq.Data(),
                                  amplitudeQT1eq.Data(), fMeanQT1[pmt] + frangeQT1 );
    TString commonSelection = Form("%s&&%s&&%s", selectionTVDCeq.Data(), selectionCFDeq.Data(), selectionQT1eq.Data());



    //QTC vs CFD:  tvdc && CFD && QT1 selection
    TH2 *CFDQTCplot = new TH2I(Form("CFDQTCplot_%i",pmt),"",kQTChistMax - kQTChistMin, kQTChistMin,kQTChistMax,cfdMax-cfdMin,cfdMin-fMeanCFD[pmt],cfdMax-fMeanCFD[pmt]);
    fTreeChain->Draw( Form("(%s - %i) : %s >> CFDQTCplot_%i", timeCFDeq.Data(), fMeanCFDmip[pmt], amplitudeQTCeq.Data(), pmt), commonSelection.Data() );

    return CFDQTCplot;
}

/*******************************************************************************************************************/
TH1 *AliT0TimeAmplCorr::CreateGraphAssay(TH2 *const sourceplot, Double_t &chisquare, const TGraph *sourcegraph)
{

    if(sourcegraph) if(fVerbose >= 2)cout <<Form("Assaying graph for \"%s\":", sourceplot->GetName())<<endl;
    if(!sourcegraph) if(fVerbose >= 2)cout <<Form("Assaying original plot \"%s\":", sourceplot->GetName())<<endl;

    //procedure parameters:
    static const Int_t kUseFuncRange = 0; //in case graph exist, use its range, graph expanded on (0, 10000), so no region by graph

    if(fVerbose >= 4)cout <<Form("Creating profile ... ")<<endl;

    TProfile *tempprofile = sourceplot->ProfileX("tempprofile", 1, -1, "o");
    if(!tempprofile) {if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::CreateGraphAssay::ERROR: profile for \"%s\" does NOT created ", sourceplot->GetName())<<endl; return 0;}

    if(fVerbose >= 4)cout <<Form("Creating X projection ... ")<<endl;
    TH1* projection = tempprofile->ProjectionX(Form("%s_assay",sourceplot->GetName()), "E");
    delete tempprofile;

    Double_t Xmin = projection->GetXaxis()->GetXmin();
    Double_t Xmax = projection->GetXaxis()->GetXmax();


    if(sourcegraph)
    {
        SubtractGraphFromHist(projection,sourcegraph);
        if(kUseFuncRange){
            Xmin = sourcegraph->GetXaxis()->GetXmin();
            Xmax = sourcegraph->GetXaxis()->GetXmax();
        }
    }

    if(fVerbose >= 4)cout <<Form("Fitting projection ... ")<<endl;

    TF1 *linearfunc = new TF1("linearfunc","pol0", Xmin, Xmax);

    projection->Fit(linearfunc, "Q", "", Xmin, Xmax);
    chisquare = linearfunc->GetChisquare();

    /***/
    //calculating no normed chi-square for linear function taking in account bin error of projection
    chisquare = 0.;
    Double_t binCntr, binVal, binErr, funcVal;
    for(Int_t projBin = 1; projBin <= projection->GetNbinsX(); projBin++)
    {
        binCntr = projection->GetBinCenter(projBin);
        binVal = projection->GetBinContent(projBin);
        binErr = projection->GetBinError(projBin);
        funcVal = linearfunc->Eval(binCntr);

        if(0.01 < binErr)
            chisquare += (funcVal-binVal)*(funcVal-binVal)/binErr/binErr;

    }
    /***/

    if(linearfunc) delete linearfunc;
    return projection;
}

/*******************************************************************************************************************/
Double_t AliT0TimeAmplCorr::TryProjFitOption(TH2 *sourceplot, TH1 *sourceprojection, TH1 **bestprojection, TGraph **bestgraph, TH1 **bestassayhist, TH1 **CHItrend,
                                             Int_t channel, TString projectionoption, TString fragoption, Int_t nslices, Int_t rebinX, Float_t minbinevX)
{
    //not set parameters
    static const Int_t kMaxFuncPower = 9;
    TString fitoption = "MOI";

    static const Int_t kMaxAttempts = 100000;

    static Double_t bestChisquare;
    static Int_t lastpmt = -1;
    static Int_t nattempt = 0;

    if(channel != lastpmt)
    {
        bestChisquare = -111;
        lastpmt = channel;
        nattempt = 0;
    };

    nattempt++;
    if(kMaxAttempts<nattempt)if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::TryProjFitOption::ERROR: too many points for trend: %i", nattempt)<<endl;

    //TH1 *currProjection;
    TGraph *currGraph;
    TH1 *currAssayHist;

    Double_t currChisquare;
    TString parString = Form("[N%i, proj:\"%s\" frag:\"%s\" nsl:%i, minEvInBin:%.2f, rX%i]",
                             nattempt, projectionoption.Data(), fragoption.Data(), nslices, minbinevX, rebinX);

    if(fVerbose >= 2)cout <<Form("\nTrying parameters: %s", parString.Data())<<endl;

    currGraph = CreateCorrectionGraph(sourceprojection, kMaxFuncPower, fitoption, fragoption, nslices);
    if(!currGraph){ if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::TryProjFitOption::ERROR: Graph absent")<<endl; return -111.;}

    currAssayHist = CreateGraphAssay(sourceplot, currChisquare, currGraph);
    if(!currAssayHist){ if(fVerbose >= 1) cout <<Form("AliT0TimeAmplCorr::TryProjFitOption::ERROR: Assay Histogram absent")<<endl; delete currGraph; return -111.;}

    //adding point to trend
    if(CHItrend)
    {
        if( !(*CHItrend) ) *CHItrend = new TH1D(Form("CHItrend_CH%d",channel),Form("Chi-square trend for CH%d",channel), kMaxAttempts,1,kMaxAttempts);
        if( nattempt <= kMaxAttempts ){
            (*CHItrend)->Fill(parString.Data(), currChisquare);
            (*CHItrend)->SetAxisRange(1,nattempt,"X");
        }
    }

    if(fVerbose >= 2)cout <<Form("RESULT CHISQUARE: %.1f",currChisquare)<<endl;

    if(  (currChisquare < bestChisquare)||(bestChisquare < 0)  )
    {
        bestChisquare = currChisquare;

        if(*bestprojection) delete *bestprojection;
        if(*bestgraph) delete *bestgraph;
        if(*bestassayhist) delete *bestassayhist;

        currGraph->SetName( Form("corgraph_CH%d", channel) );
        currGraph->SetTitle( Form("Correction graph CH%d", channel) );

        currAssayHist->SetName( Form("func_assay_wif_CH%d", channel) );
        currAssayHist->SetTitle( Form("Graph assay CH%d", channel) );

        *bestprojection = (TH1*)(  sourceprojection->Clone( Form("projection_CH%d", channel) )  );
        (*bestprojection)->SetTitle( Form("Time vs Amplitude CH%d;Amplitude;Time",channel) );

        *bestgraph = currGraph;
        *bestassayhist = currAssayHist;
    }
    else
    {
        delete currGraph;
        delete currAssayHist;
    }

    return bestChisquare;

}

/*******************************************************************************************************************/
/*!
 * \brief Make QTCCFD plots from tree in connected files
 * existing plots will be replaced, also monitoring histograms CFD, QT1, QTC, CFD_for_1_mip will be created
 * for check. Range of fit function in histograms show range of selection for creating QTCCFD plots.
 * Two procedures for mean values may be done: get values from ocdb and find values by fitting histograms
 * \param option
 * "noocdb" - do not use OCD base
 * "nomeanfit" - do not make fitting for mean values
 * "fastmeanscan" - use a part of statistic
 * \return 1 if successfully else 0
 */
Int_t AliT0TimeAmplCorr::MakeQTCCFDplots(TString option)
{
    if(  ((fIsRAWtree)? fTreeChain:fTree) == 0)
    {
        if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::MakeQTCCFDplots::ERROR: Data tree was not created")<<endl;
        return 0;
    }

    Int_t nEntries = (fIsRAWtree)? fTreeChain->GetEntries():fTree->GetEntries();

    if(nEntries < 10000)
    {
        if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::MakeQTCCFDplots::ERROR: too low events in tree: %ld", (long int)fTreeChain->GetEntries())<<endl;
        return 0;
    }

    ResetArray( &fMonHists, "recreate" );
    ResetMeanValues();

    if( !(fIsNewQTC && !fIsRAWtree) ){ //there is no need mean values for new QTC from ESD file

        if(fVerbose >= 2)cout <<Form("\n\n======================== SEARCH AVERAGES ========================")<<endl;

        Int_t meanCalculated = 0;
        if( !option.Contains("noocdb") )
        {
            meanCalculated = GetRAWmeanValuesFromOCDB();
            if( !meanCalculated ) if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::MakeQTCCFDplots::ERROR: mean values does not taken from OCDB")<<endl;
        }

        if(fIsRAWtree && !option.Contains("nomeanfit")){
            meanCalculated = GetRAWmeanValuesByFIT( option.Contains("fastmeanscan")?100000:90000000 );
            if( !meanCalculated ) if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::MakeQTCCFDplots::ERROR: mean values does not calculated")<<endl;
        }

        if( !meanCalculated && fIsRAWtree){ if(fVerbose >= 1)cout <<Form("AliT0TimeAmplCorr::MakeQTCCFDplots::ERROR: mean values does not setted")<<endl; return 0; }

    }

    if(fVerbose >= 2)cout <<Form("\n\n===================== CREATING CFD-QTC PLOTS ====================")<<endl;

    ResetArray( &fQTCCFDplot, "recreate" );

    TString pmtname;
    TH2 *currQTCCFDplot = 0;

    for(Int_t PMT=0; PMT<NPMT0; PMT++)
    {
        pmtname = Form("%s_%02d", (PMT<=11)?"C":"A", (PMT<=11)?PMT+1:PMT-11);

        if(fVerbose >= 2)cout <<Form("Making 2D CFD-QTC plots for %s with %s", pmtname.Data(), fIsRAWtree?"RAW":"ESD")<<endl;

        currQTCCFDplot = 0;
        currQTCCFDplot = (fIsRAWtree)? CreateCFDQTCplotFromRAW(PMT) : CreateCFDQTCplotFromESD(PMT);

        if(currQTCCFDplot)
        {
            currQTCCFDplot->SetName( Form("%s_QTCCFD", pmtname.Data()) );
            currQTCCFDplot->SetTitle( Form("%sQTC vs CFD %s %s, ch;QTC;CFD",fIsNewQTC?"new":"old",fIsRAWtree?"RAW":"ESD", pmtname.Data()) );
            fQTCCFDplot->Add(currQTCCFDplot);
        }

    }

    fIsPlotsLoaded = 0;
    return fQTCCFDplot->GetEntries();
}

/*******************************************************************************************************************/
/*!
 * \brief Make projection TGraph function for each created or merged plot
 * To create a correction function a 1D histograms is prepared firstly for plot with CreatePojection():\n
 *  1)QTCCFD plot smoothed if rebinoption contains "s"\n
 *  2)QTCCFD plot rebined with rebinX and rebinY values, and mean number of events in bin AverageEvInBinX is founded.\n
 *  3)Bounds bins are determined with AverageEvInBinX*minbinevedge <= (number of events in bound bin).\n
 *  4)if rebinoption contains "r", QTCCFD plot rebinned along X (amplitude) on condition\n
 *    AverageEvInBinX*minbinevX <= (number of events in new bin)\n
 *  5)1D histogram created by mean and RMS values for Y projection for each new bin\n
 * \n
 * After 1D histogram created, it divided by slices and fitted with CreateCorrectionGraph():\n
 *  1)if fragoption contains "p" maximum and minimum peaks found for 1D histogramm\n
 *  2)ranges between peaks or between first and last bins divided in slices so that the weight of each slice\n
 *    corresponded to (total weight of histogram) / maxnslices. Weight may be varied by fragoption:\n
 *    "e"-by errors; "r"-by ranges; "b"-by bins; all containing option taking in account\n
 *  3)each slice fitted by polynomials with various power, and chosen the one with the best chi-square\n
 *  4)if fragoption contains "nofit", correction function plotted by points in hist, no other procedures are done\n
 *  5)Correction function decremented from QTCCFD plot and chi-square by linear function\n
 *    is calculated for assay the correction function with CreateGraphAssay().\n
 * \n
 * QTCCFD plot are much differ from each other and there are not unique option for all, for this reason\n
 * function MakeATCorrectionForPlots(TString configfilename = "") takes range for main options from file\n
 * <configfilename> and by iterating over all options chooses the best by assay. If file is not defended
 * or is not opened, parameters ranges set by default.\n
 * Configuration file example:\n
 * \n
 * [SLICES]\n
 * 5   80	3\n
 * [REBINX]\n
 * 5   40	3\n
 * [EVINBIN]\n
 * 1    10	3\n
 * # Option for projection algorithm:\n
 * # "n" all projection bins are same and equal rebinX\n
 * # "r" - rebining projection with various bin values with evinbin in each\n
 * # "s" - smoothing plot before making projection\n
 * [REBINOPT]\n
 * r rs s\n
 * [FRAGOPT]\n
 * e r b pe pr bp perb\n
 * [PMTS]\n
 * 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23\n
 * # 0	9 10 16\n
 * \n
 * there # - comments, parameter names are placed in []:\n
 * SLICES - maxnslices; REBINX - rebinX; EVINBIN - minbinevX; REBINOPT - rebinoption; FRAGOPT - fragoption; PMTS - pmts numbers\n
 * after SLICES, REBINX, EVINBIN  parameters follow its minimum, maximum values and number of iterations\n
 * after REBINOPT, FRAGOPT, PMTS follow parameters enumeration. Values divided by space or tab\n
 * \n
 * Resulting correction function may be taken with TObjArray* GetGraphs() or wrote in file with SaveResultsInFile();\n
 * While MakeATCorrectionForPlots as well as correction function, created 1D projection histogram, function assay histograms\n
 * and chi-square trends to monitor the chi-square dependence from iterated options, this information also may be take by\n
 * "get" functions or wrote in file with  common function SaveResultsInFile();
 * \n
 * EXAMPLE:\n
 *
 * T0ampltimeCorrection->MergeQTCCFDplots( "/Plots_mip/ATC_Plots_236349.root" );
 * T0ampltimeCorrection->MergeQTCCFDplots( "/Plots_mip/ATC_Plots_236337.root" );
 * T0ampltimeCorrection->MergeQTCCFDplots( "/Plots_mip/ATC_Plots_236246.root" );
 *
 * if( T0ampltimeCorrection->MakeATCorrectionForPlots("./AT_fast_scan.corr") )
 *   T0ampltimeCorrection->SaveResultsInFile( "comment" );
 *
 * \param configfilename
 * \return
 */
Int_t AliT0TimeAmplCorr::MakeATCorrectionForPlots(TString configfilename)
{
    if(fQTCCFDplot == 0){ if(fVerbose >= 2) cout << Form("AliT0TimeAmplCorr::MakeATCorrectionForPlots::ERROR: QTCCFD plots array does NOT exist, Plots must be loaded or costructed")<<endl; return 0;}

    Int_t numPlots = fQTCCFDplot->GetEntries();

    if(numPlots <= 0){ if(fVerbose >= 2) cout << Form("AliT0TimeAmplCorr::MakeATCorrectionForPlots::ERROR: QTCCFD plots does NOT exist")<<endl; return 0;}
    else if(fVerbose >= 2) cout << Form("\n\n%i QTCCFD plots available for Making Time Amplitude Correction  ... ", numPlots)<<endl;


    //PMT correction array:
    Int_t pmtArray[NPMT0], numPMTtoCor = 0;

    //sorting out parameters; min, max, iterator, best parameter
    Int_t slicesMin = 0, slicesMax = 0, slicesI = 0, bestslices = 0;
    Int_t rebinXMin = 0, rebinXMax = 0, rebinXI = 0, bestrebinX = 0;
    Float_t evInBinMin = 0., evInBinMax = 0., evInBinI = 0., bestEvInBin = 0.;

    TObjArray *rebinOptArr = 0, *fragOptArr = 0;
    TString rebinopt, fragopt, bestrebinopt, bestfragopt;

    ResetArray( &rebinOptArr, "recreate" );
    ResetArray( &fragOptArr, "recreate" );

    //not set parameters
    static const Float_t kMinBinEvEdge = 0.001;
    static const Int_t kRebinY = 5;

    //reading values range from file
    FILE *configfile = 0;
    if(configfilename != "") configfile = fopen(configfilename.Data(),"rt");

    if(configfile != 0 )
    {
        if(fVerbose >= 2) cout << Form("Configuration file \"%s\" was opened for values range",configfilename.Data());

        char buffer[MAX_CHAR_LINE], *fgetsval; Int_t read;
        Int_t stateReading = 0, stateTableReading = 0;

        while(!feof(configfile)) {
            // read a word from the files and skip comments
            read = fscanf(configfile, "%s", buffer);
            if( !read || (read == EOF) || !strlen(buffer)) continue;
            if(buffer[0] == '#') {fgetsval = fgets(buffer, MAX_CHAR_LINE, configfile); continue; }

            if (buffer[0] == '['){

                if (strstr(buffer, "SLICES")) {
                    fgetsval = fgets(buffer , MAX_CHAR_LINE, configfile);
                    stateReading = 1;
                    stateTableReading = 1;
                    continue;}

                if (strstr(buffer, "REBINX")) {
                    fgetsval = fgets(buffer , MAX_CHAR_LINE, configfile);
                    stateReading = 2;
                    stateTableReading = 1;
                    continue;}

                if (strstr(buffer, "EVINBIN")) {
                    fgetsval = fgets(buffer , MAX_CHAR_LINE, configfile);
                    stateReading = 3;
                    stateTableReading = 1;
                    continue;}

                if (strstr(buffer, "REBINOPT")) {
                    fgetsval = fgets(buffer , MAX_CHAR_LINE, configfile);
                    stateReading = 4;
                    continue;}

                if (strstr(buffer, "FRAGOPT")) {
                    fgetsval = fgets(buffer , MAX_CHAR_LINE, configfile);
                    stateReading = 5;
                    continue;}

                if (strstr(buffer, "PMTS")) {
                    fgetsval = fgets(buffer , MAX_CHAR_LINE, configfile);
                    stateReading = 6;
                    continue;}
            } //if '['

            Int_t value;
            Float_t valuef;
            switch (stateReading) {
            //SLICES
            case 1:
                value = atoi(buffer);

                switch (stateTableReading) {
                case 1: slicesMin = value; break;
                case 2: slicesMax = value; break;
                case 3:
                    if((value < 2)||(slicesMax == slicesMin)) slicesI = (slicesMax - slicesMin + 1)*2; //only min value
                    else slicesI = (Int_t) floor( (slicesMax - slicesMin)/(Int_t)value );
                    if(slicesI < 1)slicesI = 1;
                    break;
                }

                stateTableReading++;

                break;

                //REBINX
            case 2:
                value = atoi(buffer);

                switch (stateTableReading) {
                case 1: rebinXMin = value; break;
                case 2: rebinXMax = value; break;
                case 3:
                    if((value < 2)||(rebinXMax == rebinXMin)) rebinXI = (rebinXMax - rebinXMin + 1)*2; //only min value
                    else rebinXI = (Int_t) floor( (rebinXMax - rebinXMin)/(Int_t)value );
                    if(slicesI < 1)slicesI = 1;
                    break;
                }

                stateTableReading++;

                break;

                //EVINBIN
            case 3:
                valuef = atof(buffer);

                switch (stateTableReading) {
                case 1: evInBinMin = valuef; break;
                case 2: evInBinMax = valuef; break;
                case 3:
                    if((valuef < 2.)||(evInBinMax == evInBinMin)) evInBinI = (evInBinMax - evInBinMin + 1)*2; //only min value
                    else evInBinI =  (evInBinMax - evInBinMin)/floor(valuef);
                    break;
                }

                stateTableReading++;

                break;

                //REBINOPT
            case 4:
                rebinOptArr->Add( new TObjString( buffer ) );
                break;

                //FRAGOPT
            case 5:
                fragOptArr->Add( new TObjString( buffer ) );
                break;

                //PMTS
            case 6:
                value = atoi(buffer);
                if(value < 0) for(numPMTtoCor=0;numPMTtoCor<NPMT0;numPMTtoCor++)pmtArray[numPMTtoCor] = numPMTtoCor;
                if(numPMTtoCor < NPMT0) pmtArray[numPMTtoCor++] = value;
                break;


            }


        }//while

        fclose(configfile);

    }//(configfile != 0 )
    else
    {
        if(fVerbose >= 2)
            if(configfilename == "") cout << Form("Calculating for default values range")<<endl;
            else cout << Form("Configuration file \"%s\" was NOT opened, calculating for default values range",configfilename.Data())<<endl;

        slicesMin = 20; slicesMax = 80; slicesI = 40;
        rebinXMin = 10; rebinXMax = 20; rebinXI = 10;
        evInBinMin = .01; evInBinMax = .1; evInBinI = .02;

        rebinOptArr->Add( new TObjString("") );
        rebinOptArr->Add( new TObjString("r") );
        rebinOptArr->Add( new TObjString("s") );
        rebinOptArr->Add( new TObjString("sr") );

        fragOptArr->Add( new TObjString("nofit") );
        fragOptArr->Add( new TObjString("e") );
        fragOptArr->Add( new TObjString("pe") );

        for(numPMTtoCor=0;numPMTtoCor<NPMT0;numPMTtoCor++)pmtArray[numPMTtoCor] = numPMTtoCor;
    }


    if(fVerbose >= 2)
    {
        cout << Form("Values range:")<<endl;
        cout << Form("Slice:   min %5i, max %5i, iter %5i",slicesMin, slicesMax, slicesI)<<endl;
        cout << Form("RebinX   min %5i, max %5i, iter %5i",rebinXMin, rebinXMax, rebinXI)<<endl;
        cout << Form("EvInBin: min %.3f, max %.3f, iter %.3f",evInBinMin, evInBinMax, evInBinI)<<endl;
        cout << Form("Rebin option (%i):",rebinOptArr->GetLast()+1); for(Int_t i=0; i<= rebinOptArr->GetLast();i++)
            cout<<" \""<<((TObjString*) rebinOptArr->At(i))->GetString().Data()<<"\""; cout<<endl;
        cout <<Form("Fragmentation option (%i):",fragOptArr->GetLast()+1); for(Int_t i=0; i<= fragOptArr->GetLast();i++)
            cout<<" \""<<((TObjString*) fragOptArr->At(i))->GetString().Data()<<"\""; cout<<endl;
        cout<<Form("PMTs (%i):",numPMTtoCor); for(Int_t i=0;i<numPMTtoCor;i++)cout<<" "<<pmtArray[i]<<"";cout<<endl;
    }


    Int_t nCombination = 0, currCombination = 0;

    for(Int_t rebinoptN = 0; rebinoptN <= rebinOptArr->GetLast();  rebinoptN++)
    {
        for(Int_t rebinX =  rebinXMin;    rebinX <=    rebinXMax;     rebinX += rebinXI)
            for(Float_t EvInBin =   evInBinMin;   EvInBin <=   evInBinMax;    EvInBin += evInBinI)
            {
                // for this rebin option EvInBin does not matter, so only one (max) value variation
                if( !( ((TObjString*) rebinOptArr->At(rebinoptN))->GetString().Contains("r") ) ) EvInBin = evInBinMax;

                for(Int_t fragoptN =  0;  fragoptN <=  fragOptArr->GetLast();   fragoptN++)
                {
                    for(Int_t slices =    slicesMin;    slices <=    slicesMax;     slices += slicesI)
                    {

                        nCombination++;

                    }
                }

            }
    }



    if(fVerbose >= 2) cout << Form("\nFitting iteration for one PMT: %i, n PMT: %i, Total: %i",nCombination, numPMTtoCor, nCombination*numPMTtoCor)<<endl;
    nCombination *= numPMTtoCor;

    ResetArray( &fProjection, "recreate" );
    ResetArray( &fGraphs, "recreate" );
    ResetArray( &fChiTrends, "recreate" );
    ResetArray( &fAssay, "recreate" );


    //temp variables
    TString pmtname;
    TH2 *currQTCCFDplot = 0;
    TH1 *currTrend = 0;
    TH1 *currProjection = 0;
    TH1 *bestProjection = 0;
    TGraph *bestGraph = 0;
    TH1 *bestAssayHist = 0;
    TH1 *origAssayHist = 0;
    Double_t bestChisquare, thebestChisquare = -111.;

    fChiTrends->Add( new TH1I("stat_nslices","stat of nslices",100,1,100) );
    fChiTrends->Add( new TH1I("stat_rebinX","stat of rebinX",100,1,100) );
    fChiTrends->Add( new TH1F("stat_EvInBin","stat of EvInBin",100,1,100) );
    fChiTrends->Add( new TH1I("stat_rebinopt","stat of rebinopt",100,1,100) );
    fChiTrends->Add( new TH1I("stat_fragopt","stat of fragopt",100,1,100) );

    Double_t timestart = fTimer->RealTime();
    fTimer->Continue();

    Int_t pmt;
    for(Int_t iPMT = 0; iPMT < numPMTtoCor; iPMT++)
    {
        pmt = pmtArray[iPMT];

        if((pmt < 0)||(23 < pmt))
        {
            if(fVerbose >= 2) cout << Form("WARING: Wrong PMT number :%i",pmt) << endl;
            continue;
        }

        pmtname = Form("%s_%02d", (pmt<=11)?"C":"A", (pmt<=11)?pmt+1:pmt-11);

        if(fVerbose >= 2) cout << Form("\n======================== WORKING ON %s ========================",pmtname.Data())<<endl;


        currQTCCFDplot = (TH2*)( fQTCCFDplot->FindObject( Form("%s_QTCCFD", pmtname.Data()) ) );

        if(!currQTCCFDplot)
        {
            if(fVerbose >= 2) cout << Form("WARING: There is no QTCCFD plot for \"%s\"",pmtname.Data())<<endl;
            continue;
        }


        //assaying original plot
        origAssayHist = 0;
        origAssayHist = CreateGraphAssay(currQTCCFDplot,bestChisquare); thebestChisquare = bestChisquare;
        if(origAssayHist)
        {
            fAssay->Add(origAssayHist);
            origAssayHist->SetName( Form("%s_func_assay_wof", pmtname.Data()) );
            origAssayHist->SetTitle( Form("%s Graph assay %s", pmtname.Data(),
                                          origAssayHist->GetTitle()) );
        }

        bestProjection = 0;
        bestGraph = 0;
        bestAssayHist = 0;
        currTrend = 0;

        for(Int_t rebinoptN = 0; rebinoptN <= rebinOptArr->GetLast();  rebinoptN++)
        {

            rebinopt = ((TObjString*) rebinOptArr->At(rebinoptN))->GetString();

            for(Int_t rebinX =  rebinXMin;    rebinX <=    rebinXMax;     rebinX += rebinXI){
                for(Float_t EvInBin =   evInBinMin;   EvInBin <=   evInBinMax;    EvInBin += evInBinI)
                {
                    // for this rebin option EvInBin does not matter, so only one (max) value variation
                    if(  !rebinopt.Contains("r") ) EvInBin = evInBinMax;

                    currProjection = CreatePojection(currQTCCFDplot,rebinopt, EvInBin, kMinBinEvEdge, rebinX, kRebinY);
                    if(!currProjection)
                    {
                        if(fVerbose >= 1)
                            cout <<Form("AliT0TimeAmplCorr::MakeATCorrectionForPlots::ERROR: Projection absent")<<endl;
                        continue;
                    }

                    for(Int_t fragoptN =  0;  fragoptN <=  fragOptArr->GetLast();   fragoptN++)
                    {
                        fragopt = ((TObjString*) fragOptArr->At(fragoptN))->GetString();

                        for(Int_t slices =    slicesMin;    slices <=    slicesMax;     slices += slicesI)
                        {

                            if(fVerbose >= 2){
                                currCombination++;
                                Double_t timepass = fTimer->RealTime() - timestart;
                                Double_t percents = currCombination/ (Float_t)(nCombination);
                                Double_t timeestimate = (percents < 0.1)?0:timepass/percents*( 1. - percents);
                                Int_t hoursestimate = (Int_t)floor(timeestimate/3600.);
                                Int_t minutesestimate = (Int_t)floor( (timeestimate-hoursestimate*3600.)/60. );

                                cout <<endl<<endl; fTimer->Print(); cout << Form("Iteration %i/%i (%.1f), Estimate time left: %02ih %02im", currCombination, nCombination, percents*100, hoursestimate, minutesestimate) <<endl;
                                fTimer->Continue();
                            }

                            bestChisquare = TryProjFitOption(currQTCCFDplot, currProjection, &bestProjection, &bestGraph, &bestAssayHist, &currTrend,
                                                             pmt,rebinopt, fragopt, slices,rebinX,EvInBin);

                            if(bestChisquare < 0.){ if(fVerbose >= 1)cout <<Form("WARING: current calibration missed"); continue; }

                            if(bestProjection){
                                bestProjection->SetName( Form("%s_QTCCFDprojection", pmtname.Data()) );
                                bestProjection->SetTitle( Form("QTC vs CFD %s Projection, ch;QTC;CFD", pmtname.Data()) );}

                            if(bestGraph){
                                bestGraph->SetName( Form("%s_QTCCFDgraph", pmtname.Data()) );
                                bestGraph->SetTitle( Form("%s Correction graph Chi = %.1f", pmtname.Data(), bestChisquare) );}

                            if(bestAssayHist){
                                bestAssayHist->SetName( Form("%s_func_assay_wif", pmtname.Data()) );
                                bestAssayHist->SetTitle( Form("%s [chi %.1f, slices %i, rebinX %i, EvInBin %.3f, rebin \"%s\", frag \"%s\"", pmtname.Data(), bestChisquare,
                                                              slices, rebinX, EvInBin, rebinopt.Data(), fragopt.Data() ) ); }

                            if( (bestChisquare<thebestChisquare)||(thebestChisquare < 0.) )
                            {
                                thebestChisquare = bestChisquare;
                                bestslices = slices;
                                bestrebinX = rebinX;
                                bestEvInBin = EvInBin;
                                bestrebinopt = rebinopt;
                                bestfragopt = fragopt;
                            }
                        }
                    }

                    if(currProjection) delete currProjection;

                }
            }
        }//rebin option


        if(thebestChisquare < 0.){ if(fVerbose >= 1)cout <<Form("WARING: Correction was NOT done for %s",pmtname.Data())<<endl; continue; }

        if(fVerbose >= 3)cout <<Form("\nBEST CHISQUARE: %.2f",bestChisquare)<<endl;

        ( (TH1*)fChiTrends->FindObject("stat_nslices") )->Fill(bestslices);
        ( (TH1*)fChiTrends->FindObject("stat_rebinX") )->Fill(bestrebinX);
        ( (TH1*)fChiTrends->FindObject("stat_EvInBin") )->Fill(bestEvInBin);
        ( (TH1*)fChiTrends->FindObject("stat_rebinopt") )->Fill(bestrebinopt,1);
        ( (TH1*)fChiTrends->FindObject("stat_fragopt") )->Fill(bestfragopt,1);

        currTrend->SetAxisRange(currTrend->GetMinimum(1)-100,currTrend->GetMaximum()+100,"Y");
        currTrend->SetTitle( Form("Chi-square trend for%s",pmtname.Data() ));

        fProjection->Add( bestProjection );
        fGraphs->Add( bestGraph );
        fAssay->Add( bestAssayHist );
        fChiTrends->Add( currTrend );

        if(fVerbose >= 2)cout <<Form("DONE FOR %s PMT",pmtname.Data())<<endl;
    }//PMT

    ResetArray( &rebinOptArr, "delete" );
    ResetArray( &fragOptArr, "delete" );

    return 1;
}
/*******************************************************************************************************************/
#endif // ALIT0TIMEAMPLCORR_CXX



