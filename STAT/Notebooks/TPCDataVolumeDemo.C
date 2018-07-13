///\file TPCDataVolumeDemo.C
///\ingroup Notebooks
///\notebook -nodraw
///\macro_output
///\macro_code
///\author Marian Ivanov
///\brief TPCDataVolume demo to be run to check functionality  of JUPYTER notebook TPCdataVolumeDemo.ipynb
/// =================================================================================================
/// Demo showing:
/// ------------------------------------
///    * AliExternalInfo functionality to load input data
///    * AliTreePlayer to query tree MetaData information (sude for web servewer, drawing axis,titles, ...)
///    * TStattToolkit::MakeGraph<xxx> for tree queries + drawing
/// \include TPCVolumeDemo.C.log
/// \includelineno TPCVolumeDemo.C.log
/// See more detailed description in TPCdataVolumeDemo.ipynb
/// \code
///    .L $AliRoot_SRC/STAT/Notebooks/TPCDataVolumeDemo.C
///    TPCDataVolumeDemo()
/// \endcode

void TPCDataVolumeDemo(){
  /// 1. Setup input parameters of macro
  /// =====================================
  ///  Declaring variables and setting parametersTo be exercized - study data volume for different period
  ///  PbPb data period="LHC15o";
  ///   * all 2017 data period="LHC17*";
  ///   * to use star convention (load all period according mask ) data has to be cached before)
  /// 0.) Define input variable
  AliExternalInfo info;
  TTree * treeTPC=0;
  TString period="LHC15o";   // enter here period of interest
  TString pass="pass1";      // enter here pass of interest
  TCanvas *canvasDraw = new TCanvas("canvasDraw","canvasDraw",800,500);
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptFit(1);
  Int_t entries;
  ///
  /// 2. Initialize input data - AliExternalInfo DB and define derived varaibles
  /// ==============================================
  ///\code
  ///AliExternalInfo info;
  ///treeTPC = info.GetTree("QA.TPC","LHC15o","cpass1_pass1","QA.rawTPC;QA.TPC;QA.EVS;Logbook;Logbook.detector:TPC:detector==\"TPC\"");
  ///\endcode
  //period="LHC18*"; pass="cpass1_pass1"; ///- enable this line to analyze all 2018
  info.fVerbose=7;
  treeTPC = info.GetTree("QA.TPC",period,pass,"QA.rawTPC;QA.TPC;QA.EVS;Logbook;Logbook.detector:TPC:detector==\"TPC\"");
  //treeTPC = info.GetChain("QA.TPC",period,pass,"QA.rawTPC;QA.TPC;QA.EVS;Logbook;Logbook.detector");
  treeTPC->SetMarkerStyle(21);
  treeTPC->SetAlias("normDataVolume","Logbook.detector_TPC.bytesInjectedPhysics/Logbook.detector_TPC.eventCountPhysics/1000000");
  ///
  /// 3. Query list of the variables of given class using tree functionality and AliTreePlayer Metadata functionality
  /// ================================================
  ///   *  Simple tree functionality to dump tree content (and data)
  /// 3.a) using Show
  /// \code   treeTPC->GetFriend("Logbook.detector_TPC")->Show(0); \endcode
  treeTPC->GetFriend("Logbook.detector_TPC")->Show(0);
  ///   * Simple query TPC dEdx variables
  /// \code AliTreePlayer::selectMetadata(treeTPC->GetFriend("QA.TPC"), "[class==\"TPC&&dEdx\"]",0)->Print(); \endcode
  AliTreePlayer::selectMetadata(treeTPC->GetFriend("QA.TPC"), "[class==\"TPC&&dEdx\"]",0)->Print();
  ///   * Simple query TPC efficiency variables
  ///   \code AliTreePlayer::selectMetadata(treeTPC->GetFriend("QA.TPC"), "[class==\"TPC&&Eff\"]",0)->Print(); \endcode
  AliTreePlayer::selectMetadata(treeTPC->GetFriend("QA.TPC"), "[class==\"TPC&&Eff\"]",0)->Print();
  ///   * Refined query TPC DCAz variables which are not of class Err
  ///   \code  AliTreePlayer::selectMetadata(treeTPC->GetFriend("QA.TPC"), "[class==\"TPC&&(DCAz&&!Err)\"]",0)->Print(); \endcode
  AliTreePlayer::selectMetadata(treeTPC->GetFriend("QA.TPC"), "[class==\"TPC&&(DCAz&&!Err)\"]",0)->Print();
  ///   * Select logboook time information source
  /// \code AliTreePlayer::selectMetadata(treeTPC->GetFriend("Logbook"), "[class==\"Logbook&&Time\"]",0)->Print(); \endcode
  AliTreePlayer::selectMetadata(treeTPC->GetFriend("Logbook"), "[class==\"Logbook&&Time\"]",0)->Print();
  ///   * Select logboook event stat  information source
  /// \code AliTreePlayer::selectMetadata(treeTPC->GetFriend("Logbook"), "[class==\"Logbook&&Stat\"]",0)->Print();  \endcode
  AliTreePlayer::selectMetadata(treeTPC->GetFriend("Logbook"), "[class==\"Logbook&&Stat\"]",0)->Print();
  ///
  /// 4. Make and export fits Fit dependece of data voume as function of Interaction rate
  /// ====================================================================================
  canvasDraw = new TCanvas("canvasDraw","canvasDraw",600,400);
  /// 4.a Non robust fit using the TProfile
  /// ----------------------------------------
  /// \code
  /// entries = treeTPC->Draw("normDataVolume:QA.EVS.interactionRate","","prof");
  /// treeTPC->GetHistogram()->Fit("pol1","");
  ///\endcode
  /// ![Graph of data volume as function of rate fitted by robust pol1 fit ](TPCDataVolumeDemoFig4a.png)
  treeTPC->SetAlias("normDataVolume","Logbook.detector_TPC.bytesInjectedPhysics/Logbook.detector_TPC.eventCountPhysics/1000000");
  entries = treeTPC->Draw("normDataVolume:QA.EVS.interactionRate","","prof");
  treeTPC->GetHistogram()->Fit("pol1","");
  canvasDraw->SaveAs("$AliRoot_SRC/STAT/imgdoc/TPCDataVolumeDemoFig4a.png");
  /// 4.b Creating graphs using TStatToolkit::MakeGraphErrors and fitting using robust LTM fits
  ///-------------------------------------------------------------------------------------------
  /// \code
  /// grNorm=TStatToolkit::MakeGraphErrors(treeTPC,"normDataVolume:QA.EVS.interactionRate","1",25,1);
  /// grNorm->Draw("ap");
  /// grNorm->Fit("pol1","rob=0.8","");
  ///\endcode
  /// ![Graph of data volume as function of rate fitted by robust pol1 fit ](TPCDataVolumeDemoFig4b.png)
  gStyle->SetOptFit(1);
  auto grNorm=TStatToolkit::MakeGraphErrors(treeTPC,"normDataVolume:QA.EVS.interactionRate","1",25,1);
  grNorm->Draw("ap");
  grNorm->Fit("pol1","rob=0.8","");
  treeTPC->SetAlias("normDataVolumeFit",TString::Format("(%f+%f*QA.EVS.interactionRate)",grNorm->GetFunction("pol1")->GetParameter(0), grNorm->GetFunction("pol1")->GetParameter(1)).Data());
  canvasDraw->SaveAs("$AliRoot_SRC/STAT/imgdoc/TPCDataVolumeDemoFig4b.png");
  ///
  /// 5. Correlating normalizedData volume with other observables
  /// ============================================================
  /// 5.a Draw using standard root queries
  /// \code treeTPC->Draw("normDataVolume/normDataVolumeFit:run:QA.TPC.bz","","colz"); \endcode
  /// ![Ratio of data volume to fit per run standrd tree draw query](TPCDataVolumeDemoFig5a.png)
  treeTPC->Draw("normDataVolume/normDataVolumeFit:run:QA.TPC.bz","","colz");
  canvasDraw->SaveAs("$AliRoot_SRC/STAT/imgdoc/TPCDataVolumeDemoFig5a.png");
  /// 5.b Sparse drawing  using TStatToolkit MakeGraphSparse
  /// \code TStatToolkit::MakeGraphSparse(treeTPC,"normDataVolume/normDataVolumeFit:run","1",25,1,1,0)->Draw("apl");\endcode
  /// ![Ratio of data volume to fit per run - sparse represenation](TPCDataVolumeDemoFig5b.png)
  TStatToolkit::MakeGraphSparse(treeTPC,"normDataVolume/normDataVolumeFit:run","1",25,1,1,0)->Draw("apl");
  canvasDraw->SaveAs("$AliRoot_SRC/STAT/imgdoc/TPCDataVolumeDemoFig5b.png");
  /// 5.c parse drawing  using TStatToolkit MakeGraphSparse
  /// \code treeTPC->Draw("normDataVolume:QA.EVS.interactionRate:QA.TPC.bz","","colz"); \endcode
  /// ![Normalized data volume as function of rate and Bz(color code)](TPCDataVolumeDemoFig5c.png)
  treeTPC->Draw("normDataVolume:QA.EVS.interactionRate:QA.TPC.bz","","colz");
  canvasDraw->SaveAs("$AliRoot_SRC/STAT/imgdoc/TPCDataVolumeDemoFig5c.png");
  ///
}