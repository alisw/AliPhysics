# QnCorrections framework description {#mainpage}
  \tableofcontents

\section introduction Introduction
Measuring collective phenomena in heavy-ion collisions requires an estimate of the collision symmetry plane for different flow harmonics. The orientation of each plane is estimated from the azimutal distribution of produced particles. Detectors with nonuniform azimuthal acceptance distort such distribution. The goal of the correction framework is to provide an environment that is able to implement a set of corrections that compensates nonuniform azimuthal detector responses but also flexible enough to try out new or additional correction approaches with a minimum incremental design effort.

The framework models any experimental setup as a set of detectors (AliQnCorrectionsDetector) handled by a framework manager (AliQnCorrectionsManager). Each detector needs to be configured at least on one shape (AliQnCorrectionsDetectorConfigurationBase) but a detector could have multiple configurations (AliQnCorrectionsDetectorConfigurationSet). 

A detector configuration in its simplest way might be defined as a set of cuts applied to select which of the basic potential events happening to the detector are of interest for the considered detector configuration. For instance, if you are considering a tracking detector you could decide to define your detector configuration (AliQnCorrectionsDetectorConfigurationTracks) for your flow analysis by selecting the tracks with \f$p_T\f$ between 0.2 and 20.0 GeV/c.

Detector configuration allows you to handle channelized detectors or sub-detectors (AliQnCorrectionsDetectorConfigurationChannels) by allocating a different set of channels to each configuration and allowing you to group your channels within groups to so considering them in your flow analysis. Not channelized detectors, as tracking detectors, are able to be modelled with channelized configurations if so you decide. For instance for your tracking detector you could define a channelized configuration assigning 20 channels to your 0.2 – 20.0 GeV/c \f$p_T\f$ range.

![Framework functional diagram](Framework.png "Framework functional diagram")

Once you have defined the set of detector configurations that models your setup you will proceed to define the set of corrections (AliQnCorrectionsCorrectionStepBase) you will apply on each of these detector configurations. The framework provides support for two correction types: corrections on input (to the framework) data (AliQnCorrectionsCorrectionOnInputData), such as gain equalization, and Q vector corrections (AliQnCorrectionsCorrectionOnQvector). Corrections on input data are only able to be applied to channelized detector configurations while Q vector corrections can be applied to any detector configuration.

Each channelized detector configuration has a socket where the desired corrections on input data (AliQnCorrectionsCorrectionsSetOnInputData) plug in. The available corrections on input data are defined, at design time, with information regarding its applicability order. When the framework user decides to apply a concrete set of the available corrections on input data to a concrete detector configuration such information is used at run time to 'know' in which order the corrections will be applied.

Any detector configuration has a socket where the desired Q vector corrections (AliQnCorrectionsCorrectionsSetOnQvector) plug in. The principle is the same. When a new Q vector correction is designed for its first time, information on its applicability order is incorporated during its SW coding phase. This information is the only piece that the Q vector correction should provide to the framework for being able to run on it. The Q vector correction should follow a concrete template to be able of being invoked when so is required to perform its task.

Once fully configured, the correction framework might be used in three different modes. *Calibration* is the mode by which the framework produces all the necessary information for, in a latter phase, apply the desired corrections. Calibration mode usually requires several passes before the information for all correction steps is collected. *Correct* is the mode by which the framework provides corrected Q vector to user processes according to the selected configuration. *Mixed* is the mode by which the framework apply certain selected corrections to Q vectors and builds new calibration information. In mixed mode the user selects, at configuration time, the set of corrections to apply and the calibration information to build. This mode is intended for the development of new correction approaches.

\section frameworkdef Defining the framework

\subsection expsetup The experimental setup
Before starting to model your experimental setup you first need to consider which will be the scope of your flow analysis. You need to identify the set of detectors you will be focused on and under which conditions / configurations you will make use of the information they provide for each of the events you will analyse. Then you have to identify the different event classes in which you will classify your events and which variables you will make use of to build such event classes.

\subsection correctionmanager The framework manager
AliQnCorrectionsManager framework manager is the central piece of the correction framework. It is the anchor point between the framework and the external run time environment. The data flow comes in from the external environment and gets distributed to the set of defined detectors. At configuration time the different detectors are defined and assigned to the framework manager. 

You define the framework manager in a simple way
~~~{.cxx}
  /* create the framework manager */
  AliQnCorrectionsManager *QnManager = new AliQnCorrectionsManager();
~~~
and then you configure the basic functions that control the information produced by the framework
* produce an output TTree with Q vector data: for getting the subsequent corrected Q vectors written in a TTree structure.
* produce framework QA histograms: for producing and giving the framework defined QA histograms.
* build calibration information: for providing the information needed to produce the calibration / correction parameters.

~~~{.cxx}
  /* do not fill Qn vector TTree */
  QnManager->SetShouldFillQnVectorTree(kFALSE);
  /* do not fill framework QA histograms */
  QnManager->SetShouldFillQAHistograms(kFALSE);
  /* produce calibration information */
  QnManager->SetShouldFillOutputHistograms(kTRUE);
~~~

The framework supports running a set of its instances on a concurrent scenario so that you will get results from each of the running instances. To be able to allocate the results to different processes they correspond to getting them at the end properly merged, you declare the list of processes names the framework should globally handle
~~~{.cxx}
  /* store the list of concurrent processes names */
  QnManager->SetListOfProcessesNames(procNamesList);
~~~
and then, if you have already produced correction information in a previous step, you inform the framework about the file that includes it
~~~{.cxx}
  /* transfer the TFile with correction information */
  QnManager->SetCalibrationHistogramsList(calibfile);
~~~
Of course, the framework manager holds the set of detectors but they are defined next. The detectors are addressed by an external Id defined by the user but internally they are reached using an internal address which translation is performed by the framework manager. The framework manager also owns the data container used to interchange experimental setup variables values. 

\subsection detectors Defining detectors

AliQnCorrectionsDetector mirrors the experimental setup detectors within the correction framework. They are each externally identified by an unique detector Id that is passed to the framework at detector creation time together with the detector name to be used by the framework.
~~~{.cxx}
  /* create the new detector */
  AliQnCorrectionsDetector *VZERO = new AliQnCorrectionsDetector("VZERO", VAR::kVZERO);
~~~
The detector is in charge of holding the different configurations that the user has defined on it being its main task to properly address to them the incoming data flow. Some configurations could correspond to concrete subdetectors of the proper real detector but could as well be the own detector addressed by different set of cuts. In any case, it is task of the detector configuration and not of the detector to store / handle such characterization.

Once the different detector configurations have been incorporated, the detector is includen in the framework attaching it to the framework manager.
~~~{.cxx}
  /* add the detector to the framework manager */
  QnManager->AddDetector(VZERO);
~~~

\subsection detectorconfig Detector configurations


![Framework incoming dataflow](FrameworkDataFlow.png "Framework incoming dataflow")


\subsection eventclasses Event classes
If you are interesting in doing Flow analysis most probably you will want to classify your events under different event classes. For supporting that, the correction framework defines two classes: AliQnCorrectionsEventClassVariable and AliQnCorrectionsEventClassVariablesSet.

\subsubsection eventclassesvars Event classes (EC) variables
AliQnCorrectionsEventClassVariable allows you to identify one variable you want to use to classify your events, their z vertex coordinate, for instance, and associate it a concrete binning and text label to be used when such variale is presented on the different histograms.

As you know all variables of interest whitin the correction framework own a unique Id. This Id is also passed when you instantiate a AliQnCorrectionsEventClassVariable object.

The class incorporates a wide set of constructors that basically cover the different possibilities of declaring a concrete binning: explicit, uniform and with segments of different granularity.

Examples:
* constructor form an array of pairs, where the 1st element of each pair is the upper edge of a coarse bin, and the 2nd element is the number of fine bins inside the coarse bin. The 1st element of the 1st pair is the lowest value and the 2nd element of the first pair is the number of coarse bins plus one (i.e. the total number of pairs).
~~~{.cxx}
  Double_t VtxZaxes[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  AliQnCorrectionsEventClassVariable * fpVertexZEventClass = new AliQnCorrectionsEventClassVariable(kVertexZ, VarNames[kVertexZ], VtxZaxes);
~~~
will build an event class variable associated to the kVertexZ variable and with an associated binning {-10.0, -7.0, -5.25, -3.5, -1.75, 0.0, 1.75, 3.5, 5.25, 7.0, 10.0}.
* constructor from number of bins and minimum and maximum values.
~~~{.cxx}
  AliQnCorrectionsEventClassVariable * fpCentralityEventClass = new AliQnCorrectionsEventClassVariable(kCentrality, VarNames[kCentrality], 10, 0.0, 100.0);
~~~

\subsubsection eventclassesvarset Set of EC variables
AliQnCorrectionsEventClassVariablesSet encapsulates all variables you want to be used to classify your events, their z vertex coordinate and their multiplicities, for instance, as an array that contains the QnCorrectionEventClassVariable objects.

You are free to define an event classes variables set per detector configuration. This variables set will be incorporated to the different histograms the detector configuration will need to perform the intended corrections. The variables set will be the base for supporting the needed histogram multidimensionality.

Example:
* once you have built the previous variables you define the variables set and incorporate them
~~~{.cxx}
  const Int_t nEventClassesDimensions = 2;
  AliQnCorrectionsEventClassVariablesSet fEventClasses = AliQnCorrectionsEventClassVariablesSet(nEventClassesDimensions);
  fEventClasses[0] = fpVertexZEventClass;
  fEventClasses[1] = fpCentralityEventClass;
~~~
now you are able to access to your event classes variables and their access members
  * with an iterator 
~~~{.cxx}
  TObjArrayIter *myiter = new TObjArrayIter(&fEventClasses); // requieres a pointer that's what '&' does 
  AliQnCorrectionsEventClassVariable *nextVar = (AliQnCorrectionsEventClassVariable*) myiter->Next();

  while (nextVar) {
    printf("Variable id: %d\n  name: %s\n  bins: %f", nextVar->GetVariableId(), nextVar->GetVariableLabel(), nextVar->GetBinEdge(0));
    for (Int_t bin = 1; bin < nextVar->GetNBins() + 1; bin++) {
      printf(", %f", nextVar->GetBinEdge(bin));
    }
    printf("\n");
    nextVar = (AliQnCorrectionsEventClassVariable*) myiter->Next();
  }
~~~
  * or directly
~~~{.cxx}
  for (Int_t ixvar = 0; ixvar < fEventClasses.GetEntriesFast(); ixvar++) {
    printf("Variable id: %d\n  name: %s\n  bins: %f",
        fEventClasses.At(ixvar)->GetVariableId(),
        fEventClasses.At(ixvar)->GetVariableLabel(),
        fEventClasses.At(ixvar)->GetBinEdge(0));
    for (Int_t bin = 1; bin < fEventClasses.At(ixvar)->GetNBins() + 1; bin++) {
      printf(", %f", fEventClasses.At(ixvar)->GetBinEdge(bin));
    }
    printf("\n");
  }
~~~

\subsection histograms Histograms
The correction framework makes an extensive use of profile histograms in order to provide the proper scope to the different averages needed to elaborate the intended Q vector corrections. As was previouslly described one of the key components of the correction framework is its concept of multidimensional event classes: you decide to define your event classes based on a concrete set of variables usually greather than one. This goes straight on to multidimensional profiles. 

So far ROOT only supports up to tridimensional profiles, but we wanted to give more flexibility and, due to the fact we expect not that much bins per variable you will use to define your event classes, we decided to build multidimensional profiles functionality. With this aim we provide a base class, AliQnCorrectionsHistogramBase, that should not be instantiated, and which declares the whole interface for the correction framework histogram classes, a single multidimensional histogram profile class, AliQnCorrectionsProfile, a components oriented multidimensional histogram profile class, AliQnCorrectionsProfileComponents, which provide support for both components, X and Y, for a set of selected harmonics, and a correlation components oriented multidimenional histogram profile class, AliQnCorrectionsProfileCorrelationComponents, which provide support for the correlation components, XX, XY, YX and YY, for a set of selected harmonics. 

\subsubsection histogramsbase Histograms base
AliQnCorrectionsHistogramBase is the base for the whole set of histogram classes and declares their interfaces for the basic methods. It is not oriented to be instantiated and as an additional development support will include run time error messages that will orient you when you are not doing a proper use of its descendant classes.

\subsubsection profile Basic profiles
AliQnCorrectionsProfile implements a multidimensional histogram profile. You construct one of them passing as parameter an event classes variables set. Then you ask for the actual support histograms creation and then you use it as you will expect from a conventional histogram.

Example:
* based on the event classes variables set you had created above
~~~{.cxx}
    AliQnCorrectionsProfile *myProfile = new AliQnCorrectionsProfile("AliQnCorrectionsProfile", "myProfile", fEventClasses);
    myProfile->CreateProfileHistograms(myList); // the list should collect all histograms you want them going to a file latter on
~~~
* now your code will fill a data bank structure that the framework will receive with each singular event. This data bank of course includes the actual values of the variables you identified as the event classes variables. So, when is time to fill your profile you do it, as you expected,
~~~{.cxx}
    myProfile->Fill(varContainer, weight);
~~~
* accessing the content of your profile is straight forward, as you expected,
~~~{.cxx}
    Double_t myProfileBinContent = myProfile->GetBinContent(myProfile->GetBin(varContainer));
    Double_t myProfileBinError = myProfile->GetBinError(myProfile->GetBin(varContainer));
~~~

\subsubsection compprofile Components profile
AliQnCorrectionsProfileComponents implements a multidimensional histogram profile for each of the components X and Y for each of the selected harmonics. The harmonic number is handled in the expected way, i.e., if you are interested in the second harmonic, you will always address it as the second (2) harmonic. But you have to provide some help at histogram creation time.

At construction time you provide the same information as for any other histogram profile within the framework: an event classs variable set
~~~{.cxx}
    AliQnCorrectionsProfileComponents *myProfile = new AliQnCorrectionsProfileComponents("AliQnCorrectionsProfileComponents", "myComponentsProfile", fEventClasses);
~~~
But at creation time, if you want full support for your criteria of numbering harmonics you have to provide a map to indicate so to the framework. For instance, suppose you need support for three harmonics, the even ones, 2, 4 and 6. Then you write
~~~{.cxx}
  Int_t nNoOfHarmonics = 3;
  Int_t harmonicsMap[] = {2,4,6};

  myProfile->CreateComponentsProfileHistograms(myList, nNoOfHarmonics, harmonicsMap);
~~~
Now you have access to each component in an explicit way both for filling the histograms and for getting its bin content. So, as you would expect
~~~{.cxx}
    Int_t myHarmonic = 2;
    myProfile->FillX(myHarmonic, varContainer, weight2X);
    myProfile->FillY(myHarmonic, varContainer, weight2Y);
    myHarmonic = 4;
    myProfile->FillX(myHarmonic, varContainer, weight4X);
    myProfile->FillY(myHarmonic, varContainer, weight4Y);
    myHarmonic = 6;
    myProfile->FillX(myHarmonic, varContainer, weight6X);
    myProfile->FillY(myHarmonic, varContainer, weight6Y);
~~~
The framework expects you to fill both components for the whole set of harmonics you asked support for each time you perform a fill so, if you fail doing that the framework will raise an error indicating you you are not doing what it is expected you to do.
Accessing the content of your profile is straight forward, as you expected,
~~~{.cxx}
    myHarmonic = 2;
    Double_t myProfile2XBinContent = myProfile->GetXBinContent(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2YBinContent = myProfile->GetYBinContent(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2XBinError = myProfile->GetXBinError(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2YBinError = myProfile->GetYBinError(myHarmonic, myProfile->GetBin(varContainer));
~~~
 
\subsubsection compcorrprofile Correlation components
AliQnCorrectionsProfileCorrelationComponents implements a multidimensional histogram profile for each of the correlation components XX, XY, YX and YY for each of the selected harmonics. The overall behavior matches the AliQnCorrectionsProfileComponents one so, we just include the adapted code snipets
~~~{.cxx}
    AliQnCorrectionsProfileCorrelationComponents *myProfile = 
        new AliQnCorrectionsProfileCorrelationComponents("AliQnCorrectionsProfileCorrelationComponents", "myComponentsCorrelationProfile", fEventClasses);
~~~
~~~{.cxx}
  Int_t nNoOfHarmonics = 1;
  Int_t harmonicsMap[] = {2};

  myProfile->CreateComponentsCorrelationProfileHistograms(myList, nNoOfHarmonics, harmonicsMap);
~~~
~~~{.cxx}
    Int_t myHarmonic = 2;
    myProfile->FillXX(myHarmonic, varContainer, weight2XX);
    myProfile->FillXY(myHarmonic, varContainer, weight2XY);
    myProfile->FillYX(myHarmonic, varContainer, weight2YX);
    myProfile->FillYY(myHarmonic, varContainer, weight2YY);
~~~
The framework expects you to fill both components for the whole set of harmonics you asked support for each time you perform a fill so, if you fail doing that the framework will raise an error indicating you you are not doing what it is expected you to do.
~~~{.cxx}
    myHarmonic = 2;
    Double_t myProfile2XXBinContent = myProfile->GetXXBinContent(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2XYBinContent = myProfile->GetXYBinContent(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2YXBinContent = myProfile->GetYXBinContent(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2YYBinContent = myProfile->GetYYBinContent(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2XXBinError = myProfile->GetXXBinError(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2XYBinError = myProfile->GetXYBinError(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2YXBinError = myProfile->GetYXBinError(myHarmonic, myProfile->GetBin(varContainer));
    Double_t myProfile2YYBinError = myProfile->GetYYBinError(myHarmonic, myProfile->GetBin(varContainer));
~~~


