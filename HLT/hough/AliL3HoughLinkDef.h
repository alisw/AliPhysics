// @(#) $Id$

#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliL3Hough; 
#pragma link C++ class AliL3HoughBaseTransformer; 
#pragma link C++ class AliL3HoughTransformer;
#pragma link C++ class AliL3HoughTransformerLUT;
#pragma link C++ class AliL3HoughTransformerVhdl;
#pragma link C++ class AliL3HoughTransformerNew;
#pragma link C++ class AliL3HoughTransformerRow;
#ifndef macosx
#pragma link C++ class AliL3HoughTrack;
#pragma link C++ class AliL3HoughKalmanTrack;
#endif
#pragma link C++ class AliL3HoughMaxFinder;
#pragma link C++ class AliL3HoughEval;
#pragma link C++ class AliL3Histogram;
#pragma link C++ class AliL3Histogram1D;
#pragma link C++ class AliL3HoughMerger;
#pragma link C++ class AliL3HoughIntMerger;
#pragma link C++ class AliL3HoughGlobalMerger;
#pragma link C++ class AliL3HoughDisplay;
#pragma link C++ class AliL3HoughClusterTransformer;
#pragma link C++ class AliL3HistogramAdaptive;
#ifndef macosx
#pragma link C++ class AliL3HoughTest;
#endif

#ifdef use_aliroot
#pragma link C++ class AliL3HoughTransformerGlobal;
#endif

#endif
