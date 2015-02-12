#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace AAF;

#pragma link C++ function AAF::FILTER_ESDMUON;

#pragma link C++ function AAF::FILTER_RAWMUON;
#pragma link C++ namespace AAF::RAWMUON;
#pragma link C++ function AAF::RAWMUON::CheckFile;
#pragma link C++ function AAF::RAWMUON::DisableBranches;

#pragma link C++ function AAF::FILTER_AODMUONWITHTRACKLETS;

#pragma link C++ function AAF::StageAndFilter;

#endif

