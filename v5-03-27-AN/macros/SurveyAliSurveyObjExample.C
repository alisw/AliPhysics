void SurveyAliSurveyObjExample(const char* user, const char* det, Int_t repNum, Int_t ver)
{
    // Example of use of AliSurveyObj functionality
    // user is the alien user allowing to connect to alien and open there the survey parsed file
    // repNum is the report number and ver it's version
    // As an example you could run this macro with the following line:
    // aliroot -b -q $ALICE_ROOT/macros/SurveyAliSurveyObjExample.C("myalienusername","TPC",818098,1)
    //

    AliSurveyObj *so = new AliSurveyObj();

    Int_t size = so->GetEntries();
    printf("-> %d\n", size);

    //  The survey object can be filled from local file or from the survey depot in alien.
    //  The following commented lines show an example, then we will use the second option.
    //  For this reason the alien user is required as argument.
    //  so->FillFromLocalFile("~/survey/real_data/Survey_781282_HMPID.txt");
    //  size = so->GetEntries();
    //  printf("--> %d\n", size);

    //  so->Fill("TRD", 2007, 816582, 1); 
    //  size = so->GetEntries();
    //  printf("----> %d\n", size);

    so->Fill(det, repNum, ver, user);
    size = so->GetEntries();
    printf("---> %d\n", size);

    Printf("Title: \"%s\"", so->GetReportTitle().Data());
    Printf("Date: \"%s\"", so->GetReportDate().Data());
    Printf("Detector: \"%s\"", so->GetDetector().Data());
    Printf("URL: \"%s\"", so->GetURL().Data());
    Printf("Number: \"%d\"", so->GetReportNumber());
    Printf("Version: \"%d\"", so->GetReportVersion());
    Printf("Observations: \"%s\"", so->GetObservations().Data());
    Printf("Coordinate System: \"%s\"", so->GetCoordSys().Data());
    Printf("Measurement Units: \"%s\"", so->GetUnits().Data());
    Printf("Nr Columns: \"%d\"", so->GetNrColumns());

    TObjArray *colNames = so->GetColumnNames();
    for (Int_t i = 0; i < colNames->GetEntries(); ++i)
	Printf("  Column %d --> \"%s\"", i, ((TObjString *) colNames->At(i))->GetString().Data());

    Printf("Points:");
    TObjArray *points = so->GetData();
    for (Int_t i = 0; i < points->GetEntries(); ++i)
	Printf("  Point %d --> \"%s\"", i, ((AliSurveyPoint *) points->At(i))->GetPointName().Data());

    // See "STEER/AliSurveyPoint.h" for more getters

    return;
}
