void SurveyAliSurveyObjExample()
{
  AliSurveyObj *so = new AliSurveyObj();
  
  Int_t size = so->GetEntries();
  printf("-> %d\n", size);
  
  //  so->FillFromLocalFile("~/survey/real_data/Survey_781282_HMPID.txt");
  //  size = so->GetEntries();
  //  printf("--> %d\n", size);

  so->Fill("HMPID", 2006, 781282, 1);
  size = so->GetEntries();
  printf("---> %d\n", size);

  //  so->Fill("TRD", 2007, 816582, 1); 
  //  size = so->GetEntries();
  //  printf("----> %d\n", size);

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
