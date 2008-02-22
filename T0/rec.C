 void rec()
{
  AliLog::SetModuleDebugLevel("T0", 1);
  AliReconstruction rec;
  rec.SetRunReconstruction("T0");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunTracking("");

  rec.SetNumberOfEventsPerFile(-1);
  rec.SetRunQA(kFALSE);

  //  rec.SetEventRange(0, 1000); 
 
  // rec.SetOption("T0","cosmic");
  // rec.SetInput("/home/alla/alice/testFeb08/08000019174005.10.root"); 

    rec.SetOption("T0","pdc");

  rec.Run();
}
