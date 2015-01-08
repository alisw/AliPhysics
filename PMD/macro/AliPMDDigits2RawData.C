// ----------------------------------------------------//
//                                                     //
//       This macro does digits to Raw Data            //
//                                                     //
// ----------------------------------------------------//


Int_t AliPMDDigits2RawData()
{
  AliSimulation sim;
  sim.WriteRawData("PMD");

  return 0;
}

