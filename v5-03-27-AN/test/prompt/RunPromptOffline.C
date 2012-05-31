Bool_t RunPromptOffline(Int_t run, Int_t gdcNumber, TString trg="" )
{
    //
    //  origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
    //

    //check arguments
    if (gdcNumber>999) return kFALSE;

    gROOT->LoadMacro("/local/home/daq/alisoft/macros/grp.C");
    Int_t grpRetCode = grp(run);
    if (grpRetCode<=0) return kFALSE;

    TString gdcNumberStr = "";
    gdcNumberStr += gdcNumber;
    if (gdcNumberStr.Length()==2) gdcNumberStr.Prepend("0");
    TString datasource = "mem://@aldaqpc";
    datasource.Append(gdcNumberStr);
    datasource.Append(":");  //no trg, do nothing

    //handle the low-level trigger selection
    if (trg!="")
    {
      datasource.Append("?Trigger=");
      datasource.Append(trg);
    }

    cout<<endl<<"RunPromptOffline datasource: "<<datasource<<endl<<endl;

    gROOT->LoadMacro("./rec.C");
    rec(datasource);
}
