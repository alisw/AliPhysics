Int_t readCms() 
{
// read CMS NSD data from file
//

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read CMS NSD data from file" <<endl;
cout << "Number of bins: " << binsCms <<endl;
cout << "Filename:       " << filenameCms <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

ifstream fileNsdCms;
fileNsdCms.open(filenameCms);

Int_t i = 0;
while(!fileNsdCms.eof()) {
    if(i == binsCms) break;
    // textfile content: pt:NSD_yield:(stat + syst error added linerarly)
    fileNsdCms >> centerPtCms[i] >> nsdCms[i] >> errNsdCms[i];
    
    ptCms[i]          = centerPtCms[i];
    
    // the width of the pt bins, currently hardwired....
    widthPtCms[i] = 0.1;
    if (centerPtCms[i] > 1) widthPtCms[i] = 0.2;
    
    errPtCms[i]       = widthPtCms[i] / 2.0;              
    lowPtCms[i]       = centerPtCms[i] - errPtCms[i];
    highPtCms[i]      = centerPtCms[i] + errPtCms[i];
    
    /*
    lowStatNsdCms[i]  = nsdCms[i] - statNsdCms[i];
    highStatNsdCms[i] = nsdCms[i] + statNsdCms[i];
    relStatNsdCms[i]  = statNsdCms[i] / nsdCms[i];    
    lowSystNsdCms[i]  = nsdCms[i] - systNsdCms[i];
    highSystNsdCms[i] = nsdCms[i] + systNsdCms[i];
    relSystNsdCms[i]  = systNsdCms[i] / nsdCms[i];    
    */    
    lowErrNsdCms[i]   = nsdCms[i] - errNsdCms[i];
    highErrNsdCms[i]  = nsdCms[i] + errNsdCms[i];
    relErrNsdCms[i]   = errNsdCms[i] / nsdCms[i];
    err2NsdCms[i]     = errNsdCms[i];
    lowErr2NsdCms[i]  = lowErrNsdCms[i];
    highErr2NsdCms[i] = highErrNsdCms[i];
    relErr2NsdCms[i]  = relErrNsdCms[i];
        
    cout << "ptCms[" << i << "]             = " << ptCms[i] <<endl;
    cout << "   centerPtCms[" << i << "]    = " << centerPtCms[i] <<endl;
    cout << "   widthPtCms[" << i << "]     = " << widthPtCms[i] <<endl;
    cout << "   errPtCms[" << i << "]       = " << errPtCms[i] <<endl;
    cout << "   lowPtCms[" << i << "]       = " << lowPtCms[i] <<endl;
    cout << "   highPtCms[" << i << "]      = " << highPtCms[i] <<endl;
    cout << "nsdCms[" << i << "]            = " << nsdCms[i] <<endl;
    /*
    cout << "   statNsdCms[" << i << "]     = " << statNsdCms[i] <<endl;
    cout << "   lowStatNsdCms[" << i << "]  = " << lowStatNsdCms[i] <<endl;
    cout << "   highStatNsdCms[" << i << "] = " << highStatNsdCms[i] <<endl;
    cout << "   relStatNsdCms[" << i << "]  = " << relStatNsdCms[i] <<endl;
    cout << "   systNsdCms[" << i << "]     = " << systNsdCms[i] <<endl;
    cout << "   lowSystNsdCms[" << i << "]  = " << lowSystNsdCms[i] <<endl;
    cout << "   highSystNsdCms[" << i << "] = " << highSystNsdCms[i] <<endl;
    cout << "   relSystNsdCms[" << i << "]  = " << relSystNsdCms[i] <<endl;
    */
    cout << "errNsdCms[" << i << "]         = " << errNsdCms[i] <<endl;
    cout << "   lowErrNsdCms[" << i << "]   = " << lowErrNsdCms[i] <<endl;
    cout << "   highErrNsdCms[" << i << "]  = " << highErrNsdCms[i] <<endl;
    cout << "   relErrNsdCms[" << i << "]   = " << relErrNsdCms[i] <<endl;
    cout << "err2NsdCms[" << i << "]        = " << err2NsdCms[i] <<endl;
    cout << "   lowErr2NsdCms[" << i << "]  = " << lowErr2NsdCms[i] <<endl;
    cout << "   highErr2NsdCms[" << i << "] = " << highErr2NsdCms[i] <<endl;
    cout << "   relErr2NsdCms[" << i << "]  = " << relErr2NsdCms[i] <<endl;
    cout << endl;
     
   i++;
} // while(!fileNsdCms.eof())
fileNsdCms.close();
//if (fileNsdCms) { delete fileNsdCms; }
//fileNsdCms=0;

cout << "================================" << endl;
cout << endl;
cout << "Finished reading CMS NSD data" <<endl;
cout << "Number of bins read: " << i <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return i;
}