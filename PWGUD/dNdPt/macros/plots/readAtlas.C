Int_t readAtlas() 
{
// read ATLAS NSD data from file
//

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read ATLAS NSD data from file" <<endl;
cout << "Number of bins: " << binsAtlas <<endl;
cout << "Filename:       " << filenameAtlas <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

ifstream fileNsdAtlas;
fileNsdAtlas.open(filenameAtlas);

Int_t i = 0;
while(!fileNsdAtlas.eof()) {
    if(i == binsAtlas) break;
    // textfile content: pt:pt_low:pt_high:NSD_yield:stat_error:syst_error
    fileNsdAtlas >> ptAtlas[i] >> lowPtAtlas[i] >> highPtAtlas[i] >> nsdAtlas[i] >> statNsdAtlas[i] >> systNsdAtlas[i];
    
    centerPtAtlas[i]    = (highPtAtlas[i] + lowPtAtlas[i]) / 2.0;
    widthPtAtlas[i]     = highPtAtlas[i] - lowPtAtlas[i];
    errPtAtlas[i]       = widthPtAtlas[i] / 2.0;
    lowErrPtAtlas[i] = ptAtlas[i] - lowPtAtlas[i];
    highErrPtAtlas[i] = highPtAtlas[i] - ptAtlas[i];
            
    lowStatNsdAtlas[i]  = nsdAtlas[i] - statNsdAtlas[i];
    highStatNsdAtlas[i] = nsdAtlas[i] + statNsdAtlas[i];
    relStatNsdAtlas[i]  = statNsdAtlas[i] / nsdAtlas[i];    
    lowSystNsdAtlas[i]  = nsdAtlas[i] - systNsdAtlas[i];
    highSystNsdAtlas[i] = nsdAtlas[i] + systNsdAtlas[i];
    relSystNsdAtlas[i]  = systNsdAtlas[i] / nsdAtlas[i];    
    errNsdAtlas[i]      = systNsdAtlas[i] + statNsdAtlas[i];
    lowErrNsdAtlas[i]   = nsdAtlas[i] - errNsdAtlas[i];
    highErrNsdAtlas[i]  = nsdAtlas[i] + errNsdAtlas[i];
    relErrNsdAtlas[i]   = errNsdAtlas[i] / nsdAtlas[i];
    err2NsdAtlas[i]     = sqrt(systNsdAtlas[i]*systNsdAtlas[i] + statNsdAtlas[i]*statNsdAtlas[i]);
    lowErr2NsdAtlas[i]  = nsdAtlas[i] - err2NsdAtlas[i];
    highErr2NsdAtlas[i] = nsdAtlas[i] + errNsdAtlas[i];
    relErr2NsdAtlas[i]  = err2NsdAtlas[i] / nsdAtlas[i];
        
    cout << "ptAtlas[" << i << "]             = " << ptAtlas[i] <<endl;
    cout << "   centerPtAtlas[" << i << "]    = " << centerPtAtlas[i] <<endl;
    cout << "   widthPtAtlas[" << i << "]     = " << widthPtAtlas[i] <<endl;
    cout << "   errPtAtlas[" << i << "]       = " << errPtAtlas[i] <<endl;
    cout << "   lowErrPtAtlas[" << i << "]    = " << lowErrPtAtlas[i] <<endl;
    cout << "   highErrPtAtlas[" << i << "]   = " << highErrPtAtlas[i] <<endl;
    cout << "   lowPtAtlas[" << i << "]       = " << lowPtAtlas[i] <<endl;
    cout << "   highPtAtlas[" << i << "]      = " << highPtAtlas[i] <<endl;
    cout << "nsdAtlas[" << i << "]            = " << nsdAtlas[i] <<endl;
    cout << "   statNsdAtlas[" << i << "]     = " << statNsdAtlas[i] <<endl;
    cout << "   lowStatNsdAtlas[" << i << "]  = " << lowStatNsdAtlas[i] <<endl;
    cout << "   highStatNsdAtlas[" << i << "] = " << highStatNsdAtlas[i] <<endl;
    cout << "   relStatNsdAtlas[" << i << "]  = " << relStatNsdAtlas[i] <<endl;
    cout << "   systNsdAtlas[" << i << "]     = " << systNsdAtlas[i] <<endl;
    cout << "   lowSystNsdAtlas[" << i << "]  = " << lowSystNsdAtlas[i] <<endl;
    cout << "   highSystNsdAtlas[" << i << "] = " << highSystNsdAtlas[i] <<endl;
    cout << "   relSystNsdAtlas[" << i << "]  = " << relSystNsdAtlas[i] <<endl;
    cout << "errNsdAtlas[" << i << "]         = " << errNsdAtlas[i] <<endl;
    cout << "   lowErrNsdAtlas[" << i << "]   = " << lowErrNsdAtlas[i] <<endl;
    cout << "   highErrNsdAtlas[" << i << "]  = " << highErrNsdAtlas[i] <<endl;
    cout << "   relErrNsdAtlas[" << i << "]   = " << relErrNsdAtlas[i] <<endl;
    cout << "err2NsdAtlas[" << i << "]        = " << err2NsdAtlas[i] <<endl;
    cout << "   lowErr2NsdAtlas[" << i << "]  = " << lowErr2NsdAtlas[i] <<endl;
    cout << "   highErr2NsdAtlas[" << i << "] = " << highErr2NsdAtlas[i] <<endl;
    cout << "   relErr2NsdAtlas[" << i << "]  = " << relErr2NsdAtlas[i] <<endl;
    cout << endl;
     
   i++;
} // while(!fileNsdAtlas.eof())
fileNsdAtlas.close();
//if (fileNsdAtlas) { delete fileNsdAtlas; }
//fileNsdAtlas=0;

cout << "================================" << endl;
cout << endl;
cout << "Finished reading ATLAS NSD data" <<endl;
cout << "Number of bins read: " << i <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return i;
}