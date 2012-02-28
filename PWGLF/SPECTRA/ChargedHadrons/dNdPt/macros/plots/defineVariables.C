//------------------------------------------------------------------------------
// defineVariables.C
//
// define all Variables that are needed to store data and uncertainties
//------------------------------------------------------------------------------


//
// define more arrays to store ALICE NSD data points
//     
                                             // for NSD
Double_t ptNsd2PiPtAlice[binsNsdAlice] = {0};            // pt   value (weighted)
Double_t centerPtNsd2PiPtAlice[binsNsdAlice] = {0};      //      center of bin
Double_t widthPtNsd2PiPtAlice[binsNsdAlice] = {0};       //      width
Double_t errPtNsd2PiPtAlice[binsNsdAlice] = {0};         //      error (half width)
Double_t lowPtNsd2PiPtAlice[binsNsdAlice] = {0};         //      lower edge
Double_t highPtNsd2PiPtAlice[binsNsdAlice] = {0};        //      upper edge

Double_t nsd2PiPtAlice[binsNsdAlice] = {0};                 // NSD
Double_t statNsd2PiPtAlice[binsNsdAlice] = {0};             //      stat error
Double_t lowStatNsd2PiPtAlice[binsNsdAlice] = {0};          //                  lower
Double_t highStatNsd2PiPtAlice[binsNsdAlice] = {0};         //                  upper
Double_t relStatNsd2PiPtAlice[binsNsdAlice] = {0};          //                  relative
Double_t systNsd2PiPtAlice[binsNsdAlice] = {0};             //      syst error
Double_t lowSystNsd2PiPtAlice[binsNsdAlice] = {0};          //                  lower
Double_t highSystNsd2PiPtAlice[binsNsdAlice] = {0};         //                  upper
Double_t relSystNsd2PiPtAlice[binsNsdAlice] = {0};          //                  relative
Double_t errNsd2PiPtAlice[binsNsdAlice] = {0};              //      stat+syst added linearly
Double_t lowErrNsd2PiPtAlice[binsNsdAlice] = {0};           //                  lower
Double_t highErrNsd2PiPtAlice[binsNsdAlice] = {0};          //                  upper
Double_t relErrNsd2PiPtAlice[binsNsdAlice] = {0};           //                  relative
Double_t err2Nsd2PiPtAlice[binsNsdAlice] = {0};             //      stat+syst added quadratically
Double_t lowErr2Nsd2PiPtAlice[binsNsdAlice] = {0};          //                  lower
Double_t highErr2Nsd2PiPtAlice[binsNsdAlice] = {0};         //                  upper
Double_t relErr2Nsd2PiPtAlice[binsNsdAlice] = {0};          //                  relative

//
// define more arrays to store ALICE INEL data points
//  
                                                // for INEL
Double_t ptInel2PiPtAlice[binsInelAlice] = {0};          // pt   value (weighted)
Double_t centerPtInel2PiPtAlice[binsInelAlice] = {0};    //      center of bin
Double_t widthPtInel2PiPtAlice[binsInelAlice] = {0};     //      width
Double_t errPtInel2PiPtAlice[binsInelAlice] = {0};       //      error (half width)
Double_t lowPtInel2PiPtAlice[binsInelAlice] = {0};       //      lower edge
Double_t highPtInel2PiPtAlice[binsInelAlice] = {0};      //      upper edge

Double_t inel2PiPtAlice[binsInelAlice] = {0};                // INEL
Double_t statInel2PiPtAlice[binsInelAlice] = {0};            //      stat error
Double_t lowStatInel2PiPtAlice[binsInelAlice] = {0};         //                  lower
Double_t highStatInel2PiPtAlice[binsInelAlice] = {0};        //                  upper
Double_t relStatInel2PiPtAlice[binsInelAlice] = {0};         //                  relative
Double_t systInel2PiPtAlice[binsInelAlice] = {0};            //      syst error
Double_t lowSystInel2PiPtAlice[binsInelAlice] = {0};         //                  lower
Double_t highSystInel2PiPtAlice[binsInelAlice] = {0};        //                  upper
Double_t relSystInel2PiPtAlice[binsInelAlice] = {0};         //                  relative
Double_t errInel2PiPtAlice[binsInelAlice] = {0};             //      stat+syst added linearly
Double_t lowErrInel2PiPtAlice[binsInelAlice] = {0};          //                  lower
Double_t highErrInel2PiPtAlice[binsInelAlice] = {0};         //                  upper
Double_t relErrInel2PiPtAlice[binsInelAlice] = {0};          //                  relative
Double_t err2Inel2PiPtAlice[binsInelAlice] = {0};            //      stat+syst added quadratically
Double_t lowErr2Inel2PiPtAlice[binsInelAlice] = {0};         //                  lower
Double_t highErr2Inel2PiPtAlice[binsInelAlice] = {0};        //                  upper
Double_t relErr2Inel2PiPtAlice[binsInelAlice] = {0};         //                  relative

//
// define more arrays to store ALICE YIELD data points
//  
                                                // for YIELD
Double_t ptYield2PiPtAlice[binsYieldAlice] = {0};        // pt   value (weighted)
Double_t centerPtYield2PiPtAlice[binsYieldAlice] = {0};  //      center of bin
Double_t widthPtYield2PiPtAlice[binsYieldAlice] = {0};   //      width
Double_t errPtYield2PiPtAlice[binsYieldAlice] = {0};     //      error (half width)
Double_t lowPtYield2PiPtAlice[binsYieldAlice] = {0};     //      lower edge
Double_t highPtYield2PiPtAlice[binsYieldAlice] = {0};    //      upper edge

Double_t yield2PiPtAlice[binsYieldAlice] = {0};               // Invariant Yield
Double_t statYield2PiPtAlice[binsYieldAlice] = {0};           //      stat error
Double_t lowStatYield2PiPtAlice[binsYieldAlice] = {0};        //                  lower
Double_t highStatYield2PiPtAlice[binsYieldAlice] = {0};       //                  upper
Double_t relStatYield2PiPtAlice[binsYieldAlice] = {0};        //                  relative
Double_t systYield2PiPtAlice[binsYieldAlice] = {0};           //      syst error
Double_t lowSystYield2PiPtAlice[binsYieldAlice] = {0};        //                  lower
Double_t highSystYield2PiPtAlice[binsYieldAlice] = {0};       //                  upper
Double_t relSystYield2PiPtAlice[binsYieldAlice] = {0};        //                  relative
Double_t errYield2PiPtAlice[binsYieldAlice] = {0};            //      stat+syst added linearly
Double_t lowErrYield2PiPtAlice[binsYieldAlice] = {0};         //                  lower
Double_t highErrYield2PiPtAlice[binsYieldAlice] = {0};        //                  upper
Double_t relErrYield2PiPtAlice[binsYieldAlice] = {0};         //                  relative
Double_t err2Yield2PiPtAlice[binsYieldAlice] = {0};           //      stat+syst added quadratically
Double_t lowErr2Yield2PiPtAlice[binsYieldAlice] = {0};        //                  lower
Double_t highErr2Yield2PiPtAlice[binsYieldAlice] = {0};       //                  upper
Double_t relErr2Yield2PiPtAlice[binsYieldAlice] = {0};        //                  relative






//
// define arrays to store ALICE NSD data points
//     
                                             // for NSD
Double_t ptNsdAlice[binsNsdAlice] = {0};            // pt   value (weighted)
Double_t centerPtNsdAlice[binsNsdAlice] = {0};      //      center of bin
Double_t widthPtNsdAlice[binsNsdAlice] = {0};       //      width
Double_t errPtNsdAlice[binsNsdAlice] = {0};         //      error (half width)
Double_t lowPtNsdAlice[binsNsdAlice] = {0};         //      lower edge
Double_t highPtNsdAlice[binsNsdAlice] = {0};        //      upper edge

Double_t nsdAlice[binsNsdAlice] = {0};                 // NSD
Double_t statNsdAlice[binsNsdAlice] = {0};             //      stat error
Double_t lowStatNsdAlice[binsNsdAlice] = {0};          //                  lower
Double_t highStatNsdAlice[binsNsdAlice] = {0};         //                  upper
Double_t relStatNsdAlice[binsNsdAlice] = {0};          //                  relative
Double_t systNsdAlice[binsNsdAlice] = {0};             //      syst error
Double_t lowSystNsdAlice[binsNsdAlice] = {0};          //                  lower
Double_t highSystNsdAlice[binsNsdAlice] = {0};         //                  upper
Double_t relSystNsdAlice[binsNsdAlice] = {0};          //                  relative
Double_t errNsdAlice[binsNsdAlice] = {0};              //      stat+syst added linearly
Double_t lowErrNsdAlice[binsNsdAlice] = {0};           //                  lower
Double_t highErrNsdAlice[binsNsdAlice] = {0};          //                  upper
Double_t relErrNsdAlice[binsNsdAlice] = {0};           //                  relative
Double_t err2NsdAlice[binsNsdAlice] = {0};             //      stat+syst added quadratically
Double_t lowErr2NsdAlice[binsNsdAlice] = {0};          //                  lower
Double_t highErr2NsdAlice[binsNsdAlice] = {0};         //                  upper
Double_t relErr2NsdAlice[binsNsdAlice] = {0};          //                  relative

//
// define arrays to store ALICE INEL data points
//  
                                                // for INEL
Double_t ptInelAlice[binsInelAlice] = {0};          // pt   value (weighted)
Double_t centerPtInelAlice[binsInelAlice] = {0};    //      center of bin
Double_t widthPtInelAlice[binsInelAlice] = {0};     //      width
Double_t errPtInelAlice[binsInelAlice] = {0};       //      error (half width)
Double_t lowPtInelAlice[binsInelAlice] = {0};       //      lower edge
Double_t highPtInelAlice[binsInelAlice] = {0};      //      upper edge

Double_t inelAlice[binsInelAlice] = {0};                // INEL
Double_t statInelAlice[binsInelAlice] = {0};            //      stat error
Double_t lowStatInelAlice[binsInelAlice] = {0};         //                  lower
Double_t highStatInelAlice[binsInelAlice] = {0};        //                  upper
Double_t relStatInelAlice[binsInelAlice] = {0};         //                  relative
Double_t systInelAlice[binsInelAlice] = {0};            //      syst error
Double_t lowSystInelAlice[binsInelAlice] = {0};         //                  lower
Double_t highSystInelAlice[binsInelAlice] = {0};        //                  upper
Double_t relSystInelAlice[binsInelAlice] = {0};         //                  relative
Double_t errInelAlice[binsInelAlice] = {0};             //      stat+syst added linearly
Double_t lowErrInelAlice[binsInelAlice] = {0};          //                  lower
Double_t highErrInelAlice[binsInelAlice] = {0};         //                  upper
Double_t relErrInelAlice[binsInelAlice] = {0};          //                  relative
Double_t err2InelAlice[binsInelAlice] = {0};            //      stat+syst added quadratically
Double_t lowErr2InelAlice[binsInelAlice] = {0};         //                  lower
Double_t highErr2InelAlice[binsInelAlice] = {0};        //                  upper
Double_t relErr2InelAlice[binsInelAlice] = {0};         //                  relative

//
// define arrays to store ALICE YIELD data points
//  
                                                // for YIELD
Double_t ptYieldAlice[binsYieldAlice] = {0};        // pt   value (weighted)
Double_t centerPtYieldAlice[binsYieldAlice] = {0};  //      center of bin
Double_t widthPtYieldAlice[binsYieldAlice] = {0};   //      width
Double_t errPtYieldAlice[binsYieldAlice] = {0};     //      error (half width)
Double_t lowPtYieldAlice[binsYieldAlice] = {0};     //      lower edge
Double_t highPtYieldAlice[binsYieldAlice] = {0};    //      upper edge

Double_t yieldAlice[binsYieldAlice] = {0};               // Invariant Yield
Double_t statYieldAlice[binsYieldAlice] = {0};           //      stat error
Double_t lowStatYieldAlice[binsYieldAlice] = {0};        //                  lower
Double_t highStatYieldAlice[binsYieldAlice] = {0};       //                  upper
Double_t relStatYieldAlice[binsYieldAlice] = {0};        //                  relative
Double_t systYieldAlice[binsYieldAlice] = {0};           //      syst error
Double_t lowSystYieldAlice[binsYieldAlice] = {0};        //                  lower
Double_t highSystYieldAlice[binsYieldAlice] = {0};       //                  upper
Double_t relSystYieldAlice[binsYieldAlice] = {0};        //                  relative
Double_t errYieldAlice[binsYieldAlice] = {0};            //      stat+syst added linearly
Double_t lowErrYieldAlice[binsYieldAlice] = {0};         //                  lower
Double_t highErrYieldAlice[binsYieldAlice] = {0};        //                  upper
Double_t relErrYieldAlice[binsYieldAlice] = {0};         //                  relative
Double_t err2YieldAlice[binsYieldAlice] = {0};           //      stat+syst added quadratically
Double_t lowErr2YieldAlice[binsYieldAlice] = {0};        //                  lower
Double_t highErr2YieldAlice[binsYieldAlice] = {0};       //                  upper
Double_t relErr2YieldAlice[binsYieldAlice] = {0};        //                  relative

//
// used for fitting and shifting ALICE data
//

Double_t ptAliceFit[binsAlice] = {0};
Double_t widthPtAliceFit[binsAlice] = {0};
Double_t errPtAliceFit[binsAlice] = {0};
Double_t lowPtAliceFit[binsAlice] = {0};
Double_t highPtAliceFit[binsAlice] = {0};
Double_t nsdAliceFit[binsAlice] = {0};
Double_t errNsdAliceFit[binsAlice] = {0};
Double_t statNsdAliceFit[binsAlice] = {0};
Double_t systNsdAliceFit[binsAlice] = {0};

Double_t pt2PiPtAliceFit[binsAlice] = {0};
Double_t widthPt2PiPtAliceFit[binsAlice] = {0};
Double_t errPt2PiPtAliceFit[binsAlice] = {0};
Double_t lowPt2PiPtAliceFit[binsAlice] = {0};
Double_t highPt2PiPtAliceFit[binsAlice] = {0};
Double_t nsd2PiPtAliceFit[binsAlice] = {0};
Double_t errNsd2PiPtAliceFit[binsAlice] = {0};
Double_t statNsd2PiPtAliceFit[binsAlice] = {0};
Double_t systNsd2PiPtAliceFit[binsAlice] = {0};


//
// define arrays to store ATLAS data points
//

Double_t ptAtlas[binsAtlas] = {0};              // pt   value (weighted)
Double_t centerPtAtlas[binsAtlas] = {0};        //      center of bin
Double_t widthPtAtlas[binsAtlas] = {0};         //      width
Double_t errPtAtlas[binsAtlas] = {0};           //      error (half width)
Double_t lowErrPtAtlas[binsAtlas] = {0};        //      lower error
Double_t highErrPtAtlas[binsAtlas] = {0};       //      higher error
Double_t lowPtAtlas[binsAtlas] = {0};           //      lower edge
Double_t highPtAtlas[binsAtlas] = {0};          //      upper edge

Double_t nsdAtlas[binsAtlas] = {0};             // NSD
Double_t statNsdAtlas[binsAtlas] = {0};         //      stat error
Double_t lowStatNsdAtlas[binsAtlas] = {0};      //                  lower
Double_t highStatNsdAtlas[binsAtlas] = {0};     //                  upper
Double_t relStatNsdAtlas[binsAtlas] = {0};      //                  relative
Double_t systNsdAtlas[binsAtlas] = {0};         //      syst error
Double_t lowSystNsdAtlas[binsAtlas] = {0};      //                  lower
Double_t highSystNsdAtlas[binsAtlas] = {0};     //                  upper
Double_t relSystNsdAtlas[binsAtlas] = {0};      //                  relative
Double_t errNsdAtlas[binsAtlas] = {0};          //      stat+syst added linearly
Double_t lowErrNsdAtlas[binsAtlas] = {0};       //                  lower
Double_t highErrNsdAtlas[binsAtlas] = {0};      //                  upper
Double_t relErrNsdAtlas[binsAtlas] = {0};       //                  relative
Double_t err2NsdAtlas[binsAtlas] = {0};         //      stat+syst added quadratically
Double_t lowErr2NsdAtlas[binsAtlas] = {0};      //                  lower
Double_t highErr2NsdAtlas[binsAtlas] = {0};     //                  upper
Double_t relErr2NsdAtlas[binsAtlas] = {0};      //                  relative

//
// define arrays to store CMS data points
//

Double_t ptCms[binsCms] = {0};                  // pt   value (weighted)
Double_t centerPtCms[binsCms] = {0};            //      center of bin
Double_t widthPtCms[binsCms] = {0};             //      width
Double_t errPtCms[binsCms] = {0};               //      error (half width)
Double_t lowPtCms[binsCms] = {0};               //      lower edge
Double_t highPtCms[binsCms] = {0};              //      upper edge

Double_t nsdCms[binsCms] = {0};                 // NSD
/*
Double_t statNsdCms[binsCms] = {0};             //      stat error
Double_t lowStatNsdCms[binsCms] = {0};          //                  lower
Double_t highStatNsdCms[binsCms] = {0};         //                  upper
Double_t relStatNsdCms[binsCms] = {0};          //                  relative
Double_t systNsdCms[binsCms] = {0};             //      syst error
Double_t lowSystNsdCms[binsCms] = {0};          //                  lower
Double_t highSystNsdCms[binsCms] = {0};         //                  upper
Double_t relSystNsdCms[binsCms] = {0};          //                  relative
*/
Double_t errNsdCms[binsCms] = {0};              //      stat+syst added linearly
Double_t lowErrNsdCms[binsCms] = {0};           //                  lower
Double_t highErrNsdCms[binsCms] = {0};          //                  upper
Double_t relErrNsdCms[binsCms] = {0};           //                  relative
Double_t err2NsdCms[binsCms] = {0};             //      stat+syst added quadratically
Double_t lowErr2NsdCms[binsCms] = {0};          //                  lower
Double_t highErr2NsdCms[binsCms] = {0};         //                  upper
Double_t relErr2NsdCms[binsCms] = {0};          //                  relative

//
// define arrays to store UA1 data points
//

Double_t ptUa1[binsUa1] = {0};                  // pt   value (weighted)
Double_t centerPtUa1[binsUa1] = {0};            //      center of bin
Double_t widthPtUa1[binsUa1] = {0};             //      width
Double_t errPtUa1[binsUa1] = {0};               //      error (half width)
Double_t lowPtUa1[binsUa1] = {0};               //      lower edge
Double_t highPtUa1[binsUa1] = {0};              //      upper edge

Double_t crossUa1[binsUa1] = {0};               // Invariant cross section
/*
Double_t statCrossUa1[binsUa1] = {0};           //      stat error
Double_t lowStatCrossUa1[binsUa1] = {0};        //                  lower
Double_t highStatCrossUa1[binsUa1] = {0};       //                  upper
Double_t relStatCrossUa1[binsUa1] = {0};        //                  relative
Double_t systCrossUa1[binsUa1] = {0};           //      syst error
Double_t lowSystCrossUa1[binsUa1] = {0};        //                  lower
Double_t highSystCrossUa1[binsUa1] = {0};       //                  upper
Double_t relSystCrossUa1[binsUa1] = {0};        //                  relative
*/
Double_t errCrossUa1[binsUa1] = {0};            //      stat+syst added linearly
Double_t lowErrCrossUa1[binsUa1] = {0};         //                  lower
Double_t highErrCrossUa1[binsUa1] = {0};        //                  upper
Double_t relErrCrossUa1[binsUa1] = {0};         //                  relative
Double_t err2CrossUa1[binsUa1] = {0};           //      stat+syst added quadratically
Double_t lowErr2CrossUa1[binsUa1] = {0};        //                  lower
Double_t highErr2CrossUa1[binsUa1] = {0};       //                  upper
Double_t relErr2CrossUa1[binsUa1] = {0};        //                  relative

Double_t yieldUa1[binsUa1] = {0};               // Invariant Yield
/*
Double_t statYieldUa1[binsUa1] = {0};           //      stat error
Double_t lowStatYieldUa1[binsUa1] = {0};        //                  lower
Double_t highStatYieldUa1[binsUa1] = {0};       //                  upper
Double_t relStatYieldUa1[binsUa1] = {0};        //                  relative
Double_t systYieldUa1[binsUa1] = {0};           //      syst error
Double_t lowSystYieldUa1[binsUa1] = {0};        //                  lower
Double_t highSystYieldUa1[binsUa1] = {0};       //                  upper
Double_t relSystYieldUa1[binsUa1] = {0};        //                  relative
*/
Double_t errYieldUa1[binsUa1] = {0};            //      stat+syst added linearly
Double_t lowErrYieldUa1[binsUa1] = {0};         //                  lower
Double_t highErrYieldUa1[binsUa1] = {0};        //                  upper
Double_t relErrYieldUa1[binsUa1] = {0};         //                  relative
Double_t err2YieldUa1[binsUa1] = {0};           //      stat+syst added quadratically
Double_t lowErr2YieldUa1[binsUa1] = {0};        //                  lower
Double_t highErr2YieldUa1[binsUa1] = {0};       //                  upper
Double_t relErr2YieldUa1[binsUa1] = {0};        //                  relative

//
// define arrays to store PHOJET data points
//

Double_t nEventsPhojet = 0;

Double_t centerPtPhojet[binsPhojet] = {0};     
Double_t ptPhojet[binsPhojet] =  {0};
Double_t widthPtPhojet[binsPhojet] =  {0};
Double_t errPtPhojet[binsPhojet] =  {0};
Double_t lowPtPhojet[binsPhojet] =  {0};
Double_t highPtPhojet[binsPhojet] =  {0};

Double_t inelPhojet[binsPhojet] = {0};
Double_t errInelPhojet[binsPhojet] = {0}; 
Double_t lowErrInelPhojet[binsPhojet] = {0};
Double_t highErrInelPhojet[binsPhojet] = {0};
Double_t relErrInelPhojet[binsPhojet] = {0};

Double_t centerPt2PiPtPhojet[binsPhojet] = {0};     
Double_t pt2PiPtPhojet[binsPhojet] =  {0};
Double_t widthPt2PiPtPhojet[binsPhojet] =  {0};
Double_t errPt2PiPtPhojet[binsPhojet] =  {0};
Double_t lowPt2PiPtPhojet[binsPhojet] =  {0};
Double_t highPt2PiPtPhojet[binsPhojet] =  {0};

Double_t inel2PiPtPhojet[binsPhojet] = {0};
Double_t errInel2PiPtPhojet[binsPhojet] = {0}; 
Double_t lowErrInel2PiPtPhojet[binsPhojet] = {0};
Double_t highErrInel2PiPtPhojet[binsPhojet] = {0};
Double_t relErrInel2PiPtPhojet[binsPhojet] = {0};

//
// define arrays to store PYTHIA D6T (109) data points
//

Double_t nEventsPythia109 = 0;

Double_t centerPtPythia109[binsPythia109] = {0};     
Double_t ptPythia109[binsPythia109] =  {0};
Double_t widthPtPythia109[binsPythia109] =  {0};
Double_t errPtPythia109[binsPythia109] =  {0};
Double_t lowPtPythia109[binsPythia109] =  {0};
Double_t highPtPythia109[binsPythia109] =  {0};

Double_t inelPythia109[binsPythia109] = {0};
Double_t errInelPythia109[binsPythia109] = {0}; 
Double_t lowErrInelPythia109[binsPythia109] = {0};
Double_t highErrInelPythia109[binsPythia109] = {0};
Double_t relErrInelPythia109[binsPythia109] = {0};

Double_t centerPt2PiPtPythia109[binsPythia109] = {0};     
Double_t pt2PiPtPythia109[binsPythia109] =  {0};
Double_t widthPt2PiPtPythia109[binsPythia109] =  {0};
Double_t errPt2PiPtPythia109[binsPythia109] =  {0};
Double_t lowPt2PiPtPythia109[binsPythia109] =  {0};
Double_t highPt2PiPtPythia109[binsPythia109] =  {0};

Double_t inel2PiPtPythia109[binsPythia109] = {0};
Double_t errInel2PiPtPythia109[binsPythia109] = {0}; 
Double_t lowErrInel2PiPtPythia109[binsPythia109] = {0};
Double_t highErrInel2PiPtPythia109[binsPythia109] = {0};
Double_t relErrInel2PiPtPythia109[binsPythia109] = {0};


//
// define arrays to store PYTHIA ATLAS-CSC (306) data points
//

Double_t nEventsPythia306 = 0;

Double_t centerPtPythia306[binsPythia306] = {0};     
Double_t ptPythia306[binsPythia306] =  {0};
Double_t widthPtPythia306[binsPythia306] =  {0};
Double_t errPtPythia306[binsPythia306] =  {0};
Double_t lowPtPythia306[binsPythia306] =  {0};
Double_t highPtPythia306[binsPythia306] =  {0};

Double_t inelPythia306[binsPythia306] = {0};
Double_t errInelPythia306[binsPythia306] = {0}; 
Double_t lowErrInelPythia306[binsPythia306] = {0};
Double_t highErrInelPythia306[binsPythia306] = {0};
Double_t relErrInelPythia306[binsPythia306] = {0};

Double_t centerPt2PiPtPythia306[binsPythia306] = {0};     
Double_t pt2PiPtPythia306[binsPythia306] =  {0};
Double_t widthPt2PiPtPythia306[binsPythia306] =  {0};
Double_t errPt2PiPtPythia306[binsPythia306] =  {0};
Double_t lowPt2PiPtPythia306[binsPythia306] =  {0};
Double_t highPt2PiPtPythia306[binsPythia306] =  {0};

Double_t inel2PiPtPythia306[binsPythia306] = {0};
Double_t errInel2PiPtPythia306[binsPythia306] = {0}; 
Double_t lowErrInel2PiPtPythia306[binsPythia306] = {0};
Double_t highErrInel2PiPtPythia306[binsPythia306] = {0};
Double_t relErrInel2PiPtPythia306[binsPythia306] = {0};

//
// define arrays to store PYTHIA PERUGIA0 (320) data points
//

Double_t nEventsPythia320 = 0;

Double_t centerPtPythia320[binsPythia320] = {0};     
Double_t ptPythia320[binsPythia320] =  {0};
Double_t widthPtPythia320[binsPythia320] =  {0};
Double_t errPtPythia320[binsPythia320] =  {0};
Double_t lowPtPythia320[binsPythia320] =  {0};
Double_t highPtPythia320[binsPythia320] =  {0};

Double_t inelPythia320[binsPythia320] = {0};
Double_t errInelPythia320[binsPythia320] = {0}; 
Double_t lowErrInelPythia320[binsPythia320] = {0};
Double_t highErrInelPythia320[binsPythia320] = {0};
Double_t relErrInelPythia320[binsPythia320] = {0};

Double_t centerPt2PiPtPythia320[binsPythia320] = {0};     
Double_t pt2PiPtPythia320[binsPythia320] =  {0};
Double_t widthPt2PiPtPythia320[binsPythia320] =  {0};
Double_t errPt2PiPtPythia320[binsPythia320] =  {0};
Double_t lowPt2PiPtPythia320[binsPythia320] =  {0};
Double_t highPt2PiPtPythia320[binsPythia320] =  {0};

Double_t inel2PiPtPythia320[binsPythia320] = {0};
Double_t errInel2PiPtPythia320[binsPythia320] = {0}; 
Double_t lowErrInel2PiPtPythia320[binsPythia320] = {0};
Double_t highErrInel2PiPtPythia320[binsPythia320] = {0};
Double_t relErrInel2PiPtPythia320[binsPythia320] = {0};

