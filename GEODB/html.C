{
// to run this macro, you must have the correct .rootrc file
// in your galice directory.
// The gAlice classes summary documentation go to directory html
// The gAlice classes source  documentation go to directory html/src
// The example macros documentation go to directory html/examples
   
   gSystem->Load("./libGeoAlice.so"); 
   THtml html;
   html.MakeAll(1);

}
