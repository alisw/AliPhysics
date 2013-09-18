Index: correct.C
===================================================================
--- correct.C	(revision 64235)
+++ correct.C	(working copy)
@@ -7667,6 +7667,14 @@
     Float_t leadingPtArr[] = { 2.0, 3.0, 4.0, 8.0, 15.0, 20.0 };
     Float_t assocPtArr[] =     { 0.15, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
   }
+  else if (1)
+  {
+    //Example for the Hadrons_Example wagon
+    maxLeadingPt = 3;
+    maxAssocPt = 6;
+    Float_t leadingPtArr[] =     { 3.0, 5.0, 8.0, 16.0 };
+    Float_t assocPtArr[] = {0.15, 0.3, 50.0, 0.5, 50.0, 1.0, 50.0 };
+  }
   else if (0)
   {
     //pA, trigger from all pT
@@ -7861,7 +7869,7 @@
   for (Int_t i=0; i<maxLeadingPt; i++)
     for (Int_t j=1; j<maxAssocPt; j++)
     {
-      if(1){
+      if(0){
 	if(j!=(i+1))continue;
 	Printf("\nOnly symmetric pt bins selected, leading pt: %f - %f     associated pt: %f - %f",leadingPtArr[i],leadingPtArr[i+leadingPtOffset],assocPtArr[j],assocPtArr[j+1]); 
       }
@@ -7869,6 +7877,8 @@
       gpTMin = assocPtArr[j] + 0.01;
       gpTMax = assocPtArr[j+1] - 0.01;
       
+      if(gpTMin >= gpTMax)continue;
+	
       SetupRanges(h);
       SetupRanges(hMixed);
       SetupRanges(h2);
@@ -7877,7 +7887,7 @@
       SetupRanges(hMixed3);
 //       SetupRanges(hMixed3);
 
-      if (assocPtArr[j] >= leadingPtArr[i+leadingPtOffset])
+      if(0)if (assocPtArr[j] >= leadingPtArr[i+leadingPtOffset])
 	continue;
   
       TH1* hist1 = 0;
@@ -7943,6 +7953,12 @@
       }      
       else if (1)
       {
+	// pp, MB
+	Int_t step = 8;      
+	GetSumOfRatios(h, hMixed, &hist1,  step, 0, 100, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, kTRUE, kTRUE); 
+      }      
+      else if (0)
+      {
 	// pp
 	Int_t step = 8;
       
