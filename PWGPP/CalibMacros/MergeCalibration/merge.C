void merge(const char* outputDir){
	gROOT->LoadMacro("mergeCalibObjects.C");
	IterAlien(outputDir);
	return;
}
