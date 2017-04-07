void run_FB_analysis()
{
    gROOT->ProcessLine( ".L SupplementaryClasses.cxx+");
    gROOT->ProcessLine( ".L analyse_FB_TREE.C+");

    const char *dirBase = "<path_to_output_file_with_a_tree>";

    int fileId = -1;
    int ineff = 0;
    int ptW = 0;

    analyse_FB_TREE( dirBase, fileId, -1, ineff, i, ptW );
}
