void testcvsKeyWords(){
    // This macro tests the cvs keywords Date and Revision
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    const Char_t *cvsDate="$Date$";
    const Char_t *cvsRevision="$Revision$";
    const Char_t *cvsId="$Id$";
    Char_t string[100];

    sprintf(string,"%s %s %s",cvsDate,cvsRevision,cvsId);
    printf("%s\n",string);
}
