#include "AliExternalFormatConverter.h"


int main(int argv, char* argc[]){
    if (argv != 4){
        std::cout << "Missing arguments" << std::endl;
        return(EXIT_FAILURE);
    }
    Int_t entry = atoi(argc[1]); // Event number
    const TString dirPath = argc[2]; // Path to the directory with ESDs
    const TString JSONPath = argc[3]; // Output file

    AliExternalFormatConverter Converter(dirPath);
    //Converter.LoadFiles(filePath, filenameFriends);
    // Converter.WriteCSVToFile(datFilePath, entry);
    Converter.WriteJSONToFile(JSONPath, entry);
    //Converter.WriteXMLToFile(XMLPath, entry);
    return(EXIT_SUCCESS);
}
