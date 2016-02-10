#include "AliExternalFormatConverter.h"


int main(int argv, char* argc[])
{
    if (argv != 4){
        std::cout << "Usage:\naliceToJson <event_numer> <input_path> <output_path>" << std::endl;
        return(EXIT_FAILURE);
    }
    Int_t entry = atoi(argc[1]); // Event number
    const TString dirPath = argc[2]; // Path to the directory with ESDs
    const TString path = argc[3]; // Output file

    AliExternalFormatConverter Converter(dirPath);
    //Converter.LoadFiles(filePath, filenameFriends);
    // Converter.WriteCSVToFile(datFilePath, entry);
//    Converter.SerializeAllEvents(path);
//    Converter.WriteXMLToFile(XMLPath, entry);
    Converter.WriteJSONToFile(path, entry);
    return(EXIT_SUCCESS);
}
