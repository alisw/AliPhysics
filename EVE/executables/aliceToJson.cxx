#include "AliExternalFormatConverter.h"

int8_t PROPER_ARGS_COUNT = 9;
std::string JSON = "json";
std::string XML = "xml";


int main(int argc, char* argv[]){
    if (argc != PROPER_ARGS_COUNT){
        std::cout << "Usage:\naliconvertertest -e <event_numer> -i <input_path> -o <output_path> -f <format(xml|json)>" << std::endl;
        std::cout << "Usage:\naliconvertertest -event <event_numer> -input <input_path> -output <output_path> -format <format(xml|json)>" << std::endl;
        return(EXIT_FAILURE);
    }
    Int_t entry;
    std::string char_format;
    TString dirPath;
    TString path;

    for (int8_t i = 1; i < argc - 1; i++){
        std::string arg = argv[i];
        if (arg == "-e" || arg == "-event")
            entry = atoi(argv[i + 1]);
        else if (arg == "-i" || arg == "-input")
            dirPath = argv[i + 1];
        else if (arg == "-o" || arg == "-output")
            path = argv[i + 1];
        else if (arg == "-f" || arg == "-format")
            char_format = argv[i + 1];
    }
    AliExternalFormatConverter Converter(dirPath);
    if (char_format == JSON)
        Converter.WriteJSONToFile(path, entry);
    else if (char_format == XML)
        Converter.WriteXMLToFile(path, entry);

    return(EXIT_SUCCESS);
}
