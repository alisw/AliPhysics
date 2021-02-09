#!/usr/bin/python

from os import environ as env
import os
import subprocess
from sys import argv

model_to_file = {
    "KrakowSFO": "krakow.ini",
    "BlastWave": "blastwave.ini",
    "BWAVT": "bwa.ini",
    "BWAVTDelay": "bwa.ini",
    "BWAVLinear": "bwa.ini",
    "BWAVLinearDelay": "bwa.ini",
    "BWAVLinearFormation": "bwa.ini",
    "Lhyquid3D": "lhyquid3d.ini",
    "Lhyquid2DBI": "lhyquid2dbi.ini"
}

def print_usage():
    print("Configure the Therminator2 files based on env variables.")
    print("Run this script in the directory where you run Therminator2 (where events.ini and fomodel/ are).")
    print("This script replaces variables in configuration files based on the $THERM2_PARAMS_<parameter_name> environmental variables if they exist.")
    print("It also downloads and sets the proper freeze-out xml file if the $THERM2_PARAMS_XML_PATH is set.")
    print("Script usage:")
    print("%s [-verbose]" % argv[0])


def v_print(text):
    if len(argv)>1:
        if argv[1] in ["-v", "-verbose"]:
            print(text)


def apply_config_for_file(path):
    contents = None
    temp_file = "%s.tmp" % path
    with open(path, "r") as f:
        contents = f.read()
    
    parameters = [] 
    with open(temp_file, "w") as f:
        for line in contents.splitlines():
            if not line:
                continue
            if line[0] in "#[":
                continue
            
            line = "".join(line.split()) # remove all whitespace from string
            key, value = line.split("=")
            parameters.append(key)
            write_key_value(f, key, value)

    os.remove(path)
    os.rename(temp_file, path)
    return parameters
            

def write_key_value(file, key, default_val):
    val = env.get("THERM2_PARAMS_%s" % key, default_val)
    try:
        v_print("Got non-default value for %s: %s" % (key, env["THERM2_PARAMS_%s" % key]))
    except:
        pass
    file.write("%s=%s\n" % (key, val))
    

def set_param(key, val):
    env["THERM2_PARAMS_%s" % key] = str(val)


def main():

    if len(argv)>1:
        if argv[1] not in ["-v", "-verbose"]:
            print_usage()
            exit()

    set_param("EventFileType", "text")  # the parser works with text files
    set_param("EventSubDir", "./")      # easy to handle
    apply_config_for_file("events.ini")

    model = env.get("THERM2_PARAMS_FreezeOutModel", "Lhyquid2DBI")
    custom_xml = env.get("THERM2_PARAMS_XML_PATH", -1)
    if custom_xml != -1:
        v_print("Got nonzero XML_PATH: %s" % custom_xml)
        if custom_xml[:6] == "alien:":  # if it's a GRID path, download the file
            v_print("It's a GRID path, downloading it")
            subprocess.call(["alien_cp", custom_xml, "fomodel/."])
            set_param("FreezeFile", custom_xml.split("/")[-1])    # Set the path to xml file name
            v_print("Downloaded and set XML_PATH")
        else:
            set_param("FreezeFile", custom_xml)   # local path
            v_print("It's a local path. Set.")

    model_file = "fomodel/%s" % model_to_file[model]
    model_params = apply_config_for_file(model_file)
    if model in ["Lhyquid3D", "Lhyquid2DBI"]:
        if "EventSubDir" not in model_params:
            v_print("Model is hydro, setting EventSubDir manually")
            with open(model_file, "a") as f:    # by default this option is commented out in these models
                f.write("EventSubDir=./")       # so apply_config_for_file won't change it


if __name__ == "__main__":
    main()



