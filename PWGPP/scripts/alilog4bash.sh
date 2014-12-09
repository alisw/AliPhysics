#!/usr/bin/env bash
#
#adopted for ALICE from:
#--------------------------------------------------------------------------------------------------
# log4bash - Makes logging in Bash scripting suck less
# Copyright (c) Fred Palmer
# Licensed under the MIT license
# http://github.com/fredpalmer/log4bash
#--------------------------------------------------------------------------------------------------

# Useful global variables that users may wish to reference
ALILOG_SCRIPT_ARGS="$@"
ALILOG_SCRIPT_NAME="$0"
ALILOG_SCRIPT_NAME="${SCRIPT_NAME#\./}"
ALILOG_SCRIPT_NAME="${SCRIPT_NAME##/*/}"
ALILOG_SCRIPT_BASE_DIR="$(cd "$( dirname "$0")" && pwd )"

# This should probably be the right way - didn't have time to experiment though
# declare -r ALILOG_INTERACTIVE_MODE="$([ tty --silent ] && echo on || echo off)"
# declare -r ALILOG_INTERACTIVE_MODE=$([ "$(uname)" == "Darwin" ] && echo "on" || echo "off")
declare -r ALILOG_INTERACTIVE_MODE="off"

#--------------------------------------------------------------------------------------------------
# Begin Help Section

ALILOG_HELP_TEXT=""

# This function is called in the event of an error.
# Scripts which source this script may override by defining their own "alilog_usage" function
alilog_usage() {
    echo -e "${HELP_TEXT}";
    return 1;
}

# End Help Section
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
# Begin Logging Section
if [[ "${ALILOG_INTERACTIVE_MODE}" == "off" ]]
then
    # Then we don't care about alilog colors
    declare -r ALILOG_DEFAULT_COLOR=""
    declare -r ALILOG_ERROR_COLOR=""
    declare -r ALILOG_INFO_COLOR=""
    declare -r ALILOG_SUCCESS_COLOR=""
    declare -r ALILOG_WARN_COLOR=""
    declare -r ALILOG_DEBUG_COLOR=""
else
    declare -r ALILOG_DEFAULT_COLOR="\033[0m"
    declare -r ALILOG_ERROR_COLOR="\033[1;31m"
    declare -r ALILOG_INFO_COLOR="\033[1m"
    declare -r ALILOG_SUCCESS_COLOR="\033[1;32m"
    declare -r ALILOG_WARN_COLOR="\033[1;33m"
    declare -r ALILOG_DEBUG_COLOR="\033[1;34m"
fi

# This function scrubs the output of any control characters used in colorized output
# It's designed to be piped through with text that needs scrubbing.  The scrubbed
# text will come out the other side!
prepare_alilog_for_nonterminal() {
    # Essentially this strips all the control characters for alilog colors
    sed "s/[[:cntrl:]]\[[0-9;]*m//g"
}

alilog() {
    local alilog_text="$1"
    local alilog_level="$2"
    local alilog_color="$3"

    # Default level to "info"
    [[ -z ${alilog_level} ]] && alilog_level="INFO";
    [[ -z ${alilog_color} ]] && alilog_color="${ALILOG_INFO_COLOR}";

    echo -e "${alilog_color}[$(date +"%Y-%m-%d %H:%M:%S %Z")] [${alilog_level}] ${alilog_text} ${ALILOG_DEFAULT_COLOR}";
    return 0;
}

alilog_info()      { alilog "$@"; }

alilog_speak()     {
    if type -P say >/dev/null
    then
        local easier_to_say="$1";
        case "${easier_to_say}" in
            studionowdev*)
                easier_to_say="studio now dev ${easier_to_say#studionowdev}";
                ;;
            studionow*)
                easier_to_say="studio now ${easier_to_say#studionow}";
                ;;
        esac
        say "${easier_to_say}";
    fi
    return 0;
}

alilog_success()   { alilog "$1" "SUCCESS" "${ALILOG_SUCCESS_COLOR}"; }
alilog_error()     { alilog "$1" "ERROR" "${ALILOG_ERROR_COLOR}"; alilog_speak "$1"; }
alilog_warning()   { alilog "$1" "WARNING" "${ALILOG_WARN_COLOR}"; }
alilog_debug()     { alilog "$1" "DEBUG" "${ALILOG_DEBUG_COLOR}"; }
alilog_captains()  {
    if type -P figlet >/dev/null;
    then
        figlet -f computer -w 120 "$1";
    else
        alilog "$1";
    fi
    
    alilog_speak "$1";

    return 0;
}

alilog_campfire() {
    # This function performs a campfire notification with the arguments passed to it
    if [[ -z ${CAMPFIRE_API_AUTH_TOKEN} || -z ${CAMPFIRE_NOTIFICATION_URL} ]]
    then
        alilog_warning "CAMPFIRE_API_AUTH_TOKEN and CAMPFIRE_NOTIFICATION_URL must be set in order alilog to campfire."
        return 1;
    fi

    local campfire_message="
    {
        \"message\": {
            \"type\":\"TextMessage\",
            \"body\":\"$@\"
        }
    }"

    curl                                                            \
        --write-out "\r\n"                                          \
        --user ${CAMPFIRE_API_AUTH_TOKEN}:X                         \
        --header 'Content-Type: application/json'                   \
        --data "${campfire_message}"                                \
        ${CAMPFIRE_NOTIFICATION_URL}

    return $?;
}

# End Logging Section
#--------------------------------------------------------------------------------------------------

