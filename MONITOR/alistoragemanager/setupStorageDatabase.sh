#!/bin/bash
HOST="localhost"   # IP of machine on which mysql database is located
PORT="123"
DATABASE="database"
USER="user"
PASS="pass123"
TABLE="table"
STORAGE_PATH="/some/path/"
MAX_SIZE="300000000"
MAX_OCCUPATION="100"
REMOVE_PERCENT="50"
EVENTS_IN_FILE="20"
EVENT_SERVER="localhost"   # IP of machine running alieventserver
EVENT_SERVER_USER="user"
EVENT_SERVER_PORT="124"
STORAGE_SERVER="localhost"      # IP of machine running alistorage
STORAGE_SERVER_USER="user"
STORAGE_SERVER_PORT="125"      # server thread communication port
STORAGE_CLIENT_PORT="126"      # client thread communication port
XML_SERVER_PORT="127"          # server of xml files
mysql -u root -ppassword -e "create database if not exists $DATABASE;"
mysql -u root -ppassword -e "grant ALL PRIVILEGES on $DATABASE.* to '$USER'@'$HOST' identified by '$PASS';"
mysql -u root -ppassword -e "use $DATABASE;"
mysql -u root -ppassword -e "CREATE TABLE IF NOT EXISTS $DATABASE.$TABLE(\
run_number int(6) NOT NULL,\
event_number int(6) NOT NULL,\
system text(7) DEFAULT NULL,\
multiplicity int(5) DEFAULT NULL,\
permanent tinyint(1) DEFAULT NULL,\
file_path text(100) DEFAULT NULL,\
PRIMARY KEY(run_number,event_number));"

echo "-----------------------------"
echo "Databases successfuly created"
echo "-----------------------------"
