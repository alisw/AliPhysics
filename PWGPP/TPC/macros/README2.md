# Example to test web server locally

* nginx installation
* JSROOT usage for the QA web page 
  * API https://root.cern.ch/js/latest/api.htm
  * Examples: https://root.cern.ch/js/latest/examples.htm
  
## Install nginx
* **Setting up my nginx server**
    http://nginx.org/en/docs/beginners_guide.html
* **Edit your nginx configuration file**
```
    sudo xemacs /etc/nginx/nginx.conf
``` 
  * add browsable directory /data/alice-tpc-notes/ and executable jsroot
````
    server {
                listen     localhost:90;
                location /data/alice-tpc-notes/ {
                root /;
                autoindex on;
                }
                location /jsoot/ {
                root /data/;
                }
        }
````

* **Test nginx configuration after modification**
``
    sudo nginx -t
``
* **Restart nginx**
````
  sudo service nginx restart
  or
  sudo nginx -s reload
````
* **check nginx  log files**
  * in the /var/log/nginx/error.log or access log
````
ls -alrt /var/log/nginx
````
  * in case of problem check http error code 
    https://en.wikipedia.org/wiki/List_of_HTTP_status_codes#nginx
* **setup directory permission to be accessible by www-data**
  * create new user or add existing user to www-data group
  * see https://www.digitalocean.com/community/questions/proper-permissions-for-web-server-s-directory
```
 # add existing user (miranov) to the group
 sudo usermod -a -G www-data  miranov
 # check  groups
 id miranov
```

## Browse web page

* Browse my notes data directory
````
http://localhost:90/data/alice-tpc-notes/JIRA/
````
* Load jsroot browser
  * reading data from remote file - OK
    ````
    http://localhost:90/data/jsroot/index.html?layout=grid1x2&file=http://aliqatpc.web.cern.ch/aliqatpc/histograms/report.root
    ````
  * reading data from local file absolute path
    * in one of the test example did not work properly in case file not properly closed
    * not obvious error message in browser
    ````
    http://localhost:90/data/jsroot/index.htm?layout=grid1x2file=http://localhost:90/data/alice-tpc-notes/JIRA/ATO-83/test/report.root
    ````
  * reading data from local file absolute path
    ````
    http://localhost:90/data/jsroot/index.htm?layout=grid1x2&file=../alice-tpc-notes/JIRA/ATO-83/test/report.root
    ````
  * locat jsroot - http root file
    ````
    file:///data/jsroot/index.htm?file=http://aliqatpc.web.cern.ch/aliqatpc/histograms/report.root
    ````
    
## QA consideration
 