# Example to test jsroot web server locally


* JSROOT usage for the QA web page 
  * API https://root.cern.ch/js/latest/api.htm
  * Examples: https://root.cern.ch/js/latest/examples.htm
* nginx installation  
  
### Other usefull html links    
*  Monalisa - raw run details query - used on mouse over
   * http://alimonitor.cern.ch/raw/rawrun_details.jsp?run=280897
* multi-line tooltip using title tag  - can be used for bitmask explanation
  * http://jsfiddle.net/rzea/vsp6840b/3/
* usage of datalist:
  *  https://www.w3schools.com/tags/tryit.asp?filename=tryhtml5_datalist
  *  https://www.w3schools.com/tags/tryit.asp?filename=tryhtml5_details

  
  
### Install nginx
* Needed in order to test server functionality on local computer
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
  * local jsroot - http root file
    ````
    file:///data/jsroot/index.htm?file=http://aliqatpc.web.cern.ch/aliqatpc/histograms/report.root
    ````
    
### QA consideration
 
### JSROOT problems
* marker styles 20 missing - looks like there are replaced by defaults
  * other markers have differnt definiation in jsroot and root (css style problem ?)
  * checking jsroot available marker style in interactive session (jsroot 5.3.0 ???)
  * Problem of marker size ? trying to avoid it increasing marker size from0.75->1.
  * checking source code   - JSRootPainter.js
* graph lines not shown (after some changes probelm dissapeared) 
  * http://jsroot.gsi.de/5.3.0/ - lines shown 
* graph label problems for the y axis  

version                     | y axis                 | y axis label
----------------------------|------------------------|----------------
| root                      | left and right         | shown 
|http://jsroot.gsi.de/dev/  | shown only right       | shown numbers instead of labels
|http://jsroot.gsi.de/5.3.0/| shown only right       | shown numbers instead of labels
|http://jsroot.gsi.de/5.2.4 | shown only left        | shown numbers instead of labels
 
* is the problem related to the names of graphs (not unique?) 
 