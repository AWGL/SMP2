var nano = process.argv.slice(2);

// Load node modules
var request = require(nano + '/node_modules/request');
var fs = require(nano + '/node_modules/fs');
// Load configuration files
var config = require('./config.json');
var runconfig = require('./runConfig.json')

var APISERVER = config.apiServer;
var APIVERSION = config.apiVersion;
//Set accessToken variable from the config file
var ACCESSTOKEN = config.accessToken;
//Load run-specific config files
var NUMPAIRS = runconfig.numPairs;
var PROJECTID = runconfig.projectID;
var NEGATIVECONTROL = runconfig.negativeControl;

// Obtain time at which script was launched to enable later timeout
var STARTTIME = new Date().getTime();

// Global variables required for re-calling outer iteration
var J;
var APPRES;

// Variables which may need adjusting
var TEMPLATE = "SMP2_CRUK_V2_03.15.xlsx"; //Update manually if it changes

// Variables- adjust these to the desired intervals for polling and timeout of the script
var POLLINGINTERVAL = 60000; //Change to 60000 for live
var TIMEOUT = 7200000; // 60000 is 1 minute // 7200000 is 2 hours

//temp vars
var fileID = config.get("fileIDexample");
var appResultID = config.get("appResultIDexample");

//Access appResults through projectid
//This is asynchronous- need to put in a callback to ensure that we can access the data
function appResultsByProject(cb){
    request.get(
        APISERVER + APIVERSION + "/projects/" + PROJECTID + "/appresults",
        {qs: {"access_token": ACCESSTOKEN}},
        function (error, response, body) {
            if (!error && response.statusCode === 200) {
                var projectAppResults = JSON.parse(body);
                return cb(null, projectAppResults, console.log("App results successfully retrieved"));
            }
            else if (response.statusCode !== 200) {
                return cb(Error('Response status is ' + response.statusCode + " " + body));
            }
            else if (error) {
                return cb(Error(error.message));
            }
        }
    );
}

function checkAppResultsComplete(appResults, refresh, cb) {
    console.log("Checking status of app results");
    var numComplete = 0;
    var appResultsLen = appResults.Response.Items.length;
    var appResultsArr = [];
    // See the status of all of the appSessions
    for (i = 0; i < appResultsLen; i++) {
        if (appResults.Response.Items[i].Status === "Complete") {
            numComplete += 1;
            // Store the appResults IDs which are needed for downloading the files
            appResultsArr[i] = appResults.Response.Items[i].Id;
        }
    }
    //Stop execution of the polling function after a certain time has elapsed (assume the process has failed after this time)
    if (new Date().getTime() - STARTTIME > TIMEOUT) {
        clearInterval(refresh);
        return cb(Error("Polling timed out"));
    }
    else if (appResultsLen === NUMPAIRS && numComplete === NUMPAIRS) {
        clearInterval(refresh);
        console.log("All appSessions complete");
        //In here want to call another function to kick off getting appresults, getting fileids and download of results
        return cb(null, appResultsArr);
    }
}

function iterator(appRes, j, cb) {
    var appResId = appRes[j];
    if (appRes.length === j) {
        return (console.log("Files retrieved"))
    }
    getFileIds(appResId, function(err, fileIds) {
        if (err) {
            return cb(Error(err));
        }else {
            iterFileId(fileIds, 0, j);
        }
    });
}

function iterFileId(appResFiles, i) {
    numFiles = appResFiles.Response.Items.length;
    if (i === (numFiles-1)) {
        J+=1;
        return iterator(APPRES, J);
    }
    if (i < (numFiles-1)) {
        var fileId = appResFiles.Response.Items[i].Id;
        var fileName = appResFiles.Response.Items[i].Name;
        if (fileName !== TEMPLATE && fileName !== NEGATIVECONTROL + ".bam") {
            downloadFile(fileId, fileName, function(err, result){
                if (err) {
                    throw new Error(console.log("File download failed ") + err);
                }else {
                    console.log(result);
                    iterFileId(appResFiles, i+1);}});
        }
    }
}

// Get file IDs- to retrieve xlsx, bam and bai files only
function getFileIds(appResultId, cb) {
    console.log("Getting file Ids for " + appResultId);
    request.get(
        APISERVER + APIVERSION + "/appresults/" + appResultId + "/files?SortBy=Id&Extensions=.xlsx,.bai&Offset=0&Limit=50&SortDir=Asc",
        {qs: {"access_token": ACCESSTOKEN}},
        function (error, response, body) {
            if (!error && response.statusCode === 200) {
                var appResultFiles = JSON.parse(body);
                return cb(null, appResultFiles);
            }
            else if (response.statusCode !== 200) {
                return cb(Error('Response status is ' + response.statusCode + " " + body));
            }
            else if (error) {
                return cb(Error(error.message));
            }
        }
    );
}

// Download files
function downloadFile(fileIdentifier, outFile, cb) {
    var writeFile = fs.createWriteStream(outFile);
    request.get(
        APISERVER + APIVERSION + "/files/" + fileIdentifier + "/content",
        {qs: {"access_token": ACCESSTOKEN}}).on('error', function(err) {cb(Error(err))})
        .pipe(writeFile).on('close', function(){cb(null, "Download Success " + outFile)});
s
}

// Repeatedly call the function to check if the results are complete or not
function poll(){
    var refresh = setInterval(function(){
        appResultsByProject(function(err, appRes){
            if (err) return console.log(err);
            checkAppResultsComplete(appRes, refresh, function(err, appResIds) {
                if (err) return console.log(err);
                iterator(APPRES=appResIds, J=0, function(output){
                    console.log(output)});
            });
        });
    }, POLLINGINTERVAL);
}


// Call function
poll();
