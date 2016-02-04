## Licensed under the Apache License, Version 2.0 (the 'License');
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an 'AS IS' BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.

##' Get one page of search results for a particular endpoint from
##' a GA4GH compliant server.
##'
##' In general, use higher level methods such as getReads and getVariants
##' instead.
##'
##' @param base_url the base URL of the GA4GH server to which
##' requests will be sent to. 
##' @param endpoint a endpoint such as reads, variants,
##' variantSets, etc...
##' @param method which method (GET/POST) should be used in the request to
##' the API?
##' @param body The body of the message to GET/POST to the endpoint.
##' @param query_fields The fields to be returned in the search response.
##' @param pageToken The page token. This can be NULL for the first page.
##' @return The raw response converted from JSON to an R object.
##' 
##' @examples
##' body <- list(readGroupSetIds=list("CMvnhpKTFhDnk4_9zcKO3_YB"),
##'              referenceName="22",
##'              start=16051400, end=16051500, pageToken=NULL)
##' reads <- getSearchPage("reads", body, NULL, NULL)
##' summary(reads)
##' @export
sendRequest <- function(base_url, endpoint, method,
                         body, query_fields, pageToken) {

    if (missing(base_url)) {
        stop("ERROR: baseURL not provided")
    }
    if (missing(endpoint)) {
        stop("ERROR: endpoint not provided")
    }
    if (missing(method)) {
        stop("ERROR: method not provided")
    }
    if (missing(body)) {
        stop("ERROR: body not provided")
    } else if (class(body) != "list"){
        stop("ERROR: 'body' not a list")
    } 
    if (missing(query_fields)) {
        stop("ERROR: query_fields not provided")
    }
    if (missing(pageToken)) {
        stop("ERROR: pageToken not provided")
    }
    
    if (!is.null(query_fields)) {
        if (!grepl("nextPageToken", query_fields)) {
            query_fields <- paste(query_fields, "nextPageToken",
                                  sep=",")
        }        
    }
    
    if (method == "POST") {
        request_function <- httr::POST
        endpoint_method <- "search"
    } else if (method == "GET") {
        req_function <- httr::GET
        endpoint_method <- "get"
    }

    if (stringr::str_sub(base_url, nchar(base_url)) == "/") {
        base_url <- stringr::str_sub(base_url, 1, nchar(base_url) - 1)
    }
    
    url <- paste(base_url, endpoint, endpoint_method, sep = "/")

    message("Sending ", method, " request to ", url)
    
    response <- request_function(url = url,
                                 body = body,
                                 encode = "json",
                                 httr::add_headers("Content-Type"="application/json")
                                 )

    httr::content(response)
}
