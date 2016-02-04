## Copywrite 2015 Oregon Health and Science University. 
## All rights reserved. 
##
## Originally developed by Google Inc and distributed in the
## GoogleGenomics package with this licence.
##
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

##' Get one page of reads from a GA4GH compliant server.
##'
##' In general, use the searchReads method instead.  It calls this method,
##' returning reads from all of the pages that comprise the requested
##' genomic range.
##'
##'
##' Note that the Global Alliance for Genomics and Health API uses a 0-based
##' coordinate system.  For more detail, please see GA4GH discussions such
##' as the following:
##' \itemize{
##'    \item\url{https://github.com/ga4gh/schemas/issues/168}
##'    \item\url{https://github.com/ga4gh/schemas/issues/121}
##' }
##'
##' @param base_url the base URL of the GA4GH server to which
##' requests will be sent to.
##' @param readGroupIds The read group set ID.
##' @param referenceId The chromosome.
##' @param start Start position on the chromosome in 0-based coordinates.
##' @param end End position on the chromosome in 0-based coordinates.
##' @param fields A subset of fields to retrieve.  The default (NULL) will
##'   return all fields.
##' @param pageToken The page token. This can be NULL (default) for the first
##'   page.
##' @return A two-element list is returned by the function.
##'
##'     reads: A list of R objects corresponding to the JSON objects returned
##'               by the GA4GH Reads API.
##'
##'     nextPageToken: The token to be used to retrieve the next page of
##'                    results, if applicable.
##' @family page fetch functions
##' @examples
##' readsPage <- searchReadsPage()
##' summary(readsPage)
##' summary(readsPage$reads[[1]])
##' @export
searchReadsPage <- function(base_url,
                         readGroupIds = list("MWtnLXAzLXN1YnNldDpIRzAwNTM0OkVSUjAyMDIzOA=="),
                         referenceId = "R1JDaDM4LXN1YnNldDoz",
                         start = NULL,
                         end = NULL,
                         fields = NULL,
                         pageToken = NULL) {

    if (class(readGroupIds) != "list") {
        readGroupIds <- list(readGroupIds)
    }
    
    body <- list(readGroupIds = readGroupIds,
                 referenceId = referenceId,
                 start = start, end = end,
                 pageToken = pageToken)

    results <- sendRequest(base_url = base_url,
                           endpoint = "reads",
                           method = "POST",
                           body = body,
                           query_fields = fields,
                           pageToken = pageToken)

    list(reads = results$alignments,
         nextPageToken = results$nextPageToken)
}

##' Get reads from a GA4GH compliant server. 
##'
##' This function will return all of the reads that comprise the requested
##' genomic range, iterating over paginated results as necessary.
##'
##'
##' Optionally pass a converter as appropriate for your use case.  By passing it
##' to this method, only the converted objects will be accumulated in memory. The
##' converter function should return an empty container of the desired type
##' if called without any arguments.
##'
##' @param base_url the base URL of the GA4GH server to which
##' requests will be sent to.
##' @param readGroupIds The read group set ID.
##' @param referenceId The chromosome.
##' @param start Start position on the chromosome in 0-based coordinates.
##' @param end End position on the chromosome in 0-based coordinates.
##' @param fields A subset of fields to retrieve.  The default (NULL) will
##'               return all fields.
##' @param converter A function that takes a list of read R objects and returns
##'                  them converted to the desired type.
##' @return By default, the return value is a list of R objects
##' corresponding to the JSON objects returned by the GA4GH
##' Reads API.  If a converter is passed, object(s) of the type
##' returned by the converter will be returned by this function.
searchReads <- function(base_url,
                     readGroupIds = list("MWtnLXAzLXN1YnNldDpIRzAwNTM0OkVSUjAyMDIzOA=="),
                     referenceId = "R1JDaDM4LXN1YnNldDoz",
                     start = NULL,
                     end = NULL,
                     fields = NULL,
                     converter = c) {
    pageToken <- NULL
    reads <- converter()
    repeat {
        result <- searchReadsPage(base_url = base_url,
                               readGroupIds = readGroupIds,
                               referenceId = referenceId,
                               start = start,
                               end = end,
                               fields = fields,
                               pageToken = pageToken)
        pageToken <- result$nextPageToken
        ## TODO improve performance,
        ## see https://github.com/Bioconductor/GoogleGenomics/issues/17
        reads <- c(reads, converter(result$reads))
        if (is.null(pageToken)) {
            break
        }
        message(paste("Continuing read query with the nextPageToken:", pageToken))
    }

    message("Reads are now available.")
    reads
}

## Transformation helpers
cigar_enum_map <- list(
    ALIGNMENT_MATCH = "M",
    CLIP_HARD = "H",
    CLIP_SOFT = "S",
    DELETE = "D",
    INSERT = "I",
    PAD = "P",
    SEQUENCE_MATCH = "=",
    SEQUENCE_MISMATCH = "X",
    SKIP = "N")

getCigar <- function(read) {
    paste(
        sapply(
            read$alignment$cigar,
            function(cigarPiece) {
                paste0(cigarPiece$operationLength,
                       cigar_enum_map[cigarPiece$operation])
            }),
        collapse = "")
}

getPosition <- function(read) {
    as.integer(read$alignment$position$position)
}

getReferenceName <- function(read) {
    read$alignment$position$referenceName
}

getFlags <- function(read) {
    flags <- 0

    if (isTRUE(read$numberReads == 2)) {
        flags <- flags + 1  # read_paired
    }
    if (isTRUE(read$properPlacement)) {
        flags <- flags + 2  # read_proper_pair
    }
    if (is.null(getPosition(read))) {
        flags <- flags + 4  # read_unmapped
    }
    if (is.null(read$nextMatePosition$position)) {
        flags <- flags + 8  # mate_unmapped
    }
    if (isTRUE(read$alignment$position$reverseStrand)) {
        flags <- flags + 16  # read_reverse_strand
    }
    if (isTRUE(read$nextMatePosition$reverseStrand)) {
        flags <- flags + 32  # mate_reverse_strand
    }
    if (isTRUE(read$readNumber == 0)) {
        flags <- flags + 64  # first_in_pair
    }
    if (isTRUE(read$readNumber == 1)) {
        flags <- flags + 128  # second_in_pair
    }
    if (isTRUE(read$secondaryAlignment)) {
        flags <- flags + 256  # secondary_alignment
    }
    if (isTRUE(read$failedVendorQualityChecks)) {
        flags <- flags + 512  # failed_quality_check
    }
    if (isTRUE(read$duplicateFragment)) {
        flags <- flags + 1024  # duplicate_read
    }
    if (isTRUE(read$supplementaryAlignment)) {
        flags <- flags + 2048  # supplementary_alignment
    }
    flags
}

##' Convert reads to GAlignments.
##'
##' Note that the Global Alliance for Genomics and Health API uses a 0-based
##' coordinate system.  For more detail, please see GA4GH discussions such
##' as the following:
##' \itemize{
##'    \item\url{https://github.com/ga4gh/schemas/issues/168}
##'    \item\url{https://github.com/ga4gh/schemas/issues/121}
##' }
##'
##' @param reads A list of R objects corresponding to the JSON objects
##'  returned by the GA4GH Reads API.
##' @param oneBasedCoord Convert genomic positions to 1-based coordinates.
##' @param slStyle The style for seqnames (chrN or N or...).  Default is UCSC.
##' @return \link[GenomicAlignments]{GAlignments}
##' @family reads converter functions
##' @examples
##' alignments1 <- searchReads(converter = readsToGAlignments)
##' summary(alignments1)
##' alignments2 <- readsToGAlignments(searchReads())
##' print(identical(alignments1, alignments2))
##' @export
readsToGAlignments <- function(reads, oneBasedCoord = TRUE,
                               slStyle = "UCSC") {

    if (missing(reads)) {
        return(GAlignments())
    }

    names <- sapply(reads, "[[", "fragmentName")
    cigars <- sapply(reads, getCigar)
    positions <- sapply(reads, getPosition)
    if (oneBasedCoord) {
        positions <- as.integer(positions + 1)
    }
    flags <- sapply(reads, getFlags)
    chromosomes <- sapply(reads, getReferenceName)

    isMinusStrand <- bamFlagAsBitMatrix(as.integer(flags),
                                        bitnames = "isMinusStrand")
    alignments <- GAlignments(
        seqnames = Rle(chromosomes),
        strand = strand(as.vector(ifelse(isMinusStrand, "-", "+"))),
        pos = positions, cigar = cigars, names = names, flag = flags)

    seqlevelsStyle(alignments) <- slStyle
    alignments
}
