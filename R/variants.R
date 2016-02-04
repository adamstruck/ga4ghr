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

##' Get one page of variants from a GA4GH compliant server.
##'
##' In general, use the getVariants method instead.  It calls this method,
##' returning variants from all of the pages that comprise the requested
##' genomic range.
##'
##' @param base_url the base URL of the GA4GH server to which
##' requests will be sent to.
##' @param variantSetId The variant set ID.
##' @param referenceName The chromosome.
##' @param start Start position on the chromosome in 0-based coordinates.
##' @param end End position on the chromosome in 0-based coordinates.
##' @param fields A subset of fields to retrieve.  The default (NULL) will
##'   return all fields.
##' @param pageToken The page token. This can be NULL (default) for the first
##'   page.
##' @return A two-element list is returned by the function.
##'
##'     variants: A list of R objects corresponding to the JSON objects returned
##'               by the GA4GH Variants API.
##'
##'     nextPageToken: The token to be used to retrieve the next page of
##'                    results, if applicable.
##' @family page fetch functions
##' @examples
##' variantsPage <- getVariantsPage()
##' summary(variantsPage)
##' summary(variantsPage$variants[[1]])
##' @export
searchVariantsPage <- function(base_url,
                               variantSetId = "MWtnLXAzLXN1YnNldDptdm5jYWxs",
                               referenceName = "1",
                               start = NULL,
                               end = NULL,
                               fields = NULL,
                               pageToken = NULL) {

    body <- list(variantSetId = variantSetId,
                 referenceName = referenceName,
                 start = start, end = end,
                 pageToken = pageToken)

    results <- sendRequest(base_url = base_url,
                           endpoint = "variants",
                           method = "POST",
                           body = body,
                           query_fields = fields,
                           pageToken = pageToken)

    list(variants = results$variants,
         nextPageToken = results$nextPageToken)
}

##' Get variants from a GA4GH compliant server.
##'
##' This function will return all of the variants that comprise the requested
##' genomic range, iterating over paginated results as necessary.
##'
##'
##' Optionally pass a converter as appropriate for your use case.  By passing it
##' to this method, only the converted objects will be accumulated in memory.
##' The converter function should return an empty container of the desired type
##' if called without any arguments.
##'
##' @param base_url the base URL of the GA4GH server to which
##' requests will be sent to.
##' @param variantSetId The variant set ID.
##' @param referenceName The chromosome.
##' @param start Start position on the chromosome in 0-based coordinates.
##' @param end End position on the chromosome in 0-based coordinates.
##' @param fields A subset of fields to retrieve.  The default (NULL) will
##'               return all fields.
##' @param converter A function that takes a list of variant R objects and
##'                  returns them converted to the desired type.
##' @return By default, the return value is a list of R objects
##'   corresponding to the JSON objects returned by the GA4GH
##'   Variants API.  If a converter is passed, object(s) of the type
##'   returned by the converter will be returned by this function.
##' @seealso \code{\link{searchReads}} for equivalent function for reads, and
##'   \code{\link{variantsToVRanges}} for a converter function.
##' @examples
##' variants <- searchVariants()
##' summary(variants)
##' summary(variants[[1]])
##' @export
searchVariants <- function(base_url,
                           variantSetId = "MWtnLXAzLXN1YnNldDptdm5jYWxs",
                           referenceName = "1",
                           start = NULL,
                           end = NULL,
                           fields = NULL,
                           converter = c) {
    pageToken <- NULL
    variants <- converter()
    repeat {
        result <- searchVariantsPage(base_url = base_url,
                                     variantSetId = variantSetId,
                                     referenceName = referenceName,
                                     start = start,
                                     end = end,
                                     fields = fields,
                                     pageToken = pageToken)
        pageToken <- result$nextPageToken
        # TODO improve performance,
        # see https://github.com/Bioconductor/GoogleGenomics/issues/17
        variants <- c(variants, converter(result$variants))
        if (is.null(pageToken)) {
            break
        }
        message(paste("Continuing variant query with the nextPageToken:",
                      pageToken))
    }

    message("Variants are now available.")
    variants
}

##' Convert variants to VRanges.
##'
##' Note that the Global Alliance for Genomics and Health API uses a 0-based
##' coordinate system.  For more detail, please see GA4GH discussions such
##' as the following:
##' \itemize{
##'    \item\url{https://github.com/ga4gh/schemas/issues/168}
##'    \item\url{https://github.com/ga4gh/schemas/issues/121}
##' }
##'
##' @param variants A list of R objects corresponding to the JSON objects
##'  returned by the GA4GH Variants API.
##' @param oneBasedCoord Convert genomic positions to 1-based coordinates.
##' @param slStyle The style for seqnames (chrN or N or...).  Default is UCSC.
##' @return \link[VariantAnnotation]{VRanges}
##' @family variants converter functions
##' @examples
##' variants1 <- searchVariants(converter = variantsToVRanges)
##' summary(variants1)
##' variants2 <- variantsToVRanges(searchVariants())
##' print(identical(variants1, variants2))
##' @export
variantsToVRanges <- function(variants, oneBasedCoord = TRUE,
                              slStyle = "UCSC") {
    if (missing(variants)) {
        return(VRanges())
    }

    variantsToVRangesHelper <- function(v) {
        if (TRUE == oneBasedCoord) {
            position <- as.integer(v$start) + 1
        } else {
            position <- as.integer(v$start)
        }

        ranges <- VRanges(
            seqnames = Rle(as.character(v$referenceName), 1),
            ranges = IRanges(start = position,
                             end = as.integer(v$end)),
            ref = as.character(v$referenceBases),
            alt = as.character(v$alternateBases[1]),  # TODO flatten per alt
            QUAL = as.numeric(v$quality),
            FILTER = as.character(v$filter))

        names(ranges) <- as.character(v$names[1])
        ranges
    }
    vranges <- do.call("c", lapply(variants, variantsToVRangesHelper))

    seqlevelsStyle(vranges) <- slStyle
    vranges
}

##' Convert variants to GRanges.
##'
##' Note that the Global Alliance for Genomics and Health API uses a 0-based
##' coordinate system.  For more detail, please see GA4GH discussions such
##' as the following:
##' \itemize{
##'    \item\url{https://github.com/ga4gh/schemas/issues/168}
##'    \item\url{https://github.com/ga4gh/schemas/issues/121}
##' }
##'
##' @param variants A list of R objects corresponding to the JSON objects
##'   returned by the GA4GH Variants API.
##' @param oneBasedCoord Convert genomic positions to 1-based coordinates.
##' @param slStyle The style for seqnames (chrN or N or...).  Default is UCSC.
##' @return \link[GenomicRanges]{GRanges}
##' @family variants converter functions
##' @examples
##' variants1 <- searchVariants(converter = variantsToGRanges)
##' summary(variants1)
##' variants2 <- variantsToGRanges(searchVariants())
##' print(identical(variants1, variants2))
##' @export
variantsToGRanges <- function(variants, oneBasedCoord = TRUE,
                              slStyle = "UCSC") {
    if (missing(variants)) {
        return(GRanges())
    }

    variantsToGRangesHelper <- function(v) {
        if (TRUE == oneBasedCoord) {
            position <- as.integer(v$start) + 1
        } else {
            position <- as.integer(v$start)
        }

        ranges <- GRanges(
            seqnames = Rle(as.character(v$referenceName), 1),
            ranges = IRanges(start = position,
                             end = as.integer(v$end)),
            REF = DNAStringSet(v$referenceBases),
            ALT = DNAStringSetList(v$alternateBases),
            QUAL = as.numeric(v$quality),
            FILTER = as.character(v$filter))

        names(ranges) <- as.character(v$names[1])
        ranges
    }
    granges <- do.call("c", lapply(variants, variantsToGRangesHelper))

    seqlevelsStyle(granges) <- slStyle
    granges
}
