\name{Ontime}
\alias{Ontime}
\docType{data}
\title{Ontime}
\description{Ontime data set}
\usage{data("Ontime")}
\format{
  The format is:
 int [1:500, 1:6] 11618 12889 14908 14635 13930 13930 12266 11292 14893 14524 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:500] "1" "2" "3" "4" ...
  ..$ : chr [1:6] "ORIGIN_AIRPORT_ID" "ORIGIN_AIRPORT_SEQ_ID" "ORIGIN_CITY_MARKET_ID" "DEST_AIRPORT_ID" ...
}
\details{It contains 621461 samples with 6 variables}
\source{US department of transportation}
\examples{
data(Ontime)
## maybe str(Ontime) ; plot(Ontime) ...
}
\keyword{datasets}
