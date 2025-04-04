% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{map_prep}
\alias{map_prep}
\title{Prepare map data for visualization}
\usage{
map_prep(
  ukcoast_poly,
  ukcoast_line,
  uk_poly_converted,
  uk_line_converted,
  uk_grid_converted
)
}
\arguments{
\item{ukcoast_poly}{Character string. Path to UK polygon shapefile.}

\item{ukcoast_line}{Character string. Path to UK coastline shapefile.}

\item{uk_poly_converted}{Character string. Path where the processed UK polygon will be saved.}

\item{uk_line_converted}{Character string. Path where the processed UK coastline will be saved.}

\item{uk_grid_converted}{Character string. Path where the interpolation grid will be saved.}
}
\value{
None. The function saves the processed spatial objects to the specified output files.
}
\description{
Prepares UK map data for use in OTU distribution visualization.
}
\details{
The function:
\enumerate{
\item Reads UK coastline shape files
\item Sets the correct coordinate reference system
\item Creates an interpolation grid based on the UK polygon bounding box
\item Saves the processed spatial objects as shapefiles
}
}
\examples{
map_prep(
  "data/map/ukcoast1.shp",
  "data/map/ukcoast_line.shp",
  "output/map/uk_poly.shp",
  "output/map/uk_line.shp",
  "output/map/uk_grid.shp"
)

}
