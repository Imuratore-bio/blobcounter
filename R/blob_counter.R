#' blobcounter
#'
#' @param mask_file the mask file used with trackR
#' @param tracks_file the trajectory csv output file from trackR
#'
#' @return a data frame
#' @export
#'
#' @examples
#' mask_file_input <- "C:/Users/imura/Downloads/7.19(B1)60D_mask2_1.png"
#' tracks_file_input <- 'C:/Users/imura/Downloads/sanity_tracks4.csv'
#' blob_counter(mask_file_input, tracks_file_input)

#beginning of function
blob_counter <- function(mask_file, tracks_file) {


  #out of all blobs created, how many unique ones have a position close to the centroid of a given box?

  #load packages for image processing and automatic corner detection
  requireNamespace("image.ContourDetector")
  requireNamespace("magick")
  requireNamespace("pixmap")
  requireNamespace("grDevices")
  requireNamespace("image.CornerDetectionHarris")
  requireNamespace("graphics")

  #read in tracks file
  tracks <- utils::read.csv(file = tracks_file)

  #read in mask file
  image_x <- magick::image_read(mask_file)

  #convert to magick format to use with corner detector
  x <- magick::image_convert(image_x, format = "pgm", depth = 8)

  #detect corners
  corners <- image.CornerDetectionHarris::image_harris(x)

  #for some reason the Harris corner detector results are vertically mirrored; this fixes that
  corners$y <- 1080 - corners$y

  #warning to user
  if (length(corners$x) != 16) {
    print("Warning, desired number of corners is 16. Current number of corners is:")
    print(length(corners$x))

  }

  #view identified corners
  plot(x)
  graphics::points(corners$x, corners$y, col = "red", pch = 20, lwd = 0.5)

  #define a list to store corner coordinates
  line_coords <- c()

  #define a counter variable
  l <- 0

  #loop through every pair of two corner points
  for (len in 1:length(corners$x))
  {
    for (len2 in 1:length(corners$x))

    {
      #add the two points and the length between them to a data frame
      line_coords$x1 <- append(line_coords$x1, corners$x[len])
      line_coords$y1 <- append(line_coords$y1, corners$y[len])
      line_coords$x2 <- append(line_coords$x2, corners$x[len2])
      line_coords$y2 <- append(line_coords$y2, corners$y[len2])
      #getting the length between the points
      line_coords$l <- append(line_coords$l, sqrt((corners$x[len2] - corners$x[len])^2 + (corners$y[len2] - corners$y[len])^2))
    }
  }
  #making it into a data frame
  lines <- data.frame(line_coords)

  #take the point pairs with the 16 shortest line distances between them (excluding 0 distances between identical points)
  #aka only orthogonal lines within a square
  #this can fail if some squares are bigger, with small square hypotenuses being smaller
  non_zero <- lines[lines$l>0,]
  ordered <- non_zero[order(non_zero$l),]

  mn <- pmin(ordered$x1, ordered$y1, ordered$x2, ordered$y2)
  mx <- pmax(ordered$x1, ordered$y1, ordered$x2, ordered$y2)
  int <- as.numeric(interaction(mn, mx))
  no_dupes <- ordered[match(unique(int), int),]

  no_dupes$mx <- (no_dupes$x1 + no_dupes$x2)/2
  no_dupes$my <- (no_dupes$y1 + no_dupes$y2)/2

  #remove highly similar midpoint values since there will be non-identical duplicates
  #for every combo of midpoints
  for (m in (1:length(no_dupes$mx))) {
    for (m2 in (1:length(no_dupes$mx))) {
      #see if they're similar
      if (abs(no_dupes$mx[m2]-no_dupes$mx[m]) < 5 & abs(no_dupes$my[m2]-no_dupes$my[m]) < 5) {
        #replace one with the other
        no_dupes$mx[m2] <- no_dupes$mx[m]
        no_dupes$my[m2] <- no_dupes$my[m]
      }
    }
  }

  #remove any lines where the middle area is black pixels (i.e. between squares) or exclusively white pixels (i.e. crossing the middle of a square)

  #read in file again in a different format to extract pixel hex values
  img <- png::readPNG(mask_file)
  #check to see whether it has an alpha channel
  if (length(dim(img)) > 3) {
    y <- rgb(img[,,1], img[,,2], img[,,3], alpha = img[,,4])
  } else {
    y <- rgb(img[,,1], img[,,2], img[,,3])
  }
  yg <- colorspace::desaturate(y)
  yn <- grDevices::col2rgb(yg)[1, ]/255
  dim(y) <- dim(yg) <- dim(yn) <- dim(img)[1:2]

  #needs to be transposed
  y <- t(y)

  #test whether all midpoints have at least one pure white pixel nearby
  #this is to exclude any short lines connecting points outside of squares

  #make empty variable to hold list of True/False values
  delete <- c()

  #for each set of line midpoints
  for (line in 1:length(no_dupes$mx)){
    #test whether a midpoint itself or at least one pixel in any cardinal direction going out three pixels is pure white
    delete <- c(delete,
                (grepl( "#FFFFFF", y[(round(no_dupes$mx[line])), (1080-round(no_dupes$my[line]))], fixed = TRUE) | grepl( "#FFFFFF", y[(round(no_dupes$mx[line])+3), (1080-round(no_dupes$my[line])+3)], fixed = TRUE) | grepl( "#FFFFFF", y[(round(no_dupes$mx[line])+3), (1080-round(no_dupes$my[line]))], fixed = TRUE) | grepl( "#FFFFFF", y[(round(no_dupes$mx[line])), (1080-round(no_dupes$my[line])+3)], fixed = TRUE) | grepl( "#FFFFFF", y[(round(no_dupes$mx[line])+3), (1080-round(no_dupes$my[line])-3)], fixed = TRUE) | grepl( "#FFFFFF", y[(round(no_dupes$mx[line])-3), (1080-round(no_dupes$my[line])+3)], fixed = TRUE) | grepl( "#FFFFFF", y[(round(no_dupes$mx[line])-3), (1080-round(no_dupes$my[line]))], fixed = TRUE) | grepl( "#FFFFFF", y[(round(no_dupes$mx[line])), (1080-round(no_dupes$my[line])-3)], fixed = TRUE) | grepl( "#FFFFFF", y[(round(no_dupes$mx[line])+3), (1080-round(no_dupes$my[line]))], fixed = TRUE)| grepl( "#FFFFFF", y[(round(no_dupes$mx[line])), (1080-round(no_dupes$my[line])+3)], fixed = TRUE)  | grepl( "#FFFFFF", y[(round(no_dupes$mx[line])-3), (1080-round(no_dupes$my[line])-3)], fixed = TRUE)))

  }

  #if there was at least one white pixel, keep that line, but exclude others
  no_dupes <- no_dupes[delete,]

  #make a new empty variable to hold list of True/False values
  delete <- c()

  #for each set of line midpoints
  for (line in 1:length(no_dupes$mx)){
    #test whether a midpoint itself and one pixel in every cardinal direction going out five pixels is pure white
    delete <- c(delete,
                (grepl( "#FFFFFF", y[(round(no_dupes$mx[line])), (1080-round(no_dupes$my[line]))], fixed = TRUE) & grepl( "#FFFFFF", y[(round(no_dupes$mx[line])+5), (1080-round(no_dupes$my[line])+5)], fixed = TRUE) & grepl( "#FFFFFF", y[(round(no_dupes$mx[line])+5), (1080-round(no_dupes$my[line])-5)], fixed = TRUE) & grepl( "#FFFFFF", y[(round(no_dupes$mx[line])+5), (1080-round(no_dupes$my[line]))], fixed = TRUE)& grepl( "#FFFFFF", y[(round(no_dupes$mx[line])), (1080-round(no_dupes$my[line])+5)], fixed = TRUE) & grepl( "#FFFFFF", y[(round(no_dupes$mx[line])-5), (1080-round(no_dupes$my[line]))], fixed = TRUE)& grepl( "#FFFFFF", y[(round(no_dupes$mx[line])), (1080-round(no_dupes$my[line])-5)], fixed = TRUE)  & grepl( "#FFFFFF", y[(round(no_dupes$mx[line])-5), (1080-round(no_dupes$my[line])+5)], fixed = TRUE)   & grepl( "#FFFFFF", y[(round(no_dupes$mx[line])-5), (1080-round(no_dupes$my[line])-5)], fixed = TRUE)))

  }

  #if all pixels were white, exclude that line, but keep others
  no_dupes <- no_dupes[!delete,]

  #take only the 16 shortest lines to get outlines of squares

  #reset row names since some rows were deleted before
  row.names(no_dupes) <- NULL

  #or fewer than 16 if there were fewer than 16 corners detected
  if (length(corners$x) < 16){
  dupe_limit <- 16 - ((16 - length(corners$x))*2)

  } else {
  dupe_limit <- 16
  }

  short_lines <- no_dupes[1:dupe_limit,]

  #view line segments
  plot(x)
  graphics::segments(short_lines$x1, short_lines$y1, short_lines$x2, short_lines$y2, col = "red", lwd = 0.5)

  #for selected point combos that share one point, or have an extremely similar point
  for (line1 in (1:length(short_lines$x1))) {
    for (line2 in 1:length(short_lines$x1)) {
      if (short_lines[line1,]$x1 == short_lines[line2,]$x1 & short_lines[line1,]$y1 == short_lines[line2,]$y1 & short_lines[line1,]$x2 != short_lines[line2,]$x2 & short_lines[line1,]$y2 != short_lines[line2,]$y2) {
        #take the two non-shared points and find their midpoint
        #store it as cx and cy
        #each case is defined individually
        short_lines$cx[line1] <- (short_lines$x2[line1] + short_lines$x2[line2])/2
        short_lines$cy[line1] <- (short_lines$y2[line1] + short_lines$y2[line2])/2
      }

      if (short_lines[line1,]$x1 == short_lines[line2,]$x2 & short_lines[line1,]$y1 == short_lines[line2,]$y2 & short_lines[line1,]$x2 != short_lines[line2,]$x1 & short_lines[line1,]$y2 != short_lines[line2,]$y1) {
        #store it as cx and cy
        short_lines$cx[line1] <- (short_lines$x2[line1] + short_lines$x1[line2])/2
        short_lines$cy[line1] <- (short_lines$y2[line1] + short_lines$y1[line2])/2
      }

      if (short_lines[line1,]$x2 == short_lines[line2,]$x1 & short_lines[line1,]$y2 == short_lines[line2,]$y1 & short_lines[line1,]$x1 != short_lines[line2,]$x2 & short_lines[line1,]$y1 != short_lines[line2,]$y2) {
        #store it as cx and cy
        short_lines$cx[line1] <- (short_lines$x1[line1] + short_lines$x2[line2])/2
        short_lines$cy[line1] <- (short_lines$y1[line1] + short_lines$y2[line2])/2
      }

      if (short_lines[line1,]$x2 == short_lines[line2,]$x2 & short_lines[line1,]$y2 == short_lines[line2,]$y2 & short_lines[line1,]$x1 != short_lines[line2,]$x1 & short_lines[line1,]$y1 != short_lines[line2,]$y1) {
        #store it as cx and cy
        short_lines$cx[line1] <- (short_lines$x1[line1] + short_lines$x1[line2])/2
        short_lines$cy[line1] <- (short_lines$y1[line1] + short_lines$y1[line2])/2
      }
    }
  }

  #automatically defining the centroid proximity threshold (how close does a blob need to be to a centroid to count)
  #based on the distance of half the diagonal, plus some wiggle room to account for differently sized squares
  centroid_proximity_thres <- sqrt((short_lines$x1[1] - short_lines$cx[1])^2 + (short_lines$y1[1] - short_lines$cy[1])^2 ) + 20

  #remove highly similar centroid values since there will be non-identical duplicates
  #for every combo of centroid points
  for (c in (1:length(short_lines$cx))) {
    for (c2 in (1:length(short_lines$cx))) {
      #see if they're similar
      #used to be or
      if (abs(short_lines$cx[c2]-short_lines$cx[c]) < 5 & abs(short_lines$cy[c2]-short_lines$cy[c]) < 5) {
        #replace one with the other
        short_lines$cx[c2] <- short_lines$cx[c]
        short_lines$cy[c2] <- short_lines$cy[c]
      }
    }
  }


  #exclude black centroid points

  delete <- c()

  for (line in 1:length(short_lines$cx)){
      delete <- c(delete,
                grepl("#000000", y[(round(short_lines$cx[line])), (1080-round(short_lines$cy[line]))])
    )

  }

  short_lines <- short_lines[!delete,]

  #user warnings about centroid similarity threshold
  if (length( short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx  ) != 4) {
    print("Warning, desired number of centroids is 4. Current number of centroids is:")
    print(length( short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx  ))
  }

  #view the final centroids
  plot(x)
  graphics::points(short_lines$cx, short_lines$cy, col = "red", pch = 20, lwd = 0.5)

  #variables to hold column data for output csv
  #centroid position and count of blobs that pass near it
  centroid_x_value <- c()
  centroid_y_value <- c()
  count <- c()

  #for each centroid point
  for (c in 1:length(short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx)) {
    #reset the list of used blobs
    used_blobs <- c()
    #reset count of blobs that have passed nearby
    blob_count <- 0
    #for each row in the tracks file
    for (p in 1:length(tracks$x)) {
      #if it is close to the centroid in question and not already counted
      #and closer to this centroid than to any other centroid
      if (!(tracks$track[p] %in% used_blobs) & (abs(tracks$x[p] - short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx[c]) < centroid_proximity_thres) & (abs(tracks$y[p] - short_lines[!duplicated(short_lines[c("cx","cy")]),]$cy[c]) < centroid_proximity_thres)) {

        #set this to False unless a point is found to be closer to a given centroid than to all others
        comparison_check <- FALSE

        #compare to all non-focal centroids
        for (others in 1:length(short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx)) {

          #checking distance between each point and a given centroid
          #is it closer to the focal centroid than to any other centroid?
          if( sqrt((tracks$y[p] - short_lines[!duplicated(short_lines[c("cx","cy")]),]$cy[c])^2 + (tracks$x[p] - short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx[c])^2) < sqrt((tracks$y[p] - short_lines[!duplicated(short_lines[c("cx","cy")]),]$cy[others])^2 + (tracks$x[p] - short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx[others])^2 ) ) {

            comparison_check <- TRUE

          }

        }

        #add the name to the already counted list
        if (comparison_check == TRUE) {
          used_blobs <- append(used_blobs, tracks$track[p])
          #up the blob count
          blob_count <- blob_count + 1
        }
      }

    }
    #finishing the above loops takes the longest
    #based on result of above loop append row to variables holding data for data frame
    centroid_x_value <- append(centroid_x_value, unique(short_lines$cx)[c])
    centroid_y_value <- append(centroid_y_value, unique(short_lines$cy)[c])
    count <- append(count, blob_count)
    #might be redundant
    blob_count <- 0
  }

  #create and return output data frame
  output <- data.frame(centroid_x_value, centroid_y_value, count)
  #sort by x value so it's basically displaying the results for each square from left to right
  output_sorted <- output[order(centroid_x_value),]
  return(output_sorted)

  #end of function
}
