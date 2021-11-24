#beginning of function
blob_counter <- function(mask_file, tracks_file, corner_thres = 120, centroid_similarity_thres = 10) {
  #out of all blobs created, how many unique ones have a position close to the centroid of a given box?

  #load packages for image processing and automatic corner detection
  library(image.ContourDetector)
  library(magick)
  library(pixmap)
  library(image.CornerDetectionF9)

  #read in tracks file
  tracks <- read.csv(file = tracks_file)

  #read in mask file
  image_x <- image_read(mask_file)

  #convert to pgm to use with corner detector
  x <- image_convert(image_x, format = "pgm", depth = 8)
  f <- tempfile(fileext = ".pgm")
  image_write(x, path = f, format = "pgm")
  m <- read.pnm(file = f, cellres = 1)

  #define centroid points for the four boxes in mask file

  #get corners
  #can adjust this threshold if we don't get the expected 16 corners for 4 squares
  corners <- image_detect_corners(m@grey * 255,
                                  threshold = corner_thres,
                                  suppress_non_max = TRUE)

  #warnings to user
  if (length(corners$x) > 16) {
    print("Warning, corner threshold too low, adjust up")
  }

  if (length(corners$x) < 16) {
    print("Warning, corner threshold too high, adjust down")
  }

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
  non_zero <- lines[lines$l>0,]
  ordered <- non_zero[order(non_zero$l),]

  mn <- pmin(ordered$x1, ordered$y1, ordered$x2, ordered$y2)
  mx <- pmax(ordered$x1, ordered$y1, ordered$x2, ordered$y2)
  int <- as.numeric(interaction(mn, mx))
  no_dupes <- ordered[match(unique(int), int),]

  short_lines <- no_dupes[1:16,]

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
      if (abs(short_lines$cx[c2]-short_lines$cx[c]) < centroid_similarity_thres | abs(short_lines$cy[c2]-short_lines$cy[c]) < centroid_similarity_thres) {
        #replace one with the other
        short_lines$cx[c2] <- short_lines$cx[c]
        short_lines$cy[c2] <- short_lines$cy[c]
      }
    }
  }

  #user warnings about selecting centroid similarity threshold
  #you should end up with 4 unique centroid points so if you don't you need to adjust
  if (length( short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx  ) > 4) {
    print("Warning, centroid similarity threshold too low, adjust up")
  }

  if (length( short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx  ) < 4) {
    print("Warning, centroid similarity threshold too high, adjust down")
  }

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
      if ((abs(tracks$x[p] - short_lines[!duplicated(short_lines[c("cx","cy")]),]$cx[c]) < centroid_proximity_thres) & (abs(tracks$y[p] - short_lines[!duplicated(short_lines[c("cx","cy")]),]$cy[c]) < centroid_proximity_thres) & !(tracks$track[p] %in% used_blobs)) {
        #add the name to the already counted list
        used_blobs <- append(used_blobs, tracks$track[p])
        #up the blob count
        blob_count <- blob_count + 1
      }
    }

    #based on result of above loop append row to variables holding data for data frame
    centroid_x_value <- append(centroid_x_value, unique(short_lines$cx)[c])
    centroid_y_value <- append(centroid_y_value, unique(short_lines$cy)[c])
    count <- append(count, blob_count)
  }

  #create and return output data frame
  output <- data.frame(centroid_x_value, centroid_y_value, count)
  return(output)

  #end of function
}
