#' blob_direction
#'
#' @param tracks_file the trajectory csv output file from trackR
#'
#' @return a data frame
#' @export
#'
#' @examples
#' tracks_file_input <- 'C:/Users/imura/Downloads/sanity_tracks4.csv'
#' blob_direction(tracks_file_input)

#direction and disappearance only blob counter
blob_direction <- function(tracks_file) {
  #load any packages

  #read in tracks file
  tracks <- read.csv(file = tracks_file)

  #make empty variables for data frame
  ExitLeft <- c()
  ExitRight <- c()
  minutes <- c()

  #outer loop
  for (min in 1:(round(max(tracks$frame)/1800)-1) ) {

    left_blobs <- c()
    right_blobs <- c()

    lower <- (min*1800+1)
    upper <- ((min*1800)+1801)
    #need to check here, does the blob exist between these two numbers?
    tracks_include <- subset(tracks, frame>=lower & frame<upper)

    #and if so, does it also not exist for any high frame numbers?
    #get the blobs to exclude
    tracks_exclude <- subset(tracks, frame>=upper)

    #data frame only including desired blobs
    blobs <- subset(tracks,(frame %in% tracks_include$frame) & !(frame %in% tracks_exclude$frame))

    #ORDER BLOBS BY TRACK AND THEN BY INCREASING FRAME VALUE
    blobs_order <- blobs[with(blobs, order(track, frame)),]

    #what direction is each blob going?
    #NEED TO FIX THIS TO BE FOR EVERY MINUTE
    for (b in unique(blobs_order$track)) {

      if (blobs_order[which(blobs_order$track==b),]$x[1] < blobs_order[which(blobs_order$track==b),]$x[length(blobs_order[which(blobs_order$track==b),]$x)]){
        right_blobs <- append(right_blobs, b)
      }
      if (blobs_order[which(blobs_order$track==b),]$x[1] > blobs_order[which(blobs_order$track==b),]$x[length(blobs_order[which(blobs_order$track==b),]$x)]){
        left_blobs <- append(left_blobs, b)
      }

    }

    unique_left_blobs <- length(unique(left_blobs))
    unique_right_blobs <- length(unique(right_blobs))

    #append count of passing blobs to data frame variables
    #make these data frames and record "min" as well
    ExitLeft <- append(ExitLeft, unique_left_blobs)
    ExitRight <- append(ExitRight, unique_right_blobs)
    minutes <- append(minutes, min)

    #resetting variables
    right_blobs <- c()
    left_blobs <- c()

  }

  #build output data frame
  #consider appending this to existing csv instead but keep in mind need to match/adjust length
  output <- data.frame(minutes, ExitLeft, ExitRight)

  #return data frame
  return(output)

  #end of function
}
