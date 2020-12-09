#' apply simulated annealing to 2opt algorithm. 
#'
#' @param tsp the tsp instance
#' @param dima the distance matrix
#' @param T0 initial temperature
#' @param alpha the temperature reduction factor per iteration
#'
#' @return the solution
#' @export
simulatedAnnealing <- function(tsp, dima, T0 = 1e2, alpha = 0.9) {
  # code see slide 53 from week 10
  # cur_sol <- getRandomTour(tsp = tsp)
  cur_sol <- constructNearestNeighborSolution(tsp, dima)
  cur_dist <- calculateTourDistance(cur_sol, dima)
  best_sol <- cur_sol
  best_dist <- cur_dist
  T <- T0
  cnt_non_improving_subsequent_iterations <- 0
  while(TRUE) {
    # choose move randomly: pos 1, pos 2 for swapping with 2opt
    tmp <- sort(sample(1:nrow(tsp),2, replace = FALSE))
    pos1 <- tmp[1]
    pos2 <- tmp[2]
    # apply move
    new_sol <- apply2opt(cur_sol, pos1, pos2)
    new_dist <- calculateTourDistance(new_sol, dima)
    dist <- new_dist - cur_dist
    # if new sol is better, take it!
    if (dist < 0) {
      # improving moves are always allowed
      cur_sol <- new_sol
      cur_dist <- new_dist
    } else {
      # Generate random number in [0,1]
      u <- runif(1)
      if (exp(-dist/T) > u) {
        # if Temp allows accepting bad value
        cur_sol <- new_sol
        cur_dist <- new_dist
      }
    }
    if (cur_dist < best_dist) {
      best_sol <- cur_sol
      best_dist <- cur_dist
      # reset unsuccessful iteration counter
      cnt_non_improving_subsequent_iterations <- 0
      cat("best sol: ", best_dist, "\n")
    } else {
      # increase unsuccessful iteration counter
      cnt_non_improving_subsequent_iterations <- 
        cnt_non_improving_subsequent_iterations +1
    }
    T <- T / (1 + alpha*T)
    if (cnt_non_improving_subsequent_iterations > 3000)
      break
  }
  return(best_sol)
}

#' generate a random tour
#'
#' @param tsp the tsp instance
#'
#' @return a random tour
#' @export
getRandomTour <- function(tsp) {
  tour <- sample(tsp$stop_id, size = nrow(tsp), replace = FALSE)
  return(tour)
}

#' construction heuristic to create an initial tsp-solution based on nearest
#' neighbor starting from a given location
#'
#' @param tsp the tsp instance
#' @param dima the distance matrix
#' @param start_point_id starting point, where the tsp tour should start
#'
#' @return a valid tsp-tour starting at location start_point_id (not listing the
#'   last position equal to first one, this is considered in the tour-distance
#'   calculation function
constructNearestNeighborSolution <- function(tsp, dima, start_point_id = 1) {
  t <- c()
  # startpunkt für nearest neighbor k-heuristik wählen (WICHTIG: drop = F 
  # verhindert hier, dass spalte zu vektor konvertiert wird!)
  # willkürlich start bei position = 1
  tsp_ids <- tsp$stop_id
  cur_id <- tsp$stop_id[start_point_id]
  repeat {
    # insert stop cur_id in tour t
    t <- c(t, cur_id)
    # remove cur_id from available unvisited stops in tsp_ids
    tsp_ids <- tsp_ids[tsp_ids != cur_id, drop = FALSE]
    best_stop <- -1
    best_dist <- Inf
    if (length(tsp_ids) != 0) {
      # check all left stops to find the nearest
      for (potential_next_id in tsp_ids) {
        cur_dist <- getDimaDist(fromLoc = cur_id, toLoc = potential_next_id, 
                                dima = dima)
        if (cur_dist < best_dist) {
          best_dist <- cur_dist
          best_stop <- potential_next_id
        }
      }
      cur_id <- best_stop
    } else {
      break
    }
  }
  return(t)  
}

#' berechnet eine einzelne 2opt nachbarschaft für eine gegebene TSP-Lösung indem
#' die alle Knoten zwischen Node-ID firstNodeId und secondNodeId umgekehrt werden.
#' (Siehe dazu Slide 95 aus der Vorlesung)
#'
#' @param tsp_node_sequence Vektor von TSP Knoten-IDs
#' @param firstNodeId NodeId des ersten Knotens 
#' @param secondNodeId NodeId des zweiten Knotens 
#'
#' @return eine einzelne neue TSP-Lösung, die der 2opt Nachbarschaft mit den 
#' gegebenen input-Parametern entspricht
apply2opt <- function(tsp_node_sequence, firstNodeId, secondNodeId) {
  # validiere die werte
  num_nodes <- length(tsp_node_sequence)
  ret_val <- c()
  # mühsame Fallunterscheidung in R (Subsetting nicht konsistent: a[,1:0]) ist 
  # leider nicht leer, sondern liefert immer 1 Spalte. Daher Fallunterscheidung nötig
  # a) alle Nodes vor dem startknoten übernehmen
  if (firstNodeId > 1) {
    ret_val <- c(ret_val, tsp_node_sequence[1:(firstNodeId-1)])
  }
  # b) alle knoten zwischen firstNodeId und secondNodeId umkehren
  ret_val <- c(ret_val, tsp_node_sequence[rev(firstNodeId:secondNodeId)]) 
  # c) alle knoten nach letztem knoten übernehmen
  if (secondNodeId < num_nodes) {
    ret_val <- c(ret_val, tsp_node_sequence[(secondNodeId + 1):num_nodes])
  }
  return(ret_val)
}

#' render a simulated annealing tour by first initializing the dima
#'
#' @param tsp 
#'
#' @return the tour
#' @export
renderTour <- function(tsp) {
  # build new or load existing distance-matrix
  dima <- calculateDima(tsp)
  tour <- simulatedAnnealing(tsp = tsp, dima = dima)
  return(tour)
}

#' calculate  the distance-matrix
#'
#' @param tsp the tsp instance
#'
#' @return the dima
#' @export
calculateDima <- function(tsp) {
  dima.dt <- data.table::data.table(loc_from = character(0),
                         loc_to = character(0),
                         dist = numeric(0), stringsAsFactors = FALSE)
  n <- nrow(tsp)
  # since we don't have a matrix but rather a lookup-table, we need to keep
  # track of the row-id of the 'dima'
  for (from in 1:n)  {
    fromId <- tsp[from,]$stop_id
    for (to in from:n) {
      toId <- tsp[to,]$stop_id
      lat1 <- tsp[from,]$lat
      lng1 <- tsp[from,]$lng
      lat2 <- tsp[to,]$lat
      lng2 <- tsp[to,]$lng
      result <- sqrt((lat2-lat1)^2 + (lng2-lng1)^2)
      dima.dt <- rbind(dima.dt, data.table::data.table(loc_from = fromId,
             loc_to = toId,  dist = result, stringsAsFactors = FALSE))
    }
    # now set keys on dima
    data.table::setkey(dima.dt, loc_from, loc_to)
  }
  return(dima.dt)
}

#' calculate total tour distance
#'
#' @param tour the tour as a sequence of stopIds
#' @param dima the distance matrix
#'
#' @return the total tsp-distance
#' @export
calculateTourDistance <- function(tour, dima) {
  dist <- 0.0
  for (i in 2:length(tour)) {
    a <- tour[i-1]
    b <- tour[i]
    dist <- dist + getDimaDist(fromLoc = a, toLoc = b, dima = dima)
  }
  # now return to start
  retDist <- getDimaDist(fromLoc = tour[1], toLoc = tour[length(tour)], 
                         dima = dima)
  return(dist + retDist)
}

#' get distance between two points from dima 
#'
#' @param fromLoc from location ID
#' @param toLoc to location ID
#' @param dima the distance matrix
#'
#' @return the distance
#' @export
getDimaDist <- function(fromLoc, toLoc, dima) {
  dimaEntry <- dima[loc_from == fromLoc & loc_to == toLoc]
  if (nrow(dimaEntry) != 1) {
    # if a / b lookup failed, try other way round (since we store only one
    # direction) in the distance matrix.
    dimaEntry <- dima[I(loc_from == toLoc & loc_to == fromLoc)]
    if (nrow(dimaEntry) != 1) {
      stop(
        paste0(
          "Expected to find exactly one dima entry corresponding to the given loc_from/loc_to-pair ",
          loc_from,
          loc_to,
          " but found 0 or more than 1."
        )
      )
    }
  }
  return(dimaEntry$dist)
}
# best sol:  7538.932 
#  [1] 16  2 17 30 21  0 48 31 44 18 40  7  8  9 32 50 10 51 13 12 46 25 26 27 11 24  3  5 14 42 39 38 35 34 33
# [36] 36 37  4 23 47 45 43 15 28 49 19 22 29 20 41  6  1

tsp <- read.table("berlin52.tsp", skip = 2, col.names = c("stop_id", "lat", "lng"))
tour <- renderTour(tsp)
print(tour)
