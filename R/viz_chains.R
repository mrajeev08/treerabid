# Get chains (i.e. chain membership for each case)
get_chain_id <- function(from, to, date_dt) {
  
  I_gr <- data.table(from = from[!is.na(from)], to = to[!is.na(from)])
  gr <- graph_from_data_frame(d = I_gr, vertices = unique(to), directed = TRUE)
  comps <- components(gr)
  comp_dt <- data.table(membership = comps$membership, id_case = as.numeric(names(comps$membership)))
}

# Get graph
get_graph <- function(tree) {
  
  from <- tree$id_progen
  to <- tree$id_case
  I_gr <- data.table(from = from[!is.na(from)], to = to[!is.na(from)])
  gr <- graph_from_data_frame(d = I_gr, 
                              vertices = tree[, .(id_case,
                                                  t = t - ymd("2002-01-01"), x_coord, y_coord)], 
                              directed = TRUE)
  V(gr)$membership <- components(gr)$membership
  return(gr)
}

# Get chain size & lengths
get_chain_stats <- function(from, to) {
  
  I_gr <- data.table(from = from[!is.na(from)], to = to[!is.na(from)])
  gr <- graph_from_data_frame(d = I_gr, vertices = unique(to), directed = TRUE)
  data.table(membership = seq(components(gr)$no), sizes = components(gr)$csize, 
             lengths = unlist(lapply(decompose(gr), diameter)))
}

# Chain persistence across time
get_chain_persistence <- function(from, to, date_dt) {
  
  I_gr <- data.table(from = from[!is.na(from)], to = to[!is.na(from)])
  gr <- graph_from_data_frame(d = I_gr, vertices = unique(to), directed = TRUE)
  comps <- components(gr)
  comp_dt <- data.table(membership = comps$membership, id_case = names(comps$membership))
  comp_dt[date_dt, on = "id_case"][, .(start_date = min(t), end_date = max(t)), by = "membership"]
  
}


