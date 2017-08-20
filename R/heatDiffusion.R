#' Generate heat based on log2fc and p-value
#' 
#' This function generate the initial heat for genes with the
#' given fold change and p-value vector, as the square root of
#' Log2-transformed FC times -Log10-transformed P-value.
#'
#' @param fc A numeric vector indicating fold change of expression
#' @param p A numeric vector indicating P-value of respective test
#' @return A numeric vector of intial heat with the same length as fc and p-value
#' @export
initialHeat <- function(fc = NULL, p = NULL){
	heat <- 1
	if(is.null(fc) & is.null(p)){
		warning("Both Log2FC and P are NULL. Return 1 as default.")
	} else if(is.null(fc)){
		heat <- -log10(p)
	} else if(is.null(p)){
		heat <- log2(fc)
	} else{
		heat <- -log10(p) * log2(fc)
	}
	heat <- sign(heat) * sqrt(abs(heat))
	return(heat)
}

#' Fill in the missing heat
#'
#' This function finds the vertices in the interaction network
#' whose initial heat is not given, and fill with heat using
#' specific value or function.
#' 
#' @param links A two-column matrix or data frame, with each row representing one interaction
#' @param heat A numeric vector representing initial heat of each vectex in the network
#' @param fill A numeric value, character, or function to specialize the way to fill in missing values
#' @return A numeric vector representing heat of vertices
#' @export
fillHeat <- function(links = interactions_brain, heat = vertices_heat, fill = "mean", ...){
	if(! is.vector(heat)){
		warning("'heat' is not a numeric vector. Force converted.")
		heat <- as.matrix(heat)[,1]
	}
	nodes <- sort(unique(c(links[,1], links[,2])))
	if (is.numeric(fill)){
		scores <- rep(fill, length(nodes))
	} else if (is.function(fill)){
		scores <- fill(length(nodes), ...)
	} else if (is.character(fill)){
		if (fill == "mean"){
			scores <- rep(mean(heat), length(nodes))
		} else if (fill == "median"){
			scores <- rep(median(heat), length(nodes))
		}
	} else if (is.na(fill)){
		scores <- rep(fill, length(nodes))
	}
	names(scores) <- nodes
	scores[intersect(names(scores), names(heat))] <- heat[intersect(names(scores), names(heat))]
	return(scores)
}

# get the adjacent matrix representing the graph given the edges
get.adjacent.matrix <- function(links, directed = FALSE){
	graph <- igraph::graph_from_data_frame(links, directed = directed, vertices = NULL)
	adj.mat <- igraph::as_adjacency_matrix(graph)
	return(adj.mat)
}

# get the heat influence matrix for each step of heat diffusion given the edges
get.influence.matrix <- function(links, directed = FALSE){
	mat <- get.adjacent.matrix(links, directed = directed)
	res <- apply(mat, 1, function(x){ x/sum(x) })
	return(res)
}

#' Calculate the heat diffusion matrix when reaching equilibrium
#'
#' This function calculate the equilibrium heat diffusion
#' matrix, using HotNet2 algorithm (Leiserson
#' MDM, et al. Nat Genet. 2015).
#'
#' @param links A two-column matrix or data frame, with each row representing one interaction
#' @param directed Whether the interactions are directed
#' @param beta Insulating parameter between 0 and 1
#' @return The heat diffusion matrix
#' @export
diffusionMat <- function(links = interactions_brain, directed = FALSE, beta = 0.55){
	w <- get.influence.matrix(links, directed = directed)
	i <- diag(nrow(w))
	F <- beta * solve(i - (1 - beta) * w)
	return(F)
}

#' Calculate heat after reaching equilibrium
#'
#' This function calculate the equilibrium heat of each
#' vertex in the network.
#'
#' @param diffusion The diffusion matrix
#' @param heat A numeric vector or one-column matrix representing vertices' initial heat
#' @return A numeric vector representing equilibrium heat of vectices
#' @export
equilibriumHeat <- function(diffusion, heat){
	if(is.vector(heat)){
		heat <- heat[rownames(diffusion)]
		heat <- as.matrix(heat)
	} else{
		warning("'heat' is not a numeric vector. Force converted.")
		heat <- as.matrix(heat)[rownames(diffusion),1]
	}
	E <- (diffusion %*% heat)[,1]
	names(E) <- rownames(diffusion)
	return(E)
}

