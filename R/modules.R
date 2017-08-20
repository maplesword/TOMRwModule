#' Generate the weighted igraph object
#'
#' This function generate the igraph object based on the given
#' interactions. The weight of each edge is determined as the
#' heat differences between the two interacting vertices.
#'
#' @param links A two-column matrix or data frame, with each row representing one interaction
#' @param directed Whether the interactions are directed
#' @param heat A numeric vector representing vertices' heat
#' @return A igraph object representing the interaction network. Edges are weighted accordingly
#' @export
weightedGraph <- function(links = interactions_brain, directed = FALSE, heat = NULL){
	graph <- igraph::graph_from_data_frame(links, directed=directed)
	if(is.numeric(heat) & sum(names(igraph::V(graph)) %in% names(heat)) == length(igraph::V(graph))){
		heat.diff <- apply(links, 1, function(x){
			heats <- heat[x]
			(max(heats)-min(heats))/max(abs(heats))
		})
		igraph::E(graph)$weight <- 1-(heat.diff/2)
	} else{
		igraph::E(graph)$weight <- rep(1, length(igraph::E(graph)))
	}
	
	return(graph)
}

#' Module detection based on Topological Overlap Matrix
#'
#' This function calculates TOM matrix using the implementation
#' in WGCNA package, and applied hierarchical clustering to identify
#' modules.
#'
#' @param graph The igraph object representing the interaction network
#' @param minSize The minimum requirement of identified modules
#' @return A list with each element representing vertices of one module
#' @export
TOMModules <- function(graph, minSize = 20, ...){
	adj.mat <- as.matrix(igraph::get.adjacency(graph, type="both", attr="weight"))/2
	TOM <- WGCNA::TOMsimilarity(adj.mat, TOMType="unsigned")
	dissTOM <- 1-TOM
	geneTree <- hclust(as.dist(dissTOM), method = "average")
	dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minSize, ...)
	dynamicColors <- WGCNA::labels2colors(dynamicMods)
	names(dynamicColors) <- rownames(igraph::get.adjacency(graph, type="both", attr="weight"))
	modules.TOM <- lapply(unique(dynamicColors), function(x){
		names(dynamicColors)[dynamicColors==x]
	})
	
	modules.TOM <- modules.TOM[unique(dynamicColors)!="grey"]
	names(modules.TOM) <- paste0("M", 1:length(modules.TOM))
	return(modules.TOM)
}

#' The automatic procedure to detect modules in the interaction network
#'
#' This function integrates all steps including missing value
#' imputation, heat diffusion and module detection into an
#' automatic pipeline.
#'
#' @param links A two-column matrix or data frame, with each row representing one interaction
#' @param heat A numeric vector representing initial heat of each vectex in the network
#' @param force_HD If there are no more than half of the vertices have known heat, whether the heat diffusion procedure should be applied
#' @param fill A numeric value, character, or function to specialize the way to fill in missing values
#' @param directed Whether the interactions are directed
#' @param beta Insulating parameter between 0 and 1
#' @param minSize The minimum requirement of identified modules
#' @param ... Other parameters for 'dynamicTreeCut' in TOM-based module detection
#' @return A list with each element representing vertices of one module
#' @export
heatModules <- function(links = interactions_brain, heat = NULL, force_HD = FALSE, fill = "mean", directed = FALSE, beta = 0.55, minSize = 20, ...){
	vertices <- union(links[,1], links[,2])
	if(is.null(heat)){
		graph <- weightedGraph(links = links, directed = directed, heat = NULL)
	} else{
		heat <- heat[names(heat) %in% vertices]
		if(! force_HD & length(heat) < length(vertices)/2){
			warning("Heat is assigned to less than half of vertices. Use unweighted graph for module identification.")
			graph <- weightedGraph(links = links, directed = directed, heat = NULL)
		} else{
			heat <- fillHeat(links, heat, fill = fill)
			cat("Done missing value imputation\n")
			diffusion.mat <- diffusionMat(links, directed = directed, beta = beta)
			cat("Done calculating diffusion matrix\n")
			heat.eq <- equilibriumHeat(diffusion.mat, heat)
			cat("Done calculating equilibrium heat\n")
			graph <- weightedGraph(links = links, directed = directed, heat = heat.eq)
			cat("Generated weighted graph\n")
		}
	}
	cat("Start module detection...\n")
	modules <- TOMModules(graph, minSize = minSize, ...)
	res <- list(diffusMat = diffusion.mat, weightedGraph = graph, modules = modules)
	return(res)
}
