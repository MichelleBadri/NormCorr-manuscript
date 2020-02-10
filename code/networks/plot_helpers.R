## Plot network highlighting modularity / cross edges
plot_mod_net <- function(G, main="", dropiso=TRUE) {
  if (dropiso) {
    ind <- igraph::V(G)$name[which(igraph::degree(G) < 1)]
    G <- igraph::delete.vertices(G, ind)
  }

  E(G)$weight <- ifelse(E(G)$cross, 0.1, 1)
  E(G)$color  <- ifelse(E(G)$cross, "grey", "#FF000033")
  V(G)$label  <- V(G)$cmemb
  l <- cbind(V(G)$l1, V(G)$l2)
  plot(G, vertex.shape="none", vertex.label.cex=0.8,
      vertex.label.color="black",
      main=main, layout=l, #vertex.size=4, edge.width=1,
      rescale=TRUE, ylim=c(-0.7,0.7), xlim=c(-0.7,0.7))
}


## Plot network highlighting node taxonomy / edge weights&sign
col <- unlist(yaml::read_yaml('code/helpers/fam_colors.yml'))
shapes <- c("star","square","circle","triangle", "rectangle")
plot_phy_net <- function(G, main="", dropiso=TRUE) {

  E(G)$color <- ifelse(E(G)$sign>0, "#159C0033", "#FF000033")
  V(G)$color <- col[V(G)$Rank5]
  V(G)$shape <- shapes[as.factor(V(G)$Rank2)]

  if (dropiso) {
    ind <- igraph::V(G)$name[which(igraph::degree(G) < 1)]
    G <- igraph::delete.vertices(G, ind)
  }
  l <- cbind(V(G)$l1, V(G)$l2)

  plot(G, vertex.label.cex=0.4, vertex.label=NA,
      layout=l, vertex.size=4, main=main,
      rescale=TRUE, ylim=c(-1,1), xlim=c(-1,1))
}


### Merge igraph objects for consensus net methods
### keeping or merging relevant vertex/edge attributes
p.equal <- function(x) length(unique(x))==1
gr_intersect <- function(igr_li, rseed=1000) {
  gr <- do.call('intersection', igr_li)

  ## merge edge attributes
  attr <- names(edge.attributes(gr))
  # keep only sign-consistent edges
  sattr <- grep("sign", attr, val=TRUE)
  gr <- delete.edges(
      gr,
      E(gr)[!apply(data.frame(edge.attributes(gr)[sattr]), 1, p.equal)]
    )
  E(gr)$sign <- E(gr)$sign_1
  for (a in sattr) gr <- remove.edge.attribute(gr, a)

  # average weights
  wattr <- grep("weight", attr, val=TRUE)
  E(gr)$weight <- rowMeans(data.frame(edge.attributes(gr)[wattr]))
  for (a in wattr) gr <- remove.edge.attribute(gr, a)

  # remove extra edge attributes
  cattr <- grep("cross", attr, val=TRUE)
  for (a in cattr) gr <- remove.edge.attribute(gr, a)

  ## merge vertex attributes
  attr <- names(vertex.attributes(gr))
  # remove extra vertex attributes (module membership)
  mattr <- grep("cmemb", attr, val=TRUE)
  for (a in mattr) gr <- remove.vertex.attribute(gr, a)

  ## Taxonomy should be consistent across all input graphs
  V(gr)$Rank2 <- V(gr)$Rank2_1
  V(gr)$Rank5 <- V(gr)$Rank5_1
  tattr <- grep("Rank[2,5]_|l[1-2]", attr, val=TRUE)
  for (a in tattr) gr <- remove.vertex.attribute(gr, a)

  ## Recompute layout
  set.seed(rseed)
  l <- layout.fruchterman.reingold(gr)
  V(gr)$l1 <- l[,1]
  V(gr)$l2 <- l[,2]
  gr
}



####
## adding vertex shapes: triangle and star to igraph shapes
## these will be used in plotting networks in addition to the base shapes: circles & squares
###
# functions to create triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }

  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

# generic star vertex shape, with a parameter for number of rays
mystar <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }

  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           symbols(x=x, y=y, bg=bg,
                   stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                   add=TRUE, inches=FALSE)
         })
}
# no clipping, edges will be below the vertices anyway
add_shape("star", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=5))


edge.weights <- function(community, network, weight.within = 1, weight.between = 0.1) {
  bridges <- igraph::crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights)
}
edge.color<- function(community, network, weight.within = rgb(1,0,0,alpha=0.2) #transparent red
                      , weight.between = "grey") {
  bridges <- igraph::crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights)
}

