library(igraph)
library(data.table)
#setwd("C:/Users/dplewi3/Dropbox")
setwd("~/Dropbox")

# make an example graph
# similar to making a matrix
# directed = TRUE is the default
g1 = graph( edges=c(1,2,2,3,3,1), n=3, directed=FALSE)

# quick aside, note that the syntax for this varies depending on what type of object we're converting to a graph
# so we write graph.data.frame(d, directed = FALSE)
# but graph.adjacency(m, "undirected") -- this is short for mode = "undirected"

plot(g1)

# note the difference between

typeof(g1)
# tells us the internal storage mode

class(g1)
# tells us the specific features it has

g1

# what does the same thing with n = 10 produce?
g2 = graph( edges=c(1,2,2,3,3,1), n=10, directed=FALSE)

# 7 isolates

# note in g1 the input is numeric so there are no names
V(g1)$names 
# names is generated automatically as an attribute when we create the igraph object

# can repeat the same process as in g1 with strings, and we get
g1_names = graph( edges=c("A","B","B","C","C","A"), n=3, directed=FALSE)

V(g1_names)$names 

# note the warning message--
# when the vertex names are given by the strings, R assumes no isolates
# so removing n= produces the same result
g1_names = graph( edges=c("A","B","B","C","C","A"), directed=FALSE)

# instead, add isolates manually by naming them in a separate vector
g2_names  = graph( edges=c("A","B","B","C","C","A"), isolates = c("D", "E", "F", "G", "H", "I", "J"), directed=FALSE)

plot(g2_names)

# graph_from_literal
# can enter different types of edges by using specified symbols:
# - undirected (number of - doesn't matter)
# -+ directed right (number of - doesn't matter)
# +- directed left (number of - doesn't matter)
# ++ symmetric (- in the middle don't matter)

plot(graph_from_literal(A---B, B---C))
plot(graph_from_literal(A--+B, B+--C))
plot(graph_from_literal(A+-+B, B+-+C)) 

# : incorporates sets of vertices on either side of the relationship
# make a 5-star
plot(graph_from_literal(a:b:c:d-e))
# make a clique of size 5
plot(graph_from_literal(a:b:c-c:d:e))

# to make an isolate, add a comma
plot(graph_from_literal(a:b:c-c:d:e, f))

# pull the edges of a graph using E(graph)

# without names
E(g2)

# with names 
E(g2_names)

# pull the vertices with V(graph)
# without names
V(g2)

# with names 
V(g2_names)

# can get the matrix out two ways
g2_names_mat = g2[]

# or 
g2_names_mat2 = as_adjacency_matrix(g2_names)

# note what kind of objects these are
class(g2_names_mat)
typeof(g2_names_mat)

class(g2_names_mat2)
typeof(g2_names_mat2)
# sparse matrices, which are okay, but may or may not transfer into other packages or operations

# to get these as regular matrices, need either
g2_names_mat = as.matrix(g2[])

# or
g2_names_mat2 = as_adjacency_matrix(g2_names, sparse = FALSE)

# can subset to check out and individual row or column of the matrix with the index method
g2_names[1,]
g2_names[,1]

# we can add attributes directly to the graph using the V() notation and a vector

# with strings
V(g2_names)$gender = rep(c("male", "female"), 5)

# we'll get a warning message if the vector does not match the number of nodes
V(g2_names)$gender = rep(c("male", "female", "other"), 3)
# and will recycle the remainder


# with numbers
V(g2_names)$age = sample(100, 10, replace = FALSE)

# special attributes to take note of include type and weight 

# type lets us know what type of node (type of actor) or edge (type of relationship) this is

# note, when projecting a two-mode network down to a one-mode, a type attribute is required
V(g2_names)$type = "person"

# define as an communication network
V(g2_names)$type = "talks to"

# pull all vertex attributes with
vertex_attr(g2_names)
edge_attr(g2_names)

# can also set graph-level attributes, like names
g2_names = set_graph_attr(g2_names, "name", "Communication Network")

graph_attr(g2_names)

# or any other characteristics
g2_names = set_graph_attr(g2_names, "time period", "November 2017")

# what attributes have we put in the graph
graph_attr_names(g2_names)

# and call them 
graph_attr(g2_names)

graph_attr(g2_names, "time period")

# can delete graph attributes
g2_names = delete_graph_attr(g2_names, "time period")

graph_attr(g2_names)

# as well as other attributes
g2_names = delete_edge_attr(g2_names, "type")
g2_names = delete_vertex_attr(g2_names, "age")


# suppose we have a graph with some self loops and multiple ties, and some other attributes we want to ignore
g3 = add.edges(g2_names, c("A", "A", "A", "B"))
V(g3)$gender = rep(c("male", "female"), 5)

# so say that we set the colors using...

# short aside on setting colors, using as.factor to take advantage of R's indexing--
# this avoids having to set colors manually with something like
colors = V(g3)$gender
colors[colors == "male"] = "light blue"
colors[colors == "female"] = "red"

# and then 
V(g3)$color = colors
# or
plot(g3, vertex.color = colors)

# 1+ notation achieves the same thing as as.factor here--both create a vector of 1s and 2s
identical(1+(V(g3)$gender=="male"), as.numeric(as.factor(V(g3)$gender)))

# if we wanted to figure out a value that was throwing us off, this would be a good way to check for it
which(1+(V(g3)$gender=="male") != as.numeric(as.factor(V(g3)$gender)))


# so we get a bunch of possible ways to do this with the indexing
# color is a special attribute that corresponds to what will be shown in the plot (like weight, more on this below)
V(g3)$color = c("light blue", "red")[as.factor(V(g3)$gender)]

#or adding it as a vertex attribute this way
set.vertex.attribute(g3, "color", c("light blue", "red")[as.factor(V(g3)$gender)])

# or adding it through vertex.color in the plot
plot(g3, vertex.color = c("light blue", "red")[as.factor(V(g3)$gender)])

# or creating a specific color object like above
color = c("light blue", "red")[as.factor(V(g3)$gender)]
V(g3)$color = color

# or calling the color object in the plotting function
plot(g3, vertex.color = color)
# you can start to see how all of these things work in combination with each other, can choose which option is more expedient for the situation -- for example how big is the data, do you have another indexing object that's already arranged in the right order, e.g., "genders", and so on

# okay now that we have the colors in, let's begin to combine attributes to get rid of unwanted information and concatenate some other information
plot(g3)

# remember
E(g3)$weight 

# now
g3_simple = simplify(g3, remove.multiple = TRUE, remove.loops = TRUE)
plot(g3_simple)

# but we could also specify this differently
E(g3)$weight = 1

g3_simple2 = simplify(g3, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = c(weight = "sum"))
E(g3_simple2)$weight

# set to one before because summing to a NULL will also result in a null

# calling an igraph object provides us with some information about the network
g3_simple2

# "U" indicates the graph is undirected "D" for directed
# "N" indicates that the vertices have names (will be missing otherwise)
# "W" indicates that the edges have weights (will be missing otherwise)
# "B" indicates that the graph is bipartite

# we also get information about attributes
#(g/c) - graph-level attribute
#(v/c) - vertex-level attribute
#(e/n) - edge-level attribute
# the c or n indicates whether this is a character or numeric attribute


# we can make a variety of graphs using the make_ family of functions
# helful for looking at some important properties of networks of particular types, and then checking to see the degree to which these properties hold true in the actual data

# all isolates
empty = make_empty_graph(30)
plot(empty)

# all connected
full = make_full_graph(30)
plot(full)

# everyone else connected to one
# for star, we drop the _graph
star = make_star(30)
plot(star)

# everyone connected to two
ring = make_ring(30)
plot(ring)

# note that we can turn a ring into a line with
line = delete.edges(ring, sample(30, 1))
plot(line)

# or n lines with
n = sample(30, 1) - 1
lines = delete.edges(ring, sample(30, n, replace = FALSE))
plot(lines)

# in a tree graph, each vertex is connected by exactly one path
# can make a tree with make_tree
# need to specify the number of branches from the central node, here short for children = 3 
# default is directed
tree = make_tree(30, 3, "undirected")
plot(tree)

# a forest is a tree graph with disconnected components
# get a forest similarly as a line
forest = delete.edges(tree, sample(30, n, replace = FALSE))

# small-world graph made up of tightly connected local neighborhoods and bridging shortcuts
# dim = how many dimensions should neighborhoods be laid out in, size = how big is the network on each dimension, nei = setting where the neighborhood will be defined, p = probability of being a shortcut bridge
small = sample_smallworld(dim=2, size=10, nei=1, p=0.1)

plot(small)
# can't see this, try
plot(small, layout = layout_in_circle)


# igraph also has some historical graphs, like the karate club from the book
karate = graph("Zachary") 

 plot(karate)

# rewiring graphs randomly with a certain probability (similar to producing random graphs with ergm)
ring_re = rewire(ring, each_edge(prob=0.1))
plot(ring_re)

# rewiring with adding connections, where people within a certain distance to one another make connections
ring_nei = connect.neighborhood(rn, 5)
plot(ring_nei)

# can combine graphs with %du%
comb = ring_nei %du% tree
plot(comb)
# this ignores vertex names

# if we wanted to let igraph know these were the same vertices in either graph
comb2 = union(ring_nei, tree)
plot(comb2)

# aside, if we only wanted the ties that existed in both networks, e.g., exercise 1 extra question
comb3 = intersection(ring_nei, tree)
plot(comb3)


# reading in network data from a file 

# reading in a data frame of edges
# note the differences between
nodes = read.csv("Dataset1-Media-Example-NODES.csv", header=TRUE)
typeof(nodes$media)

# and
nodes = fread("Dataset1-Media-Example-NODES.csv", header=TRUE)
typeof(nodes$media)

# so you could specify 
nodes = read.csv("C:/Users/akabo/Downloads/social network analysis/Dataset1-Media-Example-NODES.csv", header=TRUE, stringsAsFactors = FALSE)
typeof(nodes$media)

# or
nodes = read.csv("Dataset1-Media-Example-NODES.csv", header=TRUE, as.is = TRUE)
typeof(nodes$media)

links = read.csv("C:/Users/akabo/Downloads/social network analysis/Dataset1-Media-Example-EDGES.csv", header=TRUE, as.is = TRUE)

# can check for duplicate links
nrow(links) - nrow(unique(links[,c("from", "to")]))

# note that if this were an undirected graph we would also have to compare the links in the reverse order
# so might try something like this
nrow(links)*2 - nrow(rbind(unique(links[,c("from", "to")]), unique(links[,c("to", "from")])))

# but note that data frame is sorting our from column directly to the from column from above, even though we indicated it as the second column
# so 
unique(links[,c("to", "from")])
# is as intended
# but 
rbind(unique(links[,c("from", "to")]), unique(links[,c("to", "from")]))
# ignores the ordering

# sometimes the sorting could be useful, but not in this case
# note that removing the names doesn't actually work
colnames(links) = NULL
rbind(unique(links[,1:2]), unique(links[,2:1]))
# need to specify new matching names for each object for the corresponding column

# so from the original
links = read.csv("C:/Users/akabo/Downloads/social network analysis/Dataset1-Media-Example-EDGES.csv", header=TRUE, as.is = TRUE)
u = unique(links[,1:2])
u_rev = unique(links[,2:1])

colnames(u_rev) = colnames(u)
rbind(u, u_rev)

# okay

# if these were data tables this could be slightly avoided with rbindlist, see empirical exercise #2 answer key

# take a few minutes to discuss-- how could you use this method to check for reciprocity in this edge list?

# why does this look odd?

# avoiding this issue
# collapse on the other variables besides the weight in oder to retain info about the weight
# this lets us retain info about the types of links
links = aggregate(links[,3], links[,-3], sum)

# some ordering to make things easier to view if we need to
links = links[order(links$from, links$to),]

links

# name x to be the weight again
colnames(links)[4] = "weight"
# removing the row names will make sure igraph doesn't read these as vertex names
rownames(links) = NULL

# now reading in the data to igraph
# vertices = nodes automatically reads in the attributes of the nodes in the edge list by name
net1 =graph.data.frame(links)
# check edges and attributes
net1

# access different parts of the network object
E(net1)
V(net1)

# and attributes specifically
E(net1)$type
V(net1)$media

plot(net1)

# if we don't want self loops...
# note the difference between
net1 = simplify(net1, remove.multiple = FALSE, remove.loops = TRUE) 
#net1 = simplify(net, remove.multiple = TRUE, edge.attr.comb=list(weight="sum")
# in the second case, we're actually removing information about the different link types that we might want later on

# if you want to check info from igraph in a different format
as_edgelist(net1)

# or
get.edgelist(net1)

# to get a (sparse) matrix out, 
as_adjacency_matrix(net1)
get.adjacency(net1)

# can also get the attributes table back out
as_data_frame(net1, what="vertices")


# now, read in some affiliation data
nodes2 = read.csv("C:/Users/akabo/Downloads/social network analysis/Dataset2-Media-User-Example-NODES.csv", header=TRUE, as.is=TRUE)
links2 = read.csv("C:/Users/akabo/Downloads/social network analysis/Dataset2-Media-User-Example-EDGES.csv", header=TRUE, row.names=1)

# note 
dim(links2)

# so if we try 
net2 = graph.adjacency(links2)
# we get the sqaure matrix error

# instead, need to use
net2 = graph.incidence(links2)

# now check
V(net2)$type

# so igraph understands that there are two sets of nodes that don't have connections with each other, that have different types

# can get co-membership projection directly from igraph
net2.bp = bipartite.projection(net2)

# which will be two networks in a list
net2.bp

# plot co-membership using the name of each user or media type
plot(net2.bp[[1]], vertex.label=nodes2$media[is.na(nodes2$media.type)])
plot(net2.bp[[2]], vertex.label=nodes2$media[!is.na(nodes2$media.type)])


# our plot parameters for nodes all start with vertex. and the ones for edges start with edge.
# many options are listed under ?igraph.plotting

# some useful ones below

# NODES	 
# vertex.color	 Node color
# vertex.frame.color	 Node border color
# vertex.shape	 One of “none”, “circle”, “square”, “csquare”, “rectangle”
#  “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
# vertex.size	 Size of the node (default is 15)
# vertex.size2	 The second size of the node (e.g. for a rectangle)
# vertex.label	 Character vector used to label the nodes
# vertex.label.family	 Font family of the label (e.g.“Times”, “Helvetica”)
# vertex.label.font	 Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
# vertex.label.cex	 Font size (multiplication factor, device-dependent)
# vertex.label.dist	 Distance between the label and the vertex
# vertex.label.degree	 The position of the label in relation to the vertex,
#  where 0 right, “pi” is left, “pi/2” is below, and “-pi/2” is above
# EDGES	 
# edge.color	 Edge color
# edge.width	 Edge width, defaults to 1
# edge.arrow.size	 Arrow size, defaults to 1
# edge.arrow.width	 Arrow width, defaults to 1
# edge.lty	 Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”,
#  3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
# edge.label	 Character vector used to label edges
# edge.label.family	 Font family of the label (e.g.“Times”, “Helvetica”)
# edge.label.font	 Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
# edge.label.cex	 Font size for edge labels
# edge.curved	 Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
# arrow.mode	 Vector specifying whether edges should have arrows,
#  possible values: 0 no arrow, 1 back, 2 forward, 3 both
# OTHER	 
# margin	 Empty space margins around the plot, vector with length 4
# frame	 if TRUE, the plot will be framed
# main	 If set, adds a title to the plot
# sub	 If set, adds a subtitle to the plot

# take a few minutes to work and discuss: using the commands above, generate a plot for the first type in which
# all three media types have a unique color
# the node size is determined by the size of the audience
# there are no labels for nodes
# edge width is determined by weight


# adding a legend
# x and y give the positions, next the names on the legend, pch is the shape of the symbols, col is the color of the symbols, next is the color of the symbol background, next is size for the symbols, then size for the text, n = no box, ncol gives a vertical legend with 1 column

# note that the pch 21 with background similar to using pch with the col = colors--the outlined circle just lets us give them an outline
legend(x=-1.5, y=-1.1, c("Newspaper","Television", "Online News"), pch=21, col="#777777", pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1)
 
# take a few minutes to build a plot with just the media names

# we can also build a plot where the edges are colored by the type of source that they come from 
# the ends()function lets us break apart the edge origin and destinations
ends(net, es=E(net1), names=FALSE)[,1]

# so if we want just the origin
origin = ends(net1, es=E(net1), names=FALSE)[,1]

E(net1)$color = V(net1)$color[origin]

plot(net1, edge.curved=.1)  


# exploring different plot layouts

# use a sample graph generated from a preferntial attachment model (few preferred hubs that many others attach to, this model approximates the internet well)
pref = sample_pa(80) 

# so we can see everything a bit easier
V(pref)$size = 8
V(pref)$frame.color = "white"
V(pref)$color = "orange"
V(pref)$label = "" 
E(pref)$arrow.mode = 0

plot(pref)

# some layouts

# random
plot(pref, layout = layout.random)

# or 
plot(pref, layout = layout_randomly)

# you can also assign the layout as an object in advance
l = layout_in_circle
plot(pref, layout = l)

l = layout_on_sphere
plot(pref, layout = l)

# l gives the coordinates on the plot where each vertex should be

# so we could specify these on our own too
# e.g., 
l = cbind(1:vcount(pref), c(1, vcount(pref):2))
l
plot(pref, layout = l)

# a note our most used algorithm so far, fruchterman-reingold 
# this algorithm falls into the family of force-directed layouts
# another common force-directed method similar to fruchterman-reingold is kamada-kawai

# "Force-directed layouts try to get a nice-looking graph where edges are similar in length and cross each other as little as possible. They simulate the graph as a physical system. Nodes are electrically charged particles that repulse each other when they get too close. The edges act as springs that attract connected nodes closer together. As a result, nodes are evenly distributed through the chart area, and the layout is intuitive in that nodes which share more connections are closer to each other. The disadvantage of these algorithms is that they are rather slow and therefore less often used in graphs larger than ~1000 vertices. You can set the “weight” parameter which increases the attraction forces among nodes connected by heavier edges."

# fruchterman-reingold = nearby vertices attract and all vertices repel
# kamada-kawai = all vertices repel, but the strength of the repelling force is equal to the graph distance of the nodes

# compare on the same viewer
par(mfrow=c(1,2), mar = c(0,0,0,0)) 
plot(pref, layout = layout.fruchterman.reingold)
plot(pref, layout = layout.kamada.kawai)
dev.off()

# if we run fruchterman-reingold a few times, we'll notice that the layouts are different
par(mfrow=c(2,2), mar = c(0,0,0,0)) 
plot(pref, layout = layout.fruchterman.reingold)
plot(pref, layout = layout.fruchterman.reingold)
plot(pref, layout = layout.fruchterman.reingold)
plot(pref, layout = layout.fruchterman.reingold)
dev.off()

# but, there are situations where we would want the the layouts to stay the same. take a minute to discuss a few


# looking at the graph as characteristics of nodes change over time, or different types of relationships

# we can fix this by determining the layout beforehand with l 
l = layout.fruchterman.reingold(pref)
# layout_with sytax also works here

# now
par(mfrow=c(2,2), mar = c(0,0,0,0)) 
plot(pref, layout = l)
plot(pref, layout = l)
plot(pref, layout = l)
plot(pref, layout = l)


# we can also rescale the coordinates manually by turning rescale off
l = layout_with_kk(pref)

# R automatically rescales to [-1,1], so replicate this first with norm_cords
l = norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

par(mfrow=c(2,2), mar = c(0,0,0,0)) 
plot(pref, rescale=FALSE, layout=l*0.4)
plot(pref, rescale=FALSE, layout=l*0.6)
plot(pref, rescale=FALSE, layout=l*0.8)
plot(pref, rescale=FALSE, layout=l*1.0)

# lgl http://lgl.sourceforge.net is meant for large, connected graphs
# you can specify a root to be in the middle of the graph
l = lapply(1:4, function(i) layout_with_lgl(pref, root = sample(vcount(pref), 1))) 
lapply(1:4, function(i) plot(pref, layout = l[[i]]))

# to compare some more layouts, pull them from the igraph package space
layouts = grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 

# check these out at http://igraph.org/r/doc/layout_.html

# remove layouts that do not apply to our graph, e.g., it's not two-mode and has more than exactly one path in between nodes
layouts = layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]

par(mfrow=c(5,3), mar=c(1,1,1,1))
for (layout in layouts) {
  l = do.call(layout, list(pref)) 
  plot(pref, edge.arrow.mode=0, layout=l, main=layout) }

# for the original graph
for (layout in layouts) {
  l = do.call(layout, list(net1)) 
  plot(net1, edge.arrow.mode=0, layout=l, main=layout) }


# we can try to get more information out of the graph by viewing it differently, for example, taking weaker strength links out
# resetting the device parameters
dev.off()
sparse = delete_edges(net1, E(net1)[weight<mean(links$weight)])

plot(sparse) 

# next, take a few minutes and try to plot the two types of ties, hyperlink or mention, separately, on the same graphical device


# we can also plot interactively and manually with tkplot -- similar to neo4j in the visualization view

# first, open up the plot
tkid = tkplot(net1) 
#tkid is the id of the tkplot that will open

# take a few minutes to cluster the vertices together by color, then export the coordinates you chose from tkplot and import them to igraph by setting them to a layout object
l = tkplot.getcoords(tkid) 

# close the tkplot and graph the plot you made in r
tk_close(tkid, window.close = T)
plot(net1, layout=l)


# making a heatmap to represent the network, instead of a regular graph
netm = get.adjacency(net1, attr="weight", sparse=FALSE)

colnames(netm) = V(net1)$media
rownames(netm) = colnames(netm) 

colors = colorRampPalette(c("gold", "dark orange")) 

heatmap(netm[,17:1], Rowv = NA, Colv = NA, col = colors(100), scale="none", margins=c(10,10) )

# we can get some clustering on here too
heatmap(netm[,17:1], col = colors(100), scale="none", margins=c(10,10) )


# some plotting options for two-mode graphs

# not too informative
plot(net2, vertex.label = nodes2$media)

# we can use layout.bipartite
plot(net2, vertex.label = nodes2$media, layout = layout.bipartite)
# but this isn't great either

# take a few minutes to work on the bipartite make a graph that has 
# different colors for media outlets and users
# different shapes for media outlets and users
# only labels for the media outlets

# this is better now too
plot(net2, vertex.label=NA, vertex.size=7, layout=layout_as_bipartite) 

# text as nodes also helpful here
plot(net2, vertex.shape="none", vertex.label=nodes2$media, vertex.label.color=V(net2)$color, vertex.label.font=2.5,vertex.label.cex=.6, edge.color="gray70",  edge.width=2)


# some descriptive stats in the network
edge_density(net1, loops=FALSE)

# equal to (for directed)
ecount(net1)/(vcount(net1)*(vcount(net1)-1))

# mutuality/reciprocity
reciprocity(net1)
# mutual, asymmetric, and null node pairs
dyad_census(net1)
2*dyad_census(net1)$mut/ecount(net1)

# transitivity/closure
# global measure is a ratio ratio of triangles (cliques of 3, direction disregarded) to connected triples (i.e., sets of i-j-k)
# individual measure is for each node -- how many triangles is a node part of as a fraction of triples it is a part of

transitivity(net1, type="global") 
# net1 is treated as an undirected network
transitivity(as.undirected(net1, mode="collapse"))
# same as above
transitivity(net1, type="local")
triad_census(net1) 

# how many triads of each type -- only for directed networks 
# A,B,C, empty triad.
# A->B, C, triad with a single directed edge.
# A<->B, C, triad with a reciprocated connection between two vertices.
# A<-B->C, triadic out-star.
# A->B<-C triadic in-star.
# A->B->C, directed line.
# A<->B<-C.
# A<->B->C.
# A->B<-C, A->C.
# A<-B<-C, A->C.
# A<->B<->C.
# A<-B->C, A<->C.
# A->B<-C, A<->C.
# A->B->C, A<->C.
# A->B<->C, A<->C.
# A<->B<->C, A<->C,  complete triad.

# diameter, max geodesic distance in network
diameter(net1, directed=FALSE, weights=NA)
diameter(net1, directed=FALSE)

diam = get_diameter(net1, directed=TRUE)

diam

# returns a vertex sequence

# note though that when asked to behaved as a vector, a vertex sequence will produce the numeric indexes of the nodes in it
# same applies for edge sequences

class(diam)
as.vector(diam)
# okay

# make a plot coloring nodes along the diameter:

vcol = rep("gray40", vcount(net1))
vcol[diam] = "gold"
E(net1)$color = "gray80"
E(net1)[E(net1, path=diam)]$color = "orange" 

# E(net1, path=diam) finds edges along a path, here 'diam'
plot(net1, vertex.color=vcol, edge.arrow.mode=0)
# may be a little hard to see

# sizing by degree

plot(net1, vertex.size=degree(net1, mode = "all")*3)

# sizing by betweenness
plot(net1, vertex.size= betweenness(net1)/3)

# histograph of degree
deg = degree(net1, mode = "all")
hist(deg, breaks=1:vcount(net1)-1, main="Histogram of node degree")

# plot of degree distribution
deg.dist = degree_distribution(net1, cumulative=TRUE, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")

# identifying hubs and authorities

# note on h/a model
# The hubs and authorities algorithm developed by Jon Kleinberg was initially used to examine web pages. Hubs were expected to contain catalogs with a large number of outgoing links; while authorities would get many incoming links from hubs, presumably because of their high-quality relevant information.

hs = hub_score(net1, weights=NA)$vector
as = authority_score(net1, weights=NA)$vector

par(mfrow=c(1,2), mar = c(0,0,0,0))
plot(net1, vertex.size=hs*50, main="Hubs")
plot(net1, vertex.size=as*30, main="Authorities")


# incorporating information on paths

# finding shortest paths through nodes and plotting them 
news.path = shortest_paths(net1, from = V(net1)[media=="MSNBC"], to  = V(net1)[media=="New York Post"],output = "both")
# both path nodes and edges

# generate edge color variable to plot the path:
ecol = rep("gray80", ecount(net1))
ecol[unlist(news.path$epath)] = "orange"

# generate edge width variable to plot the path:
ew = rep(2, ecount(net1))
ew[unlist(news.path$epath)] = 4

# generate node color variable to plot the path:
vcol = rep("gray40", vcount(net1))
vcol[unlist(news.path$vpath)] = "gold"
dev.off()
plot(net1, vertex.color=vcol, edge.color=ecol, edge.width=ew, edge.arrow.mode=0)


# what about edges coming out of a specific node
inc.edges = incident(net1,  V(net1)[media=="Wall Street Journal"], mode="all")
# for multiple nodes, use indident edges

# set colors to plot the selected edges
ecol = rep("gray80", ecount(net1))
ecol[inc.edges]  = "orange"
vcol  = rep("grey40", vcount(net1))
vcol[V(net1)$media=="Wall Street Journal"] = "gold"
plot(net1, vertex.color=vcol, edge.color=ecol)

# what about the nodes connected to a specific node
neigh.nodes = neighbors(net1, V(net1)[media=="Wall Street Journal"], mode="out")
# set colors to plot the neighbors:
vcol[neigh.nodes] <- "#ff9d00"
plot(net1, vertex.color=vcol)


# some more operators to subset nodes
#special operators for the indexing of edge sequences: %–%, %->%, %<-%
# E(network)[X %–% Y] selects edges between vertex sets X and Y, ignoring direction
# E(network)[X %->% Y] selects edges from vertex sets X to vertex set Y
# E(network)[X %->% Y] selects edges from vertex sets Y to vertex set X

# newspaper to online
E(net1)[ V(net1)[type.label=="Newspaper"] %->% V(net1)[type.label=="Online"] ]

# cocitiation, how many shared nominations two nodes have (source points to both nodes)
cocitation(net1)

# some group and community identification

# first make undirected
net_un =as.undirected(net1, mode= "collapse", edge.attr.comb=list(weight="sum", "ignore"))

# check out cliques
cliques(net_un)

# list of cliques       
sapply(cliques(net_un), length)

# clique sizes
largest_cliques(net_un)
# cliques with max number of nodes
vcol = rep("grey80", vcount(net_un))
vcol[unlist(largest_cliques(net_un))] = "gold"
plot(as.undirected(net_un), vertex.label=V(net_un)$name, vertex.color=vcol)

# communities based on propagating labels
# note on approach
# Assigns node labels, randomizes, than replaces each vertex’s label with the label that appears most frequently among neighbors. Those steps are repeated until each vertex has the most common label of its neighbors.
# so, similar process to k-means clustering

clp = cluster_label_prop(net1)

plot(clp, net1)

# lustering on modularity, e.g., the spatial sense in which a graph can be broken down into independent clusters
cfg = cluster_fast_greedy(as.undirected(net1))

plot(cfg, as.undirected(net1))

# we can also assign our own colors based on the community detection
V(net1)$community = cfg$membership
colrs = adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(net1, vertex.color=colrs[V(net1)$community])

# finally, we can explore the degree of homophily or sorting on some variable
assortativity_nominal(net1, V(net1)$media.type, directed=FALSE)
assortativity(net1, V(net1)$audience.size, directed=FALSE)
assortativity_degree(net1, directed=FALSE)

# igraph exercise 
# there are two .rda files containing the MBA advice and trust networks from exercise 1
# this time, the data have demographic attributes attached to them

# using the network data and some attributes, produce some plots and reports to explain why the networks looks the way it does, using the techniques provided above, and whatever other techniques you find appropriate
load("trust_att.rda")
load("advice_att.rda")
