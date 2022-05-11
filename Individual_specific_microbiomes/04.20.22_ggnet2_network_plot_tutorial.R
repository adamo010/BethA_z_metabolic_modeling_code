#the goal of this script is to run through a tutorial for ggnet2
#we'll follow along with the vignette here: https://mran.microsoft.com/snapshot/2018-01-02/web/packages/ggCompNet/vignettes/examples-from-paper.html

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/") 

#these are the libraries we need: 
library(dplyr)
library(tidyr)
library(ggplot2) # needs to be version ≥ 2.1.0
library(scales)

#install some new packages:
#ggnet2
if (!require(GGally, quietly = TRUE)) {
  getFromNamespace("install_github", asNamespace("devtools"))("ggobi/ggally")
}
#geom_net
if (!require(geomnet, quietly = TRUE) ||
    packageVersion("geomnet") < "0.2.0") {
  getFromNamespace("install_github", asNamespace("devtools"))("sctyner/geomnet")
}
#ggnetwork
if (!require(ggnetwork, quietly = TRUE) ||
    packageVersion("ggnetwork") < "0.5.1") {
  getFromNamespace("install_github", asNamespace("devtools"))("briatte/ggnetwork")
}
#Load in requried libraries post-installation
library(network)
library(sna)
library(GGally)
library(geomnet)
library(ggnetwork)
library(igraph)

###########Importing data
# make data accessible for geomnet
data(madmen, package = "geomnet")

#data step: merge edges and nodes by the "from" column
MMnet <- fortify(as.edgedf(madmen$edges), madmen$vertices) #I think this is just pulling out a couple of columns?
#in my data, from_id would be microbe, and to_id would be host gene?

## Joining edge and node information by from_id and label respectively.
#I"m getting a weird error about "The first two columns of `x` must be of the same type." so I guess I'll set the type of all MMnet
str(MMnet) #all the columns are factors... is that bad?
library(dplyr)
MMnet<- MMnet %>%
  mutate(across(everything(), as.character))

# create plot
set.seed(10052016)
ggplot(data = MMnet, aes(from_id = from_id, to_id = to_id)) +
  geom_net(aes(color = Gender), layout.alg = "kamadakawai", 
           size = 2, labelon = TRUE, vjust = -0.6, ecolour = "grey60",
           directed =FALSE, fontsize = 3, ealpha = 0.5) +
  scale_colour_manual(values = c("#FF69B4", "#0099ff")) +
  xlim(c(-0.05, 1.05)) +
  theme_net() +
  theme(legend.position = "bottom")

#there's some error here that I don't know how to fix.
#let's try a different package: ggnet2
library(GGally)
library(network)

# random graph
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:10]
ggnet2(net)

#The net argument is the only compulsory argument of ggnet2.
#It can be a network object or any object that can be coerced to that class through its edgeset.constructors functions,
#such as adjacency matrixes, incidence matrixes and edge lists.

#let's see if we can do that with MMnet data
#create undirected network
#first, drop na rows
MMnet2 <- MMnet %>% drop_na()
#then, convert everything to strings
MMnet2 <- MMnet2 %>%
  mutate(across(everything(), as.character))
#NOW we can create a network:
mm.net <- network::network(MMnet2[, 1:2], directed = FALSE) #create an undirected network from the first two columns of dataframe

# create node attribute (gender)
rownames(madmen$vertices) <- madmen$vertices$label
mm.net %v% "gender" <- as.character(
  madmen$vertices[ network.vertex.names(mm.net), "Gender"]
)

# gender color palette (don't ask me, it was in the vingette)
mm.col <- c("female" = "#ff69b4", "male" = "#0099ff")

# create plot for ggnet2
set.seed(10052016)
ggnet2(mm.net, color = mm.col[ mm.net %v% "gender" ],
       label = TRUE, label.color = mm.col[ mm.net %v% "gender" ],
       size = 2, vjust = -0.6, mode = "kamadakawai", label.size = 3)

#okay, that works pretty nice.

#how can I adjust this so the line widths represent the strength of association?
#let's try a different dataset that has something like that in it. 
data(bikes, package = 'geomnet')
# data step for geomnet
tripnet <- fortify(as.edgedf(bikes$trips), bikes$stations[,c(2,1,3:5)])
#okay. Looks like it's as simple as having extra columns. The key here seems to be to create a dataframe where
#the first two columns represent the from and to nodes (labeled from_id and to_id), and the other columns can be
#named whatever. They'll be included in the network design. 
#also, NO NANs allowed. 
#also also, your to and from id columns need to be characters. Otherwise, 
#you get the dreaded "Error: The first two columns of `x` must be of the same type."

# data preparation for ggnet2 and ggnetwork
#first, drop nans
tripnet <- tripnet %>% drop_na()
#then, convert everything to strings
tripnet <- tripnet %>%
  mutate_at(c("from_id", "to_id"), as.character)
#then, create network
bikes.net <- network::network(tripnet[, 1:2], directed = TRUE) #set directed to TURE or will get duplicate edges error
#what the hell is being duplicated here? I'm getting a 'x' contains parallel edges error
#duplicated(tripnet) #NOTHING IS DUPLICATED????

# create edge attribute (number of trips)
network::set.edge.attribute(bikes.net, "n", bikes$trips[, 3 ] / 15)
# create vertex attribute for Metro Station
bikes.net %v% "station" <-  grepl("Metro", network.vertex.names(bikes.net)) #creates a Station column set as True if "Metro" in network_vertex_names
bikes.net %v% "station" <-  1 + as.integer(bikes.net %v% "station") #converts True to 2, and False to 1
rownames(bikes$stations) <- bikes$stations$name

head(ggnetwork(bikes.net))

# create node attributes (coordinates)
bikes.net %v% "lon" <-
  bikes$stations[ network.vertex.names(bikes.net), "long" ]
bikes.net %v% "lat" <-
  bikes$stations[ network.vertex.names(bikes.net), "lat" ]
bikes.col <- c("grey40", "darkorange")

# Non-geographic placement
set.seed(1232016)
ggnet2(bikes.net, mode = "fruchtermanreingold", size = 3, label = TRUE,
       vjust = -0.5, edge.size = "n", layout.exp = 1.5,
       color = bikes.col[ bikes.net %v% "station" ],
       label.color = bikes.col[ bikes.net %v% "station" ])

#Nice!let's see if we can add some aesthetics
ggnet2(bikes.net, mode = "fruchtermanreingold", size = 3, label = TRUE,
       vjust = -0.5, edge.size = "n", layout.exp = 1.5,
       color = bikes.col[ bikes.net %v% "station" ],
       label.color = bikes.col[ bikes.net %v% "station" ]) +
  geom_nodelabel_repel(aes(label = network.vertex.names, color = bikes.col), 
                       box.padding = unit(1, "lines"))

#I do not understand how this system is working. 
#let's try a different approach
set.seed(1232016)
ggplot(data = ggnetwork(bikes.net, layout.alg = "fruchtermanreingold"),
       aes(x, y, xend = xend, yend = yend)) +
  geom_edges(aes(size = n), color = "grey40") +
  geom_nodes(aes(color = factor(station)), size = 4) +
  #geom_nodetext_repel(aes(label = vertex.names, color = factor(station)),
                #vjust = -0.5) +
  #geom_nodelabel_repel(aes(label = vertex.names, color = factor(station))) +
  scale_size_continuous("Trips", breaks = c(2, 4, 6), labels = c(30, 60, 90)) +
  scale_colour_manual("Metro station", labels = c("FALSE", "TRUE"),
                      values = c("grey40", "darkorange")) +
  theme_blank() +
  theme(legend.position = "bottom", legend.box = "horizontal")

#you know what, maybe we don't worry about this too much. It shouldn't be a big issue . 

ggnet2(bikes.net, mode = "fruchtermanreingold", size = 3, label = TRUE,
       vjust = -0.5, edge.size = "n", layout.exp = 1.5,
       color = bikes.col[ bikes.net %v% "station" ],
       label.color = bikes.col[ bikes.net %v% "station" ])


#this option keeps the node labels on top of the nodes
ggplot(data= ggnetwork(bikes.net, layout.alg = "fruchtermanreingold"),
       aes(x, y, xend=xend, yend=yend)) +
  geom_edges(aes(size = n), color = "blue") + #sets edge widths to n
  geom_nodes(aes(color = factor(station)), size = 2) + #sets color to station, which is a factor we created above. 
  scale_colour_manual("Metro station", labels = c("FALSE", "TRUE"),
                      values = c("grey40", "darkorange")) +
  scale_size_continuous("Trips", breaks = c(2, 4, 6), labels = c(30, 60, 90)) +
  theme_blank() +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  geom_nodelabel(aes(color = factor(station), label = vertex.names))
  #geom_nodelabel_repel(aes(color = factor(station), label= vertex.names),
                       #box.padding = unit(1, "lines")) 

#this option scoots the nodes outward, so they're connected to (but not inside) the nodes
  ggplot(data= ggnetwork(bikes.net, layout.alg = "fruchtermanreingold"),
         aes(x, y, xend=xend, yend=yend)) +
    geom_edges(aes(size = n), color = "blue") + #sets edge widths to n
    geom_nodes(aes(color = factor(station)), size = 2) + #sets color to station, which is a factor we created above. 
    scale_colour_manual("Metro station", labels = c("FALSE", "TRUE"),
                        values = c("grey40", "darkorange")) +
    scale_size_continuous("Trips", breaks = c(2, 4, 6), labels = c(30, 60, 90)) +
    theme_blank() +
    theme(legend.position = "bottom", legend.box = "horizontal") +
    #geom_nodelabel(aes(color = factor(station), label = vertex.names))
  geom_nodelabel_repel(aes(color = factor(station), label= vertex.names),
                       box.padding = unit(1, "lines")) 
  
  







