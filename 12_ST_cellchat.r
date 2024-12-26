dev.off()
rm(list=ls())
dir()
gc()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
options(clusterProfiler.download.method = "wininet")
set.seed(1)
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)


visium = readRDS("ST.rds")
# Prepare input data for CelChat analysis
data.input = GetAssayData(visium.brain, slot = "data", assay = "SCT")
meta = data.frame(labels = Idents(visium.brain), row.names = names(Idents(visium.brain))) 
unique(meta$labels) # check the cell labels

# load spatial imaging information
# Spatial locations of spots from full (NOT high/low) resolution images are required
spatial.locs = GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters of the full resolution images 
scale.factors = jsonlite::fromJSON(txt = file.path("scalefactors_json.json"))
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, 
                     lowres = scale.factors$tissue_lowres_scalef 
)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
pathways.show <- c("THBS") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    layout = "circle")
# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    layout = "spatial", 
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 0)
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# USER can visualize this information on the spatial imaging, e.g., bigger circle indicates larger incoming signaling
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", 
                    edge.width.max = 2, alpha.image = 0.2, 
                    vertex.weight = "incoming", vertex.size.max = 3, vertex.label.cex = 0)
groupSize <- as.numeric(table(cellchat@idents))
mat <- cellchat@net$weight
par(mfrow = c(2,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
mat <- cellchat@net$count
par(mfrow = c(2,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
netAnalysis_contribution(cellchat, signaling = pathways.show)
netVisual_bubble(cellchat, sources.use = c(2,3,4), targets.use = c(1), 
                 signaling = c("THBS"), remove.isolate = FALSE)
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)




























