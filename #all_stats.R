library(metabo.pwy)
library(ggplot2)
mmu.net <- fread(system.file("extdata", "metacyc_All_reactions_of_M._musculus.csv", package = "metabo.pwy"))
mmu.cpd <- unique(unlist(str_split(mmu.net$`Substrates (reactants and products) of a reaction`, " // ")))
currency.cpd <- c("PROTON", "WATER", "ATP", "NAD", "OXYGEN-MOLECULE",
                  "NADPH", "NADP", "NADH", "NADH-P-OR-NOP", "CO-A", "Pi",
                  "Donor-H2","NAD-P-OR-NOP","Acceptor","PPI","CARBON-DIOXIDE",
                  "ADP","AMP")

# METACYC DB
db <- load.db(source.db = "metacyc")
db <- reset.ppm.range(db, 10)
#db <- filter.adduct.ions(db)
db@met.db <- db@met.db[cpd.id %in% mmu.cpd]
db@met.adducts$pos <- db@met.adducts$pos[neutral.formula %in% db@met.db$neutral.formula]
db@met.adducts$neg <- db@met.adducts$neg[neutral.formula %in% db@met.db$neutral.formula]

#
# SET2 - RP
#
s2_liver_rp <- init.msexp(name        = "Set2 - Liver - RP",
                          file.pos    = "data/for.metabo.pwy/Set2_liver_RP_pos.csv",
                          file.neg    = "data/for.metabo.pwy/Set2_liver_RP_neg.csv", 
                          sample.info = "data/for.metabo.pwy/Set2_liver_RP.sample_info.csv"
)
#s2_liver_rp <- scale.intensities(s2_liver_rp, method = "median.topn", topn = 75)
s2_liver_rp <- annotate.mz.peaks(s2_liver_rp, db)
s2_liver_rp <- peak.group.comparison(s2_liver_rp,
                                     colname = "genotype",
                                     grpA = "B6", grpB = "FVB", do.lme = F)

#
# SET2 - HILIC
#
s2_liver_hi <- init.msexp(name        = "Set2 - Liver - HI",
                          file.pos    = "data/for.metabo.pwy/Set2_liver_HILIC_pos.csv",
                          file.neg    = "data/for.metabo.pwy/Set2_liver_HILIC_neg.csv", 
                          sample.info = "data/for.metabo.pwy/Set2_liver_HILIC.sample_info.csv"
)
#s2_liver_hi <- scale.intensities(s2_liver_hi, method = "median.topn", topn = 75)
s2_liver_hi <- annotate.mz.peaks(s2_liver_hi, db)
s2_liver_hi <- peak.group.comparison(s2_liver_hi,
                                     colname = "genotype",
                                     grpA = "B6", grpB = "FVB", do.lme = F)


# ~ ~ ~ ~ ~ ~ ~ #
# Volcano  Plot #
# ~ ~ ~ ~ ~ ~ ~ #
s2.rp.p <- copy(s2_liver_rp@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$pos)
s2.rp.n <- copy(s2_liver_rp@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$neg)
s2.rp   <- rbind(s2.rp.p, s2.rp.n)[, mode := "rp"]

s2.hi.p <- copy(s2_liver_hi@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$pos)
s2.hi.n <- copy(s2_liver_hi@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$neg)
s2.hi   <- rbind(s2.hi.p, s2.hi.n)[, mode := "hi"]

s2.dt <- rbind(s2.rp, s2.hi)
s2.dt[, color := "0"]
s2.dt[wt.pval <= 0.05 & med.A > med.B, color := "1"]
s2.dt[wt.pval <= 0.05 & med.A < med.B, color := "2"]
p.vol <- ggplot(s2.dt, aes(log2(med.A/med.B), -log10(wt.pval), col = color)) + 
  geom_point(size = 0.1) + 
  geom_vline(xintercept = 0, lty = 3, col = "black") +
  geom_hline(yintercept = -log10(0.05), lty = 3, col = "black") +
  scale_color_manual(values=c("grey","blue","red")) +
  labs(y="-log10(p.value)", x = "log2(B6/FVB)") +
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p.vol


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Pathway  Enrichment #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
s2.cpd.anno   <- db@met.db[neutral.formula %in% c(s2_liver_rp@peak.annotation$neutral.formula,
                                                  s2_liver_hi@peak.annotation$neutral.formula), cpd.id]
s2.cpd.anno   <- s2.cpd.anno[!(s2.cpd.anno %in% currency.cpd)]

s2.rp.cpd     <- merge(db@met.db, s2_liver_rp@peak.annotation)
s2.rp.mz.b6   <- s2.rp[wt.pval <= 0.05 & med.A > med.B, mz.meas]
s2.rp.mz.fvb  <- s2.rp[wt.pval <= 0.05 & med.A < med.B, mz.meas]
s2.rp.cpd.b6  <- s2.rp.cpd[query.mz %in% s2.rp.mz.b6]
s2.rp.cpd.fvb <- s2.rp.cpd[query.mz %in% s2.rp.mz.fvb]

s2.hi.cpd     <- merge(db@met.db, s2_liver_hi@peak.annotation)
s2.hi.mz.b6   <- s2.hi[wt.pval <= 0.05 & med.A > med.B, mz.meas]
s2.hi.mz.fvb  <- s2.hi[wt.pval <= 0.05 & med.A < med.B, mz.meas]
s2.hi.cpd.b6  <- s2.hi.cpd[query.mz %in% s2.hi.mz.b6]
s2.hi.cpd.fvb <- s2.hi.cpd[query.mz %in% s2.hi.mz.fvb]

s2.cpd.b6     <- unique(c(s2.rp.cpd.b6$cpd.id, s2.hi.cpd.b6$cpd.id))
s2.cpd.fvb    <- unique(c(s2.rp.cpd.fvb$cpd.id, s2.hi.cpd.fvb$cpd.id))

s2.cpd.b6     <- s2.cpd.b6[!(s2.cpd.b6 %in% currency.cpd)]
s2.cpd.fvb    <- s2.cpd.fvb[!(s2.cpd.fvb %in% currency.cpd)]

isct          <- intersect(s2.cpd.b6, s2.cpd.fvb)
s2.cpd.b6     <- s2.cpd.b6[!(s2.cpd.b6 %in% isct)]
s2.cpd.fvb    <- s2.cpd.fvb[!(s2.cpd.fvb %in% isct)]

dt.rxn.pwy <- lapply(mmu.net$`Pathways of a reaction`, function(x) data.table(pwy = unlist(str_split(x, " // "))))
names(dt.rxn.pwy) <- mmu.net$Reaction
dt.rxn.pwy <- rbindlist(dt.rxn.pwy, idcol = "rxn")
dt.rxn.pwy <- dt.rxn.pwy[pwy != ""]

dt.rxn.cpd <- lapply(mmu.net$`Substrates (reactants and products) of a reaction`, function(x) data.table(cpd = unlist(str_split(x, " // "))))
names(dt.rxn.cpd) <- mmu.net$Reaction
dt.rxn.cpd <- rbindlist(dt.rxn.cpd, idcol = "rxn")
dt.rxn.cpd <- dt.rxn.cpd[cpd != ""]

dt.pwy.cpd <- merge(dt.rxn.pwy, dt.rxn.cpd)
dt.pwy.cpd <- dt.pwy.cpd[!duplicated(paste(pwy, cpd, sep = " // "))]
dt.pwy.cpd <- dt.pwy.cpd[!(cpd %in% currency.cpd)]
dt.pwy.cpd[, rxn := NULL]

# Enrichment in B6
dt.pwyenr <- copy(dt.pwy.cpd)
dt.pwyenr[, nr.cpd := .N, by = pwy]
dt.pwyenr[, cpd.sign := cpd %in% s2.cpd.b6]
dt.pwyenr[, cpd.anno := cpd %in% s2.cpd.anno]
dt.pwyenr[, nr.cpd.sign := sum(cpd.sign), by = pwy]
dt.pwyenr[, nr.cpd.anno := sum(cpd.anno), by = pwy]
dt.pwyenr <- dt.pwyenr[cpd.sign == T] # N
dt.pwyenr[, sig.cpds    := paste0(cpd, collapse = " // "), by = pwy] # N
dt.pwyenr <- dt.pwyenr[!duplicated(pwy), .(pwy, nr.cpd.sign, nr.cpd.anno, nr.cpd, sig.cpds)]
dt.pwyenr[, fisher.pval := NA_real_]
dt.pwyenr[, fisher.odds := NA_real_]

for(i in 1:nrow(dt.pwyenr)) {
  cont.mat <- matrix(c(dt.pwyenr[i, nr.cpd.anno-nr.cpd.sign], # non-sign annotated in pathway of interest
                       dt.pwyenr[i, nr.cpd.sign],
                       length(s2.cpd.anno)-length(s2.cpd.b6),
                       length(s2.cpd.b6)),
                     nrow = 2, byrow = T)
  res <- fisher.test(cont.mat, alternative = "less")
  dt.pwyenr[i, fisher.pval := res$p.value]
  dt.pwyenr[i, fisher.odds := res$estimate]
}
dt.b6up <- copy(dt.pwyenr[fisher.pval < 0.05 & nr.cpd.sign > 1])

# Enrichment in FVB
dt.pwyenr <- copy(dt.pwy.cpd)
dt.pwyenr[, nr.cpd := .N, by = pwy]
dt.pwyenr[, cpd.sign := cpd %in% s2.cpd.fvb]
dt.pwyenr[, cpd.anno := cpd %in% s2.cpd.anno]
dt.pwyenr[, nr.cpd.sign := sum(cpd.sign), by = pwy]
dt.pwyenr[, nr.cpd.anno := sum(cpd.anno), by = pwy]
dt.pwyenr <- dt.pwyenr[cpd.sign == T] # N
dt.pwyenr[, sig.cpds    := paste0(cpd, collapse = " // "), by = pwy] # N
dt.pwyenr <- dt.pwyenr[!duplicated(pwy), .(pwy, nr.cpd.sign, nr.cpd.anno, nr.cpd, sig.cpds)]
dt.pwyenr[, fisher.pval := NA_real_]
dt.pwyenr[, fisher.odds := NA_real_]

for(i in 1:nrow(dt.pwyenr)) {
  cont.mat <- matrix(c(dt.pwyenr[i, nr.cpd.anno-nr.cpd.sign], # non-sign annotated in pathway of interest
                       dt.pwyenr[i, nr.cpd.sign],
                       length(s2.cpd.anno)-length(s2.cpd.fvb),
                       length(s2.cpd.fvb)),
                     nrow = 2, byrow = T)
  res <- fisher.test(cont.mat, alternative = "less")
  dt.pwyenr[i, fisher.pval := res$p.value]
  dt.pwyenr[i, fisher.odds := res$estimate]
}
dt.fvbup <- dt.pwyenr[fisher.pval < 0.05 & nr.cpd.sign > 1]
dt.enrichment <- rbind(dt.b6up[, strain := "B6"],
                       dt.fvbup[, strain := "FVB"])

pwy.names <- fread(system.file("extdata", "metacyc_All_pathways_of_M._musculus.csv", package = "metabo.pwy"))
dt.enrichment <- merge(dt.enrichment, pwy.names, by.x = "pwy", by.y = "Pathways")

dt.enrichment <- dt.enrichment[order(strain, -fisher.pval)]
dt.enrichment$`Common-Name` <- factor(dt.enrichment$`Common-Name`, levels = dt.enrichment$`Common-Name`)
fwrite(dt.enrichment, file = "dt.enrichment.csv")

p.enr <- ggplot(dt.enrichment, aes(`Common-Name`, -log10(fisher.pval), fill = strain)) + 
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_manual(values=c("blue","red")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Pathway enrichment [-log10(p.value)]", x = "Pathway") +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())
p.enr

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Network presentation #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#

# prepare netowkr info
dt.rxn.pwy <- lapply(mmu.net$`Pathways of a reaction`, function(x) data.table(pwy = unlist(str_split(x, " // "))))
names(dt.rxn.pwy) <- mmu.net$Reaction
dt.rxn.pwy <- rbindlist(dt.rxn.pwy, idcol = "rxn")
dt.rxn.pwy <- dt.rxn.pwy[pwy != ""]

dt.rxn.cpd <- lapply(mmu.net$`Substrates (reactants and products) of a reaction`, function(x) data.table(cpd = unlist(str_split(x, " // "))))
names(dt.rxn.cpd) <- mmu.net$Reaction
dt.rxn.cpd <- rbindlist(dt.rxn.cpd, idcol = "rxn")
dt.rxn.cpd <- dt.rxn.cpd[cpd != ""]

dt.pwy.cpd <- merge(dt.rxn.pwy, dt.rxn.cpd)

# rxn edges
uniq.rxn  <- unique(dt.rxn.cpd$rxn)
rxn.edges <- lapply(uniq.rxn, function(x) {
  pcpd <- unique(dt.rxn.cpd[rxn == x, cpd])
  pcpd <- pcpd[!(pcpd %in% currency.cpd)]
  if(length(pcpd) < 2)
    return(  data.table(from      = character(0),
                        to        = character(0),
                        edge.lab  = character(0),
                        edge.type = character(0)))
  pcpd <- combn(pcpd, 2)
  data.table(from      = pcpd[1,],
             to        = pcpd[2,],
             edge.lab  = x,
             edge.type = "rxn")
})
rxn.edges <- rbindlist(rxn.edges)

# pathway edges
uniq.pwy  <- unique(dt.pwy.cpd$pwy)
pwy.edges <- lapply(uniq.pwy, function(x) {
  pcpd <- unique(dt.pwy.cpd[pwy == x, cpd])
  pcpd <- pcpd[!(pcpd %in% currency.cpd)]
  if(length(pcpd) < 2)
    return(  data.table(from      = character(0),
                        to        = character(0),
                        edge.lab  = character(0),
                        edge.type = character(0)))
  pcpd <- combn(pcpd, 2)
  data.table(from      = pcpd[1,],
             to        = pcpd[2,],
             edge.lab  = x,
             edge.type = "pwy")
})
pwy.edges <- rbindlist(pwy.edges)

pwy.net <- rbind(rxn.edges, pwy.edges)

# rxn -x- genes
bm   <- fread("genearray/mart_export.txt-1.gz") # AFFYID to
garr <- fread("genearray/Diff_gene_all_data.csv")
dt.g <- merge(garr, bm , by.x = "AFFYID", by.y = "AFFY MoGene 1 0 st v1 probe")
dt.g <- dt.g[, .(mgi.id = `MGI ID`, adj.P.Val, logFC)]
dt.g <- dt.g[order(adj.P.Val)][!duplicated(mgi.id)] # get only the most significant entry

dt.rxn.gen <- lapply(mmu.net$`Genes of a reaction`, function(x) data.table(gene = unlist(str_split(x, " // ")))) 
names(dt.rxn.gen) <- mmu.net$Reaction
dt.rxn.gen  <- rbindlist(dt.rxn.gen, idcol = "rxn")
dt.rxn.gen  <- dt.rxn.gen[gene != ""]
rel.mgi.b6 <- dt.g[adj.P.Val <= 0.05 & logFC <  0, mgi.id]
rxn.sig.b6  <- dt.rxn.gen[gene %in% rel.mgi.b6, rxn]
rel.mgi.fvb <- dt.g[adj.P.Val <= 0.05 & logFC > 0, mgi.id]
rxn.sig.fvb <- dt.rxn.gen[gene %in% rel.mgi.fvb, rxn]
rg.isct <- intersect(rxn.sig.b6, rxn.sig.fvb)
rxn.sig.b6  <- rxn.sig.b6[!(rxn.sig.b6 %in% rg.isct)]
rxn.sig.fvb <- rxn.sig.fvb[!(rxn.sig.fvb %in% rg.isct)]

# assign directions to nodes and edges
pwy.net[,edge.dir := "neutral"]
pwy.net[edge.type == "rxn" & edge.lab %in% rxn.sig.b6 , edge.dir := "B6"]
pwy.net[edge.type == "rxn" & edge.lab %in% rxn.sig.fvb, edge.dir := "FVB"]

#pwy.net[edge.type == "pwy" & edge.lab %in% dt.enrichment[strain=="B6", pwy], edge.dir := "B6"]
#pwy.net[edge.type == "pwy" & edge.lab %in% dt.enrichment[strain=="FVB", pwy], edge.dir := "FVB"]

pwy.net[,from.dir := "neutral"]
pwy.net[from %in% s2.cpd.b6, from.dir := "B6"]
pwy.net[from %in% s2.cpd.fvb, from.dir := "FVB"]

pwy.net[,to.dir := "neutral"]
pwy.net[to %in% s2.cpd.b6, to.dir := "B6"]
pwy.net[to %in% s2.cpd.fvb, to.dir := "FVB"]

# assign weights to edges
pwy.net <- pwy.net[edge.type != "pwy"]
pwy.net[, edge.weight := 0]
pwy.net[, edge.weight := edge.weight + (edge.dir != "neutral")]
pwy.net[, edge.weight := edge.weight + (from.dir != "neutral")]
pwy.net[, edge.weight := edge.weight + (to.dir != "neutral")]

pwy.net <- pwy.net[edge.weight > 1]
pwy.net <- pwy.net[!(edge.weight == 1 & edge.type == "pwy" & edge.dir != "neutral")] # if edge is supported only by pwy significance, remove it
pwy.net.nodes <- data.table(node.id  = c(pwy.net$from, pwy.net$to),
                            node.dir = c(pwy.net$from.dir, pwy.net$to.dir))
pwy.net.nodes <- pwy.net.nodes[!duplicated(node.id)]
# labels if sign
pwy.net.nodes[, label.sw := NA_character_]
pwy.net.nodes[node.dir != "neutral", label.sw := node.id]
pwy.net.nodes <- merge(pwy.net.nodes,db@met.db[,.(cpd.id, cpd.name)], by.x = "label.sw", by.y = "cpd.id"
                         ,all.x = T, sort = F)
pwy.net.nodes <- pwy.net.nodes[, .(node.id, node.dir, label.sw, cpd.name)]
pwy.net[, label.sw := NA_character_]
pwy.net[edge.dir != "neutral", label.sw := edge.lab]

library(igraph)
g <- graph.data.frame(pwy.net, vertices=pwy.net.nodes, directed=FALSE)

V(g)$size         <- 2
V(g)$label        <- NA
V(g)$frame.color  <- "white"
V(g)$color        <- ifelse(pwy.net.nodes$node.dir == "B6", "blue",
                            ifelse(pwy.net.nodes$node.dir == "FVB", "red", "grey"))
E(g)$color        <- ifelse(pwy.net$edge.dir == "B6", "blue",
                            ifelse(pwy.net$edge.dir == "FVB", "red", "grey"))
E(g)$lty          <- ifelse(pwy.net$edge.type == "pwy", 3, 1)
E(g)$width        <- ifelse(pwy.net$edge.type == "pwy", 0.5, 1)
plot(g)


V(g)$label        <- pwy.net.nodes$cpd.name
V(g)$label.family <- "sans"
V(g)$label.cex    <- 0.5
V(g)$label.color  <- "black"
#E(g)$label        <- pwy.net$label.sw
#E(g)$label.family <- "sans"

plot(g)
pdf("gfx/C_metnet.pdf",width=6,height=6)
plot(g)
dev.off()

# save also the other plots as pdf
ggsave("gfx/A_vol.pdf", p.vol, height = 3, width = 5)
ggsave("gfx/B_enr.pdf", p.enr, height = 2.5, width = 6)

# save network 
net.out <- copy(pwy.net)[,-("label.sw")]
tmp.rxngen <- merge(dt.rxn.gen, dt.g, by.x = "gene", by.y = "mgi.id")
tmp.rxngen <- tmp.rxngen[gene != ""]

net.out <- merge(net.out, tmp.rxngen, by.x = "edge.lab", by.y = "rxn", all.x = T)

fwrite(net.out, file = "network_table.csv")



# tkid <- tkplot(g)
# tkconfigure(igraph:::.tkplot.get(tkid)$canvas, "bg"="white")
# 
# l <- tkplot.getcoords(tkid)
# 
# 
# l <- layout_with_fr(g, weights = pwy.net$edge.weight)
# plot(g, l)
# 
# l <- layout_with_kk(g, weights = (pwy.net$edge.weight-0.9)/2)
# plot(g, l)
# 
# l <- layout_with_drl(g, weights = pwy.net$edge.weight)
# plot(g, l)
# 
# l <- layout_with_lgl(g)
# plot(g, l)


# 
# List for significantly abundant metabolites
#

# -> Up in FVB <-
# SET2 - RP
dt.tmp.p <- merge(s2_liver_rp@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$pos[wt.pval <= 0.05 & med.A < med.B],
                  s2_liver_rp@peak.annotation[ion.mode == "positive"],
                  by.x = "mz.meas", by.y = "query.mz")

dt.tmp.n <- merge(s2_liver_rp@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$neg[wt.pval <= 0.05 & med.A < med.B],
                  s2_liver_rp@peak.annotation[ion.mode == "negative"],
                  by.x = "mz.meas", by.y = "query.mz")

dt.met.rp <- rbind(dt.tmp.p, dt.tmp.n)
dt.met.rp <- merge(dt.met.rp, db@met.db[,.(cpd.id, cpd.name, chem.formula, cpd.charge, neutral.formula)], by = "neutral.formula")

# SET2 - HILIC
dt.tmp.p <- merge(s2_liver_hi@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$pos[wt.pval <= 0.05 & med.A < med.B],
                  s2_liver_hi@peak.annotation[ion.mode == "positive"],
                  by.x = "mz.meas", by.y = "query.mz")

dt.tmp.n <- merge(s2_liver_hi@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$neg[wt.pval <= 0.05 & med.A < med.B],
                  s2_liver_hi@peak.annotation[ion.mode == "negative"],
                  by.x = "mz.meas", by.y = "query.mz")

dt.met.hi <- rbind(dt.tmp.p, dt.tmp.n)
dt.met.hi <- merge(dt.met.hi, db@met.db[,.(cpd.id, cpd.name, chem.formula, cpd.charge, neutral.formula)], by = "neutral.formula")

dt.met.FVBup <- rbind(dt.met.rp[, LCMS.col := "RP"],
                      dt.met.hi[, LCMS.col := "HILIC"])
dt.met.FVBup <- dt.met.FVBup[order(wt.pval)]


# -> Up in B6 <-
# SET2 - RP
dt.tmp.p <- merge(s2_liver_rp@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$pos[wt.pval <= 0.05 & med.A > med.B],
                  s2_liver_rp@peak.annotation[ion.mode == "positive"],
                  by.x = "mz.meas", by.y = "query.mz")

dt.tmp.n <- merge(s2_liver_rp@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$neg[wt.pval <= 0.05 & med.A > med.B],
                  s2_liver_rp@peak.annotation[ion.mode == "negative"],
                  by.x = "mz.meas", by.y = "query.mz")

dt.met.rp <- rbind(dt.tmp.p, dt.tmp.n)
dt.met.rp <- merge(dt.met.rp, db@met.db[,.(cpd.id, cpd.name, chem.formula, cpd.charge, neutral.formula)], by = "neutral.formula")

# SET2 - HILIC
dt.tmp.p <- merge(s2_liver_hi@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$pos[wt.pval <= 0.05 & med.A > med.B],
                  s2_liver_hi@peak.annotation[ion.mode == "positive"],
                  by.x = "mz.meas", by.y = "query.mz")

dt.tmp.n <- merge(s2_liver_hi@stat.tests$`stat.pgc_genotype_B6--VS--FVB`$neg[wt.pval <= 0.05 & med.A > med.B],
                  s2_liver_hi@peak.annotation[ion.mode == "negative"],
                  by.x = "mz.meas", by.y = "query.mz")

dt.met.hi <- rbind(dt.tmp.p, dt.tmp.n)
dt.met.hi <- merge(dt.met.hi, db@met.db[,.(cpd.id, cpd.name, chem.formula, cpd.charge, neutral.formula)], by = "neutral.formula")

dt.met.B6up <- rbind(dt.met.rp[, LCMS.col := "RP"],
                      dt.met.hi[, LCMS.col := "HILIC"])
dt.met.B6up <- dt.met.B6up[order(wt.pval)]

dt.met <- rbind(dt.met.FVBup[, enr.genotype := "FVB"],
                dt.met.B6up[, enr.genotype := "B6"])

dt.met <- dt.met[, .(enr.genotype, cpd.id, cpd.name, neutral.formula, wt.pval, ion.mode.x, adduct, mz.pred, LCMS.col)]
fwrite(dt.met, "#_Table_S1.csv")


#
# significant genes in network
#

# up in FVB
tmp.Grxn <- pwy.net[edge.dir=="FVB"]
tmp.Grxn <- merge(tmp.Grxn, dt.rxn.gen, by.x = "edge.lab", by.y = "rxn")
tmp.Grxn <- tmp.Grxn[gene %in% rel.mgi.fvb]
tmp.Grxn <- tmp.Grxn[, .(involved.met = paste(from, to, sep = ",", collapse = ",")), by = c("edge.lab", "gene")]
fwrite(tmp.Grxn, "#_Table_S4_FVB.csv")

# up in B6
tmp.Grxn <- pwy.net[edge.dir=="B6"]
tmp.Grxn <- merge(tmp.Grxn, dt.rxn.gen, by.x = "edge.lab", by.y = "rxn")
tmp.Grxn <- tmp.Grxn[gene %in% rel.mgi.b6]
tmp.Grxn <- tmp.Grxn[, .(involved.met = paste(from, to, sep = ",", collapse = ",")), by = c("edge.lab", "gene")]
fwrite(tmp.Grxn, "#_Table_S4_B6.csv")

