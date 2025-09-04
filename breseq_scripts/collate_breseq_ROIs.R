setwd("E:\\Desktop_Laptop_Xfer\\Nano-Enabled_Antimicrobials_Working_Group\\R21\\Preliminary_Data_R21\\Stabryla_Unpublished")

# The parsing ----
pacman::p_load(tidyverse, rvest)

file_list <- list.files(pattern="*.html")

promote_first_row_to_colnames <- function(df) {
  colnames(df) <- df[1, ]
  df <- df[-1, ]
  return(df)
}


master_list <- list()

for(filename in file_list){
  tempdf <- read_html(filename) %>% html_table() # read html
  wrkngdf <- tempdf[-1] # Remove first iterm in list (not needed)
  
  # make list of mutations
  mutation_list <- list()
  mutation_type_list <- c() 
  for(i in 1:length(wrkngdf)) {
    mutation_type <- colnames(wrkngdf[[i]])[1]
    mutation_table <- promote_first_row_to_colnames(wrkngdf[[i]]) 
    mutation_list[[i]] <- mutation_table 
    mutation_type_list <- c(mutation_type_list, mutation_type )
  }  
  names(mutation_list) <- mutation_type_list 
  condensed_name <- str_extract(filename, "(?<=output_)[^\\-]+")
  condensed_name <- str_remove(condensed_name, "_index.html")
  master_list[[condensed_name]] <- mutation_list
}

## Bind the variants ---- 
for (name in names(master_list)) {
  # Extract the Predicted mutations dataframe
  predicted_mutations <- master_list[[name]]$`Predicted mutations`
  # Add the ID_strain column with the current name
  predicted_mutations <- predicted_mutations %>%
    mutate(ID_strain = name)
  # Update the list with the modified dataframe
  master_list[[name]]$`Predicted mutations` <- predicted_mutations
}

combined_df <- map_dfr(names(master_list), function(name) {
  master_list[[name]]$`Predicted mutations` %>%
    mutate(ID_strain = name)
})

combined_df %>%
  mutate(position = as.numeric(gsub(",", "", position))) -> 
  combined_predicted_mutations

unique(combined_predicted_mutations$ID_strain)
custom_order <- 
  c("U1", "U2", "U3", "U4", "U5", "A1", "A2", "A5", "C1", "C2", "C3", "C4", "C5")

combined_predicted_mutations$ID_strain <- 
  factor(combined_predicted_mutations$ID_strain, levels = custom_order)

# Sort the dataframe by the custom ordered factor
combined_predicted_mutations <- 
  combined_predicted_mutations[order(combined_predicted_mutations$ID_strain), ]

## Parse the mutations to get plot_variants ----

# RA = Read alignment evidence
# MC = Missing coverage evidence
# JC = New junction evidence
# UN = Unknown base evidence

# SNP = single nucleotide polymorphism
# SUB = substitution
# DEL = deletion
# INS = insertion
# MOB = mobile element
# AMP = amplification
# CON = gene conversion
# INV = inversion


ROIs <- c("32611","96910","96915","96917","639834","1080775",
          "1834467","1834959","2607200","2683945","2683967","3814231",
          "3817598","4211591")

combined_predicted_mutations <- combined_predicted_mutations %>% 
  filter(position %in% ROIs)

row <- combined_predicted_mutations[2,]
parse_row <- function(row) {
  final_colnames <- c("ID_strain", "base", "type", "delta.bases", "gene")
  
  if (str_detect(row$mutation, "→")) { # finds the SNPS
    ID_strain <- row$ID_strain
    base <- row$position
    type <- "snp"
    delta.bases <- 0
    gene <- str_remove_all(row$gene, "[\\[\\]→←]") %>% str_trim()
    row <- data.frame(ID_strain, base, type, delta.bases, gene, stringsAsFactors = FALSE)
  } else if (str_detect(row$mutation, "Δ")) { # finds the deletions
    ID_strain_v <- c(row$ID_strain, row$ID_strain)
    del.bp <- as.numeric(str_remove_all(row$mutation, "[^0-9]"))
    type_v <- c("del", "del")
    delta.bases_v <- c(del.bp, -del.bp)
    base_v <- c(row$position, row$position + del.bp)
    gene_tmp <- str_remove_all(row$gene, "[\\[\\]→←]")
    if (str_detect(gene_tmp, "/")) {
      gene_v <- c(
        str_split_1(gene_tmp, "/")[1] %>% str_trim() %>% str_replace_all(" ", ""),
        str_split_1(gene_tmp, "/")[2] %>% str_trim() %>% str_replace_all(" ", "")
      )
    } else {
      gene_v <- rep(gene_tmp, 2)
    }
    row <- data.frame(ID_strain = ID_strain_v, base = base_v, type = type_v, delta.bases = delta.bases_v, gene = gene_v, stringsAsFactors = FALSE)
  } else { # otherwise it's a structural variant
    ID_strain <- row$ID_strain
    base <- row$position
    type <- "ins"
    delta.bases <- 0
    gene_tmp <- str_remove_all(row$gene, "[\\[\\]→←]")
    gene <- str_split_1(gene_tmp, "/") %>% str_trim() %>% str_replace_all(" ", "")
    row <- data.frame(ID_strain, base, type, delta.bases, gene, stringsAsFactors = FALSE)
  }
  
  return(row)
}

parsed_rows <- combined_predicted_mutations %>%
  split(1:nrow(.)) %>%
  map_dfr(parse_row)

parsed_rows$treatment<- c(
  rep("untreated", sum(grepl("U[1-5]", parsed_rows$ID_strain)) ) ,
  rep("ampicillin", sum(grepl("A[1-5]", parsed_rows$ID_strain)) ) ,
  rep("cephalexin", sum(grepl("C[1-5]", parsed_rows$ID_strain)) ) )

parsed_rows$strain <- c(
  rep(1, sum(grepl("U1", parsed_rows$ID_strain)) ),
  rep(2, sum(grepl("U2", parsed_rows$ID_strain)) ),
  rep(3, sum(grepl("U3", parsed_rows$ID_strain)) ),
  rep(4, sum(grepl("U4", parsed_rows$ID_strain)) ),
  rep(5, sum(grepl("U5", parsed_rows$ID_strain)) ),
  
  rep(6, sum(grepl("A1", parsed_rows$ID_strain)) ),
  rep(7, sum(grepl("A2", parsed_rows$ID_strain)) ),
  rep(8, sum(grepl("A5", parsed_rows$ID_strain)) ),
  
  rep(9, sum(grepl("C1", parsed_rows$ID_strain)) ),
  rep(10, sum(grepl("C2", parsed_rows$ID_strain)) ),
  rep(11, sum(grepl("C3", parsed_rows$ID_strain)) ),
  rep(12, sum(grepl("C4", parsed_rows$ID_strain)) ),
  rep(13, sum(grepl("C5", parsed_rows$ID_strain)) ) )
  
  
parsed_rows_out <- parsed_rows %>% select(treatment, ID_strain, strain,
                                             base, type, delta.bases, gene) %>% 
  rename(mut.type = type)

parsed_rows_out <- parsed_rows_out %>% arrange(strain)

# The plotting ----
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--plotting_options"), type="character", default=NULL, 
              help="plotting options file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


raw.seq.df <- as.data.frame(parsed_rows_out)


# raw.seq.df <- raw.seq.df %>% filter(base %in% ROIs)

plotting.options <- read.csv("plotting_options.csv", 
                             header = TRUE, na.strings=c(""))


# parse the plotting options file to customize the plotting
N.bases <- as.numeric(plotting.options$Value[plotting.options$Options == "N.bases"])
section.colors <- unname(unlist(plotting.options[plotting.options$Options == "section.colors", 3:length(plotting.options)]))
section.colors <- section.colors[complete.cases(section.colors)]
bp.fill <- plotting.options$Value[plotting.options$Options == "bp.fill"]
bp.colors <- unname(unlist(plotting.options[plotting.options$Options == "bp.colors", 3:length(plotting.options)]))
bp.colors <- bp.colors[complete.cases(bp.colors)]
treatments <- unname(unlist(plotting.options[plotting.options$Options == "treatments", 3:length(plotting.options)]))
treatments <- treatments[complete.cases(treatments)]
if (treatments[1] == FALSE) {
  treatments <- unique(raw.seq.df$treatment) # default is all of the treatments. User can pass a subset of the samples if desired
}
treatment.colors <- unname(unlist(plotting.options[plotting.options$Options == "treatment.colors", 3:length(plotting.options)]))
treatment.colors <- treatment.colors[complete.cases(treatment.colors)]
if (treatment.colors[1] == FALSE) {
  library('RColorBrewer')
  treatment.colors <- brewer.pal(length(treatments),"Set3")
}
treatment.labels <- unname(unlist(plotting.options[plotting.options$Options == "treatment.labels", 3:length(plotting.options)]))
treatment.labels <- treatment.labels[complete.cases(treatment.labels)]
if (treatment.labels[1] == FALSE) {
  treatment.labels <- treatments
}
gene.color = plotting.options$Value[plotting.options$Options == "gene.color"]
gene.width = as.numeric(plotting.options$Value[plotting.options$Options == "gene.width"])
gene.section.border = plotting.options$Value[plotting.options$Options == "gene.section.border"]
mag.border = plotting.options$Value[plotting.options$Options == "mag.border"]
min.gap.length = plotting.options$Value[plotting.options$Options == "min.gap.length"]
if (min.gap.length != 'gene') {
  min.gap.length = as.numeric(min.gap.length)
}
section.buffer = as.numeric(plotting.options$Value[plotting.options$Options == "section.buffer"])
mag.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "mag.spacer"])
tab.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "tab.spacer"])
treatment.label.cex= as.numeric(plotting.options$Value[plotting.options$Options == "treatment.label.cex"])
treatment.label.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "treatment.label.spacer"])
snp.mutation.color = plotting.options$Value[plotting.options$Options == "snp.mutation.color"]
ins.mutation.color = plotting.options$Value[plotting.options$Options == "ins.mutation.color"]
del.mutation.color = plotting.options$Value[plotting.options$Options == "del.mutation.color"]
hist.spacer = as.numeric(plotting.options$Value[plotting.options$Options == "hist.spacer"])
hist.unit.height= as.numeric(plotting.options$Value[plotting.options$Options == "hist.unit.height"])
hist.axis.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "hist.axis.spacer"])
hist.ticks.length= as.numeric(plotting.options$Value[plotting.options$Options == "hist.ticks.length"])
hist.numbers.cex= as.numeric(plotting.options$Value[plotting.options$Options == "hist.numbers.cex"])
hist.numbers.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "hist.numbers.spacer"])
hist.label = plotting.options$Value[plotting.options$Options == "hist.label"]
hist.label.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "hist.label.spacer"])
hist.label.cex= as.numeric(plotting.options$Value[plotting.options$Options == "hist.label.cex"])
xmin= as.numeric(plotting.options$Value[plotting.options$Options == "xmin"])
xmax = 1+hist.unit.height+hist.spacer+0.05
ymin= as.numeric(plotting.options$Value[plotting.options$Options == "ymin"])
ymax= as.numeric(plotting.options$Value[plotting.options$Options == "ymax"])
row.spec.width= as.numeric(plotting.options$Value[plotting.options$Options == "row.spec.width"])
border.color = plotting.options$Value[plotting.options$Options == "border.color"]
treatment.border.color = plotting.options$Value[plotting.options$Options == "treatment.border.color"]
output_filename = plotting.options$Value[plotting.options$Options == "output_filename"]

#Reserve space for the plot
plot(c(xmin,xmax),c(ymin,ymax),col="white",axes=FALSE,ylab="",xlab="")

#Draw full gene
rect(0, 1-gene.width, 1, 1, col=gene.color) #xleft, ybottom, xright, ytop

#Pick the specified treatments out of the data frame
treatments <- unique(raw.seq.df$treatment)
focal.raw.seq.df <- raw.seq.df[is.element(raw.seq.df$treatment,treatments),]	

#Construct site vectors/parameters
focal.raw.seq.df.order <- focal.raw.seq.df[order(focal.raw.seq.df$base),] # orders the mutations by position
sites <- sort(unique(focal.raw.seq.df$base)) #gets all the positions for sites that will be plotted that are unique
genes <- vector(mode="character", length=length(sites)) # a vector for the gene names for each mutation
for (i in c(1:length(sites))) {
  genes[i] <- toString(focal.raw.seq.df$gene[focal.raw.seq.df$base==sites[i]][1])
}
Nsit <- length(sites) #number of variant sites that will be plotted

#Determine the number and boundaries of genetic sections to plot
sec.init <- (sites[1]-section.buffer) #grabs the first site with a mutation and subtracts the amount of bp before the variant to start plotting
if(Nsit == 1) { #if the number of variant site is one
  sec.fin <- (sites[1]+section.buffer) 
  N.sec <- 1
}
if(Nsit > 1) { #if there is more than one variant site to plot
  N.sec <- 1 # counts the number of sections to be drawn
  j <- 1
  for(i in 2:Nsit) { #loop through the number of variants starting at 2. For example, i will be 2 then 3 then 4
    if (min.gap.length != "gene") {
      if((sites[i]-sites[i-1]) > min.gap.length) { #if the current mutation is more than the space specified between variants to be plotted in the same panel. Do this loop!
        sec.init <- c(sec.init, sites[i]-section.buffer) #add to the section initiation vector the starting position for the panel
        ifelse(j==1, sec.fin <- sites[i-1] + section.buffer, sec.fin <- c(sec.fin,sites[i-1]+section.buffer)) # if on the first section initiate the sec.fin vector else append to the vector
        j <- j+1
        N.sec <- N.sec+1
      }
    } else {
      if (genes[i] != genes[i-1]) { #if the current mutation is not in the same gene as the previous variant, Do this loop!
        sec.init<-c(sec.init,sites[i]-section.buffer) #add to the section initiation vector the starting position for the panel
        ifelse(j==1,sec.fin<-sites[i-1]+section.buffer,sec.fin<-c(sec.fin,sites[i-1]+section.buffer))
        j<-j+1
        N.sec<-N.sec+1
      }
    }
  }
  sec.fin<-c(sec.fin,sites[Nsit]+section.buffer)
}	

section.colors <- rep(section.colors, ceiling(N.sec/length(section.colors)))
#Draw sections on gene
for(i in 1:N.sec) {
  rect(sec.init[i]/N.bases, 1-gene.width, sec.fin[i]/N.bases, 1, col=section.colors[i], border=gene.section.border) #xleft, ybottom, xright, ytop
}	
lines(c(0,1),c(1,1))

#Compute the lengths of each section and total length of all sections
new.run <- numeric(N.sec) # length (bp) of each section
tot <- 0 # total length (bp) of all sections
for(i in 1:N.sec) {
  new.run[i] <- (sec.fin[i]-sec.init[i])+1 #this needed a plus one
  tot <- tot + new.run[i]
}

#Put together "new" sections that are confluent on the interval [0,1], preserving relative lengths
new.sec.init <- numeric(N.sec)
new.sec.fin <- numeric(N.sec)
tot2 <- 0
for(i in 1:N.sec) {
  new.sec.init[i] <- tot2
  new.sec.fin[i] <- tot2+(new.run[i]/tot)
  tot2 <- tot2 + (new.run[i]/tot)
}

#"Magnify" the sections
for(i in 1:N.sec) {
  polygon(c(sec.init[i]/N.bases, new.sec.init[i], new.sec.fin[i], sec.fin[i]/N.bases),
          c(1-gene.width, 1-gene.width-mag.spacer, 1-gene.width-mag.spacer, 1-gene.width),
          col = section.colors[i], border = mag.border) #the x vector c(top left, bottom left, bottom right, top right), the y vector (top left, bottom left, bottom right, top right)
}
lines(c(0,1),c(1-gene.width,1-gene.width))


#Function to organize strains in a list first by number of mutations
#then by location of the first mutation (this makes the mutation table
#a bit prettier)
OrderStrains<-function(X.df) {
  Xstrains<-sort(unique(X.df$strain))
  len<-length(Xstrains)
  x<-rep(0,len)
  y<-rep(0,len)
  z<-rep(0,len)
  SO.df<-data.frame(strain=x,base=y,N.mut=z)
  for(i in 1:len) {
    SO.df$strain[i]<-X.df[X.df$strain==Xstrains[i],]$strain[1]
    SO.df$base[i]<-min(X.df[X.df$strain==Xstrains[i],]$base)
    nm<-length(X.df[X.df$strain==Xstrains[i],]$base)
    count=0
    for(j in X.df[X.df$strain==Xstrains[i],]$mut.type) {
      if(j == 'del') {
        count <- count + 1
      }
    }
    SO.df$N.mut[i] <- (nm - (count/2))
  }
  SO.df<-SO.df[order(SO.df$N.mut,SO.df$base),]
  SO.df
}

#Organize strains in all treatments and collect information on the number
#of mutations in each.	
for(i in 1:length(treatments)) { # loop through the treatments
  T.df<-OrderStrains(raw.seq.df[raw.seq.df$treatment==treatments[i],]) # pass all of the mutations for a particular treatment to the OrderStrains function
  ifelse(i==1, strains<-T.df$strain, strains<-c(strains,T.df$strain))
  ifelse(i==1, mutHist<-T.df$N.mut, mutHist<-c(mutHist,T.df$N.mut))
}

Nstr<-length(strains) #the number of strains to be plotted

# calculate the row width to fill the y-axis
if (row.spec.width == 0) {
  row.width<-(1-gene.width-mag.spacer-tab.spacer)/(Nstr+1)  #the width of each row (strain)
} else {
  row.width <- row.spec.width
}

#Draw confluent sections at row height for the magnified bar
for(i in 1:N.sec) {
  rect(new.sec.init[i], 1-gene.width-mag.spacer-row.width, new.sec.fin[i], 1-gene.width-mag.spacer, col=section.colors[i]) #xleft, ybottom, xright, ytop
}


#Put together a list of all the bases in the table
focal.bases <- c((sec.init[1]):(sec.fin[1]))
if(N.sec>1) {
  for(i in 2:N.sec) {
    nex <- c((sec.init[i]):(sec.fin[i]))
    focal.bases <- c(focal.bases, nex)
  }
}

#Draw the proper background colors for each strain
current.treatment <- raw.seq.df[raw.seq.df$strain==strains[1],]$treatment[1]
for(j in 1:length(treatments)) {
  if(current.treatment==treatments[j]) {
    current.color <- treatment.colors[j]
  }
}
row.of.last.color <- 0
for(i in 1:Nstr) {
  if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1] != current.treatment) {	
    rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), 1, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)), col=current.color, border=FALSE) #xleft, ybottom, xright, ytop
    for(j in 1:length(treatments)) {
      if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]==treatments[j]) {
        current.color<-treatment.colors[j]
      }
    }
    row.of.last.color<-(i-1)
    current.treatment<-raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]
  }
  if(i==Nstr) {
    rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), 1, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)), col=current.color, border=FALSE)
  }	
}

#draw in the tab.spacer bp rectangles
if (bp.fill != FALSE){
  lfb<-length(focal.bases)
  bp.colors <- rep(bp.colors, ceiling(lfb/length(bp.colors)))
  bp_size <- 1/lfb
  for (pos in 1:lfb) {
    rect((pos-1)*bp_size, 1-gene.width-mag.spacer-row.width-tab.spacer, pos*bp_size, 1-gene.width-mag.spacer-row.width, col= bp.colors[pos], border=NA)
  }
  rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer, 1, 1-gene.width-mag.spacer-row.width, col= NA, border='black')
}

#Draw the mutations for each strain
# strain 1, stain 2, etc
for(i in 1:Nstr) {
  mutations<-raw.seq.df[raw.seq.df$strain==strains[i],]$base
  type<-raw.seq.df[raw.seq.df$strain==strains[i],]$mut.type
  lfb<-length(focal.bases)
  for(j in 1:length(mutations)) {
    pos<-which(focal.bases==mutations[j])
    if(type[j]=='snp') {
      mut.col <- snp.mutation.color
    }  
    if(type[j]=='ins') {
      mut.col <- ins.mutation.color
    }
    if(type[j]=='del') {
      mut.col <- del.mutation.color
    }
    rect((pos-1)/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), pos/lfb,1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), col=mut.col, border=NA)
    if(type[j]=='del') {
      db <- raw.seq.df[raw.seq.df$strain==strains[i] & raw.seq.df$base==mutations[j],]$delta.bases
      print(db)
      if(any(db>0)) {
        end.pos<-which(focal.bases==(mutations[j]+db))
        rect((pos-1)/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), end.pos/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), col=mut.col, border=NA)
        #        for(p in pos:end.pos) {
        #          rect((p-1)/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), p/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), col=mut.col, border=NA)
        #        }
        if(length(raw.seq.df[raw.seq.df$strain==strains[i] & 
                             raw.seq.df$base==(mutations[j]+db),]$delta.bases) != 1) {
          print(paste("There is an error; strain ", strains[i], 
                      " should have a deletion ending at ", mutations[j]+db, sep=""))
        }
      }
    }
  }
}

#Draw boundaries around all magnified sections in the strain rows
for(i in 1:N.sec) {
  rect(new.sec.init[i], 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*Nstr), new.sec.fin[i], 1-gene.width-mag.spacer-row.width-tab.spacer, border = border.color)
}

#Draw lines between the strains in the mutation table
row.of.last.color<-0
for(i in 1:Nstr) {
  if(i!=Nstr) {	
    rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)),1, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), border = border.color)
    for(j in 1:length(treatments)) {
      if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]==treatments[j]) {
        current.color<-treatment.colors[j]
      }
    }
    row.of.last.color<-(i-1)
    current.treatment<-raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]
  }
  if(i==Nstr) {
    rect(0,1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)),1,
         1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)),border=border.color)
  }	
}

#Draw lines around the treatments 
current.treatment <- raw.seq.df[raw.seq.df$strain==strains[1],]$treatment[1]
for(j in 1:length(treatments)) {
  if(current.treatment==treatments[j]) {
    current.color<-treatment.colors[j]
  }
}
row.of.last.color<-0

for(i in 1:Nstr) {
  if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1] != current.treatment) {	
    rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), 1, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)), border=treatment.border.color)
    for(j in 1:length(treatments)) {
      if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]==treatments[j]) {
        current.color<-treatment.colors[j]
      }
    }
    row.of.last.color<-(i-1)
    current.treatment<-raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]
  }
  if(i==Nstr) {
    rect(0,1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)),1,
         1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)), border=treatment.border.color)
  }	
}


#Draw histogram to the right of the table
hist.unit.height <- hist.unit.height/max(mutHist)
for(i in 1:Nstr) {
  for(j in 1:length(treatments)) {
    if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]==treatments[j]) {
      color<-treatment.colors[j]
    }
  }	
  rect(1+hist.spacer, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), 1+hist.spacer+mutHist[i]*hist.unit.height, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)),
       col=color)
}


#Plot the treatment labels to the left of the table
Tlen <- numeric(length(treatments))
Totlen <- 0
for(i in 1:length(treatments)) {
  Tlen[i] <- length(unique(raw.seq.df[raw.seq.df$treatment==treatments[i],]$strain))
  Totlen <- (Totlen+Tlen[i])
}
Ty<-numeric(length(treatments))
RunTot<-0
if (row.spec.width == 0) {
  for(i in 1:length(treatments)) {
    Ty[i]<-(1-gene.width-mag.spacer-row.width-tab.spacer)*(1-RunTot-0.5*(Tlen[i]/Totlen))
    text(-treatment.label.spacer, Ty[i], treatment.labels[i], col="black",cex=treatment.label.cex, pos = 2)
    RunTot<-(RunTot+(Tlen[i]/Totlen))
  }
} else {
  for(i in 1:length(treatments)) {
    mutant.space <- row.spec.width * Nstr
    Ty[i]<-(1-gene.width-mag.spacer-row.spec.width-tab.spacer)-(mutant.space-(mutant.space*(1-RunTot-0.5*(Tlen[i]/Totlen))))
    text(-treatment.label.spacer, Ty[i], treatment.labels[i],  
         col="black",cex=treatment.label.cex, pos = 2)
    RunTot<-(RunTot+(Tlen[i]/Totlen))
  }
}

#Plot the histogram axis, numbers, and label
if (row.spec.width == 0) {
  arrows(1+hist.spacer,-hist.axis.spacer,1+hist.spacer+(max(mutHist)+1)*hist.unit.height,-hist.axis.spacer,length=2*hist.unit.height)
  for(i in 1:max(mutHist)) {
    lines(c(1+hist.spacer+(i)*hist.unit.height,1+hist.spacer+(i)*hist.unit.height),
          c(-hist.axis.spacer,-hist.axis.spacer-hist.ticks.length))
    text(1+hist.spacer+(i)*hist.unit.height,-hist.numbers.spacer,i,cex=hist.numbers.cex,pos=1)
  }
  text(1+hist.spacer+(max(mutHist)+1)*.5*hist.unit.height,-hist.label.spacer,hist.label,cex=hist.label.cex,pos=1)
} else {
  Bottom.row <- 1-gene.width-mag.spacer-row.spec.width-tab.spacer-mutant.space
  arrows(1+hist.spacer,Bottom.row-hist.axis.spacer,1+hist.spacer+(max(mutHist)+1)*hist.unit.height,Bottom.row-hist.axis.spacer,length=2*hist.unit.height)
  for(i in 1:max(mutHist)) {
    lines(c(1+hist.spacer+(i)*hist.unit.height,1+hist.spacer+(i)*hist.unit.height),
          c(Bottom.row-hist.axis.spacer,Bottom.row-hist.axis.spacer-hist.ticks.length))
    text(1+hist.spacer+(i)*hist.unit.height,Bottom.row-hist.numbers.spacer,i,cex=hist.numbers.cex,pos=1)
  }
  text(1+hist.spacer+(max(mutHist)+1)*.5*hist.unit.height,Bottom.row-hist.label.spacer,hist.label,cex=hist.label.cex,pos=1)
}

text((2607200/4793125)+0.03, 1.02, "frdD", cex=hist.label.cex, pos=3, srt = 45)
text((1834467/4793125)+0.03, 1.02, "acrB", cex=hist.label.cex, pos=3, srt = 45)
text((639834/4793125)+0.03 , 1.02, "marR", cex=hist.label.cex, pos=3, srt = 45)
text((96910/4793125)+0.03 , 1.02, "xerC", cex=hist.label.cex, pos=3, srt = 45)
text((3814231/4793125)+0.03 , 1.02, "iolU", cex=hist.label.cex, pos=3, srt = 45)
text((4211591/4793125)+0.03 , 1.02, "nlpD", cex=hist.label.cex, pos=3, srt = 45)
text((1080775/4793125)+0.03 , 1.02, "putP", cex=hist.label.cex, pos=3, srt = 45)




#Add size labels 
text(0, 1, paste(0, "Mb"), cex=hist.label.cex, pos=3)
text(1, 1, paste(round(N.bases/1000000, digits = 1), "Mb"),cex=hist.label.cex, pos=3)

text((2607200/4793125)+0.03, 1.02, "frdD", cex=hist.label.cex, pos=3, srt = 45)
text((1834467/4793125)+0.03, 1.02, "acrB", cex=hist.label.cex, pos=3, srt = 45)
text((639834/4793125)+0.03 , 1.02, "marR", cex=hist.label.cex, pos=3, srt = 45)
#text(96917/4793125, 1, "xerC", cex=hist.label.cex, pos=3)
