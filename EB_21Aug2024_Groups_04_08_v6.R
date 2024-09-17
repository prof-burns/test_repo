# Giving this a try in git

#### Reading in the data and checking it over ####

# Set the working directory to be where all the files are:

getwd()
setwd("~/Desktop/EB_21_Aug_2024_Endpoints/")
list.files(getwd())

# When you use setwd (short for set working directory) to the folder where all the input and output files should go
# R then knows to look there first and you can skip adding in the full filepaths all over the place.

# Use VROOM to load the ASV table
# install.packages("vroom") # Only need to do this once to install VROOM
library(vroom)
Sys.setenv("VROOM_CONNECTION_SIZE" = 400000000)
seqtab.nochim <- vroom("otu_allsamples.csv")

# It worked, but we need to change it's class to a dataframe so r knows how to interact with it
seqtab.nochim <- as.data.frame(seqtab.nochim)

# Set the row names from the first column
rownames(seqtab.nochim) <- seqtab.nochim[[1]]

# Remove the first column, since it's already there as row names 
seqtab.nochim <- seqtab.nochim[, -1]

# Now physically click on seqtab.nochim to make sure it looks the way you think it should wrt columns, rows, etc.

# Next, simply check how many rows and columns there are (with rows being samples and columns being ASVs)
dim(seqtab.nochim)

# Now switch the class back to a matrix to allow phyloseq to work with it.
seqtab.nochim <- as.matrix(seqtab.nochim)

# Read in taxonomy
file_path <- "taxa_allsamples.csv"
taxa <- read.table(file_path, header = TRUE, row.names = 1, sep = ",")

# Remember to click on the taxa in the upper right panel to open it real quick to double check that it is what we think it is!

# Check the dimensions!
dim(taxa)

# Now switch the class to a matrix to allow phyloseq to work with it.
taxa <- as.matrix(taxa)

# Read in metadata
file_path <- "All_Samples_Metadata_forEB.txt"
samdf <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t")

# Click it to check, then inspect the dimensions
dim(samdf)

# Now that all three things are loaded, are the dimensions all consistent?
# Are there the same number of ASVs in the taxa and ASV tables?
# Are there the same number of samples in the metadata table as there are in the ASV table?

# As we're loading this into phyloseq, we should update the samples names.
# Right now, they're things like 1985set2_138_AP02Week0-1.fastq.gz, but we can do better than that.
# There's a column in the metadata table called SampleID that contains a suitable name.

# We start by sorting the metadata and ASV tables (otherwise there's a chance that the names will get shuffled!)
samdf_sorted <- samdf[order(row.names(samdf)),]
seqtab.nochim_sorted <- seqtab.nochim[order(row.names(seqtab.nochim)),]

# Double check visually that this worked then update the originals
samdf <- samdf_sorted
seqtab.nochim <- seqtab.nochim_sorted

# Update the rownames in both!
row.names(samdf) <- samdf$SampleID
row.names(seqtab.nochim) <- samdf$SampleID

# Visually doublecheck that this worked!

#### Making a phyloseq object ####
library(phyloseq)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
print(ps)

#### Filtering out uninformative taxa from ps ####

# Remove all the taxa that are still garbage or uninformative
ps1 = subset_taxa(ps, Kingdom=="Bacteria")
print(ps1)

# Check the phyla as well
unique_phyla <- as.data.frame(unique(tax_table(ps1)[, "Phylum"]))

# Remove all the taxa that are not called at the level of the phylum
ps2 = subset_taxa(ps1, Phylum!="NA")
unique_phyla1 <- as.data.frame(unique(tax_table(ps2)[, "Phylum"]))
print(ps2)

ps <- ps2
ps

# How many taxa were removed? 

#### To collapse and decontam or to decontam and then collapse? ####

# Question - Is this the right spot to do this or should we run decontam first?
# Try both! What do the results look like when you do it both ways? Which should we do? 
# What makes sense and what are the jjustifications for one over the other?

# ALSO
# Unlike collapsing taxa, we need to pay close attention to which samples we use with decontam
# The dataset we're working with is 622 samples from at least 5 different projects, including the EB one.
# A "control" sample for one project is not neccessarily a good control for another project's samples. 
# Let's start by reducing the dataset to only look at the EB samples.

# For simplicity at this point, let's only work with the EB project samples.
# If you look in samdf, there's a column that let's us easily do this!

ps_EB <- prune_samples(sample_data(ps)$PROJECT=="EB", ps)

# Check that it worked
sample_data(ps_EB)

ps <- ps_EB

#### Collapsing ASVs by taxonomy ####

# Load necessary libraries
library(dplyr)
library(tibble)

# Assume your phyloseq object is named 'physeq'
physeq <- ps

# Extract OTU table, taxonomy table, and sample data
otu_table <- as.data.frame(otu_table(physeq))
tax_table <- as.data.frame(tax_table(physeq))
sample_data <- as.data.frame(sample_data(physeq))

# Make a taxa string that I can use for collapsing
tax_table$Concatenated <- apply(tax_table, 1, function(row) paste(row, collapse = ";__"))

# Transpose the ASV table so I can add the concatenated column
otu_table_t <- as.data.frame(t(otu_table))

# Merge the Concatenated column into otu_table_t
otu_table_t$Concatenated <- tax_table$Concatenated[match(rownames(otu_table_t), rownames(tax_table))]

# Reorder the columns to make it easier to see the new concatenated one
otu_table_t <- otu_table_t %>%
  dplyr::select(Concatenated, everything())

# Group by the Concatenated column and sum the values of each column
otu_table_summed <- otu_table_t %>%
  group_by(Concatenated) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

# Convert the otu_table_summed from a tibble as the next steps balk at this
otu_table_summed <- as.data.frame(otu_table_summed)

# Add row names as ASV1, ASV2, ASV3, etc.
rownames(otu_table_summed) <- paste0("ASV", seq_len(nrow(otu_table_summed)))

# Make a new taxa table
taxa_table_summed <- data.frame(
  ASV = rownames(otu_table_summed),
  Concatenated = otu_table_summed$Concatenated
)
rownames(taxa_table_summed) <- taxa_table_summed$ASV
taxa_table_summed <- taxa_table_summed[-1]

# Rebuild the taxa table
# Load the tidyverse package
library(tidyverse)

# Use the separate function to split the Concatenated column
taxa_table_summed <- taxa_table_summed %>%
  separate(Concatenated, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";__")

# Now put all this stuff back together

# Re-transpose the asv table
otu_table_summed <- as.data.frame(otu_table_summed)
rownames(otu_table_summed) <- paste0("ASV", seq_len(nrow(otu_table_summed)))
otu_table_summed_t <- t(otu_table_summed[-1])
otu_table_summed_t <- as.data.frame(otu_table_summed_t)
samdf_summed <- as.data.frame(sample_data)
taxa_summed <- as.matrix(taxa_table_summed)

ps_collapsed <- phyloseq(otu_table(otu_table_summed_t, taxa_are_rows=FALSE), 
                         sample_data(samdf_summed), 
                         tax_table(taxa_summed))

print(ps_collapsed)

# How do we know this worked properly? One way is to compare the total reads/sample
# Before and after collapsing (they shouldn't change!)

# Calculate total reads per sample
sample_sums_pre <- as.data.frame(sample_sums(ps))
sample_sums_pre

sample_sums_post <- as.data.frame(sample_sums(ps_collapsed))
sample_sums_post

# Merge these two for simplicity
merged_sums_df <- sample_sums_pre %>%
  rename(pre = 1) %>%
  bind_cols(sample_sums_post %>% rename(post = 1))

# Write them out if desired
write.table(merged_sums_df, file="samplesums.txt",quote=FALSE)

#### Decontam ####

# This assumes that you ran the collapsing first. Update the inputs as needed to run this with the uncollapsed version.

# It never hurts to make a back-up variable of the data you're working with in case you need to circle back
ps_BackUp <- ps
ps_collapsed_BackUp <- ps_collapsed

# Assuming that your phyloseq object is named "ps"

ps <- ps_collapsed

# Now to run decontam
# BiocManager::install("decontam")
library(decontam); packageVersion("decontam")

head(sample_data(ps))

# Inspect the library sizes (just for fun)
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point(size = 0.05)

# Now to use the prevalence method to identify potential contaminants

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

# How many of the ASVs are contaminants?
decontamSet_numbers <- table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)
decontamSet_numbers

# How many of the taxa/ASVs are real and how many are contaminants? 
# When you run this, you'll need to keep close track of this for later!

# Now to re-run this with a more stringent approach:
contamdf.prev01 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
which(contamdf.prev01$contaminant)

# Does the threshold matter? If not, go with the more stringent one!

# Look at these results as phyloseq opbjects
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point(size = 0.05) +
  xlab("Prevalence (Controls)") + ylab("Prevalence (Samples)")

# Make a vector to store the contaminant ASVs
contaminant_asvs <- row.names(contamdf.prev01[contamdf.prev01$contaminant == TRUE, ])

# Now to remove the identified contaminants:
badTaxa = contaminant_asvs
allTaxa = taxa_names(ps)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ex1 = prune_taxa(allTaxa, ps)
taxa_names(ex1)
# new phyloseq object with just the taxa you kept.
ex1

# Now to filter out the rest of the unassigned and the archaea (we already did this earlier, but it never hurts)
ex2 = subset_taxa(ex1, Kingdom=="Bacteria")
ex2
# Now to do some filtering and pruning
ps <- ex2

# Prune out the Control samples
ps_no_env = subset_samples(ps, Sample_or_Control != "Control")

ps <- ps_no_env
ps

# How many samples and taxa (or ASVs) do we have left and does that make sense?

#### Merging replicate samples by timepoint ####

# Extract the sample data
sample_data <- sample_data(ps)

# Create a new column in the sample data that combines patient ID and week
sample_data$PatientWeek <- paste(sample_data$Name, sample_data$Week, sep = "_")

# Update the phyloseq object with the new sample data
sample_data(ps) <- sample_data

# Merge samples based on the PatientWeek column
ps_combined <- merge_samples(ps, "PatientWeek")

# take a look to see what the results say
sample_data_combined <- sample_data(ps_combined)

# Check the counts
otu_data_combined <- as.data.frame(otu_table(ps_combined))
colSums(otu_data_combined)
rowSums(otu_data_combined)
min(rowSums(otu_data_combined))
sum(otu_data_combined)

otu_data_ps <- as.data.frame(otu_table(ps))
colSums(otu_data_ps)
rowSums(otu_data_ps)
min(rowSums(otu_data_ps))
sum(otu_data_ps)

# Seems like it works

#### Group samples for Week00-04 and Week 00-08 ####

library(dplyr)

# Extract sample data
sample_data <- sample_data(ps_combined)

# Convert to data frame for easier manipulation
sample_data_df <- as(sample_data, "data.frame")
sample_data_df$SampleID <- rownames(sample_data_df)

# Repopulate the Name and Week columns that were lost in the merging
sample_data_df <- sample_data_df %>%
  separate(SampleID, into = c("Name", "Week"), sep = "_")

sample_data_df$SampleID <- rownames(sample_data_df)

# While I'm at it, I should add a numeric week column.
# Add a numeric Week column
sample_data_df <- sample_data_df %>%
  mutate(NumericWeek = case_when(
    Week == "Week00" ~ 0,
    Week == "Week04" ~ 4,
    Week == "Week08" ~ 8,
    Week == "Week12" ~ 12,
    TRUE ~ NA_real_ # Handle unexpected Week values
  ))

# Identify patients with all Early time points (Week00 and Week04)
early_patients <- sample_data_df %>%
  group_by(Name) %>%
  filter(all(c("Week00", "Week04") %in% Week)) %>%
  pull(Name) %>%
  unique()

# Filter the sample data to include only complete patients
early_sample_data <- sample_data_df %>%
  filter(Name %in% early_patients)

# Filter the phyloseq object
physeq_early <- prune_samples(early_sample_data$SampleID, ps_combined)

# Update the sample data in the phyloseq object with the new data
sample_data(physeq_early) <- sample_data(early_sample_data)

# Optional: Verify the filtering
sample_data(physeq_early)$Name
sample_data(physeq_early)$NumericWeek

# Remove all weeks except 00 and 04
physeq_early_filtered <- subset_samples(physeq_early, NumericWeek %in% c(0, 4))

# Optional: Verify the filtering
sample_data(physeq_early_filtered)$Name
sample_data(physeq_early_filtered)$NumericWeek

physeq_early <- physeq_early_filtered

# Identify patients with all Late time points (Week00 and Week08)
late_patients <- sample_data_df %>%
  group_by(Name) %>%
  filter(all(c("Week00", "Week08") %in% Week)) %>%
  pull(Name) %>%
  unique()

# Filter the sample data to include only complete patients
late_sample_data <- sample_data_df %>%
  filter(Name %in% late_patients)

# Filter the phyloseq object
physeq_late <- prune_samples(late_sample_data$SampleID, ps_combined)

# Update the sample data in the phyloseq object with the new data
sample_data(physeq_late) <- sample_data(late_sample_data)

# Optional: Verify the filtering
sample_data(physeq_late)$Name
sample_data(physeq_late)$NumericWeek

# Remove all weeks except 00 and 04
physeq_late_filtered <- subset_samples(physeq_late, NumericWeek %in% c(0, 8))

# Optional: Verify the filtering
sample_data(physeq_late_filtered)$Name
sample_data(physeq_late_filtered)$NumericWeek

physeq_late <- physeq_late_filtered


#### Updating the metadata ####

# As per the discussion in July, we're ditching samples from patients 5 and 6
# This is because the read counts for these samples are all over the place and will
# probably mess things up.

# Also, sample AP11 for weeks 8 and 12 are compromise by Abx use
# As are sample AP14 for weeks 4, 8, and 12

# Filtering early set
Patients_to_Remove_early <- c("AP05","AP06", "AP14")
ps_early_curated <- prune_samples(!(sample_data(physeq_early)$Name %in% Patients_to_Remove_early), physeq_early)
ps_early_curated
sample_data(ps_early_curated)$Name

physeq_complete_early <- ps_early_curated

# Filtering late set
Patients_to_Remove_late <- c("AP05","AP06", "AP14", "AP11")
ps_late_curated <- prune_samples(!(sample_data(physeq_late)$Name %in% Patients_to_Remove_late), physeq_late)
ps_late_curated
sample_data(ps_late_curated)$Name

physeq_complete_late <- ps_late_curated


# We need to add in some additional bits here in order to visualize what's 
# going on with the different groups.
# Let's just start with the samples for which the staph changed or didn't change.

# The samples that (visually) showed large shifts in SA were:
# AP02, AP12, AP12, AP14, AP16
# This sucks, since AP14 was removed for Apx use (but now it makes sense)

# Let's add something that stratifies them

# Extract sample data (metadata) to edit things with new info
sample_data <- sample_data(physeq_complete)
# Convert to data frame for easier manipulation
sample_data_df <- as(sample_data, "data.frame")
# Define file path for metadata
meta_file_path <- "~/Desktop/phyloseq_metadata.csv"
# Write the metadata to CSV
write.csv(sample_data_df, file = meta_file_path, row.names = TRUE)
# Once updated, re-read it back in and add it to the ps object
sample_data_df2 <- read.csv(meta_file_path, header=TRUE, row.names = 1)

# Add in the new columns
sample_data(physeq_complete)$StaphChange <- sample_data_df2$StaphChange
sample_data(physeq_complete)$Subtype <- sample_data_df2$Subtype
sample_data(physeq_complete)


#### Filtering the ASVs/taxa by prevalence and abundance and also filtering the samples by total read count ####

# Filtering samples by read counts

### This step is redundant as all the samples have more than 1k...
# Remove samples with less than n from phyloseq object
ps_early_morethan1k <- prune_samples(sample_sums(physeq_complete_early) >= 1000, physeq_complete_early)
ps_early_morethan1k
ps_late_morethan1k <- prune_samples(sample_sums(physeq_complete_late) >= 1000, physeq_complete_late)
ps_late_morethan1k

# Remove ASVs with less than x prevalence from phyloseq object
# Taken in part from here: https://f1000research.com/articles/5-1492

# EARLY SAMPLES
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps_early_morethan1k),
               MARGIN = ifelse(taxa_are_rows(ps_early_morethan1k), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps_early_morethan1k),
                    tax_table(ps_early_morethan1k))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps_early_morethan1k),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#  Define prevalence threshold as 9% of total samples
prevalenceThreshold = 0.09 * nsamples(ps_early_morethan1k)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps4 = prune_taxa(keepTaxa, ps_early_morethan1k)
ps_early_morethan1k_filtered <- ps4
ps_early_morethan1k_filtered

# LATE SAMPLES
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps_late_morethan1k),
               MARGIN = ifelse(taxa_are_rows(ps_late_morethan1k), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps_late_morethan1k),
                    tax_table(ps_late_morethan1k))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps_late_morethan1k),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#  Define prevalence threshold as 9% of total samples
prevalenceThreshold = 0.11 * nsamples(ps_late_morethan1k)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps4 = prune_taxa(keepTaxa, ps_late_morethan1k)
ps_late_morethan1k_filtered <- ps4
ps_late_morethan1k_filtered



# Remove low abundance taxa - EARLY
physeqr <- transform_sample_counts(ps_early_morethan1k_filtered, function(x) x / sum(x))
as.data.frame(otu_table(physeqr))[1,1]
# Filter out taxa from raw count ps object with n overall abundance
# This sets a filter at the sum of a taxon across all samples being greater than 0.05% of all OTUs.
physeqrF <- filter_taxa(physeqr, function(x) mean(x) < .00005,TRUE) # 0.00005 = 0.0005%
rmtaxa <- taxa_names(physeqrF)
alltaxa <- taxa_names(ps_early_morethan1k_filtered)
myTaxa <- alltaxa[!alltaxa %in% rmtaxa]
physeqaF_early <- prune_taxa(myTaxa,ps_early_morethan1k_filtered)
physeqaF_early

otu_table(physeqaF_early)
test_taxa <- as.data.frame(tax_table(physeqaF_early))

physeqaF_early_complete <- physeqaF_early

lookie <- as.data.frame(tax_table(physeqaF_early_complete))

# Remove low abundance taxa - LATE
physeqr <- transform_sample_counts(ps_late_morethan1k_filtered, function(x) x / sum(x))
as.data.frame(otu_table(physeqr))[1,1]
# Filter out taxa from raw count ps object with n overall abundance
# This sets a filter at the sum of a taxon across all samples being greater than 0.05% of all OTUs.
physeqrF <- filter_taxa(physeqr, function(x) mean(x) < .00005,TRUE) # 0.00005 = 0.0005%
rmtaxa <- taxa_names(physeqrF)
alltaxa <- taxa_names(ps_late_morethan1k_filtered)
myTaxa <- alltaxa[!alltaxa %in% rmtaxa]
physeqaF_late <- prune_taxa(myTaxa,ps_late_morethan1k_filtered)
physeqaF_late

otu_table(physeqaF_late)
test_taxa <- as.data.frame(tax_table(physeqaF_late))

physeqaF_late_complete <- physeqaF_late


#### Alpha Diversity - Some simple stats over the timepoints ####

# First, check out the alpha rarefaction plots

library(phyloseq)
library(ggplot2)
library(plyr)

# EARLY ALPHA RAREFACTION

# Set a range of rarefaction depths
rarefaction_depths <- seq(1, min(sample_sums(physeqaF_early_complete)), by = 5000)
min(sample_sums(physeqaF_early_complete))

# Generate the alpha diversity results with an explicit SampleID column
rarefaction_results <- ldply(rarefaction_depths, function(depth) {
  rarefied_physeq <- rarefy_even_depth(physeqaF_early_complete, sample.size = depth, verbose = FALSE)
  alpha_diversity <- estimate_richness(rarefied_physeq, measures = c("Observed"))
  alpha_diversity$Depth <- depth
  alpha_diversity$SampleID <- rownames(alpha_diversity)
  alpha_diversity
})

# Create the alpha rarefaction plot with color coding by SampleID and a legend
alpha_rarefaction_plot <- ggplot(rarefaction_results, aes(x = Depth, y = Observed, group = SampleID, color = SampleID)) +
  geom_line(alpha = 0.8) +
  theme_minimal() +
  labs(x = "Sequencing Depth", y = "Observed Diversity", title = "Alpha Rarefaction Curve", color = "Sample ID") +
  theme(legend.position = "right") # Position the legend on the right side

# Display the plot
print(alpha_rarefaction_plot)

# Generate the alpha diversity results with an explicit SampleID column
rarefaction_results <- ldply(rarefaction_depths, function(depth) {
  rarefied_physeq <- rarefy_even_depth(physeqaF_early_complete, sample.size = depth, verbose = FALSE)
  alpha_diversity <- estimate_richness(rarefied_physeq, measures = c("Shannon"))
  alpha_diversity$Depth <- depth
  alpha_diversity$SampleID <- rownames(alpha_diversity)
  alpha_diversity
})

max(rarefaction_results$Shannon)

# Create the alpha rarefaction plot with color coding by SampleID and a legend
alpha_rarefaction_plot <- ggplot(rarefaction_results, aes(x = Depth, y = Shannon, group = SampleID, color = SampleID)) +
  geom_line(alpha = 0.8) +
  theme_minimal() +
  labs(x = "Sequencing Depth", y = "Shannon Diversity Index", title = "Alpha Rarefaction Curve", color = "Sample ID") +
  theme(legend.position = "right") # Position the legend on the right side

# Display the plot
print(alpha_rarefaction_plot)


# LATE ALPHA RAREFACTION

# Set a range of rarefaction depths
rarefaction_depths <- seq(1, min(sample_sums(physeqaF_late_complete)), by = 5000)

# Generate the alpha diversity results with an explicit SampleID column
rarefaction_results <- ldply(rarefaction_depths, function(depth) {
  rarefied_physeq <- rarefy_even_depth(physeqaF_late_complete, sample.size = depth, verbose = FALSE)
  alpha_diversity <- estimate_richness(rarefied_physeq, measures = c("Observed"))
  alpha_diversity$Depth <- depth
  alpha_diversity$SampleID <- rownames(alpha_diversity)
  alpha_diversity
})

# Create the alpha rarefaction plot with color coding by SampleID and a legend
alpha_rarefaction_plot <- ggplot(rarefaction_results, aes(x = Depth, y = Observed, group = SampleID, color = SampleID)) +
  geom_line(alpha = 0.8) +
  theme_minimal() +
  labs(x = "Sequencing Depth", y = "Observed Diversity", title = "Alpha Rarefaction Curve", color = "Sample ID") +
  theme(legend.position = "right") # Position the legend on the right side

# Display the plot
print(alpha_rarefaction_plot)

# Generate the alpha diversity results with an explicit SampleID column
rarefaction_results <- ldply(rarefaction_depths, function(depth) {
  rarefied_physeq <- rarefy_even_depth(physeqaF_late_complete, sample.size = depth, verbose = FALSE)
  alpha_diversity <- estimate_richness(rarefied_physeq, measures = c("Shannon"))
  alpha_diversity$Depth <- depth
  alpha_diversity$SampleID <- rownames(alpha_diversity)
  alpha_diversity
})

# Create the alpha rarefaction plot with color coding by SampleID and a legend
alpha_rarefaction_plot <- ggplot(rarefaction_results, aes(x = Depth, y = Shannon, group = SampleID, color = SampleID)) +
  geom_line(alpha = 0.8) +
  theme_minimal() +
  labs(x = "Sequencing Depth", y = "Shannon Diversity Index", title = "Alpha Rarefaction Curve", color = "Sample ID") +
  theme(legend.position = "right") # Position the legend on the right side

# Display the plot
print(alpha_rarefaction_plot)



# EARLY ALPHA DIVERSITY WITH STATS
# Plot alpha diversity over time for each subject
plot_richness(physeqaF_early_complete, x = "NumericWeek", measures = "Shannon", color = "Name")
library(phyloseq)
library(ggplot2)

# Replace 'phyloseq_object' with the name of your phyloseq object
phyloseq_object <- physeqaF_early_complete  # Example

# Calculate richness
richness <- estimate_richness(phyloseq_object, measures = "Observed")

# Extract sample data
sample_data_df <- as(sample_data(phyloseq_object), "data.frame")

# Combine richness data with sample data
richness_data <- cbind(richness, sample_data_df)

# Plot with lines connecting the same patient across different time points
ggplot(richness_data, aes(x = Week, y = Observed, group = Name, color = Name)) +
  geom_line() +  # Add lines connecting the same patient
  geom_point() + # Add points
  labs(title = "Observed Diversity Over Time",
       x = "Timepoint",
       y = "Observed Diversity") +
  theme_minimal() +
  scale_color_discrete(name = "Patient ID")  # Optional: Rename legend

# SHANNON
# Calculate richness
richness <- estimate_richness(phyloseq_object, measures = "Shannon")

# Extract sample data
sample_data_df <- as(sample_data(phyloseq_object), "data.frame")

# Combine richness data with sample data
richness_data <- cbind(richness, sample_data_df)

# Plot with lines connecting the same patient across different time points
ggplot(richness_data, aes(x = Week, y = Shannon, group = Name, color = Name)) +
  geom_line() +  # Add lines connecting the same patient
  geom_point() + # Add points
  labs(title = "Shannon Diversity Over Time",
       x = "Timepoint",
       y = "Shannon Diversity Index") +
  theme_minimal() +
  scale_color_discrete(name = "Patient ID")  # Optional: Rename legend

# LATE ALPHA DIVERSITY WITH STATS
# Plot alpha diversity over time for each subject
plot_richness(physeqaF_late_complete, x = "NumericWeek", measures = "Shannon", color = "Name")
library(phyloseq)
library(ggplot2)

# Replace 'phyloseq_object' with the name of your phyloseq object
phyloseq_object <- physeqaF_late_complete  # Example

# Calculate richness
richness <- estimate_richness(phyloseq_object, measures = "Observed")

# Extract sample data
sample_data_df <- as(sample_data(phyloseq_object), "data.frame")

# Combine richness data with sample data
richness_data <- cbind(richness, sample_data_df)

# Plot with lines connecting the same patient across different time points
ggplot(richness_data, aes(x = Week, y = Observed, group = Name, color = Name)) +
  geom_line() +  # Add lines connecting the same patient
  geom_point() + # Add points
  labs(title = "Observed Diversity Over Time",
       x = "Timepoint",
       y = "Observed Diversity") +
  theme_minimal() +
  scale_color_discrete(name = "Patient ID")  # Optional: Rename legend

# SHANNON
# Calculate richness
richness <- estimate_richness(phyloseq_object, measures = "Shannon")

# Extract sample data
sample_data_df <- as(sample_data(phyloseq_object), "data.frame")

# Combine richness data with sample data
richness_data <- cbind(richness, sample_data_df)

# Plot with lines connecting the same patient across different time points
ggplot(richness_data, aes(x = Week, y = Shannon, group = Name, color = Name)) +
  geom_line() +  # Add lines connecting the same patient
  geom_point() + # Add points
  labs(title = "Shannon Diversity Over Time",
       x = "Timepoint",
       y = "Shannon Diversity Index") +
  theme_minimal() +
  scale_color_discrete(name = "Patient ID")  # Optional: Rename legend









# Now with some stats:

# EARLY
# OBSERVED
# Calcuearly alpha diversity (Observed)
alpha_div <- estimate_richness(physeqaF_early_complete, measures = "Observed")

# Add metadata to the alpha diversity data frame
alpha_div$SampleID <- rownames(alpha_div)
metadata <- as(sample_data(physeqaF_early_complete), "data.frame")
alpha_div <- merge(alpha_div, metadata, by.x = "SampleID", by.y = "row.names")

# Ensure the data frame has the Week and Name variables
head(alpha_div)

# Reshape the data to wide format for the paired test
library(reshape2)
alpha_div_wide <- dcast(alpha_div, Name ~ Week, value.var = "Observed")

# Perform the Wilcoxon signed-rank test
wilcox_test <- wilcox.test(alpha_div_wide$`Week00`, alpha_div_wide$`Week04`, paired = TRUE)

# Print the result
print(wilcox_test)

# early Shannon
# Calculate early alpha diversity (Shannon index)
alpha_div <- estimate_richness(physeqaF_early_complete, measures = "Shannon")

# Add metadata to the alpha diversity data frame
alpha_div$SampleID <- rownames(alpha_div)
metadata <- as(sample_data(physeqaF_early_complete), "data.frame")
alpha_div <- merge(alpha_div, metadata, by.x = "SampleID", by.y = "row.names")

# Ensure the data frame has the Week and Name variables
head(alpha_div)

# Reshape the data to wide format for the paired test
library(reshape2)
alpha_div_wide <- dcast(alpha_div, Name ~ Week, value.var = "Shannon")

# Perform the Wilcoxon signed-rank test
wilcox_test <- wilcox.test(alpha_div_wide$`Week00`, alpha_div_wide$`Week04`, paired = TRUE)

# Print the result
print(wilcox_test)


# LATE
# OBSERVED
# Calculate alpha diversity (Observed)
alpha_div <- estimate_richness(physeqaF_late_complete, measures = "Observed")

# Add metadata to the alpha diversity data frame
alpha_div$SampleID <- rownames(alpha_div)
metadata <- as(sample_data(physeqaF_late_complete), "data.frame")
alpha_div <- merge(alpha_div, metadata, by.x = "SampleID", by.y = "row.names")

# Ensure the data frame has the Week and Name variables
head(alpha_div)

# Reshape the data to wide format for the paired test
library(reshape2)
alpha_div_wide <- dcast(alpha_div, Name ~ Week, value.var = "Observed")

# Perform the Wilcoxon signed-rank test
wilcox_test <- wilcox.test(alpha_div_wide$`Week00`, alpha_div_wide$`Week08`, paired = TRUE)

# Print the result
print(wilcox_test)

# LATE Shannon
# Calculate alpha diversity (Shannon index)
alpha_div <- estimate_richness(physeqaF_late_complete, measures = "Shannon")

# Add metadata to the alpha diversity data frame
alpha_div$SampleID <- rownames(alpha_div)
metadata <- as(sample_data(physeqaF_late_complete), "data.frame")
alpha_div <- merge(alpha_div, metadata, by.x = "SampleID", by.y = "row.names")

# Ensure the data frame has the Week and Name variables
head(alpha_div)

# Reshape the data to wide format for the paired test
library(reshape2)
alpha_div_wide <- dcast(alpha_div, Name ~ Week, value.var = "Shannon")

# Perform the Wilcoxon signed-rank test
wilcox_test <- wilcox.test(alpha_div_wide$`Week00`, alpha_div_wide$`Week08`, paired = TRUE)

# Print the result
print(wilcox_test)





#### Beta Diversity ####

# install.packages("vegan")
library(phyloseq)
library(vegan)

# EARLY Beta diversity
# Replace 'phyloseq_object' with your phyloseq object
phyloseq_object <- physeqaF_early_complete

# Calculate Bray-Curtis distance matrix
bray_dist <- distance(phyloseq_object, method = "bray")


# Extract sample data
metadata <- as(sample_data(phyloseq_object), "data.frame")

# Ensure that 'Name' (patient ID) and 'NumericWeeks' (timepoint) columns are present
head(metadata)  # Check the first few rows to ensure correct format

# Perform PERMANOVA accounting for repeated measures using patient ID as strata
adonis_result <- adonis2(bray_dist ~ NumericWeek, data = metadata, strata = metadata$Name, permutations = 999)

# Print the results
print(adonis_result)

# Now for the pairwise p-values
# install.packages("pairwiseAdonis")
library(vegan)
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Generate a Bray-Curtis distance matrix from phyloseq object
distance_matrix <- phyloseq::distance(phyloseq_object, method = "bray")

# Extract metadata
metadata <- as(sample_data(phyloseq_object), "data.frame")

# Perform pairwise comparisons for a factor named "Week"
pairwise_adonis_result <- pairwise.adonis(distance_matrix, metadata$Week)

# View the results
print(pairwise_adonis_result)



# Perform PCoA
ordination <- ordinate(phyloseq_object, method = "PCoA", distance = "bray")

# Extract the ordination data
ordination_data <- as.data.frame(ordination$vectors)

# Combine ordination data with metadata
ordination_data <- cbind(ordination_data, metadata)

# Plot PCoA
library(ggplot2)
library(RColorBrewer)

# Plot PCoA with colorblind-friendly palette and 95% confidence ellipses
p <- ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = as.factor(NumericWeek), shape = as.factor(NumericWeek))) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Week") +
  scale_shape_manual(values = c(16, 17, 18, 19), name = "Week") +
  stat_ellipse(level = 0.95, alpha = 0.5) +
  theme_minimal() +
  labs(title = "PCoA of Beta Diversity (Bray-Curtis)",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

print(p)

# Now by patient rather than by time

# Ensure the 'Name' column is treated as a factor for plotting
ordination_data$Name <- as.factor(ordination_data$Name)

# Plot PCoA with colorblind-friendly palette and 95% confidence ellipses
p <- ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = Name, shape = Name)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_brewer(palette = "Set3", name = "Name") +
  scale_shape_manual(values = 1:length(unique(ordination_data$Name)), name = "Name") +
 # stat_ellipse(aes(group = Name), level = 0.95, alpha = 0.5) +
  theme_minimal() +
  labs(title = "PCoA of Beta Diversity (Bray-Curtis) Grouped by Name",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

print(p)


# LATE Beta diversity
# Replace 'phyloseq_object' with your phyloseq object
phyloseq_object <- physeqaF_late_complete

# Calculate Bray-Curtis distance matrix
bray_dist <- distance(phyloseq_object, method = "bray")


# Extract sample data
metadata <- as(sample_data(phyloseq_object), "data.frame")

# Ensure that 'Name' (patient ID) and 'NumericWeeks' (timepoint) columns are present
head(metadata)  # Check the first few rows to ensure correct format

# Perform PERMANOVA accounting for repeated measures using patient ID as strata
adonis_result <- adonis2(bray_dist ~ NumericWeek, data = metadata, strata = metadata$Name, permutations = 999)

# Print the results
print(adonis_result)

# Now for the pairwise p-values
# install.packages("pairwiseAdonis")
library(vegan)
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Generate a Bray-Curtis distance matrix from phyloseq object
distance_matrix <- phyloseq::distance(phyloseq_object, method = "bray")

# Extract metadata
metadata <- as(sample_data(phyloseq_object), "data.frame")

# Perform pairwise comparisons for a factor named "Week"
pairwise_adonis_result <- pairwise.adonis(distance_matrix, metadata$Week)

# View the results
print(pairwise_adonis_result)



# Perform PCoA
ordination <- ordinate(phyloseq_object, method = "PCoA", distance = "bray")

# Extract the ordination data
ordination_data <- as.data.frame(ordination$vectors)

# Combine ordination data with metadata
ordination_data <- cbind(ordination_data, metadata)

# Plot PCoA
library(ggplot2)
library(RColorBrewer)

# Plot PCoA with colorblind-friendly palette and 95% confidence ellipses
p <- ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = as.factor(NumericWeek), shape = as.factor(NumericWeek))) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Week") +
  scale_shape_manual(values = c(16, 17, 18, 19), name = "Week") +
  stat_ellipse(level = 0.95, alpha = 0.5) +
  theme_minimal() +
  labs(title = "PCoA of Beta Diversity (Bray-Curtis)",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

print(p)

# Now by patient rather than by time

# Ensure the 'Name' column is treated as a factor for plotting
ordination_data$Name <- as.factor(ordination_data$Name)

# Plot PCoA with colorblind-friendly palette and 95% confidence ellipses
p <- ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = Name, shape = Name)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_brewer(palette = "Set3", name = "Name") +
  scale_shape_manual(values = 1:length(unique(ordination_data$Name)), name = "Name") +
  # stat_ellipse(aes(group = Name), level = 0.95, alpha = 0.5) +
  theme_minimal() +
  labs(title = "PCoA of Beta Diversity (Bray-Curtis) Grouped by Name",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

print(p)



#### Plotting some abundances ####
library(phyloseq)
library(dplyr)

# EARLY
physeq <- physeqaF_early_complete


# Step 1: Aggregate data at the desired taxonomic rank
physeq_species <- tax_glom(physeq, taxrank = "Species")

# Step 2: Transform to relative abundance
physeq_species_rel <- transform_sample_counts(physeq_species, function(x) x / sum(x))

# Step 3: Identify the top 20 most abundant taxa
top_taxa <- names(sort(taxa_sums(physeq_species_rel), decreasing = TRUE)[1:20])

# Step 4: Label taxa as "Other" if they are not in the top 20
tax_table(physeq_species_rel) <- tax_table(physeq_species_rel) %>%
  as.data.frame() %>%
  mutate(Species = ifelse(rownames(.) %in% top_taxa, Species, "Other")) %>%
  as.matrix()

# Step 5: Plot the stacked bar chart
plot_bar(physeq_species_rel, fill = "Species") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 5))  # Adjust legend appearance

# LATE species
physeq <- physeq_late_complete

# Step 1: Aggregate data at the desired taxonomic rank
physeq_species <- tax_glom(physeq, taxrank = "Species")

# Step 2: Transform to relative abundance
physeq_species_rel <- transform_sample_counts(physeq_species, function(x) x / sum(x))

# Step 3: Identify the top 20 most abundant taxa
top_taxa <- names(sort(taxa_sums(physeq_species_rel), decreasing = TRUE)[1:20])

# Step 4: Label taxa as "Other" if they are not in the top 20
tax_table(physeq_species_rel) <- tax_table(physeq_species_rel) %>%
  as.data.frame() %>%
  mutate(Species = ifelse(rownames(.) %in% top_taxa, Species, "Other")) %>%
  as.matrix()

# Step 5: Plot the stacked bar chart
plot_bar(physeq_species_rel, fill = "Species") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 5))  # Adjust legend appearance


#### MaasLin2 ####

#BiocManager::install("Maaslin2")
library(Maaslin2)

# EARLY

physeqaF <- physeqaF_early_complete

# Extract OTU table
otu_table <- as.data.frame(otu_table(physeqaF))

# Extract taxonomy table (if needed)
taxonomy_table <- as.data.frame(tax_table(physeqaF))

# Extract metadata
metadata <- as.data.frame(sample_data(physeqaF))
metadata <- as(metadata, "data.frame")
class(metadata)

# Check the data
head(otu_table)
head(taxonomy_table)
head(metadata)

# Ensure row names match
otu_table <- otu_table[match(rownames(metadata), rownames(otu_table)), ]

# Define the input and output paths
input_data <- otu_table
class(input_data)
input_metadata <- metadata
class(input_metadata)
output <- "Maaslin2_output"

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = input_data,
  input_metadata = input_metadata,
  output = output,
  fixed_effects = c("NumericWeek"), # Add your fixed effects here
  random_effects = c("Name"), # Add your random effects here (e.g., patient name)
  normalization = "NONE", # Add normalization method if needed, e.g., "TSS" for total sum scaling
  transform = "LOG" # Add transformation method if needed
)

# Now for a fancier run at this:

# Extract taxonomy table from phyloseq object and convert to data frame
taxonomy_table <- as.data.frame(tax_table(physeqaF))

# Create a combined Genus_species column
taxonomy_table$Genus_species <- paste(taxonomy_table$Genus, taxonomy_table$Species, sep = "_")

# View the taxonomy table to ensure it is correctly formatted
head(taxonomy_table)

# Extract OTU table and convert it to a data frame
otu_table <- as.data.frame(otu_table(physeqaF))
# otu_table <- t(otu_table)  # Transpose to have samples as rows

# Ensure row names of taxonomy table match column names of OTU table
taxonomy_table <- taxonomy_table[match(colnames(otu_table), rownames(taxonomy_table)), ]

# Create a named vector with ASV names as keys and Genus_species as values
asv_to_tax <- setNames(taxonomy_table$Genus_species, rownames(taxonomy_table))

# Rename the columns of the OTU table using the taxonomy information
colnames(otu_table) <- asv_to_tax[colnames(otu_table)]

# View the updated OTU table to ensure it includes taxonomy information
head(otu_table)

# Define the input and output paths
input_data <- otu_table
input_metadata <- as(metadata, "data.frame")
class(input_metadata)
output <- "Maaslin2_output_Early"

# Specify fixed effects (e.g., timepoints) and random effects (e.g., patient ID)
fixed_effects <- c("NumericWeek")  # Replace with your fixed effects
random_effects <- c("Name")         # Replace with your random effects

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = input_data,
  input_metadata = input_metadata,
  output = output,
  fixed_effects = fixed_effects,
  random_effects = random_effects,
  normalization = "CLR",  # Add normalization method if needed
  transform = "NONE",        # Add transformation method if needed
  max_significance = 0.97,
  max_pngs = 112
)

library(ggplot2)

# Load MaAsLin2 results
results <- read.csv(file.path(output, "all_results.tsv"), sep = "\t")

# View the results to check the available columns
head(results)

# LATE

physeqaF <- physeqaF_late_complete

# Extract OTU table
otu_table <- as.data.frame(otu_table(physeqaF))

# Extract taxonomy table (if needed)
taxonomy_table <- as.data.frame(tax_table(physeqaF))

# Extract metadata
metadata <- as.data.frame(sample_data(physeqaF))
metadata <- as(metadata, "data.frame")
class(metadata)

# Check the data
head(otu_table)
head(taxonomy_table)
head(metadata)

# Ensure row names match
otu_table <- otu_table[match(rownames(metadata), rownames(otu_table)), ]

# Define the input and output paths
input_data <- otu_table
class(input_data)
input_metadata <- metadata
class(input_metadata)
output <- "Maaslin2_output"

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = input_data,
  input_metadata = input_metadata,
  output = output,
  fixed_effects = c("NumericWeek"), # Add your fixed effects here
  random_effects = c("Name"), # Add your random effects here (e.g., patient name)
  normalization = "NONE", # Add normalization method if needed, e.g., "TSS" for total sum scaling
  transform = "LOG" # Add transformation method if needed
)

# Now for a fancier run at this:

# Extract taxonomy table from phyloseq object and convert to data frame
taxonomy_table <- as.data.frame(tax_table(physeqaF))

# Create a combined Genus_species column
taxonomy_table$Genus_species <- paste(taxonomy_table$Genus, taxonomy_table$Species, sep = "_")

# View the taxonomy table to ensure it is correctly formatted
head(taxonomy_table)

# Extract OTU table and convert it to a data frame
otu_table <- as.data.frame(otu_table(physeqaF))
# otu_table <- t(otu_table)  # Transpose to have samples as rows

# Ensure row names of taxonomy table match column names of OTU table
taxonomy_table <- taxonomy_table[match(colnames(otu_table), rownames(taxonomy_table)), ]

# Create a named vector with ASV names as keys and Genus_species as values
asv_to_tax <- setNames(taxonomy_table$Genus_species, rownames(taxonomy_table))

# Rename the columns of the OTU table using the taxonomy information
colnames(otu_table) <- asv_to_tax[colnames(otu_table)]

# View the updated OTU table to ensure it includes taxonomy information
head(otu_table)

# Define the input and output paths
input_data <- otu_table
input_metadata <- as(metadata, "data.frame")
class(input_metadata)
output <- "Maaslin2_output_Late"

# Specify fixed effects (e.g., timepoints) and random effects (e.g., patient ID)
fixed_effects <- c("NumericWeek")  # Replace with your fixed effects
random_effects <- c("Name")         # Replace with your random effects

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = input_data,
  input_metadata = input_metadata,
  output = output,
  fixed_effects = fixed_effects,
  random_effects = random_effects,
  normalization = "CLR",  # Add normalization method if needed
  transform = "NONE",        # Add transformation method if needed
  max_significance = 0.97,
  max_pngs = 82
)

library(ggplot2)

# Load MaAsLin2 results
results <- read.csv(file.path(output, "all_results.tsv"), sep = "\t")

# View the results to check the available columns
head(results)


#### Exporting raw data ####

# EARLY
physeq <- physeqaF_early_complete

otu_table <- as.data.frame(otu_table(physeq))
tax_table <- as.data.frame(tax_table(physeq))
sample_data <- as(sample_data(physeq), "data.frame")

export_path <- "~/Desktop"
otu_filename <- "early_otu_table.tsv"
tax_filename <- "early_tax_table.tsv"
meta_filename <- "early_meta_table.tsv"

write.table(otu_table, file.path(export_path, otu_filename), sep="\t", row.names = TRUE, quote = FALSE)
write.table(tax_table, file.path(export_path, tax_filename), sep="\t", row.names = TRUE, quote = FALSE)
write.table(sample_data, file.path(export_path, meta_filename), sep="\t", row.names = TRUE, quote = FALSE)


# LATE
physeq <- physeqaF_late_complete

otu_table <- as.data.frame(otu_table(physeq))
tax_table <- as.data.frame(tax_table(physeq))
as.data.frame(sample_data(physeq))
sample_data <- as(sample_data(physeq), "data.frame")

export_path <- "~/Desktop"
otu_filename <- "late_otu_table.tsv"
tax_filename <- "late_tax_table.tsv"
meta_filename <- "late_meta_table.tsv"

write.table(otu_table, file.path(export_path, otu_filename), sep="\t", row.names = TRUE, quote = FALSE)
write.table(tax_table, file.path(export_path, tax_filename), sep="\t", row.names = TRUE, quote = FALSE)
write.table(sample_data, file.path(export_path, meta_filename), sep="\t", row.names = TRUE, quote = FALSE)



# Alpha Diversity

# EARLY
# Replace 'phyloseq_object' with the name of your phyloseq object
phyloseq_object <- physeqaF_early_complete

# Calculate richness
richness <- estimate_richness(phyloseq_object, measures = c("Observed", "Shannon"))

# Extract sample data
sample_data_df <- as(sample_data(phyloseq_object), "data.frame")

# Combine richness data with sample data
richness_data_early <- cbind(richness, sample_data_df)

# Write this out
desktop_path <- "~/Desktop/"
filename <- "richness_data_early.tsv"

# Write the dataframe to the TSV file
write.table(richness_data_early, file.path(desktop_path, filename), sep="\t", row.names = TRUE, quote = FALSE)

# LATE

# Replace 'phyloseq_object' with the name of your phyloseq object
phyloseq_object <- physeqaF_late_complete

# Calculate richness
richness <- estimate_richness(phyloseq_object, measures = c("Observed", "Shannon"))

# Extract sample data
sample_data_df <- as(sample_data(phyloseq_object), "data.frame")

# Combine richness data with sample data
richness_data_late <- cbind(richness, sample_data_df)

# Write this out
desktop_path <- "~/Desktop/"
filename <- "richness_data_late.tsv"

# Write the dataframe to the TSV file
write.table(richness_data_late, file.path(desktop_path, filename), sep="\t", row.names = TRUE, quote = FALSE)



#### DeSeq2 ####

library(DESeq2)
packageVersion("DESeq2")
# Turns out that there's an official tool in phyloseq to feed data directly to DeSeq2

library(phyloseq)
packageVersion("phyloseq")

# EARLY
physeqaF <- physeqaF_early_complete

sample_data(physeqaF)$Week
sample_data(physeqaF)$Name

Early_DeSeq2 = phyloseq_to_deseq2(physeqaF, ~ Name + Week)
Early_DeSeq2 = DESeq(Early_DeSeq2, test="Wald", fitType="parametric")

# Investigate the results - making a table of the taxa (with names, of the significant findings)
res = results(Early_DeSeq2, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeqaF)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Make a table of all the DeSeq2 results, add the taxa names, and export it as a csv file
res_df <- as.data.frame(res)
res_df <- cbind(as(res_df, "data.frame"), as(tax_table(physeqaF)[rownames(res_df), ], "matrix"))
write.csv(res_df, file = "~/Desktop/early_deseq2_results.csv", row.names = TRUE)


# Build some plots

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
sigtab$Genus
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
# Genus order
sigtab$Genus
sigtab$Species
x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
x = sort(x, TRUE)
sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))
ggplot(sigtab, aes(x=Species, y=log2FoldChange, color=Genus)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Making some plots

library(microbiome)
# devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
library(gghalves)
library(tidyr)

pseq.rel <- microbiome::transform(physeqaF, "compositional")
select.taxa <- rownames(sigtab) #c("Akkermansia", "Dialister")

# Replace ASV IDs in the select.taxa with Genus_Species names
tax_table <- tax_table(pseq.rel)
genus_species <- paste(tax_table[, "Genus"], tax_table[, "Species"], sep = " ")
asv_to_taxa <- setNames(genus_species, rownames(tax_table))
select.taxa.labels <- asv_to_taxa[select.taxa]

group.colors <- c("brown3", "steelblue", "grey70")
p <- plot_paired_abundances(pseq.rel,
                            select.taxa = select.taxa,
                            group = "Week",
                            group.colors = group.colors,
                            dot.opacity = 0.25,
                            dot.size = 2,
                            group.order = NULL,
                            line = "Name"
)
p
}

# LATE
physeqaF <- physeqaF_late_complete

sample_data(physeqaF)$Week
sample_data(physeqaF)$Name

Late_DeSeq2 = phyloseq_to_deseq2(physeqaF, ~ Name + Week)
Late_DeSeq2 = DESeq(Late_DeSeq2, test="Wald", fitType="parametric")

# Investigate the results - making a table of the taxa (with names, of the significant findings)
res = results(Late_DeSeq2, cooksCutoff = FALSE)
alpha = 0.25
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeqaF)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Make a table of all the DeSeq2 results, add the taxa names, and export it as a csv file
res_df <- as.data.frame(res)
res_df <- cbind(as(res_df, "data.frame"), as(tax_table(physeqaF)[rownames(res_df), ], "matrix"))

write.csv(res_df, file = "~/Desktop/late_deseq2_results.csv", row.names = TRUE)



# Build some plots

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
sigtab$Genus
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
# Genus order
sigtab$Genus
sigtab$Species
x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
x = sort(x, TRUE)
sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))
ggplot(sigtab, aes(x=Species, y=log2FoldChange, color=Genus)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Making some plots

library(microbiome)
# devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
library(gghalves)
library(tidyr)

pseq.rel <- microbiome::transform(physeqaF, "compositional")
select.taxa <- rownames(sigtab) #c("Akkermansia", "Dialister")

# Replace ASV IDs in the select.taxa with Genus_Species names
tax_table <- tax_table(pseq.rel)
genus_species <- paste(tax_table[, "Genus"], tax_table[, "Species"], sep = " ")
asv_to_taxa <- setNames(genus_species, rownames(tax_table))
select.taxa.labels <- asv_to_taxa[select.taxa]

group.colors <- c("brown3", "steelblue", "grey70")
p <- plot_paired_abundances(pseq.rel,
                            select.taxa = select.taxa,
                            group = "Week",
                            group.colors = group.colors,
                            dot.opacity = 0.25,
                            dot.size = 2,
                            group.order = NULL,
                            line = "Name"
)
p
}


#### Genus - level data ####
# EARLY
# Replace 'phyloseq_object' with the name of your phyloseq object
phyloseq_object <- physeqaF_early_complete

# Collapse the taxon abundances at the genus level
EARLY_genus <- tax_glom(phyloseq_object, taxrank = "Genus")

# View the resulting object
EARLY_genus

physeq <- EARLY_genus

otu_table <- as.data.frame(otu_table(physeq))
tax_table <- as.data.frame(tax_table(physeq))

export_path <- "~/Desktop"
otu_filename <- "early_genus_otu_table.tsv"
tax_filename <- "early_genus__tax_table.tsv"

write.table(otu_table, file.path(export_path, otu_filename), sep="\t", row.names = TRUE, quote = FALSE)
write.table(tax_table, file.path(export_path, tax_filename), sep="\t", row.names = TRUE, quote = FALSE)

# LATE
# Replace 'phyloseq_object' with the name of your phyloseq object
phyloseq_object <- physeqaF_late_complete

# Collapse the taxon abundances at the genus level
LATE_genus <- tax_glom(phyloseq_object, taxrank = "Genus")

# View the resulting object
LATE_genus

physeq <- LATE_genus

otu_table <- as.data.frame(otu_table(physeq))
tax_table <- as.data.frame(tax_table(physeq))

export_path <- "~/Desktop"
otu_filename <- "late_genus_otu_table.tsv"
tax_filename <- "late_genus_tax_table.tsv"

write.table(otu_table, file.path(export_path, otu_filename), sep="\t", row.names = TRUE, quote = FALSE)
write.table(tax_table, file.path(export_path, tax_filename), sep="\t", row.names = TRUE, quote = FALSE)


#### Checking p-adjuested values ####

top_taxa_pvals <- c(
  0.57741307,
  0.418621284,
  0.884289498,
  0.06398843,
  0.00998336,
  0.863434876,
  0.49229665,
  0.995965769,
  0.61599641,
  0.642368082,
  0.611281785
)

adjusted_pvals <- p.adjust(top_taxa_pvals, method = "BH")

# View the adjusted p-values
print(adjusted_pvals)







