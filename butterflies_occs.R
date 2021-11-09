###############################################
####      Downloading occurrences of       ####
####        European butterflies           ####
###############################################


wd <- "/Users/xavi_rp/Documents/D5_FFGRCC/European_butterflies_SDMs_data"
setwd(wd)
gbif_creds <- "/Users/xavi_rp/Dropbox/GBIF_credentials/"

library(tidyr)
library(devtools)
install_github("xavi-rp/PreSPickR", 
               ref = "v2", 
               INSTALL_opts = c("--no-multiarch"))  # https://github.com/rstudio/renv/issues/162
library(PreSPickR)


## List of European butterflies (from Wiemers, et al., 2018)

# Wiemers M, Balletto E, Dincă V, Fric ZF, Lamas G, Lukhtanov V, Munguira ML, van Swaay CAM, 
# Vila R, Vliegenthart A, Wahlberg N, Verovnik R (2018) An updated checklist of the European 
# Butterflies (Lepidoptera, Papilionoidea). ZooKeys 811: 9-45. https://doi.org/10.3897/zookeys.811.28712

list.files()
butt_sps <- read.csv("28712_22.csv", header = TRUE, sep = ";")
head(butt_sps)

length(unique(butt_sps$Original.combination))

butt_sps1 <- butt_sps[butt_sps$Original.combination != "", ]

butt_sps1 <- separate(data = butt_sps1, col = Taxon, into = c("genus", "epithet", "bot", "year"), sep = " ")
# butt_sps1[c(14, 33, 47, 71, 122, 151, 162, 168, 218, 228, 232, 236, 241, 243, 247, 248, 258, 260, 262, 267), ]

butt_sps1 <- butt_sps1[, c("genus", "epithet")]
butt_sps1 <- unite(data = butt_sps1, col = species, c(genus, epithet), sep = " ", remove = TRUE)
  
head(butt_sps1)
nrow(butt_sps1)
View(butt_sps1)


## Downloading from GBIF

sp_1_key <- as.data.frame(name_backbone(name='Iphiclides podalirius'))$speciesKey

countr <- c("BE", "EL", "LT", "PT", "BG", "ES", "LU", "RO", "CZ", "FR", "HU", "SI", "DK", "HR", "MT", "SK", "DE", "IT", "NL", "FI", "EE", "CY", "AT", "SE", "IE", "LV", "PL")
countr <- sort(countr)

num_eu_occs_df <- c()
count <- 1
for(sp in butt_sps1$species){
  sp_key <- as.data.frame(name_backbone(name = sp))$usageKey
  num_eu_occs <- 0
  for(c in countr){
    num_occs <- occ_count(taxonKey = sp_key,
                          country = c,
                          from = 2000,
                          to = 2021)
    num_eu_occs <- num_eu_occs + num_occs
  }
  num_eu_occs_df <- rbind(num_eu_occs_df, data.frame(sp, sp_key, num_eu_occs))
  print(paste0(sp, " - sp ", count, ": ", num_eu_occs))
  count <- count + 1
}

write.csv(num_eu_occs_df, "Number_occs_sp_EU.csv", row.names = FALSE)


head(num_eu_occs_df)
nrow(num_eu_occs_df)
View(num_eu_occs_df)

sum(num_eu_occs_df$num_eu_occs == 0)  # 60 over 496 (total sp.) have 0 occurrences in EU-27
summary(num_eu_occs_df$num_eu_occs)

num_eu_occs_df_no0 <- num_eu_occs_df[num_eu_occs_df$num_eu_occs != 0, ]
summary(num_eu_occs_df_no0$num_eu_occs)
unique(num_eu_occs_df_no0$sp)


## Downloading the species (9) which number of occurrences are closest to the median (1483),
#  However, when cutting their extent by European coordinates, their number of occurrences can be 
#  quite reduced, or even desapeared (e.g. Pararge xiphioides, which is from the Canary Islands).
#  The number can also be much larger due to occurrences e.g. in Switzerland (e.g. Plebejus orbitulus)

head(num_eu_occs_df[order(num_eu_occs_df$num_eu_occs, decreasing = TRUE), ], 10)
nrow(num_eu_occs_df[num_eu_occs_df$num_eu_occs > 1400 & num_eu_occs_df$num_eu_occs < 1560, ])
taxons <- num_eu_occs_df[num_eu_occs_df$num_eu_occs > 1400 & num_eu_occs_df$num_eu_occs < 1560, "sp"]
taxons <- taxons[taxons != "Vanessa vulcania"]
taxons <- taxons[taxons != "Pararge xiphioides"]

t0 <- Sys.time()
GetBIF(credentials = paste0(gbif_creds, "/gbif_credentials.RData"),
       taxon_list = taxons,
       download_format = "SIMPLE_CSV",
       download_years = c(2000, 2021),
       download_coords = c(-13, 48, 35, 72), #order: xmin, xmax, ymin, ymax
       download_coords_accuracy = c(0, 100),
       rm_dupl = TRUE,
       cols2keep = c("species", "decimalLatitude", "decimalLongitude", #"elevation",
                     "gbifID",
                     "coordinateUncertaintyInMeters",
                     "countryCode", "year", 
                     #"institutionCode",	"collectionCode",
                     #"ownerInstitutionCode",
                     "datasetKey"),
       out_name = paste0("sp_records_", format(Sys.Date(), "%Y%m%d")))

Sys.time() - t0


## if GetBIF didn't manage to create/write out the data frame with presences:
taxon_dir <- getwd()
taxons <- taxons
data1 <- Prep_BIF(taxon_dir = paste0(taxon_dir, "/"),
                  taxons = taxons,
                  cols2keep = c("species", "decimalLatitude", "decimalLongitude", #"elevation",
                                "gbifID",
                                "coordinateUncertaintyInMeters",
                                "countryCode", "year", 
                                #"institutionCode",	"collectionCode",
                                #"ownerInstitutionCode",
                                "datasetKey"
                  ),
                  #cols2keep = "all",
                  rm_dupl = TRUE)

head(data1)
unique(data1$species)
table(data1$species)

print(paste0("Saving GBIF data as ", "/sp_records_20211109", ".csv"))
write.csv(data1, file = paste0("sp_records_20211109", ".csv"),
          quote = FALSE, row.names = FALSE)

