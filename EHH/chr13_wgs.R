library(rehh)
# create haplohh object from data
hh <- data2haplohh(hap_file = "chr13_wgs.hap",
                   map_file = "chr13_wgs.map",
                   allele_coding = "01",
                   min_perc_geno.mrk = 50)

# calculate EHH statistics
res <- calc_ehh(hh, mrk = "hrp3_deletion",
                include_nhaplo = TRUE,
                include_zero_values = TRUE,
                discard_integration_at_border = FALSE )

# Create EHH Plot
plot(res)

# Calculate and plot haplotype bifurcation plots 
f <- calc_furcation(hh, mrk = "hrp3_deletion")
plot(f)

# Save EHH statistics to file
write.csv(res[[3]], "chr13_wgs_ehh.csv")
