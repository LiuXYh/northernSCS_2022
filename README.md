# Main files

**RCode_..**: Main data analysis and plotting code.

**asv_table.csv**: ASV abundance table, with ASVs annotated as cyanobacteria, chloroplasts, mitochondria, and other non-heterotrophic bacterial targets removed. The copy numbers of each ASV have been corrected using the rrnDB database and normalized.

**taxonomy.csv**: ASV taxonomic annotation table, recording the classification of each ASV at the phylum, class, order, family, genus, and species levels.

**asv_seq.tre**: Phylogenetic tree file for ASVs, constructed using the representative 16S sequences of each ASV.

**env_table.csv**: Environmental factor table, including:  

    sample_id: Sample name
    
    station: Sampling station
    
    date: Sampling date
    
    time: Sampling time
    
    longitude: Sampling longitude (unit: °E)
    
    latitude: Sampling latitude (unit: °N)
    
    depth: Sampling depth (unit: m)
    
    layer: Water layer
    
    watermass: Water mass
    
    Temperature: Temperature (unit: °C)
    
    Salinity: Salinity
    
    Chla: Chlorophyll a concentration (unit: mg m⁻³)
    
    Oxygen: Dissolved oxygen concentration (unit: mg L⁻¹)
    
    pH: pH
    
    NOx: Nitrate concentration (unit: µmol L⁻¹)
    
    PO4: Phosphate concentration (unit: µmol L⁻¹)
    
    Si: Silicate concentration (unit: µmol L⁻¹)

    MLD: Mixed layer depth (unit: m)
    
    BP: Bacterial production (unit: mg C m⁻³ d⁻¹)
    
    sBP: Cell-specific bacterial production (unit: fg C cell⁻¹ d⁻¹)
    
    BA: Bacterial abundance (unit: cells mL⁻¹)
    
    HNA_LNA: Ratio of high nucleic acid (HNA) bacterial abundance to low nucleic acid (LNA) bacterial abundance
