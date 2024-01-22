####################################################

# Preprocessing the FCS files before importing them into Cyotree #
# Some FCS files have laser names and some antibody names #
# data transformation: Arcsinh transformation with co-factor 150 #
# median centring and scaling # 

####################################################

# Set a working directory
setwd("xxxx/xxxx/xxxx")

# Loading packages
library(CytoTree)
library(flowCore)
library(stringr)
library(scales)

## An example showing the preprocessing of FCS files with Antibody labels ##
## FCS files with Antibody labels ##
# Create an empty list to store the processed flowframes
fcs_from_file = list()
# Donors list
donor_list = c("028(IV)","029(IV)") 
# Get the directory of FCS files
fcs.files = list.files("~/xxxx/Ab_labels_raw_late_stage_samples", pattern=".fcs",full.names = TRUE)

# loop to alter the expression matrices columns with the selected list of markers
for (i in 1:length(fcs.files)) {
  
  # input the FCS files 
  fcs_from_file[i] = flowCore::read.FCS(fcs.files[i],transformation = F) 
  
  # List of new labels to be replaced
  recol = c(`FSC-A` = "FSC-A",
            `SSC-A` = "SSC-A",
            `FJComp-APC-A` = "CD68",
            `FJComp-APC-Cy7-A` = "Live/Dead",
            `FJComp-Alexa Fluor 700-A` = "CD33",
            `FJComp-Brilliant Violet 421-A` = "CD11b",
            `FJComp-Brilliant Violet 510-A` = "CD3",
            `FJComp-Brilliant Violet 605-A` = "CD15",
            `FJComp-FITC-A` = "HLA-DR",
            `FJComp-PE-A` = "PDL-1",
            `FJComp-PE-Cy7-A` = "CD14",
            `FJComp-PerCP-Cy5-5-A` = "CD11c",
            `FJComp-Qdot 585-A` = "CD1a",
            `FJComp-Qdot 625-A` = "CD16",
            `FJComp-Qdot 705-A` = "CD19/CD20",
            `FJComp-Qdot 800-A` = "CD45")
  
  # Rename the columns with the desired marker labels
  colnames(fcs_from_file[[i]]@exprs)[match(names(recol), colnames(fcs_from_file[[i]]@exprs))] = recol
  fcs_from_file[[i]]@exprs = fcs_from_file[[i]]@exprs[, recol]
  
  # Remove the live/dead and CD45 columns; since the samples were pre-gated for live/dead and CD45
  fcs_from_file[[i]]@exprs = fcs_from_file[[i]]@exprs[,!colnames(fcs_from_file[[i]]@exprs) %in% c("CD45", "Live/Dead")]
  
  # Arcsinh transformation with co-factor 150
  fcs_from_file[[i]]@exprs = asinh(fcs_from_file[[i]]@exprs/150) 
  
  # After arcsinh transformation, median centring and scaling
  center_scale = function(x) {
    scale(x, center = apply(x,2,median), scale = T)
  }
  fcs_from_file[[i]]@exprs = center_scale(fcs_from_file[[i]]@exprs)
  
  # Also selecting only the necessary markers in the parameters list as well
  re_select = c("FSC-A",
                "SSC-A",
                "FJComp-APC-A",
                "FJComp-Alexa Fluor 700-A",
                "FJComp-Brilliant Violet 421-A",
                "FJComp-Brilliant Violet 510-A",
                "FJComp-Brilliant Violet 605-A",
                "FJComp-FITC-A",
                "FJComp-PE-A",
                "FJComp-PE-Cy7-A",
                "FJComp-PerCP-Cy5-5-A",
                "FJComp-Qdot 585-A",
                "FJComp-Qdot 625-A",
                "FJComp-Qdot 705-A")
  fcs_from_file[[i]]@parameters@data = fcs_from_file[[i]]@parameters@data[fcs_from_file[[i]]@parameters@data$name %in% re_select,]
  
  # Save the new pre-processed FCS files 
  write.FCS(fcs_from_file[[i]],
            filename = paste0("~/xxxx/CLEAN_Mix 4_Tumor_NSCLC_with_comp_median_",donor_list[i],'.fcs'))
}

## An example showing the preprocessing of FCS files with Laser labels ##
### Laser name labels ###
# Create an empty list to store the processed flowframes
fcs_from_file = list()
# Donors list
donor_list = c("03","04","06","08","09")
# Get the directory of FCS 
fcs.files = list.files("~/xxxx/laser_labels_raw_late_stage_samples/", pattern=".fcs",full.names = TRUE)


# Alter the expression matrices columns with the selected list of markers
for (i in 1:length(fcs.files)) {
  
  # input the FCS files 
  fcs_from_file[i] = flowCore::read.FCS(fcs.files[i],transformation = F) 
  
  # Re-arrange the columns of the expression matrix
  re_order = c("FSC-A",
               "FSC-H",
               "FSC-W",
               "SSC-A",
               "SSC-H",
               "SSC-W",
               "FJComp-637 pmt C 670-30-A",
               "FJComp-637 pmt A 780-60-A",
               "FJComp-637 pmt B 710-40-A",
               "FJComp-405 pmt H 431-28-A",
               "FJComp-405 pmt G 470-20-A",
               "FJComp-405 pmt E 610-20-A",
               "FJComp-488 pmt E 524-24-A",
               "FJComp-561 pmt E 586-15-A",
               "FJComp-561 pmt A 780-60-A",
               "FJComp-488 pmt B 710-50-A",
               "FJComp-355 pmt G 379-28-A",
               "FJComp-355 pmt F 515-30-A",
               "FJComp-355 pmt B 740-35-A",
               "FJComp-355 pmt A 810-44-A",
               "Time")
  
  # Re-order the column by "re_order" list
  fcs_from_file[[i]] =  fcs_from_file[[i]][,re_order]
  
  # Rename the columns to match the samples with antibody labels
  colnames(fcs_from_file[[i]]@exprs) = c("FSC-A",
                                         "FSC-H",
                                         "FSC-W",
                                         "SSC-A",
                                         "SSC-H",
                                         "SSC-W",
                                         "FJComp-APC-A",
                                         "FJComp-APC-Cy7-A",
                                         "FJComp-Alexa Fluor 700-A",
                                         "FJComp-Brilliant Violet 421-A",
                                         "FJComp-Brilliant Violet 510-A",
                                         "FJComp-Brilliant Violet 605-A",
                                         "FJComp-FITC-A",
                                         "FJComp-PE-A",
                                         "FJComp-PE-Cy7-A",
                                         "FJComp-PerCP-Cy5-5-A",
                                         "FJComp-Qdot 585-A",
                                         "FJComp-Qdot 625-A",
                                         "FJComp-Qdot 705-A",
                                         "FJComp-Qdot 800-A" ,
                                         "Time" )
  
  # List of new labels to be replaced
  recol = c(`FSC-A` = "FSC-A",
            `SSC-A` = "SSC-A",
            `FJComp-APC-A` = "CD68",
            `FJComp-APC-Cy7-A` = "Live/Dead",
            `FJComp-Alexa Fluor 700-A` = "CD33",
            `FJComp-Brilliant Violet 421-A` = "CD11b",
            `FJComp-Brilliant Violet 510-A` = "CD3",
            `FJComp-Brilliant Violet 605-A` = "CD15",
            `FJComp-FITC-A` = "HLA-DR",
            `FJComp-PE-A` = "PDL-1",
            `FJComp-PE-Cy7-A` = "CD14",
            `FJComp-PerCP-Cy5-5-A` = "CD11c",
            `FJComp-Qdot 585-A` = "CD1a",
            `FJComp-Qdot 625-A` = "CD16",
            `FJComp-Qdot 705-A` = "CD19/CD20",
            `FJComp-Qdot 800-A` = "CD45")
  colnames(fcs_from_file[[i]]@exprs)[match(names(recol), colnames(fcs_from_file[[i]]@exprs))] = recol
  fcs_from_file[[i]]@exprs = fcs_from_file[[i]]@exprs[, recol]
  
  # Remove the live/dead and CD45 columns
  fcs_from_file[[i]]@exprs = fcs_from_file[[i]]@exprs[,!colnames(fcs_from_file[[i]]@exprs) %in% c("CD45", "Live/Dead")]

  #  Arcsinh transformation with co-factor 150
  fcs_from_file[[i]]@exprs = asinh(fcs_from_file[[i]]@exprs/150) 
  
  # Median centring and scaling
  center_scale = function(x) {
    scale(x, center = apply(x,2,median), scale = T)
  }
  fcs_from_file[[i]]@exprs = center_scale(fcs_from_file[[i]]@exprs)
  
  # Also selecting only the necessary markers in the parameters list as well
  
  # Re-order the list as per "re-order" list
  fcs_from_file[[i]]@parameters@data[["name"]] =  re_order
  
  # Change the laser labels to antibody labels 
  fcs_from_file[[i]]@parameters@data[["name"]] = c("FSC-A",
                                                   "FSC-H",
                                                   "FSC-W",
                                                   "SSC-A",
                                                   "SSC-H",
                                                   "SSC-W",
                                                   "FJComp-APC-A",
                                                   "FJComp-APC-Cy7-A",
                                                   "FJComp-Alexa Fluor 700-A",
                                                   "FJComp-Brilliant Violet 421-A",
                                                   "FJComp-Brilliant Violet 510-A",
                                                   "FJComp-Brilliant Violet 605-A",
                                                   "FJComp-FITC-A",
                                                   "FJComp-PE-A",
                                                   "FJComp-PE-Cy7-A",
                                                   "FJComp-PerCP-Cy5-5-A",
                                                   "FJComp-Qdot 585-A",
                                                   "FJComp-Qdot 625-A",
                                                   "FJComp-Qdot 705-A",
                                                   "FJComp-Qdot 800-A" ,
                                                   "Time" )

  # Also select only the necessary markers in the parameters list as well
  re_select = c("FSC-A",
                "SSC-A",
                "FJComp-APC-A",
                "FJComp-Alexa Fluor 700-A",
                "FJComp-Brilliant Violet 421-A",
                "FJComp-Brilliant Violet 510-A",
                "FJComp-Brilliant Violet 605-A",
                "FJComp-FITC-A",
                "FJComp-PE-A",
                "FJComp-PE-Cy7-A",
                "FJComp-PerCP-Cy5-5-A",
                "FJComp-Qdot 585-A",
                "FJComp-Qdot 625-A",
                "FJComp-Qdot 705-A")
  fcs_from_file[[i]]@parameters@data = fcs_from_file[[i]]@parameters@data[fcs_from_file[[i]]@parameters@data$name %in% re_select,]
  
  # Save the new pre-processed FCS files 
  write.FCS(fcs_from_file[[i]],
            filename = paste0("~/xxxx/CLEAN_Mix 4_Tumor_NSCLC_with_comp_median_",donor_list[i],'.fcs'))
}

