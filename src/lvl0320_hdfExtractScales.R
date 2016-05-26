compMetaModis <- function(hdf_filepath, granule = NULL)
{
  grnl <- c("EOS_SWATH.*EV_250.*_RefSB$", "EOS_SWATH.*EV_500.*_RefSB$", 
            "EOS_SWATH.*EV_1KM.*_RefSB$", "EOS_SWATH.*EV_1KM.*_Emissive$")
  if(!grepl("1KM.", hdf_filepath)){
    if(grepl("HKM.", hdf_filepath)){
      grnl <- grnl[1:2]
    } else {
      grnl <- grnl[1]
    }
  }
  
  info <- GDALinfo(hdf_filepath, returnScaleOffset = F)
  subds <- attr(info, "subdsmdata")

  meta <- do.call("rbind", lapply(grnl, function(x){
    sds <- unlist(strsplit(subds[grep(x, subds)], "="))[2]
    sds_info <- attr(GDALinfo(sds), "mdata")

    sds_bnds <- unlist(strsplit(unlist(strsplit(sds_info[grep("band_names", sds_info)], "="))[2], ","))
    
    sds_scls_rad <- sapply(strsplit(unlist(strsplit(sds_info[grep("radiance_scales", sds_info)], "="))[2], ", "), as.numeric)
    sds_ofst_rad <- sapply(strsplit(unlist(strsplit(sds_info[grep("radiance_offsets", sds_info)], "="))[2], ", "), as.numeric)
    
    if(x != "EOS_SWATH.*EV_1KM.*_Emissive$"){
      sds_scls_ref <- sapply(strsplit(unlist(strsplit(sds_info[grep("reflectance_scales", sds_info)], "="))[2], ", "), as.numeric)
      sds_ofst_ref <- sapply(strsplit(unlist(strsplit(sds_info[grep("reflectance_offsets", sds_info)], "="))[2], ", "), as.numeric)
    } else {
      sds_scls_ref <- NA
      sds_ofst_ref <- NA
    }
    
    data.frame(bands = sds_bnds,
               scls_rad = sds_scls_rad,
               ofst_rad = sds_ofst_rad,
               scls_ref = sds_scls_ref,
               ofst_ref = sds_ofst_ref)
  }))
  
  return(meta)
}
