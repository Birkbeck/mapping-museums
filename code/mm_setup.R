#
# title: "Mapping Museums, Andrea Ballatore"
#

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(knitr,readr,readr,plyr,dplyr,sp,vegan,ineq,data.table,foreach,R.utils,
               rgdal,maptools,pastecs,ggplot2,gplots,Hmisc,readxl,parsedate,rgeos,
               moments,car,lmtest,QuantPsyc,psych,tmap,spdep,GWmodel,reshape2,scales,
               quantmod,classInt,tmap,grid,writexl,openxlsx,shinyjs,corrplot,measurements,
               caret,rpart,rpart.plot,randomForest,party,arules,doParallel,magick,
               aspace,spatstat,stringdist,spatstat,biscale)

R.utils::gcDLLs()
#doParallel::registerDoParallel()

ALL_MUSEUMS_N = 4229

# register parallel backend
#nodes <- detectCores()
#print(paste0("activate parallel cluster on nodes=",nodes))
#cl <- makeCluster(nodes)
#registerDoParallel(cl)
#stopImplicitCluster()
#rm(nodes)


# ------------------------------------------------------------------
# Common variables 
# ------------------------------------------------------------------
# british_grid_crs = CRS("+init=epsg:27700") # outdated in 2020
british_grid_crs = sp::CRS(SRS_string = "EPSG:27700")  
ll_crs = sp::CRS(SRS_string = "EPSG:4326")  # ll_crs = CRS("+init=epsg:4326")


# ------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------
sort_df <- function( df, col, asc=T ){
  sdf = df[ with(df, order(df[,c(col)], decreasing = !asc)), ]
  return(sdf)
}

FUN_order_factor <- function( x, ordered_levels, ordered = T ){
  stopifnot(is.factor(x))
  x = factor(x, ordered = ordered, levels=ordered_levels )
  x
}

load_uk_dataset = function(file_path){
  fn = paste0('../../049-SpatialDatasets/data/uk_data/',file_path)
  stopifnot(file.exists(fn))
  ds = readRDS(fn)
  return(ds)
}

FUN_save_museum_dataset = function(ds,file_name){
  print(paste("FUN_save_museum_dataset",file_name,class(ds)))
  fn = paste0("../datasets/museums/",file_name)
  write_tsv(ds,paste0(fn,'.tsv'))
  saveRDS(ds,paste0(fn,'.rds'))
  return(fn)
}

# function to generate a file in a folder to show when it was updated
FUN_gen_last_update_file = function( outdir ){
  unlink( paste0(outdir,'/last.update-*') )
  dd = format(Sys.time(), "%Y_%m_%d")
  fn = paste0(outdir,'/last.update-',dd,".txt")
  file.create(fn)
  print(fn)
}

# clean output folder and keep readme file
FUN_clean_folder = function(outf, exclude_files = c('README.txt') ){
  ll = list.files(outf, all.files = T, recursive = T)
  ll = ll[!ll %in% exclude_files]
  ll = ll[!ll %in% c('.','..')]
  delete_list = paste0(outf,'/',ll)
  #print(delete_list)
  unlink(delete_list,recursive = F)
  FUN_gen_last_update_file(outf)
}

diffSets <- function(a,b, label_a,label_b){
  d = setdiff(a,b)
  invd = setdiff(b,a)
  print(paste0("Found in '",label_a,"' and not in '",label_b,"' [",length(d),'/',length(a),"]: ", paste(sort(d),collapse = ' ')))
  print(paste0("Found in '",label_b,"' and not in '",label_a,"' [",length(invd),'/',length(b),"]: ", paste(sort(invd),collapse = ' ')))
  return(c(length(d),length(invd)))
}

get_plot_subtitle <- function( file_path, suffix = '' ){
  path = gsub("\\.\\.\\/","",file_path)
  return(paste0("File: ", path,suffix))
}

# Simplify SpatialPolygonsDataFrame. It does not create slivers (based on rmapshaper).
#
# retain_tolerance: 1 --> no change
#                   .5 --> keep half points
#                   0 --> remove all points
#
simplify_poly <- function( sdf, retain_tolerance ){
  print( paste0("simplify with tolerance ", retain_tolerance  ))
  stopifnot(class(sdf)=="SpatialPolygonsDataFrame")
  stopifnot(retain_tolerance > 0 && retain_tolerance < 1)
  df = sdf@data
  sz = object.size(sdf)
  n = nrow(sdf)
  geom_simpl = rmapshaper::ms_simplify(sdf, keep = retain_tolerance) # does not create slivers
  #geom_simpl <- rgeos::gSimplify(sdf, tol = tolerance, topologyPreserve = T) # does create slivers
  sdf_simpl = SpatialPolygonsDataFrame( geom_simpl, df )
  sz2 = object.size(sdf_simpl)
  print(paste0(n," geometries: new size ",round(sz2/sz*100,1),"% of input"))
  return(sdf_simpl)
}


print_count_NAs <- function(x){
  n = sum( is.na(x), na.rm=FALSE)
  paste0('> NA/missing values: ',n,' (',
               round(n/length(x)*100,3),
               '%)')
}

p_summary_list_num <- function(x){
  p_summary_list(x)
  print("DISTRIBUTION")
  print(summary(x))
}

write_geojson <- function( sdf, fn, layer_name ){
  stopifnot(proj4string(sdf)==ll_crs@projargs)
  writeOGR(obj = sdf, dsn = fn, layer=layer_name, 
           driver="GeoJSON", overwrite_layer=T, delete_dsn=T)
  gzip(fn,overwrite=T)
  print(paste0('GeoJson written: ',fn,'.gz'))
}

p_summary_list <- function( x, lab, sort_by_pc = T, top_n = NA ){
  #cat('\n')
  t = table(x, useNA = "always")
  d = data.frame( VAL = t )
  d$PC = round((d$VAL.Freq / sum(d$VAL.Freq)*100),2)
  if (sort_by_pc)
    d = sort_df(d,'PC',asc = F)
  d$VAR = lab
  stopifnot(nrow(d)>0)
  if (!is.na(top_n))
    return(head(d,top_n))
  
  return(d)
  #print_count_NAs(x)
  #return(t)
}


FUN_gen_autocorrel_maps = function(sdf, fld, fout){
  fout = paste0(fout, "-", fld)
  # Map
  # P values
  sdf_map <- tm_shape(sdf) + 
    tm_fill(paste0(fld, "_locg_pval"), breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1)) + 
    tm_layout(aes.palette = list(seq = "-YlOrRd")) + 
    tm_layout(paste0("Local Gi* p value, ", fld)) +
    tm_shape(sdf) + tm_borders(col="white", lwd=.2)
  tmap_save(sdf_map, paste0(fout,"-gi_pval_map.pdf"))
  tmap_save(sdf_map, paste0(fout,"-gi_pval_map.png"))
  
  #tmap_save(lsa_map, paste0("../papers/2019-TGIS-LA_vgi_demography/figures/la_locspautoc_p-",lsa_var_name,"-v1.pdf"), width = 5, height = 5)
  # Hot and cold spot categories
  sdf_map <- tm_shape(sdf) + 
    tm_fill (paste0(fld, "_locg_fill")) + 
    tm_add_legend(
      labels = c("Hotspot 99% conf intrv", "Hotspot 95% conf intrv", "Hotspot 90% conf intrv",
                 "Not significant", "Coldspot 90% conf intrv", "Coldspot 95% conf intrv",
                 "Coldspot 99% conf intrv"),
      col = c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4")) + 
    tm_layout(paste0("Local Gi* classification, ", fld)) +
    tm_shape(sdf) + tm_borders(col="white", lwd=.2)
  tmap_save(sdf_map, paste0(fout,"-gi_map.pdf"))
  tmap_save(sdf_map, paste0(fout,"-gi_map.png"))
  
  # Map
  # P values
  sdf_map <- tm_shape(sdf) + 
    tm_fill(paste0(fld, "_loci_pval"), breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1)) + 
    tm_layout(aes.palette = list(seq = "-YlOrRd")) + 
    tm_layout(paste0("Moran's I p value, ", fld)) +
    tm_shape(sdf) + tm_borders(col="white", lwd=.2)
  tmap_save(sdf_map, paste0(fout,"-morani_pval_map.pdf"))
  tmap_save(sdf_map, paste0(fout,"-morani_pval_map.png"))
  
  # P values
  sdf_map <- tm_shape(sdf) + 
    tm_fill(paste0(fld, "_loci"), n=5) + # breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1) 
    tm_layout(aes.palette = list(seq = "-YlOrRd")) + 
    tm_layout(paste0("Local Moran's I, ", fld)) +
    tm_shape(sdf) + tm_borders(col="white", lwd=.2)
  tmap_save(sdf_map, paste0(fout,"-morani_map.pdf"))
  tmap_save(sdf_map, paste0(fout,"-morani_map.png"))
}



# Getis-Ord local Gi*
FUN_localGiStar <- function (sppoly_obj, var_name, neighbors = "queen", k = 1) {
  print(paste("FUN_localGiStar", nrow(sppoly_obj)))
  
  # Compute neightbourhood
  if (neighbors == "queen") {
    sppoly_obj_nb <- poly2nb(sppoly_obj)
  } else if (neighbors == "knn") {
    sppoly_obj_nb <- knn2nb(knearneigh(coordinates(sppoly_obj), k = k))
  } else {
    sppoly_obj_nb <- poly2nb(sppoly_obj)
  }
  
  # Get spatial weights for neighbours lists
  sppoly_obj_lw <- nb2listw(
    include.self(sppoly_obj_nb),
    style="B")
  
  # Computer Local Gi test
  localG_results <- localG(sppoly_obj@data[, var_name], sppoly_obj_lw)
  
  # Simulated p-values
  sim_G <- matrix(0, 10000, length(localG_results))
  for (i in 1:10000) {
    sim_G[i,] <- localG(sample(sppoly_obj@data[, var_name]), sppoly_obj_lw)
  }
  sim_G_pval <- ifelse(
    localG_results > 0,
    (colSums(sweep(sim_G, 2, localG_results, '>=')) + 1) / (nrow(sim_G) + 1),
    (colSums(sweep(sim_G, 2, localG_results, '<=')) + 1) / (nrow(sim_G) + 1))
  sim_G_pval <- p.adjust(sim_G_pval, method = "fdr")
  
  # Add to shapefile
  sppoly_obj$locg <- as.numeric(localG_results)
  sppoly_obj$locg_pval <- sim_G_pval
  sppoly_obj$locg_clas <- ifelse(
    sppoly_obj$locg_pval < 0.01,
    ifelse(
      sppoly_obj$locg >0, "Hotspot", "Coldspot"
    ),
    "Not_sign"
  )
  sppoly_obj$locg_nclas <- ifelse(
    sppoly_obj$locg_pval < 0.01,
    ifelse(
      sppoly_obj$locg >0, 1, -1
    ),
    0
  )
  # c("#", "#", "#", "#", "#", "#", "#")
  sppoly_obj$fill <- ifelse(
    sppoly_obj$locg_pval < 0.01,
    ifelse(
      #sppoly_obj$locg >0, "#b2182b", "#2166ac"
      sppoly_obj$locg >0, "#d73027", "#4575b4"
    ),
    ifelse(
      sppoly_obj$locg_pval < 0.05,
      ifelse(
        #sppoly_obj$locg >0, "#ef8a62", "#67a9cf"
        sppoly_obj$locg >0, "#fc8d59", "#91bfdb"
      ),
      ifelse(
        sppoly_obj$locg_pval < 0.1,
        ifelse(
          #sppoly_obj$locg >0, "#fddbc7", "#d1e5f0"
          sppoly_obj$locg >0, "#fee090", "#e0f3f8"
        ),
        #"#f7f7f7"
        "#ffffbf"
      )
    )
  )
  colnames(sppoly_obj@data)[colnames(sppoly_obj@data) == "locg"] <- paste0(var_name, "_locgistar")
  colnames(sppoly_obj@data)[colnames(sppoly_obj@data) == "locg_pval"] <- paste0(var_name, "_locg_pval")
  colnames(sppoly_obj@data)[colnames(sppoly_obj@data) == "locg_clas"] <- paste0(var_name, "_locg_clas")
  colnames(sppoly_obj@data)[colnames(sppoly_obj@data) == "locg_nclas"] <- paste0(var_name, "_locg_nclas")
  colnames(sppoly_obj@data)[colnames(sppoly_obj@data) == "fill"] <- paste0(var_name, "_locg_fill")
  
  # Return object
  sppoly_obj
  
}

# returns a df containing generic stats about the df
FUN_get_numeric_df_summary <- function(df, normStats = F){
  nums <- sapply(df, is.numeric)
  if (!any(nums)){
    print(paste0("  no numeric cols found, skipping."))
    return(NA)
  }
  sumdf = round( pastecs::stat.desc(df[,nums], norm=normStats), 4)
  sumdf = cbind( row.names(sumdf), sumdf )
  names(sumdf)[1] = 'STAT'
  return(sumdf)
}

FUN_copy_doc_uk <- function( file_path, out_folder ){
  fn = paste0('../../049-SpatialDatasets/data/uk_data_source/',file_path)
  stopifnot(file.exists(fn))
  filen = tolower(basename(fn)) #tools::file_path_sans_ext(
  new_file = paste0(out_folder,'doc-',filen)
  file.copy(fn,new_file)
  return(new_file)
}

# compare two columns in a df
# NB: COMPLEX and BRITTLE FUNCTION
FUN_gen_freq_matrix <- function( df, var1, var2, outf, doHeatmap = T, rowPercents = NA, rowPercentsOmit0 = T ){
  #print(paste("FUN_gen_freq_matrix", var1, var2, outf, doHeatmap, length(rowPercents), rowPercentsOmit0))
  stopifnot(class(df)=='data.frame')
  
  clean_var_names = function(x){
    x = gsub("\\_+", " ", x)
    x = make.names(trimws(x))
    x = gsub("\\.+", ".", x)
    x = gsub("NA\\.", "NOT_AVAIL", x)
    x = gsub("^X", "", x)
    x = gsub("\\.", " ", x)
    return(x)
  }
  
  gen_count_matrix = function(freqm,var1,var2){
    freqm$ROW_SUM = rowSums(freqm) # sum rows
    cols = names(freqm)
    cols[is.na(cols)] = "NOT_AVAIL"
    cols = as.character(cols)
    names(freqm) = clean_var_names(cols)
    rown = row.names(freqm)
    rown[is.na(rown)] = 'NOT_AVAIL'
    stopifnot( !any(is.na(rown)) )
    row.names(freqm) = clean_var_names(rown)
    rownames_df = data.frame(X=rown)
    freqmout = cbind( rownames_df, freqm )
    # data is ready (freqmout)
    freqmout[,1] = clean_var_names(freqmout[,1])
    names(freqmout) = clean_var_names(names(freqmout))
    names(freqmout)[1] = "*"
    return(freqmout)
  }
  
  # important: this calculates the DIVERGENCE by COLUMN (and not ROW)
  gen_divergence_matrix = function(freqmout, rowPercents){
    rowPercents$VAL.x = clean_var_names(rowPercents$VAL.x)
    stopifnot( length(freqmout$`*`)>0 )
    stopifnot( length(intersect(unique(freqmout$`*`), unique(rowPercents$VAL.x))) > 0)
    # merge data
    m = merge(freqmout, rowPercents, by.x='*', by.y="VAL.x", all.x=T, all.y=F)
    stopifnot(c("*","PC") %in% names(m))
    resm = m[,c('*','PC')]
    stopifnot( nrow(resm)>0, sum(resm$PC,na.rm = T)>0 )
    cols = names(m[,2:(ncol(m)-2)])
    stopifnot(length(cols)>0)
    for (c in cols){
      if(rowPercentsOmit0){
        catpc = ifelse( m[,c]==0, NA, round( m[,c] / sum(m[,c],na.rm = T) * 100, 3) )
      } else {
        catpc = round( m[,c] / sum(m[,c],na.rm = T) * 100, 3)
      }
      resm[,paste0(c,'_div')] = round( catpc - m$PC, 3)
    }
    return(resm)
  }
  
  stopifnot(var1 %in% names(df), var2 %in% names(df), nrow(df)>0)
  if (!all(is.na(rowPercents))){
    stopifnot( c('VAL.x','PC') %in% names(rowPercents) )
  }
  # Calc COUNTS
  # format values for merge (bug fix)
  df[,var1] = make.names(df[,var1])
  df[,var2] = make.names(df[,var2])
  tabfreq = table( df[,c(var2,var1)], useNA = "ifany" )
  freqm = as.data.frame.matrix(tabfreq)
  freqmout = gen_count_matrix(freqm)
  write_xlsx(freqmout, paste0(outf,'--data.xlsx'))
    
  if(doHeatmap){
    # print plain heatmap of values
    freqm$ROW_SUM = NULL # remove column
    mat = as.matrix(freqm)
    colnames(mat) = clean_var_names( colnames(mat) )
    rownames(mat) = clean_var_names( rownames(mat) )
    FUN_df_to_heatmap( mat, var1, var2, paste0(outf,'--heat.pdf') )
    rm(mat)
  }
  
  # calculate DIVERGENCES based on the rowPercents
  if (!all(is.na(rowPercents))){
    resm = gen_divergence_matrix(freqmout, rowPercents)
    write_xlsx(resm, paste0(outf,'--diverg.xlsx'))
    
    # divergence heatmap
    if (doHeatmap){
      # generate heatmap of divergence values
      resm = resm[ !is.na(resm[,2]), ] 
      stopifnot(nrow(resm)>0,ncol(resm)>0)
      mat = round(as.matrix(resm[,2:ncol(resm)]),1)
      # col labels
      colnames(mat)[1] = '[Areas %]'
      colnames(mat) = gsub( '_div', ' ', colnames(mat))
      colnames(mat) = trimws(gsub( '_', ' ', colnames(mat)))
      # row labels
      rownames(mat) = resm$`*`
      
      stopifnot( sum(!is.na(mat[,2])) > 0 ) # ensure non-NA values in the matrix
      FUN_df_to_heatmap(mat, var1, var2,
                        paste0(outf,'--diverg_heat.pdf'),
                        colPalette = c("red", 'white', "green"))
      rm(mat,rown)
    }
    rm(m,resm)
  }
  return(freqmout)
}

FUN_dfsummary_to_file <- function(df,fout){
  sink(fout)
  print("== Summary statistics ==")
  print(paste("rows=",nrow(df),"cols=",ncol(df)))
  print(summary(df))
  sink()
}

# generate a heatmap from a DF, similar to an excel spreadsheet with coloured cells
FUN_df_to_heatmap <- function( num_mat, xlab, ylab, fout, colPalette = c("white", "royalblue") ){
  print(paste("FUN_df_to_heatmap",xlab,ylab,fout))
  row_labs = rownames(num_mat)
  col_labs = colnames(num_mat)
  stopifnot(nrow(num_mat)>0,ncol(num_mat)>0, class(num_mat)=='matrix')
  # get max lengths of labels
  maxchar1 <- max( nchar(as.character(col_labs)),na.rm = T )
  maxchar2 <- max( nchar(as.character(row_labs)),na.rm = T )
  valn1 = length(unique(as.character(col_labs)))
  valn2 = length(unique(as.character(row_labs)))
  #print(paste("maxchar1",var1,maxchar2))
  #print(paste("maxchar2",var2,maxchar2))
  stopifnot( !is.na(maxchar1), !is.na(maxchar2), maxchar1>0, maxchar2>0, valn1>0, valn2>0)
  # Based on http://sebastianraschka.com/Articles/heatmaps_in_r.html
  w = h = 10
  marg_w = marg_h = 20
  pdf( fout, width=w, height=h ) 
  my_palette <- colorRampPalette(colPalette)(n = 100)
  sz = par("cex.main")
  par(cex.main=.7)
  
  calcSz = function(x,n){
    #print(paste('calcSz',x,n))
    stopifnot( is.numeric(x), !is.na(x) )
    sz = 2
    if (n>=30) return(sz*.2)
    
    if (x>=120) return(sz*.2)
    if (x>=80) return(sz*.3)
    if (x>=40) return(sz*.5)
    if (x>=20) return(sz*.7)
    return(sz)
  }
  #print(paste(maxchar1,maxchar2))
  colLabSz = calcSz(maxchar1,valn1)
  rowLabSz = calcSz(maxchar2,valn2)
  #print(paste0("heatmap.2 valn1=",valn1," valn2=",valn2," rowLabSz=",rowLabSz," colLabSz=",colLabSz))
  heatmap.2(num_mat,
            cellnote = num_mat,  # same data set for cell labels
            main = paste(xlab,"vs",ylab,"\nFile:",fout), # heat map title
            linecol = 'gray',
            cexRow = rowLabSz, 
            cexCol = colLabSz,
            notecex = 0.5,
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            xlab = xlab,
            ylab = ylab,
            key = F, # hide legend
            margins =c(marg_w,marg_h),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="none",     # only draw a row dendrogram
            srtCol=45,
            Rowv="NA",
            Colv="Rowv")            # turn off column clustering
  dev.off()    
  par(cex.main=sz)
  rm(marg_h,marg_w)
}

scaled.box.cox <- function(x,bc.lambda){
  bc.x <- ( ( x ** bc.lambda ) - 1) / bc.lambda
  s.bc.x <- scale(bc.x, center = TRUE, scale=TRUE)
  return(s.bc.x)
}

bivariate.choropleth.colors <- function(code){
  if (code == "RdBu") return(c("#e8e8e8", "#e4acac", "#c85a5a","#b0d5df", "#ad9ea5", "#985356","#64acbe", "#627f8c", "#574249"))
  if (code == "BuPu") return(c("#e8e8e8", "#ace4e4", "#5ac8c8","#dfb0d6", "#a5add3", "#5698b9", "#be64ac", "#8c62aa", "#3b4994"))
  if (code == "GnBu") return(c("#e8e8e8", "#b5c0da", "#6c83b5", "#b8d6be", "#90b2b3", "#567994", "#73ae80", "#5a9178", "#2a5a5b"))
  if (code == "PuOr") return(c("#e8e8e8", "#e4d9ac", "#c8b35a","#cbb8d7", "#c8ada0", "#af8e53","#9972af", "#976b82", "#804d36"))
  stopifnot(FALSE)
}
  
bivariate.choropleth <- function(bivmap.dataset, bivmap.vars, bivmap.labels, fout, 
                                 bivmap.style='quantile', bvColorsCode="RdBu"){
  bvColors = bivariate.choropleth.colors(bvColorsCode)
  
  bivmap.sdf <- bivmap.dataset[
    !is.na(bivmap.dataset@data[,bivmap.vars[1]]) &
    !is.na(bivmap.dataset@data[,bivmap.vars[2]]) &
    !is.infinite(bivmap.dataset@data[,bivmap.vars[1]]) &
    !is.infinite(bivmap.dataset@data[,bivmap.vars[2]])
    ,bivmap.vars]
  colnames(bivmap.sdf@data) <- c("xvar","yvar")
  # output bins
  sink(paste0(fout,"-bins.txt"))
  print("VAR X")
  print(classIntervals( bivmap.sdf@data$xvar, n=3, bivmap.style))
  print("VAR Y")
  print(classIntervals( bivmap.sdf@data$yvar, n=3, bivmap.style))
  sink()
  bivmap.sdf@data$xcat <- findCols(classIntervals( bivmap.sdf@data$xvar, n=3, bivmap.style))
  bivmap.sdf@data$ycat <- findCols(classIntervals( bivmap.sdf@data$yvar, n=3, bivmap.style))
  bivmap.sdf@data$bicat <- bivmap.sdf@data$xcat + (3 * (bivmap.sdf@data$ycat - 1))
  
  bivmap <- tm_shape(bivmap.sdf) +
    tm_fill("bicat", style="cat", palette=bvColors, border.alpha=0,
            legend.show=T, legend.hist = F) +
    tm_scale_bar(width=0.15, position=c("left","bottom")) + tm_credits("\n\n\n\n\nContains data from Office for National Statistics\nand Ordnance Survey. Crown copyright & database right 2016.",size=0.6,position=c("left","bottom")) + tm_layout(frame=FALSE) + tm_legend(scale=0.75)
  suppressWarnings(print( bivmap ))
  #tmap_save(bivmap, fout)
  print(fout)
  pdf(paste0(fout,"-legend.pdf"))
  vp <- viewport(x=.85, y=.15, width=.3, height=.3)
  pushViewport(vp)
  print(levelplot(matrix(1:9, nrow=3), axes=FALSE, col.regions=bvColors,
                  xlab=list(label=bivmap.labels[1],cex=0.7), ylab=list(label=bivmap.labels[2],cex=0.7),
                  cuts=8, colorkey=FALSE,
                  scales=list(draw=0)),
        newpage=FALSE)
  popViewport()
  dev.off()
  return(bivmap)
}

is_unique <- function(x){
  stopifnot(length(x)>0)
  return(!anyDuplicated(x))
}

# Find museums that are open in the given year
# (museums without closing date + museums that are still open, based on closing MEAN and not END,
# otherwise many ranges like 1970;2017 will be included)
FUN_get_museums_open_in_a_year <- function( all_mus, year ){
  sdf = subset( all_mus, all_mus$year_opened_END <= year &
                  (is.na(all_mus$year_closed_END) | all_mus$year_closed_MEAN >= year))
  return(sdf)
}

FUN_get_random_linetypes = function( vars, force=F ){
  stopifnot(length(vars)>=0)
  min_threshold_n = 4
  n = length(unique(vars))
  if (n <= min_threshold_n && !force){
    styles = rep(c('solid'), length(vars))
  } else {
    # http://www.sthda.com/english/wiki/ggplot2-line-types-how-to-change-line-types-of-a-graph-in-r-software
    mapvals = rep(c('solid','longdash','dotdash'), n/3+1)
    styles = mapvalues(vars, 
              from=unique(vars), 
              to=mapvals[1:n])
  }
  #stopifnot(length(styles)==n_cases)
  return(styles)
}

# get probability that a museum is open at the target year
#  updated to formula from Nick (May 2018)
#   based on MMvisNotes_AP_NL_01_05-3_with_code.docx
#   1.	open at a given time [ Bar ]  
#
# Open at a given time Logic: Suppose the time selected is t; 
#   for each museum with an opening date (fo,to) 
#   and a closing date (fc,tc), we add to the count for time t as follows:
#   If t < fo then 0
#   If fo <= t <= to then t-fo+1 / to-fo+1
#   If to < t < fc then 1
#   If fc <= t <= tc then tc-t+1 / tc-fc+1 
#   If t> tc then 0
FUN_get_probs_open_at_given_year = function( fo, to, fc, tc, t, opt='prob', prjid=NA ){
  stopifnot(opt %in% c('prob','min','max')) # min and max are to get the bounds

  stopifnot( t > 1700 && t < 2100 )
  if (is.na(fo) || is.na(to)) return(0.0) # handle cases
  stopifnot(is.numeric(fo) && is.numeric(to))
  #if (fo > to) print(paste(fo,to))
  stopifnot( fo <= to )
  
  # get bounds
  if (opt %in% c('min','max')){
    return(FUN_get_probs_open_at_given_year_bounds(fo, to, fc, tc, t, opt))
  }
  
  # normal, get probability
  if (t < fo) return(0.0)
  if (fo <= t && t <= to){
    p = (t-fo+1) / (to-fo+1)
    stopifnot( p>=0 && p<=1 )
    return(p)
  }
  bOpen = is.na(fc)
  if (bOpen){
    if (to < t) return(1.0)
  } else {
    #print(paste(fc,tc))
    stopifnot( fc <= tc || is.na(tc) )
    if (is.na(tc)) tc = fc # dirty fix for missing fc values
    if (to < t && t < fc) return(1.0)
    if (fc <= t && t <= tc){
      p = (tc-t+1) / (tc-fc+1)
      #print(paste(fo, to, fc, tc, t, p))
      stopifnot( p>=0 && p<=1 )
      return(p)
    }
    if (t > tc) return(0.0)
  }
  stopifnot(F) # this should never be reached
  return(NA) 
  
  #df = mdf@data
  #df$tmp1930 = FUN_get_openclose_probs_for_year( df, 1930 )
  #df$tmp1960 = FUN_get_openclose_probs_for_year( df, 1960 )
  #df$tmp1975 = FUN_get_openclose_probs_for_year( df, 1975 )
  #df$tmp1980 = FUN_get_openclose_probs_for_year( df, 1980 )
  #df$tmp1995 = FUN_get_openclose_probs_for_year( df, 1995 )
  #df$tmp2005 = FUN_get_openclose_probs_for_year( df, 2005 )
  #df$tmp2017 = FUN_get_openclose_probs_for_year( df, 2017 )
  #res = df[,c("year_opened_BEG","year_opened_END",
  #            "year_closed_BEG","year_closed_END",'tmp1960','tmp1980','tmp2005','tmp2017')]
  #colSums(res)
  #summary(res)
  #nrow(res)
  #write_xlsx(res,'../tmp/year_open_logic_test.xlsx')
  #rm(res,df)
}

FUN_get_probs_open_at_given_year_bounds = function( fo, to, fc, tc, t, opt ){
  stopifnot(opt %in% c('min','max')) # min and max are to get the bounds
  p = FUN_get_probs_open_at_given_year(fo, to, fc, tc, t)  
  if (p > 0 & p < 1){
    if (opt == 'min') return(0.0)
    if (opt == 'max') return(1.0)
  }
  p
}

# 1.	open at a given time
FUN_get_open_at_given_time = function(df, year, opt='prob'){
  stopifnot(opt %in% c('prob','min','max'))
  
  d = df[,c("year_opened_BEG","year_opened_END",
        "year_closed_BEG","year_closed_END")]
  year_open_prob = apply(d, 1,
      function(x){ FUN_get_probs_open_at_given_year( x[1], x[2], x[3], x[4], year, opt) }) 
  stopifnot(nrow(df)==length(year_open_prob))
  
  return( round(year_open_prob, 4) )
}

# 1.	close at a given time 
FUN_get_close_at_given_time = function(df, year, opt='prob'){
  stopifnot(opt %in% c('prob','min','max'))
  d = df[,c("year_opened_BEG","year_opened_END",
            "year_closed_BEG","year_closed_END","project_id")]
  year_close_prob = apply(d, 1,
        function(x){ FUN_get_probs_open_at_given_year( as.numeric(x[3]), as.numeric(x[4]), NA, NA, year, opt) })
  stopifnot(nrow(df)==length(year_close_prob))
  return( round(year_close_prob, 4) )
}

# 4.	openings over time (aka stats by years)
#   If t < fo then 0
#   If fo <= t <= to then 1 / to-fo+1
#   If to < t then 0
FUN_get_openings_over_time = function(df, year){
  probs = apply(df[,c("year_opened_BEG","year_opened_END")], 1, function(x){ 
    t = year
    fo = x[1]
    to = x[2]
    stopifnot(fo<=to)
    if (t<fo) return(0.0)
    if (fo <= t && t <= to) return( 1.0 / (to-fo+1))
    if (to<t) return(0.0)
  })
  stopifnot(nrow(df)==length(probs))
  return( round(probs, 4) )
}

# 5.	closings over time (aka stats by years)
#   If t < fc then 0
#   If fc <= t <= tc then 1 / tc-fc+1
#   If tc < t then 0
FUN_get_closings_over_time = function(df, year){
  probs = apply(df[,c("year_closed_BEG","year_closed_END")], 1, function(x){ 
    t = year
    fc = x[1]
    tc = x[2]
    if (is.na(fc)) return(0.0)
    if (is.na(tc)) tc=fc # fix
    stopifnot(fc<=tc)
    
    if (t<fc) return(0.0)
    if (fc <= t && t <= tc) return( 1.0 / (tc-fc+1) )
    if (tc<t) return(0.0)
  })
  stopifnot(nrow(df)==length(probs))
  return( round(probs, 4) )
}

# execute function 'funct' on subgroups of 'df', based on columns 'splitcols'
FUN_subgroups <- function( df, splitcols, funct ){
  stopifnot(nrow(df)>0)
  resdf <- foreach(block=split(df, df[,splitcols]), .combine='rbind', .errorhandling="stop") %do% {
    if (nrow(block) == 0) return(NULL)
    return( funct(block) )
  }
  return(resdf)
}

# execute 'funct' on each element in 'elements'
FUN_foreach_in <- function( elements, funct ){
  stopifnot(length(elements)>0)
  resdf <- foreach(i=elements, .combine='rbind', .errorhandling="stop") %do% {
    return( funct(i) )
  }
  return(resdf)
}

FUN_count_groups_by_vars = function( df, vars ){
  stopifnot(length(vars)>=1)
  stopifnot(nrow(df)>=1)
  df = as_tibble(df)
  df = df[,vars]
  res = df %>% group_by_all() %>% tally()
  res
}

FUN_calc_percentage_change <- function( timeseries ){
  ch = (timeseries/lag(timeseries) - 1)*100
  stopifnot(length(timeseries)==length(ch))
  return(ch)
}

FUN_calc_percentage_change_df = function(df, excl_cols){
  df1 = df
  for (c in excl_cols) df1[,c]=NULL
  res = cbind( apply(df1, 2, function(x){
    FUN_calc_percentage_change(x)
  }) )
  res = round(res,1)
  res[is.nan(res)] = NA
  res[is.infinite(res)] = NA
  res = cbind( df[,excl_cols], res )
  return(res)
}

FUN_col_exists <- function(df, colnames){
  return(all(colnames %in% colnames(df)))
}

FUN_countNotNA <- function(x){
  length(which(!is.na(x)))
}

FUN_countNA <- function(x){
  length(which(is.na(x)))
}

FUN_twoVarStats <- function(musdf, var1, var2, prefixfn){
  stopifnot(nrow(musdf)>0)
  print(paste("FUN_twoVarStats",var1,var2))
  df = as.data.frame.matrix(table(musdf[,c(var1,var2)]))
  
  if (var2 == "classification_1998"){
    names(df) = paste0("98-",names(df))
    names(df)[1] = '98-NA'
  }
  
  # get percentages
  pcdf = round(df / sum(df,na.rm=T)*100, 3)
  pcdf$ROW_TOT = rowSums(pcdf, na.rm = T)
  pcdf = rbind( pcdf, colSums(pcdf,na.rm = T) )
  
  `*` = row.names(df)
  df = cbind( `*`, df )
  pcdf = cbind( c(`*`,"COL_TOT"), pcdf )
  names(pcdf)[1]= "*"
  
  fn = paste0(prefixfn,var1,'_vs_',var2)
  FUN_df_to_heatmap(as.matrix(df[,seq(2,ncol(df))]), var2, var1, paste0(fn,"--heat.pdf"))
  write_xlsx(df,paste0(fn,".xlsx"))
  write_xlsx(pcdf,paste0(fn,"-PC.xlsx"))
  return(df)
}

get_visits_by_museum_size_mean = function(x){
  
  if (x=="small") return(5000)
  if (x=="medium") return(20000)
  if (x=="large") return(475000)
  if (x=="huge") return(1000000)
  return(0)
}

get_visits_by_museum_size_min = function(x){
  # 90% of range (exclude extreme cases)
  if (x=="small") return(1000)
  if (x=="medium") return(10500)
  if (x=="large") return(60000)
  if (x=="huge") return(1000000)
  return(0)
}

get_visits_by_museum_size_max = function(x){
  if (x=="small") return(9000)
  if (x=="medium") return(40000)
  if (x=="large") return(900000)
  if (x=="huge") return(2000000)
  return(0)
}

# gen spatial stats for museums, based on variable UNIT_GSS in the museum df.
# SLOW function.
FUN_summarise_museum_by_spatial_unit <- function( museum_df ){
  print("FUN_summarise_museum_by_spatial_unit")
  stopifnot(class(museum_df)=='data.frame', nrow(museum_df)>0, 
            "UNIT_GSS" %in% names(museum_df), is.factor(museum_df$UNIT_GSS) )
  
  cols = c("UNIT_GSS","VAR","VAL","n_museums")
  
  # ALL
  df = plyr::ddply(museum_df, .(UNIT_GSS), summarise, n_museums=length(project_id),
                   .drop=FALSE, .inform=T, .parallel=F)
  df$VAR = "all"
  df$VAL = "all"
  df = df[,c("UNIT_GSS","VAR","VAL","n_museums")]
  names(df) = cols
  
  # GOVERNANCE
  df1 <- plyr::ddply(museum_df, .(UNIT_GSS, governance_simpl), summarise, n_museums=length(project_id), .drop=FALSE)
  df1$VAR = "governance_simpl"
  df1 = df1[,c("UNIT_GSS","VAR","governance_simpl","n_museums")]
  names(df1) = cols
  
  # SIZE
  df2 <- plyr::ddply(museum_df, .(UNIT_GSS, size), summarise, n_museums=length(project_id), .drop=FALSE)
  df2$VAR = "size"
  df2 = df2[,c("UNIT_GSS","VAR","size","n_museums")]
  names(df2) = cols
  
  # SPECIAL TYPE (for geo paper)
  dfsp <- plyr::ddply(museum_df, .(UNIT_GSS, special_type), summarise, n_museums=length(project_id), .drop=FALSE)
  dfsp$VAR = "special_type"
  dfsp = dfsp[,c("UNIT_GSS","VAR","special_type","n_museums")]
  names(dfsp) = cols
  
  # VISITS (based on SIZE) min
  museum_df$size_visits_min = sapply(museum_df$size, get_visits_by_museum_size_min)
  dfv <- plyr::ddply(museum_df, .(UNIT_GSS), summarise, n_museums=sum(size_visits_min), .drop=FALSE)
  dfv$VAR = "size_visits_min"
  dfv$visit_aggr = "sum.min"
  dfv1 = dfv[,c("UNIT_GSS","VAR","visit_aggr","n_museums")]
  names(dfv1) = cols
  
  # VISITS (based on SIZE) mean
  museum_df$size_visits_mean = sapply(museum_df$size, get_visits_by_museum_size_mean)
  dfv <- plyr::ddply(museum_df, .(UNIT_GSS), summarise, n_museums=sum(size_visits_mean), .drop=FALSE)
  dfv$VAR = "size_visits_mean"
  dfv$visit_aggr = "sum.mean"
  dfv2 = dfv[,c("UNIT_GSS","VAR","visit_aggr","n_museums")]
  names(dfv2) = cols
  
  # VISITS (based on SIZE) max
  museum_df$size_visits_max = sapply(museum_df$size, get_visits_by_museum_size_max)
  dfv <- plyr::ddply(museum_df, .(UNIT_GSS), summarise, n_museums=sum(size_visits_max), .drop=FALSE)
  dfv$VAR = "size_visits_max"
  dfv$visit_aggr = "sum.max"
  dfv3 = dfv[,c("UNIT_GSS","VAR","visit_aggr","n_museums")]
  names(dfv3) = cols
  dfv = rbind(dfv1,dfv2,dfv3,dfsp)
  
  # DECADE
  museum_df$opened_decade_est = ifelse( museum_df$year_opened_MEAN<1900,
                                        "1900before",
                                        paste0(substr(as.character(museum_df$year_opened_MEAN),0,3),"0s") )
  #print(summary(as.factor(museum_df$opened_decade_est)))
  df3 <- plyr::ddply(museum_df, .(UNIT_GSS,opened_decade_est), summarise, n_museums=length(project_id), .drop=FALSE)
  df3$VAR = "opened_decade_est"
  df3 = df3[,c("UNIT_GSS","VAR","opened_decade_est","n_museums")]
  names(df3) = cols
  
  # CLASS18
  #print(summary(as.factor(museum_df$opened_decade_est)))
  df4 <- plyr::ddply(museum_df, .(UNIT_GSS,subject_matter_simpl_aggr), summarise, n_museums=length(project_id), .drop=FALSE)
  df4$VAR = "subject_matter_simpl_aggr"
  df4 = df4[,c("UNIT_GSS","VAR","subject_matter_simpl_aggr","n_museums")]
  names(df4) = cols
  
  df = rbind(df,df1,df2,df3,df4,dfv)
  df$UNIT_GSS = as.factor(as.character(df$UNIT_GSS))
  df$VAR = as.factor(df$VAR)
  df$VAL = as.factor(df$VAL)
  return(df)
}

#  binning = "fixed", "sd", "equal",spTransform "pretty", "quantile", "kmeans", "hclust", "bclust", "fisher", "jenks" or "dpih"
FUN_bin_var = function( x, binning, bin_n, roundInt = F, isolateZeros = F ){
  if (isolateZeros){
    bin_n = bin_n-1
    x = x[x>0]
    cl = classIntervals( x, style = binning, n = bin_n)
    print(cl)
    brk = c(0,cl$brks)
  } else {
    cl = classIntervals( x, style = binning, n = bin_n)
    print(cl)
    brk = cl$brks
  }
  if (roundInt) brk = round(brk)
  stopifnot(is_unique(brk)) # 'Repeated values, please choose a different binning'
  return(brk)
}

# Based on 
# http://rstudio-pubs-static.s3.amazonaws.com/6577_3b66f8d8f4984fb2807e91224defa854.html
# it changes var_name in old_geog into new_geog, using overlaps
FUN_apportion_geographies = function(old_geog, new_geog, var_name){
  print('FUN_apportion_geographies')
  print('Calc intersections...')
  overlap.geog <- gIntersection(new_geog, spTransform(old_geog,proj4string(new_geog)), 
                                byid = T)
  overlap.area <- gArea(overlap.geog, byid = T)
  
  n.overlaps <- length(overlap.area)
  n.overlaps
  
  # NOTE: the type and order of calls to R functions below ensures that
  # resulting dataframe retains IDs in character format; and the area of
  # each overlap in numeric format
  tmp <- as.character(scan(text = names(overlap.area)))
  tmp <- matrix(tmp, nrow = n.overlaps, ncol = 2, byrow = TRUE, dimnames = list(c(NULL), 
                                                                                c("new.id", "old.id")))
  tmp <- as.data.frame(tmp, stringsAsFactors = FALSE)
  overlaps <- cbind(tmp, overlap.area)
  rownames(overlaps) <- NULL
  
  tmp <- gArea(old_geog, byid = TRUE)
  old.area <- data.frame(old.id = names(tmp), old.area = tmp, stringsAsFactors = FALSE)
  head(old.area, 1)
  
  overlaps <- data.frame(overlaps, old.area = old.area[match(overlaps[, "old.id"], 
                                                             old.area[, "old.id"]), "old.area"])
  head(overlaps)
  
  overlaps$old.area.share <- overlaps$overlap.area/overlaps$old.area
  
  # Find sum of old.area.share for each old.id
  tmp <- aggregate(old.area.share ~ old.id, data = overlaps, sum)
  colnames(tmp) <- c("old.id", "sum.old.area.shares")
  # Sort into same order as spatial units are stored in old.geog@data
  old.area.shares <- tmp[match(rownames(old_geog@data), tmp[, "old.id"]), ]
  old.area.shares
  
  # Add the sum.old.area.shares specific to each *old.id* to the relevant
  # overlap areas in the overlaps dataframe
  overlaps <- data.frame(overlaps, sum.old.area.shares = old.area.shares[match(overlaps[, 
                                                                                        "old.id"], old.area.shares[, "old.id"]), "sum.old.area.shares"])
  
  # Rescale the old.area.share for each overlap by 1/sum.old.area.shares
  overlaps$old.area.share.r <- overlaps$old.area.share * (1/overlaps$sum.old.area.shares)
  
  # Show that this rescaling has worked (i.e. that sum.old.area.shares.r = 1
  # for all old.ids)
  aggregate(old.area.share.r ~ old.id, data = overlaps, sum)
  
  new_geog@data[,var_name] <- 0
  for (i in 1:n.overlaps) {
    new_geog@data[overlaps$new.id[i], var_name] <- new_geog@data[overlaps$new.id[i], 
                                                                 var_name] + (overlaps$old.area.share.r[i] * old_geog@data[overlaps$old.id[i], 
                                                                                                                           var_name])
  }
  
  return(new_geog)
}


# spatial join
FUN_calc_grid_counts = function(points_sdf, grid_sdf, checkSum=T){
  verifyclass(points_sdf,'SpatialPointsDataFrame')
  verifyclass(grid_sdf,'SpatialPolygonsDataFrame')
  grid_sdf$CELL_ID_TMP = seq(nrow(grid_sdf))
  over_df = sp::over(points_sdf, spTransform(grid_sdf,proj4string(points_sdf)))
  counts = as.data.frame(table(over_df$CELL_ID_TMP))
  names(counts)=c('CELL_ID_TMP','VAR')
  #View(counts)
  gridc_sdf = merge(grid_sdf,counts,by='CELL_ID_TMP',all.x=T)
  gridc_sdf$VAR[ is.na(gridc_sdf$VAR) ] = 0
  gridc_sdf$CELL_ID_TMP = NULL
  if (checkSum) stopifnot(nrow(points_sdf)==sum(gridc_sdf$VAR))
  return(gridc_sdf)
}

# auto correl
# https://pro.arcgis.com/en/pro-app/tool-reference/spatial-statistics/h-how-spatial-autocorrelation-moran-s-i-spatial-st.htm
# The p-value is statistically significant, and the z-score is positive.
# You may reject the null hypothesis. The spatial distribution of high values 
# and/or low values in the dataset is more spatially clustered than would be 
# expected if underlying spatial processes were random.
FUN_globalMoransI = function(shapes, var_name, neiK, removeMissing = F){
  # find neighbours
  print(paste("Moran's I: VAR =",var_name,'; neiK =',neiK))
  print(nrow(shapes))
  if (removeMissing){
    shapes = subset( shapes, !is.na(shapes@data[,var_name] ))
    print(paste('after removed',nrow(shapes)))
  }
  mgrid_nei <- spdep::knearneigh( gCentroid(shapes,byid = T), k=neiK, RANN=F )
  nei = spdep::knn2nb(mgrid_nei)
  weights = nb2listw(nei)
  #plot(nei, coordinates(mgrid))
  mt = moran.test(shapes@data[,var_name], weights) 
  print(mt)
  print(round(mt$estimate,3))
  return(mt)
}
