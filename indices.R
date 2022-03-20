
require(tuneR);
require(soundecology);
require(seewave);
require(doParallel);

source("utils.R");

MAX_FREQ_ADI_AEI <- 10000;
MAX_FREQ_BIO     <- 8000;
MAX_FREQ_NDSI    <- 11000;

PROGRESS.LIST <- list(start = "", end = "", end.ar = "", current.name = "", current.file = 0, total.file = 0, current.slice = 0, total.slice = 0, path = "")

LIST.ENV.FUNC <- c("MAX_FREQ_ADI_AEI", "MAX_FREQ_BIO", "MAX_FREQ_NDSI", 
                   "get.index", "round.file.sec", "norm.ACI", "coalesce", "get.ndsi.values", 
                   "get.power.based", "psd", "msp", "spl")
LIST.ENV.PKG  <- c("soundecology", "seewave", "tuneR")

start.progress1 <- function(path, total.file)
{
  PROGRESS.LIST$total.file <<- total.file
  PROGRESS.LIST$path <<- path
  PROGRESS.LIST$start <<- Sys.time()        
}

inc.progress1 <- function(index, name)
{
  PROGRESS.LIST$current.file <<- index
  PROGRESS.LIST$current.name <<- tools::file_path_sans_ext(basename(name)) 
  
  progress()
}

start.progress2 <- function(total)
{
  PROGRESS.LIST$total.slice <<- total
}

inc.progress2 <- function(index)
{
  PROGRESS.LIST$current.slice <<- index
  
  progress()      
}

reset.progress2 <- function()
{
  PROGRESS.LIST$current.slice <<- 0
  PROGRESS.LIST$total.slice <<- 0      
}

progress <- function(end = FALSE, end.ar = FALSE)
{
  if(end)
    PROGRESS.LIST$end = Sys.time()
  if(end.ar)
    PROGRESS.LIST$end.ar = Sys.time()
  
  cat("Start: ", format(PROGRESS.LIST$start), "\n",
      "End   : ", format(PROGRESS.LIST$end), "\n",
      "End.AR: ", format(PROGRESS.LIST$end.ar), "\n\n",
      sep = "", 
      file = PROGRESS.LIST$path);

  if(end)
    cat("End...", file = PROGRESS.LIST$path, append = TRUE)
  else
  {
    cat(PROGRESS.LIST$current.file, "of", PROGRESS.LIST$total.file, "->", 
        round(PROGRESS.LIST$current.file / PROGRESS.LIST$total.file * 100, 2), "%\n",
        "File:",  PROGRESS.LIST$current.name,
        sep = " ", 
        file = PROGRESS.LIST$path, 
        append = TRUE);
    
    if(PROGRESS.LIST$total.slice > 1)
      cat("\n Slice  :",
          PROGRESS.LIST$current.slice, "of", PROGRESS.LIST$total.slice, "->", 
          round(PROGRESS.LIST$current.slice / PROGRESS.LIST$total.slice * 100, 2), "%\n", 
          sep = " ", 
          file = PROGRESS.LIST$path, 
          append = TRUE);       
  }
}

#slice.size in seconds for each audio's slice
slice.audio <- function(audio, slice.size)
{
  slice.list  = list(audio);
  #Similar to seewave::cutw  
  slice.size2 = round(slice.size * audio@samp.rate)
  qt.slice    = round(length(audio@left) / slice.size2)

  if(coalesce(slice.size) > 0 && qt.slice > 1 && round.file.sec(audio) != slice.size)
  {
    slice.list  = list();
    
    for(i in 1:qt.slice)
    {
      first.index = floor(1 + (slice.size2 * (i - 1)))
      last.index  = floor(ifelse(i == qt.slice, length(audio@left), i * slice.size2));
      
      slice.left  = audio@left[first.index:last.index];
      
      if(audio@stereo)
      {
        slice.right = audio@right[first.index:last.index];
        slice = Wave(left = as.numeric(slice.left), right = as.numeric(slice.right), samp.rate = audio@samp.rate, bit = audio@bit);
      }else
        slice = Wave(left = as.numeric(slice.left), samp.rate = audio@samp.rate, bit = audio@bit);        
      
      slice.list = append(slice.list, list(slice));
    }
  }
  
  return(slice.list);
}

round.file.sec <- function(file)
{
  file.length = length(file) / file@samp.rate;
  int.part    = trunc(file.length);
  float.part  = file.length - int.part;
  
  if(int.part > 0 && float.part >= 0.001)
    file.length = round(file.length, 3)
  else
    file.length = int.part;
  
  return(file.length);
}

#See soundecology documentation
norm.ACI <- function(file, value)
{
  time = round.file.sec(file) / 60;
  
  if(time > 1)
    value = value / time;
  
  return(coalesce(value));
}

psd <- function(audio, wn.param, wl.param, ovl.param)
{
  cat("Calculating PSD. Please wait... \n")  
  value  = audio
  window = ftwindow(wl = wl.param, wn = wn.param)
  ovl.param = floor(ovl.param/100 * wl.param)
  
  if(audio@stereo)
    value = mono(audio, "both")
  
  value = oce::pwelch(value@left, noverlap = ovl.param, window = window, fs = value@samp.rate, plot = FALSE)
  
  return(value$spec);
}

#Terrestrial landscapes ref = 10^-6
#Underwater landscapes  ref = 20 * 10^-6
msp <- function(audio, ref, psd.value = NULL, wn.param = "hanning", wl.param = 1024, ovl.param = 10)
{
  cat("Calculating MSP. Please wait... \n")  
  value = psd.value;
  
  if(is.null(psd.value))
    value = psd(audio = audio, wn.param = wn.param, wl.param = wl.param, ovl.param = ovl.param);
  
  value = 10 * log10(value)
  value = convSPL(value , pref = ref)$p
  value = sum(value^2) / length(value)
  
  return(value)
}

#Terrestrial landscapes ref = 10^-6
#Underwater landscapes  ref = 20 * 10^-6
spl <- function(audio, ref, msp.value = NULL, wn.param = "hanning", wl.param = 1024, ovl.param = 10)
{
  cat("Calculating SPL. Please wait... \n")  
  if(is.null(msp.value))
    msp.value = msp(audio = audio, ref = ref, wn.param = wn.param, wl.param = wl.param, ovl.param =  ovl.param)

  if(msp.value == 0)
    return(msp.value)
  else
    return(10 * log10( msp.value / ref^2 ))
}

#Different from the original index (https://eprints.qut.edu.au/110634/), it implements https://arxiv.org/abs/2201.02099
bgn <- function(env, smooth.window = 5, min.threshold = -90)
{
  #It uses Hilbert envelop as input  
  env = 20 * log10(env)
  env[env < min.threshold] = min.threshold
  
  #It applies thresholding of the histogram to return a suitable value
  return(otsu.threshold(env, smooth.hist.window = smooth.window))
}

get.melfcc <- function(file, return.list = TRUE)
{
  cat("Calculating MFCC. Please wait... \n")  
  #It does not process stereo recordings. If recording is stereo, run melfcc for each channel and
  #return the mean of the channels.
  M_left = Wave(left = as.numeric(file@left), samp.rate = file@samp.rate, bit = file@bit, pcm = file@pcm);
  is.stereo = file@stereo;
  
  #==========================================================================#
  #                           TUNER
  #==========================================================================#    
  M_left    = melfcc(M_left); 
  gc()
  M_right   = NULL;
  
  if(is.stereo)
  {
    M_right = Wave(left = as.numeric(file@right), samp.rate = file@samp.rate, bit = file@bit);
    M_right = melfcc(M_right);
    gc()
  }  
  
  rm(file);  
  
  if(return.list)
  {
    values = list();
    list_M = colMeans(M_left, na.rm = TRUE);
    
    if(!is.null(M_right))
      list_M = colMeans(rbind(list_M, colMeans(M_right)), na.rm = TRUE);  
    
    rm(M_left, M_right)    
    
    for(i in 1:length(list_M))
      values[[paste("MFCC_", i, sep = "")]] = list_M[i];
    
    return(values);    
  }
  else if(is.stereo)
    return((M_left + M_right) / 2)
  else 
    return(M_left)    
}

#Copy of seewave::AR. It uses M and Ht generated by other functions to improve processing time.
get.ar <- function(data)
{
  if(!is.null(data$M) && !is.null(data$Ht))
  {
    cat("Calculating AR. Please wait... \n")  
    ar.value = (rank(data$M) * rank(data$Ht)) / nrow(data)^2
    data = cbind(data, AR = ar.value)
  }
  
  return(data)
}

#Copy of seewave::H. It improves memory usage and return all data generated to be used for other indices.
my.H <- function(wave, f = wave@samp.rate, channel = 1, wl = 1024, envt = "hil", msmooth = NULL, ksmooth = NULL, wn = "hanning", ovl = 10)
{
  cat("(my) Calculating H. Please wait... \n")  
  input <- inputw(wave = wave, f = f, channel = channel)
  wave <- input$w
  f <- input$f
  rm(input)
  spec <- meanspec(wave = wave, f = f, wl = wl, plot = FALSE, wn = wn, ovlp = ovl)
  gc() #added
  SH <- sh(spec)
  enve <- env(wave = wave, f = f, envt = envt, msmooth = msmooth,
              ksmooth = ksmooth, plot = FALSE)
  TH <- th(enve)
  z <- SH * TH
  return(list(h = z, sh = SH, th = TH, env = enve, spec = spec)) #extended results
}

get.ndsi.values <- function(file, wl.param, return.left.right)
{
  value = ndsi(file, fft_w = wl.param, bio_max = min(MAX_FREQ_NDSI, file@samp.rate / 2))
  
  if(return.left.right)
    return(list(ANTHROPHONY_LEFT  = coalesce(value$anthrophony_left), ANTHROPHONY_RIGHT = coalesce(value$anthrophony_right),
                BIOPHONY_LEFT  = coalesce(value$biophony_left), BIOPHONY_RIGHT = coalesce(value$biophony_right),
                NDSI_LEFT = coalesce(value$ndsi_left), NDSI_RIGHT = coalesce(value$ndsi_right)))
  else
    return(list(ANTHROPHONY = mean(coalesce(value$anthrophony_left), coalesce(value$anthrophony_right)),
                BIOPHONY    = mean(coalesce(value$biophony_left), coalesce(value$biophony_right)),
                NDSI        = mean(coalesce(value$ndsi_left), coalesce(value$ndsi_right))))    
}

get.entropy.based <- function(file, wn.param, wl.param, ovl.param)
{
  value.aux = my.H(file, wn = wn.param, wl = wl.param, ovl = ovl.param)
  value     = coalesce(value.aux$h)
  
  cat("Calculating RMS. Please wait... \n")    
  value.rms = rms(value.aux$env);
  
  cat("Calculating Roughness. Please wait... \n")          
  value.roughness = roughness(value.aux$spec[, 2])
  
  cat("Calculating Number of Peaks. Please wait... \n")    
  value.npeaks = 0
  
  if(!any(is.na(value.aux$spec[, 2])))
    value.npeaks = nrow(fpeaks(value.aux$spec, f = file@samp.rate, plot = FALSE))
  
  cat("Calculating M. Please wait... \n")  
  #Code similar to seewave::M function but it uses the env calculated by H function.
  value.m  = median(value.aux$env) * (2^(1 - file@bit))
  
  cat("Calculating BGN. Please wait... \n")
  value.bgn = bgn(value.aux$env)


  return(list(H = value, Hf = value.aux$sh, Ht = value.aux$th, M = value.m, RMS = value.rms, 
              ROUGHNESS = value.roughness, NPEAKS = value.npeaks, BGN = value.bgn))
}

get.power.based <- function(file, wn.param, wl.param, ovl.param, aquatic)
{
  value.aux = psd(file, wn.param, wl.param, ovl.param);
  ref = ifelse(aquatic, 10^-6, 20 * 10^-6);
  
  value.msp = msp(file, ref, value.aux)
  value.spl = spl(file, ref, value.msp)
  
  cat("Calculating SNR Please wait... \n")
  mean.aux = mean(value.aux)
  std.aux  = sd(value.aux)

  return(list(MEAN_PSD = mean.aux, STD_PSD = std.aux, MSP = value.msp, SPL = value.spl, 
              SNR = mean.aux/std.aux))
}

get.index <- function(index, file, wn.param, wl.param, ovl.param, aquatic, return.left.right)
{
  #It does not increase processing time and memory consumption
  if(wl.param > 1024)
    wl.param = 1024
  if(ovl.param > 50)
    ovl.param = 10
  
  value.list  = list()
  value.left  = NULL
  value.right = NULL
  
  if(index == "ACI")
  {
      #default j = 5
      #default max_freq = file@samp.rate / 2
      value = acoustic_complexity(file, j = ifelse(round.file.sec(file) < 5, 1, 5), fft_w = wl.param)
    
      value.left  = norm.ACI(file, value$AciTotAll_left)
      value.right = norm.ACI(file, value$AciTotAll_right)    
  }else if(index == "ADI")
  {
    value = acoustic_diversity(file, max_freq = min(MAX_FREQ_ADI_AEI, file@samp.rate / 2))
    
    value.left  = coalesce(value$adi_left)
    value.right = coalesce(value$adi_right)
  }else if(index == "AEI")
  {
    value = acoustic_evenness(file, max_freq = min(MAX_FREQ_ADI_AEI, file@samp.rate / 2))
    
    value.left  = coalesce(value$aei_left)
    value.right = coalesce(value$aei_right)
  }else if(index == "BIO")
  {
    value = bioacoustic_index(file, max_freq = min(MAX_FREQ_BIO, file@samp.rate / 2))
    
    value.left  = coalesce(value$left_area)
    value.right = coalesce(value$right_area)
  }else if(index == "RUGOSITY")
  {    
    cat("Calculating Rugosity. Please wait... \n")      
    value.left  = rugo(file@left / max(file@left))
    value.right = ifelse(file@stereo, rugo(file@right / max(file@right)), 0)  
  }else if(index == "ZCR")
  { 
    cat("Calculating ZCR. Please wait... \n")
    value.left  = zcr(file, file@samp.rate, channel = 1, wl = NULL, plot = FALSE)
    value.right = zcr(file, file@samp.rate, channel = 2, wl = NULL, plot = FALSE)    
  }else if(index == "NDSI")
    value.list = get.ndsi.values(file, wl.param, return.left.right)
  else if(index == "H")
    value.list = get.entropy.based(file, wn.param, wl.param, ovl.param)
  else if(index == "PSD")
    value.list = get.power.based(file, wn.param, wl.param, ovl.param, aquatic)  
  else if(index == "MFCC")
    value.list = get.melfcc(file)
  
  if(index %in% c("ACI", "ADI", "AEI", "BIO", "RUGOSITY", "ZCR"))
  {
    if(return.left.right)
    {
      value.list[[paste(index, "_LEFT", sep = "")]]  = value.left
      value.list[[paste(index, "_RIGHT", sep = "")]] = value.right  
    }
    else
      value.list[[index]] = mean(c(value.left, value.right), na.rm = TRUE)
  }
  
  rm(file, value.left, value.right)
  gc()
  
  return(value.list)
}

generate.indices <- function(file, wn.param = "hanning", wl.param = 1024, ovl.param = 10, aquatic = FALSE, return.left.right = TRUE,
                            return.mfcc = TRUE)
{
  list.indices =   c("ACI", "ADI", "AEI", "BIO", "NDSI",      #from soundecology package
                     "RUGOSITY",                              #from seewave package
                     "PSD"                                    #MSP, SPL, and SNR are calculated with PSD values
  )

  cat("Calculating", list.indices, ".Please wait... \n")  

  `%myinfix%` <- ifelse(TRUE, `%dopar%`, `%do%`)

  #Indices that take long period to process
  values = foreach(i = 1:length(list.indices), .combine = append, .export = LIST.ENV.FUNC, .packages = LIST.ENV.PKG)%myinfix%
  {
    get.index(list.indices[i], file, wn.param, wl.param, ovl.param, aquatic, return.left.right)
  }

  #from seewave package. Ht, Hf are calculated inside H routine and RMS, ROUGHNESS, NPEAKS, M, and BGN used H routine values.
  values = append(values, get.index("H", file, wn.param, wl.param, ovl.param, aquatic, return.left.right))
  #from seewave package
  values = append(values, get.index("ZCR", file, wn.param, wl.param, ovl.param, aquatic, return.left.right))
  #from tuneR package

  #from tuneR
  if(return.mfcc)
    values = append(values, get.index("MFCC", file, wl.param, ovl.param, aquatic, return.left.right))
  
  return(values)
}

generate.spec <- function(file, file.name, target.path, wn.param = "hanning", wl.param = 1024, ovl.param = 10,
                         img.dim = c(1366, 768), palette = spectro.colors)
{
  cat("Saving spectrogram. Please wait... \n");
  img_name = file.path(target.path, paste(file.name, ".png", sep = ""));
  
  if(!all(file@left == 0))
  {
    png(img_name, width = img.dim[1], height = img.dim[2]);
    spectro(file, wl = wl.param, ovlp = ovl.param, wn = wn.param, f = file@samp.rate, palette = palette, collevels = seq(-80, 0, 5), grid = TRUE);
    dev.off();
  }  
}

process.file <- function(path, spec.path = NULL, aquatic = FALSE, generate = "index", 
                        slice.size = 1, img.dim = c(1366, 768), palette = spectro.colors, 
                        start.parallel = FALSE)
{
  indices = NULL;  
  cl = NULL;
  
  if(start.parallel)
  {
    qt.cores = detectCores() - 1
    cl = makeCluster(qt.cores, "PSOCK")
    registerDoParallel(cl)    
  }
  
  if(file.exists(path))
  {
    name      = basename(path);
    file.type = tools::file_ext(name); 
    name      = tools::file_path_sans_ext(name); 
    file      = NULL
    log.path  = file.path(dirname(path), "empty_files.txt")
    
    # Loading recording. It converts .flac to .wav if necessary
    if(file.type == "flac")
      wav2flac(path, reverse = TRUE)
    
    if(file.type == "mp3")
      file = readMP3(file.path(dirname(path), paste(name, ".mp3", sep = "")))
    else
      file = readWave(file.path(dirname(path), paste(name, ".wav", sep = "")))    
    
    if(file.type == "flac")
      file.remove(file.path(dirname(path), paste(name, ".wav", sep = "")))
    
    if(round.file.sec(file) == 0)
      cat(name, sep = "\n", file = log.path, append = TRUE)
    else
    {
      slices = slice.audio(file, slice.size);                    
      suffix.pattern = paste("%0", nchar(length(slices)), "d", sep = "");
      file.suffix    = 0;       
      
      if(round.file.sec(slices[[1]]) > 90)
        stop("Split large recordings at most 1 1/2 minute, because routines can have high memory consumption and take for a long period.")
      
      rm(file)
  
      start.progress2(length(slices))
      
      for(i in 1:length(slices))
      {
        name.suffix = "";  
        
        if(length(slices) > 1)
        {
          inc.progress2(i)
          file.suffix = file.suffix + 1;
          name.suffix = paste("_", sprintf(suffix.pattern, file.suffix), sep = "");  
        }
        
        # Generate only indices
        if(generate == "index")
        {
          new.row = generate.indices(file = slices[[i]], aquatic = aquatic, return.left.right = FALSE);
          new.row = append(list(filename = paste(name, name.suffix, ".", file.type, sep = "")),
                           new.row)
        }  
        # Generate only spectrograms
        if(generate == "spec")
        {
          generate.spec(slices[[i]], paste(basename(name), name.suffix, sep = ""), spec.path, img.dim = img.dim, palette = palette);
          new.row = list(filename = paste(name, name.suffix, ".", file.type, sep = ""))
        }  
        
        if(is.null(indices))
          indices = data.frame(new.row, stringsAsFactors = FALSE)
        else
          indices = rbind(indices, new.row);   
        
        rm(new.row);      
      }
      
      reset.progress2()
    }
  }else
    warning(paste("File [", path, "] was not found.", sep = ""));      
  
  if(start.parallel)
  {
    stopCluster(cl)
    registerDoSEQ()   
    indices = get.ar(indices)
  }
  
  return(indices);
}

process.dir <- function(source.path, target.path = source.path, aquatic = FALSE, generate = "index", 
                       slice.size = 1, batch.size = 100, img.dim = c(1366, 768), palette = spectro.colors)
{
  if(!generate %in% c("index", "spec"))
    stop("Do you want to generate indices or spectrogramas?") 
  if(dir.exists(source.path))
  {
    list.item = list.files(source.path, full.names = TRUE, pattern = "(.flac|.wav|.mp3)");
    
    if(length(list.item) > 0)
    {
      if(is.null(target.path) || is.na(target.path))
        target.path = source.path;      
      
      spec.path = file.path(target.path, "spec"); 
      log.path  = file.path(source.path, "indices_log.txt");
      log2.path = file.path(source.path, paste("log", "_", generate, ".txt", sep = ""));
      indices.path = file.path(target.path, paste(basename(source.path), ".csv", sep = ""))
      spec.time.path = file.path(target.path, paste(basename(source.path), "_spec_time.csv", sep = ""))  
      start.index = 1;
      
      if(generate == "spec" && !dir.exists(spec.path))
        dir.create(spec.path);
      if(file.exists(log.path) && file.info(log.path)$size > 0)
      {
        log = read.table(log.path, stringsAsFactors = FALSE, sep = "\n")
        
        if(nrow(log) > 0 && ncol(log) > 0)
          start.index = grep(log[1, 1], list.item)  
        
        log = log[0, ]
        write.table(log, log.path, quote = FALSE, col.names = FALSE);
      }
      
      indices = NULL
      print.col.names = TRUE
      
      qt.cores = detectCores() - 1
      cl = makeCluster(qt.cores, "PSOCK")
      registerDoParallel(cl)
      
      if(start.index > 1)
        cat("Restarting from ", list.item[start.index], "\n")
      
      start.progress1(log2.path, length(list.item))

      for(i in start.index:length(list.item))
      {
        inc.progress1(i, list.item[i])

        if(batch.size < 2 || i %% batch.size == 1)
          cat(list.item[i], sep = "\n", file = log.path, append = FALSE);     
        
        new.row = process.file(path = list.item[i], spec.path = spec.path, 
                               aquatic = aquatic, generate = generate, 
                               slice.size = slice.size, img.dim = img.dim, palette = palette);
        
        if(is.null(indices))
          indices = data.frame(new.row, stringsAsFactors = FALSE)
        else
          indices = rbind(indices, new.row);  
        
        if(generate == "index" && (i %% batch.size == 0 || i == length(list.item)))
        {
          write.table(indices, indices.path, 
                      append = TRUE, quote = FALSE, dec = ".", sep = ",", row.names = FALSE, 
                      col.names =  (print.col.names && start.index == 1))  
          rm(indices)
          indices = NULL
          print.col.names = FALSE;
        }
      }
      
      progress(end = TRUE)

      stopCluster(cl)
      registerDoSEQ()
      
      if(file.exists(log.path) && file.info(log.path)$size > 0)
        cat("", file = log.path, append = FALSE);

      if(generate == "index")
      {
        merged.path = file.path(target.path, paste(basename(source.path), "_Indices.csv", sep = ""))
        
        indices = read.csv(indices.path)
        index   = ncol(indices)
        indices = get.ar(indices)

        write.table(indices, merged.path, quote = FALSE, dec = ".", sep = ",", row.names = FALSE)
        file.remove(indices.path)
        
        progress(end = TRUE, end.ar = TRUE)
        
        return(indices)
      }
    }else
      warning(paste("Directory [", source.path, "] is empty.", sep = ""));    
    
  }else
    warning(paste("Directory [", source.path, "] was not found.", sep = ""));    
}
