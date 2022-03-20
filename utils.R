
otsu.threshold <- function(x, qt.bins = 100, remove.outlier = FALSE, smooth.hist.window = NULL)
{
  if(remove.outlier)
    x = remove.outilier(x)
  
  brk = seq(min(x), max(x),  (max(x) - min(x)) / qt.bins )
  hist.value = hist(x, breaks = brk, plot = FALSE)
  
  mids = hist.value$mids
  hist = hist.value$counts
  
  if(!is.null(smooth.hist.window))
    hist = moving.avg(hist, smooth.hist.window)
  
  if(length(mids) > length(hist))
    mids = mids[-1]
  
  prop = hist
  
  if(sum(hist) != 1)
    prop = hist / sum(hist)
  
  omega = cumsum(prop)
  mu    = cumsum(mids * prop)
  muT   = sum(mids * prop)
  max.varB = -1
  maxk     = -1  
  
  for(k in 1:length(prop))
  {
    varB = (muT * omega[k] - mu[k])^2 / ifelse(omega[k] %in% c(0, 1), 1, (omega[k] * (1 - omega[k])))  
    
    if(varB > max.varB)
    {
      max.varB = varB
      maxk = k
    }
  }  
  
  return(mids[maxk])
}

remove.outilier <- function(x)
{
  while(TRUE)
  {
    qt  = quantile(x)
    q1  = qt[2]
    q3  = qt[4]
    iqr = q3 - q1
    
    begin_length = length(x)
    
    x = x[x >= q1 - iqr]
    x = x[x <= q3 + iqr]
    
    if(begin_length == length(x))
      break
  }
  
  return(x)  
}

coalesce <- function(value)
{
  if(is.null(value) || is.na(value))
    return(0)
  else 
    return(value)
}

moving.avg <- function(x, window = 3, padding.zero = FALSE, causal = FALSE)
{
  x.mean  = c()
  
  if(causal)
  {
    padding = window - 1
    x.new   = c(x[1:padding], x)
    
    if(padding.zero)
      x.new = c(rep(0, padding), x)
    
    for(i in (padding + 1):length(x.new))
      x.mean = c(x.mean, mean(x.new[(i - padding):i]))    
  }else
  {  
    padding = floor(window / 2)    
    x.new   = c(x[1:padding], x, x[(length(x) - padding + 1):length(x)])
    
    if(padding.zero)
      x.new = c(rep(0, padding), x, rep(0, padding))
    
    for(i in (padding + 1):(length(x.new) - padding))
      x.mean = c(x.mean, mean(x.new[(i - padding):(i + padding)]))
  }
  
  return(x.mean);
}
