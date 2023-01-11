
# compare random effect posterior levels  ---------------------------------

RanefPostComp<- function(df, comps, taxa, var){
  require(tidyverse)
  require(MCMCglmm)
  not_all_na <- function(x) any(!is.na(x))
  df$group <- df[, var]
  df$taxa <- df[, taxa]
  df<- df %>% filter(term!="(Intercept)") %>% 
    mutate(taxa=as.factor(taxa))
  LevelComps<-list()
  Levels=levels(df[, "taxa"])
  for(x in 1:length(Levels)){
    post1 <- df %>% subset(taxa==Levels[x] & group== comps[1]) %>% pull(posterior)
    post2 <- df %>% subset(taxa==Levels[x] & group== comps[2]) %>% pull(posterior)
    if(length(post2>0 ) & length(post1>0)){
      postComp <- post2-post1 %>% as.matrix() 
    }
    else{
      postComp <- rep(NA, 1000)
    }
    LevelComps[[x]] <-postComp
  }
  LevelCompsFin<- bind_cols(LevelComps) %>% as.matrix()
  colnames(LevelCompsFin) <- paste(Levels, comps[[2]], sep=".")
  LevelCompsFin %>% as.data.frame() %>% select_if(not_all_na) -> df 
  
  df1<-apply(df,2, function(a) mean(a)[1])
  df1l<-apply(df,2,function(a)  HPDinterval(as.mcmc(a))[1])
  df1u<-apply(df,2,function(a)  HPDinterval(as.mcmc(a))[2])
  df1v<-apply(df,2,function(a)  var(a)[1])
  
  df2 <- cbind(df1, df1l, df1u, df1v)
  colnames(df2) <- c("mean", "lower", "upper", "var")
  df2 <- df2 %>% as.data.frame() %>% 
    rownames_to_column(var="term") %>% 
    mutate(taxa=sapply(str_split(term, "[[.]]"), '[', 1),
           effect=var, 
           sig= ifelse(lower<0&upper<0, "*", 
                       ifelse(lower>0&upper>0, "*", ""))) 
  colnames(df2)[which(names(df2) == "taxa")] <- paste(taxa)
  return(df2)
  
}






TopAndBottom <- function(df, taxa=NULL,  threshold=NULL) {
  require(tidyverse) 
  
}
