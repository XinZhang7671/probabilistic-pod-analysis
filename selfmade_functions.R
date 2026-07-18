## 
get_pf <- function(df_assay, wth = 0){
lims <- range(df_assay$m - 3*df_assay$s, 
              df_assay$m + 3*df_assay$s)  # ±3×SE

conc_vector <- seq(lims[1], lims[2], length.out = 400)

conam <- unique(df_assay$Compound)

df_pf <- do.call('rbind',lapply(1:length(conam),function(j){
  data.frame(conc=conc_vector, 
             p_conc = sapply(conc_vector, function(conc) {
    df_c <- df_assay %>%
      filter(Compound==conam[j]) %>%
      filter(w > wth)
    1 - prod(1 - pnorm(conc, mean = df_c$m, sd = df_c$s))
  }), Compound=conam[j], Assay = df_assay$Assay[1])
}))
return(df_pf)
}

## 
plot_pao <- function(cnam, df_assay, df_pf, wth = NULL, legend = TRUE, saveplot = TRUE){
  df_comp <- df_assay %>%
    filter(Compound==cnam) 
  
  if(!is.null(wth)){
    df_comp <- df_comp%>%
    filter(w>wth)
  }
  
  n <- nrow(df_comp)
  pp = ppoints(200)
  
  dd <- do.call('rbind',lapply(1:n,function(i){

    x <- qnorm(pp,df_comp$m[i],df_comp$s[i])
    df <- data.frame(pp=pp,x=x,
                     y = dnorm(x,df_comp$m[i],df_comp$s[i]),
                     Feature = df_comp$Features[i],
                     w = df_comp$w[i])
    return(df)
  }))
  
  df_pf_comp <- df_pf %>%
    filter(Compound==cnam) %>%
    mutate(Feature = "overall") %>%
    filter(p_conc > 0.001 & p_conc < 0.999) 
  
   denom <- max(dd$y)*2
   
 p <-  ggplot(data=dd,aes(x=x, y=pp, group=Feature, col=Feature,linewidth = w)) +
#  p <-  ggplot(data=dd,aes(x=x, y=y/denom, group=Feature, col=Feature)) +
   geom_line(data=df_pf_comp,aes(x=conc,y=p_conc),col = 'darkmagenta', linewidth = 1.5, alpha = 0.5) +
   geom_line() +
   labs(title="Overall probability of at least one adverse effect at different concentrations",
         subtitle = paste(df_comp$Assay[1],"for",cnam)) +
    xlab("concentration") +
    ylab("probability") +
    ylim(0,1) 
 
  
  if(!legend){
  p <-  p + theme(legend.position="none")
  }
 if(!is.null(wth)){
 p <- p + scale_linewidth(range = c(0, 1))
 }else{
   p <- p + scale_linewidth(range = c(0.2, 0.2))
 }
 if(saveplot){
 ggsave(paste0("img/PAO_",df_comp$Assay[1],"_",cnam,".jpg"))
 }
   p
}

## 
calc_pod <- function(pinvitro, df_assay, df_pf){
df_pod <- df_pf %>%
  group_by(Compound) %>%
  filter(p_conc < pinvitro) %>%
  slice_max(conc) %>%  
  mutate(pPoD = conc, pinvitro = pinvitro) %>%
  select(Compound, Assay, pinvitro, pPoD) 

df_alt <- df_assay %>%
  group_by(Compound) %>% 
  mutate(MIN = min(m)) %>%
  mutate(P05 = quantile(m, probs = 0.05)) %>%
  select(Compound, MIN, P05, Assay) %>%
  slice(1)

df_pod %>% left_join(df_alt, by = c("Compound","Assay"))
}

## 
plot_pod <- function(cnam,df_pod,df_pf) {
  
  df_pod_c <- df_pod %>%
    filter(Compound == cnam)
    
  p <- df_pf %>% filter(Compound == cnam) %>%
    filter(p_conc > 0.001 & p_conc < 0.999) %>%
    ggplot(aes(x = conc, y = p_conc)) +
    geom_line(color = "darkmagenta", linewidth = 1, alpha = 0.5) +
    
    geom_vline(aes(xintercept = pPoD), data = df_pod_c, color = "black", linetype = "solid") +
    geom_vline(aes(xintercept = P05), data = df_pod_c, color = "blue", linetype = "solid") +
    geom_vline(aes(xintercept = MIN), data = df_pod_c, color = "red", linetype = "dashed") +
    
    labs(
      title = paste("Probability at least on effect for", cnam, "given in vitro", df_pod$Assay),
      subtitle = paste0("pPOD_",(df_pod$pinvitro*100),", min and 5th percentile over features"), 
      x = "Concentration",
      y = "Probability"
    ) +
    theme_minimal()
  ggsave(paste0("img/PAO_compare_",df_pod_c$Assay,"_",cnam,".jpg"), p)
  
  p
}

