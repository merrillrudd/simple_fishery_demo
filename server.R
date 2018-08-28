# Define server logic required to draw figure
library(shiny)

shinyServer(
  function(input, output) 
  {    

################################
## population dynamic functions
################################
    ## length at age
get_lengths <-function(Linf,k,t0,ages)
{
  Lengths_exp <- Linf*(1-exp(-k*(ages-t0)))
  return(Lengths_exp)
}

  ## selectivity at age
get_selex <- function(aselex, ages, dome){
  if(dome==1) dome_sd <- 15
  if(dome==2) dome_sd <- 8
  if(dome==3) dome_sd <- 2
  Selex_a <- 1/(1+exp(-log(19)*(ages-aselex)/((aselex+2)-aselex)))# Selectivity at age
  if(dome!=0){
    Sfull <- which(Selex_a>=0.99)[1]
    find_dome <- (Sfull+1):length(Selex_a)
    Selex_a[find_dome] <- exp((-(find_dome - Sfull)^2)/(2*dome_sd^2))
  }
  return(Selex_a)
}

  ## numbers at age
get_Na <- function(ages, aselex, M, fish, R0){
  S_a <- get_selex(aselex, ages, input$dome)
  N_a <- rep(NA, length(ages))
  N_a[1] <- R0
  for(i in 2:length(ages)){
    if(i<length(ages)) N_a[i] <- exp(-M-S_a[i-1]*fish)*N_a[i-1]
    if(i==length(ages)) N_a[i] <- (N_a[i-1]*exp(-M-fish*S_a[i-1]))/(1-exp(-M-fish*S_a[i]))
  }
  return(N_a)
}

  ## proportion mature at age
get_mature <- function(agemat, ages){
  Mat_a <- 1/(1+exp(-log(19)*(ages-agemat)/((agemat+2)-agemat)))# Selectivity at age
  return(Mat_a)
}


  ## fecundity at age
get_fec <- function(agemat, ages, Linf, k, t0, lwa, lwb){
  Mat_a <- get_mature(agemat, ages)
  L_a <- get_lengths(Linf, k, t0, ages)
  W_a <- lwa * L_a ^ lwb
  Fec_a <- Mat_a*W_a
  return(Fec_a)
}

  ## spawning biomass at age
get_SB <- function(ages, agemat, Linf, k, t0, lwa, lwb, aselex, M, fish, R0){
  N_a <- get_Na(ages, aselex, M, fish, R0)
  Fec_a <- get_fec(agemat, ages, Linf, k, t0, lwa, lwb)
  
  SB <- Fec_a*N_a
  return(SB)
}

  ## catch at age
get_catch <- function(ages, aselex, M, fish, R0, Linf,k,t0,lwa, lwb){
  N_a <- get_Na(ages, aselex, M, fish, R0)
  S_a <- get_selex(aselex, ages, input$dome)
  L_a <- get_lengths(Linf, k, t0, ages)
  W_a <- lwa * L_a ^ lwb
  Cb_a <- N_a*W_a*(1-exp(-M-S_a*fish))*fish*S_a/(M+S_a*fish)
  return(Cb_a)
}

  ## spawning biomass per recruit, either in the fished or unfished conditions
get_SBPR <- function(agemat, Linf, k, t0, lwa, lwb, aselex, M, fish, R0, unfished, text){
  ages <- c(0:input$amax)
  #SB <- get_SB(ages, agemat, Linf, k, t0, lwa, lwb, aselex, M, fish, R0)
  S_a <- get_selex(aselex, ages, input$dome)
  Mat_a <- get_mature(input$amat, ages)
  L_a <- get_lengths(Linf, k, t0, ages)
  W_a <- lwa * L_a ^ lwb
  
  Na0 <- Naf <- rep(NA, length(ages))
  Na0[1] <- Naf[1] <- 1
  for (a in 2:length(ages)) {
    if (a < length(ages)) {
      Na0[a] <- Na0[a - 1] * exp(-M)
      Naf[a] <- Naf[a - 1] * exp(-M - S_a[a - 1] * fish)
    }
    if (a == length(ages)) {
      Na0[a] <- (Na0[a - 1] * exp(-M))/(1 - exp(-M))
      Naf[a] <- (Naf[a - 1] * exp(-M - fish * S_a[a - 1]))/(1 - exp(-M - fish * S_a[a - 1]))
    }
  }
  SB0 <- sum(Na0 * Mat_a * W_a)
  SBf <- sum(Naf * Mat_a * W_a)
  
  SBPR <- SBf
  if(text==TRUE & unfished==FALSE) return(paste0("SBPR(fished) = ", round(SBPR,0)))
  if(text==TRUE & unfished==TRUE) return(paste0("SBPR(unfished) = ", round(SBPR,0)))
  
  if(text==FALSE) return(round(SBPR,0))
}

  ## spawning potential ratio
get_SPR <- function(agemat, Linf, k, t0, lwa, lwb, aselex, M, fish, R0, text=TRUE){
  ages <- c(0:input$amax)
#   SBf <- get_SB(ages, agemat, Linf, k, t0, lwa, lwb, aselex, M, fish, R0)
#   SBPRf <- sum(SBf)/R0
#   
#   SB0 <- get_SB(ages, agemat, Linf, k, t0, lwa, lwb, aselex, M, 0, R0)
#   SBPR0 <- sum(SB0)/R0

#SB <- get_SB(ages, agemat, Linf, k, t0, lwa, lwb, aselex, M, fish, R0)
S_a <- get_selex(aselex, ages, input$dome)
Mat_a <- get_mature(input$amat, ages)
L_a <- get_lengths(Linf, k, t0, ages)
W_a <- lwa * L_a ^ lwb

Na0 <- Naf <- rep(NA, length(ages))
Na0[1] <- Naf[1] <- 1
for (a in 2:length(ages)) {
  if (a < length(ages)) {
    Na0[a] <- Na0[a - 1] * exp(-M)
    Naf[a] <- Naf[a - 1] * exp(-M - S_a[a - 1] * fish)
  }
  if (a == length(ages)) {
    Na0[a] <- (Na0[a - 1] * exp(-M))/(1 - exp(-M))
    Naf[a] <- (Naf[a - 1] * exp(-M - fish * S_a[a - 1]))/(1 - exp(-M - fish * S_a[a - 1]))
  }
}
SB0 <- sum(Na0 * Mat_a * W_a)
SBf <- sum(Naf * Mat_a * W_a)

SPR <- SBf/SB0
  if(text==TRUE) return(paste0("Spawning Potential Ratio = ", round(SPR, 3)))
  if(text==FALSE) return(round(SPR, 3))
}

AgeToLengthComp <-
  function(L_a,
           S_a,
           Fish,
           ages,
           comp_sample,
           sample_type = 'catch') {
    ################################################
    ## Probability being in a length bin given age
    ################################################
    highs <- seq(1,max(L_a)*1.5,by=1)
    lows <- highs - 1
    
    lbprobs <-
      function(mnl, sdl)
        return(pnorm(highs, mnl, sdl) - pnorm(lows, mnl, sdl))
    vlprobs <- Vectorize(lbprobs, vectorize.args = c("mnl", "sdl"))
    plba <- t(vlprobs(L_a, L_a * 0.1))
    plba <- plba / rowSums(plba)
    
    ################################################
    ## Probability being in harvested at an age
    ################################################
    N_a <- get_Na(ages, input$aselex, input$M, Fish, 1)
      if(sample_type=="catch") page <- N_a * S_a
    if(sample_type!="catch") page <- N_a
    page <- page / sum(page)
    
    ################################################
    ## Probability of sampling a given length bin
    ################################################
      plb <- page %*% plba
    plb <- plb / rowSums(plb)
    
    #######################
    ## Length frequencies
    #######################
      if(is.na(sum(plb))==FALSE){
        LF <- rmultinom(n = 1,
                             size = comp_sample,
                             prob = plb)
      }
      if(is.na(sum(plb))){
        LF <- NA
      }
    
    Outs <- NULL
    Outs$plba <- plba
    Outs$plb <- plb
    Outs$page <- page
    Outs$LF <- LF
    return(Outs)
  }

  ## yield per recruit
get_YPR <- function(aselex, M, fish, Linf,k,t0,lwa, lwb, R0, text)
{
  ages <- c(0:input$amax)
  Cb_a <- get_catch(ages, input$aselex, input$M, input$fish, 1, input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79)
  YPR <- sum(Cb_a)/R0 
  if(text==TRUE) return(paste0("YPR = ", round(YPR,0)))
  if(text==FALSE) return(YPR)
}

find_Fref <- function(SBref, agemat, Linf, k, t0, lwa, lwb, aselex, M, fish, R0, unfished=FALSE, text=FALSE){
  SBcalc <- get_SBPR(input$amat, input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79, input$aselex, input$M, fish, 1, unfished=FALSE, text=FALSE)
  SBdiff <- SBcalc - SBref
  return(SBdiff)
}

################################
## output
################################

## length
output$VBGFplot <- renderPlot(
{
  ages<-c(0:input$amax)
  lengths.out<- get_lengths(input$Linf,input$k,t0=-0.01,ages)
  # plot VBGF
  par(mfrow=c(1,1),mar=c(6,6,2,2))
  plot(ages, lengths.out, col = "#0000AA",
       xlab="Age",ylab="Length (cm)",xlim=c(0,input$amax),
       ylim=c(0,input$Linf*1.1),type="l",lwd=5,
       main="Length at age", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, cex.main=2)
  
  lower <- lengths.out - (0.1*lengths.out)
  upper <- lengths.out + (0.1*lengths.out)
  polygon(x=c(ages, rev(ages)), y=c(lower,rev(upper)), col="#0000AA30", border=NA)
  abline(v=input$amat, col="red", lty=2, lwd=3)
  abline(v=input$aselex, col="black", lty=2, lwd=3)
  # legend("bottomright", legend=c("Age at 50% selectivity", "Age at 50% maturity"), col=c("black","red"), lty=c(1,2), lwd=3, cex=1.3)
}
    )

## numbers alive at age
output$NumbersAtAge <- renderPlot(
{
  ages <- c(0:input$amax)
  N_a <- get_Na(ages, input$aselex, input$M, input$fish, 1)
  N_a0 <- get_Na(ages, input$aselex, input$M, 0, 1)
  ## plot numbers at age
  par(mfrow=c(1,1),mar=c(6,6,2,2))
  plot(ages, N_a0, col = "#AAAAAA",
       xlab="Age",ylab="Numbers Alive",xlim=c(0,input$amax),
       ylim=c(0,1),type="l",lwd=5, main="Numbers at age", xaxs="i", yaxs="i",
       cex.axis=2,cex.lab=2,cex.main=2)
  if(input$fish>0) lines(ages, N_a, col="#00AA0080", lwd=5)
  # legend("topright", legend=c("Unfished", "Fished"), col=c("#AAAAAA", "#00AA0080"), lwd=5, cex=1.3)
})

output$CatchAtLength <- renderPlot({
  ages <- c(0:input$amax)
  L_a <- get_lengths(input$Linf, input$k, t0=-0.01, ages)
  S_a <- get_selex(input$aselex, ages, input$dome)
  set.seed(123)
  LF0 <- AgeToLengthComp(L_a=L_a, S_a=S_a, ages=ages, Fish=0, comp_sample=input$lsamp)
  set.seed(123)
  LF <- AgeToLengthComp(L_a=L_a, S_a=S_a, ages=ages, Fish=input$fish, comp_sample=input$lsamp)
  
  par(mfrow=c(1,1),mar=c(6,6,2,2))
  barplot(t(LF0$LF), col="#AAAAAA",xlim=c(0,nrow(LF0$LF)), ylim=c(0,max(LF0$LF)*1.5), 
          border=NA, xlab="Length (cm)", ylab="Frequency", main="Lengths in catch",
          xaxs="i",yaxs="i", cex.axis=2, cex.lab=2, cex.main=2)
  # legend("topleft", legend=c("Unfished","Fished","Length at 50% selectivity","Length at 50% maturity"), 
         # col=c("#AAAAAA", "#00AA0080", "black", "red"), bg='white', pch=c(15,15,NA,NA), lty=c(NA,NA,1,2), lwd=3, cex=1.3)
  if(input$fish>0){
    par(new=TRUE)
    barplot(t(LF$LF), col='#00AA0080', xlim=c(0,nrow(LF0$LF)), ylim=c(0,max(LF0$LF)*1.5), border=NA)
  }
  box()
  axis(1, cex.axis=2)
  SL50 <- input$Linf * (1- exp(-input$k*(input$aselex-(-0.01))))
  ML50 <- input$Linf * (1- exp(-input$k*(input$amat-(-0.01))))
  abline(v=SL50, lwd=3,lty=2)
  abline(v=ML50, col="red", lty=2, lwd=3)
  abline(v=input$Linf, col="#0000AA", lwd=3, lty=2)

})

output$plotLegend <- renderPlot({
  par(mfrow=c(1,1),mar=c(0,0,0,0))
  plot(x=1,y=1,type="n",xlim=c(1,10),ylim=c(1,10), axes=F, ann=F)
  legend("topleft", legend=c("Unfished","Fished","Selectivity","Length at \n50% selectivity","Maturity","Length at \n50% maturity", "Length-at-age", "Asymptotic length"), 
         col=c("#AAAAAA","#00AA0080","black","black","red","red", "#0000AA", "#0000AA"), lwd=3, lty=c(1,1,1,2,1,2,1,2), cex=1.3)
  
  ages <- c(0:input$amax)
  L_a <- get_lengths(input$Linf, input$k, t0=-0.01, ages)
  S_a <- get_selex(input$aselex, ages, input$dome)
  set.seed(123)
  LF <- AgeToLengthComp(L_a=L_a, S_a=S_a, ages=ages, Fish=input$fish, comp_sample=input$lsamp)
  LF_mat <- t(LF$LF)
  
  ML50 <- input$Linf * (1- exp(-input$k*(input$amat-(-0.01))))
  prop <- sum(LF_mat[1:ceiling(ML50)])/sum(LF_mat)
  
  
  if(input$fish==0)  text(x=5,y=3,"No fishing", cex=3)
  if(input$fish > 0) text(x=5,y=3, paste0(round(prop*100), "% of catch\nless than L50"), cex=3)

  
})

## weight at age
output$WeightAtAge <- renderPlot(
  {
    ages <- c(0:input$amax)
    L_a <- get_lengths(input$Linf, input$k, t0=-0.01, ages)
    W_a <- 0.0245 * L_a ^ 2.79
    plot(ages, W_a, col="navyblue", xlab="Age", ylab="Weight (g)",
         xlim=c(0, input$amax), ylim=c(0, max(W_a)*1.1), type="l", lwd=5,
         main="Weight At Age", xaxs="i", yaxs="i")
  })


## maturity at age
output$MatureAtAge <- renderPlot(
{
  ages <- c(0:input$amax)
  Mat_a <- get_mature(input$amat, ages)
  plot(ages, Mat_a, col="goldenrod",
       xlab="Age", ylab="Proportion Mature", xlim=c(0, input$amax),
       ylim=c(0, max(Mat_a)*1.1), type="l", lwd=5, main="Maturity at Age", xaxs="i", yaxs="i")
}
  )

## fecundity at age
output$FecundityAtAge <- renderPlot(
  {
    ages <- c(0:input$amax)
    Fec_a <- get_fec(input$amat, ages, input$Linf, input$k, t0=-0.01, lwa=0.0245, 2.79)
    plot(ages, Fec_a, col = "goldenrod4",
         xlab="Age", ylab="Fecundity", xlim=c(0, input$amax),
         ylim=c(0,max(Fec_a)*1.1), type="l", lwd=5, main="Fecundity at age", xaxs="i", yaxs="i")
  })


## selectivity at age
output$SelexAtAge <- renderPlot(
{
  ages <- c(0:input$amax)
  S_a <- get_selex(input$aselex, ages, input$dome)
  plot(ages, S_a, col = "forestgreen",
       xlab="Age", ylab="Selectivity", xlim=c(0, input$amax),
       ylim=c(0,max(S_a)*1.1), type="l", lwd=5, main="Selectivity at age", xaxs="i", yaxs="i")
})

## selectivity and maturity
output$SelexMature <- renderPlot(
  {
    ages <- c(0:input$amax)
    S_a <- get_selex(input$aselex, ages, input$dome)
    Mat_a <- get_mature(input$amat, ages)
    
    par(mfrow=c(1,1),mar=c(6,6,2,2))
    plot(ages, S_a, col="black", xlab="Age", ylab="Proportion vulnerable/mature", 
         xlim=c(0, input$amax), ylim=c(0, max(S_a)*1.1), type="l", lwd=4, main="Selectivity and Maturity", xaxs="i", yaxs="i",
         cex.axis=2, cex.lab=2, cex.main=2)
    lines(ages, Mat_a, col="red", lwd=4)
    # legend("bottomright", xpd=NA, legend=c("Selectivity", "Maturity"), lwd=4, lty=c(1,2), col=c("black", "red"), border=NA, cex=1.3)
  })


## spawning biomass at age
output$SpawnBioAtAge <- renderPlot(
{
  ages <- c(0:input$amax)
  SB <- get_SB(ages, input$amat, input$Linf, input$k, t0=-0.01, lwa=0.0245, 2.79, input$aselex, input$M, input$fish, 1)
  plot(ages, SB, col="tomato3", xlim=c(0, input$amax),
       ylim=c(0, max(SB)*1.1), type="l", lwd=5, main="Spawning Biomass at age", xaxs="i", yaxs="i")
}
  )

## catch at age
output$CatchAtAge <- renderPlot(
  {
    ages <- c(0:input$amax)
    Cb_a <- get_catch(ages, input$aselex, input$M, input$fish, 1, input$Linf, input$k, t0=-0.01, lwa=0.0245, 2.79)
    plot(ages, Cb_a, col="darkred", xlim=c(0, input$amax), 
         ylim=c(0, max(Cb_a)*1.1), type="l", lwd=5, main="Catch at age", xaxs="i", yaxs="i")
    
  })

## spawning biomass per recruit in fished condition
output$SBPRf <- renderText(
  {
    get_SBPR(input$amat, input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79, input$aselex, input$M, input$fish, 1, unfished=FALSE, text=TRUE)
  })


## spawning biomass per recruit in unfished condition
output$SBPR0 <- renderText(
  {
    get_SBPR(input$amat, input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79, input$aselex, input$M, 0, 1, unfished=TRUE, text=TRUE)
  })

## spawning potential ratio
output$SPR <- renderText(
  {
    get_SPR(input$amat, input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79, input$aselex, input$M, input$fish, 1)
  })


## yield per recruit
output$YPR <- renderText(
  {
    get_YPR(input$aselex, input$M, input$fish, input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79, 1, text=TRUE)
  }
)

output$YPRplot <- renderPlot(
  {
    uvec <- seq(0,2, by=0.02)
    YPRvec <- sapply(1:length(uvec), function(x) get_YPR(input$aselex, input$M, uvec[x], input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79, 1, text=FALSE))
    plot(uvec, YPRvec, pch=19, cex=1.5, col="darkred", xlim=c(min(uvec), max(uvec)), ylim=c(0, max(YPRvec)*1.1), 
         xlab="Fishing Mortality", ylab="Yield Per Recruit", xaxs="i", yaxs="i")
  })
# 
# output$SBPRplot <- renderPlot(
# {
#   uvec <- seq(0,2, by=0.02)
#   SBPRvec <- sapply(1:length(uvec), function(x) get_SBPR(input$amat, input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79, input$aselex, input$M, uvec[x], 1, unfished=FALSE, text=FALSE))
#   plot(uvec, SBPRvec, pch=19, cex=1.5, col="tomato3", xlim=c(min(uvec), max(uvec)), ylim=c(0, max(SBPRvec)*1.1),
#        xlab="Fishing Mortality", ylab="Spawning Biomass Per Recruit", xaxs="i", yaxs="i", main="Find reference point")
#   SBref <- (input$Fref/100)*SBPRvec[1] ## target percentage of spawning biomass from unfished state
#   Fval <- uniroot(find_Fref, lower=uvec[1], upper=50, SBref=SBref, agemat=input$amat, Linf=input$Linf, k=input$k, t0=-0.01, 
#                   lwa=0.0245, lwb=2.79, aselex=input$aselex, M=input$M, R0=1)$root  
#   segments(x0=0, x1=Fval, y0=SBref, y1=SBref, lty=2)
#   segments(x0=Fval, x1=Fval, y0=0, y1=SBref, lty=2)
#   points(x=Fval, y=SBref, pch=16, cex=2.5)
#   text(x=0.8, y=max(SBPRvec)*0.8, paste0("F", input$Fref, " = ", round(Fval, 2)), cex=2, font=2)
# 
# })
# 
# output$Kobe <- renderPlot(
#   {
#     plot(x=1, y=1, type="n", xlim=c(0, 1), ylim=c(0, 3), xaxs="i", yaxs="i",
#          xlab="Spawning Potential Ratio", ylab="F/Fref", main="Status relative to reference point")
#     abline(v=input$Fref/100, lty=2)
#     abline(h=1, lty=2)
#     SPRval <- get_SPR(input$amat, input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79, input$aselex, input$M, input$fish, 1, text=FALSE)
#     
#     uvec <- seq(0,2, by=0.02)
#     SBPRvec <- sapply(1:length(uvec), function(x) get_SBPR(input$amat, input$Linf, input$k, t0=-0.01, lwa=0.0245, lwb=2.79, input$aselex, input$M, uvec[x], 1, unfished=FALSE, text=FALSE))
#     SBref <- (input$Fref/100)*SBPRvec[1] ## target percentage of spawning biomass from unfished state
#     Fval <- uniroot(find_Fref, lower=uvec[1], upper=50, SBref=SBref, agemat=input$amat, Linf=input$Linf, k=input$k, t0=-0.01, 
#                     lwa=0.0245, lwb=2.79, aselex=input$aselex, M=input$M, R0=1)$root  
#     
#     polygon(x=c(0, input$Fref/100, input$Fref/100, 0), y=c(1, 1, 3, 3), col="#AA000050", border=NA)
#     polygon(x=c(input$Fref/100, 1, 1, input$Fref/100), y=c(0, 0, 1, 1), col="#00AA0050", border=NA)
#     
#     points(x=SPRval, y=input$fish/Fval, pch=19, col="blue", cex=2)
#   })


## end
}
)

