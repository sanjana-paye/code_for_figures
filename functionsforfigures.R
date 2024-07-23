figure1 <- function(d, a, b, t){
  par(mfrow = c(2, 2), mgp = c(1, 3, 10), mar = c(2, 2, 2, 2), cex.axis = 0.5, las = 1)
  ncpsurfaceplotfigure_risk(b, t)
  ncpsurfaceplotfigure_protective(b, t)
  
  
  # par(mfrow = c(1, 2), mgp = c(3, 1, 0), mar = c(5.1, 4.1, 2.1, 2.1))
  allele <- 'risk'
  
  par(mgp = c(2, 1, 0), mar = c(4, 4.5, 1, 2))
  riskvsprotectivefigure(d, a, b)
}

ncpsurfaceplotfigure_risk <- function(b, t){
  #surface plot figure
 
  # d: disease frequency values at fixed intervals
  # a: risk allele frequency values at fixed intervals
  # z: function that defines the surface
  d = seq(0, 0.99, by = 0.03)
  a = seq(0, 0.99, by = 0.03)
  #t = 1
  #b = 0.5
  z <- outer(d,a, function(d,a) (((b-d)^2)*a*t*(d*t-1))/(d*(1+b*(t-1)-d*t)*(-1+d+a+a*b*(t-1)-d*a*t)))

  z_new <- pmin(z, 1)
  
  p <- persp(d,a,z_new, theta=20, phi=15, col="darkgoldenrod1",border = 'white', lwd = 0.2, expand = 0.5, ticktype = "detailed",
             xlab=" Disease frequency", ylab=" Risk allele frequency", zlab="", cex.axis = 0.35, cex.lab = "0.5", xlim = c(0.025, .97), zlim = c(0,1))
  
  points_x <- c(0.1, 0.5)  # Example x-coordinates of points
  points_y <- c(0.01, 0.05)  # Example y-coordinates of points
  points_z <- c(0.091, 0.052)  # Example z-coordinates of points
  
  points_trans <- trans3d(points_x, points_y, points_z, pmat = p)  # Transform points to plot coordinates
  
  points(points_trans$x, points_trans$y, col = "red", pch = 16)  # Plot points on the 3D plot
  
  line_x <- c(0.1, 0.5)  # Example x-coordinates of line endpoints
  line_y <- c(0.01, 0.05)  # Example y-coordinates of line endpoints
  line_z <- c(0.091, 0.052)  # Example z-coordinates of line endpoints
  
  line_trans <- trans3d(line_x, line_y, line_z, pmat = p)  # Transform line endpoints to plot coordinates
  
  lines(line_trans$x, line_trans$y, col = "blue", lwd = 2)  # Draw the line on the plot
  
  
  legend("topleft", legend = "a)", inset= c(-0.2, -.01), bty = "n")
  mtext(expression(lambda), side = 2, line = 0.6, at = -0.05, las = 2, cex = 0.5)

} 

ncpsurfaceplotfigure_protective <- function(b, t){
  d = seq(0, 0.99, by = 0.03)
  a = seq(0, 0.99, by = 0.03)
  #t = 1
  #b = 0.5
  z <- outer(d,a, function(d,a) ((a*t*(-1+b+d)^2*(-1+t*d))/((b*(-1+t)+t*(-1+d))*(1+a*b*(-1+t)+a*t*(-1+d)-d)*d)))
  z_new <- pmin(z, 1)
  # Plots the surface
  p <- persp(d,a,z_new, theta=-20, phi=15, col="darkgoldenrod1", border = 'white', lwd = 0.2, expand = 0.5, ticktype = "detailed",
             xlab="Disease frequency", ylab="Protective allele frequency", zlab= "", cex.axis = 0.35, cex.lab = "0.5", xlim = c(0.025, .97), zlim = c(0,1))
 
  points_x <- c(0.1, 0.5)  # Example x-coordinates of points
  points_y <- c(0.01, 0.05)  # Example y-coordinates of points
  points_z <- c(0.001, 0.052)  # Example z-coordinates of points
  
  points_trans <- trans3d(points_x, points_y, points_z, pmat = p)  # Transform points to plot coordinates
  
  points(points_trans$x, points_trans$y, col = "red", pch = 16)  # Plot points on the 3D plot
  
  line_x <- c(0.1, 0.5)  # Example x-coordinates of line endpoints
  line_y <- c(0.01, 0.05)  # Example y-coordinates of line endpoints
  line_z <- c(0.001, 0.053)  # Example z-coordinates of line endpoints
  
  line_trans <- trans3d(line_x, line_y, line_z, pmat = p)  # Transform line endpoints to plot coordinates
  
  lines(line_trans$x, line_trans$y, col = "blue", lwd = 2)  # Draw the line on the plot
  
  legend("topleft", legend = "b)", inset= c(-0.2, -.01), bty = "n")
  mtext(expression(lambda), side = 2, line = 0.3, at = 0.05, las = 2, cex = 0.5)

}


#calculates ncp in haploid case
haploidncpcalculation <- function(d, a, samp.size, df, b, allele, x){
  inflationfactor <- seq(0, x, length.out = 1000)
  
  ncpdata <- vector()
  
  for (t in inflationfactor){
    
    if (allele == 'risk'){
      ncp <- ((b-d)^2*a*t*(d*t-1))/(d*(1+b*(t-1)-d*t)*(-1+d+a+a*b*(t-1)-d*a*t))
      
    }else if (allele == 'protective'){
      ncp <- (a*t*(-1+b+d)^2*(-1+t*d))/((b*(-1+t)+t*(-1+d))*(1+a*b*(-1+t)+a*t*(-1+d)-d)*d)
      
    }
    
    ncpdata[match(t, inflationfactor)] <- ncp 
  }
  plot(inflationfactor*(1/x), ncpdata, xlab = "Case Fraction", ylab = "", type = "l", las = 1, bty = "n", pch =19, xlim = c(0, 1))
  mtext(expression(lambda), side = 2, line = 3.5, las = 2)
}

#calculates ncp in diploid case
diploidncpcalculation <- function(d, a, b, h, samp.size, analysis, x, df = 2){
  inflationfactor <- seq(0.001, x, length.out = 1000)
  
  ncpdata <- vector()
  
  for (t in inflationfactor){
    
    if (analysis == 'dominant'){
      ncp <- (t*(a-1)*(a*(-2*h*(a-1)+a)^2)*((b-d)^2)*(t*d-1))/(d*(d-1+(2*h*(a-1)-a)*a*(b-1-b*t+t*d))*(-2-b*(t-1)*(2*h*(a-1)-a)*(a-1)-a+2*d+a*(a+t*d-t*a*d)+2*h*((t-1)*d+(a-2)*a*(-1+t*d))))
      
    }else if (analysis == 'recessive'){
      ncp <- (t*(a^2)*((b-d)^2)*(t*d-1))/(d*(-1+b-b*t+t*d)*(1-d+(a^2)*(-1+b-b*t+t*d)))
      
    }else{
      ncp <- -(t*a*((b-d)^2)*(-1+t*d)*(2*(h^2)*(-1+a)*(1+2*(-1+a)*a)*(1+b*(-1+t)-t*d)+a*(-1+d+(a^2)*(1+b*(-1+t)-t*d))+h*a*(-b*(-1+t)*((1-2*a)^2)+(-1+t)*d+4*(-1+a)*a*(-1+t*d))))/((d*(-1+b-b*t+t*d)*(-1-2*h*(-1+a)*a+(a^2)-b*(-1+t)*(h+2*h*(-1+a)*a-(a^2))+d-t*(a^2)*d+h*(-1+t+2*t*(-1+a)*a)*d)*(-1+d+(2*h*(-1+a)-a)*a*(-1+b-b*t+t*d))))
      
    }
    
      ncpdata[match(t, inflationfactor)] <- ncp + (df/samp.size)
  }
  plot(inflationfactor*(1/x), ncpdata,  xlab = "Case Fraction", ylab = "", type = 'l', las = 1, bty = "n", pch =19, xlim = c(0, 1), ylim = c(0, max(ncpdata)*1.1))
  mtext(expression(lambda), side = 2, line = 3.5, las = 2)
} 

#figure showing haploid and diploid cases at different penetrance levels (0.2, 0.5, 1)
haploidpenetrancefigure <- function(d, a, penetrance_levels){
  allele <- 'risk'
  par(mfrow = c(1,3), mgp = c(3.5, 1, 0), mar = c(5.1, 6, 2.1, 2.1), xpd=NA)
  haploidncpcalculation(d, a, 10000, 1, penetrance_levels[1], allele, x = (1/d))
  legend("topleft", legend = "a)", inset= c(-0.7, -.2225), bty = "n")
  par(mar = c(5.1, 4.1, 2.1, 2.1))
  haploidncpcalculation(d, a, 10000, 1, penetrance_levels[2], "risk", 1/d)
  legend("topleft", legend = "b)", inset= c(-0.6, -.2), bty = "n")
  haploidncpcalculation(d, a, 10000, 1, penetrance_levels[3], "risk", 1/d)
  legend("topleft", legend = "c)", inset= c(-0.6, -.2), bty = "n")
}

diploidpenetrancefigure <- function(d, a, penetrance_levels){
  allele <- 'unknown'
  par(mfrow = c(1,3), mgp = c(3.5, 1, 0), mar = c(5.1, 6, 2.1, 2.1), xpd=NA)
  diploidncpcalculation(d, a, penetrance_levels[1], 0.5, 10000, allele, 1/d)
  legend("topleft", legend = "d)", inset= c(-0.7, -.2225), bty = "n")
  par(mar = c(5.1, 4.1, 2.1, 2.1))
  diploidncpcalculation(d, a, penetrance_levels[2], 0.5, 10000, allele, 1/d)
  legend("topleft", legend = "e)", inset= c(-0.6, -.2225), bty = "n")
  diploidncpcalculation(d, a, 1, penetrance_levels[3], 10000, allele, 1/d)
  legend("topleft", legend = "f)", inset= c(-0.6, -.2225), bty = "n")
}

dominancefigure <- function(d, a, dominancelevels){
  allele <- "unknown"
  par(mfrow = c(2,3), mgp = c(3, 1, 0), mar = c(5.1, 6, 2.1, 2.1), xpd=NA, las = 2)
  diploidncpcalculation(d, a, 1, dominancelevels[1], 10000, allele, 1/d)
  legend("topleft", legend = "a)", inset= c(-0.7, -.2), bty = "n")
  diploidncpcalculation(d, a, 1, dominancelevels[2], 10000, allele, 1/d)
  legend("topleft", legend = "b)", inset= c(-0.7, -.2), bty = "n")
  diploidncpcalculation(d, a, 1, dominancelevels[3], 10000, allele, 1/d)
  legend("topleft", legend = "c)", inset= c(-0.7, -.2), bty = "n")
  diploidncpcalculation(d, a, 0.5, dominancelevels[1], 10000, allele, 1/d)
  legend("topleft", legend = "d)", inset= c(-0.7, -.2), bty = "n")
  diploidncpcalculation(d, a, 0.5, dominancelevels[2], 10000, allele, 1/d)
  legend("topleft", legend = "e)", inset= c(-0.7, -.2), bty = "n")
  diploidncpcalculation(d, a, 0.5, dominancelevels[3], 10000, allele, 1/d)
  legend("topleft", legend = "f)", inset= c(-0.7, -.2), bty = "n")
}

riskvsprotectivefigure <- function(d, a, penetrance){
  allele <- 'risk'
 # par(mgp = c(3, 1, 0), mar = c(3, 3, 3, 3) , axes = TRUE)
  haploidncpcalculation(d, a, 10000, 1, penetrance, allele, 1/d)
  legend("topleft", legend = "c)", inset= c(-0.38, -.1), bty = "n")
  allele <- 'protective'
  haploidncpcalculation(d, a, 10000, 1, penetrance, allele, 1/d)
  legend("topleft", legend = "d)", inset= c(-0.38, -.1), bty = "n")
}

creatediploidpopulation <- function(d, a, b, h, analysis, pop.size = 1000000){
  
  pop.genotypes <- rbinom(pop.size, 2, a)
  pop.phenotypes <- numeric(length(pop.genotypes))
  
  if (analysis == 'dominant'){
    pop.phenotypes[pop.genotypes == 2] <- rbinom(sum(pop.genotypes == 2), 1, b)
    pop.phenotypes[pop.genotypes == 1] <- rbinom(sum(pop.genotypes == 1), 1, 
                                                 (b*h+((1-h)*(-2*b*h*a-b*a^2+2*b*h*a^2+d))/((-1+a)*(-1-a+2*h*a))))
    pop.phenotypes[pop.genotypes == 0] <- rbinom(sum(pop.genotypes == 0), 1, 
                                                 ((-2*b*h*a-b*a^2+2*b*h*a^2+d)/((-1+a)*(-1-a+2*h*a))))  
    pop.genotypes[pop.genotypes == 2] <- 1
    
  }else if (analysis == 'recessive'){
    pop.phenotypes[pop.genotypes == 2] <- rbinom(sum(pop.genotypes == 2), 1, b)
    pop.phenotypes[pop.genotypes == 1] <- rbinom(sum(pop.genotypes == 1), 1, 
                                                 (b*h+((1-h)*(-2*b*h*a-b*a^2+2*b*h*a^2+d))/((-1+a)*(-1-a+2*h*a))))
    pop.phenotypes[pop.genotypes == 0] <- rbinom(sum(pop.genotypes == 0), 1, 
                                                 ((-2*b*h*a-b*a^2+2*b*h*a^2+d)/((-1+a)*(-1-a+2*h*a))))  
    pop.genotypes[pop.genotypes == 1] <- 0
    pop.genotypes[pop.genotypes== 2] <-1
    
  }else{
    #assigns phenotypes based on the three genotypes using probabilities from chi-squared tables
    pop.phenotypes[pop.genotypes == 2] <- rbinom(sum(pop.genotypes == 2), 1, b)
    pop.phenotypes[pop.genotypes == 1] <- rbinom(sum(pop.genotypes == 1), 1, 
                                                 (b*h+((1-h)*(-2*b*h*a-b*a^2+2*b*h*a^2+d))/((-1+a)*(-1-a+2*h*a))))
    pop.phenotypes[pop.genotypes == 0] <- rbinom(sum(pop.genotypes == 0), 1, 
                                                 ((-2*b*h*a-b*a^2+2*b*h*a^2+d)/((-1+a)*(-1-a+2*h*a))))  
  }
  
  pop <- data.frame(pop.genotypes, pop.phenotypes)
  return(pop)
}

powercalculation <- function(d, a, samp.size, analysis, alpha, df, b, h, x){
  
  inflation <- seq(0.0001, x, length.out = 1000)
  
  data <- vector()
  
  for (t in inflation){
    
    if (analysis == 'dominant'){
      ncp.base <- (t*(a-1)*(a*(-2*h*(a-1)+a)^2)*((b-d)^2)*(t*d-1))/
        (d*(d-1+(2*h*(a-1)-a)*a*(b-1-b*t+t*d))*(-2-b*(t-1)*(2*h*(a-1)-a)*
                                                  (a-1)-a+2*d+a*(a+t*d-t*a*d)+2*h*((t-1)*d+(a-2)*a*(-1+t*d))))
      
    }else if (analysis == 'recessive'){
      ncp.base <- (t*(a^2)*((b-d)^2)*(t*d-1))/(d*(-1+b-b*t+t*d)*(1-d+(a^2)
                                                                 *(-1+b-b*t+t*d)))
      
    }else{
      ncp.base <- -(t*a*((b-d)^2)*(-1+t*d)*(2*(h^2)*(-1+a)*(1+2*(-1+a)*a)*
                                              (1+b*(-1+t)-t*d)+a*(-1+d+(a^2)*(1+b*(-1+t)-t*d))+h*a*(-b*(-1+t)*((1-2*a)^2)+
                                                                                                      (-1+t)*d+4*(-1+a)*a*(-1+t*d))))/((d*(-1+b-b*t+t*d)*(-1-2*h*(-1+a)*a+(a^2)-b*(-1+t)*
                                                                                                                                                            (h+2*h*(-1+a)*a-(a^2))+d-t*(a^2)*d+h*(-1+t+2*t*(-1+a)*a)*d)*(-1+d+(2*h*(-1+a)-a)*a*(-1+b-b*t+t*d))))
    }
    
    ncp <- ncp.base * samp.size
    
    criticalvalue <- qchisq(1-alpha, df)
    

    power <- pchisq(criticalvalue, df, ncp, lower.tail=FALSE, log.p = FALSE)
    data[match(t, inflation)] <- power
  }
  lines(inflation*d, data, xlab = "Case Fraction", ylab = "Power", type = 'l', las = 1, bty = "n", pch =19)
  
}

ncppowersimulation <- function(d, a, b, h, samp.size, analysis, 
                               inflation = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3),
                               n.sim = 1000, pop.size = 1000000, alpha = 0.00000005, df = 2){
  #creates population
  pop <- creatediploidpopulation(d, a, b, h, analysis, pop.size)
  
  #inflation vector contains the different inflation factors that will be used
  power <- numeric(length(inflation))
  upperpower <- numeric(length(inflation))
  lowerpower <- numeric(length(inflation))
  
  
  #subsets the population into cases and controls based on phenotype
  cases <- pop[pop$pop.phenotypes == 1,]
  controls <- pop[pop$pop.phenotypes == 0,]
  
  #loops through the different inflation factors to get the chi-squared value
  for(t in inflation){
    chisq.stats <- numeric(n.sim)
    chisq.ps <- numeric(n.sim)
    case.inflation <-t
    
    #repeats for the desired number of simulations
    for(i in 1:n.sim){
      num.casesample <- case.inflation*d*samp.size
      num.controlsample <- samp.size - num.casesample
      casesample <- (cases[sample(nrow(cases), num.casesample),])
      controlsample <- (controls[sample(nrow(controls),num.controlsample),])
      
      caselist <- c(nrow(casesample[casesample$pop.genotypes == 2,]), 
                    nrow(casesample[casesample$pop.genotypes == 1,]), 
                    nrow(casesample[casesample$pop.genotypes == 0,]) )
      controllist <- c(nrow(controlsample[controlsample$pop.genotypes == 2,]), 
                       nrow(controlsample[controlsample$pop.genotypes == 1,]), 
                       nrow(controlsample[controlsample$pop.genotypes == 0,]) ) 
      
      sample <- rbind(casesample, controlsample)
      chisq.fit <- chisq.test(sample$pop.genotypes, sample$pop.phenotypes)
      
      #stores chi-squared statistic and p-value
      chisq.stats[i] <- chisq.fit$statistic
      chisq.ps[i] <- chisq.fit$p.value
    }
    
    #approximation for power
    power[match(t, inflation)] <- mean(chisq.ps <= alpha)
    
    #confidence intervals
    upperpower[match(t,inflation)] <- mean(chisq.ps <= alpha) + 
      (2*sqrt((power[match(t, inflation)])*(1-(power[match(t, inflation)])))/sqrt(n.sim))
    lowerpower[match(t,inflation)] <- mean(chisq.ps <= alpha) - 
      (2*sqrt((power[match(t, inflation)])*(1-(power[match(t, inflation)])))/sqrt(n.sim))
    
  }
  
  #plots data points for power
  plot(inflation*d, power, xlab = "Case Fraction", ylab = "Power", las = 1, bty = "n", pch =19, xlim = c(0, 1))
  for (t in inflation){
    x2 <- c(t*d, t*d)
    y2 <- c((upperpower[match(t,inflation)]), (lowerpower[match(t,inflation)]))
    lines(x2, y2)
  }
  
  #overlays values of power calculated using calculated ncp and prints case frac
  #tion at which power is greatest
  powercalculation(d, a, samp.size, analysis, alpha, df, b, h, max(inflation))
}

ncppowersimfigure <- function(d, a){
  allele <- "unknown"
  par(mfrow = c(3,3), mgp = c(2.5,0.5,0),  mar = c(4.1, 4.1, 1.2, 2.1),  cex.axis = 1, xpd=TRUE)
  ncppowersimulation(0.3, 0.1, 1, 0.1, 1450, analysis = allele, n.sim = 500)
  legend("topleft", legend = "a)", inset= c(-0.5, -0.2), bty = "n")
  ncppowersimulation(0.3, 0.1, 1, 0.5, 350, analysis = allele, n.sim = 500)
  legend("topleft", legend = "b)", inset= c(-0.5, -0.2), bty = "n")
  ncppowersimulation(0.3, 0.1, 1, .9, 100,  analysis = allele, n.sim = 500)
  legend("topleft", legend = "c)", inset= c(-0.5, -0.2), bty = "n")
  ncppowersimulation(0.3, 0.1, 0.6, .1, 8200,  analysis = allele, n.sim = 100)
  legend("topleft", legend = "d)", inset= c(-0.5, -0.2), bty = "n")
  ncppowersimulation(0.3, 0.1, 0.6, 0.5, 1700, analysis = allele, n.sim = 500)
  legend("topleft", legend = "e)", inset= c(-0.5, -0.2), bty = "n")
  ncppowersimulation(0.3, 0.1, 0.6, .9, 550,  analysis = allele, n.sim = 500)
  legend("topleft", legend = "f)", inset= c(-0.5, -0.2), bty = "n")
  ncppowersimulation(0.3, 0.1, 0.4, .1, 70000, analysis = allele,n.sim = 500)
  legend("topleft", legend = "g)", inset= c(-0.5, -0.2), bty = "n")
  ncppowersimulation(0.3, 0.1, 0.4, 0.5, 15000, analysis = allele, n.sim = 500)
  legend("topleft", legend = "h)", inset= c(-0.5, -0.2), bty = "n")
  ncppowersimulation(0.3, 0.1, 0.4, .9, 4500, analysis = allele, n.sim = 500)
  legend("topleft", legend = "i)", inset= c(-0.5, -0.2), bty = "n")
  
}


armitagepowersimulation <- function(d, a, b, h, samp.size, analysis = 'unknown', 
                                    inflation = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), n.sim = 200, 
                                    pop.size = 1000000, alpha = 0.00000005, df = 2){
  pop <- creatediploidpopulation(d, a, b, h, analysis, pop.size)
  
  powerdata <- numeric(length(inflation))
  pvalues <- numeric(length(inflation))
  
  cases <- pop[pop$pop.phenotypes == 1,]
  controls <- pop[pop$pop.phenotypes == 0,]
  
  for(t in inflation){
    armitage.stats <- numeric(n.sim)
    armitage.ps <- numeric(n.sim)
    case.inflation <-t
    
    for(i in 1:n.sim){
      num.casesample <- case.inflation*d*samp.size
      num.controlsample <- samp.size - num.casesample
      casesample <- (cases[sample(nrow(cases), num.casesample),])
      controlsample <- (controls[sample(nrow(controls),num.controlsample),])
      
      caselist <- c(nrow(casesample[casesample$pop.genotypes == 0,]),
                    nrow(casesample[casesample$pop.genotypes == 1,]), 
                    nrow(casesample[casesample$pop.genotypes == 2,]) 
      )
      controllist <- c(nrow(controlsample[controlsample$pop.genotypes == 0,]),
                       nrow(controlsample[controlsample$pop.genotypes == 1,]),
                       nrow(controlsample[controlsample$pop.genotypes == 2,]) 
      )
      
      contingencytable <- matrix(c(caselist, controllist), ncol = 3, byrow = TRUE)
      #print(contingencytable)
      armitage.fit <- CochranArmitageTest(contingencytable, alternative = "one.sided")
      
      armitage.stats[i] <- armitage.fit$statistic
      armitage.ps[i] <- armitage.fit$p.value
    }
    #data[match(t,inflation)] <- mean(armitage.stats)
    powerdata[match(t,inflation)] <- mean(armitage.ps <= alpha)
    pvalues[match(t,inflation)] <- mean(armitage.ps)
  }
  #plots the data points for ncp
  plot(inflation*d, powerdata, 
       xlab = "Case Fraction", ylab = "Power from Armitage trend test", las = 1, bty = "n", pch =19)
}

armitagetrendfigure <- function(d, a){
  #install.packages("DescTools")
  library(DescTools)
  allele <- "unknown"
  par(mfrow = c(3,3), mgp = c(2.5,0.5,0),  mar = c(4.1, 4.1, 1.2, 2.1),  cex.axis = 1, xpd=TRUE)
  armitagepowersimulation(0.3, 0.1, 1, 0.1, 1450, analysis = allele, n.sim = 500)
  legend("topleft", legend = "a)", inset= c(-0.5, -0.2), bty = "n")
  armitagepowersimulation(0.3, 0.1, 1, 0.5, 350, analysis = allele, n.sim = 500)
  legend("topleft", legend = "b)", inset= c(-0.5, -0.2), bty = "n")
  armitagepowersimulation(0.3, 0.1, 1, .9, 100,  analysis = allele, n.sim = 500)
  legend("topleft", legend = "c)", inset= c(-0.5, -0.2), bty = "n")
  armitagepowersimulation(0.3, 0.1, 0.6, .1, 8200,  analysis = allele, n.sim = 100)
  legend("topleft", legend = "d)", inset= c(-0.5, -0.2), bty = "n")
  armitagepowersimulation(0.3, 0.1, 0.6, 0.5, 1700, analysis = allele, n.sim = 500)
  legend("topleft", legend = "e)", inset= c(-0.5, -0.2), bty = "n")
  armitagepowersimulation(0.3, 0.1, 0.6, .9, 550,  analysis = allele, n.sim = 500)
  legend("topleft", legend = "f)", inset= c(-0.5, -0.2), bty = "n")
  armitagepowersimulation(0.3, 0.1, 0.4, .1, 70000, analysis = allele,n.sim = 500)
  legend("topleft", legend = "g)", inset= c(-0.5, -0.2), bty = "n")
  armitagepowersimulation(0.3, 0.1, 0.4, 0.5, 15000, analysis = allele, n.sim = 500)
  legend("topleft", legend = "h)", inset= c(-0.5, -0.2), bty = "n")
  armitagepowersimulation(0.3, 0.1, 0.4, .9, 4500, analysis = allele, n.sim = 500)
  legend("topleft", legend = "i)", inset= c(-0.5, -0.2), bty = "n")
}


logisticpowersimulation <- function(d, a, b, h, samp.size, analysis = 'unknown', inflation = 
                                      c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), n.sim = 200, pop.size = 1000000, alpha = 0.00000005, df = 2){
  #creates population
  pop <- creatediploidpopulation(d, a, b, h, analysis, pop.size)
  
  #inflation vector contains the different inflation factors that will be used
  data <- numeric(length(inflation))
  upperncp <- numeric(length(inflation))
  lowerncp <- numeric(length(inflation))
  
  #subsets the population into cases and controls based on phenotype
  cases <- pop[pop$pop.phenotypes == 1,]
  controls <- pop[pop$pop.phenotypes == 0,]
  
  #loops through the different inflation factors to run chi-squared test
  for(t in inflation){
    logistic.ps <- numeric(n.sim)
    case.inflation <-t
    
    #repeats for the desired number of simulations
    for(i in 1:n.sim){
      num.casesample <- case.inflation*d*samp.size
      num.controlsample <- samp.size - num.casesample
      casesample <- (cases[sample(nrow(cases), num.casesample),])
      controlsample <- (controls[sample(nrow(controls),num.controlsample),])
      
      sample <- rbind(casesample, controlsample)
      logistic <- glm(sample$pop.phenotypes ~ sample$pop.genotypes,family = 'binomial'( link = 'logit'))
      logistic.ps[i] <- summary(logistic)$coefficients[2,4]
      #stores chi-squared statistic and p-value
      
    }
    data[match(t,inflation)] <- mean(logistic.ps <= alpha)
  }
  #plots the data points for ncp
  plot(inflation*d, data, 
       xlab = "Case Fraction", ylab = "Power from logistic regression", las = 1, bty = "n", pch =19)
}


logisticregressionfigure <- function(d, a){
  pd <- d
  pa <- a
  allele <- "unknown"
  par(mfrow = c(3,3), mgp = c(2.5,0.5,0),  mar = c(4.1, 4.1, 1.2, 2.1),  cex.axis = 1, xpd=TRUE)
  logisticpowersimulation(0.3, 0.1, 1, 0.1, 1450, analysis = allele, n.sim = 500)
  legend("topleft", legend = "a)", inset= c(-0.5, -0.2), bty = "n")
  logisticpowersimulation(0.3, 0.1, 1, 0.5, 350, analysis = allele, n.sim = 500)
  legend("topleft", legend = "b)", inset= c(-0.5, -0.2), bty = "n")
  logisticpowersimulation(0.3, 0.1, 1, .9, 100,  analysis = allele, n.sim = 500)
  legend("topleft", legend = "c)", inset= c(-0.5, -0.2), bty = "n")
  logisticpowersimulation(0.3, 0.1, 0.6, .1, 8200,  analysis = allele, n.sim = 100)
  legend("topleft", legend = "d)", inset= c(-0.5, -0.2), bty = "n")
  logisticpowersimulation(0.3, 0.1, 0.6, 0.5, 1700, analysis = allele, n.sim = 500)
  legend("topleft", legend = "e)", inset= c(-0.5, -0.2), bty = "n")
  logisticpowersimulation(0.3, 0.1, 0.6, .9, 550,  analysis = allele, n.sim = 500)
  legend("topleft", legend = "f)", inset= c(-0.5, -0.2), bty = "n")
  logisticpowersimulation(0.3, 0.1, 0.4, .1, 70000, analysis = allele,n.sim = 500)
  legend("topleft", legend = "g)", inset= c(-0.5, -0.2), bty = "n")
  logisticpowersimulation(0.3, 0.1, 0.4, 0.5, 15000, analysis = allele, n.sim = 500)
  legend("topleft", legend = "h)", inset= c(-0.5, -0.2), bty = "n")
  logisticpowersimulation(0.3, 0.1, 0.4, .9, 4500, analysis = allele, n.sim = 500)
  legend("topleft", legend = "i)", inset= c(-0.5, -0.2), bty = "n")
  
  
}