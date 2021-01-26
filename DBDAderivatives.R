# Copied from:
# Jags-Ydich-XnomSsubj-Mbernbeta.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
plotMCMC = function( codaSamples , data , yName="y" , sName="s", compVal=0.5 , rope=NULL , 
                     compValDiff=0.0 , ropeDiff=NULL , credMass=0.95,
                     saveName=NULL , saveType="jpg" ) {
  #-----------------------------------------------------------------------------
  # N.B.: This function expects the data to be a data frame, 
  # with one component named y being a vector of integer 0,1 values,
  # and one component named s being a factor of subject identifiers.
  y = data[,yName]
  s = as.numeric(data[,sName]) # converts character to consecutive integer levels
  # Now plot the posterior:
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  Ntheta = length(grep("theta",colnames(mcmcMat)))
  openGraph(width=2.5*Ntheta,height=2.0*Ntheta)
  par( mfrow=c(Ntheta,Ntheta) )
  for ( t1Idx in 1:(Ntheta) ) {
    for ( t2Idx in (1):Ntheta ) {
      parName1 = paste0("theta[",t1Idx,"]")
      parName2 = paste0("theta[",t2Idx,"]")
      if ( t1Idx > t2Idx) {  
        # plot.new() # empty plot, advance to next
        par( mar=c(3.5,3.5,1,1) , mgp=c(2.0,0.7,0) )
        nToPlot = 700
        ptIdx = round(seq(1,chainLength,length=nToPlot))
        plot ( mcmcMat[ptIdx,parName2] , mcmcMat[ptIdx,parName1] , cex.lab=1.75 ,
               xlab=parse(text=parName2) , ylab=parse(text=parName1) , col="skyblue" )
      } else if ( t1Idx == t2Idx ) {
        par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) )
        postInfo = plotPost( mcmcMat[,parName1] , cex.lab = 1.75 , 
                             compVal=compVal , ROPE=rope , credMass=credMass, cex.main=1.5 ,
                             xlab=parse(text=parName1) , main="" )
        includeRows = ( s == t1Idx ) # identify rows of this subject in data
        dataPropor = sum(y[includeRows])/sum(includeRows) 
        points( dataPropor , 0 , pch="+" , col="red" , cex=3 )
      } else if ( t1Idx < t2Idx ) {
        par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) )
        postInfo = plotPost(mcmcMat[,parName1]-mcmcMat[,parName2] , cex.lab = 1.75 , 
                            compVal=compValDiff , ROPE=ropeDiff , credMass=credMass, cex.main=1.5 ,
                            xlab=parse(text=paste0(parName1,"-",parName2)) , main="" )
        includeRows1 = ( s == t1Idx ) # identify rows of this subject in data
        dataPropor1 = sum(y[includeRows1])/sum(includeRows1) 
        includeRows2 = ( s == t2Idx ) # identify rows of this subject in data
        dataPropor2 = sum(y[includeRows2])/sum(includeRows2) 
        points( dataPropor1-dataPropor2 , 0 , pch="+" , col="red" , cex=3 )
      }
    }
  }
  #-----------------------------------------------------------------------------  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
}


# Adapted from Jags-Ymet-Xnom1fac-MnormalHom.R
plotMCMCwithContrasts = function( codaSamples , 
                                  datFrm , yName="y" , xName="x" , contrasts=NULL , credMass=0.95,
                                  saveName=NULL , saveType="jpg" ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  y = datFrm[,yName]
  x = as.numeric(as.factor(datFrm[,xName]))
  xlevels = levels(as.factor(datFrm[,xName]))
  # Display data with posterior predictive distributions
  openGraph(width=min(10,1.25*length(xlevels)),height=5)
  par(mar=c(3,3,2,0.5)) # number of margin lines: bottom,left,top,right
  par(mgp=c(1.75,0.5,0)) # which margin lines to use for labels
  plot(-1,0, 
       xlim=c(0.1,length(xlevels)+0.1) , 
       xlab=xName , xaxt="n" , ylab=yName ,
       ylim=c(min(y)-0.2*(max(y)-min(y)),max(y)+0.2*(max(y)-min(y))) , 
       main="Data with Posterior Predictive Distrib.")
  axis( 1 , at=1:length(xlevels) , tick=FALSE , lab=xlevels )
  for ( xidx in 1:length(xlevels) ) {
    print(xidx)
    xPlotVal = xidx 
    yVals = y[ x==xidx ]
    points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) , 
            yVals , pch=1 , cex=1.5 , col="red" )
    chainSub = round(seq(1,chainLength,length=20))
    for ( chnIdx in chainSub ) {
      m = mcmcMat[chnIdx,paste("theta[",xidx,"]",sep="")]
      s = 1
      print(chnIdx)
      nu = 1000 # effectively normal instead of mcmcMat[chnIdx,"nu"]
      tlim = qt( c(0.025,0.975) , df=nu )
      yl = m+tlim[1]*s
      yh = m+tlim[2]*s
      ycomb=seq(yl,yh,length=201)
      #ynorm = dnorm(ycomb,mean=m,sd=s)
      #ynorm = 0.67*ynorm/max(ynorm)
      yt = dt( (ycomb-m)/s , df=nu )
      yt = 0.67*yt/max(yt)
      lines( xPlotVal-yt , ycomb , col="skyblue" ) 
    }
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostPred",sep=""), type=saveType)
  }
  if ( !is.null(contrasts) ) {
    if ( is.null(datFrm) | is.null(xName) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      for ( cIdx in 1:length(contrasts) ) {
        thisContrast = contrasts[[cIdx]]
        left = right = rep(FALSE,length(xlevels))
        for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
          left = left | xlevels==thisContrast[[1]][nIdx]
        }
        left = normalize(left)
        for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
          right = right | xlevels==thisContrast[[2]][nIdx]
        }
        right = normalize(right)
        contrastCoef = matrix( left-right , ncol=1 )
        postContrast = ( mcmcMat[,paste("theta[",1:length(xlevels),"]",sep="")] 
                         %*% contrastCoef )
        openGraph(height=4,width=4)
        layout(matrix(1:2,ncol=1))
        plotPost( postContrast , xlab="Difference" ,
                  main=paste0( 
                    paste(thisContrast[[1]],collapse="."), 
                    " vs ",
                    paste(thisContrast[[2]],collapse=".") ) ,
                  compVal=thisContrast$compVal , ROPE=thisContrast$ROPE, credMass=credMass )
        #plotPost( postContrast/1 , xlab="Effect Size" ,
        #         main=paste0( 
        #           paste(thisContrast[[1]],collapse="."), 
        #           "\nvs\n",
        #           paste(thisContrast[[2]],collapse=".") ) ,
        #         compVal=0.0 , 
        #         ROPE=c(-0.1,0.1) )
        if ( !is.null(saveName) ) {
          saveGraph( file=paste0(saveName, paste0( 
            paste(thisContrast[[1]],collapse=""), 
            "_v_",
            paste(thisContrast[[2]],collapse="") ) ), 
            type=saveType )
        }
      }
    }
  } # end if ( !is.null(contrasts) )
}