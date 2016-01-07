library(shiny)
library(png)
shinyServer(function(input, output) {
  

  output$myImage <- renderImage({
    appDir <- system.file("ExamplePlots", package = "Crambled")
    
    mycell<-input$mycell
    mydepth<-input$mydepth
    mywidth<-input$mywidth
    myheight<-input$myheight
    myxlim1<-input$myxlim1
    myxlim2<-input$myxlim2
    
    Dmat<-matrix(ncol=8,nrow=8,NA)
    Amat<-matrix(ncol=8,nrow=8,NA)
    Bmat<-matrix(ncol=8,nrow=8,NA)
    
    for(i in 1:5){
      for(j in i:5){
        Dmat[i,j]<-2*mydepth*(1-mycell)+mycell*(i+j-2)*mydepth
        Amat[i,j]<-(1-mycell+mycell*(i-1))/( 2*(1-mycell)+(i+j-2)*mycell) 
              
        mydens<-dbinom(0:round(Dmat[i,j]),round(Dmat[i,j]),Amat[i,j])
        myval<-(0:round(Dmat[i,j]))/round(Dmat[i,j])
        myval[myval>0.5]<-1-myval[myval>0.5]
        
        Bmat[i,j]<-sum(myval*mydens)
        
      }  
    }
    
    ABlabels<-c("0",rep("",7),"A","AB",rep("",6),"AA","AAB","AABB",rep("",5),"AAA","AAAB","AAABB","AAABBB",rep("",4),"AAAA","AAAAB","AAAABB","AAAABBB","AAAABBBB")
    
    outfile <- tempfile(fileext='.png')
    
    png(outfile, width=mywidth, height=myheight)
    plot(1,1,type="n",ylim=c(0,0.5),xlim=c(myxlim1,myxlim2),axes=F,xlab="",ylab="",main="")
    points(Dmat,Bmat,pch=21,col="black",bg="yellow")
    text(as.vector(Dmat)[1:37],as.vector(Bmat)[1:37],ABlabels,col="red",cex=0.75,pos=1)
    dev.off()
    
    inFile <- input$file1

    if (is.null(inFile)){inFile$datapath<-file.path(appDir,"HCC1395-shiny.png")}
    
    underpng<-png::readPNG(inFile$datapath)
    overpng<-png::readPNG(outfile)
    matte<-(overpng[,,1]+overpng[,,2]+overpng[,,3])<3
for(i in 1:3){
  tmpu<-underpng[,,i]
  tmpo<-overpng[,,i]
  tmpu[matte]<-tmpo[matte]
  underpng[,,i]<-tmpu
}
    png::writePNG(underpng,target=outfile)
    list(src = outfile,
         alt = paste("Image")) 
  }, deleteFile = TRUE)
})