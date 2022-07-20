Rwcs_in_poly = function(x, y, poly_x, poly_y){
  if (is.matrix(x) || is.data.frame(x)) {
    if (ncol(x) == 1) {
      x = x[, 1]
    }else if (ncol(x) == 2) {
      y = x[, 2]
      x = x[, 1]
    }
  }
  
  if (is.matrix(y) || is.data.frame(y)) {
    if (ncol(y) == 1) {
      y = y[, 1]
    }
  }
  
  if (is.matrix(poly_x) || is.data.frame(poly_x)) {
    if (ncol(poly_x) == 1) {
      poly_x = poly_x[, 1]
    }else if (ncol(poly_x) == 2) {
      poly_y = poly_x[, 2]
      poly_x = poly_x[, 1]
    }
  }
  
  if (is.matrix(poly_y) || is.data.frame(poly_y)) {
    if (ncol(poly_y) == 1) {
      poly_y = poly_y[, 1]
    }
  }
  
  return(.pnpoly(x, y, poly_x, poly_y))
}

Rwcs_overlap_poly = function(poly_x1, poly_y1, poly_x2, poly_y2, plot=FALSE){
  if(plot){
    xlim = c(min(poly_x1, poly_x2), max(poly_x1, poly_x2))
    ylim = c(min(poly_y1, poly_y2), max(poly_y1, poly_y2))
    magplot(NA, NA, xlim=xlim, ylim=ylim)
    polygon(poly_x1, poly_y1, border='red')
    polygon(poly_x2, poly_y2, border='blue')
  }
  
  edge_check = any(c(.pnpoly(poly_x1, poly_y1, poly_x2, poly_y2), .pnpoly(poly_x2, poly_y2, poly_x1, poly_y1)))
  
  if(edge_check){
    return(edge_check)
  }else{
    N1 = length(poly_x1)
    N2 = length(poly_x2)
    
    if(poly_x1[1] == poly_x1[N1] & poly_y1[1] == poly_y1[N1]){
      poly_x1 = poly_x1[-N1]
      poly_y1 = poly_y1[-N1]
      N1 = N1 - 1L
    }
    
    if(poly_x2[1] == poly_x2[N2] & poly_y2[1] == poly_y2[N2]){
      poly_x2 = poly_x2[-N2]
      poly_y2 = poly_y2[-N2]
      N2 = N2 - 1L
    }
    
    line_start_x1 = poly_x1
    line_start_y1 = poly_y1
    line_end_x1   = c(poly_x1[2:length(poly_x1)], poly_x1[1])
    line_end_y1   = c(poly_y1[2:length(poly_x1)], poly_y1[1])
    
    line_start_x2 = poly_x2
    line_start_y2 = poly_y2
    line_end_x2   = c(poly_x2[2:length(poly_x2)], poly_x2[1])
    line_end_y2   = c(poly_y2[2:length(poly_x2)], poly_y2[1])
    
    #first we use times for the repeating
    start1 = cbind(rep(line_start_x1, times=N2), rep(line_start_y1, times=N2))
    end1   = cbind(rep(line_end_x1, times=N2), rep(line_end_y1, times=N2))
    #then we use each for the repeating
    start2 = cbind(rep(line_start_x2, each=N1), rep(line_start_y2, each=N1))
    end2   = cbind(rep(line_end_x2, each=N1), rep(line_end_y2, each=N1))
    
    return(any(.line_intersect(start1=start1, end1=end1, start2=start2, end2=end2)))
  }
}

.line_intersect = function(start1 = cbind(0,0), end1 = cbind(1,1),
                          start2 = cbind(1,0), end2 = cbind(0,1)){
  #tells you whether line element intersect
  start1 = rbind(start1)
  end1 = rbind(end1)
  start2 = rbind(start2)
  end2 = rbind(end2)
  
  dx1 = end1[,1] - start1[,1]
  dx2 = end2[,1] - start2[,1]
  dy1 = end1[,2] - start1[,2]
  dy2 = end2[,2] - start2[,2]
  
  p1 = dy2*(end2[,1] - start1[,1]) - dx2*(end2[,2] - start1[,2])
  p2 = dy2*(end2[,1] - end1[,1]) -   dx2*(end2[,2] - end1[,2])
  p3 = dy1*(end1[,1] - start2[,1]) - dx1*(end1[,2] - start2[,2])
  p4 = dy1*(end1[,1] - end2[,1]) -   dx1*(end1[,2] - end2[,2])
  
  return((p1*p2 <= 0) & (p3*p4 <=0))
}

.linedist = function(start1=cbind(0,0,0), end1=cbind(1,1,0), start2=cbind(1,0,0), end2=cbind(0,1,0)){
  
  if(dim(start1)[2]!=3){stop('start1 must be Nx3 dimensions')}
  if(dim(end1)[2]!=3){stop('end1 must be Nx3 dimensions')}
  if(dim(start2)[2]!=3){stop('start2 must be Nx3 dimensions')}
  if(dim(end2)[2]!=3){stop('end2 must be Nx3 dimensions')}
  
  avec = end1 - start1
  bvec = end2 - start2
  cvec = start2 - start1
  
  cross_temp = cbind(avec[,2]*bvec[,3]-avec[,3]*bvec[,2],
                     avec[,3]*bvec[,1]-avec[,1]*bvec[,3],
                     avec[,1]*bvec[,2]-avec[,2]*bvec[,1])
  dot_temp=abs(cvec[,1]*cross_temp[,1]+cvec[,2]*cross_temp[,2]+cvec[,3]*cross_temp[,3])
  cross_temp=sqrt(cross_temp[,1]^2+cross_temp[,2]^2+cross_temp[,3]^2)
  return(as.numeric(unlist(dot_temp/cross_temp)))
}

.line_cross = function(m1, c1, m2, c2){
  #gives location of line intersection
  x_cross = (c2 - c1)/(m1 - m2)
  y_cross = (c2*m1 - c1*m2)/(m1 - m2)
  return(list(x_cross=x_cross, y_cross=y_cross))
}
