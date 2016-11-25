from rpy2.robjects import r, numpy2ri
numpy2ri.activate()

# The following R code was copied from https://github.com/meeliskull/prg/blob/master/R_package/prg/R/prg.R
r('''
library(grid)
precision_gain = function(TP,FN,FP,TN) {
  n_pos = TP+FN
  n_neg = FP+TN
  prec_gain = 1-(n_pos*FP)/(n_neg*TP)
  prec_gain[TN+FN==0] = 0
  return(prec_gain)
}
recall_gain = function(TP,FN,FP,TN) {
  n_pos = TP+FN
  n_neg = FP+TN
  rg = 1-(n_pos*FN)/(n_neg*TP)
  rg[TN+FN==0] = 1
  return(rg)
}
# create a table of segments
.create.segments = function(labels,pos_scores,neg_scores) {
  # reorder labels and pos_scores by decreasing pos_scores, using increasing neg_scores in breaking ties
  new_order = order(pos_scores,-neg_scores,decreasing=TRUE)
  labels = labels[new_order]
  pos_scores = pos_scores[new_order]
  neg_scores = neg_scores[new_order]
  # create a table of segments
  segments = data.frame(pos_score=NA,neg_score=NA,pos_count=0,neg_count=rep(0,length(labels)))
  j = 0
  for (i in seq_along(labels)) {
    if ((i==1)||(pos_scores[i-1]!=pos_scores[i])||(neg_scores[i-1]!=neg_scores[i])) {
      j = j + 1
      segments$pos_score[j] = pos_scores[i]
      segments$neg_score[j] = neg_scores[i]
    }
    segments[j,4-labels[i]] = segments[j,4-labels[i]] + 1
  }
  segments = segments[1:j,]
  return(segments)
}
.create.crossing.points = function(points,n_pos,n_neg) {
  n = n_pos+n_neg
  points$is_crossing = 0
  # introduce a crossing point at the crossing through the y-axis
  j = min(which(points$recall_gain>=0))
  if (points$recall_gain[j]>0) { # otherwise there is a point on the boundary and no need for a crossing point
    delta = points[j,,drop=FALSE]-points[j-1,,drop=FALSE]
    if (delta$TP>0) {
      alpha = (n_pos*n_pos/n-points$TP[j-1])/delta$TP
    } else {
      alpha = 0.5
    }
    new_point = points[j-1,,drop=FALSE] + alpha*delta
    new_point$precision_gain = precision_gain(new_point$TP,new_point$FN,new_point$FP,new_point$TN)
    new_point$recall_gain = 0
    new_point$is_crossing = 1
    points = rbind(points,new_point)
    points = points[order(points$index),,drop=FALSE]
  }   
  # now introduce crossing points at the crossings through the non-negative part of the x-axis
  crossings = data.frame()
  x = points$recall_gain
  y = points$precision_gain
  for (i in which((c(y,0)*c(0,y)<0)&(c(1,x)>=0))) {
    cross_x = x[i-1]+(-y[i-1])/(y[i]-y[i-1])*(x[i]-x[i-1])
    delta = points[i,,drop=FALSE]-points[i-1,,drop=FALSE]
    if (delta$TP>0) {
      alpha = (n_pos*n_pos/(n-n_neg*cross_x)-points$TP[i-1])/delta$TP
    } else {
      alpha = (n_neg/n_pos*points$TP[i-1]-points$FP[i-1])/delta$FP
    }
    new_point = points[i-1,,drop=FALSE] + alpha*delta
    new_point$precision_gain = 0
    new_point$recall_gain = recall_gain(new_point$TP,new_point$FN,new_point$FP,new_point$TN)
    new_point$is_crossing = 1
    crossings = rbind(crossings,new_point)
  }
  # add crossing points to the 'points' data frame
  points = rbind(points,crossings)
  points = points[order(points$index,points$recall_gain),2:ncol(points),drop=FALSE]
  rownames(points) = NULL
  points$in_unit_square = 1
  points$in_unit_square[points$recall_gain<0] = 0
  points$in_unit_square[points$precision_gain<0] = 0
  return(points)
}
create_prg_curve = function(labels,pos_scores,neg_scores=-pos_scores,create_crossing_points=TRUE) {
  n = length(labels)
  n_pos = sum(labels)
  n_neg = n - n_pos
  # convert negative labels into 0s
  labels = 1*(labels==1) 
  segments = .create.segments(labels,pos_scores,neg_scores)
  # calculate recall gains and precision gains for all thresholds
  points = data.frame(index=1:(1+nrow(segments)))
  points$pos_score=c(Inf,segments$pos_score)
  points$neg_score=c(-Inf,segments$neg_score)
  points$TP = c(0,cumsum(segments$pos_count))
  points$FP = c(0,cumsum(segments$neg_count))
  points$FN = n_pos-points$TP
  points$TN = n_neg-points$FP
  points$precision_gain = precision_gain(points$TP,points$FN,points$FP,points$TN)
  points$recall_gain = recall_gain(points$TP,points$FN,points$FP,points$TN)
  if (create_crossing_points) {
    points = .create.crossing.points(points,n_pos,n_neg)
  } else {
    points = points[,2:ncol(points)]
  }
  return(points)
}
calc_auprg = function(prg_curve) {
  area = 0
  for (i in 2:nrow(prg_curve)) {
    if (!is.na(prg_curve$recall_gain[i-1]) && (prg_curve$recall_gain[i-1]>=0)) {
      width = prg_curve$recall_gain[i]-prg_curve$recall_gain[i-1]
      height = (prg_curve$precision_gain[i]+prg_curve$precision_gain[i-1])/2
      area = area + width*height
    }
  }
  return(area)
}
get_recall_gain=function(prg_curve) {
return(prg_curve$recall_gain) 
}
get_precision_gain=function(prg_curve){
return(prg_curve$precision_gain)
}
''')

def get_gain_vals(labels,predictions):
    curve_name=r.create_prg_curve(labels,predictions)
    recall_gain=r.get_recall_gain(curve_name) 
    precision_gain=r.get_precision_gain(curve_name) 
    return recall_gain,precision_gain 

def auPRG(labels, predictions):
    return r.calc_auprg(r.create_prg_curve(labels, predictions))[0]

