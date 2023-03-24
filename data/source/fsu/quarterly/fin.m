%fin.m 
clear all
load est_joint_r.mat result tt
load bias_joint_r_run.mat b
table(1,1) = median(result(tt)-b(tt));

load est_fin.mat result tt
load bias_fin.mat b
table(2,1) = median(result(tt)-b(tt));

load est_post04.mat result tt
load bias_post04.mat b
table(3,1) = median(result(tt)-b(tt))