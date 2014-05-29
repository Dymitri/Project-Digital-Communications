function [ inv_map ] = rec_constell( received, map )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tmp_map=2.*round((received+1)/2)-1;
maxm=max(map(:,2));
minm=min(map(:,2));

subM = tmp_map;
subM(subM>maxm) = maxm ;
subM(subM<minm) = minm ;
inv_map(:,:) = subM ;
end

