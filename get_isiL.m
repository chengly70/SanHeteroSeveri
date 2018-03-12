function [isiL,frt]=get_isiL(tm,Volt,thres)
%
%calculates some AP features; ASSUMING tm is equally spaced, same length as
%Volt, and that max(Volt)>thres.  EXCLUDE transient volt; assuming perfectly
%periodic

%initializing outputs
isiL=0;
frt=0;

dt=tm(2)-tm(1);

indHiThres=Volt>thres; %indices where Volt>thres

ind1stX=find(diff(indHiThres)>0)+1; %indices where Volt FIRST crosses thres

isiL=mean(diff(tm(ind1stX))); %use avg ISI
frt=1/isiL;

