close all
clear all
clc

%%

%Hill_function=@(X, KM, Nh)((X.^Nh)./(KM.^Nh+X.^Nh))

%Hill_function=@(X, KM, Nh, Max, Min)Min+(Max-Min).*(1./(1+((X.^Nh)./(KM.^Nh))))


Hill_function=@(X, KM1,KM2, Nh, Max, Min)((Min+(((X.^Nh)./(KM.^Nh+X.^Nh))-Min)).*(1./(1+((X.^Nh)./(KM.^Nh))))

KM1=1
KM2=2
Nh=10
Max=10
Min=1
X=0:0.1:KM2*2


plot(X, Hill_function(X, KM1,KM2, Nh, Max, Min))
hold on
line([KM1, KM1], [0, Max])
line([0, X(end)], [Max/2, Max/2])
line([0, X(end)], [Max, Max])