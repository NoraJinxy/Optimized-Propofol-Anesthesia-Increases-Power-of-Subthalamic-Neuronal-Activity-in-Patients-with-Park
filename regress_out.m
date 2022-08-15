function [r, b, bint, stats] = regress_out(Y, conf)
% Y =[n,1], conf=[n,q], there are q confound variables
    [M,N] = size(Y);
    Y = reshape(Y, [length(Y),1]);
    X = [ones(size(conf,1),1), conf];
%     X=conf;
    [b,bint,r,rint,stats] = regress(Y, X);
    r=r+b(2);
    r = reshape(r, [M,N]);
end
