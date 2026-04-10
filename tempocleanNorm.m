function [r,f] = tempocleanNorm(LIASet,window)

r = LIASet;
N = length(r(1,:));
regressWin = window;

for ch=1:size(r,1)
    r(ch,:) = normalize(r(ch,:));
end

for in=1:regressWin:(N-regressWin)
    
    XneonL = (r(1,in:in+regressWin-1));
    XtomatoL = (r(2,in:in+regressWin-1));     
    XneonR = (r(3,in:in+regressWin-1));
    XtomatoR = (r(4,in:in+regressWin-1));

    %Perform robust multiple linear regression to remove contribution of
    %tdTomato from mNeon signals
    b = robustfit(XtomatoL,XneonL);
    YFITMNL = b(1) + b(2)*XtomatoL;
    b = robustfit(XtomatoR,XneonR);
    YFITMNR = b(1) + b(2)*XtomatoR;
    
    tempMNL = XneonL - YFITMNL;
    tempMNR = XneonR - YFITMNR;

    %Perform robust multiple linear regression to remove contribution of
    %mNeon from tdTomato signals
    b = robustfit(XneonL,XtomatoL);
    YFITTDL = b(1) + b(2)*XneonL;
    b = robustfit(XneonR,XtomatoR);
    YFITTDR = b(1) + b(2)*XneonR;
    
    tempTDL = XtomatoL - YFITTDL;
    tempTDR = XtomatoR - YFITTDR;
    
    
    j = 0;
    for i=in:in+regressWin-1
        j = j+1;
        r(1,i) = tempMNL(j);
        r(2,i) = tempTDL(j);
        r(3,i) = tempMNR(j);
        r(4,i) = tempTDR(j);
        
        f(1,i) = YFITMNL(j);
        f(2,i) = YFITMNR(j);
    end
    
end
