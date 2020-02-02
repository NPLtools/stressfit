%%Slit
%Matrix model of a slit or radial collimator (distance=0)
% 

classdef TofDetector < Slit
    properties(Constant)     
    end
    properties  
        amin,amax;  % angular range
        dtau;  % time resolution (gaussian fwhm)
    end

    methods              
        function X=TofDetector(distance,binwidth,binshape,dtau,amin,amax)
           X = X@Slit(distance, binwidth*binshape/Comp.g_uni);
           X.amin=amin;
           X.amax=amax;
           X.dtau=dtau*Comp.g_norm;
           X.id='TofDetector';
        end         
    end
end

