function [vorts] = loadVtx(start, step, fin, lostCheck)
%loads the processed vorts from the output of vort.py, and structures them
%into Matlab structs for Vortex position and index
%   start step and fin are the initial dataset, step size, and final
%   dataset to load for the vortex postions
%   lostCheck turns on the ability to test for vortex losses and ensure
%   that these vortices do not contribute to the dataset.
%   The Data should be in the format of [X Y Sign UID IsLost] 

    %vorts=zeros(1,1,2); 
    a=zeros(2,1);
    steps=start:step:fin;
   
    vtxLost.count=0;
    vtxLost.idx=[];
    
    %Check if a vortex vanishes (goes outside the boundary) over
    %the course of the simulations. If so, we can choose to discard it.
    %This is indicated by a 0 in the 5th column, whose UID is then tagged
    if lostCheck==1
        for ii=steps
            f=csvread(strcat(strcat('vort_ord_',int2str(ii)),'.csv'),0,0);
            for jj=1:size(f,1)
                if f(jj,5) == 0
                    vtxLost.idx = union(vtxLost.idx,f(jj,4));
                end
            end
        end
        vtxLost.count = length(vtxLost.idx);
    end
    
    %Arrange vortices in appropriate ordering and structuring for
    %trajectory plotting which are return by vorts data structure
    for ii=steps
        count((ii/1000))=0;
        f=csvread(strcat(strcat('vort_ord_',int2str(ii)),'.csv'),0,0);
        for jj=1:(size(f,1))
            if ( round(f(jj,1) == 0) || round(f(jj,2) == 0)) || ismember(f(jj,4),vtxLost.idx)
                0;
            else
                count(ii/1000  ) = count(ii/1000 ) +1;
                vorts(count(ii/1000 ),ii/1000 ).x = f(jj,1);
                vorts(count(ii/1000 ),ii/1000 ).y = f(jj,2);
                vorts(count(ii/1000 ),ii/1000 ).uid = round(f(jj,4));
            end
        end
    end
end