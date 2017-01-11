function [ temppos ] = Quantize_Fiber_Segments(xpos, ypos, minfsegleng)
% Quantized fiber segments for every minimum fiber segment length (minfsegleng = 5 pixels)
% For each quantized fiber length, find the straight line - fitted
% segment. For example, a curvy line is disected by many straight lines. 

fstart=1;
fend=2;
flength=length(xpos(:,1));
temppos=[]; % X1, Y1, X2, Y2, fiber length
while (fend<=flength)
    fdist=norm([xpos(fend) ypos(fend)]-[xpos(fstart) ypos(fstart)],2);
    if(fdist>minfsegleng)
        coeff=polyfit(xpos(fstart:fend),ypos(fstart:fend),1);
        if(abs(coeff(1))<1e6 && xpos(fstart)~=xpos(fend))
            regressionstart=coeff(1)*xpos(fstart)+coeff(2);
            regressionend=coeff(1)*xpos(fend)+coeff(2);
            temppos=[temppos; xpos(fstart) regressionstart xpos(fend) regressionend fdist];
        else
            temppos=[temppos; xpos(fstart) ypos(fstart) xpos(fend) ypos(fend) fdist];
        end
        fstart=fend;
        fend=fend+1;
    else
        fend=fend+1;
    end
end
if (fstart~=flength)
    fdist=norm([xpos(flength) ypos(flength)]-[xpos(fstart) ypos(fstart)],2);
    % round up process for the remaining fiber segment
    % if the fiber segment length > 2.5, we include the fiber segment
    % if not, exclude the fiber segment
    if(fdist>minfsegleng/2)
        fend=flength;
        coeff=polyfit(xpos(fstart:fend),ypos(fstart:fend),1);
        if(abs(coeff(1))<1e6 && xpos(fstart)~=xpos(fend))
            regressionstart=coeff(1)*xpos(fstart)+coeff(2);
            regressionend=coeff(1)*xpos(fend)+coeff(2);
            temppos=[temppos; xpos(fstart) regressionstart xpos(fend) regressionend fdist];
        else
            temppos=[temppos; xpos(fstart) ypos(fstart) xpos(fend) ypos(fend) fdist];
        end
    end
end

end

