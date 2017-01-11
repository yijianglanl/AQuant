function [ fcount, temp_AI, temp_ang] = Calculate_Alignment_Vector(angpos, x, y, distcheck)
fcount=0;
temp_AI=0;
temp_ang=0;
coslist=[];
for pl=1:length(angpos(:,1))
     fdist1=norm([angpos(pl,1) angpos(pl,2)]-[x y],2);
     fdist2=norm([angpos(pl,3) angpos(pl,4)]-[x y],2);
     % Accept the quantized fiber segment for the given local circle position (x, y) 
     % to calculate the local alignment vector, 
     % if any end points of the fiber segment are located inside of the
     % local circle
     if (fdist1<distcheck || fdist2<distcheck)
         fcount=fcount+1;
         fseg=[angpos(pl,3)-angpos(pl,1), angpos(pl,4)-angpos(pl,2)];
         if (fseg(1,2)<0)
            fseg=-1.*fseg;
         end
         cosangle=fseg(1,1)/norm(fseg,2);
         coslist=[coslist; cosangle];
     % Accept the quantized fiber segment for the given local circle position (x, y) 
     % to calculate the local alignment vector, 
     % if the fiber segment line pass through the local circle, although
     % two end points of the line are located outside of the local circle.
     elseif (fdist1<(distcheck+5) || fdist2<(distcheck+5))
         x1=angpos(pl,1); y1=angpos(pl,2);
         x2=angpos(pl,3); y2=angpos(pl,4);
         tempxden=(y2-y1)*(y2-y1)+(x1-x2)*(x1-x2);
         tempxnum=(x1-x2)*((x1-x2)*x-(y2-y1)*y)-(y2-y1)*(y1*(x2-x1)-(y2-y1)*x1);
         tempx=tempxnum/tempxden;
         if (tempx>x1 && tempx<x2)
             checkdist=1;
         elseif(tempx>x2 && tempx<x1)
             checkdist=1;
         else
             checkdist=0;
         end
         if (checkdist==1)
             distden=tempxden^0.5;
             distnum=abs((x2-x1)*(y1-y)-(x1-x)*(y2-y1));
             cdist=distnum/distden;
             if (cdist<distcheck)
                 fcount=fcount+1;
                 fseg=[angpos(pl,3)-angpos(pl,1), angpos(pl,4)-angpos(pl,2)];
                 if (fseg(1,2)<0)
                    fseg=-1.*fseg;
                 end
                 cosangle=fseg(1,1)/norm(fseg,2);
                 coslist=[coslist; cosangle];
             end
         end
     end
end
if(fcount>0)
    % All fiber angles are mapped into upper plane (0 - 180 degree)
    tempanghist=acos(coslist);
    
    % Fibers are mapped into the unit circle (0 - 360 degree)
    tempanghist2=tempanghist.*2;
    cosanghist=cos(tempanghist2);
    sinanghist=sin(tempanghist2);
    
    %mrv = mean resultant vector in Circular Statistics
    mrv=[sum(cosanghist) sum(sinanghist)]./length(tempanghist(:,1));
    
    %Alignment index of the local alignment vector
    temp_AI=norm(mrv,2);
    mvector=mrv./temp_AI;
    
    %Local alignment vector angle, 0 - 180 degree
    temp_ang=acos(mvector(1,1));
    if (mvector(1,2)>0)
        temp_ang=temp_ang/2;
    else
        temp_ang=(2*pi-temp_ang)/2;
    end
end

end

