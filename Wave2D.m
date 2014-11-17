% Description:
%     A function to approximate solutions to the the planar wave equation 
%     u_tt=e*(u_xx+uyy) for some scalar e on the bounds -a<x<a, -b<y<b with 
%     initial conditions u(x,y,0)=f(x,y), u_t(x,y,0)=g(x,y) and boundary value 0
% Inputs:
%     bounds:     array [a,b] for bounds -a<x<a, -b<y<b
%     deltas:     scalar dx=dy
%     tend:       ending t bound

function [approx,frames]=Wave2D(bounds,delta,tend)
    a=bounds(1); b=bounds(2); dt=delta/2; dy=delta; dx=delta;
    xn=floor(2*a/dx); yn=floor(2*b/dy); tn=floor(tend/dt); L=(dt/dx); 
    [meshx,meshy]=meshgrid(-a:dx:a,-b:dy:b);
    frames=[];
    %create initial values from initial conditions
    approx=zeros(xn+1,yn+1,tn+1); xs=-a:dx:a; ys=-b:dy:b;
    set(gcf,'Renderer','opengl');
    surf(meshx,meshy,approx(:,:,1));
    axis([-pi,pi,-pi,pi,-1,1]);
    drawnow;
    %pause(1/60);
    %use the McLaurin expansion method to find an aproximation for t1
    for j=2:yn
        for i=2:xn
            approx(i,j,2)=approx(i,j,1)+dt*f(xs(i),ys(j));
        end
    end
    surf(meshx,meshy,approx(:,:,2));
    axis([-pi,pi,-pi,pi,-1,1]);
    drawnow;
%     frames=[frames;getframe()];
%     time_msg=sprintf('Current time:\t%f\n',dt);
%     fprintf(time_msg);
    %pause(1/60);
%     min_span=abs(max(max(approx(:,:,2)))-min(min(approx(:,:,2))));
    %loop and create successive approximations until t=tend using our
    %finite differene scheme
    ti=2;
    while true
        for j=2:yn
            for i=2:xn
                approx(i,j,ti+1)=(2-4*(L^2))*approx(i,j,ti)+(L^2)*...
                    (approx(i+1,j,ti)+approx(i-1,j,ti))+(L^2)*(approx(i,j+1,ti)+approx(i,j-1,ti))...
                    -approx(i,j,ti-1);
            end
        end
        clf;
        surf(meshx,meshy,approx(:,:,ti+1), 'EdgeAlpha', .5);
        axis([-pi,pi,-pi,pi,-.5,1]);
        drawnow;
%         %frames=[frames;getframe()];
%         fprintf(repmat('\b',1,length(time_msg)));
%         cur_span=abs(max(max(approx(:,:,ti)))-min(min(approx(:,:,ti))));
%         time_msg=sprintf('Current time:\t%f\nCurrent span:\t%f\n',ti*dt,cur_span);
%         fprintf(time_msg);
%         if abs(max(max(approx(:,:,ti)))-min(min(approx(:,:,tcur))))<=min_span
%             min_span=abs(max(max(approx(:,:,ti)))-min(min(approx(:,:,ti))));
%         end
%         
        %pause(1/60);
        if ti==tn
            break;
        end
        ti=ti+1;
    end
    
    
    function sol=g(x,y)
        sol=0;
    end
    function sol=f(x,y)
        sol=cos(x)*cos(y);
    end
end
