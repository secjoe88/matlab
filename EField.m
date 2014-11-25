% Description:
%     A script to display an animated graphic of a time-dependent electric field 
%     in free space. Currently, the user may set different source charge locations,
%     and electric charges (as a function of t), and the script generates the 
%     time-dependent field in real time. Adjustments can be made for frequency
%     and charge (although all electrostatic force vectors are drawn as unit vectors).
function EField()
%create variables
    framerate=30;
    %time step
	dt=pi/framerate;
	t=0;
    %location of source charges
	Ps=zeros(6,3);Ps(1:3,3)=.5:.5:1.5; Ps(4:6,3)=-.5:-.5:-1.5;
    ps=zeros(length(Ps(:,1)),1);
    maxF=0;
    %create grid of test particles with test charge -1C
	[Qx,Qy,Qz]=genQArray(.25,10,10);q=1;
    %create initial 0 E-field
	Ex=zeros(size(Qx));Ey=Ex;Ez=Ex;
    
%Create figure handles for plots. We will update these with refreshed data on the loop
%below
    set(gcf,'Renderer','opengl','Color',[1 1 1],'CloseRequestFcn','delete(gcf)');
    %figure handle for E-Field Plot
    h=quiver3(Qx,Qy,Qz,Ex,Ey,Ez,'AutoScale','off');
    hold on;
    %figure handles for soures
    SourcePlots=zeros(length(Ps(:,1)),1);
    for i=1:length(SourcePlots)
        SourcePlots(i)=plot3(Ps(i,1),Ps(i,2),Ps(i,3),'LineStyle','none','Marker','.','HandleVisibility','off');
    end
    axis([-4,4,-4,4,-4,4]);
    set(gca,'Projection','perspective','CameraPosition',[-22.038 -54.544 36.768],...
        'GridLineStyle','none','Visible','off');
    axis vis3d;zoom(2.5);
%Generate animation label and title
    %plot invisible static markers so there is something to put on the
    %label
    plot3(50,50,50,'LineStyle','none','Marker','.','MarkerSize',40,...
        'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'Visible','off');
    plot3(50,50,50,'LineStyle','none','Marker','.','MarkerSize',40,...
        'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Visible','off');
    plot3(50,50,50,'LineStyle','none','Marker','.','MarkerSize',40,...
        'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'Visible','off');
    key=legend('Electrostatic Force in Free Space','+1C Source Charge',...
        '0C Source Charge','-1C Source Charge','Location','southeast');
    set(key,'Position',[.650 .029 .336 .106]);
    title=annotation('textbox','LineStyle','none','FontSize',14,'FontWeight',...
        'bold');
    set(title,'String','Electric Field as Would be Experienced by a 1C Point Charge in Free Space',...
        'Position',[0 .892 .185 .108]);
    set(gcf,'Position',[1 1 1280 720]);
    
%Generate t-dependent E-Field
    drawnow;
    startTime=tic;    
    while ishandle(h) %t<=60%
        ps(1:3)=cos(t*dt);ps(4:6)=-cos(t*dt);
        for i=1:numel(Qx)
            %set location of current test charge
            Q=[Qx(i),Qy(i),Qz(i)];
            %make sure test charge isnt the same as source charge
            if checkSource()
                Ex(i)=0;Ey(i)=0;Ez(i)=0;
            end
            %calculate force exerted on Q by field
            F=genForceVec(Ps,ps,Q,q);
            if norm(F)>maxF
                maxF=norm(F);
            end
            %store the force vector in its place in the vector field E
            F=.2*(F/norm(F))*abs(cos(t*dt));
            Ex(i)=F(1);Ey(i)=F(2);Ez(i)=F(3);
        end
        stopTime=toc(startTime);
        %render at 30 fps
        if stopTime>=1/framerate
            updateSources();
            set(h,'UData',Ex,'VData',Ey,'WData',Ez);
            camorbit(18/framerate,0);
            %drawnow;
            %optionally save frames to file
%             filename=sprintf('frames/frame%d.png',t);
%             saveas(gcf,filename);
            startTime=tic;
            
        end
        t=t+1;
    end
    %check if Q is a source charge return true if it is
    function bool=checkSource()
        bool=0;
        for j=1:length(Ps(:,1))
            if isequal(Ps(j,:),Q)
                bool=1;
            end
        end
    end
    %plot the sources charges and vary their color based on their electric 
    %charges
    function updateSources()
        for j=1:length(Ps(:,1))
            if ps(j)>0
                color=[abs(ps(j)) 0 1-abs(ps(j))];
            else
                color=[0 0 1-abs(ps(j))];
            end
            set(SourcePlots(j),'MarkerSize',40,...
                'MarkerFaceColor',color,'MarkerEdgeColor',color);
        end
    end
end
%Generate MxM sphere of test charges around the origin at distances
%r1+n*dr from the origin for n=0,...,N 
function [Qx,Qy,Qz]=genQArray(dr,N,M)
    %initialize test charge arrays
    Qx=[];Qy=[];Qz=[];
    for i=1:N
        R=(i*dr);
        [tx,ty,tz]=sphere(M);
        Qx=[Qx,1.5*R*tx];Qy=[Qy,1.5*R*ty];Qz=[Qz,R*tz];
    end
end
%Takes in an array Ps=[P1;P2;...;Pn] of coordinates of source charges with
%Pi=[xi,yi,zi], an array ps=[p1;p2;...;pn] of their respective electric 
%charges, coordinates Q=[Qx,Qy,Qz] of a test charge, and its electric
%charge q. Returns the force vector F that the electric Field exerts on Q.
function F=genForceVec(Ps,ps,Q,q)
    F=[0,0,0];
    for i=1:length(Ps(:,1))
        %for each P
        curP=Ps(i,:);
        %calculate the unit vector R:Q->P
        R=-(curP-Q)/norm(curP-Q);
        %use R to calculate the force FTemp that the P exerts on Q
        FTemp=(ps(i)*q)*R/(norm(curP-Q)^2);
        %add FTemp to F, the net force exerted on Q by the field
        F=F+FTemp;
    end
end






