function[movie_handle] = Gen_Pretty_DSFGraphs(Qoft,ColorMat_Matrix,THRESHOLD)
close all
figure(1)
vidfile = VideoWriter('Qoft.mp4','MPEG-4');
open(vidfile); 


for tindex = 2:3
    DQ = digraph(Qoft(:,:,tindex),'omitselfloops');
    try
        h_plot = plot(DQ,'EdgeLabel',round(DQ.Edges.Weight,2),'LineWidth',LWidths,'EdgeFontSize',20,'MarkerSize',50,'NodeFontSize',50,'ArrowPosition',0.75,'ArrowSize',75,'NodeColor',ColorMat_Matrix)
        %Fixed_Node_Label={'Y1','Y2','Y3'};
        Fixed_XData = h_plot.XData; 
        Fixed_YData = h_plot.YData; 
    catch  
    end

        
end




for tindex = 2:numel(Qoft(1,1,:))
    DQ = digraph(Qoft(:,:,tindex),'omitselfloops');
    %DQ = rmedge(DQ,find(abs(DQ.Edges.Weight)<=THRESHOLD));  
    %display(abs(DQ.Edges.Weight))
    try
        LWidths = 50.0*abs(DQ.Edges.Weight+eps)/max(abs(DQ.Edges.Weight+eps)); 
        EdgeColors = zeros(numel(LWidths),3);
        for Edge_Index  = 1:numel(DQ.Edges.Weight)
            if abs(DQ.Edges.Weight(Edge_Index)) <= THRESHOLD                
                EdgeColors(Edge_Index,:) = [1.0;1.0;1.0];
            else
                if DQ.Edges.Weight(Edge_Index)>0.0
                    EdgeColors(Edge_Index,:) = [0.0;1.0;0.0];
                end
                if DQ.Edges.Weight(Edge_Index)<0.0
                    EdgeColors(Edge_Index,:) = [1.0;0.0;0.0];
                end
            end
        end
        This_Plot_Handle = plot(DQ,'EdgeColor',EdgeColors,'LineWidth',LWidths,'EdgeFontSize',20,'MarkerSize',50,'NodeFontSize',50,'ArrowPosition',0.75,'ArrowSize',75,'NodeColor',ColorMat_Matrix);
        %This_Plot_Handle.NodeLabel = Fixed_Node_Label;
        %This_Plot_Handle.XData = Fixed_XData; 
        %This_Plot_Handle.YData = Fixed_YData; 
        writeVideo(vidfile, getframe(gcf));        
    catch e %e is an MException struct
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
    
end
movie_handle = vidfile;
close(vidfile);

