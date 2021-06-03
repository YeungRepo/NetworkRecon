function[movie_handle] = Gen_Pretty_DSFGraphs(Qoft,ColorMat_Matrix,Nodes_Names,THRESHOLD,This_Circuit_Name)
close all
figure(1)
vidfile = VideoWriter(strcat('Qoft_',This_Circuit_Name,'.mp4'),'MPEG-4');
open(vidfile); 

for tindex = 2:size(Qoft,3)
    DQ = digraph(Qoft(:,:,tindex)','omitselfloops');
    try
        LWidths = 20.0*abs(DQ.Edges.Weight+eps)/max(abs(DQ.Edges.Weight+eps)); 
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
        labelnode(This_Plot_Handle,[1 2],Nodes_Names)
        writeVideo(vidfile, getframe(gcf));   
        
    catch e %e is an MException struct
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
    
end
movie_handle = vidfile;
close(vidfile);

