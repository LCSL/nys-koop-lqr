function VID = makeVideo(phiPositions,T,nodos_borde,peli,vel_vid,dt)
%triangulamos para hacer el plot
[Xtri,Ttri] =  traingulateQuadMesh(phiPositions{1},T); 
[super,curve] = plotMesh(Xtri,Ttri,nodos_borde);
n_frames = size(phiPositions,2);
clear VID
if length(peli) > 0
%     set(gca, 'nextplot', 'replacechildren');
    VID(1) = getframe;
end  
VID(1) = getframe;
%actualizamos la figura
for i=2:vel_vid:n_frames
    Xnew = phiPositions{i};
    [super,curve] = updatePlot(super,curve,Xnew,T,nodos_borde);
    pause(dt)
    %w = waitforbuttonpress;
    if length(peli) > 0
        if vel_vid == 1
           VID(i) = getframe;
        else   
           VID(floor(i/vel_vid)+1) = getframe;
        end
    end 
end 

if length(peli) > 0
    %guardamos el video
    v = VideoWriter(peli);
    v.Quality = 95;
    v.FrameRate = 30;
    open(v)
    writeVideo(v,VID)
    close(v)
end    